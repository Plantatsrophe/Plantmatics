#!/usr/bin/env Rscript
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Script: scripts/analysis/edaphic_solus_3dep_utils.R
# Project: Plantmatics / Packera Species Delimitation Pipeline (NCU)
# Description:
#   SOLUS100 100 m lithic/soil extraction, 30 m 3DEP micro-topography (Slope,
#   TPI, McCune & Keon 2002 HLI, D8 TWI), USGS SGMC 6-class lithology, and
#   idempotent fine-grain SQLite/CSV cache orchestration (pure POSIX ERE).
# ==============================================================================

.normalize_soilgrids_prop <- function(v_raw, prop_type = "pct") {
  v <- suppressWarnings(as.numeric(v_raw)); fin <- v[is.finite(v)]
  if (length(fin) == 0L) return(rep(NA_real_, length(v)))
  med_v <- stats::median(fin, na.rm = TRUE)
  if (identical(prop_type, "ph")) {
    if (is.finite(med_v) && med_v > 100.0) v <- v / 100.0 else if (is.finite(med_v) && med_v > 14.5) v <- v / 10.0
    v[!is.finite(v) | v < 2.0 | v > 11.5] <- NA_real_
  } else if (identical(prop_type, "pct")) {
    if (is.finite(med_v) && med_v > 105.0) v <- v / 10.0
    v[!is.finite(v) | v < 0.0 | v > 100.0] <- NA_real_
  } else if (identical(prop_type, "cec")) {
    if (is.finite(med_v) && med_v > 90.0) v <- v / 10.0
    v[!is.finite(v) | v < 0.0 | v > 250.0] <- NA_real_
  } else if (identical(prop_type, "bd")) {
    if (is.finite(med_v) && med_v > 10.0) v <- v / 100.0
    v[!is.finite(v) | v < 0.1 | v > 2.65] <- NA_real_
  }
  v
}

extract_polaris_and_solus_predictors <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Extracting harmonized 0-15cm rhizosphere soil & lithic predictors (SOLUS100 / SSURGO / SoilGrids 250m; POLARIS 30m disabled)...")
  df <- flag_coordinate_uncertainty(coords_df); n <- nrow(df); if (n == 0L) return(df)
  gssurgo_dir <- file.path(env_dir, "gssurgo_30m"); solus_dir <- file.path(env_dir, "solus_100m")
  if (!isTRUE(demo)) dir.create(solus_dir, recursive = TRUE, showWarnings = FALSE)

  init_cols <- c("polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m",
                 "solus_resdept_cm", "solus_fragvol_pct", "solus_sand_pct", "solus_clay_pct", "solus_silt_pct", "solus_cec",
                 "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "soil_silt_pct", "ilr_texture_1", "ilr_texture_2",
                 "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm")
  for (col in init_cols) if (!(col %in% names(df))) df[[col]] <- rep(NA_real_, n)
  df <- sanitize_disabled_polaris_columns(df)

  lat_v <- suppressWarnings(as.numeric(df$latitude)); lon_v <- suppressWarnings(as.numeric(df$longitude))
  in_conus <- !is.na(lat_v) & !is.na(lon_v) & lat_v >= 24.5 & lat_v <= 49.5 & lon_v >= -125.0 & lon_v <= -66.5
  solus_ph_vec <- rep(NA_real_, n)

  if (!isTRUE(demo) && any(in_conus) && requireNamespace("terra", quietly = TRUE)) {
    c_idx <- which(in_conus)
    pts_v <- terra::vect(df[c_idx, c("longitude", "latitude"), drop = FALSE], geom = c("longitude", "latitude"), crs = "EPSG:4326")
    ex_sol <- c(list.files(solus_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE),
                if (dir.exists(gssurgo_dir)) list.files(gssurgo_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE) else character(0))
    if (length(ex_sol) > 0L) {
      ext_1 <- function(rp) tryCatch({ r <- terra::rast(rp); as.numeric(terra::extract(r, .project_vect_to_rast(pts_v, r))[[2L]]) }, error = function(e) rep(NA_real_, length(c_idx)))
      ext_0_15 <- function(cf, tf = identity) {
        if (length(cf) == 0L) return(rep(NA_real_, length(c_idx)))
        bn <- basename(cf)
        f05 <- cf[grepl("0_5|0-5", bn, ignore.case = TRUE)]; f515 <- cf[grepl("5_15|5-15", bn, ignore.case = TRUE)]
        f0 <- cf[grepl("(^|_)0_cm(_|\\.)", bn, ignore.case = TRUE)]; f5 <- cf[grepl("(^|_)5_cm(_|\\.)", bn, ignore.case = TRUE)]; f15 <- cf[grepl("(^|_)15_cm(_|\\.)", bn, ignore.case = TRUE)]
        if (length(f05) > 0L && length(f515) > 0L) { compute_rhizosphere_0_15cm(tf(ext_1(f05[1L])), tf(ext_1(f515[1L]))) } else if (length(f0) > 0L && length(f5) > 0L && length(f15) > 0L) { compute_rhizosphere_0_15cm(0.5 * (tf(ext_1(f0[1L])) + tf(ext_1(f5[1L]))), 0.5 * (tf(ext_1(f5[1L])) + tf(ext_1(f15[1L])))) } else { tf(ext_1(cf[1L])) }
      }
      bn_all <- basename(ex_sol)
      # Extract SOLUS100 bedrock strictly from lithic/paralithic contact layers (anylithicdpt)
      lith_f <- ex_sol[grepl("anylithicdpt", bn_all, ignore.case = TRUE)]
      if (length(lith_f) > 0L) { vr <- ext_1(lith_f[1L]); ok <- is.finite(vr) & vr >= 0.0; if (any(ok)) df$solus_resdept_cm[c_idx[ok]] <- round(pmin(200.0, pmax(0.0, vr[ok])), 1) }
      vf <- ext_0_15(ex_sol[grepl("fragvol|cfvo", bn_all, ignore.case = TRUE)])
      if (any(is.finite(vf) & vf >= 0.0)) { ok <- is.finite(vf) & vf >= 0.0; df$solus_fragvol_pct[c_idx[ok]] <- round(pmin(100.0, pmax(0.0, vf[ok])), 1) }
      vp <- ext_0_15(ex_sol[grepl("^ph1to1h2o|^ph_solus", bn_all, ignore.case = TRUE)], .rescale_solus_ph_vec)
      if (any(is.finite(vp))) solus_ph_vec[c_idx[is.finite(vp)]] <- round(vp[is.finite(vp)], 2)
      for (sp in list(c("^sandtotal", "solus_sand_pct"), c("^claytotal", "solus_clay_pct"), c("^silttotal", "solus_silt_pct"))) {
        vv <- ext_0_15(ex_sol[grepl(sp[1L], bn_all, ignore.case = TRUE)], function(z) .normalize_soilgrids_prop(z, "pct"))
        if (any(is.finite(vv))) df[[sp[2L]]][c_idx[is.finite(vv)]] <- round(vv[is.finite(vv)], 1)
      }
      vc <- ext_0_15(ex_sol[grepl("^cec7", bn_all, ignore.case = TRUE)], function(z) .normalize_soilgrids_prop(z / 10.0, "cec"))
      if (any(is.finite(vc))) df$solus_cec[c_idx[is.finite(vc)]] <- round(vc[is.finite(vc)], 2)
    }
  }

  get_sg <- function(c1, c2, cb, pt) {
    x1 <- if (c1 %in% names(df)) suppressWarnings(as.numeric(df[[c1]])) else rep(NA_real_, n)
    x2 <- if (c2 %in% names(df)) suppressWarnings(as.numeric(df[[c2]])) else if (cb %in% names(df)) suppressWarnings(as.numeric(df[[cb]])) else rep(NA_real_, n)
    .normalize_soilgrids_prop(compute_rhizosphere_0_15cm(x1, x2), prop_type = pt)
  }
  sg_ph <- get_sg("soil_ph_0_5cm", "soil_ph_5_15cm", "soil_ph", "ph"); sg_sand <- get_sg("soil_sand_0_5cm", "soil_sand_5_15cm", "soil_sand", "pct")
  sg_clay <- get_sg("soil_clay_0_5cm", "soil_clay_5_15cm", "soil_clay", "pct"); sg_silt <- get_sg("soil_silt_0_5cm", "soil_silt_5_15cm", "soil_silt", "pct")
  sg_cec <- get_sg("soil_cec_0_5cm", "soil_cec_5_15cm", "soil_cec", "cec"); sg_cfvo <- get_sg("soil_cfvo_0_5cm", "soil_cfvo_5_15cm", "soil_cfvo", "pct")
  gn <- function(cn) if (cn %in% names(df)) suppressWarnings(as.numeric(df[[cn]])) else rep(NA_real_, n)
  ss_ph <- gn("ssurgo_ph"); ss_sand <- gn("ssurgo_sand_pct"); ss_clay <- gn("ssurgo_clay_pct"); ss_cec <- gn("ssurgo_cec"); ss_bed <- gn("ssurgo_bedrock_depth_cm"); ss_out <- gn("ssurgo_rock_outcrop_pct")
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n); is_out[is.na(is_out)] <- FALSE

  if (isTRUE(demo)) {
    mr <- is.na(df$solus_resdept_cm); if (any(mr)) df$solus_resdept_cm[mr] <- ifelse(is.finite(ss_bed[mr]), round(ss_bed[mr], 1), 160.0)
    mf <- is.na(df$solus_fragvol_pct)
    if (any(mf)) { ss_sim <- ifelse(is.finite(ss_out[mf]), ss_out[mf], ifelse(is_out[mf], 40.0, 0.0)); df$solus_fragvol_pct[mf] <- round(pmin(85.0, pmax(1.0, 0.55 * ss_sim + ifelse(df$solus_resdept_cm[mf] < 100, (100 - df$solus_resdept_cm[mf]) * 0.25, 4.0))), 1) }
    mc <- is.na(df$solus_cec); if (any(mc)) df$solus_cec[mc] <- ifelse(is.finite(ss_cec[mc]), round(ss_cec[mc], 2), ifelse(is.finite(sg_cec[mc]), round(sg_cec[mc], 2), 11.5))
  }

  coalesce_num <- function(...) {
    out <- rep(NA_real_, n)
    for (v in list(...)) { vn <- suppressWarnings(as.numeric(v)); fill <- !is.finite(out) & is.finite(vn); if (any(fill)) out[fill] <- vn[fill] }
    out
  }
  df <- sanitize_disabled_polaris_columns(df)
  df$fine_soil_ph <- round(coalesce_num(solus_ph_vec, ss_ph, sg_ph), 2)
  df$fine_soil_sand <- round(coalesce_num(df$solus_sand_pct, ss_sand, sg_sand), 1)
  df$fine_soil_clay <- round(coalesce_num(df$solus_clay_pct, ss_clay, sg_clay), 1)
  df$fine_soil_cec <- round(coalesce_num(df$solus_cec, ss_cec, sg_cec), 2)
  df$fine_cfvo_pct <- round(pmin(100.0, pmax(0.0, coalesce_num(df$solus_fragvol_pct, sg_cfvo))), 1)
  bed_comb <- coalesce_num(df$solus_resdept_cm, ss_bed)
  df$fine_bedrock_depth_cm <- ifelse(is.finite(bed_comb), round(pmin(200.0, pmax(0.0, bed_comb)), 1), NA_real_)
  tex_lr <- compute_soil_texture_log_ratios(df$fine_soil_sand, df$fine_soil_clay, silt_pct = coalesce_num(df$solus_silt_pct, sg_silt))
  df$soil_silt_pct <- tex_lr$soil_silt_pct; df$ilr_texture_1 <- tex_lr$ilr_texture_1; df$ilr_texture_2 <- tex_lr$ilr_texture_2
  df
}

.compute_mccune_keon_hli <- function(lat_deg, slope_deg, aspect_deg) {
  rad <- pi / 180.0; L <- pmin(60.0, pmax(15.0, abs(as.numeric(lat_deg)))) * rad; S <- pmin(85.0, pmax(0.0, as.numeric(slope_deg))) * rad
  F_rad <- abs(pi - abs((as.numeric(aspect_deg) %% 360.0) * rad - (225.0 * rad)))
  hli <- exp(-1.467 + 1.582 * cos(L) * cos(S) - 1.500 * cos(F_rad) * sin(S) * sin(L) - 0.262 * sin(L) * sin(S) + 0.607 * sin(F_rad) * sin(S))
  ok <- is.finite(lat_deg) & is.finite(slope_deg); flat <- ok & (!is.finite(aspect_deg) | slope_deg < 1e-4)
  if (any(flat)) hli[flat] <- exp(-1.467 + 1.582 * cos(L[flat]))
  hli[!ok] <- NA_real_; pmin(1.50, pmax(0.05, hli))
}

.extract_dem_tile_metrics <- function(r_dem, coords_sub, rows_idx) {
  if (is.null(r_dem) || length(rows_idx) == 0L) return(NULL)
  tryCatch({
    pts_p <- .project_vect_to_rast(terra::vect(coords_sub, geom = c("longitude", "latitude"), crs = "EPSG:4326"), r_dem)
    ep <- terra::ext(pts_p); rxy <- terra::res(r_dem); ill <- isTRUE(tryCatch(terra::is.lonlat(r_dem), error = function(e) TRUE))
    px <- max(rxy[1L] * 65.0, if (ill) 0.018 else 1800.0); py <- max(rxy[2L] * 65.0, if (ill) 0.018 else 1800.0)
    rs <- tryCatch(terra::crop(r_dem, terra::ext(ep[1L] - px, ep[2L] + px, ep[3L] - py, ep[4L] + py)), error = function(e) NULL)
    if (is.null(rs) || terra::ncell(rs) < 9L) return(NULL)
    fd <- tryCatch(terra::terrain(rs, v = "flowdir"), error = function(e) NULL)
    fa <- if (!is.null(fd)) tryCatch(terra::flowAccumulation(fd), error = function(e) NULL) else NULL
    list(rows = rows_idx, is_lonlat = ill, res_xy = rxy, elev = as.numeric(terra::extract(rs, pts_p)[[2L]]),
         tpi = as.numeric(terra::extract(rs - terra::focal(rs, w = 5, fun = mean, na.rm = TRUE), pts_p)[[2L]]),
         slope = as.numeric(terra::extract(terra::terrain(rs, v = "slope", unit = "degrees"), pts_p)[[2L]]),
         asp = as.numeric(terra::extract(terra::terrain(rs, v = "aspect", unit = "degrees"), pts_p)[[2L]]),
         facc = if (!is.null(fa)) as.numeric(terra::extract(fa, pts_p)[[2L]]) else rep(NA_real_, length(rows_idx)))
  }, error = function(e) NULL)
}

extract_micro_topography_3dep <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Computing 30m 3DEP micro-topographic covariates (Slope, TPI, HLI, D8 Hydrological TWI)...")
  df <- as.data.frame(coords_df, stringsAsFactors = FALSE); n <- nrow(df); if (n == 0L) return(df)
  dem_dir <- file.path(env_dir, "dem_3dep_30m")
  if (!isTRUE(demo)) dir.create(dem_dir, recursive = TRUE, showWarnings = FALSE)
  for (col in c("topo_tpi_30m", "topo_slope_deg", "topo_hli_30m", "topo_twi_30m", "elevation_m")) if (!(col %in% names(df))) df[[col]] <- rep(NA_real_, n)

  lat_v <- suppressWarnings(as.numeric(df$latitude)); lon_v <- suppressWarnings(as.numeric(df$longitude)); valid_pts <- is.finite(lat_v) & is.finite(lon_v)
  stale_dem <- file.path(dem_dir, "dem_3dep_regional.tif"); if (!isTRUE(demo) && file.exists(stale_dem)) tryCatch(unlink(stale_dem), error = function(e) NULL)

  apply_tile <- function(to) {
    if (is.null(to)) return(invisible(NULL))
    rs <- to$rows; ls <- lat_v[rs]
    cw <- if (to$is_lonlat) sqrt(pmax(1.0, (to$res_xy[1L] * 111320.0 * cos(pmin(80.0, abs(ls)) * pi / 180.0)) * (to$res_xy[2L] * 110574.0))) else sqrt(pmax(1.0, to$res_xy[1L] * to$res_xy[2L]))
    ok_e <- is.finite(to$elev) & is.na(df$elevation_m[rs]); if (any(ok_e)) df$elevation_m[rs[ok_e]] <<- round(to$elev[ok_e], 1)
    ok_t <- is.finite(to$tpi); if (any(ok_t)) df$topo_tpi_30m[rs[ok_t]] <<- round(to$tpi[ok_t], 2)
    ok_s <- is.finite(to$slope)
    if (any(ok_s)) {
      sl_c <- pmax(0.05, pmin(75.0, to$slope[ok_s])); df$topo_slope_deg[rs[ok_s]] <<- round(sl_c, 2)
      hli_s <- .compute_mccune_keon_hli(ls[ok_s], sl_c, to$asp[ok_s]); ok_h <- is.finite(hli_s); if (any(ok_h)) df$topo_hli_30m[rs[ok_s][ok_h]] <<- round(hli_s[ok_h], 4)
    }
    ok_w <- is.finite(to$facc) & is.finite(to$slope)
    if (any(ok_w)) df$topo_twi_30m[rs[ok_w]] <<- round(pmin(25.0, pmax(0.5, log(pmax(1e-3, ((pmax(0.0, to$facc[ok_w]) + 1.0) * cw[ok_w]) / pmax(tan(pmax(to$slope[ok_w], 0.10) * (pi / 180.0)), 1e-4))))), 2)
  }

  if (!isTRUE(demo) && any(valid_pts) && requireNamespace("terra", quietly = TRUE)) {
    dem_files <- .list_real_fine_rasters(dem_dir, exclude_pattern = "^dem_3dep_regional\\.tif$")
    if (length(dem_files) > 0L) {
      v_idx <- which(valid_pts); stems <- toupper(sub("\\.(tif|tiff|hgt)$", "", basename(dem_files), ignore.case = TRUE))
      is_grid <- grepl("^[NS][0-9]{2}[EW][0-9]{3}$", stems); coords_xy <- df[, c("longitude", "latitude"), drop = FALSE]
      if (any(is_grid)) {
        tmap <- stats::setNames(dem_files[is_grid], stems[is_grid])
        pt_t <- sprintf("%s%02d%s%03d", ifelse(floor(lat_v[v_idx]) >= 0, "N", "S"), abs(as.integer(floor(lat_v[v_idx]))), ifelse(floor(lon_v[v_idx]) >= 0, "E", "W"), abs(as.integer(floor(lon_v[v_idx]))))
        matched <- intersect(unique(pt_t), names(tmap))
        nw <- if (length(matched) >= 2L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, length(matched)) else 1L
        run_tid <- function(tid) { rs_t <- v_idx[pt_t == tid]; .extract_dem_tile_metrics(tryCatch(terra::rast(tmap[[tid]]), error = function(e) NULL), coords_xy[rs_t, , drop = FALSE], rs_t) }
        for (to in (if (nw > 1L) packera_mclapply(matched, run_tid, mc.cores = nw) else lapply(matched, run_tid))) apply_tile(to)
      }
      rem_idx <- v_idx[is.na(df$topo_tpi_30m[v_idx])]
      for (od_p in c(dem_files[!is_grid], dem_files[is_grid])) {
        if (length(rem_idx) == 0L) break
        r_od <- tryCatch(terra::rast(od_p), error = function(e) NULL); if (is.null(r_od)) next
        ex <- tryCatch(terra::ext(r_od), error = function(e) NULL)
        in_ex <- if (!is.null(ex)) rem_idx[lon_v[rem_idx] >= ex[1L] & lon_v[rem_idx] <= ex[2L] & lat_v[rem_idx] >= ex[3L] & lat_v[rem_idx] <= ex[4L]] else rem_idx
        if (length(in_ex) > 0L) apply_tile(.extract_dem_tile_metrics(r_od, coords_xy[in_ex, , drop = FALSE], in_ex))
        rem_idx <- v_idx[is.na(df$topo_tpi_30m[v_idx])]
      }
    }
  }

  muname_v <- if ("muname" %in% names(df)) as.character(df$muname) else rep("", n)
  slope_ss <- atan(parse_ssurgo_slope_pct(muname_v) / 100.0) * (180.0 / pi)
  if (!isTRUE(demo)) {
    ms <- !is.finite(df$topo_slope_deg) & is.finite(slope_ss); if (any(ms)) df$topo_slope_deg[ms] <- round(pmax(0.15, pmin(60.0, slope_ss[ms])), 2)
    return(mask_micro_scale_predictors_by_uncertainty(df, cols = c("topo_tpi_30m", "topo_hli_30m", "topo_twi_30m"), max_uncertainty_m = 250.0))
  }

  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n); is_out[is.na(is_out)] <- FALSE
  out_p <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(NA_real_, n)
  op_v <- ifelse(is.finite(out_p), out_p, ifelse(is_out, 40.0, 0.0))
  df$topo_slope_deg <- pmax(0.15, pmin(45.0, round(ifelse(!is.na(slope_ss), slope_ss, ifelse(is_out, 4.2, 2.5)), 2)))
  is_conc <- grepl("muck|flooded|aquolls|aquepts|aquults|saprists|depression|seep", tolower(muname_v))
  df$topo_tpi_30m <- round(ifelse(is_out | op_v >= 15.0, 2.8 + 0.045 * op_v, ifelse(is_conc, -2.4, 0.0)), 2)
  df$topo_hli_30m <- round(.compute_mccune_keon_hli(ifelse(is.finite(lat_v), lat_v, 35.5), df$topo_slope_deg, ifelse(is_out, 215.0, 135.0)), 4)
  df$topo_twi_30m <- round(log(ifelse(is_conc, 450.0, ifelse(is_out, 45.0, 120.0)) / tan(pmax(df$topo_slope_deg, 0.15) * (pi / 180.0))), 2)
  mask_micro_scale_predictors_by_uncertainty(df, cols = c("topo_tpi_30m", "topo_hli_30m", "topo_twi_30m"), max_uncertainty_m = 250.0)
}

standardize_sgmc_lithology <- function(
  rocktype1 = "", rocktype2 = "", major1 = "", major2 = "", minor1 = "", unit_name = "",
  muname = "", outcrop_class = "", lat = NA_real_, lon = NA_real_, is_demo = FALSE, demo = NULL,
  generalize = "", rock_type = "", compname = "", taxsuborder = "", taxgreatgroup = "",
  reskind = "", soil_ph = NA_real_, soil_sand = NA_real_, soil_cec = NA_real_,
  bedrock_depth_cm = NA_real_, rock_outcrop_pct = NA_real_
) {
  eff_demo <- if (!is.null(demo)) isTRUE(demo) else isTRUE(is_demo)
  n <- max(length(rocktype1), length(rocktype2), length(major1), length(major2), length(minor1), length(unit_name), length(generalize), length(rock_type), length(muname), length(compname), length(taxsuborder), length(taxgreatgroup), length(reskind), length(outcrop_class), length(soil_ph), length(soil_sand), length(soil_cec), length(bedrock_depth_cm), length(rock_outcrop_pct), length(lat), length(lon), 1L)
  rs <- function(x) { if (is.null(x) || length(x) == 0L) return(rep("", n)); y <- as.character(x); y[is.na(y)] <- ""; rep_len(y, n) }
  rn <- function(x) { if (is.null(x) || length(x) == 0L) return(rep(NA_real_, n)); rep_len(suppressWarnings(as.numeric(x)), n) }

  oc_c <- tolower(trimws(rs(outcrop_class)))
  combined <- trimws(paste(tolower(trimws(paste(rs(rocktype1), rs(rocktype2), rs(major1), rs(major2), rs(minor1), rs(unit_name), rs(generalize), rs(rock_type)))), tolower(trimws(paste(rs(muname), rs(compname), rs(taxsuborder), rs(taxgreatgroup), rs(reskind), ifelse(oc_c == "non-outcrop matrix", "", oc_c))))))
  ph_v <- rn(soil_ph); sand_v <- rn(soil_sand); cec_v <- rn(soil_cec); bed_v <- rn(bedrock_depth_cm); ro_v <- rn(rock_outcrop_pct); lat_v <- rn(lat); lon_v <- rn(lon)

  m_u <- grepl("ultramafic|serpentin|peridotit|dunit|pyroxenit|soapstone|talc|ophiolit|chrome|conowingo|aldino", combined) | (oc_c == "serpentine")
  m_m <- !m_u & (grepl("mafic|diabas|gabbro|basalt|metabasalt|amphibolit|greenstone|hornblende|metavolcanic|metadiorit|norite|anorthosit|trap rock|iredell|mecklenburg|wilkes|poindexter|enon", combined) | (oc_c == "diabase/mafic glade"))
  m_c <- !(m_u | m_m) & (grepl("carbonat|limeston|dolostone|dolomit|marble|chalk|marl|calcarenit|travertin|karst|coquina|alvar|gladeville|barfield|talbott|colbert|rendolls|calciustolls|calciudolls", combined) | (oc_c == "limestone glade"))
  m_f <- !(m_u | m_m | m_c) & (grepl("felsic|crystalline|granit|granodiorit|monzonit|quartz monzonit|gneiss|pegmatit|rhyolit|aplit|dacit|migmatit|wake|louisburg|ashlar|rion|wateree|pacolet|appling|cecil|wedowee|madison", combined) | (oc_c == "granitic flatrock"))
  m_s <- !(m_u | m_m | m_c | m_f) & (grepl("siliciclastic|(^|[^[:alnum:]_])consolidated($|[^[:alnum:]_])|sandstone|quartzit|arkos|conglomerat|greywack|graywack|siltston|shale|slate|phyllit|mudston|argillit|chert|novaculit|schist|ramsey|lily|sewanee|goldston|manteo|wallen|dekab|berks|weikert|channery", combined) | (oc_c == "sandstone/siliceous"))
  m_n <- !(m_u | m_m | m_c | m_f | m_s) & grepl("unconsolidated|alluvi|colluvi|coastal plain|terrace|loess|till|outwash|glacial|dune|eolian|sandhill|muck|peat|saprist|psamment|aquult|aquoll|aquept|fluvent|wagram|lakeland|plainfield|goldsboro|norfolk|myatt|rains|lynchburg|tama", combined)

  out_cls <- character(n); out_cls[m_u] <- "Ultramafic_Serpentine"; out_cls[m_m] <- "Mafic_Igneous_Metamorphic"; out_cls[m_c] <- "Carbonate_Calcareous"
  out_cls[m_f] <- "Felsic_Granitic_Gneiss"; out_cls[m_s] <- "Siliciclastic_Consolidated"; out_cls[m_n] <- "Unconsolidated_Sedimentary"

  rem <- which(!nzchar(out_cls))
  if (length(rem) > 0L) {
    is_sh <- (is.finite(bed_v[rem]) & bed_v[rem] <= 50.0) | (is.finite(ro_v[rem]) & ro_v[rem] >= 15.0) | grepl("lithic|paralithic|bedrock|rock outcrop", combined[rem])
    ph_r <- ph_v[rem]; sa_r <- sand_v[rem]; ce_r <- cec_v[rem]; bd_r <- bed_v[rem]; ro_r <- ro_v[rem]; la_r <- lat_v[rem]; lo_r <- lon_v[rem]
    demo_fb <- ifelse(!is.na(la_r) & !is.na(lo_r), ifelse(lo_r > -81.5 & la_r < 36.8 & (lo_r > -79.2 | la_r < 33.2), "Unconsolidated_Sedimentary", ifelse(lo_r >= -84.5 & lo_r <= -77.5 & la_r >= 32.5 & la_r <= 37.5, "Felsic_Granitic_Gneiss", ifelse(lo_r >= -88.5 & lo_r < -83.0 & la_r >= 34.5 & la_r <= 38.5, "Carbonate_Calcareous", ifelse(lo_r >= -84.0 & lo_r <= -78.0 & la_r > 35.0 & la_r <= 41.0, "Siliciclastic_Consolidated", "Unconsolidated_Sedimentary")))), "Unconsolidated_Sedimentary")
    out_cls[rem] <- ifelse(is_sh & is.finite(ph_r) & ph_r >= 7.0, "Carbonate_Calcareous", ifelse(is_sh & is.finite(ph_r) & ph_r >= 6.2 & is.finite(ce_r) & ce_r >= 20.0 & (!is.finite(sa_r) | sa_r < 45.0), "Mafic_Igneous_Metamorphic", ifelse(is_sh & is.finite(ph_r) & ph_r <= 5.3 & is.finite(sa_r) & sa_r >= 60.0, "Felsic_Granitic_Gneiss", ifelse(is_sh, "Siliciclastic_Consolidated", ifelse(is.finite(bd_r) & bd_r >= 150.0 & (!is.finite(ro_r) | ro_r < 5.0) & ((is.finite(sa_r) & sa_r >= 70.0) | nzchar(combined[rem])), "Unconsolidated_Sedimentary", if (isTRUE(eff_demo)) demo_fb else "Unknown_Lithology")))))
  }
  out_cls
}

extract_sgmc_lithology <- function(df, env_dir = "data/environmental", demo = FALSE) {
  df <- as.data.frame(df, stringsAsFactors = FALSE); n <- nrow(df); if (n == 0L) return(df)
  pc <- function(cs) { for (cc in cs) if (cc %in% names(df)) return(as.character(df[[cc]])); rep("", n) }
  pn <- function(cs) { for (cc in cs) if (cc %in% names(df)) return(suppressWarnings(as.numeric(df[[cc]]))); rep(NA_real_, n) }
  r1 <- pc(c("sgmc_rocktype1", "rocktype1")); r2 <- pc(c("sgmc_rocktype2", "rocktype2")); m1 <- pc(c("sgmc_major1", "major1"))
  un <- pc(c("sgmc_unit_name", "unit_name")); gn <- pc(c("sgmc_generalize", "GENERALIZE")); rt <- pc(c("sgmc_rock_type", "ROCK_TYPE"))
  prev <- if ("sgmc_lithology_class" %in% names(df)) as.character(df$sgmc_lithology_class) else rep("", n)
  prev[is.na(prev) | tolower(trimws(prev)) %in% c("na", "unknown", "unknown_lithology")] <- ""

  lat_v <- suppressWarnings(as.numeric(df$latitude)); lon_v <- suppressWarnings(as.numeric(df$longitude)); v_pts <- !is.na(lat_v) & !is.na(lon_v)
  if (!isTRUE(demo) && any(v_pts) && requireNamespace("sf", quietly = TRUE)) {
    sg_f <- unlist(lapply(unique(c(file.path(env_dir, "sgmc"), env_dir)), function(d) if (dir.exists(d)) list.files(d, pattern = "sgmc.*\\.(gpkg|shp|geojson)$", full.names = TRUE, ignore.case = TRUE) else character(0)))
    if (length(sg_f) > 0L) tryCatch({
      sf_g <- sf::st_read(sg_f[1L], quiet = TRUE)
      hits <- sf::st_intersects(sf::st_transform(sf::st_as_sf(df[v_pts, c("longitude", "latitude"), drop = FALSE], coords = c("longitude", "latitude"), crs = 4326), sf::st_crs(sf_g)), sf_g)
      vi <- which(v_pts); cu <- toupper(names(sf_g))
      pm <- function(cs) { for (cc in cs) { j <- match(cc, cu); if (!is.na(j)) return(names(sf_g)[j]) }; NA_character_ }
      c1 <- pm(c("ROCKTYPE1", "MAJOR1")); c2 <- pm(c("ROCKTYPE2", "MAJOR2", "MINOR1")); cm <- pm(c("MAJOR1", "ROCKTYPE1"))
      cu1 <- pm(c("UNIT_NAME", "UNITNAME")); cg <- pm(c("GENERALIZE", "GENERALIZED_LITH", "GEN_LITH")); cr <- pm(c("ROCK_TYPE", "ORIG_LABEL"))
      for (k in seq_along(vi)) if (length(hits[[k]]) > 0L) {
        rg <- hits[[k]][1L]
        if (!is.na(c1)) r1[vi[k]] <- as.character(sf_g[[c1]][rg]); if (!is.na(c2)) r2[vi[k]] <- as.character(sf_g[[c2]][rg])
        if (!is.na(cm)) m1[vi[k]] <- as.character(sf_g[[cm]][rg]); if (!is.na(cu1)) un[vi[k]] <- as.character(sf_g[[cu1]][rg])
        if (!is.na(cg)) gn[vi[k]] <- as.character(sf_g[[cg]][rg]); if (!is.na(cr)) rt[vi[k]] <- as.character(sf_g[[cr]][rg])
      }
    }, error = function(e) NULL)
  }

  lith <- standardize_sgmc_lithology(
    rocktype1 = r1, rocktype2 = r2, major1 = m1, unit_name = un, muname = pc("muname"), outcrop_class = pc("outcrop_class"),
    lat = lat_v, lon = lon_v, is_demo = isTRUE(demo), generalize = ifelse(nzchar(trimws(gn)), gn, prev), rock_type = rt,
    compname = pc("compname"), taxsuborder = pc("taxsuborder"), taxgreatgroup = pc("taxgreatgroup"), reskind = pc("reskind"),
    soil_ph = pn(c("fine_soil_ph", "soil_ph")), soil_sand = pn(c("fine_soil_sand", "soil_sand")), soil_cec = pn(c("fine_soil_cec", "soil_cec")),
    bedrock_depth_cm = pn(c("fine_bedrock_depth_cm", "ssurgo_bedrock_depth_cm")), rock_outcrop_pct = pn("ssurgo_rock_outcrop_pct")
  )
  df$sgmc_lithology_class <- lith
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n); is_out[is.na(is_out)] <- FALSE
  op <- pn("ssurgo_rock_outcrop_pct"); bd <- pn("fine_bedrock_depth_cm")
  df$is_azonal_lithology <- (lith == "Ultramafic_Serpentine") | ((lith %in% c("Ultramafic_Serpentine", "Mafic_Igneous_Metamorphic", "Carbonate_Calcareous", "Felsic_Granitic_Gneiss")) & (is_out | (is.finite(op) & op >= 15.0) | (is.finite(bd) & bd <= 40.0)))
  df
}

load_and_update_fine_grain_cache <- function(
  df, env_dir = "data/environmental", ssurgo_cache_dir = "data/environmental/ssurgo_cache",
  batch_size = 75, demo = FALSE, query_ssurgo_fn = NULL, cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv")
) {
  df <- flag_coordinate_uncertainty(df); n <- nrow(df); if (n == 0L) return(df)
  df$catalogNumber <- as.character(df$catalogNumber); df$lat_lon_key <- format_lat_lon_key(df$latitude, df$longitude)
  req_core <- setdiff(CORE_FINE_GRAIN_CACHE_COLS, c("catalogNumber", "lat_lon_key", "latitude", "longitude"))

  cached_df <- NULL
  sql_c <- if (!isTRUE(demo) && exists("query_environmental_cache_sqlite", mode = "function")) tryCatch(query_environmental_cache_sqlite(catalog_numbers = df$catalogNumber, lat_lon_keys = df$lat_lon_key, cache_path = cache_path, demo = demo), error = function(e) NULL) else NULL
  if (!isTRUE(demo) && file.exists(cache_path)) {
    cached_df <- tryCatch({
      cr <- if (requireNamespace("data.table", quietly = TRUE)) as.data.frame(data.table::fread(cache_path, showProgress = FALSE, na.strings = c("", "NA")), stringsAsFactors = FALSE) else as.data.frame(readr::read_csv(cache_path, show_col_types = FALSE), stringsAsFactors = FALSE)
      cr$catalogNumber <- as.character(cr$catalogNumber); if (!("lat_lon_key" %in% names(cr))) cr$lat_lon_key <- format_lat_lon_key(cr$latitude, cr$longitude); cr
    }, error = function(e) NULL)
  }
  if (!is.null(sql_c) && nrow(sql_c) > 0L) {
    sql_c$catalogNumber <- as.character(sql_c$catalogNumber); if (!("lat_lon_key" %in% names(sql_c))) sql_c$lat_lon_key <- format_lat_lon_key(sql_c$latitude, sql_c$longitude)
    cached_df <- if (!is.null(cached_df) && nrow(cached_df) > 0L) {
      if (exists(".bind_rows_schema_safe", mode = "function")) .bind_rows_schema_safe(sql_c, cached_df) else {
        for (.cn in intersect(names(sql_c), names(cached_df))) if (!identical(typeof(sql_c[[.cn]]), typeof(cached_df[[.cn]]))) {
          if (is.logical(cached_df[[.cn]]) && is.numeric(sql_c[[.cn]])) { sql_c[[.cn]] <- ifelse(is.na(sql_c[[.cn]]), NA, sql_c[[.cn]] != 0) } else if (is.numeric(sql_c[[.cn]]) && is.numeric(cached_df[[.cn]])) { sql_c[[.cn]] <- as.numeric(sql_c[[.cn]]); cached_df[[.cn]] <- as.numeric(cached_df[[.cn]]) } else { sql_c[[.cn]] <- as.character(sql_c[[.cn]]); cached_df[[.cn]] <- as.character(cached_df[[.cn]]) }
        }
        dplyr::bind_rows(sql_c, cached_df)
      }
    } else sql_c
  }

  pol_d <- .list_real_fine_rasters(file.path(env_dir, "gssurgo_30m")); sol_d <- c(.list_real_fine_rasters(file.path(env_dir, "solus_100m")), pol_d)
  has_pol <- length(pol_d) > 0L; has_sol <- length(sol_d) > 0L; has_sol_ph <- any(grepl("^ph1to1h2o|^ph_", basename(sol_d), ignore.case = TRUE))
  stale_dem <- file.path(env_dir, "dem_3dep_30m", "dem_3dep_regional.tif"); if (!isTRUE(demo) && file.exists(stale_dem)) tryCatch(unlink(stale_dem), error = function(e) NULL)
  has_dem <- length(.list_real_fine_rasters(file.path(env_dir, "dem_3dep_30m"), exclude_pattern = "^dem_3dep_regional\\.tif$")) > 0L

  if (!is.null(cached_df) && nrow(cached_df) > 0L && any(c("fine_soil_ph", "polaris_ph_30m", "solus_resdept_cm") %in% names(cached_df))) {
    avail_c <- intersect(setdiff(FINE_GRAIN_CACHE_COLS, c("catalogNumber", "lat_lon_key", "latitude", "longitude", "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord", "is_clim_eligible_coord", "is_edaphic_eligible_coord", "edaphic_coord_weight")), names(cached_df))
    c_cat <- cached_df[!duplicated(cached_df$catalogNumber), c("catalogNumber", avail_c), drop = FALSE]
    m_i <- match(df$catalogNumber, c_cat$catalogNumber); was_in_cache <- !is.na(m_i); miss_c <- which(is.na(m_i) & !is.na(df$lat_lon_key))
    if (length(miss_c) > 0L) {
      c_ll <- cached_df[!is.na(cached_df$lat_lon_key) & !duplicated(cached_df$lat_lon_key), c("lat_lon_key", avail_c), drop = FALSE]
      ll_i <- match(df$lat_lon_key[miss_c], c_ll$lat_lon_key); was_in_cache[miss_c[!is.na(ll_i)]] <- TRUE
      for (col in avail_c) { df[[col]] <- c_cat[[col]][m_i]; if (any(!is.na(ll_i))) df[[col]][miss_c[!is.na(ll_i)]] <- c_ll[[col]][ll_i[!is.na(ll_i)]] }
    } else for (col in avail_c) df[[col]] <- c_cat[[col]][m_i]
    df <- mask_micro_scale_predictors_by_uncertainty(df, cols = c("topo_tpi_30m", "topo_hli_30m", "topo_twi_30m"), max_uncertainty_m = 250.0)

    lat_c <- suppressWarnings(as.numeric(df$latitude)); lon_c <- suppressWarnings(as.numeric(df$longitude))
    in_conus_c <- is.finite(lat_c) & lat_c >= 24.0 & lat_c <= 49.5 & is.finite(lon_c) & lon_c >= -125.0 & lon_c <= -66.0
    has_elev_col <- "elevation_m" %in% names(cached_df)
    need_dem_up <- isTRUE(has_dem) && !isTRUE(has_elev_col) && any(was_in_cache & in_conus_c)
    uncached_mask <- (!was_in_cache) | (was_in_cache & in_conus_c & need_dem_up)
    n_unc <- sum(uncached_mask)
    message(sprintf("Loaded %d / %d fine-grain edaphic records from %s (%d coordinates requiring refresh).", n - n_unc, n, cache_path, n_unc))

    if (n_unc > 0L) {
      sub_d <- df[uncached_mask, , drop = FALSE]; for (col in req_core) sub_d[[col]] <- NULL
      if (!is.null(query_ssurgo_fn)) sub_d <- query_ssurgo_fn(sub_d, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
      sub_d <- extract_uncertainty_buffered_predictors(extract_micro_topography_3dep(extract_polaris_and_solus_predictors(enrich_ssurgo_pedology_from_series(sub_d, demo = demo), env_dir = env_dir, demo = demo), env_dir = env_dir, demo = demo), env_dir = env_dir, demo = demo)
      for (col in setdiff(names(sub_d), c("catalogNumber", "lat_lon_key", "latitude", "longitude"))) { if (!(col %in% names(df))) df[[col]] <- rep(NA, n); df[[col]][uncached_mask] <- sub_d[[col]] }
    }

    has_sgmc <- "sgmc_lithology_class" %in% names(df)
    sg_emp <- if (has_sgmc) is.na(df$sgmc_lithology_class) | !nzchar(trimws(as.character(df$sgmc_lithology_class))) else rep(TRUE, n)
    need_sgmc <- sg_emp & ((!was_in_cache) | (!has_sgmc & in_conus_c))
    if (any(need_sgmc)) {
      g_sub <- extract_sgmc_lithology(df[need_sgmc, , drop = FALSE], env_dir = env_dir, demo = demo)
      for (col in c("sgmc_lithology_class", "is_azonal_lithology")) if (col %in% names(g_sub)) { if (!(col %in% names(df))) df[[col]] <- rep(NA, n); df[[col]][need_sgmc] <- g_sub[[col]] }
    }

    san_st <- FALSE
    if (!isTRUE(demo)) {
      pc_c <- intersect(c("polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m"), names(df))
      if (any(vapply(pc_c, function(cc) any(!is.na(df[[cc]])), logical(1)))) { df <- sanitize_disabled_polaris_columns(df); san_st <- TRUE }
      if (!has_pol && !has_sol_ph) { pb <- intersect(c("ph_q10_buffer", "ph_q90_buffer"), names(df)); if (any(vapply(pb, function(cc) any(!is.na(df[[cc]])), logical(1)))) { for (cc in pb) df[[cc]] <- rep(NA_real_, n); san_st <- TRUE } }
      if (!has_sol) { sc <- intersect(c("solus_resdept_cm", "solus_fragvol_pct", "solus_sand_pct", "solus_clay_pct", "solus_silt_pct", "solus_cec", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct"), names(df)); if (any(vapply(sc, function(cc) any(!is.na(df[[cc]])), logical(1)))) { for (cc in sc) df[[cc]] <- rep(NA_real_, n); san_st <- TRUE } }
      if (!has_dem) { dc <- intersect(c("topo_tpi_30m", "topo_hli_30m", "topo_twi_30m"), names(df)); if (any(vapply(dc, function(cc) any(!is.na(df[[cc]])), logical(1)))) { for (cc in dc) df[[cc]] <- rep(NA_real_, n); df <- extract_micro_topography_3dep(df, env_dir = env_dir, demo = FALSE); san_st <- TRUE } }
      if (isTRUE(san_st)) df <- extract_sgmc_lithology(extract_polaris_and_solus_predictors(df, env_dir = env_dir, demo = FALSE), env_dir = env_dir, demo = FALSE)
    }

    if (all(c("fine_soil_sand", "fine_soil_clay") %in% names(df))) {
      m_ilr <- !("ilr_texture_1" %in% names(df)) | !("ilr_texture_2" %in% names(df)) | (is.na(df$ilr_texture_1) & is.finite(suppressWarnings(as.numeric(df$fine_soil_sand))) & is.finite(suppressWarnings(as.numeric(df$fine_soil_clay))))
      if (any(m_ilr)) {
        th <- compute_soil_texture_log_ratios(df$fine_soil_sand, df$fine_soil_clay, silt_pct = if ("soil_silt_pct" %in% names(df)) df$soil_silt_pct else NULL)
        for (tc in c("soil_silt_pct", "ilr_texture_1", "ilr_texture_2")) { if (!(tc %in% names(df))) df[[tc]] <- th[[tc]] else df[[tc]][m_ilr] <- th[[tc]][m_ilr] }
      }
    }

    unc_c <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters)); has_bc <- "min_bedrock_depth_buffer_cm" %in% names(df)
    b_el <- is.finite(lat_c) & is.finite(lon_c) & in_conus_c & (df$is_edaphic_eligible_coord | is.na(unc_c) | (is.finite(unc_c) & unc_c <= 5000.0))
    need_buf <- !isTRUE(demo) && has_sol && (!has_bc || any((!was_in_cache) & b_el & is.na(df$min_bedrock_depth_buffer_cm)))
    if (isTRUE(need_buf)) {
      bm <- b_el & (if (has_bc) ((!was_in_cache) & is.na(df$min_bedrock_depth_buffer_cm)) else rep(TRUE, n))
      if (any(bm)) {
        bs <- extract_uncertainty_buffered_predictors(df[bm, , drop = FALSE], env_dir = env_dir, demo = demo)
        for (bc in c("ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer")) if (bc %in% names(bs)) { if (!(bc %in% names(df))) df[[bc]] <- rep(NA, n); df[[bc]][bm] <- bs[[bc]] }
      }
    }
    df$is_rock_outcrop <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n); df$is_rock_outcrop[is.na(df$is_rock_outcrop)] <- FALSE
    if (!san_st && n_unc == 0L && !any(need_sgmc) && !isTRUE(need_buf)) return(df)
  } else {
    if (!is.null(query_ssurgo_fn)) df <- query_ssurgo_fn(df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
    df <- extract_uncertainty_buffered_predictors(extract_sgmc_lithology(extract_micro_topography_3dep(extract_polaris_and_solus_predictors(enrich_ssurgo_pedology_from_series(df, demo = demo), env_dir = env_dir, demo = demo), env_dir = env_dir, demo = demo), env_dir = env_dir, demo = demo), env_dir = env_dir, demo = demo)
  }

  if (!isTRUE(demo)) {
    ec <- intersect(FINE_GRAIN_CACHE_COLS, names(df)); nc <- df[, ec, drop = FALSE]
    comb <- if (!is.null(cached_df) && nrow(cached_df) > 0L) {
      for (mc in setdiff(ec, names(cached_df))) cached_df[[mc]] <- NA
      cc <- intersect(ec, names(cached_df)); rb <- rbind(nc[, cc, drop = FALSE], cached_df[, cc, drop = FALSE]); rb[!duplicated(rb$catalogNumber), , drop = FALSE]
    } else nc[!duplicated(nc$catalogNumber), , drop = FALSE]
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE); tc <- paste0(cache_path, ".tmp")
    if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(comb, file = tc, na = "", quote = "auto", logical01 = FALSE) else readr::write_csv(comb, tc)
    file.rename(tc, cache_path)
    if (exists("upsert_environmental_cache_sqlite", mode = "function")) tryCatch(upsert_environmental_cache_sqlite(comb, cache_path = cache_path, demo = demo), error = function(e) NULL)
    message(sprintf("Persisted %d fine-grain edaphic & micro-topographic records to %s.", nrow(comb), cache_path))
  }
  df
}

extract_fine_grain_edaphic_metrics <- function(df, env_dir = "data/environmental", ssurgo_cache_dir = "data/environmental/ssurgo_cache", batch_size = 75, demo = FALSE, query_ssurgo_fn = NULL, cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv")) {
  load_and_update_fine_grain_cache(df = df, env_dir = env_dir, ssurgo_cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo, query_ssurgo_fn = query_ssurgo_fn, cache_path = cache_path)
}
