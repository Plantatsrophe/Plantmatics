# ==============================================================================
# Module: scripts/analysis/step05_env_raster_extract_utils.R
# Purpose: WorldClim v2.1 & SoilGrids 250 m buffered extraction, North American
#          raster cropping, SQLite cache-first hydration, and cross-border
#          candidate ecological predictor discovery.
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists("WC_BIO_LOOKUP", inherits = TRUE)) {
  WC_BIO_LOOKUP <- c(
    "1" = "bio1_temp_mean", "2" = "bio2_mean_diurnal_range", "3" = "bio3_isothermality", "4" = "bio4_temp_seasonality",
    "5" = "bio5_max_temp_warmest_month", "6" = "bio6_min_temp_coldest_month", "7" = "bio7_temp_annual_range", "8" = "bio8_mean_temp_wettest_qtr",
    "9" = "bio9_mean_temp_driest_qtr", "10" = "bio10_mean_temp_warmest_qtr", "11" = "bio11_mean_temp_coldest_qtr", "12" = "bio12_precip_annual",
    "13" = "bio13_precip_wettest_month", "14" = "bio14_precip_driest_month", "15" = "bio15_precip_seasonality", "16" = "bio16_precip_wettest_qtr",
    "17" = "bio17_precip_driest_qtr", "18" = "bio18_precip_warmest_qtr", "19" = "bio19_precip_coldest_qtr"
  )
}
VALID_CLIM_BIO_PREDICTORS <- c("bio1_temp_mean", "bio2_mean_diurnal_range", "bio4_temp_seasonality", "bio5_max_temp_warmest_month", "bio6_min_temp_coldest_month", "bio10_mean_temp_warmest_qtr", "bio11_mean_temp_coldest_qtr", "bio12_precip_annual", "bio13_precip_wettest_month", "bio14_precip_driest_month", "bio15_precip_seasonality", "bio16_precip_wettest_qtr", "bio17_precip_driest_qtr")
VALID_TOPO_PREDICTORS <- c("topo_slope_deg", "topo_tpi_30m", "topo_hli_30m", "topo_twi_30m")
NON_ENM_RAW_GEOCHEM_HYDRO_COLS <- c("soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct", "ultramafic_calcareous_gradient", "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m", "nwpl_indicator_code", "nwpl_wetland_score")

extract_rast_buffered_mean <- function(rast_obj, pts_vect, radii_m) {
  r_num <- suppressWarnings(as.numeric(radii_m))
  r_num[is.na(r_num) | !is.finite(r_num)] <- 250.0; r_num <- pmax(100.0, pmin(1000.0, r_num))
  crds <- tryCatch(terra::crds(pts_vect), error = function(e) NULL)
  inv_idx <- NULL; pts_use <- pts_vect; r_use <- r_num
  if (!is.null(crds) && is.matrix(crds) && nrow(crds) == length(r_num) && nrow(crds) > 1L) {
    u_key <- paste(round(crds[, 1L], 6L), round(crds[, 2L], 6L), round(r_num, 1L), sep = "_")
    u_keys <- unique(u_key)
    if (length(u_keys) < nrow(crds)) {
      first_idx <- match(u_keys, u_key); inv_idx <- match(u_key, u_keys)
      pts_use <- pts_vect[first_idx]; r_use <- r_num[first_idx]
    }
  }
  pt_vals <- terra::extract(rast_obj, pts_use)
  buf_vals <- tryCatch(terra::extract(rast_obj, terra::buffer(pts_use, width = r_use), fun = mean, na.rm = TRUE), error = function(e) NULL)
  if (!is.null(buf_vals)) {
    for (col_nm in intersect(setdiff(names(pt_vals), "ID"), names(buf_vals))) {
      bv <- suppressWarnings(as.numeric(buf_vals[[col_nm]])); pv <- suppressWarnings(as.numeric(pt_vals[[col_nm]]))
      pt_vals[[col_nm]] <- ifelse(is.finite(bv), bv, pv)
    }
  }
  if (!is.null(inv_idx)) { pt_vals <- pt_vals[inv_idx, , drop = FALSE]; pt_vals$ID <- seq_along(inv_idx); rownames(pt_vals) <- NULL }
  pt_vals
}

extract_soilgrids_rhizosphere_col <- function(sg_vals, var_pattern) {
  nms <- names(sg_vals)
  i0 <- grep(paste0(var_pattern, ".*(0-5|0_5)"), nms, ignore.case = TRUE)
  i5 <- grep(paste0(var_pattern, ".*(5-15|5_15)"), nms, ignore.case = TRUE)
  if (length(i0) > 0L && length(i5) > 0L) {
    v0 <- suppressWarnings(as.numeric(sg_vals[[i0[1L]]])); v5 <- suppressWarnings(as.numeric(sg_vals[[i5[1L]]]))
    if (exists("compute_rhizosphere_0_15cm", mode = "function")) return(compute_rhizosphere_0_15cm(v0, v5))
    return((5.0 * v0 + 10.0 * v5) / 15.0)
  }
  ia <- grep(var_pattern, nms, ignore.case = TRUE)
  if (length(ia) > 0L) suppressWarnings(as.numeric(sg_vals[[ia[1L]]])) else rep(NA_real_, nrow(sg_vals))
}

.ensure_wc_bio_19_cols <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  for (b_col in unname(WC_BIO_LOOKUP)) df[[b_col]] <- if (b_col %in% names(df)) suppressWarnings(as.numeric(df[[b_col]])) else rep(NA_real_, nrow(df))
  df
}

assign_wc_bio_from_extracted <- function(df, wc_vals) {
  if (is.null(wc_vals) || !is.data.frame(wc_vals)) return(.ensure_wc_bio_19_cols(df))
  val_nms <- setdiff(names(wc_vals), "ID")
  for (b_idx in seq_len(19L)) {
    tc <- WC_BIO_LOOKUP[[as.character(b_idx)]]; hit <- character(0)
    if (tc %in% val_nms) hit <- tc else {
      m1 <- grep(sprintf("(^|[._-])bio[_.-]?0?%d$", b_idx), val_nms, value = TRUE, ignore.case = TRUE)
      if (length(m1) > 0L) hit <- m1[1L] else {
        m2 <- grep(sprintf("^bio0?%d(_|$)", b_idx), val_nms, value = TRUE, ignore.case = TRUE)
        if (length(m2) > 0L) hit <- m2[1L]
      }
    }
    if (length(hit) == 0L && length(val_nms) == 19L && all(grepl("bio|lyr|band", val_nms, ignore.case = TRUE))) hit <- val_nms[b_idx]
    if (length(hit) > 0L) {
      ve <- suppressWarnings(as.numeric(wc_vals[[hit[1L]]]))
      df[[tc]] <- if (tc %in% names(df)) ifelse(is.finite(ve), ve, suppressWarnings(as.numeric(df[[tc]]))) else ve
    }
  }
  .ensure_wc_bio_19_cols(df)
}

discover_candidate_eco_predictors <- function(df, include_pheno = FALSE, min_finite = 3L, min_coverage_frac = 0.50) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(character(0))
  df_g <- if (exists("mask_micro_scale_predictors_by_uncertainty", mode = "function")) mask_micro_scale_predictors_by_uncertainty(df, max_uncertainty_m = 250.0) else df
  nms <- names(df_g); n_r <- nrow(df_g)
  min_obs <- if (n_r >= 8L) max(min(as.integer(min_finite), max(2L, n_r)), as.integer(ceiling(min_coverage_frac * n_r))) else min(as.integer(min_finite), max(2L, n_r))
  ok_c <- function(cl) {
    if (!(cl %in% nms)) return(FALSE)
    rc <- df_g[[cl]]; if (!(is.numeric(rc) || is.integer(rc) || is.logical(rc))) return(FALSE)
    v <- suppressWarnings(as.numeric(rc)); f <- is.finite(v); sum(f) >= min_obs && is.finite(stats::sd(v[f])) && stats::sd(v[f]) > 1e-8
  }
  drop_conus <- character(0); lat_g <- rep(NA_real_, n_r); lon_g <- rep(NA_real_, n_r)
  for (lc in intersect(c("latitude", "decimalLatitude"), nms)) { vl <- suppressWarnings(as.numeric(df_g[[lc]])); m <- !is.finite(lat_g) & is.finite(vl); if (any(m)) lat_g[m] <- vl[m] }
  for (lc in intersect(c("longitude", "decimalLongitude"), nms)) { vl <- suppressWarnings(as.numeric(df_g[[lc]])); m <- !is.finite(lon_g) & is.finite(vl); if (any(m)) lon_g[m] <- vl[m] }
  geo_ok <- is.finite(lat_g) & is.finite(lon_g) & (abs(lat_g) <= 90.0) & (abs(lon_g) <= 180.0)
  out_c <- geo_ok & ((lat_g > 49.5) | (lat_g < 24.5) | (lon_g < -125.0))
  if (sum(geo_ok) > 0L && (sum(out_c) / sum(geo_ok)) > 0.20) {
    for (cc in c("fine_bedrock_depth_cm", "lithic_outcrop_index")) {
      if ((cc %in% nms) && !any(is.finite(suppressWarnings(as.numeric(df_g[[cc]][out_c]))))) drop_conus <- c(drop_conus, cc)
    }
  }
  ed <- character(0)
  if (ok_c("fine_soil_ph")) ed <- c(ed, "fine_soil_ph") else if (ok_c("soil_ph")) ed <- c(ed, "soil_ph")
  if (ok_c("ilr_texture_1") && ok_c("ilr_texture_2")) ed <- c(ed, "ilr_texture_1", "ilr_texture_2") else if (ok_c("ilr_texture_1")) ed <- c(ed, "ilr_texture_1") else {
    if (ok_c("fine_soil_sand")) ed <- c(ed, "fine_soil_sand") else if (ok_c("soil_sand")) ed <- c(ed, "soil_sand")
    if (ok_c("fine_soil_clay")) ed <- c(ed, "fine_soil_clay") else if (ok_c("soil_clay")) ed <- c(ed, "soil_clay")
  }
  if (ok_c("fine_soil_cec")) ed <- c(ed, "fine_soil_cec") else if (ok_c("soil_cec")) ed <- c(ed, "soil_cec")
  if (ok_c("fine_cfvo_pct")) ed <- c(ed, "fine_cfvo_pct") else if (ok_c("soil_cfvo")) ed <- c(ed, "soil_cfvo")
  for (ex in c("fine_bedrock_depth_cm", "lithic_outcrop_index", "wetland_seepage_index")) if (!(ex %in% drop_conus) && ok_c(ex)) ed <- c(ed, ex)
  topo <- VALID_TOPO_PREDICTORS[vapply(VALID_TOPO_PREDICTORS, ok_c, logical(1))]
  clim <- c(VALID_CLIM_BIO_PREDICTORS, "summer_moisture_deficit_idx", "maritime_continentality_idx")
  clim <- clim[vapply(clim, ok_c, logical(1))]
  ph <- if (isTRUE(include_pheno) && ok_c("delta_doy")) "delta_doy" else character(0)
  out <- setdiff(unique(c(ed, topo, clim, ph)), NON_ENM_RAW_GEOCHEM_HYDRO_COLS)
  out <- out[!grepl("^polaris_", out, ignore.case = TRUE)]
  if (length(out) == 0L && min_coverage_frac > 0.0) return(discover_candidate_eco_predictors(df, include_pheno, min_finite, 0.0))
  out
}

crop_and_compress_na_rasters <- function(rast_obj = NULL, env_dir = NULL, na_bbox = c(-130.0, -50.0, 14.0, 66.0)) {
  if (!requireNamespace("terra", quietly = TRUE)) return(rast_obj)
  na_ext <- terra::ext(na_bbox[1L], na_bbox[2L], na_bbox[3L], na_bbox[4L])
  .crop_one <- function(tp) {
    if (is.null(tp) || !nzchar(tp) || !file.exists(tp) || grepl("\\.crop_tmp\\.tif$", tp, ignore.case = TRUE)) return(tp)
    rd <- tryCatch(terra::rast(tp), error = function(e) NULL)
    ex <- if (!is.null(rd)) tryCatch(terra::ext(rd), error = function(e) NULL) else NULL
    if (is.null(ex) || !(terra::xmin(ex) < -150.0 || terra::xmax(ex) > 0.0 || terra::ymin(ex) < 0.0 || terra::ymax(ex) > 75.0)) return(tp)
    tmp <- paste0(tp, ".crop_tmp.tif")
    ok <- tryCatch({
      rc <- terra::crop(rd, na_ext)
      tryCatch(terra::writeRaster(rc, filename = tmp, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")),
               error = function(e2) terra::writeRaster(rc, filename = tmp, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6")))
      TRUE
    }, error = function(e) FALSE)
    rm(rd)
    if (isTRUE(ok) && file.exists(tmp)) { unlink(tp); file.rename(tmp, tp) } else if (file.exists(tmp)) unlink(tmp)
    tp
  }
  cf <- character(0)
  if (!is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir)) cf <- c(cf, list.files(env_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE))
  if (is.character(rast_obj)) {
    cf <- c(cf, rast_obj[file.exists(rast_obj)])
  } else if (inherits(rast_obj, "SpatRaster")) {
    cf <- c(cf, tryCatch(terra::sources(rast_obj), error = function(e) character(0)))
  } else if (is.list(rast_obj)) {
    for (el in rast_obj) if (inherits(el, "SpatRaster")) cf <- c(cf, tryCatch(terra::sources(el), error = function(e) character(0)))
  }
  for (p in unique(cf[nzchar(cf) & file.exists(cf)])) tryCatch(.crop_one(p), error = function(e) NULL)
  if (inherits(rast_obj, "SpatRaster")) {
    srcs <- tryCatch(terra::sources(rast_obj), error = function(e) character(0)); srcs <- srcs[nzchar(srcs) & file.exists(srcs)]
    if (length(srcs) > 0L && all(file.exists(srcs))) rast_obj <- tryCatch(terra::rast(srcs), error = function(e) terra::crop(rast_obj, na_ext))
    ex_c <- tryCatch(terra::ext(rast_obj), error = function(e) NULL)
    if (!is.null(ex_c) && (terra::xmin(ex_c) < -135.0 || terra::xmax(ex_c) > -45.0 || terra::ymin(ex_c) < 10.0 || terra::ymax(ex_c) > 70.0)) rast_obj <- tryCatch(terra::crop(rast_obj, na_ext), error = function(e) rast_obj)
    return(rast_obj)
  } else if (is.list(rast_obj)) return(lapply(rast_obj, function(el) if (inherits(el, "SpatRaster")) crop_and_compress_na_rasters(el, NULL, na_bbox) else el))
  rast_obj
}

discover_env_rasters <- function(env_dir, recursive = TRUE, crop_global = TRUE) {
  if (is.null(env_dir) || !nzchar(env_dir) || !dir.exists(env_dir)) return(character(0))
  if (isTRUE(crop_global) && requireNamespace("terra", quietly = TRUE)) tryCatch(crop_and_compress_na_rasters(env_dir = env_dir), error = function(e) NULL)
  f_vec <- list.files(env_dir, pattern = "\\.(tif|tiff|grd)$", full.names = TRUE, recursive = isTRUE(recursive), ignore.case = TRUE)
  f_vec[!grepl("\\.crop_tmp\\.tif$", f_vec, ignore.case = TRUE)]
}

.hydrate_env_from_sqlite <- function(df, env_dir, target_cols, demo = FALSE) {
  if (isTRUE(demo) || nrow(df) == 0L || !exists("query_environmental_cache_sqlite", mode = "function")) return(list(df = df, hit = rep(FALSE, nrow(df))))
  ll <- ifelse(is.finite(suppressWarnings(as.numeric(df$latitude))) & is.finite(suppressWarnings(as.numeric(df$longitude))),
               sprintf("%.5f_%.5f", round(as.numeric(df$latitude), 5), round(as.numeric(df$longitude), 5)), NA_character_)
  sq <- tryCatch(query_environmental_cache_sqlite(as.character(df$catalogNumber), ll, file.path(env_dir, "fine_grain_edaphic_cache.csv"), demo), error = function(e) NULL)
  if (is.null(sq) || !is.data.frame(sq) || nrow(sq) == 0L) return(list(df = df, hit = rep(FALSE, nrow(df))))
  mc <- if ("catalogNumber" %in% names(sq)) match(as.character(df$catalogNumber), as.character(sq$catalogNumber)) else rep(NA_integer_, nrow(df))
  ml <- if ("lat_lon_key" %in% names(sq)) match(ll, as.character(sq$lat_lon_key)) else rep(NA_integer_, nrow(df))
  mh <- ifelse(!is.na(mc), mc, ml)
  for (hc in intersect(target_cols, names(sq))) {
    vs <- suppressWarnings(as.numeric(sq[[hc]][mh]))
    df[[hc]] <- if (hc %in% names(df)) ifelse(is.finite(suppressWarnings(as.numeric(df[[hc]]))), suppressWarnings(as.numeric(df[[hc]])), vs) else vs
  }
  list(df = .ensure_wc_bio_19_cols(df), hit = !is.na(mh))
}

extract_environmental_layers <- function(df, env_dir, demo = FALSE) {
  message("Integrating SoilGrids 250m (0-15cm Rhizosphere) & Uncertainty-Buffered WorldClim v2.1 layers...")
  for (pr_xy in list(c("latitude", "decimalLatitude"), c("longitude", "decimalLongitude"))) {
    if (!(pr_xy[1L] %in% names(df)) && pr_xy[2L] %in% names(df)) {
      df[[pr_xy[1L]]] <- suppressWarnings(as.numeric(df[[pr_xy[2L]]]))
    } else if (pr_xy[2L] %in% names(df)) {
      c1 <- suppressWarnings(as.numeric(df[[pr_xy[1L]]])); c2 <- suppressWarnings(as.numeric(df[[pr_xy[2L]]])); df[[pr_xy[1L]]] <- ifelse(is.finite(c1), c1, c2)
    }
  }
  if (!isTRUE(demo) && !is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir) && requireNamespace("terra", quietly = TRUE)) tryCatch(crop_and_compress_na_rasters(env_dir = env_dir), error = function(e) NULL)
  if (exists("flag_coordinate_uncertainty", mode = "function")) df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0) else {
    if (!("coordinateUncertaintyInMeters" %in% names(df))) df$coordinateUncertaintyInMeters <- NA_real_
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    df$high_coord_uncertainty_flag <- !is.na(u_m) & u_m > 10000.0; df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
    df$buffer_radius_m <- ifelse(is.na(u_m) | !is.finite(u_m) | u_m <= 0, 250.0, pmax(100.0, pmin(1000.0, u_m)))
  }
  .nsg <- function(v, pt) if (exists(".normalize_soilgrids_prop", mode = "function")) .normalize_soilgrids_prop(v, pt) else suppressWarnings(as.numeric(v))
  sg_spec <- list(soil_ph = c("phh2o", "ph"), soil_cec = c("cec", "cec"), soil_sand = c("sand", "pct"), soil_clay = c("clay", "pct"), soil_silt = c("silt", "pct"), soil_cfvo = c("cfvo", "pct"), soil_bulk_density = c("bdod", "bd"))
  .assign_sg <- function(d_in, sv, only_missing = FALSE) {
    for (nm in names(sg_spec)) {
      if (!only_missing || (!(nm %in% names(d_in)) && any(grepl(sg_spec[[nm]][1L], names(sv), ignore.case = TRUE)))) {
        d_in[[nm]] <- .nsg(extract_soilgrids_rhizosphere_col(sv, sg_spec[[nm]][1L]), sg_spec[[nm]][2L])
      }
    }
    d_in
  }

  req_cols <- c(unname(WC_BIO_LOOKUP), names(sg_spec))
  full_df_ref <- NULL; uncached_idx <- seq_len(nrow(df))
  is_canon_env <- !isTRUE(demo) && nrow(df) > 0L && is_canonical_prod_path(env_dir) && exists("query_environmental_cache_sqlite", mode = "function")
  if (is_canon_env) {
    hyd <- .hydrate_env_from_sqlite(df, env_dir, req_cols, demo = demo); df <- hyd$df; was_in_cache <- hyd$hit
    for (rc in setdiff(req_cols, names(df))) df[[rc]] <- rep(NA_real_, nrow(df))
    wc_c <- intersect(unname(WC_BIO_LOOKUP), names(df))
    fwc <- vapply(wc_c, function(cl) is.finite(suppressWarnings(as.numeric(df[[cl]]))), logical(nrow(df)))
    if (!is.matrix(fwc)) fwc <- matrix(fwc, nrow = nrow(df))
    has_wc <- rowSums(fwc) > 0L
    fmat <- vapply(req_cols, function(cl) is.finite(suppressWarnings(as.numeric(df[[cl]]))), logical(nrow(df)))
    if (!is.matrix(fmat)) fmat <- matrix(fmat, nrow = nrow(df))
    vll <- is.finite(suppressWarnings(as.numeric(df$latitude))) & is.finite(suppressWarnings(as.numeric(df$longitude)))
    unc_m <- vll & ((!was_in_cache & (rowSums(!fmat) > 0L)) | (was_in_cache & !has_wc))
    if (!any(unc_m) && any(vll)) { message(sprintf("Hydrated all %d records from SQLite environmental_cache (0 uncached).", sum(vll))); return(df) }
    if (sum(unc_m) > 0L && sum(unc_m) < nrow(df)) { full_df_ref <- df; uncached_idx <- which(unc_m); df <- df[uncached_idx, , drop = FALSE] }
  }
  .fin_ext <- function(sub_d) {
    out <- if (!is.null(full_df_ref)) {
      o <- full_df_ref; for (cl in setdiff(names(sub_d), c("catalogNumber", "latitude", "longitude"))) { if (!(cl %in% names(o))) o[[cl]] <- rep(NA, nrow(o)); o[[cl]][uncached_idx] <- sub_d[[cl]] }; o
    } else sub_d
    if (is_canon_env && exists("upsert_environmental_cache_sqlite", mode = "function")) {
      tryCatch({
        uc <- intersect(c("catalogNumber", "latitude", "longitude", req_cols), names(sub_d))
        if (all(c("catalogNumber", "latitude", "longitude") %in% uc)) {
          ud <- sub_d[, uc, drop = FALSE]; la_u <- suppressWarnings(as.numeric(ud$latitude)); lo_u <- suppressWarnings(as.numeric(ud$longitude))
          ud$lat_lon_key <- ifelse(is.finite(la_u) & is.finite(lo_u), sprintf("%.5f_%.5f", round(la_u, 5), round(lo_u, 5)), NA_character_)
          upsert_environmental_cache_sqlite(ud, cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv"), demo = demo)
        }
      }, error = function(e) NULL)
    }
    out
  }

  rfiles <- discover_env_rasters(env_dir, recursive = FALSE, crop_global = TRUE); extract_err <- NULL
  if (length(rfiles) > 0L && requireNamespace("terra", quietly = TRUE)) {
    ex_df <- tryCatch({
      st <- crop_and_compress_na_rasters(terra::rast(crop_and_compress_na_rasters(rfiles, env_dir = env_dir)), env_dir = env_dir)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      ex <- extract_rast_buffered_mean(st, pts, df$buffer_radius_m)
      for (v in setdiff(names(ex), "ID")) if (!(v %in% unname(WC_BIO_LOOKUP))) df[[v]] <- ex[[v]]
      df <- assign_wc_bio_from_extracted(df, ex)
      for (bv in names(sg_spec)) {
        c05 <- paste0(bv, "_0_5cm"); c515 <- paste0(bv, "_5_15cm")
        if (all(c(c05, c515) %in% names(df))) df[[bv]] <- if (exists("compute_rhizosphere_0_15cm", mode = "function")) compute_rhizosphere_0_15cm(df[[c05]], df[[c515]]) else (5.0 * as.numeric(df[[c05]]) + 10.0 * as.numeric(df[[c515]])) / 15.0
      }
      .assign_sg(df, ex, only_missing = TRUE)
    }, error = function(e) { extract_err <<- e$message; NULL })
    if (!is.null(ex_df)) return(.fin_ext(ex_df))
  }

  sub_f <- discover_env_rasters(env_dir, recursive = TRUE, crop_global = TRUE)
  wcf <- sub_f[grepl("bio[_.-]?[0-9]+", basename(sub_f), ignore.case = TRUE)]
  sgf <- sub_f[grepl("phh2o|cec|sand|clay|silt|cfvo|bdod", basename(sub_f), ignore.case = TRUE) & !grepl("solus_100m|polaris_30m", sub_f, ignore.case = TRUE)]
  if (length(wcf) > 0L && length(sgf) > 0L && requireNamespace("terra", quietly = TRUE)) {
    ex_c <- tryCatch({
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      wcv <- extract_rast_buffered_mean(crop_and_compress_na_rasters(terra::rast(crop_and_compress_na_rasters(wcf, env_dir = env_dir)), env_dir = env_dir), pts, df$buffer_radius_m)
      sgv <- extract_rast_buffered_mean(crop_and_compress_na_rasters(terra::rast(crop_and_compress_na_rasters(sgf, env_dir = env_dir)), env_dir = env_dir), pts, df$buffer_radius_m)
      .assign_sg(assign_wc_bio_from_extracted(df, wcv), sgv, FALSE)
    }, error = function(e) { extract_err <<- e$message; NULL })
    if (!is.null(ex_c)) return(.fin_ext(ex_c))
  }

  if (!isTRUE(demo) && requireNamespace("geodata", quietly = TRUE) && requireNamespace("terra", quietly = TRUE) &&
      (is_canonical_prod_path(env_dir) || identical(trimws(env_plantmatics_or_packera("PLANTMATICS_ENABLE_GEODATA_DOWNLOAD", "PACKERA_ENABLE_GEODATA_DOWNLOAD")), "TRUE"))) {
    ex_g <- tryCatch({
      dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
      wcr <- tryCatch({
        tg <- expand.grid(lon = unique(floor(as.numeric(df$longitude) / 30.0) * 30.0 + 15.0), lat = unique(floor(as.numeric(df$latitude) / 30.0) * 30.0 + 15.0))
        tg <- tg[is.finite(tg$lon) & is.finite(tg$lat), , drop = FALSE]
        if (nrow(tg) >= 1L && nrow(tg) <= 8L) {
          tl <- Filter(Negate(is.null), lapply(seq_len(nrow(tg)), function(ti) tryCatch(geodata::worldclim_tile("bio", tg$lon[ti], tg$lat[ti], path = env_dir), error = function(e) NULL)))
          if (length(tl) == 1L) tl[[1L]] else if (length(tl) > 1L) do.call(terra::merge, tl) else NULL
        } else NULL
      }, error = function(e) NULL)
      if (is.null(wcr)) wcr <- geodata::worldclim_global("bio", res = 2.5, path = env_dir)
      wcr <- crop_and_compress_na_rasters(wcr, env_dir = env_dir)
      sdl <- c("phh2o", "cec", "sand", "clay", "silt", "cfvo", "bdod")
      sgl <- crop_and_compress_na_rasters(Filter(Negate(is.null), c(lapply(sdl, function(v) tryCatch(geodata::soil_world(v, 5, path = env_dir), error = function(e) NULL)),
                                                                    lapply(sdl, function(v) tryCatch(geodata::soil_world(v, 15, path = env_dir), error = function(e) NULL)))), env_dir = env_dir)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      .assign_sg(assign_wc_bio_from_extracted(df, extract_rast_buffered_mean(wcr, pts, df$buffer_radius_m)),
                 extract_rast_buffered_mean(crop_and_compress_na_rasters(terra::rast(sgl), env_dir = env_dir), pts, df$buffer_radius_m), FALSE)
    }, error = function(e) { extract_err <<- e$message; NULL })
    if (!is.null(ex_g)) return(.fin_ext(ex_g))
  }
  if (!is.null(full_df_ref)) df <- full_df_ref

  if (!isTRUE(demo) && nrow(df) > 0L) {
    df <- .hydrate_env_from_sqlite(df, env_dir, c(req_cols, "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm"), demo = demo)$df
    has_cl <- ("bio1_temp_mean" %in% names(df)) && any(is.finite(suppressWarnings(as.numeric(df$bio1_temp_mean))))
    has_so <- any(vapply(intersect(c("soil_ph", "fine_soil_ph"), names(df)), function(cl) any(is.finite(suppressWarnings(as.numeric(df[[cl]])))), logical(1)))
    if (has_cl && has_so) {
      for (pm in list(c("soil_ph", "fine_soil_ph"), c("soil_sand", "fine_soil_sand"), c("soil_cec", "fine_soil_cec"))) {
        if (!(pm[1L] %in% names(df))) df[[pm[1L]]] <- rep(NA_real_, nrow(df))
        if (pm[2L] %in% names(df)) { cv <- suppressWarnings(as.numeric(df[[pm[1L]]])); df[[pm[1L]]] <- ifelse(is.finite(cv), cv, suppressWarnings(as.numeric(df[[pm[2L]]]))) }
      }
      if (!("soil_bulk_density" %in% names(df))) df$soil_bulk_density <- rep(NA_real_, nrow(df))
      return(df)
    }
  }
  if (!isTRUE(demo)) stop(sprintf("WorldClim v2.1 / SoilGrids 250m GeoTIFF rasters are missing or unreadable in '%s'%s. Populate '%s' or pass --demo.", env_dir, if (!is.null(extract_err)) paste0(" (", extract_err, ")") else "", env_dir))

  message("[DEMO MODE] Generating synthetic SoilGrids (0-15cm rhizosphere) & all 19 WorldClim bioclimatic layers...")
  set.seed(42L); lat <- as.numeric(df$latitude); lon <- as.numeric(df$longitude); n <- nrow(df)
  reg <- if ("regional_group" %in% names(df)) as.character(df$regional_group) else rep("Unknown", n)
  is_fl <- grepl("Flatrock", reg, ignore.case = TRUE) | (lon > -84 & lon < -79 & lat > 33 & lat < 36.5)
  is_sh <- grepl("Sandhill|Coastal", reg, ignore.case = TRUE) | (lon > -82 & lon < -75 & lat < 36.5)
  is_pr <- grepl("Midwest|Prairie", reg, ignore.case = TRUE) | (lon < -88)
  ph0 <- ifelse(is_sh, 4.6 + 0.3 * sin(lat), ifelse(is_fl, 5.1 + 0.25 * cos(lon), ifelse(is_pr, 7.0 + 0.3 * sin(lat), 6.3 + 0.3 * cos(lat)))) + stats::rnorm(n, 0, 0.15)
  ce0 <- ifelse(is_sh, 5.6 + 1.2 * abs(sin(lat)), ifelse(is_fl, 12.0 + 2.0 * cos(lat), ifelse(is_pr, 27.2 + 3.5 * sin(lon), 18.6 + 2.5 * cos(lat)))) + stats::rnorm(n, 0, 1.0)
  sa0 <- ifelse(is_sh, 83.0 + 5.0 * sin(lon), ifelse(is_fl, 53.0 + 6.0 * cos(lat), ifelse(is_pr, 24.5 + 4.5 * sin(lat), 38.5 + 5.0 * cos(lon)))) + stats::rnorm(n, 0, 3.0)
  bd0 <- ifelse(is_sh, 1.44, ifelse(is_fl, 1.28, 1.20)) + stats::rnorm(n, 0, 0.04)
  df$soil_ph <- round(pmax(3.8, pmin(8.5, (5.0 * ph0 + 10.0 * (ph0 + 0.15 + stats::rnorm(n, 0, 0.10))) / 15.0)), 2)
  df$soil_cec <- round(pmax(1.5, pmin(50.0, (5.0 * ce0 + 10.0 * (ce0 - 0.8 + stats::rnorm(n, 0, 0.8))) / 15.0)), 1)
  df$soil_sand <- round(pmax(5.0, pmin(92.0, (5.0 * sa0 + 10.0 * (sa0 - 1.5 + stats::rnorm(n, 0, 2.5))) / 15.0)), 1)
  df$soil_clay <- round(pmax(3.0, pmin(60.0, (100.0 - df$soil_sand) * 0.42 + stats::rnorm(n, 0, 2.0))), 1); df$soil_silt <- round(pmax(2.0, 100.0 - df$soil_sand - df$soil_clay), 1)
  df$soil_cfvo <- round(pmax(1.0, pmin(75.0, ifelse(is_fl, 28.0, ifelse(is_sh, 6.0, 11.0)) + stats::rnorm(n, 0, 2.5))), 1); df$soil_bulk_density <- round((5.0 * bd0 + 10.0 * (bd0 + 0.06 + stats::rnorm(n, 0, 0.04))) / 15.0, 2)
  df$bio1_temp_mean <- round(24.5 - 0.72 * (lat - 28.0) - 0.05 * abs(lon + 82) + stats::rnorm(n, 0, 0.4), 2); df$bio2_mean_diurnal_range <- round(pmax(4.5, 11.4 + 0.11 * abs(lon + 80.0) + stats::rnorm(n, 0, 0.45)), 2)
  df$bio4_temp_seasonality <- round(pmax(80.0, 520.0 + 28.5 * (lat - 30.0) + stats::rnorm(n, 0, 15.0)), 1); sa <- (df$bio4_temp_seasonality / 100.0) * 1.25
  df$bio5_max_temp_warmest_month <- round(df$bio1_temp_mean + sa + 0.55 * df$bio2_mean_diurnal_range + stats::rnorm(n, 0, 0.35), 2); df$bio6_min_temp_coldest_month <- round(df$bio1_temp_mean - sa - 0.55 * df$bio2_mean_diurnal_range + stats::rnorm(n, 0, 0.35), 2)
  df$bio7_temp_annual_range <- round(pmax(df$bio2_mean_diurnal_range + 1.0, df$bio5_max_temp_warmest_month - df$bio6_min_temp_coldest_month), 2); df$bio3_isothermality <- round(100.0 * df$bio2_mean_diurnal_range / pmax(df$bio7_temp_annual_range, 1.0), 2)
  df$bio10_mean_temp_warmest_qtr <- round(df$bio1_temp_mean + 0.95 * (df$bio4_temp_seasonality / 100.0) + stats::rnorm(n, 0, 0.3), 2); df$bio11_mean_temp_coldest_qtr <- round(df$bio1_temp_mean - 0.95 * (df$bio4_temp_seasonality / 100.0) + stats::rnorm(n, 0, 0.3), 2)
  df$bio8_mean_temp_wettest_qtr <- round(0.65 * df$bio10_mean_temp_warmest_qtr + 0.35 * df$bio1_temp_mean + stats::rnorm(n, 0, 0.6), 2); df$bio9_mean_temp_driest_qtr <- round(0.45 * df$bio11_mean_temp_coldest_qtr + 0.55 * df$bio1_temp_mean + stats::rnorm(n, 0, 0.6), 2)
  df$bio12_precip_annual <- round(pmax(120.0, 1450.0 - 15.0 * (lat - 30.0) + 12.0 * (lon + 85.0) + stats::rnorm(n, 0, 35.0)), 1); df$bio15_precip_seasonality <- round(pmax(5.0, 22.0 + 0.8 * (lat - 30.0) + stats::rnorm(n, 0, 1.5)), 1)
  df$bio13_precip_wettest_month <- round(pmax(15.0, df$bio12_precip_annual * (0.115 + 0.0012 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 6.0)), 1); df$bio14_precip_driest_month <- round(pmax(1.0, df$bio12_precip_annual * pmax(0.015, 0.068 - 0.0009 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 4.0)), 1)
  df$bio16_precip_wettest_qtr <- round(pmax(df$bio13_precip_wettest_month, df$bio12_precip_annual * (0.31 + 0.0020 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 12.0)), 1); df$bio17_precip_driest_qtr <- round(pmax(df$bio14_precip_driest_month, df$bio12_precip_annual * pmax(0.06, 0.21 - 0.0018 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 10.0)), 1)
  df$bio18_precip_warmest_qtr <- round(pmax(10.0, df$bio12_precip_annual * (0.32 + 0.02 * sin(lat)) + stats::rnorm(n, 0, 15.0)), 1); df$bio19_precip_coldest_qtr <- round(pmax(10.0, df$bio12_precip_annual * (0.22 - 0.01 * cos(lon)) + stats::rnorm(n, 0, 12.0)), 1)
  .ensure_wc_bio_19_cols(df)
}
