#!/usr/bin/env Rscript
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Script: scripts/analysis/edaphic_fine_grain_utils.R
# Project: Plantmatics / Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Core fine-grain edaphic constants, 0-15 cm rhizosphere integration,
#   isometric log-ratio (ILR) soil texture, coordinate-uncertainty gating,
#   uncertainty-buffered raster extrema, and derived habitat indices, plus
#   idempotent auto-sourcing of Step 05 edaphic/SSURGO/SOLUS/3DEP/collinearity
#   companion modules.
# ==============================================================================

.resolve_edaphic_helper_path <- function(rel_root, rel_self) {
  args_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  ofiles <- unlist(lapply(sys.frames(), function(fr) fr$ofile))
  s_dirs <- unique(dirname(c(args_f, ofiles)))
  s_dirs <- s_dirs[!is.na(s_dirs) & nzchar(s_dirs)]
  cands <- unique(c(rel_root, if (length(s_dirs) > 0L) file.path(s_dirs, rel_self) else character(0)))
  for (.p in cands) if (!is.na(.p) && nzchar(.p) && file.exists(.p)) return(.p)
  NA_character_
}

.source_edaphic_dep <- function(sentinel_fn, rel_root, rel_self) {
  if (!exists(sentinel_fn, mode = "function")) {
    .hit <- .resolve_edaphic_helper_path(rel_root, rel_self)
    if (!is.na(.hit)) source(.hit)
  }
}

.source_edaphic_dep("query_environmental_cache_sqlite", "scripts/core/db_utils.R", file.path("..", "core", "db_utils.R"))
.source_edaphic_dep("clamp_bedrock_buffer_radius_m", "scripts/analysis/step05_layer_hygiene_utils.R", "step05_layer_hygiene_utils.R")
.source_edaphic_dep("enrich_geochem_wetland_predictors", "scripts/analysis/edaphic_geochem_wetland_utils.R", "edaphic_geochem_wetland_utils.R")
.source_edaphic_dep("load_step05_macroecology_config", "scripts/analysis/step05_config_geo_utils.R", "step05_config_geo_utils.R")

SSURGO_OUTCROP_SERIES_REGEX <- "(^|[^[:alnum:]_])(rock outcrop|rock land|flatrock|glade|boulders|wake|ashlar|gladeville|barfield|goldston|manteo|ramsey|lily|chrome|conowingo)($|[^[:alnum:]_])"

SGMC_LITHOLOGY_CLASSES <- c(
  "Ultramafic_Serpentine",
  "Mafic_Igneous_Metamorphic",
  "Carbonate_Calcareous",
  "Felsic_Granitic_Gneiss",
  "Siliciclastic_Consolidated",
  "Unconsolidated_Sedimentary"
)

SSURGO_OUTCROP_CLASSES <- c(
  "Granitic flatrock",
  "Limestone glade",
  "Sandstone/Siliceous",
  "Serpentine",
  "Diabase/Mafic glade",
  "Non-outcrop matrix"
)

CORE_FINE_GRAIN_CACHE_COLS <- c(
  "catalogNumber", "lat_lon_key", "latitude", "longitude",
  "musym", "muname", "taxsuborder", "taxgreatgroup", "is_rock_outcrop", "outcrop_class",
  "ssurgo_rock_outcrop_pct", "ssurgo_lithic_series_pct", "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct", "ssurgo_clay_pct",
  "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm",
  "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m",
  "polaris_bd_30m", "polaris_ksat_30m",
  "solus_resdept_cm", "solus_fragvol_pct", "solus_sand_pct", "solus_clay_pct", "solus_silt_pct", "solus_cec",
  "topo_tpi_30m", "topo_slope_deg", "topo_hli_30m", "topo_twi_30m", "elevation_m",
  "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "soil_silt_pct",
  "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm"
)

FINE_GRAIN_CACHE_COLS <- c(
  CORE_FINE_GRAIN_CACHE_COLS,
  "sgmc_lithology_class", "is_azonal_lithology",
  "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
  "is_clim_eligible_coord", "is_edaphic_eligible_coord", "edaphic_coord_weight",
  "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm",
  "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer",
  "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
  "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m"
)

format_lat_lon_key <- function(lat, lon) {
  lat_num <- suppressWarnings(as.numeric(lat))
  lon_num <- suppressWarnings(as.numeric(lon))
  ifelse(is.na(lat_num) | is.na(lon_num), NA_character_, sprintf("%.5f_%.5f", round(lat_num, 5), round(lon_num, 5)))
}

# ------------------------------------------------------------------------------
# 1. Rhizosphere Integration (0-15 cm Trapezoidal) & ILR Soil Texture
# ------------------------------------------------------------------------------
compute_rhizosphere_0_15cm <- function(X_0_5, X_5_15) {
  x1 <- suppressWarnings(as.numeric(X_0_5))
  x2 <- suppressWarnings(as.numeric(X_5_15))
  n <- max(length(x1), length(x2))
  if (n == 0L) return(numeric(0))
  if (length(x1) != n) x1 <- rep_len(x1, n)
  if (length(x2) != n) x2 <- rep_len(x2, n)
  ok1 <- !is.na(x1) & is.finite(x1)
  ok2 <- !is.na(x2) & is.finite(x2)
  out <- rep(NA_real_, n)
  both <- ok1 & ok2
  if (any(both)) out[both] <- (5.0 * x1[both] + 10.0 * x2[both]) / 15.0
  only1 <- ok1 & !ok2; if (any(only1)) out[only1] <- x1[only1]
  only2 <- !ok1 & ok2; if (any(only2)) out[only2] <- x2[only2]
  out
}

compute_soil_texture_log_ratios <- function(sand_pct, clay_pct, silt_pct = NULL) {
  s <- suppressWarnings(as.numeric(sand_pct))
  c_v <- suppressWarnings(as.numeric(clay_pct))
  n <- max(length(s), length(c_v))
  if (n == 0L) return(data.frame(soil_silt_pct = numeric(0), ilr_texture_1 = numeric(0), ilr_texture_2 = numeric(0)))
  if (length(s) != n) s <- rep_len(s, n)
  if (length(c_v) != n) c_v <- rep_len(c_v, n)
  si_raw <- if (is.null(silt_pct)) rep(NA_real_, n) else rep_len(suppressWarnings(as.numeric(silt_pct)), n)

  valid <- is.finite(s) & is.finite(c_v) & (s >= 0.0) & (c_v >= 0.0)
  si_out <- rep(NA_real_, n); ilr1 <- rep(NA_real_, n); ilr2 <- rep(NA_real_, n)
  if (any(valid)) {
    sv <- pmax(0.5, s[valid]); cv <- pmax(0.5, c_v[valid]); siv_in <- si_raw[valid]
    siv <- ifelse(is.finite(siv_in) & siv_in >= 0.0, pmax(0.5, siv_in), pmax(0.5, 100.0 - sv - cv))
    tot <- sv + cv + siv
    sv_c <- 100.0 * (sv / tot); cv_c <- 100.0 * (cv / tot); siv_c <- 100.0 * (siv / tot)
    si_out[valid] <- round(siv_c, 2)
    ilr1[valid] <- round(sqrt(2.0 / 3.0) * log(sv_c / sqrt(siv_c * cv_c)), 4)
    ilr2[valid] <- round(sqrt(0.5) * log(siv_c / cv_c), 4)
  }
  data.frame(soil_silt_pct = si_out, ilr_texture_1 = ilr1, ilr_texture_2 = ilr2, stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Coordinate Uncertainty Gating & Shared Fine-Raster Utilities
# ------------------------------------------------------------------------------
compute_uncertainty_buffer_radius <- function(coord_unc_m, default_r = 250.0, min_r = 100.0, max_r = 1000.0) {
  u <- suppressWarnings(as.numeric(coord_unc_m))
  ifelse(is.na(u) | !is.finite(u) | u <= 0, as.numeric(default_r), pmin(as.numeric(max_r), pmax(as.numeric(min_r), u)))
}

flag_coordinate_uncertainty <- function(
  df,
  max_uncertainty_m = 10000.0,
  max_edaphic_uncertainty_m = 2000.0,
  default_r = 250.0,
  min_r = 100.0,
  max_r = 1000.0
) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0L) return(df)

  unc_raw <- if ("coordinateUncertaintyInMeters" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
  } else if ("coordinateUncertainty" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertainty))
  } else rep(NA_real_, n)
  if ("coordinateUncertainty" %in% names(df)) {
    unc_alt <- suppressWarnings(as.numeric(df$coordinateUncertainty))
    na_m <- is.na(unc_raw) & !is.na(unc_alt)
    if (any(na_m)) unc_raw[na_m] <- unc_alt[na_m]
  }

  is_t_u <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    as.logical(.identify_type_specimen_vec(df)$is_type_specimen)
  } else if ("is_type_specimen" %in% names(df)) {
    as.logical(df$is_type_specimen)
  } else rep(FALSE, n)
  is_t_u[is.na(is_t_u)] <- FALSE

  q_c <- rep(FALSE, n)
  if ("is_centroid_or_inst_coord" %in% names(df)) {
    raw_qc <- df$is_centroid_or_inst_coord
    q_c <- if (is.logical(raw_qc)) {
      !is.na(raw_qc) & raw_qc
    } else if (is.numeric(raw_qc) || is.integer(raw_qc)) {
      !is.na(raw_qc) & is.finite(raw_qc) & (raw_qc != 0)
    } else {
      tolower(trimws(as.character(raw_qc))) %in% c("true", "1", "t", "yes")
    }
    q_c[is.na(q_c)] <- FALSE
  }

  df$coordinateUncertaintyInMeters <- unc_raw
  df$high_coord_uncertainty_flag <- !is.na(unc_raw) & is.finite(unc_raw) & (unc_raw > as.numeric(max_uncertainty_m))
  df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
  df$is_clim_eligible_coord <- (!df$high_coord_uncertainty_flag & !q_c) |
                               (is_t_u & !q_c & (is.na(unc_raw) | !is.finite(unc_raw) | unc_raw <= 50000.0))

  yr_raw <- if ("year" %in% names(df)) suppressWarnings(as.numeric(df$year)) else rep(NA_real_, n)
  has_year_col <- ("year" %in% names(df)) && any(is.finite(yr_raw))
  modern_gps_fallback <- if (has_year_col) (is.finite(yr_raw) & yr_raw >= 2000) else rep(TRUE, n)
  df$is_edaphic_eligible_coord <- (is.finite(unc_raw) & unc_raw > 0 & (unc_raw <= as.numeric(max_edaphic_uncertainty_m))) |
                                  ((is.na(unc_raw) | !is.finite(unc_raw)) & modern_gps_fallback)
  df$edaphic_coord_weight <- ifelse(
    is.finite(unc_raw) & unc_raw > 0,
    ifelse(unc_raw <= 500.0, 1.0, ifelse(unc_raw <= 2000.0, 0.70, ifelse(unc_raw <= as.numeric(max_edaphic_uncertainty_m), 0.45, 0.0))),
    ifelse(modern_gps_fallback, 0.60, 0.0)
  )
  df$buffer_radius_m <- compute_uncertainty_buffer_radius(unc_raw, default_r = default_r, min_r = min_r, max_r = max_r)

  # Quarantined centroids/rooftops and > 10 km coords are edaphically ineligible regardless of type status
  bad_edaph <- q_c | df$high_coord_uncertainty_flag
  if (any(bad_edaph)) {
    df$is_edaphic_eligible_coord[bad_edaph] <- FALSE
    df$edaphic_coord_weight[bad_edaph] <- 0.0
  }
  df
}

.list_real_fine_rasters <- function(dir_path, exclude_pattern = "^dem_3dep_regional\\.tif$") {
  if (is.null(dir_path) || !nzchar(dir_path) || !dir.exists(dir_path)) return(character(0))
  fls <- list.files(dir_path, pattern = "\\.(tif|tiff|hgt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fls) == 0L) return(character(0))
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) fls <- fls[!grepl(exclude_pattern, basename(fls), ignore.case = TRUE)]
  fls[file.exists(fls) & file.info(fls)$size > 256]
}

.project_vect_to_rast <- function(v_geom, r_rast) {
  if (is.null(v_geom) || is.null(r_rast)) return(v_geom)
  r_crs <- tryCatch(terra::crs(r_rast), error = function(e) "")
  v_crs <- tryCatch(terra::crs(v_geom), error = function(e) "")
  if (nzchar(r_crs) && nzchar(v_crs) && !isTRUE(tryCatch(terra::same.crs(v_geom, r_rast), error = function(e) FALSE))) {
    v_geom <- tryCatch(terra::project(v_geom, r_crs), error = function(e) v_geom)
  }
  v_geom
}

.rescale_solus_ph_vec <- function(v_raw) {
  v <- suppressWarnings(as.numeric(v_raw))
  fin <- v[is.finite(v)]
  if (length(fin) == 0L) return(v)
  med_v <- stats::median(fin, na.rm = TRUE)
  if (is.finite(med_v) && med_v > 100.0) v <- v / 100.0 else if (is.finite(med_v) && med_v > 14.5) v <- v / 10.0
  v[!is.finite(v) | v < 2.0 | v > 11.5] <- NA_real_
  v
}

# ------------------------------------------------------------------------------
# 3. Uncertainty-Buffered Herbarium Extraction (Marcer et al. 2022)
# ------------------------------------------------------------------------------
extract_uncertainty_buffered_predictors <- function(df, env_dir = "data/environmental", demo = FALSE) {
  df <- flag_coordinate_uncertainty(df)
  n <- nrow(df)
  if (n == 0L) return(df)
  if (!("sgmc_lithology_class" %in% names(df))) df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = demo)

  df$ph_q10_buffer <- rep(NA_real_, n)
  df$ph_q90_buffer <- rep(NA_real_, n)
  df$min_bedrock_depth_buffer_cm <- rep(NA_real_, n)
  df$max_cfvo_buffer_pct <- rep(NA_real_, n)
  df$azonal_lithology_presence_buffer <- rep(NA, n)

  if (!isTRUE(demo)) {
    pol_files <- .list_real_fine_rasters(file.path(env_dir, "gssurgo_30m"))
    sol_files <- c(.list_real_fine_rasters(file.path(env_dir, "solus_100m")), pol_files)
    lat_v <- suppressWarnings(as.numeric(df$latitude))
    lon_v <- suppressWarnings(as.numeric(df$longitude))
    unc_b <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    valid_pts <- is.finite(lat_v) & is.finite(lon_v) &
                 (df$is_edaphic_eligible_coord | is.na(unc_b) | (is.finite(unc_b) & unc_b <= 5000.0))
    if (any(valid_pts) && requireNamespace("terra", quietly = TRUE) && (length(pol_files) > 0L || length(sol_files) > 0L)) {
      v_idx <- which(valid_pts)
      pts_v <- tryCatch(terra::vect(df[v_idx, c("longitude", "latitude"), drop = FALSE], geom = c("longitude", "latitude"), crs = "EPSG:4326"), error = function(e) NULL)
      buf_v <- if (!is.null(pts_v)) tryCatch(terra::buffer(pts_v, width = pmax(100.0, pmin(1000.0, as.numeric(df$buffer_radius_m[v_idx])))), error = function(e) NULL) else NULL
      buf_bd_v <- if (!is.null(pts_v)) tryCatch(terra::buffer(pts_v, width = clamp_bedrock_buffer_radius_m(df$buffer_radius_m[v_idx], min_r = 100.0, max_bedrock_r = 300.0)), error = function(e) NULL) else NULL

      if (!is.null(buf_v)) {
        ph_tifs <- c(pol_files[grepl("^ph_", basename(pol_files), ignore.case = TRUE)],
                     sol_files[grepl("^ph1to1h2o.*_5_cm_p|^ph1to1h2o|^ph_", basename(sol_files), ignore.case = TRUE)])
        if (length(ph_tifs) > 0L) {
          tryCatch({
            r_ph <- terra::rast(ph_tifs[1L]); buf_ph <- .project_vect_to_rast(buf_v, r_ph)
            df$ph_q10_buffer[v_idx] <- round(as.numeric(.rescale_solus_ph_vec(terra::extract(r_ph, buf_ph, fun = function(x) stats::quantile(x, probs = 0.10, na.rm = TRUE))[[2L]])), 2)
            df$ph_q90_buffer[v_idx] <- round(as.numeric(.rescale_solus_ph_vec(terra::extract(r_ph, buf_ph, fun = function(x) stats::quantile(x, probs = 0.90, na.rm = TRUE))[[2L]])), 2)
          }, error = function(e) NULL)
        }
        # Read strictly lithic/paralithic bedrock contact GeoTIFF (anylithicdpt_cm_p.tif)
        res_tifs <- sol_files[grepl("anylithicdpt", basename(sol_files), ignore.case = TRUE)]
        if (length(res_tifs) > 0L && !is.null(buf_bd_v)) {
          min_accum <- rep(Inf, length(v_idx))
          for (rt_p in res_tifs) {
            tryCatch({
              r_bd <- terra::rast(rt_p); buf_bd <- .project_vect_to_rast(buf_bd_v, r_bd)
              ex_m <- suppressWarnings(as.numeric(terra::extract(r_bd, buf_bd, fun = function(x) {
                xf <- x[is.finite(x)]
                if (length(xf) == 0L) NA_real_ else min(xf)
              })[[2L]]))
              ok_m <- is.finite(ex_m) & ex_m >= 0.0
              min_accum[ok_m] <- pmin(min_accum[ok_m], ex_m[ok_m])
            }, error = function(e) NULL)
          }
          min_accum[!is.finite(min_accum)] <- NA_real_
          df$min_bedrock_depth_buffer_cm[v_idx] <- round(pmin(200.0, pmax(0.0, min_accum)), 1)
        }
        cf_tifs <- sol_files[grepl("fragvol_0_cm_p|fragvol_5_cm_p|fragvol_15_cm_p|fragvol|cfvo", basename(sol_files), ignore.case = TRUE)]
        if (length(cf_tifs) > 0L) {
          max_cf <- rep(-Inf, length(v_idx))
          for (cf_p in utils::head(cf_tifs, 3L)) {
            tryCatch({
              r_cf <- terra::rast(cf_p); buf_cf <- .project_vect_to_rast(buf_v, r_cf)
              ex_c <- suppressWarnings(as.numeric(terra::extract(r_cf, buf_cf, fun = function(x) max(x, na.rm = TRUE))[[2L]]))
              ok_c <- is.finite(ex_c) & ex_c >= 0.0
              max_cf[ok_c] <- pmax(max_cf[ok_c], ex_c[ok_c])
            }, error = function(e) NULL)
          }
          max_cf[!is.finite(max_cf)] <- NA_real_
          df$max_cfvo_buffer_pct[v_idx] <- round(pmin(100.0, pmax(0.0, max_cf)), 1)
        }
      }
    }
    return(df)
  }

  # Explicit demonstration simulation mode ONLY (demo == TRUE)
  r_m <- df$buffer_radius_m
  ph_mean <- if ("fine_soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_ph)) else if ("soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph)) else rep(5.6, n)
  ph_mean[!is.finite(ph_mean)] <- 5.6
  bed_mean <- if ("fine_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$fine_bedrock_depth_cm)) else rep(150.0, n)
  bed_mean[!is.finite(bed_mean)] <- 150.0
  cfvo_mean <- if ("fine_cfvo_pct" %in% names(df)) suppressWarnings(as.numeric(df$fine_cfvo_pct)) else if ("solus_fragvol_pct" %in% names(df)) suppressWarnings(as.numeric(df$solus_fragvol_pct)) else rep(10.0, n)
  cfvo_mean[!is.finite(cfvo_mean)] <- 10.0
  out_pct <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  out_pct[!is.finite(out_pct)] <- 0.0
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n); is_out[is.na(is_out)] <- FALSE
  rough_v <- if ("topo_slope_deg" %in% names(df)) suppressWarnings(as.numeric(df$topo_slope_deg)) else rep(3.0, n); rough_v[!is.finite(rough_v)] <- 3.0
  lith_cls <- as.character(df$sgmc_lithology_class)
  is_azonal <- if ("is_azonal_lithology" %in% names(df)) as.logical(df$is_azonal_lithology) else rep(FALSE, n); is_azonal[is.na(is_azonal)] <- FALSE

  buf_scale <- pmax(0.15, log10(pmax(100.0, r_m) / 50.0))
  r_m_bed <- clamp_bedrock_buffer_radius_m(r_m, min_r = 100.0, max_bedrock_r = 300.0)
  buf_scale_bed <- pmax(0.15, log10(pmax(100.0, r_m_bed) / 50.0))
  het_ph <- pmin(1.15, pmax(0.12, (0.22 + 0.015 * rough_v + 0.006 * out_pct) * buf_scale))
  acid_pull <- ifelse(lith_cls %in% c("Felsic_Granitic_Gneiss", "Siliciclastic_Consolidated") & (is_out | out_pct >= 15.0), 0.35, 0.0)
  calc_pull <- ifelse(lith_cls %in% c("Carbonate_Calcareous", "Mafic_Igneous_Metamorphic", "Ultramafic_Serpentine") & (is_out | out_pct >= 15.0), 0.40, 0.0)

  df$ph_q10_buffer <- round(pmax(3.5, pmin(ph_mean, ph_mean - het_ph - acid_pull)), 2)
  df$ph_q90_buffer <- round(pmin(8.6, pmax(ph_mean, ph_mean + het_ph + calc_pull)), 2)
  shallow_pull <- ifelse(is_out | out_pct >= 20.0, 0.45, ifelse(out_pct >= 5.0 | is_azonal, 0.70, 0.88)) * (1.0 - 0.12 * (buf_scale_bed - 0.3))
  df$min_bedrock_depth_buffer_cm <- round(pmax(0.0, pmin(bed_mean, bed_mean * pmax(0.15, shallow_pull) - 0.25 * out_pct)), 1)
  df$max_cfvo_buffer_pct <- round(pmin(95.0, pmax(cfvo_mean, cfvo_mean + (8.0 + 0.35 * out_pct + 0.55 * rough_v) * buf_scale)), 1)
  df$azonal_lithology_presence_buffer <- as.logical(is_azonal | is_out | (out_pct >= 15.0) | (df$min_bedrock_depth_buffer_cm <= 40.0))
  df
}

# ------------------------------------------------------------------------------
# 4. Deterministic Multi-Scale Habitat & Hydrothermal Covariate Synthesis
# ------------------------------------------------------------------------------
compute_derived_habitat_predictors <- function(df, db_conn = NULL, env_dir = "data/environmental") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  n <- nrow(df); nms <- names(df)
  num_or_na <- function(col_nm) if (col_nm %in% nms) suppressWarnings(as.numeric(df[[col_nm]])) else rep(NA_real_, n)

  # 1. Smooth continuous lithic outcrop index [0, 1] from fine_bedrock_depth_cm + optional fine_cfvo_pct
  if ("fine_bedrock_depth_cm" %in% nms) {
    bd_cm <- num_or_na("fine_bedrock_depth_cm")
    cfvo_pct <- num_or_na("fine_cfvo_pct")
    ok_bd <- is.finite(bd_cm)
    if (sum(ok_bd) >= 1L) {
      bd_score <- pmin(1.0, pmax(0.0, 1.0 - log1p(pmax(bd_cm, 0.0)) / log1p(200.0)))
      cf_term <- log1p(pmax(cfvo_pct, 0.0)) / log1p(100.0)
      lith_idx <- ifelse(
        ok_bd,
        ifelse(is.finite(cfvo_pct), round(pmin(1.0, pmax(0.0, 0.80 * bd_score + 0.20 * cf_term)), 4), round(bd_score, 4)),
        NA_real_
      )
      if (n <= 2L || (sum(is.finite(lith_idx)) >= 2L && is.finite(stats::sd(lith_idx, na.rm = TRUE)) && stats::sd(lith_idx, na.rm = TRUE) > 1e-6)) {
        df$lithic_outcrop_index <- lith_idx
      }
    }
  }

  # 2. Warm-season / dry-month hydrothermal moisture deficit index
  if (all(c("bio10_mean_temp_warmest_qtr", "bio14_precip_driest_month") %in% nms)) {
    b10 <- num_or_na("bio10_mean_temp_warmest_qtr"); b14 <- num_or_na("bio14_precip_driest_month"); b17 <- num_or_na("bio17_precip_driest_qtr")
    dry_p <- ifelse(is.finite(b17), b17 / 3.0, b14); ok_sm <- is.finite(b10) & is.finite(dry_p)
    if (sum(ok_sm) >= 3L) {
      smd_idx <- ifelse(ok_sm, round(pmax(4.0, 2.5 * pmax(0.0, b10)) / pmax(dry_p + 10.0, 5.0), 4), NA_real_)
      if (sum(is.finite(smd_idx)) >= 3L && is.finite(stats::sd(smd_idx, na.rm = TRUE)) && stats::sd(smd_idx, na.rm = TRUE) > 1e-6) {
        df$summer_moisture_deficit_idx <- smd_idx
      }
    }
  }

  # 3. Maritime vs. continental thermal buffering index
  if (all(c("bio4_temp_seasonality", "bio6_min_temp_coldest_month") %in% nms)) {
    b4 <- num_or_na("bio4_temp_seasonality"); b6 <- num_or_na("bio6_min_temp_coldest_month")
    b2 <- num_or_na("bio2_mean_diurnal_range"); if (all(!is.finite(b2))) b2 <- num_or_na("bio2_diurnal_range")
    ok_mc <- is.finite(b4) & is.finite(b6)
    if (sum(ok_mc) >= 3L) {
      b4_norm <- ifelse(stats::median(abs(b4[ok_mc]), na.rm = TRUE) > 50.0, b4 / 100.0, b4)
      diurn_term <- ifelse(is.finite(b2), pmax(4.0, b2) / 12.0, 1.0)
      mc_idx <- ifelse(ok_mc, round((b6 + 20.0) / pmax(b4_norm * diurn_term, 1.0), 4), NA_real_)
      if (sum(is.finite(mc_idx)) >= 3L && is.finite(stats::sd(mc_idx, na.rm = TRUE)) && stats::sd(mc_idx, na.rm = TRUE) > 1e-6) {
        df$maritime_continentality_idx <- mc_idx
      }
    }
  }

  if (exists("enrich_geochem_wetland_predictors", mode = "function")) {
    df <- enrich_geochem_wetland_predictors(df, db_conn = db_conn, env_dir = env_dir)
  }
  df
}

# ------------------------------------------------------------------------------
# 5. Idempotent Auto-Sourcing of Companion Edaphic Sub-Modules
# ------------------------------------------------------------------------------
.source_edaphic_dep("query_ssurgo_sda", "scripts/analysis/edaphic_ssurgo_sda_utils.R", "edaphic_ssurgo_sda_utils.R")
.source_edaphic_dep("aggregate_ssurgo_tabular_to_points", "scripts/analysis/edaphic_ssurgo_sda_utils.R", "edaphic_ssurgo_sda_utils.R")
.source_edaphic_dep("extract_fine_grain_edaphic_metrics", "scripts/analysis/edaphic_solus_3dep_utils.R", "edaphic_solus_3dep_utils.R")
.source_edaphic_dep("screen_edaphic_point_mass_columns", "scripts/analysis/edaphic_collinearity_contingency_utils.R", "edaphic_collinearity_contingency_utils.R")
