if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_env_parallel.R
# Schema constants, env/fork-safe parallel helpers, SQLite pragmas, and column registries.
# ==============================================================================

PLANTMATICS_SCHEMA_USER_VERSION <- 20261127L
PACKERA_SCHEMA_USER_VERSION <- PLANTMATICS_SCHEMA_USER_VERSION

env_plantmatics_or_packera <- function(p_key, k_key = NULL, unset = "") {
  if (is.null(k_key) || !nzchar(k_key)) {
    k_key <- sub("^PLANTMATICS_", "PACKERA_", p_key)
  }
  v <- Sys.getenv(p_key, unset = "")
  if (nzchar(trimws(v))) {
    return(v)
  }
  v2 <- Sys.getenv(k_key, unset = "")
  if (nzchar(trimws(v2))) {
    return(v2)
  }
  unset
}

resolve_packera_workers <- function(max_workers = 64L, n_tasks = NULL) {
  if (identical(Sys.getenv("PLANTMATICS_IN_FORK_WORKER", unset = ""), "1")) {
    return(1L)
  }
  if (!identical(.Platform$OS.type, "unix") || !requireNamespace("parallel", quietly = TRUE)) {
    return(1L)
  }
  env_w <- suppressWarnings(as.integer(trimws(env_plantmatics_or_packera("PLANTMATICS_MAX_WORKERS", "PACKERA_MAX_WORKERS"))))
  if (is.na(env_w) || env_w < 1L) {
    det_log <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
    det_phys <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    base_cores <- if (!is.na(det_log) && det_log >= 1L) {
      as.integer(det_log)
    } else if (!is.na(det_phys) && det_phys >= 1L) {
      as.integer(det_phys)
    } else {
      1L
    }
    env_w <- min(as.integer(max_workers), base_cores)
  } else {
    env_w <- as.integer(env_w)
  }
  if (!is.null(n_tasks) && length(n_tasks) == 1L && is.finite(n_tasks)) {
    env_w <- min(env_w, max(1L, as.integer(n_tasks)))
  }
  max(1L, as.integer(env_w))
}

packera_mclapply <- function(X, FUN, mc.cores = 1L, mc.preschedule = TRUE) {
  if (identical(Sys.getenv("PLANTMATICS_IN_FORK_WORKER", unset = ""), "1")) {
    return(lapply(X, FUN))
  }
  n_x <- length(X)
  eff_cores <- max(1L, min(as.integer(mc.cores), max(1L, n_x)))
  if (eff_cores <= 1L || n_x <= 1L || !identical(.Platform$OS.type, "unix") || !requireNamespace("parallel", quietly = TRUE)) {
    return(lapply(X, FUN))
  }
  guard_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "RCPP_PARALLEL_NUM_THREADS", "PLANTMATICS_IN_FORK_WORKER")
  prev_vals <- Sys.getenv(guard_vars, unset = NA_character_)
  restore_env <- function() {
    for (vn in guard_vars) {
      pv <- prev_vals[[vn]]
      if (is.na(pv)) {
        Sys.unsetenv(vn)
      } else {
        do.call(Sys.setenv, stats::setNames(list(pv), vn))
      }
    }
  }
  on.exit(restore_env(), add = TRUE)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1",
    PLANTMATICS_IN_FORK_WORKER = "1"
  )
  res <- parallel::mclapply(X, FUN, mc.cores = eff_cores, mc.preschedule = mc.preschedule)
  restore_env()
  bad_idx <- which(vapply(res, function(el) inherits(el, "try-error"), logical(1)))
  if (length(bad_idx) > 0L) {
    res[bad_idx] <- lapply(X[bad_idx], FUN)
  }
  res
}

.get_sqlite_user_version_r <- function(conn) {
  if (is.null(conn)) {
    return(0L)
  }
  uv <- tryCatch(DBI::dbGetQuery(conn, "PRAGMA user_version;"), error = function(e) NULL)
  if (is.null(uv) || nrow(uv) == 0L || ncol(uv) == 0L) {
    return(0L)
  }
  val <- suppressWarnings(as.integer(uv[[1]][1]))
  if (is.na(val)) 0L else val
}

.apply_packera_sqlite_pragmas_r <- function(conn, read_only = FALSE) {
  if (is.null(conn)) {
    return(invisible(NULL))
  }
  tryCatch(
    {
      DBI::dbExecute(conn, "PRAGMA busy_timeout = 30000;")
      DBI::dbExecute(conn, "PRAGMA foreign_keys = ON;")
      DBI::dbExecute(conn, "PRAGMA cache_size = -262144;")
      DBI::dbExecute(conn, "PRAGMA temp_store = MEMORY;")
      tryCatch(DBI::dbExecute(conn, "PRAGMA mmap_size = 2147483648;"), error = function(e) NULL)
      if (!isTRUE(read_only)) {
        DBI::dbExecute(conn, "PRAGMA journal_mode = WAL;")
        DBI::dbExecute(conn, "PRAGMA synchronous = NORMAL;")
        tryCatch(DBI::dbExecute(conn, "PRAGMA wal_autocheckpoint = 4000;"), error = function(e) NULL)
      }
    },
    error = function(e) NULL
  )
  invisible(NULL)
}

PACKERA_DB_LOGICAL_COLS <- c(
  "is_primary_duplicate", "is_calibrated_mm", "is_bent_stem",
  "is_dissected", "lm3_gates_pass", "is_crop_clipped", "bilateral_lobe_symmetry", "efa_qc_pass", "mask_is_upright", "is_verified_anthesis",
  "is_rock_outcrop", "is_azonal_lithology", "high_coord_uncertainty_flag",
  "is_training_eligible_coord", "is_clim_eligible_coord", "is_edaphic_eligible_coord", "azonal_lithology_presence_buffer",
  "is_promoted_bronze", "is_passive",
  "is_contour_outlier", "morphometric_outlier_flag", "misidentification_flag",
  "substrate_specialist_flag", "multimodal_conflict_flag",
  "is_completed_voucher", "is_type_specimen", "has_empirical_type_anchor",
  "is_typical_morphotype", "matches_type_lithology",
  "is_confident_label_noise", "cryptic_cluster_cohesion_pass",
  "cryptic_non_clinal_pass", "hybrid_zone_cohesion_pass",
  "geno_quality_pass", "is_putative_novel_or_hybrid_query", "is_geno_admixed", "cytonuclear_discordant",
  "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass",
  "is_nwi_wetland", "is_hybrid_taxon", "is_native", "is_exotic", "is_historical_or_extirpated",
  "is_centroid_or_inst_coord", "is_open_or_hybrid_concept", "is_abaxial_pressed_leaf", "env_mcd_outlier_flag"
)

PACKERA_DB_HARMONIC_COLS <- as.vector(vapply(
  1:24,
  function(i) paste0(c("A", "B", "C", "D"), i),
  character(4)
))
PACKERA_DB_SYM_HARMONIC_COLS <- as.vector(vapply(
  1:24,
  function(i) paste0(c("A", "D"), i),
  character(2)
))
PACKERA_DB_PC_COLS <- paste0("PC", 1:15)
PACKERA_DB_WC_BIO_COLS <- c(
  "bio1_temp_mean",
  "bio2_mean_diurnal_range",
  "bio3_isothermality",
  "bio4_temp_seasonality",
  "bio5_max_temp_warmest_month",
  "bio6_min_temp_coldest_month",
  "bio7_temp_annual_range",
  "bio8_mean_temp_wettest_qtr",
  "bio9_mean_temp_driest_qtr",
  "bio10_mean_temp_warmest_qtr",
  "bio11_mean_temp_coldest_qtr",
  "bio12_precip_annual",
  "bio13_precip_wettest_month",
  "bio14_precip_driest_month",
  "bio15_precip_seasonality",
  "bio16_precip_wettest_qtr",
  "bio17_precip_driest_qtr",
  "bio18_precip_warmest_qtr",
  "bio19_precip_coldest_qtr"
)

PACKERA_DB_ENV_CACHE_COLS <- c(
  "lat_lon_key",
  "catalogNumber",
  "latitude",
  "longitude",
  "musym",
  "muname",
  "taxsuborder",
  "taxgreatgroup",
  "is_rock_outcrop",
  "outcrop_class",
  "ssurgo_rock_outcrop_pct",
  "ssurgo_lithic_series_pct",
  "ssurgo_bedrock_depth_cm",
  "ssurgo_sand_pct",
  "ssurgo_clay_pct",
  "ssurgo_ph",
  "ssurgo_om_pct",
  "ssurgo_cec",
  "ssurgo_aws_0_50cm",
  "polaris_ph_30m",
  "polaris_sand_30m",
  "polaris_clay_30m",
  "polaris_om_30m",
  "polaris_bd_30m",
  "polaris_ksat_30m",
  "solus_resdept_cm",
  "solus_fragvol_pct",
  "solus_sand_pct",
  "solus_clay_pct",
  "solus_silt_pct",
  "solus_cec",
  "topo_tpi_30m",
  "topo_slope_deg",
  "topo_hli_30m",
  "topo_twi_30m",
  "elevation_m",
  "fine_soil_ph",
  "fine_soil_sand",
  "fine_soil_clay",
  "soil_silt_pct",
  "ilr_texture_1",
  "ilr_texture_2",
  "fine_soil_cec",
  "fine_cfvo_pct",
  "fine_bedrock_depth_cm",
  "sgmc_lithology_class",
  "is_azonal_lithology",
  "buffer_radius_m",
  "high_coord_uncertainty_flag",
  "is_training_eligible_coord",
  "is_clim_eligible_coord",
  "is_edaphic_eligible_coord",
  "edaphic_coord_weight",
  "ph_q10_buffer",
  "ph_q90_buffer",
  "min_bedrock_depth_buffer_cm",
  "max_cfvo_buffer_pct",
  "azonal_lithology_presence_buffer",
  PACKERA_DB_WC_BIO_COLS,
  "soil_ph",
  "soil_cec",
  "soil_sand",
  "soil_clay",
  "soil_silt",
  "soil_cfvo",
  "soil_bulk_density",
  "soil_ph_soilgrids",
  "soil_cec_soilgrids",
  "soil_sand_soilgrids",
  "soil_clay_soilgrids",
  "soil_silt_soilgrids",
  "soil_cfvo_soilgrids",
  "soil_bulk_density_soilgrids",
  "soil_ca_mg_molar_ratio",
  "soil_ultramafic_ni_cr_index",
  "soil_carbonate_c_pct",
  "is_nwi_wetland",
  "nwi_wetland_system",
  "nwi_wetland_dist_m"
)

PACKERA_DB_SRF_COLS <- c(
  "run_id",
  "catalogNumber",
  "plant_individual_id",
  "lat_lon_key",
  "ecotype_regime",
  "lithic_outcrop_index",
  "summer_moisture_deficit_idx",
  "maritime_continentality_idx",
  "rf_predicted_taxon",
  "rf_predicted_prob",
  "rf_given_taxon_prob",
  "rf_prediction_source",
  "edaphic_best_fit_taxon",
  "spatial_block_cv_acc",
  "spatial_autocorr_range_km",
  "spatial_block_id",
  "doy",
  "expected_doy",
  "hopkins_expected",
  "delta_doy",
  "hopkins_adjusted_doy",
  "hopkins_adjusted_doy_resid",
  "anthesis_dominance_weight",
  "allochronic_bimodality_delta_bic",
  "pheno_timing_category",
  "anthesis_filter_status",
  "edaphic_z_score",
  "climatic_z_score",
  "type_edaphic_z_score",
  "nearest_type_niche_taxon",
  "nearest_type_niche_basionym",
  "nearest_type_basionym_edaphic_z",
  "matches_type_lithology",
  "primary_env_predictors",
  "primary_clim_predictors",
  "primary_edaph_predictors",
  PACKERA_DB_WC_BIO_COLS,
  "substrate_specialist_flag",
  "n_observed_modalities",
  "multimodal_concordance",
  "multimodal_conflict_flag",
  "triage_category",
  "conflict_discordance_notes",
  "env_mahal_d2",
  "env_mcd_outlier_flag"
)

.normalize_semicolon_col_r <- function(v) {
  if (is.null(v)) {
    return(NULL)
  }
  if (is.list(v)) {
    return(vapply(v, function(el) {
      if (is.null(el) || length(el) == 0L) {
        return(NA_character_)
      }
      toks <- trimws(as.character(unlist(el)))
      toks <- toks[!is.na(toks) & nzchar(toks)]
      if (length(toks) == 0L) NA_character_ else paste(toks, collapse = ";")
    }, character(1)))
  }
  as.character(v)
}
