#!/usr/bin/env Rscript
# ==============================================================================
# Script: 05_multimodal_spatial_rf.R
# Project: Plantmatics / Packera Species Delimitation & Macroecology Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
# ==============================================================================
suppressPackageStartupMessages({
  for (.pkg in c("soilDB", "elevatr", "geodata", "ranger", "adespatial", "blockCV", "sp", "adehabitatMA", "adehabitatHR", "ecospat", "digest", "spatialRF", "ENMTools", "vegan", "terra", "sf", "randomForest", "circular", "MASS", "ggridges", "dplyr", "readr", "tibble", "ggplot2", "gridExtra", "optparse")) {
    if (requireNamespace(.pkg, quietly = TRUE)) suppressMessages(library(.pkg, character.only = TRUE))
  }
})
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.step05_script_dirs <- unique(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))))
.step05_script_dirs <- .step05_script_dirs[!is.na(.step05_script_dirs) & nzchar(.step05_script_dirs)]
.source_step05_mod <- function(repo_rel, script_rel = basename(repo_rel), guard_syms = character(0), optional = FALSE) {
  if (length(guard_syms) > 0L && all(vapply(guard_syms, function(s) exists(s, mode = "function"), logical(1)))) return(invisible(TRUE))
  for (p in unique(c(repo_rel, file.path(.step05_script_dirs, script_rel)))) {
    if (!is.na(p) && nzchar(p) && file.exists(p)) { source(p, local = FALSE); return(invisible(TRUE)) }
  }
  if (!isTRUE(optional)) warning(sprintf("Step 05 helper module not found: %s", repo_rel))
  invisible(FALSE)
}
.source_step05_mod("scripts/morphometrics/pub_figures_utils.R", file.path("..", "morphometrics", "pub_figures_utils.R"), c("get_dynamic_packera_palette", "standardize_packera_taxon", "format_predictor_labels"))
.source_step05_mod("scripts/morphometrics/authority_weighting_utils.R", file.path("..", "morphometrics", "authority_weighting_utils.R"), c("compute_authority_precision_weights", "resolve_ecological_authority_weights"))
.source_step05_mod("scripts/core/db_utils.R", file.path("..", "core", "db_utils.R"), c("persist_step05_to_sqlite"))
.source_step05_mod("scripts/analysis/step05_config_geo_utils.R", "step05_config_geo_utils.R", c("load_step05_macroecology_config", "compute_haversine_dist_km_matrix"))
.source_step05_mod("scripts/analysis/step05_layer_hygiene_utils.R", "step05_layer_hygiene_utils.R", c("assemble_mem_decoupled_importance_df", "is_fine_edaphic_type_or_anchor_eligible"))
.source_step05_mod("scripts/analysis/step05_env_outlier_utils.R", "step05_env_outlier_utils.R", c("screen_within_taxon_env_mcd_outliers", "detect_r_centroid_or_inst_quarantine"))
.source_step05_mod("scripts/analysis/edaphic_geochem_wetland_utils.R", "edaphic_geochem_wetland_utils.R", c("enrich_geochem_wetland_predictors"), optional = TRUE)
.source_step05_mod("scripts/analysis/edaphic_fine_grain_utils.R", "edaphic_fine_grain_utils.R", c("flag_coordinate_uncertainty", "extract_uncertainty_buffered_predictors", "compute_derived_habitat_predictors"))
.source_step05_mod("scripts/analysis/edaphic_ssurgo_sda_utils.R", "edaphic_ssurgo_sda_utils.R", c("query_ssurgo_sda", "classify_ssurgo_outcrop", "aggregate_ssurgo_tabular_to_points"), optional = TRUE)
.source_step05_mod("scripts/analysis/edaphic_solus_3dep_utils.R", "edaphic_solus_3dep_utils.R", c("extract_polaris_and_solus_predictors", "extract_fine_grain_edaphic_metrics"), optional = TRUE)
.source_step05_mod("scripts/analysis/edaphic_collinearity_contingency_utils.R", "edaphic_collinearity_contingency_utils.R", c("screen_edaphic_point_mass_columns", "evaluate_multiclass_substrate_contingency"), optional = TRUE)
.source_step05_mod("scripts/analysis/step05_circular_hopkins_models.R", "step05_circular_hopkins_models.R", c("run_circular_multigroup_phase_test"), optional = TRUE)
.source_step05_mod("scripts/analysis/step05_allochronic_pheno_utils.R", "step05_allochronic_pheno_utils.R", c("test_allochronic_bimodality_von_mises", "compute_phenological_anomalies"), optional = TRUE)
.source_step05_mod("scripts/analysis/step05_type_niche_utils.R", "step05_type_niche_utils.R", c("compute_type_locality_ecological_control"), optional = TRUE)
.source_step05_mod("scripts/analysis/step05_phenology_ecotypes.R", "step05_phenology_ecotypes.R", c("resolve_analysis_group", "resolve_active_anchor_groups", "compute_ibe_partial_rda"))

WC_BIO_LOOKUP <- c(
  "1" = "bio1_temp_mean", "2" = "bio2_mean_diurnal_range", "3" = "bio3_isothermality",
  "4" = "bio4_temp_seasonality", "5" = "bio5_max_temp_warmest_month", "6" = "bio6_min_temp_coldest_month",
  "7" = "bio7_temp_annual_range", "8" = "bio8_mean_temp_wettest_qtr", "9" = "bio9_mean_temp_driest_qtr",
  "10" = "bio10_mean_temp_warmest_qtr", "11" = "bio11_mean_temp_coldest_qtr", "12" = "bio12_precip_annual",
  "13" = "bio13_precip_wettest_month", "14" = "bio14_precip_driest_month", "15" = "bio15_precip_seasonality",
  "16" = "bio16_precip_wettest_qtr", "17" = "bio17_precip_driest_qtr", "18" = "bio18_precip_warmest_qtr",
  "19" = "bio19_precip_coldest_qtr"
)

parse_args_robust <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  run_dir <- trimws(env_plantmatics_or_packera("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR"))
  rp <- function(rel) if (nzchar(run_dir)) file.path(run_dir, rel) else rel
  defs <- list(
    vouchers = rp("data/tables/curated_vouchers.csv"), morphometrics = rp("data/tables/morphometric_flags.csv"),
    pheno_states = rp("data/tables/voucher_phenological_states.csv"), env_dir = "data/environmental",
    output_flags = rp("data/tables/multimodal_conflict_flags.csv"), output_plot = rp("outputs/figures/fig05a_ecological_niche_and_rf_importance.png"),
    output_summary = rp("outputs/reports/multimodal_spatial_rf_summary.csv"), pheno_anomalies = rp("data/tables/phenological_anomalies.csv"),
    pheno_summary = rp("outputs/reports/phenological_anomaly_summary.csv"), pheno_plot = rp("outputs/figures/fig05b_hopkins_circular_phenology.png"),
    ssurgo_cache = "data/environmental/ssurgo_cache", ssurgo_table = rp("data/tables/ssurgo_edaphic_validation.csv"),
    contingency_table = rp("outputs/reports/micro_edaphic_outcrop_contingency.csv"), batch_size = 75L,
    permutations = 100L, n_trees = 500L, intraspecific_mode = FALSE, focal_taxon = "", demo = FALSE
  )
  raw_args <- sub("^--batch-size($|=)", "--ssurgo-batch-size\\1", cli_args)
  if (requireNamespace("optparse", quietly = TRUE)) {
    ol <- list(
      optparse::make_option(c("-v", "--vouchers"), type = "character", default = defs$vouchers),
      optparse::make_option(c("-m", "--morphometrics"), type = "character", default = defs$morphometrics),
      optparse::make_option(c("--pheno-states"), type = "character", default = defs$pheno_states),
      optparse::make_option(c("-e", "--env-dir"), type = "character", default = defs$env_dir),
      optparse::make_option(c("-f", "--output-flags"), type = "character", default = defs$output_flags),
      optparse::make_option(c("-p", "--output-plot"), type = "character", default = defs$output_plot),
      optparse::make_option(c("-s", "--output-summary"), type = "character", default = defs$output_summary),
      optparse::make_option(c("--pheno-anomalies"), type = "character", default = defs$pheno_anomalies),
      optparse::make_option(c("--pheno-summary"), type = "character", default = defs$pheno_summary),
      optparse::make_option(c("--pheno-plot"), type = "character", default = defs$pheno_plot),
      optparse::make_option(c("--ssurgo-cache"), type = "character", default = defs$ssurgo_cache),
      optparse::make_option(c("--ssurgo-table"), type = "character", default = defs$ssurgo_table),
      optparse::make_option(c("--contingency-table"), type = "character", default = defs$contingency_table),
      optparse::make_option(c("--ssurgo-batch-size"), dest = "batch_size", type = "integer", default = defs$batch_size),
      optparse::make_option(c("-k", "--permutations"), type = "integer", default = defs$permutations),
      optparse::make_option(c("-t", "--n-trees"), type = "integer", default = defs$n_trees),
      optparse::make_option(c("--intraspecific-mode"), action = "store_true", default = FALSE),
      optparse::make_option(c("--focal-taxon"), type = "character", default = ""),
      optparse::make_option(c("--demo"), action = "store_true", default = FALSE)
    )
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = ol), args = raw_args)
    names(opts) <- gsub("-", "_", names(opts)); return(opts)
  }
  fmap <- c("-v" = "vouchers", "--vouchers" = "vouchers", "-m" = "morphometrics", "--morphometrics" = "morphometrics", "--pheno-states" = "pheno_states",
            "-e" = "env_dir", "--env-dir" = "env_dir", "-f" = "output_flags", "--output-flags" = "output_flags", "-p" = "output_plot", "--output-plot" = "output_plot",
            "-s" = "output_summary", "--output-summary" = "output_summary", "--pheno-anomalies" = "pheno_anomalies", "--pheno-summary" = "pheno_summary",
            "--pheno-plot" = "pheno_plot", "--ssurgo-cache" = "ssurgo_cache", "--ssurgo-table" = "ssurgo_table", "--contingency-table" = "contingency_table",
            "--ssurgo-batch-size" = "batch_size", "-k" = "permutations", "--permutations" = "permutations", "-t" = "n_trees", "--n-trees" = "n_trees", "--focal-taxon" = "focal_taxon")
  i <- 1L
  while (i <= length(raw_args)) {
    a <- raw_args[i]
    if (a == "--demo") { defs$demo <- TRUE; i <- i + 1L; next }
    if (a == "--intraspecific-mode") { defs$intraspecific_mode <- TRUE; i <- i + 1L; next }
    if (grepl("^--[^=]+=", a)) {
      k <- sub("=.*$", "", a); v <- sub("^[^=]+=", "", a)
      if (k %in% names(fmap)) { d <- fmap[[k]]; defs[[d]] <- if (d %in% c("batch_size", "permutations", "n_trees")) as.integer(v) else v }
      i <- i + 1L; next
    }
    if (a %in% names(fmap) && i < length(raw_args)) { d <- fmap[[a]]; v <- raw_args[i + 1L]; defs[[d]] <- if (d %in% c("batch_size", "permutations", "n_trees")) as.integer(v) else v; i <- i + 2L; next }
    i <- i + 1L
  }
  defs
}

get_active_packera_taxa <- function(species_vec) {
  u_sp <- sort(unique(as.character(stats::na.omit(species_vec))))
  u_p <- u_sp[nzchar(u_sp) & u_sp != "Unknown" & grepl("^Packera[[:space:]]+", u_sp)]
  if (length(u_p) > 0L) u_p else u_sp[nzchar(u_sp) & u_sp != "Unknown"]
}

is_canonical_prod_path <- function(p) {
  if (is.null(p) || !nzchar(p)) return(FALSE)
  np <- normalizePath(p, mustWork = FALSE)
  startsWith(np, paste0(normalizePath("data", mustWork = FALSE), "/")) || startsWith(np, paste0(normalizePath("outputs", mustWork = FALSE), "/"))
}

.bind_rows_schema_safe <- function(df1, df2) {
  if (is.null(df1) || nrow(df1) == 0L) return(df2)
  if (is.null(df2) || nrow(df2) == 0L) return(df1)
  lc <- if (exists("PACKERA_DB_LOGICAL_COLS", inherits = TRUE)) get("PACKERA_DB_LOGICAL_COLS", inherits = TRUE) else character(0)
  .tb <- function(x) { if (is.logical(x)) return(x); if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), NA, x != 0)); xu <- toupper(trimws(as.character(x))); ifelse(is.na(x) | !nzchar(xu) | xu == "NA", NA, xu %in% c("TRUE", "T", "1", "YES")) }
  for (cn in intersect(names(df1), names(df2))) {
    v1 <- df1[[cn]]; v2 <- df2[[cn]]; if (identical(typeof(v1), typeof(v2)) && identical(class(v1), class(v2))) next
    if (cn %in% lc) {
      df1[[cn]] <- .tb(v1); df2[[cn]] <- .tb(v2)
    } else if (is.logical(v1) && all(is.na(v1)) && !is.logical(v2)) {
      df1[[cn]] <- if (is.character(v2) || is.factor(v2)) rep(NA_character_, nrow(df1)) else if (is.integer(v2)) rep(NA_integer_, nrow(df1)) else rep(NA_real_, nrow(df1))
    } else if (is.logical(v2) && all(is.na(v2)) && !is.logical(v1)) {
      df2[[cn]] <- if (is.character(v1) || is.factor(v1)) rep(NA_character_, nrow(df2)) else if (is.integer(v1)) rep(NA_integer_, nrow(df2)) else rep(NA_real_, nrow(df2))
    } else if (is.logical(v1) && (is.numeric(v2) || is.integer(v2)) && all(is.na(v2) | v2 %in% c(0, 1))) {
      df2[[cn]] <- .tb(v2)
    } else if (is.logical(v2) && (is.numeric(v1) || is.integer(v1)) && all(is.na(v1) | v1 %in% c(0, 1))) {
      df1[[cn]] <- .tb(v1)
    } else if ((is.numeric(v1) || is.integer(v1)) && (is.numeric(v2) || is.integer(v2)) && !is.factor(v1) && !is.factor(v2)) {
      df1[[cn]] <- suppressWarnings(as.numeric(v1)); df2[[cn]] <- suppressWarnings(as.numeric(v2))
    } else {
      df1[[cn]] <- as.character(v1); df2[[cn]] <- as.character(v2)
    }
  }
  dplyr::bind_rows(df1, df2)
}

write_csv_atomic <- function(df, target_path, demo = FALSE) {
  if (isTRUE(demo) && is_canonical_prod_path(target_path)) { message(sprintf("[DEMO MODE] Skipping write to canonical path: %s", target_path)); return(invisible(FALSE)) }
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(target_path, ".tmp")
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(df, file = tmp, na = "", quote = "auto", logical01 = FALSE) else readr::write_csv(df, tmp)
  file.rename(tmp, target_path); invisible(TRUE)
}

.split_step05_cv_domain_vars <- function(cv_cand_vars) {
  clim_cand_cv <- intersect(grep("^(bio[0-9]*_|summer_moisture_|maritime_continentality_|delta_doy|doy|hopkins_)", cv_cand_vars, value = TRUE, ignore.case = TRUE), cv_cand_vars)
  topo_cand_cv <- setdiff(intersect(grep("^(topo_|elevation($|_)|elev_)", cv_cand_vars, value = TRUE, ignore.case = TRUE), cv_cand_vars), clim_cand_cv)
  edaph_cand_cv <- setdiff(cv_cand_vars, unique(c(clim_cand_cv, topo_cand_cv)))
  list(clim_vars = clim_cand_cv, topo_vars = topo_cand_cv, edaph_vars = edaph_cand_cv)
}

synthesize_multiscale_edaphics <- function(df, ssurgo_out = "data/tables/ssurgo_edaphic_validation.csv",
                                           contingency_out = "outputs/reports/micro_edaphic_outcrop_contingency.csv", demo = FALSE) {
  message("Executing Multi-Scale Edaphic Synthesis & Fisher's Exact Outcrop Fidelity Test...")
  val_cols <- c(
    "catalogNumber", "plant_individual_id", "lat_lon_key", "latitude", "longitude", "species_standardized", "regional_group",
    "coordinateUncertaintyInMeters", "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
    "is_clim_eligible_coord", "is_edaphic_eligible_coord", "edaphic_coord_weight", "musym", "muname", "taxsuborder",
    "taxgreatgroup", "is_rock_outcrop", "outcrop_class", "ssurgo_rock_outcrop_pct", "ssurgo_lithic_series_pct",
    "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec",
    "ssurgo_aws_0_50cm", "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m",
    "polaris_ksat_30m", "solus_resdept_cm", "solus_fragvol_pct", "solus_sand_pct", "solus_clay_pct", "solus_silt_pct",
    "solus_cec", "elevation_m", "topo_slope_deg", "topo_tpi_30m", "topo_hli_30m", "topo_twi_30m", "fine_soil_ph",
    "fine_soil_sand", "fine_soil_clay", "soil_silt_pct", "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct",
    "fine_bedrock_depth_cm", "sgmc_lithology_class", "is_azonal_lithology", "soil_ca_mg_molar_ratio",
    "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct", "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m",
    "ultramafic_calcareous_gradient", "wetland_seepage_index", "soil_ph", "soil_cec", "soil_sand", "soil_clay", "soil_silt",
    "soil_cfvo", "soil_bulk_density"
  )
  val_df <- df[, intersect(val_cols, names(df)), drop = FALSE]
  sg_b <- c("soil_ph", "soil_cec", "soil_sand", "soil_clay", "soil_silt", "soil_cfvo", "soil_bulk_density")
  names(val_df)[names(val_df) %in% sg_b] <- paste0(names(val_df)[names(val_df) %in% sg_b], "_soilgrids")
  write_csv_atomic(val_df, ssurgo_out, demo = demo)
  active_tx <- get_active_packera_taxa(df$species_standardized); grp_col <- "species_standardized"; grp_lvls <- active_tx
  if (length(active_tx) < 2L && "analysis_group" %in% names(df)) {
    ag <- sort(unique(as.character(stats::na.omit(df$analysis_group)))); ag <- ag[nzchar(ag) & ag != "Unknown"]
    if (length(ag) >= 2L) { grp_col <- "analysis_group"; grp_lvls <- ag }
  }
  if (!("is_rock_outcrop" %in% names(df))) df$is_rock_outcrop <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) !is.na(df$ssurgo_rock_outcrop_pct) & as.numeric(df$ssurgo_rock_outcrop_pct) >= 15.0 else rep(FALSE, nrow(df))
  df$is_rock_outcrop[is.na(df$is_rock_outcrop)] <- FALSE
  mc <- evaluate_multiclass_substrate_contingency(df, active_groups = grp_lvls, group_col = grp_col)
  if (nrow(mc$summary_table) > 0L) write_csv_atomic(mc$summary_table, contingency_out, demo = demo)
  list(validation = val_df, contingency = mc$summary_table, multiclass_contingency = mc$multiclass_details, fisher = mc$omnibus_test)
}

left_join_by_voucher_individual <- function(base_df, mod_df, keep_cols) {
  if (is.null(mod_df) || nrow(mod_df) == 0L || !("catalogNumber" %in% names(mod_df))) return(base_df)
  base_df$catalogNumber <- as.character(base_df$catalogNumber); mod_df$catalogNumber <- as.character(mod_df$catalogNumber)
  kc <- intersect(names(mod_df), keep_cols); vc <- setdiff(kc, c("catalogNumber", "plant_individual_id"))
  ovc <- intersect(vc, names(base_df)); adc <- setdiff(vc, ovc); hmp <- "plant_individual_id" %in% kc
  if (hmp) {
    mp <- ifelse(is.na(mod_df$plant_individual_id), 0L, suppressWarnings(as.integer(mod_df$plant_individual_id))); mp[is.na(mp)] <- 0L
    mod_df$plant_individual_id <- mp; msub <- mod_df[!duplicated(paste0(mod_df$catalogNumber, "__", mp)), kc, drop = FALSE]
  } else msub <- mod_df[!duplicated(mod_df$catalogNumber), kc, drop = FALSE]
  if (("plant_individual_id" %in% names(base_df)) && hmp) {
    bp <- ifelse(is.na(base_df$plant_individual_id), 0L, suppressWarnings(as.integer(base_df$plant_individual_id))); bp[is.na(bp)] <- 0L; base_df$plant_individual_id <- bp
    bc <- table(base_df$catalogNumber); s0 <- unique(base_df$catalogNumber[bp == 0L & as.integer(bc[base_df$catalogNumber]) == 1L])
    mnz <- msub[msub$catalogNumber %in% s0 & msub$plant_individual_id != 0L, c("catalogNumber", "plant_individual_id"), drop = FALSE]
    if (nrow(mnz) > 0L) { kb <- !(base_df$catalogNumber %in% unique(mnz$catalogNumber)); eb <- dplyr::left_join(base_df[!kb, setdiff(names(base_df), "plant_individual_id"), drop = FALSE], mnz, by = "catalogNumber"); base_df <- .bind_rows_schema_safe(base_df[kb, , drop = FALSE], eb[, names(base_df), drop = FALSE]) }
    out <- dplyr::left_join(base_df, msub[, unique(c("catalogNumber", "plant_individual_id", adc)), drop = FALSE], by = c("catalogNumber", "plant_individual_id"))
    mc <- table(msub$catalogNumber); sh <- msub[msub$plant_individual_id == 0L & as.integer(mc[msub$catalogNumber]) == 1L, , drop = FALSE]
    iex <- match(paste0(out$catalogNumber, "__", out$plant_individual_id), paste0(msub$catalogNumber, "__", msub$plant_individual_id))
    ish <- if (nrow(sh) > 0L) match(out$catalogNumber, sh$catalogNumber) else rep(NA_integer_, nrow(out))
    bc2 <- table(out$catalogNumber); isc <- ifelse(out$catalogNumber %in% intersect(names(bc2[bc2 == 1L]), names(mc[mc == 1L])), match(out$catalogNumber, msub$catalogNumber), NA_integer_)
    for (cl in adc) { nm <- is.na(out[[cl]]) & is.na(iex); if (any(nm & !is.na(ish))) { h <- nm & !is.na(ish); out[[cl]][h] <- sh[[cl]][ish[h]] }; h2 <- is.na(out[[cl]]) & is.na(iex) & !is.na(isc); if (any(h2)) out[[cl]][h2] <- msub[[cl]][isc[h2]] }
    for (cl in ovc) { h1 <- is.na(out[[cl]]) & !is.na(iex); if (any(h1)) out[[cl]][h1] <- msub[[cl]][iex[h1]]; h2 <- is.na(out[[cl]]) & !is.na(ish); if (any(h2)) out[[cl]][h2] <- sh[[cl]][ish[h2]]; h3 <- is.na(out[[cl]]) & !is.na(isc); if (any(h3)) out[[cl]][h3] <- msub[[cl]][isc[h3]] }
    return(out)
  }
  jc <- unique(c("catalogNumber", if (hmp && !("plant_individual_id" %in% names(base_df))) "plant_individual_id" else character(0), adc))
  out <- dplyr::left_join(base_df, msub[, jc, drop = FALSE], by = "catalogNumber")
  if (length(ovc) > 0L) {
    ic <- if (hmp && "plant_individual_id" %in% names(out)) match(paste0(out$catalogNumber, "__", out$plant_individual_id), paste0(msub$catalogNumber, "__", msub$plant_individual_id)) else match(out$catalogNumber, msub$catalogNumber)
    for (cl in ovc) { h <- is.na(out[[cl]]) & !is.na(ic); if (any(h)) out[[cl]][h] <- msub[[cl]][ic[h]] }
  }
  out
}

.source_step05_mod("scripts/analysis/step05_env_raster_extract_utils.R", "step05_env_raster_extract_utils.R", c("discover_candidate_eco_predictors", "extract_environmental_layers"))
.source_step05_mod("scripts/analysis/step05_spatial_rf_cv_utils.R", "step05_spatial_rf_cv_utils.R", c("compute_spatial_rf_mems", "estimate_spatial_autocorrelation_block_size"))
.source_step05_mod("scripts/analysis/step05_consensus_niche_utils.R", "step05_consensus_niche_utils.R", c("execute_crossmodal_consensus", "run_warren_niche_identity_tests"))
.source_step05_mod("scripts/analysis/step05_sdm_engine_utils.R", "step05_sdm_engine_utils.R", c(".build_regional_sdm_raster_stack"), optional = TRUE)
.source_step05_mod("scripts/analysis/step05_sdm_maxent.R", "step05_sdm_maxent.R", c("fit_and_project_maxent_sdm_maps", "export_spatial_rf_figures"))

main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  set.seed(42L); opts <- parse_args_robust(cli_args)
  cfg <- if (exists("load_step05_macroecology_config", mode = "function")) load_step05_macroecology_config() else list(cor_threshold = 0.70, vif_threshold = 5.0, cumulative_importance_cutoff = 0.85)
  message("=== Starting Multimodal Spatial RF & Ecological Niche Analysis ===")
  if (!file.exists(opts$vouchers)) stop("Missing vouchers file: ", opts$vouchers)
  vdf <- readr::read_csv(opts$vouchers, show_col_types = FALSE); vdf$catalogNumber <- as.character(vdf$catalogNumber)
  sp_c <- rep("Unknown", nrow(vdf))
  for (sc in intersect(c("species_standardized", "scientificName", "species_raw"), names(vdf))) { v <- standardize_packera_taxon(vdf[[sc]]); m <- sp_c == "Unknown" & !is.na(v) & v != "Unknown"; sp_c[m] <- v[m] }
  vdf$species_standardized <- sp_c

  # Modality-Independent SQLite Hydration (zero Step 04 morphometric columns or Step 04 authority_weight)
  if (!isTRUE(opts$demo) && exists("connect_packera_db", mode = "function") && exists("is_sqlite_active_for_inputs", mode = "function") && is_sqlite_active_for_inputs(opts$vouchers, demo = isTRUE(opts$demo))) {
    con <- tryCatch(connect_packera_db(read_only = TRUE, demo = FALSE), error = function(e) NULL)
    if (!is.null(con)) {
      rid <- resolve_packera_run_id()
      rsel <- intersect(c("catalogNumber", "plant_individual_id", "species_standardized", "determiner_tier", "type_status", "is_type_specimen", "is_promoted_bronze", "spatial_balance_weight"), tryCatch(DBI::dbListFields(con, "run_active_cohort"), error = function(e) character(0)))
      rrows <- if ("catalogNumber" %in% rsel) tryCatch(DBI::dbGetQuery(con, sprintf("SELECT %s FROM run_active_cohort WHERE run_id = $1 ORDER BY catalogNumber%s;", paste(rsel, collapse = ", "), if ("plant_individual_id" %in% rsel) ", plant_individual_id" else ""), params = list(rid)), error = function(e) NULL) else NULL
      vsel <- intersect(c("typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status", "identification_year", "is_open_or_hybrid_concept", "is_centroid_or_inst_coord", "coord_quarantine_reason", "coordinateUncertainty"), tryCatch(DBI::dbListFields(con, "vouchers"), error = function(e) character(0)))
      vrows <- if (length(vsel) > 0L) tryCatch(DBI::dbGetQuery(con, sprintf("SELECT catalogNumber, %s FROM vouchers;", paste(vsel, collapse = ", "))), error = function(e) NULL) else NULL
      disconnect_packera_db(con)
      if (!is.null(rrows) && nrow(rrows) > 0L) { for (bc in intersect(c("is_type_specimen", "is_promoted_bronze"), names(rrows))) rrows[[bc]] <- .coerce_bool_flag_vec(rrows[[bc]], nrow(rrows)); vdf <- left_join_by_voucher_individual(vdf, rrows, rsel) }
      if (!is.null(vrows) && nrow(vrows) > 0L) {
        vrows$catalogNumber <- as.character(vrows$catalogNumber); vrows <- vrows[!duplicated(vrows$catalogNumber), , drop = FALSE]
        for (bc in intersect(c("is_centroid_or_inst_coord", "is_open_or_hybrid_concept"), names(vrows))) vrows[[bc]] <- .coerce_bool_flag_vec(vrows[[bc]], nrow(vrows))
        vdf <- left_join_by_voucher_individual(vdf, vrows, c("catalogNumber", vsel))
        if ("coordinateUncertainty" %in% names(vdf)) { u1 <- if ("coordinateUncertaintyInMeters" %in% names(vdf)) suppressWarnings(as.numeric(vdf$coordinateUncertaintyInMeters)) else rep(NA_real_, nrow(vdf)); u2 <- suppressWarnings(as.numeric(vdf$coordinateUncertainty)); vdf$coordinateUncertaintyInMeters <- ifelse(is.finite(u1), u1, u2) }
      }
    }
  }
  if (file.exists(opts$pheno_states)) {
    pdf <- readr::read_csv(opts$pheno_states, show_col_types = FALSE); pdf$catalogNumber <- as.character(pdf$catalogNumber)
    pc <- intersect(names(pdf), c("catalogNumber", "total_capitula", "n_anthesis", "n_fruiting", "n_bud", "n_senescent", "voucher_phenological_state", "reproductive_stage_confidence", "dominant_pappus_ratio", "corymb_head_count_est", "radiate_head_fraction", "anthesis_head_count", "fruiting_head_count", "Flower_MANY", "Seed_Fruit_MANY", "is_verified_anthesis"))
    vdf <- left_join_by_voucher_individual(vdf, pdf[!duplicated(pdf$catalogNumber), pc, drop = FALSE], pc)
  } else {
    sp_sql <- if (exists("fetch_step02_phenology_from_sqlite", mode = "function")) fetch_step02_phenology_from_sqlite(vdf$catalogNumber, opts$vouchers, demo = isTRUE(opts$demo)) else NULL
    if (is.data.frame(sp_sql) && nrow(sp_sql) > 0L) vdf <- left_join_by_voucher_individual(vdf, sp_sql, names(sp_sql))
    if (!("voucher_phenological_state" %in% names(vdf))) vdf$voucher_phenological_state <- NA_character_
    for (rc in c("total_capitula", "n_anthesis", "n_fruiting", "n_bud", "dominant_pappus_ratio")) if (!(rc %in% names(vdf))) vdf[[rc]] <- NA_real_
  }

  pm_df <- NULL
  if (!is.null(opts$morphometrics) && nzchar(trimws(opts$morphometrics)) && opts$morphometrics != "data/tables/morphometric_flags.csv" && file.exists(opts$morphometrics)) {
    pm_df <- tryCatch(readr::read_csv(opts$morphometrics, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(pm_df) && "catalogNumber" %in% names(pm_df)) {
      pm_df$catalogNumber <- as.character(pm_df$catalogNumber)
      if ("cda_predicted" %in% names(pm_df) && !("cda_predicted_taxon" %in% names(pm_df))) pm_df$cda_predicted_taxon <- pm_df$cda_predicted
      if ("is_outlier" %in% names(pm_df) && !("morphometric_outlier_flag" %in% names(pm_df))) pm_df$morphometric_outlier_flag <- as.logical(pm_df$is_outlier)
      if ("plant_individual_id" %in% names(pm_df)) vdf <- left_join_by_voucher_individual(vdf, pm_df[, c("catalogNumber", "plant_individual_id"), drop = FALSE], c("catalogNumber", "plant_individual_id"))
    }
  }
  if ("is_primary_duplicate" %in% names(vdf)) vdf <- vdf[is.na(vdf$is_primary_duplicate) | toupper(as.character(vdf$is_primary_duplicate)) %in% c("TRUE", "T", "1"), , drop = FALSE]
  for (px in list(c("latitude", "decimalLatitude"), c("longitude", "decimalLongitude"))) {
    if (!(px[1L] %in% names(vdf))) vdf[[px[1L]]] <- NA_real_; v1 <- suppressWarnings(as.numeric(vdf[[px[1L]]]))
    if (px[2L] %in% names(vdf)) { v2 <- suppressWarnings(as.numeric(vdf[[px[2L]]])); v1 <- ifelse(is.finite(v1), v1, v2) }
    vdf[[px[1L]]] <- v1
  }
  hgeo <- is.finite(vdf$latitude) & is.finite(vdf$longitude) & abs(vdf$latitude) <= 90.0 & abs(vdf$longitude) <= 180.0
  ng_df <- vdf[!hgeo, , drop = FALSE]; vdf <- vdf[hgeo, , drop = FALSE]
  if (nrow(vdf) == 0L) stop("Zero georeferenced vouchers available for Step 05 spatial ecological modeling.")
  dm <- isTRUE(opts$demo); skip_f <- identical(trimws(env_plantmatics_or_packera("PLANTMATICS_SKIP_FIGURES", "PACKERA_SKIP_FIGURES")), "1")

  vdf <- extract_environmental_layers(vdf, opts$env_dir, demo = dm)
  qsda <- if (exists("query_ssurgo_sda", mode = "function")) query_ssurgo_sda else NULL
  if (exists("load_and_update_fine_grain_cache", mode = "function")) {
    vdf <- load_and_update_fine_grain_cache(vdf, env_dir = opts$env_dir, ssurgo_cache_dir = opts$ssurgo_cache, batch_size = opts$batch_size, demo = dm, query_ssurgo_fn = qsda, cache_path = file.path(opts$env_dir, "fine_grain_edaphic_cache.csv"))
  } else if (!is.null(qsda)) {
    vdf <- qsda(vdf, cache_dir = opts$ssurgo_cache, batch_size = opts$batch_size, demo = dm)
    if (exists("extract_sgmc_lithology", mode = "function")) vdf <- extract_sgmc_lithology(vdf, env_dir = opts$env_dir, demo = dm)
    if (exists("extract_uncertainty_buffered_predictors", mode = "function")) vdf <- extract_uncertainty_buffered_predictors(vdf, env_dir = opts$env_dir, demo = dm)
  }
  synthesize_multiscale_edaphics(vdf, ssurgo_out = opts$ssurgo_table, contingency_out = opts$contingency_table, demo = dm)
  if (exists("compute_derived_habitat_predictors", mode = "function")) {
    gwc <- if (!dm && exists("connect_packera_db", mode = "function") && (nzchar(Sys.getenv("RUN_ID", "")) || nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")) || (exists("is_sqlite_active_for_inputs", mode = "function") && is_sqlite_active_for_inputs(opts$vouchers, demo = dm)))) tryCatch(connect_packera_db(read_only = TRUE, demo = FALSE), error = function(e) NULL) else NULL
    vdf <- compute_derived_habitat_predictors(vdf, db_conn = gwc, env_dir = opts$env_dir)
    if (!is.null(gwc) && exists("disconnect_packera_db", mode = "function")) tryCatch(disconnect_packera_db(gwc), error = function(e) NULL)
  }

  pheno_res <- compute_phenological_anomalies(vdf, output_csv = opts$pheno_anomalies, summary_csv = opts$pheno_summary, demo = dm)
  vdf <- pheno_res$df; if (!skip_f) export_phenological_latitudinal_figures(vdf, pheno_res$baseline, opts$pheno_plot)
  cand_eco <- discover_candidate_eco_predictors(vdf, include_pheno = FALSE, min_finite = 3L)
  for (ev in cand_eco) if (ev %in% names(vdf)) vdf[[ev]] <- suppressWarnings(as.numeric(vdf[[ev]]))
  is_intra <- isTRUE(opts$intraspecific_mode) || (length(unique(stats::na.omit(vdf$species_standardized))) < 2L)
  if (is_intra) vdf <- resolve_analysis_group(vdf, intraspecific_mode = TRUE) else vdf$analysis_group <- vdf$species_standardized
  vdf$ecotype_regime <- as.character(vdf$analysis_group)

  srf <- compute_spatial_rf_mems(vdf, cand_eco, n_trees = opts$n_trees, intraspecific_mode = is_intra, vif_res = NULL, env_dir = opts$env_dir, demo = dm)
  if (is_intra && (isTRUE(srf$ecotype_collapsed_to_unimodal) || !isTRUE(srf$is_post_selection_ecotype))) {
    sb <- stats::na.omit(trimws(as.character(vdf$species_standardized)))
    vdf$analysis_group <- rep(if (length(sb) > 0L && nzchar(sb[1L])) sb[1L] else "Intraspecific_Lineage", nrow(vdf))
    vdf$ecotype_regime <- vdf$analysis_group; attr(vdf, "is_post_selection_ecotype") <- FALSE
  }
  vif_res <- if (is.list(srf$collinearity_screening)) srf$collinearity_screening else NULL
  kept_rf <- if (!is.null(vif_res) && length(vif_res$kept_vars) > 0L) vif_res$kept_vars else cand_eco
  vdf$spatial_block_cv_acc <- round(srf$Spatial_Block_CV_Accuracy, 4); vdf$spatial_autocorr_range_km <- srf$spatial_autocorr_range_km
  vdf$rf_predicted_taxon <- srf$rf_predicted_taxon; vdf$rf_predicted_prob <- srf$rf_predicted_prob
  vdf$rf_given_taxon_prob <- srf$rf_given_taxon_prob; vdf$rf_prediction_source <- srf$rf_prediction_source; vdf$spatial_block_id <- srf$spatial_block_id
  vdf$env_mahal_d2 <- if (length(srf$env_mahal_d2) == nrow(vdf)) srf$env_mahal_d2 else rep(NA_real_, nrow(vdf))
  vdf$env_mcd_outlier_flag <- if (length(srf$env_mcd_outlier_flag) == nrow(vdf)) as.logical(srf$env_mcd_outlier_flag) else rep(FALSE, nrow(vdf))

  psel <- select_primary_predictors_from_rf(srf, vif_res, vdf, cum_imp_target = cfg$cumulative_importance_cutoff, intraspecific_mode = is_intra)
  if (!is.null(pm_df) && nrow(pm_df) > 0L) vdf <- left_join_by_voucher_individual(vdf, pm_df, intersect(names(pm_df), c("catalogNumber", "plant_individual_id", "cda_predicted_taxon", "cda_posterior_prob", "cda_posterior_margin", "morphometric_outlier_flag")))
  vdf <- execute_crossmodal_consensus(vdf, is_intra, srf, psel$primary_env_vars, psel$primary_clim_vars, psel$primary_edaph_vars, psel$rf_var_weights)
  vdf$ecotype_regime <- as.character(vdf$analysis_group); vdf$primary_env_predictors <- paste(psel$primary_env_vars, collapse = ";")
  vdf$primary_clim_predictors <- paste(psel$primary_clim_vars, collapse = ";"); vdf$primary_edaph_predictors <- paste(psel$primary_edaph_vars, collapse = ";")

  niche <- run_warren_niche_identity_tests(vdf, kept_rf, opts$permutations, is_intra, psel$primary_env_vars, psel$primary_clim_vars, psel$primary_edaph_vars, psel$rf_var_weights)
  edc <- c("Pair", "Taxon1", "Taxon2", "D_edaph", "I_edaph", "p_equiv_D_edaph", "p_bg_D_edaph")
  niche$edaphic_schoener_d <- if (length(psel$primary_edaph_vars) > 0L && is.data.frame(niche$summary) && all(edc %in% names(niche$summary))) niche$summary[, edc, drop = FALSE] else compute_edaphic_schoener_d(vdf, psel$primary_edaph_vars, psel$rf_var_weights, opts$permutations, is_intra)

  geo_vdf <- vdf
  if (nrow(ng_df) > 0L) {
    ng_df$analysis_group <- if (is_intra) rep(as.character(ng_df$species_standardized[1L]), nrow(ng_df)) else ng_df$species_standardized
    ng_df$ecotype_regime <- NA_character_; ng_df$rf_prediction_source <- "Passive_Ungeoreferenced"; ng_df$spatial_block_id <- NA_character_
    ng_df$high_coord_uncertainty_flag <- TRUE; ng_df$is_training_eligible_coord <- FALSE; ng_df$is_clim_eligible_coord <- FALSE
    ng_df$is_edaphic_eligible_coord <- FALSE; ng_df$edaphic_coord_weight <- 0.0; ng_df$env_mahal_d2 <- NA_real_; ng_df$env_mcd_outlier_flag <- FALSE
    ng_df <- execute_crossmodal_consensus(ng_df, is_intra, NULL, psel$primary_env_vars, psel$primary_clim_vars, psel$primary_edaph_vars, psel$rf_var_weights)
    vdf <- .bind_rows_schema_safe(geo_vdf, ng_df)
  }
  for (nc in c("hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic", "nearest_type_basionym_edaphic_z", "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx")) if (!(nc %in% names(vdf))) vdf[[nc]] <- rep(NA_real_, nrow(vdf))
  for (cc in c("nearest_type_niche_basionym", "ecotype_regime")) if (!(cc %in% names(vdf))) vdf[[cc]] <- rep(NA_character_, nrow(vdf))

  exp_df <- if (all(c("decimalLatitude", "decimalLongitude") %in% names(vdf))) vdf[, setdiff(names(vdf), c("latitude", "longitude")), drop = FALSE] else vdf
  write_csv_atomic(exp_df, opts$output_flags, demo = dm)
  if (!dm && exists("persist_step05_to_sqlite", mode = "function")) tryCatch(persist_step05_to_sqlite(vdf, opts$vouchers, demo = dm), error = function(e) warning("SQLite Step 05 persistence note: ", e$message))

  sdm_dir <- dirname(dirname(opts$output_plot)); if (is.null(sdm_dir) || !nzchar(sdm_dir) || identical(sdm_dir, ".")) sdm_dir <- "outputs"
  ftx <- if (!is.null(opts$focal_taxon) && nzchar(trimws(opts$focal_taxon))) standardize_packera_taxon(opts$focal_taxon) else "Packera dubia"
  sdm_res <- if (!skip_f) tryCatch(fit_and_project_maxent_sdm_maps(geo_vdf, srf, opts$env_dir, sdm_dir, kept_rf, psel$primary_clim_vars, psel$primary_edaph_vars, ftx, is_intra, dm), error = function(e) { warning("Step 05 SDM warning: ", e$message); NULL }) else NULL

  .pm <- function(pdf, cols) { if (!is.data.frame(pdf) || nrow(pdf) == 0L) return(NA_real_); for (cc in intersect(cols, names(pdf))) { v <- suppressWarnings(as.numeric(pdf[[cc]])); v <- v[is.finite(v)]; if (length(v) > 0L) return(round(mean(v), 4)) }; NA_real_ }
  oob_o <- if (!is.null(srf$OOB_Accuracy) && is.finite(srf$OOB_Accuracy)) round(as.numeric(srf$OOB_Accuracy), 4) else NA_real_
  sb_o <- if (!is.null(srf$Spatial_Block_CV_Accuracy) && is.finite(srf$Spatial_Block_CV_Accuracy)) round(as.numeric(srf$Spatial_Block_CV_Accuracy), 4) else NA_real_
  oob_gap <- if (is.finite(oob_o) && is.finite(sb_o)) round(oob_o - sb_o, 4) else NA_real_
  cvp <- if (!is.null(srf$cv_mean_fold_predictors) && is.finite(srf$cv_mean_fold_predictors)) round(as.numeric(srf$cv_mean_fold_predictors), 2) else round(as.numeric(length(kept_rf)), 2)
  cvj <- if (!is.null(srf$cv_fold_predictor_jaccard) && is.finite(srf$cv_fold_predictor_jaccard)) round(as.numeric(srf$cv_fold_predictor_jaccard), 4) else 1.0
  ns <- if (is.data.frame(niche$summary)) niche$summary else data.frame(Pair = character(0), Taxon1 = character(0), Taxon2 = character(0), Schoener_D = numeric(0), p_value_D = numeric(0), stringsAsFactors = FALSE)
  nrs <- nrow(ns); ns$CV_Mean_Fold_Predictors <- rep(cvp, nrs); ns$CV_Fold_Predictor_Jaccard <- rep(cvj, nrs)
  ns$OOB_Accuracy <- rep(oob_o, nrs); ns$Spatial_Block_CV_Accuracy <- rep(sb_o, nrs); ns$OOB_Block_Gap <- rep(oob_gap, nrs)
  ns$AUC_Block_Mean <- rep(if (is.list(sdm_res)) .pm(sdm_res$performance, c("AUC_Block_Ensemble", "AUC_Block", "AUC_Block_Mean", "AUC_Block_SpatialRF", "AUC_Block_MaxEnt")) else NA_real_, nrs)
  ns$CBI_Mean <- rep(if (is.list(sdm_res)) .pm(sdm_res$performance, c("CBI_Ensemble", "CBI", "CBI_Mean", "CBI_SpatialRF", "CBI_MaxEnt")) else NA_real_, nrs)
  if (is.list(niche)) niche$summary <- ns; write_csv_atomic(ns, opts$output_summary, demo = dm)

  if (!dm && exists("write_env_outlier_quarantine_audit", mode = "function")) {
    anc_g <- rep(FALSE, nrow(geo_vdf)); if (length(srf$anchor_idx) > 0L) anc_g[srf$anchor_idx[srf$anchor_idx >= 1L & srf$anchor_idx <= nrow(geo_vdf)]] <- TRUE
    acsv <- file.path(dirname(opts$output_summary), "environmental_outlier_quarantine_audit.csv")
    if (nrow(ng_df) > 0L) {
      ra <- if (length(srf$env_outlier_screen$env_outlier_rescued_azonal) == nrow(geo_vdf)) as.logical(srf$env_outlier_screen$env_outlier_rescued_azonal) else rep(FALSE, nrow(geo_vdf))
      write_env_outlier_quarantine_audit(vdf, list(env_mahal_d2 = vdf$env_mahal_d2, env_mcd_outlier_flag = vdf$env_mcd_outlier_flag, env_outlier_rescued_azonal = c(ra, rep(FALSE, nrow(ng_df)))), acsv, c(anc_g, rep(FALSE, nrow(ng_df))))
    } else write_env_outlier_quarantine_audit(vdf, srf$env_outlier_screen, acsv, anc_g)
  }
  if (is.data.frame(srf$ibe_summary) && nrow(srf$ibe_summary) > 0L) write_csv_atomic(srf$ibe_summary, file.path(dirname(opts$output_summary), "isolation_by_environment_summary.csv"), demo = dm)
  if (is.data.frame(srf$variance_partition) && nrow(srf$variance_partition) > 0L) write_csv_atomic(srf$variance_partition, file.path(dirname(opts$output_summary), "variance_partitioning_summary.csv"), demo = dm)
  if (!skip_f) export_spatial_rf_figures(geo_vdf, srf, niche, opts$output_plot, pheno_res = pheno_res)
  message("=== Multimodal Spatial RF & Niche Modeling Analysis Complete ===")
}

if (!interactive() && sys.nframe() == 0L) main()
