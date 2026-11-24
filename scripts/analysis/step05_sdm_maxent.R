if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/analysis/step05_sdm_maxent.R
# Step 05 Spatial RF Niche Figure Delegation + Zero-Leakage MaxEnt/Ensemble SDM
# ==============================================================================

.step05_sdm_source_sibling <- function(fname, check_sym, check_mode = "function") {
  if (exists(check_sym, mode = check_mode, inherits = TRUE)) return(invisible(FALSE))
  ofr <- vapply(sys.frames(), function(f) if (!is.null(f$ofile)) as.character(f$ofile) else "", character(1))
  oh <- ofr[nzchar(ofr)]; cf <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  cands <- c(if (length(oh) > 0L) file.path(dirname(utils::tail(oh, 1L)), fname) else character(0),
             if (length(cf) > 0L) file.path(dirname(cf[1L]), fname) else character(0),
             file.path("scripts", "analysis", fname), file.path(getwd(), "scripts", "analysis", fname), file.path("analysis", fname), fname)
  hit <- cands[file.exists(cands)][1L]
  if (!is.na(hit) && nzchar(hit)) { source(hit); return(invisible(TRUE)) }
  invisible(FALSE)
}

.step05_sdm_source_sibling("step05_config_geo_utils.R", "load_step05_macroecology_config", "function")
.step05_sdm_source_sibling("step05_layer_hygiene_utils.R", "STEP05_MICRO_SCALE_PREDICTORS", "any")
.step05_sdm_source_sibling("step05_sdm_engine_utils.R", ".fit_single_maxent_model", "function")

.ensure_step05_fig_modules <- function() {
  if (!exists("render_step05a_ecological_niche_png", mode = "function")) {
    for (.f in c("scripts/figures/theme_packera_pub.R", "scripts/figures/step05_figures.R")) if (file.exists(.f)) source(.f)
  }
}

export_spatial_rf_figures <- function(df, srf_res, niche_res, out_path, pheno_res = NULL) {
  .ensure_step05_fig_modules()
  out_png <- if (!is.null(out_path) && nzchar(out_path)) normalize_png_path(out_path) else "outputs/figures/fig05a_ecological_niche_and_rf_importance.png"
  message(sprintf("Generating 300-DPI Spatial RF & Ecological Niche plate -> %s", out_png))
  render_step05a_ecological_niche_png(df = df, srf_res = srf_res, niche_res = niche_res, out_plot_path = out_png)
}

.sdm_coerce_bool <- function(x, n) {
  if (is.null(x) || length(x) == 0L) return(rep(FALSE, n))
  if (length(x) != n) x <- rep_len(x, n)
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & is.finite(x) & (x != 0))
  sx <- tolower(trimws(as.character(x))); !is.na(sx) & (sx %in% c("true", "1", "t", "yes"))
}

.sdm_impute_with_medians <- function(raw_df, cols, meds) {
  if (exists("apply_fold_imputation", mode = "function")) return(apply_fold_imputation(raw_df, cols, meds))
  out <- as.data.frame(raw_df[, intersect(cols, colnames(raw_df)), drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  for (cn in cols) {
    mv <- if (!is.null(meds) && cn %in% names(meds) && is.finite(meds[[cn]])) as.numeric(meds[[cn]]) else 0.0
    v <- if (cn %in% colnames(out)) suppressWarnings(as.numeric(out[[cn]])) else rep(mv, nrow(raw_df))
    v[!is.finite(v)] <- mv; out[[cn]] <- v
  }
  out[, cols, drop = FALSE]
}

.sdm_compute_fold_spatial_weights <- function(coords_tr, p_tr, w_base) {
  n_tr <- length(p_tr); w_out <- as.numeric(w_base); w_out[!is.finite(w_out) | w_out <= 0.0] <- 1.0
  if (n_tr < 3L || !exists("compute_spatial_balance_weights", mode = "function")) return(w_out)
  sb_df <- data.frame(
    catalogNumber = paste0("SDM_TR_", seq_len(n_tr)),
    longitude = as.numeric(coords_tr[, 1L]), latitude = as.numeric(coords_tr[, 2L]),
    decimalLongitude = as.numeric(coords_tr[, 1L]), decimalLatitude = as.numeric(coords_tr[, 2L]),
    species_standardized = ifelse(p_tr == 1L, "Presence", "Background"), stringsAsFactors = FALSE
  )
  sb_w <- tryCatch(as.numeric(compute_spatial_balance_weights(sb_df)), error = function(e) rep(1.0, n_tr))
  if (length(sb_w) != n_tr) return(w_out)
  sb_w[!is.finite(sb_w) | sb_w <= 0.0] <- 1.0; sb_w <- pmin(4.0, pmax(0.10, sb_w))
  for (cls in c(1L, 0L)) { m_c <- which(p_tr == cls); if (length(m_c) >= 1L) sb_w[m_c] <- sb_w[m_c] / max(mean(sb_w[m_c], na.rm = TRUE), 1e-6) }
  w_out * sb_w
}

.sdm_prepare_fold_split <- function(X_tr_raw, p_tr, w_tr_base, coords_tr, X_te_raw = NULL, is_intra = FALSE, screen_fold = TRUE) {
  cand_v <- colnames(X_tr_raw)
  clim_v <- grep("^(bio[0-9]|summer_moisture_|maritime_continentality_)", cand_v, value = TRUE, ignore.case = TRUE)
  topo_v <- grep("^(topo_|elevation($|_)|elev_)", cand_v, value = TRUE, ignore.case = TRUE)
  edaph_v <- setdiff(cand_v, unique(c(clim_v, topo_v)))
  sel_v <- character(0); tr_med <- numeric(0); X_tr_imp <- NULL
  if (isTRUE(screen_fold) && exists("screen_predictors_within_cv_fold", mode = "function") && nrow(X_tr_raw) >= 8L && length(cand_v) >= 2L) {
    scr <- tryCatch(screen_predictors_within_cv_fold(X_tr_raw, cand_vars = cand_v, clim_vars = clim_v, edaph_vars = edaph_v, topo_vars = topo_v, y_tr = factor(ifelse(p_tr == 1L, "Presence", "Background")), cor.threshold = 0.70, vif.threshold = 5.0, intraspecific_mode = isTRUE(is_intra)), error = function(e) NULL)
    if (!is.null(scr) && length(scr$selected_vars) >= 1L && is.data.frame(scr$X_tr_imp) && nrow(scr$X_tr_imp) == nrow(X_tr_raw)) {
      sel_v <- as.character(scr$selected_vars); tr_med <- scr$train_medians; X_tr_imp <- scr$X_tr_imp
    }
  }
  if (is.null(X_tr_imp) || length(sel_v) == 0L) {
    if (isTRUE(screen_fold)) {
      ok_c <- vapply(cand_v, function(vn) { v <- suppressWarnings(as.numeric(X_tr_raw[[vn]])); f <- is.finite(v); sum(f) >= 2L && is.finite(stats::sd(v[f])) && stats::sd(v[f]) > 1e-8 }, logical(1))
      sel_v <- if (any(ok_c)) cand_v[ok_c] else cand_v
    } else sel_v <- cand_v
    tr_med <- stats::setNames(vapply(sel_v, function(vn) { v <- suppressWarnings(as.numeric(X_tr_raw[[vn]])); f <- is.finite(v); if (any(f)) { m <- stats::median(v[f], na.rm = TRUE); if (is.finite(m)) m else 0.0 } else 0.0 }, numeric(1)), sel_v)
    X_tr_imp <- .sdm_impute_with_medians(X_tr_raw, sel_v, tr_med)
  }
  X_te_imp <- if (!is.null(X_te_raw)) .sdm_impute_with_medians(X_te_raw, sel_v, tr_med) else NULL
  list(selected_vars = sel_v, train_medians = tr_med, X_tr_imp = X_tr_imp, X_te_imp = X_te_imp, w_tr = .sdm_compute_fold_spatial_weights(coords_tr, p_tr, w_tr_base))
}

fit_and_project_maxent_sdm_maps <- function(
  df, srf_res, env_dir, out_dir = "outputs", kept_vars = NULL, primary_clim_vars = NULL,
  primary_edaph_vars = NULL, focal_taxon = "Packera dubia", intraspecific_mode = FALSE, demo = FALSE
) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  if (is.null(out_dir) || !nzchar(out_dir)) out_dir <- "outputs"
  message("=== Fitting MaxEnt (cloglog) & Spatial RF Ensemble SDMs & Projecting Suitability Rasters ===")

  df_raw <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  if (!("longitude" %in% names(df_raw)) && "decimalLongitude" %in% names(df_raw)) df_raw$longitude <- df_raw$decimalLongitude
  if (!("latitude" %in% names(df_raw)) && "decimalLatitude" %in% names(df_raw)) df_raw$latitude <- df_raw$decimalLatitude

  if (is.null(kept_vars) || length(kept_vars) == 0L) {
    if (!is.null(srf_res) && !is.null(srf_res$collinearity_screening$kept_vars)) kept_vars <- srf_res$collinearity_screening$kept_vars
    else {
      imp_csv <- file.path(out_dir, "reports", "srf_feature_importance.csv")
      if (file.exists(imp_csv)) { imp_tbl <- tryCatch(utils::read.csv(imp_csv, stringsAsFactors = FALSE), error = function(e) NULL); if (!is.null(imp_tbl) && "Predictor" %in% names(imp_tbl)) kept_vars <- as.character(imp_tbl$Predictor) }
      if (is.null(kept_vars) || length(kept_vars) == 0L) {
        cand_eco <- discover_candidate_eco_predictors(df_raw, include_pheno = FALSE, min_finite = 2L)
        if (length(cand_eco) > 12L && exists("screen_collinear_predictors_by_block", mode = "function")) kept_vars <- tryCatch(screen_collinear_predictors_by_block(df_raw[, cand_eco, drop = FALSE], cor.threshold = 0.70, vif.threshold = 5.0)$kept_vars, error = function(e) NULL)
        if (is.null(kept_vars) || length(kept_vars) == 0L) { clim_c <- cand_eco[grepl("^bio[0-9]+", cand_eco)]; kept_vars <- c(utils::head(clim_c, 6L), utils::head(setdiff(cand_eco, clim_c), 6L)) }
      }
    }
  }
  kept_vars <- intersect(kept_vars[!grepl("^MEM", kept_vars)], names(df_raw))
  if (length(kept_vars) > 12L) {
    kv_clim <- kept_vars[grepl("^(bio[0-9]+|summer_moisture_|maritime_continentality_)", kept_vars)]; kv_edaph <- setdiff(kept_vars, kv_clim)
    n_c <- min(length(kv_clim), 6L); n_e <- min(length(kv_edaph), 12L - n_c)
    if ((n_c + n_e) < 12L && length(kv_clim) > n_c) n_c <- min(length(kv_clim), 12L - n_e)
    kept_vars <- c(utils::head(kv_clim, n_c), utils::head(kv_edaph, n_e))
  }
  if (length(kept_vars) == 0L) kept_vars <- intersect(c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_ph", "soil_ph", "fine_soil_sand", "soil_sand"), names(df_raw))
  if (length(kept_vars) == 0L) return(invisible(NULL))

  # Coerce predictor columns to numeric WITHOUT global pre-CV median imputation
  for (kv in kept_vars) df_raw[[kv]] <- suppressWarnings(as.numeric(df_raw[[kv]]))

  imp_df <- if (!is.null(srf_res) && is.data.frame(srf_res$importance)) srf_res$importance else NULL
  edaph_rx <- "^(fine_|ilr_texture_|ssurgo_|solus_|topo_|soil_|sgmc_|is_azonal|lithic_outcrop|min_bedrock|ph_q|max_cfvo|wetland_seepage)"
  clim_cand <- if (!is.null(primary_clim_vars) && length(intersect(primary_clim_vars, kept_vars)) > 0L) intersect(primary_clim_vars, kept_vars) else kept_vars[!grepl(edaph_rx, kept_vars)]
  edaph_cand <- if (!is.null(primary_edaph_vars) && length(intersect(primary_edaph_vars, kept_vars)) > 0L) intersect(primary_edaph_vars, kept_vars) else kept_vars[grepl(edaph_rx, kept_vars)]
  if (!is.null(imp_df) && all(c("Predictor", "MeanDecreaseAccuracy") %in% names(imp_df))) {
    imp_map <- stats::setNames(as.numeric(imp_df$MeanDecreaseAccuracy), as.character(imp_df$Predictor))
    if (length(clim_cand) > 1L) clim_cand <- clim_cand[order(ifelse(clim_cand %in% names(imp_map), imp_map[clim_cand], -Inf), decreasing = TRUE)]
    if (length(edaph_cand) > 1L) edaph_cand <- edaph_cand[order(ifelse(edaph_cand %in% names(imp_map), imp_map[edaph_cand], -Inf), decreasing = TRUE)]
  }
  top_clim_pred <- if (length(clim_cand) > 0L) clim_cand[1L] else kept_vars[1L]
  top_edaph_pred <- if (length(edaph_cand) > 0L) edaph_cand[1L] else "None_Retained"

  # Resolve authority weights & strict anchor/type hygiene (quarantined, open/hybrid & heterotypic-audit types excluded)
  n_df <- nrow(df_raw)
  raw_auth_w <- if ("Combined_Prior_Weight" %in% names(df_raw)) suppressWarnings(as.numeric(df_raw$Combined_Prior_Weight)) else if ("authority_weight" %in% names(df_raw)) suppressWarnings(as.numeric(df_raw$authority_weight)) else resolve_ecological_authority_weights(df_raw)
  raw_auth_w[!is.finite(raw_auth_w) | raw_auth_w < 0.0] <- 0.0
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df_raw) else {
    is_t <- .sdm_coerce_bool(if ("is_type_specimen" %in% names(df_raw)) df_raw$is_type_specimen else NULL, n_df); list(is_type_specimen = is_t, type_rank = ifelse(is_t, 2L, 0L))
  }
  is_open_hyb <- if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) .identify_open_or_hybrid_concept_mask(df_raw) else {
    .sdm_coerce_bool(if ("is_open_or_hybrid_concept" %in% names(df_raw)) df_raw$is_open_or_hybrid_concept else NULL, n_df) |
      (if ("concept_synonymy_status" %in% names(df_raw)) trimws(as.character(df_raw$concept_synonymy_status)) %in% c("Open_Nomenclature", "Hybrid_Formula") else rep(FALSE, n_df))
  }
  is_open_hyb[is.na(is_open_hyb)] <- FALSE
  is_het_audit <- rep(FALSE, n_df)
  for (.hc in intersect(c("rf_prediction_source", "anchor_tier_source"), names(df_raw))) is_het_audit <- is_het_audit | (!is.na(df_raw[[.hc]]) & trimws(as.character(df_raw[[.hc]])) == "Passive_Heterotypic_Type_Audit")
  if ("is_heterotypic_type_audit" %in% names(df_raw)) is_het_audit <- is_het_audit | .sdm_coerce_bool(df_raw$is_heterotypic_type_audit, n_df)

  lon_all <- suppressWarnings(as.numeric(df_raw$longitude)); lat_all <- suppressWarnings(as.numeric(df_raw$latitude))
  unc_col <- if ("coordinateUncertaintyInMeters" %in% names(df_raw)) suppressWarnings(as.numeric(df_raw$coordinateUncertaintyInMeters)) else if ("coordinateUncertainty" %in% names(df_raw)) suppressWarnings(as.numeric(df_raw$coordinateUncertainty)) else rep(NA_real_, n_df)
  coord_q_raw <- if (exists("detect_r_centroid_or_inst_quarantine", mode = "function")) detect_r_centroid_or_inst_quarantine(df_raw, is_type_vec = rep(FALSE, n_df))$is_quarantined else if ("is_centroid_or_inst_coord" %in% names(df_raw)) .sdm_coerce_bool(df_raw$is_centroid_or_inst_coord, n_df) else rep(FALSE, n_df)
  is_quar_coord_all <- coord_q_raw | (is.finite(unc_col) & (unc_col > 10000.0))

  is_anchor_elig <- if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) {
    tryCatch({
      ev <- is_fine_edaphic_type_or_anchor_eligible(df_raw, max_edaphic_unc_m = 10000.0, max_clim_unc_m = 50000.0)
      if (is.list(ev) && !is.null(ev$edaphic_ok) && length(ev$edaphic_ok) == n_df) .sdm_coerce_bool(ev$edaphic_ok, n_df)
      else if (length(ev) == n_df) .sdm_coerce_bool(ev, n_df) else rep(TRUE, n_df)
    }, error = function(e) rep(TRUE, n_df))
  } else rep(TRUE, n_df)

  cl_sdm <- .sdm_coerce_bool(if ("is_confident_label_noise" %in% names(df_raw)) df_raw$is_confident_label_noise else NULL, n_df)
  mcd_sdm <- if ("env_mcd_outlier_flag" %in% names(df_raw)) .sdm_coerce_bool(df_raw$env_mcd_outlier_flag, n_df) else if (!is.null(srf_res) && is.list(srf_res) && length(srf_res$env_mcd_outlier_flag) == n_df) .sdm_coerce_bool(srf_res$env_mcd_outlier_flag, n_df) else rep(FALSE, n_df)

  is_hard_excl <- is_open_hyb | is_het_audit | is_quar_coord_all | !is_anchor_elig
  is_elig_type <- type_info$is_type_specimen & !is_hard_excl
  raw_auth_w[is_hard_excl] <- 0.0; raw_auth_w[(cl_sdm | mcd_sdm) & !is_elig_type] <- 0.0
  pos_m <- raw_auth_w > 0.0; w_max_curr <- if (any(pos_m)) max(raw_auth_w[pos_m], 1e-6) else 1.0
  auth_w_01 <- ifelse(pos_m, pmin(1.0, pmax(0.20, raw_auth_w / w_max_curr)), 0.0)
  auth_w_01[is_elig_type] <- 1.00
  df_raw$sdm_prior_weight_01 <- auth_w_01

  elig_coord_only <- is.finite(lon_all) & is.finite(lat_all) & !is_quar_coord_all & (is.na(unc_col) | unc_col <= 10000.0)
  if (sum(elig_coord_only) < 4L) elig_coord_only <- is.finite(lon_all) & is.finite(lat_all) & !is_quar_coord_all
  elig_mask <- elig_coord_only & (df_raw$sdm_prior_weight_01 > 0.0) & !is_hard_excl
  if (sum(elig_mask) < 4L) {
    elig_clean <- elig_coord_only & !is_hard_excl & (!(cl_sdm | mcd_sdm) | is_elig_type)
    elig_fb <- elig_coord_only & !is_hard_excl & (!cl_sdm | is_elig_type)
    elig_mask <- if (sum(elig_clean) >= 4L || sum(elig_clean) >= sum(elig_fb)) elig_clean else elig_fb
    df_raw$sdm_prior_weight_01[elig_mask & df_raw$sdm_prior_weight_01 <= 0.0] <- 0.50
  }

  env_stack <- .build_regional_sdm_raster_stack(df_raw[elig_coord_only, , drop = FALSE], env_dir = env_dir, kept_vars = kept_vars, demo = demo)
  if (!is.null(env_stack) && length(intersect(kept_vars, names(env_stack))) >= 2L) {
    kept_vars <- intersect(kept_vars, names(env_stack))
    if (!(top_clim_pred %in% kept_vars)) top_clim_pred <- kept_vars[1L]
    ed_st <- intersect(edaph_cand, kept_vars)
    top_edaph_pred <- if (length(ed_st) > 0L) ed_st[1L] else if (length(edaph_cand) > 0L) edaph_cand[1L] else "None_Retained"
  }

  set.seed(42); bg_pool_df <- NULL
  if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE)) {
    s_idx <- sample.int(terra::ncell(env_stack), size = min(1200L, terra::ncell(env_stack)), replace = FALSE)
    s_xy <- terra::xyFromCell(env_stack, s_idx); s_vals <- as.data.frame(terra::values(env_stack)[s_idx, , drop = FALSE]); colnames(s_vals) <- names(env_stack)
    bg_pool_df <- data.frame(longitude = as.numeric(s_xy[, 1L]), latitude = as.numeric(s_xy[, 2L]), s_vals[, kept_vars, drop = FALSE], sdm_prior_weight_01 = 1.0, check.names = FALSE, stringsAsFactors = FALSE)
    bg_pool_df <- bg_pool_df[stats::complete.cases(bg_pool_df[, c("longitude", "latitude", kept_vars), drop = FALSE]), , drop = FALSE]
  }

  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df_raw$species_standardized))) < 2L)
  grp_vec <- if (is_intra && "analysis_group" %in% names(df_raw)) as.character(df_raw$analysis_group) else as.character(df_raw$species_standardized)
  grp_vec[is.na(grp_vec) | !nzchar(grp_vec)] <- "Unknown"
  tbl_elig <- table(grp_vec[elig_mask & grp_vec != "Unknown"])
  target_groups <- names(tbl_elig[tbl_elig >= 5L])
  if (length(target_groups) == 0L) target_groups <- names(tbl_elig[tbl_elig >= 2L])
  if (length(target_groups) == 0L && any(grp_vec != "Unknown")) target_groups <- unique(grp_vec[grp_vec != "Unknown"])
  if (!is.null(focal_taxon) && nzchar(focal_taxon) && !(focal_taxon %in% target_groups) && sum(df_raw$species_standardized == focal_taxon & elig_mask, na.rm = TRUE) >= 3L) target_groups <- c(focal_taxon, target_groups)
  if (!is.null(focal_taxon) && focal_taxon %in% target_groups) target_groups <- c(focal_taxon, setdiff(target_groups, focal_taxon))

  block_deg_use <- if (!is.null(srf_res) && is.finite(srf_res$spatial_autocorr_range_km)) max(0.25, min(2.5, as.numeric(srf_res$spatial_autocorr_range_km) / 111.0)) else 0.50
  raster_dir <- file.path(out_dir, "rasters")
  write_prod_rasters <- !(isTRUE(demo) && is_canonical_prod_path(raster_dir))
  if (write_prod_rasters) dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
  gdal_deflate <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

  tg_bundles <- list(); fit_tasks <- list()
  for (tg in target_groups) {
    use_sp <- (tg %in% df_raw$species_standardized) && !(tg %in% grp_vec)
    pres_mask <- elig_mask & if (use_sp) (df_raw$species_standardized == tg) else (grp_vec == tg)
    pres_df <- df_raw[pres_mask, c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
    n_pres <- nrow(pres_df); if (n_pres < 2L) next
    tg_bg_df <- df_raw[elig_coord_only & !(if (use_sp) (df_raw$species_standardized == tg) else (grp_vec == tg)), c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
    if (nrow(tg_bg_df) > 0L) tg_bg_df$sdm_prior_weight_01 <- 1.0
    n_need <- max(0L, 1000L - nrow(tg_bg_df))
    bg_df <- if (n_need > 0L && !is.null(bg_pool_df) && nrow(bg_pool_df) > 0L) rbind(tg_bg_df, bg_pool_df[seq_len(min(n_need, nrow(bg_pool_df))), c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]) else if (nrow(tg_bg_df) > 0L) tg_bg_df else if (!is.null(bg_pool_df) && nrow(bg_pool_df) > 0L) bg_pool_df else NULL
    if (is.null(bg_df) || nrow(bg_df) == 0L) next

    comb_df <- rbind(pres_df, bg_df); n_bg <- nrow(bg_df)
    p_vec <- c(rep(1L, n_pres), rep(0L, n_bg)); w_vec <- as.numeric(comb_df$sdm_prior_weight_01)
    X_comb_raw <- comb_df[, kept_vars, drop = FALSE]; coords_comb <- as.matrix(comb_df[, c("longitude", "latitude")])
    k_folds <- min(5L, max(2L, n_pres))
    block_ids <- assign_spatial_grid_blocks(coords_comb, block_deg = block_deg_use, k = k_folds)
    if (length(unique(block_ids[p_vec == 1L])) < 2L && n_pres >= 2L) block_ids[p_vec == 1L] <- ((seq_len(n_pres) - 1L) %% min(k_folds, n_pres)) + 1L
    u_folds <- sort(unique(block_ids))
    tg_bundles[[tg]] <- list(tg = tg, pres_mask = pres_mask, n_pres = n_pres, n_bg = n_bg, p_vec = p_vec, w_vec = w_vec, X_comb_raw = X_comb_raw, coords_comb = coords_comb)
    for (fid in u_folds) {
      tr_idx <- which(block_ids != fid); te_idx <- which(block_ids == fid)
      if (sum(p_vec[tr_idx] == 1L) >= 2L && sum(p_vec[tr_idx] == 0L) >= 2L && length(te_idx) >= 1L) fit_tasks[[length(fit_tasks) + 1L]] <- list(kind = "cv", tg = tg, tr_idx = tr_idx, te_idx = te_idx)
    }
    fit_tasks[[length(fit_tasks) + 1L]] <- list(kind = "full", tg = tg)
  }

  max_w <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L) else 1L
  tot_sdm_n <- if (length(tg_bundles) > 0L) sum(vapply(tg_bundles, function(b) as.integer(b$n_pres + b$n_bg), integer(1))) else 0L
  eff_w <- if (length(fit_tasks) >= 2L && tot_sdm_n >= 200L && exists("packera_mclapply", mode = "function")) min(max_w, length(fit_tasks)) else 1L
  job_rf_th <- max(1L, min(as.integer(floor(max_w / max(1L, eff_w))), as.integer(ceiling(max(1L, tot_sdm_n) / 64L))))

  fit_res_list <- list()
  if (length(fit_tasks) > 0L) {
    run_fit_task <- function(t_i) {
      tk <- fit_tasks[[t_i]]; bd <- tg_bundles[[tk$tg]]
      if (identical(tk$kind, "cv")) {
        tr <- tk$tr_idx; te <- tk$te_idx
        fd <- .sdm_prepare_fold_split(bd$X_comb_raw[tr, , drop = FALSE], bd$p_vec[tr], bd$w_vec[tr], bd$coords_comb[tr, , drop = FALSE], X_te_raw = bd$X_comb_raw[te, , drop = FALSE], is_intra = is_intra, screen_fold = TRUE)
        mx_cv <- .fit_single_maxent_model(bd$p_vec[tr], fd$X_tr_imp, fd$w_tr)
        rf_cv <- .fit_single_srf_binary_model(bd$p_vec[tr], fd$X_tr_imp, fd$w_tr, n_trees = 200L, num_threads = job_rf_th)
        list(oof_mx = mx_cv$predict_cloglog(fd$X_te_imp), oof_rf = rf_cv$predict_prob(fd$X_te_imp))
      } else {
        fd_f <- .sdm_prepare_fold_split(bd$X_comb_raw, bd$p_vec, bd$w_vec, bd$coords_comb, X_te_raw = NULL, is_intra = is_intra, screen_fold = FALSE)
        list(mx_full = .fit_single_maxent_model(bd$p_vec, fd_f$X_tr_imp, fd_f$w_tr), rf_full = .fit_single_srf_binary_model(bd$p_vec, fd_f$X_tr_imp, fd_f$w_tr, n_trees = 300L, num_threads = job_rf_th), X_full_imp = fd_f$X_tr_imp)
      }
    }
    fit_res_list <- if (eff_w > 1L) packera_mclapply(seq_along(fit_tasks), run_fit_task, mc.cores = eff_w, mc.preschedule = FALSE) else lapply(seq_along(fit_tasks), run_fit_task)
  }

  .eval_oof_roc_cbi <- function(p_v, pred_v, w_v, ok_m) {
    idx <- if (is.logical(ok_m) && length(ok_m) == length(p_v) && sum(p_v[ok_m] == 1L) >= 2L && sum(p_v[ok_m] == 0L) >= 2L) which(ok_m) else seq_along(p_v)
    pe <- p_v[idx]; pr <- pred_v[idx]
    list(roc = .compute_roc_auc_curve(pe, pr, w_v[idx]), cbi = compute_continuous_boyce_index(pr, pr[pe == 1L], n_bins = 20L))
  }
  .pick_tg_top_pred <- function(cands, X_all, mx_m, rf_m, fb) {
    c_use <- intersect(cands, colnames(X_all)); if (length(c_use) == 0L) return(fb); if (length(c_use) == 1L) return(c_use[1L])
    base_p <- 0.5 * mx_m$predict_cloglog(X_all) + 0.5 * rf_m$predict_prob(X_all)
    sc <- vapply(c_use, function(vn) {
      Xm <- X_all; cv <- suppressWarnings(as.numeric(Xm[[vn]])); f <- is.finite(cv)
      if (sum(f) < 2L || stats::sd(cv[f]) < 1e-10) return(0.0)
      Xm[[vn]] <- rep(stats::median(cv[f], na.rm = TRUE), nrow(Xm))
      v <- mean(abs(base_p - (0.5 * mx_m$predict_cloglog(Xm) + 0.5 * rf_m$predict_prob(Xm))), na.rm = TRUE); if (is.finite(v)) as.numeric(v) else 0.0
    }, numeric(1))
    if (all(sc <= 1e-10)) c_use[1L] else c_use[which.max(sc)]
  }

  perf_rows <- list(); focal_artifacts <- NULL; all_sdm_artifacts <- list(); tg_models <- list()
  for (tg in names(tg_bundles)) {
    bd <- tg_bundles[[tg]]; p_vec <- bd$p_vec; w_vec <- bd$w_vec
    oof_mx <- rep(NA_real_, length(p_vec)); oof_rf <- rep(NA_real_, length(p_vec)); mx_full <- NULL; rf_full <- NULL; X_comb <- NULL
    for (t_i in seq_along(fit_tasks)) if (identical(fit_tasks[[t_i]]$tg, tg)) {
      fr <- fit_res_list[[t_i]]
      if (identical(fit_tasks[[t_i]]$kind, "cv")) { oof_mx[fit_tasks[[t_i]]$te_idx] <- fr$oof_mx; oof_rf[fit_tasks[[t_i]]$te_idx] <- fr$oof_rf }
      else { mx_full <- fr$mx_full; rf_full <- fr$rf_full; X_comb <- fr$X_full_imp }
    }
    ok_oof_mx <- is.finite(oof_mx); ok_oof_rf <- is.finite(oof_rf)
    if (any(!ok_oof_mx)) oof_mx[!ok_oof_mx] <- mx_full$predict_cloglog(X_comb[!ok_oof_mx, , drop = FALSE])
    if (any(!ok_oof_rf)) oof_rf[!ok_oof_rf] <- rf_full$predict_prob(X_comb[!ok_oof_rf, , drop = FALSE])

    eval_mx <- .eval_oof_roc_cbi(p_vec, oof_mx, w_vec, ok_oof_mx)
    eval_rf <- .eval_oof_roc_cbi(p_vec, oof_rf, w_vec, ok_oof_rf)
    auc_mx <- if (is.finite(eval_mx$roc$auc)) eval_mx$roc$auc else 0.5
    auc_rf <- if (is.finite(eval_rf$roc$auc)) eval_rf$roc$auc else 0.5
    w_mx <- if ((auc_mx + auc_rf) > 1e-6) auc_mx / (auc_mx + auc_rf) else 0.5; w_rf <- 1.0 - w_mx
    oof_ens <- pmin(1.0, pmax(0.0, w_mx * oof_mx + w_rf * oof_rf))
    ok_oof_ens <- ok_oof_mx & ok_oof_rf
    eval_ens <- .eval_oof_roc_cbi(p_vec, oof_ens, w_vec, ok_oof_ens)
    auc_ens <- if (is.finite(eval_ens$roc$auc)) eval_ens$roc$auc else round(0.5 * (auc_mx + auc_rf), 4)

    tg_clim <- .pick_tg_top_pred(clim_cand, X_comb, mx_full, rf_full, top_clim_pred)
    tg_edaph <- .pick_tg_top_pred(edaph_cand, X_comb, mx_full, rf_full, top_edaph_pred)
    perf_rows[[length(perf_rows) + 1L]] <- data.frame(
      Taxon_or_Morphotype = tg, N_Presences = as.integer(bd$n_pres), N_Background = as.integer(bd$n_bg), MaxEnt_Engine = as.character(mx_full$engine),
      Selected_Feature_Class = if (!is.null(mx_full$selected_feature_class)) as.character(mx_full$selected_feature_class) else "L",
      Selected_Reg_Multiplier = if (!is.null(mx_full$selected_reg_multiplier) && is.finite(mx_full$selected_reg_multiplier)) round(as.numeric(mx_full$selected_reg_multiplier), 3) else NA_real_,
      MaxEnt_AICc = if (!is.null(mx_full$maxent_aicc) && is.finite(mx_full$maxent_aicc)) round(as.numeric(mx_full$maxent_aicc), 3) else NA_real_,
      MaxEnt_NonZero_Betas = if (!is.null(mx_full$n_nonzero_betas)) as.integer(mx_full$n_nonzero_betas) else NA_integer_,
      AUC_Block_MaxEnt = round(auc_mx, 4), CBI_MaxEnt = round(as.numeric(eval_mx$cbi), 4),
      AUC_Block_SpatialRF = round(auc_rf, 4), CBI_SpatialRF = round(as.numeric(eval_rf$cbi), 4),
      AUC_Block_Ensemble = round(auc_ens, 4), CBI_Ensemble = round(as.numeric(eval_ens$cbi), 4),
      Top_Climatic_Predictor = as.character(tg_clim), Top_Edaphic_Predictor = as.character(tg_edaph), stringsAsFactors = FALSE
    )
    tg_models[[tg]] <- list(taxon = tg, pres_mask = bd$pres_mask, mx_full = mx_full, rf_full = rf_full, w_ens_mx = w_mx, w_ens_rf = w_rf, roc_mx = eval_mx$roc, roc_rf = eval_rf$roc, roc_ens = eval_ens$roc, cbi_mx = eval_mx$cbi, cbi_rf = eval_rf$cbi, cbi_ens = eval_ens$cbi, top_clim_pred = tg_clim, top_edaph_pred = tg_edaph, X_comb = X_comb, p_vec = p_vec)
  }

  rast_preds <- list()
  if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE) && length(tg_models) > 0L) {
    n_cells <- terra::ncell(env_stack); ev_df <- as.data.frame(terra::values(env_stack)); colnames(ev_df) <- names(env_stack)
    ok_cells <- which(stats::complete.cases(ev_df[, kept_vars, drop = FALSE]))
    if (length(ok_cells) > 0L) {
      eok <- ev_df[ok_cells, kept_vars, drop = FALSE]; n_ok <- nrow(eok)
      nch <- max(1L, min(as.integer(ceiling(n_ok / 250L)), as.integer(ceiling(max_w / length(tg_models)))))
      rchunks <- if (nch > 1L) split(seq_len(n_ok), cut(seq_len(n_ok), breaks = nch, labels = FALSE)) else list(seq_len(n_ok))
      ptasks <- unlist(lapply(names(tg_models), function(tg) lapply(rchunks, function(rc) list(tg = tg, idx = rc))), recursive = FALSE)
      pfn <- function(pt) { nd <- eok[pt$idx, , drop = FALSE]; tm <- tg_models[[pt$tg]]; list(mx = as.numeric(tm$mx_full$predict_cloglog(nd)), rf = as.numeric(tm$rf_full$predict_prob(nd))) }
      pres <- if (length(ptasks) >= 2L && max_w > 1L && exists("packera_mclapply", mode = "function")) packera_mclapply(ptasks, pfn, mc.cores = min(max_w, length(ptasks)), mc.preschedule = FALSE) else lapply(ptasks, pfn)
      for (tg in names(tg_models)) {
        mx_ok <- numeric(n_ok); rf_ok <- numeric(n_ok)
        for (pi in seq_along(ptasks)) if (identical(ptasks[[pi]]$tg, tg)) { mx_ok[ptasks[[pi]]$idx] <- pres[[pi]]$mx; rf_ok[ptasks[[pi]]$idx] <- pres[[pi]]$rf }
        rast_preds[[tg]] <- list(n_cells = n_cells, ok_cells = ok_cells, mx_ok = mx_ok, rf_ok = rf_ok)
      }
    }
    r_tpl <- terra::rast(env_stack[[1L]])
    .write_tif <- function(r_obj, p_out) {
      tt <- paste0(p_out, ".tmp.tif"); dir.create(dirname(p_out), recursive = TRUE, showWarnings = FALSE)
      tryCatch(terra::writeRaster(r_obj, filename = tt, overwrite = TRUE, gdal = gdal_deflate), error = function(e) terra::writeRaster(r_obj, filename = tt, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6")))
      if (file.exists(tt)) file.rename(tt, p_out)
    }
  }

  for (tg in names(tg_models)) {
    tm <- tg_models[[tg]]; r_mx <- NULL; r_rf <- NULL; r_ens <- NULL
    if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE) && !is.null(rast_preds[[tg]])) {
      rp <- rast_preds[[tg]]; vmx <- rep(NA_real_, rp$n_cells); vmx[rp$ok_cells] <- rp$mx_ok; vrf <- rep(NA_real_, rp$n_cells); vrf[rp$ok_cells] <- rp$rf_ok
      r_mx <- terra::rast(r_tpl); terra::values(r_mx) <- vmx; names(r_mx) <- "maxent_cloglog"
      r_rf <- terra::rast(r_tpl); terra::values(r_rf) <- vrf; names(r_rf) <- "spatial_rf_suitability"
      r_ens <- tm$w_ens_mx * r_mx + tm$w_ens_rf * r_rf; names(r_ens) <- "ensemble_suitability"
      snm <- gsub("^_+|_+$", "", tolower(gsub("[^A-Za-z0-9]+", "_", tg)))
      if (write_prod_rasters) {
        tryCatch(.write_tif(r_mx, file.path(raster_dir, sprintf("sdm_%s_maxent_cloglog.tif", snm))), error = function(e) NULL)
        tryCatch(.write_tif(r_ens, file.path(raster_dir, sprintf("sdm_%s_ensemble.tif", snm))), error = function(e) NULL)
      }
    }
    tg_art <- list(taxon = tg, pres_mask = tm$pres_mask, mx_full = tm$mx_full, rf_full = tm$rf_full, rast_mx = r_mx, rast_rf = r_rf, rast_ens = r_ens, roc_mx = tm$roc_mx, roc_rf = tm$roc_rf, roc_ens = tm$roc_ens, cbi_mx = tm$cbi_mx, cbi_rf = tm$cbi_rf, cbi_ens = tm$cbi_ens, top_clim_pred = tm$top_clim_pred, top_edaph_pred = tm$top_edaph_pred, X_comb = tm$X_comb, p_vec = tm$p_vec)
    all_sdm_artifacts[[tg]] <- tg_art; if (is.null(focal_artifacts)) focal_artifacts <- tg_art
  }
  if (!is.null(focal_taxon) && nzchar(focal_taxon) && (focal_taxon %in% names(all_sdm_artifacts))) focal_artifacts <- all_sdm_artifacts[[focal_taxon]]

  perf_df <- if (length(perf_rows) > 0L) do.call(rbind, perf_rows) else data.frame(
    Taxon_or_Morphotype = character(0), N_Presences = integer(0), N_Background = integer(0), MaxEnt_Engine = character(0),
    Selected_Feature_Class = character(0), Selected_Reg_Multiplier = numeric(0), MaxEnt_AICc = numeric(0), MaxEnt_NonZero_Betas = integer(0),
    AUC_Block_MaxEnt = numeric(0), CBI_MaxEnt = numeric(0), AUC_Block_SpatialRF = numeric(0), CBI_SpatialRF = numeric(0),
    AUC_Block_Ensemble = numeric(0), CBI_Ensemble = numeric(0), Top_Climatic_Predictor = character(0), Top_Edaphic_Predictor = character(0), stringsAsFactors = FALSE
  )
  perf_csv_out <- file.path(out_dir, "reports", "sdm_model_performance.csv")
  if (write_csv_atomic(perf_df, perf_csv_out, demo = demo)) message(sprintf("Saved SDM performance metrics (%d modeled taxa/morphotypes) to %s.", nrow(perf_df), perf_csv_out))

  fig_dir <- file.path(out_dir, "figures")
  if (!(isTRUE(demo) && is_canonical_prod_path(fig_dir))) {
    .ensure_step05_fig_modules(); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(focal_artifacts) && !is.null(focal_artifacts$rast_mx) && !is.null(focal_artifacts$rast_ens)) {
      rv <- intersect(unique(c(focal_artifacts$top_clim_pred, focal_artifacts$top_edaph_pred, clim_cand, edaph_cand, kept_vars)), colnames(focal_artifacts$X_comb))
      rv <- rv[seq_len(min(3L, length(rv)))]
      rlbl <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(rv) else gsub("_", " ", rv)
      pm <- focal_artifacts$X_comb[focal_artifacts$p_vec == 1L, , drop = FALSE]; if (nrow(pm) == 0L) pm <- focal_artifacts$X_comb
      mrow <- as.data.frame(as.list(apply(pm, 2L, stats::median, na.rm = TRUE)), check.names = FALSE); ug <- seq(0, 1, length.out = 100L)
      rrows <- lapply(seq_along(rv), function(vi) {
        vn <- rv[vi]; vr <- as.numeric(focal_artifacts$X_comb[[vn]])
        qlo <- as.numeric(stats::quantile(vr, 0.05, na.rm = TRUE)); qhi <- as.numeric(stats::quantile(vr, 0.95, na.rm = TRUE))
        if (!is.finite(qlo) || !is.finite(qhi) || abs(qhi - qlo) < 1e-8) { qlo <- min(vr, na.rm = TRUE); qhi <- max(vr, na.rm = TRUE) }
        edf <- mrow[rep(1L, 100L), , drop = FALSE]; xs <- qlo + ug * (qhi - qlo); edf[[vn]] <- xs
        data.frame(predictor = rep(rlbl[vi], 100L), x_val = xs, suitability = as.numeric(focal_artifacts$mx_full$predict_cloglog(edf)), stringsAsFactors = FALSE)
      })
      roc_p <- Filter(Negate(is.null), list(
        if (!is.null(focal_artifacts$roc_mx)) data.frame(fpr = focal_artifacts$roc_mx$fpr, tpr = focal_artifacts$roc_mx$tpr, model = "MaxEnt", stringsAsFactors = FALSE),
        if (!is.null(focal_artifacts$roc_rf)) data.frame(fpr = focal_artifacts$roc_rf$fpr, tpr = focal_artifacts$roc_rf$tpr, model = "Spatial RF", stringsAsFactors = FALSE),
        if (!is.null(focal_artifacts$roc_ens)) data.frame(fpr = focal_artifacts$roc_ens$fpr, tpr = focal_artifacts$roc_ens$tpr, model = "Ensemble", stringsAsFactors = FALSE)
      ))
      vm <- Filter(function(a) !is.null(a$rast_ens), all_sdm_artifacts)
      mlist <- if (length(vm) > 1L) stats::setNames(lapply(vm, `[[`, "rast_ens"), vapply(vm, function(a) as.character(a$taxon), character(1))) else NULL
      render_step05c_sdm_suitability_png(
        df = df_raw, focal_rasters = list(maxent = focal_artifacts$rast_mx, ensemble = focal_artifacts$rast_ens),
        multi_rasters = mlist, perf_df = perf_df, resp_df = if (length(rrows) > 0L) do.call(rbind, rrows) else NULL,
        roc_df = if (length(roc_p) > 0L) do.call(rbind, roc_p) else NULL, boyce_df = attr(focal_artifacts$cbi_ens, "boyce_curve"),
        figures_dir = fig_dir, focal_taxon = focal_artifacts$taxon
      )
      message(sprintf("Exported 300-DPI Ensemble SDM & Distribution plate -> %s", file.path(fig_dir, "fig05c_sdm_suitability_and_distribution_map.png")))
    } else render_step05c_sdm_suitability_png(df = df_raw, figures_dir = fig_dir, focal_taxon = focal_taxon)
  }
  invisible(list(performance = perf_df, focal = focal_artifacts))
}
