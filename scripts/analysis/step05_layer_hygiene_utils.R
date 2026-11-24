# ==============================================================================
# Module: scripts/analysis/step05_layer_hygiene_utils.R
# Purpose: Step 05 layer hygiene guards (DS-801/NWI/NWPL exclusion, 300 m
#          bedrock cap, <=250 m micro-scale gate, POLARIS retirement, fine-grain
#          coordinate/type/open-hybrid eligibility, and MEM decoupling).
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

STEP05_MICRO_SCALE_PREDICTORS <- c("topo_tpi_30m", "topo_hli_30m", "topo_twi_30m", "wetland_seepage_index")
STEP05_NON_ENM_AUDIT_ONLY_COLS <- c(
  "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
  "ultramafic_calcareous_gradient", "is_nwi_wetland", "nwi_wetland_system",
  "nwi_wetland_dist_m", "nwpl_indicator_code", "nwpl_wetland_score"
)
STEP05_DISABLED_POLARIS_COLS <- c(
  "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m",
  "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m"
)

.lh_to_bool_vec <- function(x, n) {
  if (n <= 0L) return(logical(0))
  if (is.null(x) || length(x) == 0L) return(rep(FALSE, n))
  if (length(x) != n) x <- rep_len(x, n)
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & is.finite(x) & (x != 0))
  sx <- tolower(trimws(as.character(x)))
  !is.na(sx) & (sx %in% c("true", "1", "t", "yes"))
}

clamp_bedrock_buffer_radius_m <- function(buffer_radius_m, min_r = 100.0, max_bedrock_r = 300.0) {
  r <- suppressWarnings(as.numeric(buffer_radius_m))
  if (length(r) == 0L) return(numeric(0))
  bad <- !is.finite(r) | r <= 0.0
  if (any(bad)) r[bad] <- 250.0
  pmin(as.numeric(max_bedrock_r), pmax(as.numeric(min_r), r))
}

is_micro_scale_coord_eligible <- function(df, max_uncertainty_m = 250.0) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(logical(0))
  nms <- names(df)
  has_unc_m <- "coordinateUncertaintyInMeters" %in% nms
  has_unc_alt <- "coordinateUncertainty" %in% nms
  if (!has_unc_m && !has_unc_alt) return(rep(TRUE, nrow(df)))
  u <- if (has_unc_m) suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters)) else rep(NA_real_, nrow(df))
  if (has_unc_alt) {
    u_alt <- suppressWarnings(as.numeric(df$coordinateUncertainty))
    miss_u <- is.na(u) | !is.finite(u)
    if (any(miss_u)) u[miss_u] <- u_alt[miss_u]
  }
  is.finite(u) & (u > 0.0) & (u <= as.numeric(max_uncertainty_m))
}

is_fine_edaphic_type_or_anchor_eligible <- function(df, max_edaphic_unc_m = 10000.0, max_clim_unc_m = 50000.0) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(list(edaphic_ok = logical(0), clim_ok = logical(0), is_open_hyb = logical(0), is_quarantined_type = logical(0)))
  }
  n <- nrow(df); nms <- names(df)
  is_type <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    t_res <- tryCatch(.identify_type_specimen_vec(df), error = function(e) NULL)
    if (!is.null(t_res$is_type_specimen) && length(t_res$is_type_specimen) == n) {
      .lh_to_bool_vec(t_res$is_type_specimen, n)
    } else if ("is_type_specimen" %in% nms) {
      .lh_to_bool_vec(df$is_type_specimen, n)
    } else {
      rep(FALSE, n)
    }
  } else if ("is_type_specimen" %in% nms) {
    .lh_to_bool_vec(df$is_type_specimen, n)
  } else {
    rep(FALSE, n)
  }

  u <- if ("coordinateUncertaintyInMeters" %in% nms) suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters)) else rep(NA_real_, n)
  if ("coordinateUncertainty" %in% nms) {
    u_alt <- suppressWarnings(as.numeric(df$coordinateUncertainty))
    miss_u <- is.na(u) | !is.finite(u); if (any(miss_u)) u[miss_u] <- u_alt[miss_u]
  }
  is_q <- if ("is_centroid_or_inst_coord" %in% nms) {
    .lh_to_bool_vec(df$is_centroid_or_inst_coord, n)
  } else if (exists("detect_r_centroid_or_inst_quarantine", mode = "function")) {
    q_res <- tryCatch(detect_r_centroid_or_inst_quarantine(df, is_type_vec = is_type, exempt_types = FALSE), error = function(e) NULL)
    if (is.list(q_res) && length(q_res$is_quarantined) == n) .lh_to_bool_vec(q_res$is_quarantined, n) else rep(FALSE, n)
  } else rep(FALSE, n)

  edaphic_ok <- !is_q & (is.na(u) | !is.finite(u) | (u <= as.numeric(max_edaphic_unc_m)))
  clim_ok <- !is_q & (is.na(u) | !is.finite(u) | (u <= as.numeric(max_clim_unc_m)))

  is_open_hyb <- if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) {
    oh <- tryCatch(.identify_open_or_hybrid_concept_mask(df), error = function(e) rep(FALSE, n))
    if (length(oh) == n) .lh_to_bool_vec(oh, n) else rep(FALSE, n)
  } else {
    oh <- rep(FALSE, n)
    if ("is_open_or_hybrid_concept" %in% nms) oh <- oh | .lh_to_bool_vec(df$is_open_or_hybrid_concept, n)
    if ("concept_synonymy_status" %in% nms) {
      cs <- trimws(as.character(df$concept_synonymy_status))
      oh <- oh | (!is.na(cs) & (cs %in% c("Open_Nomenclature", "Hybrid_Formula")))
    }
    oh
  }
  list(
    edaphic_ok = as.logical(edaphic_ok),
    clim_ok = as.logical(clim_ok),
    is_open_hyb = as.logical(is_open_hyb),
    is_quarantined_type = as.logical(is_type & !edaphic_ok)
  )
}

mask_micro_scale_predictors_by_uncertainty <- function(df, cols = STEP05_MICRO_SCALE_PREDICTORS, max_uncertainty_m = 250.0) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  elig <- is_micro_scale_coord_eligible(df, max_uncertainty_m = max_uncertainty_m)
  target_cols <- intersect(as.character(cols), names(df))
  if (length(target_cols) == 0L || all(elig)) return(df)
  for (col_nm in target_cols) {
    v <- suppressWarnings(as.numeric(df[[col_nm]])); v[!elig] <- NA_real_; df[[col_nm]] <- v
  }
  df
}

sanitize_disabled_polaris_columns <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  for (col_nm in intersect(STEP05_DISABLED_POLARIS_COLS, names(df))) df[[col_nm]] <- rep(NA_real_, n)
  df
}

cap_diagnostic_mem_matrix <- function(MEMs, n_anchors, min_env = 4L) {
  if (is.null(MEMs) || is.null(ncol(MEMs)) || ncol(MEMs) == 0L) return(NULL)
  n_a <- max(1L, as.integer(n_anchors))
  total_b <- max(1L, as.integer(floor(n_a / 5L)))
  mem_cap <- if (n_a < 30L) 2L else min(5L, max(2L, total_b - as.integer(min_env)))
  MEMs[, seq_len(min(ncol(MEMs), as.integer(mem_cap))), drop = FALSE]
}

.extract_rf_perm_importance <- function(rf_obj) {
  if (is.null(rf_obj)) return(numeric(0))
  if (inherits(rf_obj, "ranger") && requireNamespace("ranger", quietly = TRUE)) {
    v <- tryCatch(ranger::importance(rf_obj), error = function(e) NULL)
    if (!is.null(v) && length(v) > 0L) return(stats::setNames(as.numeric(v), names(v)))
  }
  if (requireNamespace("randomForest", quietly = TRUE)) {
    m <- tryCatch(randomForest::importance(rf_obj, type = 1), error = function(e) NULL)
    if (!is.null(m) && is.matrix(m) && ncol(m) >= 1L) return(stats::setNames(as.numeric(m[, 1L]), rownames(m)))
  }
  if (!is.null(rf_obj$variable.importance)) return(stats::setNames(as.numeric(rf_obj$variable.importance), names(rf_obj$variable.importance)))
  numeric(0)
}

assemble_mem_decoupled_importance_df <- function(rf_nonspatial, rf_spatial, env_vars) {
  imp_ns <- .extract_rf_perm_importance(rf_nonspatial)
  env_clean <- setdiff(as.character(env_vars), c(STEP05_NON_ENM_AUDIT_ONLY_COLS, STEP05_DISABLED_POLARIS_COLS))
  env_clean <- env_clean[!grepl("^(MEM[0-9]+|polaris_)", env_clean)]
  env_keep <- if (length(env_clean) > 0L) intersect(env_clean, names(imp_ns)) else names(imp_ns)[!grepl("^MEM[0-9]+$", names(imp_ns))]
  env_scores <- if (length(env_keep) > 0L) imp_ns[env_keep] else numeric(0)
  imp_sp <- .extract_rf_perm_importance(rf_spatial)
  mem_names <- grep("^MEM[0-9]+$", names(imp_sp), value = TRUE)
  mem_scores <- if (length(mem_names) > 0L) imp_sp[mem_names] else numeric(0)
  all_scores <- c(env_scores, mem_scores)
  if (length(all_scores) == 0L) {
    return(data.frame(
      Predictor = character(0), MeanDecreaseAccuracy = numeric(0),
      PositiveImportance = numeric(0), Type = character(0), Group = character(0), stringsAsFactors = FALSE
    ))
  }
  pred_names <- names(all_scores); vals <- as.numeric(all_scores)
  is_mem_pred <- grepl("^MEM[0-9]+$", pred_names)
  is_fine_edaphic_pred <- grepl("^(fine_|ilr_texture_|ssurgo_|solus_|topo_|soil_|sgmc_|is_azonal|lithic_outcrop_|wetland_seepage_)", pred_names) &
    !grepl("^(polaris_|ultramafic_calcareous_)", pred_names)
  is_climate_pred <- grepl("^(bio|delta_doy|doy|summer_moisture_|maritime_continentality_)", pred_names, ignore.case = TRUE)
  imp_df <- data.frame(
    Predictor = pred_names, MeanDecreaseAccuracy = vals, PositiveImportance = pmax(0.0, vals),
    Type = ifelse(is_mem_pred, "Spatial (MEMs)", ifelse(is_fine_edaphic_pred, "Fine-Grain Edaphic & Outcrop (10-30m)", "Environmental")),
    Group = ifelse(is_mem_pred, "[C] Spatial MEMs",
                   ifelse(is_fine_edaphic_pred, "[A] Fine-Grain Edaphic & Outcrop Lithology",
                          ifelse(is_climate_pred, "[B] Macro-Climate", "[A] Fine-Grain Edaphic & Outcrop Lithology"))),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(imp_df$MeanDecreaseAccuracy, decreasing = TRUE), , drop = FALSE]
  rownames(imp_df) <- NULL
  imp_df
}
