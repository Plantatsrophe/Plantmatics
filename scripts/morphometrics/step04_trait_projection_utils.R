if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/step04_trait_projection_utils.R
# Step 04 Calibration Masking, Trait Pools, Shared Subspace & Concept Overrides
# ==============================================================================

mask_uncalibrated_step04_mm_cols <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  if (!any(c("is_calibrated_mm", "scale_source") %in% names(df))) return(df)
  cal_ok <- if (exists("is_step04_calibrated_mm", mode = "function")) is_step04_calibrated_mm(df) else rep(TRUE, nrow(df))
  if (length(cal_ok) != nrow(df) || all(cal_ok)) return(df)
  mm_cols <- c("log_centsize", "centroid_size", "area_mm2", "petiole_length_mm", "petiole_width_mm",
               "blade_length_mm", "blade_width_mm", "lamina_trace_length_mm", "cauline_blade_length_mm",
               "involucre_height_mm", "involucre_width_mm", "plant_height_mm", "plant_height_cm")
  for (mc in intersect(mm_cols, names(df))) { v <- suppressWarnings(as.numeric(df[[mc]])); v[!cal_ok] <- NA_real_; df[[mc]] <- v }
  df
}

prepare_step04_standardized_features <- function(closed_df, pca_cols, outline_only = FALSE, intraspecific_mode = FALSE) {
  pca_cols <- intersect(as.character(pca_cols), names(closed_df))
  if (is.null(closed_df) || !is.data.frame(closed_df) || nrow(closed_df) == 0L) return(list(df = closed_df, gmm_cols = pca_cols, cda_cols = pca_cols))
  n <- nrow(closed_df); closed_df <- mask_uncalibrated_step04_mm_cols(closed_df)
  if ("AR_inv" %in% names(closed_df)) closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$AR_inv))
  if (!("log1p_capitulum_count" %in% names(closed_df)) && any(c("capitula_count", "capitulum_count", "corymb_head_count_est") %in% names(closed_df))) {
    cc_col <- if ("capitula_count" %in% names(closed_df)) "capitula_count" else if ("capitulum_count" %in% names(closed_df)) "capitulum_count" else NULL
    cc_n <- if (!is.null(cc_col)) suppressWarnings(as.numeric(closed_df[[cc_col]])) else rep(NA_real_, n)
    ce_n <- if ("corymb_head_count_est" %in% names(closed_df)) suppressWarnings(as.numeric(closed_df$corymb_head_count_est)) else rep(NA_real_, n)
    cc_n[!is.finite(cc_n)] <- NA_real_; ce_n[!is.finite(ce_n)] <- NA_real_
    cc_eff <- ifelse(is.finite(cc_n) | is.finite(ce_n), pmax(cc_n, ce_n, na.rm = TRUE), NA_real_)
    closed_df$log1p_capitulum_count <- ifelse(is.finite(cc_eff) & cc_eff >= 1, round(log1p(cc_eff), 5), NA_real_)
  }
  phot_eb_cols <- unique(intersect(c(grep("^blade_cie_", names(closed_df), value = TRUE), "basal_zone_cie_l_delta", "abaxial_cie_l_sheet_delta",
                                     "cauline_cie_l_q85_sheet_delta", "involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta",
                                     "phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta"), names(closed_df)))
  if (length(phot_eb_cols) > 0L && exists("apply_eb_institution_batch_centering", mode = "function")) {
    closed_df <- apply_eb_institution_batch_centering(closed_df, target_cols = phot_eb_cols, min_inst_n = 10L)
  }
  aw_vec <- if ("authority_weight" %in% names(closed_df)) suppressWarnings(as.numeric(closed_df$authority_weight)) else rep(1.0, n)
  aw_vec[!is.finite(aw_vec)] <- 0.0; ref_mask <- aw_vec > 0.0; ref_idx <- which(ref_mask)
  if (length(ref_idx) < 3L) { ref_mask <- rep(TRUE, n); ref_idx <- seq_len(n) }

  caul_het_raw <- c("basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "cauline_cie_l_q85_sheet_delta", "basal_to_cauline_wool_delta",
                    "heteroblasty_efa_distance", "heteroblasty_ect_distance",
                    "cauline_dissection_index_med", "cauline_solidity_med", "cauline_ect_lobe_complexity", "Cauline_PC1", "Cauline_ECT_PC1")
  micro_repro_raw <- c("phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta", "involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta", "calyculus_flare_index")
  diag_phys_raw <- c("log_centsize", "petiole_length_mm", "lamina_trace_length_mm", "plant_height_cm")
  repro_arch_raw <- c("AR_inv", "radiate_head_fraction", "log1p_capitulum_count", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio", micro_repro_raw)
  topo_arch_raw <- c("ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "basal_cordate_index", "solidity", "q25_solidity", "dissection_index",
                     "terminal_lobe_area_frac", "margin_serration_index", "margin_tooth_acuteness", "base_angle_deg", "apex_angle_deg",
                     "lamina_curvature", "within_rosette_efa_dispersion", "aspect_ratio")
  phot_raw <- c("blade_cie_l_q85_sheet_delta", "blade_cie_l_sheet_delta", "basal_zone_cie_l_delta", "abaxial_cie_l_sheet_delta",
                "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta", "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_a_delta", "blade_cie_b_delta")
  for (col_nm in unique(intersect(c(repro_arch_raw, topo_arch_raw, phot_raw, caul_het_raw, diag_phys_raw), names(closed_df)))) {
    vv <- suppressWarnings(as.numeric(closed_df[[col_nm]])); closed_df[[col_nm]] <- vv
    vf <- which(is.finite(vv)); if (length(vf) < 5L) next
    vr <- intersect(vf, ref_idx); if (length(vr) < 3L) vr <- vf
    m_v <- mean(vv[vr]); s_v <- stats::sd(vv[vr]); sc_nm <- paste0(col_nm, "_std")
    if (is.finite(s_v) && s_v > 1e-6) {
      out_v <- rep(NA_real_, n); out_v[vf] <- round((vv[vf] - m_v) / s_v, 5); closed_df[[sc_nm]] <- out_v
    } else if (col_nm %in% c(repro_arch_raw, diag_phys_raw)) {
      out_v <- rep(NA_real_, n); out_v[vf] <- 0.0; closed_df[[sc_nm]] <- out_v
    }
  }
  if (isTRUE(outline_only)) return(list(df = closed_df, gmm_cols = pca_cols, cda_cols = pca_cols))

  has_cov <- function(col_nm, min_frac = 0.70) {
    if (!(col_nm %in% names(closed_df))) return(FALSE)
    vv <- suppressWarnings(as.numeric(closed_df[[col_nm]][ref_idx]))
    (mean(is.finite(vv)) >= min_frac) && (sum(is.finite(vv)) >= 5L) && isTRUE(stats::sd(vv, na.rm = TRUE) > 1e-6)
  }
  has_sheet_l <- has_cov("blade_cie_l_q85_sheet_delta_std", 0.70) || has_cov("blade_cie_l_sheet_delta_std", 0.70)
  caul_het_std <- paste0(c(caul_het_raw, micro_repro_raw), "_std")
  cand_pool <- unique(c(
    "AR_inv_std", "radiate_head_fraction_std", "log1p_capitulum_count_std", "petiole_to_blade_ratio_std", "scape_to_basal_leaf_ratio_std",
    "ECT_PC1_std", "ECT_PC2_std", "ECT_PC3_std", "ect_lobe_complexity_std", "basal_cordate_index_std",
    "solidity_std", "q25_solidity_std", "dissection_index_std", "terminal_lobe_area_frac_std",
    "margin_serration_index_std", "margin_tooth_acuteness_std", "base_angle_deg_std", "apex_angle_deg_std", "lamina_curvature_std",
    "within_rosette_efa_dispersion_std", "blade_cie_l_q85_sheet_delta_std", "blade_cie_l_sheet_delta_std",
    "basal_zone_cie_l_delta_std", "abaxial_cie_l_sheet_delta_std", "blade_cie_a_browning_delta_std", "blade_cie_b_yellowing_delta_std",
    if (!has_sheet_l) "blade_cie_l_norm_std", "aspect_ratio_std", caul_het_std
  ))
  cov_cands <- cand_pool[vapply(cand_pool, function(sc) has_cov(sc, if (sc %in% caul_het_std) 0.45 else 0.70), logical(1))]
  dis_var <- if (length(cov_cands) > 0L && exists("filter_biologically_variable_dissection_scalars", mode = "function")) {
    filter_biologically_variable_dissection_scalars(closed_df, cov_cands, active_mask = ref_mask)
  } else cov_cands

  pc_cols <- grep("^PC[0-9]+$", pca_cols, value = TRUE)
  if (length(pc_cols) >= 2L) {
    if (isTRUE(intraspecific_mode)) {
      intra_pool <- c("ECT_PC1_std", "ECT_PC2_std", "ECT_PC3_std", "petiole_to_blade_ratio_std",
                      "margin_tooth_acuteness_std", "basal_cordate_index_std", "apex_angle_deg_std",
                      "base_angle_deg_std", "blade_cie_l_q85_sheet_delta_std", "q25_solidity_std",
                      "phyllary_tip_darkening_delta_std", "phyllary_tip_anthocyanin_delta_std",
                      "involucre_basal_cie_l_delta_std", "peduncle_cie_l_sheet_delta_std", "basal_to_cauline_wool_delta_std")
      intra_cands <- intersect(intra_pool, dis_var)
      intra_cands <- intra_cands[vapply(intra_cands, function(sc) {
        vv <- round(suppressWarnings(as.numeric(closed_df[[sc]][ref_idx])), 5); vv <- vv[is.finite(vv)]
        length(vv) >= 5L && ((max(tabulate(match(vv, unique(vv)))) / length(vv)) <= 0.25)
      }, logical(1L))]
      if (length(intra_cands) > 0L && exists("screen_gmm_point_mass_columns", mode = "function")) {
        X_pm <- as.matrix(closed_df[, c(pc_cols[1:2], intra_cands), drop = FALSE])
        pm_k <- screen_gmm_point_mass_columns(X_pm, ref_mask = ref_mask, max_tie_frac = 0.25)
        intra_cands <- intersect(intra_cands, if (is.character(pm_k)) pm_k else colnames(X_pm)[pm_k])
      }
      if (length(intra_cands) > 0L && exists("prune_collinear_dissection_scalars", mode = "function")) {
        intra_cands <- prune_collinear_dissection_scalars(closed_df, base_cols = pc_cols, cand_scalars = intra_cands, active_mask = ref_mask, group_vec = NULL, max_r = 0.75)
      }
      gmm_cols <- unique(c(pc_cols, intra_cands))
    } else {
      gmm_cols <- unique(c(pc_cols, intersect(c("ECT_PC1_std", "ECT_PC2_std", "ECT_PC3_std"), dis_var)))
    }
  } else {
    gmm_cand <- intersect(c("solidity_std", "q25_solidity_std", "dissection_index_std", "ect_lobe_complexity_std"), dis_var)
    gmm_extra <- if (length(gmm_cand) > 0L && exists("prune_collinear_dissection_scalars", mode = "function")) {
      prune_collinear_dissection_scalars(closed_df, base_cols = pca_cols, cand_scalars = gmm_cand, active_mask = ref_mask, group_vec = NULL, max_r = 0.75)
    } else gmm_cand
    gmm_cols <- unique(c(pca_cols, gmm_extra))
  }
  cda_extra <- if (length(dis_var) > 0L && exists("prune_collinear_dissection_scalars", mode = "function")) {
    prune_collinear_dissection_scalars(closed_df, base_cols = pca_cols, cand_scalars = dis_var, active_mask = ref_mask, group_vec = NULL, max_r = 0.75)
  } else dis_var
  list(df = closed_df, gmm_cols = gmm_cols, cda_cols = unique(c(pca_cols, cda_extra)))
}

project_shared_subspace_step04 <- function(df, cda_res, gmm_res, pca_cols) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L || length(pca_cols) == 0L) return(list(cda_res = cda_res, gmm_res = gmm_res))
  pc_in_pca <- grep("^PC[0-9]+$", pca_cols, value = TRUE)
  guard_cols <- if (length(pc_in_pca) >= 2L) pc_in_pca[seq_len(2L)] else pca_cols[seq_len(min(2L, length(pca_cols)))]
  if (length(guard_cols) == 0L || !all(guard_cols %in% names(df))) return(list(cda_res = cda_res, gmm_res = gmm_res))
  no_efa_idx <- which(!stats::complete.cases(df[, guard_cols, drop = FALSE]))
  if (length(no_efa_idx) == 0L) return(list(cda_res = cda_res, gmm_res = gmm_res))

  n <- nrow(df); cda_res$is_active[no_efa_idx] <- FALSE
  if (!is.null(cda_res$X_proj) && is.matrix(cda_res$X_proj) && nrow(cda_res$X_proj) == n) {
    pc_xp <- grep("^PC[0-9]+$", colnames(cda_res$X_proj))
    if (length(pc_xp) > 0L) cda_res$X_proj[no_efa_idx, pc_xp] <- NA_real_
  }
  if (is.null(cda_res$self_morph_mahal_d2)) cda_res$self_morph_mahal_d2 <- rep(NA_real_, n)
  if (is.null(cda_res$self_morph_df)) cda_res$self_morph_df <- rep(NA_integer_, n)

  base_topo <- c("ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "solidity", "q25_solidity", "dissection_index",
                 "basal_cordate_index", "terminal_lobe_area_frac", "margin_serration_index", "margin_tooth_acuteness",
                 "base_angle_deg", "apex_angle_deg", "lamina_curvature", "aspect_ratio", "AR_inv", "radiate_head_fraction",
                 "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio", "blade_cie_l_q85_sheet_delta",
                 "blade_cie_l_sheet_delta", "basal_zone_cie_l_delta")
  act_idx <- which(as.logical(cda_res$is_active)); shared_cols <- character(0)
  if (length(act_idx) >= 4L) {
    for (bt in base_topo) {
      for (cn in intersect(c(paste0(bt, "_std"), bt), names(df))) {
        v_num <- suppressWarnings(as.numeric(df[[cn]]))
        if (sum(is.finite(v_num[act_idx])) >= max(3L, ceiling(0.50 * length(act_idx))) &&
            any(is.finite(v_num[no_efa_idx])) && isTRUE(stats::sd(v_num[act_idx], na.rm = TRUE) > 1e-6)) {
          shared_cols <- c(shared_cols, cn); break
        }
      }
    }
  }
  shared_cols <- unique(shared_cols)

  clear_unproj <- function(idx_v) {
    if (length(idx_v) == 0L) return(invisible(NULL))
    cda_res$predicted_taxon[idx_v] <<- NA_character_; cda_res$posterior_prob[idx_v] <<- NA_real_
    cda_res$second_taxon[idx_v] <<- NA_character_; cda_res$second_posterior_prob[idx_v] <<- NA_real_
    cda_res$posterior_margin[idx_v] <<- NA_real_
    if (!is.null(cda_res$canonical_scores)) cda_res$canonical_scores[idx_v, ] <<- NA_real_
    if (!is.null(cda_res$oof_posterior_prob)) cda_res$oof_posterior_prob[idx_v] <<- NA_real_
    cda_res$self_morph_mahal_d2[idx_v] <<- NA_real_; cda_res$self_morph_df[idx_v] <<- NA_integer_
  }
  if (length(shared_cols) < 2L || length(act_idx) < 4L) { clear_unproj(no_efa_idx); return(list(cda_res = cda_res, gmm_res = gmm_res)) }

  X_sh_raw <- as.matrix(df[, shared_cols, drop = FALSE]); storage.mode(X_sh_raw) <- "double"
  obs_counts <- as.integer(rowSums(is.finite(X_sh_raw))); p_sh <- ncol(X_sh_raw)
  mu_sh <- vapply(seq_len(p_sh), function(j) mean(X_sh_raw[act_idx, j], na.rm = TRUE), numeric(1))
  sd_sh <- vapply(seq_len(p_sh), function(j) { s <- stats::sd(X_sh_raw[act_idx, j], na.rm = TRUE); if (is.finite(s) && s > 1e-6) s else 1.0 }, numeric(1))
  X_sh <- sweep(sweep(X_sh_raw, 2L, mu_sh, "-"), 2L, sd_sh, "/")
  for (j in seq_len(p_sh)) {
    bad_j <- !is.finite(X_sh[, j])
    if (any(bad_j)) { med_j <- stats::median(X_sh[act_idx, j], na.rm = TRUE); X_sh[bad_j, j] <- if (is.finite(med_j)) med_j else 0.0 }
  }

  grp_all <- if (!is.null(cda_res$grouping_col)) as.character(cda_res$grouping_col) else if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Morphotype_K1", n)
  y_act <- grp_all[act_idx]
  g_lvls <- if (!is.null(cda_res$post_probs) && !is.null(colnames(cda_res$post_probs))) colnames(cda_res$post_probs) else sort(unique(stats::na.omit(y_act)))
  if (length(g_lvls) == 0L) g_lvls <- "Morphotype_K1"
  w_sb_sh <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight[act_idx])) else rep(1.0, length(act_idx))
  w_sb_sh[!is.finite(w_sb_sh) | w_sb_sh <= 0] <- 1.0
  w_act <- (if (!is.null(cda_res$authority_weight)) as.numeric(cda_res$authority_weight[act_idx]) else rep(1.0, length(act_idx))) * w_sb_sh
  w_act[!is.finite(w_act) | w_act <= 0.0] <- 1.0; G_sh <- length(g_lvls)
  lw_fn <- function(Sw, Xc) {
    S_r <- if (exists("compute_ledoit_wolf_cov", mode = "function")) compute_ledoit_wolf_cov(Sw, X_centered = Xc, n_obs = length(act_idx), weights = w_act) else (Sw + diag(1e-4, p_sh))
    if (is.list(S_r) && !is.null(S_r$S_reg)) S_r <- S_r$S_reg
    S_r
  }

  gm_sh <- colSums(X_sh[act_idx, , drop = FALSE] * w_act) / sum(w_act)
  if (G_sh >= 2L) {
    mu_g <- matrix(0.0, nrow = G_sh, ncol = p_sh, dimnames = list(g_lvls, shared_cols))
    X_cent_sh <- matrix(0.0, nrow = length(act_idx), ncol = p_sh)
    W_sh <- matrix(0.0, nrow = p_sh, ncol = p_sh); B_sh <- matrix(0.0, nrow = p_sh, ncol = p_sh)
    neff_g <- rep(1.0, G_sh); w_corr_sh <- 0.0
    for (ki in seq_along(g_lvls)) {
      km <- which(y_act == g_lvls[ki])
      if (length(km) > 0L) {
        sub_k <- X_sh[act_idx[km], , drop = FALSE]; wk <- w_act[km]; swk <- sum(wk)
        mk <- colSums(sub_k * wk) / swk; mu_g[ki, ] <- mk; dm <- mk - gm_sh; B_sh <- B_sh + swk * tcrossprod(dm)
        dx <- sweep(sub_k, 2L, mk, "-"); X_cent_sh[km, ] <- dx; W_sh <- W_sh + crossprod(dx, dx * wk)
        neff_g[ki] <- max(1.0, (swk^2) / max(sum(wk^2), 1e-12)); w_corr_sh <- w_corr_sh + (sum(wk^2) / swk)
      } else w_corr_sh <- w_corr_sh + 1.0
    }
    S_reg_sh <- lw_fn(W_sh / max(sum(w_act) - w_corr_sh, 1.0), X_cent_sh)
    S_inv_sh <- tryCatch(solve(S_reg_sh), error = function(e) diag(1.0 / pmax(diag(S_reg_sh), 1e-4), p_sh))
    eig_sh <- if (exists("solve_shrinkage_cda_eigen", mode = "function")) tryCatch(solve_shrinkage_cda_eigen(S_reg_sh, B_sh, max_axes = min(p_sh, G_sh - 1L)), error = function(e) NULL) else NULL
    if (!is.null(eig_sh) && !is.null(cda_res$centroids_full) && !is.null(rownames(cda_res$centroids_full)) && all(g_lvls %in% rownames(cda_res$centroids_full))) {
      c_sh <- sweep(mu_g, 2L, gm_sh, "-") %*% eig_sh$vectors
      for (aj in seq_len(min(ncol(c_sh), ncol(cda_res$centroids_full)))) {
        r_al <- suppressWarnings(stats::cor(c_sh[, aj], cda_res$centroids_full[g_lvls, aj], use = "complete.obs"))
        if (is.finite(r_al) && r_al < 0.0) eig_sh$vectors[, aj] <- -eig_sh$vectors[, aj]
      }
    }
    for (ri in no_efa_idx) {
      if (obs_counts[ri] < 2L) { clear_unproj(ri); next }
      lp <- numeric(G_sh); d2_sh_vec <- numeric(G_sh)
      for (ki in seq_along(g_lvls)) { dv <- X_sh[ri, ] - mu_g[ki, ]; d2_sh_vec[ki] <- max(0.0, as.numeric(t(dv) %*% S_inv_sh %*% dv)); lp[ki] <- log(1.0 / G_sh) - 0.5 * d2_sh_vec[ki] }
      pp <- exp(lp - max(lp)); pp <- pp / max(sum(pp), 1e-12)
      if (!is.null(cda_res$post_probs) && ncol(cda_res$post_probs) == G_sh) cda_res$post_probs[ri, ] <- pp
      best_k <- which.max(pp); oth_k <- setdiff(seq_len(G_sh), best_k); sec_k <- oth_k[which.max(pp[oth_k])]
      cda_res$predicted_taxon[ri] <- g_lvls[best_k]; cda_res$posterior_prob[ri] <- pp[best_k]
      cda_res$second_taxon[ri] <- g_lvls[sec_k]; cda_res$second_posterior_prob[ri] <- pp[sec_k]; cda_res$posterior_margin[ri] <- pp[best_k] - pp[sec_k]
      if (!is.null(cda_res$oof_posterior_prob)) cda_res$oof_posterior_prob[ri] <- pp[best_k]
      if (!is.null(cda_res$cl_class_map) && (g_lvls[best_k] %in% names(cda_res$cl_class_map))) cda_res$cl_target_threshold[ri] <- unname(cda_res$cl_class_map[g_lvls[best_k]])
      yi_sh <- match(grp_all[ri], g_lvls)
      cda_res$self_morph_mahal_d2[ri] <- if (!is.na(yi_sh)) round(d2_sh_vec[yi_sh] / (1.0 + 1.0 / max(1.0, neff_g[yi_sh])), 4) else NA_real_
      cda_res$self_morph_df[ri] <- as.integer(p_sh)
      if (!is.null(eig_sh) && !is.null(cda_res$canonical_scores)) {
        z_r <- as.numeric((X_sh[ri, ] - gm_sh) %*% eig_sh$vectors); n_fill <- min(length(z_r), ncol(cda_res$canonical_scores))
        cda_res$canonical_scores[ri, ] <- 0.0; cda_res$canonical_scores[ri, seq_len(n_fill)] <- z_r[seq_len(n_fill)]
      }
    }
  } else {
    Xc_sh <- sweep(X_sh[act_idx, , drop = FALSE], 2L, gm_sh, "-")
    denom_1 <- max(sum(w_act) - (sum(w_act^2) / max(sum(w_act), 1e-12)), 1.0)
    S_reg_sh <- lw_fn(crossprod(Xc_sh, Xc_sh * w_act) / denom_1, Xc_sh)
    S_inv_sh <- tryCatch(solve(S_reg_sh), error = function(e) diag(1.0 / pmax(diag(S_reg_sh), 1e-4), p_sh))
    eig_1 <- eigen(0.5 * (S_reg_sh + t(S_reg_sh)), symmetric = TRUE); V_1 <- eig_1$vectors[, seq_len(min(3L, p_sh)), drop = FALSE]
    if (!is.null(cda_res$canonical_scores) && nrow(cda_res$canonical_scores) == n) {
      Z_act_sh <- Xc_sh %*% V_1
      for (aj in seq_len(min(ncol(Z_act_sh), ncol(cda_res$canonical_scores)))) {
        r_al <- suppressWarnings(stats::cor(Z_act_sh[, aj], cda_res$canonical_scores[act_idx, aj], use = "complete.obs"))
        if (is.finite(r_al) && r_al < 0.0) V_1[, aj] <- -V_1[, aj]
      }
    }
    neff_1 <- max(1.0, (sum(w_act)^2) / max(sum(w_act^2), 1e-12))
    for (ri in no_efa_idx) {
      if (obs_counts[ri] < 2L) { clear_unproj(ri); next }
      dv <- X_sh[ri, ] - gm_sh; raw_d2 <- max(0.0, as.numeric(t(dv) %*% S_inv_sh %*% dv))
      cda_res$predicted_taxon[ri] <- g_lvls[1L]; cda_res$posterior_prob[ri] <- 1.0
      cda_res$second_taxon[ri] <- NA_character_; cda_res$second_posterior_prob[ri] <- 0.0; cda_res$posterior_margin[ri] <- 1.0
      if (!is.null(cda_res$post_probs) && ncol(cda_res$post_probs) == 1L) cda_res$post_probs[ri, 1L] <- 1.0
      if (!is.null(cda_res$oof_posterior_prob)) cda_res$oof_posterior_prob[ri] <- 1.0
      cda_res$self_morph_mahal_d2[ri] <- round(raw_d2 / (1.0 + 1.0 / neff_1), 4); cda_res$self_morph_df[ri] <- as.integer(p_sh)
      if (!is.null(cda_res$canonical_scores)) {
        z_r <- as.numeric(dv %*% V_1); n_fill <- min(length(z_r), ncol(cda_res$canonical_scores))
        cda_res$canonical_scores[ri, ] <- 0.0; cda_res$canonical_scores[ri, seq_len(n_fill)] <- z_r[seq_len(n_fill)]
      }
    }
  }

  if (!is.null(gmm_res) && !is.null(gmm_res$classification)) {
    cl_act <- gmm_res$classification[act_idx]; u_cl <- sort(unique(stats::na.omit(cl_act)))
    u_cl <- u_cl[vapply(u_cl, function(ck) any(!is.na(cl_act) & cl_act == ck), logical(1))]
    if (length(u_cl) >= 1L) {
      cl_meds <- do.call(rbind, lapply(u_cl, function(ck) apply(X_sh[act_idx[!is.na(cl_act) & cl_act == ck], , drop = FALSE], 2L, stats::median, na.rm = TRUE)))
      for (ri in no_efa_idx) if (obs_counts[ri] >= 2L) {
        d2_c <- rowSums(sweep(cl_meds, 2L, X_sh[ri, ], "-")^2); wc <- exp(-0.5 * (d2_c - min(d2_c))); wc <- wc / max(sum(wc), 1e-12)
        gmm_res$classification[ri] <- u_cl[which.max(wc)]
        if (!is.null(gmm_res$uncertainty)) gmm_res$uncertainty[ri] <- 1.0 - max(wc)
      }
    }
  }
  cda_res$shared_subspace <- list(no_efa_idx = no_efa_idx, act_idx = act_idx, shared_cols = shared_cols,
                                  X_sh = X_sh, obs_counts = obs_counts, S_reg_sh = S_reg_sh, S_inv_sh = S_inv_sh)
  list(cda_res = cda_res, gmm_res = gmm_res)
}

apply_open_or_hybrid_concept_audit_overrides <- function(df, cda_res, flags, triage, reasons, misid_flags, p_audit) {
  n <- if (is.null(df) || !is.data.frame(df)) 0L else nrow(df)
  if (n == 0L) return(list(cda_res = cda_res, flags = flags, triage = triage, reasons = reasons, misid_flags = misid_flags))
  is_type_v <- if (exists(".identify_type_specimen_vec", mode = "function")) as.logical(.identify_type_specimen_vec(df)$is_type_specimen) else if ("is_type_specimen" %in% names(df)) (!is.na(df$is_type_specimen) & as.logical(df$is_type_specimen)) else rep(FALSE, n)
  is_type_v[is.na(is_type_v)] <- FALSE
  oh_mask <- if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) as.logical(.identify_open_or_hybrid_concept_mask(df)) else rep(FALSE, n)
  oh_mask[is.na(oh_mask)] <- FALSE; oh_idx <- which(!is_type_v & oh_mask)
  if (length(oh_idx) == 0L) return(list(cda_res = cda_res, flags = flags, triage = triage, reasons = reasons, misid_flags = misid_flags))

  if (is.null(cda_res$anchor_tier_source) || length(cda_res$anchor_tier_source) != n) cda_res$anchor_tier_source <- rep("Passive_Projected", n)
  cda_res$anchor_tier_source[oh_idx] <- "Passive_Open_or_Hybrid_Concept"
  hyb_pat <- "(^|[[:space:]])[a-z][a-z-]+([[:space:]]*\u00d7[[:space:]]*|[[:space:]]+[xX][[:space:]]+)([A-Z][a-z]*\\.?[[:space:]]+)?[a-z][a-z-]+"
  hyb_kw <- "putative[[:space:]]+hybrid|hybrid[[:space:]]+(with|between)"
  txt_cols <- intersect(c("verbatim_taxon_concept", "scientificName", "species_raw", "verbatimIdentification", "identificationRemarks"), names(df))
  has_cs <- "concept_synonymy_status" %in% names(df)

  for (i in oh_idx) {
    is_hyb <- has_cs && !is.na(df$concept_synonymy_status[i]) && identical(trimws(as.character(df$concept_synonymy_status[i])), "Hybrid_Formula")
    if (!is_hyb && length(txt_cols) > 0L) for (tc in txt_cols) {
      tv <- trimws(as.character(df[[tc]][i]))
      if (!is.na(tv) && nzchar(tv) && (grepl(hyb_pat, tv) || grepl(hyb_kw, tv, ignore.case = TRUE))) { is_hyb <- TRUE; break }
    }
    pred_sp <- if (!is.null(cda_res$predicted_taxon) && length(cda_res$predicted_taxon) >= i && !is.na(cda_res$predicted_taxon[i])) trimws(as.character(cda_res$predicted_taxon[i])) else ""
    sec_sp <- if (!is.null(cda_res$second_taxon) && length(cda_res$second_taxon) >= i && !is.na(cda_res$second_taxon[i])) trimws(as.character(cda_res$second_taxon[i])) else ""
    post_p <- if (!is.null(cda_res$posterior_prob) && length(cda_res$posterior_prob) >= i) suppressWarnings(as.numeric(cda_res$posterior_prob[i])) else NA_real_
    raw_sp <- if (!is.null(cda_res$grouping_col) && length(cda_res$grouping_col) >= i && !is.na(cda_res$grouping_col[i]) && nzchar(trimws(as.character(cda_res$grouping_col[i])))) trimws(as.character(cda_res$grouping_col[i])) else if ("species_standardized" %in% names(df) && !is.na(df$species_standardized[i])) trimws(as.character(df$species_standardized[i])) else ""

    if (is_hyb) {
      flags[i] <- FALSE; misid_flags[i] <- FALSE
      triage[i] <- if (is.finite(post_p) && post_p < 0.80) "HIGH" else "MEDIUM"
      pred_lbl <- if (nzchar(pred_sp)) gsub(" ", "_", pred_sp) else "Unassigned"
      reasons[i] <- if (is.finite(post_p) && post_p < 0.80 && nzchar(sec_sp)) sprintf("Hybrid_Concept_Intermediate_%s_%s", pred_lbl, gsub(" ", "_", sec_sp)) else sprintf("Hybrid_Concept_Asymmetric_Affinity_to_%s", pred_lbl)
    } else {
      is_morphotype_grp <- grepl("^Morphotype_K[0-9]+$", pred_sp) || grepl("^Morphotype_K[0-9]+$", raw_sp) || isTRUE(cda_res$is_intraspecific_mode)
      self_d2 <- if (!is.null(cda_res$self_morph_mahal_d2) && length(cda_res$self_morph_mahal_d2) >= i) suppressWarnings(as.numeric(cda_res$self_morph_mahal_d2[i])) else NA_real_
      self_df <- if (!is.null(cda_res$self_morph_df) && length(cda_res$self_morph_df) >= i && is.finite(cda_res$self_morph_df[i])) as.integer(cda_res$self_morph_df[i]) else as.integer(p_audit)
      if (!is_morphotype_grp && nzchar(pred_sp) && raw_sp != pred_sp && is.finite(post_p) && post_p >= 0.90 && is.finite(self_d2) && self_d2 > stats::qchisq(0.99, df = max(1L, self_df))) {
        flags[i] <- TRUE; misid_flags[i] <- TRUE; triage[i] <- "HIGH"
        reasons[i] <- sprintf("Open_Nomenclature_Reassignment_Candidate_to_%s", gsub(" ", "_", pred_sp))
      } else {
        flags[i] <- FALSE; misid_flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Open_Nomenclature_Passive_Audit"
      }
    }
  }
  list(cda_res = cda_res, flags = flags, triage = triage, reasons = reasons, misid_flags = misid_flags)
}
