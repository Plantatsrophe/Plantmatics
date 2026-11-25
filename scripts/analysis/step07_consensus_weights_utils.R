if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists(".resolve_step07_anchor_mask", envir = environment(), mode = "function")) {
  .hc_env <- environment()
  .hc_prev <- exists(".step07_sourcing_consensus_weights", envir = .hc_env, inherits = FALSE)
  assign(".step07_sourcing_consensus_weights", TRUE, envir = .hc_env)
  tryCatch({
    .hc_dirs <- dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
                          unlist(lapply(sys.frames(), function(fr) fr$ofile))))
    .hc_dirs <- .hc_dirs[!is.na(.hc_dirs) & nzchar(.hc_dirs)]
    for (.hc_p in unique(c("scripts/analysis/step07_hopkins_consensus_utils.R", file.path(.hc_dirs, "step07_hopkins_consensus_utils.R")))) {
      if (!is.na(.hc_p) && nzchar(.hc_p) && file.exists(.hc_p)) { source(.hc_p, local = .hc_env); break }
    }
  }, finally = {
    if (!.hc_prev && exists(".step07_sourcing_consensus_weights", envir = .hc_env, inherits = FALSE)) {
      rm(".step07_sourcing_consensus_weights", envir = .hc_env)
    }
  })
}

compute_cohens_kappa <- function(y_true, y_pred, default_kappa = 0.70) {
  valid <- !is.na(y_true) & !is.na(y_pred) & nzchar(as.character(y_true)) & nzchar(as.character(y_pred))
  if (sum(valid) < 4L) return(max(0.05, default_kappa))
  yt <- as.character(y_true[valid]); yp <- as.character(y_pred[valid]); lvls <- sort(unique(c(yt, yp)))
  if (length(lvls) < 2L || length(unique(yt)) < 2L) return(max(0.05, min(0.95, mean(yt == yp))))
  tab <- table(factor(yt, levels = lvls), factor(yp, levels = lvls)); n_tot <- sum(tab)
  p_o <- sum(diag(tab)) / n_tot; p_e <- sum((rowSums(tab) / n_tot) * (colSums(tab) / n_tot))
  if (!is.finite(p_e) || (1.0 - p_e) <= 1e-8) return(max(0.05, p_o))
  kv <- (p_o - p_e) / (1.0 - p_e); if (!is.finite(kv)) return(max(0.05, default_kappa))
  max(0.05, min(0.99, round(kv, 4)))
}

.compute_oof_class_precision_map <- function(y_true_vec, y_pred_vec, prior_alpha = 0.5) {
  ok <- !is.na(y_true_vec) & !is.na(y_pred_vec) & nzchar(y_true_vec) & nzchar(y_pred_vec) & (y_true_vec != "Unknown") & (y_pred_vec != "Unknown")
  if (sum(ok) < 6L) return(numeric(0))
  yt <- as.character(y_true_vec[ok]); yp <- as.character(y_pred_vec[ok]); true_taxa <- sort(unique(yt)); all_taxa <- sort(unique(c(true_taxa, yp)))
  if (length(true_taxa) < 2L || length(all_taxa) < 2L) return(numeric(0))
  K_cols <- length(all_taxa); C_mat <- matrix(0.0, nrow = length(true_taxa), ncol = K_cols, dimnames = list(true_taxa, all_taxa))
  for (k_tx in true_taxa) { idx_k <- which(yt == k_tx); n_k <- length(idx_k); for (j_tx in all_taxa) C_mat[k_tx, j_tx] <- (sum(yp[idx_k] == j_tx) + prior_alpha / K_cols) / (n_k + prior_alpha) }
  prec_out <- stats::setNames(rep(1.0, K_cols), all_taxa)
  for (j_tx in all_taxa) {
    cs <- sum(C_mat[, j_tx]); dv <- if (j_tx %in% rownames(C_mat)) C_mat[j_tx, j_tx] else (prior_alpha / K_cols)
    prec_out[[j_tx]] <- min(1.0, max(0.05, round(if (is.finite(cs) && cs > 1e-12) dv / cs else 1.0, 4)))
  }
  prec_out
}

.resolve_step07_anthesis_flag <- function(df, n = nrow(df)) {
  if ("is_active_anthesis" %in% names(df)) { v <- as.logical(df$is_active_anthesis); v[is.na(v)] <- TRUE; v } else resolve_active_anthesis_mask(df)
}

compute_modality_reliability_weights <- function(df) {
  n <- nrow(df); if (n == 0L) return(c(morph = 0.80, edaphic = 0.65, pheno = 0.55))
  sub_a <- df[.resolve_step07_anchor_mask(df), , drop = FALSE]; na_s <- nrow(sub_a); y_true <- as.character(sub_a$species_standardized)
  n_true_taxa <- length(unique(y_true[!is.na(y_true) & nzchar(y_true) & y_true != "Unknown"]))
  m_pure <- if ("morph_tier_source" %in% names(sub_a)) (is.na(sub_a$morph_tier_source) | sub_a$morph_tier_source == "Basal_Rosette") else rep(TRUE, na_s)
  y_morph <- .chr_col(sub_a, "cda_predicted_taxon", na_s)
  def_morph <- if ("loocv_accuracy" %in% names(df) && any(is.finite(as.numeric(df$loocv_accuracy)))) stats::median(as.numeric(df$loocv_accuracy), na.rm = TRUE) else 0.80
  y_edaph <- if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else .chr_col(sub_a, "rf_predicted_taxon", na_s)
  cv_med <- if ("spatial_block_cv_acc" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$spatial_block_cv_acc))))) stats::median(suppressWarnings(as.numeric(df$spatial_block_cv_acc)), na.rm = TRUE) else 0.65
  if (n_true_taxa < 2L) {
    out_m_a <- if ("morphometric_outlier_flag" %in% names(sub_a)) as.logical(sub_a$morphometric_outlier_flag[m_pure]) else rep(FALSE, sum(m_pure)); out_m_a[is.na(out_m_a)] <- FALSE
    w_morph <- round(min(0.80, max(0.30, (if (length(out_m_a) > 0L) mean(!out_m_a) else 0.80) * 0.65)), 4)
    eco_reg_a <- if ("ecotype_regime" %in% names(sub_a)) as.character(sub_a$ecotype_regime) else .chr_col(sub_a, "analysis_group", na_s)
    rf_reg_a <- if ("rf_predicted_taxon" %in% names(sub_a)) as.character(sub_a$rf_predicted_taxon) else .chr_col(sub_a, "edaphic_best_fit_taxon", na_s)
    u_eco_reg <- unique(eco_reg_a[!is.na(eco_reg_a) & nzchar(eco_reg_a) & eco_reg_a != "Unknown"])
    w_edaph <- if (length(u_eco_reg) >= 2L && any(rf_reg_a %in% u_eco_reg, na.rm = TRUE)) max(0.15, min(0.85, compute_cohens_kappa(eco_reg_a, rf_reg_a, default_kappa = 0.50))) else round(min(0.75, max(0.20, 2.0 * (cv_med - 0.50))), 4)
  } else {
    w_morph <- max(0.05, compute_cohens_kappa(y_true[m_pure], y_morph[m_pure], default_kappa = def_morph))
    w_edaph <- max(0.05, compute_cohens_kappa(y_true, y_edaph, default_kappa = cv_med))
  }
  anth_a <- .resolve_step07_anthesis_flag(sub_a, na_s); hop_col_a <- intersect(c("hopkins_adjusted_doy", "doy"), names(sub_a))[1L]
  doy_a <- if (!is.na(hop_col_a)) .num_col(sub_a, hop_col_a, na_s) else rep(NA_real_, na_s)
  ddoy_a <- .coalesce_step07_ddoy_resid(sub_a, na_s)
  ok_ph_cls <- which(anth_a & is.finite(doy_a) & doy_a >= 1.0 & doy_a <= 366.0 & !is.na(y_true) & nzchar(y_true) & y_true != "Unknown")
  tx_ph_valid <- if (length(ok_ph_cls) >= 4L) names(which(table(y_true[ok_ph_cls]) >= 2L)) else character(0)
  ok_ph_cls <- ok_ph_cls[y_true[ok_ph_cls] %in% tx_ph_valid]
  w_pheno <- if (length(tx_ph_valid) >= 2L && length(ok_ph_cls) >= 4L) {
    th_ph <- 2.0 * pi * (doy_a[ok_ph_cls] / 365.25); yt_ph <- y_true[ok_ph_cls]; yp_ph <- character(length(ok_ph_cls))
    for (ii in seq_along(ok_ph_cls)) {
      best_d <- Inf; best_tx <- tx_ph_valid[1L]
      for (tx_k in tx_ph_valid) {
        jj <- which(yt_ph == tx_k & seq_along(ok_ph_cls) != ii); if (length(jj) == 0L) jj <- which(yt_ph == tx_k)
        mu_k <- atan2(mean(sin(th_ph[jj])), mean(cos(th_ph[jj]))); r_k <- max(1e-4, min(0.9995, sqrt(mean(sin(th_ph[jj]))^2 + mean(cos(th_ph[jj]))^2)))
        d_z <- abs(atan2(sin(th_ph[ii] - mu_k), cos(th_ph[ii] - mu_k))) / max(2.0 * pi * (6.0 / 365.25), sqrt(-2.0 * log(r_k)))
        if (d_z < best_d) { best_d <- d_z; best_tx <- tx_k }
      }
      yp_ph[ii] <- best_tx
    }
    max(0.05, min(0.95, compute_cohens_kappa(yt_ph, yp_ph, default_kappa = 0.15)))
  } else if (sum(is.finite(ddoy_a)) >= 4L) 0.20 else 0.15
  w_geno <- if (exists("compute_genomic_kappa_weight", mode = "function")) compute_genomic_kappa_weight(df) else NA_real_
  if (is.finite(w_geno)) c(morph = w_morph, edaphic = w_edaph, pheno = w_pheno, geno = w_geno) else c(morph = w_morph, edaphic = w_edaph, pheno = w_pheno)
}

compute_probabilistic_consensus_scores <- function(df) {
  n <- nrow(df); if (n == 0L) return(df)
  prop_sum_attr <- attr(df, "genomic_propagation_summary"); w_vec <- compute_modality_reliability_weights(df)
  df$w_morph_kappa <- w_vec[["morph"]]; df$w_edaphic_kappa <- w_vec[["edaphic"]]; df$w_pheno_kappa <- w_vec[["pheno"]]
  has_geno_row <- .has_step07_geno_row_mask(df); w_geno_base <- if ("geno" %in% names(w_vec) && is.finite(w_vec[["geno"]])) w_vec[["geno"]] else NA_real_
  sub_a <- df[.resolve_step07_anchor_mask(df), , drop = FALSE]
  y_true_a <- .chr_col(sub_a, "species_standardized", nrow(sub_a))
  m_pure_a <- if ("morph_tier_source" %in% names(sub_a)) (is.na(sub_a$morph_tier_source) | sub_a$morph_tier_source == "Basal_Rosette") else rep(TRUE, nrow(sub_a))
  y_morph_a <- .chr_col(sub_a, "cda_predicted_taxon", nrow(sub_a))
  y_edaph_a <- if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else .chr_col(sub_a, "rf_predicted_taxon", nrow(sub_a))
  prec_morph <- .compute_oof_class_precision_map(y_true_a[m_pure_a], y_morph_a[m_pure_a]); prec_edaph <- .compute_oof_class_precision_map(y_true_a, y_edaph_a)
  lookup_prec <- function(pmap, pred_tx) if (length(pmap) == 0L || is.na(pred_tx) || !nzchar(pred_tx) || !(pred_tx %in% names(pmap))) 1.0 else as.numeric(pmap[[pred_tx]])

  p_conf_out <- numeric(n); p_score_out <- numeric(n); n_mod_out <- integer(n); df_nms_c <- names(df)
  w_morph_sx <- numeric(n); w_edaph_sx <- numeric(n); w_pheno_sx <- numeric(n); w_geno_sx <- numeric(n)
  given_sp_cv <- as.character(df$species_standardized)
  morph_sp_cv <- .chr_col(df, "cda_predicted_taxon", n); morph_sp_cv[!is.na(morph_sp_cv) & !nzchar(morph_sp_cv)] <- NA_character_
  cda_p_cv <- .num_col(df, "cda_posterior_prob", n); is_mo_cv <- .lgl_col(df, "morphometric_outlier_flag", n); has_geo_cv <- is.finite(.num_col(df, "decimalLatitude", n))
  edaph_sp_cv <- if ("edaphic_best_fit_taxon" %in% df_nms_c) as.character(df$edaphic_best_fit_taxon) else .chr_col(df, "rf_predicted_taxon", n)
  edaph_sp_cv[!is.na(edaph_sp_cv) & !nzchar(edaph_sp_cv)] <- NA_character_
  rf_p_cv <- .num_col(df, "rf_predicted_prob", n); rf_g_cv <- .num_col(df, "rf_given_taxon_prob", n)
  ez_cv <- .num_col(df, "edaphic_z_score", n); tez_cv <- .num_col(df, "type_edaphic_z_score", n); ddoy_r_cv <- .coalesce_step07_ddoy_resid(df, n)
  ded_cv <- if ("edaphic_niche_overlap_d" %in% df_nms_c) .num_col(df, "edaphic_niche_overlap_d", n) else .num_col(df, "edaphic_schoener_d", n)
  is_anth_cv <- if ("is_active_anthesis" %in% df_nms_c) .lgl_col(df, "is_active_anthesis", n) else rep(TRUE, n)
  w_morph_base <- w_vec[["morph"]]; w_edaph_base <- w_vec[["edaphic"]]; w_pheno_base <- w_vec[["pheno"]]
  is_single_sp_cons <- length(unique(given_sp_cv[!is.na(given_sp_cv) & nzchar(given_sp_cv) & given_sp_cv != "Unknown"])) < 2L
  eco_reg_cv <- if ("ecotype_regime" %in% df_nms_c) as.character(df$ecotype_regime) else .chr_col(df, "analysis_group", n)
  sd2_cv <- .num_col(df, "self_morph_mahal_d2", n); sdf_cv <- .num_col(df, "self_morph_df", n, 6.0); sdf_cv[!is.finite(sdf_cv) | sdf_cv < 1.0] <- 6.0
  geno_sp_cv <- .chr_col(df, "geno_predicted_taxon", n); geno_p_cv <- .num_col(df, "geno_posterior_prob", n); geno_qmax_cv <- .num_col(df, "q_max", n); deme_w_cv <- .num_col(df, "geno_deme_weight_scale", n)
  geno_tier_cv <- if ("geno_linkage_tier" %in% df_nms_c) as.character(df$geno_linkage_tier) else .chr_col(df, "linkage_tier", n)
  geno_adm_cv <- if ("is_geno_admixed" %in% df_nms_c) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
  geno_qry_cv <- if ("is_putative_novel_or_hybrid_query" %in% df_nms_c) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  geno_link_scale <- ifelse(geno_tier_cv == "G1_Direct_Voucher", 1.00, ifelse(geno_tier_cv == "G1_Exsiccatae_Isoduplicate", 0.95,
                     ifelse(geno_tier_cv == "G2_Population_Voucher_Proxy", 0.90, ifelse(geno_tier_cv == "G2_Conditional_Deme_Hull", ifelse(is.finite(deme_w_cv) & deme_w_cv > 0.0, pmin(0.75, pmax(0.45, deme_w_cv)), 0.75), 0.00))))
  geno_link_scale[is.na(geno_link_scale)] <- 0.0; run_has_geno <- any(has_geno_row) && is.finite(w_geno_base)
  df$w_geno_kappa <- if (run_has_geno) ifelse(has_geno_row, round(w_geno_base * geno_link_scale, 4), 0.0) else NA_real_
  t0_abstain_cv <- if ("tier0_abstention_flag" %in% df_nms_c) suppressWarnings(as.integer(df$tier0_abstention_flag)) else rep(NA_integer_, n)
  t0_cal_tax_cv <- trimws(.chr_col(df, "tier0_calibrated_taxon", n)); t0_cal_prob_cv <- .num_col(df, "tier0_calibrated_prob", n)
  lopo_raw <- if (!is.null(prop_sum_attr) && !is.null(prop_sum_attr$lopo_cv_kappa)) suppressWarnings(as.numeric(prop_sum_attr$lopo_cv_kappa)) else NA_real_
  lopo_cv_kappa <- min(1.0, max(0.0, if (is.finite(lopo_raw)) lopo_raw else 0.60))
  known_species_pool <- unique(given_sp_cv[!is.na(given_sp_cv) & nzchar(trimws(given_sp_cv)) & given_sp_cv != "Unknown"])

  for (i in seq_len(n)) {
    given_sp <- given_sp_cv[i]; morph_sp <- morph_sp_cv[i]; cda_p <- cda_p_cv[i]; is_mo <- is_mo_cv[i]
    has_geo <- has_geo_cv[i]; edaph_sp <- edaph_sp_cv[i]; rf_p <- rf_p_cv[i]; rf_g <- rf_g_cv[i]
    ez_i <- ez_cv[i]; tez_i <- tez_cv[i]; ded_i <- ded_cv[i]; ddoy_r <- ddoy_r_cv[i]; is_anth <- is_anth_cv[i]
    morph_is_subcluster <- !is.na(morph_sp) && grepl("^Morphotype_K", morph_sp); edaph_is_ecotype <- !is.na(edaph_sp) && grepl("^Ecotype_", edaph_sp)

    p_m <- if (is_single_sp_cons || morph_is_subcluster) {
      if (is_mo) 0.72 else if (is.finite(sd2_cv[i])) min(0.85, max(0.05, stats::pchisq(sd2_cv[i], df = sdf_cv[i]) * 0.45)) else if (!is.na(morph_sp) || "gmm_cluster" %in% df_nms_c) 0.08 else NA_real_
    } else if (!is.na(morph_sp)) {
      if (morph_sp != given_sp) max(if (is.finite(cda_p)) cda_p else 0.75, if (is_mo) 0.80 else 0.55) else if (is_mo) 0.65 else if (is.finite(cda_p)) max(0.0, 1.0 - cda_p) else 0.10
    } else NA_real_
    p_e <- if (is_single_sp_cons || edaph_is_ecotype) {
      if (has_geo && (is.finite(ez_i) || is.finite(tez_i))) { z_u <- if (is.finite(ez_i)) ez_i else tez_i; min(0.92, max(0.04, stats::pchisq(z_u^2, df = 1) * 0.55)) }
      else if (edaph_is_ecotype && !is.na(eco_reg_cv[i]) && nzchar(eco_reg_cv[i])) { if (edaph_sp == eco_reg_cv[i]) 0.08 else 0.35 } else NA_real_
    } else if (!is.na(edaph_sp)) {
      if (edaph_sp != given_sp) { if (is.finite(rf_p)) max(0.55, min(0.98, rf_p)) else if (is.finite(ez_i) && ez_i > 2.0) min(0.95, 0.55 + 0.10 * (ez_i - 2.0)) else if (is.finite(ded_i) && ded_i < 0.35) 0.80 else 0.65 }
      else { if (is.finite(rf_g)) max(0.05, min(0.85, 1.0 - rf_g)) else if (is.finite(rf_p)) max(0.05, min(0.85, 1.0 - rf_p)) else if (is.finite(ez_i) && ez_i > 2.5) min(0.75, 0.25 * ez_i) else 0.10 }
    } else if (has_geo && (is.finite(ez_i) || is.finite(tez_i))) { z_u <- if (is.finite(ez_i)) ez_i else tez_i; min(0.95, max(0.05, stats::pchisq(z_u^2, df = 1) * 0.85)) } else NA_real_
    p_p <- if (is_anth && is.finite(ddoy_r)) min(1.0, abs(ddoy_r) / 30.0) else NA_real_
    p_g <- if (has_geno_row[i] && is.finite(w_geno_base) && geno_link_scale[i] > 0.0) {
      g_sp_i <- geno_sp_cv[i]; g_p_i <- geno_p_cv[i]; g_qm_i <- geno_qmax_cv[i]
      if (geno_adm_cv[i] || geno_qry_cv[i]) { if (is.finite(g_qm_i)) max(if (is_single_sp_cons) 0.35 else 0.40, min(0.92, 1.0 - g_qm_i)) else if (is_single_sp_cons) 0.72 else 0.65 }
      else if (!is.na(g_sp_i) && g_sp_i != given_sp) { if (is.finite(g_p_i)) max(0.55, min(0.99, g_p_i)) else if (is_single_sp_cons) 0.08 else 0.80 }
      else if (is_single_sp_cons) { if (is.finite(g_qm_i) && g_qm_i < 0.75) 0.55 else 0.08 }
      else { if (is.finite(g_p_i)) max(0.02, min(0.65, 1.0 - g_p_i)) else 0.08 }
    } else NA_real_
    has_t0_soft_i <- (geno_link_scale[i] == 0.0) && !is.na(t0_abstain_cv[i]) && (t0_abstain_cv[i] == 0L) &&
      is.finite(t0_cal_prob_cv[i]) && (t0_cal_prob_cv[i] >= 0.80) && !is.na(t0_cal_tax_cv[i]) && nzchar(t0_cal_tax_cv[i]) && (t0_cal_tax_cv[i] != "Sensu_Lato_DNA_Required")
    w_t0_i <- if (has_t0_soft_i && lopo_cv_kappa > 0.0) 0.50 * lopo_cv_kappa else NA_real_
    p_t0_i <- if (is.finite(w_t0_i)) {
      s_t0 <- if ((grepl(" ", t0_cal_tax_cv[i]) || (t0_cal_tax_cv[i] %in% known_species_pool)) && t0_cal_tax_cv[i] != given_sp) (1.0 - t0_cal_prob_cv[i]) else t0_cal_prob_cv[i]
      min(1.0, max(0.0, 1.0 - s_t0))
    } else NA_real_
    w_i_vec <- c(morph = w_morph_base * (0.5 + 0.5 * lookup_prec(prec_morph, morph_sp)), edaphic = w_edaph_base * (0.5 + 0.5 * lookup_prec(prec_edaph, edaph_sp)),
                 pheno = w_pheno_base, geno = if (has_geno_row[i] && is.finite(w_geno_base)) w_geno_base * geno_link_scale[i] else NA_real_, t0_soft = w_t0_i)
    p_vec <- c(morph = p_m, edaphic = p_e, pheno = p_p, geno = p_g, t0_soft = p_t0_i)
    avail <- !is.na(p_vec) & is.finite(p_vec) & !is.na(w_i_vec) & is.finite(w_i_vec); n_mod_out[i] <- sum(avail)
    p_conf_i <- if (any(avail)) sum(w_i_vec[avail] * p_vec[avail]) / sum(w_i_vec[avail]) else 0.0
    w4_vec <- pmax(0.0, w_i_vec[c("morph", "edaphic", "pheno", "geno")]); av4 <- avail[c("morph", "edaphic", "pheno", "geno")]; w4_sx <- rep(0.0, 4L)
    if (any(av4) && sum(w4_vec[av4]) > 0.0) { w4_rnd <- round(w4_vec[av4] / sum(w4_vec[av4]), 6L); im <- which.max(w4_rnd); w4_rnd[im] <- round(w4_rnd[im] + (1.0 - sum(w4_rnd)), 6L); w4_sx[av4] <- w4_rnd }
    w_morph_sx[i] <- w4_sx[1L]; w_edaph_sx[i] <- w4_sx[2L]; w_pheno_sx[i] <- w4_sx[3L]; w_geno_sx[i] <- w4_sx[4L]
    p_conf_out[i] <- round(min(1.0, max(0.0, p_conf_i)), 4); p_score_out[i] <- round(1.0 - p_conf_out[i], 4)
  }
  df$w_morph_simplex <- w_morph_sx; df$w_edaphic_simplex <- w_edaph_sx; df$w_pheno_simplex <- w_pheno_sx; df$w_geno_simplex <- w_geno_sx
  df$n_available_modalities <- n_mod_out; df$probabilistic_consensus_conflict <- p_conf_out; df$probabilistic_consensus_score <- p_score_out
  if (!is.null(prop_sum_attr)) attr(df, "genomic_propagation_summary") <- prop_sum_attr; df
}
