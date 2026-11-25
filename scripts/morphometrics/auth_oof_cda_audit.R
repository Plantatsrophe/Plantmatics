if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_oof_cda_audit.R
# Step 03 basal OTU median reconciliation and Step 04 OOF CDA/self-ellipsoid audit helpers
# ==============================================================================

# ==============================================================================
# LEVEL 1: WITHIN-PLANT BASAL ROSETTE & CAULINE FALLBACK OTU MEDIAN RECONCILIATION
# ==============================================================================

# Reconciles Step 03 specimen-level OTU basal medians (unweighted across clean basal rosette
# leaves; cauline leaves excluded whenever >= 1 basal leaf exists) and 75th-percentile petiole/blade
# architecture traits. Short-circuits when aggregate_step03_otu_medians_fast() already populated all fields.
reconcile_step03_basal_otu_medians <- function(specimen_df, closed_leaf_df, all_harm_cols = character(0), sym_names = character(0)) {
  if (is.null(specimen_df) || is.null(closed_leaf_df) || nrow(specimen_df) == 0L || nrow(closed_leaf_df) == 0L) return(specimen_df)
  if (isTRUE(attr(specimen_df, "step03_fast_aggregated"))) return(specimen_df)
  if (!exists("aggregate_step03_otu_medians_fast", mode = "function") && file.exists("scripts/morphometrics/fourier_parallel_utils.R")) source("scripts/morphometrics/fourier_parallel_utils.R")
  orig_cols <- names(closed_leaf_df)
  for (mc in c("scientificName", "species_raw")) if (!(mc %in% orig_cols)) closed_leaf_df[[mc]] <- NA_character_
  if (!("plant_individual_id" %in% orig_cols)) closed_leaf_df$plant_individual_id <- 0L
  use_harm <- intersect(unique(c(all_harm_cols, sym_names)), orig_cols)
  fast_df <- aggregate_step03_otu_medians_fast(closed_leaf_df, use_harm, intersect(sym_names, use_harm))
  if (nrow(fast_df) == 0L) return(specimen_df)
  pid_fn <- function(d) if ("plant_individual_id" %in% names(d)) ifelse(is.na(d$plant_individual_id), 0L, as.integer(d$plant_individual_id)) else rep(0L, nrow(d))
  idx <- match(paste(as.character(specimen_df$catalogNumber), pid_fn(specimen_df), sep = "__"), paste(as.character(fast_df$catalogNumber), pid_fn(fast_df), sep = "__"))
  skip_c <- c("catalogNumber", "plant_individual_id", setdiff(c("scientificName", "species_raw", "determiner_tier", "authority_weight"), orig_cols))
  for (cl in setdiff(intersect(names(specimen_df), names(fast_df)), skip_c)) {
    vals <- fast_df[[cl]][idx]; ok_v <- !is.na(vals); if (any(ok_v)) specimen_df[[cl]][ok_v] <- vals[ok_v]
  }
  specimen_df
}

# ==============================================================================
# STEP 04 MODULAR OOF CDA & AUDIT HELPERS
# ==============================================================================
.compute_oof_cda_and_self_ellipsoid <- function(
  X_proj,
  active_idx,
  y_a,
  w_a,
  g_levels,
  group_means,
  S_reg,
  grouping_col_all,
  predicted_taxon_all,
  posterior_prob_all,
  coords_df = NULL,
  intraspecific_mode = FALSE,
  deme_radius_km = 10.0
) {
  n_all <- nrow(X_proj)
  p <- max(1L, ncol(X_proj))
  G <- length(g_levels)
  self_d2 <- rep(NA_real_, n_all)
  self_df <- rep(as.integer(p), n_all)
  oof_prob <- as.numeric(posterior_prob_all)
  cl_thr_vec <- rep(NA_real_, n_all)
  cl_map <- stats::setNames(rep(NA_real_, G), g_levels)

  S_inv <- tryCatch(
    solve(S_reg),
    error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-5), nrow = p, ncol = p)
  )

  # 1. Leave-one-out within-taxon Mahalanobis tolerance distance (Zapata & Jimenez 2012)
  #    Depends strictly on active_idx so passive Tier_3_Bronze rows have zero influence on active anchors.
  sum_w_g <- stats::setNames(numeric(G), g_levels)
  sum_wx_g <- matrix(0.0, nrow = G, ncol = p, dimnames = list(g_levels, colnames(X_proj)))
  neff_g <- stats::setNames(rep(1.0, G), g_levels)
  n_act_g <- stats::setNames(integer(G), g_levels)
  act_pos_map <- rep(NA_integer_, n_all)
  act_pos_map[active_idx] <- seq_along(active_idx)

  for (k in seq_len(G)) {
    gk <- g_levels[k]
    mk <- which(y_a == gk)
    n_act_g[k] <- length(mk)
    if (length(mk) > 0L) {
      wk <- pmax(1e-6, as.numeric(w_a[mk]))
      sum_w_g[k] <- sum(wk)
      sum_wx_g[k, ] <- colSums(X_proj[active_idx[mk], , drop = FALSE] * wk)
      neff_g[k] <- max(1.0, (sum(wk)^2) / max(sum(wk^2), 1e-12))
    } else {
      sum_wx_g[k, ] <- group_means[k, ]
      sum_w_g[k] <- 1.0
    }
  }

  grp_chr <- as.character(grouping_col_all)
  for (i in seq_len(n_all)) {
    gi <- grp_chr[i]
    if (is.na(gi) || !(gi %in% g_levels)) next
    ki <- match(gi, g_levels)
    ap <- act_pos_map[i]
    if (!is.na(ap) && n_act_g[ki] >= 3L) {
      wi <- max(1e-6, as.numeric(w_a[ap]))
      rem_w <- sum_w_g[ki] - wi
      mu_i <- if (rem_w > 1e-8) (sum_wx_g[ki, ] - wi * X_proj[i, ]) / rem_w else group_means[ki, ]
      n_eff_use <- max(1.0, neff_g[ki] - 1.0)
    } else {
      mu_i <- group_means[ki, ]
      n_eff_use <- neff_g[ki]
    }
    dv <- X_proj[i, ] - mu_i
    raw_d2 <- max(0.0, as.numeric(t(dv) %*% S_inv %*% dv))
    self_d2[i] <- round(raw_d2 / (1.0 + 1.0 / max(1.0, n_eff_use)), 4)
  }

  # 2. Deme-blocked Out-of-Fold CDA Posteriors & Confident Learning Class Thresholds (Northcutt et al. 2021)
  if (!isTRUE(intraspecific_mode) && G >= 2L && length(active_idx) >= (G + 2L) && !all(grepl("^Morphotype_K", g_levels))) {
    N_a <- length(active_idx)
    X_a <- X_proj[active_idx, , drop = FALSE]
    y_fac <- factor(y_a, levels = g_levels)
    sub_coords <- if (!is.null(coords_df) && is.data.frame(coords_df) && nrow(coords_df) >= max(active_idx)) {
      coords_df[active_idx, , drop = FALSE]
    } else {
      NULL
    }
    fold_vec <- build_taxon_blind_gathering_deme_folds(sub_coords, n_obs = N_a, deme_radius_km = deme_radius_km)
    for (k in seq_len(G)) {
      ik <- which(as.integer(y_fac) == k)
      if (length(unique(fold_vec[ik])) < 2L && length(ik) >= 2L) {
        fold_vec[ik] <- max(fold_vec, na.rm = TRUE) + seq_along(ik)
      }
    }
    u_folds <- sort(unique(fold_vec))
    P_oof_act <- matrix(1.0 / G, nrow = N_a, ncol = G, dimnames = list(NULL, g_levels))
    priors <- rep(1.0 / G, G)

    for (fid in u_folds) {
      te <- which(fold_vec == fid)
      tr <- which(fold_vec != fid)
      if (length(tr) < G) tr <- setdiff(seq_len(N_a), te[1L])
      X_tr <- X_a[tr, , drop = FALSE]
      y_tr <- y_fac[tr]
      w_tr <- pmax(1e-6, as.numeric(w_a[tr]))
      mu_f <- matrix(0.0, nrow = G, ncol = p)
      W_f <- matrix(0.0, nrow = p, ncol = p)
      Xc_f <- matrix(0.0, nrow = nrow(X_tr), ncol = p)
      w_corr <- 0.0
      for (k in seq_len(G)) {
        ik <- which(as.integer(y_tr) == k)
        if (length(ik) == 0L) {
          mu_f[k, ] <- group_means[k, ]
          next
        }
        wk <- w_tr[ik]; sw <- sum(wk)
        w_corr <- w_corr + (sum(wk^2) / sw)
        mk <- colSums(X_tr[ik, , drop = FALSE] * wk) / sw
        mu_f[k, ] <- mk
        dk <- sweep(X_tr[ik, , drop = FALSE], 2L, mk, "-")
        Xc_f[ik, ] <- dk
        W_f <- W_f + crossprod(dk, dk * wk)
      }
      S_wf <- W_f / max(sum(w_tr) - w_corr, 1.0)
      S_rf <- compute_ledoit_wolf_cov(S_wf, X_centered = Xc_f, n_obs = nrow(X_tr), weights = w_tr)
      if (is.list(S_rf) && !is.null(S_rf$S_reg)) S_rf <- S_rf$S_reg
      S_if <- tryCatch(solve(S_rf), error = function(e) S_inv)
      for (ti in seq_along(te)) {
        obs_a <- te[ti]
        diff_i <- sweep(mu_f, 2L, X_a[obs_a, ], "-")
        d2_i <- pmax(0.0, rowSums((diff_i %*% S_if) * diff_i))
        ll_i <- log(priors) - 0.5 * d2_i
        pp_i <- exp(ll_i - max(ll_i))
        P_oof_act[obs_a, ] <- pp_i / max(sum(pp_i), 1e-12)
      }
    }

    for (k in seq_len(G)) {
      ik <- which(as.integer(y_fac) == k)
      if (length(ik) > 0L) {
        wk <- pmax(1e-6, as.numeric(w_a[ik]))
        tk_raw <- sum(wk * P_oof_act[ik, k]) / sum(wk)
        cl_map[k] <- round(max(0.52, min(0.86, tk_raw)), 4)
      } else {
        cl_map[k] <- 0.68
      }
    }

    pred_chr <- as.character(predicted_taxon_all)
    for (pos_a in seq_len(N_a)) {
      ri <- active_idx[pos_a]
      pr_tx <- pred_chr[ri]
      if (!is.na(pr_tx) && (pr_tx %in% g_levels)) {
        oof_prob[ri] <- round(P_oof_act[pos_a, pr_tx], 4)
      }
    }
    hit_pr <- match(pred_chr, g_levels)
    ok_pr <- which(!is.na(hit_pr))
    if (length(ok_pr) > 0L) {
      cl_thr_vec[ok_pr] <- cl_map[hit_pr[ok_pr]]
    }
  }

  list(
    self_morph_mahal_d2 = self_d2,
    self_morph_df = self_df,
    oof_posterior_prob = round(oof_prob, 4),
    cl_target_threshold = cl_thr_vec,
    cl_class_map = cl_map
  )
}

.format_step04_anchor_audit_reason <- function(i, df, cda_res, raw_sp, pred_sp, post_p, is_act) {
  if (is.na(pred_sp) || !nzchar(pred_sp) || !is.finite(post_p)) {
    return(list(flag = FALSE, triage = "LOW", reason = "Topology_Only_Unprojected_Multi_Island"))
  }
  ats <- if (!is.null(cda_res$anchor_tier_source) && length(cda_res$anchor_tier_source) >= i) as.character(cda_res$anchor_tier_source[i]) else ""
  det_t <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier[i]) else ""
  mts <- if ("morph_tier_source" %in% names(df)) as.character(df$morph_tier_source[i]) else "Basal_Rosette"
  is_typ <- (!is.na(ats) && ats == "Tier_0_Type_Anchor") || ("is_type_specimen" %in% names(df) && isTRUE(as.logical(df$is_type_specimen[i])))
  is_promo <- (!is.na(ats) && ats == "Tier_3_Bronze_Promoted") || (!is.null(cda_res$is_promoted_bronze) && isTRUE(as.logical(cda_res$is_promoted_bronze[i])))

  # Self-tolerance ellipsoid protection (Zapata & Jimenez 2012): do not flag a specimen for interspecific
  # morphological discordance when it sits inside its recorded taxon's 95% Mahalanobis tolerance region.
  is_morphotype_grp <- grepl("^Morphotype_K", raw_sp) || grepl("^Morphotype_K", pred_sp)
  self_d2_i <- if (!is.null(cda_res$self_morph_mahal_d2) && length(cda_res$self_morph_mahal_d2) >= i) cda_res$self_morph_mahal_d2[i] else NA_real_
  self_df_i <- if (!is.null(cda_res$self_morph_df) && length(cda_res$self_morph_df) >= i && !is.na(cda_res$self_morph_df[i])) as.integer(cda_res$self_morph_df[i]) else 3L
  in_self_95 <- !is_morphotype_grp && is.finite(self_d2_i) && (self_d2_i <= stats::qchisq(0.95, df = max(1L, self_df_i)))
  oof_p_i <- if (!is.null(cda_res$oof_posterior_prob) && length(cda_res$oof_posterior_prob) >= i && is.finite(cda_res$oof_posterior_prob[i])) cda_res$oof_posterior_prob[i] else post_p
  cl_thr_i <- if (!is.null(cda_res$cl_target_threshold) && length(cda_res$cl_target_threshold) >= i && is.finite(cda_res$cl_target_threshold[i])) cda_res$cl_target_threshold[i] else NA_real_

  if (is_act) {
    act_lbl <- if (is_typ) "Type_Anchor" else if (is_promo) "Promoted_Bronze_Anchor" else if (!is.na(det_t) && det_t == "Tier_2_Silver") "Tier2_Silver_Anchor" else "Tier1_Gold_Anchor"
    disc_pref <- if (is_typ) "Type_Anchor" else if (is_promo) "Promoted_Bronze" else if (!is.na(det_t) && det_t == "Tier_2_Silver") "Tier_2_Silver" else "Tier_1_Gold"
    if (raw_sp == pred_sp) return(list(flag = FALSE, triage = "CLEAN", reason = paste0("Verified_", act_lbl)))
    if (in_self_95) return(list(flag = FALSE, triage = "LOW", reason = paste0(disc_pref, "_Within_Self_Ellipsoid")))
    act_cut <- if (is.finite(cl_thr_i)) max(0.80, cl_thr_i) else 0.85
    if (oof_p_i >= act_cut && post_p >= 0.80) return(list(flag = TRUE, triage = "HIGH", reason = sprintf("%s_Morphological_Discordance_to_%s", disc_pref, gsub(" ", "_", pred_sp))))
    return(list(flag = FALSE, triage = "MEDIUM", reason = paste0(disc_pref, "_Borderline_Variant")))
  }
  cong_lbl <- if (!is.na(mts) && mts == "Cauline_Only_Fallback") {
    "Congruent_Cauline_Fallback_Determination"
  } else if (!is.na(det_t) && det_t %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver")) {
    "Congruent_Passive_Specialist_Determination"
  } else {
    "Congruent_Bronze_Determination"
  }
  if (raw_sp == pred_sp && post_p >= 0.75) return(list(flag = FALSE, triage = "CLEAN", reason = cong_lbl))
  if (raw_sp != pred_sp && in_self_95) return(list(flag = FALSE, triage = "LOW", reason = "Within_Self_Ellipsoid_Overlap_Protected"))
  pas_cut <- if (is.finite(cl_thr_i)) max(0.65, min(0.78, cl_thr_i)) else 0.70
  if (raw_sp != pred_sp && post_p >= pas_cut) return(list(flag = TRUE, triage = "HIGH", reason = sprintf("%s_Discordant_to_%s", gsub(" ", "_", raw_sp), gsub(" ", "_", pred_sp))))
  if (post_p < 0.60) return(list(flag = FALSE, triage = "MEDIUM", reason = "Morphological_Intermediate_Ambiguous"))
  if (raw_sp != pred_sp) return(list(flag = TRUE, triage = "MEDIUM", reason = "Candidate_Discordance_Moderate_Posterior"))
  list(flag = FALSE, triage = "LOW", reason = "Minor_Uncertainty")
}
