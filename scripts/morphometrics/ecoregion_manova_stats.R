if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/ecoregion_manova_stats.R
# Weighted parametric MANOVA (Pillai/Wilks) and pairwise Hotelling T^2 helpers
# ==============================================================================

.resolve_default_manova_features <- function(df, feature_cols = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  cand <- if (!is.null(feature_cols) && length(feature_cols) > 0L) {
    as.character(feature_cols)
  } else {
    c(
      paste0("PC", seq_len(5)),
      paste0("ECT_PC", seq_len(3)),
      "AR_inv_std",
      "petiole_to_blade_ratio_std",
      "scape_to_basal_leaf_ratio_std",
      "aspect_ratio_std",
      "solidity_std"
    )
  }
  avail <- intersect(unique(cand), names(df))
  if (length(avail) == 0L) avail <- intersect(c(paste0("PC", seq_len(5)), paste0("ECT_PC", seq_len(3))), names(df))
  avail
}

.fit_weighted_manova_stats <- function(Y_mat, grp_vec, w_vec) {
  grp_fac <- droplevels(factor(as.character(grp_vec)))
  g_k <- nlevels(grp_fac)
  n_obs <- nrow(Y_mat)
  df_res <- n_obs - g_k
  if (g_k < 2L || df_res < 1L || ncol(Y_mat) < 1L) return(NULL)

  w_pos <- suppressWarnings(as.numeric(w_vec))
  w_pos[!is.finite(w_pos) | w_pos <= 0] <- 1.0
  w_pos <- w_pos / max(mean(w_pos), 1e-12)

  sds <- apply(Y_mat, 2, stats::sd, na.rm = TRUE)
  keep_sd <- is.finite(sds) & (sds > 1e-7)
  if (!any(keep_sd)) return(NULL)
  Y_use <- Y_mat[, keep_sd, drop = FALSE]

  p_cap <- max(1L, df_res - 1L)
  if (ncol(Y_use) > p_cap) {
    Y_use <- Y_use[, seq_len(p_cap), drop = FALSE]
  }

  lm_pre <- tryCatch(stats::lm(Y_use ~ grp_fac, weights = w_pos), error = function(e) NULL)
  if (is.null(lm_pre)) return(NULL)
  res_mat <- as.matrix(stats::residuals(lm_pre)) * sqrt(w_pos)
  qr_e <- qr(res_mat, tol = 1e-6)
  if (qr_e$rank < 1L) return(NULL)
  if (qr_e$rank < ncol(Y_use)) {
    piv <- sort(qr_e$pivot[seq_len(qr_e$rank)])
    Y_use <- Y_use[, piv, drop = FALSE]
  }

  p_dim <- ncol(Y_use)
  s_rank <- min(p_dim, g_k - 1L)

  if (p_dim >= 2L) {
    m_fit <- tryCatch(stats::manova(Y_use ~ grp_fac, weights = w_pos), error = function(e) NULL)
    if (is.null(m_fit)) return(NULL)
    s_pil <- tryCatch(summary(m_fit, test = "Pillai")$stats, error = function(e) NULL)
    s_wlk <- tryCatch(summary(m_fit, test = "Wilks")$stats, error = function(e) NULL)
    s_hl  <- tryCatch(summary(m_fit, test = "Hotelling-Lawley")$stats, error = function(e) NULL)
    if (is.null(s_pil) || is.null(s_wlk) || nrow(s_pil) < 1L) return(NULL)

    pillai_v <- pmax(0.0, as.numeric(s_pil[1, "Pillai"]))
    wilks_l  <- pmin(1.0, pmax(0.0, as.numeric(s_wlk[1, "Wilks"])))
    hl_u     <- if (!is.null(s_hl) && "Hotelling-Lawley" %in% colnames(s_hl)) {
      pmax(0.0, as.numeric(s_hl[1, "Hotelling-Lawley"]))
    } else {
      pillai_v / max(1.0 - min(pillai_v, 1.0 - 1e-12), 1e-12)
    }
    approx_f <- pmax(0.0, as.numeric(s_pil[1, "approx F"]))
    num_df   <- as.numeric(s_pil[1, "num Df"])
    den_df   <- as.numeric(s_pil[1, "den Df"])
    p_val    <- pmin(1.0, pmax(0.0, as.numeric(s_pil[1, "Pr(>F)"])))
  } else {
    aov_1 <- tryCatch(stats::anova(stats::lm(Y_use[, 1] ~ grp_fac, weights = w_pos)), error = function(e) NULL)
    if (is.null(aov_1) || nrow(aov_1) < 2L) return(NULL)
    ss_eff <- pmax(0.0, as.numeric(aov_1[1, "Sum Sq"]))
    ss_res <- pmax(1e-12, as.numeric(aov_1[2, "Sum Sq"]))
    ss_tot <- ss_eff + ss_res
    pillai_v <- ss_eff / ss_tot
    wilks_l  <- ss_res / ss_tot
    hl_u     <- ss_eff / ss_res
    approx_f <- pmax(0.0, as.numeric(aov_1[1, "F value"]))
    num_df   <- as.numeric(aov_1[1, "Df"])
    den_df   <- as.numeric(aov_1[2, "Df"])
    p_val    <- pmin(1.0, pmax(0.0, as.numeric(aov_1[1, "Pr(>F)"])))
  }

  p_eta2 <- min(1.0, max(0.0, pillai_v / max(1L, s_rank)))
  hot_t2 <- pmax(0.0, df_res * hl_u)
  mahal_d <- NA_real_
  if (g_k == 2L) {
    lvls <- levels(grp_fac)
    w1_sum <- sum(w_pos[grp_fac == lvls[1]])
    w2_sum <- sum(w_pos[grp_fac == lvls[2]])
    mahal_d <- sqrt(max(0.0, hot_t2 * (w1_sum + w2_sum) / max(w1_sum * w2_sum, 1e-12)))
  }

  list(
    n_groups = g_k,
    n_obs = n_obs,
    n_features = p_dim,
    active_features = colnames(Y_use),
    pillai_trace = pillai_v,
    wilks_lambda = wilks_l,
    hotelling_t2 = hot_t2,
    mahalanobis_d = mahal_d,
    approx_F = approx_f,
    num_df = num_df,
    den_df = den_df,
    partial_eta2 = p_eta2,
    p_value = p_val
  )
}

.eval_manova_factor_and_pairwise <- function(
  spec, df, avail_cols, n, has_positive_anchors, auth_w, is_t1,
  has_pcs, is_pas_row, is_cln_row, is_type_row, quar_mask, sort_group_levels
) {
  g_raw <- spec$vec
  min_n <- as.integer(spec$min_n)
  geo_ok_spec <- if (identical(spec$name, "regional_group")) !quar_mask else rep(TRUE, n)
  if (length(unique(g_raw[!is.na(g_raw) & geo_ok_spec])) < 2L) return(NULL)

  # Enforce strict Tier_1/Tier_2 authority anchor split when positive anchors exist
  if (has_positive_anchors) {
    cand_mask <- !is.na(g_raw) & has_pcs & is_t1 & (auth_w > 0) & geo_ok_spec
    tab_c <- table(g_raw[cand_mask])
    act_lvls <- sort_group_levels(names(tab_c[tab_c >= min_n]))
    if (length(act_lvls) < 2L) {
      cand_alt <- !is.na(g_raw) & has_pcs & (auth_w > 0) & geo_ok_spec
      tab_alt <- table(g_raw[cand_alt])
      if (sum(tab_alt >= min_n) >= 2L) {
        cand_mask <- cand_alt
        act_lvls <- sort_group_levels(names(tab_alt[tab_alt >= min_n]))
      }
    }
    w_use <- auth_w
  } else {
    cand_mask <- !is.na(g_raw) & has_pcs & is_t1 & !is_pas_row & !is_cln_row & geo_ok_spec
    tab_c <- table(g_raw[cand_mask])
    act_lvls <- sort_group_levels(names(tab_c[tab_c >= min_n]))
    if (length(act_lvls) < 2L) {
      cand_alt <- !is.na(g_raw) & has_pcs & !is_pas_row & !is_cln_row & geo_ok_spec
      tab_alt <- table(g_raw[cand_alt])
      if (sum(tab_alt >= min_n) >= 2L) {
        cand_mask <- cand_alt
        act_lvls <- sort_group_levels(names(tab_alt[tab_alt >= min_n]))
      }
    }
    w_use <- rep(1.0, n)
  }

  if (length(act_lvls) < 2L) return(NULL)
  is_train <- cand_mask & (g_raw %in% act_lvls) & geo_ok_spec
  if (has_positive_anchors && exists("compute_spatial_balance_weights", mode = "function")) w_use <- auth_w * compute_spatial_balance_weights(df, active_mask = is_train, recompute = TRUE)

  X_full <- as.matrix(df[, avail_cols, drop = FALSE])
  storage.mode(X_full) <- "double"
  for (j in seq_len(ncol(X_full))) {
    na_j <- !is.finite(X_full[, j])
    if (any(na_j)) {
      tr_ok <- is_train & !na_j
      fill_v <- if (any(tr_ok)) stats::weighted.mean(X_full[tr_ok, j], w = pmax(w_use[tr_ok], 1e-6)) else 0.0
      X_full[na_j, j] <- fill_v
    }
  }

  sel_cols <- if (exists("select_morphometric_subspace", mode = "function")) {
    suppressMessages(select_morphometric_subspace(X_full[is_train, , drop = FALSE], group_vec = g_raw[is_train], allow_supervised_f = FALSE))
  } else {
    n_min_g <- min(as.integer(table(g_raw[is_train])))
    p_max_wt <- if (n_min_g >= 6L) max(2L, floor(n_min_g / 3L)) else ncol(X_full)
    avail_cols[seq_len(min(length(avail_cols), p_max_wt))]
  }
  X_sub <- X_full[, sel_cols, drop = FALSE]
  p_sub <- ncol(X_sub)

  # Per-group 95% chi-square Mahalanobis outlier pruning on active anchors (exempting nomenclatural types)
  for (lvl in act_lvls) {
    idx_g <- which(is_train & g_raw == lvl)
    if (length(idx_g) >= 5L) {
      sub_g <- X_sub[idx_g, , drop = FALSE]
      d2_g <- tryCatch(stats::mahalanobis(sub_g, colMeans(sub_g), stats::cov(sub_g) + diag(1e-4, p_sub)), error = function(e) rep(0.0, length(idx_g)))
      drop_g <- !is.na(d2_g) & (d2_g > stats::qchisq(0.95, df = p_sub)) & !is_type_row[idx_g]
      if (sum(!drop_g) >= max(3L, min_n)) is_train[idx_g[drop_g]] <- FALSE
    }
  }

  tab_post <- table(g_raw[is_train])
  act_lvls <- sort_group_levels(names(tab_post[tab_post >= min_n]))
  if (length(act_lvls) < 2L) return(NULL)
  is_train <- is_train & (g_raw %in% act_lvls)

  act_idx <- which(is_train)
  Y_act <- X_sub[act_idx, , drop = FALSE]
  g_act <- factor(g_raw[act_idx], levels = act_lvls)
  w_act <- w_use[act_idx]

  omni <- .fit_weighted_manova_stats(Y_act, g_act, w_act)
  if (is.null(omni)) return(NULL)

  is_post_sel_gmm <- identical(spec$name, "gmm_cluster")
  grp_row <- data.frame(
    grouping_factor = spec$name,
    factor_label = spec$label,
    n_groups = as.integer(omni$n_groups),
    n_active_vouchers = as.integer(omni$n_obs),
    n_features = as.integer(omni$n_features),
    active_features = paste(omni$active_features, collapse = ";"),
    pillai_trace = round(omni$pillai_trace, 6),
    wilks_lambda = round(omni$wilks_lambda, 6),
    approx_F = if (is_post_sel_gmm) NA_real_ else round(omni$approx_F, 5),
    num_df = if (is_post_sel_gmm) NA_real_ else round(omni$num_df, 3),
    den_df = if (is_post_sel_gmm) NA_real_ else round(omni$den_df, 3),
    partial_eta2 = round(omni$partial_eta2, 6),
    p_value = if (is_post_sel_gmm) NA_real_ else signif(omni$p_value, 6),
    q_value = NA_real_,
    stringsAsFactors = FALSE
  )

  comb_mat <- utils::combn(length(act_lvls), 2L)
  n_pairs <- ncol(comb_mat)
  eval_one_pair <- function(p_idx) {
    g1 <- act_lvls[comb_mat[1L, p_idx]]
    g2 <- act_lvls[comb_mat[2L, p_idx]]
    m_pair <- (g_act == g1) | (g_act == g2)
    n1 <- sum(g_act == g1)
    n2 <- sum(g_act == g2)
    pw_stat <- .fit_weighted_manova_stats(
      Y_act[m_pair, omni$active_features, drop = FALSE],
      droplevels(g_act[m_pair]),
      w_act[m_pair]
    )
    if (is.null(pw_stat)) return(NULL)
    data.frame(
      grouping_factor = spec$name,
      group_1 = g1,
      group_2 = g2,
      contrast = paste(g1, "vs", g2),
      n_1 = as.integer(n1),
      n_2 = as.integer(n2),
      n_features = as.integer(pw_stat$n_features),
      hotelling_t2 = round(pw_stat$hotelling_t2, 5),
      mahalanobis_d = round(pw_stat$mahalanobis_d, 5),
      pillai_trace = round(pw_stat$pillai_trace, 6),
      wilks_lambda = round(pw_stat$wilks_lambda, 6),
      approx_F = if (is_post_sel_gmm) NA_real_ else round(pw_stat$approx_F, 5),
      num_df = if (is_post_sel_gmm) NA_real_ else round(pw_stat$num_df, 3),
      den_df = if (is_post_sel_gmm) NA_real_ else round(pw_stat$den_df, 3),
      partial_eta2 = round(pw_stat$partial_eta2, 6),
      p_value = if (is_post_sel_gmm) NA_real_ else signif(pw_stat$p_value, 6),
      q_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  raw_pairs <- if (n_pairs >= 40L && exists("packera_mclapply", mode = "function")) {
    n_pw_cores <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(min(32L, n_pairs), n_pairs) else 8L
    packera_mclapply(seq_len(n_pairs), eval_one_pair, mc.cores = n_pw_cores)
  } else {
    lapply(seq_len(n_pairs), eval_one_pair)
  }
  pw_fac_df <- NULL
  pair_list_factor <- Filter(Negate(is.null), raw_pairs)
  if (length(pair_list_factor) > 0L) {
    pw_fac_df <- do.call(rbind, pair_list_factor)
    fin_pw <- is.finite(pw_fac_df$p_value)
    if (any(fin_pw)) pw_fac_df$q_value[fin_pw] <- signif(stats::p.adjust(pw_fac_df$p_value[fin_pw], method = "BH"), 6)
  }
  list(group_row = grp_row, pairwise_df = pw_fac_df)
}
