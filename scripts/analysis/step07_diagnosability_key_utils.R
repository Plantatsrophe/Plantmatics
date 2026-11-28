if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# Module: scripts/analysis/step07_diagnosability_key_utils.R
# Pure POSIX ERE statistical helpers for Stage 6 Candidate Lineage Diagnosability,
# 75%-Rule Spatially Blocked Shrinkage CDA, Pairwise Niche Equivalence,
# Dichotomous Key Character Ranking, and 300-DPI PNG Visualization.

.nz_s06d <- function(x) {
  if (is.null(x)) return(character(0))
  s <- trimws(as.character(x))
  s[is.na(s) | toupper(s) %in% c("", "NA", "NAN", "NONE", "NULL")] <- ""
  s
}

.clean_regime_lbl <- function(reg_tok, is_ss = FALSE) {
  r <- gsub("^Ecotype_[0-9]+_", "", reg_tok)
  r <- gsub("(lithic_outcrop_index|topo_slope_deg)_High", "Piedmont_Upland_Outcrop", r)
  r <- gsub("(lithic_outcrop_index|topo_slope_deg)_Low", "Coastal_Plain_Sand_Matrix", r)
  r <- gsub("[^A-Za-z0-9_]+", "_", r)
  if (isTRUE(is_ss) && !grepl("^Sensu_Stricto_", r)) paste0("Sensu_Stricto_", r) else r
}

hydrate_ecotypes_from_stage4 <- function(df, s4_mc_path) {
  eco_cur <- if ("ecotype_regime" %in% names(df)) .nz_s06d(df$ecotype_regime) else rep("", nrow(df))
  if (any(grepl("^Ecotype_", eco_cur)) || !file.exists(s4_mc_path)) return(df)
  s4 <- tryCatch(utils::read.csv(s4_mc_path, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (nrow(s4) == 0L || !("ecotype_regime" %in% names(s4))) return(df)
  s4_eco <- .nz_s06d(s4$ecotype_regime)
  regs <- sort(unique(s4_eco[grepl("^Ecotype_", s4_eco)]))
  if (length(regs) < 2L) return(df)
  s4_map <- stats::setNames(s4_eco[grepl("^Ecotype_", s4_eco)], as.character(s4$catalogNumber[grepl("^Ecotype_", s4_eco)]))
  lat_v <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, nrow(df))
  lon_v <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, nrow(df))
  has_ll <- is.finite(lat_v) & is.finite(lon_v)
  new_eco <- ifelse(as.character(df$catalogNumber) %in% names(s4_map), s4_map[as.character(df$catalogNumber)], "")
  f_cols <- intersect(c("topo_slope_deg", "lithic_outcrop_index", "fine_soil_sand", "fine_soil_ph"), intersect(names(df), names(s4)))
  if (length(f_cols) >= 2L && any(has_ll & !nzchar(new_eco))) {
    M4 <- as.matrix(sapply(f_cols, function(cl) suppressWarnings(as.numeric(s4[[cl]]))))
    ok4 <- grepl("^Ecotype_", s4_eco) & is.finite(rowSums(M4))
    if (sum(ok4) >= 8L) {
      m_c <- colMeans(M4[ok4, , drop = FALSE]); s_c <- apply(M4[ok4, , drop = FALSE], 2L, stats::sd); s_c[!is.finite(s_c) | s_c < 1e-6] <- 1.0
      M4s <- scale(M4[ok4, , drop = FALSE], center = m_c, scale = s_c)
      meds <- do.call(rbind, lapply(regs, function(rg) apply(M4s[s4_eco[ok4] == rg, , drop = FALSE], 2L, stats::median)))
      rownames(meds) <- regs
      M3 <- as.matrix(sapply(f_cols, function(cl) suppressWarnings(as.numeric(df[[cl]]))))
      fill_idx <- which(has_ll & !nzchar(new_eco) & is.finite(rowSums(M3)))
      if (length(fill_idx) > 0L) {
        M3s <- scale(M3[fill_idx, , drop = FALSE], center = m_c, scale = s_c)
        d_mat <- sapply(regs, function(rg) rowSums(sweep(M3s, 2L, meds[rg, ], "-")^2))
        if (is.null(dim(d_mat))) d_mat <- matrix(d_mat, nrow = 1L)
        new_eco[fill_idx] <- regs[max.col(-d_mat, ties.method = "first")]
      }
    }
  }
  df$ecotype_regime <- ifelse(nzchar(new_eco), new_eco, df$ecotype_regime)
  # Anchor higher lithic_outcrop_index regime as sensu stricto when Forty Acre Rock outcrop regime is present
  if ("lithic_outcrop_index" %in% names(df)) {
    lv <- suppressWarnings(as.numeric(df$lithic_outcrop_index))
    med_by_r <- vapply(regs, function(rg) stats::median(lv[new_eco == rg], na.rm = TRUE), numeric(1))
    if (all(is.finite(med_by_r)) && max(med_by_r) >= 0.06 && max(med_by_r) > 1.8 * min(med_by_r)) {
      df$sensu_stricto_ecotype_regime <- names(which.max(med_by_r))[1L]
    }
  }
  df
}

assign_stage06_candidate_lineages <- function(df, min_n = 5L) {
  n <- nrow(df)
  if (n == 0L) return(list(lineage = character(0), ss_lineage = "", valid_lineages = character(0), n_lineages = 0L))
  sub_v <- if ("taxonomic_treatment_subtype" %in% names(df)) .nz_s06d(df$taxonomic_treatment_subtype) else rep("", n)
  sec_v <- if ("biogeographic_sector" %in% names(df)) .nz_s06d(df$biogeographic_sector) else rep("", n)
  eco_v <- if ("ecotype_regime" %in% names(df)) .nz_s06d(df$ecotype_regime) else rep("", n)
  ss_v <- if ("sensu_stricto_ecotype_regime" %in% names(df)) .nz_s06d(df$sensu_stricto_ecotype_regime) else rep("", n)
  gmm_v <- if ("gmm_cluster" %in% names(df)) .nz_s06d(df$gmm_cluster) else rep("", n)
  lat_v <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, n)
  has_any_disj <- any(nzchar(sec_v) & sec_v != "Continuous_Range")
  has_any_eco <- any(grepl("^Ecotype_", eco_v))

  ss_tok <- ""
  nz_ss <- ss_v[nzchar(ss_v)]
  if (length(nz_ss) > 0L) ss_tok <- names(sort(table(nz_ss), decreasing = TRUE))[1L]

  lin <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    if (sub_v[i] == "Putative_Cryptic_Species" && nzchar(gmm_v[i])) {
      lin[i] <- paste0("Cryptic_Morphotype_K", gmm_v[i])
    } else if (nzchar(sec_v[i]) && sec_v[i] != "Continuous_Range") {
      lin[i] <- sec_v[i]
    } else if (grepl("^Ecotype_", eco_v[i])) {
      lin[i] <- .clean_regime_lbl(eco_v[i], is_ss = (nzchar(ss_tok) && eco_v[i] == ss_tok))
    } else if (has_any_disj && !has_any_eco && is.finite(lat_v[i])) {
      lin[i] <- "Sensu_Stricto_Cis_Mississippi_Core"
    }
  }
  tab_l <- table(lin[!is.na(lin)])
  valid_l <- names(tab_l[tab_l >= min_n])
  if (length(valid_l) < 2L) {
    tab_g <- table(gmm_v[nzchar(gmm_v)])
    valid_g <- names(tab_g[tab_g >= min_n])
    if (length(valid_g) >= 2L) {
      lin <- ifelse(gmm_v %in% valid_g, paste0("Morphotype_K", gmm_v), NA_character_)
      tab_l <- table(lin[!is.na(lin)])
      valid_l <- names(tab_l[tab_l >= min_n])
    }
  }
  lin[!(lin %in% valid_l)] <- NA_character_
  ss_lin <- if (nzchar(ss_tok)) .clean_regime_lbl(ss_tok, is_ss = TRUE) else ""
  if (!nzchar(ss_lin) || !(ss_lin %in% valid_l)) {
    ss_lin <- if (length(valid_l) >= 1L) valid_l[1L] else ""
  }
  list(lineage = lin, ss_lineage = ss_lin, valid_lineages = sort(valid_l), n_lineages = length(valid_l))
}

select_stage06_morph_features <- function(df) {
  pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
  if (length(pc_cols) > 6L) pc_cols <- pc_cols[seq_len(6L)]
  ect_cols <- intersect(c("ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity"), names(df))
  micro_cols <- c("phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta", "involucre_basal_cie_l_delta",
                  "peduncle_cie_l_sheet_delta", "calyculus_flare_index", "cauline_cie_l_q85_sheet_delta", "basal_to_cauline_wool_delta")
  arch_cols <- intersect(c(
    "aspect_ratio", "solidity", "dissection_index", "basal_cordate_index",
    "terminal_lobe_area_frac", "margin_serration_index", "margin_tooth_acuteness",
    "base_angle_deg", "apex_angle_deg", "lamina_curvature", "petiole_to_blade_ratio",
    "blade_cie_l_q85_sheet_delta", "AR_inv", "log1p_capitulum_count", "radiate_head_fraction", micro_cols
  ), names(df))
  cand <- unique(c(pc_cols, ect_cols, arch_cols))
  keep <- character(0)
  for (cl in cand) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    fin <- v[is.finite(v)]
    min_cov <- if (cl %in% c(micro_cols, "AR_inv", "radiate_head_fraction", "log1p_capitulum_count")) 0.30 else 0.55
    if (length(fin) >= 12L && stats::sd(fin) > 1e-6 && (length(fin) / max(1L, nrow(df))) >= min_cov) {
      keep <- c(keep, cl)
    }
  }
  keep
}

.fit_shrinkage_lda_equal_prior <- function(X_tr, y_tr) {
  classes <- sort(unique(y_tr))
  G <- length(classes); p <- ncol(X_tr); N <- nrow(X_tr)
  mu_mat <- matrix(0.0, nrow = G, ncol = p, dimnames = list(classes, colnames(X_tr)))
  S_W <- matrix(0.0, nrow = p, ncol = p)
  for (g in classes) {
    idx <- which(y_tr == g)
    Xg <- X_tr[idx, , drop = FALSE]
    mu_g <- colMeans(Xg)
    mu_mat[g, ] <- mu_g
    Xc <- sweep(Xg, 2L, mu_g, "-")
    S_W <- S_W + crossprod(Xc)
  }
  S_W <- S_W / max(1.0, N - G)
  lam <- min(0.55, max(0.08, p / max(1.0, N - G)))
  d_diag <- diag(S_W)
  d_diag[!is.finite(d_diag) | d_diag < 1e-6] <- 1.0
  Sigma <- (1.0 - lam) * S_W + lam * diag(d_diag, nrow = p)
  Sigma_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  list(classes = classes, mu = mu_mat, Sigma_inv = Sigma_inv)
}

.predict_shrinkage_lda <- function(model, X_te) {
  G <- length(model$classes); n <- nrow(X_te)
  log_p <- matrix(0.0, nrow = n, ncol = G, dimnames = list(NULL, model$classes))
  for (k in seq_len(G)) {
    g <- model$classes[k]
    diff_k <- sweep(X_te, 2L, model$mu[g, ], "-")
    d2_k <- rowSums((diff_k %*% model$Sigma_inv) * diff_k)
    log_p[, k] <- -0.5 * pmax(0.0, d2_k)
  }
  mx <- apply(log_p, 1L, max)
  P <- exp(sweep(log_p, 1L, mx, "-"))
  P <- sweep(P, 1L, pmax(1e-12, rowSums(P)), "/")
  pred <- model$classes[max.col(P, ties.method = "first")]
  list(pred = pred, prob = P, max_prob = apply(P, 1L, max))
}

run_candidate_lineage_oof_cda <- function(df, lineage_vec, morph_cols, seed = 42L) {
  n <- nrow(df)
  valid_lin <- sort(unique(stats::na.omit(lineage_vec)))
  G <- length(valid_lin)
  if (G < 2L || length(morph_cols) < 2L) {
    return(list(available = FALSE, global_oof_acc = NA_real_, recall_by_lin = stats::setNames(rep(NA_real_, G), valid_lin),
                oof_pred = rep(NA_character_, n), oof_prob = rep(NA_real_, n), cv1 = rep(NA_real_, n), cv2 = rep(NA_real_, n), pair_df = data.frame()))
  }
  X_raw <- as.matrix(df[, morph_cols, drop = FALSE]); storage.mode(X_raw) <- "double"; X_unimp <- X_raw
  for (j in seq_len(ncol(X_raw))) {
    col_v <- X_raw[, j]; med_j <- stats::median(col_v[is.finite(col_v)])
    col_v[!is.finite(col_v)] <- if (is.finite(med_j)) med_j else 0.0; X_raw[, j] <- col_v
  }
  tier_v <- if ("determiner_tier" %in% names(df)) .nz_s06d(df$determiner_tier) else rep("Tier_2_Silver", n)
  is_type_v <- if ("is_type_specimen" %in% names(df)) (!is.na(df$is_type_specimen) & as.logical(df$is_type_specimen)) else rep(FALSE, n)
  is_basal_v <- if ("morph_tier_source" %in% names(df)) !(.nz_s06d(df$morph_tier_source) %in% c("Cauline_Only_Fallback", "Transitional_Passive")) else rep(TRUE, n)

  gs_mask <- !is.na(lineage_vec) & is_basal_v & (tier_v %in% c("Tier_0_Type", "Tier_1_Gold", "Tier_2_Silver") | is_type_v)
  train_mask <- gs_mask
  for (g in valid_lin) {
    if (sum(train_mask & lineage_vec == g, na.rm = TRUE) < 8L) {
      train_mask[!is.na(lineage_vec) & lineage_vec == g & is_basal_v] <- TRUE
    }
  }
  lat_v <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, n)
  lon_v <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, n)
  st_v <- if ("stateProvince" %in% names(df)) .nz_s06d(df$stateProvince) else rep("UNK", n)
  inst_v <- if ("institutionCode" %in% names(df)) .nz_s06d(df$institutionCode) else rep("INST", n)
  blk_v <- ifelse(is.finite(lat_v) & is.finite(lon_v),
                  paste0("GEO_", round(lat_v * 2.0) / 2.0, "_", round(lon_v * 2.0) / 2.0),
                  paste0("ADM_", st_v, "_", inst_v))

  tr_idx <- which(train_mask)
  set.seed(seed)
  K_folds <- 5L
  fold_assign <- integer(length(tr_idx))
  for (g in valid_lin) {
    loc_g <- which(lineage_vec[tr_idx] == g)
    blks_g <- unique(blk_v[tr_idx[loc_g]])
    if (length(blks_g) >= 2L) {
      blk_f <- stats::setNames(((seq_along(blks_g) - 1L) %% min(K_folds, length(blks_g))) + 1L, sample(blks_g))
      fold_assign[loc_g] <- blk_f[blk_v[tr_idx[loc_g]]]
    } else {
      fold_assign[loc_g] <- ((seq_along(loc_g) - 1L) %% min(K_folds, length(loc_g))) + 1L
    }
  }

  oof_pred_tr <- character(length(tr_idx))
  oof_prob_tr <- numeric(length(tr_idx))
  oof_mat_tr <- matrix(NA_real_, nrow = length(tr_idx), ncol = G, dimnames = list(NULL, valid_lin))
  for (f in sort(unique(fold_assign))) {
    in_te <- which(fold_assign == f); in_tr <- which(fold_assign != f)
    if (length(in_tr) < G + 2L || length(unique(lineage_vec[tr_idx[in_tr]])) < G) in_tr <- seq_along(tr_idx)
    X_f_tr <- X_unimp[tr_idx[in_tr], , drop = FALSE]; X_f_te <- X_unimp[tr_idx[in_te], , drop = FALSE]
    for (j in seq_len(ncol(X_f_tr))) {
      mj <- stats::median(X_f_tr[is.finite(X_f_tr[, j]), j]); if (!is.finite(mj)) mj <- 0.0
      X_f_tr[!is.finite(X_f_tr[, j]), j] <- mj; X_f_te[!is.finite(X_f_te[, j]), j] <- mj
    }
    m_f <- colMeans(X_f_tr); s_f <- apply(X_f_tr, 2L, stats::sd); s_f[!is.finite(s_f) | s_f < 1e-6] <- 1.0
    mdl_f <- .fit_shrinkage_lda_equal_prior(scale(X_f_tr, center = m_f, scale = s_f), lineage_vec[tr_idx[in_tr]])
    pr_f <- .predict_shrinkage_lda(mdl_f, scale(X_f_te, center = m_f, scale = s_f))
    oof_pred_tr[in_te] <- pr_f$pred; oof_prob_tr[in_te] <- pr_f$max_prob; oof_mat_tr[in_te, colnames(pr_f$prob)] <- pr_f$prob
  }

  repro_mc <- intersect(c("phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta", "involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta"), colnames(X_unimp))
  has_repro_tr <- if (length(repro_mc) > 0L) rowSums(is.finite(X_unimp[tr_idx, repro_mc, drop = FALSE])) >= 1L else rep(FALSE, length(tr_idx))
  repro_recall <- stats::setNames(vapply(valid_lin, function(g) {
    ig <- which(lineage_vec[tr_idx] == g & has_repro_tr)
    if (length(ig) >= 4L) mean(oof_pred_tr[ig] == g) else NA_real_
  }, numeric(1L)), valid_lin)
  repro_oof_acc <- if (all(is.finite(repro_recall))) mean(repro_recall) else NA_real_

  recall_by_lin <- stats::setNames(numeric(G), valid_lin)
  for (g in valid_lin) {
    ig <- which(lineage_vec[tr_idx] == g)
    recall_by_lin[g] <- if (length(ig) > 0L) mean(oof_pred_tr[ig] == g) else NA_real_
  }
  global_oof_acc <- mean(recall_by_lin, na.rm = TRUE)

  # Full model for canonical axes and ungeoreferenced / Bronze projection
  m_all <- colMeans(X_raw[tr_idx, , drop = FALSE])
  s_all <- apply(X_raw[tr_idx, , drop = FALSE], 2L, stats::sd); s_all[!is.finite(s_all) | s_all < 1e-6] <- 1.0
  X_all_s <- scale(X_raw, center = m_all, scale = s_all)
  mdl_full <- .fit_shrinkage_lda_equal_prior(X_all_s[tr_idx, , drop = FALSE], lineage_vec[tr_idx])
  pr_all <- .predict_shrinkage_lda(mdl_full, X_all_s)
  oof_pred_full <- pr_all$pred; oof_prob_full <- pr_all$max_prob
  oof_pred_full[tr_idx] <- oof_pred_tr; oof_prob_full[tr_idx] <- oof_prob_tr

  # Canonical variates CV1, CV2 from between-group scatter in whitened space
  B_cov <- stats::cov(mdl_full$mu)
  eig_c <- tryCatch(eigen(mdl_full$Sigma_inv %*% B_cov), error = function(e) NULL)
  cv1 <- rep(NA_real_, n); cv2 <- rep(NA_real_, n)
  if (!is.null(eig_c)) {
    V_re <- Re(eig_c$vectors)
    cv1 <- as.numeric(X_all_s %*% V_re[, 1L])
    cv2 <- if (ncol(V_re) >= 2L && G >= 3L) as.numeric(X_all_s %*% V_re[, 2L]) else as.numeric(X_all_s[, 1L])
  }

  # Pairwise Mahalanobis D, Patten-Unitt IQR separation, and pairwise OOF accuracy
  pair_rows <- list()
  pairs <- utils::combn(valid_lin, 2L, simplify = FALSE)
  for (pr in pairs) {
    ga <- pr[1L]; gb <- pr[2L]
    ia <- tr_idx[lineage_vec[tr_idx] == ga]; ib <- tr_idx[lineage_vec[tr_idx] == gb]
    diff_mu <- mdl_full$mu[ga, ] - mdl_full$mu[gb, ]
    d2_raw <- max(0.0, as.numeric(t(diff_mu) %*% mdl_full$Sigma_inv %*% diff_mu))
    p_dim <- ncol(X_all_s)
    d2_bc <- max(0.0, d2_raw - p_dim * (1.0 / max(2L, length(ia)) + 1.0 / max(2L, length(ib))))
    w_dir <- as.numeric(mdl_full$Sigma_inv %*% diff_mu)
    z_a <- as.numeric(X_all_s[ia, , drop = FALSE] %*% w_dir)
    z_b <- as.numeric(X_all_s[ib, , drop = FALSE] %*% w_dir)
    if (stats::median(z_a) < stats::median(z_b)) { tmp_z <- z_a; z_a <- z_b; z_b <- tmp_z }
    iqr_sep <- as.numeric(stats::quantile(z_a, 0.25, na.rm = TRUE)) > as.numeric(stats::quantile(z_b, 0.75, na.rm = TRUE))
    pu_d75_90 <- (mean(z_a) - stats::qt(0.75, df = max(2L, length(z_a) - 1L)) * stats::sd(z_a)) -
      (mean(z_b) + stats::qt(0.90, df = max(2L, length(z_b) - 1L)) * stats::sd(z_b))
    loc_a <- which(lineage_vec[tr_idx] == ga); loc_b <- which(lineage_vec[tr_idx] == gb)
    pa_vs_b <- mean(oof_mat_tr[loc_a, ga] >= oof_mat_tr[loc_a, gb], na.rm = TRUE)
    pb_vs_a <- mean(oof_mat_tr[loc_b, gb] > oof_mat_tr[loc_b, ga], na.rm = TRUE)
    pw_oof <- 0.5 * (pa_vs_b + pb_vs_a)
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      lineage_a = ga, lineage_b = gb, n_a = length(ia), n_b = length(ib),
      mahalanobis_d = round(sqrt(d2_bc), 3), pairwise_oof_acc = round(pw_oof, 4),
      patten_unitt_iqr_sep = isTRUE(iqr_sep), patten_unitt_d75_90 = round(pu_d75_90, 3),
      stringsAsFactors = FALSE
    )
  }
  list(available = TRUE, global_oof_acc = global_oof_acc, recall_by_lin = recall_by_lin,
       repro_oof_acc = repro_oof_acc, repro_recall_by_lin = repro_recall,
       oof_pred = oof_pred_full, oof_prob = oof_prob_full, cv1 = cv1, cv2 = cv2,
       pair_df = do.call(rbind, pair_rows))
}

.schoener_d_2d <- function(M, g_vec, ga, gb, bins = 64L) {
  ok <- is.finite(rowSums(M)) & !is.na(g_vec) & (g_vec %in% c(ga, gb))
  if (sum(ok & g_vec == ga) < 4L || sum(ok & g_vec == gb) < 4L) return(NA_real_)
  Ms <- scale(M[ok, , drop = FALSE])
  Ms[, !is.finite(colSums(Ms))] <- 0.0
  pc <- tryCatch(stats::prcomp(Ms, rank. = 2L)$x, error = function(e) Ms[, seq_len(min(2L, ncol(Ms))), drop = FALSE])
  if (ncol(pc) == 1L) pc <- cbind(pc[, 1L], stats::rnorm(nrow(pc), 0, 0.01))
  r1 <- range(pc[, 1L]); r2 <- range(pc[, 2L])
  pad1 <- max(0.25, 0.15 * diff(r1)); pad2 <- max(0.25, 0.15 * diff(r2))
  lims <- c(r1[1L] - pad1, r1[2L] + pad1, r2[1L] - pad2, r2[2L] + pad2)
  g_ok <- g_vec[ok]
  ka <- MASS::kde2d(pc[g_ok == ga, 1L], pc[g_ok == ga, 2L], n = bins, lims = lims)$z
  kb <- MASS::kde2d(pc[g_ok == gb, 1L], pc[g_ok == gb, 2L], n = bins, lims = lims)$z
  p_a <- as.numeric(ka) / max(1e-12, sum(ka)); p_b <- as.numeric(kb) / max(1e-12, sum(kb))
  round(1.0 - 0.5 * sum(abs(p_a - p_b)), 4)
}

augment_pairwise_niche_and_pheno <- function(df, lineage_vec, pair_df) {
  if (nrow(pair_df) == 0L) return(pair_df)
  clim_cols <- grep("^bio[0-9]+($|_)", names(df), value = TRUE)
  edaph_cols <- intersect(c("lithic_outcrop_index", "fine_soil_sand", "fine_soil_ph", "ilr_texture_1", "ilr_texture_2", "topo_slope_deg", "wetland_seepage_index"), names(df))
  M_clim <- if (length(clim_cols) >= 2L) as.matrix(sapply(clim_cols, function(cl) suppressWarnings(as.numeric(df[[cl]])))) else NULL
  M_edaph <- if (length(edaph_cols) >= 2L) as.matrix(sapply(edaph_cols, function(cl) suppressWarnings(as.numeric(df[[cl]])))) else NULL
  ed_ok <- if ("is_edaphic_eligible_coord" %in% names(df)) (!is.na(df$is_edaphic_eligible_coord) & as.logical(df$is_edaphic_eligible_coord)) else rep(TRUE, nrow(df))
  hop_v <- if ("hopkins_adjusted_doy" %in% names(df)) suppressWarnings(as.numeric(df$hopkins_adjusted_doy)) else rep(NA_real_, nrow(df))
  anth_v <- if ("is_active_anthesis" %in% names(df)) (!is.na(df$is_active_anthesis) & as.logical(df$is_active_anthesis)) else rep(TRUE, nrow(df))
  lat_v <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, nrow(df))
  lon_v <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, nrow(df))

  d_clim <- numeric(nrow(pair_df)); d_edaph <- numeric(nrow(pair_df))
  d_doy <- numeric(nrow(pair_df)); gap_km <- numeric(nrow(pair_df)); bar_type <- character(nrow(pair_df))
  for (i in seq_len(nrow(pair_df))) {
    ga <- pair_df$lineage_a[i]; gb <- pair_df$lineage_b[i]
    d_clim[i] <- if (!is.null(M_clim)) .schoener_d_2d(M_clim, lineage_vec, ga, gb) else NA_real_
    lin_ed <- ifelse(ed_ok, lineage_vec, NA_character_)
    d_edaph[i] <- if (!is.null(M_edaph)) .schoener_d_2d(M_edaph, lin_ed, ga, gb) else NA_real_
    ha <- hop_v[!is.na(lineage_vec) & lineage_vec == ga & anth_v & is.finite(hop_v)]
    hb <- hop_v[!is.na(lineage_vec) & lineage_vec == gb & anth_v & is.finite(hop_v)]
    if (length(ha) >= 3L && length(hb) >= 3L) {
      tha <- atan2(mean(sin(2 * pi * ha / 365.25)), mean(cos(2 * pi * ha / 365.25)))
      thb <- atan2(mean(sin(2 * pi * hb / 365.25)), mean(cos(2 * pi * hb / 365.25)))
      d_doy[i] <- round(abs(atan2(sin(tha - thb), cos(tha - thb))) * 365.25 / (2 * pi), 1)
    } else d_doy[i] <- NA_real_
    ia <- which(!is.na(lineage_vec) & lineage_vec == ga & is.finite(lat_v) & is.finite(lon_v))
    ib <- which(!is.na(lineage_vec) & lineage_vec == gb & is.finite(lat_v) & is.finite(lon_v))
    if (length(ia) >= 1L && length(ib) >= 1L) {
      rad <- pi / 180.0
      dlat <- outer(lat_v[ia] * rad, lat_v[ib] * rad, "-")
      dlon <- outer(lon_v[ia] * rad, lon_v[ib] * rad, "-")
      a_h <- sin(dlat / 2)^2 + outer(cos(lat_v[ia] * rad), cos(lat_v[ib] * rad), "*") * sin(dlon / 2)^2
      min_d <- min(6371.0 * 2.0 * asin(pmin(1.0, sqrt(pmax(0.0, a_h)))))
      gap_km[i] <- round(min_d, 1)
      bar_type[i] <- if (min_d >= 200.0) "Allopatric_Geographic_Barrier" else "Parapatric_or_Mosaic_Contact"
    } else {
      gap_km[i] <- NA_real_; bar_type[i] <- "Unresolved_Contact"
    }
  }
  pair_df$d_clim <- d_clim; pair_df$d_edaph <- d_edaph
  pair_df$delta_hopkins_doy_days <- d_doy; pair_df$min_interlineage_gap_km <- gap_km
  pair_df$contact_or_barrier_mode <- bar_type
  pair_df
}

build_dichotomous_key_table <- function(df, lineage_vec) {
  valid_lin <- sort(unique(stats::na.omit(lineage_vec)))
  if (length(valid_lin) < 2L) return(data.frame())
  is_cal <- if ("is_calibrated_mm" %in% names(df)) (!is.na(df$is_calibrated_mm) & as.logical(df$is_calibrated_mm)) else rep(FALSE, nrow(df))
  mm_traits <- c("blade_length_mm", "blade_width_mm", "petiole_length_mm", "involucre_height_mm", "involucre_width_mm", "plant_height_cm")
  ratio_traits <- c("aspect_ratio", "petiole_to_blade_ratio", "solidity", "dissection_index", "basal_cordate_index",
                    "base_angle_deg", "apex_angle_deg", "margin_serration_index", "corymb_head_count_est", "AR_inv",
                    "radiate_head_fraction", "blade_cie_l_q85_sheet_delta", "cauline_cie_l_q85_sheet_delta",
                    "basal_to_cauline_wool_delta", "phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta",
                    "involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta", "calyculus_flare_index")
  cand_t <- intersect(c(mm_traits, ratio_traits), names(df))
  rows_k <- list()
  for (tr in cand_t) {
    v <- suppressWarnings(as.numeric(df[[tr]]))
    if (tr %in% mm_traits) v[!is_cal] <- NA_real_
    ok <- !is.na(lineage_vec) & is.finite(v)
    cnt_by_g <- table(factor(lineage_vec[ok], levels = valid_lin))
    if (any(cnt_by_g < 4L)) next
    kw <- tryCatch(stats::kruskal.test(v[ok] ~ factor(lineage_vec[ok])), error = function(e) NULL)
    if (is.null(kw)) next
    N_ok <- sum(ok); G_ok <- length(valid_lin)
    eta_h <- max(0.0, (as.numeric(kw$statistic) - G_ok + 1.0) / max(1.0, N_ok - G_ok))
    q25_v <- numeric(G_ok); q75_v <- numeric(G_ok); summ_parts <- character(G_ok)
    for (k in seq_along(valid_lin)) {
      vg <- v[ok & lineage_vec == valid_lin[k]]
      qs <- as.numeric(stats::quantile(vg, c(0.25, 0.50, 0.75), na.rm = TRUE))
      q25_v[k] <- qs[1L]; q75_v[k] <- qs[3L]
      summ_parts[k] <- sprintf("%s: %.2f [%.2f to %.2f] (n=%d)", valid_lin[k], qs[2L], qs[1L], qs[3L], length(vg))
    }
    non_ov <- FALSE
    for (a in seq_len(G_ok - 1L)) for (b in (a + 1L):G_ok) {
      if (q25_v[a] > q75_v[b] || q25_v[b] > q75_v[a]) non_ov <- TRUE
    }
    rows_k[[length(rows_k) + 1L]] <- data.frame(
      trait_name = tr, is_physical_mm = (tr %in% mm_traits), n_measured = N_ok,
      kw_statistic = round(as.numeric(kw$statistic), 2), kw_p_value = as.numeric(kw$p.value),
      eta_sq = round(eta_h, 4), has_non_overlapping_iqr = non_ov,
      lineage_iqr_summary = paste(summ_parts, collapse = " | "), stringsAsFactors = FALSE
    )
  }
  if (length(rows_k) == 0L) return(data.frame())
  kdf <- do.call(rbind, rows_k)
  kdf$kw_q_value <- signif(stats::p.adjust(kdf$kw_p_value, method = "BH"), 4)
  kdf$kw_p_value <- signif(kdf$kw_p_value, 4)
  kdf[order(-as.integer(kdf$has_non_overlapping_iqr), -kdf$eta_sq), , drop = FALSE]
}

render_lineage_diagnosability_png <- function(df, lineage_vec, oof_res, key_df, focal_taxon, out_png) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !exists("save_packera_png_300", mode = "function")) return(invisible(NULL))
  valid_lin <- sort(unique(stats::na.omit(lineage_vec)))
  if (length(valid_lin) < 2L || !isTRUE(oof_res$available)) return(invisible(NULL))
  pal_cols <- if (exists("TOL_OKABE_24")) unname(TOL_OKABE_24[seq_along(valid_lin)]) else c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(valid_lin)]
  names(pal_cols) <- valid_lin
  shp_vals <- stats::setNames(c(0L, 1L, 2L, 5L, 6L)[((seq_along(valid_lin) - 1L) %% 5L) + 1L], valid_lin)

  ok_m <- !is.na(lineage_vec) & is.finite(oof_res$cv1) & is.finite(oof_res$cv2)
  df_m <- data.frame(cv1 = oof_res$cv1[ok_m], cv2 = oof_res$cv2[ok_m], lineage = factor(lineage_vec[ok_m], levels = valid_lin))
  p1 <- ggplot2::ggplot(df_m, ggplot2::aes(x = .data$cv1, y = .data$cv2, color = .data$lineage, shape = .data$lineage)) +
    ggplot2::geom_point(size = 1.9, stroke = 0.75, alpha = 0.85) +
    ggplot2::stat_ellipse(level = 0.75, linewidth = 0.7, linetype = "dashed") +
    ggplot2::scale_color_manual(values = pal_cols, name = "Candidate Lineage") +
    ggplot2::scale_shape_manual(values = shp_vals, name = "Candidate Lineage") +
    ggplot2::labs(title = sprintf("A. Supervised Morphospace (%s)", focal_taxon),
                  subtitle = sprintf("Spatially blocked OOF accuracy = %.1f%% (dashed = 75%% contour)", 100 * oof_res$global_oof_acc),
                  x = "Canonical Variate 1", y = "Canonical Variate 2") +
    theme_packera_pub(base_size = 9.5)

  sand_v <- if ("fine_soil_sand" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_sand)) else rep(NA_real_, nrow(df))
  lith_v <- if ("lithic_outcrop_index" %in% names(df)) suppressWarnings(as.numeric(df$lithic_outcrop_index)) else rep(NA_real_, nrow(df))
  ok_e <- !is.na(lineage_vec) & is.finite(sand_v) & is.finite(lith_v)
  p2 <- if (sum(ok_e) >= 8L) {
    df_e <- data.frame(sand = sand_v[ok_e], lith = lith_v[ok_e], lineage = factor(lineage_vec[ok_e], levels = valid_lin))
    ggplot2::ggplot(df_e, ggplot2::aes(x = .data$sand, y = .data$lith, color = .data$lineage, shape = .data$lineage)) +
      ggplot2::geom_point(size = 1.9, stroke = 0.75, alpha = 0.85) +
      ggplot2::scale_color_manual(values = pal_cols, guide = "none") +
      ggplot2::scale_shape_manual(values = shp_vals, guide = "none") +
      ggplot2::labs(title = "B. 30 m Substrate Partitioning", subtitle = "Rhizosphere sand (%) vs. lithic outcrop index (0 to 1)",
                    x = "Fine Soil Sand (%)", y = "Lithic Outcrop Index (0 to 1)") +
      theme_packera_pub(base_size = 9.5)
  } else p1

  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- p1 + p2 + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")
    save_packera_png_300(comb, out_png, width = 12.5, height = 5.6, dpi = 300L)
  } else {
    save_packera_png_300(p1, out_png, width = 7.5, height = 5.6, dpi = 300L)
  }
}

