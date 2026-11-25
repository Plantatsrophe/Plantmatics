if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_shrinkage_cda.R
# Gathering/deme-blocked shrinkage CDA cross-validation and whitened eigen-solver
# ==============================================================================

compute_weighted_shrinkage_cda_loocv <- function(
  X_a,
  y_a,
  w_a = NULL,
  coords_df = NULL,
  deme_radius_km = 10.0,
  X_cand = NULL,
  allow_supervised_f = TRUE,
  trim_fold_outliers = TRUE
) {
  # Purely morphological gathering/deme-grouped cross-validation:
  #   - Zero geographic coordinate priors (preserves Padial et al. 2010 modality independence)
  #   - Holds out entire gatherings / 10 km spatial localities together to prevent pseudoreplication leakage
  #   - Imputes non-finite candidate features strictly from training-fold medians, re-estimates
  #     the feature subspace (when X_cand is supplied), optionally trims 95% chi-sq training outliers
  #     (exempting nomenclatural types), and fits Kish-N_eff Ledoit-Wolf covariance inside each fold
  X_a <- as.matrix(X_a)
  storage.mode(X_a) <- "double"
  if (is.null(colnames(X_a)) && ncol(X_a) >= 1L) {
    colnames(X_a) <- paste0("V", seq_len(ncol(X_a)))
  }
  y_f <- droplevels(as.factor(y_a))
  g_levels <- levels(y_f)
  K <- length(g_levels)
  N_a <- nrow(X_a)
  p <- ncol(X_a)
  if (N_a < K + 2L || K < 2L || p < 1L) return(NA_real_)

  use_cand <- FALSE
  if (!is.null(X_cand)) {
    X_cand <- as.matrix(X_cand)
    storage.mode(X_cand) <- "double"
    if (is.matrix(X_cand) && nrow(X_cand) == N_a && ncol(X_cand) >= 1L) {
      if (is.null(colnames(X_cand))) colnames(X_cand) <- paste0("C", seq_len(ncol(X_cand)))
      use_cand <- TRUE
    }
  }

  if (is.null(w_a) || length(w_a) != N_a) {
    w_a <- rep(1.0, N_a)
  } else {
    w_a <- suppressWarnings(as.numeric(w_a))
    w_a[!is.finite(w_a) | w_a < 0.0] <- 0.0
    if (sum(w_a) <= 1e-12) w_a <- rep(1.0, N_a)
  }

  fold_vec <- build_taxon_blind_gathering_deme_folds(coords_df, n_obs = N_a, deme_radius_km = deme_radius_km)
  # Ensure holding out a spatial deme never wipes out an entire class in small cohorts; fall back to gathering/LOO for singleton-deme taxa
  for (k in seq_len(K)) {
    idx_k <- which(as.integer(y_f) == k)
    if (length(unique(fold_vec[idx_k])) < 2L && length(idx_k) >= 2L) {
      max_f <- max(fold_vec, na.rm = TRUE)
      fold_vec[idx_k] <- max_f + seq_along(idx_k)
    }
  }
  u_folds <- sort(unique(fold_vec))
  priors <- rep(1.0 / K, K)
  pred_cv <- character(N_a)

  for (fid in u_folds) {
    test_idx <- which(fold_vec == fid)
    train_idx <- which(fold_vec != fid)
    if (length(train_idx) < K) train_idx <- setdiff(seq_len(N_a), test_idx[1L])

    y_tr_raw <- y_a[train_idx]
    y_tr <- y_f[train_idx]
    w_tr <- w_a[train_idx]
    med_tr_map <- NULL
    is_type_tr <- if (isTRUE(trim_fold_outliers) && !is.null(coords_df) && is.data.frame(coords_df) &&
                      nrow(coords_df) == N_a && exists(".identify_type_specimen_vec", mode = "function")) {
      tv <- as.logical(.identify_type_specimen_vec(coords_df[train_idx, , drop = FALSE])$is_type_specimen)
      tv[is.na(tv)] <- FALSE
      tv
    } else {
      rep(FALSE, length(train_idx))
    }

    used_cand_sel <- FALSE
    if (use_cand) {
      X_tr_cand <- X_cand[train_idx, , drop = FALSE]
      X_te_cand <- X_cand[test_idx, , drop = FALSE]
      med_tr_cand <- numeric(ncol(X_tr_cand))
      names(med_tr_cand) <- colnames(X_tr_cand)
      for (jc in seq_len(ncol(X_tr_cand))) {
        v_tr <- X_tr_cand[, jc]
        fin_tr <- is.finite(v_tr)
        m_jc <- if (any(fin_tr)) stats::median(v_tr[fin_tr], na.rm = TRUE) else 0.0
        if (!is.finite(m_jc)) m_jc <- 0.0
        med_tr_cand[jc] <- m_jc
        if (any(!fin_tr)) X_tr_cand[!fin_tr, jc] <- m_jc
        bad_te <- !is.finite(X_te_cand[, jc])
        if (any(bad_te)) X_te_cand[bad_te, jc] <- m_jc
      }
      aux_tr <- if (!is.null(coords_df) && is.data.frame(coords_df) && nrow(coords_df) == N_a) {
        coords_df[train_idx, , drop = FALSE]
      } else {
        NULL
      }
      if (isTRUE(trim_fold_outliers)) {
        keep_pre <- .compute_fold_outlier_keep_mask(X_tr_cand, y_tr, is_type_tr)
        if (any(!keep_pre) && sum(keep_pre) >= K) {
          X_tr_cand <- X_tr_cand[keep_pre, , drop = FALSE]
          y_tr_raw <- y_tr_raw[keep_pre]
          y_tr <- y_tr[keep_pre]
          w_tr <- w_tr[keep_pre]
          is_type_tr <- is_type_tr[keep_pre]
          if (!is.null(aux_tr)) aux_tr <- aux_tr[keep_pre, , drop = FALSE]
        }
      }
      sel_cols_f <- if (exists("select_morphometric_subspace", mode = "function")) {
        tryCatch(
          suppressMessages(select_morphometric_subspace(
            X_tr_cand,
            group_vec = y_tr_raw,
            aux_df = aux_tr,
            allow_supervised_f = isTRUE(allow_supervised_f)
          )),
          error = function(e) character(0)
        )
      } else {
        character(0)
      }
      sel_cols_f <- intersect(as.character(sel_cols_f), colnames(X_cand))
      if (length(sel_cols_f) >= 1L) {
        X_tr_f <- X_tr_cand[, sel_cols_f, drop = FALSE]
        X_te_f <- X_te_cand[, sel_cols_f, drop = FALSE]
        med_tr_map <- med_tr_cand[sel_cols_f]
        used_cand_sel <- TRUE
      }
    }
    if (!used_cand_sel) {
      y_tr <- y_f[train_idx]
      w_tr <- w_a[train_idx]
      if (isTRUE(trim_fold_outliers) && !is.null(coords_df) && is.data.frame(coords_df) &&
          nrow(coords_df) == N_a && exists(".identify_type_specimen_vec", mode = "function")) {
        tv <- as.logical(.identify_type_specimen_vec(coords_df[train_idx, , drop = FALSE])$is_type_specimen)
        tv[is.na(tv)] <- FALSE
        is_type_tr <- tv
      } else {
        is_type_tr <- rep(FALSE, length(train_idx))
      }
      X_tr_f <- X_a[train_idx, , drop = FALSE]
      X_te_f <- X_a[test_idx, , drop = FALSE]
    }

    if (any(!is.finite(X_tr_f)) || any(!is.finite(X_te_f))) {
      for (jc in seq_len(ncol(X_tr_f))) {
        fin_j <- is.finite(X_tr_f[, jc])
        m_j <- if (any(fin_j)) stats::median(X_tr_f[fin_j, jc], na.rm = TRUE) else 0.0
        if (!is.finite(m_j)) m_j <- 0.0
        if (any(!fin_j)) X_tr_f[!fin_j, jc] <- m_j
        bad_te_j <- !is.finite(X_te_f[, jc])
        if (any(bad_te_j)) X_te_f[bad_te_j, jc] <- m_j
      }
    }

    p_f <- ncol(X_tr_f)

    if (isTRUE(trim_fold_outliers)) {
      keep_tr <- .compute_fold_outlier_keep_mask(X_tr_f, y_tr, is_type_tr)
      if (any(!keep_tr) && sum(keep_tr) >= K) {
        X_tr_f <- X_tr_f[keep_tr, , drop = FALSE]
        y_tr <- y_tr[keep_tr]
        w_tr <- w_tr[keep_tr]
      }
    }

    mu_k_f <- matrix(0.0, nrow = K, ncol = p_f)
    W_mat_f <- matrix(0.0, nrow = p_f, ncol = p_f)
    X_c_f <- matrix(0.0, nrow = nrow(X_tr_f), ncol = p_f)
    w_tr_eff <- w_tr
    w_corr_f <- 0.0

    for (k in seq_len(K)) {
      ik <- which(as.integer(y_tr) == k)
      if (length(ik) == 0L) {
        # Fallback only if class had a single total anchor
        ik_all <- which(as.integer(y_f) == k)
        if (use_cand && !is.null(med_tr_map) && all(colnames(X_tr_f) %in% colnames(X_cand))) {
          fb_m <- X_cand[ik_all, colnames(X_tr_f), drop = FALSE]
          for (jc in seq_len(ncol(fb_m))) {
            bj <- !is.finite(fb_m[, jc])
            if (any(bj)) fb_m[bj, jc] <- med_tr_map[jc]
          }
          mu_k_f[k, ] <- colMeans(fb_m)
        } else {
          mu_k_f[k, ] <- colMeans(X_a[ik_all, , drop = FALSE], na.rm = TRUE)
          mu_k_f[k, !is.finite(mu_k_f[k, ])] <- 0.0
        }
        next
      }
      wk <- w_tr_eff[ik]; sw_k <- sum(wk)
      if (sw_k <= 1e-12) { wk <- rep(1.0, length(ik)); sw_k <- sum(wk); w_tr_eff[ik] <- wk }
      w_corr_f <- w_corr_f + (sum(wk^2) / sw_k)
      Xk <- X_tr_f[ik, , drop = FALSE]
      mk <- colSums(Xk * wk) / sw_k
      mu_k_f[k, ] <- mk
      dk <- sweep(Xk, 2L, mk, "-")
      X_c_f[ik, ] <- dk
      W_mat_f <- W_mat_f + crossprod(dk, dk * wk)
    }

    denom_f <- max(sum(w_tr_eff) - w_corr_f, 1.0)
    S_w_f <- W_mat_f / denom_f
    S_reg_f <- compute_ledoit_wolf_cov(S_w_f, X_centered = X_c_f, n_obs = nrow(X_tr_f), weights = w_tr_eff)
    if (is.list(S_reg_f) && !is.null(S_reg_f$S_reg)) S_reg_f <- S_reg_f$S_reg
    S_inv_f <- tryCatch(solve(S_reg_f), error = function(e) diag(1.0 / pmax(diag(S_reg_f), 1e-6), nrow = p_f, ncol = p_f))

    for (ti in seq_along(test_idx)) {
      obs_i <- test_idx[ti]
      diff_i <- sweep(mu_k_f, 2L, X_te_f[ti, ], "-")
      d2_i <- rowSums((diff_i %*% S_inv_f) * diff_i)
      ll_i <- -0.5 * d2_i + log(priors)
      pred_cv[obs_i] <- g_levels[which.max(ll_i)]
    }
  }
  round(mean(pred_cv == as.character(y_f), na.rm = TRUE), 4)
}

solve_shrinkage_cda_eigen <- function(S_reg, B, max_axes = NULL) {
  S_mat <- 0.5 * (as.matrix(S_reg) + t(as.matrix(S_reg)))
  B_mat <- 0.5 * (as.matrix(B) + t(as.matrix(B)))
  p <- ncol(S_mat)
  k_keep <- if (is.null(max_axes)) p else max(1L, min(as.integer(max_axes), p))

  chol_R <- tryCatch(base::chol(S_mat), error = function(e) NULL)
  if (is.null(chol_R)) {
    jitter_s <- max(1e-6, 1e-4 * mean(diag(S_mat)))
    S_mat <- S_mat + diag(jitter_s, p)
    chol_R <- tryCatch(base::chol(S_mat), error = function(e) NULL)
  }

  if (!is.null(chol_R)) {
    # Symmetric whitened generalized eigensystem: M = R^{-T} B R^{-1}
    left_whiten <- matrix(backsolve(chol_R, B_mat, transpose = TRUE), nrow = p, ncol = p)
    M_whiten <- t(matrix(backsolve(chol_R, t(left_whiten), transpose = TRUE), nrow = p, ncol = p))
    M_sym <- 0.5 * (M_whiten + t(M_whiten))
    eig_sym <- eigen(M_sym, symmetric = TRUE)
    vals <- pmax(0.0, as.numeric(eig_sym$values)[seq_len(k_keep)])
    U_sub <- eig_sym$vectors[, seq_len(k_keep), drop = FALSE]
    V_sub <- matrix(backsolve(chol_R, U_sub), nrow = p, ncol = k_keep)
  } else {
    eig_ns <- eigen(solve(S_mat, B_mat))
    real_idx <- which(abs(Im(eig_ns$values)) < 1e-6)
    if (length(real_idx) == 0L) real_idx <- seq_len(p)
    r_vals <- Re(eig_ns$values[real_idx])
    r_vecs <- Re(eig_ns$vectors[, real_idx, drop = FALSE])
    ord <- order(r_vals, decreasing = TRUE)
    take <- ord[seq_len(min(k_keep, length(ord)))]
    vals <- pmax(0.0, r_vals[take])
    V_sub <- r_vecs[, take, drop = FALSE]
  }

  for (j in seq_len(ncol(V_sub))) {
    v <- as.numeric(V_sub[, j])
    sc <- sqrt(max(as.numeric(crossprod(v, S_mat %*% v)), 1e-12))
    v <- v / sc
    dom_idx <- which.max(abs(v))
    if (length(dom_idx) == 1L && is.finite(v[dom_idx]) && v[dom_idx] < 0.0) {
      v <- -v
    }
    V_sub[, j] <- v
  }
  rownames(V_sub) <- colnames(S_mat)
  colnames(V_sub) <- paste0("Can", seq_len(ncol(V_sub)))
  list(values = vals, vectors = V_sub)
}
