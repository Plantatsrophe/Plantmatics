if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# Module: scripts/analysis/step05_sdm_engine_utils.R
# Step 05 SDM Engine: Boyce Index, ROC/AUC, IPP-Cloglog MaxEnt, Binary Spatial RF, and Regional SpatRaster Stack Builder

.STEP05_SDM_WC_BIO_FALLBACK <- c(
  "1" = "bio1_temp_mean", "2" = "bio2_mean_diurnal_range", "3" = "bio3_isothermality", "4" = "bio4_temp_seasonality",
  "5" = "bio5_max_temp_warmest_month", "6" = "bio6_min_temp_coldest_month", "7" = "bio7_temp_annual_range", "8" = "bio8_mean_temp_wettest_qtr",
  "9" = "bio9_mean_temp_driest_qtr", "10" = "bio10_mean_temp_warmest_qtr", "11" = "bio11_mean_temp_coldest_qtr", "12" = "bio12_precip_annual",
  "13" = "bio13_precip_wettest_month", "14" = "bio14_precip_driest_month", "15" = "bio15_precip_seasonality", "16" = "bio16_precip_wettest_qtr",
  "17" = "bio17_precip_driest_qtr", "18" = "bio18_precip_warmest_qtr", "19" = "bio19_precip_coldest_qtr"
)

.boyce_pe_bins <- function(s_all, s_pres, lo_v, hi_v) {
  n_p <- length(s_pres); n_a <- length(s_all)
  vapply(seq_along(lo_v), function(i) {
    e_i <- sum(s_all >= lo_v[i] & s_all <= hi_v[i]) / n_a
    if (e_i > 0) (sum(s_pres >= lo_v[i] & s_pres <= hi_v[i]) / n_p) / e_i else NA_real_
  }, numeric(1))
}

compute_continuous_boyce_index <- function(suitability_all, suitability_pres, n_bins = 20L) {
  s_all <- suppressWarnings(as.numeric(suitability_all)); s_all <- s_all[is.finite(s_all)]
  s_pres <- suppressWarnings(as.numeric(suitability_pres)); s_pres <- s_pres[is.finite(s_pres)]
  empty_c <- data.frame(suitability = numeric(0), pe_ratio = numeric(0), stringsAsFactors = FALSE)
  if (length(s_pres) < 2L || length(s_all) < 4L) { res <- NA_real_; attr(res, "boyce_curve") <- empty_c; return(res) }
  s_min <- min(c(s_all, s_pres), na.rm = TRUE); s_max <- max(c(s_all, s_pres), na.rm = TRUE)
  if (!is.finite(s_min) || !is.finite(s_max) || (s_max - s_min) < 1e-8) { res <- 0.0; attr(res, "boyce_curve") <- empty_c; return(res) }
  n_p <- length(s_pres); n_b <- if (n_p < 25L) min(max(6L, n_p), as.integer(n_bins)) else max(5L, as.integer(n_bins))
  win_w <- max((s_max - s_min) / if (n_p < 25L) 4.0 else 10.0, 1e-4)
  centers <- seq(s_min + 0.5 * win_w, s_max - 0.5 * win_w, length.out = n_b)
  pe_ratio <- .boyce_pe_bins(s_all, s_pres, centers - 0.5 * win_w, centers + 0.5 * win_w)
  ok <- is.finite(pe_ratio) & is.finite(centers)
  if (sum(ok) < 5L || stats::sd(pe_ratio[ok]) < 1e-10) {
    u_all <- unique(round(s_all, 6))
    if (length(u_all) >= 6L) {
      n_q <- min(12L, max(5L, length(u_all) %/% 4L))
      q_e <- as.numeric(stats::quantile(s_all, probs = seq(0.0, 1.0, length.out = n_q + 2L), na.rm = TRUE))
      q_lo <- q_e[seq_len(n_q)]; q_hi <- q_e[seq_len(n_q) + 2L]; q_ctr <- 0.5 * (q_lo + q_hi)
      q_pe <- .boyce_pe_bins(s_all, s_pres, q_lo, q_hi); q_ok <- is.finite(q_pe) & is.finite(q_ctr)
      if (sum(q_ok) >= 4L && stats::sd(q_pe[q_ok]) >= 1e-10) { centers <- q_ctr; pe_ratio <- q_pe; ok <- q_ok }
    }
  }
  curve_df <- data.frame(suitability = centers[ok], pe_ratio = pe_ratio[ok], stringsAsFactors = FALSE)
  if (sum(ok) < 3L || stats::sd(pe_ratio[ok]) < 1e-10) { res <- 0.0; attr(res, "boyce_curve") <- curve_df; return(res) }
  cbi <- suppressWarnings(stats::cor(centers[ok], pe_ratio[ok], method = "spearman", use = "complete.obs"))
  res <- round(if (is.finite(cbi)) as.numeric(cbi) else 0.0, 4); attr(res, "boyce_curve") <- curve_df; res
}

.compute_roc_auc_curve <- function(y01, probs, wts = NULL) {
  y_v <- suppressWarnings(as.numeric(y01)); p_v <- suppressWarnings(as.numeric(probs))
  w_v <- if (is.null(wts) || length(wts) != length(y_v)) rep(1.0, length(y_v)) else suppressWarnings(as.numeric(wts))
  w_v[!is.finite(w_v) | w_v <= 0] <- 1.0; ok <- is.finite(y_v) & is.finite(p_v) & (y_v %in% c(0, 1))
  if (sum(y_v[ok] == 1) < 1L || sum(y_v[ok] == 0) < 1L) return(list(auc = NA_real_, fpr = c(0, 1), tpr = c(0, 1)))
  y_v <- y_v[ok]; p_v <- p_v[ok]; w_v <- w_v[ok]; ord <- order(p_v, decreasing = TRUE)
  y_s <- y_v[ord]; p_s <- p_v[ord]; w_s <- w_v[ord]; pos_w <- sum(w_s[y_s == 1]); neg_w <- sum(w_s[y_s == 0])
  if (pos_w <= 0 || neg_w <= 0) return(list(auc = 0.5, fpr = c(0, 1), tpr = c(0, 1)))
  dups <- !duplicated(c(Inf, p_s), fromLast = TRUE)
  tpr <- c(0, cumsum(ifelse(y_s == 1, w_s, 0)) / pos_w)[dups]; fpr <- c(0, cumsum(ifelse(y_s == 0, w_s, 0)) / neg_w)[dups]
  if (utils::tail(fpr, 1) < 1.0 || utils::tail(tpr, 1) < 1.0) { fpr <- c(fpr, 1.0); tpr <- c(tpr, 1.0) }
  list(auc = round(min(1.0, max(0.0, sum(0.5 * (tpr[-1] + tpr[-length(tpr)]) * diff(fpr)))), 4), fpr = fpr, tpr = tpr)
}

.extract_glmnet_betas <- function(fit_obj, s_val, col_nms) {
  co_m <- stats::coef(fit_obj, s = s_val); co_v <- stats::setNames(as.numeric(co_m), rownames(co_m))[col_nms]
  co_v[!is.finite(co_v)] <- 0.0; stats::setNames(co_v, col_nms)
}

.score_ipp_betas <- function(mm, betas, p_ipp, w_bg, n_pres, compute_aicc = TRUE) {
  eta_all <- as.numeric(mm %*% betas); eta_pres <- eta_all[p_ipp == 1L]; eta_bg <- eta_all[p_ipp == 0L]
  eta_shift <- max(eta_all[is.finite(eta_all)]); if (!is.finite(eta_shift)) eta_shift <- 0.0
  logL_ipp <- sum(log(pmax(exp(eta_pres - eta_shift) / max(1e-300, sum(w_bg * exp(eta_bg - eta_shift))), 1e-300)))
  k_nz <- as.integer(sum(abs(betas) > 1e-8))
  aicc_val <- if (isTRUE(compute_aicc) && n_pres > 5L && k_nz >= 1L && k_nz < (n_pres - 1L)) -2.0 * logL_ipp + 2.0 * k_nz + (2.0 * k_nz * (k_nz + 1.0)) / (n_pres - k_nz - 1.0) else if (isTRUE(compute_aicc)) Inf else NA_real_
  list(eta_bg = eta_bg, loglik_ipp = logL_ipp, k_nz = k_nz, aicc = aicc_val)
}

.fit_single_maxent_model <- function(p_vec, X_df, w_vec = NULL, rm_grid = c(0.5, 1.0, 2.0, 3.0, 4.0)) {
  if (requireNamespace("glmnet", quietly = TRUE)) {
    tryCatch(glmnet::glmnet.control(factory = TRUE), error = function(e) NULL)
    on.exit(tryCatch(glmnet::glmnet.control(factory = TRUE), error = function(e) NULL), add = TRUE)
  }
  p_num <- as.integer(as.numeric(p_vec) > 0.5); X_mat <- as.data.frame(X_df, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(X_mat) == 0L) stop("X_df has zero columns")
  cn_init <- colnames(X_mat)
  if (is.null(cn_init) || any(!nzchar(cn_init))) {
    cn_fill <- paste0("V", seq_len(ncol(X_mat))); if (!is.null(cn_init)) cn_fill[nzchar(cn_init)] <- cn_init[nzchar(cn_init)]
    colnames(X_mat) <- make.unique(cn_fill)
  }
  pred_cols <- colnames(X_mat); n_rows <- nrow(X_mat)
  for (cn in pred_cols) { v <- suppressWarnings(as.numeric(X_mat[[cn]])); mv <- stats::median(v[is.finite(v)], na.rm = TRUE); v[!is.finite(v)] <- if (is.finite(mv)) mv else 0.0; X_mat[[cn]] <- v }
  w_use <- if (is.null(w_vec) || length(w_vec) != n_rows) rep(1.0, n_rows) else as.numeric(w_vec)
  w_use[!is.finite(w_use) | w_use <= 0] <- 1.0; n_pres <- sum(p_num == 1L)
  has_nonunif_w <- (length(w_use[p_num == 1L]) >= 2L) && (stats::sd(w_use[p_num == 1L], na.rm = TRUE) > 1e-4)

  # IPP quadrature augmentation (Phillips et al. 2017, addsamplestobackground = TRUE)
  pres_idx <- which(p_num == 1L); bg_idx <- which(p_num == 0L)
  X_pres <- X_mat[pres_idx, , drop = FALSE]; X_bg <- X_mat[bg_idx, , drop = FALSE]
  w_pres_norm <- w_use[pres_idx] / max(mean(w_use[pres_idx], na.rm = TRUE), 1e-6)
  X_bg_aug <- if (nrow(X_pres) > 0L) {
    bg_k <- if (nrow(X_bg) > 0L) do.call(paste, c(lapply(X_bg, function(z) round(as.numeric(z), 6)), sep = "\r")) else character(0)
    pr_k <- do.call(paste, c(lapply(X_pres, function(z) round(as.numeric(z), 6)), sep = "\r"))
    miss_bg <- !(pr_k %in% bg_k) & !duplicated(pr_k); if (any(miss_bg)) rbind(X_bg, X_pres[miss_bg, , drop = FALSE]) else X_bg
  } else X_bg
  X_ipp <- rbind(X_pres, X_bg_aug); n_bg_ipp <- nrow(X_bg_aug)
  p_ipp <- c(rep(1L, nrow(X_pres)), rep(0L, n_bg_ipp)); w_ipp <- c(w_pres_norm, rep(100.0, n_bg_ipp))

  .align_nd <- function(new_df) {
    nd <- as.data.frame(new_df, check.names = FALSE, stringsAsFactors = FALSE); if (nrow(nd) == 0L) return(nd); cn <- colnames(nd)
    if (is.null(cn) || (any(!nzchar(cn)) && ncol(nd) == length(pred_cols))) colnames(nd) <- pred_cols[seq_len(ncol(nd))]
    else if (!all(pred_cols %in% cn) && ncol(nd) == length(pred_cols) && all(grepl("^V[0-9]+$", cn))) colnames(nd) <- pred_cols
    nd
  }
  vmin_ipp <- vapply(pred_cols, function(cn) min(as.numeric(X_ipp[[cn]]), na.rm = TRUE), numeric(1))
  vrng_ipp <- pmax(vapply(pred_cols, function(cn) max(as.numeric(X_ipp[[cn]]), na.rm = TRUE), numeric(1)) - vmin_ipp, 1e-6)
  med_ipp <- vapply(pred_cols, function(cn) stats::median(as.numeric(X_mat[[cn]]), na.rm = TRUE), numeric(1))
  hinge_knots <- seq(1.0 / 11.0, 10.0 / 11.0, length.out = 10L)

  .build_basis <- function(df_in, cl_mode = "LQ") {
    cl_lc <- tolower(paste(as.character(cl_mode), collapse = "")); nr <- nrow(df_in); np <- length(pred_cols)
    U <- matrix(0.0, nrow = nr, ncol = np, dimnames = list(NULL, paste0("l__", pred_cols)))
    for (j in seq_len(np)) {
      v <- if (pred_cols[j] %in% names(df_in)) suppressWarnings(as.numeric(df_in[[pred_cols[j]]])) else rep(med_ipp[j], nr)
      v[!is.finite(v)] <- med_ipp[j]; U[, j] <- pmin(1.0, pmax(0.0, (v - vmin_ipp[j]) / vrng_ipp[j]))
    }
    mats <- list(U); ftypes <- rep("l", np)
    if (grepl("q", cl_lc, fixed = TRUE)) { Uq <- U * U; colnames(Uq) <- paste0("q__", pred_cols); mats[[length(mats) + 1L]] <- Uq; ftypes <- c(ftypes, rep("q", np)) }
    if (grepl("h", cl_lc, fixed = TRUE)) {
      nk <- length(hinge_knots); H <- matrix(0.0, nrow = nr, ncol = 2L * nk * np); hn <- character(2L * nk * np); cp <- 1L
      for (j in seq_len(np)) for (ki in seq_len(nk)) {
        tk <- hinge_knots[ki]; uj <- U[, j]; cn <- pred_cols[j]
        H[, cp] <- pmax(0.0, (uj - tk) / (1.0 - tk)); hn[cp] <- sprintf("hf_%s_%02d", cn, ki); cp <- cp + 1L
        H[, cp] <- pmax(0.0, (tk - uj) / tk); hn[cp] <- sprintf("hr_%s_%02d", cn, ki); cp <- cp + 1L
      }
      colnames(H) <- hn; mats[[length(mats) + 1L]] <- H; ftypes <- c(ftypes, rep("h", 2L * nk * np))
    }
    B_out <- do.call(cbind, mats); attr(B_out, "ftypes") <- ftypes; B_out
  }

  fc_cands <- if (n_pres >= 15L) c("L", "LQ", "LQH") else if (n_pres >= 8L) c("L", "LQ") else c("L")
  rm_vals <- sort(unique(as.numeric(rm_grid[is.finite(rm_grid) & rm_grid > 0.0])), decreasing = TRUE)
  if (length(rm_vals) == 0L) rm_vals <- c(4.0, 3.0, 2.0, 1.0, 0.5)

  if (requireNamespace("glmnet", quietly = TRUE) && n_pres >= 3L && n_bg_ipp >= 4L && length(pred_cols) >= 1L) {
    cand_list <- list(); l_cache <- NULL; w_tot <- sum(w_ipp); sum_p_eff <- sum(w_ipp[p_ipp == 1L]); w_bg <- rep(1.0, sum(p_ipp == 0L))
    reg_lin <- stats::approx(c(0, 10, 30, 100), c(1.00, 1.00, 0.20, 0.05), xout = n_pres, rule = 2)$y
    reg_quad <- stats::approx(c(0, 10, 17, 30, 100), c(1.30, 1.30, 0.80, 0.50, 0.05), xout = n_pres, rule = 2)$y
    for (cl_mode in fc_cands) {
      tryCatch({
        glmnet::glmnet.control(factory = TRUE); mm_raw <- .build_basis(X_ipp, cl_mode = cl_mode); ft_raw <- attr(mm_raw, "ftypes")
        b_mu_a <- colMeans(mm_raw); b_sd_a <- apply(mm_raw, 2L, stats::sd)
        kb <- which(is.finite(b_mu_a) & is.finite(b_sd_a) & b_sd_a > 1e-8); if (length(kb) < 1L) stop("zero-var")
        b_mu <- stats::setNames(b_mu_a[kb], colnames(mm_raw)[kb]); b_sd <- stats::setNames(b_sd_a[kb], colnames(mm_raw)[kb])
        mm <- sweep(sweep(mm_raw[, kb, drop = FALSE], 2L, b_mu, "-"), 2L, b_sd, "/")
        c_mult <- ifelse(ft_raw[kb] == "l", reg_lin, ifelse(ft_raw[kb] == "q", reg_quad, 0.50))
        mm_pr <- mm[p_ipp == 1L, , drop = FALSE]; wp_n <- w_ipp[p_ipp == 1L] / max(sum(w_ipp[p_ipp == 1L]), 1e-12)
        mu_p <- as.numeric(crossprod(wp_n, mm_pr)); var_p <- as.numeric(crossprod(wp_n, sweep(mm_pr, 2L, mu_p, "-")^2))
        pf_base <- c_mult * sqrt(pmax(var_p, 0.0) + (1.0 / max(n_pres, 1L))) / sqrt(max(n_pres, 1L))
        pf_base[!is.finite(pf_base) | pf_base <= 1e-6] <- 1e-3; pf_fit <- pf_base / max(mean(pf_base), 1e-12); mm_fit <- mm
        if (ncol(mm_fit) == 1L) { mm_fit <- cbind(mm_fit, `__zero_pad__` = 0.0); pf_fit <- c(pf_fit, 1.0) }
        if (cl_mode == "L") l_cache <- list(mm = mm, mm_fit = mm_fit, pf_base = pf_base, pf_fit = pf_fit, b_mu = b_mu, b_sd = b_sd)
        lam_grid <- rm_vals * mean(pf_base) * (sum_p_eff / max(1e-12, w_tot))
        lam_seq <- sort(unique(c(exp(seq(log(max(max(lam_grid) * 20.0, 1e-10)), log(max(min(lam_grid) * 0.99, 1e-12)), length.out = 50L)), lam_grid)), decreasing = TRUE)
        fit_gl <- suppressWarnings(glmnet::glmnet(x = mm_fit, y = as.factor(p_ipp), family = "binomial", alpha = 1.0, standardize = FALSE, penalty.factor = pf_fit, lambda = lam_seq, weights = w_ipp, thresh = 1e-5, maxit = 1500L))
        if (is.null(fit_gl$beta) || ncol(fit_gl$beta) < 1L) stop("empty path")
        for (idx_rm in seq_along(rm_vals)) {
          bc <- .extract_glmnet_betas(fit_gl, lam_grid[idx_rm], colnames(mm)); sc <- .score_ipp_betas(mm, bc, p_ipp, w_bg, n_pres, compute_aicc = TRUE)
          cand_list[[length(cand_list) + 1L]] <- list(cl_mode = cl_mode, rm = rm_vals[idx_rm], betas_all = bc, b_mu = b_mu, b_sd = b_sd, k_nz = sc$k_nz, loglik_ipp = sc$loglik_ipp, aicc = sc$aicc, eta_bg = sc$eta_bg, is_ridge = FALSE)
        }
      }, error = function(e) NULL)
    }
    best <- NULL
    if (length(cand_list) > 0L) {
      aicc_all <- vapply(cand_list, `[[`, numeric(1), "aicc"); fin_i <- which(is.finite(aicc_all))
      if (n_pres > 5L && length(fin_i) > 0L) {
        ni <- fin_i[(aicc_all[fin_i] - min(aicc_all[fin_i])) < 0.5]
        best <- cand_list[[ni[order(vapply(cand_list[ni], `[[`, integer(1), "k_nz"), -vapply(cand_list[ni], `[[`, numeric(1), "rm"), aicc_all[ni])[1L]]]]
      } else {
        l_i <- which(vapply(cand_list, function(cnd) identical(cnd$cl_mode, "L"), logical(1)))
        rk_fn <- function(r) { m <- match(round(r, 4), c(2.0, 1.0, 0.5, 3.0, 4.0)); if (is.na(m)) 99L else m }
        lp <- l_i[vapply(cand_list[l_i], function(cnd) cnd$k_nz >= 1L && cnd$k_nz <= max(1L, n_pres - 2L), logical(1))]
        la <- l_i[vapply(cand_list[l_i], function(cnd) cnd$k_nz >= 1L, logical(1))]
        if (length(lp) > 0L) best <- cand_list[[lp[order(vapply(cand_list[lp], function(x) rk_fn(x$rm), integer(1)), vapply(cand_list[lp], `[[`, integer(1), "k_nz"))[1L]]]]
        else if (length(la) > 0L) best <- cand_list[[la[order(vapply(cand_list[la], `[[`, integer(1), "k_nz"), vapply(cand_list[la], function(x) rk_fn(x$rm), integer(1)))[1L]]]]
        if (!is.null(best)) best$aicc <- NA_real_
      }
    }
    if (is.null(best) && !is.null(l_cache)) {
      tryCatch({
        lam_r1 <- mean(l_cache$pf_base) * (sum_p_eff / max(1e-12, w_tot))
        lam_rs <- sort(unique(c(exp(seq(log(max(lam_r1 * 20.0, 1e-10)), log(max(lam_r1 * 0.99, 1e-12)), length.out = 25L)), lam_r1)), decreasing = TRUE)
        fit_r <- suppressWarnings(glmnet::glmnet(x = l_cache$mm_fit, y = as.factor(p_ipp), family = "binomial", alpha = 0.0, standardize = FALSE, penalty.factor = l_cache$pf_fit, lambda = lam_rs, weights = w_ipp, thresh = 1e-5, maxit = 1500L))
        if (!is.null(fit_r$beta) && ncol(fit_r$beta) >= 1L) {
          br <- .extract_glmnet_betas(fit_r, lam_r1, colnames(l_cache$mm)); sc_r <- .score_ipp_betas(l_cache$mm, br, p_ipp, w_bg, n_pres, compute_aicc = FALSE)
          best <- list(cl_mode = "L", rm = 1.0, betas_all = br, b_mu = l_cache$b_mu, b_sd = l_cache$b_sd, k_nz = sc_r$k_nz, loglik_ipp = sc_r$loglik_ipp, aicc = NA_real_, eta_bg = sc_r$eta_bg, is_ridge = TRUE)
        }
      }, error = function(e) NULL)
    }
    if (!is.null(best)) {
      bnz <- best$betas_all[is.finite(best$betas_all) & abs(best$betas_all) > 1e-8]
      bgs <- max(best$eta_bg[is.finite(best$eta_bg)]); if (!is.finite(bgs)) bgs <- 0.0
      wbe <- w_bg * exp(best$eta_bg - bgs); sbr <- max(1e-300, sum(wbe)); qb <- wbe / sbr
      ent <- -sum(qb * log(pmax(qb, 1e-16))); if (!is.finite(ent)) ent <- 0.0
      alp <- -log(sbr) - bgs; cl_s <- best$cl_mode; bmu_s <- best$b_mu; bsd_s <- best$b_sd
      eng_b <- if (has_nonunif_w) { if (isTRUE(best$is_ridge)) "weighted_maxnet_ridge" else "weighted_maxnet_glmnet" } else { if (isTRUE(best$is_ridge)) "maxnet_ipp_ridge" else "maxnet_ipp_glmnet" }
      return(list(
        engine = sprintf("%s (%s_rm%.2g/cloglog)", eng_b, cl_s, best$rm),
        model = list(cl_opt = cl_s, rm = best$rm, betas = bnz, b_mu = bmu_s, b_sd = bsd_s, entropy = ent, alpha = alp, aicc = best$aicc, loglik_ipp = best$loglik_ipp),
        selected_feature_class = as.character(cl_s), selected_reg_multiplier = as.numeric(best$rm), maxent_aicc = as.numeric(best$aicc), maxent_loglik_ipp = as.numeric(best$loglik_ipp), n_nonzero_betas = as.integer(best$k_nz),
        predict_cloglog = function(new_df) {
          nd <- .align_nd(new_df); if (nrow(nd) == 0L) return(numeric(0))
          if (length(bnz) == 0L) { p0 <- rep(1.0 - exp(-exp(ent + alp)), nrow(nd)); p0[!is.finite(p0)] <- 0.5; return(pmin(1.0 - 1e-6, pmax(1e-6, p0))) }
          Bn <- .build_basis(nd, cl_mode = cl_s); ub <- intersect(names(bnz), colnames(Bn)); if (length(ub) == 0L) return(rep(0.5, nrow(nd)))
          Bs <- sweep(sweep(Bn[, ub, drop = FALSE], 2L, bmu_s[ub], "-"), 2L, bsd_s[ub], "/")
          pr <- 1.0 - exp(-exp(pmin(10.0, pmax(-20.0, ent + as.numeric(Bs %*% bnz[ub]) + alp)))); pr[!is.finite(pr)] <- 1e-6; pmin(1.0 - 1e-6, pmax(1e-6, pr))
        }
      ))
    }
  }

  # Secondary fallback: CV-regularized IPP cloglog GLM
  cl_fb <- if (n_pres >= 8L) "LQ" else "L"; Phi <- .build_basis(X_ipp, cl_mode = cl_fb)
  if (requireNamespace("glmnet", quietly = TRUE) && ncol(Phi) >= 2L && n_pres >= 4L && n_bg_ipp >= 4L) {
    glmnet::glmnet.control(factory = TRUE); lam_fb <- (10^seq(3.0, 0.0, length.out = 40L)) * (sum(w_ipp[p_ipp == 1L]) / sum(w_ipp))
    nf_fb <- min(5L, max(3L, as.integer(floor(min(n_pres, n_bg_ipp) / 2L)))); fid <- integer(length(p_ipp))
    fid[p_ipp == 1L] <- rep_len(seq_len(nf_fb), sum(p_ipp == 1L)); fid[p_ipp == 0L] <- rep_len(seq_len(nf_fb), sum(p_ipp == 0L))
    cv_fit <- tryCatch(suppressWarnings(glmnet::cv.glmnet(x = Phi, y = as.factor(p_ipp), family = "binomial", standardize = TRUE, weights = w_ipp, alpha = 0.85, lambda = lam_fb, foldid = fid, thresh = 1e-5, maxit = 1000L)), error = function(e) NULL)
    if (!is.null(cv_fit)) {
      s_u <- if (is.finite(cv_fit$lambda.1se)) cv_fit$lambda.1se else cv_fit$lambda.min
      m_ep <- stats::median(as.numeric(stats::predict(cv_fit, newx = Phi, s = s_u, type = "link"))[p_ipp == 1L], na.rm = TRUE)
      c_ent <- -if (is.finite(m_ep)) m_ep else 0.0; co_fb <- as.numeric(stats::coef(cv_fit, s = s_u))[-1L]
      return(list(
        engine = sprintf("glmnet_maxent (%s/cloglog)", cl_fb), model = cv_fit, selected_feature_class = cl_fb, selected_reg_multiplier = 1.0, maxent_aicc = NA_real_, maxent_loglik_ipp = NA_real_,
        n_nonzero_betas = as.integer(sum(is.finite(co_fb) & abs(co_fb) > 1e-8)),
        predict_cloglog = function(new_df) {
          nd <- .align_nd(new_df); if (nrow(nd) == 0L) return(numeric(0)); Bn <- .build_basis(nd, cl_mode = cl_fb)
          eta <- tryCatch(as.numeric(stats::predict(cv_fit, newx = Bn, s = s_u, type = "link")) + c_ent, error = function(e) rep(0.0, nrow(Bn)))
          pr <- 1.0 - exp(-exp(pmin(10.0, pmax(-15.0, eta)))); pr[!is.finite(pr)] <- 1e-6; pmin(1.0 - 1e-6, pmax(1e-6, pr))
        }
      ))
    }
  }

  # Tertiary fallback: ridge-stabilized weighted cloglog GLM
  X_des <- cbind(Intercept = 1.0, Phi); w_glm <- c(w_pres_norm, rep(2.0 * max(sum(w_pres_norm), 1e-6) / max(n_bg_ipp, 1L), n_bg_ipp))
  Wd <- w_glm / max(mean(w_glm), 1e-6); Pp <- diag(c(0.0, rep(0.25, ncol(Phi))), nrow = ncol(X_des), ncol = ncol(X_des))
  bh <- as.numeric(tryCatch(solve(crossprod(X_des * sqrt(Wd)) + Pp, crossprod(X_des, Wd * ifelse(p_ipp == 1L, 0.85, -1.25))), error = function(e) rep(0.0, ncol(X_des))))
  m_ep <- stats::median(as.numeric(X_des %*% bh)[p_ipp == 1L], na.rm = TRUE); c_ent <- -if (is.finite(m_ep)) m_ep else 0.0
  list(
    engine = sprintf("weighted_glm_maxent (%s/cloglog)", cl_fb), model = list(beta = bh, c_ent = c_ent), selected_feature_class = cl_fb, selected_reg_multiplier = 1.0,
    maxent_aicc = NA_real_, maxent_loglik_ipp = NA_real_, n_nonzero_betas = as.integer(sum(is.finite(bh[-1L]) & abs(bh[-1L]) > 1e-8)),
    predict_cloglog = function(new_df) {
      nd <- .align_nd(new_df); if (nrow(nd) == 0L) return(numeric(0))
      pr <- 1.0 - exp(-exp(pmin(10.0, pmax(-15.0, as.numeric(cbind(Intercept = 1.0, .build_basis(nd, cl_mode = cl_fb)) %*% bh) + c_ent))))
      pr[!is.finite(pr)] <- 1e-6; pmin(1.0 - 1e-6, pmax(1e-6, pr))
    }
  )
}

.fit_single_srf_binary_model <- function(p_vec, X_df, w_vec = NULL, n_trees = 250L, num_threads = NULL) {
  p_num <- as.integer(as.numeric(p_vec) > 0.5); X_mat <- as.data.frame(X_df); pred_cols <- colnames(X_mat); n_rows <- nrow(X_mat)
  w_use <- if (is.null(w_vec) || length(w_vec) != n_rows) rep(1.0, n_rows) else as.numeric(w_vec); w_use[!is.finite(w_use) | w_use <= 0] <- 1.0
  w_rf <- w_use * ifelse(p_num == 1L, n_rows / (2.0 * max(sum(p_num == 1L), 1L)), n_rows / (2.0 * max(sum(p_num == 0L), 1L)))
  y_fac <- factor(ifelse(p_num == 1L, "Presence", "Background"), levels = c("Background", "Presence"))
  if (requireNamespace("ranger", quietly = TRUE) && sum(p_num == 1L) >= 2L && sum(p_num == 0L) >= 2L) {
    rf_th <- if (!is.null(num_threads) && length(num_threads) == 1L && is.finite(num_threads)) max(1L, as.integer(num_threads)) else if (identical(Sys.getenv("PLANTMATICS_IN_FORK_WORKER", unset = ""), "1")) 1L else if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, max(2L, as.integer(ceiling(n_rows / 16L)))) else 1L
    fit_rg <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = y_fac, X_mat, check.names = FALSE), num.trees = n_trees, probability = TRUE, case.weights = w_rf, num.threads = rf_th, seed = 42), error = function(e) NULL)
    if (!is.null(fit_rg)) return(list(model = fit_rg, predict_prob = function(new_df) {
      nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]; if (nrow(nd) == 0L) return(numeric(0))
      for (cn in pred_cols) { v <- suppressWarnings(as.numeric(nd[[cn]])); v[!is.finite(v)] <- stats::median(X_mat[[cn]], na.rm = TRUE); nd[[cn]] <- v }
      p_th <- if (identical(Sys.getenv("PLANTMATICS_IN_FORK_WORKER", unset = ""), "1")) 1L else rf_th
      pm <- tryCatch(stats::predict(fit_rg, data = nd, num.threads = p_th)$predictions, error = function(e) NULL)
      pr <- if (is.matrix(pm) && "Presence" %in% colnames(pm)) as.numeric(pm[, "Presence"]) else rep(0.5, nrow(nd)); pr[!is.finite(pr)] <- 0.0; pmin(1.0, pmax(0.0, pr))
    }))
  }
  if (requireNamespace("randomForest", quietly = TRUE) && sum(p_num == 1L) >= 2L && sum(p_num == 0L) >= 2L) {
    fit_rf <- tryCatch(randomForest::randomForest(x = X_mat, y = y_fac, ntree = n_trees), error = function(e) NULL)
    if (!is.null(fit_rf)) return(list(model = fit_rf, predict_prob = function(new_df) {
      pm <- tryCatch(stats::predict(fit_rf, newdata = as.data.frame(new_df)[, pred_cols, drop = FALSE], type = "prob"), error = function(e) NULL)
      pr <- if (is.matrix(pm) && "Presence" %in% colnames(pm)) as.numeric(pm[, "Presence"]) else rep(0.5, nrow(new_df)); pr[!is.finite(pr)] <- 0.0; pmin(1.0, pmax(0.0, pr))
    }))
  }
  mx_fb <- .fit_single_maxent_model(p_vec, X_df, w_vec); list(model = mx_fb$model, predict_prob = mx_fb$predict_cloglog)
}

.build_regional_sdm_raster_stack <- function(df, env_dir, kept_vars, ncols = 110L, nrows = 90L, demo = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE)) return(NULL)
  lon_v <- suppressWarnings(as.numeric(df$longitude)); lat_v <- suppressWarnings(as.numeric(df$latitude))
  ok_ll <- is.finite(lon_v) & is.finite(lat_v); if (!any(ok_ll)) return(NULL)
  lon_rng <- if (sum(ok_ll) >= 20L) as.numeric(stats::quantile(lon_v[ok_ll], probs = c(0.005, 0.995), na.rm = TRUE)) else range(lon_v[ok_ll], na.rm = TRUE)
  lat_rng <- if (sum(ok_ll) >= 20L) as.numeric(stats::quantile(lat_v[ok_ll], probs = c(0.005, 0.995), na.rm = TRUE)) else range(lat_v[ok_ll], na.rm = TRUE)
  pad_x <- max(0.75, 0.10 * diff(lon_rng)); pad_y <- max(0.75, 0.10 * diff(lat_rng)); med_lat <- stats::median(lat_v[ok_ll], na.rm = TRUE)
  xmin <- max(-130.0, lon_rng[1] - pad_x); xmax <- min(-50.0, lon_rng[2] + pad_x)
  ymin <- max(if (is.finite(med_lat) && med_lat >= 28.0) 24.5 else 14.0, lat_rng[1] - pad_y); ymax <- min(66.0, lat_rng[2] + pad_y)
  if (xmax - xmin < 1.0) { xmin <- xmin - 0.75; xmax <- xmax + 0.75 }; if (ymax - ymin < 1.0) { ymin <- ymin - 0.75; ymax <- ymax + 0.75 }

  study_ext <- terra::ext(xmin, xmax, ymin, ymax); tmpl_rast <- terra::rast(ext = study_ext, ncols = ncols, nrows = nrows, crs = "EPSG:4326")
  wc_lkp <- if (exists("WC_BIO_LOOKUP", inherits = TRUE)) get("WC_BIO_LOOKUP", inherits = TRUE) else .STEP05_SDM_WC_BIO_FALLBACK

  .crop_resample <- function(r_src, nm_out, scale_div = 1.0) {
    tryCatch({
      rs <- if (!isTRUE(tryCatch(terra::same.crs(r_src, tmpl_rast), error = function(e) FALSE))) terra::project(terra::crop(r_src, terra::project(terra::as.polygons(study_ext, crs = "EPSG:4326"), terra::crs(r_src))), tmpl_rast, method = "bilinear") else terra::resample(terra::crop(r_src, study_ext), tmpl_rast, method = "bilinear")
      if (scale_div != 1.0) rs <- rs / scale_div; names(rs) <- nm_out; rs
    }, error = function(e) NULL)
  }
  .sg_div <- function(r_lyr, thr) { vm <- stats::median(terra::values(r_lyr), na.rm = TRUE); if (is.finite(vm) && vm > thr) 10.0 else 1.0 }
  .load_solus_trapezoid <- function(tif_files, prop_prefix, nm_out, scale_div = 1.0) {
    bns <- tolower(basename(tif_files))
    f0 <- tif_files[grepl(paste0("^", prop_prefix, "_0_cm"), bns)]; f5 <- tif_files[grepl(paste0("^", prop_prefix, "_5_cm"), bns)]; f15 <- tif_files[grepl(paste0("^", prop_prefix, "_15_cm"), bns)]
    if (length(f0) > 0L && length(f5) > 0L && length(f15) > 0L) {
      r0 <- tryCatch(.crop_resample(terra::rast(f0[1L])[[1L]], nm_out, scale_div), error = function(e) NULL)
      r5 <- tryCatch(.crop_resample(terra::rast(f5[1L])[[1L]], nm_out, scale_div), error = function(e) NULL)
      r15 <- tryCatch(.crop_resample(terra::rast(f15[1L])[[1L]], nm_out, scale_div), error = function(e) NULL)
      if (!is.null(r0) && !is.null(r5) && !is.null(r15)) {
        v0 <- as.numeric(terra::values(r0)[, 1L]); v5 <- as.numeric(terra::values(r5)[, 1L]); v15 <- as.numeric(terra::values(r15)[, 1L])
        v05 <- 0.5 * (v0 + v5); v515 <- 0.5 * (v5 + v15)
        v_int <- if (exists("compute_rhizosphere_0_15cm", mode = "function")) compute_rhizosphere_0_15cm(v05, v515) else ((5.0 * v05) + (10.0 * v515)) / 15.0
        r_out <- terra::rast(tmpl_rast); terra::values(r_out) <- v_int; names(r_out) <- nm_out; return(r_out)
      }
    }
    fb_f <- c(f5, f0, f15); if (length(fb_f) > 0L) tryCatch(.crop_resample(terra::rast(fb_f[1L])[[1L]], nm_out, scale_div), error = function(e) NULL) else NULL
  }

  loaded_layers <- list(); raw_aux <- list()
  if (!is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir)) {
    tif_files <- if (exists("discover_env_rasters", mode = "function")) discover_env_rasters(env_dir, recursive = TRUE, crop_global = TRUE) else {
      lf <- list.files(env_dir, pattern = "\\.(tif|tiff|grd)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE); lf[!grepl("\\.crop_tmp\\.tif$", lf, ignore.case = TRUE)]
    }
    if (length(tif_files) > 0L) {
      solus_specs <- list(
        list(p = "ph1to1h2o", v = "fine_soil_ph", aux = NULL, div = 100.0, need = "fine_soil_ph" %in% kept_vars),
        list(p = "cec7", v = "fine_soil_cec", aux = NULL, div = 10.0, need = "fine_soil_cec" %in% kept_vars),
        list(p = "fragvol", v = "fine_cfvo_pct", aux = "cfvo", div = 1.0, need = any(c("fine_cfvo_pct", "lithic_outcrop_index") %in% kept_vars)),
        list(p = "sandtotal", v = "fine_soil_sand", aux = "sand", div = 1.0, need = any(c("ilr_texture_1", "ilr_texture_2", "fine_soil_sand") %in% kept_vars)),
        list(p = "claytotal", v = "fine_soil_clay", aux = "clay", div = 1.0, need = any(c("ilr_texture_1", "ilr_texture_2", "fine_soil_clay") %in% kept_vars)),
        list(p = "silttotal", v = "soil_silt_pct", aux = "silt", div = 1.0, need = any(c("ilr_texture_1", "ilr_texture_2", "soil_silt_pct") %in% kept_vars))
      )
      for (sp in solus_specs) if (isTRUE(sp$need) && !is.null(rs <- .load_solus_trapezoid(tif_files, sp$p, sp$v, scale_div = sp$div))) {
        if (sp$v %in% kept_vars && !(sp$v %in% names(loaded_layers))) loaded_layers[[sp$v]] <- rs
        if (!is.null(sp$aux) && is.null(raw_aux[[sp$aux]])) { ra <- rs; names(ra) <- sp$aux; raw_aux[[sp$aux]] <- ra }
      }
      if (any(c("fine_bedrock_depth_cm", "lithic_outcrop_index") %in% kept_vars)) {
        f_bed <- tif_files[grepl("^anylithicdpt", tolower(basename(tif_files)))]
        if (length(f_bed) > 0L && !is.null(rb <- tryCatch(.crop_resample(terra::rast(f_bed[1L])[[1L]], "fine_bedrock_depth_cm", 1.0), error = function(e) NULL))) {
          if ("fine_bedrock_depth_cm" %in% kept_vars) loaded_layers[["fine_bedrock_depth_cm"]] <- rb; raw_aux[["bedrock"]] <- rb
        }
      }
      for (tf in tif_files) {
        bn <- tolower(basename(tf)); if (grepl("^((ph1to1h2o|cec7|fragvol|sandtotal|claytotal|silttotal)_(0|5|15)_cm|anylithicdpt)", bn)) next
        r_t <- tryCatch(terra::rast(tf), error = function(e) NULL); if (is.null(r_t)) next
        for (li in seq_along(names(r_t))) {
          lnm <- names(r_t)[li]; mapped_nm <- lnm; sdiv <- 1.0; tag <- paste(bn, tolower(lnm))
          if (!(mapped_nm %in% kept_vars)) {
            mb <- regmatches(lnm, regexec("bio[_.-]?([0-9]{1,2})$", tolower(lnm)))[[1]]
            if (length(mb) == 2L && mb[2] %in% names(wc_lkp)) mapped_nm <- wc_lkp[[mb[2]]]
            else if (grepl("phh2o.*5[-_]15", tag)) { mapped_nm <- if ("fine_soil_ph" %in% kept_vars && !("fine_soil_ph" %in% names(loaded_layers))) "fine_soil_ph" else "soil_ph"; sdiv <- .sg_div(r_t[[li]], 14.5) }
            else if (grepl("cec.*5[-_]15", tag)) { mapped_nm <- if ("fine_soil_cec" %in% kept_vars && !("fine_soil_cec" %in% names(loaded_layers))) "fine_soil_cec" else "soil_cec"; sdiv <- .sg_div(r_t[[li]], 90.0) }
            else if (grepl("cfvo.*5[-_]15", tag)) {
              mapped_nm <- if ("fine_cfvo_pct" %in% kept_vars && !("fine_cfvo_pct" %in% names(loaded_layers))) "fine_cfvo_pct" else "soil_cfvo"; sdiv <- .sg_div(r_t[[li]], 105.0)
              if (is.null(raw_aux[["cfvo"]])) raw_aux[["cfvo"]] <- .crop_resample(r_t[[li]], "cfvo", sdiv)
            } else for (txk in c("sand", "clay", "silt")) if (grepl(paste0(txk, ".*5[-_]15"), tag) && is.null(raw_aux[[txk]])) raw_aux[[txk]] <- .crop_resample(r_t[[li]], txk, .sg_div(r_t[[li]], 105.0))
          }
          if (mapped_nm %in% kept_vars && !(mapped_nm %in% names(loaded_layers)) && !is.null(r_sub <- .crop_resample(r_t[[li]], mapped_nm, scale_div = sdiv))) {
            loaded_layers[[mapped_nm]] <- r_sub
            if (mapped_nm == "fine_bedrock_depth_cm" && is.null(raw_aux[["bedrock"]])) raw_aux[["bedrock"]] <- r_sub
            if (mapped_nm == "fine_cfvo_pct" && is.null(raw_aux[["cfvo"]])) raw_aux[["cfvo"]] <- r_sub
          }
        }
      }
    }
  }

  if (any(c("ilr_texture_1", "ilr_texture_2") %in% setdiff(kept_vars, names(loaded_layers))) && !is.null(raw_aux[["sand"]]) && !is.null(raw_aux[["clay"]])) {
    sv <- terra::values(raw_aux[["sand"]])[, 1L]; cv <- terra::values(raw_aux[["clay"]])[, 1L]; siv <- if (!is.null(raw_aux[["silt"]])) terra::values(raw_aux[["silt"]])[, 1L] else NULL
    tex_g <- if (exists("compute_soil_texture_log_ratios", mode = "function")) compute_soil_texture_log_ratios(sv, cv, silt_pct = siv) else {
      ok_t <- is.finite(sv) & is.finite(cv) & sv >= 0.0 & cv >= 0.0; i1 <- rep(NA_real_, length(sv)); i2 <- i1
      if (any(ok_t)) {
        sa <- pmax(0.5, sv[ok_t]); cl <- pmax(0.5, cv[ok_t]); si <- if (!is.null(siv)) ifelse(is.finite(siv[ok_t]) & siv[ok_t] >= 0.0, pmax(0.5, siv[ok_t]), pmax(0.5, 100.0 - sa - cl)) else pmax(0.5, 100.0 - sa - cl)
        tot <- sa + cl + si; sa <- 100.0 * sa / tot; cl <- 100.0 * cl / tot; si <- 100.0 * si / tot
        i1[ok_t] <- round(sqrt(2.0 / 3.0) * log(sa / sqrt(si * cl)), 4L); i2[ok_t] <- round(sqrt(0.5) * log(si / cl), 4L)
      }
      data.frame(ilr_texture_1 = i1, ilr_texture_2 = i2)
    }
    for (tn in intersect(c("ilr_texture_1", "ilr_texture_2"), setdiff(kept_vars, names(loaded_layers)))) { r_il <- terra::rast(tmpl_rast); terra::values(r_il) <- tex_g[[tn]]; names(r_il) <- tn; loaded_layers[[tn]] <- r_il }
  }

  # Derive locked [0, 1] lithic_outcrop_index identical to voucher formula
  if ("lithic_outcrop_index" %in% kept_vars && !("lithic_outcrop_index" %in% names(loaded_layers)) && !is.null(raw_aux[["bedrock"]])) {
    bd_v <- pmin(pmax(as.numeric(terra::values(raw_aux[["bedrock"]])[, 1L]), 0.0), 200.0)
    bd_score <- pmin(1.0, pmax(0.0, 1.0 - log1p(pmax(bd_v, 0.0)) / log1p(200.0)))
    loi_v <- if (!is.null(raw_aux[["cfvo"]])) {
      cfvo_v <- pmin(pmax(as.numeric(terra::values(raw_aux[["cfvo"]])[, 1L]), 0.0), 100.0)
      ifelse(is.finite(cfvo_v), round(pmin(1.0, pmax(0.0, 0.80 * bd_score + 0.20 * (log1p(pmax(cfvo_v, 0.0)) / log1p(100.0)))), 4L), round(bd_score, 4L))
    } else round(bd_score, 4L)
    r_loi <- terra::rast(tmpl_rast); terra::values(r_loi) <- loi_v; names(r_loi) <- "lithic_outcrop_index"; loaded_layers[["lithic_outcrop_index"]] <- r_loi
  }

  if (!isTRUE(demo)) {
    if (length(loaded_layers) == 0L) { warning("Skipping Step 05 continuous SDM raster map projection: zero real GeoTIFF layers matched kept_vars in env_dir (IDW raster fabrication is disabled when demo == FALSE)."); return(NULL) }
    bio_nms <- names(loaded_layers)[grepl("^bio[0-9]+_", names(loaded_layers))]
    if (length(bio_nms) > 0L) {
      land_ref <- loaded_layers[[bio_nms[1L]]]
      for (nm_l in setdiff(names(loaded_layers), bio_nms)) {
        rl_f <- tryCatch({ rf <- terra::focal(loaded_layers[[nm_l]], w = 5, fun = mean, na.rm = TRUE, na.policy = "only"); names(rf) <- nm_l; terra::mask(rf, land_ref) }, error = function(e) NULL)
        if (!is.null(rl_f)) loaded_layers[[nm_l]] <- rl_f
      }
    }
    cov_counts <- vapply(loaded_layers, function(rl) sum(is.finite(terra::values(rl)[, 1L])), integer(1))
    loaded_layers <- loaded_layers[intersect(kept_vars, names(cov_counts)[cov_counts >= floor(0.70 * max(cov_counts, 1L))])]
    if (length(loaded_layers) < 2L) { warning(sprintf("Skipping Step 05 continuous SDM raster map projection: only %d raster layer(s) cover >=70%% of study domain in '%s' (demo == FALSE).", length(loaded_layers), env_dir)); return(NULL) }
    return(terra::rast(loaded_layers))
  }

  land_mask_rast <- if (length(loaded_layers) > 0L) loaded_layers[[1L]] else NULL; miss_vars <- setdiff(kept_vars, names(loaded_layers))
  if (length(miss_vars) > 0L) {
    u_idx <- which(ok_ll)[!duplicated(paste(round(lon_v[ok_ll], 4), round(lat_v[ok_ll], 4), sep = "_"))]
    if (length(u_idx) > 450L) u_idx <- u_idx[unique(round(seq(1L, length(u_idx), length.out = 450L)))]
    V_mat <- matrix(0.0, nrow = length(u_idx), ncol = length(miss_vars), dimnames = list(NULL, miss_vars))
    for (j in seq_along(miss_vars)) { vv <- if (miss_vars[j] %in% names(df)) suppressWarnings(as.numeric(df[[miss_vars[j]]][u_idx])) else rep(0.0, length(u_idx)); mv <- stats::median(vv[is.finite(vv)], na.rm = TRUE); vv[!is.finite(vv)] <- if (is.finite(mv)) mv else 0.0; V_mat[, j] <- vv }
    g_xy <- terra::xyFromCell(tmpl_rast, seq_len(terra::ncell(tmpl_rast)))
    W <- 1.0 / ((outer(g_xy[, 1L], lon_v[u_idx], "-")^2 + outer(g_xy[, 2L], lat_v[u_idx], "-")^2 + 0.0025)^1.15)
    G_mat <- (W / pmax(rowSums(W), 1e-12)) %*% V_mat
    for (j in seq_along(miss_vars)) {
      r_idw <- terra::rast(tmpl_rast); terra::values(r_idw) <- as.numeric(G_mat[, j])
      if (!is.null(land_mask_rast)) r_idw <- tryCatch(terra::mask(r_idw, land_mask_rast), error = function(e) r_idw)
      names(r_idw) <- miss_vars[j]; loaded_layers[[miss_vars[j]]] <- r_idw
    }
  }
  terra::rast(loaded_layers[kept_vars])
}
