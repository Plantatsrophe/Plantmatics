if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_deme_ledoit_cov.R
# Ledoit-Wolf shrinkage covariance, spatial deme folds, and fold outlier trimming
# ==============================================================================

compute_ledoit_wolf_cov <- function(
  S_w,
  X_centered = NULL,
  n_obs = NULL,
  min_lambda = 1e-5,
  max_lambda = 0.25,
  weights = NULL,
  min_shrinkage = NULL,
  max_shrinkage = NULL
) {
  if (!is.null(min_shrinkage)) min_lambda <- min_shrinkage
  if (!is.null(max_shrinkage)) max_lambda <- max_shrinkage
  S_mat <- as.matrix(S_w)
  p <- ncol(S_mat)
  dn <- dimnames(S_mat)
  S_mat[!is.finite(S_mat)] <- 0.0
  S_mat <- 0.5 * (S_mat + t(S_mat))
  min_l <- max(0.0, as.numeric(min_lambda))
  max_l <- max(min_l, min(1.0, as.numeric(max_lambda)))

  tr_s <- sum(diag(S_mat))
  mu_s <- tr_s / max(1L, p)
  if (!is.finite(mu_s) || mu_s <= 1e-12) {
    mu_s <- 1e-4
  }
  if (p <= 1L) {
    out_1d <- matrix(max(as.numeric(S_mat[1L, 1L]), mu_s * min_l, 1e-8), 1L, 1L, dimnames = dn)
    attr(out_1d, "shrinkage_lambda") <- min_l
    return(out_1d)
  }

  target_mat <- diag(mu_s, p)
  diff_target <- S_mat - target_mat
  delta2 <- sum(diff_target * diff_target) / p

  lambda_raw <- min_l
  n_eff_used <- NA_real_
  if (is.finite(delta2) && delta2 > 1e-16) {
    if (!is.null(X_centered)) {
      Xc_in <- as.matrix(X_centered)
      valid_rows <- apply(Xc_in, 1L, function(r) all(is.finite(r)))
      Xc <- Xc_in[valid_rows, , drop = FALSE]
      n_raw <- nrow(Xc)
      w_sub <- if (!is.null(weights) && length(weights) == nrow(Xc_in)) suppressWarnings(as.numeric(weights))[valid_rows] else NULL
      use_w <- !is.null(w_sub) && all(is.finite(w_sub)) && sum(w_sub > 0.0) >= 2L
      if (use_w) {
        w_pos <- pmax(0.0, as.numeric(w_sub))
        n_eff <- (sum(w_pos)^2) / max(sum(w_pos^2), 1e-12)
        n_c <- max(2.0, n_eff)
        w_norm <- w_pos / sum(w_pos)
      } else {
        n_c <- as.numeric(n_raw)
        w_norm <- rep(1.0 / max(1L, n_raw), n_raw)
      }
      n_eff_used <- n_c
      if (n_raw >= 2L && ncol(Xc) == p) {
        emp_tr <- sum((Xc * Xc) * w_norm)
        if (is.finite(emp_tr) && emp_tr > 1e-12 && tr_s > 1e-12) {
          Xc <- Xc * sqrt(tr_s / emp_tr)
        }
        q_i <- rowSums(Xc * Xc)
        s_norm2 <- sum(S_mat * S_mat)
        cross_i <- rowSums((Xc %*% S_mat) * Xc)
        u_w <- pmax(0.0, n_c * (w_norm^2))
        u_w <- u_w / max(sum(u_w), 1e-12)
        b_bar2 <- sum(u_w * pmax(0.0, q_i * q_i - 2.0 * cross_i + s_norm2)) / (p * n_c)
        beta2 <- min(delta2, max(0.0, b_bar2))
        lambda_raw <- beta2 / delta2
      }
    } else {
      w_ok <- !is.null(weights) && length(weights) >= 2L && sum(is.finite(weights) & weights > 0.0) >= 2L
      n_eff <- if (w_ok) {
        w_pos <- pmax(0.0, as.numeric(weights[is.finite(weights)]))
        max(2.0, (sum(w_pos)^2) / max(sum(w_pos^2), 1e-12))
      } else if (!is.null(n_obs) && is.finite(as.numeric(n_obs))) {
        max(2.0, as.numeric(n_obs))
      } else {
        max(2.0, as.numeric(p) + 2.0)
      }
      n_eff_used <- n_eff
      tr_s2 <- sum(S_mat * S_mat)
      denom_oas <- (n_eff + 1.0 - 2.0 / p) * (tr_s2 - (tr_s * tr_s) / p)
      if (is.finite(denom_oas) && denom_oas > 1e-16) {
        num_oas <- (1.0 - 2.0 / p) * tr_s2 + (tr_s * tr_s)
        lambda_raw <- num_oas / denom_oas
      } else {
        lambda_raw <- min(0.50, p / max(1.0, as.numeric(n_eff) * 2.0))
      }
    }
  }

  lambda <- max(min_l, min(max_l, if (is.finite(lambda_raw)) lambda_raw else min_l))
  S_reg <- (1.0 - lambda) * S_mat + lambda * target_mat
  S_reg <- 0.5 * (S_reg + t(S_reg))
  dimnames(S_reg) <- dn
  attr(S_reg, "shrinkage_lambda") <- lambda
  attr(S_reg, "n_eff") <- n_eff_used
  S_reg
}


build_taxon_blind_gathering_deme_folds <- function(coords_df, n_obs, deme_radius_km = 10.0) {
  if (is.null(coords_df) || !is.data.frame(coords_df) || nrow(coords_df) != n_obs || n_obs < 2L) {
    return(seq_len(n_obs))
  }
  if ("cv_deme_id" %in% names(coords_df) && any(!is.na(coords_df$cv_deme_id) & nzchar(as.character(coords_df$cv_deme_id)))) {
    cv_raw <- as.character(coords_df$cv_deme_id)
    cv_raw[is.na(cv_raw) | !nzchar(cv_raw)] <- paste0("ISO_", seq_len(n_obs))[is.na(cv_raw) | !nzchar(cv_raw)]
    return(match(cv_raw, unique(cv_raw)))
  }
  cat_v <- if ("catalogNumber" %in% names(coords_df)) trimws(as.character(coords_df$catalogNumber)) else as.character(seq_len(n_obs))
  cat_v[is.na(cat_v) | !nzchar(cat_v)] <- as.character(seq_len(n_obs))[is.na(cat_v) | !nzchar(cat_v)]
  exs_v <- if ("exsiccatae_key" %in% names(coords_df)) trimws(as.character(coords_df$exsiccatae_key)) else rep("", n_obs)
  exs_v[is.na(exs_v) | tolower(exs_v) %in% c("na", "none", "unknown")] <- ""
  grp_key <- ifelse(nzchar(exs_v), paste0("EXS__", exs_v), paste0("CAT__", cat_v))
  fold_id <- match(grp_key, unique(grp_key))

  lat_v <- if ("decimalLatitude" %in% names(coords_df)) suppressWarnings(as.numeric(coords_df$decimalLatitude)) else rep(NA_real_, n_obs)
  lon_v <- if ("decimalLongitude" %in% names(coords_df)) suppressWarnings(as.numeric(coords_df$decimalLongitude)) else rep(NA_real_, n_obs)
  ok_c <- which(is.finite(lat_v) & is.finite(lon_v) & abs(lat_v) <= 90.0 & abs(lon_v) <= 180.0)
  if (length(ok_c) >= 2L && is.finite(deme_radius_km) && deme_radius_km > 0.0) {
    # Complete-linkage (diameter-bounded) taxon-blind 10 km spatial locality grouping (prevents single-linkage chaining)
    lat_r <- lat_v[ok_c] * (pi / 180.0); lon_r <- lon_v[ok_c] * (pi / 180.0)
    clusters <- list()
    for (ii in seq_along(ok_c)) {
      placed <- FALSE
      for (cl_idx in seq_along(clusters)) {
        members <- clusters[[cl_idx]]
        dlat <- lat_r[members] - lat_r[ii]; dlon <- lon_r[members] - lon_r[ii]
        hav <- sin(dlat * 0.5)^2 + cos(lat_r[ii]) * cos(lat_r[members]) * (sin(dlon * 0.5)^2)
        dkm <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hav))))
        if (all(is.finite(dkm) & dkm <= deme_radius_km)) {
          clusters[[cl_idx]] <- c(members, ii); placed <- TRUE; break
        }
      }
      if (!placed) clusters[[length(clusters) + 1L]] <- ii
    }
    for (cl_idx in seq_along(clusters)) {
      rows_cl <- ok_c[clusters[[cl_idx]]]
      if (length(rows_cl) >= 2L) {
        old_f <- unique(fold_id[rows_cl])
        fold_id[fold_id %in% old_f] <- min(old_f)
      }
    }
  }
  match(fold_id, unique(fold_id))
}

assign_spatial_gathering_demes <- function(df, deme_radius_km = 10.0) {
  # Feature-preserving, taxon-blind spatial/gathering deme metadata assignment (ZERO feature shrinkage!)
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  if (n == 0L) return(df)
  if ("PC1" %in% names(df) && !("PC1_raw_individual" %in% names(df))) df$PC1_raw_individual <- df$PC1
  if ("PC2" %in% names(df) && !("PC2_raw_individual" %in% names(df))) df$PC2_raw_individual <- df$PC2
  fold_ids <- build_taxon_blind_gathering_deme_folds(df, n_obs = n, deme_radius_km = deme_radius_km)
  df$cv_deme_id <- sprintf("LOC_DEME_%04d", as.integer(fold_ids))
  df$population_deme_id <- df$cv_deme_id
  tab_d <- table(df$population_deme_id)
  df$deme_n_specimens <- as.integer(tab_d[df$population_deme_id])
  df
}

.compute_fold_outlier_keep_mask <- function(X_mat, y_fac, is_type_vec = NULL) {
  n_r <- nrow(X_mat)
  p_m <- ncol(X_mat)
  keep <- rep(TRUE, n_r)
  if (n_r < 5L || p_m < 1L) return(keep)
  if (is.null(is_type_vec) || length(is_type_vec) != n_r) is_type_vec <- rep(FALSE, n_r)

  use_cols <- which(vapply(seq_len(p_m), function(j) {
    v <- X_mat[, j]
    fin <- is.finite(v)
    sum(fin) >= 5L && (stats::mad(v[fin], na.rm = TRUE) > 1e-8 || stats::var(v[fin], na.rm = TRUE) > 1e-10)
  }, logical(1)))
  if (length(use_cols) == 0L) return(keep)

  X_d <- X_mat[, use_cols, drop = FALSE]
  p_d <- ncol(X_d)
  K_n <- nlevels(y_fac)

  core_resid_list <- vector("list", K_n)
  for (k_p in seq_len(K_n)) {
    ikp <- which(as.integer(y_fac) == k_p)
    if (length(ikp) >= 2L) {
      Xk_p <- X_d[ikp, , drop = FALSE]
      mk_p <- if (length(ikp) >= 5L && exists(".compute_robust_taxon_medoid", mode = "function")) {
        .compute_robust_taxon_medoid(Xk_p)
      } else {
        apply(Xk_p, 2L, stats::median, na.rm = TRUE)
      }
      diff_p <- sweep(Xk_p, 2L, mk_p, "-")
      if (length(ikp) >= 5L) {
        mads_p <- apply(Xk_p, 2L, stats::mad, na.rm = TRUE)
        sds_p <- apply(Xk_p, 2L, stats::sd, na.rm = TRUE)
        sc_p <- ifelse(is.finite(mads_p) & mads_p > 1e-8, mads_p, ifelse(is.finite(sds_p) & sds_p > 1e-8, sds_p, 1.0))
        r_p <- rowSums(sweep(diff_p, 2L, sc_p, "/")^2)
        n_core_p <- max(min(length(ikp), p_d + 2L), as.integer(floor(0.75 * length(ikp))))
        diff_p <- diff_p[order(r_p, na.last = TRUE)[seq_len(min(length(ikp), n_core_p))], , drop = FALSE]
      }
      core_resid_list[[k_p]] <- diff_p
    }
  }
  core_resid_list <- Filter(Negate(is.null), core_resid_list)
  S_pool_f <- NULL
  if (length(core_resid_list) > 0L) {
    R_p <- do.call(rbind, core_resid_list)
    if (nrow(R_p) >= 3L) {
      cov_p <- tryCatch(stats::cov(R_p), error = function(e) NULL)
      if (!is.null(cov_p) && is.matrix(cov_p)) {
        cov_p[!is.finite(cov_p)] <- 0.0
        reg_p <- max(1e-4, 0.05 * max(mean(diag(cov_p), na.rm = TRUE), 1e-4))
        S_pool_f <- 0.5 * (cov_p + t(cov_p)) + diag(reg_p, nrow = p_d, ncol = p_d)
      }
    }
  }

  q95_thr <- stats::qchisq(0.95, df = p_d)
  for (k in seq_len(K_n)) {
    ik_tr <- which(as.integer(y_fac) == k)
    if (length(ik_tr) >= 5L) {
      sub_k <- X_d[ik_tr, , drop = FALSE]
      mu_k <- if (exists(".compute_robust_taxon_medoid", mode = "function")) {
        .compute_robust_taxon_medoid(sub_k)
      } else {
        apply(sub_k, 2L, stats::median, na.rm = TRUE)
      }
      mads_k <- apply(sub_k, 2L, stats::mad, na.rm = TRUE)
      sds_k <- apply(sub_k, 2L, stats::sd, na.rm = TRUE)
      sc_k <- ifelse(is.finite(mads_k) & mads_k > 1e-8, mads_k, ifelse(is.finite(sds_k) & sds_k > 1e-8, sds_k, 1.0))
      r_k <- rowSums(sweep(sweep(sub_k, 2L, mu_k, "-"), 2L, sc_k, "/")^2)
      n_core_k <- max(min(length(ik_tr), p_d + 2L), as.integer(floor(0.75 * length(ik_tr))))
      sub_k_core <- sub_k[order(r_k, na.last = TRUE)[seq_len(min(length(ik_tr), n_core_k))], , drop = FALSE]
      cov_k <- tryCatch(stats::cov(sub_k_core), error = function(e) matrix(0.0, nrow = p_d, ncol = p_d))
      cov_k[!is.finite(cov_k)] <- 0.0
      if (!is.null(S_pool_f)) {
        w_self <- if (nrow(sub_k_core) >= (p_d + 3L)) 0.75 else 0.50
        cov_blend <- w_self * cov_k + (1.0 - w_self) * S_pool_f
      } else {
        cov_blend <- cov_k
      }
      reg_k <- max(1e-4, 0.05 * max(mean(diag(cov_blend), na.rm = TRUE), 1e-4))
      S_k <- 0.5 * (cov_blend + t(cov_blend)) + diag(reg_k, nrow = p_d, ncol = p_d)
      d2_k <- tryCatch(
        stats::mahalanobis(sub_k, center = mu_k, cov = S_k),
        error = function(e) rep(0.0, length(ik_tr))
      )
      drop_k <- is.finite(d2_k) & (d2_k > q95_thr) & !is_type_vec[ik_tr]
      if (any(drop_k) && sum(!drop_k) >= 2L) {
        keep[ik_tr[drop_k]] <- FALSE
      }
    }
  }
  keep
}
