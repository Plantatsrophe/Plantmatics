if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_subspace_mcd_cov.R
# Morphometric subspace selection, robust medoid, FastMCD, and pooled covariance
# ==============================================================================

.select_morphometric_subspace_core <- function(X, group_vec = NULL, min_p = 2L, max_p = NULL, aux_df = NULL, allow_supervised_f = FALSE) {
  X <- as.matrix(X)
  M <- ncol(X); N <- nrow(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(M))
  if (M <= min_p) return(colnames(X))
  valid_g <- if (!is.null(group_vec)) stats::na.omit(as.character(group_vec)) else character(0)
  valid_g <- valid_g[nzchar(valid_g) & !(valid_g %in% c("Unassigned_Packera", "Unknown"))]
  g_tab <- table(valid_g); g_tab_multi <- g_tab[g_tab >= 2L]
  N_min <- if (length(g_tab_multi) >= 2L) min(as.integer(g_tab_multi)) else N
  P_max <- if (N_min >= 6L) max(min_p, floor(N_min / 3)) else M
  if (length(g_tab_multi) >= 4L) {
    df_w_multi <- max(1L, as.integer(N) - length(g_tab_multi))
    n_q25 <- as.integer(stats::quantile(as.numeric(g_tab_multi), probs = 0.25, names = FALSE))
    eff_n <- max(N_min, min(n_q25, floor(df_w_multi / max(2L, length(g_tab_multi)))))
    P_max <- max(P_max, min(M, min(16L, floor(eff_n / 3L))))
  }
  if (isTRUE(allow_supervised_f) && length(g_tab_multi) >= 2L) {
    # Pooled within-group Ledoit-Wolf CDA estimates a single regularized W across df_W = N - G anchors,
    # so a single rare endemic with small N_min must not collapse 13-species CDA down to 2 PCs.
    df_pooled <- max(1L, as.integer(N) - length(g_tab_multi))
    P_max_pooled <- max(min_p, min(16L, floor(df_pooled / 3L)))
    P_max <- max(P_max, min(M, P_max_pooled))
  }
  if (!is.null(max_p) && length(max_p) == 1L && is.finite(as.numeric(max_p)) && as.numeric(max_p) >= 1) {
    P_max <- min(P_max, as.integer(max_p))
  }

  g_fac <- NULL
  if (isTRUE(allow_supervised_f) && !is.null(group_vec) && length(group_vec) == N && length(g_tab_multi) >= 2L) {
    g_chr <- as.character(group_vec)
    g_chr[!(g_chr %in% names(g_tab_multi))] <- NA_character_
    g_fac <- droplevels(factor(g_chr))
    if (nlevels(g_fac) < 2L) g_fac <- NULL
  }

  pc_idx <- grep("^PC[0-9]+$", colnames(X))
  non_pc_idx <- setdiff(seq_len(M), pc_idx)
  kept_pc_cols <- pc_idx
  if (length(pc_idx) >= 3L) {
    v_obs <- apply(X[, pc_idx, drop = FALSE], 2, stats::var, na.rm = TRUE)
    v_prop <- v_obs / max(sum(v_obs, na.rm = TRUE), 1e-12)
    k_m <- length(pc_idx)
    bs_prop <- vapply(seq_len(k_m), function(j) sum(1 / (j:k_m)) / k_m, numeric(1))
    horn_thr <- (1 + sqrt(k_m / max(N, k_m + 1)))^2 / k_m
    keep_pc <- which(is.finite(v_prop) & (v_prop >= pmin(bs_prop * 0.55, horn_thr * 0.55)))
    if (length(keep_pc) < min_p) keep_pc <- sort(unique(c(seq_len(min(min_p, k_m)), keep_pc)))
    rem_pc <- setdiff(which(is.finite(v_prop) & v_prop >= 0.004), keep_pc)
    if (length(rem_pc) > 0L) {
      rescued <- integer(0)
      if (!is.null(aux_df) && is.data.frame(aux_df) && nrow(aux_df) == N) {
        dis_aux_cols <- intersect(c("solidity", "q25_solidity", "dissection_index"), names(aux_df))
        for (rp in rem_pc) {
          v_p <- X[, pc_idx[rp]]
          hit_r <- any(vapply(dis_aux_cols, function(ac) {
            va <- suppressWarnings(as.numeric(aux_df[[ac]]))
            ok <- is.finite(v_p) & is.finite(va)
            if (sum(ok) >= 8L && stats::sd(v_p[ok]) > 1e-8 && stats::sd(va[ok]) > 1e-8) {
              abs(suppressWarnings(stats::cor(v_p[ok], va[ok]))) >= 0.25
            } else FALSE
          }, logical(1)))
          if (isTRUE(hit_r)) rescued <- c(rescued, rp)
        }
      }
      if (!is.null(g_fac)) {
        p_rem <- vapply(rem_pc, function(rp) {
          v_p <- X[, pc_idx[rp]]
          ok <- is.finite(v_p) & !is.na(g_fac)
          if (sum(ok) < 8L || nlevels(droplevels(g_fac[ok])) < 2L) return(1.0)
          pv <- tryCatch(summary(stats::aov(v_p[ok] ~ droplevels(g_fac[ok])))[[1]][["Pr(>F)"]][1L], error = function(e) 1.0)
          if (is.finite(pv)) as.numeric(pv) else 1.0
        }, numeric(1))
        rescued <- union(rescued, rem_pc[stats::p.adjust(p_rem, method = "BH") < 0.01])
      }
      if (length(rescued) > 0L) keep_pc <- sort(unique(c(keep_pc, rescued)))
    }
    if (length(keep_pc) < min_p) keep_pc <- seq_len(min(min_p, k_m))
    kept_pc_cols <- pc_idx[keep_pc]
  }

  cand_idx <- c(kept_pc_cols, non_pc_idx)
  if (length(cand_idx) <= P_max) return(colnames(X)[cand_idx])
  if (!is.null(g_fac)) {
    core_n <- min(length(kept_pc_cols), min(2L, P_max))
    core_idx <- if (core_n > 0L) kept_pc_cols[seq_len(core_n)] else integer(0)
    pool_idx <- setdiff(cand_idx, core_idx)
    if (length(pool_idx) > 0L && length(core_idx) < P_max) {
      f_pool <- vapply(pool_idx, function(ci) {
        v_c <- X[, ci]; ok <- is.finite(v_c) & !is.na(g_fac)
        if (sum(ok) < 6L || nlevels(droplevels(g_fac[ok])) < 2L) return(0.0)
        fv <- tryCatch(summary(stats::aov(v_c[ok] ~ droplevels(g_fac[ok])))[[1]][["F value"]][1L], error = function(e) 0.0)
        if (is.finite(fv)) as.numeric(fv) else 0.0
      }, numeric(1))
      pick_pool <- pool_idx[order(f_pool, decreasing = TRUE)][seq_len(min(length(pool_idx), P_max - length(core_idx)))]
      final_idx <- cand_idx[cand_idx %in% c(core_idx, pick_pool)]
      return(colnames(X)[final_idx])
    }
  }
  colnames(X)[cand_idx[seq_len(P_max)]]
}

select_morphometric_subspace <- function(X, group_vec = NULL, min_p = 2L, max_p = NULL, aux_df = NULL, allow_supervised_f = FALSE) {
  sel <- .select_morphometric_subspace_core(X, group_vec = group_vec, min_p = min_p, max_p = max_p, aux_df = aux_df, allow_supervised_f = allow_supervised_f)
  message(sprintf("Subspace Selection (Horn Parallel Analysis / Broken-Stick & Williams & Titus 1988 P <= floor(N_min / 3)): retained P = %d / %d features.", length(sel), ncol(as.matrix(X))))
  sel
}

.prepare_anchor_feature_matrix <- function(df, feature_matrix = NULL, max_p = NULL, group_vec = NULL) {
  n <- nrow(df)
  if (is.null(feature_matrix)) {
    auto_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
    if (length(auto_cols) == 0L) auto_cols <- .discover_numeric_sorted_cols(names(df), "^[AD][0-9]+$")
    if (length(auto_cols) == 0L) return(NULL)
    feature_matrix <- df[, auto_cols, drop = FALSE]
  }
  X <- as.matrix(feature_matrix)
  storage.mode(X) <- "double"
  if (nrow(X) != n || ncol(X) == 0L) return(NULL)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))

  primary_cols <- grep("^PC[0-9]+$|^[AD]_?[0-9]+$", colnames(X))
  if (length(primary_cols) == 0L) primary_cols <- seq_len(ncol(X))
  has_feat <- apply(X[, primary_cols, drop = FALSE], 1, function(r) all(is.finite(r)))
  if (sum(has_feat) < 2L) {
    has_feat <- apply(X, 1, function(r) any(is.finite(r)))
  }

  keep_cols <- which(vapply(seq_len(ncol(X)), function(j) {
    v <- X[has_feat, j]
    v <- v[is.finite(v)]
    length(v) >= 3L && stats::var(v) > 1e-10
  }, logical(1)))
  if (length(keep_cols) == 0L) return(list(X = NULL, has_feat = has_feat))

  X_use <- X[, keep_cols, drop = FALSE]
  for (j in seq_len(ncol(X_use))) {
    bad_j <- !is.finite(X_use[, j])
    if (any(bad_j)) {
      med_j <- stats::median(X_use[!bad_j, j], na.rm = TRUE)
      X_use[bad_j, j] <- if (is.finite(med_j)) med_j else 0.0
    }
  }
  if (!any(grepl("^PC[0-9]+$", colnames(X_use))) && all(grepl("^[ABCD]_?[0-9]+$", colnames(X_use)))) {
    h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(X_use))))
    if (!any(is.na(h_nums))) X_use <- X_use[, order(h_nums, colnames(X_use)), drop = FALSE]
  }
  dyn_cap <- if (!is.null(max_p) && length(max_p) == 1L && is.finite(as.numeric(max_p)) && as.numeric(max_p) >= 1) {
    as.integer(max_p)
  } else {
    max(2L, min(ncol(X_use), floor(max(sum(has_feat), 4L) / 4L), 15L))
  }
  X_use <- X_use[, seq_len(min(ncol(X_use), dyn_cap)), drop = FALSE]
  list(X = X_use, has_feat = has_feat)
}

.compute_robust_taxon_medoid <- function(X_sub, S_inv = NULL) {
  X_sub <- as.matrix(X_sub)
  m <- nrow(X_sub)
  p <- ncol(X_sub)
  if (m == 0L) return(rep(0.0, p))
  if (m == 1L) return(as.numeric(X_sub[1, ]))
  if (m <= 4L) return(as.numeric(apply(X_sub, 2, stats::median, na.rm = TRUE)))
  if (m <= 600L) {
    if (is.null(S_inv)) {
      mads <- apply(X_sub, 2, stats::mad, na.rm = TRUE)
      sds <- apply(X_sub, 2, stats::sd, na.rm = TRUE)
      sc <- ifelse(is.finite(mads) & mads > 1e-8, mads, ifelse(is.finite(sds) & sds > 1e-8, sds, 1.0))
      Z <- sweep(X_sub, 2, sc, "/")
      dmat <- as.matrix(stats::dist(Z, method = "euclidean"))
    } else {
      # Pairwise Mahalanobis distances under precision matrix S_inv
      G <- X_sub %*% S_inv %*% t(X_sub)
      dg <- diag(G)
      d2mat <- outer(dg, dg, "+") - 2.0 * G
      d2mat[d2mat < 0.0] <- 0.0
      dmat <- as.matrix(sqrt(d2mat))
    }
    row_d <- rowSums(dmat, na.rm = TRUE)
  } else {
    ref_idx <- unique(round(seq(1L, m, length.out = 400L)))
    if (is.null(S_inv)) {
      mads <- apply(X_sub, 2, stats::mad, na.rm = TRUE)
      sds <- apply(X_sub, 2, stats::sd, na.rm = TRUE)
      sc <- ifelse(is.finite(mads) & mads > 1e-8, mads, ifelse(is.finite(sds) & sds > 1e-8, sds, 1.0))
      Z <- sweep(X_sub, 2, sc, "/")
      Z_ref <- Z[ref_idx, , drop = FALSE]
      d2mat <- outer(rowSums(Z * Z), rowSums(Z_ref * Z_ref), "+") - 2.0 * tcrossprod(Z, Z_ref)
    } else {
      XS <- X_sub %*% S_inv
      X_ref <- X_sub[ref_idx, , drop = FALSE]
      d_all <- rowSums(XS * X_sub)
      d2mat <- outer(d_all, d_all[ref_idx], "+") - 2.0 * tcrossprod(XS, X_ref)
    }
    d2mat[d2mat < 0.0] <- 0.0
    row_d <- rowSums(sqrt(d2mat), na.rm = TRUE)
  }
  k_core <- max(2L, ceiling(0.60 * m))
  core_idx <- order(row_d)[seq_len(min(m, k_core))]
  as.numeric(colMeans(X_sub[core_idx, , drop = FALSE]))
}

.mcd_consistency_factor <- function(alpha, p) {
  alpha <- min(0.99, max(0.50, as.numeric(alpha)))
  denom <- stats::pchisq(stats::qchisq(alpha, df = p), df = p + 2)
  if (is.finite(denom) && denom > 1e-8) (alpha / denom) else 1.0
}

.fit_fast_mcd_cov <- function(X_mat, q_used = NULL, seed = 42L, allow_parallel = TRUE) {
  X_mat <- as.matrix(X_mat)
  n_m <- nrow(X_mat)
  p_m <- ncol(X_mat)
  if (n_m < max(12L, p_m + 2L) || !requireNamespace("MASS", quietly = TRUE)) return(NULL)
  if (is.null(q_used)) {
    q_used <- min(n_m - 1L, max(floor(0.80 * n_m), floor((n_m + p_m + 1L) / 2L)))
  }
  if (n_m <= 250L) {
    return(tryCatch(MASS::cov.rob(X_mat, method = "mcd", quantile.used = q_used), error = function(e) NULL))
  }
  # Partitioned parallel FastMCD (Rousseeuw & Van Driessen 1999) + 2 full-sample C-steps across all n_m rows
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({ if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv) }, add = TRUE)
  set.seed(as.integer(seed))
  q_frac <- min(0.95, max(0.50, q_used / n_m))
  n_blk <- min(32L, max(2L, as.integer(ceiling(n_m / 120L))))
  blk_cap <- min(n_m, n_blk * 350L)
  pool_idx <- if (blk_cap < n_m) round(seq(1L, n_m, length.out = blk_cap)) else seq_len(n_m)
  perm_idx <- pool_idx[sample.int(length(pool_idx))]
  blocks <- split(perm_idx, rep_len(seq_len(n_blk), length(perm_idx)))
  fit_block <- function(b_idx) {
    Xb <- X_mat[b_idx, , drop = FALSE]
    nb <- nrow(Xb)
    qb <- min(nb - 1L, max(floor(q_frac * nb), floor((nb + p_m + 1L) / 2L)))
    tryCatch(MASS::cov.rob(Xb, method = "mcd", quantile.used = qb), error = function(e) NULL)
  }
  cand_fits <- if (isTRUE(allow_parallel) && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) {
    packera_mclapply(blocks, fit_block, mc.cores = resolve_packera_workers(64L, n_blk))
  } else {
    lapply(blocks, fit_block)
  }
  cand_fits <- Filter(function(f) !is.null(f) && !is.null(f$cov) && !is.null(f$center), cand_fits)
  if (length(cand_fits) == 0L) return(NULL)
  jit_diag <- diag(1e-6, nrow = p_m, ncol = p_m)
  best_ld <- Inf; best_res <- NULL
  for (cf in cand_fits) {
    c_cur <- as.numeric(cf$center)
    s_cur <- as.matrix(cf$cov) + jit_diag
    ok_c <- TRUE
    for (c_step in 1:2) {
      d2 <- tryCatch(stats::mahalanobis(X_mat, center = c_cur, cov = s_cur), error = function(e) NULL)
      if (is.null(d2) || sum(is.finite(d2)) < q_used) { ok_c <- FALSE; break }
      inl <- order(d2, na.last = TRUE)[seq_len(q_used)]
      X_inl <- X_mat[inl, , drop = FALSE]
      c_cur <- colMeans(X_inl)
      s_raw <- stats::cov(X_inl)
      s_cur <- s_raw + jit_diag
    }
    if (!ok_c) next
    ld <- tryCatch(as.numeric(determinant(s_cur, logarithm = TRUE)$modulus), error = function(e) Inf)
    if (is.finite(ld) && ld < best_ld) {
      best_ld <- ld
      best_res <- list(center = c_cur, cov = s_raw)
    }
  }
  if (!is.null(best_res)) best_res else cand_fits[[1L]]
}

.estimate_pooled_within_cov <- function(X_use, taxa_vec, det_tier, valid_pool) {
  p <- ncol(X_use)
  u_pool_tx <- setdiff(unique(taxa_vec[valid_pool]), "Unassigned_Packera")
  build_resid <- function(tx) {
    idx_gold <- which(valid_pool & (taxa_vec == tx) & (det_tier == "Tier_1_Gold"))
    idx_gs <- which(valid_pool & (taxa_vec == tx) & (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")))
    idx_tx <- if (length(idx_gs) >= 3L) idx_gs else which(valid_pool & (taxa_vec == tx))
    if (length(idx_tx) < 3L) return(NULL)
    X_t <- X_use[idx_tx, , drop = FALSE]
    med_t <- if (length(idx_gold) >= 3L) {
      .compute_robust_taxon_medoid(X_use[idx_gold, , drop = FALSE])
    } else {
      .compute_robust_taxon_medoid(X_t)
    }
    diff_t <- sweep(X_t, 2, med_t, "-")
    nz_mask <- rowSums(diff_t^2) > 1e-12
    if (sum(nz_mask) >= 2L) diff_t <- diff_t[nz_mask, , drop = FALSE]
    if (length(idx_gs) < 3L && nrow(diff_t) >= 5L) {
      r_norm <- sqrt(rowSums(diff_t^2))
      keep_r <- r_norm <= stats::quantile(r_norm, probs = 0.75, na.rm = TRUE)
      if (sum(keep_r) >= 3L) diff_t <- diff_t[keep_r, , drop = FALSE]
    }
    df_scale <- sqrt(nrow(X_t) / max(1.0, nrow(X_t) - 1.0))
    diff_t * df_scale
  }
  resid_blocks <- if (length(u_pool_tx) >= 2L && sum(valid_pool) >= 120L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) {
    packera_mclapply(u_pool_tx, build_resid, mc.cores = resolve_packera_workers(64L, length(u_pool_tx)))
  } else {
    lapply(u_pool_tx, build_resid)
  }
  resid_blocks <- Filter(Negate(is.null), resid_blocks)
  if (length(resid_blocks) > 0L) {
    R_all <- do.call(rbind, resid_blocks)
    if (nrow(R_all) >= (p + 2L)) {
      if (nrow(R_all) >= 25L && requireNamespace("MASS", quietly = TRUE)) {
        q_used <- min(nrow(R_all) - 1L, max(floor(0.80 * nrow(R_all)), floor((nrow(R_all) + p + 1L) / 2L)))
        rob <- .fit_fast_mcd_cov(R_all, q_used = q_used)
        if (!is.null(rob) && !is.null(rob$cov)) {
          c_mcd <- .mcd_consistency_factor(q_used / nrow(R_all), p)
          return(c_mcd * rob$cov + diag(1e-5, nrow = p, ncol = p))
        }
      }
      return(stats::cov(R_all) + diag(1e-5, nrow = p, ncol = p))
    }
  }
  v_diag <- apply(X_use[valid_pool, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  v_diag[!is.finite(v_diag) | v_diag <= 1e-8] <- 1.0
  diag(v_diag, nrow = p, ncol = p)
}

.compute_type_shrinkage_gamma <- function(n_eff_type, n_ref, type_prior_mult = 2.5, max_gamma = 0.35) {
  n_t <- max(0.0, as.numeric(n_eff_type))
  if (!is.finite(n_t) || n_t <= 0.0) return(0.0)
  n_r <- max(1.0, as.numeric(n_ref))
  min(as.numeric(max_gamma), (type_prior_mult * n_t) / (n_r + type_prior_mult * n_t))
}

.count_effective_type_gatherings <- function(df, row_idx) {
  k_len <- length(row_idx)
  if (k_len == 0L) return(0.0)
  if (!is.null(df) && "exsiccatae_key" %in% names(df)) {
    ek <- trimws(as.character(df$exsiccatae_key[row_idx]))
    valid_ek <- !is.na(ek) & nzchar(ek) & !(tolower(ek) %in% c("na", "none", "null", "unknown"))
    n_distinct_ek <- length(unique(ek[valid_ek]))
    n_no_ek <- sum(!valid_ek)
    return(as.numeric(max(1L, n_distinct_ek + n_no_ek)))
  }
  as.numeric(max(1L, k_len))
}
