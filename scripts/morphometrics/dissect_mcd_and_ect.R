if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/dissect_mcd_and_ect.R
# Dissection roughness, regime-stratified MCD screening, and Euler Characteristic Transform (ECT)
# ==============================================================================

compute_dissection_aware_roughness <- function(coe, solidity_vec = NULL, pct_asymmetry_vec = NULL) {
  coe <- as.matrix(coe)
  n_rows <- nrow(coe)
  if (n_rows == 0L || ncol(coe) == 0L) {
    return(list(harmonic_roughness = numeric(0), is_contour_outlier = logical(0)))
  }
  all_h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(coe))))
  max_h <- if (any(is.finite(all_h_nums))) max(all_h_nums, na.rm = TRUE) else 12L
  # Measure high-frequency segmentation jaggedness above the biological pinna-lobe band (n = 4..8)
  min_high_n <- max(9L, min(14L, max_h - 3L))
  if (max_h < 9L) min_high_n <- max(6L, max_h - 1L)

  low_h_cols <- colnames(coe)[!is.na(all_h_nums) & all_h_nums >= 1L & all_h_nums <= 3L]
  high_h_cols <- colnames(coe)[!is.na(all_h_nums) & all_h_nums >= min_high_n]
  high_n <- as.numeric(gsub("[^0-9]", "", high_h_cols))

  low_power <- if (length(low_h_cols) > 0L) rowSums(coe[, low_h_cols, drop = FALSE]^2, na.rm = TRUE) else rep(1.0, n_rows)
  high_power <- if (length(high_h_cols) > 0L) {
    rowSums(sweep(coe[, high_h_cols, drop = FALSE]^2, 2L, (high_n / 2.2)^2, "*"), na.rm = TRUE)
  } else {
    rep(0.0, n_rows)
  }
  harmonic_roughness <- high_power / pmax(low_power, 1e-12)

  pct_asym <- if (!is.null(pct_asymmetry_vec) && length(pct_asymmetry_vec) == n_rows) {
    as.numeric(pct_asymmetry_vec)
  } else {
    rep(0.0, n_rows)
  }
  sol <- if (!is.null(solidity_vec) && length(solidity_vec) == n_rows) {
    suppressWarnings(as.numeric(solidity_vec))
  } else {
    rep(1.0, n_rows)
  }
  sol[!is.finite(sol)] <- 1.0

  # Bilaterally symmetric deeply dissected leaves naturally carry finite energy in n >= 9;
  # scale roughness ceiling smoothly up to 0.28 only when bilateral symmetry passes (pct_asymmetry <= 0.20)
  rough_ceil <- ifelse(
    sol < 0.72 & is.finite(pct_asym) & pct_asym <= 0.20,
    pmin(0.28, 0.15 + 0.24 * (0.72 - sol)),
    0.15
  )
  is_out <- if (n_rows >= 10L) (pct_asym > 0.22) | (harmonic_roughness > rough_ceil) else rep(FALSE, n_rows)
  is_out[is.na(is_out)] <- FALSE
  if (all(is_out)) is_out <- rep(FALSE, n_rows)
  list(harmonic_roughness = harmonic_roughness, is_contour_outlier = is_out)
}

screen_regime_stratified_mcd_outliers <- function(
  pc_scores,
  fit_pos = NULL,
  solidity_vec = NULL,
  is_dissected_vec = NULL,
  k_mcd = 5L,
  p_cutoff = 0.001,
  seed = 42L,
  immune_mask = NULL
) {
  pc_scores <- as.matrix(pc_scores)
  n <- nrow(pc_scores)
  mcd_out <- rep(FALSE, n)
  if (n < 15L || !requireNamespace("MASS", quietly = TRUE)) return(mcd_out)
  k_eff <- max(2L, min(as.integer(k_mcd), ncol(pc_scores)))
  pc_sub <- pc_scores[, seq_len(k_eff), drop = FALSE]

  fit_mask <- rep(FALSE, n)
  if (!is.null(fit_pos) && length(fit_pos) > 0L) {
    valid_fp <- fit_pos[!is.na(fit_pos) & fit_pos >= 1L & fit_pos <= n]
    fit_mask[valid_fp] <- TRUE
  } else {
    fit_mask[] <- TRUE
  }
  if (sum(fit_mask) < 20L) return(mcd_out)

  sol <- if (!is.null(solidity_vec) && length(solidity_vec) == n) suppressWarnings(as.numeric(solidity_vec)) else rep(0.90, n)
  sol[!is.finite(sol)] <- 0.90
  is_dis <- if (!is.null(is_dissected_vec) && length(is_dissected_vec) == n) {
    (!is.na(is_dissected_vec) & as.logical(is_dissected_vec)) | (sol < 0.72)
  } else {
    sol < 0.72
  }

  # Partition into morphological regimes so entire-leaf dominance (>50%) never collapses dissected modes (Hubert & Debruyne 2010)
  regimes <- ifelse(!is_dis, "Entire", "Dissected")
  if (sum(is_dis & (sol < 0.50)) >= 60L && sum(is_dis & (sol >= 0.50)) >= 60L) {
    regimes[is_dis & (sol < 0.50)] <- "Deep_Pinnatifid"
    regimes[is_dis & (sol >= 0.50)] <- "Moderate_Lobed"
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)

  for (reg_nm in unique(regimes)) {
    reg_idx <- which(regimes == reg_nm)
    reg_fit <- which((regimes == reg_nm) & fit_mask)
    if (length(reg_fit) < max(12L, k_eff + 5L)) reg_fit <- reg_idx
    k_reg <- min(k_eff, max(2L, floor((length(reg_fit) - 4L) / 2L)))
    if (length(reg_fit) < (k_reg + 5L) || length(reg_idx) < (k_reg + 5L)) next
    X_reg_all <- pc_sub[reg_idx, seq_len(k_reg), drop = FALSE]
    X_reg_fit <- pc_sub[reg_fit, seq_len(k_reg), drop = FALSE]
    v_ok <- apply(X_reg_fit, 2L, stats::var, na.rm = TRUE) > 1e-10
    if (sum(v_ok) < 2L) next
    X_reg_all <- X_reg_all[, v_ok, drop = FALSE]
    X_reg_fit <- X_reg_fit[, v_ok, drop = FALSE]
    df_reg <- ncol(X_reg_fit)
    var_floor <- pmax(1e-4, 0.05 * apply(X_reg_fit, 2L, stats::var, na.rm = TRUE))
    var_floor[!is.finite(var_floor)] <- 1e-4
    if (nrow(X_reg_fit) >= 30L) {
      q_used <- min(nrow(X_reg_fit) - 1L, max(floor(0.88 * nrow(X_reg_fit)), floor((nrow(X_reg_fit) + df_reg + 1L) / 2L)))
      if (nrow(X_reg_fit) > 1000L && exists(".fit_fast_mcd_cov", mode = "function")) {
        mcd_fit <- tryCatch(.fit_fast_mcd_cov(X_reg_fit, q_used = q_used, seed = as.integer(seed)), error = function(e) NULL)
      } else if (nrow(X_reg_fit) > 1000L) {
        set.seed(as.integer(seed))
        ord_r <- order(X_reg_fit[, 1L])
        sub_r <- sort(ord_r[unique(round(seq(1, nrow(X_reg_fit), length.out = 800L)))])
        q_sub <- min(length(sub_r) - 1L, max(floor((q_used / nrow(X_reg_fit)) * length(sub_r)), floor((length(sub_r) + df_reg + 1L) / 2L)))
        sub_fit <- tryCatch(MASS::cov.mcd(X_reg_fit[sub_r, , drop = FALSE], quantile.used = q_sub), error = function(e) NULL)
        if (!is.null(sub_fit) && !is.null(sub_fit$cov)) {
          ctr_c <- sub_fit$center
          cov_c <- sub_fit$cov + diag(1e-7, df_reg)
          for (c_it in 1:2) {
            d2_c <- tryCatch(stats::mahalanobis(X_reg_fit, center = ctr_c, cov = cov_c), error = function(e) NULL)
            if (is.null(d2_c) || any(!is.finite(d2_c))) break
            h_idx <- order(d2_c)[seq_len(q_used)]
            X_h <- X_reg_fit[h_idx, , drop = FALSE]
            ctr_c <- colMeans(X_h)
            cov_c <- stats::cov(X_h) + diag(1e-7, df_reg)
          }
          mcd_fit <- list(center = ctr_c, cov = cov_c)
        } else {
          mcd_fit <- NULL
        }
      } else {
        set.seed(as.integer(seed))
        mcd_fit <- tryCatch(MASS::cov.mcd(X_reg_fit, quantile.used = q_used), error = function(e) NULL)
      }
    } else {
      mcd_fit <- NULL
      q_used <- nrow(X_reg_fit)
    }
    if (!is.null(mcd_fit) && !is.null(mcd_fit$cov)) {
      alpha_q <- min(0.999, max(0.50, q_used / nrow(X_reg_fit)))
      denom_c <- stats::pgamma(stats::qchisq(alpha_q, df = df_reg) / 2.0, shape = df_reg / 2.0 + 1.0)
      c_fac <- if (is.finite(denom_c) && denom_c > 1e-8) (alpha_q / denom_c) else 1.0
      m_center <- mcd_fit$center
      cov_reg <- (c_fac * mcd_fit$cov) + diag(var_floor, nrow = df_reg, ncol = df_reg)
    } else {
      m_center <- colMeans(X_reg_fit, na.rm = TRUE)
      cov_reg <- stats::cov(X_reg_fit) + diag(var_floor, nrow = df_reg, ncol = df_reg)
    }
    d2 <- tryCatch(
      stats::mahalanobis(X_reg_all, center = m_center, cov = cov_reg),
      error = function(e) rep(0.0, nrow(X_reg_all))
    )
    thr <- stats::qchisq(1.0 - p_cutoff, df = df_reg)
    out_r <- is.finite(d2) & (d2 > thr)
    if (any(out_r) && sum(!out_r) >= 5L) {
      mcd_out[reg_idx[out_r]] <- TRUE
    }
  }
  if (!is.null(immune_mask) && length(immune_mask) == nrow(pc_scores)) {
    mcd_out[which(!is.na(immune_mask) & as.logical(immune_mask))] <- FALSE
  }
  if (sum(!mcd_out) < 5L) mcd_out[] <- FALSE
  mcd_out
}

read_numpy_ect_r <- function(npy_path, num_thresh = NULL, num_dirs = NULL) {
  if (is.null(npy_path) || is.na(npy_path) || !nzchar(trimws(as.character(npy_path))) || !file.exists(npy_path)) {
    return(NULL)
  }
  fsize <- file.info(npy_path)$size
  if (is.na(fsize) || fsize < 80L) return(NULL)
  con <- file(npy_path, "rb")
  on.exit(close(con), add = TRUE)
  # Support full 32x360 float32/float64 LM3 ECT matrices (up to 512 KB)
  raw_b <- readBin(con, what = "raw", n = min(as.integer(fsize), 524288L))
  if (length(raw_b) < 16L || as.integer(raw_b[1L]) != 147L) return(NULL)
  if (rawToChar(raw_b[2:6]) != "NUMPY") return(NULL)
  ver_maj <- as.integer(raw_b[7L])
  if (ver_maj == 1L) {
    hlen <- as.integer(raw_b[9L]) + 256L * as.integer(raw_b[10L])
    data_start <- 10L + hlen
  } else {
    hlen <- as.integer(raw_b[9L]) + 256L * as.integer(raw_b[10L]) + 65536L * as.integer(raw_b[11L]) + 16777216L * as.integer(raw_b[12L])
    data_start <- 12L + hlen
  }
  if (data_start >= length(raw_b)) return(NULL)
  hdr_str <- rawToChar(raw_b[seq.int(if (ver_maj == 1L) 11L else 13L, data_start)])
  is_f8 <- grepl("<f8|float64", hdr_str)
  is_fortran <- grepl("'fortran_order':[[:space:]]*True", hdr_str)
  elem_sz <- if (is_f8) 8L else 4L
  nr_use <- if (!is.null(num_thresh) && is.finite(num_thresh)) as.integer(num_thresh) else 32L
  nc_use <- if (!is.null(num_dirs) && is.finite(num_dirs)) as.integer(num_dirs) else 360L
  sh_m <- regmatches(hdr_str, regexec("'shape'[[:space:]]*:[[:space:]]*\\([[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*\\)", hdr_str))[[1L]]
  if (length(sh_m) >= 3L) {
    nr_hdr <- suppressWarnings(as.integer(sh_m[2L]))
    nc_hdr <- suppressWarnings(as.integer(sh_m[3L]))
    if (is.finite(nr_hdr) && nr_hdr >= 4L && is.finite(nc_hdr) && nc_hdr >= 4L) {
      nr_use <- nr_hdr
      nc_use <- nc_hdr
    }
  }
  n_need <- as.integer(nr_use) * as.integer(nc_use)
  payload <- raw_b[seq.int(data_start + 1L, length(raw_b))]
  if (length(payload) < n_need * elem_sz) return(NULL)
  vals <- readBin(payload, what = "numeric", n = n_need, size = elem_sz, endian = "little")
  if (length(vals) != n_need || any(!is.finite(vals))) return(NULL)
  matrix(vals, nrow = as.integer(nr_use), ncol = as.integer(nc_use), byrow = !is_fortran)
}

.resample_ect_directions_r <- function(ect_mat, target_nr = 32L, target_nc = 360L) {
  if (is.null(ect_mat) || !is.matrix(ect_mat)) return(NULL)
  if (nrow(ect_mat) == target_nr && ncol(ect_mat) == target_nc) return(ect_mat)
  # Periodic linear interpolation along phi in [0, 2*pi) preserving Column 1 = phi = 0 (lamina_tip)
  if (nrow(ect_mat) != target_nr) {
    x_old <- seq(0.0, 1.0, length.out = nrow(ect_mat))
    x_new <- seq(0.0, 1.0, length.out = target_nr)
    ect_r <- matrix(0.0, nrow = target_nr, ncol = ncol(ect_mat))
    for (j in seq_len(ncol(ect_mat))) {
      ect_r[, j] <- stats::approx(x_old, ect_mat[, j], xout = x_new, rule = 2)$y
    }
    ect_mat <- ect_r
  }
  if (ncol(ect_mat) != target_nc) {
    nc_old <- ncol(ect_mat)
    phi_ext <- seq(0.0, 1.0, length.out = nc_old + 1L)
    phi_new <- seq(0.0, 1.0 - 1.0 / target_nc, length.out = target_nc)
    ect_c <- matrix(0.0, nrow = target_nr, ncol = target_nc)
    for (i in seq_len(target_nr)) {
      y_ext <- c(ect_mat[i, ], ect_mat[i, 1L])
      ect_c[i, ] <- stats::approx(phi_ext, y_ext, xout = phi_new, rule = 2)$y
    }
    ect_mat <- ect_c
  }
  ect_mat
}

compute_cycle_ect_r <- function(coo, num_thresh = 32L, num_dirs = 360L) {
  if (is.null(coo) || !is.matrix(coo) || nrow(coo) < 6L || ncol(coo) < 2L) return(NULL)
  xy <- matrix(as.numeric(coo[, 1:2]), ncol = 2L)
  if (any(!is.finite(xy))) return(NULL)
  # Orient tip-up: when coo is already horizontally aligned by Step 03 (+X = lamina_tip or -X = lamina_base),
  # use the exact horizontal midrib axis; otherwise project from vertex 1 (lamina_base) to apical midline cap
  if (isTRUE(attr(coo, "lm3_xax_aligned"))) {
    cx0 <- mean(xy[, 1L]); p_base <- xy[1L, ]
    u_tip <- if (p_base[1L] <= cx0) c(1.0, 0.0) else c(-1.0, 0.0)
    v_len <- 1.0
    v_axis <- u_tip
  } else {
    p_base <- xy[1L, ]
    d_base <- sqrt((xy[, 1L] - p_base[1L])^2 + (xy[, 2L] - p_base[2L])^2)
    cap_i <- which(d_base >= 0.985 * max(d_base, na.rm = TRUE))
    v_axis <- c(stats::median(xy[cap_i, 1L]) - p_base[1L], stats::median(xy[cap_i, 2L]) - p_base[2L])
  }
  v_len <- sqrt(sum(v_axis^2))
  if (!is.finite(v_len) || v_len <= 1e-8) {
    u_tip <- c(1.0, 0.0)
  } else {
    u_tip <- v_axis / v_len
  }
  # In Momocs Cartesian space (+X = lamina_tip, +Y = anatomical left, -Y = anatomical right),
  # 90-deg clockwise from u_tip = (ux, uy) is u_right = (uy, -ux), matching Step 02 tip_up_directions().
  u_right <- c(u_tip[2L], -u_tip[1L])
  cx <- mean(xy[, 1L]); cy <- mean(xy[, 2L])
  xy_c <- cbind(xy[, 1L] - cx, xy[, 2L] - cy)
  radii <- sqrt(xy_c[, 1L]^2 + xy_c[, 2L]^2)
  r_max <- max(radii, na.rm = TRUE)
  if (!is.finite(r_max) || r_max <= 1e-8) return(NULL)
  xy_n <- xy_c / r_max
  # Sweep phi in [0, 2*pi) starting at phi = 0 (lamina_tip) and ascending clockwise toward u_right
  phis <- seq(0.0, 2.0 * pi * (1.0 - 1.0 / num_dirs), length.out = num_dirs)
  thresh <- seq(-1.0, 1.0, length.out = num_thresh)
  n_v <- nrow(xy_n)
  next_idx <- c(seq.int(2L, n_v), 1L)
  ect_mat <- matrix(0.0, nrow = num_thresh, ncol = num_dirs)
  for (j in seq_len(num_dirs)) {
    phi_j <- phis[j]
    dir_vec <- cos(phi_j) * u_tip + sin(phi_j) * u_right
    h_v <- as.numeric(xy_n[, 1L] * dir_vec[1L] + xy_n[, 2L] * dir_vec[2L])
    h_e <- pmax(h_v, h_v[next_idx])
    # Exact Euler characteristic of EmbeddedGraph cycle: chi(phi_j, t_i) = #V(h <= t_i) - #E(max(h_u, h_v) <= t_i)
    n_v_le <- findInterval(thresh, sort(h_v))
    n_e_le <- findInterval(thresh, sort(h_e))
    ect_mat[, j] <- as.numeric(n_v_le - n_e_le)
  }
  ect_mat
}

.summarize_and_symmetrize_ect_matrix <- function(ect_mat) {
  if (is.null(ect_mat) || !is.matrix(ect_mat) || nrow(ect_mat) < 4L || ncol(ect_mat) < 4L) {
    return(NULL)
  }
  nr <- nrow(ect_mat); nc <- ncol(ect_mat)
  # Step 02 tip_up_directions() and compute_cycle_ect_r() place phi = 0 (lamina_tip)
  # strictly at Column 1; never re-roll onto wider lateral axes on reniform (W > L) blades!
  # Column 1 is phi = 0 (North/Tip); Columns 2..nc reflect via nc + 2 - j
  refl_cols <- c(1L, seq.int(nc, 2L, by = -1L))
  ect_refl <- ect_mat[, refl_cols, drop = FALSE]
  ect_sym <- 0.5 * (ect_mat + ect_refl)
  tot_norm2 <- sum(ect_mat^2, na.rm = TRUE)
  asym_val <- if (is.finite(tot_norm2) && tot_norm2 > 1e-8) {
    min(1.0, sum((ect_mat - ect_refl)^2, na.rm = TRUE) / (2.0 * tot_norm2))
  } else {
    0.0
  }
  # Total variation along each direction column j: convex shape has TV_j = 2; extra lobes/sinuses increase TV_j and max chi > 1
  tv_cols <- colSums(abs(diff(rbind(0.0, ect_mat, 0.0))), na.rm = TRUE)
  excess_tv <- mean(pmax(0.0, (tv_cols - 2.0) / 2.0), na.rm = TRUE)
  multi_comp <- mean(pmax(0.0, ect_mat - 1.0), na.rm = TRUE)
  lobe_comp <- round(as.numeric(excess_tv + multi_comp), 4)
  half_nc <- floor(nc / 2L) + 1L
  sym_half <- ect_sym[, seq_len(half_nc), drop = FALSE]
  sym_vec <- as.numeric(sign(sym_half) * log1p(abs(sym_half)))
  list(
    ect_mat = ect_mat,
    ect_sym = ect_sym,
    sym_vec = sym_vec,
    lobe_complexity = lobe_comp,
    bilateral_asym = round(as.numeric(asym_val), 5)
  )
}

compute_leaf_ect_features <- function(efa_df, coo_list = NULL, contours_dir = NULL) {
  if (is.null(efa_df) || nrow(efa_df) == 0L) return(efa_df)
  if (!exists("extract_step03_leaf_ect_vectors", mode = "function")) {
    for (p in c("scripts/morphometrics/step03_ect_traits_utils.R", "step03_ect_traits_utils.R")) {
      if (file.exists(p)) { source(p); break }
    }
  }
  ect_vec_res <- extract_step03_leaf_ect_vectors(
    contours_list = coo_list,
    fac_df = efa_df,
    masks_dir = contours_dir,
    num_thresh = 32L,
    num_dirs = 360L
  )
  ect_mats <- ect_vec_res$ect_sym_mats
  sym_vecs <- ect_vec_res$sym_vecs
  lobe_comp <- ect_vec_res$complexity
  bilat_asym <- ect_vec_res$bilateral_asym

  efa_df$ect_lobe_complexity <- lobe_comp
  efa_df$ect_lobe_symmetry <- ect_vec_res$lobe_sym
  efa_df$ect_bilateral_asym <- bilat_asym
  efa_df$ECT_PC1 <- NA_real_
  efa_df$ECT_PC2 <- NA_real_
  efa_df$ECT_PC3 <- NA_real_

  valid_idx <- which(ect_vec_res$valid)
  if (length(valid_idx) >= 2L) {
    M_ect <- ect_vec_res$ect_mat[valid_idx, , drop = FALSE]
    storage.mode(M_ect) <- "double"
    M_ect[!is.finite(M_ect)] <- 0.0
    lz <- if ("leaf_zone" %in% names(efa_df)) tolower(trimws(as.character(efa_df$leaf_zone[valid_idx]))) else rep("basal", length(valid_idx))
    out_v <- if ("is_contour_outlier" %in% names(efa_df)) as.logical(efa_df$is_contour_outlier[valid_idx]) else rep(FALSE, length(valid_idx))
    typ_v <- if ("is_type_specimen" %in% names(efa_df)) (!is.na(efa_df$is_type_specimen[valid_idx]) & as.logical(efa_df$is_type_specimen[valid_idx])) else rep(FALSE, length(valid_idx))
    out_v[is.na(out_v) | typ_v] <- FALSE
    wt_v <- if ("authority_weight" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$authority_weight[valid_idx])) else rep(1.0, length(valid_idx))
    wt_v[!is.finite(wt_v) | wt_v < 0.0] <- 0.0
    sb_v <- if ("spatial_balance_weight" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$spatial_balance_weight[valid_idx])) else rep(1.0, length(valid_idx))
    sb_v[!is.finite(sb_v) | sb_v <= 0.0] <- 1.0
    eff_wt_v <- wt_v * sb_v

    bas_pos <- which((is.na(lz) | lz == "basal") & !out_v)
    if (length(bas_pos) < 2L) bas_pos <- which(is.na(lz) | lz == "basal")
    if (length(bas_pos) < 2L) bas_pos <- seq_along(valid_idx)
    ref_pos <- bas_pos[wt_v[bas_pos] > 0.0]
    if (length(ref_pos) < 2L) ref_pos <- bas_pos

    M_fit <- M_ect[ref_pos, , drop = FALSE]
    w_fit <- if (any(eff_wt_v[ref_pos] > 0.0)) eff_wt_v[ref_pos] else wt_v[ref_pos]
    if (sum(w_fit, na.rm = TRUE) <= 1e-12) w_fit <- rep(1.0, nrow(M_fit))
    w_norm <- w_fit / sum(w_fit)
    col_v <- apply(M_fit, 2L, stats::var, na.rm = TRUE)
    act_c <- which(is.finite(col_v) & col_v > 1e-10)
    if (length(act_c) >= 1L) {
      mu_fit <- as.numeric(crossprod(w_norm, M_fit[, act_c, drop = FALSE]))
      X_fit <- sweep(M_fit[, act_c, drop = FALSE], 2L, mu_fit, "-")
      Z_fit <- sweep(X_fit, 1L, sqrt(w_norm), "*")
      svd_e <- tryCatch(svd(Z_fit, nu = 0L, nv = min(3L, nrow(Z_fit), ncol(Z_fit))), error = function(e) NULL)
      if (!is.null(svd_e) && !is.null(svd_e$v) && ncol(svd_e$v) >= 1L) {
        X_all <- sweep(M_ect[, act_c, drop = FALSE], 2L, mu_fit, "-")
        sc_all <- X_all %*% svd_e$v
        sc_ref <- sc_all[ref_pos, 1L]
        lc_ref <- lobe_comp[valid_idx[ref_pos]]
        # Align ECT_PC1 sign positively with ect_lobe_complexity strictly on active basal reference anchors
        if (sum(is.finite(lc_ref)) >= 3L && stats::sd(sc_ref, na.rm = TRUE) > 1e-8 && stats::sd(lc_ref, na.rm = TRUE) > 1e-8) {
          r_sign <- suppressWarnings(stats::cor(sc_ref, lc_ref, use = "complete.obs"))
          if (is.finite(r_sign) && r_sign < 0.0) {
            sc_all[, 1L] <- -sc_all[, 1L]
            svd_e$v[, 1L] <- -svd_e$v[, 1L]
          }
        }
        efa_df$ECT_PC1[valid_idx] <- round(sc_all[, 1L], 6)
        if (ncol(sc_all) >= 2L) efa_df$ECT_PC2[valid_idx] <- round(sc_all[, 2L], 6) else efa_df$ECT_PC2[valid_idx] <- 0.0
        if (ncol(sc_all) >= 3L) efa_df$ECT_PC3[valid_idx] <- round(sc_all[, 3L], 6) else efa_df$ECT_PC3[valid_idx] <- 0.0
        denom_w <- max(1e-8, 1.0 - sum(w_norm^2))
        tot_wvar <- sum(colSums(Z_fit^2, na.rm = TRUE)) / denom_w
        ev_vals <- (svd_e$d^2) / denom_w
        attr(efa_df, "ect_basal_var_exp") <- round(ev_vals / max(tot_wvar, 1e-12) * 100.0, 2)
      }
    }
    attr(efa_df, "ect_sym_dim") <- ncol(M_ect)
  }
  attr(efa_df, "ect_sym_mats") <- ect_mats
  attr(efa_df, "ect_sym_vecs") <- sym_vecs
  attr(efa_df, "ect_coo_list") <- coo_list
  efa_df
}
