# ==============================================================================
# Module: scripts/morphometrics/step03_ect_traits_utils.R
# Step 03 Point-1-Preserving Winding, OTU-Anchored ECT PCA & Step 02 Trait Agg
# ==============================================================================

.step03_ect_traits_ensure_deps <- function() {
  if (!exists("discover_step03_photometry_cols", mode = "function") || !exists("STEP03_LEAF_ARCH_MED_COLS")) {
    for (p in c("scripts/morphometrics/step03_config_scale_utils.R", "step03_config_scale_utils.R")) if (file.exists(p)) { source(p); break }
  }
  if (!exists("read_numpy_ect_r", mode = "function") || !exists("compute_cycle_ect_r", mode = "function")) {
    for (p in c("scripts/morphometrics/dissected_leaf_utils.R", "dissected_leaf_utils.R")) if (file.exists(p)) { source(p); break }
  }
}

.step03_parse_abaxial_vec <- function(v, n) {
  if (is.null(v) || length(v) == 0L) return(rep(FALSE, n))
  v <- rep_len(v, n)
  if (is.logical(v)) return(!is.na(v) & v)
  if (is.numeric(v)) return(is.finite(v) & (v == 1))
  s <- tolower(trimws(as.character(v)))
  !is.na(s) & (s %in% c("true", "t", "1", "1.0"))
}

.step03_composite_keys <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  pids <- if ("plant_individual_id" %in% names(df)) { p <- suppressWarnings(as.integer(df$plant_individual_id)); p[is.na(p)] <- 0L; p } else rep(0L, n)
  paste(trimws(as.character(df$catalogNumber)), pids, sep = "\r")
}

ensure_clockwise_keep_start_r <- function(coo) {
  if (is.null(coo)) return(coo)
  attrs <- attributes(coo)
  if (is.data.frame(coo)) coo <- as.matrix(coo)
  if (!is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 3L) return(coo)
  storage.mode(coo) <- "double"
  x <- coo[, 1L]; y <- coo[, 2L]
  if (any(!is.finite(x)) || any(!is.finite(y))) return(coo)
  sa2 <- sum(x * c(y[-1L], y[1L]) - c(x[-1L], x[1L]) * y)
  if (is.finite(sa2) && sa2 > 0.0) {
    out <- coo[c(1L, seq.int(nrow(coo), 2L, by = -1L)), , drop = FALSE]
    colnames(out) <- colnames(coo)
    for (ak in setdiff(names(attrs), c("dim", "dimnames"))) attr(out, ak) <- attrs[[ak]]
    return(out)
  }
  coo
}

extract_step03_leaf_ect_vectors <- function(contours_list, fac_df, masks_dir = NULL, num_thresh = 16L, num_dirs = 36L) {
  .step03_ect_traits_ensure_deps()
  n_leaves <- if (!is.null(fac_df) && is.data.frame(fac_df)) nrow(fac_df) else length(contours_list)
  if (is.null(n_leaves) || n_leaves == 0L) {
    return(list(ect_mat = matrix(numeric(0), 0L, 0L), complexity = numeric(0), lobe_sym = numeric(0),
      bilateral_asym = numeric(0), valid = logical(0), ect_sym_mats = list(), sym_vecs = list(), ect_sym_dim = 0L))
  }
  ect_dirs <- character(0)
  if (!is.null(masks_dir) && length(masks_dir) >= 1L && !is.na(masks_dir[1L]) && nzchar(trimws(as.character(masks_dir[1L])))) {
    md_clean <- sub("/+$", "", trimws(as.character(masks_dir[1L])))
    ect_dirs <- c(ect_dirs, md_clean, file.path(dirname(md_clean), "ect"))
  }
  ect_dirs <- unique(c(ect_dirs, "data/ect")); ect_dirs <- ect_dirs[dir.exists(ect_dirs)]

  raw_mats <- vector("list", n_leaves); ect_mats <- vector("list", n_leaves); sym_vecs <- vector("list", n_leaves)
  complexity <- rep(NA_real_, n_leaves); lobe_sym <- rep(NA_real_, n_leaves)
  bilat_asym <- rep(NA_real_, n_leaves); valid <- rep(FALSE, n_leaves)
  has_df <- !is.null(fac_df) && is.data.frame(fac_df) && nrow(fac_df) == n_leaves

  for (i in seq_len(n_leaves)) {
    cat_i <- if (has_df && "catalogNumber" %in% names(fac_df)) trimws(as.character(fac_df$catalogNumber[i])) else ""
    pid_i <- if (has_df && "plant_individual_id" %in% names(fac_df) && !is.na(fac_df$plant_individual_id[i])) as.integer(fac_df$plant_individual_id[i]) else 0L
    lid_i <- if (has_df && "leaf_id" %in% names(fac_df) && !is.na(fac_df$leaf_id[i])) as.integer(fac_df$leaf_id[i]) else i
    sid_i <- if (has_df && "shape_id" %in% names(fac_df) && !is.na(fac_df$shape_id[i])) trimws(as.character(fac_df$shape_id[i])) else ""
    ep_i <- if (has_df && "ect_path" %in% names(fac_df)) as.character(fac_df$ect_path[i]) else NA_character_
    cand_paths <- if (!is.na(ep_i) && nzchar(trimws(ep_i))) trimws(ep_i) else character(0)
    for (ed in ect_dirs) {
      if (nzchar(sid_i)) cand_paths <- c(cand_paths, file.path(ed, sprintf("%s_ect.npy", sid_i)))
      if (nzchar(cat_i)) {
        if (pid_i > 0L) cand_paths <- c(cand_paths, file.path(ed, sprintf("%s_p%d_leaf%d_ect.npy", cat_i, pid_i, lid_i)))
        cand_paths <- c(cand_paths, file.path(ed, sprintf("%s_p0_leaf%d_ect.npy", cat_i, lid_i)), file.path(ed, sprintf("%s_leaf%d_ect.npy", cat_i, lid_i)))
      }
    }
    for (cp in unique(cand_paths)) {
      if (file.exists(cp)) {
        m_try <- read_numpy_ect_r(cp)
        if (!is.null(m_try)) { raw_mats[[i]] <- m_try; break }
      }
    }
  }

  loaded_dims <- Filter(Negate(is.null), lapply(raw_mats, dim))
  target_nr <- if (length(loaded_dims) > 0L) as.integer(stats::median(vapply(loaded_dims, `[`, integer(1L), 1L))) else as.integer(num_thresh)
  target_nc <- if (length(loaded_dims) > 0L) as.integer(stats::median(vapply(loaded_dims, `[`, integer(1L), 2L))) else as.integer(num_dirs)

  for (i in seq_len(n_leaves)) {
    mat_i <- raw_mats[[i]]
    if (is.null(mat_i) && !is.null(contours_list) && length(contours_list) > 0L) {
      sid_i <- if (has_df && "shape_id" %in% names(fac_df) && !is.na(fac_df$shape_id[i])) as.character(fac_df$shape_id[i]) else if (has_df) rownames(fac_df)[i] else NULL
      coo_i <- if (!is.null(sid_i) && !is.null(names(contours_list)) && (sid_i %in% names(contours_list))) contours_list[[sid_i]] else if (i <= length(contours_list)) contours_list[[i]] else NULL
      mat_i <- compute_cycle_ect_r(coo_i, num_thresh = target_nr, num_dirs = target_nc)
    }
    if (!is.null(mat_i)) {
      mat_i <- .resample_ect_directions_r(mat_i, target_nr = target_nr, target_nc = target_nc)
      s_res <- .summarize_and_symmetrize_ect_matrix(mat_i)
      if (!is.null(s_res)) {
        sv <- as.numeric(s_res$sym_vec); sv[!is.finite(sv)] <- 0.0
        ect_mats[[i]] <- s_res$ect_sym; sym_vecs[[i]] <- sv
        complexity[i] <- s_res$lobe_complexity; bilat_asym[i] <- s_res$bilateral_asym
        lobe_sym[i] <- round(max(0.0, min(1.0, 1.0 - as.numeric(s_res$bilateral_asym))), 5)
        valid[i] <- TRUE
      }
    }
  }

  v_idx <- which(valid)
  m_dim <- if (length(v_idx) > 0L) length(sym_vecs[[v_idx[1L]]]) else as.integer(target_nr * (floor(target_nc / 2L) + 1L))
  ect_mat <- matrix(NA_real_, nrow = n_leaves, ncol = m_dim); storage.mode(ect_mat) <- "double"
  for (vi in v_idx) ect_mat[vi, ] <- sym_vecs[[vi]]
  if (has_df) {
    rn <- if ("shape_id" %in% names(fac_df)) as.character(fac_df$shape_id) else rownames(fac_df)
    if (!is.null(rn) && length(rn) == n_leaves) rownames(ect_mat) <- rn
  }
  list(
    ect_mat = ect_mat, complexity = as.numeric(complexity), lobe_sym = as.numeric(lobe_sym),
    bilateral_asym = as.numeric(bilat_asym), valid = as.logical(valid),
    ect_sym_mats = ect_mats, sym_vecs = sym_vecs, ect_sym_dim = as.integer(m_dim)
  )
}

fit_and_project_otu_ect_pca <- function(leaf_ect_mat, leaf_df, specimen_df, active_otu_mask = NULL) {
  n_leaf <- if (is.data.frame(leaf_df)) nrow(leaf_df) else 0L
  n_spec <- if (is.data.frame(specimen_df)) nrow(specimen_df) else 0L
  ect_src_list <- is.list(leaf_ect_mat) && !is.data.frame(leaf_ect_mat)
  ect_m <- if (ect_src_list) as.matrix(leaf_ect_mat$ect_mat) else as.matrix(leaf_ect_mat)
  storage.mode(ect_m) <- "double"
  m_dim <- ncol(ect_m)

  leaf_valid <- if (ect_src_list && !is.null(leaf_ect_mat$valid) && length(leaf_ect_mat$valid) == n_leaf) {
    !is.na(leaf_ect_mat$valid) & as.logical(leaf_ect_mat$valid)
  } else if (n_leaf > 0L && nrow(ect_m) == n_leaf && m_dim > 0L) {
    apply(ect_m, 1L, function(r) any(is.finite(r)))
  } else rep(FALSE, n_leaf)

  leaf_comp <- if (ect_src_list && !is.null(leaf_ect_mat$complexity) && length(leaf_ect_mat$complexity) == n_leaf) {
    as.numeric(leaf_ect_mat$complexity)
  } else if ("ect_lobe_complexity" %in% names(leaf_df)) {
    suppressWarnings(as.numeric(leaf_df$ect_lobe_complexity))
  } else rep(NA_real_, n_leaf)
  miss_lc <- leaf_valid & !is.finite(leaf_comp)
  if (any(miss_lc) && nrow(ect_m) == n_leaf && m_dim > 0L) {
    leaf_comp[miss_lc] <- round(rowMeans(abs(ect_m[miss_lc, , drop = FALSE]), na.rm = TRUE), 4)
  }
  leaf_lsym <- if (ect_src_list && !is.null(leaf_ect_mat$lobe_sym) && length(leaf_ect_mat$lobe_sym) == n_leaf) {
    as.numeric(leaf_ect_mat$lobe_sym)
  } else if ("ect_lobe_symmetry" %in% names(leaf_df)) {
    suppressWarnings(as.numeric(leaf_df$ect_lobe_symmetry))
  } else if ("ect_bilateral_asym" %in% names(leaf_df)) {
    ba <- suppressWarnings(as.numeric(leaf_df$ect_bilateral_asym))
    ifelse(is.finite(ba), round(pmax(0.0, pmin(1.0, 1.0 - ba)), 5), NA_real_)
  } else rep(NA_real_, n_leaf)
  leaf_basym <- if (ect_src_list && !is.null(leaf_ect_mat$bilateral_asym) && length(leaf_ect_mat$bilateral_asym) == n_leaf) {
    as.numeric(leaf_ect_mat$bilateral_asym)
  } else if ("ect_bilateral_asym" %in% names(leaf_df)) {
    suppressWarnings(as.numeric(leaf_df$ect_bilateral_asym))
  } else ifelse(is.finite(leaf_lsym), round(pmax(0.0, pmin(1.0, 1.0 - leaf_lsym)), 5), NA_real_)

  for (col in c("ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "ect_lobe_symmetry")) {
    specimen_df[[col]] <- rep(NA_real_, n_spec); leaf_df[[col]] <- rep(NA_real_, n_leaf)
  }
  leaf_df$ect_lobe_complexity <- ifelse(leaf_valid, leaf_comp, NA_real_)
  leaf_df$ect_lobe_symmetry <- ifelse(leaf_valid, leaf_lsym, NA_real_)
  if (any(is.finite(leaf_basym)) || ("ect_bilateral_asym" %in% names(leaf_df))) leaf_df$ect_bilateral_asym <- ifelse(leaf_valid, leaf_basym, NA_real_)
  if (ect_src_list) {
    if (!is.null(leaf_ect_mat$sym_vecs)) attr(leaf_df, "ect_sym_vecs") <- leaf_ect_mat$sym_vecs
    if (!is.null(leaf_ect_mat$ect_sym_mats)) attr(leaf_df, "ect_sym_mats") <- leaf_ect_mat$ect_sym_mats
    if (!is.null(leaf_ect_mat$ect_sym_dim)) attr(leaf_df, "ect_sym_dim") <- leaf_ect_mat$ect_sym_dim
  } else if (m_dim > 0L) attr(leaf_df, "ect_sym_dim") <- as.integer(m_dim)

  otu_ect_mat <- matrix(NA_real_, nrow = n_spec, ncol = max(0L, m_dim)); storage.mode(otu_ect_mat) <- "double"
  otu_valid <- rep(FALSE, n_spec)
  if (n_spec == 0L || n_leaf == 0L || m_dim == 0L) {
    attr(specimen_df, "ect_otu_projected") <- TRUE
    return(list(specimen_df = specimen_df, leaf_df = leaf_df, otu_ect_mat = otu_ect_mat, otu_valid = otu_valid, loadings = matrix(0.0, nrow = max(0L, m_dim), ncol = 0L), center = numeric(max(0L, m_dim))))
  }

  grp_l <- split(seq_len(n_leaf), .step03_composite_keys(leaf_df))
  s_keys <- .step03_composite_keys(specimen_df)
  lz_l <- if ("leaf_zone" %in% names(leaf_df)) tolower(trimws(as.character(leaf_df$leaf_zone))) else rep("basal", n_leaf)
  lz_l[is.na(lz_l) | !(lz_l %in% c("basal", "cauline"))] <- "basal"
  mts_s <- if ("morph_tier_source" %in% names(specimen_df)) trimws(as.character(specimen_df$morph_tier_source)) else rep(NA_character_, n_spec)
  otu_comp <- rep(NA_real_, n_spec); otu_lsym <- rep(NA_real_, n_spec)

  for (i in seq_len(n_spec)) {
    idx_all <- grp_l[[s_keys[i]]]
    if (is.null(idx_all) || length(idx_all) == 0L) next
    idx_val <- idx_all[leaf_valid[idx_all]]
    if (length(idx_val) == 0L) next
    want_bas <- is.na(mts_s[i]) || !nzchar(mts_s[i]) || identical(mts_s[i], "Basal_Rosette")
    use_idx <- if (want_bas && ("leaf_zone" %in% names(leaf_df))) { bas_v <- idx_val[lz_l[idx_val] == "basal"]; if (length(bas_v) >= 1L) bas_v else idx_val } else idx_val
    med_v <- if (length(use_idx) == 1L) as.numeric(ect_m[use_idx, ]) else apply(ect_m[use_idx, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
    if (!any(is.finite(med_v))) next
    med_v[!is.finite(med_v)] <- 0.0
    otu_ect_mat[i, ] <- med_v; otu_valid[i] <- TRUE
    lc_u <- leaf_comp[use_idx]; lc_u <- lc_u[is.finite(lc_u)]
    otu_comp[i] <- if (length(lc_u) > 0L) round(stats::median(lc_u), 4) else round(mean(abs(med_v), na.rm = TRUE), 4)
    ls_u <- leaf_lsym[use_idx]; ls_u <- ls_u[is.finite(ls_u)]
    otu_lsym[i] <- if (length(ls_u) > 0L) round(stats::median(ls_u), 5) else NA_real_
  }
  specimen_df$ect_lobe_complexity <- otu_comp; specimen_df$ect_lobe_symmetry <- otu_lsym

  is_basal_otu <- is.na(mts_s) | !nzchar(mts_s) | (mts_s == "Basal_Rosette")
  act_mask <- if (!is.null(active_otu_mask) && length(active_otu_mask) == n_spec) (!is.na(active_otu_mask) & as.logical(active_otu_mask)) else rep(TRUE, n_spec)
  det_tier <- if ("determiner_tier" %in% names(specimen_df)) tolower(trimws(as.character(specimen_df$determiner_tier))) else rep("", n_spec)
  is_bronze <- !is.na(det_tier) & (det_tier == "tier_3_bronze")
  auth_w <- if ("authority_weight" %in% names(specimen_df)) { aw <- suppressWarnings(as.numeric(specimen_df$authority_weight)); aw[!is.finite(aw) | aw < 0.0] <- 0.0; aw } else ifelse(is_bronze, 0.0, 1.0)
  spat_w <- if ("spatial_balance_weight" %in% names(specimen_df)) { sw <- suppressWarnings(as.numeric(specimen_df$spatial_balance_weight)); sw[!is.finite(sw) | sw <= 0.0] <- 1.0; sw } else rep(1.0, n_spec)

  w_eff <- auth_w * spat_w
  w_eff[!(otu_valid & is_basal_otu & act_mask)] <- 0.0
  non_bronze_elig <- otu_valid & is_basal_otu & act_mask & !is_bronze & (w_eff > 0.0)
  if (sum(non_bronze_elig) >= 2L) {
    w_eff[is_bronze] <- 0.0; fit_idx <- which(non_bronze_elig)
  } else {
    promo_elig <- otu_valid & is_basal_otu & act_mask & (w_eff > 0.0)
    if (sum(promo_elig) >= 2L) {
      fit_idx <- which(promo_elig)
    } else {
      fit_idx <- which(otu_valid & is_basal_otu)
      if (length(fit_idx) < 2L) fit_idx <- which(otu_valid)
      w_eff[fit_idx] <- ifelse(auth_w[fit_idx] > 0.0, auth_w[fit_idx] * spat_w[fit_idx], spat_w[fit_idx])
    }
  }

  loadings <- matrix(0.0, nrow = m_dim, ncol = 0L)
  mu_full <- if (any(otu_valid)) colMeans(otu_ect_mat[otu_valid, , drop = FALSE], na.rm = TRUE) else numeric(m_dim)
  if (length(fit_idx) >= 2L) {
    X_fit <- otu_ect_mat[fit_idx, , drop = FALSE]
    w_fit <- w_eff[fit_idx]; if (sum(w_fit, na.rm = TRUE) <= 1e-12) w_fit <- rep(1.0, length(fit_idx))
    w_norm <- w_fit / sum(w_fit)
    mu_full <- as.numeric(crossprod(w_norm, X_fit))
    col_v <- apply(X_fit, 2L, stats::var, na.rm = TRUE)
    act_c <- which(is.finite(col_v) & col_v > 1e-10)
    if (length(act_c) >= 1L) {
      mu_act <- mu_full[act_c]
      Z_fit <- sweep(sweep(X_fit[, act_c, drop = FALSE], 2L, mu_act, "-"), 1L, sqrt(w_norm), "*")
      svd_e <- tryCatch(svd(Z_fit, nu = 0L, nv = min(3L, nrow(Z_fit), ncol(Z_fit))), error = function(e) NULL)
      if (!is.null(svd_e) && !is.null(svd_e$v) && ncol(svd_e$v) >= 1L) {
        v_mat <- svd_e$v; otu_v_idx <- which(otu_valid)
        sc_otu <- sweep(otu_ect_mat[otu_v_idx, act_c, drop = FALSE], 2L, mu_act, "-") %*% v_mat
        sc_f1 <- sc_otu[match(fit_idx, otu_v_idx), 1L]; cp_f <- otu_comp[fit_idx]
        if (sum(is.finite(sc_f1) & is.finite(cp_f)) >= 2L && stats::sd(sc_f1, na.rm = TRUE) > 1e-10 && stats::sd(cp_f, na.rm = TRUE) > 1e-10) {
          r_sg <- suppressWarnings(stats::cor(sc_f1, cp_f, use = "complete.obs"))
          if (is.finite(r_sg) && r_sg < 0.0) { v_mat[, 1L] <- -v_mat[, 1L]; sc_otu[, 1L] <- -sc_otu[, 1L] }
        }
        loadings <- matrix(0.0, nrow = m_dim, ncol = ncol(v_mat)); loadings[act_c, ] <- v_mat
        colnames(loadings) <- paste0("ECT_PC", seq_len(ncol(v_mat)))
        specimen_df$ECT_PC1[otu_v_idx] <- round(sc_otu[, 1L], 6)
        specimen_df$ECT_PC2[otu_v_idx] <- if (ncol(sc_otu) >= 2L) round(sc_otu[, 2L], 6) else 0.0
        specimen_df$ECT_PC3[otu_v_idx] <- if (ncol(sc_otu) >= 3L) round(sc_otu[, 3L], 6) else 0.0
        lv_idx <- which(leaf_valid)
        if (length(lv_idx) > 0L) {
          X_lv <- ect_m[lv_idx, act_c, drop = FALSE]; X_lv[!is.finite(X_lv)] <- 0.0
          sc_lv <- sweep(X_lv, 2L, mu_act, "-") %*% v_mat
          leaf_df$ECT_PC1[lv_idx] <- round(sc_lv[, 1L], 6)
          leaf_df$ECT_PC2[lv_idx] <- if (ncol(sc_lv) >= 2L) round(sc_lv[, 2L], 6) else 0.0
          leaf_df$ECT_PC3[lv_idx] <- if (ncol(sc_lv) >= 3L) round(sc_lv[, 3L], 6) else 0.0
        }
        denom_w <- max(1e-8, 1.0 - sum(w_norm^2))
        var_exp <- round(((svd_e$d^2) / denom_w) / max(sum(colSums(Z_fit^2, na.rm = TRUE)) / denom_w, 1e-12) * 100.0, 2)
        attr(specimen_df, "ect_basal_var_exp") <- var_exp; attr(leaf_df, "ect_basal_var_exp") <- var_exp
      }
    }
  }
  if (ncol(loadings) == 0L) {
    if (any(otu_valid)) specimen_df[otu_valid, c("ECT_PC1", "ECT_PC2", "ECT_PC3")] <- 0.0
    if (any(leaf_valid)) leaf_df[leaf_valid, c("ECT_PC1", "ECT_PC2", "ECT_PC3")] <- 0.0
  }
  attr(specimen_df, "ect_otu_projected") <- TRUE
  list(specimen_df = specimen_df, leaf_df = leaf_df, otu_ect_mat = otu_ect_mat, otu_valid = otu_valid, loadings = loadings, center = mu_full)
}

aggregate_step02_traits_by_otu <- function(leaf_df, grp_indices, prim_zones, cal_flags) {
  .step03_ect_traits_ensure_deps()
  n_g <- length(grp_indices)
  arch_cols <- if (exists("STEP03_LEAF_ARCH_MED_COLS")) STEP03_LEAF_ARCH_MED_COLS else c(
    "base_angle_deg", "apex_angle_deg", "lamina_curvature", "basal_cordate_index",
    "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac"
  )
  adax_phot_cols <- setdiff(discover_step03_photometry_cols(leaf_df), c("basal_zone_cie_l_delta", "abaxial_cie_l_sheet_delta"))

  out <- data.frame(row.names = seq_len(n_g))
  for (ac in arch_cols) out[[ac]] <- rep(NA_real_, n_g)
  out$lamina_trace_length_mm <- rep(NA_real_, n_g); out$abaxial_leaf_fraction <- rep(NA_real_, n_g)
  out$is_abaxial_pressed_leaf <- rep(FALSE, n_g); out$abaxial_cie_l_sheet_delta <- rep(NA_real_, n_g)
  out$basal_zone_cie_l_delta <- rep(NA_real_, n_g)
  for (pc in adax_phot_cols) out[[pc]] <- rep(NA_real_, n_g)
  if (n_g == 0L || is.null(leaf_df) || nrow(leaf_df) == 0L) return(out)

  n_l <- nrow(leaf_df); has_lz <- "leaf_zone" %in% colnames(leaf_df)
  lz_all <- if (has_lz) tolower(trimws(as.character(leaf_df$leaf_zone))) else rep("basal", n_l)
  prim_z <- if (is.null(prim_zones)) rep("basal", n_g) else tolower(trimws(as.character(rep_len(prim_zones, n_g))))
  cal_f <- if (is.null(cal_flags)) rep(FALSE, n_g) else (!is.na(rep_len(cal_flags, n_g)) & as.logical(rep_len(cal_flags, n_g)))

  pull_v <- function(cn) if (cn %in% names(leaf_df)) suppressWarnings(as.numeric(leaf_df[[cn]])) else NULL
  med4 <- function(x) { xf <- x[is.finite(x)]; if (length(xf) > 0L) round(stats::median(xf), 4) else NA_real_ }
  arch_vecs <- stats::setNames(lapply(arch_cols, pull_v), arch_cols)
  trace_v <- pull_v("lamina_trace_length_mm")
  abax_all <- .step03_parse_abaxial_vec(if ("is_abaxial_pressed_leaf" %in% names(leaf_df)) leaf_df$is_abaxial_pressed_leaf else NULL, n_l)
  abax_l_v <- pull_v("abaxial_cie_l_sheet_delta"); basal_l_v <- pull_v("basal_zone_cie_l_delta")
  blade_lsd_v <- pull_v("blade_cie_l_sheet_delta")
  adax_vecs <- stats::setNames(lapply(adax_phot_cols, pull_v), adax_phot_cols)

  for (g in seq_len(n_g)) {
    idx <- grp_indices[[g]]
    if (length(idx) == 0L) next
    p_idx <- if ((prim_z[g] %in% c("basal", "basal_rosette")) && has_lz) {
      b_sub <- idx[!is.na(lz_all[idx]) & lz_all[idx] == "basal"]; if (length(b_sub) > 0L) b_sub else idx
    } else idx

    for (ac in arch_cols) if (!is.null(arch_vecs[[ac]])) out[[ac]][g] <- med4(arch_vecs[[ac]][p_idx])
    if (isTRUE(cal_f[g]) && !is.null(trace_v)) {
      tf <- trace_v[p_idx]; tf <- tf[is.finite(tf)]
      if (length(tf) > 0L) out$lamina_trace_length_mm[g] <- round(as.numeric(stats::quantile(tf, probs = 0.75, type = 7, names = FALSE)), 4)
    }

    abax <- abax_all[p_idx]; out$abaxial_leaf_fraction[g] <- round(mean(abax), 4); out$is_abaxial_pressed_leaf[g] <- isTRUE(any(abax))
    if (!is.null(abax_l_v)) { v_ab <- med4(abax_l_v[p_idx[abax]]); out$abaxial_cie_l_sheet_delta[g] <- if (is.finite(v_ab)) v_ab else med4(abax_l_v[p_idx]) }
    if (!is.finite(out$abaxial_cie_l_sheet_delta[g]) && any(abax) && !is.null(blade_lsd_v)) out$abaxial_cie_l_sheet_delta[g] <- med4(blade_lsd_v[p_idx[abax]])
    if (!is.null(basal_l_v)) out$basal_zone_cie_l_delta[g] <- med4(basal_l_v[p_idx])
    for (pc in adax_phot_cols) { pv <- adax_vecs[[pc]][p_idx]; v_ad <- med4(pv[!abax]); out[[pc]][g] <- if (is.finite(v_ad)) v_ad else med4(pv) }
  }
  out
}
