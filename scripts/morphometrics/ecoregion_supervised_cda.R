if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/ecoregion_supervised_cda.R
# Track 2 supervised EPA Level III Ecoregion canonical discriminant analysis (CDA)
# ==============================================================================

run_supervised_ecoregion_cda <- function(
  df,
  feature_cols,
  gmm_res = NULL,
  output_summary_csv = "outputs/reports/ecoregion_supervised_cda_summary.csv",
  permutations = 999,
  seed = 42,
  min_anchors_per_ecoregion = 3L
) {
  n <- nrow(df)
  quar_mask <- .is_quarantined_ecoregion_coord(df)
  not_quar <- !quar_mask
  rg <- if ("regional_group" %in% names(df)) standardize_epa_l3_label_r(df$regional_group) else rep(NA_character_, n)

  # If regional_group was absent or contained only legacy non-EPA heuristic strings, infer via polygon lookup
  if (all(is.na(rg))) {
    lat_col <- intersect(c("decimalLatitude"), names(df))
    lon_col <- intersect(c("decimalLongitude"), names(df))
    if ((length(lat_col) > 0L && length(lon_col) > 0L) || any(c("stateProvince", "county") %in% names(df))) {
      lat_v <- if (length(lat_col) > 0L) suppressWarnings(as.numeric(df[[lat_col[1]]])) else rep(NA_real_, n)
      lon_v <- if (length(lon_col) > 0L) suppressWarnings(as.numeric(df[[lon_col[1]]])) else rep(NA_real_, n)
      if (any(quar_mask)) {
        lat_v[quar_mask] <- NA_real_
        lon_v[quar_mask] <- NA_real_
      }
      st_v <- if ("stateProvince" %in% names(df)) df$stateProvince else NULL
      co_v <- if ("county" %in% names(df)) df$county else NULL
      rg <- tryCatch(infer_epa_level3_ecoregion(lat_v, lon_v, state_province = st_v, county = co_v), error = function(e) rg)
    }
  }

  pc_req <- grep("^PC[0-9]+$", feature_cols, value = TRUE)
  pc_req <- pc_req[vapply(pc_req, function(cn) any(is.finite(suppressWarnings(as.numeric(df[[cn]])))), logical(1L))]
  ect_req <- grep("^ECT_PC[0-9]+(_std)?$", feature_cols, value = TRUE)
  ect_req <- ect_req[vapply(ect_req, function(cn) any(is.finite(suppressWarnings(as.numeric(df[[cn]])))), logical(1L))]
  req_cols <- if (length(pc_req) > 0L) pc_req else if (length(ect_req) > 0L) ect_req else feature_cols
  has_pcs <- if (length(req_cols) > 0L) stats::complete.cases(df[, req_cols, drop = FALSE]) else rep(FALSE, n)
  auth_w <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n)
  auth_w[is.na(auth_w) | !is.finite(auth_w)] <- 0.0
  is_pas_row <- if ("is_passive" %in% names(df)) { p_v <- as.logical(df$is_passive); !is.na(p_v) & p_v } else rep(FALSE, n)
  is_cln_row <- if ("is_confident_label_noise" %in% names(df)) { c_v <- as.logical(df$is_confident_label_noise); !is.na(c_v) & c_v } else rep(FALSE, n)
  auth_w[is_pas_row | is_cln_row] <- 0.0
  is_t1 <- if ("morph_tier_source" %in% names(df)) (!is.na(df$morph_tier_source) & df$morph_tier_source == "Basal_Rosette") else rep(TRUE, n)

  # Strictly prefer unquarantined Basal_Rosette Gold/Silver anchors (>= 3 per EPA Level III Ecoregion), with unquarantined fallback only when Gold/Silver < 2 groups
  cand_mask <- !is.na(rg) & has_pcs & is_t1 & (auth_w > 0) & not_quar
  tab_cand <- table(rg[cand_mask])
  active_ecos <- sort(names(tab_cand[tab_cand >= min_anchors_per_ecoregion]))
  if (length(active_ecos) < 2L) {
    cand_t1 <- !is.na(rg) & has_pcs & is_t1 & !is_cln_row & not_quar
    tab_t1 <- table(rg[cand_t1])
    if (sum(tab_t1 >= min_anchors_per_ecoregion) >= 2L) {
      cand_mask <- cand_t1
      active_ecos <- sort(names(tab_t1[tab_t1 >= min_anchors_per_ecoregion]))
    } else {
      cand_all <- !is.na(rg) & has_pcs & !is_cln_row & not_quar
      tab_all <- table(rg[cand_all])
      if (sum(tab_all >= min_anchors_per_ecoregion) >= 2L) {
        cand_mask <- cand_all
        active_ecos <- sort(names(tab_all[tab_all >= min_anchors_per_ecoregion]))
      }
    }
  }

  # Unsupervised Morphotype_K* x Supervised EPA_Ecoregion contingency statistics (excluding NA ecoregions)
  gmm_cl <- if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n) {
    as.integer(gmm_res$classification)
  } else if ("gmm_cluster" %in% names(df)) {
    as.integer(df$gmm_cluster)
  } else {
    rep(NA_integer_, n)
  }
  morph_lbl <- ifelse(is.na(gmm_cl), NA_character_, paste0("Morphotype_K", gmm_cl))
  gmm_opt_k <- length(unique(stats::na.omit(gmm_cl)))
  cont_mask <- !is.na(morph_lbl) & !is.na(rg)
  ari_val <- NA_real_
  fisher_p <- NA_real_
  if (sum(cont_mask) >= 4L) {
    ari_val <- compute_ari_pure(morph_lbl[cont_mask], rg[cont_mask])
    cont_tab <- table(morph_lbl[cont_mask], rg[cont_mask])
    if (nrow(cont_tab) >= 2L && ncol(cont_tab) >= 2L) {
      set.seed(seed)
      use_exact <- (nrow(cont_tab) <= 3L && ncol(cont_tab) <= 3L && sum(cont_tab) <= 180L)
      fisher_p <- tryCatch(
        if (use_exact) stats::fisher.test(cont_tab, workspace = 2e6)$p.value
        else stats::fisher.test(cont_tab, simulate.p.value = TRUE, B = max(2000L, as.integer(permutations)))$p.value,
        error = function(e) tryCatch(stats::fisher.test(cont_tab, simulate.p.value = TRUE, B = 2000L)$p.value, error = function(e2) NA_real_)
      )
    }
  }

  pred_eco <- rep(NA_character_, n)
  post_max <- rep(NA_real_, n)
  can1_out <- rep(NA_real_, n)
  can2_out <- rep(NA_real_, n)
  post_mat <- NULL
  loocv_acc <- NA_real_
  perm_f <- NA_real_
  perm_r2 <- NA_real_
  perm_p <- NA_real_
  var_pct <- c(100.0, 0.0)
  struct_load <- NULL
  centroids_2d <- NULL
  is_train <- rep(FALSE, n)

  if (length(active_ecos) >= 2L) {
    is_train <- cand_mask & (rg %in% active_ecos) & not_quar
    X_full <- as.matrix(df[, feature_cols, drop = FALSE])
    storage.mode(X_full) <- "double"
    for (j in seq_len(ncol(X_full))) {
      bad_j <- !is.finite(X_full[, j])
      if (any(bad_j)) {
        fin_tr_j <- is.finite(X_full[is_train, j])
        fill_v <- if (any(fin_tr_j)) stats::median(X_full[is_train, j][fin_tr_j], na.rm = TRUE) else 0.0
        if (!is.finite(fill_v)) fill_v <- 0.0
        X_full[bad_j, j] <- fill_v
      }
    }

    sel_cols <- if (exists("select_morphometric_subspace", mode = "function")) {
      select_morphometric_subspace(X_full[is_train, , drop = FALSE], group_vec = rg[is_train], allow_supervised_f = FALSE)
    } else {
      n_min_g <- min(as.integer(table(rg[is_train])))
      feature_cols[seq_len(max(1L, min(length(feature_cols), floor(n_min_g / 3L))))]
    }
    X_use <- X_full[, sel_cols, drop = FALSE]
    p_dim <- ncol(X_use)

    # Per-ecoregion 95% chi-square Mahalanobis trimming on active training anchors (exempting nomenclatural types)
    is_type_eco <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, n)
    is_type_eco[is.na(is_type_eco)] <- FALSE
    for (eco_g in active_ecos) {
      idx_g <- which(is_train & rg == eco_g)
      if (length(idx_g) >= 5L) {
        sub_g <- X_use[idx_g, , drop = FALSE]
        d2_g <- tryCatch(stats::mahalanobis(sub_g, colMeans(sub_g), stats::cov(sub_g) + diag(1e-4, p_dim)), error = function(e) rep(0, length(idx_g)))
        drop_g <- !is.na(d2_g) & (d2_g > stats::qchisq(0.95, df = p_dim)) & !is_type_eco[idx_g]
        if (sum(!drop_g) >= min_anchors_per_ecoregion) is_train[idx_g[drop_g]] <- FALSE
      }
    }

    active_idx <- which(is_train)
    X_a <- X_use[active_idx, , drop = FALSE]
    grp_a <- droplevels(factor(rg[active_idx], levels = active_ecos))
    g_n <- nlevels(grp_a)
    priors_u <- rep(1.0 / g_n, g_n)

    if (g_n >= 2L) {
      sds_a <- apply(X_a, 2, stats::sd, na.rm = TRUE)
      keep_sd <- is.finite(sds_a) & (sds_a > 1e-7)
      if (!any(keep_sd)) keep_sd[1] <- TRUE
      X_a_fit <- X_a[, keep_sd, drop = FALSE]
      X_proj_fit <- X_use[, keep_sd, drop = FALSE]
      p_ef <- ncol(X_a_fit)

      sp_w_eco <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df, active_mask = is_train, recompute = TRUE) else rep(1.0, n)
      w_eco_a <- auth_w[active_idx] * sp_w_eco[active_idx]
      if (!any(is.finite(w_eco_a) & w_eco_a > 0.0)) w_eco_a <- rep(1.0, length(active_idx))

      if (exists("compute_weighted_shrinkage_cda_loocv", mode = "function")) {
        X_cand_eco <- as.matrix(df[active_idx, feature_cols, drop = FALSE])
        storage.mode(X_cand_eco) <- "double"
        loocv_acc <- tryCatch(
          compute_weighted_shrinkage_cda_loocv(
            X_a_fit, grp_a, w_a = w_eco_a,
            coords_df = df[active_idx, , drop = FALSE],
            deme_radius_km = 10.0,
            X_cand = X_cand_eco,
            allow_supervised_f = FALSE
          ),
          error = function(e) NA_real_
        )
      }

      eco_lvls <- levels(grp_a)
      gm_e <- colSums(X_a_fit * w_eco_a) / sum(w_eco_a)
      B_e <- matrix(0.0, p_ef, p_ef); W_e <- matrix(0.0, p_ef, p_ef)
      Xc_e <- matrix(0.0, nrow(X_a_fit), p_ef)
      mu_e_mat <- matrix(0.0, g_n, p_ef, dimnames = list(eco_lvls, colnames(X_a_fit)))
      for (ki in seq_len(g_n)) {
        km <- which(grp_a == eco_lvls[ki])
        wk <- w_eco_a[km]; if (sum(wk) <= 0.0) wk <- rep(1.0, length(km))
        mk <- colSums(X_a_fit[km, , drop = FALSE] * wk) / sum(wk)
        mu_e_mat[ki, ] <- mk
        dm <- mk - gm_e; B_e <- B_e + sum(wk) * tcrossprod(dm)
        dx <- sweep(X_a_fit[km, , drop = FALSE], 2L, mk, "-"); Xc_e[km, ] <- dx
        W_e <- W_e + crossprod(dx, dx * wk)
      }
      denom_e <- max(sum(w_eco_a) - sum(vapply(eco_lvls, function(eg) { wk <- w_eco_a[grp_a == eg]; if (sum(wk) > 0) sum(wk^2) / sum(wk) else 1.0 }, numeric(1))), 1.0)
      S_w_e <- W_e / denom_e
      S_reg_e <- if (exists("compute_ledoit_wolf_cov", mode = "function")) {
        compute_ledoit_wolf_cov(S_w_e, X_centered = Xc_e, n_obs = nrow(X_a_fit), weights = w_eco_a)
      } else {
        S_w_e + diag(1e-4, p_ef)
      }
      k_ef <- max(1L, min(p_ef, g_n - 1L))
      eig_e <- solve_shrinkage_cda_eigen(S_reg_e, B_e, max_axes = k_ef)
      Z_full_e <- sweep(X_proj_fit, 2L, gm_e, "-") %*% eig_e$vectors
      cent_full_e <- sweep(mu_e_mat, 2L, gm_e, "-") %*% eig_e$vectors

      lp_mat <- matrix(0.0, nrow = n, ncol = g_n, dimnames = list(NULL, eco_lvls))
      for (ki in seq_len(g_n)) {
        dz <- sweep(Z_full_e, 2L, cent_full_e[ki, ], "-")
        lp_mat[, ki] <- log(priors_u[ki]) - 0.5 * rowSums(dz^2)
      }
      post_mat <- exp(lp_mat - apply(lp_mat, 1L, max, na.rm = TRUE))
      post_mat <- post_mat / pmax(rowSums(post_mat, na.rm = TRUE), 1e-12)
      pred_eco <- eco_lvls[max.col(post_mat, ties.method = "first")]
      post_max <- round(apply(post_mat, 1L, max, na.rm = TRUE), 4)

      Z_raw <- Z_full_e
      if (ncol(Z_raw) == 1L) {
        v1 <- as.numeric(Z_raw[, 1])
        src2 <- if ("PC2" %in% colnames(X_full) && isTRUE(stats::sd(X_full[, "PC2"], na.rm = TRUE) > 1e-7)) {
          X_full[, "PC2"]
        } else if ("ECT_PC2" %in% colnames(X_full) && isTRUE(stats::sd(X_full[, "ECT_PC2"], na.rm = TRUE) > 1e-7)) {
          X_full[, "ECT_PC2"]
        } else if (ncol(X_full) >= 2L) {
          X_full[, 2]
        } else {
          X_full[, 1]
        }
        v2 <- as.numeric(src2) - mean(src2, na.rm = TRUE)
        b1 <- v1 - mean(v1, na.rm = TRUE)
        v2 <- v2 - (sum(v2 * b1, na.rm = TRUE) / max(sum(b1^2, na.rm = TRUE), 1e-12)) * b1
        Z_raw <- cbind(Can1 = v1, Can2 = v2)
      }
      can1_out <- round(as.numeric(Z_raw[, 1]), 5)
      can2_out <- round(as.numeric(Z_raw[, 2]), 5)

      ev_e <- pmax(eig_e$values, 0.0)
      if (length(ev_e) >= 2L) {
        var_pct <- round(100.0 * ev_e[1:2] / max(sum(ev_e), 1e-12), 2)
      } else {
        var_pct <- c(100.0, 0.0)
      }
      Z_act_2d <- cbind(can1_out[active_idx], can2_out[active_idx])
      struct_load <- tryCatch(stats::cor(X_a_fit, Z_act_2d), error = function(e) NULL)
      centroids_2d <- do.call(rbind, lapply(eco_lvls, function(eg) {
        m_e <- active_idx[grp_a == eg]
        w_m <- w_eco_a[grp_a == eg]
        c(stats::weighted.mean(can1_out[m_e], w = w_m, na.rm = TRUE), stats::weighted.mean(can2_out[m_e], w = w_m, na.rm = TRUE))
      }))
      rownames(centroids_2d) <- eco_lvls
      colnames(centroids_2d) <- c("Can1", "Can2")
    }

    # Euclidean PERMANOVA across active EPA Level III Ecoregions
    set.seed(seed)
    if (g_n >= 2L && exists("fast_euclidean_permanova_trace", mode = "function")) {
      ad_fast <- tryCatch(
        fast_euclidean_permanova_trace(X_a, grp_a, permutations = as.integer(permutations), seed = seed),
        error = function(e) NULL
      )
      if (!is.null(ad_fast) && is.finite(ad_fast$R2)) {
        perm_r2 <- round(as.numeric(ad_fast$R2), 4)
        perm_f <- round(as.numeric(ad_fast$F), 4)
        perm_p <- signif(as.numeric(ad_fast$p_value), 4)
      }
    } else if (requireNamespace("vegan", quietly = TRUE) && g_n >= 2L) {
      ad_res <- tryCatch(
        vegan::adonis2(X_a ~ grp_a, permutations = as.integer(permutations), method = "euclidean"),
        error = function(e) NULL
      )
      if (!is.null(ad_res) && nrow(ad_res) >= 1L) {
        perm_r2 <- round(as.numeric(ad_res$R2[1]), 4)
        perm_f <- round(as.numeric(ad_res$F[1]), 4)
        perm_p <- signif(as.numeric(ad_res[["Pr(>F)"]][1]), 4)
      }
    }
  } else if (length(active_ecos) == 1L) {
    pred_eco <- rep(active_ecos[1], n)
  }

  all_present_ecos <- sort(unique(rg[!is.na(rg)]))
  n_na_total <- sum(is.na(rg))
  if (length(all_present_ecos) > 0L) {
    summary_df <- do.call(rbind, lapply(all_present_ecos, function(eg) {
      data.frame(
        ecoregion = eg,
        is_active_cda_group = eg %in% active_ecos,
        n_total_vouchers = sum(!is.na(rg) & rg == eg),
        n_active_anchors = sum(is_train & !is.na(rg) & rg == eg),
        n_unassigned_na_projected = n_na_total,
        loocv_accuracy = loocv_acc,
        permanova_f = perm_f,
        permanova_r2 = perm_r2,
        permanova_p_value = perm_p,
        unsupervised_gmm_k = gmm_opt_k,
        ari_morphotype_vs_ecoregion = ari_val,
        fisher_exact_p_value = signif(fisher_p, 4),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    summary_df <- data.frame(
      ecoregion = NA_character_,
      is_active_cda_group = FALSE,
      n_total_vouchers = 0L,
      n_active_anchors = 0L,
      n_unassigned_na_projected = n_na_total,
      loocv_accuracy = NA_real_,
      permanova_f = NA_real_,
      permanova_r2 = NA_real_,
      permanova_p_value = NA_real_,
      unsupervised_gmm_k = gmm_opt_k,
      ari_morphotype_vs_ecoregion = NA_real_,
      fisher_exact_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(output_summary_csv) && nzchar(output_summary_csv)) {
    dir.create(dirname(output_summary_csv), recursive = TRUE, showWarnings = FALSE)
    tmp_sum <- paste0(output_summary_csv, ".", Sys.getpid(), ".tmp")
    utils::write.csv(summary_df, file = tmp_sum, row.names = FALSE)
    file.rename(tmp_sum, output_summary_csv)
  }

  list(
    ecoregion_vec = rg,
    active_ecoregions = active_ecos,
    is_active = is_train,
    ecoregion_cda_predicted = pred_eco,
    ecoregion_cda_posterior_prob = post_max,
    ecoregion_can1 = can1_out,
    ecoregion_can2 = can2_out,
    post_probs = post_mat,
    loocv_accuracy = loocv_acc,
    permanova_f = perm_f,
    permanova_r2 = perm_r2,
    permanova_p_value = perm_p,
    ari_morphotype_vs_ecoregion = ari_val,
    fisher_exact_p_value = fisher_p,
    variance_pct = var_pct,
    structure_loadings = struct_load,
    centroids = centroids_2d,
    summary_table = summary_df
  )
}
