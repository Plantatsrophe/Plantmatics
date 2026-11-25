if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step04_restore_b.R
# Step 04 raw Step 03 restoration from leaf_efa_harmonics (part B: photometry, PCs, and coordinator).
# ==============================================================================

.restore_raw_step03_photometry_and_pcs_r <- function(state) {
  df <- state$df
  leaf_df <- state$leaf_df
  n_df <- state$n_df
  prim_idx <- state$prim_idx
  matched_rows <- state$matched_rows
  df_is_type <- state$df_is_type
  harm_cols <- state$harm_cols
  restored_mask <- state$restored_mask
  has_fin_in_leaf <- function(cn) .step04_leaf_has_fin_r(leaf_df, cn)
  agg_med_col <- function(cn, digits = NULL, pos_only = FALSE, cal_vec = NULL, subset_vec = NULL) {
    .step04_agg_med_col_r(leaf_df, prim_idx, matched_rows, n_df, cn, digits = digits, pos_only = pos_only, cal_vec = cal_vec, subset_vec = subset_vec)
  }
  has_abax_col <- "is_abaxial_pressed_leaf" %in% names(leaf_df)
  l_abax <- if (has_abax_col) {
    !is.na(leaf_df$is_abaxial_pressed_leaf) & as.logical(leaf_df$is_abaxial_pressed_leaf)
  } else {
    rep(FALSE, nrow(leaf_df))
  }
  l_adax <- !l_abax

  if (has_fin_in_leaf("basal_zone_cie_l_delta")) {
    bz_vals <- agg_med_col("basal_zone_cie_l_delta", digits = 4L)
    ok_bz <- matched_rows & is.finite(bz_vals)
    if (any(ok_bz)) {
      if (!("basal_zone_cie_l_delta" %in% names(df))) df$basal_zone_cie_l_delta <- rep(NA_real_, n_df)
      df$basal_zone_cie_l_delta[ok_bz] <- bz_vals[ok_bz]
      restored_mask[["basal_zone_cie_l_delta"]] <- ok_bz
    }
  }

  phot_cols <- unique(c(
    intersect(c(
      "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
      "blade_cie_l_q85_sheet_delta",
      "blade_cie_a_delta", "blade_cie_b_delta",
      "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta"
    ), names(leaf_df)),
    grep("^blade_cie_", names(leaf_df), value = TRUE)
  ))
  phot_cols <- phot_cols[vapply(phot_cols, has_fin_in_leaf, logical(1L))]
  for (pcol in phot_cols) {
    vals <- agg_med_col(pcol, digits = 4L, subset_vec = l_adax)
    fb_all <- agg_med_col(pcol, digits = 4L)
    use_fb <- matched_rows & !is.finite(vals) & is.finite(fb_all)
    if (any(use_fb)) vals[use_fb] <- fb_all[use_fb]
    ok_p <- matched_rows & is.finite(vals)
    if (any(ok_p)) {
      if (!(pcol %in% names(df))) df[[pcol]] <- rep(NA_real_, n_df)
      df[[pcol]][ok_p] <- vals[ok_p]
      restored_mask[[pcol]] <- ok_p
    }
  }

  if (has_abax_col) {
    abax_frac <- rep(NA_real_, n_df)
    abax_flag <- rep(FALSE, n_df)
    for (i in which(matched_rows)) {
      idx <- prim_idx[[i]]
      if (length(idx) > 0L) {
        abax_frac[i] <- round(mean(l_abax[idx]), 4)
        abax_flag[i] <- any(l_abax[idx])
      }
    }
    if (!("abaxial_leaf_fraction" %in% names(df))) df$abaxial_leaf_fraction <- rep(NA_real_, n_df)
    df$abaxial_leaf_fraction[matched_rows] <- abax_frac[matched_rows]
    restored_mask[["abaxial_leaf_fraction"]] <- matched_rows
    if (!("is_abaxial_pressed_leaf" %in% names(df))) df$is_abaxial_pressed_leaf <- rep(FALSE, n_df)
    df$is_abaxial_pressed_leaf[matched_rows] <- abax_flag[matched_rows]
    restored_mask[["is_abaxial_pressed_leaf"]] <- matched_rows
  }
  if (has_fin_in_leaf("abaxial_cie_l_sheet_delta")) {
    abax_l_vals <- agg_med_col("abaxial_cie_l_sheet_delta", digits = 4L, subset_vec = if (has_abax_col) l_abax else NULL)
    if (!("abaxial_cie_l_sheet_delta" %in% names(df))) df$abaxial_cie_l_sheet_delta <- rep(NA_real_, n_df)
    df$abaxial_cie_l_sheet_delta[matched_rows] <- abax_l_vals[matched_rows]
    restored_mask[["abaxial_cie_l_sheet_delta"]] <- matched_rows
  }

  pca_sym_cols <- grep("^[AD][0-9]+$", harm_cols, value = TRUE)
  if (length(pca_sym_cols) == 0L) {
    pca_sym_cols <- harm_cols
  } else {
    a_c <- grep("^A[0-9]+$", pca_sym_cols, value = TRUE)
    d_c <- grep("^D[0-9]+$", pca_sym_cols, value = TRUE)
    if (length(a_c) > 1L) a_c <- a_c[order(as.integer(sub("^A", "", a_c)))]
    if (length(d_c) > 1L) d_c <- d_c[order(as.integer(sub("^D", "", d_c)))]
    pca_sym_cols <- c(a_c, d_c)
  }

  leaf_pc_cols <- grep("^PC[0-9]+$", names(leaf_df), value = TRUE)
  leaf_pc_cols <- leaf_pc_cols[vapply(leaf_pc_cols, has_fin_in_leaf, logical(1L))]
  if (length(leaf_pc_cols) > 1L) {
    leaf_pc_cols <- leaf_pc_cols[order(as.integer(sub("^PC", "", leaf_pc_cols)))]
  }
  df_pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
  df_pc_cols <- df_pc_cols[vapply(df_pc_cols, function(cn) any(is.finite(suppressWarnings(as.numeric(df[[cn]])))), logical(1L))]
  if (length(df_pc_cols) > 1L) {
    df_pc_cols <- df_pc_cols[order(as.integer(sub("^PC", "", df_pc_cols)))]
  }

  k_target <- if (length(leaf_pc_cols) >= 2L) {
    length(leaf_pc_cols)
  } else if (length(df_pc_cols) >= 2L) {
    length(df_pc_cols)
  } else {
    min(5L, max(2L, length(pca_sym_cols)))
  }
  k_target <- max(2L, min(15L, as.integer(k_target)))
  target_pcs <- paste0("PC", seq_len(k_target))

  leaf_pc_med <- NULL
  if (length(leaf_pc_cols) >= 2L) {
    k_leaf <- min(k_target, length(leaf_pc_cols))
    leaf_pc_med <- matrix(NA_real_, nrow = n_df, ncol = k_target, dimnames = list(NULL, target_pcs))
    for (pj in seq_len(k_leaf)) {
      leaf_pc_med[, pj] <- agg_med_col(leaf_pc_cols[pj], digits = 6L)
    }
  }

  if (length(pca_sym_cols) >= 1L) {
    X_spec <- as.matrix(df[, pca_sym_cols, drop = FALSE])
    storage.mode(X_spec) <- "double"
    otu_has_efa <- matched_rows & stats::complete.cases(X_spec)
    pc_lm <- NULL
    pc_cov <- NULL

    if (any(otu_has_efa) && !is.null(leaf_pc_med)) {
      all_single_leaf <- all(vapply(prim_idx[otu_has_efa], length, integer(1L)) <= 1L)
      if (isTRUE(all_single_leaf) && all(is.finite(leaf_pc_med[otu_has_efa, seq_len(min(2L, ncol(leaf_pc_med))), drop = FALSE]))) {
        pc_lm <- leaf_pc_med[otu_has_efa, , drop = FALSE]
      } else {
        k_leaf <- min(k_target, length(leaf_pc_cols))
        lv_ok <- stats::complete.cases(leaf_df[, c(pca_sym_cols, leaf_pc_cols[seq_len(k_leaf)]), drop = FALSE])
        if (sum(lv_ok) > length(pca_sym_cols) + 2L) {
          lv_h <- as.matrix(leaf_df[lv_ok, pca_sym_cols, drop = FALSE])
          storage.mode(lv_h) <- "double"
          lv_v <- apply(lv_h, 2L, stats::var, na.rm = TRUE)
          fit_h <- names(lv_v[is.finite(lv_v) & lv_v > 1e-10])
          if (length(fit_h) >= 1L && sum(lv_ok) > length(fit_h) + 2L) {
            X_des <- cbind(1.0, lv_h[, fit_h, drop = FALSE])
            Y_pcs <- as.matrix(leaf_df[lv_ok, leaf_pc_cols[seq_len(k_leaf)], drop = FALSE])
            storage.mode(Y_pcs) <- "double"
            qr_x <- tryCatch(qr(X_des), error = function(e) NULL)
            if (!is.null(qr_x) && qr_x$rank == ncol(X_des)) {
              fit_b <- tryCatch(stats::lm.fit(X_des, Y_pcs), error = function(e) NULL)
              if (!is.null(fit_b) && all(is.finite(fit_b$coefficients)) &&
                  max(abs(fit_b$residuals), na.rm = TRUE) < 1e-4) {
                X_sp_des <- cbind(1.0, X_spec[otu_has_efa, fit_h, drop = FALSE])
                proj_lm <- X_sp_des %*% fit_b$coefficients
                pc_lm <- matrix(0.0, nrow = sum(otu_has_efa), ncol = k_target)
                for (pj in seq_len(k_leaf)) {
                  pc_lm[, pj] <- round(proj_lm[, pj], 6)
                }
              }
            }
          }
        }
      }
    }

    if (sum(otu_has_efa) >= 2L) {
      w_auth <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(0.0, n_df)
      if ("authority_weight" %in% names(leaf_df) && any(is.finite(suppressWarnings(as.numeric(leaf_df$authority_weight))))) {
        l_aw <- agg_med_col("authority_weight")
        use_l <- matched_rows & is.finite(l_aw)
        w_auth[use_l] <- l_aw[use_l]
      }
      w_auth[!is.finite(w_auth) | w_auth < 0.0] <- 0.0

      w_sb <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, n_df)
      if ("spatial_balance_weight" %in% names(leaf_df) && any(is.finite(suppressWarnings(as.numeric(leaf_df$spatial_balance_weight))))) {
        l_sb <- agg_med_col("spatial_balance_weight")
        use_sb <- matched_rows & is.finite(l_sb) & (l_sb > 0.0)
        w_sb[use_sb] <- l_sb[use_sb]
      }
      w_sb[!is.finite(w_sb) | w_sb <= 0.0] <- 1.0

      is_tier1_ok <- if ("morph_tier_source" %in% names(df)) {
        is.na(df$morph_tier_source) | !(df$morph_tier_source %in% c("Cauline_Only_Fallback", "Transitional_Passive"))
      } else {
        rep(TRUE, n_df)
      }
      not_out_otu <- if ("is_contour_outlier" %in% names(df)) {
        is.na(df$is_contour_outlier) | !as.logical(df$is_contour_outlier) | df_is_type
      } else {
        rep(TRUE, n_df)
      }
      fit_mask <- otu_has_efa & (w_auth > 0.0) & is_tier1_ok & not_out_otu
      fit_idx <- which(fit_mask)
      w_eff <- w_auth * w_sb

      rot_mat <- NULL
      cen_vec <- NULL
      if (length(fit_idx) >= 2L && sum(w_eff[fit_idx], na.rm = TRUE) > 0.0) {
        cv <- apply(X_spec[fit_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE)
        nz_cols <- names(cv[is.finite(cv) & cv > 1e-10])
        if (length(nz_cols) >= 1L) {
          wt_cov <- tryCatch(stats::cov.wt(X_spec[fit_idx, nz_cols, drop = FALSE], wt = w_eff[fit_idx], center = TRUE), error = function(e) NULL)
          if (!is.null(wt_cov)) {
            eig <- tryCatch(eigen(wt_cov$cov, symmetric = TRUE), error = function(e) NULL)
            if (!is.null(eig)) {
              rot_mat <- eig$vectors
              cen_vec <- wt_cov$center
            }
          }
        }
      }
      if (is.null(rot_mat)) {
        all_idx <- which(otu_has_efa)
        cv <- apply(X_spec[all_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE)
        nz_cols <- names(cv[is.finite(cv) & cv > 1e-10])
        if (length(nz_cols) >= 1L) {
          pr <- tryCatch(stats::prcomp(X_spec[all_idx, nz_cols, drop = FALSE], center = TRUE, scale. = FALSE), error = function(e) NULL)
          if (!is.null(pr)) {
            rot_mat <- pr$rotation
            cen_vec <- pr$center
          }
        }
      }
      if (!is.null(rot_mat) && !is.null(cen_vec)) {
        sc_raw <- sweep(X_spec[otu_has_efa, nz_cols, drop = FALSE], 2L, cen_vec, "-") %*% rot_mat
        pc_cov <- matrix(0.0, nrow = sum(otu_has_efa), ncol = k_target)
        n_copy <- min(k_target, ncol(sc_raw))
        for (pj in seq_len(n_copy)) {
          sc_col <- sc_raw[, pj]
          ref_col <- NULL
          raw_ind_nm <- sprintf("PC%d_raw_individual", pj)
          if (raw_ind_nm %in% names(df)) {
            cand_r <- suppressWarnings(as.numeric(df[[raw_ind_nm]][otu_has_efa]))
            if (sum(is.finite(cand_r)) >= 2L) ref_col <- cand_r
          }
          if (is.null(ref_col) && !is.null(pc_lm) && pj <= ncol(pc_lm)) {
            cand_lm <- pc_lm[, pj]
            if (sum(is.finite(cand_lm)) >= 2L) ref_col <- cand_lm
          }
          if (is.null(ref_col) && !is.null(leaf_pc_med) && pj <= ncol(leaf_pc_med)) {
            cand_l <- leaf_pc_med[otu_has_efa, pj]
            if (sum(is.finite(cand_l)) >= 2L) ref_col <- cand_l
          }
          if (!is.null(ref_col)) {
            ok_s <- is.finite(sc_col) & is.finite(ref_col)
            if (sum(ok_s) >= 2L) {
              dp <- sum((sc_col[ok_s] - mean(sc_col[ok_s])) * (ref_col[ok_s] - mean(ref_col[ok_s])))
              if (is.finite(dp) && dp < 0.0) sc_col <- -sc_col
            }
          }
          pc_cov[, pj] <- round(sc_col, 6)
        }
      }
    }

    pc_restored_mat <- if (!is.null(pc_lm) && !is.null(pc_cov)) {
      if (isTRUE(max(abs(pc_cov - pc_lm), na.rm = TRUE) < 5e-5)) pc_cov else pc_lm
    } else if (!is.null(pc_lm)) {
      pc_lm
    } else {
      pc_cov
    }

    if (!is.null(pc_restored_mat)) {
      for (pj in seq_len(k_target)) {
        p_nm <- target_pcs[pj]
        if (!(p_nm %in% names(df))) df[[p_nm]] <- rep(NA_real_, n_df)
        df[[p_nm]][otu_has_efa] <- pc_restored_mat[, pj]
        df[[p_nm]][matched_rows & !otu_has_efa] <- NA_real_
        restored_mask[[p_nm]] <- matched_rows
      }
    } else if (any(matched_rows & !otu_has_efa)) {
      for (p_nm in intersect(target_pcs, names(df))) {
        df[[p_nm]][matched_rows & !otu_has_efa] <- NA_real_
        restored_mask[[p_nm]] <- matched_rows & !otu_has_efa
      }
    }
  } else if (!is.null(leaf_pc_med)) {
    for (pj in seq_len(min(k_target, length(leaf_pc_cols)))) {
      p_nm <- target_pcs[pj]
      ok_p <- matched_rows & is.finite(leaf_pc_med[, pj])
      if (any(ok_p)) {
        if (!(p_nm %in% names(df))) df[[p_nm]] <- rep(NA_real_, n_df)
        df[[p_nm]][ok_p] <- leaf_pc_med[ok_p, pj]
        restored_mask[[p_nm]] <- ok_p
      }
    }
  }

  list(df = df, restored_mask = restored_mask)
}

.restore_raw_step03_from_leaf_harmonics_db <- function(conn, rid, df) {
  state <- .restore_raw_step03_harmonics_and_scale_r(conn, rid, df)
  if (!isTRUE(state$matched)) {
    return(list(df = state$df, restored_mask = state$restored_mask))
  }
  .restore_raw_step03_photometry_and_pcs_r(state)
}
