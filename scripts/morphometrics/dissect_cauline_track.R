if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/dissect_cauline_track.R
# Separate cauline OTU track aggregation and heteroblasty divergence metrics
# ==============================================================================

compute_separate_cauline_otu_track <- function(efa_df, specimen_df, sym_names = NULL) {
  if (is.null(specimen_df) || nrow(specimen_df) == 0L) return(specimen_df)
  n_spec <- nrow(specimen_df)
  cauline_leaf_count <- integer(n_spec)
  cauline_solidity_med <- rep(NA_real_, n_spec)
  cauline_dissection_index_med <- rep(NA_real_, n_spec)
  cauline_aspect_ratio_med <- rep(NA_real_, n_spec)
  cauline_blade_length_mm <- rep(NA_real_, n_spec)
  basal_vs_cauline_length_ratio <- rep(NA_real_, n_spec)
  basal_vs_cauline_solidity_delta <- rep(NA_real_, n_spec)
  cauline_cie_l_q85_sheet_delta <- rep(NA_real_, n_spec)
  basal_to_cauline_wool_delta <- rep(NA_real_, n_spec)
  heteroblasty_efa_distance <- rep(NA_real_, n_spec)
  Cauline_PC1 <- rep(NA_real_, n_spec)
  Cauline_PC2 <- rep(NA_real_, n_spec)
  Cauline_PC3 <- rep(NA_real_, n_spec)
  ECT_PC1 <- rep(NA_real_, n_spec)
  ECT_PC2 <- rep(NA_real_, n_spec)
  ECT_PC3 <- rep(NA_real_, n_spec)
  Cauline_ECT_PC1 <- rep(NA_real_, n_spec)
  Cauline_ECT_PC2 <- rep(NA_real_, n_spec)
  ect_lobe_complexity <- rep(NA_real_, n_spec)
  cauline_ect_lobe_complexity <- rep(NA_real_, n_spec)
  heteroblasty_ect_distance <- rep(NA_real_, n_spec)

  if (!is.null(efa_df) && nrow(efa_df) > 0L) {
    if (is.null(sym_names) || length(sym_names) == 0L) {
      sym_names <- grep("^[AD][0-9]+$", names(efa_df), value = TRUE)
    }
    sym_names <- intersect(as.character(sym_names), names(efa_df))
    ect_sym_vecs <- attr(efa_df, "ect_sym_vecs")
    pids_leaf <- if ("plant_individual_id" %in% names(efa_df)) ifelse(is.na(efa_df$plant_individual_id), 0L, as.integer(efa_df$plant_individual_id)) else rep(0L, nrow(efa_df))
    pids_spec <- if ("plant_individual_id" %in% names(specimen_df)) ifelse(is.na(specimen_df$plant_individual_id), 0L, as.integer(specimen_df$plant_individual_id)) else rep(0L, n_spec)
    leaf_keys <- paste0(trimws(as.character(efa_df$catalogNumber)), "__", pids_leaf)
    spec_keys <- paste0(trimws(as.character(specimen_df$catalogNumber)), "__", pids_spec)
    idx_by_spec <- split(seq_len(nrow(efa_df)), leaf_keys)

    med_col <- function(df_sub, cn, digits = 4L, pos_only = FALSE) {
      if (is.null(df_sub) || nrow(df_sub) == 0L || !(cn %in% names(df_sub))) return(NA_real_)
      v <- suppressWarnings(as.numeric(df_sub[[cn]]))
      v <- if (pos_only) v[is.finite(v) & v > 0.0] else v[is.finite(v)]
      if (length(v) == 0L) NA_real_ else round(stats::median(v), digits)
    }
    row_cal_mask <- function(df_sub) {
      if (is.null(df_sub) || nrow(df_sub) == 0L) return(logical(0))
      ok <- if ("is_calibrated_mm" %in% names(df_sub)) {
        v_c <- df_sub$is_calibrated_mm
        if (is.character(v_c) || is.factor(v_c)) (!is.na(v_c) & (toupper(trimws(as.character(v_c))) %in% c("1", "TRUE", "T", "YES"))) else (!is.na(v_c) & as.logical(v_c))
      } else {
        rep(TRUE, nrow(df_sub))
      }
      if ("scale_source" %in% names(df_sub)) {
        ss_v <- tolower(trimws(as.character(df_sub$scale_source)))
        ss_v[is.na(ss_v)] <- ""
        ok <- ok & !(ss_v %in% c("mp_estimated", "predicted_from_megapixels", "tick_lattice_unboxed", "uncalibrated"))
      }
      ok
    }
    spec_cal_ok <- row_cal_mask(specimen_df)

    caul_harm_list <- vector("list", n_spec)
    caul_ect_list <- vector("list", n_spec)
    for (i in seq_len(n_spec)) {
      r_idx <- idx_by_spec[[spec_keys[i]]]
      if (length(r_idx) == 0L) next
      grp <- efa_df[r_idx, , drop = FALSE]
      keep_sub <- rep(TRUE, nrow(grp))

      is_type_grp <- if ("is_type_specimen" %in% names(grp)) (!is.na(grp$is_type_specimen) & as.logical(grp$is_type_specimen)) else rep(FALSE, nrow(grp))
      if ("is_contour_outlier" %in% names(grp) && any((!as.logical(grp$is_contour_outlier[keep_sub])) | is_type_grp[keep_sub], na.rm = TRUE)) {
        keep_sub <- keep_sub & ((!is.na(grp$is_contour_outlier) & !as.logical(grp$is_contour_outlier)) | is_type_grp)
      }
      grp <- grp[keep_sub, , drop = FALSE]
      r_idx_kept <- r_idx[keep_sub]
      lz <- if ("leaf_zone" %in% names(grp)) tolower(trimws(as.character(grp$leaf_zone))) else rep("basal", nrow(grp))
      lz[is.na(lz) | !(lz %in% c("basal", "cauline"))] <- "basal"
      basal_rows <- grp[lz == "basal", , drop = FALSE]
      cauline_rows <- grp[lz == "cauline", , drop = FALSE]
      bas_orig_idx <- r_idx_kept[lz == "basal"]
      caul_orig_idx <- r_idx_kept[lz == "cauline"]
      n_bas <- nrow(basal_rows)
      n_caul <- nrow(cauline_rows)
      cauline_leaf_count[i] <- as.integer(n_caul)
      if ("basal_leaf_count" %in% names(specimen_df)) specimen_df$basal_leaf_count[i] <- as.integer(n_bas)
      if ("leaf_count" %in% names(specimen_df)) specimen_df$leaf_count[i] <- as.integer(n_bas + n_caul)
      if (n_bas == 0L && n_caul >= 1L && ("morph_tier_source" %in% names(specimen_df))) {
        specimen_df$morph_tier_source[i] <- "Cauline_Only_Fallback"
      }

      prim_rows <- if (n_bas >= 1L) basal_rows else cauline_rows
      if (nrow(prim_rows) >= 1L) {
        ECT_PC1[i] <- med_col(prim_rows, "ECT_PC1", 6L)
        ECT_PC2[i] <- med_col(prim_rows, "ECT_PC2", 6L)
        ECT_PC3[i] <- med_col(prim_rows, "ECT_PC3", 6L)
        ect_lobe_complexity[i] <- med_col(prim_rows, "ect_lobe_complexity", 4L)
        if ("n_lamina_components" %in% names(prim_rows) && any(!is.na(prim_rows$n_lamina_components) & as.integer(prim_rows$n_lamina_components) > 1L)) {
          if ("solidity" %in% names(specimen_df)) specimen_df$solidity[i] <- med_col(prim_rows, "solidity", 4L)
          if ("q25_solidity" %in% names(specimen_df) && "solidity" %in% names(prim_rows)) {
            sv_p <- suppressWarnings(as.numeric(prim_rows$solidity)); sv_p <- sv_p[is.finite(sv_p)]
            if (length(sv_p) >= 1L) specimen_df$q25_solidity[i] <- round(as.numeric(stats::quantile(sv_p, probs = 0.25, names = FALSE)), 4)
          }
          if ("dissection_index" %in% names(specimen_df)) specimen_df$dissection_index[i] <- med_col(prim_rows, "dissection_index", 4L)
          if ("dissected_leaf_fraction" %in% names(specimen_df) && "is_dissected" %in% names(prim_rows)) {
            dv_p <- as.logical(prim_rows$is_dissected); dv_p <- dv_p[!is.na(dv_p)]
            if (length(dv_p) >= 1L) specimen_df$dissected_leaf_fraction[i] <- round(mean(dv_p), 4)
          }
        }
      }

      if (n_caul >= 1L) {
        cauline_solidity_med[i] <- med_col(cauline_rows, "solidity", 4L)
        cauline_dissection_index_med[i] <- med_col(cauline_rows, "dissection_index", 4L)
        cauline_aspect_ratio_med[i] <- med_col(cauline_rows, "aspect_ratio", 4L)
        cauline_ect_lobe_complexity[i] <- med_col(cauline_rows, "ect_lobe_complexity", 4L)
        caul_cal <- cauline_rows[row_cal_mask(cauline_rows), , drop = FALSE]
        bas_cal <- if (n_bas >= 1L) basal_rows[row_cal_mask(basal_rows), , drop = FALSE] else basal_rows
        c_bl <- med_col(caul_cal, "blade_length_mm", 3L, pos_only = TRUE)
        c_cs <- med_col(caul_cal, "centroid_size", 3L, pos_only = TRUE)
        c_cs_px <- med_col(cauline_rows, "centroid_size_px", 3L, pos_only = TRUE)
        c_apx <- med_col(cauline_rows, "area_px", 1L, pos_only = TRUE)
        cauline_blade_length_mm[i] <- if (is.finite(c_bl)) c_bl else NA_real_

        b_bl <- if (n_bas >= 1L) med_col(bas_cal, "blade_length_mm", 4L, pos_only = TRUE) else NA_real_
        if (!is.finite(b_bl) && isTRUE(spec_cal_ok[i]) && ("blade_length_mm" %in% names(specimen_df))) {
          v_b <- suppressWarnings(as.numeric(specimen_df$blade_length_mm[i]))
          if (is.finite(v_b) && v_b > 0.0) b_bl <- v_b
        }
        b_cs <- if (n_bas >= 1L) med_col(bas_cal, "centroid_size", 4L, pos_only = TRUE) else NA_real_
        if (!is.finite(b_cs) && isTRUE(spec_cal_ok[i]) && ("centroid_size" %in% names(specimen_df))) {
          v_c <- suppressWarnings(as.numeric(specimen_df$centroid_size[i]))
          if (is.finite(v_c) && v_c > 0.0) b_cs <- v_c
        }
        b_cs_px <- if (n_bas >= 1L) med_col(basal_rows, "centroid_size_px", 3L, pos_only = TRUE) else NA_real_
        b_apx <- if (n_bas >= 1L) med_col(basal_rows, "area_px", 1L, pos_only = TRUE) else NA_real_
        if (is.finite(b_bl) && is.finite(c_bl) && c_bl > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_bl / c_bl, 4)
        } else if (is.finite(b_cs) && is.finite(c_cs) && c_cs > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_cs / c_cs, 4)
        } else if (is.finite(b_cs_px) && is.finite(c_cs_px) && c_cs_px > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_cs_px / c_cs_px, 4)
        } else if (is.finite(b_apx) && is.finite(c_apx) && c_apx > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(sqrt(b_apx / c_apx), 4)
        }

        b_sol <- if ("solidity" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$solidity[i])) else med_col(basal_rows, "solidity", 4L)
        if (is.finite(b_sol) && is.finite(cauline_solidity_med[i])) {
          basal_vs_cauline_solidity_delta[i] <- round(b_sol - cauline_solidity_med[i], 4)
        }

        if ("blade_cie_l_q85_sheet_delta" %in% names(cauline_rows)) {
          c_abax <- if ("is_abaxial_pressed_leaf" %in% names(cauline_rows)) (!is.na(cauline_rows$is_abaxial_pressed_leaf) & as.logical(cauline_rows$is_abaxial_pressed_leaf)) else rep(FALSE, nrow(cauline_rows))
          c_adax_rows <- if (any(!c_abax)) cauline_rows[!c_abax, , drop = FALSE] else cauline_rows
          c_lq85 <- med_col(c_adax_rows, "blade_cie_l_q85_sheet_delta", 4L)
          if (!is.finite(c_lq85)) c_lq85 <- med_col(cauline_rows, "blade_cie_l_q85_sheet_delta", 4L)
          cauline_cie_l_q85_sheet_delta[i] <- c_lq85
          b_lq85 <- if ("blade_cie_l_q85_sheet_delta" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$blade_cie_l_q85_sheet_delta[i])) else if (n_bas >= 1L) med_col(basal_rows, "blade_cie_l_q85_sheet_delta", 4L) else NA_real_
          if (is.finite(c_lq85) && is.finite(b_lq85)) {
            basal_to_cauline_wool_delta[i] <- round(c_lq85 - b_lq85, 4)
          }
        }

        if (length(sym_names) >= 1L) {
          caul_mat_i <- as.matrix(cauline_rows[, sym_names, drop = FALSE])
          storage.mode(caul_mat_i) <- "double"
          if (any(complete.cases(caul_mat_i))) {
            caul_med_h <- apply(caul_mat_i, 2L, stats::median, na.rm = TRUE)
            caul_med_h[!is.finite(caul_med_h)] <- 0.0
            caul_harm_list[[i]] <- caul_med_h
            if (n_bas >= 1L) {
              bas_mat_i <- as.matrix(basal_rows[, sym_names, drop = FALSE])
              storage.mode(bas_mat_i) <- "double"
              if (any(complete.cases(bas_mat_i))) {
                bas_med_h <- apply(bas_mat_i, 2L, stats::median, na.rm = TRUE)
                bas_med_h[!is.finite(bas_med_h)] <- 0.0
                d_efa <- sqrt(sum((bas_med_h - caul_med_h)^2, na.rm = TRUE))
                if (is.finite(d_efa)) {
                  heteroblasty_efa_distance[i] <- round(d_efa, 6)
                }
              }
            }
          }
        }

        if (!is.null(ect_sym_vecs) && length(caul_orig_idx) >= 1L) {
          c_vecs <- Filter(Negate(is.null), ect_sym_vecs[caul_orig_idx])
          if (length(c_vecs) >= 1L) {
            c_em <- do.call(rbind, c_vecs)
            c_med_e <- apply(c_em, 2L, stats::median, na.rm = TRUE)
            c_med_e[!is.finite(c_med_e)] <- 0.0
            caul_ect_list[[i]] <- c_med_e
            if (n_bas >= 1L && length(bas_orig_idx) >= 1L) {
              b_vecs <- Filter(Negate(is.null), ect_sym_vecs[bas_orig_idx])
              if (length(b_vecs) >= 1L) {
                b_em <- do.call(rbind, b_vecs)
                b_med_e <- apply(b_em, 2L, stats::median, na.rm = TRUE)
                b_med_e[!is.finite(b_med_e)] <- 0.0
                d_ect <- sqrt(mean((b_med_e - c_med_e)^2, na.rm = TRUE))
                if (is.finite(d_ect)) heteroblasty_ect_distance[i] <- round(d_ect, 6)
              }
            }
          }
        } else if (all(c("ECT_PC1", "ECT_PC2", "ECT_PC3") %in% names(cauline_rows))) {
          c_epc <- c(
            med_col(cauline_rows, "ECT_PC1", 6L),
            med_col(cauline_rows, "ECT_PC2", 6L),
            med_col(cauline_rows, "ECT_PC3", 6L)
          )
          c_epc[!is.finite(c_epc)] <- 0.0
          caul_ect_list[[i]] <- c_epc
          if (n_bas >= 1L) {
            b_epc <- c(
              med_col(basal_rows, "ECT_PC1", 6L),
              med_col(basal_rows, "ECT_PC2", 6L),
              med_col(basal_rows, "ECT_PC3", 6L)
            )
            b_epc[!is.finite(b_epc)] <- 0.0
            ect_dim_ref <- attr(efa_df, "ect_sym_dim")
            ect_dim_use <- if (!is.null(ect_dim_ref) && is.finite(ect_dim_ref) && ect_dim_ref >= 3L) as.numeric(ect_dim_ref) else 5792.0
            d_ect <- sqrt(sum((b_epc - c_epc)^2, na.rm = TRUE) / ect_dim_use)
            if (is.finite(d_ect)) heteroblasty_ect_distance[i] <- round(d_ect, 6)
          }
        }
      }
    }

    auth_spec_v <- if ("authority_weight" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$authority_weight)) else rep(1.0, n_spec)
    auth_spec_v[!is.finite(auth_spec_v) | auth_spec_v < 0.0] <- 0.0
    sb_spec_v <- if ("spatial_balance_weight" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$spatial_balance_weight)) else rep(1.0, n_spec)
    sb_spec_v[!is.finite(sb_spec_v) | sb_spec_v <= 0.0] <- 1.0
    eff_spec_v <- auth_spec_v * sb_spec_v

    caul_spec_idx <- which(cauline_leaf_count >= 1L & !vapply(caul_harm_list, is.null, logical(1)))
    if (length(caul_spec_idx) >= 1L && length(sym_names) >= 1L) {
      M_caul <- do.call(rbind, caul_harm_list[caul_spec_idx])
      storage.mode(M_caul) <- "double"
      ref_c_pos <- which(auth_spec_v[caul_spec_idx] > 0.0)
      if (length(ref_c_pos) < 3L) ref_c_pos <- seq_len(nrow(M_caul))
      for (jc in seq_len(ncol(M_caul))) {
        bad_c <- !is.finite(M_caul[, jc])
        if (any(bad_c)) {
          ok_r <- intersect(ref_c_pos, which(!bad_c))
          M_caul[bad_c, jc] <- if (length(ok_r) > 0L) stats::median(M_caul[ok_r, jc], na.rm = TRUE) else if (any(!bad_c)) stats::median(M_caul[!bad_c, jc], na.rm = TRUE) else 0.0
        }
      }
      c_scores <- NULL
      if (length(ref_c_pos) >= 3L) {
        M_ref_c <- M_caul[ref_c_pos, , drop = FALSE]
        c_vars <- apply(M_ref_c, 2L, stats::var, na.rm = TRUE)
        act_c <- which(is.finite(c_vars) & c_vars > 1e-12)
        if (length(act_c) >= 1L) {
          w_c <- eff_spec_v[caul_spec_idx[ref_c_pos]]
          if (sum(w_c, na.rm = TRUE) <= 1e-12) w_c <- rep(1.0, length(ref_c_pos))
          wt_c <- stats::cov.wt(M_ref_c[, act_c, drop = FALSE], wt = w_c, center = TRUE)
          eig_c <- eigen(wt_c$cov, symmetric = TRUE)
          c_scores <- sweep(M_caul[, act_c, drop = FALSE], 2L, wt_c$center, "-") %*% eig_c$vectors
        }
      }
      if (is.null(c_scores)) {
        all_sym_m <- as.matrix(efa_df[, sym_names, drop = FALSE])
        storage.mode(all_sym_m) <- "double"
        all_sym_m[!is.finite(all_sym_m)] <- 0.0
        auth_leaf_v <- if ("authority_weight" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$authority_weight)) else rep(1.0, nrow(efa_df))
        ref_l_pos <- which(is.finite(auth_leaf_v) & auth_leaf_v > 0.0)
        if (length(ref_l_pos) < 2L) ref_l_pos <- seq_len(nrow(all_sym_m))
        if (length(ref_l_pos) >= 2L) {
          c_vars <- apply(all_sym_m[ref_l_pos, , drop = FALSE], 2L, stats::var, na.rm = TRUE)
          act_c <- which(is.finite(c_vars) & c_vars > 1e-12)
          if (length(act_c) >= 1L) {
            mu_all <- colMeans(all_sym_m[ref_l_pos, act_c, drop = FALSE])
            eig_all <- eigen(stats::cov(all_sym_m[ref_l_pos, act_c, drop = FALSE]), symmetric = TRUE)
            c_scores <- sweep(M_caul[, act_c, drop = FALSE], 2L, mu_all, "-") %*% eig_all$vectors
          }
        }
      }
      if (!is.null(c_scores) && nrow(c_scores) == length(caul_spec_idx)) {
        if (ncol(c_scores) >= 1L) {
          pc1_c <- as.numeric(c_scores[, 1L])
          di_c <- cauline_dissection_index_med[caul_spec_idx]
          ok_sc <- is.finite(pc1_c) & is.finite(di_c)
          if (sum(ok_sc) >= 2L && stats::sd(pc1_c[ok_sc]) > 1e-10 && stats::sd(di_c[ok_sc]) > 1e-10) {
            r_sc <- suppressWarnings(stats::cor(pc1_c[ok_sc], di_c[ok_sc]))
            if (is.finite(r_sc) && r_sc < 0.0) c_scores[, 1L] <- -c_scores[, 1L]
          }
        }
        Cauline_PC1[caul_spec_idx] <- if (ncol(c_scores) >= 1L) round(c_scores[, 1L], 6) else 0.0
        Cauline_PC2[caul_spec_idx] <- if (ncol(c_scores) >= 2L) round(c_scores[, 2L], 6) else 0.0
        Cauline_PC3[caul_spec_idx] <- if (ncol(c_scores) >= 3L) round(c_scores[, 3L], 6) else 0.0
      }
    }

    caul_ect_idx <- which(cauline_leaf_count >= 1L & !vapply(caul_ect_list, is.null, logical(1)))
    if (length(caul_ect_idx) >= 1L) {
      M_cect <- do.call(rbind, caul_ect_list[caul_ect_idx])
      storage.mode(M_cect) <- "double"
      M_cect[!is.finite(M_cect)] <- 0.0
      ce_scores <- NULL
      ref_ce_pos <- which(auth_spec_v[caul_ect_idx] > 0.0)
      if (length(ref_ce_pos) < 2L) ref_ce_pos <- seq_len(nrow(M_cect))
      if (length(ref_ce_pos) >= 2L) {
        M_ref_ce <- M_cect[ref_ce_pos, , drop = FALSE]
        ce_vars <- apply(M_ref_ce, 2L, stats::var, na.rm = TRUE)
        act_ce <- which(is.finite(ce_vars) & ce_vars > 1e-12)
        if (length(act_ce) >= 1L) {
          w_ce <- eff_spec_v[caul_ect_idx[ref_ce_pos]]
          if (sum(w_ce, na.rm = TRUE) <= 1e-12) w_ce <- rep(1.0, length(ref_ce_pos))
          w_ce_n <- w_ce / sum(w_ce)
          mu_ce <- as.numeric(crossprod(w_ce_n, M_ref_ce[, act_ce, drop = FALSE]))
          X_ce_ref <- sweep(M_ref_ce[, act_ce, drop = FALSE], 2L, mu_ce, "-")
          Z_ce_ref <- sweep(X_ce_ref, 1L, sqrt(w_ce_n), "*")
          svd_ce <- tryCatch(svd(Z_ce_ref, nu = 0L, nv = min(2L, nrow(Z_ce_ref), ncol(Z_ce_ref))), error = function(e) NULL)
          if (!is.null(svd_ce) && !is.null(svd_ce$v) && ncol(svd_ce$v) >= 1L) {
            ce_scores <- sweep(M_cect[, act_ce, drop = FALSE], 2L, mu_ce, "-") %*% svd_ce$v
            denom_ce <- max(1e-8, 1.0 - sum(w_ce_n^2))
            ev_ce <- (svd_ce$d^2) / denom_ce
            tot_ce <- sum(colSums(Z_ce_ref^2, na.rm = TRUE)) / denom_ce
            attr(specimen_df, "ect_cauline_var_exp") <- round(ev_ce / max(tot_ce, 1e-12) * 100.0, 2)
          }
        }
      }
      if (!is.null(ce_scores) && nrow(ce_scores) == length(caul_ect_idx)) {
        if (ncol(ce_scores) >= 1L) {
          epc1_c <- as.numeric(ce_scores[, 1L])
          elc_c <- cauline_ect_lobe_complexity[caul_ect_idx]
          ok_ce <- is.finite(epc1_c) & is.finite(elc_c)
          if (sum(ok_ce) >= 2L && stats::sd(epc1_c[ok_ce]) > 1e-10 && stats::sd(elc_c[ok_ce]) > 1e-10) {
            r_ce <- suppressWarnings(stats::cor(epc1_c[ok_ce], elc_c[ok_ce]))
            if (is.finite(r_ce) && r_ce < 0.0) ce_scores[, 1L] <- -ce_scores[, 1L]
          }
        }
        Cauline_ECT_PC1[caul_ect_idx] <- if (ncol(ce_scores) >= 1L) round(ce_scores[, 1L], 6) else 0.0
        Cauline_ECT_PC2[caul_ect_idx] <- if (ncol(ce_scores) >= 2L) round(ce_scores[, 2L], 6) else 0.0
      } else {
        Cauline_ECT_PC1[caul_ect_idx] <- 0.0
        Cauline_ECT_PC2[caul_ect_idx] <- 0.0
      }
    }
  }

  specimen_df$cauline_leaf_count <- as.integer(cauline_leaf_count)
  specimen_df$cauline_solidity_med <- as.numeric(cauline_solidity_med)
  specimen_df$cauline_dissection_index_med <- as.numeric(cauline_dissection_index_med)
  specimen_df$cauline_aspect_ratio_med <- as.numeric(cauline_aspect_ratio_med)
  specimen_df$cauline_blade_length_mm <- as.numeric(cauline_blade_length_mm)
  specimen_df$basal_vs_cauline_length_ratio <- as.numeric(basal_vs_cauline_length_ratio)
  specimen_df$basal_vs_cauline_solidity_delta <- as.numeric(basal_vs_cauline_solidity_delta)
  specimen_df$cauline_cie_l_q85_sheet_delta <- as.numeric(cauline_cie_l_q85_sheet_delta)
  specimen_df$basal_to_cauline_wool_delta <- as.numeric(basal_to_cauline_wool_delta)
  specimen_df$heteroblasty_efa_distance <- as.numeric(heteroblasty_efa_distance)
  specimen_df$Cauline_PC1 <- as.numeric(Cauline_PC1)
  specimen_df$Cauline_PC2 <- as.numeric(Cauline_PC2)
  specimen_df$Cauline_PC3 <- as.numeric(Cauline_PC3)
  if (!isTRUE(attr(specimen_df, "ect_otu_projected"))) {
    specimen_df$ECT_PC1 <- as.numeric(ECT_PC1)
    specimen_df$ECT_PC2 <- as.numeric(ECT_PC2)
    specimen_df$ECT_PC3 <- as.numeric(ECT_PC3)
    specimen_df$ect_lobe_complexity <- as.numeric(ect_lobe_complexity)
  }
  specimen_df$Cauline_ECT_PC1 <- as.numeric(Cauline_ECT_PC1)
  specimen_df$Cauline_ECT_PC2 <- as.numeric(Cauline_ECT_PC2)
  specimen_df$cauline_ect_lobe_complexity <- as.numeric(cauline_ect_lobe_complexity)
  specimen_df$heteroblasty_ect_distance <- as.numeric(heteroblasty_ect_distance)
  if (is.null(attr(specimen_df, "ect_basal_var_exp"))) {
    attr(specimen_df, "ect_basal_var_exp") <- attr(efa_df, "ect_basal_var_exp")
  }
  specimen_df
}
