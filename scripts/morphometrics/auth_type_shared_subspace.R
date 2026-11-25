if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_type_shared_subspace.R
# Shared-subspace type control and OSR re-evaluation for multi-island/no-EFA vouchers
# ==============================================================================

.eval_type_control_shared_subspace <- function(
  df, sh, no_efa_idx, n, G, active_taxa, taxa_vec, ss_type_rows_map,
  det_tier, is_t1_pure, is_noise_vec, auth_w, emp_type_all,
  basionym_vec_all, n_bas, basionym_entries, gmm_cl, given_idx,
  typical_cl_list, has_efa_row, type_rank, eff_df_vec,
  type_morph_mahal_d2, nearest_type_morph_taxon, nearest_type_basionym,
  nearest_type_basionym_d2, type_morph_concordance, is_typical_morphotype, osr_res
) {
  compute_pool_centroid <- function(rows_sub, X_src, S_inv_src, prefer_efa = FALSE) {
    .compute_type_pool_centroid(rows_sub, X_src, S_inv_src, df, has_efa_row, type_rank, prefer_efa = prefer_efa)
  }
  has_sh_mat <- !is.null(sh$X_sh) && is.matrix(sh$X_sh) && nrow(sh$X_sh) == n && ncol(sh$X_sh) >= 2L
  if (isTRUE(has_sh_mat)) {
    X_sh <- as.matrix(sh$X_sh)
    storage.mode(X_sh) <- "double"
    p_sh <- ncol(X_sh)
    S_inv_sh <- if (!is.null(sh$S_inv_sh) && is.matrix(sh$S_inv_sh) && nrow(sh$S_inv_sh) == p_sh && ncol(sh$S_inv_sh) == p_sh) {
      sh$S_inv_sh
    } else if (!is.null(sh$S_reg_sh) && is.matrix(sh$S_reg_sh) && nrow(sh$S_reg_sh) == p_sh && ncol(sh$S_reg_sh) == p_sh) {
      tryCatch(solve(sh$S_reg_sh), error = function(e) diag(1.0 / pmax(diag(sh$S_reg_sh), 1e-4), nrow = p_sh, ncol = p_sh))
    } else {
      diag(1.0, nrow = p_sh, ncol = p_sh)
    }
    obs_sh <- if (!is.null(sh$obs_counts) && length(sh$obs_counts) == n) {
      as.integer(sh$obs_counts)
    } else if (!is.null(sh$obs_counts) && length(sh$obs_counts) == length(no_efa_idx)) {
      tmp_obs <- rep(as.integer(p_sh), n)
      tmp_obs[no_efa_idx] <- as.integer(sh$obs_counts)
      tmp_obs
    } else {
      rep(as.integer(p_sh), n)
    }
    has_sh_obs <- !is.na(obs_sh) & (obs_sh >= 2L)

    type_centroids_sh <- matrix(0.0, nrow = G, ncol = p_sh, dimnames = list(active_taxa, colnames(X_sh)))
    n_eff_sh_vec <- stats::setNames(rep(1.0, G), active_taxa)
    anchor_rows_sh_list <- stats::setNames(vector("list", G), active_taxa)
    sh_act_idx <- if (!is.null(sh$act_idx)) as.integer(sh$act_idx) else integer(0)

    for (k in seq_len(G)) {
      tx_k <- active_taxa[k]
      tx_rows <- which(taxa_vec == tx_k)
      ss_r <- ss_type_rows_map[[tx_k]]
      if (length(ss_r) > 0L) {
        use_r <- if (any(has_sh_obs[ss_r])) ss_r[has_sh_obs[ss_r]] else ss_r
        cp_sh <- compute_pool_centroid(use_r, X_src = X_sh, S_inv_src = S_inv_sh, prefer_efa = FALSE)
        type_centroids_sh[k, ] <- cp_sh$mu
        n_eff_sh_vec[tx_k] <- cp_sh$n_eff
        anchor_rows_sh_list[[tx_k]] <- cp_sh$pool
      } else {
        gold_sh <- tx_rows[(det_tier[tx_rows] == "Tier_1_Gold") & is_t1_pure[tx_rows] & !is_noise_vec[tx_rows] & has_sh_obs[tx_rows]]
        if (length(gold_sh) == 0L) gold_sh <- tx_rows[(det_tier[tx_rows] == "Tier_1_Gold") & !is_noise_vec[tx_rows] & has_sh_obs[tx_rows]]
        if (length(gold_sh) > 0L) {
          max_w_sh <- max(auth_w[gold_sh], na.rm = TRUE)
          top_sh <- gold_sh[auth_w[gold_sh] >= 0.85 * max_w_sh]
          n_eff_sh_vec[tx_k] <- as.numeric(max(1L, length(top_sh)))
          type_centroids_sh[k, ] <- .compute_robust_taxon_medoid(X_sh[top_sh, , drop = FALSE], S_inv = S_inv_sh)
        } else {
          act_tx_sh <- intersect(tx_rows, sh_act_idx)
          act_tx_sh <- act_tx_sh[has_sh_obs[act_tx_sh]]
          pool_fb <- if (length(act_tx_sh) > 0L) act_tx_sh else if (any(has_sh_obs[tx_rows])) tx_rows[has_sh_obs[tx_rows]] else tx_rows
          n_eff_sh_vec[tx_k] <- as.numeric(max(1L, length(pool_fb)))
          type_centroids_sh[k, ] <- .compute_robust_taxon_medoid(X_sh[pool_fb, , drop = FALSE], S_inv = S_inv_sh)
        }
        anchor_rows_sh_list[[tx_k]] <- integer(0)
      }
    }

    use_emp_bas_sh <- any(grepl("^Morphotype_K", active_taxa)) && length(emp_type_all) > 0L
    bas_sh_list <- list()
    if (use_emp_bas_sh) {
      u_emp_bas <- unique(basionym_vec_all[emp_type_all])
      bas_sh_list <- lapply(u_emp_bas, function(eb_nm) {
        eb_rows <- emp_type_all[basionym_vec_all[emp_type_all] == eb_nm]
        use_eb <- if (any(has_sh_obs[eb_rows])) eb_rows[has_sh_obs[eb_rows]] else eb_rows
        cp_eb <- compute_pool_centroid(use_eb, X_src = X_sh, S_inv_src = S_inv_sh, prefer_efa = FALSE)
        list(basionym = eb_nm, mu = cp_eb$mu, pool = cp_eb$pool, n_eff = cp_eb$n_eff)
      })
    } else if (n_bas > 0L) {
      bas_sh_list <- lapply(basionym_entries, function(be) {
        if (length(be$rows) > 0L) {
          use_b <- if (any(has_sh_obs[be$rows])) be$rows[has_sh_obs[be$rows]] else be$rows
          cp_b <- compute_pool_centroid(use_b, X_src = X_sh, S_inv_src = S_inv_sh, prefer_efa = FALSE)
          list(basionym = be$basionym, mu = cp_b$mu, pool = cp_b$pool, n_eff = cp_b$n_eff)
        } else {
          k_m <- match(be$taxon, active_taxa)
          list(
            basionym = be$basionym,
            mu = if (!is.na(k_m)) type_centroids_sh[k_m, ] else rep(0.0, p_sh),
            pool = integer(0),
            n_eff = if (!is.na(k_m)) n_eff_sh_vec[active_taxa[k_m]] else max(1.0, be$n_eff)
          )
        }
      })
    }

    sp_raw_sh <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else as.character(taxa_vec)
    sp_taxa_sh <- active_taxa[grepl("^Packera[[:space:]]+", active_taxa)]
    sp_cols_sh <- match(sp_taxa_sh, active_taxa)
    core_sh_rows <- which(!is.na(gmm_cl) & (gmm_cl == 1L) & has_sh_obs)
    if (length(core_sh_rows) < 2L) core_sh_rows <- which(has_sh_obs)
    if (length(core_sh_rows) < 1L) core_sh_rows <- seq_len(n)
    mu_core_sh <- colMeans(X_sh[core_sh_rows, , drop = FALSE])
    scale_core_sh <- 1.0 + (1.0 / max(1.0, length(core_sh_rows)))
    chisq_typ_90_sh <- stats::qchisq(0.90, df = max(1L, p_sh))

    for (ri in no_efa_idx) {
      if (is.na(obs_sh[ri]) || obs_sh[ri] < 2L) {
        type_morph_mahal_d2[ri] <- NA_real_
        nearest_type_morph_taxon[ri] <- NA_character_
        nearest_type_basionym[ri] <- NA_character_
        nearest_type_basionym_d2[ri] <- NA_real_
        type_morph_concordance[ri] <- NA_real_
        is_typical_morphotype[ri] <- NA
        osr_res$min_congener_mahal_d2[ri] <- NA_real_
        osr_res$min_congener_taxon[ri] <- NA_character_
        osr_res$open_set_novelty_p[ri] <- NA_real_
        osr_res$open_set_df[ri] <- NA_integer_
        next
      }
      eff_df_vec[ri] <- as.integer(p_sh)
      d2_sh_tx <- numeric(G)
      for (k in seq_len(G)) {
        dv_k <- X_sh[ri, ] - type_centroids_sh[k, ]
        raw_k <- as.numeric(t(dv_k) %*% S_inv_sh %*% dv_k)
        if (!is.finite(raw_k) || raw_k < 0.0) raw_k <- 0.0
        scale_k <- 1.0 + (1.0 / max(1.0, n_eff_sh_vec[active_taxa[k]]))
        d2_sh_tx[k] <- if (ri %in% anchor_rows_sh_list[[active_taxa[k]]]) raw_k else (raw_k / scale_k)
      }
      if (G >= 1L) {
        nearest_type_morph_taxon[ri] <- active_taxa[which.min(d2_sh_tx)]
      }
      gi <- given_idx[ri]
      if (!is.na(gi)) {
        type_morph_mahal_d2[ri] <- round(d2_sh_tx[gi], 4)
        if (G >= 2L) {
          ll_sh <- -0.5 * d2_sh_tx
          lk_sh <- exp(ll_sh - max(ll_sh, na.rm = TRUE))
          lk_sh <- lk_sh / max(sum(lk_sh, na.rm = TRUE), 1e-12)
          type_morph_concordance[ri] <- round(lk_sh[gi], 4)
        } else {
          type_morph_concordance[ri] <- round(1.0 - stats::pchisq(type_morph_mahal_d2[ri], df = max(1L, p_sh)), 4)
        }
        tx_ri <- active_taxa[gi]
        typ_cls_ri <- typical_cl_list[[tx_ri]]
        has_cl_ri <- !is.na(gmm_cl[ri])
        has_d2_ri <- is.finite(type_morph_mahal_d2[ri])
        in_typ_cl <- has_cl_ri && (length(typ_cls_ri) > 0L) && (gmm_cl[ri] %in% typ_cls_ri)
        in_d2_90 <- has_d2_ri && (type_morph_mahal_d2[ri] <= chisq_typ_90_sh)
        if (in_typ_cl || in_d2_90) {
          is_typical_morphotype[ri] <- TRUE
        } else if (has_cl_ri || has_d2_ri) {
          is_typical_morphotype[ri] <- FALSE
        }
      }
      if (length(bas_sh_list) >= 1L) {
        d2_b_sh <- vapply(bas_sh_list, function(be_sh) {
          dv_b <- X_sh[ri, ] - be_sh$mu
          raw_b <- as.numeric(t(dv_b) %*% S_inv_sh %*% dv_b)
          if (!is.finite(raw_b) || raw_b < 0.0) raw_b <- 0.0
          if (ri %in% be_sh$pool) raw_b else (raw_b / (1.0 + (1.0 / max(1.0, be_sh$n_eff))))
        }, numeric(1L))
        b_best <- which.min(d2_b_sh)
        nearest_type_basionym[ri] <- bas_sh_list[[b_best]]$basionym
        nearest_type_basionym_d2[ri] <- round(d2_b_sh[b_best], 4)
      }
      if (length(sp_taxa_sh) >= 2L) {
        d2_sp_ri <- d2_sh_tx[sp_cols_sh]
        min_all_ri <- min(d2_sp_ri, na.rm = TRUE)
        if (!is.finite(min_all_ri) || min_all_ri < 0.0) min_all_ri <- 0.0
        osr_res$open_set_novelty_p[ri] <- round(1.0 - stats::pchisq(min_all_ri, df = max(1L, p_sh)), 6)
        osr_res$open_set_df[ri] <- as.integer(p_sh)
        yi_sp <- if (!is.na(sp_raw_sh[ri]) && grepl("^Packera[[:space:]]+", sp_raw_sh[ri])) sp_raw_sh[ri] else taxa_vec[ri]
        oth_idx <- which(sp_taxa_sh != yi_sp)
        if (length(oth_idx) >= 1L) {
          d_sub <- d2_sp_ri[oth_idx]
          if (any(is.finite(d_sub))) {
            best_o <- oth_idx[which.min(d_sub)]
            osr_res$min_congener_mahal_d2[ri] <- round(as.numeric(d2_sp_ri[best_o]), 4)
            osr_res$min_congener_taxon[ri] <- sp_taxa_sh[best_o]
          }
        }
      } else if (is.na(osr_res$min_congener_taxon[ri])) {
        dv_c <- X_sh[ri, ] - mu_core_sh
        raw_c <- as.numeric(t(dv_c) %*% S_inv_sh %*% dv_c)
        if (!is.finite(raw_c) || raw_c < 0.0) raw_c <- 0.0
        d2_c_sh <- raw_c / scale_core_sh
        osr_res$open_set_novelty_p[ri] <- round(1.0 - stats::pchisq(d2_c_sh, df = max(1L, p_sh)), 6)
        osr_res$open_set_df[ri] <- as.integer(p_sh)
      }
    }
  } else {
    type_morph_mahal_d2[no_efa_idx] <- NA_real_
    nearest_type_morph_taxon[no_efa_idx] <- NA_character_
    nearest_type_basionym[no_efa_idx] <- NA_character_
    nearest_type_basionym_d2[no_efa_idx] <- NA_real_
    type_morph_concordance[no_efa_idx] <- NA_real_
    is_typical_morphotype[no_efa_idx] <- NA
    osr_res$min_congener_mahal_d2[no_efa_idx] <- NA_real_
    osr_res$min_congener_taxon[no_efa_idx] <- NA_character_
    osr_res$open_set_novelty_p[no_efa_idx] <- NA_real_
    osr_res$open_set_df[no_efa_idx] <- NA_integer_
  }
  list(
    eff_df_vec = eff_df_vec,
    type_morph_mahal_d2 = type_morph_mahal_d2,
    nearest_type_morph_taxon = nearest_type_morph_taxon,
    nearest_type_basionym = nearest_type_basionym,
    nearest_type_basionym_d2 = nearest_type_basionym_d2,
    type_morph_concordance = type_morph_concordance,
    is_typical_morphotype = is_typical_morphotype,
    osr_res = osr_res
  )
}
