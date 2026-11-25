if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_adaptive_anchors.R
# Adaptive reference anchor selection and out-of-fold confident learning pruning
# ==============================================================================

select_adaptive_reference_anchors <- function(
  df,
  feature_matrix = NULL,
  taxon_col = "species_standardized",
  min_anchors_per_taxon = 5L,
  mahal_q_gold_silver = 0.95,
  mahal_q_bronze_promo = 0.80
) {
  if (is.null(df) || nrow(df) == 0L) {
    return(list(
      active_mask = logical(0),
      is_type_anchor = logical(0),
      is_gold_anchor = logical(0),
      is_silver_anchor = logical(0),
      is_promoted_bronze = logical(0),
      is_confident_label_noise = logical(0),
      anchor_tier_source = character(0)
    ))
  }

  n <- nrow(df)
  min_anc <- max(1L, as.integer(min_anchors_per_taxon))
  type_info <- .identify_type_specimen_vec(df)
  type_rank <- type_info$type_rank
  det_tier <- .standardize_determiner_tier_vec(df)
  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  is_t1_pure <- .identify_tier1_pure_mask(df)
  if (sum(is_t1_pure) == 0L) is_t1_pure <- rep(TRUE, n)

  feat_prep <- .prepare_anchor_feature_matrix(df, feature_matrix = feature_matrix, max_p = NULL, group_vec = taxa_vec)
  X_use <- if (!is.null(feat_prep)) feat_prep$X else NULL
  has_feat <- if (!is.null(feat_prep)) feat_prep$has_feat else rep(TRUE, n)
  valid_pool <- is_t1_pure & has_feat

  extralimital_raw <- is_extralimital_or_conflicted_non_type(df)
  extralimital_mask <- rep(FALSE, n)
  if (any(extralimital_raw)) {
    for (tx_chk in unique(taxa_vec[!is.na(taxa_vec) & nzchar(taxa_vec)])) {
      if (sum(valid_pool & (taxa_vec == tx_chk) & !extralimital_raw) >= 2L) {
        extralimital_mask[(taxa_vec == tx_chk) & extralimital_raw] <- TRUE
      }
    }
  }
  is_open_hyb_non_type <- .identify_open_or_hybrid_concept_mask(df) & !(type_rank >= 1L)
  is_het_syn_non_type <- .is_heterotypic_synonym_mask(df) & !(type_rank >= 1L)
  anchor_pool <- valid_pool & !extralimital_mask & !is_open_hyb_non_type

  is_type_valid <- (type_rank >= 1L) & valid_pool
  is_nb_type_valid <- as.logical(type_info$is_name_bearing_type) & valid_pool
  is_nb_type_valid[is.na(is_nb_type_valid)] <- FALSE
  is_gold <- (((det_tier == "Tier_1_Gold") & anchor_pool) | (type_rank >= 1L & valid_pool))
  is_silver <- (det_tier == "Tier_2_Silver") & !is_gold & anchor_pool
  is_gs <- is_gold | is_silver
  is_promoted_bronze <- rep(FALSE, n)

  S_pool <- if (!is.null(X_use) && sum(valid_pool) >= 4L) {
    .estimate_pooled_within_cov(X_use, taxa_vec, det_tier, valid_pool)
  } else {
    NULL
  }
  S_pool_inv <- if (!is.null(S_pool)) {
    tryCatch(solve(S_pool), error = function(e) diag(1.0 / pmax(diag(S_pool), 1e-6), ncol(S_pool)))
  } else {
    NULL
  }

  q25_v <- if ("q25_solidity" %in% names(df)) suppressWarnings(as.numeric(df$q25_solidity)) else rep(NA_real_, n)
  sol_v <- if ("solidity" %in% names(df)) suppressWarnings(as.numeric(df$solidity)) else rep(NA_real_, n)
  dlf_v <- if ("dissected_leaf_fraction" %in% names(df)) suppressWarnings(as.numeric(df$dissected_leaf_fraction)) else rep(NA_real_, n)
  dis_b <- if ("is_dissected" %in% names(df)) (!is.na(df$is_dissected) & as.logical(df$is_dissected)) else rep(FALSE, n)
  eff_sol <- ifelse(is.finite(q25_v), q25_v, sol_v)
  is_dis_mode <- (is.finite(eff_sol) & eff_sol < 0.72) | (is.finite(dlf_v) & dlf_v >= 0.50) | (!is.finite(eff_sol) & dis_b)
  lqw_col <- if ("leaf_quality_weight" %in% names(df)) suppressWarnings(as.numeric(df$leaf_quality_weight)) else NULL

  unique_taxa <- unique(taxa_vec[valid_pool])
  par_taxa <- length(unique_taxa) >= 2L && n >= 80L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")
  eval_one_taxon <- function(tx) {
    .eval_taxon_anchor_pruning_and_promo(
      tx = tx, df = df, X_use = X_use, S_pool = S_pool, S_pool_inv = S_pool_inv,
      type_info = type_info, is_gs = is_gs, taxa_vec = taxa_vec,
      is_type_valid = is_type_valid, is_nb_type_valid = is_nb_type_valid,
      is_het_syn_non_type = is_het_syn_non_type, is_gold = is_gold,
      is_silver = is_silver, is_dis_mode = is_dis_mode, eff_sol = eff_sol,
      par_taxa = par_taxa, mahal_q_gold_silver = mahal_q_gold_silver,
      mahal_q_bronze_promo = mahal_q_bronze_promo, min_anc = min_anc,
      anchor_pool = anchor_pool, det_tier = det_tier, lqw_col = lqw_col
    )
  }

  tx_results <- if (isTRUE(par_taxa)) {
    packera_mclapply(unique_taxa, eval_one_taxon, mc.cores = resolve_packera_workers(64L, length(unique_taxa)))
  } else {
    lapply(unique_taxa, eval_one_taxon)
  }
  all_pruned <- unlist(lapply(tx_results, `[[`, "pruned"), use.names = FALSE)
  all_pruned_het_type <- unlist(lapply(tx_results, `[[`, "pruned_het_type"), use.names = FALSE)
  all_promo <- unlist(lapply(tx_results, `[[`, "promo"), use.names = FALSE)
  passive_het_type_mask <- rep(FALSE, n)
  if (length(all_pruned_het_type) > 0L) {
    passive_het_type_mask[all_pruned_het_type] <- TRUE
  }
  if (length(all_pruned) > 0L) {
    all_pruned <- setdiff(all_pruned, setdiff(which(is_nb_type_valid), all_pruned_het_type))
    if (length(all_pruned) > 0L) {
      is_gs[all_pruned] <- FALSE
      is_gold[all_pruned] <- FALSE
      is_silver[all_pruned] <- FALSE
    }
  }
  if (length(all_promo) > 0L) {
    is_promoted_bronze[all_promo] <- TRUE
  }

  active_mask <- is_gs | is_promoted_bronze

  # Step 3: Confident Learning Out-of-Fold Anchor Pruning (Northcutt et al. 2021)
  is_confident_label_noise <- rep(FALSE, n)
  bad_cl_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  cl_act_idx <- which(active_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_cl_tx))
  taxa_levels <- sort(unique(taxa_vec[cl_act_idx]))
  if (length(taxa_levels) >= 2L && !is.null(X_use) && ncol(X_use) >= 1L) {
    K_cl <- length(taxa_levels)
    p_cl <- min(ncol(X_use), max(1L, floor(length(cl_act_idx) / 3L)), 10L)
    X_cl <- X_use[, seq_len(p_cl), drop = FALSE]
    idx_by_tx <- lapply(taxa_levels, function(tx_k) which(active_mask & (taxa_vec == tx_k)))
    names(idx_by_tx) <- taxa_levels
    if (sum(vapply(idx_by_tx, length, integer(1L)) >= 2L) >= 2L) {
      S_cl <- if (!is.null(S_pool)) {
        S_pool[seq_len(p_cl), seq_len(p_cl), drop = FALSE]
      } else {
        stats::cov(X_cl[cl_act_idx, , drop = FALSE]) + diag(1e-4, nrow = p_cl, ncol = p_cl)
      }
      S_cl <- 0.75 * S_cl + 0.25 * diag(pmax(diag(S_cl), 1e-4), nrow = p_cl, ncol = p_cl)
      S_cl_inv <- tryCatch(solve(S_cl), error = function(e) diag(1.0 / pmax(diag(S_cl), 1e-4), nrow = p_cl, ncol = p_cl))
      mu_cl_mat <- matrix(0.0, nrow = K_cl, ncol = p_cl, dimnames = list(taxa_levels, NULL))
      for (k in seq_len(K_cl)) {
        ik <- idx_by_tx[[k]]
        if (length(ik) > 0L) {
          gold_ik_strict <- ik[is_gold[ik] & !is_het_syn_non_type[ik]]
          gold_ik <- if (length(gold_ik_strict) >= 3L) gold_ik_strict else ik[is_gold[ik]]
          ref_ik <- if (length(gold_ik) >= 3L) gold_ik else ik
          mu_cl_mat[k, ] <- .compute_robust_taxon_medoid(X_cl[ref_ik, , drop = FALSE], S_inv = S_cl_inv)
        }
      }
      P_cl <- matrix(0.0, nrow = length(cl_act_idx), ncol = K_cl, dimnames = list(NULL, taxa_levels))
      y_cl_idx <- match(taxa_vec[cl_act_idx], taxa_levels)
      for (ii in seq_along(cl_act_idx)) {
        r_i <- cl_act_idx[ii]
        yk <- y_cl_idx[ii]
        mu_use <- mu_cl_mat
        ik <- idx_by_tx[[yk]]
        if (length(ik) >= 3L) {
          ik_loo <- setdiff(ik, r_i)
          gold_loo_strict <- ik_loo[is_gold[ik_loo] & !is_het_syn_non_type[ik_loo]]
          gold_loo <- if (length(gold_loo_strict) >= 3L) gold_loo_strict else ik_loo[is_gold[ik_loo]]
          ref_loo <- if (length(gold_loo) >= 3L) gold_loo else ik_loo
          mu_use[yk, ] <- .compute_robust_taxon_medoid(X_cl[ref_loo, , drop = FALSE], S_inv = S_cl_inv)
        }
        diff_i <- sweep(mu_use, 2L, X_cl[r_i, ], "-")
        d2_i <- rowSums((diff_i %*% S_cl_inv) * diff_i)
        ll_i <- -0.5 * d2_i / max(1.0, 0.5 * p_cl)
        ll_i <- ll_i - max(ll_i, na.rm = TRUE)
        p_i <- exp(ll_i)
        P_cl[ii, ] <- p_i / max(sum(p_i, na.rm = TRUE), 1e-12)
      }
      t_j <- numeric(K_cl)
      for (k in seq_len(K_cl)) {
        m_k <- which(y_cl_idx == k)
        t_j[k] <- if (length(m_k) > 0L) mean(P_cl[m_k, k], na.rm = TRUE) else 0.5
      }
      for (k in seq_len(K_cl)) {
        m_k <- which(y_cl_idx == k)
        if (length(m_k) <= 3L) next
        cand_ii <- integer(0)
        for (ii in m_k) {
          r_i <- cl_act_idx[ii]
          if (is_nb_type_valid[r_i]) next
          p_self <- P_cl[ii, k]
          oth_k <- setdiff(seq_len(K_cl), k)
          best_oth <- oth_k[which.max(P_cl[ii, oth_k])]
          p_oth <- P_cl[ii, best_oth]
          if (is.finite(p_self) && is.finite(p_oth) &&
              (p_self < 0.45 * t_j[k]) && (p_oth >= max(0.55, t_j[best_oth]))) {
            cand_ii <- c(cand_ii, ii)
          }
        }
        if (length(cand_ii) > 0L) {
          max_drop_k <- max(0L, length(m_k) - 3L)
          if (max_drop_k > 0L) {
            ord_drop <- cand_ii[order(P_cl[cand_ii, k], decreasing = FALSE)]
            drop_rows <- cl_act_idx[ord_drop[seq_len(min(length(ord_drop), max_drop_k))]]
            is_confident_label_noise[drop_rows] <- TRUE
          }
        }
      }
      if (any(is_confident_label_noise)) {
        active_mask[is_confident_label_noise] <- FALSE
        is_gold[is_confident_label_noise] <- FALSE
        is_silver[is_confident_label_noise] <- FALSE
        is_promoted_bronze[is_confident_label_noise] <- FALSE
      }
    }
  }

  if (any(extralimital_mask)) {
    active_mask[extralimital_mask] <- FALSE
    is_gold[extralimital_mask] <- FALSE
    is_silver[extralimital_mask] <- FALSE
    is_promoted_bronze[extralimital_mask] <- FALSE
  }
  if (any(is_open_hyb_non_type)) {
    active_mask[is_open_hyb_non_type] <- FALSE
    is_gold[is_open_hyb_non_type] <- FALSE
    is_silver[is_open_hyb_non_type] <- FALSE
    is_promoted_bronze[is_open_hyb_non_type] <- FALSE
  }

  anchor_tier_source <- ifelse(
    is_confident_label_noise, "Passive_Confident_Noise_Demoted",
    ifelse(
      passive_het_type_mask, "Passive_Heterotypic_Type_Audit",
      ifelse(
        extralimital_mask, "Passive_Extralimital_Audit",
        ifelse(
          is_open_hyb_non_type, "Passive_Open_or_Hybrid_Concept",
          ifelse(
            is_type_valid & active_mask, "Tier_0_Type_Anchor",
            ifelse(
              is_gold, "Tier_1_Gold_Anchor",
              ifelse(
                is_silver, "Tier_2_Silver_Anchor",
                ifelse(is_promoted_bronze, "Tier_3_Bronze_Promoted", "Passive_Projected")
              )
            )
          )
        )
      )
    )
  )

  list(
    active_mask = as.logical(active_mask),
    is_type_anchor = as.logical(is_type_valid & active_mask),
    is_gold_anchor = as.logical(is_gold),
    is_silver_anchor = as.logical(is_silver),
    is_promoted_bronze = as.logical(is_promoted_bronze),
    is_confident_label_noise = as.logical(is_confident_label_noise),
    anchor_tier_source = as.character(anchor_tier_source)
  )
}
