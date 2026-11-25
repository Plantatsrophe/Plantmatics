if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_anchor_pool_helpers.R
# Per-taxon Mahalanobis anchor pruning and adaptive Bronze promotion helper
# ==============================================================================

.eval_taxon_anchor_pruning_and_promo <- function(
  tx, df, X_use, S_pool, S_pool_inv, type_info,
  is_gs, taxa_vec, is_type_valid, is_nb_type_valid, is_het_syn_non_type,
  is_gold, is_silver, is_dis_mode, eff_sol, par_taxa,
  mahal_q_gold_silver, mahal_q_bronze_promo, min_anc,
  anchor_pool, det_tier, lqw_col
) {
  pruned_out <- integer(0)
  pruned_het_type_out <- integer(0)
  promo_out <- integer(0)
  idx_gs <- which(is_gs & (taxa_vec == tx))
  idx_gs_clean <- idx_gs
  gs_type_local <- is_type_valid[idx_gs]
  gs_nb_type_local <- is_nb_type_valid[idx_gs]
  gs_paratype_local <- gs_type_local & !gs_nb_type_local
  gs_bas_local <- .canonicalize_type_basionym_vec(df[idx_gs, , drop = FALSE], rep(tx, length(idx_gs)))
  gs_homo_type_local <- gs_type_local & .is_homotypic_sensu_stricto_basionym(
    gs_bas_local, rep(tx, length(idx_gs)), df = df[idx_gs, , drop = FALSE]
  )
  gs_het_type_local <- gs_type_local & !gs_homo_type_local
  # Prefer homotypic sensu stricto name-bearing types so divergent heterotypic types never pull a species' type-medoid target.
  gs_anchor_type_local <- if (any(gs_nb_type_local & gs_homo_type_local)) {
    gs_nb_type_local & gs_homo_type_local
  } else if (any(gs_nb_type_local)) {
    gs_nb_type_local
  } else {
    gs_type_local
  }
  gs_het_local <- is_het_syn_non_type[idx_gs]
  gs_gold_local <- is_gold[idx_gs]
  gs_silver_local <- is_silver[idx_gs]
  n_gold_k <- sum(gs_gold_local)
  dis_gold_n <- sum(is_dis_mode[idx_gs[gs_gold_local]], na.rm = TRUE)
  dis_gold_frac <- if (n_gold_k > 0L) (dis_gold_n / n_gold_k) else if (length(idx_gs) > 0L) mean(is_dis_mode[idx_gs], na.rm = TRUE) else 0.0
  if (!is.finite(dis_gold_frac)) dis_gold_frac <- 0.0
  gold_sol_pool <- if (n_gold_k > 0L) eff_sol[idx_gs[gs_gold_local]] else eff_sol[idx_gs]
  dis_Mode_pool <- if (n_gold_k > 0L) is_dis_mode[idx_gs[gs_gold_local]] else is_dis_mode[idx_gs]
  ent_sol_vals <- gold_sol_pool[!dis_Mode_pool & is.finite(gold_sol_pool)]
  dis_sol_vals <- gold_sol_pool[dis_Mode_pool & is.finite(gold_sol_pool)]
  ent_med_sol <- if (length(ent_sol_vals) > 0L) stats::median(ent_sol_vals, na.rm = TRUE) else NA_real_
  dis_med_sol <- if (length(dis_sol_vals) > 0L) stats::median(dis_sol_vals, na.rm = TRUE) else NA_real_
  has_bimodal_dis_cluster <- (dis_gold_n >= 5L) && is.finite(dis_med_sol) &&
    ((dis_med_sol <= 0.65) || (is.finite(ent_med_sol) && ((ent_med_sol - dis_med_sol) >= 0.18)))
  min_dis_gold_frac <- if (has_bimodal_dis_cluster) 0.15 else 0.25
  allow_dis_split <- isTRUE(if (n_gold_k > 0L) (dis_gold_frac >= min_dis_gold_frac) else (dis_gold_frac >= 0.30))
  idx_type_k <- idx_gs[gs_anchor_type_local]
  n_type_k <- length(idx_type_k)
  n_eff_type_k <- if (n_type_k > 0L) .count_effective_type_gatherings(df, idx_type_k) else 0.0
  gamma_k <- if (n_type_k > 0L) .compute_type_shrinkage_gamma(n_eff_type_k, length(idx_gs)) else 0.0

  # Step 1: Within-taxon robust MCD / medoid Mahalanobis pruning of Gold + Silver (chi-sq_0.95, p)
  # Two-stage Gold/Silver seeding & heteroblasty gate: when a taxon has >= 5 Gold+Type anchors,
  # estimate initial reference scatter on Gold+Type first and screen Silver rows at min(mahal_q_gold_silver, 0.90).
  if (!is.null(X_use) && length(idx_gs) >= 5L) {
    p_use <- min(ncol(X_use), length(idx_gs) - 2L)
    if (p_use >= 1L) {
      X_gs_raw <- X_use[idx_gs, seq_len(p_use), drop = FALSE]
      v_ok <- apply(X_gs_raw, 2, stats::var, na.rm = TRUE) > 1e-10
      if (any(v_ok)) {
        col_sub <- which(v_ok)
        X_gs <- X_gs_raw[, col_sub, drop = FALSE]
        p_eff <- ncol(X_gs)
        extreme_type_local <- rep(FALSE, nrow(X_gs))
        if (any(gs_anchor_type_local) && sum(!gs_anchor_type_local) >= 5L) {
          nt_strict <- gs_gold_local & !gs_type_local & !gs_het_local
          nt_mask <- if (sum(nt_strict) >= 5L) nt_strict else if (sum(gs_gold_local & !gs_type_local) >= 5L) (gs_gold_local & !gs_type_local) else !gs_anchor_type_local
          X_nt <- X_gs[nt_mask, , drop = FALSE]
          mu_nt <- .compute_robust_taxon_medoid(X_nt)
          cov_nt_raw <- stats::cov(X_nt)
          cov_nt_raw[!is.finite(cov_nt_raw)] <- 0.0
          reg_nt <- max(1e-4, 0.05 * max(mean(diag(cov_nt_raw), na.rm = TRUE), 1e-4))
          S_nt <- cov_nt_raw + diag(reg_nt, nrow = p_eff, ncol = p_eff)
          type_pos <- which(gs_anchor_type_local)
          d2_type_prelim <- tryCatch(
            stats::mahalanobis(X_gs[type_pos, , drop = FALSE], center = mu_nt, cov = S_nt),
            error = function(e) rep(0.0, length(type_pos))
          )
          extreme_type_local[type_pos] <- is.finite(d2_type_prelim) & (d2_type_prelim > stats::qchisq(0.99, df = p_eff))
        }
        non_ext_type_pos <- which(gs_anchor_type_local & !extreme_type_local)
        if (n_type_k > 0L && length(non_ext_type_pos) == 0L && sum(!gs_anchor_type_local) >= 5L) {
          mu_type_gs <- NULL
          gamma_k <- 0.0
        } else if (length(non_ext_type_pos) > 0L) {
          mu_type_gs <- .compute_robust_taxon_medoid(X_gs[non_ext_type_pos, , drop = FALSE])
        } else {
          mu_type_gs <- NULL
        }
        gs_dis_pos <- if (allow_dis_split) which(is_dis_mode[idx_gs]) else integer(0)
        gs_ent_pos <- if (allow_dis_split) which(!is_dis_mode[idx_gs]) else seq_len(nrow(X_gs))
        S_sub_pool <- if (!is.null(S_pool)) S_pool[col_sub, col_sub, drop = FALSE] else (stats::cov(X_gs) + diag(1e-4, nrow = p_eff, ncol = p_eff))
        d2_gs <- rep(NA_real_, nrow(X_gs))
        if (allow_dis_split && length(gs_dis_pos) >= 3L && length(gs_ent_pos) >= 3L) {
          drop_gs <- rep(FALSE, nrow(X_gs))
          for (pos_m in list(gs_ent_pos, gs_dis_pos)) {
            Xm <- X_gs[pos_m, , drop = FALSE]
            # Prefer non-heterotypic Gold+Type anchors for mode medoid when >= 3 exist.
            pos_m_gold_strict <- pos_m[gs_gold_local[pos_m] & !extreme_type_local[pos_m] & !gs_het_local[pos_m] & !gs_het_type_local[pos_m]]
            pos_m_gold <- if (length(pos_m_gold_strict) >= 3L) pos_m_gold_strict else pos_m[gs_gold_local[pos_m] & !extreme_type_local[pos_m]]
            has_gold_med_m <- length(pos_m_gold) >= 3L
            use_gold_m <- length(pos_m_gold) >= max(5L, p_eff + 3L)
            pos_m_fit <- if (use_gold_m) pos_m_gold else pos_m[!extreme_type_local[pos_m]]
            if (length(pos_m_fit) < 3L) pos_m_fit <- pos_m
            Xm_fit <- X_gs[pos_m_fit, , drop = FALSE]
            rob_m <- NULL; q_m <- nrow(Xm_fit)
            if (nrow(Xm_fit) >= max(12L, 5L * p_eff) && requireNamespace("MASS", quietly = TRUE)) {
              q_m <- min(nrow(Xm_fit) - 1L, max(floor(0.80 * nrow(Xm_fit)), floor((nrow(Xm_fit) + p_eff + 1L) / 2L)))
              rob_m <- .fit_fast_mcd_cov(Xm_fit, q_used = q_m, allow_parallel = !par_taxa)
            }
            cent_raw_m <- if (has_gold_med_m) {
              if (!is.null(rob_m) && use_gold_m) rob_m$center else .compute_robust_taxon_medoid(X_gs[pos_m_gold, , drop = FALSE])
            } else if (!is.null(rob_m)) {
              rob_m$center
            } else {
              .compute_robust_taxon_medoid(Xm_fit)
            }
            cent_m <- cent_raw_m
            cov_shift_m <- matrix(0.0, nrow = p_eff, ncol = p_eff)
            if (!is.null(mu_type_gs) && gamma_k > 0.0) {
              sub_types <- pos_m[gs_anchor_type_local[pos_m] & !extreme_type_local[pos_m]]
              mu_sub <- if (length(sub_types) > 0L) .compute_robust_taxon_medoid(X_gs[sub_types, , drop = FALSE]) else mu_type_gs
              cent_m <- (1.0 - gamma_k) * cent_raw_m + gamma_k * mu_sub
              d_mu_m <- cent_m - cent_raw_m
              cov_shift_m <- tcrossprod(d_mu_m) + diag(sum(d_mu_m^2) / max(1.0, p_eff), nrow = p_eff, ncol = p_eff)
            }
            cov_m <- (if (!is.null(rob_m)) {
              .mcd_consistency_factor(q_m / nrow(Xm_fit), p_eff) * rob_m$cov + diag(1e-5, nrow = p_eff, ncol = p_eff)
            } else if (nrow(Xm_fit) >= (p_eff + 3L)) {
              0.65 * stats::cov(Xm_fit) + 0.35 * S_sub_pool + diag(1e-5, nrow = p_eff, ncol = p_eff)
            } else {
              S_sub_pool + diag(1e-4, nrow = p_eff, ncol = p_eff)
            }) + cov_shift_m
            d2_m <- tryCatch(stats::mahalanobis(Xm, center = cent_m, cov = cov_m), error = function(e) rep(0.0, nrow(Xm)))
            d2_gs[pos_m] <- d2_m
            # Screen Paratypes strictly alongside standard Gold anchors while keeping name-bearing types immune.
            q_thr_m <- ifelse(
              gs_paratype_local[pos_m],
              stats::qchisq(mahal_q_gold_silver, df = p_eff),
              ifelse(
                gs_silver_local[pos_m] & has_gold_med_m,
                stats::qchisq(min(mahal_q_gold_silver, 0.90), df = p_eff),
                stats::qchisq(mahal_q_gold_silver, df = p_eff)
              )
            )
            dr_m <- is.finite(d2_m) & (d2_m > q_thr_m)
            if (any(dr_m) && sum(!dr_m) >= 2L) drop_gs[pos_m[dr_m]] <- TRUE
          }
        } else {
          # Prefer non-heterotypic Gold+Type anchors for overall taxon medoid when >= 3 exist.
          gold_core_strict <- which(gs_gold_local & !extreme_type_local & !gs_het_local & !gs_het_type_local)
          gold_core_pos <- if (length(gold_core_strict) >= 3L) gold_core_strict else which(gs_gold_local & !extreme_type_local)
          has_gold_medoid <- length(gold_core_pos) >= 3L
          use_gold_core <- length(gold_core_pos) >= max(5L, p_eff + 3L)
          X_gs_fit <- if (use_gold_core) {
            X_gs[gold_core_pos, , drop = FALSE]
          } else if (any(extreme_type_local) && sum(!extreme_type_local) >= 5L) {
            X_gs[!extreme_type_local, , drop = FALSE]
          } else {
            X_gs
          }
          rob_gs <- NULL; q_mcd <- nrow(X_gs_fit)
          if (nrow(X_gs_fit) >= max(12L, 5L * p_eff) && requireNamespace("MASS", quietly = TRUE)) {
            q_mcd <- min(nrow(X_gs_fit) - 1L, max(floor(0.80 * nrow(X_gs_fit)), floor((nrow(X_gs_fit) + p_eff + 1L) / 2L)))
            rob_gs <- .fit_fast_mcd_cov(X_gs_fit, q_used = q_mcd, allow_parallel = !par_taxa)
          }
          center_raw_gs <- if (has_gold_medoid) {
            if (!is.null(rob_gs) && use_gold_core) rob_gs$center else .compute_robust_taxon_medoid(X_gs[gold_core_pos, , drop = FALSE])
          } else if (!is.null(rob_gs)) {
            rob_gs$center
          } else {
            .compute_robust_taxon_medoid(X_gs_fit)
          }
          center_gs <- center_raw_gs
          cov_shift_gs <- matrix(0.0, nrow = p_eff, ncol = p_eff)
          if (!is.null(mu_type_gs) && gamma_k > 0.0) {
            center_gs <- (1.0 - gamma_k) * center_raw_gs + gamma_k * mu_type_gs
            d_mu_gs <- center_gs - center_raw_gs
            cov_shift_gs <- tcrossprod(d_mu_gs) + diag(sum(d_mu_gs^2) / max(1.0, p_eff), nrow = p_eff, ncol = p_eff)
          }
          cov_gs <- (if (!is.null(rob_gs)) {
            .mcd_consistency_factor(q_mcd / nrow(X_gs_fit), p_eff) * rob_gs$cov + diag(1e-5, nrow = p_eff, ncol = p_eff)
          } else if (nrow(X_gs_fit) >= (p_eff + 3L)) {
            0.75 * stats::cov(X_gs_fit) + 0.25 * S_sub_pool + diag(1e-4, nrow = p_eff, ncol = p_eff)
          } else {
            0.50 * stats::cov(X_gs_fit) + 0.50 * S_sub_pool + diag(1e-4, nrow = p_eff, ncol = p_eff)
          }) + cov_shift_gs
          d2_gs <- tryCatch(stats::mahalanobis(X_gs, center = center_gs, cov = cov_gs), error = function(e) rep(0.0, nrow(X_gs)))
          q_thr_gs <- ifelse(
            gs_paratype_local,
            stats::qchisq(mahal_q_gold_silver, df = p_eff),
            ifelse(
              gs_silver_local & has_gold_medoid,
              stats::qchisq(min(mahal_q_gold_silver, 0.90), df = p_eff),
              stats::qchisq(mahal_q_gold_silver, df = p_eff)
            )
          )
          drop_gs <- is.finite(d2_gs) & (d2_gs > q_thr_gs)
        }
        # Nomenclatural type pruning immunity: preserve homotypic name-bearing/non-paratype types unconditionally;
        # allow pruning strongly divergent heterotypic types (> chi-sq_0.995) when >= 5 sensu stricto anchors exist.
        is_explicit_paratype_local <- grepl("paratype", type_info$type_status_clean[idx_gs], ignore.case = TRUE)
        is_homo_type_immune <- gs_homo_type_local & !is_explicit_paratype_local
        if (any(is_homo_type_immune)) {
          drop_gs[is_homo_type_immune] <- FALSE
        }
        if (any(gs_het_type_local)) {
          n_ss_core <- sum((gs_gold_local & !gs_het_local & !gs_het_type_local) | gs_homo_type_local)
          het_prune <- gs_het_type_local & (n_ss_core >= 5L) & is.finite(d2_gs) & (d2_gs > stats::qchisq(0.995, df = p_eff))
          drop_gs[gs_het_type_local] <- het_prune[gs_het_type_local]
        }
        if (any(drop_gs) && sum(!drop_gs) >= 3L) {
          pruned_out <- idx_gs[drop_gs]
          pruned_het_type_out <- idx_gs[drop_gs & gs_het_type_local]
          idx_gs_clean <- idx_gs[!drop_gs]
        }
      }
    }
  }

  # Step 2: Adaptive Anchor Promotion for Gold/Silver-poor taxa (proportional eff_min_anc cap)
  n_gs <- length(idx_gs_clean)
  n_gold_clean <- sum(is_gold[idx_gs_clean])
  n_tx_pool <- sum(anchor_pool & (taxa_vec == tx))
  eff_min_anc <- if (n_gs == 0L) {
    min(min_anc, n_tx_pool)
  } else if (n_gs < 2L) {
    min(2L, n_tx_pool)
  } else if (n_tx_pool >= min_anc) {
    base_anc <- min(min_anc, as.integer(floor(0.5 * n_tx_pool)))
    if (n_gold_clean == 0L) max(base_anc, min(n_gs + 3L, as.integer(floor(0.6 * n_tx_pool)))) else base_anc
  } else {
    n_gs
  }
  if (tx != "Unassigned_Packera" && n_gs < eff_min_anc) {
    idx_br <- which((det_tier == "Tier_3_Bronze") & anchor_pool & (taxa_vec == tx))
    n_needed <- eff_min_anc - n_gs
    if (length(idx_br) > 0L && n_needed > 0L) {
      if (is.null(X_use) || ncol(X_use) == 0L) {
        ord_q <- if (!is.null(lqw_col)) {
          order(lqw_col[idx_br], decreasing = TRUE, na.last = TRUE)
        } else {
          seq_along(idx_br)
        }
        promo_out <- idx_br[ord_q[seq_len(min(length(idx_br), n_needed))]]
      } else {
        p_eff <- ncol(X_use)
        ref_ids <- if (n_gs >= 1L) idx_gs_clean else idx_br
        X_ref <- X_use[ref_ids, , drop = FALSE]
        medoid_tx <- .compute_robust_taxon_medoid(X_ref, S_inv = S_pool_inv)
        if (n_type_k > 0L && gamma_k > 0.0) {
          mu_type_tx <- .compute_robust_taxon_medoid(X_use[idx_gs[gs_anchor_type_local], , drop = FALSE], S_inv = S_pool_inv)
          medoid_tx <- (1.0 - gamma_k) * medoid_tx + gamma_k * mu_type_tx
        }
        base_cov_tx <- if (n_gs >= (p_eff + 3L)) {
          0.75 * stats::cov(X_use[idx_gs_clean, , drop = FALSE]) + 0.25 * S_pool + diag(1e-5, nrow = p_eff, ncol = p_eff)
        } else if (!is.null(S_pool)) {
          S_pool
        } else {
          diag(pmax(apply(X_use, 2, stats::var, na.rm = TRUE), 1e-4), nrow = p_eff, ncol = p_eff)
        }
        pred_scale <- 1.0 + (1.0 / max(1.0, nrow(X_ref)))
        cov_tx <- pred_scale * base_cov_tx
        d2_br <- tryCatch(
          stats::mahalanobis(X_use[idx_br, , drop = FALSE], center = medoid_tx, cov = cov_tx),
          error = function(e) rowSums(sweep(X_use[idx_br, , drop = FALSE], 2, medoid_tx, "-")^2)
        )
        ref_dis <- if (allow_dis_split) ref_ids[is_dis_mode[ref_ids]] else integer(0)
        ref_ent <- if (allow_dis_split) ref_ids[!is_dis_mode[ref_ids]] else ref_ids
        if (allow_dis_split && length(ref_dis) >= 2L && length(ref_ent) >= 2L) {
          med_dis <- .compute_robust_taxon_medoid(X_use[ref_dis, , drop = FALSE], S_inv = S_pool_inv)
          med_ent <- .compute_robust_taxon_medoid(X_use[ref_ent, , drop = FALSE], S_inv = S_pool_inv)
          d2_dis <- tryCatch(stats::mahalanobis(X_use[idx_br, , drop = FALSE], center = med_dis, cov = cov_tx), error = function(e) d2_br)
          d2_ent <- tryCatch(stats::mahalanobis(X_use[idx_br, , drop = FALSE], center = med_ent, cov = cov_tx), error = function(e) d2_br)
          d2_br <- pmin(d2_br, d2_dis, d2_ent, na.rm = TRUE)
        }
        ord_br <- order(d2_br, na.last = NA)
        q80_thr <- stats::qchisq(mahal_q_bronze_promo, df = p_eff)
        elig_ord <- ord_br[d2_br[ord_br] <= q80_thr]
        promo_ids <- idx_br[elig_ord[seq_len(min(length(elig_ord), n_needed))]]

        # Minimum-viability safeguard for 0/1-anchor rare taxa: allow up to 2 total anchors if within chi-sq_0.95
        min_viable <- min(2L, eff_min_anc)
        if ((n_gs + length(promo_ids)) < min_viable && length(ord_br) > length(promo_ids)) {
          q95_thr <- stats::qchisq(0.95, df = p_eff)
          rem_ord <- setdiff(ord_br[d2_br[ord_br] <= q95_thr], elig_ord)
          if (length(rem_ord) > 0L) {
            n_extra <- min(length(rem_ord), min_viable - (n_gs + length(promo_ids)))
            promo_ids <- c(promo_ids, idx_br[rem_ord[seq_len(n_extra)]])
          }
        }
        if (length(promo_ids) > 0L) promo_out <- promo_ids
      }
    }
  }
  list(pruned = pruned_out, pruned_het_type = pruned_het_type_out, promo = promo_out)
}
