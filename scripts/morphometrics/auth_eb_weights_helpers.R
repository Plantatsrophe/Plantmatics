if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_eb_weights_helpers.R
# Empirical Bayes tier/determiner blending and Dawid-Skene EM helpers
# ==============================================================================

.eb_blend_tier_and_determiner <- function(
  conc_vec, obs_mask, pi_vec, eval_mask, tier_cell, alpha0, taxa_vec, det_id
) {
  safe_logit <- function(p) {
    pc <- pmax(1e-5, pmin(1.0 - 1e-5, p))
    log(pc / (1.0 - pc))
  }
  theta_tier <- pi_vec
  for (tc in unique(tier_cell[eval_mask])) {
    idx_t_all <- which(eval_mask & (tier_cell == tc))
    idx_t_obs <- idx_t_all[obs_mask[idx_t_all]]
    pi_t <- pi_vec[idx_t_all[1L]]
    n_t <- length(idx_t_obs)
    s_t <- if (n_t > 0L) sum(conc_vec[idx_t_obs]) else 0.0
    theta_tier[idx_t_all] <- (alpha0 * pi_t + s_t) / (alpha0 + n_t)
  }

  theta_taxon_tier <- pi_vec
  taxon_tier_cell <- paste0(taxa_vec, "||", tier_cell)
  for (ck in unique(taxon_tier_cell[eval_mask])) {
    idx_all <- which(eval_mask & (taxon_tier_cell == ck))
    idx_obs <- idx_all[obs_mask[idx_all]]
    pi_t <- pi_vec[idx_all[1L]]
    n_c <- length(idx_obs)
    s_c <- if (n_c > 0L) sum(conc_vec[idx_obs]) else 0.0
    theta_taxon_tier[idx_all] <- (alpha0 * pi_t + s_c) / (alpha0 + n_c)
  }

  theta_cell_blend <- 0.50 * theta_tier + 0.50 * theta_taxon_tier
  theta_blend <- theta_cell_blend
  known_dets <- unique(det_id[eval_mask & !is.na(det_id)])
  if (length(known_dets) > 0L) {
    for (d_key in known_dets) {
      idx_d_all <- which(eval_mask & !is.na(det_id) & (det_id == d_key))
      idx_d_obs <- idx_d_all[obs_mask[idx_d_all]]
      n_d <- length(idx_d_obs)
      if (n_d == 0L) next
      pi_bar_d <- mean(pi_vec[idx_d_all])
      theta_d <- (alpha0 * pi_bar_d + sum(conc_vec[idx_d_obs])) / (alpha0 + n_d)
      if (length(unique(pi_vec[idx_d_all])) > 1L) {
        delta_d <- safe_logit(theta_d) - safe_logit(pi_bar_d)
        theta_d_i <- stats::plogis(safe_logit(pi_vec[idx_d_all]) + delta_d)
      } else {
        theta_d_i <- rep(theta_d, length(idx_d_all))
      }
      lambda_d <- n_d / (0.5 * alpha0 + n_d)
      theta_blend[idx_d_all] <- (1.0 - lambda_d) * theta_cell_blend[idx_d_all] + lambda_d * theta_d_i
    }
  }
  pmax(0.05, pmin(0.995, theta_blend))
}

.eb_build_dawid_skene_likelihood_matrix <- function(
  n, G, valid_taxa, y_idx, eval_mask, extralimital_ds,
  eps_err, pred_sources, bad_tx, L_pc, det_tier
) {
  has_ds_Row <- rep(FALSE, n)
  L_mat <- NULL

  if (G >= 2L) {
    L_mat <- matrix(1.0 / G, nrow = n, ncol = G, dimnames = list(NULL, valid_taxa))
    ds_rows <- which(eval_mask & !is.na(y_idx) & !extralimital_ds)
    if (length(ds_rows) == 0L) ds_rows <- which(eval_mask & !is.na(y_idx))
    for (i in ds_rows) {
      yi <- y_idx[i]
      log_odds_vec <- numeric(0)
      wt_vec <- numeric(0)
      off_accum <- rep(0.0, G)
      off_count <- 0.0

      if (is.finite(eps_err[i])) {
        p_self_e <- max(0.005, min(0.995, 1.0 - eps_err[i]))
        lo_e <- log(p_self_e / (1.0 - p_self_e))
        log_odds_vec <- c(log_odds_vec, lo_e)
        wt_vec <- c(wt_vec, if (lo_e < 0.0) 2.2 else 1.4)
      }

      has_exog_disc_i <- (is.finite(eps_err[i]) && eps_err[i] > 0.35)
      if (!has_exog_disc_i && length(pred_sources) > 0L) {
        for (ps_chk in pred_sources) {
          if (isTRUE(ps_chk$valid[i]) && !is.na(ps_chk$taxon[i]) && !(ps_chk$taxon[i] %in% bad_tx) && (ps_chk$taxon[i] != valid_taxa[yi])) {
            has_exog_disc_i <- TRUE
            break
          }
        }
      }

      if (!is.null(L_pc) && all(is.finite(L_pc[i, ]))) {
        p_self_pc <- max(0.002, min(0.998, L_pc[i, yi]))
        odds_pc <- (G - 1.0) * p_self_pc / max(1.0 - p_self_pc, 1e-6)
        lo_pc <- log(max(odds_pc, 1e-4))
        # Guard peripheral Tier_1_Gold ecotypes against circular morphology-only down-weighting when exogenous signals do not conflict
        if (!has_exog_disc_i && (det_tier[i] == "Tier_1_Gold") && (lo_pc < 0.0)) {
          lo_pc_use <- max(lo_pc, -0.30)
          wt_pc_use <- 0.45
        } else {
          lo_pc_use <- lo_pc
          wt_pc_use <- if (lo_pc < 0.0) 1.5 else 1.0
        }
        log_odds_vec <- c(log_odds_vec, lo_pc_use)
        wt_vec <- c(wt_vec, wt_pc_use)
        off_v <- L_pc[i, ]
        off_v[yi] <- 0.0
        if (sum(off_v) > 1e-12) {
          off_accum <- off_accum + (off_v / sum(off_v))
          off_count <- off_count + 1.0
        }
      }

      if (length(pred_sources) > 0L) {
        for (ps in pred_sources) {
          if (!isTRUE(ps$valid[i])) next
          ptx <- ps$taxon[i]
          if (is.na(ptx) || (ptx %in% bad_tx)) next
          c_m <- ps$conf[i]
          if (ptx == valid_taxa[yi]) {
            lo_m <- log(c_m / max(1.0 - c_m, 1e-4))
            log_odds_vec <- c(log_odds_vec, lo_m)
            wt_vec <- c(wt_vec, 1.1)
          } else {
            lo_m <- log(max(1.0 - c_m, 0.02) / c_m)
            log_odds_vec <- c(log_odds_vec, lo_m)
            wt_vec <- c(wt_vec, 1.8)
            pk_idx <- match(ptx, valid_taxa)
            if (!is.na(pk_idx) && pk_idx != yi) {
              v_m <- rep((1.0 - c_m) / max(1L, G - 1L), G)
              v_m[yi] <- 0.0
              v_m[pk_idx] <- c_m
              if (sum(v_m) > 1e-12) {
                off_accum <- off_accum + (v_m / sum(v_m))
                off_count <- off_count + 1.2
              }
            }
          }
        }
      }

      if (length(log_odds_vec) > 0L) {
        has_ds_Row[i] <- TRUE
        comb_odds <- exp(sum(wt_vec * log_odds_vec) / sum(wt_vec))
        B_i <- if (comb_odds <= 1.0) {
          max(0.005, comb_odds)
        } else {
          1.0 + 0.75 * log(comb_odds) / (log(comb_odds) + 1.5)
        }
        l_self <- B_i / (B_i + G - 1.0)
        if (off_count > 0.0 && sum(off_accum) > 1e-12) {
          q_off <- off_accum / sum(off_accum)
        } else {
          q_off <- rep(1.0 / (G - 1.0), G)
          q_off[yi] <- 0.0
        }
        L_mat[i, ] <- (1.0 - l_self) * q_off
        L_mat[i, yi] <- l_self
      }
    }
  }
  list(L_mat = L_mat, has_ds_Row = has_ds_Row)
}

.eb_run_dawid_skene_em_update <- function(
  w_eb, n, G, valid_taxa, y_idx, L_mat, has_ds_Row,
  pi_vec, tier_ceil, alpha0, eval_mask, eps_err,
  extralimital_ds, blend_tier_and_determiner
) {
  # 2-Step Dawid-Skene EM Update (Dawid & Skene 1979) with Beta/Dirichlet Empirical Bayes Shrinkage
  T_mat <- matrix(0.0, nrow = n, ncol = G, dimnames = list(NULL, valid_taxa))
  ds_obs_idx <- which(has_ds_Row)
  for (i in ds_obs_idx) {
    yi <- y_idx[i]
    l_yi <- L_mat[i, yi]
    B_init <- (G - 1.0) * l_yi / max(1.0 - l_yi, 1e-8)
    r0 <- (pi_vec[i] * B_init) / max(pi_vec[i] * B_init + (1.0 - pi_vec[i]), 1e-12)
    r0 <- max(0.005, min(0.995, r0))
    off_l <- L_mat[i, ]
    off_l[yi] <- 0.0
    off_s <- sum(off_l)
    if (off_s > 1e-12) {
      T_mat[i, ] <- (1.0 - r0) * (off_l / off_s)
    } else {
      T_mat[i, ] <- (1.0 - r0) / (G - 1.0)
    }
    T_mat[i, yi] <- r0
  }

  theta_curr <- pi_vec
  for (em_step in 1:2) {
    # M-step: determiner + (species_standardized, determiner_tier) Beta-Binomial shrinkage & confusion matrix Q
    t_self_prev <- numeric(n)
    t_self_prev[ds_obs_idx] <- T_mat[cbind(ds_obs_idx, y_idx[ds_obs_idx])]
    conc_m <- pmin(t_self_prev, tier_ceil)
    theta_curr <- blend_tier_and_determiner(conc_m, has_ds_Row)

    Q_off <- matrix(1.0 / (G - 1.0), nrow = G, ncol = G)
    diag(Q_off) <- 0.0
    if (G > 2L) {
      prior_q <- alpha0 / (G - 1.0)
      for (k in seq_len(G)) {
        counts_k <- rep(prior_q, G)
        counts_k[k] <- 0.0
        for (m in seq_len(G)) {
          if (m == k) next
          m_rows <- ds_obs_idx[y_idx[ds_obs_idx] == m]
          if (length(m_rows) > 0L) {
            counts_k[m] <- counts_k[m] + sum(T_mat[m_rows, k])
          }
        }
        Q_off[k, ] <- counts_k / max(sum(counts_k), 1e-12)
      }
    }

    # E-step: update latent true-taxon posteriors T_{i, k}^{(s)}
    for (i in ds_obs_idx) {
      yi <- y_idx[i]
      th_i <- theta_curr[i]
      u_vec <- (1.0 - th_i) * Q_off[, yi] * L_mat[i, ]
      u_vec[yi] <- th_i * L_mat[i, yi]
      u_sum <- sum(u_vec)
      if (is.finite(u_sum) && u_sum > 1e-15) {
        T_mat[i, ] <- u_vec / u_sum
      }
    }
  }

  t_self_final <- numeric(n)
  t_self_final[ds_obs_idx] <- T_mat[cbind(ds_obs_idx, y_idx[ds_obs_idx])]
  conc_final <- pmin(t_self_final, tier_ceil)
  w_eb[eval_mask] <- theta_curr[eval_mask]
  p_prior <- pmax(0.05, pmin(0.995, theta_curr[ds_obs_idx]))
  l_obs <- pmax(0.02, pmin(0.995, conc_final[ds_obs_idx]))
  l_alt <- pmax(0.02, 1.0 - l_obs)
  post_bayes <- (p_prior * l_obs) / pmax(1e-8, p_prior * l_obs + (1.0 - p_prior) * l_alt)
  w_eb[ds_obs_idx] <- pmax(0.05, pmin(0.995, 0.35 * p_prior + 0.65 * post_bayes))

  # Handle any evaluated rows with unknown taxon string that still carry explicit error rates
  rem_err <- which(eval_mask & !has_ds_Row & is.finite(eps_err) & !extralimital_ds)
  if (length(rem_err) > 0L) {
    fb_obs <- eval_mask & is.finite(eps_err) & !extralimital_ds
    fb_conc <- pmin(pmax(0.0, 1.0 - eps_err), tier_ceil)
    fb_theta <- blend_tier_and_determiner(fb_conc, fb_obs)
    p_pr_r <- pmax(0.05, pmin(0.995, fb_theta[rem_err]))
    l_ob_r <- pmax(0.02, pmin(0.995, fb_conc[rem_err]))
    l_al_r <- pmax(0.02, 1.0 - l_ob_r)
    pb_r <- (p_pr_r * l_ob_r) / pmax(1e-8, p_pr_r * l_ob_r + (1.0 - p_pr_r) * l_al_r)
    w_eb[rem_err] <- pmax(0.05, pmin(0.995, 0.35 * p_pr_r + 0.65 * pb_r))
  }
  w_eb
}
