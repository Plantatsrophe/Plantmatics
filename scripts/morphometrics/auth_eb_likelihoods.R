if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_eb_likelihoods.R
# Determiner normalization, LOO PC likelihoods, and specimen discordance signals
# ==============================================================================

.normalize_determiner_id_vec <- function(df) {
  n <- nrow(df)
  raw_det <- rep("", n)
  invalid_raw <- c("na", "n/a", "none", "unknown", "null", "anonymous", "anon", "unspecified", "not recorded", "s.n.", "sn", "indet", "indet.")
  for (col in intersect(c("identifiedBy"), names(df))) {
    v <- as.character(df[[col]])
    v[is.na(v)] <- ""
    v <- trimws(v)
    need <- !nzchar(raw_det) & nzchar(v) & !(tolower(v) %in% invalid_raw)
    if (any(need)) raw_det[need] <- v[need]
  }
  out <- rep(NA_character_, n)
  has_raw <- nzchar(raw_det)
  if (!any(has_raw)) return(out)

  raw_sub <- raw_det[has_raw]
  u_raw <- unique(raw_sub)
  s <- tolower(u_raw)
  for (iter in 1:2) {
    cand <- sub(
      "^(det\\.?|determined[[:space:]]+by|annot\\.?|annotated[[:space:]]+by|ver\\.?|verified[[:space:]]+by|rev\\.?|revised[[:space:]]+by|conf\\.?|confirmed[[:space:]]+by|identified[[:space:]]+by)[:;,[:space:]]+",
      "",
      s
    )
    cand_norm <- trimws(gsub("[[:punct:][:space:]]+", " ", cand))
    keep_strip <- nchar(cand_norm) >= 3L
    if (any(keep_strip)) s[keep_strip] <- cand[keep_strip]
  }
  s <- gsub("(^|[^A-Za-z0-9_])(1[789][0-9]{2}|20[0-9]{2})([^A-Za-z0-9_]|$)", "\\1 \\3", s)
  s <- gsub("(^|[^A-Za-z0-9_])(1[789][0-9]{2}|20[0-9]{2})([^A-Za-z0-9_]|$)", "\\1 \\3", s)
  s <- gsub("[[:punct:][:space:]]+", " ", s)
  s <- trimws(s)
  invalid_norm <- c("", "na", "n a", "none", "unknown", "null", "anonymous", "anon", "unspecified", "not recorded", "s n", "indet", "unidentified")
  s[s %in% invalid_norm | nchar(s) < 2L] <- NA_character_
  out[has_raw] <- s[match(raw_sub, u_raw)]
  out
}

.compute_loo_pc_likelihoods <- function(df, taxa_vec, det_tier, is_promoted, eval_mask, valid_taxa, prior_vec) {
  n <- nrow(df)
  G <- length(valid_taxa)
  if (G < 2L) return(NULL)

  use_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
  if (length(use_cols) > 0L) {
    use_cols <- use_cols[vapply(use_cols, function(col) sum(is.finite(suppressWarnings(as.numeric(df[[col]])))) >= 3L, logical(1))]
  }
  if (length(use_cols) == 0L) return(NULL)

  X_pc <- do.call(cbind, lapply(use_cols, function(col) suppressWarnings(as.numeric(df[[col]]))))
  colnames(X_pc) <- use_cols
  finite_pc <- apply(X_pc, 1L, function(r) all(is.finite(r)))
  eval_pc <- eval_mask & finite_pc & (taxa_vec %in% valid_taxa)
  if (sum(eval_pc) < 3L) return(NULL)
  if (ncol(X_pc) > 2L) {
    p_cap <- max(2L, min(ncol(X_pc), floor(max(sum(eval_pc), 4L) / 4L), 15L))
    X_pc <- X_pc[, seq_len(p_cap), drop = FALSE]
  }
  p <- ncol(X_pc)

  is_high_tier <- (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted
  ref_idx_by_k <- vector("list", G)
  names(ref_idx_by_k) <- valid_taxa
  taxa_with_ref <- logical(G)

  for (k in seq_len(G)) {
    tx <- valid_taxa[k]
    cand_hi <- which(eval_pc & (taxa_vec == tx) & is_high_tier)
    cand_all <- which(eval_pc & (taxa_vec == tx))
    cand_k <- if (length(cand_hi) >= 2L) cand_hi else cand_all
    if (length(cand_k) >= 5L) {
      Xk <- X_pc[cand_k, , drop = FALSE]
      med_k <- apply(Xk, 2L, stats::median)
      mad_k <- pmax(apply(Xk, 2L, stats::mad), 1e-3)
      r_d2 <- rowSums(sweep(sweep(Xk, 2L, med_k, "-"), 2L, mad_k, "/")^2)
      inlier_k <- r_d2 <= stats::qchisq(0.975, df = p)
      min_keep <- max(2L, ceiling(0.55 * length(cand_k)))
      if (sum(inlier_k) >= min_keep) {
        cand_k <- cand_k[inlier_k]
      } else {
        cand_k <- cand_k[order(r_d2)[seq_len(min_keep)]]
      }
    }
    ref_idx_by_k[[k]] <- cand_k
    taxa_with_ref[k] <- length(cand_k) >= 1L
  }
  if (sum(taxa_with_ref) < 2L) return(NULL)

  sum_w_k <- numeric(G)
  sum_wx_k <- matrix(0.0, nrow = G, ncol = p)
  in_ref_mat <- matrix(FALSE, nrow = n, ncol = G)
  res_list <- list()
  n_ref_total <- 0L
  w_ref_total <- 0.0

  for (k in which(taxa_with_ref)) {
    idx_k <- ref_idx_by_k[[k]]
    uk <- pmax(0.05, prior_vec[idx_k])
    Xk <- X_pc[idx_k, , drop = FALSE]
    sw <- sum(uk)
    swx <- colSums(Xk * uk)
    sum_w_k[k] <- sw
    sum_wx_k[k, ] <- swx
    in_ref_mat[idx_k, k] <- TRUE
    mk <- swx / sw
    diff_k <- sweep(Xk, 2L, mk, "-")
    res_list[[length(res_list) + 1L]] <- diff_k * sqrt(uk)
    n_ref_total <- n_ref_total + length(idx_k)
    w_ref_total <- w_ref_total + sw
  }

  X_res <- do.call(rbind, res_list)
  denom_w <- max(w_ref_total - sum(taxa_with_ref), 1.0)
  S_w <- crossprod(X_res) / denom_w
  if (sum(diag(S_w)) <= 1e-10) {
    S_w <- S_w + diag(1e-3, p)
  }
  S_reg <- compute_ledoit_wolf_cov(S_w, X_centered = X_res, n_obs = n_ref_total, min_lambda = 1e-4, max_lambda = 0.25)
  chol_R <- tryCatch(base::chol(S_reg), error = function(e) base::chol(S_reg + diag(1e-4, p)))

  obs_idx <- which(eval_mask & finite_pc)
  if (length(obs_idx) == 0L) return(NULL)

  d2_mat <- matrix(Inf, nrow = length(obs_idx), ncol = G)
  X_obs <- X_pc[obs_idx, , drop = FALSE]
  u_obs <- pmax(0.05, prior_vec[obs_idx])

  for (k in which(taxa_with_ref)) {
    mu_mat <- matrix(sum_wx_k[k, ], nrow = length(obs_idx), ncol = p, byrow = TRUE) / sum_w_k[k]
    is_in_k <- in_ref_mat[obs_idx, k]
    if (any(is_in_k) && length(ref_idx_by_k[[k]]) >= 2L) {
      rem_w <- sum_w_k[k] - u_obs[is_in_k]
      can_loo <- rem_w > 1e-8
      if (any(can_loo)) {
        loo_rows <- which(is_in_k)[can_loo]
        num_loo <- sweep(-X_obs[loo_rows, , drop = FALSE] * u_obs[loo_rows], 2L, sum_wx_k[k, ], "+")
        mu_mat[loo_rows, ] <- num_loo / rem_w[can_loo]
      }
    }
    diff_obs <- X_obs - mu_mat
    whiten_diff <- matrix(backsolve(chol_R, t(diff_obs), transpose = TRUE), nrow = p, ncol = nrow(diff_obs))
    d2_mat[, k] <- colSums(whiten_diff * whiten_diff)
  }

  temp_scale <- max(1.0, 0.5 * p)
  log_lik <- -0.5 * (d2_mat / temp_scale)
  row_max <- apply(log_lik, 1L, max)
  lik_sub <- exp(sweep(log_lik, 1L, row_max, "-"))
  lik_sub <- sweep(lik_sub, 1L, pmax(rowSums(lik_sub), 1e-15), "/")

  L_pc <- matrix(NA_real_, nrow = n, ncol = G, dimnames = list(NULL, valid_taxa))
  L_pc[obs_idx, ] <- lik_sub
  L_pc
}

.extract_specimen_discordance_signals <- function(
  df,
  eval_mask = NULL,
  taxon_col = "scientificName",
  pc_cols = NULL,
  target_modality = c("all", "morph", "vision", "ecological")
) {
  target_modality <- match.arg(target_modality, c("all", "morph", "vision", "ecological"))
  n <- nrow(df)
  if (is.null(eval_mask) || is.character(eval_mask) || length(eval_mask) != n) {
    eval_mask <- rep(TRUE, n)
  }
  err_cols <- intersect(
    c("dawid_skene_error_rate", "label_error_rate", "error_rate_k", "c_error", "loo_discordance", "loo_error_rate"),
    names(df)
  )
  eps_err <- rep(NA_real_, n)
  if (length(err_cols) > 0L) {
    for (ec in err_cols) {
      need <- eval_mask & !is.finite(eps_err)
      if (!any(need)) break
      v <- suppressWarnings(as.numeric(df[[ec]]))
      ok <- need & is.finite(v)
      if (any(ok)) eps_err[ok] <- pmax(0.0, pmin(1.0, v[ok]))
    }
  }

  apply_bool_flag <- function(cur_eps, flag_vec, pos_rate = 0.80, neg_rate = 0.05) {
    fv <- as.logical(flag_vec)
    ok_f <- eval_mask & !is.na(fv)
    if (!any(ok_f)) return(cur_eps)
    only_f <- ok_f & !is.finite(cur_eps)
    if (any(only_f)) {
      cur_eps[only_f] <- ifelse(fv[only_f], pos_rate, neg_rate)
    }
    both_pos <- ok_f & is.finite(cur_eps) & fv
    if (any(both_pos)) {
      cur_eps[both_pos] <- pmax(cur_eps[both_pos], 1.0 - (1.0 - cur_eps[both_pos]) * 0.35)
    }
    cur_eps
  }

  if ("exsiccatae_taxonomic_conflict" %in% names(df)) {
    ex_conf <- as.logical(df$exsiccatae_taxonomic_conflict)
    if ("exsiccatae_gold_taxon" %in% names(df)) {
      g_tx <- trimws(as.character(df$exsiccatae_gold_taxon))
      s_tx <- .resolve_authority_taxon_vec(df)
      if (exists("standardize_packera_taxon", mode = "function")) {
        has_g <- !is.na(g_tx) & nzchar(g_tx) & !(tolower(g_tx) %in% c("na", "none", "unknown"))
        if (any(has_g)) g_tx[has_g] <- as.character(standardize_packera_taxon(g_tx[has_g]))
      }
      dt_ex <- .standardize_determiner_tier_vec(df)
      agrees_gold <- !is.na(ex_conf) & ex_conf & !is.na(g_tx) & nzchar(g_tx) & (g_tx == s_tx)
      # Specialist sheet itself suffers zero conflict penalty; elevated Silver duplicate gets mild 0.18 uncertainty
      ex_disagree <- ifelse(agrees_gold, FALSE, ex_conf)
      eps_err <- apply_bool_flag(eps_err, ex_disagree, pos_rate = 0.85, neg_rate = 0.05)
      elev_silv <- eval_mask & agrees_gold & (dt_ex == "Tier_2_Silver")
      if (any(elev_silv)) {
        eps_err[elev_silv] <- ifelse(is.finite(eps_err[elev_silv]), pmax(eps_err[elev_silv], 0.18), 0.18)
      }
    } else {
      eps_err <- apply_bool_flag(eps_err, ex_conf, pos_rate = 0.85, neg_rate = 0.05)
    }
  }
  if (target_modality != "vision" && "is_label_corrupted" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$is_label_corrupted, pos_rate = 0.85, neg_rate = 0.05)
  }
  if (!(target_modality %in% c("morph", "ecological")) && "morphometric_outlier_flag" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$morphometric_outlier_flag, pos_rate = 0.75, neg_rate = 0.05)
  }
  if (target_modality != "ecological" && "multimodal_conflict_flag" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$multimodal_conflict_flag, pos_rate = 0.78, neg_rate = 0.05)
  }
  eps_err
}
