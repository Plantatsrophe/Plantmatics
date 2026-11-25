#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/plasticity_mitigation_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Strictly modality-independent phenotypic plasticity mitigation utilities for
#   Step 04 unsupervised GMM clustering and Canonical Discriminant Analysis:
#     1. Specimen-intrinsic capitulum phenophase covariate ingestion from Step 02
#        (`voucher_phenological_states.csv`) and within-taxon ontogenetic/fruiting
#        residualization free of geographic coordinates, collection dates, or
#        ecoregion labels (preserving Padial et al. 2010 modality independence).
#     2. Ray & Lindsay (2005) / Hennig (2010) multivariate Gaussian ridgeline
#        density unimodality diagnosis & skewed-tail component merging in pure
#        morphological space, allowing genuine G = 1 unimodal species nulls while
#        preserving true multimodal morphological discontinuities.
# ==============================================================================

merge_voucher_phenophase_covariates <- function(df, harmonics_path = NULL, vouchers_path = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  n <- nrow(df)
  in_paths <- stats::na.omit(c(
    if (!is.null(harmonics_path) && nzchar(harmonics_path)) harmonics_path else NA_character_,
    if (!is.null(vouchers_path) && nzchar(vouchers_path)) vouchers_path else NA_character_
  ))
  cand_dirs <- unique(dirname(in_paths))
  ph_df <- NULL
  for (d_cand in cand_dirs) {
    ph_csv <- file.path(d_cand, "voucher_phenological_states.csv")
    if (file.exists(ph_csv) && isTRUE(file.info(ph_csv)$size > 0)) {
      ph_df <- tryCatch(utils::read.csv(ph_csv, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(ph_df) && nrow(ph_df) > 0L) break else ph_df <- NULL
    }
  }
  if (is.null(ph_df) && exists("is_sqlite_active_for_inputs", mode = "function") &&
      isTRUE(is_sqlite_active_for_inputs(in_paths)) && exists("connect_packera_db", mode = "function")) {
    conn <- tryCatch(connect_packera_db(read_only = TRUE), error = function(e) NULL)
    if (!is.null(conn)) {
      on.exit(if (exists("disconnect_packera_db", mode = "function")) disconnect_packera_db(conn) else tryCatch(DBI::dbDisconnect(conn), error = function(e) NULL), add = TRUE)
      if (isTRUE(tryCatch(DBI::dbExistsTable(conn, "voucher_reproductive_phenology"), error = function(e) FALSE))) {
        avail_f <- tryCatch(DBI::dbListFields(conn, "voucher_reproductive_phenology"), error = function(e) character(0))
        want_f <- intersect(c("catalogNumber", "plant_individual_id", "n_anthesis", "n_bud", "n_fruiting", "n_senescent", "total_capitula", "corymb_head_count_est", "radiate_head_fraction", "phenological_state", "is_verified_anthesis"), avail_f)
        if ("catalogNumber" %in% want_f && length(want_f) >= 2L) {
          q_cols <- paste(vapply(want_f, function(cl) as.character(DBI::dbQuoteIdentifier(conn, cl)), character(1)), collapse = ", ")
          ph_df <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT %s FROM voucher_reproductive_phenology;", q_cols)), error = function(e) NULL)
        }
      }
    }
  }
  if (!is.null(ph_df) && nrow(ph_df) > 0L && "catalogNumber" %in% names(ph_df)) {
    ph_pid <- if ("plant_individual_id" %in% names(ph_df)) ifelse(is.na(ph_df$plant_individual_id), 0L, as.integer(ph_df$plant_individual_id)) else rep(0L, nrow(ph_df))
    df_pid <- if ("plant_individual_id" %in% names(df)) ifelse(is.na(df$plant_individual_id), 0L, as.integer(df$plant_individual_id)) else rep(0L, n)
    k_ph <- paste(trimws(as.character(ph_df$catalogNumber)), ph_pid, sep = "__")
    k_df <- paste(trimws(as.character(df$catalogNumber)), df_pid, sep = "__")
    n_pids_sheet <- tapply(ph_pid, trimws(as.character(ph_df$catalogNumber)), function(z) length(unique(z)))
    single_pid_sheet <- !(trimws(as.character(df$catalogNumber)) %in% names(n_pids_sheet[n_pids_sheet > 1L]))
    m_exact <- match(k_df, k_ph); m_sheet <- match(paste(trimws(as.character(df$catalogNumber)), 0L, sep = "__"), k_ph)
    m_use <- ifelse(!is.na(m_exact), m_exact, ifelse(single_pid_sheet, m_sheet, NA_integer_))
    ph_cols <- intersect(c("n_anthesis", "n_bud", "n_fruiting", "n_senescent", "total_capitula", "corymb_head_count_est", "radiate_head_fraction", "phenological_state", "is_verified_anthesis"), names(ph_df))
    for (pc_nm in ph_cols) {
      if (!(pc_nm %in% names(df))) df[[pc_nm]] <- ph_df[[pc_nm]][m_use]
      else { miss_m <- is.na(df[[pc_nm]]) | !nzchar(trimws(as.character(df[[pc_nm]]))); if (any(miss_m)) df[[pc_nm]][miss_m] <- ph_df[[pc_nm]][m_use[miss_m]] }
    }
  }
  df
}

residualize_ontogenetic_phenophase_plasticity <- function(
    df,
    target_cols,
    ref_mask = NULL,
    intraspecific_mode = FALSE,
    min_partial_r2 = 0.04
) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 12L) return(df)
  n <- nrow(df)
  repro_excl <- c("AR_inv", "AR_inv_std", "log1p_capitulum_count", "log1p_capitulum_count_std", "radiate_head_fraction", "radiate_head_fraction_std")
  target_cols <- setdiff(intersect(as.character(target_cols), names(df)), repro_excl)
  if (length(target_cols) == 0L) return(df)

  auth_w <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n)
  auth_w[!is.finite(auth_w) | auth_w < 0.0] <- 0.0
  sb_w <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, n)
  sb_w[!is.finite(sb_w) | sb_w <= 0.0] <- 1.0
  eff_w <- auth_w * sb_w

  if (is.null(ref_mask) || length(ref_mask) != n) {
    is_t1 <- if ("morph_tier_source" %in% names(df)) (!is.na(df$morph_tier_source) & df$morph_tier_source == "Basal_Rosette") else rep(TRUE, n)
    ref_mask <- (auth_w > 0.0) & is_t1
  }
  ref_idx <- which(as.logical(ref_mask))
  if (length(ref_idx) < 12L) return(df)

  # Strictly specimen-intrinsic phenophase covariate (zero geographic coordinates, DOY, ecoregion, or size backdoors)
  cov_list <- list()
  n_fr <- rep(0.0, n)
  has_fr <- FALSE
  for (fc in intersect(c("n_fruiting", "n_senescent"), names(df))) {
    vv <- suppressWarnings(as.numeric(df[[fc]]))
    if (any(is.finite(vv) & vv > 0.0)) {
      has_fr <- TRUE
      n_fr <- n_fr + ifelse(is.finite(vv) & vv > 0.0, vv, 0.0)
    }
  }
  tot_cap <- if ("total_capitula" %in% names(df)) suppressWarnings(as.numeric(df$total_capitula)) else if ("capitula_count" %in% names(df)) suppressWarnings(as.numeric(df$capitula_count)) else rep(NA_real_, n)
  if ("corymb_head_count_est" %in% names(df)) {
    corymb_cap <- suppressWarnings(as.numeric(df$corymb_head_count_est))
    tot_cap <- ifelse(is.finite(tot_cap) | is.finite(corymb_cap), pmax(tot_cap, corymb_cap, na.rm = TRUE), NA_real_)
  }
  if (has_fr) {
    has_cap_obs <- is.finite(tot_cap) & (tot_cap >= 1.0)
    fr_frac <- ifelse(has_cap_obs, pmin(1.0, pmax(0.0, n_fr / pmax(tot_cap, n_fr))), NA_real_)
    ok_fr_ref <- ref_idx[is.finite(fr_frac[ref_idx])]
    if (length(ok_fr_ref) >= 12L && stats::sd(fr_frac[ok_fr_ref], na.rm = TRUE) > 0.05) {
      cov_list[["fruiting_frac"]] <- fr_frac
    }
  }
  if (length(cov_list) == 0L) return(df)

  sp_raw <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Packera_sp", n)
  sp_raw[is.na(sp_raw) | !nzchar(trimws(sp_raw))] <- "Packera_sp"

  for (col_nm in target_cols) {
    yv <- suppressWarnings(as.numeric(df[[col_nm]]))
    if (sum(is.finite(yv[ref_idx])) < 12L) next
    is_pc_col <- grepl("^PC[0-9]+$", col_nm)
    if (is_pc_col && !isTRUE(intraspecific_mode)) next

    C_full <- do.call(cbind, cov_list)
    fin_sub <- is.finite(yv[ref_idx]) & stats::complete.cases(C_full[ref_idx, , drop = FALSE])
    if (sum(fin_sub) < 12L) next

    sub_rows <- ref_idx[fin_sub]
    sub_sp <- droplevels(factor(sp_raw[sub_rows]))
    sp_cnts <- table(sub_sp)
    keep_sp <- names(sp_cnts[sp_cnts >= 3L])
    ok_sp <- as.character(sub_sp) %in% keep_sp
    if (sum(ok_sp) < 12L) next

    sub_rows <- sub_rows[ok_sp]
    sub_sp <- droplevels(sub_sp[ok_sp])
    w_sub <- eff_w[sub_rows]
    if (sum(w_sub) <= 1e-12) w_sub <- rep(1.0, length(sub_rows))

    # Require homogeneous within-taxon phenophase slopes before pooling across taxa
    if (nlevels(sub_sp) >= 2L && sum(table(sub_sp) >= 5L) >= 2L) {
      sp_ge5 <- names(which(table(sub_sp) >= 5L))
      ii_ok <- as.character(sub_sp) %in% sp_ge5
      df_int <- data.frame(y = yv[sub_rows[ii_ok]], fruiting_frac = C_full[sub_rows[ii_ok], 1L], sub_sp = droplevels(sub_sp[ii_ok]), w = w_sub[ii_ok])
      an_int <- tryCatch(stats::anova(stats::lm(y ~ fruiting_frac * sub_sp, data = df_int, weights = w)), error = function(e) NULL)
      p_int <- if (!is.null(an_int) && "fruiting_frac:sub_sp" %in% rownames(an_int)) suppressWarnings(as.numeric(an_int["fruiting_frac:sub_sp", "Pr(>F)"])) else NA_real_
      if (is.finite(p_int) && p_int <= 0.05) next
    }

    w_sp_mean <- function(v_in) {
      unname(vapply(split(seq_along(v_in), sub_sp), function(ii) stats::weighted.mean(v_in[ii], w = w_sub[ii], na.rm = TRUE), numeric(1))[as.character(sub_sp)])
    }
    y_dm <- as.numeric(yv[sub_rows] - w_sp_mean(yv[sub_rows]))
    ss_y <- sum(w_sub * (y_dm^2))
    if (!is.finite(ss_y) || ss_y <= 1e-8) next

    C_dm <- apply(C_full[sub_rows, , drop = FALSE], 2L, function(cv) as.numeric(cv - w_sp_mean(cv)))
    if (!is.matrix(C_dm)) C_dm <- matrix(C_dm, ncol = ncol(C_full))
    keep_c <- which(apply(C_dm, 2L, stats::sd) > 1e-5)
    if (length(keep_c) == 0L) next
    C_dm <- C_dm[, keep_c, drop = FALSE]
    C_eval <- C_full[, keep_c, drop = FALSE]

    fit_w <- tryCatch(stats::lm.wfit(x = C_dm, y = y_dm, w = w_sub), error = function(e) NULL)
    if (is.null(fit_w)) next
    beta_w <- fit_w$coefficients
    beta_w[!is.finite(beta_w)] <- 0.0
    rss_w <- sum(w_sub * (fit_w$residuals^2))
    r2_w <- max(0.0, 1.0 - (rss_w / ss_y))
    df1 <- max(1L, fit_w$rank)
    df2 <- length(sub_rows) - nlevels(sub_sp) - df1
    if (df2 < 5L) next
    f_stat <- ((ss_y - rss_w) / df1) / max(rss_w / df2, 1e-12)
    p_val <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)

    if (is.finite(r2_w) && r2_w >= min_partial_r2 && is.finite(p_val) && p_val < 0.05) {
      fin_all <- which(is.finite(yv) & stats::complete.cases(C_eval))
      if (length(fin_all) > 0L) {
        delta_vals <- as.numeric(C_eval[fin_all, , drop = FALSE] %*% beta_w)
        nz_adj <- which(is.finite(delta_vals) & abs(delta_vals) > 0.0)
        if (length(nz_adj) > 0L) {
          df[[col_nm]][fin_all[nz_adj]] <- round(yv[fin_all[nz_adj]] - delta_vals[nz_adj], 6)
        }
      }
    }
  }
  df
}

.compute_ridgeline_valley_ratio <- function(pi_j, mu_j, Sig_j, pi_l, mu_l, Sig_l, n_grid = 101L) {
  p <- length(mu_j)
  reg_j <- Sig_j + diag(1e-5, nrow = p, ncol = p)
  reg_l <- Sig_l + diag(1e-5, nrow = p, ncol = p)
  inv_j <- tryCatch(solve(reg_j), error = function(e) diag(x = 1.0 / pmax(diag(as.matrix(reg_j)), 1e-5), nrow = p, ncol = p))
  inv_l <- tryCatch(solve(reg_l), error = function(e) diag(x = 1.0 / pmax(diag(as.matrix(reg_l)), 1e-5), nrow = p, ncol = p))
  ldet_j <- as.numeric(determinant(reg_j, logarithm = TRUE)$modulus)
  ldet_l <- as.numeric(determinant(reg_l, logarithm = TRUE)$modulus)

  w_j <- pi_j / max(pi_j + pi_l, 1e-12)
  w_l <- 1.0 - w_j
  alphas <- seq(0.0, 1.0, length.out = n_grid)
  dens_v <- numeric(n_grid)

  for (idx_a in seq_along(alphas)) {
    al <- alphas[idx_a]
    A_al <- (1.0 - al) * inv_j + al * inv_l
    rhs <- (1.0 - al) * (inv_j %*% mu_j) + al * (inv_l %*% mu_l)
    x_al <- tryCatch(as.numeric(solve(A_al, rhs)), error = function(e) (1.0 - al) * mu_j + al * mu_l)
    dj <- x_al - mu_j
    dl <- x_al - mu_l
    log_fj <- -0.5 * (p * log(2.0 * pi) + ldet_j + sum(dj * as.numeric(inv_j %*% dj)))
    log_fl <- -0.5 * (p * log(2.0 * pi) + ldet_l + sum(dl * as.numeric(inv_l %*% dl)))
    m_log <- max(log_fj, log_fl)
    dens_v[idx_a] <- exp(m_log) * (w_j * exp(log_fj - m_log) + w_l * exp(log_fl - m_log))
  }

  is_peak <- logical(n_grid)
  is_peak[1L] <- dens_v[1L] >= dens_v[2L]
  is_peak[n_grid] <- dens_v[n_grid] >= dens_v[n_grid - 1L]
  if (n_grid >= 3L) {
    for (ii in 2L:(n_grid - 1L)) {
      is_peak[ii] <- (dens_v[ii] >= dens_v[ii - 1L]) && (dens_v[ii] >= dens_v[ii + 1L]) &&
        ((dens_v[ii] > dens_v[ii - 1L]) || (dens_v[ii] > dens_v[ii + 1L]))
    }
  }
  peak_idx <- which(is_peak)
  if (length(peak_idx) <= 1L) {
    return(list(n_modes = 1L, ridge_ratio = 1.0))
  }
  p_left <- peak_idx[1L]
  p_right <- peak_idx[length(peak_idx)]
  if (p_right - p_left <= 1L) {
    return(list(n_modes = 1L, ridge_ratio = 1.0))
  }
  f_saddle <- min(dens_v[p_left:p_right])
  f_small_peak <- min(dens_v[p_left], dens_v[p_right])
  ratio <- if (f_small_peak > 1e-300) min(1.0, f_saddle / f_small_peak) else 1.0
  list(n_modes = length(peak_idx), ridge_ratio = as.numeric(ratio))
}

merge_plastic_gmm_components_ridgeline <- function(
    X_eval, X_core, probs_eval, probs_core,
    intraspecific_mode = FALSE, ridge_intraspecific_max = 0.55, ridge_interspecific_max = 0.90
) {
  if (is.null(probs_eval) || !is.matrix(probs_eval) || ncol(probs_eval) < 2L) {
    return(list(best_k = 1L, classification = rep(1L, nrow(X_eval)), uncertainty = rep(0.0, nrow(X_eval)),
                probs_eval = matrix(1.0, nrow = nrow(X_eval), ncol = 1L), merged_any = FALSE))
  }
  Z_core <- as.matrix(probs_core); Z_eval <- as.matrix(probs_eval); p <- ncol(X_core); merged_any <- FALSE
  ridge_thr <- if (isTRUE(intraspecific_mode)) as.numeric(ridge_intraspecific_max) else as.numeric(ridge_interspecific_max)
  nk_0 <- colSums(Z_core, na.rm = TRUE); valid_0 <- which(nk_0 > 1e-6)
  if (length(valid_0) < ncol(Z_core)) {
    if (length(valid_0) == 0L) valid_0 <- 1L
    Z_core <- sweep(Z_core[, valid_0, drop = FALSE], 1L, pmax(rowSums(Z_core[, valid_0, drop = FALSE]), 1e-12), "/")
    Z_eval <- sweep(Z_eval[, valid_0, drop = FALSE], 1L, pmax(rowSums(Z_eval[, valid_0, drop = FALSE]), 1e-12), "/")
    merged_any <- TRUE; nk_0 <- colSums(Z_core, na.rm = TRUE)
  }
  G_0 <- ncol(Z_core)
  if (G_0 < 2L) {
    return(list(best_k = 1L, classification = rep(1L, nrow(Z_eval)), uncertainty = rep(0.0, nrow(Z_eval)),
                probs_eval = matrix(1.0, nrow = nrow(Z_eval), ncol = 1L), merged_any = merged_any))
  }
  # Precompute pairwise ridgelines on initial components so merging never chains inflated pooled covariances.
  pis_0 <- nk_0 / sum(nk_0)
  mus_0 <- lapply(seq_len(G_0), function(k) colSums(X_core * Z_core[, k]) / max(sum(Z_core[, k]), 1e-12))
  Sigs_0 <- lapply(seq_len(G_0), function(k) {
    wk <- Z_core[, k]; Xc_k <- sweep(X_core, 2L, mus_0[[k]], "-"); Sk <- crossprod(Xc_k * sqrt(wk)) / max(sum(wk), 1.0)
    Sk + diag(0.02 * mean(diag(as.matrix(Sk))) + 1e-4, nrow = p, ncol = p)
  })
  nm_0 <- matrix(2L, nrow = G_0, ncol = G_0); rr_0 <- matrix(0.0, nrow = G_0, ncol = G_0)
  for (j in seq_len(G_0 - 1L)) for (l in (j + 1L):G_0) {
    rl <- .compute_ridgeline_valley_ratio(pis_0[j], mus_0[[j]], Sigs_0[[j]], pis_0[l], mus_0[[l]], Sigs_0[[l]])
    nm_0[j, l] <- nm_0[l, j] <- as.integer(rl$n_modes); rr_0[j, l] <- rr_0[l, j] <- as.numeric(rl$ridge_ratio)
  }
  grps <- as.list(seq_len(G_0))
  repeat {
    if (length(grps) < 2L) break
    best_pair <- NULL; best_ratio <- -1.0
    for (a in seq_len(length(grps) - 1L)) for (b in (a + 1L):length(grps)) {
      nm_sub <- nm_0[grps[[a]], grps[[b]], drop = FALSE]; rr_sub <- rr_0[grps[[a]], grps[[b]], drop = FALSE]
      pair_rr <- max(rr_sub)
      if ((any(nm_sub <= 1L) || (pair_rr >= ridge_thr)) && (pair_rr > best_ratio)) {
        best_ratio <- pair_rr; best_pair <- c(a, b)
      }
    }
    if (is.null(best_pair)) break
    grps[[best_pair[1L]]] <- c(grps[[best_pair[1L]]], grps[[best_pair[2L]]]); grps[[best_pair[2L]]] <- NULL
    merged_any <- TRUE
  }
  Z_out <- do.call(cbind, lapply(grps, function(idx) rowSums(Z_eval[, idx, drop = FALSE])))
  best_k <- ncol(Z_out)
  cls_out <- if (best_k == 1L) rep(1L, nrow(Z_out)) else as.integer(max.col(Z_out, ties.method = "first"))
  unc_out <- if (best_k == 1L) rep(0.0, nrow(Z_out)) else round(pmax(0.0, 1.0 - apply(Z_out, 1L, max)), 6)
  list(best_k = as.integer(best_k), classification = cls_out, uncertainty = unc_out, probs_eval = Z_out, merged_any = merged_any)
}

# Empirical Bayes institutional location-offset shrinkage accounting for biological taxon & ecoregion structure
apply_eb_institution_batch_centering <- function(df, target_cols, min_inst_n = 10L, allow_shape_centering = FALSE) {
  if (!is.data.frame(df) || nrow(df) < 10L || !("institutionCode" %in% names(df))) return(df)
  cand_cols <- intersect(unique(as.character(target_cols)), names(df))
  if (length(cand_cols) == 0L) return(df)
  inst_all <- trimws(as.character(df$institutionCode))
  inst_all[is.na(df$institutionCode) | !nzchar(inst_all)] <- NA_character_
  cent_all <- if ("is_centroid_or_inst_coord" %in% names(df)) {
    cv_a <- df$is_centroid_or_inst_coord
    if (is.logical(cv_a)) !is.na(cv_a) & cv_a else if (is.numeric(cv_a)) !is.na(cv_a) & cv_a != 0 else tolower(trimws(as.character(cv_a))) %in% c("true", "1", "t", "yes")
  } else rep(FALSE, nrow(df))
  ref_mask <- !is.na(inst_all)
  if ("authority_weight" %in% names(df)) ref_mask <- ref_mask & !is.na(df$authority_weight) & (suppressWarnings(as.numeric(df$authority_weight)) > 0)
  if ("is_passive" %in% names(df)) ref_mask <- ref_mask & (is.na(df$is_passive) | !as.logical(df$is_passive))
  if ("morph_tier_source" %in% names(df)) ref_mask <- ref_mask & (is.na(df$morph_tier_source) | !(df$morph_tier_source %in% c("Cauline_Only_Fallback", "Transitional_Passive")))
  if (sum(ref_mask, na.rm = TRUE) < 10L || length(unique(inst_all[ref_mask])) < 2L) ref_mask <- !is.na(inst_all)
  if (sum(ref_mask, na.rm = TRUE) < 10L || length(unique(inst_all[ref_mask])) < 2L) return(df)

  is_regionally_confounded <- FALSE
  inst_ref <- inst_all[ref_mask]; cent_ref <- cent_all[ref_mask]; phot_reg_col <- NULL
  for (g_col in intersect(c("regional_group", "epa_ecoregion_l3"), names(df))) {
    rg_chk <- trimws(as.character(df[[g_col]][ref_mask]))
    ok_ph <- !cent_ref & !is.na(inst_ref) & !is.na(df[[g_col]][ref_mask]) & nzchar(rg_chk) & !(tolower(rg_chk) %in% c("na", "unknown", "unknown_ecoregion"))
    if (length(unique(rg_chk[ok_ph])) >= 2L) { phot_reg_col <- g_col; break }
  }
  for (g_col in intersect(c("regional_group", "epa_ecoregion_l3", "stateProvince"), names(df))) {
    reg_raw <- trimws(as.character(df[[g_col]][ref_mask]))
    ok_rg <- !cent_ref & !is.na(inst_ref) & !is.na(df[[g_col]][ref_mask]) & nzchar(reg_raw)
    if (length(unique(reg_raw[ok_rg])) >= 2L) {
      tab_rg <- table(inst_ref[ok_rg], reg_raw[ok_rg])
      tab_rg <- tab_rg[rowSums(tab_rg) >= 2L, , drop = FALSE]; tab_rg <- tab_rg[, colSums(tab_rg) >= 1L, drop = FALSE]
      if (nrow(tab_rg) >= 2L && ncol(tab_rg) >= 2L) {
        chi_v <- suppressWarnings(stats::chisq.test(tab_rg, simulate.p.value = FALSE)$statistic)
        cramers_v <- sqrt(as.numeric(chi_v) / max(1.0, sum(tab_rg) * (min(nrow(tab_rg), ncol(tab_rg)) - 1L)))
        dom_share <- stats::weighted.mean(apply(tab_rg, 1L, max) / pmax(rowSums(tab_rg), 1), w = rowSums(tab_rg))
        is_regionally_confounded <- (is.finite(cramers_v) && cramers_v >= 0.45) || (is.finite(dom_share) && dom_share >= 0.75)
      }
      break
    }
  }

  sp_all <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else if ("scientificName" %in% names(df)) as.character(df$scientificName) else rep("Unassigned", nrow(df))
  sp_all[is.na(sp_all) | !nzchar(trimws(sp_all))] <- "Unassigned"
  w_all <- (if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, nrow(df))) *
    (if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, nrow(df)))
  w_all[!is.finite(w_all) | w_all <= 0] <- 1.0

  pre_eta <- attr(df, "pre_eb_batch_partial_eta_sq")
  if (is.null(pre_eta) || !is.numeric(pre_eta)) pre_eta <- stats::setNames(numeric(0), character(0))

  for (col in cand_cols) {
    yv <- suppressWarnings(as.numeric(df[[col]]))
    ok_all <- is.finite(yv) & !is.na(inst_all)
    r_idx <- which(ref_mask & is.finite(yv))
    if (length(r_idx) < 10L || stats::var(yv[r_idx]) <= 1e-12) next
    r_inst <- inst_all[r_idx]; r_sp <- sp_all[r_idx]; r_y <- yv[r_idx]; r_w <- w_all[r_idx]
    u_inst <- unique(r_inst); if (length(u_inst) < 2L) next

    # Compute exact pre-EB marginal ANCOVA partial eta-squared matching audit_institutional_batch_effects
    fit_pre <- tryCatch({
      sub_p <- data.frame(y = r_y, sp = droplevels(factor(r_sp)), inst = droplevels(factor(r_inst)))
      fml_p <- if (nlevels(sub_p$sp) >= 2L) stats::as.formula("y ~ sp + inst") else stats::as.formula("y ~ inst")
      stats::anova(stats::lm(fml_p, data = sub_p, weights = r_w))
    }, error = function(e) NULL)
    if (!is.null(fit_pre) && "inst" %in% rownames(fit_pre) && "Residuals" %in% rownames(fit_pre)) {
      ii <- which(rownames(fit_pre) == "inst")[1L]; ri <- which(rownames(fit_pre) == "Residuals")[1L]
      pre_eta[[col]] <- as.numeric(fit_pre$`Sum Sq`[ii] / max(fit_pre$`Sum Sq`[ii] + fit_pre$`Sum Sq`[ri], 1e-12))
    }
    is_phot <- grepl("^(blade_cie_|basal_zone_cie_l_delta|abaxial_cie_l_sheet_delta|blade_.*browning)", col)
    if (!is_phot && !(isTRUE(allow_shape_centering) && !isTRUE(is_regionally_confounded))) next
    r_bio <- r_sp
    if (is_phot && !is.null(phot_reg_col)) {
      rg_r <- trimws(as.character(df[[phot_reg_col]][r_idx]))
      ok_st <- !cent_all[r_idx] & !is.na(df[[phot_reg_col]][r_idx]) & nzchar(rg_r) & !(tolower(rg_r) %in% c("na", "unknown", "unknown_ecoregion"))
      st_lbl <- paste(r_sp[ok_st], rg_r[ok_st], sep = "||")
      min_st_n <- max(4L, as.integer(ceiling(0.5 * min_inst_n)))
      sh_st <- names(Filter(function(v) sum(table(v) >= min_st_n) >= 2L, split(r_inst[ok_st], st_lbl)))
      if (length(sh_st) > 0L) {
        keep_s <- which(ok_st)[st_lbl %in% sh_st]
        if (length(keep_s) >= 6L && length(unique(r_inst[keep_s])) >= 2L) {
          r_inst <- r_inst[keep_s]; r_bio <- paste(r_sp[keep_s], rg_r[keep_s], sep = "||"); r_y <- r_y[keep_s]; r_w <- r_w[keep_s]; u_inst <- unique(r_inst)
        } else if (isTRUE(is_regionally_confounded)) next
      } else if (isTRUE(is_regionally_confounded)) next
    } else if (is_phot && isTRUE(is_regionally_confounded)) next

    # Estimate marginal stratum/species means and zero-mean institutional location shifts iteratively
    d_raw <- stats::setNames(rep(0.0, length(u_inst)), u_inst)
    mu_sp <- stats::setNames(rep(mean(r_y), length(unique(r_bio))), unique(r_bio))
    for (it in seq_len(6L)) {
      for (sp_u in names(mu_sp)) {
        si <- which(r_bio == sp_u)
        mu_sp[sp_u] <- stats::weighted.mean(r_y[si] - d_raw[r_inst[si]], w = r_w[si])
      }
      e_it <- r_y - mu_sp[r_bio]
      for (ih in u_inst) {
        hi <- which(r_inst == ih)
        d_raw[ih] <- if (length(hi) >= 2L) stats::weighted.mean(e_it[hi], w = r_w[hi]) else 0.0
      }
      w_mean_d <- sum(vapply(u_inst, function(ih) sum(r_w[r_inst == ih]) * d_raw[ih], numeric(1))) / max(sum(r_w), 1e-12)
      d_raw <- d_raw - w_mean_d
    }

    res_w <- r_y - mu_sp[r_bio] - d_raw[r_inst]
    s2_in <- sum(r_w * res_w^2) / max(1.0, sum(r_w) - length(unique(r_bio)) - length(u_inst))
    n_by_i <- vapply(u_inst, function(ih) sum(r_inst == ih), integer(1))
    elig_v <- which(n_by_i >= min_inst_n)
    if (length(elig_v) < 2L) elig_v <- which(n_by_i >= 2L)
    if (length(elig_v) < 2L || !is.finite(s2_in) || s2_in <= 0) next
    se2_v <- s2_in / pmax(as.numeric(n_by_i[elig_v]), 1.0)
    tau2 <- max(0.0, stats::var(d_raw[elig_v]) - mean(se2_v, na.rm = TRUE))
    if (!is.finite(tau2) || tau2 <= 1e-9) next

    delta_eb <- stats::setNames(numeric(length(u_inst)), u_inst)
    for (ih in u_inst) {
      nh <- n_by_i[ih]
      if (nh >= 2L) {
        eb_w <- tau2 / (tau2 + (s2_in / nh))
        if (nh < min_inst_n) eb_w <- eb_w * (nh / min_inst_n)
        delta_eb[ih] <- d_raw[ih] * eb_w
      }
    }
    w_eb_mean <- sum(vapply(u_inst, function(ih) sum(r_w[r_inst == ih]) * delta_eb[ih], numeric(1))) / max(sum(r_w), 1e-12)
    delta_eb <- delta_eb - w_eb_mean
    adj_mask <- ok_all & (inst_all %in% names(delta_eb))
    yv[adj_mask] <- round(yv[adj_mask] - delta_eb[inst_all[adj_mask]], 5)
    df[[col]] <- yv
  }
  attr(df, "pre_eb_batch_partial_eta_sq") <- pre_eta
  df
}
