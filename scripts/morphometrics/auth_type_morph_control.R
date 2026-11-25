if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_type_morph_control.R
# Nomenclatural type-specimen morphological control and heterotypic audit
# ==============================================================================

compute_type_specimen_morph_control <- function(
  df,
  X_mat = NULL,
  S_reg = NULL,
  cda_res = NULL,
  gmm_res = NULL,
  flags = NULL,
  triage = NULL,
  reasons = NULL,
  taxon_col = "species_standardized"
) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L) {
    return(list(
      type_status = character(0),
      is_type_specimen = logical(0),
      type_rank = integer(0),
      has_empirical_type_anchor = logical(0),
      type_morph_mahal_d2 = numeric(0),
      nearest_type_morph_taxon = character(0),
      nearest_type_basionym = character(0),
      nearest_type_basionym_d2 = numeric(0),
      type_morph_concordance = numeric(0),
      is_typical_morphotype = logical(0),
      min_congener_mahal_d2 = numeric(0),
      min_congener_taxon = character(0),
      open_set_novelty_p = numeric(0),
      open_set_df = integer(0),
      type_centroids = NULL,
      basionym_type_centroids = NULL,
      flags = flags,
      triage = triage,
      reasons = reasons
    ))
  }

  t_info <- .identify_type_specimen_vec(df)
  type_rank <- t_info$type_rank
  is_type_specimen <- t_info$is_type_specimen
  type_status <- if ("type_status" %in% names(df)) {
    ts <- trimws(as.character(df$type_status))
    ifelse(!is.na(ts) & nzchar(ts), ts, t_info$type_status_clean)
  } else {
    t_info$type_status_clean
  }

  if (is.null(X_mat) && !is.null(cda_res) && !is.null(cda_res$X_proj)) {
    X_mat <- cda_res$X_proj
  }
  if (is.null(X_mat)) {
    cand_cols <- if (!is.null(cda_res) && !is.null(cda_res$selected_features)) {
      intersect(as.character(cda_res$selected_features), names(df))
    } else {
      character(0)
    }
    if (length(cand_cols) == 0L) {
      cand_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
    }
    if (length(cand_cols) > 0L) {
      X_mat <- as.matrix(df[, cand_cols, drop = FALSE])
    } else if (!is.null(cda_res) && !is.null(cda_res$canonical_scores)) {
      X_mat <- as.matrix(cda_res$canonical_scores)
    } else {
      X_mat <- matrix(0.0, nrow = n, ncol = 2L)
    }
  }
  sh <- if (!is.null(cda_res) && is.list(cda_res) && !is.null(cda_res$shared_subspace) && is.list(cda_res$shared_subspace)) {
    cda_res$shared_subspace
  } else {
    NULL
  }
  no_efa_idx <- if (!is.null(sh) && !is.null(sh$no_efa_idx)) {
    v_ne <- as.integer(sh$no_efa_idx)
    unique(v_ne[!is.na(v_ne) & v_ne >= 1L & v_ne <= n])
  } else {
    integer(0)
  }
  has_efa_row <- rep(TRUE, n)
  if (length(no_efa_idx) > 0L) has_efa_row[no_efa_idx] <- FALSE

  X_mat <- as.matrix(X_mat)
  storage.mode(X_mat) <- "double"
  p <- ncol(X_mat)
  for (j in seq_len(p)) {
    bad_j <- !is.finite(X_mat[, j])
    if (any(bad_j)) {
      ok_j <- !bad_j & has_efa_row
      if (!any(ok_j)) ok_j <- !bad_j
      m_j <- mean(X_mat[ok_j, j], na.rm = TRUE)
      X_mat[bad_j, j] <- if (is.finite(m_j)) m_j else 0.0
    }
  }

  if (is.null(S_reg) && !is.null(cda_res) && !is.null(cda_res$S_reg)) {
    S_reg <- cda_res$S_reg
  }
  if (is.null(S_reg) || !is.matrix(S_reg) || nrow(S_reg) != p || ncol(S_reg) != p) {
    S_raw <- if (n >= 3L) stats::cov(X_mat) else diag(1.0, nrow = p, ncol = p)
    S_raw[!is.finite(S_raw)] <- 0.0
    S_reg <- compute_ledoit_wolf_cov(S_raw + diag(1e-5, nrow = p, ncol = p), X_centered = scale(X_mat, center = TRUE, scale = FALSE), n_obs = n)
  }
  S_inv <- tryCatch(
    solve(S_reg),
    error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-5), nrow = p, ncol = p)
  )

  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  if (all(taxa_vec == "Unassigned_Packera") && !is.null(cda_res) && !is.null(cda_res$grouping_col)) {
    taxa_vec <- as.character(cda_res$grouping_col)
  }
  det_tier <- .standardize_determiner_tier_vec(df)
  is_t1_pure <- .identify_tier1_pure_mask(df)
  if (sum(is_t1_pure) == 0L) is_t1_pure <- rep(TRUE, n)
  auth_w <- if ("authority_weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$authority_weight))
  } else if (!is.null(cda_res) && !is.null(cda_res$authority_weight)) {
    suppressWarnings(as.numeric(cda_res$authority_weight))
  } else {
    rep(1.0, n)
  }
  auth_w[!is.finite(auth_w) | auth_w <= 0.0] <- 1.0

  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  active_taxa <- sort(unique(taxa_vec[!is.na(taxa_vec) & !(taxa_vec %in% bad_tx)]))
  if (length(active_taxa) == 0L) {
    active_taxa <- sort(unique(stats::na.omit(taxa_vec)))
  }
  G <- length(active_taxa)

  gmm_cl <- if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n) {
    as.integer(gmm_res$classification)
  } else if ("gmm_cluster" %in% names(df)) {
    suppressWarnings(as.integer(df$gmm_cluster))
  } else {
    rep(1L, n)
  }
  u_clusters <- sort(unique(stats::na.omit(gmm_cl)))
  cl_centroids <- if (length(u_clusters) > 0L) {
    cent_mat <- matrix(0.0, nrow = length(u_clusters), ncol = p, dimnames = list(as.character(u_clusters), colnames(X_mat)))
    for (ci in seq_along(u_clusters)) {
      c_val <- u_clusters[ci]
      c_rows <- which(!is.na(gmm_cl) & (gmm_cl == c_val))
      c_rows_efa <- c_rows[has_efa_row[c_rows]]
      if (length(c_rows_efa) > 0L) c_rows <- c_rows_efa
      cent_mat[ci, ] <- colMeans(X_mat[c_rows, , drop = FALSE])
    }
    cent_mat
  } else {
    NULL
  }

  basionym_vec_all <- .canonicalize_type_basionym_vec(df, taxa_vec)
  is_noise_vec <- if (!is.null(df) && "is_confident_label_noise" %in% names(df)) {
    as.logical(df$is_confident_label_noise)
  } else {
    rep(FALSE, n)
  }
  is_noise_vec[is.na(is_noise_vec)] <- FALSE
  compute_pool_centroid <- function(rows_sub, X_src = X_mat, S_inv_src = S_inv, prefer_efa = TRUE) {
    .compute_type_pool_centroid(rows_sub, X_src, S_inv_src, df, has_efa_row, type_rank, prefer_efa = prefer_efa)
  }
  cent_ctx <- .build_taxon_type_basionym_centroids(
    df = df, X_mat = X_mat, S_inv = S_inv, active_taxa = active_taxa, taxa_vec = taxa_vec,
    is_type_specimen = is_type_specimen, is_t1_pure = is_t1_pure, has_efa_row = has_efa_row,
    type_rank = type_rank, det_tier = det_tier, is_noise_vec = is_noise_vec, auth_w = auth_w,
    basionym_vec_all = basionym_vec_all, gmm_cl = gmm_cl, u_clusters = u_clusters, cl_centroids = cl_centroids
  )
  type_centroids <- cent_ctx$type_centroids; has_emp_map <- cent_ctx$has_emp_map
  has_emp_homo_map <- cent_ctx$has_emp_homo_map; ss_type_rows_map <- cent_ctx$ss_type_rows_map
  typical_cl_map <- cent_ctx$typical_cl_map; typical_cl_list <- cent_ctx$typical_cl_list
  n_eff_type_vec <- cent_ctx$n_eff_type_vec; type_anchor_rows_list <- cent_ctx$type_anchor_rows_list
  basionym_entries <- cent_ctx$basionym_entries

  n_bas <- length(basionym_entries)
  bas_keys <- if (n_bas > 0L) {
    raw_bk <- vapply(basionym_entries, `[[`, character(1L), "basionym")
    make.unique(raw_bk, sep = "__")
  } else {
    character(0)
  }
  basionym_type_centroids <- if (n_bas > 0L) {
    b_mat <- do.call(rbind, lapply(basionym_entries, `[[`, "mu"))
    rownames(b_mat) <- bas_keys
    colnames(b_mat) <- colnames(X_mat)
    attr(b_mat, "accepted_taxon") <- vapply(basionym_entries, `[[`, character(1L), "taxon")
    attr(b_mat, "basionym_name") <- vapply(basionym_entries, `[[`, character(1L), "basionym")
    attr(b_mat, "is_homotypic") <- vapply(basionym_entries, `[[`, logical(1L), "is_homotypic")
    b_mat
  } else {
    type_centroids
  }

  d2_bas_mat <- if (n_bas > 0L) matrix(NA_real_, nrow = n, ncol = n_bas, dimnames = list(NULL, bas_keys)) else NULL
  if (n_bas > 0L) {
    for (b_idx in seq_len(n_bas)) {
      be <- basionym_entries[[b_idx]]
      diff_b <- sweep(X_mat, 2L, be$mu, "-")
      raw_b <- rowSums((diff_b %*% S_inv) * diff_b)
      raw_b[!is.finite(raw_b) | raw_b < 0.0] <- 0.0
      scale_b <- 1.0 + (1.0 / max(1.0, be$n_eff))
      sc_b <- raw_b / scale_b
      if (length(be$pool) > 0L) sc_b[be$pool] <- raw_b[be$pool]
      d2_bas_mat[, b_idx] <- sc_b
    }
  }

  emp_type_all <- which(is_type_specimen)
  if (any(grepl("^Morphotype_K", active_taxa)) && length(emp_type_all) > 0L) {
    u_emp_bas <- unique(basionym_vec_all[emp_type_all])
    emp_bas_d2 <- matrix(NA_real_, nrow = n, ncol = length(u_emp_bas), dimnames = list(NULL, u_emp_bas))
    emp_bas_mu <- matrix(0.0, nrow = length(u_emp_bas), ncol = p, dimnames = list(u_emp_bas, colnames(X_mat)))
    for (ebi in seq_along(u_emp_bas)) {
      eb_nm <- u_emp_bas[ebi]
      eb_rows <- emp_type_all[basionym_vec_all[emp_type_all] == eb_nm]
      cp_eb <- compute_pool_centroid(eb_rows)
      emp_bas_mu[ebi, ] <- cp_eb$mu
      diff_eb <- sweep(X_mat, 2L, cp_eb$mu, "-")
      raw_eb <- rowSums((diff_eb %*% S_inv) * diff_eb)
      raw_eb[!is.finite(raw_eb) | raw_eb < 0.0] <- 0.0
      sc_eb <- raw_eb / (1.0 + (1.0 / max(1.0, cp_eb$n_eff)))
      if (length(cp_eb$pool) > 0L) sc_eb[cp_eb$pool] <- raw_eb[cp_eb$pool]
      emp_bas_d2[, ebi] <- sc_eb
    }
    basionym_type_centroids <- emp_bas_mu
    eb_near_idx <- max.col(-emp_bas_d2, ties.method = "first")
    nearest_type_basionym <- u_emp_bas[eb_near_idx]
    nearest_type_basionym_d2 <- round(emp_bas_d2[cbind(seq_len(n), eb_near_idx)], 4)
  } else if (!is.null(d2_bas_mat) && ncol(d2_bas_mat) >= 1L) {
    b_near_idx <- max.col(-d2_bas_mat, ties.method = "first")
    bas_names_clean <- vapply(basionym_entries, `[[`, character(1L), "basionym")
    nearest_type_basionym <- bas_names_clean[b_near_idx]
    nearest_type_basionym_d2 <- round(d2_bas_mat[cbind(seq_len(n), b_near_idx)], 4)
  } else {
    nearest_type_basionym <- rep(NA_character_, n)
    nearest_type_basionym_d2 <- rep(NA_real_, n)
  }

  d2_mat <- matrix(NA_real_, nrow = n, ncol = G, dimnames = list(NULL, active_taxa))
  for (k in seq_len(G)) {
    tx_k <- active_taxa[k]
    diff_k <- sweep(X_mat, 2L, type_centroids[k, ], "-")
    raw_d2 <- rowSums((diff_k %*% S_inv) * diff_k)
    raw_d2[!is.finite(raw_d2) | raw_d2 < 0.0] <- 0.0
    scale_k <- 1.0 + (1.0 / max(1.0, n_eff_type_vec[tx_k]))
    scaled_d2 <- raw_d2 / scale_k
    anchor_rows <- type_anchor_rows_list[[tx_k]]
    if (length(anchor_rows) > 0L) {
      scaled_d2[anchor_rows] <- raw_d2[anchor_rows]
    }
    d2_mat[, k] <- scaled_d2
  }

  given_idx <- match(taxa_vec, active_taxa)
  nearest_idx <- if (G >= 1L) max.col(-d2_mat, ties.method = "first") else rep(NA_integer_, n)
  nearest_type_morph_taxon <- if (G >= 1L) active_taxa[nearest_idx] else rep(NA_character_, n)

  type_morph_mahal_d2 <- rep(NA_real_, n)
  valid_g <- which(!is.na(given_idx))
  if (length(valid_g) > 0L) {
    type_morph_mahal_d2[valid_g] <- round(d2_mat[cbind(valid_g, given_idx[valid_g])], 4)
  }

  type_morph_concordance <- rep(NA_real_, n)
  if (G >= 2L && length(valid_g) > 0L) {
    log_lik <- -0.5 * d2_mat
    max_ll <- apply(log_lik, 1L, max, na.rm = TRUE)
    lik_mat <- exp(log_lik - max_ll)
    lik_mat <- lik_mat / pmax(rowSums(lik_mat, na.rm = TRUE), 1e-12)
    type_morph_concordance[valid_g] <- round(lik_mat[cbind(valid_g, given_idx[valid_g])], 4)
  } else if (G == 1L && length(valid_g) > 0L) {
    type_morph_concordance[valid_g] <- round(1.0 - stats::pchisq(type_morph_mahal_d2[valid_g], df = max(1L, p)), 4)
  }

  has_empirical_type_anchor <- ifelse(!is.na(given_idx), unname(has_emp_map[active_taxa[given_idx]]), FALSE)
  is_typical_morphotype <- rep(NA, n)
  chisq_typ_90 <- stats::qchisq(0.90, df = max(1L, p))
  if (length(valid_g) > 0L) {
    for (ii in valid_g) {
      tx_ii <- active_taxa[given_idx[ii]]
      typ_cls_ii <- typical_cl_list[[tx_ii]]
      has_cl_ii <- !is.na(gmm_cl[ii])
      has_d2_ii <- is.finite(type_morph_mahal_d2[ii])
      in_typ_cl <- has_cl_ii && (length(typ_cls_ii) > 0L) && (gmm_cl[ii] %in% typ_cls_ii)
      in_d2_90 <- has_d2_ii && (type_morph_mahal_d2[ii] <= chisq_typ_90)
      if (in_typ_cl || in_d2_90) {
        is_typical_morphotype[ii] <- TRUE
      } else if (has_cl_ii || has_d2_ii) {
        is_typical_morphotype[ii] <- FALSE
      }
    }
  }

  osr_res <- .compute_open_set_recognition_metrics(df, X_mat, S_inv, d2_mat, active_taxa, taxa_vec, gmm_cl)

  eff_df_vec <- rep(as.integer(max(1L, p)), n)
  if (length(no_efa_idx) > 0L) {
    sh_upd <- .eval_type_control_shared_subspace(
      df = df, sh = sh, no_efa_idx = no_efa_idx, n = n, G = G, active_taxa = active_taxa,
      taxa_vec = taxa_vec, ss_type_rows_map = ss_type_rows_map, det_tier = det_tier,
      is_t1_pure = is_t1_pure, is_noise_vec = is_noise_vec, auth_w = auth_w,
      emp_type_all = emp_type_all, basionym_vec_all = basionym_vec_all, n_bas = n_bas,
      basionym_entries = basionym_entries, gmm_cl = gmm_cl, given_idx = given_idx,
      typical_cl_list = typical_cl_list, has_efa_row = has_efa_row, type_rank = type_rank,
      eff_df_vec = eff_df_vec, type_morph_mahal_d2 = type_morph_mahal_d2,
      nearest_type_morph_taxon = nearest_type_morph_taxon, nearest_type_basionym = nearest_type_basionym,
      nearest_type_basionym_d2 = nearest_type_basionym_d2, type_morph_concordance = type_morph_concordance,
      is_typical_morphotype = is_typical_morphotype, osr_res = osr_res
    )
    eff_df_vec <- sh_upd$eff_df_vec; type_morph_mahal_d2 <- sh_upd$type_morph_mahal_d2
    nearest_type_morph_taxon <- sh_upd$nearest_type_morph_taxon
    nearest_type_basionym <- sh_upd$nearest_type_basionym
    nearest_type_basionym_d2 <- sh_upd$nearest_type_basionym_d2
    type_morph_concordance <- sh_upd$type_morph_concordance
    is_typical_morphotype <- sh_upd$is_typical_morphotype; osr_res <- sh_upd$osr_res
  }

  if (!is.null(flags) && !is.null(triage) && !is.null(reasons) && any(is_type_specimen)) {
    grp_vec <- if (!is.null(cda_res) && !is.null(cda_res$grouping_col)) as.character(cda_res$grouping_col) else taxa_vec
    pred_vec <- if (!is.null(cda_res) && !is.null(cda_res$predicted_taxon)) as.character(cda_res$predicted_taxon) else taxa_vec
    post_vec <- if (!is.null(cda_res) && !is.null(cda_res$posterior_prob)) suppressWarnings(as.numeric(cda_res$posterior_prob)) else rep(1.0, n)
    for (i in which(is_type_specimen)) {
      flags[i] <- FALSE
      r_sp <- if (!is.na(grp_vec[i]) && nzchar(grp_vec[i])) grp_vec[i] else taxa_vec[i]
      p_sp <- if (!is.na(pred_vec[i]) && nzchar(pred_vec[i])) pred_vec[i] else r_sp
      pp <- post_vec[i]
      tx_i <- taxa_vec[i]
      has_homo_ref_i <- !is.na(tx_i) && (
        ((tx_i %in% names(has_emp_homo_map)) && isTRUE(has_emp_homo_map[tx_i])) ||
        ((tx_i %in% names(n_eff_type_vec)) && is.finite(n_eff_type_vec[tx_i]) && n_eff_type_vec[tx_i] >= 2.0 && sum(!is_type_specimen & taxa_vec == tx_i) >= 3L)
      )
      is_het_type_i <- !.is_homotypic_sensu_stricto_basionym(
        basionym_vec_all[i], tx_i, df = df[i, , drop = FALSE]
      )
      eff_df_i <- max(1L, as.integer(eff_df_vec[i]))
      if (r_sp != p_sp && is.finite(pp) && pp >= 0.85) {
        triage[i] <- "CRITICAL"
        reasons[i] <- sprintf("Type_Specimen_Taxonomic_Overlap_with_%s", gsub(" ", "_", p_sp))
      } else if (isTRUE(is_het_type_i) && isTRUE(has_homo_ref_i) &&
                 is.finite(type_morph_mahal_d2[i]) &&
                 type_morph_mahal_d2[i] > stats::qchisq(0.99, df = eff_df_i)) {
        triage[i] <- "HIGH"
        reasons[i] <- sprintf("Heterotypic_Type_Divergence_%s", gsub("[^A-Za-z0-9]+", "_", basionym_vec_all[i]))
      } else {
        triage[i] <- "CLEAN"
        reasons[i] <- "Verified_Tier1_Gold_Anchor"
      }
    }
  }

  list(
    type_status = as.character(type_status),
    is_type_specimen = as.logical(is_type_specimen),
    type_rank = as.integer(type_rank),
    has_empirical_type_anchor = as.logical(has_empirical_type_anchor),
    type_morph_mahal_d2 = as.numeric(type_morph_mahal_d2),
    nearest_type_morph_taxon = as.character(nearest_type_morph_taxon),
    nearest_type_basionym = as.character(nearest_type_basionym),
    nearest_type_basionym_d2 = as.numeric(nearest_type_basionym_d2),
    type_morph_concordance = as.numeric(type_morph_concordance),
    is_typical_morphotype = as.logical(is_typical_morphotype),
    min_congener_mahal_d2 = as.numeric(osr_res$min_congener_mahal_d2),
    min_congener_taxon = as.character(osr_res$min_congener_taxon),
    open_set_novelty_p = as.numeric(osr_res$open_set_novelty_p),
    open_set_df = as.integer(osr_res$open_set_df),
    type_centroids = type_centroids,
    basionym_type_centroids = basionym_type_centroids,
    flags = flags,
    triage = triage,
    reasons = reasons
  )
}
