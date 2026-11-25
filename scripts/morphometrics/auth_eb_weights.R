if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_eb_weights.R
# Hierarchical Beta-Binomial Empirical Bayes & Dawid-Skene authority weights
# ==============================================================================

compute_empirical_bayes_authority_weights <- function(
  df,
  mode = c("active_passive", "continuous_prior_posterior"),
  active_mask = NULL,
  prior_alpha = 5.0,
  anchor_info = NULL,
  target_modality = c("all", "morph", "vision", "ecological"),
  taxon_col = "species_standardized",
  ...
) {
  mode <- match.arg(mode)
  target_modality <- match.arg(target_modality, c("all", "morph", "vision", "ecological"))
  if (is.null(df) || nrow(df) == 0L) return(numeric(0))
  n <- nrow(df)
  if (is.null(active_mask) && !is.null(anchor_info) && is.list(anchor_info) && !is.null(anchor_info$active_mask)) {
    active_mask <- anchor_info$active_mask
  }
  alpha0 <- max(1e-6, as.numeric(prior_alpha))

  type_info <- .identify_type_specimen_vec(df)
  type_rank <- type_info$type_rank
  is_type_mask <- type_info$is_type_specimen
  det_tier <- .standardize_determiner_tier_vec(df)
  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  det_id <- .normalize_determiner_id_vec(df)

  is_promoted <- rep(FALSE, n)
  if (!is.null(active_mask) && length(active_mask) == n) {
    act_in <- as.logical(active_mask)
    act_in[is.na(act_in)] <- FALSE
    is_promoted <- (det_tier == "Tier_3_Bronze") & act_in
  } else if ("is_promoted_bronze" %in% names(df)) {
    is_promoted <- as.logical(df$is_promoted_bronze)
    is_promoted[is.na(is_promoted)] <- FALSE
  } else if ("anchor_tier_source" %in% names(df)) {
    is_promoted <- as.character(df$anchor_tier_source) == "Tier_3_Bronze_Promoted"
    is_promoted[is.na(is_promoted)] <- FALSE
  }

  act_mask <- if (!is.null(active_mask) && length(active_mask) == n) {
    m <- as.logical(active_mask)
    m[is.na(m)] <- FALSE
    m
  } else {
    is_type_mask | (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted
  }
  if (mode == "active_passive" && target_modality %in% c("all", "morph") && "morph_tier_source" %in% names(df)) {
    mts_eb <- trimws(as.character(df$morph_tier_source))
    act_mask[!is.na(mts_eb) & (mts_eb %in% c("Cauline_Only_Fallback", "Transitional_Passive", "No_Leaves")) & !is_type_mask] <- FALSE
  }

  extralimital_raw_ds <- is_extralimital_or_conflicted_non_type(df)
  extralimital_ds <- rep(FALSE, n)
  if (any(extralimital_raw_ds)) {
    for (tx_chk in unique(taxa_vec[!is.na(taxa_vec) & nzchar(taxa_vec)])) {
      if (sum(act_mask & (taxa_vec == tx_chk) & !extralimital_raw_ds) >= 2L) {
        extralimital_ds[(taxa_vec == tx_chk) & extralimital_raw_ds] <- TRUE
      }
    }
  }
  is_open_hyb_non_type_ds <- .identify_open_or_hybrid_concept_mask(df) & !is_type_mask
  act_mask[extralimital_ds | is_open_hyb_non_type_ds] <- FALSE

  if (mode == "active_passive" && !any(act_mask)) {
    return(rep(0.0, n))
  }
  eval_mask <- if (mode == "active_passive") act_mask else rep(TRUE, n)

  tier_cell <- ifelse(
    type_rank == 2L, "Tier_0_Primary_Type",
    ifelse(
      type_rank == 1L, "Tier_0_Secondary_Type",
      ifelse(
        det_tier == "Tier_1_Gold", "Tier_1_Gold",
        ifelse(
          det_tier == "Tier_2_Silver", "Tier_2_Silver",
          ifelse(is_promoted, "Tier_3_Bronze_Promoted", "Tier_3_Bronze")
        )
      )
    )
  )

  prior_map <- c(
    "Tier_0_Primary_Type" = 0.999,
    "Tier_0_Secondary_Type" = 0.985,
    "Tier_1_Gold" = 0.95,
    "Tier_2_Silver" = 0.85,
    "Tier_3_Bronze_Promoted" = 0.65,
    "Tier_3_Bronze" = 0.35
  )
  pi_vec <- unname(prior_map[tier_cell])
  tier_ceil <- pi_vec + 0.35 * (1.0 - pi_vec)

  blend_tier_and_determiner <- function(conc_vec, obs_mask) {
    .eb_blend_tier_and_determiner(conc_vec, obs_mask, pi_vec, eval_mask, tier_cell, alpha0, taxa_vec, det_id)
  }

  eps_err <- .extract_specimen_discordance_signals(df, eval_mask = eval_mask, target_modality = target_modality)

  is_art <- if ("is_artifact_driven" %in% names(df)) {
    av <- as.logical(df$is_artifact_driven)
    av[is.na(av)] <- FALSE
    av
  } else {
    rep(FALSE, n)
  }

  safe_std_tx <- function(v) {
    tx_v <- trimws(as.character(v))
    if (exists("standardize_packera_taxon", mode = "function")) {
      tx_v <- as.character(standardize_packera_taxon(tx_v))
    } else if (exists("standardize_taxon_name", mode = "function")) {
      std_v <- as.character(standardize_taxon_name(tx_v))
      tx_v <- ifelse(!is.na(std_v) & nzchar(std_v), std_v, tx_v)
    }
    tx_v
  }

  pred_sources <- list()
  if ("exsiccatae_gold_taxon" %in% names(df)) {
    p_ex <- safe_std_tx(df$exsiccatae_gold_taxon)
    valid_ex <- !is.na(p_ex) & nzchar(p_ex) & !(tolower(p_ex) %in% c("na", "none", "unknown", "unassigned_packera"))
    if (any(valid_ex)) {
      pred_sources[["exsiccatae"]] <- list(taxon = p_ex, conf = rep(0.92, n), valid = valid_ex)
    }
  }
  if (!(target_modality %in% c("morph", "ecological")) && "cda_predicted_taxon" %in% names(df)) {
    p_cda <- safe_std_tx(df$cda_predicted_taxon)
    c_cda <- if ("cda_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(0.82, n)
    c_cda[!is.finite(c_cda)] <- 0.82
    pred_sources[["cda"]] <- list(taxon = p_cda, conf = pmax(0.55, pmin(0.92, c_cda)), valid = rep(TRUE, n))
  }
  vis_col <- if ("vision_predicted_label" %in% names(df)) "vision_predicted_label" else if ("predicted_label" %in% names(df)) "predicted_label" else NULL
  if (target_modality != "vision" && !is.null(vis_col)) {
    p_vis <- safe_std_tx(df[[vis_col]])
    pred_sources[["vision"]] <- list(taxon = p_vis, conf = rep(0.82, n), valid = !is_art)
  }
  if (target_modality != "ecological") {
    if ("rf_predicted_taxon" %in% names(df)) {
      p_rf <- safe_std_tx(df$rf_predicted_taxon)
      pred_sources[["rf"]] <- list(taxon = p_rf, conf = rep(0.78, n), valid = rep(TRUE, n))
    }
    if ("edaphic_best_fit_taxon" %in% names(df)) {
      p_ed <- safe_std_tx(df$edaphic_best_fit_taxon)
      pred_sources[["edaphic"]] <- list(taxon = p_ed, conf = rep(0.75, n), valid = rep(TRUE, n))
    }
  }

  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  valid_taxa <- sort(unique(taxa_vec[eval_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_tx)]))
  G <- length(valid_taxa)
  exog_cols <- intersect(
    c("exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon", "vision_predicted_label", "is_label_corrupted", "rf_predicted_taxon", "edaphic_best_fit_taxon"),
    names(df)
  )
  has_exog_signals <- length(exog_cols) > 0L && any(vapply(exog_cols, function(cl) any(!is.na(df[[cl]]) & nzchar(trimws(as.character(df[[cl]])))), logical(1)))
  use_loo_pc <- (G >= 2L) && !(target_modality %in% c("morph", "ecological"))
  L_pc <- if (use_loo_pc) .compute_loo_pc_likelihoods(df, taxa_vec, det_tier, is_promoted, eval_mask, valid_taxa, pi_vec) else NULL

  y_idx <- match(taxa_vec, valid_taxa)
  ds_lik <- .eb_build_dawid_skene_likelihood_matrix(
    n = n, G = G, valid_taxa = valid_taxa, y_idx = y_idx, eval_mask = eval_mask,
    extralimital_ds = extralimital_ds, eps_err = eps_err, pred_sources = pred_sources,
    bad_tx = bad_tx, L_pc = L_pc, det_tier = det_tier
  )
  L_mat <- ds_lik$L_mat
  has_ds_Row <- ds_lik$has_ds_Row

  w_eb <- pi_vec

  if (G >= 2L && any(has_ds_Row)) {
    w_eb <- .eb_run_dawid_skene_em_update(
      w_eb = w_eb, n = n, G = G, valid_taxa = valid_taxa, y_idx = y_idx,
      L_mat = L_mat, has_ds_Row = has_ds_Row, pi_vec = pi_vec, tier_ceil = tier_ceil,
      alpha0 = alpha0, eval_mask = eval_mask, eps_err = eps_err,
      extralimital_ds = extralimital_ds, blend_tier_and_determiner = blend_tier_and_determiner
    )
  } else {
    # Graceful single-taxon or feature-absent degradation: Beta-Binomial tier + determiner shrinkage
    eps_fb <- eps_err
    if (length(pred_sources) > 0L) {
      for (i in which(eval_mask & !is.finite(eps_fb) & !extralimital_ds)) {
        tx_i <- taxa_vec[i]
        if (is.na(tx_i) || (tx_i %in% bad_tx)) next
        mism <- numeric(0)
        for (ps in pred_sources) {
          if (!isTRUE(ps$valid[i])) next
          ptx <- ps$taxon[i]
          if (!is.na(ptx) && !(ptx %in% bad_tx)) {
            mism <- c(mism, if (ptx == tx_i) (1.0 - ps$conf[i]) else ps$conf[i])
          }
        }
        if (length(mism) > 0L) eps_fb[i] <- mean(mism)
      }
    }
    obs_fb <- eval_mask & is.finite(eps_fb) & !extralimital_ds
    if (!any(obs_fb)) obs_fb <- eval_mask & is.finite(eps_fb)
    if (any(obs_fb)) {
      conc_fb <- pmin(pmax(0.0, 1.0 - eps_fb), tier_ceil)
      theta_fb <- blend_tier_and_determiner(conc_fb, obs_fb)
      w_eb[eval_mask] <- theta_fb[eval_mask]
      w_eb[obs_fb] <- pmax(0.05, pmin(0.995, 0.85 * theta_fb[obs_fb] + 0.15 * conc_fb[obs_fb]))
    }
  }

  w_eb[!is.finite(w_eb)] <- pi_vec[!is.finite(w_eb)]
  w_eb <- pmax(0.05, pmin(0.995, w_eb))
  if (any(is_type_mask)) {
    w_eb[is_type_mask] <- pi_vec[is_type_mask]
  }

  aw_cfg <- .load_authority_weight_config()
  is_nb_sec <- (type_rank == 1L) & as.logical(type_info$is_name_bearing_type)
  is_nb_sec[is.na(is_nb_sec)] <- FALSE
  base_prec <- ifelse(
    type_rank == 2L, aw_cfg$tier_0_type,
    ifelse(
      is_nb_sec, aw_cfg$tier_0_syntype * (w_eb / 0.985),
      ifelse(
        type_rank == 1L, (aw_cfg$tier_1_gold + 0.15) * (w_eb / 0.985),
        ifelse(
          det_tier == "Tier_1_Gold", aw_cfg$tier_1_gold * (w_eb / 0.95),
          ifelse(
            det_tier == "Tier_2_Silver", aw_cfg$tier_2_silver * (w_eb / 0.85),
            ifelse(
              is_promoted, aw_cfg$tier_3_bronze_promoted * (w_eb / 0.65),
              aw_cfg$tier_3_bronze * (w_eb / 0.35)
            )
          )
        )
      )
    )
  )
  base_prec[!is.finite(base_prec) | base_prec <= 0.0] <- 0.05

  if (mode == "active_passive") {
    w_out <- rep(0.0, n)
    if (any(act_mask)) {
      cap_res <- .apply_stratified_silver_pool_cap(
        w_vec = base_prec,
        taxa_vec = taxa_vec,
        det_tier = det_tier,
        is_type_mask = is_type_mask,
        is_promoted = is_promoted,
        act_mask = act_mask,
        silver_max_share = aw_cfg$silver_max_pool_share,
        min_type_gold = aw_cfg$min_type_gold_for_cap
      )
      w_act <- cap_res$weights[act_mask]
      mean_act <- mean(w_act, na.rm = TRUE)
      w_out[act_mask] <- if (is.finite(mean_act) && mean_act > 1e-12) round(w_act / mean_act, 6) else 1.0
    }
    return(w_out)
  }

  max_w <- max(base_prec, na.rm = TRUE)
  w_cont <- if (is.finite(max_w) && max_w > 1e-12) (3.0 * base_prec / max_w) else rep(1.0, n)
  w_cont[!is.finite(w_cont)] <- 1.0
  if (any(type_rank == 2L)) {
    w_cont[type_rank == 2L] <- 3.0
  }
  w_cont <- pmax(0.5, pmin(3.0, w_cont))
  cap_cont <- .apply_stratified_silver_pool_cap(
    w_vec = w_cont,
    taxa_vec = taxa_vec,
    det_tier = det_tier,
    is_type_mask = is_type_mask,
    is_promoted = is_promoted,
    act_mask = act_mask,
    silver_max_share = aw_cfg$silver_max_pool_share,
    min_type_gold = aw_cfg$min_type_gold_for_cap
  )
  w_final <- cap_cont$weights
  if (any(cap_cont$capped_mask)) {
    w_final[cap_cont$capped_mask] <- pmax(0.01, w_final[cap_cont$capped_mask])
    w_final[!cap_cont$capped_mask] <- pmax(0.5, pmin(3.0, w_final[!cap_cont$capped_mask]))
  }
  w_final[extralimital_ds | is_open_hyb_non_type_ds] <- 0.0
  round(w_final, 6)
}
