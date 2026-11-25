if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_gmm_mclust_fit.R
# Authority weight config, Stratified Silver Pool Cap, and parallel Mclust GMM fitting
# ==============================================================================

.AUTH_WEIGHT_CFG_ENV <- new.env(parent = emptyenv())
.load_authority_weight_config <- function(cfg_path = "config/config.yaml") {
  if (exists("cached", envir = .AUTH_WEIGHT_CFG_ENV, inherits = FALSE)) {
    return(get("cached", envir = .AUTH_WEIGHT_CFG_ENV, inherits = FALSE))
  }
  defs <- list(
    silver_max_pool_share = 0.40,
    min_type_gold_for_cap = 5L,
    tier_0_type = 3.0,
    tier_0_syntype = 2.5,
    tier_1_gold = 2.0,
    tier_2_silver = 0.75,
    tier_3_bronze_promoted = 0.50,
    tier_3_bronze = 0.35
  )
  if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    aw <- if (is.list(cfg) && is.list(cfg$thresholds) && is.list(cfg$thresholds$authority_weights)) {
      cfg$thresholds$authority_weights
    } else if (is.list(cfg) && is.list(cfg$authority_weights)) {
      cfg$authority_weights
    } else {
      NULL
    }
    if (is.list(aw)) {
      if (!is.null(aw$silver_max_pool_share)) defs$silver_max_pool_share <- min(0.85, max(0.05, as.numeric(aw$silver_max_pool_share)))
      if (!is.null(aw$min_type_gold_for_cap)) defs$min_type_gold_for_cap <- max(1L, as.integer(aw$min_type_gold_for_cap))
      if (!is.null(aw$tier_0_type)) defs$tier_0_type <- max(0.5, as.numeric(aw$tier_0_type))
      if (!is.null(aw$tier_0_syntype)) defs$tier_0_syntype <- max(0.5, as.numeric(aw$tier_0_syntype))
      if (!is.null(aw$tier_1_gold)) defs$tier_1_gold <- max(0.2, as.numeric(aw$tier_1_gold))
      if (!is.null(aw$tier_2_silver)) defs$tier_2_silver <- max(0.05, as.numeric(aw$tier_2_silver))
      if (!is.null(aw$tier_3_bronze_promoted)) defs$tier_3_bronze_promoted <- max(0.05, as.numeric(aw$tier_3_bronze_promoted))
      if (!is.null(aw$tier_3_bronze)) defs$tier_3_bronze <- max(0.01, as.numeric(aw$tier_3_bronze))
    }
  }
  assign("cached", defs, envir = .AUTH_WEIGHT_CFG_ENV)
  defs
}

.apply_stratified_silver_pool_cap <- function(
  w_vec,
  taxa_vec,
  det_tier,
  is_type_mask,
  is_promoted = NULL,
  act_mask = NULL,
  silver_max_share = 0.40,
  min_type_gold = 5L
) {
  n <- length(w_vec)
  capped_mask <- rep(FALSE, n)
  if (n == 0L || !is.finite(silver_max_share) || silver_max_share <= 0.0 || silver_max_share >= 1.0) {
    return(list(weights = w_vec, capped_mask = capped_mask))
  }
  if (is.null(act_mask)) act_mask <- rep(TRUE, n)
  if (is.null(is_promoted)) is_promoted <- rep(FALSE, n)
  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  u_tx <- unique(taxa_vec[act_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_tx)])
  ratio_cap <- silver_max_share / max(1e-6, 1.0 - silver_max_share)
  for (tx in u_tx) {
    idx_tx_act <- which(act_mask & (taxa_vec == tx))
    idx_tg <- idx_tx_act[is_type_mask[idx_tx_act] | (det_tier[idx_tx_act] == "Tier_1_Gold")]
    idx_sv <- idx_tx_act[!is_type_mask[idx_tx_act] & (det_tier[idx_tx_act] == "Tier_2_Silver")]
    if (length(idx_tg) >= as.integer(min_type_gold) && length(idx_sv) >= 1L) {
      w_tg_sum <- sum(w_vec[idx_tg], na.rm = TRUE)
      w_sv_sum <- sum(w_vec[idx_sv], na.rm = TRUE)
      max_sv_sum <- ratio_cap * w_tg_sum
      if (is.finite(w_sv_sum) && is.finite(max_sv_sum) && max_sv_sum > 0.0 && w_sv_sum > (max_sv_sum + 1e-12)) {
        w_vec[idx_sv] <- w_vec[idx_sv] * (max_sv_sum / w_sv_sum)
        capped_mask[idx_sv] <- TRUE
        idx_pr <- idx_tx_act[is_promoted[idx_tx_act]]
        if (length(idx_pr) >= 1L) {
          sv_min <- min(w_vec[idx_sv], na.rm = TRUE)
          w_vec[idx_pr] <- pmin(w_vec[idx_pr], 0.85 * sv_min)
          capped_mask[idx_pr] <- TRUE
        }
      }
    }
  }
  list(weights = w_vec, capped_mask = capped_mask)
}

compute_authority_precision_weights <- function(df, tier_col = "determiner_tier", active_mask = NULL, target_modality = "morph") {
  if (is.logical(tier_col) && is.null(active_mask)) {
    active_mask <- tier_col
    tier_col <- "determiner_tier"
  }
  compute_empirical_bayes_authority_weights(df, mode = "active_passive", active_mask = active_mask, target_modality = target_modality)
}

resolve_ecological_authority_weights <- function(df, allow_bronze_prior = FALSE) {
  n <- if (is.null(df) || !is.data.frame(df)) 0L else nrow(df)
  if (n == 0L) return(numeric(0))
  w_eb <- compute_empirical_bayes_authority_weights(
    df,
    mode = if (isTRUE(allow_bronze_prior)) "continuous_prior_posterior" else "active_passive",
    target_modality = "ecological"
  )
  if (!isTRUE(allow_bronze_prior) && "is_confident_label_noise" %in% names(df)) {
    cl_noise <- as.logical(df$is_confident_label_noise)
    cl_noise[is.na(cl_noise)] <- FALSE
    ti <- .identify_type_specimen_vec(df)
    w_eb[cl_noise & !ti$is_type_specimen] <- 0.0
  }
  round(w_eb, 4)
}

screen_gmm_point_mass_columns <- function(X_mat, ref_mask = NULL, max_tie_frac = 0.25) {
  X_mat <- as.matrix(X_mat)
  cn <- colnames(X_mat)
  if (is.null(cn) || ncol(X_mat) <= 2L) return(seq_len(ncol(X_mat)))
  ref_idx <- if (!is.null(ref_mask) && length(ref_mask) == nrow(X_mat) && sum(ref_mask, na.rm = TRUE) >= 10L) {
    which(as.logical(ref_mask))
  } else {
    seq_len(nrow(X_mat))
  }
  if (length(ref_idx) < 30L) return(cn)
  keep <- vapply(seq_len(ncol(X_mat)), function(j) {
    if (grepl("^PC[0-9]+$", cn[j])) return(TRUE)
    v <- X_mat[ref_idx, j]
    v <- round(v[is.finite(v)], 5)
    if (length(v) < 10L) return(FALSE)
    (max(tabulate(match(v, unique(v)))) / length(v)) <= max_tie_frac
  }, logical(1L))
  if (sum(keep) < 2L) return(cn)
  cn[keep]
}

.max_column_tie_count <- function(X_mat) {
  X_mat <- as.matrix(X_mat)
  if (nrow(X_mat) == 0L || ncol(X_mat) == 0L) return(0L)
  max(vapply(seq_len(ncol(X_mat)), function(j) {
    v <- round(X_mat[is.finite(X_mat[, j]), j], 5)
    if (length(v) == 0L) 0L else max(tabulate(match(v, unique(v))))
  }, integer(1L)))
}

run_bounded_clustvarsel <- function(X_mat, core_idx, max_k = 8L, timeout_sec = 20.0, seed = 42L) {
  if (!requireNamespace("clustvarsel", quietly = TRUE)) return(NULL)
  X_mat <- as.matrix(X_mat)
  sub_idx <- if (length(core_idx) > 350L) {
    core_idx[unique(round(seq(1L, length(core_idx), length.out = 350L)))]
  } else {
    core_idx
  }
  X_sub <- X_mat[sub_idx, , drop = FALSE]
  p_dim <- ncol(X_sub)
  if (nrow(X_sub) < 10L || p_dim < 3L) return(NULL)
  has_ties <- .max_column_tie_count(X_sub) > p_dim
  G_cvs <- seq_len(min(as.integer(max_k), 5L))

  run_cvs_inner <- function() {
    set.seed(as.integer(seed))
    if (has_ties) {
      clustvarsel::clustvarsel(
        X_sub, G = G_cvs, search = "headlong",
        hcModel = "VVI", emModels2 = c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
        itermax = 100, parallel = FALSE, verbose = FALSE
      )
    } else {
      clustvarsel::clustvarsel(X_sub, G = G_cvs, itermax = 100, parallel = FALSE, verbose = FALSE)
    }
  }

  cvs_res <- if (.Platform$OS.type == "unix") {
    jb <- parallel::mcparallel(tryCatch(run_cvs_inner(), error = function(e) NULL), silent = TRUE)
    col_res <- suppressWarnings(parallel::mccollect(jb, wait = FALSE, timeout = timeout_sec))
    if (is.null(col_res)) {
      tryCatch(tools::pskill(jb$pid, signal = 9L), error = function(e) NULL)
      tryCatch(suppressWarnings(parallel::mccollect(jb, wait = FALSE, timeout = 0.2)), error = function(e) NULL)
      NULL
    } else {
      val <- col_res[[1L]]
      if (inherits(val, "try-error")) NULL else val
    }
  } else {
    tryCatch(run_cvs_inner(), error = function(e) NULL)
  }

  if (is.null(cvs_res) || length(cvs_res$subset) < 2L) return(NULL)
  sub_pos <- as.integer(cvs_res$subset)
  if (has_ties && !is.null(colnames(X_mat))) {
    lead_pc <- which(colnames(X_mat) == "PC1")
    if (length(lead_pc) > 0L && !(lead_pc[1L] %in% sub_pos)) {
      sub_pos <- c(lead_pc[1L], sub_pos)
    }
  }
  cvs_res$subset <- sort(unique(sub_pos))
  cvs_res
}

fit_parallel_mclust_gmm <- function(X_core, G_vec, modelNames = NULL, max_hc_subset = 1500L) {
  if (!requireNamespace("mclust", quietly = TRUE)) return(NULL)
  if (!("package:mclust" %in% search())) {
    tryCatch(suppressPackageStartupMessages(library(mclust, quietly = TRUE, warn.conflicts = FALSE)), error = function(e) NULL)
  }
  X_core <- as.matrix(X_core)
  n_obs <- nrow(X_core)
  p_dim <- ncol(X_core)
  G_int <- as.integer(G_vec)
  has_ties <- .max_column_tie_count(X_core) > p_dim
  ctrl <- mclust::emControl(itmax = c(180L, 180L), tol = c(1e-5, sqrt(.Machine$double.eps)))

  if (n_obs <= 40L && !has_ties) {
    return(mclust::Mclust(X_core, G = G_int, modelNames = modelNames, control = ctrl, verbose = FALSE))
  }

  prior_supported <- if (p_dim == 1L) {
    c("E", "V")
  } else {
    c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV")
  }
  cand_models <- if (!is.null(modelNames) && length(modelNames) > 0L) {
    as.character(modelNames)
  } else {
    as.character(mclust::mclust.options("emModelNames"))
  }
  if (p_dim == 1L) cand_models <- intersect(cand_models, c("E", "V"))
  if (has_ties) cand_models <- intersect(cand_models, prior_supported)
  if (length(cand_models) == 0L) cand_models <- if (p_dim == 1L) c("E", "V") else c("EII", "VII", "EEI", "VEI", "EVI", "VVI")

  prior_arg <- if (has_ties) mclust::priorControl() else NULL
  n_sub <- min(n_obs, max(600L, as.integer(max_hc_subset)))
  init_arg <- if (n_obs > n_sub) {
    ord1 <- order(X_core[, 1L])
    list(subset = sort(ord1[unique(round(seq(1, n_obs, length.out = n_sub)))]))
  } else {
    NULL
  }

  grid_df <- expand.grid(G = G_int, modelName = cand_models, stringsAsFactors = FALSE)
  n_tasks <- nrow(grid_df)
  fit_one_cell <- function(idx) {
    g_i <- grid_df$G[idx]
    mn_i <- grid_df$modelName[idx]
    tryCatch(
      mclust::Mclust(
        X_core, G = g_i, modelNames = mn_i,
        prior = prior_arg, control = ctrl, initialization = init_arg, verbose = FALSE
      ),
      error = function(e) NULL
    )
  }

  use_par <- exists("packera_mclapply", mode = "function") && n_tasks >= 2L
  fits <- if (use_par) {
    n_w <- if (exists("resolve_packera_workers", mode = "function")) {
      resolve_packera_workers(64L, n_tasks)
    } else {
      min(8L, n_tasks)
    }
    packera_mclapply(seq_len(n_tasks), fit_one_cell, mc.cores = n_w)
  } else {
    lapply(seq_len(n_tasks), fit_one_cell)
  }

  valid_idx <- which(vapply(fits, function(f) !is.null(f) && !is.null(f$bic) && is.finite(f$bic[1L]), logical(1L)))
  if (length(valid_idx) == 0L) {
    return(mclust::Mclust(
      X_core, G = G_int, modelNames = intersect(cand_models, prior_supported),
      prior = mclust::priorControl(), control = ctrl, initialization = init_arg, verbose = FALSE
    ))
  }

  bic_mat <- matrix(
    NA_real_,
    nrow = length(G_int),
    ncol = length(cand_models),
    dimnames = list(as.character(G_int), cand_models)
  )
  for (idx in valid_idx) {
    r_nm <- as.character(grid_df$G[idx])
    c_nm <- as.character(grid_df$modelName[idx])
    bic_mat[r_nm, c_nm] <- as.numeric(fits[[idx]]$bic[1L])
  }
  bics <- vapply(valid_idx, function(idx) as.numeric(fits[[idx]]$bic[1L]), numeric(1L))
  best_fit <- fits[[valid_idx[which.max(bics)]]]
  best_fit$BIC <- bic_mat
  best_fit
}
