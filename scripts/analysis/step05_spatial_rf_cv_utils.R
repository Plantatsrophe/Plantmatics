# ==============================================================================
# Module: scripts/analysis/step05_spatial_rf_cv_utils.R
# Purpose: Haversine spatial autocorrelation range, terrestrial accessible
#          background sampling, fold-internal spatial balance & class weighting,
#          small-n EB calibration, and MEM-decoupled Spatial Block-CV RF.
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

compute_inverse_class_frequency_weights <- function(y_vec) {
  y_chr <- as.character(y_vec); tbl <- table(y_chr)
  w_map <- length(y_chr) / (max(length(tbl), 1L) * pmax(as.numeric(tbl), 1.0)); names(w_map) <- names(tbl)
  w_out <- as.numeric(w_map[y_chr]); w_out[is.na(w_out) | !is.finite(w_out)] <- 1.0
  w_out / max(mean(w_out, na.rm = TRUE), 1e-8)
}

.haversine_km_safe <- function(ca, cb = NULL) {
  if (exists("compute_haversine_dist_km_matrix", mode = "function")) return(compute_haversine_dist_km_matrix(ca, cb))
  ca <- as.matrix(ca); cb <- if (is.null(cb)) ca else as.matrix(cb); r <- pi / 180.0
  a <- sin(outer(ca[, 2L] * r, cb[, 2L] * r, "-") * 0.5)^2 + outer(cos(ca[, 2L] * r), cos(cb[, 2L] * r), "*") * sin(outer(ca[, 1L] * r, cb[, 1L] * r, "-") * 0.5)^2
  d <- 2.0 * 6371.0 * asin(sqrt(pmin(pmax(a, 0.0), 1.0))); dim(d) <- c(nrow(ca), nrow(cb)); if (is.null(cb)) diag(d) <- 0.0; d
}

estimate_spatial_autocorrelation_block_size <- function(coords, X_env, default_km = 55.0) {
  n <- if (is.null(coords)) 0L else nrow(coords)
  if (n < 10L) return(list(range_km = default_km, range_m = default_km * 1000.0, range_deg = round(default_km / 111.32, 4)))
  if (n > 2500L) { si <- unique(round(seq(1L, n, length.out = 2500L))); coords <- coords[si, , drop = FALSE]; X_env <- X_env[si, , drop = FALSE] }
  r_km <- tryCatch({
    d_m <- .haversine_km_safe(coords); d_v <- d_m[upper.tri(d_m)]; if (length(d_v) < 5L) return(default_km)
    pc1 <- tryCatch(stats::prcomp(X_env, center = TRUE, scale. = TRUE)$x[, 1L], error = function(e) X_env[, 1L])
    v_v <- outer(pc1, pc1, function(a, b) 0.5 * (a - b)^2)[upper.tri(d_m)]; sill <- stats::var(pc1, na.rm = TRUE)
    if (!is.finite(sill) || sill <= 1e-8) return(default_km)
    brks <- unique(stats::quantile(d_v, probs = seq(0.0, 0.8, length.out = 9L), na.rm = TRUE)); if (length(brks) < 3L) return(default_km)
    bins <- cut(d_v, breaks = brks, include.lowest = TRUE)
    gb <- tapply(v_v, bins, mean, na.rm = TRUE); db <- tapply(d_v, bins, mean, na.rm = TRUE)
    ri <- which(gb >= 0.85 * sill)
    pmax(22.0, pmin(278.0, if (length(ri) > 0L) as.numeric(db[ri[1L]]) else as.numeric(stats::median(d_v, na.rm = TRUE) * 0.25)))
  }, error = function(e) default_km)
  if (!is.finite(r_km) || is.na(r_km)) r_km <- default_km
  list(range_km = round(r_km, 2), range_m = round(r_km * 1000.0, 0), range_deg = round(r_km / 111.32, 4))
}

assign_spatial_grid_blocks <- function(coords, block_deg = 0.50, k = 5L) {
  lon <- as.numeric(coords[, 1L]); lat <- as.numeric(coords[, 2L]); bd <- if (is.finite(block_deg) && block_deg > 0.05) block_deg else 0.50
  cids <- paste0(floor((lon - min(lon, na.rm = TRUE)) / bd), "_", floor((lat - min(lat, na.rm = TRUE)) / bd)); uc <- unique(cids)
  if (length(uc) < 2L) return(rep(seq_len(k), length.out = length(lon)))
  cf <- ((seq_along(uc) - 1L) %% k) + 1L; names(cf) <- uc; as.integer(cf[cids])
}

.is_na_terrestrial_coord <- function(lon, lat) {
  ok <- is.finite(lon) & is.finite(lat) & (lat >= 24.4) & (lat <= 72.0) & (lon >= -168.0) & (lon <= -52.5)
  atl <- ifelse(lat < 27.5, -80.0, ifelse(lat < 30.8, -81.1, ifelse(lat < 33.8, -81.1 + (lat - 30.8) * 1.10,
         ifelse(lat < 36.6, -77.8 + (lat - 33.8) * 0.85, ifelse(lat < 41.2, -75.4 + (lat - 36.6) * 1.15, ifelse(lat < 45.0, -70.1 + (lat - 41.2) * 2.15, -52.5))))))
  gulf <- ifelse(lon > -82.5, 24.5, ifelse(lon > -84.5, 28.8, ifelse(lon > -89.5, 29.8, ifelse(lon > -94.0, 29.2, 26.0))))
  ok & !((lon > -81.5) & (lat < 45.0) & (lon > atl)) & !((lon > -97.5) & (lon < -80.5) & (lat < gulf))
}

.generate_accessible_background_env <- function(coords_pres, X_pres, n_bg = 500L, seed = 42L, env_dir = NULL, demo = FALSE, block_km = 180.0) {
  os <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(as.integer(seed)); on.exit({ if (!is.null(os)) assign(".Random.seed", os, envir = .GlobalEnv) }, add = TRUE)
  p_c <- ncol(X_pres); pred_nms <- colnames(X_pres); lon_p <- as.numeric(coords_pres[, 1L]); lat_p <- as.numeric(coords_pres[, 2L])
  max_km <- max(150.0, if (is.finite(block_km)) as.numeric(block_km) else 180.0); min_km <- 5.0
  nc <- max(as.integer(n_bg) * 6L, 600L); si <- sample.int(length(lon_p), nc, replace = TRUE)
  br <- stats::runif(nc, 0.0, 2.0 * pi); rk <- sqrt(stats::runif(nc, min_km^2, max_km^2))
  cand <- cbind(longitude = lon_p[si] + (rk * sin(br)) / pmax(20.0, 111.32 * cos(lat_p[si] * pi / 180.0)), latitude = lat_p[si] + (rk * cos(br)) / 111.32)
  md_km <- apply(.haversine_km_safe(cand, cbind(lon_p, lat_p)), 1L, min, na.rm = TRUE)
  is_land <- .is_na_terrestrial_coord(cand[, 1L], cand[, 2L]) | (md_km <= 20.0)
  km <- which(is_land & (md_km >= min_km) & (md_km <= max_km))
  if (length(km) < 40L) km <- which(is_land & (md_km <= max_km))
  if (length(km) < 20L) km <- which(md_km <= max_km)
  km_pool <- if (length(km) > as.integer(n_bg) * 3L) sample(km, as.integer(n_bg) * 3L) else km
  bg_pool <- cand[km_pool, , drop = FALSE]
  if (!isTRUE(demo) && !is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir) &&
      exists(".build_regional_sdm_raster_stack", mode = "function") && requireNamespace("terra", quietly = TRUE)) {
    bs <- tryCatch(.build_regional_sdm_raster_stack(data.frame(longitude = lon_p, latitude = lat_p), env_dir = env_dir, kept_vars = pred_nms, demo = FALSE), error = function(e) NULL)
    if (!is.null(bs) && all(pred_nms %in% names(bs))) {
      pb <- tryCatch(terra::vect(as.data.frame(bg_pool), geom = c("longitude", "latitude"), crs = "EPSG:4326"), error = function(e) NULL)
      eb <- if (!is.null(pb)) tryCatch(terra::extract(bs[[pred_nms]], pb), error = function(e) NULL) else NULL
      if (!is.null(eb)) {
        Xr <- as.matrix(eb[, pred_nms, drop = FALSE]); ok <- which(rowSums(!is.finite(Xr)) == 0L)
        if (length(ok) >= 20L) {
          if (length(ok) > as.integer(n_bg)) ok <- sample(ok, as.integer(n_bg))
          return(list(coords = bg_pool[ok, , drop = FALSE], X_bg = Xr[ok, , drop = FALSE]))
        }
      }
    }
  }
  if (nrow(bg_pool) > as.integer(n_bg)) bg_c <- bg_pool[sample.int(nrow(bg_pool), as.integer(n_bg)), , drop = FALSE] else bg_c <- bg_pool
  m_n <- nrow(bg_c)
  W <- exp(-(.haversine_km_safe(bg_c, cbind(lon_p, lat_p)) + 15.0) / 95.0); W <- sweep(W, 1L, pmax(rowSums(W), 1e-8), "/")
  cr <- tryCatch(stats::cov(X_pres), error = function(e) NULL)
  if (is.null(cr) || any(!is.finite(cr))) {
    s <- apply(X_pres, 2L, stats::sd, na.rm = TRUE); s[!is.finite(s) | s <= 1e-6] <- 1.0; cr <- diag(s^2, nrow = p_c, ncol = p_c)
  } else {
    cr <- 0.85 * cr + 0.15 * diag(pmax(diag(cr), 1e-6), nrow = p_c, ncol = p_c)
  }
  eg <- eigen(cr, symmetric = TRUE); Lh <- eg$vectors %*% diag(sqrt(pmax(eg$values, 1e-8)), nrow = p_c, ncol = p_c)
  X_bg <- 0.55 * (W %*% X_pres) + 0.45 * sweep(matrix(0.0, m_n, p_c), 2L, colMeans(X_pres), "+") + 1.15 * (matrix(stats::rnorm(m_n * p_c), nrow = m_n, ncol = p_c) %*% t(Lh))
  colnames(X_bg) <- colnames(X_pres)
  for (j in seq_len(p_c)) {
    ql <- min(X_pres[, j], na.rm = TRUE) - 1.5 * stats::sd(X_pres[, j], na.rm = TRUE); qh <- max(X_pres[, j], na.rm = TRUE) + 1.5 * stats::sd(X_pres[, j], na.rm = TRUE)
    if (grepl("pct$|sand|clay", colnames(X_pres)[j], ignore.case = TRUE)) {
      ql <- max(0.0, ql); qh <- min(100.0, qh)
    } else if (grepl("ph", colnames(X_pres)[j], ignore.case = TRUE)) {
      ql <- max(3.2, ql); qh <- min(8.8, qh)
    } else if (grepl("bedrock|bdricm|resdept", colnames(X_pres)[j], ignore.case = TRUE)) {
      ql <- max(0.0, ql); qh <- min(200.0, qh)
    }
    X_bg[, j] <- pmin(qh, pmax(ql, X_bg[, j]))
  }
  list(coords = bg_c, X_bg = X_bg)
}

.compute_fast_mems <- function(coords_m, max_mems = 5L, try_adespatial = TRUE) {
  n_m <- nrow(coords_m); ok_r <- is.finite(coords_m[, 1L]) & is.finite(coords_m[, 2L])
  if (!all(ok_r)) {
    for (jj in 1:2) { mj <- stats::median(coords_m[ok_r, jj], na.rm = TRUE); coords_m[!ok_r, jj] <- if (is.finite(mj)) mj else 0.0 }
  }
  if (isTRUE(try_adespatial) && n_m <= 900L && all(ok_r) && requireNamespace("adespatial", quietly = TRUE)) {
    ma <- tryCatch({ mm <- as.matrix(adespatial::dbmem(coords_m, MEM.autocor = "positive")); if (ncol(mm) >= 1L) { k <- min(ncol(mm), as.integer(max_mems)); o <- mm[, seq_len(k), drop = FALSE]; colnames(o) <- paste0("MEM", seq_len(k)); o } else NULL }, error = function(e) NULL)
    if (!is.null(ma)) return(ma)
  }
  D <- .haversine_km_safe(coords_m); knn <- min(8L, n_m); th <- max(apply(D, 1L, function(r) sort(r, partial = knn)[knn]))
  W <- exp(-D / (th + 1e-5)) * (D <= th); diag(W) <- 0.0; rw <- rowMeans(W)
  M <- sweep(sweep(W, 1L, rw, "-"), 2L, rw, "-") + mean(rw)
  eg <- if (n_m >= 20L && requireNamespace("RSpectra", quietly = TRUE)) tryCatch(RSpectra::eigs_sym(M, k = min(max(6L, as.integer(max_mems) + 2L), n_m - 1L), which = "LA"), error = function(e) eigen(M, symmetric = TRUE)) else eigen(M, symmetric = TRUE)
  pos <- which(eg$values > 1e-4); if (length(pos) == 0L) pos <- 1L
  k <- min(length(pos), as.integer(max_mems)); o <- eg$vectors[, pos[seq_len(k)], drop = FALSE]; if (!all(ok_r)) o[!ok_r, ] <- 0.0; colnames(o) <- paste0("MEM", seq_len(k)); o
}

.calibrate_small_n_rf_probs <- function(prob_mat, y_train, min_n_calib = 10L) {
  if (is.null(prob_mat) || !is.matrix(prob_mat) || nrow(prob_mat) == 0L || ncol(prob_mat) < 2L) return(prob_mat)
  cn <- colnames(prob_mat); tbl <- table(as.character(y_train))
  nk <- vapply(cn, function(cl) if (cl %in% names(tbl)) as.numeric(tbl[[cl]]) else 0.0, numeric(1))
  if (all(nk >= as.numeric(min_n_calib))) return(prob_mat)
  pr <- (nk + 0.5) / (max(1.0, sum(nk)) + 0.5 * length(cn))
  ak <- ifelse(nk < as.numeric(min_n_calib), 2.0 * (as.numeric(min_n_calib) - nk) / as.numeric(min_n_calib), 0.0)
  out <- prob_mat
  for (j in seq_along(cn)) if (ak[j] > 0.0) out[, j] <- (nk[j] * prob_mat[, j] + ak[j] * pr[j]) / max(nk[j] + ak[j], 1e-8)
  rs <- rowSums(out, na.rm = TRUE); ok <- is.finite(rs) & rs > 1e-12
  if (any(ok)) out[ok, ] <- sweep(out[ok, , drop = FALSE], 1L, rs[ok], "/")
  out
}

.fold_spatial_balance_wts <- function(df_all, active_bool) {
  if (!exists("compute_spatial_balance_weights", mode = "function")) return(rep(1.0, nrow(df_all)))
  d_sb <- df_all; d_sb$authority_weight <- NULL
  if (all(c("active_mask", "recompute") %in% names(formals(compute_spatial_balance_weights)))) return(compute_spatial_balance_weights(d_sb, active_mask = active_bool, recompute = TRUE))
  compute_spatial_balance_weights(d_sb)
}

.fit_rf_trio <- function(X_e, X_a, y_f, w_f, n_trees, mns, max_d, sf, threads) {
  if (requireNamespace("ranger", quietly = TRUE)) {
    de <- data.frame(y = y_f, as.data.frame(X_e), check.names = FALSE); da <- data.frame(y = y_f, as.data.frame(X_a), check.names = FALSE)
    r1 <- ranger::ranger(y ~ ., data = de, num.trees = n_trees, importance = "permutation", probability = FALSE, min.node.size = mns, max.depth = max_d, sample.fraction = sf, case.weights = w_f, num.threads = threads, seed = 42L)
    r2 <- ranger::ranger(y ~ ., data = da, num.trees = n_trees, importance = "permutation", probability = FALSE, min.node.size = mns, max.depth = max_d, sample.fraction = sf, case.weights = w_f, num.threads = threads, seed = 42L)
    r3 <- ranger::ranger(y ~ ., data = de, num.trees = n_trees, probability = TRUE, min.node.size = mns, max.depth = max_d, sample.fraction = sf, case.weights = w_f, num.threads = threads, seed = 42L)
    list(rf_nonspatial = r1, rf_spatial = r2, rf_prob_anchor = r3, oob_ns = 1.0 - r1$prediction.error, oob_sp = 1.0 - r2$prediction.error)
  } else {
    r1 <- randomForest::randomForest(x = X_e, y = y_f, ntree = n_trees, importance = TRUE)
    r2 <- randomForest::randomForest(x = X_a, y = y_f, ntree = n_trees, importance = TRUE)
    list(rf_nonspatial = r1, rf_spatial = r2, rf_prob_anchor = r1, oob_ns = 1.0 - r1$err.rate[n_trees, "OOB"], oob_sp = 1.0 - r2$err.rate[n_trees, "OOB"])
  }
}

.fit_predict_fold_rf <- function(X_tr_df, y_tr, X_te_df, w_tr, n_trees, mns, max_d, sf, threads, seed) {
  if (requireNamespace("ranger", quietly = TRUE)) {
    fit <- ranger::ranger(y ~ ., data = data.frame(y = y_tr, X_tr_df, check.names = FALSE), num.trees = max(100L, floor(n_trees / 2L)), probability = TRUE, min.node.size = mns, max.depth = max_d, sample.fraction = sf, case.weights = w_tr, num.threads = threads, seed = seed)
    pr <- stats::predict(fit, data = X_te_df, num.threads = threads)$predictions
    if (!is.matrix(pr)) pr <- matrix(pr, nrow = nrow(X_te_df), dimnames = list(NULL, names(pr)))
    return(pr)
  }
  set.seed(seed); fit <- randomForest::randomForest(x = X_tr_df, y = y_tr, ntree = max(100L, floor(n_trees / 2L)))
  stats::predict(fit, newdata = X_te_df, type = "prob")
}

.resolve_rf_source_tags <- function(n, anchor_idx, heldout_anc, elig_coord, is_open_hyb, is_het_type, is_q_type, has_geo) {
  src <- rep("Passive_Projection_Bronze", n); src[!elig_coord] <- "Passive_Projection_High_Uncertainty"
  src[is_q_type] <- "Passive_Type_Uncertain_Coord"; src[is_het_type] <- "Passive_Heterotypic_Type_Audit"
  src[is_open_hyb] <- "Passive_Open_or_Hybrid_Concept"; src[!has_geo] <- "Passive_Ungeoreferenced"
  if (length(anchor_idx) > 0L) { src[anchor_idx] <- "Spatial_Block_OOF_Anchor"; fb <- anchor_idx[!heldout_anc]; if (length(fb) > 0L) src[fb] <- "Fallback_OOB_Anchor" }
  src
}

.build_varpart_df <- function(imp_df, varpart_obj = NULL, include_adj_r2 = TRUE) {
  gn <- c("[A] Fine-Grain Edaphic & Outcrop Lithology", "[B] Macro-Climate", "[C] Spatial MEMs")
  gs <- vapply(gn, function(g) sum(imp_df$PositiveImportance[imp_df$Group == g], na.rm = TRUE), numeric(1))
  ar2 <- rep(NA_real_, 3L)
  if (!is.null(varpart_obj) && !is.null(varpart_obj$part$indfract)) {
    ac <- intersect(c("Adj.R.square", "Adj.R.sq"), colnames(varpart_obj$part$indfract))[1L]
    if (!is.na(ac)) { ar2 <- round(as.numeric(varpart_obj$part$indfract[1:3, ac]), 4); vp <- pmax(0.0, as.numeric(varpart_obj$part$indfract[1:3, ac])); if (sum(vp, na.rm = TRUE) > 0.0) { gs <- vp; names(gs) <- gn } }
  }
  tp <- sum(gs, na.rm = TRUE)
  vpct <- if (tp <= 0.0) round(rep(100.0 / 3.0, 3L), 2) else { v <- round((as.numeric(gs) / tp) * 100.0, 2); v[which.max(v)] <- round(v[which.max(v)] + (100.0 - sum(v)), 2); v }
  if (isTRUE(include_adj_r2)) {
    data.frame(Component = gn, Variance_Pct = vpct, Adj_R2 = ar2, MEM_Dominance_Flag = (gn == "[C] Spatial MEMs") & (vpct > 50.0), stringsAsFactors = FALSE)
  } else data.frame(Component = gn, Variance_Pct = vpct, stringsAsFactors = FALSE)
}

.assemble_srf_return <- function(trio, imp_df, vp_df, vp_obj, MEMs, sb_cv, ac_km, preds, p_prob, g_prob, src, fmat, fblk, cscr, cv_st, mi, mp, mcd, post_sel, anc_i, pas_i, extra = list()) {
  c(list(
    rf_spatial = trio$rf_spatial, rf_nonspatial = trio$rf_nonspatial, rf_prob_anchor = trio$rf_prob_anchor,
    importance = imp_df, variance_partition = vp_df, varpart_model = vp_obj, mems = MEMs,
    OOB_Accuracy = trio$oob_ns, OOB_Error = 1.0 - trio$oob_ns, NonSpatial_OOB_Accuracy = trio$oob_ns,
    Spatial_MEM_OOB_Accuracy = trio$oob_sp, Spatial_Block_CV_Accuracy = sb_cv, spatial_autocorr_range_km = ac_km,
    rf_predicted_taxon = preds, rf_predicted_prob = round(p_prob, 4), rf_given_taxon_prob = round(g_prob, 4),
    rf_prediction_source = src, oof_prob_matrix = fmat, full_prob_mat = fmat, block_cv_acc = sb_cv,
    predicted_taxon = preds, predicted_prob = round(p_prob, 4), spatial_block_id = fblk,
    collinearity_screening = cscr, cv_mean_fold_predictors = round(cv_st$mean_p, 2), cv_fold_predictor_jaccard = round(cv_st$jaccard, 4),
    residual_moran_i = round(mi, 4), residual_moran_p = round(mp, 4), env_mahal_d2 = mcd$env_mahal_d2,
    env_mcd_outlier_flag = mcd$env_mcd_outlier_flag, env_outlier_rescued_azonal = mcd$env_outlier_rescued_azonal,
    env_outlier_screen = mcd, is_post_selection_ecotype = post_sel, anchor_idx = anc_i, passive_idx = pas_i
  ), extra)
}

compute_spatial_rf_mems <- function(df, env_vars, n_trees = 500, intraspecific_mode = FALSE, vif_res = NULL, env_dir = NULL, demo = FALSE) {
  message("Generating Moran's Eigenvector Maps (MEM1..MEM5) & fitting Authority-Anchored Spatial Random Forests...")
  env_vars_in <- env_vars
  cfg <- if (exists("load_step05_macroecology_config", mode = "function")) load_step05_macroecology_config() else list(cor_threshold = 0.70, vif_threshold = 5.0, min_node_size = 12L, max_depth = 6L, sample_fraction = 0.70, ecotype_min_fold_jaccard = 0.60)
  mns_cfg <- if (!is.null(cfg$min_node_size) && is.finite(cfg$min_node_size)) as.integer(round(cfg$min_node_size)) else 12L
  md_cfg <- if (!is.null(cfg$max_depth) && is.finite(cfg$max_depth)) max(1L, as.integer(round(cfg$max_depth))) else 6L
  sf_cfg <- if (!is.null(cfg$sample_fraction) && is.finite(cfg$sample_fraction)) min(1.0, max(0.05, as.numeric(cfg$sample_fraction))) else 0.70
  eco_min_j <- if (!is.null(cfg$ecotype_min_fold_jaccard) && is.finite(cfg$ecotype_min_fold_jaccard)) as.numeric(cfg$ecotype_min_fold_jaccard) else 0.60
  .cmns <- function(m) max(2L, min(mns_cfg, max(2L, as.integer(floor(as.numeric(m) / 3.0)))))
  .csf <- function(m) max(sf_cfg, min(1.0, 2.0 / max(1.0, as.numeric(m))))

  os <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(42L); on.exit({ if (!is.null(os)) assign(".Random.seed", os, envir = .GlobalEnv) }, add = TRUE)
  if (!("analysis_group" %in% names(df)) || all(is.na(df$analysis_group) | !nzchar(trimws(as.character(df$analysis_group))))) df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
  df_in <- df
  if (exists("flag_coordinate_uncertainty", mode = "function") && !("is_training_eligible_coord" %in% names(df))) df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0)
  coords <- as.matrix(df[, c("longitude", "latitude")]); n <- nrow(coords); has_geo <- is.finite(coords[, 1L]) & is.finite(coords[, 2L])
  MEMs <- .compute_fast_mems(coords, max_mems = 5L, try_adespatial = TRUE)
  if (exists("mask_micro_scale_predictors_by_uncertainty", mode = "function")) df <- mask_micro_scale_predictors_by_uncertainty(df, max_uncertainty_m = 250.0)
  if (exists("sanitize_disabled_polaris_columns", mode = "function")) df <- sanitize_disabled_polaris_columns(df)
  env_vars <- setdiff(intersect(env_vars, names(df)), NON_ENM_RAW_GEOCHEM_HYDRO_COLS); env_vars <- env_vars[!grepl("^polaris_", env_vars, ignore.case = TRUE)]
  X_env_raw <- as.matrix(df[, env_vars, drop = FALSE]); storage.mode(X_env_raw) <- "double"; y <- factor(df$analysis_group)

  fe <- if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) is_fine_edaphic_type_or_anchor_eligible(df, 10000.0, 50000.0) else list(edaphic_ok = rep(TRUE, n), clim_ok = rep(TRUE, n), is_open_hyb = rep(FALSE, n), is_quarantined_type = rep(FALSE, n))
  elig_c <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, n); elig_c[is.na(elig_c)] <- TRUE; elig_c <- elig_c & fe$edaphic_ok & has_geo
  elig_e <- if ("is_edaphic_eligible_coord" %in% names(df)) as.logical(df$is_edaphic_eligible_coord) else if ("coordinateUncertaintyInMeters" %in% names(df)) { u2 <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters)); is.finite(u2) & u2 > 0.0 & u2 <= 2000.0 } else elig_c
  elig_e[is.na(elig_e)] <- FALSE; elig_e <- elig_e & fe$edaphic_ok
  det_t <- if (exists(".standardize_determiner_tier_vec", mode = "function")) .standardize_determiner_tier_vec(df) else if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_2_Silver", n)
  t_inf <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else { tf <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, n); tf[is.na(tf)] <- FALSE; list(is_type_specimen = tf) }
  is_prom <- if ("is_promoted_bronze" %in% names(df)) { b <- as.logical(df$is_promoted_bronze); b[is.na(b)] <- FALSE; b } else rep(FALSE, n)
  cl_n <- if ("is_confident_label_noise" %in% names(df)) { b <- as.logical(df$is_confident_label_noise); b[is.na(b)] <- FALSE; b } else rep(FALSE, n)
  is_oh <- as.logical(fe$is_open_hyb); is_qt <- as.logical(fe$is_quarantined_type)
  is_ht <- rep(FALSE, n)
  if ("is_heterotypic_type_audit" %in% names(df)) is_ht <- is_ht | (!is.na(df$is_heterotypic_type_audit) & as.logical(df$is_heterotypic_type_audit))
  if ("anchor_tier_source" %in% names(df)) is_ht <- is_ht | (!is.na(df$anchor_tier_source) & as.character(df$anchor_tier_source) == "Passive_Heterotypic_Type_Audit")
  is_tc <- t_inf$is_type_specimen & !is_qt & !is_ht
  is_gs <- ((det_t %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_prom | is_tc) & (!cl_n | is_tc) & !is_oh & !is_ht & !is_qt
  if (exists(".is_step05_extralimital_non_type", mode = "function")) is_gs <- is_gs & !.is_step05_extralimital_non_type(df, y, is_gs & elig_c & !is.na(y), is_tc)
  elig_coord <- if (any(grepl("^(fine_|topo_|ilr_)", env_vars, ignore.case = TRUE)) && sum(is_gs & elig_c & elig_e & !is.na(y)) >= 6L) (elig_c & elig_e) else elig_c
  miss_g <- if (isTRUE(intraspecific_mode) && "is_missing_gmm_cluster" %in% names(df)) { m <- as.logical(df$is_missing_gmm_cluster); m[is.na(m)] <- FALSE; m } else rep(FALSE, n)
  min_k <- min(2L, max(1L, length(unique(stats::na.omit(y)))))
  anc_m <- is_gs & elig_coord & !is.na(y) & !miss_g
  if (sum(anc_m) < 4L || length(unique(y[anc_m])) < min_k) anc_m <- is_gs & elig_c & !is.na(y) & !miss_g
  if (sum(anc_m) < 4L || length(unique(y[anc_m])) < min_k) anc_m <- is_gs & has_geo & !is.na(y)
  if (!any(is_gs)) { cfb <- !cl_n & !is_oh & !is_ht & !is_qt & has_geo & !is.na(y); anc_m <- cfb & elig_coord & !miss_g; if (sum(anc_m) < min(2L, n)) anc_m <- cfb & elig_c; if (sum(anc_m) == 0L) anc_m <- has_geo & !is.na(y) }
  anc_i <- which(anc_m); pas_i <- setdiff(seq_len(n), anc_i); y_anc <- droplevels(y[anc_i])

  .imp_m <- function(Mt, Mr) { for (jj in seq_len(ncol(Mt))) { nj <- !is.finite(Mt[, jj]); if (any(nj)) { mj <- stats::median(Mr[is.finite(Mr[, jj]), jj], na.rm = TRUE); if (!is.finite(mj)) mj <- stats::median(Mt[is.finite(Mt[, jj]), jj], na.rm = TRUE); Mt[nj, jj] <- if (is.finite(mj)) mj else 0.0 } }; Mt }
  if (ncol(X_env_raw) > 1L && length(anc_i) >= 8L) {
    mf <- max(3L, as.integer(ceiling(0.50 * length(anc_i))))
    fok <- colnames(X_env_raw)[vapply(seq_len(ncol(X_env_raw)), function(jj) { vj <- X_env_raw[anc_i, jj]; fj <- is.finite(vj); sum(fj) >= mf && is.finite(stats::sd(vj[fj])) && stats::sd(vj[fj]) > 1e-8 }, logical(1))]
    if (length(fok) >= 1L && length(fok) < ncol(X_env_raw)) {
      env_vars <- fok; X_env_raw <- X_env_raw[, env_vars, drop = FALSE]
      if (is.list(vif_res)) { if (!is.null(vif_res$kept_vars)) vif_res$kept_vars <- intersect(vif_res$kept_vars, env_vars); if (!is.null(vif_res$selected_variables)) vif_res$selected_variables <- intersect(vif_res$selected_variables, env_vars) }
    }
  }
  X_cov_ok <- X_env_raw; cv_cand <- colnames(X_cov_ok)
  dom_sp <- if (exists(".split_step05_cv_domain_vars", mode = "function")) .split_step05_cv_domain_vars(cv_cand) else {
    cl_c <- intersect(grep("^(bio[0-9]*_|summer_moisture_|maritime_continentality_|delta_doy|doy|hopkins_)", cv_cand, value = TRUE, ignore.case = TRUE), cv_cand)
    tp_c <- setdiff(intersect(grep("^(topo_|elevation($|_)|elev_)", cv_cand, value = TRUE, ignore.case = TRUE), cv_cand), cl_c)
    list(clim_vars = cl_c, topo_vars = tp_c, edaph_vars = setdiff(cv_cand, unique(c(cl_c, tp_c))))
  }
  ag_vif <- resolve_active_anchor_groups(df, intraspecific_mode = intraspecific_mode)
  cscr <- if (is.list(vif_res) && (!is.null(vif_res$kept_vars) || !is.null(vif_res$selected_variables))) {
    vif_res
  } else if (exists("screen_collinear_predictors_by_block", mode = "function")) {
    screen_collinear_predictors_by_block(X_env_raw[anc_i, , drop = FALSE], cor.threshold = cfg$cor_threshold, vif.threshold = cfg$vif_threshold, anchor_groups = ag_vif[anc_i])
  } else if (exists("screen_collinear_predictors", mode = "function")) {
    screen_collinear_predictors(X_env_raw[anc_i, , drop = FALSE], cor.threshold = cfg$cor_threshold, vif.threshold = cfg$vif_threshold, anchor_groups = ag_vif[anc_i])
  } else {
    list(selected_variables = colnames(X_env_raw), kept_vars = colnames(X_env_raw), pruned_variables = character(0))
  }
  kv <- if (!is.null(cscr$kept_vars)) cscr$kept_vars else cscr$selected_variables
  if (length(intersect(colnames(X_env_raw), kv)) >= 1L) { env_vars <- intersect(colnames(X_env_raw), unique(kv)); X_env_raw <- X_env_raw[, env_vars, drop = FALSE] }
  epv <- enforce_epv_predictor_budget(df = X_env_raw[anc_i, , drop = FALSE], selected_vars = env_vars, n_anchors = length(anc_i), anchor_groups = ag_vif[anc_i], n_mems_planned = 0L, intraspecific_mode = intraspecific_mode, quiet = FALSE)
  if (exists("cap_diagnostic_mem_matrix", mode = "function")) MEMs <- cap_diagnostic_mem_matrix(MEMs, n_anchors = length(anc_i), min_env = 4L)
  if (length(intersect(as.character(epv$selected_vars), colnames(X_env_raw))) >= 1L) { env_vars <- intersect(as.character(epv$selected_vars), colnames(X_env_raw)); X_env_raw <- X_env_raw[, env_vars, drop = FALSE] }
  if (is.list(cscr)) {
    cscr$kept_vars <- env_vars; cscr$selected_variables <- env_vars; cscr$selected_vars <- env_vars
    for (fl in c("edaphic_variables", "selected_edaphic", "climatic_variables", "selected_climatic", "topo_variables", "selected_topo")) if (!is.null(cscr[[fl]])) cscr[[fl]] <- intersect(cscr[[fl]], env_vars)
    cscr$epv_budget <- list(n_anchors = length(anc_i), total_epv_budget = epv$total_budget, env_cap = length(env_vars), mem_cap = if (!is.null(MEMs)) ncol(MEMs) else 0L)
  }

  X_env <- .imp_m(X_env_raw, X_env_raw[anc_i, , drop = FALSE]); X_all <- cbind(X_env, MEMs)
  mcd <- if (exists("screen_within_taxon_env_mcd_outliers", mode = "function")) {
    screen_within_taxon_env_mcd_outliers(df, X_env, y, anc_m, is_type_vec = is_tc, p_cutoff = 0.005, min_n_group = 12L, seed = 42L)
  } else {
    list(env_mahal_d2 = rep(NA_real_, n), env_mcd_outlier_flag = rep(FALSE, n), env_outlier_rescued_azonal = rep(FALSE, n), raw_mcd_exceed_flag = rep(FALSE, n), is_type_specimen = t_inf$is_type_specimen)
  }
  if (any(anc_m & mcd$env_mcd_outlier_flag)) {
    cm <- anc_m & !mcd$env_mcd_outlier_flag; yc <- as.character(y); pt <- table(yc[anc_m])
    for (ck in names(pt[pt >= 2L])) {
      kk <- which(cm & !is.na(yc) & (yc == ck))
      if (length(kk) < 2L) {
        dk <- which(anc_m & !cm & !is.na(yc) & (yc == ck))
        if (length(dk) > 0L) {
          ts <- ifelse(is_tc[dk], 4.0, ifelse(det_t[dk] == "Tier_1_Gold", 3.0, ifelse(det_t[dk] == "Tier_2_Silver", 2.0, ifelse(is_prom[dk], 1.5, 1.0))))
          d2 <- mcd$env_mahal_d2[dk]; d2[!is.finite(d2)] <- Inf; cm[dk[order(-ts, d2, dk)[seq_len(min(2L - length(kk), length(dk)))]]] <- TRUE
        }
      }
    }
    if (sum(cm) >= 4L && !identical(cm, anc_m)) { anc_m <- cm; anc_i <- which(anc_m); pas_i <- setdiff(seq_len(n), anc_i); y_anc <- droplevels(y[anc_i]); X_env <- .imp_m(X_env_raw, X_env_raw[anc_i, , drop = FALSE]); X_all <- cbind(X_env, MEMs) }
  }
  adir <- if (exists("env_plantmatics_or_packera", mode = "function")) trimws(env_plantmatics_or_packera("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR")) else ""
  if (!isTRUE(demo) && nzchar(adir) && exists("write_env_outlier_quarantine_audit", mode = "function")) tryCatch(write_env_outlier_quarantine_audit(df, mcd, file.path(adir, "outputs", "reports", "environmental_outlier_quarantine_audit.csv"), final_anchor_mask = anc_m), error = function(e) NULL)

  .cv_st <- function(fs, fb) {
    ne <- Filter(function(v) length(v) > 0L, fs); if (length(ne) == 0L) return(list(mean_p = as.numeric(fb), jaccard = 1.0))
    jc <- if (length(ne) < 2L) 1.0 else { pv <- numeric(0); for (a in seq_len(length(ne) - 1L)) for (b in seq(a + 1L, length(ne))) { u <- length(union(ne[[a]], ne[[b]])); pv <- c(pv, if (u > 0L) length(intersect(ne[[a]], ne[[b]])) / u else 1.0) }; mean(pv) }
    list(mean_p = mean(vapply(ne, length, integer(1))), jaccard = jc)
  }
  auth_wts <- ifelse(anc_m, pmax(resolve_ecological_authority_weights(df), 0.25), 0.0)
  rf_th <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, max(1L, as.integer(ceiling(n / 16L)))) else 1L
  fine_ed_v <- intersect(env_vars, grep("^(fine_|ilr_texture_|ssurgo_|solus_|topo_|soil_|sgmc_|is_azonal|lithic_outcrop_|wetland_seepage_)", env_vars, value = TRUE))
  mac_cl_v <- intersect(env_vars, grep("^(bio|delta_doy|doy|summer_moisture_|maritime_continentality_)", env_vars, value = TRUE, ignore.case = TRUE))
  fsets <- list()

  # K = 1 Single-Lineage Presence vs Terrestrial Accessible Background (M)
  if (length(levels(y_anc)) < 2L) {
    stx <- if (length(levels(y_anc)) == 1L) levels(y_anc)[1L] else as.character(stats::na.omit(df$species_standardized))[1L]
    if (is.na(stx) || !nzchar(stx)) stx <- "Packera dubia"
    aci <- estimate_spatial_autocorrelation_block_size(coords[anc_i, , drop = FALSE], X_env[anc_i, , drop = FALSE])
    bg <- .generate_accessible_background_env(coords[anc_i, , drop = FALSE], X_env[anc_i, , drop = FALSE], min(600L, max(120L, 2L * length(anc_i))), 42L, env_dir, demo, aci$range_km)
    cpb <- rbind(coords[anc_i, , drop = FALSE], bg$coords); Xepb <- rbind(X_env[anc_i, , drop = FALSE], bg$X_bg)
    ypb <- factor(c(rep("Presence", length(anc_i)), rep("Background", nrow(bg$X_bg))), levels = c("Presence", "Background"))
    Xapb <- cbind(Xepb, .compute_fast_mems(cpb, max_mems = ncol(MEMs), try_adespatial = FALSE))
    wp <- auth_wts[anc_i] * (if ("edaphic_coord_weight" %in% names(df)) pmax(0.25, as.numeric(df$edaphic_coord_weight[anc_i])) else rep(1.0, length(anc_i)))
    wp <- wp / max(mean(wp, na.rm = TRUE), 1e-8); wpb <- c(wp, rep(sum(wp) / max(1L, nrow(bg$X_bg)), nrow(bg$X_bg)))
    tpb <- .fit_rf_trio(Xepb, Xapb, ypb, wpb, n_trees, .cmns(nrow(Xepb)), md_cfg, .csf(nrow(Xepb)), rf_th)
    fpb <- assign_spatial_grid_blocks(cpb, block_deg = aci$range_deg, k = 5L); oof_p <- rep(NA_real_, length(fpb)); pbo <- as.character(ypb)
    Xpb_c <- rbind(X_cov_ok[anc_i, env_vars, drop = FALSE], bg$X_bg[, env_vars, drop = FALSE]); upf <- sort(unique(fpb))
    cw <- if (length(upf) >= 2L && nrow(Xpb_c) >= 120L && rf_th > 1L && exists("packera_mclapply", mode = "function")) min(rf_th, length(upf)) else 1L; pft <- max(1L, as.integer(floor(rf_th / cw)))
    ev_pb <- function(kf) {
      ti <- which(fpb != kf); te <- which(fpb == kf); if (!(length(ti) >= 8L && length(te) >= 1L && length(unique(pbo[ti])) == 2L)) return(NULL)
      fs <- screen_predictors_within_cv_fold(as.data.frame(Xpb_c[ti, , drop = FALSE], check.names = FALSE), env_vars, intersect(dom_sp$clim_vars, env_vars), intersect(dom_sp$edaph_vars, env_vars), intersect(dom_sp$topo_vars, env_vars), factor(pbo[ti]), cfg$cor_threshold, cfg$vif_threshold, intraspecific_mode)
      pr <- .fit_predict_fold_rf(fs$X_tr_imp, droplevels(ypb[ti]), apply_fold_imputation(as.data.frame(Xpb_c[te, , drop = FALSE], check.names = FALSE), fs$selected_vars, fs$train_medians), wpb[ti], n_trees, .cmns(length(ti)), md_cfg, .csf(length(ti)), pft, 42L + as.integer(kf))
      list(key = paste0("pb_", kf), selected_vars = fs$selected_vars, te_i = te, pres_p = if (is.matrix(pr) && "Presence" %in% colnames(pr)) pr[, "Presence"] else rep(NA_real_, length(te)))
    }
    pbr <- if (cw > 1L) packera_mclapply(upf, ev_pb, mc.cores = cw) else lapply(upf, ev_pb)
    for (fr in pbr) if (!is.null(fr)) { fsets[[fr$key]] <- fr$selected_vars; oof_p[fr$te_i] <- fr$pres_p }
    ho_pb <- is.finite(oof_p)
    if (any(!ho_pb) && is.matrix(tpb$rf_prob_anchor$predictions) && "Presence" %in% colnames(tpb$rf_prob_anchor$predictions)) oof_p[!ho_pb] <- tpb$rf_prob_anchor$predictions[!ho_pb, "Presence"]
    ok_p <- which(ho_pb); sb_cv <- if (length(ok_p) > 0L) sum(wpb[ok_p] * (ifelse(oof_p[ok_p] >= 0.50, "Presence", "Background") == pbo[ok_p]), na.rm = TRUE) / sum(wpb[ok_p]) else tpb$oob_ns
    if (!is.finite(sb_cv)) sb_cv <- tpb$oob_ns
    suit <- rep(NA_real_, n); suit[anc_i] <- oof_p[seq_along(anc_i)]
    if (length(pas_i) > 0L) { pp <- stats::predict(tpb$rf_prob_anchor, data = as.data.frame(X_env[pas_i, , drop = FALSE]), num.threads = rf_th)$predictions; suit[pas_i] <- if (is.matrix(pp) && "Presence" %in% colnames(pp)) pp[, "Presence"] else 0.75 }
    suit[!is.finite(suit)] <- 0.75; fmat <- matrix(round(suit, 4), nrow = n, ncol = 1L, dimnames = list(NULL, stx))
    fblk <- rep(NA_character_, n); fblk[anc_i] <- as.character(fpb[seq_along(anc_i)])
    imp_df <- assemble_mem_decoupled_importance_df(tpb$rf_nonspatial, tpb$rf_spatial, env_vars)
    ibe_df <- if (exists("compute_isolation_by_environment", mode = "function")) tryCatch(compute_isolation_by_environment(df[anc_i, , drop = FALSE], fine_ed_v, mac_cl_v, MEMs[anc_i, , drop = FALSE]), error = function(e) NULL) else NULL
    return(.assemble_srf_return(tpb, imp_df, .build_varpart_df(imp_df, NULL, FALSE), NULL, MEMs, sb_cv, aci$range_km, rep(stx, n), suit, suit,
      .resolve_rf_source_tags(n, anc_i, ho_pb[seq_along(anc_i)], elig_coord, is_oh, is_ht, is_qt, has_geo), fmat, fblk, cscr, .cv_st(fsets, length(env_vars)), 0.0, 1.0, mcd, FALSE, anc_i, pas_i,
      list(single_lineage_mode = TRUE, ibe_summary = ibe_df, accessible_background = bg)))
  }

  # K >= 2 Multi-Taxon / Ecotype Spatial Block-CV RF
  sb_full <- .fold_spatial_balance_wts(df, anc_m); inv_w <- rep(0.0, n); inv_w[anc_i] <- compute_inverse_class_frequency_weights(y_anc)
  cwts <- auth_wts * inv_w * sb_full; cwts[anc_i] <- pmin(4.0, pmax(0.25, cwts[anc_i] / max(mean(cwts[anc_i], na.rm = TRUE), 1e-8)))
  cwts[anc_i] <- cwts[anc_i] / mean(cwts[anc_i], na.rm = TRUE); cwts[!anc_m] <- 0.0
  X_ea <- X_env[anc_i, , drop = FALSE]; X_aa <- X_all[anc_i, , drop = FALSE]; w_anc <- cwts[anc_i]
  tm <- .fit_rf_trio(X_ea, X_aa, y_anc, w_anc, n_trees, .cmns(nrow(X_ea)), md_cfg, .csf(nrow(X_ea)), rf_th)
  aci <- estimate_spatial_autocorrelation_block_size(coords[anc_i, , drop = FALSE], X_ea); n_anc <- length(anc_i); fids <- NULL
  if (requireNamespace("blockCV", quietly = TRUE) && requireNamespace("sf", quietly = TRUE) && n_anc >= 10L) {
    fids <- tryCatch({
      psf <- sf::st_as_sf(data.frame(longitude = coords[anc_i, 1L], latitude = coords[anc_i, 2L], y = y_anc), coords = c("longitude", "latitude"), crs = 4326)
      tryCatch(suppressMessages(suppressWarnings(blockCV::cv_spatial(x = psf, column = "y", size = aci$range_m, k = 5, selection = "random", iteration = 10, extend = 0.5, seed = 42L, progress = FALSE, report = FALSE, plot = FALSE)))$folds_ids,
               error = function(e) suppressMessages(suppressWarnings(blockCV::cv_spatial(x = psf, column = "y", k = 5, selection = "random", iteration = 10, extend = 0.5, seed = 42L, progress = FALSE, report = FALSE, plot = FALSE)))$folds_ids)
    }, error = function(e) NULL)
  }
  if (is.null(fids) || length(unique(stats::na.omit(fids))) < 2L) fids <- assign_spatial_grid_blocks(coords[anc_i, , drop = FALSE], aci$range_deg, min(5L, max(2L, floor(n_anc / 2L))))
  ctab <- table(y_anc); mcls <- names(ctab[ctab >= 2L])
  if (length(mcls) >= 2L) {
    uf <- sort(unique(stats::na.omit(fids)))
    if (length(uf) < 2L || any(vapply(uf, function(fid) !all(mcls %in% as.character(y_anc[fids != fid])), logical(1)))) {
      kst <- min(5L, max(2L, min(as.integer(ctab[mcls])))); fids <- integer(n_anc)
      for (ci in seq_along(levels(y_anc))) {
        ic <- which(as.character(y_anc) == levels(y_anc)[ci])
        if (length(ic) == 1L) fids[ic] <- ((ci - 1L) %% kst) + 1L else if (length(ic) > 1L) { so <- order(coords[anc_i[ic], 1L] + 0.618034 * coords[anc_i[ic], 2L]); fids[ic[so]] <- ((ceiling(seq_along(so) * kst / length(so)) + ci - 2L) %% kst) + 1L }
      }
    }
  }

  clv <- levels(y_anc); oof_m <- matrix(NA_real_, nrow = n_anc, ncol = length(clv), dimnames = list(NULL, clv))
  oof_p <- rep(NA_character_, n_anc); oof_r <- rep(NA_character_, n_anc); psel <- is_post_selection_ecotype_groups(y_anc, intraspecific_mode, df)
  Xca <- as.data.frame(X_cov_ok[anc_i, cv_cand, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE); ufs <- sort(unique(stats::na.omit(fids)))
  cvw <- if (length(ufs) >= 2L && n_anc >= 80L && rf_th > 1L && exists("packera_mclapply", mode = "function")) min(rf_th, length(ufs)) else 1L; fth <- max(1L, as.integer(floor(rf_th / cvw)))
  ev_f <- function(fid) {
    tr <- which(fids != fid); te <- which(fids == fid); if (!(length(tr) >= 4L && length(te) >= 1L && length(unique(y_anc[tr])) >= 2L)) return(NULL)
    yt <- droplevels(y_anc[tr])
    fs <- screen_predictors_within_cv_fold(Xca[tr, , drop = FALSE], cv_cand, dom_sp$clim_vars, dom_sp$edaph_vars, dom_sp$topo_vars, yt, cfg$cor_threshold, cfg$vif_threshold, intraspecific_mode)
    Xte <- apply_fold_imputation(Xca[te, , drop = FALSE], fs$selected_vars, fs$train_medians)
    or_f <- if (isTRUE(psel)) compute_fold_reclustered_ecotype_ref(fs$X_tr_imp, yt, Xte, seed = 42L + as.integer(fid)) else as.character(y_anc[te])
    fam <- rep(FALSE, n); fam[anc_i[tr]] <- TRUE
    wtf <- auth_wts[anc_i[tr]] * .fold_spatial_balance_wts(df, fam)[anc_i[tr]] * compute_inverse_class_frequency_weights(yt)
    wtf <- pmin(4.0, pmax(0.25, wtf / max(mean(wtf, na.rm = TRUE), 1e-8)))
    prm <- .calibrate_small_n_rf_probs(.fit_predict_fold_rf(fs$X_tr_imp, yt, Xte, wtf, n_trees, .cmns(length(tr)), md_cfg, .csf(length(tr)), fth, 42L + as.integer(fid)), yt, 10L)
    list(key = as.character(fid), selected_vars = fs$selected_vars, te = te, oof_ref = or_f, pr_mat = prm, preds = colnames(prm)[max.col(prm, ties.method = "first")])
  }
  fres <- if (cvw > 1L) packera_mclapply(ufs, ev_f, mc.cores = cvw) else lapply(ufs, ev_f)
  for (fr in fres) if (!is.null(fr)) { fsets[[fr$key]] <- fr$selected_vars; oof_r[fr$te] <- fr$oof_ref; oof_m[fr$te, ] <- 0.0; for (cl in intersect(colnames(fr$pr_mat), clv)) oof_m[fr$te, cl] <- fr$pr_mat[, cl]; oof_p[fr$te] <- fr$preds }
  ho <- !is.na(oof_p) & !is.na(oof_r); mo <- which(is.na(oof_p))
  if (length(mo) > 0L && requireNamespace("ranger", quietly = TRUE) && is.matrix(tm$rf_prob_anchor$predictions)) {
    ob <- .calibrate_small_n_rf_probs(tm$rf_prob_anchor$predictions, y_anc, 10L)
    for (cl in intersect(colnames(ob), clv)) oof_m[mo, cl] <- ob[mo, cl]
    oof_p[mo] <- colnames(ob)[max.col(ob[mo, , drop = FALSE], ties.method = "first")]
  }
  okc <- which(ho); sb_cv <- if (length(okc) > 0L) sum(w_anc[okc] * (oof_p[okc] == oof_r[okc])) / sum(w_anc[okc]) else tm$oob_ns
  if (!is.finite(sb_cv)) sb_cv <- tm$oob_ns

  fmat <- matrix(NA_real_, nrow = n, ncol = length(clv), dimnames = list(NULL, clv)); fmat[anc_i, ] <- oof_m
  fpred <- rep(NA_character_, n); fpred[anc_i] <- oof_p; fblk <- rep(NA_character_, n); fblk[anc_i] <- as.character(fids)
  if (length(pas_i) > 0L) {
    pp <- if (requireNamespace("ranger", quietly = TRUE)) { pr <- stats::predict(tm$rf_prob_anchor, data = as.data.frame(X_env[pas_i, , drop = FALSE]), num.threads = rf_th)$predictions; if (!is.matrix(pr)) matrix(pr, nrow = length(pas_i), dimnames = list(NULL, names(pr))) else pr } else stats::predict(tm$rf_prob_anchor, newdata = X_env[pas_i, , drop = FALSE], type = "prob")
    pp <- .calibrate_small_n_rf_probs(pp, y_anc, 10L); for (cl in intersect(colnames(pp), clv)) fmat[pas_i, cl] <- pp[, cl]; fpred[pas_i] <- colnames(pp)[max.col(pp, ties.method = "first")]
  }
  p_pr <- apply(fmat, 1L, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
  g_pr <- vapply(seq_len(n), function(i) { gi <- as.character(y[i]); if (!is.na(gi) && gi %in% colnames(fmat)) fmat[i, gi] else NA_real_ }, numeric(1))
  imp_df <- assemble_mem_decoupled_importance_df(tm$rf_nonspatial, tm$rf_spatial, env_vars)
  vp_obj <- if (requireNamespace("vegan", quietly = TRUE) && length(fine_ed_v) >= 1L && length(mac_cl_v) >= 1L) tryCatch(vegan::varpart(stats::model.matrix(~ y_anc - 1), as.data.frame(X_ea[, fine_ed_v, drop = FALSE]), as.data.frame(X_ea[, mac_cl_v, drop = FALSE]), as.data.frame(MEMs[anc_i, , drop = FALSE])), error = function(e) NULL) else NULL
  osp <- if (requireNamespace("ranger", quietly = TRUE)) as.character(tm$rf_spatial$predictions) else as.character(tm$rf_spatial$predicted)
  rv <- as.numeric(osp != as.character(y_anc)) - mean(osp != as.character(y_anc), na.rm = TRUE); rv[!is.finite(rv)] <- 0.0; mi <- 0.0; mp <- 1.0
  if (n_anc >= 5L && sum(rv^2) > 1e-10) {
    W <- 1.0 / (.haversine_km_safe(coords[anc_i, , drop = FALSE]) + 5.0); diag(W) <- 0.0
    if (sum(W) > 1e-8) { mi <- as.numeric((n_anc / sum(W)) * (t(rv) %*% W %*% rv) / sum(rv^2)); mp <- 2.0 * stats::pnorm(-abs((mi + 1.0 / (n_anc - 1.0)) / max(sqrt(2.0 / max(n_anc, 4L)), 1e-6))) }
  }
  cv_st <- .cv_st(fsets, length(env_vars))
  clj <- if (length(okc) > 0L) mean(vapply(unique(c(as.character(oof_r[okc]), as.character(y_anc[okc]))), function(cl) { u <- sum(y_anc[okc] == cl | oof_r[okc] == cl); if (u > 0L) sum(y_anc[okc] == cl & oof_r[okc] == cl) / u else 0.0 }, numeric(1)), na.rm = TRUE) else 0.0
  ej <- max(c(cv_st$jaccard, min(clj, sb_cv)), na.rm = TRUE)
  if (isTRUE(intraspecific_mode) && isTRUE(psel) && (!is.finite(ej) || ej < eco_min_j || !is.finite(cv_st$jaccard) || cv_st$jaccard < 0.35 || !is.finite(sb_cv) || sb_cv < 0.70)) {
    sp_fb <- if ("species_standardized" %in% names(df_in)) stats::na.omit(trimws(as.character(df_in$species_standardized))) else character(0)
    dk1 <- df_in; dk1$analysis_group <- rep(if (length(sp_fb) > 0L && nzchar(sp_fb[1L])) sp_fb[1L] else "Intraspecific_Lineage", nrow(dk1)); dk1$ecotype_regime <- dk1$analysis_group; attr(dk1, "is_post_selection_ecotype") <- FALSE
    k1r <- compute_spatial_rf_mems(dk1, env_vars_in, n_trees, intraspecific_mode, NULL, env_dir, demo); k1r$is_post_selection_ecotype <- FALSE; k1r$ecotype_collapsed_to_unimodal <- TRUE; return(k1r)
  }
  .assemble_srf_return(tm, imp_df, .build_varpart_df(imp_df, vp_obj, TRUE), vp_obj, MEMs, sb_cv, aci$range_km, fpred, p_pr, g_pr,
    .resolve_rf_source_tags(n, anc_i, ho, elig_coord, is_oh, is_ht, is_qt, has_geo), fmat, fblk, cscr, cv_st, mi, mp, mcd, psel, anc_i, pas_i)
}
