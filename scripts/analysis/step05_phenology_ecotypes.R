# ==============================================================================
# Step 05 Module: Ecotype Regime Discovery, Anchor Resolution & IBE Coordinator
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.step05_ohits <- Filter(nzchar, vapply(sys.frames(), function(f) if (!is.null(f$ofile)) as.character(f$ofile) else "", character(1)))
.step05_cfile <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.step05_resolve_cands <- function(fname) {
  cands <- c(file.path("scripts", "analysis", fname), file.path(getwd(), "scripts", "analysis", fname), file.path("analysis", fname), fname)
  if (length(.step05_ohits) > 0L) cands <- c(file.path(dirname(utils::tail(.step05_ohits, 1L)), fname), cands)
  if (length(.step05_cfile) > 0L) cands <- c(file.path(dirname(.step05_cfile[1L]), fname), cands)
  cands
}
for (.m in list(
  c("step05_circular_hopkins_models.R", "select_circular_hopkins_baseline_model"),
  c("step05_config_geo_utils.R",        "load_step05_macroecology_config"),
  c("step05_allochronic_pheno_utils.R", "test_allochronic_bimodality_von_mises"),
  c("step05_type_niche_utils.R",        "compute_type_locality_ecological_control")
)) {
  if (!exists(.m[2L], mode = "function")) {
    .hit <- .step05_resolve_cands(.m[1L])[file.exists(.step05_resolve_cands(.m[1L]))][1L]
    if (!is.na(.hit) && nzchar(.hit)) source(.hit)
  }
}
if (!exists(".to_bool_s05", mode = "function")) .to_bool_s05 <- function(x) if (is.logical(x)) (!is.na(x) & x) else if (is.numeric(x) || is.integer(x)) (!is.na(x) & is.finite(x) & x != 0) else tolower(trimws(as.character(x))) %in% c("true", "1", "t", "yes")

.is_step05_extralimital_non_type <- function(df, group_vec = NULL, clean_base_mask = NULL, is_type_vec = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(logical(0))
  n <- nrow(df)
  is_t <- if (!is.null(is_type_vec) && length(is_type_vec) == n) as.logical(is_type_vec)
    else if (exists(".identify_type_specimen_vec", mode = "function")) as.logical(.identify_type_specimen_vec(df)$is_type_specimen)
    else if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, n)
  is_t[is.na(is_t)] <- FALSE
  us <- if ("usda_county_status" %in% names(df)) tolower(trimws(as.character(df$usda_county_status))) else rep("", n); us[is.na(us)] <- ""
  ex_c <- if ("exsiccatae_taxonomic_conflict" %in% names(df)) .to_bool_s05(df$exsiccatae_taxonomic_conflict) else rep(FALSE, n)
  raw_m <- (!is_t) & ((us %in% c("out_of_range", "natureserve_exotic_or_waif")) | ex_c)
  if (!any(raw_m)) return(raw_m)
  g_v <- if (!is.null(group_vec) && length(group_vec) == n) as.character(group_vec) else if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("ALL", n)
  bm <- if (!is.null(clean_base_mask) && length(clean_base_mask) == n) as.logical(clean_base_mask) else rep(TRUE, n); bm[is.na(bm)] <- FALSE
  out_m <- rep(FALSE, n)
  for (tx in unique(g_v[!is.na(g_v) & nzchar(g_v)])) if (sum(bm & !is.na(g_v) & (g_v == tx) & !raw_m) >= 2L) out_m[!is.na(g_v) & (g_v == tx) & raw_m] <- TRUE
  out_m
}

.step05_anchor_flags <- function(df, min_edaph_n = 6L) {
  n <- nrow(df)
  geo_el <- if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) tryCatch(is_fine_edaphic_type_or_anchor_eligible(df, max_edaphic_unc_m = 10000.0, max_clim_unc_m = 50000.0), error = function(e) NULL) else NULL
  elig <- if ("is_training_eligible_coord" %in% names(df)) { e <- as.logical(df$is_training_eligible_coord); e[is.na(e)] <- TRUE; e } else rep(TRUE, n)
  for (uc in intersect(c("coordinateUncertaintyInMeters", "coordinateUncertainty"), names(df))) { u <- suppressWarnings(as.numeric(df[[uc]])); elig <- elig & (is.na(u) | u <= 10000.0) }
  if (!is.null(geo_el) && !is.null(geo_el$edaphic_ok) && length(geo_el$edaphic_ok) == n && sum(as.logical(geo_el$edaphic_ok), na.rm = TRUE) >= min_edaph_n) elig <- elig & !is.na(geo_el$edaphic_ok) & as.logical(geo_el$edaphic_ok)

  det_t <- if ("determiner_tier" %in% names(df)) ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier)) else rep("Tier_1_Gold", n)
  is_t <- if ("is_type_specimen" %in% names(df)) .to_bool_s05(df$is_type_specimen) else rep(FALSE, n)
  is_prom <- if ("is_promoted_bronze" %in% names(df)) .to_bool_s05(df$is_promoted_bronze) else rep(FALSE, n)
  cl_noise <- if ("is_confident_label_noise" %in% names(df)) .to_bool_s05(df$is_confident_label_noise) else rep(FALSE, n)
  quar <- if ("is_centroid_or_inst_coord" %in% names(df)) .to_bool_s05(df$is_centroid_or_inst_coord) else rep(FALSE, n)
  if (exists("detect_r_centroid_or_inst_quarantine", mode = "function")) { qd <- as.logical(detect_r_centroid_or_inst_quarantine(df, is_type_vec = rep(FALSE, n))$is_quarantined); qd[is.na(qd)] <- FALSE; quar <- quar | qd }
  if (!is.null(geo_el) && !is.null(geo_el$is_quarantined_type) && length(geo_el$is_quarantined_type) == n) quar <- quar | (!is.na(geo_el$is_quarantined_type) & as.logical(geo_el$is_quarantined_type))
  open_h <- if (!is.null(geo_el) && !is.null(geo_el$is_open_hyb) && length(geo_el$is_open_hyb) == n) (!is.na(geo_el$is_open_hyb) & as.logical(geo_el$is_open_hyb))
    else if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) { m <- as.logical(.identify_open_or_hybrid_concept_mask(df)); m[is.na(m)] <- FALSE; m }
    else if ("is_open_or_hybrid_concept" %in% names(df)) .to_bool_s05(df$is_open_or_hybrid_concept) else rep(FALSE, n)
  list(elig = elig, det_t = det_t, is_t = is_t, is_prom = is_prom, cl_noise = cl_noise, quar = quar, open_h = open_h)
}

resolve_analysis_group <- function(df, intraspecific_mode = FALSE, primary_env_vars = NULL, rf_var_weights = NULL, allow_unimodal_null = NULL) {
  if (exists("compute_derived_habitat_predictors", mode = "function")) df <- compute_derived_habitat_predictors(df)
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2L)
  n_rows <- nrow(df); df$is_missing_gmm_cluster <- rep(FALSE, n_rows)
  if (!is_intra || n_rows == 0L) {
    df$analysis_group <- if (n_rows == 0L) character(0) else as.character(df$species_standardized)
    df$ecotype_regime <- df$analysis_group; attr(df, "is_post_selection_ecotype") <- FALSE; return(df)
  }
  base_taxon <- as.character(stats::na.omit(df$species_standardized))
  base_taxon <- if (length(base_taxon) > 0L && nzchar(trimws(base_taxon[1L]))) trimws(base_taxon[1L]) else "Intraspecific_Lineage"
  allow_k1 <- if (is.null(allow_unimodal_null)) TRUE else isTRUE(allow_unimodal_null)
  ret_k1 <- function(msg = FALSE) {
    if (isTRUE(msg)) message(sprintf("[Independent Step 05 ENM] Retaining K=1 unimodal ecological envelope for '%s'.", base_taxon))
    df$analysis_group <- rep(base_taxon, n_rows); df$ecotype_regime <- df$analysis_group
    attr(df, "is_post_selection_ecotype") <- FALSE; df
  }

  exclude_cols <- c(
    "catalogNumber", "plant_individual_id", "latitude", "longitude", "decimalLatitude", "decimalLongitude",
    "coordinateUncertaintyInMeters", "coordinateUncertainty", "buffer_radius_m", "year", "month", "day", "doy",
    "hopkins_adjusted_doy", "hopkins_expected_doy", "hopkins_adjusted_doy_resid", "delta_doy", "delta_doy_resid",
    "pheno_theta", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic", "authority_weight",
    "edaphic_coord_weight", "spatial_balance_weight", "Combined_Prior_Weight", "case_wts", "gmm_cluster",
    "gmm_uncertainty", "mahalanobis_d2", "mahalanobis_pval", "cda_posterior_prob", "cda_posterior_margin",
    "loocv_accuracy", "spatial_block_cv_acc", "spatial_autocorr_range_km", "rf_predicted_prob",
    "rf_given_taxon_prob", "multimodal_concordance", "edaphic_z_score", "climatic_z_score",
    "type_edaphic_z_score", "type_morph_mahal_d2", "self_morph_mahal_d2", "min_congener_mahal_d2",
    "open_set_novelty_p", "n_observed_modalities", "total_capitula", "n_anthesis", "n_fruiting", "n_bud",
    "n_senescent", "reproductive_stage_confidence", "dominant_pappus_ratio", "radiate_head_fraction",
    "pixels_per_mm", "duplicate_count", "is_rock_outcrop", "is_azonal_lithology",
    "azonal_lithology_presence_buffer", "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm",
    "min_bedrock_depth_buffer_cm", "ph_q10_buffer", "ph_q90_buffer", "max_cfvo_buffer_pct",
    "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
    "ultramafic_calcareous_gradient", "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m",
    "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m",
    "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m", "nwpl_indicator_code", "nwpl_wetland_score"
  )

  df_rag <- df
  if (exists("mask_micro_scale_predictors_by_uncertainty", mode = "function")) df_rag <- mask_micro_scale_predictors_by_uncertainty(df_rag, max_uncertainty_m = 250.0)
  if (exists("sanitize_disabled_polaris_columns", mode = "function")) df_rag <- sanitize_disabled_polaris_columns(df_rag)

  fl <- .step05_anchor_flags(df_rag, min_edaph_n = 6L)
  cl_nt <- fl$cl_noise & !fl$is_t
  rag_grp_v <- if ("species_standardized" %in% names(df_rag)) ifelse(!is.na(df_rag$species_standardized) & nzchar(trimws(as.character(df_rag$species_standardized))), as.character(df_rag$species_standardized), base_taxon) else rep(base_taxon, n_rows)
  ext_m <- .is_step05_extralimital_non_type(df_rag, group_vec = rag_grp_v, clean_base_mask = fl$elig, is_type_vec = fl$is_t)
  anc_mask <- fl$elig & (grepl("Tier_0|Tier_1|Gold|Tier_2|Silver", fl$det_t, ignore.case = TRUE) | fl$is_t) & !fl$quar & !fl$open_h & !cl_nt & !ext_m
  if (sum(anc_mask) < 8L) anc_mask <- fl$elig & !fl$quar & !fl$open_h & !cl_nt
  if (sum(anc_mask) < 8L) anc_mask <- fl$elig & !fl$quar & !fl$open_h
  if (sum(anc_mask) < 6L) anc_mask <- rep(TRUE, n_rows)
  anc_idx <- which(anc_mask)

  cand_env <- if (!is.null(primary_env_vars) && length(primary_env_vars) > 0L) setdiff(intersect(as.character(primary_env_vars), names(df_rag)), exclude_cols)
    else if (exists("discover_candidate_eco_predictors", mode = "function")) setdiff(discover_candidate_eco_predictors(df_rag, include_pheno = FALSE, min_finite = 3L), exclude_cols)
    else character(0)
  if (length(cand_env) == 0L) {
    cand_env <- names(df_rag)[vapply(names(df_rag), function(nm) {
      if (nm %in% exclude_cols || grepl("^polaris_", nm, ignore.case = TRUE) ||
          grepl("^(PC[0-9]+|ECT_|Cauline_|can[0-9]+|sym_|asym_|MEM[0-9]+|basal_|petiole_|blade_|plant_height|involucre_|capitula_)", nm)) return(FALSE)
      is.numeric(df_rag[[nm]]) || is.integer(df_rag[[nm]])
    }, logical(1))]
  }
  cand_env <- cand_env[vapply(cand_env, function(cl) {
    v <- suppressWarnings(as.numeric(df_rag[[cl]])); fin <- is.finite(v); va <- v[anc_idx]; fa <- is.finite(va)
    sum(fin) >= max(3L, ceiling(0.40 * n_rows)) && sum(fa) >= max(4L, ceiling(0.75 * length(anc_idx))) &&
      is.finite(stats::sd(v[fin])) && stats::sd(v[fin]) > 1e-8 && is.finite(stats::sd(va[fa])) && stats::sd(va[fa]) > 1e-8
  }, logical(1))]
  if (length(cand_env) == 0L || n_rows < 6L) return(ret_k1())

  E_raw <- matrix(0.0, nrow = n_rows, ncol = length(cand_env), dimnames = list(NULL, cand_env))
  for (j in seq_along(cand_env)) {
    vj <- suppressWarnings(as.numeric(df_rag[[cand_env[j]]])); va <- vj[anc_idx]
    mj <- stats::median(va[is.finite(va)], na.rm = TRUE)
    if (!is.finite(mj)) mj <- stats::median(vj[is.finite(vj)], na.rm = TRUE)
    if (!is.finite(mj)) mj <- 0.0
    vj[!is.finite(vj)] <- mj; E_raw[, j] <- vj
  }

  pm_keep <- NULL
  if (exists("screen_edaphic_point_mass_columns", mode = "function")) {
    pm_res <- tryCatch(screen_edaphic_point_mass_columns(E_raw, ref_mask = anc_mask, max_tie_frac = 0.25), error = function(e) NULL)
    if (is.character(pm_res)) pm_keep <- intersect(pm_res, cand_env)
    else if (is.logical(pm_res) && length(pm_res) == length(cand_env)) pm_keep <- cand_env[!is.na(pm_res) & pm_res]
    else if (is.numeric(pm_res) && length(pm_res) > 0L) pm_keep <- cand_env[pm_res[pm_res >= 1L & pm_res <= length(cand_env)]]
    else if (is.matrix(pm_res) || is.data.frame(pm_res)) pm_keep <- intersect(colnames(pm_res), cand_env)
    else if (is.list(pm_res)) for (kn in c("kept_cols", "kept_vars", "keep_cols", "kept", "selected_cols", "selected_variables", "cols")) {
      if (!is.null(pm_res[[kn]]) && is.character(pm_res[[kn]])) { pm_keep <- intersect(pm_res[[kn]], cand_env); break }
    }
  }
  if (is.null(pm_keep)) {
    tie_ok <- vapply(seq_along(cand_env), function(jj) {
      va <- round(E_raw[anc_idx, jj], 5); va <- va[is.finite(va)]
      length(va) >= 4L && ((max(tabulate(match(va, unique(va)))) / length(va)) <= 0.25)
    }, logical(1))
    pm_keep <- cand_env[tie_ok]
  }
  if (length(pm_keep) == 0L) return(ret_k1())
  cand_env <- pm_keep; E_raw <- E_raw[, cand_env, drop = FALSE]
  mu_v <- colMeans(E_raw[anc_idx, , drop = FALSE])
  sd_v <- pmax(apply(E_raw[anc_idx, , drop = FALSE], 2L, stats::sd), 1e-6)
  Z_full <- sweep(sweep(E_raw, 2L, mu_v, "-"), 2L, sd_v, "/"); Z_full[!is.finite(Z_full)] <- 0.0

  assign_eco_domain <- function(vn) {
    vl <- tolower(vn)
    if (grepl("bedrock|outcrop|lithic|cfvo|fragvol|resdept", vl)) return("Lithic_Outcrop")
    if (grepl("ph|cec|ilr_texture|sand|clay|silt|aws|solus|ssurgo|soil|wetland|seepage", vl)) return("Rhizosphere_Pedology")
    if (grepl("topo_|slope|tpi|twi|hli|elev", vl)) return("Microtopography")
    "Hydrothermal_Maritime_Climate"
  }
  dom_full <- vapply(cand_env, assign_eco_domain, character(1))
  prio_mech <- c("lithic_outcrop_index", "wetland_seepage_index", "fine_bedrock_depth_cm", "fine_soil_ph", "ilr_texture_1",
                 "summer_moisture_deficit_idx", "maritime_continentality_idx", "topo_slope_deg", "topo_twi_30m", "bio4_temp_seasonality", "bio15_precip_seasonality")
  keep_idx <- integer(0)
  for (dm in sort(unique(dom_full))) {
    dm_idx <- which(dom_full == dm); if (length(dm_idx) <= 1L) { keep_idx <- c(keep_idx, dm_idx); next }
    sp_sc <- vapply(dm_idx, function(jj) {
      q <- stats::quantile(Z_full[anc_idx, jj], c(0.10, 0.25, 0.75, 0.90), na.rm = TRUE, names = FALSE)
      0.5 * (q[3L] - q[2L]) + 0.5 * (q[4L] - q[1L])
    }, numeric(1))
    sel_dm <- integer(0)
    for (jj in dm_idx[order(-as.integer(cand_env[dm_idx] %in% prio_mech), -sp_sc, dm_idx)]) {
      if (length(sel_dm) == 0L) { sel_dm <- jj; next }
      cors_j <- vapply(sel_dm, function(ss) {
        rp <- suppressWarnings(stats::cor(Z_full[anc_idx, jj], Z_full[anc_idx, ss], use = "complete.obs", method = "pearson"))
        rs <- suppressWarnings(stats::cor(Z_full[anc_idx, jj], Z_full[anc_idx, ss], use = "complete.obs", method = "spearman"))
        max(c(if (is.finite(rp)) abs(rp) else 0.0, if (is.finite(rs)) abs(rs) else 0.0))
      }, numeric(1))
      if (all(cors_j <= 0.70)) sel_dm <- c(sel_dm, jj)
    }
    keep_idx <- c(keep_idx, sort(sel_dm))
  }
  keep_idx <- sort(unique(keep_idx))
  cand_env <- cand_env[keep_idx]; E_raw <- E_raw[, keep_idx, drop = FALSE]; Z_std <- Z_full[, keep_idx, drop = FALSE]; dom_vec <- dom_full[keep_idx]

  dom_counts <- table(dom_vec)
  col_dom_wts <- vapply(dom_vec, function(dm) 1.0 / sqrt(max(1.0, as.numeric(dom_counts[[dm]]))), numeric(1))
  if (!is.null(rf_var_weights) && length(rf_var_weights) > 0L && !is.null(names(rf_var_weights))) {
    hit_w <- intersect(cand_env, names(rf_var_weights))
    if (length(hit_w) > 0L) {
      rw <- stats::setNames(rep(stats::median(as.numeric(rf_var_weights), na.rm = TRUE), length(cand_env)), cand_env)
      rw[hit_w] <- pmax(as.numeric(rf_var_weights[hit_w]), 1e-4); col_dom_wts <- col_dom_wts * sqrt(rw / max(mean(rw), 1e-6))
    }
  }
  Z_bal <- sweep(Z_std, 2L, col_dom_wts, "*"); Z_anc <- Z_bal[anc_idx, , drop = FALSE]; n_anc <- nrow(Z_anc)

  proj_pca <- function(M_anc, M_all = NULL) {
    if (ncol(M_anc) >= 2L && nrow(M_anc) >= 6L) {
      pe <- tryCatch(stats::prcomp(M_anc, center = TRUE, scale. = FALSE), error = function(e) NULL)
      if (!is.null(pe) && ncol(pe$rotation) >= 1L) {
        ev <- pe$sdev^2; kp <- which((cumsum(ev) / max(sum(ev), 1e-12)) >= 0.85)[1L]; if (is.na(kp)) kp <- ncol(pe$rotation)
        rot <- pe$rotation[, seq_len(max(min(2L, ncol(pe$rotation)), min(8L, kp, ncol(pe$rotation)))), drop = FALSE]
        return(list(S_anc = sweep(M_anc, 2L, pe$center, "-") %*% rot, S_all = if (!is.null(M_all)) sweep(M_all, 2L, pe$center, "-") %*% rot else NULL))
      }
    }
    list(S_anc = M_anc, S_all = M_all)
  }
  pp <- proj_pca(Z_anc, Z_bal); S_anc <- pp$S_anc; S_all <- pp$S_all
  tot_ss <- sum(sweep(S_anc, 2L, colMeans(S_anc), "-")^2)
  cfg_s05 <- if (exists("load_step05_macroecology_config", mode = "function")) load_step05_macroecology_config() else list(ecotype_min_silhouette = 0.40, ecotype_min_fold_jaccard = 0.60)
  eco_min_sil <- if (!is.null(cfg_s05$ecotype_min_silhouette) && is.finite(cfg_s05$ecotype_min_silhouette)) as.numeric(cfg_s05$ecotype_min_silhouette) else 0.40
  eco_min_jacc <- if (!is.null(cfg_s05$ecotype_min_fold_jaccard) && is.finite(cfg_s05$ecotype_min_fold_jaccard)) as.numeric(cfg_s05$ecotype_min_fold_jaccard) else 0.60
  min_cl_n <- if (is_intra && n_anc >= 40L) max(15L, as.integer(floor(0.10 * n_anc))) else if (n_anc >= 20L) 4L else 3L
  max_k_cand <- min(4L, as.integer(floor(n_anc / max(min_cl_n, 2L))))
  if (!is.finite(tot_ss) || tot_ss <= 1e-8 || max_k_cand < 2L) return(ret_k1())

  .old_seed_rag <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(42L)
  on.exit({ if (!is.null(.old_seed_rag)) assign(".Random.seed", .old_seed_rag, envir = .GlobalEnv) }, add = TRUE)

  d_anc_mat <- as.matrix(stats::dist(S_anc))
  calc_sil <- function(cl_i, dm = d_anc_mat) {
    uc <- sort(unique(cl_i)); if (length(uc) < 2L) return(0.0)
    mean(vapply(seq_along(cl_i), function(ii) {
      si <- which(cl_i == cl_i[ii] & seq_along(cl_i) != ii); ai <- if (length(si) > 0L) mean(dm[ii, si]) else 0.0
      bi <- min(vapply(setdiff(uc, cl_i[ii]), function(ck) mean(dm[ii, cl_i == ck]), numeric(1)))
      (bi - ai) / max(ai, bi, 1e-9)
    }, numeric(1)), na.rm = TRUE)
  }

  d_s <- ncol(S_anc); cov_s <- if (n_anc > 2L) stats::cov(S_anc) else diag(1.0, nrow = d_s, ncol = d_s)
  eig_s <- eigen(0.5 * (cov_s + t(cov_s)), symmetric = TRUE)
  L_s <- eig_s$vectors %*% diag(sqrt(pmax(eig_s$values, 0.0)), nrow = d_s, ncol = d_s)
  mk_draw <- function(Sb) list(S_b = Sb, dist_b = as.matrix(stats::dist(Sb)), tot_b = sum(sweep(Sb, 2L, colMeans(Sb), "-")^2))
  gauss_draws <- lapply(seq_len(20L), function(bi) mk_draw(matrix(stats::rnorm(n_anc * d_s), nrow = n_anc, ncol = d_s) %*% t(L_s)))
  perm_draws <- if (ncol(Z_anc) >= 2L) lapply(seq_len(20L), function(bi) { Zp <- Z_anc; for (jj in seq_len(ncol(Zp))) Zp[, jj] <- Zp[sample.int(n_anc), jj]; mk_draw(proj_pca(Zp)$S_anc) }) else list()
  eval_null_q95 <- function(dlist, k_tr) {
    if (length(dlist) == 0L) return(0.0)
    vals <- vapply(dlist, function(nd) {
      kmb <- tryCatch(stats::kmeans(nd$S_b, centers = k_tr, nstart = 10L, iter.max = 40L), error = function(e) NULL)
      if (is.null(kmb)) return(NA_real_)
      tb <- table(kmb$cluster); if (length(tb) < k_tr || any(tb < min_cl_n)) return(NA_real_)
      calc_sil(kmb$cluster, nd$dist_b) + 0.35 * (1.0 - kmb$tot.withinss / max(nd$tot_b, 1e-9))
    }, numeric(1))
    fv <- vals[is.finite(vals)]; if (length(fv) > 0L) as.numeric(stats::quantile(fv, 0.95, na.rm = TRUE, names = FALSE)) else 0.0
  }

  sub_jacc <- function(S_mat, ref_cl, k_tr, n_folds = 5L) {
    n_obs <- nrow(S_mat); if (n_obs < max(6L, k_tr + 2L)) return(0.0)
    folds_j <- integer(n_obs); for (ck in sort(unique(ref_cl))) { ix <- which(ref_cl == ck); folds_j[ix] <- ((seq_along(ix) - 1L) %% n_folds) + 1L }
    fv <- vapply(seq_len(n_folds), function(fi) {
      sub_i <- which(folds_j != fi); if (length(sub_i) < max(k_tr + 2L, 6L)) return(NA_real_)
      km_f <- tryCatch(stats::kmeans(S_mat[sub_i, , drop = FALSE], centers = k_tr, nstart = 15L, iter.max = 40L), error = function(e) NULL)
      if (is.null(km_f)) return(NA_real_)
      rc <- ref_cl[sub_i]; sc <- km_f$cluster; ur <- sort(unique(rc)); us <- sort(unique(sc))
      jm <- matrix(0.0, nrow = length(ur), ncol = length(us))
      for (ri in seq_along(ur)) { a <- rc == ur[ri]; for (si in seq_along(us)) { b <- sc == us[si]; uab <- sum(a | b); if (uab > 0L) jm[ri, si] <- sum(a & b) / uab } }
      used_s <- integer(0); jb <- numeric(length(ur))
      for (ri in order(apply(jm, 1L, max), decreasing = TRUE)) { av <- setdiff(seq_along(us), used_s); if (length(av) == 0L) next; bs <- av[which.max(jm[ri, av])[1L]]; used_s <- c(used_s, bs); jb[ri] <- jm[ri, bs] }
      mean(jb, na.rm = TRUE)
    }, numeric(1))
    fv <- fv[is.finite(fv)]; if (length(fv) > 0L) mean(fv) else 0.0
  }

  active_doms <- unique(dom_vec); cand_fits <- list()
  for (k_try in 2L:max_k_cand) {
    km <- tryCatch(stats::kmeans(S_anc, centers = k_try, nstart = 30L, iter.max = 60L), error = function(e) NULL)
    if (is.null(km) || length(table(km$cluster)) < k_try || any(table(km$cluster) < min_cl_n)) next
    sil_k <- calc_sil(km$cluster, d_anc_mat); var_red_k <- 1.0 - (km$tot.withinss / max(tot_ss, 1e-9))
    fold_jacc_k <- sub_jacc(S_anc, km$cluster, k_try)
    null_q95_combo <- max(eval_null_q95(gauss_draws, k_try), eval_null_q95(perm_draws, k_try))

    u_cls <- sort(unique(km$cluster))
    sigma_pool <- crossprod(S_anc - km$centers[km$cluster, , drop = FALSE]) / max(n_anc - length(u_cls), 1L) + diag(0.08, nrow = d_s, ncol = d_s)
    pair_min_mahal_d <- Inf; pair_min_best_d <- Inf; var_max_pair_d <- rep(0.0, length(cand_env))
    for (ia in seq_len(length(u_cls) - 1L)) for (ib in (ia + 1L):length(u_cls)) {
      ca <- u_cls[ia]; cb <- u_cls[ib]; diff_mu <- as.numeric(km$centers[ca, ] - km$centers[cb, ])
      mab <- tryCatch(sqrt(max(0.0, as.numeric(crossprod(diff_mu, solve(sigma_pool, diff_mu))))), error = function(e) 0.0)
      if (is.finite(mab) && mab < pair_min_mahal_d) pair_min_mahal_d <- mab
      best_ab_d <- 0.0
      for (j in seq_along(cand_env)) {
        ga <- Z_std[anc_idx[km$cluster == ca], j]; gb <- Z_std[anc_idx[km$cluster == cb], j]
        va <- if (length(ga) > 1L) stats::var(ga, na.rm = TRUE) else 0.0; vb <- if (length(gb) > 1L) stats::var(gb, na.rm = TRUE) else 0.0
        dab <- abs(mean(ga, na.rm = TRUE) - mean(gb, na.rm = TRUE)) / sqrt(0.5 * (ifelse(is.finite(va), va, 0.0) + ifelse(is.finite(vb), vb, 0.0)) + 0.04)
        if (is.finite(dab)) { if (dab > best_ab_d) best_ab_d <- dab; if (dab > var_max_pair_d[j]) var_max_pair_d[j] <- dab }
      }
      if (best_ab_d < pair_min_best_d) pair_min_best_d <- best_ab_d
    }

    kw_q <- stats::p.adjust(vapply(seq_along(cand_env), function(j) {
      vj <- Z_std[anc_idx, j]; if (stats::sd(vj) <= 1e-6) return(1.0)
      tryCatch(suppressWarnings(as.numeric(stats::kruskal.test(vj ~ factor(km$cluster))$p.value)), error = function(e) 1.0)
    }, numeric(1)), method = "BH")
    eta2_vec <- numeric(length(cand_env)); f_stats <- numeric(length(cand_env))
    for (j in seq_along(cand_env)) {
      vj <- Z_std[anc_idx, j]; if (stats::sd(vj) <= 1e-6) next
      eta2_vec[j] <- sum(vapply(u_cls, function(cc) { gc <- vj[km$cluster == cc]; length(gc) * ((mean(gc, na.rm = TRUE) - mean(vj, na.rm = TRUE))^2) }, numeric(1))) / max(sum((vj - mean(vj, na.rm = TRUE))^2), 1e-9)
      f_stats[j] <- tryCatch(suppressWarnings(as.numeric(summary(stats::aov(vj ~ factor(km$cluster)))[[1]][["F value"]][1L])), error = function(e) 0.0)
    }
    f_stats[!is.finite(f_stats)] <- 0.0
    min_p_cl <- min(as.numeric(table(km$cluster))) / max(as.numeric(n_anc), 1.0)
    eta2_thr_k <- max(0.10, min(0.20, 0.20 * (4.0 * min_p_cl * (1.0 - min_p_cl)) / 0.75))
    sig_mask_k <- is.finite(kw_q) & (kw_q < 0.005) & (eta2_vec >= eta2_thr_k) & (var_max_pair_d >= 1.60)
    n_sig_doms <- length(unique(dom_vec[sig_mask_k]))
    pass_eco <- is.finite(sil_k) && sil_k >= eco_min_sil && is.finite(fold_jacc_k) && fold_jacc_k >= eco_min_jacc &&
      is.finite(var_red_k) && var_red_k >= (0.30 + 0.05 * (k_try - 2L)) && is.finite(pair_min_mahal_d) && pair_min_mahal_d >= 2.05 &&
      is.finite(pair_min_best_d) && pair_min_best_d >= 1.60 && (sil_k + 0.35 * var_red_k) > null_q95_combo &&
      sum(sig_mask_k) >= (if (length(cand_env) >= 4L) 2L else 1L) && n_sig_doms >= min(2L, length(active_doms))
    score_k <- sil_k + 0.35 * var_red_k + 0.05 * min(n_sig_doms, 4L) + 0.02 * min(pair_min_mahal_d, 6.0) - 0.02 * (k_try - 2L)
    cand_fits[[as.character(k_try)]] <- list(k = k_try, km = km, sil = sil_k, var_red = var_red_k, f_stats = stats::setNames(f_stats, cand_env), n_sig_doms = n_sig_doms, pass = pass_eco, score = score_k)
  }

  valid_fits <- Filter(function(x) isTRUE(x$pass), cand_fits)
  if (length(valid_fits) == 0L) { if (allow_k1 || length(cand_fits) == 0L) return(ret_k1(msg = TRUE)); valid_fits <- cand_fits }
  best_fit <- valid_fits[[which.max(vapply(valid_fits, function(x) x$score, numeric(1)))]]
  k_opt <- best_fit$k; centers_pc <- best_fit$km$centers
  d2_all <- matrix(0.0, nrow = n_rows, ncol = k_opt)
  for (kc in seq_len(k_opt)) d2_all[, kc] <- rowSums(sweep(S_all, 2L, centers_pc[kc, ], "-")^2)
  raw_cl <- max.col(-d2_all, ties.method = "first")
  top_var <- names(which.max(best_fit$f_stats))[1L]; if (is.na(top_var) || !nzchar(top_var)) top_var <- cand_env[1L]
  clean_tok <- function(vn) gsub("[^A-Za-z0-9_]", "_", vn)
  ord_cl <- order(tapply(E_raw[anc_idx, top_var], raw_cl[anc_idx], mean, na.rm = TRUE), decreasing = FALSE); lbl_map <- character(k_opt)
  if (k_opt == 2L) {
    lbl_map[ord_cl[1L]] <- sprintf("Ecotype_1_%s_Low", clean_tok(top_var)); lbl_map[ord_cl[2L]] <- sprintf("Ecotype_2_%s_High", clean_tok(top_var))
  } else {
    used_s <- character(0)
    for (rk in seq_len(k_opt)) {
      oc <- ord_cl[rk]; zc <- colMeans(Z_std[anc_idx[raw_cl[anc_idx] == oc], , drop = FALSE], na.rm = TRUE); zc[!is.finite(zc)] <- 0.0; sig <- NULL
      for (cj in order(abs(zc), decreasing = TRUE)) { st <- paste0(clean_tok(cand_env[cj]), "_", if (zc[cj] >= 0.0) "High" else "Low"); if (!(st %in% used_s)) { sig <- st; used_s <- c(used_s, st); break } }
      if (is.null(sig)) sig <- paste0(clean_tok(top_var), "_", if (rk == 1L) "Low" else if (rk == k_opt) "High" else paste0("Mid", rk - 1L))
      lbl_map[oc] <- sprintf("Ecotype_%d_%s", rk, sig)
    }
  }
  resolved_lbls <- lbl_map[raw_cl]; tbl_res <- table(resolved_lbls)
  message(sprintf("[Independent Step 05 ENM] Discovered K=%d ecological regimes for '%s' (sil=%.3f, var_red=%.1f%%, sig_domains=%d): %s",
                  k_opt, base_taxon, best_fit$sil, 100.0 * best_fit$var_red, best_fit$n_sig_doms, paste(names(tbl_res), as.integer(tbl_res), sep = "=", collapse = ", ")))
  df$analysis_group <- as.character(resolved_lbls); df$ecotype_regime <- as.character(resolved_lbls)
  attr(df, "is_post_selection_ecotype") <- isTRUE(k_opt >= 2L); df
}

resolve_active_anchor_groups <- function(df, intraspecific_mode = FALSE) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(character(0))
  n_r <- nrow(df)
  is_intra <- isTRUE(intraspecific_mode) || (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  g_raw <- if (is_intra && "analysis_group" %in% names(df)) as.character(df$analysis_group)
    else if ("species_standardized" %in% names(df)) as.character(df$species_standardized)
    else if ("analysis_group" %in% names(df)) as.character(df$analysis_group) else rep(NA_character_, n_r)
  fl <- .step05_anchor_flags(df, min_edaph_n = 4L)
  miss_gmm <- if (is_intra && "is_missing_gmm_cluster" %in% names(df)) .to_bool_s05(df$is_missing_gmm_cluster) else rep(FALSE, n_r)
  tier_ok <- grepl("Tier_1|Gold|Tier_2|Silver", fl$det_t, ignore.case = TRUE) | fl$is_prom | fl$is_t
  valid_g <- !is.na(g_raw) & nzchar(g_raw) & (g_raw != "Unknown") & ((!fl$cl_noise) | fl$is_t) & !fl$quar & !fl$open_h
  ext_m <- .is_step05_extralimital_non_type(df, group_vec = g_raw, clean_base_mask = (tier_ok & fl$elig & valid_g), is_type_vec = fl$is_t)
  base_anc <- tier_ok & fl$elig & valid_g & !ext_m
  is_anc <- base_anc & !miss_gmm
  if (sum(is_anc) < 4L) is_anc <- base_anc
  if (!any(tier_ok) && sum(is_anc) == 0L) { is_anc <- fl$elig & valid_g & !miss_gmm; if (sum(is_anc) == 0L) is_anc <- valid_g }
  ifelse(is_anc, g_raw, NA_character_)
}

compute_isolation_by_environment <- function(df, edaph_vars = character(0), clim_vars = character(0), mem_mat = NULL, n_perm = 99L) {
  n <- nrow(df)
  empty_ibe <- data.frame(Modality = character(0), N_Specimens = integer(0), Pure_Edaphic_AdjR2 = numeric(0), Pure_Climate_AdjR2 = numeric(0),
                          Pure_Space_MEM_AdjR2 = numeric(0), Total_Explained_AdjR2 = numeric(0), Partial_RDA_F = numeric(0),
                          Partial_RDA_p = numeric(0), IBE_Inference = character(0), stringsAsFactors = FALSE)
  if (is.null(df) || n < 8L) return(empty_ibe)
  clean_blk <- function(cols) {
    if (length(cols) == 0L) return(NULL)
    M <- matrix(0.0, nrow = n, ncol = length(cols), dimnames = list(NULL, cols)); keep <- logical(length(cols))
    for (j in seq_along(cols)) {
      v <- suppressWarnings(as.numeric(df[[cols[j]]])); fin <- is.finite(v)
      if (sum(fin) >= 5L && stats::sd(v[fin]) > 1e-8) { v[!fin] <- stats::median(v[fin], na.rm = TRUE); M[, j] <- as.numeric(scale(v)); keep[j] <- TRUE }
    }
    if (!any(keep)) NULL else M[, keep, drop = FALSE]
  }
  X_ed <- clean_blk(intersect(as.character(edaph_vars), names(df)))
  X_cl <- clean_blk(intersect(setdiff(as.character(clim_vars), "delta_doy"), names(df)))
  if (is.null(X_ed) && is.null(X_cl)) return(empty_ibe)
  if (is.null(mem_mat) || nrow(as.matrix(mem_mat)) != n) {
    if (all(c("longitude", "latitude") %in% names(df))) {
      lo <- suppressWarnings(as.numeric(df$longitude)); la <- suppressWarnings(as.numeric(df$latitude))
      lo[!is.finite(lo)] <- stats::median(lo, na.rm = TRUE); la[!is.finite(la)] <- stats::median(la, na.rm = TRUE)
      mem_mat <- scale(cbind(MEM_Lon = lo, MEM_Lat = la))
    } else mem_mat <- matrix(0.0, nrow = n, ncol = 1L, dimnames = list(NULL, "MEM1"))
  } else { mem_mat <- scale(as.matrix(mem_mat)); mem_mat[!is.finite(mem_mat)] <- 0.0 }

  eval_blk <- function(Y_raw, mod_name) {
    if (is.null(Y_raw) || ncol(Y_raw) == 0L) return(NULL)
    Y <- scale(as.matrix(Y_raw)); Y[!is.finite(Y)] <- 0.0; X_env <- cbind(X_ed, X_cl)
    adj_r2 <- function(Xm) {
      if (is.null(Xm) || ncol(Xm) == 0L) return(0.0)
      qr_x <- qr(cbind(1.0, Xm)); p_k <- qr_x$rank - 1L; if (p_k < 1L || (n - qr_x$rank) < 2L) return(0.0)
      tss <- sum(sweep(Y, 2L, colMeans(Y), "-")^2); if (tss <= 1e-10) return(0.0)
      1.0 - (1.0 - max(0.0, 1.0 - sum(qr.resid(qr_x, Y)^2) / tss)) * ((n - 1.0) / max(1.0, n - p_k - 1.0))
    }
    r2_all <- adj_r2(cbind(X_ed, X_cl, mem_mat))
    pure_ed <- max(0.0, r2_all - adj_r2(cbind(X_cl, mem_mat))); pure_cl <- max(0.0, r2_all - adj_r2(cbind(X_ed, mem_mat)))
    pure_sp <- max(0.0, r2_all - adj_r2(cbind(X_ed, X_cl)))
    nw_ibe <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L) else 1L
    f_obs <- NA_real_; p_perm <- NA_real_
    if (requireNamespace("vegan", quietly = TRUE) && ncol(X_env) >= 1L && (n - ncol(X_env) - ncol(mem_mat)) >= 4L) {
      rr <- tryCatch({ an <- suppressMessages(suppressWarnings(vegan::anova.cca(vegan::rda(Y ~ X_env + Condition(mem_mat)), permutations = as.integer(max(19L, n_perm)), parallel = nw_ibe))); list(F_val = as.numeric(an$F[1L]), p_val = as.numeric(an$`Pr(>F)`[1L])) }, error = function(e) NULL)
      if (!is.null(rr)) { f_obs <- rr$F_val; p_perm <- rr$p_val }
    }
    if (!is.finite(p_perm)) {
      q0 <- qr(cbind(1.0, mem_mat)); q1 <- qr(cbind(1.0, mem_mat, X_env))
      rss0 <- sum(qr.resid(q0, Y)^2); rss1 <- sum(qr.resid(q1, Y)^2)
      df1 <- max(1L, q1$rank - q0$rank); df2 <- max(1L, n - q1$rank); f_obs <- ((rss0 - rss1) / df1) / max(rss1 / df2, 1e-8)
      Y_r0 <- qr.resid(q0, Y); perms_i <- replicate(n_perm, sample.int(n), simplify = FALSE)
      ev_p <- function(b) { Yp <- Y_r0[perms_i[[b]], , drop = FALSE]; ((sum(qr.resid(q0, Yp)^2) - sum(qr.resid(q1, Yp)^2)) / df1) / max(sum(qr.resid(q1, Yp)^2) / df2, 1e-8) }
      f_null <- if (n_perm >= 99L && nw_ibe > 1L && exists("packera_mclapply", mode = "function")) unlist(packera_mclapply(seq_len(n_perm), ev_p, mc.cores = min(nw_ibe, n_perm)), use.names = FALSE) else vapply(seq_len(n_perm), ev_p, numeric(1))
      p_perm <- (sum(f_null >= f_obs) + 1.0) / (n_perm + 1.0)
    }
    inf_lbl <- if (is.finite(p_perm) && p_perm < 0.01 && (pure_ed + pure_cl) >= 0.05) "Significant Isolation by Environment (Adaptive Divergence Candidate)"
      else if (is.finite(p_perm) && p_perm < 0.05) "Weak Environmental Cline (Plasticity / Gradual Ecotone)" else "Null IBE (Panmictic / Facultative Niche Breadth)"
    data.frame(Modality = mod_name, N_Specimens = as.integer(n), Pure_Edaphic_AdjR2 = round(pure_ed, 4), Pure_Climate_AdjR2 = round(pure_cl, 4),
               Pure_Space_MEM_AdjR2 = round(pure_sp, 4), Total_Explained_AdjR2 = round(max(0.0, r2_all), 4),
               Partial_RDA_F = round(f_obs, 4), Partial_RDA_p = round(p_perm, 4), IBE_Inference = inf_lbl, stringsAsFactors = FALSE)
  }

  rows_out <- list()
  Ym <- clean_blk(unique(c(grep("^(PC|ECT_PC)[1-4]$", names(df), value = TRUE), intersect(c("aspect_ratio", "petiole_to_blade_ratio", "solidity", "dissection_index"), names(df)))))
  if (!is.null(Ym)) rows_out[[length(rows_out) + 1L]] <- eval_blk(Ym, "Foliar_Morphology_IBE")
  Yp <- clean_blk(intersect("delta_doy", names(df)))
  if (!is.null(Yp)) rows_out[[length(rows_out) + 1L]] <- eval_blk(Yp, "Anthesis_Phenology_IBE")
  if (length(rows_out) == 0L) empty_ibe else do.call(rbind, rows_out)
}

compute_geographic_disjunction_audit <- function(df, min_gap_km = 250.0, min_gap_ratio = 2.5) {
  n <- if (is.data.frame(df)) nrow(df) else 0L
  empty_sum <- data.frame(Taxon = character(0), N_Georef_Clean = integer(0), Max_MST_Gap_km = numeric(0), Q90_MST_Edge_km = numeric(0),
                          Gap_Ratio = numeric(0), N_Primary_Sector = integer(0), N_Disjunct_Sector = integer(0), Disjunction_Sector = character(0),
                          Has_Major_Disjunction = logical(0), stringsAsFactors = FALSE)
  if (n == 0L) return(list(sector_vec = character(0), gap_km_vec = numeric(0), summary_df = empty_sum))
  la <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else if ("latitude" %in% names(df)) suppressWarnings(as.numeric(df$latitude)) else rep(NA_real_, n)
  lo <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else if ("longitude" %in% names(df)) suppressWarnings(as.numeric(df$longitude)) else rep(NA_real_, n)
  has_xy <- is.finite(la) & is.finite(lo) & !(abs(la) < 1e-4 & abs(lo) < 1e-4)
  fl <- .step05_anchor_flags(df, min_edaph_n = 4L)
  sp_v <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unknown", n)
  st_v <- if ("stateProvince" %in% names(df)) toupper(trimws(ifelse(is.na(df$stateProvince), "", as.character(df$stateProvince)))) else rep("", n)
  sec_vec <- rep(NA_character_, n); gap_vec <- rep(NA_real_, n); sum_rows <- list()
  for (tx in sort(unique(sp_v[!is.na(sp_v) & nzchar(sp_v) & sp_v != "Unknown"]))) {
    sp_idx <- which(sp_v == tx); g_idx <- sp_idx[has_xy[sp_idx] & !fl$quar[sp_idx] & !fl$open_h[sp_idx] & !fl$cl_noise[sp_idx]]
    ng <- length(g_idx); if (ng < 12L) next
    D_km <- compute_haversine_dist_km_matrix(cbind(longitude = lo[g_idx], latitude = la[g_idx])); diag(D_km) <- Inf
    in_t <- logical(ng); in_t[1L] <- TRUE; mu <- integer(ng - 1L); mv <- integer(ng - 1L); mw <- numeric(ng - 1L)
    bd <- D_km[1L, ]; bp <- rep(1L, ng)
    for (ei in seq_len(ng - 1L)) {
      cd <- which(!in_t); pk <- cd[which.min(bd[cd])[1L]]; mu[ei] <- bp[pk]; mv[ei] <- pk; mw[ei] <- bd[pk]; in_t[pk] <- TRUE
      rm_i <- which(!in_t); if (length(rm_i) > 0L) { im <- rm_i[D_km[pk, rm_i] < bd[rm_i]]; if (length(im) > 0L) { bd[im] <- D_km[pk, im]; bp[im] <- pk } }
    }
    min_cl <- max(4L, as.integer(ceiling(0.04 * ng))); best_e <- NA_integer_; comp_a <- integer(0)
    for (ei in order(mw, decreasing = TRUE)) {
      vis <- logical(ng); q <- mu[ei]; vis[q] <- TRUE; hd <- 1L
      while (hd <= length(q)) {
        cur <- q[hd]; hd <- hd + 1L
        nb <- c(mv[mu == cur & seq_along(mu) != ei], mu[mv == cur & seq_along(mv) != ei]); nb <- nb[!vis[nb]]
        if (length(nb) > 0L) { vis[nb] <- TRUE; q <- c(q, nb) }
      }
      ca <- which(vis); if (min(length(ca), ng - length(ca)) >= min_cl) { best_e <- ei; comp_a <- ca; break }
    }
    if (is.na(best_e)) next
    comp_b <- setdiff(seq_len(ng), comp_a)
    if (length(comp_a) <= length(comp_b)) { min_c <- comp_a; maj_c <- comp_b } else { min_c <- comp_b; maj_c <- comp_a }
    q90_w <- as.numeric(stats::quantile(mw, probs = 0.90, names = FALSE)); g_km <- mw[best_e]; g_rat <- g_km / max(q90_w, 10.0)
    is_disj <- is.finite(g_km) && (g_km >= min_gap_km) && (g_rat >= min_gap_ratio)
    sec_lbl <- if (is_disj) {
      if (mean(lo[g_idx[min_c]]) < -90.5 && mean(lo[g_idx[maj_c]]) > -89.0) "Trans_Mississippi_Disjunct" else "Allopatric_Disjunct_Sector"
    } else "Continuous_Range"
    gap_vec[sp_idx] <- round(g_km, 2)
    if (is_disj) {
      sec_vec[g_idx[min_c]] <- sec_lbl; min_st <- unique(st_v[g_idx[min_c]][nzchar(st_v[g_idx[min_c]])]); maj_st <- unique(st_v[g_idx[maj_c]][nzchar(st_v[g_idx[maj_c]])])
      ug_sp <- sp_idx[!has_xy[sp_idx] & nzchar(st_v[sp_idx]) & (st_v[sp_idx] %in% setdiff(min_st, maj_st))]
      if (length(ug_sp) > 0L) sec_vec[ug_sp] <- sec_lbl
    }
    sum_rows[[length(sum_rows) + 1L]] <- data.frame(
      Taxon = tx, N_Georef_Clean = ng, Max_MST_Gap_km = round(g_km, 2), Q90_MST_Edge_km = round(q90_w, 2), Gap_Ratio = round(g_rat, 3),
      N_Primary_Sector = length(maj_c), N_Disjunct_Sector = length(min_c), Disjunction_Sector = sec_lbl, Has_Major_Disjunction = is_disj, stringsAsFactors = FALSE)
  }
  list(sector_vec = sec_vec, gap_km_vec = gap_vec, summary_df = if (length(sum_rows) > 0L) do.call(rbind, sum_rows) else empty_sum)
}
