if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.step07_safe_bool <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & (x != 0))
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}
.num_col <- function(df, col_nm, n = nrow(df), default = NA_real_) if (col_nm %in% names(df)) suppressWarnings(as.numeric(df[[col_nm]])) else rep(default, n)
.chr_col <- function(df, col_nm, n = nrow(df), default = NA_character_) if (col_nm %in% names(df)) as.character(df[[col_nm]]) else rep(default, n)
.lgl_col <- function(df, col_nm, n = nrow(df)) if (col_nm %in% names(df)) (!is.na(df[[col_nm]]) & as.logical(df[[col_nm]])) else rep(FALSE, n)

.coalesce_step07_ddoy_resid <- function(df, n = nrow(df)) {
  out <- rep(NA_real_, n)
  for (dc in c("delta_doy_resid", "hopkins_adjusted_doy_resid", "delta_doy")) if (dc %in% names(df)) { dv <- suppressWarnings(as.numeric(df[[dc]])); fm <- !is.finite(out) & is.finite(dv); if (any(fm)) out[fm] <- dv[fm] }
  out
}

.has_step07_geno_row_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  if (exists("has_active_genomic_row_mask", mode = "function")) return(has_active_genomic_row_mask(df))
  nms <- names(df); if (!("geno_predicted_taxon" %in% nms)) return(rep(FALSE, n))
  pred <- as.character(df$geno_predicted_taxon); pred_ok <- !is.na(pred) & nzchar(trimws(pred)) & (pred != "Unknown")
  tier_col <- if ("geno_linkage_tier" %in% nms) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% nms) as.character(df$linkage_tier) else rep("G1_Direct_Voucher", n)
  tier_ok <- tier_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G2_Conditional_Deme_Hull")
  qc_ok <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(TRUE, n)
  prob_ok <- if ("geno_posterior_prob" %in% nms) is.finite(suppressWarnings(as.numeric(df$geno_posterior_prob))) else rep(TRUE, n)
  pred_ok & tier_ok & qc_ok & prob_ok
}

.step07_type_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  is_t <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df)$is_type_specimen else rep(FALSE, n); is_t[is.na(is_t)] <- FALSE
  if ("is_type_specimen" %in% names(df)) { tv <- suppressWarnings(as.logical(df$is_type_specimen)); tv[is.na(tv)] <- FALSE; is_t <- is_t | tv }
  if ("determiner_tier" %in% names(df)) is_t <- is_t | (trimws(as.character(df$determiner_tier)) == "Tier_0_Type_Anchor")
  is_t
}

.is_step07_centroid_coord <- function(df) { n <- nrow(df); if (n == 0L || !("is_centroid_or_inst_coord" %in% names(df))) return(rep(FALSE, n)); .step07_safe_bool(df$is_centroid_or_inst_coord) & !.step07_type_mask(df) }

.circ_median_doy_step07 <- function(dv) {
  v <- dv[is.finite(dv) & dv >= 1.0 & dv <= 366.0]; if (length(v) == 0L) return(NA_real_)
  th <- 2.0 * pi * (v / 365.25); mu <- atan2(mean(sin(th)), mean(cos(th))); if (mu < 0.0) mu <- mu + 2.0 * pi
  off <- atan2(sin(th - mu), cos(th - mu)) * (365.25 / (2.0 * pi))
  ((mu * (365.25 / (2.0 * pi)) + stats::median(off) - 1.0) %% 365.25) + 1.0
}

resolve_active_anthesis_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  ver_vec <- if ("is_verified_anthesis" %in% names(df)) suppressWarnings(as.logical(df$is_verified_anthesis)) else rep(NA, n)
  na_cnt <- if ("n_anthesis" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis)) else if ("n_anthesis_heads" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis_heads)) else rep(NA_real_, n)
  state_cols <- intersect(c("voucher_phenological_state", "phenological_state", "phenophase", "head_state"), names(df))
  is_anthesis <- if (length(state_cols) > 0L) {
    st_raw <- tolower(trimws(as.character(df[[state_cols[1L]]]))); unk_mask <- is.na(st_raw) | !nzchar(st_raw) | st_raw == "unknown"
    ia <- (!unk_mask) & (st_raw %in% c("anthesis", "flowering", "active_anthesis", "open_flower") | grepl("anthesis|flower", st_raw))
    ia[!is.na(st_raw) & (st_raw %in% c("bud", "fruit", "sterile", "vegetative", "fruiting", "senescent") | grepl("pre[_-]?anthesis|post[_-]?anthesis|non[_-]?anthesis", st_raw))] <- FALSE
    if (any(unk_mask)) {
      has_ver <- unk_mask & !is.na(ver_vec); if (any(has_ver)) ia[has_ver] <- ver_vec[has_ver]
      rem_unk <- unk_mask & is.na(ver_vec) & is.finite(na_cnt); if (any(rem_unk)) ia[rem_unk] <- na_cnt[rem_unk] >= 1.0
    }
    ia
  } else if ("is_verified_anthesis" %in% names(df)) {
    ia <- ver_vec; if (any(is.na(ia)) && any(is.finite(na_cnt))) ia[is.na(ia) & is.finite(na_cnt)] <- na_cnt[is.na(ia) & is.finite(na_cnt)] >= 1.0
    if (all(is.na(ia)) && ("is_active_anthesis" %in% names(df))) { ia_fb <- suppressWarnings(as.logical(df$is_active_anthesis)); ia_fb[is.na(ia_fb)] <- FALSE; ia_fb } else { ia[is.na(ia)] <- FALSE; ia }
  } else if ("is_active_anthesis" %in% names(df)) { ia <- suppressWarnings(as.logical(df$is_active_anthesis)); ia[is.na(ia)] <- FALSE; ia
  } else if ("is_anthesis" %in% names(df)) { ia <- as.logical(df$is_anthesis); ia[is.na(ia)] <- FALSE; ia
  } else if ("n_anthesis" %in% names(df)) { !is.na(na_cnt) & is.finite(na_cnt) & na_cnt >= 1.0
  } else if ("n_anthesis_heads" %in% names(df)) { nh <- suppressWarnings(as.numeric(df$n_anthesis_heads)); !is.na(nh) & is.finite(nh) & nh >= 1.0
  } else rep(TRUE, n)
  if ("anthesis_filter_status" %in% names(df)) { afs <- tolower(trimws(as.character(df$anthesis_filter_status))); is_anthesis[!is.na(afs) & (afs == "late_fruiting_senescent_excluded")] <- FALSE }
  if ("anthesis_dominance_weight" %in% names(df)) { adw <- suppressWarnings(as.numeric(df$anthesis_dominance_weight)); is_anthesis[is.finite(adw) & (adw < 0.20)] <- FALSE }
  defs <- .load_step07_pheno_config_defaults(); doy_v <- .num_col(df, "doy", n)
  if (all(!is.finite(doy_v)) && ("hopkins_adjusted_doy" %in% names(df))) doy_v <- .num_col(df, "hopkins_adjusted_doy", n)
  has_d <- is.finite(doy_v) & (doy_v >= 1.0) & (doy_v <= 366.0)
  if (any(has_d)) {
    is_anthesis[has_d & (doy_v < defs$min_flowering_doy)] <- FALSE
    coh_med <- .circ_median_doy_step07(doy_v[is_anthesis & has_d])
    sp_col <- if ("species_standardized" %in% names(df)) trimws(as.character(df$species_standardized)) else rep("", n)
    med_by_row <- rep(coh_med, n); ok_sp <- !is.na(sp_col) & nzchar(sp_col) & sp_col != "Unknown"
    for (tx in unique(sp_col[ok_sp])) { m_tx <- ok_sp & (sp_col == tx) & is_anthesis & has_d; if (sum(m_tx) >= 5L && is.finite(tx_med <- .circ_median_doy_step07(doy_v[m_tx]))) med_by_row[ok_sp & (sp_col == tx)] <- tx_med }
    is_anthesis[is.finite(med_by_row) & (med_by_row < 170.0) & has_d & ((doy_v < defs$min_flowering_doy) | (doy_v > defs$max_flowering_doy) | (doy_v > 240.0))] <- FALSE
  }
  is_anthesis
}

.tier0_geno_ctx <- function(df) {
  n <- nrow(df); nms <- names(df)
  tier_g <- if ("geno_linkage_tier" %in% nms) trimws(as.character(df$geno_linkage_tier)) else if ("linkage_tier" %in% nms) trimws(as.character(df$linkage_tier)) else rep("", n)
  qc_g <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(FALSE, n)
  adm_g <- if ("is_geno_admixed" %in% nms) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
  qry_g <- if ("is_putative_novel_or_hybrid_query" %in% nms) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  pred_g <- trimws(.chr_col(df, "geno_predicted_taxon", n, "")); pred_g[is.na(pred_g)] <- ""
  sp_g <- trimws(.chr_col(df, "species_standardized", n, "")); sp_g[is.na(sp_g)] <- ""
  prob_g <- .num_col(df, "geno_posterior_prob", n); marg_g <- .num_col(df, "geno_margin", n)
  g1_ok <- (tier_g %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc_g & nzchar(pred_g) & (pred_g != "Unknown") & is.finite(prob_g) & (prob_g >= 0.85)
  list(g1_ok = g1_ok, adm_g = adm_g, qry_g = qry_g, pred_g = pred_g, sp_g = sp_g, marg_g = marg_g)
}

resolve_step07_ecological_authority_weights <- function(df, allow_bronze_prior = TRUE) {
  n <- nrow(df); if (n == 0L) return(numeric(0))
  w_out <- if (exists("resolve_ecological_authority_weights", mode = "function")) resolve_ecological_authority_weights(df, allow_bronze_prior = allow_bronze_prior) else rep(1.0, n)
  if (isTRUE(allow_bronze_prior) && length(w_out) == n) w_out <- round(pmax(0.01, w_out), 4)
  if (any(.has_step07_geno_row_mask(df))) {
    is_type_m <- .step07_type_mask(df); gx <- .tier0_geno_ctx(df)
    w_out[!is_type_m & gx$g1_ok & !gx$adm_g & !gx$qry_g & (is.na(gx$marg_g) | gx$marg_g >= 0.40) & (gx$pred_g == gx$sp_g)] <- 4.0
    w_out[!is_type_m & gx$g1_ok & (gx$pred_g != gx$sp_g)] <- 0.0
  }
  w_out
}

.load_step07_pheno_config_defaults <- function(cfg_path = "config/config.yaml") {
  defs <- list(tau_ecotype = 14.0, tau_cryptic = 18.0, min_flowering_doy = 30.0, max_flowering_doy = 220.0)
  if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    if (is.list(cfg) && is.list(cfg$synthesis)) {
      if (!is.null(cfg$synthesis$ecotype_ddoy_cutoff_days)) defs$tau_ecotype <- as.numeric(cfg$synthesis$ecotype_ddoy_cutoff_days)
      if (!is.null(cfg$synthesis$cryptic_ddoy_cutoff_days)) defs$tau_cryptic <- as.numeric(cfg$synthesis$cryptic_ddoy_cutoff_days)
    }
    ph_c <- if (is.list(cfg) && is.list(cfg$macroecology) && is.list(cfg$macroecology$phenology)) cfg$macroecology$phenology else if (is.list(cfg) && is.list(cfg$phenology)) cfg$phenology else NULL
    if (is.list(ph_c)) { mn <- suppressWarnings(as.numeric(ph_c$min_flowering_doy)); mx <- suppressWarnings(as.numeric(ph_c$max_flowering_doy)); if (length(mn) == 1L && is.finite(mn)) defs$min_flowering_doy <- mn; if (length(mx) == 1L && is.finite(mx)) defs$max_flowering_doy <- mx }
  }
  defs
}

.circular_sd_days_from_residuals <- function(ddoy_days_vec) {
  v <- ddoy_days_vec[is.finite(ddoy_days_vec)]; if (length(v) < 2L) return(NA_real_)
  phi <- 2.0 * pi * (v / 365.25); r_clamped <- max(1e-6, min(0.999999, sqrt(mean(sin(phi))^2 + mean(cos(phi))^2)))
  (365.25 / (2.0 * pi)) * sqrt(-2.0 * log(r_clamped))
}

.is_step07_extralimital_non_type <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0)); nms <- names(df); is_type_m <- .step07_type_mask(df)
  u_stat <- if ("usda_county_status" %in% nms) tolower(trimws(ifelse(is.na(df$usda_county_status), "", as.character(df$usda_county_status)))) else rep("", n)
  exs_c <- if ("exsiccatae_taxonomic_conflict" %in% nms) { rc <- df$exsiccatae_taxonomic_conflict; if (is.logical(rc)) (!is.na(rc) & rc) else (!is.na(suppressWarnings(as.integer(rc))) & suppressWarnings(as.integer(rc)) != 0L) } else rep(FALSE, n)
  pas_aud <- if ("anchor_tier_source" %in% nms) (trimws(ifelse(is.na(df$anchor_tier_source), "", as.character(df$anchor_tier_source))) == "Passive_Extralimital_Audit") else rep(FALSE, n)
  flag <- !is_type_m & ((u_stat %in% c("out_of_range", "natureserve_exotic_or_waif")) | exs_c | pas_aud)
  if (any(flag) && ("species_standardized" %in% nms)) {
    sp_v <- trimws(ifelse(is.na(df$species_standardized), "", as.character(df$species_standardized))); noise_m <- .lgl_col(df, "is_confident_label_noise", n)
    for (tx in unique(sp_v[flag & nzchar(sp_v)])) if (sum(!flag[which(sp_v == tx & !noise_m)]) < 2L) flag[sp_v == tx] <- FALSE
  }
  flag
}

.resolve_step07_anchor_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  tier_col <- if ("determiner_tier" %in% names(df)) {
    if (exists(".standardize_determiner_tier_vec", mode = "function")) .standardize_determiner_tier_vec(df)
    else ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier))
  } else rep("Tier_1_Gold", n)
  ats_v <- trimws(.chr_col(df, "anchor_tier_source", n, "")); ats_v[is.na(ats_v)] <- ""
  css_v <- trimws(.chr_col(df, "concept_synonymy_status", n, "")); css_v[is.na(css_v)] <- ""
  quarantined_concept <- (ats_v %in% c("Passive_Heterotypic_Type_Audit", "Passive_Type_Uncertain_Coord", "Passive_Open_or_Hybrid_Concept")) |
    (css_v == "Heterotypic_Synonym") | (if ("is_open_or_hybrid_concept" %in% names(df)) .step07_safe_bool(df$is_open_or_hybrid_concept) else rep(FALSE, n))
  is_type_m <- .step07_type_mask(df); is_clean_type_m <- is_type_m & !quarantined_concept
  promo_mask <- (.lgl_col(df, "is_promoted_bronze", n) | is_clean_type_m) & !quarantined_concept
  noise_mask <- .lgl_col(df, "is_confident_label_noise", n); passive_mask <- .lgl_col(df, "is_passive", n)
  if ("authority_weight" %in% names(df)) { aw_v <- .num_col(df, "authority_weight", n); if (any(is.finite(aw_v) & aw_v > 0.0)) passive_mask <- passive_mask | (is.finite(aw_v) & aw_v <= 0.0) }
  anchor_mask <- ((((tier_col %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze_Promoted")) | promo_mask) & !noise_mask & !passive_mask & !.is_step07_extralimital_non_type(df)) | (is_clean_type_m & !noise_mask)) & !quarantined_concept
  if (any(.has_step07_geno_row_mask(df))) {
    gx <- .tier0_geno_ctx(df); g1_base <- gx$g1_ok & !gx$adm_g
    t0_contra <- !is_clean_type_m & g1_base & (gx$pred_g != gx$sp_g); t0_conc <- g1_base & !gx$qry_g & (gx$pred_g == gx$sp_g) & !noise_mask & !quarantined_concept
    anchor_mask <- (anchor_mask | t0_conc) & !t0_contra & !quarantined_concept
    if (sum(anchor_mask, na.rm = TRUE) < 4L) {
      clean_fb <- !noise_mask & !passive_mask & !t0_contra & !quarantined_concept
      anchor_mask <- if (sum(clean_fb, na.rm = TRUE) >= 1L) clean_fb else if (sum(!noise_mask & !t0_contra & !quarantined_concept, na.rm = TRUE) >= 1L) (!noise_mask & !t0_contra & !quarantined_concept) else if (sum(!noise_mask & !t0_contra, na.rm = TRUE) >= 1L) (!noise_mask & !t0_contra) else if (sum(!noise_mask, na.rm = TRUE) >= 1L) !noise_mask else rep(TRUE, n)
    }
  } else if (sum(anchor_mask, na.rm = TRUE) < 4L) {
    clean_fb <- !noise_mask & !passive_mask & !quarantined_concept
    anchor_mask <- if (sum(clean_fb, na.rm = TRUE) >= 1L) clean_fb else if (sum(!noise_mask & !quarantined_concept, na.rm = TRUE) >= 1L) (!noise_mask & !quarantined_concept) else if (sum(!noise_mask, na.rm = TRUE) >= 1L) !noise_mask else rep(TRUE, n)
  }
  anchor_mask
}

.compute_taxon_pheno_thresholds <- function(df) {
  n <- nrow(df); defs <- .load_step07_pheno_config_defaults()
  if (n == 0L) return(list(tau_ecotype = numeric(0), tau_cryptic = numeric(0), sigma_by_taxon = numeric(0)))
  sp_vec <- .chr_col(df, "species_standardized", n, "Unknown")
  is_anth <- if ("is_active_anthesis" %in% names(df)) as.logical(df$is_active_anthesis) else resolve_active_anthesis_mask(df); is_anth[is.na(is_anth)] <- TRUE
  ddoy_vec <- .coalesce_step07_ddoy_resid(df, n); valid_ph <- is_anth & is.finite(ddoy_vec)
  mo_m <- .lgl_col(df, "morphometric_outlier_flag", n); noise_m <- .lgl_col(df, "is_confident_label_noise", n); clean_anch <- .resolve_step07_anchor_mask(df) & !mo_m & !noise_m
  coh_idx <- which(clean_anch & valid_ph); if (length(coh_idx) < 5L) coh_idx <- which(!mo_m & !noise_m & valid_ph); if (length(coh_idx) < 5L) coh_idx <- which(!noise_m & valid_ph)
  sigma_cohort <- if (length(coh_idx) >= 5L) .circular_sd_days_from_residuals(ddoy_vec[coh_idx]) else NA_real_
  u_taxa <- unique(sp_vec[!is.na(sp_vec) & nzchar(sp_vec)])
  tau_eco_map <- stats::setNames(rep(defs$tau_ecotype, length(u_taxa)), u_taxa); tau_cry_map <- stats::setNames(rep(defs$tau_cryptic, length(u_taxa)), u_taxa); sigma_map <- stats::setNames(rep(sigma_cohort, length(u_taxa)), u_taxa)
  if (is.finite(sigma_cohort)) { tau_eco_map[] <- min(28.0, max(10.0, 2.0 * sigma_cohort)); tau_cry_map[] <- min(35.0, max(14.0, 2.5 * sigma_cohort)) }
  for (tx in u_taxa) {
    idx_t <- which(sp_vec == tx); tx_pool <- idx_t[clean_anch[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) < 5L) tx_pool <- idx_t[!mo_m[idx_t] & !noise_m[idx_t] & valid_ph[idx_t]]; if (length(tx_pool) < 5L) tx_pool <- idx_t[!noise_m[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) >= 5L && is.finite(sig_tx <- .circular_sd_days_from_residuals(ddoy_vec[tx_pool]))) { sigma_map[[tx]] <- sig_tx; tau_eco_map[[tx]] <- min(28.0, max(10.0, 2.0 * sig_tx)); tau_cry_map[[tx]] <- min(35.0, max(14.0, 2.5 * sig_tx)) }
  }
  fallback_eco <- if (is.finite(sigma_cohort)) min(28.0, max(10.0, 2.0 * sigma_cohort)) else defs$tau_ecotype
  fallback_cry <- if (is.finite(sigma_cohort)) min(35.0, max(14.0, 2.5 * sigma_cohort)) else defs$tau_cryptic
  list(tau_ecotype = as.numeric(ifelse(sp_vec %in% names(tau_eco_map), tau_eco_map[sp_vec], fallback_eco)),
       tau_cryptic = as.numeric(ifelse(sp_vec %in% names(tau_cry_map), tau_cry_map[sp_vec], fallback_cry)), sigma_by_taxon = sigma_map)
}

compute_hopkins_adjusted_phenology <- function(df) {
  n <- nrow(df); if (n == 0L) return(df)
  is_anthesis <- resolve_active_anthesis_mask(df); df$is_active_anthesis <- is_anthesis
  doy_raw <- .num_col(df, "doy", n); has_valid_doy <- is.finite(doy_raw) & (doy_raw >= 1.0) & (doy_raw <= 366.0)
  anchor_mask <- .resolve_step07_anchor_mask(df); mo_m <- .lgl_col(df, "morphometric_outlier_flag", n); noise_m <- .lgl_col(df, "is_confident_label_noise", n)
  clean_anth_doy <- is_anthesis & has_valid_doy & anchor_mask & !mo_m & !noise_m
  doy_pool <- doy_raw[clean_anth_doy]; if (length(doy_pool) == 0L) doy_pool <- doy_raw[is_anthesis & has_valid_doy & !noise_m]; if (length(doy_pool) == 0L) doy_pool <- doy_raw[has_valid_doy]
  if (length(doy_pool) > 0L) { th_pool <- 2.0 * pi * (doy_pool / 365.25); ref_theta <- atan2(mean(sin(th_pool)), mean(cos(th_pool))); if (ref_theta < 0.0) ref_theta <- ref_theta + 2.0 * pi; ref_doy <- ref_theta * (365.25 / (2.0 * pi)) } else { ref_doy <- 120.0; ref_theta <- 2.0 * pi * (ref_doy / 365.25) }
  doy_vec <- ifelse(has_valid_doy, doy_raw, ref_doy)
  is_cent_h <- .is_step07_centroid_coord(df); lat_vec <- .num_col(df, "decimalLatitude", n)
  elev_col <- intersect(c("elevation_m", "elevation", "minimumElevationInMeters"), names(df)); elev_vec <- if (length(elev_col) > 0L) .num_col(df, elev_col[1L], n) else rep(NA_real_, n)
  mat_col <- intersect(c("bio1_temp_mean", "bio1", "mat_c"), names(df)); mat_vec <- if (length(mat_col) > 0L) .num_col(df, mat_col[1L], n) else rep(NA_real_, n)
  lat_vec[is_cent_h] <- NA_real_; elev_vec[is_cent_h] <- NA_real_; mat_vec[is_cent_h] <- NA_real_
  pick_pool <- function(v) { p <- v[clean_anth_doy & is.finite(v)]; if (length(p) == 0L) p <- v[is_anthesis & is.finite(v)]; if (length(p) == 0L) p <- v[is.finite(v)]; p }
  ref_lat <- if (length(lp <- pick_pool(lat_vec)) > 0L) stats::median(lp, na.rm = TRUE) else 35.0
  ref_elev <- if (length(ep <- pick_pool(elev_vec)) > 0L) stats::median(ep, na.rm = TRUE) else 150.0
  ref_mat <- if (length(mp <- pick_pool(mat_vec)) > 0L) stats::median(mp, na.rm = TRUE) else 14.5
  lat_shift <- ifelse(is.finite(lat_vec), 4.0 * (lat_vec - ref_lat), 0.0); elev_shift <- ifelse(is.finite(elev_vec), (4.0 / 121.92) * (elev_vec - ref_elev), 0.0); mat_shift <- ifelse(is.finite(mat_vec), -1.5 * (mat_vec - ref_mat), 0.0)
  has_geo_cov <- is.finite(lat_vec) | is.finite(elev_vec); has_mat_cov <- is.finite(mat_vec); has_any_cov <- has_geo_cov | has_mat_cov
  theo_shift <- ifelse(has_geo_cov & has_mat_cov, 0.70 * (lat_shift + elev_shift) + 0.30 * mat_shift, ifelse(has_geo_cov, lat_shift + elev_shift, ifelse(has_mat_cov, mat_shift, 0.0)))
  theta_obs <- 2.0 * pi * (doy_vec / 365.25); theta_theo <- ref_theta + 2.0 * pi * (theo_shift / 365.25); sin_theo <- sin(theta_theo); cos_theo <- cos(theta_theo)
  fit_mask <- is_anthesis & has_valid_doy & has_any_cov & anchor_mask & !mo_m & !noise_m
  if (sum(fit_mask) < 5L) fit_mask <- is_anthesis & has_valid_doy & has_any_cov & anchor_mask & !noise_m
  if (sum(fit_mask) < 5L) fit_mask <- is_anthesis & has_valid_doy & has_any_cov & !noise_m
  cov_raw <- list(decimalLatitude = lat_vec, elevation = elev_vec, bio1_temp_mean = mat_vec)
  cov_terms <- names(cov_raw)[vapply(cov_raw, function(v) { vf <- v[fit_mask & is.finite(v)]; length(vf) >= 3L && is.finite(stats::sd(vf)) && stats::sd(vf) > 1e-6 }, logical(1))]
  alpha_eb <- if (length(cov_terms) > 0L) min(1.0, max(0.0, (sum(fit_mask) - 5.0) / 15.0)) else 0.0; theta_blend <- theta_theo; theta_blend_ref <- ref_theta
  if (alpha_eb > 0.0 && length(cov_terms) > 0L) {
    w_auth <- resolve_step07_ecological_authority_weights(df, allow_bronze_prior = TRUE); if (length(w_auth) != n) w_auth <- rep(1.0, n); w_auth[!is.finite(w_auth) | w_auth < 0.0] <- 0.25
    w_spb <- .num_col(df, "spatial_balance_weight", n); w_spb[!is.finite(w_spb) | w_spb <= 0.0] <- 1.0
    w_anth <- .num_col(df, "anthesis_dominance_weight", n); w_anth[!is.finite(w_anth) | w_anth <= 0.0] <- 1.0
    w_hop <- w_auth * w_spb * pmax(0.20, w_anth); if (!any(w_hop[fit_mask] > 0.0)) w_hop[fit_mask] <- 1.0
    hop_df <- data.frame(sin_doy = sin(theta_obs), cos_doy = cos(theta_obs), decimalLatitude = ifelse(is.finite(lat_vec), lat_vec, ref_lat),
                         elevation = ifelse(is.finite(elev_vec), elev_vec, ref_elev), bio1_temp_mean = ifelse(is.finite(mat_vec), mat_vec, ref_mat), .w_hop = w_hop)
    ref_df <- data.frame(decimalLatitude = ref_lat, elevation = ref_elev, bio1_temp_mean = ref_mat); rhs_str <- paste(cov_terms, collapse = " + ")
    lm_sin <- tryCatch(stats::lm(stats::as.formula(paste("sin_doy ~", rhs_str)), data = hop_df[fit_mask, , drop = FALSE], weights = .w_hop), error = function(e) NULL)
    lm_cos <- tryCatch(stats::lm(stats::as.formula(paste("cos_doy ~", rhs_str)), data = hop_df[fit_mask, , drop = FALSE], weights = .w_hop), error = function(e) NULL)
    if (!is.null(lm_sin) && !is.null(lm_cos)) {
      pred_s <- tryCatch(as.numeric(stats::predict(lm_sin, newdata = hop_df)), error = function(e) rep(NA_real_, n)); pred_c <- tryCatch(as.numeric(stats::predict(lm_cos, newdata = hop_df)), error = function(e) rep(NA_real_, n))
      ref_s <- tryCatch(as.numeric(stats::predict(lm_sin, newdata = ref_df)), error = function(e) NA_real_); ref_c <- tryCatch(as.numeric(stats::predict(lm_cos, newdata = ref_df)), error = function(e) NA_real_)
      if (all(is.finite(pred_s)) && all(is.finite(pred_c)) && is.finite(ref_s) && is.finite(ref_c)) {
        th_e <- atan2(pred_s, pred_c); th_er <- atan2(ref_s, ref_c)
        theta_blend <- atan2(alpha_eb * sin(th_e) + (1.0 - alpha_eb) * sin_theo, alpha_eb * cos(th_e) + (1.0 - alpha_eb) * cos_theo)
        theta_blend_ref <- atan2(alpha_eb * sin(th_er) + (1.0 - alpha_eb) * sin(ref_theta), alpha_eb * cos(th_er) + (1.0 - alpha_eb) * cos(ref_theta))
      }
    }
  }
  delta_theta_env <- ifelse(has_any_cov, atan2(sin(theta_blend - theta_blend_ref), cos(theta_blend - theta_blend_ref)), 0.0)
  theta_adj <- theta_obs - delta_theta_env; shift_days <- delta_theta_env * (365.25 / (2.0 * pi))
  adj_doy_calc <- ifelse(has_valid_doy, round(((doy_vec - shift_days - 1.0) %% 365.25) + 1.0, 2), NA_real_)
  adj_in <- .num_col(df, "hopkins_adjusted_doy", n); has_valid_adj <- is.finite(adj_in) & (adj_in >= 1.0) & (adj_in <= 366.0); has_pheno_doy <- has_valid_doy | has_valid_adj
  df$hopkins_adjusted_doy <- ifelse(has_valid_adj, round(adj_in, 2), adj_doy_calc)
  df$sin_doy <- ifelse(has_valid_doy, round(sin(theta_obs), 4), NA_real_); df$cos_doy <- ifelse(has_valid_doy, round(cos(theta_obs), 4), NA_real_)
  theta_adj_use <- ifelse(is.finite(df$hopkins_adjusted_doy), 2.0 * pi * (df$hopkins_adjusted_doy / 365.25), theta_adj)
  circ_resid_days <- numeric(n); sp_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unassigned_Packera", n)
  t0_hop <- if (any(.has_step07_geno_row_mask(df))) { if (exists("get_tier0_rectified_species_and_evictions", mode = "function")) get_tier0_rectified_species_and_evictions(df) else if (exists("resolve_tier0_anchor_overrides", mode = "function")) resolve_tier0_anchor_overrides(df) else NULL } else NULL
  for (tx in unique(sp_vec[!is.na(sp_vec)])) {
    idx_tx <- which(sp_vec == tx)
    if (!is.null(t0_hop)) {
      ev_tx <- t0_hop$evict_from_given_anchor & (!t0_hop$is_tier0_dna_anchor | t0_hop$rectified_taxon != tx)
      rect_tx <- which(t0_hop$is_tier0_dna_anchor & !is.na(t0_hop$rectified_taxon) & (t0_hop$rectified_taxon == tx) & is_anthesis & has_pheno_doy & !noise_m)
      anth_tx <- unique(c(idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & anchor_mask[idx_tx] & !mo_m[idx_tx] & !noise_m[idx_tx] & !ev_tx[idx_tx]], rect_tx))
      if (length(anth_tx) < 2L) anth_tx <- unique(c(idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & !noise_m[idx_tx] & !ev_tx[idx_tx]], rect_tx))
      if (length(anth_tx) == 0L) anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx]]
      w_h <- if (length(anth_tx) > 0L) ifelse(t0_hop$is_tier0_dna_anchor[anth_tx] & (t0_hop$rectified_taxon[anth_tx] == tx), ifelse(is.finite(t0_hop$tier0_weight[anth_tx]), t0_hop$tier0_weight[anth_tx], 4.0), 1.0) else numeric(0)
      theta_tx_c <- if (length(anth_tx) > 0L) atan2(sum(w_h * sin(theta_adj_use[anth_tx]), na.rm = TRUE), sum(w_h * cos(theta_adj_use[anth_tx]), na.rm = TRUE)) else ref_theta
    } else {
      anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & anchor_mask[idx_tx] & !mo_m[idx_tx] & !noise_m[idx_tx]]
      if (length(anth_tx) < 2L) anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & !noise_m[idx_tx]]
      if (length(anth_tx) == 0L) anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx]]
      theta_tx_c <- if (length(anth_tx) > 0L) atan2(mean(sin(theta_adj_use[anth_tx]), na.rm = TRUE), mean(cos(theta_adj_use[anth_tx]), na.rm = TRUE)) else ref_theta
    }
    d_th <- theta_adj_use[idx_tx] - theta_tx_c; circ_resid_days[idx_tx] <- atan2(sin(d_th), cos(d_th)) * (365.25 / (2.0 * pi))
  }
  circ_resid_gated <- ifelse(is_anthesis & has_pheno_doy, round(circ_resid_days, 2), NA_real_)
  ddoy_in <- .coalesce_step07_ddoy_resid(df, n); ddoy_use <- ifelse(is_anthesis & is.finite(ddoy_in), ddoy_in, circ_resid_gated)
  ddoy_use[!is_anthesis] <- NA_real_; df$delta_doy_resid <- ddoy_use; df$delta_doy <- ddoy_use; df
}

if (!exists("compute_probabilistic_consensus_scores", envir = environment(), inherits = FALSE, mode = "function") && !isTRUE(get0(".step07_sourcing_consensus_weights", envir = environment(), inherits = FALSE, ifnotfound = FALSE))) {
  .cw_dirs <- dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile))))
  for (.cw_p in unique(c("scripts/analysis/step07_consensus_weights_utils.R", file.path(.cw_dirs[!is.na(.cw_dirs) & nzchar(.cw_dirs)], "step07_consensus_weights_utils.R")))) {
    if (!is.na(.cw_p) && nzchar(.cw_p) && file.exists(.cw_p)) { source(.cw_p, local = environment()); break }
  }
}
