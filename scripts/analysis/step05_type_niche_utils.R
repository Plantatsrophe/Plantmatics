# ==============================================================================
# Step 05 Module: Primary Predictor Selection, Specimen Z & Type Niche Controls
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.to_bool_s05 <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & is.finite(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "1", "t", "yes")
}

.resolve_sub_wts_s05 <- function(vs, rf_var_weights, unhit_mult = 0.25) {
  if (length(vs) == 0L) return(numeric(0))
  wo <- stats::setNames(rep(1.0 / length(vs), length(vs)), vs)
  if (!is.null(rf_var_weights) && length(rf_var_weights) > 0L && !is.null(names(rf_var_weights))) {
    hit <- intersect(vs, names(rf_var_weights))
    if (length(hit) > 0L) {
      mp <- stats::median(as.numeric(rf_var_weights[rf_var_weights > 0]), na.rm = TRUE)
      if (!is.finite(mp) || mp <= 0.0) mp <- 1.0 / length(vs)
      wr <- if (unhit_mult > 0.0) stats::setNames(rep(unhit_mult * mp, length(vs)), vs) else wo
      wh <- suppressWarnings(as.numeric(rf_var_weights[hit])); wh[!is.finite(wh) | wh <= 0.0] <- max(unhit_mult * mp, 1e-4)
      wr[hit] <- wh; wo <- wr / sum(wr)
    }
  }
  wo
}

.fallback_is_homotypic_s05 <- function(bas_vec, acc_vec, df = NULL) {
  n <- length(bas_vec); out <- rep(TRUE, n)
  cs <- if (!is.null(df) && "concept_synonymy_status" %in% names(df)) tolower(trimws(as.character(df$concept_synonymy_status))) else rep("", n)
  cs[is.na(cs)] <- ""
  stem_ep <- function(nm) {
    toks <- unlist(strsplit(tolower(trimws(gsub("[^[:alpha:][:space:]]", " ", as.character(nm)))), "[[:space:]]+"))
    toks <- toks[nzchar(toks)]
    if (length(toks) < 2L) return("")
    rk_i <- which(toks %in% c("var", "subsp", "ssp", "f", "forma"))
    ep <- if (length(rk_i) > 0L && rk_i[1L] < length(toks)) toks[rk_i[1L] + 1L] else toks[2L]
    sub("(orum|arum|ensis|ense|ibus|ium|ius|ia|is|us|um|er|ra|rum|ae|ii|a|e|i|o)$", "", ep)
  }
  for (i in seq_len(n)) {
    if (grepl("heterotypic", cs[i])) { out[i] <- FALSE; next }
    if (grepl("homotypic|accepted_exact", cs[i])) { out[i] <- TRUE; next }
    b_i <- trimws(as.character(bas_vec[i])); a_i <- trimws(as.character(acc_vec[i]))
    if (!is.na(b_i) && nzchar(b_i) && !is.na(a_i) && nzchar(a_i)) {
      sb <- stem_ep(b_i); sa <- stem_ep(a_i)
      if (nzchar(sb) && nzchar(sa) && !identical(sb, sa) && !startsWith(sb, sa) && !startsWith(sa, sb)) out[i] <- FALSE
    }
  }
  out
}

.fallback_topotype_rows_s05 <- function(df, cand_rows, tx_k) {
  if (is.null(df) || length(cand_rows) == 0L || !("stateProvince" %in% names(df))) return(integer(0))
  lc <- if ("locus_classicus_text" %in% names(df)) { v <- trimws(as.character(df$locus_classicus_text[cand_rows])); v <- v[!is.na(v) & nzchar(v)]; if (length(v) > 0L) v[1L] else "" } else ""
  if (!nzchar(lc) && identical(tx_k, "Packera dubia")) lc <- "Flat Rock, Kershaw County, South Carolina"
  if (!nzchar(lc)) return(integer(0))
  st_map <- c("alabama"="al","arkansas"="ar","florida"="fl","georgia"="ga","kentucky"="ky","louisiana"="la","maryland"="md","mississippi"="ms","north carolina"="nc","south carolina"="sc","tennessee"="tn","texas"="tx","virginia"="va","west virginia"="wv")
  lc_l <- tolower(gsub("[[:space:]]+", " ", lc)); lc_st <- ""
  for (sn in names(st_map)) if (grepl(paste0("(^|[^[:alpha:]])", sn, "($|[^[:alpha:]])"), lc_l)) { lc_st <- st_map[[sn]]; break }
  mc <- regmatches(lc_l, regexpr("[a-z][a-z.-]*[[:space:]]+(county|parish)", lc_l))
  lc_c <- if (length(mc) > 0L && nzchar(mc[1L])) sub("[[:space:]]+(county|parish)$", "", mc[1L]) else ""
  sr <- tolower(trimws(as.character(df$stateProvince[cand_rows]))); sr[is.na(sr)] <- ""
  sn_v <- ifelse(nchar(sr) == 2L, sr, ifelse(sr %in% names(st_map), st_map[sr], sr))
  cr <- if ("county" %in% names(df)) sub("[[:space:]]+(county|parish|co\\.?)$", "", tolower(trimws(as.character(df$county[cand_rows])))) else rep("", length(cand_rows))
  cr[is.na(cr)] <- ""; qv <- if ("is_centroid_or_inst_coord" %in% names(df)) .to_bool_s05(df$is_centroid_or_inst_coord[cand_rows]) else rep(FALSE, length(cand_rows))
  ok_s <- !qv & nzchar(lc_st) & (sn_v == lc_st); if (!any(ok_s)) return(integer(0))
  ch <- if (nzchar(lc_c)) cand_rows[ok_s & (cr == lc_c)] else integer(0)
  if (length(ch) >= 2L) ch else if (sum(ok_s) >= 2L) cand_rows[ok_s] else ch
}

select_primary_predictors_from_rf <- function(rf_res = NULL, vif_res = NULL, df_env = NULL, cum_imp_target = NULL, intraspecific_mode = FALSE) {
  if (is.null(cum_imp_target) || !is.finite(cum_imp_target)) {
    cum_imp_target <- if (exists("load_step05_macroecology_config", mode = "function")) {
      v <- tryCatch(load_step05_macroecology_config()$cumulative_importance_cutoff, error = function(e) NA_real_)
      if (is.finite(v) && v > 0.0 && v <= 1.0) as.numeric(v) else 0.85
    } else 0.85
  }
  clim_pat <- "^(bio[0-9]*_|summer_moisture_|maritime_continentality_|delta_doy$|doy$|hopkins_|sin_doy$|cos_doy$)"
  topo_pat <- "^(topo_|elevation($|_)|elev_)"

  cand_vars <- if (!is.null(vif_res) && is.list(vif_res) && length(if (!is.null(vif_res$kept_vars)) vif_res$kept_vars else vif_res$selected_variables) > 0L)
    as.character(if (!is.null(vif_res$kept_vars)) vif_res$kept_vars else vif_res$selected_variables)
    else if (!is.null(rf_res) && is.list(rf_res) && !is.null(rf_res$importance)) as.character(rf_res$importance$Predictor)
    else if (!is.null(df_env) && is.data.frame(df_env)) discover_candidate_eco_predictors(df_env, include_pheno = FALSE, min_finite = 3L) else character(0)

  audit_excl <- c("soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct", "ultramafic_calcareous_gradient",
                  "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m",
                  "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m", "nwpl_indicator_code", "nwpl_wetland_score")
  cand_vars <- setdiff(unique(cand_vars[!is.na(cand_vars) & nzchar(cand_vars) & !grepl("^MEM[0-9]*$|^polaris_", cand_vars, ignore.case = TRUE)]), audit_excl)
  if (!is.null(df_env) && is.data.frame(df_env)) {
    cand_vars <- intersect(cand_vars, names(df_env))
    cand_vars <- cand_vars[vapply(cand_vars, function(cl) {
      v <- suppressWarnings(as.numeric(df_env[[cl]])); fin <- is.finite(v)
      sum(fin) >= min(3L, max(2L, nrow(df_env))) && is.finite(stats::sd(v[fin])) && stats::sd(v[fin]) > 1e-8
    }, logical(1))]
  }
  if (length(cand_vars) == 0L) return(list(primary_env_vars = character(0), primary_clim_vars = character(0), primary_edaph_vars = character(0), primary_topo_vars = character(0), rf_var_weights = numeric(0), selection_method = "empty"))

  mda_map <- stats::setNames(rep(NA_real_, length(cand_vars)), cand_vars)
  if (!is.null(rf_res) && is.list(rf_res) && !is.null(rf_res$importance) && all(c("Predictor", "MeanDecreaseAccuracy") %in% names(rf_res$importance))) {
    hm <- match(cand_vars, as.character(rf_res$importance$Predictor)); ok_m <- !is.na(hm)
    mda_map[ok_m] <- suppressWarnings(as.numeric(rf_res$importance$MeanDecreaseAccuracy[hm[ok_m]]))
  }
  pos_vars <- names(mda_map[is.finite(mda_map) & mda_map > 0.0]); sel_method <- "rf_permutation_85pct"
  if (length(pos_vars) >= 1L) {
    ord_pos <- pos_vars[order(mda_map[pos_vars], decreasing = TRUE)]; pv <- mda_map[ord_pos]
    k_cut <- which((cumsum(pv) / max(sum(pv), 1e-12)) >= cum_imp_target)[1L]; if (is.na(k_cut)) k_cut <- length(ord_pos)
    k_keep <- min(length(ord_pos), max(min(2L, length(ord_pos)), as.integer(k_cut)))
    selected_vars <- ord_pos[seq_len(k_keep)]; raw_weights <- pv[seq_len(k_keep)]
  } else {
    sel_method <- "univariate_bh_fdr_fallback"
    anc_grp <- if (!is.null(df_env) && is.data.frame(df_env)) resolve_active_anchor_groups(df_env, intraspecific_mode = intraspecific_mode) else character(0)
    p_vals <- stats::setNames(rep(1.0, length(cand_vars)), cand_vars)
    sep_scores <- if (!is.null(df_env) && is.data.frame(df_env) && exists("compute_anchor_separation_scores", mode = "function")) compute_anchor_separation_scores(df_env[, cand_vars, drop = FALSE], anchor_groups = anc_grp) else stats::setNames(seq(length(cand_vars), 1L), cand_vars)
    if (!is.null(df_env) && is.data.frame(df_env) && length(anc_grp) == nrow(df_env)) {
      for (cl in cand_vars) {
        v <- suppressWarnings(as.numeric(df_env[[cl]])); ok <- is.finite(v) & !is.na(anc_grp) & nzchar(anc_grp)
        if (sum(ok) >= 4L && length(unique(anc_grp[ok])) >= 2L) {
          fg <- factor(anc_grp[ok])
          pv <- suppressWarnings(min(c(tryCatch(suppressWarnings(as.numeric(stats::kruskal.test(v[ok] ~ fg)$p.value)), error = function(e) NA_real_),
                                       tryCatch(suppressWarnings(as.numeric(summary(stats::aov(v[ok] ~ fg))[[1]][["Pr(>F)"]][1L])), error = function(e) NA_real_)), na.rm = TRUE))
          if (is.finite(pv)) p_vals[[cl]] <- pv
        }
      }
    }
    p_adj <- stats::p.adjust(p_vals, method = "BH"); sig_vars <- names(p_adj[is.finite(p_adj) & p_adj < 0.05])
    if (length(sig_vars) >= 1L) selected_vars <- sig_vars[order(sep_scores[sig_vars], decreasing = TRUE)]
    else { sel_method <- "anchor_separation_top_fallback"; ord_s <- names(sort(sep_scores, decreasing = TRUE)); selected_vars <- ord_s[seq_len(min(length(ord_s), max(2L, min(6L, length(ord_s)))))] }
    raw_weights <- stats::setNames(pmax(as.numeric(sep_scores[selected_vars]), 1e-4), selected_vars)
  }

  primary_clim <- grep(clim_pat, selected_vars, value = TRUE, ignore.case = TRUE)
  primary_topo <- setdiff(grep(topo_pat, selected_vars, value = TRUE, ignore.case = TRUE), primary_clim)
  primary_edaph <- setdiff(selected_vars, unique(c(primary_clim, primary_topo)))
  cand_clim_pos <- intersect(grep(clim_pat, cand_vars, value = TRUE, ignore.case = TRUE), pos_vars)
  if (length(primary_clim) == 0L && length(cand_clim_pos) > 0L) { bc <- cand_clim_pos[which.max(mda_map[cand_clim_pos])][1L]; primary_clim <- bc; selected_vars <- unique(c(selected_vars, bc)); raw_weights[[bc]] <- max(mda_map[[bc]], 1e-4) }
  cand_ed_pos <- intersect(setdiff(cand_vars, unique(c(grep(clim_pat, cand_vars, value = TRUE, ignore.case = TRUE), grep(topo_pat, cand_vars, value = TRUE, ignore.case = TRUE)))), pos_vars)
  if (length(primary_edaph) == 0L && length(cand_ed_pos) > 0L) { be <- cand_ed_pos[which.max(mda_map[cand_ed_pos])][1L]; primary_edaph <- be; selected_vars <- unique(c(selected_vars, be)); raw_weights[[be]] <- max(mda_map[[be]], 1e-4) }
  raw_weights <- raw_weights[selected_vars]; raw_weights[!is.finite(raw_weights) | raw_weights <= 0.0] <- 1e-4
  list(primary_env_vars = selected_vars, primary_clim_vars = primary_clim, primary_edaph_vars = primary_edaph,
       primary_topo_vars = primary_topo, rf_var_weights = stats::setNames(raw_weights / sum(raw_weights), selected_vars), selection_method = sel_method)
}

compute_specimen_edaphic_z <- function(df, primary_edaph_vars = NULL, primary_clim_vars = NULL, primary_env_vars = NULL, rf_var_weights = NULL, intraspecific_mode = FALSE) {
  n <- if (is.null(df) || !is.data.frame(df)) 0L else nrow(df)
  if (n == 0L) return(list(edaphic_z_score = numeric(0), climatic_z_score = numeric(0), env_z_score = numeric(0), primary_edaph_vars = character(0), primary_clim_vars = character(0), primary_env_vars = character(0), eval_specimen_z = function(idx, target_tx, mode = "edaphic") 0.0))
  is_intra <- isTRUE(intraspecific_mode) || (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  if (is_intra && !("analysis_group" %in% names(df))) df <- resolve_analysis_group(df, intraspecific_mode = TRUE, primary_env_vars = primary_env_vars, rf_var_weights = rf_var_weights)
  grp_vec <- if (is_intra && "analysis_group" %in% names(df)) as.character(df$analysis_group) else as.character(df$species_standardized)
  active_taxa <- sort(unique(grp_vec[!is.na(grp_vec) & nzchar(grp_vec) & grp_vec != "Unknown"]))
  clim_pat <- "^(bio[0-9]*_|summer_moisture_|maritime_continentality_)"
  cand_all <- discover_candidate_eco_predictors(df, include_pheno = FALSE, min_finite = 3L)

  if (is.null(primary_env_vars) || length(primary_env_vars) == 0L) { primary_env_vars <- unique(c(as.character(primary_edaph_vars), as.character(primary_clim_vars))); if (length(primary_env_vars) == 0L) primary_env_vars <- cand_all }
  primary_env_vars <- intersect(unique(as.character(primary_env_vars)), names(df))
  primary_env_vars <- primary_env_vars[vapply(primary_env_vars, function(cl) sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= min(3L, max(2L, n)), logical(1))]
  primary_clim_vars <- if (is.null(primary_clim_vars) || length(primary_clim_vars) == 0L) { c0 <- grep(clim_pat, primary_env_vars, value = TRUE, ignore.case = TRUE); if (length(c0) > 0L) c0 else grep(clim_pat, cand_all, value = TRUE, ignore.case = TRUE) } else intersect(unique(as.character(primary_clim_vars)), names(df))
  primary_edaph_vars <- if (is.null(primary_edaph_vars) || length(primary_edaph_vars) == 0L) { e0 <- setdiff(primary_env_vars, grep(clim_pat, primary_env_vars, value = TRUE, ignore.case = TRUE)); if (length(e0) > 0L) e0 else setdiff(cand_all, grep(clim_pat, cand_all, value = TRUE, ignore.case = TRUE)) } else intersect(unique(as.character(primary_edaph_vars)), names(df))
  all_z_vars <- unique(c(primary_env_vars, primary_edaph_vars, primary_clim_vars))
  all_z_vars <- all_z_vars[vapply(all_z_vars, function(cl) (cl %in% names(df)) && sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= min(3L, max(2L, n)), logical(1))]

  auth_w <- if (exists("resolve_ecological_authority_weights", mode = "function")) resolve_ecological_authority_weights(df) else rep(1.0, n)
  auth_w[!is.finite(auth_w) | auth_w < 0.0] <- 0.0
  spat_w <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df) else rep(1.0, n)
  spat_w[!is.finite(spat_w) | spat_w <= 0.0] <- 1.0

  elig_coord <- if ("is_training_eligible_coord" %in% names(df)) { e <- as.logical(df$is_training_eligible_coord); e[is.na(e)] <- TRUE; e } else rep(TRUE, n)
  det_tier <- if ("determiner_tier" %in% names(df)) ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier)) else rep("Tier_1_Gold", n)
  is_t_z <- if ("is_type_specimen" %in% names(df)) .to_bool_s05(df$is_type_specimen) else rep(FALSE, n)
  is_prom <- if ("is_promoted_bronze" %in% names(df)) .to_bool_s05(df$is_promoted_bronze) else rep(FALSE, n)
  cl_z <- if ("is_confident_label_noise" %in% names(df)) .to_bool_s05(df$is_confident_label_noise) else rep(FALSE, n)
  open_z <- if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) { m <- as.logical(.identify_open_or_hybrid_concept_mask(df)); m[is.na(m)] <- FALSE; m } else rep(FALSE, n)
  cent_z <- if ("is_centroid_or_inst_coord" %in% names(df)) .to_bool_s05(df$is_centroid_or_inst_coord) else rep(FALSE, n)

  tier_ok_z <- (grepl("Tier_1|Gold|Tier_2|Silver", det_tier, ignore.case = TRUE) | is_prom | is_t_z) & (!cl_z | is_t_z) & !open_z & !cent_z
  row_wts <- ifelse(tier_ok_z, pmax(auth_w, 0.25), 0.0) * spat_w
  is_anchor <- tier_ok_z & elig_coord & (grp_vec %in% active_taxa)
  if (!any(tier_ok_z)) { is_anchor <- (!cl_z) & (!open_z) & (!cent_z) & elig_coord & (grp_vec %in% active_taxa); row_wts <- spat_w }

  w_mean_z <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) mean(x[is.finite(x)], na.rm = TRUE) else sum(x[ok] * w[ok]) / sum(w[ok]) }
  w_sd_z <- function(x, w, fl) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (sum(ok) < 2L) return(max(stats::sd(x[is.finite(x)], na.rm = TRUE), fl, na.rm = TRUE))
    wm <- sum(x[ok] * w[ok]) / sum(w[ok]); dn <- sum(w[ok]) - (sum(w[ok]^2) / sum(w[ok]))
    max(sqrt(max(if (dn > 1e-6) sum(w[ok] * (x[ok] - wm)^2) / dn else stats::var(x[ok]), 0.0, na.rm = TRUE)), fl, na.rm = TRUE)
  }
  scale_floors <- stats::setNames(vapply(all_z_vars, function(cl) {
    va <- suppressWarnings(as.numeric(df[[cl]])); vf <- va[is_anchor]; vf <- vf[is.finite(vf)]; if (length(vf) < 2L) vf <- va[is.finite(va)]
    if (length(vf) < 2L) return(1e-3)
    dg <- max(c(stats::sd(vf), stats::IQR(vf) / 1.349), na.rm = TRUE); if (!is.finite(dg) || dg <= 1e-8) 1e-3 else max(0.15 * dg, 1e-4)
  }, numeric(1)), all_z_vars)

  z_profiles <- list()
  if (length(active_taxa) > 0L && length(all_z_vars) > 0L) {
    for (tx in active_taxa) {
      idx_tx <- which(is_anchor & grp_vec == tx); if (length(idx_tx) == 0L) idx_tx <- which(tier_ok_z & grp_vec == tx)
      if (length(idx_tx) == 0L) idx_tx <- which((!cl_z) & elig_coord & grp_vec == tx); if (length(idx_tx) == 0L) idx_tx <- which(grp_vec == tx)
      w_tx <- pmax(row_wts[idx_tx], 1e-4)
      z_profiles[[tx]] <- list(
        mu = stats::setNames(vapply(all_z_vars, function(cl) w_mean_z(suppressWarnings(as.numeric(df[[cl]][idx_tx])), w_tx), numeric(1)), all_z_vars),
        sd = stats::setNames(vapply(all_z_vars, function(cl) w_sd_z(suppressWarnings(as.numeric(df[[cl]][idx_tx])), w_tx, scale_floors[[cl]]), numeric(1)), all_z_vars)
      )
    }
    for (cl in all_z_vars) {
      bar_sd <- mean(vapply(z_profiles, function(pr) pr$sd[[cl]], numeric(1)), na.rm = TRUE)
      if (is.finite(bar_sd) && bar_sd > 0) for (tx in active_taxa) z_profiles[[tx]]$sd[[cl]] <- max(0.5 * z_profiles[[tx]]$sd[[cl]] + 0.5 * bar_sd, scale_floors[[cl]])
    }
  }

  edaph_use <- if (length(primary_edaph_vars) > 0L) intersect(primary_edaph_vars, all_z_vars) else intersect(primary_env_vars, all_z_vars)
  clim_use <- intersect(primary_clim_vars, all_z_vars)
  env_use <- if (length(primary_env_vars) > 0L) intersect(primary_env_vars, all_z_vars) else all_z_vars
  w_edaph <- .resolve_sub_wts_s05(edaph_use, rf_var_weights, 0.25)
  w_clim <- .resolve_sub_wts_s05(clim_use, rf_var_weights, 0.25)
  w_env <- .resolve_sub_wts_s05(env_use, rf_var_weights, 0.25)

  eval_specimen_z <- function(idx, target_tx, mode = "edaphic") {
    na_fb <- if (identical(mode, "climatic")) NA_real_ else 0.0
    if (is.na(target_tx) || !(as.character(target_tx) %in% names(z_profiles))) return(na_fb)
    vars_m <- if (identical(mode, "climatic")) clim_use else if (identical(mode, "env")) env_use else edaph_use
    wts_m <- if (identical(mode, "climatic")) w_clim else if (identical(mode, "env")) w_env else w_edaph
    if (length(vars_m) == 0L) return(na_fb)
    pr <- z_profiles[[as.character(target_tx)]]
    xr <- vapply(vars_m, function(cl) suppressWarnings(as.numeric(df[[cl]][idx])), numeric(1))
    z2 <- ((xr - pr$mu[vars_m]) / pmax(pr$sd[vars_m], scale_floors[vars_m], 1e-4))^2
    ok <- is.finite(z2) & is.finite(wts_m) & (wts_m > 0); if (!any(ok)) na_fb else sqrt(sum(wts_m[ok] * z2[ok]) / sum(wts_m[ok]))
  }
  list(edaphic_z_score = round(vapply(seq_len(n), function(i) eval_specimen_z(i, grp_vec[i], "edaphic"), numeric(1)), 4),
       climatic_z_score = round(vapply(seq_len(n), function(i) eval_specimen_z(i, grp_vec[i], "climatic"), numeric(1)), 4),
       env_z_score = round(vapply(seq_len(n), function(i) eval_specimen_z(i, grp_vec[i], "env"), numeric(1)), 4),
       primary_edaph_vars = edaph_use, primary_clim_vars = clim_use, primary_env_vars = env_use, rf_var_weights = w_env, eval_specimen_z = eval_specimen_z)
}

compute_type_locality_ecological_control <- function(df, primary_env_vars = NULL, primary_edaph_vars = NULL, rf_var_weights = NULL, intraspecific_mode = FALSE, max_radius_km = 150.0) {
  n <- if (is.null(df) || !is.data.frame(df)) 0L else nrow(df)
  if (n == 0L) {
    return(list(type_edaphic_z_score = numeric(0), nearest_type_niche_taxon = character(0), matches_type_lithology = logical(0),
                nearest_type_niche_basionym = character(0), nearest_type_basionym_edaphic_z = numeric(0), is_heterotypic_type_audit = logical(0),
                type_niche_status = character(0), type_niche_audit_status = character(0), rf_prediction_source_hint = character(0),
                type_niche_profiles = list(), basionym_niche_profiles = list()))
  }

  t_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else {
    ts_raw <- if ("type_status" %in% names(df)) as.character(df$type_status) else rep("", n); ts_raw[is.na(ts_raw)] <- ""
    is_t <- grepl("(^|[^[:alnum:]_])(holotype|lectotype|neotype|epitype|isotype|isolectotype|isoneotype|syntype|isosyntype|paratype|isoparatype|type|cotype)s?($|[^[:alnum:]_])", ts_raw, ignore.case = TRUE) &
            !grepl("not a type|non-type|false type", ts_raw, ignore.case = TRUE)
    if ("is_type_specimen" %in% names(df)) is_t <- is_t | .to_bool_s05(df$is_type_specimen)
    list(type_rank = ifelse(is_t, 2L, 0L), is_type_specimen = is_t)
  }
  is_type_vec <- as.logical(t_info$is_type_specimen); is_type_vec[is.na(is_type_vec)] <- FALSE
  is_intra <- isTRUE(intraspecific_mode) || (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  grp_vec <- if (is_intra && "analysis_group" %in% names(df)) as.character(df$analysis_group)
    else if ("species_standardized" %in% names(df)) as.character(df$species_standardized)
    else if ("analysis_group" %in% names(df)) as.character(df$analysis_group) else rep("Unknown", n)
  active_taxa <- sort(unique(grp_vec[!is.na(grp_vec) & nzchar(grp_vec) & grp_vec != "Unknown"]))
  if (length(active_taxa) == 0L) active_taxa <- "Unknown"

  if ((is.null(primary_env_vars) || length(primary_env_vars) == 0L) && "primary_env_predictors" %in% names(df)) {
    r0 <- stats::na.omit(as.character(df$primary_env_predictors)); if (length(r0) > 0L && nzchar(r0[1L])) primary_env_vars <- trimws(unlist(strsplit(r0[1L], ";", fixed = TRUE)))
  }
  if ((is.null(primary_edaph_vars) || length(primary_edaph_vars) == 0L) && "primary_edaph_predictors" %in% names(df)) {
    r1 <- stats::na.omit(as.character(df$primary_edaph_predictors)); if (length(r1) > 0L && nzchar(r1[1L])) primary_edaph_vars <- trimws(unlist(strsplit(r1[1L], ";", fixed = TRUE)))
  }
  cand_all <- discover_candidate_eco_predictors(df, include_pheno = FALSE, min_finite = 2L)
  env_use <- if (!is.null(primary_env_vars) && length(primary_env_vars) > 0L) intersect(as.character(primary_env_vars), names(df)) else cand_all
  if (length(env_use) == 0L) env_use <- intersect(c("fine_soil_ph", "ilr_texture_1", "ilr_texture_2", "fine_soil_sand", "fine_soil_cec", "soil_ph", "soil_sand", "soil_cec", "bio1_temp_mean", "bio12_precip_annual"), names(df))
  clim_pat <- "^(bio[0-9]*_|summer_moisture_|maritime_continentality_)"
  edaph_use <- if (!is.null(primary_edaph_vars) && length(primary_edaph_vars) > 0L) intersect(as.character(primary_edaph_vars), names(df)) else setdiff(env_use, grep(clim_pat, env_use, value = TRUE, ignore.case = TRUE))
  if (length(edaph_use) == 0L) edaph_use <- env_use
  all_vars <- unique(c(env_use, edaph_use))
  all_vars <- all_vars[vapply(all_vars, function(cl) (cl %in% names(df)) && sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= 1L, logical(1))]

  lat_v <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else if ("latitude" %in% names(df)) suppressWarnings(as.numeric(df$latitude)) else rep(NA_real_, n)
  lon_v <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else if ("longitude" %in% names(df)) suppressWarnings(as.numeric(df$longitude)) else rep(NA_real_, n)
  has_xy <- is.finite(lat_v) & is.finite(lon_v)
  has_unc_col <- ("coordinateUncertaintyInMeters" %in% names(df)) || ("coordinateUncertainty" %in% names(df))
  unc_v <- if ("coordinateUncertaintyInMeters" %in% names(df)) suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters)) else if ("coordinateUncertainty" %in% names(df)) suppressWarnings(as.numeric(df$coordinateUncertainty)) else rep(NA_real_, n)
  yr_tc <- if ("year" %in% names(df)) suppressWarnings(as.numeric(df$year)) else rep(NA_real_, n)
  modern_gps_tc <- if (("year" %in% names(df)) && any(is.finite(yr_tc))) (is.finite(yr_tc) & yr_tc >= 2000) else !has_unc_col

  geo_elig <- if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) tryCatch(is_fine_edaphic_type_or_anchor_eligible(df, max_edaphic_unc_m = 10000.0, max_clim_unc_m = 50000.0), error = function(e) NULL) else NULL
  raw_cent <- if ("is_centroid_or_inst_coord" %in% names(df)) .to_bool_s05(df$is_centroid_or_inst_coord) else rep(FALSE, n)
  if (exists("detect_r_centroid_or_inst_quarantine", mode = "function")) {
    qa <- as.logical(detect_r_centroid_or_inst_quarantine(df, is_type_vec = rep(FALSE, n))$is_quarantined); qa[is.na(qa)] <- FALSE; raw_cent <- raw_cent | qa
  }
  over_10km <- is.finite(unc_v) & (unc_v > 10000.0)
  is_quar_type <- is_type_vec & (raw_cent | over_10km | (any(has_xy) & !has_xy))
  if (!is.null(geo_elig) && !is.null(geo_elig$is_quarantined_type) && length(geo_elig$is_quarantined_type) == n) is_quar_type <- is_quar_type | (is_type_vec & !is.na(geo_elig$is_quarantined_type) & as.logical(geo_elig$is_quarantined_type))
  is_open_hyb <- if (!is.null(geo_elig) && !is.null(geo_elig$is_open_hyb) && length(geo_elig$is_open_hyb) == n) (!is.na(geo_elig$is_open_hyb) & as.logical(geo_elig$is_open_hyb))
    else if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) { m <- as.logical(.identify_open_or_hybrid_concept_mask(df)); m[is.na(m)] <- FALSE; m } else rep(FALSE, n)
  is_type_georef_ok <- is_type_vec & !is_quar_type

  is_precise_edaph <- !raw_cent & !over_10km & ((is.finite(unc_v) & unc_v > 0.0 & unc_v <= 2000.0) | (!is.finite(unc_v) & modern_gps_tc))
  is_precise_clim <- !raw_cent & (!is.finite(unc_v) | (unc_v >= 0.0 & unc_v <= 10000.0))
  reg_v <- if ("regional_group" %in% names(df)) trimws(as.character(df$regional_group)) else rep("", n)
  reg_v[is.na(reg_v) | tolower(reg_v) %in% c("na", "unknown", "unknown_ecoregion")] <- ""
  det_t <- if ("determiner_tier" %in% names(df)) ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier)) else rep("Tier_1_Gold", n)
  is_gold_v <- (grepl("tier_?[01]|gold|type", det_t, ignore.case = TRUE) | is_type_georef_ok) & !raw_cent & !over_10km & !is_open_hyb
  is_anc_v <- (is_gold_v | grepl("tier_?2|silver", det_t, ignore.case = TRUE) | (if ("is_promoted_bronze" %in% names(df)) .to_bool_s05(df$is_promoted_bronze) else rep(FALSE, n))) & !raw_cent & !over_10km & !is_open_hyb

  bas_all <- if (exists(".canonicalize_type_basionym_vec", mode = "function")) .canonicalize_type_basionym_vec(df, fallback_taxa = grp_vec)
    else if ("typified_basionym" %in% names(df)) { tb <- trimws(as.character(df$typified_basionym)); ifelse(!is.na(tb) & nzchar(tb), tb, grp_vec) } else grp_vec
  is_homo_all <- if (exists(".is_homotypic_sensu_stricto_basionym", mode = "function")) .is_homotypic_sensu_stricto_basionym(bas_all, grp_vec, df = df)
    else .fallback_is_homotypic_s05(bas_all, grp_vec, df = df)
  is_homo_all[is.na(is_homo_all)] <- TRUE

  haversine_km <- function(la1, lo1, la2, lo2) {
    r <- pi / 180.0; a <- sin((la2 - la1) * r / 2.0)^2 + cos(la1 * r) * cos(la2 * r) * sin((lo2 - lo1) * r / 2.0)^2
    6371.0 * (2.0 * atan2(sqrt(pmax(0.0, a)), sqrt(pmax(0.0, 1.0 - a))))
  }
  clean_scale_mask <- !raw_cent & !over_10km & !(is_type_vec & !is_homo_all)
  scale_floors <- stats::setNames(vapply(all_vars, function(cl) {
    v <- suppressWarnings(as.numeric(df[[cl]])); vf <- v[clean_scale_mask & is.finite(v)]; if (length(vf) < 2L) vf <- v[is.finite(v)]
    if (length(vf) < 2L) return(1.0)
    dg <- max(c(stats::mad(vf, na.rm = TRUE), stats::sd(vf, na.rm = TRUE), stats::IQR(vf, na.rm = TRUE) / 1.349), na.rm = TRUE)
    if (!is.finite(dg) || dg <= 1e-8) 1.0 else max(0.20 * dg, 1e-4)
  }, numeric(1)), all_vars)

  lith_col <- if ("sgmc_lithology_class" %in% names(df)) "sgmc_lithology_class" else if ("Dominant_Lithology_Class" %in% names(df)) "Dominant_Lithology_Class" else NULL
  lith_vec <- if (!is.null(lith_col)) trimws(as.character(df[[lith_col]])) else rep("", n)
  lith_vec[is.na(lith_vec) | tolower(lith_vec) %in% c("na", "unknown", "none", "unmapped")] <- ""
  edaph_v <- intersect(edaph_use, all_vars); env_v <- intersect(env_use, all_vars)
  w_ed <- .resolve_sub_wts_s05(edaph_v, rf_var_weights, 0.0); w_en <- .resolve_sub_wts_s05(env_v, rf_var_weights, 0.0)

  build_profile <- function(tx_rows, focal_types, excl_rows = integer(0), fixed_sd = NULL) {
    tx_clean <- setdiff(tx_rows, excl_rows)
    if (length(focal_types) > 0L) {
      gc <- setdiff(tx_clean[is_gold_v[tx_clean]], focal_types); near_g <- integer(0)
      if (length(gc) > 0L) {
        t_regs <- unique(reg_v[focal_types[nzchar(reg_v[focal_types])]])
        gc_e <- gc[is_precise_edaph[gc]]; if (length(gc_e) == 0L) gc_e <- gc
        n_xy <- integer(0); t_xy <- focal_types[has_xy[focal_types]]
        c_xy <- gc[has_xy[gc] & is_precise_edaph[gc]]; if (length(c_xy) == 0L) c_xy <- gc[has_xy[gc]]
        if (length(t_xy) > 0L && length(c_xy) > 0L) {
          md <- rep(Inf, length(c_xy)); for (ti in t_xy) md <- pmin(md, haversine_km(lat_v[ti], lon_v[ti], lat_v[c_xy], lon_v[c_xy]))
          n_xy <- c_xy[is.finite(md) & (md <= max_radius_km)]
        }
        gc_no_xy <- if (length(t_xy) > 0L) gc_e[!has_xy[gc_e]] else gc_e
        n_eco <- if (length(t_regs) > 0L) gc_no_xy[reg_v[gc_no_xy] %in% t_regs] else integer(0)
        if (length(n_eco) == 0L && length(t_xy) == 0L && length(t_regs) > 0L) n_eco <- gc[reg_v[gc] %in% t_regs]
        near_g <- unique(c(n_eco, n_xy))
      }
      loc_rows <- c(focal_types, near_g)
    } else {
      atx <- tx_clean[is_anc_v[tx_clean]]; loc_rows <- if (length(atx) >= 1L) atx else tx_clean
    }
    mu_tx <- stats::setNames(numeric(length(all_vars)), all_vars); sd_tx <- stats::setNames(numeric(length(all_vars)), all_vars)
    for (cl in all_vars) {
      cv <- suppressWarnings(as.numeric(df[[cl]])); is_ed <- cl %in% edaph_use
      tp <- if (is_ed) focal_types[is_precise_edaph[focal_types]] else focal_types[is_precise_clim[focal_types]]
      lp <- if (is_ed) loc_rows[is_precise_edaph[loc_rows]] else loc_rows[is_precise_clim[loc_rows]]
      if (length(lp) == 0L || !any(is.finite(cv[lp]))) lp <- loc_rows
      tv <- cv[tp]; tv <- tv[is.finite(tv)]; if (length(tv) == 0L && length(focal_types) > 0L) { tf <- cv[focal_types]; tv <- tf[is.finite(tf)] }
      lv <- cv[lp]; lv <- lv[is.finite(lv)]; xa <- cv[tx_clean]; xa <- xa[is.finite(xa)]
      mu_tx[[cl]] <- if (length(tv) > 0L && length(lv) > 0L && length(loc_rows) > length(focal_types)) 0.65 * stats::median(tv) + 0.35 * stats::median(lv)
        else if (length(tv) > 0L) stats::median(tv) else if (length(lv) > 0L) stats::median(lv) else if (length(xa) > 0L) stats::median(xa) else 0.0
      if (!is.finite(mu_tx[[cl]])) mu_tx[[cl]] <- 0.0
      if (!is.null(fixed_sd) && cl %in% names(fixed_sd) && is.finite(fixed_sd[[cl]])) {
        sd_tx[[cl]] <- max(fixed_sd[[cl]], scale_floors[[cl]])
      } else {
        xanc <- cv[tx_clean[is_anc_v[tx_clean]]]; xanc <- xanc[is.finite(xanc)]
        px <- if (length(lv) >= 2L) lv else if (length(xanc) >= 2L) xanc else xa
        mv <- if (length(px) >= 2L) stats::mad(px, center = mu_tx[[cl]], na.rm = TRUE) else NA_real_
        sv <- if (length(px) >= 2L) stats::sd(px, na.rm = TRUE) else NA_real_
        sd_tx[[cl]] <- max(c(if (is.finite(mv) && mv > 1e-6) mv else sv, scale_floors[[cl]]), na.rm = TRUE)
        if (!is.finite(sd_tx[[cl]]) || sd_tx[[cl]] <= 1e-6) sd_tx[[cl]] <- scale_floors[[cl]]
      }
    }
    tl_r <- focal_types[is_precise_edaph[focal_types]]; lith_ref <- unique(lith_vec[tl_r[nzchar(lith_vec[tl_r])]])
    if (length(lith_ref) == 0L) {
      lc <- lith_vec[loc_rows[is_precise_edaph[loc_rows]]]; lc <- lc[nzchar(lc)]
      if (length(lc) == 0L) { lc <- lith_vec[loc_rows]; lc <- lc[nzchar(lc)] }
      if (length(lc) == 0L) { lc <- lith_vec[focal_types]; lc <- lc[nzchar(lc)] }
      if (length(lc) > 0L) { tb <- table(lc); lith_ref <- names(tb[tb == max(tb)]) }
    }
    list(mu = mu_tx, sd = sd_tx, lith_ref = lith_ref, anchor_type_rows = focal_types)
  }

  eval_prof <- function(idx, pr, vars_s, wts_s) {
    if (length(vars_s) == 0L || is.null(pr)) return(list(z_rms = 0.0, sum_z2 = 0.0, p_eff = 0L))
    xv <- vapply(vars_s, function(cl) suppressWarnings(as.numeric(df[[cl]][idx])), numeric(1))
    z2 <- ((xv - pr$mu[vars_s]) / pmax(pr$sd[vars_s], scale_floors[vars_s], 1e-4))^2
    ok <- is.finite(z2) & is.finite(wts_s) & (wts_s > 0)
    if (!any(ok)) return(list(z_rms = 0.0, sum_z2 = 0.0, p_eff = 0L))
    list(z_rms = sqrt(sum(wts_s[ok] * z2[ok]) / sum(wts_s[ok])), sum_z2 = sum(z2[ok]), p_eff = as.integer(sum(ok)))
  }

  type_niche_status <- rep("Non_Type_Specimen", n); type_niche_audit_status <- rep("Non_Type", n)
  is_het_audit <- rep(FALSE, n); rf_source_hint <- rep(NA_character_, n)
  type_niche_status[is_quar_type] <- "Passive_Type_Uncertain_Coord"
  type_niche_audit_status[is_quar_type] <- "Passive_Type_Uncertain_Coord"
  rf_source_hint[is_quar_type] <- "Passive_Type_Uncertain_Coord"

  type_profiles <- list()
  audit_vars <- if (length(edaph_v) > 0L) edaph_v else env_v
  audit_wts <- if (length(edaph_v) > 0L) w_ed else w_en

  for (tx in active_taxa) {
    tx_rows <- which(grp_vec == tx); elig_t <- tx_rows[is_type_georef_ok[tx_rows]]
    homo_t <- elig_t[is_homo_all[elig_t]]; het_t <- elig_t[!is_homo_all[elig_t]]
    q_tx <- tx_rows[is_quar_type[tx_rows]]; excl_ss <- c(q_tx, het_t)
    clean_nt <- setdiff(tx_rows[is_anc_v[tx_rows]], c(elig_t, q_tx))
    topo_t <- if (length(homo_t) == 0L && length(clean_nt) >= 2L) {
      if (exists(".resolve_locus_classicus_topotype_rows", mode = "function")) .resolve_locus_classicus_topotype_rows(df, clean_nt, tx) else .fallback_topotype_rows_s05(df, clean_nt, tx)
    } else integer(0)
    if (length(homo_t) > 0L) {
      pr_ss <- build_profile(tx_rows, homo_t, excl_rows = excl_ss)
      type_niche_status[homo_t] <- "Confirmed_Type_Control"; type_niche_audit_status[homo_t] <- "Confirmed_Type_Control"
    } else if (length(topo_t) >= 2L) {
      pr_ss <- build_profile(tx_rows, topo_t, excl_rows = excl_ss); pr_ss$anchor_type_rows <- integer(0)
    } else if (length(het_t) > 0L) {
      pr_ss <- if (length(clean_nt) >= 2L) build_profile(tx_rows, integer(0), excl_rows = excl_ss) else build_profile(tx_rows, het_t, excl_rows = q_tx)
    } else pr_ss <- build_profile(tx_rows, integer(0), excl_rows = q_tx)

    if (length(het_t) > 0L) {
      for (hi in het_t) {
        ev_h <- eval_prof(hi, pr_ss, audit_vars, audit_wts); pe <- ev_h$p_eff
        if (pe >= 1L) {
          c_chi <- stats::qchisq(0.99, df = pe); c_rms <- sqrt(c_chi / pe)
          if ((is.finite(ev_h$z_rms) && ev_h$z_rms > c_rms) || (is.finite(ev_h$sum_z2) && ev_h$sum_z2 > c_chi)) {
            is_het_audit[hi] <- TRUE; type_niche_status[hi] <- "Type_Synonymy_Candidate"
            type_niche_audit_status[hi] <- "Type_Synonymy_Candidate"; rf_source_hint[hi] <- "Passive_Heterotypic_Type_Audit"
          } else {
            type_niche_status[hi] <- "Confirmed_Type_Control"; type_niche_audit_status[hi] <- "Concordant_Heterotypic_Type"
          }
        }
      }
    }
    type_profiles[[tx]] <- pr_ss
  }

  basionym_profiles <- list(); elig_type_all <- which(is_type_georef_ok)
  for (bn in unique(bas_all[elig_type_all])) {
    if (is.na(bn) || !nzchar(bn)) next
    b_rows <- elig_type_all[bas_all[elig_type_all] == bn]; tx_bn <- grp_vec[b_rows[1L]]; tx_r <- which(grp_vec == tx_bn)
    basionym_profiles[[bn]] <- if (all(is_homo_all[b_rows])) build_profile(tx_r, b_rows, excl_rows = tx_r[is_quar_type[tx_r] | (is_type_georef_ok[tx_r] & !is_homo_all[tx_r])])
      else build_profile(b_rows, b_rows, excl_rows = integer(0), fixed_sd = if (!is.null(type_profiles[[tx_bn]])) type_profiles[[tx_bn]]$sd else NULL)
  }
  for (tx in active_taxa) if (!any(grp_vec[elig_type_all] == tx) && !(tx %in% names(basionym_profiles))) basionym_profiles[[tx]] <- type_profiles[[tx]]
  bas_names <- names(basionym_profiles)

  type_edaphic_z <- numeric(n); nearest_type_niche <- character(n); matches_lith <- rep(TRUE, n)
  nearest_basionym <- rep(NA_character_, n); nearest_bas_z <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    g_sp <- grp_vec[i]; if (is.na(g_sp) || !(g_sp %in% active_taxa)) g_sp <- active_taxa[1L]
    z_ed <- vapply(active_taxa, function(tx) eval_prof(i, type_profiles[[tx]], edaph_v, w_ed)$z_rms, numeric(1))
    z_en <- vapply(active_taxa, function(tx) eval_prof(i, type_profiles[[tx]], env_v, w_en)$z_rms, numeric(1))
    gp <- match(g_sp, active_taxa); type_edaphic_z[i] <- round(z_ed[if (is.na(gp)) 1L else gp], 4)
    nearest_type_niche[i] <- active_taxa[which.min(0.60 * z_ed + 0.40 * z_en)][1L]
    if (length(bas_names) > 0L) {
      zb_ed <- vapply(bas_names, function(bn) eval_prof(i, basionym_profiles[[bn]], edaph_v, w_ed)$z_rms, numeric(1))
      zb_en <- vapply(bas_names, function(bn) eval_prof(i, basionym_profiles[[bn]], env_v, w_en)$z_rms, numeric(1))
      bb <- which.min(0.60 * zb_ed + 0.40 * zb_en)[1L]; nearest_basionym[i] <- bas_names[bb]; nearest_bas_z[i] <- round(zb_ed[bb], 4)
    }
    pr_g <- type_profiles[[g_sp]]
    matches_lith[i] <- if (!is.null(pr_g) && length(pr_g$lith_ref) > 0L && nzchar(lith_vec[i])) lith_vec[i] %in% pr_g$lith_ref else (type_edaphic_z[i] <= 2.5)
  }
  list(type_edaphic_z_score = type_edaphic_z, nearest_type_niche_taxon = nearest_type_niche, matches_type_lithology = as.logical(matches_lith),
       nearest_type_niche_basionym = nearest_basionym, nearest_type_basionym_edaphic_z = nearest_bas_z, is_heterotypic_type_audit = is_het_audit,
       type_niche_status = type_niche_status, type_niche_audit_status = type_niche_audit_status, rf_prediction_source_hint = rf_source_hint,
       type_niche_profiles = type_profiles, basionym_niche_profiles = basionym_profiles)
}
