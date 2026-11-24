#!/usr/bin/env Rscript
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Script: scripts/analysis/edaphic_collinearity_contingency_utils.R
# Project: Plantmatics / Packera Species Delimitation Pipeline (NCU)
# Description:
#   Point-mass censoring screening, Pearson/Spearman block-partitioned
#   collinearity & VIF pruning, EPV budgeting, CV-fold imputation, and
#   multi-class SGMC/SSURGO/NWI-hydroperiod substrate contingency testing.
# ==============================================================================

screen_edaphic_point_mass_columns <- function(X_mat, ref_mask = NULL, max_tie_frac = 0.25) {
  X_d <- as.data.frame(X_mat, check.names = FALSE); cols <- colnames(X_d)
  if (length(cols) == 0L) return(character(0))
  n_r <- nrow(X_d)
  ref_idx <- if (is.null(ref_mask)) seq_len(n_r) else if (is.logical(ref_mask)) which(!is.na(ref_mask) & ref_mask) else { iv <- suppressWarnings(as.integer(ref_mask)); iv[is.finite(iv) & iv >= 1L & iv <= n_r] }
  if (length(ref_idx) < 16L) return(cols)
  tie_f <- vapply(cols, function(cl) {
    v <- suppressWarnings(as.numeric(X_d[[cl]][ref_idx])); vf <- v[is.finite(v)]
    if (length(vf) == 0L) 1.0 else max(as.numeric(table(round(vf, 4L)))) / length(vf)
  }, numeric(1))
  names(tie_f) <- cols
  kept <- cols[is.finite(tie_f) & tie_f <= as.numeric(max_tie_frac)]
  if (length(kept) < 2L && length(cols) >= 2L) cols[order(tie_f, seq_along(cols))[seq_len(min(2L, length(cols)))]] else if (length(kept) == 0L && length(cols) == 1L) cols else kept
}

.pairwise_max_abs_cor <- function(X_num) {
  cp <- abs(suppressWarnings(stats::cor(X_num, use = "pairwise.complete.obs", method = "pearson")))
  cs <- abs(suppressWarnings(stats::cor(X_num, use = "pairwise.complete.obs", method = "spearman")))
  cp[!is.finite(cp)] <- 0.0; cs[!is.finite(cs)] <- 0.0
  cm <- pmax(cp, cs); diag(cm) <- 0.0; cm
}

compute_anchor_separation_scores <- function(X_df, anchor_groups = NULL) {
  X_d <- as.data.frame(X_df, check.names = FALSE); cols <- colnames(X_d)
  if (length(cols) == 0L) return(numeric(0))
  n_rows <- nrow(X_d); has_g <- !is.null(anchor_groups) && (length(anchor_groups) == n_rows)
  g_all <- if (has_g) trimws(as.character(anchor_groups)) else character(0)
  scores <- vapply(cols, function(cl) {
    v <- suppressWarnings(as.numeric(X_d[[cl]])); fin <- is.finite(v); ff <- if (n_rows > 0L) sum(fin) / n_rows else 0.0
    if (sum(fin) < 3L || !is.finite(stats::sd(v[fin])) || stats::sd(v[fin]) <= 1e-8) return(0.0)
    if (has_g) {
      ok <- fin & !is.na(g_all) & nzchar(g_all) & !(tolower(g_all) %in% c("unknown", "na", "none"))
      if (sum(ok) >= 4L) {
        vg <- names(table(g_all[ok])[table(g_all[ok]) >= 2L])
        if (length(vg) >= 2L) {
          mg <- ok & (g_all %in% vg); f_g <- factor(g_all[mg]); v_g <- v[mg]
          kw <- tryCatch(suppressWarnings(as.numeric(stats::kruskal.test(v_g ~ f_g)$statistic[1L])), error = function(e) NA_real_)
          fs <- tryCatch(suppressWarnings(as.numeric(summary(stats::aov(v_g ~ f_g))[[1L]][["F value"]][1L])), error = function(e) NA_real_)
          if (is.finite(kw) && kw >= 0.0) return(kw + (if (is.finite(fs) && fs >= 0.0) 1e-4 * fs else 0.0) + 1e-6 * ff)
          if (is.finite(fs) && fs >= 0.0) return(fs + 1e-6 * ff)
        }
      }
    }
    ff
  }, numeric(1))
  names(scores) <- cols; scores
}

screen_collinear_predictors <- function(X_env, cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = NULL, var_scores = NULL, quiet = FALSE) {
  X_df <- as.data.frame(X_env, check.names = FALSE); all_c <- colnames(X_df)
  mk_ret <- function(sel, sc) list(selected_variables = sel, selected_vars = sel, kept_vars = sel, pruned_variables = setdiff(all_c, sel), anchor_scores = sc)
  if (length(all_c) <= 1L) return(mk_ret(all_c, if (length(all_c) == 1L) stats::setNames(1.0, all_c) else numeric(0)))

  sds <- vapply(X_df, function(v) stats::sd(suppressWarnings(as.numeric(v)), na.rm = TRUE), numeric(1))
  keep_c <- names(sds[is.finite(sds) & sds > 1e-8])
  if (length(keep_c) == 0L) return(mk_ret(all_c[1L], stats::setNames(rep(0.0, length(all_c)), all_c)))
  if (length(keep_c) == 1L) return(mk_ret(keep_c, stats::setNames(1.0, keep_c)))

  sc_sub <- if (!is.null(var_scores) && all(keep_c %in% names(var_scores))) var_scores[keep_c] else compute_anchor_separation_scores(X_df[, keep_c, drop = FALSE], anchor_groups = anchor_groups)
  sc_sub[!is.finite(sc_sub)] <- 0.0; ranked <- keep_c[order(sc_sub[keep_c], decreasing = TRUE)]
  X_imp <- as.data.frame(lapply(ranked, function(cl) { v <- suppressWarnings(as.numeric(X_df[[cl]])); if (any(!is.finite(v))) { m <- stats::median(v[is.finite(v)], na.rm = TRUE); v[!is.finite(v)] <- if (is.finite(m)) m else 0.0 }; v }), stringsAsFactors = FALSE)
  colnames(X_imp) <- ranked

  sel_v <- NULL
  if (requireNamespace("spatialRF", quietly = TRUE) && ncol(X_imp) >= 3L) {
    ac <- tryCatch(spatialRF::auto_cor(x = X_imp, cor.threshold = cor.threshold, preference.order = ranked), error = function(e) tryCatch(spatialRF::auto_cor(x = X_imp, cor.threshold = cor.threshold), error = function(e2) NULL))
    X_ac <- if (!is.null(ac) && length(ac$selected.variables) >= 1L) X_imp[, ac$selected.variables, drop = FALSE] else X_imp
    if (ncol(X_ac) >= 2L) {
      av <- tryCatch(spatialRF::auto_vif(x = X_ac, vif.threshold = vif.threshold, preference.order = intersect(ranked, colnames(X_ac))), error = function(e) tryCatch(spatialRF::auto_vif(x = X_ac, vif.threshold = vif.threshold), error = function(e2) NULL))
      if (!is.null(av) && length(av$selected.variables) >= 1L) sel_v <- av$selected.variables
    } else if (ncol(X_ac) == 1L) sel_v <- colnames(X_ac)
  }

  # Enforce pmax(|r_pearson|, |r_spearman|) <= cor.threshold and VIF <= vif.threshold
  cmat <- .pairwise_max_abs_cor(X_imp)
  need_a <- is.null(sel_v) || length(sel_v) == 0L || (length(sel_v) >= 2L && is.finite(max(cmat[sel_v, sel_v, drop = FALSE], na.rm = TRUE)) && max(cmat[sel_v, sel_v, drop = FALSE], na.rm = TRUE) > (cor.threshold + 1e-6))
  if (need_a) {
    cand <- ranked
    while (length(cand) > 1L && max(cmat[cand, cand, drop = FALSE], na.rm = TRUE) > cor.threshold) {
      sc_m <- cmat[cand, cand, drop = FALSE]; ij <- which(sc_m == max(sc_m, na.rm = TRUE), arr.ind = TRUE)[1L, ]
      v1 <- cand[ij[1L]]; v2 <- cand[ij[2L]]; s1 <- sc_sub[[v1]]; s2 <- sc_sub[[v2]]
      cand <- setdiff(cand, if (is.finite(s1) && is.finite(s2) && abs(s1 - s2) > 1e-9) (if (s1 < s2) v1 else v2) else (if (mean(sc_m[v1, ]) >= mean(sc_m[v2, ])) v1 else v2))
    }
    while (length(cand) > 2L) {
      R_s <- suppressWarnings(stats::cor(X_imp[, cand, drop = FALSE], use = "pairwise.complete.obs")); R_s[!is.finite(R_s)] <- 0.0; diag(R_s) <- 1.0
      vifs <- tryCatch(diag(solve(R_s + diag(1e-6, length(cand)))), error = function(e) rep(1.0, length(cand))); names(vifs) <- cand
      mv <- max(vifs, na.rm = TRUE); if (!is.finite(mv) || mv <= vif.threshold) break
      ov <- names(vifs[vifs >= (mv - 1e-9)])
      cand <- setdiff(cand, if (length(ov) > 1L) ov[order(sc_sub[ov])][1L] else ov[1L])
    }
    sel_v <- cand
  }
  if (!isTRUE(quiet)) message(sprintf("Collinearity Screening (Dormann et al. 2013; |r| <= %.2f, VIF <= %.1f): retained %d / %d predictors.", cor.threshold, vif.threshold, length(sel_v), ncol(X_df)))
  mk_ret(sel_v, sc_sub)
}

screen_collinear_predictors_by_block <- function(X_env, clim_vars = NULL, edaph_vars = NULL, cor.threshold = 0.70, vif.threshold = 5.0, topo_vars = NULL, anchor_groups = NULL, quiet = FALSE) {
  X_df <- as.data.frame(X_env, check.names = FALSE); all_c <- colnames(X_df)
  if (length(all_c) == 0L) return(list(selected_variables = character(0), selected_vars = character(0), kept_vars = character(0), selected_edaphic = character(0), selected_climatic = character(0), selected_topo = character(0), edaphic_variables = character(0), climatic_variables = character(0), topo_variables = character(0), pruned_variables = character(0), anchor_scores = numeric(0)))

  clim_vars <- if (is.null(clim_vars)) grep("^(bio[0-9]*_|summer_moisture_|maritime_continentality_|delta_doy|doy|hopkins_|sin_doy|cos_doy)", all_c, value = TRUE, ignore.case = TRUE) else intersect(clim_vars, all_c)
  topo_vars <- if (is.null(topo_vars)) setdiff(grep("^(topo_|elevation($|_)|elev_)", all_c, value = TRUE, ignore.case = TRUE), unique(c(clim_vars, if (is.null(edaph_vars)) character(0) else edaph_vars))) else setdiff(intersect(topo_vars, all_c), clim_vars)
  edaph_vars <- if (is.null(edaph_vars)) setdiff(all_c, unique(c(clim_vars, topo_vars))) else setdiff(intersect(edaph_vars, all_c), unique(c(clim_vars, topo_vars)))

  sc_all <- compute_anchor_separation_scores(X_df, anchor_groups = anchor_groups)
  run_b <- function(vs) if (length(vs) > 0L) screen_collinear_predictors(X_df[, vs, drop = FALSE], cor.threshold = cor.threshold, vif.threshold = vif.threshold, anchor_groups = anchor_groups, var_scores = sc_all, quiet = quiet)$selected_variables else character(0)
  sel_c <- run_b(clim_vars); sel_e <- run_b(edaph_vars); sel_t <- run_b(topo_vars)
  sel_j <- unique(c(sel_e, sel_t, sel_c)); cross_pr <- character(0)

  if (length(sel_j) >= 2L && is.finite(vif.threshold)) {
    Xj <- do.call(cbind, lapply(sel_j, function(cl) { v <- suppressWarnings(as.numeric(X_df[[cl]])); f <- is.finite(v); m <- if (any(f)) stats::median(v[f], na.rm = TRUE) else 0.0; v[!f] <- if (is.finite(m)) m else 0.0; v }))
    colnames(Xj) <- sel_j; storage.mode(Xj) <- "double"
    sds_j <- apply(Xj, 2L, stats::sd); cand <- colnames(Xj)[is.finite(sds_j) & sds_j > 1e-8]; if (length(cand) == 0L) cand <- sel_j[1L]

    if (length(cand) >= 2L) {
      pb <- stats::setNames(rep(0.0, length(cand)), cand); pb[intersect(cand, sel_c)] <- 0.15; pb[intersect(cand, sel_e)] <- 0.18; pb[intersect(cand, c("wetland_seepage_index", "lithic_outcrop_index"))] <- 0.25
      jpref <- as.numeric(sc_all[cand]); jpref[!is.finite(jpref)] <- 0.0; names(jpref) <- cand; jpref <- jpref + pb[cand]
      has_mg <- FALSE; Xr <- NULL
      if (!is.null(anchor_groups) && length(anchor_groups) == nrow(Xj)) {
        gc <- trimws(as.character(anchor_groups)); okg <- !is.na(gc) & nzchar(gc) & !(tolower(gc) %in% c("unknown", "na", "none"))
        vg <- names(table(gc[okg])[table(gc[okg]) >= 2L])
        if (length(vg) >= 2L) {
          has_mg <- TRUE; ivg <- which(okg & (gc %in% vg)); Xr <- Xj[ivg, cand, drop = FALSE]; gs <- gc[ivg]
          for (gn in vg) { gm <- which(gs == gn); Xr[gm, ] <- if (length(gm) >= 2L) sweep(Xr[gm, , drop = FALSE], 2L, colMeans(Xr[gm, , drop = FALSE]), "-") else 0.0 }
        }
      }
      while (length(cand) >= 2L) {
        cmj <- .pairwise_max_abs_cor(Xj[, cand, drop = FALSE]); if (max(cmj, na.rm = TRUE) <= cor.threshold) break
        cmr <- if (isTRUE(has_mg) && !is.null(Xr)) .pairwise_max_abs_cor(Xr[, cand, drop = FALSE]) else NULL
        hp <- which(cmj > cor.threshold & upper.tri(cmj), arr.ind = TRUE); if (nrow(hp) == 0L) break
        dropped <- FALSE
        for (pi in order(cmj[hp], decreasing = TRUE)) {
          v1 <- cand[hp[pi, 1L]]; v2 <- cand[hp[pi, 2L]]
          if (!is.null(cmr) && (v1 %in% c(sel_c, sel_e)) && (v2 %in% c(sel_c, sel_e)) && is.finite(cmr[v1, v2]) && cmr[v1, v2] <= cor.threshold) next
          cand <- setdiff(cand, if (jpref[[v1]] < jpref[[v2]]) v1 else v2); dropped <- TRUE; break
        }
        if (!dropped) break
      }
      while (length(cand) >= 2L) {
        Rs <- suppressWarnings(stats::cor(Xj[, cand, drop = FALSE], use = "pairwise.complete.obs")); Rs[!is.finite(Rs)] <- 0.0; diag(Rs) <- 1.0
        vifs <- tryCatch(diag(solve(Rs + diag(1e-6, length(cand)))), error = function(e) rep(1.0, length(cand))); names(vifs) <- cand
        ov <- names(vifs[is.finite(vifs) & vifs > vif.threshold]); if (length(ov) == 0L) break
        prot <- character(0)
        if (isTRUE(has_mg) && !is.null(Xr)) {
          Rr <- suppressWarnings(stats::cor(Xr[, cand, drop = FALSE], use = "pairwise.complete.obs")); Rr[!is.finite(Rr)] <- 0.0; diag(Rr) <- 1.0
          rv <- tryCatch(diag(solve(Rr + diag(1e-6, length(cand)))), error = function(e) rep(Inf, length(cand))); names(rv) <- cand
          for (vv in ov) if ((vv %in% c(sel_c, sel_e)) && is.finite(rv[[vv]]) && rv[[vv]] <= vif.threshold) prot <- c(prot, vv)
        }
        drp <- setdiff(ov, prot); if (length(drp) == 0L) break
        top_o <- drp[vifs[drp] >= (max(vifs[drp], na.rm = TRUE) - 0.25)]
        cand <- setdiff(cand, top_o[order(jpref[top_o])][1L])
      }
    }
    cross_pr <- setdiff(sel_j, cand); sel_j <- cand; sel_e <- intersect(sel_e, sel_j); sel_t <- intersect(sel_t, sel_j); sel_c <- intersect(sel_c, sel_j)
  }
  if (!isTRUE(quiet)) message(sprintf("Partitioned 3-Block + Cross-Block Collinearity Pruning (VIF <= %.1f): retained %d Edaphic-Lithic, %d Topographic, and %d Climatic%s.", vif.threshold, length(sel_e), length(sel_t), length(sel_c), if (length(cross_pr) > 0L) sprintf(" [cross-block pruned: %s]", paste(cross_pr, collapse = ", ")) else ""))
  list(selected_variables = sel_j, selected_vars = sel_j, kept_vars = sel_j, selected_edaphic = sel_e, selected_climatic = sel_c, selected_topo = sel_t,
       edaphic_variables = sel_e, climatic_variables = sel_c, topo_variables = sel_t, pruned_variables = setdiff(all_c, sel_j), anchor_scores = sc_all)
}

enforce_epv_predictor_budget <- function(df, env_vars = NULL, selected_vars = NULL, n_anchors = nrow(df), anchor_groups = NULL, n_mems_planned = 0L, min_env = 4L, max_env = 12L, intraspecific_mode = FALSE, quiet = FALSE) {
  df <- as.data.frame(df, check.names = FALSE)
  ban <- if (exists("NON_ENM_RAW_GEOCHEM_HYDRO_COLS", inherits = TRUE)) get("NON_ENM_RAW_GEOCHEM_HYDRO_COLS", inherits = TRUE) else c("soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct", "ultramafic_calcareous_gradient", "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m", "nwpl_indicator_code", "nwpl_wetland_score")
  cand_in <- if (!is.null(selected_vars)) selected_vars else if (!is.null(env_vars)) env_vars else colnames(df)
  env_c <- setdiff(unique(intersect(as.character(cand_in), colnames(df))), ban)
  env_c <- env_c[!grepl("^polaris_", env_c, ignore.case = TRUE)]

  n_anc <- max(1L, as.integer(n_anchors)); tot_b <- max(1L, as.integer(floor(n_anc / 5L))); n_mp <- max(0L, as.integer(n_mems_planned))
  max_mems <- if (n_mp > 0L) min(n_mp, as.integer(if (n_anc < 30L) 2L else min(5L, max(2L, tot_b - as.integer(min_env))))) else 0L
  env_cap <- as.integer(min(as.integer(max_env), max(as.integer(min_env), tot_b - max_mems)))

  kept <- env_c
  if (length(env_c) > env_cap && env_cap >= 1L) {
    ss <- compute_anchor_separation_scores(df[, env_c, drop = FALSE], anchor_groups = anchor_groups); ss[!is.finite(ss)] <- 0.0
    pt <- intersect(env_c, c("wetland_seepage_index", "lithic_outcrop_index")); if (length(pt) > 0L) ss[pt] <- ss[pt] + 0.02
    kept <- names(sort(ss[env_c], decreasing = TRUE))[seq_len(env_cap)]
    if (env_cap >= 2L) {
      ac <- grep("^(bio[0-9]*_|summer_moisture_|maritime_continentality_)", env_c, value = TRUE, ignore.case = TRUE)
      ae <- setdiff(env_c, unique(c(ac, grep("^(topo_|elevation($|_)|elev_)", env_c, value = TRUE, ignore.case = TRUE))))
      pick_sw <- function(pv, prot = character(0)) {
        dp <- pv[order(ss[pv])]; bc <- c(c = sum(pv %in% ac), e = sum(pv %in% ae), t = sum(!(pv %in% c(ac, ae))))
        dc <- dp[vapply(dp, function(v) bc[[if (v %in% ac) "c" else if (v %in% ae) "e" else "t"]] >= 2L, logical(1))]
        if (length(dc) > 0L) dc[1L] else { np <- setdiff(dp, prot); if (length(np) > 0L) np[1L] else dp[1L] }
      }
      if (length(ac) > 0L && !any(kept %in% ac)) kept <- unique(c(setdiff(kept, pick_sw(kept)), names(sort(ss[ac], decreasing = TRUE))[1L]))
      if (length(ae) > 0L && !any(kept %in% ae)) kept <- unique(c(setdiff(kept, pick_sw(kept, ac)), names(sort(ss[ae], decreasing = TRUE))[1L]))
    }
    kept <- kept[order(ss[kept], decreasing = TRUE)]
  }
  if (!isTRUE(quiet) && length(kept) >= 1L && length(kept) < length(env_c)) message(sprintf("[EPV Predictor Budget] Capped predictors for N_anchors=%d (budget=%d, MEMs=%d): retained %d / %d.", n_anc, tot_b, max_mems, length(kept), length(env_c)))
  list(env_vars = kept, selected_vars = kept, max_mems = max_mems, total_budget = tot_b)
}

apply_fold_imputation <- function(X_te_raw_df, selected_vars, train_medians) {
  X_te <- as.data.frame(X_te_raw_df, check.names = FALSE); sv <- as.character(selected_vars); n_te <- nrow(X_te)
  if (length(sv) == 0L) return(X_te[, character(0), drop = FALSE])
  ol <- lapply(sv, function(vn) {
    mv <- if (!is.null(train_medians) && vn %in% names(train_medians)) suppressWarnings(as.numeric(train_medians[[vn]])) else NA_real_
    if (length(mv) != 1L || !is.finite(mv)) mv <- 0.0
    cv <- if (vn %in% colnames(X_te)) suppressWarnings(as.numeric(X_te[[vn]])) else rep(NA_real_, n_te)
    if (length(cv) != n_te) cv <- rep(mv, n_te); if (n_te > 0L) cv[!is.finite(cv)] <- mv; cv
  })
  names(ol) <- sv; out <- as.data.frame(ol, stringsAsFactors = FALSE, check.names = FALSE); colnames(out) <- sv; out
}

screen_predictors_within_cv_fold <- function(X_tr_raw_df, cand_vars, clim_vars = character(0), edaph_vars = character(0), topo_vars = character(0), y_tr = NULL, cor.threshold = 0.70, vif.threshold = 5.0, intraspecific_mode = FALSE) {
  X_tr <- as.data.frame(X_tr_raw_df, check.names = FALSE); cp <- intersect(as.character(cand_vars), colnames(X_tr))
  if (length(cp) == 0L) return(list(selected_vars = character(0), train_medians = numeric(0), X_tr_imp = X_tr[, character(0), drop = FALSE]))
  chk <- function(mk) vapply(cp, function(vn) { v <- suppressWarnings(as.numeric(X_tr[[vn]])); f <- is.finite(v); sum(f) >= mk && is.finite(stats::sd(v[f])) && stats::sd(v[f]) > 1e-8 }, logical(1))
  okc <- cp[chk(max(2L, as.integer(ceiling(0.50 * nrow(X_tr)))))]; if (length(okc) < 2L) okc <- cp[chk(2L)]; if (length(okc) == 0L) okc <- cp[1L]
  med_tr <- vapply(okc, function(vn) { v <- suppressWarnings(as.numeric(X_tr[[vn]])); f <- is.finite(v); m <- if (any(f)) stats::median(v[f], na.rm = TRUE) else 0.0; if (is.finite(m)) m else 0.0 }, numeric(1)); names(med_tr) <- okc
  Xi <- apply_fold_imputation(X_tr[, okc, drop = FALSE], okc, med_tr)
  ci <- intersect(as.character(clim_vars), colnames(Xi)); ei <- intersect(as.character(edaph_vars), colnames(Xi)); ti <- intersect(as.character(topo_vars), colnames(Xi))
  ub <- (length(ci) + length(ei) + length(ti)) > 0L
  scr <- screen_collinear_predictors_by_block(Xi, clim_vars = if (ub) ci else NULL, edaph_vars = if (ub) ei else NULL, cor.threshold = cor.threshold, vif.threshold = vif.threshold, topo_vars = if (ub) ti else NULL, anchor_groups = y_tr, quiet = TRUE)
  ecap <- max(1L, as.integer(floor(nrow(Xi) / 5L)))
  epv <- enforce_epv_predictor_budget(Xi, selected_vars = if (length(scr$selected_vars) > 0L) scr$selected_vars else colnames(Xi), n_anchors = nrow(Xi), anchor_groups = y_tr, n_mems_planned = 0L, min_env = min(4L, ecap), max_env = min(12L, ecap), intraspecific_mode = intraspecific_mode, quiet = TRUE)
  sv <- intersect(as.character(epv$selected_vars), colnames(Xi))
  if (length(sv) == 0L && ncol(Xi) >= 1L) sv <- colnames(Xi)[which.max(vapply(colnames(Xi), function(cn) { s <- stats::sd(Xi[[cn]], na.rm = TRUE); if (is.finite(s)) s else 0.0 }, numeric(1)))]
  if (nrow(Xi) >= 5L && length(sv) > ecap) sv <- sv[seq_len(ecap)]
  list(selected_vars = sv, train_medians = med_tr, X_tr_imp = as.data.frame(Xi[, sv, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE))
}

compute_fisher_2x2_safe <- function(a, b, c0, d) {
  a <- max(0L, as.integer(round(a))); b <- max(0L, as.integer(round(b))); c0 <- max(0L, as.integer(round(c0))); d <- max(0L, as.integer(round(d)))
  fr <- tryCatch(stats::fisher.test(matrix(c(a, b, c0, d), nrow = 2L, byrow = TRUE)), error = function(e) list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = 1.0))
  ha_or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c0 + 0.5)); ha_se <- sqrt(1.0 / (a + 0.5) + 1.0 / (b + 0.5) + 1.0 / (c0 + 0.5) + 1.0 / (d + 0.5))
  ro <- suppressWarnings(as.numeric(fr$estimate)); rl <- if (length(fr$conf.int) >= 1L) suppressWarnings(as.numeric(fr$conf.int[1L])) else NA_real_
  rh <- if (length(fr$conf.int) >= 2L) suppressWarnings(as.numeric(fr$conf.int[2L])) else NA_real_; rp <- suppressWarnings(as.numeric(fr$p.value))
  list(odds_ratio = if (is.finite(ro) && ro > 0) ro else ha_or, ci_lower = if (is.finite(rl) && rl >= 0) rl else exp(log(ha_or) - 1.96 * ha_se),
       ci_upper = if (is.finite(rh)) rh else exp(log(ha_or) + 1.96 * ha_se), p_value = if (is.finite(rp)) max(0.0, min(1.0, rp)) else 1.0, fisher_obj = fr)
}

.parse_bool_flag <- function(x, n) {
  if (is.null(x)) return(rep(FALSE, n))
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & is.finite(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "1", "t", "yes")
}

evaluate_multiclass_substrate_contingency <- function(df, active_groups = NULL, group_col = "species_standardized") {
  d <- as.data.frame(df, stringsAsFactors = FALSE); n_d <- nrow(d)
  if (!(group_col %in% names(d))) stop(sprintf("Column '%s' not found in dataframe for substrate contingency analysis.", group_col))
  grp_all <- as.character(d[[group_col]])
  if (is.null(active_groups) || length(active_groups) == 0L) { active_groups <- sort(unique(stats::na.omit(grp_all))); active_groups <- active_groups[nzchar(active_groups) & active_groups != "Unknown"] }
  in_m <- !is.na(grp_all) & (grp_all %in% active_groups)
  is_t2 <- if ("is_type_specimen" %in% names(d)) .parse_bool_flag(d$is_type_specimen, n_d) else rep(FALSE, n_d)
  if ("is_confident_label_noise" %in% names(d)) in_m <- in_m & (!.parse_bool_flag(d$is_confident_label_noise, n_d) | is_t2)

  # Exclude open/hybrid concepts and quarantined centroid/rooftop/fine-edaphic-ineligible coordinates unconditionally
  oh <- if (exists(".identify_open_or_hybrid_concept_mask", mode = "function")) as.logical(.identify_open_or_hybrid_concept_mask(d)) else rep(FALSE, n_d)
  oh[is.na(oh)] <- FALSE
  if ("is_open_or_hybrid_concept" %in% names(d)) oh <- oh | .parse_bool_flag(d$is_open_or_hybrid_concept, n_d)
  cq <- rep(FALSE, n_d)
  if (exists("detect_r_centroid_or_inst_quarantine", mode = "function")) {
    cr <- detect_r_centroid_or_inst_quarantine(d, is_type_vec = rep(FALSE, n_d))
    if (is.list(cr) && length(cr$is_quarantined) == n_d) cq <- cq | (!is.na(cr$is_quarantined) & as.logical(cr$is_quarantined))
  }
  if ("is_centroid_or_inst_coord" %in% names(d)) cq <- cq | .parse_bool_flag(d$is_centroid_or_inst_coord, n_d)
  if (exists("is_fine_edaphic_type_or_anchor_eligible", mode = "function")) { fe <- tryCatch(is_fine_edaphic_type_or_anchor_eligible(d, 10000.0, 50000.0), error = function(e) NULL); if (is.list(fe) && length(fe$edaphic_ok) == n_d) cq <- cq | !as.logical(fe$edaphic_ok) }
  if ("is_fine_edaphic_type_or_anchor_eligible" %in% names(d)) { fe <- d$is_fine_edaphic_type_or_anchor_eligible; cq <- cq | (!is.na(fe) & !.parse_bool_flag(fe, n_d)) }
  in_m <- in_m & !oh & !cq

  eq <- if (exists(".is_step05_extralimital_non_type", mode = "function")) { ev <- .is_step05_extralimital_non_type(d, group_vec = grp_all, clean_base_mask = in_m, is_type_vec = is_t2); if (length(ev) == n_d) (!is.na(ev) & as.logical(ev)) else rep(FALSE, n_d) } else rep(FALSE, n_d)
  if ("determiner_tier" %in% names(d)) {
    dt <- as.character(d$determiner_tier); dt[is.na(dt)] <- "Tier_3_Bronze"
    pr <- if ("is_promoted_bronze" %in% names(d)) .parse_bool_flag(d$is_promoted_bronze, n_d) else rep(FALSE, n_d)
    am <- in_m & (grepl("Tier_1|Gold|Tier_2|Silver", dt, ignore.case = TRUE) | pr | is_t2)
    in_m <- if (sum(am & !eq) >= 6L) am & !eq else if (sum(am) >= 6L) am else if (sum(in_m & !eq) >= 6L) in_m & !eq else in_m
  } else if (sum(in_m & !eq) >= 6L) in_m <- in_m & !eq

  post_sel <- isTRUE((group_col != "species_standardized") && exists("is_post_selection_ecotype_groups", mode = "function") && isTRUE(is_post_selection_ecotype_groups(active_groups, TRUE, d)))
  tdf <- d[in_m, , drop = FALSE]; n_tot <- nrow(tdf)
  if (n_tot == 0L) return(list(summary_table = data.frame(), multiclass_details = data.frame(), omnibus_test = NULL))

  grp_v <- as.character(tdf[[group_col]])
  is_rk <- if ("is_rock_outcrop" %in% names(tdf)) .parse_bool_flag(tdf$is_rock_outcrop, n_tot) else rep(FALSE, n_tot)
  is_az <- if ("is_azonal_lithology" %in% names(tdf)) .parse_bool_flag(tdf$is_azonal_lithology, n_tot) else rep(FALSE, n_tot)
  if ("ssurgo_rock_outcrop_pct" %in% names(tdf)) { rp <- suppressWarnings(as.numeric(tdf$ssurgo_rock_outcrop_pct)); is_rk <- is_rk | (is.finite(rp) & rp >= 15.0) }
  is_oa <- is_rk | is_az

  sgmc_v <- if ("sgmc_lithology_class" %in% names(tdf)) as.character(tdf$sgmc_lithology_class) else rep(NA_character_, n_tot)
  sgmc_v[is.na(sgmc_v) | !nzchar(trimws(sgmc_v))] <- NA_character_; p_sgmc <- unique(c(SGMC_LITHOLOGY_CLASSES, stats::na.omit(unique(sgmc_v))))
  out_v <- if ("outcrop_class" %in% names(tdf)) as.character(tdf$outcrop_class) else rep("Non-outcrop matrix", n_tot)
  out_v[is.na(out_v) | !nzchar(trimws(out_v))] <- "Non-outcrop matrix"; p_out <- unique(c(SSURGO_OUTCROP_CLASSES, unique(out_v)))
  sub_v <- if ("taxsuborder" %in% names(tdf)) trimws(as.character(tdf$taxsuborder)) else rep(NA_character_, n_tot)
  sub_v[is.na(sub_v) | !nzchar(sub_v) | tolower(sub_v) == "na"] <- NA_character_; p_sub <- sort(unique(stats::na.omit(sub_v)))

  # 4) USFWS NWI Hydroperiod Regime (Saturated_Seep_Fen_Pocosin, Seasonally_Flooded_Bottomland, Upland_Xeric_Mesic)
  nw_w <- if ("is_nwi_wetland" %in% names(tdf)) suppressWarnings(as.numeric(tdf$is_nwi_wetland)) else rep(NA_real_, n_tot)
  nw_d <- if ("nwi_wetland_dist_m" %in% names(tdf)) suppressWarnings(as.numeric(tdf$nwi_wetland_dist_m)) else rep(NA_real_, n_tot)
  nw_s <- if ("nwi_wetland_system" %in% names(tdf)) trimws(as.character(tdf$nwi_wetland_system)) else rep(NA_character_, n_tot)
  nw_s[is.na(nw_s) | !nzchar(nw_s) | tolower(nw_s) %in% c("na", "null")] <- NA_character_
  gg_v <- if ("taxgreatgroup" %in% names(tdf)) trimws(as.character(tdf$taxgreatgroup)) else rep("", n_tot)
  sg_txt <- paste(ifelse(is.na(sub_v), "", sub_v), ifelse(is.na(gg_v), "", gg_v))
  is_sat_sub <- grepl("aquod|saprist|hemist|fibrist", sg_txt, ignore.case = TRUE)
  has_wet <- any(is.finite(nw_w)) || any(is.finite(nw_d)) || any(!is.na(nw_s)) || any(is_sat_sub)
  nw_sc <- ifelse(is.na(nw_s), "", nw_s)
  is_sat <- is_sat_sub | grepl("^(PSS|PFO|PEM)[0-9/]*((EM|FO|SS)[0-9/]*)?[BEFGH][a-z]*$", nw_sc, ignore.case = TRUE) | grepl("saturated|semiperm|semi-perm|fen|pocosin|bog|seep", nw_sc, ignore.case = TRUE)
  is_wet_a <- (!is.na(nw_w) & nw_w == 1) | (is.finite(nw_d) & nw_d <= 30.0) | grepl("^(PSS|PFO|PEM)", nw_sc, ignore.case = TRUE)
  hydro_v <- ifelse(is_sat, "Saturated_Seep_Fen_Pocosin", ifelse(is_wet_a | (has_wet & grepl("aquult|aquept|aquoll|aqualf|aquent|fluvent", sg_txt, ignore.case = TRUE)), "Seasonally_Flooded_Bottomland", "Upland_Xeric_Mesic"))
  p_hydro <- c("Saturated_Seep_Fen_Pocosin", "Seasonally_Flooded_Bottomland", "Upland_Xeric_Mesic")

  n_bg <- max(120L, 2L * n_tot); single_g <- length(active_groups) < 2L
  bg_map <- list(
    SGMC_Lithology = c(Ultramafic_Serpentine = 0.02, Mafic_Igneous_Metamorphic = 0.05, Carbonate_Calcareous = 0.18, Felsic_Granitic_Gneiss = 0.11, Siliciclastic_Consolidated = 0.36, Unconsolidated_Sedimentary = 0.28),
    SSURGO_Outcrop = c("Granitic flatrock" = 0.03, "Limestone glade" = 0.03, "Sandstone/Siliceous" = 0.02, Serpentine = 0.005, "Diabase/Mafic glade" = 0.015, "Non-outcrop matrix" = 0.90),
    NWI_Hydroperiod_Regime = c(Saturated_Seep_Fen_Pocosin = 0.03, Seasonally_Flooded_Bottomland = 0.08, Upland_Xeric_Mesic = 0.89)
  )

  detail_rows <- list()
  eval_dom <- function(tx, in_tx, vv, clss, dom, def_bg, req_a = FALSE, s_req_a = FALSE) {
    vi <- in_tx & !is.na(vv); vo <- (!in_tx) & !is.na(vv)
    if (!any(vi) || length(clss) == 0L) return(invisible(NULL))
    for (cls in clss) {
      a <- sum(vi & vv == cls); b <- sum(vi & vv != cls)
      p_bg <- if (dom %in% names(bg_map) && cls %in% names(bg_map[[dom]])) bg_map[[dom]][[cls]] else if (dom == "SSURGO_Drainage_Suborder" && grepl("^Aqu|Sapr|Hist", cls, ignore.case = TRUE)) 0.12 else def_bg
      c0 <- if (!single_g && any(vo)) sum(vo & vv == cls) else max(1L, as.integer(round(n_bg * p_bg)))
      d0 <- if (!single_g && any(vo)) sum(vo & vv != cls) else max(1L, n_bg - c0)
      if ((a + c0) > 0L && (!req_a || a > 0L) && (!s_req_a || !single_g || a > 0L)) {
        ft <- compute_fisher_2x2_safe(a, b, c0, d0)
        detail_rows[[length(detail_rows) + 1L]] <<- data.frame(Taxon = tx, Substrate_Domain = dom, Substrate_Class = cls, In_Class_Count = a, Out_Class_Count = b,
                                                               Comparator_In_Count = c0, Comparator_Out_Count = d0, Class_Frequency_Pct = round(100.0 * a / max(1L, a + b), 2),
                                                               Odds_Ratio = round(ft$odds_ratio, 3), CI_Lower = round(ft$ci_lower, 3), CI_Upper = round(ft$ci_upper, 3), P_Value = ft$p_value, stringsAsFactors = FALSE)
      }
    }
  }
  for (tx in active_groups) {
    in_tx <- grp_v == tx
    eval_dom(tx, in_tx, sgmc_v, p_sgmc, "SGMC_Lithology", 0.05)
    eval_dom(tx, in_tx, out_v, p_out, "SSURGO_Outcrop", 0.02, s_req_a = TRUE)
    eval_dom(tx, in_tx, sub_v, p_sub, "SSURGO_Drainage_Suborder", 0.25, req_a = TRUE)
    if (has_wet) eval_dom(tx, in_tx, hydro_v, p_hydro, "NWI_Hydroperiod_Regime", 0.05, s_req_a = TRUE)
  }

  mc_df <- if (length(detail_rows) > 0L) do.call(rbind, detail_rows) else data.frame(Taxon = character(0), Substrate_Domain = character(0), Substrate_Class = character(0), In_Class_Count = integer(0), Out_Class_Count = integer(0), Comparator_In_Count = integer(0), Comparator_Out_Count = integer(0), Class_Frequency_Pct = numeric(0), Odds_Ratio = numeric(0), CI_Lower = numeric(0), CI_Upper = numeric(0), P_Value = numeric(0), P_Value_BH = numeric(0), stringsAsFactors = FALSE)
  if (nrow(mc_df) > 0L) {
    if (isTRUE(post_sel)) {
      mc_df$P_Value_BH <- rep(NA_real_, nrow(mc_df)); no <- mc_df$Substrate_Domain != "SSURGO_Outcrop"
      if (any(no)) mc_df$P_Value_BH[no] <- stats::p.adjust(mc_df$P_Value[no], method = "BH")
      mc_df$P_Value[!no] <- NA_real_; mc_df$P_Value_BH[!no] <- NA_real_
    } else mc_df$P_Value_BH <- stats::p.adjust(mc_df$P_Value, method = "BH")
  }

  smean <- function(sd, cn) { if (!(cn %in% names(sd))) return(NA_real_); v <- suppressWarnings(as.numeric(sd[[cn]])); if (all(is.na(v))) NA_real_ else round(mean(v, na.rm = TRUE), 2) }
  sp_rows <- list(); raw_bp <- stats::setNames(numeric(length(active_groups)), active_groups)
  for (tx in active_groups) {
    in_tx <- grp_v == tx; sub_t <- tdf[in_tx, , drop = FALSE]; so <- sum(in_tx & is_oa); sm <- sum(in_tx & !is_oa)
    oo <- if (!single_g && any(!in_tx)) sum((!in_tx) & is_oa) else max(1L, as.integer(round(n_bg * 0.08)))
    om <- if (!single_g && any(!in_tx)) sum((!in_tx) & !is_oa) else max(1L, n_bg - oo)
    ft_b <- compute_fisher_2x2_safe(so, sm, oo, om); raw_bp[[tx]] <- ft_b$p_value

    tx_mc <- mc_df[mc_df$Taxon == tx & mc_df$In_Class_Count > 0L, , drop = FALSE]; tx_sg <- tx_mc[tx_mc$Substrate_Domain == "SGMC_Lithology", , drop = FALSE]
    dom_l <- "Unclassified"; top_l <- "Unclassified"
    if (nrow(tx_sg) > 0L) {
      dom_l <- as.character(tx_sg$Substrate_Class[order(tx_sg$In_Class_Count, tx_sg$Odds_Ratio, decreasing = TRUE)[1L]])
      sig_s <- tx_sg[!is.na(tx_sg$P_Value_BH) & tx_sg$P_Value_BH < 0.05 & tx_sg$Odds_Ratio > 1.0, , drop = FALSE]
      if (nrow(sig_s) > 0L) {
        top_l <- as.character(sig_s$Substrate_Class[order(sig_s$Odds_Ratio, sig_s$In_Class_Count, decreasing = TRUE)[1L]]); dom_l <- top_l
      } else {
        ol <- order(tx_sg$Odds_Ratio, tx_sg$In_Class_Count, decreasing = TRUE); top_l <- if (tx_sg$Odds_Ratio[ol[1L]] > 1.0) as.character(tx_sg$Substrate_Class[ol[1L]]) else dom_l
      }
    } else if (any(in_tx)) { tb <- sort(table(out_v[in_tx]), decreasing = TRUE); if (length(tb) > 0L) { dom_l <- names(tb)[1L]; top_l <- dom_l } }

    dom_s <- dom_l; top_or <- round(ft_b$odds_ratio, 3); best_q <- NA_real_
    if (nrow(tx_mc) > 0L) {
      sig_m <- tx_mc[!is.na(tx_mc$P_Value_BH) & tx_mc$P_Value_BH < 0.05 & tx_mc$Odds_Ratio > 1.0, , drop = FALSE]
      pm <- if (nrow(sig_m) > 0L) sig_m else { nm <- tx_mc[tx_mc$Substrate_Class != "Non-outcrop matrix" & tx_mc$Odds_Ratio > 1.0, , drop = FALSE]; if (nrow(nm) > 0L) nm else tx_mc }
      bi <- order(pm$Odds_Ratio, pm$In_Class_Count, decreasing = TRUE)[1L]
      dom_s <- as.character(pm$Substrate_Class[bi]); top_or <- round(as.numeric(pm$Odds_Ratio[bi]), 3); best_q <- as.numeric(pm$P_Value_BH[bi])
    }
    pv_o <- if (isTRUE(post_sel)) NA_real_ else signif(ft_b$p_value, 5)
    sp_rows[[tx]] <- data.frame(
      Taxon = tx, Outcrop_Count = so, Matrix_Count = sm, Total_Vouchers = so + sm, Outcrop_Affinity_Pct = round(100.0 * so / max(1L, so + sm), 2),
      Odds_Ratio = round(ft_b$odds_ratio, 3), CI_Lower = round(ft_b$ci_lower, 3), CI_Upper = round(ft_b$ci_upper, 3),
      P_Value = pv_o, P_Value_BH = pv_o, Outcrop_Fisher_P = pv_o, Top_Substrate_P_Value_BH = signif(ifelse(is.finite(best_q), best_q, 1.0), 5),
      Dominant_Lithology_Class = dom_l, Dominant_Substrate_Class = dom_s, Top_Enriched_Lithology = top_l, Top_Substrate_Odds_Ratio = top_or,
      Mean_Rock_Outcrop_Pct = smean(sub_t, "ssurgo_rock_outcrop_pct"), Mean_Bedrock_Depth_cm = smean(sub_t, "fine_bedrock_depth_cm"),
      Mean_Fine_Sand_Pct = smean(sub_t, "fine_soil_sand"), Mean_Fine_pH = smean(sub_t, "fine_soil_ph"), Mean_AWS_0_50cm = smean(sub_t, "ssurgo_aws_0_50cm"),
      Mean_Topo_TPI_30m = smean(sub_t, "topo_tpi_30m"), Mean_Topo_TWI_30m = smean(sub_t, "topo_twi_30m"), stringsAsFactors = FALSE
    )
  }

  sp_df <- do.call(rbind, sp_rows); rownames(sp_df) <- NULL
  if (!isTRUE(post_sel)) { bbh <- stats::p.adjust(raw_bp, method = "BH"); sp_df$P_Value_BH <- vapply(sp_df$Taxon, function(tx) signif(if (is.finite(bbh[[tx]])) bbh[[tx]] else 1.0, 5), numeric(1)) }

  omni_f <- NULL; omni_p <- 1.0; p_bin_o <- if (isTRUE(post_sel)) NA_real_ else as.numeric(sp_df$Outcrop_Fisher_P[1L])
  if (!single_g) {
    ttab <- function(lv) {
      ok <- !is.na(lv) & nzchar(trimws(as.character(lv))); if (sum(ok) < 2L) return(NULL)
      tm <- table(factor(grp_v[ok], levels = active_groups), as.character(lv[ok])); tm <- tm[rowSums(tm) > 0L, colSums(tm) > 0L, drop = FALSE]
      if (nrow(tm) < 2L || ncol(tm) < 2L) return(NULL)
      tryCatch(if (nrow(tm) == 2L && ncol(tm) == 2L) stats::fisher.test(tm) else stats::fisher.test(tm, simulate.p.value = TRUE, B = 5000),
               error = function(e) tryCatch(stats::chisq.test(tm, simulate.p.value = TRUE, B = 5000), error = function(e2) NULL))
    }
    if (isTRUE(post_sel)) {
      cs <- Filter(function(ft) !is.null(ft) && is.finite(ft$p.value), list(ttab(sgmc_v), ttab(sub_v), if (has_wet) ttab(hydro_v) else NULL))
      if (length(cs) > 0L) { pc <- vapply(cs, function(ft) as.numeric(ft$p.value), numeric(1)); omni_f <- cs[[which.min(pc)]]; omni_p <- min(pc) }
      p_bin_o <- NA_real_
    } else {
      fb <- ttab(ifelse(is_oa, "Outcrop_Azonal", "Matrix")); p_bin_o <- if (!is.null(fb) && is.finite(fb$p.value)) as.numeric(fb$p.value) else 1.0
      fm <- ttab(if (any(!is.na(sgmc_v))) ifelse(!is.na(sgmc_v), sgmc_v, out_v) else out_v); pm <- if (!is.null(fm) && is.finite(fm$p.value)) as.numeric(fm$p.value) else 1.0
      omni_p <- min(p_bin_o, pm, na.rm = TRUE); omni_f <- if (!is.null(fm) && pm <= p_bin_o) fm else fb
    }
  } else omni_p <- if (isTRUE(post_sel)) NA_real_ else as.numeric(sp_df$P_Value[1L])

  ti <- order(sp_df$Top_Substrate_Odds_Ratio, sp_df$Odds_Ratio, decreasing = TRUE)[1L]; bi <- order(sp_df$Odds_Ratio, sp_df$Top_Substrate_Odds_Ratio, decreasing = TRUE)[1L]
  to_n <- sum(is_oa); tm_n <- sum(!is_oa); o_or <- max(sp_df$Odds_Ratio[bi], sp_df$Top_Substrate_Odds_Ratio[ti], na.rm = TRUE)
  qp <- c(sp_df$P_Value_BH, omni_p); qp <- qp[is.finite(qp)]
  sum_r <- data.frame(
    Taxon = "Omnibus_Substrate_Contingency_Across_Active_Groups", Outcrop_Count = to_n, Matrix_Count = tm_n, Total_Vouchers = to_n + tm_n,
    Outcrop_Affinity_Pct = round(100.0 * to_n / max(1L, to_n + tm_n), 2), Odds_Ratio = round(o_or, 3), CI_Lower = round(sp_df$CI_Lower[bi], 3),
    CI_Upper = round(max(sp_df$CI_Upper[bi], o_or, na.rm = TRUE), 3), P_Value = if (is.finite(omni_p)) signif(omni_p, 5) else NA_real_,
    P_Value_BH = if (length(qp) > 0L) signif(min(qp), 5) else NA_real_, Outcrop_Fisher_P = if (isTRUE(post_sel)) NA_real_ else signif(if (is.finite(p_bin_o)) p_bin_o else omni_p, 5),
    Top_Substrate_P_Value_BH = signif(if (any(is.finite(sp_df$Top_Substrate_P_Value_BH))) as.numeric(sp_df$Top_Substrate_P_Value_BH[ti]) else 1.0, 5),
    Dominant_Lithology_Class = sp_df$Dominant_Lithology_Class[ti], Dominant_Substrate_Class = sp_df$Dominant_Substrate_Class[ti],
    Top_Enriched_Lithology = sp_df$Top_Enriched_Lithology[ti], Top_Substrate_Odds_Ratio = sp_df$Top_Substrate_Odds_Ratio[ti],
    Mean_Rock_Outcrop_Pct = smean(tdf, "ssurgo_rock_outcrop_pct"), Mean_Bedrock_Depth_cm = smean(tdf, "fine_bedrock_depth_cm"),
    Mean_Fine_Sand_Pct = smean(tdf, "fine_soil_sand"), Mean_Fine_pH = smean(tdf, "fine_soil_ph"), Mean_AWS_0_50cm = smean(tdf, "ssurgo_aws_0_50cm"),
    Mean_Topo_TPI_30m = smean(tdf, "topo_tpi_30m"), Mean_Topo_TWI_30m = smean(tdf, "topo_twi_30m"), stringsAsFactors = FALSE
  )
  rep_df <- rbind(sum_r, sp_df); rownames(rep_df) <- NULL
  if (is.null(omni_f)) omni_f <- list(p.value = omni_p, estimate = o_or)
  list(summary_table = rep_df, multiclass_details = mc_df, omnibus_test = omni_f)
}

evaluate_outcrop_contingency <- evaluate_multiclass_substrate_contingency
