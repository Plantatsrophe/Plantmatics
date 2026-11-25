if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/dissect_scalars_permanova.R
# Proximity-bounded sheet consolidation, scalar variance/collinearity filters, and trace PERMANOVA
# ==============================================================================

consolidate_fragmented_sheet_individuals <- function(efa_df, max_dist_mm = 35.0, max_rel_span = 0.12) {
  if (is.null(efa_df) || !is.data.frame(efa_df) || nrow(efa_df) < 2L) return(efa_df)
  if (!all(c("catalogNumber", "plant_individual_id") %in% names(efa_df))) return(efa_df)
  n_r <- nrow(efa_df)
  cat_s03 <- as.character(efa_df$catalogNumber)
  pid_s03 <- ifelse(is.na(efa_df$plant_individual_id), 0L, suppressWarnings(as.integer(efa_df$plant_individual_id)))
  pid_s03[is.na(pid_s03)] <- 0L
  lz_s03 <- if ("leaf_zone" %in% names(efa_df)) tolower(trimws(as.character(efa_df$leaf_zone))) else rep("basal", n_r)
  lz_s03[is.na(lz_s03) | !(lz_s03 %in% c("basal", "cauline"))] <- "basal"
  cy_s03 <- if ("cy_rel_plant_bbox" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$cy_rel_plant_bbox)) else rep(NA_real_, n_r)
  dc_s03 <- if ("dist_from_crown_norm" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$dist_from_crown_norm)) else rep(NA_real_, n_r)
  cx_px <- if ("cx_px" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$cx_px)) else if ("leaf_cx_px" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$leaf_cx_px)) else rep(NA_real_, n_r)
  if (("cx_px" %in% names(efa_df)) && ("leaf_cx_px" %in% names(efa_df))) {
    alt_cx <- suppressWarnings(as.numeric(efa_df$leaf_cx_px))
    cx_px[!is.finite(cx_px)] <- alt_cx[!is.finite(cx_px)]
  }
  cy_px <- if ("cy_px" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$cy_px)) else if ("leaf_cy_px" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$leaf_cy_px)) else rep(NA_real_, n_r)
  if (("cy_px" %in% names(efa_df)) && ("leaf_cy_px" %in% names(efa_df))) {
    alt_cy <- suppressWarnings(as.numeric(efa_df$leaf_cy_px))
    cy_px[!is.finite(cy_px)] <- alt_cy[!is.finite(cy_px)]
  }
  if (all(c("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2") %in% names(efa_df))) {
    bx1 <- suppressWarnings(as.numeric(efa_df$bbox_x1))
    by1 <- suppressWarnings(as.numeric(efa_df$bbox_y1))
    bx2 <- suppressWarnings(as.numeric(efa_df$bbox_x2))
    by2 <- suppressWarnings(as.numeric(efa_df$bbox_y2))
    bb_ok <- is.finite(bx1) & is.finite(by1) & is.finite(bx2) & is.finite(by2)
    need_cx <- !is.finite(cx_px) & bb_ok
    need_cy <- !is.finite(cy_px) & bb_ok
    if (any(need_cx)) cx_px[need_cx] <- 0.5 * (bx1[need_cx] + bx2[need_cx])
    if (any(need_cy)) cy_px[need_cy] <- 0.5 * (by1[need_cy] + by2[need_cy])
  }
  ppm <- if ("pixels_per_mm" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$pixels_per_mm)) else rep(NA_real_, n_r)
  ss_s03 <- if ("scale_source" %in% names(efa_df)) tolower(trimws(as.character(efa_df$scale_source))) else rep("", n_r)
  ss_s03[is.na(ss_s03)] <- ""
  cal_flag <- if ("is_calibrated_mm" %in% names(efa_df)) {
    v_c <- efa_df$is_calibrated_mm
    if (is.character(v_c) || is.factor(v_c)) (!is.na(v_c) & (toupper(trimws(as.character(v_c))) %in% c("1", "TRUE", "T", "YES"))) else (!is.na(v_c) & as.logical(v_c))
  } else {
    rep(TRUE, n_r)
  }
  is_cal <- cal_flag & !(ss_s03 %in% c("mp_estimated", "predicted_from_megapixels", "tick_lattice_unboxed", "uncalibrated"))
  cx_sh <- if ("cx_rel_sheet" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$cx_rel_sheet)) else rep(NA_real_, n_r)
  cy_sh <- if ("cy_rel_sheet" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$cy_rel_sheet)) else rep(NA_real_, n_r)

  by_cat_s03 <- split(seq_len(n_r), cat_s03)
  reassigned_s03 <- 0L
  rezoned_strays_s03 <- 0L
  for (c_idx_vec in by_cat_s03) {
    u_pids <- sort(unique(pid_s03[c_idx_vec]))
    if (length(u_pids) < 2L) next
    bas_counts <- vapply(u_pids, function(p) sum(pid_s03[c_idx_vec] == p & lz_s03[c_idx_vec] == "basal"), integer(1))
    tot_counts <- vapply(u_pids, function(p) sum(pid_s03[c_idx_vec] == p), integer(1))
    hub_mask <- bas_counts >= 2L
    kept_pids <- if (any(hub_mask)) u_pids[hub_mask] else u_pids[which.max(bas_counts * 1000L + tot_counts)]
    primary_pid <- kept_pids[which.max(bas_counts[match(kept_pids, u_pids)] * 1000L + tot_counts[match(kept_pids, u_pids)])]
    prim_rows <- c_idx_vec[pid_s03[c_idx_vec] == primary_pid]
    prim_bas_rows <- c_idx_vec[pid_s03[c_idx_vec] == primary_pid & lz_s03[c_idx_vec] == "basal"]
    prim_cy_fin <- cy_s03[prim_bas_rows][is.finite(cy_s03[prim_bas_rows])]
    min_prim_bas_cy <- if (length(prim_cy_fin) > 0L) min(prim_cy_fin) else NA_real_
    med_prim_bas_cy <- if (length(prim_cy_fin) > 0L) stats::median(prim_cy_fin) else NA_real_
    has_sheet_cal <- any(is_cal[c_idx_vec]) && any(is.finite(ppm[c_idx_vec]) & ppm[c_idx_vec] > 0)
    cal_ppm_v <- ppm[c_idx_vec][is_cal[c_idx_vec] & is.finite(ppm[c_idx_vec]) & ppm[c_idx_vec] > 0]
    ppm_sheet <- if (length(cal_ppm_v) > 0L) stats::median(cal_ppm_v) else if (has_sheet_cal) stats::median(ppm[c_idx_vec][is.finite(ppm[c_idx_vec]) & ppm[c_idx_vec] > 0]) else NA_real_
    for (k_p in seq_along(u_pids)) {
      p_cur <- u_pids[k_p]
      if (p_cur %in% kept_pids) next
      if (bas_counts[k_p] == 0L || tot_counts[k_p] <= 1L) {
        hit_rows <- c_idx_vec[pid_s03[c_idx_vec] == p_cur]
        h_px_ok <- hit_rows[is.finite(cx_px[hit_rows]) & is.finite(cy_px[hit_rows])]
        p_px_ok <- prim_rows[is.finite(cx_px[prim_rows]) & is.finite(cy_px[prim_rows])]
        h_sh_ok <- hit_rows[is.finite(cx_sh[hit_rows]) & is.finite(cy_sh[hit_rows])]
        p_sh_ok <- prim_rows[is.finite(cx_sh[prim_rows]) & is.finite(cy_sh[prim_rows])]
        h_cy_ok <- hit_rows[is.finite(cy_sh[hit_rows])]
        p_cy_ok <- prim_rows[is.finite(cy_sh[prim_rows])]
        can_merge <- if (length(h_px_ok) > 0L && length(p_px_ok) > 0L && has_sheet_cal && is.finite(ppm_sheet) && ppm_sheet > 0) {
          dx <- outer(cx_px[h_px_ok], cx_px[p_px_ok], "-")
          dy <- outer(cy_px[h_px_ok], cy_px[p_px_ok], "-")
          min_d_mm <- min(sqrt(dx^2 + dy^2) / ppm_sheet)
          is.finite(min_d_mm) && (min_d_mm <= max_dist_mm)
        } else if (length(h_sh_ok) > 0L && length(p_sh_ok) > 0L) {
          dx <- outer(cx_sh[h_sh_ok], cx_sh[p_sh_ok], "-")
          dy <- outer(cy_sh[h_sh_ok], cy_sh[p_sh_ok], "-")
          min_d_rel <- min(sqrt(dx^2 + dy^2))
          is.finite(min_d_rel) && (min_d_rel <= max_rel_span)
        } else if (length(h_px_ok) > 0L && length(p_px_ok) > 0L) {
          norm_span <- max(diff(range(cx_px[c_idx_vec], na.rm = TRUE)), diff(range(cy_px[c_idx_vec], na.rm = TRUE)), 1000.0)
          dx <- outer(cx_px[h_px_ok], cx_px[p_px_ok], "-")
          dy <- outer(cy_px[h_px_ok], cy_px[p_px_ok], "-")
          min_d_rel <- min(sqrt(dx^2 + dy^2)) / norm_span
          is.finite(min_d_rel) && (min_d_rel <= max_rel_span)
        } else if (length(h_cy_ok) > 0L && length(p_cy_ok) > 0L) {
          min_d_cy <- min(abs(outer(cy_sh[h_cy_ok], cy_sh[p_cy_ok], "-")))
          is.finite(min_d_cy) && (min_d_cy <= max(max_rel_span, 0.25))
        } else {
          FALSE
        }
        if (!isTRUE(can_merge)) next
        pid_s03[hit_rows] <- primary_pid
        reassigned_s03 <- reassigned_s03 + length(hit_rows)
        for (hr in hit_rows) {
          if (lz_s03[hr] == "basal" && length(prim_bas_rows) >= 2L) {
            is_upper_stray <- (is.finite(dc_s03[hr]) && dc_s03[hr] > 0.24) ||
              (is.finite(cy_s03[hr]) && is.finite(min_prim_bas_cy) && (cy_s03[hr] < (min_prim_bas_cy - 0.12) || (med_prim_bas_cy >= 0.60 && cy_s03[hr] < 0.45)))
            if (isTRUE(is_upper_stray)) {
              lz_s03[hr] <- "cauline"
              rezoned_strays_s03 <- rezoned_strays_s03 + 1L
            }
          }
        }
      }
    }
  }
  if (reassigned_s03 > 0L) {
    efa_df$plant_individual_id <- pid_s03
    if ("leaf_zone" %in% names(efa_df) && rezoned_strays_s03 > 0L) efa_df$leaf_zone <- lz_s03
    if ("leaf_zone_routed" %in% names(efa_df) && rezoned_strays_s03 > 0L) efa_df$leaf_zone_routed <- ifelse(lz_s03 == "cauline", "cauline", as.character(efa_df$leaf_zone_routed))
    message(sprintf("[Step 03 Multi-Plant Hygiene] Consolidated %d cauline-split / single-leaf stray rows onto primary basal individuals (re-zoned %d upper-stem strays to cauline).", reassigned_s03, rezoned_strays_s03))
  }
  efa_df
}

filter_biologically_variable_dissection_scalars <- function(df, cand_scalars, active_mask = NULL) {
  cand_scalars <- intersect(as.character(cand_scalars), names(df))
  if (length(cand_scalars) == 0L || nrow(df) < 5L) return(character(0))
  act <- if (!is.null(active_mask) && length(active_mask) == nrow(df) && sum(active_mask, na.rm = TRUE) >= 5L) {
    which(as.logical(active_mask))
  } else {
    seq_len(nrow(df))
  }
  cauline_het_scalars <- c(
    "heteroblasty_efa_distance", "basal_vs_cauline_length_ratio",
    "basal_vs_cauline_solidity_delta", "cauline_dissection_index_med",
    "cauline_solidity_med", "Cauline_PC1",
    "Cauline_ECT_PC1", "cauline_ect_lobe_complexity", "heteroblasty_ect_distance"
  )
  keep <- character(0)
  for (sc in cand_scalars) {
    raw_nm <- sub("_std$", "", sc)
    only_std <- grepl("_std$", sc) && !(raw_nm %in% names(df))
    raw_v <- if (raw_nm %in% names(df)) suppressWarnings(as.numeric(df[[raw_nm]][act])) else suppressWarnings(as.numeric(df[[sc]][act]))
    raw_v <- raw_v[is.finite(raw_v)]
    if (length(raw_v) < 5L) next
    if ((raw_nm %in% cauline_het_scalars) && (length(raw_v) < ceiling(0.50 * length(act)))) next
    sd_v <- stats::sd(raw_v)
    q_span <- as.numeric(diff(stats::quantile(raw_v, probs = c(0.10, 0.90), names = FALSE)))
    ok <- if (only_std) {
      (is.finite(sd_v) && sd_v >= 0.25) && (is.finite(q_span) && q_span >= 0.50)
    } else if (raw_nm %in% c("solidity", "q25_solidity", "cauline_solidity_med")) {
      (is.finite(sd_v) && sd_v >= 0.045) || any(raw_v < 0.75) || (is.finite(q_span) && q_span >= 0.10)
    } else if (raw_nm %in% c("dissection_index", "cauline_dissection_index_med")) {
      (is.finite(sd_v) && sd_v >= 0.12) || any(raw_v > 1.35)
    } else if (raw_nm %in% c("within_rosette_efa_dispersion", "heteroblasty_efa_distance", "heteroblasty_ect_distance")) {
      (is.finite(sd_v) && sd_v >= 0.015) && (sum(raw_v > 0) >= 5L)
    } else if (raw_nm %in% c("basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "Cauline_PC1",
                             "ECT_PC1", "ECT_PC2", "ECT_PC3", "Cauline_ECT_PC1", "ect_lobe_complexity", "cauline_ect_lobe_complexity")) {
      is.finite(sd_v) && sd_v >= 0.015
    } else if (raw_nm == "basal_cordate_index") {
      (is.finite(sd_v) && sd_v >= 0.025) || (is.finite(q_span) && q_span >= 0.06)
    } else if (raw_nm == "terminal_lobe_area_frac") {
      (is.finite(sd_v) && sd_v >= 0.04) && (sum(raw_v < 0.95) >= 5L)
    } else if (raw_nm == "margin_serration_index") {
      (is.finite(sd_v) && sd_v >= 0.008) || (is.finite(q_span) && q_span >= 0.02)
    } else if (raw_nm %in% c("margin_tooth_acuteness", "radiate_head_fraction")) {
      (is.finite(sd_v) && sd_v >= 0.08) && (is.finite(q_span) && q_span >= 0.10)
    } else if (raw_nm %in% c("base_angle_deg", "apex_angle_deg")) {
      is.finite(sd_v) && sd_v >= 4.0
    } else if (raw_nm == "lamina_curvature") {
      is.finite(sd_v) && sd_v >= 0.35
    } else if (raw_nm == "blade_cie_l_norm" && max(abs(raw_v), na.rm = TRUE) <= 2.0) {
      is.finite(sd_v) && sd_v >= 0.03
    } else if (grepl("(cie_|blade_cie_)", raw_nm)) {
      is.finite(sd_v) && sd_v >= 1.25
    } else {
      (is.finite(sd_v) && sd_v >= 0.02) && (is.finite(q_span) && q_span >= 0.04)
    }
    if (isTRUE(ok)) keep <- c(keep, sc)
  }
  keep
}

prune_collinear_dissection_scalars <- function(
  df,
  base_cols,
  cand_scalars,
  active_mask = NULL,
  group_vec = NULL,
  max_r = 0.75
) {
  base_cols <- intersect(as.character(base_cols), names(df))
  cand_scalars <- setdiff(intersect(as.character(cand_scalars), names(df)), base_cols)
  if (length(cand_scalars) == 0L || nrow(df) < 5L) return(character(0))
  act <- if (!is.null(active_mask) && length(active_mask) == nrow(df) && sum(active_mask, na.rm = TRUE) >= 5L) {
    which(as.logical(active_mask))
  } else {
    seq_len(nrow(df))
  }

  tier3_traits <- c(
    "ECT_PC1", "ECT_PC2", "ECT_PC3", "basal_cordate_index",
    "q25_solidity", "solidity", "dissection_index", "terminal_lobe_area_frac",
    "AR_inv", "radiate_head_fraction", "petiole_to_blade_ratio"
  )
  tier2_traits <- c(
    "margin_serration_index", "margin_tooth_acuteness",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "abaxial_cie_l_sheet_delta", "base_angle_deg", "apex_angle_deg",
    "lamina_curvature", "log1p_capitulum_count", "scape_to_basal_leaf_ratio"
  )

  # Rank candidate scalars empirically by univariate ANOVA F on active reference anchors (MorphoTools2 protocol),
  # falling back in unsupervised/single-taxon runs to deterministic biological tier + unstandardized relative dispersion.
  f_scores <- vapply(cand_scalars, function(sc) {
    v <- suppressWarnings(as.numeric(df[[sc]][act]))
    fin <- is.finite(v)
    if (sum(fin) < 5L || stats::sd(v[fin]) <= 1e-8) return(-1.0)
    if (!is.null(group_vec) && length(group_vec) == nrow(df)) {
      g <- droplevels(factor(as.character(group_vec[act][fin])))
      if (nlevels(g) >= 2L && min(table(g)) >= 2L) {
        fv <- tryCatch(summary(stats::aov(v[fin] ~ g))[[1]][["F value"]][1L], error = function(e) NA_real_)
        if (is.finite(fv)) return(as.numeric(fv))
      }
    }
    raw_nm <- sub("_std$", "", sc)
    tier_priority <- if (raw_nm %in% tier3_traits) {
      3.0
    } else if (raw_nm %in% tier2_traits) {
      2.0
    } else {
      1.0
    }
    rel_disp <- NA_real_
    if (raw_nm %in% names(df) && raw_nm != sc) {
      raw_v <- suppressWarnings(as.numeric(df[[raw_nm]][act]))
      raw_v <- raw_v[is.finite(raw_v)]
      if (length(raw_v) >= 5L) {
        rel_disp <- min(50.0, stats::IQR(raw_v, na.rm = TRUE) / max(stats::median(abs(raw_v), na.rm = TRUE), 0.10))
      }
    }
    if (!is.finite(rel_disp)) {
      rel_disp <- min(50.0, stats::var(v[fin]))
    }
    (tier_priority * 100.0) + rel_disp
  }, numeric(1))

  valid_cands <- cand_scalars[f_scores >= 0.0]
  if (length(valid_cands) == 0L) return(character(0))
  ord_cands <- valid_cands[order(f_scores[valid_cands], decreasing = TRUE)]

  accepted <- base_cols
  kept_scalars <- character(0)
  for (sc in ord_cands) {
    v_sc <- suppressWarnings(as.numeric(df[[sc]][act]))
    if (sum(is.finite(v_sc)) < 5L || stats::sd(v_sc, na.rm = TRUE) <= 1e-8) next
    too_corr <- FALSE
    for (ac in accepted) {
      v_ac <- suppressWarnings(as.numeric(df[[ac]][act]))
      ok <- is.finite(v_sc) & is.finite(v_ac)
      if (sum(ok) >= 5L && stats::sd(v_ac[ok]) > 1e-8 && stats::sd(v_sc[ok]) > 1e-8) {
        r_val <- suppressWarnings(abs(stats::cor(v_sc[ok], v_ac[ok])))
        if (is.finite(r_val) && r_val >= max_r) {
          too_corr <- TRUE
          break
        }
      }
    }
    if (!too_corr) {
      accepted <- c(accepted, sc)
      kept_scalars <- c(kept_scalars, sc)
    }
  }
  kept_scalars
}

fast_euclidean_permanova_trace <- function(X_mat, grp_vec, permutations = 999L, seed = 42L) {
  X <- as.matrix(X_mat)
  storage.mode(X) <- "double"
  g_fac <- droplevels(factor(as.character(grp_vec)))
  ok <- stats::complete.cases(X) & !is.na(g_fac)
  X <- X[ok, , drop = FALSE]
  g_fac <- droplevels(g_fac[ok])
  N <- nrow(X)
  g <- nlevels(g_fac)
  if (N <= g || g < 2L || ncol(X) < 1L) {
    return(list(F = NA_real_, R2 = NA_real_, p_value = NA_real_, df_group = max(0L, g - 1L), df_res = max(0L, N - g)))
  }

  X_c <- scale(X, center = TRUE, scale = FALSE)
  ss_tot <- sum(X_c^2)
  if (!is.finite(ss_tot) || ss_tot <= 1e-12) {
    return(list(F = 0.0, R2 = 0.0, p_value = 1.0, df_group = g - 1L, df_res = N - g))
  }

  g_int <- as.integer(g_fac)
  n_vec <- as.numeric(tabulate(g_int, nbins = g))
  grp_sums <- rowsum(X_c, g_int, reorder = TRUE)
  ss_between <- sum(rowSums(grp_sums^2) / n_vec)
  ss_within <- max(0.0, ss_tot - ss_between)
  df1 <- as.numeric(g - 1L)
  df2 <- as.numeric(N - g)
  r2_obs <- ss_between / ss_tot
  f_obs <- if (ss_within > 1e-12) (ss_between / df1) / (ss_within / df2) else Inf

  B <- max(0L, as.integer(permutations))
  p_val <- NA_real_
  if (B >= 1L) {
    set.seed(as.integer(seed))
    tol <- ss_between - 1e-12 * max(1.0, ss_tot)
    if (N * B >= 250000L && exists("packera_mclapply", mode = "function")) {
      perm_list <- replicate(B, sample.int(N), simplify = FALSE)
      n_w <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(16L, B) else 4L
      chunks <- split(seq_len(B), cut(seq_len(B), breaks = min(n_w, B), labels = FALSE))
      cnts <- packera_mclapply(chunks, function(idx_b) {
        c_loc <- 0L
        for (bi in idx_b) {
          sums_b <- rowsum(X_c[perm_list[[bi]], , drop = FALSE], g_int, reorder = TRUE)
          if (sum(rowSums(sums_b^2) / n_vec) >= tol) c_loc <- c_loc + 1L
        }
        c_loc
      }, mc.cores = length(chunks))
      ge_cnt <- sum(unlist(cnts))
    } else {
      ge_cnt <- 0L
      for (b in seq_len(B)) {
        g_perm <- sample.int(N)
        sums_b <- rowsum(X_c[g_perm, , drop = FALSE], g_int, reorder = TRUE)
        ss_b <- sum(rowSums(sums_b^2) / n_vec)
        if (ss_b >= tol) ge_cnt <- ge_cnt + 1L
      }
    }
    p_val <- (1.0 + ge_cnt) / (1.0 + B)
  }

  list(
    F = as.numeric(f_obs),
    R2 = as.numeric(r2_obs),
    p_value = as.numeric(p_val),
    SumOfSqs = as.numeric(ss_between),
    Residual_SS = as.numeric(ss_within),
    Total_SS = as.numeric(ss_tot),
    df_group = as.integer(df1),
    df_res = as.integer(df2)
  )
}
