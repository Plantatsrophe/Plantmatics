if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step04_restore_a.R
# Step 04 raw Step 03 restoration from leaf_efa_harmonics (part A: matching, harmonics, and scale calibration).
# ==============================================================================

.step04_leaf_has_fin_r <- function(leaf_df, cn) {
  (cn %in% names(leaf_df)) && any(is.finite(suppressWarnings(as.numeric(leaf_df[[cn]]))))
}

.step04_agg_med_col_r <- function(leaf_df, prim_idx, matched_rows, n_df, cn, digits = NULL, pos_only = FALSE, cal_vec = NULL, subset_vec = NULL) {
  out <- rep(NA_real_, n_df)
  if (!(cn %in% names(leaf_df))) return(out)
  v_all <- suppressWarnings(as.numeric(leaf_df[[cn]]))
  for (i in which(matched_rows)) {
    idx <- prim_idx[[i]]
    if (!is.null(cal_vec)) idx <- idx[cal_vec[idx]]
    if (!is.null(subset_vec)) idx <- idx[subset_vec[idx]]
    if (length(idx) == 0L) next
    vf <- v_all[idx]
    vf <- if (isTRUE(pos_only)) vf[is.finite(vf) & vf > 0.0] else vf[is.finite(vf)]
    if (length(vf) > 0L) {
      mv <- stats::median(vf)
      out[i] <- if (!is.null(digits)) round(mv, digits) else mv
    }
  }
  out
}

.restore_raw_step03_harmonics_and_scale_r <- function(conn, rid, df) {
  empty_res <- list(matched = FALSE, df = df, restored_mask = list())
  if (is.null(conn) || is.null(df) || !is.data.frame(df) || nrow(df) == 0L ||
      !db_table_exists_safe(conn, "leaf_efa_harmonics")) {
    return(empty_res)
  }
  leaf_df <- tryCatch(
    DBI::dbGetQuery(conn, "SELECT * FROM leaf_efa_harmonics WHERE run_id = $1;", params = list(rid)),
    error = function(e) NULL
  )
  if (is.null(leaf_df) || nrow(leaf_df) == 0L) {
    return(empty_res)
  }
  leaf_df <- coerce_sqlite_logical_cols(leaf_df)
  n_df <- nrow(df)

  norm_pid <- function(d) {
    if ("plant_individual_id" %in% names(d)) {
      v <- suppressWarnings(as.integer(d$plant_individual_id))
      v[is.na(v)] <- 0L
      v
    } else {
      rep(0L, nrow(d))
    }
  }
  df_keys <- paste(trimws(as.character(df$catalogNumber)), norm_pid(df), sep = "\r")
  l_keys <- paste(trimws(as.character(leaf_df$catalogNumber)), norm_pid(leaf_df), sep = "\r")
  leaf_groups <- split(seq_len(nrow(leaf_df)), l_keys)

  type_mask_fn <- function(d) {
    if (exists(".identify_type_specimen_vec", mode = "function")) {
      res_t <- tryCatch(as.logical(.identify_type_specimen_vec(d)$is_type_specimen), error = function(e) NULL)
      if (!is.null(res_t) && length(res_t) == nrow(d)) {
        res_t[is.na(res_t)] <- FALSE
        return(res_t)
      }
    }
    n_r <- nrow(d)
    b_t <- if ("is_type_specimen" %in% names(d)) {
      v <- as.logical(d$is_type_specimen)
      !is.na(v) & v
    } else {
      rep(FALSE, n_r)
    }
    raw_t <- if ("type_status" %in% names(d)) {
      trimws(as.character(d$type_status))
    } else if ("typeStatus" %in% names(d)) {
      trimws(as.character(d$typeStatus))
    } else {
      rep("", n_r)
    }
    raw_t[is.na(raw_t)] <- ""
    neg_m <- grepl("not a type|non-type|non type|false type|no type|unverified|type locality|type species|type collection", raw_t, ignore.case = TRUE)
    pos_m <- !neg_m & grepl("(^|[^a-z])(holotype|lectotype|neotype|epitype|isotype|isolectotype|isoneotype|syntype|isosyntype|paratype|isoparatype|type|cotype)s?([^a-z]|$)", raw_t, ignore.case = TRUE)
    pos_m | (b_t & !neg_m)
  }

  df_is_type <- type_mask_fn(df)
  l_is_type <- type_mask_fn(leaf_df)
  l_out <- if ("is_contour_outlier" %in% names(leaf_df)) {
    !is.na(leaf_df$is_contour_outlier) & as.logical(leaf_df$is_contour_outlier)
  } else {
    rep(FALSE, nrow(leaf_df))
  }
  l_lz <- if ("leaf_zone" %in% names(leaf_df)) {
    tolower(trimws(as.character(leaf_df$leaf_zone)))
  } else {
    rep("basal", nrow(leaf_df))
  }

  prim_idx <- vector("list", n_df)
  matched_rows <- logical(n_df)
  for (i in seq_len(n_df)) {
    idx <- leaf_groups[[df_keys[i]]]
    if (is.null(idx) || length(idx) == 0L) next
    matched_rows[i] <- TRUE
    row_type <- l_is_type[idx] | isTRUE(df_is_type[i])
    eff_out <- l_out[idx] & !row_type
    if (any(!eff_out)) {
      idx <- idx[!eff_out]
    }
    lz_i <- l_lz[idx]
    is_basal_primary <- !is.na(lz_i) & (lz_i == "basal")
    if (any(is_basal_primary)) {
      idx <- idx[is_basal_primary]
    }
    prim_idx[[i]] <- idx
  }
  if (!any(matched_rows)) {
    return(empty_res)
  }

  restored_mask <- list()
  has_fin_in_leaf <- function(cn) .step04_leaf_has_fin_r(leaf_df, cn)
  agg_med_col <- function(cn, digits = NULL, pos_only = FALSE, cal_vec = NULL, subset_vec = NULL) {
    .step04_agg_med_col_r(leaf_df, prim_idx, matched_rows, n_df, cn, digits = digits, pos_only = pos_only, cal_vec = cal_vec, subset_vec = subset_vec)
  }
  df_has_asym <- any(grepl("^[BC][0-9]+$", names(df)))
  allowed_harm_pool <- if (df_has_asym) PACKERA_DB_HARMONIC_COLS else PACKERA_DB_SYM_HARMONIC_COLS
  harm_cols <- intersect(allowed_harm_pool, names(leaf_df))
  harm_cols <- harm_cols[vapply(harm_cols, has_fin_in_leaf, logical(1L))]
  for (hc in harm_cols) {
    vals <- agg_med_col(hc)
    if (!(hc %in% names(df))) df[[hc]] <- rep(NA_real_, n_df)
    df[[hc]][matched_rows] <- vals[matched_rows]
    restored_mask[[hc]] <- matched_rows
  }

  for (spec in list(
    list(col = "ECT_PC1", d = 6L), list(col = "ECT_PC2", d = 6L), list(col = "ECT_PC3", d = 6L),
    list(col = "ect_lobe_complexity", d = 4L),
    list(col = "solidity", d = 4L), list(col = "dissection_index", d = 4L), list(col = "aspect_ratio", d = 4L),
    list(col = "base_angle_deg", d = 4L), list(col = "apex_angle_deg", d = 4L), list(col = "lamina_curvature", d = 4L),
    list(col = "basal_cordate_index", d = 4L), list(col = "margin_serration_index", d = 4L),
    list(col = "margin_tooth_acuteness", d = 4L), list(col = "terminal_lobe_area_frac", d = 4L)
  )) {
    cn <- spec$col
    if (has_fin_in_leaf(cn)) {
      vals <- agg_med_col(cn, digits = spec$d)
      ok_v <- matched_rows & is.finite(vals)
      if (any(ok_v)) {
        if (!(cn %in% names(df))) df[[cn]] <- rep(NA_real_, n_df)
        df[[cn]][ok_v] <- vals[ok_v]
        restored_mask[[cn]] <- ok_v
      }
    }
  }

  if (has_fin_in_leaf("solidity")) {
    q25_vals <- rep(NA_real_, n_df)
    sol_all <- suppressWarnings(as.numeric(leaf_df$solidity))
    for (i in which(matched_rows)) {
      sf <- sol_all[prim_idx[[i]]]
      sf <- sf[is.finite(sf)]
      if (length(sf) > 0L) {
        q25_vals[i] <- round(as.numeric(stats::quantile(sf, probs = 0.25, na.rm = TRUE, names = FALSE)), 4)
      }
    }
    ok_q <- matched_rows & is.finite(q25_vals)
    if (any(ok_q)) {
      if (!("q25_solidity" %in% names(df))) df$q25_solidity <- rep(NA_real_, n_df)
      df$q25_solidity[ok_q] <- q25_vals[ok_q]
      restored_mask[["q25_solidity"]] <- ok_q
    }
  }

  uncal_ss_set <- c("mp_estimated", "predicted_from_megapixels", "tick_lattice_unboxed", "uncalibrated")
  l_cal_raw <- if ("is_calibrated_mm" %in% names(leaf_df)) {
    as.logical(leaf_df$is_calibrated_mm)
  } else {
    rep(NA, nrow(leaf_df))
  }
  ss_l <- if ("scale_source" %in% names(leaf_df)) {
    s_tmp <- tolower(trimws(as.character(leaf_df$scale_source)))
    s_tmp[is.na(s_tmp) | s_tmp %in% c("na", "null")] <- ""
    s_tmp
  } else {
    rep("", nrow(leaf_df))
  }
  l_cal <- ifelse(is.na(l_cal_raw), !(ss_l %in% uncal_ss_set), l_cal_raw) & !(ss_l %in% uncal_ss_set)
  d_cal <- if ("is_calibrated_mm" %in% names(df)) {
    !is.na(df$is_calibrated_mm) & as.logical(df$is_calibrated_mm)
  } else {
    rep(TRUE, n_df)
  }
  if ("scale_source" %in% names(df)) {
    ss_d <- tolower(trimws(as.character(df$scale_source)))
    ss_d[is.na(ss_d) | ss_d %in% c("na", "null")] <- ""
    d_cal <- d_cal & !(ss_d %in% uncal_ss_set)
  }
  if (("scale_source" %in% names(leaf_df)) || ("is_calibrated_mm" %in% names(leaf_df))) {
    ss_rest <- rep(NA_character_, n_df)
    ok_ss <- rep(FALSE, n_df)
    ok_cal_meta <- rep(FALSE, n_df)
    for (i in which(matched_rows)) {
      idx <- prim_idx[[i]]
      if (length(idx) == 0L) next
      has_cal_meta_i <- any(!is.na(l_cal_raw[idx]) | nzchar(ss_l[idx]))
      if (has_cal_meta_i) {
        ok_cal_meta[i] <- TRUE
        has_cal_i <- any(l_cal[idx])
        d_cal[i] <- d_cal[i] & has_cal_i
      }
      if (any(nzchar(ss_l[idx]))) {
        ok_ss[i] <- TRUE
        if (d_cal[i]) {
          cal_ss <- ss_l[idx[l_cal[idx]]]
          cal_ss <- cal_ss[nzchar(cal_ss)]
          ss_rest[i] <- if (length(cal_ss) > 0L) cal_ss[1L] else "ruler_calibrated"
        } else if ("tick_lattice_unboxed" %in% ss_l[idx]) {
          ss_rest[i] <- "tick_lattice_unboxed"
        } else if (any(ss_l[idx] %in% c("mp_estimated", "predicted_from_megapixels"))) {
          ss_rest[i] <- "mp_estimated"
        } else {
          ss_rest[i] <- "uncalibrated"
        }
      }
    }
    if (any(ok_ss)) {
      if (!("scale_source" %in% names(df))) df$scale_source <- rep(NA_character_, n_df)
      df$scale_source[ok_ss] <- ss_rest[ok_ss]
      restored_mask[["scale_source"]] <- ok_ss
    }
    if (any(ok_cal_meta)) {
      if (!("is_calibrated_mm" %in% names(df))) df$is_calibrated_mm <- rep(NA, n_df)
      df$is_calibrated_mm[ok_cal_meta] <- d_cal[ok_cal_meta]
      restored_mask[["is_calibrated_mm"]] <- ok_cal_meta
    }
  }
  if (has_fin_in_leaf("centroid_size") || has_fin_in_leaf("log_centsize")) {
    cs_vals <- if (has_fin_in_leaf("centroid_size")) agg_med_col("centroid_size", digits = 4L, pos_only = TRUE, cal_vec = l_cal) else rep(NA_real_, n_df)
    lcs_fb <- if (has_fin_in_leaf("log_centsize")) agg_med_col("log_centsize", digits = 6L, cal_vec = l_cal) else rep(NA_real_, n_df)
    cs_vals[!d_cal] <- NA_real_
    lcs_vals <- ifelse(is.finite(cs_vals) & cs_vals > 0.0, round(log10(cs_vals), 6), ifelse(d_cal, lcs_fb, NA_real_))
    if (has_fin_in_leaf("centroid_size")) {
      if (!("centroid_size" %in% names(df))) df$centroid_size <- rep(NA_real_, n_df)
      df$centroid_size[matched_rows] <- cs_vals[matched_rows]
      restored_mask[["centroid_size"]] <- matched_rows
    }
    if (!("log_centsize" %in% names(df))) df$log_centsize <- rep(NA_real_, n_df)
    df$log_centsize[matched_rows] <- lcs_vals[matched_rows]
    restored_mask[["log_centsize"]] <- matched_rows
  }

  has_ltl_px <- has_fin_in_leaf("lamina_trace_length_px")
  has_ltl_mm <- has_fin_in_leaf("lamina_trace_length_mm")
  if (has_ltl_mm || has_ltl_px) {
    if (!("lamina_trace_length_mm" %in% names(leaf_df))) {
      leaf_df$lamina_trace_length_mm <- rep(NA_real_, nrow(leaf_df))
    }
    ltl_mm_l <- suppressWarnings(as.numeric(leaf_df$lamina_trace_length_mm))
    if (has_ltl_px && ("pixels_per_mm" %in% names(leaf_df))) {
      l_ppm <- suppressWarnings(as.numeric(leaf_df$pixels_per_mm))
      l_px <- suppressWarnings(as.numeric(leaf_df$lamina_trace_length_px))
      fill_ltl <- !is.finite(ltl_mm_l) & l_cal & is.finite(l_ppm) & (l_ppm > 0.0) & is.finite(l_px)
      if (any(fill_ltl)) {
        ltl_mm_l[fill_ltl] <- round(l_px[fill_ltl] / l_ppm[fill_ltl], 4)
      }
    }
    ltl_mm_l[!l_cal] <- NA_real_
    leaf_df$lamina_trace_length_mm <- ltl_mm_l
    ltl_vals <- rep(NA_real_, n_df)
    for (i in which(matched_rows & d_cal)) {
      idx <- prim_idx[[i]]
      idx <- idx[l_cal[idx]]
      if (length(idx) == 0L) next
      vf <- ltl_mm_l[idx]
      vf <- vf[is.finite(vf) & vf > 0.0]
      if (length(vf) > 0L) {
        ltl_vals[i] <- round(as.numeric(stats::quantile(vf, probs = 0.75, type = 7L, names = FALSE)), 4)
      }
    }
    if (!("lamina_trace_length_mm" %in% names(df))) df$lamina_trace_length_mm <- rep(NA_real_, n_df)
    df$lamina_trace_length_mm[matched_rows] <- ltl_vals[matched_rows]
    restored_mask[["lamina_trace_length_mm"]] <- matched_rows
  }
  list(
    matched = TRUE,
    df = df,
    leaf_df = leaf_df,
    n_df = n_df,
    prim_idx = prim_idx,
    matched_rows = matched_rows,
    df_is_type = df_is_type,
    harm_cols = harm_cols,
    restored_mask = restored_mask
  )
}
