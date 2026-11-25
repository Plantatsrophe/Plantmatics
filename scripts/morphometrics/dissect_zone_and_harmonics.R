if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/dissect_zone_and_harmonics.R
# Positional basal/cauline zone partitioning, routing, and EFA harmonic calibration
# ==============================================================================

compute_polygon_area_perim_pure <- function(coo) {
  if (is.null(coo) || !is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 3L) {
    return(list(area = NA_real_, perim = NA_real_, circularity = NA_real_, solidity = NA_real_))
  }
  x <- as.numeric(coo[, 1L]); y <- as.numeric(coo[, 2L])
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    return(list(area = NA_real_, perim = NA_real_, circularity = NA_real_, solidity = NA_real_))
  }
  xn <- c(x[-1L], x[1L]); yn <- c(y[-1L], y[1L])
  area <- 0.5 * abs(sum(x * yn - xn * y))
  perim <- sum(sqrt((xn - x)^2 + (yn - y)^2))
  circ <- if (is.finite(perim) && perim > 1e-8 && is.finite(area) && area > 0) {
    min(1.0, (4.0 * pi * area) / (perim^2))
  } else {
    NA_real_
  }
  sol <- NA_real_
  if (is.finite(area) && area > 1e-8) {
    ch <- tryCatch(grDevices::chull(x, y), error = function(e) integer(0))
    if (length(ch) >= 3L) {
      hx <- x[ch]; hy <- y[ch]
      ha <- 0.5 * abs(sum(hx * c(hy[-1L], hy[1L]) - c(hx[-1L], hx[1L]) * hy))
      if (is.finite(ha) && ha > 1e-8) sol <- min(1.0, area / ha)
    }
  }
  list(area = area, perim = perim, circularity = circ, solidity = sol)
}

compute_dissection_index_from_coo <- function(coo = NULL, area_val = NULL, perim_val = NULL) {
  a <- if (!is.null(area_val) && length(area_val) == 1L) suppressWarnings(as.numeric(area_val)) else NA_real_
  p <- if (!is.null(perim_val) && length(perim_val) == 1L) suppressWarnings(as.numeric(perim_val)) else NA_real_
  if ((!is.finite(a) || a <= 1e-8 || !is.finite(p) || p <= 1e-8) && !is.null(coo)) {
    gp <- compute_polygon_area_perim_pure(coo)
    a <- gp$area; p <- gp$perim
  }
  if (!is.finite(a) || a <= 1e-8 || !is.finite(p) || p <= 1e-8) return(NA_real_)
  round(max(1.0, p / (2.0 * sqrt(pi * a))), 4)
}

classify_and_partition_leaf_zones <- function(
  fac_df,
  min_cy_plant = 0.40,
  min_cy_sheet = 0.42,
  fallback_band = 0.15
) {
  if (is.null(fac_df) || nrow(fac_df) == 0L) return(fac_df)
  n <- nrow(fac_df)
  if (!("plant_individual_id" %in% names(fac_df))) fac_df$plant_individual_id <- 0L
  fac_df$plant_individual_id <- ifelse(
    is.na(fac_df$plant_individual_id),
    0L,
    suppressWarnings(as.integer(fac_df$plant_individual_id))
  )
  fac_df$plant_individual_id[is.na(fac_df$plant_individual_id)] <- 0L

  raw_zone <- if ("leaf_zone" %in% names(fac_df)) tolower(trimws(as.character(fac_df$leaf_zone))) else rep(NA_character_, n)
  raw_zone[!(raw_zone %in% c("basal", "cauline"))] <- NA_character_
  zone_conf <- if ("leaf_zone_confidence" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$leaf_zone_confidence)) else rep(NA_real_, n)
  d_crown <- if ("dist_from_crown_norm" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$dist_from_crown_norm)) else rep(NA_real_, n)
  cy_p <- if ("cy_rel_plant_bbox" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$cy_rel_plant_bbox)) else rep(NA_real_, n)
  cy_s <- if ("cy_rel_sheet" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$cy_rel_sheet)) else rep(NA_real_, n)

  leaf_zone <- raw_zone
  is_fallback_vec <- rep(FALSE, n)
  is_caul_only_vec <- rep(FALSE, n)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)

  for (g_idx in idx_by_grp) {
    lam_ok <- g_idx
    if (length(lam_ok) == 0L) next
    all_explicit_cauline <- all(!is.na(raw_zone[lam_ok]) & raw_zone[lam_ok] == "cauline")

    for (ii in g_idx) {
      if (!is.na(leaf_zone[ii])) {
        if (!is.finite(zone_conf[ii])) zone_conf[ii] <- 0.90
        next
      }
      if (is.finite(d_crown[ii])) {
        s_b <- if (d_crown[ii] <= 0.18) 0.82 else if (d_crown[ii] > 0.24) 0.20 else 0.50
        if (d_crown[ii] > 0.18 && d_crown[ii] <= 0.24) {
          if (is.finite(cy_p[ii])) {
            s_b <- s_b + (if (cy_p[ii] >= min_cy_plant) 0.16 else -0.16)
          } else if (is.finite(cy_s[ii])) {
            s_b <- s_b + (if (cy_s[ii] >= min_cy_sheet) 0.16 else -0.16)
          }
        }
        p_b <- min(1.0, max(0.0, s_b))
        leaf_zone[ii] <- if (p_b >= 0.50) "basal" else "cauline"
        zone_conf[ii] <- round(max(0.55, abs(p_b - 0.50) * 2.0), 4)
      } else if (is.finite(cy_p[ii]) || is.finite(cy_s[ii])) {
        in_prim <- if (is.finite(cy_p[ii])) (cy_p[ii] >= min_cy_plant) else (cy_s[ii] >= min_cy_sheet)
        leaf_zone[ii] <- if (isTRUE(in_prim)) "basal" else "cauline"
        zone_conf[ii] <- if (is.finite(zone_conf[ii])) zone_conf[ii] else 0.80
      } else {
        leaf_zone[ii] <- "basal"
        zone_conf[ii] <- if (is.finite(zone_conf[ii])) zone_conf[ii] else 0.75
      }
    }

    has_basal <- any(leaf_zone[lam_ok] == "basal", na.rm = TRUE)
    if (!has_basal) {
      has_dc <- is.finite(d_crown[lam_ok])
      crown_verified_cauline_only <- any(has_dc) && (min(d_crown[lam_ok[has_dc]], na.rm = TRUE) > 0.28)
      if (all_explicit_cauline || crown_verified_cauline_only) {
        is_caul_only_vec[lam_ok] <- TRUE
        is_fallback_vec[lam_ok] <- TRUE
      } else {
        has_cy_p <- is.finite(cy_p[lam_ok])
        has_cy_s <- !has_cy_p & is.finite(cy_s[lam_ok])
        if (any(has_cy_p)) {
          c_max <- max(cy_p[lam_ok[has_cy_p]], na.rm = TRUE)
          fb_hit <- has_cy_p & (cy_p[lam_ok] >= (c_max - fallback_band))
        } else if (any(has_cy_s)) {
          c_max <- max(cy_s[lam_ok[has_cy_s]], na.rm = TRUE)
          fb_hit <- has_cy_s & (cy_s[lam_ok] >= (c_max - fallback_band))
        } else if (any(has_dc)) {
          d_min <- min(d_crown[lam_ok[has_dc]], na.rm = TRUE)
          fb_hit <- has_dc & (d_crown[lam_ok] <= (d_min + fallback_band))
        } else {
          fb_hit <- seq_along(lam_ok) == 1L
        }
        if (any(fb_hit)) {
          leaf_zone[lam_ok[fb_hit]] <- "basal"
          is_fallback_vec[lam_ok[fb_hit]] <- TRUE
        }
      }
    }
  }

  leaf_zone_routed <- ifelse(
    !is.na(leaf_zone) & leaf_zone == "basal" & (!is.finite(zone_conf) | zone_conf >= 0.65),
    "basal",
    ifelse(
      !is.na(leaf_zone) & leaf_zone == "cauline" & (!is.finite(zone_conf) | zone_conf >= 0.65),
      "cauline",
      "ambiguous_transitional"
    )
  )
  for (g_idx in idx_by_grp) {
    has_hi_basal <- any(leaf_zone_routed[g_idx] == "basal", na.rm = TRUE)
    bas_cand <- g_idx[!is.na(leaf_zone[g_idx]) & leaf_zone[g_idx] == "basal"]
    if (!has_hi_basal && length(bas_cand) > 0L && !all(is_caul_only_vec[g_idx])) {
      sc_c <- ifelse(is.finite(zone_conf[bas_cand]), zone_conf[bas_cand], 0.50)
      leaf_zone_routed[bas_cand[which.max(sc_c)]] <- "basal"
    }
  }

  fac_df$leaf_zone <- leaf_zone
  fac_df$leaf_zone_routed <- leaf_zone_routed
  fac_df$leaf_zone_confidence <- zone_conf
  fac_df$dist_from_crown_norm <- d_crown
  fac_df$is_basal_fallback <- is_fallback_vec
  fac_df$is_cauline_only_individual <- is_caul_only_vec
  fac_df
}

filter_basal_lamina_contours <- function(
  fac_df,
  coo_list,
  min_cy_plant = 0.40,
  min_cy_sheet = 0.42,
  fallback_band = 0.15,
  retain_cauline_tagged = FALSE
) {
  if (is.null(fac_df) || nrow(fac_df) == 0L || length(coo_list) == 0L) {
    return(list(fac_df = fac_df, coo_list = coo_list, stats = list()))
  }
  n <- nrow(fac_df)
  di_vec <- if ("dissection_index" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$dissection_index)) else rep(NA_real_, n)
  sol_vec <- if ("solidity" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$solidity)) else rep(NA_real_, n)
  circ_vec <- if ("circularity" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$circularity)) else rep(NA_real_, n)

  miss_geom <- which(!is.finite(di_vec) | !is.finite(sol_vec) | !is.finite(circ_vec))
  if (length(miss_geom) > 0L) {
    for (i in miss_geom) {
      gp <- compute_polygon_area_perim_pure(coo_list[[i]])
      if (!is.finite(sol_vec[i])) sol_vec[i] <- round(gp$solidity, 4)
      if (!is.finite(circ_vec[i])) circ_vec[i] <- round(gp$circularity, 4)
      if (!is.finite(di_vec[i])) di_vec[i] <- compute_dissection_index_from_coo(area_val = gp$area, perim_val = gp$perim)
    }
  }

  fac_df$solidity <- sol_vec
  fac_df$circularity <- circ_vec
  fac_df$dissection_index <- di_vec
  raw_dis <- if ("is_dissected" %in% names(fac_df)) as.logical(fac_df$is_dissected) else rep(NA, n)
  fac_df$is_dissected <- ifelse(
    !is.na(raw_dis),
    raw_dis | (is.finite(sol_vec) & sol_vec < 0.72),
    (is.finite(sol_vec) & sol_vec < 0.72)
  )

  fac_df <- classify_and_partition_leaf_zones(
    fac_df,
    min_cy_plant = min_cy_plant,
    min_cy_sheet = min_cy_sheet,
    fallback_band = fallback_band
  )
  is_fallback_vec <- as.logical(fac_df$is_basal_fallback)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)
  n_fb_individuals <- sum(vapply(idx_by_grp, function(g_idx) any(is_fallback_vec[g_idx] & !fac_df$is_cauline_only_individual[g_idx]), logical(1)))

  basal_keep <- (fac_df$leaf_zone == "basal") | is_fallback_vec
  for (g_idx in idx_by_grp) {
    if (!any(basal_keep[g_idx]) && length(g_idx) > 0L) {
      basal_keep[g_idx[1L]] <- TRUE
    }
  }
  if (!any(basal_keep)) basal_keep <- rep(TRUE, n)

  n_area_drop <- 0L
  n_caul_drop <- sum(!basal_keep)
  keep_mask <- if (isTRUE(retain_cauline_tagged)) rep(TRUE, n) else basal_keep

  if (n_caul_drop > 0L || n_fb_individuals > 0L) {
    message(sprintf(
      paste0(
        "Positional Homology & Basal Lamina Gate: classified %d / %d primary basal blades ",
        "(cauline-partitioned = %d; lowest-lamina fallback individuals = %d)."
      ),
      sum(basal_keep), n, n_caul_drop, n_fb_individuals
    ))
  }
  fac_sub <- fac_df[keep_mask, , drop = FALSE]
  coo_sub <- coo_list[keep_mask]
  rownames(fac_sub) <- names(coo_sub)
  list(
    fac_df = fac_sub,
    coo_list = coo_sub,
    stats = list(
      n_input = n,
      n_retained = sum(keep_mask),
      n_basal_retained = sum(basal_keep),
      n_cauline_excluded = n_caul_drop,
      n_subfloor_excluded = n_area_drop,
      n_fallback_individuals = n_fb_individuals
    )
  )
}

apply_otu_confidence_routing_tier <- function(specimen_df, closed_leaf_df) {
  if (is.null(specimen_df) || nrow(specimen_df) == 0L || is.null(closed_leaf_df) || nrow(closed_leaf_df) == 0L) return(specimen_df)
  if (!("morph_tier_source" %in% names(specimen_df)) || !("leaf_zone_confidence" %in% names(closed_leaf_df))) return(specimen_df)
  s_keys <- paste0(trimws(as.character(specimen_df$catalogNumber)), "__", ifelse(is.na(specimen_df$plant_individual_id), 0L, as.integer(specimen_df$plant_individual_id)))
  l_keys <- paste0(trimws(as.character(closed_leaf_df$catalogNumber)), "__", ifelse(is.na(closed_leaf_df$plant_individual_id), 0L, as.integer(closed_leaf_df$plant_individual_id)))
  lz_l <- if ("leaf_zone" %in% names(closed_leaf_df)) tolower(trimws(as.character(closed_leaf_df$leaf_zone))) else rep("basal", nrow(closed_leaf_df))
  zc_l <- suppressWarnings(as.numeric(closed_leaf_df$leaf_zone_confidence))
  grp_l <- split(seq_len(nrow(closed_leaf_df)), l_keys)
  n_trans <- 0L
  for (i in seq_len(nrow(specimen_df))) {
    if (isTRUE(specimen_df$morph_tier_source[i] == "Cauline_Only_Fallback")) next
    idx_l <- grp_l[[s_keys[i]]]
    if (is.null(idx_l) || length(idx_l) == 0L) next
    bas_i <- idx_l[is.na(lz_l[idx_l]) | lz_l[idx_l] == "basal"]
    if (length(bas_i) >= 1L && all(is.finite(zc_l[bas_i]) & zc_l[bas_i] < 0.65)) {
      specimen_df$morph_tier_source[i] <- "Transitional_Passive"
      n_trans <- n_trans + 1L
    }
  }
  if (n_trans > 0L) message(sprintf("[Step 03 Confidence Routing] Flagged %d OTUs with zero high-confidence (conf >= 0.65) basal leaves as Transitional_Passive.", n_trans))
  specimen_df
}

.extract_harmonic_opt_from_cal <- function(cal_res, min_h = 6L, max_h = 24L, target_power = 0.99) {
  if (is.null(cal_res) || is.null(cal_res$q) || length(cal_res$q) == 0L) return(NULL)
  q_mat <- as.matrix(cal_res$q)
  q_vec <- if ("50%" %in% rownames(q_mat)) as.numeric(q_mat["50%", ]) else apply(q_mat, 2, stats::median, na.rm = TRUE)
  if (!any(is.finite(q_vec))) return(NULL)
  if (max(q_vec, na.rm = TRUE) <= 1.05) q_vec <- q_vec * 100.0
  h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(q_mat))))
  if (length(h_nums) != length(q_vec) || any(is.na(h_nums))) h_nums <- seq_along(q_vec)
  tp <- if (is.null(target_power) || length(target_power) == 0L) 0.99 else suppressWarnings(as.numeric(target_power[1L]))
  thr_pct <- if (is.finite(tp) && tp <= 1.0) tp * 100.0 else if (is.finite(tp)) tp else 99.0
  hit_idx <- which(is.finite(q_vec) & q_vec >= thr_pct)
  pick_i <- if (length(hit_idx) > 0L) hit_idx[1L] else which.max(h_nums)
  h_val <- max(as.integer(min_h), min(as.integer(max_h), as.integer(h_nums[pick_i])))
  col_m <- match(h_val, h_nums)
  cum_p <- as.numeric(q_vec[if (!is.na(col_m)) col_m else pick_i])
  list(h_opt = as.integer(h_val), cum_pct = cum_p)
}

calibrate_dissection_aware_harmonics <- function(
  out_obj = NULL,
  harm_mat = NULL,
  solidity_vec = NULL,
  nb_h = 12L,
  min_h = 6L,
  max_h = 24L,
  target_power = 0.99
) {
  fb_h <- max(as.integer(min_h), min(as.integer(max_h), as.integer(if (is.null(nb_h) || is.na(nb_h)) 12L else nb_h)))
  h_opt <- fb_h
  cum_pct <- NA_real_

  if (!is.null(out_obj) && requireNamespace("Momocs", quietly = TRUE) && !is.null(out_obj$coo)) {
    n_tot <- length(out_obj$coo)
    sol <- if (!is.null(solidity_vec) && length(solidity_vec) == n_tot) {
      suppressWarnings(as.numeric(solidity_vec))
    } else if (!is.null(out_obj$fac) && "solidity" %in% names(out_obj$fac)) {
      suppressWarnings(as.numeric(out_obj$fac$solidity))
    } else {
      rep(NA_real_, n_tot)
    }

    strat_idx <- if (n_tot > 240L && sum(is.finite(sol)) >= 20L) {
      qs <- stats::quantile(sol[is.finite(sol)], probs = c(0, 0.25, 0.50, 0.75, 1.0), na.rm = TRUE, names = FALSE)
      bins <- findInterval(ifelse(is.finite(sol), sol, stats::median(sol, na.rm = TRUE)), unique(qs), rightmost.closed = TRUE)
      unlist(lapply(sort(unique(bins)), function(b) {
        ib <- which(bins == b)
        if (length(ib) <= 60L) ib else ib[round(seq(1, length(ib), length.out = 60L))]
      }), use.names = FALSE)
    } else if (n_tot > 240L) {
      round(seq(1, n_tot, length.out = 240L))
    } else {
      seq_len(n_tot)
    }
    strat_idx <- sort(unique(pmax(1L, pmin(n_tot, as.integer(strat_idx)))))

    cal_all <- tryCatch({
      sub_all <- if (length(strat_idx) < n_tot) Momocs::slice(out_obj, strat_idx) else out_obj
      Momocs::calibrate_harmonicpower_efourier(sub_all, nb.h = as.integer(max_h), plot = FALSE)
    }, error = function(e) NULL)
    parsed_all <- .extract_harmonic_opt_from_cal(cal_all, min_h = min_h, max_h = max_h, target_power = target_power)
    if (!is.null(parsed_all)) {
      h_opt <- parsed_all$h_opt
      cum_pct <- parsed_all$cum_pct
    }

    # Ensure dissected / lowest-solidity quartile also achieves target cumulative harmonic power
    dis_idx <- which(is.finite(sol) & sol < 0.72)
    if (length(dis_idx) >= 5L) {
      if (length(dis_idx) > 80L) {
        ord_d <- dis_idx[order(sol[dis_idx])]
        dis_idx <- sort(ord_d[round(seq(1, length(ord_d), length.out = 80L))])
      }
      cal_dis <- tryCatch({
        sub_dis <- Momocs::slice(out_obj, dis_idx)
        Momocs::calibrate_harmonicpower_efourier(sub_dis, nb.h = as.integer(max_h), plot = FALSE)
      }, error = function(e) NULL)
      parsed_dis <- .extract_harmonic_opt_from_cal(cal_dis, min_h = min_h, max_h = max_h, target_power = target_power)
      if (!is.null(parsed_dis) && parsed_dis$h_opt > h_opt) {
        h_opt <- parsed_dis$h_opt
        cum_pct <- min(c(cum_pct, parsed_dis$cum_pct), na.rm = TRUE)
      }
    }
  }

  list(h_opt = as.integer(h_opt), cum_pct = as.numeric(cum_pct))
}
