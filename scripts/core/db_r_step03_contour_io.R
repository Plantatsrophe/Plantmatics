if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step03_contour_io.R
# Step 03 manifest reconciliation, parallel contour decoding, and fetch_step03_contours_from_db coordinator.
# ==============================================================================

.decode_step03_contour_rows_r <- function(
  rows,
  vouchers_path = NULL,
  manifest_path = NULL,
  active_catalogs = NULL,
  align_fn = NULL,
  centsize_fn = NULL
) {
  # If explicit manifest_path exists on disk with rows that are NOT yet in SQLite, defer to disk fallback;
  # otherwise build O(1) lookup tables for positional homology & cohort-adaptive blade area columns
  eff_mpath <- if (!is.null(manifest_path) && nzchar(manifest_path) && file.exists(manifest_path)) {
    manifest_path
  } else if (!is.null(vouchers_path) && nzchar(vouchers_path) && file.exists(file.path(dirname(vouchers_path), "extracted_leaf_manifest.csv"))) {
    file.path(dirname(vouchers_path), "extracted_leaf_manifest.csv")
  } else {
    ""
  }
  m_hdr <- NULL
  m_idx_3k <- NULL
  m_idx_2k <- NULL
  if (nzchar(eff_mpath) && file.exists(eff_mpath)) {
    m_hdr <- if (requireNamespace("data.table", quietly = TRUE)) {
      tryCatch(data.table::fread(eff_mpath, showProgress = FALSE, data.table = FALSE), error = function(e) NULL)
    } else {
      tryCatch(as.data.frame(readr::read_csv(eff_mpath, show_col_types = FALSE, progress = FALSE)), error = function(e) NULL)
    }
    if (!is.null(m_hdr) && nrow(m_hdr) > 0L && "catalogNumber" %in% names(m_hdr)) {
      if (!is.null(manifest_path) && nzchar(manifest_path) && identical(eff_mpath, manifest_path)) {
        disk_cats <- unique(trimws(as.character(m_hdr$catalogNumber)))
        if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
          disk_cats <- intersect(disk_cats, trimws(as.character(active_catalogs)))
        }
        if (length(disk_cats) > 0L && length(intersect(disk_cats, unique(rows$catalogNumber))) == 0L) {
          return(NULL)
        }
      }
      if ("plant_individual_id" %in% names(m_hdr)) {
        mc_all_v <- trimws(as.character(m_hdr$catalogNumber))
        mp_all_v <- suppressWarnings(as.integer(m_hdr$plant_individual_id))
        mp_all_v[is.na(mp_all_v)] <- 0L
        m_cat_pids <- unique(paste(mc_all_v, mp_all_v, sep = "__"))
        r_cats_v <- trimws(as.character(rows$catalogNumber))
        r_pids_v <- if ("plant_individual_id" %in% names(rows)) suppressWarnings(as.integer(rows$plant_individual_id)) else rep(0L, nrow(rows))
        r_pids_v[is.na(r_pids_v)] <- 0L
        keep_ramet_m <- !(r_cats_v %in% mc_all_v) | (paste(r_cats_v, r_pids_v, sep = "__") %in% m_cat_pids)
        if (!all(keep_ramet_m)) {
          rows <- rows[keep_ramet_m, , drop = FALSE]
          if (nrow(rows) == 0L) return(NULL)
        }
      }
      if ("leaf_id" %in% names(m_hdr)) {
        mc_v <- trimws(as.character(m_hdr$catalogNumber))
        ml_v <- suppressWarnings(as.integer(m_hdr$leaf_id))
        mp_v <- if ("plant_individual_id" %in% names(m_hdr)) ifelse(is.na(m_hdr$plant_individual_id), 0L, suppressWarnings(as.integer(m_hdr$plant_individual_id))) else rep(0L, nrow(m_hdr))
        mp_v[is.na(mp_v)] <- 0L
        k3 <- paste(mc_v, mp_v, ml_v, sep = "__")
        k2 <- paste(mc_v, ml_v, sep = "__")
        u3 <- !duplicated(k3, fromLast = TRUE)
        m_idx_3k <- which(u3)
        names(m_idx_3k) <- k3[u3]
        u2 <- !duplicated(k2, fromLast = TRUE)
        m_idx_2k <- which(u2)
        names(m_idx_2k) <- k2[u2]
      }
    }
  }

  n_rows <- nrow(rows)
  keep_qc <- evaluate_step03_lm3_efa_qc(rows, allow_multi_component_ect = TRUE) & !is.na(rows$contour_xy_json) & nzchar(as.character(rows$contour_xy_json))
  n_excl_qc <- sum(!keep_qc, na.rm = TRUE)
  if (n_excl_qc > 0L) {
    message(sprintf("Step 03 SQLite EFA QC filter: excluded %d clump/multi-island contours (retained in SQLite extracted_leaves for audit/ECT).", n_excl_qc))
  }
  if (!any(keep_qc)) {
    return(NULL)
  }
  if (!all(keep_qc)) {
    rows <- rows[keep_qc, , drop = FALSE]
    n_rows <- nrow(rows)
  }

  uncal_sources <- c("mp_estimated", "predicted_from_megapixels", "tick_lattice_unboxed", "uncalibrated")
  mh_rows <- rep(NA_integer_, n_rows)
  if (!is.null(m_hdr) && (!is.null(m_idx_3k) || !is.null(m_idx_2k))) {
    rc_v <- trimws(as.character(rows$catalogNumber))
    rp_v <- ifelse(is.na(rows$plant_individual_id), 0L, suppressWarnings(as.integer(rows$plant_individual_id)))
    rp_v[is.na(rp_v)] <- 0L
    rl_v <- suppressWarnings(as.integer(rows$leaf_id))
    if (!is.null(m_idx_3k)) {
      mh_rows <- unname(m_idx_3k[paste(rc_v, rp_v, rl_v, sep = "__")])
    }
    if (any(is.na(mh_rows)) && !is.null(m_idx_2k)) {
      miss_m <- is.na(mh_rows)
      mh_rows[miss_m] <- unname(m_idx_2k[paste(rc_v[miss_m], rl_v[miss_m], sep = "__")])
    }
  }

  ppm_vec <- suppressWarnings(as.numeric(rows$pixels_per_mm))
  raw_ss_col <- if ("scale_source" %in% names(rows)) tolower(trimws(as.character(rows$scale_source))) else rep("", n_rows)
  raw_ss_col[is.na(raw_ss_col) | raw_ss_col %in% c("na", "null")] <- ""
  if (!is.null(m_hdr) && "scale_source" %in% names(m_hdr) && any(!nzchar(raw_ss_col) & !is.na(mh_rows))) {
    fill_ss <- which(!nzchar(raw_ss_col) & !is.na(mh_rows))
    m_ss_raw <- tolower(trimws(as.character(m_hdr$scale_source[mh_rows[fill_ss]])))
    m_ss_raw[is.na(m_ss_raw) | m_ss_raw %in% c("na", "null")] <- ""
    nz_ss <- nzchar(m_ss_raw)
    if (any(nz_ss)) {
      raw_ss_col[fill_ss[nz_ss]] <- m_ss_raw[nz_ss]
    }
  }

  raw_cal_col <- if ("is_calibrated_mm" %in% names(rows)) toupper(trimws(as.character(rows$is_calibrated_mm))) else rep(NA_character_, n_rows)
  raw_cal_col[is.na(raw_cal_col) | !nzchar(raw_cal_col) | raw_cal_col %in% c("NA", "NULL")] <- NA_character_
  if (!is.null(m_hdr) && "is_calibrated_mm" %in% names(m_hdr) && any(is.na(raw_cal_col) & !is.na(mh_rows))) {
    fill_cal <- which(is.na(raw_cal_col) & !is.na(mh_rows))
    m_cal_raw <- toupper(trimws(as.character(m_hdr$is_calibrated_mm[mh_rows[fill_cal]])))
    m_cal_raw[is.na(m_cal_raw) | !nzchar(m_cal_raw) | m_cal_raw %in% c("NA", "NULL")] <- NA_character_
    ok_cal <- !is.na(m_cal_raw)
    if (any(ok_cal)) {
      raw_cal_col[fill_cal[ok_cal]] <- m_cal_raw[ok_cal]
    }
  }

  cal_explicit_false <- !is.na(raw_cal_col) & (raw_cal_col %in% c("0", "0.0", "FALSE", "F", "NO"))
  cal_explicit_true <- !is.na(raw_cal_col) & (raw_cal_col %in% c("1", "1.0", "TRUE", "T", "YES"))
  cal_flag_vec <- ifelse(
    raw_ss_col %in% uncal_sources | cal_explicit_false,
    FALSE,
    ifelse(is.na(raw_cal_col), TRUE, cal_explicit_true)
  )
  is_cal_vec <- !is.na(ppm_vec) & is.finite(ppm_vec) & (ppm_vec > 0) & cal_flag_vec & !(raw_ss_col %in% uncal_sources)
  db_ar_vec <- if ("aspect_ratio" %in% names(rows)) suppressWarnings(as.numeric(rows$aspect_ratio)) else rep(NA_real_, n_rows)
  db_upright_vec <- rep(NA, n_rows)
  if ("mask_is_upright" %in% names(rows)) {
    raw_up <- toupper(trimws(as.character(rows$mask_is_upright)))
    db_upright_vec[raw_up %in% c("1", "1.0", "TRUE", "T", "YES")] <- TRUE
    db_upright_vec[raw_up %in% c("0", "0.0", "FALSE", "F", "NO")] <- FALSE
  }
  if (!is.null(m_hdr) && "mask_is_upright" %in% names(m_hdr) && any(is.na(db_upright_vec))) {
    fill_u <- which(is.na(db_upright_vec) & !is.na(mh_rows))
    if (length(fill_u) > 0L) {
      m_up_raw <- toupper(trimws(as.character(m_hdr$mask_is_upright[mh_rows[fill_u]])))
      db_upright_vec[fill_u[m_up_raw %in% c("1", "1.0", "TRUE", "T", "YES")]] <- TRUE
      db_upright_vec[fill_u[m_up_raw %in% c("0", "0.0", "FALSE", "F", "NO")]] <- FALSE
    }
  }
  db_area_vec <- suppressWarnings(as.numeric(rows$area_px))
  db_sol_vec <- suppressWarnings(as.numeric(rows$solidity))
  js_vec <- as.character(rows$contour_xy_json)
  has_momocs <- requireNamespace("Momocs", quietly = TRUE)

  process_contour_chunk <- function(idx_vec) {
    m_len <- length(idx_vec)
    c_list <- vector("list", m_len)
    v_ok <- logical(m_len)
    v_sol <- rep(NA_real_, m_len)
    v_di <- rep(NA_real_, m_len)
    v_cs <- rep(NA_real_, m_len)
    v_ar <- rep(NA_real_, m_len)
    v_area <- rep(NA_real_, m_len)
    v_circ <- rep(NA_real_, m_len)

    for (k in seq_len(m_len)) {
      i <- idx_vec[k]
      coo <- tryCatch(jsonlite::fromJSON(js_vec[i], simplifyMatrix = TRUE), error = function(e) NULL)
      if (is.null(coo) || !is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 10L) next
      coo <- matrix(as.numeric(coo[, 1:2]), ncol = 2L)
      if (any(!is.finite(coo))) next
      colnames(coo) <- c("x", "y")
      was_tipup_db <- FALSE
      if (exists("align_contour_horizontal_xax", mode = "function")) {
        coo <- align_contour_horizontal_xax(coo, is_upright = db_upright_vec[i])
        was_tipup_db <- isTRUE(attr(coo, "lm3_from_tipup"))
      }
      if (exists("ensure_clockwise_keep_start_r", mode = "function")) {
        coo <- ensure_clockwise_keep_start_r(coo)
      } else if (nrow(coo) >= 3L) {
        sa2 <- sum(coo[, 1L] * c(coo[-1L, 2L], coo[1L, 2L]) - c(coo[-1L, 1L], coo[1L, 1L]) * coo[, 2L])
        if (is.finite(sa2) && sa2 > 0.0) {
          coo <- coo[c(1L, seq.int(nrow(coo), 2L, by = -1L)), , drop = FALSE]
        }
      }

      x_next <- c(coo[-1L, 1L], coo[1L, 1L])
      y_next <- c(coo[-1L, 2L], coo[1L, 2L])
      area_px <- 0.5 * abs(sum(coo[, 1L] * y_next - x_next * coo[, 2L]))
      perim_px <- sum(sqrt((x_next - coo[, 1L])^2 + (y_next - coo[, 2L])^2))
      circ_val <- if (is.finite(perim_px) && perim_px > 1e-6) round(min(1.0, (4 * pi * area_px) / (perim_px^2)), 4) else NA_real_
      di_val <- if (is.finite(perim_px) && is.finite(area_px) && area_px > 1e-6) round(max(1.0, perim_px / (2.0 * sqrt(pi * area_px))), 4) else NA_real_

      if (is_cal_vec[i]) {
        coo <- coo / ppm_vec[i]
        cs_val <- if (!is.null(centsize_fn)) {
          centsize_fn(coo)
        } else {
          cx <- mean(coo[, 1L])
          cy <- mean(coo[, 2L])
          sqrt(sum((coo[, 1L] - cx)^2 + (coo[, 2L] - cy)^2))
        }
      } else {
        cs_val <- NA_real_
      }

      if (!is.null(align_fn)) {
        coo <- align_fn(coo, circularity = circ_val, has_momocs = has_momocs)
      }
      attr(coo, "lm3_xax_aligned") <- TRUE
      attr(coo, "lm3_from_tipup") <- was_tipup_db

      db_ar_i <- db_ar_vec[i]
      if (is.finite(db_ar_i) && db_ar_i > 0.0) {
        ar_val <- round(as.numeric(db_ar_i), 4)
      } else if (was_tipup_db && diff(range(coo[, 2L], na.rm = TRUE)) > 1e-6) {
        ar_val <- round(as.numeric(diff(range(coo[, 1L], na.rm = TRUE)) / diff(range(coo[, 2L], na.rm = TRUE))), 4)
      } else {
        lw <- if (has_momocs) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
        if (is.null(lw) || length(lw) < 2L || is.na(lw[2]) || lw[2] <= 1e-6) {
          coo_centered <- scale(coo, center = TRUE, scale = FALSE)
          svd_rot <- coo_centered %*% svd(coo_centered)$v
          lw <- c(diff(range(svd_rot[, 1], na.rm = TRUE)), diff(range(svd_rot[, 2], na.rm = TRUE)))
        }
        ar_val <- round(as.numeric(lw[1] / lw[2]), 4)
      }

      db_a <- db_area_vec[i]
      eff_area <- if (!is.na(db_a) && is.finite(db_a) && db_a > 0) db_a else round(area_px, 1)
      sol_val <- db_sol_vec[i]
      if (is.na(sol_val) || !is.finite(sol_val)) {
        ch_idx <- tryCatch(grDevices::chull(coo[, 1L], coo[, 2L]), error = function(e) integer(0))
        if (length(ch_idx) >= 3L) {
          hc <- coo[ch_idx, , drop = FALSE]
          h_area <- 0.5 * abs(sum(hc[, 1L] * c(hc[-1L, 2L], hc[1L, 2L]) - c(hc[-1L, 1L], hc[1L, 1L]) * hc[, 2L]))
          if (h_area > 1e-6) sol_val <- round(min(1.0, area_px / h_area), 4)
        }
      }

      v_ok[k] <- TRUE
      c_list[[k]] <- coo
      v_sol[k] <- sol_val
      v_di[k] <- di_val
      v_cs[k] <- as.numeric(cs_val)
      v_ar[k] <- as.numeric(ar_val)
      v_area[k] <- as.numeric(eff_area)
      v_circ[k] <- as.numeric(circ_val)
    }

    list(
      ok = v_ok, coo = c_list, sol = v_sol, di = v_di,
      cs = v_cs, ar = v_ar, area = v_area, circ = v_circ
    )
  }

  n_workers <- if (n_rows >= 160L) resolve_packera_workers(32L, ceiling(n_rows / 80L)) else 1L
  if (n_workers <= 1L) {
    merged <- process_contour_chunk(seq_len(n_rows))
  } else {
    chunk_ids <- split(seq_len(n_rows), cut(seq_len(n_rows), breaks = n_workers, labels = FALSE))
    chunk_out <- packera_mclapply(chunk_ids, process_contour_chunk, mc.cores = n_workers)
    merged <- list(
      ok = unlist(lapply(chunk_out, `[[`, "ok"), use.names = FALSE),
      coo = unlist(lapply(chunk_out, `[[`, "coo"), recursive = FALSE, use.names = FALSE),
      sol = unlist(lapply(chunk_out, `[[`, "sol"), use.names = FALSE),
      di = unlist(lapply(chunk_out, `[[`, "di"), use.names = FALSE),
      cs = unlist(lapply(chunk_out, `[[`, "cs"), use.names = FALSE),
      ar = unlist(lapply(chunk_out, `[[`, "ar"), use.names = FALSE),
      area = unlist(lapply(chunk_out, `[[`, "area"), use.names = FALSE),
      circ = unlist(lapply(chunk_out, `[[`, "circ"), use.names = FALSE)
    )
  }

  ok_idx <- which(merged$ok)
  if (length(ok_idx) == 0L) {
    return(NULL)
  }
  list(
    rows = rows,
    ok_idx = ok_idx,
    merged = merged,
    m_hdr = m_hdr,
    mh_rows = mh_rows,
    is_cal_vec = is_cal_vec,
    ppm_vec = ppm_vec,
    raw_ss_col = raw_ss_col,
    db_upright_vec = db_upright_vec
  )
}

fetch_step03_contours_from_db <- function(
  vouchers_path = NULL,
  manifest_path = NULL,
  contours_dir = NULL,
  active_catalogs = NULL,
  align_fn = NULL,
  centsize_fn = NULL,
  run_id = NULL
) {
  rows <- .query_step03_contour_rows_db(
    vouchers_path = vouchers_path,
    manifest_path = manifest_path,
    contours_dir = contours_dir,
    active_catalogs = active_catalogs,
    run_id = run_id
  )
  if (is.null(rows) || nrow(rows) == 0L) {
    return(NULL)
  }
  decoded <- .decode_step03_contour_rows_r(
    rows = rows,
    vouchers_path = vouchers_path,
    manifest_path = manifest_path,
    active_catalogs = active_catalogs,
    align_fn = align_fn,
    centsize_fn = centsize_fn
  )
  if (is.null(decoded)) {
    return(NULL)
  }
  .assemble_step03_contour_meta_df_r(decoded)
}
