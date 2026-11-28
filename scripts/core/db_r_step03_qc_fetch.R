if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step03_qc_fetch.R
# Step 03 LM3 EFA QC evaluator and SQLite contour row query helper.
# ==============================================================================

# Reproducible Step 03 EFA inclusion rule on read-only LeafMachine3 SQLite / manifest records:
#   - Respects explicit `efa_qc_pass` when present, and evaluates LeafMachine3's native `bilateral_dice`
#     & `si_a` with relaxed bounds on genuinely dissected Packera leaves (`dice >= 0.45`, `|si_a| < 0.48`)
#     while excluding fused overlapping clumps (`UNRESOLVED_CLUMP`) and unbridged multi-island pinnae
#     (`n_lamina_components > 1`, which remain in SQLite for multi-component ECT).
evaluate_step03_lm3_efa_qc <- function(df, allow_multi_component_ect = FALSE) {
  if (is.null(df) || nrow(df) == 0L) {
    return(logical(0))
  }
  n <- nrow(df)
  keep <- rep(TRUE, n)

  n_cc_vec <- if ("n_lamina_components" %in% names(df)) suppressWarnings(as.integer(df$n_lamina_components)) else rep(1L, n)
  is_multi_cc <- !is.na(n_cc_vec) & (n_cc_vec > 1L)
  evf_vec <- if ("enclosed_void_frac" %in% names(df)) suppressWarnings(as.numeric(df$enclosed_void_frac)) else rep(0.0, n)
  is_enc_void <- !is.na(evf_vec) & is.finite(evf_vec) & (evf_vec > 0.020)

  ect_only_codes <- c("MULTI_COMPONENT_PINNAE", "ENCLOSED_INTERLOBE_VOID")
  raw_qfr <- if ("qc_flag_reason" %in% names(df)) toupper(trimws(as.character(df$qc_flag_reason))) else rep("", n)
  raw_qfr[is.na(raw_qfr)] <- ""
  is_ect_only_flag <- vapply(raw_qfr, function(s) {
    if (!nzchar(s) || s %in% c("PASS", "OK", "NONE", "NA")) {
      return(FALSE)
    }
    toks <- trimws(unlist(strsplit(s, "[;,]")))
    toks <- toks[nzchar(toks) & !(toks %in% c("PASS", "OK", "NONE", "NA"))]
    length(toks) > 0L && all(toks %in% ect_only_codes)
  }, logical(1))

  if ("efa_qc_pass" %in% names(df)) {
    raw_eqc <- tolower(trimws(as.character(df$efa_qc_pass)))
    eqc_fail <- !is.na(raw_eqc) & raw_eqc %in% c("0", "0.0", "false", "f", "no")
    if (isTRUE(allow_multi_component_ect)) {
      keep[eqc_fail & !(is_ect_only_flag | ((is_multi_cc | is_enc_void) & !nzchar(raw_qfr)))] <- FALSE
    } else {
      keep[eqc_fail] <- FALSE
    }
  }
  if ("qc_flag_reason" %in% names(df)) {
    bad_qfr <- nzchar(raw_qfr) & !(raw_qfr %in% c("PASS", "OK", "NONE", "NA"))
    if (isTRUE(allow_multi_component_ect)) {
      keep[bad_qfr & !is_ect_only_flag] <- FALSE
    } else {
      keep[bad_qfr] <- FALSE
    }
  }
  if ("detector_class_id" %in% names(df)) {
    dcid <- suppressWarnings(as.integer(df$detector_class_id))
    has_exp_pass <- if ("efa_qc_pass" %in% names(df)) {
      tolower(trimws(as.character(df$efa_qc_pass))) %in% c("1", "1.0", "true", "t", "yes")
    } else {
      rep(FALSE, n)
    }
    keep[!is.na(dcid) & (dcid == 1L) & !has_exp_pass] <- FALSE
  }
  if ("is_crop_clipped" %in% names(df)) {
    raw_cc <- tolower(trimws(as.character(df$is_crop_clipped)))
    keep[!is.na(raw_cc) & raw_cc %in% c("1", "1.0", "true", "t", "yes")] <- FALSE
  }
  if ("hole_area_frac" %in% names(df)) {
    haf <- suppressWarnings(as.numeric(df$hole_area_frac))
    keep[!is.na(haf) & is.finite(haf) & (haf > 0.035)] <- FALSE
  }

  sol <- if ("solidity" %in% names(df)) suppressWarnings(as.numeric(df$solidity)) else rep(NA_real_, n)
  sia_abs <- if ("si_a" %in% names(df)) abs(suppressWarnings(as.numeric(df$si_a))) else rep(NA_real_, n)
  bls_vec <- if ("bilateral_lobe_symmetry" %in% names(df)) {
    raw_bls <- tolower(trimws(as.character(df$bilateral_lobe_symmetry)))
    ifelse(is.na(raw_bls) | !nzchar(raw_bls) | raw_bls == "na", NA, raw_bls %in% c("1", "1.0", "true", "t", "yes"))
  } else {
    rep(NA, n)
  }
  # Grant alternate-pinna bilateral_dice relaxation ONLY when bilateral_lobe_symmetry is verified
  is_verified_diss <- ((!is.na(sol) & is.finite(sol) & (sol < 0.72)) &
    ifelse(!is.na(bls_vec), bls_vec, is.na(sia_abs) | (sia_abs < 0.28))) |
    (is_multi_cc & ifelse(!is.na(bls_vec), bls_vec, is.na(sia_abs) | (sia_abs < 0.28)))

  if ("lm3_gates_pass" %in% names(df)) {
    raw_gp <- tolower(trimws(as.character(df$lm3_gates_pass)))
    gp_fail <- !is.na(raw_gp) & raw_gp %in% c("0", "0.0", "false", "f", "no")
    keep[!is_verified_diss & gp_fail] <- FALSE
  }
  if ("bilateral_dice" %in% names(df)) {
    bd <- suppressWarnings(as.numeric(df$bilateral_dice))
    min_bd <- ifelse(is_verified_diss, 0.58, 0.70)
    keep[!is.na(bd) & is.finite(bd) & (bd < min_bd)] <- FALSE
  }
  if ("si_a" %in% names(df)) {
    max_sia <- ifelse(!is.na(sol) & is.finite(sol) & (sol < 0.72), 0.28, 0.32)
    keep[!is.na(sia_abs) & is.finite(sia_abs) & (sia_abs >= max_sia)] <- FALSE
  }
  if (!isTRUE(allow_multi_component_ect)) {
    if (any(is_multi_cc)) keep[is_multi_cc] <- FALSE
    if (any(is_enc_void)) keep[is_enc_void] <- FALSE
  }
  keep
}

.query_step03_contour_rows_db <- function(
  vouchers_path = NULL,
  manifest_path = NULL,
  contours_dir = NULL,
  active_catalogs = NULL,
  run_id = NULL
) {
  if (!is_sqlite_active_for_inputs(c(vouchers_path, manifest_path, contours_dir))) {
    return(NULL)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(NULL)
  }

  conn <- connect_packera_db(read_only = TRUE)
  if (is.null(conn)) {
    return(NULL)
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "extracted_leaves")) {
    return(NULL)
  }

  rid <- resolve_packera_run_id(run_id)
  has_cohort <- FALSE
  if (db_table_exists_safe(conn, "run_active_cohort")) {
    cnt_df <- tryCatch(
      DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM run_active_cohort WHERE run_id = $1;", params = list(rid)),
      error = function(e) data.frame(n = 0L)
    )
    has_cohort <- isTRUE(cnt_df$n[1] > 0L)
  }

  el_fields <- tryCatch(DBI::dbListFields(conn, "extracted_leaves"), error = function(e) character(0))
  has_pi_tbl <- db_table_exists_safe(conn, "plant_individuals")
  has_vrp_tbl <- db_table_exists_safe(conn, "voucher_reproductive_phenology")
  has_v_tbl <- db_table_exists_safe(conn, "vouchers")
  pi_fields <- if (has_pi_tbl) tryCatch(DBI::dbListFields(conn, "plant_individuals"), error = function(e) character(0)) else character(0)
  vrp_fields <- if (has_vrp_tbl) tryCatch(DBI::dbListFields(conn, "voucher_reproductive_phenology"), error = function(e) character(0)) else character(0)
  v_fields <- if (has_v_tbl) tryCatch(DBI::dbListFields(conn, "vouchers"), error = function(e) character(0)) else character(0)

  sel_opt <- function(alias, fields, col, fallback = "NULL") {
    if (col %in% fields) sprintf("%s.%s", alias, col) else sprintf("%s AS %s", fallback, col)
  }
  sel_coal2 <- function(col, a1, f1, a2, f2, fallback = "NULL") {
    in1 <- col %in% f1
    in2 <- col %in% f2
    if (in1 && in2) {
      sprintf("COALESCE(%s.%s, %s.%s) AS %s", a1, col, a2, col, col)
    } else if (in1) {
      sprintf("%s.%s AS %s", a1, col, col)
    } else if (in2) {
      sprintf("%s.%s AS %s", a2, col, col)
    } else {
      sprintf("%s AS %s", fallback, col)
    }
  }
  sel_vrp <- function(col) {
    if (col %in% vrp_fields) {
      sprintf("CASE WHEN vrp.catalogNumber IS NOT NULL THEN vrp.%s ELSE vrp0.%s END AS %s", col, col, col)
    } else {
      sprintf("NULL AS %s", col)
    }
  }

  sel_ect <- sel_opt("el", el_fields, "ect_path")
  sel_zone <- sel_opt("el", el_fields, "leaf_zone", "'basal'")
  sel_zone_conf <- sel_opt("el", el_fields, "leaf_zone_confidence")
  sel_dist_cr <- sel_opt("el", el_fields, "dist_from_crown_norm")
  sel_cy_rel <- sel_opt("el", el_fields, "cy_rel_plant_bbox")
  sel_arch <- sel_opt("el", el_fields, "archetype_score")
  sel_bdice <- sel_opt("el", el_fields, "bilateral_dice")
  sel_sia <- sel_opt("el", el_fields, "si_a")
  sel_lgates <- sel_opt("el", el_fields, "lm3_gates_pass")
  sel_ncc <- sel_opt("el", el_fields, "n_lamina_components", "1")
  sel_scale_src <- sel_opt("el", el_fields, "scale_source")
  sel_upright <- sel_opt("el", el_fields, "mask_is_upright")
  sel_amm2 <- sel_opt("el", el_fields, "area_mm2")
  sel_bl_lnorm <- sel_opt("el", el_fields, "blade_cie_l_norm")
  sel_bl_lsd <- sel_opt("el", el_fields, "blade_cie_l_sd")
  sel_bl_ldelta <- sel_opt("el", el_fields, "blade_cie_l_sheet_delta")
  sel_bl_adelta <- sel_opt("el", el_fields, "blade_cie_a_delta")
  sel_bl_bdelta <- sel_opt("el", el_fields, "blade_cie_b_delta")
  sel_bw_mm <- sel_opt("el", el_fields, "blade_width_mm")
  sel_pw_mm <- sel_coal2("petiole_width_mm", "el", el_fields, "pi", pi_fields)
  sel_dcid <- sel_opt("el", el_fields, "detector_class_id", "0")
  sel_cclip <- sel_opt("el", el_fields, "is_crop_clipped", "0")
  sel_haf <- sel_opt("el", el_fields, "hole_area_frac", "0.0")
  sel_evf <- sel_opt("el", el_fields, "enclosed_void_frac", "0.0")
  sel_bls <- sel_opt("el", el_fields, "bilateral_lobe_symmetry")
  sel_eqc <- sel_opt("el", el_fields, "efa_qc_pass")
  sel_qfr <- sel_opt("el", el_fields, "qc_flag_reason")

  sel_el_step02_extra <- paste(
    c(
      sel_opt("el", el_fields, "bbox_x1"),
      sel_opt("el", el_fields, "bbox_y1"),
      sel_opt("el", el_fields, "bbox_x2"),
      sel_opt("el", el_fields, "bbox_y2"),
      sel_opt("el", el_fields, "base_angle_deg"),
      sel_opt("el", el_fields, "apex_angle_deg"),
      sel_opt("el", el_fields, "lamina_curvature"),
      sel_opt("el", el_fields, "lamina_trace_length_px"),
      sel_opt("el", el_fields, "basal_cordate_index"),
      sel_opt("el", el_fields, "margin_serration_index"),
      sel_opt("el", el_fields, "margin_tooth_acuteness"),
      sel_opt("el", el_fields, "terminal_lobe_area_frac"),
      sel_opt("el", el_fields, "blade_cie_l_q85_sheet_delta"),
      sel_opt("el", el_fields, "basal_zone_cie_l_delta"),
      sel_opt("el", el_fields, "is_abaxial_pressed_leaf", "0"),
      sel_opt("el", el_fields, "abaxial_cie_l_sheet_delta")
    ),
    collapse = ", "
  )

  sel_pi_extra <- paste(
    c(
      sel_opt("pi", pi_fields, "rosette_centroid_x"),
      sel_opt("pi", pi_fields, "rosette_centroid_y")
    ),
    collapse = ", "
  )

  sel_vrp_all <- paste(
    vapply(
      c(
        "AR_inv", "capitula_count", "involucre_height_px", "involucre_width_px",
        "involucre_height_mm", "involucre_width_mm",
        "corymb_head_count_est", "radiate_head_fraction",
        "phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta",
        "involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta", "calyculus_flare_index"
      ),
      sel_vrp,
      character(1L)
    ),
    collapse = ", "
  )

  sel_v_extra <- paste(
    c(
      sel_opt("v", v_fields, "institutionCode"),
      sel_opt("v", v_fields, "scientificName"),
      sel_opt("v", v_fields, "species_raw"),
      sel_opt("v", v_fields, "species_standardized"),
      sel_opt("v", v_fields, "recordedBy"),
      sel_opt("v", v_fields, "recordNumber"),
      sel_opt("v", v_fields, "identifiedBy"),
      sel_opt("v", v_fields, "determiner_tier"),
      sel_opt("v", v_fields, "type_status"),
      sel_opt("v", v_fields, "county"),
      sel_opt("v", v_fields, "stateProvince"),
      sel_opt("v", v_fields, "locality"),
      sel_opt("v", v_fields, "decimalLatitude"),
      sel_opt("v", v_fields, "decimalLongitude"),
      sel_opt("v", v_fields, "coordinateUncertainty"),
      sel_opt("v", v_fields, "year"),
      sel_opt("v", v_fields, "month"),
      sel_opt("v", v_fields, "day"),
      sel_opt("v", v_fields, "eventDate"),
      sel_opt("v", v_fields, "doy"),
      sel_opt("v", v_fields, "regional_group"),
      sel_opt("v", v_fields, "infraspecific_taxon"),
      sel_opt("v", v_fields, "is_primary_duplicate"),
      sel_opt("v", v_fields, "duplicate_count"),
      sel_opt("v", v_fields, "spatial_balance_weight"),
      sel_opt("v", v_fields, "typified_basionym"),
      sel_opt("v", v_fields, "verbatim_taxon_concept"),
      sel_opt("v", v_fields, "concept_synonymy_status"),
      sel_opt("v", v_fields, "identification_year"),
      sel_opt("v", v_fields, "is_open_or_hybrid_concept", "0"),
      sel_opt("v", v_fields, "georef_qc_status"),
      sel_opt("v", v_fields, "is_centroid_or_inst_coord", "0"),
      sel_opt("v", v_fields, "coord_quarantine_reason"),
      sel_opt("v", v_fields, "usda_symbol"),
      sel_opt("v", v_fields, "usda_county_status"),
      sel_opt("v", v_fields, "exsiccatae_key"),
      sel_opt("v", v_fields, "exsiccatae_taxonomic_conflict", "0"),
      sel_opt("v", v_fields, "exsiccatae_gold_taxon")
    ),
    collapse = ", "
  )

  v_join_sql <- if (has_v_tbl) "LEFT JOIN vouchers v ON v.catalogNumber = el.catalogNumber " else ""
  pi_join_sql <- if (has_pi_tbl) "LEFT JOIN plant_individuals pi ON pi.catalogNumber = el.catalogNumber AND pi.plant_individual_id = el.plant_individual_id " else ""
  vrp_join_sql <- if (has_vrp_tbl) {
    paste0(
      "LEFT JOIN voucher_reproductive_phenology vrp ON vrp.catalogNumber = el.catalogNumber AND vrp.plant_individual_id = el.plant_individual_id ",
      "LEFT JOIN voucher_reproductive_phenology vrp0 ON vrp0.catalogNumber = el.catalogNumber AND vrp0.plant_individual_id = 0 "
    )
  } else ""

  base_sql <- paste0(
    "SELECT ",
    "  el.catalogNumber, el.plant_individual_id, el.leaf_id, el.shape_id, ",
    "  el.solidity, el.aspect_ratio, el.area_px, ",
    "  ", sel_ect, ", ", sel_zone, ", ", sel_zone_conf, ", ", sel_dist_cr, ", ", sel_cy_rel, ", ",
    "  ", sel_arch, ", ", sel_bdice, ", ", sel_sia, ", ", sel_lgates, ", ", sel_ncc, ", ", sel_scale_src, ", ",
    "  ", sel_upright, ", ", sel_amm2, ", ", sel_bl_lnorm, ", ", sel_bl_lsd, ", ",
    sel_bl_ldelta, ", ", sel_bl_adelta, ", ", sel_bl_bdelta, ", ", sel_bw_mm, ", ", sel_pw_mm, ", ",
    "  ", sel_dcid, ", ", sel_cclip, ", ", sel_haf, ", ", sel_evf, ", ", sel_bls, ", ", sel_eqc, ", ", sel_qfr, ", ",
    "  ", sel_el_step02_extra, ", ",
    "  ", sel_coal2("pixels_per_mm", "el", el_fields, "v", v_fields), ", ",
    "  ", sel_coal2("is_calibrated_mm", "el", el_fields, "v", v_fields), ", ",
    "  ", sel_coal2("petiole_length_mm", "el", el_fields, "pi", pi_fields), ", ",
    "  ", sel_coal2("blade_length_mm", "el", el_fields, "pi", pi_fields), ", ",
    "  ", sel_coal2("petiole_to_blade_ratio", "el", el_fields, "pi", pi_fields), ", ",
    "  ", sel_pi_extra, ", ",
    "  ", sel_vrp_all, ", ",
    "  ", sel_v_extra, ", ",
    "  el.contour_xy_json ",
    "FROM extracted_leaves el ",
    v_join_sql,
    pi_join_sql,
    vrp_join_sql
  )

  contour_pred <- "el.contour_xy_json IS NOT NULL AND el.contour_xy_json != ''"

  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) == 0L) {
      return(NULL)
    }
    q_cats <- paste(vapply(cats, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
    rac_cats <- if (has_cohort) {
      tryCatch(trimws(as.character(DBI::dbGetQuery(conn, "SELECT DISTINCT catalogNumber FROM run_active_cohort WHERE run_id = $1;", params = list(rid))$catalogNumber)), error = function(e) character(0))
    } else character(0)
    if (has_cohort && length(rac_cats) > 0L && all(cats %in% rac_cats)) {
      sql <- paste0(
        base_sql,
        "INNER JOIN run_active_cohort rac ON rac.catalogNumber = el.catalogNumber AND rac.plant_individual_id = el.plant_individual_id ",
        sprintf("WHERE rac.run_id = $1 AND %s AND el.catalogNumber IN (%s) ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;", contour_pred, q_cats)
      )
      rows <- tryCatch(DBI::dbGetQuery(conn, sql, params = list(rid)), error = function(e) NULL)
    } else {
      sql <- paste0(
        base_sql,
        sprintf("WHERE %s AND el.catalogNumber IN (%s) ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;", contour_pred, q_cats)
      )
      rows <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
    }
  } else if (has_cohort) {
    sql <- paste0(
      base_sql,
      "INNER JOIN run_active_cohort rac ON rac.catalogNumber = el.catalogNumber AND rac.plant_individual_id = el.plant_individual_id ",
      "WHERE rac.run_id = $1 AND ", contour_pred, " ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;"
    )
    rows <- tryCatch(DBI::dbGetQuery(conn, sql, params = list(rid)), error = function(e) NULL)
  } else {
    sql <- paste0(base_sql, "WHERE ", contour_pred, " ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;")
    rows <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  }

  # Close DBI handle immediately before any parallel fork
  disconnect_packera_db(conn)
  conn <- NULL

  if (is.null(rows) || nrow(rows) == 0L) {
    return(NULL)
  }
  rows
}
