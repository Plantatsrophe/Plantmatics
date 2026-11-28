if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step03_persist.R
# Step 03 contour metadata assembly, dynamic morphometric normalization, stale-row pruning, and SQLite persistence.
# ==============================================================================

.assemble_step03_contour_meta_df_r <- function(decoded) {
  rows <- decoded$rows
  ok_idx <- decoded$ok_idx
  merged <- decoded$merged
  m_hdr <- decoded$m_hdr
  mh_rows <- decoded$mh_rows
  is_cal_vec <- decoded$is_cal_vec
  ppm_vec <- decoded$ppm_vec
  raw_ss_col <- decoded$raw_ss_col
  db_upright_vec <- decoded$db_upright_vec
  sub_rows <- rows[ok_idx, , drop = FALSE]
  n_kept <- length(ok_idx)
  cat_vec <- as.character(sub_rows$catalogNumber)
  pid_vec <- ifelse(is.na(sub_rows$plant_individual_id), 0L, as.integer(sub_rows$plant_individual_id))
  lid_vec <- ifelse(is.na(sub_rows$leaf_id), ok_idx, as.integer(sub_rows$leaf_id))
  raw_sid_vec <- as.character(sub_rows$shape_id)
  sid_vec <- ifelse(!is.na(raw_sid_vec) & nzchar(raw_sid_vec), raw_sid_vec, sprintf("%s_p%d_leaf%d", cat_vec, pid_vec, lid_vec))
  if (anyDuplicated(sid_vec) > 0L) {
    seen_sids <- new.env(hash = TRUE, parent = emptyenv())
    for (k in seq_len(n_kept)) {
      s_cur <- sid_vec[k]
      if (exists(s_cur, envir = seen_sids, inherits = FALSE)) {
        s_cur <- sprintf("%s_%d", s_cur, ok_idx[k])
        sid_vec[k] <- s_cur
      }
      assign(s_cur, TRUE, envir = seen_sids)
    }
  }

  m_hit_vec <- mh_rows[ok_idx]
  pull_manifest_num <- function(cn, fb_vec = rep(NA_real_, n_kept)) {
    out_v <- fb_vec
    if (!is.null(m_hdr) && cn %in% names(m_hdr) && any(!is.na(m_hit_vec))) {
      hit_mask <- !is.na(m_hit_vec)
      cand_v <- suppressWarnings(as.numeric(m_hdr[[cn]][m_hit_vec[hit_mask]]))
      fin_m <- is.finite(cand_v)
      if (any(fin_m)) {
        hit_pos <- which(hit_mask)
        out_v[hit_pos[fin_m]] <- cand_v[fin_m]
      }
    }
    out_v
  }

  db_cy_rel_vec <- if ("cy_rel_plant_bbox" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$cy_rel_plant_bbox)) else rep(NA_real_, n_kept)
  db_cy_rel_vec[!is.finite(db_cy_rel_vec)] <- NA_real_
  db_zone_conf_vec <- if ("leaf_zone_confidence" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$leaf_zone_confidence)) else rep(NA_real_, n_kept)
  db_zone_conf_vec[!is.finite(db_zone_conf_vec)] <- NA_real_
  db_dist_cr_vec <- if ("dist_from_crown_norm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$dist_from_crown_norm)) else rep(NA_real_, n_kept)
  db_dist_cr_vec[!is.finite(db_dist_cr_vec)] <- NA_real_
  db_zone_vec <- if ("leaf_zone" %in% names(sub_rows)) trimws(as.character(sub_rows$leaf_zone)) else rep("basal", n_kept)
  db_zone_vec[is.na(db_zone_vec) | !nzchar(db_zone_vec)] <- "basal"
  if (!is.null(m_hdr) && "leaf_zone" %in% names(m_hdr) && any(!is.na(m_hit_vec))) {
    hit_m_z <- !is.na(m_hit_vec)
    cand_z <- trimws(as.character(m_hdr[["leaf_zone"]][m_hit_vec[hit_m_z]]))
    valid_z <- !is.na(cand_z) & nzchar(cand_z)
    if (any(valid_z)) {
      hit_pos_z <- which(hit_m_z)
      db_zone_vec[hit_pos_z[valid_z]] <- cand_z[valid_z]
    }
  }
  db_ect_vec <- if ("ect_path" %in% names(sub_rows)) trimws(as.character(sub_rows$ect_path)) else rep(NA_character_, n_kept)
  db_ect_vec[is.na(db_ect_vec) | !nzchar(db_ect_vec)] <- NA_character_
  if (!is.null(m_hdr) && "ect_path" %in% names(m_hdr) && any(!is.na(m_hit_vec))) {
    hit_m_e <- !is.na(m_hit_vec) & is.na(db_ect_vec)
    if (any(hit_m_e)) {
      cand_e <- trimws(as.character(m_hdr[["ect_path"]][m_hit_vec[hit_m_e]]))
      valid_e <- !is.na(cand_e) & nzchar(cand_e)
      if (any(valid_e)) {
        hit_pos_e <- which(hit_m_e)
        db_ect_vec[hit_pos_e[valid_e]] <- cand_e[valid_e]
      }
    }
  }
  sol_kept <- merged$sol[ok_idx]
  eff_diss_vec <- is.finite(sol_kept) & (sol_kept < 0.72)
  cs_kept <- merged$cs[ok_idx]
  log_cs_kept <- ifelse(!is.na(cs_kept) & cs_kept > 0, round(log10(cs_kept), 6), NA_real_)
  arch_kept <- if ("archetype_score" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$archetype_score)) else rep(NA_real_, n_kept)

  cal_kept <- is_cal_vec[ok_idx]
  ppm_kept <- ppm_vec[ok_idx]
  ss_kept <- raw_ss_col[ok_idx]
  resolved_ss_kept <- ifelse(
    cal_kept,
    ifelse(nzchar(ss_kept), ss_kept, "ruler_calibrated"),
    ifelse(
      ss_kept == "tick_lattice_unboxed",
      "tick_lattice_unboxed",
      ifelse(
        ss_kept %in% c("mp_estimated", "predicted_from_megapixels") | (is.finite(ppm_kept) & ppm_kept > 0 & !nzchar(ss_kept)),
        "mp_estimated",
        "uncalibrated"
      )
    )
  )

  pull_db_then_manifest_num <- function(cn) {
    out_v <- if (cn %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows[[cn]])) else rep(NA_real_, n_kept)
    out_v[!is.finite(out_v)] <- NA_real_
    if (!is.null(m_hdr) && cn %in% names(m_hdr) && any(!is.na(m_hit_vec) & !is.finite(out_v))) {
      fill_mask <- !is.na(m_hit_vec) & !is.finite(out_v)
      cand_v <- suppressWarnings(as.numeric(m_hdr[[cn]][m_hit_vec[fill_mask]]))
      fin_m <- is.finite(cand_v)
      if (any(fin_m)) {
        fill_pos <- which(fill_mask)
        out_v[fill_pos[fin_m]] <- cand_v[fin_m]
      }
    }
    out_v
  }
  pull_sub_chr <- function(cn) {
    if (cn %in% names(sub_rows)) as.character(sub_rows[[cn]]) else rep(NA_character_, n_kept)
  }
  pull_sub_int <- function(cn) {
    if (cn %in% names(sub_rows)) suppressWarnings(as.integer(sub_rows[[cn]])) else rep(NA_integer_, n_kept)
  }
  pull_sub_num <- function(cn) {
    if (cn %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows[[cn]])) else rep(NA_real_, n_kept)
  }

  bbox_x1_v <- suppressWarnings(as.integer(round(pull_db_then_manifest_num("bbox_x1"))))
  bbox_y1_v <- suppressWarnings(as.integer(round(pull_db_then_manifest_num("bbox_y1"))))
  bbox_x2_v <- suppressWarnings(as.integer(round(pull_db_then_manifest_num("bbox_x2"))))
  bbox_y2_v <- suppressWarnings(as.integer(round(pull_db_then_manifest_num("bbox_y2"))))
  cx_px_v <- pull_db_then_manifest_num("cx_px")
  cy_px_v <- pull_db_then_manifest_num("cy_px")
  miss_cx <- !is.finite(cx_px_v) & is.finite(bbox_x1_v) & is.finite(bbox_x2_v)
  if (any(miss_cx)) cx_px_v[miss_cx] <- 0.5 * (as.numeric(bbox_x1_v[miss_cx]) + as.numeric(bbox_x2_v[miss_cx]))
  miss_cy <- !is.finite(cy_px_v) & is.finite(bbox_y1_v) & is.finite(bbox_y2_v)
  if (any(miss_cy)) cy_px_v[miss_cy] <- 0.5 * (as.numeric(bbox_y1_v[miss_cy]) + as.numeric(bbox_y2_v[miss_cy]))

  ltl_px_v <- pull_db_then_manifest_num("lamina_trace_length_px")
  ltl_mm_v <- ifelse(
    cal_kept & is.finite(ppm_kept) & ppm_kept > 0 & is.finite(ltl_px_v),
    round(ltl_px_v / ppm_kept, 4),
    NA_real_
  )

  abax_flag_v <- rep(NA, n_kept)
  if ("is_abaxial_pressed_leaf" %in% names(sub_rows)) {
    raw_ab <- toupper(trimws(as.character(sub_rows$is_abaxial_pressed_leaf)))
    abax_flag_v[raw_ab %in% c("1", "1.0", "TRUE", "T", "YES")] <- TRUE
    abax_flag_v[raw_ab %in% c("0", "0.0", "FALSE", "F", "NO")] <- FALSE
  }
  if (!is.null(m_hdr) && "is_abaxial_pressed_leaf" %in% names(m_hdr) && any(is.na(abax_flag_v) & !is.na(m_hit_vec))) {
    fill_ab <- which(is.na(abax_flag_v) & !is.na(m_hit_vec))
    m_ab_raw <- toupper(trimws(as.character(m_hdr$is_abaxial_pressed_leaf[m_hit_vec[fill_ab]])))
    abax_flag_v[fill_ab[m_ab_raw %in% c("1", "1.0", "TRUE", "T", "YES")]] <- TRUE
    abax_flag_v[fill_ab[m_ab_raw %in% c("0", "0.0", "FALSE", "F", "NO")]] <- FALSE
  }
  abax_flag_v[is.na(abax_flag_v)] <- FALSE

  coords_list <- merged$coo[ok_idx]
  names(coords_list) <- sid_vec
  meta_df <- data.frame(
    shape_id = sid_vec,
    catalogNumber = cat_vec,
    plant_individual_id = pid_vec,
    leaf_id = lid_vec,
    n_lamina_components = if ("n_lamina_components" %in% names(sub_rows)) pmax(1L, ifelse(is.na(suppressWarnings(as.integer(sub_rows$n_lamina_components))), 1L, suppressWarnings(as.integer(sub_rows$n_lamina_components)))) else rep(1L, n_kept),
    enclosed_void_frac = if ("enclosed_void_frac" %in% names(sub_rows)) ifelse(is.na(suppressWarnings(as.numeric(sub_rows$enclosed_void_frac))), 0.0, suppressWarnings(as.numeric(sub_rows$enclosed_void_frac))) else rep(0.0, n_kept),
    qc_flag_reason = if ("qc_flag_reason" %in% names(sub_rows)) ifelse(is.na(sub_rows$qc_flag_reason) | !nzchar(trimws(as.character(sub_rows$qc_flag_reason))), "PASS", trimws(as.character(sub_rows$qc_flag_reason))) else rep("PASS", n_kept),
    solidity = sol_kept,
    dissection_index = merged$di[ok_idx],
    archetype_score = arch_kept,
    pixels_per_mm = ppm_kept,
    is_calibrated_mm = cal_kept,
    scale_source = resolved_ss_kept,
    mask_is_upright = ifelse(
      is.na(db_upright_vec[ok_idx]),
      vapply(coords_list, function(c_el) isTRUE(attr(c_el, "lm3_from_tipup")), logical(1L)),
      as.logical(db_upright_vec[ok_idx])
    ),
    bbox_x1 = bbox_x1_v,
    bbox_y1 = bbox_y1_v,
    bbox_x2 = bbox_x2_v,
    bbox_y2 = bbox_y2_v,
    cx_px = cx_px_v,
    cy_px = cy_px_v,
    base_angle_deg = pull_db_then_manifest_num("base_angle_deg"),
    apex_angle_deg = pull_db_then_manifest_num("apex_angle_deg"),
    lamina_curvature = pull_db_then_manifest_num("lamina_curvature"),
    lamina_trace_length_px = ltl_px_v,
    lamina_trace_length_mm = ltl_mm_v,
    basal_cordate_index = pull_db_then_manifest_num("basal_cordate_index"),
    margin_serration_index = pull_db_then_manifest_num("margin_serration_index"),
    margin_tooth_acuteness = pull_db_then_manifest_num("margin_tooth_acuteness"),
    terminal_lobe_area_frac = pull_db_then_manifest_num("terminal_lobe_area_frac"),
    area_mm2 = {
      db_am <- if ("area_mm2" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$area_mm2)) else rep(NA_real_, n_kept)
      db_am[!is.finite(db_am)] <- NA_real_
      am_m <- pull_manifest_num("area_mm2", db_am)
      out_am <- ifelse(is.finite(am_m), am_m, ifelse(cal_kept & is.finite(ppm_kept) & ppm_kept > 0, round(merged$area[ok_idx] / (ppm_kept^2), 2), NA_real_))
      out_am[!cal_kept] <- NA_real_
      out_am
    },
    petiole_length_mm = {
      out_pl <- suppressWarnings(as.numeric(sub_rows$petiole_length_mm))
      out_pl[!cal_kept] <- NA_real_
      out_pl
    },
    petiole_width_mm = {
      out_pw <- pull_manifest_num("petiole_width_mm", if ("petiole_width_mm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$petiole_width_mm)) else rep(NA_real_, n_kept))
      out_pw[!cal_kept] <- NA_real_
      out_pw
    },
    blade_length_mm = {
      out_bl <- suppressWarnings(as.numeric(sub_rows$blade_length_mm))
      out_bl[!cal_kept] <- NA_real_
      out_bl
    },
    blade_width_mm = {
      out_bw <- pull_manifest_num("blade_width_mm", if ("blade_width_mm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_width_mm)) else rep(NA_real_, n_kept))
      out_bw[!cal_kept] <- NA_real_
      out_bw
    },
    blade_cie_l_norm = pull_manifest_num("blade_cie_l_norm", if ("blade_cie_l_norm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_cie_l_norm)) else rep(NA_real_, n_kept)),
    blade_cie_l_sd = pull_manifest_num("blade_cie_l_sd", if ("blade_cie_l_sd" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_cie_l_sd)) else rep(NA_real_, n_kept)),
    blade_cie_l_sheet_delta = pull_manifest_num("blade_cie_l_sheet_delta", if ("blade_cie_l_sheet_delta" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_cie_l_sheet_delta)) else rep(NA_real_, n_kept)),
    blade_cie_a_delta = pull_manifest_num("blade_cie_a_delta", if ("blade_cie_a_delta" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_cie_a_delta)) else rep(NA_real_, n_kept)),
    blade_cie_b_delta = pull_manifest_num("blade_cie_b_delta", if ("blade_cie_b_delta" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_cie_b_delta)) else rep(NA_real_, n_kept)),
    blade_cie_l_q85_sheet_delta = pull_db_then_manifest_num("blade_cie_l_q85_sheet_delta"),
    basal_zone_cie_l_delta = pull_db_then_manifest_num("basal_zone_cie_l_delta"),
    is_abaxial_pressed_leaf = as.logical(abax_flag_v),
    abaxial_cie_l_sheet_delta = pull_db_then_manifest_num("abaxial_cie_l_sheet_delta"),
    petiole_to_blade_ratio = suppressWarnings(as.numeric(sub_rows$petiole_to_blade_ratio)),
    is_dissected = eff_diss_vec,
    ect_path = db_ect_vec,
    leaf_zone = db_zone_vec,
    leaf_zone_confidence = pull_manifest_num("leaf_zone_confidence", db_zone_conf_vec),
    dist_from_crown_norm = pull_manifest_num("dist_from_crown_norm", db_dist_cr_vec),
    cy_rel_plant_bbox = pull_manifest_num("cy_rel_plant_bbox", db_cy_rel_vec),
    cy_rel_sheet = pull_manifest_num("cy_rel_sheet"),
    rosette_centroid_x = pull_sub_num("rosette_centroid_x"),
    rosette_centroid_y = pull_sub_num("rosette_centroid_y"),
    AR_inv = pull_sub_num("AR_inv"),
    capitula_count = pull_sub_num("capitula_count"),
    involucre_height_px = pull_sub_num("involucre_height_px"),
    involucre_width_px = pull_sub_num("involucre_width_px"),
    involucre_height_mm = pull_sub_num("involucre_height_mm"),
    involucre_width_mm = pull_sub_num("involucre_width_mm"),
    corymb_head_count_est = pull_sub_int("corymb_head_count_est"),
    radiate_head_fraction = pull_sub_num("radiate_head_fraction"),
    phyllary_tip_darkening_delta = pull_sub_num("phyllary_tip_darkening_delta"),
    phyllary_tip_anthocyanin_delta = pull_sub_num("phyllary_tip_anthocyanin_delta"),
    involucre_basal_cie_l_delta = pull_sub_num("involucre_basal_cie_l_delta"),
    peduncle_cie_l_sheet_delta = pull_sub_num("peduncle_cie_l_sheet_delta"),
    calyculus_flare_index = pull_sub_num("calyculus_flare_index"),
    institutionCode = pull_sub_chr("institutionCode"),
    scientificName = pull_sub_chr("scientificName"),
    species_raw = pull_sub_chr("species_raw"),
    species_standardized = pull_sub_chr("species_standardized"),
    recordedBy = pull_sub_chr("recordedBy"),
    recordNumber = pull_sub_chr("recordNumber"),
    identifiedBy = pull_sub_chr("identifiedBy"),
    determiner_tier = pull_sub_chr("determiner_tier"),
    type_status = pull_sub_chr("type_status"),
    county = pull_sub_chr("county"),
    stateProvince = pull_sub_chr("stateProvince"),
    locality = pull_sub_chr("locality"),
    decimalLatitude = pull_sub_num("decimalLatitude"),
    decimalLongitude = pull_sub_num("decimalLongitude"),
    coordinateUncertainty = pull_sub_num("coordinateUncertainty"),
    year = pull_sub_int("year"),
    month = pull_sub_int("month"),
    day = pull_sub_int("day"),
    eventDate = pull_sub_chr("eventDate"),
    doy = pull_sub_num("doy"),
    regional_group = pull_sub_chr("regional_group"),
    infraspecific_taxon = pull_sub_chr("infraspecific_taxon"),
    is_primary_duplicate = pull_sub_int("is_primary_duplicate"),
    duplicate_count = pull_sub_int("duplicate_count"),
    spatial_balance_weight = pull_sub_num("spatial_balance_weight"),
    typified_basionym = pull_sub_chr("typified_basionym"),
    verbatim_taxon_concept = pull_sub_chr("verbatim_taxon_concept"),
    concept_synonymy_status = pull_sub_chr("concept_synonymy_status"),
    identification_year = pull_sub_int("identification_year"),
    is_open_or_hybrid_concept = pull_sub_int("is_open_or_hybrid_concept"),
    georef_qc_status = pull_sub_chr("georef_qc_status"),
    is_centroid_or_inst_coord = pull_sub_int("is_centroid_or_inst_coord"),
    coord_quarantine_reason = pull_sub_chr("coord_quarantine_reason"),
    usda_symbol = pull_sub_chr("usda_symbol"),
    usda_county_status = pull_sub_chr("usda_county_status"),
    exsiccatae_key = pull_sub_chr("exsiccatae_key"),
    exsiccatae_taxonomic_conflict = pull_sub_int("exsiccatae_taxonomic_conflict"),
    exsiccatae_gold_taxon = pull_sub_chr("exsiccatae_gold_taxon"),
    centroid_size = cs_kept,
    log_centsize = log_cs_kept,
    aspect_ratio = merged$ar[ok_idx],
    area_px = merged$area[ok_idx],
    circularity = merged$circ[ok_idx],
    mask_source = paste0("sqlite://extracted_leaves/", sid_vec),
    stringsAsFactors = FALSE
  )
  meta_df <- coerce_sqlite_logical_cols(meta_df)
  rownames(meta_df) <- sid_vec
  list(
    coords_list = coords_list,
    meta_df = meta_df
  )
}

.normalize_morph_dynamic_cols_r <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(df)
  }
  for (ic in intersect(c("plant_individual_id", "leaf_id", "efd_harmonic_opt", "efd_k95_pcs", "basal_leaf_count", "cauline_leaf_count", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"), names(df))) {
    df[[ic]] <- suppressWarnings(as.integer(df[[ic]]))
  }
  if ("is_abaxial_pressed_leaf" %in% names(df)) {
    df <- coerce_sqlite_logical_cols(df, cols = "is_abaxial_pressed_leaf")
  }
  if ("selected_morph_features" %in% names(df)) {
    df[["selected_morph_features"]] <- .normalize_semicolon_col_r(df[["selected_morph_features"]])
  }
  dyn_num_cols <- intersect(
    c(
      PACKERA_DB_HARMONIC_COLS, PACKERA_DB_PC_COLS, "within_rosette_efa_dispersion",
      "dissection_index", "q25_solidity", "dissected_leaf_fraction",
      "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox",
      "cauline_solidity_med", "cauline_dissection_index_med", "cauline_aspect_ratio_med",
      "cauline_blade_length_mm", "basal_vs_cauline_length_ratio",
      "basal_vs_cauline_solidity_delta", "cauline_cie_l_q85_sheet_delta", "basal_to_cauline_wool_delta", "heteroblasty_efa_distance",
      "Cauline_PC1", "Cauline_PC2", "Cauline_PC3",
      "ECT_PC1", "ECT_PC2", "ECT_PC3",
      "Cauline_ECT_PC1", "Cauline_ECT_PC2",
      "ect_lobe_complexity", "cauline_ect_lobe_complexity",
      "ect_bilateral_asym", "heteroblasty_ect_distance",
      "area_mm2", "petiole_width_mm", "blade_width_mm",
      "base_angle_deg", "apex_angle_deg", "lamina_curvature",
      "lamina_trace_length_px", "lamina_trace_length_mm",
      "basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness",
      "terminal_lobe_area_frac",
      "blade_cie_l_norm", "blade_cie_l_sd",
      "blade_cie_l_sheet_delta", "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
      "blade_cie_a_delta", "blade_cie_b_delta",
      "abaxial_leaf_fraction", "abaxial_cie_l_sheet_delta",
      "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta",
      "corymb_head_count_est", "radiate_head_fraction", "anthesis_head_count", "fruiting_head_count",
      "phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta",
      "involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta", "calyculus_flare_index",
      "Boot_LRT_Stat_prev", "Boot_LRT_P_prev",
      "log1p_capitulum_count", "log1p_capitulum_count_std", "nearest_type_basionym_d2", "spatial_balance_weight"
    ),
    names(df)
  )
  for (nc in dyn_num_cols) {
    df[[nc]] <- suppressWarnings(as.numeric(df[[nc]]))
  }
  df
}

.prune_stale_run_rows_r <- function(conn, table_name, rid, df) {
  if (is.null(conn) || is.null(df) || nrow(df) == 0L || !("catalogNumber" %in% names(df))) {
    return(invisible(0L))
  }
  cats <- trimws(as.character(df$catalogNumber))
  pids <- if ("plant_individual_id" %in% names(df)) ifelse(is.na(df$plant_individual_id), 0L, suppressWarnings(as.integer(df$plant_individual_id))) else rep(0L, nrow(df))
  pids[is.na(pids)] <- 0L
  keep_keys <- unique(paste(cats, pids, sep = "\r"))
  ex <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT rowid AS _rid, catalogNumber, plant_individual_id FROM %s WHERE run_id = $1;", table_name), params = list(rid)), error = function(e) NULL)
  if (is.null(ex) || nrow(ex) == 0L) {
    return(invisible(0L))
  }
  ex_pids <- ifelse(is.na(ex$plant_individual_id), 0L, suppressWarnings(as.integer(ex$plant_individual_id)))
  ex_pids[is.na(ex_pids)] <- 0L
  ex_keys <- paste(trimws(as.character(ex$catalogNumber)), ex_pids, sep = "\r")
  drop_rids <- as.integer(ex[["_rid"]][!(ex_keys %in% keep_keys)])
  if (length(drop_rids) > 0L) {
    for (s in seq.int(1L, length(drop_rids), by = 500L)) {
      chunk <- drop_rids[s:min(length(drop_rids), s + 499L)]
      tryCatch(DBI::dbExecute(conn, sprintf("DELETE FROM %s WHERE rowid IN (%s);", table_name, paste(chunk, collapse = ","))), error = function(e) NULL)
    }
  }
  invisible(length(drop_rids))
}

persist_step03_to_sqlite <- function(
  efa_df,
  specimen_df,
  vouchers_path = NULL,
  run_id = NULL
) {
  if (!is_sqlite_active_for_inputs(vouchers_path)) {
    return(invisible(FALSE))
  }
  conn <- connect_packera_db(read_only = FALSE)
  if (is.null(conn)) {
    return(invisible(FALSE))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  if (!is.null(efa_df) && nrow(efa_df) > 0L) {
    ind_out <- .normalize_morph_dynamic_cols_r(as.data.frame(efa_df, stringsAsFactors = FALSE))
    ind_out$run_id <- rid
    if (!("plant_individual_id" %in% names(ind_out))) ind_out$plant_individual_id <- 0L
    ind_out$plant_individual_id[is.na(ind_out$plant_individual_id)] <- 0L
    tryCatch(DBI::dbExecute(conn, "DELETE FROM leaf_efa_harmonics WHERE run_id = $1;", params = list(rid)), error = function(e) NULL)
    db_upsert_dataframe(
      conn,
      table_name = "leaf_efa_harmonics",
      df = ind_out,
      pk_cols = c("run_id", "catalogNumber", "plant_individual_id", "leaf_id"),
      preserve_existing_on_null = FALSE
    )
  }

  if (!is.null(specimen_df) && nrow(specimen_df) > 0L) {
    otu_out <- .normalize_morph_dynamic_cols_r(as.data.frame(specimen_df, stringsAsFactors = FALSE))
    otu_out$run_id <- rid
    if (!("plant_individual_id" %in% names(otu_out))) otu_out$plant_individual_id <- 0L
    otu_out$plant_individual_id[is.na(otu_out$plant_individual_id)] <- 0L
    if (!("is_passive" %in% names(otu_out)) && ("authority_weight" %in% names(otu_out))) {
      aw_v <- suppressWarnings(as.numeric(otu_out$authority_weight))
      otu_out$is_passive <- !is.finite(aw_v) | (aw_v <= 0.0)
    }
    tryCatch(DBI::dbExecute(conn, "DELETE FROM otu_morphometrics WHERE run_id = $1;", params = list(rid)), error = function(e) NULL)
    db_upsert_dataframe(
      conn,
      table_name = "otu_morphometrics",
      df = otu_out,
      pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
      preserve_existing_on_null = FALSE
    )
  }
  invisible(TRUE)
}
