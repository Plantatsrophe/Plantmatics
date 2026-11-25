if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step04_persist.R
# Step 04 OTU morphometrics SQLite joins, idempotent rerun cleanup, and persistence.
# ==============================================================================

fetch_step04_morphometrics_from_db <- function(
  vouchers_path = NULL,
  harmonics_path = NULL,
  active_catalogs = NULL,
  run_id = NULL
) {
  if (!is_sqlite_active_for_inputs(c(vouchers_path, harmonics_path))) {
    return(NULL)
  }
  conn <- connect_packera_db(read_only = TRUE)
  if (is.null(conn)) {
    return(NULL)
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "otu_morphometrics")) {
    return(NULL)
  }

  rid <- resolve_packera_run_id(run_id)
  cnt <- tryCatch(
    DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM otu_morphometrics WHERE run_id = $1;", params = list(rid))$n[1],
    error = function(e) 0L
  )
  if (is.na(cnt) || cnt == 0L) {
    return(NULL)
  }

  om_cols <- tryCatch(DBI::dbListFields(conn, "otu_morphometrics"), error = function(e) character(0))
  pi_cols <- if (db_table_exists_safe(conn, "plant_individuals")) tryCatch(DBI::dbListFields(conn, "plant_individuals"), error = function(e) character(0)) else character(0)
  vrp_cols <- if (db_table_exists_safe(conn, "voucher_reproductive_phenology")) tryCatch(DBI::dbListFields(conn, "voucher_reproductive_phenology"), error = function(e) character(0)) else character(0)
  v_cols <- if (db_table_exists_safe(conn, "vouchers")) tryCatch(DBI::dbListFields(conn, "vouchers"), error = function(e) character(0)) else character(0)

  extra_v_cols <- intersect(
    c(
      "usda_symbol", "usda_county_status", "exsiccatae_key",
      "exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon",
      "is_centroid_or_inst_coord", "coord_quarantine_reason",
      "typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status",
      "identification_year", "is_open_or_hybrid_concept"
    ),
    v_cols
  )

  dyn_select_exprs <- character(0)
  dyn_joined_cols <- character(0)
  used_pi <- FALSE
  used_vrp <- FALSE
  used_v <- length(extra_v_cols) > 0L

  add_coalesce_col <- function(col_nm, order_tbls) {
    exprs <- character(0)
    for (tb in order_tbls) {
      if (tb == "om" && (col_nm %in% om_cols)) {
        exprs <- c(exprs, if (col_nm == "scale_source") "NULLIF(TRIM(om.scale_source), '')" else paste0("om.", col_nm))
      } else if (tb == "pi" && (col_nm %in% pi_cols)) {
        exprs <- c(exprs, if (col_nm == "scale_source") "NULLIF(TRIM(pi.scale_source), '')" else paste0("pi.", col_nm))
        used_pi <<- TRUE
      } else if (tb == "vrp" && (col_nm %in% vrp_cols)) {
        exprs <- c(exprs, sprintf("CASE WHEN vrp.catalogNumber IS NOT NULL THEN vrp.%s ELSE vrp0.%s END", col_nm, col_nm))
        used_vrp <<- TRUE
      } else if (tb == "v" && (col_nm %in% v_cols)) {
        exprs <- c(exprs, if (col_nm == "scale_source") "NULLIF(TRIM(v.scale_source), '')" else paste0("v.", col_nm))
        used_v <<- TRUE
      }
    }
    if (length(exprs) == 0L) return(invisible(NULL))
    sel_e <- if (length(exprs) == 1L) exprs[1L] else sprintf("COALESCE(%s)", paste(exprs, collapse = ", "))
    dyn_select_exprs <<- c(dyn_select_exprs, sprintf("%s AS %s_joined", sel_e, col_nm))
    dyn_joined_cols <<- c(dyn_joined_cols, col_nm)
    invisible(NULL)
  }

  for (cn_pi in c("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "scale_source", "is_calibrated_mm", "blade_cie_l_sheet_delta", "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta")) {
    add_coalesce_col(cn_pi, c("om", "pi", "vrp", "v"))
  }
  for (cn_arch in c(
    "lamina_trace_length_mm", "basal_cordate_index", "base_angle_deg", "apex_angle_deg",
    "lamina_curvature", "margin_serration_index", "margin_tooth_acuteness",
    "terminal_lobe_area_frac", "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "abaxial_leaf_fraction", "abaxial_cie_l_sheet_delta"
  )) {
    add_coalesce_col(cn_arch, c("om", "pi"))
  }
  add_coalesce_col("AR_inv", c("om", "vrp"))
  for (cn_vrp in c("capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm", "corymb_head_count_est", "radiate_head_fraction", "anthesis_head_count", "fruiting_head_count")) {
    add_coalesce_col(cn_vrp, c("vrp", "pi", "om", "v"))
  }

  extra_v_sql <- if (length(extra_v_cols) > 0L) {
    paste0(", ", paste0("v.", extra_v_cols, " AS ", extra_v_cols, "_v_joined", collapse = ", "))
  } else ""
  dyn_sql_part <- if (length(dyn_select_exprs) > 0L) {
    paste0(", ", paste(dyn_select_exprs, collapse = ", "))
  } else ""
  pi_join <- if (isTRUE(used_pi)) {
    "LEFT JOIN plant_individuals pi ON pi.catalogNumber = om.catalogNumber AND pi.plant_individual_id = om.plant_individual_id "
  } else ""
  vrp_join <- if (isTRUE(used_vrp)) {
    paste0(
      "LEFT JOIN voucher_reproductive_phenology vrp ON vrp.catalogNumber = om.catalogNumber AND vrp.plant_individual_id = om.plant_individual_id ",
      "LEFT JOIN voucher_reproductive_phenology vrp0 ON vrp0.catalogNumber = om.catalogNumber AND vrp0.plant_individual_id = 0 "
    )
  } else ""
  v_join <- if (isTRUE(used_v)) {
    "LEFT JOIN vouchers v ON v.catalogNumber = om.catalogNumber "
  } else ""

  sql <- paste0(
    "SELECT om.*", dyn_sql_part, extra_v_sql, " ",
    "FROM otu_morphometrics om ",
    pi_join, vrp_join, v_join,
    "WHERE om.run_id = $1 ",
    "ORDER BY om.catalogNumber, om.plant_individual_id;"
  )
  df <- tryCatch(DBI::dbGetQuery(conn, sql, params = list(rid)), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) {
    return(NULL)
  }

  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    df <- df[trimws(as.character(df$catalogNumber)) %in% cats, , drop = FALSE]
    if (nrow(df) == 0L) {
      return(NULL)
    }
  }

  for (cn in dyn_joined_cols) {
    jc <- paste0(cn, "_joined")
    if (jc %in% names(df)) {
      df[[cn]] <- df[[jc]]
      df[[jc]] <- NULL
    }
  }
  for (vcol in extra_v_cols) {
    vj <- paste0(vcol, "_v_joined")
    if (vj %in% names(df)) {
      df[[vcol]] <- df[[vj]]
      df[[vj]] <- NULL
    }
  }
  df$run_id <- NULL
  df <- coerce_sqlite_logical_cols(df)
  has_prev_s04 <- ("cda_predicted_taxon" %in% names(df) && any(!is.na(df$cda_predicted_taxon))) ||
    ("gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster)))
  if (isTRUE(has_prev_s04)) {
    restored_db <- .restore_raw_step03_from_leaf_harmonics_db(conn, rid, df)
    df <- restored_db$df
    restored_map <- restored_db$restored_mask
    if (!is.null(harmonics_path) && nzchar(harmonics_path) && file.exists(harmonics_path)) {
      hm_raw <- tryCatch(
        if (requireNamespace("data.table", quietly = TRUE)) {
          as.data.frame(data.table::fread(harmonics_path, showProgress = FALSE))
        } else {
          utils::read.csv(harmonics_path, stringsAsFactors = FALSE, check.names = FALSE)
        },
        error = function(e) NULL
      )
      if (!is.null(hm_raw) && "catalogNumber" %in% names(hm_raw) && nrow(hm_raw) > 0L) {
        if (!("plant_individual_id" %in% names(hm_raw))) hm_raw$plant_individual_id <- 0L
        hm_raw$plant_individual_id[is.na(hm_raw$plant_individual_id)] <- 0L
        k_db <- paste0(trimws(as.character(df$catalogNumber)), "||", as.integer(df$plant_individual_id))
        k_hm <- paste0(trimws(as.character(hm_raw$catalogNumber)), "||", as.integer(hm_raw$plant_individual_id))
        if (!any(duplicated(k_hm))) {
          m_idx <- match(k_db, k_hm)
          ok_m <- which(!is.na(m_idx))
          if (length(ok_m) > 0L) {
            raw_s03_cols <- intersect(
              c(
                PACKERA_DB_HARMONIC_COLS, PACKERA_DB_PC_COLS,
                "ECT_PC1", "ECT_PC2", "ECT_PC3", "Cauline_PC1", "Cauline_PC2", "Cauline_PC3",
                "Cauline_ECT_PC1", "Cauline_ECT_PC2", "aspect_ratio", "solidity", "q25_solidity",
                "dissection_index", "dissected_leaf_fraction", "within_rosette_efa_dispersion",
                "heteroblasty_efa_distance", "heteroblasty_ect_distance", "ect_lobe_complexity",
                "cauline_ect_lobe_complexity",
                "base_angle_deg", "apex_angle_deg", "lamina_curvature", "lamina_trace_length_mm",
                "basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness",
                "terminal_lobe_area_frac",
                "blade_cie_l_norm", "blade_cie_l_sd",
                "blade_cie_l_sheet_delta", "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
                "blade_cie_a_delta", "blade_cie_b_delta",
                "is_abaxial_pressed_leaf", "abaxial_leaf_fraction", "abaxial_cie_l_sheet_delta",
                "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta",
                "centroid_size", "log_centsize", "scale_source"
              ),
              names(hm_raw)
            )
            for (rc in raw_s03_cols) {
              fill_pos <- if (rc %in% names(restored_map)) {
                ok_m[!restored_map[[rc]][ok_m]]
              } else {
                ok_m
              }
              if (length(fill_pos) > 0L) {
                if (!(rc %in% names(df))) {
                  df[[rc]] <- hm_raw[[rc]][rep(NA_integer_, nrow(df))]
                }
                df[[rc]][fill_pos] <- hm_raw[[rc]][m_idx[fill_pos]]
              }
            }
            df <- coerce_sqlite_logical_cols(df)
          }
        }
      }
    }
    stale_s04_cols <- c(
      "gmm_cluster", "gmm_uncertainty", "cda_predicted_taxon", "cda_posterior_prob",
      "cda_oof_posterior_prob", "cl_target_threshold", "self_morph_mahal_d2", "self_morph_df",
      "cda_second_taxon", "cda_second_posterior_prob", "cda_posterior_margin",
      "loocv_accuracy", "can1", "can2", "can3", "ecoregion_cda_predicted",
      "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2",
      "morphometric_outlier_flag", "misidentification_flag", "triage_priority",
      "discordance_reason", "has_empirical_type_anchor", "type_morph_mahal_d2",
      "nearest_type_morph_taxon", "nearest_type_basionym", "nearest_type_basionym_d2",
      "type_morph_concordance", "is_typical_morphotype",
      "min_congener_mahal_d2", "min_congener_taxon", "open_set_novelty_p", "open_set_df",
      "is_confident_label_noise", "selected_morph_features"
    )
    df <- df[, setdiff(names(df), stale_s04_cols), drop = FALSE]
  }
  df
}

persist_step04_to_sqlite <- function(
  out_df,
  vouchers_path = NULL,
  run_id = NULL,
  focal_taxon = NULL,
  intraspecific_mode = FALSE
) {
  if (!is_sqlite_active_for_inputs(vouchers_path) || is.null(out_df) || nrow(out_df) == 0L) {
    return(invisible(FALSE))
  }
  conn <- connect_packera_db(read_only = FALSE)
  if (is.null(conn)) {
    return(invisible(FALSE))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid, focal_taxon = focal_taxon, intraspecific_mode = intraspecific_mode)

  df_db <- .normalize_morph_dynamic_cols_r(as.data.frame(out_df, stringsAsFactors = FALSE))
  df_db$run_id <- rid
  if (!("plant_individual_id" %in% names(df_db))) df_db$plant_individual_id <- 0L
  df_db$plant_individual_id[is.na(df_db$plant_individual_id)] <- 0L

  db_upsert_dataframe(
    conn,
    table_name = "otu_morphometrics",
    df = df_db,
    pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
    preserve_existing_on_null = TRUE
  )
  .prune_stale_run_rows_r(conn, "otu_morphometrics", rid, df_db)

  if (db_table_exists_safe(conn, "run_active_cohort") && db_table_exists_safe(conn, "vouchers")) {
    # Update authority_weight and is_promoted_bronze on existing cohort rows without re-inserting pruned ramets
    existing_rac <- tryCatch(
      DBI::dbGetQuery(conn, "SELECT catalogNumber, COALESCE(plant_individual_id, 0) AS plant_individual_id FROM run_active_cohort WHERE run_id = $1;", params = list(rid)),
      error = function(e) data.frame()
    )
    if (nrow(existing_rac) > 0L) {
      ex_keys <- paste(trimws(as.character(existing_rac$catalogNumber)), as.integer(existing_rac$plant_individual_id), sep = "__")
      df_keys <- paste(trimws(as.character(df_db$catalogNumber)), as.integer(df_db$plant_individual_id), sep = "__")
      coh_sub <- df_db[df_keys %in% ex_keys, , drop = FALSE]
    } else {
      v_cats <- tryCatch(DBI::dbGetQuery(conn, "SELECT catalogNumber FROM vouchers;")$catalogNumber, error = function(e) character(0))
      coh_sub <- df_db[df_db$catalogNumber %in% v_cats, , drop = FALSE]
    }
    if (nrow(coh_sub) > 0L) {
      coh_cols <- intersect(
        c(
          "run_id", "catalogNumber", "plant_individual_id", "species_standardized",
          "determiner_tier", "type_status", "is_type_specimen",
          "is_promoted_bronze", "authority_weight", "spatial_balance_weight"
        ),
        names(coh_sub)
      )
      db_upsert_dataframe(
        conn,
        table_name = "run_active_cohort",
        df = coh_sub[, coh_cols, drop = FALSE],
        pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null = TRUE
      )
    }
  }
  invisible(TRUE)
}
