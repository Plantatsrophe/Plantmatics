if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_path_conn.R
# SQLite path/run_id resolution, legacy DB migration, Python schema bootstrap, and connection lifecycle.
# ==============================================================================

.get_env_run_dir_r <- function() {
  trimws(env_plantmatics_or_packera("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR"))
}

.get_env_db_path_r <- function() {
  trimws(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"))
}

.migrated_legacy_sqlite_env_r <- new.env(parent = emptyenv())

.migrate_legacy_sqlite_file_r <- function(target_path) {
  if (is.null(target_path) || identical(target_path, ":memory:")) {
    return(target_path)
  }
  if (identical(basename(target_path), "plantmatics.sqlite") && !file.exists(target_path)) {
    legacy_p <- file.path(dirname(target_path), "packera_pipeline.sqlite")
    if (file.exists(legacy_p)) {
      dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
      if (isTRUE(tryCatch(file.rename(legacy_p, target_path), error = function(e) FALSE))) {
        for (sfx in c("-wal", "-shm")) {
          leg_sc <- paste0(legacy_p, sfx)
          if (file.exists(leg_sc)) {
            tryCatch(file.rename(leg_sc, paste0(target_path, sfx)), error = function(e) FALSE)
          }
        }
        if (!exists(target_path, envir = .migrated_legacy_sqlite_env_r, inherits = FALSE)) {
          assign(target_path, TRUE, envir = .migrated_legacy_sqlite_env_r)
          message(sprintf("Migrated legacy SQLite database %s -> %s", legacy_p, target_path))
        }
      }
    }
  }
  target_path
}

resolve_packera_run_id <- function(run_id = NULL) {
  if (!is.null(run_id) && nzchar(trimws(as.character(run_id)[1]))) {
    return(trimws(as.character(run_id)[1]))
  }
  env_rid <- trimws(Sys.getenv("RUN_ID", unset = ""))
  if (nzchar(env_rid)) {
    return(env_rid)
  }
  env_rdir <- .get_env_run_dir_r()
  if (nzchar(env_rdir)) {
    raw_rdir <- gsub("\\\\", "/", trimws(env_rdir))
    raw_rdir <- sub("^\\./+", "", sub("/+$", "", raw_rdir))
    wd_abs <- sub("/+$", "", normalizePath(getwd(), winslash = "/", mustWork = FALSE))
    runs_cand <- if (file.exists("runs")) {
      normalizePath("runs", winslash = "/", mustWork = FALSE)
    } else {
      file.path(wd_abs, "runs")
    }
    runs_abs <- sub("/+$", "", runs_cand)
    cand_path <- if (grepl("^(/|~|[A-Za-z]:/)", raw_rdir)) raw_rdir else file.path(wd_abs, raw_rdir)
    norm_rdir <- sub("/+$", "", normalizePath(cand_path, winslash = "/", mustWork = FALSE))
    if (nzchar(runs_abs) && startsWith(norm_rdir, paste0(runs_abs, "/"))) {
      return(substring(norm_rdir, nchar(runs_abs) + 2L))
    }
    if (!file.exists(cand_path) && startsWith(raw_rdir, "runs/")) {
      return(sub("^runs/+", "", raw_rdir))
    }
    if (grepl("/runs/.+", norm_rdir)) {
      return(sub("^.*/runs/", "", norm_rdir))
    }
    return(basename(norm_rdir))
  }
  "canonical"
}

resolve_packera_db_path <- function(
  db_path = NULL,
  demo = FALSE
) {
  canonical_abs <- normalizePath("data/plantmatics.sqlite", mustWork = FALSE)
  legacy_abs <- normalizePath("data/packera_pipeline.sqlite", mustWork = FALSE)
  if (isTRUE(demo)) {
    if (!is.null(db_path) && nzchar(trimws(db_path))) {
      if (identical(trimws(db_path), ":memory:")) {
        return(":memory:")
      }
      cand_abs <- normalizePath(trimws(db_path), mustWork = FALSE)
      if (!identical(cand_abs, canonical_abs) && !identical(cand_abs, legacy_abs)) {
        return(cand_abs)
      }
    }
    env_rdir <- .get_env_run_dir_r()
    if (nzchar(env_rdir) && dir.exists(env_rdir)) {
      return(file.path(normalizePath(env_rdir, mustWork = FALSE), "demo_ephemeral.sqlite"))
    }
    return(":memory:")
  }

  if (!is.null(db_path) && nzchar(trimws(db_path))) {
    raw <- trimws(db_path)
    if (identical(raw, ":memory:")) {
      return(":memory:")
    }
    return(.migrate_legacy_sqlite_file_r(normalizePath(raw, mustWork = FALSE)))
  }
  env_db <- .get_env_db_path_r()
  if (nzchar(env_db)) {
    if (identical(env_db, ":memory:")) {
      return(":memory:")
    }
    return(.migrate_legacy_sqlite_file_r(normalizePath(env_db, mustWork = FALSE)))
  }
  .migrate_legacy_sqlite_file_r(canonical_abs)
}

resolve_plantmatics_db_path <- resolve_packera_db_path

is_sqlite_active_for_inputs <- function(input_paths = character(0), demo = FALSE) {
  if (isTRUE(demo)) {
    return(FALSE)
  }
  env_db <- .get_env_db_path_r()
  if (nzchar(env_db)) {
    return(TRUE)
  }

  valid_paths <- as.character(input_paths)
  valid_paths <- valid_paths[!is.na(valid_paths) & nzchar(trimws(valid_paths))]
  if (length(valid_paths) == 0L) {
    return(TRUE)
  }

  ws_root <- normalizePath(getwd(), mustWork = FALSE)
  ws_prefix <- paste0(ws_root, .Platform$file.sep)
  for (p in valid_paths) {
    np <- normalizePath(p, mustWork = FALSE)
    if (!identical(np, ws_root) && !startsWith(np, ws_prefix)) {
      return(FALSE)
    }
  }
  TRUE
}

.run_python_db_init_schema <- function(eff_path) {
  if (identical(eff_path, ":memory:")) {
    return(invisible(FALSE))
  }
  py_p <- strsplit(Sys.getenv("PYTHONPATH", unset = ""), .Platform$path.sep, fixed = TRUE)[[1L]]
  py_venv_pp <- if (length(py_p) > 0L && nzchar(py_p[1L]) && file.exists(file.path(py_p[1L], ".venv/bin/python"))) file.path(py_p[1L], ".venv/bin/python") else character(0)
  py_cands <- unique(c(
    if (file.exists(".venv/bin/python")) ".venv/bin/python" else character(0),
    py_venv_pp,
    Sys.which("python3")
  ))
  py_cands <- py_cands[nzchar(py_cands)]
  if (length(py_cands) == 0L) {
    return(invisible(FALSE))
  }
  py_code <- "import sys; from scripts.core.db import get_db_connection; c = get_db_connection(sys.argv[1], init_schema=True); c.close()"
  for (py_bin in py_cands) {
    rc <- tryCatch(
      system2(
        py_bin,
        args = c("-c", shQuote(py_code), shQuote(eff_path)),
        stdout = FALSE,
        stderr = FALSE
      ),
      error = function(e) 1L
    )
    if (identical(as.integer(rc), 0L)) {
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

connect_packera_db <- function(
  db_path = .get_env_db_path_r(),
  read_only = FALSE,
  demo = FALSE
) {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) {
    return(NULL)
  }
  eff_path <- resolve_packera_db_path(
    db_path = if (!is.null(db_path) && length(db_path) == 1L && nzchar(trimws(db_path))) db_path else NULL,
    demo = demo
  )

  if (isTRUE(read_only) && !identical(eff_path, ":memory:")) {
    if (!file.exists(eff_path)) {
      return(NULL)
    }
    conn <- tryCatch(
      DBI::dbConnect(RSQLite::SQLite(), dbname = eff_path, flags = RSQLite::SQLITE_RO),
      error = function(e) NULL
    )
    if (is.null(conn)) {
      return(NULL)
    }
    .apply_packera_sqlite_pragmas_r(conn, read_only = TRUE)
    return(conn)
  }

  if (!identical(eff_path, ":memory:")) {
    dir.create(dirname(eff_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(eff_path)) {
      # Initialize full schema via Python scripts.core.db if creating a fresh DB file from R
      .run_python_db_init_schema(eff_path)
    }
  }

  conn <- tryCatch(
    DBI::dbConnect(RSQLite::SQLite(), dbname = eff_path),
    error = function(e) NULL
  )
  if (is.null(conn)) {
    return(NULL)
  }
  .apply_packera_sqlite_pragmas_r(conn, read_only = FALSE)

  cur_ver <- .get_sqlite_user_version_r(conn)
  v_fields_chk <- if (db_table_exists_safe(conn, "vouchers")) tryCatch(DBI::dbListFields(conn, "vouchers"), error = function(e) character(0)) else c("typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status", "identification_year", "is_open_or_hybrid_concept")
  el_fields_chk <- if (db_table_exists_safe(conn, "extracted_leaves")) tryCatch(DBI::dbListFields(conn, "extracted_leaves"), error = function(e) character(0)) else c("leaf_zone", "ect_path", "bbox_x1", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "is_abaxial_pressed_leaf")
  vrp_fields_chk <- if (db_table_exists_safe(conn, "voucher_reproductive_phenology")) tryCatch(DBI::dbListFields(conn, "voucher_reproductive_phenology"), error = function(e) character(0)) else c("corymb_head_count_est", "radiate_head_fraction")
  leh_fields_chk <- if (db_table_exists_safe(conn, "leaf_efa_harmonics")) tryCatch(DBI::dbListFields(conn, "leaf_efa_harmonics"), error = function(e) character(0)) else c("scale_source", "bbox_x1", "lamina_trace_length_mm", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "is_abaxial_pressed_leaf")
  om_fields_chk <- if (db_table_exists_safe(conn, "otu_morphometrics")) tryCatch(DBI::dbListFields(conn, "otu_morphometrics"), error = function(e) character(0)) else c("cauline_leaf_count", "ECT_PC1", "scale_source", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "abaxial_leaf_fraction")
  ec_fields_chk <- if (db_table_exists_safe(conn, "environmental_cache")) tryCatch(DBI::dbListFields(conn, "environmental_cache"), error = function(e) character(0)) else c("ssurgo_lithic_series_pct")
  srf_fields_chk <- if (db_table_exists_safe(conn, "spatial_rf_predictions")) tryCatch(DBI::dbListFields(conn, "spatial_rf_predictions"), error = function(e) character(0)) else c(
    "ecotype_regime", "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx",
    "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight",
    "allochronic_bimodality_delta_bic", "nearest_type_niche_basionym", "nearest_type_basionym_edaphic_z",
    "env_mahal_d2", "env_mcd_outlier_flag"
  )
  gs_fields_chk <- if (db_table_exists_safe(conn, "genomic_samples")) tryCatch(DBI::dbListFields(conn, "genomic_samples"), error = function(e) character(0)) else c(
    "exsiccatae_key", "genbank_accession", "coordinateUncertainty",
    "is_centroid_or_inst_coord", "coord_quarantine_reason",
    "typified_basionym", "concept_synonymy_status", "is_open_or_hybrid_concept"
  )
  ct_fields_chk <- if (db_table_exists_safe(conn, "curation_triage")) tryCatch(DBI::dbListFields(conn, "curation_triage"), error = function(e) character(0)) else c(
    "w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex",
    "geno_deme_weight_scale", "portal_metadata_lag_flag", "portal_lag_trigger_gate",
    "is_marginal_range_extension"
  )
  mv_fields_chk <- if (db_table_exists_safe(conn, "v_integrative_triage_master")) tryCatch(DBI::dbListFields(conn, "v_integrative_triage_master"), error = function(e) character(0)) else character(0)
  v_has_concept <- all(c("typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status", "identification_year", "is_open_or_hybrid_concept") %in% v_fields_chk)
  el_has_zone <- all(c("leaf_zone", "ect_path", "bbox_x1", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "is_abaxial_pressed_leaf") %in% el_fields_chk)
  vrp_has_latest <- all(c("corymb_head_count_est", "radiate_head_fraction") %in% vrp_fields_chk)
  leh_has_latest <- all(c("scale_source", "bbox_x1", "lamina_trace_length_mm", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "is_abaxial_pressed_leaf") %in% leh_fields_chk)
  om_has_cauline <- all(c("cauline_leaf_count", "ECT_PC1", "scale_source", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "abaxial_leaf_fraction") %in% om_fields_chk)
  ec_has_latest <- all(c("ssurgo_lithic_series_pct") %in% ec_fields_chk)
  srf_has_latest <- all(c(
    "ecotype_regime", "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx",
    "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight",
    "allochronic_bimodality_delta_bic", "nearest_type_niche_basionym", "nearest_type_basionym_edaphic_z",
    "env_mahal_d2", "env_mcd_outlier_flag"
  ) %in% srf_fields_chk)
  gs_has_latest <- all(c(
    "exsiccatae_key", "genbank_accession", "coordinateUncertainty",
    "is_centroid_or_inst_coord", "coord_quarantine_reason",
    "typified_basionym", "concept_synonymy_status", "is_open_or_hybrid_concept"
  ) %in% gs_fields_chk)
  ct_has_latest <- all(c(
    "w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex",
    "geno_deme_weight_scale", "portal_metadata_lag_flag", "portal_lag_trigger_gate",
    "is_marginal_range_extension"
  ) %in% ct_fields_chk)
  mv_has_latest <- identical(eff_path, ":memory:") || all(c(
    "is_passive", "scale_source", "AR_inv_std", "reproductive_trait_status",
    "corymb_head_count_est", "radiate_head_fraction",
    "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "abaxial_leaf_fraction",
    "ecoregion_can1", "ecoregion_can2", "w_geno_kappa",
    "geno_predicted_taxon", "geno_linkage_tier", "is_geno_admixed",
    "coalescent_gdi", "hybrid_generation_class", "cytonuclear_discordant",
    "public_plastome_match", "public_its_match", "waxy_lineage_match",
    "public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index",
    "genomic_propagation_tier", "cluster_propagated_gdi",
    "cytology_linkage_tier", "voucher_sporophytic_2n", "voucher_ploidy_level",
    "taxon_reported_cytotypes", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass",
    "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
    "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m",
    "hydro_edaphic_niche_regime", "wetland_fidelity_concordance", "basionym_synonymy_type",
    "w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex",
    "geno_deme_weight_scale", "portal_metadata_lag_flag", "portal_lag_trigger_gate",
    "is_marginal_range_extension",
    "backbone_basionym_name", "protologue_url", "natureserve_global_rank",
    "nwpl_indicator_code", "nwpl_wetland_score",
    "is_centroid_or_inst_coord", "coord_quarantine_reason",
    "typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status",
    "identification_year", "is_open_or_hybrid_concept",
    "env_mahal_d2", "env_mcd_outlier_flag",
    "ssurgo_lithic_series_pct",
    "ecotype_regime", "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx",
    "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight",
    "allochronic_bimodality_delta_bic", "nearest_type_niche_basionym", "nearest_type_basionym_edaphic_z"
  ) %in% mv_fields_chk)
  if (cur_ver >= PLANTMATICS_SCHEMA_USER_VERSION &&
    v_has_concept && el_has_zone && vrp_has_latest && leh_has_latest && om_has_cauline &&
    ec_has_latest && srf_has_latest && gs_has_latest && ct_has_latest && mv_has_latest &&
    (identical(eff_path, ":memory:") || (
      db_table_exists_safe(conn, "v_integrative_triage_master") &&
        db_table_exists_safe(conn, "genomic_samples") &&
        db_table_exists_safe(conn, "genetic_lineage_assignments") &&
        db_table_exists_safe(conn, "cytological_records") &&
        db_table_exists_safe(conn, "taxon_cytology_summary") &&
        db_table_exists_safe(conn, "taxon_backbone_profiles") &&
        db_table_exists_safe(conn, "natureserve_subnational_ranges")
    ))) {
    return(conn)
  }

  schema_migrated <- .migrate_dynamic_sqlite_schema_r(conn)
  if (!identical(eff_path, ":memory:")) {
    view_needs_refresh <- schema_migrated || (cur_ver < PLANTMATICS_SCHEMA_USER_VERSION) ||
      !db_table_exists_safe(conn, "genomic_samples") ||
      !db_table_exists_safe(conn, "genetic_lineage_assignments") ||
      !db_table_exists_safe(conn, "cytological_records") ||
      !db_table_exists_safe(conn, "taxon_cytology_summary") ||
      !db_table_exists_safe(conn, "taxon_backbone_profiles") ||
      !db_table_exists_safe(conn, "natureserve_subnational_ranges")
    if (!view_needs_refresh && db_table_exists_safe(conn, "v_integrative_triage_master")) {
      v_cols <- tryCatch(DBI::dbListFields(conn, "v_integrative_triage_master"), error = function(e) character(0))
      req_v_cols <- c(
        "PC15", "within_rosette_efa_dispersion", "selected_morph_features", "climatic_z_score",
        "bio19_precip_coldest_qtr", "dissection_index", "q25_solidity",
        "dissected_leaf_fraction", "basal_leaf_count", "cauline_leaf_count",
        "basal_vs_cauline_length_ratio", "heteroblasty_efa_distance", "is_type_specimen",
        "ECT_PC1", "Cauline_ECT_PC1", "ect_lobe_complexity", "heteroblasty_ect_distance",
        "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance",
        "is_typical_morphotype", "min_congener_mahal_d2", "min_congener_taxon",
        "open_set_novelty_p", "is_confident_label_noise",
        "type_edaphic_z_score", "nearest_type_niche_taxon", "matches_type_lithology",
        "cryptic_cluster_cohesion_pass", "cryptic_non_clinal_pass", "resurrection_basionym",
        "hybrid_zone_cohesion_pass", "hybrid_parent_1", "hybrid_parent_2",
        "cda_oof_posterior_prob", "self_morph_mahal_d2", "open_set_df",
        "is_passive", "scale_source",
        "AR_inv_std", "reproductive_trait_status", "corymb_head_count_est", "radiate_head_fraction",
        "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "abaxial_leaf_fraction",
        "ecoregion_can1", "ecoregion_can2",
        "w_geno_kappa", "geno_predicted_taxon", "geno_linkage_tier", "is_geno_admixed",
        "coalescent_gdi", "hybrid_generation_class", "cytonuclear_discordant",
        "public_plastome_match", "public_its_match", "waxy_lineage_match",
        "public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index",
        "genomic_propagation_tier", "cluster_propagated_gdi",
        "cytology_linkage_tier", "voucher_sporophytic_2n", "voucher_ploidy_level",
        "taxon_reported_cytotypes", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass",
        "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
        "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m",
        "hydro_edaphic_niche_regime", "wetland_fidelity_concordance", "basionym_synonymy_type",
        "w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex",
        "geno_deme_weight_scale", "portal_metadata_lag_flag", "portal_lag_trigger_gate",
        "is_marginal_range_extension",
        "backbone_basionym_name", "protologue_url", "natureserve_global_rank",
        "nwpl_indicator_code", "nwpl_wetland_score",
        "is_centroid_or_inst_coord", "coord_quarantine_reason",
        "typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status",
        "identification_year", "is_open_or_hybrid_concept",
        "env_mahal_d2", "env_mcd_outlier_flag",
        "ssurgo_lithic_series_pct",
        "ecotype_regime", "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx",
        "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight",
        "allochronic_bimodality_delta_bic", "nearest_type_niche_basionym", "nearest_type_basionym_edaphic_z"
      )
      if (length(v_cols) > 0L && !all(req_v_cols %in% v_cols)) {
        view_needs_refresh <- TRUE
      }
    }
    if (isTRUE(view_needs_refresh)) {
      disconnect_packera_db(conn)
      .run_python_db_init_schema(eff_path)
      conn <- tryCatch(
        DBI::dbConnect(RSQLite::SQLite(), dbname = eff_path),
        error = function(e) NULL
      )
      if (!is.null(conn)) {
        .apply_packera_sqlite_pragmas_r(conn, read_only = FALSE)
      }
    } else if (!is.null(conn) && db_table_exists_safe(conn, "v_integrative_triage_master")) {
      tryCatch(
        DBI::dbExecute(conn, sprintf("PRAGMA user_version = %d;", PLANTMATICS_SCHEMA_USER_VERSION)),
        error = function(e) NULL
      )
    }
  }
  conn
}

disconnect_packera_db <- function(conn) {
  if (!is.null(conn)) {
    is_val <- tryCatch(DBI::dbIsValid(conn), error = function(e) FALSE)
    if (isTRUE(is_val)) {
      tryCatch(suppressWarnings(DBI::dbDisconnect(conn)), error = function(e) NULL)
    }
  }
  invisible(NULL)
}
