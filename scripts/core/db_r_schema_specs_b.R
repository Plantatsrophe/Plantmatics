if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_schema_specs_b.R
# Dynamic SQLite column migration specs (part B: cytology through NatureServe) and schema migrator.
# ==============================================================================

.dynamic_sqlite_table_specs_b_r <- function(table_name) {
  if (identical(table_name, "cytological_records")) {
    return(c(
      source_db = "TEXT",
      source_record_id = "TEXT",
      original_taxon_name = "TEXT",
      basionym_name = "TEXT",
      species_standardized = "TEXT",
      catalogNumber = "TEXT",
      exsiccatae_key = "TEXT",
      cytology_linkage_tier = "TEXT",
      raw_count_str = "TEXT",
      gametophytic_n = "INTEGER",
      sporophytic_2n = "INTEGER",
      base_x = "INTEGER",
      ploidy_level = "INTEGER",
      cytotype_label = "TEXT",
      b_chromosomes = "INTEGER",
      is_aneuploid_or_odd_ploid = "INTEGER",
      is_approximate = "INTEGER",
      c_value_1c_pg = "REAL",
      c_value_2c_pg = "REAL",
      breeding_system = "TEXT",
      genbank_accession = "TEXT",
      collector_surname = "TEXT",
      record_number_clean = "TEXT",
      institution_code = "TEXT",
      voucher_citation = "TEXT",
      locality_text = "TEXT",
      reference_citation = "TEXT"
    ))
  }
  if (identical(table_name, "taxon_cytology_summary")) {
    return(c(
      n_cytology_records = "INTEGER",
      n_direct_cytovouchers = "INTEGER",
      modal_sporophytic_2n = "INTEGER",
      modal_base_x = "INTEGER",
      modal_ploidy_level = "INTEGER",
      reported_cytotypes = "TEXT",
      is_multi_cytotype_taxon = "INTEGER",
      has_dysploid_series = "INTEGER",
      has_odd_polyploid_or_apomixis = "INTEGER",
      mean_2c_pg = "REAL",
      breeding_system_summary = "TEXT"
    ))
  }
  if (identical(table_name, "spatial_rf_predictions")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_WC_BIO_COLS)), PACKERA_DB_WC_BIO_COLS),
      c(
        rf_prediction_source = "TEXT",
        n_observed_modalities = "INTEGER",
        climatic_z_score = "REAL",
        primary_env_predictors = "TEXT",
        primary_clim_predictors = "TEXT",
        primary_edaph_predictors = "TEXT",
        type_edaphic_z_score = "REAL",
        nearest_type_niche_taxon = "TEXT",
        matches_type_lithology = "INTEGER",
        ecotype_regime = "TEXT",
        lithic_outcrop_index = "REAL",
        summer_moisture_deficit_idx = "REAL",
        maritime_continentality_idx = "REAL",
        hopkins_adjusted_doy = "REAL",
        hopkins_adjusted_doy_resid = "REAL",
        anthesis_dominance_weight = "REAL",
        allochronic_bimodality_delta_bic = "REAL",
        nearest_type_niche_basionym = "TEXT",
        nearest_type_basionym_edaphic_z = "REAL",
        env_mahal_d2 = "REAL",
        env_mcd_outlier_flag = "INTEGER"
      )
    ))
  }
  if (identical(table_name, "environmental_cache")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_WC_BIO_COLS)), PACKERA_DB_WC_BIO_COLS),
      c(
        ssurgo_lithic_series_pct = "REAL",
        solus_sand_pct = "REAL",
        solus_clay_pct = "REAL",
        solus_silt_pct = "REAL",
        solus_cec = "REAL",
        topo_hli_30m = "REAL",
        elevation_m = "REAL",
        ilr_texture_1 = "REAL",
        ilr_texture_2 = "REAL",
        fine_cfvo_pct = "REAL",
        is_clim_eligible_coord = "INTEGER",
        is_edaphic_eligible_coord = "INTEGER",
        edaphic_coord_weight = "REAL",
        soil_clay = "REAL",
        soil_silt = "REAL",
        soil_cfvo = "REAL",
        soil_clay_soilgrids = "REAL",
        soil_silt_soilgrids = "REAL",
        soil_cfvo_soilgrids = "REAL",
        soil_ca_mg_molar_ratio = "REAL",
        soil_ultramafic_ni_cr_index = "REAL",
        soil_carbonate_c_pct = "REAL",
        is_nwi_wetland = "INTEGER",
        nwi_wetland_system = "TEXT",
        nwi_wetland_dist_m = "REAL"
      )
    ))
  }
  if (identical(table_name, "genomic_samples")) {
    return(c(
      hdplot_retained_loci = "INTEGER",
      inferred_cytotype = "TEXT",
      nquire_delta_logl = "REAL",
      plastid_clade = "TEXT",
      its_ribotype = "TEXT",
      cytonuclear_discordant = "INTEGER",
      public_plastome_match = "TEXT",
      public_its_match = "TEXT",
      waxy_lineage_match = "TEXT",
      exsiccatae_key = "TEXT",
      genbank_accession = "TEXT",
      coordinateUncertainty = "REAL",
      is_centroid_or_inst_coord = "INTEGER DEFAULT 0",
      coord_quarantine_reason = "TEXT",
      typified_basionym = "TEXT",
      concept_synonymy_status = "TEXT",
      is_open_or_hybrid_concept = "INTEGER DEFAULT 0"
    ))
  }
  if (identical(table_name, "genetic_lineage_assignments")) {
    return(c(
      hybrid_index_h = "REAL",
      interspecific_ancestry_q12 = "REAL",
      hybrid_generation_class = "TEXT",
      coancestry_cluster = "TEXT",
      f_branch_max = "REAL",
      f_branch_donor_taxon = "TEXT",
      snaq_gamma_minor = "REAL",
      coalescent_gdi = "REAL",
      bpp_introgression_bf = "REAL",
      mmrr_beta_ibd = "REAL",
      mmrr_beta_ibe = "REAL",
      mmrr_p_ibe = "REAL",
      public_backbone_clade = "TEXT",
      public_phytop_ils_index = "REAL",
      public_phytop_ih_index = "REAL"
    ))
  }
  if (identical(table_name, "taxon_backbone_profiles")) {
    return(c(
      target_genus = "TEXT DEFAULT 'Packera'",
      accepted_scientific_name = "TEXT",
      accepted_authorship = "TEXT",
      basionym_name = "TEXT",
      basionym_authorship = "TEXT",
      homotypic_synonyms_json = "TEXT",
      heterotypic_synonyms_json = "TEXT",
      is_hybrid_taxon = "INTEGER DEFAULT 0",
      hybrid_formula = "TEXT",
      publication_citation = "TEXT",
      publication_year = "INTEGER",
      protologue_url = "TEXT",
      locus_classicus_text = "TEXT",
      tdwg_native_codes_json = "TEXT",
      tdwg_introduced_codes_json = "TEXT",
      natureserve_element_id = "TEXT",
      natureserve_global_rank = "TEXT",
      nwpl_indicator_code = "TEXT",
      nwpl_wetland_score = "REAL",
      backbone_source = "TEXT",
      updated_at = "TEXT"
    ))
  }
  if (identical(table_name, "natureserve_subnational_ranges")) {
    return(c(
      s_rank = "TEXT",
      is_native = "INTEGER DEFAULT 1",
      is_exotic = "INTEGER DEFAULT 0",
      is_historical_or_extirpated = "INTEGER DEFAULT 0"
    ))
  }
  character(0)
}

.dynamic_sqlite_table_specs_r <- function(table_name) {
  c(.dynamic_sqlite_table_specs_a_r(table_name), .dynamic_sqlite_table_specs_b_r(table_name))
}

.migrate_dynamic_sqlite_table_r <- function(conn, table_name) {
  if (is.null(conn) || !db_table_exists_safe(conn, table_name)) {
    return(FALSE)
  }
  specs <- .dynamic_sqlite_table_specs_r(table_name)
  if (length(specs) == 0L) {
    return(FALSE)
  }
  existing <- tryCatch(DBI::dbListFields(conn, table_name), error = function(e) character(0))
  if (length(existing) == 0L) {
    return(FALSE)
  }
  missing_cols <- setdiff(names(specs), existing)
  if (length(missing_cols) == 0L) {
    return(FALSE)
  }

  q_tbl <- as.character(DBI::dbQuoteIdentifier(conn, table_name))
  added_any <- FALSE
  for (cn in missing_cols) {
    col_type <- specs[[cn]]
    q_col <- as.character(DBI::dbQuoteIdentifier(conn, cn))
    ok <- tryCatch(
      {
        DBI::dbExecute(conn, sprintf("ALTER TABLE %s ADD COLUMN %s %s;", q_tbl, q_col, col_type))
        TRUE
      },
      error = function(e) FALSE
    )
    if (isTRUE(ok)) added_any <- TRUE
  }

  added_any
}

.migrate_dynamic_sqlite_schema_r <- function(conn) {
  if (is.null(conn)) {
    return(FALSE)
  }
  changed <- FALSE
  for (tbl in c("vouchers", "extracted_leaves", "voucher_reproductive_phenology", "run_active_cohort", "leaf_efa_harmonics", "otu_morphometrics", "spatial_rf_predictions", "curation_triage", "environmental_cache", "genomic_samples", "genetic_lineage_assignments", "cytological_records", "taxon_cytology_summary", "taxon_backbone_profiles", "natureserve_subnational_ranges")) {
    if (isTRUE(.migrate_dynamic_sqlite_table_r(conn, tbl))) {
      changed <- TRUE
    }
  }
  changed
}
