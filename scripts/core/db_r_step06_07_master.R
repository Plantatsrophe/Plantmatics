if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/core/db_r_step06_07_master.R
# Steps 06 and 07 master view queries (v_integrative_triage_master) and SQLite persistence.
# ==============================================================================

fetch_run_master_view <- function(conn, run_id = resolve_packera_run_id()) {
  if (is.null(conn)) {
    return(NULL)
  }
  rid <- resolve_packera_run_id(run_id)
  res <- tryCatch(
    DBI::dbGetQuery(
      conn,
      "SELECT * FROM v_integrative_triage_master WHERE run_id = $1 ORDER BY catalogNumber, plant_individual_id;",
      params = list(rid)
    ),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0L) {
    return(NULL)
  }
  res$run_id <- NULL
  coerce_sqlite_logical_cols(res)
}

try_fetch_master_view_for_cohort <- function(
  input_paths = character(0),
  active_catalogs = NULL,
  run_id = NULL,
  demo = FALSE,
  vouchers_df = NULL,
  vouchers_path = NULL,
  require_step04 = FALSE,
  require_step05 = FALSE,
  require_step07 = FALSE,
  ...
) {
  eff_paths <- unique(c(input_paths, vouchers_path))
  eff_paths <- eff_paths[!is.na(eff_paths) & nzchar(eff_paths)]
  if (!is_sqlite_active_for_inputs(eff_paths, demo = demo)) {
    return(NULL)
  }
  conn <- connect_packera_db(read_only = TRUE, demo = demo)
  if (is.null(conn)) {
    return(NULL)
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  master_df <- fetch_run_master_view(conn, run_id = rid)
  if (is.null(master_df) || nrow(master_df) == 0L) {
    return(NULL)
  }

  if (is.null(active_catalogs) && !is.null(vouchers_df) && "catalogNumber" %in% names(vouchers_df)) {
    active_catalogs <- vouchers_df$catalogNumber
  }
  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) > 0L) {
      if (!all(cats %in% trimws(as.character(master_df$catalogNumber)))) {
        return(NULL)
      }
      master_df <- master_df[trimws(as.character(master_df$catalogNumber)) %in% cats, , drop = FALSE]
    }
  }
  if (nrow(master_df) == 0L) {
    return(NULL)
  }

  has_m <- ("cda_predicted_taxon" %in% names(master_df) && any(!is.na(master_df$cda_predicted_taxon))) ||
    ("PC1" %in% names(master_df) && any(!is.na(master_df$PC1))) ||
    ("ECT_PC1" %in% names(master_df) && any(is.finite(suppressWarnings(as.numeric(master_df$ECT_PC1))))) ||
    ("can1" %in% names(master_df) && any(is.finite(suppressWarnings(as.numeric(master_df$can1)))))
  has_s <- ("rf_predicted_taxon" %in% names(master_df) && any(!is.na(master_df$rf_predicted_taxon))) ||
    ("edaphic_best_fit_taxon" %in% names(master_df) && any(!is.na(master_df$edaphic_best_fit_taxon)))

  if (isTRUE(require_step04) && !has_m) {
    return(NULL)
  }
  if ((isTRUE(require_step05) || isTRUE(require_step07)) && !has_s) {
    return(NULL)
  }

  if (!isTRUE(require_step07)) {
    stale_ct_cols <- c(
      "w_morph_kappa", "w_edaphic_kappa", "w_pheno_kappa", "w_geno_kappa",
      "w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex",
      "geno_deme_weight_scale", "portal_metadata_lag_flag", "portal_lag_trigger_gate",
      "is_marginal_range_extension", "n_available_modalities",
      "probabilistic_consensus_conflict", "probabilistic_consensus_score",
      "delta_doy_resid", "taxonomic_status_call", "recommended_determination",
      "taxonomic_treatment_broad", "taxonomic_treatment_subtype",
      "synthesis_triage_priority", "synthesis_triage_action", "synthesis_rationale",
      "cryptic_cluster_cohesion_pass", "cryptic_non_clinal_pass", "resurrection_basionym",
      "hybrid_zone_cohesion_pass", "hybrid_parent_1", "hybrid_parent_2",
      "is_completed_voucher",
      "genomic_propagation_tier", "tier0_calibrated_taxon", "tier0_calibrated_prob",
      "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi",
      "cytology_linkage_tier", "voucher_sporophytic_2n", "voucher_ploidy_level",
      "taxon_reported_cytotypes", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass",
      "hydro_edaphic_niche_regime", "wetland_fidelity_concordance", "basionym_synonymy_type"
    )
    # Strip stale ct.delta_doy_resid while explicitly preserving Step 05 hopkins_adjusted_doy_resid
    stale_ct_cols <- setdiff(stale_ct_cols, "hopkins_adjusted_doy_resid")
    master_df <- master_df[, setdiff(names(master_df), stale_ct_cols), drop = FALSE]
  }
  if (!isTRUE(require_step05) && !isTRUE(require_step07)) {
    stale_srf_cols <- c(
      "ecotype_regime", "lithic_outcrop_index", "summer_moisture_deficit_idx",
      "maritime_continentality_idx", "rf_predicted_taxon", "rf_predicted_prob",
      "rf_given_taxon_prob", "rf_prediction_source",
      "edaphic_best_fit_taxon", "spatial_block_cv_acc", "spatial_autocorr_range_km",
      "spatial_block_id", "expected_doy", "hopkins_expected", "hopkins_adjusted_doy",
      "hopkins_adjusted_doy_resid", "anthesis_dominance_weight",
      "allochronic_bimodality_delta_bic", "delta_doy", "pheno_timing_category",
      "anthesis_filter_status", "edaphic_z_score",
      "climatic_z_score", "type_edaphic_z_score", "nearest_type_niche_taxon",
      "nearest_type_niche_basionym", "nearest_type_basionym_edaphic_z",
      "matches_type_lithology", "primary_env_predictors", "primary_clim_predictors",
      "primary_edaph_predictors", "substrate_specialist_flag", "n_observed_modalities",
      "multimodal_concordance", "multimodal_conflict_flag", "triage_category",
      "conflict_discordance_notes", "env_mahal_d2", "env_mcd_outlier_flag"
    )
    master_df <- master_df[, setdiff(names(master_df), stale_srf_cols), drop = FALSE]
  }
  master_df
}

persist_step07_to_sqlite <- function(
  triage_df,
  vouchers_path = NULL,
  run_id = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(triage_df) || nrow(triage_df) == 0L) {
    return(invisible(FALSE))
  }
  if (!is_sqlite_active_for_inputs(vouchers_path, demo = demo)) {
    return(invisible(FALSE))
  }

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) {
    return(invisible(FALSE))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  ct_df <- as.data.frame(triage_df, stringsAsFactors = FALSE)
  ct_df$run_id <- rid
  if (!("plant_individual_id" %in% names(ct_df))) ct_df$plant_individual_id <- 0L
  ct_df$plant_individual_id[is.na(ct_df$plant_individual_id)] <- 0L
  for (.bcol in intersect(c("cryptic_cluster_cohesion_pass", "cryptic_non_clinal_pass", "hybrid_zone_cohesion_pass", "is_completed_voucher"), names(ct_df))) {
    ct_df[[.bcol]] <- ifelse(is.na(ct_df[[.bcol]]), NA_integer_, as.integer(as.logical(ct_df[[.bcol]])))
  }
  for (.fcol in intersect(c("tier0_abstention_flag", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass", "portal_metadata_lag_flag", "is_marginal_range_extension"), names(ct_df))) {
    v_abs <- ct_df[[.fcol]]
    if (is.logical(v_abs)) {
      ct_df[[.fcol]] <- as.integer(v_abs)
    } else if (is.character(v_abs) || is.factor(v_abs)) {
      s_abs <- trimws(as.character(v_abs))
      out_abs <- suppressWarnings(as.integer(s_abs))
      u_abs <- toupper(s_abs)
      out_abs[is.na(out_abs) & u_abs %in% c("TRUE", "T", "YES", "Y", "1.0")] <- 1L
      out_abs[is.na(out_abs) & u_abs %in% c("FALSE", "F", "NO", "N", "0.0")] <- 0L
      ct_df[[.fcol]] <- out_abs
    } else {
      v_num <- suppressWarnings(as.numeric(v_abs))
      out_abs <- rep(NA_integer_, length(v_num))
      ok_abs <- !is.na(v_num) & is.finite(v_num)
      out_abs[ok_abs] <- as.integer(round(v_num[ok_abs]))
      ct_df[[.fcol]] <- out_abs
    }
  }
  for (.icol in intersect(c("voucher_sporophytic_2n", "voucher_ploidy_level"), names(ct_df))) {
    v_i <- ct_df[[.icol]]
    v_num <- if (is.factor(v_i)) suppressWarnings(as.numeric(as.character(v_i))) else suppressWarnings(as.numeric(v_i))
    out_i <- rep(NA_integer_, length(v_num))
    ok_i <- !is.na(v_num) & is.finite(v_num)
    out_i[ok_i] <- as.integer(round(v_num[ok_i]))
    ct_df[[.icol]] <- out_i
  }
  for (.fcol in intersect(c("w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex", "geno_deme_weight_scale", "tier0_calibrated_prob", "deme_mahal_d2", "cluster_propagated_gdi", "wetland_fidelity_concordance"), names(ct_df))) {
    v_f <- ct_df[[.fcol]]
    v_num <- if (is.factor(v_f)) suppressWarnings(as.numeric(as.character(v_f))) else suppressWarnings(as.numeric(v_f))
    v_num[!is.finite(v_num)] <- NA_real_
    ct_df[[.fcol]] <- v_num
  }
  for (.scol in intersect(c("genomic_propagation_tier", "tier0_calibrated_taxon", "cytology_linkage_tier", "taxon_reported_cytotypes", "hydro_edaphic_niche_regime", "basionym_synonymy_type", "portal_lag_trigger_gate"), names(ct_df))) {
    ct_df[[.scol]] <- as.character(ct_df[[.scol]])
  }

  db_upsert_dataframe(
    conn,
    table_name = "curation_triage",
    df = ct_df,
    pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
    preserve_existing_on_null = FALSE
  )
  .prune_stale_run_rows_r(conn, "curation_triage", rid, ct_df)
  invisible(TRUE)
}

persist_step07_complexes_to_sqlite <- function(
  complexes_df,
  vouchers_path = NULL,
  run_id = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(complexes_df) || nrow(complexes_df) == 0L) {
    return(invisible(FALSE))
  }
  if (!is_sqlite_active_for_inputs(vouchers_path, demo = demo)) {
    return(invisible(FALSE))
  }

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) {
    return(invisible(FALSE))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  tc_df <- as.data.frame(complexes_df, stringsAsFactors = FALSE)
  tc_df$run_id <- rid
  if ("integrative_congruence" %in% names(tc_df)) {
    tc_df$integrative_congruence <- ifelse(is.na(tc_df$integrative_congruence), NA_integer_, as.integer(as.logical(tc_df$integrative_congruence)))
  }

  db_upsert_dataframe(
    conn,
    table_name = "taxon_complexes",
    df = tc_df,
    pk_cols = c("run_id", "species_standardized"),
    preserve_existing_on_null = FALSE
  )
  DBI::dbExecute(
    conn,
    sprintf(
      "DELETE FROM taxon_complexes WHERE run_id = %s AND species_standardized NOT IN (%s);",
      as.character(DBI::dbQuoteString(conn, rid)),
      paste(vapply(as.character(tc_df$species_standardized), function(s) as.character(DBI::dbQuoteString(conn, s)), character(1)), collapse = ", ")
    )
  )
  invisible(TRUE)
}

persist_step06_genomics_to_sqlite <- function(
  assignments_df,
  run_id = NULL,
  db_path = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(assignments_df) || nrow(assignments_df) == 0L) {
    return(invisible(FALSE))
  }
  conn <- connect_packera_db(db_path = db_path, read_only = FALSE, demo = FALSE)
  if (is.null(conn)) {
    return(invisible(FALSE))
  }
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  gla_df <- as.data.frame(assignments_df, stringsAsFactors = FALSE)
  gla_df$run_id <- rid
  if ("is_geno_admixed" %in% names(gla_df)) {
    gla_df$is_geno_admixed <- ifelse(
      is.na(gla_df$is_geno_admixed),
      NA_integer_,
      as.integer(as.logical(gla_df$is_geno_admixed))
    )
  }
  if ("public_backbone_clade" %in% names(gla_df)) {
    gla_df$public_backbone_clade <- as.character(gla_df$public_backbone_clade)
  }
  for (.ncol in intersect(c("public_phytop_ils_index", "public_phytop_ih_index"), names(gla_df))) {
    gla_df[[.ncol]] <- suppressWarnings(as.numeric(gla_df[[.ncol]]))
  }

  db_upsert_dataframe(
    conn,
    table_name = "genetic_lineage_assignments",
    df = gla_df,
    pk_cols = c("run_id", "genetic_sample_id"),
    preserve_existing_on_null = FALSE
  )
  kept_ids <- as.character(gla_df$genetic_sample_id[!is.na(gla_df$genetic_sample_id) & nzchar(trimws(as.character(gla_df$genetic_sample_id)))])
  if (length(kept_ids) > 0L) {
    DBI::dbExecute(
      conn,
      sprintf(
        "DELETE FROM genetic_lineage_assignments WHERE run_id = %s AND genetic_sample_id NOT IN (%s);",
        as.character(DBI::dbQuoteString(conn, rid)),
        paste(vapply(kept_ids, function(s) as.character(DBI::dbQuoteString(conn, s)), character(1)), collapse = ", ")
      )
    )
  }

  if (db_table_exists_safe(conn, "genomic_samples") && ("genetic_sample_id" %in% names(gla_df))) {
    gs_sync_cols <- intersect(
      c(
        "cytonuclear_discordant",
        "public_plastome_match",
        "public_its_match",
        "waxy_lineage_match",
        "exsiccatae_key",
        "genbank_accession",
        "coordinateUncertainty",
        "is_centroid_or_inst_coord",
        "coord_quarantine_reason",
        "typified_basionym",
        "concept_synonymy_status",
        "is_open_or_hybrid_concept"
      ),
      names(gla_df)
    )
    if (length(gs_sync_cols) > 0L) {
      gs_id_raw <- trimws(as.character(gla_df$genetic_sample_id))
      keep_gs <- !is.na(gs_id_raw) & nzchar(gs_id_raw)
      if (any(keep_gs)) {
        gs_df <- gla_df[keep_gs, c("genetic_sample_id", gs_sync_cols), drop = FALSE]
        gs_df$genetic_sample_id <- gs_id_raw[keep_gs]
        for (.bcol in intersect(
          c("cytonuclear_discordant", "is_centroid_or_inst_coord", "is_open_or_hybrid_concept"),
          names(gs_df)
        )) {
          .bv <- gs_df[[.bcol]]
          if (is.logical(.bv)) {
            gs_df[[.bcol]] <- ifelse(is.na(.bv), NA_integer_, ifelse(.bv, 1L, 0L))
          } else if (is.numeric(.bv) || is.integer(.bv)) {
            gs_df[[.bcol]] <- ifelse(is.na(.bv) | is.nan(.bv) | is.infinite(.bv), NA_integer_, ifelse(.bv != 0, 1L, 0L))
          } else {
            .bs <- toupper(trimws(as.character(.bv)))
            .bo <- rep(NA_integer_, length(.bs))
            .bo[.bs %in% c("TRUE", "T", "1", "1.0", "YES", "Y")] <- 1L
            .bo[.bs %in% c("FALSE", "F", "0", "0.0", "NO", "N")] <- 0L
            gs_df[[.bcol]] <- .bo
          }
        }
        if ("coordinateUncertainty" %in% names(gs_df)) {
          gs_df$coordinateUncertainty <- suppressWarnings(as.numeric(gs_df$coordinateUncertainty))
        }
        for (.scol in intersect(
          c(
            "public_plastome_match",
            "public_its_match",
            "waxy_lineage_match",
            "exsiccatae_key",
            "genbank_accession",
            "coord_quarantine_reason",
            "typified_basionym",
            "concept_synonymy_status"
          ),
          names(gs_df)
        )) {
          gs_df[[.scol]] <- as.character(gs_df[[.scol]])
        }
        db_upsert_dataframe(
          conn,
          table_name = "genomic_samples",
          df = gs_df,
          pk_cols = c("genetic_sample_id"),
          preserve_existing_on_null = TRUE
        )
      }
    }
  }
  invisible(TRUE)
}

persist_step06_to_sqlite <- persist_step06_genomics_to_sqlite
