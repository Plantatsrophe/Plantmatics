"""
Idempotent SQLite schema initialization, column upgrades, and master view refresh helpers.
"""

from __future__ import annotations

import sqlite3
from typing import Set

from scripts.core.db_schema_ddl_tier_b import SCHEMA_DDL
from scripts.core.db_schema_master_view_b import MASTER_VIEW_DDL
from scripts.core.db_schema_migration_specs import (
    PACKERA_SCHEMA_USER_VERSION,
    PLANTMATICS_SCHEMA_USER_VERSION,
    _DYNAMIC_TABLE_MIGRATION_SPECS,
    _MASTER_VIEW_REQUIRED_TABLES,
    _VOUCHERS_MIGRATION_COLUMNS,
)
from scripts.core.db_schema_step02_rebuild import (
    _rebuild_tier_a_step02_tables_if_needed,
    _tier_a_step02_needs_rebuild,
)

def _migrate_vouchers_schema_columns(conn: sqlite3.Connection) -> None:
    """Idempotently add optional exsiccatae conflict audit columns to existing `vouchers` tables."""
    existing_cols = {
        str(r[1]) for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()
    }
    if not existing_cols:
        return
    for col_name, col_ddl in _VOUCHERS_MIGRATION_COLUMNS:
        if col_name not in existing_cols:
            conn.execute(f"ALTER TABLE vouchers ADD COLUMN {col_name} {col_ddl};")
            existing_cols.add(col_name)

def _migrate_dynamic_schema_columns(
    conn: sqlite3.Connection,
    recreate_view: bool = True,
) -> None:
    """
    Idempotently upgrade existing SQLite databases in-place:
      1. Inspects `PRAGMA table_info(...)` on `leaf_efa_harmonics`, `otu_morphometrics`,
         `spatial_rf_predictions`, and `environmental_cache`.
      2. Issues `ALTER TABLE ... ADD COLUMN ...` for any missing dynamic harmonic (`A1..D24`),
         PCA (`PC1..PC15`), adaptive EFD metadata (`efd_harmonic_opt`,
         `efd_k95_pcs`, `selected_morph_features`), dissection scalars, 19-variable WorldClim
         (`bio1..bio19`), `climatic_z_score`, or `primary_*_predictors` columns.
      3. Optionally drops and recreates `v_integrative_triage_master` whenever all base tables exist.
    """
    existing_tables: Set[str] = {
        str(r[0])
        for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table';"
        ).fetchall()
    }
    if not existing_tables:
        return

    conn.execute("DROP VIEW IF EXISTS v_integrative_triage_master;")

    for table_name, col_specs in _DYNAMIC_TABLE_MIGRATION_SPECS.items():
        if table_name not in existing_tables:
            continue
        existing_cols = {
            str(r[1]) for r in conn.execute(f"PRAGMA table_info({table_name});").fetchall()
        }
        if not existing_cols:
            continue
        for col_name, col_ddl in col_specs:
            if col_name not in existing_cols:
                conn.execute(f"ALTER TABLE {table_name} ADD COLUMN {col_name} {col_ddl};")
                existing_cols.add(col_name)



    _rebuild_tier_a_step02_tables_if_needed(conn, existing_tables)

    if (
        "genomic_samples" not in existing_tables
        or "genetic_lineage_assignments" not in existing_tables
        or "cytological_records" not in existing_tables
        or "taxon_cytology_summary" not in existing_tables
        or "taxon_backbone_profiles" not in existing_tables
        or "natureserve_subnational_ranges" not in existing_tables
    ):
        conn.executescript(SCHEMA_DDL)
        # Ensure master view stays dropped until subsequent table rebuilds finish.
        conn.execute("DROP VIEW IF EXISTS v_integrative_triage_master;")
        existing_tables.update(
            {
                "genomic_samples",
                "genetic_lineage_assignments",
                "cytological_records",
                "taxon_cytology_summary",
                "taxon_backbone_profiles",
                "natureserve_subnational_ranges",
            }
        )

    if "curation_triage" in existing_tables:
        ct_cols_now = {
            str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()
        }
        for _ct_col, _ct_type in _DYNAMIC_TABLE_MIGRATION_SPECS["curation_triage"]:
            if ct_cols_now and _ct_col not in ct_cols_now:
                conn.execute(f"ALTER TABLE curation_triage ADD COLUMN {_ct_col} {_ct_type};")
                ct_cols_now.add(_ct_col)
        cur_uv_row = conn.execute("PRAGMA user_version;").fetchone()
        cur_uv = int(cur_uv_row[0]) if (cur_uv_row and cur_uv_row[0] is not None) else 0
        canon_ct_cols = (
            "run_id", "catalogNumber", "plant_individual_id",
            "w_morph_kappa", "w_edaphic_kappa", "w_pheno_kappa", "w_geno_kappa",
            "n_available_modalities", "probabilistic_consensus_conflict",
            "probabilistic_consensus_score", "hopkins_adjusted_doy",
            "delta_doy_resid", "taxonomic_status_call", "recommended_determination",
            "taxonomic_treatment_broad", "taxonomic_treatment_subtype",
            "synthesis_triage_priority", "synthesis_triage_action",
            "synthesis_rationale", "cryptic_cluster_cohesion_pass",
            "cryptic_non_clinal_pass", "resurrection_basionym",
            "hybrid_zone_cohesion_pass", "hybrid_parent_1", "hybrid_parent_2",
            "is_completed_voucher", "genomic_propagation_tier",
            "tier0_calibrated_taxon", "tier0_calibrated_prob",
            "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi",
            "cytology_linkage_tier", "voucher_sporophytic_2n",
            "voucher_ploidy_level", "taxon_reported_cytotypes",
            "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass",
            "hydro_edaphic_niche_regime", "wetland_fidelity_concordance",
            "basionym_synonymy_type", "w_morph_simplex", "w_edaphic_simplex",
            "w_pheno_simplex", "w_geno_simplex", "geno_deme_weight_scale",
            "portal_metadata_lag_flag", "portal_lag_trigger_gate",
            "is_marginal_range_extension",
        )
        if ct_cols_now and (
            ct_cols_now != set(canon_ct_cols) or cur_uv < PLANTMATICS_SCHEMA_USER_VERSION
        ):
            ct_pre_rebuild_cols = {
                str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()
            }
            shared_ct = [c for c in canon_ct_cols if c in ct_pre_rebuild_cols]
            conn.executescript(
                """
                DROP TABLE IF EXISTS _curation_triage_v20261127;
                CREATE TABLE _curation_triage_v20261127 (
                    run_id TEXT NOT NULL,
                    catalogNumber TEXT NOT NULL,
                    plant_individual_id INTEGER NOT NULL DEFAULT 0,
                    w_morph_kappa REAL,
                    w_edaphic_kappa REAL,
                    w_pheno_kappa REAL,
                    w_geno_kappa REAL,
                    n_available_modalities INTEGER,
                    probabilistic_consensus_conflict REAL,
                    probabilistic_consensus_score REAL,
                    hopkins_adjusted_doy REAL,
                    delta_doy_resid REAL,
                    taxonomic_status_call TEXT,
                    recommended_determination TEXT,
                    taxonomic_treatment_broad TEXT,
                    taxonomic_treatment_subtype TEXT,
                    synthesis_triage_priority TEXT,
                    synthesis_triage_action TEXT,
                    synthesis_rationale TEXT,
                    cryptic_cluster_cohesion_pass INTEGER,
                    cryptic_non_clinal_pass INTEGER,
                    resurrection_basionym TEXT,
                    hybrid_zone_cohesion_pass INTEGER,
                    hybrid_parent_1 TEXT,
                    hybrid_parent_2 TEXT,
                    is_completed_voucher INTEGER,
                    genomic_propagation_tier TEXT,
                    tier0_calibrated_taxon TEXT,
                    tier0_calibrated_prob REAL,
                    tier0_abstention_flag INTEGER,
                    deme_mahal_d2 REAL,
                    cluster_propagated_gdi REAL,
                    cytology_linkage_tier TEXT,
                    voucher_sporophytic_2n INTEGER,
                    voucher_ploidy_level INTEGER,
                    taxon_reported_cytotypes TEXT,
                    is_multi_cytotype_taxon INTEGER,
                    cytotype_reproductive_barrier_pass INTEGER,
                    hydro_edaphic_niche_regime TEXT,
                    wetland_fidelity_concordance REAL,
                    basionym_synonymy_type TEXT,
                    w_morph_simplex REAL,
                    w_edaphic_simplex REAL,
                    w_pheno_simplex REAL,
                    w_geno_simplex REAL,
                    geno_deme_weight_scale REAL,
                    portal_metadata_lag_flag INTEGER,
                    portal_lag_trigger_gate TEXT,
                    is_marginal_range_extension INTEGER,
                    PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
                    FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
                );
                """
            )
            if shared_ct:
                cols_csv = ", ".join(shared_ct)
                conn.execute(
                    f"INSERT OR REPLACE INTO _curation_triage_v20261127 ({cols_csv}) SELECT {cols_csv} FROM curation_triage;"
                )
            conn.executescript(
                """
                DROP TABLE curation_triage;
                ALTER TABLE _curation_triage_v20261127 RENAME TO curation_triage;
                """
            )

    if recreate_view and set(_MASTER_VIEW_REQUIRED_TABLES).issubset(existing_tables):
        conn.executescript(MASTER_VIEW_DDL)

def migrate_database_schema(conn: sqlite3.Connection) -> None:
    """Idempotently migrate existing SQLite schemas to PLANTMATICS_SCHEMA_USER_VERSION and refresh the master view."""
    _migrate_vouchers_schema_columns(conn)
    _migrate_dynamic_schema_columns(conn, recreate_view=False)
    conn.executescript(SCHEMA_DDL)
    conn.execute(f"PRAGMA user_version = {PLANTMATICS_SCHEMA_USER_VERSION};")
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.commit()

def init_db_schema(conn: sqlite3.Connection) -> None:
    """Create or idempotently migrate all Tier A and Tier B tables, indexes, triggers, R*Tree, and master view."""
    cur_ver_row = conn.execute("PRAGMA user_version;").fetchone()
    cur_ver = int(cur_ver_row[0]) if (cur_ver_row and cur_ver_row[0] is not None) else 0
    if cur_ver == PLANTMATICS_SCHEMA_USER_VERSION:
        existing_objs = {
            str(r[0])
            for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type IN ('table', 'view');"
            ).fetchall()
        }
        required_objs = set(_MASTER_VIEW_REQUIRED_TABLES) | {
            "taxon_synonyms",
            "usda_county_ranges",
            "extracted_leaves",
            "cytological_records",
            "taxon_cytology_summary",
            "taxon_backbone_profiles",
            "natureserve_subnational_ranges",
            "pipeline_runs",
            "leaf_efa_harmonics",
            "taxon_complexes",
            "v_integrative_triage_master",
        }
        if required_objs.issubset(existing_objs):
            v_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()}
            if (
                "usda_symbol" not in v_cols
                or "usda_county_status" not in v_cols
                or "is_centroid_or_inst_coord" not in v_cols
                or "coord_quarantine_reason" not in v_cols
                or "typified_basionym" not in v_cols
                or "verbatim_taxon_concept" not in v_cols
                or "concept_synonymy_status" not in v_cols
                or "identification_year" not in v_cols
                or "is_open_or_hybrid_concept" not in v_cols
            ):
                _migrate_vouchers_schema_columns(conn)
            el_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(extracted_leaves);").fetchall()}
            leh_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(leaf_efa_harmonics);").fetchall()}
            om_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(otu_morphometrics);").fetchall()}
            ec_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(environmental_cache);").fetchall()}
            srf_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(spatial_rf_predictions);").fetchall()}
            gs_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genomic_samples);").fetchall()}
            ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
            mv_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
            req_mv_cols = {"is_passive", "scale_source", "AR_inv_std", "log1p_capitulum_count", "corymb_head_count_est", "radiate_head_fraction", "blade_cie_l_norm", "blade_cie_l_sheet_delta", "blade_cie_a_delta", "blade_cie_b_delta", "basal_cordate_index", "blade_cie_l_q85_sheet_delta", "abaxial_leaf_fraction", "nearest_type_basionym", "reproductive_trait_status", "ecoregion_can1", "ecoregion_can2", "hybrid_zone_cohesion_pass", "hybrid_parent_1", "hybrid_parent_2", "cda_oof_posterior_prob", "self_morph_mahal_d2", "open_set_df", "cryptic_non_clinal_pass", "resurrection_basionym", "w_geno_kappa", "geno_predicted_taxon", "geno_linkage_tier", "is_geno_admixed", "coalescent_gdi", "hybrid_generation_class", "cytonuclear_discordant", "public_plastome_match", "public_its_match", "waxy_lineage_match", "public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index", "genomic_propagation_tier", "cluster_propagated_gdi", "cytology_linkage_tier", "voucher_sporophytic_2n", "voucher_ploidy_level", "taxon_reported_cytotypes", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass", "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct", "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m", "hydro_edaphic_niche_regime", "wetland_fidelity_concordance", "basionym_synonymy_type", "w_morph_simplex", "w_edaphic_simplex", "w_pheno_simplex", "w_geno_simplex", "geno_deme_weight_scale", "portal_metadata_lag_flag", "portal_lag_trigger_gate", "is_marginal_range_extension", "backbone_basionym_name", "protologue_url", "natureserve_global_rank", "nwpl_indicator_code", "nwpl_wetland_score", "is_centroid_or_inst_coord", "coord_quarantine_reason", "typified_basionym", "verbatim_taxon_concept", "concept_synonymy_status", "identification_year", "is_open_or_hybrid_concept", "env_mahal_d2", "env_mcd_outlier_flag", "ssurgo_lithic_series_pct", "ecotype_regime", "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx", "hopkins_adjusted_doy", "hopkins_adjusted_doy_resid", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic", "nearest_type_niche_basionym", "nearest_type_basionym_edaphic_z", "phyllary_tip_darkening_delta", "basal_to_cauline_wool_delta"}
            if (
                "leaf_zone" not in el_cols
                or "area_mm2" not in el_cols
                or "ect_path" not in el_cols
                or "mask_is_upright" not in el_cols
                or "petiole_width_px" not in el_cols
                or "petiole_width_mm" not in el_cols
                or "efa_qc_pass" not in el_cols
                or "enclosed_void_frac" not in el_cols
                or "bilateral_lobe_symmetry" not in el_cols
                or "bbox_x1" not in el_cols
                or "basal_cordate_index" not in el_cols
                or "blade_cie_l_q85_sheet_delta" not in el_cols
                or "is_abaxial_pressed_leaf" not in el_cols
                or _tier_a_step02_needs_rebuild(conn)
                or "scale_source" not in leh_cols
                or "bbox_x1" not in leh_cols
                or "lamina_trace_length_mm" not in leh_cols
                or "basal_cordate_index" not in leh_cols
                or "blade_cie_l_q85_sheet_delta" not in leh_cols
                or "is_abaxial_pressed_leaf" not in leh_cols
                or "cauline_leaf_count" not in om_cols
                or "ECT_PC1" not in om_cols
                or "blade_width_mm" not in om_cols
                or "scale_source" not in om_cols
                or "is_promoted_bronze" not in om_cols
                or "within_rosette_efa_dispersion" not in om_cols
                or "basal_cordate_index" not in om_cols
                or "blade_cie_l_q85_sheet_delta" not in om_cols
                or "abaxial_leaf_fraction" not in om_cols
                or "ssurgo_lithic_series_pct" not in ec_cols
                or "env_mahal_d2" not in srf_cols
                or "env_mcd_outlier_flag" not in srf_cols
                or "ecotype_regime" not in srf_cols
                or "lithic_outcrop_index" not in srf_cols
                or "summer_moisture_deficit_idx" not in srf_cols
                or "maritime_continentality_idx" not in srf_cols
                or "hopkins_adjusted_doy" not in srf_cols
                or "hopkins_adjusted_doy_resid" not in srf_cols
                or "anthesis_dominance_weight" not in srf_cols
                or "allochronic_bimodality_delta_bic" not in srf_cols
                or "nearest_type_niche_basionym" not in srf_cols
                or "nearest_type_basionym_edaphic_z" not in srf_cols
                or "exsiccatae_key" not in gs_cols
                or "genbank_accession" not in gs_cols
                or "coordinateUncertainty" not in gs_cols
                or "is_centroid_or_inst_coord" not in gs_cols
                or "coord_quarantine_reason" not in gs_cols
                or "typified_basionym" not in gs_cols
                or "concept_synonymy_status" not in gs_cols
                or "is_open_or_hybrid_concept" not in gs_cols
                or "w_morph_simplex" not in ct_cols
                or "w_edaphic_simplex" not in ct_cols
                or "w_pheno_simplex" not in ct_cols
                or "w_geno_simplex" not in ct_cols
                or "geno_deme_weight_scale" not in ct_cols
                or "portal_metadata_lag_flag" not in ct_cols
                or "portal_lag_trigger_gate" not in ct_cols
                or "is_marginal_range_extension" not in ct_cols
                or not req_mv_cols.issubset(mv_cols)
            ):
                _migrate_dynamic_schema_columns(conn, recreate_view=True)
            conn.execute("PRAGMA foreign_keys = ON;")
            return

    migrate_database_schema(conn)

__all__ = (
    "_migrate_vouchers_schema_columns", "_migrate_dynamic_schema_columns",
    "migrate_database_schema", "init_db_schema",
)
