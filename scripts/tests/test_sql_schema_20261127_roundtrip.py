#!/usr/bin/env python3
"""
Relational contract test for consolidated 20261118 -> 20261127 SQLite schema
migration, backbone/geochemistry/wetland/cytology/triage round-trip through
v_integrative_triage_master, and Tier A preservation across purge_runs_only.
"""

from __future__ import annotations

import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    ensure_pipeline_run,
    get_db_connection,
    persist_otu_morphometrics,
    populate_run_active_cohort,
    purge_runs_only,
    upsert_rows_generic,
    upsert_vouchers_records,
)
from scripts.core.db_schema import (
    PACKERA_SCHEMA_USER_VERSION,
    PLANTMATICS_SCHEMA_USER_VERSION,
)


class TestSQLSchema20261127Roundtrip(unittest.TestCase):
    """Integration test verifying 20261127 schema migration and master view roundtrip."""

    def test_backbone_geochem_wetland_schema_20261127_roundtrip(self) -> None:
        """Verify 20261118 -> 20261127 consolidated migration creates/migrates all schema tables/columns, roundtrips through v_integrative_triage_master, and preserves Tier A tables across purge_runs_only."""
        self.assertEqual(PLANTMATICS_SCHEMA_USER_VERSION, 20261127)
        self.assertEqual(PACKERA_SCHEMA_USER_VERSION, 20261127)

        expected_om_cols = (
            "cda_oof_posterior_prob", "cl_target_threshold",
            "self_morph_mahal_d2", "self_morph_df", "open_set_df",
        )
        expected_backbone_cols = {
            "species_standardized", "target_genus", "accepted_scientific_name",
            "accepted_authorship", "basionym_name", "basionym_authorship",
            "homotypic_synonyms_json", "heterotypic_synonyms_json", "is_hybrid_taxon",
            "hybrid_formula", "publication_citation", "publication_year",
            "protologue_url", "locus_classicus_text", "tdwg_native_codes_json",
            "tdwg_introduced_codes_json", "natureserve_element_id",
            "natureserve_global_rank", "nwpl_indicator_code", "nwpl_wetland_score",
            "backbone_source", "updated_at",
        }
        expected_ns_cols = {
            "species_standardized", "country_code", "subnation_code",
            "s_rank", "is_native", "is_exotic", "is_historical_or_extirpated",
        }
        expected_env_cols = (
            "soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
            "is_nwi_wetland", "nwi_wetland_system", "nwi_wetland_dist_m",
        )
        expected_ct_cols = (
            "hybrid_zone_cohesion_pass", "hybrid_parent_1", "hybrid_parent_2",
            "cryptic_non_clinal_pass", "resurrection_basionym", "genomic_propagation_tier",
            "tier0_calibrated_taxon", "tier0_calibrated_prob", "tier0_abstention_flag",
            "deme_mahal_d2", "cluster_propagated_gdi", "cytology_linkage_tier",
            "voucher_sporophytic_2n", "voucher_ploidy_level", "taxon_reported_cytotypes",
            "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass",
            "hydro_edaphic_niche_regime", "wetland_fidelity_concordance", "basionym_synonymy_type",
        )
        expected_mv_aliases = (
            "backbone_basionym_name", "protologue_url", "natureserve_global_rank",
            "nwpl_indicator_code", "nwpl_wetland_score",
        )

        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp)
            db_p = ws / "data" / "plantmatics.sqlite"
            db_p.parent.mkdir(parents=True, exist_ok=True)
            raw_conn = sqlite3.connect(db_p)
            raw_conn.executescript(
                """
                PRAGMA user_version = 20261118;
                CREATE TABLE pipeline_runs (
                    run_id TEXT PRIMARY KEY, run_dir TEXT, start_timestamp TEXT,
                    end_timestamp TEXT, git_commit_sha TEXT, cli_arguments_json TEXT,
                    focal_taxon TEXT, intraspecific_mode INTEGER DEFAULT 0
                );
                INSERT INTO pipeline_runs (run_id, start_timestamp) VALUES ('r_27', '2026-11-14T00:00:00Z');
                CREATE TABLE environmental_cache (lat_lon_key TEXT PRIMARY KEY, catalogNumber TEXT);
                INSERT INTO environmental_cache (lat_lon_key, catalogNumber) VALUES ('35.50000_-82.50000', 'NCU-OLD-26');
                CREATE TABLE curation_triage (
                    run_id TEXT NOT NULL, catalogNumber TEXT NOT NULL,
                    plant_individual_id INTEGER NOT NULL DEFAULT 0,
                    taxonomic_status_call TEXT, taxonomic_treatment_subtype TEXT,
                    recommended_determination TEXT, cryptic_cluster_cohesion_pass INTEGER,
                    is_completed_voucher INTEGER,
                    PRIMARY KEY (run_id, catalogNumber, plant_individual_id)
                );
                INSERT INTO curation_triage VALUES (
                    'r_27', 'NCU-OLD-26', 0, 'Species_Confirmed', 'Confirmed_Species', 'Packera aurea', 0, 1
                );
                """
            )
            raw_conn.commit()
            raw_conn.close()

            conn = get_db_connection(db_p, init_schema=True)
            try:
                self.assertEqual(int(conn.execute("PRAGMA user_version;").fetchone()[0]), 20261127)
                tables = {str(r[0]) for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()}
                for req_tbl in ("taxon_backbone_profiles", "natureserve_subnational_ranges", "cytological_records", "taxon_cytology_summary"):
                    self.assertIn(req_tbl, tables)

                tbp_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(taxon_backbone_profiles);").fetchall()}
                ns_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(natureserve_subnational_ranges);").fetchall()}
                om_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(otu_morphometrics);").fetchall()}
                env_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(environmental_cache);").fetchall()}
                ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
                mv_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}

                self.assertEqual(len(expected_backbone_cols), 22)
                self.assertTrue(expected_backbone_cols.issubset(tbp_cols), f"Missing backbone cols: {expected_backbone_cols - tbp_cols}")
                self.assertTrue(expected_ns_cols.issubset(ns_cols), f"Missing natureserve cols: {expected_ns_cols - ns_cols}")
                for col in (*expected_om_cols, *expected_env_cols, *expected_ct_cols, *expected_mv_aliases):
                    self.assertIn(col, mv_cols, f"v_integrative_triage_master missing {col}")
                for col in expected_om_cols:
                    self.assertIn(col, om_cols, f"otu_morphometrics missing {col}")
                for col in expected_env_cols:
                    self.assertIn(col, env_cols, f"environmental_cache missing {col}")
                for col in expected_ct_cols:
                    self.assertIn(col, ct_cols, f"curation_triage missing {col}")

                old_ct = conn.execute(
                    "SELECT recommended_determination FROM curation_triage WHERE run_id='r_27' AND catalogNumber='NCU-OLD-26';"
                ).fetchone()
                self.assertEqual(old_ct[0], "Packera aurea")

                upsert_rows_generic(
                    conn,
                    table_name="taxon_backbone_profiles",
                    rows=[{
                        "species_standardized": "Packera serpenticola", "target_genus": "Packera",
                        "accepted_scientific_name": "Packera serpenticola",
                        "accepted_authorship": "Boufford, Kartesz, Shi & R.Zhou",
                        "basionym_name": "Packera serpenticola",
                        "basionym_authorship": "Boufford, Kartesz, Shi & R.Zhou",
                        "homotypic_synonyms_json": "[]", "heterotypic_synonyms_json": "[]",
                        "is_hybrid_taxon": 0, "hybrid_formula": "",
                        "publication_citation": "Syst. Bot. 39(3): 1027", "publication_year": 2014,
                        "protologue_url": "https://www.biodiversitylibrary.org/page/123456",
                        "locus_classicus_text": "Buck Creek Serpentine Barrens, Clay Co., NC",
                        "tdwg_native_codes_json": '["NCA"]', "tdwg_introduced_codes_json": "[]",
                        "natureserve_element_id": "ELEMENT_GLOBAL.2.999999",
                        "natureserve_global_rank": "G1", "nwpl_indicator_code": "UPL",
                        "nwpl_wetland_score": 5.0, "backbone_source": "COL_WCVP_2232+BHL+NatureServe",
                        "updated_at": "2026-11-14T00:00:00Z",
                    }],
                    pk_cols=("species_standardized",),
                )
                upsert_rows_generic(
                    conn,
                    table_name="natureserve_subnational_ranges",
                    rows=[{
                        "species_standardized": "Packera serpenticola", "country_code": "US",
                        "subnation_code": "NC", "s_rank": "S1", "is_native": 1,
                        "is_exotic": 0, "is_historical_or_extirpated": 0,
                    }],
                    pk_cols=("species_standardized", "country_code", "subnation_code"),
                )
                upsert_rows_generic(
                    conn,
                    table_name="cytological_records",
                    rows=[{
                        "cytology_id": "CYTO-001", "source_db": "CCDB",
                        "original_taxon_name": "Packera serpenticola", "basionym_name": "Packera serpenticola",
                        "species_standardized": "Packera serpenticola", "catalogNumber": "NCU-BB-27",
                        "cytology_linkage_tier": "C1_Direct_Cytovoucher",
                        "sporophytic_2n": 44, "base_x": 22, "ploidy_level": 2,
                    }],
                    pk_cols=("cytology_id",),
                )
                upsert_rows_generic(
                    conn,
                    table_name="taxon_cytology_summary",
                    rows=[{
                        "species_standardized": "Packera serpenticola", "n_cytology_records": 1,
                        "n_direct_cytovouchers": 1, "modal_sporophytic_2n": 44, "modal_base_x": 22,
                        "modal_ploidy_level": 2, "reported_cytotypes": "2x", "is_multi_cytotype_taxon": 0,
                    }],
                    pk_cols=("species_standardized",),
                )
                vrows = [{
                    "catalogNumber": "NCU-BB-27", "institutionCode": "NCU",
                    "scientificName": "Packera serpenticola", "species_standardized": "Packera serpenticola",
                    "decimalLatitude": 35.08333, "decimalLongitude": -83.61667,
                    "determiner_tier": "Tier_1_Gold", "qc_status": "PASS",
                }]
                upsert_vouchers_records(conn, vrows)
                ensure_pipeline_run(conn, run_id="r_27")
                populate_run_active_cohort(conn, run_id="r_27", voucher_rows=vrows, clear_existing=False)
                persist_otu_morphometrics(
                    conn,
                    run_id="r_27",
                    otu_rows=[{
                        "catalogNumber": "NCU-BB-27", "plant_individual_id": 0,
                        "species_standardized": "Packera serpenticola", "determiner_tier": "Tier_1_Gold",
                        "cda_predicted_taxon": "Packera serpenticola", "cda_posterior_prob": 0.92,
                        "cda_oof_posterior_prob": 0.89, "cl_target_threshold": 0.74,
                        "self_morph_mahal_d2": 3.25, "self_morph_df": 4,
                        "min_congener_mahal_d2": 29.4, "min_congener_taxon": "Packera obovata",
                        "open_set_novelty_p": 0.65, "open_set_df": 14,
                        "nearest_type_basionym": "Packera serpenticola", "nearest_type_basionym_d2": 1.15,
                    }],
                )
                upsert_rows_generic(
                    conn,
                    table_name="environmental_cache",
                    rows=[{
                        "lat_lon_key": "35.08333_-83.61667", "catalogNumber": "NCU-BB-27",
                        "soil_ca_mg_molar_ratio": 0.42, "soil_ultramafic_ni_cr_index": 2.85,
                        "soil_carbonate_c_pct": 0.05, "is_nwi_wetland": 0,
                        "nwi_wetland_system": "Upland", "nwi_wetland_dist_m": 650.0,
                    }],
                    pk_cols=("lat_lon_key",),
                )
                upsert_rows_generic(
                    conn,
                    table_name="curation_triage",
                    rows=[{
                        "run_id": "r_27", "catalogNumber": "NCU-BB-27", "plant_individual_id": 0,
                        "taxonomic_status_call": "Species_Confirmed", "recommended_determination": "Packera serpenticola",
                        "hybrid_zone_cohesion_pass": 0, "hybrid_parent_1": None, "hybrid_parent_2": None,
                        "cryptic_non_clinal_pass": 0, "resurrection_basionym": "Senecio earlei",
                        "genomic_propagation_tier": "Tier0_Projected_Calibrated",
                        "tier0_calibrated_taxon": "Packera serpenticola", "tier0_calibrated_prob": 0.94,
                        "tier0_abstention_flag": 0, "deme_mahal_d2": 1.85, "cluster_propagated_gdi": 0.84,
                        "cytology_linkage_tier": "C1_Direct_Cytovoucher", "voucher_sporophytic_2n": 44,
                        "voucher_ploidy_level": 2, "taxon_reported_cytotypes": "2x, 4x",
                        "is_multi_cytotype_taxon": 1, "cytotype_reproductive_barrier_pass": 1,
                        "hydro_edaphic_niche_regime": "Ultramafic_Serpentine_Barren",
                        "wetland_fidelity_concordance": 1.0, "basionym_synonymy_type": "Heterotypic",
                    }],
                    pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
                )
                mv_row = conn.execute(
                    """
                    SELECT backbone_basionym_name, protologue_url, natureserve_global_rank,
                           nwpl_indicator_code, nwpl_wetland_score,
                           soil_ca_mg_molar_ratio, soil_ultramafic_ni_cr_index, soil_carbonate_c_pct,
                           is_nwi_wetland, nwi_wetland_system, nwi_wetland_dist_m,
                           cda_oof_posterior_prob, cl_target_threshold, self_morph_mahal_d2, self_morph_df, open_set_df,
                           hybrid_zone_cohesion_pass, cryptic_non_clinal_pass, resurrection_basionym,
                           genomic_propagation_tier, tier0_calibrated_taxon, tier0_calibrated_prob,
                           tier0_abstention_flag, deme_mahal_d2, cluster_propagated_gdi,
                           cytology_linkage_tier, voucher_sporophytic_2n, voucher_ploidy_level,
                           taxon_reported_cytotypes, is_multi_cytotype_taxon, cytotype_reproductive_barrier_pass,
                           hydro_edaphic_niche_regime, wetland_fidelity_concordance, basionym_synonymy_type
                    FROM v_integrative_triage_master
                    WHERE run_id = 'r_27' AND catalogNumber = 'NCU-BB-27';
                    """
                ).fetchone()
                self.assertEqual(
                    tuple(mv_row),
                    (
                        "Packera serpenticola", "https://www.biodiversitylibrary.org/page/123456",
                        "G1", "UPL", 5.0, 0.42, 2.85, 0.05, 0, "Upland", 650.0,
                        0.89, 0.74, 3.25, 4, 14, 0, 0, "Senecio earlei",
                        "Tier0_Projected_Calibrated", "Packera serpenticola", 0.94, 0, 1.85, 0.84,
                        "C1_Direct_Cytovoucher", 44, 2, "2x, 4x", 1, 1,
                        "Ultramafic_Serpentine_Barren", 1.0, "Heterotypic",
                    ),
                )
            finally:
                conn.close()

            purge_runs_only(workspace_root=ws, db_path=db_p, vacuum=False)
            conn2 = get_db_connection(db_p, init_schema=False)
            try:
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM taxon_backbone_profiles;").fetchone()[0], 1)
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM natureserve_subnational_ranges;").fetchone()[0], 1)
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM cytological_records;").fetchone()[0], 1)
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM taxon_cytology_summary;").fetchone()[0], 1)
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM environmental_cache;").fetchone()[0], 2)
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM otu_morphometrics;").fetchone()[0], 0)
                self.assertEqual(conn2.execute("SELECT COUNT(*) FROM curation_triage;").fetchone()[0], 0)
            finally:
                conn2.close()


if __name__ == "__main__":
    unittest.main()
