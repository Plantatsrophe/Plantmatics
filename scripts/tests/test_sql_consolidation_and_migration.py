#!/usr/bin/env python3
"""
Relational contract tests for SQLite phantom pid=0 stub healing & false-split
consolidation, legacy extracted_leaves CHECK constraint rebuilding, and atomic
packera_pipeline.sqlite -> plantmatics.sqlite WAL migration across Python & R.
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from scripts.tests._r_session import rscript_available, run_warm_r_cmd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    consolidate_false_multi_plant_splits_db,
    get_db_connection,
    populate_run_active_cohort,
    purge_segments_and_runs,
    record_voucher_segmentation_transaction,
    upsert_vouchers_records,
)


class TestSQLConsolidationAndMigration(unittest.TestCase):
    """Integration tests for DB split consolidation, schema rebuilds, and legacy migration."""

    def test_phantom_pid0_stub_cleanup_and_false_split_db_consolidation(self) -> None:
        """Verify SQLite heals phantom pid=0 stubs and cauline/1-leaf splits in place without GPU re-segmentation while preserving true multi-rosette sheets."""
        with tempfile.TemporaryDirectory() as tmp_d:
            db_file = Path(tmp_d) / "heal_splits.sqlite"
            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                v_rows = [
                    {"catalogNumber": "CAT-PHANTOM", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "CAT-CAUL-SPLIT", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver"},
                    {"catalogNumber": "CAT-TRUE-MULTI", "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
                ]
                upsert_vouchers_records(conn, v_rows)
                with conn:
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('CAT-PHANTOM', 0);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-PHANTOM', 1, 900.0, 2400.0);")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-PHANTOM', 1, 1, 'CAT-PHANTOM_p1_leaf1', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-PHANTOM', 1, 2, 'CAT-PHANTOM_p1_leaf2', '[[0,0],[1,0],[0,1]]', 'basal');")

                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-CAUL-SPLIT', 0, 950.0, 2500.0);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y, plant_height_cm) VALUES ('CAT-CAUL-SPLIT', 1, 960.0, 1100.0, 42.5);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-CAUL-SPLIT', 2, 1350.0, 2480.0);")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 0, 1, 'CAT-CAUL-SPLIT_p0_leaf1', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 0, 2, 'CAT-CAUL-SPLIT_p0_leaf2', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 1, 3, 'CAT-CAUL-SPLIT_p1_leaf3', '[[0,0],[1,0],[0,1]]', 'cauline');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 1, 4, 'CAT-CAUL-SPLIT_p1_leaf4', '[[0,0],[1,0],[0,1]]', 'cauline');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 2, 5, 'CAT-CAUL-SPLIT_p2_leaf5', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count, total_capitula, n_anthesis, voucher_phenological_state, is_verified_anthesis) VALUES ('CAT-CAUL-SPLIT', 0, 0, 0, 0, 'sterile', 0);")
                    conn.execute("INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count, total_capitula, n_anthesis, voucher_phenological_state, is_verified_anthesis) VALUES ('CAT-CAUL-SPLIT', 1, 6, 6, 5, 'anthesis', 1);")

                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-TRUE-MULTI', 0, 350.0, 2400.0);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-TRUE-MULTI', 1, 1650.0, 2420.0);")
                    for pid_m in (0, 1):
                        for lid_m in (1, 2):
                            eff_l = pid_m * 10 + lid_m
                            conn.execute(
                                "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES (?, ?, ?, ?, '[[0,0],[1,0],[0,1]]', 'basal');",
                                ("CAT-TRUE-MULTI", pid_m, eff_l, f"CAT-TRUE-MULTI_p{pid_m}_leaf{eff_l}"),
                            )

                n_coh = populate_run_active_cohort(conn, run_id="heal_run", voucher_rows=v_rows, clear_existing=True)
                self.assertEqual(n_coh, 4)

                phantom_pids = [r[0] for r in conn.execute("SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = 'CAT-PHANTOM';").fetchall()]
                self.assertEqual(phantom_pids, [1])

                split_pids = [r[0] for r in conn.execute("SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = 'CAT-CAUL-SPLIT';").fetchall()]
                self.assertEqual(split_pids, [0])
                split_leaves = conn.execute("SELECT COUNT(*) FROM extracted_leaves WHERE catalogNumber = 'CAT-CAUL-SPLIT' AND plant_individual_id = 0;").fetchone()[0]
                self.assertEqual(int(split_leaves), 5)
                vrp_merged = conn.execute("SELECT capitula_count, n_anthesis, voucher_phenological_state, is_verified_anthesis FROM voucher_reproductive_phenology WHERE catalogNumber = 'CAT-CAUL-SPLIT' AND plant_individual_id = 0;").fetchone()
                self.assertEqual((int(vrp_merged[0]), int(vrp_merged[1]), str(vrp_merged[2]), int(vrp_merged[3])), (6, 5, "anthesis", 1))
                h_cm = conn.execute("SELECT plant_height_cm FROM plant_individuals WHERE catalogNumber = 'CAT-CAUL-SPLIT' AND plant_individual_id = 0;").fetchone()[0]
                self.assertAlmostEqual(float(h_cm), 42.5)

                true_multi_pids = sorted(r[0] for r in conn.execute("SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = 'CAT-TRUE-MULTI';").fetchall())
                self.assertEqual(true_multi_pids, [0, 1])

                stats_again = consolidate_false_multi_plant_splits_db(conn)
                self.assertEqual(stats_again, {"phantom_zero_stubs_removed": 0, "false_split_individuals_merged": 0})
            finally:
                conn.close()

    def test_legacy_extracted_leaves_check_constraint_rebuilt_on_init_and_purge(self) -> None:
        """Verify pre-LM3 extracted_leaves CHECK constraints and retired NOT NULL columns are stripped by init_db_schema, savepoint transactions, and purge_segments_and_runs."""
        leg_col = "_".join(("assigned", "tier"))
        leg_v1 = "_".join(("Tier", "1", "Direct"))
        leg_v2 = "_".join(("Tier", "2", "Reflected"))

        def _seed_legacy_db(db_file: Path) -> None:
            init_conn = get_db_connection(db_file, init_schema=True)
            init_conn.close()
            raw = sqlite3.connect(db_file)
            raw.executescript(
                f"""
                PRAGMA foreign_keys = OFF;
                PRAGMA user_version = 20261123;
                INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, qc_status) VALUES
                    ('WILLI_82182', 'Packera dubia', 'Packera dubia', 'PASS'),
                    ('GA202225', 'Packera anonyma', 'Packera anonyma', 'PENDING'),
                    ('VSC0067189', 'Packera aurea', 'Packera aurea', 'PENDING');
                INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('WILLI_82182', 0);
                DROP TABLE IF EXISTS extracted_leaves;
                CREATE TABLE extracted_leaves (
                    catalogNumber TEXT NOT NULL,
                    plant_individual_id INTEGER NOT NULL DEFAULT 0,
                    leaf_id INTEGER NOT NULL,
                    shape_id TEXT UNIQUE,
                    {leg_col} TEXT NOT NULL CHECK({leg_col} IN ('{leg_v1}', '{leg_v2}')),
                    solidity REAL,
                    contour_xy_json TEXT NOT NULL,
                    PRIMARY KEY (catalogNumber, plant_individual_id, leaf_id),
                    FOREIGN KEY (catalogNumber, plant_individual_id)
                        REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
                );
                INSERT INTO extracted_leaves (
                    catalogNumber, plant_individual_id, leaf_id, shape_id, {leg_col}, solidity, contour_xy_json
                ) VALUES (
                    'WILLI_82182', 0, 1, 'WILLI_82182_p0_leaf1', '{leg_v1}', 0.84, '[[0,0],[10,0],[5,8]]'
                );
                """
            )
            raw.commit()
            raw.close()

        with tempfile.TemporaryDirectory() as tmp:
            db_a = Path(tmp) / "inplace" / "data" / "plantmatics.sqlite"
            db_a.parent.mkdir(parents=True, exist_ok=True)
            _seed_legacy_db(db_a)

            conn_a = get_db_connection(db_a, init_schema=True)
            try:
                el_cols_a = {str(r[1]) for r in conn_a.execute("PRAGMA table_info(extracted_leaves);").fetchall()}
                self.assertNotIn(leg_col, el_cols_a)
                self.assertIn("efa_qc_pass", el_cols_a)
                kept = conn_a.execute("SELECT solidity FROM extracted_leaves WHERE catalogNumber = 'WILLI_82182';").fetchone()
                self.assertIsNotNone(kept)
                self.assertAlmostEqual(float(kept[0]), 0.84, places=3)

                with conn_a:
                    for idx, cat in enumerate(("WILLI_82182", "GA202225", "VSC0067189")):
                        sp = f"sp_v_{idx}"
                        conn_a.execute(f"SAVEPOINT {sp};")
                        record_voucher_segmentation_transaction(
                            conn=conn_a,
                            catalog_number=cat,
                            extracted_records=[{
                                "catalogNumber": cat, "plant_individual_id": 0, "leaf_id": 1,
                                "solidity": 0.79, "area_px": 1200.0, "pixels_per_mm": 11.5,
                                "is_calibrated_mm": 1, "scale_source": "ruler_calibrated",
                                "leaf_zone": "basal", "efa_qc_pass": 1, "qc_flag_reason": "PASS",
                                "contour_xy_json": "[[0,0],[12,0],[12,9],[0,9]]",
                            }],
                            commit=False,
                        )
                        conn_a.execute(f"RELEASE SAVEPOINT {sp};")
                self.assertEqual(int(conn_a.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0]), 3)
            finally:
                conn_a.close()

            ws_b = Path(tmp) / "purge_ws"
            (ws_b / "data" / "tables").mkdir(parents=True, exist_ok=True)
            (ws_b / "data" / "tables" / "curated_vouchers.csv").write_text("catalogNumber\nWILLI_82182\n", encoding="utf-8")
            db_b = ws_b / "data" / "plantmatics.sqlite"
            _seed_legacy_db(db_b)

            summary = purge_segments_and_runs(workspace_root=ws_b, db_path=db_b, purge_runs=True, vacuum=True)
            self.assertEqual(summary["purged_db_rows"], 2)

            conn_b = get_db_connection(db_b, init_schema=False)
            try:
                el_cols_b = {str(r[1]) for r in conn_b.execute("PRAGMA table_info(extracted_leaves);").fetchall()}
                self.assertNotIn(leg_col, el_cols_b)
                ddl_b = str(conn_b.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='extracted_leaves';").fetchone()[0])
                self.assertNotIn("CHECK", ddl_b.upper())
                record_voucher_segmentation_transaction(
                    conn=conn_b,
                    catalog_number="GA202225",
                    extracted_records=[{"catalogNumber": "GA202225", "plant_individual_id": 0, "leaf_id": 1, "solidity": 0.81, "contour_xy_json": "[[0,0],[10,0],[5,7]]"}],
                    commit=True,
                )
                self.assertEqual(int(conn_b.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0]), 1)
            finally:
                conn_b.close()

    def test_legacy_packera_pipeline_sqlite_auto_migrates_to_plantmatics(self) -> None:
        """Verify legacy data/packera_pipeline.sqlite (+ WAL/SHM sidecars) atomically migrates to plantmatics.sqlite."""

        def _seed_legacy_wal_workspace(ws: Path, wal_catalog: str, base_catalog: str) -> tuple[Path, Path]:
            data_dir = ws / "data"
            scratch = ws / "_legacy_seed_scratch"
            data_dir.mkdir(parents=True, exist_ok=True)
            scratch.mkdir(parents=True, exist_ok=True)

            src_db = scratch / "packera_pipeline.sqlite"
            conn = get_db_connection(db_path=src_db, init_schema=True, workspace_root=ws)
            try:
                conn.executescript(
                    """
                    DROP VIEW IF EXISTS v_integrative_triage_master;
                    DROP TABLE IF EXISTS genomic_samples;
                    DROP TABLE IF EXISTS genetic_lineage_assignments;
                    CREATE TABLE genomic_samples (
                        genetic_sample_id TEXT PRIMARY KEY, catalogNumber TEXT,
                        linkage_tier TEXT, species_standardized TEXT,
                        geno_quality_pass INTEGER, mapped_reads INTEGER
                    );
                    CREATE TABLE genetic_lineage_assignments (
                        run_id TEXT NOT NULL, genetic_sample_id TEXT NOT NULL,
                        catalogNumber TEXT, geno_predicted_taxon TEXT,
                        PRIMARY KEY (run_id, genetic_sample_id)
                    );
                    """
                )
                conn.execute(
                    "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier) VALUES (?, ?, ?, ?);",
                    (base_catalog, "Packera dubia", "Packera dubia", "Tier_1_Gold"),
                )
                conn.execute("PRAGMA user_version = 20261123;")
                conn.commit()
            finally:
                conn.close()

            wal_conn = sqlite3.connect(str(src_db))
            try:
                wal_conn.execute("PRAGMA journal_mode = WAL;")
                wal_conn.execute("PRAGMA wal_autocheckpoint = 0;")
                wal_conn.execute(
                    "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier) VALUES (?, ?, ?, ?);",
                    (wal_catalog, "Packera schweinitziana", "Packera schweinitziana", "Tier_1_Gold"),
                )
                wal_conn.commit()
                for suffix in ("", "-wal", "-shm"):
                    src_f = Path(f"{src_db}{suffix}")
                    self.assertTrue(src_f.is_file(), f"Missing seeded WAL sidecar: {src_f}")
                    if suffix:
                        self.assertGreater(src_f.stat().st_size, 0, f"Empty WAL sidecar: {src_f}")
                    shutil.copy2(src_f, data_dir / f"packera_pipeline.sqlite{suffix}")
            finally:
                wal_conn.close()
                shutil.rmtree(scratch, ignore_errors=True)

            legacy_db = data_dir / "packera_pipeline.sqlite"
            target_db = data_dir / "plantmatics.sqlite"
            self.assertTrue(legacy_db.is_file())
            self.assertFalse(target_db.exists())
            return legacy_db, target_db

        new_gs_cols = {"public_plastome_match", "public_its_match", "waxy_lineage_match"}
        new_gla_cols = {"public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index"}
        new_ct_prop_cols = {
            "genomic_propagation_tier", "tier0_calibrated_taxon", "tier0_calibrated_prob",
            "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi",
        }

        # 1. Python db*.py atomic migration + WAL replay + 20261127 schema upgrade
        with tempfile.TemporaryDirectory() as td_py:
            ws_py = Path(td_py)
            legacy_db, target_db = _seed_legacy_wal_workspace(ws_py, wal_catalog="NCU_LEGACY_WAL_02", base_catalog="NCU_LEGACY_01")
            env_clean = os.environ.copy()
            for k in ("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH", "PACKERA_SQLITE_DB", "PACKERA_DEMO_MODE"):
                env_clean.pop(k, None)
            with mock.patch.dict(os.environ, env_clean, clear=True):
                conn = get_db_connection(workspace_root=ws_py, init_schema=True)
                try:
                    self.assertFalse(legacy_db.exists())
                    self.assertFalse(Path(f"{legacy_db}-wal").exists())
                    self.assertFalse(Path(f"{legacy_db}-shm").exists())
                    self.assertTrue(target_db.is_file())
                    self.assertEqual(int(conn.execute("PRAGMA user_version;").fetchone()[0]), 20261127)

                    gs_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genomic_samples);").fetchall()}
                    gla_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genetic_lineage_assignments);").fetchall()}
                    ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
                    view_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
                    self.assertTrue(new_gs_cols.issubset(gs_cols))
                    self.assertTrue(new_gla_cols.issubset(gla_cols))
                    self.assertTrue(new_ct_prop_cols.issubset(ct_cols))
                    self.assertTrue((new_gs_cols | new_gla_cols | new_ct_prop_cols).issubset(view_cols))
                    v_cats = {str(r[0]) for r in conn.execute("SELECT catalogNumber FROM vouchers;").fetchall()}
                    self.assertEqual(v_cats, {"NCU_LEGACY_01", "NCU_LEGACY_WAL_02"})
                finally:
                    conn.close()

        # 2. R db_utils.R (and sourced db*.R helpers) resolve_packera_db_path() + connect_packera_db() migration
        if rscript_available():
            with tempfile.TemporaryDirectory() as td_r:
                ws_r = Path(td_r)
                legacy_r, target_r = _seed_legacy_wal_workspace(ws_r, wal_catalog="NCU_R_WAL_02", base_catalog="NCU_R_LEGACY_01")
                r_code = f"""
source({json.dumps(str(PROJECT_ROOT / "scripts" / "core" / "db_utils.R"))})
local({{
  resolved_p <- resolve_packera_db_path()
  cat("RESOLVED_DB=", resolved_p, "\\n", sep = "")
  stopifnot(grepl("plantmatics\\\\.sqlite$", resolved_p))
  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  on.exit(disconnect_packera_db(conn), add = TRUE)
  stopifnot(!is.null(conn))
  uv <- as.integer(DBI::dbGetQuery(conn, "PRAGMA user_version;")[[1L]][1L])
  stopifnot(uv == 20261127L)
  gs_cols <- DBI::dbListFields(conn, "genomic_samples")
  gla_cols <- DBI::dbListFields(conn, "genetic_lineage_assignments")
  ct_cols <- DBI::dbListFields(conn, "curation_triage")
  mv_cols <- DBI::dbListFields(conn, "v_integrative_triage_master")
  req_cols <- c(
    "public_plastome_match", "public_its_match", "waxy_lineage_match",
    "public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index",
    "genomic_propagation_tier", "tier0_calibrated_taxon", "tier0_calibrated_prob",
    "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi"
  )
  stopifnot(all(c("public_plastome_match", "public_its_match", "waxy_lineage_match") %in% gs_cols))
  stopifnot(all(c("public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index") %in% gla_cols))
  stopifnot(all(c(
    "genomic_propagation_tier", "tier0_calibrated_taxon", "tier0_calibrated_prob",
    "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi"
  ) %in% ct_cols))
  stopifnot(all(req_cols %in% mv_cols))
  vs <- DBI::dbGetQuery(conn, "SELECT catalogNumber FROM vouchers ORDER BY catalogNumber;")$catalogNumber
  stopifnot(all(c("NCU_R_LEGACY_01", "NCU_R_WAL_02") %in% vs))
}})
"""
                proc = run_warm_r_cmd(["Rscript", "-e", r_code], stage_key="step03", cwd=ws_r, env={"PYTHONPATH": str(PROJECT_ROOT)})
                self.assertEqual(proc.returncode, 0, f"R legacy DB migration failed:\n{proc.stdout}\n{proc.stderr}")
                self.assertFalse(legacy_r.exists())
                self.assertFalse(Path(f"{legacy_r}-wal").exists())
                self.assertFalse(Path(f"{legacy_r}-shm").exists())
                self.assertTrue(target_r.is_file())


if __name__ == "__main__":
    unittest.main()
