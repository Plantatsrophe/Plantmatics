#!/usr/bin/env python3
"""
Relational contract & cross-language handshake tests for SQLite WAL DDL,
FOREIGN KEY cascades, R*Tree spatial index, DuckDB attachment, --demo
isolation, and Step 03 EFA execution from SQLite contours with zero CSV reads.
"""

from __future__ import annotations

import csv
import importlib.util
import math
import os
import shutil
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Tuple

from scripts.tests._r_session import rscript_available, run_warm_r_cmd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    DEFAULT_SQLITE_REL_PATH,
    get_db_connection,
    get_duckdb_connection,
    ingest_existing_csvs_to_db,
    query_environmental_bbox,
    resolve_db_path,
    upsert_environmental_cache_records,
    upsert_vouchers_records,
)


def _make_closed_leaf_contour(
    n_pts: int = 120,
    rx: float = 95.0,
    ry: float = 58.0,
    taper: float = 0.14,
    serration: float = 1.8,
) -> List[Tuple[float, float]]:
    """Generate a smooth closed 120-vertex botanical lamina contour."""
    pts: List[Tuple[float, float]] = []
    for idx in range(n_pts):
        theta = 2.0 * math.pi * (idx / float(n_pts))
        r_mod = 1.0 + (serration / 100.0) * math.cos(8.0 * theta)
        x = 220.0 + rx * r_mod * math.cos(theta)
        width_scale = 1.0 + taper * math.cos(theta)
        y = 220.0 + ry * r_mod * width_scale * math.sin(theta)
        pts.append((round(x, 4), round(y, 4)))
    return pts


class TestSQLDDLWALAndContours(unittest.TestCase):
    """Integration tests verifying SQLite WAL DDL, R*Tree triggers, and zero-CSV Step 03 EFA."""

    def test_ddl_schema_wal_foreign_keys_and_rtree_triggers(self) -> None:
        """Verify WAL mode, Tier A/B DDL, FK cascade rules, R*Tree spatial index, DuckDB, and --demo isolation."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_dir = Path(tmp_d)
            db_file = tmp_dir / "packera_test.sqlite"

            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                journal_mode = str(conn.execute("PRAGMA journal_mode;").fetchone()[0]).lower()
                fk_enabled = int(conn.execute("PRAGMA foreign_keys;").fetchone()[0])
                self.assertEqual(journal_mode, "wal")
                self.assertEqual(fk_enabled, 1)

                obj_rows = conn.execute(
                    "SELECT name, type FROM sqlite_master WHERE type IN ('table', 'view');"
                ).fetchall()
                existing_objects = {str(r[0]): str(r[1]) for r in obj_rows}

                expected_tables = {
                    "taxon_synonyms", "vouchers", "plant_individuals", "extracted_leaves",
                    "voucher_reproductive_phenology", "environmental_cache", "environmental_rtree",
                    "pipeline_runs", "run_active_cohort", "leaf_efa_harmonics",
                    "otu_morphometrics", "spatial_rf_predictions", "curation_triage",
                }
                for tbl in expected_tables:
                    self.assertIn(tbl, existing_objects, f"Missing SQLite table: {tbl}")
                self.assertEqual(existing_objects.get("v_integrative_triage_master"), "view")

                with self.assertRaises(sqlite3.IntegrityError):
                    with conn:
                        conn.execute(
                            "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('ORPHAN-1', 0);"
                        )
                with self.assertRaises(sqlite3.IntegrityError):
                    with conn:
                        conn.execute(
                            "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id) VALUES ('ORPHAN-1', 0, 1);"
                        )

                upsert_vouchers_records(
                    conn,
                    [{
                        "catalogNumber": "NCU-FK-01", "institutionCode": "NCU",
                        "scientificName": "Packera aurea", "species_standardized": "Packera aurea",
                        "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.9132, "decimalLongitude": -79.0558,
                    }],
                )
                with conn:
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU-FK-01', 0);")
                    conn.execute(
                        "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, contour_xy_json, contour_n_pts) "
                        "VALUES ('NCU-FK-01', 0, 1, '[[0,0],[1,0],[1,1]]', 3);"
                    )
                with conn:
                    conn.execute("DELETE FROM vouchers WHERE catalogNumber = 'NCU-FK-01';")
                rem_pi = conn.execute("SELECT COUNT(*) FROM plant_individuals WHERE catalogNumber = 'NCU-FK-01';").fetchone()[0]
                rem_el = conn.execute("SELECT COUNT(*) FROM extracted_leaves WHERE catalogNumber = 'NCU-FK-01';").fetchone()[0]
                self.assertEqual(int(rem_pi), 0)
                self.assertEqual(int(rem_el), 0)

                upsert_vouchers_records(
                    conn,
                    [
                        {"catalogNumber": "NCU-GEO-NC", "scientificName": "Packera aurea"},
                        {"catalogNumber": "NCU-GEO-ME", "scientificName": "Packera paupercula"},
                    ],
                )
                upsert_environmental_cache_records(
                    conn,
                    [
                        {"catalogNumber": "NCU-GEO-NC", "latitude": 35.91320, "longitude": -79.05580, "fine_soil_ph": 5.6, "sgmc_lithology_class": "Felsic_Crystalline"},
                        {"catalogNumber": "NCU-GEO-ME", "latitude": 44.80160, "longitude": -68.77120, "fine_soil_ph": 6.2, "sgmc_lithology_class": "Carbonate"},
                    ],
                    demo=False,
                )
                nc_hits = query_environmental_bbox(conn, min_lon=-80.5, max_lon=-78.0, min_lat=35.0, max_lat=36.5)
                self.assertEqual(len(nc_hits), 1)
                self.assertEqual(nc_hits[0]["catalogNumber"], "NCU-GEO-NC")
                self.assertAlmostEqual(float(nc_hits[0]["fine_soil_ph"]), 5.6, places=3)

                n_demo_env = upsert_environmental_cache_records(
                    conn,
                    [{"catalogNumber": "DEMO-1", "latitude": 30.0, "longitude": -80.0, "fine_soil_ph": 7.0}],
                    demo=True,
                )
                self.assertEqual(n_demo_env, 0)
                canonical_resolved = (PROJECT_ROOT / DEFAULT_SQLITE_REL_PATH).resolve()
                demo_resolved = resolve_db_path(db_path=canonical_resolved, workspace_root=PROJECT_ROOT, demo=True)
                self.assertNotEqual(str(demo_resolved), str(canonical_resolved))
            finally:
                conn.close()

            if rscript_available():
                r_expr = (
                    "source('scripts/core/db_utils.R'); "
                    f"df <- query_environmental_cache_sqlite(db_path = '{db_file.as_posix()}', "
                    "bbox = c(-80.5, -78.0, 35.0, 36.5)); "
                    "cat('R_BBOX_N=', nrow(df), 'CAT=', as.character(df$catalogNumber[1]), sep='')"
                )
                r_proc = run_warm_r_cmd(["Rscript", "-e", r_expr], stage_key="step03", cwd=PROJECT_ROOT)
                self.assertEqual(r_proc.returncode, 0, msg=r_proc.stderr)
                self.assertIn("R_BBOX_N=1CAT=NCU-GEO-NC", r_proc.stdout)

            if importlib.util.find_spec("duckdb") is not None:
                dconn = get_duckdb_connection(db_path=db_file, read_only=True)
                try:
                    cnt = dconn.execute("SELECT COUNT(*) FROM vouchers;").fetchone()[0]
                    self.assertEqual(int(cnt), 2)
                finally:
                    dconn.close()

    def test_contour_migration_and_step03_efa_zero_csv_reads(self) -> None:
        """Migrate CSVs + 120-vertex contours to SQLite, delete contour CSVs, and run Step 03 EFA purely from SQLite."""
        if not rscript_available():
            self.skipTest("Rscript not installed or fast mode enabled.")
        rscript_bin = "Rscript"

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d) / "ws"
            tables_dir = ws / "data" / "tables"
            contours_dir = ws / "data" / "contours"
            reports_dir = ws / "outputs" / "reports"
            figures_dir = ws / "outputs" / "figures"
            for d in (tables_dir, contours_dir, reports_dir, figures_dir):
                d.mkdir(parents=True, exist_ok=True)

            db_path = ws / "data" / "plantmatics.sqlite"
            run_id = "test_efa_sql_zero_csv"
            voucher_rows: List[Dict[str, Any]] = []
            manifest_rows: List[Dict[str, Any]] = []

            specs: List[Tuple[str, str, List[int]]] = [
                ("NCU-A1", "Packera aurea", [0, 1]),
                ("NCU-A2", "Packera aurea", [0]),
                ("NCU-A3", "Packera aurea", [0]),
                ("NCU-A4", "Packera aurea", [0]),
                ("NCU-A5", "Packera aurea", [0]),
                ("NCU-B1", "Packera obovata", [0]),
                ("NCU-B2", "Packera obovata", [0]),
                ("NCU-B3", "Packera obovata", [0]),
                ("NCU-B4", "Packera obovata", [0]),
                ("NCU-B5", "Packera obovata", [0]),
            ]

            for v_idx, (cat, taxon, pids) in enumerate(specs):
                voucher_rows.append({
                    "catalogNumber": cat, "institutionCode": "NCU" if v_idx % 2 == 0 else "US",
                    "scientificName": taxon, "species_raw": taxon, "species_standardized": taxon,
                    "identifiedBy": "T.M. Barkley", "determiner_tier": "Tier_1_Gold", "type_status": "",
                    "county": "Orange", "stateProvince": "North Carolina",
                    "decimalLatitude": round(35.5 + 0.05 * v_idx, 5), "decimalLongitude": round(-79.5 + 0.05 * v_idx, 5),
                    "coordinateUncertainty": 250.0, "year": 2019, "month": 4, "day": 18, "eventDate": "2019-04-18",
                    "regional_group": "Piedmont", "image_path": f"data/raw_vouchers/{cat}.jpg",
                    "exsiccatae_key": f"barkley_{v_idx}", "is_primary_duplicate": "True", "duplicate_count": 1,
                })
                is_aurea = taxon == "Packera aurea"
                for pid in pids:
                    for local_lid in (1, 2):
                        leaf_id = pid * 10 + local_lid
                        c_file = contours_dir / f"{cat}_p{pid}_leaf{leaf_id}.csv"
                        rx = (92.0 + 2.0 * local_lid + 4.0 * pid) if is_aurea else (108.0 + 2.0 * local_lid)
                        ry = (74.0 + 1.5 * local_lid - 3.0 * pid) if is_aurea else (48.0 + 1.5 * local_lid)
                        taper = 0.08 if is_aurea else -0.18
                        pts = _make_closed_leaf_contour(120, rx=rx, ry=ry, taper=taper, serration=1.2)
                        with open(c_file, "w", encoding="utf-8", newline="") as cf:
                            w = csv.writer(cf)
                            w.writerow(["x", "y"])
                            for px, py in pts:
                                w.writerow([px, py])

                        manifest_rows.append({
                            "catalogNumber": cat, "plant_individual_id": pid, "leaf_id": leaf_id,
                            "solidity": 0.91, "circularity": 0.68 if is_aurea else 0.54,
                            "extent": 0.73, "archetype_score": 0.94, "area_px": 18500.0,
                            "pixels_per_mm": 10.0, "is_calibrated_mm": "True",
                            "petiole_length_px": 280.0 if is_aurea else 95.0,
                            "petiole_length_mm": (28.0 + 5.0 * pid) if is_aurea else 9.5,
                            "blade_length_mm": 32.0 if is_aurea else 41.0,
                            "petiole_to_blade_ratio": round(((28.0 + 5.0 * pid) / 32.0) if is_aurea else (9.5 / 41.0), 4),
                            "is_dissected": "False", "contour_path": str(c_file),
                        })

            vouchers_csv = tables_dir / "curated_vouchers.csv"
            with open(vouchers_csv, "w", encoding="utf-8", newline="") as vf:
                vw = csv.DictWriter(vf, fieldnames=list(voucher_rows[0].keys()))
                vw.writeheader()
                vw.writerows(voucher_rows)

            manifest_csv = tables_dir / "extracted_leaf_manifest.csv"
            with open(manifest_csv, "w", encoding="utf-8", newline="") as mf:
                mw = csv.DictWriter(mf, fieldnames=list(manifest_rows[0].keys()))
                mw.writeheader()
                mw.writerows(manifest_rows)

            mig_counts = ingest_existing_csvs_to_db(
                workspace_root=ws, run_id=run_id, db_path=db_path,
                tables_dir=tables_dir, contours_dir=contours_dir,
            )
            self.assertEqual(mig_counts["vouchers"], 10)
            self.assertEqual(mig_counts["plant_individuals"], 11)
            self.assertEqual(mig_counts["extracted_leaves"], 22)
            self.assertEqual(mig_counts["run_active_cohort"], 11)

            shutil.rmtree(contours_dir)
            contours_dir.mkdir(parents=True, exist_ok=True)
            self.assertEqual(list(contours_dir.glob("*.csv")), [])

            out_otu_csv = tables_dir / "leaf_efa_harmonics.csv"
            out_ind_csv = tables_dir / "leaf_efa_harmonics_individual.csv"
            env = os.environ.copy()
            env["PLANTMATICS_DB_PATH"] = str(db_path)
            env["RUN_ID"] = run_id
            env["PLANTMATICS_RUN_DIR"] = str(ws)

            cmd = [
                rscript_bin, str(PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R"),
                "--input", str(contours_dir), "--contours-dir", str(contours_dir),
                "--manifest", str(manifest_csv), "--vouchers", str(vouchers_csv),
                "--output", str(out_otu_csv), "--output-individual", str(out_ind_csv),
            ]
            proc = run_warm_r_cmd(cmd, stage_key="step03", cwd=PROJECT_ROOT, env=env, skip_figures=True)
            combined_log = proc.stdout + "\n" + proc.stderr
            self.assertEqual(proc.returncode, 0, f"Step 03 failed with empty contours_dir:\n{combined_log}")
            self.assertIn("from SQLite extracted_leaves", combined_log)
            self.assertTrue(out_otu_csv.is_file())

            with open(out_otu_csv, "r", encoding="utf-8", newline="") as f:
                otu_rows = list(csv.DictReader(f))
            self.assertEqual(len(otu_rows), 11)
            otu_keys = {(r["catalogNumber"], int(float(r["plant_individual_id"]))) for r in otu_rows}
            self.assertIn(("NCU-A1", 0), otu_keys)
            self.assertIn(("NCU-A1", 1), otu_keys)
            for r in otu_rows:
                self.assertEqual(r["morph_tier_source"], "Basal_Rosette")
                self.assertTrue(math.isfinite(float(r["PC1"])))

            conn = get_db_connection(db_path=db_path, read_only=True)
            try:
                n_leaf_db = conn.execute("SELECT COUNT(*) FROM leaf_efa_harmonics WHERE run_id = ?;", (run_id,)).fetchone()[0]
                n_otu_db = conn.execute("SELECT COUNT(*) FROM otu_morphometrics WHERE run_id = ?;", (run_id,)).fetchone()[0]
                self.assertEqual(int(n_leaf_db), 22)
                self.assertEqual(int(n_otu_db), 11)
            finally:
                conn.close()


if __name__ == "__main__":
    unittest.main()
