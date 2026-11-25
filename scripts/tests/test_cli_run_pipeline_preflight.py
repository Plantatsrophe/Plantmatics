#!/usr/bin/env python3
"""
Unit tests for run_pipeline.sh preflight cohort filtering:
  - --gold-only and --determiner-tier stratification
  - SQLite curated_vouchers.csv reconciliation and missing-taxon fail-fast
  - --exclude-misidentified and --triage-source (CSV and pure-SQLite) filtering
"""

import os
import re
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _load_pipeline_sh_corpus() -> str:
    """Concatenate run_pipeline.sh and any scripts/pipeline/sh/*.sh helper files."""
    parts = [(PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")]
    sh_dir = PROJECT_ROOT / "scripts" / "pipeline" / "sh"
    if sh_dir.is_dir():
        for p in sorted(sh_dir.glob("*.sh")):
            parts.append(p.read_text(encoding="utf-8"))
    return "\n".join(parts)


def _load_preflight_subset_filter_body() -> str:
    """Load run_preflight_subset_filter and its helpers from sh/ modules or run_pipeline.sh."""
    sh_dir = PROJECT_ROOT / "scripts" / "pipeline" / "sh"
    parts = [p for p in (sh_dir / "pipe_preflight_cohort_b.sh", sh_dir / "pipe_preflight_cohort_a.sh") if p.is_file()]
    if parts:
        cleaned = [
            re.sub(r"(?m)^\s*(?:source|\.)\s+.*pipe_preflight_cohort_[ab]\.sh.*\n?", "", p.read_text(encoding="utf-8"))
            for p in parts
        ]
        return "\n".join(cleaned)
    sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
    fn_start = sh_src.index("run_preflight_subset_filter() {")
    fn_end = sh_src.index("\nrun_preflight_subset_filter\n", fn_start)
    return sh_src[fn_start:fn_end]


def _find_preflight_and_step02_positions() -> tuple[int, int]:
    """Locate preflight call and Step 02 invocation positions in run_pipeline.sh."""
    sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
    pf_m = re.search(r"(?m)^\s*run_preflight_subset_filter\b(?!\s*\(\))", sh_src)
    if not pf_m:
        pf_m = re.search(r"pipe_preflight_cohort_a\.sh", sh_src)
    preflight_pos = pf_m.start() if pf_m else -1

    step02_pos = sh_src.find("02_segment_and_extract.py")
    if step02_pos < 0:
        sh_dir = PROJECT_ROOT / "scripts" / "pipeline" / "sh"
        if sh_dir.is_dir():
            for p in sorted(sh_dir.glob("*.sh")):
                txt = p.read_text(encoding="utf-8")
                idx = txt.find("02_segment_and_extract.py")
                if idx < 0:
                    continue
                fn_matches = list(re.finditer(r"(?m)^\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(\)\s*\{", txt[:idx]))
                if fn_matches:
                    fn_name = fn_matches[-1].group(1)
                    call_m = re.search(rf"(?m)^\s*{re.escape(fn_name)}\b(?!\s*\(\))", sh_src)
                    if call_m:
                        step02_pos = call_m.start()
                        break
                src_m = re.search(re.escape(p.name), sh_src)
                if src_m:
                    step02_pos = src_m.start()
                    break
    return preflight_pos, step02_pos


def _run_preflight_bash(
    run_dir: Path,
    root_data: Path,
    bash_vars: dict[str, str],
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Execute run_preflight_subset_filter in an isolated bash subshell against root_data."""
    fn_body = _load_preflight_subset_filter_body()
    lines = [
        "set -euo pipefail",
        f'PROJECT_ROOT="{PROJECT_ROOT}"',
        f'RUN_DIR="{run_dir}"',
        f'PYTHON_MAIN="{sys.executable}"',
    ]
    for k, v in bash_vars.items():
        lines.append(f'{k}="{v}"')
    lines.append(fn_body)
    script = "\n".join(lines)
    script = script.replace('"${PROJECT_ROOT}/data"', f'"{root_data}"')
    script = script.replace("'${PROJECT_ROOT}/data'", f"'{root_data}'")
    script += "\nrun_preflight_subset_filter\n"
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("RUN_ID", "PLANTMATICS_DB_PATH", "PLANTMATICS_DEMO_MODE", "PACKERA_DB_PATH", "PACKERA_DEMO_MODE")
    }
    return subprocess.run(["bash", "-c", script], cwd=str(PROJECT_ROOT), capture_output=True, text=True, env=env, check=check)


class TestCLIRunPipelinePreflight(unittest.TestCase):
    """Test suite for run_pipeline.sh preflight cohort filtering and SQLite sync."""

    def test_run_pipeline_sh_gold_only_and_determiner_tier_filtering(self):
        """Verify run_pipeline.sh --gold-only and --determiner-tier filter curated_vouchers.csv, Step 02 manifests, and cropped_patches strictly to Tier_1_Gold before Steps 02-07."""
        sh_corpus = _load_pipeline_sh_corpus()
        self.assertIn("--gold-only)", sh_corpus)
        self.assertIn("--determiner-tier)", sh_corpus)
        self.assertIn("--determiner-tier=*)", sh_corpus)
        preflight_pos, step02_pos = _find_preflight_and_step02_positions()
        self.assertGreaterEqual(preflight_pos, 0, "Could not locate run_preflight_subset_filter call")
        self.assertGreaterEqual(step02_pos, 0, "Could not locate Step 02 segmentation call")
        self.assertLess(preflight_pos, step02_pos)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            root_data = tmp_path / "root_data"
            run_dir = tmp_path / "runs" / "test_gold_only"
            tables_dir = run_dir / "data" / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)
            (root_data / "cropped_patches" / "capitula").mkdir(parents=True, exist_ok=True)

            pd.DataFrame({
                "catalogNumber": ["NCU001", "UNCC_49308", "US002", "NY003"],
                "scientificName": ["Packera dubia", "Packera anonyma", "Packera paupercula", "Packera plattensis"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"],
                "authority_weight": [3.0, 3.0, 1.5, 0.5],
            }).to_csv(tables_dir / "curated_vouchers.csv", index=False)

            core_tbls = (
                "extracted_leaf_manifest.csv",
                "failed_qc_vouchers.csv",
                "voucher_reproductive_metrics.csv",
                "voucher_phenological_states.csv",
            )
            for tbl_name in core_tbls:
                pd.DataFrame({
                    "catalogNumber": ["NCU001", "UNCC_49308", "US002", "NY003"],
                    "plant_individual_id": [0, 0, 0, 0],
                }).to_csv(tables_dir / tbl_name, index=False)

            (root_data / "cropped_patches" / "capitula" / "NCU001_cap1.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "capitula" / "US002_cap1.jpg").write_bytes(b"jpg")

            res = _run_preflight_bash(
                run_dir,
                root_data,
                {
                    "FOCAL_TAXON": "",
                    "EXCLUDE_TAXON_CSV": "",
                    "INTRASPECIFIC_MODE": "false",
                    "GOLD_ANCHORS": "30",
                    "GOLD_ONLY": "true",
                    "DETERMINER_TIER": "Tier_1_Gold",
                    "REUSE_CANONICAL_DATA": "true",
                },
                check=True,
            )
            self.assertIn("Tier_1_Gold", res.stdout)

            out_vouchers = pd.read_csv(tables_dir / "curated_vouchers.csv")
            self.assertEqual(set(out_vouchers["determiner_tier"].unique()), {"Tier_1_Gold"})
            self.assertEqual(set(out_vouchers["catalogNumber"]), {"NCU001", "UNCC_49308"})
            for tbl_name in core_tbls:
                tdf = pd.read_csv(tables_dir / tbl_name)
                self.assertEqual(set(tdf["catalogNumber"]), {"NCU001", "UNCC_49308"}, f"Failed filter on {tbl_name}")

            link_p = run_dir / "data" / "cropped_patches" / "capitula"
            self.assertTrue(link_p.is_symlink(), "Expected directory-level symlink for capitula")
            self.assertEqual(link_p.resolve(), (root_data / "cropped_patches" / "capitula").resolve())

            db_file = root_data / "plantmatics.sqlite"
            self.assertTrue(db_file.is_file())
            with sqlite3.connect(db_file) as conn:
                cohort_rows = conn.execute(
                    "SELECT catalogNumber, determiner_tier FROM run_active_cohort WHERE run_id = ?;",
                    (run_dir.name,),
                ).fetchall()
            self.assertEqual({r[0] for r in cohort_rows}, {"NCU001", "UNCC_49308"})
            self.assertEqual({r[1] for r in cohort_rows}, {"Tier_1_Gold"})

    def test_reconcile_curated_vouchers_csv_with_db_restores_requested_taxa(self):
        """Verify SQLite vouchers missing from curated_vouchers.csv are appended (taxon-restricted), preserving existing rows byte-for-byte."""
        from scripts.core.db import get_db_connection, reconcile_curated_vouchers_csv_with_db, upsert_vouchers_records

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            csv_path = tmp_path / "curated_vouchers.csv"
            original = (
                "catalogNumber,scientificName,determiner_tier,is_primary_duplicate,duplicate_count\n"
                "A1,Packera anonyma,Tier_1_Gold,True,1\n"
            )
            csv_path.write_text(original, encoding="utf-8")

            conn = get_db_connection(db_path=tmp_path / "test.sqlite", init_schema=True)
            try:
                upsert_vouchers_records(conn, [
                    {"catalogNumber": "A1", "scientificName": "Packera anonyma", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "O1", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold",
                     "species_standardized": "Packera obovata", "is_primary_duplicate": 1, "duplicate_count": 2},
                    {"catalogNumber": "P1", "scientificName": "Packera paupercula", "determiner_tier": "Tier_3_Bronze",
                     "species_standardized": "Packera paupercula"},
                ])
                self.assertEqual(reconcile_curated_vouchers_csv_with_db(conn, csv_path, taxa={"Packera obovata"}), 1)
                self.assertEqual(reconcile_curated_vouchers_csv_with_db(conn, csv_path, taxa={"Packera obovata"}), 0)
                self.assertEqual(reconcile_curated_vouchers_csv_with_db(conn, csv_path, taxa=set()), 0)
            finally:
                conn.close()

            text = csv_path.read_text(encoding="utf-8")
            self.assertTrue(text.startswith(original))
            self.assertEqual(text[len(original):], "O1,Packera obovata,Tier_1_Gold,True,2\n")

    def test_run_pipeline_sh_preflight_restores_sqlite_taxa_and_fails_on_missing_taxa(self):
        """Verify --taxon restores requested-taxon vouchers present in SQLite but missing from curated_vouchers.csv, and fails fast on taxa with zero vouchers."""
        from scripts.core.db import get_db_connection, upsert_vouchers_records

        def _run_preflight(run_dir: Path, root_data: Path, custom_taxa_csv: str):
            return _run_preflight_bash(
                run_dir,
                root_data,
                {
                    "FOCAL_TAXON": "",
                    "EXCLUDE_TAXON_CSV": "",
                    "INTRASPECIFIC_MODE": "false",
                    "GOLD_ANCHORS": "30",
                    "DETERMINER_TIER": "",
                    "REUSE_CANONICAL_DATA": "true",
                    "CUSTOM_TAXA_CSV": custom_taxa_csv,
                    "SEED_DATA_DIR": str(root_data),
                    "START_STEP": "3",
                },
                check=False,
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            root_data = tmp_path / "root_data"
            run_dir = tmp_path / "runs" / "test_sqlite_taxa_sync"
            tables_dir = run_dir / "data" / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)
            (root_data / "cropped_patches" / "capitula").mkdir(parents=True, exist_ok=True)

            pd.DataFrame({
                "catalogNumber": ["A1", "D1"],
                "scientificName": ["Packera anonyma", "Packera dubia"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold"],
                "is_primary_duplicate": [True, True],
                "duplicate_count": [1, 1],
            }).to_csv(tables_dir / "curated_vouchers.csv", index=False)

            conn = get_db_connection(db_path=root_data / "plantmatics.sqlite", init_schema=True)
            try:
                upsert_vouchers_records(conn, [
                    {"catalogNumber": "A1", "scientificName": "Packera anonyma", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "D1", "scientificName": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "O1", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold",
                     "species_standardized": "Packera obovata", "is_primary_duplicate": 1, "duplicate_count": 1},
                    {"catalogNumber": "P1", "scientificName": "Packera paupercula", "determiner_tier": "Tier_1_Gold",
                     "species_standardized": "Packera paupercula", "is_primary_duplicate": 1, "duplicate_count": 1},
                ])
            finally:
                conn.close()

            res = _run_preflight(run_dir, root_data, "Packera anonyma,Packera dubia,Packera obovata")
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertIn("[Preflight: SQLite Voucher Sync] Restored 1 vouchers", res.stdout)
            out_vouchers = pd.read_csv(tables_dir / "curated_vouchers.csv")
            self.assertEqual(set(out_vouchers["catalogNumber"]), {"A1", "D1", "O1"})

            res_missing = _run_preflight(run_dir, root_data, "Packera anonyma,Packera aurea")
            self.assertNotEqual(res_missing.returncode, 0)
            self.assertIn("Packera aurea", res_missing.stderr)
            self.assertIn("zero vouchers", res_missing.stderr)

    def test_run_pipeline_sh_exclude_misidentified_and_triage_source_filtering(self):
        """Verify run_pipeline.sh --exclude-misidentified (with default ALL and explicit --triage-source) filters out Misidentification_Reassign while retaining valid subtypes across curated_vouchers.csv, Step 02 manifests, cropped_patches, and SQLite."""
        sh_corpus = _load_pipeline_sh_corpus()
        self.assertIn("--exclude-misidentified)", sh_corpus)
        self.assertIn("--triage-source)", sh_corpus)
        self.assertIn("--triage-source=*)", sh_corpus)
        self.assertIn("Misidentification_Reassign", sh_corpus)

        all_cats = ["CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT", "CAT_MISID_1", "CAT_MISID_2"]
        expected_keep_all = {"CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT"}
        expected_keep_run_a = {"CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT", "CAT_MISID_2"}
        core_tbls = (
            "extracted_leaf_manifest.csv",
            "failed_qc_vouchers.csv",
            "voucher_reproductive_metrics.csv",
            "voucher_phenological_states.csv",
        )

        for mode, triage_src, expected_keep in [("all_default", "ALL", expected_keep_all), ("single_run", "run_A", expected_keep_run_a)]:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                root_data = tmp_path / "root_data"
                runs_root = tmp_path / "runs"
                run_dir = runs_root / f"test_{mode}"
                tables_dir = run_dir / "data" / "tables"
                tables_dir.mkdir(parents=True, exist_ok=True)

                run_a_tbl = runs_root / "run_A" / "data" / "tables"
                run_b_tbl = runs_root / "run_B" / "data" / "tables"
                run_a_tbl.mkdir(parents=True, exist_ok=True)
                run_b_tbl.mkdir(parents=True, exist_ok=True)

                pd.DataFrame({
                    "catalogNumber": ["CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT", "CAT_MISID_1"],
                    "taxonomic_treatment_subtype": [
                        "Confirmed_Species", "Hybrid_Introgressant", "Subspecies_Ecotype",
                        "Ecophenotypic_Plasticity", "Putative_Cryptic_Species", "Misidentification_Reassign",
                    ],
                }).to_csv(run_a_tbl / "triage_queue.csv", index=False)
                pd.DataFrame({
                    "catalogNumber": ["CAT_CONF", "CAT_MISID_2"],
                    "taxonomic_treatment_subtype": ["Confirmed_Species", "Misidentification_Reassign"],
                }).to_csv(run_b_tbl / "triage_queue.csv", index=False)

                pdir = root_data / "cropped_patches" / "capitula"
                pdir.mkdir(parents=True, exist_ok=True)
                for c in all_cats:
                    (pdir / f"{c}_cap1.jpg").write_bytes(b"jpg")

                pd.DataFrame({
                    "catalogNumber": all_cats,
                    "scientificName": ["Packera dubia"] * len(all_cats),
                    "determiner_tier": ["Tier_1_Gold"] * len(all_cats),
                    "authority_weight": [3.0] * len(all_cats),
                }).to_csv(tables_dir / "curated_vouchers.csv", index=False)
                for tbl_name in core_tbls:
                    pd.DataFrame({"catalogNumber": all_cats, "plant_individual_id": [0] * len(all_cats)}).to_csv(tables_dir / tbl_name, index=False)

                res = _run_preflight_bash(
                    run_dir,
                    root_data,
                    {
                        "FOCAL_TAXON": "",
                        "EXCLUDE_TAXON_CSV": "",
                        "INTRASPECIFIC_MODE": "false",
                        "GOLD_ANCHORS": "30",
                        "GOLD_ONLY": "false",
                        "DETERMINER_TIER": "",
                        "EXCLUDE_MISIDENTIFIED": "true",
                        "TRIAGE_SOURCE": triage_src,
                        "REUSE_CANONICAL_DATA": "true",
                    },
                    check=True,
                )
                self.assertIn("[Preflight: Exclude Misidentified]", res.stdout)
                self.assertEqual(set(pd.read_csv(tables_dir / "curated_vouchers.csv")["catalogNumber"]), expected_keep)
                for tbl_name in core_tbls:
                    self.assertEqual(set(pd.read_csv(tables_dir / tbl_name)["catalogNumber"]), expected_keep, f"Failed filter on {tbl_name} in {mode}")

                link_p = run_dir / "data" / "cropped_patches" / "capitula"
                self.assertTrue(link_p.is_symlink(), "Expected directory-level symlink for capitula")
                self.assertEqual(link_p.resolve(), pdir.resolve())

                db_file = root_data / "plantmatics.sqlite"
                self.assertTrue(db_file.is_file())
                with sqlite3.connect(db_file) as conn:
                    cohort_cats = {r[0] for r in conn.execute("SELECT catalogNumber FROM run_active_cohort WHERE run_id = ?;", (run_dir.name,)).fetchall()}
                self.assertEqual(cohort_cats, expected_keep)

        # Pure-SQLite --exclude-misidentified (with zero CSV triage reports on disk)
        with tempfile.TemporaryDirectory() as tmpdir_sql:
            tmp_sql = Path(tmpdir_sql)
            root_data_sql = tmp_sql / "root_data"
            root_data_sql.mkdir(parents=True, exist_ok=True)
            run_dir_sql = tmp_sql / "runs" / "test_sql_only"
            tables_dir_sql = run_dir_sql / "data" / "tables"
            tables_dir_sql.mkdir(parents=True, exist_ok=True)

            from scripts.core.db import ensure_pipeline_run, get_db_connection, upsert_rows_generic, upsert_vouchers_records

            conn_seed = get_db_connection(db_path=root_data_sql / "plantmatics.sqlite", init_schema=True)
            try:
                upsert_vouchers_records(conn_seed, [{"catalogNumber": c, "scientificName": "Packera dubia", "determiner_tier": "Tier_1_Gold"} for c in all_cats])
                ensure_pipeline_run(conn_seed, run_id="hist_sql_run")
                upsert_rows_generic(
                    conn_seed,
                    table_name="curation_triage",
                    rows=[
                        {"run_id": "hist_sql_run", "catalogNumber": "CAT_CONF", "plant_individual_id": 0, "taxonomic_treatment_subtype": "Confirmed_Species"},
                        {"run_id": "hist_sql_run", "catalogNumber": "CAT_MISID_1", "plant_individual_id": 0, "taxonomic_treatment_subtype": "Misidentification_Reassign"},
                    ],
                    pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
                )
            finally:
                conn_seed.close()

            pd.DataFrame({
                "catalogNumber": all_cats,
                "scientificName": ["Packera dubia"] * len(all_cats),
                "determiner_tier": ["Tier_1_Gold"] * len(all_cats),
                "authority_weight": [3.0] * len(all_cats),
            }).to_csv(tables_dir_sql / "curated_vouchers.csv", index=False)

            res_sql = _run_preflight_bash(
                run_dir_sql,
                root_data_sql,
                {
                    "FOCAL_TAXON": "",
                    "EXCLUDE_TAXON_CSV": "",
                    "INTRASPECIFIC_MODE": "false",
                    "GOLD_ANCHORS": "30",
                    "GOLD_ONLY": "false",
                    "DETERMINER_TIER": "",
                    "EXCLUDE_MISIDENTIFIED": "true",
                    "TRIAGE_SOURCE": "hist_sql_run",
                    "REUSE_CANONICAL_DATA": "false",
                },
                check=True,
            )
            self.assertIn("SQLite runs=1", res_sql.stdout)
            self.assertEqual(set(pd.read_csv(tables_dir_sql / "curated_vouchers.csv")["catalogNumber"]), expected_keep_run_a)


if __name__ == "__main__":
    unittest.main()
