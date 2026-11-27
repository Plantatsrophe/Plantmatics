#!/usr/bin/env python3
"""
Unit tests for compute flags (-j/--threads, --gpu-batch-size), help hygiene,
--purge-segments / --purge-runs-only execution, focal complex dry-run CLI,
and Python/Bash CLI flag & environment export parity.
"""

import argparse
import importlib
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main

_HEAVY_MOD_NAMES = (
    "pandas", "numpy", "cv2", "PIL", "PIL.Image", "scipy",
    "scipy.ndimage", "scipy.spatial", "scipy.spatial.distance",
    "sklearn", "sklearn.cluster", "sklearn.linear_model",
    "sklearn.calibration", "sklearn.model_selection", "sklearn.metrics",
    "sklearn.preprocessing", "matplotlib", "matplotlib.pyplot",
    "matplotlib.patches", "matplotlib.gridspec", "pygbif", "aiohttp", "tqdm",
    "torch", "torch.nn", "torch.nn.functional", "torch.utils",
    "torch.utils.data", "torchvision", "torchvision.transforms",
)


def _fake_heavy_modules() -> dict:
    return {k: MagicMock() for k in _HEAVY_MOD_NAMES if k not in sys.modules}


class TestCLIPurgeAndComputeFlags(unittest.TestCase):
    """Test suite for compute flags, purge commands, focal CLI, and Python/Bash parity."""

    def test_cli_threads_and_gpu_batch_size_parsing_and_env_export(self):
        """Verify -j/--threads and --gpu-batch-size parsing, config resolution, and subprocess forwarding."""
        from scripts.core.config import resolve_cpu_threads, resolve_gpu_batch_size

        parser = main.build_parser()
        args_seg = parser.parse_args(["segment", "-j", "12", "--gpu-batch-size", "24"])
        self.assertEqual(args_seg.threads, 12)
        self.assertEqual(args_seg.gpu_batch_size, 24)

        args_alias = parser.parse_args(["synthesis", "--threads", "6", "--batch-size", "16", "--ssurgo-batch-size", "50"])
        self.assertEqual(args_alias.threads, 6)
        self.assertEqual(args_alias.gpu_batch_size, 16)
        self.assertEqual(args_alias.ssurgo_batch_size, 50)

        with patch.dict(os.environ, {}, clear=False):
            for env_k in ("PLANTMATICS_MAX_WORKERS", "PLANTMATICS_GPU_BATCH_SIZE", "PACKERA_MAX_WORKERS", "PACKERA_GPU_BATCH_SIZE"):
                os.environ.pop(env_k, None)
            main._apply_explicit_db_arg(args_seg)
            self.assertEqual(os.environ.get("PACKERA_MAX_WORKERS"), "12")
            self.assertEqual(os.environ.get("PACKERA_GPU_BATCH_SIZE"), "24")
            self.assertEqual(resolve_cpu_threads(None, default=4), 12)
            self.assertEqual(resolve_gpu_batch_size(None, default=16), 24)
            self.assertEqual(resolve_cpu_threads(2, default=4), 2)
            self.assertEqual(resolve_gpu_batch_size(8, default=16), 8)

        with patch.dict(sys.modules, _fake_heavy_modules()):
            mod_01 = importlib.import_module("scripts.data_prep.01_voucher_harvester")
            p01 = mod_01.parse_args(["--taxon", "Packera dubia", "--geographic-scope", "north_america", "--max-records", "250"])
            self.assertEqual(p01.taxa, ["Packera dubia"])
            self.assertEqual(p01.geographic_scope, "north_america")
            self.assertEqual(p01.max_records, 250)
            p01_legacy = mod_01.parse_args(["--taxa", "Packera dubia", "--include-western", "--max-records-per-taxon", "120", "--overwrite"])
            self.assertEqual(p01_legacy.taxa, ["Packera dubia"])
            self.assertTrue(p01_legacy.include_western)
            self.assertEqual(p01_legacy.max_records, 120)
            self.assertTrue(p01_legacy.force)

            mod_02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
            with patch.object(sys, "argv", ["02_segment_and_extract.py", "-j", "10", "--gpu-batch-size", "20"]):
                p02 = mod_02.parse_args()
                self.assertEqual(p02.threads, 10)
                self.assertEqual(p02.gpu_batch_size, 20)

    def test_subcommand_help_hides_suppressed_internal_and_legacy_flags(self):
        """Verify --help output exposes consolidated flags and hides internal/retired aliases."""
        parser = main.build_parser()
        subparsers_action = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
        synth_help = subparsers_action.choices["synthesis"].format_help()
        for flag in ("--threads", "--gpu-batch-size", "--ssurgo-batch-size"):
            self.assertIn(flag, synth_help)
        for hidden_flag in ("--ssurgo-table", "--contingency-table", "--synthesis-report"):
            self.assertNotIn(hidden_flag, synth_help)

        run_all_help = subparsers_action.choices["run-all"].format_help()
        for flag in ("--taxon", "--threads", "--gpu-batch-size"):
            self.assertIn(flag, run_all_help)
        for hidden_flag in ("--overwrite", "--include-western", "--allometry-report"):
            self.assertNotIn(hidden_flag, run_all_help)

    def test_run_pipeline_sh_help_and_compute_flags(self):
        """Verify ./run_pipeline.sh --help exits 0 and documents -j/--threads and --gpu-batch-size."""
        res = subprocess.run(
            ["bash", str(PROJECT_ROOT / "run_pipeline.sh"), "--help"],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(res.returncode, 0)
        for active_flag in (
            "-j, --threads", "--gpu-batch-size", "--taxon", "--skip-harvest",
            "--resegment", "--purge-segments", "--purge-old-data", "--gold-only", "--auto-neighbors",
        ):
            self.assertIn(active_flag, res.stdout)

    def test_purge_segments_cli_and_purge_segments_and_runs_execution(self):
        """Verify --purge-segments / --purge-old-data CLI parsing and purge_segments_and_runs execution."""
        from scripts.core.db import get_db_connection, purge_segments_and_runs

        parser = main.build_parser()
        args_top = parser.parse_args(["--purge-segments"])
        self.assertTrue(args_top.purge_segments)
        self.assertIsNone(args_top.subcommand)

        args_top_alias = parser.parse_args(["--purge-old-data", "--no-purge-runs"])
        self.assertTrue(args_top_alias.purge_segments)
        self.assertFalse(args_top_alias.purge_runs)

        self.assertEqual(parser.parse_args(["purge-segments"]).subcommand, "purge-segments")
        args_sub_alias = parser.parse_args(["purge-old-data", "--no-purge-runs"])
        self.assertEqual(args_sub_alias.subcommand, "purge-old-data")
        self.assertFalse(args_sub_alias.purge_runs)

        self.assertTrue(parser.parse_args(["segment", "--purge-segments"]).purge_segments)
        self.assertTrue(parser.parse_args(["run-all", "--purge-old-data"]).purge_segments)

        with patch.dict(sys.modules, _fake_heavy_modules()):
            mod_02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
            with patch.object(sys, "argv", ["02_segment_and_extract.py", "--purge-segments"]):
                self.assertTrue(mod_02.parse_args().purge_segments)

        with tempfile.TemporaryDirectory() as tmpdir:
            ws = Path(tmpdir)
            for rel, payload in (
                ("data/raw_vouchers/NCU0001.webp", b"WEBP"),
                ("data/masks/NCU0001_leaf_1.png", b"PNG"),
                ("data/ect/NCU0001_leaf_1_ect.npy", b"NPY"),
                ("data/cropped_patches/rosette/NCU0001_leaf_1.webp", b"WEBP"),
            ):
                p = ws / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(payload)
            for rel, txt in (
                ("data/tables/curated_vouchers.csv", "catalogNumber\nNCU0001\n"),
                ("data/tables/extracted_leaf_traits.csv", "catalogNumber\nNCU0001\n"),
                ("data/contours/NCU0001_leaf_1.csv", "x,y\n0,0\n"),
                ("runs/20260101_000000/run.log", "log"),
            ):
                p = ws / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(txt, encoding="utf-8")

            db_file = ws / "data" / "plantmatics.sqlite"
            conn = get_db_connection(db_file, init_schema=True)
            conn.execute("INSERT INTO vouchers (catalogNumber, scientificName, qc_status) VALUES ('NCU0001', 'Packera aurea', 'PASS')")
            conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU0001', 1)")
            conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, contour_xy_json) VALUES ('NCU0001', 1, 'L1', '[[0,0]]')")
            conn.commit()
            conn.close()

            summary = purge_segments_and_runs(workspace_root=ws, db_path=db_file, purge_runs=True, vacuum=True)
            self.assertGreaterEqual(summary["removed_files"], 5)
            self.assertEqual(summary["purged_db_rows"], 2)

            self.assertTrue((ws / "data" / "raw_vouchers" / "NCU0001.webp").exists())
            self.assertTrue((ws / "data" / "tables" / "curated_vouchers.csv").exists())
            self.assertFalse((ws / "data" / "tables" / "extracted_leaf_traits.csv").exists())
            self.assertEqual(list((ws / "data" / "contours").glob("*.csv")), [])
            self.assertEqual(list((ws / "data" / "ect").glob("*.npy")), [])
            self.assertTrue((ws / "data" / "masks").is_dir())
            self.assertEqual(list((ws / "data" / "masks").glob("*")), [])
            self.assertEqual(list((ws / "runs").iterdir()), [])

            conn_check = sqlite3.connect(str(db_file))
            try:
                self.assertEqual(conn_check.execute("SELECT catalogNumber, qc_status FROM vouchers").fetchone(), ("NCU0001", None))
                self.assertEqual(conn_check.execute("SELECT COUNT(*) FROM plant_individuals").fetchone()[0], 0)
                self.assertEqual(conn_check.execute("SELECT COUNT(*) FROM extracted_leaves").fetchone()[0], 0)
            finally:
                conn_check.close()

    def test_purge_runs_only_preserves_step01_and_step02_segmentation(self) -> None:
        """Verify --purge-runs-only / purge_runs_only wipes Step 03-07 run outputs & Tier B tables while preserving Step 01 & Step 02 segmentation."""
        from scripts.core.db import get_db_connection, purge_runs_only
        from scripts.pipeline.resolve_focal_species_complex import build_four_stage_commands

        parser = main.build_parser()
        args_top = parser.parse_args(["--purge-runs-only"])
        self.assertTrue(args_top.purge_runs_only)
        self.assertFalse(args_top.purge_segments)
        self.assertEqual(parser.parse_args(["purge-runs-only"]).subcommand, "purge-runs-only")

        stages_skip = build_four_stage_commands(
            run_script=PROJECT_ROOT / "run_pipeline.sh",
            focal_taxon="Packera dubia",
            prefix="dubia",
            skip_harvest=True,
        )
        self.assertIn("--skip-harvest", stages_skip[0][2])

        with tempfile.TemporaryDirectory() as tmpdir:
            ws = Path(tmpdir)
            for rel, payload in (
                ("data/raw_vouchers/NCU0001.webp", b"WEBP"),
                ("data/masks/NCU0001_leaf_1.png", b"PNG"),
                ("data/ect/NCU0001_leaf_1_ect.npy", b"NPY"),
                ("data/cropped_patches/capitula/NCU0001_cap_1.webp", b"WEBP"),
                ("runs/dubia_04_gold_thin_audit/outputs/figures/fig05a.png", b"PNG"),
            ):
                p = ws / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(payload)
            (ws / "data" / "contours").mkdir(parents=True, exist_ok=True)
            (ws / "data" / "contours" / "NCU0001_leaf_1.csv").write_text("x,y\n0,0\n", encoding="utf-8")

            tables = ws / "data" / "tables"
            tables.mkdir(parents=True, exist_ok=True)
            s01_02_tbls = (
                "curated_vouchers.csv", "extracted_leaf_manifest.csv", "failed_qc_vouchers.csv",
                "voucher_reproductive_metrics.csv", "voucher_phenological_states.csv",
            )
            s03_07_tbls = ("leaf_efa_harmonics.csv", "morphometric_flags.csv", "multimodal_conflict_flags.csv", "triage_queue.csv")
            for tbl in (*s01_02_tbls, *s03_07_tbls):
                (tables / tbl).write_text("catalogNumber\nNCU0001\n", encoding="utf-8")

            env_dir = ws / "data" / "environmental"
            env_dir.mkdir(parents=True, exist_ok=True)
            (env_dir / "fine_grain_edaphic_cache.csv").write_text(
                "catalogNumber,lat_lon_key,fine_soil_ph,min_bedrock_depth_buffer_cm\nNCU0001,35.90000_-79.05000,5.1,162.4\n",
                encoding="utf-8",
            )
            (ws / "outputs" / "reports").mkdir(parents=True, exist_ok=True)
            (ws / "outputs" / "reports" / "old_summary.md").write_text("old", encoding="utf-8")

            db_file = ws / "data" / "plantmatics.sqlite"
            conn = get_db_connection(db_file, init_schema=True)
            conn.execute("INSERT INTO vouchers (catalogNumber, scientificName, qc_status, pixels_per_mm, is_calibrated_mm) VALUES ('NCU0001', 'Packera dubia', 'PASS', 24.5, 1);")
            conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU0001', 1);")
            conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, leaf_zone, contour_xy_json) VALUES ('NCU0001', 1, 'L1', 'basal', '[[0,0]]');")
            conn.execute("INSERT INTO environmental_cache (lat_lon_key, catalogNumber, fine_soil_ph, min_bedrock_depth_buffer_cm) VALUES ('35.90000_-79.05000', 'NCU0001', 5.1, 162.4);")
            conn.execute("INSERT INTO pipeline_runs (run_id, run_dir) VALUES ('dubia_04_gold_thin_audit', 'runs/dubia_04_gold_thin_audit');")
            conn.execute("INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id) VALUES ('dubia_04_gold_thin_audit', 'NCU0001', 1);")
            conn.execute("INSERT INTO otu_morphometrics (run_id, catalogNumber, plant_individual_id, gmm_cluster) VALUES ('dubia_04_gold_thin_audit', 'NCU0001', 1, '1');")
            conn.execute("INSERT INTO spatial_rf_predictions (run_id, catalogNumber, plant_individual_id, rf_predicted_taxon, ecotype_regime) VALUES ('dubia_04_gold_thin_audit', 'NCU0001', 1, 'Packera dubia', 'Ecotype_1_lithic_outcrop_index_High');")
            conn.execute("INSERT INTO curation_triage (run_id, catalogNumber, plant_individual_id, taxonomic_treatment_subtype) VALUES ('dubia_04_gold_thin_audit', 'NCU0001', 1, 'Confirmed_Species');")
            conn.commit()
            conn.close()

            summary = purge_runs_only(workspace_root=ws, db_path=db_file, refresh_env_buffers=True, vacuum=True)
            self.assertGreaterEqual(summary["removed_files"], 6)
            self.assertEqual(summary["purged_db_rows"], 5)

            for keep_rel in (
                "data/raw_vouchers/NCU0001.webp", "data/masks/NCU0001_leaf_1.png",
                "data/contours/NCU0001_leaf_1.csv", "data/ect/NCU0001_leaf_1_ect.npy",
                "data/cropped_patches/capitula/NCU0001_cap_1.webp",
            ):
                self.assertTrue((ws / keep_rel).is_file())
            for tbl in s01_02_tbls:
                self.assertTrue((tables / tbl).is_file(), f"Step 01/02 table {tbl} must be preserved!")
            for tbl in s03_07_tbls:
                self.assertFalse((tables / tbl).exists(), f"Step 03-07 table {tbl} must be purged!")
            self.assertEqual(list((ws / "runs").iterdir()), [])
            self.assertEqual(list((ws / "outputs" / "reports").iterdir()), [])

            conn_chk = sqlite3.connect(str(db_file))
            try:
                self.assertEqual(conn_chk.execute("SELECT catalogNumber, qc_status, pixels_per_mm FROM vouchers;").fetchone(), ("NCU0001", "PASS", 24.5))
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM plant_individuals;").fetchone()[0], 1)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0], 1)
                for tier_b_tbl in ("pipeline_runs", "run_active_cohort", "otu_morphometrics", "spatial_rf_predictions", "curation_triage"):
                    self.assertEqual(conn_chk.execute(f"SELECT COUNT(*) FROM {tier_b_tbl};").fetchone()[0], 0)
                self.assertEqual(
                    conn_chk.execute("SELECT fine_soil_ph, min_bedrock_depth_buffer_cm FROM environmental_cache WHERE catalogNumber='NCU0001';").fetchone(),
                    (5.1, None),
                )
            finally:
                conn_chk.close()

            res_sh = subprocess.run(
                ["bash", str(PROJECT_ROOT / "run_pipeline.sh"), "--purge-runs-only", "--db", str(db_file)],
                cwd=ws,
                env={**os.environ, "PLANTMATICS_WORKSPACE_ROOT": str(ws), "PACKERA_WORKSPACE_ROOT": str(ws)},
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(res_sh.returncode, 0, f"run_pipeline.sh --purge-runs-only failed:\n{res_sh.stderr}\n{res_sh.stdout}")
            self.assertIn("Step 01 vouchers & Step 02 segments preserved", res_sh.stdout)

    def test_resolve_focal_species_complex_dry_run_and_cli_flags(self) -> None:
        """Verify resolve_focal_species_complex.py and resolve_packera_dubia_complex.sh construct all 4 ordered stages."""
        from scripts.pipeline.resolve_focal_species_complex import build_four_stage_commands

        stages = build_four_stage_commands(
            run_script=PROJECT_ROOT / "run_pipeline.sh",
            focal_taxon="Packera dubia",
            prefix="dubia",
            geographic_scope="eastern_na",
            gold_anchors=30,
            spatial_thin_km=15.0,
            threads=12,
            gpu_batch_size=16,
        )
        self.assertEqual(len(stages), 4)
        self.assertEqual(stages[0][1], "dubia_01_eastern_discovery")
        self.assertEqual(stages[1][1], "dubia_02_complex_triage")
        self.assertIn("--discover-complexes", stages[1][2])
        self.assertEqual(stages[2][1], "dubia_03_intraspecific_clean")
        self.assertIn("--intraspecific-mode", stages[2][2])
        self.assertIn("--exclude-misidentified", stages[2][2])
        self.assertEqual(stages[3][1], "dubia_04_gold_thin_audit")
        self.assertIn("--spatial-thin-km", stages[3][2])

        proc = subprocess.run(
            ["bash", str(PROJECT_ROOT / "resolve_packera_dubia_complex.sh"), "--dry-run", "--threads", "8"],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertRegex(proc.stdout, r"dubia_\d{8}_[0-9a-f]{6}/01_eastern_discovery")
        for sub_st in (
            "01_eastern_discovery",
            "02_complex_triage",
            "03_intraspecific_clean",
            "04_spatial_thin_field_audit",
            "05_monographer_gold_audit",
            "06_lineage_diagnosability",
        ):
            self.assertIn(sub_st, proc.stdout)

    def test_python_and_bash_cli_flag_parity(self) -> None:
        """Verify Python main.build_parser() and Bash run_pipeline.sh + pipe_cli_args.sh maintain flag and dual-env export parity."""
        parser = main.build_parser()
        for subcmd in ("harvest", "synthesis", "run-all"):
            ns = parser.parse_args([
                subcmd,
                "--gold-only",
                "--determiner-tier", "Tier_1_Gold",
                "--exclude-misidentified",
                "--triage-source", "run_prev",
                "--auto-neighbors",
            ])
            self.assertEqual(ns.subcommand, subcmd)
            self.assertTrue(getattr(ns, "gold_only", False))
            self.assertEqual(str(getattr(ns, "determiner_tier", "")).strip(), "Tier_1_Gold")
            self.assertTrue(getattr(ns, "exclude_misidentified", False))
            triage_val = getattr(ns, "triage_source", None)
            if isinstance(triage_val, (list, tuple)):
                self.assertIn("run_prev", triage_val)
            else:
                self.assertEqual(triage_val, "run_prev")
            self.assertTrue(getattr(ns, "auto_neighbors", False))

        args_env = parser.parse_args(["segment", "-j", "6", "--gpu-batch-size", "12", "--db", ":memory:"])
        with patch.dict(os.environ, {}, clear=False):
            for k in (
                "PLANTMATICS_MAX_WORKERS", "PACKERA_MAX_WORKERS",
                "PLANTMATICS_GPU_BATCH_SIZE", "PACKERA_GPU_BATCH_SIZE",
                "PLANTMATICS_DB_PATH", "PACKERA_DB_PATH",
            ):
                os.environ.pop(k, None)
            main._apply_explicit_db_arg(args_env)
            self.assertEqual(os.environ.get("PLANTMATICS_MAX_WORKERS"), "6")
            self.assertEqual(os.environ.get("PACKERA_MAX_WORKERS"), "6")
            self.assertEqual(os.environ.get("PLANTMATICS_GPU_BATCH_SIZE"), "12")
            self.assertEqual(os.environ.get("PACKERA_GPU_BATCH_SIZE"), "12")
            self.assertEqual(os.environ.get("PLANTMATICS_DB_PATH"), ":memory:")
            self.assertEqual(os.environ.get("PACKERA_DB_PATH"), ":memory:")

        sh_dir = PROJECT_ROOT / "scripts" / "pipeline" / "sh"
        cli_args_path = sh_dir / "pipe_cli_args.sh"
        sh_files = [PROJECT_ROOT / "run_pipeline.sh"]
        if sh_dir.is_dir():
            sh_files.extend(sorted(sh_dir.glob("*.sh")))
        sh_corpus = "\n".join(p.read_text(encoding="utf-8") for p in sh_files if p.is_file())
        cli_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        if cli_args_path.is_file():
            cli_src += "\n" + cli_args_path.read_text(encoding="utf-8")

        for flag in (
            "--purge-runs-only",
            "-j",
            "--threads",
            "--skip-nwi-network",
            "--harvest-full-continental-range",
            "--no-harvest-full-continental-range",
        ):
            self.assertRegex(
                cli_src,
                rf"(?m)^\s*[^#\n]*(?<![A-Za-z0-9_-]){re.escape(flag)}(?=[=)|*\s])[^#\n]*\)",
                f"Missing Bash CLI flag handler for {flag}",
            )

        for suffix in ("MAX_WORKERS", "GPU_BATCH_SIZE", "RUN_DIR", "DB_PATH"):
            for prefix in ("PLANTMATICS", "PACKERA"):
                self.assertRegex(
                    sh_corpus,
                    rf"\bexport\s+{prefix}_{suffix}\b",
                    f"Missing Bash export for {prefix}_{suffix}",
                )


if __name__ == "__main__":
    unittest.main()
