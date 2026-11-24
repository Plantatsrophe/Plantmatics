#!/usr/bin/env python3
"""
===============================================================================
Script: rescore_step02_qc_and_zones.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Fast CPU post-Step-02 re-scoring utility that refreshes Tier A leaf QC flags
    (`detector_class_id`, `is_crop_clipped`, `hole_area_frac`, `enclosed_void_frac`,
    `bilateral_lobe_symmetry`, `efa_qc_pass`, `qc_flag_reason`), read-only Packera
    margin/base/terminal-lobe traits, void-aware `solidity` / `circularity`,
    calibrated CIELAB photometry plus rosette abaxial-leaf flags, and anatomical
    `basal` vs. `cauline` leaf zones across `extracted_leaves` in SQLite and
    `data/tables/extracted_leaf_manifest.csv` without mutating masks or contours
    or re-running GPU inference.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig, resolve_max_workers
from scripts.core.db import get_db_connection
from scripts.core.io_utils import atomic_write_text
from scripts.vision.leaf_geometry_utils import compute_lm3_shape_metrics
from scripts.vision.lm3_compat import rotate_image
from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES
from scripts.vision.lm3_inference_engines import (
    LM3ScaleCalibrator,
    evaluate_leaf_damage_and_fold_qc,
    solve_ruler_tick_lattice_px_per_mm,
)
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical
from scripts.vision.step02_leaf_traits_utils import (
    evaluate_rosette_abaxial_leaf_flags,
    extract_packera_leaf_margin_and_base_traits,
    measure_enclosed_interlobe_void_frac,
)
from scripts.vision.step02_photometry_corymb_utils import (
    compute_calibrated_blade_photometry,
    estimate_sheet_paper_cie_l_ref,
)
from scripts.vision.step02_rescore_photometry import (
    _photo_sheet_worker,
    recompute_leaf_photometry_and_abaxial_flags,
)
from scripts.vision.step02_rescore_sqlite_sync import (
    persist_rescored_leaves_to_sqlite,
    rescore_anatomical_leaf_zones,
    rescore_tier_a_leaves_impl,
    sync_rescored_csv_snapshots,
)
from scripts.vision.step02_rescore_workers import (
    _align_rotated_patch_to_mask,
    _calibrate_voucher_ruler_worker,
    _eval_single_leaf_qc_worker,
    _resolve_crop_file,
    _resolve_mask_file,
    _resolve_voucher_image_file,
    run_leaf_qc_evaluation,
    run_ruler_calibration_backfill,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("RescoreStep02")


def rescore_tier_a_leaves(
    db_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    manifest_csv: Optional[Path] = None,
    segment_partial_leaves: Optional[bool] = None,
    recompute_photometry: bool = False,
    workers: Optional[int] = None,
) -> Dict[str, int]:
    """Re-evaluates read-only leaf QC gates, Packera traits, ruler scales, photometry, and basal/cauline zones."""
    return rescore_tier_a_leaves_impl(
        db_path=db_path,
        data_dir=data_dir,
        manifest_csv=manifest_csv,
        segment_partial_leaves=segment_partial_leaves,
        recompute_photometry=recompute_photometry,
        workers=workers,
    )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Re-score Step 02 leaf QC gates and anatomical basal/cauline zones in SQLite + CSV without GPU inference."
    )
    parser.add_argument("--db", type=Path, default=None, help="Optional SQLite database path (default: data/plantmatics.sqlite).")
    parser.add_argument("--data-dir", type=Path, default=PROJECT_ROOT / "data", help="Tier A data directory containing masks/, contours/, and tables/.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional extracted_leaf_manifest.csv path override.")
    parser.add_argument("--allow-partial-leaves", action="store_true", default=False, help="Allow detector_class_id == 1 (Leaf_PARTIAL) to pass QC.")
    parser.add_argument("--recompute-photometry", action="store_true", default=False, help="Recompute blade_cie_l_norm/sd where crops or sheet-matched masks exist.")
    parser.add_argument("--workers", "-j", "--threads", dest="workers", type=int, default=None, help="Worker thread/process count cap (default: auto up to 64).")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    rescore_tier_a_leaves(
        db_path=args.db,
        data_dir=args.data_dir,
        manifest_csv=args.manifest,
        segment_partial_leaves=True if args.allow_partial_leaves else None,
        recompute_photometry=bool(args.recompute_photometry),
        workers=args.workers,
    )


__all__ = [
    "PROJECT_ROOT",
    "logger",
    "ProcessPoolExecutor",
    "ThreadPoolExecutor",
    "PipelineConfig",
    "resolve_max_workers",
    "get_db_connection",
    "atomic_write_text",
    "compute_lm3_shape_metrics",
    "rotate_image",
    "_CALIBRATED_SCALE_SOURCES",
    "LM3ScaleCalibrator",
    "evaluate_leaf_damage_and_fold_qc",
    "solve_ruler_tick_lattice_px_per_mm",
    "classify_leaf_zone_anatomical",
    "evaluate_rosette_abaxial_leaf_flags",
    "extract_packera_leaf_margin_and_base_traits",
    "measure_enclosed_interlobe_void_frac",
    "compute_calibrated_blade_photometry",
    "estimate_sheet_paper_cie_l_ref",
    "_resolve_mask_file",
    "_resolve_voucher_image_file",
    "_resolve_crop_file",
    "_align_rotated_patch_to_mask",
    "_eval_single_leaf_qc_worker",
    "_calibrate_voucher_ruler_worker",
    "_photo_sheet_worker",
    "run_ruler_calibration_backfill",
    "run_leaf_qc_evaluation",
    "recompute_leaf_photometry_and_abaxial_flags",
    "rescore_anatomical_leaf_zones",
    "persist_rescored_leaves_to_sqlite",
    "sync_rescored_csv_snapshots",
    "rescore_tier_a_leaves",
    "parse_args",
    "main",
]


if __name__ == "__main__":
    main()
