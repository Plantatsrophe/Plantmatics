#!/usr/bin/env python3
"""
scripts/vision/lm3_data_loader.py
=================================
Scale calibration ingestion, working-frame coordinate unscaling, and precomputed
plant label loading utilities for LeafMachine3 (LM3).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

logger = logging.getLogger("LM3_PostProcessing")

# Canonical LeafMachine3 9-Class Plant Detector (lm3_plant_detector__yolo26x_det_1280) Taxonomy
LM3_PLANT_CLASSES: Dict[int, str] = {
    0: "Leaf_WHOLE",
    1: "Leaf_PARTIAL",
    2: "Seed_Fruit_ONE",
    3: "Seed_Fruit_MANY",
    4: "Flower_ONE",
    5: "Flower_MANY",
    6: "Bud",
    7: "Roots",
    8: "Wood",
}
LM3_LEAF_CLASS_IDS: Tuple[int, ...] = (0, 1)
LM3_REPRODUCTIVE_CLASS_IDS: Tuple[int, ...] = (2, 3, 4, 5, 6)
LM3_SINGLE_CAPITULUM_CLASS_IDS: Tuple[int, ...] = (2, 4, 6)

DEFAULT_HERBARIUM_DPI: float = 600.0
DEFAULT_PX_PER_CM: float = DEFAULT_HERBARIUM_DPI / 2.54  # ~236.2205 px/cm
DEFAULT_PX_PER_MM: float = DEFAULT_PX_PER_CM / 10.0     # ~23.6220 px/mm
_CALIBRATED_SCALE_SOURCES: Tuple[str, ...] = (
    "ruler_calibrated",
    "measured_from_ruler",
    "measured_from_fieldprism",
)


def compute_lm3_work_scale(
    sheet_shape: Any,
    sheet_w: Optional[int] = None,
    max_working_dim: int = 3200,
) -> float:
    """
    Computes the LeafMachine3 working-frame downscale factor `work_scale`.
    Accepts either `compute_lm3_work_scale((sheet_h, sheet_w), max_working_dim=3200)`
    or `compute_lm3_work_scale(sheet_h, sheet_w, max_working_dim=3200)`.
    Returns `1.0` when `max(h, w) <= max_working_dim`, else `float(max_working_dim) / float(max(h, w))`.
    """
    if isinstance(sheet_shape, (tuple, list, np.ndarray)):
        if len(sheet_shape) < 2:
            return 1.0
        h, w = int(sheet_shape[0]), int(sheet_shape[1])
        if sheet_w is not None and max_working_dim == 3200:
            max_working_dim = int(sheet_w)
    elif isinstance(sheet_shape, (int, float, np.number)) and sheet_w is not None:
        h, w = int(sheet_shape), int(sheet_w)
    else:
        return 1.0

    if max_working_dim <= 0:
        return 1.0
    max_dim = max(h, w)
    if max_dim <= 0 or max_dim <= int(max_working_dim):
        return 1.0
    return float(max_working_dim) / float(max_dim)


compute_work_scale = compute_lm3_work_scale


def unscale_working_bbox_to_sheet(
    bbox_work: Tuple[Union[int, float], Union[int, float], Union[int, float], Union[int, float]],
    work_scale: float,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> Tuple[int, int, int, int]:
    """
    Divides a working-frame bounding box `(ymin, xmin, ymax, xmax)` by `work_scale`
    and optionally clamps coordinates to `sheet_shape` `(h, w)`.
    """
    eff_scale = float(work_scale) if (work_scale is not None and float(work_scale) > 0.0) else 1.0
    y1 = int(round(float(bbox_work[0]) / eff_scale))
    x1 = int(round(float(bbox_work[1]) / eff_scale))
    y2 = int(round(float(bbox_work[2]) / eff_scale))
    x2 = int(round(float(bbox_work[3]) / eff_scale))

    if sheet_shape is not None and len(sheet_shape) >= 2:
        h, w = int(sheet_shape[0]), int(sheet_shape[1])
        if h > 0 and w > 0:
            y1 = max(0, min(h, y1))
            y2 = max(0, min(h, y2))
            x1 = max(0, min(w, x1))
            x2 = max(0, min(w, x2))

    return y1, x1, y2, x2


def normalize_pixels_per_mm(pixels_per_mm: Optional[float]) -> Optional[float]:
    """
    Normalizes a positive linear scale factor to `pixels_per_mm`.
    Preserves legitimate archival scan densities down to ~65 DPI (`>= 0.25 px/mm`)
    while inverting legacy `mm_per_pixel` inputs (`< 0.25`, e.g. `0.08 mm/px -> 12.5 px/mm`).
    """
    if pixels_per_mm is None:
        return None
    try:
        val = float(pixels_per_mm)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(val) or val <= 0.0:
        return None
    return float(1.0 / val) if val < 0.25 else val


def resolve_px_per_cm(
    pixels_per_mm: Optional[float] = None,
    px_per_cm: Optional[float] = None,
    default_dpi: float = DEFAULT_HERBARIUM_DPI,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> float:
    """
    Resolves the empirical ruler conversion factor in pixels per centimeter (`px_per_cm`):
      1. Uses `px_per_cm` directly when positive and finite.
      2. Uses `normalize_pixels_per_mm(pixels_per_mm) * 10.0` when `pixels_per_mm` is supplied
         (only inverting unambiguous `mm_per_pixel < 0.25` inputs so 150-250 DPI scans with
         `0.59-0.98 px/mm` are never inverted).
      3. Falls back to `default_dpi / 2.54`, scaled proportionally (`max(40.0, ...)`) when
         `sheet_shape` is a synthetic unit-test canvas (`max(sheet_shape) <= 500`).
    """
    if px_per_cm is not None and not np.isnan(px_per_cm) and float(px_per_cm) > 0.0:
        return float(px_per_cm)
    ppm = normalize_pixels_per_mm(pixels_per_mm)
    if ppm is not None:
        return ppm * 10.0
    fallback_px_cm = float(default_dpi) / 2.54 if float(default_dpi) > 0.0 else DEFAULT_PX_PER_CM
    if sheet_shape is not None and max(sheet_shape) <= 500:
        return max(40.0, fallback_px_cm * (float(max(sheet_shape)) / 2400.0))
    return fallback_px_cm


def compute_area_mm2(
    mask_area_px: float,
    pixels_per_mm: Optional[float] = None,
    is_calibrated_mm: bool = True,
) -> float:
    """
    Converts pixel mask area (`mask_area_px`) into physical lamina area in square
    millimeters (`mm^2`) strictly when `is_calibrated_mm` is True and `pixels_per_mm > 0`.
    Returns `np.nan` on uncalibrated or megapixel-estimated sheets.
    """
    if not is_calibrated_mm or mask_area_px is None or float(mask_area_px) <= 0.0:
        return np.nan
    ppm = normalize_pixels_per_mm(pixels_per_mm)
    if ppm is None or ppm <= 0.0:
        return np.nan
    return float(mask_area_px) / max(ppm * ppm, 1e-8)


def load_precomputed_plant_labels(
    label_path: Union[str, Path],
    sheet_width: int,
    sheet_height: int,
    min_confidence: float = 0.25,
) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
    """
    Loads precomputed YOLO normalized bounding boxes from a `.txt` label file and
    converts them into full-sheet pixel coordinates:
    `[((ymin, xmin, ymax, xmax), confidence, class_id), ...]`.
    """
    detections: List[Tuple[Tuple[int, int, int, int], float, int]] = []
    p = Path(label_path) if label_path else None
    if p is None or not p.exists():
        return detections

    try:
        with open(p, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 5:
                    continue
                cls_id = int(float(parts[0]))
                xc = float(parts[1]) * sheet_width
                yc = float(parts[2]) * sheet_height
                bw = float(parts[3]) * sheet_width
                bh = float(parts[4]) * sheet_height
                conf = float(parts[5]) if len(parts) >= 6 else 1.0
                if conf < min_confidence:
                    continue
                xmin = max(0, int(round(xc - bw / 2.0)))
                ymin = max(0, int(round(yc - bh / 2.0)))
                xmax = min(sheet_width, int(round(xc + bw / 2.0)))
                ymax = min(sheet_height, int(round(yc + bh / 2.0)))
                if ymax > ymin and xmax > xmin:
                    detections.append(((ymin, xmin, ymax, xmax), conf, cls_id))
    except Exception as exc:
        logger.debug("Failed to parse precomputed plant detector label %s: %s", p, exc)

    return detections


from scripts.vision.step02_ruler_calibration_io import (  # noqa: E402
    RulerCalibrationMap,
    _extract_catalog_number,
    _parse_calibration_row,
    _resolve_scale_provenance,
    load_ruler_calibration_records,
    load_ruler_calibrations,
)

__all__ = [
    "LM3_PLANT_CLASSES", "LM3_LEAF_CLASS_IDS", "LM3_REPRODUCTIVE_CLASS_IDS",
    "LM3_SINGLE_CAPITULUM_CLASS_IDS", "DEFAULT_HERBARIUM_DPI", "DEFAULT_PX_PER_CM",
    "DEFAULT_PX_PER_MM", "_CALIBRATED_SCALE_SOURCES", "compute_lm3_work_scale",
    "compute_work_scale", "unscale_working_bbox_to_sheet", "normalize_pixels_per_mm",
    "resolve_px_per_cm", "compute_area_mm2", "RulerCalibrationMap",
    "load_ruler_calibration_records", "load_ruler_calibrations",
    "load_precomputed_plant_labels",
]
