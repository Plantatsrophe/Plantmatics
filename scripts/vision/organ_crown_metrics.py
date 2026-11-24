#!/usr/bin/env python3
"""
===============================================================================
Module: organ_crown_metrics.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Thin backward-compatible facade for Step 02 multi-plant DBSCAN spatial
    clustering, 31-keypoint rosette crown / petiole / bent-scape linear trait
    estimation, and 4-case vascular-topology basal vs. cauline leaf zone
    classification on LeafMachine3 outputs.
===============================================================================
"""

from __future__ import annotations

from typing import Tuple

import cv2
import numpy as np

from scripts.vision.lm3_compat import DBSCAN


def _mask_basal_apical_fallback(mask: np.ndarray) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Minimal contour endpoint fallback `(p_apex, p_base)` when bare 2D masks lack LM3 pose keypoints."""
    if mask is None or mask.size == 0:
        return (0, 0), (0, 0)
    m_u8 = (mask > 0).astype(np.uint8) * 255
    cnts, _ = cv2.findContours(m_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        return (0, 0), (0, 0)
    pts = max(cnts, key=cv2.contourArea).reshape(-1, 2)
    base_idx = int(np.argmax(pts[:, 1]))
    pb = (int(pts[base_idx, 0]), int(pts[base_idx, 1]))
    apex_idx = int(np.argmax(np.hypot(pts[:, 0] - pb[0], pts[:, 1] - pb[1])))
    pa = (int(pts[apex_idx, 0]), int(pts[apex_idx, 1]))
    return pa, pb


from scripts.vision.step02_leaf_zone_classifier import classify_leaf_zone_anatomical
from scripts.vision.step02_petiole_scape_metrics import approximate_petiole_and_scape_metrics
from scripts.vision.step02_plant_cluster_merge import cluster_plant_individuals
from scripts.vision.step02_rosette_crown_pose import (
    _extract_xy_point,
    _parse_basal_leaf_pose,
)
from scripts.vision.step02_rosette_crown_solver import (
    _estimate_single_plant_crown,
    estimate_rosette_crown_point,
)

__all__ = [
    "DBSCAN",
    "_mask_basal_apical_fallback",
    "cluster_plant_individuals",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
    "classify_leaf_zone_anatomical",
]
