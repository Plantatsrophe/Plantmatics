#!/usr/bin/env python3
"""
Step 02 early capitulum phenology biomarkers and upright crop orientation:
pappus plume extraction, active yellow corolla detection, involucre posture
analysis, and peduncle/biomarker-guided upright reorientation.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

import cv2
import numpy as np


def extract_biomarker_pappus_plume(
    crop_rgb: np.ndarray,
    upper_pct: float = 0.45,
    min_l_star: float = 75.0,
    max_s: float = 0.18,
    min_laplacian_var: float = 40.0,
    min_edge_density: float = 0.03,
    plume_threshold: float = 0.18,
) -> Tuple[bool, float, Dict[str, float]]:
    """
    Biomarker A: Expanded Capillary Pappus Plume (Fruit Indicator).

    Extracts high-luminance, low-saturation pixels in the upper ROI of the capitulum
    dome, confirmed by high-frequency texture (Laplacian variance / Canny edge density)
    to differentiate fine capillary bristles from plain white mounting paper.

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image crop of a detected capitulum.
    upper_pct : float
        Fraction of upper vertical height representing the apical dome (default 0.45).
    min_l_star : float
        Minimum CIELAB lightness (L* >= 75.0) for bright pappus bristles.
    max_s : float
        Maximum HSV saturation (S <= 0.18) for neutral white/cream bristles.
    min_laplacian_var : float
        Threshold for Laplacian variance confirming fibrous bristle micro-texture.
    min_edge_density : float
        Threshold for Canny edge density confirming fine structural fibers.
    plume_threshold : float
        Threshold ratio of white fibrous pixels in upper ROI to assert plume presence (18%).

    Returns
    -------
    Tuple[bool, float, Dict[str, float]]
        - has_pappus_plume (bool)
        - white_fibrous_ratio (float, 0.0 to 1.0)
        - diagnostics (dict containing intermediate ratios and texture metrics)
    """
    h, w = crop_rgb.shape[:2]
    if h < 3 or w < 3:
        return False, 0.0, {
            "white_ratio": 0.0,
            "white_fibrous_ratio": 0.0,
            "laplacian_var": 0.0,
            "edge_density": 0.0,
        }

    upper_h = max(1, int(round(h * upper_pct)))
    upper_roi = crop_rgb[:upper_h, :]
    n_upper = max(1, upper_roi.shape[0] * upper_roi.shape[1])

    # Convert upper dome to LAB and HSV
    upper_lab = cv2.cvtColor(upper_roi, cv2.COLOR_RGB2LAB)
    upper_hsv = cv2.cvtColor(upper_roi, cv2.COLOR_RGB2HSV)
    upper_gray = cv2.cvtColor(upper_roi, cv2.COLOR_RGB2GRAY)

    # L* in OpenCV uint8 is scaled: L* = L_cv2 * (100.0 / 255.0)
    l_star = upper_lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    s_norm = upper_hsv[:, :, 1].astype(np.float32) / 255.0

    white_mask = (l_star >= min_l_star) & (s_norm <= max_s)
    white_pixel_ratio = float(np.count_nonzero(white_mask) / n_upper)

    # High-frequency texture analysis
    laplacian = cv2.Laplacian(upper_gray, cv2.CV_64F)
    laplacian_var = float(laplacian.var())

    edges = cv2.Canny(upper_gray, 50, 150)
    edge_count = np.count_nonzero(edges)
    edge_density = float(edge_count / n_upper)

    has_texture = (laplacian_var >= min_laplacian_var) or (edge_density >= min_edge_density)

    # Differentiate textured white bristles from uniform white mounting paper:
    # Dilate edges slightly (3x3 kernel) to encompass fibrous bristle shafts and gradients
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
    dilated_edges = cv2.dilate(edges, kernel) > 0
    fibrous_white_mask = white_mask & (dilated_edges | (np.abs(laplacian) > 15.0))
    fibrous_white_count = np.count_nonzero(fibrous_white_mask)
    white_fibrous_ratio = float(fibrous_white_count / n_upper)

    has_pappus_plume = bool((white_fibrous_ratio > plume_threshold) and has_texture)

    diagnostics = {
        "white_ratio": round(white_pixel_ratio, 4),
        "white_fibrous_ratio": round(white_fibrous_ratio, 4),
        "laplacian_var": round(laplacian_var, 2),
        "edge_density": round(edge_density, 4),
    }

    return has_pappus_plume, float(round(white_fibrous_ratio, 4)), diagnostics


def extract_biomarker_yellow_corollas(
    crop_rgb: np.ndarray,
    h_min_deg: float = 18.0,
    h_max_deg: float = 45.0,
    s_min: float = 0.32,
    v_min: float = 0.35,
    corolla_threshold: float = 0.12,
) -> Tuple[bool, float, Dict[str, float]]:
    """
    Biomarker B: Active Yellow Corollas / Ligules (Anthesis Indicator).

    Isolates the bright golden-yellow hue band across the capitulum disc and ray margins:
      18° <= Hue <= 45°, Saturation >= 0.32, Value >= 0.35.
    Saturated yellow area > 12% indicates fresh, active anthesis.

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image crop of a detected capitulum.
    h_min_deg : float
        Minimum hue in degrees (default 18.0°).
    h_max_deg : float
        Maximum hue in degrees (default 45.0°).
    s_min : float
        Minimum saturation (default 0.32).
    v_min : float
        Minimum brightness value (default 0.35).
    corolla_threshold : float
        Proportion of saturated yellow pixels to confirm fresh corollas (default 0.12).

    Returns
    -------
    Tuple[bool, float, Dict[str, float]]
        - has_fresh_yellow_corollas (bool)
        - yellow_ratio (float, 0.0 to 1.0)
        - diagnostics (dict containing yellow ratio and thresholds)
    """
    h, w = crop_rgb.shape[:2]
    n_total = max(1, h * w)

    if h < 3 or w < 3:
        return False, 0.0, {"yellow_ratio": 0.0}

    crop_hsv = cv2.cvtColor(crop_rgb, cv2.COLOR_RGB2HSV)

    # In OpenCV uint8, H is in [0, 179], representing 0° to 360° (H_deg = H_cv2 * 2.0)
    h_deg = crop_hsv[:, :, 0].astype(np.float32) * 2.0
    s_norm = crop_hsv[:, :, 1].astype(np.float32) / 255.0
    v_norm = crop_hsv[:, :, 2].astype(np.float32) / 255.0

    yellow_mask = (
        (h_deg >= h_min_deg)
        & (h_deg <= h_max_deg)
        & (s_norm >= s_min)
        & (v_norm >= v_min)
    )

    yellow_pixels = np.count_nonzero(yellow_mask)
    yellow_ratio = float(yellow_pixels / n_total)
    has_fresh_yellow_corollas = bool(yellow_ratio > corolla_threshold)

    diagnostics = {
        "yellow_ratio": round(yellow_ratio, 4),
        "yellow_pixel_count": int(yellow_pixels),
    }

    return has_fresh_yellow_corollas, float(round(yellow_ratio, 4)), diagnostics


def extract_biomarker_involucre_posture(
    crop_rgb: np.ndarray,
    lower_pct: float = 0.45,
    min_aspect_ratio: float = 1.0,
    max_reflexed_ar: float = 0.75,
    roughness_thresh: float = 1.15,
) -> Tuple[bool, bool, Dict[str, float]]:
    """
    Biomarker C: Involucre Posture (Lower Contour Profile).

    Segments the herbaceous lower involucre from background mounting paper and
    characterizes the basal profile:
      - Erect involucre: Smooth, tapered U-shape (H_inv / W_inv >= 1.0, low roughness).
      - Reflexed involucre: Flattened or flared downward (H_inv / W_inv <= 0.75, high boundary roughness).

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image crop of a detected capitulum.
    lower_pct : float
        Vertical split point starting the lower involucre (default 0.45).
    min_aspect_ratio : float
        Minimum H_inv / W_inv for an erect cylindrical/campanulate involucre (1.0).
    max_reflexed_ar : float
        Maximum H_inv / W_inv for reflexed/flaring phyllaries (0.75).
    roughness_thresh : float
        Contour perimeter / convex hull perimeter threshold for jagged reflexed phyllaries (1.15).

    Returns
    -------
    Tuple[bool, bool, Dict[str, float]]
        - erect_involucre (bool)
        - reflexed_involucre (bool)
        - diagnostics (dict with H_inv, W_inv, aspect_ratio, boundary_roughness)
    """
    h, w = crop_rgb.shape[:2]
    if h < 4 or w < 4:
        return True, False, {
            "h_inv_px": 0.0,
            "w_inv_px": 0.0,
            "involucre_aspect_ratio": 1.0,
            "boundary_roughness": 1.0,
        }

    lower_y = int(round(h * lower_pct))
    lower_roi = crop_rgb[lower_y:, :]
    lh, lw = lower_roi.shape[:2]
    if lh < 2 or lw < 2:
        return True, False, {
            "h_inv_px": float(lh),
            "w_inv_px": float(lw),
            "involucre_aspect_ratio": 1.0,
            "boundary_roughness": 1.0,
        }

    lower_lab = cv2.cvtColor(lower_roi, cv2.COLOR_RGB2LAB)
    lower_hsv = cv2.cvtColor(lower_roi, cv2.COLOR_RGB2HSV)

    l_star = lower_lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    s_norm = lower_hsv[:, :, 1].astype(np.float32) / 255.0
    h_deg = lower_hsv[:, :, 0].astype(np.float32) * 2.0

    # Plant tissue segmentation: filter out high-luminance neutral mounting paper
    is_paper = (l_star >= 75.0) & (s_norm <= 0.20)
    is_plant = ~is_paper

    # Clean morphological structure
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    plant_mask = cv2.morphologyEx(is_plant.astype(np.uint8) * 255, cv2.MORPH_CLOSE, kernel)

    contours, _ = cv2.findContours(plant_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    if contours and cv2.contourArea(max(contours, key=cv2.contourArea)) >= 10:
        main_c = max(contours, key=cv2.contourArea)
        bx, by, bw, bh = cv2.boundingRect(main_c)
        h_inv = float(bh)
        w_inv = float(max(1, bw))
        ar_inv = h_inv / w_inv

        peri = float(cv2.arcLength(main_c, True))
        hull = cv2.convexHull(main_c)
        hull_peri = float(cv2.arcLength(hull, True))
        roughness = float(peri / max(1e-6, hull_peri))

        erect_involucre = bool(ar_inv >= min_aspect_ratio or (ar_inv > max_reflexed_ar and roughness < roughness_thresh))
        reflexed_involucre = bool(ar_inv <= max_reflexed_ar and roughness >= roughness_thresh)
    else:
        # Fallback to lower ROI bounding aspect ratio
        h_inv = float(lh)
        w_inv = float(max(1, lw))
        ar_inv = h_inv / w_inv
        roughness = 1.0
        erect_involucre = bool(ar_inv >= 0.85)
        reflexed_involucre = bool(ar_inv <= 0.65)

    diagnostics = {
        "h_inv_px": round(h_inv, 2),
        "w_inv_px": round(w_inv, 2),
        "involucre_aspect_ratio": round(ar_inv, 4),
        "boundary_roughness": round(roughness, 4),
    }

    return erect_involucre, reflexed_involucre, diagnostics


def orient_capitulum_crop_upright(
    crop_rgb: np.ndarray,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
) -> np.ndarray:
    """
    Reorients a capitulum RGB crop so its basal-to-distal axis points upward.

    Uses `peduncle_dir_xy` when supplied and falls back to crop biomarker centroids
    between distal florets or pappus bristles and the herbaceous involucre cup.
    """
    if (
        crop_rgb is None
        or not isinstance(crop_rgb, np.ndarray)
        or crop_rgb.size == 0
        or crop_rgb.ndim < 3
        or crop_rgb.shape[2] < 3
        or min(crop_rgb.shape[:2]) < 5
    ):
        return crop_rgb

    h, w = crop_rgb.shape[:2]
    dx: Optional[float] = None
    dy: Optional[float] = None

    if peduncle_dir_xy is not None and len(peduncle_dir_xy) >= 2:
        try:
            px = float(peduncle_dir_xy[0])
            py = float(peduncle_dir_xy[1])
        except (TypeError, ValueError):
            px, py = float("nan"), float("nan")
        if np.isfinite(px) and np.isfinite(py):
            nrm = math.hypot(px, py)
            if nrm >= 1e-3:
                dx = px / nrm
                dy = py / nrm

    if dx is None or dy is None:
        rgb3 = crop_rgb[:, :, :3]
        lab = cv2.cvtColor(rgb3, cv2.COLOR_RGB2LAB)
        hsv = cv2.cvtColor(rgb3, cv2.COLOR_RGB2HSV)
        gray = cv2.cvtColor(rgb3, cv2.COLOR_RGB2GRAY)

        l_star = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
        b_star = lab[:, :, 2].astype(np.float32) - 128.0
        h_deg = hsv[:, :, 0].astype(np.float32) * 2.0
        s_norm = hsv[:, :, 1].astype(np.float32) / 255.0
        v_norm = hsv[:, :, 2].astype(np.float32) / 255.0

        non_paper = ~((l_star >= 75.0) & (s_norm <= 0.20))
        floret_m = (
            (h_deg >= 10.0)
            & (h_deg <= 50.0)
            & (s_norm >= 0.18)
            & (v_norm >= 0.28)
            & (b_star >= 8.0)
        )
        edges = cv2.Canny(gray, 50, 150)
        dil_edges = cv2.dilate(edges, cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))) > 0
        pappus_m = (l_star >= 75.0) & (s_norm <= 0.18) & dil_edges
        distal_m = floret_m | pappus_m
        cup_m = non_paper & (~floret_m)

        if np.count_nonzero(distal_m) >= 8 and np.count_nonzero(cup_m) >= 8:
            ys_d, xs_d = np.where(distal_m)
            ys_c, xs_c = np.where(cup_m)
            cy_d, cx_d = float(np.mean(ys_d)), float(np.mean(xs_d))
            cy_c, cx_c = float(np.mean(ys_c)), float(np.mean(xs_c))
            vx = cx_d - cx_c
            vy = cy_d - cy_c
            if math.hypot(vx, vy) >= 3.0:
                dx = vx
                dy = vy

    if dx is None or dy is None:
        return crop_rgb

    theta_vec = math.degrees(math.atan2(dy, dx))
    rot_deg = ((theta_vec + 90.0 + 180.0) % 360.0) - 180.0
    if abs(rot_deg) <= 25.0:
        return crop_rgb

    rad = math.radians(rot_deg)
    abs_cos = abs(math.cos(rad))
    abs_sin = abs(math.sin(rad))
    new_w = max(1, int(round(h * abs_sin + w * abs_cos)))
    new_h = max(1, int(round(h * abs_cos + w * abs_sin)))
    rot_mat = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), rot_deg, 1.0)
    rot_mat[0, 2] += (new_w / 2.0) - (w / 2.0)
    rot_mat[1, 2] += (new_h / 2.0) - (h / 2.0)
    return cv2.warpAffine(
        crop_rgb,
        rot_mat,
        (new_w, new_h),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=(240, 240, 240),
    )


__all__ = [
    "extract_biomarker_pappus_plume",
    "extract_biomarker_yellow_corollas",
    "extract_biomarker_involucre_posture",
    "orient_capitulum_crop_upright",
]
