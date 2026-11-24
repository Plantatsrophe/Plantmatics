#!/usr/bin/env python3
"""
Step 02 LeafMachine3 data structures, native-DPI crop/reprojection utilities,
and lateral single-capitulum involucre metric extraction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np


@dataclass
class DetectedInstance:
    """Represents a segmented botanical instance on a herbarium sheet."""
    catalog_number: str
    leaf_id: int
    bbox: Tuple[int, int, int, int]  # (ymin, xmin, ymax, xmax) in full-sheet coordinates
    mask: np.ndarray                 # Binary uint8 lamina mask (H, W)
    score: float
    class_id: int
    plant_individual_id: int = 0
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    pixels_per_mm: Optional[float] = None
    is_calibrated_mm: bool = False
    crop_offset: Tuple[int, int] = (0, 0)
    mask_origin: Tuple[int, int] = (0, 0)
    trim_meta: Optional[Dict[str, Any]] = None
    pcd_bbox: Optional[Tuple[int, int, int, int]] = None
    circularity: float = 0.0
    extent: float = 0.0
    aspect_ratio: float = 0.0
    area_mm2: float = np.nan
    mask_path: Optional[str] = None
    contour_path: Optional[str] = None
    petiole_mask: Optional[np.ndarray] = None
    hole_mask: Optional[np.ndarray] = None
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None
    work_scale: float = 1.0
    cf_source: str = "measured_from_ruler"
    archetype_score: float = 0.0
    petiole_width_px: Optional[float] = None
    petiole_width_mm: Optional[float] = None
    blade_length_px: Optional[float] = None
    blade_width_px: Optional[float] = None
    ect_path: Optional[str] = None
    blade_cie_l_norm: float = np.nan
    blade_cie_l_sd: float = np.nan
    blade_cie_l_sheet_delta: float = np.nan
    blade_cie_a_delta: float = np.nan
    blade_cie_b_delta: float = np.nan
    base_angle_deg: float = float("nan")
    apex_angle_deg: float = float("nan")
    lamina_curvature: float = float("nan")
    lamina_trace_length_px: float = float("nan")
    bbox_x1: Optional[int] = None
    bbox_y1: Optional[int] = None
    bbox_x2: Optional[int] = None
    bbox_y2: Optional[int] = None
    blade_cie_l_q85_sheet_delta: float = float("nan")
    basal_zone_cie_l_delta: float = float("nan")
    is_abaxial_pressed_leaf: bool = False
    abaxial_cie_l_sheet_delta: float = float("nan")


@dataclass
class _VoucherFinalizeBundle:
    """Thread-safe per-voucher post-processing result bundle for wave-pipelined execution."""
    catalog_number: str
    image_path: str
    extracted_records: List[Dict[str, Any]]
    failure_record: Optional[Dict[str, Any]]
    repro_record: Optional[Dict[str, Any]] = None
    pheno_record: Optional[Dict[str, Any]] = None
    repro_by_pid: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    pheno_by_pid: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    crowns_by_pid: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    rosette_centroids_by_pid: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    by_plant_petiole: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    contour_coords_by_leaf: Dict[int, List[List[float]]] = field(default_factory=dict)


def crop_leaf_bbox_with_padding(
    image_bgr: np.ndarray,
    bbox: Tuple[int, int, int, int],
    pad_ratio: float = 0.0,
) -> Tuple[np.ndarray, Tuple[int, int, int, int]]:
    """
    Crops a candidate leaf bounding box (ymin, xmin, ymax, xmax) from the native-DPI
    voucher sheet with optional isotropic padding (default 0.0 for tight LeafMachine3
    Stage 1 -> Stage 2 crops), clamped to sheet boundaries.
    Returns (crop_bgr, (crop_y0, crop_x0, crop_y1, crop_x1)).
    """
    sheet_h, sheet_w = image_bgr.shape[:2]
    ymin, xmin, ymax, xmax = [int(v) for v in bbox]
    box_h = max(0, ymax - ymin)
    box_w = max(0, xmax - xmin)
    pad_h = int(round(box_h * pad_ratio))
    pad_w = int(round(box_w * pad_ratio))

    crop_y0 = max(0, ymin - pad_h)
    crop_x0 = max(0, xmin - pad_w)
    crop_y1 = min(sheet_h, ymax + pad_h)
    crop_x1 = min(sheet_w, xmax + pad_w)

    if crop_y1 <= crop_y0 or crop_x1 <= crop_x0:
        return np.zeros((0, 0, 3), dtype=image_bgr.dtype), (crop_y0, crop_x0, crop_y0, crop_x0)

    return image_bgr[crop_y0:crop_y1, crop_x0:crop_x1].copy(), (crop_y0, crop_x0, crop_y1, crop_x1)


def reproject_crop_mask_to_sheet(
    crop_mask: np.ndarray,
    crop_box: Tuple[int, int, int, int],
    sheet_shape: Tuple[int, int],
) -> Tuple[np.ndarray, Tuple[int, int, int, int]]:
    """
    Projects a binary instance mask (crop_h, crop_w) from a padded Stage 1 crop
    back into full-sheet coordinates (sheet_h, sheet_w) and returns
    (sheet_mask, full_sheet_bbox).
    """
    sheet_h, sheet_w = int(sheet_shape[0]), int(sheet_shape[1])
    crop_y0, crop_x0, crop_y1, crop_x1 = [int(v) for v in crop_box]
    crop_y0 = max(0, min(sheet_h, crop_y0))
    crop_x0 = max(0, min(sheet_w, crop_x0))
    crop_y1 = max(crop_y0, min(sheet_h, crop_y1))
    crop_x1 = max(crop_x0, min(sheet_w, crop_x1))

    target_h = crop_y1 - crop_y0
    target_w = crop_x1 - crop_x0
    sheet_mask = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
    if target_h <= 0 or target_w <= 0 or crop_mask is None or crop_mask.size == 0:
        return sheet_mask, (crop_y0, crop_x0, crop_y1, crop_x1)

    mask_bin = (crop_mask > 0).astype(np.uint8) * 255
    if mask_bin.shape[0] != target_h or mask_bin.shape[1] != target_w:
        mask_bin = cv2.resize(mask_bin, (target_w, target_h), interpolation=cv2.INTER_NEAREST)

    sheet_mask[crop_y0:crop_y1, crop_x0:crop_x1] = mask_bin

    ys, xs = np.where(mask_bin > 0)
    if len(ys) > 0 and len(xs) > 0:
        ymin_s = crop_y0 + int(np.min(ys))
        ymax_s = crop_y0 + int(np.max(ys)) + 1
        xmin_s = crop_x0 + int(np.min(xs))
        xmax_s = crop_x0 + int(np.max(xs)) + 1
        full_bbox = (ymin_s, xmin_s, min(sheet_h, ymax_s), min(sheet_w, xmax_s))
    else:
        full_bbox = (crop_y0, crop_x0, crop_y1, crop_x1)

    return sheet_mask, full_bbox


def compute_capitulum_metrics(
    mask: Optional[np.ndarray],
    bbox: Tuple[int, int, int, int],
    pixels_per_mm: Optional[float] = None,
    min_ar: float = 0.25,
    max_ar: float = 4.5,
    crop_bgr: Optional[np.ndarray] = None,
    class_id: Optional[int] = None,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
    is_calibrated_mm: bool = True,
) -> Optional[Dict[str, Any]]:
    """
    Computes herbaceous involucre height (H_inv), width (W_inv), and aspect ratio (AR_inv = H/W)
    strictly for single-capitulum detections (`Seed_Fruit_ONE = 2`, `Flower_ONE = 4`, `Bud = 6`).

    Anatomical & Geometric Safeguards:
      1. Rejects multi-head corymbose cluster detections (`Seed_Fruit_MANY = 3`, `Flower_MANY = 5`).
      2. When `crop_bgr` is supplied, isolates the green/olive/anthocyanic herbaceous phyllary cup
         by masking out sheet paper background, yellow/gold ray/disc corollas, and bright white
         capillary pappus plumes in HSV/RGB space.
      3. Resolves longitudinal `involucre_height_px` vs. transverse `involucre_width_px` along the
         anatomical peduncle-to-apex vector (from phyllary cup centroid to distal floret/pappus
         centroid when present, else `peduncle_dir_xy`, else upward sheet axis), avoiding gravity
         axis inversions on tilted or nodding peduncles.
      4. Uses wide artifact-only bounds (`0.25 <= AR_inv <= 4.5`) so genuine campanulate/hemispheric
         fruiting involucres and slender cylindrical buds are never censored.
    """
    if class_id is not None and int(class_id) not in (2, 4, 6):
        return None

    eff_mask: Optional[np.ndarray] = None
    if mask is not None and isinstance(mask, np.ndarray) and mask.size > 0:
        eff_mask = (mask > 0).astype(np.uint8) * 255

    axis_ux, axis_uy = 0.0, -1.0
    has_anatomical_axis = False
    if peduncle_dir_xy is not None and len(peduncle_dir_xy) >= 2:
        px_d, py_d = float(peduncle_dir_xy[0]), float(peduncle_dir_xy[1])
        p_norm = float(np.hypot(px_d, py_d))
        if np.isfinite(p_norm) and p_norm > 1e-6:
            axis_ux, axis_uy = px_d / p_norm, py_d / p_norm
            has_anatomical_axis = True

    if crop_bgr is not None and isinstance(crop_bgr, np.ndarray) and crop_bgr.ndim == 3 and crop_bgr.shape[0] >= 5 and crop_bgr.shape[1] >= 5:
        ch, cw = crop_bgr.shape[:2]
        hsv = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2HSV)
        H_ch = hsv[:, :, 0]
        S_ch = hsv[:, :, 1]
        V_ch = hsv[:, :, 2]

        border_pts = np.concatenate([
            crop_bgr[0, :, :],
            crop_bgr[-1, :, :],
            crop_bgr[:, 0, :],
            crop_bgr[:, -1, :],
        ], axis=0).astype(np.float32)
        bg_med = np.median(border_pts, axis=0)
        bg_dist = np.linalg.norm(crop_bgr.astype(np.float32) - bg_med[None, None, :], axis=2)
        paper_mask = (bg_dist < 18.0) & (S_ch <= 42)

        fg_plant = ~paper_mask
        if eff_mask is not None and eff_mask.shape[:2] == (ch, cw) and np.any(eff_mask > 0):
            fg_plant = fg_plant & (eff_mask > 0)

        floret_mask = fg_plant & (H_ch >= 10) & (H_ch <= 36) & (S_ch >= 76) & (V_ch >= 85)
        pappus_mask = fg_plant & (~floret_mask) & (V_ch >= 158) & (S_ch <= 48)

        # Resolve distal half along anatomical axis and exclude distal aged ochre/amber corollas on non-Bud heads
        fg_ys_all, fg_xs_all = np.where(fg_plant)
        if fg_xs_all.size > 0 and (class_id is None or int(class_id) in (2, 4)):
            prelim_distal = floret_mask | pappus_mask
            prelim_cup = fg_plant & (~prelim_distal) & (V_ch <= 198)
            dist_ax_ux, dist_ax_uy = axis_ux, axis_uy
            if (
                not has_anatomical_axis
                and int(np.count_nonzero(prelim_distal)) >= 8
                and int(np.count_nonzero(prelim_cup)) >= 12
            ):
                py_c, px_c = np.where(prelim_cup)
                py_d, px_d = np.where(prelim_distal)
                dx_p = float(np.mean(px_d) - np.mean(px_c))
                dy_p = float(np.mean(py_d) - np.mean(py_c))
                dn_p = float(np.hypot(dx_p, dy_p))
                if dn_p >= 2.0:
                    dist_ax_ux, dist_ax_uy = dx_p / dn_p, dy_p / dn_p
            proj_fg = fg_xs_all.astype(np.float32) * dist_ax_ux + fg_ys_all.astype(np.float32) * dist_ax_uy
            proj_mid = 0.5 * (float(np.min(proj_fg)) + float(np.max(proj_fg)))
            yy_g, xx_g = np.ogrid[:ch, :cw]
            distal_half = (xx_g.astype(np.float32) * dist_ax_ux + yy_g.astype(np.float32) * dist_ax_uy) >= proj_mid
            aged_ochre_mask = (
                fg_plant
                & distal_half
                & (H_ch >= 8)
                & (H_ch <= 36)
                & (S_ch >= 38)
                & (S_ch <= 75)
                & (V_ch >= 75)
            )
            floret_mask = floret_mask | aged_ochre_mask
            pappus_mask = fg_plant & (~floret_mask) & (V_ch >= 158) & (S_ch <= 48)

        involucre_cand = fg_plant & (~floret_mask) & (~pappus_mask) & (V_ch <= 198)

        min_inv_px = max(12, int(0.04 * ch * cw))
        if int(np.count_nonzero(involucre_cand)) >= min_inv_px:
            eff_mask = involucre_cand.astype(np.uint8) * 255
            distal_mask = floret_mask | pappus_mask
            if int(np.count_nonzero(distal_mask)) >= 8:
                iy_inv, ix_inv = np.where(involucre_cand)
                iy_dst, ix_dst = np.where(distal_mask)
                dx_a = float(np.mean(ix_dst) - np.mean(ix_inv))
                dy_a = float(np.mean(iy_dst) - np.mean(iy_inv))
                d_norm = float(np.hypot(dx_a, dy_a))
                if d_norm >= 2.0:
                    axis_ux, axis_uy = dx_a / d_norm, dy_a / d_norm
                    has_anatomical_axis = True

        # Collect viewpoint symmetry cues for top-down radial disc detection
        n_floret_px = int(np.count_nonzero(floret_mask))
        n_inv_px = int(np.count_nonzero(involucre_cand))
        topdown_radial_4q = False
        topdown_radial_3q = False
        if fg_xs_all.size >= 24 and n_floret_px >= 12 and n_inv_px >= 12:
            fg_cx = float(np.mean(fg_xs_all))
            fg_cy = float(np.mean(fg_ys_all))
            fg_w = float(int(fg_xs_all.max()) - int(fg_xs_all.min()) + 1)
            fg_h = float(int(fg_ys_all.max()) - int(fg_ys_all.min()) + 1)
            head_diam = max(fg_w, fg_h, 1.0)
            fg_bbox_ar = min(fg_w, fg_h) / head_diam

            fl_y, fl_x = np.where(floret_mask)
            inv_y, inv_x = np.where(involucre_cand)
            fl_cx, fl_cy = float(np.mean(fl_x)), float(np.mean(fl_y))
            inv_cx, inv_cy = float(np.mean(inv_x)), float(np.mean(inv_y))

            fl_offset_frac = float(np.hypot(fl_cx - fg_cx, fl_cy - fg_cy)) / head_diam
            inv_fl_offset_frac = float(np.hypot(fl_cx - inv_cx, fl_cy - inv_cy)) / head_diam

            def _occupied_quadrants(xs_a: np.ndarray, ys_a: np.ndarray, ref_x: float, ref_y: float) -> int:
                dx_q = xs_a.astype(np.float32) - ref_x
                dy_q = ys_a.astype(np.float32) - ref_y
                min_q = max(4, int(0.06 * xs_a.size))
                return int(
                    (int(np.count_nonzero((dx_q >= 0.0) & (dy_q >= 0.0))) >= min_q)
                    + (int(np.count_nonzero((dx_q < 0.0) & (dy_q >= 0.0))) >= min_q)
                    + (int(np.count_nonzero((dx_q >= 0.0) & (dy_q < 0.0))) >= min_q)
                    + (int(np.count_nonzero((dx_q < 0.0) & (dy_q < 0.0))) >= min_q)
                )

            n_inv_quads = _occupied_quadrants(inv_x, inv_y, fl_cx, fl_cy)
            n_fl_quads = _occupied_quadrants(fl_x, fl_y, inv_cx, inv_cy)
            r_fl_q75 = float(np.percentile(np.hypot(fl_x - fg_cx, fl_y - fg_cy), 75))
            r_inv_q75 = float(np.percentile(np.hypot(inv_x - fg_cx, inv_y - fg_cy), 75))
            outer_lig_ring_4q = bool(n_fl_quads == 4 and r_fl_q75 > 1.15 * max(r_inv_q75, 1.0))
            outer_lig_ring_3q = bool(n_fl_quads >= 3 and r_fl_q75 > 1.15 * max(r_inv_q75, 1.0))

            if (
                fg_bbox_ar >= 0.72
                and fl_offset_frac <= 0.18
                and inv_fl_offset_frac <= 0.18
                and (n_inv_quads == 4 or outer_lig_ring_4q)
            ):
                topdown_radial_4q = True
            if fl_offset_frac <= 0.22 and (n_inv_quads >= 3 or outer_lig_ring_3q):
                topdown_radial_3q = True
    else:
        topdown_radial_4q = False
        topdown_radial_3q = False

    h_inv_px = 0.0
    w_inv_px = 0.0
    contours = []
    if eff_mask is not None and eff_mask.size > 0:
        contours, _ = cv2.findContours(eff_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours and cv2.contourArea(max(contours, key=cv2.contourArea)) >= 10:
        main_contour = max(contours, key=cv2.contourArea)
        rect = cv2.minAreaRect(main_contour)
        box = cv2.boxPoints(rect)
        v1 = box[1] - box[0]
        v2 = box[2] - box[1]
        len1 = float(np.linalg.norm(v1))
        len2 = float(np.linalg.norm(v2))
        norm1 = max(len1, 1e-6)
        norm2 = max(len2, 1e-6)
        proj1 = abs(float(v1[0]) * axis_ux + float(v1[1]) * axis_uy) / norm1
        proj2 = abs(float(v2[0]) * axis_ux + float(v2[1]) * axis_uy) / norm2
        if proj1 >= proj2:
            h_inv_px, w_inv_px = len1, len2
        else:
            h_inv_px, w_inv_px = len2, len1
    else:
        ymin, xmin, ymax, xmax = bbox
        dy_box = float(max(ymax - ymin, 1))
        dx_box = float(max(xmax - xmin, 1))
        if has_anatomical_axis and abs(axis_ux) > abs(axis_uy):
            h_inv_px, w_inv_px = dx_box, dy_box
        else:
            h_inv_px, w_inv_px = dy_box, dx_box

    if w_inv_px <= 1e-6:
        return None

    aspect_ratio = h_inv_px / w_inv_px
    is_lateral_profile = not (topdown_radial_4q or (topdown_radial_3q and aspect_ratio < 0.52))
    if not is_lateral_profile:
        return {
            "involucre_height_px": np.nan,
            "involucre_width_px": np.nan,
            "involucre_height_mm": np.nan,
            "involucre_width_mm": np.nan,
            "AR_inv": np.nan,
            "is_lateral_profile": False,
        }

    if not (min_ar <= aspect_ratio <= max_ar):
        return None

    use_cal = bool(is_calibrated_mm and pixels_per_mm is not None and not np.isnan(pixels_per_mm) and float(pixels_per_mm) > 0.0)
    h_inv_mm = (h_inv_px / float(pixels_per_mm)) if use_cal else np.nan
    w_inv_mm = (w_inv_px / float(pixels_per_mm)) if use_cal else np.nan

    return {
        "involucre_height_px": round(h_inv_px, 2),
        "involucre_width_px": round(w_inv_px, 2),
        "involucre_height_mm": round(h_inv_mm, 3) if not np.isnan(h_inv_mm) else np.nan,
        "involucre_width_mm": round(w_inv_mm, 3) if not np.isnan(w_inv_mm) else np.nan,
        "AR_inv": round(aspect_ratio, 4),
        "is_lateral_profile": True,
    }
