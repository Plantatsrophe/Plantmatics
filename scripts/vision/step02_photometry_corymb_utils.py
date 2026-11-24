#!/usr/bin/env python3
"""
Step 02 calibrated CIELAB blade photometry (narrow-segment safe, wool-safe strap
exclusion, midrib corridor exclusion) and corymb cluster capitulum deconvolution
utilities (`scripts/vision/step02_photometry_corymb_utils.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Tuple

import cv2
import numpy as np


def resolve_effective_px_per_mm(
    pixels_per_mm: Optional[float],
    sheet_shape_hw: Tuple[int, int],
    fallback_px_per_mm: float = 10.0,
) -> Tuple[float, bool]:
    """Resolve active `(px/mm, is_empirical)` using verified ruler calibration or LM3 megapixel prior."""
    if pixels_per_mm is not None:
        try:
            ppm = float(pixels_per_mm)
            if np.isfinite(ppm) and ppm > 0.0:
                return (ppm, True)
        except (TypeError, ValueError):
            pass
    h, w = (max(0, int(sheet_shape_hw[0])), max(0, int(sheet_shape_hw[1]))) if (sheet_shape_hw and len(sheet_shape_hw) >= 2) else (0, 0)
    mp = (float(h) * float(w)) / 1e6
    return (float(2.7112084 * math.sqrt(mp)), False) if mp >= 0.5 else (float(fallback_px_per_mm), False)


def estimate_sheet_paper_cie_l_ref_and_ab(image_bgr: Optional[np.ndarray]) -> Tuple[float, float, float]:
    """Estimate mounting-paper CIE L*, a*, b* reference from an interior 3%..15% border band."""
    if image_bgr is None or not isinstance(image_bgr, np.ndarray) or image_bgr.ndim != 3 or image_bgr.shape[2] < 3 or image_bgr.size == 0:
        return (90.0, 0.0, 2.0)
    sample = image_bgr[::8, ::8, :3]
    if sample.size == 0:
        sample = image_bgr[:, :, :3]
    lab = cv2.cvtColor(sample, cv2.COLOR_BGR2LAB)
    l_cie = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    a_c, b_c = lab[:, :, 1].astype(np.float32) - 128.0, lab[:, :, 2].astype(np.float32) - 128.0
    chroma = np.hypot(a_c, b_c)

    sh, sw = l_cie.shape
    if sh >= 24 and sw >= 24:
        r0, r1 = max(1, int(round(0.03 * sh))), max(3, int(round(0.15 * sh)))
        c0, c1 = max(1, int(round(0.03 * sw))), max(3, int(round(0.15 * sw)))
        yy, xx = np.ogrid[:sh, :sw]
        border_band = ((yy >= r0) & (yy < sh - r0) & (xx >= c0) & (xx < sw - c0)) & ~((yy >= r1) & (yy < sh - r1) & (xx >= c1) & (xx < sw - c1))
    else:
        border_band = np.ones((sh, sw), dtype=bool)
    paper_mask = border_band & (l_cie >= 70.0) & (l_cie <= 99.0) & (chroma <= 12.0)
    if int(np.count_nonzero(paper_mask)) < 32:
        paper_mask = (l_cie >= 65.0) & (l_cie <= 99.0) & (chroma <= 18.0)
    if int(np.count_nonzero(paper_mask)) >= 16:
        l_ref, a_ref, b_ref = float(np.percentile(l_cie[paper_mask], 90.0)), float(np.median(a_c[paper_mask])), float(np.median(b_c[paper_mask]))
    else:
        l_ref, a_ref, b_ref = float(np.percentile(l_cie, 85.0)), 0.0, 2.0
    return (float(np.clip(l_ref, 65.0, 99.0)), a_ref, b_ref)


def estimate_sheet_paper_cie_l_ref(image_bgr: Optional[np.ndarray]) -> float:
    """Estimate mounting-paper CIE L* white-point reference (0..100 scale)."""
    return estimate_sheet_paper_cie_l_ref_and_ab(image_bgr)[0]


def _cie_l_to_rel_y(l_val: float) -> float:
    return ((l_val + 16.0) / 116.0) ** 3 if l_val > 8.0 else l_val / 903.3


def _parse_kp_xy(raw: Any) -> Optional[Tuple[float, float]]:
    if raw is None:
        return None
    try:
        if isinstance(raw, Mapping):
            conf = raw.get("conf", raw.get("confidence"))
            x, y = float(raw["x"]), float(raw["y"])
        elif isinstance(raw, (tuple, list, np.ndarray)) and len(raw) >= 2:
            conf = raw[2] if len(raw) >= 3 else None
            x, y = float(raw[0]), float(raw[1])
        else:
            return None
        if conf is not None and float(conf) < 0.25:
            return None
        return (x, y) if (math.isfinite(x) and math.isfinite(y)) else None
    except (KeyError, TypeError, ValueError):
        return None


def _component_intrinsic_dog_mad(l_cie: np.ndarray, comp: np.ndarray, sigma_fine: float, sigma_coarse: float, scale_norm: float) -> float:
    """Measure edge-bleed-free DoG MAD texture strictly inside a candidate strap component."""
    ys, xs = np.where(comp)
    pad, (h, w) = max(2, int(math.ceil(2.0 * sigma_coarse))), l_cie.shape[:2]
    y0, y1 = max(0, int(ys.min()) - pad), min(h, int(ys.max()) + pad + 1)
    x0, x1 = max(0, int(xs.min()) - pad), min(w, int(xs.max()) + pad + 1)
    comp_sub, l_sub = comp[y0:y1, x0:x1], l_cie[y0:y1, x0:x1]
    filled = np.where(comp_sub, l_sub, float(np.median(l_sub[comp_sub]))).astype(np.float32)
    bp = (cv2.GaussianBlur(filled, (0, 0), sigmaX=sigma_fine, sigmaY=sigma_fine) - cv2.GaussianBlur(filled, (0, 0), sigmaX=sigma_coarse, sigmaY=sigma_coarse))[comp_sub]
    return 1.4826 * float(np.median(np.abs(bp - float(np.median(bp))))) * scale_norm


def compute_calibrated_blade_photometry(
    sheet_bgr: Optional[np.ndarray],
    blade_mask: Optional[np.ndarray],
    mask_origin: Tuple[int, int] = (0, 0),
    paper_l_ref: Optional[float] = None,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
    pixels_per_mm: Optional[float] = None,
    keypoints_31: Optional[Mapping[str, Any]] = None,
) -> Dict[str, float]:
    """Compute calibrated lamina CIELAB photometry with narrow-segment erosion, midrib masking, and wool-safe strap exclusion."""
    nan_out: Dict[str, float] = {
        "blade_cie_l_norm": np.nan, "blade_cie_l_sd": np.nan,
        "blade_cie_l_sheet_delta": np.nan, "blade_cie_l_q85_sheet_delta": np.nan,
        "basal_zone_cie_l_delta": np.nan, "blade_cie_a_delta": np.nan, "blade_cie_b_delta": np.nan,
    }
    if (
        sheet_bgr is None or blade_mask is None
        or not isinstance(sheet_bgr, np.ndarray) or not isinstance(blade_mask, np.ndarray)
        or sheet_bgr.ndim != 3 or sheet_bgr.shape[2] < 3 or blade_mask.ndim < 2 or blade_mask.size == 0
    ):
        return nan_out

    sh, sw = sheet_bgr.shape[:2]
    mh, mw = blade_mask.shape[:2]
    ox, oy = int(mask_origin[0]), int(mask_origin[1])
    x0, y0, x1, y1 = max(0, ox), max(0, oy), min(sw, ox + mw), min(sh, oy + mh)
    if x1 <= x0 or y1 <= y0:
        return nan_out

    mx0, my0 = x0 - ox, y0 - oy
    mx1, my1 = mx0 + (x1 - x0), my0 + (y1 - y0)
    crop_bgr = sheet_bgr[y0:y1, x0:x1, :3]
    lam_slice = blade_mask[my0:my1, mx0:mx1]
    lam = (lam_slice[:, :, 0] > 0) if lam_slice.ndim == 3 else (lam_slice > 0)

    for excl_mask in (hole_mask, petiole_mask):
        if excl_mask is not None and isinstance(excl_mask, np.ndarray) and excl_mask.shape[:2] == (mh, mw):
            ex_slice = excl_mask[my0:my1, mx0:mx1]
            lam = lam & ~((ex_slice[:, :, 0] > 0) if ex_slice.ndim == 3 else (ex_slice > 0))
    if int(np.count_nonzero(lam)) < 9:
        return nan_out

    eff_px_per_mm, _ = resolve_effective_px_per_mm(pixels_per_mm, sheet_bgr.shape[:2])
    dist = cv2.distanceTransform(lam.astype(np.uint8), cv2.DIST_L2, 3)
    max_d = float(dist.max())
    if max_d <= 0.0:
        return nan_out

    margin_px = min(float(np.clip(0.40 * eff_px_per_mm, 2.0, 12.0)), max(1.0, 0.22 * max_d))
    interior = dist >= margin_px
    if int(np.count_nonzero(interior)) < 16:
        interior = dist >= 1.0
    if int(np.count_nonzero(interior)) < 9:
        interior = lam

    if isinstance(keypoints_31, Mapping):
        mv_items = [
            (int(ks[8:]), (pt[0] - mx0, pt[1] - my0))
            for k, v in keypoints_31.items()
            for ks in (str(k),)
            if ks.startswith("midvein_") and ks[8:].isdigit()
            for pt in (_parse_kp_xy(v),)
            if pt is not None
        ]
        if len(mv_items) >= 2:
            mv_items.sort(key=lambda item: item[0])
            pts_mv = np.round(np.array([p for _, p in mv_items], dtype=np.float64)).astype(np.int32).reshape(-1, 1, 2)
            midrib_u8 = np.zeros(lam.shape, dtype=np.uint8)
            cv2.polylines(midrib_u8, [pts_mv], False, 255, thickness=max(1, int(round(0.18 * eff_px_per_mm))))
            lat_interior = interior & (midrib_u8 == 0)
            if int(np.count_nonzero(lat_interior)) >= 24:
                interior = lat_interior

    auto_l, paper_a, paper_b = estimate_sheet_paper_cie_l_ref_and_ab(sheet_bgr)
    eff_paper_l = auto_l
    if paper_l_ref is not None:
        try:
            p_val = float(paper_l_ref)
            if np.isfinite(p_val) and p_val >= 50.0:
                eff_paper_l = p_val
        except (TypeError, ValueError):
            pass

    lab = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2LAB)
    l_cie = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    a_c, b_c = lab[:, :, 1].astype(np.float32) - 128.0, lab[:, :, 2].astype(np.float32) - 128.0
    chroma = np.hypot(a_c, b_c)

    sigma_fine = float(np.clip(0.08 * eff_px_per_mm, 0.6, 4.0))
    sigma_coarse = float(max(2.5 * sigma_fine, np.clip(0.35 * eff_px_per_mm, 2.2, 16.0)))
    bp_map = cv2.GaussianBlur(l_cie, (0, 0), sigmaX=sigma_fine, sigmaY=sigma_fine) - cv2.GaussianBlur(l_cie, (0, 0), sigmaX=sigma_coarse, sigmaY=sigma_coarse)
    scale_norm = 90.0 / max(eff_paper_l, 50.0)

    # Wool-safe connected-component strap exclusion (retains fibrous white tomentum with DoG MAD > 1.8).
    if eff_paper_l >= 72.0:
        cand_strap = lam & (l_cie >= (eff_paper_l - 3.5)) & (chroma <= 6.0)
        if int(np.count_nonzero(cand_strap)) >= 8:
            rim_mask = lam & (dist <= max(2.0, margin_px + 1.5))
            n_cc, cc_labels, cc_stats, _ = cv2.connectedComponentsWithStats(cand_strap.astype(np.uint8), connectivity=8)
            strap_mask = np.zeros(lam.shape, dtype=bool)
            for lbl in range(1, n_cc):
                if int(cc_stats[lbl, cv2.CC_STAT_AREA]) < 8:
                    continue
                comp = cc_labels == lbl
                comp_dog_mad = _component_intrinsic_dog_mad(l_cie, comp, sigma_fine, sigma_coarse, scale_norm)
                if comp_dog_mad > 1.8:
                    continue
                matches_paper = abs(float(np.median(a_c[comp])) - paper_a) <= 2.2 and abs(float(np.median(b_c[comp])) - paper_b) <= 2.8 and comp_dog_mad <= 1.4
                is_cross = False
                if not matches_paper:
                    ys_c, xs_c = np.where(comp)
                    (_, (rw, rh), _) = cv2.minAreaRect(np.column_stack([xs_c, ys_c]).astype(np.float32))
                    long_d, short_d = max(float(rw), float(rh)), max(min(float(rw), float(rh)), 1e-3)
                    if (long_d / short_d) >= 2.8:
                        rim_hit = comp & rim_mask
                        n_rim, _, _, _ = cv2.connectedComponentsWithStats(rim_hit.astype(np.uint8), connectivity=8)
                        if (n_rim - 1) >= 2:
                            is_cross = True
                        elif np.any(rim_hit):
                            ry, rx = np.where(rim_hit)
                            is_cross = math.hypot(float(rx.max() - rx.min()), float(ry.max() - ry.min())) >= 0.5 * long_d
                if matches_paper or is_cross:
                    strap_mask |= comp
            non_strap = interior & (~strap_mask)
            if int(np.count_nonzero(non_strap)) >= 16:
                interior = non_strap
    if int(np.count_nonzero(interior)) < 9:
        return nan_out

    l_int = l_cie[interior]
    med_l, q85_l = float(np.median(l_int)), float(np.percentile(l_int, 85.0))
    med_a, med_b = float(np.median(a_c[interior])), float(np.median(b_c[interior]))
    y_ratio = max(0.0, _cie_l_to_rel_y(med_l) / max(_cie_l_to_rel_y(eff_paper_l), 1e-4))
    l_norm = float(np.clip(1.16 * (y_ratio ** (1.0 / 3.0)) - 0.16 if y_ratio > (6.0 / 29.0) ** 3 else 9.03296296 * y_ratio, 0.0, 1.25))
    bp_vals = bp_map[interior]
    mad_sd = 1.4826 * float(np.median(np.abs(bp_vals - float(np.median(bp_vals))))) * scale_norm

    iy, ix = np.where(interior)
    base_pt = (_parse_kp_xy(keypoints_31.get("lamina_base")) or _parse_kp_xy(keypoints_31.get("base_center"))) if isinstance(keypoints_31, Mapping) else None
    if base_pt is not None:
        d_base = np.hypot(ix.astype(np.float32) - float(base_pt[0] - mx0), iy.astype(np.float32) - float(base_pt[1] - my0))
        basal_sel = d_base <= float(np.percentile(d_base, 20.0))
    else:
        basal_sel = iy >= float(np.percentile(iy.astype(np.float32), 80.0))
    basal_delta = float(np.median(l_cie[iy[basal_sel], ix[basal_sel]])) - eff_paper_l if int(np.count_nonzero(basal_sel)) >= 8 else (med_l - eff_paper_l)

    return {
        "blade_cie_l_norm": round(l_norm, 4),
        "blade_cie_l_sd": round(max(0.0, mad_sd), 4),
        "blade_cie_l_sheet_delta": round(med_l - eff_paper_l, 4),
        "blade_cie_l_q85_sheet_delta": round(q85_l - eff_paper_l, 4),
        "basal_zone_cie_l_delta": round(basal_delta, 4),
        "blade_cie_a_delta": round(med_a - paper_a, 4),
        "blade_cie_b_delta": round(med_b - paper_b, 4),
    }


def _extract_item_bbox(item: Dict[str, Any]) -> Optional[Tuple[float, float, float, float]]:
    bb = item.get("bbox")
    return (float(bb[0]), float(bb[1]), float(bb[2]), float(bb[3])) if (bb is not None and len(bb) >= 4) else None


def estimate_corymb_cluster_head_count(
    cluster_items: List[Dict[str, Any]],
    single_head_items: List[Dict[str, Any]],
    pixels_per_mm: Optional[float] = None,
    sheet_shape_hw: Tuple[int, int] = (0, 0),
) -> Tuple[int, int]:
    """Estimate additional capitula inside `Seed_Fruit_MANY` (3) and `Flower_MANY` (5) corymb clusters."""
    clusters = [it for it in (cluster_items or []) if isinstance(it, dict) and it.get("class_id") is not None and int(it["class_id"]) in (3, 5)]
    if not clusters:
        return (0, 0)
    singles = [it for it in (single_head_items or []) if isinstance(it, dict) and it.get("class_id") is not None and int(it["class_id"]) in (2, 4, 6)]
    eff_px_per_mm, _ = resolve_effective_px_per_mm(pixels_per_mm, sheet_shape_hw)

    single_areas, single_centers = [], []
    for item in singles:
        s_mask, s_bb, area_val = item.get("mask"), _extract_item_bbox(item), 0.0
        if isinstance(s_mask, np.ndarray) and s_mask.size > 0:
            area_val = float(np.count_nonzero(s_mask > 0))
        if area_val <= 0.0 and s_bb is not None:
            area_val = 0.55 * max(0.0, s_bb[2] - s_bb[0]) * max(0.0, s_bb[3] - s_bb[1])
        if area_val > 0.0:
            single_areas.append(area_val)
        cxy = item.get("center_xy")
        if cxy is not None and len(cxy) >= 2:
            single_centers.append((float(cxy[0]), float(cxy[1])))
        elif s_bb is not None:
            single_centers.append((0.5 * (s_bb[1] + s_bb[3]), 0.5 * (s_bb[0] + s_bb[2])))

    ref_head_area = float(np.median(single_areas)) if single_areas else (28.0 * (eff_px_per_mm ** 2))
    extra_heads = 0
    for cl in clusters:
        cl_bb = _extract_item_bbox(cl)
        if cl_bb is not None:
            ymin, xmin, ymax, xmax = cl_bb
            bbox_area = max(0.0, ymax - ymin) * max(0.0, xmax - xmin)
            n_overlap = sum(1 for cx, cy in single_centers if xmin <= cx <= xmax and ymin <= cy <= ymax)
        else:
            bbox_area, n_overlap = 0.0, 0
        cl_mask, cl_crop, fg_area = cl.get("mask"), cl.get("crop_bgr"), None
        if isinstance(cl_mask, np.ndarray) and cl_mask.size > 0:
            nz = int(np.count_nonzero(cl_mask > 0))
            if nz > 0:
                fg_area = float(nz)
        if fg_area is None and isinstance(cl_crop, np.ndarray) and cl_crop.ndim == 3 and cl_crop.shape[2] >= 3 and cl_crop.size > 0:
            crop3 = cl_crop[:, :, :3]
            l_star = cv2.cvtColor(crop3, cv2.COLOR_BGR2LAB)[:, :, 0].astype(np.float32) * (100.0 / 255.0)
            s_norm = cv2.cvtColor(crop3, cv2.COLOR_BGR2HSV)[:, :, 1].astype(np.float32) / 255.0
            fg_area = float(np.count_nonzero((l_star < 78.0) | (s_norm > 0.16)))
        if fg_area is None:
            fg_area = 0.55 * bbox_area
        extra_heads += max(0, int(np.clip(round(fg_area / max(ref_head_area, 16.0)), 2, 35)) - n_overlap)
    return (int(extra_heads), int(len(clusters)))


__all__ = [
    "resolve_effective_px_per_mm",
    "estimate_sheet_paper_cie_l_ref_and_ab",
    "estimate_sheet_paper_cie_l_ref",
    "compute_calibrated_blade_photometry",
    "estimate_corymb_cluster_head_count",
]
