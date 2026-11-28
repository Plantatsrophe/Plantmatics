#!/usr/bin/env python3
"""
Step 02 lateral capitulum, phyllary-tip, calyculus, and subtending-peduncle
optical micro-character extraction (`scripts/vision/step02_capitulum_microchars.py`).

Operates on native-DPI single-capitulum crops (`Seed_Fruit_ONE = 2`,
`Flower_ONE = 4`, `Bud = 6`) and subtending full-sheet peduncle corridors.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

from scripts.vision.step02_photometry_corymb_utils import (
    _component_intrinsic_dog_mad,
    estimate_sheet_paper_cie_l_ref_and_ab,
    resolve_effective_px_per_mm,
)

MICROCHAR_REPRO_KEYS: Tuple[str, ...] = (
    "phyllary_tip_darkening_delta",
    "phyllary_tip_anthocyanin_delta",
    "involucre_basal_cie_l_delta",
    "peduncle_cie_l_sheet_delta",
    "calyculus_flare_index",
)


def _nan_microchars() -> Dict[str, float]:
    return {k: np.nan for k in MICROCHAR_REPRO_KEYS}


def _strip_smooth_paper_straps(
    fg_mask: np.ndarray,
    l_cie: np.ndarray,
    a_c: np.ndarray,
    b_c: np.ndarray,
    eff_paper_l: float,
    paper_a: float,
    paper_b: float,
    sigma_fine: float,
    sigma_coarse: float,
    scale_norm: float,
) -> np.ndarray:
    """Remove smooth white paper mounting straps while retaining fibrous white wool (`DoG MAD > 1.8`)."""
    if eff_paper_l < 70.0 or int(np.count_nonzero(fg_mask)) < 12:
        return fg_mask
    chroma = np.hypot(a_c, b_c)
    cand_strap = fg_mask & (l_cie >= (eff_paper_l - 4.0)) & (chroma <= 7.0)
    if int(np.count_nonzero(cand_strap)) < 8:
        return fg_mask
    n_cc, cc_labels, cc_stats, _ = cv2.connectedComponentsWithStats(cand_strap.astype(np.uint8), connectivity=8)
    strap_mask = np.zeros(fg_mask.shape, dtype=bool)
    for lbl in range(1, n_cc):
        if int(cc_stats[lbl, cv2.CC_STAT_AREA]) < 6:
            continue
        comp = cc_labels == lbl
        dog_mad = _component_intrinsic_dog_mad(l_cie, comp, sigma_fine, sigma_coarse, scale_norm)
        if dog_mad > 1.8:
            continue
        matches_paper = (
            abs(float(np.median(a_c[comp])) - paper_a) <= 2.5
            and abs(float(np.median(b_c[comp])) - paper_b) <= 3.2
            and dog_mad <= 1.45
        )
        if matches_paper:
            strap_mask |= comp
    cleaned = fg_mask & (~strap_mask)
    return cleaned if int(np.count_nonzero(cleaned)) >= 8 else fg_mask


def extract_single_capitulum_microchars(
    crop_bgr: Optional[np.ndarray],
    sheet_bgr: Optional[np.ndarray],
    bbox: Tuple[int, int, int, int],
    mask: Optional[np.ndarray] = None,
    class_id: Optional[int] = None,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
    pixels_per_mm: Optional[float] = None,
    is_calibrated_mm: bool = False,
    paper_lab_ref: Optional[Tuple[float, float, float]] = None,
) -> Dict[str, float]:
    """
    Extract lateral involucre tip pigmentation, basal wool, subtending peduncle wool,
    and calyculus flare metrics for single capitula (`class_id` in `{2, 4, 6}`).
    """
    out = _nan_microchars()
    if class_id is not None and int(class_id) not in (2, 4, 6):
        return out
    if (
        crop_bgr is None
        or not isinstance(crop_bgr, np.ndarray)
        or crop_bgr.ndim != 3
        or crop_bgr.shape[2] < 3
        or crop_bgr.shape[0] < 10
        or crop_bgr.shape[1] < 10
    ):
        return out

    ch, cw = crop_bgr.shape[:2]
    ref_shape = sheet_bgr.shape[:2] if (isinstance(sheet_bgr, np.ndarray) and sheet_bgr.ndim >= 2) else (ch, cw)
    eff_px_per_mm, _ = resolve_effective_px_per_mm(pixels_per_mm if is_calibrated_mm else None, ref_shape)
    eff_paper_l, paper_a, paper_b = (
        paper_lab_ref
        if (paper_lab_ref is not None and len(paper_lab_ref) >= 3)
        else estimate_sheet_paper_cie_l_ref_and_ab(sheet_bgr if isinstance(sheet_bgr, np.ndarray) else crop_bgr)
    )

    sigma_fine = float(np.clip(0.08 * eff_px_per_mm, 0.6, 3.5))
    sigma_coarse = float(max(2.5 * sigma_fine, np.clip(0.32 * eff_px_per_mm, 2.0, 12.0)))
    scale_norm = 90.0 / max(eff_paper_l, 50.0)

    crop3 = crop_bgr[:, :, :3]
    hsv = cv2.cvtColor(crop3, cv2.COLOR_BGR2HSV)
    lab = cv2.cvtColor(crop3, cv2.COLOR_BGR2LAB)
    H_ch, S_ch, V_ch = hsv[:, :, 0], hsv[:, :, 1], hsv[:, :, 2]
    l_cie = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    a_c = lab[:, :, 1].astype(np.float32) - 128.0
    b_c = lab[:, :, 2].astype(np.float32) - 128.0
    chroma = np.hypot(a_c, b_c)
    dog_abs = np.abs(
        cv2.GaussianBlur(l_cie, (0, 0), sigmaX=sigma_fine, sigmaY=sigma_fine)
        - cv2.GaussianBlur(l_cie, (0, 0), sigmaX=sigma_coarse, sigmaY=sigma_coarse)
    ) * scale_norm

    border_pts = np.concatenate([crop3[0, :, :], crop3[-1, :, :], crop3[:, 0, :], crop3[:, -1, :]], axis=0).astype(np.float32)
    bg_med = np.median(border_pts, axis=0)
    bg_dist = np.linalg.norm(crop3.astype(np.float32) - bg_med[None, None, :], axis=2)
    l_blur5 = cv2.blur(l_cie, (5, 5))
    local_l_sd = np.sqrt(np.maximum(0.0, cv2.blur(l_cie * l_cie, (5, 5)) - l_blur5 * l_blur5)) * scale_norm
    paper_mask = (bg_dist < 16.5) & (S_ch <= 38) & (dog_abs <= 1.65) & (local_l_sd <= 1.85) & (l_cie >= max(68.0, eff_paper_l - 11.0))
    fg_plant = ~paper_mask
    if mask is not None and isinstance(mask, np.ndarray) and mask.shape[:2] == (ch, cw) and np.any(mask > 0):
        fg_plant = fg_plant & (mask > 0)
    if int(np.count_nonzero(fg_plant)) < 24:
        return out

    axis_ux, axis_uy = 0.0, -1.0
    has_axis = False
    if peduncle_dir_xy is not None and len(peduncle_dir_xy) >= 2:
        px_d, py_d = float(peduncle_dir_xy[0]), float(peduncle_dir_xy[1])
        p_norm = math.hypot(px_d, py_d)
        if math.isfinite(p_norm) and p_norm > 1e-6:
            axis_ux, axis_uy = px_d / p_norm, py_d / p_norm
            has_axis = True

    floret_mask = fg_plant & (H_ch >= 10) & (H_ch <= 36) & (S_ch >= 76) & (V_ch >= 85)
    pappus_raw = fg_plant & (~floret_mask) & (V_ch >= 158) & (S_ch <= 48)
    fg_ys, fg_xs = np.where(fg_plant)
    prelim_distal = floret_mask | pappus_raw
    prelim_cup = fg_plant & (~prelim_distal) & (V_ch <= 198)
    if int(np.count_nonzero(prelim_distal)) >= 8 and int(np.count_nonzero(prelim_cup)) >= 12:
        cy_c, cx_c = np.where(prelim_cup)
        dy_d, dx_d = np.where(prelim_distal)
        vx, vy = float(np.mean(dx_d) - np.mean(cx_c)), float(np.mean(dy_d) - np.mean(cy_c))
        vn = math.hypot(vx, vy)
        if vn >= 2.5:
            ux_loc, uy_loc = vx / vn, vy / vn
            if not has_axis or (ux_loc * axis_ux + uy_loc * axis_uy) < 0.25:
                axis_ux, axis_uy = ux_loc, uy_loc
                has_axis = True

    yy_g, xx_g = np.ogrid[:ch, :cw]
    proj_grid = xx_g.astype(np.float32) * axis_ux + yy_g.astype(np.float32) * axis_uy
    proj_fg = proj_grid[fg_plant]
    proj_mid = 0.5 * (float(np.min(proj_fg)) + float(np.max(proj_fg)))
    distal_half = proj_grid >= proj_mid
    if class_id is None or int(class_id) in (2, 4):
        aged_ochre = fg_plant & distal_half & (H_ch >= 8) & (H_ch <= 36) & (S_ch >= 38) & (S_ch <= 75) & (V_ch >= 75)
        floret_mask = floret_mask | aged_ochre
    distal_pappus = fg_plant & distal_half & (~floret_mask) & (V_ch >= 158) & (S_ch <= 48)
    involucre_herb = fg_plant & (~floret_mask) & (~distal_pappus) & (V_ch <= 202)
    min_inv_px = max(18, int(0.04 * ch * cw))
    if int(np.count_nonzero(involucre_herb)) < min_inv_px:
        return out

    # Reject top-down radial disc views where florets sit concentrically inside involucre
    if int(np.count_nonzero(floret_mask)) >= 12:
        fl_y, fl_x = np.where(floret_mask)
        in_y, in_x = np.where(involucre_herb)
        fg_w = float(int(fg_xs.max()) - int(fg_xs.min()) + 1)
        fg_h = float(int(fg_ys.max()) - int(fg_ys.min()) + 1)
        diam = max(fg_w, fg_h, 1.0)
        if (
            (min(fg_w, fg_h) / diam) >= 0.72
            and math.hypot(float(np.mean(fl_x) - np.mean(in_x)), float(np.mean(fl_y) - np.mean(in_y))) / diam <= 0.18
        ):
            return out

    cup_total = _strip_smooth_paper_straps(
        fg_plant & (~floret_mask) & (~distal_pappus),
        l_cie, a_c, b_c, eff_paper_l, paper_a, paper_b, sigma_fine, sigma_coarse, scale_norm,
    )
    iy, ix = np.where(involucre_herb)
    ty, tx = np.where(cup_total)
    p_long = ix.astype(np.float32) * axis_ux + iy.astype(np.float32) * axis_uy
    p_orth = -ix.astype(np.float32) * axis_uy + iy.astype(np.float32) * axis_ux
    t_long = tx.astype(np.float32) * axis_ux + ty.astype(np.float32) * axis_uy
    t_orth = -tx.astype(np.float32) * axis_uy + ty.astype(np.float32) * axis_ux
    p0 = min(float(np.percentile(p_long, 2.0)), float(np.percentile(t_long, 2.0))) if t_long.size else float(np.percentile(p_long, 2.0))
    p1 = float(np.percentile(p_long, 98.0))
    h_span = p1 - p0
    w_span = float(np.percentile(t_orth if t_orth.size >= 8 else p_orth, 95.0) - np.percentile(t_orth if t_orth.size >= 8 else p_orth, 5.0))
    if h_span < 9.0 or w_span < 5.0:
        return out

    # 1. Phyllary tip darkening and anthocyanic reddening (eroded 1 px from background edge)
    eroded = cv2.erode(involucre_herb.astype(np.uint8), np.ones((3, 3), dtype=np.uint8)) > 0
    cup_core = eroded if int(np.count_nonzero(eroded)) >= 16 else involucre_herb
    cy_i, cx_i = np.where(cup_core)
    s_core = (cx_i.astype(np.float32) * axis_ux + cy_i.astype(np.float32) * axis_uy - p0) / h_span
    keel_mask = (s_core >= 0.30) & (s_core <= 0.70)
    tip_mask = (s_core >= 0.78) & (s_core <= 1.02)
    if int(np.count_nonzero(keel_mask)) >= 8 and int(np.count_nonzero(tip_mask)) >= 5:
        kl_vals = l_cie[cy_i[keel_mask], cx_i[keel_mask]]
        ka_vals = a_c[cy_i[keel_mask], cx_i[keel_mask]]
        kc_vals = chroma[cy_i[keel_mask], cx_i[keel_mask]]
        non_wool_k = ~((kl_vals >= (eff_paper_l - 8.5)) & (kc_vals <= 8.5))
        if int(np.count_nonzero(non_wool_k)) >= 6:
            kl_vals, ka_vals = kl_vals[non_wool_k], ka_vals[non_wool_k]
        tl_vals = l_cie[cy_i[tip_mask], cx_i[tip_mask]]
        ta_vals = a_c[cy_i[tip_mask], cx_i[tip_mask]]
        out["phyllary_tip_darkening_delta"] = round(float(np.median(kl_vals) - np.percentile(tl_vals, 25.0)), 4)
        out["phyllary_tip_anthocyanin_delta"] = round(float(np.percentile(ta_vals, 75.0) - np.median(ka_vals)), 4)

    # 2. Proximal involucre basal wool (retaining proximal flocculent white wool, stripping paper straps)
    s_grid = (proj_grid - p0) / h_span
    prox_clean = cup_total & (s_grid >= -0.05) & (s_grid <= 0.34)
    if int(np.count_nonzero(prox_clean)) >= 8:
        out["involucre_basal_cie_l_delta"] = round(float(np.percentile(l_cie[prox_clean], 85.0) - eff_paper_l), 4)

    # 3. Calyculus basal flare on sharp lateral heads (>= 13.5 px/mm)
    s_total = (t_long - p0) / h_span
    if eff_px_per_mm >= 13.5 and h_span >= 24.0 and w_span >= 12.0 and t_long.size >= 16:
        basal_widths: List[float] = []
        mid_widths: List[float] = []
        for b0, b1, target_l in (
            (0.06, 0.14, basal_widths), (0.14, 0.22, basal_widths), (0.22, 0.30, basal_widths),
            (0.38, 0.48, mid_widths), (0.48, 0.58, mid_widths), (0.58, 0.68, mid_widths),
        ):
            sel_b = (s_total >= b0) & (s_total < b1)
            if int(np.count_nonzero(sel_b)) >= 4:
                orth_b = t_orth[sel_b]
                target_l.append(float(np.percentile(orth_b, 95.0) - np.percentile(orth_b, 5.0)))
        if basal_widths and mid_widths:
            w_mid_med = float(np.median(mid_widths))
            if w_mid_med >= 6.0:
                flare = float(np.max(basal_widths)) / w_mid_med
                out["calyculus_flare_index"] = round(float(np.clip(flare, 0.35, 2.50)), 4)

    # 4. Subtending peduncle wool reflectance along -axis_u on full-sheet image
    if isinstance(sheet_bgr, np.ndarray) and sheet_bgr.ndim == 3 and sheet_bgr.shape[2] >= 3 and t_long.size >= 8:
        sh, sw = sheet_bgr.shape[:2]
        ymin, xmin, _ymax, _xmax = [int(v) for v in bbox]
        prox_pts = s_total <= 0.16
        if int(np.count_nonzero(prox_pts)) >= 3:
            orth_med = float(np.median(t_orth[prox_pts]))
            px_s = float(xmin) + (p0 * axis_ux - orth_med * axis_uy)
            py_s = float(ymin) + (p0 * axis_uy + orth_med * axis_ux)
            l_ped_max = float(np.clip(3.2 * eff_px_per_mm, 8.0, 1.6 * h_span))
            hw_ped = float(np.clip(0.30 * w_span, 3.0, 1.4 * eff_px_per_mm))
            end_x, end_y = px_s - axis_ux * l_ped_max, py_s - axis_uy * l_ped_max
            pad_r = int(math.ceil(hw_ped + 4.0))
            rx0, ry0 = max(0, int(math.floor(min(px_s, end_x))) - pad_r), max(0, int(math.floor(min(py_s, end_y))) - pad_r)
            rx1, ry1 = min(sw, int(math.ceil(max(px_s, end_x))) + pad_r + 1), min(sh, int(math.ceil(max(py_s, end_y))) + pad_r + 1)
            if (rx1 - rx0) >= 6 and (ry1 - ry0) >= 6:
                roi_bgr = sheet_bgr[ry0:ry1, rx0:rx1, :3]
                ry_g, rx_g = np.ogrid[ry0:ry1, rx0:rx1]
                dx_p, dy_p = rx_g.astype(np.float32) - px_s, ry_g.astype(np.float32) - py_s
                u_ped = -(dx_p * axis_ux + dy_p * axis_uy)
                v_ped = np.abs(-dx_p * axis_uy + dy_p * axis_ux)
                corridor = (u_ped >= 1.5) & (u_ped <= l_ped_max) & (v_ped <= hw_ped)
                if int(np.count_nonzero(corridor)) >= 12:
                    lab_p = cv2.cvtColor(roi_bgr, cv2.COLOR_BGR2LAB)
                    l_p = lab_p[:, :, 0].astype(np.float32) * (100.0 / 255.0)
                    a_p, b_p = lab_p[:, :, 1].astype(np.float32) - 128.0, lab_p[:, :, 2].astype(np.float32) - 128.0
                    lp_b5 = cv2.blur(l_p, (5, 5))
                    lp_sd = np.sqrt(np.maximum(0.0, cv2.blur(l_p * l_p, (5, 5)) - lp_b5 * lp_b5)) * scale_norm
                    dog_p = np.abs(cv2.GaussianBlur(l_p, (0, 0), sigmaX=sigma_fine, sigmaY=sigma_fine) - cv2.GaussianBlur(l_p, (0, 0), sigmaX=sigma_coarse, sigmaY=sigma_coarse)) * scale_norm
                    is_paper_p = (l_p >= (eff_paper_l - 4.5)) & (np.hypot(a_p, b_p) <= 7.5) & (dog_p <= 1.65) & (lp_sd <= 1.85)
                    fg_ped = _strip_smooth_paper_straps(corridor & (~is_paper_p), l_p, a_p, b_p, eff_paper_l, paper_a, paper_b, sigma_fine, sigma_coarse, scale_norm)
                    if int(np.count_nonzero(fg_ped)) >= 10:
                        out["peduncle_cie_l_sheet_delta"] = round(float(np.percentile(l_p[fg_ped], 80.0) - eff_paper_l), 4)

    return out


def summarize_capitulum_microchars(cap_items_sub: Sequence[Dict[str, Any]]) -> Dict[str, float]:
    """
    Aggregate head-level micro-characters to `(catalogNumber, plant_individual_id)`.
    Prioritizes pre-senescent `Flower_ONE (4)` and `Bud (6)` heads over `Seed_Fruit_ONE (2)`
    for phyllary-tip pigmentation so fruiting senescence never masks fresh tip color.
    """
    out = _nan_microchars()
    single_items = [
        it for it in (cap_items_sub or ())
        if isinstance(it, dict) and isinstance(it.get("microchars"), dict) and int(it.get("class_id", -1)) in (2, 4, 6)
    ]
    if not single_items:
        return out

    fresh_items = [it for it in single_items if int(it.get("class_id", -1)) in (4, 6)]

    def _med_from(pool: Sequence[Dict[str, Any]], key: str) -> float:
        vals = [
            float(it["microchars"][key])
            for it in pool
            if it["microchars"].get(key) is not None and np.isfinite(float(it["microchars"][key]))
        ]
        return round(float(np.median(vals)), 4) if vals else np.nan

    for tip_key in ("phyllary_tip_darkening_delta", "phyllary_tip_anthocyanin_delta"):
        v_fresh = _med_from(fresh_items, tip_key)
        out[tip_key] = v_fresh if np.isfinite(v_fresh) else _med_from(single_items, tip_key)

    for key in ("involucre_basal_cie_l_delta", "peduncle_cie_l_sheet_delta", "calyculus_flare_index"):
        out[key] = _med_from(single_items, key)

    return out


__all__ = [
    "MICROCHAR_REPRO_KEYS",
    "extract_single_capitulum_microchars",
    "summarize_capitulum_microchars",
]

