#!/usr/bin/env python3
"""
Step 02 late capitulum phenology biomarkers: oxidized aged-ochre/amber corolla
detection for historical herbarium sheets and peripheral radiate ray-ligule
discrimination separating radiate from discoid heads.
"""

from __future__ import annotations

from typing import Dict, Tuple

import cv2
import numpy as np


def extract_biomarker_aged_ochre_corollas(
    crop_rgb: np.ndarray,
    h_min_deg: float = 10.0,
    h_max_deg: float = 48.0,
    s_min: float = 0.18,
    s_max: float = 0.32,
    v_min: float = 0.28,
    min_b_star_excess: float = 8.0,
    ochre_threshold: float = 0.10,
) -> Tuple[bool, float, Dict[str, float]]:
    """
    Biomarker B2: Aged Ochre and Amber Corollas (Herbarium Anthesis Indicator).

    Detects oxidized ochre, amber, and browned ligules or disc corollas in the
    upper capitulum dome by thresholding HSV hue/saturation/value and CIE b*
    excess over border mounting paper while excluding bright achromatic pappus plumes.
    """
    if (
        crop_rgb is None
        or not isinstance(crop_rgb, np.ndarray)
        or crop_rgb.size == 0
        or crop_rgb.ndim < 3
        or crop_rgb.shape[2] < 3
    ):
        return False, 0.0, {
            "ochre_ratio": 0.0,
            "ochre_pixel_count": 0,
            "paper_b_star": 0.0,
        }

    h, w = crop_rgb.shape[:2]
    if h < 3 or w < 3:
        return False, 0.0, {
            "ochre_ratio": 0.0,
            "ochre_pixel_count": 0,
            "paper_b_star": 0.0,
        }

    rgb3 = crop_rgb[:, :, :3]
    crop_lab = cv2.cvtColor(rgb3, cv2.COLOR_RGB2LAB)
    crop_hsv = cv2.cvtColor(rgb3, cv2.COLOR_RGB2HSV)

    l_full = crop_lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    b_full = crop_lab[:, :, 2].astype(np.float32) - 128.0
    h_full = crop_hsv[:, :, 0].astype(np.float32) * 2.0
    s_full = crop_hsv[:, :, 1].astype(np.float32) / 255.0
    v_full = crop_hsv[:, :, 2].astype(np.float32) / 255.0

    # Estimate border mounting paper CIE b* baseline from outermost 3 px border
    border_mask = np.zeros((h, w), dtype=bool)
    bw = min(3, max(1, h // 4), max(1, w // 4))
    border_mask[:bw, :] = True
    border_mask[-bw:, :] = True
    border_mask[:, :bw] = True
    border_mask[:, -bw:] = True
    paper_border = border_mask & (l_full >= 65.0) & (s_full <= 0.22)
    b_paper = float(np.median(b_full[paper_border])) if np.any(paper_border) else 0.0

    upper_h = max(1, int(round(h * 0.58)))
    l_star = l_full[:upper_h, :]
    b_star = b_full[:upper_h, :]
    h_deg = h_full[:upper_h, :]
    s_norm = s_full[:upper_h, :]
    v_norm = v_full[:upper_h, :]

    is_achromatic_pappus = (l_star >= 75.0) & (s_norm <= 0.18)
    ochre_mask = (
        (h_deg >= h_min_deg)
        & (h_deg <= h_max_deg)
        & (s_norm >= s_min)
        & (s_norm < max(s_max, s_min + 1e-4))
        & (v_norm >= v_min)
        & ((b_star - b_paper) >= min_b_star_excess)
        & (~is_achromatic_pappus)
    )

    ochre_pixels = int(np.count_nonzero(ochre_mask))
    n_upper = max(1, upper_h * w)
    ochre_ratio = float(ochre_pixels / n_upper)
    has_aged_ochre = bool(ochre_ratio >= ochre_threshold)

    diagnostics = {
        "ochre_ratio": round(ochre_ratio, 4),
        "ochre_pixel_count": ochre_pixels,
        "paper_b_star": round(b_paper, 2),
    }

    return has_aged_ochre, float(round(ochre_ratio, 4)), diagnostics


def extract_biomarker_radiate_ligules(
    oriented_rgb: np.ndarray,
    min_floret_yellow_saturation: float = 0.32,
) -> Tuple[bool, Dict[str, float]]:
    """
    Detects peripheral yellow or aged-ochre ray ligules extending outside the
    compact central disc core / involucre cup, separating radiate from discoid heads.
    """
    if (
        oriented_rgb is None
        or not isinstance(oriented_rgb, np.ndarray)
        or oriented_rgb.size == 0
        or oriented_rgb.ndim < 3
        or oriented_rgb.shape[2] < 3
    ):
        return False, {"outer_corolla_frac": 0.0, "corolla_to_cup_width_ratio": 0.0}

    h, w = oriented_rgb.shape[:2]
    if h < 8 or w < 8:
        return False, {"outer_corolla_frac": 0.0, "corolla_to_cup_width_ratio": 0.0}

    rgb3 = oriented_rgb[:, :, :3]
    lab = cv2.cvtColor(rgb3, cv2.COLOR_RGB2LAB)
    hsv = cv2.cvtColor(rgb3, cv2.COLOR_RGB2HSV)

    l_star = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    b_star = lab[:, :, 2].astype(np.float32) - 128.0
    h_deg = hsv[:, :, 0].astype(np.float32) * 2.0
    s_norm = hsv[:, :, 1].astype(np.float32) / 255.0
    v_norm = hsv[:, :, 2].astype(np.float32) / 255.0

    border_mask = np.zeros((h, w), dtype=bool)
    bw = min(3, max(1, h // 4), max(1, w // 4))
    border_mask[:bw, :] = True
    border_mask[-bw:, :] = True
    border_mask[:, :bw] = True
    border_mask[:, -bw:] = True
    paper_border = border_mask & (l_star >= 65.0) & (s_norm <= 0.22)
    b_paper = float(np.median(b_star[paper_border])) if np.any(paper_border) else 0.0

    is_paper = (l_star >= 75.0) & (s_norm <= 0.20)
    is_pappus = (l_star >= 75.0) & (s_norm <= 0.18)
    non_paper = ~is_paper
    if int(np.count_nonzero(non_paper)) < 24:
        return False, {"outer_corolla_frac": 0.0, "corolla_to_cup_width_ratio": 0.0}

    fresh_yellow = (
        (h_deg >= 18.0)
        & (h_deg <= 52.0)
        & (s_norm >= max(0.28, float(min_floret_yellow_saturation) - 0.02))
        & (v_norm >= 0.35)
    )
    aged_ochre = (
        (h_deg >= 10.0)
        & (h_deg <= 48.0)
        & (s_norm >= 0.18)
        & (s_norm < max(float(min_floret_yellow_saturation), 0.19))
        & (v_norm >= 0.28)
        & ((b_star - b_paper) >= 8.0)
        & (~is_pappus)
    )
    corolla_mask = non_paper & (fresh_yellow | aged_ochre)
    n_cor = int(np.count_nonzero(corolla_mask))
    if n_cor < 16:
        return False, {"outer_corolla_frac": 0.0, "corolla_to_cup_width_ratio": 0.0}

    fg_ys, fg_xs = np.where(non_paper)
    fg_w = max(1.0, float(int(fg_xs.max()) - int(fg_xs.min()) + 1))
    fg_h = max(1.0, float(int(fg_ys.max()) - int(fg_ys.min()) + 1))
    fg_cx = 0.5 * float(int(fg_xs.min()) + int(fg_xs.max()))
    fg_cy = 0.5 * float(int(fg_ys.min()) + int(fg_ys.max()))

    cor_ys, cor_xs = np.where(corolla_mask)
    cor_cy = float(np.mean(cor_ys))
    inv_tissue = non_paper & (~corolla_mask) & (~is_pappus)
    lower_y = int(round(h * 0.48))
    inv_lower = inv_tissue[lower_y:, :]

    yy_grid, xx_grid = np.ogrid[:h, :w]
    has_lateral_rays = False
    has_radial_rays = False
    outer_cor_frac = 0.0
    width_ratio = 0.0

    is_lateral_cup = False
    if int(np.count_nonzero(inv_lower)) >= 12:
        il_y, il_x = np.where(inv_lower)
        cup_x0 = float(np.percentile(il_x, 5))
        cup_x1 = float(np.percentile(il_x, 95))
        cup_w = max(1.0, cup_x1 - cup_x0)
        cup_cx = 0.5 * (cup_x0 + cup_x1)
        cup_cy = float(lower_y) + float(np.mean(il_y))
        if cup_w >= 4.0 and (cup_cy - cor_cy) >= 0.12 * fg_h:
            is_lateral_cup = True
            cor_x0 = float(np.percentile(cor_xs, 2))
            cor_x1 = float(np.percentile(cor_xs, 98))
            cor_w = max(1.0, cor_x1 - cor_x0)
            width_ratio = cor_w / cup_w
            flank_r = max(0.60 * cup_w, min(0.30 * float(w), 0.28 * fg_w))
            flank_px = int(np.count_nonzero(corolla_mask & (np.abs(xx_grid.astype(np.float32) - cup_cx) >= flank_r)))
            outer_cor_frac = float(flank_px) / float(n_cor)
            has_lateral_rays = bool(width_ratio >= 1.45 and flank_px >= 12 and outer_cor_frac >= 0.14)

    if not is_lateral_cup:
        r_map = np.hypot(xx_grid.astype(np.float32) - fg_cx, yy_grid.astype(np.float32) - fg_cy)
        r_head = max(1.0, float(np.percentile(r_map[non_paper], 98)))
        outer_cor_px = int(np.count_nonzero(corolla_mask & (r_map >= 0.55 * r_head)))
        outer_cor_frac = float(outer_cor_px) / float(n_cor)
        all_cor_pts = np.column_stack([cor_xs, cor_ys]).astype(np.int32)
        hull_area = float(cv2.contourArea(cv2.convexHull(all_cor_pts))) if len(all_cor_pts) >= 3 else float(n_cor)
        cor_solidity = float(n_cor) / max(hull_area, 1.0)
        inner_fg = non_paper & (r_map < 0.45 * r_head)
        inner_fg_n = int(np.count_nonzero(inner_fg))
        inner_non_cor_frac = (
            float(np.count_nonzero(inner_fg & (~corolla_mask))) / float(inner_fg_n)
            if inner_fg_n > 0
            else 0.0
        )
        has_radial_rays = bool(
            outer_cor_px >= 16
            and outer_cor_frac >= 0.20
            and (cor_solidity < 0.88 or inner_non_cor_frac >= 0.20)
        )

    return bool(has_lateral_rays or has_radial_rays), {
        "outer_corolla_frac": round(float(outer_cor_frac), 4),
        "corolla_to_cup_width_ratio": round(float(width_ratio), 4),
    }


__all__ = [
    "extract_biomarker_aged_ochre_corollas",
    "extract_biomarker_radiate_ligules",
]
