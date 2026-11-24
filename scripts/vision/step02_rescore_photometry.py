#!/usr/bin/env python3
"""
Step 02 rescore photometry: calibrated CIELAB blade photometry recomputation,
algebraic sheet-delta backfills, and within-rosette abaxial pressed-leaf flags
for `scripts/vision/rescore_step02_qc_and_zones.py`.
"""

from __future__ import annotations

import logging
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.lm3_compat import rotate_image
from scripts.vision.step02_leaf_traits_utils import evaluate_rosette_abaxial_leaf_flags
from scripts.vision.step02_photometry_corymb_utils import (
    compute_calibrated_blade_photometry,
    estimate_sheet_paper_cie_l_ref,
)
from scripts.vision.step02_rescore_workers import (
    _align_rotated_patch_to_mask,
    _resolve_crop_file,
    _resolve_hook,
    _resolve_mask_file,
    _resolve_voucher_image_file,
)

logger = logging.getLogger("RescoreStep02")
_ORIG_GLOBALS: Dict[str, Any] = {}

PHOTOMETRY_FLOAT_COLS: Tuple[str, ...] = (
    "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "blade_cie_a_delta", "blade_cie_b_delta", "abaxial_cie_l_sheet_delta",
)
PHOTOMETRY_UPDATE_KEYS: Tuple[str, ...] = (
    "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "blade_cie_a_delta", "blade_cie_b_delta",
)

PhotoRecord = Tuple[
    int, Any, Any, Any, Optional[float], bool, Optional[float], Optional[float],
    Optional[Tuple[int, int, int, int]], Optional[float],
]
PhotoOutput = Tuple[
    int, Optional[float], Optional[float], Optional[float], Optional[float],
    Optional[float], Optional[float], Optional[float],
]


def _photo_sheet_worker(item: Tuple[str, List[PhotoRecord], Path]) -> List[PhotoOutput]:
    """Compute calibrated CIELAB blade photometry (or algebraic sheet-delta backfill) for one voucher sheet."""
    cat_w, recs_w, eff_data_dir = item
    resolve_img = _resolve_hook("_resolve_voucher_image_file", _resolve_voucher_image_file, _ORIG_GLOBALS)
    resolve_mask = _resolve_hook("_resolve_mask_file", _resolve_mask_file, _ORIG_GLOBALS)
    resolve_crop = _resolve_hook("_resolve_crop_file", _resolve_crop_file, _ORIG_GLOBALS)
    align_patch = _resolve_hook("_align_rotated_patch_to_mask", _align_rotated_patch_to_mask, _ORIG_GLOBALS)
    rot_img = _resolve_hook("rotate_image", rotate_image, _ORIG_GLOBALS)
    est_paper = _resolve_hook("estimate_sheet_paper_cie_l_ref", estimate_sheet_paper_cie_l_ref, _ORIG_GLOBALS)
    comp_photo = _resolve_hook("compute_calibrated_blade_photometry", compute_calibrated_blade_photometry, _ORIG_GLOBALS)

    img_f = resolve_img(recs_w[0][1], cat_w, eff_data_dir)
    img_bgr = cv2.imread(str(img_f), cv2.IMREAD_COLOR) if img_f is not None else None
    sh_hw = img_bgr.shape[:2] if (img_bgr is not None and img_bgr.ndim == 3 and img_bgr.size > 0) else None
    paper_ref: Optional[float] = est_paper(img_bgr) if sh_hw is not None else None
    out_w: List[PhotoOutput] = []

    for pos_i, _, m_raw, c_raw, px_mm_i, is_up_i, ex_ln_i, ex_lsd_i, bbox_i, midrib_ang_i in recs_w:
        m_f = resolve_mask(m_raw, eff_data_dir)
        m_u8 = cv2.imread(str(m_f), cv2.IMREAD_GRAYSCALE) if m_f is not None else None
        if m_u8 is None or int(np.count_nonzero(m_u8)) == 0:
            out_w.append((pos_i, None, None, None, None, None, None, None))
            continue
        mh_mw = m_u8.shape[:2]
        c_f = resolve_crop(c_raw, eff_data_dir)
        c_bgr = cv2.imread(str(c_f), cv2.IMREAD_COLOR) if c_f is not None else None
        res_p: Optional[Dict[str, float]] = None
        try:
            if c_bgr is not None and c_bgr.ndim == 3 and c_bgr.size > 0:
                m_use = (
                    m_u8
                    if c_bgr.shape[:2] == mh_mw
                    else cv2.resize(m_u8, (c_bgr.shape[1], c_bgr.shape[0]), interpolation=cv2.INTER_NEAREST)
                )
                res_p = comp_photo(c_bgr, m_use, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i)
            elif img_bgr is not None and sh_hw is not None and bbox_i is not None:
                sh_h, sh_w = int(sh_hw[0]), int(sh_hw[1])
                bx1 = max(0, min(sh_w - 1, bbox_i[0]))
                by1 = max(0, min(sh_h - 1, bbox_i[1]))
                bx2 = max(bx1 + 1, min(sh_w, bbox_i[2]))
                by2 = max(by1 + 1, min(sh_h, bbox_i[3]))
                crop_b = img_bgr[by1:by2, bx1:bx2]
                if crop_b.size > 0:
                    if is_up_i and midrib_ang_i is not None and math.isfinite(midrib_ang_i):
                        ang_cw = (-(midrib_ang_i + 90.0)) % 360.0
                        if 0.25 < ang_cw < 359.75:
                            crop_b = rot_img(crop_b, ang_cw, bg=255, nearest=False)
                    aligned_bgr = align_patch(crop_b, mh_mw, pad_val=255)
                    res_p = comp_photo(aligned_bgr, m_u8, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i)
            elif img_bgr is not None and sh_hw is not None:
                if mh_mw == sh_hw:
                    res_p = comp_photo(img_bgr, m_u8, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i)
                else:
                    ar_m = float(mh_mw[1]) / float(max(mh_mw[0], 1))
                    ar_s = float(sh_hw[1]) / float(max(sh_hw[0], 1))
                    is_scaled_sheet_mask = (not is_up_i) or (
                        abs(ar_m - ar_s) <= 0.08 * max(ar_s, 1e-6)
                        and min(mh_mw) >= 0.25 * min(sh_hw)
                    )
                    if is_scaled_sheet_mask:
                        m_res = cv2.resize(m_u8, (sh_hw[1], sh_hw[0]), interpolation=cv2.INTER_NEAREST)
                        res_p = comp_photo(img_bgr, m_res, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i)
                    elif ex_lsd_i is None and ex_ln_i is not None and 0.0 <= ex_ln_i <= 1.25 and paper_ref is not None:
                        y_p = ((paper_ref + 16.0) / 116.0) ** 3 if paper_ref > 8.0 else paper_ref / 903.3
                        y_rat = ((ex_ln_i + 0.16) / 1.16) ** 3 if ex_ln_i > 0.08 else (ex_ln_i / 9.03296296)
                        y_bld = max(0.0, y_rat * max(y_p, 1e-4))
                        med_l_rec = (
                            116.0 * (y_bld ** (1.0 / 3.0)) - 16.0
                            if y_bld > (6.0 / 29.0) ** 3
                            else 903.3 * y_bld
                        )
                        lsd_rec = round(float(med_l_rec - paper_ref), 4)
                        out_w.append((pos_i, None, None, lsd_rec if math.isfinite(lsd_rec) else None, None, None, None, None))
                        continue
        except Exception:
            res_p = None
        if not isinstance(res_p, dict):
            out_w.append((pos_i, None, None, None, None, None, None, None))
            continue

        def _fin_opt(k_name: str) -> Optional[float]:
            raw_v = res_p.get(k_name)
            if raw_v is None:
                return None
            fv = float(raw_v)
            return fv if math.isfinite(fv) else None

        out_w.append((
            pos_i,
            _fin_opt("blade_cie_l_norm"),
            _fin_opt("blade_cie_l_sd"),
            _fin_opt("blade_cie_l_sheet_delta"),
            _fin_opt("blade_cie_l_q85_sheet_delta"),
            _fin_opt("basal_zone_cie_l_delta"),
            _fin_opt("blade_cie_a_delta"),
            _fin_opt("blade_cie_b_delta"),
        ))
    return out_w


def recompute_leaf_photometry_and_abaxial_flags(
    el_df: pd.DataFrame,
    eff_data_dir: Path,
    recompute_photometry: bool,
    n_workers: int,
) -> Tuple[Dict[str, List[Optional[float]]], int, int]:
    """Recompute calibrated CIELAB blade photometry and evaluate rosette abaxial leaf flags in place."""
    log = _resolve_hook("logger", logger, _ORIG_GLOBALS)
    thread_pool = _resolve_hook("ThreadPoolExecutor", ThreadPoolExecutor, _ORIG_GLOBALS)
    eval_abaxial = _resolve_hook("evaluate_rosette_abaxial_leaf_flags", evaluate_rosette_abaxial_leaf_flags, _ORIG_GLOBALS)
    photo_worker = _resolve_hook("_photo_sheet_worker", _photo_sheet_worker, _ORIG_GLOBALS)

    for _pcol in PHOTOMETRY_FLOAT_COLS:
        if _pcol not in el_df.columns:
            el_df[_pcol] = np.nan
    if "is_abaxial_pressed_leaf" not in el_df.columns:
        el_df["is_abaxial_pressed_leaf"] = 0

    n_rows = len(el_df)
    upd_photo: Dict[str, List[Optional[float]]] = {k: [None] * n_rows for k in PHOTOMETRY_UPDATE_KEYS}
    n_photo_upd = 0
    n_sheet_delta_rec = 0
    prior_abax_mask = pd.to_numeric(el_df["is_abaxial_pressed_leaf"], errors="coerce").fillna(0).astype(int) == 1
    null_sheet_delta_mask = pd.to_numeric(el_df["blade_cie_l_sheet_delta"], errors="coerce").isna() & (~prior_abax_mask)

    if recompute_photometry or bool(null_sheet_delta_mask.any()):
        by_cat_photo: Dict[str, List[PhotoRecord]] = {}
        for pos_p, (_idx_p, row_p) in enumerate(el_df.iterrows()):
            if not recompute_photometry and not bool(null_sheet_delta_mask.iloc[pos_p]):
                continue
            cat_p = str(row_p["catalogNumber"]).strip()
            is_cal_p = int(row_p.get("is_calibrated_mm") or 0) == 1
            px_raw_p = row_p.get("pixels_per_mm")
            px_mm_p = (
                float(px_raw_p)
                if (
                    is_cal_p
                    and px_raw_p is not None
                    and pd.notna(px_raw_p)
                    and math.isfinite(float(px_raw_p))
                    and float(px_raw_p) > 0.0
                )
                else None
            )
            is_up_p = True if pd.isna(row_p.get("mask_is_upright")) else bool(int(row_p.get("mask_is_upright") or 0))
            ex_ln_raw = row_p.get("blade_cie_l_norm")
            ex_ln_p = (
                float(ex_ln_raw)
                if (ex_ln_raw is not None and pd.notna(ex_ln_raw) and math.isfinite(float(ex_ln_raw)))
                else None
            )
            ex_lsd_raw = row_p.get("blade_cie_l_sheet_delta")
            ex_lsd_p = (
                float(ex_lsd_raw)
                if (ex_lsd_raw is not None and pd.notna(ex_lsd_raw) and math.isfinite(float(ex_lsd_raw)))
                else None
            )
            bbox_xyxy_p: Optional[Tuple[int, int, int, int]] = None
            if all(c in row_p and pd.notna(row_p.get(c)) for c in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")):
                try:
                    bx1_p, by1_p = int(round(float(row_p["bbox_x1"]))), int(round(float(row_p["bbox_y1"])))
                    bx2_p, by2_p = int(round(float(row_p["bbox_x2"]))), int(round(float(row_p["bbox_y2"])))
                    if bx2_p > bx1_p and by2_p > by1_p:
                        bbox_xyxy_p = (bx1_p, by1_p, bx2_p, by2_p)
                except (TypeError, ValueError):
                    bbox_xyxy_p = None
            ma_raw = row_p.get("midrib_angle_deg")
            midrib_ang_p = (
                float(ma_raw)
                if (ma_raw is not None and pd.notna(ma_raw) and math.isfinite(float(ma_raw)))
                else None
            )
            by_cat_photo.setdefault(cat_p, []).append((
                pos_p,
                row_p.get("image_path"),
                row_p.get("mask_path"),
                row_p.get("crop_path"),
                px_mm_p,
                is_up_p,
                ex_ln_p,
                ex_lsd_p,
                bbox_xyxy_p,
                midrib_ang_p,
            ))

        photo_items = [(cat_w, recs_w, eff_data_dir) for cat_w, recs_w in by_cat_photo.items()]
        if len(photo_items) >= 4 and n_workers > 1:
            with thread_pool(max_workers=n_workers) as tpool:
                photo_batches = list(tpool.map(photo_worker, photo_items))
        else:
            photo_batches = [photo_worker(it) for it in photo_items]

        col_locs = {k: el_df.columns.get_loc(k) for k in PHOTOMETRY_UPDATE_KEYS}
        for batch in photo_batches:
            for pos_i, ln_v, ls_v, lsd_v, lq85_v, bas_l_v, a_v, b_v in batch:
                vals = (ln_v, ls_v, lsd_v, lq85_v, bas_l_v, a_v, b_v)
                did_upd = False
                for col_k, val_k in zip(PHOTOMETRY_UPDATE_KEYS, vals):
                    if val_k is not None:
                        upd_photo[col_k][pos_i] = val_k
                        el_df.iat[pos_i, col_locs[col_k]] = val_k
                        did_upd = True
                if did_upd:
                    if any(v is not None for v in (ln_v, ls_v, lq85_v, bas_l_v, a_v, b_v)):
                        n_photo_upd += 1
                    else:
                        n_sheet_delta_rec += 1
        log.info(
            f"CPU photometry recompute updated {n_photo_upd} / {n_rows} leaves where crop, bbox, or sheet-matched mask existed "
            f"(+ {n_sheet_delta_rec} algebraic sheet-delta backfills)."
        )

    # Evaluate within-rosette flipped abaxial basal leaves across rescored records
    abax_records: List[Dict[str, Any]] = []
    for pos_a, (_, r_a) in enumerate(el_df.iterrows()):
        lsd_a = r_a.get("blade_cie_l_sheet_delta")
        if (lsd_a is None or pd.isna(lsd_a)) and bool(prior_abax_mask.iloc[pos_a]) and pd.notna(r_a.get("abaxial_cie_l_sheet_delta")):
            lsd_a = r_a.get("abaxial_cie_l_sheet_delta")
        abax_records.append({
            "catalogNumber": str(r_a["catalogNumber"]).strip(),
            "plant_individual_id": int(r_a["plant_individual_id"]) if pd.notna(r_a.get("plant_individual_id")) else 0,
            "leaf_zone": str(r_a.get("leaf_zone") or "basal"),
            "blade_cie_l_norm": r_a.get("blade_cie_l_norm"),
            "blade_cie_l_sheet_delta": lsd_a,
            "blade_cie_a_delta": r_a.get("blade_cie_a_delta"),
            "blade_cie_b_delta": r_a.get("blade_cie_b_delta"),
        })
    eval_abaxial(abax_records)
    col_ln_idx = el_df.columns.get_loc("blade_cie_l_norm")
    col_lsd_idx = el_df.columns.get_loc("blade_cie_l_sheet_delta")
    el_df["is_abaxial_pressed_leaf"] = [1 if rec_a.get("is_abaxial_pressed_leaf") else 0 for rec_a in abax_records]
    el_df["abaxial_cie_l_sheet_delta"] = [rec_a.get("abaxial_cie_l_sheet_delta", np.nan) for rec_a in abax_records]
    for pos_a, rec_a in enumerate(abax_records):
        if rec_a.get("is_abaxial_pressed_leaf"):
            el_df.iat[pos_a, col_ln_idx] = np.nan
            el_df.iat[pos_a, col_lsd_idx] = np.nan
            upd_photo["blade_cie_l_norm"][pos_a] = None
            upd_photo["blade_cie_l_sheet_delta"][pos_a] = None

    return upd_photo, n_photo_upd, n_sheet_delta_rec


_ORIG_GLOBALS.update(globals())
