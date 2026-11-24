#!/usr/bin/env python3
"""
Step 02 rescore orchestration: per-voucher anatomical leaf-zone re-scoring,
atomic SQLite + CSV Tier A sync, and `rescore_tier_a_leaves_impl` for
`scripts/vision/rescore_step02_qc_and_zones.py`.
"""

from __future__ import annotations

import logging
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig, resolve_max_workers
from scripts.core.db import get_db_connection
from scripts.core.io_utils import atomic_write_text
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical
from scripts.vision.step02_rescore_photometry import recompute_leaf_photometry_and_abaxial_flags
from scripts.vision.step02_rescore_workers import (
    _resolve_hook,
    run_leaf_qc_evaluation,
    run_ruler_calibration_backfill,
)

logger = logging.getLogger("RescoreStep02")
_ORIG_GLOBALS: Dict[str, Any] = {}

_SYNC_COLS: Tuple[str, ...] = (
    "plant_individual_id", "solidity", "circularity", "detector_class_id",
    "is_crop_clipped", "hole_area_frac", "enclosed_void_frac", "bilateral_lobe_symmetry",
    "efa_qc_pass", "qc_flag_reason", "leaf_zone", "leaf_zone_confidence",
    "dist_from_crown_norm", "pixels_per_mm", "is_calibrated_mm", "scale_source",
    "area_mm2", "blade_length_mm", "blade_width_mm", "petiole_length_mm",
    "petiole_width_mm", "basal_cordate_index", "margin_serration_index",
    "margin_tooth_acuteness", "terminal_lobe_area_frac", "blade_cie_l_norm",
    "blade_cie_l_sd", "blade_cie_l_sheet_delta", "blade_cie_l_q85_sheet_delta",
    "basal_zone_cie_l_delta", "blade_cie_a_delta", "blade_cie_b_delta",
    "is_abaxial_pressed_leaf", "abaxial_cie_l_sheet_delta",
)
_STR_SYNC_COLS = frozenset({"qc_flag_reason", "leaf_zone", "scale_source"})
_BOOL_SYNC_COLS = frozenset({
    "is_crop_clipped", "bilateral_lobe_symmetry", "efa_qc_pass",
    "is_calibrated_mm", "is_abaxial_pressed_leaf",
})
_INT_SYNC_COLS = frozenset({"detector_class_id", "plant_individual_id"})
_NON_FLOAT_SYNC_COLS = _STR_SYNC_COLS | _BOOL_SYNC_COLS | _INT_SYNC_COLS

_EXTRACT_LEAVES_SQL = """
    SELECT el.*, pi.crown_x, pi.crown_y, pi.rosette_centroid_x, pi.rosette_centroid_y,
           pi.plant_height_px, pi.is_bent_stem, v.image_path,
           v.pixels_per_mm AS voucher_px_per_mm, v.is_calibrated_mm AS voucher_is_cal,
           COALESCE(vrp.capitula_count, vrp0.capitula_count, 0) AS eff_capitula_count
    FROM extracted_leaves el
    LEFT JOIN vouchers v ON v.catalogNumber = el.catalogNumber
    LEFT JOIN plant_individuals pi
        ON pi.catalogNumber = el.catalogNumber AND pi.plant_individual_id = el.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp
        ON vrp.catalogNumber = el.catalogNumber AND vrp.plant_individual_id = el.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp0
        ON vrp0.catalogNumber = el.catalogNumber AND vrp0.plant_individual_id = 0
    ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;
"""

_UPDATE_EXTRACTED_LEAVES_SQL = """
    UPDATE extracted_leaves
    SET solidity = COALESCE(?, solidity), circularity = COALESCE(?, circularity),
        detector_class_id = ?, is_crop_clipped = ?,
        hole_area_frac = ?, enclosed_void_frac = ?,
        bilateral_lobe_symmetry = ?, efa_qc_pass = ?,
        qc_flag_reason = ?, leaf_zone = ?,
        leaf_zone_confidence = COALESCE(?, leaf_zone_confidence),
        dist_from_crown_norm = COALESCE(?, dist_from_crown_norm),
        pixels_per_mm = COALESCE(?, pixels_per_mm), is_calibrated_mm = ?,
        scale_source = CASE WHEN ? != '' THEN ? ELSE scale_source END,
        area_mm2 = COALESCE(?, area_mm2), blade_length_mm = COALESCE(?, blade_length_mm),
        blade_width_mm = COALESCE(?, blade_width_mm), petiole_length_mm = COALESCE(?, petiole_length_mm),
        petiole_width_mm = COALESCE(?, petiole_width_mm), basal_cordate_index = COALESCE(?, basal_cordate_index),
        margin_serration_index = COALESCE(?, margin_serration_index),
        margin_tooth_acuteness = COALESCE(?, margin_tooth_acuteness),
        terminal_lobe_area_frac = COALESCE(?, terminal_lobe_area_frac),
        blade_cie_l_norm = CASE WHEN ? = 1 THEN NULL ELSE COALESCE(?, blade_cie_l_norm) END,
        blade_cie_l_sd = COALESCE(?, blade_cie_l_sd),
        blade_cie_l_sheet_delta = CASE WHEN ? = 1 THEN NULL ELSE COALESCE(?, blade_cie_l_sheet_delta) END,
        blade_cie_l_q85_sheet_delta = COALESCE(?, blade_cie_l_q85_sheet_delta),
        basal_zone_cie_l_delta = COALESCE(?, basal_zone_cie_l_delta),
        blade_cie_a_delta = COALESCE(?, blade_cie_a_delta),
        blade_cie_b_delta = COALESCE(?, blade_cie_b_delta),
        is_abaxial_pressed_leaf = ?, abaxial_cie_l_sheet_delta = ?
    WHERE catalogNumber = ? AND plant_individual_id = ? AND leaf_id = ?;
"""


def _opt_float(val: Any) -> Optional[float]:
    if val is None or pd.isna(val):
        return None
    fv = float(val)
    return fv if math.isfinite(fv) else None


def rescore_anatomical_leaf_zones(el_df: pd.DataFrame) -> List[Tuple[float, float, str, int]]:
    """Re-evaluate basal vs. cauline leaf zones and sibling half-fold screens per plant individual."""
    classify_zone = _resolve_hook("classify_leaf_zone_anatomical", classify_leaf_zone_anatomical, _ORIG_GLOBALS)
    crown_backfills: List[Tuple[float, float, str, int]] = []
    for (_cat, _pid), grp in el_df.groupby(["catalogNumber", "plant_individual_id"], sort=False):
        if len(grp) == 0:
            continue
        r0 = grp.iloc[0]
        cx_cr = float(r0["crown_x"]) if pd.notna(r0.get("crown_x")) else None
        cy_cr = float(r0["crown_y"]) if pd.notna(r0.get("crown_y")) else None
        if cx_cr is None or cy_cr is None:
            rx = float(r0["rosette_centroid_x"]) if pd.notna(r0.get("rosette_centroid_x")) else None
            ry = float(r0["rosette_centroid_y"]) if pd.notna(r0.get("rosette_centroid_y")) else None
            if rx is not None and ry is not None:
                cx_cr, cy_cr = rx, ry
                crown_backfills.append((cx_cr, cy_cr, str(_cat), int(_pid)))

        has_caps = int(r0.get("eff_capitula_count") or 0) > 0
        med_bld = float(grp["blade_length_px"].median()) if grp["blade_length_px"].notna().any() else 90.0
        p_height = (
            float(r0["plant_height_px"])
            if pd.notna(r0.get("plant_height_px")) and float(r0["plant_height_px"]) > 20.0
            else (max(320.0, 3.8 * max(med_bld, 40.0)) if len(grp) >= 2 else None)
        )
        cx_eval, cy_eval = (cx_cr if cx_cr is not None else 1500.0), (cy_cr if cy_cr is not None else 2200.0)

        # Reconstruct anatomical frame whenever >= 2 leaves exist (flowering OR sterile)
        if len(grp) >= 2 and p_height is not None:
            leaf_recs: List[Dict[str, Any]] = []
            n_prox_pet = 0
            for _, r_l in grp.iterrows():
                lid_i = int(r_l["leaf_id"])
                bld_px = float(r_l["blade_length_px"]) if pd.notna(r_l.get("blade_length_px")) and float(r_l["blade_length_px"]) > 0 else 80.0
                bld_w_px = float(r_l["blade_width_px"]) if pd.notna(r_l.get("blade_width_px")) and float(r_l["blade_width_px"]) > 0 else np.nan
                pet_px = float(r_l["petiole_length_px"]) if pd.notna(r_l.get("petiole_length_px")) and float(r_l["petiole_length_px"]) >= 0 else 0.0
                pet_w_px = float(r_l["petiole_width_px"]) if pd.notna(r_l.get("petiole_width_px")) and float(r_l["petiole_width_px"]) >= 0 else np.nan
                pbr_v = float(r_l["petiole_to_blade_ratio"]) if pd.notna(r_l.get("petiole_to_blade_ratio")) else (pet_px / max(bld_px, 1.0))
                dc_prev = float(r_l["dist_from_crown_norm"]) if pd.notna(r_l.get("dist_from_crown_norm")) else np.nan
                cy_rel = float(r_l["cy_rel_plant_bbox"]) if pd.notna(r_l.get("cy_rel_plant_bbox")) else np.nan
                s_est = max(0.0, min(1.0, 1.0 - cy_rel)) if np.isfinite(cy_rel) else (max(0.0, min(1.0, dc_prev)) if np.isfinite(dc_prev) else 0.0)
                y_c = cy_eval - s_est * p_height
                y_base, y_apex = min(cy_eval, y_c + 0.5 * bld_px), y_c - 0.5 * bld_px
                y_tip = min(cy_eval, y_base + min(pet_px, 0.45 * p_height))
                has_pet_m = bool(pet_px > 0.0 or pbr_v >= 0.12)
                if has_pet_m and s_est <= 0.45:
                    n_prox_pet += 1
                leaf_recs.append({
                    "leaf_id": lid_i, "plant_individual_id": int(_pid),
                    "p_base_sheet": (cx_eval, y_base), "p_apex_sheet": (cx_eval, y_apex), "p_petiole_tip": (cx_eval, y_tip),
                    "blade_length_px": bld_px, "blade_width_px": bld_w_px,
                    "petiole_length_px": pet_px, "petiole_width_px": pet_w_px,
                    "petiole_to_blade_ratio": pbr_v, "has_petiole_mask": has_pet_m,
                    "cy_rel_plant_bbox": cy_rel if np.isfinite(cy_rel) else (1.0 - s_est),
                    "bilateral_dice": float(r_l["bilateral_dice"]) if pd.notna(r_l.get("bilateral_dice")) else np.nan,
                    "si_a": float(r_l["si_a"]) if pd.notna(r_l.get("si_a")) else np.nan,
                    "lm3_gates_pass": bool(int(r_l["lm3_gates_pass"]) != 0) if pd.notna(r_l.get("lm3_gates_pass")) else True,
                })
            has_rosette_root = bool((cx_cr is not None and cy_cr is not None and n_prox_pet >= 1) or n_prox_pet >= 2)
            cap_box = [(cy_eval - p_height - 20.0, cx_eval - 20.0, cy_eval - p_height + 20.0, cx_eval + 20.0)] if has_caps else []
            root_box = [(cy_eval - 20.0, cx_eval - 20.0, cy_eval + 20.0, cx_eval + 20.0)] if has_rosette_root else []
            z_res = classify_zone(
                leaf_recs, crowns_by_pid={int(_pid): (cx_eval, cy_eval)}, capitula_detections=cap_box,
                root_boxes=root_box, is_bent_stem=bool(int(r0.get("is_bent_stem") or 0)),
            )
            for _, r_l in grp.iterrows():
                lid_i = int(r_l["leaf_id"])
                if lid_i in z_res:
                    mask_r = (el_df["catalogNumber"] == _cat) & (el_df["plant_individual_id"] == _pid) & (el_df["leaf_id"] == lid_i)
                    el_df.loc[mask_r, "leaf_zone"] = z_res[lid_i]["leaf_zone"]
                    el_df.loc[mask_r, "leaf_zone_confidence"] = z_res[lid_i]["leaf_zone_confidence"]
                    el_df.loc[mask_r, "dist_from_crown_norm"] = z_res[lid_i]["dist_from_crown_norm"]

        # Within-rosette sibling aspect-ratio half-fold screen on basal leaves
        grp_after = el_df.loc[grp.index]
        bas_grp = grp_after[grp_after["leaf_zone"].astype(str).str.lower() == "basal"]
        if len(bas_grp) >= 3:
            clean_b = bas_grp[bas_grp["efa_qc_pass"] == 1]
            ref_b = clean_b if len(clean_b) >= 2 else bas_grp
            med_ar, med_area = float(ref_b["aspect_ratio"].median()), float(ref_b["area_px"].median())
            for idx_b, rb in bas_grp.iterrows():
                ar_b = float(rb["aspect_ratio"]) if pd.notna(rb.get("aspect_ratio")) else med_ar
                area_b = float(rb["area_px"]) if pd.notna(rb.get("area_px")) else med_area
                rms_b = float(rb["_min_rms"])
                bd_b = float(rb["bilateral_dice"]) if pd.notna(rb.get("bilateral_dice")) else 1.0
                if (
                    ar_b > 1.85 * max(med_ar, 0.5)
                    and area_b < 0.65 * max(med_area, 1.0)
                    and ((np.isfinite(rms_b) and rms_b < 0.028) or bd_b < 0.78)
                ):
                    toks = [t for t in str(el_df.at[idx_b, "qc_flag_reason"]).split(";") if t and t != "PASS"]
                    if "LONGITUDINAL_MARGIN_CREASE" not in toks:
                        toks.append("LONGITUDINAL_MARGIN_CREASE")
                    el_df.at[idx_b, "efa_qc_pass"] = 0
                    el_df.at[idx_b, "qc_flag_reason"] = ";".join(toks)
    return crown_backfills


def persist_rescored_leaves_to_sqlite(
    conn: Any,
    el_df: pd.DataFrame,
    upd_photo: Dict[str, List[Optional[float]]],
    crown_backfills: List[Tuple[float, float, str, int]],
    ruler_cal_scales: Dict[str, float],
    unboxed_scales: Dict[str, float],
) -> None:
    """Persist rescored QC, traits, ruler scales, photometry, and leaf zones atomically to SQLite."""
    sql_rows: List[Tuple[Any, ...]] = []
    for pos_i, (_, r_u) in enumerate(el_df.iterrows()):
        ss_val, is_abax_i = str(r_u.get("scale_source") or ""), int(r_u.get("is_abaxial_pressed_leaf") or 0)
        sql_rows.append((
            _opt_float(r_u["solidity"]), _opt_float(r_u["circularity"]),
            int(r_u["detector_class_id"]), int(r_u["is_crop_clipped"]),
            float(r_u["hole_area_frac"]), float(r_u["enclosed_void_frac"]),
            int(r_u["bilateral_lobe_symmetry"]), int(r_u["efa_qc_pass"]),
            str(r_u["qc_flag_reason"]), str(r_u["leaf_zone"] or "basal"),
            _opt_float(r_u.get("leaf_zone_confidence")), _opt_float(r_u.get("dist_from_crown_norm")),
            _opt_float(r_u.get("pixels_per_mm")), int(r_u.get("is_calibrated_mm") or 0),
            ss_val, ss_val,
            _opt_float(r_u.get("area_mm2")), _opt_float(r_u.get("blade_length_mm")),
            _opt_float(r_u.get("blade_width_mm")), _opt_float(r_u.get("petiole_length_mm")),
            _opt_float(r_u.get("petiole_width_mm")), _opt_float(r_u.get("basal_cordate_index")),
            _opt_float(r_u.get("margin_serration_index")), _opt_float(r_u.get("margin_tooth_acuteness")),
            _opt_float(r_u.get("terminal_lobe_area_frac")),
            is_abax_i, _opt_float(upd_photo["blade_cie_l_norm"][pos_i]),
            _opt_float(upd_photo["blade_cie_l_sd"][pos_i]),
            is_abax_i, _opt_float(upd_photo["blade_cie_l_sheet_delta"][pos_i]),
            _opt_float(upd_photo["blade_cie_l_q85_sheet_delta"][pos_i]),
            _opt_float(upd_photo["basal_zone_cie_l_delta"][pos_i]),
            _opt_float(upd_photo["blade_cie_a_delta"][pos_i]),
            _opt_float(upd_photo["blade_cie_b_delta"][pos_i]),
            is_abax_i, _opt_float(r_u.get("abaxial_cie_l_sheet_delta")),
            str(r_u["catalogNumber"]), int(r_u["plant_individual_id"]), int(r_u["leaf_id"]),
        ))

    with conn:
        conn.executemany(_UPDATE_EXTRACTED_LEAVES_SQL, sql_rows)
        if crown_backfills:
            conn.executemany(
                "UPDATE plant_individuals SET crown_x = COALESCE(crown_x, ?), crown_y = COALESCE(crown_y, ?) "
                "WHERE catalogNumber = ? AND plant_individual_id = ?;",
                crown_backfills,
            )
        if ruler_cal_scales:
            v_cal_rows = [(float(px), str(cat)) for cat, px in ruler_cal_scales.items()]
            conn.executemany("UPDATE vouchers SET pixels_per_mm = ?, is_calibrated_mm = 1 WHERE catalogNumber = ?;", v_cal_rows)
        if unboxed_scales:
            v_rows = [(float(px), str(cat)) for cat, px in unboxed_scales.items() if cat not in ruler_cal_scales]
            if v_rows:
                conn.executemany("UPDATE vouchers SET pixels_per_mm = ?, is_calibrated_mm = 0 WHERE catalogNumber = ?;", v_rows)


def sync_rescored_csv_snapshots(
    el_df: pd.DataFrame,
    eff_manifest: Path,
    ruler_cal_scales: Dict[str, float],
    unboxed_scales: Dict[str, float],
) -> None:
    """Synchronize rescored voucher scale and leaf QC/zone columns into companion CSV snapshots."""
    log = _resolve_hook("logger", logger, _ORIG_GLOBALS)
    write_atomic = _resolve_hook("atomic_write_text", atomic_write_text, _ORIG_GLOBALS)

    if ruler_cal_scales or unboxed_scales:
        curated_csv = eff_manifest.parent / "curated_vouchers.csv"
        if curated_csv.is_file() and curated_csv.stat().st_size > 0:
            c_df = pd.read_csv(curated_csv, low_memory=False)
            if "catalogNumber" in c_df.columns:
                for ci in c_df.index:
                    c_key = str(c_df.at[ci, "catalogNumber"]).strip()
                    if c_key in ruler_cal_scales:
                        upd_v = (float(ruler_cal_scales[c_key]), True, "ruler_calibrated")
                    elif c_key in unboxed_scales:
                        upd_v = (float(unboxed_scales[c_key]), False, "tick_lattice_unboxed")
                    else:
                        continue
                    for col_v, val_v in zip(("pixels_per_mm", "is_calibrated_mm", "scale_source"), upd_v):
                        if col_v in c_df.columns:
                            c_df.at[ci, col_v] = val_v
                write_atomic(curated_csv, c_df.to_csv(index=False))

    if eff_manifest.is_file() and eff_manifest.stat().st_size > 0:
        m_df = pd.read_csv(eff_manifest, low_memory=False)
        if {"catalogNumber", "leaf_id"}.issubset(m_df.columns):
            if "plant_individual_id" not in m_df.columns:
                m_df["plant_individual_id"] = 0
            key_to_idx: Dict[Tuple[str, int, int], int] = {}
            cat_leaf_to_idx: Dict[Tuple[str, int], int] = {}
            for r_i, mr in m_df.iterrows():
                cat_m, lid_m = str(mr["catalogNumber"]).strip(), int(mr["leaf_id"])
                pid_m = int(mr["plant_individual_id"]) if pd.notna(mr["plant_individual_id"]) else 0
                key_to_idx[(cat_m, pid_m, lid_m)] = int(r_i)
                cat_leaf_to_idx.setdefault((cat_m, lid_m), int(r_i))
            for col_s in _SYNC_COLS:
                if col_s not in m_df.columns:
                    m_df[col_s] = (
                        "" if col_s in _STR_SYNC_COLS else
                        (col_s == "efa_qc_pass" if col_s in _BOOL_SYNC_COLS else (0 if col_s in _INT_SYNC_COLS else np.nan))
                    )
                elif col_s not in _NON_FLOAT_SYNC_COLS:
                    m_df[col_s] = pd.to_numeric(m_df[col_s], errors="coerce")
            for _, r_u in el_df.iterrows():
                cat_u, lid_u = str(r_u["catalogNumber"]).strip(), int(r_u["leaf_id"])
                mi = key_to_idx.get((cat_u, int(r_u["plant_individual_id"]), lid_u))
                if mi is None:
                    mi = cat_leaf_to_idx.get((cat_u, lid_u))
                if mi is not None:
                    for col_s in _SYNC_COLS:
                        if col_s in _BOOL_SYNC_COLS:
                            m_df.at[mi, col_s] = bool(r_u[col_s])
                        elif col_s in _INT_SYNC_COLS:
                            m_df.at[mi, col_s] = int(r_u[col_s])
                        elif col_s in _STR_SYNC_COLS:
                            m_df.at[mi, col_s] = str(r_u.get(col_s) or "")
                        elif col_s in r_u:
                            m_df.at[mi, col_s] = r_u[col_s]
            write_atomic(eff_manifest, m_df.to_csv(index=False))
            log.info(f"Synchronized updated QC and zone columns to {eff_manifest}")


def rescore_tier_a_leaves_impl(
    db_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    manifest_csv: Optional[Path] = None,
    segment_partial_leaves: Optional[bool] = None,
    recompute_photometry: bool = False,
    workers: Optional[int] = None,
) -> Dict[str, int]:
    """Re-evaluates read-only leaf QC gates, Packera traits, ruler scales, photometry, and basal/cauline zones."""
    proj_root = _resolve_hook("PROJECT_ROOT", PROJECT_ROOT, _ORIG_GLOBALS)
    log = _resolve_hook("logger", logger, _ORIG_GLOBALS)
    cfg_cls = _resolve_hook("PipelineConfig", PipelineConfig, _ORIG_GLOBALS)
    resolve_wk = _resolve_hook("resolve_max_workers", resolve_max_workers, _ORIG_GLOBALS)
    get_conn = _resolve_hook("get_db_connection", get_db_connection, _ORIG_GLOBALS)

    cfg = cfg_cls.from_yaml()
    eff_data_dir = (Path(data_dir) if data_dir is not None else (proj_root / "data")).resolve()
    eff_manifest = Path(manifest_csv).resolve() if manifest_csv is not None else (eff_data_dir / "tables" / "extracted_leaf_manifest.csv")
    seg_partial = bool(segment_partial_leaves) if segment_partial_leaves is not None else bool(getattr(cfg.segmentation, "segment_partial_leaves", False))

    conn = get_conn(db_path=db_path, init_schema=True)
    try:
        el_df = pd.read_sql_query(_EXTRACT_LEAVES_SQL, conn)
        if el_df.empty:
            log.info("No rows in extracted_leaves to re-score.")
            return {
                "total": 0, "efa_qc_pass": 0, "flagged": 0, "basal": 0, "cauline": 0,
                "ruler_calibrated_vouchers": 0, "tick_lattice_unboxed_vouchers": 0, "photometry_recomputed": 0,
            }

        n_rows = len(el_df)
        n_workers = max(1, int(workers)) if (workers is not None and int(workers) >= 1) else resolve_wk(default=os.cpu_count() or 8, cap=64)
        ruler_cal_scales, unboxed_scales = run_ruler_calibration_backfill(el_df, eff_data_dir, n_workers)
        run_leaf_qc_evaluation(el_df, eff_data_dir, seg_partial, n_workers)
        crown_backfills = rescore_anatomical_leaf_zones(el_df)
        upd_photo, n_photo_upd, n_sheet_delta_rec = recompute_leaf_photometry_and_abaxial_flags(
            el_df, eff_data_dir, recompute_photometry, n_workers
        )
        persist_rescored_leaves_to_sqlite(conn, el_df, upd_photo, crown_backfills, ruler_cal_scales, unboxed_scales)
        sync_rescored_csv_snapshots(el_df, eff_manifest, ruler_cal_scales, unboxed_scales)

        n_pass = int((el_df["efa_qc_pass"] == 1).sum())
        n_flag = int(n_rows - n_pass)
        n_bas = int((el_df["leaf_zone"].astype(str).str.lower() == "basal").sum())
        n_caul = int((el_df["leaf_zone"].astype(str).str.lower() == "cauline").sum())
        log.info(
            f"Re-scored {n_rows} leaves: efa_qc_pass={n_pass}, flagged={n_flag}, "
            f"basal={n_bas}, cauline={n_caul}, ruler_calibrated={len(ruler_cal_scales)}, "
            f"tick_lattice_unboxed={len(unboxed_scales)}."
        )
        return {
            "total": n_rows, "efa_qc_pass": n_pass, "flagged": n_flag,
            "basal": n_bas, "cauline": n_caul,
            "ruler_calibrated_vouchers": len(ruler_cal_scales),
            "tick_lattice_unboxed_vouchers": len(unboxed_scales),
            "photometry_recomputed": n_photo_upd,
            "sheet_delta_recovered": n_sheet_delta_rec,
        }
    finally:
        conn.close()


_ORIG_GLOBALS.update(globals())
