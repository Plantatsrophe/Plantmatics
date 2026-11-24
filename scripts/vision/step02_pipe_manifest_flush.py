#!/usr/bin/env python3
"""Step 02 atomic CSV manifest flushing, SQLite/CSV resumption, image prefetch, and skip helpers."""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

import cv2
import numpy as np
import pandas as pd

from scripts.core.db import (
    get_completed_qc_vouchers,
    get_db_connection,
    is_canonical_or_explicit_db_context,
    load_step02_records_from_db,
    upsert_vouchers_records,
)
from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES
from scripts.vision.step02_pipe_stage1_candidates import _resolve_scale_source

logger = logging.getLogger("SegmentAndExtract")

_MANIFEST_COLS: List[str] = [
    "catalogNumber", "plant_individual_id", "leaf_id",
    "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
    "solidity", "circularity", "extent", "aspect_ratio", "area_px", "area_mm2",
    "archetype_score", "bilateral_dice", "si_a", "lm3_gates_pass",
    "n_lamina_components", "detector_class_id", "is_crop_clipped",
    "hole_area_frac", "enclosed_void_frac", "bilateral_lobe_symmetry",
    "efa_qc_pass", "qc_flag_reason",
    "base_angle_deg", "apex_angle_deg", "lamina_curvature", "lamina_trace_length_px",
    "basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac",
    "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "blade_cie_a_delta", "blade_cie_b_delta",
    "is_abaxial_pressed_leaf", "abaxial_cie_l_sheet_delta",
    "scale_source", "pixels_per_mm", "is_calibrated_mm",
    "mask_path", "mask_is_upright", "contour_path", "ect_path", "leaf_zone",
]


def _flatten_pid_records(
    base_records: Dict[str, Dict[str, Any]], by_pid_records: Dict[str, Dict[int, Dict[str, Any]]]
) -> List[Dict[str, Any]]:
    """Expand per-voucher and per-plant_individual_id record maps into a flat row list."""
    rows: List[Dict[str, Any]] = []
    for c_k, base_rec in base_records.items():
        pid_map = by_pid_records.get(c_k)
        if pid_map and len(pid_map) > 1:
            for p_k in sorted(pid_map.keys()):
                rows.append(pid_map[p_k])
        elif pid_map and len(pid_map) == 1:
            rows.append(next(iter(pid_map.values())))
        else:
            rows.append(base_rec)
    return rows


def flush_pipeline_manifests(
    pipe: Any,
    all_extracted_records: List[Dict[str, Any]],
    all_failed_records: List[Dict[str, Any]],
    manifest_path: Path,
    failed_qc_path: Path,
    repro_manifest_path: Path,
    pheno_manifest_path: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build and atomically persist all 4 canonical Step 02 output tables."""
    extracted_df = pd.DataFrame(all_extracted_records)
    if extracted_df.empty:
        extracted_df = pd.DataFrame(columns=_MANIFEST_COLS + [
            "midrib_angle_deg", "petiole_length_px", "blade_length_px",
            "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
            "capitula_count", "involucre_height_px", "involucre_width_px",
            "involucre_height_mm", "involucre_width_mm", "AR_inv",
        ])
    else:
        cal_src_set = {str(s).strip().lower() for s in _CALIBRATED_SCALE_SOURCES}
        if "scale_source" not in extracted_df.columns:
            px_s = extracted_df.get("pixels_per_mm", pd.Series([np.nan] * len(extracted_df)))
            cal_s = extracted_df["is_calibrated_mm"] if "is_calibrated_mm" in extracted_df.columns else pd.Series([True] * len(extracted_df))
            cf_s = extracted_df["cf_source"] if "cf_source" in extracted_df.columns else pd.Series([""] * len(extracted_df))
            extracted_df["scale_source"] = [
                _resolve_scale_source(px, bool(ic) and (str(cf).strip().lower() in cal_src_set or not str(cf).strip()), cf)
                for px, ic, cf in zip(px_s, cal_s, cf_s)
            ]
        if "is_calibrated_mm" not in extracted_df.columns:
            extracted_df["is_calibrated_mm"] = [
                bool(pd.notna(px) and float(px) > 0 and str(src).strip().lower() in cal_src_set)
                for px, src in zip(extracted_df.get("pixels_per_mm", pd.Series([np.nan] * len(extracted_df))), extracted_df["scale_source"])
            ]
        if "plant_individual_id" not in extracted_df.columns:
            extracted_df["plant_individual_id"] = 0
        for pet_col in ("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio"):
            extracted_df[pet_col] = pd.to_numeric(extracted_df[pet_col], errors="coerce") if pet_col in extracted_df.columns else np.nan
        for c in _MANIFEST_COLS:
            if c not in extracted_df.columns:
                extracted_df[c] = False if c == "is_abaxial_pressed_leaf" else np.nan
        extracted_df = extracted_df[_MANIFEST_COLS + [c for c in extracted_df.columns if c not in _MANIFEST_COLS]]

    failed_df = pd.DataFrame(all_failed_records)
    if failed_df.empty:
        failed_df = pd.DataFrame(columns=["catalogNumber", "image_path", "failure_reason", "details"])

    repro_cols = [
        "catalogNumber", "plant_individual_id", "capitula_count", "corymb_head_count_est",
        "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm",
        "AR_inv", "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem", "radiate_head_fraction",
    ]
    repro_df = pd.DataFrame(_flatten_pid_records(pipe.voucher_repro_records, pipe.voucher_repro_by_pid))
    if repro_df.empty:
        repro_df = pd.DataFrame(columns=repro_cols)
    else:
        if "plant_individual_id" not in repro_df.columns:
            repro_df["plant_individual_id"] = 0
        repro_df["corymb_head_count_est"] = (
            pd.to_numeric(repro_df["corymb_head_count_est"], errors="coerce").fillna(0).astype(int)
            if "corymb_head_count_est" in repro_df.columns else 0
        )
        for col_f, def_v in (("plant_height_cm", np.nan), ("scape_to_basal_leaf_ratio", np.nan), ("is_bent_stem", False), ("radiate_head_fraction", np.nan)):
            if col_f not in repro_df.columns:
                repro_df[col_f] = def_v
        repro_df = repro_df[[c for c in repro_cols if c in repro_df.columns]]

    pheno_cols = [
        "catalogNumber", "plant_individual_id", "total_capitula", "n_anthesis", "n_fruiting",
        "n_bud", "voucher_phenological_state", "dominant_pappus_ratio", "radiate_head_fraction",
    ]
    pheno_df = pd.DataFrame(_flatten_pid_records(pipe.voucher_pheno_records, pipe.voucher_pheno_by_pid))
    if pheno_df.empty:
        pheno_df = pd.DataFrame(columns=pheno_cols)
    else:
        if "plant_individual_id" not in pheno_df.columns:
            pheno_df["plant_individual_id"] = 0
        if "radiate_head_fraction" not in pheno_df.columns:
            pheno_df["radiate_head_fraction"] = np.nan
        pheno_df = pheno_df[[c for c in pheno_cols if c in pheno_df.columns]]

    pipe._atomic_save_csv(extracted_df, manifest_path)
    pipe._atomic_save_csv(failed_df, failed_qc_path)
    pipe._atomic_save_csv(repro_df, repro_manifest_path)
    pipe._atomic_save_csv(pheno_df, pheno_manifest_path)
    return extracted_df, failed_df, repro_df, pheno_df


def load_resumption_state(
    pipe: Any,
    df_vouchers: pd.DataFrame,
    manifest_path: Path,
    failed_qc_path: Path,
    repro_manifest_path: Path,
    pheno_manifest_path: Path,
    base_is_completed_fn: Optional[Callable[..., Any]] = None,
) -> Dict[str, Any]:
    """Initialize SQLite connection, preload prior completed records, and determine wave vs. sequential mode."""
    existing_extracted_by_cat: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    existing_failed_by_cat: Dict[str, Dict[str, Any]] = {}
    existing_repro_by_cat: Dict[str, Dict[str, Any]] = {}
    existing_pheno_by_cat: Dict[str, Dict[str, Any]] = {}
    existing_repro_by_pid: Dict[str, Dict[int, Dict[str, Any]]] = defaultdict(dict)
    existing_pheno_by_pid: Dict[str, Dict[int, Dict[str, Any]]] = defaultdict(dict)

    use_db = (pipe.db_path is not None) or is_canonical_or_explicit_db_context(pipe.output_dir)
    db_conn = None
    if use_db:
        try:
            db_conn = get_db_connection(db_path=pipe.db_path, init_schema=True)
            if not df_vouchers.empty:
                upsert_vouchers_records(db_conn, df_vouchers.to_dict(orient="records"))
        except Exception as db_init_exc:
            logger.warning(f"SQLite initialization skipped ({db_init_exc}); falling back to CSV manifests.")
            db_conn = None

    completed_in_db: Set[str] = set()
    if not pipe.force:
        if db_conn is not None:
            try:
                target_cats = df_vouchers["catalogNumber"].astype(str).str.strip().tolist()
                completed_in_db = get_completed_qc_vouchers(db_conn, catalog_filter=target_cats)
                if completed_in_db:
                    try:
                        db_ext, db_fail, db_rep, db_phe, db_rep_pid, db_phe_pid = load_step02_records_from_db(
                            db_conn, list(completed_in_db), include_by_pid=True
                        )
                    except TypeError:
                        db_ext, db_fail, db_rep, db_phe = load_step02_records_from_db(db_conn, list(completed_in_db))
                        db_rep_pid, db_phe_pid = getattr(db_rep, "by_pid", {}) or {}, getattr(db_phe, "by_pid", {}) or {}
                    for k_cat, recs in db_ext.items():
                        existing_extracted_by_cat[k_cat].extend(recs)
                    existing_failed_by_cat.update(db_fail)
                    existing_repro_by_cat.update(db_rep)
                    existing_pheno_by_cat.update(db_phe)
                    for src_map, dst_map in ((db_rep_pid, existing_repro_by_pid), (db_phe_pid, existing_pheno_by_pid)):
                        for k_key, p_val in (src_map or {}).items():
                            if isinstance(k_key, tuple) and len(k_key) == 2 and isinstance(p_val, dict):
                                dst_map[str(k_key[0]).strip()][int(k_key[1])] = dict(p_val)
                            elif isinstance(p_val, dict):
                                dst_map[str(k_key).strip()].update(
                                    {int(pk): dict(pv) for pk, pv in p_val.items() if isinstance(pv, dict)}
                                )
            except Exception as db_load_exc:
                logger.warning(f"Could not query completed vouchers from SQLite: {db_load_exc}")

        if manifest_path.exists() and manifest_path.stat().st_size > 0:
            try:
                prev_manifest = pd.read_csv(manifest_path)
                if not prev_manifest.empty and "catalogNumber" in prev_manifest.columns:
                    for rec in prev_manifest.to_dict(orient="records"):
                        c_key = str(rec["catalogNumber"]).strip()
                        if c_key not in existing_extracted_by_cat:
                            existing_extracted_by_cat[c_key].append(rec)
            except Exception as e:
                logger.warning(f"Could not load previous manifest at {manifest_path}: {e}")

        if failed_qc_path.exists() and failed_qc_path.stat().st_size > 0:
            try:
                prev_failed = pd.read_csv(failed_qc_path)
                if not prev_failed.empty and "catalogNumber" in prev_failed.columns:
                    for rec in prev_failed.to_dict(orient="records"):
                        existing_failed_by_cat.setdefault(str(rec["catalogNumber"]).strip(), rec)
            except Exception as e:
                logger.warning(f"Could not load previous failed QC table at {failed_qc_path}: {e}")

        for csv_p, base_d, pid_d, label in (
            (repro_manifest_path, existing_repro_by_cat, existing_repro_by_pid, "reproductive"),
            (pheno_manifest_path, existing_pheno_by_cat, existing_pheno_by_pid, "phenology"),
        ):
            if csv_p.exists() and csv_p.stat().st_size > 0:
                try:
                    prev_df = pd.read_csv(csv_p)
                    if not prev_df.empty and "catalogNumber" in prev_df.columns:
                        db_cats = set(base_d.keys())
                        for rec in prev_df.to_dict(orient="records"):
                            c_r = str(rec["catalogNumber"]).strip()
                            base_d.setdefault(c_r, rec)
                            if c_r not in db_cats:
                                pid_r = int(rec["plant_individual_id"]) if pd.notna(rec.get("plant_individual_id")) else 0
                                pid_d[c_r].setdefault(pid_r, rec)
                except Exception as e:
                    logger.warning(f"Could not load previous {label} table at {csv_p}: {e}")

    total_vouchers = len(df_vouchers)
    completed_on_disk = (
        completed_in_db | (pipe._index_completed_vouchers() if len(completed_in_db) < total_vouchers else set())
        if not pipe.force else set()
    )
    custom_process_voucher = "process_voucher" in pipe.__dict__
    custom_is_completed = "is_voucher_completed" in pipe.__dict__ or (
        base_is_completed_fn is not None and getattr(pipe.is_voucher_completed, "__func__", None) is not base_is_completed_fn
    )
    custom_hooks = (
        custom_process_voucher
        or "predict" in getattr(pipe.stage1_detector, "__dict__", {})
        or "predict_batch" in getattr(pipe.stage1_detector, "__dict__", {})
        or "predict" in getattr(pipe.model_engine, "__dict__", {})
        or "predict_batch" in getattr(pipe.model_engine, "__dict__", {})
    )
    return {
        "db_conn": db_conn, "existing_extracted_by_cat": existing_extracted_by_cat,
        "existing_failed_by_cat": existing_failed_by_cat, "existing_repro_by_cat": existing_repro_by_cat,
        "existing_pheno_by_cat": existing_pheno_by_cat, "existing_repro_by_pid": existing_repro_by_pid,
        "existing_pheno_by_pid": existing_pheno_by_pid, "completed_on_disk": completed_on_disk,
        "custom_process_voucher": custom_process_voucher, "custom_is_completed": custom_is_completed,
        "use_wave_engine": (not custom_hooks) and (total_vouchers > 4),
    }


def prefetch_voucher_data(pipe: Any, cat: str, img_p: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Prefetch voucher image via `cv2.imread` and optionally precompute Stage 1 letterbox."""
    if not img_p.exists():
        for ext in (".webp", ".jpg", ".jpeg", ".png"):
            cand_p = img_p.with_suffix(ext)
            if cand_p.exists():
                img_p = cand_p
                break
    if not img_p.exists():
        return None, None
    img = cv2.imread(str(img_p))
    if img is None:
        return None, None
    s1_lb: Optional[np.ndarray] = None
    if (
        hasattr(pipe.stage1_detector, "prepare_letterbox")
        and getattr(pipe.stage1_detector, "model", None) is not None
        and (pipe.force or cat not in getattr(pipe.stage1_detector, "cached_detections", {}))
    ):
        try:
            s1_lb = pipe.stage1_detector.prepare_letterbox(img)
        except Exception:
            s1_lb = None
    return img, s1_lb


def should_skip_voucher(pipe: Any, cat_num: str, state: Dict[str, Any], manifest_path: Path, total_vouchers: int) -> bool:
    """Return True if `cat_num` already has completed artifacts or prior failed-only record."""
    if pipe.force:
        return False
    has_disk_artifacts = (cat_num in state["completed_on_disk"]) or (
        pipe.is_voucher_completed(cat_num) if state["custom_is_completed"] else False
    )
    has_manifest_entry = (cat_num in state["existing_extracted_by_cat"]) or (cat_num in state["existing_failed_by_cat"])
    is_prev_failed_only = (cat_num in state["existing_failed_by_cat"]) and (cat_num not in state["existing_extracted_by_cat"]) and (total_vouchers > 20)
    return bool(
        (has_disk_artifacts and (has_manifest_entry or state["custom_process_voucher"] or (not manifest_path.exists() and total_vouchers <= 20)))
        or is_prev_failed_only
    )


def retire_skipped_voucher(
    pipe: Any, cat_num: str, state: Dict[str, Any],
    all_extracted_records: List[Dict[str, Any]], all_failed_records: List[Dict[str, Any]],
    processed_count: int, skipped_count: int, total_vouchers: int,
) -> int:
    """Retire a completed skipped voucher into active manifests and pipeline state dicts."""
    skipped_count += 1
    logger.info(
        f"Skipping completed voucher: {cat_num} (artifacts exist in contours/masks). "
        f"[Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}]"
    )
    if cat_num in state["existing_extracted_by_cat"]:
        all_extracted_records.extend(state["existing_extracted_by_cat"][cat_num])
    elif cat_num in state["existing_failed_by_cat"]:
        all_failed_records.append(state["existing_failed_by_cat"][cat_num])
    if cat_num in state["existing_repro_by_cat"]:
        pipe.voucher_repro_records[cat_num] = state["existing_repro_by_cat"][cat_num]
    if cat_num in state["existing_pheno_by_cat"]:
        pipe.voucher_pheno_records[cat_num] = state["existing_pheno_by_cat"][cat_num]
    if cat_num in state["existing_repro_by_pid"] and state["existing_repro_by_pid"][cat_num]:
        pipe.voucher_repro_by_pid[cat_num] = dict(state["existing_repro_by_pid"][cat_num])
    if cat_num in state["existing_pheno_by_pid"] and state["existing_pheno_by_pid"][cat_num]:
        pipe.voucher_pheno_by_pid[cat_num] = dict(state["existing_pheno_by_pid"][cat_num])
    return skipped_count
