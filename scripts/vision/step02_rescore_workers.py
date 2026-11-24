#!/usr/bin/env python3
"""
Step 02 rescore workers: path resolvers, patch alignment, picklable CPU QC &
ruler-calibration workers, and batch execution helpers for `rescore_step02_qc_and_zones.py`.
"""

from __future__ import annotations

import inspect, logging, sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2, numpy as np, pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.leaf_geometry_utils import compute_lm3_shape_metrics
from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES
from scripts.vision.lm3_inference_engines import LM3ScaleCalibrator, evaluate_leaf_damage_and_fold_qc, solve_ruler_tick_lattice_px_per_mm
from scripts.vision.step02_leaf_traits_utils import extract_packera_leaf_margin_and_base_traits, measure_enclosed_interlobe_void_frac

logger = logging.getLogger("RescoreStep02")
_FACADE_MOD = "scripts.vision.rescore_step02_qc_and_zones"
_ORIG_GLOBALS: Dict[str, Any] = {}

_QC_INPUT_COLS: Tuple[str, ...] = (
    "catalogNumber", "leaf_id", "mask_path", "solidity", "circularity",
    "bilateral_dice", "si_a", "lm3_gates_pass", "n_lamina_components",
    "detector_class_id", "hole_area_frac", "enclosed_void_frac", "is_crop_clipped",
    "mask_is_upright", "pixels_per_mm", "is_calibrated_mm", "scale_source",
    "voucher_px_per_mm", "voucher_is_cal", "basal_cordate_index",
    "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac",
)
_QC_RESULT_COLS: Tuple[str, ...] = (
    "solidity", "circularity", "enclosed_void_frac", "hole_area_frac",
    "is_crop_clipped", "detector_class_id", "bilateral_lobe_symmetry",
    "efa_qc_pass", "qc_flag_reason", "basal_cordate_index",
    "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac",
)


def _resolve_hook(name: str, local_val: Any, orig_globals: Optional[Dict[str, Any]] = None) -> Any:
    """Resolve facade or local test-patched symbol without introducing circular imports."""
    orig_map = orig_globals if orig_globals is not None else _ORIG_GLOBALS
    if local_val is not orig_map.get(name, local_val):
        return local_val
    mod = sys.modules.get(_FACADE_MOD)
    return getattr(mod, name, local_val) if mod is not None else local_val


def _resolve_mask_file(raw_mask_path: Any, data_dir: Path) -> Optional[Path]:
    s = str(raw_mask_path or "").strip()
    if not s or s.lower() in {"nan", "none", "na"}:
        return None
    proj_root, p = _resolve_hook("PROJECT_ROOT", PROJECT_ROOT), Path(s)
    return next((c for c in (p, data_dir / "masks" / p.name, proj_root / s, proj_root / "data" / "masks" / p.name) if c.is_file()), None)


def _resolve_voucher_image_file(raw_image_path: Any, catalog_number: str, data_dir: Path) -> Optional[Path]:
    proj_root = _resolve_hook("PROJECT_ROOT", PROJECT_ROOT)
    s, candidates = str(raw_image_path or "").strip(), []
    if s and s.lower() not in {"nan", "none", "na"}:
        p = Path(s)
        candidates.extend([p, proj_root / s, data_dir / "raw_vouchers" / p.name])
    cat_safe = str(catalog_number or "").strip()
    if cat_safe:
        for ext in (".webp", ".jpg", ".jpeg", ".png"):
            candidates.extend([data_dir / "raw_vouchers" / f"{cat_safe}{ext}", proj_root / "data" / "raw_vouchers" / f"{cat_safe}{ext}"])
    return next((c for c in candidates if c.is_file()), None)


def _resolve_crop_file(raw_crop_path: Any, data_dir: Path) -> Optional[Path]:
    s = str(raw_crop_path or "").strip()
    if not s or s.lower() in {"nan", "none", "na"}:
        return None
    proj_root, p = _resolve_hook("PROJECT_ROOT", PROJECT_ROOT), Path(s)
    return next((c for c in (p, proj_root / s, data_dir / "crops" / p.name, proj_root / "data" / "crops" / p.name) if c.is_file()), None)


def _align_rotated_patch_to_mask(patch_bgr: np.ndarray, target_hw: Tuple[int, int], pad_val: int = 255) -> np.ndarray:
    """Align a rotated sheet crop to an upright mask by 1:1 center-crop/pad or scale-matched resize."""
    th, tw = int(target_hw[0]), int(target_hw[1])
    ph, pw = int(patch_bgr.shape[0]), int(patch_bgr.shape[1])
    if (ph, pw) == (th, tw):
        return patch_bgr
    if th <= ph + 20 and tw <= pw + 20 and ph <= max(int(2.2 * th), th + 48) and pw <= max(int(2.2 * tw), tw + 48):
        canvas = np.full((th, tw, patch_bgr.shape[2]), pad_val, dtype=patch_bgr.dtype)
        sy0, sx0 = max(0, (ph - th) // 2), max(0, (pw - tw) // 2)
        dy0, dx0 = max(0, (th - ph) // 2), max(0, (tw - pw) // 2)
        copy_h, copy_w = min(ph - sy0, th - dy0), min(pw - sx0, tw - dx0)
        if copy_h > 0 and copy_w > 0:
            canvas[dy0 : dy0 + copy_h, dx0 : dx0 + copy_w] = patch_bgr[sy0 : sy0 + copy_h, sx0 : sx0 + copy_w]
        return canvas
    return cv2.resize(patch_bgr, (tw, th), interpolation=cv2.INTER_AREA)


def _eval_single_leaf_qc_worker(task: Tuple[Dict[str, Any], str, bool]) -> Dict[str, Any]:
    """Read-only evaluation of enclosed voids, Packera margin/base traits, solidity/circ, and EFA QC gates."""
    row, data_dir_str, seg_partial = task
    resolve_mask = _resolve_hook("_resolve_mask_file", _resolve_mask_file)
    measure_void = _resolve_hook("measure_enclosed_interlobe_void_frac", measure_enclosed_interlobe_void_frac)
    extract_traits = _resolve_hook("extract_packera_leaf_margin_and_base_traits", extract_packera_leaf_margin_and_base_traits)
    compute_shape = _resolve_hook("compute_lm3_shape_metrics", compute_lm3_shape_metrics)
    eval_qc = _resolve_hook("evaluate_leaf_damage_and_fold_qc", evaluate_leaf_damage_and_fold_qc)

    m_file = resolve_mask(row.get("mask_path"), Path(data_dir_str))
    mask_u8: Optional[np.ndarray] = None
    if m_file is not None:
        m_raw = cv2.imread(str(m_file), cv2.IMREAD_GRAYSCALE)
        if m_raw is not None and np.count_nonzero(m_raw) > 0:
            mask_u8 = ((m_raw > 0).astype(np.uint8)) * 255

    def _fval(key: str, default: float) -> float:
        v = row.get(key)
        return float(v) if v is not None and pd.notna(v) else default

    def _ival(key: str, default: int) -> int:
        v = row.get(key)
        return int(v) if v is not None and pd.notna(v) else default

    sol_val, circ_val = _fval("solidity", 0.85), _fval("circularity", float("nan"))
    bd_val, sia_val = _fval("bilateral_dice", 0.85), _fval("si_a", 0.0)
    gp_val, ncc_val, dcid_val = bool(_ival("lm3_gates_pass", 1)), _ival("n_lamina_components", 1), _ival("detector_class_id", 0)
    haf_val, clipped_val, evf_val = _fval("hole_area_frac", 0.0), bool(_ival("is_crop_clipped", 0)), _fval("enclosed_void_frac", 0.0)
    packera_traits: Dict[str, float] = {k: _fval(k, float("nan")) for k in ("basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac")}

    if mask_u8 is not None:
        hm, wm = mask_u8.shape[:2]
        edge_allow = max(6, int(round(0.06 * max(hm, wm))))
        if any(int(np.count_nonzero(b)) > edge_allow for b in (mask_u8[0, :], mask_u8[-1, :], mask_u8[:, 0], mask_u8[:, -1])):
            clipped_val = True
        evf_val = float(measure_void(mask_u8, hole_mask=None))
        packera_traits = extract_traits(mask_u8, oriented_points=None)
        shape_m = compute_shape(mask_u8)
        if shape_m["area_px"] > 0.0:
            sol_val, circ_val = round(float(shape_m["solidity"]), 4), round(float(shape_m["circularity"]), 4)

    qc_res = eval_qc(
        oriented_mask_u8=mask_u8, solidity=sol_val, bilateral_dice=bd_val, si_a=sia_val,
        lm3_gates_pass=gp_val, n_lamina_components=ncc_val, detector_class_id=dcid_val,
        segment_partial_leaves=seg_partial, is_crop_clipped=clipped_val,
        hole_area_frac=haf_val, enclosed_void_frac=evf_val, is_geometric_fallback_pose=True,
    )
    return {
        "solidity": sol_val, "circularity": circ_val,
        "enclosed_void_frac": round(float(evf_val), 4), "hole_area_frac": round(float(haf_val), 4),
        "is_crop_clipped": 1 if clipped_val else 0, "detector_class_id": dcid_val,
        "bilateral_lobe_symmetry": 1 if qc_res["bilateral_lobe_symmetry"] else 0,
        "efa_qc_pass": 1 if qc_res["efa_qc_pass"] else 0,
        "qc_flag_reason": str(qc_res["qc_flag_reason"]), "min_rms": float(qc_res["min_margin_rms_norm"]),
        "basal_cordate_index": packera_traits.get("basal_cordate_index", float("nan")),
        "margin_serration_index": packera_traits.get("margin_serration_index", float("nan")),
        "margin_tooth_acuteness": packera_traits.get("margin_tooth_acuteness", float("nan")),
        "terminal_lobe_area_frac": packera_traits.get("terminal_lobe_area_frac", float("nan")),
    }


def _calibrate_voucher_ruler_worker(task: Tuple[str, str, Optional[float]]) -> Optional[Tuple[str, float]]:
    """Run CPU unboxed margin metric tick-lattice scan for one uncalibrated voucher sheet."""
    cat, img_path_str, mp_prior_mm = task
    img = cv2.imread(img_path_str, cv2.IMREAD_COLOR)
    if img is None or img.size == 0:
        return None
    solve_fn = _resolve_hook("solve_ruler_tick_lattice_px_per_mm", solve_ruler_tick_lattice_px_per_mm)
    if "verified_ruler_crop" in inspect.signature(solve_fn).parameters:
        px_per_mm = solve_fn(img, mp_prior_px_per_mm=mp_prior_mm, verified_ruler_crop=False)
    else:
        px_per_mm = solve_fn(img, mp_prior_px_per_mm=mp_prior_mm)
    if px_per_mm is None or not np.isfinite(float(px_per_mm)) or float(px_per_mm) <= 1.0:
        return None
    return (cat, round(float(px_per_mm), 4))


def _apply_scale_updates_and_normalize(el_df: pd.DataFrame, ruler_cal_scales: Dict[str, float], unboxed_scales: Dict[str, float], cal_src_set: set[str]) -> None:
    """Apply newly solved ruler/unboxed scales onto `el_df` and normalize `scale_source` provenance."""
    def _num_at(idx_c: Any, col: str) -> Any:
        return pd.to_numeric(el_df.at[idx_c, col], errors="coerce") if col in el_df.columns else np.nan

    if ruler_cal_scales:
        cat_series = el_df["catalogNumber"].astype(str).str.strip()
        for idx_c in el_df.index[cat_series.isin(ruler_cal_scales)]:
            new_px = float(ruler_cal_scales[str(el_df.at[idx_c, "catalogNumber"]).strip()])
            el_df.at[idx_c, "pixels_per_mm"], el_df.at[idx_c, "is_calibrated_mm"], el_df.at[idx_c, "scale_source"] = new_px, 1, "ruler_calibrated"
            if "voucher_px_per_mm" in el_df.columns:
                el_df.at[idx_c, "voucher_px_per_mm"] = new_px
            if "voucher_is_cal" in el_df.columns:
                el_df.at[idx_c, "voucher_is_cal"] = 1
            for px_col, mm_col, denom, dec, allow_zero in (
                ("area_px", "area_mm2", new_px ** 2, 2, False),
                ("blade_length_px", "blade_length_mm", new_px, 3, False),
                ("petiole_length_px", "petiole_length_mm", new_px, 3, True),
                ("petiole_width_px", "petiole_width_mm", new_px, 4, False),
            ):
                val = _num_at(idx_c, px_col)
                if pd.notna(val) and (float(val) >= 0.0 if allow_zero else float(val) > 0.0):
                    el_df.at[idx_c, mm_col] = round(float(val) / denom, dec)
            blpx, bwpx = _num_at(idx_c, "blade_length_px"), _num_at(idx_c, "blade_width_px")
            if (pd.isna(bwpx) or float(bwpx) <= 0.0) and pd.notna(blpx) and float(blpx) > 0.0 and "aspect_ratio" in el_df.columns:
                ar_v = _num_at(idx_c, "aspect_ratio")
                if pd.notna(ar_v) and float(ar_v) > 1e-6:
                    bwpx = float(blpx) / float(ar_v)
            if pd.notna(bwpx) and float(bwpx) > 0.0:
                el_df.at[idx_c, "blade_width_mm"] = round(float(bwpx) / new_px, 3)

    if unboxed_scales:
        cat_series = el_df["catalogNumber"].astype(str).str.strip()
        for idx_c in el_df.index[cat_series.isin(unboxed_scales) & (~cat_series.isin(ruler_cal_scales))]:
            el_df.at[idx_c, "pixels_per_mm"] = float(unboxed_scales[str(el_df.at[idx_c, "catalogNumber"]).strip()])
            el_df.at[idx_c, "is_calibrated_mm"], el_df.at[idx_c, "scale_source"] = 0, "tick_lattice_unboxed"

    is_cal_arr = pd.to_numeric(el_df["is_calibrated_mm"], errors="coerce").fillna(0).astype(int)
    px_arr = pd.to_numeric(el_df["pixels_per_mm"], errors="coerce")
    ss_low = el_df["scale_source"].fillna("").astype(str).str.strip().str.lower()
    is_cal_arr = ((is_cal_arr == 1) & px_arr.notna() & (px_arr > 0.0) & (ss_low.isin(cal_src_set) | (ss_low == ""))).astype(int)
    el_df["is_calibrated_mm"] = is_cal_arr
    el_df.loc[(is_cal_arr == 1) & (~ss_low.isin(cal_src_set)), "scale_source"] = "ruler_calibrated"
    el_df.loc[ss_low == "tick_lattice_unboxed", "scale_source"] = "tick_lattice_unboxed"
    el_df.loc[ss_low == "tick_lattice_unboxed", "is_calibrated_mm"] = 0
    mp_fix = (el_df["is_calibrated_mm"] == 0) & px_arr.notna() & (px_arr > 0.0) & (~ss_low.isin(cal_src_set | {"mp_estimated", "tick_lattice_unboxed"}))
    el_df.loc[mp_fix, "scale_source"] = "mp_estimated"


def run_ruler_calibration_backfill(el_df: pd.DataFrame, eff_data_dir: Path, n_workers: int) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Backfill CPU archival ruler calibration and unboxed tick-lattice scales onto `el_df`."""
    log, resolve_img = _resolve_hook("logger", logger), _resolve_hook("_resolve_voucher_image_file", _resolve_voucher_image_file)
    cal_cls, cal_worker = _resolve_hook("LM3ScaleCalibrator", LM3ScaleCalibrator), _resolve_hook("_calibrate_voucher_ruler_worker", _calibrate_voucher_ruler_worker)
    proc_pool, thread_pool = _resolve_hook("ProcessPoolExecutor", ProcessPoolExecutor), _resolve_hook("ThreadPoolExecutor", ThreadPoolExecutor)

    cal_src_set = {str(s).strip().lower() for s in _CALIBRATED_SCALE_SOURCES}
    for col_d, val_d in (("scale_source", ""), ("is_calibrated_mm", 0), ("pixels_per_mm", np.nan)):
        if col_d not in el_df.columns:
            el_df[col_d] = val_d

    ruler_tasks: List[Tuple[str, str, Optional[float]]] = []
    seen_cats: set[str] = set()
    for _, rv in el_df[["catalogNumber", "image_path", "voucher_px_per_mm", "voucher_is_cal"]].drop_duplicates("catalogNumber").iterrows():
        cat_s = str(rv["catalogNumber"]).strip()
        if not cat_s or cat_s in seen_cats:
            continue
        seen_cats.add(cat_s)
        is_cal_v = int(rv["voucher_is_cal"]) if pd.notna(rv.get("voucher_is_cal")) else 0
        cat_mask = el_df["catalogNumber"].astype(str).str.strip() == cat_s
        l_ss = el_df.loc[cat_mask, "scale_source"].fillna("").astype(str).str.strip().str.lower()
        l_cal = pd.to_numeric(el_df.loc[cat_mask, "is_calibrated_mm"], errors="coerce").fillna(0).astype(int)
        if is_cal_v != 0 and bool(((l_cal == 1) & (l_ss.isin(cal_src_set) | (l_ss == ""))).any()):
            continue
        img_f = resolve_img(rv.get("image_path"), cat_s, eff_data_dir)
        if img_f is not None:
            v_px = float(rv["voucher_px_per_mm"]) if pd.notna(rv.get("voucher_px_per_mm")) else None
            ruler_tasks.append((cat_s, str(img_f), v_px))

    ruler_cal_scales: Dict[str, float] = {}
    unboxed_scales: Dict[str, float] = {}

    def _run_unboxed_tasks(tasks: List[Tuple[str, str, Optional[float]]]) -> None:
        if len(tasks) >= 8 and n_workers > 1:
            with proc_pool(max_workers=n_workers) as pool:
                hits = list(pool.map(cal_worker, tasks, chunksize=8))
        else:
            hits = [cal_worker(t) for t in tasks]
        for h in hits:
            if h is not None:
                unboxed_scales[h[0]] = h[1]

    if ruler_tasks:
        calibrator: Optional[Any] = None
        try:
            calibrator = cal_cls(device="cpu", use_cf_predicted_by_mp=False)
        except Exception as cal_exc:
            log.warning(f"LM3ScaleCalibrator CPU init skipped ({cal_exc}); using unboxed tick-lattice scan.")

        if calibrator is not None:
            def _cal_sheet_task(task: Tuple[str, str, Optional[float]]) -> Tuple[str, str, Optional[float], bool, Tuple[str, str, Optional[float]]]:
                cat_t, img_p_t, _v_px_t = task
                img_c = cv2.imread(img_p_t, cv2.IMREAD_COLOR)
                if img_c is None or img_c.size == 0:
                    return (cat_t, "uncalibrated", None, False, task)
                try:
                    res_c = calibrator.calibrate_sheet(img_c, catalog_number=cat_t, use_cf_predicted_by_mp=False)
                    ss_c = str(res_c.get("scale_source") or "").strip().lower()
                    cal_c, px_c = bool(res_c.get("is_calibrated_mm")), res_c.get("pixels_per_mm")
                    if px_c is not None and np.isfinite(float(px_c)) and float(px_c) > 1.0:
                        return (cat_t, ss_c, round(float(px_c), 4), cal_c, task)
                    return (cat_t, ss_c, None, False, task)
                except Exception:
                    return (cat_t, "error", None, False, task)

            if (
                len(ruler_tasks) >= 8 and n_workers > 1 and getattr(calibrator, "archival_session", None) is None
                and "calibrate_sheet" not in getattr(calibrator, "__dict__", {}) and "_detect_ruler_boxes" not in getattr(calibrator, "__dict__", {})
                and getattr(type(calibrator), "calibrate_sheet", None) is cal_cls.calibrate_sheet
            ):
                _run_unboxed_tasks(ruler_tasks)
            else:
                if len(ruler_tasks) >= 8 and n_workers > 1:
                    with thread_pool(max_workers=max(1, min(32, (n_workers + 1) // 2))) as tpool:
                        cal_hits = list(tpool.map(_cal_sheet_task, ruler_tasks))
                else:
                    cal_hits = [_cal_sheet_task(t) for t in ruler_tasks]
                fallback_tasks: List[Tuple[str, str, Optional[float]]] = []
                for cat_r, ss_r, px_r, is_cal_r, orig_task in cal_hits:
                    if ss_r == "ruler_calibrated" and is_cal_r and px_r is not None:
                        ruler_cal_scales[cat_r] = px_r
                    elif ss_r == "tick_lattice_unboxed" and px_r is not None:
                        unboxed_scales[cat_r] = px_r
                    elif ss_r == "error":
                        fallback_tasks.append(orig_task)
                for fb_t in fallback_tasks:
                    hit_fb = cal_worker(fb_t)
                    if hit_fb is not None:
                        unboxed_scales[hit_fb[0]] = hit_fb[1]
        else:
            _run_unboxed_tasks(ruler_tasks)
        log.info(f"CPU scale backfill: ruler_calibrated={len(ruler_cal_scales)}, tick_lattice_unboxed={len(unboxed_scales)} / {len(ruler_tasks)} uncalibrated vouchers.")

    _apply_scale_updates_and_normalize(el_df, ruler_cal_scales, unboxed_scales, cal_src_set)
    return ruler_cal_scales, unboxed_scales


def run_leaf_qc_evaluation(el_df: pd.DataFrame, eff_data_dir: Path, seg_partial: bool, n_workers: int) -> None:
    """Evaluate read-only leaf QC gates and Packera traits across `el_df` in parallel."""
    qc_worker, proc_pool = _resolve_hook("_eval_single_leaf_qc_worker", _eval_single_leaf_qc_worker), _resolve_hook("ProcessPoolExecutor", ProcessPoolExecutor)
    avail_qc_cols = [c for c in _QC_INPUT_COLS if c in el_df.columns]
    qc_tasks = [(rec, str(eff_data_dir), seg_partial) for rec in el_df[avail_qc_cols].to_dict(orient="records")]
    if len(el_df) >= 64 and n_workers > 1:
        with proc_pool(max_workers=n_workers) as pool:
            qc_results = list(pool.map(qc_worker, qc_tasks, chunksize=64))
    else:
        qc_results = [qc_worker(t) for t in qc_tasks]
    for col in _QC_RESULT_COLS:
        el_df[col] = [r[col] for r in qc_results]
    el_df["_min_rms"] = [r["min_rms"] for r in qc_results]


_ORIG_GLOBALS.update(globals())
