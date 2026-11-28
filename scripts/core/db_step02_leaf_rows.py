"""Row-preparation helpers for atomic Step 02 SQLite segmentation transactions."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from scripts.core.db_column_constants import _to_bool_int, _to_float, _to_int, _to_str

_EXTRACTED_LEAVES_INSERT_COLS: Tuple[str, ...] = (
    "catalogNumber", "plant_individual_id", "leaf_id", "shape_id",
    "solidity", "circularity", "extent", "aspect_ratio",
    "area_px", "area_mm2", "midrib_angle_deg", "petiole_length_px",
    "petiole_length_mm", "petiole_width_px", "petiole_width_mm", "blade_length_px",
    "blade_length_mm", "blade_width_px", "blade_width_mm", "petiole_to_blade_ratio",
    "pixels_per_mm", "is_calibrated_mm", "scale_source", "mask_path",
    "mask_is_upright", "crop_path", "contour_path", "ect_path",
    "contour_n_pts", "contour_xy_json", "leaf_zone", "leaf_zone_confidence",
    "dist_from_crown_norm", "cy_rel_plant_bbox", "archetype_score", "bilateral_dice",
    "si_a", "lm3_gates_pass", "n_lamina_components", "detector_class_id",
    "is_crop_clipped", "hole_area_frac", "enclosed_void_frac", "bilateral_lobe_symmetry",
    "efa_qc_pass", "qc_flag_reason", "blade_cie_l_norm", "blade_cie_l_sd",
    "blade_cie_l_sheet_delta", "blade_cie_a_delta", "blade_cie_b_delta", "bbox_x1",
    "bbox_y1", "bbox_x2", "bbox_y2", "base_angle_deg",
    "apex_angle_deg", "lamina_curvature", "lamina_trace_length_px", "basal_cordate_index",
    "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "is_abaxial_pressed_leaf", "abaxial_cie_l_sheet_delta",
)

_UNCALIBRATED_SCALE_SOURCES = frozenset({"mp_estimated", "predicted_from_megapixels", "uncalibrated"})


def _p75_leaf_metric(leaves: Sequence[Mapping[str, Any]], key: str) -> Optional[float]:
    vals_clean = sorted(v for v in (_to_float(x.get(key)) for x in leaves) if v is not None)
    if not vals_clean:
        return None
    return vals_clean[min(len(vals_clean) - 1, int(0.75 * (len(vals_clean) - 1)))]


def _resolve_voucher_calibration_and_qc(
    extracted_records: Sequence[Mapping[str, Any]],
    failure_record: Optional[Mapping[str, Any]],
    image_path: Optional[str],
) -> Tuple[Optional[float], Optional[int], str, Optional[str], Optional[str], Optional[str]]:
    """Determine voucher-level calibration & QC status (preserving is_calibrated_mm == 0 for mp_estimated)."""
    ppm_val: Optional[float] = None
    is_cal_val: Optional[int] = None
    for r in extracted_records:
        p = _to_float(r.get("pixels_per_mm"))
        c_flag = _to_bool_int(r.get("is_calibrated_mm"))
        s_src = (_to_str(r.get("scale_source")) or _to_str(r.get("cf_source")) or "").strip().lower()
        is_mp_est = s_src in _UNCALIBRATED_SCALE_SOURCES
        if p is not None and p > 0:
            ppm_val = p
            is_cal_val = (0 if is_mp_est else c_flag) if c_flag is not None else (0 if is_mp_est else 1)
            break
        if c_flag is not None and is_cal_val is None:
            is_cal_val = 0 if is_mp_est else c_flag
    if ppm_val is None and failure_record:
        p_fail = _to_float(failure_record.get("pixels_per_mm"))
        c_fail = _to_bool_int(failure_record.get("is_calibrated_mm"))
        s_fail = (_to_str(failure_record.get("scale_source")) or _to_str(failure_record.get("cf_source")) or "").strip().lower()
        is_mp_fail = s_fail in _UNCALIBRATED_SCALE_SOURCES
        if p_fail is not None and p_fail > 0:
            ppm_val = p_fail
            is_cal_val = 0 if is_mp_fail else (c_fail if c_fail is not None else 1)
        elif is_cal_val is None and c_fail is not None:
            is_cal_val = 0 if is_mp_fail else c_fail

    qc_status = "PASS" if len(extracted_records) > 0 else "FAILED_QC"
    fail_reason = _to_str(failure_record.get("failure_reason")) if failure_record else None
    fail_details = _to_str(failure_record.get("details") or failure_record.get("failure_details")) if failure_record else None
    eff_img_path = _to_str(image_path) or (_to_str(failure_record.get("image_path")) if failure_record else None)
    return ppm_val, is_cal_val, qc_status, fail_reason, fail_details, eff_img_path


def _build_step02_plant_and_repro_rows(
    cat: str,
    extracted_records: Sequence[Mapping[str, Any]],
    ppm_val: Optional[float],
    repro: Mapping[str, Any],
    pheno: Mapping[str, Any],
    repro_by_pid: Optional[Mapping[int, Mapping[str, Any]]],
    pheno_by_pid: Optional[Mapping[int, Mapping[str, Any]]],
    crowns: Mapping[int, Tuple[float, float]],
    centroids: Mapping[int, Tuple[float, float]],
    petiole_by_pid: Mapping[int, Mapping[str, Any]],
) -> Tuple[List[Tuple[Any, ...]], List[Tuple[Any, ...]]]:
    """Assemble per-individual `plant_individuals` and `voucher_reproductive_phenology` tuples."""
    pids: Set[int] = set()
    for r in extracted_records:
        pids.add(_to_int(r.get("plant_individual_id"), default=0) or 0)
    for mapping in (crowns, centroids, repro_by_pid or {}, pheno_by_pid or {}):
        for k_pid in mapping.keys():
            pids.add(int(k_pid))
    if not pids:
        pids.add(0)

    plant_rows: List[Tuple[Any, ...]] = []
    repro_pheno_rows: List[Tuple[Any, ...]] = []
    repro_pid_map: Dict[int, Mapping[str, Any]] = {int(k): v for k, v in (repro_by_pid or {}).items() if v is not None}
    pheno_pid_map: Dict[int, Mapping[str, Any]] = {int(k): v for k, v in (pheno_by_pid or {}).items() if v is not None}
    is_multi_pid = len([p for p in pids if p != 0]) > 1 or (len(pids) > 1 and 0 not in pids)

    leaves_by_pid: Dict[int, List[Mapping[str, Any]]] = {pid: [] for pid in pids}
    for r in extracted_records:
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        leaves_by_pid.setdefault(pid, []).append(r)

    for pid in sorted(pids):
        rc_xy = centroids.get(pid)
        cr_xy = crowns.get(pid)
        p_pet_info = dict(petiole_by_pid.get(pid) or {})
        pid_repro: Mapping[str, Any] = repro_pid_map.get(pid) if pid in repro_pid_map else ({} if (is_multi_pid and repro_pid_map) else repro)
        pid_pheno: Mapping[str, Any] = pheno_pid_map.get(pid) if pid in pheno_pid_map else ({} if (is_multi_pid and pheno_pid_map) else pheno)

        p_height_cm = _to_float(pid_repro.get("plant_height_cm"))
        if p_height_cm is None:
            p_height_cm = _to_float(p_pet_info.get("plant_height_cm"))
        p_height_mm = (p_height_cm * 10.0) if p_height_cm is not None else _to_float(pid_repro.get("plant_height_mm"))
        p_height_px = _to_float(pid_repro.get("plant_height_px"))
        if p_height_px is None:
            p_height_px = _to_float(p_pet_info.get("scape_length_px"))
        if p_height_px is None and p_height_mm is not None and ppm_val:
            p_height_px = p_height_mm * ppm_val
        scape_ratio = _to_float(pid_repro.get("scape_to_basal_leaf_ratio"))
        if scape_ratio is None:
            scape_ratio = _to_float(p_pet_info.get("scape_to_basal_leaf_ratio"))
        is_bent = _to_bool_int(
            pid_repro.get("is_bent_stem") if pid_repro.get("is_bent_stem") is not None else p_pet_info.get("is_bent_stem"),
            default=0,
        )

        pet_mm = _to_float(p_pet_info.get("petiole_length_mm"))
        bld_mm = _to_float(p_pet_info.get("blade_length_mm"))
        p2b = _to_float(p_pet_info.get("petiole_to_blade_ratio"))
        pet_px = _to_float(p_pet_info.get("petiole_length_px"))
        bld_px = _to_float(p_pet_info.get("blade_length_px"))
        pet_w_px = _to_float(p_pet_info.get("petiole_width_px"))
        pet_w_mm = _to_float(p_pet_info.get("petiole_width_mm"))

        pid_leaves = leaves_by_pid.get(pid) or []
        basal_pid_leaves = [
            x for x in pid_leaves
            if str(x.get("leaf_zone") or "basal").strip().lower() != "cauline"
        ]
        summary_pid_leaves = basal_pid_leaves if basal_pid_leaves else pid_leaves
        if summary_pid_leaves:
            if pet_mm is None:
                pet_mm = _p75_leaf_metric(summary_pid_leaves, "petiole_length_mm")
            if pet_w_px is None:
                pet_w_px = _p75_leaf_metric(summary_pid_leaves, "petiole_width_px")
            if pet_w_mm is None:
                pet_w_mm = _p75_leaf_metric(summary_pid_leaves, "petiole_width_mm")
            if bld_mm is None:
                bld_mm = _p75_leaf_metric(summary_pid_leaves, "blade_length_mm")
            if p2b is None:
                p2b = _p75_leaf_metric(summary_pid_leaves, "petiole_to_blade_ratio")

        plant_rows.append((
            cat, pid,
            _to_float(rc_xy[0]) if rc_xy else None,
            _to_float(rc_xy[1]) if rc_xy else None,
            _to_float(cr_xy[0]) if cr_xy else None,
            _to_float(cr_xy[1]) if cr_xy else None,
            p_height_px, p_height_mm, p_height_cm, is_bent,
            pet_px, pet_mm, pet_w_px, pet_w_mm, bld_px, bld_mm,
            p2b, scape_ratio, None, None,
        ))

        if pid_repro or pid_pheno:
            ar_inv = _to_float(pid_repro.get("AR_inv"))
            n_fr_ing = _to_int(pid_pheno.get("n_fruiting"))
            ph_state = _to_str(pid_pheno.get("voucher_phenological_state"))
            is_ver_anth = _to_bool_int(pid_pheno.get("is_verified_anthesis"))
            if is_ver_anth is None and ph_state:
                low_st = ph_state.lower()
                is_ver_anth = 1 if ("anthesis" in low_st and "pre" not in low_st) else 0
            corymb_est_raw = (
                pid_repro.get("corymb_head_count_est")
                if pid_repro.get("corymb_head_count_est") is not None
                else pid_pheno.get("corymb_head_count_est")
            )
            radiate_frac_raw = (
                pid_repro.get("radiate_head_fraction")
                if pid_repro.get("radiate_head_fraction") is not None
                else pid_pheno.get("radiate_head_fraction")
            )
            repro_pheno_rows.append((
                cat, pid,
                _to_int(pid_repro.get("capitula_count")),
                _to_float(pid_repro.get("involucre_height_px")),
                _to_float(pid_repro.get("involucre_width_px")),
                _to_float(pid_repro.get("involucre_height_mm")),
                _to_float(pid_repro.get("involucre_width_mm")),
                ar_inv,
                _to_int(pid_pheno.get("total_capitula"), default=_to_int(pid_repro.get("capitula_count"))),
                _to_int(pid_pheno.get("n_anthesis")),
                _to_int(pid_pheno.get("n_bud")),
                n_fr_ing,
                _to_int(pid_pheno.get("n_senescent")),
                ph_state,
                is_ver_anth,
                _to_float(pid_pheno.get("reproductive_stage_confidence")),
                _to_float(pid_pheno.get("dominant_pappus_ratio")),
                _to_int(corymb_est_raw, default=0) or 0,
                _to_float(radiate_frac_raw),
                _to_float(pid_repro.get("phyllary_tip_darkening_delta")),
                _to_float(pid_repro.get("phyllary_tip_anthocyanin_delta")),
                _to_float(pid_repro.get("involucre_basal_cie_l_delta")),
                _to_float(pid_repro.get("peduncle_cie_l_sheet_delta")),
                _to_float(pid_repro.get("calyculus_flare_index")),
            ))
    return plant_rows, repro_pheno_rows


def _build_step02_leaf_rows(
    cat: str,
    extracted_records: Sequence[Mapping[str, Any]],
    ppm_val: Optional[float],
    coords_map: Mapping[int, Sequence[Sequence[float]]],
    out_base: Path,
    workspace_root: Path,
    resolve_contour_fn: Callable[..., Optional[Path]],
    read_contour_fn: Callable[[Path], Optional[List[List[float]]]],
) -> List[Tuple[Any, ...]]:
    """Resolve real 120-vertex contours and build canonical `extracted_leaves` insert tuples."""
    leaf_rows: List[Tuple[Any, ...]] = []
    for r in extracted_records:
        leaf_id = _to_int(r.get("leaf_id"), default=1) or 1
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0

        xy_pts: Optional[List[List[float]]] = None
        if r.get("contour_xy_json"):
            try:
                parsed = json.loads(str(r["contour_xy_json"]))
                if isinstance(parsed, list) and len(parsed) >= 3:
                    xy_pts = [[round(float(pt[0]), 4), round(float(pt[1]), 4)] for pt in parsed]
            except (ValueError, TypeError, KeyError, IndexError):
                xy_pts = None
        if xy_pts is None and leaf_id in coords_map:
            xy_pts = [[round(float(pt[0]), 4), round(float(pt[1]), 4)] for pt in coords_map[leaf_id]]
        if xy_pts is None:
            c_file = resolve_contour_fn(_to_str(r.get("contour_path")), output_dir=out_base, workspace_root=workspace_root)
            if c_file is not None:
                xy_pts = read_contour_fn(c_file)
        if not xy_pts:
            # Skip unbacked leaf rows lacking real contour vertices (Rule 2: Zero Silent Synthetic Data)
            continue

        contour_json = json.dumps(xy_pts, separators=(",", ":"))
        crop_p = _to_str(r.get("crop_path"))
        shape_id = _to_str(r.get("shape_id")) or f"{cat}_p{pid}_leaf{leaf_id}"
        l_ppm = _to_float(r.get("pixels_per_mm"))
        if l_ppm is None:
            l_ppm = ppm_val
        l_src_raw = (_to_str(r.get("scale_source")) or _to_str(r.get("cf_source")) or "").strip().lower()
        l_is_mp = l_src_raw in _UNCALIBRATED_SCALE_SOURCES
        l_cal = _to_bool_int(r.get("is_calibrated_mm"))
        if l_is_mp:
            l_cal = 0
        elif l_cal is None:
            l_cal = 1 if (l_ppm is not None and l_ppm > 0) else 0

        area_px_val = _to_float(r.get("area_px"))
        area_mm2_val = (
            _to_float(r.get("area_mm2"))
            if r.get("area_mm2") is not None
            else (round(float(area_px_val) / (float(l_ppm) ** 2), 2) if (l_cal == 1 and l_ppm is not None and l_ppm > 0 and area_px_val is not None) else None)
        )
        scale_src_val = _to_str(r.get("scale_source")) or (
            "mp_estimated" if (l_cal == 0 and l_ppm is not None and l_ppm > 0)
            else ("ruler_calibrated" if l_cal == 1 else "uncalibrated")
        )

        leaf_rows.append((
            cat, pid, leaf_id, shape_id,
            _to_float(r.get("solidity")), _to_float(r.get("circularity")),
            _to_float(r.get("extent")), _to_float(r.get("aspect_ratio")),
            area_px_val, area_mm2_val, _to_float(r.get("midrib_angle_deg")),
            _to_float(r.get("petiole_length_px")), _to_float(r.get("petiole_length_mm")),
            _to_float(r.get("petiole_width_px")), _to_float(r.get("petiole_width_mm")),
            _to_float(r.get("blade_length_px")), _to_float(r.get("blade_length_mm")),
            _to_float(r.get("blade_width_px")), _to_float(r.get("blade_width_mm")),
            _to_float(r.get("petiole_to_blade_ratio")), l_ppm, l_cal, scale_src_val,
            _to_str(r.get("mask_path")), _to_bool_int(r.get("mask_is_upright"), default=1),
            crop_p, _to_str(r.get("contour_path")), _to_str(r.get("ect_path")),
            len(xy_pts), contour_json, _to_str(r.get("leaf_zone")) or "basal",
            _to_float(r.get("leaf_zone_confidence")), _to_float(r.get("dist_from_crown_norm")),
            _to_float(r.get("cy_rel_plant_bbox")), _to_float(r.get("archetype_score")),
            _to_float(r.get("bilateral_dice")), _to_float(r.get("si_a")),
            _to_bool_int(r.get("lm3_gates_pass"), default=1),
            max(1, int(_to_int(r.get("n_lamina_components")) or 1)),
            int(_to_int(r.get("detector_class_id"), default=0) or 0),
            _to_bool_int(r.get("is_crop_clipped"), default=0),
            _to_float(r.get("hole_area_frac")) if r.get("hole_area_frac") is not None else 0.0,
            _to_float(r.get("enclosed_void_frac")) if r.get("enclosed_void_frac") is not None else 0.0,
            _to_bool_int(r.get("bilateral_lobe_symmetry"), default=0),
            _to_bool_int(r.get("efa_qc_pass"), default=1),
            _to_str(r.get("qc_flag_reason")) or "PASS",
            _to_float(r.get("blade_cie_l_norm")), _to_float(r.get("blade_cie_l_sd")),
            _to_float(r.get("blade_cie_l_sheet_delta")), _to_float(r.get("blade_cie_a_delta")),
            _to_float(r.get("blade_cie_b_delta")),
            _to_int(r.get("bbox_x1")), _to_int(r.get("bbox_y1")),
            _to_int(r.get("bbox_x2")), _to_int(r.get("bbox_y2")),
            _to_float(r.get("base_angle_deg")), _to_float(r.get("apex_angle_deg")),
            _to_float(r.get("lamina_curvature")), _to_float(r.get("lamina_trace_length_px")),
            _to_float(r.get("basal_cordate_index")), _to_float(r.get("margin_serration_index")),
            _to_float(r.get("margin_tooth_acuteness")), _to_float(r.get("terminal_lobe_area_frac")),
            _to_float(r.get("blade_cie_l_q85_sheet_delta")), _to_float(r.get("basal_zone_cie_l_delta")),
            _to_bool_int(r.get("is_abaxial_pressed_leaf"), default=0),
            _to_float(r.get("abaxial_cie_l_sheet_delta")),
        ))
    return leaf_rows


__all__ = (
    "_EXTRACTED_LEAVES_INSERT_COLS",
    "_p75_leaf_metric",
    "_resolve_voucher_calibration_and_qc",
    "_build_step02_plant_and_repro_rows",
    "_build_step02_leaf_rows",
)
