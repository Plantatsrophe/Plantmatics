#!/usr/bin/env python3
"""Step 02 leaf orientation, morphology/void/margin QC, contour/ECT export, and crown-refreshed capitula."""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

from scripts.vision.capitulum_phenology_classifier import classify_single_capitulum
from scripts.vision.leaf_geometry_utils import (
    compute_lm3_shape_metrics,
    estimate_rosette_crown_point,
    export_standardized_contour,
    save_leaf_mask_png,
)
from scripts.vision.lm3_data_loader import DEFAULT_HERBARIUM_DPI, compute_area_mm2, resolve_px_per_cm
from scripts.vision.lm3_inference_engines import (
    DetectedInstance,
    compute_capitulum_metrics,
    compute_leaf_ect_descriptor,
    evaluate_leaf_damage_and_fold_qc,
    orient_and_score_leaf_lm3,
)
from scripts.vision.step02_capitulum_microchars import extract_single_capitulum_microchars
from scripts.vision.step02_leaf_traits_utils import (
    extract_packera_leaf_margin_and_base_traits,
    measure_enclosed_interlobe_void_frac,
)
from scripts.vision.step02_photometry_corymb_utils import compute_calibrated_blade_photometry

_PHOTOMETRY_KEYS: Tuple[str, ...] = (
    "blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta", "blade_cie_a_delta", "blade_cie_b_delta",
)


def _compute_plant_bbox_rel_y_and_crown_hints(
    instances: Sequence[DetectedInstance],
    inst_geom_map: Dict[int, Tuple[float, float, float, float, float]],
    capitulum_detections: Sequence[Any],
    root_detections: Sequence[Tuple[int, int, int, int]],
    upper_cauline_pts: Sequence[Tuple[float, float]],
    rosette_centroids_by_pid: Dict[int, Optional[Tuple[float, float]]],
    basal_stats_by_pid: Dict[int, Tuple[float, float, float, float]],
    sheet_h: int,
    sheet_w: int,
) -> Tuple[Dict[int, float], Dict[int, Optional[Tuple[float, float]]]]:
    """Compute per-leaf relative y within its plant bounding box and per-plant crown tie-break hints."""
    pid_ymins: Dict[int, List[float]] = {}
    pid_ymaxs: Dict[int, List[float]] = {}
    pid_xs: Dict[int, List[float]] = {}
    for inst_b in instances:
        p_id = int(inst_b.plant_individual_id)
        pid_ymins.setdefault(p_id, []).append(float(inst_b.bbox[0]))
        pid_ymaxs.setdefault(p_id, []).append(float(inst_b.bbox[2]))
        pid_xs.setdefault(p_id, []).extend([float(inst_b.bbox[1]), float(inst_b.bbox[3])])

    if len(pid_ymins) == 1:
        sole_pid = next(iter(pid_ymins.keys()))
        for box_b in [cb for cb, _, _, _ in capitulum_detections] + list(root_detections):
            pid_ymins[sole_pid].append(float(box_b[0]))
            pid_ymaxs[sole_pid].append(float(box_b[2]))
        for _uc_x, uc_y in upper_cauline_pts:
            pid_ymins[sole_pid].append(float(uc_y))
            pid_ymaxs[sole_pid].append(float(uc_y))
    elif len(pid_ymins) > 1:
        pid_centers_xy = {
            k: (float(np.median(pid_xs[k])), float(np.median(pid_ymins[k] + pid_ymaxs[k]))) for k in pid_ymins.keys()
        }
        max_assoc_d = 0.35 * max(sheet_h, sheet_w)
        for box_b in [cb for cb, _, _, _ in capitulum_detections] + list(root_detections):
            cx_b, cy_b = 0.5 * (float(box_b[1]) + float(box_b[3])), 0.5 * (float(box_b[0]) + float(box_b[2]))
            best_k = min(pid_centers_xy.keys(), key=lambda k: math.hypot(cx_b - pid_centers_xy[k][0], cy_b - pid_centers_xy[k][1]))
            if math.hypot(cx_b - pid_centers_xy[best_k][0], cy_b - pid_centers_xy[best_k][1]) <= max_assoc_d:
                pid_ymins[best_k].append(float(box_b[0]))
                pid_ymaxs[best_k].append(float(box_b[2]))

    inst_cy_rel_plant_bbox: Dict[int, float] = {}
    for inst_b in instances:
        p_id = int(inst_b.plant_individual_id)
        ymin_p, ymax_p = float(min(pid_ymins.get(p_id, [0.0]))), float(max(pid_ymaxs.get(p_id, [float(sheet_h)])))
        inst_cy_rel_plant_bbox[inst_b.leaf_id] = max(0.0, min(1.0, (inst_geom_map[inst_b.leaf_id][1] - ymin_p) / max(1.0, ymax_p - ymin_p)))

    crown_hints_sheet: Dict[int, Optional[Tuple[float, float]]] = {}
    for pid_h, rc_h in rosette_centroids_by_pid.items():
        best_root_pt: Optional[Tuple[float, float]] = None
        if root_detections:
            ref_xy = rc_h if rc_h is not None else (basal_stats_by_pid[pid_h][0], basal_stats_by_pid[pid_h][1])
            r_cands = [(0.5 * float(rb[1] + rb[3]), 0.5 * float(rb[0] + rb[2])) for rb in root_detections]
            nearest_r = min(r_cands, key=lambda rxy: math.hypot(rxy[0] - ref_xy[0], rxy[1] - ref_xy[1]))
            if len(rosette_centroids_by_pid) == 1 or math.hypot(nearest_r[0] - ref_xy[0], nearest_r[1] - ref_xy[1]) <= 0.35 * max(sheet_h, sheet_w):
                best_root_pt = nearest_r
        crown_hints_sheet[pid_h] = best_root_pt if best_root_pt is not None else rc_h
    return inst_cy_rel_plant_bbox, crown_hints_sheet


def _apply_blade_photometry_to_inst(
    inst: DetectedInstance, img_bgr: np.ndarray, sheet_paper_l_ref: float, pixels_per_mm: Optional[float], is_calibrated_mm: bool,
) -> None:
    """Compute wool-preserving 7-key CIELAB lamina photometry and assign attributes onto `inst`."""
    try:
        phot = compute_calibrated_blade_photometry(
            img_bgr, inst.mask, mask_origin=inst.mask_origin, hole_mask=inst.hole_mask,
            petiole_mask=inst.petiole_mask, paper_l_ref=sheet_paper_l_ref,
            pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
            keypoints_31=inst.keypoints_31 if isinstance(inst.keypoints_31, dict) else None,
        )
    except Exception:
        phot = {}
    for k in _PHOTOMETRY_KEYS:
        setattr(inst, k, float(phot[k]) if (isinstance(phot, dict) and phot.get(k) is not None) else np.nan)


def _refresh_cap_item_with_crown(
    c_it: Dict[str, Any], crown_xy: Tuple[float, float], pixels_per_mm: Optional[float], is_calibrated_mm: bool,
    pheno_params: Dict[str, Any], sheet_bgr: Optional[np.ndarray] = None,
    paper_lab_ref: Optional[Tuple[float, float, float]] = None,
) -> None:
    """Re-evaluate involucre metrics, micro-characters, and phenology using crown-anchored peduncle direction."""
    ped_dir = (c_it["center_xy"][0] - float(crown_xy[0]), c_it["center_xy"][1] - float(crown_xy[1]))
    m_ref = compute_capitulum_metrics(
        c_it["mask"], c_it["bbox"], pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
        crop_bgr=c_it["crop_bgr"], class_id=c_it["class_id"], peduncle_dir_xy=ped_dir, is_calibrated_mm=is_calibrated_mm,
    )
    if m_ref is not None:
        m_ref["score"], m_ref["class_id"] = c_it["score"], c_it["class_id"]
    c_it["metrics"] = m_ref
    c_bgr_r = c_it.get("crop_bgr")
    if c_bgr_r is not None and c_bgr_r.size > 0 and c_bgr_r.shape[0] >= 5 and c_bgr_r.shape[1] >= 5:
        c_it["head_pheno"] = classify_single_capitulum(
            cv2.cvtColor(c_bgr_r, cv2.COLOR_BGR2RGB),
            min_pappus_luminance=pheno_params["min_pappus_lum"], max_pappus_saturation=pheno_params["max_pappus_sat"],
            min_pappus_area_ratio=pheno_params["min_pappus_area"], min_floret_yellow_saturation=pheno_params["min_floret_sat"],
            min_floret_yellow_ratio=pheno_params["min_floret_rat"], class_id=int(c_it["class_id"]), peduncle_dir_xy=ped_dir,
        )
        c_it["microchars"] = extract_single_capitulum_microchars(
            crop_bgr=c_bgr_r, sheet_bgr=sheet_bgr, bbox=c_it["bbox"], mask=c_it["mask"],
            class_id=int(c_it["class_id"]), peduncle_dir_xy=ped_dir,
            pixels_per_mm=pixels_per_mm if is_calibrated_mm else None, is_calibrated_mm=is_calibrated_mm,
            paper_lab_ref=paper_lab_ref,
        )


refresh_cap_item_with_crown = _refresh_cap_item_with_crown


def export_leaf_mask_contour_and_ect(
    pipe: Any,
    catalog_number: str,
    inst: DetectedInstance,
    p_base_pt: Optional[Tuple[int, int]],
    p_apex_pt: Optional[Tuple[int, int]],
    angle_cw: Optional[float],
    persisted_mask_u8: np.ndarray,
    ori_pts_map: Dict[str, Tuple[float, float]],
    is_upright_mask: bool,
    ori_hole_by_leaf: Dict[int, Optional[np.ndarray]],
    contours_by_leaf: Dict[int, List[List[float]]],
) -> None:
    """Export upright PNG mask, 120-vertex tip-up contour CSV, and 32x360 ECT .npy for one leaf."""
    inst.mask_path = save_leaf_mask_png(persisted_mask_u8, catalog_number, inst.leaf_id, pipe.output_dir)
    if is_upright_mask:
        c_res = export_standardized_contour(
            persisted_mask_u8, catalog_number, inst.leaf_id, pipe.output_dir,
            petiole_attachment_pt=(ori_pts_map.get("lamina_base") or ori_pts_map.get("base_center")),
            p_apex_pt=(ori_pts_map.get("lamina_tip") or ori_pts_map.get("apex_center")),
            return_coords=True, already_oriented=True,
        )
    else:
        c_res = export_standardized_contour(
            persisted_mask_u8, catalog_number, inst.leaf_id, pipe.output_dir,
            petiole_attachment_pt=p_base_pt, p_apex_pt=p_apex_pt,
            hole_mask=inst.hole_mask, return_coords=True, keypoints_31=inst.keypoints_31,
        )
    if isinstance(c_res, tuple):
        inst.contour_path, c_xy = c_res
        if c_xy:
            contours_by_leaf[int(inst.leaf_id)] = c_xy
    else:
        inst.contour_path = c_res

    if getattr(pipe, "ect_dir", None) is not None and getattr(pipe, "ect_enabled", True):
        ect_out = pipe.ect_dir / f"{catalog_number}_leaf{inst.leaf_id}_ect.npy"
        ori_hm = ori_hole_by_leaf.get(int(inst.leaf_id))
        if not is_upright_mask and ori_hm is None:
            ori_hm = inst.hole_mask
        if (
            not bool(getattr(pipe, "ect_include_holes", True)) and ori_hm is not None
            and isinstance(ori_hm, np.ndarray) and ori_hm.shape[:2] == persisted_mask_u8.shape[:2]
        ):
            ect_mask_in = ((persisted_mask_u8 > 0) | (ori_hm > 0)).astype(np.uint8) * 255
        else:
            ect_mask_in = persisted_mask_u8
        try:
            ect_mat = compute_leaf_ect_descriptor(
                ect_mask_in, num_dirs=getattr(pipe, "ect_num_dirs", 360), num_thresh=getattr(pipe, "ect_num_thetas", 32),
                output_path=ect_out, hole_mask=(None if is_upright_mask else inst.hole_mask),
                petiole_mask=(None if is_upright_mask else inst.petiole_mask), angle_cw=(None if is_upright_mask else angle_cw),
                include_petiole=bool(getattr(pipe, "ect_include_petiole", False)),
                include_holes=bool(getattr(pipe, "ect_include_holes", True)), already_oriented=is_upright_mask,
            )
        except TypeError:
            ect_mat = compute_leaf_ect_descriptor(ect_mask_in, output_path=ect_out)
        if ect_mat is not None:
            inst.ect_path = str(ect_out)


def _resolve_leaf_base_apex(
    inst: DetectedInstance, crown_local: Optional[Tuple[float, float]]
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Resolve local-window `(p_base_use, p_apex_use)` from pose landmarks, petiole anchor, or contour fallback."""
    p_base_use: Optional[Tuple[int, int]] = None
    p_apex_use: Optional[Tuple[int, int]] = None
    if (
        isinstance(inst.keypoints_31, dict) and "lamina_base" in inst.keypoints_31 and "lamina_tip" in inst.keypoints_31
        and len(inst.keypoints_31["lamina_base"]) >= 2 and len(inst.keypoints_31["lamina_tip"]) >= 2
    ):
        lb_i, lt_i = inst.keypoints_31["lamina_base"], inst.keypoints_31["lamina_tip"]
        p_base_use = (int(round(float(lb_i[0]))), int(round(float(lb_i[1]))))
        p_apex_use = (int(round(float(lt_i[0]))), int(round(float(lt_i[1]))))
    elif isinstance(inst.trim_meta, dict):
        if inst.trim_meta.get("p_base") is not None:
            p_base_use = inst.trim_meta["p_base"]
        if inst.trim_meta.get("p_apex") is not None:
            p_apex_use = inst.trim_meta["p_apex"]

    if p_base_use is None or p_apex_use is None:
        m_u8 = (inst.mask > 0).astype(np.uint8) * 255
        cnts_fb, _ = cv2.findContours(m_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if cnts_fb:
            pts_fb = max(cnts_fb, key=cv2.contourArea).reshape(-1, 2)
            if p_base_use is None:
                if inst.petiole_mask is not None and int(np.count_nonzero(inst.petiole_mask)) > 0:
                    pys, pxs = np.where(inst.petiole_mask > 0)
                    bi = int(np.argmin(np.hypot(pts_fb[:, 0] - float(np.mean(pxs)), pts_fb[:, 1] - float(np.mean(pys)))))
                elif crown_local is not None:
                    bi = int(np.argmin(np.hypot(pts_fb[:, 0] - crown_local[0], pts_fb[:, 1] - crown_local[1])))
                else:
                    bi = int(np.argmax(pts_fb[:, 1]))
                p_base_use = (int(pts_fb[bi, 0]), int(pts_fb[bi, 1]))
            ai = int(np.argmax(np.hypot(pts_fb[:, 0] - p_base_use[0], pts_fb[:, 1] - p_base_use[1])))
            p_apex_use = p_apex_use or (int(pts_fb[ai, 0]), int(pts_fb[ai, 1]))
        else:
            p_base_use, p_apex_use = p_base_use or (0, 0), p_apex_use or (0, 0)
    return p_base_use, p_apex_use


def orient_and_score_voucher_leaves(
    pipe: Any,
    img_bgr: np.ndarray,
    instances: Sequence[DetectedInstance],
    inst_geom_map: Dict[int, Tuple[float, float, float, float, float]],
    capitulum_detections: Sequence[Any],
    root_detections: Sequence[Tuple[int, int, int, int]],
    upper_cauline_pts: Sequence[Tuple[float, float]],
    rosette_centroids_by_pid: Dict[int, Optional[Tuple[float, float]]],
    basal_stats_by_pid: Dict[int, Tuple[float, float, float, float]],
    sheet_h: int,
    sheet_w: int,
    is_high_dpi: bool,
    pixels_per_mm: Optional[float],
    is_calibrated_mm: bool,
    cf_source: str,
    sheet_paper_l_ref: float,
    spec_mask_attempted: bool,
    spec_mask_sheet: Optional[np.ndarray],
) -> Dict[str, Any]:
    """Run Step 6a read-only LeafMachine3 orientation, void/margin/damage QC, pose packaging, and crown solving."""
    px_per_cm_sheet = (
        resolve_px_per_cm(pixels_per_mm=pixels_per_mm, default_dpi=DEFAULT_HERBARIUM_DPI, sheet_shape=(sheet_h, sheet_w))
        if is_calibrated_mm else None
    )
    inst_cy_rel_plant_bbox, crown_hints_sheet = _compute_plant_bbox_rel_y_and_crown_hints(
        instances, inst_geom_map, capitulum_detections, root_detections, upper_cauline_pts,
        rosette_centroids_by_pid, basal_stats_by_pid, sheet_h, sheet_w,
    )
    raw_px_by_leaf_id: Dict[int, float] = {}
    for inst_pre in instances:
        rpx = float(np.count_nonzero(inst_pre.mask))
        raw_px_by_leaf_id[inst_pre.leaf_id] = rpx
        inst_pre.area_mm2 = compute_area_mm2(rpx, pixels_per_mm=pixels_per_mm, is_calibrated_mm=is_calibrated_mm)

    ori_hole_by_leaf: Dict[int, Optional[np.ndarray]] = {}
    candidate_extractions: List[Any] = []
    seg_partial_allow = bool(getattr(pipe.config.segmentation, "segment_partial_leaves", False))

    for inst in instances:
        if int(np.count_nonzero(inst.mask)) <= 0:
            continue
        _apply_blade_photometry_to_inst(inst, img_bgr, sheet_paper_l_ref, pixels_per_mm, is_calibrated_mm)
        ox_l, oy_l = inst.mask_origin
        pid_hint_s = crown_hints_sheet.get(int(inst.plant_individual_id))
        crown_local = (float(pid_hint_s[0]) - float(ox_l), float(pid_hint_s[1]) - float(oy_l)) if pid_hint_s is not None else None
        p_base_use, p_apex_use = _resolve_leaf_base_apex(inst, crown_local)

        inst.cf_source = cf_source
        lm3_bm_raw = orient_and_score_leaf_lm3(
            inst.mask, keypoints_31=inst.keypoints_31, p_base=p_base_use, p_apex=p_apex_use,
            offset_xy=(0, 0), hole_mask=inst.hole_mask, petiole_mask=inst.petiole_mask, crown_hint_xy=crown_local,
        )
        angle_cw_val = float(lm3_bm_raw.get("angle_cw", 0.0)) if lm3_bm_raw.get("success") else None
        if angle_cw_val is not None:
            inst.midrib_angle_deg = float((((-90.0 - float(angle_cw_val)) + 180.0) % 360.0) - 180.0)
        else:
            inst.midrib_angle_deg = float(math.degrees(math.atan2(float(p_apex_use[1] - p_base_use[1]), float(p_apex_use[0] - p_base_use[0]))))

        lm3_morph = compute_lm3_shape_metrics(inst.mask, angle_cw=angle_cw_val)
        inst.solidity, inst.circularity, inst.extent = float(lm3_morph["solidity"]), float(lm3_morph["circularity"]), float(lm3_morph["extent"])
        bm_lw, bm_bl_px, bm_bw_px = lm3_bm_raw.get("aspect_ratio", np.nan), lm3_bm_raw.get("blade_length_px"), lm3_bm_raw.get("blade_width_px")
        if lm3_bm_raw.get("success") and bm_lw is not None and np.isfinite(float(bm_lw)) and float(bm_lw) > 0.0:
            inst.aspect_ratio = float(bm_lw)
            if bm_bl_px is not None and np.isfinite(float(bm_bl_px)) and float(bm_bl_px) > 0.0:
                inst.blade_length_px = float(bm_bl_px)
            if bm_bw_px is not None and np.isfinite(float(bm_bw_px)) and float(bm_bw_px) > 0.0:
                inst.blade_width_px = float(bm_bw_px)
        else:
            inst.aspect_ratio = float(lm3_morph["aspect_ratio"])
        inst.archetype_score = float(lm3_bm_raw.get("archetype_score", 0.0))
        inst.bilateral_dice = float(lm3_bm_raw.get("bilateral_dice", 1.0))  # type: ignore[attr-defined]
        inst.si_a = float(lm3_bm_raw.get("si_a", 0.0))  # type: ignore[attr-defined]
        inst.lm3_gates_pass = bool(lm3_bm_raw.get("gates_pass", True))  # type: ignore[attr-defined]
        lm3_tm_sub = lm3_bm_raw.get("trim_meta") if isinstance(lm3_bm_raw.get("trim_meta"), dict) else {}
        if not isinstance(inst.trim_meta, dict):
            inst.trim_meta = {}
        for nm_k in ("base_angle_deg", "apex_angle_deg", "lamina_curvature", "lamina_trace_length_px"):
            nm_raw = lm3_bm_raw.get(nm_k, lm3_tm_sub.get(nm_k, inst.trim_meta.get(nm_k, np.nan)))
            try:
                nm_f = float(nm_raw) if nm_raw is not None and np.isfinite(float(nm_raw)) else np.nan
            except (TypeError, ValueError):
                nm_f = np.nan
            inst.trim_meta[nm_k] = nm_f
            setattr(inst, nm_k, nm_f)

        ori_m = lm3_bm_raw.get("oriented_mask")
        if isinstance(ori_m, np.ndarray) and ori_m.size > 0 and int(np.count_nonzero(ori_m)) > 0:
            persisted_mask_u8 = (ori_m > 0).astype(np.uint8) * 255
            ori_pts_map: Dict[str, Tuple[float, float]] = dict(lm3_bm_raw.get("oriented_points") or {})
            is_upright_mask = True
        else:
            persisted_mask_u8 = (inst.mask > 0).astype(np.uint8) * 255
            if bool(getattr(pipe, "subtract_holes", True)) and isinstance(inst.hole_mask, np.ndarray) and inst.hole_mask.shape[:2] == persisted_mask_u8.shape[:2]:
                persisted_mask_u8 = np.where((persisted_mask_u8 > 0) & (inst.hole_mask == 0), 255, 0).astype(np.uint8)
            ori_pts_map, is_upright_mask = {}, False

        n_cc_lbls, _, cc_stats, _ = cv2.connectedComponentsWithStats((inst.mask > 0).astype(np.uint8), connectivity=8)
        n_sig_comps = sum(1 for cid in range(1, n_cc_lbls) if int(cc_stats[cid, cv2.CC_STAT_AREA]) >= 32) if n_cc_lbls > 1 else 1
        inst.n_lamina_components = max(1, n_sig_comps)  # type: ignore[attr-defined]

        ori_hole_m = lm3_bm_raw.get("oriented_hole_mask") if is_upright_mask else inst.hole_mask
        lam_fg_px = float(max(np.count_nonzero(persisted_mask_u8 > 0), 1))
        hole_px = 0.0
        if isinstance(ori_hole_m, np.ndarray) and ori_hole_m.shape[:2] == persisted_mask_u8.shape[:2]:
            cnts_h_hull, _ = cv2.findContours(persisted_mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if cnts_h_hull:
                hull_mask_h = np.zeros_like(persisted_mask_u8)
                cv2.fillPoly(hull_mask_h, [cv2.convexHull(np.vstack([c.reshape(-1, 1, 2) for c in cnts_h_hull]))], 255)
                hole_px = float(np.count_nonzero((ori_hole_m > 0) & (hull_mask_h > 0)))
            else:
                hole_px = float(np.count_nonzero(ori_hole_m > 0))
        hole_area_frac = hole_px / max(lam_fg_px + hole_px, 1.0)
        enclosed_void_frac = float(measure_enclosed_interlobe_void_frac(persisted_mask_u8, hole_mask=ori_hole_m if isinstance(ori_hole_m, np.ndarray) else None))
        packera_traits = extract_packera_leaf_margin_and_base_traits(persisted_mask_u8, oriented_points=ori_pts_map if ori_pts_map else None)
        ori_hole_by_leaf[int(inst.leaf_id)] = ori_hole_m if isinstance(ori_hole_m, np.ndarray) else None
        trim_d = inst.trim_meta if isinstance(inst.trim_meta, dict) else {}
        is_clipped = bool(trim_d.get("is_crop_clipped", False))
        if not is_clipped and inst.mask.shape[:2] == (sheet_h, sheet_w):
            sh_allow = max(8, int(round(0.08 * max(inst.bbox[2] - inst.bbox[0], inst.bbox[3] - inst.bbox[1], 1))))
            is_clipped = bool(
                int(np.count_nonzero(inst.mask[0, :] > 0)) > sh_allow or int(np.count_nonzero(inst.mask[-1, :] > 0)) > sh_allow
                or int(np.count_nonzero(inst.mask[:, 0] > 0)) > sh_allow or int(np.count_nonzero(inst.mask[:, -1] > 0)) > sh_allow
            )
        qc_eval = evaluate_leaf_damage_and_fold_qc(
            oriented_mask_u8=persisted_mask_u8, solidity=float(inst.solidity),
            bilateral_dice=float(inst.bilateral_dice), si_a=float(inst.si_a),  # type: ignore[attr-defined]
            lm3_gates_pass=bool(inst.lm3_gates_pass), n_lamina_components=int(inst.n_lamina_components),  # type: ignore[attr-defined]
            detector_class_id=int(inst.class_id), segment_partial_leaves=seg_partial_allow,
            is_crop_clipped=is_clipped, hole_area_frac=hole_area_frac, enclosed_void_frac=enclosed_void_frac,
            oriented_points=ori_pts_map, keypoints_31=inst.keypoints_31 if isinstance(inst.keypoints_31, dict) else None,
            is_geometric_fallback_pose=bool(trim_d.get("is_geometric_fallback_pose", inst.keypoints_31 is None)),
        )
        inst.detector_class_id, inst.is_crop_clipped = int(inst.class_id), bool(is_clipped)  # type: ignore[attr-defined]
        inst.hole_area_frac, inst.enclosed_void_frac = round(float(hole_area_frac), 4), round(float(enclosed_void_frac), 4)  # type: ignore[attr-defined]
        for tk in ("basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac"):
            setattr(inst, tk, packera_traits.get(tk, np.nan))
        inst.bilateral_lobe_symmetry = bool(qc_eval["bilateral_lobe_symmetry"])  # type: ignore[attr-defined]
        inst.efa_qc_pass = bool(qc_eval["efa_qc_pass"])  # type: ignore[attr-defined]
        inst.qc_flag_reason = str(qc_eval["qc_flag_reason"])  # type: ignore[attr-defined]
        inst.min_margin_rms_norm = float(qc_eval["min_margin_rms_norm"])  # type: ignore[attr-defined]
        candidate_extractions.append((inst, p_base_use, p_apex_use, angle_cw_val, persisted_mask_u8, ori_pts_map, is_upright_mask))

    pose_inputs: List[Dict[str, Any]] = []
    for inst, pb, pa, _ang_cw, _pmask, _opts, _is_up in candidate_extractions:
        ox, oy = inst.mask_origin
        pb_use, pa_use = pb or (0, 0), pa or (0, 0)
        pb_sheet, pa_sheet = (float(pb_use[0] + ox), float(pb_use[1] + oy)), (float(pa_use[0] + ox), float(pa_use[1] + oy))
        trim_info = inst.trim_meta or {}
        geod_px = float(trim_info.get("petiole_length_px", 0.0) or trim_info.get("geodesic_petiole_skeleton_px", 0.0) or trim_info.get("trim_distance_px", 0.0) or 0.0)
        pose_inputs.append({
            "leaf_id": inst.leaf_id, "plant_individual_id": inst.plant_individual_id,
            "p_base": pb_sheet, "p_apex": pa_sheet, "p_base_sheet": pb_sheet, "p_apex_sheet": pa_sheet,
            "mask_origin": (0.0, 0.0), "bbox": inst.bbox, "mask": inst.mask,
            "area_px": raw_px_by_leaf_id.get(inst.leaf_id, float(np.count_nonzero(inst.mask))),
            "area_mm2": inst.area_mm2, "blade_length_px": inst.blade_length_px, "blade_width_px": inst.blade_width_px,
            "petiole_width_px": inst.petiole_width_px, "petiole_width_mm": inst.petiole_width_mm,
            "aspect_ratio": inst.aspect_ratio, "solidity": inst.solidity,
            "geodesic_petiole_px": geod_px, "geodesic_petiole_skeleton_px": geod_px,
            "cy_rel_plant_bbox": inst_cy_rel_plant_bbox.get(inst.leaf_id, np.nan),
            "keypoints_31": inst.keypoints_31, "petiole_mask": inst.petiole_mask, "local_mask_origin": (ox, oy),
        })

    crowns_by_pid = estimate_rosette_crown_point(
        pose_inputs, root_boxes=list(root_detections), rosette_centroids=rosette_centroids_by_pid, return_by_plant=True,
    )
    if (
        not spec_mask_attempted and spec_mask_sheet is None and len(pose_inputs) >= 2 and is_high_dpi
        and getattr(pipe, "specimen_segmenter", None) is not None
        and getattr(pipe.specimen_segmenter, "model", None) is not None
    ):
        spec_mask_attempted = True
        try:
            spec_mask_sheet = pipe.specimen_segmenter.predict_mask(img_bgr, return_sheet_size=True)
        except Exception:
            spec_mask_sheet = None

    return {
        "px_per_cm_sheet": px_per_cm_sheet, "inst_cy_rel_plant_bbox": inst_cy_rel_plant_bbox,
        "raw_px_by_leaf_id": raw_px_by_leaf_id, "ori_hole_by_leaf": ori_hole_by_leaf,
        "candidate_extractions": candidate_extractions, "pose_inputs": pose_inputs,
        "crowns_by_pid": crowns_by_pid, "spec_mask_attempted": spec_mask_attempted, "spec_mask_sheet": spec_mask_sheet,
    }
