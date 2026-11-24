#!/usr/bin/env python3
"""Step 02 capitula partitioning, petiole/scape/leaf-zone readout, leaf assembly, and `process_voucher_impl`."""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import cv2
import numpy as np

from scripts.vision.capitulum_phenology_classifier import classify_voucher_phenology
from scripts.vision.leaf_geometry_utils import approximate_petiole_and_scape_metrics
from scripts.vision.lm3_inference_engines import DetectedInstance, _VoucherFinalizeBundle
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical
from scripts.vision.step02_leaf_traits_utils import evaluate_rosette_abaxial_leaf_flags
from scripts.vision.step02_pipe_leaf_readout_a import (
    export_leaf_mask_contour_and_ect,
    orient_and_score_voucher_leaves,
    refresh_cap_item_with_crown,
)
from scripts.vision.step02_pipe_repro_clustering import (
    build_clustered_instances_and_basal_stats,
    route_capitula_phenology,
    summarize_repro_subset,
)
from scripts.vision.step02_pipe_stage1_candidates import (
    _resolve_scale_source,
    load_voucher_image_or_failure,
    prepare_voucher_scale_stage1_and_specimen,
)
from scripts.vision.step02_pipe_stage2_inference import run_stage2_leaf_and_pose_inference


def _partition_capitula_by_plant(
    catalog_number: str,
    unique_leaf_pids: List[int],
    pose_inputs: List[Dict[str, Any]],
    raw_cap_items: List[Dict[str, Any]],
    crowns_by_pid: Any,
    spec_mask_sheet: Optional[np.ndarray],
    sheet_h: int,
    sheet_w: int,
    pixels_per_mm: Optional[float],
    is_calibrated_mm: bool,
    pheno_params: Dict[str, Any],
) -> Tuple[Dict[int, List[Tuple[int, int, int, int]]], Dict[int, List[Dict[str, Any]]], Dict[str, Any], Dict[str, Any]]:
    """Assign capitula across `plant_individual_id`, refresh crown-aligned metrics, and update voucher summaries."""
    caps_by_pid_boxes: Dict[int, List[Tuple[int, int, int, int]]] = {pid: [] for pid in unique_leaf_pids}
    cap_items_by_pid: Dict[int, List[Dict[str, Any]]] = {pid: [] for pid in unique_leaf_pids}

    if len(unique_leaf_pids) <= 1:
        only_pid = unique_leaf_pids[0] if unique_leaf_pids else 0
        cr_only = crowns_by_pid.get(only_pid) if isinstance(crowns_by_pid, dict) else None
        for c_it in raw_cap_items:
            if cr_only is not None:
                refresh_cap_item_with_crown(c_it, cr_only, pixels_per_mm, is_calibrated_mm, pheno_params)
            caps_by_pid_boxes.setdefault(only_pid, []).append(c_it["bbox"])
            cap_items_by_pid.setdefault(only_pid, []).append(c_it)
        repro_summary, voucher_extra_heads = summarize_repro_subset(
            catalog_number, raw_cap_items, pixels_per_mm, is_calibrated_mm, sheet_h, sheet_w, pid_val=only_pid
        )
        pheno_pid_v = int(only_pid)
    else:
        cc_labels_sheet: Optional[np.ndarray] = None
        pid_cc_sets: Dict[int, Set[int]] = {pid: set() for pid in unique_leaf_pids}
        if isinstance(spec_mask_sheet, np.ndarray) and spec_mask_sheet.shape[:2] == (sheet_h, sheet_w) and np.any(spec_mask_sheet > 0):
            _, cc_labels_sheet = cv2.connectedComponents((spec_mask_sheet > 0).astype(np.uint8), connectivity=8)
            for p_in in pose_inputs:
                px_i = int(np.clip(round(p_in["p_base"][0]), 0, sheet_w - 1))
                py_i = int(np.clip(round(p_in["p_base"][1]), 0, sheet_h - 1))
                lbl_p = int(cc_labels_sheet[py_i, px_i])
                if lbl_p > 0:
                    pid_cc_sets[int(p_in["plant_individual_id"])].add(lbl_p)
        max_link_dist = 0.55 * float(max(sheet_h, sheet_w))
        for c_it in raw_cap_items:
            cx_c, cy_c = c_it["center_xy"]
            assigned_pid: Optional[int] = None
            if cc_labels_sheet is not None:
                lbl_c = int(cc_labels_sheet[int(np.clip(round(cy_c), 0, sheet_h - 1)), int(np.clip(round(cx_c), 0, sheet_w - 1))])
                if lbl_c > 0:
                    matching_pids = [pid for pid, s_lbls in pid_cc_sets.items() if lbl_c in s_lbls]
                    if len(matching_pids) == 1:
                        assigned_pid = matching_pids[0]
            if assigned_pid is None:
                best_d = float("inf")
                for pid in unique_leaf_pids:
                    pid_pts = [p["p_base"] for p in pose_inputs if int(p["plant_individual_id"]) == pid] + [
                        p["p_apex"] for p in pose_inputs if int(p["plant_individual_id"]) == pid
                    ]
                    cr_p = crowns_by_pid.get(pid) if isinstance(crowns_by_pid, dict) else None
                    if cr_p is not None:
                        pid_pts.append((float(cr_p[0]), float(cr_p[1])))
                    if not pid_pts:
                        continue
                    min_leaf_d = min(math.hypot(cx_c - px_p, cy_c - py_p) for px_p, py_p in pid_pts)
                    eff_d = min_leaf_d * (1.35 if (cr_p is not None and cy_c > float(cr_p[1]) + 40.0) else 1.0)
                    if eff_d < best_d:
                        best_d, assigned_pid = eff_d, pid
                if best_d > max_link_dist:
                    assigned_pid = None
            if assigned_pid is not None:
                cr_ass = crowns_by_pid.get(assigned_pid) if isinstance(crowns_by_pid, dict) else None
                if cr_ass is not None:
                    refresh_cap_item_with_crown(c_it, cr_ass, pixels_per_mm, is_calibrated_mm, pheno_params)
                caps_by_pid_boxes[assigned_pid].append(c_it["bbox"])
                cap_items_by_pid[assigned_pid].append(c_it)
        repro_summary, voucher_extra_heads = summarize_repro_subset(
            catalog_number, raw_cap_items, pixels_per_mm, is_calibrated_mm, sheet_h, sheet_w, pid_val=0
        )
        pheno_pid_v = 0

    pheno_rec = dict(classify_voucher_phenology(
        catalog_number,
        [it["head_pheno"] for it in raw_cap_items if it.get("head_pheno") is not None],
        min_anthesis_heads_for_flowering=pheno_params["min_anth_heads"],
    ))
    pheno_rec["plant_individual_id"] = pheno_pid_v
    pheno_rec["total_capitula"] = int(pheno_rec.get("total_capitula", 0)) + int(voucher_extra_heads)
    pheno_rec["corymb_head_count_est"] = int(voucher_extra_heads)
    pheno_rec.setdefault("radiate_head_fraction", np.nan)
    repro_summary["radiate_head_fraction"] = pheno_rec.get("radiate_head_fraction", np.nan)
    return caps_by_pid_boxes, cap_items_by_pid, repro_summary, pheno_rec


def _screen_basal_rosette_half_folds(
    candidate_extractions: Sequence[Any], leaf_zones_by_id: Dict[int, Dict[str, Any]], raw_px_by_leaf_id: Dict[int, float],
) -> None:
    """Screen basal rosettes (>= 3 basal leaves) for half-folded blades (`LONGITUDINAL_MARGIN_CREASE`)."""
    basal_by_pid_qc: Dict[int, List[DetectedInstance]] = defaultdict(list)
    for inst_q, *_ in candidate_extractions:
        lz_q = str((leaf_zones_by_id.get(int(inst_q.leaf_id)) or {}).get("leaf_zone", "basal") or "basal").lower()
        if lz_q == "basal" and np.isfinite(float(inst_q.aspect_ratio)) and float(inst_q.aspect_ratio) > 0.0:
            basal_by_pid_qc[int(inst_q.plant_individual_id)].append(inst_q)
    for _pid_q, b_list in basal_by_pid_qc.items():
        if len(b_list) < 3:
            continue
        clean_pool = [b for b in b_list if bool(getattr(b, "efa_qc_pass", True))]
        ref_pool = clean_pool if len(clean_pool) >= 2 else b_list
        med_ar = float(np.median([float(b.aspect_ratio) for b in ref_pool]))
        med_area = float(np.median([raw_px_by_leaf_id.get(b.leaf_id, 1.0) for b in ref_pool]))
        for b_inst in b_list:
            b_rms, b_bd = float(getattr(b_inst, "min_margin_rms_norm", np.nan)), float(getattr(b_inst, "bilateral_dice", 1.0))
            b_area = float(raw_px_by_leaf_id.get(b_inst.leaf_id, med_area))
            if (
                float(b_inst.aspect_ratio) > 1.85 * max(med_ar, 0.5)
                and b_area < 0.65 * max(med_area, 1.0)
                and ((np.isfinite(b_rms) and b_rms < 0.028) or b_bd < 0.78)
            ):
                toks = [t for t in str(getattr(b_inst, "qc_flag_reason", "PASS") or "PASS").split(";") if t and t != "PASS"]
                if "LONGITUDINAL_MARGIN_CREASE" not in toks:
                    toks.append("LONGITUDINAL_MARGIN_CREASE")
                b_inst.efa_qc_pass = False  # type: ignore[attr-defined]
                b_inst.qc_flag_reason = ";".join(toks)  # type: ignore[attr-defined]


def _build_extracted_leaf_record(
    catalog_number: str, inst: DetectedInstance, is_upright_mask: bool, raw_px_by_leaf_id: Dict[int, float],
    leaf_zones_by_id: Dict[int, Dict[str, Any]], per_leaf_petiole: Dict[int, Dict[str, Any]],
    repro_by_pid_out: Dict[int, Dict[str, Any]], repro_summary: Dict[str, Any],
    inst_cy_rel_plant_bbox: Dict[int, float], rosette_bbox_fallback_used: bool,
) -> Dict[str, Any]:
    """Assemble the canonical Step 02 leaf manifest dictionary for one `DetectedInstance`."""
    leaf_pet = per_leaf_petiole.get(int(inst.leaf_id), {})
    inst_area_px_val = raw_px_by_leaf_id.get(inst.leaf_id, np.nan)
    lz_info = leaf_zones_by_id.get(int(inst.leaf_id), {}) if isinstance(leaf_zones_by_id, dict) else {}
    leaf_zone_str = str(lz_info.get("leaf_zone", "basal") or "basal")
    leaf_zone_conf, dist_crown_norm = lz_info.get("leaf_zone_confidence", np.nan), lz_info.get("dist_from_crown_norm", np.nan)
    pet_len_px, pet_len_mm = leaf_pet.get("petiole_length_px", np.nan), leaf_pet.get("petiole_length_mm", np.nan)

    if inst.blade_length_px is not None and np.isfinite(float(inst.blade_length_px)) and float(inst.blade_length_px) > 0.0:
        blade_len_px = round(float(inst.blade_length_px), 4)
        blade_len_mm = (
            round(float(inst.blade_length_px) / float(inst.pixels_per_mm), 3)
            if (inst.is_calibrated_mm and inst.pixels_per_mm is not None and float(inst.pixels_per_mm) > 0.0)
            else np.nan
        )
    else:
        blade_len_px, blade_len_mm = leaf_pet.get("blade_length_px", np.nan), leaf_pet.get("blade_length_mm", np.nan)

    if (
        pet_len_px is not None and np.isfinite(float(pet_len_px)) and float(pet_len_px) >= 0.0
        and blade_len_px is not None and np.isfinite(float(blade_len_px)) and float(blade_len_px) > 0.0
    ):
        pet_ratio = round(float(pet_len_px) / float(blade_len_px), 4)
    else:
        pet_ratio = leaf_pet.get("petiole_to_blade_ratio", np.nan)

    pid_repro_rec = repro_by_pid_out.get(int(inst.plant_individual_id), repro_summary)
    pet_w_px = round(float(inst.petiole_width_px), 4) if (inst.petiole_width_px is not None and np.isfinite(inst.petiole_width_px)) else leaf_pet.get("petiole_width_px", np.nan)
    pet_w_mm = round(float(inst.petiole_width_mm), 4) if (inst.petiole_width_mm is not None and np.isfinite(inst.petiole_width_mm)) else leaf_pet.get("petiole_width_mm", np.nan)
    if inst.blade_width_px is not None and np.isfinite(inst.blade_width_px) and inst.blade_width_px > 0.0:
        blade_w_px = round(float(inst.blade_width_px), 4)
    elif (
        blade_len_px is not None and np.isfinite(float(blade_len_px)) and float(blade_len_px) > 0.0
        and np.isfinite(float(inst.aspect_ratio)) and float(inst.aspect_ratio) > 1e-6
    ):
        blade_w_px = round(float(blade_len_px) / float(inst.aspect_ratio), 4)
    else:
        blade_w_px = np.nan
    blade_w_mm = (
        round(float(blade_w_px) / float(inst.pixels_per_mm), 3)
        if (np.isfinite(blade_w_px) and inst.is_calibrated_mm and inst.pixels_per_mm is not None and float(inst.pixels_per_mm) > 0.0)
        else np.nan
    )
    bb_y1, bb_x1, bb_y2, bb_x2 = [int(v) for v in inst.bbox]
    tm_out = inst.trim_meta if isinstance(inst.trim_meta, dict) else {}

    return {
        "catalogNumber": catalog_number, "plant_individual_id": inst.plant_individual_id, "leaf_id": inst.leaf_id,
        "bbox_x1": bb_x1, "bbox_y1": bb_y1, "bbox_x2": bb_x2, "bbox_y2": bb_y2,
        "solidity": round(float(inst.solidity), 4),
        "circularity": round(float(inst.circularity), 4) if np.isfinite(inst.circularity) else np.nan,
        "extent": round(float(inst.extent), 4) if np.isfinite(inst.extent) else np.nan,
        "aspect_ratio": round(float(inst.aspect_ratio), 4) if np.isfinite(inst.aspect_ratio) else np.nan,
        "area_px": round(float(inst_area_px_val), 1) if np.isfinite(inst_area_px_val) else np.nan,
        "area_mm2": round(float(inst.area_mm2), 2) if np.isfinite(inst.area_mm2) else np.nan,
        "archetype_score": round(float(inst.archetype_score), 4),
        "bilateral_dice": round(float(getattr(inst, "bilateral_dice", 1.0)), 4),
        "si_a": round(float(getattr(inst, "si_a", 0.0)), 4), "lm3_gates_pass": bool(getattr(inst, "lm3_gates_pass", True)),
        "n_lamina_components": int(getattr(inst, "n_lamina_components", 1)),
        "detector_class_id": int(getattr(inst, "detector_class_id", inst.class_id)),
        "is_crop_clipped": bool(getattr(inst, "is_crop_clipped", False)),
        "hole_area_frac": round(float(getattr(inst, "hole_area_frac", 0.0)), 4),
        "enclosed_void_frac": round(float(getattr(inst, "enclosed_void_frac", 0.0)), 4),
        "bilateral_lobe_symmetry": bool(getattr(inst, "bilateral_lobe_symmetry", False)),
        "efa_qc_pass": bool(getattr(inst, "efa_qc_pass", True)), "qc_flag_reason": str(getattr(inst, "qc_flag_reason", "PASS") or "PASS"),
        "base_angle_deg": tm_out.get("base_angle_deg", getattr(inst, "base_angle_deg", np.nan)),
        "apex_angle_deg": tm_out.get("apex_angle_deg", getattr(inst, "apex_angle_deg", np.nan)),
        "lamina_curvature": tm_out.get("lamina_curvature", getattr(inst, "lamina_curvature", np.nan)),
        "lamina_trace_length_px": tm_out.get("lamina_trace_length_px", getattr(inst, "lamina_trace_length_px", np.nan)),
        "basal_cordate_index": getattr(inst, "basal_cordate_index", np.nan), "margin_serration_index": getattr(inst, "margin_serration_index", np.nan),
        "margin_tooth_acuteness": getattr(inst, "margin_tooth_acuteness", np.nan), "terminal_lobe_area_frac": getattr(inst, "terminal_lobe_area_frac", np.nan),
        "blade_cie_l_norm": getattr(inst, "blade_cie_l_norm", np.nan), "blade_cie_l_sd": getattr(inst, "blade_cie_l_sd", np.nan),
        "blade_cie_l_sheet_delta": getattr(inst, "blade_cie_l_sheet_delta", np.nan), "blade_cie_l_q85_sheet_delta": getattr(inst, "blade_cie_l_q85_sheet_delta", np.nan),
        "basal_zone_cie_l_delta": getattr(inst, "basal_zone_cie_l_delta", np.nan),
        "blade_cie_a_delta": getattr(inst, "blade_cie_a_delta", np.nan), "blade_cie_b_delta": getattr(inst, "blade_cie_b_delta", np.nan),
        "is_abaxial_pressed_leaf": False, "abaxial_cie_l_sheet_delta": np.nan,
        "midrib_angle_deg": round(float(inst.midrib_angle_deg), 2),
        "pixels_per_mm": round(float(inst.pixels_per_mm), 4) if inst.pixels_per_mm else np.nan,
        "is_calibrated_mm": bool(inst.is_calibrated_mm), "cf_source": inst.cf_source,
        "scale_source": _resolve_scale_source(inst.pixels_per_mm, bool(inst.is_calibrated_mm), inst.cf_source),
        "mask_path": inst.mask_path, "mask_is_upright": bool(is_upright_mask),
        "contour_path": inst.contour_path, "ect_path": inst.ect_path,
        "petiole_length_px": pet_len_px, "petiole_width_px": pet_w_px,
        "blade_length_px": blade_len_px, "blade_width_px": blade_w_px,
        "petiole_length_mm": pet_len_mm, "petiole_width_mm": pet_w_mm,
        "blade_length_mm": blade_len_mm, "blade_width_mm": blade_w_mm,
        "petiole_to_blade_ratio": pet_ratio,
        "capitula_count": pid_repro_rec["capitula_count"], "log1p_capitulum_count": pid_repro_rec.get("log1p_capitulum_count", np.nan),
        "involucre_height_px": pid_repro_rec["involucre_height_px"], "involucre_width_px": pid_repro_rec["involucre_width_px"],
        "involucre_height_mm": pid_repro_rec["involucre_height_mm"], "involucre_width_mm": pid_repro_rec["involucre_width_mm"],
        "AR_inv": pid_repro_rec["AR_inv"], "leaf_zone": leaf_zone_str,
        "leaf_zone_confidence": round(float(leaf_zone_conf), 4) if leaf_zone_conf is not None and np.isfinite(leaf_zone_conf) else np.nan,
        "dist_from_crown_norm": round(float(dist_crown_norm), 4) if dist_crown_norm is not None and np.isfinite(dist_crown_norm) else np.nan,
        "cy_rel_plant_bbox": round(float(inst_cy_rel_plant_bbox.get(inst.leaf_id, np.nan)), 4),
        "rosette_bbox_fallback_used": bool(rosette_bbox_fallback_used),
    }


def process_voucher_impl(
    pipe: Any, catalog_number: str, image_path: Path, taxon: Optional[str] = None,
    preloaded_img: Optional[np.ndarray] = None, precomputed_s1_letterbox: Optional[np.ndarray] = None,
    precomputed_s1_detections: Optional[List[Any]] = None, precomputed_s2_payload: Optional[Tuple[Any, ...]] = None,
    precomputed_single_stage_detections: Optional[Union[List[Any], Tuple[List[Any], bool]]] = None,
    _return_bundle: bool = False, precomputed_spec_mask: Optional[np.ndarray] = None,
) -> Any:
    """Process a single herbarium voucher specimen sheet using LeafMachine3 ONNX engines."""
    if taxon is None:
        taxon = pipe.vouchers_taxa.get(catalog_number, None)

    def _finish(
        ext_recs: List[Dict[str, Any]], fail_rec: Optional[Dict[str, Any]],
        repro_rec: Optional[Dict[str, Any]] = None, pheno_rec: Optional[Dict[str, Any]] = None,
        crowns_out: Optional[Dict[int, Tuple[float, float]]] = None, centroids_out: Optional[Dict[int, Tuple[float, float]]] = None,
        petiole_out: Optional[Dict[int, Dict[str, Any]]] = None, contours_out: Optional[Dict[int, List[List[float]]]] = None,
        repro_pid_out: Optional[Dict[int, Dict[str, Any]]] = None, pheno_pid_out: Optional[Dict[int, Dict[str, Any]]] = None,
    ) -> Any:
        bundle = _VoucherFinalizeBundle(
            catalog_number=catalog_number, image_path=str(image_path),
            extracted_records=ext_recs, failure_record=fail_rec, repro_record=repro_rec, pheno_record=pheno_rec,
            repro_by_pid=repro_pid_out or {}, pheno_by_pid=pheno_pid_out or {},
            crowns_by_pid=crowns_out or {}, rosette_centroids_by_pid=centroids_out or {},
            by_plant_petiole=petiole_out or {}, contour_coords_by_leaf=contours_out or {},
        )
        if _return_bundle:
            return bundle
        pipe._apply_finalize_bundle(bundle)
        return bundle.extracted_records, bundle.failure_record

    img_bgr, image_path, load_fail = load_voucher_image_or_failure(catalog_number, image_path, preloaded_img=preloaded_img)
    if load_fail is not None or img_bgr is None:
        return _finish([], load_fail, pheno_rec=classify_voucher_phenology(catalog_number, []))
    if pipe.force:
        pipe._purge_voucher_disk_artifacts(catalog_number)

    s1_ctx = prepare_voucher_scale_stage1_and_specimen(
        pipe, catalog_number, img_bgr, precomputed_s1_letterbox=precomputed_s1_letterbox,
        precomputed_s1_detections=precomputed_s1_detections, precomputed_spec_mask=precomputed_spec_mask,
    )
    sheet_h, sheet_w, is_high_dpi = s1_ctx["sheet_h"], s1_ctx["sheet_w"], s1_ctx["is_high_dpi"]
    pixels_per_mm, cf_source = s1_ctx["pixels_per_mm"], s1_ctx["cf_source"]

    s2_ctx = run_stage2_leaf_and_pose_inference(
        pipe, catalog_number, image_path, img_bgr, sheet_h, sheet_w, is_high_dpi,
        pixels_per_mm, cf_source, s1_ctx["stage1_detections"], s1_ctx["spec_mask_sheet"],
        precomputed_s2_payload=precomputed_s2_payload,
        precomputed_single_stage_detections=precomputed_single_stage_detections,
    )
    if s2_ctx["early_exit"] is not None:
        fail_r, emp_rep, emp_phe = s2_ctx["early_exit"]
        return _finish([], fail_r, repro_rec=emp_rep, pheno_rec=emp_phe)

    leaf_detections, capitulum_detections = s2_ctx["leaf_detections"], s2_ctx["capitulum_detections"]
    root_detections, upper_cauline_pts = s2_ctx["root_detections"], s2_ctx["upper_cauline_pts"]
    rosette_bbox_fallback_used, is_calibrated_mm = s2_ctx["rosette_bbox_fallback_used"], s2_ctx["is_calibrated_mm"]

    repro_ctx = route_capitula_phenology(
        pipe, catalog_number, img_bgr, sheet_h, sheet_w, capitulum_detections, pixels_per_mm, is_calibrated_mm
    )
    if not leaf_detections:
        failure_record = {
            "catalogNumber": catalog_number, "image_path": str(image_path), "failure_reason": "no_leaves_detected",
            "details": f"Model produced 0 candidate leaf segmentations ({repro_ctx['capitula_count']} capitula detected)",
            "pixels_per_mm": round(float(pixels_per_mm), 4) if is_calibrated_mm else np.nan,
            "is_calibrated_mm": bool(is_calibrated_mm),
        }
        return _finish([], failure_record, repro_rec=repro_ctx["repro_summary"], pheno_rec=repro_ctx["pheno_rec"])

    instances, inst_geom_map, rosette_centroids_by_pid, basal_stats_by_pid = build_clustered_instances_and_basal_stats(
        pipe, catalog_number, leaf_detections, root_detections, sheet_h, sheet_w,
        pixels_per_mm, is_calibrated_mm, s1_ctx["spec_mask_sheet"],
    )
    ori_ctx = orient_and_score_voucher_leaves(
        pipe, img_bgr, instances, inst_geom_map, capitulum_detections, root_detections,
        upper_cauline_pts, rosette_centroids_by_pid, basal_stats_by_pid, sheet_h, sheet_w, is_high_dpi,
        pixels_per_mm, is_calibrated_mm, cf_source, s1_ctx["sheet_paper_l_ref"],
        s1_ctx["spec_mask_attempted"], s1_ctx["spec_mask_sheet"],
    )
    pose_inputs, crowns_by_pid, spec_mask_sheet = ori_ctx["pose_inputs"], ori_ctx["crowns_by_pid"], ori_ctx["spec_mask_sheet"]
    unique_leaf_pids = sorted({int(p["plant_individual_id"]) for p in pose_inputs})
    pheno_params = repro_ctx["pheno_params"]

    caps_by_pid_boxes, cap_items_by_pid, repro_summary, pheno_rec = _partition_capitula_by_plant(
        catalog_number, unique_leaf_pids, pose_inputs, repro_ctx["raw_cap_items"],
        crowns_by_pid, spec_mask_sheet, sheet_h, sheet_w, pixels_per_mm, is_calibrated_mm, pheno_params,
    )
    cap_boxes_all = [b for b, _, _, _ in capitulum_detections]
    petiole_scape_res = approximate_petiole_and_scape_metrics(
        basal_leaves=pose_inputs, p_crown=crowns_by_pid, px_per_cm=ori_ctx["px_per_cm_sheet"],
        is_calibrated_mm=is_calibrated_mm, upper_cauline_detections=upper_cauline_pts,
        root_boxes=root_detections, capitula_detections=cap_boxes_all, capitula_by_pid=caps_by_pid_boxes,
    )
    repro_summary["plant_height_cm"] = petiole_scape_res.get("plant_height_cm", np.nan)
    repro_summary["scape_to_basal_leaf_ratio"] = petiole_scape_res.get("scape_to_basal_leaf_ratio", np.nan)
    repro_summary["is_bent_stem"] = bool(petiole_scape_res.get("is_bent_stem", False))

    leaf_zones_by_id = classify_leaf_zone_anatomical(
        pose_inputs, crowns_by_pid=crowns_by_pid if isinstance(crowns_by_pid, dict) else None,
        rosette_centroids_by_pid=rosette_centroids_by_pid, basal_stats_by_pid=basal_stats_by_pid,
        capitula_detections=cap_boxes_all, root_boxes=root_detections,
        is_bent_stem=bool(repro_summary.get("is_bent_stem", False)),
        sheet_shape=(sheet_h, sheet_w), specimen_mask=spec_mask_sheet,
    )
    has_cauline_zone = False
    for p_in in pose_inputs:
        lz = str((leaf_zones_by_id.get(int(p_in["leaf_id"])) or {}).get("leaf_zone", "basal") or "basal").lower()
        p_in["leaf_zone"] = lz
        if lz == "cauline":
            has_cauline_zone = True
    if has_cauline_zone:
        petiole_scape_res = approximate_petiole_and_scape_metrics(
            basal_leaves=pose_inputs, p_crown=crowns_by_pid, px_per_cm=ori_ctx["px_per_cm_sheet"],
            is_calibrated_mm=is_calibrated_mm, upper_cauline_detections=upper_cauline_pts,
            root_boxes=root_detections, capitula_detections=cap_boxes_all, capitula_by_pid=caps_by_pid_boxes,
        )
        if np.isfinite(petiole_scape_res.get("scape_to_basal_leaf_ratio", np.nan)):
            repro_summary["scape_to_basal_leaf_ratio"] = petiole_scape_res.get("scape_to_basal_leaf_ratio", np.nan)

    by_plant_petiole = petiole_scape_res.get("by_plant", {}) or {}
    per_leaf_petiole = {int(r["leaf_id"]): r for r in (petiole_scape_res.get("leaf_metrics") or []) if isinstance(r, dict) and r.get("leaf_id") is not None}
    centroids_out = {int(k): (float(v[0]), float(v[1])) for k, v in rosette_centroids_by_pid.items() if v is not None}
    crowns_out = {int(k): (float(v[0]), float(v[1])) for k, v in crowns_by_pid.items() if isinstance(v, (tuple, list)) and len(v) >= 2} if isinstance(crowns_by_pid, dict) else {}
    petiole_out = {int(k): dict(v) for k, v in by_plant_petiole.items() if isinstance(v, dict)} if isinstance(by_plant_petiole, dict) else {}

    repro_by_pid_out: Dict[int, Dict[str, Any]] = {}
    pheno_by_pid_out: Dict[int, Dict[str, Any]] = {}
    for pid_r in unique_leaf_pids:
        p_caps = cap_items_by_pid.get(pid_r) or []
        r_sub, sub_extra_heads = summarize_repro_subset(
            catalog_number, p_caps, pixels_per_mm, is_calibrated_mm, sheet_h, sheet_w, pid_val=pid_r
        )
        pet_sub = petiole_out.get(pid_r, {})
        if len(unique_leaf_pids) <= 1:
            r_sub["plant_height_cm"], r_sub["scape_to_basal_leaf_ratio"], r_sub["is_bent_stem"] = (
                repro_summary["plant_height_cm"], repro_summary["scape_to_basal_leaf_ratio"], repro_summary["is_bent_stem"]
            )
        else:
            r_sub["plant_height_cm"] = pet_sub.get("plant_height_cm", np.nan)
            r_sub["scape_to_basal_leaf_ratio"] = pet_sub.get("scape_to_basal_leaf_ratio", np.nan)
            r_sub["is_bent_stem"] = bool(pet_sub.get("is_bent_stem", False))
        ph_sub = dict(classify_voucher_phenology(
            catalog_number, [it["head_pheno"] for it in p_caps if it.get("head_pheno") is not None],
            min_anthesis_heads_for_flowering=pheno_params["min_anth_heads"],
        ))
        ph_sub["plant_individual_id"] = int(pid_r)
        ph_sub["total_capitula"] = int(ph_sub.get("total_capitula", 0)) + int(sub_extra_heads)
        ph_sub["corymb_head_count_est"] = int(sub_extra_heads)
        ph_sub.setdefault("radiate_head_fraction", np.nan)
        r_sub["radiate_head_fraction"] = ph_sub.get("radiate_head_fraction", np.nan)
        repro_by_pid_out[pid_r], pheno_by_pid_out[pid_r] = r_sub, ph_sub

    _screen_basal_rosette_half_folds(ori_ctx["candidate_extractions"], leaf_zones_by_id, ori_ctx["raw_px_by_leaf_id"])
    contours_by_leaf: Dict[int, List[List[float]]] = {}
    extracted_records: List[Dict[str, Any]] = []
    for inst, pb_pt, pa_pt, angle_cw, persisted_mask_u8, ori_pts_map, is_upright_mask in ori_ctx["candidate_extractions"]:
        export_leaf_mask_contour_and_ect(
            pipe, catalog_number, inst, pb_pt, pa_pt, angle_cw, persisted_mask_u8,
            ori_pts_map, is_upright_mask, ori_ctx["ori_hole_by_leaf"], contours_by_leaf,
        )
        extracted_records.append(_build_extracted_leaf_record(
            catalog_number, inst, is_upright_mask, ori_ctx["raw_px_by_leaf_id"],
            leaf_zones_by_id, per_leaf_petiole, repro_by_pid_out, repro_summary,
            ori_ctx["inst_cy_rel_plant_bbox"], rosette_bbox_fallback_used,
        ))
    if extracted_records:
        evaluate_rosette_abaxial_leaf_flags(extracted_records)

    failure_record = None if extracted_records else {
        "catalogNumber": catalog_number, "image_path": str(image_path), "failure_reason": "no_leaves_detected",
        "details": "LeafMachine3 produced 0 non-empty leaf masks",
        "pixels_per_mm": round(float(pixels_per_mm), 4) if is_calibrated_mm else np.nan,
        "is_calibrated_mm": bool(is_calibrated_mm),
    }
    return _finish(
        extracted_records, failure_record, repro_rec=repro_summary, pheno_rec=pheno_rec,
        crowns_out=crowns_out, centroids_out=centroids_out, petiole_out=petiole_out,
        contours_out=contours_by_leaf, repro_pid_out=repro_by_pid_out, pheno_pid_out=pheno_by_pid_out,
    )
