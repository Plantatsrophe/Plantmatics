#!/usr/bin/env python3
"""Step 02 reproductive capitula phenology routing, DBSCAN clustering, and basal rosette stats."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

from scripts.vision.capitulum_phenology_classifier import classify_single_capitulum, classify_voucher_phenology
from scripts.vision.leaf_geometry_utils import cluster_plant_individuals
from scripts.vision.lm3_inference_engines import DetectedInstance, compute_capitulum_metrics
from scripts.vision.step02_capitulum_microchars import (
    MICROCHAR_REPRO_KEYS,
    extract_single_capitulum_microchars,
    summarize_capitulum_microchars,
)
from scripts.vision.step02_photometry_corymb_utils import (
    estimate_corymb_cluster_head_count,
    estimate_sheet_paper_cie_l_ref_and_ab,
)


def summarize_repro_subset(
    catalog_number: str,
    cap_items_sub: List[Dict[str, Any]],
    pixels_per_mm: Optional[float],
    is_calibrated_mm: bool,
    sheet_h: int,
    sheet_w: int,
    pid_val: int = 0,
) -> Tuple[Dict[str, Any], int]:
    """Summarize involucre geometry and corymb head count for a voucher or plant individual."""
    valid_metrics_sub = [
        it["metrics"] for it in cap_items_sub
        if it.get("metrics") is not None and int(it.get("class_id", -1)) in (2, 4, 6)
    ]
    cluster_items_sub = [it for it in cap_items_sub if int(it.get("class_id", -1)) in (3, 5)]
    single_items_sub = [it for it in cap_items_sub if int(it.get("class_id", -1)) in (2, 4, 6)]
    extra_cluster_heads, _ = estimate_corymb_cluster_head_count(
        cluster_items_sub, single_items_sub,
        pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
        sheet_shape_hw=(sheet_h, sheet_w),
    )
    c_cnt = len(valid_metrics_sub) + int(extra_cluster_heads)

    def _med_key(key: str) -> float:
        vals = [float(c[key]) for c in valid_metrics_sub if c.get(key) is not None and np.isfinite(float(c[key]))]
        return float(np.median(vals)) if vals else np.nan

    if valid_metrics_sub:
        mh_px, mw_px = _med_key("involucre_height_px"), _med_key("involucre_width_px")
        mh_mm, mw_mm = _med_key("involucre_height_mm"), _med_key("involucre_width_mm")
        m_ar = _med_key("AR_inv")
    else:
        mh_px = mw_px = mh_mm = mw_mm = m_ar = np.nan

    mc_summary = summarize_capitulum_microchars(cap_items_sub)
    rec = {
        "catalogNumber": catalog_number, "plant_individual_id": int(pid_val),
        "capitula_count": int(c_cnt), "corymb_head_count_est": int(extra_cluster_heads),
        "log1p_capitulum_count": round(math.log1p(int(c_cnt)), 5) if int(c_cnt) >= 1 else np.nan,
        "involucre_height_px": round(mh_px, 2) if not np.isnan(mh_px) else np.nan,
        "involucre_width_px": round(mw_px, 2) if not np.isnan(mw_px) else np.nan,
        "involucre_height_mm": round(mh_mm, 3) if not np.isnan(mh_mm) else np.nan,
        "involucre_width_mm": round(mw_mm, 3) if not np.isnan(mw_mm) else np.nan,
        "AR_inv": round(m_ar, 4) if not np.isnan(m_ar) else np.nan,
        "plant_height_cm": np.nan, "scape_to_basal_leaf_ratio": np.nan,
        "is_bent_stem": False, "radiate_head_fraction": np.nan,
        **mc_summary,
    }
    return rec, int(extra_cluster_heads)


def route_capitula_phenology(
    pipe: Any,
    catalog_number: str,
    img_bgr: np.ndarray,
    sheet_h: int,
    sheet_w: int,
    capitulum_detections: Sequence[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]],
    pixels_per_mm: Optional[float],
    is_calibrated_mm: bool,
) -> Dict[str, Any]:
    """Route reproductive detections (classes 2..6) to the phenology classifier and involucre metrics."""
    pheno_cfg = getattr(getattr(pipe.config, "thresholds", None), "capitulum_phenology", None)
    pheno_params = {
        "min_pappus_lum": float(getattr(pheno_cfg, "min_pappus_luminance", 75.0) if pheno_cfg is not None else 75.0),
        "max_pappus_sat": float(getattr(pheno_cfg, "max_pappus_saturation", 0.18) if pheno_cfg is not None else 0.18),
        "min_pappus_area": float(getattr(pheno_cfg, "min_pappus_area_ratio", 0.18) if pheno_cfg is not None else 0.18),
        "min_floret_sat": float(getattr(pheno_cfg, "min_floret_yellow_saturation", 0.32) if pheno_cfg is not None else 0.32),
        "min_floret_rat": float(getattr(pheno_cfg, "min_floret_yellow_ratio", 0.12) if pheno_cfg is not None else 0.12),
        "min_anth_heads": int(getattr(pheno_cfg, "min_anthesis_heads_for_flowering", 1) or 1) if pheno_cfg is not None else 1,
    }
    valid_capitula: List[Dict[str, Any]] = []
    head_pheno_results: List[Dict[str, Any]] = []
    raw_cap_items: List[Dict[str, Any]] = []
    paper_lab_ref = estimate_sheet_paper_cie_l_ref_and_ab(img_bgr) if capitulum_detections else None

    for cap_idx, (bbox, mask, score, class_id) in enumerate(capitulum_detections, 1):
        ymin, xmin, ymax, xmax = [int(v) for v in bbox]
        y0, y1, x0, x1 = max(0, ymin), min(sheet_h, max(0, ymax)), max(0, xmin), min(sheet_w, max(0, xmax))
        crop_bgr = img_bgr[y0:y1, x0:x1]
        local_cap_mask: Optional[np.ndarray] = None
        if mask is not None and isinstance(mask, np.ndarray) and mask.size > 0:
            local_cap_mask = mask[y0:y1, x0:x1] if mask.shape[:2] == (sheet_h, sheet_w) else mask

        head_pheno_item: Optional[Dict[str, Any]] = None
        if crop_bgr.size > 0:
            if bool(getattr(pipe, "save_debug_crops", False)):
                pipe.capitula_patches_dir.mkdir(parents=True, exist_ok=True)
                cap_patch_path = pipe.capitula_patches_dir / f"{catalog_number}_cap{cap_idx}.jpg"
                try:
                    cv2.imwrite(str(cap_patch_path), crop_bgr, [int(cv2.IMWRITE_JPEG_QUALITY), 85])
                except TypeError:
                    cv2.imwrite(str(cap_patch_path), crop_bgr)

            if crop_bgr.shape[0] >= 5 and crop_bgr.shape[1] >= 5:
                head_pheno_item = classify_single_capitulum(
                    cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGB),
                    min_pappus_luminance=pheno_params["min_pappus_lum"],
                    max_pappus_saturation=pheno_params["max_pappus_sat"],
                    min_pappus_area_ratio=pheno_params["min_pappus_area"],
                    min_floret_yellow_saturation=pheno_params["min_floret_sat"],
                    min_floret_yellow_ratio=pheno_params["min_floret_rat"],
                    class_id=int(class_id),
                )
                head_pheno_results.append(head_pheno_item)

        metrics = compute_capitulum_metrics(
            local_cap_mask if local_cap_mask is not None else mask, bbox,
            pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
            crop_bgr=crop_bgr if crop_bgr.size > 0 else None,
            class_id=int(class_id), is_calibrated_mm=is_calibrated_mm,
        )
        if metrics is not None:
            metrics["score"] = score
            metrics["class_id"] = int(class_id)
            valid_capitula.append(metrics)

        mc_item = extract_single_capitulum_microchars(
            crop_bgr=crop_bgr if crop_bgr.size > 0 else None, sheet_bgr=img_bgr, bbox=(y0, x0, y1, x1),
            mask=local_cap_mask if local_cap_mask is not None else mask, class_id=int(class_id),
            pixels_per_mm=pixels_per_mm if is_calibrated_mm else None, is_calibrated_mm=is_calibrated_mm,
            paper_lab_ref=paper_lab_ref,
        )
        raw_cap_items.append({
            "bbox": bbox, "center_xy": (0.5 * float(xmin + xmax), 0.5 * float(ymin + ymax)),
            "mask": local_cap_mask if local_cap_mask is not None else mask,
            "crop_bgr": crop_bgr if crop_bgr.size > 0 else None,
            "score": float(score), "class_id": int(class_id), "metrics": metrics,
            "head_pheno": head_pheno_item, "microchars": mc_item,
        })

    repro_summary, voucher_extra_heads = summarize_repro_subset(
        catalog_number, raw_cap_items, pixels_per_mm, is_calibrated_mm, sheet_h, sheet_w, pid_val=0
    )
    pheno_rec = dict(classify_voucher_phenology(
        catalog_number, head_pheno_results, min_anthesis_heads_for_flowering=pheno_params["min_anth_heads"],
    ))
    pheno_rec["plant_individual_id"] = 0
    pheno_rec["total_capitula"] = int(pheno_rec.get("total_capitula", 0)) + int(voucher_extra_heads)
    pheno_rec["corymb_head_count_est"] = int(voucher_extra_heads)
    pheno_rec.setdefault("radiate_head_fraction", np.nan)
    repro_summary["radiate_head_fraction"] = pheno_rec.get("radiate_head_fraction", np.nan)

    return {
        "pheno_params": pheno_params, "valid_capitula": valid_capitula, "raw_cap_items": raw_cap_items,
        "repro_summary": repro_summary, "voucher_extra_heads": voucher_extra_heads,
        "pheno_rec": pheno_rec, "capitula_count": int(repro_summary["capitula_count"]),
        "paper_lab_ref": paper_lab_ref, "microchar_keys": MICROCHAR_REPRO_KEYS,
    }


def _inst_sheet_geom(inst_obj: DetectedInstance) -> Tuple[float, float, float, float, float]:
    """Compute sheet-space `(cx, cy, area, aspect_ratio, solidity)` for a `DetectedInstance`."""
    ox_i, oy_i = inst_obj.mask_origin
    m_u8 = inst_obj.mask if inst_obj.mask.dtype == np.uint8 else ((inst_obj.mask > 0).astype(np.uint8) * 255)
    bx_i, by_i, bw_i, bh_i = cv2.boundingRect(m_u8)
    if bw_i > 0 and bh_i > 0:
        sub_i = m_u8[by_i : by_i + bh_i, bx_i : bx_i + bw_i]
        mom_i = cv2.moments(sub_i, binaryImage=True)
        area_i = float(mom_i["m00"])
    else:
        sub_i, mom_i, area_i = None, None, 0.0
    if sub_i is not None and mom_i is not None and area_i > 0.0:
        cx_i = float(mom_i["m10"] / area_i) + float(bx_i + ox_i)
        cy_i = float(mom_i["m01"] / area_i) + float(by_i + oy_i)
        cnts_i, _ = cv2.findContours(sub_i, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if cnts_i:
            main_c = max(cnts_i, key=cv2.contourArea)
            rw, rh = cv2.minAreaRect(main_c)[1]
            ar_i = float(max(rw, rh)) / max(float(min(rw, rh)), 1.0)
            hull_a = float(cv2.contourArea(cv2.convexHull(main_c)))
            sol_i = float(area_i / hull_a) if hull_a > 0 else 1.0
        else:
            ar_i, sol_i = float(max(bw_i, bh_i)) / max(float(min(bw_i, bh_i)), 1.0), 1.0
    else:
        cx_i = 0.5 * float(inst_obj.bbox[1] + inst_obj.bbox[3])
        cy_i = 0.5 * float(inst_obj.bbox[0] + inst_obj.bbox[2])
        bh_b = max(1.0, float(inst_obj.bbox[2] - inst_obj.bbox[0]))
        bw_b = max(1.0, float(inst_obj.bbox[3] - inst_obj.bbox[1]))
        area_i, ar_i, sol_i = bh_b * bw_b, max(bh_b, bw_b) / max(min(bh_b, bw_b), 1.0), 1.0
    return cx_i, cy_i, area_i, ar_i, sol_i


def build_clustered_instances_and_basal_stats(
    pipe: Any,
    catalog_number: str,
    leaf_detections: Sequence[Any],
    root_detections: Sequence[Tuple[int, int, int, int]],
    sheet_h: int,
    sheet_w: int,
    pixels_per_mm: Optional[float],
    is_calibrated_mm: bool,
    spec_mask_sheet: Optional[np.ndarray],
) -> Tuple[
    List[DetectedInstance],
    Dict[int, Tuple[float, float, float, float, float]],
    Dict[int, Optional[Tuple[float, float]]],
    Dict[int, Tuple[float, float, float, float]],
]:
    """Construct `DetectedInstance` objects, cluster plant individuals, and compute basal rosette stats."""
    instances: List[DetectedInstance] = []
    for idx, (bbox, mask, score, class_id, crop_offset, mask_origin, _trimmed, trim_meta, _, _) in enumerate(leaf_detections, 1):
        is_tm = isinstance(trim_meta, dict)
        pw_px_raw = trim_meta.get("petiole_width_px") if is_tm else None
        pw_px = float(pw_px_raw) if (pw_px_raw is not None and np.isfinite(float(pw_px_raw)) and float(pw_px_raw) > 0.0) else None
        pw_mm = float(pw_px) / float(pixels_per_mm) if (pw_px is not None and is_calibrated_mm and pixels_per_mm and float(pixels_per_mm) > 0.0) else None
        instances.append(DetectedInstance(
            catalog_number=catalog_number, leaf_id=idx, bbox=bbox, mask=mask, score=score, class_id=class_id,
            pixels_per_mm=pixels_per_mm, is_calibrated_mm=is_calibrated_mm, crop_offset=crop_offset,
            mask_origin=mask_origin, trim_meta=trim_meta,
            petiole_mask=(trim_meta.get("petiole_mask") if is_tm else None),
            hole_mask=(trim_meta.get("hole_mask") if is_tm else None),
            keypoints_31=(trim_meta.get("keypoints_31") if is_tm else None),
            petiole_width_px=pw_px, petiole_width_mm=pw_mm,
        ))

    # 4. Stem-Aware Spatial Clustering (bridged Specimen_Segmenter mask + vertical shoot continuity)
    instances = cluster_plant_individuals(
        instances, sheet_w, sheet_h, specimen_mask=spec_mask_sheet,
        root_boxes=list(root_detections), pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
    )

    # 5. Basal Rosette Zone Statistics per plant_individual_id
    inst_geom_map = {inst.leaf_id: _inst_sheet_geom(inst) for inst in instances}
    lower_frac_cfg = float(getattr(getattr(pipe.config, "segmentation", None), "rosette_plant_bbox_lower_fraction", 0.60) or 0.60)
    basal_y_quantile = float(np.clip(1.0 - lower_frac_cfg, 0.05, 0.90))
    all_cx = np.array([inst_geom_map[i.leaf_id][0] for i in instances], dtype=np.float64)
    all_cy = np.array([inst_geom_map[i.leaf_id][1] for i in instances], dtype=np.float64)
    all_area = np.array([inst_geom_map[i.leaf_id][2] for i in instances], dtype=np.float64)
    sheet_b_mask = all_cy >= float(np.quantile(all_cy, basal_y_quantile)) if len(all_cy) > 0 else np.array([], dtype=bool)
    if len(sheet_b_mask) > 0 and not np.any(sheet_b_mask):
        sheet_b_mask = np.ones(len(all_cy), dtype=bool)
    sheet_basal_cx = float(np.median(all_cx[sheet_b_mask])) if len(all_cx) > 0 else (0.5 * sheet_w)
    sheet_basal_cy = float(np.median(all_cy[sheet_b_mask])) if len(all_cy) > 0 else (0.75 * sheet_h)
    sheet_med_area = float(np.median(all_area[sheet_b_mask])) if len(all_area) > 0 else 1.0
    sheet_radii = np.hypot(all_cx[sheet_b_mask] - sheet_basal_cx, all_cy[sheet_b_mask] - sheet_basal_cy) if len(all_cx) > 0 else np.array([0.0])
    sheet_rosette_rad = max(float(np.median(sheet_radii)), 0.04 * float(min(sheet_h, sheet_w))) if len(instances) >= 2 else 0.0

    clusters_by_pid: Dict[int, List[DetectedInstance]] = defaultdict(list)
    for inst in instances:
        clusters_by_pid[inst.plant_individual_id].append(inst)
    rosette_centroids_by_pid: Dict[int, Optional[Tuple[float, float]]] = {}
    basal_stats_by_pid: Dict[int, Tuple[float, float, float, float]] = {}

    for pid, pid_instances in sorted(clusters_by_pid.items()):
        p_cx = np.array([inst_geom_map[i.leaf_id][0] for i in pid_instances], dtype=np.float64)
        p_cy = np.array([inst_geom_map[i.leaf_id][1] for i in pid_instances], dtype=np.float64)
        p_area = np.array([inst_geom_map[i.leaf_id][2] for i in pid_instances], dtype=np.float64)
        b_mask = p_cy >= float(np.quantile(p_cy, basal_y_quantile))
        if not np.any(b_mask):
            b_mask = np.ones(len(pid_instances), dtype=bool)
        basal_cx_med, basal_cy_med, med_area = float(np.median(p_cx[b_mask])), float(np.median(p_cy[b_mask])), float(np.median(p_area[b_mask]))
        radii = np.hypot(p_cx[b_mask] - basal_cx_med, p_cy[b_mask] - basal_cy_med)
        med_rosette_radius = max(float(np.median(radii)), 0.04 * float(min(sheet_h, sheet_w))) if len(pid_instances) >= 2 else 0.0
        if (
            len(pid_instances) == 1 and len(instances) >= 2
            and (p_cy[0] < (sheet_basal_cy - 0.15 * float(sheet_h)) or p_area[0] < (0.55 * sheet_med_area))
        ):
            basal_stats_by_pid[pid] = (sheet_basal_cx, sheet_basal_cy, sheet_med_area, sheet_rosette_rad)
        else:
            basal_stats_by_pid[pid] = (basal_cx_med, basal_cy_med, med_area, med_rosette_radius)
        rosette_centroids_by_pid[pid] = (basal_cx_med, basal_cy_med) if len(pid_instances) >= 2 else None

    return instances, inst_geom_map, rosette_centroids_by_pid, basal_stats_by_pid
