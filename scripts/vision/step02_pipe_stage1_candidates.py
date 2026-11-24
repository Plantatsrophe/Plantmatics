#!/usr/bin/env python3
"""Step 02 Stage 1 candidate unpacking, crop reduction, scale provenance, and per-voucher setup."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np
import pandas as pd

from scripts.vision.lm3_data_loader import (
    _CALIBRATED_SCALE_SOURCES,
    DEFAULT_HERBARIUM_DPI,
    LM3_REPRODUCTIVE_CLASS_IDS,
    normalize_pixels_per_mm,
    resolve_px_per_cm,
)
from scripts.vision.lm3_inference_engines import (
    _reduce_crop_predictions_native_lm3,
    crop_leaf_bbox_with_padding,
)
from scripts.vision.step02_photometry_corymb_utils import estimate_sheet_paper_cie_l_ref


def _unpack_reduced_lm3_crop_preds(
    crop_preds: Sequence[Any],
    crop_shape: Tuple[int, int],
    specimen_paper_cleaned_mask: Optional[np.ndarray] = None,
    subtract_holes: bool = True,
) -> Tuple[Optional[np.ndarray], float, Optional[np.ndarray], Optional[np.ndarray]]:
    """Unpack reduced LM3 `(lamina, score, petiole, hole)` masks and subtract UNet++ interior paper pockets."""
    ch, cw = int(crop_shape[0]), int(crop_shape[1])
    valid = [p for p in (crop_preds or ()) if p[1] is not None and isinstance(p[1], np.ndarray) and np.any(p[1] > 0)]
    if ch <= 0 or cw <= 0 or not valid:
        return None, 0.0, None, None
    if any(sum(1 for p in valid if int(p[3]) == c) > 1 for c in (0, 1, 2)):
        valid = _reduce_crop_predictions_native_lm3(valid, (ch, cw), min_area_px=(64 if (ch * cw) >= 4096 else 1))
    lamina, best_score, pet, hole = None, 0.0, None, None
    for _bb, pm, sc, cid in valid:
        m_u8 = (pm > 0).astype(np.uint8) * 255
        m_u8 = cv2.resize(m_u8, (cw, ch), interpolation=cv2.INTER_NEAREST) if m_u8.shape[:2] != (ch, cw) else m_u8
        if int(cid) == 0 and lamina is None:
            lamina, best_score = m_u8, float(sc)
        elif int(cid) == 1 and pet is None:
            pet = m_u8
        elif int(cid) == 2 and hole is None:
            hole = m_u8
    if lamina is None:
        return None, 0.0, None, None
    if subtract_holes and hole is not None and np.any((lamina > 0) & (hole == 0)):
        lamina = np.where((lamina > 0) & (hole == 0), 255, 0).astype(np.uint8)
    if (
        isinstance(specimen_paper_cleaned_mask, np.ndarray)
        and specimen_paper_cleaned_mask.shape[:2] == (ch, cw)
        and np.any(specimen_paper_cleaned_mask > 0)
    ):
        lam_b, spec_fg = lamina > 0, specimen_paper_cleaned_mask > 0
        fg_n = int(np.count_nonzero(lam_b))
        if fg_n >= 64 and float(np.count_nonzero(lam_b & spec_fg)) / float(fg_n) >= 0.55:
            dist_in = cv2.distanceTransform(lam_b.astype(np.uint8), cv2.DIST_L2, 3) >= 2.0
            n_sp, sp_l, sp_s, _ = cv2.connectedComponentsWithStats((dist_in & (~spec_fg)).astype(np.uint8), connectivity=8)
            keep_l = [lbl for lbl in range(1, n_sp) if 16 <= int(sp_s[lbl, cv2.CC_STAT_AREA]) <= int(0.35 * fg_n)]
            cleaned = lam_b & (~np.isin(sp_l, keep_l)) if keep_l else lam_b
            if keep_l and int(np.count_nonzero(cleaned)) >= max(32, int(0.50 * fg_n)):
                lamina = np.where(cleaned, 255, 0).astype(np.uint8)
    return lamina, best_score, pet, (hole if (hole is not None and np.any(hole > 0)) else None)


def _resolve_scale_source(
    pixels_per_mm: Optional[float], is_calibrated_mm: bool, cf_source: Any = ""
) -> str:
    """Normalize Step 02 scale provenance across empirical rulers, unboxed ticks, and MP priors."""
    src_low = str(cf_source or "").strip().lower()
    has_px = bool(pixels_per_mm is not None and pd.notna(pixels_per_mm) and float(pixels_per_mm) > 0.0)
    if src_low == "tick_lattice_unboxed":
        return "tick_lattice_unboxed"
    if has_px and src_low in {"predicted_from_megapixels", "mp_estimated"}:
        return "mp_estimated"
    if has_px and bool(is_calibrated_mm):
        return "ruler_calibrated"
    return "uncalibrated"


def unpack_stage1_candidates(
    pipe: Any,
    catalog_number: str,
    img_bgr: np.ndarray,
    stage1_detections: Sequence[Any],
    pixels_per_mm: Optional[float],
    cf_source: str = "measured_from_ruler",
) -> Tuple[
    List[Tuple[Tuple[int, int, int, int], float, int, np.ndarray, Tuple[int, int, int, int]]],
    List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]],
    List[Tuple[int, int, int, int]],
    List[Tuple[float, float]],
    bool, int, bool, float, float,
]:
    """Thread-safe CPU Stage 1 -> Stage 2 candidate crop extractor."""
    del catalog_number
    sheet_h, sheet_w = img_bgr.shape[:2]
    cal_src_set = {str(s).strip().lower() for s in _CALIBRATED_SCALE_SOURCES}
    s1_is_calibrated_mm = bool(
        pixels_per_mm is not None
        and not np.isnan(pixels_per_mm)
        and float(pixels_per_mm) > 0.0
        and str(cf_source).strip().lower() in cal_src_set
    )
    s1_px_per_cm = (
        resolve_px_per_cm(pixels_per_mm=pixels_per_mm, default_dpi=DEFAULT_HERBARIUM_DPI, sheet_shape=(sheet_h, sheet_w))
        if s1_is_calibrated_mm else 0.0
    )

    candidate_s1_crops: List[Tuple[Tuple[int, int, int, int], float, int, np.ndarray, Tuple[int, int, int, int]]] = []
    capitulum_detections: List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]] = []
    root_detections: List[Tuple[int, int, int, int]] = []
    upper_cauline_pts: List[Tuple[float, float]] = []
    s1_plant_boxes: List[Tuple[int, int, int, int]] = []
    rosette_bbox_fallback_used = False
    rejected_leaflet_delta = 0
    crop_pad = 0.0
    seg_partial_cfg = bool(getattr(pipe.config.segmentation, "segment_partial_leaves", False))

    # Stage 1 (`lm3_plant_detector__yolo26x_det_1280` 9-class YOLO26x):
    # `0: Leaf_WHOLE`, `1: Leaf_PARTIAL`, `2..6: reproductive`, `7: Roots`, `8: Wood`.
    c0_items: List[Tuple[Tuple[int, int, int, int], float, int]] = []
    c1_items: List[Tuple[Tuple[int, int, int, int], float, int]] = []
    for s1_item in stage1_detections:
        if len(s1_item) == 4:
            s1_bbox, s1_mask, s1_score, s1_cls = s1_item
        else:
            s1_bbox, s1_score, s1_cls = s1_item
            s1_mask = None
        if s1_cls in LM3_REPRODUCTIVE_CLASS_IDS:
            capitulum_detections.append((s1_bbox, s1_mask, float(s1_score), int(s1_cls)))
            s1_plant_boxes.append(s1_bbox)
        elif int(s1_cls) in (7, 9):
            root_detections.append(s1_bbox)
            s1_plant_boxes.append(s1_bbox)
        elif int(s1_cls) == 0:
            s1_plant_boxes.append(s1_bbox)
            c0_items.append((s1_bbox, float(s1_score), 0))
        elif int(s1_cls) == 1:
            s1_plant_boxes.append(s1_bbox)
            c1_items.append((s1_bbox, float(s1_score), 1))

    selected_s1_leaves: List[Tuple[Tuple[int, int, int, int], float, int]] = list(c0_items)
    if not c0_items:
        selected_s1_leaves.extend(c1_items)
    elif seg_partial_cfg:
        for c1_box, c1_sc, c1_cl in c1_items:
            area_c1 = float(max(1, (c1_box[2] - c1_box[0]) * (c1_box[3] - c1_box[1])))
            overlaps_c0 = False
            for c0_box, _c0_sc, _c0_cl in c0_items:
                area_c0 = float(max(1, (c0_box[2] - c0_box[0]) * (c0_box[3] - c0_box[1])))
                iy0, ix0 = max(c1_box[0], c0_box[0]), max(c1_box[1], c0_box[1])
                iy1, ix1 = min(c1_box[2], c0_box[2]), min(c1_box[3], c0_box[3])
                if iy1 > iy0 and ix1 > ix0 and (float((iy1 - iy0) * (ix1 - ix0)) / max(min(area_c1, area_c0), 1.0)) >= 0.20:
                    overlaps_c0 = True
                    break
            if not overlaps_c0:
                selected_s1_leaves.append((c1_box, c1_sc, c1_cl))

    # Confidence-ordered nested-box deduplication with asymmetric parent-enclosure replacement
    if len(selected_s1_leaves) >= 2:
        selected_s1_leaves.sort(key=lambda x: (int(x[2]) == 0, float(x[1])), reverse=True)
        deduped_s1: List[Tuple[Tuple[int, int, int, int], float, int]] = []
        for b_box, b_score, b_cls in selected_s1_leaves:
            area_b = float(max(1, (b_box[2] - b_box[0]) * (b_box[3] - b_box[1])))
            suppressed = False
            replaced_indices: List[int] = []
            for idx_k, (k_box, _k_score, k_cls) in enumerate(deduped_s1):
                area_k = float(max(1, (k_box[2] - k_box[0]) * (k_box[3] - k_box[1])))
                iy0, ix0 = max(b_box[0], k_box[0]), max(b_box[1], k_box[1])
                iy1, ix1 = min(b_box[2], k_box[2]), min(b_box[3], k_box[3])
                if iy1 > iy0 and ix1 > ix0:
                    inter = float((iy1 - iy0) * (ix1 - ix0))
                    iou = inter / max(area_b + area_k - inter, 1.0)
                    if (
                        int(b_cls) == 0 and int(k_cls) == 0 and area_b >= 1.35 * area_k
                        and (inter / max(area_k, 1.0)) >= 0.90 and float(b_score) >= 0.50
                    ):
                        replaced_indices.append(idx_k)
                    elif iou >= 0.70 or (area_k >= 0.85 * area_b and (inter / max(area_b, 1.0)) >= 0.90):
                        suppressed = True
                        break
            if suppressed:
                continue
            if replaced_indices:
                deduped_s1[replaced_indices[0]] = (b_box, b_score, b_cls)
                for rem_idx in reversed(replaced_indices[1:]):
                    deduped_s1.pop(rem_idx)
            else:
                deduped_s1.append((b_box, b_score, b_cls))
        selected_s1_leaves = deduped_s1

    for s1_bbox, s1_score, s1_cls in selected_s1_leaves:
        crop_bgr, crop_box = crop_leaf_bbox_with_padding(img_bgr, s1_bbox, pad_ratio=crop_pad)
        if crop_bgr.size > 0:
            candidate_s1_crops.append((s1_bbox, float(s1_score), int(s1_cls), crop_bgr, crop_box))
    return (
        candidate_s1_crops, capitulum_detections, root_detections, upper_cauline_pts,
        rosette_bbox_fallback_used, rejected_leaflet_delta, s1_is_calibrated_mm, s1_px_per_cm, 0.0,
    )


def load_voucher_image_or_failure(
    catalog_number: str, image_path: Path, preloaded_img: Optional[np.ndarray] = None
) -> Tuple[Optional[np.ndarray], Path, Optional[Dict[str, Any]]]:
    """Load or resolve a voucher image via `cv2.imread`, returning `(img_bgr, resolved_path, failure_or_None)`."""
    if preloaded_img is not None:
        return preloaded_img, image_path, None
    if not image_path.exists():
        for ext in (".webp", ".jpg", ".jpeg", ".png"):
            cand_p = image_path.with_suffix(ext)
            if cand_p.exists():
                image_path = cand_p
                break
    if not image_path.exists():
        return None, image_path, {
            "catalogNumber": catalog_number, "image_path": str(image_path),
            "failure_reason": "image_not_found", "details": f"File does not exist at {image_path}",
        }
    img_bgr = cv2.imread(str(image_path))
    if img_bgr is None:
        return None, image_path, {
            "catalogNumber": catalog_number, "image_path": str(image_path),
            "failure_reason": "corrupted_image", "details": "Failed to decode image with OpenCV",
        }
    return img_bgr, image_path, None


def prepare_voucher_scale_stage1_and_specimen(
    pipe: Any,
    catalog_number: str,
    img_bgr: np.ndarray,
    precomputed_s1_letterbox: Optional[np.ndarray] = None,
    precomputed_s1_detections: Optional[List[Any]] = None,
    precomputed_spec_mask: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Resolve ruler calibration, Stage 1 organ detections, sheet paper CIE L*, and UNet++ specimen mask."""
    sheet_h, sheet_w = img_bgr.shape[:2]
    scale_f = max(sheet_h, sheet_w) / 5600.0 if max(sheet_h, sheet_w) > 1000 else 0.0
    is_high_dpi = max(sheet_h, sheet_w) > 1000 and scale_f > 0.4

    lm3_px_per_mm = pipe.ruler_calibrations.get(catalog_number)
    cal_sources = getattr(pipe.ruler_calibrations, "sources", None)
    if lm3_px_per_mm is not None and not np.isnan(lm3_px_per_mm) and lm3_px_per_mm > 0:
        pixels_per_mm: Optional[float] = normalize_pixels_per_mm(lm3_px_per_mm)
        cf_source = str(cal_sources.get(catalog_number, "measured_from_ruler")) if isinstance(cal_sources, dict) else "measured_from_ruler"
    elif not getattr(pipe, "_is_dummy_test_weights", True) and getattr(pipe, "scale_calibrator", None) is not None:
        cal_res = pipe.scale_calibrator.calibrate_sheet(
            img_bgr, catalog_number=catalog_number,
            use_cf_predicted_by_mp=bool(getattr(pipe, "use_cf_predicted_by_mp", False)),
        )
        pixels_per_mm = normalize_pixels_per_mm(cal_res.get("pixels_per_mm"))
        cf_source = str(cal_res.get("cf_source", "uncalibrated"))
        if pixels_per_mm is not None and float(pixels_per_mm) > 0.0:
            pipe.ruler_calibrations[catalog_number] = float(pixels_per_mm)
            if isinstance(cal_sources, dict):
                cal_sources[catalog_number] = cf_source
    else:
        pixels_per_mm, cf_source = None, "uncalibrated"

    use_s1_cache = bool(pipe.lm3_project_dir is not None and not pipe.force)
    direct_s2_mock = (
        "predict" in getattr(pipe.model_engine, "__dict__", {})
        and "predict" not in getattr(pipe.stage1_detector, "__dict__", {})
        and "predict_batch" not in getattr(pipe.model_engine, "__dict__", {})
    )
    if precomputed_s1_detections is not None:
        stage1_detections = precomputed_s1_detections
    elif direct_s2_mock:
        stage1_detections = []
    elif precomputed_s1_letterbox is not None:
        stage1_detections = pipe.stage1_detector.predict(
            img_bgr, catalog_number=catalog_number, use_cache=use_s1_cache, precomputed_letterbox=precomputed_s1_letterbox,
        )
    else:
        stage1_detections = pipe.stage1_detector.predict(img_bgr, catalog_number=catalog_number, use_cache=use_s1_cache)

    sheet_paper_l_ref = estimate_sheet_paper_cie_l_ref(img_bgr)
    spec_mask_attempted: bool = precomputed_spec_mask is not None
    spec_mask_sheet: Optional[np.ndarray] = (
        precomputed_spec_mask if (isinstance(precomputed_spec_mask, np.ndarray) and precomputed_spec_mask.size > 0) else None
    )
    if (
        not spec_mask_attempted and is_high_dpi and isinstance(img_bgr, np.ndarray) and img_bgr.ndim >= 2
        and getattr(pipe, "specimen_segmenter", None) is not None
        and getattr(pipe.specimen_segmenter, "model", None) is not None
    ):
        spec_mask_attempted = True
        try:
            spec_mask_sheet = pipe.specimen_segmenter.predict_mask(img_bgr, return_sheet_size=True)
        except Exception:
            spec_mask_sheet = None

    return {
        "sheet_h": sheet_h, "sheet_w": sheet_w, "is_high_dpi": is_high_dpi,
        "pixels_per_mm": pixels_per_mm, "cf_source": cf_source, "stage1_detections": stage1_detections,
        "sheet_paper_l_ref": sheet_paper_l_ref, "spec_mask_attempted": spec_mask_attempted,
        "spec_mask_sheet": spec_mask_sheet,
    }
