#!/usr/bin/env python3
"""Step 02 Stage 2 crop mask reprojection, keypoint shifting, and YOLO26x-seg/pose inference."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np

from scripts.vision.capitulum_phenology_classifier import classify_voucher_phenology
from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES, LM3_REPRODUCTIVE_CLASS_IDS
from scripts.vision.lm3_inference_engines import _extract_native_lm3_petiole_anchor
from scripts.vision.step02_pipe_stage1_candidates import _unpack_reduced_lm3_crop_preds


def _reproject_crop_mask_to_local_win(
    cmask: Optional[np.ndarray],
    crop_box_yx: Tuple[int, int, int, int],
    win_box_yx: Tuple[int, int, int, int],
) -> Optional[np.ndarray]:
    """Reproject a crop-space binary mask into a padded local sheet window."""
    if cmask is None:
        return None
    cy0_c, cx0_c, _cy1_c, _cx1_c = crop_box_yx
    wy0_w, wx0_w, wy1_w, wx1_w = win_box_yx
    wh, ww = wy1_w - wy0_w, wx1_w - wx0_w
    if wh <= 0 or ww <= 0:
        return None
    out_m = np.zeros((wh, ww), dtype=np.uint8)
    ch_m, cw_m = cmask.shape[:2]
    iy0, ix0 = max(cy0_c, wy0_w), max(cx0_c, wx0_w)
    iy1, ix1 = min(cy0_c + ch_m, wy1_w), min(cx0_c + cw_m, wx1_w)
    if iy1 > iy0 and ix1 > ix0:
        out_m[iy0 - wy0_w : iy1 - wy0_w, ix0 - wx0_w : ix1 - wx0_w] = cmask[
            iy0 - cy0_c : iy1 - cy0_c, ix0 - cx0_c : ix1 - cx0_c
        ]
    return out_m if np.any(out_m > 0) else None


def _shift_crop_keypoints_to_local_win(
    kpts: Optional[Any],
    crop_box_yx: Tuple[int, int, int, int],
    win_box_yx: Tuple[int, int, int, int],
) -> Optional[Dict[str, Tuple[float, float, float]]]:
    """Shift crop-space 31-keypoint coordinates into local sheet window coordinates."""
    if not kpts or not isinstance(kpts, dict):
        return None
    cy0_c, cx0_c, _, _ = crop_box_yx
    wy0_w, wx0_w, _, _ = win_box_yx
    dx_w, dy_w = float(cx0_c - wx0_w), float(cy0_c - wy0_w)
    shifted: Dict[str, Tuple[float, float, float]] = {}
    for k_name, vals in kpts.items():
        if isinstance(vals, (tuple, list, np.ndarray)) and len(vals) >= 2:
            conf_v = float(vals[2]) if len(vals) >= 3 else 1.0
            shifted[str(k_name)] = (float(vals[0]) + dx_w, float(vals[1]) + dy_w, conf_v)
    return shifted or None


def _run_batched_crop_pose_for_voucher(
    pipe: Any,
    candidate_s1_crops: Sequence[Any],
    all_crop_preds: Sequence[Any],
    spec_mask_sheet: Optional[np.ndarray],
    sheet_h: int,
    sheet_w: int,
) -> List[Optional[Any]]:
    """Run batched YOLO26x-pose landmark inference across Stage 1 crops containing lamina masks."""
    if not (
        candidate_s1_crops and all_crop_preds
        and getattr(pipe, "pose_engine", None) is not None
        and getattr(pipe.pose_engine, "model", None) is not None
        and "predict_landmarks" not in getattr(pipe.pose_engine, "__dict__", {})
        and hasattr(pipe.pose_engine, "predict_batch")
    ):
        return []

    all_crop_pose_preds: List[Optional[Any]] = [None] * len(candidate_s1_crops)
    p_idxs: List[int] = []
    p_crops: List[np.ndarray] = []
    p_lamina_masks: List[Optional[np.ndarray]] = []
    p_petiole_masks: List[Optional[np.ndarray]] = []
    for idx_c, (c_item, c_preds) in enumerate(zip(candidate_s1_crops, all_crop_preds)):
        if not c_preds or not any(int(p_cls) == 0 and np.any(pm_u8 > 0) for _pb, pm_u8, _p_sc, p_cls in c_preds):
            continue
        cy0_p, cx0_p, cy1_p, cx1_p = c_item[4]
        spec_crop_p = (
            spec_mask_sheet[cy0_p:cy1_p, cx0_p:cx1_p]
            if (spec_mask_sheet is not None and spec_mask_sheet.shape[:2] == (sheet_h, sheet_w))
            else None
        )
        bm_c, _, pm_c, _ = _unpack_reduced_lm3_crop_preds(
            c_preds, c_item[3].shape[:2],
            specimen_paper_cleaned_mask=spec_crop_p,
            subtract_holes=bool(getattr(pipe, "subtract_holes", True)),
        )
        p_idxs.append(idx_c)
        p_crops.append(c_item[3])
        p_lamina_masks.append(bm_c)
        p_petiole_masks.append(pm_c)

    if p_crops:
        pose_out = pipe.pose_engine.predict_batch(
            p_crops, lamina_masks=p_lamina_masks, petiole_masks=p_petiole_masks,
            batch_size=max(8, pipe.stage2_batch_size),
        )
        for idx_c, p_res in zip(p_idxs, pose_out):
            all_crop_pose_preds[idx_c] = p_res
    return all_crop_pose_preds


def _run_single_stage_fallback(
    pipe: Any,
    img_bgr: np.ndarray,
    sheet_h: int,
    sheet_w: int,
    precomputed_single_stage_detections: Optional[Union[List[Any], Tuple[List[Any], bool]]],
) -> Tuple[List[Any], List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]], List[Tuple[int, int, int, int]]]:
    """Direct single-stage fallback when Stage 1 detector is uninitialized or directly mocked."""
    if precomputed_single_stage_detections is not None:
        raw_detections = (
            precomputed_single_stage_detections[0]
            if isinstance(precomputed_single_stage_detections, tuple)
            else precomputed_single_stage_detections
        )
    else:
        raw_detections = pipe.model_engine.predict(img_bgr)

    capitulum_detections: List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]] = []
    root_detections: List[Tuple[int, int, int, int]] = []
    raw_s2_leaves: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]] = []
    raw_s2_petioles: List[Tuple[Tuple[int, int, int, int], np.ndarray, float]] = []
    for bbox, mask, score, cls_id in raw_detections:
        if cls_id in LM3_REPRODUCTIVE_CLASS_IDS:
            capitulum_detections.append((bbox, mask, float(score), int(cls_id)))
        elif int(cls_id) in (7, 9):
            root_detections.append(bbox)
        elif int(cls_id) == 0:
            raw_s2_leaves.append((bbox, mask, float(score), int(cls_id)))
        elif int(cls_id) == 1:
            raw_s2_petioles.append((bbox, mask, float(score)))

    leaf_detections: List[Any] = []
    for bbox, mask, score, cls_id in raw_s2_leaves:
        sheet_mask, full_bbox = mask, bbox
        fg_area = int(np.count_nonzero(sheet_mask))
        if fg_area <= 0:
            continue
        native_pet_meta: Optional[Dict[str, Any]] = None
        if raw_s2_petioles and sheet_mask.shape[:2] == (sheet_h, sheet_w):
            ly0, lx0, ly1, lx1 = [int(v) for v in full_bbox]
            prox = max(32, int(0.35 * max(ly1 - ly0, lx1 - lx0, 1)))
            nearby_pets = [
                (pscore, pbbox, pmask)
                for pbbox, pmask, pscore in raw_s2_petioles
                if pmask is not None and pmask.shape[:2] == sheet_mask.shape[:2]
                and not (
                    int(pbbox[2]) < ly0 - prox or int(pbbox[0]) > ly1 + prox
                    or int(pbbox[3]) < lx0 - prox or int(pbbox[1]) > lx1 + prox
                )
            ]
            if nearby_pets:
                nearby_pets.sort(key=lambda item: item[0], reverse=True)
                _, pbbox_best, pmask_best = nearby_pets[0]
                jy0, jx0 = max(0, min(ly0, int(pbbox_best[0])) - 8), max(0, min(lx0, int(pbbox_best[1])) - 8)
                jy1, jx1 = min(sheet_h, max(ly1, int(pbbox_best[2])) + 8), min(sheet_w, max(lx1, int(pbbox_best[3])) + 8)
                if jy1 > jy0 and jx1 > jx0:
                    native_pet_meta = _extract_native_lm3_petiole_anchor(
                        sheet_mask[jy0:jy1, jx0:jx1], pmask_best[jy0:jy1, jx0:jx1],
                    )
                    if native_pet_meta is not None:
                        p_bc, p_ac = native_pet_meta["p_base_crop"], native_pet_meta["p_apex_crop"]
                        native_pet_meta["p_base"] = (int(jx0 + p_bc[0]), int(jy0 + p_bc[1]))
                        native_pet_meta["p_apex"] = (int(jx0 + p_ac[0]), int(jy0 + p_ac[1]))

        leaf_detections.append((
            full_bbox, sheet_mask, float(score), int(cls_id),
            (int(full_bbox[1]), int(full_bbox[0])), (0, 0), False, native_pet_meta, 0.0, fg_area,
        ))
    return leaf_detections, capitulum_detections, root_detections


def run_stage2_leaf_and_pose_inference(
    pipe: Any,
    catalog_number: str,
    image_path: Path,
    img_bgr: np.ndarray,
    sheet_h: int,
    sheet_w: int,
    is_high_dpi: bool,
    pixels_per_mm: Optional[float],
    cf_source: str,
    stage1_detections: Sequence[Any],
    spec_mask_sheet: Optional[np.ndarray],
    precomputed_s2_payload: Optional[Tuple[Any, ...]] = None,
    precomputed_single_stage_detections: Optional[Union[List[Any], Tuple[List[Any], bool]]] = None,
) -> Dict[str, Any]:
    """Run Stage 2 YOLO26x-seg + YOLO26x-pose crop inference, reprojection, and zero-detection check."""
    capitulum_detections: List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]] = []
    root_detections: List[Tuple[int, int, int, int]] = []
    upper_cauline_pts: List[Tuple[float, float]] = []
    rosette_bbox_fallback_used: bool = False
    leaf_detections: List[Any] = []

    if stage1_detections:
        all_crop_pose_preds: List[Optional[Any]] = []
        if precomputed_s2_payload is not None:
            if len(precomputed_s2_payload) >= 11:
                (
                    candidate_s1_crops, capitulum_detections, root_detections, upper_cauline_pts,
                    rosette_bbox_fallback_used, _rej_lflt, _s1_cal, _s1_px_cm, _s1_min_cm2,
                    all_crop_preds, all_crop_pose_preds,
                ) = precomputed_s2_payload[:11]
            else:
                (
                    candidate_s1_crops, capitulum_detections, root_detections, upper_cauline_pts,
                    rosette_bbox_fallback_used, _rej_lflt, _s1_cal, _s1_px_cm, _s1_min_cm2, all_crop_preds,
                ) = precomputed_s2_payload[:10]
        else:
            (
                candidate_s1_crops, capitulum_detections, root_detections, upper_cauline_pts,
                rosette_bbox_fallback_used, _rej_lflt, _s1_cal, _s1_px_cm, _s1_min_cm2,
            ) = pipe._unpack_stage1_candidates(
                catalog_number=catalog_number, img_bgr=img_bgr,
                stage1_detections=stage1_detections, pixels_per_mm=pixels_per_mm, cf_source=cf_source,
            )
            all_crop_preds = pipe.model_engine.predict_batch([c[3] for c in candidate_s1_crops]) if candidate_s1_crops else []
            all_crop_pose_preds = _run_batched_crop_pose_for_voucher(
                pipe, candidate_s1_crops, all_crop_preds, spec_mask_sheet, sheet_h, sheet_w
            )

        for crop_idx_s2, ((_s1_bbox, _s1_score, s1_cls, crop_bgr, crop_box), crop_preds) in enumerate(
            zip(candidate_s1_crops, all_crop_preds)
        ):
            if not crop_preds:
                continue
            cy0, cx0, cy1, cx1 = crop_box
            spec_crop_s2 = (
                spec_mask_sheet[cy0:cy1, cx0:cx1]
                if (spec_mask_sheet is not None and spec_mask_sheet.shape[:2] == (sheet_h, sheet_w))
                else None
            )
            best_mask, best_score, best_petiole_mask, best_hole_mask = _unpack_reduced_lm3_crop_preds(
                crop_preds, crop_bgr.shape[:2],
                specimen_paper_cleaned_mask=spec_crop_s2,
                subtract_holes=bool(getattr(pipe, "subtract_holes", True)),
            )
            if best_mask is None:
                continue

            keypoints_31_crop: Optional[Dict[str, Tuple[float, float, float]]] = None
            petiole_width_px_crop: Optional[float] = None
            is_geom_fb_pose_crop = True
            pre_pose_crop = all_crop_pose_preds[crop_idx_s2] if crop_idx_s2 < len(all_crop_pose_preds) else None
            if pre_pose_crop is not None and bool(getattr(pipe, "do_apply_landmarks", True)):
                keypoints_31_crop = getattr(pre_pose_crop, "keypoints_31", None) or None
                petiole_width_px_crop = getattr(pre_pose_crop, "petiole_width_px", None)
                is_geom_fb_pose_crop = bool(getattr(pre_pose_crop, "is_geometric_fallback", False))
            elif (
                bool(getattr(pipe, "do_apply_landmarks", True))
                and getattr(pipe, "pose_engine", None) is not None
                and (getattr(pipe.pose_engine, "model", None) is not None or "predict_landmarks" in getattr(pipe.pose_engine, "__dict__", {}))
            ):
                try:
                    pose_res = pipe.pose_engine.predict_landmarks(
                        crop_bgr, lamina_mask=best_mask, petiole_mask=best_petiole_mask
                    )
                    keypoints_31_crop = pose_res.keypoints_31 or None
                    petiole_width_px_crop = pose_res.petiole_width_px
                    is_geom_fb_pose_crop = bool(getattr(pose_res, "is_geometric_fallback", False))
                except Exception:
                    keypoints_31_crop = None

            m_raw = best_mask
            fg_area = int(np.count_nonzero(m_raw))
            if fg_area <= 0:
                continue
            bx_c, by_c, bw_c, bh_c = cv2.boundingRect(m_raw)
            if bw_c == 0 or bh_c == 0:
                continue
            sub_raw = m_raw[by_c : by_c + bh_c, bx_c : bx_c + bw_c]
            cnts_m, _ = cv2.findContours(sub_raw, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if not cnts_m:
                continue
            hull_m_area = float(cv2.contourArea(cv2.convexHull(max(cnts_m, key=cv2.contourArea))))
            m_solidity = float(fg_area) / hull_m_area if hull_m_area > 0 else 1.0

            ymin_s, ymax_s = max(0, cy0 + by_c), min(sheet_h, cy0 + by_c + bh_c)
            xmin_s, xmax_s = max(0, cx0 + bx_c), min(sheet_w, cx0 + bx_c + bw_c)
            if ymax_s <= ymin_s or xmax_s <= xmin_s:
                continue
            full_bbox = (ymin_s, xmin_s, ymax_s, xmax_s)

            margin = max(24, int(0.35 * max(bw_c, bh_c)))
            win_y0, win_x0 = max(0, ymin_s - margin), max(0, xmin_s - margin)
            win_y1, win_x1 = min(sheet_h, ymax_s + margin), min(sheet_w, xmax_s + margin)
            local_mask = np.zeros((win_y1 - win_y0, win_x1 - win_x0), dtype=np.uint8)
            dst_y0, dst_x0 = ymin_s - win_y0, xmin_s - win_x0
            eff_h, eff_w = min(bh_c, win_y1 - ymin_s), min(bw_c, win_x1 - xmin_s)
            local_mask[dst_y0 : dst_y0 + eff_h, dst_x0 : dst_x0 + eff_w] = sub_raw[:eff_h, :eff_w]

            native_pet_meta = _extract_native_lm3_petiole_anchor(m_raw, best_petiole_mask, keypoints_31=keypoints_31_crop)
            if native_pet_meta is not None:
                p_bc, p_ac = native_pet_meta["p_base_crop"], native_pet_meta["p_apex_crop"]
                native_pet_meta["p_base"] = (int(cx0 + p_bc[0] - win_x0), int(cy0 + p_bc[1] - win_y0))
                native_pet_meta["p_apex"] = (int(cx0 + p_ac[0] - win_x0), int(cy0 + p_ac[1] - win_y0))
            else:
                native_pet_meta = {}
            win_box = (win_y0, win_x0, win_y1, win_x1)
            native_pet_meta["petiole_mask"] = _reproject_crop_mask_to_local_win(best_petiole_mask, crop_box, win_box)
            native_pet_meta["hole_mask"] = _reproject_crop_mask_to_local_win(best_hole_mask, crop_box, win_box)
            native_pet_meta["keypoints_31"] = _shift_crop_keypoints_to_local_win(keypoints_31_crop, crop_box, win_box)
            native_pet_meta["is_geometric_fallback_pose"] = bool(is_geom_fb_pose_crop)
            native_pet_meta["petiole_width_px"] = petiole_width_px_crop

            ch_m, cw_m = m_raw.shape[:2]
            clip_mask = m_raw > 0
            if best_petiole_mask is not None and best_petiole_mask.shape[:2] == (ch_m, cw_m):
                pet_dil_edge = cv2.dilate((best_petiole_mask > 0).astype(np.uint8), np.ones((9, 9), dtype=np.uint8), iterations=1) > 0
                clip_mask = clip_mask & (~pet_dil_edge)
            max_edge_allow = max(8, int(round(0.08 * float(max(bw_c, bh_c)))))
            sheet_edge_allow = max(6, int(round(0.05 * float(max(bw_c, bh_c)))))
            native_pet_meta["is_crop_clipped"] = bool(
                ch_m >= 16 and cw_m >= 16 and (
                    int(np.count_nonzero(clip_mask[0, :])) > (sheet_edge_allow if cy0 <= 2 else max_edge_allow)
                    or int(np.count_nonzero(clip_mask[-1, :])) > (sheet_edge_allow if cy1 >= sheet_h - 2 else max_edge_allow)
                    or int(np.count_nonzero(clip_mask[:, 0])) > (sheet_edge_allow if cx0 <= 2 else max_edge_allow)
                    or int(np.count_nonzero(clip_mask[:, -1])) > (sheet_edge_allow if cx1 >= sheet_w - 2 else max_edge_allow)
                )
            )
            leaf_detections.append((
                full_bbox, local_mask, float(best_score), int(s1_cls),
                (crop_box[1], crop_box[0]), (win_x0, win_y0), False, native_pet_meta, m_solidity, fg_area,
            ))
    elif not is_high_dpi or pipe.stage1_detector.model is None or "predict" in pipe.model_engine.__dict__:
        leaf_detections, capitulum_detections, root_detections = _run_single_stage_fallback(
            pipe, img_bgr, sheet_h, sheet_w, precomputed_single_stage_detections
        )

    cal_src_set = {str(s).strip().lower() for s in _CALIBRATED_SCALE_SOURCES}
    is_calibrated_mm = bool(
        pixels_per_mm is not None and not np.isnan(pixels_per_mm)
        and float(pixels_per_mm) > 0.0 and str(cf_source).strip().lower() in cal_src_set
    )

    early_exit = None
    if not leaf_detections and not capitulum_detections:
        failure_record = {
            "catalogNumber": catalog_number, "image_path": str(image_path),
            "failure_reason": "no_leaves_detected", "details": "Model produced 0 candidate leaf segmentations",
            "pixels_per_mm": round(float(pixels_per_mm), 4) if is_calibrated_mm else np.nan,
            "is_calibrated_mm": bool(is_calibrated_mm),
        }
        empty_repro = {
            "catalogNumber": catalog_number, "plant_individual_id": 0, "capitula_count": 0,
            "corymb_head_count_est": 0, "involucre_height_px": np.nan, "involucre_width_px": np.nan,
            "involucre_height_mm": np.nan, "involucre_width_mm": np.nan, "AR_inv": np.nan,
            "plant_height_cm": np.nan, "scape_to_basal_leaf_ratio": np.nan,
            "is_bent_stem": False, "radiate_head_fraction": np.nan,
        }
        pheno_cfg = getattr(getattr(pipe.config, "thresholds", None), "capitulum_phenology", None)
        min_anth_heads = int(getattr(pheno_cfg, "min_anthesis_heads_for_flowering", 1) or 1)
        empty_pheno = dict(classify_voucher_phenology(catalog_number, [], min_anthesis_heads_for_flowering=min_anth_heads))
        empty_pheno["plant_individual_id"] = 0
        empty_pheno["corymb_head_count_est"] = 0
        empty_pheno.setdefault("radiate_head_fraction", np.nan)
        early_exit = (failure_record, empty_repro, empty_pheno)

    return {
        "leaf_detections": leaf_detections, "capitulum_detections": capitulum_detections,
        "root_detections": root_detections, "upper_cauline_pts": upper_cauline_pts,
        "rosette_bbox_fallback_used": rosette_bbox_fallback_used,
        "is_calibrated_mm": is_calibrated_mm, "early_exit": early_exit,
    }
