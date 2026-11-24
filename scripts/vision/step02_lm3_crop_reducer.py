#!/usr/bin/env python3
"""
Step 02 LeafMachine3 Stage 2 native crop-prediction reducer separating per-instance
lamina masks from petiole and hole masks without synthetic bridging pixels.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np


def _reduce_crop_predictions_native_lm3(
    raw_dets: Sequence[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]],
    crop_shape: Tuple[int, int],
    min_area_px: int = 64,
) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
    """
    Reduces raw per-instance YOLO26x-seg crop predictions to canonical LeafMachine3
    masks matching `YoloLeafSegmenter.build_instances` and `Reporter._build_mask_products`
    while preserving distinct `cls_id == 0` (`Leaf`) instances so touching neighbor leaves
    are never fused together:
      - Unions `cls_id == 1` (`Petiole`) into its largest 8-connected component (`area >= min_area_px`)
        and `cls_id == 2` (`Hole`) into all 8-connected components (`area >= min_area_px`).
      - Processes each `cls_id == 0` (`Leaf`) instance individually, retaining its primary
        8-connected component (`area >= min_area_px`) plus all interior unbridged lateral
        pinna islands (`area >= min_area_px`) and any component touching the segmented
        `petiole` corridor (`cls_id == 1`) without distance truncation (while rejecting
        border-clipped neighbor bleed), subtracting `hole_mask_u8` as
        `lamina = silhouette & ~hole` without adding any synthetic bridging pixels.
    """
    ch, cw = int(crop_shape[0]), int(crop_shape[1])
    if not raw_dets or ch <= 0 or cw <= 0:
        return []
    min_area = max(1, int(min_area_px))

    raw_leaf_dets: List[Tuple[float, np.ndarray]] = []
    union_by_cls: Dict[int, np.ndarray] = {}
    max_score_by_cls: Dict[int, float] = {}
    for _bbox, mask_arr, score, cls_id in raw_dets:
        cid = int(cls_id)
        if mask_arr is None or mask_arr.size == 0:
            continue
        m_u8 = (np.asarray(mask_arr) > 0).astype(np.uint8)
        if m_u8.shape[:2] != (ch, cw):
            m_u8 = cv2.resize(m_u8, (cw, ch), interpolation=cv2.INTER_NEAREST)
        if not np.any(m_u8):
            continue
        if cid == 0:
            raw_leaf_dets.append((float(score), m_u8))
        elif cid in (1, 2):
            if cid not in union_by_cls:
                union_by_cls[cid] = m_u8 > 0
                max_score_by_cls[cid] = float(score)
            else:
                union_by_cls[cid] |= (m_u8 > 0)
                if float(score) > max_score_by_cls[cid]:
                    max_score_by_cls[cid] = float(score)

    if not raw_leaf_dets and not union_by_cls:
        return []

    class_masks_u8: Dict[int, np.ndarray] = {}

    if 1 in union_by_cls:
        union_u8 = union_by_cls[1].astype(np.uint8)
        n_cc, labels, stats, _ = cv2.connectedComponentsWithStats(union_u8, connectivity=8)
        if n_cc > 1:
            best_lbl = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
            best_area = int(stats[best_lbl, cv2.CC_STAT_AREA])
            if best_area >= min_area:
                comp_u8 = (labels == best_lbl).astype(np.uint8)
                cnts, _ = cv2.findContours(comp_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                polys = [c.reshape(-1, 1, 2).astype(np.int32) for c in cnts if c.reshape(-1, 2).shape[0] >= 3]
                if polys:
                    raster_u8 = np.zeros((ch, cw), dtype=np.uint8)
                    cv2.fillPoly(raster_u8, polys, 255)
                    if int(np.count_nonzero(raster_u8)) >= min_area:
                        class_masks_u8[1] = raster_u8

    if 2 in union_by_cls:
        union_hole_u8 = union_by_cls[2].astype(np.uint8)
        n_cc, labels, stats, _ = cv2.connectedComponentsWithStats(union_hole_u8, connectivity=8)
        hole_polys: List[np.ndarray] = []
        for lbl in range(1, n_cc):
            if int(stats[lbl, cv2.CC_STAT_AREA]) < min_area:
                continue
            comp_u8 = (labels == lbl).astype(np.uint8)
            cnts, _ = cv2.findContours(comp_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            for c in cnts:
                if c.reshape(-1, 2).shape[0] >= 3:
                    hole_polys.append(c.reshape(-1, 1, 2).astype(np.int32))
        if hole_polys:
            hole_mask_u8 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.fillPoly(hole_mask_u8, hole_polys, 255)
            if int(np.count_nonzero(hole_mask_u8)) >= min_area:
                class_masks_u8[2] = hole_mask_u8

    pet_ref_u8 = class_masks_u8.get(1)
    pet_corridor_bool: Optional[np.ndarray] = None
    if pet_ref_u8 is not None and int(np.count_nonzero(pet_ref_u8)) > 0:
        kern_pet_adj = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
        pet_corridor_bool = cv2.dilate((pet_ref_u8 > 0).astype(np.uint8), kern_pet_adj, iterations=1) > 0
    hole_sub_u8 = class_masks_u8.get(2)
    cx_mid = 0.5 * float(cw)
    cy_mid = 0.5 * float(ch)
    crop_diag = max(1.0, math.hypot(float(cw), float(ch)))

    ranked_leaf_entries: List[Tuple[float, Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]] = []
    for l_score, l_u8 in raw_leaf_dets:
        n_cc, labels, stats, _ = cv2.connectedComponentsWithStats(l_u8, connectivity=8)
        if n_cc <= 1:
            continue
        best_lbl = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
        best_area = int(stats[best_lbl, cv2.CC_STAT_AREA])
        if best_area < min_area:
            continue
        keep_lbls: List[int] = [best_lbl]
        if n_cc > 2:
            for cid_sub in range(1, n_cc):
                if cid_sub == best_lbl:
                    continue
                sub_area = int(stats[cid_sub, cv2.CC_STAT_AREA])
                if sub_area < min_area:
                    continue
                sub_bool = labels == cid_sub
                sx0 = int(stats[cid_sub, cv2.CC_STAT_LEFT])
                sy0 = int(stats[cid_sub, cv2.CC_STAT_TOP])
                sx1 = sx0 + int(stats[cid_sub, cv2.CC_STAT_WIDTH])
                sy1 = sy0 + int(stats[cid_sub, cv2.CC_STAT_HEIGHT])
                touches_crop_edge = (sx0 <= 1) or (sy0 <= 1) or (sx1 >= cw - 1) or (sy1 >= ch - 1)
                touches_pet_corridor = bool(
                    pet_corridor_bool is not None and np.any(sub_bool & pet_corridor_bool)
                )
                if (not touches_crop_edge) or touches_pet_corridor:
                    keep_lbls.append(cid_sub)
        raster_u8 = (np.isin(labels, keep_lbls).astype(np.uint8)) * 255
        if hole_sub_u8 is not None:
            raster_u8[hole_sub_u8 > 0] = 0
            n_cc_h, labels_h, stats_h, _ = cv2.connectedComponentsWithStats((raster_u8 > 0).astype(np.uint8), connectivity=8)
            if n_cc_h <= 1:
                continue
            valid_h_lbls = [
                lbl_h for lbl_h in range(1, n_cc_h) if int(stats_h[lbl_h, cv2.CC_STAT_AREA]) >= min_area
            ]
            if not valid_h_lbls:
                continue
            raster_u8 = np.isin(labels_h, valid_h_lbls).astype(np.uint8) * 255
        fg_count = int(np.count_nonzero(raster_u8))
        if fg_count < min_area:
            continue
        ys, xs = np.where(raster_u8 > 0)
        if ys.size == 0 or xs.size == 0:
            continue
        ymin = int(ys.min())
        ymax = min(ch, int(ys.max()) + 1)
        xmin = int(xs.min())
        xmax = min(cw, int(xs.max()) + 1)
        if ymax <= ymin or xmax <= xmin:
            continue
        mx_l = float(np.mean(xs))
        my_l = float(np.mean(ys))
        dist_norm = math.hypot(mx_l - cx_mid, my_l - cy_mid) / (0.65 * crop_diag)
        cent_w = max(0.35, 1.0 - dist_norm)
        focal_rank = float(l_score) * cent_w * (float(fg_count) ** 0.35)
        ranked_leaf_entries.append((focal_rank, ((ymin, xmin, ymax, xmax), raster_u8, float(l_score), 0)))

    ranked_leaf_entries.sort(key=lambda item: item[0], reverse=True)
    reduced: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]] = [
        entry for _rank, entry in ranked_leaf_entries
    ]

    for cid in (1, 2):
        if cid not in class_masks_u8:
            continue
        m_u8 = class_masks_u8[cid]
        ys, xs = np.where(m_u8 > 0)
        if ys.size == 0 or xs.size == 0:
            continue
        ymin = int(ys.min())
        ymax = min(ch, int(ys.max()) + 1)
        xmin = int(xs.min())
        xmax = min(cw, int(xs.max()) + 1)
        if ymax <= ymin or xmax <= xmin:
            continue
        reduced.append(((ymin, xmin, ymax, xmax), m_u8, float(max_score_by_cls[cid]), int(cid)))

    return reduced
