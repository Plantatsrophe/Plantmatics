#!/usr/bin/env python3
"""
Step 02 root-crown deduplication, leaf attachment resolution, and initial
DBSCAN / Specimen_Segmenter component assignment (`scripts/vision/step02_plant_cluster_attach.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import cv2
import numpy as np

from scripts.vision.lm3_compat import DBSCAN
from scripts.vision.step02_rosette_crown_pose import (
    _extract_xy_point,
    _parse_basal_leaf_pose,
)


def deduplicate_root_crowns(
    root_boxes: Optional[Sequence[Any]],
    max_dim: float,
    cal_ppm: float,
    has_cal_px: bool,
) -> Tuple[List[Tuple[float, float]], List[Optional[Tuple[float, float, float, float]]]]:
    """Extract and union-find deduplicate optional root crowns into 2D distinct centers and merged boxes."""
    raw_root_entries: List[Tuple[Tuple[float, float], Optional[Tuple[float, float, float, float]]]] = []
    if root_boxes:
        for rb in root_boxes:
            pt_r = _extract_xy_point(rb)
            if pt_r is None or not (np.isfinite(pt_r[0]) and np.isfinite(pt_r[1])):
                continue
            r_box: Optional[Tuple[float, float, float, float]] = None
            raw_b = (
                rb[0]
                if (isinstance(rb, (tuple, list)) and len(rb) >= 2 and isinstance(rb[0], (tuple, list, np.ndarray)))
                else rb
            )
            if isinstance(raw_b, (tuple, list, np.ndarray)) and len(raw_b) == 4:
                try:
                    by0, bx0, by1, bx1 = [float(v) for v in raw_b]
                    if all(np.isfinite(v) for v in (by0, bx0, by1, bx1)):
                        r_box = (min(by0, by1), min(bx0, bx1), max(by0, by1), max(bx0, bx1))
                except (TypeError, ValueError):
                    r_box = None
            raw_root_entries.append(((float(pt_r[0]), float(pt_r[1])), r_box))

    root_sep_thresh = max((45.0 * cal_ppm) if has_cal_px else 0.0, 0.16 * max_dim)
    distinct_roots: List[Tuple[float, float]] = []
    distinct_root_boxes: List[Optional[Tuple[float, float, float, float]]] = []
    if not raw_root_entries:
        return distinct_roots, distinct_root_boxes

    n_r = len(raw_root_entries)
    r_parent = list(range(n_r))

    def _r_find(x: int) -> int:
        while r_parent[x] != x:
            r_parent[x] = r_parent[r_parent[x]]
            x = r_parent[x]
        return x

    for i_r in range(n_r):
        for j_r in range(i_r + 1, n_r):
            d_r = math.hypot(
                raw_root_entries[i_r][0][0] - raw_root_entries[j_r][0][0],
                raw_root_entries[i_r][0][1] - raw_root_entries[j_r][0][1],
            )
            if d_r <= root_sep_thresh:
                ri_p, rj_p = _r_find(i_r), _r_find(j_r)
                if ri_p != rj_p:
                    r_parent[max(ri_p, rj_p)] = min(ri_p, rj_p)

    r_groups: Dict[int, List[int]] = {}
    for i_r in range(n_r):
        r_groups.setdefault(_r_find(i_r), []).append(i_r)
    for grp_key in sorted(r_groups.keys()):
        g_idxs = r_groups[grp_key]
        mx_r = float(np.mean([raw_root_entries[k][0][0] for k in g_idxs]))
        my_r = float(np.mean([raw_root_entries[k][0][1] for k in g_idxs]))
        distinct_roots.append((mx_r, my_r))
        g_boxes = [raw_root_entries[k][1] for k in g_idxs if raw_root_entries[k][1] is not None]
        if g_boxes:
            distinct_root_boxes.append((
                min(b[0] for b in g_boxes),
                min(b[1] for b in g_boxes),
                max(b[2] for b in g_boxes),
                max(b[3] for b in g_boxes),
            ))
        else:
            distinct_root_boxes.append(None)

    return distinct_roots, distinct_root_boxes


def resolve_leaf_attach_point(inst: Any, fallback_xy: Tuple[float, float]) -> Tuple[float, float]:
    """Resolve sheet-space petiole-tip or base attachment coordinate for a leaf instance."""
    fb = (float(fallback_xy[0]), float(fallback_xy[1]))
    meta = getattr(inst, "trim_meta", None)

    def _meta_get(key: str) -> Any:
        if meta is None:
            return None
        return meta.get(key) if isinstance(meta, Mapping) else getattr(meta, key, None)

    origin = getattr(inst, "mask_origin", None)
    ox, oy = (
        (float(origin[0]), float(origin[1]))
        if origin is not None and len(origin) >= 2
        else (0.0, 0.0)
    )
    in_sheet = bool(
        getattr(inst, "keypoints_in_sheet_space", False)
        or getattr(inst, "coords_in_sheet_space", False)
        or _meta_get("coords_in_sheet_space")
    )

    def _coerce_pt(raw: Any, force_sheet: bool = False) -> Optional[Tuple[float, float]]:
        if raw is None:
            return None
        pt: Optional[Tuple[float, float]] = None
        if (
            isinstance(raw, (tuple, list, np.ndarray))
            and len(raw) in (2, 3)
            and all(isinstance(v, (int, float, np.number)) for v in raw[:2])
        ):
            pt = (float(raw[0]), float(raw[1]))
        else:
            pt = _extract_xy_point(raw)
        if pt is None or not (np.isfinite(pt[0]) and np.isfinite(pt[1])):
            return None
        px, py = float(pt[0]), float(pt[1])
        if force_sheet or in_sheet or (ox == 0.0 and oy == 0.0):
            return (px, py)
        d_raw = math.hypot(px - fb[0], py - fb[1])
        d_off = math.hypot((px + ox) - fb[0], (py + oy) - fb[1])
        if d_off < d_raw:
            return (px + ox, py + oy)
        return (px, py)

    b_ymin, b_xmin, b_ymax, b_xmax = [float(v) for v in getattr(inst, "bbox", (fb[1], fb[0], fb[1], fb[0]))]
    pad_b = max(128.0, 2.0 * max(b_ymax - b_ymin, b_xmax - b_xmin))

    def _within_bbox_reach(pt_xy: Tuple[float, float]) -> bool:
        return (
            (b_xmin - pad_b) <= pt_xy[0] <= (b_xmax + pad_b)
            and (b_ymin - pad_b) <= pt_xy[1] <= (b_ymax + pad_b)
        )

    for raw_cand, f_sheet in (
        (getattr(inst, "p_petiole_tip_sheet", None), True),
        (_meta_get("p_petiole_tip_sheet"), True),
        (getattr(inst, "p_petiole_tip", None), False),
        (_meta_get("p_petiole_tip"), False),
        (getattr(inst, "p_base_sheet", None), True),
        (_meta_get("p_base_sheet"), True),
        (getattr(inst, "p_base", None), False),
        (_meta_get("p_base"), False),
    ):
        cand = _coerce_pt(raw_cand, force_sheet=f_sheet)
        if cand is not None and _within_bbox_reach(cand):
            return cand

    if getattr(inst, "keypoints_31", None) is not None or ox != 0.0 or oy != 0.0:
        parsed_l = _parse_basal_leaf_pose(inst)
        if parsed_l is not None:
            cand_p = parsed_l.get("p_petiole_tip") or parsed_l.get("p_base")
            if cand_p is not None and np.isfinite(cand_p[0]) and np.isfinite(cand_p[1]):
                pt_cand = (float(cand_p[0]), float(cand_p[1]))
                if _within_bbox_reach(pt_cand):
                    return pt_cand
    return fb


def assign_initial_dbscan_and_components(
    instances: List[Any],
    centroids: List[List[float]],
    coords: np.ndarray,
    sheet_width: int,
    sheet_height: int,
    max_dim: float,
    eps: float,
    specimen_mask: Optional[np.ndarray],
    has_spec_mask: bool,
    distinct_roots: List[Tuple[float, float]],
    distinct_root_boxes: List[Optional[Tuple[float, float, float, float]]],
) -> Tuple[List[int], List[Optional[Tuple[int, int]]], List[Optional[int]]]:
    """Run initial DBSCAN, map Specimen_Segmenter components, and split multi-root straddling clusters."""
    has_multi_roots = len(distinct_roots) >= 2

    def _leaf_attach(idx_leaf: int) -> Tuple[float, float]:
        return resolve_leaf_attach_point(instances[idx_leaf], (float(coords[idx_leaf, 0]), float(coords[idx_leaf, 1])))

    clustering = DBSCAN(eps=eps, min_samples=1).fit(coords)
    raw_labels = [int(lbl) for lbl in clustering.labels_]
    uniq_init = sorted(set(raw_labels))
    init_map = {lbl: idx for idx, lbl in enumerate(uniq_init)}
    inst_clusters = [init_map[lbl] for lbl in raw_labels]

    leaf_part_keys: List[Optional[Tuple[int, int]]] = [None] * len(instances)
    inst_root_owner: List[Optional[int]] = [None] * len(instances)

    if has_spec_mask and specimen_mask is not None:
        k_bridge = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
        bridged = cv2.dilate((specimen_mask > 0).astype(np.uint8), k_bridge, iterations=1)
        _, cc_map = cv2.connectedComponents(bridged, connectivity=8)

        roots_by_cid: Dict[int, List[int]] = {}
        if has_multi_roots:
            pad_r = max(16, int(round(0.03 * max_dim)))
            for r_idx, (rx, ry) in enumerate(distinct_roots):
                xi_r = int(np.clip(round(rx), 0, sheet_width - 1))
                yi_r = int(np.clip(round(ry), 0, sheet_height - 1))
                r_cid = int(cc_map[yi_r, xi_r])
                if r_cid == 0:
                    rb_box = distinct_root_boxes[r_idx]
                    if rb_box is not None:
                        ry0 = max(0, min(yi_r - pad_r, int(math.floor(rb_box[0])) - 8))
                        ry1 = min(sheet_height, max(yi_r + pad_r + 1, int(math.ceil(rb_box[2])) + 8))
                        rx0 = max(0, min(xi_r - pad_r, int(math.floor(rb_box[1])) - 8))
                        rx1 = min(sheet_width, max(xi_r + pad_r + 1, int(math.ceil(rb_box[3])) + 8))
                    else:
                        ry0 = max(0, yi_r - pad_r)
                        ry1 = min(sheet_height, yi_r + pad_r + 1)
                        rx0 = max(0, xi_r - pad_r)
                        rx1 = min(sheet_width, xi_r + pad_r + 1)
                    nz_r = cc_map[ry0:ry1, rx0:rx1]
                    nz_r = nz_r[nz_r > 0]
                    if nz_r.size > 0:
                        vals_r, cnts_r = np.unique(nz_r, return_counts=True)
                        r_cid = int(vals_r[int(np.argmax(cnts_r))])
                if r_cid > 0:
                    roots_by_cid.setdefault(r_cid, []).append(r_idx)

        for idx_i, (cx_i, cy_i) in enumerate(centroids):
            xi = int(np.clip(round(cx_i), 0, sheet_width - 1))
            yi = int(np.clip(round(cy_i), 0, sheet_height - 1))
            cid = int(cc_map[yi, xi])
            if cid == 0:
                b_ymin, b_xmin, b_ymax, b_xmax = instances[idx_i].bbox
                y0 = max(0, min(yi - 8, int(math.floor(b_ymin))))
                y1 = min(sheet_height, max(yi + 9, int(math.ceil(b_ymax))))
                x0 = max(0, min(xi - 8, int(math.floor(b_xmin))))
                x1 = min(sheet_width, max(xi + 9, int(math.ceil(b_xmax))))
                nz = cc_map[y0:y1, x0:x1]
                nz = nz[nz > 0]
                if nz.size > 0:
                    vals, cnts = np.unique(nz, return_counts=True)
                    cid = int(vals[int(np.argmax(cnts))])
            if cid > 0:
                comp_roots = roots_by_cid.get(cid, []) if has_multi_roots else []
                if len(comp_roots) >= 2:
                    ax_i, ay_i = _leaf_attach(idx_i)
                    best_r = min(
                        comp_roots,
                        key=lambda ri: math.hypot(ax_i - distinct_roots[ri][0], ay_i - distinct_roots[ri][1]),
                    )
                    leaf_part_keys[idx_i] = (cid, best_r)
                    inst_root_owner[idx_i] = best_r
                elif len(comp_roots) == 1:
                    leaf_part_keys[idx_i] = (cid, comp_roots[0])
                    inst_root_owner[idx_i] = comp_roots[0]
                else:
                    leaf_part_keys[idx_i] = (cid, -1)

    if has_multi_roots:
        for idx_i in range(len(instances)):
            if inst_root_owner[idx_i] is None:
                ax_i, ay_i = _leaf_attach(idx_i)
                dists_r = [math.hypot(ax_i - rx, ay_i - ry) for rx, ry in distinct_roots]
                best_ri = int(np.argmin(dists_r))
                if dists_r[best_ri] <= 0.28 * max_dim:
                    inst_root_owner[idx_i] = best_ri
        init_by_cl: Dict[int, List[int]] = {}
        for idx_i, cl_i in enumerate(inst_clusters):
            init_by_cl.setdefault(cl_i, []).append(idx_i)
        split_keys: List[Tuple[int, int]] = []
        for idx_i, cl_i in enumerate(inst_clusters):
            owners_in_cl = {inst_root_owner[j] for j in init_by_cl[cl_i] if inst_root_owner[j] is not None}
            if len(owners_in_cl) >= 2:
                own = inst_root_owner[idx_i]
                if own is None:
                    ax_i, ay_i = _leaf_attach(idx_i)
                    own = min(
                        owners_in_cl,
                        key=lambda ri: math.hypot(ax_i - distinct_roots[ri][0], ay_i - distinct_roots[ri][1]),
                    )
                    inst_root_owner[idx_i] = own
                split_keys.append((cl_i, int(own)))
            else:
                split_keys.append((cl_i, -1))
        uniq_split = sorted(set(split_keys))
        split_map = {sk: idx for idx, sk in enumerate(uniq_split)}
        inst_clusters = [split_map[sk] for sk in split_keys]

    return inst_clusters, leaf_part_keys, inst_root_owner
