#!/usr/bin/env python3
"""
Step 02 specimen-mask merging, bolting-scape shoot-column merging, single-leaf
rosette absorption, and `cluster_plant_individuals` entrypoint
(`scripts/vision/step02_plant_cluster_merge.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from scripts.vision.step02_plant_cluster_attach import (
    assign_initial_dbscan_and_components,
    deduplicate_root_crowns,
    resolve_leaf_attach_point,
)


def _merge_plant_clusters(
    instances: List[Any], coords: np.ndarray, inst_clusters: List[int],
    leaf_part_keys: List[Optional[Tuple[int, int]]], inst_root_owner: List[Optional[int]],
    distinct_roots: List[Tuple[float, float]], specimen_mask: Optional[np.ndarray],
    has_spec_mask: bool, sheet_width: int, sheet_height: int,
    max_dim: float, cal_ppm: float, has_cal_px: bool,
) -> None:
    """Execute 3-pass cluster merging and write normalized `plant_individual_id` onto `instances`."""
    n_cl = len(set(inst_clusters))
    sw, sh = float(max(sheet_width, 1)), float(max(sheet_height, 1))
    has_multi_roots = len(distinct_roots) >= 2
    parent = list(range(n_cl))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a: int, b: int) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    def _resolve_root_anchor(idxs_m: List[int], fallback_xy: Tuple[float, float]) -> Optional[int]:
        if not has_multi_roots:
            return None
        assigned = [inst_root_owner[j] for j in idxs_m if inst_root_owner[j] is not None]
        if assigned:
            vals_a, cnts_a = np.unique(assigned, return_counts=True)
            return int(vals_a[int(np.argmax(cnts_a))])
        dists = [math.hypot(fallback_xy[0] - rx, fallback_xy[1] - ry) for rx, ry in distinct_roots]
        best_i = int(np.argmin(dists))
        return best_i if dists[best_i] <= 0.28 * max_dim else None

    # 1. Merge leaf clusters sharing the same Specimen_Segmenter component (respecting multi-root Voronoi partitions)
    if has_spec_mask:
        comp_owner: Dict[Tuple[int, int], int] = {}
        for idx_i, p_key in enumerate(leaf_part_keys):
            if p_key is None:
                continue
            if p_key in comp_owner:
                _union(inst_clusters[idx_i], comp_owner[p_key])
            else:
                comp_owner[p_key] = inst_clusters[idx_i]

    # 2. Merge vertically aligned bolting-scape clusters in the same shoot column
    cl_members: Dict[int, List[int]] = {}
    for idx_i, cid_i in enumerate(inst_clusters):
        cl_members.setdefault(_find(cid_i), []).append(idx_i)

    cl_keys = sorted(cl_members.keys())
    cl_med_xy: Dict[int, Tuple[float, float]] = {
        ck: (float(np.median(coords[idxs, 0])), float(np.median(coords[idxs, 1])))
        for ck, idxs in cl_members.items()
    }
    cl_root_anchor = {ck: _resolve_root_anchor(cl_members[ck], cl_med_xy[ck]) for ck in cl_keys}
    diag_s = math.hypot(sw, sh)

    for i_pos, ck_a in enumerate(cl_keys):
        for ck_b in cl_keys[i_pos + 1 :]:
            if _find(ck_a) == _find(ck_b):
                continue
            ra_idx, rb_idx = cl_root_anchor.get(ck_a), cl_root_anchor.get(ck_b)
            if ra_idx is not None and rb_idx is not None and ra_idx != rb_idx:
                continue
            xa, ya = cl_med_xy[ck_a]
            xb, yb = cl_med_xy[ck_b]
            pts_a, pts_b = coords[cl_members[ck_a]], coords[cl_members[ck_b]]
            d_min_ab = float(np.min(np.hypot(pts_a[:, 0:1] - pts_b[:, 0:1].T, pts_a[:, 1:2] - pts_b[:, 1:2].T)))
            has_bridge = False
            if has_spec_mask and specimen_mask is not None:
                ia_x, ia_y = int(np.clip(round(xa), 0, sheet_width - 1)), int(np.clip(round(ya), 0, sheet_height - 1))
                ib_x, ib_y = int(np.clip(round(xb), 0, sheet_width - 1)), int(np.clip(round(yb), 0, sheet_height - 1))
                n_chk = 16
                hits_sb = sum(
                    1
                    for s_i in range(n_chk + 1)
                    if specimen_mask[
                        int(np.clip(round(ia_y + (ib_y - ia_y) * (s_i / n_chk)), 0, sheet_height - 1)),
                        int(np.clip(round(ia_x + (ib_x - ia_x) * (s_i / n_chk)), 0, sheet_width - 1)),
                    ] > 0
                )
                has_bridge = hits_sb >= int(0.65 * (n_chk + 1))
            v_gap_close = abs(ya - yb) <= 0.14 * sh or d_min_ab <= 0.22 * diag_s
            if has_spec_mask:
                if abs(xa - xb) <= 0.20 * sw and (v_gap_close or (has_bridge and abs(ya - yb) <= 0.55 * sh)):
                    _union(ck_a, ck_b)
            elif abs(xa - xb) <= 0.12 * sw and abs(ya - yb) <= 0.60 * sh:
                _union(ck_a, ck_b)

    # 3. Absorb isolated 1-leaf outliers within rosette/petiole reach of a multi-leaf cluster
    merged_members: Dict[int, List[int]] = {}
    for idx_i, cid_i in enumerate(inst_clusters):
        merged_members.setdefault(_find(cid_i), []).append(idx_i)
    multi_clusters = [ck for ck, idxs in merged_members.items() if len(idxs) >= 2]
    if multi_clusters:
        absorb_rad = (35.0 * cal_ppm) if has_cal_px else (0.12 * max_dim)
        attach_coords = np.array(
            [resolve_leaf_attach_point(instances[i], (float(coords[i, 0]), float(coords[i, 1]))) for i in range(len(instances))],
            dtype=np.float32,
        )
        merged_root_anchor = {
            ck: _resolve_root_anchor(idxs_m, (float(np.median(coords[idxs_m, 0])), float(np.median(coords[idxs_m, 1]))))
            for ck, idxs_m in merged_members.items()
        }

        def _has_mask_bridge(p1: Sequence[float], p2: Sequence[float]) -> bool:
            if not has_spec_mask or specimen_mask is None:
                return False
            n_chk = 16
            x1, y1, x2, y2 = float(p1[0]), float(p1[1]), float(p2[0]), float(p2[1])
            hits = 0
            for s_i in range(n_chk + 1):
                frac = float(s_i) / float(n_chk)
                sx = int(np.clip(round(x1 + (x2 - x1) * frac), 0, sheet_width - 1))
                sy = int(np.clip(round(y1 + (y2 - y1) * frac), 0, sheet_height - 1))
                if np.any(specimen_mask[max(0, sy - 1) : min(sheet_height, sy + 2), max(0, sx - 1) : min(sheet_width, sx + 2)] > 0):
                    hits += 1
            return (float(hits) / float(n_chk + 1)) >= 0.70

        for ck_s, idxs_s in list(merged_members.items()):
            if len(idxs_s) != 1:
                continue
            ra_s, pt_s = merged_root_anchor.get(ck_s), attach_coords[idxs_s[0]]
            best_mc: Optional[int] = None
            best_attach_d = float("inf")
            best_nearest_pt: Optional[np.ndarray] = None
            for ck_m in multi_clusters:
                ra_m = merged_root_anchor.get(ck_m)
                if ra_s is not None and ra_m is not None and ra_s != ra_m:
                    continue
                pts_m = attach_coords[merged_members[ck_m]]
                dists_m = np.hypot(pts_m[:, 0] - pt_s[0], pts_m[:, 1] - pt_s[1])
                min_idx = int(np.argmin(dists_m))
                d_m = float(dists_m[min_idx])
                if d_m < best_attach_d:
                    best_attach_d, best_mc, best_nearest_pt = d_m, ck_m, pts_m[min_idx]
            has_bridge_to_best = (
                _has_mask_bridge(pt_s, best_nearest_pt)
                if (best_mc is not None and best_nearest_pt is not None and has_spec_mask)
                else False
            )
            if best_mc is not None and (has_bridge_to_best or best_attach_d <= absorb_rad):
                _union(ck_s, best_mc)

    ordered_roots: List[int] = []
    root_to_pid: Dict[int, int] = {}
    for idx_i, cid_i in enumerate(inst_clusters):
        r_id = _find(cid_i)
        if r_id not in root_to_pid:
            root_to_pid[r_id] = len(ordered_roots)
            ordered_roots.append(r_id)
        instances[idx_i].plant_individual_id = root_to_pid[r_id]


def cluster_plant_individuals(
    instances: List[Any],
    sheet_width: int,
    sheet_height: int,
    eps_ratio: float = 0.15,
    specimen_mask: Optional[np.ndarray] = None,
    root_boxes: Optional[Sequence[Any]] = None,
    pixels_per_mm: Optional[float] = None,
) -> List[Any]:
    """
    Stem-aware spatial clustering across detected leaves on a herbarium sheet.
    Starts from DBSCAN centroid neighbourhoods and merges clusters connected by:
      1. Shared Specimen_Segmenter foreground plant silhouette components,
      2. Vertical bolting-scape alignment within a single shoot column (unless
         anchored to distinct horizontally separated root crowns), or
      3. Splayed single-leaf outliers within rosette reach of a multi-leaf plant.
    """
    if not instances:
        return instances
    if len(instances) == 1:
        instances[0].plant_individual_id = 0
        return instances

    centroids = [[(inst.bbox[1] + inst.bbox[3]) / 2.0, (inst.bbox[0] + inst.bbox[2]) / 2.0] for inst in instances]
    coords = np.array(centroids, dtype=np.float32)
    max_dim = float(max(sheet_width, sheet_height, 1))
    eps = eps_ratio * max_dim

    has_cal_px, cal_ppm = False, 0.0
    if pixels_per_mm is not None:
        try:
            cal_ppm = float(pixels_per_mm)
            has_cal_px = bool(np.isfinite(cal_ppm) and cal_ppm > 0.0)
        except (TypeError, ValueError):
            has_cal_px = False

    has_spec_mask = bool(
        isinstance(specimen_mask, np.ndarray)
        and specimen_mask.ndim == 2
        and specimen_mask.shape == (sheet_height, sheet_width)
        and np.any(specimen_mask > 0)
    )

    distinct_roots, distinct_root_boxes = deduplicate_root_crowns(root_boxes, max_dim, cal_ppm, has_cal_px)
    inst_clusters, leaf_part_keys, inst_root_owner = assign_initial_dbscan_and_components(
        instances, centroids, coords, sheet_width, sheet_height, max_dim, eps,
        specimen_mask, has_spec_mask, distinct_roots, distinct_root_boxes,
    )

    if len(set(inst_clusters)) > 1:
        _merge_plant_clusters(
            instances, coords, inst_clusters, leaf_part_keys, inst_root_owner, distinct_roots,
            specimen_mask, has_spec_mask, sheet_width, sheet_height, max_dim, cal_ppm, has_cal_px,
        )
        return instances

    for inst, cid_i in zip(instances, inst_clusters):
        inst.plant_individual_id = cid_i
    return instances
