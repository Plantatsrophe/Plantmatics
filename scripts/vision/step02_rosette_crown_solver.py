#!/usr/bin/env python3
"""
Step 02 weighted least-squares rosette caudex crown point solver
(`scripts/vision/step02_rosette_crown_solver.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from scripts.vision.step02_rosette_crown_pose import (
    _extract_xy_point,
    _parse_basal_leaf_pose,
)


def _estimate_single_plant_crown(
    poses: List[Dict[str, Any]],
    root_pts: List[Tuple[float, float]],
    other_plant_centers: Optional[Sequence[Tuple[float, float]]] = None,
) -> Tuple[float, float]:
    """
    Compute the caudex crown point `p_crown = (x_crown, y_crown)` for a single plant individual
    via weighted least-squares intersection of inward basal-leaf midrib/petiole rays anchored by
    LeafMachine3 `petiole_0..4` / `petiole_tip` distal endpoints, root detections, and/or the
    innermost `p_base` coordinates.
    """
    if not poses:
        if root_pts:
            rx = float(np.median([p[0] for p in root_pts]))
            ry = float(np.median([p[1] for p in root_pts]))
            return (rx, ry)
        return (0.0, 0.0)

    bases = np.array([p["p_base"] for p in poses], dtype=np.float64)
    apices = np.array([p["p_apex"] for p in poses], dtype=np.float64)
    vecs = bases - apices  # Inward midrib vector: apex -> base -> caudex crown
    lengths = np.hypot(vecs[:, 0], vecs[:, 1])

    # Incorporate aligned LM3 ray directions and distal petiole_tip endpoints when available
    petiole_tips: List[Tuple[float, float]] = []
    ray_weights = np.ones(len(poses), dtype=np.float64)
    ray_origins = bases.copy()
    for idx_p, p in enumerate(poses):
        if p.get("has_lm3_petiole_pose"):
            ray_weights[idx_p] = 2.5
        ptip = p.get("p_petiole_tip")
        if ptip is not None and len(ptip) == 2:
            petiole_tips.append((float(ptip[0]), float(ptip[1])))
            ray_origins[idx_p] = np.array([float(ptip[0]), float(ptip[1])], dtype=np.float64)
        rdir = p.get("ray_dir")
        if rdir is not None and len(rdir) == 2 and p.get("has_lm3_petiole_pose"):
            rnorm = math.hypot(float(rdir[0]), float(rdir[1]))
            if rnorm > 1e-6:
                eff_len = max(float(lengths[idx_p]), 25.0)
                vecs[idx_p] = np.array([float(rdir[0]) / rnorm, float(rdir[1]) / rnorm], dtype=np.float64) * eff_len
                lengths[idx_p] = eff_len

    valid_mask = lengths > 1e-4

    mean_base = np.median(bases, axis=0)
    mean_len = float(np.median(lengths[valid_mask])) if np.any(valid_mask) else 50.0

    # Filter root points to those associated with this plant individual (and closer to this plant than to any other plant)
    near_roots: List[Tuple[float, float]] = []
    if root_pts:
        max_root_dist = max(4.0 * mean_len, 450.0)
        for rp in root_pts:
            d_self = float(np.min(np.hypot(bases[:, 0] - rp[0], bases[:, 1] - rp[1])))
            d_center = math.hypot(rp[0] - float(mean_base[0]), rp[1] - float(mean_base[1]))
            if other_plant_centers:
                d_other = min(math.hypot(rp[0] - oc[0], rp[1] - oc[1]) for oc in other_plant_centers)
                if d_center > 1.05 * d_other and d_self > d_other:
                    continue
            if min(d_self, d_center) <= max_root_dist:
                near_roots.append(rp)
        if not near_roots and len(root_pts) == 1 and not other_plant_centers:
            near_roots = list(root_pts)

    has_petiole_tip_anchor = len(petiole_tips) > 0
    pet_tip_anchor = (
        np.array(
            [np.median([p[0] for p in petiole_tips]), np.median([p[1] for p in petiole_tips])],
            dtype=np.float64,
        )
        if has_petiole_tip_anchor
        else None
    )

    # Determine innermost p_base / petiole_6 / root coordinate anchor
    if near_roots:
        root_anchor = np.array(
            [np.median([p[0] for p in near_roots]), np.median([p[1] for p in near_roots])],
            dtype=np.float64,
        )
        dists_to_root = np.hypot(bases[:, 0] - root_anchor[0], bases[:, 1] - root_anchor[1])
        inner_idx = int(np.argmin(dists_to_root))
        innermost_base = bases[inner_idx]
        if pet_tip_anchor is not None:
            p_anchor = 0.50 * pet_tip_anchor + 0.35 * root_anchor + 0.15 * innermost_base
        else:
            p_anchor = 0.70 * root_anchor + 0.30 * innermost_base
        has_root_anchor = True
    elif pet_tip_anchor is not None:
        has_root_anchor = True
        innermost_base = pet_tip_anchor.copy()
        p_anchor = pet_tip_anchor.copy()
    else:
        has_root_anchor = False
        if np.any(valid_mask):
            unit_dirs = vecs[valid_mask] / lengths[valid_mask, None]
            mean_dir = np.sum(unit_dirs, axis=0)
            norm_dir = float(np.linalg.norm(mean_dir))
            if norm_dir > 0.25 and len(bases) >= 2:
                mean_dir = mean_dir / norm_dir
                proj_inward = (bases - mean_base) @ mean_dir
                lat_vec = np.array([-mean_dir[1], mean_dir[0]], dtype=np.float64)
                lat_dist = np.abs((bases - mean_base) @ lat_vec)
                scores = proj_inward - 0.35 * lat_dist
                k_inner = max(1, int(math.ceil(0.5 * len(bases))))
                top_idx = np.argsort(scores)[-k_inner:]
                innermost_base = np.median(bases[top_idx], axis=0)
            else:
                dists_center = np.hypot(bases[:, 0] - mean_base[0], bases[:, 1] - mean_base[1])
                innermost_base = bases[int(np.argmin(dists_center))]
        else:
            innermost_base = mean_base
        p_anchor = innermost_base.copy()

    if int(np.sum(valid_mask)) < 2:
        if has_root_anchor:
            if has_petiole_tip_anchor:
                return float(p_anchor[0]), float(p_anchor[1])
            if int(np.sum(valid_mask)) == 1:
                idx0 = int(np.where(valid_mask)[0][0])
                u0 = vecs[idx0] / lengths[idx0]
                b0 = ray_origins[idx0]
                t_root = float(np.dot(p_anchor - b0, u0))
                if t_root >= 0.0:
                    proj_pt = b0 + min(t_root, 2.5 * mean_len) * u0
                    crown_pt = 0.65 * proj_pt + 0.35 * p_anchor
                    return float(crown_pt[0]), float(crown_pt[1])
            return float(p_anchor[0]), float(p_anchor[1])
        return float(innermost_base[0]), float(innermost_base[1])

    # Weighted least-squares intersection of 2D inward midrib/petiole lines:
    # Each line i has normal projector P_i = w_i * (I - u_i u_i^T) passing through ray_origins[i].
    A_mat = np.zeros((2, 2), dtype=np.float64)
    b_vec = np.zeros(2, dtype=np.float64)
    valid_indices = np.where(valid_mask)[0]
    unit_list = []
    for idx_i in valid_indices:
        u_i = vecs[idx_i] / lengths[idx_i]
        w_i = float(ray_weights[idx_i])
        unit_list.append(u_i)
        P_i = w_i * (np.eye(2, dtype=np.float64) - np.outer(u_i, u_i))
        A_mat += P_i
        b_vec += P_i @ ray_origins[idx_i]

    eigvals = np.linalg.eigvalsh(A_mat)
    min_eig = float(np.min(eigvals))
    max_reach = max(2.5 * float(np.max(lengths[valid_mask])), 60.0)

    if min_eig > 1e-3:
        reg_w = (0.65 * float(len(petiole_tips))) if has_petiole_tip_anchor else (0.35 if has_root_anchor else 1e-4)
        p_ls = np.linalg.solve(A_mat + reg_w * np.eye(2, dtype=np.float64), b_vec + reg_w * p_anchor)
        u_arr = np.array(unit_list, dtype=np.float64)
        b_arr = bases[valid_indices]
        t_vals = np.sum((p_ls[None, :] - b_arr) * u_arr, axis=1)
        dist_from_inner = float(np.linalg.norm(p_ls - innermost_base))

        if float(np.median(t_vals)) >= -0.15 * mean_len and dist_from_inner <= max_reach:
            return float(p_ls[0]), float(p_ls[1])
        if dist_from_inner > max_reach and float(np.median(t_vals)) > 0.0:
            direction = (p_ls - innermost_base) / max(dist_from_inner, 1e-6)
            clamped = innermost_base + max_reach * direction
            if has_root_anchor:
                clamped = 0.5 * clamped + 0.5 * p_anchor
            return float(clamped[0]), float(clamped[1])

    # Fallback when rays are nearly parallel or divergent: anchor by roots / petiole_6 / innermost p_base
    return float(p_anchor[0]), float(p_anchor[1])


def estimate_rosette_crown_point(
    basal_leaves: Any = None,
    root_boxes: Optional[Sequence[Any]] = None,
    p_bases: Optional[Sequence[Tuple[float, float]]] = None,
    p_apices: Optional[Sequence[Tuple[float, float]]] = None,
    return_by_plant: bool = False,
    **kwargs: Any,
) -> Union[Tuple[float, float], Dict[int, Tuple[float, float]]]:
    """
    Compute the caudex crown point `p_crown` per `plant_individual_id` via weighted least-squares
    intersection of inward basal-leaf midrib/petiole rays (`p_apex -> p_base -> petiole_6`)
    anchored by LeafMachine3 `petiole_1..6` pose skeletons, `roots` detections, and/or the
    innermost `p_base` coordinates.

    Returns:
      - `Tuple[float, float]` `(x_crown, y_crown)` when invoked on a single plant's leaf sequence
        (with `return_by_plant=False`).
      - `Dict[int, Tuple[float, float]]` mapping `plant_individual_id -> (x_crown, y_crown)` when
        `basal_leaves` is a dict keyed by `plant_individual_id` or when `return_by_plant=True`.
    """
    raw_roots = (
        root_boxes
        if root_boxes is not None
        else (
            kwargs.get("root_detections")
            or kwargs.get("roots")
            or kwargs.get("root_pts")
            or kwargs.get("stage1_roots")
        )
    )

    if isinstance(basal_leaves, dict) and not (
        "p_base" in basal_leaves or "mask" in basal_leaves or "keypoints_31" in basal_leaves
    ):
        out_by_pid: Dict[int, Tuple[float, float]] = {}
        for pid_k, leaves_k in basal_leaves.items():
            pt_k = estimate_rosette_crown_point(
                leaves_k,
                root_boxes=raw_roots,
                return_by_plant=False,
            )
            out_by_pid[int(pid_k)] = (float(pt_k[0]), float(pt_k[1]))
        return out_by_pid

    poses: List[Dict[str, Any]] = []
    if p_bases is not None and p_apices is not None:
        for pb, pa in zip(p_bases, p_apices):
            parsed = _parse_basal_leaf_pose((pb, pa))
            if parsed is not None:
                poses.append(parsed)
    elif isinstance(basal_leaves, (dict, np.ndarray)):
        parsed = _parse_basal_leaf_pose(basal_leaves)
        if parsed is not None:
            poses.append(parsed)
    elif basal_leaves is not None:
        for item in basal_leaves:
            parsed = _parse_basal_leaf_pose(item)
            if parsed is not None:
                poses.append(parsed)

    ref_pt: Optional[Tuple[float, float]] = None
    if poses:
        ref_pt = (
            float(np.median([p["p_base"][0] for p in poses])),
            float(np.median([p["p_base"][1] for p in poses])),
        )

    root_pts: List[Tuple[float, float]] = []
    if raw_roots is not None:
        if isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 2 and all(
            isinstance(v, (int, float, np.number)) for v in raw_roots
        ):
            root_pts.append((float(raw_roots[0]), float(raw_roots[1])))
        elif isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 4 and all(
            isinstance(v, (int, float, np.number)) for v in raw_roots
        ):
            pt = _extract_xy_point(raw_roots, ref_pt=ref_pt)
            if pt is not None:
                root_pts.append(pt)
        else:
            for r_item in raw_roots:
                if isinstance(r_item, (tuple, list)) and len(r_item) >= 3 and isinstance(r_item[0], (tuple, list, np.ndarray)):
                    cls_id = int(r_item[-1]) if isinstance(r_item[-1], (int, float, np.integer)) else 7
                    if cls_id not in (7, 8, 9):
                        continue
                elif isinstance(r_item, dict) and "class_id" in r_item and int(r_item["class_id"]) not in (7, 8, 9):
                    continue
                pt = _extract_xy_point(r_item, ref_pt=ref_pt)
                if pt is not None:
                    root_pts.append(pt)

    if return_by_plant:
        by_pid: Dict[int, List[Dict[str, Any]]] = {}
        for p in poses:
            by_pid.setdefault(int(p["plant_individual_id"]), []).append(p)
        if not by_pid:
            return {0: _estimate_single_plant_crown([], root_pts)}
        pid_centers: Dict[int, Tuple[float, float]] = {
            pid_k: (
                float(np.median([p["p_base"][0] for p in pid_poses])),
                float(np.median([p["p_base"][1] for p in pid_poses])),
            )
            for pid_k, pid_poses in by_pid.items()
        }
        return {
            pid_k: _estimate_single_plant_crown(
                pid_poses,
                root_pts,
                other_plant_centers=[c for k, c in pid_centers.items() if k != pid_k],
            )
            for pid_k, pid_poses in sorted(by_pid.items())
        }

    return _estimate_single_plant_crown(poses, root_pts)
