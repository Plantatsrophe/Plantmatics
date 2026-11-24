#!/usr/bin/env python3
"""
Step 02 reusable leaf-zone geometry helpers: organ-coordinate extraction,
base/apex/blade/petiole resolution, per-PID organ scoping, vertical insertion
metrics, and archetype confidence scoring (`scripts/vision/step02_leaf_zone_geometry.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import cv2
import numpy as np

from scripts.vision.step02_rosette_crown_pose import _extract_xy_point, _parse_basal_leaf_pose


def _safe_finite(val: Any) -> Optional[float]:
    if val is None:
        return None
    try:
        fval = float(val)
        return fval if np.isfinite(fval) else None
    except (TypeError, ValueError):
        return None


def _safe_nonneg_float(val: Any) -> Optional[float]:
    fval = _safe_finite(val)
    return fval if (fval is not None and fval >= 0.0) else None


def _coerce_finite_xy(val: Any) -> Optional[Tuple[float, float]]:
    if val is not None and isinstance(val, (tuple, list, np.ndarray)) and len(val) >= 2:
        if np.isfinite(val[0]) and np.isfinite(val[1]):
            return (float(val[0]), float(val[1]))
    return None


def _parse_xy_map(raw_map: Optional[Mapping[int, Any]]) -> Dict[int, Tuple[float, float]]:
    out: Dict[int, Tuple[float, float]] = {}
    if raw_map:
        for k, v in raw_map.items():
            if v is not None and len(v) >= 2 and np.isfinite(v[0]) and np.isfinite(v[1]):
                out[int(k)] = (float(v[0]), float(v[1]))
    return out


def parse_zone_inputs(
    crowns_by_pid: Optional[Mapping[int, Tuple[float, float]]],
    rosette_centroids_by_pid: Optional[Mapping[int, Optional[Tuple[float, float]]]],
    basal_stats_by_pid: Optional[Mapping[int, Tuple[float, float, float, float]]],
    capitula_detections: Optional[Sequence[Any]],
    root_boxes: Optional[Sequence[Any]],
    specimen_mask: Optional[np.ndarray],
    p_crown: Optional[Tuple[float, float]],
    p_corymb: Optional[Tuple[float, float]],
    kwargs: Mapping[str, Any],
) -> Tuple[Any, ...]:
    """Normalize per-PID crown/rosette/stat maps, sheet-space root/capitula points, and specimen components."""
    crowns_map = _parse_xy_map(crowns_by_pid)
    explicit_crown = p_crown if p_crown is not None else kwargs.get("p_crown")
    if explicit_crown is not None and len(explicit_crown) >= 2 and np.isfinite(explicit_crown[0]) and np.isfinite(explicit_crown[1]):
        crowns_map.setdefault(0, (float(explicit_crown[0]), float(explicit_crown[1])))

    rosette_map = _parse_xy_map(rosette_centroids_by_pid)
    stats_map: Dict[int, Tuple[float, float, float, float]] = {}
    if basal_stats_by_pid:
        for k, v in basal_stats_by_pid.items():
            if v is not None and len(v) >= 4:
                stats_map[int(k)] = (float(v[0]), float(v[1]), float(v[2]), float(v[3]))

    raw_roots = (
        root_boxes if root_boxes is not None
        else (kwargs.get("root_detections") or kwargs.get("roots") or kwargs.get("root_pts") or kwargs.get("stage1_roots") or [])
    )
    root_pts_all: List[Tuple[float, float]] = []
    if isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 2 and all(isinstance(v, (int, float, np.number)) for v in raw_roots):
        root_pts_all.append((float(raw_roots[0]), float(raw_roots[1])))
    elif isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 4 and all(isinstance(v, (int, float, np.number)) for v in raw_roots):
        pt = _extract_xy_point(raw_roots)
        if pt is not None and np.isfinite(pt[0]) and np.isfinite(pt[1]):
            root_pts_all.append((float(pt[0]), float(pt[1])))
    elif raw_roots is not None:
        for r_item in raw_roots:
            if isinstance(r_item, (tuple, list)) and len(r_item) >= 3 and isinstance(r_item[0], (tuple, list, np.ndarray)):
                cls_id = int(r_item[-1]) if isinstance(r_item[-1], (int, float, np.integer)) else 7
                if cls_id not in (7, 8, 9):
                    continue
            elif isinstance(r_item, dict) and "class_id" in r_item and int(r_item["class_id"]) not in (7, 8, 9):
                continue
            pt = _extract_xy_point(r_item)
            if pt is not None and np.isfinite(pt[0]) and np.isfinite(pt[1]):
                root_pts_all.append((float(pt[0]), float(pt[1])))

    raw_caps = capitula_detections or kwargs.get("capitula_boxes") or kwargs.get("capitula_pts") or []
    cap_pts_all: List[Tuple[float, float]] = []
    for c_item in raw_caps:
        pt = _extract_xy_point(c_item)
        if pt is not None and np.isfinite(pt[0]) and np.isfinite(pt[1]):
            cap_pts_all.append((float(pt[0]), float(pt[1])))
    explicit_corymb = p_corymb if p_corymb is not None else kwargs.get("p_corymb")
    if explicit_corymb is not None and len(explicit_corymb) >= 2 and np.isfinite(explicit_corymb[0]) and np.isfinite(explicit_corymb[1]):
        cap_pts_all.append((float(explicit_corymb[0]), float(explicit_corymb[1])))

    if specimen_mask is None:
        specimen_mask = kwargs.get("specimen_paper_cleaned_mask")
    comp_labels: Optional[np.ndarray] = None
    spec_bin: Optional[np.ndarray] = None
    if isinstance(specimen_mask, np.ndarray) and specimen_mask.ndim == 2 and specimen_mask.size > 0:
        spec_bin = (specimen_mask > 0).astype(np.uint8)
        if np.count_nonzero(spec_bin) > 0:
            k_bridge = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9))
            _, comp_labels = cv2.connectedComponents(cv2.dilate(spec_bin, k_bridge, iterations=1), connectivity=8)

    return crowns_map, rosette_map, stats_map, root_pts_all, cap_pts_all, spec_bin, comp_labels


def comp_id_at(comp_labels: Optional[np.ndarray], pt: Tuple[float, float]) -> int:
    """Look up the dominant Specimen_Segmenter connected-component label at or within 6 px of `pt`."""
    if comp_labels is None:
        return 0
    h_c, w_c = comp_labels.shape[:2]
    ix, iy = int(round(float(pt[0]))), int(round(float(pt[1])))
    if 0 <= iy < h_c and 0 <= ix < w_c:
        cid = int(comp_labels[iy, ix])
        if cid > 0:
            return cid
    patch = comp_labels[max(0, iy - 6) : min(h_c, iy + 7), max(0, ix - 6) : min(w_c, ix + 7)]
    nz = patch[patch > 0]
    if nz.size > 0:
        vals, counts = np.unique(nz, return_counts=True)
        return int(vals[int(np.argmax(counts))])
    return 0


def resolve_leaf_base_and_apex(rec: Mapping[str, Any], lid: int) -> Dict[str, Any]:
    """Resolve sheet-space `p_base`, `p_apex`, `p_tip`, inward `ray_dir`, and parsed LM3 pose for a leaf record."""
    rec_d = dict(rec)
    rec_d.setdefault("leaf_id", lid)
    parsed = _parse_basal_leaf_pose(rec_d)
    p_base = parsed.get("p_base") if parsed is not None else None
    p_apex = parsed.get("p_apex") if parsed is not None else None
    p_tip = parsed.get("p_petiole_tip") if parsed is not None else None
    ray_dir = parsed.get("ray_dir") if parsed is not None else None
    has_pose_pet = bool(parsed.get("has_lm3_petiole_pose", False)) if parsed is not None else False
    pose_poly_px = float(parsed.get("pose_petiole_polyline_px", 0.0) or 0.0) if parsed is not None else 0.0

    for key in ("p_base_sheet", "p_base"):
        if p_base is None:
            p_base = _coerce_finite_xy(rec.get(key))
    for key in ("p_apex_sheet", "p_apex"):
        if p_apex is None:
            p_apex = _coerce_finite_xy(rec.get(key))
    if p_tip is None:
        p_tip = _coerce_finite_xy(rec.get("p_petiole_tip"))
    if ray_dir is None:
        cand_rd = _coerce_finite_xy(rec.get("ray_dir"))
        if cand_rd is not None:
            rn = math.hypot(cand_rd[0], cand_rd[1])
            if rn > 1e-6:
                ray_dir = (cand_rd[0] / rn, cand_rd[1] / rn)

    if p_base is None:
        cx_val, cy_val = rec.get("cx", rec.get("centroid_x")), rec.get("cy", rec.get("centroid_y"))
        if cx_val is not None and cy_val is not None and np.isfinite(float(cx_val)) and np.isfinite(float(cy_val)):
            p_base = (float(cx_val), float(cy_val))
        else:
            bbox = rec.get("bbox")
            if bbox is not None and len(bbox) == 4:
                ymin, xmin, ymax, xmax = [float(v) for v in bbox]
                p_base = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))
            else:
                p_base = (0.0, 0.0)

    if p_apex is None:
        p_apex = p_base
    if p_tip is None:
        p_tip = p_base

    if ray_dir is None:
        dx_tb, dy_tb = float(p_tip[0] - p_base[0]), float(p_tip[1] - p_base[1])
        norm_tb = math.hypot(dx_tb, dy_tb)
        if norm_tb > 3.0:
            ray_dir = (dx_tb / norm_tb, dy_tb / norm_tb)
        else:
            dx_ab, dy_ab = float(p_base[0] - p_apex[0]), float(p_base[1] - p_apex[1])
            norm_ab = math.hypot(dx_ab, dy_ab)
            ray_dir = (dx_ab / norm_ab, dy_ab / norm_ab) if norm_ab > 1e-6 else (0.0, 1.0)

    return {
        "p_base": p_base, "p_apex": p_apex, "p_tip": p_tip, "ray_dir": ray_dir,
        "has_pose_pet": has_pose_pet, "pose_poly_px": pose_poly_px, "parsed_pose": parsed,
    }


def resolve_leaf_blade_length_px(rec: Mapping[str, Any], geom: Mapping[str, Any]) -> float:
    """Resolve positive blade length in pixels from record metadata, base-apex chord, or bbox span."""
    for key in ("blade_length_px", "length_px"):
        fval = _safe_nonneg_float(rec.get(key))
        if fval is not None and fval > 0.0:
            return fval
    pb, pa = geom["p_base"], geom["p_apex"]
    d = float(math.hypot(float(pa[0]) - float(pb[0]), float(pa[1]) - float(pb[1])))
    if np.isfinite(d) and d > 0.0:
        return d
    bbox = rec.get("bbox")
    if bbox is not None and len(bbox) == 4:
        try:
            ymin, xmin, ymax, xmax = [float(v) for v in bbox]
            span = max(abs(ymax - ymin), abs(xmax - xmin))
            if np.isfinite(span) and span > 0.0:
                return float(span)
        except (TypeError, ValueError):
            pass
    return 50.0


def resolve_intrinsic_petiole_ratio(
    rec: Mapping[str, Any], geom: Mapping[str, Any], blade_px: float, dist_crown_px: float
) -> Tuple[float, bool]:
    """Evaluate leaf-intrinsic petiole ratio and topological petiole presence without conflating crown distance."""
    pet_mask = rec.get("petiole_mask")
    has_pet_mask = bool(rec.get("has_petiole_mask", False)) or (
        isinstance(pet_mask, np.ndarray) and np.count_nonzero(pet_mask) > 0
    )
    pose_poly = max(
        _safe_nonneg_float(rec.get("pose_petiole_polyline_px")) or 0.0,
        _safe_nonneg_float(geom.get("pose_poly_px")) or 0.0,
    )
    geo_skel = _safe_nonneg_float(rec.get("geodesic_petiole_skeleton_px")) or 0.0
    p_tip, p_base = geom["p_tip"], geom["p_base"]
    tip_base_px = float(math.hypot(float(p_tip[0]) - float(p_base[0]), float(p_tip[1]) - float(p_base[1])))

    direct_pet_px = max(pose_poly, geo_skel, tip_base_px if tip_base_px > 3.0 else 0.0)
    if direct_pet_px > 0.0:
        return (float(direct_pet_px / blade_px) if blade_px > 1e-6 else 0.0), True

    far_from_crown = dist_crown_px > 1.75 * max(blade_px, 1.0) and not has_pet_mask
    r_val = _safe_nonneg_float(rec.get("petiole_to_blade_ratio"))
    if r_val is not None and r_val > 0.0:
        if far_from_crown:
            crown_ratio_est = dist_crown_px / max(blade_px, 1.0)
            if abs(r_val - crown_ratio_est) <= 0.25 * max(crown_ratio_est, 1.0):
                return 0.0, False
        return r_val, True

    p_val = _safe_nonneg_float(rec.get("petiole_length_px"))
    if p_val is not None and p_val > 0.0:
        if far_from_crown and abs(p_val - dist_crown_px) <= 5.0:
            return 0.0, False
        return (float(p_val / blade_px) if blade_px > 1e-6 else 0.0), True

    pet_mm = _safe_nonneg_float(rec.get("petiole_length_mm"))
    if pet_mm is not None and pet_mm > 0.0:
        blade_mm = _safe_nonneg_float(rec.get("blade_length_mm"))
        if blade_mm is not None and blade_mm > 1e-6:
            r_mm = float(pet_mm / blade_mm)
            if far_from_crown and r_mm > 1.75:
                return 0.0, False
            return r_mm, True
        return 0.18, True

    return (0.18 if has_pet_mask else 0.0), has_pet_mask


def build_pid_geometry_context(
    records_by_pid: Mapping[int, List[Tuple[int, Mapping[str, Any]]]], comp_labels: Optional[np.ndarray]
) -> Dict[str, Dict[int, Any]]:
    """Precompute per-PID leaf geometries, blade lengths, anchor clouds, component sets, and spans."""
    pid_geoms: Dict[int, List[Dict[str, Any]]] = {}
    pid_blade_lengths: Dict[int, List[float]] = {}
    pid_anchor_pts: Dict[int, List[Tuple[float, float]]] = {}
    pid_comp_sets: Dict[int, Set[int]] = {}
    pid_med_blade: Dict[int, float] = {}
    pid_span: Dict[int, float] = {}

    for pid, pid_items in records_by_pid.items():
        g_list = [resolve_leaf_base_and_apex(rec, lid) for lid, rec in pid_items]
        b_unsorted = [resolve_leaf_blade_length_px(rec, g_list[i]) for i, (_, rec) in enumerate(pid_items)]
        b_lens = sorted(b_unsorted)
        med_b = (
            float(np.median(b_lens[len(b_lens) // 2 :])) if len(b_lens) >= 4
            else (float(np.median(b_lens)) if b_lens else 50.0)
        )
        pts_cloud: List[Tuple[float, float]] = []
        comp_set: Set[int] = set()
        for g in g_list:
            pts_cloud.extend((g["p_tip"], g["p_base"]))
            for pt_cand in (g["p_tip"], g["p_base"], g["p_apex"]):
                cid = comp_id_at(comp_labels, pt_cand)
                if cid > 0:
                    comp_set.add(cid)
        pid_geoms[pid], pid_blade_lengths[pid], pid_med_blade[pid] = g_list, b_unsorted, med_b
        pid_anchor_pts[pid], pid_comp_sets[pid] = pts_cloud, comp_set
        pid_span[pid] = (
            max(float(np.ptp([p[0] for p in pts_cloud])), float(np.ptp([p[1] for p in pts_cloud])), 2.5 * med_b, 160.0)
            if pts_cloud
            else max(2.5 * med_b, 160.0)
        )

    return {
        "pid_geoms": pid_geoms, "pid_blade_lengths": pid_blade_lengths,
        "pid_anchor_pts": pid_anchor_pts, "pid_comp_sets": pid_comp_sets,
        "pid_med_blade": pid_med_blade, "pid_span": pid_span,
    }


def scope_organ_pts_to_pid(
    organ_pts: Sequence[Tuple[float, float]], target_pid: int, records_by_pid: Mapping[int, Any],
    pid_ctx: Mapping[str, Mapping[int, Any]], comp_labels: Optional[np.ndarray],
    max_mult_span: float, max_mult_blade: float, min_abs_px: float,
) -> List[Tuple[float, float]]:
    """Restrict sheet-level root or capitula organ points to those belonging to `target_pid`."""
    if not organ_pts:
        return []
    if len(records_by_pid) == 1:
        return list(organ_pts)

    pid_anchor_pts, pid_comp_sets = pid_ctx["pid_anchor_pts"], pid_ctx["pid_comp_sets"]

    def _min_dist(pt: Tuple[float, float], pid_k: int) -> float:
        pts_t = pid_anchor_pts.get(pid_k, [])
        return min(math.hypot(pt[0] - p[0], pt[1] - p[1]) for p in pts_t) if pts_t else float("inf")

    scoped: List[Tuple[float, float]] = []
    max_allow = max(
        max_mult_span * pid_ctx["pid_span"].get(target_pid, 200.0),
        max_mult_blade * pid_ctx["pid_med_blade"].get(target_pid, 50.0),
        min_abs_px,
    )
    t_comps = pid_comp_sets.get(target_pid, set())
    other_pids = [k for k in records_by_pid.keys() if k != target_pid]

    for op in organ_pts:
        d_self = _min_dist(op, target_pid)
        d_other = min((_min_dist(op, ok) for ok in other_pids), default=float("inf"))
        op_cid = comp_id_at(comp_labels, op)
        if op_cid > 0 and t_comps:
            in_self_comp = op_cid in t_comps
            in_other_comp = any(op_cid in pid_comp_sets.get(ok, set()) for ok in other_pids)
            if in_other_comp and not in_self_comp:
                continue
            if in_self_comp and not in_other_comp and d_self <= 1.5 * max_allow:
                scoped.append(op)
                continue
        if d_self <= max_allow and d_self <= 1.05 * d_other:
            scoped.append(op)
    return scoped


def compute_pid_vertical_metrics(
    pid_items: Sequence[Tuple[int, Mapping[str, Any]]], g_list: Sequence[Mapping[str, Any]],
    roots_pid: Sequence[Tuple[float, float]], caps_pid: Sequence[Tuple[float, float]],
    c_pt: Tuple[float, float], plant_bbox: Optional[Tuple[float, float, float, float]], n_pids: int,
) -> List[float]:
    """Compute per-leaf relative vertical coordinate `cy_rel_plant_bbox` across the individual's vertical span."""
    ys_span_pid: List[float] = []
    leaf_cy_vals: List[float] = []
    for idx_y, (_, rec_y) in enumerate(pid_items):
        g_y = g_list[idx_y]
        pby, pay, pty = float(g_y["p_base"][1]), float(g_y["p_apex"][1]), float(g_y["p_tip"][1])
        ys_span_pid.extend((pby, pay, pty))
        bbox_y = rec_y.get("bbox")
        if bbox_y is not None and len(bbox_y) == 4:
            by0, by1 = _safe_finite(bbox_y[0]), _safe_finite(bbox_y[2])
            if by0 is not None and by1 is not None:
                ys_span_pid.extend((by0, by1))
        cy_cand = _safe_finite(rec_y.get("cy"))
        if cy_cand is None:
            cy_cand = _safe_finite(rec_y.get("centroid_y"))
        leaf_cy_vals.append(cy_cand if cy_cand is not None else 0.5 * (pby + pay))

    ys_span_pid.extend(float(rp_y[1]) for rp_y in roots_pid)
    ys_span_pid.extend(float(cp_y[1]) for cp_y in caps_pid)
    ys_span_pid.append(float(c_pt[1]))
    if plant_bbox is not None and len(plant_bbox) == 4 and n_pids == 1:
        pb0, pb2 = _safe_finite(plant_bbox[0]), _safe_finite(plant_bbox[2])
        if pb0 is not None and pb2 is not None:
            ys_span_pid.extend((pb0, pb2))

    ymin_pid = min(ys_span_pid) if ys_span_pid else 0.0
    ymax_pid = max(ys_span_pid) if ys_span_pid else 0.0
    span_y_pid = max(ymax_pid - ymin_pid, 1.0)

    cy_rel_by_leaf: List[float] = []
    for idx_y, (_, rec_y) in enumerate(pid_items):
        cy_rel_f = _safe_finite(rec_y.get("cy_rel_plant_bbox"))
        if cy_rel_f is not None:
            cy_rel_by_leaf.append(float(np.clip(cy_rel_f, 0.0, 1.0)))
        else:
            cy_rel_by_leaf.append(
                float(np.clip((leaf_cy_vals[idx_y] - ymin_pid) / span_y_pid, 0.0, 1.0))
                if (ymax_pid - ymin_pid) > 1e-6
                else 0.5
            )
    return cy_rel_by_leaf


def compute_leaf_archetype_score(rec: Mapping[str, Any], zone_conf: float) -> float:
    """Compute composite morphological/pose/zone archetype confidence score in `[0.0, 1.0]`."""
    bd_f = _safe_finite(rec.get("bilateral_dice"))
    bd_v = float(np.clip(bd_f, 0.0, 1.0)) if bd_f is not None else 0.85
    sia_f = _safe_finite(rec.get("si_a"))
    sia_v = abs(sia_f) if sia_f is not None else 0.05

    kps_raw = rec.get("keypoints_31")
    kp_confs: List[float] = []
    if isinstance(kps_raw, Mapping):
        for kp_v in kps_raw.values():
            c_val = kp_v.get("conf", kp_v.get("confidence")) if isinstance(kp_v, Mapping) else (
                kp_v[2] if isinstance(kp_v, (tuple, list, np.ndarray)) and len(kp_v) >= 3 else None
            )
            cf = _safe_finite(c_val)
            if cf is not None:
                kp_confs.append(float(np.clip(cf, 0.0, 1.0)))
    elif isinstance(kps_raw, np.ndarray) and kps_raw.ndim == 2 and kps_raw.shape[1] >= 3:
        col_c = kps_raw[:, 2].astype(np.float64, copy=False)
        fin_c = col_c[np.isfinite(col_c)]
        if fin_c.size > 0:
            kp_confs = [float(np.clip(c, 0.0, 1.0)) for c in fin_c]

    has_pet_fallback = bool(rec.get("has_petiole_mask")) or (
        isinstance(rec.get("petiole_mask"), np.ndarray) and np.count_nonzero(rec.get("petiole_mask")) > 0
    )
    kp_conf_mean = float(np.mean(kp_confs)) if kp_confs else (0.90 if has_pet_fallback else 0.75)
    raw_gates = rec.get("lm3_gates_pass", True)
    qc_ok = 1.0 if (True if raw_gates is None else bool(raw_gates)) else 0.40
    return float(
        np.clip((0.40 * bd_v + 0.20 * (1.0 - min(1.0, sia_v)) + 0.20 * kp_conf_mean + 0.20 * zone_conf) * qc_ok, 0.0, 1.0)
    )
