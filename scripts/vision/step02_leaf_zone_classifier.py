#!/usr/bin/env python3
"""
Step 02 anatomical leaf-zone classifier (`basal` vs. `cauline`): stem-polyline
projection, inward petiole ray tracing, positional scoring, and
`classify_leaf_zone_anatomical` entrypoint (`scripts/vision/step02_leaf_zone_classifier.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from scripts.vision.step02_leaf_zone_geometry import (
    build_pid_geometry_context,
    compute_leaf_archetype_score,
    compute_pid_vertical_metrics,
    parse_zone_inputs,
    resolve_intrinsic_petiole_ratio,
    scope_organ_pts_to_pid,
)
from scripts.vision.step02_rosette_crown_solver import _estimate_single_plant_crown


def project_along_stem_polyline(
    pt: Tuple[float, float], p0: Tuple[float, float], p1: Tuple[float, float], pb_mid: Optional[Tuple[float, float]] = None
) -> float:
    """Project `pt` onto the straight or V-bent piecewise stem axis `p0 -> (pb_mid) -> p1` in `[0.0, 1.0]`."""
    if pb_mid is not None:
        seg1_dx, seg1_dy = float(pb_mid[0] - p0[0]), float(pb_mid[1] - p0[1])
        seg2_dx, seg2_dy = float(p1[0] - pb_mid[0]), float(p1[1] - pb_mid[1])
        l1, l2 = math.hypot(seg1_dx, seg1_dy), math.hypot(seg2_dx, seg2_dy)
        tot_l = l1 + l2
        if tot_l > 25.0 and l1 > 10.0 and l2 > 10.0:
            t1 = max(0.0, min(1.0, ((pt[0] - p0[0]) * seg1_dx + (pt[1] - p0[1]) * seg1_dy) / (l1 * l1)))
            d1 = math.hypot(pt[0] - (p0[0] + t1 * seg1_dx), pt[1] - (p0[1] + t1 * seg1_dy))
            t2 = max(0.0, min(1.0, ((pt[0] - pb_mid[0]) * seg2_dx + (pt[1] - pb_mid[1]) * seg2_dy) / (l2 * l2)))
            d2 = math.hypot(pt[0] - (pb_mid[0] + t2 * seg2_dx), pt[1] - (pb_mid[1] + t2 * seg2_dy))
            s_poly = (t1 * l1) / tot_l if d1 <= d2 else (l1 + t2 * l2) / tot_l
            return max(0.0, min(1.0, float(s_poly)))

    dx_s, dy_s = float(p1[0] - p0[0]), float(p1[1] - p0[1])
    l_s = math.hypot(dx_s, dy_s)
    if l_s <= 20.0:
        return 0.0
    return max(0.0, min(1.0, float(((pt[0] - p0[0]) * dx_s + (pt[1] - p0[1]) * dy_s) / (l_s * l_s))))


def _resolve_pid_poles_and_hub_radius(
    pid: int, g_list: Sequence[Mapping[str, Any]], roots_pid: List[Tuple[float, float]],
    caps_pid: List[Tuple[float, float]], other_centers: List[Tuple[float, float]], med_blade_px: float,
    crowns_map: Mapping[int, Tuple[float, float]], rosette_map: Mapping[int, Tuple[float, float]],
    stats_map: Mapping[int, Tuple[float, float, float, float]],
) -> Tuple[Tuple[float, float], Optional[Tuple[float, float]], float]:
    """Resolve basal pole `c_pt`, optional apical corymb pole `cap_pt`, and basal hub radius `r_hub`."""
    has_roots, has_caps = len(roots_pid) > 0, len(caps_pid) > 0
    poses_pid = [g["parsed_pose"] for g in g_list if g.get("parsed_pose") is not None]
    p_hub = (
        _estimate_single_plant_crown(poses_pid, roots_pid, other_plant_centers=other_centers)
        if poses_pid
        else ((float(np.median([g["p_tip"][0] for g in g_list])), float(np.median([g["p_tip"][1] for g in g_list]))) if g_list else (0.0, 0.0))
    )
    cap_pt = (float(np.median([p[0] for p in caps_pid])), float(np.median([p[1] for p in caps_pid]))) if has_caps else None

    if has_roots:
        root_med = (float(np.median([p[0] for p in roots_pid])), float(np.median([p[1] for p in roots_pid])))
        if pid in crowns_map and math.hypot(crowns_map[pid][0] - root_med[0], crowns_map[pid][1] - root_med[1]) <= max(2.0 * med_blade_px, 240.0):
            c_pt = crowns_map[pid]
        else:
            c_pt = _estimate_single_plant_crown(poses_pid, roots_pid, other_plant_centers=other_centers) if poses_pid else root_med
    elif has_caps and cap_pt is not None:
        cand_c = crowns_map.get(pid) or rosette_map.get(pid) or (stats_map[pid][:2] if pid in stats_map else None)
        dists_to_cap = [math.hypot(g["p_tip"][0] - cap_pt[0], g["p_tip"][1] - cap_pt[1]) for g in g_list]
        max_d_cap = max(dists_to_cap) if dists_to_cap else 0.0
        far_poses = [g["parsed_pose"] for i_g, g in enumerate(g_list) if g.get("parsed_pose") is not None and dists_to_cap[i_g] >= 0.55 * max_d_cap]
        far_hub = _estimate_single_plant_crown(far_poses, [], other_plant_centers=other_centers) if far_poses else p_hub
        c_pt = (float(cand_c[0]), float(cand_c[1])) if (cand_c is not None and math.hypot(cand_c[0] - cap_pt[0], cand_c[1] - cap_pt[1]) >= 0.45 * max_d_cap) else far_hub
    else:
        cand_c = crowns_map.get(pid)
        c_pt = (float(cand_c[0]), float(cand_c[1])) if (cand_c is not None and math.hypot(cand_c[0] - p_hub[0], cand_c[1] - p_hub[1]) <= max(1.8 * med_blade_px, 200.0)) else p_hub

    stat_rad = 0.0
    if pid in stats_map:
        st_tup = stats_map[pid]
        if len(st_tup) >= 4 and np.isfinite(st_tup[3]) and float(st_tup[3]) > 10.0:
            stat_rad = float(st_tup[3])
        elif len(st_tup) == 3 and np.isfinite(st_tup[2]) and 10.0 < float(st_tup[2]) <= 2500.0:
            stat_rad = float(st_tup[2])

    dists_tip_c = [float(math.hypot(g["p_tip"][0] - c_pt[0], g["p_tip"][1] - c_pt[1])) for g in g_list]
    low_dists = sorted(dists_tip_c)[: max(1, min(3, int(math.ceil(0.30 * len(dists_tip_c)))))]
    est_rad = float(np.median(low_dists)) * 1.20 if low_dists else 0.0
    if stat_rad > 10.0:
        rosette_rad = max(min(stat_rad, max(est_rad * 1.35, 2.35 * med_blade_px, 125.0)) if has_caps else stat_rad, 1.15 * med_blade_px, 45.0)
    else:
        rosette_rad = max(est_rad, 1.25 * med_blade_px, 45.0)
    r_hub = min(max(rosette_rad, 1.25 * med_blade_px, 50.0), max(2.35 * med_blade_px, 125.0)) if has_caps else max(rosette_rad, 1.35 * med_blade_px, 55.0)
    return c_pt, cap_pt, r_hub


def _trace_pid_inward_rays(
    pid_items: Sequence[Tuple[int, Mapping[str, Any]]], g_list: Sequence[Mapping[str, Any]],
    blade_lengths: Sequence[float], med_blade_px: float, plant_span: float,
    c_pt: Tuple[float, float], cap_pt: Optional[Tuple[float, float]], r_hub: float,
    has_roots: bool, has_caps: bool, pid_has_explicit_crown: bool,
    spec_bin: Optional[np.ndarray], is_bent_stem: bool, explicit_p_bend: Optional[Tuple[float, float]],
) -> Dict[str, Any]:
    """Trace leaf inward petiole rays toward `c_pt`, evaluate tuft convergence, PCA stem axis, and V-bend arch."""
    traced_att_pts: List[Tuple[float, float]] = []
    ray_crown_miss_px: List[float] = []
    ray_converges_to_crown: List[bool] = []
    intrinsic_has_pet: List[bool] = []
    effective_pet_ratios: List[float] = []
    effective_has_pet: List[bool] = []

    prelim_scape_dx = float(cap_pt[0] - c_pt[0]) if cap_pt is not None else 0.0
    prelim_scape_dy = float(cap_pt[1] - c_pt[1]) if cap_pt is not None else 0.0
    prelim_scape_len = math.hypot(prelim_scape_dx, prelim_scape_dy)

    for idx_i, (_, rec_i) in enumerate(pid_items):
        g_i, blade_px_i = g_list[idx_i], blade_lengths[idx_i]
        p_base_i, p_tip_i, u_in_i = g_i["p_base"], g_i["p_tip"], g_i["ray_dir"]

        raw_d_crown = math.hypot(p_tip_i[0] - c_pt[0], p_tip_i[1] - c_pt[1])
        pet_ratio_i, has_pet_i = resolve_intrinsic_petiole_ratio(rec_i, g_i, blade_px_i, raw_d_crown)
        intrinsic_has_pet.append(has_pet_i)

        max_ray_reach = 3.8 * max(blade_px_i, med_blade_px, 50.0)
        t_proj = float(c_pt[0] - p_tip_i[0]) * u_in_i[0] + float(c_pt[1] - p_tip_i[1]) * u_in_i[1]
        t_clamped = max(0.0, min(max_ray_reach, t_proj))
        ray_closest_pt = (p_tip_i[0] + t_clamped * u_in_i[0], p_tip_i[1] + t_clamped * u_in_i[1])
        d_perp_crown = math.hypot(ray_closest_pt[0] - c_pt[0], ray_closest_pt[1] - c_pt[1])

        mask_ray_ok = True
        if spec_bin is not None and t_clamped > 12.0:
            h_s, w_s = spec_bin.shape[:2]
            n_steps = max(4, min(24, int(round(t_clamped / 12.0))))
            hits, valid_steps = 0, 0
            for s_idx in range(1, n_steps + 1):
                frac = float(s_idx) / float(n_steps)
                sx, sy = int(round(p_tip_i[0] + frac * t_clamped * u_in_i[0])), int(round(p_tip_i[1] + frac * t_clamped * u_in_i[1]))
                if 0 <= sy < h_s and 0 <= sx < w_s:
                    valid_steps += 1
                    if np.any(spec_bin[max(0, sy - 5) : min(h_s, sy + 6), max(0, sx - 5) : min(w_s, sx + 6)] > 0):
                        hits += 1
            if valid_steps >= 4 and (float(hits) / float(valid_steps)) < 0.35:
                mask_ray_ok = False

        is_aerial_stem_node = False
        if has_caps and prelim_scape_len > 60.0 and raw_d_crown > 1.35 * r_hub:
            s_prelim = ((p_tip_i[0] - c_pt[0]) * prelim_scape_dx + (p_tip_i[1] - c_pt[1]) * prelim_scape_dy) / (prelim_scape_len * prelim_scape_len)
            if s_prelim > (0.22 if has_pet_i else 0.16):
                is_aerial_stem_node = True

        converges = (
            raw_d_crown <= r_hub
            or (
                t_proj >= -0.15 * blade_px_i
                and raw_d_crown <= max_ray_reach
                and d_perp_crown <= 1.05 * r_hub
                and mask_ray_ok
                and not is_aerial_stem_node
                and (
                    not has_caps
                    or (has_pet_i and pet_ratio_i >= 0.30 and raw_d_crown <= 1.95 * r_hub)
                    or (has_pet_i and raw_d_crown <= 1.55 * r_hub)
                    or raw_d_crown <= 1.25 * r_hub
                )
            )
        )

        if (has_roots or pid_has_explicit_crown) and converges and t_clamped > 0.0 and d_perp_crown < raw_d_crown:
            p_att_i = ray_closest_pt
            traced_pet_px = math.hypot(p_base_i[0] - p_att_i[0], p_base_i[1] - p_att_i[1])
            if blade_px_i > 1e-6 and traced_pet_px > pet_ratio_i * blade_px_i:
                pet_ratio_i = float(traced_pet_px / blade_px_i)
                if pet_ratio_i >= 0.15:
                    has_pet_i = True
        else:
            p_att_i = p_tip_i

        traced_att_pts.append(p_att_i)
        ray_crown_miss_px.append(d_perp_crown)
        ray_converges_to_crown.append(converges)
        effective_pet_ratios.append(pet_ratio_i)
        effective_has_pet.append(has_pet_i)

    n_leaves_pid = len(pid_items)
    conv_pet_vals = [effective_pet_ratios[k] for k in range(n_leaves_pid) if ray_converges_to_crown[k] and effective_pet_ratios[k] > 0.0]
    if n_leaves_pid >= 4 and len(conv_pet_vals) >= max(3, n_leaves_pid - 1) and float(np.median(conv_pet_vals)) >= 0.75:
        med_conv_pet = float(np.median(conv_pet_vals))
        for k in range(n_leaves_pid):
            if not ray_converges_to_crown[k] and blade_lengths[k] >= 0.75 * med_blade_px and (ray_crown_miss_px[k] <= 0.85 * r_hub or med_conv_pet >= 0.85):
                ray_converges_to_crown[k] = True
                traced_att_pts[k] = c_pt
                traced_pet_k = math.hypot(g_list[k]["p_base"][0] - c_pt[0], g_list[k]["p_base"][1] - c_pt[1])
                if blade_lengths[k] > 1e-6:
                    effective_pet_ratios[k] = max(effective_pet_ratios[k], float(traced_pet_k / blade_lengths[k]))
                    effective_has_pet[k] = True

    att_span_px, att_elongation = 0.0, 1.0
    stem_axis_vec: Tuple[float, float] = (0.0, -1.0)
    if n_leaves_pid >= 2:
        att_arr = np.asarray(traced_att_pts, dtype=np.float64)
        centered = att_arr - np.mean(att_arr, axis=0, keepdims=True)
        eigvals_pca, eigvecs_pca = np.linalg.eigh((centered.T @ centered) / float(max(n_leaves_pid - 1, 1)))
        att_elongation = math.sqrt(float(max(eigvals_pca[-1], 0.0)) / max(float(max(eigvals_pca[0], 1.0)), 1.0))
        axis_u = eigvecs_pca[:, -1]
        att_span_px = float(np.ptp(centered @ axis_u))
        out_vecs = np.asarray([(-g["ray_dir"][0], -g["ray_dir"][1]) for g in g_list], dtype=np.float64)
        if float(np.dot(np.mean(out_vecs, axis=0), axis_u)) < 0.0:
            axis_u = -axis_u
        stem_axis_vec = (float(axis_u[0]), float(axis_u[1]))

    p_bend = explicit_p_bend
    if p_bend is None and cap_pt is not None and n_leaves_pid >= 2:
        chord_dx, chord_dy = float(cap_pt[0] - c_pt[0]), float(cap_pt[1] - c_pt[1])
        chord_len = math.hypot(chord_dx, chord_dy)
        if is_bent_stem or chord_len < 0.55 * plant_span:
            best_arch_d, best_arch_pt = 0.0, None
            for pt_a in traced_att_pts:
                if chord_len > 20.0:
                    ux_c, uy_c = chord_dx / chord_len, chord_dy / chord_len
                    cl_c = max(0.0, min(chord_len, (pt_a[0] - c_pt[0]) * ux_c + (pt_a[1] - c_pt[1]) * uy_c))
                    d_arch = math.hypot(pt_a[0] - (c_pt[0] + cl_c * ux_c), pt_a[1] - (c_pt[1] + cl_c * uy_c))
                else:
                    d_arch = math.hypot(pt_a[0] - c_pt[0], pt_a[1] - c_pt[1])
                if d_arch > best_arch_d:
                    best_arch_d, best_arch_pt = d_arch, pt_a
            if best_arch_pt is not None and best_arch_d >= max(0.32 * plant_span, 1.5 * r_hub):
                p_bend = best_arch_pt

    dists_att_crown = [float(math.hypot(pt_a[0] - c_pt[0], pt_a[1] - c_pt[1])) for pt_a in traced_att_pts]
    n_converging = sum(1 for c_flag in ray_converges_to_crown if c_flag)
    return {
        "traced_att_pts": traced_att_pts, "dists_att_crown": dists_att_crown,
        "max_dist_att": max(dists_att_crown) if dists_att_crown else 0.0,
        "ray_converges_to_crown": ray_converges_to_crown, "n_converging": n_converging,
        "frac_converging": float(n_converging) / float(max(n_leaves_pid, 1)),
        "intrinsic_has_pet": intrinsic_has_pet, "effective_pet_ratios": effective_pet_ratios,
        "effective_has_pet": effective_has_pet, "att_span_px": att_span_px,
        "att_elongation": att_elongation, "stem_axis_vec": stem_axis_vec, "p_bend": p_bend,
    }


def _score_single_leaf_zone(
    idx_item: int, g_item: Mapping[str, Any], ray_info: Mapping[str, Any],
    c_pt: Tuple[float, float], cap_pt: Optional[Tuple[float, float]], r_hub: float,
    flags: Mapping[str, bool],
) -> Tuple[str, float, float]:
    """Compute `(zone_label, zone_conf, dist_norm)` across all 4 anatomical organ-presence cases."""
    traced_att_pts = ray_info["traced_att_pts"]
    p_att = traced_att_pts[idx_item]
    dist_crown_px = ray_info["dists_att_crown"][idx_item]
    max_dist_att = ray_info["max_dist_att"]
    converges_i = ray_info["ray_converges_to_crown"][idx_item]
    pet_ratio = ray_info["effective_pet_ratios"][idx_item]
    has_pet = ray_info["effective_has_pet"][idx_item]
    p_bend = ray_info["p_bend"]

    if flags["is_detached_sterile_rosette"]:
        dist_norm = min(0.24, dist_crown_px / max(3.5 * r_hub, 150.0))
        p_basal = max(0.76, min(0.98, 0.92 - 0.35 * dist_norm + (0.06 if has_pet else 0.0)))
    elif flags["is_rooted_sterile_tuft"] and converges_i:
        dist_norm = min(0.22, dist_crown_px / max(3.5 * r_hub, 150.0))
        p_basal = max(0.78, min(0.98, 0.94 - 0.35 * dist_norm + (0.04 if has_pet else 0.0)))
    elif flags["is_clipped_upper_flowering_stem"]:
        dist_norm = max(0.55, min(1.0, 0.55 + 0.45 * project_along_stem_polyline(p_att, c_pt, cap_pt, p_bend))) if cap_pt is not None else 0.75
        p_basal = max(0.04, min(0.32, 0.35 * (1.0 - dist_norm)))
    elif flags["has_caps"] and cap_pt is not None:
        p_tip_raw = g_item["p_tip"]
        raw_d_tip = float(math.hypot(p_tip_raw[0] - c_pt[0], p_tip_raw[1] - c_pt[1]))
        s_stem_tip = project_along_stem_polyline(p_tip_raw, c_pt, cap_pt, p_bend)
        scape_span = max(math.hypot(cap_pt[0] - c_pt[0], cap_pt[1] - c_pt[1]), 1.5 * r_hub, 100.0)
        dist_norm = max(0.0, min(1.0, 0.78 * s_stem_tip + 0.22 * min(1.0, raw_d_tip / scape_span)))

        is_long_pet_rosette = bool(has_pet and pet_ratio >= 0.30 and s_stem_tip <= 0.22 and (converges_i or raw_d_tip <= 1.95 * r_hub))
        if converges_i and s_stem_tip <= 0.18 and (raw_d_tip <= 1.15 * r_hub or (has_pet and raw_d_tip <= 1.55 * r_hub)):
            dist_norm = min(dist_norm, 0.18)
            p_basal = max(0.68, min(0.98, 0.90 - 0.55 * dist_norm + min(0.10, 0.15 * pet_ratio)))
        elif is_long_pet_rosette:
            dist_norm = min(dist_norm, 0.24)
            p_basal = max(0.64, min(0.95, 0.82 - 0.35 * dist_norm + min(0.10, 0.12 * pet_ratio)))
        elif s_stem_tip <= 0.14 and raw_d_tip <= 1.10 * r_hub:
            p_basal = max(0.62, min(0.95, 0.84 - 0.60 * dist_norm + (0.08 if has_pet else 0.0)))
        elif s_stem_tip <= 0.18 and has_pet and raw_d_tip <= 1.45 * r_hub:
            p_basal = max(0.56, min(0.88, 0.78 - 0.55 * dist_norm + min(0.12, 0.15 * pet_ratio)))
        else:
            p_basal = max(0.04, min(0.46, 0.42 * (1.0 - dist_norm) + (0.03 if (has_pet and s_stem_tip < 0.24) else 0.0)))
    elif flags["has_roots"] and not flags["has_caps"]:
        dist_norm = max(0.0, min(1.0, dist_crown_px / max(max_dist_att, 2.5 * r_hub, 120.0)))
        if converges_i or dist_crown_px <= 1.25 * r_hub or (dist_norm <= 0.42 and has_pet):
            dist_norm = min(dist_norm, 0.30)
            p_basal = max(0.65, min(0.98, 0.88 - 0.45 * dist_norm + (0.08 if has_pet else 0.0)))
        else:
            p_basal = max(0.06, min(0.46, 0.44 * (1.0 - dist_norm)))
    else:
        radial_c4 = max(0.0, min(1.0, dist_crown_px / max(max_dist_att, 2.5 * r_hub, 120.0)))
        sv = ray_info["stem_axis_vec"]
        proj_vals = [(pt_a[0] - c_pt[0]) * sv[0] + (pt_a[1] - c_pt[1]) * sv[1] for pt_a in traced_att_pts]
        p_min = min(proj_vals) if proj_vals else 0.0
        p_ptp = max(max(proj_vals) - p_min, 2.5 * r_hub, 120.0) if proj_vals else 120.0
        s_axis = max(0.0, min(1.0, (proj_vals[idx_item] - p_min) / p_ptp))
        dist_norm = max(radial_c4, s_axis) if ray_info["n_converging"] >= 2 else s_axis
        if flags["is_sessile_stem_fragment"] and not flags["pid_has_explicit_crown"]:
            dist_norm = max(0.55, min(1.0, 0.55 + 0.45 * s_axis))
            p_basal = max(0.04, min(0.42, 0.40 * (1.0 - dist_norm)))
        elif converges_i or dist_crown_px <= 1.15 * r_hub or (dist_norm <= 0.32 and has_pet):
            dist_norm = min(dist_norm, 0.28)
            p_basal = max(0.64, min(0.96, 0.86 - 0.45 * dist_norm + (0.06 if has_pet else 0.0)))
        else:
            p_basal = max(0.06, min(0.46, 0.44 * (1.0 - dist_norm)))

    p_basal = max(0.0, min(1.0, float(p_basal)))
    return ("basal" if p_basal >= 0.50 else "cauline"), round(float(max(p_basal, 1.0 - p_basal)), 4), dist_norm


def classify_leaf_zone_anatomical(
    leaf_records: Sequence[Mapping[str, Any]],
    crowns_by_pid: Optional[Mapping[int, Tuple[float, float]]] = None,
    rosette_centroids_by_pid: Optional[Mapping[int, Optional[Tuple[float, float]]]] = None,
    basal_stats_by_pid: Optional[Mapping[int, Tuple[float, float, float, float]]] = None,
    capitula_detections: Optional[Sequence[Any]] = None,
    root_boxes: Optional[Sequence[Any]] = None,
    is_bent_stem: bool = False,
    sheet_shape: Optional[Tuple[int, int]] = None,
    specimen_mask: Optional[np.ndarray] = None,
    p_crown: Optional[Tuple[float, float]] = None,
    p_corymb: Optional[Tuple[float, float]] = None,
    plant_bbox: Optional[Tuple[float, float, float, float]] = None,
    **kwargs: Any,
) -> Dict[int, Dict[str, Any]]:
    """
    Post-segmentation anatomical leaf zone classifier (`basal` vs. `cauline`)
    operating per `plant_individual_id` via vascular attachment topology and
    organ-pole anchoring across all 4 organ-presence states.

    Never drops or rejects any leaf. Returns `int(leaf_id) -> {leaf_zone,
    leaf_zone_confidence, dist_from_crown_norm, cy_rel_plant_bbox, archetype_score}`.
    """
    _ = sheet_shape
    if not leaf_records:
        return {}

    crowns_map, rosette_map, stats_map, root_pts_all, cap_pts_all, spec_bin, comp_labels = parse_zone_inputs(
        crowns_by_pid, rosette_centroids_by_pid, basal_stats_by_pid,
        capitula_detections, root_boxes, specimen_mask, p_crown, p_corymb, kwargs,
    )

    records_by_pid: Dict[int, List[Tuple[int, Mapping[str, Any]]]] = {}
    for idx, rec in enumerate(leaf_records):
        lid = int(rec["leaf_id"]) if rec.get("leaf_id") is not None else (idx + 1)
        pid = int(rec["plant_individual_id"]) if rec.get("plant_individual_id") is not None else 0
        records_by_pid.setdefault(pid, []).append((lid, rec))

    pid_ctx = build_pid_geometry_context(records_by_pid, comp_labels)
    results: Dict[int, Dict[str, Any]] = {}
    explicit_p_bend: Optional[Tuple[float, float]] = kwargs.get("p_bend")

    for pid, pid_items in records_by_pid.items():
        g_list = pid_ctx["pid_geoms"][pid]
        blade_lengths = pid_ctx["pid_blade_lengths"][pid]
        med_blade_px = pid_ctx["pid_med_blade"][pid]
        plant_span = pid_ctx["pid_span"][pid]

        roots_pid = scope_organ_pts_to_pid(root_pts_all, pid, records_by_pid, pid_ctx, comp_labels, 1.65, 4.5, 380.0)
        caps_pid = scope_organ_pts_to_pid(cap_pts_all, pid, records_by_pid, pid_ctx, comp_labels, 1.85, 5.5, 450.0)
        has_roots, has_caps = len(roots_pid) > 0, len(caps_pid) > 0
        other_centers = [
            (float(np.median([p[0] for p in pid_ctx["pid_anchor_pts"][ok]])), float(np.median([p[1] for p in pid_ctx["pid_anchor_pts"][ok]])))
            for ok in records_by_pid.keys()
            if ok != pid and pid_ctx["pid_anchor_pts"].get(ok)
        ]

        c_pt, cap_pt, r_hub = _resolve_pid_poles_and_hub_radius(
            pid, g_list, roots_pid, caps_pid, other_centers, med_blade_px, crowns_map, rosette_map, stats_map
        )
        pid_has_explicit_crown = bool(pid in crowns_map or pid in rosette_map)
        ray_info = _trace_pid_inward_rays(
            pid_items, g_list, blade_lengths, med_blade_px, plant_span, c_pt, cap_pt, r_hub,
            has_roots, has_caps, pid_has_explicit_crown, spec_bin, is_bent_stem, explicit_p_bend,
        )
        cy_rel_by_leaf = compute_pid_vertical_metrics(
            pid_items, g_list, roots_pid, caps_pid, c_pt, plant_bbox, len(records_by_pid)
        )

        n_leaves_pid = len(pid_items)
        dists_att_crown = ray_info["dists_att_crown"]
        max_dist_att = ray_info["max_dist_att"]
        att_span_px, att_elongation = ray_info["att_span_px"], ray_info["att_elongation"]
        frac_converging = ray_info["frac_converging"]
        intrinsic_has_pet = ray_info["intrinsic_has_pet"]
        n_pet_pid = sum(1 for hp in ray_info["effective_has_pet"] if hp)

        is_clipped_upper_flowering_stem = bool(
            has_caps and not has_roots
            and sum(1 for i in range(n_leaves_pid) if intrinsic_has_pet[i] and dists_att_crown[i] <= 1.85 * r_hub) == 0
        )
        is_sessile_stem_fragment = False
        is_detached_sterile_rosette = False
        if not has_roots and not has_caps:
            is_sessile_stem_fragment = n_pet_pid == 0 and n_leaves_pid >= 2 and att_span_px > 1.15 * med_blade_px
            if not is_sessile_stem_fragment and (
                (n_leaves_pid <= 2 and max_dist_att <= 2.2 * r_hub and (n_pet_pid >= 1 or att_span_px <= 1.15 * med_blade_px))
                or (frac_converging >= 0.85 and max_dist_att <= 1.85 * r_hub and n_pet_pid >= 1)
                or (att_span_px < 1.85 * med_blade_px and (n_pet_pid >= 1 or frac_converging >= 0.65))
                or (att_elongation < 2.4 and max_dist_att <= 1.85 * r_hub and n_pet_pid >= 1)
            ):
                is_detached_sterile_rosette = True

        flags = {
            "has_roots": has_roots,
            "has_caps": has_caps,
            "pid_has_explicit_crown": pid_has_explicit_crown,
            "is_detached_sterile_rosette": is_detached_sterile_rosette,
            "is_rooted_sterile_tuft": bool(has_roots and not has_caps and (frac_converging >= 0.70 or att_span_px < 2.2 * r_hub)),
            "is_clipped_upper_flowering_stem": is_clipped_upper_flowering_stem,
            "is_sessile_stem_fragment": is_sessile_stem_fragment,
        }

        for idx_item, (lid, rec) in enumerate(pid_items):
            zone_label, zone_conf, dist_norm = _score_single_leaf_zone(
                idx_item, g_list[idx_item], ray_info, c_pt, cap_pt, r_hub, flags
            )
            results[int(lid)] = {
                "leaf_zone": zone_label,
                "leaf_zone_confidence": zone_conf,
                "dist_from_crown_norm": round(float(dist_norm), 4),
                "cy_rel_plant_bbox": round(float(cy_rel_by_leaf[idx_item]), 4),
                "archetype_score": round(float(compute_leaf_archetype_score(rec, zone_conf)), 4),
            }

    return results
