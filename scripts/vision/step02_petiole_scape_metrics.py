#!/usr/bin/env python3
"""
Step 02 empirical basal petiole length, petiole width, and straight/bent
flowering scape height approximations (`scripts/vision/step02_petiole_scape_metrics.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from scripts.vision.step02_rosette_crown_pose import (
    _extract_xy_point,
    _parse_basal_leaf_pose,
)
from scripts.vision.step02_rosette_crown_solver import estimate_rosette_crown_point


def approximate_petiole_and_scape_metrics(
    basal_leaves: Any = None,
    p_crown: Optional[Any] = None,
    px_per_cm: Optional[float] = None,
    p_corymb: Optional[Tuple[float, float]] = None,
    p_bend: Optional[Tuple[float, float]] = None,
    upper_cauline_detections: Optional[Sequence[Any]] = None,
    stem_skeleton_mask: Optional[np.ndarray] = None,
    root_boxes: Optional[Sequence[Any]] = None,
    capitula_detections: Optional[Sequence[Any]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    Compute empirical approximations for basal petiole length, petiole width, and flowering plant height:
      1. Per-leaf:
         - When `keypoints_31` has valid `petiole_0..4` / `petiole_tip` (or `petiole_1..6`)
           keypoints (`confidence >= 0.25` via `_parse_basal_leaf_pose`) or a segmented
           petiole skeleton, honors LM3's direct `max(d_pose_polyline, geodesic_petiole_skeleton_px)`
           measurement without adding unsegmented crown-gap continuations, falling back to
           `||p_base - p_crown||_2` only when both direct petiole measurements are absent.
         - `blade_length_px = ||p_apex - p_base||_2`
         - `petiole_to_blade_ratio = petiole_length_px / blade_length_px`
         - `petiole_length_mm`, `blade_length_mm`, and `petiole_width_mm` (from LM3 Module 13
           `petiole_width` or `width_1`/`width_2` when available, scaled via `px_per_mm = px_per_cm / 10.0`)
      2. Specimen-level outer-rosette aggregation (75th percentile across basal leaves per specimen):
         - `petiole_length_mm`, `blade_length_mm`, `petiole_to_blade_ratio`, `petiole_width_mm`
      3. Per-voucher flowering scape height:
         - Traces the piecewise path from `p_crown` through any upper stem-bend apex (`p_bend`
           from upper cauline detections / specimen stem skeleton when bent) to the flowering
           capitula centroid (`p_corymb`), exporting `plant_height_cm`, `scape_to_basal_leaf_ratio`,
           and `is_bent_stem`.
    """
    # Resolve scale conversion strictly when ruler-calibrated (`is_calibrated_mm` is not False)
    is_cal_arg = kwargs.get("is_calibrated_mm", True)
    eff_px_per_cm: Optional[float] = None
    if is_cal_arg is not False:
        if px_per_cm is not None and not np.isnan(px_per_cm) and float(px_per_cm) > 0.0:
            eff_px_per_cm = float(px_per_cm)
        else:
            ppm = kwargs.get("pixels_per_mm")
            if ppm is not None and not np.isnan(ppm) and float(ppm) > 0.0:
                eff_px_per_cm = float(ppm) * 10.0
    px_per_mm = (eff_px_per_cm / 10.0) if (eff_px_per_cm is not None and eff_px_per_cm > 0.0) else None

    # Parse basal leaf poses
    if basal_leaves is None and kwargs.get("leaf_records") is not None:
        basal_leaves = kwargs["leaf_records"]
    poses: List[Dict[str, Any]] = []
    if isinstance(basal_leaves, (dict, np.ndarray)):
        parsed = _parse_basal_leaf_pose(basal_leaves)
        if parsed is not None:
            poses.append(parsed)
    elif basal_leaves is not None:
        for item in basal_leaves:
            parsed = _parse_basal_leaf_pose(item)
            if parsed is not None:
                poses.append(parsed)
    elif kwargs.get("p_base") is not None and kwargs.get("p_apex") is not None:
        parsed = _parse_basal_leaf_pose({
            "p_base": kwargs["p_base"],
            "p_apex": kwargs["p_apex"],
            "geodesic_petiole_skeleton_px": kwargs.get("geodesic_petiole_skeleton_px", 0.0),
            "blade_length_px": kwargs.get("blade_length_px"),
            "keypoints_31": kwargs.get("keypoints_31"),
            "petiole_width_mm": kwargs.get("petiole_width_mm"),
            "petiole_width_px": kwargs.get("petiole_width_px"),
        })
        if parsed is not None:
            poses.append(parsed)

    # Resolve p_crown per plant_individual_id and primary voucher crown
    crowns_by_pid: Dict[int, Tuple[float, float]] = {}
    if isinstance(p_crown, dict):
        for k, v in p_crown.items():
            if v is not None and len(v) == 2:
                crowns_by_pid[int(k)] = (float(v[0]), float(v[1]))
    elif p_crown is not None and isinstance(p_crown, (tuple, list, np.ndarray)) and len(p_crown) == 2:
        pt_c = (float(p_crown[0]), float(p_crown[1]))
        for p in poses:
            crowns_by_pid[int(p["plant_individual_id"])] = pt_c
        crowns_by_pid.setdefault(0, pt_c)

    if not crowns_by_pid:
        est_crowns = estimate_rosette_crown_point(
            poses,
            root_boxes=root_boxes or kwargs.get("root_detections") or kwargs.get("roots"),
            return_by_plant=True,
        )
        if isinstance(est_crowns, dict):
            crowns_by_pid.update(est_crowns)
        else:
            crowns_by_pid[0] = est_crowns

    primary_crown = crowns_by_pid.get(0) or next(iter(crowns_by_pid.values()), (0.0, 0.0))

    # Compute per-leaf metrics
    leaf_metrics: List[Dict[str, Any]] = []
    for idx_l, p in enumerate(poses):
        pid = int(p["plant_individual_id"])
        c_pt = crowns_by_pid.get(pid, primary_crown)
        pb = p["p_base"]
        pa = p["p_apex"]
        geodesic_px = float(p["geodesic_petiole_skeleton_px"])
        crown_dist_px = float(math.hypot(pb[0] - c_pt[0], pb[1] - c_pt[1]))

        blade_px = float(p["blade_length_px"])
        if blade_px <= 0.0:
            blade_px = float(math.hypot(pa[0] - pb[0], pa[1] - pb[1]))

        is_caul_leaf = str(p.get("leaf_zone", "basal") or "basal").strip().lower() == "cauline"
        base_poly_px = float(p.get("pose_petiole_polyline_px", 0.0) or 0.0)
        p_tip = p.get("p_petiole_tip")
        d_pose_polyline = base_poly_px

        crown_fallback_px = (
            crown_dist_px
            if (not is_caul_leaf and blade_px > 0.0 and crown_dist_px <= 1.75 * blade_px)
            else 0.0
        )
        direct_petiole_px = max(float(d_pose_polyline), float(geodesic_px))
        if direct_petiole_px > 0.0:
            petiole_px = direct_petiole_px
        elif is_caul_leaf:
            petiole_px = 0.0
        elif crown_fallback_px > 0.0:
            petiole_px = crown_fallback_px
        else:
            petiole_px = np.nan
        pet_to_blade = (petiole_px / blade_px) if (np.isfinite(petiole_px) and blade_px > 1e-6) else np.nan
        petiole_mm = (
            (petiole_px / px_per_mm)
            if (np.isfinite(petiole_px) and px_per_mm is not None and px_per_mm > 0.0)
            else np.nan
        )
        blade_mm = (blade_px / px_per_mm) if (px_per_mm is not None and px_per_mm > 0.0) else np.nan

        pw_px = p.get("petiole_width_px")
        pw_mm = p.get("petiole_width_mm")
        if pw_mm is None and pw_px is not None and px_per_mm is not None and px_per_mm > 0.0:
            pw_mm = float(pw_px) / px_per_mm
        elif pw_px is None and pw_mm is not None and px_per_mm is not None and px_per_mm > 0.0:
            pw_px = float(pw_mm) * px_per_mm

        leaf_metrics.append({
            "leaf_id": p["leaf_id"] if p["leaf_id"] is not None else (idx_l + 1),
            "plant_individual_id": pid,
            "leaf_zone": "cauline" if is_caul_leaf else "basal",
            "p_base": pb,
            "p_apex": pa,
            "p_petiole_tip": p_tip,
            "p_crown": c_pt,
            "pose_petiole_polyline_px": round(float(d_pose_polyline), 4),
            "geodesic_petiole_skeleton_px": round(max(geodesic_px, float(d_pose_polyline)), 4),
            "crown_distance_px": round(crown_dist_px, 4),
            "petiole_length_px": round(float(petiole_px), 4) if np.isfinite(petiole_px) else np.nan,
            "blade_length_px": round(blade_px, 4),
            "petiole_to_blade_ratio": round(float(pet_to_blade), 4) if np.isfinite(pet_to_blade) else np.nan,
            "petiole_length_mm": round(float(petiole_mm), 3) if np.isfinite(petiole_mm) else np.nan,
            "blade_length_mm": round(float(blade_mm), 3) if np.isfinite(blade_mm) else np.nan,
            "petiole_width_px": round(float(pw_px), 4) if (pw_px is not None and np.isfinite(float(pw_px))) else np.nan,
            "petiole_width_mm": round(float(pw_mm), 3) if (pw_mm is not None and np.isfinite(float(pw_mm))) else np.nan,
        })

    def _agg_p75(records: List[Dict[str, Any]], key: str, digits: int = 4) -> float:
        vals = [float(r[key]) for r in records if r.get(key) is not None and np.isfinite(float(r[key]))]
        if not vals:
            return np.nan
        return round(float(np.percentile(vals, 75)), digits)

    # Aggregate 75th-percentile outer-rosette metrics per specimen (plant_individual_id) and overall
    by_plant_summary: Dict[int, Dict[str, Any]] = {}
    pids = sorted({int(r["plant_individual_id"]) for r in leaf_metrics})
    for pid in pids:
        sub_all = [r for r in leaf_metrics if int(r["plant_individual_id"]) == pid]
        sub_bas = [r for r in sub_all if r.get("leaf_zone") != "cauline"]
        sub_r = sub_bas if sub_bas else sub_all
        pid_cr = crowns_by_pid.get(pid, primary_crown)
        pid_pet_p75 = _agg_p75(sub_r, "petiole_length_px", 4)
        pid_bld_p75 = _agg_p75(sub_r, "blade_length_px", 4)
        pid_scape_px = np.nan
        pid_height_cm = np.nan
        pid_scape_ratio = np.nan
        pid_bent = False
        caps_by_pid = kwargs.get("capitula_by_pid")
        if isinstance(caps_by_pid, dict):
            pid_caps = caps_by_pid.get(pid) or []
            pid_cap_pts = [pt for c_it in pid_caps if (pt := _extract_xy_point(c_it, ref_pt=pid_cr)) is not None]
            if pid_cap_pts:
                pid_cor = (
                    float(np.median([p[0] for p in pid_cap_pts])),
                    float(np.median([p[1] for p in pid_cap_pts])),
                )
                pid_scape_px = float(math.hypot(pid_cor[0] - pid_cr[0], pid_cor[1] - pid_cr[1]))
                if eff_px_per_cm is not None and eff_px_per_cm > 0.0:
                    pid_height_cm = round(float(pid_scape_px / eff_px_per_cm), 3)
                pid_ref = (
                    (pid_pet_p75 + pid_bld_p75)
                    if (np.isfinite(pid_pet_p75) and np.isfinite(pid_bld_p75) and (pid_pet_p75 + pid_bld_p75) > 1e-6)
                    else (pid_bld_p75 if (np.isfinite(pid_bld_p75) and pid_bld_p75 > 1e-6) else np.nan)
                )
                if np.isfinite(pid_ref) and pid_ref > 1e-6:
                    pid_scape_ratio = round(float(pid_scape_px / pid_ref), 4)
        by_plant_summary[pid] = {
            "p_crown": pid_cr,
            "petiole_length_px": pid_pet_p75,
            "blade_length_px": pid_bld_p75,
            "petiole_length_mm": _agg_p75(sub_r, "petiole_length_mm", 3),
            "blade_length_mm": _agg_p75(sub_r, "blade_length_mm", 3),
            "petiole_to_blade_ratio": _agg_p75(sub_r, "petiole_to_blade_ratio", 4),
            "petiole_width_px": _agg_p75(sub_r, "petiole_width_px", 4),
            "petiole_width_mm": _agg_p75(sub_r, "petiole_width_mm", 3),
            "scape_length_px": round(float(pid_scape_px), 2) if np.isfinite(pid_scape_px) else np.nan,
            "plant_height_cm": pid_height_cm,
            "scape_to_basal_leaf_ratio": pid_scape_ratio,
            "is_bent_stem": pid_bent,
        }

    bas_all = [r for r in leaf_metrics if r.get("leaf_zone") != "cauline"]
    agg_pool = bas_all if bas_all else leaf_metrics
    petiole_px_p75 = _agg_p75(agg_pool, "petiole_length_px", 4)
    blade_px_p75 = _agg_p75(agg_pool, "blade_length_px", 4)
    petiole_mm_p75 = _agg_p75(agg_pool, "petiole_length_mm", 3)
    blade_mm_p75 = _agg_p75(agg_pool, "blade_length_mm", 3)
    pet_to_blade_p75 = _agg_p75(agg_pool, "petiole_to_blade_ratio", 4)
    petiole_width_px_p75 = _agg_p75(agg_pool, "petiole_width_px", 4)
    petiole_width_mm_p75 = _agg_p75(agg_pool, "petiole_width_mm", 3)

    # Resolve flowering capitula centroid p_corymb
    resolved_corymb: Optional[Tuple[float, float]] = None
    if p_corymb is not None and len(p_corymb) == 2:
        resolved_corymb = (float(p_corymb[0]), float(p_corymb[1]))
    else:
        raw_caps = capitula_detections or kwargs.get("capitula_boxes") or kwargs.get("capitula_pts")
        if raw_caps:
            cap_pts = []
            for c_item in raw_caps:
                pt = _extract_xy_point(c_item, ref_pt=primary_crown)
                if pt is not None:
                    cap_pts.append(pt)
            if cap_pts:
                resolved_corymb = (
                    float(np.median([p[0] for p in cap_pts])),
                    float(np.median([p[1] for p in cap_pts])),
                )

    # Trace piecewise scape path from p_crown through any upper stem-bend apex (p_bend) to p_corymb
    is_bent_stem = False
    resolved_bend: Optional[Tuple[float, float]] = None
    scape_length_px = np.nan
    plant_height_cm = np.nan
    scape_to_basal_leaf_ratio = np.nan

    if resolved_corymb is not None:
        straight_px = float(math.hypot(resolved_corymb[0] - primary_crown[0], resolved_corymb[1] - primary_crown[1]))
        explicit_bend = p_bend or kwargs.get("bend_apex") or kwargs.get("stem_bend_pt")
        cand_bends: List[Tuple[Tuple[float, float], bool]] = []
        if explicit_bend is not None:
            pt_b = _extract_xy_point(explicit_bend, ref_pt=primary_crown)
            if pt_b is not None:
                cand_bends.append((pt_b, True))

        raw_cauline = (
            upper_cauline_detections
            or kwargs.get("cauline_detections")
            or kwargs.get("cauline_boxes")
            or kwargs.get("upper_cauline_pts")
            or kwargs.get("stem_points")
        )
        if raw_cauline:
            for c_item in raw_cauline:
                pt_b = _extract_xy_point(c_item, ref_pt=primary_crown)
                if pt_b is not None:
                    cand_bends.append((pt_b, False))

        skel_mask = stem_skeleton_mask if stem_skeleton_mask is not None else kwargs.get("specimen_mask")
        if isinstance(skel_mask, np.ndarray) and np.count_nonzero(skel_mask) > 0:
            ys_s, xs_s = np.where(skel_mask > 0)
            min_y_idx = int(np.argmin(ys_s))
            top_y = float(ys_s[min_y_idx])
            top_x = float(np.median(xs_s[ys_s <= top_y + 3]))
            cand_bends.append(((top_x, top_y), False))
            if straight_px > 1e-4:
                chord_u = np.array(
                    [resolved_corymb[0] - primary_crown[0], resolved_corymb[1] - primary_crown[1]],
                    dtype=np.float64,
                ) / straight_px
                pts_rel = np.column_stack([xs_s - primary_crown[0], ys_s - primary_crown[1]])
                perp_d = np.abs(pts_rel[:, 0] * (-chord_u[1]) + pts_rel[:, 1] * chord_u[0])
                max_p_idx = int(np.argmax(perp_d))
                cand_bends.append(((float(xs_s[max_p_idx]), float(ys_s[max_p_idx])), False))

        best_piecewise_px = straight_px
        for bend_pt, is_explicit in cand_bends:
            d1 = float(math.hypot(bend_pt[0] - primary_crown[0], bend_pt[1] - primary_crown[1]))
            d2 = float(math.hypot(resolved_corymb[0] - bend_pt[0], resolved_corymb[1] - bend_pt[1]))
            pw_len = d1 + d2
            if straight_px > 1e-4:
                ux = (resolved_corymb[0] - primary_crown[0]) / straight_px
                uy = (resolved_corymb[1] - primary_crown[1]) / straight_px
                perp_dist = abs((bend_pt[0] - primary_crown[0]) * (-uy) + (bend_pt[1] - primary_crown[1]) * ux)
            else:
                perp_dist = d1

            is_upper_fold = bend_pt[1] < (min(primary_crown[1], resolved_corymb[1]) - 8.0)
            is_lateral_bend = (
                bend_pt[1] < primary_crown[1]
                and perp_dist >= max(15.0, 0.12 * straight_px)
                and pw_len >= 1.04 * straight_px
            )
            if (is_explicit and (pw_len > straight_px + 1.0 or perp_dist > 2.0)) or is_upper_fold or is_lateral_bend:
                if pw_len >= best_piecewise_px:
                    best_piecewise_px = pw_len
                    resolved_bend = bend_pt
                    is_bent_stem = True

        scape_length_px = best_piecewise_px if is_bent_stem else straight_px
        if eff_px_per_cm is not None and eff_px_per_cm > 0.0:
            plant_height_cm = round(float(scape_length_px / eff_px_per_cm), 3)

        if np.isfinite(petiole_px_p75) and np.isfinite(blade_px_p75) and (petiole_px_p75 + blade_px_p75) > 1e-6:
            basal_ref_px = float(petiole_px_p75 + blade_px_p75)
        elif np.isfinite(blade_px_p75) and blade_px_p75 > 1e-6:
            basal_ref_px = float(blade_px_p75)
        else:
            basal_ref_px = np.nan

        if np.isfinite(basal_ref_px) and basal_ref_px > 1e-6:
            scape_to_basal_leaf_ratio = round(float(scape_length_px / basal_ref_px), 4)

    return {
        "p_crown": primary_crown,
        "crowns_by_plant": crowns_by_pid,
        "p_corymb": resolved_corymb,
        "p_bend": resolved_bend,
        "is_bent_stem": bool(is_bent_stem),
        "scape_length_px": round(float(scape_length_px), 4) if np.isfinite(scape_length_px) else np.nan,
        "plant_height_cm": plant_height_cm,
        "scape_to_basal_leaf_ratio": scape_to_basal_leaf_ratio,
        "petiole_length_px": petiole_px_p75,
        "blade_length_px": blade_px_p75,
        "petiole_length_mm": petiole_mm_p75,
        "blade_length_mm": blade_mm_p75,
        "petiole_to_blade_ratio": pet_to_blade_p75,
        "petiole_width_px": petiole_width_px_p75,
        "petiole_width_mm": petiole_width_mm_p75,
        "leaf_metrics": leaf_metrics,
        "by_plant": by_plant_summary,
    }
