#!/usr/bin/env python3
"""
Step 02 rosette leaf point extraction and LeafMachine3 31-keypoint pose parsing
(`scripts/vision/step02_rosette_crown_pose.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np


def _call_mask_fallback(mask: np.ndarray) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Invoke the bare-mask contour endpoint fallback lazily without circular imports."""
    from scripts.vision.organ_crown_metrics import _mask_basal_apical_fallback

    return _mask_basal_apical_fallback(mask)


def _safe_pos_float(val: Any) -> Optional[float]:
    if val is None:
        return None
    try:
        fval = float(val)
    except (TypeError, ValueError):
        return None
    return fval if (np.isfinite(fval) and fval > 0.0) else None


def _extract_xy_point(item: Any, ref_pt: Optional[Tuple[float, float]] = None) -> Optional[Tuple[float, float]]:
    """Convert a 2-tuple point, 4-tuple bbox, or dict into a sheet-space (x, y) coordinate."""
    if item is None:
        return None
    if isinstance(item, dict):
        if "x" in item and "y" in item:
            return float(item["x"]), float(item["y"])
        if "centroid" in item and item["centroid"] is not None:
            c = item["centroid"]
            return float(c[0]), float(c[1])
        if "p_bend" in item and item["p_bend"] is not None:
            b = item["p_bend"]
            return float(b[0]), float(b[1])
        if "bbox" in item and item["bbox"] is not None:
            return _extract_xy_point(item["bbox"], ref_pt=ref_pt)
        return None
    if hasattr(item, "bbox") and not isinstance(item, (tuple, list, np.ndarray)):
        return _extract_xy_point(getattr(item, "bbox"), ref_pt=ref_pt)
    if isinstance(item, (tuple, list, np.ndarray)):
        seq = list(item)
        if len(seq) >= 3 and isinstance(seq[0], (tuple, list, np.ndarray)):
            cls_id = int(seq[-1]) if isinstance(seq[-1], (int, float, np.integer)) else None
            if cls_id is not None and cls_id not in (0, 1, 3, 4, 5, 6, 7, 8, 9, 10):
                return None
            return _extract_xy_point(seq[0], ref_pt=ref_pt)
        if len(seq) == 2 and all(isinstance(v, (int, float, np.number)) for v in seq):
            return float(seq[0]), float(seq[1])
        if len(seq) == 4 and all(isinstance(v, (int, float, np.number)) for v in seq):
            v0, v1, v2, v3 = float(seq[0]), float(seq[1]), float(seq[2]), float(seq[3])
            # Canonical repository bbox order is strictly (ymin, xmin, ymax, xmax) -> (cx, cy) = (0.5*(v1+v3), 0.5*(v0+v2))
            return (0.5 * (v1 + v3), 0.5 * (v0 + v2))
    return None


def _coerce_keypoints_dict(keypoints_31: Any) -> Optional[Dict[str, Any]]:
    """Normalize a 31-keypoint ndarray or mapping into a named keypoint dictionary."""
    if (
        isinstance(keypoints_31, np.ndarray)
        and keypoints_31.ndim == 2
        and keypoints_31.shape[0] >= 8
        and keypoints_31.shape[1] >= 2
    ):
        def _row_tup(r: np.ndarray) -> Tuple[float, float, float]:
            return (float(r[0]), float(r[1]), float(r[2]) if len(r) >= 3 else 1.0)

        kp_dict: Dict[str, Any] = {
            "lamina_tip": _row_tup(keypoints_31[0]),
            "lamina_base": _row_tup(keypoints_31[1]),
        }
        for pi in range(min(6, keypoints_31.shape[0] - 2)):
            kp_dict[f"petiole_{pi + 1}"] = _row_tup(keypoints_31[2 + pi])
        for mi in range(min(15, max(0, keypoints_31.shape[0] - 8))):
            kp_dict[f"midvein_{mi + 1}"] = _row_tup(keypoints_31[8 + mi])
        return kp_dict
    return keypoints_31 if (isinstance(keypoints_31, dict) and len(keypoints_31) > 0) else None


def _apply_keypoints_31_pose(
    item: Any, keypoints_31: Dict[str, Any], kp_in_sheet_space: bool, ox: float, oy: float, state: Dict[str, Any]
) -> None:
    """Refine leaf base/apex, petiole polyline, inward ray direction, and width in `state` from LM3 keypoints."""
    if kp_in_sheet_space:
        kp_ox, kp_oy = 0.0, 0.0
    elif isinstance(item, dict) and item.get("local_mask_origin") is not None and len(item["local_mask_origin"]) == 2:
        kp_ox, kp_oy = float(item["local_mask_origin"][0]), float(item["local_mask_origin"][1])
    else:
        kp_ox, kp_oy = ox, oy

    p_base = state["p_base"]
    if kp_ox == 0.0 and kp_oy == 0.0 and not kp_in_sheet_space and isinstance(item, dict) and p_base is not None:
        raw_lb = keypoints_31.get("lamina_base")
        if isinstance(raw_lb, (tuple, list, np.ndarray)) and len(raw_lb) >= 2:
            rx, ry = float(raw_lb[0]), float(raw_lb[1])
        elif isinstance(raw_lb, dict):
            rx, ry = float(raw_lb.get("x", 0.0)), float(raw_lb.get("y", 0.0))
        else:
            rx, ry = 0.0, 0.0
        if item.get("p_base_sheet") is not None and (rx != 0.0 or ry != 0.0) and math.hypot(p_base[0] - rx, p_base[1] - ry) > 64.0:
            kp_ox, kp_oy = p_base[0] - rx, p_base[1] - ry

    def _extract_kp(name: str, min_conf: float = 0.25) -> Optional[Tuple[float, float, float]]:
        val = keypoints_31.get(name)
        if val is None:
            return None
        if isinstance(val, dict):
            x, y, c = float(val.get("x", 0.0)), float(val.get("y", 0.0)), float(val.get("conf", val.get("confidence", 1.0)))
        elif isinstance(val, (tuple, list, np.ndarray)) and len(val) >= 2:
            x, y, c = float(val[0]), float(val[1]), float(val[2]) if len(val) >= 3 else 1.0
        else:
            return None
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(c)) or c < min_conf or (x == 0.0 and y == 0.0):
            return None
        return (x + kp_ox, y + kp_oy, c)

    lb_kp, lt_kp = _extract_kp("lamina_base", 0.25), _extract_kp("lamina_tip", 0.25)
    mv_range = range(15) if "midvein_0" in keypoints_31 else range(1, 16)
    mv_kps = [kp for kp in (_extract_kp(f"midvein_{i}", 0.25) for i in mv_range) if kp is not None]

    pet_spec = (
        [(i + 1, f"petiole_{i}") for i in range(5)] + [(6, "petiole_tip")]
        if "petiole_0" in keypoints_31
        else [(i, f"petiole_{i}") for i in range(1, 7)] + [(7, "petiole_tip")]
    )
    pet_kps: List[Tuple[int, Tuple[float, float, float]]] = []
    for idx_p, key_p in pet_spec:
        pkp = _extract_kp(key_p, 0.25)
        if pkp is not None and (not pet_kps or math.hypot(pkp[0] - pet_kps[-1][1][0], pkp[1] - pet_kps[-1][1][1]) > 1e-4):
            pet_kps.append((idx_p, pkp))

    if lb_kp is not None:
        state["p_base"] = (lb_kp[0], lb_kp[1])
    elif pet_kps and state["p_base"] is None:
        state["p_base"] = (pet_kps[0][1][0], pet_kps[0][1][1])

    if lt_kp is not None:
        state["p_apex"] = (lt_kp[0], lt_kp[1])
    elif mv_kps and state["p_apex"] is None and state["p_base"] is not None:
        pb_cur = state["p_base"]
        far_mv = max(mv_kps, key=lambda kp: math.hypot(kp[0] - pb_cur[0], kp[1] - pb_cur[1]))
        state["p_apex"] = (far_mv[0], far_mv[1])

    if pet_kps:
        state["has_lm3_petiole_pose"] = True
        p_petiole_tip = (pet_kps[-1][1][0], pet_kps[-1][1][1])
        state["p_petiole_tip"] = p_petiole_tip

        poly_pts: List[Tuple[float, float]] = [(lb_kp[0], lb_kp[1])] if lb_kp is not None else []
        for _, pkp in pet_kps:
            pt_xy = (pkp[0], pkp[1])
            if not poly_pts or math.hypot(pt_xy[0] - poly_pts[-1][0], pt_xy[1] - poly_pts[-1][1]) > 1e-4:
                poly_pts.append(pt_xy)
        if len(poly_pts) >= 2:
            state["pose_petiole_polyline_px"] = float(
                sum(math.hypot(poly_pts[k + 1][0] - poly_pts[k][0], poly_pts[k + 1][1] - poly_pts[k][1]) for k in range(len(poly_pts) - 1))
            )

        spine_pts: List[Tuple[float, float]] = []
        if lt_kp is not None:
            spine_pts.append((lt_kp[0], lt_kp[1]))
        spine_pts.extend((mkp[0], mkp[1]) for mkp in mv_kps)
        if lb_kp is not None:
            spine_pts.append((lb_kp[0], lb_kp[1]))
        spine_pts.extend((pkp[0], pkp[1]) for _, pkp in pet_kps)

        if len(spine_pts) >= 2 and state["p_base"] is not None:
            pet_vec_x, pet_vec_y = p_petiole_tip[0] - state["p_base"][0], p_petiole_tip[1] - state["p_base"][1]
            pet_vec_len = math.hypot(pet_vec_x, pet_vec_y)
            lp = cv2.fitLine(np.asarray(spine_pts, dtype=np.float32), cv2.DIST_L2, 0, 0.01, 0.01)
            vx_s, vy_s = float(lp[0][0]), float(lp[1][0])
            if ((spine_pts[-1][0] - spine_pts[0][0]) * vx_s + (spine_pts[-1][1] - spine_pts[0][1]) * vy_s) < 0.0:
                vx_s, vy_s = -vx_s, -vy_s
            if pet_vec_len > 2.0:
                ux_pet, uy_pet = pet_vec_x / pet_vec_len, pet_vec_y / pet_vec_len
                bx_dir, by_dir = 0.65 * ux_pet + 0.35 * vx_s, 0.65 * uy_pet + 0.35 * vy_s
                bnorm = math.hypot(bx_dir, by_dir)
                state["ray_dir"] = (bx_dir / bnorm, by_dir / bnorm) if bnorm > 1e-6 else (ux_pet, uy_pet)
            else:
                state["ray_dir"] = (vx_s, vy_s)

    if state["petiole_width_mm"] is None:
        state["petiole_width_mm"] = _safe_pos_float(keypoints_31.get("petiole_width_mm"))
    if state["petiole_width_px"] is None:
        state["petiole_width_px"] = _safe_pos_float(keypoints_31.get("petiole_width_px", keypoints_31.get("petiole_width")))
    if state["petiole_width_px"] is None and state["petiole_width_mm"] is None:
        w1_kp = _extract_kp("width_left", 0.25) or _extract_kp("width_1", 0.25)
        w2_kp = _extract_kp("width_right", 0.25) or _extract_kp("width_2", 0.25)
        if w1_kp is not None and w2_kp is not None:
            w_dist = math.hypot(w1_kp[0] - w2_kp[0], w1_kp[1] - w2_kp[1])
            if w_dist > 0.0:
                state["petiole_width_px"] = float(w_dist)


def _parse_basal_leaf_pose(item: Any) -> Optional[Dict[str, Any]]:
    """
    Extract standardized `(p_base, p_apex, plant_individual_id, geodesic_petiole_skeleton_px,
    blade_length_px)` and LeafMachine3 31-keypoint pose attributes (`p_petiole_tip`,
    `pose_petiole_polyline_px`, `ray_dir`, `petiole_width_px`, `petiole_width_mm`) from a leaf record.
    """
    if item is None:
        return None

    pid, ox, oy = 0, 0.0, 0.0
    geodesic_px = 0.0
    blade_px: Optional[float] = None
    leaf_id: Optional[int] = None
    leaf_zone = "basal"
    keypoints_31: Any = None
    kp_in_sheet_space = False
    state: Dict[str, Any] = {
        "p_base": None, "p_apex": None, "p_petiole_tip": None, "ray_dir": None,
        "has_lm3_petiole_pose": False, "pose_petiole_polyline_px": 0.0,
        "petiole_width_px": None, "petiole_width_mm": None,
    }

    if isinstance(item, dict):
        raw_pid = item.get("plant_individual_id")
        pid = int(raw_pid) if raw_pid is not None else 0
        leaf_id = item.get("leaf_id")
        leaf_zone = str(item.get("leaf_zone", "basal") or "basal").strip().lower()
        origin = item.get("mask_origin") or item.get("origin")
        if origin is not None and len(origin) == 2 and not item.get("coords_in_sheet_space", False):
            ox, oy = float(origin[0]), float(origin[1])

        keypoints_31 = item.get("keypoints_31")
        kp_in_sheet_space = bool(item.get("keypoints_in_sheet_space", item.get("coords_in_sheet_space", False)))

        if item.get("p_base_sheet") is not None:
            pb = item["p_base_sheet"]
            state["p_base"] = (float(pb[0]), float(pb[1]))
        elif item.get("p_base") is not None:
            pb = item["p_base"]
            state["p_base"] = (float(pb[0]) + ox, float(pb[1]) + oy)

        if item.get("p_apex_sheet") is not None:
            pa = item["p_apex_sheet"]
            state["p_apex"] = (float(pa[0]), float(pa[1]))
        elif item.get("p_apex") is not None:
            pa = item["p_apex"]
            state["p_apex"] = (float(pa[0]) + ox, float(pa[1]) + oy)

        if item.get("p_petiole_tip") is not None:
            ptip = item["p_petiole_tip"]
            state["p_petiole_tip"] = (
                float(ptip[0]) + (0.0 if kp_in_sheet_space else ox),
                float(ptip[1]) + (0.0 if kp_in_sheet_space else oy),
            )
        if item.get("ray_dir") is not None and len(item["ray_dir"]) == 2:
            rd = item["ray_dir"]
            state["ray_dir"] = (float(rd[0]), float(rd[1]))
        if item.get("has_lm3_petiole_pose") is not None:
            state["has_lm3_petiole_pose"] = bool(item["has_lm3_petiole_pose"])
        if item.get("pose_petiole_polyline_px") is not None:
            state["pose_petiole_polyline_px"] = max(0.0, float(item["pose_petiole_polyline_px"] or 0.0))

        state["petiole_width_mm"] = _safe_pos_float(item.get("petiole_width_mm"))
        state["petiole_width_px"] = _safe_pos_float(item.get("petiole_width_px", item.get("petiole_width")))

        if (state["p_base"] is None or state["p_apex"] is None) and isinstance(item.get("mask"), np.ndarray):
            pa_m, pb_m = _call_mask_fallback(item["mask"])
            if state["p_base"] is None:
                state["p_base"] = (float(pb_m[0]) + ox, float(pb_m[1]) + oy)
            if state["p_apex"] is None:
                state["p_apex"] = (float(pa_m[0]) + ox, float(pa_m[1]) + oy)

        trim_info = item.get("trim_meta") or item.get("petiole_trim_info") or {}
        geodesic_px = float(
            item.get("geodesic_petiole_skeleton_px", item.get("trim_distance_px", trim_info.get("trim_distance_px", 0.0))) or 0.0
        )
        if item.get("blade_length_px") is not None and float(item["blade_length_px"]) > 0:
            blade_px = float(item["blade_length_px"])

    elif isinstance(item, np.ndarray) and item.ndim == 2:
        pa_m, pb_m = _call_mask_fallback(item)
        state["p_base"], state["p_apex"] = (float(pb_m[0]), float(pb_m[1])), (float(pa_m[0]), float(pa_m[1]))

    elif isinstance(item, (tuple, list)) and len(item) >= 2:
        first, second = item[0], item[1]
        if (
            isinstance(first, (tuple, list, np.ndarray))
            and isinstance(second, (tuple, list, np.ndarray))
            and len(first) == 2
            and len(second) == 2
        ):
            state["p_base"] = (float(first[0]), float(first[1]))
            state["p_apex"] = (float(second[0]), float(second[1]))
            if len(item) >= 3 and isinstance(item[2], (int, float, np.number)):
                geodesic_px = float(item[2])
            if len(item) >= 4 and isinstance(item[3], (int, np.integer)):
                pid = int(item[3])

    elif hasattr(item, "p_base") or hasattr(item, "mask") or hasattr(item, "keypoints_31"):
        raw_pid = getattr(item, "plant_individual_id", None)
        pid = int(raw_pid) if raw_pid is not None else 0
        leaf_id = getattr(item, "leaf_id", None)
        leaf_zone = str(getattr(item, "leaf_zone", "basal") or "basal").strip().lower()
        origin = getattr(item, "mask_origin", (0, 0)) or (0, 0)
        ox, oy = float(origin[0]), float(origin[1])
        keypoints_31 = getattr(item, "keypoints_31", None)
        kp_in_sheet_space = bool(getattr(item, "keypoints_in_sheet_space", False))
        state["petiole_width_mm"] = _safe_pos_float(getattr(item, "petiole_width_mm", None))
        state["petiole_width_px"] = _safe_pos_float(getattr(item, "petiole_width_px", getattr(item, "petiole_width", None)))
        pb, pa = getattr(item, "p_base", None), getattr(item, "p_apex", None)
        if (pb is None or pa is None) and getattr(item, "mask", None) is not None:
            pa, pb = _call_mask_fallback(getattr(item, "mask"))
        if pb is not None and pa is not None:
            state["p_base"] = (float(pb[0]) + ox, float(pb[1]) + oy)
            state["p_apex"] = (float(pa[0]) + ox, float(pa[1]) + oy)
        trim_info = getattr(item, "trim_meta", None) or {}
        geodesic_px = float(getattr(item, "geodesic_petiole_skeleton_px", trim_info.get("trim_distance_px", 0.0)) or 0.0)

    kp_dict = _coerce_keypoints_dict(keypoints_31)
    if kp_dict is not None:
        _apply_keypoints_31_pose(item, kp_dict, kp_in_sheet_space, ox, oy, state)

    p_base, p_apex = state["p_base"], state["p_apex"]
    if p_base is None or p_apex is None:
        return None

    chord_len = math.hypot(p_apex[0] - p_base[0], p_apex[1] - p_base[1])
    if blade_px is None or blade_px <= 0.0:
        blade_px = chord_len
    ray_dir = state["ray_dir"]
    if ray_dir is None and chord_len > 1e-6:
        ray_dir = ((p_base[0] - p_apex[0]) / chord_len, (p_base[1] - p_apex[1]) / chord_len)

    pose_poly = max(0.0, float(state["pose_petiole_polyline_px"]))
    return {
        "leaf_id": leaf_id,
        "plant_individual_id": pid,
        "leaf_zone": leaf_zone if leaf_zone in ("basal", "cauline") else "basal",
        "p_base": p_base,
        "p_apex": p_apex,
        "p_petiole_tip": state["p_petiole_tip"],
        "ray_dir": ray_dir,
        "has_lm3_petiole_pose": bool(state["has_lm3_petiole_pose"]),
        "pose_petiole_polyline_px": pose_poly,
        "geodesic_petiole_skeleton_px": max(0.0, max(float(geodesic_px), pose_poly)),
        "blade_length_px": max(0.0, float(blade_px)),
        "petiole_width_px": state["petiole_width_px"],
        "petiole_width_mm": state["petiole_width_mm"],
    }
