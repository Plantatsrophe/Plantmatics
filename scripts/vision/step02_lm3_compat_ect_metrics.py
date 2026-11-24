#!/usr/bin/env python3
"""
Pure NumPy / OpenCV Euler Characteristic Transform (ECT), radial ECT overlay,
leaf product rendering, petiole measurement, and landmark metric fallbacks for
`scripts/vision/lm3_compat.py` when native `leafmachine3` is absent.
"""

from __future__ import annotations

import math
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

from scripts.vision.step02_lm3_compat_geometry import (
    _fb_extract_contour,
    _fb_resample_closed,
    _fb_rotate_image,
)


def _fb_ect_product_for(include_petiole: bool = False, include_holes: bool = True) -> str:
    if include_petiole and include_holes:
        return "laminapet"
    if include_petiole and not include_holes:
        return "silhouettepet"
    if not include_petiole and include_holes:
        return "lamina"
    return "silhouette"


def _fb_render_leaf_products(
    _rgb: np.ndarray,
    masks_dict: Dict[str, Optional[np.ndarray]],
    angle_cw: Optional[float] = None,
    want: Optional[Sequence[str]] = None,
) -> Dict[str, Optional[np.ndarray]]:
    keys = list(want) if want is not None else list(masks_dict.keys())
    out: Dict[str, Optional[np.ndarray]] = {}
    ang = float(angle_cw) if angle_cw is not None else 0.0
    for k in keys:
        m = masks_dict.get(k)
        if m is None:
            out[k] = None
            continue
        u8 = (np.asarray(m) > 0).astype(np.uint8) * 255
        if abs(ang % 360.0) > 0.25:
            u8 = _fb_rotate_image(u8, ang, bg=0, nearest=True)
        out[k] = u8 > 0
    return out


def _fb_compute_ect(
    mask_u8: np.ndarray,
    num_dirs: int = 32,
    bound_radius: float = 1.0,
    want_simple: bool = False,
    num_thresh: int = 32,
) -> Optional[Any]:
    u8 = (np.asarray(mask_u8) > 0).astype(np.uint8) * 255
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return None
    cnts_raw, hier = cv2.findContours(u8, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    if not cnts_raw or hier is None:
        return None
    cycles: List[Tuple[float, np.ndarray]] = []
    for idx_c, c_raw in enumerate(cnts_raw):
        c_xy = c_raw.reshape(-1, 2).astype(np.float64)
        if len(c_xy) < 3 or abs(float(cv2.contourArea(c_raw))) < 4.0:
            continue
        sign_c = 1.0 if int(hier[0][idx_c][3]) < 0 else -1.0
        n_res = min(256, max(48, len(c_xy)))
        cycles.append((sign_c, _fb_resample_closed(c_xy, n_points=n_res)))
    if not cycles:
        cnt_xy = _fb_extract_contour(u8)
        if cnt_xy is None or len(cnt_xy) < 3:
            return None
        cycles = [(1.0, _fb_resample_closed(cnt_xy, n_points=min(256, max(64, len(cnt_xy)))))]

    outer_stack = np.vstack([pts_c for s_c, pts_c in cycles if s_c > 0] or [cycles[0][1]])
    centroid = np.mean(outer_stack, axis=0)
    centered_outer = outer_stack - centroid
    max_r = float(np.max(np.hypot(centered_outer[:, 0], centered_outer[:, 1])))
    scale = (float(bound_radius) * 0.95 / max_r) if max_r > 1e-9 else 1.0
    outline_norm = ((cycles[0][1] - centroid) * scale).astype(np.float32)

    n_dirs = max(4, int(num_dirs))
    n_thr = max(4, int(num_thresh))
    thetas = np.linspace(0.0, 2.0 * np.pi, n_dirs, endpoint=False, dtype=np.float32)
    thresholds = np.linspace(-float(bound_radius), float(bound_radius), n_thr, dtype=np.float32)

    # Multi-component polygonal Euler Characteristic Transform: sum_c sign_c * (#V_c(<=t) - #E_c(<=t))
    dirs = np.stack([np.cos(thetas), np.sin(thetas)], axis=1)  # (D, 2)
    ect_acc = np.zeros((n_thr, n_dirs), dtype=np.float32)
    for sign_c, pts_c in cycles:
        c_norm = ((pts_c - centroid) * scale).astype(np.float32)
        heights = c_norm @ dirs.T  # (N_c, D)
        heights_next = np.roll(heights, -1, axis=0)
        edge_max = np.maximum(heights, heights_next)
        v_le = (heights[None, :, :] <= thresholds[:, None, None]).sum(axis=1)
        e_le = (edge_max[None, :, :] <= thresholds[:, None, None]).sum(axis=1)
        ect_acc += np.float32(sign_c) * (v_le - e_le).astype(np.float32)
    ect_mat = ect_acc.astype(np.float32)
    return SimpleNamespace(
        ect_matrix=ect_mat,
        thetas=thetas,
        thresholds=thresholds,
        bound_radius=float(bound_radius),
        outline_norm=outline_norm,
        want_simple=bool(want_simple),
    )


def _fb_tip_up_directions(ect_raw: np.ndarray, thetas: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mat = np.asarray(ect_raw, dtype=np.float32)
    th = np.asarray(thetas, dtype=np.float32)
    if mat.ndim != 2 or len(th) != mat.shape[1] or len(th) == 0:
        return mat, th
    # Tip-up direction in image coordinates is (0, -1) -> angle 3*pi/2 (or -pi/2)
    target = 1.5 * math.pi
    diffs = np.abs(np.angle(np.exp(1j * (th - target))))
    shift = int(np.argmin(diffs))
    return np.roll(mat, -shift, axis=1), np.roll(th, -shift)


def _fb_visual_log(mat: np.ndarray) -> np.ndarray:
    arr = np.asarray(mat, dtype=np.float32)
    return (np.sign(arr) * np.log1p(np.abs(arr))).astype(np.float32)


def _fb_render_radial_ect_overlay(
    ect_vis: np.ndarray,
    thetas: np.ndarray,
    thresholds: np.ndarray,
    bound_radius: float,
    outline_norm: Optional[np.ndarray] = None,
    cmap: str = "viridis",
    size: int = 256,
) -> np.ndarray:
    canvas = np.full((size, size, 3), 248, dtype=np.uint8)
    mat = np.asarray(ect_vis, dtype=np.float32)
    if mat.size > 0:
        yy, xx = np.mgrid[0:size, 0:size].astype(np.float32)
        cx = cy = (size - 1) * 0.5
        rad_max = size * 0.44
        dx = xx - cx
        dy = yy - cy
        rr = np.hypot(dx, dy) / max(rad_max, 1.0)
        inside = rr <= 1.0
        ang = np.mod(np.arctan2(dy, dx), 2.0 * np.pi)
        t_idx = np.clip((rr * (mat.shape[0] - 1)).astype(np.int32), 0, mat.shape[0] - 1)
        d_idx = np.clip((ang / (2.0 * np.pi) * mat.shape[1]).astype(np.int32), 0, mat.shape[1] - 1)
        vals = mat[t_idx, d_idx]
        vmin, vmax = float(np.min(vals)), float(np.max(vals))
        norm_u8 = (
            np.round((vals - vmin) / max(vmax - vmin, 1e-6) * 255.0).astype(np.uint8)
            if vmax > vmin
            else np.full_like(t_idx, 128, dtype=np.uint8)
        )
        cv_cmap = getattr(cv2, f"COLORMAP_{str(cmap).upper()}", cv2.COLORMAP_VIRIDIS)
        colored = cv2.applyColorMap(norm_u8, cv_cmap)
        canvas[inside] = colored[inside]
    if outline_norm is not None and len(outline_norm) >= 3:
        cx = cy = (size - 1) * 0.5
        rad_max = size * 0.44 / max(float(bound_radius), 1e-6)
        pts_px = np.round(np.asarray(outline_norm, dtype=np.float32) * rad_max + np.array([cx, cy])).astype(np.int32)
        cv2.polylines(canvas, [pts_px.reshape(-1, 1, 2)], True, (255, 255, 255), 2, cv2.LINE_AA)
    return canvas


def _fb_measure_petiole(
    pet_u8: np.ndarray,
    lam_u8: np.ndarray,
    pet_chain: Sequence[Tuple[float, float]],
) -> Optional[Any]:
    p_bin = (np.asarray(pet_u8) > 0).astype(np.uint8)
    l_bin = (np.asarray(lam_u8) > 0).astype(np.uint8)
    if p_bin.size == 0 or int(np.count_nonzero(p_bin)) == 0:
        return None
    arc_len = 0.0
    for i in range(len(pet_chain) - 1):
        arc_len += math.hypot(
            float(pet_chain[i + 1][0]) - float(pet_chain[i][0]),
            float(pet_chain[i + 1][1]) - float(pet_chain[i][1]),
        )
    dt = cv2.distanceTransform(p_bin, cv2.DIST_L2, 3)
    h, w = p_bin.shape[:2]
    radii: List[float] = []
    for px, py in pet_chain:
        ix = max(0, min(w - 1, int(round(float(px)))))
        iy = max(0, min(h - 1, int(round(float(py)))))
        if dt[iy, ix] > 0.0:
            radii.append(2.0 * float(dt[iy, ix]))
    if not radii:
        nz = dt[dt > 0]
        width_px = float(2.0 * np.median(nz)) if nz.size > 0 else 1.0
    else:
        width_px = float(np.median(radii))
    touches = bool(
        l_bin.shape == p_bin.shape
        and np.any(cv2.dilate(p_bin, np.ones((3, 3), dtype=np.uint8)) & l_bin)
    )
    return SimpleNamespace(width_px=width_px, length_px=arc_len, touches_leaf=touches)


def _fb_compute_measurements(pts_2d: Dict[str, Tuple[float, float]]) -> Optional[Any]:
    if not pts_2d:
        return None
    lt = pts_2d.get("lamina_tip") or pts_2d.get("midvein_0")
    lb = pts_2d.get("lamina_base") or pts_2d.get("midvein_14")
    tip_base_len = (
        math.hypot(float(lt[0]) - float(lb[0]), float(lt[1]) - float(lb[1]))
        if (lt is not None and lb is not None)
        else None
    )
    mv_pts: List[Tuple[float, float]] = []
    if lt is not None:
        mv_pts.append((float(lt[0]), float(lt[1])))
    for i in range(15):
        mk = f"midvein_{i}"
        if mk in pts_2d and pts_2d[mk] is not None:
            pt = (float(pts_2d[mk][0]), float(pts_2d[mk][1]))
            if not mv_pts or math.hypot(pt[0] - mv_pts[-1][0], pt[1] - mv_pts[-1][1]) > 1e-6:
                mv_pts.append(pt)
    if lb is not None:
        pt_b = (float(lb[0]), float(lb[1]))
        if not mv_pts or math.hypot(pt_b[0] - mv_pts[-1][0], pt_b[1] - mv_pts[-1][1]) > 1e-6:
            mv_pts.append(pt_b)
    lam_trace = (
        sum(
            math.hypot(mv_pts[i + 1][0] - mv_pts[i][0], mv_pts[i + 1][1] - mv_pts[i][1])
            for i in range(len(mv_pts) - 1)
        )
        if len(mv_pts) >= 2
        else tip_base_len
    )
    pet_pts: List[Tuple[float, float]] = []
    if lb is not None:
        pet_pts.append((float(lb[0]), float(lb[1])))
    for k in ("petiole_0", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_tip"):
        if k in pts_2d and pts_2d[k] is not None:
            pt = (float(pts_2d[k][0]), float(pts_2d[k][1]))
            if not pet_pts or math.hypot(pt[0] - pet_pts[-1][0], pt[1] - pet_pts[-1][1]) > 1e-6:
                pet_pts.append(pt)
    pet_trace = (
        sum(
            math.hypot(pet_pts[i + 1][0] - pet_pts[i][0], pet_pts[i + 1][1] - pet_pts[i][1])
            for i in range(len(pet_pts) - 1)
        )
        if len(pet_pts) >= 2
        else None
    )
    w1 = pts_2d.get("width_1") or pts_2d.get("width_left")
    w2 = pts_2d.get("width_2") or pts_2d.get("width_right")
    leaf_w = (
        math.hypot(float(w2[0]) - float(w1[0]), float(w2[1]) - float(w1[1]))
        if (w1 is not None and w2 is not None)
        else None
    )

    def _interior_angle_deg(
        vertex: Optional[Tuple[float, float]],
        left: Optional[Tuple[float, float]],
        right: Optional[Tuple[float, float]],
    ) -> Optional[float]:
        if vertex is None or left is None or right is None:
            return None
        v1x = float(left[0]) - float(vertex[0])
        v1y = float(left[1]) - float(vertex[1])
        v2x = float(right[0]) - float(vertex[0])
        v2y = float(right[1]) - float(vertex[1])
        n1 = math.hypot(v1x, v1y)
        n2 = math.hypot(v2x, v2y)
        if n1 <= 1e-9 or n2 <= 1e-9:
            return None
        dot = v1x * v2x + v1y * v2y
        cos_val = float(np.clip(dot / (n1 * n2), -1.0, 1.0))
        return float(math.degrees(math.acos(cos_val)))

    a_left = pts_2d.get("apex_left") or pts_2d.get("apex_1")
    a_right = pts_2d.get("apex_right") or pts_2d.get("apex_3")
    a_center = pts_2d.get("apex_center") or pts_2d.get("apex_2") or pts_2d.get("lamina_tip") or lt
    if a_left is not None and a_right is not None and a_center is not None:
        apex_ang = _interior_angle_deg(a_center, a_left, a_right)
    elif lt is not None and w1 is not None and w2 is not None:
        apex_ang = _interior_angle_deg(lt, w1, w2)
    else:
        apex_ang = None

    b_left = pts_2d.get("base_left") or pts_2d.get("base_1")
    b_right = pts_2d.get("base_right") or pts_2d.get("base_3")
    b_center = pts_2d.get("base_center") or pts_2d.get("base_2") or pts_2d.get("lamina_base") or lb
    if b_left is not None and b_right is not None and b_center is not None:
        base_ang = _interior_angle_deg(b_center, b_left, b_right)
    elif lb is not None and w1 is not None and w2 is not None:
        base_ang = _interior_angle_deg(lb, w1, w2)
    else:
        base_ang = None

    curv = (lam_trace / tip_base_len) if (lam_trace and tip_base_len and tip_base_len > 1e-6) else 1.0
    return SimpleNamespace(
        lamina_trace_length=lam_trace,
        lamina_tip_base_length=tip_base_len,
        petiole_trace_length=pet_trace,
        leaf_width=leaf_w,
        apex_angle=apex_ang,
        base_angle=base_ang,
        lamina_trace_length_px=lam_trace,
        lamina_tip_base_length_px=tip_base_len,
        petiole_trace_length_px=pet_trace,
        leaf_width_px=leaf_w,
        apex_angle_deg=apex_ang,
        base_angle_deg=base_ang,
        lamina_curvature=curv,
    )


__all__ = [
    "_fb_ect_product_for", "_fb_render_leaf_products", "_fb_compute_ect",
    "_fb_tip_up_directions", "_fb_visual_log", "_fb_render_radial_ect_overlay",
    "_fb_measure_petiole", "_fb_compute_measurements",
]
