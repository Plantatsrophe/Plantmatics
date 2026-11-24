#!/usr/bin/env python3
"""
Pure NumPy / OpenCV geometry, contour, mask, affine orientation, and bilateral
leaf measurement fallbacks for `scripts/vision/lm3_compat.py` when native
`leafmachine3` is absent from the active Python interpreter.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np


def _fb_signed_area(pts: np.ndarray) -> float:
    xy = np.asarray(pts, dtype=np.float64).reshape(-1, 2)
    if len(xy) < 3:
        return 0.0
    x = xy[:, 0]
    y = xy[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _fb_extract_contour(mask: np.ndarray) -> Optional[np.ndarray]:
    if mask is None or mask.size == 0:
        return None
    u8 = (np.asarray(mask) > 0).astype(np.uint8) * 255
    cnts, _ = cv2.findContours(u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not cnts:
        return None
    best = max(cnts, key=cv2.contourArea)
    pts = best.reshape(-1, 2).astype(np.float64)
    return pts if len(pts) >= 3 else None


def _fb_prepare_mask(
    mask: np.ndarray,
    min_area: int = 1,
    smooth_kernel: int = 0,
    pad_px: int = 0,
    fill_holes: bool = True,
    largest_only: bool = True,
    **_kwargs: Any,
) -> np.ndarray:
    u8 = (np.asarray(mask) > 0).astype(np.uint8)
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return u8.astype(bool)
    num_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(u8, connectivity=8)
    if num_lbl <= 1:
        return np.zeros_like(u8, dtype=bool)
    areas = stats[1:, cv2.CC_STAT_AREA]
    valid_ids = [1 + idx for idx, a in enumerate(areas) if int(a) >= int(min_area)]
    if not valid_ids:
        valid_ids = [1 + int(np.argmax(areas))]
    if largest_only:
        best_id = max(valid_ids, key=lambda cid: int(stats[cid, cv2.CC_STAT_AREA]))
        sel_u8 = np.where(labels == best_id, 255, 0).astype(np.uint8)
    else:
        sel_u8 = np.where(np.isin(labels, valid_ids), 255, 0).astype(np.uint8)
    if fill_holes:
        cnts, _ = cv2.findContours(sel_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        filled = np.zeros_like(sel_u8)
        if cnts:
            cv2.drawContours(filled, cnts, -1, 255, thickness=cv2.FILLED)
        sel_u8 = filled
    if smooth_kernel and int(smooth_kernel) > 1:
        k = int(smooth_kernel) | 1
        kern = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
        sel_u8 = cv2.morphologyEx(sel_u8, cv2.MORPH_CLOSE, kern)
    if pad_px and int(pad_px) > 0:
        p = int(pad_px)
        sel_u8 = np.pad(sel_u8, ((p, p), (p, p)), mode="constant", constant_values=0)
    return sel_u8 > 0


def _fb_resample_closed(pts: np.ndarray, n_points: int = 120) -> np.ndarray:
    xy = np.asarray(pts, dtype=np.float64).reshape(-1, 2)
    if len(xy) < 2 or n_points <= 0:
        return xy
    closed = np.vstack([xy, xy[0:1]])
    diffs = np.diff(closed, axis=0)
    seg_lens = np.hypot(diffs[:, 0], diffs[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total = float(cum[-1])
    if total <= 1e-9:
        return np.repeat(xy[0:1], n_points, axis=0)
    targets = np.linspace(0.0, total, int(n_points), endpoint=False)
    rx = np.interp(targets, cum, closed[:, 0])
    ry = np.interp(targets, cum, closed[:, 1])
    return np.column_stack([rx, ry])


def _fb_polygon_morphology(contour_xy: np.ndarray, method: str = "minarearect") -> Optional[Any]:
    pts = np.asarray(contour_xy, dtype=np.float32).reshape(-1, 1, 2)
    if len(pts) < 3:
        return None
    area = float(abs(cv2.contourArea(pts)))
    perim = float(cv2.arcLength(pts, True))
    hull = cv2.convexHull(pts)
    hull_area = float(abs(cv2.contourArea(hull))) if hull is not None and len(hull) >= 3 else 0.0
    solidity = (area / hull_area) if hull_area > 0.0 else 0.0
    circularity = (4.0 * math.pi * area / (perim * perim)) if perim > 0.0 else 0.0
    rect = cv2.minAreaRect(pts)
    rw, rh = float(rect[1][0]), float(rect[1][1])
    dim_max = max(rw, rh)
    dim_min = min(rw, rh)
    box_pts = cv2.boxPoints(rect).astype(np.float64)
    bx, by, bw, bh = cv2.boundingRect(pts.astype(np.int32))
    extent_rect = (area / max(dim_max * dim_min, 1.0)) if dim_max * dim_min > 0 else 0.0
    return SimpleNamespace(
        area_px=area,
        perimeter_px=perim,
        convex_area_px=hull_area,
        convex_hull_area=hull_area,
        solidity=min(max(solidity, 0.0), 1.0),
        convexity=min(max(solidity, 0.0), 1.0),
        circularity=min(max(circularity, 0.0), 1.0),
        extent=min(max(extent_rect, 0.0), 1.0),
        dim_max=dim_max,
        dim_min=dim_min,
        rotated_bbox=box_pts,
        axis_bbox=(bx, by, bw, bh),
        method=method,
    )


def _fb_components(mask: np.ndarray, min_area: int = 1, **_kwargs: Any) -> List[np.ndarray]:
    u8 = (np.asarray(mask) > 0).astype(np.uint8)
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return []
    num_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(u8, connectivity=8)
    comps: List[np.ndarray] = []
    for cid in range(1, num_lbl):
        if int(stats[cid, cv2.CC_STAT_AREA]) >= int(min_area):
            comps.append(labels == cid)
    return comps


def _fb_largest_component(mask: np.ndarray) -> np.ndarray:
    u8 = (np.asarray(mask) > 0).astype(np.uint8)
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return np.zeros_like(u8, dtype=bool)
    num_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(u8, connectivity=8)
    if num_lbl <= 1:
        return np.zeros_like(u8, dtype=bool)
    best = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
    return labels == best


def _fb_outer_polygons(mask: np.ndarray) -> List[np.ndarray]:
    u8 = (np.asarray(mask) > 0).astype(np.uint8) * 255
    cnts, _ = cv2.findContours(u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    return [c.reshape(-1, 2).astype(np.float64) for c in cnts if len(c) >= 3]


def _fb_polygon_mask(polygons: Sequence[np.ndarray], shape: Tuple[int, int]) -> np.ndarray:
    h, w = int(shape[0]), int(shape[1])
    canvas = np.zeros((h, w), dtype=np.uint8)
    valid = [np.round(np.asarray(p)).astype(np.int32).reshape(-1, 1, 2) for p in polygons if len(p) >= 3]
    if valid:
        cv2.fillPoly(canvas, valid, 1)
    return canvas > 0


def _fb_mask_bbox(mask: np.ndarray) -> Optional[Tuple[int, int, int, int]]:
    if mask is None or mask.size == 0:
        return None
    ys, xs = np.where(np.asarray(mask) > 0)
    if len(xs) == 0:
        return None
    return int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1


def _fb_crop_to_box(
    arr: Optional[np.ndarray],
    bbox: Tuple[int, int, int, int],
    pad: int = 0,
    **_kwargs: Any,
) -> Optional[np.ndarray]:
    if arr is None:
        return None
    x0, y0, x1, y1 = [int(v) for v in bbox]
    p = max(0, int(pad))
    h, w = arr.shape[:2]
    x0 = max(0, min(w, x0 - p))
    x1 = max(x0, min(w, x1 + p))
    y0 = max(0, min(h, y0 - p))
    y1 = max(y0, min(h, y1 + p))
    return arr[y0:y1, x0:x1].copy()


def _fb_oriented_affine(shape_hw: Tuple[int, int], angle_cw: float) -> Tuple[np.ndarray, Tuple[int, int]]:
    h, w = int(shape_hw[0]), int(shape_hw[1])
    ang = float(angle_cw) % 360.0
    if abs(ang) <= 1e-6 or abs(ang - 360.0) <= 1e-6:
        return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64), (h, w)
    cx, cy = (w - 1) * 0.5, (h - 1) * 0.5
    # In image coords (x right, y down), OpenCV positive angle is CCW, so -angle_cw is CW
    mat = cv2.getRotationMatrix2D((cx, cy), -ang, 1.0)
    cos_v = abs(float(mat[0, 0]))
    sin_v = abs(float(mat[0, 1]))
    new_w = max(1, int(math.ceil(w * cos_v + h * sin_v)))
    new_h = max(1, int(math.ceil(w * sin_v + h * cos_v)))
    mat[0, 2] += (new_w - 1) * 0.5 - cx
    mat[1, 2] += (new_h - 1) * 0.5 - cy
    return mat.astype(np.float64), (new_h, new_w)


def _fb_apply_affine(mat: np.ndarray, pts: np.ndarray) -> np.ndarray:
    arr = np.asarray(pts, dtype=np.float64)
    if arr.size == 0:
        return arr.reshape(0, 2)
    m = np.asarray(mat, dtype=np.float64)
    xy = arr.reshape(-1, 2)
    out = (xy @ m[:, :2].T) + m[:, 2]
    return out.reshape(arr.shape)


def _fb_rotate_image(
    img: np.ndarray,
    angle_cw: float,
    bg: Any = 0,
    nearest: bool = True,
) -> np.ndarray:
    h, w = int(img.shape[0]), int(img.shape[1])
    mat, (new_h, new_w) = _fb_oriented_affine((h, w), angle_cw)
    interp = cv2.INTER_NEAREST if nearest else cv2.INTER_LINEAR
    border_val = (int(bg),) * (img.shape[2] if img.ndim == 3 else 1) if np.isscalar(bg) else tuple(bg)
    return cv2.warpAffine(
        img,
        mat,
        (new_w, new_h),
        flags=interp,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=border_val,
    )


def _fb_compute_orientation(
    keypoints_or_mask: Any = None,
    p_base: Optional[Tuple[float, float]] = None,
    p_apex: Optional[Tuple[float, float]] = None,
    **_kwargs: Any,
) -> Any:
    if isinstance(keypoints_or_mask, dict):
        kp = keypoints_or_mask
        if p_base is None:
            for k in ("lamina_base", "base_center", "midvein_14"):
                if k in kp and kp[k] is not None and len(kp[k]) >= 2:
                    p_base = (float(kp[k][0]), float(kp[k][1]))
                    break
        if p_apex is None:
            for k in ("lamina_tip", "apex_center", "midvein_0"):
                if k in kp and kp[k] is not None and len(kp[k]) >= 2:
                    p_apex = (float(kp[k][0]), float(kp[k][1]))
                    break
    if p_base is not None and p_apex is not None:
        dx = float(p_apex[0] - p_base[0])
        dy = float(p_apex[1] - p_base[1])
        if math.hypot(dx, dy) > 1e-6:
            angle_cw = (-math.degrees(math.atan2(dx, -dy))) % 360.0
            return SimpleNamespace(
                angle_cw=angle_cw,
                angle_deg=math.degrees(math.atan2(dy, dx)),
                p_base=p_base,
                p_apex=p_apex,
                success=True,
                method="landmarks",
            )
    if isinstance(keypoints_or_mask, np.ndarray) and keypoints_or_mask.ndim >= 2:
        pts = _fb_extract_contour(keypoints_or_mask)
        if pts is not None and len(pts) >= 3:
            m = cv2.moments((keypoints_or_mask > 0).astype(np.uint8), binaryImage=True)
            if m["m00"] > 0:
                theta = 0.5 * math.atan2(2.0 * m["mu11"], m["mu20"] - m["mu02"])
                vx, vy = math.cos(theta), math.sin(theta)
                angle_cw = (-math.degrees(math.atan2(vx, -vy))) % 360.0
                return SimpleNamespace(
                    angle_cw=angle_cw,
                    angle_deg=math.degrees(theta),
                    p_base=None,
                    p_apex=None,
                    success=True,
                    method="pca_moments",
                )
    return SimpleNamespace(angle_cw=0.0, angle_deg=-90.0, p_base=p_base, p_apex=p_apex, success=False, method="fallback")


def _fb_length_width_from_box(rotated_bbox: np.ndarray, angle_cw: float) -> Tuple[float, float]:
    pts = np.asarray(rotated_bbox, dtype=np.float64).reshape(-1, 2)
    rad = math.radians(float(angle_cw))
    cos_a, sin_a = math.cos(rad), math.sin(rad)
    xr = pts[:, 0] * cos_a - pts[:, 1] * sin_a
    yr = pts[:, 0] * sin_a + pts[:, 1] * cos_a
    length_px = float(np.max(yr) - np.min(yr))
    width_px = float(np.max(xr) - np.min(xr))
    return max(length_px, 0.0), max(width_px, 0.0)


@dataclass
class _FbOrientedLeaf:
    specimen_id: int
    leaf_id: int
    detection_id: int
    stem: str
    silhouette: np.ndarray
    lamina: np.ndarray
    holes: np.ndarray
    petiole: Optional[np.ndarray]
    kpts: Dict[str, Tuple[float, float]]
    kpt_conf: Dict[str, float]
    angle_cw: float
    det_box: Tuple[int, int, int, int]
    crop_shape: Tuple[int, int]
    truncated: bool = False


def _fb_measure_leaf(oleaf: Any, min_score: float = 0.70) -> Optional[Any]:
    lam = np.asarray(oleaf.lamina, dtype=bool)
    if lam.size == 0 or int(np.count_nonzero(lam)) == 0:
        return None
    h, w = lam.shape[:2]
    num_lbl, _, stats, _ = cv2.connectedComponentsWithStats(lam.astype(np.uint8), connectivity=8)
    total_fg = float(np.count_nonzero(lam))
    largest_fg = float(np.max(stats[1:, cv2.CC_STAT_AREA])) if num_lbl > 1 else total_fg
    largest_frac = largest_fg / max(total_fg, 1.0)

    kpts = getattr(oleaf, "kpts", None) or {}
    mid_xs: List[float] = []
    for k in ("lamina_tip", "lamina_base"):
        if k in kpts and kpts[k] is not None:
            mid_xs.append(float(kpts[k][0]))
    for i in range(15):
        mk = f"midvein_{i}"
        if mk in kpts and kpts[mk] is not None:
            mid_xs.append(float(kpts[mk][0]))

    ys_fg, xs_fg = np.where(lam)
    x_axis = float(np.median(mid_xs)) if mid_xs else float(np.mean(xs_fg))

    # Reflect foreground pixels across vertical midrib axis x = x_axis
    xs_refl = np.round(2.0 * x_axis - xs_fg.astype(np.float64)).astype(np.int64)
    in_bounds = (xs_refl >= 0) & (xs_refl < w)
    refl_mask = np.zeros_like(lam, dtype=bool)
    refl_mask[ys_fg[in_bounds], xs_refl[in_bounds]] = True

    inter = float(np.count_nonzero(lam & refl_mask))
    denom = float(np.count_nonzero(lam)) + float(np.count_nonzero(refl_mask))
    dice = (2.0 * inter / denom) if denom > 0.0 else 0.0

    area_left = float(np.count_nonzero(xs_fg < x_axis))
    area_right = float(np.count_nonzero(xs_fg > x_axis))
    lr_sum = area_left + area_right
    si_a = ((area_right - area_left) / lr_sum) if lr_sum > 0.0 else 0.0

    archetype_score = float(
        np.clip(0.65 * dice + 0.35 * (1.0 - min(1.0, abs(si_a))) * largest_frac, 0.0, 1.0)
    )
    gates_pass = bool(dice >= 0.70 and abs(si_a) < 0.35 and largest_frac >= 0.85)
    is_archetypal = bool(gates_pass and archetype_score >= float(min_score))
    return SimpleNamespace(
        dice=dice,
        si_a=si_a,
        largest_frac=largest_frac,
        archetype_score=archetype_score,
        gates_pass=gates_pass,
        is_archetypal=is_archetypal,
    )


__all__ = [
    "_fb_signed_area", "_fb_extract_contour", "_fb_prepare_mask", "_fb_resample_closed",
    "_fb_polygon_morphology", "_fb_components", "_fb_largest_component", "_fb_outer_polygons",
    "_fb_polygon_mask", "_fb_mask_bbox", "_fb_crop_to_box", "_fb_oriented_affine",
    "_fb_apply_affine", "_fb_rotate_image", "_fb_compute_orientation",
    "_fb_length_width_from_box", "_FbOrientedLeaf", "_fb_measure_leaf",
]
