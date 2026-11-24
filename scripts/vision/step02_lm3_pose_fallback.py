#!/usr/bin/env python3
"""
Step 02 LeafMachine3 pose result container, orthogonal petiole-width estimation,
and deterministic geometric fallback for 31-keypoint `mid15_pet5` skeletons.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np

import scripts.vision.lm3_compat  # noqa: F401
from scripts.vision.step02_lm3_petiole_kpts import normalize_lm3_keypoints_dict

_KEYPOINT_NAMES_31: Tuple[str, ...] = (
    "lamina_tip", "apex_1", "apex_2", "apex_3",
    "midvein_1", "midvein_2", "midvein_3", "midvein_4", "midvein_5",
    "midvein_6", "midvein_7", "midvein_8", "midvein_9", "midvein_10",
    "midvein_11", "midvein_12", "midvein_13", "midvein_14", "midvein_15",
    "base_1", "base_2", "base_3", "lamina_base",
    "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_5", "petiole_6",
    "width_1", "width_2",
)


class LM3PoseResult(dict):
    """
    Dictionary mapping all 31 canonical `mid15_pet5` keypoint names (`lamina_tip`, `apex_1..3`,
    `midvein_1..15`, `base_1..3`, `lamina_base`, `petiole_1..6`, `width_1..2`) to
    `(x, y, confidence)` tuples while also exposing `.petiole_width_px` and `.keypoints_31`.
    """

    def __init__(
        self,
        keypoints_31: Dict[str, Tuple[float, float, float]],
        petiole_width_px: Optional[float] = None,
        is_geometric_fallback: bool = False,
    ) -> None:
        super().__init__(keypoints_31)
        self.petiole_width_px: Optional[float] = (
            float(petiole_width_px) if (petiole_width_px is not None and np.isfinite(petiole_width_px)) else None
        )
        self.keypoints_31: Dict[str, Tuple[float, float, float]] = dict(keypoints_31)
        self.is_geometric_fallback: bool = bool(is_geometric_fallback)

    def __getitem__(self, key: str) -> Any:
        if key == "petiole_width_px":
            return self.petiole_width_px
        if key == "keypoints_31":
            return self.keypoints_31
        return super().__getitem__(key)

    def get(self, key: str, default: Any = None) -> Any:
        if key == "petiole_width_px":
            return self.petiole_width_px if self.petiole_width_px is not None else default
        if key == "keypoints_31":
            return self.keypoints_31
        return super().get(key, default)


def estimate_petiole_width_px(
    keypoints_31: Dict[str, Tuple[float, float, float]],
    petiole_mask: Optional[np.ndarray] = None,
    lamina_mask: Optional[np.ndarray] = None,
) -> Optional[float]:
    """Estimates petiole width in pixels by delegating to `leafmachine3.core.petiole.measure_petiole`."""
    if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
        pet_u8 = (petiole_mask > 0).astype(np.uint8)
        lam_u8 = (
            (lamina_mask > 0).astype(np.uint8)
            if (lamina_mask is not None and lamina_mask.size > 0 and lamina_mask.shape[:2] == pet_u8.shape[:2])
            else None
        )
        norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
        petiole_pts: List[Tuple[float, float]] = []
        for kp_name in ("lamina_base", "petiole_0", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_tip"):
            kp = norm_kps.get(kp_name)
            if kp is not None and len(kp) >= 3 and float(kp[2]) >= 0.20:
                pt_xy = (float(kp[0]), float(kp[1]))
                if not petiole_pts or math.hypot(pt_xy[0] - petiole_pts[-1][0], pt_xy[1] - petiole_pts[-1][1]) > 1e-3:
                    petiole_pts.append(pt_xy)
        if len(petiole_pts) < 2:
            ys_p, xs_p = np.where(pet_u8 > 0)
            if len(xs_p) >= 2:
                pts_cov = np.column_stack((xs_p.astype(np.float64), ys_p.astype(np.float64)))
                mean_pt = pts_cov.mean(axis=0)
                centered = pts_cov - mean_pt
                _, _, vh = np.linalg.svd(centered, full_matrices=False)
                proj = centered @ vh[0]
                i_min = int(np.argmin(proj))
                i_max = int(np.argmax(proj))
                p0 = (float(pts_cov[i_min, 0]), float(pts_cov[i_min, 1]))
                p1 = (float(pts_cov[i_max, 0]), float(pts_cov[i_max, 1]))
                if math.hypot(p1[0] - p0[0], p1[1] - p0[1]) > 1e-3:
                    petiole_pts = [p0, p1]
        if lam_u8 is not None and len(petiole_pts) >= 2:
            try:
                from leafmachine3.core import petiole as petiole_metrics

                res = petiole_metrics.measure_petiole(pet_u8, lam_u8, petiole_pts)
                if res is not None and res.width_px is not None and np.isfinite(res.width_px) and res.width_px > 0:
                    return round(float(res.width_px), 3)
            except Exception:
                pass

        dt = cv2.distanceTransform(pet_u8, cv2.DIST_L2, 3)
        sampled_radii: List[float] = []
        h_p, w_p = pet_u8.shape[:2]
        for pt_x, pt_y in petiole_pts:
            ix = max(0, min(w_p - 1, int(round(pt_x))))
            iy = max(0, min(h_p - 1, int(round(pt_y))))
            r_val = float(dt[iy, ix])
            if r_val > 0.0:
                sampled_radii.append(2.0 * r_val)
        if sampled_radii:
            return round(float(np.median(sampled_radii)), 3)
        nz_dt = dt[dt > 0]
        if nz_dt.size > 0:
            return round(float(2.0 * np.percentile(nz_dt, 75)), 3)

    if lamina_mask is not None and lamina_mask.size > 0 and int(np.count_nonzero(lamina_mask)) > 0:
        lam_u8 = (lamina_mask > 0).astype(np.uint8)
        dt_lam = cv2.distanceTransform(lam_u8, cv2.DIST_L2, 3)
        ys_l, _ = np.where(lam_u8 > 0)
        if ys_l.size > 0:
            y_cut = float(np.percentile(ys_l, 85))
            basal_dt = dt_lam[(lam_u8 > 0) & (np.arange(lam_u8.shape[0])[:, None] >= y_cut)]
            if basal_dt.size > 0:
                sw = float(2.0 * np.median(basal_dt))
                if sw > 0.0:
                    return round(sw, 3)

    return None


def geometric_fallback_keypoints_31(
    crop_bgr: Optional[np.ndarray],
    lamina_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    offset_xy: Tuple[float, float] = (0.0, 0.0),
    crown_hint_xy: Optional[Tuple[float, float]] = None,
) -> LM3PoseResult:
    """
    Derives the 31 canonical `mid15_pet5` keypoints deterministically from `lamina_mask`
    and `petiole_mask` (or `crop_bgr`) when ONNX weights are absent in unit tests,
    anchoring reniform/cordate (`W > L`) blades along the petiole-to-centroid axis
    rather than rotating `90 deg` onto lateral lobes.
    """
    ox, oy = float(offset_xy[0]), float(offset_xy[1])
    mask_use = lamina_mask
    if (mask_use is None or mask_use.size == 0) and crop_bgr is not None and crop_bgr.size > 0:
        gray = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2GRAY) if crop_bgr.ndim == 3 else crop_bgr
        _, mask_use = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        if int(np.count_nonzero(mask_use)) == 0:
            mask_use = (gray > 0).astype(np.uint8) * 255

    if mask_use is None or mask_use.size == 0 or int(np.count_nonzero(mask_use)) == 0:
        empty_kps = {k: (ox, oy, 0.0) for k in _KEYPOINT_NAMES_31}
        return LM3PoseResult(empty_kps, petiole_width_px=None, is_geometric_fallback=True)

    ys, xs = np.where(mask_use > 0)
    pts = np.column_stack((xs.astype(np.float64), ys.astype(np.float64)))
    centroid = np.mean(pts, axis=0)
    centered = pts - centroid
    cov = np.cov(centered, rowvar=False) if len(pts) >= 3 else np.eye(2)
    eigvals, eigvecs = np.linalg.eigh(cov)
    major_vec = eigvecs[:, int(np.argmax(eigvals))]

    proj_major = centered @ major_vec
    p_min = centroid + float(np.min(proj_major)) * major_vec
    p_max = centroid + float(np.max(proj_major)) * major_vec

    # Orient base -> tip along petiole attachment axis (preserving reniform W > L blades),
    # else crown_hint_xy proximity, else apex/base taper & mount orientation
    if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
        pys, pxs = np.where(petiole_mask > 0)
        pet_cent = np.array([float(np.mean(pxs)), float(np.mean(pys))], dtype=np.float64)
        junc_pt = pet_cent
        if petiole_mask.shape[:2] == mask_use.shape[:2]:
            kern = np.ones((5, 5), dtype=np.uint8)
            pet_dil = cv2.dilate((petiole_mask > 0).astype(np.uint8), kern, iterations=1)
            jys, jxs = np.where((pet_dil > 0) & (mask_use > 0))
            if len(jxs) > 0:
                junc_pt = np.array([float(np.mean(jxs)), float(np.mean(jys))], dtype=np.float64)
        v_from_pet = centroid - junc_pt
        v_norm = float(np.linalg.norm(v_from_pet))
        if v_norm > 1.0:
            u_pet = v_from_pet / v_norm
            e0 = eigvecs[:, 0]
            e1 = eigvecs[:, 1]
            d0 = abs(float(np.dot(e0, u_pet)))
            d1 = abs(float(np.dot(e1, u_pet)))
            best_e = e0 if d0 >= d1 else e1
            if max(d0, d1) >= 0.9063:
                axis_dir = best_e if float(np.dot(best_e, u_pet)) >= 0.0 else -best_e
            else:
                axis_dir = u_pet
            proj_ax = centered @ axis_dir
            base_pt = centroid + float(np.min(proj_ax)) * axis_dir
            tip_pt = centroid + float(np.max(proj_ax)) * axis_dir
        elif np.linalg.norm(p_min - pet_cent) <= np.linalg.norm(p_max - pet_cent):
            base_pt, tip_pt = p_min, p_max
        else:
            base_pt, tip_pt = p_max, p_min
    elif crown_hint_xy is not None:
        cr_local = np.array([float(crown_hint_xy[0]) - ox, float(crown_hint_xy[1]) - oy], dtype=np.float64)
        if np.linalg.norm(p_min - cr_local) <= np.linalg.norm(p_max - cr_local):
            base_pt, tip_pt = p_min, p_max
        else:
            base_pt, tip_pt = p_max, p_min
    else:
        dt_fb = cv2.distanceTransform((mask_use > 0).astype(np.uint8), cv2.DIST_L2, 5)
        h_fb, w_fb = dt_fb.shape[:2]
        span_fb = float(np.linalg.norm(p_max - p_min))
        win_fb = max(2, min(6, int(round(0.015 * span_fb))))
        r_min_vals: List[float] = []
        r_max_vals: List[float] = []
        for frac_fb in (0.01, 0.02, 0.03, 0.04, 0.05):
            pt_a = p_min + frac_fb * (p_max - p_min)
            pt_b = p_max + frac_fb * (p_min - p_max)
            xa, ya = int(round(pt_a[0])), int(round(pt_a[1]))
            xb, yb = int(round(pt_b[0])), int(round(pt_b[1]))
            pa = dt_fb[max(0, ya - win_fb) : min(h_fb, ya + win_fb + 1), max(0, xa - win_fb) : min(w_fb, xa + win_fb + 1)]
            pb = dt_fb[max(0, yb - win_fb) : min(h_fb, yb + win_fb + 1), max(0, xb - win_fb) : min(w_fb, xb + win_fb + 1)]
            if pa.size > 0:
                r_min_vals.append(float(np.max(pa)))
            if pb.size > 0:
                r_max_vals.append(float(np.max(pb)))
        r_min_m = float(np.mean(r_min_vals)) if r_min_vals else 0.0
        r_max_m = float(np.mean(r_max_vals)) if r_max_vals else 0.0
        if abs(r_min_m - r_max_m) > 0.35:
            base_pt, tip_pt = (p_min, p_max) if r_min_m < r_max_m else (p_max, p_min)
        elif p_min[1] >= p_max[1]:
            base_pt, tip_pt = p_min, p_max
        else:
            base_pt, tip_pt = p_max, p_min

    axis_vec = tip_pt - base_pt
    axis_len = max(float(np.linalg.norm(axis_vec)), 1.0)
    axis_unit = axis_vec / axis_len
    perp_unit = np.array([-axis_unit[1], axis_unit[0]], dtype=np.float64)

    proj_perp = (pts - centroid) @ perp_unit
    w1_pt = centroid + float(np.min(proj_perp)) * perp_unit
    w2_pt = centroid + float(np.max(proj_perp)) * perp_unit
    half_w = max(1.0, 0.5 * float(np.linalg.norm(w2_pt - w1_pt)))

    kps: Dict[str, Tuple[float, float, float]] = {}
    kps["lamina_tip"] = (round(float(tip_pt[0] + ox), 2), round(float(tip_pt[1] + oy), 2), 0.92)
    kps["apex_1"] = (
        round(float(tip_pt[0] - 0.06 * axis_vec[0] - 0.18 * half_w * perp_unit[0] + ox), 2),
        round(float(tip_pt[1] - 0.06 * axis_vec[1] - 0.18 * half_w * perp_unit[1] + oy), 2),
        0.88,
    )
    kps["apex_2"] = (round(float(tip_pt[0] + ox), 2), round(float(tip_pt[1] + oy), 2), 0.90)
    kps["apex_3"] = (
        round(float(tip_pt[0] - 0.06 * axis_vec[0] + 0.18 * half_w * perp_unit[0] + ox), 2),
        round(float(tip_pt[1] - 0.06 * axis_vec[1] + 0.18 * half_w * perp_unit[1] + oy), 2),
        0.88,
    )

    # midvein_1 (near apex) .. midvein_15 (near base)
    for idx_m in range(1, 16):
        frac_from_tip = float(idx_m) / 16.0
        m_pt = tip_pt - frac_from_tip * axis_vec
        kps[f"midvein_{idx_m}"] = (round(float(m_pt[0] + ox), 2), round(float(m_pt[1] + oy), 2), 0.90)

    kps["base_1"] = (
        round(float(base_pt[0] + 0.06 * axis_vec[0] - 0.18 * half_w * perp_unit[0] + ox), 2),
        round(float(base_pt[1] + 0.06 * axis_vec[1] - 0.18 * half_w * perp_unit[1] + oy), 2),
        0.88,
    )
    kps["base_2"] = (round(float(base_pt[0] + ox), 2), round(float(base_pt[1] + oy), 2), 0.90)
    kps["base_3"] = (
        round(float(base_pt[0] + 0.06 * axis_vec[0] + 0.18 * half_w * perp_unit[0] + ox), 2),
        round(float(base_pt[1] + 0.06 * axis_vec[1] + 0.18 * half_w * perp_unit[1] + oy), 2),
        0.88,
    )
    kps["lamina_base"] = (round(float(base_pt[0] + ox), 2), round(float(base_pt[1] + oy), 2), 0.92)

    # Petiole 1..6 trace
    if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
        pys, pxs = np.where(petiole_mask > 0)
        p_pts = np.column_stack((pxs.astype(np.float64), pys.astype(np.float64)))
        dists_from_base = np.linalg.norm(p_pts - base_pt[None, :], axis=1)
        far_pet_pt = p_pts[int(np.argmax(dists_from_base))]
        for idx_p in range(1, 7):
            t_p = float(idx_p) / 6.0
            p_coord = base_pt + t_p * (far_pet_pt - base_pt)
            kps[f"petiole_{idx_p}"] = (round(float(p_coord[0] + ox), 2), round(float(p_coord[1] + oy), 2), 0.86)
    else:
        for idx_p in range(1, 7):
            t_p = float(idx_p) / 6.0
            p_coord = base_pt - (0.12 * t_p) * axis_vec
            kps[f"petiole_{idx_p}"] = (round(float(p_coord[0] + ox), 2), round(float(p_coord[1] + oy), 2), 0.50)

    kps["width_1"] = (round(float(w1_pt[0] + ox), 2), round(float(w1_pt[1] + oy), 2), 0.90)
    kps["width_2"] = (round(float(w2_pt[0] + ox), 2), round(float(w2_pt[1] + oy), 2), 0.90)

    pet_w = estimate_petiole_width_px(kps, petiole_mask=petiole_mask, lamina_mask=mask_use)
    return LM3PoseResult(kps, petiole_width_px=pet_w, is_geometric_fallback=True)
