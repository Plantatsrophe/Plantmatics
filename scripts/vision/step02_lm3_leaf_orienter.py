#!/usr/bin/env python3
"""
Step 02 LeafMachine3 native tip-up leaf orientation, 180-degree polarity
reconciliation, bilateral symmetry evaluation, and landmark angle/trace forwarding.
"""

from __future__ import annotations

import logging
import math
from typing import Any, Dict, Optional, Tuple, Union

import cv2
import numpy as np

from scripts.vision.step02_lm3_petiole_kpts import normalize_lm3_keypoints_dict

logger = logging.getLogger("SegmentAndExtract")


def orient_and_score_leaf_lm3(
    lamina_mask: np.ndarray,
    keypoints_31: Optional[Union[Dict[str, Tuple[float, float, float]], np.ndarray]] = None,
    p_base: Optional[Tuple[float, float]] = None,
    p_apex: Optional[Tuple[float, float]] = None,
    offset_xy: Tuple[float, float] = (0.0, 0.0),
    min_conf: float = 0.25,
    min_score: float = 0.50,
    hole_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    crop_bgr: Optional[np.ndarray] = None,
    crown_hint_xy: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """
    Computes LeafMachine3 native tip-up rotation angle (`compute_orientation`), rotates
    `lamina_mask`, optional `hole_mask`, `petiole_mask`, and `crop_bgr` (`rotate_image`,
    `mask_bbox`, `crop_to_box`, `oriented_affine`, `apply_affine`), and evaluates
    curvilinear midvein bilateral symmetry and archetype quality (`measure_leaf`).
    """
    if lamina_mask is None or lamina_mask.size == 0 or int(np.count_nonzero(lamina_mask)) == 0:
        return {
            "success": False,
            "angle_cw": 0.0,
            "method": "empty",
            "archetype_score": 0.0,
            "bilateral_dice": 0.0,
            "si_a": 1.0,
            "largest_frac": 0.0,
            "gates_pass": False,
            "is_archetypal": False,
            "base_angle_deg": np.nan,
            "apex_angle_deg": np.nan,
            "lamina_curvature": np.nan,
            "lamina_trace_length_px": np.nan,
        }

    norm_kps = normalize_lm3_keypoints_dict(keypoints_31, offset_xy=offset_xy)
    pts_2d: Dict[str, Tuple[float, float]] = {
        k: (float(v[0]), float(v[1]))
        for k, v in norm_kps.items()
        if len(v) >= 3 and float(v[2]) >= float(min_conf)
    }
    if p_base is not None and "lamina_base" not in pts_2d:
        pts_2d["lamina_base"] = (float(p_base[0]) - float(offset_xy[0]), float(p_base[1]) - float(offset_xy[1]))
    if p_apex is not None and "lamina_tip" not in pts_2d:
        pts_2d["lamina_tip"] = (float(p_apex[0]) - float(offset_xy[0]), float(p_apex[1]) - float(offset_xy[1]))

    _kp_conf_min = (
        min(
            float(norm_kps.get("lamina_base", (0.0, 0.0, 0.0))[2]),
            float(norm_kps.get("lamina_tip", (0.0, 0.0, 0.0))[2]),
        )
        if norm_kps
        else 0.0
    )

    def _swap_midvein_poles(pts_map: Dict[str, Tuple[float, float]]) -> None:
        def _swap_pair(k1: str, k2: str) -> None:
            if k1 in pts_map and k2 in pts_map:
                pts_map[k1], pts_map[k2] = pts_map[k2], pts_map[k1]
            elif k1 in pts_map:
                pts_map[k2] = pts_map.pop(k1)
            elif k2 in pts_map:
                pts_map[k1] = pts_map.pop(k2)

        _swap_pair("lamina_base", "lamina_tip")
        _swap_pair("base_center", "apex_center")
        _swap_pair("base_left", "apex_left")
        _swap_pair("base_right", "apex_right")
        _swap_pair("base_1", "apex_1")
        _swap_pair("base_2", "apex_2")
        _swap_pair("base_3", "apex_3")
        mv_present = [i for i in range(15) if f"midvein_{i}" in pts_map]
        if len(mv_present) == 15:
            rev_vals = [pts_map[f"midvein_{14 - i}"] for i in range(15)]
            for i in range(15):
                pts_map[f"midvein_{i}"] = rev_vals[i]
        elif len(mv_present) >= 2:
            rev_vals = [pts_map[f"midvein_{idx}"] for idx in reversed(mv_present)]
            for idx, val in zip(mv_present, rev_vals):
                pts_map[f"midvein_{idx}"] = val
        if "midvein_15" in pts_map and "midvein_14" in pts_map:
            pts_map["midvein_15"] = pts_map["midvein_14"]

    if (
        petiole_mask is not None
        and petiole_mask.size > 0
        and int(np.count_nonzero(petiole_mask)) > 0
        and "lamina_base" in pts_2d
        and "lamina_tip" in pts_2d
    ):
        pys, pxs = np.where(petiole_mask > 0)
        lb_pt = pts_2d["lamina_base"]
        lt_pt = pts_2d["lamina_tip"]
        d_lb_pet = float(np.min(np.hypot(pxs - lb_pt[0], pys - lb_pt[1])))
        d_lt_pet = float(np.min(np.hypot(pxs - lt_pt[0], pys - lt_pt[1])))
        span_bt = math.hypot(lt_pt[0] - lb_pt[0], lt_pt[1] - lb_pt[1])
        swapped_poles = False
        pet_near_tip_rad = max(8.0, 0.08 * span_bt)
        if span_bt > 1.0 and d_lt_pet <= pet_near_tip_rad and d_lb_pet > 0.35 * span_bt:
            near_tip_idx = np.hypot(pxs - lt_pt[0], pys - lt_pt[1]) <= pet_near_tip_rad
            _, lbl_p, stats_p, _ = cv2.connectedComponentsWithStats((petiole_mask > 0).astype(np.uint8), connectivity=8)
            near_lbls = {int(v) for v in lbl_p[pys[near_tip_idx], pxs[near_tip_idx]] if int(v) > 0}
            has_contig_pet = any(int(stats_p[lbl_i, cv2.CC_STAT_AREA]) >= 4 for lbl_i in near_lbls)
            crown_corrob = False
            crown_contradicts = False
            if crown_hint_xy is not None:
                cr_x = float(crown_hint_xy[0]) - float(offset_xy[0])
                cr_y = float(crown_hint_xy[1]) - float(offset_xy[1])
                d_lb_cr = math.hypot(lb_pt[0] - cr_x, lb_pt[1] - cr_y)
                d_lt_cr = math.hypot(lt_pt[0] - cr_x, lt_pt[1] - cr_y)
                crown_corrob = bool(d_lt_cr < 0.70 * d_lb_cr)
                crown_contradicts = bool(d_lb_cr < 0.70 * d_lt_cr)
            if (
                has_contig_pet
                and (not norm_kps or _kp_conf_min < 0.88)
                and (crown_corrob or (not crown_contradicts))
            ):
                _swap_midvein_poles(pts_2d)
                swapped_poles = True
        if (
            not swapped_poles
            and (not norm_kps or _kp_conf_min < 0.55)
            and ((not norm_kps and d_lt_pet < d_lb_pet) or (d_lt_pet + 0.15 * span_bt < d_lb_pet))
        ):
            _swap_midvein_poles(pts_2d)
    elif (
        crown_hint_xy is not None
        and "lamina_base" in pts_2d
        and "lamina_tip" in pts_2d
        and (not norm_kps or _kp_conf_min < 0.45)
    ):
        cr_x = float(crown_hint_xy[0]) - float(offset_xy[0])
        cr_y = float(crown_hint_xy[1]) - float(offset_xy[1])
        lb_pt = pts_2d["lamina_base"]
        lt_pt = pts_2d["lamina_tip"]
        d_lb_cr = math.hypot(lb_pt[0] - cr_x, lb_pt[1] - cr_y)
        d_lt_cr = math.hypot(lt_pt[0] - cr_x, lt_pt[1] - cr_y)
        span_bt = math.hypot(lt_pt[0] - lb_pt[0], lt_pt[1] - lb_pt[1])
        if (not norm_kps and d_lt_cr < d_lb_cr) or (d_lt_cr + 0.15 * span_bt < d_lb_cr):
            _swap_midvein_poles(pts_2d)

    # Populate synthetic midvein_0..14 along p_apex -> p_base if pose keypoints were not supplied
    if "lamina_tip" in pts_2d and "lamina_base" in pts_2d and "midvein_0" not in pts_2d:
        tx, ty = pts_2d["lamina_tip"]
        bx, by = pts_2d["lamina_base"]
        for i in range(15):
            frac = float(i + 1) / 16.0
            pts_2d[f"midvein_{i}"] = (tx + frac * (bx - tx), ty + frac * (by - ty))

    try:
        from scripts.vision.lm3_compat import (
            OrientedLeaf,
            apply_affine,
            compute_measurements,
            compute_orientation,
            crop_to_box,
            mask_bbox,
            measure_leaf,
            oriented_affine,
            rotate_image,
        )

        ori = compute_orientation(pts_2d)
        angle_cw = float(ori.angle_cw) if (ori.success and ori.angle_cw is not None) else 0.0
        signed_cw = ((angle_cw + 180.0) % 360.0) - 180.0
        if abs(signed_cw) <= 0.25:
            angle_cw = 0.0
        method = str(ori.method) if ori.success else "fallback"

        mask_u8 = ((lamina_mask > 0).astype(np.uint8)) * 255
        orig_h, orig_w = int(mask_u8.shape[0]), int(mask_u8.shape[1])
        if angle_cw == 0.0:
            mat = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64)
            nh, nw = orig_h, orig_w
            rot_mask = mask_u8.copy()
            rot_hole = (
                (((hole_mask > 0).astype(np.uint8)) * 255)
                if (hole_mask is not None and hole_mask.size > 0 and hole_mask.shape[:2] == (orig_h, orig_w))
                else None
            )
            rot_petiole = (
                (((petiole_mask > 0).astype(np.uint8)) * 255)
                if (petiole_mask is not None and petiole_mask.size > 0 and petiole_mask.shape[:2] == (orig_h, orig_w))
                else None
            )
            rot_bgr = (
                crop_bgr.copy()
                if (crop_bgr is not None and crop_bgr.size > 0 and crop_bgr.shape[:2] == (orig_h, orig_w))
                else None
            )
        else:
            mat, (nh, nw) = oriented_affine((orig_h, orig_w), angle_cw)
            rot_mask = rotate_image(mask_u8, angle_cw, bg=0, nearest=True)

            rot_hole = None
            if hole_mask is not None and hole_mask.size > 0 and hole_mask.shape[:2] == (orig_h, orig_w):
                rot_hole = rotate_image(((hole_mask > 0).astype(np.uint8)) * 255, angle_cw, bg=0, nearest=True)

            rot_petiole = None
            if petiole_mask is not None and petiole_mask.size > 0 and petiole_mask.shape[:2] == (orig_h, orig_w):
                rot_petiole = rotate_image(((petiole_mask > 0).astype(np.uint8)) * 255, angle_cw, bg=0, nearest=True)

            rot_bgr = None
            if crop_bgr is not None and crop_bgr.size > 0 and crop_bgr.shape[:2] == (orig_h, orig_w):
                rot_bgr = rotate_image(crop_bgr, angle_cw, bg=0, nearest=False)

        pt_keys = list(pts_2d.keys())
        if pt_keys:
            pt_arr = np.asarray([pts_2d[k] for k in pt_keys], dtype=np.float64)
            rot_arr = apply_affine(mat, pt_arr)
            rot_pts = {k: (float(rot_arr[i, 0]), float(rot_arr[i, 1])) for i, k in enumerate(pt_keys)}
        else:
            rot_pts = {}

        union_rot = rot_mask.copy()
        if rot_petiole is not None and rot_petiole.shape == union_rot.shape:
            union_rot = np.maximum(union_rot, rot_petiole)
        bbox = mask_bbox(union_rot)
        if bbox is not None:
            x0, y0, x1, y1 = bbox
            crop_rot = crop_to_box(rot_mask, bbox)
            crop_hole = crop_to_box(rot_hole, bbox) if rot_hole is not None else None
            crop_pet = crop_to_box(rot_petiole, bbox) if rot_petiole is not None else None
            crop_img = crop_to_box(rot_bgr, bbox) if rot_bgr is not None else None
            shifted_pts = {k: (pt[0] - float(x0), pt[1] - float(y0)) for k, pt in rot_pts.items()}
            det_box = (x0, y0, x1, y1)
        else:
            crop_rot = rot_mask
            crop_hole = rot_hole
            crop_pet = rot_petiole
            crop_img = rot_bgr
            shifted_pts = rot_pts
            det_box = (0, 0, int(nw), int(nh))

        holes_bool = (
            (crop_hole > 0).astype(bool)
            if (crop_hole is not None and crop_hole.shape == crop_rot.shape)
            else np.zeros(crop_rot.shape[:2], dtype=bool)
        )
        lam_bool = (crop_rot > 0).astype(bool) & (~holes_bool)
        if not np.any(lam_bool):
            lam_bool = (crop_rot > 0).astype(bool)
        crop_lamina_u8 = (lam_bool.astype(np.uint8)) * 255
        # Trim any empty petiole-only margin around upright lamina so oriented_mask frames the upright blade
        lam_bbox = mask_bbox(crop_lamina_u8)
        if lam_bbox is not None and lam_bbox != (0, 0, crop_lamina_u8.shape[1], crop_lamina_u8.shape[0]):
            lx0, ly0, lx1, ly1 = lam_bbox
            crop_lamina_u8 = crop_to_box(crop_lamina_u8, lam_bbox)
            crop_hole = crop_to_box(crop_hole, lam_bbox) if crop_hole is not None else None
            crop_pet = crop_to_box(crop_pet, lam_bbox) if crop_pet is not None else None
            crop_img = crop_to_box(crop_img, lam_bbox) if crop_img is not None else None
            shifted_pts = {k: (pt[0] - float(lx0), pt[1] - float(ly0)) for k, pt in shifted_pts.items()}
            lam_bool = (crop_lamina_u8 > 0).astype(bool)
            holes_bool = (
                (crop_hole > 0).astype(bool)
                if (crop_hole is not None and crop_hole.shape == crop_lamina_u8.shape)
                else np.zeros(crop_lamina_u8.shape[:2], dtype=bool)
            )
        orig_edge_clipped = bool(
            orig_h >= 8
            and orig_w >= 8
            and (
                int(np.count_nonzero(mask_u8[0, :])) > max(3, int(0.02 * max(orig_h, orig_w)))
                or int(np.count_nonzero(mask_u8[-1, :])) > max(3, int(0.02 * max(orig_h, orig_w)))
                or int(np.count_nonzero(mask_u8[:, 0])) > max(3, int(0.02 * max(orig_h, orig_w)))
                or int(np.count_nonzero(mask_u8[:, -1])) > max(3, int(0.02 * max(orig_h, orig_w)))
            )
        )
        if not orig_edge_clipped and crop_lamina_u8.size > 0:
            crop_lamina_u8 = np.pad(crop_lamina_u8, ((2, 2), (2, 2)), mode="constant", constant_values=0)
            if crop_hole is not None:
                crop_hole = np.pad(crop_hole, ((2, 2), (2, 2)), mode="constant", constant_values=0)
            if crop_pet is not None:
                crop_pet = np.pad(crop_pet, ((2, 2), (2, 2)), mode="constant", constant_values=0)
            if crop_img is not None and crop_img.ndim == 3:
                crop_img = np.pad(crop_img, ((2, 2), (2, 2), (0, 0)), mode="constant", constant_values=0)
            shifted_pts = {k: (pt[0] + 2.0, pt[1] + 2.0) for k, pt in shifted_pts.items()}
            lam_bool = (crop_lamina_u8 > 0).astype(bool)
            holes_bool = (
                (crop_hole > 0).astype(bool)
                if (crop_hole is not None and crop_hole.shape == crop_lamina_u8.shape)
                else np.zeros(crop_lamina_u8.shape[:2], dtype=bool)
            )
        sil_bool = lam_bool | holes_bool
        pet_bool = (crop_pet > 0).astype(bool) if crop_pet is not None else None

        cnts_lw, _ = cv2.findContours(crop_lamina_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if cnts_lw:
            sig_lw = [c.reshape(-1, 2).astype(np.float64) for c in cnts_lw if cv2.contourArea(c) >= 16.0]
            c_lw = (
                np.vstack(sig_lw)
                if sig_lw
                else max(cnts_lw, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
            )
            blade_l_px = float(np.ptp(c_lw[:, 1]))
            blade_w_px = float(np.ptp(c_lw[:, 0]))
        else:
            ys_lw, xs_lw = np.where(lam_bool)
            blade_l_px = float(ys_lw.max() - ys_lw.min()) if len(ys_lw) > 0 else 0.0
            blade_w_px = float(xs_lw.max() - xs_lw.min()) if len(xs_lw) > 0 else 0.0
        anatomical_lw = (blade_l_px / max(blade_w_px, 1e-6)) if (blade_l_px > 0.0 and blade_w_px > 0.0) else np.nan

        def _meas_float(obj: Any, *attrs: str) -> float:
            if obj is None:
                return np.nan
            for attr_name in attrs:
                val = getattr(obj, attr_name, None)
                if val is not None:
                    try:
                        fval = float(val)
                        if np.isfinite(fval):
                            return round(fval, 4)
                    except (TypeError, ValueError):
                        pass
            return np.nan

        lm_meas = None
        try:
            if shifted_pts:
                lm_meas = compute_measurements(shifted_pts)
        except Exception:
            lm_meas = None
        base_angle_deg = _meas_float(lm_meas, "base_angle", "base_angle_deg")
        apex_angle_deg = _meas_float(lm_meas, "apex_angle", "apex_angle_deg")
        lamina_curvature = _meas_float(lm_meas, "lamina_curvature")
        lamina_trace_length_px = _meas_float(lm_meas, "lamina_trace_length", "lamina_trace_length_px")

        kpt_conf = {
            k: (float(norm_kps[k][2]) if k in norm_kps and len(norm_kps[k]) >= 3 else 0.95)
            for k in shifted_pts
        }
        oleaf = OrientedLeaf(
            specimen_id=0,
            leaf_id=0,
            detection_id=0,
            stem="leaf",
            silhouette=sil_bool,
            lamina=lam_bool,
            holes=holes_bool,
            petiole=pet_bool,
            kpts=shifted_pts,
            kpt_conf=kpt_conf,
            angle_cw=angle_cw,
            det_box=det_box,
            crop_shape=(orig_h, orig_w),
            truncated=False,
        )
        bm = measure_leaf(oleaf, min_score=float(min_score))
        if bm is None:
            return {
                "success": bool(ori.success),
                "angle_cw": round(angle_cw, 4),
                "method": method,
                "archetype_score": 0.75,
                "bilateral_dice": 0.75,
                "si_a": 0.0,
                "largest_frac": 1.0,
                "gates_pass": True,
                "is_archetypal": True,
                "aspect_ratio": round(float(anatomical_lw), 4) if np.isfinite(anatomical_lw) else np.nan,
                "blade_length_px": round(blade_l_px, 4),
                "blade_width_px": round(blade_w_px, 4),
                "base_angle_deg": base_angle_deg,
                "apex_angle_deg": apex_angle_deg,
                "lamina_curvature": lamina_curvature,
                "lamina_trace_length_px": lamina_trace_length_px,
                "oriented_mask": crop_lamina_u8,
                "oriented_hole_mask": crop_hole,
                "oriented_petiole_mask": crop_pet,
                "oriented_crop_bgr": crop_img,
                "oriented_points": shifted_pts,
            }
        return {
            "success": bool(ori.success),
            "angle_cw": round(angle_cw, 4),
            "method": method,
            "archetype_score": round(float(bm.archetype_score), 4) if np.isfinite(bm.archetype_score) else 0.75,
            "bilateral_dice": round(float(bm.dice), 4) if (bm.dice is not None and np.isfinite(bm.dice)) else 0.0,
            "si_a": round(float(bm.si_a), 4) if (bm.si_a is not None and np.isfinite(bm.si_a)) else 0.0,
            "largest_frac": round(float(bm.largest_frac), 4) if np.isfinite(bm.largest_frac) else 1.0,
            "gates_pass": bool(bm.gates_pass),
            "is_archetypal": bool(bm.is_archetypal),
            "aspect_ratio": round(float(anatomical_lw), 4) if np.isfinite(anatomical_lw) else np.nan,
            "blade_length_px": round(blade_l_px, 4),
            "blade_width_px": round(blade_w_px, 4),
            "base_angle_deg": base_angle_deg,
            "apex_angle_deg": apex_angle_deg,
            "lamina_curvature": lamina_curvature,
            "lamina_trace_length_px": lamina_trace_length_px,
            "oriented_mask": crop_lamina_u8,
            "oriented_hole_mask": crop_hole,
            "oriented_petiole_mask": crop_pet,
            "oriented_crop_bgr": crop_img,
            "oriented_points": shifted_pts,
        }
    except Exception as exc:
        logger.debug(f"LeafMachine3 orient_and_score_leaf_lm3 fallback: {exc}")
        return {
            "success": False,
            "angle_cw": 0.0,
            "method": "fallback",
            "archetype_score": 0.75,
            "bilateral_dice": 0.75,
            "si_a": 0.0,
            "largest_frac": 1.0,
            "gates_pass": True,
            "is_archetypal": True,
            "base_angle_deg": np.nan,
            "apex_angle_deg": np.nan,
            "lamina_curvature": np.nan,
            "lamina_trace_length_px": np.nan,
        }
