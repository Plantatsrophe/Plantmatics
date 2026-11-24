#!/usr/bin/env python3
"""
Step 02 LeafMachine3 read-only petiole/keypoint anchoring and 31-keypoint
canonical-to-native dictionary normalization utilities.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np

import scripts.vision.lm3_compat  # noqa: F401


def _extract_native_lm3_petiole_anchor(
    leaf_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None,
    min_kp_conf: float = 0.25,
    lamina_mask: Optional[np.ndarray] = None,
    offset_xy: Tuple[float, float] = (0.0, 0.0),
) -> Optional[Dict[str, Any]]:
    """
    Uses LeafMachine3's explicit Stage 2 `petiole_mask` (`class_id == 1`) and/or
    31-keypoint `mid15_pet5` (`lamina_base`, `petiole_0..4`, `petiole_tip`, `lamina_tip`)
    skeleton to locate the anatomical petiole-lamina junction `p_base_crop`, distal lamina
    apex `p_apex_crop`, native landmark metrics (`leafmachine3.core.landmark_metrics.compute_measurements`),
    and orthogonal petiole measurements (`leafmachine3.core.petiole.measure_petiole`).
    Falls back to the basal-most mask contour point when both `keypoints_31` and
    `petiole_mask` are absent.
    """
    if leaf_mask is None and lamina_mask is not None:
        leaf_mask = lamina_mask
    if leaf_mask is None or leaf_mask.size == 0 or int(np.count_nonzero(leaf_mask)) == 0:
        return None

    bx_l, by_l, bw_l, bh_l = cv2.boundingRect(leaf_mask)
    if bw_l == 0 or bh_l == 0:
        return None
    cnts_l, _ = cv2.findContours(
        leaf_mask[by_l : by_l + bh_l, bx_l : bx_l + bw_l],
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_NONE,
        offset=(bx_l, by_l),
    )
    if not cnts_l:
        return None
    leaf_pts = max(cnts_l, key=cv2.contourArea).reshape(-1, 2)
    if len(leaf_pts) == 0:
        return None

    p_base_crop: Optional[Tuple[int, int]] = None
    p_apex_crop: Optional[Tuple[int, int]] = None
    petiole_len_px: float = 0.0
    petiole_width_px: Optional[float] = None
    touches_leaf: bool = False
    used_native_lm3 = False
    lm_metrics_dict: Dict[str, Any] = {}
    pet_chain: List[Tuple[float, float]] = []

    norm_kps = normalize_lm3_keypoints_dict(keypoints_31, offset_xy=offset_xy) if keypoints_31 is not None else {}

    # 1. Try 31-keypoint pose skeleton (`lamina_base`, `petiole_0..4`, `petiole_tip`, `lamina_tip`)
    if norm_kps:
        pts_2d: Dict[str, Tuple[float, float]] = {
            k: (float(v[0]), float(v[1]))
            for k, v in norm_kps.items()
            if len(v) >= 3 and float(v[2]) >= float(min_kp_conf)
        }
        try:
            from leafmachine3.core import landmark_metrics

            lm_res = landmark_metrics.compute_measurements(pts_2d)
            if lm_res is not None:
                lam_tr = getattr(lm_res, "lamina_trace_length", getattr(lm_res, "lamina_trace_length_px", None))
                lam_tb = getattr(lm_res, "lamina_tip_base_length", getattr(lm_res, "lamina_tip_base_length_px", None))
                pet_tr = getattr(lm_res, "petiole_trace_length", getattr(lm_res, "petiole_trace_length_px", None))
                lw_px = getattr(lm_res, "leaf_width", getattr(lm_res, "leaf_width_px", None))
                ap_ang = getattr(lm_res, "apex_angle", getattr(lm_res, "apex_angle_deg", None))
                bs_ang = getattr(lm_res, "base_angle", getattr(lm_res, "base_angle_deg", None))
                lm_curv = getattr(lm_res, "lamina_curvature", None)
                lm_metrics_dict = {
                    "lamina_trace_length_px": lam_tr,
                    "lamina_tip_base_length_px": lam_tb,
                    "petiole_trace_length_px": pet_tr,
                    "leaf_width_px": lw_px,
                    "apex_angle_deg": ap_ang,
                    "base_angle_deg": bs_ang,
                    "lamina_curvature": lm_curv,
                }
                if pet_tr is not None and np.isfinite(pet_tr):
                    petiole_len_px = max(petiole_len_px, float(pet_tr))
        except Exception:
            pass

        lb_kp = (
            norm_kps.get("lamina_base")
            or norm_kps.get("base_center")
            or norm_kps.get("base_2")
            or norm_kps.get("midvein_14")
        )
        if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= 0.25:
            dists_lb = np.hypot(leaf_pts[:, 0] - float(lb_kp[0]), leaf_pts[:, 1] - float(lb_kp[1]))
            min_lb_idx = int(np.argmin(dists_lb))
            if float(dists_lb[min_lb_idx]) <= max(36.0, 0.35 * float(max(bw_l, bh_l))):
                p_base_crop = (int(leaf_pts[min_lb_idx, 0]), int(leaf_pts[min_lb_idx, 1]))
            else:
                p_base_crop = (int(round(float(lb_kp[0]))), int(round(float(lb_kp[1]))))
            used_native_lm3 = True

        lt_kp = (
            norm_kps.get("lamina_tip")
            or norm_kps.get("apex_center")
            or norm_kps.get("apex_2")
            or norm_kps.get("midvein_0")
        )
        if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= 0.25:
            dists_lt = np.hypot(leaf_pts[:, 0] - float(lt_kp[0]), leaf_pts[:, 1] - float(lt_kp[1]))
            min_lt_idx = int(np.argmin(dists_lt))
            if float(dists_lt[min_lt_idx]) <= max(36.0, 0.35 * float(max(bw_l, bh_l))):
                p_apex_crop = (int(leaf_pts[min_lt_idx, 0]), int(leaf_pts[min_lt_idx, 1]))
            else:
                p_apex_crop = (int(round(float(lt_kp[0]))), int(round(float(lt_kp[1]))))

        if p_base_crop is not None:
            pet_chain.append((float(p_base_crop[0]), float(p_base_crop[1])))
        elif lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= min_kp_conf:
            pet_chain.append((float(lb_kp[0]), float(lb_kp[1])))
        for native_pk in ("petiole_0", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_tip"):
            pkp = norm_kps.get(native_pk)
            if pkp is not None and len(pkp) >= 3 and float(pkp[2]) >= min_kp_conf:
                pt_xy = (float(pkp[0]), float(pkp[1]))
                if not pet_chain or math.hypot(pt_xy[0] - pet_chain[-1][0], pt_xy[1] - pet_chain[-1][1]) > 1e-3:
                    pet_chain.append(pt_xy)
        if len(pet_chain) >= 2 and petiole_len_px <= 0.0:
            arc_len = 0.0
            for i_c in range(len(pet_chain) - 1):
                arc_len += math.hypot(
                    pet_chain[i_c + 1][0] - pet_chain[i_c][0],
                    pet_chain[i_c + 1][1] - pet_chain[i_c][1],
                )
            petiole_len_px = max(petiole_len_px, float(arc_len))

    # 2. Delegate orthogonal petiole measurement directly to LeafMachine3 `measure_petiole`
    if (
        petiole_mask is not None
        and petiole_mask.size > 0
        and petiole_mask.shape[:2] == leaf_mask.shape[:2]
        and int(np.count_nonzero(petiole_mask)) > 0
    ):
        pet_u8 = (petiole_mask > 0).astype(np.uint8)
        lam_u8 = (leaf_mask > 0).astype(np.uint8)
        if p_base_crop is None:
            pys, pxs = np.where(pet_u8 > 0)
            pet_cx, pet_cy = float(np.mean(pxs)), float(np.mean(pys))
            kern_j = np.ones((5, 5), dtype=np.uint8)
            jys_c, jxs_c = np.where((cv2.dilate(pet_u8, kern_j, iterations=1) > 0) & (lam_u8 > 0))
            if len(jxs_c) > 0:
                j_cx, j_cy = float(np.median(jxs_c)), float(np.median(jys_c))
                j_idx = int(np.argmin(np.hypot(leaf_pts[:, 0] - j_cx, leaf_pts[:, 1] - j_cy)))
            else:
                lys_c, lxs_c = np.where(lam_u8 > 0)
                lam_cx = float(np.mean(lxs_c)) if len(lxs_c) > 0 else float(np.mean(leaf_pts[:, 0]))
                lam_cy = float(np.mean(lys_c)) if len(lys_c) > 0 else float(np.mean(leaf_pts[:, 1]))
                rx_p, ry_p = lam_cx - pet_cx, lam_cy - pet_cy
                rn_p = math.hypot(rx_p, ry_p)
                if rn_p > 1.0:
                    ux_p, uy_p = rx_p / rn_p, ry_p / rn_p
                    dx_lp = leaf_pts[:, 0] - pet_cx
                    dy_lp = leaf_pts[:, 1] - pet_cy
                    along_p = dx_lp * ux_p + dy_lp * uy_p
                    perp_p = np.abs(dx_lp * (-uy_p) + dy_lp * ux_p)
                    corridor = (perp_p <= max(8.0, 0.06 * float(max(bw_l, bh_l)))) & (along_p <= rn_p)
                    if np.any(corridor):
                        c_idxs = np.where(corridor)[0]
                        min_al = float(np.min(along_p[c_idxs]))
                        tie_pts = leaf_pts[c_idxs[along_p[c_idxs] <= (min_al + 1.5)]]
                        t_cx, t_cy = float(np.median(tie_pts[:, 0])), float(np.median(tie_pts[:, 1]))
                        j_idx = int(np.argmin(np.hypot(leaf_pts[:, 0] - t_cx, leaf_pts[:, 1] - t_cy)))
                    else:
                        j_idx = int(np.argmin(np.hypot(dx_lp, dy_lp)))
                else:
                    j_idx = int(np.argmin(np.hypot(leaf_pts[:, 0] - pet_cx, leaf_pts[:, 1] - pet_cy)))
            p_base_crop = (int(leaf_pts[j_idx, 0]), int(leaf_pts[j_idx, 1]))
            used_native_lm3 = True
            if len(pet_chain) < 2:
                far_i = int(np.argmax(np.hypot(pxs - p_base_crop[0], pys - p_base_crop[1])))
                pet_chain = [(float(p_base_crop[0]), float(p_base_crop[1])), (float(pxs[far_i]), float(pys[far_i]))]
                petiole_len_px = max(petiole_len_px, math.hypot(pet_chain[1][0] - pet_chain[0][0], pet_chain[1][1] - pet_chain[0][1]))
        if len(pet_chain) >= 2:
            try:
                from leafmachine3.core import petiole as petiole_metrics

                pet_res = petiole_metrics.measure_petiole(pet_u8, lam_u8, pet_chain)
                if pet_res is not None:
                    touches_leaf = bool(pet_res.touches_leaf)
                    if pet_res.width_px is not None and np.isfinite(pet_res.width_px) and pet_res.width_px > 0:
                        petiole_width_px = round(float(pet_res.width_px), 3)
                    if pet_res.length_px is not None and np.isfinite(pet_res.length_px) and pet_res.length_px > 0:
                        petiole_len_px = float(pet_res.length_px)
            except Exception:
                pass

    def _select_midline_apex(base_xy: Tuple[int, int]) -> Tuple[int, int]:
        ys_m, xs_m = np.where(leaf_mask > 0)
        if len(xs_m) > 0:
            cx_m, cy_m = float(np.mean(xs_m)), float(np.mean(ys_m))
        else:
            cx_m, cy_m = float(np.mean(leaf_pts[:, 0])), float(np.mean(leaf_pts[:, 1]))
        vx_a = cx_m - float(base_xy[0])
        vy_a = cy_m - float(base_xy[1])
        vn_a = math.hypot(vx_a, vy_a)
        if vn_a > 1.0:
            ux_a, uy_a = vx_a / vn_a, vy_a / vn_a
            proj_a = (leaf_pts[:, 0] - float(base_xy[0])) * ux_a + (leaf_pts[:, 1] - float(base_xy[1])) * uy_a
            max_p = float(np.max(proj_a))
            cap_pts = leaf_pts[proj_a >= (max_p - 1.5)]
            if len(cap_pts) > 0:
                cap_cx = float(np.median(cap_pts[:, 0]))
                cap_cy = float(np.median(cap_pts[:, 1]))
                a_i = int(np.argmin(np.hypot(cap_pts[:, 0] - cap_cx, cap_pts[:, 1] - cap_cy)))
                return (int(cap_pts[a_i, 0]), int(cap_pts[a_i, 1]))
        d_fb = np.hypot(leaf_pts[:, 0] - float(base_xy[0]), leaf_pts[:, 1] - float(base_xy[1]))
        a_fb = int(np.argmax(d_fb))
        return (int(leaf_pts[a_fb, 0]), int(leaf_pts[a_fb, 1]))

    if p_base_crop is not None:
        if p_apex_crop is None:
            p_apex_crop = _select_midline_apex(p_base_crop)
        out_anchor: Dict[str, Any] = {
            "is_trimmed": False,
            "native_lm3_petiole": bool(used_native_lm3),
            "trim_distance_px": round(float(petiole_len_px), 2),
            "petiole_length_px": round(float(petiole_len_px), 2),
            "geodesic_petiole_skeleton_px": round(float(petiole_len_px), 2),
            "petiole_width_px": petiole_width_px,
            "touches_leaf": bool(touches_leaf),
            "p_base_crop": p_base_crop,
            "p_apex_crop": p_apex_crop,
            "p_base": p_base_crop,
            "p_apex": p_apex_crop,
            "landmark_measurements": lm_metrics_dict,
        }
        out_anchor.update({k: v for k, v in lm_metrics_dict.items() if v is not None})
        return out_anchor

    # 3. Fallback to basal-most midline contour vertex when both keypoints and petiole_mask are absent
    y_max_lp = float(np.max(leaf_pts[:, 1]))
    bot_cap = leaf_pts[leaf_pts[:, 1] >= (y_max_lp - 1.5)]
    bot_cx = float(np.median(bot_cap[:, 0])) if len(bot_cap) > 0 else float(np.mean(leaf_pts[:, 0]))
    b_sub_i = int(np.argmin(np.abs(bot_cap[:, 0] - bot_cx))) if len(bot_cap) > 0 else int(np.argmax(leaf_pts[:, 1]))
    p_base_crop = (int(bot_cap[b_sub_i, 0]), int(bot_cap[b_sub_i, 1])) if len(bot_cap) > 0 else (int(leaf_pts[b_sub_i, 0]), int(leaf_pts[b_sub_i, 1]))
    p_apex_crop = _select_midline_apex(p_base_crop)
    out_fallback: Dict[str, Any] = {
        "is_trimmed": False,
        "native_lm3_petiole": False,
        "morphological_fallback": True,
        "trim_distance_px": 0.0,
        "petiole_length_px": 0.0,
        "geodesic_petiole_skeleton_px": 0.0,
        "petiole_width_px": petiole_width_px,
        "touches_leaf": False,
        "p_base_crop": p_base_crop,
        "p_apex_crop": p_apex_crop,
        "p_base": p_base_crop,
        "p_apex": p_apex_crop,
        "landmark_measurements": lm_metrics_dict,
    }
    out_fallback.update({k: v for k, v in lm_metrics_dict.items() if v is not None})
    return out_fallback


_LM3_CANONICAL_TO_NATIVE_KPT_MAP: Dict[str, str] = {
    "lamina_tip": "lamina_tip",
    "apex_1": "apex_left",
    "apex_2": "apex_center",
    "apex_3": "apex_right",
    **{f"midvein_{i}": f"midvein_{i - 1}" for i in range(1, 16)},
    "base_1": "base_left",
    "base_2": "base_center",
    "base_3": "base_right",
    "lamina_base": "lamina_base",
    **{f"petiole_{i}": f"petiole_{i - 1}" for i in range(1, 6)},
    "petiole_6": "petiole_tip",
    "width_1": "width_left",
    "width_2": "width_right",
}


def normalize_lm3_keypoints_dict(
    keypoints_31: Optional[Union[Dict[str, Tuple[float, float, float]], np.ndarray]],
    offset_xy: Tuple[float, float] = (0.0, 0.0),
) -> Dict[str, Tuple[float, float, float]]:
    """
    Normalizes a 31-keypoint dictionary or `(31, 3)` array so both `LM3_KEYPOINT_NAMES_31`
    (`midvein_1..15`, `apex_1..3`, `base_1..3`, `petiole_1..6`, `width_1..2`) and LeafMachine3
    native `KPT_NAMES` (`midvein_0..14`, `apex_left/center/right`, `base_left/center/right`,
    `petiole_0..4`, `petiole_tip`, `width_left/right`) are populated. Applies optional
    `(ox, oy)` subtraction when converting sheet-space keypoints back to crop coordinates.
    """
    ox, oy = float(offset_xy[0]), float(offset_xy[1])
    out: Dict[str, Tuple[float, float, float]] = {}
    if keypoints_31 is None:
        return out
    if isinstance(keypoints_31, dict):
        for k, val in keypoints_31.items():
            if not isinstance(val, (tuple, list, np.ndarray)) or len(val) < 2:
                continue
            vx = float(val[0]) - ox
            vy = float(val[1]) - oy
            vc = float(val[2]) if len(val) >= 3 else 1.0
            out[str(k)] = (vx, vy, vc)
        src_kpts = dict(out)
        # 1. Resolve 1-based (`midvein_1..15`) vs 0-based (`midvein_0..14`) midvein indexing without collision
        if "midvein_0" not in src_kpts and "midvein_1" in src_kpts:
            for i in range(1, 16):
                ck = f"midvein_{i}"
                if ck in src_kpts:
                    out[f"midvein_{i - 1}"] = src_kpts[ck]
        elif "midvein_0" in src_kpts and "midvein_15" not in src_kpts:
            for i in range(1, 16):
                ck = f"midvein_{i}"
                nk = f"midvein_{i - 1}"
                if nk in src_kpts and ck not in src_kpts:
                    out[ck] = src_kpts[nk]
        # 2. Resolve 1-based (`petiole_1..6`) vs 0-based (`petiole_0..4`, `petiole_tip`) petiole indexing
        if "petiole_0" not in src_kpts and "petiole_1" in src_kpts:
            for i in range(1, 6):
                ck = f"petiole_{i}"
                if ck in src_kpts:
                    out[f"petiole_{i - 1}"] = src_kpts[ck]
            if "petiole_6" in src_kpts:
                out["petiole_tip"] = src_kpts["petiole_6"]
            elif "petiole_5" in src_kpts and "petiole_tip" not in src_kpts:
                out["petiole_tip"] = src_kpts["petiole_5"]
        elif "petiole_0" in src_kpts and "petiole_6" not in src_kpts:
            for i in range(1, 6):
                ck = f"petiole_{i}"
                nk = f"petiole_{i - 1}"
                if nk in src_kpts and ck not in src_kpts:
                    out[ck] = src_kpts[nk]
            if "petiole_tip" in src_kpts and "petiole_6" not in src_kpts:
                out["petiole_6"] = src_kpts["petiole_tip"]
        # 3. Map non-overlapping canonical <-> native aliases
        _NON_INDEXED_ALIASES = (
            ("lamina_tip", "apex_center"),
            ("apex_2", "apex_center"),
            ("apex_1", "apex_left"),
            ("apex_3", "apex_right"),
            ("lamina_base", "base_center"),
            ("base_2", "base_center"),
            ("base_1", "base_left"),
            ("base_3", "base_right"),
            ("petiole_6", "petiole_tip"),
            ("width_1", "width_left"),
            ("width_2", "width_right"),
        )
        for canon_k, native_k in _NON_INDEXED_ALIASES:
            if canon_k in src_kpts and native_k not in out:
                out[native_k] = src_kpts[canon_k]
            if native_k in src_kpts and canon_k not in out:
                out[canon_k] = src_kpts[native_k]
        if "lamina_tip" in out and "apex_2" not in out:
            out["apex_2"] = out["lamina_tip"]
        if "lamina_base" in out and "base_2" not in out:
            out["base_2"] = out["lamina_base"]
        return out

    arr = np.asarray(keypoints_31, dtype=np.float64)
    if arr.ndim == 2 and arr.shape[0] >= 2 and arr.shape[1] >= 2:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        # Support both LM3_KEYPOINT_NAMES_31 order (lamina_base at index 22) and
        # compact array order (0: lamina_tip, 1: lamina_base, 2..7: petiole, 8..30: midvein)
        c22 = float(arr[22, 2]) if (arr.shape[0] > 22 and arr.shape[1] >= 3) else 0.0
        c4 = float(arr[4, 2]) if (arr.shape[0] > 4 and arr.shape[1] >= 3) else 0.0
        if arr.shape[0] >= 23 and c22 > 0.0 and c4 > 0.0:
            for idx, canon_k in enumerate(_lm3_eng.LM3_KEYPOINT_NAMES_31):
                if idx >= arr.shape[0]:
                    break
                vx = float(arr[idx, 0]) - ox
                vy = float(arr[idx, 1]) - oy
                vc = float(arr[idx, 2]) if arr.shape[1] >= 3 else 1.0
                out[canon_k] = (vx, vy, vc)
                native_k = _LM3_CANONICAL_TO_NATIVE_KPT_MAP.get(canon_k)
                if native_k:
                    out[native_k] = (vx, vy, vc)
        else:
            vx0, vy0 = float(arr[0, 0]) - ox, float(arr[0, 1]) - oy
            vc0 = float(arr[0, 2]) if arr.shape[1] >= 3 else 1.0
            out["lamina_tip"] = (vx0, vy0, vc0)
            vx1, vy1 = float(arr[1, 0]) - ox, float(arr[1, 1]) - oy
            vc1 = float(arr[1, 2]) if arr.shape[1] >= 3 else 1.0
            out["lamina_base"] = (vx1, vy1, vc1)
            for p_i in range(min(6, max(0, arr.shape[0] - 2))):
                row_idx = 2 + p_i
                vxp, vyp = float(arr[row_idx, 0]) - ox, float(arr[row_idx, 1]) - oy
                vcp = float(arr[row_idx, 2]) if arr.shape[1] >= 3 else 1.0
                canon_pk = f"petiole_{p_i + 1}"
                out[canon_pk] = (vxp, vyp, vcp)
                native_pk = _LM3_CANONICAL_TO_NATIVE_KPT_MAP.get(canon_pk)
                if native_pk:
                    out[native_pk] = (vxp, vyp, vcp)
            mv_rows = arr[8:] if arr.shape[0] > 8 else np.empty((0, arr.shape[1]))
            if len(mv_rows) > 0:
                idxs = np.linspace(0, len(mv_rows) - 1, min(15, len(mv_rows))).round().astype(int)
                for m_i, r_idx in enumerate(idxs):
                    vxm, vym = float(mv_rows[r_idx, 0]) - ox, float(mv_rows[r_idx, 1]) - oy
                    vcm = float(mv_rows[r_idx, 2]) if arr.shape[1] >= 3 else 1.0
                    # Ensure tip -> base ordering along midvein_0..14
                    out[f"midvein_{m_i}"] = (vxm, vym, vcm)
            if "midvein_0" in out and "midvein_14" in out:
                d0_tip = math.hypot(out["midvein_0"][0] - vx0, out["midvein_0"][1] - vy0)
                d14_tip = math.hypot(out["midvein_14"][0] - vx0, out["midvein_14"][1] - vy0)
                if d0_tip > d14_tip:
                    rev_mv = [out[f"midvein_{i}"] for i in range(14, -1, -1) if f"midvein_{i}" in out]
                    for i, val in enumerate(rev_mv):
                        out[f"midvein_{i}"] = val
                out["midvein_15"] = out["midvein_14"]
    return out
