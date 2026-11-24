#!/usr/bin/env python3
"""
Step 02 LeafMachine3 morphological damage/fold QC gatekeeper and UNet++ v9
whole-specimen paper-removal segmenter (`LM3SpecimenSegmenter`).
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np

from scripts.vision.step02_lm3_onnx_runtime import _prepare_working_image
from scripts.vision.step02_lm3_petiole_kpts import normalize_lm3_keypoints_dict

logger = logging.getLogger("SegmentAndExtract")


def _count_half_width_lobes(half_widths: np.ndarray, max_hw: float) -> int:
    """Counts prominent lateral lobe maxima along a smoothed 1D margin half-width profile."""
    if half_widths.size < 9 or max_hw <= 3.0:
        return 0
    k_win = max(3, min(9, (int(half_widths.size // 12) * 2) + 1))
    kernel = np.ones(k_win, dtype=np.float64) / float(k_win)
    sm = np.convolve(half_widths, kernel, mode="same")
    prom_min = max(1.8, 0.10 * float(max_hw))
    peaks: List[int] = []
    for i in range(2, len(sm) - 2):
        if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1] and (sm[i] > sm[i - 2] or sm[i] > sm[i + 2]):
            left_trough = float(np.min(sm[max(0, i - max(4, len(sm) // 6)) : i + 1]))
            right_trough = float(np.min(sm[i : min(len(sm), i + max(5, len(sm) // 6))]))
            if (float(sm[i]) - max(left_trough, right_trough)) >= prom_min or (
                (float(sm[i]) - min(left_trough, right_trough)) >= 1.35 * prom_min
            ):
                if not peaks or (i - peaks[-1]) >= max(3, len(sm) // 14):
                    peaks.append(i)
    return len(peaks)


def evaluate_leaf_damage_and_fold_qc(
    oriented_mask_u8: Optional[np.ndarray],
    solidity: float = 1.0,
    bilateral_dice: float = 1.0,
    si_a: float = 0.0,
    lm3_gates_pass: bool = True,
    n_lamina_components: int = 1,
    detector_class_id: int = 0,
    segment_partial_leaves: bool = False,
    is_crop_clipped: bool = False,
    hole_area_frac: float = 0.0,
    enclosed_void_frac: float = 0.0,
    oriented_points: Optional[Dict[str, Tuple[float, float]]] = None,
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None,
    is_geometric_fallback_pose: bool = False,
) -> Dict[str, Any]:
    """
    Evaluates morphological and anatomical leaf QC gates on an upright LeafMachine3 lamina mask:
      1. Stage 1 `Leaf_PARTIAL` (`detector_class_id == 1` when `segment_partial_leaves == False`)
      2. Crop/sheet border lamina truncation (`is_crop_clipped`)
      3. Severe interior or marginal insect/herbivory perforation (`hole_area_frac > 0.035`)
      4. 31-keypoint midvein pose fold/bend (`FOLDED_MIDVEIN_TRACE`, `APICAL_TIP_FOLD`, `LOW_LANDMARK_CONFIDENCE`)
      5. Longitudinal margin fold crease (`LONGITUDINAL_MARGIN_CREASE`, ruler-straight folded edge vs curved free margin)
      6. Verified `bilateral_lobe_symmetry` closing the `solidity < 0.72` unilateral-tear loophole
      7. Unopened broad inter-lobe overlap loops (`ENCLOSED_INTERLOBE_VOID`, `enclosed_void_frac > 0.020`)
         and multi-island pinnae (`MULTI_COMPONENT_PINNAE`, `n_lamina_components > 1`).
    """
    reasons: List[str] = []
    bilateral_lobe_sym = False
    min_margin_rms_norm = np.nan
    crease_rms_ratio = 1.0

    if int(detector_class_id) == 1 and not bool(segment_partial_leaves):
        reasons.append("LEAF_PARTIAL_FRAGMENT")
    if bool(is_crop_clipped):
        reasons.append("CROP_BORDER_CLIPPED")
    if np.isfinite(hole_area_frac) and float(hole_area_frac) > 0.035:
        reasons.append("SEVERE_HERBIVORY_HOLES")

    # Inspect upright blade margin profiles for bilateral lobing and longitudinal fold creases
    if oriented_mask_u8 is not None and oriented_mask_u8.size > 0 and int(np.count_nonzero(oriented_mask_u8)) >= 32:
        fg = oriented_mask_u8 > 0
        ys_all, xs_all = np.where(fg)
        if ys_all.size >= 16:
            y_min, y_max = int(ys_all.min()), int(ys_all.max())
            x_min, x_max = int(xs_all.min()), int(xs_all.max())
            blade_l = float(max(y_max - y_min, 1))
            blade_w = float(max(x_max - x_min, 1))

            pts_map = oriented_points or {}
            mid_xs = [
                float(pts_map[k][0])
                for k in ("lamina_tip", "lamina_base", *[f"midvein_{i}" for i in range(15)])
                if k in pts_map and pts_map[k] is not None
            ]
            x_axis = float(np.median(mid_xs)) if mid_xs else float(np.mean(xs_all))

            # Scanlines over [15%, 85%] of upright blade length
            y0_s = int(round(y_min + 0.15 * blade_l))
            y1_s = int(round(y_min + 0.85 * blade_l))
            if (y1_s - y0_s) >= 12 and blade_w >= 10.0:
                left_x_list: List[float] = []
                right_x_list: List[float] = []
                scan_ys: List[float] = []
                for y_cur in range(y0_s, y1_s + 1):
                    row_xs = np.where(fg[y_cur, :])[0]
                    if row_xs.size >= 2:
                        scan_ys.append(float(y_cur))
                        left_x_list.append(float(row_xs[0]))
                        right_x_list.append(float(row_xs[-1]))

                if len(scan_ys) >= 12:
                    sy_arr = np.asarray(scan_ys, dtype=np.float64)
                    lx_arr = np.asarray(left_x_list, dtype=np.float64)
                    rx_arr = np.asarray(right_x_list, dtype=np.float64)

                    # Linear chord-deviation RMS on left vs right margins ("taco" longitudinal fold detector)
                    p_l = np.polyfit(sy_arr, lx_arr, 1)
                    p_r = np.polyfit(sy_arr, rx_arr, 1)
                    rms_l = float(np.sqrt(np.mean((lx_arr - np.polyval(p_l, sy_arr)) ** 2))) / blade_w
                    rms_r = float(np.sqrt(np.mean((rx_arr - np.polyval(p_r, sy_arr)) ** 2))) / blade_w
                    min_margin_rms_norm = min(rms_l, rms_r)
                    max_margin_rms_norm = max(rms_l, rms_r)
                    crease_rms_ratio = max_margin_rms_norm / max(min_margin_rms_norm, 1e-4)
                    if (
                        blade_l >= 28.0
                        and min_margin_rms_norm < 0.016
                        and max_margin_rms_norm >= 0.045
                        and crease_rms_ratio >= 3.2
                    ):
                        reasons.append("LONGITUDINAL_MARGIN_CREASE")

                    # Bilateral multi-lobe verification for dissected blades (`solidity < 0.72`)
                    hw_left = np.maximum(0.0, x_axis - lx_arr)
                    hw_right = np.maximum(0.0, rx_arr - x_axis)
                    max_hw_l = float(np.max(hw_left)) if hw_left.size else 0.0
                    max_hw_r = float(np.max(hw_right)) if hw_right.size else 0.0
                    span_balance = min(max_hw_l, max_hw_r) / max(max_hw_l, max_hw_r, 1e-6)
                    lobes_l = _count_half_width_lobes(hw_left, max_hw_l)
                    lobes_r = _count_half_width_lobes(hw_right, max_hw_r)
                    if (
                        lobes_l >= 2
                        and lobes_r >= 2
                        and span_balance >= 0.58
                        and abs(float(si_a)) < 0.28
                        and float(bilateral_dice) >= 0.55
                    ):
                        bilateral_lobe_sym = True
                    elif (
                        int(n_lamina_components) > 1
                        and span_balance >= 0.55
                        and abs(float(si_a)) < 0.28
                    ):
                        bilateral_lobe_sym = True

            # Check 31-keypoint midvein trace tortuosity & apical dog-ear fold
            if keypoints_31 and isinstance(keypoints_31, dict):
                norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
                lt_kp = norm_kps.get("lamina_tip")
                lb_kp = norm_kps.get("lamina_base")
                mv_confs = [
                    float(norm_kps[f"midvein_{i}"][2])
                    for i in range(15)
                    if f"midvein_{i}" in norm_kps and len(norm_kps[f"midvein_{i}"]) >= 3
                ]
                if not is_geometric_fallback_pose:
                    tip_c = float(lt_kp[2]) if (lt_kp is not None and len(lt_kp) >= 3) else 1.0
                    base_c = float(lb_kp[2]) if (lb_kp is not None and len(lb_kp) >= 3) else 1.0
                    med_mv_c = float(np.median(mv_confs)) if mv_confs else 1.0
                    if tip_c < 0.30 or base_c < 0.30 or med_mv_c < 0.30:
                        reasons.append("LOW_LANDMARK_CONFIDENCE")

                # Midvein polyline trace vs straight tip-to-base chord
                chain_pts: List[Tuple[float, float]] = []
                if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= 0.25:
                    chain_pts.append((float(lt_kp[0]), float(lt_kp[1])))
                for i in range(15):
                    mk = f"midvein_{i}"
                    if mk in norm_kps and len(norm_kps[mk]) >= 3 and float(norm_kps[mk][2]) >= 0.25:
                        pt_m = (float(norm_kps[mk][0]), float(norm_kps[mk][1]))
                        if not chain_pts or math.hypot(pt_m[0] - chain_pts[-1][0], pt_m[1] - chain_pts[-1][1]) > 1e-3:
                            chain_pts.append(pt_m)
                if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= 0.25:
                    pt_b = (float(lb_kp[0]), float(lb_kp[1]))
                    if not chain_pts or math.hypot(pt_b[0] - chain_pts[-1][0], pt_b[1] - chain_pts[-1][1]) > 1e-3:
                        chain_pts.append(pt_b)
                if len(chain_pts) >= 4:
                    chord_px = math.hypot(chain_pts[-1][0] - chain_pts[0][0], chain_pts[-1][1] - chain_pts[0][1])
                    if chord_px >= 20.0:
                        trace_px = sum(
                            math.hypot(chain_pts[j + 1][0] - chain_pts[j][0], chain_pts[j + 1][1] - chain_pts[j][1])
                            for j in range(len(chain_pts) - 1)
                        )
                        ux_c = (chain_pts[-1][0] - chain_pts[0][0]) / chord_px
                        uy_c = (chain_pts[-1][1] - chain_pts[0][1]) / chord_px
                        max_lat = max(
                            abs((px_i - chain_pts[0][0]) * (-uy_c) + (py_i - chain_pts[0][1]) * ux_c)
                            for px_i, py_i in chain_pts[1:-1]
                        )
                        if (trace_px / chord_px) > 1.22 or (max_lat / chord_px) > 0.22:
                            reasons.append("FOLDED_MIDVEIN_TRACE")

                # Apical dog-ear fold: lamina_tip landmark displaced far below top upright blade margin
                if (
                    not is_geometric_fallback_pose
                    and oriented_points
                    and "lamina_tip" in oriented_points
                    and lt_kp is not None
                    and len(lt_kp) >= 3
                    and float(lt_kp[2]) >= 0.35
                    and blade_l >= 28.0
                ):
                    oy_tip = float(oriented_points["lamina_tip"][1])
                    if (oy_tip - float(y_min)) > 0.14 * blade_l:
                        # Exempt genuinely symmetric obcordate/emarginate apical notches
                        if float(bilateral_dice) < 0.84 or abs(float(si_a)) >= 0.14:
                            reasons.append("APICAL_TIP_FOLD")

    # Symmetry & clump/tear evaluation (closing the `solidity < 0.72` unilateral-damage loophole)
    is_verified_dissected = (float(solidity) < 0.72 and bilateral_lobe_sym) or (int(n_lamina_components) > 1 and bilateral_lobe_sym)
    min_bd = 0.58 if is_verified_dissected else 0.70
    max_sia = 0.28 if (float(solidity) < 0.72 or is_verified_dissected) else 0.32
    if float(bilateral_dice) < min_bd or abs(float(si_a)) >= max_sia or (not is_verified_dissected and not bool(lm3_gates_pass)):
        if float(solidity) < 0.72 and not bilateral_lobe_sym:
            reasons.append("UNILATERAL_ASYMMETRY_OR_TEAR")
        else:
            reasons.append("UNRESOLVED_CLUMP")

    if np.isfinite(enclosed_void_frac) and float(enclosed_void_frac) > 0.020:
        reasons.append("ENCLOSED_INTERLOBE_VOID")
    if int(n_lamina_components) > 1:
        reasons.append("MULTI_COMPONENT_PINNAE")

    dedup_reasons: List[str] = []
    for r_code in reasons:
        if r_code not in dedup_reasons:
            dedup_reasons.append(r_code)

    efa_qc_pass = len(dedup_reasons) == 0
    return {
        "bilateral_lobe_symmetry": bool(bilateral_lobe_sym),
        "efa_qc_pass": bool(efa_qc_pass),
        "qc_flag_reason": "PASS" if efa_qc_pass else ";".join(dedup_reasons),
        "min_margin_rms_norm": round(float(min_margin_rms_norm), 5) if np.isfinite(min_margin_rms_norm) else np.nan,
        "crease_rms_ratio": round(float(crease_rms_ratio), 3),
    }


class LM3SpecimenSegmenter:
    """
    LeafMachine3 Stage 2b ONNX Whole-Specimen Segmenter
    (`lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain`).
    Delegates directly to `leafmachine3.inference.specimen_segmenter` (`prep_letterbox_imagenet`,
    `post_letterbox_imagenet`) and `leafmachine3.core.paper_removal.run_paperclean`
    (`METHOD_LLR_HYST_KEEPBLOBS` in CIE L*a*b* space) to produce the cleaned plant-vs-paper
    binary mask (`specimen_paper_cleaned_mask`).
    """

    MODEL_NAME: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.50,
        img_size: int = 1024,
        max_working_dim: int = 3200,
        paperclean: bool = True,
        providers: Optional[Sequence[str]] = None,
    ) -> None:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        self.device = device if (device == "cuda" and _lm3_eng.torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.img_size = int(img_size)
        self.max_working_dim = int(max_working_dim)
        self.paperclean = bool(paperclean)
        self.weights_path, self.model_name, self.providers, self.session = _lm3_eng._load_onnx_session(
            weights_path=weights_path,
            default_model_name=self.MODEL_NAME,
            device=self.device,
            providers=providers,
        )
        self.model = self.session

    def predict_mask(
        self,
        image_bgr: np.ndarray,
        max_working_dim: Optional[int] = None,
        return_sheet_size: bool = False,
    ) -> Optional[np.ndarray]:
        """
        Runs UNet++ v9 whole-specimen segmentation + LLR hysteresis paper cleanup on `image_bgr`.
        Returns a uint8 binary mask (`0` or `255`) at working resolution (or full-sheet resolution
        when `return_sheet_size=True`). Returns `None` when ONNX weights are unavailable.
        """
        import scripts.vision.lm3_inference_engines as _lm3_eng

        active_sess = self.session if self.session is not None else self.model
        if active_sess is None or image_bgr is None or image_bgr.size == 0:
            return None
        try:
            from leafmachine3.core.paper_removal import METHOD_LLR_HYST_KEEPBLOBS, run_paperclean
            from leafmachine3.inference.specimen_segmenter import (
                _sigmoid,
                post_letterbox_imagenet,
                prep_letterbox_imagenet,
            )

            sh, sw = image_bgr.shape[:2]
            w_dim = int(max_working_dim) if max_working_dim is not None else self.max_working_dim
            work_img, _ = _prepare_working_image(image_bgr, max_working_dim=w_dim)
            wh, ww = work_img.shape[:2]
            inp_tensor, meta = prep_letterbox_imagenet(work_img, self.img_size)
            input_name = active_sess.get_inputs()[0].name
            output_name = active_sess.get_outputs()[0].name
            with _lm3_eng._ONNX_GPU_LOCK:
                raw_out = active_sess.run([output_name], {input_name: inp_tensor})[0]
            prob_sq = _sigmoid(raw_out[0, 0])
            raw_mask = post_letterbox_imagenet(prob_sq, meta, wh, ww, self.conf_thres)
            if self.paperclean:
                cleaned_mask, _, _ = run_paperclean(work_img, raw_mask, METHOD_LLR_HYST_KEEPBLOBS)
            else:
                cleaned_mask = raw_mask
            mask_u8 = (cleaned_mask > 0).astype(np.uint8) * 255
            if return_sheet_size and (wh != sh or ww != sw):
                mask_u8 = cv2.resize(mask_u8, (sw, sh), interpolation=cv2.INTER_NEAREST)
            return mask_u8
        except Exception as exc:
            logger.debug(f"LM3SpecimenSegmenter inference fallback: {exc}")
            return None

    def __call__(self, image_bgr: np.ndarray, **kwargs: Any) -> Optional[np.ndarray]:
        return self.predict_mask(image_bgr, **kwargs)
