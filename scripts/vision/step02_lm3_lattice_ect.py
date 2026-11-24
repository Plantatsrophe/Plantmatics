#!/usr/bin/env python3
"""
Step 02 LeafMachine3 CPU millimeter ruler tick-lattice solver and
Module 18 Euler Characteristic Transform (ECT) descriptor computation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import cv2
import numpy as np

_MP_PRIOR_SLOPE_PX_PER_CM: float = 27.112084


def _acorr_unbiased_1d(signal: np.ndarray, max_lag: int) -> np.ndarray:
    """Compute normalized unbiased 1D autocorrelation up to max_lag using FFT."""
    n = int(signal.size)
    if n < max_lag + 8:
        return np.zeros(0, dtype=np.float64)
    x = signal.astype(np.float64) - float(np.mean(signal))
    denom_zero = float(np.dot(x, x))
    if denom_zero <= 1e-9:
        return np.zeros(0, dtype=np.float64)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    fx = np.fft.rfft(x, n=nfft)
    raw = np.fft.irfft(fx * np.conj(fx), n=nfft)[: max_lag + 1]
    counts = np.arange(n, n - max_lag - 1, -1, dtype=np.float64)
    scale = (denom_zero / float(n)) * np.maximum(counts, 1.0)
    return np.clip(raw / scale, -1.0, 1.0)


def _parabolic_peak_lag(ac: np.ndarray, k_int: int) -> float:
    """Refine discrete autocorrelation peak position via 3-point quadratic fit."""
    if k_int <= 1 or k_int >= int(ac.size) - 1:
        return float(k_int)
    denom = 2.0 * (float(ac[k_int - 1]) - 2.0 * float(ac[k_int]) + float(ac[k_int + 1]))
    if denom >= -1e-9:
        return float(k_int)
    delta = (float(ac[k_int - 1]) - float(ac[k_int + 1])) / denom
    return float(k_int) + float(np.clip(delta, -0.5, 0.5))


def _scan_margin_strips_for_mm_lattice(
    gray_roi: np.ndarray,
    horizontal: bool,
    p_min: float,
    p_max: float,
    verified_ruler_crop: bool = False,
) -> Optional[tuple[float, int, float]]:
    """Scan parallel strips across a margin ROI for strict hierarchical 1-5-10 mm ticks."""
    roi = gray_roi if horizontal else gray_roi.T
    rh, rw = roi.shape
    min_rw = int(max(180.0, 14.0 * p_min)) if verified_ruler_crop else int(max(180.0, 18.0 * p_max))
    if rh < 16 or rw < min_rw:
        return None
    p_lo = int(np.floor(p_min))
    p_hi = int(np.ceil(p_max))
    min_p_lo = 3 if verified_ruler_crop else 4
    if p_lo < min_p_lo or p_hi <= p_lo + 1:
        return None
    strip_h = max(6, min(14, rh // 10))
    step_y = max(3, strip_h // 2)
    win_w = rw if rw < 2400 else max(1200, int(0.45 * rw))
    step_x = rw if win_w == rw else max(400, int(0.25 * rw))
    bg_half = max(3, int(round(0.45 * 0.5 * (p_min + p_max))))

    min_win_rw = int(max(120.0, 12.0 * p_min)) if verified_ruler_crop else int(15.0 * p_max)
    candidates: list[tuple[float, int, int, int, float]] = []
    for x0 in range(0, max(1, rw - win_w + 1), step_x):
        x1 = min(rw, x0 + win_w)
        if x1 - x0 < min_win_rw:
            continue
        for y0 in range(0, max(1, rh - strip_h + 1), step_y):
            patch = roi[y0 : y0 + strip_h, x0:x1]
            prof = np.mean(patch, axis=0)
            cs = np.cumsum(np.pad(prof, (1, 0), mode="constant"))
            idx = np.arange(prof.size)
            lo_idx = np.maximum(0, idx - bg_half)
            hi_idx = np.minimum(prof.size, idx + bg_half + 1)
            local_bg = (cs[hi_idx] - cs[lo_idx]) / np.maximum(1, hi_idx - lo_idx)
            # Test dominant polarity (dark ticks on light scale or bright ticks on dark scale)
            pos_sig = np.maximum(0.0, prof - local_bg)
            neg_sig = np.maximum(0.0, local_bg - prof)
            sig = neg_sig if float(np.sum(neg_sig * neg_sig)) >= float(np.sum(pos_sig * pos_sig)) else pos_sig
            if float(np.std(sig)) < 1.5 or float(np.mean(sig > 0.5 * np.max(sig))) > 0.45:
                continue
            max_lag = min(sig.size // 4, max(p_hi * 4 + 4, 48))
            ac = _acorr_unbiased_1d(sig, max_lag=max_lag)
            if verified_ruler_crop:
                p_hi_eff = min(p_hi, ac.size - 3)
                if p_hi_eff < p_lo:
                    continue
            else:
                if ac.size <= p_hi + 2:
                    continue
                p_hi_eff = p_hi
            search = ac[p_lo : p_hi_eff + 1]
            rel_idx = int(np.argmax(search))
            k = p_lo + rel_idx
            peak_val = float(ac[k])
            if verified_ruler_crop and k >= 2 * p_lo:
                for div in (3, 2):
                    c_sub = int(round(k / float(div)))
                    s_lo, s_hi = max(p_lo, c_sub - 1), min(p_hi_eff, c_sub + 1)
                    if s_lo <= s_hi:
                        sub_k = s_lo + int(np.argmax(ac[s_lo : s_hi + 1]))
                        if (
                            1 < sub_k < ac.size - 2
                            and ac[sub_k] > ac[sub_k - 1]
                            and ac[sub_k] > ac[sub_k + 1]
                            and float(ac[sub_k]) >= max(0.24, 0.72 * peak_val)
                        ):
                            k = sub_k
                            peak_val = float(ac[sub_k])
                            break
            if peak_val < 0.24 or k <= 1 or k >= ac.size - 2:
                continue
            # Anti-plateau and local concavity guard
            if not (ac[k] > ac[k - 1] and ac[k] > ac[k + 1]):
                continue
            p_sub = _parabolic_peak_lag(ac, k)
            if p_sub < p_min or p_sub > p_max:
                continue
            # Require coherent harmonics at 2*p and 3*p across >= 15 cycles
            if 3 * k + 3 >= ac.size:
                continue
            h2 = float(np.max(ac[max(2, 2 * k - 2) : 2 * k + 3]))
            h3 = float(np.max(ac[max(3, 3 * k - 3) : 3 * k + 4]))
            if h2 < 0.16 or h3 < 0.12:
                continue
            candidates.append((p_sub, x0, x1, y0, peak_val))

    if len(candidates) < 3:
        return None
    # Group nearby spatial windows and keep tightest cluster with >= 3 parallel strips
    best_cluster: Optional[tuple[float, int, float]] = None
    for x0_ref in sorted({c[1] for c in candidates}):
        sub = [c for c in candidates if abs(c[1] - x0_ref) <= max(64, int(0.10 * rw))]
        if len(sub) < 3:
            continue
        periods = np.array([c[0] for c in sub], dtype=np.float64)
        med_p = float(np.median(periods))
        agree = [c for c in sub if abs(c[0] - med_p) <= 0.25]
        if len(agree) < 3:
            continue
        agree_p = np.array([c[0] for c in agree], dtype=np.float64)
        if float(np.std(agree_p)) > 0.12:
            continue
        # Decimal (5x/10x) vs imperial (4x/8x) check across integrated tick band
        y_min = max(0, min(c[3] for c in agree) - strip_h)
        y_max = min(rh, max(c[3] for c in agree) + 2 * strip_h)
        x1_ref = max(c[2] for c in agree)
        full_patch = roi[y_min:y_max, x0_ref:x1_ref]
        f_prof = np.mean(full_patch, axis=0)
        f_cs = np.cumsum(np.pad(f_prof, (1, 0), mode="constant"))
        f_idx = np.arange(f_prof.size)
        f_lo = np.maximum(0, f_idx - bg_half)
        f_hi = np.minimum(f_prof.size, f_idx + bg_half + 1)
        f_bg = (f_cs[f_hi] - f_cs[f_lo]) / np.maximum(1, f_hi - f_lo)
        f_pos = np.maximum(0.0, f_prof - f_bg)
        f_neg = np.maximum(0.0, f_bg - f_prof)
        f_sig = f_neg if float(np.sum(f_neg * f_neg)) >= float(np.sum(f_pos * f_pos)) else f_pos
        p_est = float(np.median(agree_p))
        need_lag = int(np.ceil(10.5 * p_est)) + 2
        if f_sig.size > need_lag + 16:
            f_ac = _acorr_unbiased_1d(f_sig, max_lag=need_lag)
            if f_ac.size > int(round(10.0 * p_est)) + 2:
                def _win_peak(mult: int) -> float:
                    c_lag = int(round(mult * p_est))
                    tol = max(1, int(round(0.06 * c_lag)))
                    return float(np.max(f_ac[max(1, c_lag - tol) : min(f_ac.size, c_lag + tol + 1)]))

                h4, h5, h8, h10 = _win_peak(4), _win_peak(5), _win_peak(8), _win_peak(10)
                if (h4 + h8) > 1.06 * (h5 + h10) and h4 > h5:
                    continue
                if h10 >= 0.14:
                    c10 = int(round(10.0 * p_est))
                    tol10 = max(1, int(round(0.05 * c10)))
                    sub10_slice = f_ac[max(1, c10 - tol10) : min(f_ac.size, c10 + tol10 + 1)]
                    k10 = max(1, c10 - tol10) + int(np.argmax(sub10_slice))
                    if 1 < k10 < f_ac.size - 1 and f_ac[k10] > f_ac[k10 - 1] and f_ac[k10] > f_ac[k10 + 1]:
                        p_refined = _parabolic_peak_lag(f_ac, k10) / 10.0
                        if abs(p_refined - p_est) <= 0.15:
                            p_est = p_refined
        score = float(np.mean([c[4] for c in agree])) * float(len(agree))
        if best_cluster is None or score > best_cluster[2]:
            best_cluster = (p_est, len(agree), score)
    return best_cluster


def solve_ruler_tick_lattice_px_per_mm(
    image_bgr: np.ndarray,
    mp_prior_px_per_mm: Optional[float] = None,
    verified_ruler_crop: bool = False,
) -> Optional[float]:
    """
    Detect a metric millimeter tick lattice in outer sheet margins or verified
    ruler crops on CPU. Enforces non-plateau parabolic peaks, >=3 consistent
    parallel strips (std <= 0.12 px), decimal 5x/10x vs imperial 4x/8x
    discrimination, and anatomical megapixel/sheet-width bounds.
    """
    if image_bgr is None or not isinstance(image_bgr, np.ndarray) or image_bgr.ndim < 2 or image_bgr.size == 0:
        return None
    sh, sw = int(image_bgr.shape[0]), int(image_bgr.shape[1])
    if verified_ruler_crop:
        if min(sh, sw) < 16 or max(sh, sw) < 180:
            return None
    elif sh < 400 or sw < 400:
        return None
    if image_bgr.ndim == 3 and image_bgr.shape[2] >= 3:
        gray = (
            0.114 * image_bgr[:, :, 0].astype(np.float32)
            + 0.587 * image_bgr[:, :, 1].astype(np.float32)
            + 0.299 * image_bgr[:, :, 2].astype(np.float32)
        )
    else:
        gray = image_bgr[:, :, 0].astype(np.float32) if image_bgr.ndim == 3 else image_bgr.astype(np.float32)

    if mp_prior_px_per_mm is not None and np.isfinite(float(mp_prior_px_per_mm)) and float(mp_prior_px_per_mm) > 0.5:
        mp_est = float(mp_prior_px_per_mm)
    else:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        slope = float(
            getattr(getattr(_lm3_eng, "LM3ScaleCalibrator", None), "MP_PRIOR_SLOPE_PX_PER_CM", _MP_PRIOR_SLOPE_PX_PER_CM)
        )
        mp_est = (slope * float(np.sqrt((sh * sw) / 1e6))) / 10.0

    if verified_ruler_crop:
        p_min = max(3.5, 0.45 * mp_est)
        p_max = max(p_min + 0.5, 2.20 * mp_est)
        rois: List[Tuple[np.ndarray, bool]] = [
            (gray, True),
            (gray, False),
        ]
        if sh >= 400 and sw >= 400:
            m_h = max(32, int(0.18 * sh))
            m_w = max(32, int(0.18 * sw))
            rois.extend([
                (gray[:m_h, :], True),
                (gray[sh - m_h :, :], True),
                (gray[:, :m_w], False),
                (gray[:, sw - m_w :], False),
            ])
    else:
        sheet_floor = 0.85 * (float(min(sh, sw)) / 292.0) if min(sh, sw) >= 1200 else 0.0
        p_min = max(4.0, 0.80 * mp_est, sheet_floor)
        p_max = max(p_min + 0.5, 1.25 * mp_est)
        m_h = max(32, int(0.18 * sh))
        m_w = max(32, int(0.18 * sw))
        rois = [
            (gray[:m_h, :], True),
            (gray[sh - m_h :, :], True),
            (gray[:, :m_w], False),
            (gray[:, sw - m_w :], False),
        ]

    best: Optional[tuple[float, int, float]] = None
    for roi, horiz in rois:
        hit = _scan_margin_strips_for_mm_lattice(
            roi,
            horizontal=horiz,
            p_min=p_min,
            p_max=p_max,
            verified_ruler_crop=verified_ruler_crop,
        )
        if hit is not None and (best is None or hit[2] > best[2]):
            best = hit
    if best is None:
        return None
    return round(float(best[0]), 4)


def compute_leaf_ect_descriptor(
    mask: np.ndarray,
    num_dirs: Union[int, str, Path, None] = 360,
    num_thresh: int = 32,
    output_path: Optional[Union[str, Path]] = None,
    hole_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    angle_cw: Optional[float] = None,
    include_petiole: bool = False,
    include_holes: bool = True,
    already_oriented: bool = False,
) -> Optional[np.ndarray]:
    """
    Computes the LeafMachine3 Module 18 Euler Characteristic Transform (ECT) matrix
    (default `32 x 360` per `config/LM3_settings.yaml`) by delegating directly to
    `leafmachine3.core.ect_compute.compute_ect`, `ect_product_for`, and
    `leafmachine3.reporting.leaf_products.render_leaf_products`.
    """
    if mask is None or mask.size == 0 or int(np.count_nonzero(mask)) == 0:
        return None

    # Support legacy positional call `compute_leaf_ect_descriptor(mask, output_path)`
    if isinstance(num_dirs, (str, Path)) or num_dirs is None:
        if output_path is None and isinstance(num_dirs, (str, Path)):
            output_path = num_dirs
        eff_num_dirs = 360
    else:
        eff_num_dirs = int(num_dirs)
    eff_num_thresh = int(num_thresh)

    import scripts.vision.lm3_inference_engines as _lm3_eng

    _lm3_eng._ensure_ect_module_loaded()
    from leafmachine3.core import ect_compute
    from leafmachine3.reporting import leaf_products

    mask_u8 = ((mask > 0).astype(np.uint8)) * 255
    if not already_oriented and (angle_cw is not None or hole_mask is not None or petiole_mask is not None):
        h_m, w_m = mask_u8.shape[:2]
        dummy_rgb = np.zeros((h_m, w_m, 3), dtype=np.uint8)
        lam_bool = mask_u8 > 0
        hole_bool = (
            (hole_mask > 0)
            if (hole_mask is not None and hole_mask.size > 0 and hole_mask.shape[:2] == (h_m, w_m))
            else None
        )
        pet_bool = (
            (petiole_mask > 0)
            if (petiole_mask is not None and petiole_mask.size > 0 and petiole_mask.shape[:2] == (h_m, w_m))
            else None
        )
        if hole_bool is not None:
            lam_bool = lam_bool & (~hole_bool)
            sil_bool = lam_bool | hole_bool
        else:
            sil_bool = lam_bool
        prod_name = ect_compute.ect_product_for(bool(include_petiole), bool(include_holes))
        masks_dict: Dict[str, Optional[np.ndarray]] = {
            "lamina": lam_bool,
            "laminapet": (lam_bool | pet_bool) if pet_bool is not None else None,
            "silhouette": sil_bool,
            "silhouettepet": (sil_bool | pet_bool) if pet_bool is not None else None,
            "hole": hole_bool,
        }
        prods = leaf_products.render_leaf_products(
            dummy_rgb,
            masks_dict,
            angle_cw=float(angle_cw) if angle_cw is not None else None,
            want={prod_name},
        )
        selected = prods.get(prod_name)
        if selected is not None and selected.size > 0 and int(np.count_nonzero(selected)) > 0:
            mask_u8 = ((selected > 0).astype(np.uint8)) * 255

    res = ect_compute.compute_ect(mask_u8, num_dirs=eff_num_dirs, bound_radius=1.0, want_simple=False)
    if res is None or getattr(res, "ect_matrix", None) is None:
        return None

    ect_raw = np.asarray(res.ect_matrix, dtype=np.float32)
    thetas = getattr(res, "thetas", None)
    if thetas is not None and len(thetas) == ect_raw.shape[1]:
        ect_mat = np.asarray(ect_compute.tip_up_directions(ect_raw, thetas)[0], dtype=np.float32)
    else:
        ect_mat = ect_raw
    if ect_mat.shape != (eff_num_thresh, eff_num_dirs):
        ect_mat = cv2.resize(ect_mat, (eff_num_dirs, eff_num_thresh), interpolation=cv2.INTER_LINEAR).astype(
            np.float32
        )

    if output_path is not None:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        np.save(str(out_p), ect_mat)

    return ect_mat
