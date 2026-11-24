#!/usr/bin/env python3
"""
Read-only Step 02 leaf traits: enclosed interlobe void quantification, Packera
margin/base/terminal-lobe descriptors, and within-rosette abaxial pressed-leaf
flagging (`scripts/vision/step02_leaf_traits_utils.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import cv2
import numpy as np


def _safe_finite_float(val: Any) -> Optional[float]:
    if val is None:
        return None
    try:
        f = float(val)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _parse_point_xy(points: Optional[Mapping[str, Any]], key: str) -> Optional[Tuple[float, float]]:
    if not isinstance(points, Mapping) or key not in points:
        return None
    raw = points.get(key)
    if isinstance(raw, Mapping):
        conf = _safe_finite_float(raw.get("conf", raw.get("confidence")))
        x, y = _safe_finite_float(raw.get("x")), _safe_finite_float(raw.get("y"))
    elif isinstance(raw, (tuple, list, np.ndarray)) and len(raw) >= 2:
        conf = _safe_finite_float(raw[2]) if len(raw) >= 3 else None
        x, y = _safe_finite_float(raw[0]), _safe_finite_float(raw[1])
    else:
        return None
    if conf is not None and conf < 0.25:
        return None
    return (x, y) if (x is not None and y is not None) else None


def measure_enclosed_interlobe_void_frac(
    mask_u8: Optional[np.ndarray],
    hole_mask: Optional[np.ndarray] = None,
) -> float:
    """Measure enclosed non-hole background void area relative to lamina area without mutating `mask_u8`."""
    if mask_u8 is None or not isinstance(mask_u8, np.ndarray) or mask_u8.ndim < 2 or mask_u8.size == 0:
        return 0.0
    fg = (mask_u8[:, :, 0] > 0) if mask_u8.ndim == 3 else (mask_u8 > 0)
    fg_area = int(np.count_nonzero(fg))
    if fg_area < 16:
        return 0.0

    h, w = fg.shape[:2]
    num_bg, bg_labels, bg_stats, _ = cv2.connectedComponentsWithStats((~fg).astype(np.uint8), connectivity=8)
    if num_bg <= 1:
        return 0.0

    min_void_px = max(16, int(round(0.003 * float(fg_area))))
    hm_bool: Optional[np.ndarray] = None
    if hole_mask is not None and isinstance(hole_mask, np.ndarray) and hole_mask.ndim >= 2 and hole_mask.shape[:2] == (h, w):
        hm_bool = (hole_mask[:, :, 0] > 0) if hole_mask.ndim == 3 else (hole_mask > 0)

    unopened_void_px = 0
    for lbl in range(1, num_bg):
        x, y = int(bg_stats[lbl, cv2.CC_STAT_LEFT]), int(bg_stats[lbl, cv2.CC_STAT_TOP])
        cw, ch = int(bg_stats[lbl, cv2.CC_STAT_WIDTH]), int(bg_stats[lbl, cv2.CC_STAT_HEIGHT])
        comp_area = int(bg_stats[lbl, cv2.CC_STAT_AREA])
        if x == 0 or y == 0 or (x + cw) >= w or (y + ch) >= h or comp_area < min_void_px:
            continue
        if hm_bool is not None and float(np.count_nonzero((bg_labels == lbl) & hm_bool)) / float(max(comp_area, 1)) > 0.30:
            continue
        unopened_void_px += comp_area
    return round(float(unopened_void_px) / float(max(fg_area, 1)), 4)


def _savgol_quadratic_baseline(vals: np.ndarray, win: int) -> np.ndarray:
    """Compute a degree-2 Savitzky-Golay macro-outline baseline with polynomial edge interpolation."""
    n = int(vals.size)
    if n < 5:
        return vals.astype(np.float64, copy=True)
    w = min(max(5, int(win)), n)
    w = (w - 1 if w > 5 else min(n, 5)) if (w % 2 == 0) else w
    if w < 5 or w > n:
        return vals.astype(np.float64, copy=True)
    half = w // 2
    t = np.arange(-half, half + 1, dtype=np.float64)
    kernel = np.linalg.pinv(np.column_stack([np.ones_like(t), t, t * t]))[0]
    out = np.empty(n, dtype=np.float64)
    out[half : n - half] = np.convolve(vals, kernel, mode="valid")
    idx_w = np.arange(w, dtype=np.float64)
    out[:half] = np.polyval(np.polyfit(idx_w, vals[:w], 2), np.arange(half, dtype=np.float64))
    out[n - half :] = np.polyval(np.polyfit(idx_w, vals[n - w :], 2), np.arange(w - half, w, dtype=np.float64))
    return out


def _crossing_distance(res: np.ndarray, peak_idx: int, level: float, direction: int, max_steps: int) -> float:
    """Linearly interpolate distance from `peak_idx` to first crossing below `level`."""
    n, prev_val, idx = int(res.size), float(res[peak_idx]), peak_idx
    for steps in range(max_steps):
        nxt = idx + direction
        if nxt < 0 or nxt >= n:
            return float(steps)
        cur_val = float(res[nxt])
        if cur_val <= level:
            denom = prev_val - cur_val
            frac = (prev_val - level) / denom if abs(denom) > 1e-9 else 0.0
            return float(steps) + float(np.clip(frac, 0.0, 1.0))
        prev_val, idx = cur_val, nxt
    return float(max_steps)


def _side_trough(res: np.ndarray, peak_idx: int, peak_val: float, direction: int, span: int) -> float:
    n, low = int(res.size), peak_val
    stop = max(-1, peak_idx - span - 1) if direction < 0 else min(n, peak_idx + span + 1)
    for k in range(peak_idx + direction, stop, direction):
        v = float(res[k])
        if v < low:
            low = v
        if v > peak_val:
            break
    return low


def _measure_side_tooth_acuteness(res: np.ndarray, prom_thresh: float, search_span: int) -> List[float]:
    """Detect outward marginal teeth in high-pass residual `res` and score apex acuteness."""
    n = int(res.size)
    if n < 5:
        return []
    raw_peaks: List[int] = []
    i = 1
    while i < n - 1:
        if res[i] > res[i - 1]:
            j = i
            while j + 1 < n - 1 and res[j + 1] == res[i]:
                j += 1
            if res[j + 1] < res[i]:
                raw_peaks.append((i + j) // 2)
            i = j + 1
        else:
            i += 1
    if not raw_peaks:
        return []

    peaks: List[int] = []
    for p in raw_peaks:
        if peaks and (p - peaks[-1]) <= 2:
            if res[p] > res[peaks[-1]]:
                peaks[-1] = p
        else:
            peaks.append(p)

    scores: List[float] = []
    span = max(4, int(search_span))
    for p in peaks:
        peak_val = float(res[p])
        if peak_val < prom_thresh:
            continue
        trough = max(_side_trough(res, p, peak_val, -1, span), _side_trough(res, p, peak_val, +1, span))
        prom = peak_val - trough
        if prom < (1.25 * prom_thresh):
            continue
        lvl_50, lvl_80 = trough + 0.50 * prom, trough + 0.80 * prom
        w_50 = max(1.0, _crossing_distance(res, p, lvl_50, -1, span) + _crossing_distance(res, p, lvl_50, +1, span))
        w_80 = max(0.5, _crossing_distance(res, p, lvl_80, -1, span) + _crossing_distance(res, p, lvl_80, +1, span))
        k_q = max(1, int(round(0.25 * w_50)))
        d_q = max(0.0, peak_val - 0.5 * (float(res[max(0, p - k_q)]) + float(res[min(n - 1, p + k_q)]))) / max(prom, 1e-6)
        taper_ratio = max(0.0, (w_50 / w_80) - 1.0)
        scores.append(float(np.clip(4.0 * d_q * taper_ratio, 0.0, 5.0)))
    return scores


def extract_packera_leaf_margin_and_base_traits(
    upright_mask_u8: Optional[np.ndarray],
    oriented_points: Optional[Mapping[str, Tuple[float, float]]] = None,
) -> Dict[str, float]:
    """Extract read-only Packera lamina base, margin serration/acuteness, and terminal-lobe traits."""
    nan_out: Dict[str, float] = {
        "basal_cordate_index": np.nan, "margin_serration_index": np.nan,
        "margin_tooth_acuteness": np.nan, "terminal_lobe_area_frac": np.nan,
    }
    if upright_mask_u8 is None or not isinstance(upright_mask_u8, np.ndarray) or upright_mask_u8.ndim < 2 or upright_mask_u8.size == 0:
        return nan_out
    fg = (upright_mask_u8[:, :, 0] > 0) if upright_mask_u8.ndim == 3 else (upright_mask_u8 > 0)
    fg_area = int(np.count_nonzero(fg))
    if fg_area < 32:
        return nan_out

    ys, xs = np.where(fg)
    y_min, y_max, x_min, x_max = int(ys.min()), int(ys.max()), int(xs.min()), int(xs.max())
    blade_l, blade_w = max(float(y_max - y_min), 1.0), max(float(x_max - x_min), 1.0)
    sub = fg[y_min : y_max + 1, x_min : x_max + 1]
    bh, bw = sub.shape

    # 1. Basal cordate index (>0 cordate, ~0 truncate, <0 cuneate)
    pt_base = _parse_point_xy(oriented_points, "lamina_base") or _parse_point_xy(oriented_points, "base_center")
    if pt_base is not None and (x_min - 0.5 * blade_w) <= pt_base[0] <= (x_max + 0.5 * blade_w) and (y_min - 0.5 * blade_l) <= pt_base[1] <= (y_max + 0.5 * blade_l):
        x_base, y_base = pt_base
    else:
        x_mid = 0.5 * (float(x_min) + float(x_max))
        x_base, half_band, r_floor = x_mid, max(1.0, 0.05 * blade_w), int(math.floor(0.35 * blade_l))
        bot_rows = [
            float(y_min + int(nz_r[-1]))
            for c in range(bw)
            if abs(float(x_min + c) - x_mid) <= half_band
            for nz_r in (np.flatnonzero(sub[:, c]),)
            if nz_r.size > 0 and int(nz_r[-1]) >= r_floor
        ]
        y_base = float(np.median(bot_rows)) if bot_rows else float(y_max)

    r_start_45 = min(bh - 1, max(0, int(math.floor(0.55 * blade_l))))
    lat_cols = np.abs((np.arange(bw, dtype=np.float64) + float(x_min)) - x_base) >= (0.10 * blade_w)
    max_basal_lobe_y: Optional[float] = None
    if np.any(lat_cols):
        lat_ys = np.where(sub[r_start_45:, :][:, lat_cols])[0]
        if lat_ys.size > 0:
            max_basal_lobe_y = float(y_min + r_start_45 + int(lat_ys.max()))

    if max_basal_lobe_y is not None and max_basal_lobe_y > (y_base + 1.0):
        bci_raw = (max_basal_lobe_y - y_base) / blade_l
    else:
        y_at_25, thresh_w = float(y_max), 0.25 * blade_w
        for r in range(bh - 1, -1, -1):
            nz_c = np.flatnonzero(sub[r])
            if nz_c.size > 0 and float(int(nz_c[-1]) - int(nz_c[0])) >= thresh_w:
                y_at_25 = float(y_min + r)
                break
        bci_raw = (y_at_25 - float(y_max)) / blade_l
    basal_cordate_index = round(float(np.clip(bci_raw, -0.50, 0.50)), 4)

    # 2. Margin serration index & tooth acuteness over central 70% of blade
    r0, r1 = min(bh - 1, max(0, int(math.ceil(0.15 * blade_l)))), min(bh - 1, int(math.floor(0.85 * blade_l)))
    r1 = max(r0, r1)
    n_span = r1 - r0 + 1
    x_left, x_right = np.full(n_span, np.nan, dtype=np.float64), np.full(n_span, np.nan, dtype=np.float64)
    for idx, r in enumerate(range(r0, r1 + 1)):
        nz_c = np.flatnonzero(sub[r])
        if nz_c.size > 0:
            x_left[idx], x_right[idx] = float(nz_c[0]), float(nz_c[-1])

    valid_idx = np.flatnonzero(np.isfinite(x_left))
    margin_serration_index, margin_tooth_acuteness = np.nan, np.nan
    if valid_idx.size >= 12:
        if valid_idx.size < n_span:
            all_idx = np.arange(n_span, dtype=np.float64)
            x_left = np.interp(all_idx, valid_idx.astype(np.float64), x_left[valid_idx])
            x_right = np.interp(all_idx, valid_idx.astype(np.float64), x_right[valid_idx])
        win = max(9, int(round(0.18 * n_span)))
        win = win + 1 if (win % 2 == 0) else win
        if win > n_span:
            win = n_span if (n_span % 2 == 1) else (n_span - 1)
        res_left = (-x_left) - _savgol_quadratic_baseline(-x_left, win)
        res_right = x_right - _savgol_quadratic_baseline(x_right, win)
        all_res = np.concatenate([res_left, res_right])
        margin_serration_index = round(float(np.sqrt(np.mean(all_res * all_res)) / blade_w), 4)
        prom_thresh = max(0.75, 0.006 * blade_w)
        tooth_scores = _measure_side_tooth_acuteness(res_left, prom_thresh, win) + _measure_side_tooth_acuteness(res_right, prom_thresh, win)
        margin_tooth_acuteness = round(float(np.median(tooth_scores)), 4) if tooth_scores else 0.0

    # 3. Terminal lobe area fraction (1.0 on unlobed blades; <1.0 above deepest inter-lobe sinus)
    w_row = np.count_nonzero(sub, axis=1).astype(np.float64)
    r_lo = min(bh - 1, max(1, int(math.floor(0.25 * blade_l))))
    r_hi = min(bh - 2, max(r_lo, int(math.ceil(0.75 * blade_l))))
    terminal_lobe_area_frac = 1.0
    if r_hi >= r_lo and bh >= 6:
        best_r, best_w = None, float("inf")
        for r in range(r_lo, r_hi + 1):
            w_cur = float(w_row[r])
            w_ap = float(np.max(w_row[:r])) if r > 0 else 0.0
            w_bs = float(np.max(w_row[r + 1 :])) if (r + 1) < bh else 0.0
            if w_ap >= max(3.0, 1.35 * w_cur) and w_bs >= max(3.0, 1.25 * w_cur) and w_cur < best_w:
                best_w, best_r = w_cur, r
        if best_r is not None:
            terminal_lobe_area_frac = round(float(np.clip(float(np.count_nonzero(sub[: best_r + 1, :])) / float(fg_area), 0.0, 1.0)), 4)

    return {
        "basal_cordate_index": basal_cordate_index,
        "margin_serration_index": margin_serration_index,
        "margin_tooth_acuteness": margin_tooth_acuteness,
        "terminal_lobe_area_frac": terminal_lobe_area_frac,
    }


def evaluate_rosette_abaxial_leaf_flags(records: Sequence[MutableMapping[str, Any]]) -> None:
    """Flag flipped abaxial-pressed basal leaves within each `(catalogNumber, plant_individual_id)` rosette."""
    groups: Dict[Tuple[str, int], List[MutableMapping[str, Any]]] = {}
    for rec in records:
        rec["is_abaxial_pressed_leaf"] = False
        rec["abaxial_cie_l_sheet_delta"] = np.nan
        cat = str(rec.get("catalogNumber", ""))
        try:
            pid = int(rec.get("plant_individual_id", 0) or 0)
        except (TypeError, ValueError):
            pid = 0
        if str(rec.get("leaf_zone", "")).strip().lower() == "basal" and _safe_finite_float(rec.get("blade_cie_l_sheet_delta")) is not None:
            groups.setdefault((cat, pid), []).append(rec)

    for basal_finite in groups.values():
        if len(basal_finite) < 3:
            continue
        sibling_median = float(np.median([float(r["blade_cie_l_sheet_delta"]) for r in basal_finite]))
        for rec in basal_finite:
            l_delta = float(rec["blade_cie_l_sheet_delta"])
            a_delta, b_delta = _safe_finite_float(rec.get("blade_cie_a_delta")), _safe_finite_float(rec.get("blade_cie_b_delta"))
            if a_delta is not None and b_delta is not None and math.hypot(a_delta, b_delta) <= 18.0 and l_delta >= (sibling_median + 15.0):
                rec["is_abaxial_pressed_leaf"] = True
                rec["abaxial_cie_l_sheet_delta"] = round(l_delta, 4)
                rec["blade_cie_l_norm"] = np.nan
                rec["blade_cie_l_sheet_delta"] = np.nan


__all__ = [
    "measure_enclosed_interlobe_void_frac",
    "extract_packera_leaf_margin_and_base_traits",
    "evaluate_rosette_abaxial_leaf_flags",
]
