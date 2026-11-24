#!/usr/bin/env python3
"""
Step 02 whole-sheet LeafMachine3 diagnostic overlay exporter, 2D medoid label
repulsion helper, and hollow-marker Solidity vs. Aspect Ratio zone panel renderer.
"""

from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.step02_fig_style_assets import (
    _parse_bool_val,
    _resolve_asset_path,
    _resolve_facade_hook,
    _standardize_packera_species_py,
    format_icn_taxon_mathtext,
)

logger = logging.getLogger("SegmentAndExtract")

_ORIG_RESOLVE_ASSET_PATH = _resolve_asset_path


def export_lm3_sheet_diagnostic_overlays(
    extracted_df: pd.DataFrame,
    repro_df: Optional[pd.DataFrame],
    pheno_df: Optional[pd.DataFrame],
    figures_dir: Union[str, Path],
    vouchers_df: Optional[pd.DataFrame] = None,
    max_sheets: int = 24,
) -> List[Path]:
    """
    Export capped whole-sheet LeafMachine3 diagnostic overlays (`lm3_sheet_overlays/*.png`)
    at 300 DPI for priority sheets (Nomenclatural Types, DBSCAN multi-plant ramets K >= 2,
    basal+cauline heteroblastic sheets, and ruler-calibrated exemplars) whenever full scans exist.
    """
    resolve_fn = _resolve_facade_hook("_resolve_asset_path", _resolve_asset_path, _ORIG_RESOLVE_ASSET_PATH)
    out_dir = Path(figures_dir)
    sheet_dir = out_dir / "lm3_sheet_overlays"
    if extracted_df is None or len(extracted_df) == 0 or "catalogNumber" not in extracted_df.columns:
        return []

    df = extracted_df.copy()
    df["cat_clean"] = df["catalogNumber"].astype(str).str.strip()

    v_map_img: Dict[str, str] = {}
    v_map_tier: Dict[str, str] = {}
    v_map_type: Dict[str, bool] = {}
    if vouchers_df is not None and len(vouchers_df) > 0 and "catalogNumber" in vouchers_df.columns:
        for _, vr in vouchers_df.iterrows():
            c_k = str(vr.get("catalogNumber", "")).strip()
            if not c_k:
                continue
            if pd.notna(vr.get("image_path")):
                v_map_img[c_k] = str(vr.get("image_path")).strip()
            if pd.notna(vr.get("determiner_tier")):
                v_map_tier[c_k] = str(vr.get("determiner_tier")).strip()
            v_map_type[c_k] = _parse_bool_val(vr.get("is_type_specimen", False)) or (
                str(vr.get("type_status", "")).strip().lower() not in {"", "none", "nan", "na", "not a type"}
            )

    pheno_map: Dict[str, str] = {}
    if pheno_df is not None and len(pheno_df) > 0 and "catalogNumber" in pheno_df.columns:
        st_col = next(
            (c for c in ("voucher_phenological_state", "phenological_state", "pheno_state") if c in pheno_df.columns),
            None,
        )
        if st_col is not None:
            for _, pr in pheno_df.iterrows():
                pheno_map[str(pr["catalogNumber"]).strip()] = str(pr.get(st_col, "")).strip()

    # Score candidate sheets so Types, multi-plant DBSCAN ramets, and heteroblastic sheets come first
    candidates: List[Tuple[int, str]] = []
    for cat_id, grp in df.groupby("cat_clean", sort=False):
        if not cat_id:
            continue
        img_p = resolve_fn(v_map_img.get(cat_id, grp.iloc[0].get("image_path", "")), out_dir)
        if img_p is None:
            for ext in (".webp", ".jpg", ".jpeg", ".png", ".tif"):
                cand_raw = resolve_fn(f"raw_vouchers/{cat_id}{ext}", out_dir)
                if cand_raw is not None:
                    img_p = cand_raw
                    break
        if img_p is None or not img_p.exists():
            continue

        is_t = v_map_type.get(cat_id, False) or any(
            _parse_bool_val(x) for x in grp.get("is_type_specimen", pd.Series(dtype=object))
        )
        pids = set(
            int(x)
            for x in pd.to_numeric(grp.get("plant_individual_id", pd.Series([1])), errors="coerce").dropna()
            if int(x) > 0
        )
        zones = set(str(z).strip().lower() for z in grp.get("leaf_zone", pd.Series(dtype=object)).dropna())
        is_cal = any(_parse_bool_val(x) for x in grp.get("is_calibrated_mm", pd.Series(dtype=object)))
        prio = (100 if is_t else 0) + (50 if len(pids) >= 2 else 0) + (25 if {"basal", "cauline"}.issubset(zones) else 0) + (10 if is_cal else 0) + min(len(grp), 9)
        candidates.append((prio, cat_id))

    if not candidates:
        return []

    candidates.sort(key=lambda item: (-item[0], item[1]))
    selected_cats = [cat_id for _, cat_id in candidates[: int(max_sheets)]]
    sheet_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for cat_id in selected_cats:
        grp = df[df["cat_clean"] == cat_id]
        img_p = resolve_fn(v_map_img.get(cat_id, grp.iloc[0].get("image_path", "")), out_dir)
        if img_p is None:
            for ext in (".webp", ".jpg", ".jpeg", ".png", ".tif"):
                cand_raw = resolve_fn(f"raw_vouchers/{cat_id}{ext}", out_dir)
                if cand_raw is not None:
                    img_p = cand_raw
                    break
        if img_p is None or not img_p.exists():
            continue
        bgr = cv2.imread(str(img_p), cv2.IMREAD_COLOR)
        if bgr is None or bgr.size == 0:
            continue

        H_orig, W_orig = bgr.shape[:2]
        max_dim = 1800
        disp_scale = min(1.0, float(max_dim) / float(max(H_orig, W_orig, 1)))
        if disp_scale < 1.0:
            canvas = cv2.resize(
                bgr,
                (max(1, int(round(W_orig * disp_scale))), max(1, int(round(H_orig * disp_scale)))),
                interpolation=cv2.INTER_AREA,
            )
        else:
            canvas = bgr.copy()

        # Draw per-ramet envelopes when bbox coordinates exist
        b_cols = ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")
        has_bboxes = all(c in grp.columns for c in b_cols)
        ramet_boxes: Dict[int, List[Tuple[int, int, int, int]]] = {}

        for _, lrow in grp.iterrows():
            if not has_bboxes or any(pd.isna(lrow.get(c)) for c in b_cols):
                continue
            x1 = int(round(float(lrow["bbox_x1"]) * disp_scale))
            y1 = int(round(float(lrow["bbox_y1"]) * disp_scale))
            x2 = int(round(float(lrow["bbox_x2"]) * disp_scale))
            y2 = int(round(float(lrow["bbox_y2"]) * disp_scale))
            lz = str(lrow.get("leaf_zone", "basal")).strip().lower()
            gates_ok = _parse_bool_val(lrow.get("lm3_gates_pass", True), default=True)
            # BGR colors: Basal=#B27200 (blue), Cauline=#739E00 (green), Flagged=#888888
            col_bgr = (178, 114, 0) if lz != "cauline" else (115, 158, 0)
            if not gates_ok:
                col_bgr = (130, 130, 130)
            cv2.rectangle(canvas, (x1, y1), (x2, y2), col_bgr, 2, cv2.LINE_AA)
            pid = int(pd.to_numeric(lrow.get("plant_individual_id", 1), errors="coerce") or 1)
            ramet_boxes.setdefault(pid, []).append((x1, y1, x2, y2))
            lbl = f"L{int(lrow.get('leaf_id', 1))}:{lz[:3].upper()}"
            cv2.putText(canvas, lbl, (x1 + 3, max(14, y1 - 4)), cv2.FONT_HERSHEY_SIMPLEX, 0.42, col_bgr, 1, cv2.LINE_AA)

        if len(ramet_boxes) >= 2:
            for pid, boxes in sorted(ramet_boxes.items()):
                rx1 = max(2, min(b[0] for b in boxes) - 10)
                ry1 = max(2, min(b[1] for b in boxes) - 24)
                rx2 = min(canvas.shape[1] - 2, max(b[2] for b in boxes) + 10)
                ry2 = min(canvas.shape[0] - 2, max(b[3] for b in boxes) + 10)
                cv2.rectangle(canvas, (rx1, ry1), (rx2, ry2), (0, 94, 213), 2, cv2.LINE_AA)
                tag_txt = f"Ramet #{pid}"
                (tw, th), _ = cv2.getTextSize(tag_txt, cv2.FONT_HERSHEY_SIMPLEX, 0.48, 1)
                by_top = max(4, ry1 - th - 10)
                cv2.rectangle(canvas, (rx1, by_top), (rx1 + tw + 12, by_top + th + 8), (0, 94, 213), -1)
                cv2.putText(
                    canvas, tag_txt, (rx1 + 6, by_top + th + 2),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.48, (255, 255, 255), 1, cv2.LINE_AA,
                )

        # Two-line non-clipping header banner summarizing voucher metadata & LM3 calibration
        sp_raw = grp.iloc[0].get("species_standardized", grp.iloc[0].get("scientificName", ""))
        sp_lbl = _standardize_packera_species_py(sp_raw, grp.iloc[0].get("infraspecific_taxon", None)).replace("×", "x")
        px_mm = pd.to_numeric(grp.iloc[0].get("pixels_per_mm", np.nan), errors="coerce")
        sc_src = str(grp.iloc[0].get("scale_source", "uncalibrated")).strip()
        cal_txt = f"{float(px_mm):.2f} px/mm ({sc_src})" if pd.notna(px_mm) and float(px_mm) > 0 else sc_src
        ph_txt = pheno_map.get(cat_id, "unscored")
        tier_txt = v_map_tier.get(cat_id, str(grp.iloc[0].get("determiner_tier", "")))
        banner_h = 64
        banner = np.full((banner_h, canvas.shape[1], 3), 248, dtype=np.uint8)
        cv2.line(banner, (0, banner_h - 1), (canvas.shape[1], banner_h - 1), (190, 195, 202), 1)
        line1 = f"{cat_id}  |  {sp_lbl}  |  Authority: {tier_txt}"
        line2 = f"Leaves: {len(grp)} (DBSCAN Ramets: {max(1, len(ramet_boxes))})  |  Phenology: {ph_txt}  |  Scale: {cal_txt}"
        fscale = min(0.50, max(0.36, canvas.shape[1] / 1750.0))
        cv2.putText(banner, line1, (14, 25), cv2.FONT_HERSHEY_SIMPLEX, fscale, (20, 24, 28), 1, cv2.LINE_AA)
        cv2.putText(banner, line2, (14, 50), cv2.FONT_HERSHEY_SIMPLEX, fscale * 0.94, (55, 65, 81), 1, cv2.LINE_AA)
        composed = np.vstack([banner, canvas])

        out_overlay = sheet_dir / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', cat_id)}_lm3_sheet_overlay.png"
        tmp_ov = sheet_dir / f".{out_overlay.stem}_{os.getpid()}.png"
        cv2.imwrite(str(tmp_ov), composed)
        tmp_ov.replace(out_overlay)
        written.append(out_overlay)

    if written:
        logger.info(f"Exported {len(written)} whole-sheet LeafMachine3 diagnostic overlays to {sheet_dir}")
    return written


def _repel_medoid_annotations(
    ax: Any,
    points_xy_label_col: List[Tuple[float, float, str, str]],
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> None:
    """Place direct centroid/medoid labels inside plot bounds with 2D iterative box repulsion."""
    if not points_xy_label_col:
        return
    x_span = max(xlim[1] - xlim[0], 1e-6)
    y_span = max(ylim[1] - ylim[0], 1e-6)
    n_lbl = len(points_xy_label_col)
    pos = np.array([[p[0], p[1] + 0.035 * y_span] for p in points_xy_label_col], dtype=float)
    widths = np.array([min(0.42, max(0.16, 0.0125 * len(re.sub(r"[$\\{}]|it", "", p[2])))) * x_span for p in points_xy_label_col])
    box_h = 0.072 * y_span

    for _ in range(80):
        for i in range(n_lbl):
            pos[i, 0] = np.clip(pos[i, 0], xlim[0] + 0.52 * widths[i], xlim[1] - 0.52 * widths[i])
            pos[i, 1] = np.clip(pos[i, 1], ylim[0] + 0.65 * box_h, ylim[1] - 0.65 * box_h)
        for i in range(n_lbl):
            for j in range(i + 1, n_lbl):
                dx = pos[i, 0] - pos[j, 0]
                dy = pos[i, 1] - pos[j, 1]
                req_x = 0.54 * (widths[i] + widths[j])
                req_y = 1.05 * box_h
                if abs(dx) < req_x and abs(dy) < req_y:
                    push_y = 0.52 * (req_y - abs(dy)) * (1.0 if dy >= 0 else -1.0)
                    pos[i, 1] += push_y
                    pos[j, 1] -= push_y
                    push_x = 0.18 * (req_x - abs(dx)) * (1.0 if dx >= 0 else -1.0)
                    pos[i, 0] += push_x
                    pos[j, 0] -= push_x

    for idx, (x0, y0, lbl, col) in enumerate(points_xy_label_col):
        xk = float(np.clip(pos[idx, 0], xlim[0] + 0.52 * widths[idx], xlim[1] - 0.52 * widths[idx]))
        yk = float(np.clip(pos[idx, 1], ylim[0] + 0.65 * box_h, ylim[1] - 0.65 * box_h))
        ax.annotate(
            lbl, xy=(x0, y0), xytext=(xk, yk), ha="center", va="center", fontsize=7.2, color="#1f2937",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=col, alpha=0.90, linewidth=0.75),
            arrowprops=dict(arrowstyle="-", color=col, lw=0.65, shrinkA=2, shrinkB=3), zorder=8,
        )


_ORIG_REPEL_MEDOID = _repel_medoid_annotations


def _draw_solidity_vs_ar_zone_panel(
    ax: Any,
    zone_df: pd.DataFrame,
    active_taxa: List[str],
    style_map: Dict[str, Tuple[str, str]],
    panel_title: str,
    show_legend: bool = True,
) -> None:
    """Render hollow-marker Solidity vs. Aspect Ratio panel with reniform L/W < 1.0 dynamic scaling."""
    sol_vals = pd.to_numeric(zone_df["solidity_clean"], errors="coerce")
    ar_vals = pd.to_numeric(zone_df["ar_clean"], errors="coerce")
    ok = sol_vals.notna() & ar_vals.notna() & np.isfinite(sol_vals) & np.isfinite(ar_vals)
    sub = zone_df.loc[ok].copy()
    if len(sub) == 0:
        ax.text(0.5, 0.5, "No leaves in zone", ha="center", va="center", transform=ax.transAxes, fontsize=9, color="#6b7280")
        ax.set_title(panel_title, fontsize=10, fontweight="bold", loc="left")
        return

    x_min = max(0.15, float(sub["solidity_clean"].min()) - 0.05)
    x_max = min(1.02, float(sub["solidity_clean"].max()) + 0.03)
    if x_max - x_min < 0.20:
        x_min = max(0.15, x_max - 0.25)

    min_ar = float(sub["ar_clean"].min())
    max_ar = float(sub["ar_clean"].max())
    # Preserve reniform/orbicular blades (L/W < 1.0) without clipping at 1.0
    y_min = max(0.30, min(0.92, min_ar - 0.12)) if min_ar < 1.0 else max(0.75, min_ar - 0.18)
    y_max = max(y_min + 0.8, min(10.0, max_ar + 0.35))

    # Subtle morphological regime guides (Dissected vs. Entire/Undivided; Reniform L/W = 1.0)
    if x_min < 0.75 < x_max:
        ax.axvspan(x_min, 0.75, color="#f3f4f6", alpha=0.55, zorder=0)
        ax.axvline(0.75, color="#9ca3af", linestyle="--", linewidth=0.8, zorder=1)
        ax.text(
            x_min + 0.01 * (x_max - x_min), y_max - 0.05 * (y_max - y_min),
            "Lobed / Pinnatifid", fontsize=7.5, color="#4b5563", va="top",
        )
    if y_min < 1.0 < y_max:
        ax.axhline(1.0, color="#9ca3af", linestyle=":", linewidth=0.8, zorder=1)

    n_pts = len(sub)
    alpha_pt = 0.82 if n_pts <= 150 else (0.68 if n_pts <= 600 else 0.52)
    medoid_labels: List[Tuple[float, float, str, str]] = []

    for tx in active_taxa:
        tx_rows = sub[sub["species_clean"] == tx]
        if len(tx_rows) == 0:
            continue
        col, mk = style_map[tx]
        lbl_math = f"{format_icn_taxon_mathtext(tx, abbreviate_genus=(len(active_taxa) > 6))} (n={len(tx_rows)})"
        ax.scatter(
            tx_rows["solidity_clean"], tx_rows["ar_clean"], s=36, marker=mk,
            facecolors="none", edgecolors=col, linewidths=1.0, alpha=alpha_pt, label=lbl_math, zorder=4,
        )
        med_x = float(tx_rows["solidity_clean"].median())
        med_y = float(tx_rows["ar_clean"].median())
        medoid_labels.append((med_x, med_y, format_icn_taxon_mathtext(tx, abbreviate_genus=True), col))

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("Lamina Solidity (Area / Convex Hull Area)", fontsize=9)
    ax.set_ylabel("Anatomical Blade Aspect Ratio (L / W)", fontsize=9)
    ax.set_title(panel_title, fontsize=10, fontweight="bold", loc="left")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)

    if len(active_taxa) > 10:
        repel_fn = _resolve_facade_hook("_repel_medoid_annotations", _repel_medoid_annotations, _ORIG_REPEL_MEDOID)
        repel_fn(ax, medoid_labels, (x_min, x_max), (y_min, y_max))
    elif show_legend and len(active_taxa) >= 1:
        ncol_leg = 2 if len(active_taxa) > 5 else 1
        ax.legend(
            loc="best", fontsize=7.3, frameon=True, facecolor="white",
            framealpha=0.88, edgecolor="#d1d5db", ncol=ncol_leg,
        )


__all__ = [
    "export_lm3_sheet_diagnostic_overlays",
    "_repel_medoid_annotations",
    "_draw_solidity_vs_ar_zone_panel",
]
