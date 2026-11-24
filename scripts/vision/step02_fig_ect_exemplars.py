#!/usr/bin/env python3
"""
Step 02 native LeafMachine3 silhouette/cutout card renderer, polar radial ECT
overlay renderer, and taxon-stratified exemplar cache exporter.
"""

from __future__ import annotations

import json
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import scripts.vision.lm3_compat  # noqa: F401
from scripts.vision.lm3_inference_engines import _ensure_ect_module_loaded
from scripts.vision.step02_fig_style_assets import (
    _is_valid_smooth_closed_mask,
    _load_mask_or_contour_u8,
    _parse_bool_val,
    _resolve_facade_hook,
    _standardize_packera_species_py,
    get_combinatorial_taxon_style,
)

logger = logging.getLogger("SegmentAndExtract")
_ORIG_IS_VALID_MASK, _ORIG_LOAD_MASK = _is_valid_smooth_closed_mask, _load_mask_or_contour_u8


def _render_lm3_silhouette_card(
    mask_u8: np.ndarray,
    rgb_patch: Optional[np.ndarray],
    accent_hex: str = "#0072B2",
    pixels_per_mm: Optional[float] = None,
    is_calibrated: bool = False,
    card_size: int = 320,
    prefer_rgb: bool = True,
) -> np.ndarray:
    """Render square RGB card with masked RGB cutout or silhouette, midrib vector, and 10 mm scale bar."""
    card = np.full((card_size, card_size, 3), 255, dtype=np.uint8)
    h, w = mask_u8.shape[:2]
    if h == 0 or w == 0 or int(np.count_nonzero(mask_u8)) == 0:
        return card

    scale = float(int(card_size * 0.80)) / float(max(h, w, 1))
    nh, nw = max(4, int(round(h * scale))), max(4, int(round(w * scale)))
    m_res = cv2.resize(mask_u8, (nw, nh), interpolation=cv2.INTER_NEAREST)
    ox, oy = (card_size - nw) // 2, (card_size - nh) // 2
    acc_rgb = tuple(int(accent_hex.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
    fg = m_res > 0

    patch_canvas = card[oy : oy + nh, ox : ox + nw]
    if prefer_rgb and rgb_patch is not None and rgb_patch.shape[:2] == (h, w):
        patch_canvas[fg] = cv2.resize(rgb_patch, (nw, nh), interpolation=cv2.INTER_AREA)[fg]
    else:
        patch_canvas[fg] = tuple(int(0.22 * c + 0.78 * 235) for c in acc_rgb)

    # Crisp LeafMachine3 contour outline
    cnts, _ = cv2.findContours((fg.astype(np.uint8) * 255), cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    if cnts:
        shifted = [c + np.array([[[ox, oy]]], dtype=np.int32) for c in cnts]
        cv2.drawContours(card, shifted, -1, acc_rgb, thickness=2, lineType=cv2.LINE_AA)

    # Tip-up midrib axis and basal lamina attachment point (p_base at bottom pole of upright mask)
    ys, xs = np.where(fg)
    if len(ys) >= 10:
        y_top, y_bot = int(ys.min()), int(ys.max())
        top_band = xs[ys <= y_top + max(2, int(0.06 * (y_bot - y_top + 1)))]
        bot_band = xs[ys >= y_bot - max(2, int(0.06 * (y_bot - y_top + 1)))]
        x_tip = int(round(float(np.median(top_band)))) + ox if len(top_band) else (ox + nw // 2)
        x_base = int(round(float(np.median(bot_band)))) + ox if len(bot_band) else (ox + nw // 2)
        y_tip, y_base = y_top + oy, y_bot + oy

        # Draw dashed midrib line and hollow p_base / p_apex markers
        n_seg = 14
        for s_i in range(0, n_seg, 2):
            t0, t1 = s_i / n_seg, min(1.0, (s_i + 1) / n_seg)
            p0 = (int(round(x_base + t0 * (x_tip - x_base))), int(round(y_base + t0 * (y_tip - y_base))))
            p1 = (int(round(x_base + t1 * (x_tip - x_base))), int(round(y_base + t1 * (y_tip - y_base))))
            cv2.line(card, p0, p1, (75, 85, 99), 1, cv2.LINE_AA)
        cv2.circle(card, (x_base, y_base), 5, (213, 94, 0), 2, cv2.LINE_AA)
        cv2.circle(card, (x_tip, y_tip), 4, acc_rgb, 1, cv2.LINE_AA)

    # Physical 10 mm scale bar on ruler-calibrated specimens
    if is_calibrated and pixels_per_mm is not None and np.isfinite(pixels_per_mm) and float(pixels_per_mm) > 0:
        bar_px, bar_mm = int(round(10.0 * float(pixels_per_mm) * scale)), 10
        if bar_px > int(card_size * 0.45):
            bar_px, bar_mm = int(round(5.0 * float(pixels_per_mm) * scale)), 5
        if 12 <= bar_px <= int(card_size * 0.60):
            bx0, bx1, by = card_size - bar_px - 14, card_size - 14, card_size - 14
            cv2.line(card, (bx0, by), (bx1, by), (30, 30, 30), 2, cv2.LINE_AA)
            cv2.line(card, (bx0, by - 3), (bx0, by + 3), (30, 30, 30), 1, cv2.LINE_AA)
            cv2.line(card, (bx1, by - 3), (bx1, by + 3), (30, 30, 30), 1, cv2.LINE_AA)
            cv2.putText(card, f"{bar_mm} mm", (max(6, bx0 - 2), by - 6), cv2.FONT_HERSHEY_SIMPLEX, 0.36, (30, 30, 30), 1, cv2.LINE_AA)

    return card


def _render_lm3_ect_rgb(mask_u8: np.ndarray, n_dirs: int = 160, size: int = 320) -> Optional[np.ndarray]:
    """Delegate polar radial ECT + contour rendering directly to LeafMachine3 `render_radial_ect_overlay`."""
    _ensure_ect_module_loaded()
    from leafmachine3.core.ect_compute import compute_ect
    from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log

    ect_res = compute_ect(mask_u8, num_dirs=int(n_dirs), bound_radius=1.0, want_simple=False)
    if ect_res is None or getattr(ect_res, "ect_matrix", None) is None:
        return None
    bgr = render_radial_ect_overlay(
        visual_log(ect_res.ect_matrix), ect_res.thetas, ect_res.thresholds,
        ect_res.bound_radius, getattr(ect_res, "outline_norm", None), cmap="viridis",
    )
    if bgr is None or bgr.size == 0:
        return None
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    return cv2.resize(rgb, (size, size), interpolation=cv2.INTER_AREA) if rgb.shape[:2] != (size, size) else rgb


_ORIG_RENDER_SIL, _ORIG_RENDER_ECT = _render_lm3_silhouette_card, _render_lm3_ect_rgb


def export_taxon_ect_exemplars(
    extracted_df: Optional[pd.DataFrame] = None,
    figures_dir: Union[str, Path] = "outputs/figures",
    zone: str = "basal",
    n_dirs: int = 160,
    max_taxa: Optional[int] = None,
    vouchers_csv: Optional[Union[str, Path]] = None,
    tables_dir: Optional[Union[str, Path]] = None,
    db_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Export authentic LeafMachine3 radial ECT overlays (`ect_overlays/<zone>_<slug>.png`)
    and native LeafMachine3 silhouette/RGB exemplar badges (`lm3_assets/<zone>_<slug>_*.png`)
    for ALL active taxa/ranks (unlimited by default when `max_taxa=None`).
    """
    is_valid_mask_fn = _resolve_facade_hook("_is_valid_smooth_closed_mask", _is_valid_smooth_closed_mask, _ORIG_IS_VALID_MASK)
    load_mask_fn = _resolve_facade_hook("_load_mask_or_contour_u8", _load_mask_or_contour_u8, _ORIG_LOAD_MASK)
    render_ect_fn = _resolve_facade_hook("_render_lm3_ect_rgb", _render_lm3_ect_rgb, _ORIG_RENDER_ECT)
    render_sil_fn = _resolve_facade_hook("_render_lm3_silhouette_card", _render_lm3_silhouette_card, _ORIG_RENDER_SIL)

    fig_dir = Path(figures_dir).resolve()
    overlays_dir, lm3_assets_dir = fig_dir / "ect_overlays", fig_dir / "lm3_assets"
    overlays_dir.mkdir(parents=True, exist_ok=True)
    lm3_assets_dir.mkdir(parents=True, exist_ok=True)
    zone_clean = str(zone).lower().strip()

    if extracted_df is None or len(extracted_df) == 0:
        cand_tables: List[Path] = []
        t_names = ("extracted_leaf_manifest.csv", "leaf_efa_harmonics_individual.csv", "leaf_efa_harmonics.csv")
        if tables_dir is not None:
            cand_tables.extend(Path(tables_dir) / n for n in t_names)
        for base in (fig_dir.parent.parent / "data" / "tables", fig_dir.parent / "data" / "tables", PROJECT_ROOT / "data" / "tables"):
            cand_tables.extend(base / n for n in t_names)
        for ct in cand_tables:
            if ct.exists() and ct.stat().st_size > 0:
                try:
                    extracted_df = pd.read_csv(ct, low_memory=False)
                    break
                except Exception:
                    pass

    if extracted_df is None or len(extracted_df) == 0:
        logger.warning(f"export_taxon_ect_exemplars: No extracted leaf data found to export for zone '{zone_clean}'.")
        return {}

    df_work = extracted_df.copy()
    if vouchers_csv is not None and Path(vouchers_csv).exists() and "catalogNumber" in df_work.columns:
        try:
            v_df = pd.read_csv(vouchers_csv, dtype=str)
            if "catalogNumber" in v_df.columns:
                for col_v in ("species_standardized", "scientificName", "infraspecific_taxon", "image_path"):
                    if col_v in v_df.columns and col_v not in df_work.columns:
                        v_map = dict(zip(v_df["catalogNumber"].astype(str).str.strip(), v_df[col_v]))
                        df_work[col_v] = df_work["catalogNumber"].astype(str).str.strip().map(v_map)
        except Exception:
            pass

    if "species_clean" not in df_work.columns:
        sp_col = next((c for c in ("species_standardized", "scientificName", "species_name", "species_raw", "taxon") if c in df_work.columns), None)
        inf_col = "infraspecific_taxon" if "infraspecific_taxon" in df_work.columns else None
        if sp_col is not None:
            df_work["species_clean"] = (
                [_standardize_packera_species_py(s, inf) for s, inf in zip(df_work[sp_col], df_work[inf_col])]
                if inf_col is not None
                else df_work[sp_col].apply(_standardize_packera_species_py)
            )
        else:
            df_work["species_clean"] = "Unassigned_Packera"

    if "leaf_zone" in df_work.columns:
        lz_norm = df_work["leaf_zone"].astype(str).str.strip().str.lower()
        zc_num = (
            pd.to_numeric(df_work["leaf_zone_confidence"], errors="coerce").fillna(1.0)
            if "leaf_zone_confidence" in df_work.columns
            else pd.Series(1.0, index=df_work.index)
        )
        pbr_num = (
            pd.to_numeric(df_work["petiole_to_blade_ratio"], errors="coerce")
            if "petiole_to_blade_ratio" in df_work.columns
            else pd.Series(np.nan, index=df_work.index)
        )
        if zone_clean == "cauline":
            df_zone = df_work[lz_norm == "cauline"].copy()
            if len(df_zone) >= 2:
                hi_c = df_zone[(zc_num.loc[df_zone.index] >= 0.72) & (pbr_num.loc[df_zone.index].isna() | (pbr_num.loc[df_zone.index] <= 0.30))]
                if len(hi_c) >= 1:
                    df_zone = hi_c.copy()
        else:
            df_zone = df_work[lz_norm.isin({"basal", "", "nan", "none"})].copy()
            if len(df_zone) == 0:
                df_zone = df_work.copy()
            if len(df_zone) >= 2:
                hi_b = df_zone[(zc_num.loc[df_zone.index] >= 0.72) & (pbr_num.loc[df_zone.index].isna() | (pbr_num.loc[df_zone.index] >= 0.15))]
                if len(hi_b) >= 1:
                    df_zone = hi_b.copy()
    else:
        df_zone = df_work.iloc[0:0].copy() if zone_clean == "cauline" else df_work.copy()

    if len(df_zone) == 0:
        with open(overlays_dir / f"manifest_{zone_clean}.json", "w", encoding="utf-8") as f_m:
            json.dump({}, f_m, indent=2)
        logger.info(f"Exported 0 authentic LeafMachine3 radial ECT overlays to {overlays_dir}")
        return {}

    sp_counts = df_zone["species_clean"].value_counts()
    valid_sp = [s for s in sp_counts.index if s not in {"Unassigned_Packera", "Unknown", ""} and pd.notna(s)]
    target_sp = valid_sp[: int(max_taxa)] if (max_taxa is not None and int(max_taxa) > 0) else valid_sp
    style_map = get_combinatorial_taxon_style(target_sp)

    manifest: Dict[str, Dict[str, Any]] = {}
    for sp_name in target_sp:
        sp_df = df_zone[df_zone["species_clean"] == sp_name]
        if len(sp_df) == 0:
            continue

        # Filter to QC-passing, clean closed/unclipped masks when available
        if "efa_qc_pass" in sp_df.columns:
            eqc_flags = [_parse_bool_val(v, default=True) for v in sp_df["efa_qc_pass"]]
            if any(eqc_flags):
                sp_df = sp_df[eqc_flags]
        valid_mask_flags = [is_valid_mask_fn(r.get("mask_path", ""), fig_dir) for _, r in sp_df.iterrows()]
        sp_clean_df = sp_df[valid_mask_flags] if any(valid_mask_flags) else sp_df

        best_row = None
        if "ECT_PC1" in sp_clean_df.columns and "ECT_PC2" in sp_clean_df.columns:
            pc1 = pd.to_numeric(sp_clean_df["ECT_PC1"], errors="coerce")
            pc2 = pd.to_numeric(sp_clean_df["ECT_PC2"], errors="coerce")
            valid_pcs = pc1.notna() & pc2.notna()
            if valid_pcs.sum() >= 2:
                sub_pc = sp_clean_df[valid_pcs]
                dist_pc = np.sqrt(
                    (pc1[valid_pcs] - float(pc1[valid_pcs].median())) ** 2
                    + (pc2[valid_pcs] - float(pc2[valid_pcs].median())) ** 2
                )
                is_type_s = (
                    sub_pc["is_type_specimen"].apply(_parse_bool_val)
                    if "is_type_specimen" in sub_pc.columns
                    else pd.Series(False, index=sub_pc.index)
                )
                type_hits = sub_pc[is_type_s & (dist_pc <= 1.5 * float(dist_pc.median() + 1e-6))]
                best_row = type_hits.iloc[0] if len(type_hits) > 0 else sub_pc.loc[dist_pc.idxmin()]

        if best_row is None:
            sol_s = pd.to_numeric(sp_clean_df.get("solidity", 0.88), errors="coerce").fillna(0.88)
            ar_s = pd.to_numeric(sp_clean_df.get("aspect_ratio", 2.0), errors="coerce").fillna(2.0)
            d_geom = (sol_s - float(sol_s.median())) ** 2 + (ar_s - float(ar_s.median())) ** 2
            best_row = sp_clean_df.loc[d_geom.idxmin()]

        cat_num = str(best_row.get("catalogNumber", "Unknown"))
        leaf_id = int(best_row.get("leaf_id", 1)) if pd.notna(best_row.get("leaf_id", 1)) else 1
        taxon_slug = re.sub(r"[^A-Za-z0-9_]+", "_", sp_name.lower()).strip("_")
        out_ect_png = overlays_dir / f"{zone_clean}_{taxon_slug}.png"
        out_sil_png = lm3_assets_dir / f"{zone_clean}_{taxon_slug}_silhouette.png"
        out_rgb_png = lm3_assets_dir / f"{zone_clean}_{taxon_slug}_cutout.png"

        mask_u8, rgb_patch = load_mask_fn(best_row, fig_dir)
        wrote_ect = False
        if mask_u8 is not None and int(np.count_nonzero(mask_u8)) > 0:
            ect_rgb = render_ect_fn(mask_u8, n_dirs=n_dirs, size=320)
            if ect_rgb is not None:
                cv2.imwrite(str(out_ect_png), cv2.cvtColor(ect_rgb, cv2.COLOR_RGB2BGR))
                wrote_ect = True

            acc_col = style_map.get(sp_name, ("#0072B2", "o"))[0]
            px_mm = pd.to_numeric(best_row.get("pixels_per_mm", np.nan), errors="coerce")
            px_arg = float(px_mm) if pd.notna(px_mm) else None
            is_cal = _parse_bool_val(best_row.get("is_calibrated_mm", False)) and str(
                best_row.get("scale_source", "")
            ).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated"}
            sil_card = render_sil_fn(mask_u8, None, accent_hex=acc_col, pixels_per_mm=px_arg, is_calibrated=is_cal, card_size=320, prefer_rgb=False)
            cv2.imwrite(str(out_sil_png), cv2.cvtColor(sil_card, cv2.COLOR_RGB2BGR))
            if rgb_patch is not None:
                cut_card = render_sil_fn(mask_u8, rgb_patch, accent_hex=acc_col, pixels_per_mm=px_arg, is_calibrated=is_cal, card_size=320, prefer_rgb=True)
                cv2.imwrite(str(out_rgb_png), cv2.cvtColor(cut_card, cv2.COLOR_RGB2BGR))
        elif out_ect_png.exists() and out_ect_png.stat().st_size > 500:
            wrote_ect = True

        if not wrote_ect:
            continue

        comp_val = best_row.get("ect_lobe_complexity", best_row.get("cauline_ect_lobe_complexity", np.nan))
        comp_float = float(comp_val) if pd.notna(comp_val) else np.nan
        manifest[sp_name] = {
            "image_path": str(out_ect_png),
            "silhouette_path": str(out_sil_png) if out_sil_png.exists() else None,
            "cutout_path": str(out_rgb_png) if out_rgb_png.exists() else None,
            "catalogNumber": cat_num,
            "leaf_id": leaf_id,
            "taxon": sp_name,
            "taxon_slug": taxon_slug,
            "zone": zone_clean,
            "is_type": bool(_parse_bool_val(best_row.get("is_type_specimen", False))),
            "ect_lobe_complexity": round(comp_float, 4) if np.isfinite(comp_float) else None,
            "aspect_ratio": round(float(best_row.get("aspect_ratio", np.nan)), 4) if pd.notna(best_row.get("aspect_ratio")) else None,
            "solidity": round(float(best_row.get("solidity", np.nan)), 4) if pd.notna(best_row.get("solidity")) else None,
        }

    for m_dir in (overlays_dir, lm3_assets_dir):
        with open(m_dir / f"manifest_{zone_clean}.json", "w", encoding="utf-8") as f_m:
            json.dump(manifest, f_m, indent=2)

    logger.info(f"Exported {len(manifest)} authentic LeafMachine3 radial ECT overlays to {overlays_dir}")
    return manifest


export_taxon_lm3_exemplars = export_taxon_ect_exemplars

__all__ = [
    "_render_lm3_silhouette_card", "_render_lm3_ect_rgb",
    "export_taxon_ect_exemplars", "export_taxon_lm3_exemplars",
]
