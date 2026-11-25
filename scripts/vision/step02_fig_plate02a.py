#!/usr/bin/env python3
"""
Step 02 figure cohort preparation, legacy artifact cleanup, atomic saver hook
resolution, quantile exemplar selection, and Plate 02A (`fig02a_lm3_organ_and_ect_exemplars.png`).
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.step02_fig_ect_exemplars import _render_lm3_ect_rgb, _render_lm3_silhouette_card, export_taxon_ect_exemplars
from scripts.vision.step02_fig_style_assets import (
    _atomic_save_matplotlib_fig, _discover_downstream_table, _is_valid_smooth_closed_mask,
    _load_mask_or_contour_u8, _parse_bool_val, _resolve_asset_path, _resolve_facade_hook,
    _standardize_packera_species_py, format_icn_taxon_mathtext, get_combinatorial_taxon_style,
)

logger = logging.getLogger("SegmentAndExtract")

_ORIG_ATOMIC_SAVE, _ORIG_DISCOVER_TABLE = _atomic_save_matplotlib_fig, _discover_downstream_table
_ORIG_IS_VALID_MASK, _ORIG_LOAD_MASK = _is_valid_smooth_closed_mask, _load_mask_or_contour_u8
_ORIG_RESOLVE_ASSET_PATH, _ORIG_EXPORT_ECT = _resolve_asset_path, export_taxon_ect_exemplars
_ORIG_RENDER_ECT, _ORIG_RENDER_SIL = _render_lm3_ect_rgb, _render_lm3_silhouette_card

_LEGACY_STEP02_STEMS: Tuple[str, ...] = (
    "fig02a_segmentation_qc_distributions", "fig02b_silhouette_montage",
    "fig02b_silhouette_montage_supplement", "fig02c_cauline_silhouette_montage",
    "fig02c_cauline_silhouette_montage_supplement", "fig02d_basal_ect_radial_montage",
    "fig02e_cauline_ect_radial_montage",
)


def _cleanup_legacy_step02_artifacts(out_dir: Path) -> None:
    """Remove legacy multi-page or duplicate Step 02 PDF/PNG remnants from prior pipeline versions."""
    for legacy_stem in _LEGACY_STEP02_STEMS:
        for ext in (".pdf", ".png"):
            rem = out_dir / f"{legacy_stem}{ext}"
            if rem.exists():
                rem.unlink(missing_ok=True)


def _resolve_atomic_fig_saver() -> Any:
    """Resolve atomic figure saver, honoring monkeypatches on Step 02 pipeline or figures facades."""
    for mod_name in ("scripts.pipeline.02_segment_and_extract", "scripts.vision.step02_figures", "scripts.vision.step02_fig_style_assets"):
        mod = sys.modules.get(mod_name)
        if mod is not None:
            patched = getattr(mod, "_atomic_save_matplotlib_fig", None)
            if patched is not None and patched is not _ORIG_ATOMIC_SAVE:
                return patched
    return _atomic_save_matplotlib_fig


def _is_exemplar_eligible(
    row: pd.Series,
    out_dir: Path,
    misid_cats: Set[str],
    outlier_cats: Set[str],
    outlier_leaf_keys: Set[Tuple[str, int]],
) -> bool:
    """Check whether a segmented leaf passes downstream outlier, misid, LM3 gate, and closed-mask checks."""
    c_s = str(row.get("cat_clean", "")).strip()
    lid = int(float(row.get("leaf_id", 1))) if pd.notna(row.get("leaf_id", 1)) else 1
    if c_s in misid_cats or c_s in outlier_cats or (c_s, lid) in outlier_leaf_keys:
        return False
    if "lm3_gates_pass" in row and pd.notna(row.get("lm3_gates_pass")):
        if not _parse_bool_val(row.get("lm3_gates_pass", True), default=True):
            return False
    is_valid_fn = _resolve_facade_hook("_is_valid_smooth_closed_mask", _is_valid_smooth_closed_mask, _ORIG_IS_VALID_MASK)
    return is_valid_fn(row.get("mask_path", ""), out_dir)


def _pick_distinct_quantile_rows(
    cand_df: pd.DataFrame,
    quantiles: Sequence[float],
    window: float = 0.10,
) -> List[Optional[pd.Series]]:
    """Select distinct authority- and calibration-prioritized leaves near target shape quantiles."""
    if len(cand_df) == 0:
        return [None] * len(quantiles)
    ranked = cand_df.assign(
        _ar=cand_df["ar_clean"].fillna(cand_df["ar_clean"].median() if cand_df["ar_clean"].notna().any() else 1.8),
        _sol=cand_df["solidity_clean"].fillna(0.85),
    ).sort_values(by=["_ar", "_sol"], ascending=[True, False])
    n_r = len(ranked)
    q_pos = [p / max(1, n_r - 1) if n_r > 1 else 0.5 for p in range(n_r)]
    meta_scores: List[Tuple[int, int, int, float]] = []
    for p in range(n_r):
        r_p = ranked.iloc[p]
        is_t = 1 if _parse_bool_val(r_p.get("is_type_specimen", False)) else 0
        t_str = str(r_p.get("determiner_tier", "")).strip()
        t_rk = 3 if ("Tier_1" in t_str or "Gold" in t_str) else (2 if ("Tier_2" in t_str or "Silver" in t_str) else (1 if ("Tier_3" in t_str or "Bronze" in t_str) else 0))
        is_c = 1 if (
            _parse_bool_val(r_p.get("is_calibrated_mm", False))
            and str(r_p.get("scale_source", "")).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated"}
        ) else 0
        lqw_v = pd.to_numeric(r_p.get("leaf_quality_weight", r_p.get("archetype_score", 0.5)), errors="coerce")
        meta_scores.append((is_t, t_rk, is_c, round(float(lqw_v), 4) if pd.notna(lqw_v) else 0.5))

    used_idx: Set[Any] = set()
    picked: List[Optional[pd.Series]] = []
    for q in quantiles:
        avail = [p for p in range(n_r) if ranked.index[p] not in used_idx]
        if not avail:
            picked.append(None)
            continue
        win = [p for p in avail if abs(q_pos[p] - q) <= window]
        if win:
            chosen_pos = max(win, key=lambda p: (*meta_scores[p], -abs(q_pos[p] - q)))
        else:
            chosen_pos = min(
                avail,
                key=lambda p: (round(abs(q_pos[p] - q), 6), -meta_scores[p][0], -meta_scores[p][1], -meta_scores[p][2], -meta_scores[p][3]),
            )
        used_idx.add(ranked.index[chosen_pos])
        picked.append(ranked.iloc[chosen_pos])
    return picked


def prepare_step02_figure_context(
    extracted_df: pd.DataFrame,
    out_dir: Path,
    vouchers_df: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    """Clean stale legacy artifacts, enrich Step 02 manifests with voucher/QC flags, and build shared plate context."""
    discover_fn = _resolve_facade_hook("_discover_downstream_table", _discover_downstream_table, _ORIG_DISCOVER_TABLE)
    _cleanup_legacy_step02_artifacts(out_dir)

    df = extracted_df.copy()
    if "catalogNumber" not in df.columns:
        df["catalogNumber"] = [f"SPEC_{i:04d}" for i in range(len(df))]
    df["cat_clean"] = df["catalogNumber"].astype(str).str.strip()
    active_cats: Set[str] = set(df["cat_clean"].unique()) - {"", "nan", "None"}

    # Enrich taxon, infraspecific rank, ruler scale, and scan paths from curated_vouchers.csv if available
    if vouchers_df is None:
        v_cand = discover_fn("curated_vouchers.csv", out_dir, active_cats)
        if v_cand is not None:
            try:
                vouchers_df = pd.read_csv(v_cand, low_memory=False)
            except Exception:
                vouchers_df = None

    if vouchers_df is not None and len(vouchers_df) > 0 and "catalogNumber" in vouchers_df.columns:
        v_sub = vouchers_df.drop_duplicates(subset=["catalogNumber"]).copy()
        v_sub["cat_clean"] = v_sub["catalogNumber"].astype(str).str.strip()
        # Restrict active species strictly to vouchers present in curated_vouchers.csv when matched
        matched_vcats = set(v_sub["cat_clean"]) & active_cats
        if matched_vcats:
            df = df[df["cat_clean"].isin(matched_vcats)].copy()
            active_cats = matched_vcats
        for v_col in ("species_standardized", "scientificName", "infraspecific_taxon", "determiner_tier", "is_type_specimen", "type_status", "image_path"):
            if v_col in v_sub.columns and v_col not in df.columns:
                df[v_col] = df["cat_clean"].map(dict(zip(v_sub["cat_clean"], v_sub[v_col])))

    sp_col = next((c for c in ("species_standardized", "scientificName", "species_name", "species_raw", "taxon") if c in df.columns), None)
    inf_col = "infraspecific_taxon" if "infraspecific_taxon" in df.columns else None
    if sp_col is not None:
        df["species_clean"] = (
            [_standardize_packera_species_py(s, inf) for s, inf in zip(df[sp_col], df[inf_col])]
            if inf_col is not None
            else df[sp_col].apply(_standardize_packera_species_py)
        )
    else:
        df["species_clean"] = "Unassigned_Packera"

    df["solidity_clean"] = pd.to_numeric(df.get("solidity", df.get("sol_val", np.nan)), errors="coerce")
    df["ar_clean"] = pd.to_numeric(df.get("aspect_ratio", df.get("ar_val", np.nan)), errors="coerce")
    df["zone_clean"] = df["leaf_zone"].astype(str).str.strip().str.lower() if "leaf_zone" in df.columns else pd.Series("basal", index=df.index)
    df["zone_clean"] = df["zone_clean"].where(df["zone_clean"].isin({"basal", "cauline"}), "basal")

    # Discover downstream Step 03 contour outlier flags and Step 04 misidentification flags so exemplars are 100% clean
    outlier_leaf_keys: Set[Tuple[str, int]] = set()
    outlier_cats: Set[str] = set()
    misid_cats: Set[str] = set()

    ind_efa_p = discover_fn("leaf_efa_harmonics_individual.csv", out_dir, active_cats)
    if ind_efa_p is not None:
        try:
            ind_df = pd.read_csv(ind_efa_p, low_memory=False)
            if {"catalogNumber", "is_contour_outlier"}.issubset(ind_df.columns):
                for _, r_i in ind_df.iterrows():
                    if _parse_bool_val(r_i.get("is_contour_outlier", False)):
                        c_s = str(r_i["catalogNumber"]).strip()
                        if "leaf_id" in ind_df.columns and pd.notna(r_i.get("leaf_id")):
                            outlier_leaf_keys.add((c_s, int(float(r_i["leaf_id"]))))
                        else:
                            outlier_cats.add(c_s)
            if "leaf_quality_weight" in ind_df.columns and "leaf_quality_weight" not in df.columns and "catalogNumber" in ind_df.columns:
                lqw_map: Dict[Tuple[str, int], float] = {}
                for _, r_w in ind_df.iterrows():
                    c_w = str(r_w.get("catalogNumber", "")).strip()
                    lid_w = int(float(r_w.get("leaf_id", 1))) if pd.notna(r_w.get("leaf_id", 1)) else 1
                    w_v = pd.to_numeric(r_w.get("leaf_quality_weight", None), errors="coerce")
                    if c_w and pd.notna(w_v):
                        lqw_map[(c_w, lid_w)] = float(w_v)
                        lqw_map.setdefault((c_w, -1), float(w_v))
                if lqw_map:
                    df["leaf_quality_weight"] = [
                        lqw_map.get(
                            (str(r.get("cat_clean", "")).strip(), int(float(r.get("leaf_id", 1))) if pd.notna(r.get("leaf_id", 1)) else 1),
                            lqw_map.get((str(r.get("cat_clean", "")).strip(), -1), float("nan")),
                        )
                        for _, r in df.iterrows()
                    ]
        except Exception:
            pass

    otu_efa_p = discover_fn("leaf_efa_harmonics.csv", out_dir, active_cats)
    if otu_efa_p is not None:
        try:
            otu_df = pd.read_csv(otu_efa_p, low_memory=False)
            if {"catalogNumber", "is_contour_outlier"}.issubset(otu_df.columns):
                for _, r_o in otu_df.iterrows():
                    if _parse_bool_val(r_o.get("is_contour_outlier", False)):
                        outlier_cats.add(str(r_o["catalogNumber"]).strip())
        except Exception:
            pass

    mflags_p = discover_fn("morphometric_flags.csv", out_dir, active_cats)
    if mflags_p is not None:
        try:
            mf_df = pd.read_csv(mflags_p, low_memory=False)
            if "catalogNumber" in mf_df.columns:
                for _, r_m in mf_df.iterrows():
                    if _parse_bool_val(r_m.get("misidentification_flag", False)) or _parse_bool_val(r_m.get("is_confident_label_noise", False)):
                        misid_cats.add(str(r_m["catalogNumber"]).strip())
        except Exception:
            pass

    df["exemplar_ok"] = [_is_exemplar_eligible(r, out_dir, misid_cats, outlier_cats, outlier_leaf_keys) for _, r in df.iterrows()]
    sp_counts = df["species_clean"].value_counts()
    valid_taxa = [s for s in sp_counts.index if s not in {"Unassigned_Packera", "Unknown", ""} and pd.notna(s)]
    if not valid_taxa:
        valid_taxa = [str(s) for s in sp_counts.index if str(s).strip()]
    active_taxa = sorted(valid_taxa)

    return {
        "df": df, "vouchers_df": vouchers_df, "active_taxa": active_taxa,
        "style_map": get_combinatorial_taxon_style(active_taxa),
        "has_cauline_any": bool((df["zone_clean"] == "cauline").any()),
        "misid_cats": misid_cats, "outlier_cats": outlier_cats, "outlier_leaf_keys": outlier_leaf_keys,
    }


def _filter_high_conf_zone(pool_df: pd.DataFrame) -> pd.DataFrame:
    """Prefer high-confidence zone assignments (>= 0.72) when at least two candidate leaves remain."""
    if len(pool_df) >= 2 and "leaf_zone_confidence" in pool_df.columns:
        hi = pool_df[pd.to_numeric(pool_df["leaf_zone_confidence"], errors="coerce").fillna(1.0) >= 0.72]
        if len(hi) >= 2:
            return hi
    return pool_df


def render_plate02a_and_ect_cache(
    ctx: Dict[str, Any],
    out_dir: Path,
    save_fig_fn: Optional[Any] = None,
) -> Path:
    """Export native LM3/ECT exemplar caches and render Plate 02A (`fig02a_lm3_organ_and_ect_exemplars.png`)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    df: pd.DataFrame = ctx["df"]
    active_taxa: List[str] = ctx["active_taxa"]
    style_map: Dict[str, Tuple[str, str]] = ctx["style_map"]
    has_cauline_any: bool = ctx["has_cauline_any"]
    misid_cats: Set[str] = ctx["misid_cats"]
    outlier_cats: Set[str] = ctx["outlier_cats"]

    export_ect_fn = _resolve_facade_hook("export_taxon_ect_exemplars", export_taxon_ect_exemplars, _ORIG_EXPORT_ECT)
    resolve_fn = _resolve_facade_hook("_resolve_asset_path", _resolve_asset_path, _ORIG_RESOLVE_ASSET_PATH)
    load_mask_fn = _resolve_facade_hook("_load_mask_or_contour_u8", _load_mask_or_contour_u8, _ORIG_LOAD_MASK)
    render_ect_fn = _resolve_facade_hook("_render_lm3_ect_rgb", _render_lm3_ect_rgb, _ORIG_RENDER_ECT)
    render_sil_fn = _resolve_facade_hook("_render_lm3_silhouette_card", _render_lm3_silhouette_card, _ORIG_RENDER_SIL)
    saver = save_fig_fn if save_fig_fn is not None else _resolve_atomic_fig_saver()

    # 1. Export Native LeafMachine3 Asset Cache (`ect_overlays/` & `lm3_assets/`)
    clean_for_ect = df[df["exemplar_ok"]].copy() if df["exemplar_ok"].any() else df
    try:
        export_ect_fn(extracted_df=clean_for_ect, figures_dir=out_dir, zone="basal", max_taxa=None)
        if has_cauline_any:
            export_ect_fn(extracted_df=clean_for_ect, figures_dir=out_dir, zone="cauline", max_taxa=None)
    except Exception as exc_ect:
        logger.warning(f"Could not export taxon LM3/ECT exemplars: {exc_ect}")

    # 2. Plate 02A: dynamic rows = ALL active taxa/ranks (zero [:8] truncation)
    n_rows_a = max(1, len(active_taxa))
    n_cols_a = 5 if has_cauline_any else 4
    fig_w_a = 14.6 if has_cauline_any else 12.4
    fig_h_a = max(3.4, min(95.0, 1.72 * n_rows_a + 0.85))
    fig_a = plt.figure(figsize=(fig_w_a, fig_h_a), dpi=300, layout="constrained")
    gs_a = GridSpec(n_rows_a, n_cols_a, figure=fig_a, wspace=0.08, hspace=0.18)

    col_headers = (
        ["20th %ile (Compact)", "50th %ile (Medoid)", "80th %ile (Dissected)", "LM3 Radial ECT"]
        if not has_cauline_any
        else ["Basal Medoid (LM3 Cutout)", "Basal 80th %ile (Dissected)", "Basal Radial ECT (LM3)", "Cauline Medoid (LM3)", "Cauline Radial ECT (LM3)"]
    )

    for r_idx, sp_name in enumerate(active_taxa):
        sp_all = df[df["species_clean"] == sp_name]
        sp_clean = sp_all[sp_all["exemplar_ok"]]
        if len(sp_clean) == 0:
            n_missing = sum(1 for _, r_chk in sp_all.iterrows() if resolve_fn(r_chk.get("mask_path", ""), out_dir) is None)
            if n_missing > 0:
                logger.warning(f"fig02a: {n_missing} / {len(sp_all)} {sp_name} extracted leaves have no resolvable mask/contour asset on disk.")
        col_hex, _ = style_map[sp_name]

        basal_pool = sp_clean[sp_clean["zone_clean"] == "basal"]
        if len(basal_pool) == 0 and len(sp_all[sp_all["zone_clean"] == "basal"]) > 0:
            # Only fall back if masks were synthetic/empty path (never admit flagged outliers/clipped masks)
            sp_basal_only = sp_all[sp_all["zone_clean"] == "basal"]
            empty_mask_flags = pd.Series(
                [not bool(str(pv or "").strip()) for pv in sp_basal_only.get("mask_path", pd.Series("", index=sp_basal_only.index))],
                index=sp_basal_only.index, dtype=bool,
            )
            basal_pool = sp_basal_only[(~sp_basal_only["cat_clean"].isin(misid_cats | outlier_cats)) & empty_mask_flags]
        basal_pool = _filter_high_conf_zone(basal_pool)
        cauline_pool = _filter_high_conf_zone(sp_clean[sp_clean["zone_clean"] == "cauline"])

        if not has_cauline_any:
            q_rows = _pick_distinct_quantile_rows(basal_pool, [0.24, 0.50, 0.80])
            med_for_ect = q_rows[1] if q_rows[1] is not None else (q_rows[0] if q_rows[0] is not None else q_rows[2])
            tile_specs = [
                ("leaf", q_rows[0], False, "20th %ile (Compact)"), ("leaf", q_rows[1], True, "50th %ile (Medoid)"),
                ("leaf", q_rows[2], False, "80th %ile (Dissected)"), ("ect", med_for_ect, False, "LM3 Radial ECT"),
            ]
        else:
            b_rows = _pick_distinct_quantile_rows(basal_pool, [0.50, 0.80])
            c_rows = _pick_distinct_quantile_rows(cauline_pool, [0.50])
            b_med, c_med = b_rows[0], c_rows[0]
            tile_specs = [
                ("leaf", b_med, True, "Basal Medoid"), ("leaf", b_rows[1], False, "Basal 80th %ile"),
                ("ect", b_med, False, "Basal Radial ECT"), ("leaf", c_med, True, "Cauline Medoid"), ("ect", c_med, False, "Cauline Radial ECT"),
            ]

        for c_idx, (kind, row_sel, prefer_rgb, sub_tag) in enumerate(tile_specs):
            ax = fig_a.add_subplot(gs_a[r_idx, c_idx])
            ax.set_xticks([])
            ax.set_yticks([])
            for sp_edge in ax.spines.values():
                sp_edge.set_edgecolor("#e5e7eb")
                sp_edge.set_linewidth(0.8)

            if c_idx == 0:
                wrapped_sp = format_icn_taxon_mathtext(sp_name, abbreviate_genus=True)
                wrapped_sp = re.sub(r"\s+(subsp\.|var\.|f\.)\s+", r"\n\1 ", wrapped_sp)
                ax.set_ylabel(f"{wrapped_sp}\n(n = {len(sp_all)} leaves)", fontsize=8.2, color="#111827", rotation=0, ha="right", va="center", labelpad=12)

            top_tag = col_headers[c_idx] if r_idx == 0 else sub_tag
            if row_sel is None:
                ax.set_facecolor("#f9fafb")
                ax.text(0.5, 0.5, "No clean organ\nin cohort", ha="center", va="center", fontsize=7.5, color="#9ca3af", transform=ax.transAxes)
                ax.set_title(f"{top_tag}\nN/A", fontsize=7.4, color="#6b7280", pad=3)
                continue

            cat_lbl = str(row_sel.get("cat_clean", row_sel.get("catalogNumber", "")))
            lid_lbl = int(float(row_sel.get("leaf_id", 1))) if pd.notna(row_sel.get("leaf_id", 1)) else 1
            mask_u8, rgb_patch = load_mask_fn(row_sel, out_dir)
            if mask_u8 is None or int(np.count_nonzero(mask_u8)) == 0:
                ax.set_facecolor("#f9fafb")
                ax.text(0.5, 0.5, f"{cat_lbl}\n(mask unavailable)", ha="center", va="center", fontsize=7.5, color="#9ca3af", transform=ax.transAxes)
                ax.set_title(f"{top_tag}\n{cat_lbl} (L{lid_lbl})", fontsize=7.4, color="#374151", pad=3)
                continue

            if kind == "ect":
                ect_card = render_ect_fn(mask_u8, n_dirs=160, size=320)
                if ect_card is not None:
                    ax.imshow(ect_card)
                else:
                    ax.text(0.5, 0.5, "ECT unavailable", ha="center", va="center", fontsize=7.5, transform=ax.transAxes)
            else:
                px_mm = pd.to_numeric(row_sel.get("pixels_per_mm", np.nan), errors="coerce")
                is_cal = _parse_bool_val(row_sel.get("is_calibrated_mm", False)) and str(row_sel.get("scale_source", "")).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated", "tick_lattice_unboxed"}
                is_type_row = _parse_bool_val(row_sel.get("is_type_specimen", False)) or (
                    str(row_sel.get("type_status", "")).strip().lower() not in {"", "none", "nan", "na", "not a type"}
                )
                try:
                    leaf_card = render_sil_fn(
                        mask_u8, rgb_patch, accent_hex=col_hex, pixels_per_mm=float(px_mm) if pd.notna(px_mm) else None,
                        is_calibrated=is_cal, card_size=320, prefer_rgb=prefer_rgb, is_type=is_type_row,
                    )
                except TypeError:
                    leaf_card = render_sil_fn(
                        mask_u8, rgb_patch, accent_hex=col_hex, pixels_per_mm=float(px_mm) if pd.notna(px_mm) else None,
                        is_calibrated=is_cal, card_size=320, prefer_rgb=prefer_rgb,
                    )
                ax.imshow(leaf_card)

            ar_v = pd.to_numeric(row_sel.get("ar_clean", np.nan), errors="coerce")
            sol_v = pd.to_numeric(row_sel.get("solidity_clean", np.nan), errors="coerce")
            stat_s = f"L/W={float(ar_v):.2f}, Sol={float(sol_v):.2f}" if (pd.notna(ar_v) and pd.notna(sol_v)) else f"L{lid_lbl}"
            ax.set_title(f"{top_tag}\n{cat_lbl} ({stat_s})", fontsize=7.3, color="#1f2937", pad=3)

    fig_a.suptitle("Step 02: Native LeafMachine3 Organ Silhouettes, Pose Landmarks & Radial ECT Exemplars", fontsize=11.5, fontweight="bold")
    fig02a_png = saver(fig_a, out_dir / "fig02a_lm3_organ_and_ect_exemplars", dpi=300)
    plt.close(fig_a)
    return fig02a_png[1] if isinstance(fig02a_png, tuple) else Path(fig02a_png)


__all__ = [
    "_LEGACY_STEP02_STEMS", "_cleanup_legacy_step02_artifacts", "_resolve_atomic_fig_saver",
    "_is_exemplar_eligible", "_pick_distinct_quantile_rows",
    "prepare_step02_figure_context", "render_plate02a_and_ect_cache",
]
