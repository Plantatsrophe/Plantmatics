#!/usr/bin/env python3
"""
Step 02 Plate 02B (`fig02b_segmentation_qc_and_phenophase.png` Panels A-E),
capped whole-sheet overlay dispatch, and top-level Step 02 figure orchestration.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.step02_fig_plate02a import (
    _resolve_atomic_fig_saver,
    prepare_step02_figure_context,
    render_plate02a_and_ect_cache,
)
from scripts.vision.step02_fig_sheet_overlays import (
    _draw_solidity_vs_ar_zone_panel,
    export_lm3_sheet_diagnostic_overlays,
)
from scripts.vision.step02_fig_style_assets import (
    _parse_bool_val,
    _resolve_facade_hook,
    format_icn_taxon_mathtext,
)

logger = logging.getLogger("SegmentAndExtract")
_ORIG_DRAW_ZONE_PANEL = _draw_solidity_vs_ar_zone_panel
_ORIG_EXPORT_SHEET_OVERLAYS = export_lm3_sheet_diagnostic_overlays


def render_plate02b_qc_and_phenophase(
    ctx: Dict[str, Any],
    failed_df: pd.DataFrame,
    repro_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    out_dir: Path,
    save_fig_fn: Optional[Any] = None,
) -> Path:
    """Render Plate 02B (`fig02b_segmentation_qc_and_phenophase.png`) with Panels A-E at 300 DPI."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    df: pd.DataFrame = ctx["df"]
    active_taxa: List[str] = ctx["active_taxa"]
    style_map: Dict[str, Tuple[str, str]] = ctx["style_map"]
    has_cauline_any: bool = ctx["has_cauline_any"]

    draw_zone_fn = _resolve_facade_hook("_draw_solidity_vs_ar_zone_panel", _draw_solidity_vs_ar_zone_panel, _ORIG_DRAW_ZONE_PANEL)
    saver = save_fig_fn if save_fig_fn is not None else _resolve_atomic_fig_saver()

    # Panel A (A1 Basal / A2 Cauline hollow scatter), Panel B (Scale & QC),
    # Panel C (Involucral AR_inv hollow-jitter boxplots), Panel D (Phenophases), optional Panel E
    extra_h = max(0.0, 0.24 * max(0, len(active_taxa) - 8))
    if has_cauline_any:
        fig_b = plt.figure(figsize=(14.4, 13.2 + 1.5 * extra_h), dpi=300, layout="constrained")
        gs_b = GridSpec(3, 2, figure=fig_b, height_ratios=[1.0, 1.0 + extra_h / 5.5, 0.95 + extra_h / 5.5])
        draw_zone_fn(
            fig_b.add_subplot(gs_b[0, 0]), df[df["zone_clean"] == "basal"], active_taxa, style_map,
            "A1. Basal Rosette Homology Space (Hollow Markers)", show_legend=True,
        )
        draw_zone_fn(
            fig_b.add_subplot(gs_b[0, 1]), df[df["zone_clean"] == "cauline"], active_taxa, style_map,
            "A2. Cauline Stem Homology Space (Hollow Markers)", show_legend=False,
        )
        ax_qc, ax_ar_inv = fig_b.add_subplot(gs_b[1, 0]), fig_b.add_subplot(gs_b[1, 1])
        ax_pheno, ax_het = fig_b.add_subplot(gs_b[2, 0]), fig_b.add_subplot(gs_b[2, 1])
    else:
        fig_b = plt.figure(figsize=(14.2, 9.6 + extra_h), dpi=300, layout="constrained")
        gs_b = GridSpec(2, 2, figure=fig_b, height_ratios=[1.0, 1.0 + extra_h / 5.2])
        draw_zone_fn(
            fig_b.add_subplot(gs_b[0, 0]), df, active_taxa, style_map,
            "A. Basal Rosette Solidity vs. Blade Aspect Ratio (L / W)", show_legend=True,
        )
        ax_qc, ax_ar_inv = fig_b.add_subplot(gs_b[0, 1]), fig_b.add_subplot(gs_b[1, 0])
        ax_pheno, ax_het = fig_b.add_subplot(gs_b[1, 1]), None

    y_pos = np.arange(len(active_taxa))
    pct_ruler, pct_mp, pct_uncal, qc_pass_pct, counts_lbl = [], [], [], [], []
    for tx in active_taxa:
        sub_tx = df[df["species_clean"] == tx]
        n_t = max(1, len(sub_tx))
        sc_raw = sub_tx.get("scale_source", pd.Series("", index=sub_tx.index)).astype(str).str.strip().str.lower()
        cal_b = sub_tx["is_calibrated_mm"].apply(_parse_bool_val) if "is_calibrated_mm" in sub_tx.columns else pd.Series(False, index=sub_tx.index)
        is_mp = sc_raw.isin({"mp_estimated", "predicted_from_megapixels"})
        is_rul = (cal_b & (~is_mp)) | sc_raw.isin({"measured_from_ruler", "measured_from_fieldprism", "ruler_calibrated"})
        is_unc = ~(is_rul | is_mp)
        pct_ruler.append(100.0 * float(is_rul.sum()) / n_t)
        pct_mp.append(100.0 * float(is_mp.sum()) / n_t)
        pct_uncal.append(100.0 * float(is_unc.sum()) / n_t)
        qc_pass_pct.append(100.0 * float(sub_tx["exemplar_ok"].sum()) / n_t)
        counts_lbl.append(f"{format_icn_taxon_mathtext(tx, abbreviate_genus=(len(active_taxa) > 4))} (n={len(sub_tx)})")

    b1, b2, b3 = np.array(pct_ruler), np.array(pct_mp), np.array(pct_uncal)
    ax_qc.barh(y_pos, b1, height=0.62, color="#0072B2", alpha=0.85, label="Ruler / FieldPrism (mm)")
    ax_qc.barh(y_pos, b2, left=b1, height=0.62, color="#E69F00", alpha=0.85, label="Megapixel Prior")
    ax_qc.barh(y_pos, b3, left=b1 + b2, height=0.62, color="#d1d5db", alpha=0.90, label="Uncalibrated (px)")
    ax_qc.scatter(
        qc_pass_pct, y_pos, s=42, marker="D", facecolors="none", edgecolors="#111827",
        linewidths=1.2, zorder=5, label="Clean EFA/LM3 QC %",
    )
    ax_qc.set_yticks(y_pos)
    ax_qc.set_yticklabels(counts_lbl, fontsize=7.8)
    ax_qc.invert_yaxis()
    ax_qc.set_xlim(0, 108)
    ax_qc.set_xlabel("Percentage of Extracted Leaves (%)", fontsize=9)
    n_fail_sheets = len(failed_df) if failed_df is not None else 0
    ax_qc.set_title(f"B. Ruler Calibration & Clean LM3 Gate Yield (Failed: {n_fail_sheets})", fontsize=10, fontweight="bold", loc="left", pad=36)
    ax_qc.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=7.0, frameon=False)
    ax_qc.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)

    # Panel C: Involucral Aspect Ratio (AR_inv) hollow-jitter boxplots
    rep_work = pd.DataFrame()
    if repro_df is not None and len(repro_df) > 0 and "AR_inv" in repro_df.columns:
        rep_work = repro_df.copy()
        if "species_clean" not in rep_work.columns and "catalogNumber" in rep_work.columns:
            rep_work["species_clean"] = rep_work["catalogNumber"].astype(str).str.strip().map(dict(zip(df["cat_clean"], df["species_clean"])))
        rep_work["AR_inv_num"] = pd.to_numeric(rep_work["AR_inv"], errors="coerce")
        rep_work = rep_work[rep_work["AR_inv_num"].notna() & rep_work["species_clean"].isin(active_taxa)]

    if len(rep_work) > 0:
        taxa_with_rep = [tx for tx in active_taxa if (rep_work["species_clean"] == tx).any()]
        box_data = [rep_work.loc[rep_work["species_clean"] == tx, "AR_inv_num"].to_numpy() for tx in taxa_with_rep]
        pos_r = np.arange(len(taxa_with_rep))
        bp = ax_ar_inv.boxplot(
            box_data, positions=pos_r, orientation="horizontal", widths=0.52,
            patch_artist=True, showfliers=False, medianprops=dict(color="#111827", linewidth=1.4),
        )
        for patch_b in bp["boxes"]:
            patch_b.set_facecolor("none")
            patch_b.set_edgecolor("#4b5563")
            patch_b.set_linewidth(0.95)
        rng_j = np.random.default_rng(42)
        for idx_r, tx in enumerate(taxa_with_rep):
            vals_r = box_data[idx_r]
            col_r, mk_r = style_map[tx]
            jit = rng_j.uniform(-0.16, 0.16, size=len(vals_r))
            ax_ar_inv.scatter(
                vals_r, np.full(len(vals_r), idx_r, dtype=float) + jit, s=28, marker=mk_r,
                facecolors="none", edgecolors=col_r, linewidths=0.95, alpha=0.82, zorder=4,
            )
        ax_ar_inv.axvline(1.0, color="#9ca3af", linestyle=":", linewidth=0.8)
        ax_ar_inv.set_yticks(pos_r)
        ax_ar_inv.set_yticklabels(
            [f"{format_icn_taxon_mathtext(tx, abbreviate_genus=True)} (n={len(box_data[i])})" for i, tx in enumerate(taxa_with_rep)],
            fontsize=7.6,
        )
        ax_ar_inv.invert_yaxis()
        ax_ar_inv.set_xlabel("Involucral Aspect Ratio (Height / Width)", fontsize=8.8)
        ax_ar_inv.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)
    else:
        ax_ar_inv.text(0.5, 0.5, "No measurable capitula\nin cohort", ha="center", va="center", fontsize=8.5, color="#6b7280", transform=ax_ar_inv.transAxes)
        ax_ar_inv.set_xticks([])
        ax_ar_inv.set_yticks([])
    ax_ar_inv.set_title("C. Capitulum Involucre Aspect Ratio", fontsize=10, fontweight="bold", loc="left")

    # Panel D: Determinate Cyme Phenophase Proportions
    ph_work = pd.DataFrame()
    if pheno_df is not None and len(pheno_df) > 0:
        ph_work = pheno_df.copy()
        st_col = next((c for c in ("voucher_phenological_state", "phenological_state", "pheno_state") if c in ph_work.columns), None)
        if st_col is not None:
            if "species_clean" not in ph_work.columns and "catalogNumber" in ph_work.columns:
                ph_work["species_clean"] = ph_work["catalogNumber"].astype(str).str.strip().map(dict(zip(df["cat_clean"], df["species_clean"])))
            ph_work["state_norm"] = ph_work[st_col].astype(str).str.strip().str.lower()
            ph_work = ph_work[ph_work["species_clean"].isin(active_taxa)]
        else:
            ph_work = pd.DataFrame()

    if len(ph_work) > 0:
        taxa_ph = [tx for tx in active_taxa if (ph_work["species_clean"] == tx).any()]
        pos_p = np.arange(len(taxa_ph))
        cats_ph = [
            ("Anthesis", lambda s: s.str.contains("anthesis|flower|bloom", regex=True), "#E69F00"),
            ("Bud", lambda s: s.str.contains("bud", regex=True) & ~s.str.contains("anthesis", regex=True), "#56B4E9"),
            ("Fruiting", lambda s: s.str.contains("fruit|pappus", regex=True) & ~s.str.contains("anthesis", regex=True), "#009E73"),
            ("Sterile / Senescent", lambda s: ~(s.str.contains("anthesis|flower|bloom|bud|fruit|pappus", regex=True)), "#9ca3af"),
        ]
        left_acc = np.zeros(len(taxa_ph), dtype=float)
        for lbl_p, pred_fn, col_p in cats_ph:
            vals_p = [
                100.0 * float(pred_fn(s_tx).sum()) / max(1, len(s_tx))
                for s_tx in (ph_work.loc[ph_work["species_clean"] == tx, "state_norm"] for tx in taxa_ph)
            ]
            arr_p = np.array(vals_p, dtype=float)
            ax_pheno.barh(pos_p, arr_p, left=left_acc, height=0.60, color=col_p, alpha=0.88, label=lbl_p)
            left_acc += arr_p
        ax_pheno.set_yticks(pos_p)
        ax_pheno.set_yticklabels([format_icn_taxon_mathtext(tx, abbreviate_genus=True) for tx in taxa_ph], fontsize=7.6)
        ax_pheno.invert_yaxis()
        ax_pheno.set_xlim(0, 100)
        ax_pheno.set_xlabel("Voucher Proportion (%)", fontsize=8.8)
        ax_pheno.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=6.9, frameon=False)
        ax_pheno.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)
    else:
        ax_pheno.text(0.5, 0.5, "No capitulum phenology\nrecords in cohort", ha="center", va="center", fontsize=8.5, color="#6b7280", transform=ax_pheno.transAxes)
        ax_pheno.set_xticks([])
        ax_pheno.set_yticks([])
    ax_pheno.set_title("D. Cyme Phenophase Stratification", fontsize=10, fontweight="bold", loc="left", pad=36)

    if ax_het is not None:
        # Panel E: Basal vs. Cauline Organ Counts & Heteroblastic Solidity Shift per Taxon
        b_cnts = [int(((df["species_clean"] == tx) & (df["zone_clean"] == "basal")).sum()) for tx in active_taxa]
        c_cnts = [int(((df["species_clean"] == tx) & (df["zone_clean"] == "cauline")).sum()) for tx in active_taxa]
        pos_h = np.arange(len(active_taxa))
        ax_het.barh(pos_h - 0.17, b_cnts, height=0.32, color="#0072B2", alpha=0.86, label="Basal Rosette Leaves")
        ax_het.barh(pos_h + 0.17, c_cnts, height=0.32, color="#009E73", alpha=0.86, label="Cauline Stem Leaves")
        ax_het.set_yticks(pos_h)
        ax_het.set_yticklabels([format_icn_taxon_mathtext(tx, abbreviate_genus=True) for tx in active_taxa], fontsize=7.6)
        ax_het.invert_yaxis()
        ax_het.set_xlabel("Extracted Organ Count (n)", fontsize=8.8)
        ax_het.set_title("E. Basal Rosette vs. Cauline Stem Yield by Taxon", fontsize=10, fontweight="bold", loc="left")
        ax_het.legend(loc="lower right", fontsize=7.4, frameon=True, facecolor="white", framealpha=0.9)
        ax_het.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)

    fig_b.suptitle("Step 02: Segmentation Quality Control, Ruler Calibration & Reproductive Traits", fontsize=11.5, fontweight="bold")
    fig02b_png = saver(fig_b, out_dir / "fig02b_segmentation_qc_and_phenophase", dpi=300)
    plt.close(fig_b)
    return fig02b_png[1] if isinstance(fig02b_png, tuple) else Path(fig02b_png)


def generate_step02_publication_figures(
    extracted_df: pd.DataFrame,
    failed_df: pd.DataFrame,
    repro_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    figures_dir: Union[str, Path],
    vouchers_df: Optional[pd.DataFrame] = None,
) -> Dict[str, Path]:
    """Generate the two canonical Step 02 300-DPI PNG plates + native LM3 asset and sheet-overlay caches."""
    out_dir = Path(figures_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved_paths: Dict[str, Path] = {}
    if extracted_df is None or len(extracted_df) == 0:
        return saved_paths

    save_fig_fn = _resolve_atomic_fig_saver()
    ctx = prepare_step02_figure_context(extracted_df, out_dir, vouchers_df=vouchers_df)
    saved_paths["fig02a_png"] = render_plate02a_and_ect_cache(ctx, out_dir, save_fig_fn=save_fig_fn)
    saved_paths["fig02b_png"] = render_plate02b_qc_and_phenophase(ctx, failed_df, repro_df, pheno_df, out_dir, save_fig_fn=save_fig_fn)

    # Capped Whole-Sheet LeafMachine3 Diagnostic Overlays (`lm3_sheet_overlays/`)
    export_sheets_fn = _resolve_facade_hook("export_lm3_sheet_diagnostic_overlays", export_lm3_sheet_diagnostic_overlays, _ORIG_EXPORT_SHEET_OVERLAYS)
    try:
        overlays_written = export_sheets_fn(
            extracted_df=ctx["df"], repro_df=repro_df, pheno_df=pheno_df,
            figures_dir=out_dir, vouchers_df=ctx["vouchers_df"], max_sheets=24,
        )
        if overlays_written:
            saved_paths["sheet_overlays_dir"] = out_dir / "lm3_sheet_overlays"
    except Exception as exc_sh:
        logger.warning(f"Whole-sheet LeafMachine3 diagnostic overlay export note: {exc_sh}")

    logger.info(f"Step 02 publication figures exported (300-DPI PNG) to: {out_dir}")
    return saved_paths


_ORIG_GEN_STEP02_FIGS = generate_step02_publication_figures
generate_step02_publication_figures_impl = generate_step02_publication_figures


def regenerate_step02_figures_from_tables(
    output_dir: Union[str, Path],
    vouchers_csv: Optional[Union[str, Path]] = None,
) -> Dict[str, Path]:
    """Reload Step 02 + downstream QC tables and regenerate 300-DPI PNG Step 02 plates."""
    out_data_dir = Path(output_dir)
    tables_dir = out_data_dir / "tables" if (out_data_dir / "tables").exists() else out_data_dir
    figures_dir = out_data_dir.parent / "outputs" / "figures"

    def _safe_read_csv(p: Path) -> pd.DataFrame:
        if p.exists() and p.stat().st_size > 0:
            try:
                return pd.read_csv(p, low_memory=False)
            except Exception as e_csv:
                logger.warning(f"Failed to read {p}: {e_csv}")
        return pd.DataFrame()

    extracted_df = _safe_read_csv(tables_dir / "extracted_leaf_manifest.csv")
    failed_df = _safe_read_csv(tables_dir / "failed_qc_vouchers.csv")
    repro_df = _safe_read_csv(tables_dir / "voucher_reproductive_metrics.csv")
    pheno_df = _safe_read_csv(tables_dir / "voucher_phenological_states.csv")
    vouchers_df: Optional[pd.DataFrame] = None

    if vouchers_csv is not None and Path(vouchers_csv).exists() and len(extracted_df) > 0:
        try:
            vouchers_df = pd.read_csv(vouchers_csv, low_memory=False)
            if "catalogNumber" in vouchers_df.columns and "catalogNumber" in extracted_df.columns:
                active_vcats = set(vouchers_df["catalogNumber"].astype(str).str.strip())
                extracted_df = extracted_df[extracted_df["catalogNumber"].astype(str).str.strip().isin(active_vcats)].copy()
        except Exception:
            vouchers_df = None

    gen_fn = _resolve_facade_hook("generate_step02_publication_figures", generate_step02_publication_figures, _ORIG_GEN_STEP02_FIGS)
    return gen_fn(
        extracted_df=extracted_df, failed_df=failed_df, repro_df=repro_df,
        pheno_df=pheno_df, figures_dir=figures_dir, vouchers_df=vouchers_df,
    )


regenerate_step02_figures_from_tables_impl = regenerate_step02_figures_from_tables

__all__ = [
    "render_plate02b_qc_and_phenophase", "generate_step02_publication_figures",
    "generate_step02_publication_figures_impl", "regenerate_step02_figures_from_tables",
    "regenerate_step02_figures_from_tables_impl",
]
