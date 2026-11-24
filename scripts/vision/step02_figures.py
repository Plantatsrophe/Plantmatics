#!/usr/bin/env python3
"""
===============================================================================
Module: step02_figures.py
Project: Packera Species Delimitation & Integrative Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Backward-compatible facade and CLI entrypoint for the Step 02 300-DPI PNG
    publication figure suite and native LeafMachine3 visual asset exporters:
      1. fig02a_lm3_organ_and_ect_exemplars.png:
         Dynamic taxon-stratified organ montage delegating directly to native
         LeafMachine3 `render_leaf_products` and `render_radial_ect_overlay`
         across Basal Rosette and Cauline Stem tracks.
      2. fig02b_segmentation_qc_and_phenophase.png:
         Multi-panel segmentation QC, ruler calibration, basal-vs-cauline
         solidity x aspect-ratio homology separation (unfilled hollow markers,
         reniform L/W < 1.0 support), involucral AR_inv distributions, and
         determinate cyme phenophase proportions.
      3. lm3_assets/ & ect_overlays/:
         Cached native LeafMachine3 exemplar silhouettes and polar ECT PNGs.
      4. lm3_sheet_overlays/:
         Capped whole-sheet LeafMachine3 diagnostic overlays for priority sheets.
===============================================================================
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import scripts.vision.lm3_compat  # noqa: F401
from scripts.vision.lm3_inference_engines import _ensure_ect_module_loaded  # noqa: F401
from scripts.vision.step02_fig_ect_exemplars import (
    _render_lm3_ect_rgb,
    _render_lm3_silhouette_card,
    export_taxon_ect_exemplars,
    export_taxon_lm3_exemplars,
)
from scripts.vision.step02_fig_plate02a import (
    _is_exemplar_eligible,
    _pick_distinct_quantile_rows,
    _resolve_atomic_fig_saver,
    prepare_step02_figure_context,
    render_plate02a_and_ect_cache,
)
from scripts.vision.step02_fig_plate02b import (
    generate_step02_publication_figures,
    regenerate_step02_figures_from_tables,
    render_plate02b_qc_and_phenophase,
)
from scripts.vision.step02_fig_sheet_overlays import (
    _draw_solidity_vs_ar_zone_panel,
    _repel_medoid_annotations,
    export_lm3_sheet_diagnostic_overlays,
)
from scripts.vision.step02_fig_style_assets import (
    HOLLOW_MARKERS_5,
    ICN_ROMAN_TOKENS,
    TOL_OKABE_12,
    _ASSET_ANCHOR_DIRS,
    _atomic_save_matplotlib_fig,
    _discover_downstream_table,
    _is_valid_smooth_closed_mask,
    _load_mask_or_contour_u8,
    _parse_bool_val,
    _resolve_asset_path,
    _standardize_packera_species_py,
    format_icn_taxon_mathtext,
    get_combinatorial_taxon_style,
)

logger = logging.getLogger("SegmentAndExtract")

__all__ = [
    "_atomic_save_matplotlib_fig",
    "_standardize_packera_species_py",
    "format_icn_taxon_mathtext",
    "get_combinatorial_taxon_style",
    "_parse_bool_val",
    "_ASSET_ANCHOR_DIRS",
    "_resolve_asset_path",
    "_is_valid_smooth_closed_mask",
    "_discover_downstream_table",
    "_load_mask_or_contour_u8",
    "TOL_OKABE_12",
    "HOLLOW_MARKERS_5",
    "ICN_ROMAN_TOKENS",
    "generate_step02_publication_figures",
    "regenerate_step02_figures_from_tables",
    "export_taxon_ect_exemplars",
    "export_taxon_lm3_exemplars",
    "export_lm3_sheet_diagnostic_overlays",
    "_render_lm3_silhouette_card",
    "_render_lm3_ect_rgb",
    "_repel_medoid_annotations",
    "_draw_solidity_vs_ar_zone_panel",
    "_resolve_atomic_fig_saver",
    "_is_exemplar_eligible",
    "_pick_distinct_quantile_rows",
    "prepare_step02_figure_context",
    "render_plate02a_and_ect_cache",
    "render_plate02b_qc_and_phenophase",
]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 02 300-DPI PNG Figures and LeafMachine3 Exemplar Exporter")
    parser.add_argument("--export-ect-exemplars", action="store_true", help="Export authentic LeafMachine3 radial ECT & silhouette overlays")
    parser.add_argument("--output-dir", type=str, default="outputs/figures", help="Output directory for figures/overlays")
    parser.add_argument("--zone", type=str, default="all", choices=["basal", "cauline", "all"], help="Leaf zone (basal, cauline, or all)")
    parser.add_argument("--vouchers-csv", type=str, default="", help="Path to curated_vouchers.csv")
    parser.add_argument("--manifest-csv", type=str, default="", help="Path to extracted_leaf_manifest.csv")
    parser.add_argument("--regenerate-figures-only", action="store_true", help="Regenerate Step 02 publication figures")
    args = parser.parse_args()

    if args.export_ect_exemplars:
        out_fdir = Path(args.output_dir)
        v_csv = args.vouchers_csv or None
        ext_in = pd.read_csv(args.manifest_csv, low_memory=False) if (args.manifest_csv and Path(args.manifest_csv).exists()) else None
        zones = ["basal", "cauline"] if args.zone == "all" else [args.zone]
        for z in zones:
            export_taxon_ect_exemplars(extracted_df=ext_in, figures_dir=out_fdir, zone=z, max_taxa=None, vouchers_csv=v_csv)
    elif args.regenerate_figures_only:
        regenerate_step02_figures_from_tables(output_dir=args.output_dir, vouchers_csv=args.vouchers_csv or None)
