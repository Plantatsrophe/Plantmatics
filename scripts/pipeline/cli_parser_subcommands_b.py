"""CLI Group B subcommands (`morphometrics`, `genomics`, `skim-organellar`, `download-public-*`, `synthesis`, `run-all`)."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from scripts.pipeline.cli_parser_subcommands_a import (
    _add_common_subparser_args,
    _add_curation_parity_args,
)


def register_subcommands_b(subparsers: Any) -> None:
    """Register morphometrics, genomics, public downloaders, synthesis, and run-all subcommands."""
    # Subcommand: morphometrics
    p_morph = subparsers.add_parser("morphometrics", help="Phase 3: Fourier EFA, GMM, & passive CDA in R")
    _add_common_subparser_args(p_morph, include_gpu=False)
    p_morph.add_argument("--input", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--manifest", type=Path, default=None, help="Extracted leaf manifest CSV")
    p_morph.add_argument("--contours-dir", type=Path, default=None, help="Input contours directory")
    p_morph.add_argument("--vouchers", type=Path, default=None, help="Curated vouchers CSV")
    p_morph.add_argument("--reproductive", type=Path, default=None, help="Voucher reproductive metrics CSV")
    p_morph.add_argument("--harmonics-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--flags-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--plot-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--report-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--allometry-report", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--allometry-pairwise", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--allometry-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--ecoregion-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--output-batch-audit", type=Path, default=None, help=argparse.SUPPRESS)
    p_morph.add_argument("--harmonics", type=int, default=None, help="Number of Fourier harmonics (k)")
    p_morph.add_argument("--num-pcs", type=int, default=None, help="Number of PCA dimensions")
    p_morph.add_argument("--max-k", type=int, default=None, help="Max GMM mixture components")
    p_morph.add_argument("--focal-taxon", type=str, default=None, help="Optional focal Packera taxon for targeted morphometric CDA")
    p_morph.add_argument("--intraspecific-mode", action="store_true", default=False, help="Run single-taxon intraspecific delimitation mode")

    # Subcommand: genomics (Step 06 standalone phylogenomics)
    p_geno = subparsers.add_parser("genomics", help="Step 06: Independent genomic ingestion & phylogenomics")
    _add_common_subparser_args(p_geno, include_gpu=False)
    p_geno.add_argument("--genomics-dir", type=Path, default=None, help="Directory holding staged genomic files (default: data/genomics)")

    # Subcommand: skim-organellar (offline off-target chloroplast/ITS skimming)
    p_skim = subparsers.add_parser("skim-organellar", help="Offline skimming of off-target plastid/ITS/waxy reads into data/genomics/public_organellar_and_loci/")
    p_skim.add_argument("--fastq-source", type=Path, default=None, help="Optional .zip, .tar.gz, or directory containing raw FASTQ reads")
    p_skim.add_argument("--genomics-dir", type=Path, default=None, help="Destination data/genomics/ directory")
    p_skim.add_argument("--nuclear-tree", type=Path, default=None, help="Optional nuclear Newick tree for cytonuclear comparison")
    p_skim.add_argument("--out-dir", type=Path, default=None, help="Optional output directory for skimmed FASTAs and cytonuclear CSV (default: <genomics-dir>/public_organellar_and_loci)")
    p_skim.add_argument("-j", "--threads", "--workers", dest="threads", type=int, default=None, help="Maximum CPU worker processes for parallel FASTQ skimming (default: auto).")

    # Subcommand: download-public-genomics (public NCBI/SRA/Dryad/BOLD staging into relational data/genomics/ subdirectories)
    p_dl_pub = subparsers.add_parser("download-public-genomics", help="Download public Packera sequences, trees, and SRA manifests into relational data/genomics/ subdirectories")
    p_dl_pub.add_argument("--genomics-dir", type=Path, default=None, help="Destination data/genomics/ directory")
    p_dl_pub.add_argument("--sources", nargs="+", default=None, help="Public sources to query (ncbi, dryad, bold, paraloss)")
    p_dl_pub.add_argument("--genus", type=str, default=None, help="Target genus for dynamic ENA taxon-tree SRA discovery")
    p_dl_pub.add_argument("--family", type=str, default=None, help="Target family for Compositae-ParaLoss-1272 metadata")
    p_dl_pub.add_argument("--dryad-zip-dir", type=Path, default=None, help="Directory holding manually downloaded Dryad .zip archives")
    p_dl_pub.add_argument("--skim-sra-accessions", nargs="+", default=None, help="SRA/ENA run accessions to stream-skim via skim_offtarget_organellar")
    p_dl_pub.add_argument("--force", action="store_true", default=False, help="Overwrite existing non-empty public_packera_* files")

    # Subcommand: download-public-cytology (public CCDB/ICNA/FNA/label staging into flat data/cytology/)
    p_dl_cyto = subparsers.add_parser("download-public-cytology", help="Download and cross-link public Packera cytological records into flat data/cytology/")
    p_dl_cyto.add_argument("--cytology-dir", type=Path, default=None, help="Flat destination data/cytology/ directory")
    p_dl_cyto.add_argument("--sources", nargs="+", default=None, help="Public sources to query (ccdb, icna, fna, labels)")
    p_dl_cyto.add_argument("--local-cytology-csv", type=Path, default=None, help="Optional local CSV/TSV file with curated chromosome counts or C-values")
    p_dl_cyto.add_argument("--sync-db-only", action="store_true", default=False, help="Skip network downloads and re-sync existing flat CSVs and label remarks into SQLite")
    p_dl_cyto.add_argument("--no-db", action="store_true", default=False, help="Write flat CSV/JSON artifacts without updating SQLite")
    p_dl_cyto.add_argument("--force", action="store_true", default=False, help="Overwrite existing non-empty public_packera_cytology_* files")

    # Subcommand: download-taxonomic-backbone
    p_dl_backbone = subparsers.add_parser("download-taxonomic-backbone", help="Download dynamic taxonomic backbone profiles and synonymy into flat data/taxonomy/")
    p_dl_backbone.add_argument("--taxonomy-dir", type=Path, default=None, help="Flat destination data/taxonomy/ directory")
    p_dl_backbone.add_argument("--genus", type=str, default=None, help="Target genus for taxonomic backbone discovery")
    p_dl_backbone.add_argument("--family", type=str, default=None, help="Target botanical family filter")
    p_dl_backbone.add_argument("--sources", nargs="+", default=None, help="Optional taxonomic backbone sources")
    p_dl_backbone.add_argument("--sync-only", "--sync-db-only", dest="sync_only", action="store_true", default=False, help="Skip network downloads and sync cached backbone into SQLite")
    p_dl_backbone.add_argument("--db", type=Path, default=None, help="Optional override path to SQLite database")
    p_dl_backbone.add_argument("--no-db", action="store_true", default=False, help="Write flat CSV/JSON artifacts without updating SQLite")
    p_dl_backbone.add_argument("--force", action="store_true", default=False, help="Overwrite existing taxonomic backbone cache files")

    # Subcommand: synthesis
    p_synth = subparsers.add_parser("synthesis", aliases=["analyze"], help="Steps 05-07: Spatial RF, Phylogenomics, & Triage Dashboard")
    _add_common_subparser_args(p_synth, include_gpu=True)
    _add_curation_parity_args(p_synth)
    p_synth.add_argument("--skip-nwi-network", action="store_true", default=False, help="Reuse cached NWI points without REST queries")
    p_synth.add_argument("--vouchers", type=Path, default=None, help="Curated vouchers CSV")
    p_synth.add_argument("--manifest", dest="manifest_csv", type=Path, default=None, help="Extracted leaf manifest CSV")
    p_synth.add_argument("--manifest-csv", dest="manifest_csv", type=Path, default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    p_synth.add_argument("--morphometrics", dest="morphometrics", type=Path, default=None, help="Morphometrics flags CSV")
    p_synth.add_argument("--morphometrics-csv", dest="morphometrics", type=Path, default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    p_synth.add_argument("--pheno-states", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--export-figures", action="store_true", default=True, help="Export diagnostic figures")
    p_synth.add_argument("--no-export-figures", dest="export_figures", action="store_false", help="Disable figure export")
    p_synth.add_argument("--env-dir", type=Path, default=None, help="Directory containing environmental rasters")
    p_synth.add_argument("--skip-env-download", action="store_true", default=False, help="Skip preflight environmental layer download check")
    p_synth.add_argument("--force-env-download", action="store_true", default=False, help="Force re-download of SOLUS100, SGMC, and 3DEP layers")
    p_synth.add_argument("--skip-cytology-download", action="store_true", default=False, help="Skip preflight public cytology download and sync check")
    p_synth.add_argument("--force-cytology-download", action="store_true", default=False, help="Force fresh public cytology download before Steps 06-07")
    p_synth.add_argument("--skip-backbone-download", action="store_true", default=False, help="Skip preflight taxonomic backbone download and sync check")
    p_synth.add_argument("--force-backbone-download", action="store_true", default=False, help="Force fresh taxonomic backbone download before synthesis")
    p_synth.add_argument("--ssurgo-batch-size", type=int, default=None, help="Batch size for USDA SDA SSURGO queries")
    p_synth.add_argument("--ssurgo-table", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--contingency-table", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--pheno-anomalies", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--pheno-summary", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--pheno-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--conflict-flags", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--spatial-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--sdm-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--spatial-summary", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--permutations", type=int, default=999, help="Warren's identity test permutations (default: 999)")
    p_synth.add_argument("--n-trees", type=int, default=500, help="Random forest trees")
    p_synth.add_argument("--gmm-summary", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--output-queue", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--synthesis-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--synthesis-report", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--synthesis-html", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--min-affinity", type=float, default=None, help="Minimum multimodal affinity threshold for Step 07 complex neighbors (default: 0.03)")
    p_synth.add_argument("--output-complexes", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--output-affinity", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--output-neighbors", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--output-complex-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_synth.add_argument("--focal-taxon", type=str, default=None, help="Optional focal Packera taxon for targeted revision synthesis")
    p_synth.add_argument("--intraspecific-mode", action="store_true", default=False, help="Run single-taxon intraspecific delimitation mode")

    # Subcommand: run-all
    p_all = subparsers.add_parser("run-all", help="Execute complete pipeline (harvest -> segment -> morphometrics -> synthesis)")
    _add_common_subparser_args(p_all, include_gpu=True)
    _add_curation_parity_args(p_all)
    p_all.add_argument("--skip-nwi-network", action="store_true", default=False, help="Reuse cached NWI points without REST queries")
    p_all.add_argument(
        "--taxon",
        dest="taxa",
        nargs="+",
        default=None,
        help="Target taxonomic binomials (omit to harvest all Packera species genus-wide)",
    )
    p_all.add_argument(
        "--taxa",
        dest="taxa",
        nargs="+",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    p_all.add_argument(
        "--geographic-scope",
        choices=["eastern_na", "north_america", "global"],
        default=None,
        help="Geographic scope for GBIF harvesting (eastern_na, north_america, global)",
    )
    p_all.add_argument("--countries", nargs="+", default=None, help="ISO-2 country codes for GBIF harvesting")
    p_all.add_argument(
        "--expand-synonyms",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Toggle dynamic GBIF Backbone unaccepted Packera + Senecio/Cineraria synonym unification (default: True)",
    )
    p_all.add_argument(
        "--usda-county-filter",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Toggle USDA PLANTS county range filtering during Step 01 harvest (default: from config.yaml)",
    )
    p_all.add_argument(
        "--usda-county-buffer-km",
        type=float,
        default=None,
        help="Geodesic buffer in km around USDA PLANTS documented counties (default: from config.yaml)",
    )
    p_all.add_argument(
        "--usda-exempt-gold",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Toggle Tier_1_Gold specialist exemption during USDA county range filtering",
    )
    p_all.add_argument(
        "--harvest-full-continental-range",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Toggle full US + CA continental range harvesting for verified Eastern NA resident species",
    )
    p_all.add_argument(
        "--spatial-thin-km",
        type=float,
        default=None,
        help="Minimum spatial grid thinning distance in km per species (default: 0.0, disabled)",
    )
    p_all.add_argument(
        "--max-per-state",
        type=int,
        default=None,
        help="Maximum vouchers retained per (species, stateProvince)",
    )
    p_all.add_argument(
        "--max-per-region",
        type=int,
        default=None,
        help="Maximum vouchers retained per (species, regional_group)",
    )
    p_all.add_argument("--include-western", action="store_true", default=False, help=argparse.SUPPRESS)
    p_all.add_argument("--out-dir", type=Path, default=None, help="Output directory for raw vouchers")
    p_all.add_argument("--output-csv", type=Path, default=None, help="Curated vouchers output CSV")
    p_all.add_argument("--output-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--ecoregion-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--sdm-plot", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--max-records", type=int, default=None, help="Max records per taxon")
    p_all.add_argument("--max-uncertainty", type=float, default=None, help="Max georeference uncertainty (m)")
    p_all.add_argument("--min-megapixels", type=float, default=None, help="Min resolution in Megapixels")
    p_all.add_argument("--min-file-size-kb", type=float, default=None, help="Min file size in KB")
    p_all.add_argument("--concurrency", type=int, default=None, help="Download concurrency")
    p_all.add_argument("--download-images", action="store_true", default=False, help="Download specimen images")
    p_all.add_argument(
        "--purge-segments",
        "--purge-old-data",
        dest="purge_segments",
        action="store_true",
        default=False,
        help="Purge all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables before running.",
    )
    p_all.add_argument(
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-harvesting and re-segmentation, overwriting existing artifacts.",
    )
    p_all.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    p_all.add_argument("--vouchers", type=Path, default=None, help="Input curated vouchers CSV")
    p_all.add_argument("--weights", type=Path, default=None, help="LeafMachine3 model weights checkpoint (.onnx / .pt)")
    p_all.add_argument("--output-dir", type=Path, default=None, help="Root output directory")
    p_all.add_argument("--device", type=str, default=None, choices=["cuda", "cpu"], help="Inference device")
    p_all.add_argument("--score-thresh", type=float, default=None, help="PCD detection score threshold")
    p_all.add_argument("--limit", type=int, default=None, help="Limit number of vouchers to process")
    p_all.add_argument("--input", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--manifest", type=Path, default=None, help="Extracted leaf manifest CSV")
    p_all.add_argument("--contours-dir", type=Path, default=None, help="Input contours directory")
    p_all.add_argument("--harmonics-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--flags-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--plot-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--report-out", type=Path, default=None, help=argparse.SUPPRESS)
    p_all.add_argument("--harmonics", type=int, default=None, help="Number of Fourier harmonics (k)")
    p_all.add_argument("--num-pcs", type=int, default=None, help="Number of PCA dimensions")
    p_all.add_argument("--max-k", type=int, default=None, help="Max GMM mixture components")
    p_all.add_argument("--ssurgo-batch-size", type=int, default=None, help="Batch size for USDA SDA SSURGO queries")
    p_all.add_argument("--skip-env-download", action="store_true", default=False, help="Skip preflight environmental layer download check")
    p_all.add_argument("--force-env-download", action="store_true", default=False, help="Force re-download of SOLUS100, SGMC, and 3DEP layers")
    p_all.add_argument("--skip-cytology-download", action="store_true", default=False, help="Skip preflight public cytology download and sync check")
    p_all.add_argument("--force-cytology-download", action="store_true", default=False, help="Force fresh public cytology download before Steps 06-07")
    p_all.add_argument("--skip-backbone-download", action="store_true", default=False, help="Skip preflight taxonomic backbone download and sync check")
    p_all.add_argument("--force-backbone-download", action="store_true", default=False, help="Force fresh taxonomic backbone download before synthesis")
    p_all.add_argument("--export-figures", action="store_true", default=True, help="Export diagnostic figures")
    p_all.add_argument("--no-export-figures", dest="export_figures", action="store_false", help="Disable figure export")
    p_all.add_argument("--focal-taxon", type=str, default=None, help="Optional focal Packera taxon for targeted revision synthesis")
    p_all.add_argument("--intraspecific-mode", action="store_true", default=False, help="Run single-taxon intraspecific delimitation mode")
    p_all.add_argument("--min-affinity", type=float, default=None, help="Minimum multimodal affinity threshold for Step 07 complex neighbors (default: 0.03)")
