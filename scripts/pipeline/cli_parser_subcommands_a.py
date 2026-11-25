"""CLI parser helpers and Group A subcommands (`purge*`, `check-env`, `download-*`, `harvest`, `segment`)."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any


class _CommaAppendAction(argparse.Action):
    """Accumulate single or repeated `--triage-source` values as a comma-separated string."""

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: Any, option_string: str | None = None) -> None:
        val = str(values).strip()
        if not val:
            return
        cur = getattr(namespace, self.dest, None)
        setattr(namespace, self.dest, val if (not cur or str(cur).upper() == "ALL") else f"{cur},{val}")


def _add_common_subparser_args(subparser: argparse.ArgumentParser, include_gpu: bool = False) -> None:
    """Attach shared run-scoping, database, demo, and compute flags to a subparser."""
    subparser.add_argument("--run-id", type=str, default=None, help="Run identifier for isolated run directory under runs/<run_id>/")
    subparser.add_argument("--db", type=Path, default=argparse.SUPPRESS, help="Optional override path to SQLite database")
    subparser.add_argument("--demo", action="store_true", default=False, help="Run in ephemeral test isolation mode (PACKERA_DEMO_MODE=1).")
    subparser.add_argument("-j", "--threads", "--workers", dest="threads", type=int, default=None, help="Maximum CPU worker threads across local Python and R stages (default: auto).")
    if include_gpu:
        subparser.add_argument("--gpu-batch-size", dest="gpu_batch_size", type=int, default=None, help="GPU inference batch size for LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime) (default: stage-optimal).")
        subparser.add_argument("--batch-size", dest="gpu_batch_size", type=int, default=argparse.SUPPRESS, help=argparse.SUPPRESS)
        subparser.add_argument("--wave-size", dest="wave_size", type=int, default=None, help="Number of herbarium sheets processed per cross-voucher GPU wave (default: 16 or --gpu-batch-size).")
        subparser.add_argument("--prefetch-workers", dest="prefetch_workers", type=int, default=None, help="Optional override for Stage A disk image prefetch worker thread count.")
        subparser.add_argument("--prep-workers", dest="prep_workers", type=int, default=None, help="Optional override for Stage B2/B3 CPU crop prep and prototype mask decode worker thread count.")
        subparser.add_argument("--post-workers", dest="post_workers", type=int, default=None, help="Optional override for Stage C contour/ECT/rosette post-processing worker thread count.")


def _add_curation_parity_args(subparser: argparse.ArgumentParser) -> None:
    """Attach run_pipeline.sh curation & complex-neighbor parity flags to a subparser."""
    subparser.add_argument("--gold-only", action="store_true", default=False, help="Restrict active cohort to Nomenclatural Types and Tier_1_Gold vouchers")
    subparser.add_argument("--determiner-tier", type=str, default=None, help="Restrict cohort to determiner authority tier (gold|silver|bronze or Tier_1_Gold)")
    subparser.add_argument("--exclude-misidentified", action="store_true", default=False, help="Prune vouchers flagged as Misidentification_Reassign in prior runs")
    subparser.add_argument("--triage-source", dest="triage_source", action=_CommaAppendAction, default=None, help="Source run ID(s) or ALL for --exclude-misidentified")
    subparser.add_argument("--auto-neighbors", action="store_true", default=False, help="Append high-affinity neighbor taxa from Step 07 complex network output")


def register_subcommands_a(subparsers: Any) -> None:
    """Register purge, diagnostics, provisioning, harvest, and segment subcommands."""
    # Subcommand: purge-segments
    p_purge = subparsers.add_parser(
        "purge-segments",
        aliases=["purge-old-data"],
        help="Purge all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables",
    )
    _add_common_subparser_args(p_purge, include_gpu=False)
    p_purge.add_argument(
        "--no-purge-runs",
        dest="purge_runs",
        action="store_false",
        default=True,
        help="Preserve historical runs/ directories when purging Step 02 segments.",
    )

    # Subcommand: purge-runs-only
    p_purge_runs = subparsers.add_parser(
        "purge-runs-only",
        aliases=["purge-run-outputs"],
        help="Purge only Step 03 to 07 run outputs and Tier B SQLite tables while preserving Step 01 & Step 02 segmentation",
    )
    _add_common_subparser_args(p_purge_runs, include_gpu=False)
    p_purge_runs.add_argument("--output-dir", type=Path, default=None, help="Optional root data directory override")

    # Subcommand: purge (--runs-only / --segments)
    p_purge_cmd = subparsers.add_parser(
        "purge",
        help="Purge pipeline artifacts (--runs-only for Tier B runs or --segments for Step 02 to 07)",
    )
    _add_common_subparser_args(p_purge_cmd, include_gpu=False)
    p_purge_cmd.add_argument(
        "--runs-only",
        action="store_true",
        default=False,
        help="Purge Tier B run outputs and genetic_lineage_assignments while preserving Tier A vouchers, segments, and genomic_samples",
    )
    p_purge_cmd.add_argument(
        "--segments",
        action="store_true",
        default=False,
        help="Purge Step 02 to 07 segments, masks, contours, and run directories",
    )
    p_purge_cmd.add_argument(
        "--no-purge-runs",
        dest="purge_runs",
        action="store_false",
        default=True,
        help="Preserve historical runs/ directories when purging Step 02 segments.",
    )
    p_purge_cmd.add_argument("--output-dir", type=Path, default=None, help="Optional root data directory override")

    # Subcommand: check-env
    p_check = subparsers.add_parser(
        "check-env",
        help="Run comprehensive preflight environment & dependency diagnostics",
    )
    p_check.add_argument(
        "--strict",
        action="store_true",
        default=False,
        help="Treat warnings (such as missing CUDA) as critical failures.",
    )

    # Subcommand: download-weights
    p_dl = subparsers.add_parser(
        "download-weights",
        help="Download and verify LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime) model weights from remote repository",
    )
    p_dl.add_argument(
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-download even if model weights file exists locally.",
    )
    p_dl.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )

    # Subcommand: download-env
    p_env = subparsers.add_parser(
        "download-env",
        help="Verify or download SOLUS100 (100 m), USGS SGMC, and 3DEP 30 m DEM environmental layers",
    )
    _add_common_subparser_args(p_env, include_gpu=False)
    p_env.add_argument("--env-dir", type=Path, default=None, help="Target environmental directory (default: data/environmental)")
    p_env.add_argument("--vouchers", type=Path, default=None, help="Optional curated_vouchers.csv to scope 30 m DEM tiles")
    p_env.add_argument("--genus", type=str, default=None, help="Target genus for USACE NWPL wetland ratings")
    p_env.add_argument("--skip-nwi-network", action="store_true", default=False, help="Reuse cached NWI points without REST queries")
    p_env.add_argument("--force", action="store_true", default=False, help="Force re-download even if cached layers exist")

    # Subcommand: harvest
    p_harvest = subparsers.add_parser("harvest", help="Phase 1: Ingest GBIF occurrences & download images")
    _add_common_subparser_args(p_harvest, include_gpu=False)
    _add_curation_parity_args(p_harvest)
    p_harvest.add_argument(
        "--taxon",
        dest="taxa",
        nargs="+",
        default=None,
        help="Target taxonomic binomials (omit to harvest all Packera species genus-wide)",
    )
    p_harvest.add_argument(
        "--taxa",
        dest="taxa",
        nargs="+",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    p_harvest.add_argument(
        "--geographic-scope",
        choices=["eastern_na", "north_america", "global"],
        default=None,
        help="Geographic scope for GBIF harvesting (eastern_na, north_america, global)",
    )
    p_harvest.add_argument("--countries", nargs="+", default=None, help="ISO-2 country codes for GBIF harvesting")
    p_harvest.add_argument(
        "--expand-synonyms",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Toggle dynamic GBIF Backbone unaccepted Packera + Senecio/Cineraria synonym unification (default: True)",
    )
    p_harvest.add_argument(
        "--usda-county-filter",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Toggle USDA PLANTS county range filtering during Step 01 harvest (default: from config.yaml)",
    )
    p_harvest.add_argument(
        "--usda-county-buffer-km",
        type=float,
        default=None,
        help="Geodesic buffer in km around USDA PLANTS documented counties (default: from config.yaml)",
    )
    p_harvest.add_argument(
        "--usda-exempt-gold",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Toggle Tier_1_Gold specialist exemption during USDA county range filtering",
    )
    p_harvest.add_argument(
        "--harvest-full-continental-range",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Toggle full US + CA continental range harvesting for verified Eastern NA resident species",
    )
    p_harvest.add_argument(
        "--spatial-thin-km",
        type=float,
        default=None,
        help="Minimum spatial grid thinning distance in km per species (default: 0.0, disabled)",
    )
    p_harvest.add_argument(
        "--max-per-state",
        type=int,
        default=None,
        help="Maximum vouchers retained per (species, stateProvince)",
    )
    p_harvest.add_argument(
        "--max-per-region",
        type=int,
        default=None,
        help="Maximum vouchers retained per (species, regional_group)",
    )
    p_harvest.add_argument("--include-western", action="store_true", default=False, help=argparse.SUPPRESS)
    p_harvest.add_argument("--out-dir", type=Path, default=None, help="Output directory for raw vouchers")
    p_harvest.add_argument("--output-csv", type=Path, default=None, help="Curated vouchers output CSV")
    p_harvest.add_argument("--max-records", type=int, default=None, help="Max records per taxon")
    p_harvest.add_argument("--max-uncertainty", type=float, default=None, help="Max georeference uncertainty (m)")
    p_harvest.add_argument("--min-megapixels", type=float, default=None, help="Min resolution in Megapixels")
    p_harvest.add_argument("--min-file-size-kb", type=float, default=None, help="Min file size in KB")
    p_harvest.add_argument("--concurrency", type=int, default=None, help="Download concurrency")
    p_harvest.add_argument("--download-images", action="store_true", default=False, help="Download specimen images")
    p_harvest.add_argument(
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-download and overwrite existing cached voucher images.",
    )
    p_harvest.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )

    # Subcommand: segment
    p_segment = subparsers.add_parser("segment", help="Phase 2: Segment leaves & export contours")
    _add_common_subparser_args(p_segment, include_gpu=True)
    p_segment.add_argument("--vouchers", type=Path, default=None, help="Input curated vouchers CSV")
    p_segment.add_argument("--weights", type=Path, default=None, help="LeafMachine3 model weights checkpoint (.onnx / .pt)")
    p_segment.add_argument("--output-dir", type=Path, default=None, help="Root output directory")
    p_segment.add_argument("--device", type=str, default=None, choices=["cuda", "cpu"], help="Inference device")
    p_segment.add_argument("--score-thresh", type=float, default=None, help="LeafMachine3 leaf segmentation & Stage 1 confidence threshold")
    p_segment.add_argument("--limit", type=int, default=None, help="Limit number of vouchers to process")
    p_segment.add_argument(
        "--purge-segments",
        "--purge-old-data",
        dest="purge_segments",
        action="store_true",
        default=False,
        help="Purge all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables before segmenting.",
    )
    p_segment.add_argument(
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-segmentation and overwrite existing contour CSVs and masks.",
    )
    p_segment.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
