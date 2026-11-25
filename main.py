#!/usr/bin/env python3
"""Unified Plantmatics CLI facade orchestrating Steps 01 to 07 across SQLite WAL and CSV snapshots."""

from __future__ import annotations

import argparse
import importlib
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Canonical Momocs GitHub remediation command referenced by check-env diagnostics
MOMOCS_GITHUB_REMEDIATION = 'Rscript -e \'remotes::install_github("MomX/Momocs")\''

from scripts.core.artifact_manager import ensure_model_weights
from scripts.core.config import (
    PipelineConfig,
    ensure_run_curated_vouchers_snapshot,
    generate_default_run_id,
    get_active_run_dir,
    initialize_run_directory,
    resolve_canonical_run_id,
    write_run_metadata,
)
from scripts.core.db import purge_runs_only, purge_segments_and_runs
from scripts.pipeline.cli_env_doctor import (
    check_environment,
    run_check_env,
    run_download_weights,
)
from scripts.pipeline.cli_parser_subcommands_a import (
    _add_common_subparser_args,
    _add_curation_parity_args,
    register_subcommands_a,
)
from scripts.pipeline.cli_parser_subcommands_b import register_subcommands_b
from scripts.pipeline.cli_run_context import (
    _COHORT_SYNC_CACHE,
    _apply_explicit_db_arg,
    _attach_run_log_handler,
    _db_has_extracted_contours,
    _resolve_active_run_context,
    _resolve_focal_and_intraspecific,
    _sync_run_cohort_in_db,
    check_gpu_availability,
    logger,
    verify_dir_has_files,
    verify_file_exists,
)
from scripts.pipeline.cli_stage_genomics_dl import (
    run_download_public_cytology_cli,
    run_download_public_genomics_cli,
    run_download_taxonomic_backbone_cli,
    run_genomics_cli,
    run_genomics_stream,
    run_skim_organellar_cli,
)
from scripts.pipeline.cli_stage_harvest_seg import (
    get_lm3_python_executable,
    run_harvest,
    run_segment,
)
from scripts.pipeline.cli_stage_morph_layers import (
    _maybe_ensure_cytology_records,
    _maybe_ensure_environmental_layers,
    _maybe_ensure_taxonomic_backbone,
    run_download_env,
    run_morphometrics,
)
from scripts.pipeline.cli_stage_synthesis import (
    run_all,
    run_purge_runs_only,
    run_purge_segments,
    run_synthesis,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="main.py",
        description="Packera Species Delimitation & Morphometrics Production Pipeline.",
    )
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=None,
        help="Path to YAML configuration file (default: config/config.yaml).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose debugging output.",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Optional run identifier for isolated run directory under runs/<run_id>/ (default: run_YYYYMMDD_HHMMSS).",
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Optional override path to SQLite database (default: data/plantmatics.sqlite).",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        default=False,
        help="Run in ephemeral test isolation mode (PACKERA_DEMO_MODE=1).",
    )
    parser.add_argument(
        "-j", "--threads", "--workers",
        dest="threads",
        type=int,
        default=None,
        help="Maximum CPU worker threads across local Python and R stages (default: auto).",
    )
    parser.add_argument(
        "--gpu-batch-size",
        dest="gpu_batch_size",
        type=int,
        default=None,
        help="GPU inference batch size for LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime) (default: stage-optimal).",
    )
    parser.add_argument(
        "--batch-size",
        dest="gpu_batch_size",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--purge-segments",
        "--purge-old-data",
        dest="purge_segments",
        action="store_true",
        default=False,
        help="Purge all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables while preserving Step 01 vouchers and raw images.",
    )
    parser.add_argument(
        "--purge-runs-only",
        "--purge-run-outputs",
        dest="purge_runs_only",
        action="store_true",
        default=False,
        help="Purge only Step 03 to 07 run outputs and Tier B SQLite tables while preserving Step 01 vouchers and Step 02 segmentation.",
    )
    parser.add_argument(
        "--no-purge-runs",
        dest="purge_runs",
        action="store_false",
        default=True,
        help="When purging segments, preserve historical runs/ directories.",
    )

    subparsers = parser.add_subparsers(dest="subcommand", required=False, help="Pipeline phase to execute")
    register_subcommands_a(subparsers)
    register_subcommands_b(subparsers)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    _apply_explicit_db_arg(args)

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if args.subcommand == "purge":
        if getattr(args, "runs_only", False) and not getattr(args, "segments", False):
            run_purge_runs_only(args)
        elif getattr(args, "segments", False):
            run_purge_segments(args)
        else:
            parser.error("purge requires --runs-only or --segments")
        return

    if args.subcommand in ("purge-runs-only", "purge-run-outputs"):
        run_purge_runs_only(args)
        return

    if getattr(args, "purge_runs_only", False) and not args.subcommand:
        run_purge_runs_only(args)
        return

    if args.subcommand in ("purge-segments", "purge-old-data"):
        run_purge_segments(args)
        return

    if getattr(args, "purge_segments", False) and not args.subcommand:
        run_purge_segments(args)
        return

    if not args.subcommand:
        parser.error("the following arguments are required: subcommand")

    if args.subcommand == "check-env":
        run_check_env(args)
        return

    if getattr(args, "purge_segments", False) and args.subcommand != "segment":
        run_purge_segments(args)

    cfg = PipelineConfig.from_yaml(args.config)

    dispatch = {
        "check-env": run_check_env,
        "download-weights": run_download_weights,
        "download-env": run_download_env,
        "purge-segments": run_purge_segments,
        "purge-old-data": run_purge_segments,
        "purge-runs-only": run_purge_runs_only,
        "purge-run-outputs": run_purge_runs_only,
        "harvest": run_harvest,
        "segment": run_segment,
        "morphometrics": run_morphometrics,
        "genomics": run_genomics_cli,
        "skim-organellar": run_skim_organellar_cli,
        "download-public-genomics": run_download_public_genomics_cli,
        "download-public-cytology": run_download_public_cytology_cli,
        "download-taxonomic-backbone": run_download_taxonomic_backbone_cli,
        "synthesis": run_synthesis,
        "analyze": run_synthesis,
        "run-all": run_all,
    }

    handler = dispatch.get(args.subcommand)
    if handler:
        handler(args, cfg)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
