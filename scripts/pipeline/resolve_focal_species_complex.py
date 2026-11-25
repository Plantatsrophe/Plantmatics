#!/usr/bin/env python3
"""
scripts/pipeline/resolve_focal_species_complex.py
=================================================
Four-stage iterative species-complex discovery and intraspecific delimitation
orchestrator for Packera dubia (and any focal North American Packera taxon).

Executes four isolated, run-scoped passes through ./run_pipeline.sh backed by
shared Tier A SQLite caches (data/plantmatics.sqlite), then compiles a
cross-stage quantitative comparison CSV and a structured taxonomic adjudication
report in outputs/reports/:

  Stage 1 (<prefix>_01_eastern_discovery):
      Broad regional harvest, LeafMachine3 organ segmentation, and genus-level
      empirical taxon-complex network discovery.
  Stage 2 (<prefix>_02_complex_triage):
      Targeted focal-taxon boundary audit against auto-discovered high-affinity
      congeneric neighbors to separate Stage 1 Misidentification_Reassign sheets
      from Stage 2 Hybrid_Introgressant and Type_Synonymy_Candidate signals.
  Stage 3 (<prefix>_03_intraspecific_clean):
      Purified intraspecific mode (--intraspecific-mode --exclude-misidentified)
      dedicating 100% of EFA/ECT variance, plasticity ridgeline merging, 0 to 15 cm
      rhizosphere/outcrop contingency, and Hopkins circular phenology to internal
      lineage structure within the focal species sensu stricto.
  Stage 4 (<prefix>_04_gold_thin_audit):
      Fast Step 03+ spatially thinned sensitivity check (--spatial-thin-km) to
      confirm that candidate morphotypes or ecotypes survive collection-density
      de-clumping.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon
from scripts.pipeline.step07_focal_cohort_utils import (
    _extract_stage_metrics,
    _fetch_complex_rows_for_run,
    _fetch_triage_rows_for_run,
    _is_truthy,
    _read_csv_rows,
    _safe_float,
    _stage_run_ids,
    default_prefix_for_taxon,
    slugify_taxon,
)
from scripts.pipeline.step07_focal_stage_runners import (
    build_four_stage_commands,
    compile_complex_resolution_summary,
)

__all__ = [
    "PROJECT_ROOT",
    "standardize_packera_taxon",
    "slugify_taxon",
    "default_prefix_for_taxon",
    "_stage_run_ids",
    "_read_csv_rows",
    "_safe_float",
    "_is_truthy",
    "_fetch_triage_rows_for_run",
    "_fetch_complex_rows_for_run",
    "_extract_stage_metrics",
    "build_four_stage_commands",
    "compile_complex_resolution_summary",
    "build_parser",
    "main",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Execute and synthesize the 4-stage species-complex resolution workflow for Packera dubia or another focal Packera taxon."
    )
    parser.add_argument(
        "--focal-taxon",
        default="Packera dubia",
        help="Focal Packera binomial to resolve (default: 'Packera dubia').",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Run ID prefix for the four stages (default: derived from focal epithet, e.g. 'dubia').",
    )
    parser.add_argument(
        "--discovery-run-id",
        default=None,
        help="Optional existing Stage 1 discovery run ID to reuse as Stage 2 triage source.",
    )
    parser.add_argument(
        "--parent-run-id",
        default=None,
        help="Optional parent monograph run ID placing outputs under runs/<PARENT_RUN_ID>/<Species_Slug>/.",
    )
    parser.add_argument(
        "--geographic-scope",
        default="eastern_na",
        choices=("eastern_na", "north_america", "global"),
        help="Geographic scope for Stage 1 baseline discovery and downstream stages (default: eastern_na).",
    )
    parser.add_argument(
        "--gold-anchors",
        type=int,
        default=30,
        help="Target Tier 1 Gold / Silver anchors per neighbor taxon in Stage 2 (default: 30).",
    )
    parser.add_argument(
        "--spatial-thin-km",
        type=float,
        default=15.0,
        help="Spatial thinning distance in km for the Stage 4 sensitivity check (default: 15.0).",
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=None,
        help="CPU worker thread limit forwarded to run_pipeline.sh.",
    )
    parser.add_argument(
        "--gpu-batch-size",
        type=int,
        default=None,
        help="GPU inference batch size for Stage 1 LeafMachine3 segmentation.",
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Optional custom SQLite WAL database path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for cross-run comparison CSV and synthesis Markdown report (default: outputs/reports).",
    )
    parser.add_argument(
        "--gold-only",
        action="store_true",
        default=False,
        help="Propagate --gold-only across Stages 1 to 3 (Stage 4 always enforces --gold-only).",
    )
    parser.add_argument(
        "--skip-harvest",
        action="store_true",
        default=False,
        help="Reuse existing canonical Step 01 vouchers and Step 02 segments in Stage 1 without GBIF harvesting.",
    )
    parser.add_argument(
        "--skip-discovery",
        action="store_true",
        default=False,
        help="Skip Stage 1 when an Eastern NA discovery baseline has already been harvested and segmented.",
    )
    parser.add_argument(
        "--summarize-only",
        action="store_true",
        default=False,
        help="Skip pipeline execution and compile the cross-stage comparison CSV and Markdown report from existing runs.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print the exact four-stage run_pipeline.sh commands without executing them.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    focal_std = standardize_packera_taxon(args.focal_taxon) or args.focal_taxon.strip()
    parent_raw = (args.parent_run_id or "").strip().strip("/")
    if parent_raw.startswith("runs/"):
        parent_raw = parent_raw[len("runs/"):].strip("/")
    if parent_raw:
        prefix = f"{parent_raw}/{slugify_taxon(focal_std)}/"
        output_dir = args.output_dir if args.output_dir is not None else (
            PROJECT_ROOT / "runs" / parent_raw / slugify_taxon(focal_std)
        )
    else:
        prefix = args.prefix.strip() or default_prefix_for_taxon(focal_std)
        output_dir = args.output_dir
    run_script = PROJECT_ROOT / "run_pipeline.sh"

    stages = build_four_stage_commands(
        run_script=run_script,
        focal_taxon=focal_std,
        prefix=prefix,
        discovery_run_id=args.discovery_run_id,
        geographic_scope=args.geographic_scope,
        gold_anchors=args.gold_anchors,
        spatial_thin_km=args.spatial_thin_km,
        threads=args.threads,
        gpu_batch_size=args.gpu_batch_size,
        db_path=args.db,
        skip_discovery=args.skip_discovery,
        skip_harvest=args.skip_harvest,
        gold_only=args.gold_only,
    )

    if args.dry_run:
        for st_label, rid, cmd in stages:
            print(f"[{st_label}] ({rid}): {' '.join(cmd)}")
        return 0

    if not args.summarize_only:
        for st_label, rid, cmd in stages:
            print(f"\n=== Executing {st_label} ({rid}) ===", flush=True)
            proc = subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=False)
            if proc.returncode != 0:
                print(
                    f"[ERROR] {st_label} ({rid}) exited with non-zero status {proc.returncode}.",
                    file=sys.stderr,
                )
                return proc.returncode

    csv_out, md_out = compile_complex_resolution_summary(
        focal_taxon=focal_std,
        prefix=prefix,
        discovery_run_id=args.discovery_run_id,
        workspace_root=PROJECT_ROOT,
        db_path=args.db,
        output_dir=output_dir,
    )
    print(f"[Complex Resolution Summary] Comparison CSV : {csv_out}")
    print(f"[Complex Resolution Summary] Synthesis Report: {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
