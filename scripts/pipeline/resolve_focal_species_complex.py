#!/usr/bin/env python3
"""
scripts/pipeline/resolve_focal_species_complex.py
=================================================
Six-stage iterative species-complex discovery, misidentification purge,
unconfounded spatial/monographer sensitivity audits, and supervised candidate
lineage diagnosability orchestrator for Packera dubia (and any focal North
American Packera taxon).

Nests all six stage folders, cross-stage comparison CSVs, synthesis Markdown
report, and JSON run manifest inside a single unique timestamped parent run
directory `runs/<prefix>_YYYYMMDD_<hex6>/` (matching resolve_eastern_na_packera.sh):

  01_eastern_discovery:
      Broad regional harvest, LeafMachine3 organ segmentation, and genus-level
      empirical taxon-complex network discovery.
  02_complex_triage:
      Targeted focal-taxon boundary audit against auto-discovered complex
      neighbors to purge Misidentification_Reassign sheets and isolate hybrids.
  03_intraspecific_clean:
      Purified intraspecific mode (--intraspecific-mode --exclude-misidentified)
      across Gold + Silver + Bronze vouchers to discover candidate internal
      morphotypes, 30 m soil regimes, and Prim MST geographic disjunctions.
  04_spatial_thin_field_audit:
      Unconfounded 15 km habitat-stratified spatial thinning across Tier_1_Gold
      and Tier_2_Silver vouchers to verify geographic de-clumping with modern
      GPS field data intact.
  05_monographer_gold_audit:
      Strict Tier_1_Gold + nomenclatural types + G1 sequenced baseline at 15 km
      thinning to evaluate historical monographer concept stability.
  06_lineage_diagnosability:
      Supervised multi-group hypothesis test promoting Stage 3 candidate
      lineages to explicit operational groups to test the 75% morphological
      diagnosability rule, pairwise climatic/edaphic niche identity, dichotomous
      key traits, and WCVP/IPNI basionym coverage.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

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
from scripts.pipeline.step07_focal_parent_run import (
    _six_stage_run_ids,
    generate_focal_parent_run_id,
    resolve_focal_parent_run_id,
    stage_is_complete_focal,
    update_latest_focal_symlink,
    write_focal_run_manifest,
)
from scripts.pipeline.step07_focal_stage_runners import (
    build_four_stage_commands,
    build_six_stage_commands,
    compile_complex_resolution_summary,
)

__all__ = [
    "PROJECT_ROOT",
    "standardize_packera_taxon",
    "slugify_taxon",
    "default_prefix_for_taxon",
    "_stage_run_ids",
    "_six_stage_run_ids",
    "generate_focal_parent_run_id",
    "resolve_focal_parent_run_id",
    "_read_csv_rows",
    "_safe_float",
    "_is_truthy",
    "_fetch_triage_rows_for_run",
    "_fetch_complex_rows_for_run",
    "_extract_stage_metrics",
    "build_four_stage_commands",
    "build_six_stage_commands",
    "compile_complex_resolution_summary",
    "build_parser",
    "main",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Execute and synthesize the 6-stage species-complex resolution workflow inside a single parent run directory."
    )
    parser.add_argument("--focal-taxon", default="Packera dubia", help="Focal Packera binomial to resolve (default: 'Packera dubia').")
    parser.add_argument("--prefix", default="", help="Base prefix for unique parent run ID '<prefix>_YYYYMMDD_<hex6>' (default: focal epithet, e.g. 'dubia').")
    parser.add_argument("--parent-run-id", default=None, help="Explicit parent run directory under runs/<PARENT_RUN_ID>/ (or monograph parent).")
    parser.add_argument("--discovery-run-id", default=None, help="Optional existing Stage 1 discovery run ID to reuse as Stage 2 triage source.")
    parser.add_argument("--geographic-scope", default="eastern_na", choices=("eastern_na", "north_america", "global"), help="Geographic scope (default: eastern_na).")
    parser.add_argument("--gold-anchors", type=int, default=30, help="Target Tier 1 Gold / Silver anchors per neighbor taxon in Stage 2 (default: 30).")
    parser.add_argument("--spatial-thin-km", type=float, default=15.0, help="Spatial thinning distance in km for Stages 4 and 5 (default: 15.0).")
    parser.add_argument("-j", "--threads", type=int, default=None, help="CPU worker thread limit forwarded to run_pipeline.sh.")
    parser.add_argument("--gpu-batch-size", type=int, default=None, help="GPU inference batch size for Stage 1 LeafMachine3 segmentation.")
    parser.add_argument("--db", type=Path, default=None, help="Optional custom SQLite WAL database path.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional extra directory to mirror summary CSV and Markdown report (default: parent run directory).")
    for flag, help_txt in (
        ("--resume", "Resume the most recent parent run pointed to by runs/latest_<epithet>, skipping completed stages."),
        ("--force", "Re-execute stages even when a completed triage_queue.csv checkpoint already exists."),
        ("--gold-only", "Propagate --gold-only across all stages."),
        ("--skip-harvest", "Reuse canonical Step 01 vouchers and Step 02 segments in Stage 1 without GBIF harvesting."),
        ("--skip-discovery", "Skip Stage 1 when an Eastern NA discovery baseline has already been segmented."),
        ("--summarize-only", "Skip pipeline execution and compile the cross-stage summary from existing stage folders."),
        ("--dry-run", "Print the six-stage run_pipeline.sh commands and parent directory without executing or writing files."),
    ):
        parser.add_argument(flag, action="store_true", default=False, help=help_txt)
    return parser


def main(argv: Optional[Sequence[str]] = None, workspace_root: Optional[Path] = None) -> int:
    raw_argv = list(argv) if argv is not None else sys.argv[1:]
    args = build_parser().parse_args(argv)
    ws = (workspace_root or PROJECT_ROOT).resolve()
    focal_std = standardize_packera_taxon(args.focal_taxon) or args.focal_taxon.strip()

    parent_run_id, unique_id, stage_prefix, parent_out_dir = resolve_focal_parent_run_id(
        workspace_root=ws,
        focal_taxon=focal_std,
        prefix=args.prefix,
        parent_run_id=args.parent_run_id,
        resume=args.resume,
    )
    run_script = ws / "run_pipeline.sh"

    stages = build_six_stage_commands(
        run_script=run_script,
        focal_taxon=focal_std,
        prefix=stage_prefix,
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
        print(f"[Dry-Run Parent Directory] runs/{stage_prefix.rstrip('/')}")
        for st_label, rid, cmd in stages:
            print(f"[{st_label}] ({rid}): {' '.join(cmd)}")
        return 0

    parent_out_dir.mkdir(parents=True, exist_ok=True)
    if not args.summarize_only:
        update_latest_focal_symlink(ws, focal_std, parent_run_id)

    stage_statuses: Dict[str, str] = {}
    if args.skip_discovery:
        r1_eff = _six_stage_run_ids(stage_prefix, discovery_run_id=args.discovery_run_id)[0]
        stage_statuses[r1_eff] = "skipped_discovery_flag"

    if not args.summarize_only:
        for st_label, rid, cmd in stages:
            if not args.force and stage_is_complete_focal(ws, rid):
                print(f"[Checkpoint Skip] {st_label} ({rid}) already has triage_queue.csv.", flush=True)
                stage_statuses[rid] = "skipped_checkpoint"
                continue
            print(f"\n=== Executing {st_label} ({rid}) ===", flush=True)
            proc = subprocess.run(cmd, cwd=str(ws), check=False)
            stage_statuses[rid] = "completed" if proc.returncode == 0 else f"failed_{proc.returncode}"
            if proc.returncode != 0:
                print(
                    f"[ERROR] {st_label} ({rid}) exited with non-zero status {proc.returncode}.",
                    file=sys.stderr,
                )
                return proc.returncode
    else:
        stage_statuses.update({rid: "summarize_only" for _, rid, _ in stages})

    csv_out, md_out = compile_complex_resolution_summary(
        focal_taxon=focal_std,
        prefix=stage_prefix,
        discovery_run_id=args.discovery_run_id,
        workspace_root=ws,
        db_path=args.db,
        output_dir=parent_out_dir,
    )
    r1, r2, r3, r4, r5, r6 = _six_stage_run_ids(stage_prefix, discovery_run_id=args.discovery_run_id)
    manifest_out = write_focal_run_manifest(
        output_dir=parent_out_dir,
        parent_run_id=parent_run_id,
        unique_id=unique_id,
        focal_taxon=focal_std,
        cli_arguments=raw_argv,
        stage_run_ids={
            "stage_1_discovery": r1,
            "stage_2_complex_triage": r2,
            "stage_3_intraspecific_clean": r3,
            "stage_4_spatial_thin_field_audit": r4,
            "stage_5_monographer_gold_audit": r5,
            "stage_6_lineage_diagnosability": r6,
        },
        stage_statuses=stage_statuses,
        workspace_root=ws,
    )
    if args.output_dir is not None and args.output_dir.resolve() != parent_out_dir.resolve():
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for src in (csv_out, md_out, manifest_out):
            shutil.copy2(src, args.output_dir / src.name)
    print(f"[Complex Resolution Summary] Parent Run Dir  : {parent_out_dir}")
    print(f"[Complex Resolution Summary] Comparison CSV  : {csv_out}")
    print(f"[Complex Resolution Summary] Synthesis Report: {md_out}")
    print(f"[Complex Resolution Summary] Run Manifest    : {manifest_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
