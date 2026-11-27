#!/usr/bin/env python3
"""
scripts/pipeline/step07_focal_stage_runners.py
==============================================
Six-stage (and backward-compatible four-stage) command builders plus cross-stage
comparison CSV and Markdown synthesis compilers for focal Packera species complexes.
"""

from __future__ import annotations

import csv
import os
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon
from scripts.pipeline.step07_focal_cohort_utils import (
    _extract_stage_metrics,
    _stage_run_ids,
    slugify_taxon,
)
from scripts.pipeline.step07_focal_parent_run import _six_stage_run_ids
from scripts.pipeline.step07_focal_report_writer import render_complex_resolution_markdown


def _common_stage_tail(db_path: Optional[Path], threads: Optional[int]) -> List[str]:
    tail: List[str] = []
    if db_path is not None:
        tail.extend(["--db", str(db_path)])
    if threads is not None and threads > 0:
        tail.extend(["--threads", str(threads)])
    return tail


def build_six_stage_commands(
    *,
    run_script: Path,
    focal_taxon: str,
    prefix: str,
    discovery_run_id: Optional[str] = None,
    geographic_scope: str = "eastern_na",
    gold_anchors: int = 30,
    spatial_thin_km: float = 15.0,
    threads: Optional[int] = None,
    gpu_batch_size: Optional[int] = None,
    db_path: Optional[Path] = None,
    skip_discovery: bool = False,
    skip_harvest: bool = False,
    gold_only: bool = False,
) -> List[Tuple[str, str, List[str]]]:
    """Return ordered `(stage_label, run_id, argv)` tuples for the six-stage resolution funnel."""
    common_tail = _common_stage_tail(db_path, threads)
    gold_flag = ["--gold-only"] if gold_only else []
    r1, r2, r3, r4, r5, r6 = _six_stage_run_ids(prefix, discovery_run_id=discovery_run_id)

    stages: List[Tuple[str, str, List[str]]] = []
    if not skip_discovery:
        cmd1 = [
            "bash",
            str(run_script),
            "--run-id",
            r1,
            *(["--skip-harvest"] if skip_harvest else []),
            "--geographic-scope",
            geographic_scope,
            *gold_flag,
            *common_tail,
        ]
        if gpu_batch_size is not None and gpu_batch_size > 0:
            cmd1.extend(["--gpu-batch-size", str(gpu_batch_size)])
        stages.append(("Stage_1_Discovery", r1, cmd1))

    cmd2 = [
        "bash",
        str(run_script),
        "--run-id",
        r2,
        "--skip-harvest",
        "--geographic-scope",
        geographic_scope,
        "--focal-taxon",
        focal_taxon,
        "--discover-complexes",
        "--triage-source",
        r1,
        "--gold-anchors",
        str(gold_anchors),
        *gold_flag,
        "--skip-backbone-download",
        *common_tail,
    ]
    stages.append(("Stage_2_Complex_Triage", r2, cmd2))

    cmd3 = [
        "bash",
        str(run_script),
        "--run-id",
        r3,
        "--skip-harvest",
        "--geographic-scope",
        geographic_scope,
        "--focal-taxon",
        focal_taxon,
        "--intraspecific-mode",
        "--exclude-misidentified",
        "--triage-source",
        r2,
        *gold_flag,
        "--skip-backbone-download",
        *common_tail,
    ]
    stages.append(("Stage_3_Intraspecific_Clean", r3, cmd3))

    cmd4 = [
        "bash",
        str(run_script),
        "--run-id",
        r4,
        "--start-step",
        "3",
        "--geographic-scope",
        geographic_scope,
        "--focal-taxon",
        focal_taxon,
        "--intraspecific-mode",
        "--exclude-misidentified",
        "--triage-source",
        r2,
        "--spatial-thin-km",
        f"{spatial_thin_km:g}",
        *(["--gold-only"] if gold_only else ["--determiner-tier", "gold_silver"]),
        "--skip-backbone-download",
        *common_tail,
    ]
    stages.append(("Stage_4_Spatial_Thin_Field_Audit", r4, cmd4))

    cmd5 = [
        "bash",
        str(run_script),
        "--run-id",
        r5,
        "--start-step",
        "3",
        "--geographic-scope",
        geographic_scope,
        "--focal-taxon",
        focal_taxon,
        "--intraspecific-mode",
        "--exclude-misidentified",
        "--triage-source",
        r2,
        "--spatial-thin-km",
        f"{spatial_thin_km:g}",
        "--gold-only",
        "--skip-backbone-download",
        *common_tail,
    ]
    stages.append(("Stage_5_Monographer_Gold_Audit", r5, cmd5))

    cmd6 = [
        "bash",
        str(run_script),
        "--run-id",
        r6,
        "--start-step",
        "7",
        "--geographic-scope",
        geographic_scope,
        "--focal-taxon",
        focal_taxon,
        "--lineage-diagnosability-source",
        r3,
        "--skip-backbone-download",
        *common_tail,
    ]
    stages.append(("Stage_6_Lineage_Diagnosability", r6, cmd6))
    return stages


def build_four_stage_commands(
    *,
    run_script: Path,
    focal_taxon: str,
    prefix: str,
    discovery_run_id: Optional[str] = None,
    geographic_scope: str = "eastern_na",
    gold_anchors: int = 30,
    spatial_thin_km: float = 15.0,
    threads: Optional[int] = None,
    gpu_batch_size: Optional[int] = None,
    db_path: Optional[Path] = None,
    skip_discovery: bool = False,
    skip_harvest: bool = False,
    gold_only: bool = False,
) -> List[Tuple[str, str, List[str]]]:
    """Backward-compatible four-stage builder preserved for legacy 4-stage callers."""
    six = build_six_stage_commands(
        run_script=run_script,
        focal_taxon=focal_taxon,
        prefix=prefix,
        discovery_run_id=discovery_run_id,
        geographic_scope=geographic_scope,
        gold_anchors=gold_anchors,
        spatial_thin_km=spatial_thin_km,
        threads=threads,
        gpu_batch_size=gpu_batch_size,
        db_path=db_path,
        skip_discovery=skip_discovery,
        skip_harvest=skip_harvest,
        gold_only=gold_only,
    )
    _, _, _, r4_leg = _stage_run_ids(prefix, discovery_run_id=discovery_run_id)
    head = [st for st in six if st[0] in {"Stage_1_Discovery", "Stage_2_Complex_Triage", "Stage_3_Intraspecific_Clean"}]
    s5_cmd = list(next(st[2] for st in six if st[0] == "Stage_5_Monographer_Gold_Audit"))
    rid_idx = s5_cmd.index("--run-id") + 1
    s5_cmd[rid_idx] = r4_leg
    return [*head, ("Stage_4_Spatially_Thinned", r4_leg, s5_cmd)]


def _row_otu_key(row: Dict[str, Any]) -> Tuple[str, str]:
    return (
        str(row.get("catalogNumber") or "").strip(),
        str(row.get("plant_individual_id") or "0").strip(),
    )


def _resolve_active_stage_specs(
    prefix: str, discovery_run_id: Optional[str], runs_root: Path
) -> List[Tuple[str, str]]:
    r1, r2, r3, r4, r5, r6 = _six_stage_run_ids(prefix, discovery_run_id=discovery_run_id)
    _, _, _, r4_legacy = _stage_run_ids(prefix, discovery_run_id=discovery_run_id)
    has_six = any((runs_root / rid).exists() for rid in (r4, r5, r6))
    has_legacy_only = not has_six and (runs_root / r4_legacy).exists()
    if has_legacy_only:
        return [
            ("Stage_1_Discovery", r1),
            ("Stage_2_Complex_Triage", r2),
            ("Stage_3_Intraspecific_Clean", r3),
            ("Stage_4_Spatially_Thinned", r4_legacy),
        ]
    return [
        ("Stage_1_Discovery", r1),
        ("Stage_2_Complex_Triage", r2),
        ("Stage_3_Intraspecific_Clean", r3),
        ("Stage_4_Spatial_Thin_Field_Audit", r4),
        ("Stage_5_Monographer_Gold_Audit", r5),
        ("Stage_6_Lineage_Diagnosability", r6),
    ]


def compile_complex_resolution_summary(
    *,
    focal_taxon: str,
    prefix: str,
    discovery_run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    """Compile cross-run comparison CSVs and Markdown report across all active stages."""
    focal_std = standardize_packera_taxon(focal_taxon) or focal_taxon.strip()
    slug = slugify_taxon(focal_std)
    runs_root = workspace_root / "runs"
    out_dir = output_dir or (workspace_root / "outputs" / "reports")
    out_dir.mkdir(parents=True, exist_ok=True)

    if db_path is not None:
        resolved_db = db_path
    else:
        cand_new = workspace_root / "data" / "plantmatics.sqlite"
        cand_old = workspace_root / "data" / "packera_pipeline.sqlite"
        resolved_db = cand_new if cand_new.is_file() else (cand_old if cand_old.is_file() else cand_new)
    conn: Optional[sqlite3.Connection] = None
    if resolved_db.is_file():
        try:
            conn = sqlite3.connect(f"file:{resolved_db}?mode=ro", uri=True)
        except sqlite3.Error:
            conn = None

    stage_specs = _resolve_active_stage_specs(prefix, discovery_run_id, runs_root)
    summary_rows: List[Dict[str, Any]] = []
    stage_triage_map: Dict[str, List[Dict[str, Any]]] = {}
    try:
        for st_name, rid in stage_specs:
            m, eval_rows = _extract_stage_metrics(
                stage_name=st_name,
                run_id=rid,
                focal_std=focal_std,
                runs_root=runs_root,
                conn=conn,
            )
            summary_rows.append(m)
            stage_triage_map[st_name] = eval_rows
    finally:
        if conn is not None:
            conn.close()

    clean_rows = stage_triage_map.get("Stage_3_Intraspecific_Clean", [])
    s2_eval_rows = stage_triage_map.get("Stage_2_Complex_Triage", [])
    misid_calls = {"Misidentification_Reassign", "Misidentification_Severe"}
    s2_misid_keys = {
        _row_otu_key(r)
        for r in s2_eval_rows
        if str(r.get("taxonomic_treatment_subtype", "")).strip() in misid_calls
        or str(r.get("taxonomic_status_call", "")).strip() in misid_calls
    }
    s2_all_keys = {_row_otu_key(r) for r in s2_eval_rows}
    s2_misid_sheets = {c for c, _ in s2_misid_keys if c}
    s3_all_keys = {_row_otu_key(r) for r in clean_rows}
    s3_sheets = {c for c, _ in s3_all_keys if c}
    if clean_rows:
        n_purged_before_stage3_otus = len(s2_misid_keys - s3_all_keys)
        n_purged_before_stage3_sheets = len(s2_misid_sheets - s3_sheets)
    else:
        s2_non_misid_sheets = {c for c, _ in (s2_all_keys - s2_misid_keys) if c}
        n_purged_before_stage3_otus = len(s2_misid_keys)
        n_purged_before_stage3_sheets = len(s2_misid_sheets - s2_non_misid_sheets)

    for r_sum in summary_rows:
        if r_sum.get("stage") in {"Stage_2_Complex_Triage", "Stage_3_Intraspecific_Clean"}:
            r_sum["n_purged_before_stage3_otus"] = n_purged_before_stage3_otus
            r_sum["n_purged_before_stage3_sheets"] = n_purged_before_stage3_sheets
        else:
            r_sum["n_purged_before_stage3_otus"] = ""
            r_sum["n_purged_before_stage3_sheets"] = ""

    fieldnames = list(summary_rows[0].keys())
    csv_primary = out_dir / f"{slug}_Stage_Comparison.csv"
    csv_targets = [
        csv_primary,
        out_dir / f"{slug}_Six_Stage_Comparison.csv",
        out_dir / f"{slug}_Four_Stage_Comparison.csv",
        out_dir / f"{slug}_complex_stage_comparison.csv",
    ]
    for target_csv in csv_targets:
        tmp_csv = target_csv.with_suffix(".csv.tmp")
        with open(tmp_csv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_rows)
            fh.flush()
            os.fsync(fh.fileno())
        tmp_csv.replace(target_csv)

    md_out = render_complex_resolution_markdown(
        focal_std=focal_std,
        summary_rows=summary_rows,
        stage_triage_map=stage_triage_map,
        runs_root=runs_root,
        n_purged_before_stage3_otus=n_purged_before_stage3_otus,
        n_purged_before_stage3_sheets=n_purged_before_stage3_sheets,
        n_s2_misid_otus=len(s2_misid_keys),
        n_s2_misid_sheets=len(s2_misid_sheets),
        out_md_path=out_dir / f"{slug}_Complex_Resolution_Report.md",
    )
    return csv_primary, md_out
