#!/usr/bin/env python3
"""
scripts/pipeline/step07_focal_stage_runners.py
==============================================
Four-stage command builders and cross-stage comparison CSV / taxonomic
adjudication Markdown report compilers for focal Packera species complexes.
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
    NA_TOKENS,
    _extract_stage_metrics,
    _stage_run_ids,
    slugify_taxon,
)


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
    """Return ordered (stage_label, run_id, argv) tuples for the four-stage funnel."""
    common_tail: List[str] = []
    if db_path is not None:
        common_tail.extend(["--db", str(db_path)])
    if threads is not None and threads > 0:
        common_tail.extend(["--threads", str(threads)])
    gold_flag = ["--gold-only"] if gold_only else []

    r1, r2, r3, r4 = _stage_run_ids(prefix, discovery_run_id=discovery_run_id)

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
        "--gold-only",
        "--skip-backbone-download",
        *common_tail,
    ]
    stages.append(("Stage_4_Spatially_Thinned", r4, cmd4))
    return stages


def _nz_token(val: Any) -> str:
    s = str(val or "").strip()
    return "" if s.upper() in NA_TOKENS else s


def _row_otu_key(row: Dict[str, Any]) -> Tuple[str, str]:
    return (
        str(row.get("catalogNumber") or "").strip(),
        str(row.get("plant_individual_id") or "0").strip(),
    )


def compile_complex_resolution_summary(
    *,
    focal_taxon: str,
    prefix: str,
    discovery_run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    """Compile cross-run comparison CSV and Markdown synthesis report across Stages 1 to 4."""
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

    r1, r2, r3, r4 = _stage_run_ids(prefix, discovery_run_id=discovery_run_id)
    stage_specs = [
        ("Stage_1_Discovery", r1),
        ("Stage_2_Complex_Triage", r2),
        ("Stage_3_Intraspecific_Clean", r3),
        ("Stage_4_Spatially_Thinned", r4),
    ]

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

    s1, s2, s3, s4 = summary_rows
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
    csv_out = out_dir / f"{slug}_Four_Stage_Comparison.csv"
    for target_csv in (csv_out, out_dir / f"{slug}_complex_stage_comparison.csv"):
        tmp_csv = target_csv.with_suffix(".csv.tmp")
        with open(tmp_csv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_rows)
            fh.flush()
            os.fsync(fh.fileno())
        tmp_csv.replace(target_csv)

    type_controls = [
        r
        for r in clean_rows
        if str(r.get("taxonomic_treatment_subtype", "")).strip() in {"Confirmed_Type_Control", "Type_Synonymy_Candidate"}
        or str(r.get("taxonomic_status_call", "")).strip() in {"Confirmed_Type_Control", "Type_Synonymy_Candidate"}
    ]
    pull_subtypes = {
        "Putative_Cryptic_Species",
        "Hybrid_Introgressant",
        "Subspecies_Ecotype",
        "Type_Synonymy_Candidate",
        "Misidentification_Reassign",
        "Misidentification_Severe",
    }
    priority_pulls = [
        r
        for r in clean_rows
        if str(r.get("taxonomic_treatment_subtype", "")).strip() in pull_subtypes
        or str(r.get("taxonomic_status_call", "")).strip() in misid_calls
        or str(r.get("synthesis_triage_priority", "")).strip() == "CRITICAL"
    ]
    clean_keys = set(s3_all_keys)
    for r in s2_eval_rows:
        sub_s2 = str(r.get("taxonomic_treatment_subtype", "")).strip()
        call_s2 = str(r.get("taxonomic_status_call", "")).strip()
        key_s2 = _row_otu_key(r)
        if (sub_s2 in misid_calls or call_s2 in misid_calls) and key_s2 not in clean_keys:
            r_copy = dict(r)
            r_copy["_purged_stage2"] = True
            priority_pulls.append(r_copy)
            clean_keys.add(key_s2)

    prio_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 3, "LOW": 4, "RESOLVED": 5}
    priority_pulls.sort(
        key=lambda r: (
            1 if r.get("_purged_stage2") else prio_rank.get(str(r.get("synthesis_triage_priority", "")).strip(), 4),
            0 if str(r.get("taxonomic_treatment_subtype", "")).strip() != "Subspecies_Ecotype" else 2,
            str(r.get("catalogNumber", "")),
        )
    )

    eco_counts: Dict[str, int] = {}
    disj_counts: Dict[str, int] = {}
    n_ungeoref_clean = 0
    for r_c in clean_rows:
        er_v = _nz_token(r_c.get("ecotype_regime"))
        bs_v = _nz_token(r_c.get("biogeographic_sector"))
        if er_v.startswith("Ecotype_"):
            eco_counts[er_v] = eco_counts.get(er_v, 0) + 1
        elif not er_v:
            n_ungeoref_clean += 1
        if bs_v and bs_v != "Continuous_Range":
            disj_counts[bs_v] = disj_counts.get(bs_v, 0) + 1

    if int(s3["n_cohesive_cryptic_otus"]) > 0:
        verdict = (
            f"We detected a candidate cryptic species inside *{focal_std}* ({s3['n_cohesive_cryptic_otus']} purified "
            f"plants in shape K = {s3['gmm_optimal_k'] or 'N/A'}) that separates cleanly from the type baseline."
        )
    elif int(s3["n_subspecies_ecotype"]) > 0:
        verdict = (
            f"*{focal_std}* holds together as one species (0 cryptic species), while splitting internally into "
            f"well-supported habitat or geographic varieties ({s3['n_subspecies_ecotype']} plants flagged at infraspecific rank; "
            f"shape K = {s3['gmm_optimal_k'] or 1})."
        )
    else:
        verdict = (
            f"*{focal_std}* resolves as one cohesive species after we removed {n_purged_before_stage3_otus} misidentified "
            f"plants ({n_purged_before_stage3_sheets} full sheets) before Stage 3 (shape K = {s3['gmm_optimal_k'] or 1}; "
            f"cryptic species = 0)."
        )

    eco_line = ""
    if eco_counts or disj_counts:
        parts_e = [f"`{k}` ({v} GPS-verified plants)" for k, v in sorted(eco_counts.items())]
        parts_e.extend(f"`{k}` ({v} plants)" for k, v in sorted(disj_counts.items()))
        eco_line = (
            f"* **Habitat and Geographic Groups (Stage 3):** {'; '.join(parts_e)}. "
            f"{n_ungeoref_clean} sheets without GPS coordinates were kept unassigned for 30 m soil modeling so variable label handwriting never biases habitat groups."
        )

    thin_note = (
        f"Stage 4 habitat-stratified 15 km thinning retained {s4['n_focal_otus']} *{focal_std}* plants "
        f"(shape K = {s4['gmm_optimal_k'] or 'N/A'}; varieties/ecotypes = {s4['n_subspecies_ecotype']}; "
        f"cryptic species = {s4['n_cohesive_cryptic_otus']})."
    )
    has_cyto_line = bool(
        s3.get("reported_cytotypes")
        or s2.get("reported_cytotypes")
        or s1.get("reported_cytotypes")
        or max(int(s.get("n_cytovouchers_linked") or 0) for s in (s1, s2, s3, s4)) > 0
    )

    md_lines = [
        f"# Species Complex Summary Report: *{focal_std}*",
        "",
        "## 1. Plain-Language Taxonomic Summary",
        "",
        f"* **Bottom-Line Call:** {verdict}",
        f"* **Compared Complex Neighbors (Stages 1 and 2):** {s2['complex_neighbors'] or s1['complex_neighbors'] or 'None above affinity threshold'}",
        (
            f"* **Misidentified Sheets Removed Before Stage 3:** "
            f"{n_purged_before_stage3_otus} plants across {n_purged_before_stage3_sheets} full herbarium sheets "
            f"({len(s2_misid_keys)} Stage 2 misidentifications flagged across {len(s2_misid_sheets)} sheets)."
        ),
        *([eco_line] if eco_line else []),
        (
            f"* **Contact-Zone Hybrids (Stage 2):** {s2['n_hybrid_introgressant']} plants ({s2['hybrid_parent_pairs']})."
            if s2.get("hybrid_parent_pairs")
            else f"* **Contact-Zone Hybrids (Stage 2):** {s2['n_hybrid_introgressant']} plants."
        ),
        f"* **Spatial Thinning Check (Stage 4):** {thin_note}",
        *(
            [
                f"* **Reported Ploidy and Cytovouchers:** "
                f"{(s3['reported_cytotypes'] or s2['reported_cytotypes'] or s1['reported_cytotypes'] or 'Unassigned')} "
                f"(linked cytovouchers in Stage 3: {s3['n_cytovouchers_linked']})."
            ]
            if has_cyto_line
            else []
        ),
        "",
        "## 2. Four-Stage Numbers at a Glance",
        "",
        "| Stage | Focal Plants | Shape K (Raw to Final) | Soil CV Accuracy | Misid. Removed | Hybrids | Varieties or Ecotypes | Plastic Outliers | Cryptic Species | Soil Overlap (D) | Flowering Shift (Days) |",
        "| :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in summary_rows:
        md_lines.append(
            f"| {r['stage']} | {r['n_focal_otus']} | {r['gmm_raw_k'] or 'N/A'} to {r['gmm_optimal_k'] or 'N/A'} | "
            f"{r['spatial_block_cv_acc'] or 'N/A'} | {r['n_misidentification_reassign']} | "
            f"{r['n_hybrid_introgressant']} | {r['n_subspecies_ecotype']} | "
            f"{r['n_ecophenotypic_plasticity']} | {r['n_cohesive_cryptic_otus']} | "
            f"{r['min_d_edaph'] or 'N/A'} | {r['max_abs_delta_doy_days'] or 'N/A'} |"
        )

    md_lines.extend(["", "## 3. Type Specimens and Naming Anchors", ""])
    if type_controls:
        for tr in type_controls[:10]:
            cat = tr.get("catalogNumber", "Unknown")
            inst = tr.get("institutionCode", "")
            cl = _nz_token(tr.get("gmm_cluster")) or "unclustered"
            t_sub = _nz_token(tr.get("taxonomic_treatment_subtype")) or "Confirmed_Type_Control"
            bas = (
                _nz_token(tr.get("nearest_type_basionym"))
                or _nz_token(tr.get("typified_basionym"))
                or _nz_token(tr.get("backbone_basionym_name"))
                or _nz_token(tr.get("basionym_name"))
                or _nz_token(tr.get("species_raw"))
                or focal_std
            )
            md_lines.append(f"* `{cat}` ({inst}): basionym *{bas}* in shape group `{cl}` (`{t_sub}`).")
        if focal_std == "Packera dubia":
            md_lines.append(
                "* Primary *locus classicus* note: Andre Michaux's original *Senecio tomentosus* Michx. specimen from Flat Rock "
                "(Forty Acre Rock area, Lancaster and Kershaw Counties, South Carolina; `P-MICH`, Paris) awaits digital imaging, "
                "so verified Forty Acre Rock granite flatrock topotypes anchor typical *Packera dubia* while NYBG *Senecio alabamensis* "
                "Britton types (`NY00232445`, `NY00232447`) anchor the Coastal Plain sand lineage."
            )
    else:
        md_lines.append("* No explicit type specimen was present in Stage 3; verify physical type images before formal naming changes.")

    md_lines.extend(["", "## 4. Herbarium Sheets to Review First", ""])
    if priority_pulls:
        for pr in priority_pulls[:25]:
            cat = pr.get("catalogNumber", "Unknown")
            pid = pr.get("plant_individual_id", 0)
            inst = pr.get("institutionCode", "")
            sub = pr.get("taxonomic_treatment_subtype", "")
            rec = pr.get("recommended_determination", "")
            rat = str(pr.get("synthesis_rationale", "")).strip()
            s2_tag = " [Stage 2 purged]" if pr.get("_purged_stage2") else ""
            md_lines.append(f"* `{cat}` (plant `{pid}`, {inst}){s2_tag}: `{sub}` to *{rec}*. {rat}")
    else:
        md_lines.append("* Zero non-confirming vouchers remained in Stage 3 after Stage 2 cleanup.")

    md_out = out_dir / f"{slug}_Complex_Resolution_Report.md"
    tmp_md = md_out.with_suffix(".md.tmp")
    with open(tmp_md, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(md_lines) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    tmp_md.replace(md_out)
    return csv_out, md_out
