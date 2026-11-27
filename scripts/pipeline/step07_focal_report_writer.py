#!/usr/bin/env python3
"""
scripts/pipeline/step07_focal_report_writer.py
==============================================
Plain-language Markdown adjudication report compiler for the six-stage focal
species-complex workflow, integrating Stage 4/5 unconfounded sensitivity audits
and Stage 6 supervised candidate-lineage diagnosability and key-trait outputs.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from scripts.pipeline.step07_focal_cohort_utils import NA_TOKENS, _read_csv_rows


def _nz(val: Any) -> str:
    s = str(val or "").strip()
    return "" if s.upper() in NA_TOKENS else s


def _row_otu_key_rep(row: Mapping[str, Any]) -> tuple[str, str]:
    return (str(row.get("catalogNumber") or "").strip(), str(row.get("plant_individual_id") or "0").strip())


def render_complex_resolution_markdown(
    *, focal_std: str, summary_rows: Sequence[Mapping[str, Any]], stage_triage_map: Mapping[str, List[Dict[str, Any]]],
    runs_root: Path, n_purged_before_stage3_otus: int, n_purged_before_stage3_sheets: int,
    n_s2_misid_otus: int, n_s2_misid_sheets: int, out_md_path: Path,
) -> Path:
    by_stage: Dict[str, Mapping[str, Any]] = {str(r["stage"]): r for r in summary_rows}
    s1 = by_stage.get("Stage_1_Discovery", summary_rows[0])
    s2 = by_stage.get("Stage_2_Complex_Triage", summary_rows[min(1, len(summary_rows) - 1)])
    s3 = by_stage.get("Stage_3_Intraspecific_Clean", summary_rows[min(2, len(summary_rows) - 1)])
    s4 = by_stage.get("Stage_4_Spatial_Thin_Field_Audit") or by_stage.get("Stage_4_Spatially_Thinned")
    s5 = by_stage.get("Stage_5_Monographer_Gold_Audit")
    s6 = by_stage.get("Stage_6_Lineage_Diagnosability")

    clean_rows = stage_triage_map.get("Stage_3_Intraspecific_Clean", [])
    s2_eval_rows = stage_triage_map.get("Stage_2_Complex_Triage", [])
    misid_calls = {"Misidentification_Reassign", "Misidentification_Severe"}
    s3_all_keys = {_row_otu_key_rep(r) for r in clean_rows}
    type_controls = [
        r for r in clean_rows
        if _nz(r.get("taxonomic_treatment_subtype")) in {"Confirmed_Type_Control", "Type_Synonymy_Candidate"}
        or _nz(r.get("taxonomic_status_call")) in {"Confirmed_Type_Control", "Type_Synonymy_Candidate"}
    ]
    pull_subtypes = {"Putative_Cryptic_Species", "Hybrid_Introgressant", "Subspecies_Ecotype",
                     "Type_Synonymy_Candidate", "Misidentification_Reassign", "Misidentification_Severe"}
    priority_pulls: List[Dict[str, Any]] = [
        dict(r) for r in clean_rows
        if _nz(r.get("taxonomic_treatment_subtype")) in pull_subtypes
        or _nz(r.get("taxonomic_status_call")) in misid_calls
        or _nz(r.get("synthesis_triage_priority")) == "CRITICAL"
    ]
    clean_keys = set(s3_all_keys)
    for r in s2_eval_rows:
        key_s2 = _row_otu_key_rep(r)
        if (_nz(r.get("taxonomic_treatment_subtype")) in misid_calls or _nz(r.get("taxonomic_status_call")) in misid_calls) and key_s2 not in clean_keys:
            priority_pulls.append({**r, "_purged_stage2": True})
            clean_keys.add(key_s2)

    prio_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 3, "LOW": 4, "RESOLVED": 5}
    priority_pulls.sort(key=lambda r: (
        1 if r.get("_purged_stage2") else prio_rank.get(_nz(r.get("synthesis_triage_priority")), 4),
        0 if _nz(r.get("taxonomic_treatment_subtype")) != "Subspecies_Ecotype" else 2,
        str(r.get("catalogNumber", "")),
    ))

    eco_counts: Dict[str, int] = {}
    disj_counts: Dict[str, int] = {}
    n_ungeoref_clean = 0
    for r_c in clean_rows:
        er_v = _nz(r_c.get("ecotype_regime"))
        bs_v = _nz(r_c.get("biogeographic_sector"))
        if er_v.startswith("Ecotype_"):
            eco_counts[er_v] = eco_counts.get(er_v, 0) + 1
        elif not er_v:
            n_ungeoref_clean += 1
        if bs_v and bs_v != "Continuous_Range":
            disj_counts[bs_v] = disj_counts.get(bs_v, 0) + 1

    if int(s3.get("n_cohesive_cryptic_otus") or 0) > 0:
        verdict = (
            f"We detected a candidate cryptic species inside *{focal_std}* ({s3['n_cohesive_cryptic_otus']} purified "
            f"plants in shape K = {s3.get('gmm_optimal_k') or 'N/A'}) that separates cleanly from the type baseline."
        )
    elif int(s3.get("n_subspecies_ecotype") or 0) > 0:
        verdict = (
            f"*{focal_std}* holds together as one species (0 cryptic species), while splitting internally into "
            f"well-supported habitat or geographic lineages ({s3['n_subspecies_ecotype']} plants flagged at infraspecific rank; "
            f"shape K = {s3.get('gmm_optimal_k') or 1})."
        )
    else:
        verdict = (
            f"*{focal_std}* resolves as one cohesive species after we removed {n_purged_before_stage3_otus} misidentified "
            f"plants ({n_purged_before_stage3_sheets} full sheets) before Stage 3 (shape K = {s3.get('gmm_optimal_k') or 1}; "
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

    thin_lines: List[str] = []
    if s4 is not None and s5 is not None:
        thin_lines.append(
            f"* **Field GPS Spatial Thinning Audit (Stage 4, Gold and Silver at 15 km):** "
            f"Retained {s4['n_focal_otus']} *{focal_std}* plants (shape K = {s4.get('gmm_optimal_k') or 'N/A'}; "
            f"varieties or ecotypes = {s4['n_subspecies_ecotype']}; soil overlap D = {s4.get('min_d_edaph') or 'N/A'}; "
            f"cryptic species = {s4['n_cohesive_cryptic_otus']})."
        )
        thin_lines.append(
            f"* **Strict Monographer Baseline (Stage 5, Gold-Only at 15 km):** "
            f"Retained {s5['n_focal_otus']} *{focal_std}* plants (shape K = {s5.get('gmm_optimal_k') or 'N/A'}; "
            f"varieties or ecotypes = {s5['n_subspecies_ecotype']}; soil overlap D = {s5.get('min_d_edaph') or 'N/A'}; "
            f"cryptic species = {s5['n_cohesive_cryptic_otus']})."
        )
    elif s4 is not None:
        thin_lines.append(
            f"* **Spatial Thinning Check (Stage 4):** Stage 4 habitat-stratified 15 km thinning retained {s4['n_focal_otus']} "
            f"*{focal_std}* plants (shape K = {s4.get('gmm_optimal_k') or 'N/A'}; varieties/ecotypes = {s4['n_subspecies_ecotype']}; "
            f"cryptic species = {s4['n_cohesive_cryptic_otus']})."
        )

    s6_diag_rows: List[Dict[str, str]] = []
    s6_pair_rows: List[Dict[str, str]] = []
    s6_key_rows: List[Dict[str, str]] = []
    if s6 is not None and _nz(s6.get("run_id")):
        s6_rep_dir = runs_root / str(s6["run_id"]) / "outputs" / "reports"
        s6_diag_rows = _read_csv_rows(s6_rep_dir / "candidate_lineage_diagnosability.csv")
        s6_pair_rows = _read_csv_rows(s6_rep_dir / "candidate_lineage_pairwise_tests.csv")
        s6_key_rows = _read_csv_rows(s6_rep_dir / "dichotomous_key_character_table.csv")

    s6_bullet = ""
    if s6_diag_rows:
        oof_acc = _nz(s6.get("oof_lineage_morph_acc")) or _nz(s6_diag_rows[0].get("global_oof_morph_accuracy")) or "N/A"
        lin_calls = "; ".join(
            f"`{r.get('candidate_lineage', '')}` -> `{r.get('recommended_rank', '')}` (held-out morph recall = {r.get('oof_morph_recall', 'N/A')})"
            for r in s6_diag_rows
            if _nz(r.get("candidate_lineage"))
        )
        s6_bullet = (
            f"* **Candidate Lineage Diagnosability (Stage 6, 75% Rule Check):** "
            f"Global spatially blocked morphological accuracy = {oof_acc}. {lin_calls}."
        )

    all_stages_chk = [r for r in summary_rows]
    has_cyto_line = bool(
        any(_nz(r.get("reported_cytotypes")) for r in all_stages_chk)
        or max((int(r.get("n_cytovouchers_linked") or 0) for r in all_stages_chk), default=0) > 0
    )
    rep_cyto = next((_nz(r.get("reported_cytotypes")) for r in (s3, s2, s1) if _nz(r.get("reported_cytotypes"))), "Unassigned")

    n_st_word = "Six-Stage" if len(summary_rows) >= 6 else "Stage"
    md_lines: List[str] = [
        f"# Species Complex Summary Report: *{focal_std}*",
        "",
        "## 1. Plain-Language Taxonomic Summary",
        "",
        f"* **Bottom-Line Call:** {verdict}",
        f"* **Compared Complex Neighbors (Stages 1 and 2):** {s2.get('complex_neighbors') or s1.get('complex_neighbors') or 'None above affinity threshold'}",
        (
            f"* **Misidentified Sheets Removed Before Stage 3:** "
            f"{n_purged_before_stage3_otus} plants across {n_purged_before_stage3_sheets} full herbarium sheets "
            f"({n_s2_misid_otus} Stage 2 misidentifications flagged across {n_s2_misid_sheets} sheets)."
        ),
        *([eco_line] if eco_line else []),
        (
            f"* **Contact-Zone Hybrids (Stage 2):** {s2['n_hybrid_introgressant']} plants ({s2['hybrid_parent_pairs']})."
            if _nz(s2.get("hybrid_parent_pairs"))
            else f"* **Contact-Zone Hybrids (Stage 2):** {s2['n_hybrid_introgressant']} plants."
        ),
        *thin_lines,
        *([s6_bullet] if s6_bullet else []),
        *(
            [
                f"* **Reported Ploidy and Cytovouchers:** "
                f"{rep_cyto} (linked cytovouchers in Stage 3: {s3.get('n_cytovouchers_linked', 0)})."
            ]
            if has_cyto_line
            else []
        ),
        "",
        f"## 2. {n_st_word} Numbers at a Glance",
        "",
        "| Stage | Focal Plants | Shape K (Raw to Final) | Soil or Lineage CV Acc | Misid. Removed | Hybrids | Varieties or Ecotypes | Plastic Outliers | Cryptic Species | Soil Overlap (D) | Flowering Shift (Days) |",
        "| :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in summary_rows:
        cv_show = _nz(r.get("oof_lineage_morph_acc")) or _nz(r.get("spatial_block_cv_acc")) or "N/A"
        md_lines.append(
            f"| {r['stage']} | {r['n_focal_otus']} | {r.get('gmm_raw_k') or 'N/A'} to {r.get('gmm_optimal_k') or 'N/A'} | "
            f"{cv_show} | {r['n_misidentification_reassign']} | "
            f"{r['n_hybrid_introgressant']} | {r['n_subspecies_ecotype']} | "
            f"{r['n_ecophenotypic_plasticity']} | {r['n_cohesive_cryptic_otus']} | "
            f"{r.get('min_d_edaph') or 'N/A'} | {r.get('max_abs_delta_doy_days') or 'N/A'} |"
        )

    if s6_diag_rows:
        md_lines.extend([
            "",
            "## 3. Stage 6 Candidate Lineage Diagnosability and Key Characters",
            "",
            "| Candidate Lineage | Plants (Gold+Silver) | Morph Recall (75% Gate) | Median Sand % | Median Soil pH | Outcrop Index | Available Basionym | Recommended Rank |",
            "| :--- | ---: | ---: | ---: | ---: | ---: | :--- | :--- |",
        ])
        for dr in s6_diag_rows:
            md_lines.append(
                f"| `{dr.get('candidate_lineage', '')}` | {dr.get('n_plants', 0)} ({dr.get('n_gold_silver', 0)}) | "
                f"{dr.get('oof_morph_recall', 'N/A')} | {dr.get('median_soil_sand_pct', 'N/A')} | "
                f"{dr.get('median_soil_ph', 'N/A')} | {dr.get('median_lithic_outcrop_index', 'N/A')} | "
                f"{dr.get('available_basionym', 'Unassigned')} | `{dr.get('recommended_rank', '')}` |"
            )
        if s6_pair_rows:
            md_lines.extend(["", "**Pairwise Lineage Divergence:**", ""])
            for pr_w in s6_pair_rows[:6]:
                md_lines.append(
                    f"* `{pr_w.get('lineage_a', '')}` vs. `{pr_w.get('lineage_b', '')}`: "
                    f"held-out morph accuracy = {pr_w.get('pairwise_oof_acc', 'N/A')} "
                    f"(Mahalanobis D = {pr_w.get('mahalanobis_d', 'N/A')}, IQR separated = {pr_w.get('patten_unitt_iqr_sep', 'FALSE')}); "
                    f"soil overlap D_edaph = {pr_w.get('d_edaph', 'N/A')}; climate overlap D_clim = {pr_w.get('d_clim', 'N/A')}; "
                    f"flowering shift = {pr_w.get('delta_hopkins_doy_days', 'N/A')} days."
                )
        if s6_key_rows:
            md_lines.extend(["", "**Top Ranked Herbarium Characters Across Candidate Lineages:**", ""])
            for kr in s6_key_rows[:5]:
                md_lines.append(
                    f"* `{kr.get('trait_name', '')}` (effect size eta_sq = {kr.get('eta_sq', 'N/A')}, "
                    f"BH q = {kr.get('kw_q_value', 'N/A')}, non-overlapping IQR = {kr.get('has_non_overlapping_iqr', 'FALSE')}): "
                    f"{kr.get('lineage_iqr_summary', '')}"
                )

    sec_type_num = 4 if s6_diag_rows else 3
    sec_pull_num = sec_type_num + 1
    md_lines.extend(["", f"## {sec_type_num}. Type Specimens and Naming Anchors", ""])
    if type_controls or focal_std == "Packera dubia":
        for tr in type_controls[:10]:
            cat = tr.get("catalogNumber", "Unknown")
            inst = tr.get("institutionCode", "")
            cl = _nz(tr.get("gmm_cluster")) or "unclustered"
            t_sub = _nz(tr.get("taxonomic_treatment_subtype")) or "Confirmed_Type_Control"
            bas = (
                _nz(tr.get("nearest_type_basionym"))
                or _nz(tr.get("typified_basionym"))
                or _nz(tr.get("backbone_basionym_name"))
                or _nz(tr.get("basionym_name"))
                or _nz(tr.get("species_raw"))
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

    md_lines.extend(["", f"## {sec_pull_num}. Herbarium Sheets to Review First", ""])
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

    out_md_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_md = out_md_path.with_suffix(".md.tmp")
    with open(tmp_md, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(md_lines) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    tmp_md.replace(out_md_path)
    return out_md_path
