#!/usr/bin/env python3
"""
scripts/pipeline/step07_focal_cohort_utils.py
=============================================
Cohort naming, CSV/SQLite triage loaders, and quantitative cross-stage metric
extraction helpers for the four-stage focal species-complex orchestrator.
"""

from __future__ import annotations

import csv, math, re, sqlite3, statistics, sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon

NA_TOKENS: Set[str] = {"", "NA", "NAN", "NONE", "NULL"}


def slugify_taxon(taxon_name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9]+", "_", taxon_name.strip()).strip("_")
    return clean or "Packera_focal"


def default_prefix_for_taxon(taxon_name: str) -> str:
    parts = taxon_name.strip().split()
    if len(parts) >= 2 and parts[0].lower() == "packera":
        return re.sub(r"[^a-z0-9]+", "_", parts[1].lower()).strip("_")
    return slugify_taxon(taxon_name).lower()


def _stage_run_ids(prefix: str, discovery_run_id: Optional[str] = None) -> Tuple[str, str, str, str]:
    """Return ordered (r1, r2, r3, r4) run IDs for flat ('dubia') or nested ('parent/slug/') prefixes."""
    clean = prefix.strip()
    sep = "" if clean.endswith("/") else "_"
    r1 = f"{clean}{sep}01_eastern_discovery"
    r2 = f"{clean}{sep}02_complex_triage"
    r3 = f"{clean}{sep}03_intraspecific_clean"
    r4 = f"{clean}{sep}04_gold_thin_audit"
    if discovery_run_id is not None and str(discovery_run_id).strip():
        r1 = str(discovery_run_id).strip()
    return r1, r2, r3, r4


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def _safe_float(val: Any) -> Optional[float]:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s.upper() in NA_TOKENS:
        return None
    try:
        out = float(s)
        return out if math.isfinite(out) else None
    except ValueError:
        return None


def _first_float(row: Mapping[str, Any], cols: Sequence[str]) -> Optional[float]:
    for col in cols:
        if (val := _safe_float(row.get(col))) is not None:
            return val
    return None


def _is_truthy(val: Any) -> bool:
    return val is not None and str(val).strip().upper() in {"1", "TRUE", "T", "YES"}


def _fetch_triage_rows_for_run(
    conn: Optional[sqlite3.Connection], run_id: str, run_dir: Path
) -> List[Dict[str, Any]]:
    for cand_name in ("triage_queue_complex_full.csv", "triage_queue.csv"):
        if csv_rows := _read_csv_rows(run_dir / "data" / "tables" / cand_name):
            return csv_rows
    if conn is None:
        return []
    for query in (
        "SELECT * FROM v_integrative_triage_master WHERE run_id = ?;",
        "SELECT * FROM curation_triage WHERE run_id = ?;",
    ):
        try:
            cur = conn.execute(query, (run_id,))
            cols = [d[0] for d in cur.description]
            if rows := [dict(zip(cols, r)) for r in cur.fetchall()]:
                return rows
        except sqlite3.Error:
            continue
    return []


def _fetch_complex_rows_for_run(
    conn: Optional[sqlite3.Connection], run_id: str, run_dir: Path
) -> List[Dict[str, Any]]:
    if csv_rows := _read_csv_rows(run_dir / "outputs" / "reports" / "empirical_taxon_complexes.csv"):
        return csv_rows
    if conn is None:
        return []
    try:
        cur = conn.execute("SELECT * FROM taxon_complexes WHERE run_id = ?;", (run_id,))
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]
    except sqlite3.Error:
        return []


def _parse_gmm_metrics(gmm_rows: Sequence[Mapping[str, Any]], clusters: Set[str]) -> Dict[str, Any]:
    raw_k: Optional[int] = None
    opt_k: Optional[int] = None
    max_unmerged: Optional[float] = None
    two_delta_by_k: Dict[int, float] = {}
    delta_icl_by_k: Dict[int, float] = {}
    boot_p_by_k: Dict[int, float] = {}
    best_bic_val: Optional[float] = None
    best_bic_k: Optional[int] = None
    is_opt_k: Optional[int] = None

    for gr in gmm_rows:
        if raw_k is None and (rkv := _first_float(gr, ("optimal_K_raw", "optimal_k_raw", "optimal_K", "optimal_k", "G_opt", "best_K"))) is not None:
            raw_k = int(round(rkv))
        if opt_k is None and (kv := _first_float(gr, ("optimal_K_merged", "optimal_k_merged", "optimal_K", "optimal_k", "G_opt", "best_K"))) is not None:
            opt_k = int(round(kv))
        k_row = _safe_float(gr.get("K"))
        k_int = int(round(k_row)) if k_row is not None else None
        if is_opt_k is None and k_int is not None and _is_truthy(gr.get("is_optimal")):
            is_opt_k = k_int
        bic_v = _safe_float(gr.get("BIC"))
        if bic_v is not None and k_int is not None and (best_bic_val is None or bic_v > best_bic_val):
            best_bic_val, best_bic_k = bic_v, k_int
        if (bv := _first_float(gr, ("Two_Delta_BIC", "two_delta_BIC", "two_delta_bic", "delta_bic_2x", "bayes_factor_2lnB", "delta_BIC"))) is not None:
            if k_int is not None and k_int not in two_delta_by_k:
                two_delta_by_k[k_int] = bv
            if k_row is None or k_row >= 2.0:
                max_unmerged = bv if max_unmerged is None else max(max_unmerged, bv)
        if (iv := _first_float(gr, ("Delta_ICL_vs_Null", "delta_icl_vs_null"))) is not None and k_int is not None:
            delta_icl_by_k[k_int] = iv
        if (bpv := _first_float(gr, ("Boot_LRT_P_K1_vs_K2", "boot_lrt_p_k1_vs_k2"))) is not None:
            if k_int is not None:
                boot_p_by_k[k_int] = bpv
            elif -1 not in boot_p_by_k:
                boot_p_by_k[-1] = bpv

    if raw_k is None:
        raw_k = is_opt_k if is_opt_k is not None else best_bic_k
    if opt_k is None:
        opt_k = raw_k if raw_k is not None else (len(clusters) if clusters else None)
    if raw_k is None and opt_k is not None:
        raw_k = opt_k

    two_delta_at_opt: Optional[float] = None
    if opt_k is not None:
        if opt_k <= 1:
            two_delta_at_opt = 0.0
        elif opt_k in two_delta_by_k:
            two_delta_at_opt = two_delta_by_k[opt_k]
        elif max_unmerged is not None:
            two_delta_at_opt = max_unmerged

    delta_icl = next((delta_icl_by_k[ck] for ck in (raw_k, opt_k, 2) if ck is not None and ck in delta_icl_by_k), None)
    if delta_icl is None and delta_icl_by_k:
        delta_icl = delta_icl_by_k[sorted(delta_icl_by_k.keys())[-1]]

    boot_lrt_p = boot_p_by_k[2] if 2 in boot_p_by_k else (next(iter(boot_p_by_k.values())) if boot_p_by_k else None)
    return {
        "raw_k": raw_k, "opt_k": opt_k, "two_delta_bic_at_opt": two_delta_at_opt,
        "max_unmerged_two_delta_bic": max_unmerged, "delta_icl_vs_null": delta_icl,
        "boot_lrt_p_k1_vs_k2": boot_lrt_p,
    }


def _row_involves_focal(row_map: Mapping[str, Any], focal_std: str) -> bool:
    for tcol in ("Taxon1", "Taxon2", "taxon1", "taxon2"):
        tv = str(row_map.get(tcol) or "").strip()
        if tv and standardize_packera_taxon(tv) == focal_std:
            return True
    pair_raw = str(row_map.get("Pair") or row_map.get("pair") or "").strip()
    if pair_raw:
        for part in re.split(r"\s+vs\.?\s+|\s*[|/&]\s*|__", pair_raw, flags=re.IGNORECASE):
            p_clean = part.replace("_", " ").strip()
            if p_clean and standardize_packera_taxon(p_clean) == focal_std:
                return True
    return False


def _format_hybrid_parent_pairs(eval_rows: Sequence[Mapping[str, Any]]) -> str:
    counts: Counter[str] = Counter()
    for r in eval_rows:
        if str(r.get("taxonomic_treatment_subtype", "")).strip() != "Hybrid_Introgressant":
            continue
        p1, p2 = str(r.get("hybrid_parent_1") or "").strip(), str(r.get("hybrid_parent_2") or "").strip()
        if p1 and p2:
            pair_key = " x ".join(sorted([p1, p2]))
        else:
            rec = str(r.get("recommended_determination") or "")
            pair_key = " x ".join(sorted(p.strip() for p in rec.split(" x ", 1))) if " x " in rec else ""
        if pair_key:
            counts[pair_key] += 1
    return "; ".join(f"{pair} (n={cnt})" for pair, cnt in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def _read_genomics_csv(run_dir: Path, fname: str) -> List[Dict[str, str]]:
    for sub in ("outputs/genomics", "outputs/reports", "data/tables"):
        if rows_g := _read_csv_rows(run_dir / sub / fname):
            return rows_g
    return []


def _collect_genomics_and_cytology(
    eval_rows: Sequence[Mapping[str, Any]], run_dir: Path, focal_std: str
) -> Dict[str, Any]:
    gdi_vals: List[float] = []
    backbone_clades: List[str] = []
    n_tier0_abstained = 0
    reported_cytotypes_set: Set[str] = set()
    n_cytovouchers_linked = 0
    for r in eval_rows:
        if (gv := _first_float(r, ("coalescent_gdi", "cluster_propagated_gdi"))) is not None and math.isfinite(gv):
            gdi_vals.append(gv)
        bc = str(r.get("public_backbone_clade") or "").strip()
        if bc and bc.upper() not in NA_TOKENS:
            backbone_clades.append(bc)
        ab_raw = r.get("tier0_abstention_flag")
        if (isinstance(ab_raw, bool) and ab_raw) or (
            not isinstance(ab_raw, bool)
            and ab_raw is not None
            and (str(ab_raw).strip() in {"1", "1.0"} or _safe_float(ab_raw) == 1.0)
        ):
            n_tier0_abstained += 1
        cy_raw = str(r.get("taxon_reported_cytotypes") or "").strip()
        if cy_raw and cy_raw.upper() not in NA_TOKENS:
            for tok in re.split(r"[;,]+", cy_raw):
                if (t_clean := tok.strip()) and t_clean.upper() not in NA_TOKENS:
                    reported_cytotypes_set.add(t_clean)
        if str(r.get("cytology_linkage_tier") or "").strip() in {"C1_Direct_Cytovoucher", "C2_Exsiccatae_Cytoduplicate"}:
            n_cytovouchers_linked += 1

    if not gdi_vals:
        for gr in _read_genomics_csv(run_dir, "coalescent_gdi_delimitation.csv"):
            if str(gr.get("gdi_source") or "").strip() not in ("", "nextrad_msc", "hybseq_supermatrix_msc"):
                continue
            if focal_std and standardize_packera_taxon(str(gr.get("taxon") or "")) != focal_std:
                continue
            if (gv := _safe_float(gr.get("coalescent_gdi"))) is not None and math.isfinite(gv):
                gdi_vals.append(gv)

    if not backbone_clades:
        for pr in _read_genomics_csv(run_dir, "public_backbone_concordance.csv"):
            if focal_std and standardize_packera_taxon(str(pr.get("taxon") or pr.get("normalized_taxon") or "")) != focal_std:
                continue
            bc = str(pr.get("public_backbone_clade") or "").strip()
            if bc and bc.upper() not in NA_TOKENS:
                backbone_clades.append(bc)

    gm_rows = _read_genomics_csv(run_dir, "genomic_model_metrics.csv")
    return {
        "median_coalescent_gdi": float(statistics.median(gdi_vals)) if gdi_vals else None,
        "public_backbone_clade": (
            sorted(Counter(backbone_clades).items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
            if backbone_clades
            else ""
        ),
        "active_genomic_tier": str(gm_rows[0].get("active_genomic_tier") or "").strip() if gm_rows else "",
        "n_tier0_abstained": n_tier0_abstained,
        "reported_cytotypes": "; ".join(sorted(reported_cytotypes_set)),
        "n_cytovouchers_linked": n_cytovouchers_linked,
    }


def _collect_complex_neighbors(complex_rows: Sequence[Mapping[str, Any]], run_dir: Path, focal_std: str) -> List[str]:
    neighbors: List[str] = []
    focal_cid: Optional[str] = None
    for cr in complex_rows:
        if standardize_packera_taxon(str(cr.get("species_standardized") or "")) == focal_std:
            focal_cid = str(cr.get("consensus_complex_id") or "").strip()
            nn = standardize_packera_taxon(str(cr.get("nearest_neighbor_taxon") or ""))
            if nn and nn != focal_std:
                neighbors.append(nn)
    if focal_cid:
        for cr in complex_rows:
            sp = standardize_packera_taxon(str(cr.get("species_standardized") or ""))
            if str(cr.get("consensus_complex_id") or "").strip() == focal_cid and sp and sp != focal_std and sp not in neighbors:
                neighbors.append(sp)
    neigh_txt = run_dir / "outputs" / "reports" / "focal_auto_neighbors.txt"
    if neigh_txt.is_file():
        for ln in neigh_txt.read_text(encoding="utf-8").splitlines():
            sp = standardize_packera_taxon(ln.strip())
            if sp and sp != focal_std and sp not in neighbors:
                neighbors.append(sp)
    return neighbors


def _extract_stage_metrics(
    *, stage_name: str, run_id: str, focal_std: str, runs_root: Path, conn: Optional[sqlite3.Connection]
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    run_dir = runs_root / run_id
    triage_rows = _fetch_triage_rows_for_run(conn, run_id, run_dir)
    complex_rows = _fetch_complex_rows_for_run(conn, run_id, run_dir)
    gmm_rows: List[Dict[str, str]] = []
    for gmm_cand in (
        run_dir / "outputs" / "reports" / "gmm_bic_model_selection.csv",
        run_dir / "data" / "tables" / "gmm_bic_model_selection.csv",
        run_dir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv",
    ):
        if gmm_rows := _read_csv_rows(gmm_cand):
            break
    srf_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "multimodal_spatial_rf_summary.csv")
    srf_flag_rows = _read_csv_rows(run_dir / "data" / "tables" / "multimodal_conflict_flags.csv")
    warren_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "warren_niche_identity_summary.csv")
    batch_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "institutional_batch_effect_audit.csv")

    focal_rows = [
        r for r in triage_rows
        if standardize_packera_taxon(str(r.get("species_standardized") or r.get("scientificName") or "")) == focal_std
    ]
    eval_rows = focal_rows if focal_rows else triage_rows
    non_focal_rows = (
        [r for r in triage_rows if standardize_packera_taxon(str(r.get("species_standardized") or r.get("scientificName") or "")) != focal_std]
        if focal_rows
        else []
    )

    def _count_sub(subtype: str, *, pool: Optional[List[Dict[str, Any]]] = None) -> int:
        rows_chk = pool if pool is not None else eval_rows
        if subtype in {"Confirmed_Type_Control", "Type_Synonymy_Candidate"}:
            return sum(
                1 for r in rows_chk
                if str(r.get("taxonomic_treatment_subtype", "")).strip() == subtype
                or str(r.get("taxonomic_status_call", "")).strip() == subtype
            )
        return sum(1 for r in rows_chk if str(r.get("taxonomic_treatment_subtype", "")).strip() == subtype)

    clusters = {s for r in eval_rows if (s := str(r.get("gmm_cluster", "")).strip()) not in {"", "NA", "None"}}
    taxa_present = {s for r in triage_rows if (s := str(r.get("species_standardized", "")).strip()) not in {"", "NA", "None"}}
    gmm = _parse_gmm_metrics(gmm_rows, clusters)

    srf_scope = [sr for sr in srf_rows if _row_involves_focal(sr, focal_std)] or srf_rows
    d_edaph_vals = [v for sr in srf_scope if (v := _safe_float(sr.get("D_edaph"))) is not None]
    if not d_edaph_vals:
        w_pool = warren_rows or srf_rows
        w_scope = [wr for wr in w_pool if _row_involves_focal(wr, focal_std)] or w_pool
        d_edaph_vals = [v for wr in w_scope if (v := _safe_float(wr.get("Schoener_D"))) is not None]
    min_d_edaph = min(d_edaph_vals) if d_edaph_vals else None

    cv_acc_vals: List[float] = []
    for src_rows in (srf_rows, srf_flag_rows, eval_rows):
        cv_acc_vals = [
            v for r_cv in src_rows
            if (v := _first_float(r_cv, ("spatial_block_cv_acc", "Spatial_Block_CV_Acc", "spatial_block_cv_accuracy"))) is not None
        ]
        if cv_acc_vals:
            break
    spatial_block_cv_acc = float(statistics.mean(cv_acc_vals)) if cv_acc_vals else None

    doy_vals = [abs(v) for r in eval_rows if (v := _safe_float(r.get("delta_doy_resid"))) is not None]
    max_abs_doy = max(doy_vals) if doy_vals else None
    eta_vals = [v for br in batch_rows if (v := _first_float(br, ("partial_eta_sq", "R_squared"))) is not None]
    max_eta = max(eta_vals) if eta_vals else None
    cohesive_cryptic = sum(
        1 for r in eval_rows
        if str(r.get("taxonomic_treatment_subtype", "")).strip() == "Putative_Cryptic_Species"
        and _is_truthy(r.get("cryptic_cluster_cohesion_pass"))
    )
    neighbors = _collect_complex_neighbors(complex_rows, run_dir, focal_std)
    geno_cyto = _collect_genomics_and_cytology(eval_rows, run_dir, focal_std)

    def _fmt(val: Optional[float], spec: str) -> str:
        return format(val, spec) if val is not None else ""

    metrics: Dict[str, Any] = {
        "stage": stage_name,
        "run_id": run_id,
        "n_total_otus": len(triage_rows),
        "n_focal_otus": len(focal_rows),
        "n_taxa_present": len(taxa_present),
        "gmm_raw_k": gmm["raw_k"] if gmm["raw_k"] is not None else "",
        "gmm_optimal_k": gmm["opt_k"] if gmm["opt_k"] is not None else "",
        "two_delta_bic_at_optimal_k": _fmt(gmm["two_delta_bic_at_opt"], ".2f"),
        "max_unmerged_two_delta_bic": _fmt(gmm["max_unmerged_two_delta_bic"], ".2f"),
        "max_two_delta_bic": _fmt(gmm["max_unmerged_two_delta_bic"], ".2f"),
        "delta_icl_vs_null": _fmt(gmm["delta_icl_vs_null"], ".2f"),
        "boot_lrt_p_k1_vs_k2": _fmt(gmm["boot_lrt_p_k1_vs_k2"], ".4f"),
        "spatial_block_cv_acc": _fmt(spatial_block_cv_acc, ".4f"),
        "n_confirmed_species": _count_sub("Confirmed_Species"),
        "n_confirmed_type_control": _count_sub("Confirmed_Type_Control"),
        "n_type_synonymy_candidate": _count_sub("Type_Synonymy_Candidate"),
        "n_neighbor_type_controls": _count_sub("Confirmed_Type_Control", pool=non_focal_rows),
        "n_misidentification_reassign": _count_sub("Misidentification_Reassign"),
        "n_hybrid_introgressant": _count_sub("Hybrid_Introgressant"),
        "hybrid_parent_pairs": _format_hybrid_parent_pairs(eval_rows),
        "n_subspecies_ecotype": _count_sub("Subspecies_Ecotype"),
        "n_ecophenotypic_plasticity": _count_sub("Ecophenotypic_Plasticity"),
        "n_putative_cryptic_species": _count_sub("Putative_Cryptic_Species"),
        "n_cohesive_cryptic_otus": cohesive_cryptic,
        "min_d_edaph": _fmt(min_d_edaph, ".3f"),
        "max_abs_delta_doy_days": _fmt(max_abs_doy, ".1f"),
        "max_batch_partial_eta_sq": _fmt(max_eta, ".4f"),
        "median_coalescent_gdi": geno_cyto["median_coalescent_gdi"],
        "public_backbone_clade": geno_cyto["public_backbone_clade"],
        "active_genomic_tier": geno_cyto["active_genomic_tier"],
        "n_tier0_abstained": geno_cyto["n_tier0_abstained"],
        "reported_cytotypes": geno_cyto["reported_cytotypes"],
        "n_cytovouchers_linked": geno_cyto["n_cytovouchers_linked"],
        "complex_neighbors": "; ".join(neighbors),
    }
    return metrics, eval_rows
