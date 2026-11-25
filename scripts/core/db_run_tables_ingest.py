"""
Run-scoped Step 03-07 CSV table ingestion into Tier B SQLite tables.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sqlite3
from typing import Any, Dict, List, Union

from scripts.core.db_column_constants import _to_bool_int, _to_float, _to_int, _to_str
from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated, upsert_rows_generic
from scripts.core.db_genomics_cytology_persist import persist_genetic_lineage_assignments
from scripts.core.db_morph_and_rf_persist import (
    persist_leaf_efa_harmonics,
    persist_otu_morphometrics,
    persist_spatial_rf_predictions,
)
from scripts.core.db_step02_loader import ensure_pipeline_run

def _read_csv_dicts(csv_path: Path) -> List[Dict[str, str]]:
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        return []
    with open(csv_path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f))

def ingest_existing_run_tables_into_db(
    conn: sqlite3.Connection,
    tables_dir: Union[str, Path],
    run_id: str = "canonical",
    demo: bool = False,
) -> Dict[str, int]:
    """
    Ingest Step 03-07 run-scoped CSV outputs (`leaf_efa_harmonics_individual.csv`,
    `leaf_efa_harmonics.csv`, `morphometric_flags.csv`,
    `multimodal_conflict_flags.csv`, `triage_queue.csv`)
    into Tier B SQLite tables with full dynamic `A1..D24`, `PC1..PC15`,
    and 19-variable `bio1..bio19` support.
    """
    if demo:
        return {}
    t_dir = Path(tables_dir).resolve()
    eff_run = str(run_id).strip() or "canonical"
    ensure_pipeline_run(conn, run_id=eff_run, run_dir=t_dir.parent)
    run_counts: Dict[str, int] = {}

    ind_efa = _read_csv_dicts(t_dir / "leaf_efa_harmonics_individual.csv")
    if ind_efa:
        run_counts["leaf_efa_harmonics"] = persist_leaf_efa_harmonics(
            conn, eff_run, ind_efa, demo=False
        )

    step03_otu = _read_csv_dicts(t_dir / "leaf_efa_harmonics.csv")
    step04_otu = _read_csv_dicts(t_dir / "morphometric_flags.csv")
    n_otu = 0
    if step03_otu:
        n_otu = persist_otu_morphometrics(
            conn, eff_run, step03_otu, demo=False, preserve_existing_on_null=True
        )
    if step04_otu:
        n_otu = persist_otu_morphometrics(
            conn, eff_run, step04_otu, demo=False, preserve_existing_on_null=True
        )
    if step03_otu or step04_otu:
        run_counts["otu_morphometrics"] = n_otu

    srf_rows = _read_csv_dicts(t_dir / "multimodal_conflict_flags.csv")
    if srf_rows:
        run_counts["spatial_rf_predictions"] = persist_spatial_rf_predictions(
            conn, eff_run, srf_rows, demo=False, also_upsert_env_cache=False
        )

    tri_path = t_dir / "triage_queue.csv"
    tri_rows = _read_csv_dicts(tri_path)
    if tri_rows:
        _ensure_dynamic_table_migrated(conn, "curation_triage")
        norm_tri: List[Dict[str, Any]] = []
        for r in tri_rows:
            cat = _to_str(r.get("catalogNumber"))
            if not cat:
                continue
            d_tri: Dict[str, Any] = {
                "run_id": eff_run,
                "catalogNumber": cat,
                "plant_individual_id": _to_int(r.get("plant_individual_id"), default=0) or 0,
            }
            for sc in (
                "taxonomic_status_call",
                "recommended_determination",
                "taxonomic_treatment_broad",
                "taxonomic_treatment_subtype",
                "synthesis_triage_priority",
                "synthesis_triage_action",
                "synthesis_rationale",
                "resurrection_basionym",
                "hybrid_parent_1",
                "hybrid_parent_2",
                "genomic_propagation_tier",
                "tier0_calibrated_taxon",
                "cytology_linkage_tier",
                "taxon_reported_cytotypes",
                "hydro_edaphic_niche_regime",
                "basionym_synonymy_type",
                "portal_lag_trigger_gate",
            ):
                if sc in r:
                    d_tri[sc] = _to_str(r.get(sc))
            for ic in (
                "n_available_modalities",
                "tier0_abstention_flag",
                "voucher_sporophytic_2n",
                "voucher_ploidy_level",
            ):
                if ic in r:
                    iv = _to_int(r.get(ic))
                    if iv is None and ic == "tier0_abstention_flag":
                        iv = _to_bool_int(r.get(ic))
                    d_tri[ic] = iv
            for bc in (
                "cryptic_cluster_cohesion_pass",
                "cryptic_non_clinal_pass",
                "hybrid_zone_cohesion_pass",
                "is_completed_voucher",
                "is_multi_cytotype_taxon",
                "cytotype_reproductive_barrier_pass",
                "portal_metadata_lag_flag",
                "is_marginal_range_extension",
            ):
                if bc in r:
                    d_tri[bc] = _to_bool_int(r.get(bc))
            for fc in (
                "w_morph_kappa",
                "w_edaphic_kappa",
                "w_pheno_kappa",
                "w_geno_kappa",
                "w_morph_simplex",
                "w_edaphic_simplex",
                "w_pheno_simplex",
                "w_geno_simplex",
                "geno_deme_weight_scale",
                "probabilistic_consensus_conflict",
                "probabilistic_consensus_score",
                "hopkins_adjusted_doy",
                "delta_doy_resid",
                "tier0_calibrated_prob",
                "deme_mahal_d2",
                "cluster_propagated_gdi",
                "wetland_fidelity_concordance",
            ):
                if fc in r:
                    d_tri[fc] = _to_float(r.get(fc))
            norm_tri.append(d_tri)
        run_counts["curation_triage"] = upsert_rows_generic(
            conn, "curation_triage", norm_tri, ("run_id", "catalogNumber", "plant_individual_id")
        )

    gla_rows = _read_csv_dicts(t_dir / "genomic_assignments.csv")
    if gla_rows:
        run_counts["genetic_lineage_assignments"] = persist_genetic_lineage_assignments(
            conn, eff_run, gla_rows, demo=False
        )

    return run_counts

__all__ = (
    "_read_csv_dicts", "ingest_existing_run_tables_into_db",
)
