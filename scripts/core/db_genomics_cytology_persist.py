"""
Tier A genomics, cytology, taxonomic backbone, NatureServe, and Tier B lineage persistence.
"""

from __future__ import annotations

import sqlite3
from typing import Any, Dict, List, Mapping, Sequence

from scripts.core.db_column_constants import _to_bool_int, _to_float, _to_int, _to_str
from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated, upsert_rows_generic
from scripts.core.db_step02_loader import ensure_pipeline_run

def upsert_genomic_samples(
    conn: sqlite3.Connection,
    rows: Sequence[Mapping[str, Any]],
    *,
    preserve_existing_on_null: bool = True,
) -> int:
    """Upsert Tier A genomic sample metadata rows into `genomic_samples`."""
    if not rows:
        return 0
    str_cols = (
        "genetic_sample_id", "catalogNumber", "project_id", "population_id",
        "replicate_pair_id", "taxon_raw", "species_standardized", "marker_platform",
        "stateProvince", "county", "eventDate", "linkage_tier",
        "inferred_cytotype", "plastid_clade", "its_ribotype",
        "public_plastome_match", "public_its_match", "waxy_lineage_match",
        "exsiccatae_key", "genbank_accession", "coord_quarantine_reason",
        "typified_basionym", "concept_synonymy_status",
    )
    int_cols = ("ploidy_level", "hq_reads", "mapped_reads", "n_called_loci", "hdplot_retained_loci")
    bool_cols = (
        "is_putative_novel_or_hybrid_query",
        "geno_quality_pass",
        "cytonuclear_discordant",
        "is_centroid_or_inst_coord",
        "is_open_or_hybrid_concept",
    )
    zero_default_bool_cols = {
        "is_putative_novel_or_hybrid_query",
        "is_centroid_or_inst_coord",
        "is_open_or_hybrid_concept",
    }
    float_cols = (
        "pct_repetitive", "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
        "missing_genotype_rate", "observed_het", "nquire_delta_logl",
    )

    norm_rows: List[Dict[str, Any]] = []
    for raw in rows:
        gid = _to_str(raw.get("genetic_sample_id"))
        if not gid:
            continue
        item: Dict[str, Any] = {"genetic_sample_id": gid}
        for col in str_cols:
            if col != "genetic_sample_id":
                item[col] = _to_str(raw.get(col))
        for col in int_cols:
            item[col] = _to_int(raw.get(col))
        for col in bool_cols:
            bval = _to_bool_int(raw.get(col))
            item[col] = 0 if (col in zero_default_bool_cols and bval is None) else bval
        for col in float_cols:
            item[col] = _to_float(raw.get(col))
        norm_rows.append(item)

    return upsert_rows_generic(
        conn,
        "genomic_samples",
        norm_rows,
        ("genetic_sample_id",),
        preserve_existing_on_null=preserve_existing_on_null,
    )

def upsert_cytological_records(
    conn: sqlite3.Connection,
    rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    *,
    preserve_existing_on_null: bool = False,
) -> int:
    """Upsert Tier A cytological voucher and literature records into `cytological_records`."""
    if demo or not rows:
        return 0
    str_cols = (
        "source_record_id",
        "basionym_name",
        "catalogNumber",
        "exsiccatae_key",
        "cytology_linkage_tier",
        "raw_count_str",
        "cytotype_label",
        "breeding_system",
        "genbank_accession",
        "collector_surname",
        "record_number_clean",
        "institution_code",
        "voucher_citation",
        "locality_text",
        "reference_citation",
    )
    int_cols = ("gametophytic_n", "sporophytic_2n", "base_x", "ploidy_level")
    float_cols = ("c_value_1c_pg", "c_value_2c_pg")

    norm_rows: List[Dict[str, Any]] = []
    for raw in rows:
        cid = _to_str(raw.get("cytology_id"))
        sp_std = _to_str(raw.get("species_standardized"))
        if not cid or not sp_std:
            continue
        item: Dict[str, Any] = {
            "cytology_id": cid,
            "source_db": _to_str(raw.get("source_db")) or "local",
            "original_taxon_name": _to_str(raw.get("original_taxon_name")) or sp_std,
            "species_standardized": sp_std,
            "b_chromosomes": _to_int(raw.get("b_chromosomes"), default=0) or 0,
            "is_aneuploid_or_odd_ploid": _to_bool_int(raw.get("is_aneuploid_or_odd_ploid")) or 0,
            "is_approximate": _to_bool_int(raw.get("is_approximate")) or 0,
        }
        for col in str_cols:
            item[col] = _to_str(raw.get(col))
        for col in int_cols:
            item[col] = _to_int(raw.get(col))
        for col in float_cols:
            item[col] = _to_float(raw.get(col))
        norm_rows.append(item)

    if not norm_rows:
        return 0
    return upsert_rows_generic(
        conn,
        "cytological_records",
        norm_rows,
        ("cytology_id",),
        preserve_existing_on_null=preserve_existing_on_null,
    )

def upsert_taxon_cytology_summary(
    conn: sqlite3.Connection,
    rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    *,
    preserve_existing_on_null: bool = False,
) -> int:
    """Upsert Tier A taxon-level cytological summaries into `taxon_cytology_summary`."""
    if demo or not rows:
        return 0
    norm_rows: List[Dict[str, Any]] = []
    for raw in rows:
        sp_std = _to_str(raw.get("species_standardized"))
        if not sp_std:
            continue
        item: Dict[str, Any] = {
            "species_standardized": sp_std,
            "n_cytology_records": _to_int(raw.get("n_cytology_records"), default=0) or 0,
            "n_direct_cytovouchers": _to_int(raw.get("n_direct_cytovouchers"), default=0) or 0,
            "modal_sporophytic_2n": _to_int(raw.get("modal_sporophytic_2n")),
            "modal_base_x": _to_int(raw.get("modal_base_x")),
            "modal_ploidy_level": _to_int(raw.get("modal_ploidy_level")),
            "reported_cytotypes": _to_str(raw.get("reported_cytotypes")),
            "is_multi_cytotype_taxon": _to_bool_int(raw.get("is_multi_cytotype_taxon")) or 0,
            "has_dysploid_series": _to_bool_int(raw.get("has_dysploid_series")) or 0,
            "has_odd_polyploid_or_apomixis": _to_bool_int(raw.get("has_odd_polyploid_or_apomixis")) or 0,
            "mean_2c_pg": _to_float(raw.get("mean_2c_pg")),
            "breeding_system_summary": _to_str(raw.get("breeding_system_summary")),
        }
        norm_rows.append(item)

    if not norm_rows:
        return 0
    return upsert_rows_generic(
        conn,
        "taxon_cytology_summary",
        norm_rows,
        ("species_standardized",),
        preserve_existing_on_null=preserve_existing_on_null,
    )

def upsert_taxon_backbone_profiles(
    conn: sqlite3.Connection,
    rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    *,
    preserve_existing_on_null: bool = True,
) -> int:
    """Upsert Tier A taxonomic, protologue, NatureServe rank, and NWPL profiles into `taxon_backbone_profiles`."""
    if demo or not rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "taxon_backbone_profiles")
    str_cols = (
        "accepted_scientific_name",
        "accepted_authorship",
        "basionym_name",
        "basionym_authorship",
        "homotypic_synonyms_json",
        "heterotypic_synonyms_json",
        "hybrid_formula",
        "publication_citation",
        "protologue_url",
        "locus_classicus_text",
        "tdwg_native_codes_json",
        "tdwg_introduced_codes_json",
        "natureserve_element_id",
        "natureserve_global_rank",
        "nwpl_indicator_code",
        "backbone_source",
        "updated_at",
    )
    has_any_tg = any(_to_str(raw.get("target_genus")) for raw in rows)
    norm_rows: List[Dict[str, Any]] = []
    for raw in rows:
        sp_std = _to_str(raw.get("species_standardized"))
        if not sp_std:
            continue
        item: Dict[str, Any] = {"species_standardized": sp_std}
        if has_any_tg:
            item["target_genus"] = _to_str(raw.get("target_genus")) or "Packera"
        for col in str_cols:
            if col in raw or not preserve_existing_on_null:
                item[col] = _to_str(raw.get(col))
        if "is_hybrid_taxon" in raw:
            bval = _to_bool_int(raw.get("is_hybrid_taxon"))
            item["is_hybrid_taxon"] = (0 if bval is None and not preserve_existing_on_null else bval)
        elif not preserve_existing_on_null:
            item["is_hybrid_taxon"] = 0
        if "publication_year" in raw or not preserve_existing_on_null:
            item["publication_year"] = _to_int(raw.get("publication_year"))
        if "nwpl_wetland_score" in raw or not preserve_existing_on_null:
            item["nwpl_wetland_score"] = _to_float(raw.get("nwpl_wetland_score"))
        norm_rows.append(item)

    if not norm_rows:
        return 0
    return upsert_rows_generic(
        conn,
        "taxon_backbone_profiles",
        norm_rows,
        ("species_standardized",),
        preserve_existing_on_null=preserve_existing_on_null,
    )

def upsert_natureserve_subnational_ranges(
    conn: sqlite3.Connection,
    rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    *,
    preserve_existing_on_null: bool = False,
) -> int:
    """Upsert Tier A NatureServe subnational US state and Canadian province rows into `natureserve_subnational_ranges`."""
    if demo or not rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "natureserve_subnational_ranges")
    norm_rows: List[Dict[str, Any]] = []
    for raw in rows:
        sp_std = _to_str(raw.get("species_standardized"))
        c_code = (_to_str(raw.get("country_code")) or "").upper()
        sub_code = (_to_str(raw.get("subnation_code")) or "").upper()
        if not sp_std or not c_code or not sub_code:
            continue
        is_nat = _to_bool_int(raw.get("is_native"))
        is_exo = _to_bool_int(raw.get("is_exotic"))
        is_hist = _to_bool_int(raw.get("is_historical_or_extirpated"))
        item: Dict[str, Any] = {
            "species_standardized": sp_std,
            "country_code": c_code,
            "subnation_code": sub_code,
            "s_rank": _to_str(raw.get("s_rank")),
            "is_native": 1 if is_nat is None else is_nat,
            "is_exotic": 0 if is_exo is None else is_exo,
            "is_historical_or_extirpated": 0 if is_hist is None else is_hist,
        }
        norm_rows.append(item)

    if not norm_rows:
        return 0
    return upsert_rows_generic(
        conn,
        "natureserve_subnational_ranges",
        norm_rows,
        ("species_standardized", "country_code", "subnation_code"),
        preserve_existing_on_null=preserve_existing_on_null,
    )

def persist_genetic_lineage_assignments(
    conn: sqlite3.Connection,
    run_id: str,
    rows: Sequence[Mapping[str, Any]],
    *,
    demo: bool = False,
    prune_stale: bool = True,
) -> int:
    """Persist Tier B run-scoped genomic lineage assignments and optionally prune stale rows for `run_id`."""
    if demo or not rows:
        return 0
    eff_run = str(run_id).strip() or "canonical"
    ensure_pipeline_run(conn, run_id=eff_run)

    str_cols = (
        "genetic_sample_id", "catalogNumber", "population_id", "geno_predicted_taxon",
        "genetic_cluster", "q_admixture_json", "geno_parent_1", "geno_parent_2",
        "nearest_tree_sister_taxon", "hybrid_generation_class", "coancestry_cluster",
        "f_branch_donor_taxon", "public_backbone_clade",
    )
    float_cols = (
        "geno_posterior_prob", "geno_margin", "q_max",
        "patristic_dist_to_dubia", "replicate_pair_concordance",
        "hybrid_index_h", "interspecific_ancestry_q12", "f_branch_max",
        "snaq_gamma_minor", "coalescent_gdi", "bpp_introgression_bf",
        "mmrr_beta_ibd", "mmrr_beta_ibe", "mmrr_p_ibe",
        "public_phytop_ils_index", "public_phytop_ih_index",
    )
    norm_rows: List[Dict[str, Any]] = []
    kept_ids: List[str] = []
    for raw in rows:
        gid = _to_str(raw.get("genetic_sample_id"))
        if not gid:
            continue
        kept_ids.append(gid)
        item: Dict[str, Any] = {"run_id": eff_run, "genetic_sample_id": gid}
        for col in str_cols:
            if col != "genetic_sample_id":
                item[col] = _to_str(raw.get(col))
        for col in float_cols:
            item[col] = _to_float(raw.get(col))
        item["is_geno_admixed"] = _to_bool_int(raw.get("is_geno_admixed"))
        norm_rows.append(item)

    if not norm_rows:
        return 0
    written = upsert_rows_generic(
        conn,
        "genetic_lineage_assignments",
        norm_rows,
        ("run_id", "genetic_sample_id"),
        preserve_existing_on_null=False,
    )
    if prune_stale and kept_ids:
        placeholders = ", ".join("?" for _ in kept_ids)
        conn.execute(
            f"DELETE FROM genetic_lineage_assignments WHERE run_id = ? AND genetic_sample_id NOT IN ({placeholders});",
            [eff_run, *kept_ids],
        )
        conn.commit()
    return written

__all__ = (
    "upsert_genomic_samples", "upsert_cytological_records", "upsert_taxon_cytology_summary",
    "upsert_taxon_backbone_profiles", "upsert_natureserve_subnational_ranges",
    "persist_genetic_lineage_assignments",
)
