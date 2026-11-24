"""step01_corpus_authority_discovery.py - Layer 3 empirical specialist and regional herbarium discovery."""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.step01_determiner_surnames import (
    PACKERA_SPECIALISTS,
    _AMBIGUOUS_COMMON_SURNAMES,
    _INVALID_INSTITUTION_CODES,
    _has_non_person_determination_token,
    _is_explicit_collector_reference,
    _is_verified_type_status,
    normalize_determiner_surnames,
)
from scripts.core.step01_exsiccatae_fingerprint import (
    _resolve_row_standardized_packera,
    normalize_collector,
    normalize_event_date,
    normalize_record_number,
)


def discover_corpus_authorities(
    df: pd.DataFrame,
    min_gold_determinations: Optional[int] = None,
    min_gold_herbaria: int = 3,
    min_herbarium_vouchers: int = 10,
) -> Tuple[Set[str], Set[str]]:
    """
    Layer 3 Authority Engine: Discovers empirical taxonomic specialists and regional
    research herbaria directly from the assembled Darwin Core voucher corpus.

    Rules:
      - Dynamic Specialist Discovery (`dynamic_specialists`):
        1. Determiner (`identifiedBy` only, never falling back to `recordedBy` or
           matching `_NON_PERSON_DETERMINATION_TOKENS`) has identified >= `min_gold_determinations`
           specimens across >= `min_gold_herbaria` distinct `institutionCode`s (and in full
           corpora with >= 100 rows, either >= 4 distinct Packera species or >= 1 type specimen), OR
        2. Determiner has explicitly annotated >= 1 verified type specimen (giving +2
           weight to type determinations, excluding collector-reference fallbacks) AND
           reaches weighted count >= 3 across >= 2 distinct `institutionCode`s, AND
        3. Bare surnames in `_AMBIGUOUS_COMMON_SURNAMES` are excluded unless already in
           `PACKERA_SPECIALISTS` or backed by >= 1 verified type determination, AND
        4. Determiner has cross-duplicate disagreement rate <= 15% across shared
           `exsiccatae_key` duplicate groups when >= 2 duplicate groups exist.
      - Dynamic Regional Research Herbarium Discovery (`dynamic_herbaria`):
        1. `institutionCode` holds >= 1 verified type specimen for target Packera taxa
           AND >= 5 total sheets in the corpus, OR
        2. `institutionCode` ranks in the top quartile of institutional holdings
           (>= `min_herbarium_vouchers` curated sheets) with >= 70% determiner-attributed
           records for any harvested taxon.

    Returns:
        Tuple[Set[str], Set[str]]: (dynamic_specialists, dynamic_herbaria)
    """
    dynamic_specialists: Set[str] = set()
    dynamic_herbaria: Set[str] = set()
    if df is None or df.empty:
        return dynamic_specialists, dynamic_herbaria

    is_large_corpus = len(df) >= 100
    eff_min_det = (15 if is_large_corpus else 5) if min_gold_determinations is None else int(min_gold_determinations)
    eff_min_herbs = max(1, int(min_gold_herbaria))
    eff_min_herb_vouchers = max(1, int(min_herbarium_vouchers))

    id_col = (
        "_raw_identified_by"
        if "_raw_identified_by" in df.columns
        else ("identifiedBy" if "identifiedBy" in df.columns else None)
    )
    sp_col = (
        "species_standardized"
        if "species_standardized" in df.columns
        else ("scientificName" if "scientificName" in df.columns else ("species_raw" if "species_raw" in df.columns else None))
    )
    type_col = (
        "_raw_type_status"
        if "_raw_type_status" in df.columns
        else ("type_status" if "type_status" in df.columns else ("typeStatus" if "typeStatus" in df.columns else None))
    )
    inst_col = "institutionCode" if "institutionCode" in df.columns else None
    has_exsiccatae_col = "exsiccatae_key" in df.columns

    # Track (determiner_surname, standardized_taxon) -> [count, set_of_institutions]
    det_taxon_counts: Dict[Tuple[str, str], int] = {}
    det_taxon_insts: Dict[Tuple[str, str], Set[str]] = {}
    det_total_counts: Dict[str, int] = {}
    det_total_insts: Dict[str, Set[str]] = {}
    det_species: Dict[str, Set[str]] = {}
    det_type_hits: Dict[str, int] = {}
    det_type_annotated: Set[str] = set()

    # Track institutional holdings, type holdings, and determiner-attribution rates
    inst_total_sheets: Dict[str, int] = {}
    inst_type_sheets: Dict[str, int] = {}
    inst_taxon_totals: Dict[Tuple[str, str], int] = {}
    inst_taxon_attributed: Dict[Tuple[str, str], int] = {}

    from scripts.core.harvester import _extract_canonical_binomial, standardize_packera_taxon

    # Track duplicate groups for cross-duplicate disagreement auditing (<= 15% threshold)
    group_entries: Dict[str, List[Tuple[str, List[str], bool]]] = {}

    for _, row in df.iterrows():
        raw_sp = row.get(sp_col) if sp_col else None
        std_sp = standardize_packera_taxon(raw_sp)
        if not std_sp and raw_sp is not None and not pd.isna(raw_sp):
            s_tmp = str(raw_sp).strip()
            if s_tmp.startswith("Packera "):
                std_sp = _extract_canonical_binomial(s_tmp, {"Packera"})
        if not std_sp:
            std_sp = _resolve_row_standardized_packera(row)

        raw_inst = row.get(inst_col) if inst_col else ""
        inst = str(raw_inst if not pd.isna(raw_inst) else "").strip().upper()
        valid_inst = inst not in _INVALID_INSTITUTION_CODES

        raw_id = row.get(id_col) if id_col else ""
        is_explicit_col_ref = _is_explicit_collector_reference(raw_id)
        is_non_person_id = _has_non_person_determination_token(raw_id)
        has_explicit_det_str = (
            raw_id is not None
            and not (pd is not None and pd.isna(raw_id))
            and bool(str(raw_id).strip())
            and not is_explicit_col_ref
            and not is_non_person_id
        )
        # Never fall back to recordedBy when counting determinations in discover_corpus_authorities
        surnames = normalize_determiner_surnames(raw_id) if has_explicit_det_str else []
        has_attributed_det = len(surnames) > 0
        raw_rec_by = row.get("recordedBy", "") if "recordedBy" in row else ""
        rec_col_norm = normalize_collector("" if (raw_rec_by is None or pd.isna(raw_rec_by)) else str(raw_rec_by))
        cross_collector_surnames = [
            s_n for s_n in surnames
            if not (rec_col_norm and rec_col_norm != "unknown" and s_n == rec_col_norm)
        ]

        raw_type = row.get(type_col) if type_col else None
        is_type, _ = _is_verified_type_status(raw_type)

        if is_type:
            if valid_inst and inst != "UNKNOWN" and std_sp:
                inst_type_sheets[inst] = inst_type_sheets.get(inst, 0) + 1
            if has_explicit_det_str:
                for s_name in cross_collector_surnames:
                    det_type_annotated.add(s_name)
                    det_type_hits[s_name] = det_type_hits.get(s_name, 0) + 1

        if valid_inst:
            inst_total_sheets[inst] = inst_total_sheets.get(inst, 0) + 1
            if std_sp:
                it_key = (inst, std_sp)
                inst_taxon_totals[it_key] = inst_taxon_totals.get(it_key, 0) + 1
                if has_attributed_det:
                    inst_taxon_attributed[it_key] = inst_taxon_attributed.get(it_key, 0) + 1

        if std_sp and cross_collector_surnames:
            det_weight = 2 if (is_type and has_explicit_det_str) else 1
            for s_name in cross_collector_surnames:
                dt_key = (s_name, std_sp)
                det_taxon_counts[dt_key] = det_taxon_counts.get(dt_key, 0) + det_weight
                det_total_counts[s_name] = det_total_counts.get(s_name, 0) + det_weight
                det_species.setdefault(s_name, set()).add(std_sp)
                if valid_inst:
                    det_taxon_insts.setdefault(dt_key, set()).add(inst)
                    det_total_insts.setdefault(s_name, set()).add(inst)

        dup_key = ""
        if has_exsiccatae_col:
            raw_ek = row.get("exsiccatae_key")
            if raw_ek is not None and not pd.isna(raw_ek):
                dup_key = str(raw_ek).strip()
        if not dup_key and "recordNumber" in row:
            r_num = row.get("recordNumber", "")
            norm_num = normalize_record_number("" if pd.isna(r_num) else r_num)
            if norm_num:
                r_by = row.get("recordedBy", "")
                norm_col = normalize_collector("" if pd.isna(r_by) else str(r_by))
                ev_dt = row.get("eventDate", "")
                norm_dt = normalize_event_date(
                    "" if pd.isna(ev_dt) else str(ev_dt),
                    year=row.get("year"),
                    month=row.get("month"),
                    day=row.get("day"),
                )
                dup_key = f"{norm_col}_{norm_num}_{norm_dt}"
        if dup_key and std_sp:
            pre_conf = bool(
                "exsiccatae_taxonomic_conflict" in row
                and not pd.isna(row.get("exsiccatae_taxonomic_conflict"))
                and bool(row.get("exsiccatae_taxonomic_conflict"))
            )
            group_entries.setdefault(dup_key, []).append((std_sp, surnames, pre_conf))

    # Build cross-duplicate conflict counts per determiner before promotion
    det_dup_groups: Dict[str, Set[str]] = {}
    det_conflicted_groups: Dict[str, Set[str]] = {}
    for dup_key, entries in group_entries.items():
        if len(entries) < 2:
            continue
        distinct_group_taxa = {e[0] for e in entries if e[0]}
        is_conflicted_group = (len(distinct_group_taxa) >= 2) or any(e[2] for e in entries)
        for _, entry_surnames, _ in entries:
            for s_name in entry_surnames:
                det_dup_groups.setdefault(s_name, set()).add(dup_key)
                if is_conflicted_group:
                    det_conflicted_groups.setdefault(s_name, set()).add(dup_key)

    def _passes_specialist_gates(s_name: str) -> bool:
        type_hits = det_type_hits.get(s_name, 0)
        if s_name in _AMBIGUOUS_COMMON_SURNAMES and s_name not in PACKERA_SPECIALISTS and type_hits < 1:
            return False
        dup_groups = len(det_dup_groups.get(s_name, set()))
        conflicted_groups = len(det_conflicted_groups.get(s_name, set()))
        if dup_groups >= 2 and (conflicted_groups / float(dup_groups)) > 0.15:
            return False
        if is_large_corpus and type_hits < 1 and len(det_species.get(s_name, set())) < 4:
            return False
        return True

    for dt_key, count in det_taxon_counts.items():
        s_name = dt_key[0]
        if not _passes_specialist_gates(s_name):
            continue
        inst_set = det_taxon_insts.get(dt_key, set())
        if (count >= eff_min_det and len(inst_set) >= eff_min_herbs) or (
            s_name in det_type_annotated and count >= 3 and len(inst_set) >= 2
        ):
            dynamic_specialists.add(s_name)

    for s_name, tot_count in det_total_counts.items():
        if s_name in dynamic_specialists or not _passes_specialist_gates(s_name):
            continue
        tot_insts = det_total_insts.get(s_name, set())
        if (tot_count >= eff_min_det and len(tot_insts) >= eff_min_herbs) or (
            s_name in det_type_annotated and tot_count >= 3 and len(tot_insts) >= 2
        ):
            dynamic_specialists.add(s_name)

    if inst_total_sheets:
        for inst, n_types in inst_type_sheets.items():
            if n_types >= 1 and inst_total_sheets.get(inst, 0) >= 5:
                dynamic_herbaria.add(inst)
        q75_holdings = float(np.percentile(list(inst_total_sheets.values()), 75))
        for inst, total_n in inst_total_sheets.items():
            if total_n >= eff_min_herb_vouchers and total_n >= q75_holdings:
                min_sp_sheets = min(3, total_n)
                for (i_code, sp_name), sp_total in inst_taxon_totals.items():
                    if i_code == inst and sp_total >= min_sp_sheets:
                        n_attr = inst_taxon_attributed.get((i_code, sp_name), 0)
                        if (n_attr / float(sp_total)) >= 0.70:
                            dynamic_herbaria.add(inst)
                            break

    return dynamic_specialists, dynamic_herbaria
