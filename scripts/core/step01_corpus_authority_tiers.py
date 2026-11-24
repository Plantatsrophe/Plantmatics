"""step01_corpus_authority_tiers.py - 3-layer determiner tier scoring and corpus-level authority application."""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.config import MAJOR_HERBARIA_CODES
from scripts.core.step01_concept_authority_utils import apply_concept_and_segregate_recency_caps
from scripts.core.step01_corpus_authority_discovery import discover_corpus_authorities
from scripts.core.step01_determiner_surnames import (
    GENERIC_DETERMINER_TOKENS,
    PACKERA_SPECIALISTS,
    _AMBIGUOUS_COMMON_SURNAMES,
    _INVALID_INSTITUTION_CODES,
    _fold_ascii,
    _has_core_monographer_match,
    _has_non_person_determination_token,
    _has_regional_or_historical_specialist_match,
    _has_specialist_match,
    _is_explicit_collector_reference,
    _is_missing_identified_by,
    _should_fallback_to_collector,
    classify_nomenclatural_type_rank,
    normalize_determiner_surnames,
)
from scripts.core.step01_exsiccatae_fingerprint import (
    _as_bool_flag,
    _resolve_row_standardized_packera,
    normalize_collector,
)


def parse_determiner_tier(
    type_status_raw: Optional[str],
    identified_by_raw: Optional[str],
    recorded_by_raw: Optional[str],
    history_raw: Optional[str],
    institution_code_raw: Optional[str],
    locality_raw: Optional[str],
    habitat_raw: Optional[str],
    dynamic_specialists: Optional[Iterable[str]] = None,
    dynamic_herbaria: Optional[Iterable[str]] = None,
    is_open_or_hybrid: bool = False,
) -> Tuple[str, str, str]:
    """
    Evaluates Darwin Core fields to assign a taxonomic determination credibility tier
    using the 3-Layer Dynamic Expert & Herbarium Authority Engine:
      - Layer 1: Genus-wide Packera / Senecioneae monographer patterns (`SPECIALIST_PATTERNS`).
      - Layer 2 & 3: Dynamic protologue authors and corpus-derived multi-institution
        specialists (`dynamic_specialists`) and regional research herbaria (`dynamic_herbaria`).

    Tiers:
      - Tier_1_Gold: Primary/secondary nomenclatural types, monographers, protologue authors,
        or corpus-verified multi-herbarium specialists.
      - Tier_2_Silver: Determinations from major or corpus-qualified research herbaria with
        complete ecological locality and attributed determiner.
      - Tier_3_Bronze: Unverified candidate collections, general collectors, or missing determiner.

    Returns:
        Tuple[str, str, str]: (determiner_tier, type_status_clean, identified_by_combined)
    """
    type_status = "" if (type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw))) else str(type_status_raw).strip()
    identified_by = "" if (identified_by_raw is None or (isinstance(identified_by_raw, float) and math.isnan(identified_by_raw))) else str(identified_by_raw).strip()
    recorded_by = "" if (recorded_by_raw is None or (isinstance(recorded_by_raw, float) and math.isnan(recorded_by_raw))) else str(recorded_by_raw).strip()
    history = "" if (history_raw is None or (isinstance(history_raw, float) and math.isnan(history_raw))) else str(history_raw).strip()
    institution = "" if (institution_code_raw is None or (isinstance(institution_code_raw, float) and math.isnan(institution_code_raw))) else str(institution_code_raw).strip().upper()
    locality = "" if (locality_raw is None or (isinstance(locality_raw, float) and math.isnan(locality_raw))) else str(locality_raw).strip()
    habitat = "" if (habitat_raw is None or (isinstance(habitat_raw, float) and math.isnan(habitat_raw))) else str(habitat_raw).strip()

    type_rank, type_status_clean = classify_nomenclatural_type_rank(type_status)
    is_type = type_rank > 0

    is_non_person_det = _has_non_person_determination_token(identified_by)
    has_collector_ref = False if is_non_person_det else _is_explicit_collector_reference(identified_by)
    has_explicit_det = bool(
        not is_non_person_det
        and identified_by
        and not _is_missing_identified_by(identified_by)
        and not has_collector_ref
    )
    fallback_to_collector = (
        False
        if (is_non_person_det and type_rank <= 0)
        else _should_fallback_to_collector(identified_by)
    )

    if has_collector_ref and recorded_by:
        identified_by_combined = recorded_by
    elif has_explicit_det:
        identified_by_combined = identified_by
    else:
        identified_by_combined = history if history else (recorded_by if recorded_by else identified_by)

    # Quarantine non-type open-nomenclature and hybrid-formula sheets strictly to Tier_3_Bronze
    if is_open_or_hybrid and not is_type:
        return "Tier_3_Bronze", type_status_clean, identified_by_combined

    # Match specialists strictly on current determination (`identifiedBy`), only consulting
    # `verbatimIdentificationHistory` when `identifiedBy` is absent/collector-ref (so an older
    # superseded specialist annotation in history never promotes a newer non-specialist override).
    eff_id_for_specialist = "" if is_non_person_det else identified_by
    eff_hist_for_specialist = (
        ""
        if (has_explicit_det or _has_non_person_determination_token(history))
        else history
    )

    is_core_monographer = False
    is_regional_specialist_silver = False
    if not (is_non_person_det and not eff_hist_for_specialist and type_rank <= 0):
        explicit_det_text = f"{eff_id_for_specialist} | {eff_hist_for_specialist}".strip(" |")
        if explicit_det_text and _has_core_monographer_match(explicit_det_text):
            is_core_monographer = True
        elif explicit_det_text and (
            _has_regional_or_historical_specialist_match(explicit_det_text)
            or _has_specialist_match(explicit_det_text)
        ):
            is_regional_specialist_silver = True
        elif fallback_to_collector and recorded_by:
            if _has_core_monographer_match(recorded_by):
                is_core_monographer = True
            elif _has_regional_or_historical_specialist_match(recorded_by) or _has_specialist_match(recorded_by):
                # Regional flora collectors on blank-identifiedBy field labels route to Silver
                is_regional_specialist_silver = True

    # Corpus-discovered multi-herbarium specialists cap at Tier_2_Silver (never Tier_1_Gold),
    # and apply ONLY to explicit determinations/unoverridden history (never collector fallback).
    is_dynamic_corpus_specialist_silver = False
    if not is_core_monographer and not is_regional_specialist_silver and dynamic_specialists and (has_explicit_det or eff_hist_for_specialist):
        dyn_spec_set: Set[str] = set()
        for item in dynamic_specialists:
            if item is None:
                continue
            s_tok = str(item).strip().lower()
            if len(s_tok) >= 3 and s_tok not in _AMBIGUOUS_COMMON_SURNAMES:
                dyn_spec_set.add(s_tok)
                dyn_spec_set.add(_fold_ascii(s_tok))
        if dyn_spec_set:
            det_sources = [eff_id_for_specialist if has_explicit_det else "", eff_hist_for_specialist]
            det_surnames: Set[str] = set()
            for src in det_sources:
                det_surnames.update(normalize_determiner_surnames(src))
            if det_surnames & dyn_spec_set:
                is_dynamic_corpus_specialist_silver = True
            else:
                det_audit_raw = " | ".join(det_sources).lower()
                det_audit_ascii = _fold_ascii(det_audit_raw)
                for spec_tok in dyn_spec_set:
                    pat = rf"(?<![a-zA-ZÀ-ÿ]){re.escape(spec_tok)}(?![a-zA-ZÀ-ÿ])"
                    if re.search(pat, det_audit_raw) or re.search(pat, det_audit_ascii):
                        is_dynamic_corpus_specialist_silver = True
                        break

    if is_type or is_core_monographer:
        return "Tier_1_Gold", type_status_clean, identified_by_combined

    if is_regional_specialist_silver or is_dynamic_corpus_specialist_silver:
        return "Tier_2_Silver", type_status_clean, identified_by_combined

    if "[exsiccatae " in identified_by.lower():
        return "Tier_2_Silver", type_status_clean, identified_by_combined

    dyn_herb_set = {str(h).strip().upper() for h in dynamic_herbaria if h} if dynamic_herbaria else set()
    is_major_herbarium = (institution in MAJOR_HERBARIA_CODES) or (institution in dyn_herb_set)
    has_rich_locality = len(locality) > 5 or len(habitat) > 3
    eff_silver_det = "" if is_non_person_det else (recorded_by if (has_collector_ref and recorded_by) else identified_by)
    silver_low = eff_silver_det.lower()
    is_student_or_intern = bool(re.search(r"\b(?:student|undergrad(?:uate)?|intern|class|volunteer)\b", silver_low))
    has_determiner = bool(
        eff_silver_det
        and not is_student_or_intern
        and silver_low not in GENERIC_DETERMINER_TOKENS
        and not _has_non_person_determination_token(eff_silver_det)
        and len(normalize_determiner_surnames(eff_silver_det)) > 0
    )

    if is_major_herbarium and has_rich_locality and has_determiner:
        return "Tier_2_Silver", type_status_clean, identified_by_combined

    return "Tier_3_Bronze", type_status_clean, identified_by_combined


def apply_corpus_authority_tiers(
    df: pd.DataFrame,
    dynamic_specialists: Optional[Set[str]] = None,
    dynamic_herbaria: Optional[Set[str]] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Applies the 3-Layer Dynamic Authority Engine across an assembled DataFrame:
      - Discovers corpus-level multi-herbarium specialists and dynamic regional herbaria.
      - Collects per-taxon protologue authors from Darwin Core `scientificName` / `species_raw`.
      - Detects undated institutional batch transfers (>= 25 sheets at 1 institution across >= 10
        collectors) and demotes non-type historical sheets collected by others from `Tier_1_Gold`
        to `Tier_2_Silver` for out-of-sample screening, while exempting recognized *Packera*
        specialists who annotate across multiple herbaria or multiple species in a batch.
      - Re-evaluates `determiner_tier` across all rows cleanly without downgrading pre-scored
        Silver rows whose `locality` / `habitat` columns were already omitted upon CSV persistence.
    """
    if df is None or df.empty:
        return df

    out = df.copy()
    if dynamic_specialists is None or dynamic_herbaria is None:
        disc_specs, disc_herbs = discover_corpus_authorities(out)
        eff_specs = disc_specs if dynamic_specialists is None else (set(dynamic_specialists) | disc_specs)
        eff_herbs = disc_herbs if dynamic_herbaria is None else (set(dynamic_herbaria) | disc_herbs)
    else:
        eff_specs = set(dynamic_specialists)
        eff_herbs = set(dynamic_herbaria)

    from scripts.core.harvester import standardize_packera_taxon

    # Tally (institutionCode, identifiedBy) batch statistics to detect undated institutional batch transfers
    batch_voucher_counts: Dict[Tuple[str, str], int] = {}
    batch_collectors: Dict[Tuple[str, str], Set[str]] = {}
    batch_dates: Dict[Tuple[str, str], Set[str]] = {}
    det_institutions: Dict[str, Set[str]] = {}
    surname_institutions: Dict[str, Set[str]] = {}
    key_species_counts: Dict[Tuple[str, str], Dict[str, int]] = {}

    for _, row in out.iterrows():
        raw_names_row = [
            row.get("scientificName"),
            row.get("species_raw"),
        ]
        std_sp = None
        for rn in raw_names_row:
            std_sp = standardize_packera_taxon(rn)
            if std_sp:
                break
        if not std_sp:
            std_sp = _resolve_row_standardized_packera(row)

        inst_r = row.get("institutionCode")
        inst_key = "" if (inst_r is None or pd.isna(inst_r)) else str(inst_r).strip().upper()
        id_r = row.get("_raw_identified_by") if "_raw_identified_by" in out.columns else row.get("identifiedBy")
        id_str_clean = "" if (id_r is None or pd.isna(id_r)) else str(id_r).strip()
        id_key = re.sub(r"\s+", " ", id_str_clean.lower()) if id_str_clean else ""
        if (
            inst_key
            and inst_key not in _INVALID_INSTITUTION_CODES
            and id_key
            and not _is_missing_identified_by(id_key)
            and not _is_explicit_collector_reference(id_key)
            and not _has_non_person_determination_token(id_key)
        ):
            b_key = (inst_key, id_key)
            batch_voucher_counts[b_key] = batch_voucher_counts.get(b_key, 0) + 1
            det_institutions.setdefault(id_key, set()).add(inst_key)
            for s_n in normalize_determiner_surnames(id_str_clean):
                surname_institutions.setdefault(s_n, set()).add(inst_key)
            sp_val = str(std_sp or "").strip().lower()
            if sp_val:
                sp_map = key_species_counts.setdefault(b_key, {})
                sp_map[sp_val] = sp_map.get(sp_val, 0) + 1
            rec_r = row.get("recordedBy")
            col_norm = normalize_collector("" if (rec_r is None or pd.isna(rec_r)) else str(rec_r))
            if col_norm and col_norm != "unknown":
                batch_collectors.setdefault(b_key, set()).add(col_norm)
            dt_id_r = row.get("_raw_date_identified") if "_raw_date_identified" in out.columns else row.get("dateIdentified")
            dt_id_str = "" if (dt_id_r is None or pd.isna(dt_id_r)) else str(dt_id_r).strip()
            if dt_id_str and dt_id_str.lower() not in {"none", "nan", "null", "na", "n/a", "unknown"}:
                batch_dates.setdefault(b_key, set()).add(dt_id_str)

    flagged_batch_keys: Set[Tuple[str, str]] = set()
    for b_key, v_cnt in batch_voucher_counts.items():
        if v_cnt < 25 or len(batch_collectors.get(b_key, set())) < 10:
            continue
        _inst_k, _id_k = b_key
        if _has_specialist_match(_id_k):
            inst_set_for_det = set(det_institutions.get(_id_k, set()))
            for _sn in normalize_determiner_surnames(_id_k):
                if _sn in PACKERA_SPECIALISTS:
                    inst_set_for_det |= surname_institutions.get(_sn, set())
            n_inst_for_det = len(inst_set_for_det)
            sp_counts = key_species_counts.get(b_key, {})
            n_species_in_batch = len(sp_counts)
            max_sp_share = (
                (max(sp_counts.values()) / max(1, sum(sp_counts.values())))
                if sp_counts
                else 1.0
            )
            if n_inst_for_det >= 2 or (n_species_in_batch >= 2 and max_sp_share < 0.85):
                continue
        flagged_batch_keys.add(b_key)

    has_determiner_cols = any(
        c in out.columns
        for c in (
            "_raw_identified_by",
            "identifiedBy",
            "recordedBy",
            "_raw_type_status",
            "type_status",
            "typeStatus",
            "_raw_history",
            "verbatimIdentificationHistory",
        )
    )
    if not has_determiner_cols and "determiner_tier" in out.columns:
        return out

    has_locality_cols = any(
        c in out.columns for c in ("_raw_locality", "locality", "verbatimLocality", "_raw_habitat", "habitat")
    )

    new_tiers: List[str] = []
    new_types: List[str] = []
    update_type_col = "type_status" in out.columns

    for _, row in out.iterrows():
        t_raw = row.get("_raw_type_status") if "_raw_type_status" in out.columns else row.get("type_status", row.get("typeStatus"))
        id_raw = row.get("_raw_identified_by") if "_raw_identified_by" in out.columns else row.get("identifiedBy")
        rec_raw = row.get("recordedBy")
        hist_raw = row.get("_raw_history") if "_raw_history" in out.columns else row.get("verbatimIdentificationHistory")
        inst_raw = row.get("institutionCode")
        loc_raw = (
            row.get("_raw_locality")
            if "_raw_locality" in out.columns
            else (row.get("locality") if "locality" in out.columns else row.get("verbatimLocality"))
        )
        hab_raw = row.get("_raw_habitat") if "_raw_habitat" in out.columns else row.get("habitat")

        raw_sci = row.get("scientificName") or row.get("species_raw")
        std_sp = standardize_packera_taxon(raw_sci) or _resolve_row_standardized_packera(row)
        c_stat_r = str(row.get("concept_synonymy_status", "") if not pd.isna(row.get("concept_synonymy_status", "")) else "").strip()
        row_open_hyb = _as_bool_flag(row.get("is_open_or_hybrid_concept")) or (c_stat_r in {"Open_Nomenclature", "Hybrid_Formula"})

        tier_eval, type_clean, _ = parse_determiner_tier(
            type_status_raw=t_raw,
            identified_by_raw=id_raw,
            recorded_by_raw=rec_raw,
            history_raw=hist_raw,
            institution_code_raw=inst_raw,
            locality_raw=loc_raw,
            habitat_raw=hab_raw,
            dynamic_specialists=eff_specs,
            dynamic_herbaria=eff_herbs,
            is_open_or_hybrid=row_open_hyb,
        )

        prev_tier = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
        ex_gold = str(row.get("exsiccatae_gold_taxon", "") if not pd.isna(row.get("exsiccatae_gold_taxon", "")) else "").strip()
        is_exsiccatae_transfer = "[exsiccatae " in str(id_raw or "").lower()
        t_rank_row, _ = classify_nomenclatural_type_rank(t_raw)
        if row_open_hyb and t_rank_row == 0:
            tier_eval = "Tier_3_Bronze"
        elif is_exsiccatae_transfer:
            tier_eval = "Tier_2_Silver"
        elif prev_tier in {"Tier_1_Gold", "Tier_2_Silver"} and ex_gold and std_sp == ex_gold and tier_eval == "Tier_3_Bronze":
            tier_eval = prev_tier
        elif not has_locality_cols and prev_tier == "Tier_2_Silver" and tier_eval == "Tier_3_Bronze":
            tier_eval = "Tier_2_Silver"

        # Demote undated (or batch-identical-dated) institutional batch transfers on historical sheets
        # collected by others (type_rank == 0) from Tier_1_Gold to Tier_2_Silver
        if tier_eval == "Tier_1_Gold" and not (ex_gold and std_sp == ex_gold):
            t_rank, _ = classify_nomenclatural_type_rank(t_raw)
            if t_rank == 0:
                inst_key = "" if (inst_raw is None or pd.isna(inst_raw)) else str(inst_raw).strip().upper()
                id_str_row = "" if (id_raw is None or pd.isna(id_raw)) else str(id_raw).strip()
                id_key = re.sub(r"\s+", " ", id_str_row.lower()) if id_str_row else ""
                b_key = (inst_key, id_key)
                if b_key in flagged_batch_keys:
                    dt_id_r = row.get("_raw_date_identified") if "_raw_date_identified" in out.columns else row.get("dateIdentified")
                    dt_id_str = "" if (dt_id_r is None or pd.isna(dt_id_r)) else str(dt_id_r).strip()
                    has_embedded_year = bool(id_str_row and re.search(r"\b(19|20)\d{2}\b", id_str_row))
                    has_date = bool(dt_id_str and re.search(r"\b(19|20)\d{2}\b", dt_id_str)) or has_embedded_year
                    has_blank_or_batch_date = (
                        not has_date
                        or (not has_embedded_year and len(batch_dates.get(b_key, set())) <= 1)
                    )
                    if has_blank_or_batch_date:
                        det_surnames = set(normalize_determiner_surnames(id_raw))
                        rec_str = "" if (rec_raw is None or pd.isna(rec_raw)) else str(rec_raw)
                        rec_surnames = set(normalize_determiner_surnames(rec_str))
                        rec_primary = normalize_collector(rec_str)
                        if rec_primary and rec_primary != "unknown":
                            rec_surnames.add(rec_primary)
                        rec_ascii = _fold_ascii(rec_str)
                        collector_matches_det = bool(det_surnames & rec_surnames) or any(
                            bool(re.search(rf"(?<![a-zA-ZÀ-ÿ]){re.escape(ds)}(?![a-zA-ZÀ-ÿ])", rec_ascii))
                            for ds in det_surnames
                        )
                        if not collector_matches_det:
                            tier_eval = "Tier_2_Silver"

        new_tiers.append(tier_eval)
        prev_type = str(row.get("type_status", "") if not pd.isna(row.get("type_status", "")) else "").strip()
        new_types.append(type_clean if (type_clean != "None" or not prev_type) else prev_type)

    out["determiner_tier"] = new_tiers
    if update_type_col:
        out["type_status"] = new_types
    return apply_concept_and_segregate_recency_caps(out, db_path=db_path, workspace_dir=workspace_dir)
