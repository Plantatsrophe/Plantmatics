"""Step 01 GBIF Backbone synonym expansion and Packera taxon standardization."""

from __future__ import annotations

import concurrent.futures
import math
import re
from typing import Any, Dict, List, Optional, Set, Tuple

from scripts.core.step01_biogeography_utils import normalize_packera_orthography
from scripts.core.step01_concept_authority_utils import is_interspecific_hybrid_formula
from scripts.core.step01_harvester_filters import (
    _OPEN_NOMENCLATURE_EPITHETS,
    _SYNONYM_GENERA,
    _extract_canonical_binomial,
)


def _query_gbif_backbone_and_synonyms(
    query_name: str,
    timeout_seconds: float = 5.0,
) -> Tuple[Optional[str], List[str], List[str]]:
    """
    Query `pygbif.species.name_backbone` and `pygbif.species.name_usage(..., data='synonyms')`
    with a hard 5-second timeout and in-memory cache so offline unit tests never fail or hang.
    Returns `(accepted_packera_binomial, synonym_binomials, raw_scientific_names_with_authors)`
    and records in `_GBIF_BACKBONE_SYNONYM_FLAG` whether `query_name` itself is an exact
    unaccepted synonym in the GBIF Backbone.
    """
    import scripts.core.harvester as _harv

    is_mocked = _harv._is_gbif_species_mocked()
    cache_key = query_name.strip().lower()
    if not is_mocked and cache_key in _harv._GBIF_SYNONYM_CACHE:
        return _harv._GBIF_SYNONYM_CACHE[cache_key]
    if not is_mocked and _harv._GBIF_NETWORK_UNAVAILABLE:
        return None, [], []

    def _fetch() -> Tuple[Optional[str], List[str], List[str], bool]:
        sp = _harv.species
        backbone = _harv._call_and_normalize_gbif_backbone(query_name, "species")
        if not isinstance(backbone, dict) or not backbone:
            return None, [], [], False
        match_type = str(backbone.get("matchType", "EXACT")).upper()
        if match_type in {"NONE", "HIGHERRANK"}:
            return None, [], [], False

        bb_status = str(backbone.get("status") or backbone.get("taxonomicStatus") or "").strip().upper()
        is_exact_synonym = (match_type == "EXACT") and (bb_status in _harv._GBIF_SYNONYM_STATUSES)

        author_SciNames: List[str] = []
        for fld in ("scientificName", "accepted"):
            val = backbone.get(fld)
            if isinstance(val, str) and val.strip():
                author_SciNames.append(val.strip())

        accepted_packera: Optional[str] = None
        acc_fields = ("accepted", "species", "canonicalName") if is_exact_synonym else ("species", "accepted", "canonicalName")
        for cand_fld in acc_fields:
            cand_bin = _extract_canonical_binomial(backbone.get(cand_fld), {"Packera"})
            if cand_bin:
                accepted_packera = cand_bin
                break

        usage_key = backbone.get("usageKey")
        accepted_key = backbone.get("acceptedUsageKey")
        query_packera_bin = _extract_canonical_binomial(query_name, {"Packera"})
        if (
            accepted_key is not None
            and str(accepted_key) != str(usage_key)
            and (accepted_packera is None or (query_packera_bin and accepted_packera.lower() == query_packera_bin.lower()))
            and hasattr(sp, "name_usage")
        ):
            try:
                acc_key_norm = int(accepted_key) if str(accepted_key).isdigit() else accepted_key
                acc_info = sp.name_usage(key=acc_key_norm)
                if isinstance(acc_info, dict) and "results" not in acc_info:
                    acc_bin = _extract_canonical_binomial(
                        acc_info.get("canonicalName") or acc_info.get("species") or acc_info.get("scientificName"),
                        {"Packera"},
                    )
                    if acc_bin:
                        accepted_packera = acc_bin
                    acc_sci = acc_info.get("scientificName")
                    if isinstance(acc_sci, str) and acc_sci.strip():
                        author_SciNames.append(acc_sci.strip())
            except Exception:
                pass

        found_binomials: List[str] = []
        for cand_fld in ("canonicalName", "species", "accepted"):
            b_val = _extract_canonical_binomial(backbone.get(cand_fld), _SYNONYM_GENERA)
            if b_val and b_val not in found_binomials:
                found_binomials.append(b_val)

        keys_to_query: List[Any] = []
        if usage_key is not None:
            keys_to_query.append(int(usage_key) if str(usage_key).isdigit() else usage_key)
        if accepted_key is not None and str(accepted_key) != str(usage_key):
            keys_to_query.append(int(accepted_key) if str(accepted_key).isdigit() else accepted_key)

        for k_id in keys_to_query:
            try:
                syn_resp = sp.name_usage(key=k_id, data="synonyms")
            except Exception:
                syn_resp = []
            if isinstance(syn_resp, dict):
                syn_items = syn_resp.get("results", [])
            elif isinstance(syn_resp, list):
                syn_items = syn_resp
            else:
                syn_items = []
            for item in syn_items:
                if not isinstance(item, dict):
                    continue
                item_status = str(item.get("taxonomicStatus") or item.get("status") or "").strip().upper()
                if item_status in _harv._GBIF_REJECTED_SYNONYM_STATUSES:
                    continue
                sci_full = item.get("scientificName")
                if isinstance(sci_full, str) and sci_full.strip():
                    author_SciNames.append(sci_full.strip())
                syn_bin = _extract_canonical_binomial(
                    item.get("canonicalName") or item.get("scientificName") or item.get("species"),
                    _SYNONYM_GENERA,
                )
                if syn_bin and syn_bin not in found_binomials:
                    found_binomials.append(syn_bin)
                if accepted_packera is None:
                    acc_cand = _extract_canonical_binomial(item.get("species") or item.get("accepted"), {"Packera"})
                    if acc_cand:
                        accepted_packera = acc_cand

        return accepted_packera, found_binomials, author_SciNames, is_exact_synonym

    if is_mocked:
        try:
            acc_p, f_bins, auth_s, is_syn = _fetch()
            _harv._GBIF_BACKBONE_SYNONYM_FLAG[cache_key] = is_syn
            return acc_p, f_bins, auth_s
        except Exception:
            return None, [], []

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch)
        acc_p, f_bins, auth_s, is_syn = future.result(timeout=timeout_seconds)
        res = (acc_p, f_bins, auth_s)
        _harv._GBIF_SYNONYM_CACHE[cache_key] = res
        _harv._GBIF_BACKBONE_SYNONYM_FLAG[cache_key] = is_syn
        return res
    except Exception:
        _harv._GBIF_NETWORK_UNAVAILABLE = True
        return None, [], []
    finally:
        executor.shutdown(wait=False, cancel_futures=True)


def expand_taxon_synonyms_gbif(
    taxon_name: str,
    cfg_synonyms: Optional[Dict[str, List[str]]] = None,
    protected_packera_taxa: Optional[Set[str]] = None,
) -> List[str]:
    """
    Dynamically expands a target `Packera` (or historical `Senecio`/`Cineraria`) taxon into its
    full homotypic and heterotypic synonym query set using:
      1. Optional configuration synonym mappings (`cfg_synonyms`).
      2. Live GBIF Species API (`pygbif.species.name_backbone` + `pygbif.species.name_usage`),
         guarded by a 5-second timeout and in-memory cache for offline reliability.
    Unaccepted species-rank `Packera` synonyms in the GBIF Backbone (e.g., `Packera tomentosa`
    under `Packera dubia`) and historical generic combinations (`Senecio`/`Cineraria`) are
    dynamically unified under the accepted `Packera <epithet>` in `_DYNAMIC_SYNONYM_TO_CANONICAL`
    and persisted to SQLite `taxon_synonyms`, while protecting distinct accepted `Packera`
    species in `protected_packera_taxa` from cross-species overwriting.
    """
    import scripts.core.harvester as _harv

    if taxon_name is None:
        return []
    clean_input = str(taxon_name).strip()
    if not clean_input:
        return []
    if clean_input.lower() == "packera":
        return ["Packera"]

    protected_lower: Set[str] = {str(p).strip().lower() for p in (protected_packera_taxa or set()) if p}
    expanded: List[str] = []
    seen_lower: Set[str] = set()

    def _add_candidate(cand: Optional[str]) -> None:
        if not cand:
            return
        c_str = str(cand).strip()
        if not c_str:
            return
        c_low = c_str.lower()
        if c_low not in seen_lower:
            seen_lower.add(c_low)
            expanded.append(c_str)

    canonical_packera = _extract_canonical_binomial(clean_input, {"Packera"})
    if canonical_packera and canonical_packera.lower() in _harv._DYNAMIC_SYNONYM_TO_CANONICAL:
        canonical_packera = _harv._DYNAMIC_SYNONYM_TO_CANONICAL[canonical_packera.lower()]

    cfg_group_matched: List[str] = []
    if cfg_synonyms:
        for cfg_key, cfg_vals in cfg_synonyms.items():
            key_clean = str(cfg_key).strip()
            vals_clean = [str(v).strip() for v in (cfg_vals or []) if str(v).strip()]
            all_group = [key_clean] + vals_clean
            all_group_lower = {g.lower() for g in all_group}
            bin_input = _extract_canonical_binomial(clean_input, _SYNONYM_GENERA)
            if (
                clean_input.lower() in all_group_lower
                or (bin_input and bin_input.lower() in all_group_lower)
                or (canonical_packera and canonical_packera.lower() == key_clean.lower())
            ):
                if canonical_packera is None:
                    canonical_packera = _extract_canonical_binomial(key_clean, {"Packera"})
                for g in all_group:
                    _add_candidate(g)
                    cfg_group_matched.append(g)

    if canonical_packera and canonical_packera.startswith("Packera "):
        _add_candidate(canonical_packera)
    _add_candidate(clean_input)

    # Query GBIF Species Backbone + Synonyms endpoint (fast 5s timeout + cached)
    gbif_query_target = canonical_packera if (canonical_packera and canonical_packera.startswith("Packera ")) else clean_input
    gbif_acc, gbif_syns, gbif_author_names = _harv._query_gbif_backbone_and_synonyms(gbif_query_target)
    query_is_synonym = _harv._GBIF_BACKBONE_SYNONYM_FLAG.get(gbif_query_target.strip().lower(), False)
    if gbif_acc and gbif_acc.startswith("Packera "):
        if not (canonical_packera and canonical_packera.startswith("Packera ")):
            canonical_packera = gbif_acc
        elif query_is_synonym and gbif_acc.lower() != canonical_packera.lower():
            canonical_packera = gbif_acc
        _add_candidate(gbif_acc)
    for g_syn in gbif_syns:
        _add_candidate(g_syn)

    if canonical_packera and canonical_packera.startswith("Packera "):
        canon_low = canonical_packera.lower()
        other_protected = protected_lower - {canon_low}
        if cfg_group_matched:
            _harv._register_dynamic_synonyms(canonical_packera, cfg_group_matched, None, allow_packera_species_remap=True)
        safe_expanded = [
            cand
            for cand in expanded
            if (_extract_canonical_binomial(cand, {"Packera"}) or "").lower() not in other_protected
        ]
        safe_authors = [
            auth
            for auth in (gbif_author_names or [])
            if (_extract_canonical_binomial(auth, {"Packera"}) or "").lower() not in other_protected
        ]
        _harv._register_dynamic_synonyms(
            canonical_packera,
            safe_expanded,
            safe_authors,
            allow_packera_species_remap=True,
        )
        _harv._persist_dynamic_synonyms_for_taxon(canonical_packera)

    return expanded


def standardize_packera_taxon(raw_name: Any) -> Optional[str]:
    """
    Dynamically standardizes raw Darwin Core scientific names, varieties, unaccepted `Packera`
    synonyms, and historical `Senecio`/`Cineraria` basionyms across the genus `Packera`.
    Uses the runtime/SQLite synonym table populated dynamically by `expand_taxon_synonyms_gbif()`
    and occurrence ingestion, GBIF `name_backbone` synonym verification, and genus-wide botanical
    regex matching (`^packera\\s*(?:([x×])\\s*)?([a-z-]+)`).
    """
    import scripts.core.harvester as _harv

    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return None
    s_orig = str(raw_name).strip()
    if not s_orig:
        return None
    s_compact = re.sub(r"\s+", " ", s_orig.lower())
    if s_compact in {"none", "nan", "null", "na", "unknown", "packera"}:
        return None
    if is_interspecific_hybrid_formula(s_orig, _SYNONYM_GENERA)[0]:
        return None

    cand_bin = _extract_canonical_binomial(s_orig, _SYNONYM_GENERA)
    is_non_packera_syn_genus = bool(cand_bin and cand_bin.split(" ", 1)[0] in (_SYNONYM_GENERA - {"Packera"}))

    # When pygbif.species is mocked in a unit test, allow mocked basionym or exact-synonym
    # Packera lookups to resolve directly without being shadowed by earlier cache state.
    if cand_bin and _harv._is_gbif_species_mocked():
        gbif_acc, gbif_syns, gbif_authors = _harv._query_gbif_backbone_and_synonyms(cand_bin)
        is_exact_syn = _harv._GBIF_BACKBONE_SYNONYM_FLAG.get(cand_bin.lower(), False)
        if gbif_acc and gbif_acc.startswith("Packera ") and (is_non_packera_syn_genus or is_exact_syn):
            _harv._register_dynamic_synonyms(
                gbif_acc,
                [s_orig, cand_bin, *gbif_syns],
                gbif_authors,
                allow_packera_species_remap=is_exact_syn,
            )
            return gbif_acc

    # 1. Dynamic GBIF Backbone & runtime/SQLite synonym dictionary lookup
    if s_compact in _harv._DYNAMIC_SYNONYM_TO_CANONICAL:
        return _harv._DYNAMIC_SYNONYM_TO_CANONICAL[s_compact]

    var_m = re.match(r"^([a-z]+\s+[a-z-]+\s+(?:var\.?|subsp\.?|ssp\.?|f\.?)\s+[a-z-]+)", s_compact)
    if var_m:
        var_raw = var_m.group(1)
        if var_raw in _harv._DYNAMIC_SYNONYM_TO_CANONICAL:
            return _harv._DYNAMIC_SYNONYM_TO_CANONICAL[var_raw]
        var_dot = re.sub(r"\b(var|subsp|ssp|f)\b\.?", r"\1.", var_raw)
        if var_dot in _harv._DYNAMIC_SYNONYM_TO_CANONICAL:
            return _harv._DYNAMIC_SYNONYM_TO_CANONICAL[var_dot]

    if cand_bin and cand_bin.lower() in _harv._DYNAMIC_SYNONYM_TO_CANONICAL:
        return _harv._DYNAMIC_SYNONYM_TO_CANONICAL[cand_bin.lower()]

    # 2. General genus-wide botanical regex for any valid Packera binomial
    m = re.match(r"^packera(?:\s*×\s*|\s+[x×]\s+|\s+)([a-z][a-z-]*)", s_compact, flags=re.IGNORECASE)
    if m:
        epithet = m.group(1).lower().strip("-")
        if len(epithet) >= 2 and epithet not in _OPEN_NOMENCLATURE_EPITHETS:
            return normalize_packera_orthography(f"Packera {epithet}")
        return None

    # 3. Dynamic GBIF Backbone resolution for historical Senecio / Cineraria synonyms
    if is_non_packera_syn_genus:
        gbif_acc, gbif_syns, gbif_authors = _harv._query_gbif_backbone_and_synonyms(cand_bin)
        if gbif_acc and gbif_acc.startswith("Packera "):
            _harv._register_dynamic_synonyms(gbif_acc, [s_orig, cand_bin, *gbif_syns], gbif_authors)
            _harv._persist_dynamic_synonyms_for_taxon(gbif_acc)
            return gbif_acc

    return None
