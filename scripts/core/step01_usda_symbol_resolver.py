"""
Module: step01_usda_symbol_resolver.py
Description:
    USDA PLANTS master checklist (`plantlst.txt`) parsing, Virginia independent city
    and East Texas county alias augmentation, algorithmic ICN generic transfers, and
    synonym-aware USDA accepted symbol resolution.
"""

from __future__ import annotations

import csv
import io
import logging
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

VA_INDEPENDENT_CITIES: Set[str] = {
    "chesapeake",
    "colonial heights",
    "hampton",
    "hopewell",
    "newport news",
    "norfolk",
    "petersburg",
    "poquoson",
    "portsmouth",
    "suffolk",
    "virginia beach",
    "williamsburg",
}
_EAST_TX_PINEY_WOODS_COUNTIES: Set[str] = {
    "smith",
    "upshur",
    "gregg",
    "harrison",
    "marion",
    "cass",
    "panola",
    "rusk",
    "cherokee",
    "nacogdoches",
    "anderson",
    "wood",
    "camp",
    "morris",
    "titus",
}


def normalize_virginia_independent_city(
    county: Optional[str],
    state: Optional[str] = "VA",
) -> Optional[str]:
    """
    Normalize Virginia independent cities lacking `"City"` when no Virginia county
    shares the same name, preserving Title Case `"<Name> City"`.
    """
    import scripts.core.usda_plants_ranges as _usda

    if county is None:
        return None
    raw_clean = re.sub(r"\s+", " ", str(county).strip())
    if not raw_clean:
        return None
    st_norm = _usda._normalize_us_state_name(state) if state is not None and str(state).strip() else None
    if st_norm == "virginia":
        stripped = re.sub(r"\s*\(independent\s+city\)\s*$", "", raw_clean, flags=re.IGNORECASE).strip()
        stripped = re.sub(r"\s+independent\s+city$", "", stripped, flags=re.IGNORECASE).strip()
        low = re.sub(r"[^a-z0-9]+", " ", stripped.lower()).strip()
        base_low = re.sub(r"\s+city$", "", low).strip()
        if base_low in VA_INDEPENDENT_CITIES:
            return f"{base_low.title()} City"
    return raw_clean


def _expand_virginia_city_aliases(counties: Set[str]) -> Set[str]:
    """Treat `'virginia|<city>'` and `'virginia|<city> city'` equivalently and link Chesapeake to Norfolk/Virginia Beach."""
    expanded: Set[str] = set(counties)
    for ck in list(counties):
        if not ck.startswith("virginia|"):
            continue
        co_part = ck.split("|", 1)[1].strip().lower()
        base_city = re.sub(r"\s+city$", "", co_part).strip()
        if base_city in VA_INDEPENDENT_CITIES:
            expanded.add(f"virginia|{base_city}")
            expanded.add(f"virginia|{base_city} city")
        if base_city in {"norfolk", "virginia beach"}:
            expanded.add("virginia|chesapeake")
            expanded.add("virginia|chesapeake city")
    return expanded


def _augment_documented_counties(
    canonical_packera: str,
    counties: Set[str],
    documented_states: Set[str],
) -> Set[str]:
    """Apply Virginia independent city aliases and documented East Texas Piney Woods augmentation."""
    augmented = _expand_virginia_city_aliases(counties)
    canon_str = str(canonical_packera or "").strip()
    if canon_str in {"Packera dubia", "Packera tomentosa"} and "texas" in documented_states:
        for tx_co in _EAST_TX_PINEY_WOODS_COUNTIES:
            augmented.add(f"texas|{tx_co}")
    return augmented


def _parse_plantlst_text(text: str) -> Dict[str, Any]:
    """
    Parse USDA PLANTS `plantlst.txt` (`"Symbol","Synonym Symbol","Scientific Name with Author",...`)
    into species-rank and infraspecific binomial -> accepted USDA `Symbol` lookup tables.
    """
    import scripts.core.usda_plants_ranges as _usda

    species_binom_to_symbols: Dict[str, List[str]] = {}
    infra_binom_to_symbols: Dict[str, List[str]] = {}
    symbol_to_accepted_name: Dict[str, str] = {}
    syn_symbol_to_accepted: Dict[str, List[str]] = {}
    all_symbols: Set[str] = set()

    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        if not isinstance(row, dict):
            continue
        sym = str(row.get("Symbol") or "").strip().upper()
        syn_sym = str(row.get("Synonym Symbol") or "").strip().upper()
        sci_auth = str(
            row.get("Scientific Name with Author")
            or row.get("Scientific Name with Authors")
            or row.get("Scientific Name")
            or ""
        ).strip()
        if not sym or not sci_auth:
            continue
        all_symbols.add(sym)
        if syn_sym:
            syn_bucket = syn_symbol_to_accepted.setdefault(syn_sym, [])
            if sym not in syn_bucket:
                syn_bucket.append(sym)
        binom, is_infra = _usda._extract_clean_binomial(sci_auth)
        if not binom:
            continue
        if not syn_sym and not is_infra and sym not in symbol_to_accepted_name:
            symbol_to_accepted_name[sym] = sci_auth
        target_map = infra_binom_to_symbols if is_infra else species_binom_to_symbols
        bucket = target_map.setdefault(binom, [])
        if sym not in bucket:
            # Prioritize accepted rows (empty Synonym Symbol) at front of list
            if not syn_sym:
                bucket.insert(0, sym)
            else:
                bucket.append(sym)

    return {
        "species_map": species_binom_to_symbols,
        "infra_map": infra_binom_to_symbols,
        "symbol_names": symbol_to_accepted_name,
        "syn_symbol_map": syn_symbol_to_accepted,
        "all_symbols": all_symbols,
    }


def load_or_fetch_plantlst_index(
    cache_dir: Optional[Path] = None,
    allow_network: bool = True,
    timeout: float = 20.0,
    logger: Optional[logging.Logger] = None,
) -> Optional[Dict[str, Any]]:
    """Load cached `plantlst.txt` index or download once from USDA PLANTS if allowed."""
    import scripts.core.usda_plants_ranges as _usda

    eff_dir = Path(cache_dir) if cache_dir is not None else _usda.DEFAULT_USDA_CACHE_DIR
    plantlst_path = eff_dir / "plantlst.txt"

    if not plantlst_path.is_file() and allow_network and not _usda._USDA_NETWORK_UNAVAILABLE:
        if logger:
            logger.info(f"Downloading USDA PLANTS master checklist to cache: {plantlst_path} ...")
        fetched = _usda._http_get_text(_usda.USDA_PLANTLST_URL, timeout=timeout)
        if fetched and "Symbol" in fetched[:256]:
            try:
                _usda.atomic_write_text(plantlst_path, fetched)
            except Exception as exc:
                if logger:
                    logger.debug(f"Could not cache plantlst.txt at {plantlst_path}: {exc}")
                return _usda._parse_plantlst_text(fetched)

    if not plantlst_path.is_file() or plantlst_path.stat().st_size == 0:
        return None

    mtime = plantlst_path.stat().st_mtime
    cache_key = (str(plantlst_path.resolve()), mtime)
    if cache_key in _usda._PLANTLST_INDEX_CACHE:
        return _usda._PLANTLST_INDEX_CACHE[cache_key]

    try:
        text = plantlst_path.read_text(encoding="utf-8-sig", errors="replace")
        parsed = _usda._parse_plantlst_text(text)
        _usda._PLANTLST_INDEX_CACHE[cache_key] = parsed
        return parsed
    except Exception:
        return None


def _generate_algorithmic_generic_transfers(binom_lower: str) -> List[str]:
    """
    Generate species-agnostic ICN Art. 23.5 generic transfer combinations across
    `Packera` (feminine), `Senecio` (masculine), and `Cineraria` (feminine) without
    any hardcoded taxon dictionary.
    """
    parts = binom_lower.split(" ", 1)
    if len(parts) != 2:
        return []
    _gen, ep = parts
    ep_variants: List[str] = [ep]
    if ep.endswith("ia") and len(ep) > 4:
        ep_variants.append(ep[:-2] + "ius")
    elif ep.endswith("a") and len(ep) > 3:
        ep_variants.append(ep[:-1] + "us")
    elif ep.endswith("us") and len(ep) > 4:
        ep_variants.append(ep[:-2] + "a")
    elif ep.endswith("is") and len(ep) > 4:
        ep_variants.append(ep[:-2] + "e")
    elif ep.endswith("e") and len(ep) > 4:
        ep_variants.append(ep[:-1] + "is")
    out: List[str] = []
    for g_cand in ("packera", "senecio", "cineraria"):
        for ep_v in ep_variants:
            comb = f"{g_cand} {ep_v}"
            if comb != binom_lower and comb not in out:
                out.append(comb)
    return out


def resolve_usda_symbols_for_taxon(
    canonical_packera: str,
    synonym_names: Optional[Iterable[str]] = None,
    cache_dir: Optional[Path] = None,
    db_path: Optional[Union[str, Path]] = None,
    allow_network: bool = True,
    logger: Optional[logging.Logger] = None,
) -> List[str]:
    """
    Dynamically resolve one or more USDA PLANTS accepted `Symbol` codes for
    `canonical_packera` without hardcoded taxon dictionaries, combining:
      1. Caller-supplied `synonym_names`, in-memory `_DYNAMIC_SYNONYM_TO_CANONICAL`,
         and persisted SQLite `taxon_synonyms` mappings.
      2. Algorithmic ICN Art. 23.5 generic transfers across `Packera`/`Senecio`/`Cineraria`.
      3. Dynamic GBIF Species Backbone synonym expansion when checklist symbols are not yet
         matched (e.g., resolving post-2006 combinations through their older checklist synonyms).
    """
    import scripts.core.usda_plants_ranges as _usda

    candidates: List[str] = []
    seen_cand: Set[str] = set()
    symbol_candidates: List[str] = []

    def _add_cand(raw: Any) -> None:
        txt_r = str(raw or "").strip()
        if re.fullmatch(r"[A-Za-z]{2,8}\d*", txt_r) and " " not in txt_r:
            sym_up = txt_r.upper()
            if sym_up not in symbol_candidates:
                symbol_candidates.append(sym_up)
        b, _ = _usda._extract_clean_binomial(raw)
        if b and b not in seen_cand:
            seen_cand.add(b)
            candidates.append(b)

    _add_cand(canonical_packera)
    for syn in synonym_names or ():
        _add_cand(syn)

    canon_clean_low = str(canonical_packera or "").strip().lower()
    syn_pairs: List[Tuple[str, str]] = []
    try:
        from scripts.core import harvester as _harv_mod
        syn_pairs.extend(
            (str(k).strip().lower(), str(v).strip().lower())
            for k, v in list(_harv_mod._DYNAMIC_SYNONYM_TO_CANONICAL.items())
        )
    except Exception:
        _harv_mod = None  # type: ignore[assignment]

    if db_path is not None and str(db_path) != ":memory:" and Path(db_path).is_file():
        try:
            conn_syn = _usda.get_db_connection(db_path=db_path, read_only=True)
            try:
                syn_pairs.extend(
                    (str(k).strip().lower(), str(v).strip().lower())
                    for k, v in _usda.load_taxon_synonyms(conn_syn).items()
                )
            finally:
                conn_syn.close()
        except Exception:
            pass

    target_canon_set: Set[str] = {canon_clean_low} if canon_clean_low else set()
    for raw_low, mapped_canon_low in syn_pairs:
        if raw_low == canon_clean_low and mapped_canon_low:
            target_canon_set.add(mapped_canon_low)
    for raw_low, mapped_canon_low in syn_pairs:
        if mapped_canon_low in target_canon_set or raw_low in target_canon_set:
            _add_cand(raw_low)
            if mapped_canon_low:
                _add_cand(mapped_canon_low)

    for base_c in list(candidates):
        for gen_comb in _generate_algorithmic_generic_transfers(base_c):
            _add_cand(gen_comb)

    if not candidates and not symbol_candidates:
        return []

    eff_dir = Path(cache_dir) if cache_dir is not None else _usda.DEFAULT_USDA_CACHE_DIR
    symbols: List[str] = []

    def _match_against_plantlst(idx_obj: Dict[str, Any]) -> List[str]:
        matched: List[str] = []
        sp_map: Dict[str, List[str]] = idx_obj["species_map"]
        infra_map: Dict[str, List[str]] = idx_obj["infra_map"]
        syn_sym_map: Dict[str, List[str]] = idx_obj.get("syn_symbol_map", {})
        known_syms: Set[str] = idx_obj.get("all_symbols", set())
        for cand in candidates:
            for sym in sp_map.get(cand, ()):
                if sym not in matched:
                    matched.append(sym)
        if not matched:
            for cand in candidates:
                for sym in infra_map.get(cand, ()):
                    if sym not in matched:
                        matched.append(sym)
        for s_cand in symbol_candidates:
            for acc_sym in syn_sym_map.get(s_cand, ()):
                if acc_sym not in matched:
                    matched.append(acc_sym)
            if s_cand in known_syms and s_cand not in matched:
                matched.append(s_cand)
        return matched

    # 1. Check local / downloaded plantlst.txt index first (O(1), zero per-species API latency)
    idx = _usda.load_or_fetch_plantlst_index(cache_dir=eff_dir, allow_network=allow_network, logger=logger)
    if idx is not None:
        symbols = _match_against_plantlst(idx)
        if symbols:
            return symbols

    # 1b. Dynamic GBIF Backbone synonym query if direct + generic-transfer candidates did not match
    if _harv_mod is not None and (allow_network or _harv_mod._is_gbif_species_mocked()):
        pre_n = len(candidates)
        gbif_acc, gbif_syns, gbif_auths = _harv_mod._query_gbif_backbone_and_synonyms(canonical_packera)
        for gs in (gbif_acc, *(gbif_syns or []), *(gbif_auths or [])):
            _add_cand(gs)
        for base_c in list(candidates[pre_n:]):
            for gen_comb in _generate_algorithmic_generic_transfers(base_c):
                _add_cand(gen_comb)
        if idx is not None and len(candidates) > pre_n:
            symbols = _match_against_plantlst(idx)
            if symbols:
                return symbols

    # 2. Fallback: POST /api/PlantSearch per candidate binomial
    if allow_network and not _usda._USDA_NETWORK_UNAVAILABLE:
        for cand in candidates:
            title_cand = f"{cand.split(' ', 1)[0].capitalize()} {cand.split(' ', 1)[1]}"
            resp = _usda._http_post_json(
                _usda.USDA_PLANT_SEARCH_API_URL,
                {"SearchCriterion": title_cand, "SearchType": 1},
            )
            if not resp:
                continue
            results = resp.get("PlantSearchResults") or []
            if not isinstance(results, list):
                continue
            for item in results:
                if not isinstance(item, dict):
                    continue
                sym = str(item.get("AcceptedSymbol") or item.get("Symbol") or "").strip().upper()
                sci = str(item.get("ScientificName") or "").strip()
                rank = str(item.get("Rank") or "").strip().lower()
                item_binom, is_infra = _usda._extract_clean_binomial(sci)
                if sym and item_binom in seen_cand and not is_infra and rank in {"", "species"}:
                    if sym not in symbols:
                        symbols.append(sym)
            if symbols:
                break

    return symbols
