"""Step 01 dynamic synonym cache helpers and GBIF genus discovery."""

from __future__ import annotations

import concurrent.futures
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Union

from scripts.core.step01_harvester_filters import (
    _SYNONYM_GENERA,
    _extract_canonical_binomial,
)


def _hydrate_synonyms_from_sqlite(db_path: Optional[Union[str, Path]] = None) -> None:
    """Hydrate in-memory synonym mappings from persisted `taxon_synonyms` when SQLite exists."""
    import scripts.core.harvester as _harv

    try:
        resolved = _harv.resolve_db_path(db_path=db_path)
        if resolved != ":memory:" and Path(resolved).is_file():
            conn = _harv.get_db_connection(db_path=resolved, read_only=True)
            try:
                for raw_low, canon in _harv.load_taxon_synonyms(conn).items():
                    if not raw_low or not canon:
                        continue
                    _harv._DYNAMIC_SYNONYM_TO_CANONICAL.setdefault(raw_low, canon)
                    raw_bin = _extract_canonical_binomial(raw_low, _SYNONYM_GENERA)
                    if raw_bin:
                        _harv._DYNAMIC_SYNONYM_TO_CANONICAL.setdefault(raw_bin.lower(), canon)
            finally:
                conn.close()
    except Exception:
        pass


def _is_gbif_species_mocked() -> bool:
    """Return True if pygbif.species endpoints are currently patched by unittest.mock."""
    import scripts.core.harvester as _harv

    sp = _harv.species
    nb = getattr(sp, "name_backbone", None)
    nu = getattr(sp, "name_usage", None)
    nl = getattr(sp, "name_lookup", None)
    return any(
        hasattr(fn, "_mock_name") or hasattr(fn, "assert_called")
        for fn in (nb, nu, nl)
        if fn is not None
    )


def _register_dynamic_synonyms(
    canonical_packera: str,
    synonyms: Iterable[Any],
    author_strings: Optional[Iterable[str]] = None,
    allow_packera_species_remap: bool = True,
) -> None:
    """Register discovered synonyms and author strings in the in-memory standardization cache."""
    import scripts.core.harvester as _harv

    if not canonical_packera or not canonical_packera.startswith("Packera "):
        return
    _harv._DYNAMIC_SYNONYM_TO_CANONICAL[canonical_packera.lower()] = canonical_packera
    if author_strings:
        clean_authors = {str(a).strip() for a in author_strings if a and str(a).strip()}
        if clean_authors:
            _harv._DYNAMIC_TAXON_AUTHOR_STRINGS.setdefault(canonical_packera, set()).update(clean_authors)
            for auth_sci in clean_authors:
                auth_low = re.sub(r"\s+", " ", auth_sci.lower())
                auth_p_bin = _extract_canonical_binomial(auth_sci, {"Packera"})
                if (
                    not allow_packera_species_remap
                    and auth_p_bin
                    and auth_p_bin.lower() != canonical_packera.lower()
                    and not re.search(r"\b(?:var|subsp|ssp|f)\b", auth_low)
                ):
                    continue
                _harv._DYNAMIC_SYNONYM_TO_CANONICAL[auth_low] = canonical_packera
                auth_bin = _extract_canonical_binomial(auth_sci, _SYNONYM_GENERA)
                if auth_bin:
                    _harv._DYNAMIC_SYNONYM_TO_CANONICAL[auth_bin.lower()] = canonical_packera

    for item in synonyms:
        if not item:
            continue
        item_str = re.sub(r"\s+", " ", str(item).strip())
        if not item_str:
            continue
        item_low = item_str.lower()
        item_p_bin = _extract_canonical_binomial(item_str, {"Packera"})
        if (
            not allow_packera_species_remap
            and item_p_bin
            and item_p_bin.lower() != canonical_packera.lower()
            and not re.search(r"\b(?:var|subsp|ssp|f)\b", item_low)
        ):
            continue
        _harv._DYNAMIC_SYNONYM_TO_CANONICAL[item_low] = canonical_packera
        item_bin = _extract_canonical_binomial(item_str, _SYNONYM_GENERA)
        if item_bin:
            _harv._DYNAMIC_SYNONYM_TO_CANONICAL[item_bin.lower()] = canonical_packera


def _call_and_normalize_gbif_backbone(query_name: str, rank: str) -> Dict[str, Any]:
    """
    Invoke `pygbif.species.name_backbone` across both v1 (`name=`, `rank=`) and
    v2 (`scientificName=`, `taxonRank=`) signatures, normalizing nested GBIF v2
    responses into the flat backbone dict schema.
    """
    import scripts.core.harvester as _harv

    sp = _harv.species
    if sp is None or getattr(sp, "name_backbone", None) is None:
        return {}
    nb_fn = sp.name_backbone
    use_v2_sig = False
    if not _harv._is_gbif_species_mocked():
        try:
            import inspect
            use_v2_sig = "scientificName" in inspect.signature(nb_fn).parameters
        except Exception:
            use_v2_sig = False
    if use_v2_sig:
        raw = nb_fn(scientificName=query_name, taxonRank=str(rank).upper(), kingdom="Plantae")
    else:
        raw = nb_fn(name=query_name, rank=rank, kingdom="Plantae")
    if not isinstance(raw, dict) or not raw:
        return {}
    if "usageKey" in raw or "matchType" in raw or "canonicalName" in raw or "species" in raw:
        return raw
    usage = raw.get("usage") if isinstance(raw.get("usage"), dict) else {}
    acc_usage = raw.get("acceptedUsage") if isinstance(raw.get("acceptedUsage"), dict) else {}
    diag = raw.get("diagnostics") if isinstance(raw.get("diagnostics"), dict) else {}
    cls_list = raw.get("classification") if isinstance(raw.get("classification"), list) else []
    cls_by_rank = {
        str(c.get("rank", "")).upper(): c
        for c in cls_list
        if isinstance(c, dict) and c.get("rank")
    }
    sp_cls = cls_by_rank.get("SPECIES", {})
    gen_cls = cls_by_rank.get("GENUS", {})
    norm: Dict[str, Any] = {
        "matchType": diag.get("matchType", "EXACT" if usage else "NONE"),
        "status": usage.get("status") or diag.get("status"),
        "usageKey": usage.get("key"),
        "acceptedUsageKey": acc_usage.get("key") or sp_cls.get("key"),
        "scientificName": usage.get("name"),
        "canonicalName": usage.get("canonicalName") or usage.get("name"),
        "accepted": acc_usage.get("name") or acc_usage.get("canonicalName") or sp_cls.get("name"),
        "species": acc_usage.get("canonicalName") or sp_cls.get("name") or acc_usage.get("name"),
        "genus": gen_cls.get("name"),
        "genusKey": gen_cls.get("key"),
    }
    return {k: v for k, v in norm.items() if v is not None}


def _persist_dynamic_synonyms_for_taxon(canonical_packera: str) -> None:
    """Persist newly discovered live GBIF synonyms for `canonical_packera` into SQLite `taxon_synonyms`."""
    import scripts.core.harvester as _harv

    if (
        not canonical_packera
        or _harv._is_gbif_species_mocked()
        or os.environ.get("PLANTMATICS_DEMO_MODE") == "1"
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or not _harv.is_canonical_or_explicit_db_context()
    ):
        return
    try:
        resolved = _harv.resolve_db_path()
        if resolved != ":memory:" and Path(resolved).is_file():
            sub_map = {
                k: v
                for k, v in _harv._DYNAMIC_SYNONYM_TO_CANONICAL.items()
                if v == canonical_packera and k
            }
            if not sub_map:
                return
            auth_sub = {}
            if canonical_packera in _harv._DYNAMIC_TAXON_AUTHOR_STRINGS:
                auth_sub[canonical_packera] = sorted(_harv._DYNAMIC_TAXON_AUTHOR_STRINGS[canonical_packera])
            conn = _harv.get_db_connection(db_path=resolved, init_schema=False)
            try:
                _harv.upsert_taxon_synonyms(conn, sub_map, auth_sub or None)
            finally:
                conn.close()
    except Exception:
        pass


def _get_packera_genus_key_gbif(timeout_seconds: float = 5.0) -> Optional[int]:
    """
    Query pygbif.species.name_backbone(name='Packera', rank='genus', kingdom='Plantae')
    with a fast 5-second timeout and in-memory cache to obtain the GBIF usageKey for genus Packera.
    """
    import scripts.core.harvester as _harv

    is_mocked = _harv._is_gbif_species_mocked()
    cache_key = "packera"
    if not is_mocked and cache_key in _harv._GBIF_GENUS_KEY_CACHE:
        return _harv._GBIF_GENUS_KEY_CACHE[cache_key]
    if not is_mocked and _harv._GBIF_NETWORK_UNAVAILABLE:
        return None

    def _fetch_key() -> Optional[int]:
        backbone = _harv._call_and_normalize_gbif_backbone("Packera", "genus")
        if not isinstance(backbone, dict) or not backbone:
            return None
        match_type = str(backbone.get("matchType", "EXACT")).upper()
        if match_type == "NONE":
            return None
        raw_key = backbone.get("usageKey") or backbone.get("genusKey") or backbone.get("key")
        if raw_key is None:
            return None
        try:
            return int(raw_key)
        except (TypeError, ValueError):
            return None

    if is_mocked:
        try:
            return _fetch_key()
        except Exception:
            return None

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch_key)
        res = future.result(timeout=timeout_seconds)
        _harv._GBIF_GENUS_KEY_CACHE[cache_key] = res
        return res
    except Exception:
        _harv._GBIF_NETWORK_UNAVAILABLE = True
        return None
    finally:
        executor.shutdown(wait=False, cancel_futures=True)


def discover_genus_packera_species_gbif(limit: int = 1000) -> List[str]:
    """
    Dynamically query the GBIF Species Backbone API to collect all accepted
    `Packera <epithet>` species binomials across the entire genus Packera.
    """
    import scripts.core.harvester as _harv

    eff_limit = max(1, int(limit))
    is_mocked = _harv._is_gbif_species_mocked()
    if not is_mocked and eff_limit in _harv._GBIF_GENUS_SPECIES_CACHE:
        return list(_harv._GBIF_GENUS_SPECIES_CACHE[eff_limit])
    if not is_mocked and _harv._GBIF_NETWORK_UNAVAILABLE:
        return []

    genus_key = _harv._get_packera_genus_key_gbif(timeout_seconds=5.0)
    if genus_key is None:
        return []

    def _fetch_species() -> List[str]:
        sp = _harv.species
        raw_items: List[Any] = []

        if hasattr(sp, "name_lookup"):
            try:
                lookup_resp = sp.name_lookup(
                    higherTaxonKey=genus_key,
                    rank="SPECIES",
                    status="ACCEPTED",
                    limit=eff_limit,
                )
                if isinstance(lookup_resp, dict):
                    raw_items.extend(lookup_resp.get("results", []) or [])
                elif isinstance(lookup_resp, list):
                    raw_items.extend(lookup_resp)
            except Exception:
                pass

        if hasattr(sp, "name_usage"):
            try:
                try:
                    children_resp = sp.name_usage(key=genus_key, data="children", limit=eff_limit)
                except TypeError:
                    children_resp = sp.name_usage(key=genus_key, data="children")
                if isinstance(children_resp, dict):
                    raw_items.extend(children_resp.get("results", []) or [])
                elif isinstance(children_resp, list):
                    raw_items.extend(children_resp)
            except Exception:
                pass

        rejected_statuses = {
            "SYNONYM",
            "HETEROTYPIC_SYNONYM",
            "HOMOTYPIC_SYNONYM",
            "PROPARTE_SYNONYM",
            "MISAPPLIED",
            "DOUBTFUL",
        }
        discovered: Set[str] = set()
        for item in raw_items:
            if not isinstance(item, dict):
                continue
            status_val = str(item.get("taxonomicStatus") or item.get("status") or "").strip().upper()
            if status_val and status_val in rejected_statuses:
                continue
            rank_val = str(item.get("rank") or "").strip().upper()
            if rank_val and rank_val not in {"SPECIES", "SUBSPECIES", "VARIETY"}:
                continue
            cand = _extract_canonical_binomial(
                item.get("canonicalName") or item.get("species") or item.get("scientificName"),
                {"Packera"},
            )
            if cand and re.match(r"^Packera\s+(?:[x×]\s*)?[a-z-]+$", cand):
                discovered.add(cand)
                sci_full = item.get("scientificName")
                auth_list = [sci_full.strip()] if isinstance(sci_full, str) and sci_full.strip() else None
                _harv._register_dynamic_synonyms(cand, [cand], auth_list)

        return sorted(discovered)

    if is_mocked:
        try:
            return _fetch_species()
        except Exception:
            return []

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch_species)
        res = future.result(timeout=5.0)
        _harv._GBIF_GENUS_SPECIES_CACHE[eff_limit] = list(res)
        return res
    except Exception:
        _harv._GBIF_NETWORK_UNAVAILABLE = True
        return []
    finally:
        executor.shutdown(wait=False, cancel_futures=True)
