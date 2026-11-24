"""
Module: step01_usda_distribution_io.py
Description:
    SQLite (`usda_county_ranges`) persistence, CSV cache management, and USDA PLANTS
    distribution endpoint retrieval for canonical taxa.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import re
import sqlite3
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union
import urllib.parse

from scripts.core.step01_usda_csv_parser import (  # noqa: F401
    _convert_usda_json_records_to_csv,
    _normalize_usda_distribution_payload,
    _slice_csv_from_header,
    parse_usda_distribution_csv,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _load_taxon_range_from_sqlite(
    conn: sqlite3.Connection,
    canonical_packera: str,
) -> Optional[Dict[str, Any]]:
    """Load persisted USDA county distribution for `canonical_packera` from SQLite if present."""
    import scripts.core.usda_plants_ranges as _usda

    try:
        cur = conn.execute(
            """
            SELECT usda_symbol, state_norm, county_norm, fips_code, native_status, is_state_only
            FROM usda_county_ranges
            WHERE canonical_packera = ?;
            """,
            (canonical_packera,),
        )
        fetched = cur.fetchall()
    except sqlite3.Error:
        return None

    if not fetched:
        return None

    counties: Set[str] = set()
    state_only_states: Set[str] = set()
    documented_states: Set[str] = set()
    symbols: List[str] = []

    for r in fetched:
        sym = str(r[0] or "").strip().upper()
        st_norm = str(r[1] or "").strip().lower()
        co_norm = str(r[2] or "").strip().lower()
        is_st_only = bool(int(r[5] or 0))
        if sym and sym not in symbols:
            symbols.append(sym)
        if not st_norm:
            continue
        documented_states.add(st_norm)
        if is_st_only or not co_norm:
            state_only_states.add(st_norm)
        else:
            counties.add(f"{st_norm}|{co_norm}")

    state_only_states = {
        st for st in state_only_states
        if not any(ck.startswith(f"{st}|") for ck in counties)
    }
    counties = _usda._augment_documented_counties(canonical_packera, counties, documented_states)
    return {
        "canonical_packera": canonical_packera,
        "symbols": symbols,
        "counties": counties,
        "state_only_states": state_only_states,
        "documented_states": documented_states,
    }


def _persist_taxon_range_to_sqlite(
    conn: sqlite3.Connection,
    canonical_packera: str,
    rows: Sequence[Mapping[str, Any]],
) -> int:
    """Persist normalized USDA county presence rows for `canonical_packera` into `usda_county_ranges`."""
    if not rows:
        return 0
    now_iso = datetime.now(timezone.utc).isoformat()
    payload: List[Tuple[str, str, str, str, Optional[str], Optional[str], int, str]] = []
    for r in rows:
        st_norm = str(r.get("state_norm") or "").strip().lower()
        co_norm = str(r.get("county_norm") or "").strip().lower()
        sym = str(r.get("usda_symbol") or "").strip().upper()
        if not st_norm or not sym:
            continue
        payload.append((
            canonical_packera,
            sym,
            st_norm,
            co_norm,
            str(r.get("fips_code") or "").strip() or None,
            str(r.get("native_status") or "").strip() or None,
            1 if r.get("is_state_only") else 0,
            now_iso,
        ))
    if not payload:
        return 0

    sql = """
    INSERT INTO usda_county_ranges (
        canonical_packera, usda_symbol, state_norm, county_norm,
        fips_code, native_status, is_state_only, updated_at
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(canonical_packera, state_norm, county_norm) DO UPDATE SET
        usda_symbol = excluded.usda_symbol,
        fips_code = COALESCE(excluded.fips_code, usda_county_ranges.fips_code),
        native_status = COALESCE(excluded.native_status, usda_county_ranges.native_status),
        is_state_only = excluded.is_state_only,
        updated_at = excluded.updated_at;
    """
    with conn:
        conn.executemany(sql, payload)
    return len(payload)


def _fetch_usda_symbol_distribution_csv(sym: str, timeout: float = 15.0) -> Optional[str]:
    """
    Download county distribution CSV for USDA symbol `sym` using:
      1. `GET /api/PlantProfile?symbol=<SYM>` (extracts numeric `Id`/`MasterId` or inline records,
         and short-circuits if mocked `_http_get_text` returns CSV directly).
      2. `POST /api/PlantProfile/getDownloadDistributionDocumentation` with `MasterId` + `Symbol`.
      3. `GET /api/PlantProfileDataDownload?symbol=<SYM>` legacy fallback.
    """
    import scripts.core.usda_plants_ranges as _usda

    master_id: Optional[int] = None
    prof_url = f"{_usda.USDA_PLANT_PROFILE_API_URL}?{urllib.parse.urlencode({'symbol': sym})}"
    prof_raw = _usda._http_get_text(prof_url, timeout=timeout)
    if prof_raw:
        direct_norm = _usda._normalize_usda_distribution_payload(prof_raw, default_symbol=sym)
        if direct_norm:
            return direct_norm
        try:
            prof_obj = json.loads(prof_raw)
            if isinstance(prof_obj, dict):
                raw_id = (
                    prof_obj.get("Id")
                    or prof_obj.get("MasterId")
                    or prof_obj.get("PlantId")
                    or prof_obj.get("id")
                    or prof_obj.get("masterId")
                )
                if raw_id is not None and str(raw_id).strip().isdigit():
                    master_id = int(str(raw_id).strip())
        except Exception:
            pass

    if master_id is None and not _usda._USDA_NETWORK_UNAVAILABLE:
        srch = _usda._http_post_json(_usda.USDA_PLANT_SEARCH_API_URL, {"SearchCriterion": sym, "SearchType": 0}, timeout=10.0)
        if srch and isinstance(srch.get("PlantSearchResults"), list):
            for item in srch["PlantSearchResults"]:
                if isinstance(item, dict):
                    item_sym = str(item.get("Symbol") or item.get("AcceptedSymbol") or "").strip().upper()
                    raw_id = item.get("Id") or item.get("MasterId") or item.get("PlantId")
                    if item_sym == sym and raw_id is not None and str(raw_id).strip().isdigit():
                        master_id = int(str(raw_id).strip())
                        break

    if not _usda._USDA_NETWORK_UNAVAILABLE:
        post_payload: Dict[str, Any] = {
            "Field": "Symbol",
            "SortBy": "sortSciName",
            "Offset": None,
            "MasterId": master_id,
            "Symbol": sym,
            "Value": sym,
        }
        docs_raw = _usda._http_post_text(_usda.USDA_DISTRIBUTION_DOCS_API_URL, payload=post_payload, timeout=timeout)
        docs_norm = _usda._normalize_usda_distribution_payload(docs_raw, default_symbol=sym)
        if docs_norm:
            return docs_norm

    if not _usda._USDA_NETWORK_UNAVAILABLE:
        dl_url = f"{_usda.USDA_DISTRIBUTION_CSV_API_URL}?{urllib.parse.urlencode({'symbol': sym})}"
        legacy_raw = _usda._http_get_text(dl_url, timeout=timeout)
        legacy_norm = _usda._normalize_usda_distribution_payload(legacy_raw, default_symbol=sym)
        if legacy_norm:
            return legacy_norm

    return None


def fetch_or_load_usda_taxon_range(
    canonical_packera: str,
    synonym_names: Optional[Iterable[str]] = None,
    cache_dir: Optional[Path] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Path = PROJECT_ROOT,
    allow_network: bool = True,
    force_refresh: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """
    Load or download the unioned USDA PLANTS county range for `canonical_packera`
    across all mapped USDA symbols and synonyms, persisting to CSV cache and SQLite.
    """
    import scripts.core.usda_plants_ranges as _usda

    eff_cache_dir = Path(cache_dir) if cache_dir is not None else (Path(workspace_dir) / "data" / "environmental" / "usda_plants")
    dist_dir = eff_cache_dir / "distributions"

    # 1. Check SQLite cache first unless force_refresh=True
    use_sqlite = db_path is not None or (
        os.environ.get("PLANTMATICS_DEMO_MODE") != "1"
        and os.environ.get("PACKERA_DEMO_MODE") != "1"
        and _usda.is_canonical_or_explicit_db_context(eff_cache_dir, workspace_root=workspace_dir)
    )
    resolved_db: Optional[Union[str, Path]] = None
    if use_sqlite:
        try:
            resolved_db = _usda.resolve_db_path(db_path=db_path, workspace_root=workspace_dir)
        except Exception:
            resolved_db = None

    if not force_refresh and resolved_db is not None and resolved_db != ":memory:" and Path(resolved_db).is_file():
        try:
            conn = _usda.get_db_connection(db_path=resolved_db, read_only=True, workspace_root=workspace_dir)
            try:
                cached_sql = _usda._load_taxon_range_from_sqlite(conn, canonical_packera)
                if cached_sql is not None and (cached_sql["counties"] or cached_sql["state_only_states"]):
                    return cached_sql
            finally:
                conn.close()
        except Exception:
            pass

    # 2. Check taxon-named local CSV override/cache `<Packera_epithet>.csv` if present
    slug_csv = dist_dir / f"{re.sub(r'[^A-Za-z0-9_-]+', '_', canonical_packera.strip())}.csv"
    if not force_refresh and slug_csv.is_file() and slug_csv.stat().st_size > 0:
        parsed_slug = _usda.parse_usda_distribution_csv(
            slug_csv.read_text(encoding="utf-8-sig", errors="replace"),
            default_symbol=canonical_packera.upper(),
        )
        syms = sorted({str(r["usda_symbol"]) for r in parsed_slug["rows"] if r.get("usda_symbol")})
        aug_slug_counties = _usda._augment_documented_counties(
            canonical_packera,
            set(parsed_slug["counties"]),
            set(parsed_slug["documented_states"]),
        )
        result_slug = {
            "canonical_packera": canonical_packera,
            "symbols": syms,
            "counties": aug_slug_counties,
            "state_only_states": parsed_slug["state_only_states"],
            "documented_states": parsed_slug["documented_states"],
        }
        if resolved_db is not None and parsed_slug["rows"]:
            try:
                conn_w = _usda.get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=workspace_dir)
                try:
                    _usda._persist_taxon_range_to_sqlite(conn_w, canonical_packera, parsed_slug["rows"])
                finally:
                    conn_w.close()
            except Exception:
                pass
        return result_slug

    # 3. Resolve USDA symbol(s) via checklist / PlantSearch API
    symbols = _usda.resolve_usda_symbols_for_taxon(
        canonical_packera=canonical_packera,
        synonym_names=synonym_names,
        cache_dir=eff_cache_dir,
        db_path=resolved_db,
        allow_network=allow_network,
        logger=logger,
    )
    if not symbols:
        return {
            "canonical_packera": canonical_packera,
            "symbols": [],
            "counties": set(),
            "state_only_states": set(),
            "documented_states": set(),
        }

    all_counties: Set[str] = set()
    all_state_only: Set[str] = set()
    all_doc_states: Set[str] = set()
    all_rows: List[Dict[str, Any]] = []

    for sym in symbols:
        sym_csv_path = dist_dir / f"{sym}.csv"
        csv_text: Optional[str] = None
        if not force_refresh and sym_csv_path.is_file() and sym_csv_path.stat().st_size > 0:
            try:
                csv_text = sym_csv_path.read_text(encoding="utf-8-sig", errors="replace")
            except Exception:
                csv_text = None

        if csv_text is not None:
            probe_cached = _usda.parse_usda_distribution_csv(csv_text, default_symbol=sym)
            if not probe_cached["counties"] and not probe_cached["state_only_states"] and allow_network:
                csv_text = None

        if csv_text is None and allow_network and not _usda._USDA_NETWORK_UNAVAILABLE:
            csv_text = _usda._fetch_usda_symbol_distribution_csv(sym, timeout=15.0)
            if csv_text:
                try:
                    _usda.atomic_write_text(sym_csv_path, csv_text)
                except Exception as exc:
                    if logger:
                        logger.debug(f"Could not write cached USDA CSV {sym_csv_path}: {exc}")

        if not csv_text:
            continue

        parsed = _usda.parse_usda_distribution_csv(csv_text, default_symbol=sym)
        all_counties.update(parsed["counties"])
        all_state_only.update(parsed["state_only_states"])
        all_doc_states.update(parsed["documented_states"])
        all_rows.extend(parsed["rows"])

    # Any state with >=1 county record across synonyms is county-resolved rather than state-only
    all_state_only = {
        st for st in all_state_only
        if not any(ck.startswith(f"{st}|") for ck in all_counties)
    }
    filtered_rows = [
        r for r in all_rows
        if not (r.get("is_state_only") and r.get("state_norm") not in all_state_only)
    ]

    if resolved_db is not None and filtered_rows:
        try:
            conn_w = _usda.get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=workspace_dir)
            try:
                _usda._persist_taxon_range_to_sqlite(conn_w, canonical_packera, filtered_rows)
            finally:
                conn_w.close()
        except Exception as db_exc:
            if logger:
                logger.debug(f"SQLite usda_county_ranges persistence skipped: {db_exc}")

    all_counties = _usda._augment_documented_counties(canonical_packera, all_counties, all_doc_states)

    return {
        "canonical_packera": canonical_packera,
        "symbols": symbols,
        "counties": all_counties,
        "state_only_states": all_state_only,
        "documented_states": all_doc_states,
    }
