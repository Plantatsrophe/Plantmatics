"""
Module: step01_usda_natureserve_utils.py
Description:
    NatureServe subnational state/province distribution lookup and voucher evaluation
    helpers for USDA PLANTS range fallback and Canadian occurrences.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import sqlite3
from typing import Any, Dict, Mapping, Optional, Sequence, Set, Tuple, Union

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _resolve_voucher_country_and_subnation(
    row: Mapping[str, Any],
    st_raw: Optional[str],
    st_norm: Optional[str],
) -> Tuple[Optional[str], Optional[str]]:
    """Normalize voucher `(country_code, subnation_code)` for NatureServe fallback checks."""
    import scripts.core.usda_plants_ranges as _usda

    raw_c = ""
    for col in ("countryCode", "country"):
        val = row.get(col)
        if val is not None and not pd.isna(val) and str(val).strip():
            raw_c = str(val).strip().upper()
            break
    c_code: Optional[str] = None
    if raw_c in {"US", "USA", "U.S.", "U.S.A.", "UNITED STATES", "UNITED STATES OF AMERICA"}:
        c_code = "US"
    elif raw_c in {"CA", "CAN", "CANADA"}:
        c_code = "CA"

    raw_st_clean = str(st_raw or "").strip()
    if re.match(r"^(?:US|CA)-([A-Za-z]{2})$", raw_st_clean, flags=re.IGNORECASE):
        c_prefix, sub_tok = raw_st_clean.split("-", 1)
        return (c_code or c_prefix.upper()), sub_tok.upper()

    if c_code == "CA" and raw_st_clean.upper() in {"CA", "CAN", "CANADA"}:
        return "CA", None

    st_norm_dwc = _usda._normalize_us_state_name(raw_st_clean)
    eff_full = st_norm_dwc or st_norm
    if not eff_full:
        return c_code, None
    if eff_full in _usda._CA_FULL_TO_ABBR or _usda._is_canadian_province_label(raw_st_clean):
        ca_full = _usda._CANADIAN_PROV_ABBR_TO_FULL.get(eff_full, eff_full)
        sub_ca = _usda._CA_FULL_TO_ABBR.get(ca_full)
        if sub_ca and (c_code is None or c_code == "CA"):
            return "CA", sub_ca
    if eff_full in _usda._US_FULL_TO_ABBR:
        sub_us = _usda._US_FULL_TO_ABBR[eff_full]
        if c_code is None or c_code == "US":
            return "US", sub_us
    return c_code, None


def _load_natureserve_ranges_for_taxa(
    taxa: Sequence[str],
    *,
    conn: Optional[sqlite3.Connection] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Path = PROJECT_ROOT,
    cache_dir: Optional[Path] = None,
) -> Dict[str, Dict[Tuple[str, str], Dict[str, Any]]]:
    """Read subnational state/province distributions from SQLite `natureserve_subnational_ranges`."""
    import scripts.core.usda_plants_ranges as _usda

    clean_taxa = sorted({str(t).strip() for t in taxa if str(t or "").strip()})
    if not clean_taxa:
        return {}

    own_conn = False
    active_conn = conn
    if active_conn is None:
        eff_cache_dir = Path(cache_dir) if cache_dir is not None else (Path(workspace_dir) / "data" / "environmental" / "usda_plants")
        use_sqlite = db_path is not None or (
            os.environ.get("PLANTMATICS_DEMO_MODE") != "1"
            and os.environ.get("PACKERA_DEMO_MODE") != "1"
            and _usda.is_canonical_or_explicit_db_context(eff_cache_dir, workspace_root=workspace_dir)
        )
        if not use_sqlite:
            return {}
        try:
            resolved_db = _usda.resolve_db_path(db_path=db_path, workspace_root=workspace_dir)
        except Exception:
            return {}
        if resolved_db == ":memory:" or not Path(resolved_db).is_file():
            return {}
        try:
            active_conn = _usda.get_db_connection(db_path=resolved_db, read_only=True, workspace_root=workspace_dir)
            own_conn = True
        except Exception:
            return {}

    out: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]] = {}
    try:
        tbl_chk = active_conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'natureserve_subnational_ranges';"
        ).fetchone()
        if not tbl_chk:
            return {}
        placeholders = ", ".join("?" for _ in clean_taxa)
        rows = active_conn.execute(
            f"""
            SELECT species_standardized, country_code, subnation_code, s_rank,
                   is_native, is_exotic, is_historical_or_extirpated
            FROM natureserve_subnational_ranges
            WHERE species_standardized IN ({placeholders});
            """,
            tuple(clean_taxa),
        ).fetchall()
        for r in rows:
            sp = str(r[0] or "").strip()
            c_code = str(r[1] or "").strip().upper()
            sub_code = str(r[2] or "").strip().upper()
            if not sp or not c_code or not sub_code:
                continue
            out.setdefault(sp, {})[(c_code, sub_code)] = {
                "s_rank": str(r[3] or "").strip(),
                "is_native": bool(int(r[4])) if r[4] is not None else True,
                "is_exotic": bool(int(r[5])) if r[5] is not None else False,
                "is_historical_or_extirpated": bool(int(r[6])) if r[6] is not None else False,
            }
    except sqlite3.Error:
        return {}
    finally:
        if own_conn and active_conn is not None:
            active_conn.close()
    return out


def _evaluate_natureserve_voucher(
    c_code: Optional[str],
    sub_code: Optional[str],
    ns_map: Mapping[Tuple[str, str], Mapping[str, Any]],
    extra_us_states: Optional[Set[str]] = None,
) -> str:
    """Evaluate a US/CA voucher against NatureServe subnational records + 1-ring state adjacency."""
    import scripts.core.usda_plants_ranges as _usda

    if c_code not in {"US", "CA"} or not sub_code or not ns_map:
        return "non_us_or_unresolved"
    rec = ns_map.get((c_code, sub_code))
    if rec is not None:
        if rec.get("is_exotic") and not rec.get("is_native"):
            return "natureserve_exotic_or_waif"
        if rec.get("is_exotic"):
            return "natureserve_exotic_or_waif"
        return "natureserve_state_in_range"

    mapped_subs = {sub for (_, sub) in ns_map.keys()}
    for st_full in extra_us_states or ():
        abbr = _usda._US_FULL_TO_ABBR.get(st_full)
        if abbr:
            mapped_subs.add(abbr)
    if _usda._SUBNATION_ADJACENCY.get(sub_code, set()) & mapped_subs:
        return "buffered_adjacent"
    return "out_of_range"
