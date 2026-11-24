"""
Module: step01_usda_csv_parser.py
Description:
    Normalization and schema-agnostic CSV/JSON parsing for USDA NRCS PLANTS
    distribution exports.
"""

from __future__ import annotations

import csv
import io
import json
import re
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple


def _slice_csv_from_header(raw_text: str) -> str:
    """Strip any leading metadata/title lines prior to the actual CSV header row."""
    cleaned = raw_text.lstrip("\ufeff").strip()
    if not cleaned:
        return ""
    lines = cleaned.splitlines()
    for idx, ln in enumerate(lines[:12]):
        low = ln.lower()
        if "," in low and any(h in low for h in ("state", "county", "fips", "symbol")):
            return "\n".join(lines[idx:])
    return cleaned


def _convert_usda_json_records_to_csv(records: Sequence[Mapping[str, Any]], default_symbol: str = "") -> str:
    """Serialize USDA JSON distribution dicts into canonical 7-column CSV text."""
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["Symbol", "Scientific Name", "Common Name", "State", "County", "FIPS Code", "Distribution Status"])
    wrote = 0
    for rec in records:
        if not isinstance(rec, Mapping):
            continue
        norm = {re.sub(r"[\s_]+", " ", str(k or "").strip().lower()): v for k, v in rec.items()}
        sym = str(norm.get("symbol") or norm.get("acceptedsymbol") or default_symbol).strip().upper()
        sci = str(norm.get("scientific name") or norm.get("scientificname") or "").strip()
        com = str(norm.get("common name") or norm.get("commonname") or "").strip()
        st = str(
            norm.get("state")
            or norm.get("state and province")
            or norm.get("statename")
            or norm.get("state name")
            or norm.get("stateprovince")
            or norm.get("location")
            or ""
        ).strip()
        co = str(
            norm.get("county")
            or norm.get("county name")
            or norm.get("countyname")
            or ""
        ).strip()
        st_fip = str(norm.get("state fip") or norm.get("statefip") or "").strip()
        co_fip = str(norm.get("county fip") or norm.get("countyfip") or "").strip()
        fips = str(norm.get("fips code") or norm.get("fips") or norm.get("fipscode") or "").strip()
        if not fips and (st_fip or co_fip):
            fips = f"{st_fip.zfill(2)}{co_fip.zfill(3)}" if (st_fip and co_fip) else (co_fip or st_fip)
        status = str(
            norm.get("distribution status")
            or norm.get("distributionstatus")
            or norm.get("native status")
            or norm.get("nativestatus")
            or norm.get("status")
            or ""
        ).strip()
        if st:
            writer.writerow([sym, sci, com, st, co, fips, status])
            wrote += 1
    return buf.getvalue() if wrote > 0 else ""


def _normalize_usda_distribution_payload(raw_text: Optional[str], default_symbol: str = "") -> Optional[str]:
    """Accept raw CSV text, JSON-wrapped CSV text, or JSON distribution objects/arrays from USDA PLANTS."""
    if not raw_text:
        return None
    s = raw_text.lstrip("\ufeff").strip()
    if not s:
        return None

    if s[0] in ('{', '[', '"') and (s[0] != '"' or s.endswith('"')):
        try:
            obj = json.loads(s)
        except Exception:
            obj = None
        if isinstance(obj, str):
            return _normalize_usda_distribution_payload(obj, default_symbol=default_symbol)
        if isinstance(obj, list):
            csv_conv = _convert_usda_json_records_to_csv(obj, default_symbol=default_symbol)
            return csv_conv or None
        if isinstance(obj, dict):
            for str_key in ("Csv", "CsvData", "csv", "csvData", "FileContents", "Content", "Data", "data"):
                val = obj.get(str_key)
                if isinstance(val, str) and "," in val:
                    res = _normalize_usda_distribution_payload(val, default_symbol=default_symbol)
                    if res:
                        return res
            for list_key in (
                "PlantsDistributionResults",
                "plantsDistributionResults",
                "DistributionDocumentation",
                "Distributions",
                "distributions",
                "Counties",
                "counties",
                "Data",
                "data",
                "Results",
                "results",
            ):
                val = obj.get(list_key)
                if isinstance(val, list) and val:
                    csv_conv = _convert_usda_json_records_to_csv(val, default_symbol=default_symbol)
                    if csv_conv:
                        return csv_conv
            return None

    sliced = _slice_csv_from_header(s)
    parsed_probe = parse_usda_distribution_csv(sliced, default_symbol=default_symbol)
    if parsed_probe["counties"] or parsed_probe["state_only_states"]:
        return sliced
    return None


def parse_usda_distribution_csv(
    csv_text: str,
    default_symbol: str = "",
) -> Dict[str, Any]:
    """
    Parse USDA PLANTS distribution CSV across both USDA PLANTS 3.0 export schemas:
      - `Symbol,Scientific Name,Common Name,State and Province,County Name,FIPS Code,Distribution Status`
      - `Symbol,Synonym Symbol,Scientific Name,State FIP,State,County FIP,County,Native Status`

    Distinguishes county-resolved states from state-only (`No County Data`) states:
    a state is marked `state_only` ONLY when USDA records non-excluded presence for
    that state with blank `County` and zero county-level rows in that state.
    """
    import scripts.core.usda_plants_ranges as _usda

    counties: Set[str] = set()
    counties_by_state: Dict[str, Set[str]] = {}
    present_us_states: Set[str] = set()
    db_rows: List[Dict[str, Any]] = []

    if not csv_text or not csv_text.strip():
        return {
            "counties": counties,
            "state_only_states": set(),
            "documented_states": set(),
            "rows": db_rows,
        }

    sliced_csv = _slice_csv_from_header(csv_text)
    reader = csv.DictReader(io.StringIO(sliced_csv))
    raw_parsed_rows: List[Tuple[str, str, Optional[str], str, str]] = []

    for row in reader:
        if not isinstance(row, dict):
            continue
        norm_row = {
            re.sub(r"[\s_]+", " ", str(k or "").strip().lower()): str(v or "").strip()
            for k, v in row.items()
            if k is not None
        }
        sym = (norm_row.get("symbol") or default_symbol).strip().upper()
        raw_state = (
            norm_row.get("state")
            or norm_row.get("state and province")
            or norm_row.get("stateprovince")
            or norm_row.get("state name")
            or norm_row.get("statename")
            or ""
        )
        raw_county = (
            norm_row.get("county")
            or norm_row.get("county name")
            or norm_row.get("countyname")
            or ""
        )
        st_fip = norm_row.get("state fip") or norm_row.get("statefip") or ""
        co_fip = norm_row.get("county fip") or norm_row.get("countyfip") or ""
        fips = norm_row.get("fips code") or norm_row.get("fips") or norm_row.get("fipscode") or ""
        if not fips and (st_fip or co_fip):
            fips = f"{st_fip.zfill(2)}{co_fip.zfill(3)}" if (st_fip and co_fip) else (co_fip or st_fip)
        native_status = (
            norm_row.get("native status")
            or norm_row.get("distribution status")
            or norm_row.get("nativestatus")
            or norm_row.get("distributionstatus")
            or norm_row.get("status")
            or ""
        )

        low_status = native_status.lower()
        if any(tok in low_status for tok in _usda._EXCLUDED_STATUS_TOKENS):
            continue

        st_norm = _usda._normalize_us_state_name(raw_state)
        if not st_norm or st_norm not in _usda._US_FULL_STATE_NAMES:
            continue

        present_us_states.add(st_norm)
        co_norm = _usda._normalize_us_county_name(raw_county, state_norm=st_norm) if raw_county else None
        raw_parsed_rows.append((sym, st_norm, co_norm, fips, native_status))
        if co_norm:
            ck = f"{st_norm}|{co_norm}"
            counties.add(ck)
            counties_by_state.setdefault(st_norm, set()).add(co_norm)

    state_only_states: Set[str] = {
        st for st in present_us_states if not counties_by_state.get(st)
    }
    documented_states: Set[str] = set(counties_by_state.keys()) | state_only_states
    counties = _usda._expand_virginia_city_aliases(counties)

    seen_pk: Set[Tuple[str, str]] = set()
    for sym, st_norm, co_norm, fips, native_status in raw_parsed_rows:
        if co_norm:
            pk = (st_norm, co_norm)
            if pk in seen_pk:
                continue
            seen_pk.add(pk)
            db_rows.append({
                "usda_symbol": sym,
                "state_norm": st_norm,
                "county_norm": co_norm,
                "fips_code": fips,
                "native_status": native_status,
                "is_state_only": 0,
            })
        elif st_norm in state_only_states:
            pk = (st_norm, "")
            if pk in seen_pk:
                continue
            seen_pk.add(pk)
            db_rows.append({
                "usda_symbol": sym,
                "state_norm": st_norm,
                "county_norm": "",
                "fips_code": fips,
                "native_status": native_status,
                "is_state_only": 1,
            })

    return {
        "counties": counties,
        "state_only_states": state_only_states,
        "documented_states": documented_states,
        "rows": db_rows,
    }
