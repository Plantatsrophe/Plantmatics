"""
Module: usda_plants_ranges.py
Project: Packera Species Delimitation & Morphometrics Pipeline (NCU Herbarium)
Description:
    Facade for USDA NRCS PLANTS county distribution fetching, synonym-aware symbol
    resolution, SQLite/CSV caching, spatial county adjacency, and Step 01 range filtering.

    Official Public-Domain Endpoints:
      1. Master Checklist (`plantlst.txt`):
         https://plants.sc.egov.usda.gov/DocumentLibrary/Txt/plantlst.txt
      2. Symbol Search (`POST`):
         https://plantsservices.sc.egov.usda.gov/api/PlantSearch
      3. Per-Symbol Distribution CSV (`GET`):
         https://plantsservices.sc.egov.usda.gov/api/PlantProfileDataDownload?symbol=<SYMBOL>
"""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any, Dict, Mapping, Optional, Set, Tuple
import urllib.error
import urllib.request

from scripts.core.config import resolve_target_genus  # noqa: F401
from scripts.core.db import (  # noqa: F401
    get_db_connection,
    is_canonical_or_explicit_db_context,
    load_taxon_synonyms,
    resolve_db_path,
)
from scripts.core.epa_ecoregions import (  # noqa: F401
    DEFAULT_EPA_L3_GEOJSON,
    DEFAULT_US_COUNTIES_GEOJSON,
    _CANADIAN_PROV_ABBR_TO_FULL,
    _US_STATE_ABBR_TO_FULL,
    _ensure_epa_level3_cache_files,
    _is_canadian_province_label,
    _normalize_us_county_name,
    _normalize_us_state_name,
    load_epa_level3_spatial_index,
)
from scripts.core.io_utils import atomic_write_text  # noqa: F401
from scripts.core.taxonomy_authority import (  # noqa: F401
    classify_nomenclatural_type_rank,
    haversine_distance_meters,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_USDA_CACHE_DIR = PROJECT_ROOT / "data" / "environmental" / "usda_plants"
USDA_PLANTLST_URL = "https://plants.sc.egov.usda.gov/DocumentLibrary/Txt/plantlst.txt"
USDA_PLANT_SEARCH_API_URL = "https://plantsservices.sc.egov.usda.gov/api/PlantSearch"
USDA_PLANT_PROFILE_API_URL = "https://plantsservices.sc.egov.usda.gov/api/PlantProfile"
USDA_DISTRIBUTION_DOCS_API_URL = (
    "https://plantsservices.sc.egov.usda.gov/api/PlantProfile/getDownloadDistributionDocumentation"
)
USDA_DISTRIBUTION_CSV_API_URL = "https://plantsservices.sc.egov.usda.gov/api/PlantProfileDataDownload"
USDA_USER_AGENT = "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"

_US_FULL_STATE_NAMES: Set[str] = set(_US_STATE_ABBR_TO_FULL.values())
_US_FULL_TO_ABBR: Dict[str, str] = {full: abbr.upper() for abbr, full in _US_STATE_ABBR_TO_FULL.items()}
_CA_FULL_TO_ABBR: Dict[str, str] = {
    "ontario": "ON", "quebec": "QC", "new brunswick": "NB", "nova scotia": "NS",
    "prince edward island": "PE", "newfoundland and labrador": "NL", "manitoba": "MB",
    "saskatchewan": "SK", "alberta": "AB", "british columbia": "BC",
    "yukon": "YT", "northwest territories": "NT", "nunavut": "NU",
}
_SUBNATION_EDGE_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("AL", "FL"), ("AL", "GA"), ("AL", "MS"), ("AL", "TN"),
    ("AK", "YT"), ("AK", "BC"),
    ("AZ", "CA"), ("AZ", "CO"), ("AZ", "NM"), ("AZ", "NV"), ("AZ", "UT"),
    ("AR", "LA"), ("AR", "MO"), ("AR", "MS"), ("AR", "OK"), ("AR", "TN"), ("AR", "TX"),
    ("CA", "NV"), ("CA", "OR"),
    ("CO", "KS"), ("CO", "NE"), ("CO", "NM"), ("CO", "OK"), ("CO", "UT"), ("CO", "WY"),
    ("CT", "MA"), ("CT", "NY"), ("CT", "RI"),
    ("DE", "MD"), ("DE", "NJ"), ("DE", "PA"), ("DC", "MD"), ("DC", "VA"),
    ("FL", "GA"), ("GA", "NC"), ("GA", "SC"), ("GA", "TN"),
    ("ID", "MT"), ("ID", "NV"), ("ID", "OR"), ("ID", "UT"), ("ID", "WA"), ("ID", "WY"), ("ID", "BC"), ("ID", "AB"),
    ("IL", "IN"), ("IL", "IA"), ("IL", "KY"), ("IL", "MO"), ("IL", "WI"),
    ("IN", "KY"), ("IN", "MI"), ("IN", "OH"), ("IA", "MN"), ("IA", "MO"), ("IA", "NE"), ("IA", "SD"), ("IA", "WI"),
    ("KS", "MO"), ("KS", "NE"), ("KS", "OK"),
    ("KY", "MO"), ("KY", "OH"), ("KY", "TN"), ("KY", "VA"), ("KY", "WV"),
    ("LA", "MS"), ("LA", "TX"),
    ("ME", "NH"), ("ME", "NB"), ("ME", "QC"),
    ("MD", "PA"), ("MD", "VA"), ("MD", "WV"),
    ("MA", "NH"), ("MA", "NY"), ("MA", "RI"), ("MA", "VT"),
    ("MI", "MN"), ("MI", "OH"), ("MI", "WI"), ("MI", "ON"),
    ("MN", "ND"), ("MN", "SD"), ("MN", "WI"), ("MN", "MB"), ("MN", "ON"),
    ("MS", "TN"), ("MO", "NE"), ("MO", "OK"), ("MO", "TN"),
    ("MT", "ND"), ("MT", "SD"), ("MT", "WY"), ("MT", "AB"), ("MT", "BC"), ("MT", "SK"),
    ("NE", "SD"), ("NE", "WY"), ("NV", "OR"), ("NV", "UT"),
    ("NH", "VT"), ("NH", "QC"), ("NJ", "NY"), ("NJ", "PA"),
    ("NM", "OK"), ("NM", "TX"), ("NM", "UT"),
    ("NY", "PA"), ("NY", "VT"), ("NY", "ON"), ("NY", "QC"),
    ("NC", "SC"), ("NC", "TN"), ("NC", "VA"),
    ("ND", "SD"), ("ND", "MB"), ("ND", "SK"),
    ("OH", "PA"), ("OH", "WV"), ("OH", "ON"), ("OK", "TX"),
    ("OR", "WA"), ("PA", "WV"), ("PA", "ON"),
    ("SD", "WY"), ("TN", "VA"), ("UT", "WY"),
    ("VT", "QC"), ("VA", "WV"), ("WA", "BC"),
    ("YT", "BC"), ("YT", "NT"), ("BC", "AB"), ("BC", "NT"),
    ("AB", "SK"), ("AB", "NT"), ("SK", "MB"), ("SK", "NT"), ("SK", "NU"),
    ("MB", "ON"), ("MB", "NU"), ("ON", "QC"), ("QC", "NB"), ("QC", "NL"), ("NB", "NS"), ("NB", "PE"), ("NS", "PE"),
)
_SUBNATION_ADJACENCY: Dict[str, Set[str]] = {}
for _s_a, _s_b in _SUBNATION_EDGE_PAIRS:
    _SUBNATION_ADJACENCY.setdefault(_s_a, set()).add(_s_b)
    _SUBNATION_ADJACENCY.setdefault(_s_b, set()).add(_s_a)
_EXCLUDED_STATUS_TOKENS = ("excluded", "(x)", "false report", "erroneous")

# Process-level caches
_PLANTLST_INDEX_CACHE: Dict[Tuple[str, float], Dict[str, Any]] = {}
_COUNTY_SPATIAL_GRAPH_CACHE: Dict[Tuple[str, float], Dict[str, Any]] = {}
_USDA_NETWORK_UNAVAILABLE: bool = False


def reset_usda_network_state() -> None:
    """Reset cached network-unavailability state (used by unit tests)."""
    global _USDA_NETWORK_UNAVAILABLE
    _USDA_NETWORK_UNAVAILABLE = False


def _extract_clean_binomial(raw_name: Any) -> Tuple[Optional[str], bool]:
    """
    Extract lowercase `'genus epithet'` binomial and whether the string specifies
    an infraspecific rank (`var.`, `subsp.`, `ssp.`, `f.`).
    """
    if raw_name is None:
        return None, False
    txt = str(raw_name).strip()
    if not txt:
        return None, False
    m = re.match(r"^([A-Za-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", txt)
    if not m:
        return None, False
    genus = m.group(1).lower()
    epithet = m.group(2).lower().strip("-")
    if len(epithet) < 2 or epithet in {"sp", "spp", "cf", "aff", "indet", "hybrid", "unknown", "na", "x"}:
        return None, False
    is_infraspecific = bool(re.search(r"\b(?:var|subsp|ssp|f)\.?\s+[a-z-]{2,}", txt.lower()))
    return f"{genus} {epithet}", is_infraspecific


def _http_get_text(url: str, timeout: float = 12.0) -> Optional[str]:
    """Perform an unauthenticated HTTP GET returning decoded UTF-8 text or None on error."""
    global _USDA_NETWORK_UNAVAILABLE
    if _USDA_NETWORK_UNAVAILABLE:
        return None
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": USDA_USER_AGENT,
            "Accept": "application/json, text/csv, text/plain, */*",
        },
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if getattr(resp, "status", 200) != 200:
                return None
            raw_bytes = resp.read()
            return raw_bytes.decode("utf-8-sig", errors="replace")
    except urllib.error.HTTPError:
        return None
    except Exception:
        _USDA_NETWORK_UNAVAILABLE = True
        return None


def _http_post_text(url: str, payload: Mapping[str, Any], timeout: float = 15.0) -> Optional[str]:
    """Perform an unauthenticated HTTP POST with JSON payload returning raw UTF-8 text (CSV or JSON)."""
    global _USDA_NETWORK_UNAVAILABLE
    if _USDA_NETWORK_UNAVAILABLE:
        return None
    body_bytes = json.dumps(dict(payload)).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=body_bytes,
        headers={
            "User-Agent": USDA_USER_AGENT,
            "Content-Type": "application/json;charset=UTF-8",
            "Accept": "text/csv, application/json, text/plain, */*",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if getattr(resp, "status", 200) != 200:
                return None
            return resp.read().decode("utf-8-sig", errors="replace")
    except urllib.error.HTTPError:
        return None
    except Exception:
        _USDA_NETWORK_UNAVAILABLE = True
        return None


def _http_post_json(url: str, payload: Mapping[str, Any], timeout: float = 10.0) -> Optional[Dict[str, Any]]:
    """Perform an unauthenticated HTTP POST with JSON payload returning parsed dict or None."""
    txt = _http_post_text(url, payload=payload, timeout=timeout)
    if not txt:
        return None
    try:
        parsed = json.loads(txt)
        return parsed if isinstance(parsed, dict) else None
    except Exception:
        return None


from scripts.core.step01_county_spatial_graph import (  # noqa: E402, F401
    _ensure_standalone_us_counties_geojson,
    _is_within_geodesic_buffer,
    load_us_county_spatial_graph,
    resolve_voucher_us_county,
)
from scripts.core.step01_usda_csv_parser import (  # noqa: E402, F401
    _convert_usda_json_records_to_csv,
    _normalize_usda_distribution_payload,
    _slice_csv_from_header,
    parse_usda_distribution_csv,
)
from scripts.core.step01_usda_distribution_io import (  # noqa: E402, F401
    _fetch_usda_symbol_distribution_csv,
    _load_taxon_range_from_sqlite,
    _persist_taxon_range_to_sqlite,
    fetch_or_load_usda_taxon_range,
)
from scripts.core.step01_usda_natureserve_utils import (  # noqa: E402, F401
    _evaluate_natureserve_voucher,
    _load_natureserve_ranges_for_taxa,
    _resolve_voucher_country_and_subnation,
)
from scripts.core.step01_usda_symbol_resolver import (  # noqa: E402, F401
    VA_INDEPENDENT_CITIES,
    _EAST_TX_PINEY_WOODS_COUNTIES,
    _augment_documented_counties,
    _expand_virginia_city_aliases,
    _generate_algorithmic_generic_transfers,
    _parse_plantlst_text,
    load_or_fetch_plantlst_index,
    normalize_virginia_independent_city,
    resolve_usda_symbols_for_taxon,
)
from scripts.core.step01_usda_voucher_filter import (  # noqa: E402, F401
    filter_vouchers_by_usda_county_range,
)
