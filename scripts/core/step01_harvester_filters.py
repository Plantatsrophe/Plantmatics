"""Step 01 canonical binomial, infraspecific, cultivation, and geographic bbox filters."""

from __future__ import annotations

import math
import re
from typing import Any, Dict, Optional, Set, Tuple

from scripts.core.step01_biogeography_utils import resolve_synonym_genera
from scripts.core.step01_concept_authority_utils import is_interspecific_hybrid_formula

_SYNONYM_GENERA: Set[str] = resolve_synonym_genera()
_OPEN_NOMENCLATURE_EPITHETS: Set[str] = {
    "sp", "spp", "cf", "aff", "indet", "hybrid", "unknown", "na", "x", "×",
}


def _extract_canonical_binomial(
    raw_name: Any,
    allowed_genera: Optional[Set[str]] = None,
) -> Optional[str]:
    """Extract canonical 'Genus epithet' binomial from a Darwin Core scientific name string."""
    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return None
    txt = str(raw_name).strip()
    if not txt:
        return None
    if is_interspecific_hybrid_formula(txt, allowed_genera=allowed_genera)[0]:
        return None
    m = re.match(r"^([A-Za-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", txt)
    if not m:
        return None
    genus = m.group(1).capitalize()
    epithet = m.group(2).lower().strip("-")
    if len(epithet) < 2 or epithet in _OPEN_NOMENCLATURE_EPITHETS:
        return None
    if allowed_genera is not None and genus not in allowed_genera:
        return None
    return f"{genus} {epithet}"


def _extract_infraspecific_taxon(raw_name: Any) -> str:
    """Extract normalized infraspecific rank + epithet ('var. gracilis', 'subsp. ...'), stripping autonyms."""
    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return ""
    txt = str(raw_name).strip()
    if not txt:
        return ""
    m = re.search(
        r"\b(subsp|ssp|var|forma|f|nothosubsp|nothovar)\.?\s+([a-z][a-z-]{1,})\b",
        txt,
        flags=re.IGNORECASE,
    )
    if not m:
        return ""
    rk_raw = m.group(1).lower()
    ep = m.group(2).lower().strip("-")
    if ep in _OPEN_NOMENCLATURE_EPITHETS:
        return ""
    m_sp = re.match(r"^[A-Za-z]+(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", txt)
    if m_sp and ep == m_sp.group(1).lower().strip("-"):
        return ""
    rk_map = {
        "ssp": "subsp.",
        "subsp": "subsp.",
        "var": "var.",
        "forma": "f.",
        "f": "f.",
        "nothosubsp": "nothosubsp.",
        "nothovar": "nothovar.",
    }
    return f"{rk_map.get(rk_raw, rk_raw + '.')} {ep}"


_CULTIVATED_STATUS_TOKENS: Set[str] = {
    "cultivated", "captive", "managed", "garden", "greenhouse", "horticultural",
}
_CULTIVATED_TEXT_RE = re.compile(
    r"(?<!not\s)(?<!non-)(?<!escaped\sfrom\s)\b"
    r"(?:cultivated|in\s+cultivation|grown\s+in\s+cultivation|greenhouse|glasshouse|"
    r"botanical\s+garden|botanic\s+garden|arboretum|hortus\s+botanicus|experimental\s+garden)\b",
    re.IGNORECASE,
)
_CA_PROVINCE_ALIASES_TO_KEY: Dict[str, str] = {
    "on": "ontario", "ont": "ontario", "ontario": "ontario",
    "qc": "quebec", "pq": "quebec", "que": "quebec", "quebec": "quebec", "québec": "quebec",
    "nb": "new brunswick", "new brunswick": "new brunswick", "nouveau-brunswick": "new brunswick",
    "ns": "nova scotia", "nova scotia": "nova scotia", "nouvelle-écosse": "nova scotia", "nouvelle-ecosse": "nova scotia",
    "pe": "prince edward island", "pei": "prince edward island", "prince edward island": "prince edward island",
    "nl": "newfoundland and labrador", "nf": "newfoundland and labrador", "newfoundland": "newfoundland and labrador", "labrador": "newfoundland and labrador", "newfoundland and labrador": "newfoundland and labrador",
    "mb": "manitoba", "man": "manitoba", "manitoba": "manitoba",
    "sk": "saskatchewan", "sask": "saskatchewan", "saskatchewan": "saskatchewan",
    "ab": "alberta", "alta": "alberta", "alberta": "alberta",
    "bc": "british columbia", "british columbia": "british columbia", "colombie-britannique": "british columbia",
    "yt": "yukon", "yukon": "yukon", "yukon territory": "yukon",
    "nt": "northwest territories", "nwt": "northwest territories", "northwest territories": "northwest territories",
    "nu": "nunavut", "nunavut": "nunavut",
}

# (lat_min, lat_max, lon_min, lon_max) envelopes for offline CoordinateCleaner state/province verification
_NA_STATE_PROVINCE_BBOXES: Dict[str, Tuple[float, float, float, float]] = {
    "alabama": (30.1, 35.1, -88.5, -84.8),
    "alaska": (51.2, 71.5, -179.9, -129.9),
    "arizona": (31.3, 37.1, -114.9, -109.0),
    "arkansas": (33.0, 36.6, -94.7, -89.6),
    "california": (32.5, 42.1, -124.5, -114.1),
    "colorado": (36.9, 41.1, -109.1, -102.0),
    "connecticut": (40.9, 42.1, -73.8, -71.7),
    "delaware": (38.4, 39.9, -75.8, -75.0),
    "district of columbia": (38.7, 39.0, -77.2, -76.9),
    "florida": (24.4, 31.1, -87.7, -79.9),
    "georgia": (30.3, 35.1, -85.7, -80.8),
    "hawaii": (18.9, 22.3, -160.3, -154.8),
    "idaho": (41.9, 49.1, -117.3, -111.0),
    "illinois": (36.9, 42.6, -91.6, -87.4),
    "indiana": (37.7, 41.8, -88.1, -84.7),
    "iowa": (40.3, 43.6, -96.7, -90.1),
    "kansas": (36.9, 40.1, -102.1, -94.5),
    "kentucky": (36.4, 39.2, -89.6, -81.9),
    "louisiana": (28.9, 33.1, -94.1, -88.8),
    "maine": (43.0, 47.5, -71.1, -66.9),
    "maryland": (37.9, 39.8, -79.5, -75.0),
    "massachusetts": (41.2, 42.9, -73.6, -69.9),
    "michigan": (41.6, 48.3, -90.5, -82.1),
    "minnesota": (43.4, 49.4, -97.3, -89.4),
    "mississippi": (30.1, 35.1, -91.7, -88.0),
    "missouri": (35.9, 40.7, -95.8, -89.0),
    "montana": (44.3, 49.1, -116.1, -104.0),
    "nebraska": (39.9, 43.1, -104.1, -95.3),
    "nevada": (35.0, 42.1, -120.1, -114.0),
    "new hampshire": (42.6, 45.4, -72.6, -70.6),
    "new jersey": (38.9, 41.4, -75.6, -73.8),
    "new mexico": (31.3, 37.1, -109.1, -103.0),
    "new york": (40.4, 45.1, -79.8, -71.8),
    "north carolina": (33.8, 36.6, -84.4, -75.4),
    "north dakota": (45.9, 49.1, -104.1, -96.5),
    "ohio": (38.4, 42.0, -84.9, -80.5),
    "oklahoma": (33.6, 37.1, -103.1, -94.4),
    "oregon": (41.9, 46.3, -124.6, -116.4),
    "pennsylvania": (39.7, 42.3, -80.6, -74.6),
    "rhode island": (41.1, 42.1, -71.9, -71.1),
    "south carolina": (32.0, 35.3, -83.4, -78.5),
    "south dakota": (42.4, 46.0, -104.1, -96.4),
    "tennessee": (34.9, 36.7, -90.4, -81.6),
    "texas": (25.8, 36.6, -106.7, -93.5),
    "utah": (36.9, 42.1, -114.1, -109.0),
    "vermont": (42.7, 45.1, -73.5, -71.4),
    "virginia": (36.5, 39.5, -83.7, -75.2),
    "washington": (45.5, 49.1, -124.8, -116.9),
    "west virginia": (37.2, 40.7, -82.7, -77.7),
    "wisconsin": (42.4, 47.1, -92.9, -86.7),
    "wyoming": (40.9, 45.1, -111.1, -104.0),
    "ontario": (41.6, 56.9, -95.2, -74.3),
    "quebec": (44.9, 62.6, -79.8, -57.1),
    "new brunswick": (44.5, 48.1, -69.1, -63.7),
    "nova scotia": (43.3, 47.1, -66.4, -59.6),
    "prince edward island": (45.9, 47.1, -64.5, -61.9),
    "newfoundland and labrador": (46.5, 60.4, -67.8, -52.6),
    "manitoba": (49.0, 60.1, -102.1, -88.9),
    "saskatchewan": (49.0, 60.1, -110.1, -101.3),
    "alberta": (49.0, 60.1, -120.1, -110.0),
    "british columbia": (48.2, 60.1, -139.1, -114.0),
    "yukon": (60.0, 69.7, -141.0, -123.8),
    "northwest territories": (60.0, 78.8, -136.5, -102.0),
    "nunavut": (51.6, 83.2, -120.7, -61.0),
}


def _normalize_ca_province_bbox_key(raw_state_prov: Optional[str]) -> Optional[str]:
    if not raw_state_prov:
        return None
    cleaned = re.sub(r"[^a-zà-ÿ\s-]", "", str(raw_state_prov).strip().lower())
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return _CA_PROVINCE_ALIASES_TO_KEY.get(cleaned)


_FATAL_GBIF_GEO_ISSUES: Set[str] = {
    "ZERO_COORDINATE",
    "COORDINATE_INVALID",
    "COORDINATE_OUT_OF_RANGE",
    "COUNTRY_COORDINATE_MISMATCH",
}


def _is_cultivated_occurrence(rec: Dict[str, Any]) -> bool:
    """Return True if Darwin Core establishment or locality metadata indicates a cultivated/garden plant."""
    for fld in ("establishmentMeans", "degreeOfEstablishment"):
        val = str(rec.get(fld) or "").strip().lower()
        if val and any(tok in val for tok in _CULTIVATED_STATUS_TOKENS) and "escaped" not in val and "naturalized" not in val:
            return True
    for txt_fld in ("locality", "verbatimLocality", "habitat", "occurrenceRemarks"):
        txt = str(rec.get(txt_fld) or "").strip()
        if txt and _CULTIVATED_TEXT_RE.search(txt):
            return True
    return False
