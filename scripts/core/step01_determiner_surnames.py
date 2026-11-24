"""step01_determiner_surnames.py - Determiner surname normalization, specialist regex matching,
and ICN Art. 9.6 nomenclatural type rank/sub-rank classification for Step 01."""

from __future__ import annotations

import math
import re
from typing import Any, Iterable, List, Optional, Set, Tuple
import unicodedata

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

from scripts.core.config import PRIMARY_TYPE_STATUSES, SECONDARY_TYPE_STATUSES, SPECIALIST_PATTERNS
from scripts.core.step01_concept_authority_utils import (
    CORE_PACKERA_GOLD_REGEXES,
    CORE_PACKERA_GOLD_SURNAMES,
    REGIONAL_AND_HISTORICAL_SILVER_REGEXES,
    REGIONAL_AND_HISTORICAL_SILVER_SURNAMES,
)

GENERIC_DETERMINER_TOKENS: Set[str] = {
    "",
    "unknown",
    "anonymous",
    "none",
    "nan",
    "null",
    "na",
    "n/a",
    "staff",
    "student",
    "students",
    "herbarium",
    "curator",
    "volunteer",
    "unassigned",
    "unspecified",
    "not",
    "evident",
    "indet",
    "indetermined",
    "undetermined",
    "various",
    "class",
    "admin",
    "database",
    "digitizer",
    "illegible",
    "unrecorded",
}

_INVALID_INSTITUTION_CODES: Set[str] = {
    "",
    "UNKNOWN",
    "UNKNOWN_INST",
    "NONE",
    "NAN",
    "NULL",
    "NA",
}


def _fold_ascii(text: str) -> str:
    """Fold Unicode diacritics to base ASCII lowercase characters (e.g., 'Löve' -> 'love')."""
    norm = unicodedata.normalize("NFKD", text)
    return "".join(ch for ch in norm if not unicodedata.combining(ch)).lower()


_MISSING_IDENTIFIED_BY_TOKENS: Set[str] = {
    "", "none", "nan", "null", "na", "n/a", "<na>",
    "unknown", "unassigned", "unspecified", "unrecorded", "not recorded", "not evident",
}
_NON_PERSON_DETERMINATION_TOKENS: Tuple[str, ...] = (
    "nomenclatural adjustment",
    "nomenclature update",
    "nomenclatural update",
    "taxonomic update",
    "taxonomic adjustment",
    "synonymy update",
    "synonymy",
    "folder update",
    "herbarium folder",
    "database import",
    "batch update",
    "flora of ",
    "fna vol",
    "fna ",
    "usda plants",
    "bonap",
    "checklist",
    "kartesz",
    "[exsiccatae ",
    "sernec",
    "gbif",
    "idigbio",
    "unverified",
    "unknown",
    "anonymous",
    "leg. ign",
    "s.n.",
)
_FLORA_CITATION_REGEX = re.compile(
    r"\b(?:fide|sensu|sec\.|according\s+to|per|after|following)\s+[A-Z]"
    r"|\b(?:FNA|BONAP|USDA)\b"
    r"|^[A-Z][a-z]+(?:\s+(?:et\s+al\.?|&\s*[A-Z][a-z]+))?\s*\(?\b(?:19|20)\d{2}[a-z]?\)?$",
    re.IGNORECASE,
)
_AMBIGUOUS_COMMON_SURNAMES: Set[str] = {
    "anderson", "atwood", "baker", "bell", "britton", "brooks", "brown", "bush", "campbell", "carter",
    "chapman", "clark", "collins", "cook", "cooper", "cox", "davis", "diamond", "dorn",
    "duncan", "estes", "evans", "fisher", "freeman", "fuller", "golden", "gray", "green", "hall", "heller", "hill",
    "horn", "howard", "jackson", "james", "jeffrey", "johnson", "jones", "kelly", "king",
    "lee", "lewis", "long", "love", "martin", "miller", "mitchell", "moore", "morgan",
    "morris", "murphy", "nelson", "packer", "parker", "perry", "phillips", "powell", "price",
    "reed", "richardson", "roberts", "robinson", "rogers", "ross", "russell",
    "scott", "smith", "stewart", "taylor", "thomas", "thompson", "turner", "walker",
    "ward", "watson", "weber", "welsh", "white", "williams", "wilson", "wood", "wright", "young",
}
CORE_PACKERA_MONOGRAPHERS: Set[str] = set(CORE_PACKERA_GOLD_SURNAMES)
REGIONAL_FLORA_SPECIALISTS: Set[str] = set(REGIONAL_AND_HISTORICAL_SILVER_SURNAMES)
PACKERA_SPECIALISTS: Set[str] = CORE_PACKERA_MONOGRAPHERS | REGIONAL_FLORA_SPECIALISTS
_CORE_QUALIFIED_SPECIALIST_REGEXES: Tuple[re.Pattern, ...] = CORE_PACKERA_GOLD_REGEXES
_REGIONAL_QUALIFIED_SPECIALIST_REGEXES: Tuple[re.Pattern, ...] = REGIONAL_AND_HISTORICAL_SILVER_REGEXES
_QUALIFIED_SPECIALIST_REGEXES: Tuple[re.Pattern, ...] = (
    _CORE_QUALIFIED_SPECIALIST_REGEXES + _REGIONAL_QUALIFIED_SPECIALIST_REGEXES
)
_COLLECTOR_REFERENCE_REGEX = re.compile(
    r"^(?:(?:id(?:entified)?|det(?:ermined)?|annot(?:ated)?|ver(?:ified)?|conf(?:irmed)?)\.?\s*(?:by\s+)?)?"
    r"(?:the\s+)?(?:collector(?:s)?|same(?:\s+as\s+collector(?:s)?)?)\.?$",
    re.IGNORECASE,
)


def _has_core_monographer_match(text: str) -> bool:
    """Return True if `text` matches a core Packera / Senecioneae monographer."""
    if not text or not str(text).strip():
        return False
    s_str = str(text).strip()
    if _has_non_person_determination_token(s_str):
        return False
    surnames = set(normalize_determiner_surnames(s_str))
    if surnames & CORE_PACKERA_MONOGRAPHERS:
        return True
    return any(q_rx.search(s_str) for q_rx in _CORE_QUALIFIED_SPECIALIST_REGEXES)


def _has_regional_or_historical_specialist_match(text: str) -> bool:
    """Return True if `text` matches a regional flora or historical Senecioneae specialist."""
    if not text or not str(text).strip():
        return False
    s_str = str(text).strip()
    if _has_non_person_determination_token(s_str):
        return False
    surnames = set(normalize_determiner_surnames(s_str))
    if surnames & REGIONAL_FLORA_SPECIALISTS:
        return True
    return any(q_rx.search(s_str) for q_rx in _REGIONAL_QUALIFIED_SPECIALIST_REGEXES)


def _has_specialist_match(text: str, extra_patterns: Optional[Iterable[str]] = None) -> bool:
    """
    Returns True if `text` matches a recognized North American Packera / Senecioneae
    taxonomic specialist via normalized surnames (`PACKERA_SPECIALISTS`), initial-qualified
    specialist regexes (`_QUALIFIED_SPECIALIST_REGEXES`), `SPECIALIST_PATTERNS`, or optional
    `extra_patterns`.
    """
    if not text or not str(text).strip():
        return False
    s_str = str(text).strip()
    if _has_non_person_determination_token(s_str):
        return False
    surnames = set(normalize_determiner_surnames(s_str))
    if surnames & PACKERA_SPECIALISTS:
        return True
    for q_rx in _QUALIFIED_SPECIALIST_REGEXES:
        if q_rx.search(s_str):
            return True
    for pat in SPECIALIST_PATTERNS:
        if re.search(pat, s_str, re.IGNORECASE):
            return True
    if extra_patterns:
        for pat in extra_patterns:
            if pat and re.search(pat, s_str, re.IGNORECASE):
                return True
    return False


def _has_non_person_determination_token(val: Optional[Any]) -> bool:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return False
    if pd is not None and pd.isna(val):
        return False
    s_raw = str(val).strip()
    s_low = s_raw.lower()
    if not s_low:
        return False
    if any(tok in s_low for tok in _NON_PERSON_DETERMINATION_TOKENS):
        return True
    return bool(_FLORA_CITATION_REGEX.search(s_raw))


def _is_missing_identified_by(val: Optional[Any]) -> bool:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return True
    return str(val).strip().lower() in _MISSING_IDENTIFIED_BY_TOKENS


def _is_explicit_collector_reference(val: Optional[Any]) -> bool:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return False
    s = re.sub(r"\s+", " ", str(val).strip())
    return bool(s and _COLLECTOR_REFERENCE_REGEX.match(s))


def _should_fallback_to_collector(identified_by_val: Optional[Any]) -> bool:
    if _has_non_person_determination_token(identified_by_val):
        return False
    return _is_missing_identified_by(identified_by_val) or _is_explicit_collector_reference(identified_by_val)


def normalize_determiner_surnames(identified_by_raw: Optional[Any]) -> List[str]:
    """
    Normalizes Darwin Core `identifiedBy` strings into clean ASCII-folded determiner
    surnames, filtering out dates, prefixes, initials, and generic non-person strings.
    """
    if (
        identified_by_raw is None
        or (isinstance(identified_by_raw, float) and math.isnan(identified_by_raw))
        or _is_missing_identified_by(identified_by_raw)
        or _is_explicit_collector_reference(identified_by_raw)
        or _has_non_person_determination_token(identified_by_raw)
    ):
        return []
    raw = str(identified_by_raw).strip()
    if not raw:
        return []
    raw_low = raw.lower()
    if raw_low in GENERIC_DETERMINER_TOKENS or any(
        bad in raw_low for bad in ("herbarium staff", "not evident", "unrecorded", "anonymous", "digitiz")
    ):
        return []

    # Strip annotation prefixes, affirmation marks, parenthetical notes, dates, and years
    cleaned = re.sub(r"\([^)]*\)", " ", raw)
    cleaned = re.sub(
        r"\b(?:det\.?\s*(?:by)?|determined\s+by|annot\.?\s*(?:by)?|annotated\s+by|ver\.?\s*(?:by)?|verified\s+by|conf\.?\s*(?:by)?|confirmed\s+by|rev\.?\s*(?:by)?)\b",
        " ",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(
        r"[,\s]*(?:\d{1,2}\s+)?\b(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:t(?:ember)?)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\.?\s*\d{0,4}\b",
        " ",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(r"[,\s]*\b\d{4}(?:[-/]\d{1,2}(?:[-/]\d{1,2})?)?\b", " ", cleaned)
    cleaned = re.sub(r"\b\d+\b", " ", cleaned)
    cleaned = cleaned.replace("!", " ").strip(" ,;|-")
    if not cleaned:
        return []

    parts = re.split(r"\s*(?:&|;|\band\b|\bwith\b|\bet\s+al\.?)\s*", cleaned, flags=re.IGNORECASE)
    surnames: List[str] = []
    for part in parts:
        p = part.strip(" ,;|-.")
        if not p:
            continue
        p = re.sub(r"\b(?:dr|prof|mr|mrs|ms|jr|sr|ii|iii|iv)\b\.?", " ", p, flags=re.IGNORECASE).strip(" ,.")
        if not p:
            continue
        if "," in p:
            cand_raw = p.split(",", 1)[0].strip()
            cand_tokens = [re.sub(r"[^A-Za-zÀ-ÿ-]", "", t) for t in cand_raw.split()]
            cand_tokens = [t for t in cand_tokens if len(t) > 1]
            surname = _fold_ascii(cand_tokens[-1]) if cand_tokens else ""
        else:
            # Handle compact initials like 'Á.Löve' or 'D.K.Trock'
            p_spaced = re.sub(r"\.", ". ", p)
            tokens = [re.sub(r"[^A-Za-zÀ-ÿ-]", "", t) for t in p_spaced.split()]
            non_init = [t for t in tokens if len(t) > 1]
            surname = _fold_ascii(non_init[-1]) if non_init else ""

        surname = surname.strip("-")
        if (
            len(surname) >= 3
            and surname not in GENERIC_DETERMINER_TOKENS
            and surname not in {"collector", "collectors", "same"}
            and not any(tok in GENERIC_DETERMINER_TOKENS for tok in _fold_ascii(p).split() if tok in {"staff", "student", "students", "herbarium", "unknown", "anonymous"})
        ):
            if surname not in surnames:
                surnames.append(surname)

    return surnames


_HOLOTYPE_TIER_STATUSES: Set[str] = {"HOLOTYPE", "LECTOTYPE", "NEOTYPE", "EPITYPE"}
_ISOTYPE_TIER_STATUSES: Set[str] = {"ISOTYPE", "ISOLECTOTYPE", "ISONEOTYPE"}
_SYNTYPE_TIER_STATUSES: Set[str] = {"SYNTYPE", "ISOSYNTYPE"}
_PARATYPE_TIER_STATUSES: Set[str] = {"PARATYPE", "ISOPARATYPE", "TYPE", "COTYPE"}
_NEGATED_TYPE_RE = re.compile(r"\b(?:not\s+a?\s*type|non[-\s]?type)\b", re.IGNORECASE)


def classify_nomenclatural_type_rank(type_status_raw: Optional[Any]) -> Tuple[int, str]:
    """
    Classify a Darwin Core `typeStatus` value into a 3-tier nomenclatural type rank:
      - (2, cleaned_label) for PRIMARY_TYPE_STATUSES (Holotype, Lectotype, Neotype,
        Epitype, Isotype, Isolectotype, Isoneotype; `Tier_0_Primary_Type`)
      - (1, cleaned_label) for SECONDARY_TYPE_STATUSES (Syntype, Isosyntype,
        Paratype, Isoparatype, Type, Cotype; `Tier_0_Secondary_Type`)
      - (0, "None") for non-types or negated strings.
    """
    if type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw)):
        return 0, "None"
    if pd is not None and pd.isna(type_status_raw):
        return 0, "None"
    type_status = str(type_status_raw).strip()
    if not type_status or type_status.upper() in {
        "NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE",
        "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN",
    }:
        return 0, "None"
    if _NEGATED_TYPE_RE.search(type_status):
        return 0, "None"
    for prim_type in sorted(PRIMARY_TYPE_STATUSES, key=len, reverse=True):
        if re.search(rf"\b{prim_type}\b", type_status, re.IGNORECASE):
            return 2, prim_type.title()
    for sec_type in sorted(SECONDARY_TYPE_STATUSES, key=len, reverse=True):
        if re.search(rf"\b{sec_type}\b", type_status, re.IGNORECASE):
            return 1, sec_type.title()
    return 0, "None"


def classify_nomenclatural_type_subrank(type_status_raw: Optional[Any]) -> float:
    """
    Classify a Darwin Core `typeStatus` value into an ICN Art. 9.6 fine-grained float sub-rank:
      - 2.2: Holotype-tier (HOLOTYPE, LECTOTYPE, NEOTYPE, EPITYPE)
      - 2.0: Isotype-tier (ISOTYPE, ISOLECTOTYPE, ISONEOTYPE)
      - 1.5: Syntype-tier (SYNTYPE, ISOSYNTYPE; jointly name-bearing under ICN Art. 9.6)
      - 1.0: Paratype-tier (PARATYPE, ISOPARATYPE, TYPE, COTYPE; non-name-bearing)
      - 0.0: Non-types or negated strings.
    """
    if type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw)):
        return 0.0
    if pd is not None and pd.isna(type_status_raw):
        return 0.0
    type_status = str(type_status_raw).strip()
    if not type_status or type_status.upper() in {
        "NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE",
        "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN",
    }:
        return 0.0
    if _NEGATED_TYPE_RE.search(type_status):
        return 0.0
    for tier_val, tier_set in (
        (2.2, _HOLOTYPE_TIER_STATUSES),
        (2.0, _ISOTYPE_TIER_STATUSES),
        (1.5, _SYNTYPE_TIER_STATUSES),
        (1.0, _PARATYPE_TIER_STATUSES),
    ):
        for tok in sorted(tier_set, key=len, reverse=True):
            if re.search(rf"\b{tok}\b", type_status, re.IGNORECASE):
                return tier_val
    return 0.0


def _is_verified_type_status(type_status_raw: Optional[Any]) -> Tuple[bool, str]:
    """Return (is_type, cleaned_type_label) for a Darwin Core typeStatus value."""
    t_rank, t_lbl = classify_nomenclatural_type_rank(type_status_raw)
    return (t_rank >= 1), t_lbl
