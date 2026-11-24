#!/usr/bin/env python3
"""Chromosome-count parsing, herbarium label extraction, and offline-safe Packera cytology taxon resolution."""
from __future__ import annotations

import os
import re
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import resolve_db_path
from scripts.core.harvester import _DYNAMIC_SYNONYM_TO_CANONICAL, standardize_packera_taxon

try:
    from scripts.data_prep.step06_synonym_utils import CANONICAL_SENECIO_PACKERA_SYNONYMS
except ImportError:
    from step06_synonym_utils import CANONICAL_SENECIO_PACKERA_SYNONYMS  # type: ignore[no-redef]

_X22: Dict[int, int] = {22 * k: k for k in range(2, 9)}
_X23: Dict[int, int] = {23 * k: k for k in range(2, 9)}
_OLD_WORLD_2N = {10, 20, 30, 40, 60, 80}
_DB_SYN_CACHE: Dict[str, Dict[str, str]] = {}

_B_CHR_RE = re.compile(r"(?:(?:\+|,|\s)\s*(?:(\d+)\s*[-–]\s*)?(\d+)|(\d+))\s*[Bb]\b")
_APPROX_RE = re.compile(r"\b(?:ca\.?|c\.?|approx\.?)\b|[~?]", re.I)
_MEIOTIC_RE = re.compile(r"(\d+)\s*(IV|III|II|I)s?\b")
_MULT_MAP: Dict[str, int] = {"IV": 4, "III": 3, "II": 2, "I": 1}
_SPLIT_RE = re.compile(r"\s*[;]\s*|\s+\bor\b\s+|\s*,\s*(?![+\s]*(?:\d+\s*[-–]\s*)?\d+\s*[Bb]\b)", re.I)

_NON_CNT = (
    r"(?!\s*(?:mm|cm|dm|m|km|in|ft|plants?|specimens?|sheets?|samples?|individuals?|"
    r"vouchers?|stems?|leaves|heads?|flowers?|collections?|sites?|populations?|plots?|"
    r"reps?|years?|branched?)\b)"
)
_RE_2N = re.compile(
    r"\b2n\s*=\s*(?:(?:ca\.?|c\.?|approx\.?|~)\s*)?\d+\b"
    r"(?:\s*(?:[-–]\s*\d+\b|\+\s*(?:\d+\s*[-–]\s*)?\d+\s*[Bb]\b|\s*[,;]\s*\d+\b" + _NON_CNT + r"))*"
    + _NON_CNT
)
_RE_N = re.compile(
    r"\bn\s*=\s*(?:(?:ca\.?|c\.?|approx\.?|~)\s*)?\d+"
    r"(?:\s*(?:IV|III|II|I)\b(?:\s*\+\s*\d+\s*(?:IV|III|II|I)\b)*"
    r"|\b(?:\s*(?:\+\s*(?:\d+\s*[-–]\s*)?\d+\s*[Bb]\b|\s*[,;]\s*\d+\b" + _NON_CNT + r"))*)"
    + _NON_CNT
)
_RE_CHR = re.compile(
    r"\bchromosome\s+(?:count|number|voucher)\b[^.;\n]{0,80}?\b"
    r"(?:2n\s*=\s*|n\s*=\s*)?(?:ca\.?\s*)?\d+(?:\s*(?:IV|III|II|I)\b|\b)" + _NON_CNT,
    re.I,
)
_RE_CVAL = re.compile(r"\b(?:1C|2C)\s*(?:=|:)?\s*\d+(?:\.\d+)?\s*pg\b", re.I)


def _default_db_synonyms() -> Dict[str, str]:
    dbp = resolve_db_path(workspace_root=PROJECT_ROOT)
    ck = f"{os.environ.get('PLANTMATICS_DB_PATH', '')}|{os.environ.get('PACKERA_DB_PATH', '')}|{dbp}"
    if ck not in _DB_SYN_CACHE:
        _DB_SYN_CACHE[ck] = {}
        if Path(dbp).exists():
            try:
                with sqlite3.connect(str(dbp)) as conn:
                    _DB_SYN_CACHE[ck] = {
                        str(k).strip().lower(): str(v).strip()
                        for k, v in conn.execute(
                            "SELECT raw_name_lower, canonical_packera FROM taxon_synonyms WHERE canonical_packera LIKE 'Packera %'"
                        )
                        if k and v
                    }
            except Exception:
                pass
    return _DB_SYN_CACHE[ck]


def _classify_cytotype(
    n_val: Optional[int], sp_2n: int, b_chr: int, approx: int, raw_str: str, force_odd: bool = False,
) -> Dict[str, Any]:
    bx, pl = (22, _X22[sp_2n]) if sp_2n in _X22 else ((23, _X23[sp_2n]) if sp_2n in _X23 else (None, None))
    is_odd = int(force_odd or pl is None or pl in (3, 5, 7) or sp_2n % 2 == 1)
    n_eff = None if (force_odd or sp_2n % 2 == 1) else (sp_2n // 2 if (pl is not None and n_val is None) else n_val)
    return {
        "raw_count_str": raw_str, "gametophytic_n": n_eff, "sporophytic_2n": sp_2n,
        "base_x": bx, "ploidy_level": pl, "cytotype_label": f"{pl}x" if pl is not None else "aneuploid",
        "b_chromosomes": b_chr, "is_aneuploid_or_odd_ploid": is_odd, "is_approximate": approx,
    }


def _parse_single_text(text: str, default_mode: Optional[str] = None) -> List[Dict[str, Any]]:
    raw_orig = " ".join(str(text or "").split())
    if not raw_orig:
        return []
    cleaned = raw_orig
    for _ in range(4):
        cleaned = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", " ", cleaned)
    if not cleaned.strip(" -?"):
        return []
    g_approx, out, inherited = int(bool(_APPROX_RE.search(cleaned))), [], default_mode
    for piece in (p.strip() for p in _SPLIT_RE.split(cleaned) if p.strip()):
        approx, bm = int(bool(g_approx or _APPROX_RE.search(piece))), _B_CHR_RE.search(piece)
        b_chr = int(bm.group(2) or bm.group(3) or 0) if bm else 0
        piece = _APPROX_RE.sub("", _B_CHR_RE.sub("", piece) if bm else piece).strip()
        mode_cur = inherited
        if re.match(r"^\s*2n\b", piece, re.I):
            mode_cur = inherited = "2n"
            piece = re.sub(r"^\s*2n\s*=?\s*", "", piece, flags=re.I).strip()
        elif re.match(r"^\s*n\b", piece, re.I):
            mode_cur = inherited = "n"
            piece = re.sub(r"^\s*n\s*=?\s*", "", piece, flags=re.I).strip()
        m_terms = _MEIOTIC_RE.findall(piece)
        if m_terms:
            sp_2n = sum(int(c) * _MULT_MAP[s] for c, s in m_terms)
            has_un = any(s in ("I", "III") for _, s in m_terms)
            if 8 <= sp_2n <= 250:
                out.append(_classify_cytotype(None if has_un or sp_2n % 2 else sp_2n // 2, sp_2n, b_chr, approx, raw_orig, has_un))
            continue
        rng = re.match(r"^\s*(\d+)\s*[-–]\s*(\d+)\s*$", piece)
        for val in ([int(rng.group(1)), int(rng.group(2))] if rng else [int(x) for x in re.findall(r"\b\d+\b", piece)]):
            m = mode_cur or ("n" if val < 36 else "2n")
            if m == "n" and 4 <= val <= 125:
                out.append(_classify_cytotype(val, val * 2, b_chr, approx or int(bool(rng)), raw_orig))
            elif m == "2n" and 8 <= val <= 250:
                out.append(_classify_cytotype(None if val % 2 else val // 2, val, b_chr, approx or int(bool(rng)), raw_orig))
    return out


def parse_packera_chromosome_record(
    raw_count: Any, gametophytic_raw: Any = None, sporophytic_raw: Any = None,
) -> List[Dict[str, Any]]:
    """Parse raw chromosome strings, gametophytic n, or sporophytic 2n counts into normalized cytotype dicts."""
    raw_s, g_s, s_s = (" ".join(str(x or "").split()) for x in (raw_count, gametophytic_raw, sporophytic_raw))
    if raw_s:
        parsed = _parse_single_text(raw_s)
        if parsed:
            return parsed
    g_p, s_p = (_parse_single_text(g_s, "n") if g_s else []), (_parse_single_text(s_s, "2n") if s_s else [])
    raw_fb = raw_s or "; ".join(p for p, v in ((f"n={g_s}", g_p), (f"2n={s_s}", s_p)) if v)
    if len(g_p) == 1 and len(s_p) == 1 and g_p[0]["sporophytic_2n"] == s_p[0]["sporophytic_2n"]:
        return [dict(s_p[0], gametophytic_n=g_p[0]["gametophytic_n"], raw_count_str=raw_fb)]
    seen, uniq = set(), []
    for rec in s_p + g_p:
        k = (rec["sporophytic_2n"], rec["b_chromosomes"])
        if k not in seen:
            seen.add(k)
            rec["raw_count_str"] = raw_fb
            uniq.append(rec)
    return uniq


def _strip_authors(name: str) -> str:
    s = re.split(r"\s+\bas\b\s+", re.sub(r"\s*\[[^\]]*\]", " ", str(name or "")), maxsplit=1, flags=re.I)[0]
    for _ in range(4):
        s = re.sub(r"\([^()]*\)", " ", s)
    toks = [t for t in s.replace("×", " ").split() if t]
    out, i = toks[:2], 2
    while i < len(toks) - 1 and toks[i].lower() in ("var.", "subsp.", "ssp.", "f."):
        out.extend([toks[i].lower(), toks[i + 1]])
        i += 2
    return " ".join(out)


def resolve_packera_cytology_taxon(
    raw_name: Any, synonym_map: Optional[Mapping[str, str]] = None,
    sporophytic_2n: Optional[int] = None, allow_standardize_fallback: bool = True,
) -> Optional[Dict[str, Optional[str]]]:
    """Resolve a Packera or Senecio/Cineraria basionym taxon name using canonical synonyms before DB overlays."""
    orig = " ".join(str(raw_name or "").split())
    clean = _strip_authors(orig)
    words = clean.split()
    if len(words) < 2:
        return None
    smap: Dict[str, str] = {
        str(k).strip().lower(): v
        for k, v in CANONICAL_SENECIO_PACKERA_SYNONYMS.items()
        if isinstance(v, str) and v.startswith("Packera ")
    }
    for src in (_default_db_synonyms(), _DYNAMIC_SYNONYM_TO_CANONICAL, synonym_map or {}):
        smap.update({str(k).strip().lower(): v for k, v in src.items() if isinstance(v, str) and v.startswith("Packera ")})
    binom, is_sen = f"{words[0].lower()} {words[1].lower()}", words[0].lower() in ("senecio", "cineraria")
    if re.search(r"\bledeb\.?\b", orig, re.I) or ((is_sen or "tomentos" in orig.lower()) and sporophytic_2n is not None and int(sporophytic_2n) in _OLD_WORLD_2N):
        return None
    if clean.lower().startswith("senecio aureus var. pauperculus"):
        return {"original_taxon_name": orig, "basionym_name": clean, "species_standardized": "Packera paupercula"}
    for cand in (orig.lower(), clean.lower(), binom):
        if cand in smap:
            return {"original_taxon_name": orig, "basionym_name": clean if is_sen else None, "species_standardized": smap[cand]}
    if words[0].lower() == "packera":
        std = standardize_packera_taxon(clean) or f"Packera {words[1].lower()}"
        return {"original_taxon_name": orig, "basionym_name": None, "species_standardized": std} if std.startswith("Packera ") else None
    if is_sen and allow_standardize_fallback and not (sporophytic_2n is not None and int(sporophytic_2n) in _OLD_WORLD_2N):
        std = standardize_packera_taxon(clean)
        if std and std.startswith("Packera "):
            return {"original_taxon_name": orig, "basionym_name": clean, "species_standardized": std}
    return None


def extract_label_cytology_mentions(text: Any) -> List[str]:
    """Extract non-overlapping chromosome count and 1C/2C flow cytometry mentions from label prose."""
    s = " ".join(str(text or "").split())
    spans: List[Tuple[int, int, str]] = []
    for rx in (_RE_2N, _RE_N, _RE_CHR, _RE_CVAL):
        for m in rx.finditer(s):
            a, b = m.span()
            if not any(not (b <= ea or a >= eb) for ea, eb, _ in spans):
                spans.append((a, b, m.group(0).strip()))
    return list(dict.fromkeys(t for _, _, t in sorted(spans, key=lambda x: x[0]))) if s else []
