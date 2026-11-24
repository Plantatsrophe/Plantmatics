"""step01_exsiccatae_fingerprint.py - Darwin Core exsiccatae normalization, geodesic distance,
image path/resolution scoring, and voucher helper predicates for Step 01."""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any, List, Mapping, Optional, Union

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]
try:
    from PIL import Image
except ImportError:
    Image = None  # type: ignore[assignment]

if Image is not None:
    Image.MAX_IMAGE_PIXELS = None

from scripts.core.config import PROJECT_ROOT


def normalize_collector(recorded_by: Optional[str]) -> str:
    """
    Normalizes botanical collector names for exsiccatae collection event fingerprinting:
    lowercased, non-alphanumeric characters stripped, initials removed, primary surname extracted.
    e.g., 'J. Brandon Fuller' -> 'fuller', 'Fuller, J. B.' -> 'fuller', 'A. Cronquist' -> 'cronquist'.
    """
    if not recorded_by or not isinstance(recorded_by, str):
        return "unknown"
    s = str(recorded_by).strip()
    if not s or s.lower() in {"none", "nan", "null", "unknown", "anonymous", "not evident"}:
        return "unknown"

    # Separate multiple collectors (e.g. 'J. B. Fuller & R. R. Kowal') -> retain primary collector
    primary = re.split(r"\s*(?:&|;|\band\b|\bwith\b|\bet\s+al\.?)\s*", s, flags=re.IGNORECASE)[0].strip()

    # Handle 'Fuller, J. Brandon' format
    if "," in primary:
        surname_part = primary.split(",")[0].strip()
        cleaned = re.sub(r"[^a-zA-Z0-9]", "", surname_part).lower()
        if cleaned:
            return cleaned

    # Tokenize and remove punctuation
    tokens = [re.sub(r"[^a-zA-Z0-9]", "", t).lower() for t in primary.split()]
    tokens = [t for t in tokens if t]

    if not tokens:
        return "unknown"

    # Remove single-letter initials (e.g. 'j', 'b')
    non_initials = [t for t in tokens if len(t) > 1]
    if not non_initials:
        # All tokens were initials (e.g. 'J. B.')
        return "".join(tokens)

    # In 'First Middle Last', the surname is the last token
    return non_initials[-1]


def normalize_record_number(record_number: Optional[Any]) -> str:
    """
    Normalizes collection number: strips boilerplate prefixes ('#', 'No.', 'Coll.'),
    preserves compound/year-prefixed numeric sequences ('2022-014' -> '2022-014') and
    ICN Art. 8.2 split-gathering letter suffixes ('1042b' -> '1042b', '1042-A' -> '1042a').
    """
    if record_number is None:
        return ""
    s = str(record_number).strip()
    if not s or s.lower() in {"none", "nan", "null", "s.n.", "sn", "s. n.", "unnumbered", "sine numero"}:
        return ""
    s_clean = re.sub(r"^(?:coll(?:ection)?\.?\s*)?(?:no\.?|nr\.?|num(?:ber)?\.?|#)\s*", "", s, flags=re.IGNORECASE).strip()
    if re.fullmatch(
        r"[A-Za-z]{2,}\d+(?:[-_./][A-Za-z0-9]+)+|[A-Za-z]{2,}\d+[A-Za-z]{2,}[A-Za-z0-9_./-]*",
        s_clean,
    ):
        return re.sub(r"[./\s_]+", "-", s_clean).lower()
    match = re.search(r"(\d+(?:[-./]\d+)*)(?:[-\s]*([a-zA-Z])(?![a-zA-Z0-9]))?", s_clean)
    if not match:
        return ""
    num_part = re.sub(r"[./]+", "-", match.group(1))
    sfx_part = (match.group(2) or "").lower()
    return f"{num_part}{sfx_part}"


def normalize_event_date(
    event_date_raw: Optional[str] = None,
    year: Optional[Any] = None,
    month: Optional[Any] = None,
    day: Optional[Any] = None,
) -> str:
    """
    Standardizes collection event date to ISO format YYYY-MM-DD.
    Evaluates explicit year, month, day components or parses date string.
    """
    try:
        if year is not None and month is not None and day is not None:
            if not (pd.isna(year) or pd.isna(month) or pd.isna(day)):
                y, m, d = int(year), int(month), int(day)
                if 1700 <= y <= 2100 and 1 <= m <= 12 and 1 <= d <= 31:
                    return f"{y:04d}-{m:02d}-{d:02d}"
    except (ValueError, TypeError):
        pass

    if event_date_raw and isinstance(event_date_raw, str):
        s = event_date_raw.strip()
        match_iso = re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", s)
        if match_iso:
            y, m, d = int(match_iso.group(1)), int(match_iso.group(2)), int(match_iso.group(3))
            return f"{y:04d}-{m:02d}-{d:02d}"
        match_y = re.search(r"\b(1[789]\d{2}|20\d{2})\b", s)
        if match_y:
            return f"{match_y.group(1)}-00-00"

    try:
        if year is not None and not pd.isna(year):
            y = int(year)
            if 1700 <= y <= 2100:
                return f"{y:04d}-00-00"
    except (ValueError, TypeError):
        pass

    return "unknown_date"


def generate_exsiccatae_key(
    collector: Optional[str],
    number: Optional[Any],
    date: Optional[str] = None,
    year: Optional[Any] = None,
    month: Optional[Any] = None,
    day: Optional[Any] = None,
) -> str:
    """
    Constructs a standardized collection event signature (exsiccatae fingerprint):
    f"{norm_collector}_{norm_number}_{norm_date}".
    """
    norm_collector = normalize_collector(collector)
    norm_number = normalize_record_number(number)
    norm_date = normalize_event_date(date, year=year, month=month, day=day)
    return f"{norm_collector}_{norm_number}_{norm_date}"


def haversine_distance_meters(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Calculates geodesic distance between two coordinate pairs using the Haversine formula."""
    r_earth = 6371000.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    a = math.sin(delta_phi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2.0) ** 2
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return r_earth * c


def resolve_voucher_image_path(
    image_path: Union[str, Path],
    workspace_dir: Optional[Path] = None,
) -> Path:
    """
    Resolves a voucher image path on disk by checking the given path, its `.webp` sibling,
    and its `.jpg`/`.jpeg`/`.png` siblings (both relative to `workspace_dir` if provided
    and directly), returning the existing file path on disk (or the `.webp` path by default
    if none exists yet).
    """
    if image_path is None:
        return Path("")
    raw_str = str(image_path).strip()
    if not raw_str:
        return Path("")
    p = Path(raw_str)
    ws_root = Path(workspace_dir) if workspace_dir is not None else None

    base_candidates: List[Path] = []
    if p.is_absolute():
        base_candidates.append(p)
    else:
        if ws_root is not None:
            base_candidates.append(ws_root / p)
        base_candidates.append(p)
        base_candidates.append(PROJECT_ROOT / p)

    for base_p in base_candidates:
        if base_p.is_file():
            return base_p
        for ext in (".webp", ".jpg", ".jpeg", ".png"):
            sib = base_p.with_suffix(ext)
            if sib != base_p and sib.is_file():
                return sib

    default_target = base_candidates[0] if base_candidates else p
    if default_target.suffix.lower() in {".jpg", ".jpeg", ".png", ".webp"}:
        return default_target.with_suffix(".webp")
    return default_target


def get_image_resolution_score(image_path_str: Optional[str], workspace_dir: Optional[Path] = None) -> float:
    """Computes a numeric image resolution score (megapixels or pixel count, fallback to byte size)."""
    if not image_path_str or not isinstance(image_path_str, str) or not image_path_str.strip():
        return 0.0
    resolved_path = resolve_voucher_image_path(image_path_str, workspace_dir=workspace_dir)
    if resolved_path.is_file():
        try:
            with Image.open(resolved_path) as img:
                w, h = img.size
                return float(w * h)
        except Exception:
            try:
                return float(resolved_path.stat().st_size)
            except Exception:
                return 0.0
    return 0.0


def _resolve_row_standardized_packera(row: Mapping[str, Any]) -> str:
    """Resolve standardized canonical 'Packera <epithet>' binomial from a voucher row."""
    from scripts.core.harvester import _extract_canonical_binomial, standardize_packera_taxon
    for col in ("species_standardized", "scientificName", "species_raw"):
        if col not in row:
            continue
        raw_val = row.get(col)
        if raw_val is None or (isinstance(raw_val, float) and math.isnan(raw_val)) or pd.isna(raw_val):
            continue
        s_txt = str(raw_val).strip()
        if not s_txt:
            continue
        std = standardize_packera_taxon(s_txt)
        if std:
            return std
        if s_txt.startswith("Packera "):
            fb = _extract_canonical_binomial(s_txt, {"Packera"})
            if fb:
                return fb
    return ""


def _is_extralimital_status(row_or_status: Any) -> bool:
    """Return True if a voucher row carries an out-of-range or exotic/waif USDA/NatureServe flag."""
    if isinstance(row_or_status, dict) or (pd is not None and isinstance(row_or_status, pd.Series)):
        raw_st = row_or_status.get("usda_county_status")
    else:
        raw_st = row_or_status
    if raw_st is None or (pd is not None and pd.isna(raw_st)):
        return False
    return str(raw_st).strip().lower() in {"out_of_range", "natureserve_exotic_or_waif"}


def _as_bool_flag(val: Any) -> bool:
    if val is None or (pd is not None and pd.isna(val)):
        return False
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in {"true", "1", "t", "yes"}
