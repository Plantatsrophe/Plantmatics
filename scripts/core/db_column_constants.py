"""
Core SQLite database path constants, shared migration tracking set, column sets, and scalar coercers.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Optional, Set, Tuple

from scripts.core.db_schema_env_cols import (
    ENVIRONMENTAL_BOOL_COLS,
    ENVIRONMENTAL_CACHE_COLUMNS,
    ENVIRONMENTAL_TEXT_COLS,
    HARMONIC_COLS,
    PC_COLS,
    SYM_HARMONIC_COLS,
    WC_BIO_COLUMNS,
)

logger = logging.getLogger("PackeraDB")
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SQLITE_REL_PATH = Path("data") / "plantmatics.sqlite"
LEGACY_SQLITE_REL_PATH = Path("data") / "packera_pipeline.sqlite"
_MIGRATED_LEGACY_SQLITE_PATHS: Set[str] = set()

DWC_VOUCHER_COLUMNS: Tuple[str, ...] = (
    "catalogNumber",
    "institutionCode",
    "scientificName",
    "species_raw",
    "recordedBy",
    "recordNumber",
    "identifiedBy",
    "determiner_tier",
    "type_status",
    "county",
    "stateProvince",
    "decimalLatitude",
    "decimalLongitude",
    "coordinateUncertainty",
    "year",
    "month",
    "day",
    "eventDate",
    "regional_group",
    "image_path",
    "exsiccatae_key",
    "is_primary_duplicate",
    "duplicate_count",
)

def _is_missing(val: Any) -> bool:
    if val is None:
        return True
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return True
    if isinstance(val, str) and val.strip().lower() in {"", "na", "nan", "null", "none"}:
        return True
    return False

def _to_str(val: Any) -> Optional[str]:
    if _is_missing(val):
        return None
    s = str(val).strip()
    return s if s else None

def _to_float(val: Any) -> Optional[float]:
    if _is_missing(val):
        return None
    if isinstance(val, bool):
        return 1.0 if val else 0.0
    try:
        f = float(val)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except (TypeError, ValueError):
        return None

def _to_int(val: Any, default: Optional[int] = None) -> Optional[int]:
    if _is_missing(val):
        return default
    if isinstance(val, bool):
        return 1 if val else 0
    try:
        return int(round(float(val)))
    except (TypeError, ValueError):
        return default

def _to_bool_int(val: Any, default: Optional[int] = None) -> Optional[int]:
    if _is_missing(val):
        return default
    if isinstance(val, bool):
        return 1 if val else 0
    s = str(val).strip().upper()
    if s in {"TRUE", "T", "1", "1.0", "YES", "Y"}:
        return 1
    if s in {"FALSE", "F", "0", "0.0", "NO", "N"}:
        return 0
    return default

def format_lat_lon_key(lat: Any, lon: Any) -> Optional[str]:
    f_lat = _to_float(lat)
    f_lon = _to_float(lon)
    if f_lat is None or f_lon is None:
        return None
    return f"{f_lat:.5f}_{f_lon:.5f}"

__all__ = (
    "logger", "PROJECT_ROOT", "DEFAULT_SQLITE_REL_PATH", "LEGACY_SQLITE_REL_PATH",
    "_MIGRATED_LEGACY_SQLITE_PATHS", "DWC_VOUCHER_COLUMNS", "HARMONIC_COLS", "SYM_HARMONIC_COLS",
    "PC_COLS", "WC_BIO_COLUMNS", "ENVIRONMENTAL_CACHE_COLUMNS", "ENVIRONMENTAL_TEXT_COLS",
    "ENVIRONMENTAL_BOOL_COLS", "_is_missing", "_to_str", "_to_float", "_to_int", "_to_bool_int",
    "format_lat_lon_key",
)
