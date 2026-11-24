#!/usr/bin/env python3
"""
Step 02 LeafMachine3 ruler calibration SQLite/CSV ingestion, scale provenance
normalization, and provenance-preserving `RulerCalibrationMap` container.
"""

from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger("LM3_PostProcessing")


def _extract_catalog_number(raw_identifier: Any) -> str:
    """Normalizes a filename, stem, or specimen_id into a clean catalogNumber."""
    if raw_identifier is None:
        return ""
    s = str(raw_identifier).strip()
    if not s or s.lower() in {"nan", "none", "null"}:
        return ""
    stem = Path(s).stem
    cat = stem.split("__")[0].strip()
    return cat


def _resolve_scale_provenance(
    eff_source: str,
    raw_scale_src: str,
    empirical_sources: Any,
) -> Tuple[str, str, bool]:
    """Normalize `(cf_source, scale_source, is_empirical_ruler)` across row schemas."""
    src_norm = str(eff_source or "").strip().lower()
    raw_norm = str(raw_scale_src or "").strip().lower()
    emp_norm = {str(s).strip().lower() for s in (empirical_sources or ())}
    if raw_norm == "tick_lattice_unboxed" or src_norm == "tick_lattice_unboxed":
        return "tick_lattice_unboxed", "tick_lattice_unboxed", False
    if src_norm in emp_norm or raw_norm in emp_norm:
        out_cf = str(eff_source or "").strip()
        if not out_cf or src_norm not in emp_norm:
            out_cf = raw_norm if raw_norm in emp_norm else "ruler_calibrated"
        return out_cf, "ruler_calibrated", True
    if raw_norm == "mp_estimated" or src_norm in {"mp_estimated", "predicted_from_megapixels"}:
        return "predicted_from_megapixels", "mp_estimated", False
    return str(eff_source or "").strip(), "mp_estimated", False


def _parse_calibration_row(
    row: Mapping[str, Any],
    min_confidence: float = 0.50,
) -> Optional[Dict[str, Any]]:
    """Parses a single SQLite or CSV calibration record into a normalized dict."""
    from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES, normalize_pixels_per_mm

    raw_id = (
        row.get("catalogNumber")
        or row.get("catalog_number")
        or row.get("specimen_id")
        or row.get("image_name")
        or row.get("filename")
        or row.get("stem")
    )
    cat = _extract_catalog_number(raw_id)
    if not cat:
        return None

    _EMPIRICAL_RULER_SOURCES = set(_CALIBRATED_SCALE_SOURCES)
    raw_scale_src = str(row.get("scale_source") or "").strip().lower()
    cf_source = str(row.get("cf_source") or row.get("conversion_determination") or "").strip()
    if raw_scale_src == "tick_lattice_unboxed" or cf_source.lower() == "tick_lattice_unboxed":
        cf_source = "tick_lattice_unboxed"
    elif not cf_source:
        if raw_scale_src == "mp_estimated":
            cf_source = "predicted_from_megapixels"
        elif raw_scale_src == "uncalibrated":
            cf_source = "uncalibrated"
        elif raw_scale_src == "ruler_calibrated":
            cf_source = "ruler_calibrated"
    if cf_source.lower() == "uncalibrated" or raw_scale_src == "uncalibrated":
        return None
    work_scale_raw = row.get("work_scale", 1.0)
    try:
        work_scale = float(work_scale_raw) if pd.notna(work_scale_raw) and float(work_scale_raw) > 0 else 1.0
    except (ValueError, TypeError):
        work_scale = 1.0

    # Check ruler classifier confidence when present
    conf_raw = row.get(
        "ruler_class_confidence",
        row.get("ruler_confidence", row.get("confidence", np.nan)),
    )
    if pd.notna(conf_raw) and cf_source.lower() not in {
        "measured_from_fieldprism",
        "predicted_from_megapixels",
        "tick_lattice_unboxed",
    }:
        try:
            conf_val = float(conf_raw)
            if conf_val < 0:
                return None
            if conf_val > 1.0:
                conf_val /= 100.0
            if conf_val < min_confidence:
                return None
        except (ValueError, TypeError):
            pass

    # 1. Explicit mm_per_pixel / mm_per_px or direct pixels_per_mm column if provided
    mm_per_px_raw = row.get("mm_per_pixel", row.get("mm_per_px", np.nan))
    direct_px_mm = row.get("pixels_per_mm", np.nan)
    if pd.notna(mm_per_px_raw) and pd.isna(direct_px_mm):
        try:
            mpx = float(mm_per_px_raw)
            if np.isfinite(mpx) and mpx > 0.0:
                direct_px_mm = 1.0 / mpx
        except (ValueError, TypeError):
            pass
    if pd.notna(direct_px_mm):
        try:
            val_mm = normalize_pixels_per_mm(float(direct_px_mm))
            if val_mm is not None and val_mm > 0:
                eff_direct_src, out_scale_src, is_emp = _resolve_scale_provenance(
                    cf_source or "measured_from_ruler",
                    raw_scale_src,
                    _EMPIRICAL_RULER_SOURCES,
                )
                return {
                    "catalog_number": cat,
                    "pixels_per_mm": val_mm,
                    "cf_px_per_cm": val_mm * 10.0,
                    "fp_cf_px_per_cm": None,
                    "cf_source": eff_direct_src,
                    "scale_source": out_scale_src,
                    "work_scale": work_scale,
                    "is_empirical_ruler": is_emp,
                }
        except (ValueError, TypeError):
            pass

    # 2. LeafMachine3 cf_px_per_cm / fp_cf_px_per_cm / mp_cf_px_per_cm
    cf_px_cm_raw = row.get("cf_px_per_cm", np.nan)
    fp_px_cm_raw = row.get("fp_cf_px_per_cm", np.nan)
    mp_px_cm_raw = row.get("mp_cf_px_per_cm", np.nan)

    fp_val: Optional[float] = None
    if pd.notna(fp_px_cm_raw):
        try:
            if float(fp_px_cm_raw) > 0:
                fp_val = float(fp_px_cm_raw)
        except (ValueError, TypeError):
            pass

    cf_val: Optional[float] = None
    for cand_raw in (cf_px_cm_raw, fp_px_cm_raw, mp_px_cm_raw):
        if pd.notna(cand_raw):
            try:
                f_val = float(cand_raw)
                if f_val > 0:
                    cf_val = f_val
                    break
            except (ValueError, TypeError):
                pass

    if cf_val is not None and cf_val > 0:
        eff_source = cf_source
        if not eff_source:
            if fp_val is not None and abs(cf_val - fp_val) < 1e-4:
                eff_source = "measured_from_fieldprism"
            else:
                eff_source = "measured_from_ruler"
        eff_source, out_scale_src, is_emp = _resolve_scale_provenance(
            eff_source,
            raw_scale_src,
            _EMPIRICAL_RULER_SOURCES,
        )
        return {
            "catalog_number": cat,
            "pixels_per_mm": cf_val / 10.0,
            "cf_px_per_cm": cf_val,
            "fp_cf_px_per_cm": fp_val,
            "cf_source": eff_source,
            "scale_source": out_scale_src,
            "work_scale": work_scale,
            "is_empirical_ruler": is_emp,
        }

    # 3. CSV fallback columns (conversion_mean, predicted_conversion_factor_cm)
    success = str(row.get("ruler_success", "False")).strip().lower() == "true"
    conv_mean = row.get("conversion_mean", np.nan)
    pred_conv = row.get("predicted_conversion_factor_cm", np.nan)

    conversion_val: Optional[float] = None
    inferred_source = cf_source or "measured_from_ruler"
    try:
        if pd.notna(conv_mean) and float(conv_mean) > 0:
            conversion_val = float(conv_mean)
            inferred_source = cf_source or "measured_from_ruler"
        elif pd.notna(pred_conv) and float(pred_conv) > 0:
            conversion_val = float(pred_conv)
            inferred_source = cf_source or "predicted_from_megapixels"
    except (ValueError, TypeError):
        pass

    if (success and conversion_val and conversion_val > 0) or (conversion_val and conversion_val > 10.0):
        inferred_source, out_scale_src, is_emp = _resolve_scale_provenance(
            inferred_source,
            raw_scale_src,
            _EMPIRICAL_RULER_SOURCES,
        )
        return {
            "catalog_number": cat,
            "pixels_per_mm": conversion_val / 10.0,
            "cf_px_per_cm": conversion_val,
            "fp_cf_px_per_cm": fp_val,
            "cf_source": inferred_source,
            "scale_source": out_scale_src,
            "work_scale": work_scale,
            "is_empirical_ruler": is_emp,
        }

    return None


def load_ruler_calibration_records(
    lm3_dir: Path,
    min_confidence: float = 0.50,
) -> Dict[str, Dict[str, Any]]:
    """
    Loads detailed LeafMachine3 ruler / FieldPrism / MP calibration records from
    `<run_name>.sqlite` (`specimen` and `ruler_CF_lattice` tables) and
    `Data/specimen_summary.csv` / `Data/ruler_conversion_factor.csv`.
    """
    records: Dict[str, Dict[str, Any]] = {}
    if not lm3_dir or not Path(lm3_dir).exists():
        return records

    lm3_path = Path(lm3_dir)

    def _merge_rec(parsed: Optional[Dict[str, Any]]) -> None:
        if not parsed:
            return
        cat = parsed["catalog_number"]
        existing = records.get(cat)
        if existing is None:
            records[cat] = parsed
            return
        if parsed.get("is_empirical_ruler") and not existing.get("is_empirical_ruler"):
            records[cat] = parsed

    sqlite_files = sorted(
        set(list(lm3_path.glob("*.sqlite")) + list(lm3_path.glob("**/*.sqlite")))
    )
    for db_file in sqlite_files:
        try:
            conn = sqlite3.connect(f"file:{db_file.as_posix()}?mode=ro", uri=True, timeout=5.0)
            conn.row_factory = sqlite3.Row
            try:
                tables = {
                    str(r[0])
                    for r in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table';"
                    ).fetchall()
                }
                for tbl in ("specimen", "ruler_CF_lattice", "ruler_conversion_factor", "ruler"):
                    if tbl in tables:
                        for row in conn.execute(f"SELECT * FROM {tbl};").fetchall():
                            _merge_rec(_parse_calibration_row(dict(row), min_confidence=min_confidence))
            finally:
                conn.close()
        except Exception as exc:
            logger.debug("Error reading LM3 SQLite calibrations from %s: %s", db_file, exc)

    csv_files = (
        list(lm3_path.glob("Data/specimen_summary.csv"))
        + list(lm3_path.glob("**/Data/specimen_summary.csv"))
        + list(lm3_path.glob("Data/ruler_conversion_factor.csv"))
        + list(lm3_path.glob("**/Data/ruler_conversion_factor.csv"))
        + list(lm3_path.glob("**/Data/Ruler/*.csv"))
        + list(lm3_path.glob("**/Data/Batch/Ruler/*.csv"))
        + list(lm3_path.glob("**/Archival_Components/Ruler_Info/*.csv"))
    )
    seen_csvs = set()
    for rf in csv_files:
        rf_res = rf.resolve()
        if rf_res in seen_csvs:
            continue
        seen_csvs.add(rf_res)
        try:
            df = pd.read_csv(rf_res)
            if df.empty:
                continue
            for _, row in df.iterrows():
                _merge_rec(_parse_calibration_row(row.to_dict(), min_confidence=min_confidence))
        except Exception as exc:
            logger.debug("Error reading ruler CSV %s: %s", rf_res, exc)

    return records


class RulerCalibrationMap(dict):
    """
    Dictionary mapping `catalogNumber -> pixels_per_mm` (`float`) while preserving
    empirical vs. megapixel-prior provenance in `.sources` (`cf_source`) and `.records`
    so `predicted_from_megapixels` estimates are never promoted to `is_calibrated_mm = True`.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.sources: Dict[str, str] = {}
        self.records: Dict[str, Dict[str, Any]] = {}
        super().__init__(*args, **kwargs)

    def update(self, *args: Any, **kwargs: Any) -> None:  # type: ignore[override]
        if args:
            other = args[0]
            other_sources = getattr(other, "sources", None)
            other_records = getattr(other, "records", None)
            if isinstance(other_sources, dict):
                from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES

                empirical_sources = set(_CALIBRATED_SCALE_SOURCES)
                for cat_k, src_v in other_sources.items():
                    prev_src = str(self.sources.get(cat_k, "")).strip().lower()
                    prev_emp = prev_src in empirical_sources
                    new_emp = str(src_v or "").strip().lower() in empirical_sources
                    if prev_emp and not new_emp:
                        continue
                    self.sources[str(cat_k)] = str(src_v)
                    if isinstance(other, dict) and cat_k in other:
                        super().__setitem__(cat_k, other[cat_k])
                    if isinstance(other_records, dict) and cat_k in other_records:
                        self.records[str(cat_k)] = dict(other_records[cat_k])
                for cat_k, val_v in (other.items() if isinstance(other, dict) else other):
                    if str(cat_k) not in self.sources:
                        super().__setitem__(cat_k, val_v)
                if kwargs:
                    super().update(**kwargs)
                return
        super().update(*args, **kwargs)

    def copy(self) -> "RulerCalibrationMap":
        dup = RulerCalibrationMap(super().copy())
        dup.sources = dict(self.sources)
        dup.records = {k: dict(v) for k, v in self.records.items()}
        return dup


def load_ruler_calibrations(
    lm3_dir: Path,
    min_confidence: float = 0.50,
) -> Dict[str, float]:
    """
    Fast loader for specimen scale calibrations (`pixels_per_mm`) across LeafMachine3
    SQLite (`<run_name>.sqlite`) and CSV (`Data/specimen_summary.csv`,
    `Data/ruler_conversion_factor.csv`) outputs, preserving per-voucher `cf_source`
    provenance (`measured_from_ruler`, `measured_from_fieldprism`, `predicted_from_megapixels`).
    """
    records = load_ruler_calibration_records(lm3_dir=lm3_dir, min_confidence=min_confidence)
    calibrations = RulerCalibrationMap()
    for cat, rec in records.items():
        if rec.get("pixels_per_mm") is not None and float(rec["pixels_per_mm"]) > 0.0:
            calibrations[cat] = float(rec["pixels_per_mm"])
            calibrations.sources[cat] = str(rec.get("cf_source", "measured_from_ruler"))
            calibrations.records[cat] = dict(rec)
    logger.info("Loaded %d specimen ruler scale calibrations from LM3.", len(calibrations))
    return calibrations


__all__ = [
    "_extract_catalog_number",
    "_resolve_scale_provenance",
    "_parse_calibration_row",
    "load_ruler_calibration_records",
    "RulerCalibrationMap",
    "load_ruler_calibrations",
]
