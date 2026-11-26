"""
Step 01 taxon synonym and Darwin Core voucher SQLite persistence plus CSV reconciliation.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import io
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

from scripts.core.db_column_constants import (
    DWC_VOUCHER_COLUMNS,
    _is_missing,
    _to_bool_int,
    _to_float,
    _to_int,
    _to_str,
)
from scripts.core.db_conn_pragmas import get_table_columns, upsert_rows_generic
from scripts.core.step01_concept_authority_utils import parse_voucher_taxonomic_concept

def upsert_taxon_synonyms(
    conn: sqlite3.Connection,
    synonym_to_canonical: Mapping[str, str],
    author_strings_by_canonical: Optional[Mapping[str, Iterable[str]]] = None,
) -> int:
    """Persist discovered GBIF / configuration taxon synonym mappings into `taxon_synonyms`."""
    if not synonym_to_canonical:
        return 0
    now_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    authors_map: Dict[str, str] = {}
    if author_strings_by_canonical:
        for canon, auths in author_strings_by_canonical.items():
            clean_list = sorted({str(a).strip() for a in (auths or []) if str(a).strip()})
            if clean_list:
                authors_map[str(canon).strip()] = json.dumps(clean_list, ensure_ascii=False)

    payload: List[Tuple[str, str, Optional[str], str]] = []
    for raw_low, canon in synonym_to_canonical.items():
        r_key = re.sub(r"\s+", " ", str(raw_low or "").strip().lower())
        c_val = str(canon or "").strip()
        if not r_key or not c_val:
            continue
        payload.append((r_key, c_val, authors_map.get(c_val), now_iso))

    if not payload:
        return 0

    sql = """
    INSERT INTO taxon_synonyms (raw_name_lower, canonical_packera, protologue_authors_json, updated_at)
    VALUES (?, ?, ?, ?)
    ON CONFLICT(raw_name_lower) DO UPDATE SET
        canonical_packera = excluded.canonical_packera,
        protologue_authors_json = COALESCE(excluded.protologue_authors_json, taxon_synonyms.protologue_authors_json),
        updated_at = excluded.updated_at;
    """
    with conn:
        conn.executemany(sql, payload)
    return len(payload)

def load_taxon_synonyms(conn: sqlite3.Connection) -> Dict[str, str]:
    """Load persisted `raw_name_lower -> canonical_packera` dictionary from `taxon_synonyms`."""
    cur = conn.execute("SELECT raw_name_lower, canonical_packera FROM taxon_synonyms;")
    return {str(r[0]): str(r[1]) for r in cur.fetchall() if r[0] and r[1]}

def _infer_species_standardized_fallback(sci_name: Optional[str], sp_raw: Optional[str], syn_map: Mapping[str, str]) -> Optional[str]:
    for cand in (sci_name, sp_raw):
        if not cand:
            continue
        low = re.sub(r"\s+", " ", cand.strip().lower())
        if low in syn_map:
            return syn_map[low]
        m = re.match(r"^([A-Za-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", cand.strip())
        if m:
            genus = m.group(1).capitalize()
            ep = m.group(2).lower().strip("-")
            binom_low = f"{genus.lower()} {ep}"
            if binom_low in syn_map:
                return syn_map[binom_low]
            if genus == "Packera" and len(ep) >= 2:
                return f"Packera {ep}"
    return None

def upsert_vouchers_records(
    conn: sqlite3.Connection,
    records: Sequence[Mapping[str, Any]],
) -> int:
    """
    UPSERT Darwin Core voucher rows into `vouchers`, preserving existing Step 02
    calibration and QC columns when not explicitly provided in `records`.
    """
    if not records:
        return 0
    syn_map = load_taxon_synonyms(conn)
    tbls = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()}
    bb_map = {str(r[0]): {"basionym_name": r[1], "homotypic_synonyms_json": r[2]} for r in conn.execute("SELECT species_standardized, basionym_name, homotypic_synonyms_json FROM taxon_backbone_profiles;").fetchall() if r[0]} if "taxon_backbone_profiles" in tbls else {}
    v_cols = get_table_columns(conn, "vouchers")
    ex_conc = {str(r[0]): (str(r[1] or "").strip(), str(r[2] or "").strip()) for r in conn.execute("SELECT catalogNumber, concept_synonymy_status, typified_basionym FROM vouchers;").fetchall() if r[0]} if {"catalogNumber", "concept_synonymy_status", "typified_basionym"}.issubset(v_cols) else {}
    norm_rows: List[Dict[str, Any]] = []

    for rec in records:
        cat = _to_str(rec.get("catalogNumber"))
        if not cat:
            continue
        sci = _to_str(rec.get("scientificName")) or _to_str(rec.get("species_raw"))
        sp_raw = _to_str(rec.get("species_raw")) or sci
        dec_lat = _to_float(rec.get("decimalLatitude"))
        dec_lon = _to_float(rec.get("decimalLongitude"))
        ident = _to_str(rec.get("identifiedBy"))
        sp_std = _to_str(rec.get("species_standardized")) or _infer_species_standardized_fallback(sci, sp_raw, syn_map)

        row_dict: Dict[str, Any] = {
            "catalogNumber": cat,
            "institutionCode": _to_str(rec.get("institutionCode")),
            "scientificName": sci,
            "species_raw": sp_raw,
            "recordedBy": _to_str(rec.get("recordedBy")),
            "recordNumber": _to_str(rec.get("recordNumber")),
            "identifiedBy": ident,
            "determiner_tier": _to_str(rec.get("determiner_tier")),
            "type_status": _to_str(rec.get("type_status")),
            "county": _to_str(rec.get("county")),
            "stateProvince": _to_str(rec.get("stateProvince")),
            "locality": _to_str(rec.get("locality")) or _to_str(" ; ".join(x for x in (str(rec.get("_raw_locality") or "").strip(), str(rec.get("habitat") or rec.get("_raw_habitat") or "").strip()) if x)),
            "decimalLatitude": dec_lat,
            "decimalLongitude": dec_lon,
            "coordinateUncertainty": _to_float(rec.get("coordinateUncertainty")),
            "year": _to_int(rec.get("year")),
            "month": _to_int(rec.get("month")),
            "day": _to_int(rec.get("day")),
            "eventDate": _to_str(rec.get("eventDate")),
            "doy": _to_float(rec.get("doy")),
            "regional_group": _to_str(rec.get("regional_group")),
            "image_path": _to_str(rec.get("image_path")),
            "exsiccatae_key": _to_str(rec.get("exsiccatae_key")),
            "species_standardized": sp_std,
            "pixels_per_mm": _to_float(rec.get("pixels_per_mm")),
            "is_calibrated_mm": _to_bool_int(rec.get("is_calibrated_mm")),
        }
        if "is_primary_duplicate" in rec and not _is_missing(rec.get("is_primary_duplicate")):
            row_dict["is_primary_duplicate"] = _to_bool_int(rec.get("is_primary_duplicate"), default=1)
        if "duplicate_count" in rec and not _is_missing(rec.get("duplicate_count")):
            row_dict["duplicate_count"] = _to_int(rec.get("duplicate_count"), default=1)
        if "exsiccatae_taxonomic_conflict" in rec and not _is_missing(rec.get("exsiccatae_taxonomic_conflict")):
            row_dict["exsiccatae_taxonomic_conflict"] = _to_bool_int(rec.get("exsiccatae_taxonomic_conflict"), default=0)
        if "exsiccatae_gold_taxon" in rec and not _is_missing(rec.get("exsiccatae_gold_taxon")):
            row_dict["exsiccatae_gold_taxon"] = _to_str(rec.get("exsiccatae_gold_taxon"))
        if "usda_symbol" in rec and not _is_missing(rec.get("usda_symbol")):
            row_dict["usda_symbol"] = _to_str(rec.get("usda_symbol"))
        if "usda_county_status" in rec and not _is_missing(rec.get("usda_county_status")):
            row_dict["usda_county_status"] = _to_str(rec.get("usda_county_status"))
        if "infraspecific_taxon" in rec and not _is_missing(rec.get("infraspecific_taxon")):
            row_dict["infraspecific_taxon"] = _to_str(rec.get("infraspecific_taxon"))
        if "georef_qc_status" in rec and not _is_missing(rec.get("georef_qc_status")):
            row_dict["georef_qc_status"] = _to_str(rec.get("georef_qc_status"))
        if "is_centroid_or_inst_coord" in rec and not _is_missing(rec.get("is_centroid_or_inst_coord")):
            row_dict["is_centroid_or_inst_coord"] = _to_bool_int(rec.get("is_centroid_or_inst_coord"), default=0)
        if "coord_quarantine_reason" in rec and not _is_missing(rec.get("coord_quarantine_reason")):
            row_dict["coord_quarantine_reason"] = _to_str(rec.get("coord_quarantine_reason"))
        if "typified_basionym" in rec and not _is_missing(rec.get("typified_basionym")):
            row_dict["typified_basionym"] = _to_str(rec.get("typified_basionym"))
        if "verbatim_taxon_concept" in rec and not _is_missing(rec.get("verbatim_taxon_concept")):
            row_dict["verbatim_taxon_concept"] = _to_str(rec.get("verbatim_taxon_concept"))
        if "concept_synonymy_status" in rec and not _is_missing(rec.get("concept_synonymy_status")):
            row_dict["concept_synonymy_status"] = _to_str(rec.get("concept_synonymy_status"))
        if "identification_year" in rec and not _is_missing(rec.get("identification_year")):
            row_dict["identification_year"] = _to_int(rec.get("identification_year"))
        if "is_open_or_hybrid_concept" in rec and not _is_missing(rec.get("is_open_or_hybrid_concept")):
            row_dict["is_open_or_hybrid_concept"] = _to_bool_int(rec.get("is_open_or_hybrid_concept"), default=0)
        ts_clean = str(row_dict.get("type_status") or "").strip()
        is_t_row = bool(ts_clean and ts_clean.upper() not in {"NONE", "NOT A TYPE", "NON-TYPE"})
        prev_cs, prev_tb = ex_conc.get(cat, ("", ""))
        need_conc = bool((sci or sp_raw or is_t_row) and "concept_synonymy_status" not in row_dict and (not prev_cs or (is_t_row and "typified_basionym" not in row_dict and not prev_tb)))
        if need_conc:
            ci = parse_voucher_taxonomic_concept(rec, sp_std, is_t_row, {"Packera", "Senecio", "Cineraria"}, bb_map)
            row_dict.setdefault("concept_synonymy_status", ci["concept_synonymy_status"])
            row_dict.setdefault("verbatim_taxon_concept", ci["verbatim_taxon_concept"])
            row_dict.setdefault("is_open_or_hybrid_concept", 1 if ci["is_open_or_hybrid_concept"] else 0)
            if ci["typified_basionym"]:
                row_dict.setdefault("typified_basionym", ci["typified_basionym"])
            if ci["identification_year"] is not None:
                row_dict.setdefault("identification_year", int(ci["identification_year"]))
        if "spatial_balance_weight" in rec and not _is_missing(rec.get("spatial_balance_weight")):
            row_dict["spatial_balance_weight"] = _to_float(rec.get("spatial_balance_weight"))
        qc_val = _to_str(rec.get("qc_status"))
        if qc_val in {"PASS", "FAILED_QC", "PENDING"}:
            row_dict["qc_status"] = qc_val
        if "failure_reason" in rec and not _is_missing(rec.get("failure_reason")):
            row_dict["failure_reason"] = _to_str(rec.get("failure_reason"))
        if "failure_details" in rec and not _is_missing(rec.get("failure_details")):
            row_dict["failure_details"] = _to_str(rec.get("failure_details"))
        norm_rows.append(row_dict)

    return upsert_rows_generic(
        conn,
        table_name="vouchers",
        rows=norm_rows,
        pk_cols=("catalogNumber",),
        preserve_existing_on_null=True,
    )

_CSV_BOOL_VOUCHER_COLUMNS = frozenset({
    "is_primary_duplicate", "exsiccatae_taxonomic_conflict", "is_calibrated_mm",
    "is_centroid_or_inst_coord", "is_open_or_hybrid_concept",
})

def _format_voucher_csv_value(col: str, val: Any) -> Any:
    """Render a SQLite `vouchers` value in the pandas-style `curated_vouchers.csv` convention."""
    if _is_missing(val):
        return ""
    if col in _CSV_BOOL_VOUCHER_COLUMNS and not isinstance(val, str):
        b = _to_bool_int(val)
        return "" if b is None else ("True" if b else "False")
    return val

def reconcile_curated_vouchers_csv_with_db(
    conn: sqlite3.Connection,
    csv_path: Union[str, Path],
    taxa: Optional[Iterable[str]] = None,
    taxon_resolver: Optional[Any] = None,
) -> int:
    """
    Append canonical SQLite `vouchers` rows that are missing from a
    `curated_vouchers.csv` snapshot (matched on `catalogNumber`).

    Taxon-restricted harvests (`--taxon` / `--focal-taxon`) UPSERT into `vouchers`
    but are intentionally never promoted to the root CSV snapshot, so the snapshot
    can silently lack requested taxa. When `taxa` is given, only vouchers whose
    resolved taxon is in `taxa` are restored; `taxon_resolver(scientificName,
    species_raw, species_standardized)` resolves each row (defaulting to the
    stored `species_standardized`).

    Existing CSV rows are preserved byte-for-byte; appended rows are projected onto
    the existing CSV header (or the DwC voucher schema when the CSV does not exist)
    and written atomically (`.tmp` + `fsync` + `replace`). Returns the number of
    appended rows.
    """
    target = Path(csv_path)
    db_cols = get_table_columns(conn, "vouchers")
    if "catalogNumber" not in db_cols:
        return 0

    taxa_set = {str(t).strip() for t in (taxa or []) if t is not None and str(t).strip()}
    restrict_taxa = taxa is not None
    if restrict_taxa and not taxa_set:
        return 0

    original_text = ""
    header: List[str] = []
    existing_cats: Set[str] = set()
    if target.is_file() and target.stat().st_size > 0:
        with open(target, "r", encoding="utf-8", newline="") as fh:
            original_text = fh.read()
        reader = csv.DictReader(io.StringIO(original_text, newline=""))
        header = list(reader.fieldnames or [])
        if header and "catalogNumber" not in header:
            raise ValueError(f"{target} is missing the required 'catalogNumber' column.")
        for row in reader:
            cat = str(row.get("catalogNumber") or "").strip()
            if cat:
                existing_cats.add(cat)
    if not header:
        original_text = ""
        header = [c for c in DWC_VOUCHER_COLUMNS if c in db_cols]

    taxon_cols = [c for c in ("scientificName", "species_raw", "species_standardized") if c in db_cols]
    select_cols: List[str] = []
    for c in ["catalogNumber"] + [h for h in header if h in db_cols] + taxon_cols:
        if c not in select_cols:
            select_cols.append(c)

    cur = conn.execute(f"SELECT {', '.join(select_cols)} FROM vouchers ORDER BY catalogNumber;")
    new_rows: List[List[Any]] = []
    for rec in cur.fetchall():
        rec_map = dict(zip(select_cols, tuple(rec)))
        cat = str(rec_map.get("catalogNumber") or "").strip()
        if not cat or cat in existing_cats:
            continue
        if restrict_taxa:
            if taxon_resolver is not None:
                resolved = taxon_resolver(
                    rec_map.get("scientificName"),
                    rec_map.get("species_raw"),
                    rec_map.get("species_standardized"),
                )
            else:
                resolved = _to_str(rec_map.get("species_standardized"))
            if not resolved or str(resolved).strip() not in taxa_set:
                continue
        existing_cats.add(cat)
        new_rows.append([_format_voucher_csv_value(c, rec_map.get(c)) for c in header])

    if not new_rows:
        return 0

    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_name = tempfile.mkstemp(dir=str(target.parent), prefix=f".{target.stem}_", suffix=".csv.tmp")
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(tmp_fd, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, lineterminator="\n")
            if original_text:
                f.write(original_text)
                if not original_text.endswith(("\n", "\r")):
                    f.write("\n")
            else:
                writer.writerow(header)
            writer.writerows(new_rows)
            f.flush()
            os.fsync(f.fileno())
        tmp_path.replace(target)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
        raise
    return len(new_rows)

def get_completed_qc_vouchers(
    conn: sqlite3.Connection,
    catalog_filter: Optional[Iterable[str]] = None,
) -> Set[str]:
    """Query completed vouchers in O(1) via indexed `qc_status IN ('PASS', 'FAILED_QC')`."""
    cur = conn.execute(
        "SELECT catalogNumber FROM vouchers WHERE qc_status IN ('PASS', 'FAILED_QC');"
    )
    done = {str(r[0]).strip() for r in cur.fetchall() if r[0]}
    if catalog_filter is not None:
        allowed = {str(c).strip() for c in catalog_filter}
        done &= allowed
    return done

__all__ = (
    "upsert_taxon_synonyms", "load_taxon_synonyms", "_infer_species_standardized_fallback",
    "upsert_vouchers_records", "_CSV_BOOL_VOUCHER_COLUMNS", "_format_voucher_csv_value",
    "reconcile_curated_vouchers_csv_with_db", "get_completed_qc_vouchers",
)
