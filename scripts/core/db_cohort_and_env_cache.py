"""
Run-scoped active cohort population, triage misidentification queries, and R*Tree environmental cache.
"""

from __future__ import annotations

from pathlib import Path
import re
import sqlite3
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from scripts.core.db_cohort_consolidation import consolidate_false_multi_plant_splits_db
from scripts.core.db_column_constants import (
    ENVIRONMENTAL_BOOL_COLS,
    ENVIRONMENTAL_CACHE_COLUMNS,
    ENVIRONMENTAL_TEXT_COLS,
    PROJECT_ROOT,
    _is_missing,
    _to_bool_int,
    _to_float,
    _to_int,
    _to_str,
    format_lat_lon_key,
)
from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated, upsert_rows_generic
from scripts.core.db_step02_loader import ensure_pipeline_run
from scripts.core.db_vouchers_synonyms import upsert_vouchers_records

def populate_run_active_cohort(
    conn: sqlite3.Connection,
    run_id: str,
    voucher_rows: Sequence[Mapping[str, Any]],
    clear_existing: bool = True,
) -> int:
    """
    Populate `run_active_cohort` for `run_id`, expanding every retained voucher
    across all of its `(catalogNumber, plant_individual_id)` rows in `plant_individuals`
    (or `plant_individual_id = 0` for unsegmented/historical sheets).
    """
    rid = str(run_id).strip() or "canonical"
    ensure_pipeline_run(conn, run_id=rid)
    _tier_b_tables = (
        "leaf_efa_harmonics",
        "otu_morphometrics",
        "spatial_rf_predictions",
        "curation_triage",
    )
    if not voucher_rows:
        if clear_existing:
            with conn:
                conn.execute("DELETE FROM run_active_cohort WHERE run_id = ?;", (rid,))
                for _tbl in _tier_b_tables:
                    conn.execute(f"DELETE FROM {_tbl} WHERE run_id = ?;", (rid,))
        return 0

    _ensure_dynamic_table_migrated(conn, "run_active_cohort")
    upsert_vouchers_records(conn, voucher_rows)
    consolidate_false_multi_plant_splits_db(conn)

    cats = [str(r.get("catalogNumber", "")).strip() for r in voucher_rows if str(r.get("catalogNumber", "")).strip()]
    pids_by_cat: Dict[str, List[int]] = {}
    chunk_size = 500
    for i in range(0, len(cats), chunk_size):
        chunk = cats[i : i + chunk_size]
        ph = ", ".join("?" for _ in chunk)
        cur = conn.execute(
            f"SELECT catalogNumber, plant_individual_id FROM plant_individuals WHERE catalogNumber IN ({ph});",
            chunk,
        )
        for row in cur.fetchall():
            pids_by_cat.setdefault(str(row[0]).strip(), []).append(int(row[1]))

    cohort_payload: List[Tuple[Any, ...]] = []
    seen_keys: Set[Tuple[str, int]] = set()
    for r in voucher_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        sp_std = _to_str(r.get("species_standardized"))
        det_tier = _to_str(r.get("determiner_tier"))
        type_stat = _to_str(r.get("type_status"))
        is_type_sp = _to_bool_int(r.get("is_type_specimen"))
        if is_type_sp is None:
            is_type_sp = (
                1
                if (
                    type_stat
                    and type_stat.upper()
                    not in {
                        "NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE",
                        "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN",
                    }
                    and not re.search(r"\b(?:not\s+a?\s*type|non[-\s]?type)\b", type_stat, re.IGNORECASE)
                )
                else 0
            )
        is_prom = _to_bool_int(r.get("is_promoted_bronze"), default=0) or 0
        auth_wt = _to_float(r.get("authority_weight"))
        sp_wt = _to_float(r.get("spatial_balance_weight"))
        if sp_wt is None:
            sp_wt = 1.0

        explicit_pid = _to_int(r.get("plant_individual_id"))
        target_pids = [explicit_pid] if explicit_pid is not None else sorted(set(pids_by_cat.get(cat, [0])))
        for pid in target_pids:
            if (cat, pid) in seen_keys:
                continue
            seen_keys.add((cat, pid))
            cohort_payload.append((rid, cat, pid, sp_std, det_tier, type_stat, is_type_sp, is_prom, auth_wt, sp_wt))

    with conn:
        if clear_existing:
            conn.execute("DELETE FROM run_active_cohort WHERE run_id = ?;", (rid,))
        conn.executemany(
            """
            INSERT INTO run_active_cohort (
                run_id, catalogNumber, plant_individual_id,
                species_standardized, determiner_tier,
                type_status, is_type_specimen,
                is_promoted_bronze, authority_weight, spatial_balance_weight
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id, catalogNumber, plant_individual_id) DO UPDATE SET
                species_standardized = COALESCE(excluded.species_standardized, run_active_cohort.species_standardized),
                determiner_tier = COALESCE(excluded.determiner_tier, run_active_cohort.determiner_tier),
                type_status = COALESCE(excluded.type_status, run_active_cohort.type_status),
                is_type_specimen = COALESCE(excluded.is_type_specimen, run_active_cohort.is_type_specimen),
                is_promoted_bronze = COALESCE(excluded.is_promoted_bronze, run_active_cohort.is_promoted_bronze),
                authority_weight = COALESCE(excluded.authority_weight, run_active_cohort.authority_weight),
                spatial_balance_weight = COALESCE(excluded.spatial_balance_weight, run_active_cohort.spatial_balance_weight);
            """,
            cohort_payload,
        )
        if clear_existing:
            for _tbl in _tier_b_tables:
                conn.execute(
                    f"""
                    DELETE FROM {_tbl}
                    WHERE run_id = ?
                      AND NOT EXISTS (
                          SELECT 1 FROM run_active_cohort rac
                          WHERE rac.run_id = {_tbl}.run_id
                            AND rac.catalogNumber = {_tbl}.catalogNumber
                            AND rac.plant_individual_id = {_tbl}.plant_individual_id
                      );
                    """,
                    (rid,),
                )
            conn.execute(
                """
                UPDATE otu_morphometrics
                SET species_standardized = (
                    SELECT rac.species_standardized
                    FROM run_active_cohort rac
                    WHERE rac.run_id = otu_morphometrics.run_id
                      AND rac.catalogNumber = otu_morphometrics.catalogNumber
                      AND rac.plant_individual_id = otu_morphometrics.plant_individual_id
                )
                WHERE run_id = ?
                  AND EXISTS (
                      SELECT 1
                      FROM run_active_cohort rac
                      WHERE rac.run_id = otu_morphometrics.run_id
                        AND rac.catalogNumber = otu_morphometrics.catalogNumber
                        AND rac.plant_individual_id = otu_morphometrics.plant_individual_id
                        AND rac.species_standardized IS NOT NULL
                        AND TRIM(rac.species_standardized) <> ''
                  );
                """,
                (rid,),
            )
    return len(cohort_payload)

def normalize_run_id_token(token: str, workspace_root: Path = PROJECT_ROOT) -> str:
    """
    Normalize a `--triage-source` token (flat ID, nested ID, `runs/`-prefixed path,
    or path to a run dir / `triage_queue.csv`) to its canonical `run_id` string
    while preserving internal `/` separators in nested run identifiers.
    """
    cleaned = str(token or "").strip().replace("\\", "/")
    if not cleaned:
        return ""
    while cleaned.startswith("./"):
        cleaned = cleaned[2:]
    cleaned = cleaned.rstrip("/")
    for suffix in ("/data/tables/triage_queue.csv", "/data"):
        if cleaned.endswith(suffix):
            cleaned = cleaned[: -len(suffix)]
            break
    cleaned = cleaned.rstrip("/")
    if not cleaned:
        return ""

    if Path(cleaned).is_absolute():
        tp = Path(cleaned).resolve()
        runs_root = (Path(workspace_root) / "runs").resolve()
        if runs_root in tp.parents:
            return tp.relative_to(runs_root).as_posix()
        return tp.name

    while cleaned.startswith("./"):
        cleaned = cleaned[2:]
    if cleaned.startswith("runs/"):
        return cleaned[len("runs/"):].strip("/")
    return cleaned.strip("/")

def query_misidentified_vouchers_from_db(
    conn: sqlite3.Connection,
    triage_sources: Optional[Sequence[str]] = None,
    exclude_run_id: Optional[str] = None,
    allowed_subtypes: Optional[Set[str]] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> Tuple[Set[str], int]:
    """
    Query `curation_triage` directly in SQLite to identify vouchers flagged as
    `Misidentification_Reassign` (or outside `allowed_subtypes`) for `--exclude-misidentified`.
    Returns `(flagged_catalog_numbers, matched_run_count)`.
    """
    allowed = allowed_subtypes or {
        "Confirmed_Species",
        "Confirmed_Type_Control",
        "Type_Synonymy_Candidate",
        "Hybrid_Introgressant",
        "Subspecies_Ecotype",
        "Ecophenotypic_Plasticity",
        "Putative_Cryptic_Species",
    }
    tokens = [str(t).strip() for t in (triage_sources or ["ALL"]) if str(t).strip()]
    if not tokens:
        tokens = ["ALL"]

    where_clauses: List[str] = []
    params: List[Any] = []
    eff_exclude = (
        normalize_run_id_token(exclude_run_id, workspace_root=workspace_root)
        if exclude_run_id
        else ""
    )
    if eff_exclude:
        where_clauses.append("run_id != ?")
        params.append(eff_exclude)

    if not any(t.upper() in {"ALL", "*"} for t in tokens):
        run_ids: List[str] = []
        seen_run_ids: Set[str] = set()
        for t in tokens:
            norm_tok = normalize_run_id_token(t, workspace_root=workspace_root)
            if norm_tok and norm_tok not in seen_run_ids:
                seen_run_ids.add(norm_tok)
                run_ids.append(norm_tok)
        if not run_ids:
            return set(), 0
        ph = ", ".join("?" for _ in run_ids)
        where_clauses.append(f"run_id IN ({ph})")
        params.extend(run_ids)

    where_sql = (" WHERE " + " AND ".join(where_clauses)) if where_clauses else ""
    cur_runs = conn.execute(f"SELECT COUNT(DISTINCT run_id) FROM curation_triage{where_sql};", params)
    n_runs = int((cur_runs.fetchone() or [0])[0] or 0)
    if n_runs == 0:
        return set(), 0

    cur = conn.execute(
        f"SELECT catalogNumber, taxonomic_treatment_subtype FROM curation_triage{where_sql};",
        params,
    )
    flagged: Set[str] = set()
    for row in cur.fetchall():
        cat = str(row[0] or "").strip()
        sub = str(row[1] or "").strip()
        if not cat or not sub:
            continue
        if sub == "Misidentification_Reassign" or sub not in allowed:
            flagged.add(cat)
    return flagged, n_runs

def _normalize_semicolon_list_str(val: Any) -> Optional[str]:
    if val is None:
        return None
    if isinstance(val, (list, tuple, set)):
        items = [str(x).strip() for x in val if not _is_missing(x) and str(x).strip()]
        return ";".join(items) if items else None
    return _to_str(val)

def upsert_environmental_cache_records(
    conn: sqlite3.Connection,
    records: Sequence[Mapping[str, Any]],
    demo: bool = False,
) -> int:
    """
    UPSERT fine-grain SSURGO / POLARIS / SOLUS100 / 3DEP / SGMC / 19-variable WorldClim
    rows into `environmental_cache` (automatically updating `environmental_rtree` via
    triggers). Strictly no-ops when `demo=True`.
    """
    if demo or not records:
        return 0
    _ensure_dynamic_table_migrated(conn, "environmental_cache")
    norm_rows: List[Dict[str, Any]] = []
    for r in records:
        lat = _to_float(r.get("latitude") if "latitude" in r else r.get("decimalLatitude"))
        lon = _to_float(r.get("longitude") if "longitude" in r else r.get("decimalLongitude"))
        ll_key = _to_str(r.get("lat_lon_key")) or format_lat_lon_key(lat, lon)
        if not ll_key:
            continue
        d: Dict[str, Any] = {
            "lat_lon_key": ll_key,
            "catalogNumber": _to_str(r.get("catalogNumber")),
            "latitude": lat,
            "longitude": lon,
        }
        for col in ENVIRONMENTAL_CACHE_COLUMNS:
            if col in {"lat_lon_key", "catalogNumber", "latitude", "longitude"}:
                continue
            if col not in r:
                continue
            if col in ENVIRONMENTAL_TEXT_COLS:
                d[col] = _to_str(r.get(col))
            elif col in ENVIRONMENTAL_BOOL_COLS:
                d[col] = _to_bool_int(r.get(col))
            else:
                d[col] = _to_float(r.get(col))
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="environmental_cache",
        rows=norm_rows,
        pk_cols=("lat_lon_key",),
        preserve_existing_on_null=True,
    )

def query_environmental_bbox(
    conn: sqlite3.Connection,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
) -> List[Dict[str, Any]]:
    """Query `environmental_cache` in O(log N) using the `environmental_rtree` spatial index."""
    cur = conn.execute(
        """
        SELECT ec.*
        FROM environmental_cache ec
        INNER JOIN environmental_rtree rt ON ec.rowid = rt.id
        WHERE rt.min_lon >= ? AND rt.max_lon <= ?
          AND rt.min_lat >= ? AND rt.max_lat <= ?;
        """,
        (float(min_lon), float(max_lon), float(min_lat), float(max_lat)),
    )
    out: List[Dict[str, Any]] = []
    for r in cur.fetchall():
        out.append(dict(r))
    return out

__all__ = (
    "populate_run_active_cohort", "normalize_run_id_token", "query_misidentified_vouchers_from_db",
    "_normalize_semicolon_list_str", "upsert_environmental_cache_records",
    "query_environmental_bbox", "_ensure_dynamic_table_migrated",
)
