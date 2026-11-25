"""
Database and filesystem purge utilities (`purge_segments_and_runs`, `purge_runs_only`).
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
import sqlite3
from typing import Dict, List, Optional, Union

from scripts.core.db_column_constants import PROJECT_ROOT, logger
from scripts.core.db_conn_pragmas import get_db_connection, get_table_columns, resolve_db_path
from scripts.core.db_schema import init_db_schema

def purge_segments_and_runs(
    workspace_root: Optional[Union[str, Path]] = None,
    db_path: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    purge_runs: bool = True,
    vacuum: bool = True,
) -> Dict[str, int]:
    """
    Purges all Step 02 to 07 segmentation artifacts (masks, contours, ECT matrices,
    cropped patches, non-Step-01 tables, figures, reports, and run directories) and
    clears Step 02 to 07 SQLite tables while preserving Step 01 harvested vouchers,
    raw voucher imagery, taxon synonyms, `taxon_backbone_profiles`,
    `natureserve_subnational_ranges`, `cytological_records`,
    `taxon_cytology_summary`, `data/taxonomy/`, `data/cytology/`, and `environmental_cache`.
    """
    import shutil as _shutil

    ws = Path(workspace_root).resolve() if workspace_root is not None else PROJECT_ROOT
    data_root = Path(output_dir).resolve() if output_dir is not None else (ws / "data")

    target_dirs = [
        data_root / "contours",
        data_root / "masks",
        data_root / "ect",
        data_root / "cropped_patches",
    ]
    if data_root == (ws / "data").resolve():
        target_dirs.extend([
            ws / "outputs" / "figures",
            ws / "outputs" / "reports",
        ])
        if purge_runs:
            target_dirs.append(ws / "runs")

    removed_files = 0
    for d in target_dirs:
        if d.exists() or d.is_symlink():
            if d.is_symlink() or d.is_file():
                try:
                    d.unlink(missing_ok=True)
                    removed_files += 1
                except OSError:
                    pass
            else:
                removed_files += sum(1 for p in d.rglob("*") if p.is_file() or p.is_symlink())
                _shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True, exist_ok=True)


    for sub in ("capitula", "rosettes"):
        (data_root / "cropped_patches" / sub).mkdir(parents=True, exist_ok=True)

    tables_dir = data_root / "tables"
    if tables_dir.exists():
        for f in tables_dir.iterdir():
            if f.is_file() and f.name != "curated_vouchers.csv":
                try:
                    f.unlink(missing_ok=True)
                    removed_files += 1
                except OSError:
                    pass

    resolved_db = resolve_db_path(db_path=db_path, workspace_root=ws)
    purged_rows = 0
    if resolved_db != ":memory:" and Path(resolved_db).exists():
        conn = get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=ws)
        try:
            conn.execute("PRAGMA foreign_keys = OFF;")
            purge_tables = (
                "taxon_complexes",
                "genetic_lineage_assignments",
                "curation_triage",
                "spatial_rf_predictions",
                "otu_morphometrics",
                "leaf_efa_harmonics",
                "run_active_cohort",
                "pipeline_runs",
                "voucher_reproductive_phenology",
                "extracted_leaves",
                "plant_individuals",
            )
            conn.execute("DROP VIEW IF EXISTS v_integrative_triage_master;")
            for tbl in purge_tables:
                try:
                    cnt_row = conn.execute(f"SELECT COUNT(*) FROM {tbl};").fetchone()
                    if cnt_row and cnt_row[0]:
                        purged_rows += int(cnt_row[0])
                except sqlite3.Error:
                    pass
                try:
                    conn.execute(f"DROP TABLE IF EXISTS {tbl};")
                except sqlite3.Error:
                    pass
            try:
                conn.execute(
                    "UPDATE vouchers SET qc_status = NULL, failure_reason = NULL, failure_details = NULL;"
                )
            except sqlite3.Error:
                pass
            conn.execute("PRAGMA user_version = 0;")
            conn.commit()
            init_db_schema(conn)
            if vacuum:
                try:
                    conn.execute("VACUUM;")
                    conn.commit()
                except sqlite3.Error:
                    pass
        finally:
            conn.close()

    summary = {"removed_files": int(removed_files), "purged_db_rows": int(purged_rows)}
    logger.info("Purged Step 02-07 segments and run artifacts: %s", summary)
    return summary

def purge_runs_only(
    workspace_root: Optional[Union[str, Path]] = None,
    db_path: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    refresh_env_buffers: bool = True,
    vacuum: bool = True,
) -> Dict[str, int]:
    """
    Purges ONLY Step 03 to 07 run outputs (`runs/*`, `outputs/figures/*`,
    `outputs/reports/*`, `outputs/rasters/*`, non-Step-01/02 tables in `data/tables/`,
    and Tier B SQLite run/morphometrics/ENM/triage tables) while 100% preserving
    Step 01 harvested vouchers, raw voucher images, and Step 02 segmentation assets
    (`data/masks`, `data/contours`, `data/ect`, `data/cropped_patches`,
    `extracted_leaf_manifest.csv`, `voucher_reproductive_metrics.csv`,
    `voucher_phenological_states.csv`, `failed_qc_vouchers.csv`, `data/taxonomy/`, `data/cytology/`,
    and SQLite `vouchers`, `plant_individuals`, `extracted_leaves`,
    `voucher_reproductive_phenology`, `taxon_backbone_profiles`,
    `natureserve_subnational_ranges`, `cytological_records`, and `taxon_cytology_summary`).

    When `refresh_env_buffers=True`, resets cached `min_bedrock_depth_buffer_cm` to NULL/NA
    in `environmental_cache` and `fine_grain_edaphic_cache.csv` so Step 05 recomputes true
    minimum bedrock depth within coordinate-uncertainty buffers from local rasters without
    re-querying network APIs.
    """
    import shutil as _shutil

    ws = Path(workspace_root).resolve() if workspace_root is not None else PROJECT_ROOT
    data_root = Path(output_dir).resolve() if output_dir is not None else (ws / "data")

    target_dirs: List[Path] = []
    if data_root == (ws / "data").resolve():
        target_dirs.extend([
            ws / "runs",
            ws / "outputs" / "figures",
            ws / "outputs" / "reports",
            ws / "outputs" / "rasters",
        ])

    removed_files = 0
    for d in target_dirs:
        if d.exists() or d.is_symlink():
            if d.is_symlink() or d.is_file():
                try:
                    d.unlink(missing_ok=True)
                    removed_files += 1
                except OSError:
                    pass
            else:
                removed_files += sum(1 for p in d.rglob("*") if p.is_file() or p.is_symlink())
                _shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True, exist_ok=True)

    preserved_step01_02_tables = {
        "curated_vouchers.csv",
        "extracted_leaf_manifest.csv",
        "failed_qc_vouchers.csv",
        "voucher_reproductive_metrics.csv",
        "voucher_phenological_states.csv",
    }
    tables_dir = data_root / "tables"
    if tables_dir.exists():
        for f in tables_dir.iterdir():
            if f.is_file() and f.name not in preserved_step01_02_tables:
                try:
                    f.unlink(missing_ok=True)
                    removed_files += 1
                except OSError:
                    pass

    if refresh_env_buffers:
        env_cache_csv = data_root / "environmental" / "fine_grain_edaphic_cache.csv"
        if env_cache_csv.is_file() and env_cache_csv.stat().st_size > 0:
            try:
                with open(env_cache_csv, "r", encoding="utf-8", newline="") as fh:
                    reader = csv.DictReader(fh)
                    fieldnames = list(reader.fieldnames or [])
                    rows_c = list(reader)
                if "min_bedrock_depth_buffer_cm" in fieldnames and rows_c:
                    for rc in rows_c:
                        rc["min_bedrock_depth_buffer_cm"] = ""
                    tmp_c = env_cache_csv.with_suffix(".csv.tmp")
                    with open(tmp_c, "w", encoding="utf-8", newline="") as out_fh:
                        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, lineterminator="\n")
                        writer.writeheader()
                        writer.writerows(rows_c)
                        out_fh.flush()
                        os.fsync(out_fh.fileno())
                    tmp_c.replace(env_cache_csv)
            except Exception:
                pass

    resolved_db = resolve_db_path(db_path=db_path, workspace_root=ws)
    purged_rows = 0
    if resolved_db != ":memory:" and Path(resolved_db).exists():
        conn = get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=ws)
        try:
            conn.execute("PRAGMA foreign_keys = OFF;")
            tier_b_run_tables = (
                "taxon_complexes",
                "genetic_lineage_assignments",
                "curation_triage",
                "spatial_rf_predictions",
                "otu_morphometrics",
                "leaf_efa_harmonics",
                "run_active_cohort",
                "pipeline_runs",
            )
            conn.execute("DROP VIEW IF EXISTS v_integrative_triage_master;")
            for tbl in tier_b_run_tables:
                try:
                    cnt_row = conn.execute(f"SELECT COUNT(*) FROM {tbl};").fetchone()
                    if cnt_row and cnt_row[0]:
                        purged_rows += int(cnt_row[0])
                except sqlite3.Error:
                    pass
                try:
                    conn.execute(f"DROP TABLE IF EXISTS {tbl};")
                except sqlite3.Error:
                    pass
            if refresh_env_buffers:
                try:
                    env_cols = get_table_columns(conn, "environmental_cache")
                    if "min_bedrock_depth_buffer_cm" in env_cols:
                        conn.execute("UPDATE environmental_cache SET min_bedrock_depth_buffer_cm = NULL;")
                except sqlite3.Error:
                    pass
            conn.commit()
            init_db_schema(conn)
            if vacuum:
                try:
                    conn.execute("VACUUM;")
                    conn.commit()
                except sqlite3.Error:
                    pass
        finally:
            conn.close()

    summary = {"removed_files": int(removed_files), "purged_db_rows": int(purged_rows)}
    logger.info("Purged Step 03-07 run outputs only (Step 01 & Step 02 segmentation preserved): %s", summary)
    return summary

__all__ = (
    "purge_segments_and_runs", "purge_runs_only",
)
