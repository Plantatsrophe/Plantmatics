"""
Step 02 contour file readers and atomic per-voucher SQLite segmentation transaction persistence.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sqlite3
from typing import Any, List, Mapping, Optional, Sequence, Tuple, Union

from scripts.core.db_column_constants import PROJECT_ROOT, _to_float
from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated
from scripts.core.db_step02_leaf_rows import (
    _EXTRACTED_LEAVES_INSERT_COLS,
    _build_step02_leaf_rows,
    _build_step02_plant_and_repro_rows,
    _resolve_voucher_calibration_and_qc,
)

def read_contour_csv_xy(contour_path: Union[str, Path]) -> Optional[List[List[float]]]:
    """
    Read a 120-vertex contour CSV from disk and return compact `[[x, y], ...]`
    coordinates (`x, y` pixel contour vertices so Step 03
    metric scale calibration `coo / resolved_ppm` and `area_px` remain exact).
    """
    p = Path(contour_path)
    if not p.is_file() or p.stat().st_size == 0:
        return None
    try:
        with open(p, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            fnames = set(reader.fieldnames or [])
            if {"x", "y"}.issubset(fnames):
                xk, yk = "x", "y"

            else:
                return None
            pts: List[List[float]] = []
            for row in reader:
                fx = _to_float(row.get(xk))
                fy = _to_float(row.get(yk))
                if fx is not None and fy is not None:
                    pts.append([round(fx, 4), round(fy, 4)])
            return pts if len(pts) >= 3 else None
    except OSError:
        return None

def _resolve_contour_file(
    raw_contour_path: Optional[str],
    output_dir: Optional[Path] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> Optional[Path]:
    if not raw_contour_path:
        return None
    cp = Path(raw_contour_path)
    if cp.is_file():
        return cp
    fname = cp.name
    candidates: List[Path] = []
    if output_dir is not None:
        candidates.extend([
            output_dir / "contours" / fname,
            output_dir / "data" / "contours" / fname,
            output_dir.parent / raw_contour_path,
        ])
    candidates.extend([
        workspace_root / raw_contour_path,
        workspace_root / "data" / "contours" / fname,
    ])
    for cand in candidates:
        if cand.is_file():
            return cand
    return None

def record_voucher_segmentation_transaction(
    conn: sqlite3.Connection,
    catalog_number: str,
    extracted_records: Sequence[Mapping[str, Any]],
    failure_record: Optional[Mapping[str, Any]] = None,
    repro_record: Optional[Mapping[str, Any]] = None,
    pheno_record: Optional[Mapping[str, Any]] = None,
    repro_by_pid: Optional[Mapping[int, Mapping[str, Any]]] = None,
    pheno_by_pid: Optional[Mapping[int, Mapping[str, Any]]] = None,
    crowns_by_pid: Optional[Mapping[int, Tuple[float, float]]] = None,
    rosette_centroids_by_pid: Optional[Mapping[int, Tuple[float, float]]] = None,
    by_plant_petiole: Optional[Mapping[int, Mapping[str, Any]]] = None,
    contour_coords_by_leaf: Optional[Mapping[int, Sequence[Sequence[float]]]] = None,
    output_dir: Optional[Path] = None,
    image_path: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
    commit: bool = True,
) -> None:
    """
    Persist a single voucher's Step 02 segmentation results in one atomic SQLite
    transaction across `vouchers`, `plant_individuals`, `extracted_leaves`, and
    `voucher_reproductive_phenology`.
    """
    cat = str(catalog_number).strip()
    if not cat:
        return

    repro = dict(repro_record or {})
    pheno = dict(pheno_record or {})
    crowns = dict(crowns_by_pid or {})
    centroids = dict(rosette_centroids_by_pid or {})
    petiole_by_pid = dict(by_plant_petiole or {})
    coords_map = dict(contour_coords_by_leaf or {})

    ppm_val, is_cal_val, qc_status, fail_reason, fail_details, eff_img_path = (
        _resolve_voucher_calibration_and_qc(extracted_records, failure_record, image_path)
    )
    out_base = Path(output_dir) if output_dir is not None else (workspace_root / "data")

    plant_rows, repro_pheno_rows = _build_step02_plant_and_repro_rows(
        cat=cat,
        extracted_records=extracted_records,
        ppm_val=ppm_val,
        repro=repro,
        pheno=pheno,
        repro_by_pid=repro_by_pid,
        pheno_by_pid=pheno_by_pid,
        crowns=crowns,
        centroids=centroids,
        petiole_by_pid=petiole_by_pid,
    )
    leaf_rows = _build_step02_leaf_rows(
        cat=cat,
        extracted_records=extracted_records,
        ppm_val=ppm_val,
        coords_map=coords_map,
        out_base=out_base,
        workspace_root=workspace_root,
        resolve_contour_fn=_resolve_contour_file,
        read_contour_fn=read_contour_csv_xy,
    )

    if extracted_records and not leaf_rows:
        # All extracted records lacked readable contours
        qc_status = "FAILED_QC"
        fail_reason = fail_reason or "MISSING_CONTOUR_COORDINATES"

    _ensure_dynamic_table_migrated(conn, "plant_individuals")
    _ensure_dynamic_table_migrated(conn, "extracted_leaves")
    _ensure_dynamic_table_migrated(conn, "voucher_reproductive_phenology")

    def _apply_sql() -> None:
        conn.execute(
            """
            INSERT INTO vouchers (catalogNumber, image_path, pixels_per_mm, is_calibrated_mm, qc_status, failure_reason, failure_details)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(catalogNumber) DO UPDATE SET
                image_path = COALESCE(vouchers.image_path, excluded.image_path),
                pixels_per_mm = COALESCE(excluded.pixels_per_mm, vouchers.pixels_per_mm),
                is_calibrated_mm = COALESCE(excluded.is_calibrated_mm, vouchers.is_calibrated_mm),
                qc_status = excluded.qc_status,
                failure_reason = excluded.failure_reason,
                failure_details = excluded.failure_details;
            """,
            (cat, eff_img_path, ppm_val, is_cal_val, qc_status, fail_reason, fail_details),
        )
        conn.execute("DELETE FROM extracted_leaves WHERE catalogNumber = ?;", (cat,))
        conn.execute("DELETE FROM voucher_reproductive_phenology WHERE catalogNumber = ?;", (cat,))
        conn.execute("DELETE FROM plant_individuals WHERE catalogNumber = ?;", (cat,))

        conn.executemany(
            """
            INSERT INTO plant_individuals (
                catalogNumber, plant_individual_id,
                rosette_centroid_x, rosette_centroid_y, crown_x, crown_y,
                plant_height_px, plant_height_mm, plant_height_cm, is_bent_stem,
                petiole_length_px, petiole_length_mm, petiole_width_px, petiole_width_mm,
                blade_length_px, blade_length_mm,
                petiole_to_blade_ratio, scape_to_basal_leaf_ratio,
                rosette_patch_path, rosette_mask_path
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """,
            plant_rows,
        )
        if repro_pheno_rows:
            conn.executemany(
                """
                INSERT INTO voucher_reproductive_phenology (
                    catalogNumber, plant_individual_id,
                    capitula_count, involucre_height_px, involucre_width_px,
                    involucre_height_mm, involucre_width_mm, AR_inv,
                    total_capitula, n_anthesis, n_bud, n_fruiting, n_senescent,
                    voucher_phenological_state, is_verified_anthesis,
                    reproductive_stage_confidence, dominant_pappus_ratio,
                    corymb_head_count_est, radiate_head_fraction,
                    phyllary_tip_darkening_delta, phyllary_tip_anthocyanin_delta,
                    involucre_basal_cie_l_delta, peduncle_cie_l_sheet_delta, calyculus_flare_index
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                """,
                repro_pheno_rows,
            )
        if leaf_rows:
            el_cols_info = conn.execute("PRAGMA table_info(extracted_leaves);").fetchall()
            insert_cols = list(_EXTRACTED_LEAVES_INSERT_COLS)
            known_set = set(insert_cols)
            extra_vals: List[Any] = []
            for _cid, cname, ctype, notnull, dflt_val, _pk in el_cols_info:
                if cname not in known_set and int(notnull or 0) == 1 and dflt_val is None:
                    insert_cols.append(cname)
                    if "INT" in str(ctype or "").upper():
                        extra_vals.append(0)
                    elif "REAL" in str(ctype or "").upper():
                        extra_vals.append(0.0)
                    else:
                        extra_vals.append("legacy")
            rows_to_insert = [tuple(r) + tuple(extra_vals) for r in leaf_rows] if extra_vals else leaf_rows
            placeholders = ", ".join(["?"] * len(insert_cols))
            conn.executemany(
                f"INSERT INTO extracted_leaves ({', '.join(insert_cols)}) VALUES ({placeholders});",
                rows_to_insert,
            )

    if commit:
        with conn:
            _apply_sql()
    else:
        _apply_sql()

__all__ = (
    "read_contour_csv_xy", "_resolve_contour_file", "record_voucher_segmentation_transaction",
)
