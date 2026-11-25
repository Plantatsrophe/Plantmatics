"""
Canonical Tier-A Step 02 SQLite table rebuild utilities stripping legacy columns and CHECK constraints.
"""

from __future__ import annotations

import sqlite3
from typing import Set, Tuple

_CANONICAL_PLANT_INDIVIDUALS_COLS: Tuple[str, ...] = (
    "catalogNumber", "plant_individual_id", "rosette_centroid_x", "rosette_centroid_y", "crown_x",
    "crown_y", "plant_height_px", "plant_height_mm", "plant_height_cm", "is_bent_stem",
    "petiole_length_px", "petiole_length_mm", "petiole_width_px", "petiole_width_mm",
    "blade_length_px", "blade_length_mm", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio",
    "rosette_patch_path", "rosette_mask_path",
)

_CANONICAL_EXTRACTED_LEAVES_COLS: Tuple[str, ...] = (
    "catalogNumber", "plant_individual_id", "leaf_id", "shape_id", "solidity", "circularity",
    "extent", "aspect_ratio", "area_px", "area_mm2", "midrib_angle_deg", "petiole_length_px",
    "petiole_length_mm", "petiole_width_px", "petiole_width_mm", "blade_length_px",
    "blade_length_mm", "blade_width_px", "blade_width_mm", "petiole_to_blade_ratio",
    "pixels_per_mm", "is_calibrated_mm", "scale_source", "mask_path", "mask_is_upright",
    "crop_path", "contour_path", "ect_path", "contour_n_pts", "contour_xy_json", "leaf_zone",
    "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox", "archetype_score",
    "bilateral_dice", "si_a", "lm3_gates_pass", "n_lamina_components", "detector_class_id",
    "is_crop_clipped", "hole_area_frac", "enclosed_void_frac", "bilateral_lobe_symmetry",
    "efa_qc_pass", "qc_flag_reason", "blade_cie_l_norm", "blade_cie_l_sd",
    "blade_cie_l_sheet_delta", "blade_cie_a_delta", "blade_cie_b_delta", "bbox_x1", "bbox_y1",
    "bbox_x2", "bbox_y2", "base_angle_deg", "apex_angle_deg", "lamina_curvature",
    "lamina_trace_length_px", "basal_cordate_index", "margin_serration_index",
    "margin_tooth_acuteness", "terminal_lobe_area_frac", "blade_cie_l_q85_sheet_delta",
    "basal_zone_cie_l_delta", "is_abaxial_pressed_leaf", "abaxial_cie_l_sheet_delta",
)

_CANONICAL_VOUCHER_REPRO_COLS: Tuple[str, ...] = (
    "catalogNumber", "plant_individual_id", "capitula_count", "involucre_height_px",
    "involucre_width_px", "involucre_height_mm", "involucre_width_mm", "AR_inv", "total_capitula",
    "n_anthesis", "n_bud", "n_fruiting", "n_senescent", "voucher_phenological_state",
    "is_verified_anthesis", "reproductive_stage_confidence", "dominant_pappus_ratio",
    "corymb_head_count_est", "radiate_head_fraction",
)

def _table_needs_canonical_rebuild(
    conn: sqlite3.Connection,
    table_name: str,
    canonical_cols: Tuple[str, ...],
) -> bool:
    cols_now = {
        str(r[1]) for r in conn.execute(f"PRAGMA table_info({table_name});").fetchall()
    }
    if not cols_now:
        return False
    ddl_row = conn.execute(
        "SELECT sql FROM sqlite_master WHERE type = 'table' AND name = ?;",
        (table_name,),
    ).fetchone()
    ddl_sql = str(ddl_row[0] or "") if ddl_row else ""
    return cols_now != set(canonical_cols) or ("CHECK" in ddl_sql.upper())

def _tier_a_step02_needs_rebuild(conn: sqlite3.Connection) -> bool:
    return (
        _table_needs_canonical_rebuild(conn, "plant_individuals", _CANONICAL_PLANT_INDIVIDUALS_COLS)
        or _table_needs_canonical_rebuild(conn, "extracted_leaves", _CANONICAL_EXTRACTED_LEAVES_COLS)
        or _table_needs_canonical_rebuild(conn, "voucher_reproductive_phenology", _CANONICAL_VOUCHER_REPRO_COLS)
    )

def _rebuild_tier_a_step02_tables_if_needed(
    conn: sqlite3.Connection,
    existing_tables: Set[str],
) -> None:
    """
    Rebuild Tier-A Step 02 tables (`plant_individuals`, `extracted_leaves`,
    `voucher_reproductive_phenology`) whenever an older on-disk SQLite database
    retains retired columns or obsolete `CHECK` constraints that `ALTER TABLE ADD COLUMN`
    and `CREATE TABLE IF NOT EXISTS` cannot strip. Uses individual `conn.execute()`
    statements rather than `executescript()` so active savepoints remain intact.
    """
    specs = (
        (
            "plant_individuals",
            _CANONICAL_PLANT_INDIVIDUALS_COLS,
            """
            CREATE TABLE _plant_individuals_v20261121 (
                catalogNumber TEXT NOT NULL,
                plant_individual_id INTEGER NOT NULL DEFAULT 0,
                rosette_centroid_x REAL,
                rosette_centroid_y REAL,
                crown_x REAL,
                crown_y REAL,
                plant_height_px REAL,
                plant_height_mm REAL,
                plant_height_cm REAL,
                is_bent_stem INTEGER,
                petiole_length_px REAL,
                petiole_length_mm REAL,
                petiole_width_px REAL,
                petiole_width_mm REAL,
                blade_length_px REAL,
                blade_length_mm REAL,
                petiole_to_blade_ratio REAL,
                scape_to_basal_leaf_ratio REAL,
                rosette_patch_path TEXT,
                rosette_mask_path TEXT,
                PRIMARY KEY (catalogNumber, plant_individual_id),
                FOREIGN KEY (catalogNumber) REFERENCES vouchers(catalogNumber) ON DELETE CASCADE
            );
            """,
            (),
        ),
        (
            "extracted_leaves",
            _CANONICAL_EXTRACTED_LEAVES_COLS,
            """
            CREATE TABLE _extracted_leaves_v20261121 (
                catalogNumber TEXT NOT NULL,
                plant_individual_id INTEGER NOT NULL DEFAULT 0,
                leaf_id INTEGER NOT NULL,
                shape_id TEXT UNIQUE,
                solidity REAL,
                circularity REAL,
                extent REAL,
                aspect_ratio REAL,
                area_px REAL,
                area_mm2 REAL,
                midrib_angle_deg REAL,
                petiole_length_px REAL,
                petiole_length_mm REAL,
                petiole_width_px REAL,
                petiole_width_mm REAL,
                blade_length_px REAL,
                blade_length_mm REAL,
                blade_width_px REAL,
                blade_width_mm REAL,
                petiole_to_blade_ratio REAL,
                pixels_per_mm REAL,
                is_calibrated_mm INTEGER,
                scale_source TEXT,
                mask_path TEXT,
                mask_is_upright INTEGER DEFAULT 1,
                crop_path TEXT,
                contour_path TEXT,
                ect_path TEXT,
                contour_n_pts INTEGER DEFAULT 120,
                contour_xy_json TEXT NOT NULL,
                leaf_zone TEXT DEFAULT 'basal',
                leaf_zone_confidence REAL,
                dist_from_crown_norm REAL,
                cy_rel_plant_bbox REAL,
                archetype_score REAL,
                bilateral_dice REAL,
                si_a REAL,
                lm3_gates_pass INTEGER,
                n_lamina_components INTEGER DEFAULT 1,
                detector_class_id INTEGER DEFAULT 0,
                is_crop_clipped INTEGER DEFAULT 0,
                hole_area_frac REAL DEFAULT 0.0,
                enclosed_void_frac REAL DEFAULT 0.0,
                bilateral_lobe_symmetry INTEGER DEFAULT 0,
                efa_qc_pass INTEGER DEFAULT 1,
                qc_flag_reason TEXT DEFAULT 'PASS',
                blade_cie_l_norm REAL,
                blade_cie_l_sd REAL,
                blade_cie_l_sheet_delta REAL,
                blade_cie_a_delta REAL,
                blade_cie_b_delta REAL,
                bbox_x1 INTEGER,
                bbox_y1 INTEGER,
                bbox_x2 INTEGER,
                bbox_y2 INTEGER,
                base_angle_deg REAL,
                apex_angle_deg REAL,
                lamina_curvature REAL,
                lamina_trace_length_px REAL,
                basal_cordate_index REAL,
                margin_serration_index REAL,
                margin_tooth_acuteness REAL,
                terminal_lobe_area_frac REAL,
                blade_cie_l_q85_sheet_delta REAL,
                basal_zone_cie_l_delta REAL,
                is_abaxial_pressed_leaf INTEGER DEFAULT 0,
                abaxial_cie_l_sheet_delta REAL,
                PRIMARY KEY (catalogNumber, plant_individual_id, leaf_id),
                FOREIGN KEY (catalogNumber, plant_individual_id)
                    REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
            );
            """,
            ("CREATE INDEX IF NOT EXISTS idx_extracted_leaves_cat ON extracted_leaves(catalogNumber);",),
        ),
        (
            "voucher_reproductive_phenology",
            _CANONICAL_VOUCHER_REPRO_COLS,
            """
            CREATE TABLE _voucher_reproductive_phenology_v20261121 (
                catalogNumber TEXT NOT NULL,
                plant_individual_id INTEGER NOT NULL DEFAULT 0,
                capitula_count INTEGER,
                involucre_height_px REAL,
                involucre_width_px REAL,
                involucre_height_mm REAL,
                involucre_width_mm REAL,
                AR_inv REAL,
                total_capitula INTEGER,
                n_anthesis INTEGER,
                n_bud INTEGER,
                n_fruiting INTEGER,
                n_senescent INTEGER,
                voucher_phenological_state TEXT,
                is_verified_anthesis INTEGER,
                reproductive_stage_confidence REAL,
                dominant_pappus_ratio REAL,
                corymb_head_count_est INTEGER DEFAULT 0,
                radiate_head_fraction REAL,
                PRIMARY KEY (catalogNumber, plant_individual_id),
                FOREIGN KEY (catalogNumber, plant_individual_id)
                    REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
            );
            """,
            (),
        ),
    )

    to_rebuild = [
        (tname, ccols, ddl, idxs)
        for tname, ccols, ddl, idxs in specs
        if tname in existing_tables and _table_needs_canonical_rebuild(conn, tname, ccols)
    ]
    if not to_rebuild:
        return

    toggled_fk = False
    if not conn.in_transaction:
        fk_row = conn.execute("PRAGMA foreign_keys;").fetchone()
        if fk_row and int(fk_row[0]) == 1:
            conn.execute("PRAGMA foreign_keys = OFF;")
            toggled_fk = True

    try:
        for tname, ccols, create_tmp_ddl, post_indexes in to_rebuild:
            tmp_name = f"_{tname}_v20261121"
            cols_now = {
                str(r[1]) for r in conn.execute(f"PRAGMA table_info({tname});").fetchall()
            }
            shared_cols = [c for c in ccols if c in cols_now]
            if tname == "extracted_leaves":
                conn.execute("DROP INDEX IF EXISTS idx_extracted_leaves_cat;")
            conn.execute(f"DROP TABLE IF EXISTS {tmp_name};")
            conn.execute(create_tmp_ddl)
            if shared_cols:
                cols_csv = ", ".join(shared_cols)
                select_exprs = ", ".join(
                    "COALESCE(contour_xy_json, '[]') AS contour_xy_json"
                    if c == "contour_xy_json"
                    else c
                    for c in shared_cols
                )
                where_clause = ""
                if not toggled_fk:
                    if tname == "plant_individuals" and "vouchers" in existing_tables:
                        where_clause = " WHERE catalogNumber IN (SELECT catalogNumber FROM vouchers)"
                    elif tname in {"extracted_leaves", "voucher_reproductive_phenology"} and "plant_individuals" in existing_tables:
                        where_clause = (
                            " WHERE EXISTS (SELECT 1 FROM plant_individuals pi "
                            f"WHERE pi.catalogNumber = {tname}.catalogNumber "
                            f"AND pi.plant_individual_id = COALESCE({tname}.plant_individual_id, 0))"
                        )
                conn.execute(
                    f"INSERT OR REPLACE INTO {tmp_name} ({cols_csv}) SELECT {select_exprs} FROM {tname}{where_clause};"
                )
            conn.execute(f"DROP TABLE {tname};")
            conn.execute(f"ALTER TABLE {tmp_name} RENAME TO {tname};")
            for idx_sql in post_indexes:
                conn.execute(idx_sql)
    finally:
        if toggled_fk:
            conn.execute("PRAGMA foreign_keys = ON;")

__all__ = (
    "_CANONICAL_PLANT_INDIVIDUALS_COLS", "_CANONICAL_EXTRACTED_LEAVES_COLS",
    "_CANONICAL_VOUCHER_REPRO_COLS", "_table_needs_canonical_rebuild",
    "_tier_a_step02_needs_rebuild", "_rebuild_tier_a_step02_tables_if_needed",
)
