"""
Step 02 SQLite resume loader (`load_step02_records_from_db`) and run registration.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

from scripts.core.db_column_constants import _to_str
from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated

class _PidAwareDict(dict):
    """Dict keyed by catalogNumber that also exposes per-(catalogNumber, pid) rows via `.by_pid`."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.by_pid: Dict[Tuple[str, int], Dict[str, Any]] = {}

def load_step02_records_from_db(
    conn: sqlite3.Connection,
    catalog_numbers: Optional[Sequence[str]] = None,
    include_by_pid: bool = False,
) -> Any:
    """
    Fetch existing Step 02 leaf, failed-QC, reproductive, and phenology records
    from SQLite for `catalog_numbers` (or all stored vouchers when `None`) when
    resuming Step 02. Both `repro_by_cat` and `pheno_by_cat` expose `.by_pid`
    mappings keyed by `(catalogNumber, plant_individual_id)`; passing
    `include_by_pid=True` returns a 6-tuple including those mappings explicitly.
    """
    extracted_by_cat: Dict[str, List[Dict[str, Any]]] = {}
    failed_by_cat: Dict[str, Dict[str, Any]] = {}
    repro_by_cat: _PidAwareDict = _PidAwareDict()
    pheno_by_cat: _PidAwareDict = _PidAwareDict()
    repro_by_pid = repro_by_cat.by_pid
    pheno_by_pid = pheno_by_cat.by_pid

    if catalog_numbers is None:
        cat_set: Set[str] = set()
        for q_sql in (
            "SELECT DISTINCT catalogNumber FROM extracted_leaves;",
            "SELECT DISTINCT catalogNumber FROM voucher_reproductive_phenology;",
            "SELECT DISTINCT catalogNumber FROM vouchers WHERE qc_status IN ('PASS', 'FAILED_QC');",
        ):
            for r in conn.execute(q_sql).fetchall():
                if r[0]:
                    cat_set.add(str(r[0]).strip())
        cats = sorted(cat_set)
    else:
        cats = [str(c).strip() for c in catalog_numbers if str(c).strip()]

    if not cats:
        if include_by_pid:
            return extracted_by_cat, failed_by_cat, repro_by_cat, pheno_by_cat, repro_by_pid, pheno_by_pid
        return extracted_by_cat, failed_by_cat, repro_by_cat, pheno_by_cat

    _ensure_dynamic_table_migrated(conn, "plant_individuals")
    _ensure_dynamic_table_migrated(conn, "extracted_leaves")
    _ensure_dynamic_table_migrated(conn, "voucher_reproductive_phenology")
    chunk_size = 500
    for i in range(0, len(cats), chunk_size):
        chunk = cats[i : i + chunk_size]
        ph = ", ".join("?" for _ in chunk)

        cur_l = conn.execute(
            f"""
            SELECT
                el.catalogNumber, el.plant_individual_id, el.leaf_id, el.solidity,
                el.circularity, el.extent, el.aspect_ratio,
                el.area_px, el.area_mm2,
                el.contour_path, el.ect_path, el.mask_path, el.mask_is_upright,
                el.midrib_angle_deg, el.pixels_per_mm, el.is_calibrated_mm, el.scale_source,
                el.leaf_zone, el.leaf_zone_confidence, el.dist_from_crown_norm, el.cy_rel_plant_bbox,
                el.archetype_score, el.bilateral_dice, el.si_a, el.lm3_gates_pass, el.n_lamina_components,
                el.detector_class_id, el.is_crop_clipped, el.hole_area_frac, el.enclosed_void_frac,
                el.bilateral_lobe_symmetry, el.efa_qc_pass, el.qc_flag_reason,
                el.blade_cie_l_norm, el.blade_cie_l_sd,
                el.blade_cie_l_sheet_delta, el.blade_cie_a_delta, el.blade_cie_b_delta,
                el.bbox_x1, el.bbox_y1, el.bbox_x2, el.bbox_y2,
                el.base_angle_deg, el.apex_angle_deg, el.lamina_curvature, el.lamina_trace_length_px,
                el.basal_cordate_index, el.margin_serration_index, el.margin_tooth_acuteness, el.terminal_lobe_area_frac,
                el.blade_cie_l_q85_sheet_delta, el.basal_zone_cie_l_delta,
                el.is_abaxial_pressed_leaf, el.abaxial_cie_l_sheet_delta,
                el.petiole_length_px,
                el.blade_length_px,
                el.blade_width_px, el.blade_width_mm,
                el.petiole_width_px,
                el.petiole_width_mm,
                el.petiole_length_mm,
                el.blade_length_mm,
                el.petiole_to_blade_ratio,
                vrp.capitula_count, vrp.involucre_height_px, vrp.involucre_width_px,
                vrp.involucre_height_mm, vrp.involucre_width_mm,
                vrp.AR_inv, vrp.corymb_head_count_est, vrp.radiate_head_fraction
            FROM extracted_leaves el
            LEFT JOIN plant_individuals pi
                ON pi.catalogNumber = el.catalogNumber
               AND pi.plant_individual_id = el.plant_individual_id
            LEFT JOIN voucher_reproductive_phenology vrp
                ON vrp.catalogNumber = el.catalogNumber
               AND vrp.plant_individual_id = el.plant_individual_id
            WHERE el.catalogNumber IN ({ph})
            ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;
            """,
            chunk,
        )
        for row in cur_l.fetchall():
            d = dict(row)
            d["is_calibrated_mm"] = bool(d.get("is_calibrated_mm"))
            raw_up = d.get("mask_is_upright")
            d["mask_is_upright"] = True if raw_up is None else bool(raw_up)
            if d.get("lm3_gates_pass") is not None:
                d["lm3_gates_pass"] = bool(d.get("lm3_gates_pass"))
            if d.get("is_crop_clipped") is not None:
                d["is_crop_clipped"] = bool(d.get("is_crop_clipped"))
            if d.get("bilateral_lobe_symmetry") is not None:
                d["bilateral_lobe_symmetry"] = bool(d.get("bilateral_lobe_symmetry"))
            if d.get("efa_qc_pass") is not None:
                d["efa_qc_pass"] = bool(d.get("efa_qc_pass"))
            if d.get("is_abaxial_pressed_leaf") is not None:
                d["is_abaxial_pressed_leaf"] = bool(d.get("is_abaxial_pressed_leaf"))
            if "cf_source" not in d or d.get("cf_source") is None:
                s_src = str(d.get("scale_source") or "").strip().lower()
                if s_src == "mp_estimated":
                    d["cf_source"] = "predicted_from_megapixels"
                elif s_src == "ruler_calibrated":
                    d["cf_source"] = "measured_from_ruler"
                else:
                    d["cf_source"] = "uncalibrated"
            extracted_by_cat.setdefault(str(d["catalogNumber"]), []).append(d)

        cur_f = conn.execute(
            f"""
            SELECT catalogNumber, image_path, failure_reason, failure_details AS details
            FROM vouchers
            WHERE qc_status = 'FAILED_QC' AND catalogNumber IN ({ph});
            """,
            chunk,
        )
        for row in cur_f.fetchall():
            d = dict(row)
            failed_by_cat[str(d["catalogNumber"])] = d

        cur_rp = conn.execute(
            f"""
            SELECT
                vrp.catalogNumber, vrp.plant_individual_id,
                vrp.capitula_count, vrp.involucre_height_px, vrp.involucre_width_px,
                vrp.involucre_height_mm, vrp.involucre_width_mm,
                vrp.AR_inv, vrp.corymb_head_count_est, vrp.radiate_head_fraction,
                pi.plant_height_cm, pi.scape_to_basal_leaf_ratio, pi.is_bent_stem,
                vrp.total_capitula, vrp.n_anthesis, vrp.n_fruiting,
                vrp.n_bud, vrp.n_senescent, vrp.voucher_phenological_state,
                vrp.is_verified_anthesis, vrp.reproductive_stage_confidence,
                vrp.dominant_pappus_ratio
            FROM voucher_reproductive_phenology vrp
            LEFT JOIN plant_individuals pi
                ON pi.catalogNumber = vrp.catalogNumber
               AND pi.plant_individual_id = vrp.plant_individual_id
            WHERE vrp.catalogNumber IN ({ph})
            ORDER BY vrp.catalogNumber ASC, vrp.plant_individual_id ASC;
            """,
            chunk,
        )
        for row in cur_rp.fetchall():
            d = dict(row)
            c = str(d["catalogNumber"])
            pid_i = int(d.get("plant_individual_id") or 0)
            repro_item = {
                "catalogNumber": c,
                "plant_individual_id": pid_i,
                "capitula_count": d.get("capitula_count"),
                "involucre_height_px": d.get("involucre_height_px"),
                "involucre_width_px": d.get("involucre_width_px"),
                "involucre_height_mm": d.get("involucre_height_mm"),
                "involucre_width_mm": d.get("involucre_width_mm"),
                "AR_inv": d.get("AR_inv"),
                "corymb_head_count_est": int(d.get("corymb_head_count_est") or 0),
                "radiate_head_fraction": d.get("radiate_head_fraction"),
                "plant_height_cm": d.get("plant_height_cm"),
                "scape_to_basal_leaf_ratio": d.get("scape_to_basal_leaf_ratio"),
                "is_bent_stem": bool(d.get("is_bent_stem")),
            }
            repro_by_pid[(c, pid_i)] = repro_item
            if c not in repro_by_cat:
                repro_by_cat[c] = repro_item
            if d.get("voucher_phenological_state") is not None:
                pheno_item = {
                    "catalogNumber": c,
                    "plant_individual_id": pid_i,
                    "total_capitula": d.get("total_capitula"),
                    "n_anthesis": d.get("n_anthesis"),
                    "n_fruiting": d.get("n_fruiting"),
                    "n_bud": d.get("n_bud"),
                    "n_senescent": d.get("n_senescent"),
                    "voucher_phenological_state": d.get("voucher_phenological_state"),
                    "is_verified_anthesis": bool(d.get("is_verified_anthesis")) if d.get("is_verified_anthesis") is not None else None,
                    "reproductive_stage_confidence": d.get("reproductive_stage_confidence"),
                    "dominant_pappus_ratio": d.get("dominant_pappus_ratio"),
                    "corymb_head_count_est": int(d.get("corymb_head_count_est") or 0),
                    "radiate_head_fraction": d.get("radiate_head_fraction"),
                }
                pheno_by_pid[(c, pid_i)] = pheno_item
                if c not in pheno_by_cat:
                    pheno_by_cat[c] = pheno_item

    if include_by_pid:
        return extracted_by_cat, failed_by_cat, repro_by_cat, pheno_by_cat, repro_by_pid, pheno_by_pid
    return extracted_by_cat, failed_by_cat, repro_by_cat, pheno_by_cat

def ensure_pipeline_run(
    conn: sqlite3.Connection,
    run_id: str,
    run_dir: Optional[Union[str, Path]] = None,
    start_timestamp: Optional[str] = None,
    end_timestamp: Optional[str] = None,
    git_commit_sha: Optional[str] = None,
    cli_arguments: Optional[Sequence[str]] = None,
    focal_taxon: Optional[str] = None,
    intraspecific_mode: Optional[bool] = None,
) -> None:
    """Register or update a run record inside `pipeline_runs`."""
    rid = str(run_id).strip() or "canonical"
    rdir_str = str(run_dir) if run_dir is not None else None
    cli_json = json.dumps(list(cli_arguments)) if cli_arguments is not None else None
    now_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    im_ins = 1 if intraspecific_mode else 0
    im_opt = None if intraspecific_mode is None else (1 if intraspecific_mode else 0)

    with conn:
        conn.execute(
            """
            INSERT INTO pipeline_runs (
                run_id, run_dir, start_timestamp, end_timestamp,
                git_commit_sha, cli_arguments_json, focal_taxon, intraspecific_mode
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id) DO UPDATE SET
                run_dir = COALESCE(excluded.run_dir, pipeline_runs.run_dir),
                start_timestamp = COALESCE(pipeline_runs.start_timestamp, excluded.start_timestamp),
                end_timestamp = COALESCE(excluded.end_timestamp, pipeline_runs.end_timestamp),
                git_commit_sha = COALESCE(excluded.git_commit_sha, pipeline_runs.git_commit_sha),
                cli_arguments_json = COALESCE(excluded.cli_arguments_json, pipeline_runs.cli_arguments_json),
                focal_taxon = COALESCE(excluded.focal_taxon, pipeline_runs.focal_taxon),
                intraspecific_mode = CASE WHEN ? IS NOT NULL THEN ? ELSE pipeline_runs.intraspecific_mode END;
            """,
            (
                rid,
                rdir_str,
                start_timestamp or now_iso,
                end_timestamp,
                git_commit_sha,
                cli_json,
                _to_str(focal_taxon),
                im_ins,
                im_opt,
                im_opt,
            ),
        )

__all__ = (
    "_PidAwareDict", "load_step02_records_from_db", "ensure_pipeline_run",
)
