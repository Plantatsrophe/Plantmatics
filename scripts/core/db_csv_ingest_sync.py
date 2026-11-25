"""
Bidirectional bootstrap synchronization of existing CSV manifests and contours into SQLite.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from scripts.core.db_cohort_and_env_cache import (
    populate_run_active_cohort,
    upsert_environmental_cache_records,
)
from scripts.core.db_cohort_consolidation import consolidate_false_multi_plant_splits_db
from scripts.core.db_column_constants import PROJECT_ROOT, _to_float, _to_int, _to_str, logger
from scripts.core.db_conn_pragmas import get_db_connection, upsert_rows_generic
from scripts.core.db_run_tables_ingest import _read_csv_dicts, ingest_existing_run_tables_into_db
from scripts.core.db_step02_contours_tx import (
    read_contour_csv_xy,
    record_voucher_segmentation_transaction,
)
from scripts.core.db_step02_loader import ensure_pipeline_run
from scripts.core.db_vouchers_synonyms import upsert_taxon_synonyms, upsert_vouchers_records

def ingest_existing_csvs_to_db(
    workspace_root: Union[str, Path] = PROJECT_ROOT,
    run_id: str = "canonical",
    db_path: Optional[Union[str, Path]] = None,
    tables_dir: Optional[Union[str, Path]] = None,
    contours_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, int]:
    """
    Bootstrap or synchronize `data/plantmatics.sqlite` (migrating legacy
    `data/packera_pipeline.sqlite` if present) from existing CSV manifests,
    120-vertex contour CSVs (`data/contours/*.csv`), and `fine_grain_edaphic_cache.csv`.
    """
    ws = Path(workspace_root).resolve()
    t_dir = Path(tables_dir).resolve() if tables_dir is not None else (ws / "data" / "tables")
    c_dir = Path(contours_dir).resolve() if contours_dir is not None else (ws / "data" / "contours")
    env_csv = ws / "data" / "environmental" / "fine_grain_edaphic_cache.csv"

    conn = get_db_connection(db_path=db_path, read_only=False, demo=False, init_schema=True, workspace_root=ws)
    counts: Dict[str, int] = {}

    try:
        # 1. Seed taxon_synonyms from harvester if importable
        try:
            from scripts.core.harvester import _DYNAMIC_SYNONYM_TO_CANONICAL, _DYNAMIC_TAXON_AUTHOR_STRINGS

            if _DYNAMIC_SYNONYM_TO_CANONICAL:
                counts["taxon_synonyms"] = upsert_taxon_synonyms(
                    conn, _DYNAMIC_SYNONYM_TO_CANONICAL, _DYNAMIC_TAXON_AUTHOR_STRINGS
                )
        except Exception:
            counts["taxon_synonyms"] = 0

        # 2. Vouchers (`curated_vouchers.csv`)
        vouchers_rows = _read_csv_dicts(t_dir / "curated_vouchers.csv")
        counts["vouchers"] = upsert_vouchers_records(conn, vouchers_rows)

        # 3. Step 02 tables + 120-vertex contours
        failed_rows = _read_csv_dicts(t_dir / "failed_qc_vouchers.csv")
        manifest_path = t_dir / "extracted_leaf_manifest.csv"
        leaf_manifest_rows = _read_csv_dicts(manifest_path)
        repro_rows = _read_csv_dicts(t_dir / "voucher_reproductive_metrics.csv")
        pheno_rows = _read_csv_dicts(t_dir / "voucher_phenological_states.csv")

        repro_by_cat: Dict[str, Dict[str, Any]] = {}
        repro_by_cat_pid: Dict[str, Dict[int, Dict[str, Any]]] = {}
        for r in repro_rows:
            c_key = str(r.get("catalogNumber", "")).strip()
            if not c_key:
                continue
            p_key = _to_int(r.get("plant_individual_id"), default=0) or 0
            repro_by_cat.setdefault(c_key, r)
            repro_by_cat_pid.setdefault(c_key, {})[p_key] = r
        pheno_by_cat: Dict[str, Dict[str, Any]] = {}
        pheno_by_cat_pid: Dict[str, Dict[int, Dict[str, Any]]] = {}
        for r in pheno_rows:
            c_key = str(r.get("catalogNumber", "")).strip()
            if not c_key:
                continue
            p_key = _to_int(r.get("plant_individual_id"), default=0) or 0
            pheno_by_cat.setdefault(c_key, r)
            pheno_by_cat_pid.setdefault(c_key, {})[p_key] = r
        failed_by_cat = {str(r.get("catalogNumber", "")).strip(): r for r in failed_rows if r.get("catalogNumber")}
        leaves_by_cat: Dict[str, List[Dict[str, Any]]] = {}
        for lr in leaf_manifest_rows:
            cat = str(lr.get("catalogNumber", "")).strip()
            if cat:
                leaves_by_cat.setdefault(cat, []).append(lr)

        all_step02_cats = sorted(
            set(leaves_by_cat.keys())
            | set(failed_by_cat.keys())
            | set(repro_by_cat.keys())
            | set(pheno_by_cat.keys())
        )

        # Index contour files by filename for O(1) resolution
        contour_index: Dict[str, Path] = {}
        for scan_cdir in (c_dir, ws / "data" / "contours"):
            if scan_cdir.is_dir():
                with os.scandir(scan_cdir) as it:
                    for entry in it:
                        if entry.is_file() and entry.name.endswith(".csv"):
                            contour_index.setdefault(entry.name, Path(entry.path))

        migrated_leaves = 0
        for cat in all_step02_cats:
            cat_leaves = leaves_by_cat.get(cat, [])
            coords_by_leaf: Dict[int, List[List[float]]] = {}
            valid_cat_leaves: List[Dict[str, Any]] = []
            for lr in cat_leaves:
                lid = _to_int(lr.get("leaf_id"), default=1) or 1
                raw_cp = _to_str(lr.get("contour_path"))
                c_file: Optional[Path] = None
                if raw_cp:
                    cand = Path(raw_cp)
                    if cand.is_file():
                        c_file = cand
                    elif cand.name in contour_index:
                        c_file = contour_index[cand.name]
                if c_file is None:
                    fallback_name = f"{cat}_leaf{lid}.csv"
                    c_file = contour_index.get(fallback_name)
                if c_file is not None:
                    pts = read_contour_csv_xy(c_file)
                    if pts:
                        coords_by_leaf[lid] = pts
                        valid_cat_leaves.append(lr)

            fail_rec = failed_by_cat.get(cat)
            if not valid_cat_leaves and fail_rec is None and not cat_leaves:
                # Voucher only appeared in reproductive/phenology table; ensure stub individual without marking FAILED_QC
                with conn:
                    conn.execute(
                        "INSERT OR IGNORE INTO vouchers (catalogNumber, qc_status) VALUES (?, 'PENDING');",
                        (cat,),
                    )
                    conn.execute(
                        "INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 0);",
                        (cat,),
                    )
                if cat in repro_by_cat or cat in pheno_by_cat:
                    rep = repro_by_cat.get(cat, {})
                    phe = pheno_by_cat.get(cat, {})
                    ar_val = _to_float(rep.get("AR_inv"))
                    upsert_rows_generic(
                        conn,
                        "voucher_reproductive_phenology",
                        [{
                            "catalogNumber": cat,
                            "plant_individual_id": _to_int(rep.get("plant_individual_id"), default=0) or 0,
                            "capitula_count": _to_int(rep.get("capitula_count")),
                            "involucre_height_px": _to_float(rep.get("involucre_height_px")),
                            "involucre_width_px": _to_float(rep.get("involucre_width_px")),
                            "involucre_height_mm": _to_float(rep.get("involucre_height_mm")),
                            "involucre_width_mm": _to_float(rep.get("involucre_width_mm")),
                            "AR_inv": ar_val,
                            "total_capitula": _to_int(phe.get("total_capitula")),
                            "n_anthesis": _to_int(phe.get("n_anthesis")),
                            "n_bud": _to_int(phe.get("n_bud")),
                            "n_fruiting": _to_int(phe.get("n_fruiting")),
                            "voucher_phenological_state": _to_str(phe.get("voucher_phenological_state")),
                            "dominant_pappus_ratio": _to_float(phe.get("dominant_pappus_ratio")),
                            "corymb_head_count_est": _to_int(
                                rep.get("corymb_head_count_est")
                                if rep.get("corymb_head_count_est") is not None
                                else phe.get("corymb_head_count_est"),
                                default=0,
                            ) or 0,
                            "radiate_head_fraction": _to_float(
                                rep.get("radiate_head_fraction")
                                if rep.get("radiate_head_fraction") is not None
                                else phe.get("radiate_head_fraction")
                            ),
                        }],
                        pk_cols=("catalogNumber", "plant_individual_id"),
                    )
                continue

            record_voucher_segmentation_transaction(
                conn=conn,
                catalog_number=cat,
                extracted_records=valid_cat_leaves,
                failure_record=fail_rec,
                repro_record=repro_by_cat.get(cat),
                pheno_record=pheno_by_cat.get(cat),
                repro_by_pid=repro_by_cat_pid.get(cat),
                pheno_by_pid=pheno_by_cat_pid.get(cat),
                contour_coords_by_leaf=coords_by_leaf,
                output_dir=t_dir.parent,
                workspace_root=ws,
            )
            migrated_leaves += len(valid_cat_leaves)

        # Ensure every voucher in `vouchers` has at least one row in `plant_individuals`
        # without injecting phantom `pid = 0` stubs onto vouchers that already hold `pid != 0`
        with conn:
            conn.execute(
                """
                INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id)
                SELECT v.catalogNumber, 0
                FROM vouchers v
                WHERE NOT EXISTS (
                    SELECT 1 FROM plant_individuals pi
                    WHERE pi.catalogNumber = v.catalogNumber
                );
                """
            )
        consolidate_false_multi_plant_splits_db(conn)

        counts["extracted_leaves"] = migrated_leaves
        cur_pi = conn.execute("SELECT COUNT(*) FROM plant_individuals;")
        counts["plant_individuals"] = int((cur_pi.fetchone() or [0])[0])

        # 4. Environmental cache (`data/environmental/fine_grain_edaphic_cache.csv`)
        if env_csv.is_file():
            env_rows = _read_csv_dicts(env_csv)
            counts["environmental_cache"] = upsert_environmental_cache_records(conn, env_rows, demo=False)
        else:
            counts["environmental_cache"] = 0

        # 5. Register active run & cohort (`run_id`)
        eff_run = str(run_id).strip() or "canonical"
        ensure_pipeline_run(conn, run_id=eff_run, run_dir=t_dir.parent)
        if vouchers_rows:
            counts["run_active_cohort"] = populate_run_active_cohort(
                conn, run_id=eff_run, voucher_rows=vouchers_rows, clear_existing=True
            )
        else:
            cur_all_v = conn.execute("SELECT catalogNumber, species_standardized, determiner_tier FROM vouchers;")
            v_stub = [dict(r) for r in cur_all_v.fetchall()]
            counts["run_active_cohort"] = populate_run_active_cohort(
                conn, run_id=eff_run, voucher_rows=v_stub, clear_existing=True
            )

        # 6. Optional downstream run-scoped tables if present in `t_dir`
        counts.update(
            ingest_existing_run_tables_into_db(
                conn=conn, tables_dir=t_dir, run_id=eff_run, demo=False
            )
        )

        logger.info("SQLite migration summary (%s): %s", eff_run, counts)
        return counts
    finally:
        conn.close()

__all__ = (
    "ingest_existing_csvs_to_db",
)
