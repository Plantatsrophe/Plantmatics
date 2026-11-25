#!/usr/bin/env python3
"""
Relational contract tests for multi-plant (catalogNumber, plant_individual_id)
retention, non-georeferenced historic voucher preservation across Step 05/07,
and v_integrative_triage_master Step 04 parity & Hopkins precedence.
"""

from __future__ import annotations

import csv
import os
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from scripts.tests._r_session import rscript_available, run_warm_r_cmd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    ensure_pipeline_run,
    get_db_connection,
    populate_run_active_cohort,
    upsert_environmental_cache_records,
    upsert_rows_generic,
    upsert_vouchers_records,
)


def _mk_voucher(
    cat: str,
    inst: str,
    sp: str,
    tier: str,
    lat: Optional[float],
    lon: Optional[float],
    unc: Optional[float],
    yr: int,
    mo: int,
    dy: int,
    doy: int,
    st: str,
    reg: str,
) -> Dict[str, Any]:
    return {
        "catalogNumber": cat, "institutionCode": inst, "scientificName": sp,
        "species_raw": sp, "species_standardized": sp, "determiner_tier": tier,
        "decimalLatitude": lat, "latitude": lat, "decimalLongitude": lon, "longitude": lon,
        "coordinateUncertainty": unc, "year": yr, "month": mo, "day": dy,
        "eventDate": f"{yr:04d}-{mo:02d}-{dy:02d}", "doy": doy,
        "stateProvince": st, "regional_group": reg,
    }


class TestSQLMultiplantAndMasterView(unittest.TestCase):
    """Integration tests for multi-plant/non-georeferenced retention and master view contracts."""

    def test_multi_plant_and_non_georeferenced_retention_across_step06_and_step07(self) -> None:
        """Verify multi-plant (catalogNumber, plant_individual_id) & non-georeferenced retention across Step 05, Step 07, and v_integrative_triage_master."""
        if not rscript_available():
            self.skipTest("Rscript not installed or fast mode enabled.")
        rscript_bin = "Rscript"

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d) / "ws_master"
            tables_dir = ws / "data" / "tables"
            env_dir = ws / "data" / "environmental"
            reports_dir = ws / "outputs" / "reports"
            figures_dir = ws / "outputs" / "figures"
            for d in (tables_dir, env_dir, reports_dir, figures_dir):
                d.mkdir(parents=True, exist_ok=True)

            db_path = ws / "data" / "plantmatics.sqlite"
            run_id = "test_master_multi_nogeo"

            voucher_specs: List[Dict[str, Any]] = [
                _mk_voucher("MULTI-001", "NCU", "Packera aurea", "Tier_2_Silver", 35.91000, -79.05000, 200.0, 2020, 4, 20, 111, "North Carolina", "Piedmont"),
                _mk_voucher("HIST-NOGEO-99", "NY", "Packera obovata", "Tier_3_Bronze", None, None, None, 1898, 5, 4, 124, "Virginia", ""),
            ]
            for i in range(1, 6):
                voucher_specs.append(
                    _mk_voucher(
                        f"ANC-AUR-{i}", "NCU" if i % 2 == 0 else "DUKE", "Packera aurea", "Tier_1_Gold",
                        round(35.60 + 0.08 * i, 5), round(-79.40 + 0.08 * i, 5), 150.0, 2018, 4, 15 + i, 105 + i, "North Carolina", "Piedmont",
                    )
                )
                voucher_specs.append(
                    _mk_voucher(
                        f"ANC-OBO-{i}", "US" if i % 2 == 0 else "MO", "Packera obovata", "Tier_1_Gold",
                        round(37.20 + 0.08 * i, 5), round(-80.40 + 0.08 * i, 5), 150.0, 2018, 4, 22 + i, 112 + i, "Virginia", "Ridge and Valley",
                    )
                )

            vouchers_csv = tables_dir / "curated_vouchers.csv"
            with open(vouchers_csv, "w", encoding="utf-8", newline="") as vf:
                vw = csv.DictWriter(vf, fieldnames=list(voucher_specs[0].keys()))
                vw.writeheader()
                for row in voucher_specs:
                    vw.writerow({k: ("" if v is None else v) for k, v in row.items()})

            otu_specs: List[Tuple[str, int, str, str]] = [
                ("MULTI-001", 0, "Packera aurea", "Packera aurea"),
                ("MULTI-001", 1, "Packera aurea", "Packera obovata"),
                ("HIST-NOGEO-99", 0, "Packera obovata", "Packera obovata"),
            ]
            for i in range(1, 6):
                otu_specs.append((f"ANC-AUR-{i}", 0, "Packera aurea", "Packera aurea"))
                otu_specs.append((f"ANC-OBO-{i}", 0, "Packera obovata", "Packera obovata"))

            conn = get_db_connection(db_path=db_path, init_schema=True)
            try:
                upsert_vouchers_records(conn, voucher_specs)
                with conn:
                    for cat, pid, _, _ in otu_specs:
                        conn.execute(
                            "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, petiole_length_mm, blade_length_mm, petiole_to_blade_ratio, plant_height_mm) "
                            "VALUES (?, ?, ?, ?, ?, ?);",
                            (cat, pid, 30.0 if pid == 0 else 11.0, 34.0 if pid == 0 else 42.0, round((30.0 / 34.0) if pid == 0 else (11.0 / 42.0), 4), 360.0 + 25.0 * pid),
                        )
                        conn.execute(
                            "INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count, n_anthesis, n_bud, n_fruiting, voucher_phenological_state) "
                            "VALUES (?, ?, 12, 8, 2, 2, 'anthesis');",
                            (cat, pid),
                        )

                ensure_pipeline_run(conn, run_id=run_id, run_dir=ws)
                n_coh = populate_run_active_cohort(conn, run_id=run_id, voucher_rows=voucher_specs, clear_existing=True)
                self.assertEqual(n_coh, 13)

                otu_rows: List[Dict[str, Any]] = []
                env_rows: List[Dict[str, Any]] = []
                v_by_cat = {r["catalogNumber"]: r for r in voucher_specs}

                for idx_o, (cat, pid, given_sp, true_sp) in enumerate(otu_specs):
                    is_aur = true_sp == "Packera aurea"
                    second_sp = "Packera obovata" if is_aur else "Packera aurea"
                    otu_rows.append({
                        "run_id": run_id, "catalogNumber": cat, "plant_individual_id": pid,
                        "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3,
                        "PC1": (-1.6 + 0.05 * idx_o) if is_aur else (1.7 + 0.05 * idx_o),
                        "PC2": 0.12 * ((idx_o % 3) - 1), "PC3": 0.05 * ((idx_o % 2) - 0.5),
                        "can1": (-2.4 + 0.06 * idx_o) if is_aur else (2.5 + 0.06 * idx_o),
                        "can2": 0.15 * ((idx_o % 3) - 1), "can3": 0.05,
                        "gmm_cluster": 1 if is_aur else 2, "gmm_uncertainty": 0.04,
                        "cda_predicted_taxon": true_sp, "cda_posterior_prob": 0.94,
                        "cda_second_taxon": second_sp, "cda_second_posterior_prob": 0.06,
                        "cda_posterior_margin": 0.88, "loocv_accuracy": 0.95,
                        "mahalanobis_d2": 1.4 if given_sp == true_sp else 14.2,
                        "morphometric_outlier_flag": 0 if given_sp == true_sp else 1,
                        "misidentification_flag": 0 if given_sp == true_sp else 1,
                    })
                    v_rec = v_by_cat[cat]
                    if v_rec["latitude"] is not None and pid == 0:
                        env_rows.append({
                            "catalogNumber": cat, "latitude": v_rec["latitude"], "longitude": v_rec["longitude"],
                            "fine_soil_ph": 5.4 + (0.05 * idx_o if is_aur else 1.4 + 0.04 * idx_o),
                            "fine_soil_sand": 42.0 + (1.2 * idx_o if is_aur else 18.0 + 0.9 * idx_o),
                            "fine_soil_clay": 24.0 - 0.5 * idx_o, "fine_bedrock_depth_cm": 135.0 if is_aur else 48.0,
                            "ssurgo_rock_outcrop_pct": 2.0 if is_aur else 22.0,
                            "elevation_m": 180.0 + 5.0 * idx_o if is_aur else 460.0 + 8.0 * idx_o,
                            "topo_slope_deg": 4.5 if is_aur else 16.0, "topo_tpi_30m": -0.8 if is_aur else 2.4,
                            "topo_hli_30m": 0.72 if is_aur else 0.94, "topo_twi_30m": 9.8 if is_aur else 5.2,
                            "sgmc_lithology_class": "Felsic_Crystalline" if is_aur else "Carbonate",
                            "musym": "WkB" if is_aur else "GlC",
                            "muname": "Wake-Wateree complex" if is_aur else "Gladeville-Rock outcrop",
                            "outcrop_class": "Granitic flatrock" if is_aur else "Limestone glade",
                            "is_rock_outcrop": False if is_aur else True,
                            "soil_ph": 5.4 + (0.05 * idx_o if is_aur else 1.4 + 0.04 * idx_o),
                            "soil_cec": 11.5 if is_aur else 22.0,
                            "soil_sand": 42.0 + (1.2 * idx_o if is_aur else 18.0 + 0.9 * idx_o),
                            "soil_bulk_density": 1.35 if is_aur else 1.28,
                            "bio1_temp_mean": 14.8 if is_aur else 11.6,
                            "bio4_temp_seasonality": 740.0 + 5.0 * idx_o,
                            "bio12_precip_annual": 1180.0 if is_aur else 980.0,
                            "bio15_precip_seasonality": 14.0 + 0.3 * idx_o,
                            "bio18_precip_warmest_qtr": 345.0 if is_aur else 290.0,
                            "bio19_precip_coldest_qtr": 265.0 if is_aur else 215.0,
                            "is_training_eligible_coord": True,
                        })

                upsert_rows_generic(conn, "otu_morphometrics", otu_rows, ("run_id", "catalogNumber", "plant_individual_id"))
                upsert_environmental_cache_records(conn, env_rows, demo=False)

                morph_csv = tables_dir / "morphometric_flags.csv"
                with open(morph_csv, "w", encoding="utf-8", newline="") as mf:
                    mw = csv.DictWriter(mf, fieldnames=list(otu_rows[0].keys()))
                    mw.writeheader()
                    mw.writerows(otu_rows)

                view_rows_pre = conn.execute(
                    "SELECT catalogNumber, plant_individual_id, cda_predicted_taxon, fine_soil_ph "
                    "FROM v_integrative_triage_master WHERE run_id = ? ORDER BY catalogNumber, plant_individual_id;",
                    (run_id,),
                ).fetchall()
                self.assertEqual(len(view_rows_pre), 13)
                pre_by_key = {(str(r[0]), int(r[1])): dict(r) for r in view_rows_pre}
                self.assertEqual(pre_by_key[("MULTI-001", 0)]["cda_predicted_taxon"], "Packera aurea")
                self.assertEqual(pre_by_key[("MULTI-001", 1)]["cda_predicted_taxon"], "Packera obovata")
                self.assertIn(("HIST-NOGEO-99", 0), pre_by_key)
                self.assertEqual(pre_by_key[("HIST-NOGEO-99", 0)]["cda_predicted_taxon"], "Packera obovata")
                self.assertIsNone(pre_by_key[("HIST-NOGEO-99", 0)]["fine_soil_ph"])
            finally:
                conn.close()

            env = os.environ.copy()
            env["PLANTMATICS_DB_PATH"] = str(db_path)
            env["RUN_ID"] = run_id
            env["PLANTMATICS_RUN_DIR"] = str(ws)

            conflict_csv = tables_dir / "multimodal_conflict_flags.csv"
            srf_summary_csv = reports_dir / "multimodal_spatial_rf_summary.csv"
            cmd_step05 = [
                rscript_bin, str(PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"),
                "--vouchers", str(vouchers_csv), "--morphometrics", str(tables_dir / "morphometric_flags.csv"),
                "--env-dir", str(env_dir), "--ssurgo-table", str(tables_dir / "ssurgo_edaphic_validation.csv"),
                "--contingency-table", str(reports_dir / "micro_edaphic_outcrop_contingency.csv"),
                "--pheno-anomalies", str(tables_dir / "phenological_anomalies.csv"),
                "--pheno-summary", str(reports_dir / "phenological_anomaly_summary.csv"),
                "--pheno-plot", str(figures_dir / "fig05b_hopkins_circular_phenology.png"),
                "--output-flags", str(conflict_csv), "--output-summary", str(srf_summary_csv),
                "--output-plot", str(figures_dir / "fig05a_ecological_niche_and_rf_importance.png"),
                "--permutations", "9", "--n-trees", "50",
            ]
            proc05 = run_warm_r_cmd(cmd_step05, stage_key="step05", cwd=PROJECT_ROOT, env=env, skip_figures=True)
            self.assertEqual(proc05.returncode, 0, f"Step 05 failed:\n{proc05.stdout}\n{proc05.stderr}")
            self.assertTrue(conflict_csv.is_file())

            (tables_dir / "morphometric_flags.csv").unlink()
            conflict_csv.unlink()

            triage_queue_csv = tables_dir / "triage_queue.csv"
            rev_report_md = reports_dir / "Packera_Taxonomic_Revision_Summary.md"
            rev_html = reports_dir / "Triage_Dashboard_Interactive.html"
            rev_fig_png = figures_dir / "fig07a_integrative_taxonomic_synthesis.png"
            cmd_step07 = [
                rscript_bin, str(PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"),
                "--vouchers", str(vouchers_csv), "--morphometrics", str(tables_dir / "morphometric_flags.csv"),
                "--multimodal-flags", str(conflict_csv), "--niche-summary", str(srf_summary_csv),
                "--output-queue", str(triage_queue_csv), "--output-plot", str(rev_fig_png),
                "--output-report", str(rev_report_md), "--output-html", str(rev_html),
            ]
            proc07 = run_warm_r_cmd(cmd_step07, stage_key="step07", cwd=PROJECT_ROOT, env=env, skip_figures=True)
            combined07 = proc07.stdout + "\n" + proc07.stderr
            self.assertEqual(proc07.returncode, 0, f"Step 07 failed:\n{combined07}")
            self.assertIn("v_integrative_triage_master", combined07)
            self.assertTrue(triage_queue_csv.is_file())
            self.assertTrue(rev_html.is_file())

            with open(triage_queue_csv, "r", encoding="utf-8", newline="") as tf:
                triage_rows = list(csv.DictReader(tf))
            self.assertEqual(len(triage_rows), 13)
            by_key = {(r["catalogNumber"], int(float(r["plant_individual_id"]))): r for r in triage_rows}

            self.assertIn(("MULTI-001", 0), by_key)
            self.assertIn(("MULTI-001", 1), by_key)
            self.assertEqual(by_key[("MULTI-001", 0)]["recommended_determination"], "Packera aurea")
            self.assertEqual(by_key[("MULTI-001", 1)]["recommended_determination"], "Packera obovata")
            self.assertEqual(by_key[("MULTI-001", 1)]["taxonomic_treatment_subtype"], "Misidentification_Reassign")

            self.assertIn(("HIST-NOGEO-99", 0), by_key)
            nogeo = by_key[("HIST-NOGEO-99", 0)]
            self.assertIn(str(nogeo["is_completed_voucher"]).strip().upper(), {"TRUE", "1"})
            self.assertTrue(str(nogeo["recommended_determination"]).startswith("Packera obovata"))

            conn_post = get_db_connection(db_path=db_path, read_only=True)
            try:
                n_ct = conn_post.execute("SELECT COUNT(*) FROM curation_triage WHERE run_id = ?;", (run_id,)).fetchone()[0]
                self.assertEqual(int(n_ct), 13)
                nogeo_db = conn_post.execute(
                    "SELECT recommended_determination, is_completed_voucher FROM v_integrative_triage_master "
                    "WHERE run_id = ? AND catalogNumber = 'HIST-NOGEO-99' AND plant_individual_id = 0;",
                    (run_id,),
                ).fetchone()
                self.assertIsNotNone(nogeo_db)
                self.assertTrue(str(nogeo_db[0]).startswith("Packera obovata"))
                self.assertEqual(int(nogeo_db[1]), 1)
            finally:
                conn_post.close()

    def test_master_view_step04_parity_and_upstream_hopkins_precedence(self) -> None:
        """Verify v_integrative_triage_master exposes Step 04 parity columns and fresh Step 05 hopkins_adjusted_doy overrides stale curation_triage."""
        conn = get_db_connection(":memory:", init_schema=True)
        try:
            mv_cols = {r[1] for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
            om_cols = {r[1] for r in conn.execute("PRAGMA table_info(otu_morphometrics);").fetchall()}
            self.assertIn("scale_source", om_cols)
            for req in ("is_passive", "scale_source", "AR_inv_std", "reproductive_trait_status", "ecoregion_can1", "ecoregion_can2"):
                self.assertIn(req, mv_cols)

            conn.execute("INSERT INTO pipeline_runs (run_id) VALUES ('run_hopkins_test');")
            conn.execute(
                "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier) "
                "VALUES ('NCU-HOP-1', 'Packera aurea', 'Packera aurea', 'Tier_3_Bronze');"
            )
            conn.execute(
                "INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier, is_promoted_bronze, authority_weight) "
                "VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 'Packera aurea', 'Tier_3_Bronze', 0, 0.0);"
            )
            conn.execute(
                "INSERT INTO otu_morphometrics (run_id, catalogNumber, plant_individual_id, is_passive, scale_source, AR_inv_std, reproductive_trait_status, ecoregion_can1, ecoregion_can2) "
                "VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 1, 'ruler_calibrated', 1.25, 'Measured', -0.45, 0.85);"
            )
            conn.execute(
                "INSERT INTO spatial_rf_predictions (run_id, catalogNumber, plant_individual_id, hopkins_adjusted_doy) "
                "VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 134.5);"
            )
            conn.execute(
                "INSERT INTO curation_triage (run_id, catalogNumber, plant_individual_id, hopkins_adjusted_doy) "
                "VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 99.0);"
            )
            row = conn.execute(
                "SELECT is_passive, scale_source, AR_inv_std, reproductive_trait_status, ecoregion_can1, ecoregion_can2, hopkins_adjusted_doy "
                "FROM v_integrative_triage_master WHERE run_id = 'run_hopkins_test' AND catalogNumber = 'NCU-HOP-1';"
            ).fetchone()
            self.assertIsNotNone(row)
            self.assertEqual(int(row[0]), 1)
            self.assertEqual(str(row[1]), "ruler_calibrated")
            self.assertAlmostEqual(float(row[2]), 1.25)
            self.assertEqual(str(row[3]), "Measured")
            self.assertAlmostEqual(float(row[4]), -0.45)
            self.assertAlmostEqual(float(row[5]), 0.85)
            self.assertAlmostEqual(float(row[6]), 134.5)
        finally:
            conn.close()


if __name__ == "__main__":
    unittest.main()
