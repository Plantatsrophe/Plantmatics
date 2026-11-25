#!/usr/bin/env python3
"""Warm-R tests for Step 04 intraspecific allometry, ecoregion CDA, and dynamic non-default taxa."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


class TestStep04AllometryAndEcoregionCDA(unittest.TestCase):
    """Warm-R tests for Step 04 intraspecific allometry, ecoregion CDA, and dynamic non-default taxa."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "04_gmm_morphotools.R")

    def test_dual_intraspecific_allometry_morphotype_and_supervised_epa_ecoregion_cda(self):
        """Verify single-taxon (--intraspecific-mode) allometry, partial AR_inv & Cauline_Only_Fallback retention, Track 1 Morphotype_K* CDA, and Track 2 EPA Level III Ecoregion CDA."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(202610)
        rng_ar = np.random.default_rng(777)
        eco_specs = [
            ("Piedmont", 35.9132, -79.0558, "North Carolina", "Orange", -3.2, -1.6),
            ("Middle_Atlantic_Coastal_Plain", 34.2257, -77.9447, "North Carolina", "New Hanover", 3.2, -1.4),
            ("South_Central_Plains", 31.6035, -94.6555, "Texas", "Nacogdoches", 0.2, 3.0),
        ]
        rows = []
        vouch_rows = []
        idx_counter = 0
        for eco_name, lat, lon, st_name, co_name, c1, c2 in eco_specs:
            for j in range(12):
                cat = f"DUAL_{idx_counter:03d}"
                idx_counter += 1
                tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                is_t2_fb = j == 11
                cs = float(24.0 + 1.1 * j)
                log_cs = float(np.log(cs))
                rec = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 0 if is_t2_fb else 3,
                    "morph_tier_source": "Cauline_Only_Fallback" if is_t2_fb else "Basal_Rosette",
                    "species_raw": "Packera anonyma",
                    "scientificName": "Packera anonyma",
                    "determiner_tier": tier,
                    "PC1": float(c1 + rng.normal(0.0, 0.22)),
                    "PC2": float(c2 + rng.normal(0.0, 0.22)),
                    "PC3": float(rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)),
                    "PC5": float(rng.normal(0.0, 0.12)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.9 + 0.25 * (c1 > 0) + rng.normal(0.0, 0.05)),
                    "solidity": float(0.88 + rng.normal(0.0, 0.02)),
                    "archetype_score": 0.94,
                    "AR_inv": np.nan if j % 5 == 0 else float(1.15 + 0.08 * c1 + rng_ar.normal(0.0, 0.015)),
                }
                for k in range(1, 13):
                    rec[f"A{k}"] = float((c1 * 0.06) / k + rng.normal(0.0, 0.003))
                    rec[f"B{k}"] = 0.0
                    rec[f"C{k}"] = 0.0
                    rec[f"D{k}"] = float((c2 * 0.06) / k + rng.normal(0.0, 0.003))
                rows.append(rec)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": "Packera anonyma",
                    "determiner_tier": tier,
                    "institutionCode": "NCU" if j % 2 == 0 else "TEX",
                    "stateProvince": st_name,
                    "county": co_name,
                    "decimalLatitude": lat if j < 10 else np.nan,
                    "decimalLongitude": lon if j < 10 else np.nan,
                    "latitude": lat if j < 10 else np.nan,
                    "longitude": lon if j < 10 else np.nan,
                    "regional_group": eco_name,
                    "doy": 118 + j,
                })

        for u_i in range(4):
            cat = f"DUAL_NA_{u_i:02d}"
            rec_na = dict(rows[u_i])
            rec_na.update({
                "catalogNumber": cat,
                "PC1": float(-2.9 + 0.2 * u_i),
                "PC2": float(-1.5 + 0.1 * u_i),
            })
            rows.append(rec_na)
            vouch_rows.append({
                "catalogNumber": cat,
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold",
                "institutionCode": "NCU",
                "stateProvince": "North Carolina",
                "county": "Wake",
                "decimalLatitude": np.nan,
                "decimalLongitude": np.nan,
                "latitude": np.nan,
                "longitude": np.nan,
                "regional_group": np.nan,
                "doy": 122,
            })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            fig_dir = tmp_p / "figures"
            rep_dir = tmp_p / "reports"
            fig_dir.mkdir(parents=True, exist_ok=True)
            rep_dir.mkdir(parents=True, exist_ok=True)

            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = rep_dir / "gmm_bayes_factors_summary.csv"
            allom_csv = rep_dir / "allometry.csv"
            allom_pdf = fig_dir / "fig04a_foliar_allometry_audit.png"
            eco_sum_csv = rep_dir / "ecoregion_supervised_cda_summary.csv"
            cda_pdf = fig_dir / "fig04b_gmm_cda_morphometrics.png"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(rep_dir / "batch_audit.csv"),
                    "--allometry-report", str(allom_csv),
                    "--allometry-pairwise", str(rep_dir / "pairwise.csv"),
                    "--allometry-plot", str(allom_pdf),
                    "--focal-taxon", "Packera anonyma",
                    "--intraspecific-mode",
                    "--permutations", "29",
                ],
                stage_key="step04",
                cwd=PROJECT_ROOT,
                skip_figures=False,
            )
            self.assertEqual(proc.returncode, 0, f"Dual/Intraspecific Step 04 failed:\n{proc.stderr}\n{proc.stdout}")

            # 0. Single-taxon allometry audit CSV & PDF + zero dropped Tier_1 / Cauline_Only_Fallback / partial-AR_inv vouchers
            self.assertTrue(allom_csv.exists(), "Allometry audit CSV must be exported in single-taxon mode")
            self.assertTrue(allom_pdf.exists() and allom_pdf.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), 40, "All 40 vouchers (including Cauline_Only_Fallback and partial AR_inv) must be retained")
            self.assertFalse(flags_df["gmm_cluster"].isna().any(), "No voucher should have NA gmm_cluster when AR_inv is partially missing")
            self.assertFalse(flags_df["cda_predicted_taxon"].isna().any(), "No voucher should have NA cda_predicted_taxon")

            # 1. Track 1: Unsupervised GMM BIC table evaluates G = 1:eff_k with Two_Delta_BIC vs G = 1
            self.assertTrue(bic_csv.exists())
            bic_df = pd.read_csv(bic_csv)
            self.assertEqual(int(bic_df["K"].iloc[0]), 1)
            self.assertAlmostEqual(float(bic_df["Two_Delta_BIC"].iloc[0]), 0.0, places=5)
            self.assertAlmostEqual(
                float(bic_df["Two_Delta_BIC"].iloc[2]),
                round(float(bic_df["BIC"].iloc[2] - bic_df["BIC"].iloc[0]), 2),
                places=2,
            )

            # 2. Track 1: cda_predicted_taxon is strictly Morphotype_K1..Kn (never regional_group)
            pred_morphs = set(flags_df["cda_predicted_taxon"].dropna().astype(str))
            self.assertGreaterEqual(len(pred_morphs), 2)
            for m_lbl in pred_morphs:
                self.assertRegex(m_lbl, r"^Morphotype_K[0-9]+$")

            # 3. Track 2: Supervised EPA Level III Ecoregion CDA columns in morphometric_flags.csv
            for req_col in ("ecoregion_cda_predicted", "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2"):
                self.assertIn(req_col, flags_df.columns)
                self.assertFalse(flags_df[req_col].isna().any(), f"{req_col} must be populated for active + NA-projected rows")
            pred_ecos = set(flags_df["ecoregion_cda_predicted"].dropna().astype(str))
            self.assertEqual(pred_ecos, {"Piedmont", "Middle_Atlantic_Coastal_Plain", "South_Central_Plains"})
            na_rows = flags_df[flags_df["catalogNumber"].str.startswith("DUAL_NA_")]
            self.assertEqual(len(na_rows), 4)
            self.assertTrue(na_rows["regional_group"].isna().all())
            self.assertTrue((na_rows["ecoregion_cda_predicted"] == "Piedmont").all())

            # 4. Track 2 summary report & fig04c_supervised_ecoregion_cda.{pdf,png}
            self.assertTrue(eco_sum_csv.exists(), "ecoregion_supervised_cda_summary.csv must be generated")
            sum_df = pd.read_csv(eco_sum_csv)
            self.assertEqual(set(sum_df["ecoregion"]), {"Piedmont", "Middle_Atlantic_Coastal_Plain", "South_Central_Plains"})
            for s_col in (
                "n_total_vouchers",
                "n_active_anchors",
                "loocv_accuracy",
                "permanova_r2",
                "permanova_p_value",
                "ari_morphotype_vs_ecoregion",
                "fisher_exact_p_value",
            ):
                self.assertIn(s_col, sum_df.columns)
            self.assertGreater(float(sum_df["loocv_accuracy"].iloc[0]), 0.85)
            self.assertGreater(float(sum_df["permanova_r2"].iloc[0]), 0.50)
            self.assertGreater(float(sum_df["ari_morphotype_vs_ecoregion"].iloc[0]), 0.50)
            self.assertLess(float(sum_df["fisher_exact_p_value"].iloc[0]), 0.01)

            for stem in (
                "fig04a_foliar_allometry_audit",
                "fig04b_gmm_cda_morphometrics",
                "fig04c_ecoregion_and_batch_audit",
            ):
                png_f = fig_dir / f"{stem}.png"
                self.assertTrue(png_f.exists() and png_f.stat().st_size > 10000, f"Missing/empty {png_f}")
            self.assertEqual(list(fig_dir.glob("*.pdf")), [])

    def test_dynamic_non_default_taxa_cda_second_taxon_and_multi_plant_petioles(self):
        """Verify dynamic non-default Packera species discovery, singleton passive projection, cda_second_taxon export, and multi-plant (catalogNumber, plant_individual_id) petiole aggregation."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(2026)
        non_default_specs = [
            ("Packera schweinitziana", "Senecio schweinitzianus", (-4.5, -2.5, 0.5)),
            ("Packera serpenticola", "Packera eburnea", (4.5, -2.0, -0.5)),
            ("Packera millefolium", "Senecio millefolium", (-1.5, 4.0, -2.0)),
            ("Packera aurea", "Senecio aureus", (2.5, 3.5, 2.5)),
        ]
        rows = []
        vouch_rows = []
        for idx, (canon_sp, raw_syn, (c1, c2, c3)) in enumerate(non_default_specs):
            for j in range(10):
                is_multi = idx in (0, 1) and j == 0
                cat = "MULTI_SHEET_01" if is_multi else f"ENA_{idx}_{j:02d}"
                pid = (idx + 1) if is_multi else 0
                tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 8 else "Tier_3_Bronze")
                cs = float(28.0 + 1.2 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": pid,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": raw_syn if j % 2 == 0 else canon_sp,
                    "scientificName": canon_sp,
                    "determiner_tier": tier,
                    "PC1": float(c1 + rng.normal(0.0, 0.15)),
                    "PC2": float(c2 + rng.normal(0.0, 0.15)),
                    "PC3": float(c3 + rng.normal(0.0, 0.12)),
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.8 + 0.25 * idx + rng.normal(0.0, 0.04)),
                    "solidity": float(0.42 if canon_sp == "Packera millefolium" else 0.89),
                    "archetype_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((c1 * 0.06) / k + rng.normal(0.0, 0.003))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((c2 * 0.06) / k + rng.normal(0.0, 0.003))
                rows.append(row)
                if not any(v["catalogNumber"] == cat for v in vouch_rows):
                    vouch_rows.append({
                        "catalogNumber": cat,
                        "species_raw": row["species_raw"],
                        "determiner_tier": tier,
                        "institutionCode": "NCU" if len(vouch_rows) % 2 == 0 else "NY",
                        "stateProvince": "North Carolina",
                        "latitude": 35.6,
                        "longitude": -82.5,
                        "doy": 132,
                    })

        # Add a singleton rare outgroup row (n = 1) to verify passive projection without singular W
        singleton_row = dict(rows[2])
        singleton_row.update({
            "catalogNumber": "ENA_SINGLETON_01",
            "plant_individual_id": 0,
            "species_raw": "Senecio glabellus",
            "scientificName": "Packera glabella",
            "determiner_tier": "Tier_1_Gold",
            "PC1": 0.0,
            "PC2": 0.0,
            "PC3": 0.0,
        })
        rows.append(singleton_row)
        vouch_rows.append({
            "catalogNumber": "ENA_SINGLETON_01",
            "species_raw": "Senecio glabellus",
            "determiner_tier": "Tier_1_Gold",
            "institutionCode": "NCU",
            "stateProvince": "North Carolina",
            "latitude": 35.6,
            "longitude": -82.5,
            "doy": 132,
        })

        man_rows = [
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 1, "leaf_id": 1, "petiole_length_mm": 32.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": round(32.0 / 28.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 1, "leaf_id": 2, "petiole_length_mm": 32.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": round(32.0 / 28.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 2, "leaf_id": 1, "petiole_length_mm": 12.0, "blade_length_mm": 30.0, "petiole_to_blade_ratio": round(12.0 / 30.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 2, "leaf_id": 2, "petiole_length_mm": 12.0, "blade_length_mm": 30.0, "petiole_to_blade_ratio": round(12.0 / 30.0, 4)},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            man_csv = tmp_p / "extracted_leaf_manifest.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            allom_csv = tmp_p / "allometry_significance_audit.csv"
            allom_pdf = tmp_p / "fig04a_foliar_allometry_audit.png"
            cda_pdf = tmp_p / "fig04b_gmm_cda_morphometrics.png"

            df_in = pd.DataFrame(rows)
            df_in.to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)
            pd.DataFrame(man_rows).to_csv(man_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(allom_csv),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
                    "--allometry-plot", str(allom_pdf),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 failed on non-default Packera taxa:\n{proc.stderr}\n{proc.stdout}")

            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), len(df_in), "All vouchers including multi-plant individuals and singleton outgroup must be retained.")

            expected_four = {
                "Packera schweinitziana",
                "Packera serpenticola",
                "Packera millefolium",
                "Packera aurea",
            }
            pred_taxa = set(flags_df["cda_predicted_taxon"].dropna().unique())
            self.assertEqual(pred_taxa, expected_four)
            self.assertTrue(
                pred_taxa.isdisjoint({"Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"}),
                f"Default 4-species names must not appear when absent from dataset: {pred_taxa}",
            )

            # Verify cda_second_taxon export
            self.assertIn("cda_second_taxon", flags_df.columns)
            self.assertTrue(flags_df["cda_second_taxon"].notna().all())
            self.assertTrue((flags_df["cda_second_taxon"] != flags_df["cda_predicted_taxon"]).all())
            self.assertTrue(set(flags_df["cda_second_taxon"].unique()).issubset(expected_four))

            # Verify high classification accuracy on the 4 well-separated species
            core_rows = flags_df[flags_df["species_standardized"].isin(expected_four)]
            acc = float((core_rows["species_standardized"] == core_rows["cda_predicted_taxon"]).mean())
            self.assertGreaterEqual(acc, 0.95, f"Expected >=95% CDA classification accuracy on 4 distinct taxa, got {acc:.3f}")

            # Verify singleton Packera glabella row was projected passively with finite 3D canonical scores
            sing_out = flags_df[flags_df["catalogNumber"] == "ENA_SINGLETON_01"]
            self.assertEqual(len(sing_out), 1)
            self.assertEqual(sing_out["species_standardized"].iloc[0], "Packera glabella")
            self.assertTrue(bool(sing_out["is_passive"].iloc[0]))
            for col in ("can1", "can2", "can3", "cda_posterior_prob"):
                self.assertTrue(np.isfinite(float(sing_out[col].iloc[0])), f"Singleton {col} must be finite.")

            # Verify multi-plant sheet individuals (pid=1, pid=2) are never collapsed and retain distinct petioles
            multi_rows = flags_df[flags_df["catalogNumber"] == "MULTI_SHEET_01"].sort_values("plant_individual_id")
            self.assertEqual(len(multi_rows), 2, "Multi-plant sheet individuals (pid=1, pid=2) must not be collapsed.")
            self.assertAlmostEqual(float(multi_rows.iloc[0]["petiole_length_mm"]), 32.0, places=2)
            self.assertAlmostEqual(float(multi_rows.iloc[1]["petiole_length_mm"]), 12.0, places=2)
            self.assertAlmostEqual(float(multi_rows.iloc[0]["petiole_to_blade_ratio"]), round(32.0 / 28.0, 4), places=3)
            self.assertAlmostEqual(float(multi_rows.iloc[1]["petiole_to_blade_ratio"]), round(12.0 / 30.0, 4), places=3)


if __name__ == "__main__":
    unittest.main()
