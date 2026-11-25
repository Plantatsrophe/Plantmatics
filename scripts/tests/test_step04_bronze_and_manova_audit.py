#!/usr/bin/env python3
"""Warm-R tests for Step 04 Bronze invariance, multivariate MANOVA, >5 PC subspace, and BH batch audit."""

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


class TestStep04BronzeAndManovaAudit(unittest.TestCase):
    """Warm-R tests for Step 04 Bronze invariance, multivariate MANOVA, >5 PC subspace, and BH batch audit."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "04_gmm_morphotools.R")

    def test_tier3_bronze_invariance_and_multivariate_manova_audit(self):
        """Verify Tier_3_Bronze zero-influence invariance across allometry, GMM, CDA, and MANOVA, plus batch audit R^2 < 0.05 and CDA subtitles."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(314)
        specs = [
            ("Packera anonyma", "Piedmont", 35.9132, -79.0558, "NC", -5.5, -1.5, 11),
            ("Packera dubia", "Middle_Atlantic_Coastal_Plain", 34.2257, -77.9447, "NC", 5.5, -1.5, 11),
            ("Packera paupercula", "South_Central_Plains", 31.6035, -94.6555, "TX", 0.0, 5.0, 11),
        ]
        rows = []
        vouch_rows = []
        for s_idx, (sp, eco_name, lat0, lon0, st_code, c1, c2, n_v) in enumerate(specs):
            for j in range(n_v):
                cat = f"MANOVA_{s_idx}_{j:02d}"
                det_tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                cs = float(28.0 + 0.8 * j)
                log_cs = float(np.log(cs))
                pc1 = float(c1 + 0.25 * (log_cs - 3.4) + rng.normal(0.0, 0.25))
                pc2 = float(c2 + rng.normal(0.0, 0.25))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": det_tier,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": float(0.4 * s_idx + rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(2.0 + 0.3 * s_idx + rng.normal(0.0, 0.05)),
                    "solidity": float(0.89 + 0.02 * s_idx),
                    "archetype_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((-1.25 - 0.06 * pc1 + 0.03 * log_cs) / k + rng.normal(0.0, 0.003))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((0.46 + 0.05 * pc2 + 0.02 * log_cs) / k + rng.normal(0.0, 0.003))
                rows.append(row)
                vouch_rows.append(
                    {
                        "catalogNumber": cat,
                        "species_raw": sp,
                        "determiner_tier": det_tier,
                        "institutionCode": "NCU" if j % 2 == 0 else "BRIT",
                        "stateProvince": st_code,
                        "county": "",
                        "regional_group": eco_name,
                        "decimalLatitude": float(lat0 + rng.normal(0.0, 0.03)),
                        "decimalLongitude": float(lon0 + rng.normal(0.0, 0.03)),
                        "doy": 115 + j,
                    }
                )

        # Add 2 Blue_Ridge Gold vouchers (< 3 threshold -> excluded from regional_group MANOVA) and 2 NA-ecoregion vouchers
        for extra_idx, (rg_val, lat_v, lon_v) in enumerate(
            [
                ("Blue_Ridge", 35.5951, -82.5515),
                ("Blue_Ridge", 35.6100, -82.5400),
                ("", np.nan, np.nan),
                ("", np.nan, np.nan),
            ]
        ):
            cat = f"MANOVA_EXTRA_{extra_idx:02d}"
            row = dict(rows[0])
            row["catalogNumber"] = cat
            row["determiner_tier"] = "Tier_1_Gold"
            row["PC1"] = float(-5.4 + 0.1 * extra_idx)
            row["PC2"] = float(-1.4 - 0.1 * extra_idx)
            rows.append(row)
            vouch_rows.append(
                {
                    "catalogNumber": cat,
                    "species_raw": "Packera anonyma",
                    "determiner_tier": "Tier_1_Gold",
                    "institutionCode": "NCU",
                    "stateProvince": "NC",
                    "county": "",
                    "regional_group": rg_val,
                    "decimalLatitude": lat_v,
                    "decimalLongitude": lon_v,
                    "doy": 120,
                }
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)

            def _run_case(sub_name: str, perturb_bronze: bool = False):
                case_dir = tmp_p / sub_name
                fig_dir = case_dir / "outputs" / "figures"
                rep_dir = case_dir / "outputs" / "reports"
                tbl_dir = case_dir / "data" / "tables"
                for d in (fig_dir, rep_dir, tbl_dir):
                    d.mkdir(parents=True, exist_ok=True)
                df_h = pd.DataFrame(rows).copy()
                df_v = pd.DataFrame(vouch_rows).copy()
                if perturb_bronze:
                    b_mask = df_h["determiner_tier"] == "Tier_3_Bronze"
                    df_h.loc[b_mask, ["PC1", "PC2", "PC3"]] += 45.0
                    df_h.loc[b_mask, "log_centsize"] += 2.5
                    for k in range(1, 13):
                        df_h.loc[b_mask, f"A{k}"] += 0.5
                        df_h.loc[b_mask, f"D{k}"] -= 0.5
                    df_v.loc[df_v["determiner_tier"] == "Tier_3_Bronze", "regional_group"] = "Blue_Ridge"
                harm_csv = tbl_dir / "leaf_efa_harmonics.csv"
                vouch_csv = tbl_dir / "curated_vouchers.csv"
                flags_csv = tbl_dir / "morphometric_flags.csv"
                bic_csv = rep_dir / "gmm_bayes_factors_summary.csv"
                allom_csv = rep_dir / "allometry.csv"
                batch_csv = rep_dir / "batch_audit.csv"
                cda_pdf = fig_dir / "fig04b_gmm_cda_morphometrics.png"
                df_h.to_csv(harm_csv, index=False)
                df_v.to_csv(vouch_csv, index=False)
                proc = run_warm_r_cmd(
                    [
                        "Rscript",
                        str(self.r_script_path),
                        "--vouchers", str(vouch_csv),
                        "--harmonics", str(harm_csv),
                        "--output-flags", str(flags_csv),
                        "--output-report", str(bic_csv),
                        "--output-plot", str(cda_pdf),
                        "--output-batch-audit", str(batch_csv),
                        "--allometry-report", str(allom_csv),
                        "--allometry-pairwise", str(rep_dir / "pairwise.csv"),
                        "--allometry-plot", str(fig_dir / "fig04a_foliar_allometry_audit.png"),
                        "--correct-allometry", "true",
                        "--permutations", "19",
                    ],
                    stage_key="step04",
                    cwd=PROJECT_ROOT,
                    skip_figures=perturb_bronze,
                )
                self.assertEqual(proc.returncode, 0, f"Step 04 run failed ({sub_name}):\n{proc.stderr}\n{proc.stdout}")
                return tbl_dir, rep_dir, fig_dir

            tbl_base, rep_base, fig_base = _run_case("base", perturb_bronze=False)
            flags1 = pd.read_csv(tbl_base / "morphometric_flags.csv")
            bic1 = pd.read_csv(rep_base / "gmm_bayes_factors_summary.csv")
            allom1 = pd.read_csv(rep_base / "allometry.csv")
            batch1 = pd.read_csv(rep_base / "batch_audit.csv")

            # 1. Verify Empirical Bayes active/passive authority weights
            bronze_rows = flags1[flags1["determiner_tier"] == "Tier_3_Bronze"]
            active_rows = flags1[flags1["authority_weight"] > 0.0]
            gold_rows = active_rows[active_rows["determiner_tier"] == "Tier_1_Gold"]
            silver_rows = active_rows[active_rows["determiner_tier"] == "Tier_2_Silver"]
            self.assertTrue((bronze_rows["authority_weight"] == 0.0).all())
            self.assertTrue((bronze_rows["is_passive"] == True).all())
            self.assertAlmostEqual(float(active_rows["authority_weight"].mean()), 1.0, places=4)
            self.assertGreater(float(gold_rows["authority_weight"].mean()), float(silver_rows["authority_weight"].mean()))

            # 2. Verify institutional batch effect report R^2 < 0.05 and allometry ANOVA schema
            for resp_pc in ("PC1", "PC2"):
                b_row = batch1[(batch1["Response"] == resp_pc) & (batch1["Term"] == "institutionCode")]
                self.assertEqual(len(b_row), 1)
                self.assertLess(float(b_row["R_squared"].iloc[0]), 0.05)
                self.assertIn("Pass", str(b_row["Batch_Effect_Status"].iloc[0]))
            for req_term in ("log_centsize", "scientificName", "Residual"):
                self.assertIn(req_term, allom1["Term"].tolist())

            # 3. Verify omnibus & pairwise MANOVA tables and CDA PDF subtitles
            grp_df = pd.read_csv(rep_base / "manova_group_significance.csv")
            pw_df = pd.read_csv(rep_base / "manova_pairwise_posthoc.csv")
            req_grp_cols = {
                "grouping_factor", "pillai_trace", "wilks_lambda", "approx_F",
                "num_df", "den_df", "partial_eta2", "p_value", "q_value",
            }
            self.assertTrue(req_grp_cols.issubset(set(grp_df.columns)))
            self.assertEqual(set(grp_df["grouping_factor"]), {"species_standardized", "regional_group"})
            for _, r in grp_df.iterrows():
                self.assertGreater(float(r["pillai_trace"]), 0.20)
                self.assertGreater(float(r["wilks_lambda"]), 0.0)
                self.assertLess(float(r["wilks_lambda"]), 0.90)
                self.assertGreater(float(r["approx_F"]), 3.0)
                self.assertGreater(float(r["partial_eta2"]), 0.10)
                self.assertLess(float(r["p_value"]), 0.01)
                self.assertGreaterEqual(float(r["q_value"]) + 1e-12, float(r["p_value"]))

            rg_row = grp_df[grp_df["grouping_factor"] == "regional_group"].iloc[0]
            self.assertEqual(int(rg_row["n_groups"]), 3)

            sp_pw = pw_df[pw_df["grouping_factor"] == "species_standardized"]
            rg_pw = pw_df[pw_df["grouping_factor"] == "regional_group"]
            cl_pw = pw_df[pw_df["grouping_factor"] == "gmm_cluster"]
            self.assertEqual(len(sp_pw), 3)
            self.assertEqual(len(rg_pw), 3)
            self.assertEqual(len(cl_pw), 0, "Post-selection gmm_cluster p-values must be excluded from pipeline MANOVA")
            self.assertNotIn("Blue_Ridge", set(rg_pw["group_1"]).union(set(rg_pw["group_2"])))

            for _, pr in pw_df.iterrows():
                self.assertGreater(float(pr["hotelling_t2"]), 0.0)
                self.assertAlmostEqual(float(pr["pillai_trace"] + pr["wilks_lambda"]), 1.0, places=4)
                self.assertAlmostEqual(float(pr["partial_eta2"]), float(pr["pillai_trace"]), places=4)
                self.assertGreaterEqual(float(pr["q_value"]) + 1e-12, float(pr["p_value"]))

            self.assertTrue((fig_base / "fig04b_gmm_cda_morphometrics.png").exists())
            self.assertTrue((fig_base / "fig04c_ecoregion_and_batch_audit.png").exists())
            self.assertEqual(list(fig_base.glob("*.pdf")), [])
            sp_row = grp_df[grp_df["grouping_factor"] == "species_standardized"].iloc[0]
            self.assertGreater(float(sp_row["approx_F"]), 0.0)
            self.assertGreater(float(rg_row["approx_F"]), 0.0)

            # 4. Perturb Tier_3_Bronze drastically and assert strict invariance across allometry, GMM BIC, CDA, and MANOVA
            tbl_pert, rep_pert, _ = _run_case("perturbed_bronze", perturb_bronze=True)
            flags2 = pd.read_csv(tbl_pert / "morphometric_flags.csv")
            bic2 = pd.read_csv(rep_pert / "gmm_bayes_factors_summary.csv")
            allom2 = pd.read_csv(rep_pert / "allometry.csv")

            np.testing.assert_allclose(
                allom1.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                allom2.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                atol=1e-6,
                err_msg="Tier_3_Bronze perturbation must have zero influence on reference allometry ANOVA.",
            )
            np.testing.assert_allclose(
                bic1.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                bic2.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                atol=1e-6,
                err_msg="Tier_3_Bronze perturbation must have zero influence on reference GMM BIC table.",
            )
            act1 = flags1[flags1["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])].sort_values("catalogNumber")
            act2 = flags2[flags2["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])].sort_values("catalogNumber")
            for col in ["can1", "can2", "cda_posterior_prob", "PC1", "PC2"]:
                np.testing.assert_allclose(
                    act1[col].to_numpy(),
                    act2[col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter reference {col}.",
                )
            pd.testing.assert_frame_equal(grp_df, pd.read_csv(rep_pert / "manova_group_significance.csv"))
            pd.testing.assert_frame_equal(pw_df, pd.read_csv(rep_pert / "manova_pairwise_posthoc.csv"))

    def test_dynamic_subspace_beyond_five_pcs_and_multipc_bh_batch_audit(self):
        """Verify Task 1C dynamic subspace selection (>5 PCs when informative), selected_morph_features export, and multi-feature BH institutional batch audit."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(20260431)
        rows = []
        vouch_rows = []
        # 2 taxa x 24 Gold anchors each (N_min = 24 -> Williams & Titus ceiling = 8 >= 6)
        taxa = ["Packera cymbalaria", "Packera hyperborea"]
        for s_idx, sp in enumerate(taxa):
            for j in range(24):
                cat = f"HI_DIM_{s_idx}_{j:02d}"
                cs = float(18.0 + 0.3 * j)
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": "Tier_1_Gold",
                    "foliar_variance": 0.015,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": float(np.log(cs)),
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.55 + rng.normal(0.0, 0.02)),
                    "solidity": 0.91,
                    "archetype_score": 0.95,
                }
                # 6 orthogonal high-variance signal PCs (PC1..PC6) + 2 low-noise PCs (PC7..PC8)
                for p in range(1, 7):
                    val = float(rng.normal(2.5 if (s_idx == 0 and p == 1) else (-2.5 if (s_idx == 1 and p == 1) else 0.0), 4.0))
                    row[f"PC{p}"] = val
                for p in range(7, 9):
                    val = float(rng.normal(0.0, 0.03))
                    row[f"PC{p}"] = val
                for k in range(1, 13):
                    row[f"A{k}"] = float(0.05 / k + rng.normal(0.0, 0.002))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float(0.04 / k + rng.normal(0.0, 0.002))
                rows.append(row)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "determiner_tier": "Tier_1_Gold",
                    "institutionCode": "NCU" if j % 2 == 0 else "ALA",
                    "stateProvince": "Alaska",
                    "decimalLatitude": 63.4,
                    "decimalLongitude": -148.9,
                    "regional_group": "Alaska_Range",
                    "doy": 210,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            batch_csv = tmp_p / "batch_audit.csv"
            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(tmp_p / "bic.csv"),
                    "--output-plot", str(tmp_p / "fig04b_gmm_cda_morphometrics.png"),
                    "--output-batch-audit", str(batch_csv),
                    "--allometry-report", str(tmp_p / "allom.csv"),
                    "--allometry-pairwise", str(tmp_p / "pw.csv"),
                    "--allometry-plot", str(tmp_p / "fig04a_foliar_allometry_audit.png"),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 dynamic >5 PC test failed:\n{proc.stderr}\n{proc.stdout}")
            flags_df = pd.read_csv(flags_csv)
            batch_df = pd.read_csv(batch_csv)

            self.assertIn("selected_morph_features", flags_df.columns)
            sel_feats = [t for t in str(flags_df["selected_morph_features"].iloc[0]).split(";") if t]
            self.assertGreaterEqual(len(sel_feats), 6, f"Expected >= 6 dynamically retained morphometric features without max_p=5 cap, got {sel_feats}")

            inst_rows = batch_df[batch_df["Term"] == "institutionCode"]
            self.assertGreaterEqual(len(inst_rows), len(sel_feats))
            self.assertIn("p_value_bh", batch_df.columns)
            self.assertTrue(inst_rows["p_value_bh"].notna().all())


if __name__ == "__main__":
    unittest.main()
