#!/usr/bin/env python3
"""Warm-R tests for Step 03 EFA QC, adaptive anchor promotion, and dynamic harmonic calibration."""

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


class TestStep03EFAAnchorsAndHarmonics(unittest.TestCase):
    """Warm-R tests for Step 03 EFA QC, adaptive anchor promotion, and dynamic harmonic calibration."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "03_fourier_extractor.R")

    def test_tier1_only_efa_unweighted_median_and_tier3_passive_pca_invariance(self):
        """Verify only Tier 1 intact contours enter EFA, OTU aggregation uses unweighted median, and Tier_3_Bronze has zero influence on reference PCA."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        def write_contour_csv(path: Path, cat: str, lid: int, a: float, b: float, d2_mod: float = 0.0):
            x = a * np.cos(t)
            y = b * np.sin(t) * (1.0 + d2_mod * np.cos(t))
            pd.DataFrame({
                "catalogNumber": cat,
                "leaf_id": lid,
                "point_index": np.arange(120),
                "x": x + 200.0,
                "y": y + 200.0,
            }).to_csv(path, index=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = [
                {"catalogNumber": "G1", "species_raw": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "G2", "species_raw": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "G3", "species_raw": "Packera anonyma", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "S1", "species_raw": "Packera anonyma", "determiner_tier": "Tier_2_Silver"},
                {"catalogNumber": "S2", "species_raw": "Packera paupercula", "determiner_tier": "Tier_2_Silver"},
                {"catalogNumber": "B1", "species_raw": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
                {"catalogNumber": "T2_ONLY", "species_raw": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            write_contour_csv(contours_dir / "G1_leaf1.csv", "G1", 1, 90.0, 35.0, -0.15)
            write_contour_csv(contours_dir / "G1_leaf2.csv", "G1", 2, 90.0, 45.0, 0.00)
            write_contour_csv(contours_dir / "G1_leaf3.csv", "G1", 3, 90.0, 60.0, 0.20)
            write_contour_csv(contours_dir / "G1_leaf4.csv", "G1", 4, 90.0, 85.0, 0.45)

            write_contour_csv(contours_dir / "G2_leaf1.csv", "G2", 1, 85.0, 42.0, -0.10)
            write_contour_csv(contours_dir / "G3_leaf1.csv", "G3", 1, 95.0, 30.0, 0.18)
            write_contour_csv(contours_dir / "S1_leaf1.csv", "S1", 1, 92.0, 32.0, 0.12)
            write_contour_csv(contours_dir / "S2_leaf1.csv", "S2", 1, 100.0, 28.0, -0.22)
            write_contour_csv(contours_dir / "B1_leaf1.csv", "B1", 1, 88.0, 40.0, 0.05)
            write_contour_csv(contours_dir / "T2_ONLY_leaf1.csv", "T2_ONLY", 1, 80.0, 50.0, 0.10)
            write_contour_csv(contours_dir / "EXCL_PLAT_leaf1.csv", "EXCL_PLAT", 1, 85.0, 45.0, 0.08)

            manifest_rows = [
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.95, "archetype_score": 0.95, "contour_path": str(contours_dir / "G1_leaf1.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 2, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.90, "archetype_score": 0.90, "contour_path": str(contours_dir / "G1_leaf2.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 3, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.85, "archetype_score": 0.85, "contour_path": str(contours_dir / "G1_leaf3.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 4, "efa_qc_pass": False, "qc_flag_reason": "UNRESOLVED_CLUMP", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.92, "archetype_score": 0.40, "contour_path": str(contours_dir / "G1_leaf4.csv")},
                {"catalogNumber": "G2", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.93, "archetype_score": 0.93, "contour_path": str(contours_dir / "G2_leaf1.csv")},
                {"catalogNumber": "G3", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.91, "archetype_score": 0.91, "contour_path": str(contours_dir / "G3_leaf1.csv")},
                {"catalogNumber": "S1", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.92, "archetype_score": 0.92, "contour_path": str(contours_dir / "S1_leaf1.csv")},
                {"catalogNumber": "S2", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.90, "archetype_score": 0.90, "contour_path": str(contours_dir / "S2_leaf1.csv")},
                {"catalogNumber": "B1", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.91, "archetype_score": 0.91, "contour_path": str(contours_dir / "B1_leaf1.csv")},
                {"catalogNumber": "T2_ONLY", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": False, "qc_flag_reason": "UNRESOLVED_CLUMP", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.89, "archetype_score": 0.35, "contour_path": str(contours_dir / "T2_ONLY_leaf1.csv")},
                {"catalogNumber": "EXCL_PLAT", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.94, "archetype_score": 0.94, "contour_path": str(contours_dir / "EXCL_PLAT_leaf1.csv")},
            ]
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            def run_step03(tag: str):
                out_spec = tmp_p / f"leaf_efa_harmonics_{tag}.csv"
                out_ind = tmp_p / f"leaf_efa_harmonics_ind_{tag}.csv"
                proc = run_warm_r_cmd([
                        "Rscript",
                        str(self.r_script_path),
                        "--contours-dir", str(contours_dir),
                        "--masks-dir", str(tmp_p / "masks"),
                        "--vouchers", str(vouchers_csv),
                        "--manifest", str(manifest_csv),
                        "--output", str(out_spec),
                        "--output-individual", str(out_ind),
                    ], stage_key='step03', cwd=PROJECT_ROOT, skip_figures=True)
                self.assertEqual(proc.returncode, 0, f"Step 03 Rscript failed ({tag}):\n{proc.stderr}\n{proc.stdout}")
                return pd.read_csv(out_spec), pd.read_csv(out_ind)

            spec1, ind1 = run_step03("run1")

            # 1. Only QC-passing contours present in vouchers_df enter leaf_efa_harmonics.csv
            self.assertEqual(len(ind1[ind1["catalogNumber"] == "G1"]), 3, "UNRESOLVED_CLUMP leaf 4 must be excluded by Step 03 EFA QC.")
            self.assertNotIn("T2_ONLY", set(spec1["catalogNumber"]), "UNRESOLVED_CLUMP-only specimen must not enter EFA table.")
            self.assertNotIn("EXCL_PLAT", set(ind1["catalogNumber"]), "Excluded taxon absent from vouchers_df must not be ingested into individual EFA table.")
            self.assertNotIn("EXCL_PLAT", set(spec1["catalogNumber"]), "Excluded taxon absent from vouchers_df must not enter specimen EFA table.")
            pub_utils_content = "\n".join(
                p.read_text(encoding="utf-8")
                for p in (
                    PROJECT_ROOT / "scripts" / "morphometrics" / "pub_figures_utils.R",
                    PROJECT_ROOT / "scripts" / "morphometrics" / "pub_taxon_label_utils.R",
                )
                if p.exists()
            )
            self.assertNotIn("if (length(disp_tx) < 4) disp_tx <- unique(c(disp_tx, focal_4))[1:4]", pub_utils_content)

            for col in ("catalogNumber", "plant_individual_id", "pixels_per_mm", "is_calibrated_mm", "solidity", "archetype_score", "centroid_size", "log_centsize"):
                self.assertIn(col, spec1.columns)
                self.assertIn(col, ind1.columns)

            # 2. Unweighted coordinate-wise median per harmonic coefficient for G1
            g1_ind = ind1[ind1["catalogNumber"] == "G1"]
            g1_spec = spec1[spec1["catalogNumber"] == "G1"].iloc[0]
            for h_col in ["A1", "A2", "A3", "D1", "D2", "D3"]:
                expected_med = float(np.median(g1_ind[h_col].to_numpy()))
                self.assertAlmostEqual(float(g1_spec[h_col]), expected_med, places=5)

            # 3. Tier_3_Bronze has authority_weight == 0.0 and zero influence on reference PCA
            b1_row = spec1[spec1["catalogNumber"] == "B1"].iloc[0]
            self.assertEqual(float(b1_row["authority_weight"]), 0.0)

            # Perturb B1's contour drastically and re-run Step 03
            write_contour_csv(contours_dir / "B1_leaf1.csv", "B1", 1, 130.0, 22.0, -0.38)
            spec2, ind2 = run_step03("run2")

            ref1 = spec1[spec1["catalogNumber"] != "B1"].sort_values("catalogNumber")
            ref2 = spec2[spec2["catalogNumber"] != "B1"].sort_values("catalogNumber")
            for pc_col in ["PC1", "PC2", "PC3"]:
                np.testing.assert_allclose(
                    ref1[pc_col].to_numpy(),
                    ref2[pc_col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter reference {pc_col} scores.",
                )
            ref_ind1 = ind1[ind1["catalogNumber"] != "B1"].sort_values(["catalogNumber", "leaf_id"])
            ref_ind2 = ind2[ind2["catalogNumber"] != "B1"].sort_values(["catalogNumber", "leaf_id"])
            for pc_col in ["PC1", "PC2", "PC3"]:
                np.testing.assert_allclose(
                    ref_ind1[pc_col].to_numpy(),
                    ref_ind2[pc_col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter individual reference {pc_col} scores.",
                )

    def test_zero_gold_silver_adaptive_anchor_promotion_and_step01_step02_handoff(self):
        """Verify Step 03 promotes clean Bronze anchors for 0-Gold/Silver species and preserves Step 01/02 metadata and petiole/dissection traits."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        def write_contour_csv(path: Path, cat: str, lid: int, a: float, b: float, d2_mod: float = 0.0):
            x = a * np.cos(t)
            y = b * np.sin(t) * (1.0 + d2_mod * np.cos(t))
            pd.DataFrame({
                "catalogNumber": cat,
                "leaf_id": lid,
                "point_index": np.arange(120),
                "x": x + 200.0,
                "y": y + 200.0,
            }).to_csv(path, index=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            voucher_defs = [
                ("D_G1", "Packera dubia", "Tier_1_Gold", "Holotype"),
                ("D_G2", "Packera dubia", "Tier_1_Gold", ""),
                ("D_B1", "Packera dubia", "Tier_3_Bronze", ""),
                ("A_S1", "Packera anonyma", "Tier_2_Silver", ""),
                ("A_S2", "Packera anonyma", "Tier_2_Silver", ""),
                # Rare endemic Packera serpenticola with ZERO Gold and ZERO Silver vouchers (only Bronze)
                ("R_B1", "Packera serpenticola", "Tier_3_Bronze", ""),
                ("R_B2", "Packera serpenticola", "Tier_3_Bronze", ""),
                ("R_B3", "Packera serpenticola", "Tier_3_Bronze", ""),
            ]
            vouchers = [
                {
                    "catalogNumber": cat,
                    "taxon": sp,
                    "scientificName": sp,
                    "species_raw": sp,
                    "institutionCode": "NCU",
                    "determiner_tier": tier,
                    "authority_tier": tier,
                    "type_status": t_stat,
                    "stateProvince": "North Carolina",
                    "county": "Macon",
                    "locality": "Whiteside Mountain granitic dome",
                    "decimalLatitude": 35.081,
                    "decimalLongitude": -83.138,
                    "coordinateUncertainty": 120.0,
                    "year": 2021,
                    "month": 5,
                    "day": 14,
                    "eventDate": "2021-05-14",
                }
                for cat, sp, tier, t_stat in voucher_defs
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            specs = [
                ("D_G1", 92.0, 38.0, -0.14, False),
                ("D_G2", 90.0, 40.0, -0.12, False),
                ("D_B1", 89.0, 41.0, -0.10, False),
                ("A_S1", 96.0, 28.0, 0.16, True),
                ("A_S2", 94.0, 30.0, 0.14, True),
                ("R_B1", 72.0, 58.0, 0.28, False),
                ("R_B2", 74.0, 56.0, 0.26, False),
                ("R_B3", 73.0, 57.0, 0.27, False),
            ]
            manifest_rows = []
            for idx, (cat, a_v, b_v, d2_v, is_diss) in enumerate(specs):
                cpath = contours_dir / f"{cat}_leaf1.csv"
                write_contour_csv(cpath, cat, 1, a_v, b_v, d2_v)
                manifest_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.42 if is_diss else 0.93,
                    "archetype_score": 0.94,
                    "petiole_length_mm": 18.5 + idx,
                    "blade_length_mm": 24.0 + idx,
                    "petiole_to_blade_ratio": round((18.5 + idx) / (24.0 + idx), 4),
                    "is_dissected": is_diss,
                    "contour_path": str(cpath),
                })
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            out_spec = tmp_p / "leaf_efa_harmonics.csv"
            out_ind = tmp_p / "leaf_efa_harmonics_ind.csv"
            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--contours-dir", str(contours_dir),
                    "--masks-dir", str(tmp_p / "masks"),
                    "--vouchers", str(vouchers_csv),
                    "--manifest", str(manifest_csv),
                    "--output", str(out_spec),
                    "--output-individual", str(out_ind),
                ], stage_key='step03', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 03 failed on 0-Gold/Silver & handoff run:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            spec_df = pd.read_csv(out_spec)
            self.assertEqual(len(spec_df), 8)

            # 1. Adaptive Anchor Promotion & Empirical Bayes active weight ordering: Gold > Silver > Promoted Bronze > 0
            self.assertEqual(float(spec_df.loc[spec_df["catalogNumber"] == "D_B1", "authority_weight"].iloc[0]), 0.0)
            rare_rows = spec_df[spec_df["catalogNumber"].isin(["R_B1", "R_B2", "R_B3"])]
            self.assertGreaterEqual(int((rare_rows["authority_weight"] > 0.0).sum()), 2)

            w_gold = float(spec_df.loc[spec_df["catalogNumber"] == "D_G1", "authority_weight"].iloc[0])
            w_silver = float(spec_df.loc[spec_df["catalogNumber"] == "A_S1", "authority_weight"].iloc[0])
            w_promo = float(rare_rows.loc[rare_rows["authority_weight"] > 0.0, "authority_weight"].iloc[0])
            self.assertGreater(w_gold, w_silver)
            self.assertGreater(w_silver, w_promo)
            self.assertGreater(w_promo, 0.0)
            active_w = spec_df.loc[spec_df["authority_weight"] > 0.0, "authority_weight"]
            self.assertAlmostEqual(float(active_w.mean()), 1.0, places=4)

            # 2. Step 01 locality/date metadata & Step 02 petiole/dissection trait handoff
            for df_chk in (ind_df, spec_df):
                for req_col in (
                    "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "is_dissected",
                    "locality", "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
                    "year", "month", "day", "eventDate", "type_status",
                ):
                    self.assertIn(req_col, df_chk.columns)
            row0 = spec_df[spec_df["catalogNumber"] == "D_G1"].iloc[0]
            self.assertAlmostEqual(float(row0["petiole_length_mm"]), 18.5, places=3)
            self.assertAlmostEqual(float(row0["blade_length_mm"]), 24.0, places=3)
            self.assertEqual(str(row0["locality"]), "Whiteside Mountain granitic dome")
            self.assertEqual(str(row0["eventDate"]), "2021-05-14")
            row_diss = spec_df[spec_df["catalogNumber"] == "A_S1"].iloc[0]
            self.assertTrue(bool(row_diss["is_dissected"]))

    def test_dynamic_harmonic_calibration_and_dynamic_pc_export(self):
        """Verify Task 1B dynamic H_opt in [6, 24] harmonic calibration and dynamic PC1..PC{k_export} export (>5 PCs when k_95 > 5)."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("Momocs::efourier(out_obj, nb.h = h_opt, norm = FALSE, start = FALSE)", content)
        self.assertIn("efd_harmonic_opt", content)
        self.assertIn("efd_k95_pcs", content)

        if not rscript_available():
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)
        rng = np.random.default_rng(20260430)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = []
            manifest_rows = []
            # Create 18 outlines with 8 orthogonal symmetric harmonic variation modes so k_95 >= 6
            freqs = [2, 3, 4, 5, 6, 7, 8, 9]
            for i in range(18):
                cat = f"DYN_H_{i:02d}"
                sp = "Packera millefolium" if i < 9 else "Packera schweinitziana"
                coefs = rng.normal(0.0, 0.055, size=len(freqs))
                x = 85.0 * np.cos(t) + sum(85.0 * coefs[m] * np.cos(freqs[m] * t) for m in range(len(freqs)))
                y = 42.0 * np.sin(t) + sum(42.0 * coefs[m] * np.sin(freqs[m] * t) for m in range(len(freqs)))
                cpath = contours_dir / f"{cat}_leaf1.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": 1,
                    "point_index": np.arange(120),
                    "x": x + 250.0,
                    "y": y + 250.0,
                }).to_csv(cpath, index=False)
                vouchers.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": "Tier_1_Gold",
                })
                manifest_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.62,
                    "archetype_score": 0.94,
                    "is_dissected": True,
                    "contour_path": str(cpath),
                })

            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            out_spec = tmp_p / "leaf_efa_harmonics.csv"
            out_ind = tmp_p / "leaf_efa_harmonics_ind.csv"
            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--contours-dir", str(contours_dir),
                    "--masks-dir", str(tmp_p / "masks"),
                    "--vouchers", str(vouchers_csv),
                    "--manifest", str(manifest_csv),
                    "--output", str(out_spec),
                    "--output-individual", str(out_ind),
                ], stage_key='step03', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Dynamic H_opt Step 03 failed:\n{proc.stderr}\n{proc.stdout}")
            spec_df = pd.read_csv(out_spec)
            self.assertIn("efd_harmonic_opt", spec_df.columns)
            self.assertIn("efd_k95_pcs", spec_df.columns)
            h_opt_val = int(spec_df["efd_harmonic_opt"].iloc[0])
            k95_val = int(spec_df["efd_k95_pcs"].iloc[0])
            self.assertGreaterEqual(h_opt_val, 6)
            self.assertLessEqual(h_opt_val, 24)
            self.assertGreaterEqual(k95_val, 6, f"Expected k_95 >= 6 on 8-mode outlines, got {k95_val}")
            self.assertIn("PC6", spec_df.columns)


if __name__ == "__main__":
    unittest.main()
