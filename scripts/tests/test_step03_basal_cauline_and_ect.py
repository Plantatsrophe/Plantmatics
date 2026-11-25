#!/usr/bin/env python3
"""Warm-R tests for Step 03 tip-up alignment, basal/cauline tracks, multi-island immunity, and ECT weights."""

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


class TestStep03BasalCaulineAndECT(unittest.TestCase):
    """Warm-R tests for Step 03 tip-up alignment, basal/cauline tracks, multi-island immunity, and ECT weights."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "03_fourier_extractor.R")

    def test_vertical_tip_up_contour_horizontal_xax_alignment_before_rm_asym(self) -> None:
        """
        Verify align_contour_horizontal_xax in fourier_parallel_utils.R and 03_fourier_extractor.R
        rotates a vertically oriented tip-up symmetric leaf contour onto the horizontal X-axis
        before Momocs::rm_asym() so pct_asymmetry < 0.05 and is_contour_outlier == False.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = []
            manifest_rows = []
            for idx in range(6):
                cat = f"VERT_TIPUP_{idx:02d}"
                sp = "Packera dubia" if idx < 3 else "Packera anonyma"
                # Vertically oriented tip-up bilaterally symmetric leaf:
                # point 0 at lamina_base (x=250, y=340), point 60 at lamina_tip (x=250, y=160)
                # lateral width along X-axis symmetric around x=250
                a_len = 90.0 + idx * 1.5
                b_wid = 36.0 + (idx % 2) * 3.0
                y = 250.0 + a_len * np.cos(t)
                x = 250.0 + b_wid * np.sin(t) * (1.0 - 0.18 * np.cos(t))
                cpath = contours_dir / f"{cat}_leaf1.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": 1,
                    "point_index": np.arange(120),
                    "x": x,
                    "y": y,
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
                    "solidity": 0.94,
                    "archetype_score": 0.96,
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
            self.assertEqual(proc.returncode, 0, f"Step 03 vertical tip-up test failed:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            self.assertEqual(len(ind_df), 6)
            self.assertTrue(
                (ind_df["pct_asymmetry"].astype(float) < 0.05).all(),
                f"Expected pct_asymmetry < 0.05 on vertically oriented symmetric leaves, got {ind_df['pct_asymmetry'].tolist()}",
            )
            outlier_flags = ind_df["is_contour_outlier"].astype(str).str.upper()
            self.assertTrue(
                (outlier_flags == "FALSE").all(),
                f"Expected zero contour outliers on vertically oriented symmetric leaves, got {outlier_flags.tolist()}",
            )

    def test_step03_preserves_basal_and_cauline_tracks_with_pure_basal_otu_medians(self) -> None:
        """
        Verify that when a specimen has both leaf_zone = "basal" and leaf_zone = "cauline" leaves in Step 03:
          - Both basal and cauline leaves are preserved in leaf_efa_harmonics_ind.csv with their leaf_zone tags.
          - The specimen's primary OTU harmonics (A1..D_H, solidity) in specimen_efa_medians.csv match the
            pure basal leaf outlines (w_cauline = 0.0, uncontaminated by the cauline leaf), while the separate
            cauline track columns (cauline_leaf_count, cauline_solidity_med, cauline_dissection_index_med,
            basal_vs_cauline_length_ratio, heteroblasty_efa_distance) are populated.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_mixed_csv = tmp_p / "manifest_mixed.csv"
            manifest_basal_only_csv = tmp_p / "manifest_basal_only.csv"

            def write_contour(cat: str, lid: int, r_vec: np.ndarray) -> Path:
                x = 250.0 + r_vec * np.cos(t)
                y = 250.0 + r_vec * np.sin(t)
                p = contours_dir / f"{cat}_leaf{lid}.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": lid,
                    "point_index": np.arange(120),
                    "x": x,
                    "y": y,
                }).to_csv(p, index=False)
                return p

            # Basal rosette outline (broad ovate, solidity ~ 0.94) vs dissected cauline outline (deep lobes, solidity ~ 0.48)
            r_basal_1 = 80.0 * (1.0 + 0.12 * np.cos(t))
            r_basal_2 = 82.0 * (1.0 + 0.14 * np.cos(t))
            r_cauline = 55.0 * (1.0 + 0.46 * np.cos(6 * t))

            vouchers = [
                {"catalogNumber": "V_HETERO_1", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF_2", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF_3", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF_4", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            p_h1_b1 = write_contour("V_HETERO_1", 1, r_basal_1)
            p_h1_b2 = write_contour("V_HETERO_1", 2, r_basal_2)
            p_h1_c3 = write_contour("V_HETERO_1", 3, r_cauline)

            basal_rows = [
                {
                    "catalogNumber": "V_HETERO_1",
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "leaf_zone": "basal",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.94,
                    "dissection_index": 1.08,
                    "blade_length_mm": 42.0,
                    "petiole_length_mm": 28.0,
                    "petiole_to_blade_ratio": 0.67,
                    "cy_rel_plant_bbox": 0.86,
                    "cy_rel_sheet": 0.80,
                    "area_mm2": 420.0,
                    "archetype_score": 0.95,
                    "contour_path": str(p_h1_b1),
                },
                {
                    "catalogNumber": "V_HETERO_1",
                    "plant_individual_id": 0,
                    "leaf_id": 2,
                    "leaf_zone": "basal",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.92,
                    "dissection_index": 1.10,
                    "blade_length_mm": 40.0,
                    "petiole_length_mm": 25.0,
                    "petiole_to_blade_ratio": 0.625,
                    "cy_rel_plant_bbox": 0.84,
                    "cy_rel_sheet": 0.78,
                    "area_mm2": 400.0,
                    "archetype_score": 0.95,
                    "contour_path": str(p_h1_b2),
                },
            ]
            for idx_ref, cat_ref in enumerate(["V_REF_2", "V_REF_3", "V_REF_4"], start=2):
                p_ref = write_contour(cat_ref, 1, 78.0 * (1.0 + (0.10 + 0.03 * idx_ref) * np.cos(2 * t)))
                basal_rows.append({
                    "catalogNumber": cat_ref,
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "leaf_zone": "basal",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.91,
                    "dissection_index": 1.12,
                    "blade_length_mm": 38.0,
                    "petiole_length_mm": 20.0,
                    "petiole_to_blade_ratio": 0.53,
                    "cy_rel_plant_bbox": 0.82,
                    "cy_rel_sheet": 0.76,
                    "area_mm2": 360.0,
                    "archetype_score": 0.94,
                    "contour_path": str(p_ref),
                })

            cauline_row = {
                "catalogNumber": "V_HETERO_1",
                "plant_individual_id": 0,
                "leaf_id": 3,
                "leaf_zone": "cauline",
                "pixels_per_mm": 10.0,
                "is_calibrated_mm": True,
                "solidity": 0.48,
                "dissection_index": 2.35,
                "blade_length_mm": 20.5,
                "petiole_length_mm": 0.8,
                "petiole_to_blade_ratio": 0.039,
                "cy_rel_plant_bbox": 0.20,
                "cy_rel_sheet": 0.22,
                "area_mm2": 110.0,
                "archetype_score": 0.90,
                "contour_path": str(p_h1_c3),
            }

            pd.DataFrame(basal_rows).to_csv(manifest_basal_only_csv, index=False)
            pd.DataFrame(basal_rows + [cauline_row]).to_csv(manifest_mixed_csv, index=False)

            figs_b_dir = tmp_p / "figs_b"
            figs_m_dir = tmp_p / "figs_m"
            figs_b_dir.mkdir(parents=True, exist_ok=True)
            figs_m_dir.mkdir(parents=True, exist_ok=True)

            out_spec_basal = tmp_p / "spec_basal_only.csv"
            out_ind_basal = tmp_p / "ind_basal_only.csv"
            proc_b = run_warm_r_cmd(
                ["Rscript", str(self.r_script_path), "--contours-dir", str(contours_dir), "--masks-dir", str(tmp_p / "masks"), "--vouchers", str(vouchers_csv), "--manifest", str(manifest_basal_only_csv), "--output", str(out_spec_basal), "--output-individual", str(out_ind_basal), "--figures-dir", str(figs_b_dir)],
                stage_key="step03", cwd=PROJECT_ROOT, skip_figures=False,
            )
            self.assertEqual(proc_b.returncode, 0, f"Basal-only Step 03 failed:\n{proc_b.stderr}\n{proc_b.stdout}")

            out_spec_mix = tmp_p / "spec_mixed.csv"
            out_ind_mix = tmp_p / "ind_mixed.csv"
            proc_m = run_warm_r_cmd(
                ["Rscript", str(self.r_script_path), "--contours-dir", str(contours_dir), "--masks-dir", str(tmp_p / "masks"), "--vouchers", str(vouchers_csv), "--manifest", str(manifest_mixed_csv), "--output", str(out_spec_mix), "--output-individual", str(out_ind_mix), "--figures-dir", str(figs_m_dir)],
                stage_key="step03", cwd=PROJECT_ROOT, skip_figures=False,
            )
            self.assertEqual(proc_m.returncode, 0, f"Mixed basal+cauline Step 03 failed:\n{proc_m.stderr}\n{proc_m.stdout}")

            # Verify Step 03 exports canonical 300-DPI PNG plates only (zero .pdf files)
            self.assertTrue((figs_b_dir / "fig03a_efa_morphospace_and_harmonics.png").exists())
            self.assertTrue((figs_b_dir / "fig03b_topological_ect_and_heteroblasty.png").exists())
            self.assertEqual(list(figs_b_dir.glob("*.pdf")), [])
            self.assertTrue((figs_m_dir / "fig03a_efa_morphospace_and_harmonics.png").exists())
            self.assertTrue((figs_m_dir / "fig03b_topological_ect_and_heteroblasty.png").exists())
            self.assertEqual(list(figs_m_dir.glob("*.pdf")), [])

            ind_mix = pd.read_csv(out_ind_mix)
            spec_basal = pd.read_csv(out_spec_basal)
            spec_mix = pd.read_csv(out_spec_mix)

            # 1. Both basal and cauline leaves are preserved in leaf_efa_harmonics_ind.csv with leaf_zone tags
            h1_ind = ind_mix.loc[ind_mix["catalogNumber"] == "V_HETERO_1"].sort_values("leaf_id")
            self.assertEqual(h1_ind["leaf_id"].astype(int).tolist(), [1, 2, 3])
            self.assertEqual(
                h1_ind["leaf_zone"].astype(str).str.lower().tolist(),
                ["basal", "basal", "cauline"],
            )

            # 2. Primary OTU harmonics (A1..D6, solidity) in spec_mix match pure basal leaf outlines (w_cauline = 0.0)
            row_b = spec_basal.loc[spec_basal["catalogNumber"] == "V_HETERO_1"].iloc[0]
            row_m = spec_mix.loc[spec_mix["catalogNumber"] == "V_HETERO_1"].iloc[0]
            self.assertEqual(int(row_m["basal_leaf_count"]), 2)
            self.assertAlmostEqual(float(row_m["solidity"]), float(row_b["solidity"]), places=4)
            for h_col in ["A1", "A2", "A3", "D1", "D2", "D3"]:
                if h_col in row_m and h_col in row_b:
                    self.assertAlmostEqual(
                        float(row_m[h_col]),
                        float(row_b[h_col]),
                        places=5,
                        msg=f"Primary OTU harmonic {h_col} must match pure basal leaves without cauline contamination.",
                    )

            # 3. Separate cauline track columns (EFA & ECT) are populated on the mixed specimen
            self.assertEqual(int(row_m["cauline_leaf_count"]), 1)
            self.assertAlmostEqual(float(row_m["cauline_solidity_med"]), 0.48, places=2)
            self.assertGreater(float(row_m["cauline_dissection_index_med"]), 1.8)
            self.assertGreater(float(row_m["basal_vs_cauline_length_ratio"]), 1.5)
            self.assertGreater(float(row_m["heteroblasty_efa_distance"]), 0.01)
            self.assertIn("ECT_PC1", ind_mix.columns)
            self.assertIn("ect_lobe_complexity", ind_mix.columns)
            self.assertIn("ECT_PC1", spec_mix.columns)
            self.assertIn("Cauline_ECT_PC1", spec_mix.columns)
            self.assertIn("ect_lobe_complexity", spec_mix.columns)
            self.assertIn("cauline_ect_lobe_complexity", spec_mix.columns)
            self.assertIn("heteroblasty_ect_distance", spec_mix.columns)
            self.assertTrue(np.isfinite(float(row_m["ECT_PC1"])))
            self.assertTrue(np.isfinite(float(row_m["ect_lobe_complexity"])))
            self.assertTrue(np.isfinite(float(row_m["cauline_ect_lobe_complexity"])))
            self.assertGreater(float(row_m["heteroblasty_ect_distance"]), 0.0)

    def test_step03_cordate_upright_multi_island_rm_asym_type_immunity_and_ect_weights(self) -> None:
        """
        Live Rscript integration test verifying:
          1. Cordate wider-than-long mask_is_upright == TRUE preserves exact 0.000-deg horizontal midrib.
          2. Multi-island (n_lamina_components > 1) sibling leaves do not disable Momocs::rm_asym() on closed leaves,
             and multi-island-only OTUs retain PC1 == NA rather than 0.0.
          3. Nomenclatural type leaves (is_type_specimen == TRUE) are immune to is_contour_outlier.
          4. compute_cycle_ect_r uses clockwise u_right = (uy, -ux) and ECT_PC1..3 excludes Tier_3_Bronze (wt=0).
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            r_test = tmp_p / "run_step03_zero_trust.R"
            r_test.write_text(
                f"""
                if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
                source("{PROJECT_ROOT / 'scripts/morphometrics/pub_figures_utils.R'}")
                source("{PROJECT_ROOT / 'scripts/morphometrics/authority_weighting_utils.R'}")
                source("{PROJECT_ROOT / 'scripts/morphometrics/dissected_leaf_utils.R'}")
                source("{PROJECT_ROOT / 'scripts/morphometrics/fourier_parallel_utils.R'}")
                source("{PROJECT_ROOT / 'scripts/morphometrics/03_fourier_extractor.R'}")

                # 1. Cordate wider-than-long upright contour (basal auricles 45% below sinus)
                th <- seq(0, 2 * pi * (1 - 1/120), length.out = 120)
                aur <- exp(-((abs(((th + pi) %% (2 * pi)) - pi) - 0.72)^2) / 0.10)
                coo_cord <- cbind(250 - 115 * sin(th), 235 + 50 * cos(th) + 62 * aur)
                al_cord <- align_contour_horizontal_xax(coo_cord, is_upright = TRUE)
                stopifnot(isTRUE(attr(al_cord, "lm3_from_tipup")))
                stopifnot(abs(al_cord[1L, 2L] - al_cord[61L, 2L]) < 1e-9)

                # 2. Clockwise u_right check: Step 02 tip-up contour with 3 lobes strictly on anatomical right (x_img > 250)
                x_img_r <- ifelse(sin(th) >= 0, 250 - 35 * sin(th), 250 - (35 + 55 * (sin(4 * th)^2)) * sin(th))
                y_img_r <- 250 + 90 * cos(th)
                coo_img_r <- cbind(x_img_r, y_img_r)
                al_right <- align_contour_horizontal_xax(coo_img_r, is_upright = TRUE)
                ect_r <- compute_cycle_ect_r(al_right, num_thresh = 16L, num_dirs = 36L)
                stopifnot(!is.null(ect_r), ncol(ect_r) == 36L)
                # Verify exact 0.000-diff parity against Step 02 image-space tip_up_directions (North=(0,-1), East=(1,0)):
                xy_c2 <- sweep(coo_img_r, 2L, colMeans(coo_img_r), "-")
                xy_n2 <- xy_c2 / max(sqrt(rowSums(xy_c2^2)))
                phis2 <- seq(0.0, 2.0 * pi * (1.0 - 1.0 / 36L), length.out = 36L)
                thr2 <- seq(-1.0, 1.0, length.out = 16L)
                nxt2 <- c(2:nrow(xy_n2), 1L)
                ect_step02 <- matrix(0.0, nrow = 16L, ncol = 36L)
                for (jj in seq_len(36L)) {{
                  hv2 <- as.numeric(xy_n2[, 1L] * sin(phis2[jj]) - xy_n2[, 2L] * cos(phis2[jj]))
                  he2 <- pmax(hv2, hv2[nxt2])
                  ect_step02[, jj] <- as.numeric(findInterval(thr2, sort(hv2)) - findInterval(thr2, sort(he2)))
                }}
                stopifnot(max(abs(ect_r - ect_step02)) == 0)
                stopifnot(sum(ect_r[10:14, 10L]) > sum(ect_r[10:14, 28L]))

                # 3. Type immunity + multi-island sibling OutCoe rm_asym invariance
                coe_m <- matrix(0.0, nrow = 12, ncol = 16)
                colnames(coe_m) <- as.vector(outer(c("A", "B", "C", "D"), 1:4, paste0))
                set.seed(11)
                for (i in 1:11) {{
                  coe_m[i, c("A1","D1","A2","D2","A3","D3","A4","D4")] <- c(-1.2 + 0.05*i, 0.55 - 0.02*i, rnorm(6, 0, 0.02))
                  coe_m[i, c("B1","C1","B2","C2")] <- c(0.08, 0.14 * ((-1)^i), 0.03, -0.04)
                }}
                # Row 1: holotype with high roughness/asymmetry -> must stay is_contour_outlier == FALSE
                coe_m[1, c("A4","D4","B1","C1")] <- c(0.45, 0.45, 0.35, 0.35)
                # Row 12: multi-island leaf (NA EFA harmonics)
                coe_m[12, ] <- NA_real_
                fac_t <- data.frame(
                  catalogNumber = paste0("CAT_", 1:12),
                  plant_individual_id = 0L,
                  leaf_id = 1L,
                  scientificName = rep(c("Packera dubia", "Packera aurea"), each = 6L),
                  determiner_tier = c("Tier_1_Gold", rep("Tier_1_Gold", 5), rep("Tier_1_Gold", 4), "Tier_3_Bronze", "Tier_1_Gold"),
                  type_status = c("holotype", rep("", 11)),
                  leaf_zone = "basal",
                  solidity = c(rep(0.88, 11), 0.38),
                  stringsAsFactors = FALSE
                )
                res_s <- extract_symmetric_fourier_pcs(list(coe = coe_m, fac = fac_t), fac_df = fac_t, num_pcs = 3)
                stopifnot(isFALSE(res_s$df$is_contour_outlier[1L]))
                stopifnot(!("C1" %in% colnames(res_s$sym_mat)))
                stopifnot(is.na(res_s$df$PC1[12L]))
                cat("STEP03_ZERO_TRUST_OK\\n")
                """,
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_test)], stage_key="step03", cwd=PROJECT_ROOT)
            self.assertEqual(proc.returncode, 0, f"Step 03 R zero-trust test failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertIn("STEP03_ZERO_TRUST_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
