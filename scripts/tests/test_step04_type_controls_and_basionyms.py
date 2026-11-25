#!/usr/bin/env python3
"""Warm-R tests for Step 04 adaptive anchor promotion, type Mahalanobis immunity, and basionym centroids."""

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


class TestStep04TypeControlsAndBasionyms(unittest.TestCase):
    """Warm-R tests for Step 04 adaptive anchor promotion, type Mahalanobis immunity, and basionym centroids."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "04_gmm_morphotools.R")

    def test_zero_gold_species_adaptive_anchor_promotion_and_empirical_bayes_weights(self):
        """Verify multi-species GMM and CDA succeed via Adaptive Anchor Promotion when species have 0 Tier_1_Gold vouchers."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(202604)
        rows = []
        vouch_rows = []

        species_specs = [
            ("Packera dubia", (-5.0, -1.8, 0.4), ["Tier_1_Gold"] * 5 + ["Tier_2_Silver"] * 4 + ["Tier_3_Bronze"] * 3),
            ("Packera anonyma", (5.0, -1.5, -0.4), ["Tier_2_Silver"] * 3 + ["Tier_3_Bronze"] * 7),
            ("Packera serpenticola", (0.0, 4.8, 0.8), ["Tier_3_Bronze"] * 10),
        ]

        for s_idx, (sp, (c1, c2, c3), tiers) in enumerate(species_specs):
            for j, tier in enumerate(tiers):
                cat = f"ADAPT_{s_idx}_{j:02d}"
                pc1 = float(c1 + rng.normal(0.0, 0.16))
                pc2 = float(c2 + rng.normal(0.0, 0.16))
                pc3 = float(c3 + rng.normal(0.0, 0.12))
                if sp == "Packera dubia" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = 5.05, -1.45, -0.38
                elif sp == "Packera anonyma" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = -18.0, 16.0, 9.0
                elif sp == "Packera serpenticola" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = 19.0, -17.0, -8.0

                cs = float(26.0 + 0.9 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": tier,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": pc3,
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.9 + 0.2 * s_idx + rng.normal(0.0, 0.04)),
                    "solidity": 0.90,
                    "archetype_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((pc1 * 0.06) / k + rng.normal(0.0, 0.002))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((pc2 * 0.06) / k + rng.normal(0.0, 0.002))
                rows.append(row)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "determiner_tier": tier,
                    "institutionCode": "NCU" if j % 2 == 0 else "DUKE",
                    "stateProvince": "NC",
                    "county": "Orange",
                    "decimalLatitude": 35.91,
                    "decimalLongitude": -79.05,
                    "regional_group": "Piedmont",
                    "doy": 115 + j,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            cda_pdf = tmp_p / "fig04b_gmm_cda_morphometrics.png"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(tmp_p / "allometry.csv"),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
                    "--allometry-plot", str(tmp_p / "fig04a_foliar_allometry_audit.png"),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 failed with 0-Gold species:\n{proc.stderr}\n{proc.stdout}")

            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), len(rows))
            self.assertEqual(
                set(flags_df["cda_predicted_taxon"].dropna().unique()),
                {"Packera dubia", "Packera anonyma", "Packera serpenticola"},
            )
            for req_col in (
                "min_congener_mahal_d2",
                "min_congener_taxon",
                "open_set_novelty_p",
                "is_confident_label_noise",
            ):
                self.assertIn(req_col, flags_df.columns, f"Missing {req_col} in morphometric_flags.csv")
                self.assertTrue(flags_df[req_col].notna().all(), f"Unexpected NA in {req_col}")

            anon_df = flags_df[flags_df["species_standardized"] == "Packera anonyma"]
            anon_promo = anon_df[(anon_df["determiner_tier"] == "Tier_3_Bronze") & (anon_df["authority_weight"] > 0.0)]
            self.assertGreaterEqual(len(anon_promo), 3, "Expected clean Bronze anchors promoted for 0-Gold Packera anonyma")
            anon_outlier = flags_df[flags_df["catalogNumber"] == "ADAPT_1_09"].iloc[0]
            self.assertEqual(float(anon_outlier["authority_weight"]), 0.0)
            self.assertTrue(bool(anon_outlier["is_passive"]))

            serp_df = flags_df[flags_df["species_standardized"] == "Packera serpenticola"]
            serp_promo = serp_df[serp_df["authority_weight"] > 0.0]
            self.assertGreaterEqual(len(serp_promo), 5, "Expected clean Bronze anchors promoted for 0-Gold/0-Silver species")
            serp_outlier = flags_df[flags_df["catalogNumber"] == "ADAPT_2_09"].iloc[0]
            self.assertEqual(float(serp_outlier["authority_weight"]), 0.0)
            self.assertTrue(bool(serp_outlier["is_passive"]))

            dub_misid = flags_df[flags_df["catalogNumber"] == "ADAPT_0_11"].iloc[0]
            self.assertEqual(str(dub_misid["cda_predicted_taxon"]), "Packera anonyma")
            self.assertEqual(str(dub_misid["discordance_reason"]), "Packera_dubia_Discordant_to_Packera_anonyma")

    def test_nomenclatural_type_mahalanobis_immunity_eb_pinning_and_typical_morphotype_control(self):
        """Verify select_adaptive_reference_anchors() never prunes an intact Holotype/Lectotype during 95% chi^2 FastMCD screening, pins EB weights, and exports nearest_type_morph_taxon and is_typical_morphotype."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        r_code = """
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        options(rgl.useNULL = TRUE)
        source("scripts/morphometrics/04_gmm_morphotools.R")
        set.seed(202611)

        # 1. Verify .identify_type_specimen_vec() ranks primary (2L), secondary (1L), and non-types (0L)
        df_t <- data.frame(
          type_status = c("Holotype", "Isolectotype", "Paratype", "Syntype", "Not a type", NA_character_),
          stringsAsFactors = FALSE
        )
        t_id <- .identify_type_specimen_vec(df_t)
        stopifnot(identical(as.integer(t_id$type_rank), c(2L, 2L, 1L, 1L, 0L, 0L)))
        stopifnot(identical(as.logical(t_id$is_type_specimen), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)))

        # 2. Verify select_adaptive_reference_anchors() never prunes an extreme Holotype during 95% chi^2 MCD screening,
        #    while pruning an equally extreme non-type Tier_1_Gold outlier
        df_anc <- data.frame(
          catalogNumber = c("HOLO_EXTREME", "GOLD_EXTREME_OUTLIER", paste0("GOLD_CORE_", 1:18), paste0("OBO_GOLD_", 1:12)),
          plant_individual_id = 0L,
          species_standardized = c(rep("Packera dubia", 20), rep("Packera obovata", 12)),
          determiner_tier = "Tier_1_Gold",
          type_status = c("Holotype", "None", rep("None", 18), "Lectotype", rep("None", 11)),
          morph_tier_source = "Basal_Rosette",
          PC1 = c(-4.8, 14.5, rnorm(18, -1.8, 0.18), 4.2, rnorm(11, 3.8, 0.18)),
          PC2 = c(3.9, -12.0, rnorm(18, 0.4, 0.18), -1.1, rnorm(11, -0.9, 0.18)),
          PC3 = c(2.6, 9.5, rnorm(18, 0.1, 0.15), -0.4, rnorm(11, -0.3, 0.15)),
          stringsAsFactors = FALSE
        )
        X_anc <- as.matrix(df_anc[, c("PC1", "PC2", "PC3")])
        anc_out <- select_adaptive_reference_anchors(df_anc, feature_matrix = X_anc, taxon_col = "species_standardized")
        # Holotype at row 1 MUST remain active and marked Tier_0_Type_Anchor; non-type extreme outlier at row 2 MUST be pruned
        stopifnot(isTRUE(anc_out$active_mask[1]))
        stopifnot(identical(anc_out$anchor_tier_source[1], "Tier_0_Type_Anchor"))
        stopifnot(!isTRUE(anc_out$active_mask[2]))

        # 3. Verify compute_empirical_bayes_authority_weights() pins primary types to 3.0 in continuous mode,
        #    gives type anchors >= 1.4x higher active weight than ordinary Gold anchors, and supports target_modality="morph"
        w_act <- compute_empirical_bayes_authority_weights(
          df_anc, mode = "active_passive", active_mask = anc_out$active_mask, target_modality = "morph"
        )
        w_cont <- compute_empirical_bayes_authority_weights(
          df_anc, mode = "continuous_prior_posterior", active_mask = anc_out$active_mask, target_modality = "morph"
        )
        stopifnot(w_act[1] > 1.35 * stats::median(w_act[3:20]))
        stopifnot(abs(w_cont[1] - 3.0) < 1e-6)

        # 3b. Verify Kish's N_eff in compute_ledoit_wolf_cov(..., weights = w) and two-sample scaling in compute_type_specimen_morph_control()
        X_lw <- scale(X_anc[anc_out$active_mask, , drop = FALSE], center = TRUE, scale = FALSE)
        S_emp_lw <- stats::cov(X_lw)
        lw_unw <- compute_ledoit_wolf_cov(S_emp_lw, X_centered = X_lw, n_obs = nrow(X_lw))
        w_skew <- c(100.0, rep(0.05, nrow(X_lw) - 1L))
        lw_skew <- compute_ledoit_wolf_cov(S_emp_lw, X_centered = X_lw, n_obs = nrow(X_lw), weights = w_skew)
        stopifnot(attr(lw_skew, "n_eff") < attr(lw_unw, "n_eff"))
        stopifnot(attr(lw_skew, "shrinkage_lambda") >= attr(lw_unw, "shrinkage_lambda"))

        # Verify two-sample (1 + 1/n_eff_type)^(-1) scaling in compute_type_specimen_morph_control()
        df_ts_1 <- df_anc[c(1, 3:10, 21:28), ]
        X_ts_1 <- as.matrix(df_ts_1[, c("PC1", "PC2", "PC3")])
        tc_1 <- compute_type_specimen_morph_control(df_ts_1, X_ts_1, S_reg = diag(3), taxon_col = "species_standardized")
        raw_d2_row2 <- sum((X_ts_1[2, ] - X_ts_1[1, ])^2)
        # Single Holotype (n_eff_type = 1) -> scaled by 1 / (1 + 1/1) = 0.5
        stopifnot(abs(tc_1$type_morph_mahal_d2[2] - round(0.5 * raw_d2_row2, 4)) < 1e-4)
        stopifnot(abs(tc_1$type_morph_mahal_d2[1]) < 1e-8)

        # 4. Verify audit_misidentifications() exports type_morph_mahal_d2, nearest_type_morph_taxon,
        #    type_morph_concordance, is_typical_morphotype, and protects type specimens from routine misidentification
        df_anc$gmm_cluster <- c(2L, 3L, rep(2L, 18), rep(1L, 12))
        df_anc$gmm_uncertainty <- 0.04
        cda_res <- run_canonical_discriminant_analysis(
          df_anc, feature_cols = c("PC1", "PC2", "PC3"), target_taxa = c("Packera dubia", "Packera obovata")
        )
        aud_df <- audit_misidentifications(
          df_anc, cda_res, list(classification = df_anc$gmm_cluster, uncertainty = df_anc$gmm_uncertainty)
        )
        req_cols <- c("is_type_specimen", "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance", "is_typical_morphotype", "has_empirical_type_anchor")
        stopifnot(all(req_cols %in% names(aud_df)))
        # P. dubia Holotype is in gmm_cluster 2 -> cluster 2 is typical morphotype for P. dubia
        stopifnot(isTRUE(aud_df$is_typical_morphotype[1]))
        stopifnot(!isTRUE(aud_df$is_typical_morphotype[2]))
        stopifnot(!isTRUE(aud_df$morphometric_outlier_flag[1]))
        stopifnot(identical(aud_df$nearest_type_morph_taxon[21], "Packera obovata"))
        cat("STEP04_TYPE_CONTROL_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step04.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Step 04 type control R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP04_TYPE_CONTROL_OK", proc.stdout)

    def test_basionym_type_centroids_log1p_capitula_and_blade_cie_l_vestiture(self):
        """Verify basionym-level type centroid preservation (Cineraria dubia vs Senecio alabamensis), sterile-safe log1p_capitulum_count, and eroded-blade CIE L* vestiture integration."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            v_csv = tmp_p / "vouchers.csv"
            h_csv = tmp_p / "harmonics.csv"
            r_csv = tmp_p / "repro.csv"
            f_csv = tmp_p / "flags.csv"
            rep_csv = tmp_p / "bic.csv"
            allom_csv = tmp_p / "allometry.csv"
            batch_csv = tmp_p / "batch.csv"
            plot_png = tmp_p / "plot.png"

            rng = np.random.default_rng(20261112)
            v_rows = []
            h_rows = []
            r_rows = []
            for idx in range(24):
                is_m2 = idx >= 12
                cat = f"DUB_BAS_{idx:02d}"
                if idx == 0:
                    sp_raw = "Cineraria dubia Spreng."
                    t_stat = "Holotype"
                elif idx == 12:
                    sp_raw = "Senecio alabamensis Britton"
                    t_stat = "Holotype"
                else:
                    sp_raw = "Packera dubia"
                    t_stat = ""
                pc1 = float(rng.normal(3.2 if is_m2 else -3.2, 0.22))
                pc2 = float(rng.normal(-1.2 if is_m2 else 1.2, 0.22))
                cap_cnt = 0 if idx in (5, 17) else (int(rng.integers(24, 42)) if is_m2 else int(rng.integers(4, 9)))
                l_norm = float(rng.normal(0.41 if is_m2 else 0.74, 0.025))
                l_sd = float(rng.normal(1.05 if is_m2 else 3.85, 0.18))
                v_rows.append({
                    "catalogNumber": cat,
                    "institutionCode": "NCU" if idx % 2 == 0 else "US",
                    "scientificName": sp_raw,
                    "species_raw": sp_raw,
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": t_stat,
                    "decimalLatitude": 35.2 + 0.08 * idx,
                    "decimalLongitude": -80.5 + 0.08 * idx,
                    "stateProvince": "North Carolina",
                    "regional_group": "Southeastern_Plains" if is_m2 else "Piedmont",
                })
                h_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "morph_tier_source": "Basal_Rosette",
                    "basal_leaf_count": 3,
                    "scale_source": "ruler_calibrated",
                    "is_calibrated_mm": True,
                    "pixels_per_mm": 23.6,
                    "centroid_size": 45.0,
                    "log_centsize": float(np.log10(45.0)),
                    "A1": 1.0, "D1": 0.5, "A2": 0.1 * pc1, "D2": 0.1 * pc2,
                    "A3": 0.02, "D3": 0.02, "A4": 0.01, "D4": 0.01,
                    "PC1": pc1, "PC2": pc2, "PC3": float(rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)), "PC5": float(rng.normal(0.0, 0.10)),
                    "solidity": 0.86,
                    "dissection_index": 1.22,
                    "blade_cie_l_norm": round(l_norm, 4),
                    "blade_cie_l_sd": round(l_sd, 4),
                })
                r_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "capitula_count": cap_cnt,
                    "AR_inv": np.nan if cap_cnt == 0 else 1.28,
                })

            pd.DataFrame(v_rows).to_csv(v_csv, index=False)
            pd.DataFrame(h_rows).to_csv(h_csv, index=False)
            pd.DataFrame(r_rows).to_csv(r_csv, index=False)

            proc = run_warm_r_cmd(
                [
                    "Rscript", str(self.r_script_path),
                    "--vouchers", str(v_csv),
                    "--harmonics", str(h_csv),
                    "--reproductive", str(r_csv),
                    "--output-flags", str(f_csv),
                    "--output-report", str(rep_csv),
                    "--allometry-report", str(allom_csv),
                    "--output-batch-audit", str(batch_csv),
                    "--output-plot", str(plot_png),
                    "--intraspecific-mode",
                    "--focal-taxon", "Packera dubia",
                    "--permutations", "19",
                ],
                stage_key="step04",
                cwd=PROJECT_ROOT,
                skip_figures=True,
            )
            self.assertEqual(proc.returncode, 0, f"Step 04 basionym/vestiture run failed:\n{proc.stderr}\n{proc.stdout}")
            flags = pd.read_csv(f_csv)
            self.assertIn("nearest_type_basionym", flags.columns)
            self.assertIn("nearest_type_basionym_d2", flags.columns)
            self.assertIn("log1p_capitulum_count", flags.columns)
            self.assertIn("blade_cie_l_norm", flags.columns)
            self.assertIn("blade_cie_l_sd", flags.columns)

            sterile_rows = flags[flags["catalogNumber"].isin(["DUB_BAS_05", "DUB_BAS_17"])]
            self.assertTrue(sterile_rows["log1p_capitulum_count"].isna().all())
            fert_rows = flags[~flags["catalogNumber"].isin(["DUB_BAS_05", "DUB_BAS_17"])]
            self.assertTrue(fert_rows["log1p_capitulum_count"].notna().all())

            g1 = flags[flags["catalogNumber"].isin([f"DUB_BAS_{i:02d}" for i in range(12)])]
            g2 = flags[flags["catalogNumber"].isin([f"DUB_BAS_{i:02d}" for i in range(12, 24)])]
            self.assertTrue((g1["nearest_type_basionym"] == "Cineraria dubia").all())
            self.assertTrue((g2["nearest_type_basionym"] == "Senecio alabamensis").all())
            t_cin = flags.loc[flags["catalogNumber"] == "DUB_BAS_00", "nearest_type_basionym_d2"].iloc[0]
            t_ala = flags.loc[flags["catalogNumber"] == "DUB_BAS_12", "nearest_type_basionym_d2"].iloc[0]
            self.assertAlmostEqual(float(t_cin), 0.0, places=4)
            self.assertAlmostEqual(float(t_ala), 0.0, places=4)


if __name__ == "__main__":
    unittest.main()
