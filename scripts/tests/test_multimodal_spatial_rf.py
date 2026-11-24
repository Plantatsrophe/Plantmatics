#!/usr/bin/env python3
"""
===============================================================================
Script: test_multimodal_spatial_rf.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Core Step 05 Spatial RF integration and CLI contract tests
    (scripts/analysis/05_multimodal_spatial_rf.R):
    - Species-agnostic Packera taxonomic standardization in Python & R
    - Organic PCA-env intraspecific ecotypes & dynamic Spatial RF figure output
    - Composite-key (catalogNumber, plant_individual_id) multi-plant joins
    - Spatial Block-CV OOF RF posterior integration in edaphic_best_fit_taxon
    - End-to-end authority-anchored CLI run, zero synthetic cache, & Step 07 handoff
===============================================================================
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon
from scripts.tests._r_session import rscript_available, run_warm_r_cmd


class TestMultimodalSpatialRF(unittest.TestCase):
    """Core Spatial RF integration and CLI contract suite for Step 05."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"
        cls.util_path = PROJECT_ROOT / "scripts" / "analysis" / "edaphic_fine_grain_utils.R"
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.pca_env_dir = cls._class_tmp_path / "pca_env_plots"
        cls.pca_env_dir.mkdir(parents=True, exist_ok=True)
        cls.pca_out_png = cls.pca_env_dir / "fig05a_ecological_niche_and_rf_importance.png"
        cls.r_unit_proc = None

        if rscript_available():
            r_batch_code = f"""
            if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/05_multimodal_spatial_rf.R")

            # 1. Taxonomic standardization
            local({{
              stopifnot(standardize_packera_taxon("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve") == "Packera streptanthifolia")
              stopifnot(standardize_packera_taxon("P. sanguisorbae (DC.) C.Jeffrey") == "Packera sanguisorbae")
              stopifnot(standardize_packera_taxon("Packera x memmingeri") == "Packera memmingeri")
              cat("STD_TAXON_OK\\n")
            }})

            # 2. Organic PCA-env intraspecific ecotypes & dynamic Spatial RF plots
            local({{
              set.seed(123)
              n <- 28
              df <- data.frame(
                catalogNumber = sprintf("INTRA_%03d", seq_len(n)),
                species_standardized = rep("Packera aurea", n),
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver"), length.out = n),
                latitude = c(seq(34.0, 35.4, length.out = 14), seq(40.0, 41.8, length.out = 14)),
                longitude = c(seq(-83.0, -81.6, length.out = 14), seq(-76.0, -74.2, length.out = 14)),
                coordinateUncertaintyInMeters = rep(200.0, n),
                fine_bedrock_depth_cm = c(rnorm(14, mean = 25, sd = 3), rnorm(14, mean = 175, sd = 5)),
                bio1_temp_mean = c(rnorm(14, mean = 16.5, sd = 0.4), rnorm(14, mean = 8.2, sd = 0.4)),
                fine_soil_ph = rnorm(n, mean = 5.8, sd = 0.05),
                fine_soil_sand = rnorm(n, mean = 45.0, sd = 0.5),
                topo_twi_30m = c(rnorm(14, mean = 4.5, sd = 0.3), rnorm(14, mean = 9.5, sd = 0.3)),
                doy = rep(130, n),
                delta_doy = rnorm(n, mean = 0, sd = 2),
                stringsAsFactors = FALSE
              )

              resolved <- resolve_analysis_group(df, intraspecific_mode = TRUE)
              u_grps <- sort(unique(resolved$analysis_group))
              stopifnot(length(u_grps) == 2)
              stopifnot(!any(grepl("Acidic_Sand_Coastal|Upland_Piedmont_Outcrop", u_grps)))
              stopifnot(all(grepl("^Ecotype_[12]_[A-Za-z0-9_]+_(Low|High)$", u_grps)))

              df_flat <- data.frame(
                catalogNumber = sprintf("FLAT_%02d", 1:8),
                species_standardized = rep("Packera aurea", 8),
                fine_soil_ph = rep(6.0, 8),
                stringsAsFactors = FALSE
              )
              res_flat <- resolve_analysis_group(df_flat, intraspecific_mode = TRUE)
              stopifnot(identical(unique(res_flat$analysis_group), "Packera aurea"))

              env_vars <- c("fine_bedrock_depth_cm", "bio1_temp_mean", "fine_soil_ph", "fine_soil_sand", "topo_twi_30m")
              srf_res <- compute_spatial_rf_mems(resolved, env_vars, n_trees = 50, intraspecific_mode = TRUE)
              imp_env <- srf_res$importance[!grepl("^MEM[0-9]+$", srf_res$importance$Predictor), , drop = FALSE]
              imp_env <- imp_env[order(imp_env$MeanDecreaseAccuracy, decreasing = TRUE), , drop = FALSE]
              stopifnot(nrow(imp_env) >= 1L)
              stopifnot(imp_env$Predictor[1L] %in% c("fine_bedrock_depth_cm", "bio1_temp_mean", "topo_twi_30m"))
              imp_ed_only <- imp_env[imp_env$Group == "[A] Fine-Grain Edaphic & Outcrop Lithology", , drop = FALSE]
              stopifnot(nrow(imp_ed_only) >= 1L && imp_ed_only$Predictor[1L] %in% c("fine_bedrock_depth_cm", "topo_twi_30m"))
              pred_lbl <- format_predictor_labels("fine_bedrock_depth_cm")
              stopifnot(nzchar(pred_lbl) && pred_lbl != "fine_bedrock_depth_cm")

              niche_res <- run_warren_niche_identity_tests(resolved, env_vars, n_perm = 10, intraspecific_mode = TRUE)
              export_spatial_rf_figures(resolved, srf_res, niche_res, "{cls.pca_out_png.as_posix()}")
              cat("ORGANIC_PCA_ENV_AND_PLOTS_OK\\n")
            }})

            # 3. Composite key (catalogNumber, plant_individual_id) multi-plant join & expanded morph/vision columns
            local({{
              vouchers_df <- data.frame(
                catalogNumber = c("NCU_MULTI", "NCU_SINGLE"),
                species_standardized = c("Packera aurea", "Packera obovata"),
                latitude = c(35.6, 36.1),
                longitude = c(-82.5, -79.1),
                stringsAsFactors = FALSE
              )
              morph_df <- data.frame(
                catalogNumber = c("NCU_MULTI", "NCU_MULTI", "NCU_SINGLE"),
                plant_individual_id = c(1L, 2L, 1L),
                morph_tier_source = rep("Basal_Rosette", 3),
                basal_leaf_count = c(3L, 2L, 4L),
                gmm_cluster = c(1L, 2L, 1L),
                gmm_uncertainty = c(0.04, 0.18, 0.03),
                cda_predicted_taxon = c("Packera aurea", "Packera obovata", "Packera obovata"),
                cda_posterior_prob = c(0.94, 0.88, 0.96),
                cda_second_taxon = c("Packera obovata", "Packera aurea", "Packera aurea"),
                cda_second_posterior_prob = c(0.06, 0.12, 0.04),
                cda_posterior_margin = c(0.88, 0.76, 0.92),
                PC1 = c(-1.4, 1.3, 1.5),
                PC2 = c(0.2, -0.3, -0.1),
                petiole_length_mm = c(28.0, 11.5, 10.0),
                blade_length_mm = c(25.0, 29.0, 30.0),
                petiole_to_blade_ratio = c(1.12, 0.3966, 0.3333),
                plant_height_cm = c(42.0, 35.0, 38.0),
                scape_to_basal_leaf_ratio = c(7.9, 8.6, 9.5),
                is_bent_stem = c(FALSE, TRUE, FALSE),
                capitula_count = c(12L, 8L, 10L),
                involucre_height_mm = c(6.5, 5.8, 6.0),
                involucre_width_mm = c(5.2, 4.9, 5.0),
                AR_inv = c(1.25, 1.18, 1.20),
                ecoregion_cda_predicted = c("Blue_Ridge", "Piedmont", "Piedmont"),
                ecoregion_cda_posterior_prob = c(0.84, 0.79, 0.91),
                stringsAsFactors = FALSE
              )
              vis_df <- data.frame(
                catalogNumber = c("NCU_MULTI", "NCU_SINGLE"),
                plant_individual_id = c(0L, 0L),
                vision_predicted_label = c("Packera aurea", "Packera obovata"),
                confidence_predicted_class = c(0.91, 0.95),
                vision_second_prob = c(0.09, 0.05),
                vision_margin = c(0.82, 0.90),
                c_error = c(0.08, 0.04),
                is_label_corrupted = c(FALSE, FALSE),
                vision_ood_score = c(0.21, 0.72),
                vision_is_ood = c(FALSE, TRUE),
                foreground_attention_ratio = c(0.83, 0.51),
                is_artifact_driven = c(FALSE, TRUE),
                surface_orientation = c("adaxial", "abaxial"),
                mean_luminance = c(54.2, 71.8),
                mean_saturation = c(0.32, 0.14),
                tomentum_taxon_discordance = c(FALSE, TRUE),
                stringsAsFactors = FALSE
              )

              m_joined <- left_join_by_voucher_individual(vouchers_df, morph_df, names(morph_df))
              mv_joined <- left_join_by_voucher_individual(m_joined, vis_df, names(vis_df))
              stopifnot(nrow(mv_joined) == 3L)
              multi_sub <- mv_joined[mv_joined$catalogNumber == "NCU_MULTI", ]
              stopifnot(identical(sort(as.integer(multi_sub$plant_individual_id)), c(1L, 2L)))
              stopifnot(all(c("cda_second_taxon", "PC1", "petiole_to_blade_ratio", "plant_height_cm",
                              "AR_inv", "ecoregion_cda_predicted", "vision_ood_score", "vision_is_ood",
                              "foreground_attention_ratio", "is_artifact_driven", "surface_orientation",
                              "tomentum_taxon_discordance") %in% names(mv_joined)))
              stopifnot(all(multi_sub$foreground_attention_ratio == 0.83))
              cat("STEP06_COMPOSITE_HANDOFF_OK\\n")
            }})

            # 4. Spatial Block-CV OOF RF posterior integration in edaphic_best_fit_taxon
            local({{
              set.seed(106)
              n6 <- 24
              sp6 <- rep(c("Packera aurea", "Packera obovata"), each = n6 / 2)
              df6 <- data.frame(
                catalogNumber = sprintf("S6_%02d", seq_len(n6)),
                plant_individual_id = 1L,
                species_standardized = sp6,
                analysis_group = sp6,
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver"), length.out = n6),
                latitude = c(seq(35.0, 36.1, length.out = 12), seq(36.2, 37.3, length.out = 12)),
                longitude = c(seq(-82.5, -81.4, length.out = 12), seq(-80.0, -78.9, length.out = 12)),
                coordinateUncertaintyInMeters = rep(150.0, n6),
                fine_soil_ph = c(rnorm(12, 5.2, 0.15), rnorm(12, 7.1, 0.15)),
                fine_soil_sand = c(rnorm(12, 65.0, 3.0), rnorm(12, 28.0, 3.0)),
                fine_bedrock_depth_cm = c(rnorm(12, 120.0, 8.0), rnorm(12, 35.0, 6.0)),
                ssurgo_rock_outcrop_pct = c(rep(2.0, 12), rep(28.0, 12)),
                topo_twi_30m = c(rnorm(12, 9.2, 0.3), rnorm(12, 4.8, 0.3)),
                doy = rep(120, n6),
                stringsAsFactors = FALSE
              )
              df6$fine_soil_ph[1] <- 7.2
              df6$fine_soil_sand[1] <- 26.0
              df6$fine_bedrock_depth_cm[1] <- 33.0
              df6$topo_twi_30m[1] <- 4.7
              df6$determiner_tier[1] <- "Tier_3_Bronze"

              srf6 <- compute_spatial_rf_mems(
                df6,
                c("fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "topo_twi_30m"),
                n_trees = 50
              )
              stopifnot(!is.null(srf6$oof_prob_matrix))
              stopifnot(nrow(srf6$oof_prob_matrix) == nrow(df6))

              cons_mahal_only <- execute_crossmodal_consensus(df6, intraspecific_mode = FALSE, rf_res = NULL)
              cons_pooled <- execute_crossmodal_consensus(df6, intraspecific_mode = FALSE, rf_res = srf6)
              stopifnot("edaphic_best_fit_taxon" %in% names(cons_pooled))
              stopifnot("edaphic_z_score" %in% names(cons_pooled))
              stopifnot(cons_mahal_only$edaphic_best_fit_taxon[1] == "Packera obovata")
              stopifnot(cons_pooled$edaphic_best_fit_taxon[1] == "Packera obovata")

              df_mid <- df6
              df_mid$fine_soil_ph[2] <- mean(df6$fine_soil_ph)
              df_mid$fine_soil_sand[2] <- mean(df6$fine_soil_sand)
              df_mid$fine_bedrock_depth_cm[2] <- mean(df6$fine_bedrock_depth_cm)
              df_mid$ssurgo_rock_outcrop_pct[2] <- mean(df6$ssurgo_rock_outcrop_pct)
              df_mid$topo_twi_30m[2] <- mean(df6$topo_twi_30m)
              synth_rf <- list(
                oof_prob_matrix = matrix(
                  c(rep(c(0.90, 0.10), 1), c(0.03, 0.97), rep(c(0.90, 0.10), 10), rep(c(0.10, 0.90), 12)),
                  nrow = n6,
                  ncol = 2,
                  byrow = TRUE,
                  dimnames = list(NULL, c("Packera aurea", "Packera obovata"))
                ),
                spatial_block_cv_acc = 0.92
              )
              cons_rf_tiebreak <- execute_crossmodal_consensus(df_mid, intraspecific_mode = FALSE, rf_res = synth_rf)
              stopifnot(cons_rf_tiebreak$edaphic_best_fit_taxon[2] == "Packera obovata")
              cat("STEP06_RF_POSTERIOR_POOLING_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step05_core_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = run_warm_r_cmd(
                ["Rscript", str(r_batch_file)],
                stage_key="step05",
                cwd=PROJECT_ROOT,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"
        self.util_path = PROJECT_ROOT / "scripts" / "analysis" / "edaphic_fine_grain_utils.R"

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 05 R unit suite failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_taxonomic_standardization(self):
        """Verify generic species-agnostic standardization across North American Packera binomials in Python and R."""
        self.assertEqual(standardize_packera_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera anonyma (Alph.Wood) W.A.Weber"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. savannarum"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve"), "Packera streptanthifolia")
        self.assertEqual(standardize_packera_taxon("Packera sanguisorbae (DC.) C.Jeffrey"), "Packera sanguisorbae")
        self.assertEqual(standardize_packera_taxon("Packera x memmingeri"), "Packera memmingeri")
        self.assertIn("STD_TAXON_OK", self._require_batched_r())

    def test_organic_pca_env_composite_join_and_rf_pooling(self):
        """Verify organic PCA-env ecotypes, multi-plant composite join, and Spatial Block-CV OOF RF posterior pooling."""
        out = self._require_batched_r()
        for marker in (
            "ORGANIC_PCA_ENV_AND_PLOTS_OK",
            "STEP06_COMPOSITE_HANDOFF_OK",
            "STEP06_RF_POSTERIOR_POOLING_OK",
        ):
            self.assertIn(marker, out)
        self.assertTrue(self.pca_out_png.exists(), "fig05a_ecological_niche_and_rf_importance.png must be generated.")
        self.assertEqual(list(self.pca_env_dir.glob("*.pdf")), [])

    def test_authority_anchored_spatial_block_cv_zero_synthetic_cache_and_step07_contract(self):
        """Verify authority-anchored Spatial RF training, zero synthetic disk caching under --demo, passive projection, and Step 07 compatibility."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vouchers_csv, morph_csv = tmp / "curated_vouchers.csv", tmp / "morphometric_flags.csv"
            pheno_states_csv = tmp / "voucher_phenological_states.csv"
            out_flags, out_summary = tmp / "multimodal_conflict_flags.csv", tmp / "multimodal_spatial_rf_summary.csv"
            out_ssurgo, out_contingency = tmp / "ssurgo_edaphic_validation.csv", tmp / "micro_edaphic_outcrop_contingency.csv"
            out_pheno_anom, out_pheno_sum = tmp / "phenological_anomalies.csv", tmp / "phenological_anomaly_summary.csv"
            out_plot = tmp / "fig05a_ecological_niche_and_rf_importance.png"
            out_pheno_plot = tmp / "fig05b_hopkins_circular_phenology.png"

            fixture_taxa = ["Packera aurea", "Packera millefolium", "Packera obovata", "Packera schweinitziana"]
            rows = []
            for idx in range(32):
                sp = fixture_taxa[idx % len(fixture_taxa)]
                lat, lon = 33.0 + (idx % 8) * 1.05, -84.0 + (idx % 6) * 0.85
                doy = int(round(-30.0 + 3.9 * lat))
                tier = "Tier_1_Gold" if idx % 4 in (0, 1) else ("Tier_2_Silver" if idx % 4 == 2 else "Tier_3_Bronze")
                coord_unc = 16000.0 if idx == 31 else float(120 + (idx % 5) * 150)
                rows.append({
                    "catalogNumber": f"ANC_{idx:03d}", "institutionCode": "NCU", "scientificName": sp,
                    "species_raw": sp, "determiner_tier": tier, "latitude": lat, "longitude": lon,
                    "decimalLatitude": lat, "decimalLongitude": lon, "coordinateUncertaintyInMeters": coord_unc,
                    "eventDate": f"2022-04-{(idx % 20) + 1:02d}", "year": 2022, "month": 4,
                    "day": (idx % 20) + 1, "doy": doy, "reproductiveCondition": "fl",
                    "regional_group": "Piedmont_Granite_Flatrocks" if sp == "Packera anonyma" else "Southeastern_Coastal_Plain",
                    "is_primary_duplicate": True,
                })
            pd.DataFrame(rows).to_csv(vouchers_csv, index=False)
            pd.DataFrame([
                {
                    "catalogNumber": r["catalogNumber"], "plant_individual_id": 1,
                    "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3,
                    "cda_predicted_taxon": r["species_raw"], "cda_posterior_prob": 0.92,
                    "can1": 1.1, "can2": -0.3, "can3": 0.2, "gmm_cluster": (i % 4) + 1,
                    "gmm_uncertainty": 0.05, "morphometric_outlier_flag": False, "misidentification_flag": False,
                }
                for i, r in enumerate(rows)
            ]).to_csv(morph_csv, index=False)
            pd.DataFrame([
                {
                    "catalogNumber": r["catalogNumber"], "total_capitula": 4, "n_anthesis": 3,
                    "n_fruiting": 1, "n_bud": 0, "voucher_phenological_state": "anthesis",
                    "dominant_pappus_ratio": 0.10,
                }
                for r in rows
            ]).to_csv(pheno_states_csv, index=False)

            cmd = [
                "Rscript", str(self.r_script_path),
                "--vouchers", str(vouchers_csv), "--morphometrics", str(morph_csv),
                "--pheno-states", str(pheno_states_csv), "--env-dir", str(tmp / "env"),
                "--ssurgo-cache", str(tmp / "ssurgo_cache"), "--ssurgo-table", str(out_ssurgo),
                "--contingency-table", str(out_contingency), "--pheno-anomalies", str(out_pheno_anom),
                "--pheno-summary", str(out_pheno_sum), "--pheno-plot", str(out_pheno_plot),
                "--output-flags", str(out_flags), "--output-summary", str(out_summary),
                "--output-plot", str(out_plot), "--permutations", "10", "--n-trees", "50", "--demo",
            ]
            res = run_warm_r_cmd(cmd, stage_key="step05", cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(res.returncode, 0, f"Step 05 Rscript failed: {res.stderr}\n{res.stdout}")

            cache_dir = tmp / "ssurgo_cache"
            if cache_dir.exists():
                cached_files = list(cache_dir.glob("*.rds")) + list(cache_dir.glob("*.csv"))
                self.assertEqual(cached_files, [], "Synthetic SSURGO records must never be cached via saveRDS/write_csv.")

            flags_df = pd.read_csv(out_flags)
            self.assertEqual(len(flags_df), 32)
            for col in (
                "multimodal_conflict_flag", "triage_category", "edaphic_best_fit_taxon",
                "rf_predicted_taxon", "rf_predicted_prob", "rf_given_taxon_prob",
                "rf_prediction_source", "spatial_block_cv_acc", "sgmc_lithology_class",
                "fine_bedrock_depth_cm", "fine_cfvo_pct", "ilr_texture_1", "ilr_texture_2",
                "topo_hli_30m", "n_observed_modalities", "ph_q10_buffer", "ph_q90_buffer",
                "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer",
            ):
                self.assertIn(col, flags_df.columns, f"Missing column '{col}' in multimodal_conflict_flags.csv")

            gold_silver_rows = flags_df[
                flags_df["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])
                & (~flags_df["high_coord_uncertainty_flag"])
            ]
            bronze_rows = flags_df[
                (flags_df["determiner_tier"] == "Tier_3_Bronze")
                & (~flags_df["high_coord_uncertainty_flag"])
            ]
            high_unc_rows = flags_df[flags_df["catalogNumber"] == "ANC_031"]
            self.assertTrue((gold_silver_rows["rf_prediction_source"] == "Spatial_Block_OOF_Anchor").all())
            self.assertTrue((bronze_rows["rf_prediction_source"] == "Passive_Projection_Bronze").all())
            self.assertTrue((high_unc_rows["rf_prediction_source"] == "Passive_Projection_High_Uncertainty").all())

            sum_df = pd.read_csv(out_summary)
            for col in ("D_clim", "I_clim", "D_edaph", "I_edaph", "D_joint", "I_joint", "Schoener_D", "Warren_D", "Warren_I"):
                self.assertIn(col, sum_df.columns, f"Missing column '{col}' in multimodal_spatial_rf_summary.csv")

            step07_script = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"
            out_queue, out_report = tmp / "triage_queue.csv", tmp / "revision_summary.md"
            out_triage_plot = tmp / "fig07a_integrative_taxonomic_synthesis.png"
            cmd07 = [
                "Rscript", str(step07_script),
                "--vouchers", str(vouchers_csv), "--morphometrics", str(morph_csv),
                "--multimodal-flags", str(out_flags), "--niche-summary", str(out_summary),
                "--output-queue", str(out_queue), "--output-report", str(out_report),
                "--output-plot", str(out_triage_plot),
            ]
            res07 = run_warm_r_cmd(cmd07, stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(res07.returncode, 0, f"Step 07 execution failed: {res07.stderr}\n{res07.stdout}")
            self.assertTrue(out_queue.exists(), "Step 07 triage_queue.csv must be generated from Step 05 outputs.")
            self.assertEqual(len(pd.read_csv(out_queue)), 32)


if __name__ == "__main__":
    unittest.main()
