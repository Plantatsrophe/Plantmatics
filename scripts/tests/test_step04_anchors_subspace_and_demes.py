#!/usr/bin/env python3
"""Warm-R tests for Step 04 batched R helpers, dissection pruning, subspace rescue, and gathering demes."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


class TestStep04AnchorsSubspaceAndDemes(unittest.TestCase):
    """Warm-R tests for Step 04 batched R helpers, dissection pruning, subspace rescue, and gathering demes."""

    def test_batched_r_helpers_epa_manova_dawid_skene_and_k_full_cda(self):
        """Verify batched R unit checks for EPA Level III helper, single-species MANOVA collinearity, Dawid-Skene EM, and K_full (>3 axes) CDA."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            grp_out = tmp_p / "manova_group_significance.csv"
            pw_out = tmp_p / "manova_pairwise_posthoc.csv"
            r_batch_code = f"""
            if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
            source("scripts/morphometrics/pub_figures_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/morphometrics/04_gmm_morphotools.R")

            # Block 1: EPA Level III Ecoregion spatial helper agreement (offline + on-disk county_norm GeoJSON schema)
            local({{
              res <- infer_epa_level3_ecoregion(
                lat = c(35.9132, 34.2257, 32.5404, 31.6035, NA, NA, NA),
                lon = c(-79.0558, -77.9447, -82.9037, -94.6555, NA, NA, NA),
                state_province = c("NC", "NC", "GA", "TX", "North Carolina", "North Carolina", "NC"),
                county = c("", "", "", "", "Orange County", "Wake County", "")
              )
              off_b <- .build_offline_epa_l3_sf_r()
              eco_gj <- file.path("{tmp_p.as_posix()}", "test_eco_l3.geojson")
              cnty_gj <- file.path("{tmp_p.as_posix()}", "test_us_counties.geojson")
              sf::st_write(off_b$eco_sf, eco_gj, driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
              sf::st_write(off_b$cnty_sf[, c("county_key", "state_name", "county_norm", "centroid_lon", "centroid_lat")], cnty_gj, driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
              res_disk <- infer_epa_level3_ecoregion(
                lat = c(NA, NA),
                lon = c(NA, NA),
                state_province = c("North Carolina", "North Carolina"),
                county = c("Orange County", "Wake County"),
                ecoregion_path = eco_gj,
                county_path = cnty_gj
              )
              stopifnot(identical(res_disk, c("Piedmont", NA_character_)))
              cat("EPA_L3:", paste(ifelse(is.na(res), "NA", res), collapse = "|"), "\\n", sep = "")
            }})

            # Block 2: Direct run_multivariate_manova_audit single-species & collinearity safety
            local({{
              set.seed(17)
              df <- data.frame(
                catalogNumber = sprintf("V%02d", 1:18),
                species_standardized = rep("Packera dubia", 18),
                determiner_tier = rep("Tier_1_Gold", 18),
                morph_tier_source = rep("Basal_Rosette", 18),
                regional_group = c(rep("Piedmont", 10), rep("Blue_Ridge", 2), rep(NA_character_, 6)),
                PC1 = c(rnorm(9, -3.0, 0.3), rnorm(9, 3.0, 0.3)),
                PC2 = c(rnorm(9, 1.5, 0.3), rnorm(9, -1.5, 0.3)),
                stringsAsFactors = FALSE
              )
              df$PC3_collinear <- df$PC1 + df$PC2
              df$PC4_const <- 0.0
              gmm_res <- list(classification = c(rep(1L, 9), rep(2L, 9)))
              # Default include_gmm_clusters = FALSE omits post-selection gmm_cluster inference
              res_def <- run_multivariate_manova_audit(
                df,
                c("PC1", "PC2", "PC3_collinear", "PC4_const"),
                gmm_res,
                output_group_csv = "{grp_out.as_posix()}",
                output_pairwise_csv = "{pw_out.as_posix()}"
              )
              stopifnot(nrow(res_def$group_table) == 0L)
              # Even when descriptive gmm_cluster geometry is requested, p_value and q_value are suppressed (NA)
              res <- run_multivariate_manova_audit(
                df,
                c("PC1", "PC2", "PC3_collinear", "PC4_const"),
                gmm_res,
                output_group_csv = "{grp_out.as_posix()}",
                output_pairwise_csv = "{pw_out.as_posix()}",
                include_gmm_clusters = TRUE
              )
              stopifnot(nrow(res$group_table) == 1L)
              stopifnot(res$group_table$grouping_factor[1] == "gmm_cluster")
              stopifnot(is.na(res$group_table$p_value[1]) && is.na(res$group_table$q_value[1]) && is.na(res$group_table$approx_F[1]))
              stopifnot(all(is.na(res$pairwise_table$p_value)) && all(is.na(res$pairwise_table$q_value)))
              cat("MANOVA_DIRECT_OK\\n")
            }})

            # Block 3: Beta-Binomial Empirical Bayes active_passive & continuous_prior_posterior modes
            local({{
              df_t <- data.frame(
                species_standardized = rep("Packera aurea", 8),
                determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver", "Tier_2_Silver",
                                    "Tier_3_Bronze", "Tier_3_Bronze", "Tier_3_Bronze", "Tier_3_Bronze"),
                is_promoted_bronze = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
                c_error = c(0.02, 0.04, 0.12, 0.14, 0.25, 0.30, 0.55, 0.65),
                stringsAsFactors = FALSE
              )
              w_ap <- compute_empirical_bayes_authority_weights(df_t, mode = "active_passive")
              w_cp <- compute_empirical_bayes_authority_weights(df_t, mode = "continuous_prior_posterior")
              stopifnot(all(w_ap[7:8] == 0.0), abs(mean(w_ap[1:6]) - 1.0) < 1e-5)
              stopifnot(w_ap[1] > w_ap[3], w_ap[3] > w_ap[5], w_ap[5] > 0.0)
              stopifnot(abs(max(w_cp) - 3.0) < 1e-5, min(w_cp) >= 0.5, w_cp[1] > w_cp[3], w_cp[3] > w_cp[5], w_cp[5] > w_cp[7])
              cat("EB_OK\\n")
            }})

            # Block 4: 2-step Dawid-Skene EM determiner shrinkage, no double pruning, and K_full (>3 axes) CDA
            local({{
              set.seed(20261030)
              n_sp <- 10
              df_eb <- data.frame(
                catalogNumber = sprintf("EB_%02d", 1:(2 * n_sp)),
                plant_individual_id = 1L,
                species_standardized = rep(c("Packera aurea", "Packera obovata"), each = n_sp),
                determiner_tier = rep(c(rep("Tier_2_Silver", 8), "Tier_1_Gold", "Tier_3_Bronze"), 2),
                identifiedBy = rep(c(rep("Dr. Reliable Specialist", 4), rep("Noisy Generalist", 4), "T.M. Barkley", "Field Student"), 2),
                is_promoted_bronze = FALSE,
                PC1 = c(rnorm(4, -3.0, 0.15), rnorm(4, +3.0, 0.15), -3.1, -2.9,
                        rnorm(4, +3.0, 0.15), rnorm(4, -3.0, 0.15), +3.1, +2.9),
                PC2 = rnorm(2 * n_sp, 0.0, 0.15),
                vision_predicted_label = c(rep("Packera aurea", 4), rep("Packera obovata", 4), "Packera aurea", "Packera aurea",
                                           rep("Packera obovata", 4), rep("Packera aurea", 4), "Packera obovata", "Packera obovata"),
                c_error = c(rep(0.05, 4), rep(0.92, 4), 0.03, 0.10, rep(0.05, 4), rep(0.92, 4), 0.03, 0.10),
                stringsAsFactors = FALSE
              )

              w_ap <- compute_empirical_bayes_authority_weights(df_eb, mode = "active_passive")
              stopifnot(length(w_ap) == nrow(df_eb))
              stopifnot(w_ap[10] == 0.0, w_ap[20] == 0.0)
              stopifnot(abs(mean(w_ap[w_ap > 0]) - 1.0) < 1e-5)
              stopifnot(mean(w_ap[df_eb$identifiedBy == "Dr. Reliable Specialist"]) >
                        mean(w_ap[df_eb$identifiedBy == "Noisy Generalist"]))

              w_cp <- compute_empirical_bayes_authority_weights(df_eb, mode = "continuous_prior_posterior")
              stopifnot(all(w_cp >= 0.5 - 1e-9 & w_cp <= 3.0 + 1e-9))
              stopifnot(mean(w_cp[df_eb$identifiedBy == "Dr. Reliable Specialist"]) >
                        mean(w_cp[df_eb$identifiedBy == "Noisy Generalist"]))

              set.seed(104)
              taxa5 <- c("Packera spA", "Packera spB", "Packera spC", "Packera spD", "Packera spE")
              centers5 <- rbind(
                c(-4.0,  0.0,  0.0,  0.0),
                c( 4.0,  0.0,  0.0,  0.0),
                c( 0.0,  4.0,  0.0,  0.0),
                c( 0.0, -2.0,  3.0, -4.5),
                c( 0.0, -2.0,  3.0, +4.5)
              )
              rows5 <- list()
              idx <- 1L
              for (g_i in seq_along(taxa5)) {{
                for (j in 1:14) {{
                  rows5[[idx]] <- data.frame(
                    catalogNumber = sprintf("G%d_%02d", g_i, j),
                    plant_individual_id = 1L,
                    species_standardized = taxa5[g_i],
                    determiner_tier = "Tier_1_Gold",
                    morph_tier_source = "Basal_Rosette",
                    PC1 = centers5[g_i, 1] + rnorm(1, 0, 0.18),
                    PC2 = centers5[g_i, 2] + rnorm(1, 0, 0.18),
                    PC3 = centers5[g_i, 3] + rnorm(1, 0, 0.18),
                    PC4 = centers5[g_i, 4] + rnorm(1, 0, 0.18),
                    stringsAsFactors = FALSE
                  )
                  idx <- idx + 1L
                }}
              }}
              df5 <- dplyr::bind_rows(rows5)
              anc5 <- select_adaptive_reference_anchors(
                df5,
                feature_matrix = as.matrix(df5[, c("PC1", "PC2", "PC3", "PC4")]),
                taxon_col = "species_standardized",
                mahal_q_gold_silver = 0.95
              )
              cda5 <- run_canonical_discriminant_analysis(
                df5,
                feature_cols = c("PC1", "PC2", "PC3", "PC4"),
                target_taxa = taxa5
              )
              stopifnot(identical(as.logical(cda5$is_active), as.logical(anc5$active_mask)))
              stopifnot(all(c("Can1", "Can2", "Can3") %in% colnames(cda5$canonical_scores)))
              stopifnot(cda5$num_canonical_axes == 4L)
              spD_idx <- which(df5$species_standardized == "Packera spD")
              spE_idx <- which(df5$species_standardized == "Packera spE")
              stopifnot(mean(cda5$predicted_taxon[spD_idx] == "Packera spD") >= 0.90)
              stopifnot(mean(cda5$predicted_taxon[spE_idx] == "Packera spE") >= 0.90)
              stopifnot(mean(cda5$posterior_prob[c(spD_idx, spE_idx)]) > 0.85)
              cat("PHASE4_DAWID_SKENE_AND_KFULL_OK\\n")
            }})
            """
            r_batch_file = tmp_p / "batched_step04_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            r_unit_proc = run_warm_r_cmd(["Rscript", str(r_batch_file)], stage_key="step04", cwd=PROJECT_ROOT)
            self.assertEqual(r_unit_proc.returncode, 0, f"Batched Step 04 R helper tests failed:\n{r_unit_proc.stderr}\n{r_unit_proc.stdout}")
            out = r_unit_proc.stdout

            # 1. EPA Level III Ecoregion spatial helper
            epa_line = [ln for ln in out.splitlines() if ln.startswith("EPA_L3:")][0]
            got_epa = epa_line.split("EPA_L3:", 1)[1].strip().split("|")
            self.assertEqual(
                got_epa,
                ["Piedmont", "Middle_Atlantic_Coastal_Plain", "Southeastern_Plains", "South_Central_Plains", "Piedmont", "NA", "NA"],
            )

            # 2. Direct single-species MANOVA factor gating & collinearity safety
            self.assertIn("MANOVA_DIRECT_OK", out)
            grp_out = tmp_p / "manova_group_significance.csv"
            pw_out = tmp_p / "manova_pairwise_posthoc.csv"
            self.assertTrue(grp_out.exists() and pw_out.exists())
            g_df = pd.read_csv(grp_out)
            p_df = pd.read_csv(pw_out)
            self.assertEqual(list(g_df["grouping_factor"]), ["gmm_cluster"])
            self.assertEqual(int(g_df["n_features"].iloc[0]), 2)
            self.assertTrue(pd.isna(g_df["p_value"].iloc[0]) and pd.isna(g_df["q_value"].iloc[0]))
            self.assertEqual(len(p_df), 1)
            self.assertEqual(p_df["contrast"].iloc[0], "Morphotype_K1 vs Morphotype_K2")
            self.assertGreater(float(p_df["hotelling_t2"].iloc[0]), 20.0)
            self.assertTrue(pd.isna(p_df["p_value"].iloc[0]) and pd.isna(p_df["q_value"].iloc[0]))

            # 3. Empirical Bayes & Dawid-Skene + K_full (>3 axes) CDA
            self.assertIn("EB_OK", out)
            self.assertIn("PHASE4_DAWID_SKENE_AND_KFULL_OK", out)

    def test_dissection_collinearity_pruning_lobing_pc_rescue_bimodal_anchors_and_fast_permanova(self):
        """Verify |r|<0.75 scalar pruning, PC3+ lobing rescue, bimodal heteroblastic anchor preservation, and fast PERMANOVA exactness."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        r_code = """
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        options(rgl.useNULL = TRUE)
        source("scripts/morphometrics/04_gmm_morphotools.R")
        set.seed(123)

        # 1. |r| < 0.75 collinearity pruning on active anchors
        n <- 60
        s <- runif(n, 0.35, 0.95)
        df_sc <- data.frame(
          solidity_std = as.numeric(scale(s)),
          q25_solidity_std = as.numeric(scale(s + rnorm(n, 0, 0.02))),
          dissection_index_std = as.numeric(scale(1 / sqrt(s) + rnorm(n, 0, 0.35))),
          within_rosette_efa_dispersion_std = as.numeric(scale(rnorm(n)))
        )
        kept_sc <- prune_collinear_dissection_scalars(
          df_sc,
          base_cols = character(0),
          cand_scalars = c("solidity_std", "q25_solidity_std", "dissection_index_std", "within_rosette_efa_dispersion_std"),
          active_mask = rep(TRUE, n),
          group_vec = rep(c("SpA", "SpB"), each = 30),
          max_r = 0.75
        )
        stopifnot(!all(c("solidity_std", "q25_solidity_std") %in% kept_sc))
        stopifnot("within_rosette_efa_dispersion_std" %in% kept_sc)

        # 2. Subspace rescue of lobing PC4 within floor(N_min / 3)
        n_sub <- 72
        tx <- rep(c("Packera anonyma", "Packera dubia", "Packera obovata"), each = 24)
        lobing_sig <- ifelse(tx == "Packera anonyma", 1.8, -0.9) + rnorm(n_sub, 0, 0.35)
        pc_mat <- cbind(
          PC1 = rnorm(n_sub, 0, 4.0),
          PC2 = rnorm(n_sub, 0, 3.0),
          PC3 = rnorm(n_sub, 0, 0.55),
          PC4 = 0.35 * lobing_sig,
          PC5 = rnorm(n_sub, 0, 0.25)
        )
        aux_df <- data.frame(
          dissection_index = 1.6 + 0.4 * lobing_sig,
          solidity = 0.75 - 0.12 * lobing_sig
        )
        sel_pcs <- select_morphometric_subspace(
          pc_mat,
          group_vec = tx,
          aux_df = aux_df
        )
        stopifnot("PC4" %in% sel_pcs)
        stopifnot(length(sel_pcs) <= floor(24 / 3))

        # 3. Heteroblasty-aware bimodal Gold+Silver anchor preservation
        df_bi <- rbind(
          data.frame(
            catalogNumber = paste0("E", 1:26), plant_individual_id = 0L,
            species_standardized = "Packera paupercula", determiner_tier = "Tier_1_Gold",
            morph_tier_source = "Basal_Rosette", solidity = rnorm(26, 0.90, 0.02),
            q25_solidity = rnorm(26, 0.88, 0.02),
            PC1 = rnorm(26, -1.4, 0.2), PC2 = rnorm(26, 0.3, 0.2), PC3 = rnorm(26, -0.4, 0.15)
          ),
          data.frame(
            catalogNumber = paste0("D", 1:6), plant_individual_id = 0L,
            species_standardized = "Packera paupercula", determiner_tier = "Tier_1_Gold",
            morph_tier_source = "Basal_Rosette", solidity = rnorm(6, 0.46, 0.03),
            q25_solidity = rnorm(6, 0.42, 0.03),
            PC1 = rnorm(6, 2.9, 0.2), PC2 = rnorm(6, -1.5, 0.2), PC3 = rnorm(6, 3.1, 0.2)
          )
        )
        anc_res <- select_adaptive_reference_anchors(
          df_bi,
          feature_matrix = as.matrix(df_bi[, c("PC1", "PC2", "PC3")]),
          taxon_col = "species_standardized"
        )
        stopifnot(sum(anc_res$active_mask[27:32]) >= 5L)

        # 4. Exact equivalence of O(N*P) fast_euclidean_permanova_trace vs vegan::adonis2
        Xp <- matrix(rnorm(90 * 4), nrow = 90, ncol = 4)
        gp <- factor(rep(c("E1", "E2", "E3"), each = 30))
        Xp[31:60, 1] <- Xp[31:60, 1] + 1.4
        fp <- fast_euclidean_permanova_trace(Xp, gp, permutations = 99L, seed = 7L)
        vp <- vegan::adonis2(Xp ~ gp, permutations = 99L, method = "euclidean")
        stopifnot(abs(fp$R2 - vp$R2[1]) < 1e-8)
        stopifnot(abs(fp$F - vp$F[1]) < 1e-8)
        cat("STEP04_DISSECTION_HELPERS_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step04.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Step 04 dissection R helper test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP04_DISSECTION_HELPERS_OK", proc.stdout)

    def test_spatial_gathering_demes_and_per_species_typical_morphotype(self):
        """
        Verify that:
        1. `assign_spatial_gathering_demes()` groups conspecific individuals within
           a 10 km Haversine deme without pre-shrinking feature coordinates.
        2. `compute_type_specimen_morph_control()` evaluates `is_typical_morphotype`
           per-species rather than requiring global `gmm_cluster == 1L`.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        source("scripts/morphometrics/authority_weighting_utils.R")
        source("scripts/morphometrics/plasticity_mitigation_utils.R")
        stopifnot(exists("assign_spatial_gathering_demes", mode = "function"))

        # 1. Test 10 km conspecific deme assignment without pre-modeling coordinate shrinkage
        df_deme <- data.frame(
          catalogNumber = c("SHEET-001", "SHEET-001", "SHEET-002", "SHEET-003", "SHEET-ISOLATE"),
          plant_individual_id = c(0L, 1L, 0L, 0L, 0L),
          species_standardized = rep("Packera dubia", 5L),
          determiner_tier = rep("Tier_1_Gold", 5L),
          authority_weight = rep(0.95, 5L),
          basal_leaf_count = c(2L, 2L, 4L, 4L, 3L),
          decimalLatitude = c(35.910, 35.910, 35.914, 35.918, 37.500),
          decimalLongitude = c(-79.050, -79.050, -79.054, -79.058, -81.200),
          PC1 = c(-2.60, -1.40, -2.05, -1.95, -2.10),
          PC2 = c(0.60, -0.40, 0.12, 0.08, 0.10),
          stringsAsFactors = FALSE
        )
        deme_df <- assign_spatial_gathering_demes(df_deme, deme_radius_km = 10.0)
        stopifnot("population_deme_id" %in% names(deme_df))
        stopifnot("deme_n_specimens" %in% names(deme_df))
        # First 4 rows are within 10 km of each other and share a deme; raw PC1 coordinates stay untouched
        stopifnot(length(unique(deme_df$population_deme_id[1:4])) == 1L)
        stopifnot(deme_df$population_deme_id[5L] != deme_df$population_deme_id[1L])
        stopifnot(all(deme_df$deme_n_specimens[1:4] == 4L))
        stopifnot(max(abs(deme_df$PC1 - df_deme$PC1)) == 0.0)

        # 2. Test per-species evaluation of is_typical_morphotype in compute_type_specimen_morph_control()
        df_typ <- data.frame(
          catalogNumber = c("DUB-TYPE", "DUB-GOLD-1", "DUB-GOLD-2", "OBO-TYPE", "OBO-GOLD-1", "OBO-GOLD-2"),
          plant_individual_id = rep(0L, 6L),
          species_standardized = c(rep("Packera dubia", 3L), rep("Packera obovata", 3L)),
          determiner_tier = rep("Tier_1_Gold", 6L),
          authority_weight = rep(0.98, 6L),
          type_status = c("Holotype", "None", "None", "Holotype", "None", "None"),
          is_type_specimen = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
          gmm_cluster = c(1L, 1L, 1L, 2L, 2L, 2L),
          PC1 = c(-2.20, -2.15, -2.25, 2.30, 2.25, 2.35),
          PC2 = c(0.10, 0.12, 0.08, -0.10, -0.12, -0.08),
          stringsAsFactors = FALSE
        )
        ctrl <- compute_type_specimen_morph_control(
          df = df_typ,
          X_mat = as.matrix(df_typ[, c("PC1", "PC2")]),
          taxon_col = "species_standardized"
        )
        stopifnot(all(ctrl$is_typical_morphotype[1:3]))
        stopifnot(all(ctrl$is_typical_morphotype[4:6]))
        cat("DEME_ASSIGNMENT_AND_PER_SPECIES_TYPICAL_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_deme_and_typical.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Deme assignment / typical morphotype R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("DEME_ASSIGNMENT_AND_PER_SPECIES_TYPICAL_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
