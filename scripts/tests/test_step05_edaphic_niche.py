#!/usr/bin/env python3
"""
===============================================================================
Script: test_step05_edaphic_niche.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Batched warm-R suite for Step 05 fine-grain edaphic pedology, SGMC lithology,
    uncertainty buffering, partitioned D_clim vs D_edaph PCA-env niche tests,
    multi-class substrate contingency, 19-BioClim cohort-adaptive RF gating,
    and North American raster cropping + MaxEnt/Spatial RF SDM maps.
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

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


class TestStep05EdaphicNiche(unittest.TestCase):
    """Batched warm-R tests for Step 05 pedology, niche overlap, 19-BioClim gating, and SDM projection."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._class_temp_dir = tempfile.TemporaryDirectory(prefix="step05_edaph_niche_")
        cls._tmp = Path(cls._class_temp_dir.name)
        cls.con_out = cls._tmp / "contingency.csv"
        cls.r_proc = None

        if rscript_available():
            ssurgo_out = cls._tmp / "ssurgo_val.csv"
            sdm_root = cls._tmp / "sdm_test_root"
            sdm_root.mkdir(parents=True, exist_ok=True)
            r_code = f"""
            if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
            suppressPackageStartupMessages({{ library(readr); library(dplyr); library(terra) }})
            source("scripts/core/db_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/05_multimodal_spatial_rf.R")

            # 1. 0-15 cm trapezoidal rhizosphere integration, log-ratio texture, & strict non-blended coalesce
            local({{
              v_0_5 <- c(4.5, 6.0, NA, 5.2)
              v_5_15 <- c(5.1, 6.6, 5.8, NA)
              rhizo <- compute_rhizosphere_0_15cm(v_0_5, v_5_15)
              stopifnot(abs(rhizo[1] - ((5 * 4.5 + 10 * 5.1) / 15)) < 1e-8)
              stopifnot(abs(rhizo[2] - ((5 * 6.0 + 10 * 6.6) / 15)) < 1e-8)
              stopifnot(abs(rhizo[3] - 5.8) < 1e-8 && abs(rhizo[4] - 5.2) < 1e-8)

              tex <- compute_soil_texture_log_ratios(sand_pct = c(70.0, 30.0, NA), clay_pct = c(10.0, 35.0, 20.0))
              stopifnot(all(c("soil_silt_pct", "ilr_texture_1", "ilr_texture_2") %in% names(tex)))
              stopifnot(!any(c("alr_sand_clay", "alr_silt_clay") %in% names(tex)))
              stopifnot(tex$ilr_texture_1[1] > 0 && tex$ilr_texture_1[2] < 0)
              stopifnot(is.na(tex$ilr_texture_1[3]) && is.na(tex$ilr_texture_2[3]))

              df_in <- data.frame(
                catalogNumber = c("V1", "V2", "V3"), latitude = c(35.8, 36.1, 50.5), longitude = c(-79.1, -86.4, -85.0),
                is_rock_outcrop = c(TRUE, FALSE, FALSE),
                outcrop_class = c("Granitic flatrock", "Non-outcrop matrix", "Non-outcrop matrix"),
                ssurgo_rock_outcrop_pct = c(55.0, 0.0, NA), ssurgo_bedrock_depth_cm = c(18.0, 180.0, NA),
                ssurgo_sand_pct = c(74.0, 32.0, NA), ssurgo_clay_pct = c(8.0, 28.0, NA),
                ssurgo_ph = c(4.7, 6.5, NA), ssurgo_cec = c(4.2, 19.5, NA), ssurgo_aws_0_50cm = c(2.2, 9.0, NA),
                soil_ph = c(6.9, 4.1, 5.8), soil_sand = c(20.0, 85.0, 62.0), soil_cec = c(30.0, 3.0, 14.5),
                stringsAsFactors = FALSE
              )
              df_out <- extract_polaris_and_solus_predictors(df_in, env_dir = tempdir(), demo = TRUE)
              req_cols <- c("fine_bedrock_depth_cm", "fine_cfvo_pct", "fine_soil_cec", "ilr_texture_1", "ilr_texture_2")
              stopifnot(all(req_cols %in% names(df_out)))
              stopifnot(df_out$fine_bedrock_depth_cm[1] < df_out$fine_bedrock_depth_cm[2])
              stopifnot(df_out$fine_cfvo_pct[1] > df_out$fine_cfvo_pct[2])
              stopifnot(abs(df_out$fine_soil_ph[1] - 4.7) < 1e-6 && abs(df_out$fine_soil_ph[2] - 6.5) < 1e-6 && abs(df_out$fine_soil_ph[3] - 5.8) < 1e-6)
              stopifnot(abs(df_out$fine_soil_sand[1] - 74.0) < 1e-6 && abs(df_out$fine_soil_sand[2] - 32.0) < 1e-6 && abs(df_out$fine_soil_sand[3] - 62.0) < 1e-6)
              stopifnot(abs(df_out$fine_soil_cec[1] - 4.2) < 1e-6 && abs(df_out$fine_soil_cec[2] - 19.5) < 1e-6 && abs(df_out$fine_soil_cec[3] - 14.5) < 1e-6)

              sub_pts <- data.frame(gid = 1L, catalogNumber = "V1", latitude = 35.8, longitude = -79.1, stringsAsFactors = FALSE)
              tab_res <- data.frame(
                gid = 1L, mukey = "101", musym = "WkB", muname = "Wake-Rock outcrop complex",
                brockdepmin = 18.0, aws050wta = 2.4, cokey = "101_1", compname = "Wake",
                comppct_r = 65.0, compkind = "Series", majcompflag = "Yes",
                taxsuborder = "Udipsamments", taxgreatgroup = "Quartzipsamments",
                hzdept_r = 0.0, hzdepb_r = 15.0, sandtotal_r = 76.0, claytotal_r = 8.0,
                ph1to1h2o_r = 4.8, om_r = 1.7, cec7_r = 4.1, awc_r = 0.05,
                resdept_r = 18.0, reskind = "Lithic bedrock", stringsAsFactors = FALSE
              )
              agg1 <- aggregate_ssurgo_tabular_to_points(sub_pts, tab_res)
              stopifnot(nrow(agg1) == 1 && isTRUE(agg1$is_rock_outcrop[1]) && agg1$outcrop_class[1] == "Granitic flatrock")
              cat("RHIZOSPHERE_PEDOLOGY_OK\\n")
            }})

            # 2. 6-class USGS SGMC bedrock lithology
            local({{
              stopifnot(length(SGMC_LITHOLOGY_CLASSES) == 6L)
              expected_classes <- c(
                "Ultramafic_Serpentine", "Mafic_Igneous_Metamorphic", "Carbonate_Calcareous",
                "Felsic_Granitic_Gneiss", "Siliciclastic_Consolidated", "Unconsolidated_Sedimentary"
              )
              stopifnot(identical(SGMC_LITHOLOGY_CLASSES, expected_classes))
              mapped <- standardize_sgmc_lithology(
                rocktype1 = c("Serpentinite and peridotite", "Diabase and gabbro", "Limestone and dolostone", "Biotite granite gneiss", "Quartzite and sandstone", "Unconsolidated coastal plain sand"),
                muname = c("Chrome silt loam", "Iredell loam", "Gladeville-Rock outcrop", "Wake-Rock outcrop complex", "Ramsey sandy loam", "Wagram loamy sand")
              )
              stopifnot(identical(as.character(mapped), expected_classes))
              df_geo <- data.frame(
                catalogNumber = paste0("G", 1:6),
                latitude = c(39.7, 35.9, 36.1, 35.8, 35.5, 34.2),
                longitude = c(-76.1, -78.9, -86.4, -80.2, -82.5, -78.1),
                sgmc_rocktype1 = c("Serpentinite", "Diabase dike", "Ordovician limestone", "Granitic gneiss", "Sandstone", "Alluvial sand"),
                is_rock_outcrop = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
                ssurgo_rock_outcrop_pct = c(45, 35, 50, 65, 5, 0),
                fine_bedrock_depth_cm = c(25, 35, 18, 12, 95, 200),
                stringsAsFactors = FALSE
              )
              df_geo <- extract_sgmc_lithology(df_geo, env_dir = tempdir(), demo = TRUE)
              stopifnot(all(c("sgmc_lithology_class", "is_azonal_lithology") %in% names(df_geo)))
              stopifnot(!any(c("azonal_contact_dist_m", "azonal_proximity_index") %in% names(df_geo)))
              stopifnot(all(df_geo$is_azonal_lithology[1:4]) && !df_geo$is_azonal_lithology[6])
              cat("SGMC_LITHOLOGY_OK\\n")
            }})

            # 3. Coordinate-uncertainty buffered extraction & >10000m training exclusion
            local({{
              radii <- compute_uncertainty_buffer_radius(c(NA, 25, 250, 600, 3500, 25000))
              stopifnot(identical(as.numeric(radii), c(250, 100, 250, 600, 1000, 1000)))
              df <- data.frame(
                catalogNumber = c("U1", "U2", "U3", "U4"),
                latitude = c(35.8, 35.9, 36.0, 34.5), longitude = c(-80.1, -80.2, -79.8, -78.5),
                coordinateUncertaintyInMeters = c(NA, 45, 750, 15000),
                fine_soil_ph = c(4.9, 5.1, 6.2, 5.4), fine_bedrock_depth_cm = c(22, 35, 110, 160),
                fine_cfvo_pct = c(32, 24, 12, 6), ssurgo_rock_outcrop_pct = c(55, 30, 5, 0),
                is_rock_outcrop = c(TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE
              )
              df_buf <- extract_uncertainty_buffered_predictors(df, env_dir = tempdir(), demo = TRUE)
              req <- c("buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
                       "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer")
              stopifnot(all(req %in% names(df_buf)))
              stopifnot(identical(df_buf$high_coord_uncertainty_flag, c(FALSE, FALSE, FALSE, TRUE)))
              stopifnot(identical(df_buf$is_training_eligible_coord, c(TRUE, TRUE, TRUE, FALSE)))
              stopifnot(all(df_buf$ph_q10_buffer <= df_buf$fine_soil_ph) && all(df_buf$ph_q90_buffer >= df_buf$fine_soil_ph))
              stopifnot(all(df_buf$min_bedrock_depth_buffer_cm <= df_buf$fine_bedrock_depth_cm))
              stopifnot(all(df_buf$max_cfvo_buffer_pct >= df_buf$fine_cfvo_pct) && df_buf$azonal_lithology_presence_buffer[1] == TRUE)
              cat("UNCERTAINTY_BUFFER_OK\\n")
            }})

            # 4. Block-wise collinearity pruning & partitioned D_clim vs D_edaph PCA-env niche tests
            local({{
              set.seed(42)
              n <- 40; sp <- rep(c("Packera anonyma", "Packera dubia"), each = n / 2); is_dub <- sp == "Packera dubia"
              df <- data.frame(
                catalogNumber = sprintf("V_%03d", seq_len(n)), species_standardized = sp, analysis_group = sp,
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"), length.out = n),
                latitude = ifelse(is_dub, runif(n, 33.5, 36.0), runif(n, 33.8, 36.2)),
                longitude = ifelse(is_dub, runif(n, -81.0, -78.0), runif(n, -81.2, -78.2)),
                coordinateUncertaintyInMeters = c(rep(200, n - 2), 18000, 25000),
                bio1_temp_mean = rnorm(n, 16.5, 0.8), bio4_temp_seasonality = rnorm(n, 620, 25),
                bio12_precip_annual = rnorm(n, 1250, 60), bio15_precip_seasonality = rnorm(n, 18, 2.0),
                delta_doy = ifelse(is_dub, rnorm(n, -10, 4), rnorm(n, 8, 4)),
                fine_soil_ph = ifelse(is_dub, rnorm(n, 4.6, 0.25), rnorm(n, 5.9, 0.30)),
                ilr_texture_1 = ifelse(is_dub, rnorm(n, 1.8, 0.22), rnorm(n, 0.2, 0.22)),
                ilr_texture_2 = ifelse(is_dub, rnorm(n, -0.3, 0.18), rnorm(n, 0.4, 0.18)),
                fine_bedrock_depth_cm = ifelse(is_dub, rnorm(n, 22, 6), rnorm(n, 135, 20)),
                fine_cfvo_pct = ifelse(is_dub, rnorm(n, 38, 6), rnorm(n, 9, 3)),
                topo_hli_30m = ifelse(is_dub, rnorm(n, 0.96, 0.04), rnorm(n, 0.72, 0.04)),
                stringsAsFactors = FALSE
              )
              df$bio1_collinear_dup <- df$bio1_temp_mean * 1.01 + rnorm(n, 0, 0.005)
              env_vars <- c("bio1_temp_mean", "bio1_collinear_dup", "bio4_temp_seasonality", "bio12_precip_annual",
                            "bio15_precip_seasonality", "fine_soil_ph", "ilr_texture_1", "ilr_texture_2",
                            "fine_bedrock_depth_cm", "fine_cfvo_pct", "topo_hli_30m")
              bs <- screen_collinear_predictors_by_block(df[, env_vars], cor.threshold = 0.70, vif.threshold = 5.0)
              stopifnot(length(bs$climatic_variables) >= 1 && length(bs$edaphic_variables) >= 1)
              stopifnot(!(all(c("bio1_temp_mean", "bio1_collinear_dup") %in% bs$climatic_variables)))

              niche_out <- run_warren_niche_identity_tests(df, env_vars, n_perm = 15)
              sum_df <- niche_out$summary
              req_cols <- c("Pair", "Taxon1", "Taxon2", "N1", "N2",
                            "D_clim", "I_clim", "p_equiv_D_clim", "p_bg_D_clim",
                            "D_edaph", "I_edaph", "p_equiv_D_edaph", "p_bg_D_edaph",
                            "D_joint", "I_joint", "Schoener_D", "Warren_D", "Warren_I",
                            "p_value_D", "p_value_I", "p_bg_D_joint", "Niche_Divergence_Axis", "Taxonomic_Inference")
              stopifnot(all(req_cols %in% names(sum_df)))
              exp_active_n <- sum(df$determiner_tier %in% c("Tier_1_Gold", "Tier_2_Silver") & df$coordinateUncertaintyInMeters <= 10000)
              stopifnot((sum_df$N1[1] + sum_df$N2[1]) == exp_active_n && sum_df$D_edaph[1] < sum_df$D_clim[1])
              cat("PARTITIONED_NICHE_OK\\n")
            }})

            # 5. Multi-class lithology & substrate contingency on non-dubia species
            local({{
              df <- data.frame(
                catalogNumber = sprintf("ND_%03d", 1:45),
                latitude = c(rep(36.15, 15), rep(35.80, 15), rep(42.20, 15)),
                longitude = c(rep(-81.65, 15), rep(-79.10, 15), rep(-72.80, 15)),
                species_standardized = c(rep("Packera serpenticola", 15), rep("Packera aurea", 15), rep("Packera schweinitziana", 15)),
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver", "Tier_1_Gold"), 15),
                coordinateUncertaintyInMeters = rep(150.0, 45),
                is_rock_outcrop = c(rep(TRUE, 14), FALSE, rep(FALSE, 15), rep(FALSE, 14), TRUE),
                outcrop_class = c(rep("Serpentine", 14), "Non-outcrop matrix", rep("Non-outcrop matrix", 15), rep("Non-outcrop matrix", 14), "Sandstone/Siliceous"),
                taxsuborder = c(rep("Udepts", 15), rep("Aquults", 14), "Udults", rep("Orthents", 15)),
                sgmc_lithology_class = c(rep("Ultramafic_Serpentine", 14), "Mafic_Igneous_Metamorphic", rep("Unconsolidated_Sedimentary", 15), rep("Siliciclastic_Consolidated", 14), "Carbonate_Calcareous"),
                is_azonal_lithology = c(rep(TRUE, 15), rep(FALSE, 15), rep(FALSE, 14), TRUE),
                ssurgo_rock_outcrop_pct = c(rep(55.0, 14), 5.0, rep(0.0, 15), rep(2.0, 14), 22.0),
                fine_bedrock_depth_cm = c(rep(24.0, 14), 85.0, rep(180.0, 15), rep(135.0, 14), 45.0),
                fine_soil_sand = c(rep(42.0, 15), rep(28.0, 15), rep(51.0, 15)),
                fine_soil_ph = c(rep(6.8, 15), rep(5.4, 15), rep(4.9, 15)),
                ssurgo_aws_0_50cm = c(rep(3.1, 15), rep(9.4, 15), rep(6.8, 15)),
                topo_tpi_30m = c(rep(3.6, 15), rep(-1.8, 15), rep(0.4, 15)),
                topo_twi_30m = c(rep(4.1, 15), rep(9.8, 15), rep(6.2, 15)),
                doy = rep(132, 45), stringsAsFactors = FALSE
              )
              res <- synthesize_multiscale_edaphics(df, ssurgo_out = "{ssurgo_out.as_posix()}", contingency_out = "{cls.con_out.as_posix()}", demo = FALSE)
              con_df <- res$contingency; mc_df <- res$multiclass_contingency
              stopifnot(nrow(con_df) == 4 && con_df$Taxon[1] == "Omnibus_Substrate_Contingency_Across_Active_Groups")
              stopifnot(!any(grepl("dubia|anonyma|paupercula|plattensis", con_df$Taxon, ignore.case = TRUE)))
              req_cols <- c("Taxon", "Outcrop_Count", "Matrix_Count", "Total_Vouchers", "Outcrop_Affinity_Pct",
                            "Odds_Ratio", "CI_Lower", "CI_Upper", "P_Value", "P_Value_BH",
                            "Dominant_Lithology_Class", "Dominant_Substrate_Class", "Top_Enriched_Lithology", "Top_Substrate_Odds_Ratio")
              stopifnot(all(req_cols %in% names(con_df)))
              serp_row <- con_df[con_df$Taxon == "Packera serpenticola", ]
              stopifnot(serp_row$Dominant_Lithology_Class == "Ultramafic_Serpentine" && as.numeric(serp_row$Top_Substrate_Odds_Ratio) > 10.0 && as.numeric(serp_row$P_Value_BH) < 0.01)
              aurea_row <- con_df[con_df$Taxon == "Packera aurea", ]
              stopifnot(aurea_row$Dominant_Lithology_Class == "Unconsolidated_Sedimentary" && aurea_row$Dominant_Substrate_Class %in% c("Unconsolidated_Sedimentary", "Aquults"))
              stopifnot(all(c("SGMC_Lithology", "SSURGO_Outcrop", "SSURGO_Drainage_Suborder") %in% unique(mc_df$Substrate_Domain)) && "P_Value_BH" %in% names(mc_df))

              cons_df <- execute_crossmodal_consensus(df, intraspecific_mode = FALSE)
              stopifnot("substrate_specialist_flag" %in% names(cons_df) && !("piedmont_outcrop_specialist_flag" %in% names(cons_df)))
              stopifnot(sum(cons_df$substrate_specialist_flag[cons_df$species_standardized == "Packera serpenticola"]) >= 14)
              stopifnot(sum(cons_df$substrate_specialist_flag[cons_df$species_standardized == "Packera aurea"]) == 0)
              cat("NON_DUBIA_CONTINGENCY_OK\\n")
            }})

            # 6. Cohort-adaptive 19-BioClim, coordinate-free lithology, Arctic circular phenology, & RF primary gating
            local({{
              stopifnot(length(WC_BIO_LOOKUP) == 19L && identical(unname(WC_BIO_LOOKUP), PACKERA_DB_WC_BIO_COLS))
              sg_nd <- standardize_sgmc_lithology(c("", ""), lat = c(34.7, 35.9), lon = c(-80.5, -79.1), is_demo = FALSE)
              stopifnot(all(sg_nd == "Unknown_Lithology"))
              ss_nd <- classify_ssurgo_outcrop(muname = "", compname = "", lat = 34.7, lon = -80.5, is_demo = FALSE)
              stopifnot(ss_nd$outcrop_class == "Non-outcrop matrix" && !isTRUE(ss_nd$is_rock_outcrop))
              ss_ak <- classify_ssurgo_outcrop(muname = "Rock outcrop", compname = "Wake", soil_sand = 82, soil_ph = 5.1, bedrock_depth_cm = 22, rock_outcrop_pct = 28, lat = 64.5, lon = -148.0, is_demo = FALSE)
              stopifnot(ss_ak$outcrop_class == "Granitic flatrock" && isTRUE(ss_ak$is_rock_outcrop))

              arctic_df <- data.frame(
                catalogNumber = paste0("ARC_", 1:12), species_standardized = rep(c("Packera hyperborea", "Packera contermina"), each = 6),
                determiner_tier = "Tier_1_Gold", latitude = seq(62.0, 67.5, length.out = 12), longitude = seq(-149.0, -138.0, length.out = 12),
                doy = c(228, 232, 235, 238, 241, 245, 236, 240, 244, 248, 252, 255),
                n_anthesis = c(4, 3, 5, 2, 4, 3, 3, 5, 4, 2, 0, 0), total_capitula = rep(5, 12),
                voucher_phenological_state = c(rep("Peak_Anthesis", 10), "Unknown", "Unknown"), stringsAsFactors = FALSE
              )
              ph_res <- compute_phenological_anomalies(arctic_df, demo = FALSE)
              stopifnot(sum(ph_res$df$is_verified_anthesis, na.rm = TRUE) >= 10L && all(is.finite(ph_res$df$delta_doy)))

              set.seed(42); n <- 24; sp <- rep(c("Packera hyperborea", "Packera contermina"), each = 12); base_t <- seq(-4, 4, length.out = n)
              df_tb <- data.frame(
                catalogNumber = paste0("TB_", seq_len(n)), species_standardized = sp, determiner_tier = "Tier_1_Gold",
                latitude = seq(61.0, 66.0, length.out = n), longitude = seq(-150.0, -140.0, length.out = n),
                bio10_mean_temp_warmest_qtr = base_t + ifelse(sp == "Packera hyperborea", -2.5, 2.5),
                bio1_temp_mean = base_t + ifelse(sp == "Packera hyperborea", -0.2, 0.2) + rnorm(n, 0, 0.05),
                fine_soil_cec = ifelse(sp == "Packera hyperborea", 8.0, 24.0) + rnorm(n, 0, 0.8),
                topo_tpi_30m = rnorm(n, 0, 1.0), stringsAsFactors = FALSE
              )
              df_tb$bio1_temp_mean <- 0.92 * df_tb$bio10_mean_temp_warmest_qtr + rnorm(n, 0, 0.9)
              anc_g <- resolve_active_anchor_groups(df_tb)
              cand_tb <- c("bio1_temp_mean", "bio10_mean_temp_warmest_qtr", "fine_soil_cec", "topo_tpi_30m")
              blk_tb <- screen_collinear_predictors_by_block(df_tb[, cand_tb], cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = anc_g)
              stopifnot("bio10_mean_temp_warmest_qtr" %in% blk_tb$climatic_variables)
              rf_tb <- compute_spatial_rf_mems(df_tb, blk_tb$kept_vars, n_trees = 60, vif_res = blk_tb)
              prim_tb <- select_primary_predictors_from_rf(rf_tb, blk_tb, df_tb)
              stopifnot(length(prim_tb$primary_env_vars) >= 1L)
              cons_tb <- execute_crossmodal_consensus(
                df_tb, rf_res = rf_tb, primary_env_vars = prim_tb$primary_env_vars,
                primary_clim_vars = prim_tb$primary_clim_vars, primary_edaph_vars = prim_tb$primary_edaph_vars,
                rf_var_weights = prim_tb$rf_var_weights
              )
              stopifnot(all(c("edaphic_z_score", "climatic_z_score", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors") %in% names(cons_tb)))
              stopifnot(all(is.finite(cons_tb$edaphic_z_score)))
              cat("STEP06_COHORT_ADAPTIVE_OK\\n")
            }})

            # 7. NA raster cropping/DEFLATE compression & default MaxEnt + Spatial RF SDM maps
            local({{
              tmp_root <- {repr(str(sdm_root))}
              env_dir <- file.path(tmp_root, "env"); out_dir <- file.path(tmp_root, "outputs")
              dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
              dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

              r_glob <- terra::rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90, res = 1.0, crs = "EPSG:4326")
              terra::values(r_glob) <- seq_len(terra::ncell(r_glob)) %% 50
              names(r_glob) <- "bio1_temp_mean"
              glob_tif <- file.path(env_dir, "wc2.1_bio_1.tif")
              terra::writeRaster(r_glob, glob_tif, overwrite = TRUE)

              crop_and_compress_na_rasters(env_dir = env_dir)
              ex_after <- terra::ext(terra::rast(glob_tif))
              stopifnot(terra::xmin(ex_after) >= -130.1 && terra::xmax(ex_after) <= -49.9)
              stopifnot(terra::ymin(ex_after) >= 13.9 && terra::ymax(ex_after) <= 66.1)
              stopifnot(length(list.files(env_dir, pattern = "\\\\.crop_tmp\\\\.tif$")) == 0L)

              set.seed(123); n <- 30; sp <- rep(c("Packera dubia", "Packera anonyma"), each = 15)
              df_sdm <- data.frame(
                catalogNumber = paste0("SDM_", seq_len(n)), plant_individual_id = 0L,
                species_standardized = sp, analysis_group = sp,
                determiner_tier = c("Tier_1_Gold", rep("Tier_2_Silver", n - 1)),
                type_status = c("Holotype", rep("None", n - 1)), is_type_specimen = c(TRUE, rep(FALSE, n - 1)),
                Combined_Prior_Weight = c(1.0, rep(0.85, n - 1)), coordinateUncertaintyInMeters = rep(250, n),
                latitude = ifelse(sp == "Packera dubia", seq(32.5, 36.2, length.out = 15), seq(34.0, 38.5, length.out = 15)) + rnorm(n, 0, 0.15),
                longitude = ifelse(sp == "Packera dubia", seq(-81.5, -78.2, length.out = 15), seq(-83.5, -77.5, length.out = 15)) + rnorm(n, 0, 0.15),
                bio1_temp_mean = ifelse(sp == "Packera dubia", 16.5, 13.2) + rnorm(n, 0, 0.8),
                bio12_precip_annual = ifelse(sp == "Packera dubia", 1180, 1320) + rnorm(n, 0, 45),
                fine_soil_sand = ifelse(sp == "Packera dubia", 82.0, 44.0) + rnorm(n, 0, 4.0),
                fine_soil_ph = ifelse(sp == "Packera dubia", 4.9, 6.1) + rnorm(n, 0, 0.25),
                stringsAsFactors = FALSE
              )
              kept_v <- c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_sand", "fine_soil_ph")
              srf_res <- compute_spatial_rf_mems(df_sdm, kept_v, n_trees = 50)
              df_sdm$spatial_block_id <- srf_res$spatial_block_id
              df_sdm$rf_predicted_taxon <- srf_res$rf_predicted_taxon
              df_sdm$rf_predicted_prob <- srf_res$rf_predicted_prob
              df_sdm$rf_given_taxon_prob <- srf_res$rf_given_taxon_prob

              out_empty <- file.path(tmp_root, "outputs_empty")
              suppressWarnings(fit_and_project_maxent_sdm_maps(
                df = df_sdm, srf_res = srf_res, env_dir = file.path(tmp_root, "empty_env"), out_dir = out_empty,
                kept_vars = kept_v, primary_clim_vars = c("bio1_temp_mean", "bio12_precip_annual"),
                primary_edaph_vars = c("fine_soil_sand", "fine_soil_ph"), focal_taxon = "Packera dubia",
                intraspecific_mode = FALSE, demo = FALSE
              ))
              stopifnot(file.exists(file.path(out_empty, "reports", "sdm_model_performance.csv")))
              stopifnot(!file.exists(file.path(out_empty, "rasters", "sdm_packera_dubia_maxent_cloglog.tif")))

              r_tmpl <- terra::rast(xmin = -90, xmax = -72, ymin = 28, ymax = 42, res = 0.25, crs = "EPSG:4326")
              xy_g <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
              r_b1 <- r_tmpl; terra::values(r_b1) <- 22.0 - 0.55 * (xy_g[, 2] - 28.0); names(r_b1) <- "bio1_temp_mean"
              r_b12 <- r_tmpl; terra::values(r_b12) <- 1100.0 + 18.0 * (xy_g[, 1] + 85.0); names(r_b12) <- "bio12_precip_annual"
              r_snd <- r_tmpl; terra::values(r_snd) <- pmin(92.0, pmax(15.0, 65.0 + 3.5 * (xy_g[, 1] + 80.0))); names(r_snd) <- "fine_soil_sand"
              r_ph <- r_tmpl; terra::values(r_ph) <- pmin(7.5, pmax(4.2, 5.4 - 0.12 * (xy_g[, 1] + 80.0))); names(r_ph) <- "fine_soil_ph"
              terra::writeRaster(c(r_b1, r_b12, r_snd, r_ph), file.path(env_dir, "regional_predictors.tif"), overwrite = TRUE)

              fit_and_project_maxent_sdm_maps(
                df = df_sdm, srf_res = srf_res, env_dir = env_dir, out_dir = out_dir,
                kept_vars = kept_v, primary_clim_vars = c("bio1_temp_mean", "bio12_precip_annual"),
                primary_edaph_vars = c("fine_soil_sand", "fine_soil_ph"), focal_taxon = "Packera dubia",
                intraspecific_mode = FALSE, demo = FALSE
              )
              perf_csv <- file.path(out_dir, "reports", "sdm_model_performance.csv")
              fig_png <- file.path(out_dir, "figures", "fig05c_sdm_suitability_and_distribution_map.png")
              tif_maxent <- file.path(out_dir, "rasters", "sdm_packera_dubia_maxent_cloglog.tif")
              tif_ens <- file.path(out_dir, "rasters", "sdm_packera_dubia_ensemble.tif")
              stopifnot(file.exists(perf_csv) && file.info(perf_csv)$size > 0)
              stopifnot(file.exists(fig_png) && file.info(fig_png)$size > 0)
              stopifnot(length(list.files(file.path(out_dir, "figures"), pattern = "\\\\.pdf$")) == 0L)
              stopifnot(file.exists(tif_maxent) && file.info(tif_maxent)$size > 0)
              stopifnot(file.exists(tif_ens) && file.info(tif_ens)$size > 0)
              perf_df <- readr::read_csv(perf_csv, show_col_types = FALSE)
              stopifnot("Packera dubia" %in% perf_df$Taxon_or_Morphotype)
              stopifnot(all(c("AUC_Block_MaxEnt", "CBI_MaxEnt", "AUC_Block_SpatialRF", "CBI_SpatialRF", "AUC_Block_Ensemble", "CBI_Ensemble") %in% names(perf_df)))
              cat("SDM_AND_CROP_OK\\n")
            }})
            """
            r_file = cls._tmp / "batched_step05_edaphic_niche.R"
            r_file.write_text(r_code, encoding="utf-8")
            cls.r_proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step05", cwd=PROJECT_ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def _require_r(self) -> str:
        if self.r_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_proc.returncode,
            0,
            f"Batched Step 05 edaphic/niche R suite failed:\n{self.r_proc.stderr}\n{self.r_proc.stdout}",
        )
        return self.r_proc.stdout

    def test_rhizosphere_0_15cm_sgmc_uncertainty_and_collinearity_niche_tests(self):
        """Verify 0-15 cm trapezoidal rhizosphere pedology, 6-class SGMC lithology, uncertainty buffering, and partitioned climatic vs edaphic niche tests."""
        out = self._require_r()
        for marker in ("RHIZOSPHERE_PEDOLOGY_OK", "SGMC_LITHOLOGY_OK", "UNCERTAINTY_BUFFER_OK", "PARTITIONED_NICHE_OK"):
            self.assertIn(marker, out)

    def test_multiclass_substrate_contingency_non_dubia(self):
        """Verify multi-class SGMC/SSURGO substrate contingency and CSV schema on non-dubia Packera taxa."""
        self.assertIn("NON_DUBIA_CONTINGENCY_OK", self._require_r())
        con_csv_df = pd.read_csv(self.con_out)
        self.assertEqual(con_csv_df.loc[0, "Taxon"], "Omnibus_Substrate_Contingency_Across_Active_Groups")
        self.assertIn("Dominant_Lithology_Class", con_csv_df.columns)
        self.assertIn("P_Value_BH", con_csv_df.columns)

    def test_cohort_adaptive_19_bioclim_rf_primary_gating_arctic_phenology_and_coord_free_lithology(self):
        """Verify all 19 WorldClim bio1..bio19, anchor-ranked 3-block VIF + RF primary gating, Arctic circular phenology, and coordinate-free non-demo lithology."""
        self.assertIn("STEP06_COHORT_ADAPTIVE_OK", self._require_r())

    def test_na_raster_cropping_deflate_and_maxent_spatial_rf_sdm_maps(self):
        """Verify automatic North American raster cropping/DEFLATE compression and MaxEnt + Spatial RF SDM GeoTIFF & fig05c PNG generation."""
        self.assertIn("SDM_AND_CROP_OK", self._require_r())


if __name__ == "__main__":
    unittest.main()
