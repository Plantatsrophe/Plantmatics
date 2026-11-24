#!/usr/bin/env python3
"""
===============================================================================
Script: test_step05_zero_trust_enm.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Zero-trust scientific, anti-leakage, EPV budget, and morphometric
    independence regression suite for Step 05.
===============================================================================
"""

from __future__ import annotations

import re
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


def _extract_r_function_body(src: str, fn_name: str) -> str:
    m = re.search(rf"^{re.escape(fn_name)}\s*<-\s*function\(", src, flags=re.MULTILINE)
    if not m:
        return ""
    rest = src[m.start():]
    m_next = re.search(r"\n[A-Za-z._][A-Za-z0-9._]*\s*<-\s*function\(", rest[len(fn_name):])
    return rest[: len(fn_name) + m_next.start()] if m_next else rest


class TestStep05ZeroTrustENM(unittest.TestCase):
    """Zero-trust ENM, anti-leakage, morphometric independence, and EPV budget suite."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._class_temp_dir = tempfile.TemporaryDirectory(prefix="step05_zt_enm_")
        cls._tmp = Path(cls._class_temp_dir.name)
        cls.r_proc = None

        if rscript_available():
            r_code = f"""
            if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
            suppressPackageStartupMessages({{ library(readr); library(dplyr); library(terra) }})
            source("scripts/core/db_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/morphometrics/pub_figures_utils.R")
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/05_multimodal_spatial_rf.R")
            tmp_root <- {repr(str(cls._tmp))}

            # 1. Non-circular K=1 IBE, accessible background M, & empty-raster hygiene
            local({{
              empty_env <- file.path(tmp_root, "empty_env")
              for (.sd in c("polaris_30m", "solus_100m", "dem_3dep_30m")) dir.create(file.path(empty_env, .sd), recursive = TRUE, showWarnings = FALSE)
              df_real <- data.frame(
                catalogNumber = c("R1", "R2", "R3", "R4"),
                latitude = c(35.6, 35.7, 34.2, 33.9), longitude = c(-79.1, -80.2, -78.4, -79.5),
                coordinateUncertaintyInMeters = c(120, 450, 2500, 14000),
                muname = c("Rion sandy loam, 8 to 15 percent slopes", "Mecklenburg clay loam, 2 to 8 percent slopes",
                           "Wake-Rock outcrop complex, 2 to 15 percent slopes", "Wagram loamy sand, 0 to 6 percent slopes"),
                compname = c("Rion", "Mecklenburg", "Wake", "Wagram"),
                ssurgo_ph = c(5.3, 6.2, 4.8, 5.1), ssurgo_sand_pct = c(66.0, 38.0, 76.0, 84.0),
                ssurgo_clay_pct = c(15.0, 32.0, 8.0, 6.0), ssurgo_cec = c(5.8, 19.5, 4.1, 3.6),
                ssurgo_aws_0_50cm = c(6.4, 7.2, 2.1, 4.0), ssurgo_bedrock_depth_cm = c(160.0, 110.0, 18.0, 200.0),
                ssurgo_rock_outcrop_pct = c(0.0, 0.0, 45.0, 0.0), stringsAsFactors = FALSE
              )
              cl_rion <- classify_ssurgo_outcrop(muname = df_real$muname[1], compname = "Rion", bedrock_depth_cm = 160.0, rock_outcrop_pct = 0.0, is_demo = FALSE)
              cl_meck <- classify_ssurgo_outcrop(muname = df_real$muname[2], compname = "Mecklenburg", bedrock_depth_cm = 110.0, rock_outcrop_pct = 0.0, is_demo = FALSE)
              cl_wake <- classify_ssurgo_outcrop(muname = df_real$muname[3], compname = "Wake", bedrock_depth_cm = 18.0, rock_outcrop_pct = 45.0, is_demo = FALSE)
              stopifnot(!isTRUE(cl_rion$is_rock_outcrop) && cl_rion$outcrop_class == "Non-outcrop matrix")
              stopifnot(!isTRUE(cl_meck$is_rock_outcrop) && cl_meck$outcrop_class == "Non-outcrop matrix")
              stopifnot(isTRUE(cl_wake$is_rock_outcrop) && cl_wake$outcrop_class == "Granitic flatrock")

              df_real$is_rock_outcrop <- c(FALSE, FALSE, TRUE, FALSE)
              df_real$outcrop_class <- c("Non-outcrop matrix", "Non-outcrop matrix", "Granitic flatrock", "Non-outcrop matrix")
              df_ps <- extract_polaris_and_solus_predictors(df_real, env_dir = empty_env, demo = FALSE)
              df_ps <- extract_micro_topography_3dep(df_ps, env_dir = empty_env, demo = FALSE)
              df_ps <- extract_sgmc_lithology(df_ps, env_dir = empty_env, demo = FALSE)
              df_ps <- extract_uncertainty_buffered_predictors(df_ps, env_dir = empty_env, demo = FALSE)
              stopifnot(all(is.na(df_ps$polaris_ph_30m)) && all(is.na(df_ps$solus_fragvol_pct)))
              stopifnot(all(is.na(df_ps$topo_tpi_30m)) && all(is.na(df_ps$topo_twi_30m)) && all(is.finite(df_ps$topo_slope_deg)))
              stopifnot(all(is.na(df_ps$ph_q10_buffer)) && all(is.na(df_ps$azonal_lithology_presence_buffer)))
              stopifnot(identical(df_ps$is_edaphic_eligible_coord, c(TRUE, TRUE, FALSE, FALSE)))
              stopifnot(identical(df_ps$is_clim_eligible_coord, c(TRUE, TRUE, TRUE, FALSE)))
              stopifnot(!df_ps$is_azonal_lithology[1] && !df_ps$is_azonal_lithology[2] && df_ps$is_azonal_lithology[3])

              set.seed(2026); n_dub <- 60
              base_lat <- rep(seq(32.0, 37.5, length.out = 30), each = 2) + rnorm(n_dub, 0, 0.03)
              base_lon <- rep(seq(-84.0, -77.0, length.out = 30), each = 2) + rnorm(n_dub, 0, 0.03)
              is_fl <- rep(c(FALSE, TRUE), 30)
              df_sym <- data.frame(
                catalogNumber = sprintf("DUB_%03d", seq_len(n_dub)), plant_individual_id = 1L,
                species_standardized = rep("Packera dubia", n_dub),
                determiner_tier = c(rep(c("Tier_1_Gold", "Tier_2_Silver"), 28), rep("Tier_3_Bronze", 3), "Tier_1_Gold"),
                coordinateUncertaintyInMeters = c(rep(180.0, n_dub - 1), 16000.0),
                latitude = base_lat, longitude = base_lon, gmm_cluster = rep(1L, n_dub),
                PC1 = rnorm(n_dub, 0, 1.0), PC2 = rnorm(n_dub, 0, 0.8), aspect_ratio = rnorm(n_dub, 2.2, 0.25),
                doy = round(95 + 4.0 * (base_lat - 32.0) + rnorm(n_dub, 0, 6)), delta_doy = rnorm(n_dub, 0, 5.5),
                fine_soil_ph = rnorm(n_dub, 5.3, 0.35),
                fine_soil_sand = ifelse(is_fl, rnorm(n_dub, 74, 5), rnorm(n_dub, 78, 6)),
                fine_soil_clay = rnorm(n_dub, 9.0, 2.0), fine_soil_cec = rnorm(n_dub, 4.8, 1.1),
                fine_bedrock_depth_cm = ifelse(is_fl, rnorm(n_dub, 28, 8), rnorm(n_dub, 165, 20)),
                ssurgo_rock_outcrop_pct = ifelse(is_fl, runif(n_dub, 20, 55), rep(0.0, n_dub)),
                ilr_texture_1 = rnorm(n_dub, 1.6, 0.25), ilr_texture_2 = rnorm(n_dub, 0.1, 0.20),
                is_azonal_lithology = is_fl,
                bio10_mean_temp_warmest_qtr = 27.2 - 0.55 * (base_lat - 32.0) + rnorm(n_dub, 0, 0.35),
                bio12_precip_annual = rnorm(n_dub, 1220, 75), bio15_precip_seasonality = rnorm(n_dub, 16, 2.2),
                stringsAsFactors = FALSE
              )
              cand_v <- discover_candidate_eco_predictors(df_sym, include_pheno = FALSE)
              stopifnot(!any(c("is_azonal_lithology", "ssurgo_rock_outcrop_pct") %in% cand_v))
              df_res <- resolve_analysis_group(df_sym, intraspecific_mode = TRUE)
              stopifnot(identical(unique(df_res$analysis_group), "Packera dubia") && identical(unique(resolve_plot_group_column(df_res)), "Packera dubia"))
              srf_k1 <- compute_spatial_rf_mems(df_res, cand_v, n_trees = 60, intraspecific_mode = TRUE)
              stopifnot(isTRUE(srf_k1$single_lineage_mode) && all(srf_k1$rf_predicted_taxon == "Packera dubia"))
              stopifnot(srf_k1$rf_prediction_source[n_dub] == "Passive_Projection_High_Uncertainty" && srf_k1$rf_prediction_source[n_dub - 1] == "Passive_Projection_Bronze")
              stopifnot(!is.null(srf_k1$ibe_summary) && nrow(srf_k1$ibe_summary) >= 1L && all(srf_k1$ibe_summary$Partial_RDA_p > 0.01))
              niche_k1 <- run_warren_niche_identity_tests(df_res, cand_v, n_perm = 15, intraspecific_mode = TRUE)
              stopifnot(nrow(niche_k1$summary) == 1L && niche_k1$summary$Schoener_D[1] == 1.0)
              fig_k1_png <- file.path(tmp_root, "k1_fig", "fig05a_ecological_niche_and_rf_importance.png")
              export_spatial_rf_figures(df_res, srf_k1, niche_k1, fig_k1_png)
              stopifnot(file.exists(fig_k1_png) && length(list.files(file.path(tmp_root, "k1_fig"), pattern = "\\\\.pdf$")) == 0L)
              cat("NONCIRCULAR_K1_IBE_AND_HYGIENE_OK\\n")
            }})

            # 2. Zero-trust scientific, Hopkins' Law, anti-leakage, invariance, & 30 m DEM HLI/TWI suite
            local({{
              df_w <- data.frame(
                catalogNumber = paste0("W", 1:12),
                fine_soil_ph = seq(4.5, 6.7, length.out = 12), ilr_texture_1 = seq(-1.2, 1.5, length.out = 12),
                ilr_texture_2 = seq(-0.6, 0.8, length.out = 12), fine_soil_sand = seq(25, 80, length.out = 12),
                fine_soil_clay = seq(35, 8, length.out = 12), soil_silt_pct = seq(40, 12, length.out = 12),
                fine_soil_cec = seq(4, 22, length.out = 12), fine_cfvo_pct = seq(2, 42, length.out = 12),
                fine_bedrock_depth_cm = seq(20, 180, length.out = 12), polaris_ph_30m = seq(4.5, 6.7, length.out = 12),
                ssurgo_ph = seq(4.6, 6.8, length.out = 12), ssurgo_rock_outcrop_pct = seq(0, 60, length.out = 12),
                ssurgo_aws_0_50cm = seq(2.0, 10.0, length.out = 12), ph_q10_buffer = seq(4.1, 6.2, length.out = 12),
                min_bedrock_depth_buffer_cm = seq(10, 140, length.out = 12), elevation_m = seq(100, 1200, length.out = 12),
                topo_slope_deg = seq(1.5, 24.0, length.out = 12), topo_tpi_30m = seq(-3.0, 4.0, length.out = 12),
                topo_hli_30m = seq(0.62, 1.02, length.out = 12), topo_twi_30m = seq(4.5, 11.5, length.out = 12),
                bio1_temp_mean = seq(7.0, 18.0, length.out = 12), bio3_isothermality = seq(28, 39, length.out = 12),
                bio7_temp_annual_range = seq(26, 37, length.out = 12), bio8_mean_temp_wettest_qtr = seq(15, 24, length.out = 12),
                bio9_mean_temp_driest_qtr = seq(2, 14, length.out = 12), bio12_precip_annual = seq(850, 1450, length.out = 12),
                bio18_precip_warmest_qtr = seq(240, 420, length.out = 12), bio19_precip_coldest_qtr = seq(180, 340, length.out = 12),
                delta_doy = seq(-15, 15, length.out = 12), doy = seq(95, 150, length.out = 12), stringsAsFactors = FALSE
              )
              cand_w <- discover_candidate_eco_predictors(df_w)
              forbidden_w <- c("bio3_isothermality", "bio7_temp_annual_range", "bio8_mean_temp_wettest_qtr", "bio9_mean_temp_driest_qtr",
                               "bio18_precip_warmest_qtr", "bio19_precip_coldest_qtr", "delta_doy", "doy", "elevation_m",
                               "ph_q10_buffer", "min_bedrock_depth_buffer_cm", "polaris_ph_30m", "ssurgo_ph",
                               "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm", "soil_silt_pct", "fine_soil_sand", "fine_soil_clay")
              expected_w <- c("fine_soil_ph", "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm",
                              "topo_slope_deg", "topo_tpi_30m", "topo_hli_30m", "topo_twi_30m", "bio1_temp_mean", "bio12_precip_annual")
              stopifnot(!any(forbidden_w %in% cand_w) && all(expected_w %in% cand_w))

              set.seed(77); n_h <- 24
              lat_h <- seq(34.0, 38.0, length.out = n_h); lon_h <- rep(c(-83.0, -78.0), length.out = n_h)
              elev_h <- c(rep(200.0, n_h - 2), 1200.0, 1400.0); yr_h <- rep(seq(1960, 2020, length.out = 6), 4)
              true_doy <- round(115.0 + 4.0 * (lat_h - 34.0) + 0.80 * (lon_h - median(lon_h)) + (4.0 / 121.92) * (elev_h - median(elev_h)) - 0.10 * (yr_h - median(yr_h)))
              df_hop <- data.frame(
                catalogNumber = paste0("HOP_", seq_len(n_h)), species_standardized = "Packera dubia", determiner_tier = "Tier_1_Gold",
                decimalLatitude = lat_h, decimalLongitude = lon_h, elevation_m = elev_h, year = yr_h, doy = true_doy,
                reproductiveCondition = "fl", voucher_phenological_state = "anthesis", total_capitula = 5, n_anthesis = 4, stringsAsFactors = FALSE
              )
              hop_res <- compute_phenological_anomalies(df_hop, demo = FALSE)
              stopifnot(abs(hop_res$baseline$b_elev - (4.0 / 121.92)) < 1e-6 && abs(hop_res$baseline$b_lon - 0.80) < 1e-6)
              stopifnot(hop_res$baseline$b_year < -0.04 && all(abs(hop_res$df$delta_doy[23:24]) < 3.0))

              set.seed(101); n_lk <- 80; sp_lk <- rep(c("Packera dubia", "Packera paupercula"), each = 40)
              lat_lk <- rep(seq(36.0, 47.0, length.out = 40), 2) + rnorm(n_lk, 0, 0.05)
              lon_lk <- rep(seq(-84.0, -76.0, length.out = 40), 2) + rnorm(n_lk, 0, 0.05)
              bed_sparse <- rep(NA_real_, n_lk); bed_sparse[1:3] <- c(30.0, 32.0, 34.0); bed_sparse[41:43] <- c(150.0, 155.0, 160.0)
              df_lk_sparse <- data.frame(
                catalogNumber = paste0("LK_", seq_len(n_lk)), species_standardized = sp_lk, analysis_group = sp_lk,
                determiner_tier = "Tier_1_Gold", coordinateUncertaintyInMeters = 200.0, latitude = lat_lk, longitude = lon_lk,
                bio1_temp_mean = rnorm(n_lk, 11.0, 1.0), bio12_precip_annual = rnorm(n_lk, 1050.0, 40.0),
                fine_bedrock_depth_cm = bed_sparse, stringsAsFactors = FALSE
              )
              rf_noleak <- compute_spatial_rf_mems(df_lk_sparse, c("bio1_temp_mean", "bio12_precip_annual", "fine_bedrock_depth_cm"), n_trees = 100)
              stopifnot(!("fine_bedrock_depth_cm" %in% rf_noleak$importance$Predictor) && rf_noleak$Spatial_Block_CV_Accuracy < 0.72)
              bed_55 <- rep(NA_real_, n_lk)
              bed_55[1:22] <- c(40.0 + seq(-1.1, 1.1, length.out = 12), 160.0 + seq(-0.9, 0.9, length.out = 10))
              bed_55[41:62] <- c(40.0 + seq(-0.9, 0.9, length.out = 10), 160.0 + seq(-1.1, 1.1, length.out = 12))
              df_lk_55 <- df_lk_sparse; df_lk_55$fine_bedrock_depth_cm <- bed_55
              rf_noleak_55 <- compute_spatial_rf_mems(df_lk_55, c("bio1_temp_mean", "bio12_precip_annual", "fine_bedrock_depth_cm"), n_trees = 100)
              stopifnot("fine_bedrock_depth_cm" %in% rf_noleak_55$importance$Predictor && rf_noleak_55$Spatial_Block_CV_Accuracy < 0.72)

              df_clean <- data.frame(
                catalogNumber = paste0("INV_", 1:24),
                species_standardized = rep(c("Packera dubia", "Packera anonyma"), each = 12),
                analysis_group = rep(c("Packera dubia", "Packera anonyma"), each = 12),
                determiner_tier = rep(c(rep("Tier_1_Gold", 5), rep("Tier_2_Silver", 4), rep("Tier_3_Bronze", 2), "Tier_1_Gold"), 2),
                is_confident_label_noise = rep(c(rep(FALSE, 11), TRUE), 2), coordinateUncertaintyInMeters = 180.0,
                latitude = rep(seq(34.0, 36.5, length.out = 12), 2), longitude = rep(seq(-81.0, -78.5, length.out = 12), 2),
                fine_soil_ph = c(rnorm(12, 4.9, 0.2), rnorm(12, 6.3, 0.2)), ilr_texture_1 = c(rnorm(12, 1.6, 0.2), rnorm(12, -0.3, 0.2)),
                bio1_temp_mean = c(rnorm(12, 16.0, 0.4), rnorm(12, 13.5, 0.4)),
                is_verified_anthesis = c(rep(TRUE, 10), FALSE, FALSE, rep(TRUE, 10), FALSE, FALSE),
                doy = c(rep(108, 10), NA, 260, rep(132, 10), NA, 25),
                cda_predicted_taxon = rep(c("Packera dubia", "Packera anonyma"), each = 12), cda_posterior_prob = 0.94, stringsAsFactors = FALSE
              )
              df_corrupt <- df_clean; bad_rows <- which(df_corrupt$determiner_tier == "Tier_3_Bronze" | df_corrupt$is_confident_label_noise)
              df_corrupt$fine_soil_ph[bad_rows] <- c(9.5, 2.5, 9.8, 2.2, 9.2, 2.8)
              df_corrupt$ilr_texture_1[bad_rows] <- c(-5.0, 5.0, -4.8, 4.8, -5.2, 5.2)
              df_corrupt$bio1_temp_mean[bad_rows] <- c(-10.0, 35.0, -8.0, 38.0, -12.0, 40.0)
              rf_c1 <- compute_spatial_rf_mems(df_clean, c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), n_trees = 80)
              rf_c2 <- compute_spatial_rf_mems(df_corrupt, c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), n_trees = 80)
              act_rows <- setdiff(seq_len(24), bad_rows)
              stopifnot(max(abs(rf_c1$oof_prob_matrix[act_rows, ] - rf_c2$oof_prob_matrix[act_rows, ])) < 1e-9)
              cons_c1 <- execute_crossmodal_consensus(df_clean, rf_res = rf_c1, primary_env_vars = c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), primary_edaph_vars = c("fine_soil_ph", "ilr_texture_1"), primary_clim_vars = "bio1_temp_mean")
              cons_c2 <- execute_crossmodal_consensus(df_corrupt, rf_res = rf_c2, primary_env_vars = c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), primary_edaph_vars = c("fine_soil_ph", "ilr_texture_1"), primary_clim_vars = "bio1_temp_mean")
              stopifnot(max(abs(cons_c1$edaphic_z_score[act_rows] - cons_c2$edaphic_z_score[act_rows])) < 1e-9)
              stopifnot(!any(cons_c1$triage_category[c(11, 12, 23, 24)] == "Phenological_Anomaly_Outlier") && all(cons_c1$n_observed_modalities[c(11, 12, 23, 24)] == 2L))

              dem_dir <- file.path(tmp_root, "dem_env", "dem_3dep_30m"); dir.create(dem_dir, recursive = TRUE, showWarnings = FALSE)
              r_dem <- terra::rast(xmin = -79.20, xmax = -78.80, ymin = 35.40, ymax = 35.80, nrows = 160, ncols = 160, crs = "EPSG:4326")
              xy_d <- terra::xyFromCell(r_dem, seq_len(terra::ncell(r_dem)))
              terra::values(r_dem) <- 200.0 + 1800.0 * abs(xy_d[, 1] - (-79.00)) + 400.0 * (xy_d[, 2] - 35.40)
              terra::writeRaster(r_dem, file.path(dem_dir, "N35W080.tif"), overwrite = TRUE)
              pts_dem <- data.frame(catalogNumber = c("VALLEY_FLOOR", "UPPER_FLANK_WEST_FACING", "UPPER_FLANK_EAST_FACING"),
                                    latitude = c(35.55, 35.55, 35.55), longitude = c(-79.00, -78.95, -79.05), stringsAsFactors = FALSE)
              dem_out <- extract_micro_topography_3dep(pts_dem, env_dir = file.path(tmp_root, "dem_env"), demo = FALSE)
              stopifnot(all(is.finite(dem_out$topo_hli_30m)) && all(is.finite(dem_out$topo_twi_30m)) && all(is.finite(dem_out$elevation_m)))
              stopifnot(dem_out$topo_twi_30m[1] > dem_out$topo_twi_30m[2] && dem_out$topo_tpi_30m[1] < 0.0 && dem_out$topo_hli_30m[2] > dem_out$topo_hli_30m[3])

              imp_ilr <- rf_c1$importance
              stopifnot(identical(as.character(imp_ilr$Group[imp_ilr$Predictor == "ilr_texture_1"]), "[A] Fine-Grain Edaphic & Outcrop Lithology"))
              imp_ilr_ed <- imp_ilr[!grepl("^MEM[0-9]+$", imp_ilr$Predictor) & imp_ilr$Group == "[A] Fine-Grain Edaphic & Outcrop Lithology", , drop = FALSE]
              imp_ilr_ed <- imp_ilr_ed[order(imp_ilr_ed$MeanDecreaseAccuracy, decreasing = TRUE), , drop = FALSE]
              stopifnot(nrow(imp_ilr_ed) >= 1L && imp_ilr_ed$Predictor[1L] %in% c("fine_soil_ph", "ilr_texture_1"))
              t0_mx <- proc.time()[["elapsed"]]
              mx_sep <- .fit_single_maxent_model(
                p_vec = c(rep(1L, 15), rep(0L, 15)),
                X_df = data.frame(bio1_temp_mean = c(rnorm(15, 16.5, 0.3), rnorm(15, 11.0, 0.3)), fine_soil_ph = c(rnorm(15, 4.8, 0.15), rnorm(15, 6.6, 0.15))),
                w_vec = c(1.0, rep(0.85, 14), rep(1.0, 15))
              )
              stopifnot((proc.time()[["elapsed"]] - t0_mx) < 10.0 && all(is.finite(mx_sep$predict_cloglog(data.frame(bio1_temp_mean = c(16.5, 11.0), fine_soil_ph = c(4.8, 6.6))))))
              if (requireNamespace("glmnet", quietly = TRUE)) stopifnot(glmnet::glmnet.control()$fdev > 0)
              cat("STEP05_ZERO_TRUST_SUITE_OK\\n")
            }})

            # 3. Composite-only geochem/wetland ENM inputs, cross-block VIF <= 5.0, & EPV predictor budgets
            local({{
              set.seed(202611); n_g1 <- 24L; sp_g1 <- rep(c("Packera serpenticola", "Packera aurea"), each = 12L)
              df_g1 <- data.frame(
                catalogNumber = sprintf("G1_%02d", seq_len(n_g1)), species_standardized = sp_g1, determiner_tier = "Tier_1_Gold",
                latitude = c(seq(35.0, 35.4, length.out = 12L), seq(36.1, 36.6, length.out = 12L)),
                longitude = c(seq(-83.8, -83.4, length.out = 12L), seq(-81.6, -81.1, length.out = 12L)),
                coordinateUncertaintyInMeters = 150.0, fine_soil_ph = c(rnorm(12L, 6.4, 0.2), rnorm(12L, 5.5, 0.2)),
                ilr_texture_1 = rnorm(n_g1, 0.2, 0.3), ilr_texture_2 = rnorm(n_g1, -0.1, 0.25),
                topo_twi_30m = c(rnorm(12L, 4.2, 0.4), rnorm(12L, 11.4, 0.5)),
                bio1_temp_mean = c(rnorm(12L, 13.8, 0.4), rnorm(12L, 12.6, 0.4)), bio12_precip_annual = rnorm(n_g1, 1380.0, 50.0),
                soil_ca_mg_molar_ratio = c(rnorm(12L, 0.12, 0.03), rnorm(12L, 3.8, 0.4)),
                soil_ultramafic_ni_cr_index = c(rnorm(12L, 2.6, 0.3), rnorm(12L, 0.15, 0.05)),
                soil_carbonate_c_pct = c(rnorm(12L, 0.05, 0.01), rnorm(12L, 2.2, 0.3)),
                is_nwi_wetland = c(rep(0L, 12L), rep(1L, 12L)), nwi_wetland_system = c(rep("UPL", 12L), rep("PEM1A", 12L)),
                nwi_wetland_dist_m = c(runif(12L, 220, 480), runif(12L, 0, 12)),
                nwpl_indicator_code = c(rep("FACU", 12L), rep("FACW", 12L)), nwpl_wetland_score = c(rep(4.0, 12L), rep(2.0, 12L)),
                stringsAsFactors = FALSE
              )
              df_g1 <- compute_derived_habitat_predictors(df_g1)
              cand_g1 <- discover_candidate_eco_predictors(df_g1, include_pheno = FALSE)
              raw_ban <- c("soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct",
                           "ultramafic_calcareous_gradient", "is_nwi_wetland", "nwi_wetland_system",
                           "nwi_wetland_dist_m", "nwpl_indicator_code", "nwpl_wetland_score")
              stopifnot("wetland_seepage_index" %in% cand_g1 && !any(raw_ban %in% cand_g1))
              ct_g1 <- evaluate_multiclass_substrate_contingency(df_g1)
              stopifnot(is.list(ct_g1) && nrow(ct_g1$multiclass_details) >= 1L)

              n_g2 <- 30L; sp_g2 <- rep(c("Packera dubia", "Packera anonyma"), each = 15L); s_g2 <- ifelse(sp_g2 == "Packera dubia", -1.0, 1.0)
              u1_g2 <- stats::resid(stats::lm(rnorm(n_g2) ~ s_g2)); u2_g2 <- stats::resid(stats::lm(rnorm(n_g2) ~ u1_g2 + s_g2))
              u1_g2 <- u1_g2 / stats::sd(u1_g2); u2_g2 <- u2_g2 / stats::sd(u2_g2)
              wsi_base <- u1_g2 + 0.6 * s_g2; bio1_base <- 14.0 + u2_g2 - 0.5 * s_g2
              df_g2 <- data.frame(
                fine_soil_ph = c(rnorm(15L, 4.9, 0.3), rnorm(15L, 6.2, 0.3)), wetland_seepage_index = wsi_base,
                topo_twi_30m = round(6.5 + 3.2 * wsi_base, 5), bio1_temp_mean = bio1_base,
                elevation_m = round(2500.0 - 130.0 * bio1_base, 5), bio12_precip_annual = rnorm(n_g2, 1200.0, 90.0), stringsAsFactors = FALSE
              )
              bs_g2 <- screen_collinear_predictors_by_block(df_g2, cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = sp_g2)
              stopifnot("wetland_seepage_index" %in% bs_g2$kept_vars && !("topo_twi_30m" %in% bs_g2$kept_vars))
              stopifnot("bio1_temp_mean" %in% bs_g2$kept_vars && !("elevation_m" %in% bs_g2$kept_vars))
              R_g2 <- stats::cor(df_g2[, bs_g2$kept_vars, drop = FALSE])
              stopifnot(max(diag(solve(R_g2 + diag(1e-6, length(bs_g2$kept_vars)))), na.rm = TRUE) <= 5.0 + 1e-4)

              make_orth_cohort <- function(n_anc, n_bronze = 0L) {{
                n_tot <- n_anc + n_bronze; n_half <- floor(n_anc / 2L)
                sp_anc <- c(rep("Packera serpenticola", n_half), rep("Packera anonyma", n_anc - n_half))
                sp_v <- c(sp_anc, rep("Packera anonyma", n_bronze)); tier_v <- c(rep("Tier_1_Gold", n_anc), rep("Tier_3_Bronze", n_bronze))
                s_anc <- ifelse(sp_anc == "Packera serpenticola", 1.0, -1.0)
                Q_mat <- qr.Q(qr(stats::resid(stats::lm(matrix(rnorm(n_anc * 16L), nrow = n_anc, ncol = 16L) ~ s_anc)))) * sqrt(n_anc - 1L)
                if (n_bronze > 0L) Q_mat <- rbind(Q_mat, matrix(rnorm(n_bronze * 16L), nrow = n_bronze, ncol = 16L))
                pred_nms <- c("fine_soil_ph", "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm",
                              "lithic_outcrop_index", "wetland_seepage_index", "topo_slope_deg", "topo_tpi_30m", "topo_hli_30m",
                              "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality", "summer_moisture_deficit_idx")
                out_df <- as.data.frame(Q_mat, stringsAsFactors = FALSE); names(out_df) <- pred_nms
                s_full <- c(s_anc, rep(-1.0, n_bronze))
                out_df$fine_soil_ph <- out_df$fine_soil_ph + 1.2 * s_full; out_df$bio1_temp_mean <- out_df$bio1_temp_mean - 1.2 * s_full
                out_df$catalogNumber <- sprintf("EPV_%d_%03d", n_anc, seq_len(n_tot))
                out_df$species_standardized <- sp_v; out_df$analysis_group <- sp_v; out_df$determiner_tier <- tier_v
                out_df$coordinateUncertaintyInMeters <- 200.0
                out_df$latitude <- seq(34.5, 37.8, length.out = n_tot) + rnorm(n_tot, 0, 0.05)
                out_df$longitude <- seq(-83.5, -79.0, length.out = n_tot) + rnorm(n_tot, 0, 0.05)
                out_df
              }}
              df_narrow <- make_orth_cohort(18L, n_bronze = 4L)
              cand_n <- discover_candidate_eco_predictors(df_narrow, include_pheno = FALSE); stopifnot(length(cand_n) == 16L)
              srf_n <- compute_spatial_rf_mems(df_narrow, cand_n, n_trees = 50)
              stopifnot(ncol(srf_n$mems) == 2L && identical(colnames(srf_n$mems), c("MEM1", "MEM2")))
              stopifnot(length(srf_n$collinearity_screening$kept_vars) == 4L && nrow(srf_n$importance) == 6L)
              df_med <- make_orth_cohort(45L)
              srf_m <- compute_spatial_rf_mems(df_med, discover_candidate_eco_predictors(df_med, include_pheno = FALSE), n_trees = 50)
              stopifnot(length(srf_m$collinearity_screening$kept_vars) >= 4L && length(srf_m$collinearity_screening$kept_vars) <= 9L && ncol(srf_m$mems) <= 5L)
              df_wide <- make_orth_cohort(120L)
              srf_w <- compute_spatial_rf_mems(df_wide, discover_candidate_eco_predictors(df_wide, include_pheno = FALSE), n_trees = 50)
              stopifnot(ncol(srf_w$mems) == 5L && length(srf_w$collinearity_screening$kept_vars) == 12L)
              cat("ENM_COMPOSITE_VIF_EPV_GUARDS_OK\\n")
            }})

            # 4. Step 05 morphometric independence & min() bedrock buffer (source Step 07 last inside local)
            local({{
              source("scripts/analysis/07_triage_dashboard_synthesis.R")
              sol_dir <- file.path(tmp_root, "env_indep", "solus_100m"); dir.create(sol_dir, recursive = TRUE, showWarnings = FALSE)
              r_bd <- terra::rast(xmin = -80.10, xmax = -79.90, ymin = 35.90, ymax = 36.10, nrows = 100, ncols = 100, crs = "EPSG:4326")
              terra::values(r_bd) <- 180.0; r_bd[terra::cellFromXY(r_bd, cbind(-80.001, 36.000))] <- 18.0
              terra::writeRaster(r_bd, file.path(sol_dir, "anylithicdpt_cm_p.tif"), overwrite = TRUE)
              df_out_buf <- extract_uncertainty_buffered_predictors(
                data.frame(catalogNumber = "FLATROCK_01", latitude = 36.000, longitude = -80.000, coordinateUncertaintyInMeters = 350.0, stringsAsFactors = FALSE),
                env_dir = file.path(tmp_root, "env_indep"), demo = FALSE
              )
              stopifnot(is.finite(df_out_buf$min_bedrock_depth_buffer_cm[1]) && df_out_buf$min_bedrock_depth_buffer_cm[1] <= 20.0)

              set.seed(2027); n_h <- 16; hab <- rep(c("Flatrock", "Sandhills", "CoastalDune", "WesternPrairie"), each = n_h); n_tot <- length(hab)
              lat_v <- ifelse(hab == "Flatrock", 34.65 + rnorm(n_tot, 0, 0.04), ifelse(hab == "Sandhills", 34.60 + rnorm(n_tot, 0, 0.04), ifelse(hab == "CoastalDune", 33.85 + rnorm(n_tot, 0, 0.08), 34.10 + rnorm(n_tot, 0, 0.25))))
              lon_v <- ifelse(hab == "Flatrock", -80.55 + rnorm(n_tot, 0, 0.04), ifelse(hab == "Sandhills", -80.48 + rnorm(n_tot, 0, 0.04), ifelse(hab == "CoastalDune", -78.00 + rnorm(n_tot, 0, 0.08), -94.20 + rnorm(n_tot, 0, 0.30))))
              df_dub <- data.frame(
                catalogNumber = sprintf("PD_%03d", seq_len(n_tot)), plant_individual_id = 1L,
                species_standardized = "Packera dubia", determiner_tier = "Tier_1_Gold", coordinateUncertaintyInMeters = 200.0,
                latitude = lat_v, longitude = lon_v, decimalLatitude = lat_v, decimalLongitude = lon_v,
                gmm_cluster = 1L, cda_predicted_taxon = "Morphotype_K1", morphometric_outlier_flag = FALSE,
                PC1 = rnorm(n_tot, 0, 0.6), PC2 = rnorm(n_tot, 0, 0.5),
                doy = round(102 + rnorm(n_tot, 0, 5)), delta_doy = rnorm(n_tot, 0, 4), delta_doy_resid = rnorm(n_tot, 0, 4), is_verified_anthesis = TRUE,
                fine_bedrock_depth_cm = ifelse(hab == "Flatrock", rnorm(n_tot, 22, 6), rnorm(n_tot, 185, 5)),
                min_bedrock_depth_buffer_cm = ifelse(hab == "Flatrock", rnorm(n_tot, 14, 4), rnorm(n_tot, 178, 8)),
                fine_cfvo_pct = ifelse(hab == "Flatrock", rnorm(n_tot, 38, 5), rnorm(n_tot, 6, 2)),
                ssurgo_rock_outcrop_pct = ifelse(hab == "Flatrock", runif(n_tot, 30, 65), rep(0.0, n_tot)),
                ssurgo_aws_0_50cm = ifelse(hab == "Flatrock", rnorm(n_tot, 1.8, 0.3), ifelse(hab == "Sandhills", rnorm(n_tot, 2.9, 0.4), ifelse(hab == "CoastalDune", rnorm(n_tot, 3.4, 0.5), rnorm(n_tot, 10.5, 0.9)))),
                fine_soil_ph = ifelse(hab %in% c("Flatrock", "Sandhills"), rnorm(n_tot, 4.8, 0.2), ifelse(hab == "CoastalDune", rnorm(n_tot, 6.7, 0.25), rnorm(n_tot, 6.9, 0.3))),
                fine_soil_cec = ifelse(hab == "WesternPrairie", rnorm(n_tot, 34.0, 5.0), rnorm(n_tot, 4.5, 1.2)),
                ilr_texture_1 = ifelse(hab %in% c("Sandhills", "CoastalDune"), rnorm(n_tot, 2.3, 0.2), ifelse(hab == "Flatrock", rnorm(n_tot, 1.1, 0.2), rnorm(n_tot, -0.8, 0.25))),
                ilr_texture_2 = ifelse(hab == "WesternPrairie", rnorm(n_tot, 0.9, 0.15), rnorm(n_tot, -0.2, 0.15)),
                bio4_temp_seasonality = ifelse(hab == "CoastalDune", rnorm(n_tot, 560, 18), ifelse(hab == "WesternPrairie", rnorm(n_tot, 835, 22), rnorm(n_tot, 690, 20))),
                bio6_min_temp_coldest_month = ifelse(hab == "CoastalDune", rnorm(n_tot, 4.8, 0.5), ifelse(hab == "WesternPrairie", rnorm(n_tot, -2.5, 0.6), rnorm(n_tot, 0.5, 0.5))),
                bio10_mean_temp_warmest_qtr = ifelse(hab == "WesternPrairie", rnorm(n_tot, 27.8, 0.4), rnorm(n_tot, 25.8, 0.4)),
                bio14_precip_driest_month = ifelse(hab == "WesternPrairie", rnorm(n_tot, 42, 6), rnorm(n_tot, 88, 7)), stringsAsFactors = FALSE
              )
              res_clean <- resolve_analysis_group(df_dub, intraspecific_mode = TRUE)
              u_eco <- sort(unique(res_clean$analysis_group))
              stopifnot(length(u_eco) >= 3L && all(grepl("^Ecotype_[1-4]_", u_eco)))
              df_mut <- df_dub
              df_mut$gmm_cluster <- sample(1:3, n_tot, replace = TRUE); df_mut$PC1 <- rnorm(n_tot, 50, 20)
              df_mut$cda_predicted_taxon <- sample(c("Packera anonyma", "Packera paupercula"), n_tot, replace = TRUE)
              df_mut$morphometric_outlier_flag <- sample(c(TRUE, FALSE), n_tot, replace = TRUE)
              res_mut <- resolve_analysis_group(df_mut, intraspecific_mode = TRUE)
              stopifnot(identical(res_clean$analysis_group, res_mut$analysis_group))

              cand_d <- discover_candidate_eco_predictors(res_clean, include_pheno = FALSE)
              stopifnot(all(c("lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx") %in% cand_d))
              stopifnot(identical(resolve_ecological_authority_weights(res_clean), resolve_ecological_authority_weights(res_mut)))
              srf_c <- compute_spatial_rf_mems(res_clean, cand_d, n_trees = 60, intraspecific_mode = TRUE)
              srf_m <- compute_spatial_rf_mems(res_mut, cand_d, n_trees = 60, intraspecific_mode = TRUE)
              stopifnot(max(abs(srf_c$oof_prob_matrix - srf_m$oof_prob_matrix), na.rm = TRUE) < 1e-12 && srf_c$Spatial_Block_CV_Accuracy >= 0.75)

              res_clean$rf_predicted_taxon <- srf_c$rf_predicted_taxon; res_clean$rf_predicted_prob <- srf_c$rf_predicted_prob
              res_clean$edaphic_best_fit_taxon <- srf_c$rf_predicted_taxon; res_clean$edaphic_z_score <- runif(n_tot, 0.3, 1.4)
              tri_out <- apply_taxonomic_decision_matrix(res_clean, intraspecific_mode = TRUE, gmm_two_delta_bic = 0.0)
              stopifnot(mean(tri_out$taxonomic_treatment_subtype == "Confirmed_Species") >= 0.90)
              stopifnot(tri_out$w_morph_kappa[1] >= 0.30 && tri_out$w_edaphic_kappa[1] < 0.90 && mean(tri_out$probabilistic_consensus_score) >= 0.70)
              cat("INDEPENDENT_STEP05_AND_STEP06_OK\\n")
            }})
            """
            r_file = cls._tmp / "batched_step05_zero_trust.R"
            r_file.write_text(r_code, encoding="utf-8")
            cls.r_proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step05", cwd=PROJECT_ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def _require_r(self) -> str:
        if self.r_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(self.r_proc.returncode, 0, f"Batched Step 05 zero-trust R suite failed:\n{self.r_proc.stderr}\n{self.r_proc.stdout}")
        return self.r_proc.stdout

    def test_noncircular_k1_ibe_background_similarity_and_empty_raster_hygiene(self):
        """Verify non-circular K=1 ecotype gate, IBE partial RDA, accessible M similarity, and empty-raster hygiene."""
        self.assertIn("NONCIRCULAR_K1_IBE_AND_HYGIENE_OK", self._require_r())

    def test_step05_zero_trust_scientific_and_anti_leakage_suite(self):
        """Verify physiological ENM whitelist, exact Hopkins' Law, anti-leakage imputation, Bronze/noise invariance, sterile-sheet protection, and DEM HLI/D8 TWI."""
        self.assertIn("STEP05_ZERO_TRUST_SUITE_OK", self._require_r())

    def test_enm_composite_only_cross_block_vif_and_epv_budget_guards(self):
        """Verify composite-only geochem/wetland ENM inputs, cross-block VIF <= 5.0 pruning, and anchor-count EPV predictor budgeting."""
        self.assertIn("ENM_COMPOSITE_VIF_EPV_GUARDS_OK", self._require_r())

    def test_step05_morphometric_independence_and_min_bedrock_buffer_contracts(self):
        """Static + R verification that Step 05 ENM is independent of Step 04 morphometrics and uses min() for buffered bedrock depth."""
        ana_dir = PROJECT_ROOT / "scripts" / "analysis"
        pe_src = "\n".join((ana_dir / fn).read_text(encoding="utf-8") for fn in ("step05_phenology_ecotypes.R", "step05_allochronic_pheno_utils.R", "step05_type_niche_utils.R") if (ana_dir / fn).is_file())
        ed_src = "\n".join((ana_dir / fn).read_text(encoding="utf-8") for fn in ("edaphic_fine_grain_utils.R", "edaphic_ssurgo_sda_utils.R", "edaphic_solus_3dep_utils.R", "edaphic_collinearity_contingency_utils.R", "edaphic_geochem_wetland_utils.R") if (ana_dir / fn).is_file())
        s05_src = "\n".join((ana_dir / fn).read_text(encoding="utf-8") for fn in ("05_multimodal_spatial_rf.R", "step05_env_raster_extract_utils.R", "step05_spatial_rf_cv_utils.R", "step05_consensus_niche_utils.R") if (ana_dir / fn).is_file())
        rp_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")

        rag_slice = _extract_r_function_body(pe_src, "resolve_analysis_group")
        self.assertTrue(bool(rag_slice), "resolve_analysis_group missing from Step 05 phenology/ecotype modules.")
        for bad_tok in ("Morphotype_K", "pass_pheno_gate", "symp_10km"):
            self.assertNotIn(bad_tok, rag_slice)
        for req_idx in ("lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx"):
            self.assertIn(req_idx, ed_src)

        buf_full = _extract_r_function_body(ed_src, "extract_uncertainty_buffered_predictors")
        self.assertTrue(bool(buf_full), "extract_uncertainty_buffered_predictors missing from Step 05 edaphic modules.")
        buf_slice = re.split(r"#\s*Explicit demonstration|if\s*\(\s*isTRUE\(\s*demo\s*\)\s*\)", buf_full, maxsplit=1)[0]
        self.assertIsNotNone(re.search(r"\bmin\(", buf_slice), "Buffered bedrock extractor must use min() on finite pixels.")
        self.assertNotIn("mean(x, na.rm = TRUE)", buf_slice)

        step05_sh_slice = rp_src.split("05_multimodal_spatial_rf.R")[1].split("06_ingest_genetic_data.py")[0]
        self.assertNotIn("--morphometrics", step05_sh_slice)
        self.assertNotIn("require_step04 = TRUE", s05_src)
        if rscript_available():
            self.assertIn("INDEPENDENT_STEP05_AND_STEP06_OK", self._require_r())


if __name__ == "__main__":
    unittest.main()
