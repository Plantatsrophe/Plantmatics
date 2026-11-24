#!/usr/bin/env python3
"""Batched warm-R and Python SQLite contract suite for Step 05 taxonomic clarity invariants."""

from __future__ import annotations

import json, math, re, sqlite3, sys, tempfile, unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import persist_spatial_rf_predictions
from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION, init_db_schema
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

RUN_PY, RUN_R = "run_s05_clarity_py", "run_s05_clarity_r"
SKIP_R_REASON = "Skipping warm-R Step 05 taxonomic clarity tests under fast audit or missing Rscript."
SKIP_R_COND = not rscript_available()


def _is_finite(val: object) -> bool:
    try:
        return val is not None and not isinstance(val, bool) and math.isfinite(float(val))
    except (TypeError, ValueError):
        return False


def _seed_clarity_db(tmp_p: Path) -> Path:
    db_path = tmp_p / "step05_clarity.sqlite"
    with sqlite3.connect(str(db_path)) as conn:
        init_db_schema(conn)
        conn.executemany("INSERT INTO pipeline_runs (run_id) VALUES (?);", [(RUN_PY,), (RUN_R,)])
        conn.executemany(
            "INSERT INTO vouchers (catalogNumber,institutionCode,scientificName,species_standardized,determiner_tier,type_status,typified_basionym,concept_synonymy_status,decimalLatitude,decimalLongitude,coordinateUncertainty,doy,is_primary_duplicate,is_centroid_or_inst_coord,is_open_or_hybrid_concept) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);",
            [("PY_RT_01", "NCU", "Packera dubia", "Packera dubia", "Tier_1_Gold", "Holotype", "Cineraria dubia", "Homotypic_Synonym", 35.80, -79.10, 120.0, 110.0, 1, 0, 0),
             ("R_RT_01", "NCU", "Packera dubia", "Packera dubia", "Tier_1_Gold", "Holotype", "Cineraria dubia", "Homotypic_Synonym", 35.82, -79.12, 120.0, 112.0, 1, 0, 0)],
        )
        conn.executemany("INSERT INTO plant_individuals (catalogNumber,plant_individual_id) VALUES (?,?);", [("PY_RT_01", 0), ("R_RT_01", 0)])
        conn.executemany(
            "INSERT INTO run_active_cohort (run_id,catalogNumber,plant_individual_id,species_standardized,determiner_tier,type_status,is_type_specimen,authority_weight,spatial_balance_weight) VALUES (?,?,?,?,?,?,?,?,?);",
            [(RUN_PY, "PY_RT_01", 0, "Packera dubia", "Tier_1_Gold", "Holotype", 1, 1.0, 1.0), (RUN_R, "R_RT_01", 0, "Packera dubia", "Tier_1_Gold", "Holotype", 1, 1.0, 1.0)],
        )
        conn.commit()
    return db_path


_WARM_R_SCRIPT = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE); out_json <- args[1L]; tmp_root <- args[2L]; db_path <- args[3L]
Sys.setenv(PLANTMATICS_SKIP_FIGURES = "1", PACKERA_SKIP_FIGURES = "1")
suppressPackageStartupMessages({ library(jsonlite); library(terra) })
if (!exists("compute_spatial_rf_mems", mode = "function")) for (.f in c("scripts/analysis/05_multimodal_spatial_rf.R", "scripts/analysis/step05_sdm_maxent.R", "scripts/analysis/edaphic_geochem_wetland_utils.R")) if (file.exists(.f)) source(.f)
res <- list(); empty_env <- file.path(tmp_root, "empty_env"); dir.create(empty_env, recursive = TRUE, showWarnings = FALSE)
ev_c1 <- c("fine_soil_ph", "fine_bedrock_depth_cm", "bio1_temp_mean", "bio12_precip_annual")

# 1. Concept/quarantine integration, rooftop/centroid/coarse type gates, & POSIX .norm_name_tokens()
el_df <- data.frame(
  catalogNumber = c("EL_HOLO", "EL_CENT", "EL_18KM", "EL_HYB", "EL_60KM"), species_standardized = "Packera dubia", determiner_tier = "Tier_1_Gold",
  type_status = c("Holotype", "Holotype", "Holotype", "", ""), is_type_specimen = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  concept_synonymy_status = c("Homotypic_Synonym", "Homotypic_Synonym", "Homotypic_Synonym", "Hybrid_Formula", "Accepted_Sensu_Stricto"),
  is_open_or_hybrid_concept = c(FALSE, FALSE, FALSE, TRUE, FALSE), is_centroid_or_inst_coord = c(FALSE, TRUE, FALSE, FALSE, FALSE),
  coordinateUncertaintyInMeters = c(150.0, 200.0, 18000.0, 150.0, 60000.0),
  latitude = c(35.5, 35.8, 35.7, 35.6, 35.4), longitude = c(-80.5, -79.5, -79.8, -80.1, -80.3), year = 2020L, stringsAsFactors = FALSE
)
el_out <- is_fine_edaphic_type_or_anchor_eligible(el_df); fc_el <- flag_coordinate_uncertainty(el_df)
el_roof <- data.frame(catalogNumber = "EL_ROOF", species_standardized = "Packera dubia", determiner_tier = "Tier_1_Gold", type_status = "Holotype",
  is_type_specimen = TRUE, coordinateUncertaintyInMeters = 100.0, latitude = 35.9108, longitude = -79.0505, decimalLatitude = 35.9108,
  decimalLongitude = -79.0505, locality = "Chapel Hill, NC", year = 2020L, stringsAsFactors = FALSE)
roof_el <- is_fine_edaphic_type_or_anchor_eligible(el_roof); roof_q <- detect_r_centroid_or_inst_quarantine(el_roof, is_type_vec = TRUE, exempt_types = FALSE)

set.seed(501L); n_cln <- 24L; sp_cln <- rep(c("Packera aurea", "Packera dubia"), each = 12L); is_d <- sp_cln == "Packera dubia"
lat_cln <- c(seq(35.0, 37.2, length.out = 12L), seq(34.2, 36.4, length.out = 12L)); lon_cln <- c(seq(-82.8, -80.6, length.out = 12L), seq(-80.8, -78.6, length.out = 12L))
df_base <- data.frame(
  catalogNumber = c("A_TYPE_CLEAN", sprintf("A_%02d", 2:12), "D_TYPE_CLEAN", sprintf("D_%02d", 2:12)), plant_individual_id = 0L,
  species_standardized = sp_cln, analysis_group = sp_cln, determiner_tier = "Tier_1_Gold",
  type_status = ifelse(seq_len(n_cln) %in% c(1L, 13L), "Holotype", ""), is_type_specimen = seq_len(n_cln) %in% c(1L, 13L),
  typified_basionym = ifelse(seq_len(n_cln) == 1L, "Senecio aureus", ifelse(seq_len(n_cln) == 13L, "Cineraria dubia", "")),
  concept_synonymy_status = ifelse(seq_len(n_cln) %in% c(1L, 13L), "Homotypic_Synonym", "Accepted_Sensu_Stricto"),
  is_open_or_hybrid_concept = FALSE, is_centroid_or_inst_coord = FALSE, coordinateUncertaintyInMeters = 150.0,
  latitude = lat_cln, longitude = lon_cln, decimalLatitude = lat_cln, decimalLongitude = lon_cln, locality = "Wild outcrop 5 km N of creek", year = 2020L,
  fine_soil_ph = ifelse(is_d, rnorm(n_cln, 4.75, 0.14), rnorm(n_cln, 6.20, 0.14)), fine_bedrock_depth_cm = ifelse(is_d, rnorm(n_cln, 22.0, 3.5), rnorm(n_cln, 155.0, 9.0)),
  bio1_temp_mean = ifelse(is_d, rnorm(n_cln, 16.4, 0.30), rnorm(n_cln, 12.1, 0.30)), bio12_precip_annual = ifelse(is_d, rnorm(n_cln, 1180.0, 35.0), rnorm(n_cln, 1390.0, 35.0)),
  doy = round(98.0 + 4.0 * (lat_cln - 34.0) + ifelse(is_d, -10.0, 10.0)), voucher_phenological_state = "Peak_Anthesis", n_anthesis = 8L, n_fruiting = 1L,
  total_capitula = 9L, is_verified_anthesis = TRUE, sgmc_lithology_class = ifelse(is_d, "Felsic_Granitic_Gneiss", "Unconsolidated_Sedimentary"),
  is_rock_outcrop = is_d, outcrop_class = ifelse(is_d, "Granitic flatrock", "Non-outcrop matrix"), stringsAsFactors = FALSE
)
extra_c1 <- df_base[c(14L, 1L, 13L, 1L), , drop = FALSE]; extra_c1$catalogNumber <- c("HYB_OPEN_01", "TYPE_CENT", "TYPE_15KM", "TYPE_ROOF")
extra_c1$type_status <- c("", "Holotype", "Isotype", "Holotype"); extra_c1$is_type_specimen <- c(FALSE, TRUE, TRUE, TRUE)
extra_c1$concept_synonymy_status <- c("Hybrid_Formula", "Homotypic_Synonym", "Homotypic_Synonym", "Homotypic_Synonym")
extra_c1$is_open_or_hybrid_concept <- c(TRUE, FALSE, FALSE, FALSE); extra_c1$is_centroid_or_inst_coord <- c(FALSE, TRUE, FALSE, TRUE)
extra_c1$coordinateUncertaintyInMeters <- c(150.0, 200.0, 15000.0, 100.0); extra_c1$locality <- c("Roadside", "Orange County", "Piedmont", "Chapel Hill, NC")
df_c1 <- rbind(df_base, extra_c1); df_c1_mut <- df_c1
df_c1_mut$fine_soil_ph[25L] <- 9.8; df_c1_mut$fine_bedrock_depth_cm[25L] <- 1.0; df_c1_mut$bio1_temp_mean[25L] <- -12.0; df_c1_mut$bio12_precip_annual[25L] <- 4200.0
df_c1_mut$doy[25L] <- 215.0; df_c1_mut$sgmc_lithology_class[25L] <- "Ultramafic_Serpentine"; df_c1_mut$is_rock_outcrop[25L] <- TRUE; df_c1_mut$outcrop_class[25L] <- "Serpentine"
rf_c1a <- compute_spatial_rf_mems(df_c1, ev_c1, n_trees = 25L); rf_c1b <- compute_spatial_rf_mems(df_c1_mut, ev_c1, n_trees = 25L)
sdm_c1a <- suppressWarnings(fit_and_project_maxent_sdm_maps(df_c1, rf_c1a, empty_env, file.path(tmp_root, "sdm_c1a"), kept_vars = ev_c1, focal_taxon = "Packera dubia", demo = FALSE))
sdm_c1b <- suppressWarnings(fit_and_project_maxent_sdm_maps(df_c1_mut, rf_c1b, empty_env, file.path(tmp_root, "sdm_c1b"), kept_vars = ev_c1, focal_taxon = "Packera dubia", demo = FALSE))
ph_c1a <- compute_phenological_anomalies(df_c1, demo = FALSE); ph_c1b <- compute_phenological_anomalies(df_c1_mut, demo = FALSE)
ct_c1a <- evaluate_multiclass_substrate_contingency(df_c1)$summary_table; ct_c1b <- evaluate_multiclass_substrate_contingency(df_c1_mut)$summary_table
tok_var <- .norm_name_tokens("Packera paupercula var. savannarum (R.R.Kowal) Trock"); tok_ssp <- .norm_name_tokens("Packera streptanthifolia ssp. rubricaulis (Greene) Bain")
tok_hyb1 <- .norm_name_tokens("Packera ×memmingeri (Britton ex Small) Weakley"); tok_hyb2 <- .norm_name_tokens("Packera × memmingeri")
res$c1 <- list(
  el_edaph = as.logical(el_out$edaphic_ok), el_clim = as.logical(el_out$clim_ok), el_open = as.logical(el_out$is_open_hyb), el_qtype = as.logical(el_out$is_quarantined_type),
  fc_edaph = as.logical(fc_el$is_edaphic_eligible_coord[c(1L, 2L, 3L, 5L)]), roof_blocked = !isTRUE(roof_el$edaphic_ok[1L]) && isTRUE(roof_el$is_quarantined_type[1L]) && isTRUE(roof_q$is_quarantined[1L]),
  src_hyb = as.character(rf_c1a$rf_prediction_source[25L]), src_qtypes = as.character(rf_c1a$rf_prediction_source[26:28]),
  clean_types_active = all(c(1L, 13L) %in% rf_c1a$anchor_idx) && !any(25:28 %in% rf_c1a$anchor_idx),
  rf_oof_diff = max(abs(rf_c1a$oof_prob_matrix[1:24, ] - rf_c1b$oof_prob_matrix[1:24, ]), na.rm = TRUE),
  sdm_npres_eq = identical(sdm_c1a$performance$N_Presences, sdm_c1b$performance$N_Presences), sdm_auc_diff = max(abs(sdm_c1a$performance$AUC_Block_Ensemble - sdm_c1b$performance$AUC_Block_Ensemble), na.rm = TRUE),
  hop_n_eq = identical(as.integer(ph_c1a$baseline$n), as.integer(ph_c1b$baseline$n)), hop_slope_diff = abs(as.numeric(ph_c1a$baseline$slope) - as.numeric(ph_c1b$baseline$slope)),
  hop_ddoy_diff = max(abs(ph_c1a$df$delta_doy[1:24] - ph_c1b$df$delta_doy[1:24]), na.rm = TRUE), ct_tot_eq = identical(ct_c1a$Total_Vouchers, ct_c1b$Total_Vouchers), ct_out_eq = identical(ct_c1a$Outcrop_Count, ct_c1b$Outcrop_Count),
  tok_var_ok = "packera paupercula var. savannarum" %in% tok_var, tok_ssp_ok = "packera streptanthifolia subsp. rubricaulis" %in% tok_ssp,
  tok_hyb_match = any(grepl("^packera[[:space:]]+([×x][[:space:]]*)?memmingeri$", tok_hyb1)) && length(intersect(tok_hyb1, tok_hyb2)) >= 1L
)

# 2. Basionym-anchored locus classicus control & heterotypic synonym divergence audit
het_row <- df_base[13L, , drop = FALSE]; het_row$catalogNumber <- "HET_TYPE_DIV"; het_row$typified_basionym <- "Senecio alabamensis"; het_row$concept_synonymy_status <- "Heterotypic_Synonym"
het_row$latitude <- 32.5; het_row$longitude <- -87.5; het_row$decimalLatitude <- 32.5; het_row$decimalLongitude <- -87.5; het_row$fine_soil_ph <- 8.4; het_row$fine_bedrock_depth_cm <- 198.0
het_row$bio1_temp_mean <- 21.0; het_row$bio12_precip_annual <- 950.0; het_row$sgmc_lithology_class <- "Unconsolidated_Sedimentary"; het_row$is_rock_outcrop <- FALSE; het_row$outcrop_class <- "Non-outcrop matrix"
df_het <- rbind(df_base, het_row)
tc_hom <- compute_type_locality_ecological_control(df_base, primary_env_vars = ev_c1, primary_edaph_vars = c("fine_soil_ph", "fine_bedrock_depth_cm"))
tc_het <- compute_type_locality_ecological_control(df_het, primary_env_vars = ev_c1, primary_edaph_vars = c("fine_soil_ph", "fine_bedrock_depth_cm"))
df_het$is_heterotypic_type_audit <- tc_het$is_heterotypic_type_audit; rf_het <- compute_spatial_rf_mems(df_het, ev_c1, n_trees = 25L)
cons_het <- execute_crossmodal_consensus(df_het, rf_res = rf_het, primary_env_vars = ev_c1, primary_edaph_vars = c("fine_soil_ph", "fine_bedrock_depth_cm"), primary_clim_vars = c("bio1_temp_mean", "bio12_precip_annual"))
res$c2 <- list(
  hom_z_shift = max(abs(tc_hom$type_edaphic_z_score[13:24] - tc_het$type_edaphic_z_score[13:24]), na.rm = TRUE), hom_holo_z = as.numeric(tc_het$type_edaphic_z_score[13L]),
  het_vs_hom_z = as.numeric(tc_het$type_edaphic_z_score[25L]), tc_het_st = as.character(tc_het$type_niche_status[25L]), tc_het_hint = as.character(tc_het$rf_prediction_source_hint[25L]),
  het_rf_src = as.character(rf_het$rf_prediction_source[25L]), het_audit_flag = isTRUE(tc_het$is_heterotypic_type_audit[25L]), het_conf_flag = isTRUE(cons_het$multimodal_conflict_flag[25L]),
  bas_hom_a = as.character(cons_het$nearest_type_niche_basionym[1L]), bas_hom_d = as.character(cons_het$nearest_type_niche_basionym[13L]), bas_het = as.character(cons_het$nearest_type_niche_basionym[25L]),
  bas_z_het = as.numeric(cons_het$nearest_type_basionym_edaphic_z[25L]), bas_z_hom = as.numeric(cons_het$nearest_type_basionym_edaphic_z[13L])
)

# 3. Pedology, hardpan ignore, smooth indices, point-mass screening, & coastal/Canadian NA cache
ss_tab <- data.frame(
  gid = c(1L, 2L, 3L, 3L), mukey = c("101", "201", "301", "301"), musym = c("RnC", "LoC", "WkB", "WkB"),
  muname = c("Rion sandy loam, 8 to 15 percent slopes", "Louisburg loamy sand, 6 to 15 percent slopes", "Wake-Rock outcrop complex, 2 to 15 percent slopes", "Wake-Rock outcrop complex, 2 to 15 percent slopes"),
  brockdepmin = c(160.0, 90.0, 18.0, 18.0), aws050wta = c(6.5, 4.5, 2.2, 2.2), cokey = c("101_1", "201_1", "301_1", "301_2"), compname = c("Rion", "Louisburg", "Wake", "Rock outcrop"),
  comppct_r = c(85.0, 80.0, 65.0, 25.0), compkind = c("Series", "Series", "Series", "Miscellaneous area"), majcompflag = "Yes",
  taxsuborder = c("Udults", "Udepts", "Udipsamments", NA), taxgreatgroup = c("Hapludults", "Dystrudepts", "Quartzipsamments", NA),
  hzdept_r = c(0.0, 0.0, 0.0, NA), hzdepb_r = c(15.0, 15.0, 15.0, NA), sandtotal_r = c(68.0, 74.0, 78.0, NA), claytotal_r = c(14.0, 10.0, 7.0, NA),
  ph1to1h2o_r = c(5.3, 5.1, 4.8, NA), om_r = c(1.5, 1.4, 1.8, NA), cec7_r = c(5.5, 4.8, 3.9, NA), awc_r = c(0.11, 0.09, 0.05, NA),
  resdept_r = c(NA, 90.0, 18.0, 0.0), reskind = c("", "Paralithic bedrock", "Lithic bedrock", "Lithic bedrock"), stringsAsFactors = FALSE
)
agg_ss <- aggregate_ssurgo_tabular_to_points(data.frame(gid = 1:3, catalogNumber = c("SS_RION", "SS_LOU", "SS_WAKE"), latitude = c(35.8, 35.9, 35.7), longitude = c(-79.1, -79.2, -79.3), stringsAsFactors = FALSE), ss_tab)
cand_ex <- discover_candidate_eco_predictors(data.frame(
  catalogNumber = paste0("EX_", 1:12), fine_soil_ph = seq(4.6, 6.8, length.out = 12L), fine_bedrock_depth_cm = seq(20.0, 180.0, length.out = 12L),
  bio1_temp_mean = seq(10.0, 18.0, length.out = 12L), bio12_precip_annual = seq(950.0, 1450.0, length.out = 12L),
  ssurgo_rock_outcrop_pct = seq(0.0, 65.0, length.out = 12L), ssurgo_aws_0_50cm = seq(2.0, 11.0, length.out = 12L), ssurgo_lithic_series_pct = seq(0.0, 80.0, length.out = 12L), stringsAsFactors = FALSE
), include_pheno = FALSE)

hp_env <- file.path(tmp_root, "hp_env"); dir.create(file.path(hp_env, "solus_100m"), recursive = TRUE, showWarnings = FALSE)
r_hp <- terra::rast(xmin = -80.1, xmax = -79.9, ymin = 35.9, ymax = 36.1, nrows = 40, ncols = 40, crs = "EPSG:4326")
r_lith <- r_hp; terra::values(r_lith) <- 175.0; r_pan <- r_hp; terra::values(r_pan) <- 22.0
terra::writeRaster(r_lith, file.path(hp_env, "solus_100m", "anylithicdpt_cm_p.tif"), overwrite = TRUE); terra::writeRaster(r_pan, file.path(hp_env, "solus_100m", paste0("resdept_", "all_cm_p.tif")), overwrite = TRUE)
pt_hp <- data.frame(catalogNumber = "HP_1", latitude = 36.0, longitude = -80.0, coordinateUncertaintyInMeters = 200.0, stringsAsFactors = FALSE)
hp_sol <- extract_polaris_and_solus_predictors(pt_hp, env_dir = hp_env, demo = FALSE); hp_buf <- extract_uncertainty_buffered_predictors(pt_hp, env_dir = hp_env, demo = FALSE)
loi_env <- file.path(tmp_root, "loi_env"); dir.create(file.path(loi_env, "solus_100m"), recursive = TRUE, showWarnings = FALSE)
r_lt <- terra::rast(xmin = -82.0, xmax = -78.0, ymin = 34.0, ymax = 38.0, nrows = 40, ncols = 40, crs = "EPSG:4326"); xy_l <- terra::xyFromCell(r_lt, seq_len(terra::ncell(r_lt)))
r_b1 <- r_lt; terra::values(r_b1) <- 18.0 - 0.8 * (xy_l[, 2] - 34.0); names(r_b1) <- "bio1_temp_mean"; r_b12 <- r_lt; terra::values(r_b12) <- 1100.0 + 40.0 * (xy_l[, 1] + 82.0); names(r_b12) <- "bio12_precip_annual"
r_bd_l <- r_lt; terra::values(r_bd_l) <- 20.0 + 40.0 * (xy_l[, 1] + 82.0); names(r_bd_l) <- "anylithicdpt"; r_cf_l <- r_lt; terra::values(r_cf_l) <- 60.0 - 12.0 * (xy_l[, 1] + 82.0); names(r_cf_l) <- "fragvol"
terra::writeRaster(r_b1, file.path(loi_env, "wc2.1_bio_1.tif"), overwrite = TRUE); terra::writeRaster(r_b12, file.path(loi_env, "wc2.1_bio_12.tif"), overwrite = TRUE)
terra::writeRaster(r_bd_l, file.path(loi_env, "solus_100m", "anylithicdpt_cm_p.tif"), overwrite = TRUE); terra::writeRaster(r_cf_l, file.path(loi_env, "solus_100m", "fragvol_5_cm_p.tif"), overwrite = TRUE)
lon_l <- seq(-81.5, -78.5, length.out = 8L); lat_l <- seq(34.5, 37.5, length.out = 8L); pts_l <- terra::vect(data.frame(longitude = lon_l, latitude = lat_l), geom = c("longitude", "latitude"), crs = "EPSG:4326")
df_loi <- compute_derived_habitat_predictors(data.frame(
  catalogNumber = paste0("LOI_", 1:8), longitude = lon_l, latitude = lat_l, coordinateUncertaintyInMeters = 120.0,
  fine_bedrock_depth_cm = as.numeric(terra::extract(r_bd_l, pts_l)[[2L]]), fine_cfvo_pct = as.numeric(terra::extract(r_cf_l, pts_l)[[2L]]),
  bio1_temp_mean = as.numeric(terra::extract(r_b1, pts_l)[[2L]]), bio12_precip_annual = as.numeric(terra::extract(r_b12, pts_l)[[2L]]),
  is_rock_outcrop = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE), is_nwi_wetland = c(1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L),
  nwi_wetland_system = c("PEM1A", "PFO1A", "UPL", "UPL", "PSS1C", "UPL", "UPL", "UPL"), nwi_wetland_dist_m = c(0.0, 25.0, 300.0, 450.0, 10.0, 180.0, 500.0, 600.0),
  topo_twi_30m = c(11.5, 10.0, 5.2, 4.5, 10.8, 6.0, 4.8, 4.0), stringsAsFactors = FALSE
))
st_loi <- .build_regional_sdm_raster_stack(df_loi, env_dir = loi_env, kept_vars = c("bio1_temp_mean", "bio12_precip_annual", "lithic_outcrop_index"), demo = FALSE)
r_loi_pts <- as.numeric(terra::extract(st_loi[["lithic_outcrop_index"]], pts_l)[[2L]])
exp_loi <- 0.80 * pmin(1.0, pmax(0.0, 1.0 - log1p(pmax(df_loi$fine_bedrock_depth_cm, 0.0)) / log1p(200.0))) + 0.20 * (log1p(pmax(df_loi$fine_cfvo_pct, 0.0)) / log1p(100.0))
twi_rng <- range(df_loi$topo_twi_30m); exp_wsi <- 0.70 * ifelse(grepl("^(PEM|PFO|PSS)", df_loi$nwi_wetland_system, ignore.case = TRUE), exp(-df_loi$nwi_wetland_dist_m / 90.0), 0.45 * exp(-df_loi$nwi_wetland_dist_m / 90.0)) + 0.30 * ((df_loi$topo_twi_30m - twi_rng[1L]) / diff(twi_rng))
set.seed(502L); X_pm <- cbind(fine_bedrock_depth_cm = c(rep(200.0, 26L), seq(15.0, 45.0, length.out = 14L)), fine_soil_ph = rnorm(40L, 5.8, 0.4), bio1_temp_mean = rnorm(40L, 14.0, 1.2))
pm_res <- screen_edaphic_point_mass_columns(X_pm, max_tie_frac = 0.25)
pm_kept <- if (is.logical(pm_res)) colnames(X_pm)[pm_res] else if (is.list(pm_res)) as.character(pm_res[[intersect(c("kept_cols", "kept_vars", "selected_vars", "kept"), names(pm_res))[1L]]]) else as.character(pm_res)
res_pan <- resolve_analysis_group(data.frame(
  catalogNumber = sprintf("PAN_%02d", 1:40), species_standardized = "Packera aurea", determiner_tier = "Tier_1_Gold",
  latitude = seq(35.0, 39.0, length.out = 40L), longitude = seq(-82.0, -78.0, length.out = 40L), coordinateUncertaintyInMeters = 150.0,
  fine_bedrock_depth_cm = c(rep(200.0, 26L), runif(14L, 18.0, 35.0)), is_rock_outcrop = c(rep(FALSE, 26L), rep(TRUE, 14L)), ssurgo_rock_outcrop_pct = c(rep(0.0, 26L), rep(45.0, 14L)),
  fine_soil_ph = rnorm(40L, 5.7, 0.35), ilr_texture_1 = rnorm(40L, 0.2, 0.30), ilr_texture_2 = rnorm(40L, -0.1, 0.25), bio1_temp_mean = rnorm(40L, 13.5, 0.9), bio12_precip_annual = rnorm(40L, 1250.0, 75.0), stringsAsFactors = FALSE
), intraspecific_mode = TRUE)
can_env <- file.path(tmp_root, "can_env"); for (.sd in c("dem_3dep_30m", "solus_100m")) dir.create(file.path(can_env, .sd), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(r_lith, file.path(can_env, "dem_3dep_30m", "N35W080.tif"), overwrite = TRUE); terra::writeRaster(r_lith, file.path(can_env, "solus_100m", "anylithicdpt_cm_p.tif"), overwrite = TRUE)
ss_calls <- 0L; mock_ss_na <- function(d_in, cache_dir = NULL, batch_size = 75, demo = FALSE) { ss_calls <<- ss_calls + 1L; d_in$muname <- NA_character_; d_in$is_rock_outcrop <- FALSE; d_in$outcrop_class <- "Non-outcrop matrix"; d_in$fine_soil_ph <- NA_real_; d_in$ssurgo_ph <- NA_real_; d_in$ilr_texture_1 <- NA_real_; d_in$ilr_texture_2 <- NA_real_; d_in }
df_can <- data.frame(catalogNumber = c("COAST_NA", "CAN_NA"), latitude = c(34.2104, 46.5000), longitude = c(-77.8868, -75.0000), coordinateUncertaintyInMeters = 100.0, stringsAsFactors = FALSE); can_csv <- file.path(can_env, "fine_grain_edaphic_cache.csv")
load_and_update_fine_grain_cache(df_can, env_dir = can_env, demo = FALSE, query_ssurgo_fn = mock_ss_na, cache_path = can_csv); calls_1 <- ss_calls
msg_can2 <- capture.output({ load_and_update_fine_grain_cache(df_can, env_dir = can_env, demo = FALSE, query_ssurgo_fn = mock_ss_na, cache_path = can_csv) }, type = "message")
res$c3 <- list(
  rion_pct = as.numeric(agg_ss$ssurgo_rock_outcrop_pct[1L]), lou_pct = as.numeric(agg_ss$ssurgo_rock_outcrop_pct[2L]), wake_out_pct = as.numeric(agg_ss$ssurgo_rock_outcrop_pct[3L]), wake_lith_pct = as.numeric(agg_ss$ssurgo_lithic_series_pct[3L]),
  excl_ro = !("ssurgo_rock_outcrop_pct" %in% cand_ex), excl_aws = !("ssurgo_aws_0_50cm" %in% cand_ex), keep_ph = "fine_soil_ph" %in% cand_ex,
  hp_sol_bd = as.numeric(hp_sol$fine_bedrock_depth_cm[1L]), hp_buf_bd = as.numeric(hp_buf$min_bedrock_depth_buffer_cm[1L]),
  loi_v_max = max(df_loi$lithic_outcrop_index), loi_r_max = max(r_loi_pts, na.rm = TRUE), loi_formula_err = max(abs(df_loi$lithic_outcrop_index - exp_loi)), loi_vr_err = max(abs(df_loi$lithic_outcrop_index - r_loi_pts), na.rm = TRUE),
  wsi_formula_err = max(abs(df_loi$wetland_seepage_index - exp_wsi)), pm_dropped_bd = !("fine_bedrock_depth_cm" %in% pm_kept) && ("fine_soil_ph" %in% pm_kept),
  pan_grps = unique(as.character(res_pan$analysis_group)), can_calls_1 = as.integer(calls_1), can_calls_2 = as.integer(ss_calls - calls_1), can_zero_refresh = any(grepl("0 coordinates requiring refresh", msg_can2, fixed = TRUE))
)

# 4. Haversine ocean-box exclusion, fold-safe SDM CV, fruiting gate, allochronic bimodality, & R SQLite persist
oc_env <- file.path(tmp_root, "ocean_env"); dir.create(oc_env, recursive = TRUE, showWarnings = FALSE)
r_oc <- terra::rast(xmin = -82.0, xmax = -74.0, ymin = 32.0, ymax = 38.0, nrows = 48, ncols = 64, crs = "EPSG:4326"); xy_o <- terra::xyFromCell(r_oc, seq_len(terra::ncell(r_oc))); is_land <- xy_o[, 1] <= -76.5
r_ob1 <- r_oc; terra::values(r_ob1) <- ifelse(is_land, 17.5 - 0.5 * (xy_o[, 2] - 32.0), NA_real_); names(r_ob1) <- "bio1_temp_mean"
r_ob12 <- r_oc; terra::values(r_ob12) <- ifelse(is_land, 1250.0 + 20.0 * (xy_o[, 1] + 80.0), NA_real_); names(r_ob12) <- "bio12_precip_annual"
terra::writeRaster(r_ob1, file.path(oc_env, "wc2.1_bio_1.tif"), overwrite = TRUE); terra::writeRaster(r_ob12, file.path(oc_env, "wc2.1_bio_12.tif"), overwrite = TRUE)
pres_xy <- cbind(longitude = seq(-77.4, -76.7, length.out = 10L), latitude = seq(34.0, 36.0, length.out = 10L))
pres_X <- cbind(bio1_temp_mean = seq(16.5, 15.5, length.out = 10L), bio12_precip_annual = seq(1300.0, 1315.0, length.out = 10L))
bg_oc <- .generate_accessible_background_env(pres_xy, pres_X, n_bg = 120L, seed = 42L, env_dir = oc_env, demo = FALSE)
max_near_km <- max(apply(compute_haversine_dist_km_matrix(bg_oc$coords, pres_xy), 1L, min))
set.seed(504L); n_sc <- 40L; sp_sc <- rep(c("Packera aurea", "Packera dubia"), each = 20L); b1_sc <- rnorm(n_sc, 14.0, 1.0); bd_sc <- rnorm(n_sc, 90.0, 20.0); bd_sc[seq(1L, n_sc, by = 3L)] <- NA_real_
sdm_sc <- suppressWarnings(fit_and_project_maxent_sdm_maps(data.frame(
  catalogNumber = sprintf("SC_%02d", seq_len(n_sc)), species_standardized = sp_sc, analysis_group = sp_sc, determiner_tier = "Tier_1_Gold", coordinateUncertaintyInMeters = 150.0,
  latitude = rep(seq(34.0, 38.0, length.out = 20L), 2L), longitude = rep(seq(-82.0, -78.0, length.out = 20L), 2L),
  bio1_temp_mean = b1_sc, bio1_dup = b1_sc * 1.005 + rnorm(n_sc, 0, 0.001), bio12_precip_annual = rnorm(n_sc, 1250.0, 60.0), fine_bedrock_depth_cm = bd_sc, stringsAsFactors = FALSE
), NULL, empty_env, file.path(tmp_root, "sdm_sc"), kept_vars = c("bio1_temp_mean", "bio1_dup", "bio12_precip_annual", "fine_bedrock_depth_cm"), focal_taxon = "Packera aurea", demo = FALSE))
lat_a <- seq(33.0, 37.5, length.out = 12L)
ph_anth <- compute_phenological_anomalies(data.frame(
  catalogNumber = c(sprintf("PK_%02d", 1:12), "FR_DOM_01", "FR_DOM_02"), species_standardized = "Packera dubia", determiner_tier = "Tier_1_Gold",
  decimalLatitude = c(lat_a, 35.0, 35.2), decimalLongitude = c(seq(-81.5, -78.5, length.out = 12L), -80.0, -79.8), doy = c(round(95.0 + 4.0 * (lat_a - 33.0)), 195.0, 198.0),
  n_anthesis = c(rep(8L, 12L), 1L, 1L), n_fruiting = c(rep(1L, 12L), 25L, 25L), total_capitula = c(rep(9L, 12L), 26L, 26L),
  voucher_phenological_state = c(rep("Peak_Anthesis", 12L), "anthesis", "anthesis"), reproductiveCondition = c(rep("fl", 12L), "fl; fr", "fl; fr"), stringsAsFactors = FALSE
), demo = FALSE)
set.seed(505L); bm_res <- test_allochronic_bimodality_von_mises(c(rnorm(20L, -30.0, 4.0), rnorm(20L, 30.0, 4.0))); um_res <- test_allochronic_bimodality_von_mises(rnorm(40L, 0.0, 5.0), weights = rep(1.0, 40L))
Sys.setenv(PLANTMATICS_DB_PATH = db_path, PACKERA_DB_PATH = db_path, RUN_ID = "run_s05_clarity_r")
ok_r_sql <- persist_step05_to_sqlite(data.frame(
  catalogNumber = "R_RT_01", plant_individual_id = 0L, latitude = 35.82, longitude = -79.12, ecotype_regime = "Ecotype_1_Lithic_High",
  lithic_outcrop_index = 0.7425, summer_moisture_deficit_idx = 1.85, maritime_continentality_idx = 2.42, hopkins_adjusted_doy = 109.5, hopkins_adjusted_doy_resid = -6.25,
  nearest_type_niche_basionym = "Cineraria dubia", nearest_type_basionym_edaphic_z = 0.415, anthesis_dominance_weight = 0.875, allochronic_bimodality_delta_bic = 14.32,
  rf_prediction_source = "Spatial_Block_OOF_Anchor", rf_predicted_taxon = "Packera dubia", rf_predicted_prob = 0.93, stringsAsFactors = FALSE
), vouchers_path = NULL, run_id = "run_s05_clarity_r", demo = FALSE)
Sys.unsetenv(c("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH", "RUN_ID"))
res$c4 <- list(
  bg_n = nrow(bg_oc$coords), bg_all_land = all(bg_oc$coords[, 1L] <= -76.5), bg_max_km = as.numeric(max_near_km), sdm_auc_ens = as.numeric(sdm_sc$performance$AUC_Block_Ensemble[1L]),
  has_adw = "anthesis_dominance_weight" %in% names(ph_anth$df), has_hresid = "hopkins_adjusted_doy_resid" %in% names(ph_anth$df),
  fr_ver = as.logical(ph_anth$df$is_verified_anthesis[13:14]), pk_ver = all(as.logical(ph_anth$df$is_verified_anthesis[1:12])),
  fr_adw = as.numeric(ph_anth$df$anthesis_dominance_weight[13:14]), pk_adw_min = min(as.numeric(ph_anth$df$anthesis_dominance_weight[1:12])),
  anth_base_n = as.integer(ph_anth$baseline$n), hresid_fin = all(is.finite(as.numeric(ph_anth$df$hopkins_adjusted_doy_resid[1:12]))),
  bm_dbic = as.numeric(bm_res$allochronic_bimodality_delta_bic), bm_ashman = as.numeric(bm_res$ashman_d), um_dbic = as.numeric(um_res$allochronic_bimodality_delta_bic), ok_r_sql = isTRUE(ok_r_sql)
)
jsonlite::write_json(res, out_json, auto_unbox = TRUE, na = "null", digits = 10)
"""


class TestStep05TaxonomicClarity(unittest.TestCase):
    """Batched warm-R, Python SQLite, and static contract suite for Step 05 taxonomic clarity."""

    results: dict = {}

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp_dir = tempfile.TemporaryDirectory(prefix="step05_clarity_")
        cls.tmp_p = Path(cls._tmp_dir.name); cls.db_path = _seed_clarity_db(cls.tmp_p); cls.results = {}
        if SKIP_R_COND:
            return
        r_path, out_json = cls.tmp_p / "run_step05_clarity.R", cls.tmp_p / "clarity_results.json"
        r_path.write_text(_WARM_R_SCRIPT, encoding="utf-8")
        proc = run_warm_r_cmd(["Rscript", str(r_path), str(out_json), str(cls.tmp_p), str(cls.db_path)], stage_key="step05", cwd=PROJECT_ROOT, skip_figures=True)
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(f"Batched Step 05 clarity R failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls.results = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_01_concept_quarantine_type_coord_and_posix_norm_tokens(self) -> None:
        c1 = self.results["c1"]
        self.assertEqual(([bool(x) for x in c1["el_edaph"]], [bool(x) for x in c1["el_clim"]]), ([True, False, False, True, False], [True, False, True, True, False]))
        self.assertEqual(([bool(x) for x in c1["el_open"]], [bool(x) for x in c1["el_qtype"]]), ([False, False, False, True, False], [False, True, True, False, False]))
        self.assertEqual([bool(x) for x in c1["fc_edaph"]], [True, False, False, False])
        self.assertTrue(bool(c1["roof_blocked"]) and str(c1["src_hyb"]) == "Passive_Open_or_Hybrid_Concept")
        self.assertEqual([str(x) for x in c1["src_qtypes"]], ["Passive_Type_Uncertain_Coord"] * 3)
        self.assertTrue(bool(c1["clean_types_active"]) and float(c1["rf_oof_diff"]) < 1e-12)
        self.assertTrue(bool(c1["sdm_npres_eq"]) and float(c1["sdm_auc_diff"]) < 1e-12)
        self.assertTrue(bool(c1["hop_n_eq"]) and float(c1["hop_slope_diff"]) < 1e-12 and float(c1["hop_ddoy_diff"]) < 1e-12)
        self.assertTrue(bool(c1["ct_tot_eq"]) and bool(c1["ct_out_eq"]) and bool(c1["tok_var_ok"]) and bool(c1["tok_ssp_ok"]) and bool(c1["tok_hyb_match"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_02_basionym_anchored_locus_classicus_and_heterotypic_audit(self) -> None:
        c2 = self.results["c2"]
        self.assertTrue(float(c2["hom_z_shift"]) < 0.15 and float(c2["hom_holo_z"]) < 0.50 and float(c2["het_vs_hom_z"]) > 2.50)
        self.assertEqual((str(c2["tc_het_st"]), str(c2["tc_het_hint"]), str(c2["het_rf_src"])), ("Type_Synonymy_Candidate", "Passive_Heterotypic_Type_Audit", "Passive_Heterotypic_Type_Audit"))
        self.assertTrue(bool(c2["het_audit_flag"]) and bool(c2["het_conf_flag"]))
        self.assertEqual((str(c2["bas_hom_a"]), str(c2["bas_hom_d"]), str(c2["bas_het"])), ("Senecio aureus", "Cineraria dubia", "Senecio alabamensis"))
        self.assertTrue(_is_finite(c2["bas_z_het"]) and float(c2["bas_z_het"]) < 1.0 and _is_finite(c2["bas_z_hom"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_03_pedology_layer_hygiene_smooth_indices_point_mass_and_cache(self) -> None:
        c3 = self.results["c3"]
        self.assertEqual((float(c3["rion_pct"]), float(c3["lou_pct"])), (0.0, 0.0))
        self.assertAlmostEqual(float(c3["wake_out_pct"]), 25.0, places=2); self.assertAlmostEqual(float(c3["wake_lith_pct"]), 65.0, places=2)
        self.assertTrue(bool(c3["excl_ro"]) and bool(c3["excl_aws"]) and bool(c3["keep_ph"]))
        self.assertAlmostEqual(float(c3["hp_sol_bd"]), 175.0, places=1); self.assertAlmostEqual(float(c3["hp_buf_bd"]), 175.0, places=1)
        self.assertTrue(float(c3["loi_v_max"]) <= 1.0 and float(c3["loi_r_max"]) <= 1.0 and float(c3["loi_formula_err"]) < 1e-3 and float(c3["loi_vr_err"]) < 0.04 and float(c3["wsi_formula_err"]) < 1e-3)
        self.assertTrue(bool(c3["pm_dropped_bd"]) and (c3["pan_grps"] if isinstance(c3["pan_grps"], list) else [c3["pan_grps"]]) == ["Packera aurea"])
        self.assertTrue((int(c3["can_calls_1"]), int(c3["can_calls_2"])) == (1, 0) and bool(c3["can_zero_refresh"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_04_fold_safe_sdm_haversine_ocean_exclusion_and_circular_phenology(self) -> None:
        c4 = self.results["c4"]
        self.assertTrue(int(c4["bg_n"]) >= 20 and bool(c4["bg_all_land"]) and float(c4["bg_max_km"]) <= 260.0 and _is_finite(c4["sdm_auc_ens"]) and 0.0 <= float(c4["sdm_auc_ens"]) < 0.78)
        self.assertTrue(bool(c4["has_adw"]) and bool(c4["has_hresid"]) and bool(c4["pk_ver"]) and bool(c4["hresid_fin"]) and int(c4["anth_base_n"]) == 12)
        self.assertTrue([bool(x) for x in c4["fr_ver"]] == [False, False] and all(float(w) < 0.15 for w in c4["fr_adw"]) and float(c4["pk_adw_min"]) >= 0.80)
        self.assertTrue(_is_finite(c4["bm_dbic"]) and float(c4["bm_dbic"]) > 0.0 and float(c4["bm_ashman"]) > 2.0 and _is_finite(c4["um_dbic"]) and float(c4["um_dbic"]) <= 2.0)

    def test_05_sqlite_roundtrip_python_and_r_step05_clarity_columns(self) -> None:
        cols = ("ecotype_regime", "lithic_outcrop_index", "hopkins_adjusted_doy_resid", "nearest_type_niche_basionym", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic", "nearest_type_basionym_edaphic_z")
        with sqlite3.connect(str(self.db_path)) as conn:
            uv = int(conn.execute("PRAGMA user_version;").fetchone()[0])
            self.assertEqual((uv, uv), (PLANTMATICS_SCHEMA_USER_VERSION, 20261127))
            self.assertEqual(persist_spatial_rf_predictions(conn, RUN_PY, [{
                "catalogNumber": "PY_RT_01", "plant_individual_id": 0, "latitude": 35.80, "longitude": -79.10, "ecotype_regime": "Ecotype_1_Lithic_High",
                "lithic_outcrop_index": 0.7425, "summer_moisture_deficit_idx": 1.85, "maritime_continentality_idx": 2.42, "hopkins_adjusted_doy": 109.5,
                "hopkins_adjusted_doy_resid": -6.25, "nearest_type_niche_basionym": "Cineraria dubia", "nearest_type_basionym_edaphic_z": 0.415,
                "anthesis_dominance_weight": 0.875, "allochronic_bimodality_delta_bic": 14.32, "rf_prediction_source": "Spatial_Block_OOF_Anchor",
                "rf_predicted_taxon": "Packera dubia", "rf_predicted_prob": 0.93,
            }], demo=False), 1)
            q = f"SELECT {', '.join(cols)} FROM spatial_rf_predictions WHERE run_id=? AND catalogNumber=?;"
            if not SKIP_R_COND:
                self.assertTrue(bool(self.results["c4"]["ok_r_sql"]))
            for rid, cat in [(RUN_PY, "PY_RT_01")] + ([] if SKIP_R_COND else [(RUN_R, "R_RT_01")]):
                row = conn.execute(q, (rid, cat)).fetchone()
                self.assertIsNotNone(row); self.assertEqual((row[0], row[3]), ("Ecotype_1_Lithic_High", "Cineraria dubia"))
                for got, exp in zip((row[1], row[2], row[4], row[5], row[6]), (0.7425, -6.25, 0.875, 14.32, 0.415)):
                    self.assertAlmostEqual(float(got), exp, places=4)

    def test_06_static_step05_clarity_exports_and_posix_hygiene(self) -> None:
        self.assertLess(len(Path(__file__).resolve().read_text(encoding="utf-8").splitlines()), 380)
        ana = PROJECT_ROOT / "scripts" / "analysis"
        for fname, fn in (("edaphic_collinearity_contingency_utils.R", "screen_edaphic_point_mass_columns"),
                          ("step05_layer_hygiene_utils.R", "is_fine_edaphic_type_or_anchor_eligible"),
                          ("step05_allochronic_pheno_utils.R", "test_allochronic_bimodality_von_mises")):
            p = ana / fname
            self.assertTrue(p.is_file() and bool(re.search(rf"^{re.escape(fn)}\s*<-\s*function\(", p.read_text(encoding="utf-8"), flags=re.MULTILINE)))
        for fname in ("step05_layer_hygiene_utils.R", "edaphic_geochem_wetland_utils.R", "edaphic_ssurgo_sda_utils.R", "edaphic_solus_3dep_utils.R",
                      "edaphic_collinearity_contingency_utils.R", "step05_allochronic_pheno_utils.R", "step05_type_niche_utils.R", "step05_sdm_engine_utils.R"):
            if (ana / fname).is_file():
                self.assertIsNone(re.search(r"\\\\[sdwb]", "\n".join(ln for ln in (ana / fname).read_text(encoding="utf-8").splitlines() if not ln.lstrip().startswith("#"))), f"Non-POSIX regex escape in {fname}")


if __name__ == "__main__":
    unittest.main()
