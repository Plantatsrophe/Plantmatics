#!/usr/bin/env python3
"""Warm-R tests for Step 07 queue coalescing, multi-plant ramet consensus, ANCOVA, and non-dubia precedence."""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    ensure_pipeline_run,
    get_db_connection,
    persist_otu_morphometrics,
    persist_spatial_rf_predictions,
    populate_run_active_cohort,
    upsert_vouchers_records,
)
from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _seed_sqlite_cohort(db_path: Path, run_id: str) -> None:
    """Seed minimal Tier A + Tier B rows to verify SQLite queue coalescing in Step 07."""
    with contextlib.closing(get_db_connection(db_path=db_path, init_schema=True)) as conn:
        v_rows = [
            {
                "catalogNumber": "SQL_DUB_01", "institutionCode": "NCU", "scientificName": "Packera dubia",
                "species_raw": "Packera dubia", "species_standardized": "Packera dubia",
                "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.8, "decimalLongitude": -79.0,
                "doy": 115, "stateProvince": "NC", "county": "Orange", "usda_county_status": "in_range",
                "typified_basionym": "Senecio tomentosus", "concept_synonymy_status": "Homotypic_Synonym",
                "is_open_or_hybrid_concept": 0, "qc_status": "PASS",
            },
            {
                "catalogNumber": "SQL_AUR_01", "institutionCode": "DUKE", "scientificName": "Packera aurea",
                "species_raw": "Packera aurea", "species_standardized": "Packera aurea",
                "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.9, "decimalLongitude": -79.1,
                "doy": 122, "stateProvince": "NC", "county": "Orange", "usda_county_status": "in_range",
                "typified_basionym": "Senecio aureus", "concept_synonymy_status": "Homotypic_Synonym",
                "is_open_or_hybrid_concept": 0, "qc_status": "PASS",
            },
        ]
        upsert_vouchers_records(conn, v_rows)
        with conn:
            for r in v_rows:
                conn.execute("INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 1);", (r["catalogNumber"],))
        ensure_pipeline_run(conn, run_id=run_id)
        populate_run_active_cohort(conn, run_id=run_id, voucher_rows=v_rows, clear_existing=True)
        persist_otu_morphometrics(
            conn, run_id,
            [
                {"catalogNumber": "SQL_DUB_01", "plant_individual_id": 1, "species_standardized": "Packera dubia",
                 "PC1": -2.0, "PC2": 0.1, "can1": -2.0, "can2": 0.1, "morph_tier_source": "Basal_Rosette",
                 "basal_leaf_count": 3, "cda_predicted_taxon": "Packera dubia", "cda_posterior_prob": 0.94,
                 "open_set_novelty_p": 0.62, "gmm_cluster": 1},
                {"catalogNumber": "SQL_AUR_01", "plant_individual_id": 1, "species_standardized": "Packera aurea",
                 "PC1": 2.0, "PC2": -0.1, "can1": 2.0, "can2": -0.1, "morph_tier_source": "Basal_Rosette",
                 "basal_leaf_count": 3, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.95,
                 "open_set_novelty_p": 0.71, "gmm_cluster": 1},
            ],
        )
        persist_spatial_rf_predictions(
            conn, run_id,
            [
                {"catalogNumber": "SQL_DUB_01", "plant_individual_id": 1, "rf_predicted_taxon": "Packera dubia",
                 "rf_predicted_prob": 0.91, "edaphic_best_fit_taxon": "Packera dubia", "edaphic_z_score": 0.4,
                 "hopkins_adjusted_doy": 115.0, "hopkins_adjusted_doy_resid": 1.2, "anthesis_dominance_weight": 0.88,
                 "allochronic_bimodality_delta_bic": 0.4, "lithic_outcrop_index": 0.72,
                 "nearest_type_niche_basionym": "Senecio tomentosus"},
                {"catalogNumber": "SQL_AUR_01", "plant_individual_id": 1, "rf_predicted_taxon": "Packera aurea",
                 "rf_predicted_prob": 0.93, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_z_score": 0.3,
                 "hopkins_adjusted_doy": 122.0, "hopkins_adjusted_doy_resid": -0.8, "anthesis_dominance_weight": 0.92,
                 "allochronic_bimodality_delta_bic": 0.3, "lithic_outcrop_index": 0.14,
                 "nearest_type_niche_basionym": "Senecio aureus"},
            ],
        )


_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; work_dir <- args[2L]; db_path <- args[3L]; sql_run_id <- args[4L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/analysis/07_triage_dashboard_synthesis.R")

# 1. CLI alias & CSV queue coalescing (historic ungeoreferenced, multi-plant, open_set_p_value alias, new cols, focal full export)
t_dir <- file.path(work_dir, "csv_case", "data", "tables")
r_dir <- file.path(work_dir, "csv_case", "outputs", "reports")
dir.create(t_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_dir, recursive = TRUE, showWarnings = FALSE)

v_csv <- data.frame(
  catalogNumber = c("DUB_01", "DUB_02", "DUB_03", "AUR_01", "AUR_02", "AUR_03", "NCU-HIST-2", "CAT-MULTI-SHEET", "CAT-MIXED-SHEET", "NCU-SOLO-LEAFLESS"),
  institutionCode = c("NCU", "DUKE", "US", "NCU", "DUKE", "US", "NCU", "NCU", "NCU", "NCU"),
  scientificName = c(rep("Packera dubia", 3L), rep("Packera aurea", 3L), "Packera dubia", "Packera dubia", "Packera dubia", "Packera dubia"),
  species_standardized = c(rep("Packera dubia", 3L), rep("Packera aurea", 3L), "Packera dubia", "Packera dubia", "Packera dubia", "Packera dubia"),
  determiner_tier = c(rep("Tier_1_Gold", 6L), "Tier_2_Silver", "Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"),
  decimalLatitude = c(35.80, 35.81, 35.82, 35.85, 35.86, 35.87, NA_real_, 35.83, 35.84, 35.82),
  decimalLongitude = c(-79.00, -79.01, -79.02, -79.10, -79.11, -79.12, NA_real_, -79.03, -79.04, -79.02),
  stateProvince = c(rep("NC", 6L), NA_character_, "NC", "NC", "NC"),
  county = c(rep("Orange", 6L), NA_character_, "Orange", "Orange", "Orange"),
  usda_county_status = rep("in_range", 10L),
  doy = c(114, 116, 118, 121, 123, 125, 115, 116, 117, 116),
  typified_basionym = c(rep("Senecio tomentosus", 3L), rep("Senecio aureus", 3L), "Senecio tomentosus", "Senecio tomentosus", "Senecio smallii", "Senecio tomentosus"),
  concept_synonymy_status = c(rep("Homotypic_Synonym", 8L), "Heterotypic_Synonym", "Homotypic_Synonym"),
  anchor_tier_source = c(rep("Tier_1_Gold_Core", 6L), "Tier_2_Silver_Regional", "Tier_1_Gold_Core", "Tier_2_Silver_Regional", "Tier_3_Bronze"),
  is_open_or_hybrid_concept = c(rep(FALSE, 9L), FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(v_csv, file.path(t_dir, "curated_vouchers.csv"), row.names = FALSE)

m_csv <- data.frame(
  catalogNumber = c("DUB_01", "DUB_02", "DUB_03", "AUR_01", "AUR_02", "AUR_03", "NCU-HIST-2", "CAT-MULTI-SHEET", "CAT-MULTI-SHEET", "CAT-MIXED-SHEET", "CAT-MIXED-SHEET"),
  plant_individual_id = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L),
  morph_tier_source = c(rep("Basal_Rosette", 8L), "Cauline_Only_Fallback", "Basal_Rosette", "Basal_Rosette"),
  basal_leaf_count = c(rep(3L, 8L), 0L, 4L, 3L),
  cauline_leaf_count = c(rep(1L, 8L), 1L, 1L, 1L),
  PC1 = c(-2.0, -2.1, -1.9, 2.0, 2.1, 1.9, -2.05, -1.95, 1.80, -2.0, 2.2),
  PC2 = c(0.1, -0.1, 0.0, 0.1, -0.1, 0.0, 0.05, 0.0, 0.3, 0.0, 0.1),
  can1 = c(-2.0, -2.1, -1.9, 2.0, 2.1, 1.9, -2.05, -1.95, 1.80, -2.0, 2.2),
  can2 = c(0.1, -0.1, 0.0, 0.1, -0.1, 0.0, 0.05, 0.0, 0.3, 0.0, 0.1),
  cda_predicted_taxon = c(rep("Packera dubia", 3L), rep("Packera aurea", 3L), "Packera dubia", "Packera dubia", "Packera aurea", "Packera dubia", "Packera aurea"),
  cda_posterior_prob = c(0.95, 0.94, 0.96, 0.95, 0.96, 0.94, 0.91, 0.95, 0.88, 0.95, 0.96),
  cda_posterior_margin = c(0.90, 0.88, 0.92, 0.90, 0.92, 0.88, 0.82, 0.90, 0.76, 0.90, 0.92),
  type_morph_mahal_d2 = c(1.2, 1.4, 1.1, 1.3, 1.2, 1.5, 1.6, 1.2, 19.5, 1.3, 30.0),
  self_morph_mahal_d2 = c(1.2, 1.4, 1.1, 1.3, 1.2, 1.5, 1.6, 1.2, 19.5, 1.3, 30.0),
  is_typical_morphotype = c(rep(TRUE, 8L), FALSE, TRUE, FALSE),
  morphometric_outlier_flag = c(rep(FALSE, 8L), TRUE, FALSE, TRUE),
  open_set_p_value = c(0.75, 0.72, 0.78, 0.80, 0.74, 0.77, 0.42, 0.69, 0.15, 0.70, 0.004),
  gmm_cluster = rep(1L, 11L),
  stringsAsFactors = FALSE
)
utils::write.csv(m_csv, file.path(t_dir, "morphometric_flags.csv"), row.names = FALSE)

mm_csv <- data.frame(
  catalogNumber = c("DUB_01", "DUB_02", "DUB_03", "AUR_01", "AUR_02", "AUR_03", "CAT-MULTI-SHEET", "CAT-MIXED-SHEET", "NCU-SOLO-LEAFLESS"),
  rf_predicted_taxon = c(rep("Packera dubia", 3L), rep("Packera aurea", 3L), "Packera dubia", "Packera aurea", "Packera dubia"),
  rf_predicted_prob = c(0.92, 0.91, 0.93, 0.92, 0.94, 0.90, 0.91, 0.89, 0.88),
  edaphic_best_fit_taxon = c(rep("Packera dubia", 3L), rep("Packera aurea", 3L), "Packera dubia", "Packera aurea", "Packera dubia"),
  edaphic_z_score = c(0.3, 0.4, 0.2, 0.3, 0.2, 0.5, 0.4, 3.4, 0.5),
  hopkins_adjusted_doy = c(114, 116, 118, 121, 123, 125, 116, 117, 116),
  hopkins_adjusted_doy_resid = c(-1.5, 0.2, 1.8, -1.2, 0.5, 1.1, 0.4, 2.2, 0.0),
  anthesis_dominance_weight = c(0.90, 0.85, 0.92, 0.88, 0.95, 0.87, 0.82, 0.86, 0.78),
  allochronic_bimodality_delta_bic = c(0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.5, 3.7, 0.5),
  lithic_outcrop_index = c(0.74, 0.71, 0.69, 0.12, 0.15, 0.11, 0.73, 0.18, 0.70),
  wetland_seepage_index = c(0.10, 0.12, 0.08, 0.82, 0.85, 0.79, 0.11, 0.78, 0.09),
  nearest_type_niche_basionym = c(rep("Senecio tomentosus", 3L), rep("Senecio aureus", 3L), "Senecio tomentosus", "Senecio aureus", "Senecio tomentosus"),
  stringsAsFactors = FALSE
)
utils::write.csv(mm_csv, file.path(t_dir, "multimodal_conflict_flags.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(catalogNumber = v_csv$catalogNumber, voucher_phenological_state = rep("anthesis", nrow(v_csv)), n_anthesis = rep(5L, nrow(v_csv)), stringsAsFactors = FALSE),
  file.path(t_dir, "voucher_phenological_states.csv"), row.names = FALSE
)

opts_csv <- parse_args_robust(c(
  "--vouchers", file.path(t_dir, "curated_vouchers.csv"),
  "--morphometrics", file.path(t_dir, "morphometric_flags.csv"),
  "--spatial-flags", file.path(t_dir, "multimodal_conflict_flags.csv"),
  "--output-queue", file.path(t_dir, "triage_queue.csv"),
  "--output-report", file.path(r_dir, "revision.md"),
  "--output-html", file.path(r_dir, "dashboard.html"),
  "--focal-taxon", "Packera dubia",
  "--demo"
))
q_csv <- build_triage_queue(opts_csv)
q_full <- utils::read.csv(file.path(t_dir, "triage_queue_complex_full.csv"), stringsAsFactors = FALSE)

# 2. SQLite queue coalescing path
sql_tdir <- file.path(work_dir, "sql_case", "data", "tables")
sql_rdir <- file.path(work_dir, "sql_case", "outputs", "reports")
dir.create(sql_tdir, recursive = TRUE, showWarnings = FALSE)
dir.create(sql_rdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(catalogNumber = c("SQL_DUB_01", "SQL_AUR_01"), species_standardized = c("Packera dubia", "Packera aurea"), stringsAsFactors = FALSE),
  file.path(sql_tdir, "curated_vouchers.csv"), row.names = FALSE
)
utils::write.csv(data.frame(catalogNumber = "SQL_DUB_01"), file.path(sql_tdir, "morphometric_flags.csv"), row.names = FALSE)
utils::write.csv(data.frame(catalogNumber = "SQL_DUB_01"), file.path(sql_tdir, "multimodal_conflict_flags.csv"), row.names = FALSE)
old_db <- Sys.getenv("PLANTMATICS_DB_PATH", unset = NA_character_)
old_pdb <- Sys.getenv("PACKERA_DB_PATH", unset = NA_character_)
old_rid <- Sys.getenv("RUN_ID", unset = NA_character_)
on.exit({
  if (is.na(old_db)) Sys.unsetenv("PLANTMATICS_DB_PATH") else Sys.setenv(PLANTMATICS_DB_PATH = old_db)
  if (is.na(old_pdb)) Sys.unsetenv("PACKERA_DB_PATH") else Sys.setenv(PACKERA_DB_PATH = old_pdb)
  if (is.na(old_rid)) Sys.unsetenv("RUN_ID") else Sys.setenv(RUN_ID = old_rid)
}, add = TRUE)
Sys.setenv(PLANTMATICS_DB_PATH = db_path, PACKERA_DB_PATH = db_path, RUN_ID = sql_run_id)
opts_sql <- parse_args_robust(c(
  "--vouchers", file.path(sql_tdir, "curated_vouchers.csv"),
  "--morphometrics", file.path(sql_tdir, "morphometric_flags.csv"),
  "--multimodal-flags", file.path(sql_tdir, "multimodal_conflict_flags.csv"),
  "--output-queue", file.path(sql_tdir, "triage_queue.csv"),
  "--output-report", file.path(sql_rdir, "revision.md")
))
q_sql <- build_triage_queue(opts_sql)

# 3. Marginal drop1 ANCOVA (clean vs severe herbarium batch shift)
set.seed(19L)
mk_anc <- function(us_shift) {
  n_a <- 24L
  sp <- rep(c("Packera aurea", "Packera obovata"), each = 12L)
  inst <- rep(c("NCU", "DUKE", "US"), times = 8L)
  eps_a <- rep(c(-0.20, -0.20, -0.20, 0.20, 0.20, 0.20, -0.10, -0.10, -0.10, 0.10, 0.10, 0.10), times = 2L)
  data.frame(
    catalogNumber = sprintf("ANC_%02d", seq_len(n_a)),
    species_standardized = sp,
    institutionCode = inst,
    regional_group = rep(c("Piedmont", "Mountains"), length.out = n_a),
    morph_tier_source = rep("Basal_Rosette", n_a),
    selected_morph_features = rep("PC1;PC2;blade_cie_l_norm", n_a),
    PC1 = ifelse(sp == "Packera aurea", -2.0, 2.0) + ifelse(inst == "US", us_shift, 0.0) + eps_a,
    PC2 = eps_a * 0.8,
    blade_cie_l_norm = ifelse(inst == "US", 90.0, 50.0),
    authority_weight = rep(2.0, n_a),
    spatial_balance_weight = rep(1.0, n_a),
    stringsAsFactors = FALSE
  )
}
anc_clean <- audit_institutional_batch_effects_ancova(mk_anc(0.0), min_obs = 4L)
anc_shift <- audit_institutional_batch_effects_ancova(mk_anc(5.5), min_obs = 4L)

# 4. Non-dubia taxon-agnostic 5-stage precedence & self-reassignment guard
nd_df <- data.frame(
  catalogNumber = c(paste0("AUR_", 1:5), paste0("OBO_", 1:5), "ND_MISID", "ND_ECOTYPE", "ND_PLASTIC", "ND_SELF_GUARD"),
  species_standardized = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera aurea", "Packera obovata", "Packera aurea", "Packera aurea"),
  species_raw = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera aurea", "Packera obovata var. montana", "Packera aurea", "Packera aurea"),
  determiner_tier = c(rep("Tier_1_Gold", 10L), "Tier_2_Silver", "Tier_2_Silver", "Tier_3_Bronze", "Tier_3_Bronze"),
  decimalLatitude = c(rep(35.8, 5L), rep(35.82, 5L), 35.81, 35.83, 35.80, 35.80),
  decimalLongitude = c(rep(-79.0, 5L), rep(-79.02, 5L), -79.01, -79.03, -79.00, -79.00),
  stateProvince = rep("NC", 14L), county = rep("Orange", 14L), usda_county_status = rep("in_range", 14L),
  phenological_state = rep("anthesis", 14L),
  doy = c(rep(115, 5L), rep(120, 5L), 120, 145, 138, 115),
  PC1 = c(rep(-2.0, 5L), rep(2.0, 5L), 2.1, 3.2, -1.3, -2.0),
  PC2 = c(rep(0.0, 10L), 0.1, 1.1, 0.2, 0.0),
  can1 = c(rep(-2.0, 5L), rep(2.0, 5L), 2.1, 3.2, -1.3, -2.0),
  can2 = c(rep(0.0, 10L), 0.1, 1.1, 0.2, 0.0),
  cda_predicted_taxon = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera obovata", "Packera aurea", "Packera aurea"),
  cda_posterior_prob = c(rep(0.95, 10L), 0.94, 0.62, 0.58, 0.90),
  edaphic_best_fit_taxon = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera obovata", "Packera aurea", "Packera aurea"),
  rf_predicted_taxon = c(rep("Packera aurea", 5L), rep("Packera obovata", 5L), "Packera obovata", "Packera obovata", "Packera aurea", "Packera aurea"),
  rf_predicted_prob = c(rep(0.92, 10L), 0.91, 0.60, 0.85, 0.90),
  edaphic_z_score = c(rep(0.4, 10L), 3.4, 2.4, 0.8, 0.4),
  self_morph_mahal_d2 = c(rep(1.0, 10L), 28.0, 8.0, 12.0, 1.2),
  gmm_cluster = c(rep(1L, 11L), 2L, 1L, 1L),
  morphometric_outlier_flag = c(rep(FALSE, 10L), TRUE, FALSE, TRUE, FALSE),
  morph_tier_source = rep("Basal_Rosette", 14L),
  basal_leaf_count = rep(3L, 14L),
  stringsAsFactors = FALSE
)
res_nd <- apply_taxonomic_decision_matrix(nd_df, intraspecific_mode = FALSE, opts = list(demo = TRUE))
pick_sub <- function(d, cid, pid = NULL) {
  m <- if (is.null(pid)) which(d$catalogNumber == cid) else which(d$catalogNumber == cid & d$plant_individual_id == pid)
  if (length(m) == 0L) return(list(sub = NA_character_, rec = NA_character_, rat = NA_character_))
  list(sub = as.character(d$taxonomic_treatment_subtype[m[1L]]), rec = as.character(d$recommended_determination[m[1L]]), rat = as.character(d$synthesis_rationale[m[1L]]))
}

i_hist <- which(q_csv$catalogNumber == "NCU-HIST-2")[1L]
i_mix2 <- which(q_full$catalogNumber == "CAT-MIXED-SHEET" & q_full$plant_individual_id == 2L)[1L]
write_json(list(
  spatial_alias = opts_csv$multimodal_flags,
  hist_present = !is.na(i_hist),
  hist_score = if (!is.na(i_hist)) as.numeric(q_csv$probabilistic_consensus_score[i_hist]) else NA_real_,
  hist_osr_p = if (!is.na(i_hist) && "open_set_novelty_p" %in% names(q_csv)) as.numeric(q_csv$open_set_novelty_p[i_hist]) else NA_real_,
  multi_pids = sort(as.integer(q_csv$plant_individual_id[q_csv$catalogNumber == "CAT-MULTI-SHEET"])),
  multi_r1 = pick_sub(q_csv, "CAT-MULTI-SHEET", 1L),
  multi_r2 = pick_sub(q_csv, "CAT-MULTI-SHEET", 2L),
  mixed_r1 = pick_sub(q_full, "CAT-MIXED-SHEET", 1L),
  mixed_r2 = pick_sub(q_full, "CAT-MIXED-SHEET", 2L),
  solo_leafless = pick_sub(q_csv, "NCU-SOLO-LEAFLESS"),
  coalesced_cols_present = all(c(
    "hopkins_adjusted_doy_resid", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic",
    "lithic_outcrop_index", "wetland_seepage_index", "nearest_type_niche_basionym",
    "typified_basionym", "concept_synonymy_status", "anchor_tier_source", "is_open_or_hybrid_concept"
  ) %in% names(q_full)),
  mix2_vals = list(
    hop_resid = as.numeric(q_full$hopkins_adjusted_doy_resid[i_mix2]),
    anth_wt = as.numeric(q_full$anthesis_dominance_weight[i_mix2]),
    bimod_bic = as.numeric(q_full$allochronic_bimodality_delta_bic[i_mix2]),
    lithic = as.numeric(q_full$lithic_outcrop_index[i_mix2]),
    wetland = as.numeric(q_full$wetland_seepage_index[i_mix2]),
    nt_bas = as.character(q_full$nearest_type_niche_basionym[i_mix2]),
    typ_bas = as.character(q_full$typified_basionym[i_mix2]),
    syn_st = as.character(q_full$concept_synonymy_status[i_mix2]),
    anc_src = as.character(q_full$anchor_tier_source[i_mix2])
  ),
  focal_taxa_in_filtered = I(sort(unique(as.character(q_csv$species_standardized)))),
  full_taxa_in_complex = I(sort(unique(as.character(q_full$species_standardized)))),
  sql_cats = sort(as.character(q_sql$catalogNumber)),
  sql_lithic_dub = as.numeric(q_sql$lithic_outcrop_index[q_sql$catalogNumber == "SQL_DUB_01"][1L]),
  sql_typ_bas_dub = as.character(q_sql$typified_basionym[q_sql$catalogNumber == "SQL_DUB_01"][1L]),
  anc_clean_status = as.character(anc_clean$Batch_Effect_Status[anc_clean$Response == "PC1"][1L]),
  anc_clean_eta = as.numeric(anc_clean$partial_eta_sq[anc_clean$Response == "PC1"][1L]),
  anc_shift_status = as.character(anc_shift$Batch_Effect_Status[anc_shift$Response == "PC1"][1L]),
  anc_shift_eta = as.numeric(anc_shift$partial_eta_sq[anc_shift$Response == "PC1"][1L]),
  anc_shift_qbh = as.numeric(anc_shift$p_adj_bh[anc_shift$Response == "PC1"][1L]),
  anc_responses = as.character(anc_shift$Response),
  nd_misid = pick_sub(res_nd, "ND_MISID"),
  nd_ecotype = pick_sub(res_nd, "ND_ECOTYPE"),
  nd_plastic = pick_sub(res_nd, "ND_PLASTIC"),
  nd_core = pick_sub(res_nd, "AUR_1"),
  nd_any_dubia = any(grepl("dubia", c(res_nd$recommended_determination, res_nd$synthesis_rationale), ignore.case = TRUE)),
  nd_any_self_reassign = any(res_nd$taxonomic_treatment_subtype == "Misidentification_Reassign" & res_nd$recommended_determination == res_nd$species_standardized)
), out_json, auto_unbox = TRUE, na = "null")
"""


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R queue/ANCOVA suite under fast audit or missing Rscript.")
class TestStep07QueueAndAncova(unittest.TestCase):
    """Contract tests for Step 07 CSV/SQLite queue coalescing, multi-plant sheets, ANCOVA, and non-dubia precedence."""

    res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        with tempfile.TemporaryDirectory(prefix="step07_q_anc_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            db_path = tmp / "step07_queue_test.sqlite"
            run_id = "step07_queue_sql_unit"
            _seed_sqlite_cohort(db_path, run_id)
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(
                ["Rscript", str(r_p), str(out_j), str(tmp), str(db_path), run_id],
                stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True,
            )
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))

    def test_historic_ungeoref_multi_plant_ramets_and_open_set_alias(self) -> None:
        r = self.res
        self.assertTrue(r["spatial_alias"].endswith("multimodal_conflict_flags.csv"))
        self.assertTrue(r["hist_present"])
        self.assertIsNotNone(r["hist_score"])
        self.assertGreater(float(r["hist_score"]), 0.5)
        self.assertAlmostEqual(float(r["hist_osr_p"]), 0.42, places=3)
        self.assertEqual(r["multi_pids"], [1, 2])
        self.assertEqual(r["multi_r1"]["sub"], "Confirmed_Species")
        self.assertEqual(r["multi_r2"]["sub"], "Confirmed_Species")
        self.assertEqual(r["mixed_r1"]["sub"], "Confirmed_Species")
        self.assertEqual(r["mixed_r2"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["mixed_r2"]["rec"], "Packera aurea")
        self.assertNotIn(r["solo_leafless"]["sub"], ("Hybrid_Introgressant", "Putative_Cryptic_Species"))

    def test_new_step01_step05_coalesced_columns_focal_full_and_sqlite_hydration(self) -> None:
        r = self.res
        self.assertTrue(r["coalesced_cols_present"])
        mv = r["mix2_vals"]
        self.assertAlmostEqual(float(mv["hop_resid"]), 2.2, places=2)
        self.assertAlmostEqual(float(mv["anth_wt"]), 0.86, places=2)
        self.assertAlmostEqual(float(mv["bimod_bic"]), 3.7, places=2)
        self.assertAlmostEqual(float(mv["lithic"]), 0.18, places=2)
        self.assertAlmostEqual(float(mv["wetland"]), 0.78, places=2)
        self.assertEqual(mv["nt_bas"], "Senecio aureus")
        self.assertEqual(mv["typ_bas"], "Senecio smallii")
        self.assertEqual(mv["syn_st"], "Heterotypic_Synonym")
        self.assertEqual(mv["anc_src"], "Tier_2_Silver_Regional")
        self.assertEqual(r["focal_taxa_in_filtered"], ["Packera dubia"])
        self.assertEqual(r["full_taxa_in_complex"], ["Packera aurea", "Packera dubia"])
        self.assertEqual(r["sql_cats"], ["SQL_AUR_01", "SQL_DUB_01"])
        self.assertAlmostEqual(float(r["sql_lithic_dub"]), 0.72, places=2)
        self.assertEqual(r["sql_typ_bas_dub"], "Senecio tomentosus")

    def test_marginal_drop1_ancova_fdr_and_non_dubia_five_stage_precedence(self) -> None:
        r = self.res
        self.assertEqual(r["anc_clean_status"], "Pass (<0.05)")
        self.assertLess(float(r["anc_clean_eta"]), 0.05)
        self.assertEqual(r["anc_shift_status"], "Exceeds 0.05")
        self.assertGreater(float(r["anc_shift_eta"]), 0.50)
        self.assertLess(float(r["anc_shift_qbh"]), 0.05)
        self.assertNotIn("blade_cie_l_norm", r["anc_responses"])
        self.assertEqual(r["nd_misid"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["nd_misid"]["rec"], "Packera obovata")
        self.assertEqual(r["nd_ecotype"]["sub"], "Subspecies_Ecotype")
        self.assertEqual(r["nd_plastic"]["sub"], "Ecophenotypic_Plasticity")
        self.assertEqual(r["nd_core"]["sub"], "Confirmed_Species")
        self.assertFalse(bool(r["nd_any_dubia"]))
        self.assertFalse(bool(r["nd_any_self_reassign"]))


if __name__ == "__main__":
    unittest.main()
