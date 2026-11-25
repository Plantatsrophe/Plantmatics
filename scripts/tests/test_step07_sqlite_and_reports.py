#!/usr/bin/env python3
"""Warm-R & Python tests for Step 07 SQLite round-trip, offline HTML/MD reports, fig07a ECT fallback, and 3-module focal CLI."""

from __future__ import annotations

import contextlib
import csv
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    ensure_pipeline_run,
    get_db_connection,
    ingest_existing_run_tables_into_db,
    persist_otu_morphometrics,
    persist_spatial_rf_predictions,
    populate_run_active_cohort,
    upsert_vouchers_records,
)
from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION
from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

NEW_EIGHT_TRIAGE_COLS = (
    "w_morph_simplex",
    "w_edaphic_simplex",
    "w_pheno_simplex",
    "w_geno_simplex",
    "geno_deme_weight_scale",
    "portal_metadata_lag_flag",
    "portal_lag_trigger_gate",
    "is_marginal_range_extension",
)


def _seed_roundtrip_db(db_path: Path, py_run: str, r_run: str, ect_run: str, empty_run: str, tables_dir: Path) -> None:
    """Seed Tier A + Tier B cohorts and ingest a Python Step 07 triage CSV containing all 8 newly persisted columns."""
    v_rows = [
        {"catalogNumber": "RT_A_COALESCE", "institutionCode": "NCU", "scientificName": "Packera dubia",
         "species_raw": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold",
         "decimalLatitude": 35.80, "decimalLongitude": -79.00, "stateProvince": "NC", "county": "Orange",
         "usda_county_status": "in_range", "qc_status": "PASS"},
        {"catalogNumber": "RT_B_OVERRIDE", "institutionCode": "DUKE", "scientificName": "Packera aurea",
         "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver",
         "decimalLatitude": 35.85, "decimalLongitude": -79.05, "stateProvince": "NC", "county": "Orange",
         "usda_county_status": "out_of_range", "qc_status": "PASS"},
        {"catalogNumber": "RT_ECT_ONLY", "institutionCode": "NCU", "scientificName": "Packera millefolium",
         "species_raw": "Packera millefolium", "species_standardized": "Packera millefolium", "determiner_tier": "Tier_1_Gold",
         "decimalLatitude": 35.20, "decimalLongitude": -82.50, "stateProvince": "NC", "county": "Transylvania",
         "usda_county_status": "in_range", "qc_status": "PASS"},
    ]
    with contextlib.closing(get_db_connection(db_path=db_path, init_schema=True)) as conn:
        upsert_vouchers_records(conn, v_rows)
        with conn:
            for r in v_rows:
                conn.execute("INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 1);", (r["catalogNumber"],))

        for rid in (py_run, r_run):
            ensure_pipeline_run(conn, run_id=rid)
            populate_run_active_cohort(conn, run_id=rid, voucher_rows=v_rows[:2], clear_existing=True)
            persist_otu_morphometrics(conn, rid, [
                {"catalogNumber": "RT_A_COALESCE", "plant_individual_id": 1, "species_standardized": "Packera dubia",
                 "PC1": -2.0, "PC2": 0.1, "can1": -2.0, "can2": 0.1, "cda_predicted_taxon": "Packera dubia", "cda_posterior_prob": 0.94},
                {"catalogNumber": "RT_B_OVERRIDE", "plant_individual_id": 1, "species_standardized": "Packera aurea",
                 "PC1": 2.0, "PC2": -0.1, "can1": 2.0, "can2": -0.1, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.91},
            ])
            persist_spatial_rf_predictions(conn, rid, [
                {"catalogNumber": "RT_A_COALESCE", "plant_individual_id": 1, "rf_predicted_taxon": "Packera dubia",
                 "rf_predicted_prob": 0.90, "edaphic_best_fit_taxon": "Packera dubia", "edaphic_z_score": 0.4,
                 "hopkins_adjusted_doy": 115.0, "hopkins_adjusted_doy_resid": 14.25, "anthesis_dominance_weight": 0.90},
                {"catalogNumber": "RT_B_OVERRIDE", "plant_individual_id": 1, "rf_predicted_taxon": "Packera aurea",
                 "rf_predicted_prob": 0.88, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_z_score": 0.5,
                 "hopkins_adjusted_doy": 122.0, "hopkins_adjusted_doy_resid": 14.25, "anthesis_dominance_weight": 0.85},
            ])

        # Multi-island ECT_PC1-only cohort (PC1, can1, cda_predicted_taxon all NULL)
        ensure_pipeline_run(conn, run_id=ect_run)
        populate_run_active_cohort(conn, run_id=ect_run, voucher_rows=[v_rows[2]], clear_existing=True)
        persist_otu_morphometrics(conn, ect_run, [
            {"catalogNumber": "RT_ECT_ONLY", "plant_individual_id": 1, "species_standardized": "Packera millefolium",
             "ECT_PC1": 2.75, "ECT_PC2": -0.85, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
        ])
        ensure_pipeline_run(conn, run_id=empty_run)
        populate_run_active_cohort(conn, run_id=empty_run, voucher_rows=[v_rows[2]], clear_existing=True)
        persist_otu_morphometrics(conn, empty_run, [
            {"catalogNumber": "RT_ECT_ONLY", "plant_individual_id": 1, "species_standardized": "Packera millefolium",
             "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
        ])

        tables_dir.mkdir(parents=True, exist_ok=True)
        tri_csv = tables_dir / "triage_queue.csv"
        fieldnames = [
            "catalogNumber", "plant_individual_id", "taxonomic_status_call", "recommended_determination",
            "taxonomic_treatment_broad", "taxonomic_treatment_subtype", "synthesis_triage_priority",
            "synthesis_triage_action", "synthesis_rationale", "delta_doy_resid", "w_morph_kappa", "w_edaphic_kappa",
            "w_pheno_kappa", "w_geno_kappa", *NEW_EIGHT_TRIAGE_COLS,
            "resurrection_basionym", "genomic_propagation_tier", "cytology_linkage_tier",
            "hydro_edaphic_niche_regime", "basionym_synonymy_type",
        ]
        rows_py = [
            {"catalogNumber": "RT_A_COALESCE", "plant_individual_id": 1, "taxonomic_status_call": "Species_Confirmed",
             "recommended_determination": "Packera dubia", "taxonomic_treatment_broad": "Confirmed Species",
             "taxonomic_treatment_subtype": "Confirmed_Species", "synthesis_triage_priority": "RESOLVED",
             "synthesis_triage_action": "Accept_Current_Determination", "synthesis_rationale": "Verified anchor",
             "delta_doy_resid": "", "w_morph_kappa": 0.85, "w_edaphic_kappa": 0.72, "w_pheno_kappa": 0.60, "w_geno_kappa": 0.90,
             "w_morph_simplex": 0.30, "w_edaphic_simplex": 0.25, "w_pheno_simplex": 0.15, "w_geno_simplex": 0.30,
             "geno_deme_weight_scale": 0.75, "portal_metadata_lag_flag": "FALSE", "portal_lag_trigger_gate": "",
             "is_marginal_range_extension": "TRUE", "resurrection_basionym": "",
             "genomic_propagation_tier": "G2_Conditional_Deme_Hull", "cytology_linkage_tier": "C1_Direct_Cytovoucher",
             "hydro_edaphic_niche_regime": "Siliceous_Felsic_Outcrop", "basionym_synonymy_type": "Homotypic"},
            {"catalogNumber": "RT_B_OVERRIDE", "plant_individual_id": 1, "taxonomic_status_call": "Misidentification_Severe",
             "recommended_determination": "Packera dubia", "taxonomic_treatment_broad": "Misidentification / Cryptic Species",
             "taxonomic_treatment_subtype": "Misidentification_Reassign", "synthesis_triage_priority": "HIGH",
             "synthesis_triage_action": "Reassign_Determination", "synthesis_rationale": "Gate3 lag",
             "delta_doy_resid": -9.50, "w_morph_kappa": 0.85, "w_edaphic_kappa": 0.72, "w_pheno_kappa": 0.60, "w_geno_kappa": "",
             "w_morph_simplex": 0.45, "w_edaphic_simplex": 0.35, "w_pheno_simplex": 0.20, "w_geno_simplex": 0.0,
             "geno_deme_weight_scale": "", "portal_metadata_lag_flag": "TRUE", "portal_lag_trigger_gate": "Gate3_Extralimital_Bioclimatic",
             "is_marginal_range_extension": "FALSE", "resurrection_basionym": "Senecio smallii",
             "genomic_propagation_tier": "", "cytology_linkage_tier": "C3_Taxon_Monomorphic",
             "hydro_edaphic_niche_regime": "Mesic_Upland_Matrix", "basionym_synonymy_type": "Heterotypic"},
        ]
        with open(tri_csv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows_py)
        ingest_existing_run_tables_into_db(conn=conn, tables_dir=tables_dir, run_id=py_run, demo=False)


_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; work_dir <- args[2L]; db_path <- args[3L]
py_run <- args[4L]; r_run <- args[5L]; ect_run <- args[6L]; empty_run <- args[7L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/morphometrics/pub_figures_utils.R")
source("scripts/analysis/07_triage_dashboard_synthesis.R")
source("scripts/figures/step07_figures.R")

old_db <- Sys.getenv("PLANTMATICS_DB_PATH", unset = NA_character_)
old_pdb <- Sys.getenv("PACKERA_DB_PATH", unset = NA_character_)
old_rid <- Sys.getenv("RUN_ID", unset = NA_character_)
on.exit({
  if (is.na(old_db)) Sys.unsetenv("PLANTMATICS_DB_PATH") else Sys.setenv(PLANTMATICS_DB_PATH = old_db)
  if (is.na(old_pdb)) Sys.unsetenv("PACKERA_DB_PATH") else Sys.setenv(PACKERA_DB_PATH = old_pdb)
  if (is.na(old_rid)) Sys.unsetenv("RUN_ID") else Sys.setenv(RUN_ID = old_rid)
}, add = TRUE)
Sys.setenv(PLANTMATICS_DB_PATH = db_path, PACKERA_DB_PATH = db_path)

# 1. R persist_step07_to_sqlite() & try_fetch_master_view_for_cohort() round-trip + ECT-only readiness
df_r_tri <- data.frame(
  catalogNumber = c("RT_A_COALESCE", "RT_B_OVERRIDE"),
  plant_individual_id = c(1L, 1L),
  taxonomic_status_call = c("Species_Confirmed", "Type_Synonymy_Candidate"),
  recommended_determination = c("Packera dubia", "Packera aurea"),
  taxonomic_treatment_broad = c("Confirmed Species", "Type Synonymy / Resurrection"),
  taxonomic_treatment_subtype = c("Confirmed_Species", "Type_Synonymy_Candidate"),
  synthesis_triage_priority = c("RESOLVED", "CRITICAL"),
  synthesis_triage_action = c("Accept_Current_Determination", "Evaluate_Heterotypic_Synonymy"),
  synthesis_rationale = c("R persist verified", "Heterotypic type divergence"),
  delta_doy_resid = c(NA_real_, -6.25),
  w_morph_simplex = c(0.32, 0.50), w_edaphic_simplex = c(0.28, 0.30),
  w_pheno_simplex = c(0.15, 0.20), w_geno_simplex = c(0.25, 0.00),
  geno_deme_weight_scale = c(0.66, NA_real_),
  portal_metadata_lag_flag = c(FALSE, TRUE),
  portal_lag_trigger_gate = c(NA_character_, "Gate1_Exsiccatae_Conflict"),
  is_marginal_range_extension = c(TRUE, FALSE),
  resurrection_basionym = c(NA_character_, "Senecio smallii"),
  stringsAsFactors = FALSE
)
ok_persist <- isTRUE(persist_step07_to_sqlite(df_r_tri, vouchers_path = NULL, run_id = r_run, demo = FALSE))
mv_r <- try_fetch_master_view_for_cohort(run_id = r_run, demo = FALSE, require_step04 = TRUE, require_step05 = TRUE, require_step07 = TRUE)
mv_py <- try_fetch_master_view_for_cohort(run_id = py_run, demo = FALSE, require_step04 = TRUE, require_step05 = TRUE, require_step07 = TRUE)
mv_ect <- try_fetch_master_view_for_cohort(run_id = ect_run, demo = FALSE, require_step04 = TRUE)
mv_empty <- try_fetch_master_view_for_cohort(run_id = empty_run, demo = FALSE, require_step04 = TRUE)

# 2. Offline HTML dashboard & Markdown revision report (with build_genomic_markdown_summary)
rep_dir <- file.path(work_dir, "reports")
fig_dir <- file.path(work_dir, "figures")
gen_dir <- file.path(work_dir, "genomics")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gen_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(data.frame(
  genetic_sample_id = c("GS_1", "GS_2", "GS_3"),
  catalogNumber = c("REP_1", "REP_2", NA_character_),
  species_standardized = c("Packera dubia", "Packera aurea", "Packera dubia"),
  linkage_tier = c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G3_Georeferenced_Unvouchered"),
  geno_quality_pass = c(TRUE, TRUE, TRUE),
  geno_predicted_taxon = c("Packera dubia", "Packera aurea", "Packera dubia"),
  geno_posterior_prob = c(0.96, 0.94, 0.91),
  stringsAsFactors = FALSE
), file.path(gen_dir, "genomic_assignments.csv"), row.names = FALSE)

df_rep <- data.frame(
  catalogNumber = c("REP_1", "REP_2", "REP_TYPE_SYN", "REP_ECT_1", "REP_ECT_2", "REP_ECT_3"),
  plant_individual_id = rep(1L, 6L),
  institutionCode = c("NCU", "DUKE", "NY", "NCU", "US", "DUKE"),
  recordedBy = paste0("Coll", 1:6),
  scientificName = c("Packera dubia", "Packera aurea", "Senecio smallii", "Packera millefolium", "Packera millefolium", "Packera dubia"),
  species_standardized = c("Packera dubia", "Packera aurea", "Packera dubia", "Packera millefolium", "Packera millefolium", "Packera dubia"),
  determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver", "Tier_2_Silver"),
  is_type_specimen = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
  decimalLatitude = c(35.80, 35.90, 35.85, 35.20, 35.25, 35.75),
  decimalLongitude = c(-79.00, -79.10, -79.05, -82.50, -82.55, -78.95),
  regional_group = rep("Piedmont", 6L),
  doy = c(115, 122, 118, 130, 132, 116),
  is_active_anthesis = rep(TRUE, 6L),
  delta_doy_resid = c(0.5, -1.0, 2.0, 1.5, -0.5, 0.8),
  can1 = rep(NA_real_, 6L), can2 = rep(NA_real_, 6L),
  PC1 = rep(NA_real_, 6L), PC2 = rep(NA_real_, 6L),
  ECT_PC1 = c(-2.1, 2.2, 1.9, 4.5, 4.7, -1.8),
  ECT_PC2 = c(0.2, -0.3, 0.1, 1.4, 1.6, -0.1),
  cda_predicted_taxon = c("Packera dubia", "Packera aurea", "Packera aurea", "Packera millefolium", "Packera millefolium", "Packera dubia"),
  cda_posterior_prob = c(0.94, 0.95, 0.91, 0.96, 0.92, 0.90),
  edaphic_best_fit_taxon = c("Packera dubia", "Packera aurea", "Packera aurea", "Packera millefolium", "Packera millefolium", "Packera dubia"),
  rf_predicted_taxon = c("Packera dubia", "Packera aurea", "Packera aurea", "Packera millefolium", "Packera millefolium", "Packera dubia"),
  edaphic_z_score = c(0.4, 0.3, 3.1, 0.5, 0.6, 0.4),
  fine_soil_ph = c(5.5, 6.4, 6.3, 6.8, 6.7, 5.4),
  fine_soil_sand = c(48, 35, 36, 62, 60, 50),
  fine_bedrock_depth_cm = c(40, 110, 95, 20, 25, 45),
  w_morph_kappa = rep(0.84, 6L), w_edaphic_kappa = rep(0.71, 6L), w_pheno_kappa = rep(0.58, 6L),
  probabilistic_consensus_score = c(0.93, 0.94, 0.31, 0.95, 0.91, 0.90),
  taxonomic_status_call = c("Species_Confirmed", "Species_Confirmed", "Type_Synonymy_Candidate", "Species_Confirmed", "Species_Confirmed", "Species_Confirmed"),
  taxonomic_treatment_broad = c("Confirmed Species", "Confirmed Species", "Type Synonymy / Resurrection", "Confirmed Species", "Confirmed Species", "Confirmed Species"),
  taxonomic_treatment_subtype = c("Confirmed_Species", "Confirmed_Species", "Type_Synonymy_Candidate", "Confirmed_Species", "Confirmed_Species", "Confirmed_Species"),
  recommended_determination = c("Packera dubia", "Packera aurea", "Packera dubia", "Packera millefolium", "Packera millefolium", "Packera dubia"),
  synthesis_triage_priority = c("RESOLVED", "RESOLVED", "CRITICAL", "RESOLVED", "RESOLVED", "RESOLVED"),
  synthesis_triage_action = c("Accept_Current_Determination", "Accept_Current_Determination", "Evaluate_Heterotypic_Synonymy", "Accept_Current_Determination", "Accept_Current_Determination", "Accept_Current_Determination"),
  synthesis_rationale = c("Verified", "Verified", "Heterotypic type audit", "Dissected ECT anchor", "Dissected ECT anchor", "Verified"),
  morph_tier_source = rep("Basal_Rosette", 6L),
  basal_leaf_count = rep(3L, 6L), cauline_leaf_count = rep(1L, 6L),
  stringsAsFactors = FALSE
)

seed_plates <- c(
  "fig01_cohort_biogeography_and_authority.png", "fig02a_lm3_organ_and_ect_exemplars.png",
  "fig02b_segmentation_qc_and_phenophase.png", "fig03a_efa_morphospace_and_harmonics.png",
  "fig03b_topological_ect_and_heteroblasty.png", "fig04a_foliar_allometry_audit.png",
  "fig04b_gmm_cda_morphometrics.png", "fig04c_ecoregion_and_batch_audit.png",
  "fig05a_ecological_niche_and_rf_importance.png", "fig05b_hopkins_circular_phenology.png",
  "fig05c_sdm_suitability_and_distribution_map.png", "fig06a_genomic_pcoa_tree_concordance.png",
  "fig06b_reticulation_gdi_landscape.png", "fig06c_cytonuclear_paralog_cytotype_audit.png",
  "fig07a_integrative_taxonomic_synthesis.png", "fig07b_taxon_complex_delimitation_network.png"
)
for (sp_fn in seed_plates) writeBin(as.raw(c(137L, 80L, 78L, 71L)), file.path(fig_dir, sp_fn))

html_out <- file.path(rep_dir, "triage_dashboard.html")
md_out <- file.path(rep_dir, "revision_report.md")
png_out <- file.path(fig_dir, "fig07a_ect_fallback.png")
opts_rep <- list(output_report = md_out, output_queue = file.path(gen_dir, "triage_queue.csv"),
                 output_plot = png_out, genomic_dir = gen_dir, focal_taxon = "", gmm_bic = "", gmm_summary = "")
generate_interactive_triage_html(df_rep, html_out, focal_taxon = "", opts = opts_rep)
generate_taxonomic_revision_report(df_rep, opts_rep)
render_step07a_integrative_synthesis_png(df_rep, out_plot_path = png_out)

pick_mv <- function(d, cid) {
  k <- which(d$catalogNumber == cid)[1L]
  list(
    delta_doy_resid = as.numeric(d$delta_doy_resid[k]),
    w_morph_simplex = as.numeric(d$w_morph_simplex[k]),
    w_edaphic_simplex = as.numeric(d$w_edaphic_simplex[k]),
    w_pheno_simplex = as.numeric(d$w_pheno_simplex[k]),
    w_geno_simplex = as.numeric(d$w_geno_simplex[k]),
    geno_deme_weight_scale = if (is.na(d$geno_deme_weight_scale[k])) NULL else as.numeric(d$geno_deme_weight_scale[k]),
    portal_metadata_lag_flag = as.integer(d$portal_metadata_lag_flag[k]),
    portal_lag_trigger_gate = if (is.na(d$portal_lag_trigger_gate[k])) NULL else as.character(d$portal_lag_trigger_gate[k]),
    is_marginal_range_extension = as.integer(d$is_marginal_range_extension[k])
  )
}

pdf_files <- list.files(work_dir, pattern = "\\.pdf$", recursive = TRUE, ignore.case = TRUE)
write_json(list(
  r_schema_ver = as.integer(PLANTMATICS_SCHEMA_USER_VERSION),
  ok_persist = ok_persist,
  mv_r_a = pick_mv(mv_r, "RT_A_COALESCE"),
  mv_r_b = pick_mv(mv_r, "RT_B_OVERRIDE"),
  mv_py_a = pick_mv(mv_py, "RT_A_COALESCE"),
  mv_py_b = pick_mv(mv_py, "RT_B_OVERRIDE"),
  mv_ect_nonnull = !is.null(mv_ect) && nrow(mv_ect) == 1L,
  mv_ect_pc1 = if (!is.null(mv_ect) && nrow(mv_ect) == 1L) as.numeric(mv_ect$ECT_PC1[1L]) else NA_real_,
  mv_empty_isnull = is.null(mv_empty),
  html_path = html_out, md_path = md_out, png_path = png_out,
  n_pdfs = length(pdf_files)
), out_json, auto_unbox = TRUE, na = "null")
"""


class TestStep07SqliteAndReports(unittest.TestCase):
    """Contract tests for Step 07 SQLite round-trip, offline HTML/MD reports, fig07a ECT fallback, and focal CLI exports."""

    res: Dict[str, Any] = {}
    py_db_rows: Dict[str, Dict[str, Any]] = {}
    py_user_version: int = 0
    html_text: str = ""
    md_text: str = ""
    png_bytes: bytes = b""

    @classmethod
    def setUpClass(cls) -> None:
        with tempfile.TemporaryDirectory(prefix="step07_sql_rep_") as td:
            tmp = Path(td)
            db_path = tmp / "step07_roundtrip.sqlite"
            py_run, r_run, ect_run, empty_run = "rt_py_run", "rt_r_run", "rt_ect_only", "rt_empty_morph"
            _seed_roundtrip_db(db_path, py_run, r_run, ect_run, empty_run, tmp / "py_tables")
            with contextlib.closing(get_db_connection(db_path=db_path, init_schema=False)) as conn:
                cls.py_user_version = int(conn.execute("PRAGMA user_version;").fetchone()[0])
                cols_sql = ", ".join(["catalogNumber", "delta_doy_resid", *NEW_EIGHT_TRIAGE_COLS])
                cur = conn.execute(f"SELECT {cols_sql} FROM v_integrative_triage_master WHERE run_id = ?;", (py_run,))
                cls.py_db_rows = {str(r["catalogNumber"]): dict(r) for r in cur.fetchall()}

            if SKIP_R_COND:
                return
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(
                ["Rscript", str(r_p), str(out_j), str(tmp), str(db_path), py_run, r_run, ect_run, empty_run],
                stage_key="step07", cwd=PROJECT_ROOT, skip_figures=False,
            )
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))
            cls.html_text = Path(cls.res["html_path"]).read_text(encoding="utf-8")
            cls.md_text = Path(cls.res["md_path"]).read_text(encoding="utf-8")
            cls.png_bytes = Path(cls.res["png_path"]).read_bytes()

    def test_python_sqlite_schema_version_and_eight_triage_columns_roundtrip(self) -> None:
        self.assertEqual(PLANTMATICS_SCHEMA_USER_VERSION, 20261127)
        self.assertEqual(self.py_user_version, 20261127)
        ra, rb = self.py_db_rows["RT_A_COALESCE"], self.py_db_rows["RT_B_OVERRIDE"]
        self.assertAlmostEqual(float(ra["delta_doy_resid"]), 14.25, places=3)
        self.assertAlmostEqual(float(rb["delta_doy_resid"]), -9.50, places=3)
        self.assertAlmostEqual(sum(float(ra[k]) for k in NEW_EIGHT_TRIAGE_COLS[:4]), 1.0, places=5)
        self.assertAlmostEqual(float(ra["geno_deme_weight_scale"]), 0.75, places=3)
        self.assertEqual(int(ra["portal_metadata_lag_flag"]), 0)
        self.assertEqual(int(ra["is_marginal_range_extension"]), 1)
        self.assertEqual(int(rb["portal_metadata_lag_flag"]), 1)
        self.assertEqual(rb["portal_lag_trigger_gate"], "Gate3_Extralimital_Bioclimatic")
        self.assertEqual(int(rb["is_marginal_range_extension"]), 0)

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_r_sqlite_roundtrip_coalesce_and_multi_island_ect_readiness(self) -> None:
        r = self.res
        self.assertEqual(int(r["r_schema_ver"]), 20261127)
        self.assertTrue(bool(r["ok_persist"]))
        for tag in ("mv_r_a", "mv_py_a"):
            self.assertAlmostEqual(float(r[tag]["delta_doy_resid"]), 14.25, places=3)
            self.assertEqual(int(r[tag]["is_marginal_range_extension"]), 1)
            self.assertEqual(int(r[tag]["portal_metadata_lag_flag"]), 0)
        self.assertAlmostEqual(float(r["mv_r_b"]["delta_doy_resid"]), -6.25, places=3)
        self.assertEqual(r["mv_r_b"]["portal_lag_trigger_gate"], "Gate1_Exsiccatae_Conflict")
        self.assertAlmostEqual(float(r["mv_r_a"]["geno_deme_weight_scale"]), 0.66, places=3)
        self.assertAlmostEqual(float(r["mv_py_b"]["delta_doy_resid"]), -9.50, places=3)
        self.assertTrue(bool(r["mv_ect_nonnull"]))
        self.assertAlmostEqual(float(r["mv_ect_pc1"]), 2.75, places=3)
        self.assertTrue(bool(r["mv_empty_isnull"]))

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_offline_html_markdown_genomic_summary_and_fig07a_ect_fallback(self) -> None:
        self.assertIn("Export Filtered CSV", self.html_text)
        self.assertIn('<option value="Type_Synonymy_Candidate">', self.html_text)
        self.assertIn("REP_TYPE_SYN", self.html_text)
        self.assertIn("_lm3_sheet_overlay.png", self.html_text)
        self.assertIn("replace(/[^A-Za-z0-9_.-]+/g, '_')", self.html_text)
        for ext in ("'.webp'", "'.jpg'", "'.jpeg'", "'.png'"):
            self.assertIn(ext, self.html_text)
        self.assertIn(".s-Confirmed_Type_Control{background:#EEF2FF;color:#332288;border:1px solid #C7D2FE;}", self.html_text)
        self.assertIn("sub-badge s-Confirmed_Type_Control", self.html_text)
        for plate in (
            "fig01_cohort_biogeography_and_authority.png", "fig02a_lm3_organ_and_ect_exemplars.png",
            "fig02b_segmentation_qc_and_phenophase.png", "fig03a_efa_morphospace_and_harmonics.png",
            "fig03b_topological_ect_and_heteroblasty.png", "fig04a_foliar_allometry_audit.png",
            "fig04b_gmm_cda_morphometrics.png", "fig04c_ecoregion_and_batch_audit.png",
            "fig05a_ecological_niche_and_rf_importance.png", "fig05b_hopkins_circular_phenology.png",
            "fig05c_sdm_suitability_and_distribution_map.png", "fig06a_genomic_pcoa_tree_concordance.png",
            "fig06b_reticulation_gdi_landscape.png", "fig06c_cytonuclear_paralog_cytotype_audit.png",
            "fig07a_integrative_taxonomic_synthesis.png", "fig07b_taxon_complex_delimitation_network.png",
        ):
            self.assertIn(f'href="../figures/{plate}"', self.html_text)
        self.assertIsNone(re.search(r"https?://", self.html_text), "Offline HTML dashboard must have zero external CDN URLs")
        self.assertNotIn("\u2013", self.md_text)
        self.assertNotIn("\u2014", self.md_text)
        self.assertIn("Rule A (Heterotypic Type Synonymy Candidates)", self.md_text)
        self.assertIn("Independent Step 06 Genomic & Phylogenomic Evidence Stream", self.md_text)
        self.assertIn("2 voucher-linked (`G1`/`G2`)", self.md_text)
        self.assertTrue(self.png_bytes.startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertGreater(len(self.png_bytes), 10000)
        self.assertEqual(int(self.res["n_pdfs"]), 0)

    def test_resolve_focal_species_complex_three_module_reexports(self) -> None:
        from scripts.pipeline import resolve_focal_species_complex as rfsc
        from scripts.pipeline import step07_focal_cohort_utils as cohort_mod
        from scripts.pipeline import step07_focal_stage_runners as runner_mod

        expected: List[str] = [
            "slugify_taxon", "default_prefix_for_taxon", "_stage_run_ids", "build_four_stage_commands",
            "_read_csv_rows", "_fetch_triage_rows_for_run", "_fetch_complex_rows_for_run",
            "_extract_stage_metrics", "compile_complex_resolution_summary", "build_parser", "main",
        ]
        for sym in expected:
            self.assertTrue(callable(getattr(rfsc, sym, None)), f"Missing callable re-export: {sym}")
        self.assertIs(rfsc.slugify_taxon, cohort_mod.slugify_taxon)
        self.assertIs(rfsc.build_four_stage_commands, runner_mod.build_four_stage_commands)
        self.assertIs(rfsc.compile_complex_resolution_summary, runner_mod.compile_complex_resolution_summary)


if __name__ == "__main__":
    unittest.main()
