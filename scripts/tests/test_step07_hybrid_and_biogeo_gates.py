#!/usr/bin/env python3
"""Warm-R tests for Step 07 5-gate hybrid swarms, Buck Creek sympatry, marginal range extensions, and Portal Lag Gates 1/2/3/5."""

from __future__ import annotations

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

from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; work_dir <- args[2L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/analysis/07_triage_dashboard_synthesis.R")

# 1. Conservative 5-gate hybrid evaluation (4 isolated syntopic sites >50 km apart) + Buck Creek pure-parent guard + intraspecific suppression
# Site A (35.10, -83.50): deme >= 2 + pure parent + allochronic + off-chord; Site B (35.60): singleton x-formula; Site C (36.10): singleton G1 admixed; Site D (36.60): singleton reject; Site E (41.50): allopatric reject
n_p <- 8L
p_lats <- c(35.10, 35.12, 35.60, 35.61, 36.10, 36.11, 36.60, 36.61)
p_lons <- rep(-83.50, n_p)
p_cnty <- c("Clay", "Clay", "Macon", "Macon", "Buncombe", "Buncombe", "Watauga", "Watauga")
df_hyb <- data.frame(
  catalogNumber = c(
    sprintf("DUB_%02d", 1:n_p), sprintf("AUR_%02d", 1:n_p),
    "HYB_DEME_1", "HYB_DEME_2", "HYB_X_SINGLE", "HYB_G1_ADM",
    "HYB_FAIL_ALLOPATRIC", "HYB_FAIL_ALLOCHRONIC", "HYB_FAIL_OFFCHORD", "HYB_FAIL_SINGLETON",
    "BUCK_CREEK_PURE"
  ),
  species_standardized = c(rep("Packera dubia", n_p), rep("Packera aurea", n_p), rep("Packera dubia", 9L)),
  species_raw = c(
    rep("Packera dubia", n_p), rep("Packera aurea", n_p),
    "Packera dubia", "Packera dubia", "Packera dubia x aurea", "Packera dubia",
    "Packera dubia", "Packera dubia", "Packera dubia", "Packera dubia", "Packera dubia"
  ),
  determiner_tier = c(rep("Tier_1_Gold", 2L * n_p), rep("Tier_2_Silver", 8L), "Tier_1_Gold"),
  decimalLatitude = c(
    p_lats, p_lats + 0.01,
    35.11, 35.12, 35.605, 36.105,
    41.50, 35.11, 35.11, 36.605, 35.11
  ),
  decimalLongitude = c(
    p_lons, p_lons - 0.01,
    -83.505, -83.508, -83.505, -83.505,
    -75.00, -83.505, -83.505, -83.505, -83.502
  ),
  stateProvince = c(rep("NC", 2L * n_p + 4L), "PA", rep("NC", 4L)),
  county = c(p_cnty, p_cnty, "Clay", "Clay", "Macon", "Buncombe", "Berks", "Clay", "Clay", "Watauga", "Clay"),
  usda_county_status = rep("in_range", 2L * n_p + 9L),
  recordedBy = paste0("Coll_", seq_len(2L * n_p + 9L)), recordNumber = as.character(101:(100 + 2L * n_p + 9L)),
  institutionCode = rep("NCU", 2L * n_p + 9L),
  phenological_state = rep("anthesis", 2L * n_p + 9L),
  doy = c(rep(118, n_p), rep(122, n_p), 120, 121, 120, 120, 120, 185, 120, 120, 118),
  PC1 = c(
    seq(-3.1, -2.9, length.out = n_p), seq(2.9, 3.1, length.out = n_p),
    0.0, 0.1, -0.05, 0.05,
    0.0, 0.0, 0.0, 0.0, -3.0
  ),
  PC2 = c(
    seq(-0.1, 0.1, length.out = n_p), seq(-0.1, 0.1, length.out = n_p),
    0.02, -0.02, 0.0, 0.01,
    0.0, 0.0, 18.0, 0.0, 0.02
  ),
  can1 = c(
    seq(-3.1, -2.9, length.out = n_p), seq(2.9, 3.1, length.out = n_p),
    0.0, 0.1, -0.05, 0.05,
    0.0, 0.0, 0.0, 0.0, -3.0
  ),
  can2 = c(
    seq(-0.1, 0.1, length.out = n_p), seq(-0.1, 0.1, length.out = n_p),
    0.02, -0.02, 0.0, 0.01,
    0.0, 0.0, 18.0, 0.0, 0.02
  ),
  cda_predicted_taxon = c(rep("Packera dubia", n_p), rep("Packera aurea", n_p), rep("Packera aurea", 8L), "Packera dubia"),
  cda_posterior_prob = c(rep(0.95, 2L * n_p), rep(0.51, 8L), 0.94),
  cda_second_taxon = c(rep("Packera aurea", n_p), rep("Packera dubia", n_p), rep("Packera dubia", 8L), "Packera aurea"),
  cda_second_posterior_prob = c(rep(0.05, 2L * n_p), rep(0.49, 8L), 0.06),
  cda_posterior_margin = c(rep(0.90, 2L * n_p), rep(0.02, 8L), 0.88),
  gmm_uncertainty = c(rep(0.02, 2L * n_p), rep(0.36, 8L), 0.02),
  gmm_cluster = rep(1L, 2L * n_p + 9L),
  edaphic_best_fit_taxon = c(rep("Packera dubia", n_p), rep("Packera aurea", n_p), rep("Packera dubia", 9L)),
  rf_predicted_taxon = c(rep("Packera dubia", n_p), rep("Packera aurea", n_p), rep("Packera dubia", 9L)),
  rf_predicted_prob = c(rep(0.90, 2L * n_p), rep(0.55, 8L), 0.92),
  edaphic_z_score = c(rep(0.3, 2L * n_p), rep(1.1, 8L), 0.2),
  geno_linkage_tier = c(rep(NA_character_, 2L * n_p + 3L), "G1_Direct_Voucher", rep(NA_character_, 5L)),
  geno_quality_pass = c(rep(NA, 2L * n_p + 3L), TRUE, rep(NA, 5L)),
  geno_predicted_taxon = c(rep(NA_character_, 2L * n_p + 3L), "Packera dubia", rep(NA_character_, 5L)),
  geno_posterior_prob = c(rep(NA_real_, 2L * n_p + 3L), 0.52, rep(NA_real_, 5L)),
  is_geno_admixed = c(rep(FALSE, 2L * n_p + 3L), TRUE, rep(FALSE, 5L)),
  interspecific_ancestry_q12 = c(rep(NA_real_, 2L * n_p + 3L), 0.72, rep(NA_real_, 5L)),
  hybrid_generation_class = c(rep(NA_character_, 2L * n_p + 3L), "F1_Hybrid", rep(NA_character_, 5L)),
  geno_parent_1 = c(rep(NA_character_, 2L * n_p + 3L), "Packera dubia", rep(NA_character_, 5L)),
  geno_parent_2 = c(rep(NA_character_, 2L * n_p + 3L), "Packera aurea", rep(NA_character_, 5L)),
  morph_tier_source = rep("Basal_Rosette", 2L * n_p + 9L),
  basal_leaf_count = rep(3L, 2L * n_p + 9L),
  stringsAsFactors = FALSE
)
res_hyb <- apply_taxonomic_decision_matrix(df_hyb, intraspecific_mode = FALSE, opts = list(demo = TRUE))
res_intra_sup <- apply_taxonomic_decision_matrix(df_hyb, intraspecific_mode = TRUE, opts = list(demo = TRUE))
res_single_sup <- apply_taxonomic_decision_matrix(df_hyb[df_hyb$species_standardized == "Packera dubia", , drop = FALSE], intraspecific_mode = FALSE, opts = list(demo = TRUE))

# 2. Marginal range extension + Portal Metadata Lag Gates 1, 2, 3a, 3b & Gate 5 audit CSV export
p_tdir <- file.path(work_dir, "portal_case", "data", "tables")
p_rdir <- file.path(work_dir, "portal_case", "outputs", "reports")
dir.create(p_tdir, recursive = TRUE, showWarnings = FALSE)
dir.create(p_rdir, recursive = TRUE, showWarnings = FALSE)

df_port <- data.frame(
  catalogNumber = c(
    sprintf("PA_NC_%02d", 1:5), sprintf("PD_NC_%02d", 1:5),
    "PORTAL_G1_EXS", "PORTAL_G2_BARKLEY", "PORTAL_G3_BIOCLIM", "PORTAL_G3_UNGEOREF", "PORTAL_MARG_EXT", "PORTAL_WAIF_EXOTIC"
  ),
  species_standardized = c(rep("Packera aurea", 5L), rep("Packera dubia", 5L), rep("Packera dubia", 6L)),
  determiner_tier = c(rep("Tier_1_Gold", 10L), "Tier_2_Silver", "Tier_1_Gold", "Tier_2_Silver", "Tier_2_Silver", "Tier_2_Silver", "Tier_2_Silver"),
  identifiedBy = c(rep("A. Weakley", 10L), "Unk Det", "T. M. Barkley", "Old Det", "Old Det", "Field Bot", "Old Det"),
  institutionCode = rep("NCU", 16L),
  recordedBy = paste0("CollP", 1:16), recordNumber = as.character(201:216),
  stateProvince = c(rep("NC", 5L), rep("NC", 5L), "NC", "NC", "NC", "NC", "VA", "NC"),
  county = c(rep("Orange", 5L), rep("Wake", 5L), "Orange", "Wake", "Orange", "Orange", "Mecklenburg", "Orange"),
  decimalLatitude = c(rep(35.90, 5L), rep(35.78, 5L), 35.90, 35.78, 35.90, NA_real_, 36.68, 35.90),
  decimalLongitude = c(rep(-79.05, 5L), rep(-78.64, 5L), -79.05, -78.64, -79.05, NA_real_, -78.35, -79.05),
  usda_county_status = c(rep("in_range", 10L), "in_range", "in_range", "out_of_range", "out_of_range", "out_of_range", "natureserve_exotic_or_waif"),
  exsiccatae_taxonomic_conflict = c(rep(0L, 10L), 1L, 0L, 0L, 0L, 0L, 0L),
  exsiccatae_gold_taxon = c(rep(NA_character_, 10L), "Packera aurea", "Packera dubia", NA_character_, NA_character_, NA_character_, NA_character_),
  climatic_z_score = c(rep(0.4, 10L), 0.8, 0.3, 4.9, NA_real_, 1.4, 5.2),
  rf_given_taxon_prob = c(rep(0.92, 10L), 0.15, 0.94, 0.05, NA_real_, 0.62, 0.08),
  rf_predicted_taxon = c(rep("Packera aurea", 5L), rep("Packera dubia", 5L), "Packera aurea", "Packera dubia", "Packera aurea", NA_character_, "Packera dubia", "Packera aurea"),
  rf_predicted_prob = c(rep(0.92, 10L), 0.88, 0.94, 0.91, NA_real_, 0.62, 0.90),
  edaphic_best_fit_taxon = c(rep("Packera aurea", 5L), rep("Packera dubia", 5L), "Packera aurea", "Packera dubia", "Packera aurea", NA_character_, "Packera dubia", "Packera aurea"),
  edaphic_z_score = c(rep(0.3, 10L), 2.5, 0.2, 3.8, NA_real_, 1.1, 3.9),
  phenological_state = rep("anthesis", 16L),
  doy = c(rep(115, 5L), rep(120, 5L), 116, 120, 116, 116, 121, 116),
  PC1 = c(rep(-2.0, 5L), rep(2.0, 5L), -1.8, 2.0, -1.9, -2.0, 1.9, -1.9),
  PC2 = rep(0.0, 16L),
  can1 = c(rep(-2.0, 5L), rep(2.0, 5L), -1.8, 2.0, -1.9, -2.0, 1.9, -1.9),
  can2 = rep(0.0, 16L),
  cda_predicted_taxon = c(rep("Packera aurea", 5L), rep("Packera dubia", 5L), "Packera aurea", "Packera dubia", "Packera aurea", "Packera aurea", "Packera dubia", "Packera aurea"),
  cda_posterior_prob = c(rep(0.95, 10L), 0.88, 0.96, 0.86, 0.90, 0.92, 0.89),
  self_morph_mahal_d2 = c(rep(1.0, 10L), 22.0, 0.8, 24.0, 25.0, 1.2, 23.0),
  gmm_cluster = rep(1L, 16L),
  morph_tier_source = rep("Basal_Rosette", 16L),
  basal_leaf_count = rep(3L, 16L),
  stringsAsFactors = FALSE
)
utils::write.csv(df_port, file.path(p_tdir, "curated_vouchers.csv"), row.names = FALSE)
opts_p <- parse_args_robust(c(
  "--vouchers", file.path(p_tdir, "curated_vouchers.csv"),
  "--output-queue", file.path(p_tdir, "triage_queue.csv"),
  "--output-report", file.path(p_rdir, "revision.md"),
  "--demo"
))
res_port <- build_triage_queue(opts_p)
audit_csv_p <- file.path(p_rdir, "portal_metadata_lag_audit.csv")
audit_df <- if (file.exists(audit_csv_p)) utils::read.csv(audit_csv_p, stringsAsFactors = FALSE) else data.frame()

pick_h <- function(d, cid) {
  k <- which(d$catalogNumber == cid)[1L]
  list(
    sub = as.character(d$taxonomic_treatment_subtype[k]),
    rec = as.character(d$recommended_determination[k]),
    coh = isTRUE(d$hybrid_zone_cohesion_pass[k]),
    p1 = if (is.na(d$hybrid_parent_1[k])) NULL else as.character(d$hybrid_parent_1[k]),
    p2 = if (is.na(d$hybrid_parent_2[k])) NULL else as.character(d$hybrid_parent_2[k])
  )
}
pick_p <- function(d, cid) {
  k <- which(d$catalogNumber == cid)[1L]
  list(
    sub = as.character(d$taxonomic_treatment_subtype[k]),
    rec = as.character(d$recommended_determination[k]),
    lag = isTRUE(d$portal_metadata_lag_flag[k]),
    gate = if (is.na(d$portal_lag_trigger_gate[k])) NULL else as.character(d$portal_lag_trigger_gate[k]),
    marg = isTRUE(d$is_marginal_range_extension[k]),
    rat = as.character(d$synthesis_rationale[k])
  )
}

write_json(list(
  h_deme1 = pick_h(res_hyb, "HYB_DEME_1"),
  h_deme2 = pick_h(res_hyb, "HYB_DEME_2"),
  h_x = pick_h(res_hyb, "HYB_X_SINGLE"),
  h_g1 = pick_h(res_hyb, "HYB_G1_ADM"),
  h_allo = pick_h(res_hyb, "HYB_FAIL_ALLOPATRIC"),
  h_alloch = pick_h(res_hyb, "HYB_FAIL_ALLOCHRONIC"),
  h_off = pick_h(res_hyb, "HYB_FAIL_OFFCHORD"),
  h_sing = pick_h(res_hyb, "HYB_FAIL_SINGLETON"),
  h_buck = pick_h(res_hyb, "BUCK_CREEK_PURE"),
  intra_can_hyb = isTRUE(attr(res_intra_sup, "can_evaluate_hybrids")),
  intra_n_hyb = sum(res_intra_sup$taxonomic_treatment_subtype == "Hybrid_Introgressant", na.rm = TRUE),
  single_can_hyb = isTRUE(attr(res_single_sup, "can_evaluate_hybrids")),
  single_n_hyb = sum(res_single_sup$taxonomic_treatment_subtype == "Hybrid_Introgressant", na.rm = TRUE),
  p_g1 = pick_p(res_port, "PORTAL_G1_EXS"),
  p_g2 = pick_p(res_port, "PORTAL_G2_BARKLEY"),
  p_g3b = pick_p(res_port, "PORTAL_G3_BIOCLIM"),
  p_g3u = pick_p(res_port, "PORTAL_G3_UNGEOREF"),
  p_marg = pick_p(res_port, "PORTAL_MARG_EXT"),
  p_waif = pick_p(res_port, "PORTAL_WAIF_EXOTIC"),
  audit_exists = file.exists(audit_csv_p),
  audit_cats = if ("catalogNumber" %in% names(audit_df)) as.character(audit_df$catalogNumber) else character(0)
), out_json, auto_unbox = TRUE, na = "null")
"""


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R hybrid/biogeo suite under fast audit or missing Rscript.")
class TestStep07HybridAndBiogeoGates(unittest.TestCase):
    """Contract tests for 5-gate hybrid swarms, intraspecific suppression, marginal extensions, and Portal Lag Gates."""

    res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        with tempfile.TemporaryDirectory(prefix="step07_hyb_bg_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_p), str(out_j), str(tmp)], stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True)
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))

    def test_five_gate_hybrid_swarm_pass_and_veto_conditions(self) -> None:
        r = self.res
        for key in ("h_deme1", "h_deme2", "h_x", "h_g1"):
            self.assertEqual(r[key]["sub"], "Hybrid_Introgressant", f"{key} should pass 5-gate hybrid check")
            self.assertTrue(bool(r[key]["coh"]))
            self.assertEqual({r[key]["p1"], r[key]["p2"]}, {"Packera dubia", "Packera aurea"})
        for key in ("h_allo", "h_alloch", "h_off", "h_sing", "h_buck"):
            self.assertNotEqual(r[key]["sub"], "Hybrid_Introgressant", f"{key} must fail hybrid gates")
            self.assertFalse(bool(r[key]["coh"]))
        self.assertEqual(r["h_buck"]["sub"], "Confirmed_Species")
        self.assertFalse(bool(r["intra_can_hyb"]))
        self.assertEqual(int(r["intra_n_hyb"]), 0)
        self.assertFalse(bool(r["single_can_hyb"]))
        self.assertEqual(int(r["single_n_hyb"]), 0)

    def test_marginal_range_extensions_and_portal_metadata_lag_gates_1_2_3_5(self) -> None:
        r = self.res
        self.assertEqual(r["p_g1"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["p_g1"]["rec"], "Packera aurea")
        self.assertTrue(bool(r["p_g1"]["lag"]))
        self.assertEqual(r["p_g1"]["gate"], "Gate1_Exsiccatae_Conflict")
        self.assertEqual(r["p_g2"]["sub"], "Confirmed_Species")
        self.assertFalse(bool(r["p_g2"]["lag"]))
        self.assertEqual(r["p_g3b"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["p_g3b"]["gate"], "Gate3_Extralimital_Bioclimatic")
        self.assertEqual(r["p_g3u"]["sub"], "Misidentification_Reassign")
        self.assertEqual(r["p_g3u"]["gate"], "Gate3_Ungeoref_Extralimital")
        self.assertEqual(r["p_marg"]["sub"], "Confirmed_Species")
        self.assertTrue(bool(r["p_marg"]["marg"]))
        self.assertFalse(bool(r["p_marg"]["lag"]))
        self.assertIn("Marginal Range Extension Preserved", r["p_marg"]["rat"])
        self.assertEqual(r["p_waif"]["sub"], "Misidentification_Reassign")
        self.assertFalse(bool(r["p_waif"]["marg"]))
        self.assertTrue(bool(r["audit_exists"]))
        for cid in ("PORTAL_G1_EXS", "PORTAL_G3_BIOCLIM", "PORTAL_G3_UNGEOREF"):
            self.assertIn(cid, r["audit_cats"])


if __name__ == "__main__":
    unittest.main()
