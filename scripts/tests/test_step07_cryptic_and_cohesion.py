#!/usr/bin/env python3
"""Warm-R tests for Step 07 GMM anti-inflation, allochronic bimodality, G1 exsiccatae ranking, multi-island ECT, nQuire, and cohesion."""

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

# 1. resolve_gmm_bic_stats() anti-inflation (optimal_K_merged = 1 and Boot_LRT_P_K1_vs_K2 >= 0.05)
g1_p <- file.path(work_dir, "gmm_merged1.csv")
utils::write.csv(data.frame(
  K = c(1L, 4L), BIC = c(-90.0, -112.0), Two_Delta_BIC = c(0.0, 22.0),
  Boot_LRT_P_K1_vs_K2 = c(0.45, 0.45), is_optimal = c(FALSE, TRUE),
  optimal_K = c(4L, 4L), optimal_K_merged = c(1L, 1L), stringsAsFactors = FALSE
), g1_p, row.names = FALSE)
gs_merged1 <- resolve_gmm_bic_stats(list(gmm_bic = g1_p, gmm_summary = ""))

g2_p <- file.path(work_dir, "gmm_lrt_fail.csv")
utils::write.csv(data.frame(
  K = c(1L, 2L, 3L), BIC = c(-95.0, -104.0, -110.0), Two_Delta_BIC = c(0.0, 9.0, 15.0),
  Boot_LRT_P_K1_vs_K2 = c(NA_real_, 0.38, 0.38), is_optimal = c(FALSE, FALSE, TRUE), stringsAsFactors = FALSE
), g2_p, row.names = FALSE)
gs_lrt_fail <- resolve_gmm_bic_stats(list(gmm_bic = g2_p, gmm_summary = ""))

g3_p <- file.path(work_dir, "gmm_valid2.csv")
utils::write.csv(data.frame(
  K = c(1L, 2L), BIC = c(-95.0, -113.5), Two_Delta_BIC = c(0.0, 18.5),
  Boot_LRT_P_K1_vs_K2 = c(0.002, 0.002), is_optimal = c(FALSE, TRUE),
  optimal_K = c(2L, 2L), optimal_K_merged = c(2L, 2L), stringsAsFactors = FALSE
), g3_p, row.names = FALSE)
gs_valid2 <- resolve_gmm_bic_stats(list(gmm_bic = g3_p, gmm_summary = ""))

# 2. Stage 3 allochronic_bimodality_delta_bic > 2.0 (Putative_Cryptic_Species) vs <= 2.0 (Phenological_Ecotype / Subspecies_Ecotype)
mk_allo <- function(dbic_val) {
  data.frame(
    catalogNumber = c(sprintf("AL_C1_%02d", 1:6), sprintf("AL_K2_%02d", 1:6)),
    species_standardized = rep("Packera dubia", 12L),
    determiner_tier = c(rep("Tier_1_Gold", 6L), rep("Tier_2_Silver", 6L)),
    recordedBy = c(paste0("Core", 1:6), "Smith", "Jones", "Brown", "Davis", "Miller", "Wilson"),
    recordNumber = as.character(101:112),
    institutionCode = rep(c("NCU", "DUKE", "US"), length.out = 12L),
    decimalLatitude = c(35.10, 35.12, 35.15, 35.20, 35.25, 35.30, 35.11, 35.22, 35.33, 35.14, 35.26, 35.36),
    decimalLongitude = c(-79.10, -79.15, -79.20, -79.25, -79.30, -79.35, -79.11, -79.22, -79.33, -79.14, -79.26, -79.36),
    phenological_state = rep("anthesis", 12L),
    doy = c(rep(112, 6L), rep(154, 6L)),
    delta_doy_resid = c(rep(0.5, 6L), rep(40.0, 6L)),
    allochronic_bimodality_delta_bic = c(rep(0.2, 6L), rep(dbic_val, 6L)),
    PC1 = c(seq(-2.1, -1.9, length.out = 6L), seq(3.2, 3.8, length.out = 6L)),
    PC2 = c(rep(0.0, 6L), rep(0.2, 6L)),
    can1 = c(seq(-2.1, -1.9, length.out = 6L), seq(3.2, 3.8, length.out = 6L)),
    can2 = c(rep(0.0, 6L), rep(0.2, 6L)),
    cda_predicted_taxon = rep("Packera dubia", 12L),
    cda_posterior_prob = c(rep(0.94, 6L), rep(0.42, 6L)),
    self_morph_mahal_d2 = c(rep(1.0, 6L), rep(28.0, 6L)),
    edaphic_best_fit_taxon = rep("Packera dubia", 12L),
    edaphic_niche_overlap_d = rep(0.82, 12L),
    edaphic_z_score = rep(0.4, 12L),
    gmm_cluster = c(rep(1L, 6L), rep(2L, 6L)),
    morphometric_outlier_flag = c(rep(FALSE, 6L), rep(TRUE, 6L)),
    morph_tier_source = rep("Basal_Rosette", 12L),
    basal_leaf_count = rep(3L, 12L),
    stringsAsFactors = FALSE
  )
}
res_bimod <- apply_taxonomic_decision_matrix(mk_allo(5.4), intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0, opts = list(demo = TRUE))
res_unimod <- apply_taxonomic_decision_matrix(mk_allo(0.8), intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0, opts = list(demo = TRUE))

# 3. G1_Exsiccatae_Isoduplicate loader/replicate ranking & n_g12 summary count
g_dir <- file.path(work_dir, "geno_rank")
dir.create(g_dir, recursive = TRUE, showWarnings = FALSE)
ga_rank <- data.frame(
  genetic_sample_id = c("S_G2_HI", "S_G1E_WIN", "S_P_LO", "S_P_HI", "S_HQ_LO", "S_HQ_HI", "S_SID_Z", "S_SID_A", "S_G3_FLD", "S_G4_TIP"),
  catalogNumber = c("CAT_TIER", "CAT_TIER", "CAT_PROB", "CAT_PROB", "CAT_HQ", "CAT_HQ", "CAT_SID", "CAT_SID", NA_character_, NA_character_),
  species_standardized = rep("Packera dubia", 10L),
  linkage_tier = c("G2_Population_Voucher_Proxy", "G1_Exsiccatae_Isoduplicate", "G1_Exsiccatae_Isoduplicate", "G1_Exsiccatae_Isoduplicate", "G1_Direct_Voucher", "G1_Direct_Voucher", "G1_Direct_Voucher", "G1_Direct_Voucher", "G3_Georeferenced_Unvouchered", "G4_Ungeoreferenced_Tip"),
  geno_quality_pass = rep(TRUE, 10L),
  geno_predicted_taxon = rep("Packera dubia", 10L),
  geno_posterior_prob = c(0.99, 0.91, 0.88, 0.97, 0.95, 0.95, 0.95, 0.95, 0.94, 0.93),
  hq_reads = c(500000L, 100000L, 900000L, 120000L, 100000L, 400000L, 200000L, 200000L, 300000L, 250000L),
  stringsAsFactors = FALSE
)
utils::write.csv(ga_rank, file.path(g_dir, "genomic_assignments.csv"), row.names = FALSE)
df_v_g <- data.frame(catalogNumber = c("CAT_TIER", "CAT_PROB", "CAT_HQ", "CAT_SID"), species_standardized = rep("Packera dubia", 4L), stringsAsFactors = FALSE)
opts_g <- list(output_queue = file.path(g_dir, "triage_queue.csv"), output_report = file.path(g_dir, "rep.md"), genomic_dir = g_dir)
df_g_loaded <- load_step07_genomic_metrics(opts_g, df_v_g)
md_g_lines <- build_genomic_markdown_summary(opts_g, df_g_loaded)

# 4. Dissected multi-island (PC1 = NA, ECT_PC1..3 finite) deme propagation & Tier-0 calibration
df_ect <- data.frame(
  catalogNumber = c(sprintf("MIL_%02d", 1:4), sprintf("DUB_%02d", 1:4), "MIL_UNSEQ_DEME", "MIL_UNSEQ_T0"),
  plant_individual_id = 0L,
  species_standardized = c(rep("Packera millefolium", 4L), rep("Packera dubia", 4L), "Packera millefolium", "Packera millefolium"),
  determiner_tier = c(rep("Tier_1_Gold", 8L), "Tier_3_Bronze", "Tier_3_Bronze"),
  population_deme_id = c("PM1", "PM1", "PM2", "PM2", "PD1", "PD1", "PD2", "PD2", "PM1", "PM3"),
  decimalLatitude = c(35.100, 35.102, 35.400, 35.402, 36.200, 36.202, 36.500, 36.502, 35.103, 35.900),
  decimalLongitude = c(-82.500, -82.501, -82.600, -82.601, -80.200, -80.201, -80.400, -80.401, -82.501, -82.900),
  coordinateUncertaintyInMeters = rep(150.0, 10L),
  morph_tier_source = rep("Basal_Rosette", 10L),
  basal_leaf_count = rep(3L, 10L),
  PC1 = rep(NA_real_, 10L), PC2 = rep(NA_real_, 10L),
  ECT_PC1 = c(4.8, 5.2, 4.9, 5.1, -4.8, -5.2, -4.9, -5.1, 5.05, 4.95),
  ECT_PC2 = c(0.1, -0.1, 0.2, -0.2, 0.1, -0.1, 0.2, -0.2, 0.05, -0.05),
  ECT_PC3 = c(0.05, -0.05, 0.1, -0.1, -0.05, 0.05, -0.1, 0.1, 0.0, 0.02),
  fine_soil_ph = c(rep(6.8, 4L), rep(5.4, 4L), 6.7, 6.8),
  hopkins_adjusted_doy = c(rep(130, 4L), rep(115, 4L), 132, 131),
  is_active_anthesis = rep(TRUE, 10L),
  geno_linkage_tier = c(rep("G1_Direct_Voucher", 8L), NA_character_, NA_character_),
  geno_quality_pass = c(rep(TRUE, 8L), NA, NA),
  is_geno_admixed = rep(FALSE, 10L),
  geno_predicted_taxon = c(rep("Packera millefolium", 4L), rep("Packera dubia", 4L), NA_character_, NA_character_),
  geno_posterior_prob = c(rep(0.96, 8L), NA_real_, NA_real_),
  geno_margin = c(rep(0.88, 8L), NA_real_, NA_real_),
  stringsAsFactors = FALSE
)
ect_deme <- propagate_conditional_deme_priors(df_ect)
ect_t0 <- fit_tier0_calibrated_projection(df_ect)

# 5. Sign-verified nQuire cytotype barriers (2x <= -4.0, 4x >= +4.0 vs wrong signs rejected)
mk_nq <- function(dl_2x, dl_4x) {
  data.frame(
    catalogNumber = c("NQ_2X_1", "NQ_2X_2", "NQ_4X_1", "NQ_4X_2", "NQ_4X_U"),
    species_standardized = rep("Packera dubia", 5L),
    gmm_cluster = c(1L, 1L, 2L, 2L, 2L),
    geno_linkage_tier = c("G1_Direct_Voucher", "G1_Direct_Voucher", "G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", NA_character_),
    geno_quality_pass = c(TRUE, TRUE, TRUE, TRUE, NA),
    is_geno_admixed = rep(FALSE, 5L),
    inferred_cytotype = c("2x", "2x", "4x", "4x", NA_character_),
    nquire_delta_logl = c(dl_2x, dl_2x, dl_4x, dl_4x, NA_real_),
    stringsAsFactors = FALSE
  )
}
nq_valid <- run_step07_cytology_integration(mk_nq(-7.5, 8.2), opts = list(cytological_records_df = data.frame(), taxon_cytology_summary_df = data.frame()))
nq_wrong <- run_step07_cytology_integration(mk_nq(7.5, -8.2), opts = list(cytological_records_df = data.frame(), taxon_cytology_summary_df = data.frame()))

# 6. Duplicate single-gathering cohesion block & allopatric spatial-cline block
df_dup_blk <- mk_allo(6.0)
df_dup_blk$recordedBy[7:12] <- "SameCollector"
df_dup_blk$recordNumber[7:12] <- "999"
df_dup_blk$institutionCode[7:12] <- "NCU"
df_dup_blk$decimalLatitude[7:12] <- 35.11; df_dup_blk$decimalLongitude[7:12] <- -79.11
res_dup_blk <- apply_taxonomic_decision_matrix(df_dup_blk, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0, opts = list(demo = TRUE))

df_cline <- data.frame(
  catalogNumber = sprintf("CLN_%02d", 1:12),
  species_standardized = rep("Packera dubia", 12L),
  determiner_tier = c(rep("Tier_1_Gold", 6L), rep("Tier_2_Silver", 6L)),
  recordedBy = paste0("Coll", 1:12), recordNumber = as.character(1:12),
  institutionCode = rep(c("NCU", "DUKE", "US"), 4L),
  decimalLatitude = c(33.0, 33.4, 33.8, 34.2, 34.6, 35.0, 38.0, 38.6, 39.2, 39.8, 40.4, 41.0),
  decimalLongitude = c(-82.0, -81.8, -81.6, -81.4, -81.2, -81.0, -78.0, -77.6, -77.2, -76.8, -76.4, -76.0),
  phenological_state = rep("anthesis", 12L),
  doy = rep(120, 12L), delta_doy_resid = rep(1.0, 12L),
  PC1 = seq(-3.0, 4.5, length.out = 12L),
  PC2 = rep(0.0, 12L),
  can1 = seq(-3.0, 4.5, length.out = 12L),
  can2 = rep(0.0, 12L),
  cda_predicted_taxon = rep("Packera dubia", 12L),
  cda_posterior_prob = c(rep(0.92, 6L), rep(0.45, 6L)),
  edaphic_best_fit_taxon = rep("Packera dubia", 12L),
  edaphic_niche_overlap_d = c(rep(0.85, 6L), rep(0.20, 6L)),
  edaphic_z_score = c(rep(0.4, 6L), rep(3.2, 6L)),
  gmm_cluster = c(rep(1L, 6L), rep(2L, 6L)),
  morphometric_outlier_flag = c(rep(FALSE, 6L), rep(TRUE, 6L)),
  morph_tier_source = rep("Basal_Rosette", 12L),
  basal_leaf_count = rep(3L, 12L),
  stringsAsFactors = FALSE
)
res_cline <- apply_taxonomic_decision_matrix(df_cline, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0, opts = list(demo = TRUE))

write_json(list(
  gs_merged1 = gs_merged1, gs_lrt_fail = gs_lrt_fail, gs_valid2 = gs_valid2,
  bimod_k2_call = as.character(res_bimod$taxonomic_status_call[7L]),
  bimod_k2_sub = as.character(res_bimod$taxonomic_treatment_subtype[7L]),
  unimod_k2_call = as.character(res_unimod$taxonomic_status_call[7L]),
  unimod_k2_sub = as.character(res_unimod$taxonomic_treatment_subtype[7L]),
  unimod_k2_rat = as.character(res_unimod$synthesis_rationale[7L]),
  rank_sids = stats::setNames(as.list(as.character(df_g_loaded$genetic_sample_id)), df_g_loaded$catalogNumber),
  rank_tiers = stats::setNames(as.list(as.character(df_g_loaded$geno_linkage_tier)), df_g_loaded$catalogNumber),
  md_g_text = paste(md_g_lines, collapse = "\n"),
  ect_deme_tier = as.character(ect_deme$genomic_propagation_tier[ect_deme$catalogNumber == "MIL_UNSEQ_DEME"][1L]),
  ect_deme_pred = as.character(ect_deme$geno_predicted_taxon[ect_deme$catalogNumber == "MIL_UNSEQ_DEME"][1L]),
  ect_t0_flag = as.integer(ect_t0$tier0_abstention_flag[ect_t0$catalogNumber == "MIL_UNSEQ_T0"][1L]),
  ect_t0_tax = as.character(ect_t0$tier0_calibrated_taxon[ect_t0$catalogNumber == "MIL_UNSEQ_T0"][1L]),
  nq_valid_pass = as.integer(nq_valid$cytotype_reproductive_barrier_pass[5L]),
  nq_valid_pl = as.integer(nq_valid$voucher_ploidy_level[1:4]),
  nq_wrong_pass = as.integer(nq_wrong$cytotype_reproductive_barrier_pass[5L]),
  nq_wrong_all_na = all(is.na(nq_wrong$voucher_ploidy_level)),
  dup_coh_pass = isTRUE(res_dup_blk$cryptic_cluster_cohesion_pass[7L]),
  dup_k2_sub = as.character(res_dup_blk$taxonomic_treatment_subtype[7L]),
  cline_nonclinal = isTRUE(res_cline$cryptic_non_clinal_pass[7L]),
  cline_k2_sub = as.character(res_cline$taxonomic_treatment_subtype[7L]),
  cline_k2_rat = as.character(res_cline$synthesis_rationale[7L])
), out_json, auto_unbox = TRUE, na = "null")
"""


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R cryptic/cohesion suite under fast audit or missing Rscript.")
class TestStep07CrypticAndCohesion(unittest.TestCase):
    """Contract tests for GMM anti-inflation, allochronic bimodality, genomic ranking, multi-island ECT, nQuire, and cohesion."""

    res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        with tempfile.TemporaryDirectory(prefix="step07_cry_coh_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_p), str(out_j), str(tmp)], stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True)
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.res = json.loads(out_j.read_text(encoding="utf-8"))

    def test_gmm_bic_anti_inflation_and_allochronic_bimodality_routing(self) -> None:
        r = self.res
        self.assertEqual((r["gs_merged1"]["opt_k"], r["gs_merged1"]["two_delta_bic"], r["gs_merged1"]["raw_opt_k"]), (1, 0.0, 4))
        self.assertEqual((r["gs_lrt_fail"]["opt_k"], r["gs_lrt_fail"]["two_delta_bic"], r["gs_lrt_fail"]["raw_opt_k"]), (1, 0.0, 3))
        self.assertEqual((r["gs_valid2"]["opt_k"], r["gs_valid2"]["two_delta_bic"]), (2, 18.5))
        self.assertEqual(r["bimod_k2_call"], "Putative_Cryptic_Species")
        self.assertEqual(r["bimod_k2_sub"], "Putative_Cryptic_Species")
        self.assertEqual(r["unimod_k2_call"], "Phenological_Ecotype")
        self.assertEqual(r["unimod_k2_sub"], "Subspecies_Ecotype")
        self.assertIn("Unimodal Allochronic Cline", r["unimod_k2_rat"])

    def test_g1_exsiccatae_replicate_ranking_and_markdown_g12_count(self) -> None:
        r = self.res
        self.assertEqual(r["rank_sids"]["CAT_TIER"], "S_G1E_WIN")
        self.assertEqual(r["rank_tiers"]["CAT_TIER"], "G1_Exsiccatae_Isoduplicate")
        self.assertEqual(r["rank_sids"]["CAT_PROB"], "S_P_HI")
        self.assertEqual(r["rank_sids"]["CAT_HQ"], "S_HQ_HI")
        self.assertEqual(r["rank_sids"]["CAT_SID"], "S_SID_A")
        self.assertIn("8 voucher-linked (`G1`/`G2`)", r["md_g_text"])

    def test_multi_island_ect_propagation_sign_verified_nquire_and_cohesion_clines(self) -> None:
        r = self.res
        self.assertEqual(r["ect_deme_tier"], "G2_Conditional_Deme_Hull")
        self.assertEqual(r["ect_deme_pred"], "Packera millefolium")
        self.assertEqual(int(r["ect_t0_flag"]), 0)
        self.assertEqual(r["ect_t0_tax"], "Packera millefolium")
        self.assertEqual(int(r["nq_valid_pass"]), 1)
        self.assertEqual(r["nq_valid_pl"], [2, 2, 4, 4])
        self.assertIn(r["nq_wrong_pass"], (0, None))
        self.assertTrue(bool(r["nq_wrong_all_na"]))
        self.assertFalse(bool(r["dup_coh_pass"]))
        self.assertEqual(r["dup_k2_sub"], "Subspecies_Ecotype")
        self.assertFalse(bool(r["cline_nonclinal"]))
        self.assertEqual(r["cline_k2_sub"], "Subspecies_Ecotype")
        self.assertIn("Allopatric Spatial Cline", r["cline_k2_rat"])


if __name__ == "__main__":
    unittest.main()
