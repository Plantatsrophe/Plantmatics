#!/usr/bin/env python3
"""
Warm-R contract tests for Step 07 four-tier genomic propagation (`Scenarios a-e`):
  (a) Tier 1 exsiccatae isoduplicate propagation + mixed-mount/type/multi-mount vetoes;
  (b) Tier 2 conditional deme hull + edaphic/pheno/uncertainty/sympatry vetoes + non-override;
  (c) Tier_0 DNA anchor eviction + LOPO-CV GDA & overlap/OSR abstention;
  (d) Tier 4 cluster-wide coalescent gdi governance (0.09 oversplit block & 0.84 support);
  (e) Zero-genomics strict invariance across all decision-matrix columns.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R genomic propagation contract under fast audit or missing Rscript.")
class TestStep07GenomicPropagationContract(unittest.TestCase):
    """Comprehensive Warm-R verification of the 4-tier Step 07 genomic propagation contract."""

    def test_batched_step07_genomic_propagation_four_tiers_and_invariance(self) -> None:
        """Warm-R batched contract test covering Scenarios (a) through (e) of Step 07 genomic propagation."""
        r_code = """
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        source("scripts/analysis/07_triage_dashboard_synthesis.R")
        if (!exists("run_step07_genomic_propagation", mode = "function")) {
          source("scripts/analysis/step07_genomic_propagation.R")
        }

        # ------------------------------------------------------------------
        # Scenario (a): Tier 1 Exsiccatae Propagation + Mixed-Mount/Type Veto
        # ------------------------------------------------------------------
        n_bg <- 36L
        df_exs <- data.frame(
          catalogNumber = c(sprintf("BG_%02d", 1:n_bg), "EXS_DONOR", "EXS_DUP_OK", "EXS_MIXED_VETO", "EXS_TYPE_VETO", "EXS_EXP_DONOR", "EXS_EXP_MULTI_MO", "ANC_AUR_1", "ANC_AUR_2"),
          plant_individual_id = 0L,
          species_standardized = c(rep(c("Packera dubia", "Packera anonyma"), each = 18L), rep("Packera dubia", 4L), rep("Packera aurea", 4L)),
          determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", basal_leaf_count = 3L,
          n_plant_individuals = c(rep(1L, n_bg + 5L), 2L, 1L, 1L),
          morphometric_outlier_flag = c(rep(FALSE, n_bg + 5L), TRUE, FALSE, FALSE),
          type_status = c(rep("none", n_bg + 3L), "ISOTYPE", rep("none", 4L)),
          is_type_specimen = c(rep(FALSE, n_bg + 3L), TRUE, rep(FALSE, 4L)),
          recordedBy = c(sprintf("Coll_%02d", 1:n_bg), "J. R. Radford", "Radford, John R.", "A. Collector", "J. R. Radford", "E. Weakley", "E. Weakley", "A. Ray", "B. Ray"),
          recordNumber = c(sprintf("%d", 1000L + seq_len(n_bg)), "4421", " 4421 ", "4421", "4421", "901", "901", "801", "802"),
          eventDate = c(rep("1990-05-01", n_bg), "1968-05-14", "1968-05-14T12:00:00Z", "1968-05-14", "1968-05-14", rep("1995-06-01", 4L)),
          exsiccatae_key = c(rep(NA_character_, n_bg + 4L), "explicit_ex_99", "explicit_ex_99", NA_character_, NA_character_),
          PC1 = c(seq(-1.5, 1.5, length.out = n_bg), 0.10, 0.14, 0.10, 0.12, 0.20, 0.22, -0.50, -0.45),
          PC2 = c(rep(c(-0.8, 0.8), 18L), 0.10, 0.12, 12.00, 0.11, 0.15, 0.16, -0.20, -0.25),
          cda_predicted_taxon = c(rep(c("Packera dubia", "Packera anonyma"), each = 18L), rep("Packera dubia", 4L), rep("Packera aurea", 4L)),
          cda_posterior_prob = 0.88,
          edaphic_best_fit_taxon = c(rep(c("Packera dubia", "Packera anonyma"), each = 18L), rep("Packera dubia", 4L), rep("Packera aurea", 4L)),
          rf_predicted_prob = 0.85, delta_doy_resid = 1.0, is_active_anthesis = TRUE,
          geno_linkage_tier = c(c("G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, n_bg - 2L)), "G1_Direct_Voucher", NA_character_, NA_character_, NA_character_, "G1_Direct_Voucher", NA_character_, "G1_Direct_Voucher", "G1_Direct_Voucher"),
          geno_quality_pass = c(c(TRUE, TRUE, rep(NA, n_bg - 2L)), TRUE, NA, NA, NA, TRUE, NA, TRUE, TRUE),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c(c("Packera dubia", "Packera anonyma", rep(NA_character_, n_bg - 2L)), "Packera serpenticola", NA_character_, NA_character_, NA_character_, "Packera aurea", NA_character_, "Packera aurea", "Packera aurea"),
          geno_posterior_prob = c(c(0.95, 0.94, rep(NA_real_, n_bg - 2L)), 0.96, NA_real_, NA_real_, NA_real_, 0.95, NA_real_, 0.95, 0.95),
          geno_margin = 0.85, stringsAsFactors = FALSE
        )
        df_exs$recordedBy[df_exs$catalogNumber == "EXS_MIXED_VETO"] <- "J. R. Radford"
        res_a <- propagate_exsiccatae_isoduplicates(df_exs)
        i_dup <- which(res_a$catalogNumber == "EXS_DUP_OK")
        i_mix <- which(res_a$catalogNumber == "EXS_MIXED_VETO")
        i_typ <- which(res_a$catalogNumber == "EXS_TYPE_VETO")
        i_mul <- which(res_a$catalogNumber == "EXS_EXP_MULTI_MO")
        stopifnot(identical(res_a$genomic_propagation_tier[i_dup], "G1_Exsiccatae_Isoduplicate"))
        stopifnot(identical(res_a$geno_linkage_tier[i_dup], "G1_Exsiccatae_Isoduplicate"))
        stopifnot(identical(res_a$geno_predicted_taxon[i_dup], "Packera serpenticola"))
        stopifnot(is.na(res_a$genomic_propagation_tier[i_mix]) && is.na(res_a$geno_predicted_taxon[i_mix]))
        stopifnot(is.na(res_a$genomic_propagation_tier[i_typ]) && is.na(res_a$geno_predicted_taxon[i_typ]))
        stopifnot(is.na(res_a$genomic_propagation_tier[i_mul]) && is.na(res_a$geno_predicted_taxon[i_mul]))

        sc_a <- compute_probabilistic_consensus_scores(res_a)
        wg_donor <- sc_a$w_geno_kappa[sc_a$catalogNumber == "EXS_DONOR"]
        wg_dup <- sc_a$w_geno_kappa[sc_a$catalogNumber == "EXS_DUP_OK"]
        stopifnot(is.finite(wg_donor) && is.finite(wg_dup) && wg_donor > 0)
        stopifnot(abs((wg_dup / wg_donor) - 0.95) < 1e-3)

        # ------------------------------------------------------------------
        # Scenario (b): Tier 2 Deme Hull + Buck Creek Sympatry Veto + Non-Override
        # ------------------------------------------------------------------
        df_deme <- data.frame(
          catalogNumber = c(
            sprintf("SERP_BG_%02d", 1:6), sprintf("AUR_BG_%02d", 1:6),
            "SERP_SITE1_ANC", "DEME_MATCH_OK", "DEME_VETO_EDAPH", "DEME_VETO_PHENO",
            "DEME_VETO_UNC", "SYMP_SERP_ANC", "SYMP_AUR_ANC", "DEME_VETO_SYMP"
          ),
          plant_individual_id = 0L,
          species_standardized = c(
            rep("Packera serpenticola", 6L), rep("Packera aurea", 6L),
            "Packera serpenticola", "Packera aurea", rep("Packera serpenticola", 4L),
            "Packera aurea", "Packera serpenticola"
          ),
          determiner_tier = c(rep("Tier_1_Gold", 13L), rep("Tier_3_Bronze", 4L), "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze"),
          morph_tier_source = "Basal_Rosette", basal_leaf_count = 3L,
          decimalLatitude = c(seq(35.60, 35.65, length.out = 6L), seq(35.80, 35.85, length.out = 6L), 35.100, 35.108, 35.108, 35.108, 35.108, 35.300, 35.301, 35.302),
          decimalLongitude = c(rep(-83.20, 6L), rep(-82.50, 6L), -83.800, -83.808, -83.808, -83.808, -83.808, -83.500, -83.500, -83.501),
          coordinateUncertaintyInMeters = c(rep(200, 16L), 2500, 200, 200, 200),
          is_azonal_lithology = c(rep(1L, 6L), rep(0L, 6L), 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L),
          fine_soil_ph = c(rep(7.1, 6L), rep(5.5, 6L), 7.1, 7.0, 5.4, 7.1, 7.1, 7.1, 7.1, 7.1),
          type_edaphic_z_score = 0.3,
          hopkins_adjusted_doy = c(rep(125, 6L), rep(120, 6L), 125, 132, 125, 175, 125, 125, 125, 125),
          is_active_anthesis = TRUE,
          PC1 = c(1.4, 1.6, 1.3, 1.7, 1.5, 1.5, -1.4, -1.6, -1.3, -1.7, -1.5, -1.5, 1.5, 1.55, 1.52, 1.52, 1.52, 1.5, -1.5, 0.20),
          PC2 = c(-0.4, 0.4, -0.3, 0.3, -0.2, 0.2, -0.4, 0.4, -0.3, 0.3, -0.2, 0.2, 0.0, 0.05, 0.02, 0.02, 0.02, 0.0, 0.0, 0.0),
          cda_predicted_taxon = c(rep("Packera serpenticola", 6L), rep("Packera aurea", 6L), "Packera serpenticola", "Packera aurea", rep("Packera serpenticola", 4L), "Packera aurea", "Packera serpenticola"),
          cda_posterior_prob = c(rep(0.92, 13L), 0.72, rep(0.92, 6L)),
          edaphic_best_fit_taxon = c(rep("Packera serpenticola", 6L), rep("Packera aurea", 6L), "Packera serpenticola", "Packera aurea", rep("Packera serpenticola", 4L), "Packera aurea", "Packera serpenticola"),
          rf_predicted_prob = c(rep(0.90, 13L), 0.75, rep(0.90, 6L)),
          geno_linkage_tier = c("G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 4L), "G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 4L), "G1_Direct_Voucher", rep(NA_character_, 4L), "G1_Direct_Voucher", "G1_Direct_Voucher", NA_character_),
          geno_quality_pass = c(TRUE, TRUE, rep(NA, 4L), TRUE, TRUE, rep(NA, 4L), TRUE, rep(NA, 4L), TRUE, TRUE, NA),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c("Packera serpenticola", "Packera serpenticola", rep(NA_character_, 4L), "Packera aurea", "Packera aurea", rep(NA_character_, 4L), "Packera serpenticola", rep(NA_character_, 4L), "Packera serpenticola", "Packera aurea", NA_character_),
          geno_posterior_prob = c(0.94, 0.93, rep(NA_real_, 4L), 0.94, 0.93, rep(NA_real_, 4L), 0.92, rep(NA_real_, 4L), 0.92, 0.92, NA_real_),
          geno_margin = 0.82, stringsAsFactors = FALSE
        )
        res_b <- propagate_conditional_deme_priors(df_deme)
        ib_ok <- which(res_b$catalogNumber == "DEME_MATCH_OK")
        stopifnot(identical(res_b$genomic_propagation_tier[ib_ok], "G2_Conditional_Deme_Hull"))
        stopifnot(identical(res_b$geno_linkage_tier[ib_ok], "G2_Conditional_Deme_Hull"))
        stopifnot(identical(res_b$geno_predicted_taxon[ib_ok], "Packera serpenticola"))
        stopifnot(is.finite(res_b$deme_mahal_d2[ib_ok]))
        stopifnot(is.finite(res_b$geno_posterior_prob[ib_ok]) && res_b$geno_posterior_prob[ib_ok] >= 0.65 && res_b$geno_posterior_prob[ib_ok] <= 0.90)
        for (vc in c("DEME_VETO_EDAPH", "DEME_VETO_PHENO", "DEME_VETO_UNC", "DEME_VETO_SYMP")) {
          iv <- which(res_b$catalogNumber == vc)
          stopifnot(is.na(res_b$genomic_propagation_tier[iv]) && is.na(res_b$geno_predicted_taxon[iv]))
        }
        sc_b <- compute_probabilistic_consensus_scores(res_b)
        stopifnot(abs((sc_b$w_geno_kappa[ib_ok] / sc_b$w_geno_kappa[sc_b$catalogNumber == "SERP_SITE1_ANC"]) - 0.75) < 1e-3)
        dec_b <- apply_taxonomic_decision_matrix(res_b)
        stopifnot(dec_b$taxonomic_treatment_subtype[dec_b$catalogNumber == "DEME_MATCH_OK"] != "Misidentification_Reassign")

        # ------------------------------------------------------------------
        # Scenario (c): Tier_0 Anchor Eviction + LOPO-CV GDA & Overlap Abstention
        # ------------------------------------------------------------------
        df_t0 <- data.frame(
          catalogNumber = c(sprintf("SERP_T0_%02d", 1:4), sprintf("DUB_T0_%02d", 1:4), sprintf("ANO_T0_%02d", 1:4), "MISID_GOLD_EVICT", "F1_ADMIX_EVICT", "UNSEQ_SERP", "UNSEQ_OVERLAP", "UNSEQ_OSR_ABST"),
          plant_individual_id = 0L,
          species_standardized = c(rep("Packera serpenticola", 4L), rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera dubia", 5L)),
          determiner_tier = c(rep("Tier_2_Silver", 12L), "Tier_1_Gold", "Tier_1_Gold", rep("Tier_3_Bronze", 3L)),
          morph_tier_source = "Basal_Rosette", basal_leaf_count = 3L,
          population_deme_id = c(rep(c("P1", "P2"), 6L), "P3", "P3", "P4", "P4", "P4"),
          PC1 = c(5.8, 6.2, 6.0, 5.9, -2.1, -1.9, -2.0, -2.05, -1.95, -2.05, -2.0, -2.1, 6.1, -0.5, 6.0, -2.02, 22.0),
          PC2 = c(5.9, 6.1, 6.0, 6.05, 0.1, -0.1, 0.0, 0.05, -0.05, 0.05, 0.0, -0.1, 6.0, 0.0, 6.0, 0.0, 22.0),
          hopkins_adjusted_doy = c(rep(125, 4L), rep(115, 4L), rep(116, 4L), 210, 120, 125, 115, 120),
          is_active_anthesis = TRUE,
          geno_linkage_tier = c(rep("G1_Direct_Voucher", 14L), rep(NA_character_, 3L)),
          geno_quality_pass = c(rep(TRUE, 14L), rep(NA, 3L)),
          is_geno_admixed = c(rep(FALSE, 13L), TRUE, rep(FALSE, 3L)),
          geno_predicted_taxon = c(rep("Packera serpenticola", 4L), rep("Packera dubia", 4L), rep("Packera anonyma", 4L), "Packera serpenticola", "Packera dubia", rep(NA_character_, 3L)),
          geno_posterior_prob = c(rep(0.96, 13L), 0.55, rep(NA_real_, 3L)),
          geno_margin = c(rep(0.88, 13L), 0.10, rep(NA_real_, 3L)),
          interspecific_ancestry_q12 = c(rep(0.02, 13L), 0.72, rep(NA_real_, 3L)), stringsAsFactors = FALSE
        )
        t0_res <- resolve_tier0_anchor_overrides(df_t0)
        i_mg <- which(df_t0$catalogNumber == "MISID_GOLD_EVICT")
        i_f1 <- which(df_t0$catalogNumber == "F1_ADMIX_EVICT")
        stopifnot(isTRUE(t0_res$evict_from_given_anchor[i_mg]) && isTRUE(t0_res$is_tier0_dna_anchor[i_mg]) && t0_res$tier0_weight[i_mg] == 4.0)
        stopifnot(isTRUE(t0_res$evict_from_given_anchor[i_f1]) && !isTRUE(t0_res$is_tier0_dna_anchor[i_f1]) && is.na(t0_res$tier0_weight[i_f1]))
        # Passive heterotypic type exclusion from Tier_0 DNA anchor override
        df_t0_het <- df_t0[1L, , drop = FALSE]
        df_t0_het$catalogNumber <- "HET_TYPE_G1"; df_t0_het$is_type_specimen <- TRUE; df_t0_het$type_status <- "HOLOTYPE"
        df_t0_het$anchor_tier_source <- "Passive_Heterotypic_Type_Audit"; df_t0_het$concept_synonymy_status <- "Heterotypic_Synonym"
        t0_het <- resolve_tier0_anchor_overrides(df_t0_het)
        stopifnot(!isTRUE(t0_het$is_tier0_dna_anchor[1L]), !isTRUE(t0_het$evict_from_given_anchor[1L]))
        # G1_Exsiccatae_Isoduplicate chains into conditional deme hull propagation
        df_chain <- df_deme[1:14, , drop = FALSE]
        df_chain$geno_linkage_tier[df_chain$catalogNumber == "SERP_SITE1_ANC"] <- "G1_Exsiccatae_Isoduplicate"
        res_chain <- propagate_conditional_deme_priors(df_chain)
        ic_ok <- which(res_chain$catalogNumber == "DEME_MATCH_OK")
        stopifnot(identical(res_chain$genomic_propagation_tier[ic_ok], "G2_Conditional_Deme_Hull"),
                  identical(res_chain$geno_predicted_taxon[ic_ok], "Packera serpenticola"))

        proj_c <- fit_tier0_calibrated_projection(df_t0)
        iu_s <- which(proj_c$catalogNumber == "UNSEQ_SERP")
        iu_o <- which(proj_c$catalogNumber == "UNSEQ_OVERLAP")
        iu_osr <- which(proj_c$catalogNumber == "UNSEQ_OSR_ABST")
        stopifnot(identical(proj_c$tier0_abstention_flag[iu_s], 0L), identical(proj_c$tier0_calibrated_taxon[iu_s], "Packera serpenticola"))
        stopifnot(identical(proj_c$tier0_abstention_flag[iu_o], 1L), identical(proj_c$tier0_calibrated_taxon[iu_o], "Sensu_Lato_DNA_Required"))
        stopifnot(identical(proj_c$tier0_abstention_flag[iu_osr], 1L), identical(proj_c$tier0_calibrated_taxon[iu_osr], "Sensu_Lato_DNA_Required"))

        dec_t0 <- apply_taxonomic_decision_matrix(proj_c)
        stopifnot(abs(dec_t0$delta_doy_resid[dec_t0$catalogNumber == "DUB_T0_01"]) < 3.0)
        stopifnot(grepl("[Tier_0 Abstention: Phenotypic overlap across genetically distinct lineages; DNA required]", dec_t0$synthesis_rationale[dec_t0$catalogNumber == "UNSEQ_OVERLAP"], fixed = TRUE))

        # ------------------------------------------------------------------
        # Scenario (d): Tier 4 Cluster-Wide gdi Governance (0.09 block & 0.84 support)
        # ------------------------------------------------------------------
        n_k2 <- 15L
        df_k2 <- data.frame(
          catalogNumber = c(sprintf("DUB_CORE_%02d", 1:5), sprintf("K2_ANC_%02d", 1:3), sprintf("K2_UNSEQ_%02d", 1:n_k2)),
          plant_individual_id = 0L,
          institutionCode = c(rep("NCU", 5L), rep("NCU", 3L), rep(c("NCU", "DUKE", "US"), 5L)),
          recordedBy = c(rep("CoreColl", 5L), "AncA", "AncB", "AncA", rep(c("Smith", "Jones", "Brown"), 5L)),
          recordNumber = sprintf("%d", 500L + seq_len(8L + n_k2)),
          species_standardized = "Packera dubia",
          determiner_tier = c(rep("Tier_1_Gold", 8L), rep("Tier_3_Bronze", n_k2)),
          morph_tier_source = "Basal_Rosette", basal_leaf_count = 3L,
          gmm_cluster = c(rep(1L, 5L), rep(2L, 3L + n_k2)),
          morphometric_outlier_flag = c(rep(FALSE, 5L), rep(TRUE, 3L + n_k2)),
          PC1 = c(seq(-0.3, 0.3, length.out = 5L), 2.3, 2.5, 2.4, seq(2.2, 2.6, length.out = n_k2)),
          PC2 = c(rep(0.0, 5L), 0.1, -0.1, 0.0, rep(c(0.05, -0.05, 0.0), 5L)),
          cda_predicted_taxon = "Packera dubia", cda_posterior_prob = c(rep(0.92, 5L), rep(0.45, 3L + n_k2)),
          edaphic_best_fit_taxon = "Packera dubia", edaphic_niche_overlap_d = c(rep(0.85, 5L), rep(0.15, 3L + n_k2)),
          type_edaphic_z_score = c(rep(0.4, 5L), rep(3.2, 3L + n_k2)),
          decimalLatitude = c(rep(35.80, 5L), 35.81, 35.82, 35.95, rep(c(35.81, 35.95, 36.12, 36.30, 36.48), 3L)),
          decimalLongitude = c(rep(-79.00, 5L), -79.01, -79.02, -79.15, rep(c(-79.01, -79.15, -79.35, -79.55, -79.75), 3L)),
          hopkins_adjusted_doy = c(rep(115, 5L), rep(155, 3L + n_k2)), delta_doy_resid = c(rep(1.0, 5L), rep(40.0, 3L + n_k2)),
          is_active_anthesis = TRUE, min_congener_mahal_d2 = c(rep(45.0, 5L), rep(30.0, 3L + n_k2)), open_set_p_value = 0.001,
          geno_linkage_tier = c(rep(NA_character_, 5L), rep("G1_Direct_Voucher", 3L), rep(NA_character_, n_k2)),
          geno_quality_pass = c(rep(NA, 5L), rep(TRUE, 3L), rep(NA, n_k2)), is_geno_admixed = FALSE,
          geno_predicted_taxon = c(rep(NA_character_, 5L), rep("Packera dubia", 3L), rep(NA_character_, n_k2)),
          geno_posterior_prob = c(rep(NA_real_, 5L), rep(0.95, 3L), rep(NA_real_, n_k2)),
          geno_margin = c(rep(NA_real_, 5L), rep(0.85, 3L), rep(NA_real_, n_k2)),
          coalescent_gdi = c(rep(NA_real_, 5L), rep(0.09, 3L), rep(NA_real_, n_k2)), stringsAsFactors = FALSE
        )
        prop_k2 <- run_step07_genomic_propagation(df_k2)
        unseq_k2 <- which(grepl("^K2_UNSEQ_", prop_k2$catalogNumber))
        stopifnot(all(is.finite(prop_k2$cluster_propagated_gdi[unseq_k2])), all(abs(prop_k2$cluster_propagated_gdi[unseq_k2] - 0.09) < 1e-4))
        dec_k2 <- apply_taxonomic_decision_matrix(prop_k2, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0)
        stopifnot(!any(dec_k2$taxonomic_treatment_subtype[unseq_k2] == "Putative_Cryptic_Species"))
        stopifnot(all(grepl("[Cluster Coalescent Oversplit Block: propagated gdi=0.09]", dec_k2$synthesis_rationale[unseq_k2], fixed = TRUE)))

        df_nov <- data.frame(
          catalogNumber = c(sprintf("NOV_DUB_%02d", 1:4), sprintf("NOV_ANO_%02d", 1:4), sprintf("NOV_ANC_%02d", 1:4), "NOVEL_UNSEQ_OK", "NOVEL_UNSEQ_ABST"),
          plant_individual_id = 0L,
          species_standardized = c(rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera dubia", 6L)),
          determiner_tier = c(rep("Tier_2_Silver", 12L), "Tier_3_Bronze", "Tier_3_Bronze"),
          morph_tier_source = "Basal_Rosette", basal_leaf_count = 3L,
          gmm_cluster = c(rep(1L, 8L), rep(3L, 6L)), morphometric_outlier_flag = c(rep(FALSE, 8L), rep(TRUE, 6L)),
          population_deme_id = c(rep(c("PA", "PB"), 6L), "PC", "PC"),
          PC1 = c(-3.1, -2.9, -3.0, -3.05, -2.95, -3.05, -3.0, -2.98, 4.9, 5.1, 5.0, 5.05, 5.0, -3.0),
          PC2 = c(-0.05, 0.05, 0.0, 0.02, 0.02, -0.02, 0.0, 0.04, 4.9, 5.1, 5.0, 4.95, 5.0, 0.0),
          cda_predicted_taxon = c(rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera dubia", 6L)),
          cda_posterior_prob = c(rep(0.90, 8L), rep(0.45, 6L)),
          edaphic_best_fit_taxon = "Packera dubia", edaphic_niche_overlap_d = c(rep(0.80, 8L), rep(0.55, 6L)),
          type_edaphic_z_score = 0.5, delta_doy_resid = 4.0, is_active_anthesis = TRUE,
          min_congener_mahal_d2 = c(rep(35.0, 8L), rep(30.0, 6L)), open_set_p_value = 0.001,
          geno_linkage_tier = c(rep("G1_Direct_Voucher", 12L), NA_character_, NA_character_),
          geno_quality_pass = c(rep(TRUE, 12L), NA, NA), is_geno_admixed = FALSE,
          geno_predicted_taxon = c(rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera sp. nov. A", 4L), NA_character_, NA_character_),
          geno_posterior_prob = c(rep(0.96, 12L), NA_real_, NA_real_), geno_margin = c(rep(0.88, 12L), NA_real_, NA_real_),
          coalescent_gdi = c(rep(0.50, 8L), rep(0.84, 4L), NA_real_, NA_real_), stringsAsFactors = FALSE
        )
        prop_nov <- run_step07_genomic_propagation(df_nov)
        i_nok <- which(prop_nov$catalogNumber == "NOVEL_UNSEQ_OK")
        i_nab <- which(prop_nov$catalogNumber == "NOVEL_UNSEQ_ABST")
        stopifnot(identical(prop_nov$tier0_abstention_flag[i_nok], 0L), identical(prop_nov$tier0_calibrated_taxon[i_nok], "Packera sp. nov. A"))
        stopifnot(abs(prop_nov$cluster_propagated_gdi[i_nok] - 0.84) < 1e-4, identical(prop_nov$tier0_abstention_flag[i_nab], 1L), is.na(prop_nov$cluster_propagated_gdi[i_nab]))

        dec_nov <- apply_taxonomic_decision_matrix(prop_nov, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0)
        stopifnot(identical(dec_nov$taxonomic_treatment_subtype[i_nok], "Putative_Cryptic_Species"))
        stopifnot(grepl("[Candidate Member via Cluster-Propagated gdi=0.84; confirm by targeted sequencing]", dec_nov$synthesis_rationale[i_nok], fixed = TRUE))
        stopifnot(dec_nov$taxonomic_treatment_subtype[i_nab] != "Putative_Cryptic_Species")
        stopifnot(grepl("[Tier_0 Abstention: Phenotypic overlap across genetically distinct lineages; DNA required]", dec_nov$synthesis_rationale[i_nab], fixed = TRUE))

        # ------------------------------------------------------------------
        # Scenario (e): Zero-Genomics Invariance
        # ------------------------------------------------------------------
        df_zero <- data.frame(
          catalogNumber = sprintf("Z_%02d", 1:6), plant_individual_id = 0L,
          species_standardized = rep(c("Packera dubia", "Packera anonyma"), each = 3L),
          determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", basal_leaf_count = 3L,
          PC1 = c(-2.0, -2.1, -1.9, 2.0, 2.1, 1.9), PC2 = c(0.1, -0.1, 0.0, 0.1, -0.1, 0.0),
          cda_predicted_taxon = rep(c("Packera dubia", "Packera anonyma"), each = 3L), cda_posterior_prob = 0.90,
          edaphic_best_fit_taxon = rep(c("Packera dubia", "Packera anonyma"), each = 3L), rf_predicted_prob = 0.88,
          delta_doy_resid = 2.0, is_active_anthesis = TRUE,
          geno_linkage_tier = NA_character_, geno_quality_pass = NA, geno_predicted_taxon = NA_character_, geno_posterior_prob = NA_real_,
          stringsAsFactors = FALSE
        )
        zp <- run_step07_genomic_propagation(df_zero)
        stopifnot(isFALSE(attr(zp, "genomic_propagation_summary")$active))
        for (pcol in c("genomic_propagation_tier", "tier0_calibrated_taxon", "tier0_calibrated_prob", "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi")) {
          stopifnot(all(is.na(zp[[pcol]])))
        }
        r_base <- apply_taxonomic_decision_matrix(df_zero)
        r_prop <- apply_taxonomic_decision_matrix(zp)
        for (ccol in c("probabilistic_consensus_score", "authority_weight", "taxonomic_treatment_subtype", "taxonomic_status_call", "recommended_determination", "synthesis_rationale")) {
          stopifnot(identical(r_base[[ccol]], r_prop[[ccol]]))
        }
        cat("STEP07_FOUR_TIER_PROPAGATION_ALL_OK\\n")
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_step07_genomic_prop.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step07", cwd=PROJECT_ROOT, skip_figures=True)
        self.assertEqual(proc.returncode, 0, f"Step 07 genomic propagation R suite failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP07_FOUR_TIER_PROPAGATION_ALL_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
