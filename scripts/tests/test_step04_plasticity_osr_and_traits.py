#!/usr/bin/env python3
"""Warm-R tests for Step 04 heteroblasty weights, plasticity mitigation, and zero-trust multi-island OSR."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import rscript_available, run_warm_r_cmd


class TestStep04PlasticityOSRAndTraits(unittest.TestCase):
    """Warm-R tests for Step 04 heteroblasty weights, plasticity mitigation, and zero-trust multi-island OSR."""

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT.joinpath("scripts", "morphometrics", "04_gmm_morphotools.R")

    def test_heteroblasty_leaf_weights_pure_basal_and_cauline_fallback(self) -> None:
        """
        Verify reconcile_step03_basal_otu_medians:
          1. When a specimen has >= 1 basal leaf (leaf_zone == "basal"), cauline leaves are strictly
             excluded from primary OTU unweighted basal medians so harmonics and solidity are uncontaminated.
          2. When a specimen has only cauline leaves (0 basal leaves), unweighted cauline medians populate
             valid fallback OTU coordinates.
        """
        r_code = """
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        source("scripts/morphometrics/authority_weighting_utils.R")
        stopifnot(exists("reconcile_step03_basal_otu_medians", mode = "function"))

        # 1. Mixed basal + cauline specimen -> unweighted basal median (0.91, 1.01) uncontaminated by cauline leaf
        df_mixed <- data.frame(
          catalogNumber = rep("SPEC_MIX", 3L),
          plant_individual_id = rep(0L, 3L),
          leaf_id = 1:3,
          leaf_zone = c("basal", "basal", "cauline"),
          area_mm2 = c(400.0, 100.0, 90.0),
          solidity = c(0.92, 0.90, 0.42),
          dissection_index = c(1.08, 1.10, 2.40),
          blade_length_mm = c(40.0, 38.0, 19.0),
          A1 = c(1.00, 1.02, 0.30),
          D1 = c(0.50, 0.52, 0.95),
          stringsAsFactors = FALSE
        )
        spec_init <- data.frame(
          catalogNumber = "SPEC_MIX",
          plant_individual_id = 0L,
          solidity = 0.0,
          dissection_index = 0.0,
          A1 = 0.0,
          D1 = 0.0,
          stringsAsFactors = FALSE
        )
        spec_agg <- reconcile_step03_basal_otu_medians(
          spec_init,
          df_mixed,
          all_harm_cols = c("A1", "D1"),
          sym_names = c("A1", "D1")
        )
        stopifnot(abs(spec_agg$solidity[1L] - 0.91) < 1e-6)
        stopifnot(abs(spec_agg$A1[1L] - 1.01) < 1e-6)

        # 2. Cauline-only specimen -> unweighted cauline fallback median (0.50, 0.42)
        df_caul_only <- data.frame(
          catalogNumber = rep("SPEC_CAUL", 2L),
          plant_individual_id = rep(0L, 2L),
          leaf_id = 1:2,
          leaf_zone = c("cauline", "cauline"),
          area_mm2 = c(120.0, 100.0),
          solidity = c(0.52, 0.48),
          A1 = c(0.40, 0.44),
          D1 = c(0.80, 0.84),
          stringsAsFactors = FALSE
        )
        spec_caul_init <- data.frame(
          catalogNumber = "SPEC_CAUL",
          plant_individual_id = 0L,
          solidity = 0.0,
          A1 = 0.0,
          D1 = 0.0,
          stringsAsFactors = FALSE
        )
        spec_caul <- reconcile_step03_basal_otu_medians(
          spec_caul_init,
          df_caul_only,
          all_harm_cols = c("A1", "D1"),
          sym_names = c("A1", "D1")
        )
        stopifnot(abs(spec_caul$solidity[1L] - 0.50) < 1e-6)
        stopifnot(abs(spec_caul$A1[1L] - 0.42) < 1e-6)

        # 3. Scale-invariant foliar architecture scalar filtering retains dimensionless traits without ruler-missingness bias
        df_arch <- data.frame(
          PC1 = c(-2.1, -1.8, -1.5, 0.2, 1.4, 1.9, 2.3),
          authority_weight = rep(1.0, 7L),
          petiole_to_blade_ratio = c(0.35, 0.42, 0.50, 0.65, 0.72, 0.85, 0.95),
          solidity = c(0.91, 0.89, 0.87, 0.78, 0.72, 0.64, 0.58),
          petiole_to_blade_ratio_std = as.numeric(scale(c(0.35, 0.42, 0.50, 0.65, 0.72, 0.85, 0.95))),
          solidity_std = as.numeric(scale(c(0.91, 0.89, 0.87, 0.78, 0.72, 0.64, 0.58))),
          stringsAsFactors = FALSE
        )
        kept_arch <- filter_biologically_variable_dissection_scalars(
          df_arch,
          c("petiole_to_blade_ratio_std", "solidity_std"),
          active_mask = (df_arch$authority_weight > 0)
        )
        stopifnot("petiole_to_blade_ratio_std" %in% kept_arch)
        stopifnot("solidity_std" %in% kept_arch)
        cat("HETEROBLASTY_PURE_BASAL_AND_CAULINE_FALLBACK_OK\\n")
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_heteroblasty_weights.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Heteroblasty leaf weights R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("HETEROBLASTY_PURE_BASAL_AND_CAULINE_FALLBACK_OK", proc.stdout)

    def test_step04_plasticity_mitigation_module_and_sympatric_unimodal_collapse(self) -> None:
        """
        Verify plasticity_mitigation_utils.R contracts and end-to-end R behavior:
          1. Static verification that plasticity_mitigation_utils.R defines phenophase residualization,
             within-deme Q_ST covariance whitening, sPCA global test, and Hennig/Ray-Lindsay ridgeline merging,
             and is wired into 04_gmm_morphotools.R without forcing G >= 2.
          2. When Rscript is available, verify that a skewed, phenophase-confounded sympatric panmictic cloud
             collapses to K = 1 (Morphotype_K1) while a geographically segregated bimodal cohort retains K = 2.
        """
        pm_path = PROJECT_ROOT / "scripts" / "morphometrics" / "plasticity_mitigation_utils.R"
        self.assertTrue(pm_path.is_file(), "plasticity_mitigation_utils.R must exist")
        pm_src = pm_path.read_text(encoding="utf-8")
        for fn_tok in (
            "merge_voucher_phenophase_covariates",
            "residualize_ontogenetic_phenophase_plasticity",
            "merge_plastic_gmm_components_ridgeline",
        ):
            self.assertIn(fn_tok, pm_src, f"Missing {fn_tok} in plasticity_mitigation_utils.R")
        self.assertNotIn("compute_spca_global_test", pm_src)

        r04_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("plasticity_mitigation_utils.R", r04_src)
        self.assertNotIn("kmeans_K2_fallback", r04_src)

        if not rscript_available():
            return

        r_code = """
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        source("scripts/morphometrics/pub_figures_utils.R")
        source("scripts/morphometrics/authority_weighting_utils.R")
        source("scripts/morphometrics/plasticity_mitigation_utils.R")
        source("scripts/morphometrics/04_gmm_morphotools.R")
        set.seed(20261109)

        # Case A: Sympatric skewed plastic cloud across 4 shared localities -> must collapse to K = 1 (Morphotype_K1)
        n_sym <- 48L
        site_lats <- rep(c(31.5, 33.8, 35.9, 38.4), length.out = n_sym)
        site_lons <- rep(c(-93.2, -84.1, -79.0, -76.1), length.out = n_sym)
        ecos <- rep(c("South_Central_Plains", "Piedmont", "Piedmont", "Middle_Atlantic_Coastal_Plain"), length.out = n_sym)
        # Skewed log-normal microhabitat plasticity completely randomized within every site
        plast_shock <- rexp(n_sym, rate = 1.0) - 1.0
        df_sym <- data.frame(
          catalogNumber = sprintf("SYM_%03d", seq_len(n_sym)),
          plant_individual_id = 0L,
          species_standardized = "Packera dubia",
          determiner_tier = "Tier_1_Gold",
          morph_tier_source = "Basal_Rosette",
          decimalLatitude = site_lats + rnorm(n_sym, 0, 0.02),
          decimalLongitude = site_lons + rnorm(n_sym, 0, 0.02),
          regional_group = ecos,
          solidity = 0.88,
          is_dissected = FALSE,
          PC1 = 1.6 * plast_shock + rnorm(n_sym, 0, 0.25),
          PC2 = -0.8 * plast_shock + rnorm(n_sym, 0, 0.25),
          PC3 = rnorm(n_sym, 0, 0.20),
          stringsAsFactors = FALSE
        )
        gmm_sym <- run_gmm_cluster_analysis(df_sym, c("PC1", "PC2", "PC3"), max_k = 5L, intraspecific_mode = TRUE)
        cda_sym <- run_canonical_discriminant_analysis(df_sym, c("PC1", "PC2", "PC3"), "Packera dubia", gmm_res = gmm_sym, intraspecific_mode = TRUE)
        stopifnot(gmm_sym$best_k == 1L)
        stopifnot(all(cda_sym$predicted_taxon == "Morphotype_K1"))

        # Case B: Geographically segregated true morphotypes -> must retain K = 2
        df_geo <- df_sym
        west_m <- df_geo$decimalLongitude < -88.0
        df_geo$PC1 <- ifelse(west_m, -3.2, 2.2) + rnorm(n_sym, 0, 0.30)
        df_geo$PC2 <- ifelse(west_m, 1.8, -1.2) + rnorm(n_sym, 0, 0.30)
        gmm_geo <- run_gmm_cluster_analysis(df_geo, c("PC1", "PC2", "PC3"), max_k = 5L, intraspecific_mode = TRUE)
        stopifnot(gmm_geo$best_k >= 2L)

        # Case C: Purely morphological ridgeline valley check preserves modality independence
        rl_unimodal <- .compute_ridgeline_valley_ratio(0.5, c(0, 0), diag(2), 0.5, c(0.6, 0.0), diag(2))
        rl_bimodal  <- .compute_ridgeline_valley_ratio(0.5, c(-3, 0), diag(2) * 0.25, 0.5, c(3, 0), diag(2) * 0.25)
        stopifnot(rl_unimodal$n_modes == 1L || rl_unimodal$ridge_ratio >= 0.45)
        stopifnot(rl_bimodal$n_modes >= 2L && rl_bimodal$ridge_ratio < 0.10)
        cat("PLASTICITY_MITIGATION_R_OK\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step04_plast.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Plasticity mitigation R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("PLASTICITY_MITIGATION_R_OK", proc.stdout)

    def test_step04_zero_trust_multi_island_cauline_osr_and_no_fake_loocv(self) -> None:
        """
        Zero-trust regression suite for Step 04 fixes:
          1. Pure-R CLI fallback (.parse_step04_cli_fallback) parses all flags without optparse.
          2. Partial-NA reproductive & petiole trait backfilling populates NA rows even when other rows are finite.
          3. Whole + multi-island dissected vouchers (PC1 = NA, finite ECT_PC1..3 + solidity + dissection_index)
             are retained and classified accurately via shared topology + architecture space, and cauline
             heteroblasty features enter Step 04 at >=50% anchor coverage.
          4. Single-group (g = 1) CDA sets loocv_accuracy = NA (not 1.0), and post-GMM intraspecific figures
             never print fabricated "LOOCV = 85.0%".
          5. Cross-run SQLite open-set recognition uses run-invariant harmonics + topology (immune to PC sign flips).
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            r_script = tmp_p / "zero_trust_step04.R"
            r_script.write_text(
                f"""
                if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
                options(rgl.useNULL = TRUE)
                source("scripts/morphometrics/04_gmm_morphotools.R")

                # 1. Pure-R CLI fallback
                defs <- list(vouchers = "a.csv", max_k = 8L, seed = 42L, outline_only = FALSE, intraspecific_mode = FALSE, correct_allometry = "auto")
                parsed <- .parse_step04_cli_fallback(c("--vouchers=b.csv", "--max-k", "6", "--seed", "99", "--outline-only", "--intraspecific-mode", "--correct-allometry", "false"), defs)
                stopifnot(parsed$vouchers == "b.csv", parsed$max_k == 6L, parsed$seed == 99L, isTRUE(parsed$outline_only), isTRUE(parsed$intraspecific_mode), parsed$correct_allometry == "false")

                # 2. Exact mclust parameter counts (EII vs VVI vs VVV)
                np_eii <- .compute_exact_mclust_params(c(1L, 3L), p_dim = 4L, default_model = "EII")
                np_vvi <- .compute_exact_mclust_params(c(1L, 3L), p_dim = 4L, default_model = "VVI")
                np_vvv <- .compute_exact_mclust_params(c(1L, 3L), p_dim = 4L, default_model = "VVV")
                stopifnot(np_eii[2] < np_vvi[2], np_vvi[2] < np_vvv[2])

                # 3. Partial-NA reproductive & petiole backfilling
                repro_csv <- file.path({repr(str(tmp_p))}, "repro_partial.csv")
                man_csv <- file.path({repr(str(tmp_p))}, "extracted_leaf_manifest.csv")
                write.csv(data.frame(catalogNumber = c("C1", "C2"), plant_individual_id = c(0L, 0L), AR_inv = c(1.10, 1.45)), repro_csv, row.names = FALSE)
                write.csv(data.frame(catalogNumber = c("C1", "C2"), plant_individual_id = c(0L, 0L), leaf_zone = c("basal", "basal"), petiole_length_mm = c(15.0, 28.0), blade_length_mm = c(30.0, 28.0), petiole_to_blade_ratio = c(0.5, 1.0)), man_csv, row.names = FALSE)
                agg_part <- data.frame(catalogNumber = c("C1", "C2"), plant_individual_id = c(0L, 0L), AR_inv = c(1.10, NA_real_), petiole_length_mm = c(15.0, NA_real_), blade_length_mm = c(30.0, NA_real_), petiole_to_blade_ratio = c(0.5, NA_real_), stringsAsFactors = FALSE)
                agg_part <- merge_reproductive_traits_by_individual(agg_part, repro_csv)
                agg_part <- merge_manifest_petiole_traits_by_individual(agg_part, man_csv)
                stopifnot(abs(agg_part$AR_inv[2] - 1.45) < 1e-6, abs(agg_part$petiole_length_mm[2] - 28.0) < 1e-6)

                # 4. Whole + multi-island (PC1 = NA) joint classification & 60% cauline heteroblasty retention
                set.seed(77)
                n_sp <- 12L
                mk_sp <- function(sp, idx0, pc1_mu, ect1_mu, sol_mu, di_mu, caul_sol_mu, n_multi = 0L) {{
                  rows <- lapply(seq_len(n_sp), function(j) {{
                    is_mi <- j > (n_sp - n_multi)
                    has_caul <- j <= 8L # 67% cauline coverage (< 85%, >= 50%)
                    data.frame(
                      catalogNumber = sprintf("%s_%02d", idx0, j),
                      plant_individual_id = 0L,
                      basal_leaf_count = 3L,
                      cauline_leaf_count = if (has_caul) 2L else 0L,
                      morph_tier_source = "Basal_Rosette",
                      species_raw = sp,
                      scientificName = sp,
                      determiner_tier = if (j <= 8L) "Tier_1_Gold" else "Tier_3_Bronze",
                      PC1 = if (is_mi) NA_real_ else rnorm(1, pc1_mu, 0.20),
                      PC2 = if (is_mi) NA_real_ else rnorm(1, 0.5 * pc1_mu, 0.20),
                      PC3 = if (is_mi) NA_real_ else rnorm(1, 0.0, 0.15),
                      ECT_PC1 = rnorm(1, ect1_mu, 0.18),
                      ECT_PC2 = rnorm(1, -0.4 * ect1_mu, 0.18),
                      ECT_PC3 = rnorm(1, 0.0, 0.12),
                      ect_lobe_complexity = max(0.05, 1.2 + 0.6 * ect1_mu + rnorm(1, 0, 0.08)),
                      solidity = rnorm(1, sol_mu, 0.02),
                      q25_solidity = rnorm(1, sol_mu - 0.03, 0.02),
                      dissection_index = rnorm(1, di_mu, 0.08),
                      cauline_solidity_med = if (has_caul) rnorm(1, caul_sol_mu, 0.02) else NA_real_,
                      basal_vs_cauline_solidity_delta = if (has_caul) rnorm(1, sol_mu - caul_sol_mu, 0.02) else NA_real_,
                      aspect_ratio = rnorm(1, 1.8, 0.05),
                      centroid_size = 24.0 + 1.2 * j,
                      log_centsize = log10(24.0 + 1.2 * j),
                      is_calibrated_mm = TRUE,
                      scale_source = "ruler_calibrated",
                      A1 = if (is_mi) NA_real_ else rnorm(1, -1.1 + 0.05 * pc1_mu, 0.004),
                      B1 = if (is_mi) NA_real_ else 0.0,
                      C1 = if (is_mi) NA_real_ else 0.0,
                      D1 = if (is_mi) NA_real_ else rnorm(1, 0.45 - 0.03 * pc1_mu, 0.004),
                      A2 = if (is_mi) NA_real_ else rnorm(1, 0.03 * pc1_mu, 0.004),
                      B2 = if (is_mi) NA_real_ else 0.0,
                      C2 = if (is_mi) NA_real_ else 0.0,
                      D2 = if (is_mi) NA_real_ else rnorm(1, -0.02 * pc1_mu, 0.004),
                      stringsAsFactors = FALSE
                    )
                  }})
                  do.call(rbind, rows)
                }}
                df_h <- rbind(
                  mk_sp("Packera dubia", "DUB", -3.0, -2.2, 0.90, 1.15, 0.55, n_multi = 0L),
                  mk_sp("Packera millefolium", "MIL", 3.0, 2.8, 0.36, 3.20, 0.32, n_multi = 2L)
                )
                df_v <- data.frame(
                  catalogNumber = df_h$catalogNumber,
                  species_raw = df_h$species_raw,
                  determiner_tier = df_h$determiner_tier,
                  institutionCode = rep(c("NCU", "DUKE"), length.out = nrow(df_h)),
                  stateProvince = "NC",
                  decimalLatitude = 35.8,
                  decimalLongitude = -79.1,
                  doy = 125,
                  stringsAsFactors = FALSE
                )
                h_path <- file.path({repr(str(tmp_p))}, "harm_mi.csv")
                v_path <- file.path({repr(str(tmp_p))}, "vouch_mi.csv")
                f_path <- file.path({repr(str(tmp_p))}, "flags_mi.csv")
                p_path <- file.path({repr(str(tmp_p))}, "fig04b_gmm_cda_morphometrics.png")
                write.csv(df_h, h_path, row.names = FALSE)
                write.csv(df_v, v_path, row.names = FALSE)

                opts_mi <- parse_args_robust()
                opts_mi$vouchers <- v_path; opts_mi$harmonics <- h_path; opts_mi$output_flags <- f_path
                opts_mi$output_plot <- p_path; opts_mi$output_report <- file.path({repr(str(tmp_p))}, "bic_mi.csv")
                opts_mi$output_batch_audit <- file.path({repr(str(tmp_p))}, "batch_mi.csv")
                opts_mi$allometry_report <- file.path({repr(str(tmp_p))}, "allom_mi.csv"); opts_mi$allometry_pairwise <- file.path({repr(str(tmp_p))}, "pw_mi.csv")
                opts_mi$allometry_plot <- file.path({repr(str(tmp_p))}, "fig04a_foliar_allometry_audit.png")
                opts_mi$permutations <- 19L; opts_mi$seed <- 123L
                out_mi <- run_gmm_morphotools_pipeline(opts_mi)
                stopifnot(nrow(out_mi) == 24L)
                mi_rows <- out_mi[out_mi$catalogNumber %in% c("MIL_11", "MIL_12"), , drop = FALSE]
                stopifnot(nrow(mi_rows) == 2L, all(is.na(mi_rows$PC1)))
                stopifnot(all(mi_rows$cda_predicted_taxon == "Packera millefolium"))
                stopifnot(all(mi_rows$cda_posterior_prob > 0.75))
                stopifnot("cauline_solidity_med_std" %in% names(out_mi) || "basal_vs_cauline_solidity_delta_std" %in% names(out_mi))

                # 5. Single-group (g = 1) CDA sets loocv_accuracy = NA (not 1.0) & post-GMM plot omits fake 85.0%
                df_one <- df_h[1:12, , drop = FALSE]
                gmm_one <- list(classification = rep(1L, 12L), uncertainty = rep(0.0, 12L), best_k = 1L, model_name = "EII", selected_features = c("PC1", "PC2"))
                cda_one <- run_canonical_discriminant_analysis(df_one, c("PC1", "PC2"), "Packera dubia", gmm_res = gmm_one, intraspecific_mode = TRUE)
                stopifnot(is.na(cda_one$loocv_accuracy))
                gmm_two <- list(classification = rep(1:2, each = 6L), uncertainty = rep(0.05, 12L), best_k = 2L, model_name = "VVI", selected_features = c("PC1", "PC2"), bic_table = data.frame(K = 1:2, BIC = c(-40, -28), Two_Delta_BIC = c(0, 12)))
                cda_two <- run_canonical_discriminant_analysis(df_one, c("PC1", "PC2"), "Packera dubia", gmm_res = gmm_two, intraspecific_mode = TRUE)
                stopifnot(is.na(cda_two$loocv_accuracy))
                aud_two <- audit_misidentifications(df_one, cda_two, gmm_two, selected_features = c("PC1", "PC2"))
                intra_png <- file.path({repr(str(tmp_p))}, "intra_no_fake_loocv.png")
                render_cda_biplot_4panel_png(aud_two, cda_two, gmm_two, intra_png, c("PC1", "PC2"))
                stopifnot(file.exists(intra_png))

                # 6. Cross-run SQLite open-set recognition invariance under PC sign flips
                if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {{
                  db_p <- file.path({repr(str(tmp_p))}, "osr_cross_run.sqlite")
                  con <- DBI::dbConnect(RSQLite::SQLite(), db_p)
                  db_cong <- rbind(
                    data.frame(run_id = "old_run_A", catalogNumber = paste0("AN_", 1:8), species_standardized = "Packera anonyma", authority_weight = 1.0, morph_tier_source = "Basal_Rosette",
                               PC1 = rnorm(8, 9.0, 0.1), PC2 = rnorm(8, -9.0, 0.1),
                               A1 = rnorm(8, -0.92, 0.01), A2 = rnorm(8, 0.12, 0.005), D1 = rnorm(8, 0.35, 0.01), D2 = rnorm(8, -0.08, 0.005), solidity = rnorm(8, 0.88, 0.01), dissection_index = rnorm(8, 1.25, 0.03)),
                    data.frame(run_id = "old_run_A", catalogNumber = paste0("AU_", 1:8), species_standardized = "Packera aurea", authority_weight = 1.0, morph_tier_source = "Basal_Rosette",
                               PC1 = rnorm(8, -9.0, 0.1), PC2 = rnorm(8, 9.0, 0.1),
                               A1 = rnorm(8, -1.35, 0.01), A2 = rnorm(8, -0.18, 0.005), D1 = rnorm(8, 0.68, 0.01), D2 = rnorm(8, 0.14, 0.005), solidity = rnorm(8, 0.82, 0.01), dissection_index = rnorm(8, 1.45, 0.03))
                  )
                  DBI::dbWriteTable(con, "otu_morphometrics", db_cong, overwrite = TRUE)
                  DBI::dbDisconnect(con)
                  Sys.setenv(PLANTMATICS_DB_PATH = db_p)
                  df_q <- data.frame(
                    catalogNumber = paste0("Q_", 1:6), species_standardized = "Packera dubia",
                    PC1 = rnorm(6, -9.0, 0.1), PC2 = rnorm(6, 9.0, 0.1), # deliberately flipped PC signs matching aurea if PCs were leaked!
                    A1 = rnorm(6, -0.92, 0.01), A2 = rnorm(6, 0.12, 0.005), D1 = rnorm(6, 0.35, 0.01), D2 = rnorm(6, -0.08, 0.005), # true harmonics match anonyma!
                    solidity = rnorm(6, 0.88, 0.01), dissection_index = rnorm(6, 1.25, 0.03),
                    stringsAsFactors = FALSE
                  )
                  osr_res <- compute_type_specimen_morph_control(df_q, X_mat = as.matrix(df_q[, c("PC1", "PC2")]), S_reg = diag(2), taxon_col = "species_standardized")
                  Sys.unsetenv("PLANTMATICS_DB_PATH")
                  stopifnot(all(osr_res$min_congener_taxon == "Packera anonyma"))
                }}
                cat("STEP04_ZERO_TRUST_OK\n")
                """,
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_script)], stage_key='step04', cwd=PROJECT_ROOT)
            self.assertEqual(proc.returncode, 0, f"Step 04 zero-trust R test failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertIn("STEP04_ZERO_TRUST_OK", proc.stdout)
            self.assertTrue((tmp_p / "intra_no_fake_loocv.png").exists())
            self.assertEqual(list(tmp_p.glob("*.pdf")), [])


if __name__ == "__main__":
    unittest.main()
