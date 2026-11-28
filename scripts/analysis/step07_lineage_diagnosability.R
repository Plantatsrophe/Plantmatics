#!/usr/bin/env Rscript
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# Script: scripts/analysis/step07_lineage_diagnosability.R
# Stage 06 Supervised Candidate-Lineage Diagnosability, Pairwise Niche/Barrier
# Equivalence, Dichotomous Key Character Table, and Basionym Adjudication Engine.

suppressPackageStartupMessages({
  for (pkg in c("MASS", "ggplot2", "optparse")) {
    if (requireNamespace(pkg, quietly = TRUE)) library(pkg, character.only = TRUE)
  }
})

.s06_dirs <- unique(c(
  "scripts/analysis",
  dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
            unlist(lapply(sys.frames(), function(fr) fr$ofile))))
))
.src_s06d <- function(rel) {
  for (p in unique(c(file.path("scripts", rel), file.path(.s06_dirs, basename(rel)), file.path(.s06_dirs, "..", rel)))) {
    if (!is.na(p) && nzchar(p) && file.exists(p)) { source(p); return(invisible(TRUE)) }
  }
  invisible(FALSE)
}
.src_s06d("figures/theme_packera_pub.R")
.src_s06d("core/db_utils.R")
.src_s06d("analysis/step07_diagnosability_key_utils.R")

.write_csv_atomic_s06d <- function(df_obj, path_out) {
  dir.create(dirname(path_out), recursive = TRUE, showWarnings = FALSE)
  tmp_p <- paste0(path_out, ".", Sys.getpid(), ".tmp")
  utils::write.csv(df_obj, file = tmp_p, row.names = FALSE, na = "")
  file.rename(tmp_p, path_out)
}

parse_stage06_diag_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (requireNamespace("optparse", quietly = TRUE)) {
    mo <- optparse::make_option
    ol <- list(
      mo("--source-run", type = "character", default = ""),
      mo("--run-dir", type = "character", default = ""),
      mo("--run-id", type = "character", default = ""),
      mo("--focal-taxon", type = "character", default = "Packera dubia"),
      mo("--demo", action = "store_true", default = FALSE)
    )
    opts <- optparse::parse_args(optparse::OptionParser(option_list = ol), args = args)
    names(opts) <- gsub("-", "_", names(opts))
    return(opts)
  }
  stop("optparse package is required for step07_lineage_diagnosability.R")
}

resolve_source_run_dir <- function(src_tok) {
  clean <- trimws(as.character(src_tok))
  cands <- c(clean, file.path("runs", clean), sub("^runs/", "", clean))
  for (c_p in cands) {
    if (nzchar(c_p) && dir.exists(c_p) && file.exists(file.path(c_p, "data", "tables", "triage_queue.csv"))) {
      return(normalizePath(c_p, mustWork = TRUE))
    }
  }
  stop(sprintf("Stage 06 could not locate Stage 03 source run tables for --source-run='%s'", clean))
}

adjudicate_lineage_ranks <- function(df, lin_res, oof_res, pair_df, key_df, focal_taxon) {
  valid_lin <- lin_res$valid_lineages
  ss_lin <- lin_res$ss_lineage
  if (length(valid_lin) < 2L) {
    return(data.frame(
      candidate_lineage = if (length(valid_lin) == 1L) valid_lin[1L] else gsub("[[:space:]]+", "_", focal_taxon),
      n_plants = nrow(df), n_gold_silver = nrow(df),
      global_oof_morph_accuracy = NA_real_, oof_morph_recall = NA_real_,
      reproductive_oof_morph_accuracy = NA_real_, reproductive_oof_recall = NA_real_, passes_75_rule = FALSE,
      has_non_overlapping_key_trait = FALSE, median_soil_sand_pct = NA_real_, median_soil_ph = NA_real_,
      median_lithic_outcrop_index = NA_real_, min_d_edaph_vs_others = 1.0, min_gap_km_vs_others = 0.0,
      median_coalescent_gdi = 0.0, available_basionym = focal_taxon,
      recommended_rank = "Uniform_Single_Species", stringsAsFactors = FALSE
    ))
  }
  tier_v <- if ("determiner_tier" %in% names(df)) .nz_s06d(df$determiner_tier) else rep("", nrow(df))
  sand_v <- if ("fine_soil_sand" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_sand)) else rep(NA_real_, nrow(df))
  ph_v <- if ("fine_soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_ph)) else rep(NA_real_, nrow(df))
  lith_v <- if ("lithic_outcrop_index" %in% names(df)) suppressWarnings(as.numeric(df$lithic_outcrop_index)) else rep(NA_real_, nrow(df))
  gdi_v <- if ("coalescent_gdi" %in% names(df)) suppressWarnings(as.numeric(df$coalescent_gdi)) else rep(NA_real_, nrow(df))
  is_type_v <- if ("is_type_specimen" %in% names(df)) (!is.na(df$is_type_specimen) & as.logical(df$is_type_specimen)) else rep(FALSE, nrow(df))
  typ_bas_v <- if ("typified_basionym" %in% names(df)) .nz_s06d(df$typified_basionym) else rep("", nrow(df))
  nt_bas_v <- if ("nearest_type_niche_basionym" %in% names(df)) .nz_s06d(df$nearest_type_niche_basionym) else rep("", nrow(df))
  nt_bz_v <- if ("nearest_type_basionym_edaphic_z" %in% names(df)) suppressWarnings(as.numeric(df$nearest_type_basionym_edaphic_z)) else rep(NA_real_, nrow(df))

  any_iqr_key <- nrow(key_df) > 0L && any(isTRUE(key_df$has_non_overlapping_iqr) | key_df$has_non_overlapping_iqr == TRUE)
  rows_out <- list()
  for (g in valid_lin) {
    idx <- which(!is.na(lin_res$lineage) & lin_res$lineage == g)
    n_gs <- sum(tier_v[idx] %in% c("Tier_0_Type", "Tier_1_Gold", "Tier_2_Silver") | is_type_v[idx])
    rec_g <- if (g %in% names(oof_res$recall_by_lin)) as.numeric(oof_res$recall_by_lin[g]) else NA_real_
    rep_rec_g <- if (!is.null(oof_res$repro_recall_by_lin) && (g %in% names(oof_res$repro_recall_by_lin))) as.numeric(oof_res$repro_recall_by_lin[g]) else NA_real_
    rep_acc_all <- if (!is.null(oof_res$repro_oof_acc)) as.numeric(oof_res$repro_oof_acc) else NA_real_
    pass_75 <- (is.finite(rec_g) && (rec_g >= 0.75)) || (is.finite(rep_rec_g) && (rep_rec_g >= 0.75) && is.finite(rep_acc_all) && (rep_acc_all >= 0.70))

    pw_sub <- if (nrow(pair_df) > 0L) pair_df[pair_df$lineage_a == g | pair_df$lineage_b == g, , drop = FALSE] else data.frame()
    min_de <- if (nrow(pw_sub) > 0L && any(is.finite(pw_sub$d_edaph))) min(pw_sub$d_edaph, na.rm = TRUE) else NA_real_
    min_gap <- if (nrow(pw_sub) > 0L && any(is.finite(pw_sub$min_interlineage_gap_km))) min(pw_sub$min_interlineage_gap_km, na.rm = TRUE) else 0.0
    med_gdi <- if (any(is.finite(gdi_v[idx]))) stats::median(gdi_v[idx], na.rm = TRUE) else 0.0

    # Resolve available basionym
    bas_lbl <- ""
    if (g == ss_lin) {
      bas_lbl <- paste0(focal_taxon, " (sensu stricto locus classicus)")
    } else {
      t_hits <- unique(typ_bas_v[idx][is_type_v[idx] & nzchar(typ_bas_v[idx])])
      if (length(t_hits) > 0L) {
        bas_lbl <- t_hits[1L]
      } else if (!grepl("Disjunct", g, ignore.case = TRUE)) {
        nb_ok <- nt_bas_v[idx][nzchar(nt_bas_v[idx]) & is.finite(nt_bz_v[idx]) & nt_bz_v[idx] <= 2.5]
        if (length(nb_ok) >= 3L) bas_lbl <- names(sort(table(nb_ok), decreasing = TRUE))[1L]
      }
      if (!nzchar(bas_lbl)) bas_lbl <- "None_Matched (Untypified_Lineage)"
    }

    rank_call <- if (g == ss_lin) {
      "Sensu_Stricto_Typical_Lineage"
    } else if (is.finite(med_gdi) && med_gdi >= 0.70 && pass_75) {
      "Distinct_Species_Candidate"
    } else if (pass_75 || any_iqr_key) {
      "Formal_Subspecies_or_Variety"
    } else if (is.finite(min_gap) && min_gap >= 200.0) {
      "Allopatric_Disjunct_Lineage (Morphologically_Overlapping)"
    } else if (!is.finite(min_de) || min_de < 0.65) {
      "Edaphic_Ecotype (Informal_Ecological_Race)"
    } else {
      "Continuous_Intraspecific_Variant"
    }

    rows_out[[length(rows_out) + 1L]] <- data.frame(
      candidate_lineage = g, n_plants = length(idx), n_gold_silver = n_gs,
      global_oof_morph_accuracy = round(oof_res$global_oof_acc, 4),
      oof_morph_recall = round(rec_g, 4),
      reproductive_oof_morph_accuracy = if (is.finite(rep_acc_all)) round(rep_acc_all, 4) else NA_real_,
      reproductive_oof_recall = if (is.finite(rep_rec_g)) round(rep_rec_g, 4) else NA_real_,
      passes_75_rule = pass_75,
      has_non_overlapping_key_trait = any_iqr_key,
      median_soil_sand_pct = if (any(is.finite(sand_v[idx]))) round(stats::median(sand_v[idx], na.rm = TRUE), 1) else NA_real_,
      median_soil_ph = if (any(is.finite(ph_v[idx]))) round(stats::median(ph_v[idx], na.rm = TRUE), 2) else NA_real_,
      median_lithic_outcrop_index = if (any(is.finite(lith_v[idx]))) round(stats::median(lith_v[idx], na.rm = TRUE), 4) else NA_real_,
      min_d_edaph_vs_others = if (is.finite(min_de)) round(min_de, 3) else NA_real_,
      min_gap_km_vs_others = if (is.finite(min_gap)) round(min_gap, 1) else 0.0,
      median_coalescent_gdi = round(med_gdi, 4),
      available_basionym = bas_lbl, recommended_rank = rank_call,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows_out)
}

main_stage06_diagnosability <- function() {
  opts <- parse_stage06_diag_args()
  src_dir <- resolve_source_run_dir(opts$source_run)
  out_run_dir <- opts$run_dir
  if (!nzchar(out_run_dir)) stop("--run-dir is required")

  tq_src <- utils::read.csv(file.path(src_dir, "data", "tables", "triage_queue.csv"), stringsAsFactors = FALSE)
  mc_path <- file.path(src_dir, "data", "tables", "multimodal_conflict_flags.csv")
  mf_path <- file.path(src_dir, "data", "tables", "morphometric_flags.csv")
  mc_df <- if (file.exists(mc_path)) utils::read.csv(mc_path, stringsAsFactors = FALSE) else data.frame()
  mf_df <- if (file.exists(mf_path)) utils::read.csv(mf_path, stringsAsFactors = FALSE) else data.frame()

  key_join <- function(d) paste(as.character(d$catalogNumber), ifelse(is.na(d$plant_individual_id), "0", as.character(d$plant_individual_id)), sep = "||")
  merged <- tq_src
  merged$._key <- key_join(merged)
  for (aux in list(mc_df, mf_df)) {
    if (nrow(aux) > 0L && all(c("catalogNumber", "plant_individual_id") %in% names(aux))) {
      aux$._key <- key_join(aux)
      aux <- aux[!duplicated(aux$._key), , drop = FALSE]
      add_cols <- setdiff(names(aux), names(merged))
      if (length(add_cols) > 0L) {
        m_idx <- match(merged$._key, aux$._key)
        for (cl in add_cols) merged[[cl]] <- aux[[cl]][m_idx]
      }
    }
  }

  s4_sibling_mc <- file.path(dirname(src_dir), "04_spatial_thin_field_audit", "data", "tables", "multimodal_conflict_flags.csv")
  merged <- hydrate_ecotypes_from_stage4(merged, s4_sibling_mc)
  lin_res <- assign_stage06_candidate_lineages(merged, min_n = 5L)
  morph_cols <- select_stage06_morph_features(merged)
  oof_res <- run_candidate_lineage_oof_cda(merged, lin_res$lineage, morph_cols)
  pair_df <- if (isTRUE(oof_res$available)) augment_pairwise_niche_and_pheno(merged, lin_res$lineage, oof_res$pair_df) else data.frame()
  key_df <- build_dichotomous_key_table(merged, lin_res$lineage)
  diag_df <- adjudicate_lineage_ranks(merged, lin_res, oof_res, pair_df, key_df, opts$focal_taxon)

  rank_map <- stats::setNames(diag_df$recommended_rank, diag_df$candidate_lineage)
  bas_map <- stats::setNames(diag_df$available_basionym, diag_df$candidate_lineage)
  eff_lin <- ifelse(!is.na(lin_res$lineage), lin_res$lineage,
                    ifelse(!is.na(oof_res$oof_pred), paste0("Projected_", oof_res$oof_pred), "Unassigned"))
  base_lin <- sub("^Projected_", "", eff_lin)

  tq_out <- tq_src
  tq_out$candidate_lineage <- eff_lin
  tq_out$lineage_oof_morph_prob <- round(oof_res$oof_prob, 4)
  tq_out$lineage_recommended_rank <- ifelse(base_lin %in% names(rank_map), rank_map[base_lin], "")
  tq_out$lineage_available_basionym <- ifelse(base_lin %in% names(bas_map), bas_map[base_lin], "")

  out_tables <- file.path(out_run_dir, "data", "tables")
  out_reports <- file.path(out_run_dir, "outputs", "reports")
  out_figs <- file.path(out_run_dir, "outputs", "figures")
  .write_csv_atomic_s06d(tq_out, file.path(out_tables, "triage_queue.csv"))
  .write_csv_atomic_s06d(diag_df, file.path(out_reports, "candidate_lineage_diagnosability.csv"))
  .write_csv_atomic_s06d(pair_df, file.path(out_reports, "candidate_lineage_pairwise_tests.csv"))
  .write_csv_atomic_s06d(key_df, file.path(out_reports, "dichotomous_key_character_table.csv"))
  .write_csv_atomic_s06d(diag_df[, c("candidate_lineage", "n_plants", "available_basionym", "recommended_rank"), drop = FALSE],
                         file.path(out_reports, "candidate_lineage_basionym_audit.csv"))

  # Mirror Stage 3 summary tables and record Stage 6 pairwise D_edaph + OOF morph accuracy for cross-stage table
  for (fn_cp in c("gmm_bayes_factors_summary.csv", "gmm_bic_model_selection.csv")) {
    src_f <- file.path(src_dir, "outputs", "reports", fn_cp)
    if (file.exists(src_f)) file.copy(src_f, file.path(out_reports, fn_cp), overwrite = TRUE)
  }
  min_pw_de <- if (nrow(pair_df) > 0L && any(is.finite(pair_df$d_edaph))) min(pair_df$d_edaph, na.rm = TRUE) else NA_real_
  .write_csv_atomic_s06d(
    data.frame(Taxon1 = opts$focal_taxon, Taxon2 = opts$focal_taxon,
               Spatial_Block_CV_Acc = round(oof_res$global_oof_acc, 4), D_edaph = min_pw_de, stringsAsFactors = FALSE),
    file.path(out_reports, "multimodal_spatial_rf_summary.csv")
  )

  slug <- gsub("[^A-Za-z0-9]+", "_", trimws(opts$focal_taxon))
  render_lineage_diagnosability_png(merged, lin_res$lineage, oof_res, key_df, opts$focal_taxon,
                                    file.path(out_figs, sprintf("Figure_Lineage_Diagnosability_%s.png", slug)))
  message(sprintf("[Stage 06 Diagnosability] Evaluated %d candidate lineages (global OOF morph acc = %s).",
                  lin_res$n_lineages, ifelse(is.finite(oof_res$global_oof_acc), sprintf("%.4f", oof_res$global_oof_acc), "NA")))
}

if (sys.nframe() == 0L) {
  main_stage06_diagnosability()
}

