if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/analysis/step07_genomic_report_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Step 07 genomic field triage CSV export, cytology/nQuire discordance lines,
#   and Markdown summary report synthesis for the Step 06 genomic stream.
# ==============================================================================

export_genomic_field_triage_csv <- function(opts, triage_df) {
  t_dir <- .resolve_step07_genomic_tables_dir(opts)
  r_dir <- .resolve_step07_genomic_reports_dir(opts)
  ga_path <- file.path(t_dir, "genomic_assignments.csv")
  if (!file.exists(ga_path)) ga_path <- file.path(r_dir, "genomic_assignments.csv")
  if (!file.exists(ga_path)) return(invisible(NULL))
  ga_df <- tryCatch(utils::read.csv(ga_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(ga_df) || nrow(ga_df) == 0L) return(invisible(NULL))

  ga_df$catalogNumber <- .coalesce_chr_cols_r(ga_df, c("catalogNumber", "linked_catalog_number"), default = NA_character_)
  n_g <- nrow(ga_df)
  ga_df$voucher_recommended_determination <- NA_character_
  ga_df$voucher_treatment_subtype <- NA_character_
  ga_df$voucher_triage_priority <- NA_character_
  ga_df$genomic_field_status <- rep("Genomic_Reference_Sample", n_g)

  if (!is.null(triage_df) && nrow(triage_df) > 0L && "catalogNumber" %in% names(triage_df)) {
    tri_u <- triage_df[!duplicated(triage_df$catalogNumber), , drop = FALSE]
    m <- match(ga_df$catalogNumber, trimws(as.character(tri_u$catalogNumber)))
    hit <- !is.na(m)
    if (any(hit)) {
      ga_df$voucher_recommended_determination[hit] <- as.character(tri_u$recommended_determination[m[hit]])
      ga_df$voucher_treatment_subtype[hit] <- as.character(tri_u$taxonomic_treatment_subtype[m[hit]])
      ga_df$voucher_triage_priority[hit] <- as.character(tri_u$synthesis_triage_priority[m[hit]])
    }
  }

  qc_ok <- .safe_bool_vec_r(ga_df$geno_quality_pass)
  is_q <- .safe_bool_vec_r(ga_df$is_putative_novel_or_hybrid_query)
  hyb_cls <- if ("hybrid_generation_class" %in% names(ga_df)) as.character(ga_df$hybrid_generation_class) else rep(NA_character_, n_g)
  is_div <- !is.na(hyb_cls) & (hyb_cls == "Divergent_NonHybrid_Lineage")
  is_adm <- if ("is_geno_admixed" %in% names(ga_df)) .safe_bool_vec_r(ga_df$is_geno_admixed) & !is_div else rep(FALSE, n_g)
  sp_std <- as.character(ga_df$species_standardized)
  g_pred <- .align_packera_orthography(as.character(ga_df$geno_predicted_taxon), sp_std)
  ga_df$geno_predicted_taxon <- g_pred
  g_prob <- suppressWarnings(as.numeric(ga_df$geno_posterior_prob))

  ga_df$genomic_field_status[!qc_ok] <- "Failed_Extraction_QC_Dropout"
  ga_df$genomic_field_status[qc_ok & is_q & is_div] <- "Divergent_NonHybrid_Lineage_Query"
  ga_df$genomic_field_status[qc_ok & is_q & !is_div] <- "Putative_Novel_Or_Hybrid_Query"
  ga_df$genomic_field_status[qc_ok & !is_q & is_adm] <- "Genomic_Admixed_Or_Introgressant"
  disc_mask <- qc_ok & !is_q & !is_adm & !is.na(sp_std) & nzchar(sp_std) & !is.na(g_pred) & nzchar(g_pred) & (sp_std != g_pred) & is.finite(g_prob) & (g_prob >= 0.70)
  ga_df$genomic_field_status[disc_mask] <- "Genomic_Discordant_Reassign_Candidate"
  conc_mask <- qc_ok & !is_q & !is_adm & !is.na(sp_std) & nzchar(sp_std) & !is.na(g_pred) & nzchar(g_pred) & (sp_std == g_pred)
  ga_df$genomic_field_status[conc_mask] <- "Genomic_Concordant_Anchor"

  out_p <- file.path(t_dir, "genomic_field_triage.csv")
  tmp_p <- paste0(out_p, ".tmp")
  readr::write_csv(ga_df, tmp_p, na = "")
  file.rename(tmp_p, out_p)
  invisible(out_p)
}

.build_cytology_summary_lines_r <- function(df) {
  out <- character(0)
  if (!is.data.frame(df) || nrow(df) == 0L) return(out)
  cy_lines <- attr(df, "cytology_markdown_lines")
  if (length(cy_lines) > 0L) {
    out <- c(out, cy_lines)
  } else if ("cytology_linkage_tier" %in% names(df)) {
    ct_v <- as.character(df$cytology_linkage_tier)
    n_c1 <- sum(!is.na(ct_v) & ct_v == "C1_Direct_Cytovoucher"); n_c2 <- sum(!is.na(ct_v) & ct_v == "C2_Exsiccatae_Cytoduplicate")
    mc_v <- if ("is_multi_cytotype_taxon" %in% names(df)) suppressWarnings(as.integer(df$is_multi_cytotype_taxon)) else rep(NA_integer_, nrow(df))
    sp_v <- if ("species_standardized" %in% names(df)) trimws(as.character(df$species_standardized)) else rep("", nrow(df))
    n_mc <- length(unique(sp_v[!is.na(mc_v) & mc_v == 1L & nzchar(sp_v)]))
    if (any(!is.na(ct_v)) || n_mc > 0L) out <- c(out, sprintf("- **Step 07 Cytology Coverage:** %d `C1` direct cytovouchers, %d `C2` exsiccatae cytoduplicates across %d multi-cytotype taxa.", n_c1, n_c2, n_mc))
  }
  if ("inferred_cytotype" %in% names(df)) {
    nq_pl <- suppressWarnings(as.integer(sub("x$", "", tolower(trimws(as.character(df$inferred_cytotype))))))
    nq_pl[is.na(nq_pl) | nq_pl < 2L | nq_pl > 16L] <- NA_integer_
    qc_v <- if ("geno_quality_pass" %in% names(df)) .safe_bool_vec_r(df$geno_quality_pass) else rep(TRUE, nrow(df))
    v_pl <- if ("voucher_ploidy_level" %in% names(df)) suppressWarnings(as.integer(df$voucher_ploidy_level)) else rep(NA_integer_, nrow(df))
    rep_v <- if ("taxon_reported_cytotypes" %in% names(df)) as.character(df$taxon_reported_cytotypes) else rep(NA_character_, nrow(df))
    sp_v <- if ("species_standardized" %in% names(df)) trimws(as.character(df$species_standardized)) else rep("", nrow(df))
    cat_v <- if ("catalogNumber" %in% names(df)) trimws(as.character(df$catalogNumber)) else as.character(seq_len(nrow(df)))
    disc_lbls <- character(0); n_chk <- 0L
    for (r_i in which(!is.na(nq_pl) & qc_v)) {
      lit_pl <- suppressWarnings(as.integer(sub("x$", "", unlist(strsplit(tolower(as.character(rep_v[r_i])), "[,;/ ]+"))))); lit_pl <- lit_pl[!is.na(lit_pl)]
      if (!is.na(v_pl[r_i]) || length(lit_pl) > 0L) n_chk <- n_chk + 1L
      if ((!is.na(v_pl[r_i]) && (nq_pl[r_i] != v_pl[r_i])) || (length(lit_pl) > 0L && !(nq_pl[r_i] %in% lit_pl))) {
        exp_s <- if (!is.na(v_pl[r_i])) paste0(v_pl[r_i], "x") else as.character(rep_v[r_i])
        disc_lbls <- c(disc_lbls, sprintf("%s %s (nQuire %dx vs %s)", sp_v[r_i], cat_v[r_i], nq_pl[r_i], exp_s))
      }
    }
    if (length(disc_lbls) > 0L) {
      out <- c(out, sprintf("- **Cytology vs. nQuire Discordance Audit:** %d discordant sample(s) out of %d checked (%s).", length(disc_lbls), n_chk, paste(utils::head(unique(disc_lbls), 5L), collapse = "; ")))
    } else if (n_chk > 0L) {
      out <- c(out, sprintf("- **Cytology vs. nQuire Discordance Audit:** 0 discordances across %d cytologically cross-checked sample(s).", n_chk))
    }
  }
  out
}

build_genomic_markdown_summary <- function(opts, df) {
  t_dir <- .resolve_step07_genomic_tables_dir(opts)
  r_dir <- .resolve_step07_genomic_reports_dir(opts)
  resolve_csv <- function(fname) {
    cands <- unique(c(
      r_dir, t_dir,
      if (!is.null(opts$genomic_dir) && nzchar(opts$genomic_dir)) opts$genomic_dir else character(0),
      if (!is.null(opts$out_dir) && nzchar(opts$out_dir)) opts$out_dir else character(0),
      file.path(dirname(r_dir), "genomics"), file.path(dirname(t_dir), "genomics"),
      file.path(dirname(dirname(t_dir)), "outputs", "genomics")
    ))
    wd <- normalizePath(getwd(), mustWork = FALSE)
    is_local <- nzchar(t_dir) && (!grepl("^/", t_dir) || startsWith(normalizePath(t_dir, mustWork = FALSE), wd))
    if (is_local && dir.exists("outputs/genomics")) cands <- c(cands, "outputs/genomics")
    for (d in cands) if (nzchar(d) && file.exists(p <- file.path(d, fname))) return(p)
    ""
  }
  read_csv_safe <- function(p) if (nzchar(p)) tryCatch(utils::read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  ga_path <- resolve_csv("genomic_assignments.csv")
  if (!nzchar(ga_path) && !any(has_active_genomic_row_mask(df))) return(.build_cytology_summary_lines_r(df))

  ga_df <- read_csv_safe(ga_path)
  gm_df <- read_csv_safe(resolve_csv("genomic_model_metrics.csv"))
  gc_df <- read_csv_safe(resolve_csv("genomic_multimodal_congruence_summary.csv"))
  pbc_df <- read_csv_safe(resolve_csv("public_backbone_concordance.csv"))
  pop_df <- read_csv_safe(resolve_csv("public_organellar_placement.csv"))
  tc_df <- read_csv_safe(resolve_csv("taxon_complexes.csv"))
  gdi_df <- attr(df, "genomic_gdi_delimitation"); if (is.null(gdi_df)) gdi_df <- read_csv_safe(resolve_csv("coalescent_gdi_delimitation.csv"))
  fb_df <- attr(df, "genomic_fbranch_summary"); if (is.null(fb_df)) fb_df <- read_csv_safe(resolve_csv("introgressive_fbranch_summary.csv"))
  land_df <- attr(df, "genomic_landscape_summary"); if (is.null(land_df)) land_df <- read_csv_safe(resolve_csv("landscape_mmrr_feems_summary.csv"))

  lines <- c("", "## 5. Independent Step 06 Genomic & Phylogenomic Evidence Stream")
  if (!is.null(gm_df) && nrow(gm_df) > 0L) {
    lines <- c(lines, sprintf(
      "- **Cohort QC, HDplot Paralog Filter & Shrinkage DAPC (LORO-CV):** %d/%d sequenced libraries passed QC (`hq_reads >= 50000`, `missing_genotype_rate <= 0.80`; %d reference anchors, %d novel/hybrid queries, %d HDplot-clean loci retained; LORO-CV accuracy = %.3f, Cohen's kappa = %.3f, mean `M`/`P` replicate concordance = %.4f).",
      as.integer(gm_df$n_samples_pass[1L]), as.integer(gm_df$n_samples_total[1L]),
      as.integer(gm_df$n_reference_anchors[1L]), as.integer(gm_df$n_novel_queries[1L]),
      as.integer(gm_df$n_loci_retained[1L]), suppressWarnings(as.numeric(gm_df$loro_cv_accuracy[1L])),
      suppressWarnings(as.numeric(gm_df$loro_cv_kappa[1L])), suppressWarnings(as.numeric(gm_df$mean_replicate_concordance[1L]))
    ))
  }
  if (!is.null(ga_df) && nrow(ga_df) > 0L) {
    lt_col <- .coalesce_chr_cols_r(ga_df, c("geno_linkage_tier", "linkage_tier"), default = "")
    n_g12 <- sum(lt_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy"), na.rm = TRUE)
    n_g3 <- sum(lt_col == "G3_Georeferenced_Unvouchered", na.rm = TRUE)
    n_g4 <- sum(lt_col == "G4_Ungeoreferenced_Tip", na.rm = TRUE)
    lines <- c(lines, sprintf(
      "- **Voucher Linkage Tiers:** %d voucher-linked (`G1`/`G2`), %d georeferenced field samples (`G3`), %d ungeoreferenced sequence tips (`G4`). Field and unvouchered query dispositions exported to `genomic_field_triage.csv`.",
      n_g12, n_g3, n_g4
    ))
    if ("inferred_cytotype" %in% names(ga_df)) {
      ct <- as.character(ga_df$inferred_cytotype[!is.na(ga_df$inferred_cytotype) & nzchar(ga_df$inferred_cytotype)])
      if (length(ct) > 0L) {
        ct_tab <- table(ct)
        ct_str <- paste(sprintf("%s=%d", names(ct_tab), as.integer(ct_tab)), collapse = ", ")
        hc_str <- if ("hybrid_generation_class" %in% names(ga_df)) {
          hc <- as.character(ga_df$hybrid_generation_class[!is.na(ga_df$hybrid_generation_class) & nzchar(ga_df$hybrid_generation_class)])
          if (length(hc) > 0L) sprintf("; `(h, Q12)` classes: %s", paste(sprintf("%s=%d", names(table(hc)), as.integer(table(hc))), collapse = ", ")) else ""
        } else ""
        lines <- c(lines, sprintf("- **Allele-Balance Cytotype & Hybrid Architecture:** Inferred ploidies (%s)%s.", ct_str, hc_str))
      }
    }
  }
  if (!is.null(gdi_df) && nrow(gdi_df) > 0L && "coalescent_gdi" %in% names(gdi_df)) {
    gdv <- suppressWarnings(as.numeric(gdi_df$coalescent_gdi)); gdv <- gdv[is.finite(gdv)]
    if (length(gdv) > 0L) {
      lines <- c(lines, sprintf(
        "- **Coalescent Delimitation (`SVDquartets` / `BPP` `MSci` / `gdi`):** Across %d evaluated lineages, median `gdi` = %.3f (range %.3f to %.3f; Jackson et al. 2017: `<0.20` conspecific oversplit block, `>0.70` distinct lineage).",
        nrow(gdi_df), stats::median(gdv), min(gdv), max(gdv)
      ))
    }
  }
  f_col <- if (!is.null(fb_df) && "f_branch_fb" %in% names(fb_df)) "f_branch_fb" else "f_branch"
  if (!is.null(fb_df) && nrow(fb_df) > 0L && (f_col %in% names(fb_df))) {
    f_vals <- suppressWarnings(as.numeric(fb_df[[f_col]]))
    fmax <- if (any(is.finite(f_vals))) max(f_vals, na.rm = TRUE) else NA_real_
    if (is.finite(fmax)) {
      top_fb <- fb_df[which.max(f_vals), , drop = FALSE]
      donor_s <- if ("donor_taxon" %in% names(top_fb)) as.character(top_fb$donor_taxon[1L]) else "NA"
      pol_s <- if ("polarization_source" %in% names(top_fb) && !is.na(top_fb$polarization_source[1L]) && nzchar(as.character(top_fb$polarization_source[1L]))) {
        sprintf("; polarization `%s`", as.character(top_fb$polarization_source[1L]))
      } else ""
      lines <- c(lines, sprintf("- **Tree-Guided Introgression (`f-branch`):** Maximum `f_b` = %.4f (donor `%s`%s).", fmax, donor_s, pol_s))
    }
  }
  if (!is.null(land_df) && nrow(land_df) > 0L && "mmrr_beta_ibd" %in% names(land_df)) {
    b_ibd <- suppressWarnings(as.numeric(land_df$mmrr_beta_ibd[1L]))
    b_ibe <- suppressWarnings(as.numeric(land_df$mmrr_beta_ibe[1L]))
    p_ibe <- suppressWarnings(as.numeric(land_df$mmrr_p_ibe[1L]))
    if (is.finite(b_ibd) || is.finite(b_ibe)) {
      lines <- c(lines, sprintf(
        "- **Landscape Genomics (`MMRR` / `conStruct` / `FEEMS`):** `beta_IBD` = %s, `beta_IBE` = %s (`p_IBE` = %s).",
        ifelse(is.finite(b_ibd), sprintf("%.3f", b_ibd), "NA"),
        ifelse(is.finite(b_ibe), sprintf("%.3f", b_ibe), "NA"),
        ifelse(is.finite(p_ibe), sprintf("%.4f", p_ibe), "NA")
      ))
    }
  }
  if (!is.null(gc_df) && nrow(gc_df) > 0L) {
    c_strs <- vapply(seq_len(nrow(gc_df)), function(k) {
      sprintf("%s (n=%d, Mantel r=%s)", gc_df$comparison[k], as.integer(gc_df$n_samples[k]), ifelse(is.finite(suppressWarnings(as.numeric(gc_df$mantel_r[k]))), sprintf("%.3f", as.numeric(gc_df$mantel_r[k])), "NA"))
    }, character(1))
    lines <- c(lines, paste0("- **Post-Hoc Multi-Modal Congruence:** ", paste(c_strs, collapse = "; "), "."))
  }
  tc_src <- if (!is.null(tc_df) && nrow(tc_df) > 0L) tc_df else df
  cx_col <- intersect(c("consensus_complex_id", "complex_id", "complex_label", "species_complex"), names(tc_src))[1L]
  cx_str <- ""
  if (!is.na(cx_col)) {
    cx_vals <- as.character(tc_src[[cx_col]])
    cx_vals <- cx_vals[!is.na(cx_vals) & nzchar(trimws(cx_vals)) & !grepl("^(NA|Singleton|Unassigned)$", cx_vals, ignore.case = TRUE)]
    if (length(cx_vals) > 0L) cx_str <- sprintf("; Step 07 Leiden complexes represented (%d): %s", length(unique(cx_vals)), paste(utils::head(unique(cx_vals), 5L), collapse = ", "))
  }
  if (!is.null(pbc_df) && nrow(pbc_df) > 0L) {
    n_mo <- if ("in_mo_backbone" %in% names(pbc_df)) sum(.safe_bool_vec_r(pbc_df$in_mo_backbone), na.rm = TRUE) else 0L
    n_ap <- if ("in_astralpro_backbone" %in% names(pbc_df)) sum(.safe_bool_vec_r(pbc_df$in_astralpro_backbone), na.rm = TRUE) else 0L
    cmask <- if ("mo_astralpro_concordant" %in% names(pbc_df)) !is.na(pbc_df$mo_astralpro_concordant) else rep(FALSE, nrow(pbc_df))
    n_both <- sum(cmask); n_conc <- if (n_both > 0L) sum(.safe_bool_vec_r(pbc_df$mo_astralpro_concordant[cmask])) else 0L
    jacc_v <- if ("mo_astralpro_jaccard" %in% names(pbc_df)) suppressWarnings(as.numeric(pbc_df$mo_astralpro_jaccard)) else numeric(0); jacc_v <- jacc_v[is.finite(jacc_v)]
    clades <- if ("public_backbone_clade" %in% names(pbc_df)) unique(as.character(pbc_df$public_backbone_clade[!is.na(pbc_df$public_backbone_clade) & nzchar(pbc_df$public_backbone_clade)])) else character(0)
    ils_v <- if ("public_phytop_ils_index" %in% names(pbc_df)) suppressWarnings(as.numeric(pbc_df$public_phytop_ils_index)) else numeric(0)
    ih_v <- if ("public_phytop_ih_index" %in% names(pbc_df)) suppressWarnings(as.numeric(pbc_df$public_phytop_ih_index)) else numeric(0)
    hot_str <- if (any(is.finite(ih_v)) || any(is.finite(ils_v))) {
      hi <- if (any(is.finite(ih_v))) which.max(ifelse(is.finite(ih_v), ih_v, -Inf)) else which.max(ifelse(is.finite(ils_v), ils_v, -Inf))
      htx <- if ("taxon" %in% names(pbc_df)) as.character(pbc_df$taxon[hi]) else "NA"
      sprintf("PhyTop ILS vs IH hotspot: `%s` (ILS_index=%s, IH_index=%s; median ILS=%s, max IH=%s)", htx,
              if (length(ils_v) >= hi && is.finite(ils_v[hi])) sprintf("%.3f", ils_v[hi]) else "NA",
              if (length(ih_v) >= hi && is.finite(ih_v[hi])) sprintf("%.3f", ih_v[hi]) else "NA",
              if (any(is.finite(ils_v))) sprintf("%.3f", stats::median(ils_v[is.finite(ils_v)])) else "NA",
              if (any(is.finite(ih_v))) sprintf("%.3f", max(ih_v[is.finite(ih_v)])) else "NA")
    } else "PhyTop ILS_index/IH_index not staged"
    lines <- c(lines, sprintf(
      "- **Public 111-Taxon Nuclear Backbone (MO vs. ASTRAL-Pro & PhyTop):** %d/%d taxa mapped to MO and %d to ASTRAL-Pro (%d/%d shared clades concordant, median Jaccard=%s); clades: %s; %s%s.",
      n_mo, nrow(pbc_df), n_ap, n_conc, n_both,
      if (length(jacc_v) > 0L) sprintf("%.3f", stats::median(jacc_v)) else "NA",
      if (length(clades) > 0L) paste(utils::head(clades, 4L), collapse = ", ") else "none", hot_str, cx_str
    ))
  }
  if (!is.null(pop_df) && nrow(pop_df) > 0L) {
    n_pl <- sum(as.character(pop_df$locus) == "plastid", na.rm = TRUE)
    n_it <- sum(as.character(pop_df$locus) == "its", na.rm = TRUE)
    n_wx <- sum(as.character(pop_df$locus) == "waxy", na.rm = TRUE)
    top_tx <- unique(as.character(pop_df$best_match_taxon[!is.na(pop_df$best_match_taxon) & nzchar(pop_df$best_match_taxon)]))
    med_c <- if ("kmer_containment" %in% names(pop_df)) suppressWarnings(stats::median(as.numeric(pop_df$kmer_containment), na.rm = TRUE)) else NA_real_
    lines <- c(lines, sprintf(
      "- **GenBank/BOLD Plastid-ITS-waxy Placement:** %d sample-locus matches (plastid=%d, ITS=%d, waxy=%d; median 11-mer containment=%s; top matched references: %s)%s.",
      nrow(pop_df), n_pl, n_it, n_wx, if (is.finite(med_c)) sprintf("%.3f", med_c) else "NA",
      if (length(top_tx) > 0L) paste(utils::head(top_tx, 5L), collapse = ", ") else "none",
      if (is.null(pbc_df) || nrow(pbc_df) == 0L) cx_str else ""
    ))
  } else if (!is.null(ga_df) && nrow(ga_df) > 0L && any(c("public_plastome_match", "public_its_match", "waxy_lineage_match") %in% names(ga_df))) {
    cnt_nn <- function(col) if (col %in% names(ga_df)) sum(!is.na(ga_df[[col]]) & nzchar(trimws(as.character(ga_df[[col]])))) else 0L
    n_pl <- cnt_nn("public_plastome_match"); n_it <- cnt_nn("public_its_match"); n_wx <- cnt_nn("waxy_lineage_match")
    if ((n_pl + n_it + n_wx) > 0L) {
      lines <- c(lines, sprintf("- **GenBank/BOLD Plastid-ITS-waxy Placement:** plastid=%d, ITS=%d, waxy=%d sample matches recorded in `genomic_assignments.csv`%s.", n_pl, n_it, n_wx, if (is.null(pbc_df) || nrow(pbc_df) == 0L) cx_str else ""))
    }
  }
  prop_sum <- attr(df, "genomic_propagation_summary")
  n_ex <- if (!is.null(prop_sum)) as.integer(prop_sum$n_g1_exsiccatae) else if ("genomic_propagation_tier" %in% names(df)) sum(df$genomic_propagation_tier == "G1_Exsiccatae_Isoduplicate", na.rm = TRUE) else 0L
  n_dh <- if (!is.null(prop_sum)) as.integer(prop_sum$n_g2_deme_hull) else if ("genomic_propagation_tier" %in% names(df)) sum(df$genomic_propagation_tier == "G2_Conditional_Deme_Hull", na.rm = TRUE) else 0L
  n_t0a <- if (!is.null(prop_sum)) as.integer(prop_sum$n_tier0_anchors) else 0L
  n_t0c <- if (!is.null(prop_sum)) as.integer(prop_sum$n_tier0_calibrated) else if ("tier0_calibrated_taxon" %in% names(df) && "tier0_abstention_flag" %in% names(df)) sum(!is.na(df$tier0_calibrated_taxon) & df$tier0_abstention_flag == 0L, na.rm = TRUE) else 0L
  n_t0s <- if (!is.null(prop_sum)) as.integer(prop_sum$n_tier0_abstained) else if ("tier0_abstention_flag" %in% names(df)) sum(df$tier0_abstention_flag == 1L, na.rm = TRUE) else 0L
  n_cg <- if (!is.null(prop_sum)) as.integer(prop_sum$n_cluster_gdi) else if ("cluster_propagated_gdi" %in% names(df)) sum(is.finite(df$cluster_propagated_gdi)) else 0L
  if ((n_ex + n_dh + n_t0a + n_t0c + n_t0s + n_cg) > 0L) {
    k_str <- if (!is.null(prop_sum) && is.finite(prop_sum$lopo_cv_kappa)) sprintf(", LOPO-CV kappa = %.3f", prop_sum$lopo_cv_kappa) else ""
    lines <- c(lines, sprintf("- **Step 07 Genomic Propagation & Tier_0 Calibration:** %d exsiccatae isoduplicates (`G1`), %d conditional deme priors (`G2`), %d `Tier_0` DNA anchors, %d `Tier_0`-calibrated projections (%d abstained as `Sensu_Lato_DNA_Required`), %d cluster-level `gdi` transfers%s.", n_ex, n_dh, n_t0a, n_t0c, n_t0s, n_cg, k_str))
  }
  c(lines, .build_cytology_summary_lines_r(df))
}

(function(.env) {
  req <- c(".safe_bool_vec_r", ".coalesce_chr_cols_r", ".resolve_step07_genomic_tables_dir",
           ".resolve_step07_genomic_reports_dir", "has_active_genomic_row_mask", ".align_packera_orthography")
  if (all(vapply(req, exists, logical(1L), envir = .env, inherits = TRUE, mode = "function"))) return(invisible(NULL))
  ofiles <- as.character(unlist(lapply(sys.frames(), function(fr) fr$ofile)))
  dirs <- dirname(c(rev(ofiles), sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
  dirs <- dirs[!is.na(dirs) & nzchar(dirs)]
  for (p in unique(c(file.path(dirs, "genomic_triage_utils.R"), "scripts/analysis/genomic_triage_utils.R"))) {
    if (nzchar(p) && file.exists(p)) { source(p, local = .env); break }
  }
})(environment())
