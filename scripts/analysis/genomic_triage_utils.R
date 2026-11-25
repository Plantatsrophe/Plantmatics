if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/analysis/genomic_triage_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Helper utilities bridging the independent Step 06 Genomic & Phylogenomic
#   Evidence Stream into Step 07 Multi-Evidence Synthesis without altering upstream
#   Step 03-05 modalities or legacy 3-stream scores on unsequenced vouchers.
# ==============================================================================

.safe_bool_vec_r <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  s <- toupper(trimws(as.character(x)))
  s %in% c("1", "TRUE", "T", "YES", "Y")
}

.coalesce_chr_cols_r <- function(tbl, cols, default = NA_character_) {
  n <- if (is.data.frame(tbl)) nrow(tbl) else 0L
  if (n == 0L) return(character(0))
  out <- rep(NA_character_, n)
  for (cn in cols) {
    if (!(cn %in% names(tbl))) next
    v <- trimws(as.character(tbl[[cn]]))
    fill <- (is.na(out) | !nzchar(out)) & !is.na(v) & nzchar(v)
    if (any(fill)) out[fill] <- v[fill]
  }
  if (!is.na(default)) out[is.na(out) | !nzchar(out)] <- default
  out
}

.coalesce_num_cols_r <- function(tbl, cols, default = NA_real_) {
  n <- if (is.data.frame(tbl)) nrow(tbl) else 0L
  if (n == 0L) return(numeric(0))
  out <- rep(NA_real_, n)
  for (cn in cols) {
    if (!(cn %in% names(tbl))) next
    v <- suppressWarnings(as.numeric(tbl[[cn]]))
    fill <- !is.finite(out) & is.finite(v)
    if (any(fill)) out[fill] <- v[fill]
  }
  out[!is.finite(out)] <- default
  out
}

.resolve_step07_genomic_tables_dir <- function(opts) {
  q_dir <- if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) dirname(opts$output_queue) else "outputs/tables"
  g_dirs <- c(
    if (!is.null(opts$genomic_dir) && nzchar(opts$genomic_dir)) opts$genomic_dir else character(0),
    if (!is.null(opts$genomics_dir) && nzchar(opts$genomics_dir)) opts$genomics_dir else character(0)
  )
  cands <- unique(c(q_dir, g_dirs, file.path(dirname(q_dir), "reports"), file.path(dirname(dirname(q_dir)), "outputs", "reports")))
  for (cd in cands) {
    if (nzchar(cd) && (file.exists(file.path(cd, "genomic_assignments.csv")) || file.exists(file.path(cd, "genomic_model_metrics.csv")))) {
      return(cd)
    }
  }
  q_dir
}

.resolve_step07_genomic_reports_dir <- function(opts) {
  q_dir <- if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) dirname(opts$output_queue) else "outputs/tables"
  rep_dir <- if (!is.null(opts$output_report) && nzchar(opts$output_report)) dirname(opts$output_report) else character(0)
  g_dirs <- c(
    if (!is.null(opts$genomic_dir) && nzchar(opts$genomic_dir)) opts$genomic_dir else character(0),
    if (!is.null(opts$genomics_dir) && nzchar(opts$genomics_dir)) opts$genomics_dir else character(0)
  )
  cands <- unique(c(g_dirs, rep_dir, file.path(dirname(q_dir), "reports"), file.path(dirname(dirname(q_dir)), "outputs", "reports"), q_dir))
  for (cd in cands) {
    if (nzchar(cd) && (
      file.exists(file.path(cd, "coalescent_gdi_delimitation.csv")) ||
      file.exists(file.path(cd, "landscape_mmrr_feems_summary.csv")) ||
      file.exists(file.path(cd, "genomic_assignments.csv"))
    )) {
      return(cd)
    }
  }
  .resolve_step07_genomic_tables_dir(opts)
}

has_active_genomic_row_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  nms <- names(df)
  if (!("geno_predicted_taxon" %in% nms)) return(rep(FALSE, n))
  pred <- as.character(df$geno_predicted_taxon)
  pred_ok <- !is.na(pred) & nzchar(trimws(pred)) & (pred != "Unknown")
  tier_col <- if (any(c("geno_linkage_tier", "linkage_tier") %in% nms)) {
    .coalesce_chr_cols_r(df, c("geno_linkage_tier", "linkage_tier"), default = "")
  } else rep("G1_Direct_Voucher", n)
  tier_ok <- tier_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G2_Conditional_Deme_Hull")
  qc_ok <- if ("geno_quality_pass" %in% nms) .safe_bool_vec_r(df$geno_quality_pass) else rep(TRUE, n)
  prob_ok <- if ("geno_posterior_prob" %in% nms) is.finite(suppressWarnings(as.numeric(df$geno_posterior_prob))) else rep(TRUE, n)
  pred_ok & tier_ok & qc_ok & prob_ok
}

.align_packera_orthography <- function(vec, ref_taxa) {
  k <- function(x) {
    s <- tolower(trimws(gsub("[[:space:]]+", " ", gsub("_+", " ", as.character(x)))))
    s <- ifelse(s == "packera millefolia", "packera millefolium", s)
    ifelse(s %in% c("packera memmingeri", "packera x memmingeri", "packera ×memmingeri", "packera × memmingeri", "packera xmemmingeri"), "packera xmemmingeri", s)
  }
  r <- unique(trimws(as.character(ref_taxa[!is.na(ref_taxa) & nzchar(trimws(as.character(ref_taxa)))])))
  if (length(r) == 0L || length(vec) == 0L) return(as.character(vec))
  rmap <- stats::setNames(r, k(r)); out <- as.character(vec); vk <- k(out)
  hit <- !is.na(out) & nzchar(trimws(out)) & vk %in% names(rmap)
  out[hit] <- unname(rmap[vk[hit]]); out
}

load_step07_genomic_metrics <- function(opts, df) {
  t_dir <- .resolve_step07_genomic_tables_dir(opts)
  r_dir <- .resolve_step07_genomic_reports_dir(opts)
  ga_path <- file.path(t_dir, "genomic_assignments.csv")
  if (!file.exists(ga_path)) ga_path <- file.path(r_dir, "genomic_assignments.csv")
  gm_path <- file.path(t_dir, "genomic_model_metrics.csv")
  if (!file.exists(gm_path)) gm_path <- file.path(r_dir, "genomic_model_metrics.csv")

  geno_keep_schema <- c(
    "catalogNumber", "genetic_sample_id", "geno_linkage_tier", "geno_quality_pass",
    "is_putative_novel_or_hybrid_query", "geno_predicted_taxon", "geno_posterior_prob",
    "geno_margin", "genetic_cluster", "q_max", "q_admixture_json", "is_geno_admixed",
    "geno_parent_1", "geno_parent_2", "nearest_tree_sister_taxon",
    "patristic_dist_to_dubia", "replicate_pair_concordance",
    "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl",
    "plastid_clade", "its_ribotype", "cytonuclear_discordant",
    "hybrid_index_h", "interspecific_ancestry_q12", "hybrid_generation_class",
    "coancestry_cluster", "f_branch_max", "f_branch_donor_taxon",
    "snaq_gamma_minor", "coalescent_gdi", "bpp_introgression_bf",
    "mmrr_beta_ibd", "mmrr_beta_ibe", "mmrr_p_ibe",
    "public_plastome_match", "public_its_match", "waxy_lineage_match",
    "public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index",
    "genomic_propagation_tier", "deme_mahal_d2", "geno_deme_weight_scale",
    "tier0_calibrated_taxon", "tier0_calibrated_prob", "tier0_abstention_flag", "cluster_propagated_gdi"
  )

  ga_df <- if (file.exists(ga_path)) tryCatch(utils::read.csv(ga_path, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  has_geno_cols <- ("geno_predicted_taxon" %in% names(df)) && any(!is.na(df$geno_predicted_taxon) & nzchar(as.character(df$geno_predicted_taxon)))
  missing_phase2_cols <- setdiff(geno_keep_schema, names(df))
  need_csv_hydrate <- !is.null(ga_df) && nrow(ga_df) > 0L && ("catalogNumber" %in% names(df)) && (!has_geno_cols || length(missing_phase2_cols) > 0L)

  if (need_csv_hydrate && any(c("catalogNumber", "linked_catalog_number") %in% names(ga_df))) {
    ga_cat <- .coalesce_chr_cols_r(ga_df, c("catalogNumber", "linked_catalog_number"), default = NA_character_)
    ga_df$catalogNumber <- ga_cat
    lt_vec <- .coalesce_chr_cols_r(ga_df, c("geno_linkage_tier", "linkage_tier"), default = NA_character_)
    qc_vec <- if ("geno_quality_pass" %in% names(ga_df)) .safe_bool_vec_r(ga_df$geno_quality_pass) else rep(FALSE, nrow(ga_df))
    keep_tiers <- c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy")
    pass_mask <- !is.na(ga_cat) & nzchar(ga_cat) & qc_vec & (lt_vec %in% keep_tiers)
    if (any(pass_mask)) {
      pass_ga <- ga_df[pass_mask, , drop = FALSE]
      lt_pass <- lt_vec[pass_mask]
      tier_ord <- ifelse(lt_pass == "G1_Direct_Voucher", 1L, ifelse(lt_pass == "G1_Exsiccatae_Isoduplicate", 2L, 3L))
      prob_ord <- .coalesce_num_cols_r(pass_ga, c("geno_posterior_prob", "posterior_prob"), default = -Inf)
      hq_ord <- .coalesce_num_cols_r(pass_ga, c("hq_reads", "geno_hq_reads"), default = 0.0)
      sid_ord <- if ("genetic_sample_id" %in% names(pass_ga)) ifelse(is.na(pass_ga$genetic_sample_id), "", as.character(pass_ga$genetic_sample_id)) else rep("", nrow(pass_ga))
      ord_idx <- order(tier_ord, -prob_ord, -hq_ord, sid_ord)
      pass_ga <- pass_ga[ord_idx, , drop = FALSE]
      lt_pass <- lt_pass[ord_idx]
      uniq_mask <- !duplicated(pass_ga$catalogNumber)
      pass_ga <- pass_ga[uniq_mask, , drop = FALSE]
      pass_ga$geno_linkage_tier <- lt_pass[uniq_mask]
      if (!("geno_posterior_prob" %in% names(pass_ga)) && ("posterior_prob" %in% names(pass_ga))) {
        pass_ga$geno_posterior_prob <- suppressWarnings(as.numeric(pass_ga$posterior_prob))
      }
      keep_cols <- intersect(geno_keep_schema, names(pass_ga))
      m_idx <- match(trimws(as.character(df$catalogNumber)), pass_ga$catalogNumber)
      hit <- !is.na(m_idx)
      if (any(hit)) {
        for (cn in setdiff(keep_cols, "catalogNumber")) {
          if (!(cn %in% names(df))) {
            df[[cn]] <- pass_ga[[cn]][rep(NA_integer_, nrow(df))]
            df[[cn]][hit] <- pass_ga[[cn]][m_idx[hit]]
          } else if (!has_geno_cols) {
            df[[cn]][hit] <- pass_ga[[cn]][m_idx[hit]]
          } else {
            fill_idx <- if (is.character(df[[cn]])) hit & (is.na(df[[cn]]) | !nzchar(trimws(df[[cn]]))) else hit & is.na(df[[cn]])
            if (any(fill_idx)) df[[cn]][fill_idx] <- pass_ga[[cn]][m_idx[fill_idx]]
          }
        }
      }
    }
  }

  n_df <- nrow(df)
  for (cn in c("genomic_propagation_tier", "tier0_calibrated_taxon")) {
    if (!(cn %in% names(df))) df[[cn]] <- rep(NA_character_, n_df)
  }
  for (cn in c("deme_mahal_d2", "geno_deme_weight_scale", "tier0_calibrated_prob", "cluster_propagated_gdi")) {
    if (!(cn %in% names(df))) df[[cn]] <- rep(NA_real_, n_df)
  }
  if (!("tier0_abstention_flag" %in% names(df))) df[["tier0_abstention_flag"]] <- rep(NA_integer_, n_df)

  read_opt_csv <- function(fname) {
    p <- file.path(r_dir, fname)
    if (!file.exists(p)) p <- file.path(t_dir, fname)
    if (file.exists(p)) tryCatch(utils::read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  }

  ref_sp <- if ("species_standardized" %in% names(df)) df$species_standardized else character(0)
  if ("geno_predicted_taxon" %in% names(ga_df)) ga_df$geno_predicted_taxon <- .align_packera_orthography(ga_df$geno_predicted_taxon, ref_sp)
  if ("geno_predicted_taxon" %in% names(df)) df$geno_predicted_taxon <- .align_packera_orthography(df$geno_predicted_taxon, ref_sp)
  pbc_df <- read_opt_csv("public_backbone_concordance.csv")
  if (!is.null(pbc_df) && nrow(pbc_df) > 0L && "species_standardized" %in% names(df) && "taxon" %in% names(pbc_df)) {
    mp <- match(as.character(df$species_standardized), as.character(pbc_df$taxon))
    for (bc_col in intersect(c("public_backbone_clade", "public_phytop_ils_index", "public_phytop_ih_index"), names(pbc_df))) {
      if (!(bc_col %in% names(df))) df[[bc_col]] <- pbc_df[[bc_col]][rep(NA_integer_, nrow(df))]
      need_b <- !is.na(mp) & (is.na(df[[bc_col]]) | !nzchar(trimws(as.character(df[[bc_col]]))))
      if (any(need_b)) df[[bc_col]][need_b] <- pbc_df[[bc_col]][mp[need_b]]
    }
  }
  attr(df, "genomic_assignments_df") <- ga_df
  attr(df, "genomic_model_metrics") <- if (file.exists(gm_path)) tryCatch(utils::read.csv(gm_path, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  attr(df, "genomic_gdi_delimitation") <- read_opt_csv("coalescent_gdi_delimitation.csv")
  attr(df, "genomic_fbranch_summary") <- read_opt_csv("introgressive_fbranch_summary.csv")
  attr(df, "genomic_landscape_summary") <- read_opt_csv("landscape_mmrr_feems_summary.csv")
  df
}

compute_genomic_kappa_weight <- function(df, metrics_df = NULL) {
  active_mask <- has_active_genomic_row_mask(df)
  if (!any(active_mask)) return(NA_real_)
  if (is.null(metrics_df)) metrics_df <- attr(df, "genomic_model_metrics")
  if (!is.null(metrics_df) && nrow(metrics_df) > 0L && "loro_cv_kappa" %in% names(metrics_df)) {
    k_val <- suppressWarnings(as.numeric(metrics_df$loro_cv_kappa[1L]))
    if (is.finite(k_val)) return(max(0.05, min(1.0, k_val)))
  }
  sub_g <- df[active_mask, , drop = FALSE]
  yt <- as.character(sub_g$species_standardized)
  yp <- as.character(sub_g$geno_predicted_taxon)
  ok <- !is.na(yt) & nzchar(yt) & yt != "Unknown" & !is.na(yp) & nzchar(yp)
  if (sum(ok) >= 3L && length(unique(yt[ok])) >= 2L && exists("compute_cohens_kappa", mode = "function")) {
    return(max(0.05, min(1.0, compute_cohens_kappa(yt[ok], yp[ok], default_kappa = 0.85))))
  }
  0.85
}

if (!exists("build_genomic_markdown_summary", envir = environment(), inherits = FALSE, mode = "function") ||
    !exists("export_genomic_field_triage_csv", envir = environment(), inherits = FALSE, mode = "function")) {
  (function(.env) {
    ofiles <- as.character(unlist(lapply(sys.frames(), function(fr) fr$ofile)))
    dirs <- dirname(c(rev(ofiles), sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
    dirs <- dirs[!is.na(dirs) & nzchar(dirs)]
    for (p in unique(c(file.path(dirs, "step07_genomic_report_utils.R"), "scripts/analysis/step07_genomic_report_utils.R"))) {
      if (nzchar(p) && file.exists(p)) { source(p, local = .env); break }
    }
  })(environment())
}
