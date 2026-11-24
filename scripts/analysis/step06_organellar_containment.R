# ==============================================================================
# Module: scripts/analysis/step06_organellar_containment.R
# Project: Plantmatics | Affiliation: UNC-Chapel Hill Herbarium (NCU)
# Purpose: Quality-gated skimmed plastome / nrITS / nuclear waxy 11-mer matching
#   and containment-gate screening for Step 06 public organellar placement.
# ==============================================================================

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists("screen_public_fasta_outliers", mode = "function") || !exists(".empty_rogue_audit_df", mode = "function")) {
  .oc_srcs <- c(
    "scripts/analysis/step06_rogue_sequence_qc.R",
    "step06_rogue_sequence_qc.R",
    file.path(dirname(as.character(c(
      sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
      unlist(lapply(sys.frames(), function(fr) fr$ofile))
    ))), "step06_rogue_sequence_qc.R")
  )
  for (.oc_p in .oc_srcs) if (!is.na(.oc_p) && nzchar(.oc_p) && file.exists(.oc_p)) { source(.oc_p); break }
}

.find_organellar_file <- function(dirs, fname) {
  if (exists(".find_geno_file", mode = "function")) return(.find_geno_file(dirs, fname))
  for (d in dirs) {
    if (is.null(d) || !nzchar(d) || !dir.exists(d)) next
    p <- file.path(d, fname)
    if (file.exists(p)) return(p)
    cand <- sort(list.files(d, full.names = TRUE, recursive = TRUE))
    hits <- cand[tolower(basename(cand)) == tolower(fname)]
    if (length(hits) > 0L) return(hits[1L])
  }
  ""
}

.norm_organellar_taxon <- function(raw) {
  if (exists(".normalize_public_packera_taxon", mode = "function")) return(.normalize_public_packera_taxon(raw))
  if (exists(".norm_taxon_qc", mode = "function")) return(.norm_taxon_qc(raw))
  trimws(gsub("_+", " ", ifelse(is.na(raw), "", as.character(raw))))
}

.read_fasta_named <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(data.frame(header = character(0), seq = character(0), stringsAsFactors = FALSE))
  lns <- readLines(path, warn = FALSE)
  hdr_idx <- grep("^>", lns)
  if (length(hdr_idx) == 0L) return(data.frame(header = character(0), seq = character(0), stringsAsFactors = FALSE))
  ends <- c(hdr_idx[-1L] - 1L, length(lns))
  hdrs <- sub("^>", "", trimws(lns[hdr_idx]))
  seqs <- vapply(seq_along(hdr_idx), function(i) {
    if (hdr_idx[i] + 1L > ends[i]) return("")
    gsub("[^A-Z]", "", toupper(paste(lns[(hdr_idx[i] + 1L):ends[i]], collapse = "")))
  }, character(1L))
  keep <- nzchar(hdrs) & nchar(seqs) >= 11L
  data.frame(header = hdrs[keep], seq = seqs[keep], stringsAsFactors = FALSE)
}

.kmers_11 <- function(s) {
  n <- nchar(s)
  if (is.na(n) || n < 11L) return(character(0))
  st <- seq_len(n - 10L)
  k <- substring(s, st, st + 10L)
  unique(k[!grepl("[^ACGT]", k)])
}

.revcomp_seq <- function(s) {
  ch <- rev(strsplit(s, "", fixed = TRUE)[[1L]])
  chartr("ACGTU", "TGCAA", paste(ch, collapse = ""))
}

.empty_organellar_placement_df <- function() {
  data.frame(
    genetic_sample_id = character(0), locus = character(0), best_match_taxon = character(0),
    best_ref_header = character(0), ref_source_file = character(0), n_shared_11mers = integer(0),
    kmer_containment = numeric(0), kmer_jaccard = numeric(0), containment_gate_pass = integer(0),
    stringsAsFactors = FALSE
  )
}

.match_skimmed_organellar_11mers <- function(
  geno_dir,
  sample_df,
  min_containment = 0.80,
  min_shared_11mers = 80L,
  out_dir = ""
) {
  sample_ids <- if (is.data.frame(sample_df)) {
    if ("genetic_sample_id" %in% names(sample_df)) as.character(sample_df$genetic_sample_id) else character(0)
  } else if (is.character(sample_df)) {
    as.character(sample_df)
  } else character(0)
  n_s <- length(sample_ids)
  search_dirs <- unique(Filter(
    function(d) !is.null(d) && nzchar(d) && dir.exists(d),
    c(geno_dir, out_dir, if (nzchar(out_dir)) dirname(out_dir) else "", if (nzchar(out_dir)) file.path(dirname(out_dir), "genomics") else "")
  ))
  locus_specs <- list(
    list(locus = "plastid", qfile = "skimmed_plastid_consensus.fasta",
         refs = c("public_packera_111mapping_plastomes.fasta", "public_packera_18complete_plastomes.fasta", "public_packera_plastomes_and_cpdna.fasta")),
    list(locus = "its", qfile = "skimmed_its_consensus.fasta",
         refs = c("public_packera_its_ets.fasta", "public_packera_bold_barcodes.fasta", "bain_jansen_sequences.fasta")),
    list(locus = "waxy", qfile = "skimmed_waxy_consensus.fasta",
         refs = c("public_packera_waxy_nuclear.fasta"))
  )
  rows <- list(); audit_list <- list()
  min_c <- as.numeric(min_containment); min_sh <- as.integer(min_shared_11mers)
  for (sp in locus_specs) {
    ref_list <- list()
    for (rf in sp$refs) {
      rp <- .find_organellar_file(search_dirs, rf)
      rdf <- .read_fasta_named(rp)
      if (nrow(rdf) > 0L) {
        rdf$taxon <- .norm_organellar_taxon(rdf$header)
        rdf <- rdf[!is.na(rdf$taxon) & nzchar(rdf$taxon), , drop = FALSE]
        if (nrow(rdf) > 0L) { rdf$ref_source <- rf; ref_list[[length(ref_list) + 1L]] <- rdf }
      }
    }
    if (length(ref_list) == 0L) next
    refs <- do.call(rbind, ref_list)
    r_km <- lapply(refs$seq, .kmers_11)
    qc_fa <- screen_public_fasta_outliers(refs, sp$locus, kmers_list = r_km)
    if (!is.null(qc_fa$audit_df) && nrow(qc_fa$audit_df) > 0L) audit_list[[length(audit_list) + 1L]] <- qc_fa$audit_df
    refs <- qc_fa$clean_refs; r_km <- qc_fa$clean_kmers
    qp <- .find_organellar_file(search_dirs, sp$qfile)
    qdf <- .read_fasta_named(qp)
    if (nrow(qdf) == 0L || nrow(refs) == 0L) next
    r_len <- vapply(r_km, length, integer(1L)); valid_r <- which(r_len > 0L)
    if (length(valid_r) == 0L) next
    r_sk <- lapply(seq_along(r_km), function(ri) if (r_len[ri] > 6000L) r_km[[ri]][endsWith(r_km[[ri]], "AA")] else r_km[[ri]])
    r_sk_len <- vapply(r_sk, length, integer(1L))
    for (qi in seq_len(nrow(qdf))) {
      sid <- sub("\\|.*$", "", sub("[[:space:]].*$", "", qdf$header[qi]))
      if (n_s > 0L && !(sid %in% sample_ids)) next
      qf <- .kmers_11(qdf$seq[qi]); qr <- .kmers_11(.revcomp_seq(qdf$seq[qi]))
      q_len <- length(qf); if (q_len < 3L) next
      use_q_sk <- q_len > 6000L
      qf_u <- if (use_q_sk) qf[endsWith(qf, "AA")] else qf
      qr_u <- if (use_q_sk) qr[endsWith(qr, "AA")] else qr
      q_eff <- length(qf_u); if (q_eff < 3L) { qf_u <- qf; qr_u <- qr; q_eff <- q_len; use_q_sk <- FALSE }
      best <- list(score = -1, sh = 0L, cont = 0.0, jac = 0.0, ri = NA_integer_)
      for (ri in valid_r) {
        rk <- if (use_q_sk) r_sk[[ri]] else r_km[[ri]]; rl <- if (use_q_sk) r_sk_len[ri] else r_len[ri]
        if (rl < 3L) next
        sh <- max(sum(qf_u %in% rk), sum(qr_u %in% rk)); if (sh < 3L) next
        cont <- sh / q_eff; if (cont < 0.05) next
        jac <- sh / (q_eff + rl - sh); sc <- cont + 1e-4 * jac
        if (sc > best$score) best <- list(score = sc, sh = as.integer(sh), cont = cont, jac = jac, ri = ri)
      }
      if (!is.na(best$ri)) {
        gate_pass <- as.integer(is.finite(best$cont) && best$cont >= min_c && as.integer(best$sh) >= min_sh)
        rows[[length(rows) + 1L]] <- data.frame(
          genetic_sample_id = sid, locus = sp$locus, best_match_taxon = refs$taxon[best$ri],
          best_ref_header = refs$header[best$ri], ref_source_file = refs$ref_source[best$ri],
          n_shared_11mers = as.integer(best$sh), kmer_containment = round(best$cont, 4),
          kmer_jaccard = round(best$jac, 4), containment_gate_pass = gate_pass, stringsAsFactors = FALSE
        )
      }
    }
  }
  res_df <- if (length(rows) > 0L) do.call(rbind, rows) else .empty_organellar_placement_df()
  fasta_audit <- if (length(audit_list) > 0L) do.call(rbind, audit_list) else .empty_rogue_audit_df()
  attr(res_df, "fasta_rogue_audit") <- fasta_audit
  extract_gated <- function(loc) {
    out <- rep(NA_character_, n_s); c_v <- rep(NA_real_, n_s)
    sub <- res_df[res_df$locus == loc & res_df$containment_gate_pass == 1L, , drop = FALSE]
    if (nrow(sub) > 0L && n_s > 0L) {
      m <- match(sample_ids, sub$genetic_sample_id); ok <- !is.na(m)
      out[ok] <- sub$best_match_taxon[m[ok]]; c_v[ok] <- sub$kmer_containment[m[ok]]
    }
    list(match = out, cont = c_v)
  }
  cp_r <- extract_gated("plastid"); its_r <- extract_gated("its"); wx_r <- extract_gated("waxy")
  comb_ok <- !is.na(its_r$match) & !is.na(wx_r$match) & (tolower(its_r$match) == tolower(wx_r$match)) &
    is.finite(its_r$cont) & is.finite(wx_r$cont) & (0.5 * (its_r$cont + wx_r$cont) >= 0.85)
  list(
    placement_df = res_df, public_plastome_match = cp_r$match, public_its_match = its_r$match,
    waxy_lineage_match = wx_r$match, combined_its_waxy_concordant = comb_ok, fasta_rogue_audit = fasta_audit
  )
}

.place_organellar_loci <- function(raw_geno_dir, out_dir, sample_ids, min_containment = 0.80, min_shared_11mers = 80L) {
  .match_skimmed_organellar_11mers(raw_geno_dir, sample_ids, min_containment = min_containment, min_shared_11mers = min_shared_11mers, out_dir = out_dir)$placement_df
}
