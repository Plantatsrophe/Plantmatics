if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/auth_concept_type_masks.R
# Type specimen identification, concept/hybrid quarantine, and tier masks
# ==============================================================================

is_extralimital_or_conflicted_non_type <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(logical(0))
  }
  n <- nrow(df)
  t_info <- .identify_type_specimen_vec(df)
  is_t <- as.logical(t_info$is_type_specimen)
  is_t[is.na(is_t)] <- FALSE
  usda_st <- if ("usda_county_status" %in% names(df)) {
    tolower(trimws(as.character(df$usda_county_status)))
  } else {
    rep("", n)
  }
  usda_st[is.na(usda_st)] <- ""
  is_out <- usda_st %in% c("out_of_range", "natureserve_exotic_or_waif")
  ex_conf <- if ("exsiccatae_taxonomic_conflict" %in% names(df)) {
    raw_ec <- df$exsiccatae_taxonomic_conflict
    if (is.logical(raw_ec)) {
      !is.na(raw_ec) & raw_ec
    } else if (is.numeric(raw_ec)) {
      !is.na(raw_ec) & (raw_ec != 0)
    } else {
      tolower(trimws(as.character(raw_ec))) %in% c("true", "1", "t", "yes")
    }
  } else {
    rep(FALSE, n)
  }
  (!is_t) & (is_out | ex_conf)
}

.identify_type_specimen_vec <- function(df) {
  if (is.null(df)) {
    return(list(type_rank = integer(0), is_type_specimen = logical(0), type_status_clean = character(0)))
  }
  if (!is.data.frame(df)) {
    raw_vec <- as.character(df)
    n <- length(raw_vec)
    explicit_bool <- rep(FALSE, n)
  } else {
    n <- nrow(df)
    if (n == 0L) {
      return(list(type_rank = integer(0), is_type_specimen = logical(0), type_status_clean = character(0)))
    }
    raw_vec <- if ("type_status" %in% names(df)) {
      as.character(df$type_status)
    } else if ("typeStatus" %in% names(df)) {
      as.character(df$typeStatus)
    } else {
      rep("", n)
    }
    explicit_bool <- if ("is_type_specimen" %in% names(df)) {
      b <- as.logical(df$is_type_specimen)
      b[is.na(b)] <- FALSE
      b
    } else {
      rep(FALSE, n)
    }
  }
  raw_vec[is.na(raw_vec)] <- ""
  raw_trim <- trimws(raw_vec)
  neg_mask <- grepl("not a type|non-type|non type|false type|no type|unverified|type locality|type species|type collection", raw_trim, ignore.case = TRUE)
  prim_mask <- !neg_mask & grepl("(^|[^A-Za-z0-9_])(holotype|lectotype|neotype|epitype|isotype|isolectotype|isoneotype)s?([^A-Za-z0-9_]|$)", raw_trim, ignore.case = TRUE)
  syntype_mask <- !neg_mask & !prim_mask & grepl("(^|[^A-Za-z0-9_])(syntype|isosyntype)s?([^A-Za-z0-9_]|$)", raw_trim, ignore.case = TRUE)
  sec_mask <- !neg_mask & !prim_mask & grepl("(^|[^A-Za-z0-9_])(syntype|isosyntype|paratype|isoparatype|type|cotype)s?([^A-Za-z0-9_]|$)", raw_trim, ignore.case = TRUE)

  type_rank <- integer(n)
  type_rank[sec_mask] <- 1L
  type_rank[prim_mask] <- 2L
  fb_mask <- (type_rank == 0L) & explicit_bool & !neg_mask
  if (any(fb_mask)) {
    type_rank[fb_mask] <- 1L
  }
  is_nb_type <- prim_mask | syntype_mask
  if (is.data.frame(df) && ("species_standardized" %in% names(df)) && any(is_nb_type)) {
    sp_std <- tolower(trimws(as.character(df$species_standardized)))
    sp_ep <- sub("^[a-z]+[[:space:]]+([a-z-]+).*$", "\\1", sp_std)
    has_of_binomial <- grepl("(^|[^A-Za-z0-9_])of[[:space:]]+(packera|senecio|cineraria)[[:space:]]+[a-z]{3,}", raw_trim, ignore.case = TRUE)
    if (any(is_nb_type & has_of_binomial)) {
      of_binom <- sub("(^|.*[^A-Za-z0-9_])of[[:space:]]+((packera|senecio|cineraria)[[:space:]]+[a-z-]+).*$", "\\2", tolower(raw_trim))
      of_ep <- sub("^[a-z]+[[:space:]]+([a-z-]+).*$", "\\1", of_binom)
      stem6 <- function(s) substr(gsub("(us|a|um)$", "", s), 1L, 6L)
      dyn_equiv <- rep(FALSE, n)
      if (exists("standardize_taxon_name", mode = "function")) {
        chk_idx <- which(is_nb_type & has_of_binomial & nzchar(of_binom))
        if (length(chk_idx) > 0L) {
          resolved_of <- tolower(trimws(as.character(standardize_taxon_name(of_binom[chk_idx]))))
          dyn_equiv[chk_idx] <- !is.na(resolved_of) & nzchar(resolved_of) & (resolved_of == sp_std[chk_idx])
        }
      }
      het_mismatch <- is_nb_type & has_of_binomial & nzchar(sp_ep) & nzchar(of_ep) &
        (sp_ep != of_ep) & (stem6(sp_ep) != stem6(of_ep)) & !dyn_equiv
      is_nb_type[het_mismatch] <- FALSE
    }
  }
  clean_lbl <- ifelse(type_rank >= 1L & nzchar(raw_trim), raw_trim, ifelse(type_rank >= 1L, "TYPE", "None"))
  list(
    type_rank = as.integer(type_rank),
    is_type_specimen = as.logical(type_rank >= 1L),
    is_name_bearing_type = as.logical(is_nb_type),
    type_status_clean = as.character(clean_lbl)
  )
}

# Identify non-type open-nomenclature (cf./aff./?/s.l.) or two-epithet hybrid formulas.
.identify_open_or_hybrid_concept_mask <- function(df) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L) return(logical(0))
  res <- rep(FALSE, n)
  if ("is_open_or_hybrid_concept" %in% names(df)) {
    v_oh <- df$is_open_or_hybrid_concept
    if (is.logical(v_oh)) {
      res <- res | (!is.na(v_oh) & v_oh)
    } else if (is.numeric(v_oh) || is.integer(v_oh)) {
      res <- res | (!is.na(v_oh) & v_oh != 0)
    } else {
      sv <- tolower(trimws(as.character(v_oh)))
      res <- res | (!is.na(sv) & (sv %in% c("1", "true", "t", "yes")))
    }
  }
  if ("concept_synonymy_status" %in% names(df)) {
    cs <- trimws(as.character(df$concept_synonymy_status))
    res <- res | (!is.na(cs) & (cs %in% c("Open_Nomenclature", "Hybrid_Formula")))
  }
  # Match two-epithet hybrid formulas without flagging single-epithet nothospecies (e.g. Packera ×memmingeri).
  hyb_pat <- "[A-Z][a-z]+[[:space:]]+[a-z][a-z-]+([[:space:]]*×|[[:space:]]+[xX][[:space:]]+)([A-Z][a-z]*\\.?[[:space:]]+)?[a-z][a-z-]+"
  open_pat <- "(^|[[:space:][:punct:]])(cf|aff|vel[[:space:]]+aff)\\.?($|[[:space:][:punct:]])|\\?|sensu[[:space:]]+lato|(^|[[:space:][:punct:]])s\\.l\\.($|[[:space:][:punct:]])|putative[[:space:]]+hybrid|hybrid[[:space:]]+(with|between)|intermediate[[:space:]]+(between|toward)|intergrade[[:space:]]+with"
  for (col_nm in intersect(c("verbatim_taxon_concept", "species_raw", "scientificName", "identificationQualifier", "verbatimIdentification", "identificationRemarks"), names(df))) {
    txt <- trimws(as.character(df[[col_nm]]))
    ok <- !is.na(txt) & nzchar(txt)
    if (!any(ok)) next
    hit_hyb <- grepl(hyb_pat, txt[ok])
    hit_open <- grepl(open_pat, txt[ok], ignore.case = TRUE)
    res[ok] <- res[ok] | hit_hyb | hit_open
  }
  res
}

# Identify non-type vouchers mapped from a heterotypic synonym rather than sensu stricto concept.
.is_heterotypic_synonym_mask <- function(df) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L || !("concept_synonymy_status" %in% names(df))) return(rep(FALSE, n))
  cs <- trimws(as.character(df$concept_synonymy_status))
  !is.na(cs) & (cs == "Heterotypic_Synonym")
}

.standardize_determiner_tier_vec <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  raw <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_3_Bronze", n)
  raw[is.na(raw) | !nzchar(trimws(raw))] <- "Tier_3_Bronze"
  out <- rep("Tier_3_Bronze", n)
  out[grepl("tier_?2|silver", raw, ignore.case = TRUE)] <- "Tier_2_Silver"
  out[grepl("tier_?[01]|gold|type", raw, ignore.case = TRUE)] <- "Tier_1_Gold"
  t_info <- .identify_type_specimen_vec(df)
  # Quarantine non-type open-nomenclature and hybrid-formula sheets to Tier_3_Bronze.
  is_open_hyb <- .identify_open_or_hybrid_concept_mask(df)
  if (any(is_open_hyb & !t_info$is_type_specimen)) {
    out[is_open_hyb & !t_info$is_type_specimen] <- "Tier_3_Bronze"
  }
  if (any(t_info$is_type_specimen)) {
    out[t_info$is_type_specimen] <- "Tier_1_Gold"
  }
  out
}

.resolve_authority_taxon_vec <- function(df, taxon_col = "species_standardized") {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  cand_cols <- unique(c(taxon_col, "species_standardized", "grouping_col", "scientificName", "species_raw"))
  chosen_col <- NULL
  for (cc in cand_cols) {
    if (!is.null(cc) && nzchar(cc) && (cc %in% names(df))) {
      vals <- trimws(as.character(df[[cc]]))
      if (any(!is.na(vals) & nzchar(vals))) {
        chosen_col <- cc
        break
      }
    }
  }
  if (is.null(chosen_col)) return(rep("Unassigned_Packera", n))
  tx <- trimws(as.character(df[[chosen_col]]))
  if (chosen_col %in% c("species_standardized", "scientificName", "species_raw") && !any(grepl("^Morphotype_K|^Ecotype_", tx[!is.na(tx)]))) {
    if (exists("standardize_taxon_name", mode = "function")) {
      std_tx <- as.character(standardize_taxon_name(tx))
      tx <- ifelse(!is.na(std_tx) & nzchar(std_tx), std_tx, tx)
    } else if (exists("standardize_packera_taxon", mode = "function")) {
      tx <- as.character(standardize_packera_taxon(tx))
    }
  }
  tx[is.na(tx) | !nzchar(tx) | tolower(tx) == "unknown"] <- "Unassigned_Packera"
  tx
}

.identify_tier1_pure_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  if ("morph_tier_source" %in% names(df)) {
    mts <- trimws(as.character(df$morph_tier_source))
    is_t1 <- is.na(mts) | !nzchar(mts) | !(mts %in% c("Cauline_Only_Fallback", "Transitional_Passive", "No_Leaves"))
  } else {
    is_t1 <- rep(TRUE, n)
  }
  is_t1[is.na(is_t1)] <- FALSE
  if ("is_contour_outlier" %in% names(df)) {
    out_flag <- as.logical(df$is_contour_outlier)
    out_flag[is.na(out_flag)] <- FALSE
    if (sum(is_t1 & !out_flag) >= 2L) is_t1 <- is_t1 & !out_flag
  }
  is_t1
}

.discover_numeric_sorted_cols <- function(nms, pattern) {
  hits <- grep(pattern, nms, value = TRUE)
  if (length(hits) <= 1L) return(hits)
  nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", hits)))
  hits[order(nums, hits)]
}
