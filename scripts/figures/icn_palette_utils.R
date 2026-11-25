if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/figures/icn_palette_utils.R
# 24-color colorblind-safe combinatorial scales, ICN botanical typography,
# nomenclatural type overlays, and high-richness legend helpers.
# ==============================================================================

PACKERA_TOL_OKABE_24 <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
  "#332288", "#882255", "#44AA99", "#999933", "#666666", "#117733",
  "#DDCC77", "#661100", "#6699CC", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#008080", "#88CCEE", "#AA4499"
)
PACKERA_TOL_OKABE_12 <- PACKERA_TOL_OKABE_24[1:12]
PACKERA_HOLLOW_SHAPES_5 <- c(1L, 2L, 5L, 0L, 6L)

PACKERA_ICN_ROMAN_TOKENS <- c(
  "subsp.", "subsp", "ssp.", "ssp", "var.", "var", "f.", "forma",
  "sect.", "sect", "ser.", "ser", "subgen.", "subgen",
  "x", "\u00d7", "aff.", "aff", "cf.", "cf", "sp.", "sp",
  "nov.", "nov", "comb.", "s.s.", "s.l.", "ined.", "auct.", "non", "p.p."
)

.esc_plotmath_str <- function(x) {
  s <- gsub("\\\\", "\\\\\\\\", as.character(x), fixed = TRUE)
  gsub("\n", "\\n", gsub('"', '\\"', s, fixed = TRUE), fixed = TRUE)
}

abbreviate_icn_taxon <- function(label_vec, abbreviate_genus = FALSE) {
  vapply(as.character(label_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("Unknown")
    s_clean <- trimws(s)
    if (!isTRUE(abbreviate_genus)) return(s_clean)
    if (grepl("^(Morphotype|Ecotype|Complex|Cluster|Deme)_", s_clean)) return(gsub("_", " ", s_clean))
    toks <- strsplit(s_clean, "[[:space:]]+")[[1L]]
    if (length(toks) >= 2L && grepl("^[A-Z][a-z]+$", toks[1L])) {
      toks[1L] <- paste0(substr(toks[1L], 1L, 1L), ".")
      return(paste(toks, collapse = " "))
    }
    s_clean
  }, character(1), USE.NAMES = FALSE)
}

.icn_taxon_plotmath_str <- function(s, bold_all = FALSE) {
  toks <- strsplit(trimws(as.character(s)), "[[:space:]]+")[[1L]]
  n_t <- length(toks)
  if (n_t == 0L) return("")
  wrap_it <- function(txt) sprintf('%s("%s")', if (isTRUE(bold_all)) "bolditalic" else "italic", .esc_plotmath_str(txt))
  wrap_rm <- function(txt) if (isTRUE(bold_all)) sprintf('bold("%s")', .esc_plotmath_str(txt)) else sprintf('"%s"', .esc_plotmath_str(txt))
  rank_resets <- c("subsp.", "subsp", "ssp.", "ssp", "var.", "var", "f.", "forma",
                   "sect.", "sect", "ser.", "ser", "subgen.", "subgen", "aff.", "aff", "cf.", "cf")
  parts <- character(0); in_auth <- FALSE; prev_was_x <- FALSE; i <- 1L
  while (i <= n_t) {
    tk <- toks[i]; low <- tolower(tk)
    nxt <- if (i < n_t) toks[i + 1L] else ""; nxt_low <- tolower(nxt)
    if (i == 1L) {
      parts <- c(parts, if (grepl("aceae$", tk, ignore.case = TRUE) || low %in% c("asteraceae", "compositae")) wrap_rm(tk) else wrap_it(tk))
      prev_was_x <- FALSE; i <- i + 1L; next
    }
    if (grepl("^\u00d7[a-z][a-z-]+$", tk)) {
      parts <- c(parts, sprintf('paste(%s, %s)', wrap_rm("\u00d7"), wrap_it(sub("^\u00d7", "", tk))))
      in_auth <- FALSE; prev_was_x <- FALSE; i <- i + 1L; next
    }
    if (low == "x" && grepl("^=", nxt)) {
      parts <- c(parts, wrap_rm("x")); in_auth <- TRUE; prev_was_x <- FALSE; i <- i + 1L; next
    }
    if (low %in% c("x", "\u00d7")) {
      if (nzchar(nxt) && !(nxt_low %in% PACKERA_ICN_ROMAN_TOKENS) && grepl("^[a-z][a-z-]+$", nxt)) {
        parts <- c(parts, sprintf('paste(%s, %s)', wrap_rm("\u00d7"), wrap_it(nxt)))
        in_auth <- FALSE; prev_was_x <- FALSE; i <- i + 2L; next
      }
      parts <- c(parts, wrap_rm("\u00d7")); in_auth <- FALSE; prev_was_x <- TRUE; i <- i + 1L; next
    }
    if (low %in% PACKERA_ICN_ROMAN_TOKENS) {
      rom <- if (low %in% c("subsp", "ssp", "subsp.", "ssp.")) "subsp."
             else if (low %in% c("var", "var.")) "var."
             else if (low %in% c("f", "f.", "forma")) "f." else tk
      parts <- c(parts, wrap_rm(rom))
      if (low %in% rank_resets) in_auth <- FALSE
      prev_was_x <- FALSE; i <- i + 1L; next
    }
    if (isTRUE(prev_was_x) && grepl("^[A-Z][a-z]*\\.?$", tk) && nzchar(nxt) &&
        !(nxt_low %in% PACKERA_ICN_ROMAN_TOKENS) && grepl("^[a-z][a-z-]+$", nxt)) {
      parts <- c(parts, wrap_it(tk)); in_auth <- FALSE; prev_was_x <- FALSE; i <- i + 1L; next
    }
    if (!in_auth && grepl("^[a-z][a-z-]+$", tk)) {
      parts <- c(parts, wrap_it(tk)); prev_was_x <- FALSE; i <- i + 1L; next
    }
    parts <- c(parts, wrap_rm(tk)); in_auth <- TRUE; prev_was_x <- FALSE; i <- i + 1L
  }
  paste(parts, collapse = " ~ ")
}

format_icn_label_expr <- function(tx, abbreviate_genus = FALSE, count_n = NULL) {
  s <- abbreviate_icn_taxon(tx, abbreviate_genus = abbreviate_genus)
  suffix_n <- if (!is.null(count_n) && is.finite(count_n)) sprintf(" (n=%d)", as.integer(count_n)) else ""
  if (grepl("^(Morphotype|Ecotype|Complex|Cluster|Deme|Unknown|Unassigned|Asteraceae)", s, ignore.case = TRUE) ||
      !grepl("^[A-Z][a-z.]*[[:space:]]+", s)) {
    lbl_plain <- paste0(gsub("_", " ", s), suffix_n)
    return(bquote(.(lbl_plain)))
  }
  code_str <- .icn_taxon_plotmath_str(s, bold_all = FALSE)
  if (nzchar(suffix_n)) code_str <- paste0(code_str, ' ~ "', .esc_plotmath_str(trimws(suffix_n)), '"')
  tryCatch(parse(text = code_str)[[1L]], error = function(e) paste0(s, suffix_n))
}

format_icn_title_expr <- function(prefix = "", taxon = "", suffix = "", bold_all = TRUE) {
  pre_s <- if (is.null(prefix) || length(prefix) == 0L || is.na(prefix[1L])) "" else as.character(prefix)[1L]
  suf_s <- if (is.null(suffix) || length(suffix) == 0L || is.na(suffix[1L])) "" else as.character(suffix)[1L]
  tx_s <- if (is.null(taxon) || length(taxon) == 0L || is.na(taxon[1L])) "" else trimws(as.character(taxon)[1L])
  if (!nzchar(tx_s) || grepl("^(Ecotype|Morphotype|Complex|Cluster|Deme)_", tx_s, ignore.case = TRUE) ||
      grepl("^(Unknown|Unassigned|Asteraceae)", tx_s, ignore.case = TRUE) || !grepl("^[A-Z][a-z.]*[[:space:]]+", tx_s)) {
    return(paste0(pre_s, gsub("_", " ", tx_s), suf_s))
  }
  wrap_txt <- function(txt) if (isTRUE(bold_all)) sprintf('bold("%s")', .esc_plotmath_str(txt)) else sprintf('"%s"', .esc_plotmath_str(txt))
  tx_code <- .icn_taxon_plotmath_str(tx_s, bold_all = isTRUE(bold_all))
  args_pm <- c(if (nzchar(pre_s)) wrap_txt(pre_s) else NULL, tx_code, if (nzchar(suf_s)) wrap_txt(suf_s) else NULL)
  full_code <- if (length(args_pm) == 1L) args_pm[1L] else sprintf("paste(%s)", paste(args_pm, collapse = ", "))
  tryCatch(parse(text = full_code)[[1L]], error = function(e) paste0(pre_s, tx_s, suf_s))
}

format_icn_label_vec <- function(taxa_vec, abbreviate_genus = FALSE, counts_map = NULL) {
  tx_chr <- as.character(taxa_vec)
  if (length(tx_chr) == 0L) return(expression())
  out_expr <- as.expression(lapply(tx_chr, function(tx) {
    n_v <- if (!is.null(counts_map) && !is.na(tx) && tx %in% names(counts_map)) counts_map[[tx]] else NULL
    format_icn_label_expr(tx, abbreviate_genus = abbreviate_genus, count_n = n_v)
  }))
  names(out_expr) <- tx_chr
  out_expr
}

icn_discrete_labeller <- function(abbreviate_genus = FALSE, counts_map = NULL) {
  force(abbreviate_genus); force(counts_map)
  function(x) unname(format_icn_label_vec(x, abbreviate_genus = abbreviate_genus, counts_map = counts_map))
}

get_combinatorial_palette_r <- function(groups_vec, preserve_order = FALSE) {
  lvls <- unique(as.character(stats::na.omit(groups_vec)))
  lvls <- lvls[nzchar(trimws(lvls))]
  if (!isTRUE(preserve_order)) lvls <- sort(lvls)
  if (length(lvls) == 0L) return(character(0))
  n_c <- length(PACKERA_TOL_OKABE_24)
  stats::setNames(PACKERA_TOL_OKABE_24[((seq_along(lvls) - 1L) %% n_c) + 1L], lvls)
}

get_combinatorial_hollow_shapes_r <- function(groups_vec, preserve_order = FALSE) {
  lvls <- unique(as.character(stats::na.omit(groups_vec)))
  lvls <- lvls[nzchar(trimws(lvls))]
  if (!isTRUE(preserve_order)) lvls <- sort(lvls)
  if (length(lvls) == 0L) return(integer(0))
  n_c <- length(PACKERA_TOL_OKABE_24); n_s <- length(PACKERA_HOLLOW_SHAPES_5)
  shps <- if (length(lvls) > 8L) PACKERA_HOLLOW_SHAPES_5[(((seq_along(lvls) - 1L) %/% n_c) %% n_s) + 1L] else PACKERA_HOLLOW_SHAPES_5[((seq_along(lvls) - 1L) %% n_s) + 1L]
  stats::setNames(shps, lvls)
}

classify_type_plot_tier_r <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(character(0))
  n <- nrow(df)
  truthy_col <- function(col) {
    if (!(col %in% names(df))) return(rep(FALSE, n))
    v <- df[[col]]
    if (is.logical(v)) return(!is.na(v) & v)
    if (is.numeric(v)) return(is.finite(v) & v != 0)
    tolower(trimws(as.character(v))) %in% c("true", "t", "1", "1.0", "yes")
  }
  chr_col <- function(col) {
    if (!(col %in% names(df))) return(rep("", n))
    v <- trimws(as.character(df[[col]])); v[is.na(v)] <- ""; v
  }
  is_tp <- truthy_col("is_type_specimen"); is_cent <- truthy_col("is_centroid_or_inst_coord")
  anc <- chr_col("anchor_tier_source"); subt <- chr_col("taxonomic_treatment_subtype")
  stat <- chr_col("taxonomic_status_call"); tns <- chr_col("type_niche_status")
  any_type <- is_tp | (anc == "Tier_0_Type_Anchor") | (subt %in% c("Confirmed_Type_Control", "Type_Synonymy_Candidate")) | (tns == "Passive_Type_Uncertain_Coord")
  out <- rep("None", n)
  if (!any(any_type)) return(out)
  is_unc <- any_type & ((tns == "Passive_Type_Uncertain_Coord") | is_cent)
  is_syn <- any_type & !is_unc & ((subt == "Type_Synonymy_Candidate") | (stat == "Type_Synonymy_Candidate") |
    grepl("Synonymy|Overlap", subt, ignore.case = TRUE) | grepl("Synonymy|Overlap", stat, ignore.case = TRUE))
  out[any_type] <- "Confirmed_Type_Control"; out[is_syn] <- "Type_Synonymy_Candidate"; out[is_unc] <- "Passive_Type_Uncertain_Coord"
  out
}

add_nomenclatural_type_layer_gg <- function(p, df, x_col, y_col, show_uncertain_coord = FALSE) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L || !all(c(x_col, y_col) %in% names(df))) return(p)
  xv <- suppressWarnings(as.numeric(df[[x_col]])); yv <- suppressWarnings(as.numeric(df[[y_col]]))
  tiers <- classify_type_plot_tier_r(df); ok <- is.finite(xv) & is.finite(yv) & tiers != "None"
  if (!any(ok)) return(p)
  tdf <- data.frame(tx = xv[ok], ty = yv[ok], tier = tiers[ok], stringsAsFactors = FALSE)
  c_df <- tdf[tdf$tier == "Confirmed_Type_Control", , drop = FALSE]
  s_df <- tdf[tdf$tier == "Type_Synonymy_Candidate", , drop = FALSE]
  u_df <- tdf[tdf$tier == "Passive_Type_Uncertain_Coord", , drop = FALSE]
  if (nrow(c_df) > 0L) {
    p <- p +
      ggplot2::geom_point(data = c_df, ggplot2::aes(x = .data$tx, y = .data$ty), shape = 23, fill = NA, size = 4.4, stroke = 1.1, color = "#332288", inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::geom_point(data = c_df, ggplot2::aes(x = .data$tx, y = .data$ty), shape = 5, size = 2.3, stroke = 0.9, color = "#111827", inherit.aes = FALSE, show.legend = FALSE)
  }
  if (nrow(s_df) > 0L) {
    p <- p +
      ggplot2::geom_point(data = s_df, ggplot2::aes(x = .data$tx, y = .data$ty), shape = 24, fill = NA, size = 4.2, stroke = 1.1, color = "#882255", inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::geom_point(data = s_df, ggplot2::aes(x = .data$tx, y = .data$ty), shape = 2, size = 2.1, stroke = 0.9, color = "#882255", inherit.aes = FALSE, show.legend = FALSE)
  }
  if (isTRUE(show_uncertain_coord) && nrow(u_df) > 0L) {
    p <- p +
      ggplot2::geom_point(data = u_df, ggplot2::aes(x = .data$tx, y = .data$ty), shape = 1, size = 4.5, stroke = 0.95, color = "#64748B", inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::geom_point(data = u_df, ggplot2::aes(x = .data$tx, y = .data$ty), shape = 5, size = 2.3, stroke = 0.8, color = "#64748B", inherit.aes = FALSE, show.legend = FALSE)
  }
  p
}

apply_high_richness_legend_gg <- function(p, n_groups, cutoff = 12L, ncol_bottom = 4L, hide_above = 18L) {
  ng <- suppressWarnings(as.integer(n_groups[1L]))
  if (is.na(ng) || ng <= cutoff) return(p)
  if (ng > hide_above) return(p + ggplot2::theme(legend.position = "none"))
  p + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 7.2), legend.key.size = grid::unit(0.32, "cm")) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = ncol_bottom, byrow = TRUE), shape = ggplot2::guide_legend(ncol = ncol_bottom, byrow = TRUE))
}
