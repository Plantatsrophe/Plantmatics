if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# [WORKER1_MODULAR_FACADE]
# Module: scripts/figures/theme_packera_pub.R (300-DPI PNG theme, hollow scales & ICN labels)

if (!exists("env_plantmatics_or_packera", mode = "function")) {
  env_plantmatics_or_packera <- function(p_key, k_key = NULL, unset = "") {
    if (is.null(k_key) || !nzchar(k_key)) {
      k_key <- sub("^PLANTMATICS_", "PACKERA_", p_key)
    }
    v <- Sys.getenv(p_key, unset = "")
    if (nzchar(trimws(v))) return(v)
    v2 <- Sys.getenv(k_key, unset = "")
    if (nzchar(trimws(v2))) return(v2)
    unset
  }
}

normalize_png_path <- function(out_path) {
  if (is.null(out_path) || !nzchar(trimws(as.character(out_path)))) return("")
  p <- trimws(as.character(out_path))
  if (grepl("\\.pdf$", p, ignore.case = TRUE)) {
    p <- sub("\\.pdf$", ".png", p, ignore.case = TRUE)
  } else if (!grepl("\\.png$", p, ignore.case = TRUE)) {
    p <- paste0(p, ".png")
  }
  p
}

save_packera_png_300 <- function(plot_or_fn, out_path, width = 12.5, height = 8.0, dpi = 300L, also_copy_to = NULL) {
  png_final <- normalize_png_path(out_path)
  if (!nzchar(png_final)) return(invisible(""))
  if (identical(trimws(env_plantmatics_or_packera("PLANTMATICS_SKIP_FIGURES", "PACKERA_SKIP_FIGURES")), "1")) {
    return(invisible(png_final))
  }
  fdir <- dirname(png_final)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  stem <- tools::file_path_sans_ext(basename(png_final))
  png_tmp <- file.path(fdir, sprintf(".%s_%d.png.tmp", stem, Sys.getpid()))
  stale_pdf <- file.path(fdir, paste0(stem, ".pdf"))

  open_png_device <- function(target_file) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(filename = target_file, width = width, height = height, units = "in", res = dpi, background = "white")
      return(TRUE)
    }
    ok_cairo <- tryCatch({
      grDevices::png(filename = target_file, width = width, height = height, units = "in", res = dpi, bg = "white", type = "cairo")
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(ok_cairo)) {
      grDevices::png(filename = target_file, width = width, height = height, units = "in", res = dpi, bg = "white")
    }
    TRUE
  }

  prev_dev <- grDevices::dev.cur()
  opened_dev <- NA_integer_
  ok_render <- tryCatch({
    open_png_device(png_tmp)
    opened_dev <- grDevices::dev.cur()
    if (is.function(plot_or_fn)) {
      res_p <- plot_or_fn()
      if (inherits(res_p, "gg") || inherits(res_p, "ggplot") || inherits(res_p, "patchwork")) {
        print(res_p)
      }
    } else {
      print(plot_or_fn)
    }
    TRUE
  }, error = function(e) {
    warning("save_packera_png_300 warning (", basename(png_final), "): ", e$message)
    FALSE
  }, finally = {
    if (is.finite(opened_dev) && opened_dev > 1L && opened_dev != prev_dev && opened_dev %in% grDevices::dev.list()) {
      tryCatch(grDevices::dev.off(opened_dev), error = function(e) NULL)
    }
    if (prev_dev > 1L && prev_dev %in% grDevices::dev.list()) {
      tryCatch(grDevices::dev.set(prev_dev), error = function(e) NULL)
    }
  })

  if (isTRUE(ok_render) && file.exists(png_tmp) && file.info(png_tmp)$size > 0) {
    file.rename(png_tmp, png_final)
    if (file.exists(stale_pdf)) unlink(stale_pdf, force = TRUE)
    if (!is.null(also_copy_to)) {
      for (cp_raw in as.character(also_copy_to)) {
        cp_png <- normalize_png_path(cp_raw)
        if (nzchar(cp_png) && normalizePath(cp_png, mustWork = FALSE) != normalizePath(png_final, mustWork = FALSE)) {
          dir.create(dirname(cp_png), recursive = TRUE, showWarnings = FALSE)
          cp_tmp <- paste0(cp_png, ".", Sys.getpid(), ".tmp")
          if (isTRUE(file.copy(png_final, cp_tmp, overwrite = TRUE))) file.rename(cp_tmp, cp_png)
          cp_stale_pdf <- sub("\\.png$", ".pdf", cp_png, ignore.case = TRUE)
          if (cp_stale_pdf != cp_png && file.exists(cp_stale_pdf)) unlink(cp_stale_pdf, force = TRUE)
        }
      }
    }
  } else if (file.exists(png_tmp)) {
    unlink(png_tmp, force = TRUE)
  }
  invisible(png_final)
}

resolve_taxon_or_rank_vector <- function(df, group_col = NULL) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  base_v <- if (!is.null(group_col) && group_col %in% names(df)) {
    as.character(df[[group_col]])
  } else if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("scientificName" %in% names(df)) {
    as.character(df$scientificName)
  } else {
    rep("Unknown", n)
  }
  base_v[is.na(base_v) | !nzchar(trimws(base_v))] <- "Unknown"
  base_v <- trimws(base_v)

  # Append infraspecific rank (var. / subsp. / f.) ONLY when explicitly requested via group_col = "infraspecific_taxon"
  if ("infraspecific_taxon" %in% names(df) && identical(group_col, "infraspecific_taxon")) {
    inf_v <- trimws(as.character(df$infraspecific_taxon))
    has_inf <- !is.na(inf_v) & nzchar(inf_v) & !(tolower(inf_v) %in% c("na", "none", "null", "unknown"))
    if (any(has_inf)) {
      for (i in which(has_inf)) {
        b_i <- base_v[i]
        if (!grepl("^Packera[[:space:]]+", b_i) || grepl("(^|[[:space:]])(subsp\\.|var\\.|f\\.)[[:space:]]+", b_i, ignore.case = TRUE)) next
        iv <- inf_v[i]
        if (grepl("^(subsp\\.?|ssp\\.?|var\\.?|f\\.?|forma)[[:space:]]+", iv, ignore.case = TRUE)) {
          iv_clean <- sub("^(subsp|ssp)\\.?[[:space:]]+", "subsp. ", iv, ignore.case = TRUE)
          iv_clean <- sub("^var\\.?[[:space:]]+", "var. ", iv_clean, ignore.case = TRUE)
          iv_clean <- sub("^(f|forma)\\.?[[:space:]]+", "f. ", iv_clean, ignore.case = TRUE)
          base_v[i] <- paste(b_i, iv_clean)
        } else if (grepl("^[a-z][a-z-]+$", tolower(iv))) {
          ep <- strsplit(b_i, "[[:space:]]+")[[1L]]
          if (length(ep) >= 2L && tolower(iv) != tolower(ep[2L])) {
            base_v[i] <- paste0(b_i, " var. ", tolower(iv))
          }
        }
      }
    }
  }
  base_v
}

resolve_plot_group_labels <- function(df, cda_res = NULL, gmm_res = NULL, group_col = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  if (!is.null(group_col) && is.character(group_col) && length(group_col) == 1L && group_col %in% names(df)) {
    return(resolve_taxon_or_rank_vector(df, group_col = group_col))
  }
  if (exists("resolve_plot_group_column", mode = "function")) {
    gv <- tryCatch(resolve_plot_group_column(df, cda_res = cda_res, gmm_res = gmm_res), error = function(e) NULL)
    if (is.list(gv) && !is.null(gv$col)) {
      return(resolve_taxon_or_rank_vector(df, group_col = gv$col))
    }
    if (is.character(gv) && length(gv) == nrow(df)) {
      gv_clean <- trimws(gv)
      gv_clean[is.na(gv_clean) | !nzchar(gv_clean)] <- "Unknown"
      non_unk <- gv_clean[gv_clean != "Unknown"]
      if (length(non_unk) > 0L && all(grepl("^Packera[[:space:]]+", non_unk)) && "infraspecific_taxon" %in% names(df)) {
        tmp_df <- df
        tmp_df$.__tmp_base_sp <- gv_clean
        return(resolve_taxon_or_rank_vector(tmp_df, group_col = ".__tmp_base_sp"))
      }
      return(gv_clean)
    }
  }
  resolve_taxon_or_rank_vector(df)
}

theme_packera_pub <- function(base_size = 10.2) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.8, hjust = 0, color = "#111827"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1.1, hjust = 0, color = "#4B5563"),
      axis.title = ggplot2::element_text(size = base_size - 0.3, color = "#1F2937"),
      axis.text = ggplot2::element_text(size = base_size - 1.1, color = "#1F2937"),
      panel.border = ggplot2::element_rect(color = "#9CA3AF", fill = NA, linewidth = 0.65),
      panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35, linetype = "dotted"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", color = "#D1D5DB", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.8, color = "#111827"),
      legend.background = ggplot2::element_rect(fill = "white", color = "#D1D5DB", linewidth = 0.4),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1.0),
      legend.text = ggplot2::element_text(size = base_size - 1.3),
      plot.margin = ggplot2::margin(6, 8, 6, 8, unit = "pt")
    )
}

wrap_clean_label_r <- function(txt_vec, width = 28L) {
  vapply(as.character(txt_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("")
    paste(strwrap(trimws(s), width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

add_centroid_repel_labels_gg <- function(p, cent_df, x_col, y_col, label_col, color_col = NULL, size = 2.85, allow_ggrepel = TRUE) {
  if (is.null(cent_df) || nrow(cent_df) == 0L) return(p)
  xc <- as.character(x_col)[1L]
  yc <- as.character(y_col)[1L]
  lc <- as.character(label_col)[1L]
  if (!all(c(xc, yc, lc) %in% names(cent_df))) return(p)
  xv <- suppressWarnings(as.numeric(cent_df[[xc]]))
  yv <- suppressWarnings(as.numeric(cent_df[[yc]]))
  lv <- as.character(cent_df[[lc]])
  has_col <- !is.null(color_col) && length(color_col) > 0L && nzchar(as.character(color_col)[1L]) && (as.character(color_col)[1L] %in% names(cent_df))
  cv <- if (has_col) as.character(cent_df[[as.character(color_col)[1L]]]) else rep("#111827", nrow(cent_df))
  ok <- is.finite(xv) & is.finite(yv) & !is.na(lv) & nzchar(trimws(lv))
  if (!any(ok)) return(p)
  rep_df <- data.frame(rx = xv[ok], ry = yv[ok], rlbl = lv[ok], rcol = cv[ok], stringsAsFactors = FALSE)
  is_latin <- all(grepl("^[A-Z][a-z.]*[[:space:]]+[a-z]", rep_df$rlbl))
  face_use <- if (is_latin) "italic" else "plain"
  sz_num <- as.numeric(size)[1L]
  if (isTRUE(allow_ggrepel) && requireNamespace("ggrepel", quietly = TRUE)) {
    rep_fmls <- names(formals(ggrepel::geom_label_repel))
    base_args <- list(
      data = rep_df,
      fill = grDevices::adjustcolor("white", alpha.f = 0.88),
      size = sz_num,
      fontface = face_use,
      show.legend = FALSE,
      inherit.aes = FALSE
    )
    if ("max.overlaps" %in% rep_fmls) base_args$max.overlaps <- 50
    if ("box.padding" %in% rep_fmls) base_args$box.padding <- 0.28
    if (has_col) {
      base_args$mapping <- ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl, color = .data$rcol)
    } else {
      base_args$mapping <- ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl)
      base_args$color <- "#111827"
    }
    lyr <- tryCatch(do.call(ggrepel::geom_label_repel, base_args), error = function(e) NULL)
    if (!is.null(lyr)) return(p + lyr)
  }
  if (has_col) {
    return(p + ggplot2::geom_label(
      data = rep_df,
      mapping = ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl, color = .data$rcol),
      fill = grDevices::adjustcolor("white", alpha.f = 0.86),
      size = sz_num,
      fontface = face_use,
      show.legend = FALSE,
      inherit.aes = FALSE
    ))
  }
  p + ggplot2::geom_label(
    data = rep_df,
    mapping = ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl),
    fill = grDevices::adjustcolor("white", alpha.f = 0.86),
    color = "#111827",
    size = sz_num,
    fontface = face_use,
    show.legend = FALSE,
    inherit.aes = FALSE
  )
}

.resolve_self_dir <- function(default_rel, probe_file) {
  f_frm <- as.character(unlist(lapply(sys.frames(), function(fr) fr$ofile)))
  f_cmd <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand_d <- unique(c(rev(dirname(f_frm)), dirname(f_cmd), default_rel))
  cand_d <- cand_d[!is.na(cand_d) & nzchar(cand_d)]
  for (d in cand_d) {
    if (file.exists(file.path(d, probe_file))) return(d)
  }
  default_rel
}
.load_env <- environment()
.fig_self_dir <- .resolve_self_dir("scripts/figures", "icn_palette_utils.R")
source(file.path(.fig_self_dir, "icn_palette_utils.R"), local = .load_env)
source(file.path(.resolve_self_dir("scripts/figures", "basemap_exemplar_utils.R"), "basemap_exemplar_utils.R"), local = .load_env)
