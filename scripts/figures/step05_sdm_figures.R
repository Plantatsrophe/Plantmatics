if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/figures/step05_sdm_figures.R
# Step 05 Equal-Area Projected MaxEnt/Ensemble SDM Suitability Figures:
#   fig05c_sdm_suitability_and_distribution_map.png
# ==============================================================================

.load_na_sf_5070 <- function() {
  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("maps", quietly = TRUE)) return(list(land = NULL, states = NULL))
  old_s2 <- tryCatch(sf::sf_use_s2(FALSE), error = function(e) NULL)
  on.exit(if (!is.null(old_s2)) tryCatch(sf::sf_use_s2(old_s2), error = function(e) NULL), add = TRUE)
  mk_sf <- function(db, ...) tryCatch({
    s <- sf::st_as_sf(maps::map(db, plot = FALSE, fill = TRUE, ...))
    s <- if (is.na(sf::st_crs(s))) sf::st_set_crs(s, 4326) else sf::st_transform(s, 4326)
    s <- suppressWarnings(sf::st_crop(sf::st_make_valid(s), c(xmin = -170, ymin = 18, xmax = -50, ymax = 75)))
    sf::st_transform(s, 5070L)
  }, error = function(e) NULL)
  list(land = mk_sf("world", regions = c("USA", "Canada")), states = mk_sf("state"))
}

.sdm_filter_occ <- function(d_in) {
  if (exists("filter_map_occurrences_r", mode = "function")) {
    return(filter_map_occurrences_r(d_in, lon_col = "decimalLongitude", lat_col = "decimalLatitude"))
  }
  list(clean_df = d_in, uncertain_type_df = d_in[0L, , drop = FALSE], n_quarantined_non_type = 0L)
}

.sdm_proj_xy <- function(d_in, projected = TRUE) {
  if (is.null(d_in) || nrow(d_in) == 0L) return(data.frame())
  lon <- suppressWarnings(as.numeric(d_in$decimalLongitude))
  lat <- suppressWarnings(as.numeric(d_in$decimalLatitude))
  ok <- is.finite(lon) & is.finite(lat) & abs(lat) <= 90 & abs(lon) <= 180
  if (!any(ok)) return(data.frame())
  sub <- d_in[ok, , drop = FALSE]; sub$lon_num <- lon[ok]; sub$lat_num <- lat[ok]
  sub$px <- sub$lon_num; sub$py <- sub$lat_num
  if (isTRUE(projected) && requireNamespace("sf", quietly = TRUE)) {
    xy <- tryCatch(sf::st_coordinates(sf::st_transform(sf::st_as_sf(sub, coords = c("lon_num", "lat_num"), crs = 4326), 5070L)), error = function(e) NULL)
    if (!is.null(xy) && nrow(xy) == nrow(sub)) { sub$px <- as.numeric(xy[, 1L]); sub$py <- as.numeric(xy[, 2L]) }
  }
  sub
}

.sdm_overlay_types <- function(p, clean_p, unc_p, x_col = "px", y_col = "py", facet_col = NULL) {
  cm <- intersect(names(clean_p), names(unc_p))
  t_all <- if (nrow(clean_p) > 0L && nrow(unc_p) > 0L) rbind(clean_p[, cm, drop = FALSE], unc_p[, cm, drop = FALSE]) else if (nrow(unc_p) > 0L) unc_p else clean_p
  if (nrow(t_all) == 0L || !all(c(x_col, y_col) %in% names(t_all))) return(p)
  if (is.null(facet_col)) {
    if (!exists("add_nomenclatural_type_layer_gg", mode = "function")) return(p)
    return(add_nomenclatural_type_layer_gg(p, t_all, x_col, y_col, show_uncertain_coord = TRUE))
  }
  if (!exists("classify_type_plot_tier_r", mode = "function") || !(facet_col %in% names(t_all))) return(p)
  t_all$.tier <- classify_type_plot_tier_r(t_all)
  t_all <- t_all[t_all$.tier != "None" & is.finite(suppressWarnings(as.numeric(t_all[[x_col]]))) & is.finite(suppressWarnings(as.numeric(t_all[[y_col]]))), , drop = FALSE]
  if (nrow(t_all) == 0L) return(p)
  for (sp in list(
    list(tr = "Confirmed_Type_Control", s1 = 23, f1 = NA, z1 = 3.2, c1 = "#332288", s2 = 5, z2 = 1.8, c2 = "#111827"),
    list(tr = "Type_Synonymy_Candidate", s1 = 24, f1 = NA, z1 = 3.1, c1 = "#882255", s2 = 2, z2 = 1.7, c2 = "#882255"),
    list(tr = "Passive_Type_Uncertain_Coord", s1 = 1, f1 = NA, z1 = 3.3, c1 = "#64748B", s2 = 5, z2 = 1.7, c2 = "#64748B")
  )) {
    sd_t <- t_all[t_all$.tier == sp$tr, , drop = FALSE]
    if (nrow(sd_t) > 0L) {
      p <- p + ggplot2::geom_point(data = sd_t, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]), shape = sp$s1, fill = sp$f1, size = sp$z1, stroke = 0.9, color = sp$c1, inherit.aes = FALSE, show.legend = FALSE) +
        ggplot2::geom_point(data = sd_t, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]), shape = sp$s2, size = sp$z2, stroke = 0.75, color = sp$c2, inherit.aes = FALSE, show.legend = FALSE)
    }
  }
  p
}

.sdm_raster_to_proj_df <- function(r_obj, max_cells = 45000L, target_crs = "EPSG:5070") {
  if (is.null(r_obj) || !requireNamespace("terra", quietly = TRUE)) return(list(df = data.frame(), projected = FALSE))
  tryCatch({
    r_use <- if (terra::ncell(r_obj) > max_cells) terra::aggregate(r_obj, fact = max(2L, ceiling(sqrt(terra::ncell(r_obj) / max_cells))), fun = "mean", na.rm = TRUE) else r_obj
    is_proj <- FALSE; r_out <- r_use
    if (requireNamespace("sf", quietly = TRUE)) {
      if (!nzchar(terra::crs(r_use))) terra::crs(r_use) <- "EPSG:4326"
      r_tmp <- tryCatch(terra::project(r_use, target_crs, method = "bilinear"), error = function(e) NULL)
      if (!is.null(r_tmp) && terra::ncell(r_tmp) > 0L) {
        d_chk <- as.data.frame(r_tmp, xy = TRUE, na.rm = TRUE)
        if (ncol(d_chk) >= 3L && nrow(d_chk) > 0L && any(is.finite(d_chk[[3L]]))) { r_out <- r_tmp; is_proj <- TRUE }
      }
    }
    d_r <- as.data.frame(r_out, xy = TRUE, na.rm = TRUE)
    if (ncol(d_r) < 3L || nrow(d_r) == 0L) return(list(df = data.frame(), projected = FALSE))
    names(d_r)[1:3] <- c("x", "y", "suitability")
    list(df = d_r[is.finite(d_r$suitability), c("x", "y", "suitability"), drop = FALSE], projected = is_proj)
  }, error = function(e) list(df = data.frame(), projected = FALSE))
}

render_step05c_sdm_suitability_png <- function(
  df, focal_rasters = NULL, multi_rasters = NULL, perf_df = NULL, resp_df = NULL,
  roc_df = NULL, boyce_df = NULL, figures_dir = "outputs/figures", focal_taxon = "", out_sdm_path = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(""))
  fdir <- as.character(figures_dir); dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  vdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(vdf) == 0L) return(invisible(""))
  if (!("decimalLongitude" %in% names(vdf)) && "longitude" %in% names(vdf)) vdf$decimalLongitude <- vdf$longitude
  if (!("decimalLatitude" %in% names(vdf)) && "latitude" %in% names(vdf)) vdf$decimalLatitude <- vdf$latitude
  vdf$plot_group <- resolve_plot_group_labels(vdf)
  active_groups <- sort(unique(vdf$plot_group[!is.na(vdf$plot_group) & nzchar(vdf$plot_group)]))
  pal_map <- get_combinatorial_palette_r(active_groups); shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = TRUE, counts_map = table(vdf$plot_group))
  names(lbl_exprs) <- active_groups

  na_sf <- .load_na_sf_5070()
  water_bg <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#EEF4F8", color = NA))
  .app_q <- function(s, nq) if (is.finite(nq) && nq > 0L) paste0(s, sprintf(" | %d centroid/inst. non-types excluded", as.integer(nq))) else s

  build_suitability_panel <- function(r_info, pts_sub, title_expr, sub_txt = "") {
    r_df <- if (is.list(r_info) && is.data.frame(r_info$df)) r_info$df else data.frame()
    if (nrow(r_df) > 0L) r_df <- r_df[is.finite(r_df$x) & is.finite(r_df$y) & is.finite(r_df$suitability), , drop = FALSE]
    occ_s <- .sdm_filter_occ(pts_sub); sub_u <- .app_q(sub_txt, occ_s$n_quarantined_non_type)
    is_proj <- is.list(r_info) && isTRUE(r_info$projected) && !is.null(na_sf$states)
    if (nrow(r_df) > 0L) {
      clean_p <- .sdm_proj_xy(occ_s$clean_df, projected = is_proj)
      unc_p <- .sdm_proj_xy(occ_s$uncertain_type_df, projected = is_proj)
      if (is_proj && requireNamespace("sf", quietly = TRUE)) {
        p_m <- ggplot2::ggplot()
        if (!is.null(na_sf$land)) p_m <- p_m + ggplot2::geom_sf(data = na_sf$land, fill = "#F3F4F6", color = "#CBD5E1", linewidth = 0.25, inherit.aes = FALSE)
        p_m <- p_m +
          ggplot2::geom_tile(data = r_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$suitability), alpha = 0.92, na.rm = TRUE) +
          ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Suitability") +
          ggplot2::geom_sf(data = na_sf$states, fill = NA, color = "#E5E7EB", linewidth = 0.26, inherit.aes = FALSE)
        if (!is.null(na_sf$land)) p_m <- p_m + ggplot2::geom_sf(data = na_sf$land, fill = NA, color = "#94A3B8", linewidth = 0.24, inherit.aes = FALSE)
        if (nrow(clean_p) > 0L) p_m <- p_m + ggplot2::geom_point(data = clean_p, ggplot2::aes(x = .data$px, y = .data$py), shape = 1, size = 2.0, stroke = 0.8, color = "#FF6B00", alpha = 0.88, inherit.aes = FALSE)
        p_m <- .sdm_overlay_types(p_m, clean_p, unc_p, "px", "py")
        bb_r <- c(xmin = min(r_df$x), xmax = max(r_df$x), ymin = min(r_df$y), ymax = max(r_df$y))
        if (exists("add_projected_scale_north_gg", mode = "function")) p_m <- add_projected_scale_north_gg(p_m, bb_r)
        return(p_m + ggplot2::coord_sf(crs = 5070L, xlim = range(r_df$x), ylim = range(r_df$y), expand = FALSE) + theme_packera_pub(base_size = 9.5) + water_bg + ggplot2::labs(title = title_expr, subtitle = sub_u, x = NULL, y = NULL))
      }
      cos_l <- cos(mean(range(r_df$y, na.rm = TRUE)) * pi / 180.0)
      p_m <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = r_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$suitability), alpha = 0.90, na.rm = TRUE) +
        ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Suitability")
      if (requireNamespace("maps", quietly = TRUE)) {
        st_d <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
        if (!is.null(st_d)) p_m <- p_m + ggplot2::geom_path(data = st_d, ggplot2::aes(x = .data$long, y = .data$lat, group = .data$group), color = "#E5E7EB", linewidth = 0.25)
      }
      if (nrow(clean_p) > 0L) p_m <- p_m + ggplot2::geom_point(data = clean_p, ggplot2::aes(x = .data$px, y = .data$py), shape = 1, size = 2.0, stroke = 0.8, color = "#FF6B00", alpha = 0.88)
      p_m <- .sdm_overlay_types(p_m, clean_p, unc_p, "px", "py")
      return(p_m + ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_l), xlim = range(r_df$x), ylim = range(r_df$y), expand = FALSE) + theme_packera_pub(base_size = 9.5) + water_bg + ggplot2::labs(title = title_expr, subtitle = sub_u, x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)"))
    }
    m_ctx <- build_na_projected_basemap_gg(occ_s$clean_df$decimalLongitude, occ_s$clean_df$decimalLatitude)
    clean_p <- project_lonlat_points_r(occ_s$clean_df, m_ctx)
    unc_p <- project_lonlat_points_r(occ_s$uncertain_type_df, m_ctx)
    p_m <- m_ctx$plot
    if (nrow(clean_p) > 0L) {
      p_m <- p_m + ggplot2::geom_point(data = clean_p, ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$plot_group, shape = .data$plot_group), size = 2.2, stroke = 0.85, alpha = 0.84) +
        ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
        ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group")
    }
    p_m <- .sdm_overlay_types(p_m, clean_p, unc_p, "map_x", "map_y")
    if (exists("apply_high_richness_legend_gg", mode = "function")) p_m <- apply_high_richness_legend_gg(p_m, length(active_groups))
    else if (length(active_groups) > 10L) p_m <- p_m + ggplot2::theme(legend.position = "none")
    p_m + ggplot2::labs(title = title_expr, subtitle = sub_u, x = if (isTRUE(m_ctx$projected)) NULL else "Longitude (\u00b0E)", y = if (isTRUE(m_ctx$projected)) NULL else "Latitude (\u00b0N)")
  }

  has_multi_sdm <- is.list(multi_rasters) && length(multi_rasters) >= 2L
  has_focal_sdm <- is.list(focal_rasters) && (!is.null(focal_rasters$maxent) || !is.null(focal_rasters$ensemble))
  p_multi_grid <- NULL; n_rows_multi <- 0L
  if (has_multi_sdm) {
    .fmt_sdm_strip <- function(v) wrap_clean_label_r(gsub("_", " ", abbreviate_icn_taxon(v, abbreviate_genus = TRUE)), width = 22L)
    strip_face_m <- if (any(grepl("^(Ecotype|Morphotype|Complex|Cluster)", names(multi_rasters)))) "bold" else "bold.italic"
    grid_infos <- lapply(names(multi_rasters), function(tx_nm) { ri <- .sdm_raster_to_proj_df(multi_rasters[[tx_nm]], max_cells = 18000L); if (nrow(ri$df) > 0L) ri$df$taxon_lbl <- .fmt_sdm_strip(tx_nm); ri })
    proj_m <- all(vapply(grid_infos, function(ri) nrow(ri$df) == 0L || isTRUE(ri$projected), logical(1))) && !is.null(na_sf$states)
    grid_all <- do.call(rbind, lapply(grid_infos, `[[`, "df"))
    if (!is.null(grid_all) && nrow(grid_all) > 0L) grid_all <- grid_all[is.finite(grid_all$x) & is.finite(grid_all$y) & is.finite(grid_all$suitability), , drop = FALSE]
    if (!is.null(grid_all) && nrow(grid_all) > 0L) {
      occ_all <- .sdm_filter_occ(vdf); sub_e_txt <- .app_q("Orange hollow circles mark coordinate-uncertainty-gated herbarium presences", occ_all$n_quarantined_non_type)
      cln_e <- .sdm_proj_xy(occ_all$clean_df, projected = proj_m); unc_e <- .sdm_proj_xy(occ_all$uncertain_type_df, projected = proj_m)
      if (nrow(cln_e) > 0L) { cln_e$taxon_lbl <- .fmt_sdm_strip(cln_e$plot_group); cln_e <- cln_e[cln_e$taxon_lbl %in% unique(grid_all$taxon_lbl), , drop = FALSE] }
      if (nrow(unc_e) > 0L) { unc_e$taxon_lbl <- .fmt_sdm_strip(unc_e$plot_group); unc_e <- unc_e[unc_e$taxon_lbl %in% unique(grid_all$taxon_lbl), , drop = FALSE] }
      n_sp_m <- length(unique(grid_all$taxon_lbl)); ncol_m <- min(4L, max(2L, ceiling(sqrt(n_sp_m)))); n_rows_multi <- as.integer(ceiling(n_sp_m / ncol_m))
      p_mg <- ggplot2::ggplot(grid_all, ggplot2::aes(x = .data$x, y = .data$y))
      if (proj_m && !is.null(na_sf$land)) p_mg <- p_mg + ggplot2::geom_sf(data = na_sf$land, fill = "#F3F4F6", color = "#CBD5E1", linewidth = 0.20, inherit.aes = FALSE)
      p_mg <- p_mg + ggplot2::geom_tile(ggplot2::aes(fill = .data$suitability), na.rm = TRUE)
      if (proj_m) {
        p_mg <- p_mg + ggplot2::geom_sf(data = na_sf$states, fill = NA, color = "#E5E7EB", linewidth = 0.22, inherit.aes = FALSE)
        if (!is.null(na_sf$land)) p_mg <- p_mg + ggplot2::geom_sf(data = na_sf$land, fill = NA, color = "#94A3B8", linewidth = 0.20, inherit.aes = FALSE)
      }
      if (nrow(cln_e) > 0L) p_mg <- p_mg + ggplot2::geom_point(data = cln_e, ggplot2::aes(x = .data$px, y = .data$py), shape = 1, size = 1.35, stroke = 0.65, color = "#FF6B00", alpha = 0.80, inherit.aes = FALSE)
      p_mg <- .sdm_overlay_types(p_mg, cln_e, unc_e, "px", "py", facet_col = "taxon_lbl") +
        ggplot2::facet_wrap(~ taxon_lbl, ncol = ncol_m) +
        ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Ensemble\nSuitability")
      p_multi_grid <- if (proj_m) {
        p_mg + ggplot2::coord_sf(crs = 5070L, xlim = range(grid_all$x), ylim = range(grid_all$y), expand = FALSE) + theme_packera_pub(base_size = 9.2) + water_bg +
          ggplot2::theme(strip.text = ggplot2::element_text(face = strip_face_m, size = 7.6, lineheight = 0.92)) +
          ggplot2::labs(title = sprintf("E. Comparative Multi-Taxon Ensemble Habitat Suitability (%d Modeled Taxa, EPSG:5070)", n_sp_m), subtitle = sub_e_txt, x = NULL, y = NULL)
      } else {
        cos_m <- cos(mean(range(grid_all$y, na.rm = TRUE)) * pi / 180.0)
        p_mg + ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_m), expand = FALSE) + theme_packera_pub(base_size = 9.2) + water_bg +
          ggplot2::theme(strip.text = ggplot2::element_text(face = strip_face_m, size = 7.6, lineheight = 0.92)) +
          ggplot2::labs(title = sprintf("E. Comparative Multi-Taxon Ensemble Habitat Suitability (%d Modeled Taxa)", n_sp_m), subtitle = sub_e_txt, x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)")
      }
    }
  }

  f_raw <- trimws(as.character(focal_taxon)); f_tx <- if (nzchar(f_raw)) f_raw else active_groups[1L]
  pts_f <- if (has_focal_sdm && f_tx %in% vdf$plot_group) vdf[vdf$plot_group == f_tx, , drop = FALSE]
    else if (has_focal_sdm && ("species_standardized" %in% names(vdf)) && (f_tx %in% vdf$species_standardized)) vdf[vdf$species_standardized == f_tx, , drop = FALSE] else vdf
  info_max <- if (has_focal_sdm) .sdm_raster_to_proj_df(focal_rasters$maxent) else list(df = data.frame(), projected = FALSE)
  info_ens <- if (has_focal_sdm) .sdm_raster_to_proj_df(focal_rasters$ensemble) else list(df = data.frame(), projected = FALSE)
  out_canon <- file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png")

  if (nrow(info_max$df) > 0L || nrow(info_ens$df) > 0L) {
    .mk_t <- function(pre) if (exists("format_icn_title_expr", mode = "function")) format_icn_title_expr(prefix = pre, taxon = f_tx) else paste0(pre, gsub("_", " ", f_tx))
    p_a <- build_suitability_panel(if (nrow(info_max$df) > 0L) info_max else info_ens, pts_f, .mk_t("A. MaxEnt (cloglog) Suitability: "), "L1-regularized IPP presence-background model (EPSG:5070)")
    p_b <- build_suitability_panel(if (nrow(info_ens$df) > 0L) info_ens else info_max, pts_f, .mk_t("B. Spatial RF + MaxEnt Ensemble: "), "AUC-weighted ensemble consensus surface (EPSG:5070)")

    p_c <- if (!is.null(resp_df) && is.data.frame(resp_df) && nrow(resp_df) > 0L && all(c("predictor", "x_val", "suitability") %in% names(resp_df))) {
      ggplot2::ggplot(resp_df, ggplot2::aes(x = .data$x_val, y = .data$suitability, color = .data$predictor)) +
        ggplot2::geom_line(linewidth = 0.95, show.legend = FALSE) +
        ggplot2::facet_wrap(~ predictor, scales = "free_x", ncol = min(3L, length(unique(resp_df$predictor)))) +
        ggplot2::scale_color_manual(values = rep(c("#0072B2", "#009E73", "#D55E00", "#CC79A7"), length.out = length(unique(resp_df$predictor)))) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) + theme_packera_pub(base_size = 9.2) +
        ggplot2::labs(title = "C. Marginal MaxEnt (cloglog) Ecological Response Curves", subtitle = "5th to 95th percentile gradient holding co-predictors at presence median", x = "Predictor Value", y = "Predicted Suitability")
    } else ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Marginal Response Curves")

    p_d <- if (!is.null(roc_df) && is.data.frame(roc_df) && nrow(roc_df) > 0L && all(c("fpr", "tpr") %in% names(roc_df))) {
      p_row <- if (!is.null(perf_df) && is.data.frame(perf_df) && nrow(perf_df) > 0L) {
        if ("Taxon_or_Morphotype" %in% names(perf_df) && any(perf_df$Taxon_or_Morphotype == f_tx)) perf_df[perf_df$Taxon_or_Morphotype == f_tx, , drop = FALSE][1L, , drop = FALSE] else perf_df[1L, , drop = FALSE]
      } else NULL
      .pm <- function(cs) { if (is.null(p_row)) return(NA_real_); h <- intersect(cs, names(p_row)); if (length(h) == 0L) NA_real_ else suppressWarnings(as.numeric(p_row[[h[1L]]][1L])) }
      auc_v <- .pm(c("AUC_Block_Ensemble", "auc_ensemble", "AUC_Block_MaxEnt", "auc_maxent"))
      boy_v <- .pm(c("CBI_Ensemble", "boyce_index", "CBI_MaxEnt"))
      boyce_r <- if (!is.null(boyce_df) && is.data.frame(boyce_df) && nrow(boyce_df) >= 3L && all(c("suitability", "pe_ratio") %in% names(boyce_df))) {
        suppressWarnings(stats::cor(as.numeric(boyce_df$suitability), as.numeric(boyce_df$pe_ratio), method = "spearman", use = "complete.obs"))
      } else NA_real_
      if (!("model" %in% names(roc_df))) roc_df$model <- "Ensemble"
      sub_d <- paste(c(sprintf("Ensemble Block AUC = %s", ifelse(is.finite(auc_v), sprintf("%.3f", auc_v), "N/A")),
                       sprintf("CBI = %s", ifelse(is.finite(boy_v), sprintf("%+.3f", boy_v), "N/A")),
                       if (is.finite(boyce_r)) sprintf("Boyce Spearman r = %+.3f", boyce_r) else NULL), collapse = " | ")
      p_roc <- ggplot2::ggplot(roc_df, ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$model, linetype = .data$model)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#9CA3AF") +
        ggplot2::geom_line(linewidth = 0.95) +
        ggplot2::scale_color_manual(values = c("MaxEnt" = "#0072B2", "Spatial RF" = "#009E73", "Ensemble" = "#D55E00"), name = "SDM Engine") +
        ggplot2::scale_linetype_manual(values = c("MaxEnt" = "dashed", "Spatial RF" = "dotdash", "Ensemble" = "solid"), name = "SDM Engine") +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + theme_packera_pub() +
        ggplot2::labs(title = "D. Spatial-Block Cross-Validation & Boyce Calibration", subtitle = sub_d, x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)")
      if (is.finite(boyce_r)) {
        p_roc <- p_roc + ggplot2::annotate("label", x = 0.97, y = 0.05, hjust = 1, vjust = 0, size = 2.7, fill = grDevices::adjustcolor("white", alpha.f = 0.88),
                                           color = "#1F2937", label = sprintf("Boyce cal. Spearman r = %+.3f\n(n = %d bins)", boyce_r, nrow(boyce_df)))
      }
      p_roc
    } else ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "D. Spatial Block-CV ROC & Boyce Calibration")

    if (!is.null(p_multi_grid)) {
      plate_05c <- ((p_a | p_b) / (p_c | p_d) / p_multi_grid) +
        patchwork::plot_layout(heights = c(1.22, 0.92, max(1.15, 0.72 * n_rows_multi))) +
        patchwork::plot_annotation(title = "Step 05: Ensemble Species Distribution Modeling, Spatial-Block Validation & Multi-Taxon Suitability",
                                   theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4)))
      return(save_packera_png_300(plate_05c, out_canon, width = 14.2, height = 9.4 + 2.7 * n_rows_multi, also_copy_to = out_sdm_path))
    }
    plate_05c <- ((p_a | p_b) / (p_c | p_d)) +
      patchwork::plot_layout(widths = c(1.1, 1.1), heights = c(1.25, 0.95)) +
      patchwork::plot_annotation(title = "Step 05: Ensemble Species Distribution Modeling & Spatial Block Validation",
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4)))
    return(save_packera_png_300(plate_05c, out_canon, width = 14.0, height = 9.8, also_copy_to = out_sdm_path))
  }

  p_geo <- build_suitability_panel(data.frame(), vdf, "Step 05: Equal-Area Projected Geographic Distribution (Albers EPSG:5070)", "Hollow markers stratified by taxon; scale bar & north arrow in projected CRS")
  save_packera_png_300(p_geo, out_canon, width = 11.8, height = 8.2, also_copy_to = out_sdm_path)
}
