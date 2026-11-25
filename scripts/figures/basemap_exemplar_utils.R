if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/figures/basemap_exemplar_utils.R
# Equal-area sf North America basemaps, LM3 exemplar asset loaders, and trait labels
# ==============================================================================

build_na_projected_basemap_gg <- function(lon_vec, lat_vec, target_crs = 5070L) {
  lon <- suppressWarnings(as.numeric(lon_vec))
  lat <- suppressWarnings(as.numeric(lat_vec))
  ok <- is.finite(lon) & is.finite(lat) & lat >= 10 & lat <= 75 & lon >= -170 & lon <= -50
  lon_use <- if (any(ok)) lon[ok] else c(-96.0, -75.0)
  lat_use <- if (any(ok)) lat[ok] else c(30.0, 44.0)

  pad_lon <- max(2.2, 0.14 * diff(range(lon_use)))
  pad_lat <- max(1.8, 0.14 * diff(range(lat_use)))
  xlim_ll <- c(max(-168, min(lon_use) - pad_lon), min(-52, max(lon_use) + pad_lon))
  ylim_ll <- c(max(18, min(lat_use) - pad_lat), min(72, max(lat_use) + pad_lat))

  can_sf <- requireNamespace("sf", quietly = TRUE) && requireNamespace("maps", quietly = TRUE)
  if (isTRUE(can_sf)) {
    bg_sf <- tryCatch({
      m_usa <- maps::map("state", plot = FALSE, fill = TRUE)
      sf_usa <- sf::st_as_sf(m_usa)
      if (is.na(sf::st_crs(sf_usa))) sf_usa <- sf::st_set_crs(sf_usa, 4326) else sf_usa <- sf::st_transform(sf_usa, 4326)
      sf_usa <- sf::st_make_valid(sf_usa)
      sf::st_transform(sf_usa, target_crs)
    }, error = function(e) NULL)

    bbox_proj <- tryCatch({
      pts_bb <- sf::st_sfc(
        sf::st_point(c(xlim_ll[1], ylim_ll[1])),
        sf::st_point(c(xlim_ll[2], ylim_ll[1])),
        sf::st_point(c(xlim_ll[2], ylim_ll[2])),
        sf::st_point(c(xlim_ll[1], ylim_ll[2])),
        crs = 4326
      )
      bb <- sf::st_bbox(sf::st_transform(pts_bb, target_crs))
      bb
    }, error = function(e) NULL)

    if (!is.null(bg_sf) && !is.null(bbox_proj)) {
      # Compute scale-bar length in meters rounded to nice km
      span_x_m <- as.numeric(bbox_proj["xmax"] - bbox_proj["xmin"])
      span_y_m <- as.numeric(bbox_proj["ymax"] - bbox_proj["ymin"])
      nice_km_opts <- c(25, 50, 100, 200, 250, 500, 1000)
      bar_km <- nice_km_opts[which.min(abs(nice_km_opts * 1000 - 0.20 * span_x_m))]
      bar_m <- bar_km * 1000.0
      bx1 <- as.numeric(bbox_proj["xmax"]) - 0.05 * span_x_m
      bx0 <- bx1 - bar_m
      by0 <- as.numeric(bbox_proj["ymin"]) + 0.06 * span_y_m
      nx <- as.numeric(bbox_proj["xmax"]) - 0.06 * span_x_m
      ny0 <- as.numeric(bbox_proj["ymax"]) - 0.14 * span_y_m
      ny1 <- as.numeric(bbox_proj["ymax"]) - 0.05 * span_y_m

      p_base <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = bg_sf, fill = "#FAFAFA", color = "#CBD5E1", linewidth = 0.30) +
        ggplot2::annotate("segment", x = bx0, xend = bx1, y = by0, yend = by0, linewidth = 0.85, color = "#1F2937") +
        ggplot2::annotate("text", x = 0.5 * (bx0 + bx1), y = by0 + 0.032 * span_y_m, label = paste0(bar_km, " km"), size = 2.7, color = "#1F2937") +
        ggplot2::annotate("segment", x = nx, xend = nx, y = ny0, yend = ny1, linewidth = 0.75, color = "#1F2937",
                          arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "open")) +
        ggplot2::annotate("text", x = nx, y = ny1 + 0.028 * span_y_m, label = "N", fontface = "bold", size = 2.9, color = "#1F2937") +
        ggplot2::coord_sf(
          crs = target_crs,
          xlim = c(as.numeric(bbox_proj["xmin"]), as.numeric(bbox_proj["xmax"])),
          ylim = c(as.numeric(bbox_proj["ymin"]), as.numeric(bbox_proj["ymax"])),
          expand = FALSE
        ) +
        theme_packera_pub() +
        ggplot2::labs(x = "Longitude", y = "Latitude")
      return(list(plot = p_base, projected = TRUE, crs = target_crs, xlim = xlim_ll, ylim = ylim_ll))
    }
  }

  # Conformal fallback when sf/maps projection is unavailable
  p_fb <- ggplot2::ggplot()
  if (requireNamespace("maps", quietly = TRUE)) {
    st_df <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
    if (!is.null(st_df) && nrow(st_df) > 0L) {
      p_fb <- p_fb + ggplot2::geom_polygon(
        data = st_df,
        ggplot2::aes(x = .data$long, y = .data$lat, group = .data$group),
        fill = "#FAFAFA", color = "#CBD5E1", linewidth = 0.30
      )
    }
  }
  cos_lat <- cos(mean(ylim_ll) * pi / 180.0)
  p_fb <- p_fb +
    ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_lat), xlim = xlim_ll, ylim = ylim_ll, expand = FALSE) +
    theme_packera_pub() +
    ggplot2::labs(x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)")
  list(plot = p_fb, projected = FALSE, crs = 4326L, xlim = xlim_ll, ylim = ylim_ll)
}

project_lonlat_points_r <- function(df, map_ctx, lon_col = "decimalLongitude", lat_col = "decimalLatitude") {
  if (is.null(df) || nrow(df) == 0L || !all(c(lon_col, lat_col) %in% names(df))) return(data.frame())
  lon <- suppressWarnings(as.numeric(df[[lon_col]]))
  lat <- suppressWarnings(as.numeric(df[[lat_col]]))
  ok <- is.finite(lon) & is.finite(lat) & abs(lat) <= 90 & abs(lon) <= 180
  if (!any(ok)) return(data.frame())
  sub <- df[ok, , drop = FALSE]
  sub$map_x <- lon[ok]
  sub$map_y <- lat[ok]
  if (isTRUE(map_ctx$projected) && requireNamespace("sf", quietly = TRUE)) {
    xy_p <- tryCatch({
      pts_sf <- sf::st_as_sf(sub, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
      pts_tr <- sf::st_transform(pts_sf, map_ctx$crs)
      sf::st_coordinates(pts_tr)
    }, error = function(e) NULL)
    if (!is.null(xy_p) && nrow(xy_p) == nrow(sub)) {
      sub$map_x <- as.numeric(xy_p[, 1L])
      sub$map_y <- as.numeric(xy_p[, 2L])
    }
  }
  sub
}

resolve_lm3_exemplar_manifest_r <- function(figures_dir, zone = "basal") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(list())
  z <- tolower(trimws(as.character(zone)))
  cands <- unique(c(
    file.path(figures_dir, "lm3_assets", sprintf("manifest_%s.json", z)),
    file.path(figures_dir, "ect_overlays", sprintf("manifest_%s.json", z)),
    file.path("outputs/figures/lm3_assets", sprintf("manifest_%s.json", z)),
    file.path("outputs/figures/ect_overlays", sprintf("manifest_%s.json", z))
  ))
  for (cp in cands) {
    if (file.exists(cp) && file.info(cp)$size > 2L) {
      res <- tryCatch(jsonlite::fromJSON(cp, simplifyVector = FALSE), error = function(e) NULL)
      if (is.list(res) && length(res) > 0L) return(res)
    }
  }
  list()
}

read_png_raster_grob_safe <- function(png_path) {
  if (is.null(png_path) || !nzchar(as.character(png_path)) || !file.exists(png_path)) return(NULL)
  if (!requireNamespace("png", quietly = TRUE)) return(NULL)
  arr <- tryCatch(png::readPNG(png_path), error = function(e) NULL)
  if (is.null(arr)) return(NULL)
  grid::rasterGrob(arr, interpolate = TRUE)
}

PACKERA_MORPH_TRAIT_LABEL_DICT <- c(
  "ect_lobe_complexity"            = "Topological ECT Lobe Complexity",
  "solidity"                       = "Lamina Solidity (Area / Convex Hull)",
  "q25_solidity"                   = "Lower-Quartile Lamina Solidity (Q25)",
  "dissection_index"               = "Perimeter Dissection Ratio",
  "blade_cie_l_norm"               = "Normalized Lamina Brightness (CIE L*)",
  "blade_cie_l_sd"                 = "Lamina Brightness Dispersion (SD L*)",
  "petiole_to_blade_ratio"         = "Petiole-to-Blade Length Ratio",
  "scape_to_basal_leaf_ratio"      = "Scape-to-Basal Leaf Ratio",
  "AR_inv"                         = "Involucre Aspect Ratio (H / W)",
  "aspect_ratio"                   = "Leaf Aspect Ratio (L / W)",
  "basal_cordate_index"            = "Basal Cordate Sinus Index",
  "terminal_lobe_area_frac"        = "Terminal Lobe Area Fraction",
  "margin_serration_index"         = "Margin Serration Index",
  "margin_tooth_acuteness"         = "Margin Tooth Acuteness",
  "base_angle_deg"                 = "Lamina Base Angle (deg)",
  "apex_angle_deg"                 = "Lamina Apex Angle (deg)",
  "lamina_curvature"               = "Lamina Midrib Curvature",
  "lamina_trace_length_mm"         = "Lamina Midrib Length (mm)",
  "radiate_head_fraction"          = "Radiate Capitulum Fraction",
  "blade_cie_l_q85_sheet_delta"    = "Lamina Brightness Delta (Q85 CIE L*)",
  "blade_cie_l_sheet_delta"        = "Lamina Brightness Delta (CIE L*)",
  "basal_zone_cie_l_delta"         = "Basal-Zone Brightness Delta (CIE L*)",
  "abaxial_cie_l_sheet_delta"      = "Abaxial Brightness Delta (CIE L*)",
  "blade_cie_a_browning_delta"     = "Lamina Browning Delta (CIE a*)",
  "blade_cie_b_yellowing_delta"    = "Lamina Yellowing Delta (CIE b*)",
  "log1p_capitulum_count"          = "Capitulum Count (log1p)",
  "basal_vs_cauline_length_ratio"  = "Basal-to-Cauline Blade Length Ratio",
  "basal_vs_cauline_solidity_delta"= "Basal-to-Cauline Solidity Shift",
  "heteroblasty_efa_distance"      = "Basal-to-Cauline EFA Divergence",
  "heteroblasty_ect_distance"      = "Within-Voucher Basal-to-Cauline ECT Divergence",
  "within_rosette_efa_dispersion"  = "Within-Rosette EFA Shape Dispersion",
  "cauline_heteroblasty_pc1"       = "Cauline Heteroblasty PC1",
  "log_leaf_area_mm2"              = "Leaf Lamina Area (log mm^2)",
  "foliar_variance"                = "Within-Individual Foliar Variance"
)

format_morph_trait_labels <- function(var_names) {
  vapply(as.character(var_names), function(nm) {
    if (is.na(nm) || !nzchar(trimws(nm))) return("")
    base_nm <- sub("_(std|scaled)$", "", trimws(nm))
    if (base_nm %in% names(PACKERA_MORPH_TRAIT_LABEL_DICT)) {
      return(PACKERA_MORPH_TRAIT_LABEL_DICT[[base_nm]])
    }
    if (grepl("^ECT_PC[0-9]+$", base_nm)) {
      return(sub("^ECT_PC([0-9]+)$", "Topological ECT PC\\1", base_nm))
    }
    if (grepl("^PC[0-9]+$", base_nm)) {
      return(sub("^PC([0-9]+)$", "EFA Morphometric PC\\1", base_nm))
    }
    tools::toTitleCase(gsub("_", " ", base_nm))
  }, character(1), USE.NAMES = FALSE)
}
