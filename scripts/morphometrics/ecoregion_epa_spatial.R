if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/ecoregion_epa_spatial.R
# US EPA Level III Ecoregion spatial polygons, county containment, and centroid quarantine
# ==============================================================================

.epa_l3_r_cache <- new.env(parent = emptyenv())

.LEGACY_HEURISTIC_REGIONS <- c(
  "coastal_plain_sandhills",
  "piedmont_granite_flatrocks",
  "appalachian_highlands",
  "interior_prairie_midwest",
  "other_us",
  "other_region",
  "other_state",
  "unassigned_na"
)

standardize_epa_l3_label_r <- function(raw_vec) {
  vapply(as.character(raw_vec), function(val) {
    if (is.na(val)) return(NA_character_)
    txt <- trimws(val)
    if (!nzchar(txt)) return(NA_character_)
    low <- tolower(txt)
    if (low %in% c("na", "nan", "none", "null", "unknown", .LEGACY_HEURISTIC_REGIONS)) {
      return(NA_character_)
    }
    cleaned <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", txt))
    if (!nzchar(cleaned) || tolower(cleaned) %in% .LEGACY_HEURISTIC_REGIONS) return(NA_character_)
    cleaned
  }, character(1), USE.NAMES = FALSE)
}

.norm_us_state_r <- function(st_vec) {
  abbr_map <- c(
    stats::setNames(tolower(datasets::state.name), tolower(datasets::state.abb)),
    "dc" = "district of columbia", "d.c." = "district of columbia",
    "on" = "ontario", "ont" = "ontario",
    "qc" = "quebec", "pq" = "quebec", "que" = "quebec",
    "nb" = "new brunswick", "ns" = "nova scotia",
    "pe" = "prince edward island", "pei" = "prince edward island",
    "nl" = "newfoundland and labrador", "nf" = "newfoundland and labrador", "nfld" = "newfoundland and labrador",
    "mb" = "manitoba", "man" = "manitoba",
    "sk" = "saskatchewan", "sask" = "saskatchewan",
    "ab" = "alberta", "alta" = "alberta",
    "bc" = "british columbia", "yt" = "yukon", "nt" = "northwest territories", "nu" = "nunavut"
  )
  vapply(as.character(st_vec), function(s) {
    if (is.na(s)) return(NA_character_)
    s_ascii <- iconv(s, from = "", to = "ASCII//TRANSLIT")
    if (is.na(s_ascii)) s_ascii <- s
    s_clean <- gsub("\\.$", "", gsub("[[:space:]]+", " ", tolower(trimws(s_ascii))))
    if (!nzchar(s_clean) || s_clean %in% c("na", "nan", "none", "unknown")) return(NA_character_)
    if (s_clean %in% names(abbr_map)) unname(abbr_map[[s_clean]]) else s_clean
  }, character(1), USE.NAMES = FALSE)
}

.HOMONYMOUS_INDEP_CITIES_R <- list(
  "virginia" = c("fairfax", "franklin", "richmond", "roanoke"),
  "maryland" = c("baltimore"),
  "missouri" = c("st louis")
)

.norm_us_county_r <- function(co_vec, st_norm_vec = NULL) {
  n <- length(co_vec)
  if (is.null(st_norm_vec)) st_norm_vec <- rep(NA_character_, n)
  vapply(seq_len(n), function(i) {
    raw <- as.character(co_vec[i])
    if (is.na(raw)) return(NA_character_)
    c_ascii <- iconv(raw, from = "", to = "ASCII//TRANSLIT")
    if (is.na(c_ascii)) c_ascii <- raw
    c_str <- tolower(trimws(c_ascii))
    if (!nzchar(c_str) || c_str %in% c("na", "nan", "none", "unknown", "unplaced", "indeterminate")) {
      return(NA_character_)
    }
    if (grepl(" or |/|\\?", c_str)) return(NA_character_)
    c_str <- gsub("\\[|\\]|\\(|\\)", "", c_str)
    was_explicit_city <- grepl("^city[[:space:]]+of[[:space:]]+|(^|[^A-Za-z0-9_])(independent[[:space:]]+city|city)$", c_str) && !grepl("(^|[^A-Za-z0-9_])county([^A-Za-z0-9_]|$)", c_str)
    c_str <- gsub("^(regional[[:space:]]+municipality[[:space:]]+of|united[[:space:]]+counties[[:space:]]+of|county[[:space:]]+of|city[[:space:]]+of|mrc[[:space:]]+de|mrc[[:space:]]+d'|municipalite[[:space:]]+regionale[[:space:]]+de[[:space:]]+comte[[:space:]]+de)[[:space:]]+", "", c_str)
    c_str <- gsub("[[:space:]]+(regional[[:space:]]+municipality|united[[:space:]]+counties|county[[:space:]]+municipality|census[[:space:]]+division|census[[:space:]]+area|city[[:space:]]+and[[:space:]]+borough|independent[[:space:]]+city|independent|municipality|county|parish|par\\.?|borough|district|mrc|co\\.?)$", "", c_str)
    c_str <- gsub("[[:space:]]+city$", "", c_str)
    c_str <- gsub("^(st\\.|saint)[[:space:]]+", "st ", c_str)
    c_str <- gsub("^ste\\.[[:space:]]*", "ste ", c_str)
    c_str <- trimws(gsub("[^a-z0-9]+", " ", c_str))
    if (!nzchar(c_str)) return(NA_character_)
    st_i <- as.character(st_norm_vec[i])
    if (!is.na(st_i)) {
      if (st_i == "virginia" && c_str == "nansemond") c_str <- "suffolk"
      else if (st_i == "virginia" && c_str == "princess anne") c_str <- "virginia beach"
      else if (st_i == "virginia" && c_str == "bedford") was_explicit_city <- FALSE
      else if (st_i == "georgia" && c_str == "columbus") c_str <- "muscogee"
      else if (st_i == "south dakota" && c_str == "shannon") c_str <- "oglala lakota"
      if (isTRUE(was_explicit_city) && (st_i %in% names(.HOMONYMOUS_INDEP_CITIES_R)) && (c_str %in% .HOMONYMOUS_INDEP_CITIES_R[[st_i]])) {
        c_str <- paste0(c_str, " city")
      }
    }
    c_str
  }, character(1), USE.NAMES = FALSE)
}

.build_offline_epa_l3_sf_r <- function() {
  mk_poly <- function(mat) sf::st_polygon(list(rbind(mat, mat[1, , drop = FALSE])))
  eco_df <- data.frame(
    ecoregion_label = c(
      "Blue_Ridge", "Ridge_and_Valley", "Piedmont", "Middle_Atlantic_Coastal_Plain",
      "Southeastern_Plains", "Southern_Coastal_Plain", "Mississippi_Alluvial_Plain",
      "South_Central_Plains", "East_Central_Texas_Plains", "Ouachita_Mountains",
      "Ozark_Highlands", "Interior_Plateau", "Southwestern_Appalachians",
      "Mixed_Wood_Plains", "Northern_Lakes_and_Forests", "Atlantic_Highlands",
      "Lake_Manitoba_and_Lake_Winnipeg_Plain", "Aspen_Parkland_Northern_Glaciated_Plains"
    ),
    stringsAsFactors = FALSE
  )
  eco_geoms <- sf::st_sfc(
    mk_poly(matrix(c(-84.5,34.5, -83.2,34.8, -81.6,36.1, -79.0,37.9, -78.2,38.9, -79.5,39.2, -82.8,36.8, -84.6,35.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-85.6,34.2, -84.6,35.2, -82.8,36.8, -79.5,39.2, -77.2,40.8, -78.8,41.2, -83.8,37.4, -85.9,34.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-85.2,32.5, -82.0,33.5, -79.4,34.9, -77.6,36.5, -77.1,38.9, -78.2,38.9, -79.0,37.9, -81.6,36.1, -83.2,34.8, -84.5,34.5), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-79.6,33.0, -75.4,35.2, -75.0,38.8, -76.4,39.2, -77.1,38.5, -77.6,36.5, -78.5,35.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-90.2,31.0, -85.0,30.8, -81.3,32.2, -79.6,33.0, -78.5,35.2, -79.4,34.9, -82.0,33.5, -85.2,32.5, -88.8,34.8, -90.1,33.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-85.0,29.6, -80.2,26.5, -80.8,31.2, -81.3,32.2, -85.0,30.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-91.8,29.5, -90.2,30.5, -90.1,35.0, -89.2,37.0, -90.6,37.0, -91.5,34.2, -92.0,31.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-95.8,30.2, -92.0,31.0, -91.5,34.2, -94.8,34.1, -96.0,32.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-98.2,29.0, -95.8,30.2, -96.0,32.8, -97.2,33.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-95.6,34.1, -92.4,34.3, -92.6,35.1, -95.4,35.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-94.8,35.7, -91.0,35.8, -90.5,38.2, -94.6,38.1), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-88.5,34.8, -85.5,35.2, -84.5,38.2, -87.8,38.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-86.8,33.6, -85.6,34.2, -84.2,36.8, -85.5,36.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-83.5,41.8, -71.0,45.0, -71.0,47.2, -83.5,44.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-95.2,46.0, -76.0,45.2, -76.0,50.5, -95.2,50.5), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-71.0,45.0, -64.0,45.0, -64.0,49.2, -71.0,49.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-99.5,49.0, -96.0,49.0, -96.0,53.0, -99.5,53.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-104.05,49.0, -99.5,49.0, -99.5,53.5, -104.05,53.5), ncol = 2, byrow = TRUE)),
    crs = 4326
  )
  eco_sf <- sf::st_sf(eco_df, geometry = eco_geoms)
  cnty_df <- data.frame(
    county_key = c("north carolina|orange", "north carolina|new hanover", "north carolina|wake", "texas|nacogdoches", "ontario|essex", "quebec|gatineau"),
    state_name = c("north carolina", "north carolina", "north carolina", "texas", "ontario", "quebec"),
    county_norm = c("orange", "new hanover", "wake", "nacogdoches", "essex", "gatineau"),
    centroid_lon = c(-79.10, -77.90, -78.40, -94.65, -82.80, -75.95),
    centroid_lat = c(36.05, 34.20, 35.78, 31.62, 42.15, 46.10),
    n_ecoregions = rep(NA_integer_, 6L),
    top_ecoregion = rep(NA_character_, 6L),
    top_area_share = rep(NA_real_, 6L),
    stringsAsFactors = FALSE
  )
  cnty_geoms <- sf::st_sfc(
    mk_poly(matrix(c(-79.25,35.86, -78.95,35.86, -78.95,36.24, -79.25,36.24), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-78.02,34.05, -77.78,34.05, -77.78,34.35, -78.02,34.35), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-79.00,35.50, -77.80,35.50, -77.80,36.05, -79.00,36.05), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-94.90,31.40, -94.40,31.40, -94.40,31.85, -94.90,31.85), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-83.10,41.95, -82.50,41.95, -82.50,42.32, -83.10,42.32), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-76.30,45.75, -75.65,45.75, -75.65,46.45, -76.30,46.45), ncol = 2, byrow = TRUE)),
    crs = 4326
  )
  cnty_sf <- sf::st_sf(cnty_df, geometry = cnty_geoms)
  list(eco_sf = eco_sf, cnty_sf = cnty_sf, use_live_county_intersection = TRUE)
}

load_epa_level3_spatial_r <- function(ecoregion_path = NULL, county_path = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required for EPA Level III Ecoregion lookup.")
  custom_paths <- !is.null(ecoregion_path) || !is.null(county_path)
  default_dir <- "data/environmental/epa_ecoregions_l3"
  gpkg_p <- file.path(default_dir, "us_eco_l3.gpkg")
  na_eco_p <- file.path(default_dir, "na_cec_eco_l3.geojson")
  ca_cnty_p <- file.path(default_dir, "ca_census_divisions.geojson")
  eco_json_p <- if (!is.null(ecoregion_path)) as.character(ecoregion_path) else if (file.exists(na_eco_p)) na_eco_p else file.path(default_dir, "us_eco_l3.geojson")
  cnty_json_p <- if (!is.null(county_path)) as.character(county_path) else file.path(default_dir, "us_counties.geojson")

  cache_id <- paste(
    if (is.null(ecoregion_path)) "DEFAULT_ECO" else eco_json_p,
    if (is.null(county_path)) "DEFAULT_CNTY" else cnty_json_p,
    sep = "::"
  )
  if (exists(cache_id, envir = .epa_l3_r_cache, inherits = FALSE)) {
    return(get(cache_id, envir = .epa_l3_r_cache, inherits = FALSE))
  }

  .merge_ca_cec_r <- function(eco_sf_in, cnty_sf_in) {
    if (!custom_paths && file.exists(na_eco_p) && !identical(eco_json_p, na_eco_p)) {
      na_sf <- tryCatch(sf::st_read(na_eco_p, quiet = TRUE), error = function(e) NULL)
      if (!is.null(na_sf) && nrow(na_sf) > 0L) {
        na_lbl <- if ("ecoregion_label" %in% names(na_sf)) na_sf$ecoregion_label else if ("NA_L3NAME" %in% names(na_sf)) na_sf$NA_L3NAME else if ("US_L3NAME" %in% names(na_sf)) na_sf$US_L3NAME else na_sf[[1]]
        na_sf$ecoregion_label <- standardize_epa_l3_label_r(na_lbl)
        na_sf <- sf::st_transform(na_sf[, "ecoregion_label", drop = FALSE], sf::st_crs(eco_sf_in))
        gcol_e <- attr(eco_sf_in, "sf_column")
        if (!is.null(gcol_e) && !identical(attr(na_sf, "sf_column"), gcol_e)) sf::st_geometry(na_sf) <- gcol_e
        eco_sf_in <- rbind(eco_sf_in[, "ecoregion_label", drop = FALSE], na_sf[, "ecoregion_label", drop = FALSE])
      }
    }
    if (!is.null(eco_sf_in) && !is.na(sf::st_crs(eco_sf_in)) && !isTRUE(sf::st_is_longlat(eco_sf_in))) {
      eco_sf_in <- sf::st_transform(eco_sf_in, 4326)
    }
    if (!is.null(cnty_sf_in) && !is.na(sf::st_crs(cnty_sf_in)) && !isTRUE(sf::st_is_longlat(cnty_sf_in))) {
      cnty_sf_in <- sf::st_transform(cnty_sf_in, 4326)
    }
    if (!custom_paths && file.exists(ca_cnty_p) && !is.null(cnty_sf_in)) {
      ca_sf <- tryCatch(sf::st_read(ca_cnty_p, quiet = TRUE), error = function(e) NULL)
      if (!is.null(ca_sf) && nrow(ca_sf) > 0L) {
        ca_sf <- sf::st_transform(ca_sf, sf::st_crs(cnty_sf_in))
        gcol_c <- attr(cnty_sf_in, "sf_column")
        if (!is.null(gcol_c) && !identical(attr(ca_sf, "sf_column"), gcol_c)) sf::st_geometry(ca_sf) <- gcol_c
        ca_st_c <- intersect(c("state_name", "PRNAME", "STATE_NAME", "state"), names(ca_sf))
        ca_co_c <- intersect(c("county_norm", "CDNAME", "NAME", "county"), names(ca_sf))
        if (length(ca_st_c) > 0L && length(ca_co_c) > 0L) {
          ca_st_v <- .norm_us_state_r(ca_sf[[ca_st_c[1L]]])
          ca_co_v <- .norm_us_county_r(ca_sf[[ca_co_c[1L]]], ca_st_v)
          ca_sf$state_name <- ca_st_v
          ca_sf$county_norm <- ca_co_v
          ca_sf$county_key <- paste0(ca_st_v, "|", ca_co_v)
        }
        non_geom_c <- setdiff(names(cnty_sf_in), gcol_c)
        for (miss_c in setdiff(non_geom_c, names(ca_sf))) ca_sf[[miss_c]] <- NA
        cnty_sf_in <- rbind(cnty_sf_in[, c(non_geom_c, gcol_c), drop = FALSE], ca_sf[, c(non_geom_c, gcol_c), drop = FALSE])
      }
    }
    list(eco_sf = eco_sf_in, cnty_sf = cnty_sf_in)
  }

  if (!custom_paths && file.exists(gpkg_p)) {
    eco_sf <- sf::st_read(gpkg_p, layer = "us_eco_l3", quiet = TRUE)
    cnty_sf <- tryCatch(sf::st_read(gpkg_p, layer = "us_counties", quiet = TRUE), error = function(e) NULL)
    eco_sf$ecoregion_label <- standardize_epa_l3_label_r(
      if ("ecoregion_label" %in% names(eco_sf)) eco_sf$ecoregion_label else if ("NA_L3NAME" %in% names(eco_sf)) eco_sf$NA_L3NAME else eco_sf$US_L3NAME
    )
    m_ca <- .merge_ca_cec_r(eco_sf, cnty_sf)
    bundle <- list(eco_sf = m_ca$eco_sf, cnty_sf = m_ca$cnty_sf, use_live_county_intersection = FALSE)
    assign(cache_id, bundle, envir = .epa_l3_r_cache)
    return(bundle)
  }

  if (file.exists(eco_json_p)) {
    eco_sf <- sf::st_read(eco_json_p, quiet = TRUE)
    raw_lbl <- if ("ecoregion_label" %in% names(eco_sf)) eco_sf$ecoregion_label else if ("NA_L3NAME" %in% names(eco_sf)) eco_sf$NA_L3NAME else if ("US_L3NAME" %in% names(eco_sf)) eco_sf$US_L3NAME else eco_sf[[1]]
    eco_sf$ecoregion_label <- standardize_epa_l3_label_r(raw_lbl)
    cnty_sf <- if (file.exists(cnty_json_p)) tryCatch(sf::st_read(cnty_json_p, quiet = TRUE), error = function(e) NULL) else NULL
    m_ca <- .merge_ca_cec_r(eco_sf, cnty_sf)
    bundle <- list(eco_sf = m_ca$eco_sf, cnty_sf = m_ca$cnty_sf, use_live_county_intersection = custom_paths || identical(eco_json_p, na_eco_p))
    assign(cache_id, bundle, envir = .epa_l3_r_cache)
    return(bundle)
  }

  .build_offline_epa_l3_sf_r()
}

infer_epa_level3_ecoregion <- function(
  lat,
  lon,
  state_province = NULL,
  county = NULL,
  ecoregion_path = NULL,
  county_path = NULL
) {
  n <- max(length(lat), length(lon), 1L)
  if (length(lat) == 0L && length(lon) == 0L) return(character(0))
  lat_num <- suppressWarnings(as.numeric(rep_len(lat, n)))
  lon_num <- suppressWarnings(as.numeric(rep_len(lon, n)))
  st_vec <- if (!is.null(state_province) && length(state_province) > 0L) rep_len(as.character(state_province), n) else rep(NA_character_, n)
  co_vec <- if (!is.null(county) && length(county) > 0L) rep_len(as.character(county), n) else rep(NA_character_, n)

  old_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  bundle <- load_epa_level3_spatial_r(ecoregion_path = ecoregion_path, county_path = county_path)
  eco_sf <- bundle$eco_sf
  cnty_sf <- bundle$cnty_sf
  out <- rep(NA_character_, n)

  # (a) Finite coordinates: strict point-in-polygon intersection
  has_coords <- is.finite(lat_num) & is.finite(lon_num) & (lat_num >= -90) & (lat_num <= 90) & (lon_num >= -180) & (lon_num <= 180)
  if (any(has_coords)) {
    pts_sf <- sf::st_as_sf(
      data.frame(lon = lon_num[has_coords], lat = lat_num[has_coords]),
      coords = c("lon", "lat"),
      crs = 4326
    )
    hits <- suppressWarnings(sf::st_intersects(pts_sf, eco_sf))
    out[has_coords] <- vapply(hits, function(h) {
      if (length(h) >= 1L) as.character(eco_sf$ecoregion_label[h[1]]) else NA_character_
    }, character(1))
  }

  # (b) Missing coordinates: county-containment fallback (>= 95% single-ecoregion area share)
  miss_idx <- which(!has_coords)
  if (length(miss_idx) > 0L && !is.null(cnty_sf) && nrow(cnty_sf) > 0L) {
    st_norm <- .norm_us_state_r(st_vec[miss_idx])
    co_norm <- .norm_us_county_r(co_vec[miss_idx], st_norm)
    c_st_col <- intersect(c("state_name", "STATE_NAME", "state", "PRNAME"), names(cnty_sf))
    c_st <- .norm_us_state_r(if (length(c_st_col) > 0L) cnty_sf[[c_st_col[1L]]] else if ("county_key" %in% names(cnty_sf)) sub("\\|.*$", "", as.character(cnty_sf$county_key)) else cnty_sf[[1]])
    c_co_col <- intersect(c("county_norm", "NAME", "CDNAME", "county"), names(cnty_sf))
    c_co <- .norm_us_county_r(if (length(c_co_col) > 0L) cnty_sf[[c_co_col[1L]]] else if ("county_key" %in% names(cnty_sf)) sub("^.*\\|", "", as.character(cnty_sf$county_key)) else cnty_sf[[min(2L, ncol(cnty_sf))]], c_st)
    c_keys <- paste0(c_st, "|", c_co)
    c_keys_ns <- paste0(c_st, "|", gsub(" ", "", c_co))

    for (k_i in seq_along(miss_idx)) {
      row_i <- miss_idx[k_i]
      if (is.na(st_norm[k_i]) || is.na(co_norm[k_i])) next
      q_key <- paste0(st_norm[k_i], "|", co_norm[k_i])
      q_key_ns <- paste0(st_norm[k_i], "|", gsub(" ", "", co_norm[k_i]))
      m_idx <- match(q_key, c_keys)
      if (is.na(m_idx)) m_idx <- match(q_key_ns, c_keys_ns)
      if (is.na(m_idx)) next

      c_row <- cnty_sf[m_idx, ]
      use_live <- isTRUE(bundle$use_live_county_intersection) || !("top_area_share" %in% names(c_row)) || is.na(c_row$top_area_share[1])
      if (use_live) {
        ix <- tryCatch(suppressWarnings(sf::st_intersection(c_row, eco_sf[, "ecoregion_label", drop = FALSE])), error = function(e) NULL)
        if (is.null(ix) || nrow(ix) == 0L) next
        a_parts <- as.numeric(sf::st_area(ix))
        a_by_eco <- tapply(a_parts, ix$ecoregion_label, sum, na.rm = TRUE)
        tot_a <- sum(a_by_eco, na.rm = TRUE)
        c_tot_a <- as.numeric(sf::st_area(c_row))[1]
        if (!is.finite(tot_a) || tot_a <= 0 || (is.finite(c_tot_a) && c_tot_a > 0 && (tot_a / c_tot_a) < 0.50)) next
        top_idx <- which.max(a_by_eco)
        top_eco <- names(a_by_eco)[top_idx]
        top_share <- as.numeric(a_by_eco[top_idx]) / tot_a
        n_eco <- sum((a_by_eco / tot_a) >= 0.001)
      } else {
        top_eco <- standardize_epa_l3_label_r(c_row$top_ecoregion[1])
        top_share <- suppressWarnings(as.numeric(c_row$top_area_share[1]))
        n_eco <- suppressWarnings(as.integer(c_row$n_ecoregions[1]))
      }

      if (is.na(top_eco) || !(isTRUE(n_eco == 1L) || (is.finite(top_share) && top_share >= 0.95))) next

      c_lon <- if ("centroid_lon" %in% names(c_row)) suppressWarnings(as.numeric(c_row$centroid_lon[1])) else NA_real_
      c_lat <- if ("centroid_lat" %in% names(c_row)) suppressWarnings(as.numeric(c_row$centroid_lat[1])) else NA_real_
      if (!is.finite(c_lon) || !is.finite(c_lat)) {
        pos_pt <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(c_row)))
        c_lon <- pos_pt[1, 1]; c_lat <- pos_pt[1, 2]
      }
      c_pt_sf <- sf::st_as_sf(data.frame(lon = c_lon, lat = c_lat), coords = c("lon", "lat"), crs = 4326)
      c_hit <- suppressWarnings(sf::st_intersects(c_pt_sf, eco_sf))[[1]]
      hit_lbl <- if (length(c_hit) >= 1L) as.character(eco_sf$ecoregion_label[c_hit[1]]) else NA_character_
      if (!is.na(hit_lbl) && hit_lbl == top_eco) {
        out[row_i] <- hit_lbl
      } else if (is.na(hit_lbl)) {
        out[row_i] <- top_eco
      }
    }
  }
  out
}

.is_quarantined_ecoregion_coord <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(logical(0))
  if (!("is_centroid_or_inst_coord" %in% names(df))) return(rep(FALSE, nrow(df)))
  v <- df$is_centroid_or_inst_coord
  if (is.logical(v)) return(!is.na(v) & v)
  if (is.numeric(v)) return(!is.na(v) & is.finite(v) & (v != 0))
  lv <- tolower(trimws(as.character(v)))
  !is.na(lv) & (lv %in% c("true", "1", "t", "yes"))
}
