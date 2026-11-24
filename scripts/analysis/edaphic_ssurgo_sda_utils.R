#!/usr/bin/env Rscript
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Script: scripts/analysis/edaphic_ssurgo_sda_utils.R
# Project: Plantmatics / Packera Species Delimitation Pipeline (NCU)
# Description:
#   USDA SSURGO Vector Pedology (SDA / soilDB) focal 30-50 m querying,
#   component-weighted rhizosphere aggregation, rock-outcrop vs. lithic-series
#   separation, and outcrop functional classification (pure POSIX ERE).
# ==============================================================================

SSURGO_CANONICAL_COLS <- c(
  "catalogNumber", "latitude", "longitude", "musym", "muname", "taxsuborder", "taxgreatgroup",
  "is_rock_outcrop", "outcrop_class", "ssurgo_rock_outcrop_pct", "ssurgo_lithic_series_pct",
  "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph",
  "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm"
)

.ssurgo_outcrop_series_pat <- function() {
  if (exists("SSURGO_OUTCROP_SERIES_REGEX", inherits = TRUE)) get("SSURGO_OUTCROP_SERIES_REGEX", inherits = TRUE) else "(^|[^[:alnum:]_])(rock outcrop|rock land|flatrock|glade|boulders|wake|ashlar|gladeville|barfield|goldston|manteo|ramsey|lily|chrome|conowingo)($|[^[:alnum:]_])"
}

parse_ssurgo_slope_pct <- function(muname_vec) {
  vapply(as.character(muname_vec), function(m) {
    if (is.na(m) || !nzchar(m)) return(NA_real_)
    ml <- tolower(m)
    reg_m <- regmatches(ml, regexec("([0-9]+)[[:space:]]*to[[:space:]]*([0-9]+)[[:space:]]*percent[[:space:]]*slope", ml))[[1L]]
    if (length(reg_m) == 3L) return(mean(c(as.numeric(reg_m[2L]), as.numeric(reg_m[3L]))))
    if (grepl("very steep", ml)) return(35.0)
    if (grepl("hilly|steep", ml)) return(18.0)
    if (grepl("dissected|rolling", ml)) return(12.0)
    if (grepl("undulating|gently sloping", ml)) return(4.0)
    if (grepl("nearly level|level|muck|flooded|lake plain", ml)) return(1.0)
    if (grepl("rock outcrop|flatrock|glade", ml)) return(6.0)
    NA_real_
  }, numeric(1), USE.NAMES = FALSE)
}

classify_ssurgo_outcrop <- function(
  muname = "", compname = "", taxsuborder = "", taxgreatgroup = "", compkind = "",
  lat = NA, lon = NA, is_demo = FALSE, demo = NULL, reskind = "",
  soil_ph = NA_real_, soil_sand = NA_real_, soil_cec = NA_real_,
  bedrock_depth_cm = NA_real_, rock_outcrop_pct = NA_real_
) {
  eff_demo <- if (!is.null(demo)) isTRUE(demo) else isTRUE(is_demo)
  mu_cn <- paste(muname, compname)
  desc <- tolower(trimws(paste(muname, compname, taxsuborder, taxgreatgroup, compkind, reskind)))
  is_lithic_series <- grepl(.ssurgo_outcrop_series_pat(), mu_cn, ignore.case = TRUE)
  is_misc_rock <- grepl("miscellaneous area", tolower(as.character(compkind))) && grepl("rock|rubble|outcrop|boulder", desc)

  ph_num <- suppressWarnings(as.numeric(soil_ph)[1L]); sand_num <- suppressWarnings(as.numeric(soil_sand)[1L])
  cec_num <- suppressWarnings(as.numeric(soil_cec)[1L]); bed_num <- suppressWarnings(as.numeric(bedrock_depth_cm)[1L])
  ro_num <- suppressWarnings(as.numeric(rock_outcrop_pct)[1L])
  has_shallow_phys <- (is.finite(ro_num) && ro_num >= 15.0) || (is.finite(bed_num) && bed_num <= 50.0)
  lith_pat <- "(^|[^[:alnum:]_])lithic($|[^[:alnum:]_])"
  has_lithic_word <- grepl(lith_pat, tolower(paste(muname, compname, taxsuborder, taxgreatgroup))) ||
                     (grepl(lith_pat, tolower(as.character(reskind))) && (!is.finite(bed_num) || bed_num <= 50.0))

  is_serp_txt  <- grepl("chrome|trego|conowingo|aldino|serpentine|ultramafic|peridotite|dunite|soapstone", desc)
  is_mafic_txt <- grepl("iredell|mecklenburg|wilkes|enon|poindexter|diabase|gabbro|basalt|metabasalt|amphibolite|trap rock|greenstone|hornblende|mafic", desc)
  is_lime_txt  <- grepl("gladeville|barfield|cedars|talbott|colbert|limestone|dolomite|dolostone|calcareous|alvar|karst|marble|chalk|rendolls|calciustolls|calciudolls|cedar glade|limestone glade", desc)
  is_sand_txt  <- grepl("ramsey|lily|sewanee|goldston|manteo|wallen|dekab|berks|weikert|channery|shale barren|sandstone|quartzite|siliceous|chert|novaculite|shale|argillite", desc)
  is_gran_txt  <- grepl("wake|louisburg|flatrock|granit|gneiss|felsic|catao|hempstead|ashlar|rion|wateree|pegmatite|rhyolite", desc) ||
                  (grepl("appling", desc) && grepl("flatrock", desc))

  demo_gran <- FALSE; demo_lime <- FALSE
  if (isTRUE(eff_demo)) {
    la <- suppressWarnings(as.numeric(lat)[1L]); lo <- suppressWarnings(as.numeric(lon)[1L])
    demo_gran <- grepl("dystrudepts|rock outcrop", desc) && !is.na(lo) && lo > -84.5 && lo < -77.5 && !is.na(la) && la > 32.5 && la < 37.0
    demo_lime <- grepl("hapludolls|rendolls", desc) && !is.na(lo) && lo > -88.5 && lo < -83.0 && !is.na(la) && la > 34.5 && la < 37.5
  }

  is_outcrop <- grepl("rock outcrop|rock land|flatrock|glade|shale barren|barren|rubble|boulder|weikert|talus|cliff", desc) ||
                has_lithic_word || is_lithic_series || is_misc_rock || has_shallow_phys || is_serp_txt || demo_gran || demo_lime

  outcrop_class <- if (!is_outcrop) {
    "Non-outcrop matrix"
  } else if (is_serp_txt) {
    "Serpentine"
  } else if (is_mafic_txt) {
    "Diabase/Mafic glade"
  } else if (is_lime_txt || demo_lime) {
    "Limestone glade"
  } else if (is_sand_txt) {
    "Sandstone/Siliceous"
  } else if (is_gran_txt || demo_gran) {
    "Granitic flatrock"
  } else if (is.finite(ph_num) && ph_num >= 6.8) {
    "Limestone glade"
  } else if (is.finite(cec_num) && cec_num >= 20.0 && (!is.finite(sand_num) || sand_num < 50.0)) {
    "Diabase/Mafic glade"
  } else if (is.finite(sand_num) && sand_num >= 60.0 && (!is.finite(ph_num) || ph_num <= 5.5)) {
    "Granitic flatrock"
  } else {
    "Sandstone/Siliceous"
  }
  list(is_rock_outcrop = is_outcrop, outcrop_class = outcrop_class)
}

enrich_ssurgo_pedology_from_series <- function(d, demo = FALSE) {
  n <- nrow(d)
  if (n == 0L) return(d)
  num_cols <- c("ssurgo_rock_outcrop_pct", "ssurgo_lithic_series_pct", "ssurgo_bedrock_depth_cm",
                "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm")
  for (col in num_cols) d[[col]] <- if (col %in% names(d)) suppressWarnings(as.numeric(d[[col]])) else rep(NA_real_, n)

  if (!isTRUE(demo)) {
    ok_b <- is.finite(d$ssurgo_bedrock_depth_cm); if (any(ok_b)) d$ssurgo_bedrock_depth_cm[ok_b] <- pmin(200.0, pmax(0.0, d$ssurgo_bedrock_depth_cm[ok_b]))
    ok_r <- is.finite(d$ssurgo_rock_outcrop_pct); if (any(ok_r)) d$ssurgo_rock_outcrop_pct[ok_r] <- pmin(100.0, pmax(0.0, d$ssurgo_rock_outcrop_pct[ok_r]))
    ok_l <- is.finite(d$ssurgo_lithic_series_pct); if (any(ok_l)) d$ssurgo_lithic_series_pct[ok_l] <- pmin(100.0, pmax(0.0, d$ssurgo_lithic_series_pct[ok_l]))
    mu_v <- if ("muname" %in% names(d)) ifelse(is.na(d$muname), "", trimws(as.character(d$muname))) else rep("", n)
    ms_v <- if ("musym" %in% names(d)) ifelse(is.na(d$musym), "", trimws(as.character(d$musym))) else rep("", n)
    need_ro <- is.na(d$ssurgo_rock_outcrop_pct) & (nzchar(mu_v) | nzchar(ms_v))
    if (any(need_ro)) {
      mu_lc <- tolower(mu_v[need_ro])
      d$ssurgo_rock_outcrop_pct[need_ro] <- ifelse(grepl("^rock outcrop$|^rock land$", mu_lc), 85.0, ifelse(grepl("rock outcrop|rock land", mu_lc), 45.0, 0.0))
    }
    return(d)
  }

  # Explicit demonstration mode synthesis ONLY (demo == TRUE)
  lith_pat <- "(^|[^[:alnum:]_])lithic($|[^[:alnum:]_])"
  for (i in seq_len(n)) {
    muname <- if ("muname" %in% names(d) && !is.na(d$muname[i])) trimws(as.character(d$muname[i])) else ""
    musym  <- if ("musym" %in% names(d) && !is.na(d$musym[i])) trimws(as.character(d$musym[i])) else ""
    subord <- if ("taxsuborder" %in% names(d) && !is.na(d$taxsuborder[i])) as.character(d$taxsuborder[i]) else ""
    grtgrp <- if ("taxgreatgroup" %in% names(d) && !is.na(d$taxgreatgroup[i])) as.character(d$taxgreatgroup[i]) else ""
    is_out <- if ("is_rock_outcrop" %in% names(d) && !is.na(d$is_rock_outcrop[i])) isTRUE(d$is_rock_outcrop[i]) else FALSE
    ocl    <- if ("outcrop_class" %in% names(d) && !is.na(d$outcrop_class[i])) as.character(d$outcrop_class[i]) else "Non-outcrop matrix"
    if (!nzchar(muname) && !nzchar(musym) && !nzchar(subord) && !nzchar(grtgrp) && !is_out) next

    desc <- paste(muname, subord, grtgrp); desc_lc <- tolower(desc); mu_lc <- tolower(muname)
    is_lith <- grepl(lith_pat, desc_lc) || grepl(.ssurgo_outcrop_series_pat(), desc, ignore.case = TRUE)
    if (is.na(d$ssurgo_rock_outcrop_pct[i])) d$ssurgo_rock_outcrop_pct[i] <- if (grepl("^rock outcrop$|^rock land$", mu_lc)) 85.0 else if (grepl("rock outcrop|rock land", mu_lc)) 45.0 else 0.0
    if (is.na(d$ssurgo_lithic_series_pct[i])) d$ssurgo_lithic_series_pct[i] <- if (is_lith) 55.0 else 0.0
    if (is.na(d$ssurgo_bedrock_depth_cm[i])) {
      d$ssurgo_bedrock_depth_cm[i] <- if (grepl("^rock outcrop$|^rock land$", mu_lc)) 6.0 else if (grepl("rock outcrop|flatrock|gladeville", desc_lc)) 22.0 else if (grepl(paste0(lith_pat, "|wake|ashlar|goldston|manteo|ramsey|lily|chrome|weikert"), desc_lc)) 36.0 else if (grepl("louisburg|rion|wateree|wilkes|poindexter|enon|iredell|mecklenburg|diabase|granite", desc_lc)) 64.0 else if (is_out) 38.0 else if (grepl("pacolet|appling|cecil|wedowee|madison|hapludults|kanhapludults", desc_lc)) 155.0 else 200.0
    }
    d$ssurgo_bedrock_depth_cm[i] <- pmin(200.0, pmax(0.0, d$ssurgo_bedrock_depth_cm[i]))
    est <- if (ocl == "Granitic flatrock" || grepl("rock outcrop|wake|louisburg|ashlar|rion|wateree|granite|flatrock|dystrudepts", desc_lc)) {
      c(74.0, 8.5, 4.75, 1.8, 4.2, 2.4)
    } else if (grepl("iredell|mecklenburg|wilkes|poindexter|enon|diabase|amphibolite|metabasalt|greenstone|hornblende", desc_lc)) {
      c(38.0, 34.0, 6.35, 2.2, 21.5, 5.2)
    } else if (ocl == "Limestone glade" || grepl("gladeville|barfield|talbott|limestone|dolomite|rendolls|hapludolls", desc_lc)) {
      c(26.0, 36.0, 7.25, 3.8, 26.0, 3.1)
    } else if (ocl == "Serpentine" || grepl("chrome|serpentine|ultramafic", desc_lc)) {
      c(34.0, 28.0, 6.65, 2.6, 24.0, 4.0)
    } else if (ocl == "Sandstone/Siliceous" || grepl("ramsey|lily|sewanee|sandstone|goldston|manteo|berks|weikert|channery|shale barren", desc_lc)) {
      c(68.0, 11.0, 4.85, 1.9, 5.1, 3.2)
    } else if (grepl("muck|saprists|fibrists|hemists|houghton|carbondale|lupton|tawas|wasda|deford", desc_lc)) {
      c(22.0, 14.0, 5.80, 42.0, 36.0, 16.5)
    } else if (grepl("psamments|quartzipsamments|udipsamments|lakeland|wagram|plainfield|coloma|ottokee|corolla|duckston|bonneau|ocilla|spinks|boyer|dune", desc_lc)) {
      c(85.0, 5.5, 5.05, 1.1, 3.9, 3.8)
    } else if (grepl("olls|udolls|aquolls|argiaquolls|endoaquolls|tama|pewamo|sloan|cohoctah|granby", desc_lc)) {
      c(28.0, 26.0, 6.75, 4.2, 21.5, 10.2)
    } else if (grepl("sandy clay loam|clay loam|cecil|pacolet", desc_lc)) {
      c(46.0, 28.0, 5.30, 1.6, 8.2, 7.8)
    } else if (grepl("fine sandy loam|sandy loam|loamy sand|appling|goldsboro|kempsville|slagle|myatt|tomotley|nansemond|hurlock|fallsington", desc_lc)) {
      c(66.0, 13.0, 5.20, 1.5, 6.2, 6.5)
    } else if (grepl("silt loam", desc_lc)) {
      c(24.0, 21.0, 6.10, 2.8, 14.5, 9.5)
    } else {
      c(52.0, 19.0, 5.65, 2.0, 10.5, 7.2)
    }
    if (is.na(d$ssurgo_sand_pct[i])) d$ssurgo_sand_pct[i] <- est[1L]; if (is.na(d$ssurgo_clay_pct[i])) d$ssurgo_clay_pct[i] <- est[2L]
    if (is.na(d$ssurgo_ph[i])) d$ssurgo_ph[i] <- est[3L]; if (is.na(d$ssurgo_om_pct[i])) d$ssurgo_om_pct[i] <- est[4L]
    if (is.na(d$ssurgo_cec[i])) d$ssurgo_cec[i] <- est[5L]; if (is.na(d$ssurgo_aws_0_50cm[i])) d$ssurgo_aws_0_50cm[i] <- est[6L]
  }
  d
}

build_ssurgo_buffered_sda_sql <- function(sub_pts, buffer_deg = 0.00045) {
  # Tight 30-50 m focal positional tolerance (0.00027 to 0.00045 deg); never expands to 1 km
  tol_deg <- suppressWarnings(as.numeric(buffer_deg)[1L])
  if (!is.finite(tol_deg) || tol_deg <= 0) tol_deg <- 0.00045
  tol_deg <- pmin(0.00045, pmax(0.00027, tol_deg))
  union_rows <- vapply(seq_len(nrow(sub_pts)), function(i) sprintf("SELECT %d AS gid, %.6f AS lon, %.6f AS lat, %.5f AS buf_deg", as.integer(sub_pts$gid[i]), as.numeric(sub_pts$longitude[i]), as.numeric(sub_pts$latitude[i]), tol_deg), character(1))
  sprintf(
    paste0("-- WKT ST_Buffer(ST_Point(lon, lat), %.5f) (tight 30-50m focal positional tolerance)\n",
           "WITH pts AS (%s), buffered_pts AS (SELECT gid, lon, lat, geometry::STGeomFromText('POINT(' + CAST(lon AS VARCHAR(32)) + ' ' + CAST(lat AS VARCHAR(32)) + ')', 4326).STBuffer(buf_deg) AS geom FROM pts),\n",
           "mu_hits AS (SELECT DISTINCT p.gid, mp.mukey FROM buffered_pts p INNER JOIN mupolygon mp ON mp.mupolygongeo.STIntersects(p.geom) = 1)\n",
           "SELECT h.gid, h.mukey, m.musym, m.muname, ma.brockdepmin, ma.aws050wta, c.cokey, c.compname, c.comppct_r, c.compkind, c.majcompflag, c.taxsuborder, c.taxgrtgroup AS taxgreatgroup, c.taxsubgrp, ",
           "ch.hzdept_r, ch.hzdepb_r, ch.sandtotal_r, ch.claytotal_r, ch.ph1to1h2o_r, ch.om_r, ch.cec7_r, ch.awc_r, cr.resdept_r, cr.reskind ",
           "FROM mu_hits h INNER JOIN mapunit m ON h.mukey = m.mukey LEFT JOIN muaggatt ma ON m.mukey = ma.mukey ",
           "LEFT JOIN component c ON m.mukey = c.mukey LEFT JOIN chorizon ch ON c.cokey = ch.cokey AND ch.hzdept_r < 15 LEFT JOIN corestrictions cr ON c.cokey = cr.cokey"),
    tol_deg, paste(union_rows, collapse = " UNION ALL ")
  )
}

aggregate_ssurgo_tabular_to_points <- function(sub_pts, tab_res, classify_fn = NULL, is_demo = FALSE) {
  if (is.null(classify_fn) || !is.function(classify_fn)) classify_fn <- classify_ssurgo_outcrop
  fn_formals <- names(formals(classify_fn)); has_dots <- "..." %in% fn_formals
  n_sub <- nrow(sub_pts); out_list <- vector("list", n_sub); lith_pat <- "(^|[^[:alnum:]_])lithic($|[^[:alnum:]_])"

  for (i in seq_len(n_sub)) {
    gid_i <- sub_pts$gid[i]; lat_i <- sub_pts$latitude[i]; lon_i <- sub_pts$longitude[i]; cat_i <- as.character(sub_pts$catalogNumber[i])
    rows_i <- if (!is.null(tab_res) && nrow(tab_res) > 0L && "gid" %in% names(tab_res)) tab_res[!is.na(tab_res$gid) & as.integer(tab_res$gid) == as.integer(gid_i), , drop = FALSE] else data.frame()
    if (nrow(rows_i) == 0L) {
      out_list[[i]] <- data.frame(catalogNumber = cat_i, latitude = lat_i, longitude = lon_i, musym = NA_character_, muname = NA_character_, taxsuborder = NA_character_, taxgreatgroup = NA_character_, is_rock_outcrop = FALSE, outcrop_class = "Non-outcrop matrix", ssurgo_rock_outcrop_pct = NA_real_, ssurgo_lithic_series_pct = NA_real_, ssurgo_bedrock_depth_cm = NA_real_, ssurgo_sand_pct = NA_real_, ssurgo_clay_pct = NA_real_, ssurgo_ph = NA_real_, ssurgo_om_pct = NA_real_, ssurgo_cec = NA_real_, ssurgo_aws_0_50cm = NA_real_, stringsAsFactors = FALSE)
      next
    }

    mk_chr <- ifelse(is.na(rows_i$mukey), "", as.character(rows_i$mukey))
    ck_raw <- if ("cokey" %in% names(rows_i)) as.character(rows_i$cokey) else paste0(mk_chr, "_", rows_i$compname)
    rows_i$cokey_id <- ifelse(is.na(ck_raw) | !nzchar(ck_raw), paste0("r_", seq_len(nrow(rows_i))), ck_raw)
    comp_df <- rows_i[!duplicated(rows_i$cokey_id), , drop = FALSE]
    comp_mk <- ifelse(is.na(comp_df$mukey), "", as.character(comp_df$mukey))
    comp_pct_w <- suppressWarnings(as.numeric(comp_df$comppct_r)); comp_pct_w[!is.finite(comp_pct_w) | comp_pct_w <= 0] <- 15.0

    mu_lc_all <- tolower(trimws(ifelse(is.na(comp_df$muname), "", as.character(comp_df$muname))))
    focal_mk <- comp_mk[order(ifelse(grepl("rock outcrop|rock land|flatrock|glade", mu_lc_all), 1L, 0L), ifelse(!is.na(comp_df$majcompflag) & comp_df$majcompflag == "Yes", 1L, 0L), comp_pct_w, decreasing = TRUE)[1L]]
    focal_rows <- rows_i[mk_chr == focal_mk, , drop = FALSE]; focal_comp <- comp_df[comp_mk == focal_mk, , drop = FALSE]
    f_pct_raw <- suppressWarnings(as.numeric(focal_comp$comppct_r)); f_has_pct <- is.finite(f_pct_raw) & f_pct_raw > 0; f_pct_w <- ifelse(f_has_pct, f_pct_raw, 15.0)

    cn_lc <- tolower(trimws(ifelse(is.na(focal_comp$compname), "", as.character(focal_comp$compname))))
    ck_lc <- tolower(trimws(ifelse(is.na(focal_comp$compkind), "", as.character(focal_comp$compkind))))
    mu_lc_f <- tolower(trimws(ifelse(is.na(focal_comp$muname), "", as.character(focal_comp$muname))))
    mu_exp <- any(grepl("rock outcrop|rock land|flatrock|glade", mu_lc_f))

    # Separate exposed rock outcrop miscellaneous areas from shallow/Lithic soil series
    is_rock_c <- (ck_lc == "miscellaneous area" & grepl("rock|rubble|outcrop|boulder", cn_lc)) | grepl("(^|[^[:alnum:]_])(rock outcrop|rock land)($|[^[:alnum:]_])", cn_lc)
    tsg <- if ("taxsubgrp" %in% names(focal_comp)) as.character(focal_comp$taxsubgrp) else rep("", nrow(focal_comp))
    tgg <- if ("taxgreatgroup" %in% names(focal_comp)) as.character(focal_comp$taxgreatgroup) else rep("", nrow(focal_comp))
    comp_lith_r <- vapply(seq_len(nrow(focal_comp)), function(ci) {
      sr <- focal_rows[focal_rows$cokey_id == focal_comp$cokey_id[ci], , drop = FALSE]
      rd <- if ("resdept_r" %in% names(sr)) suppressWarnings(as.numeric(sr$resdept_r)) else numeric(0)
      rk <- if ("reskind" %in% names(sr)) tolower(ifelse(is.na(sr$reskind), "", as.character(sr$reskind))) else character(0)
      any(is.finite(rd) & rd >= 0 & rd <= 50.0 & grepl(lith_pat, rk))
    }, logical(1))
    is_lith_c <- (!is_rock_c) & (grepl(lith_pat, tolower(paste(tsg, tgg))) | comp_lith_r)

    rock_pct <- if (any(is_rock_c)) { if (any(f_has_pct[is_rock_c])) pmin(100.0, sum(f_pct_raw[is_rock_c & f_has_pct], na.rm = TRUE)) else NA_real_ } else if (all(!f_has_pct) && any(grepl("rock outcrop|rock land", mu_lc_f))) NA_real_ else 0.0
    lith_pct <- if (any(is_lith_c)) { if (any(f_has_pct[is_lith_c])) pmin(100.0, sum(f_pct_raw[is_lith_c & f_has_pct], na.rm = TRUE)) else NA_real_ } else 0.0

    # Prevent < 25% minor inclusions from hijacking rep_row unless muname explicitly names Rock outcrop/flatrock/glade
    elig_out <- (is_rock_c | is_lith_c) & (mu_exp | (f_has_pct & f_pct_raw >= 25.0))
    rep_pool <- if (any(elig_out)) focal_comp[elig_out, , drop = FALSE] else focal_comp
    rep_row <- rep_pool[order(ifelse(!is.na(rep_pool$majcompflag) & rep_pool$majcompflag == "Yes", 1L, 0L), suppressWarnings(as.numeric(rep_pool$comppct_r)), decreasing = TRUE)[1L], , drop = FALSE]

    maj_ck <- focal_comp$cokey_id[mu_exp | (!f_has_pct) | (f_pct_raw >= 25.0)]
    if (length(maj_ck) == 0L) maj_ck <- focal_comp$cokey_id
    bed_rows <- focal_rows[focal_rows$cokey_id %in% maj_ck, , drop = FALSE]
    res_c <- if ("resdept_r" %in% names(bed_rows)) suppressWarnings(as.numeric(bed_rows$resdept_r)) else numeric(0)
    res_c <- res_c[is.finite(res_c) & res_c >= 0]
    brk_m <- if ("brockdepmin" %in% names(bed_rows)) suppressWarnings(as.numeric(bed_rows$brockdepmin)) else numeric(0)
    brk_m <- brk_m[is.finite(brk_m) & brk_m >= 0]
    res_d <- if (mu_exp || length(res_c) == 0L) c(res_c, brk_m) else res_c
    bed_cm <- if (length(res_d) > 0L) min(200.0, min(res_d)) else NA_real_

    hz_t <- if ("hzdept_r" %in% names(focal_rows)) suppressWarnings(as.numeric(focal_rows$hzdept_r)) else rep(0.0, nrow(focal_rows))
    hz_b <- if ("hzdepb_r" %in% names(focal_rows)) suppressWarnings(as.numeric(focal_rows$hzdepb_r)) else hz_t + 15.0
    hz_t[is.na(hz_t)] <- 0.0; hz_b[is.na(hz_b) | hz_b <= hz_t] <- pmin(15.0, hz_t[is.na(hz_b) | hz_b <= hz_t] + 10.0)
    w_c <- suppressWarnings(as.numeric(focal_rows$comppct_r)); w_c[is.na(w_c) | w_c <= 0] <- 10.0
    ov1 <- pmax(0.0, pmin(hz_b, 5.0) - pmax(hz_t, 0.0)); ov2 <- pmax(0.0, pmin(hz_b, 15.0) - pmax(hz_t, 5.0))
    w_rhizo <- function(vr) {
      v <- suppressWarnings(as.numeric(vr)); o1 <- is.finite(v) & is.finite(w_c) & ov1 > 0; o2 <- is.finite(v) & is.finite(w_c) & ov2 > 0
      res <- compute_rhizosphere_0_15cm(if (any(o1)) sum(v[o1] * w_c[o1] * ov1[o1]) / sum(w_c[o1] * ov1[o1]) else NA_real_, if (any(o2)) sum(v[o2] * w_c[o2] * ov2[o2]) / sum(w_c[o2] * ov2[o2]) else NA_real_)
      if (!is.na(res)) return(res)
      oa <- is.finite(v) & is.finite(w_c) & w_c > 0; if (!any(oa)) NA_real_ else sum(v[oa] * w_c[oa]) / sum(w_c[oa])
    }
    w_avg <- function(vr, wv) { v <- suppressWarnings(as.numeric(vr)); ok <- is.finite(v) & is.finite(wv) & wv > 0; if (!any(ok)) NA_real_ else sum(v[ok] * wv[ok]) / sum(wv[ok]) }

    sand_p <- w_rhizo(focal_rows$sandtotal_r); clay_p <- w_rhizo(focal_rows$claytotal_r); ph_v <- w_rhizo(focal_rows$ph1to1h2o_r)
    om_p <- w_rhizo(focal_rows$om_r); cec_v <- w_rhizo(focal_rows$cec7_r); aws_v <- w_avg(focal_comp$aws050wta, f_pct_w)
    rk_val <- if ("reskind" %in% names(bed_rows)) { rkv <- stats::na.omit(unique(trimws(as.character(bed_rows$reskind)))); paste(rkv[nzchar(rkv)], collapse = "; ") } else ""

    c_args <- list(muname = ifelse(!is.na(rep_row$muname[1L]), as.character(rep_row$muname[1L]), ""), compname = ifelse(!is.na(rep_row$compname[1L]), as.character(rep_row$compname[1L]), ""), taxsuborder = ifelse(!is.na(rep_row$taxsuborder[1L]), as.character(rep_row$taxsuborder[1L]), ""), taxgreatgroup = ifelse(!is.na(rep_row$taxgreatgroup[1L]), as.character(rep_row$taxgreatgroup[1L]), ""), lat = lat_i, lon = lon_i)
    x_args <- list(compkind = if ("compkind" %in% names(rep_row) && !is.na(rep_row$compkind[1L])) as.character(rep_row$compkind[1L]) else "", is_demo = isTRUE(is_demo), reskind = rk_val, soil_ph = ph_v, soil_sand = sand_p, soil_cec = cec_v, bedrock_depth_cm = bed_cm, rock_outcrop_pct = rock_pct)
    for (nm in names(x_args)) if (has_dots || (nm %in% fn_formals)) c_args[[nm]] <- x_args[[nm]]
    cl <- do.call(classify_fn, c_args)
    is_of <- isTRUE(cl$is_rock_outcrop) || (is.finite(rock_pct) && rock_pct >= 15.0) || (is.finite(lith_pct) && lith_pct >= 25.0)
    ocl_f <- if (is_of && cl$outcrop_class == "Non-outcrop matrix") {
      if (is.finite(ph_v) && ph_v >= 6.8) { "Limestone glade" } else if (is.finite(cec_v) && cec_v >= 20.0 && (!is.finite(sand_p) || sand_p < 50.0)) { "Diabase/Mafic glade" } else if (is.finite(sand_p) && sand_p >= 60.0 && (!is.finite(ph_v) || ph_v <= 5.5)) { "Granitic flatrock" } else { "Sandstone/Siliceous" }
    } else cl$outcrop_class

    out_list[[i]] <- enrich_ssurgo_pedology_from_series(data.frame(
      catalogNumber = cat_i, latitude = lat_i, longitude = lon_i, musym = as.character(rep_row$musym[1L]), muname = as.character(rep_row$muname[1L]),
      taxsuborder = as.character(rep_row$taxsuborder[1L]), taxgreatgroup = as.character(rep_row$taxgreatgroup[1L]), is_rock_outcrop = is_of, outcrop_class = ocl_f,
      ssurgo_rock_outcrop_pct = if (is.finite(rock_pct)) round(rock_pct, 2) else NA_real_, ssurgo_lithic_series_pct = if (is.finite(lith_pct)) round(lith_pct, 2) else NA_real_,
      ssurgo_bedrock_depth_cm = if (is.finite(bed_cm)) round(bed_cm, 1) else NA_real_, ssurgo_sand_pct = round(sand_p, 2), ssurgo_clay_pct = round(clay_p, 2),
      ssurgo_ph = round(ph_v, 2), ssurgo_om_pct = round(om_p, 2), ssurgo_cec = round(cec_v, 2), ssurgo_aws_0_50cm = round(aws_v, 2), stringsAsFactors = FALSE
    ), demo = isTRUE(is_demo))
  }
  dplyr::bind_rows(out_list)
}

normalize_ssurgo_cols <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  for (col in c("catalogNumber", "musym", "muname", "taxsuborder", "taxgreatgroup", "outcrop_class")) d[[col]] <- if (col %in% names(d)) as.character(d[[col]]) else if (col == "outcrop_class") "Non-outcrop matrix" else NA_character_
  for (col in c("latitude", "longitude")) d[[col]] <- if (col %in% names(d)) suppressWarnings(as.numeric(d[[col]])) else NA_real_
  d$is_rock_outcrop <- if ("is_rock_outcrop" %in% names(d)) as.logical(d$is_rock_outcrop) else FALSE; d$is_rock_outcrop[is.na(d$is_rock_outcrop)] <- FALSE
  d$outcrop_class[is.na(d$outcrop_class) | !nzchar(d$outcrop_class)] <- "Non-outcrop matrix"
  enrich_ssurgo_pedology_from_series(d)[, SSURGO_CANONICAL_COLS, drop = FALSE]
}

.hydrate_ssurgo_from_sqlite <- function(pts_df, cache_dir, demo = FALSE) {
  if (isTRUE(demo) || !exists("query_environmental_cache_sqlite", mode = "function") || nrow(pts_df) == 0L) return(list(hit_df = NULL, hit_idx = integer(0)))
  p_ll <- sprintf("%.5f_%.5f", round(as.numeric(pts_df$latitude), 5), round(as.numeric(pts_df$longitude), 5))
  sql_ss <- tryCatch(query_environmental_cache_sqlite(catalog_numbers = as.character(pts_df$catalogNumber), lat_lon_keys = p_ll, cache_path = cache_dir, demo = demo), error = function(e) NULL)
  if (is.null(sql_ss) || nrow(sql_ss) == 0L || !any(c("musym", "ssurgo_rock_outcrop_pct", "outcrop_class") %in% names(sql_ss))) return(list(hit_df = NULL, hit_idx = integer(0)))
  mc <- if ("catalogNumber" %in% names(sql_ss)) match(as.character(pts_df$catalogNumber), as.character(sql_ss$catalogNumber)) else rep(NA_integer_, nrow(pts_df))
  ml <- if ("lat_lon_key" %in% names(sql_ss)) match(p_ll, as.character(sql_ss$lat_lon_key)) else rep(NA_integer_, nrow(pts_df))
  mh <- ifelse(!is.na(mc), mc, ml); has_h <- !is.na(mh)
  if (!any(has_h)) return(list(hit_df = NULL, hit_idx = integer(0)))
  cand <- sql_ss[mh[has_h], , drop = FALSE]
  ok_s <- (if ("musym" %in% names(cand)) !is.na(cand$musym) & nzchar(as.character(cand$musym)) else rep(FALSE, nrow(cand))) | (if ("ssurgo_rock_outcrop_pct" %in% names(cand)) !is.na(cand$ssurgo_rock_outcrop_pct) else rep(FALSE, nrow(cand)))
  if (!any(ok_s)) return(list(hit_df = NULL, hit_idx = integer(0)))
  hi <- which(has_h)[ok_s]; hr <- cand[ok_s, , drop = FALSE]
  hr$catalogNumber <- as.character(pts_df$catalogNumber[hi]); hr$latitude <- pts_df$latitude[hi]; hr$longitude <- pts_df$longitude[hi]
  list(hit_df = normalize_ssurgo_cols(hr), hit_idx = hi)
}

query_ssurgo_sda_batch <- function(batch_df, cache_dir, batch_idx, max_retries = 4, demo = FALSE, batch_size = 75) {
  bh <- if (requireNamespace("digest", quietly = TRUE)) substr(digest::digest(paste(batch_df$catalogNumber, batch_df$latitude, batch_df$longitude, collapse = ";"), algo = "md5"), 1L, 12L) else sprintf("b%04d", batch_idx)
  c_rds <- file.path(cache_dir, sprintf("ssurgo_batch_%03d_%s.rds", batch_idx, bh))
  c_csv <- file.path(cache_dir, sprintf("ssurgo_batch_%03d_%s.csv", batch_idx, bh))
  if (file.exists(c_rds)) return(normalize_ssurgo_cols(readRDS(c_rds)))
  if (file.exists(c_csv)) return(normalize_ssurgo_cols(as.data.frame(readr::read_csv(c_csv, show_col_types = FALSE))))
  hyd <- .hydrate_ssurgo_from_sqlite(batch_df, cache_dir, demo = demo)
  if (!is.null(hyd$hit_df) && length(hyd$hit_idx) == nrow(batch_df)) return(hyd$hit_df)

  pts_df <- batch_df[, c("catalogNumber", "longitude", "latitude")]
  n_pts <- nrow(pts_df); query_ok <- FALSE; res_df <- NULL; delays <- c(1, 2, 4, 8)
  if (!isTRUE(demo) && requireNamespace("soilDB", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
    chunks <- split(seq_len(n_pts), ceiling(seq_len(n_pts) / max(1L, min(as.integer(batch_size), 100L))))
    c_res <- list(); all_ok <- TRUE
    for (ci in seq_along(chunks)) {
      sub_pts <- pts_df[chunks[[ci]], , drop = FALSE]; sub_pts$gid <- seq_len(nrow(sub_pts)); ch_ok <- FALSE; ch_df <- NULL
      for (att in seq_len(max_retries)) {
        tryCatch({
          bt <- tryCatch(soilDB::SDA_query(build_ssurgo_buffered_sda_sql(sub_pts, buffer_deg = 0.00045)), error = function(e) NULL)
          if (is.data.frame(bt) && nrow(bt) > 0L && "mukey" %in% names(bt)) { ch_df <- normalize_ssurgo_cols(aggregate_ssurgo_tabular_to_points(sub_pts, bt, classify_fn = classify_ssurgo_outcrop, is_demo = FALSE)); ch_ok <- TRUE; break }
          pts_sf <- sf::st_as_sf(sub_pts, coords = c("longitude", "latitude"), crs = 4326)
          sp_res <- tryCatch(soilDB::SDA_spatialQuery(tryCatch(suppressWarnings(sf::st_buffer(pts_sf, dist = 0.00045)), error = function(e) pts_sf), what = "mukey", db = "SSURGO", byFeature = TRUE), error = function(e_sp) {
            if (!grepl("numbers of columns", conditionMessage(e_sp), fixed = TRUE)) stop(e_sp)
            pl <- Filter(Negate(is.null), lapply(seq_len(nrow(sub_pts)), function(p_i) { rp <- tryCatch(soilDB::SDA_spatialQuery(pts_sf[p_i, , drop = FALSE], what = "mukey", db = "SSURGO", byFeature = FALSE), error = function(e) NULL); if (is.data.frame(rp) && nrow(rp) > 0L && "mukey" %in% names(rp)) { rp$gid <- sub_pts$gid[p_i]; rp } else NULL }))
            if (length(pl) > 0L) dplyr::bind_rows(pl) else NULL
          })
          if (is.data.frame(sp_res) && nrow(sp_res) > 0L && "mukey" %in% names(sp_res)) {
            mks <- unique(stats::na.omit(as.character(sp_res$mukey)))
            if (length(mks) > 0L) {
              sql_b <- sprintf("SELECT m.mukey, m.musym, m.muname, ma.brockdepmin, ma.aws050wta, c.cokey, c.compname, c.comppct_r, c.compkind, c.majcompflag, c.taxsuborder, c.taxgrtgroup AS taxgreatgroup, c.taxsubgrp, cr.resdept_r, cr.reskind, ch.hzdept_r, ch.hzdepb_r, ch.sandtotal_r, ch.claytotal_r, ch.ph1to1h2o_r, ch.om_r, ch.cec7_r FROM mapunit m LEFT JOIN muaggatt ma ON m.mukey = ma.mukey LEFT JOIN component c ON m.mukey = c.mukey LEFT JOIN corestrictions cr ON c.cokey = cr.cokey LEFT JOIN chorizon ch ON c.cokey = ch.cokey AND ch.hzdept_r < 20 WHERE m.mukey IN (%s)", paste0("'", mks, "'", collapse = ","))
              tr <- soilDB::SDA_query(sql_b)
              if (is.data.frame(tr) && nrow(tr) > 0L) {
                ch_df <- normalize_ssurgo_cols(aggregate_ssurgo_tabular_to_points(sub_pts, merge(sp_res[, intersect(c("gid", "mukey"), names(sp_res)), drop = FALSE], tr, by = "mukey", all.x = TRUE), classify_fn = classify_ssurgo_outcrop, is_demo = FALSE))
                ch_ok <- TRUE; break
              }
            }
          }
        }, error = function(e) message(sprintf("Batch %d chunk %d attempt %d failed: %s", batch_idx, ci, att, e$message)))
        if (att < max_retries && !ch_ok) Sys.sleep(delays[att])
      }
      if (ch_ok && !is.null(ch_df)) c_res[[ci]] <- ch_df else { all_ok <- FALSE; break }
    }
    if (all_ok && length(c_res) > 0L) { res_df <- normalize_ssurgo_cols(dplyr::bind_rows(c_res)); query_ok <- TRUE }
  }

  if (query_ok && !is.null(res_df)) {
    if (!isTRUE(demo)) { dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE); readr::write_csv(res_df, paste0(c_csv, ".tmp")); file.rename(paste0(c_csv, ".tmp"), c_csv); saveRDS(res_df, paste0(c_rds, ".tmp")); file.rename(paste0(c_rds, ".tmp"), c_rds) }
    return(res_df)
  }
  if (!isTRUE(demo)) {
    warning(sprintf("soilDB::SDA_spatialQuery() failed after %d retries for batch %d; retaining NA SSURGO attributes.", max_retries, batch_idx))
    return(normalize_ssurgo_cols(data.frame(catalogNumber = pts_df$catalogNumber, latitude = pts_df$latitude, longitude = pts_df$longitude, musym = NA_character_, muname = NA_character_, taxsuborder = NA_character_, taxgreatgroup = NA_character_, is_rock_outcrop = FALSE, outcrop_class = "Non-outcrop matrix", stringsAsFactors = FALSE)))
  }

  musyms <- character(n_pts); munames <- character(n_pts); suborders <- character(n_pts); greatgroups <- character(n_pts); is_outcrops <- logical(n_pts); outcrop_classes <- character(n_pts)
  for (i in seq_len(n_pts)) {
    lat <- pts_df$latitude[i]; lon <- pts_df$longitude[i]; reg <- if ("regional_group" %in% names(batch_df)) batch_df$regional_group[i] else "Unknown"
    row_m <- if (grepl("Flatrock", reg, ignore.case = TRUE) || (lon > -84.5 && lon < -78.5 && lat > 33.2 && lat < 36.8)) {
      if (as.integer(abs(sin(lat * 100 + lon * 100)) * 1000) %% 10L < 7L) c("RoB", "Rock outcrop-Wake complex, 2 to 10 percent slopes", "Udepts", "Lithic Dystrudepts") else c("ApB", "Appling sandy loam, 2 to 6 percent slopes", "Udults", "Typic Kanhapludults")
    } else if (grepl("Sandhill|Coastal", reg, ignore.case = TRUE) || (lon > -82.0 && lon < -75.5 && lat < 36.5)) {
      c("WaB", "Wagram sand, 0 to 6 percent slopes", "Udults", "Arenic Kandiudults")
    } else if (grepl("Midwest|Prairie", reg, ignore.case = TRUE) || (lon < -88.0)) {
      c("TaA", "Tama silt loam, 0 to 2 percent slopes", "Udolls", "Typic Argiudolls")
    } else if (grepl("Appalachian", reg, ignore.case = TRUE) || (lon >= -84.0 && lon <= -79.0 && lat >= 35.0 && lat <= 39.0)) {
      if (as.integer(abs(cos(lat * 50 + lon * 50)) * 1000) %% 10L < 3L) c("GvC", "Gladeville-Rock outcrop complex, 2 to 12 percent slopes", "Udolls", "Lithic Hapludolls") else c("CeB", "Cecil sandy clay loam, 2 to 8 percent slopes", "Udults", "Typic Hapludults")
    } else {
      c("MuB", "Generic Upland loam, 1 to 5 percent slopes", "Udults", "Typic Hapludults")
    }
    musyms[i] <- row_m[1L]; munames[i] <- row_m[2L]; suborders[i] <- row_m[3L]; greatgroups[i] <- row_m[4L]
    cl <- classify_ssurgo_outcrop(munames[i], "", suborders[i], greatgroups[i], "", lat, lon, is_demo = TRUE)
    is_outcrops[i] <- cl$is_rock_outcrop; outcrop_classes[i] <- cl$outcrop_class
  }
  normalize_ssurgo_cols(data.frame(catalogNumber = pts_df$catalogNumber, latitude = pts_df$latitude, longitude = pts_df$longitude, musym = musyms, muname = munames, taxsuborder = suborders, taxgreatgroup = greatgroups, is_rock_outcrop = is_outcrops, outcrop_class = outcrop_classes, stringsAsFactors = FALSE))
}

query_ssurgo_sda <- function(df, cache_dir = "data/environmental/ssurgo_cache", batch_size = 75, demo = FALSE) {
  message(sprintf("Initiating USDA SSURGO Vector Pedology SDA Query (cache: %s, batch_size: %d)...", cache_dir, batch_size))
  conus_df <- df[!is.na(df$latitude) & !is.na(df$longitude) & df$latitude >= 24.5 & df$latitude <= 49.5 & df$longitude >= -125.0 & df$longitude <= -66.5, , drop = FALSE]
  q_cols <- c("ssurgo_rock_outcrop_pct", "ssurgo_lithic_series_pct", "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm")
  init_empty <- function(d_in) {
    d_in$musym <- NA_character_; d_in$muname <- NA_character_; d_in$taxsuborder <- NA_character_; d_in$taxgreatgroup <- NA_character_
    d_in$is_rock_outcrop <- FALSE; d_in$outcrop_class <- "Non-outcrop matrix"; for (sc in q_cols) d_in[[sc]] <- NA_real_; d_in
  }
  if (nrow(conus_df) == 0L) return(init_empty(df))

  tryCatch({
    batch_results <- list(); rem_conus <- conus_df
    hyd <- .hydrate_ssurgo_from_sqlite(conus_df, cache_dir, demo = demo)
    if (!is.null(hyd$hit_df) && length(hyd$hit_idx) > 0L) { batch_results[[1L]] <- hyd$hit_df; rem_conus <- conus_df[-hyd$hit_idx, , drop = FALSE] }
    if (nrow(rem_conus) > 0L) {
      b_idx <- split(seq_len(nrow(rem_conus)), ceiling(seq_len(nrow(rem_conus)) / batch_size))
      run_one <- function(bi) query_ssurgo_sda_batch(rem_conus[b_idx[[bi]], , drop = FALSE], cache_dir, bi, demo = demo, batch_size = batch_size)
      rem_b <- if (!isTRUE(demo) && length(b_idx) >= 2L && exists("packera_mclapply", mode = "function")) {
        packera_mclapply(seq_along(b_idx), run_one, mc.cores = if (exists("resolve_packera_workers", mode = "function")) min(16L, resolve_packera_workers(16L, length(b_idx))) else min(4L, length(b_idx)), mc.preschedule = FALSE)
      } else lapply(seq_along(b_idx), run_one)
      batch_results <- c(batch_results, rem_b)
    }
    all_ss <- as.data.frame(dplyr::bind_rows(batch_results), stringsAsFactors = FALSE)
    keep_c <- c("catalogNumber", "musym", "muname", "taxsuborder", "taxgreatgroup", "is_rock_outcrop", "outcrop_class", q_cols)
    all_ss <- all_ss[!duplicated(all_ss$catalogNumber), intersect(names(all_ss), keep_c), drop = FALSE]
    df$catalogNumber <- as.character(df$catalogNumber); all_ss$catalogNumber <- as.character(all_ss$catalogNumber)
    for (col_rm in intersect(setdiff(keep_c, "catalogNumber"), names(df))) df[[col_rm]] <- NULL
    m_df <- dplyr::left_join(df, all_ss, by = "catalogNumber")
    m_df$is_rock_outcrop[is.na(m_df$is_rock_outcrop)] <- FALSE; m_df$outcrop_class[is.na(m_df$outcrop_class)] <- "Non-outcrop matrix"
    m_df
  }, error = function(e) { warning(sprintf("USDA SSURGO SDA query failed (%s); retaining SoilGrids 250m + WorldClim v2.1 raster covariates.", e$message)); init_empty(df) })
}
