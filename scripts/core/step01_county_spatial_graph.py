"""
Module: step01_county_spatial_graph.py
Description:
    US county polygon layer generation, 1-ring topological adjacency graph construction,
    voucher county resolution, and geodesic buffer membership queries.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List, Mapping, Optional, Set, Tuple

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    from shapely.geometry import Point, shape
    from shapely.ops import nearest_points
    from shapely.strtree import STRtree
except ImportError:
    Point = shape = nearest_points = STRtree = None  # type: ignore[assignment]


def _ensure_standalone_us_counties_geojson(cnty_geojson_path: Path) -> bool:
    """Build `us_counties.geojson` locally via R `maps` + `sf` if absent and Rscript is installed."""
    if cnty_geojson_path.is_file() and cnty_geojson_path.stat().st_size > 0:
        return True
    rscript_bin = shutil.which("Rscript")
    if not rscript_bin:
        return False
    cnty_geojson_path.parent.mkdir(parents=True, exist_ok=True)
    r_code = (
        "suppressPackageStartupMessages({library(sf); library(maps)}); sf_use_s2(FALSE);"
        "cm <- st_make_valid(st_as_sf(maps::map('county', fill=TRUE, plot=FALSE))); st_crs(cm) <- 4326;"
        "pts <- strsplit(as.character(cm$ID), ',');"
        "st_n <- vapply(pts, `[`, character(1), 1);"
        "co_n <- trimws(gsub('[^a-z0-9]+', ' ', gsub('^saint\\\\s+|^st\\\\.\\\\s*', 'st ', gsub('\\\\s+city$', '', sub(':.*$', '', vapply(pts, `[`, character(1), 2))))));"
        "ck <- paste0(st_n, '|', co_n); dk <- unique(ck[duplicated(ck)]); sm <- !(ck %in% dk); ga <- st_geometry(cm);"
        "ok <- ck[sm]; os <- st_n[sm]; oc <- co_n[sm]; og <- lapply(which(sm), function(i) st_cast(ga[[i]], 'MULTIPOLYGON'));"
        "for(d in dk){ id <- which(ck==d); ok <- c(ok, d); os <- c(os, st_n[id[1]]); oc <- c(oc, co_n[id[1]]); og[[length(og)+1L]] <- st_cast(st_make_valid(st_union(ga[id]))[[1]], 'MULTIPOLYGON'); };"
        "cd <- st_sf(county_key=ok, state_name=os, county_norm=oc, geometry=st_sfc(og, crs=4326));"
        "cp <- st_coordinates(suppressWarnings(st_point_on_surface(cd)));"
        "cd$centroid_lon <- cp[,1]; cd$centroid_lat <- cp[,2];"
        f"tc <- {json.dumps(str(cnty_geojson_path) + '.tmp.geojson')};"
        "st_write(cd, tc, driver='GeoJSON', quiet=TRUE, delete_dsn=TRUE);"
        "file.rename(tc, sub('\\\\.tmp\\\\.geojson$', '', tc));"
    )
    try:
        subprocess.run([rscript_bin, "-e", r_code], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
    except Exception:
        return False
    return cnty_geojson_path.is_file() and cnty_geojson_path.stat().st_size > 0


def load_us_county_spatial_graph(
    county_path: Optional[Path] = None,
    ecoregion_path: Optional[Path] = None,
    allow_download: bool = True,
) -> Dict[str, Any]:
    """
    Build and cache a spatial index (`STRtree`) and 1-ring topological county
    adjacency map over `us_counties.geojson`.
    """
    import scripts.core.usda_plants_ranges as _usda

    cnty_p = Path(county_path) if county_path is not None else _usda.DEFAULT_US_COUNTIES_GEOJSON
    eco_p = Path(ecoregion_path) if ecoregion_path is not None else _usda.DEFAULT_EPA_L3_GEOJSON

    if not cnty_p.exists() and county_path is None:
        if allow_download:
            _usda._ensure_epa_level3_cache_files(eco_p, cnty_p)
        if not cnty_p.exists():
            _usda._ensure_standalone_us_counties_geojson(cnty_p)

    mtime = cnty_p.stat().st_mtime if cnty_p.exists() else -1.0
    cache_key = (str(cnty_p.resolve()) if cnty_p.exists() else "OFFLINE_FALLBACK", mtime)
    if cache_key in _usda._COUNTY_SPATIAL_GRAPH_CACHE:
        return _usda._COUNTY_SPATIAL_GRAPH_CACHE[cache_key]

    idx_bundle = _usda.load_epa_level3_spatial_index(
        ecoregion_path=eco_p if eco_p.exists() else ecoregion_path,
        county_path=cnty_p if cnty_p.exists() else county_path,
        allow_download=False,
    )
    counties_dict: Dict[str, Dict[str, Any]] = dict(idx_bundle.get("counties") or {})
    # If eco_p did not exist yet, load cnty_p directly if it exists
    if cnty_p.exists() and len(counties_dict) <= 10:
        try:
            cnty_fc = json.loads(cnty_p.read_text(encoding="utf-8"))
            counties_dict = {}
            for feat in cnty_fc.get("features", []):
                props = feat.get("properties") or {}
                st_n = _usda._normalize_us_state_name(props.get("state_name") or props.get("STATE_NAME") or props.get("state"))
                co_n = _usda._normalize_us_county_name(props.get("county_norm") or props.get("NAME") or props.get("county"), st_n)
                geom_dict = feat.get("geometry")
                if st_n and co_n and geom_dict and shape is not None:
                    geom_obj = shape(geom_dict)
                    if not geom_obj.is_empty:
                        c_pt = geom_obj.representative_point()
                        counties_dict[f"{st_n}|{co_n}"] = {
                            "geom": geom_obj,
                            "centroid_lon": float(props["centroid_lon"]) if props.get("centroid_lon") is not None else float(c_pt.x),
                            "centroid_lat": float(props["centroid_lat"]) if props.get("centroid_lat") is not None else float(c_pt.y),
                        }
        except Exception:
            pass

    counties_nospace: Dict[str, str] = {
        f"{ck.split('|', 1)[0]}|{ck.split('|', 1)[1].replace(' ', '')}": ck
        for ck in counties_dict
        if "|" in ck
    }

    ck_list: List[str] = []
    geoms: List[Any] = []
    for ck, rec in counties_dict.items():
        g = rec.get("geom")
        if g is not None and not g.is_empty:
            ck_list.append(ck)
            geoms.append(g)

    tree = STRtree(geoms) if (STRtree is not None and geoms) else None
    geom_id_map = {id(g): i for i, g in enumerate(geoms)}
    neighbors: Dict[str, Set[str]] = {ck: set() for ck in ck_list}

    if tree is not None:
        for i, g_i in enumerate(geoms):
            ck_i = ck_list[i]
            # 1e-3 deg (~100 m) tolerance bridges simplified polygon vertices across county lines
            buffered_i = g_i.buffer(1e-3)
            for hit in tree.query(buffered_i):
                j = int(hit) if (np is not None and isinstance(hit, (int, np.integer))) or isinstance(hit, int) else geom_id_map.get(id(hit), -1)
                if j >= 0 and j != i and buffered_i.intersects(geoms[j]):
                    ck_j = ck_list[j]
                    neighbors[ck_i].add(ck_j)
                    neighbors[ck_j].add(ck_i)

    graph = {
        "counties": counties_dict,
        "counties_nospace": counties_nospace,
        "ck_list": ck_list,
        "geoms": geoms,
        "tree": tree,
        "geom_id_map": geom_id_map,
        "neighbors": neighbors,
        "is_complete_national_layer": len(counties_dict) >= 500,
    }
    _usda._COUNTY_SPATIAL_GRAPH_CACHE[cache_key] = graph
    return graph


def resolve_voucher_us_county(
    lat: Optional[float],
    lon: Optional[float],
    state_province: Optional[str],
    county: Optional[str],
    spatial_graph: Mapping[str, Any],
) -> Tuple[Optional[str], Optional[str], bool]:
    """
    Resolve `(state_norm, county_norm, is_non_us)` for a voucher using:
      1. Point-in-polygon against `us_counties.geojson` when finite coordinates exist
         (with a 0.02 deg coastal barrier-island snap when state matches or is blank).
      2. Normalized Darwin Core `(stateProvince, county)` fallback.
    """
    import scripts.core.usda_plants_ranges as _usda

    st_dwc = _usda._normalize_us_state_name(state_province)
    is_explicit_non_us = bool(st_dwc and st_dwc not in _usda._US_FULL_STATE_NAMES)

    has_coords = False
    lat_f = math.nan
    lon_f = math.nan
    if lat is not None and lon is not None:
        try:
            lat_f = float(lat)
            lon_f = float(lon)
            if math.isfinite(lat_f) and math.isfinite(lon_f) and -90.0 <= lat_f <= 90.0 and -180.0 <= lon_f <= 180.0:
                has_coords = True
        except (TypeError, ValueError):
            has_coords = False

    tree = (spatial_graph or {}).get("tree")
    geoms = spatial_graph.get("geoms") or []
    ck_list = spatial_graph.get("ck_list") or []
    geom_id_map = spatial_graph.get("geom_id_map") or {}

    if has_coords and tree is not None and Point is not None:
        pt = Point(lon_f, lat_f)
        for hit in tree.query(pt):
            idx = int(hit) if (np is not None and isinstance(hit, (int, np.integer))) or isinstance(hit, int) else geom_id_map.get(id(hit), -1)
            if idx >= 0 and geoms[idx].covers(pt):
                st_hit, co_hit = ck_list[idx].split("|", 1)
                return st_hit, co_hit, (st_hit not in _usda._US_FULL_STATE_NAMES)

        # Coastal / estuarine tolerance (~2 km) strictly for US / unlabeled coastal points (never snap non-US border points into US counties)
        if not is_explicit_non_us:
            pt_buf = pt.buffer(0.02)
            best_idx = -1
            best_dist = float("inf")
            for hit in tree.query(pt_buf):
                idx = int(hit) if (np is not None and isinstance(hit, (int, np.integer))) or isinstance(hit, int) else geom_id_map.get(id(hit), -1)
                if idx >= 0 and geoms[idx].intersects(pt_buf):
                    st_cand = ck_list[idx].split("|", 1)[0]
                    if st_dwc and st_dwc in _usda._US_FULL_STATE_NAMES and st_cand != st_dwc:
                        continue
                    d_val = float(geoms[idx].distance(pt))
                    if d_val < best_dist:
                        best_dist = d_val
                        best_idx = idx
            if best_idx >= 0:
                st_hit, co_hit = ck_list[best_idx].split("|", 1)
                return st_hit, co_hit, (st_hit not in _usda._US_FULL_STATE_NAMES)

    if is_explicit_non_us:
        return st_dwc, _usda._normalize_us_county_name(county, state_norm=st_dwc), True

    if has_coords and st_dwc and st_dwc in _usda._US_FULL_STATE_NAMES:
        in_us_bbox = (17.5 <= lat_f <= 71.5) and (-179.9 <= lon_f <= -64.5)
        in_conus_bbox = (st_dwc in ("Alaska", "Hawaii")) or ((24.3 <= lat_f <= 49.5) and (-125.0 <= lon_f <= -66.8))
        if not (in_us_bbox and in_conus_bbox):
            return None, None, False

    if st_dwc and st_dwc in _usda._US_FULL_STATE_NAMES:
        co_dwc = _usda._normalize_us_county_name(county, state_norm=st_dwc)
        if co_dwc:
            ns_key = f"{st_dwc}|{co_dwc.replace(' ', '')}"
            canon_ck = (spatial_graph.get("counties_nospace") or {}).get(ns_key)
            if canon_ck:
                _, co_dwc = canon_ck.split("|", 1)
        return st_dwc, co_dwc, False

    # Coordinates outside US without US state label (e.g., Canada/Mexico)
    if has_coords and spatial_graph.get("is_complete_national_layer"):
        return None, None, True

    return None, None, False


def _is_within_geodesic_buffer(
    lat_f: Optional[float],
    lon_f: Optional[float],
    voucher_ck: Optional[str],
    doc_counties: Set[str],
    doc_union_geom: Any,
    spatial_graph: Mapping[str, Any],
    buffer_km: float,
) -> bool:
    """Return True if the voucher lies within `buffer_km` of any documented county."""
    import scripts.core.usda_plants_ranges as _usda

    if buffer_km <= 0.0 or not doc_counties:
        return False
    buffer_m = float(buffer_km) * 1000.0
    counties_dict: Mapping[str, Mapping[str, Any]] = spatial_graph.get("counties") or {}

    q_lat = lat_f
    q_lon = lon_f
    if (q_lat is None or q_lon is None) and voucher_ck and voucher_ck in counties_dict:
        c_rec = counties_dict[voucher_ck]
        if c_rec.get("centroid_lat") is not None and c_rec.get("centroid_lon") is not None:
            q_lat = float(c_rec["centroid_lat"])
            q_lon = float(c_rec["centroid_lon"])

    if q_lat is None or q_lon is None:
        return False

    # Exact polygon boundary geodesic check when documented county geometries are available
    if doc_union_geom is not None and Point is not None and nearest_points is not None and not doc_union_geom.is_empty:
        pt = Point(q_lon, q_lat)
        near_geom_pt, _ = nearest_points(doc_union_geom, pt)
        dist_m = _usda.haversine_distance_meters(q_lat, q_lon, float(near_geom_pt.y), float(near_geom_pt.x))
        if dist_m <= buffer_m:
            return True

    # Centroid-to-centroid fallback for counties not in a partial/offline polygon layer
    for d_ck in doc_counties:
        c_rec = counties_dict.get(d_ck)
        if not c_rec:
            continue
        clat = c_rec.get("centroid_lat")
        clon = c_rec.get("centroid_lon")
        if clat is not None and clon is not None:
            if _usda.haversine_distance_meters(q_lat, q_lon, float(clat), float(clon)) <= buffer_m:
                return True

    return False
