"""
Module: step01_cec_ecoregion_utils.py
Description:
    CEC North American Level III ecoregion indexing, Canadian province label
    recognition, point-in-ecoregion lookup, and offline fallback polygon builders.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any, Dict, List, Optional
import urllib.request
import zipfile

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    from shapely.geometry import Point, Polygon, shape
    from shapely.strtree import STRtree
except ImportError:
    Point = Polygon = shape = STRtree = None  # type: ignore[assignment]

CEC_NA_L3_OFFICIAL_ZIP_URL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l3.zip"

_CANADIAN_PROVINCE_TOKENS = {
    "on", "ont", "ontario",
    "qc", "pq", "que", "quebec", "québec",
    "nb", "new brunswick", "nouveau-brunswick",
    "ns", "nova scotia", "nouvelle-écosse", "nouvelle-ecosse",
    "pe", "pei", "prince edward island",
    "nl", "nf", "newfoundland", "labrador", "newfoundland and labrador",
    "mb", "man", "manitoba",
    "sk", "sask", "saskatchewan",
    "ab", "alta", "alberta",
    "bc", "british columbia", "colombie-britannique",
    "yt", "yukon", "yukon territory",
    "nt", "nwt", "northwest territories",
    "nu", "nunavut",
}


def _is_canadian_province_label(state_province: Optional[str]) -> bool:
    if not state_province:
        return False
    cleaned = re.sub(r"[^a-zà-ÿ\s-]", "", str(state_province).strip().lower())
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned in _CANADIAN_PROVINCE_TOKENS


def _build_offline_epa_l3_fallback_index() -> Dict[str, Any]:
    """Deterministic in-memory multi-vertex EPA Level III polygons for isolated offline unit tests (never persisted)."""
    offline_polys = [
        ("Blue_Ridge", Polygon([(-84.5, 34.5), (-83.2, 34.8), (-81.6, 36.1), (-79.0, 37.9), (-78.2, 38.9), (-79.5, 39.2), (-82.8, 36.8), (-84.6, 35.2)])),
        ("Ridge_and_Valley", Polygon([(-85.6, 34.2), (-84.6, 35.2), (-82.8, 36.8), (-79.5, 39.2), (-77.2, 40.8), (-78.8, 41.2), (-83.8, 37.4), (-85.9, 34.8)])),
        ("Piedmont", Polygon([(-85.2, 32.5), (-82.0, 33.5), (-79.4, 34.9), (-77.6, 36.5), (-77.1, 38.9), (-78.2, 38.9), (-79.0, 37.9), (-81.6, 36.1), (-83.2, 34.8), (-84.5, 34.5)])),
        ("Middle_Atlantic_Coastal_Plain", Polygon([(-79.6, 33.0), (-75.4, 35.2), (-75.0, 38.8), (-76.4, 39.2), (-77.1, 38.5), (-77.6, 36.5), (-78.5, 35.2)])),
        ("Southeastern_Plains", Polygon([(-90.2, 31.0), (-85.0, 30.8), (-81.3, 32.2), (-79.6, 33.0), (-78.5, 35.2), (-79.4, 34.9), (-82.0, 33.5), (-85.2, 32.5), (-88.8, 34.8), (-90.1, 33.2)])),
        ("Southern_Coastal_Plain", Polygon([(-85.0, 29.6), (-80.2, 26.5), (-80.8, 31.2), (-81.3, 32.2), (-85.0, 30.8)])),
        ("Mississippi_Alluvial_Plain", Polygon([(-91.8, 29.5), (-90.2, 30.5), (-90.1, 35.0), (-89.2, 37.0), (-90.6, 37.0), (-91.5, 34.2), (-92.0, 31.0)])),
        ("South_Central_Plains", Polygon([(-95.8, 30.2), (-92.0, 31.0), (-91.5, 34.2), (-94.8, 34.1), (-96.0, 32.0)])),
        ("East_Central_Texas_Plains", Polygon([(-98.2, 29.0), (-95.8, 30.2), (-96.0, 32.8), (-97.2, 33.0)])),
        ("Ouachita_Mountains", Polygon([(-95.6, 34.1), (-92.4, 34.3), (-92.6, 35.1), (-95.4, 35.0)])),
        ("Ozark_Highlands", Polygon([(-94.8, 35.7), (-91.0, 35.8), (-90.5, 38.2), (-94.6, 38.1)])),
        ("Interior_Plateau", Polygon([(-88.5, 34.8), (-85.5, 35.2), (-84.5, 38.2), (-87.8, 38.8)])),
        ("Southwestern_Appalachians", Polygon([(-86.8, 33.6), (-85.6, 34.2), (-84.2, 36.8), (-85.5, 36.8)])),
        ("Lake_Erie_Lowland", Polygon([(-83.5, 41.6), (-78.8, 42.7), (-79.2, 44.2), (-81.9, 44.6), (-83.5, 43.2)])),
        ("Eastern_Great_Lakes_Lowlands", Polygon([(-79.2, 43.3), (-75.0, 44.2), (-71.2, 46.8), (-73.8, 47.0), (-76.8, 45.7), (-79.8, 44.6)])),
        ("Acadian_Plains_and_Hills", Polygon([(-67.8, 43.4), (-60.8, 45.2), (-64.2, 47.8), (-67.8, 46.2)])),
        ("Northern_Appalachian_and_Atlantic_Maritime_Highlands", Polygon([(-73.0, 44.5), (-64.0, 44.0), (-59.8, 46.8), (-64.5, 49.2), (-71.2, 46.8)])),
        ("Northern_Lakes_and_Forests", Polygon([(-92.5, 45.0), (-84.0, 45.2), (-78.5, 46.0), (-80.5, 48.5), (-90.0, 48.8)])),
        ("Boreal_Shield", Polygon([(-95.5, 48.8), (-80.5, 48.5), (-55.0, 49.5), (-55.0, 56.5), (-95.5, 56.5)])),
        ("Aspen_Parkland_Northern_Glaciated_Plains", Polygon([(-114.0, 49.0), (-96.0, 49.0), (-96.5, 53.0), (-114.0, 53.5)])),
        ("Boreal_Plain", Polygon([(-120.0, 53.5), (-98.0, 53.0), (-100.0, 58.5), (-120.0, 58.5)])),
        ("Canadian_Rockies", Polygon([(-121.0, 49.0), (-114.0, 49.0), (-115.0, 54.5), (-121.0, 54.5)])),
    ]
    eco_labels = [lbl for lbl, _ in offline_polys]
    eco_geoms = [geom for _, geom in offline_polys]
    offline_counties = {
        "north carolina|orange": {"geom": Polygon([(-79.25, 35.86), (-78.95, 35.86), (-78.95, 36.24), (-79.25, 36.24)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -79.10, "centroid_lat": 36.05},
        "north carolina|new hanover": {"geom": Polygon([(-78.02, 34.05), (-77.78, 34.05), (-77.78, 34.35), (-78.02, 34.35)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -77.90, "centroid_lat": 34.20},
        "north carolina|wake": {"geom": Polygon([(-79.00, 35.50), (-77.80, 35.50), (-77.80, 36.05), (-79.00, 36.05)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -78.40, "centroid_lat": 35.78},
        "texas|nacogdoches": {"geom": Polygon([(-94.90, 31.40), (-94.40, 31.40), (-94.40, 31.85), (-94.90, 31.85)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -94.65, "centroid_lat": 31.62},
        "arkansas|baxter": {"geom": Polygon([(-92.65, 36.05), (-92.10, 36.05), (-92.10, 36.50), (-92.65, 36.50)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -92.38, "centroid_lat": 36.28},
        "arkansas|garland": {"geom": Polygon([(-93.40, 34.42), (-92.90, 34.42), (-92.90, 34.78), (-93.40, 34.78)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -93.15, "centroid_lat": 34.60},
        "ontario|essex": {"geom": Polygon([(-83.12, 41.95), (-82.45, 41.95), (-82.45, 42.32), (-83.12, 42.32)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -82.78, "centroid_lat": 42.14},
        "ontario|ottawa": {"geom": Polygon([(-76.35, 45.10), (-75.35, 45.10), (-75.35, 45.52), (-76.35, 45.52)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -75.70, "centroid_lat": 45.32},
        "quebec|gatineau": {"geom": Polygon([(-76.10, 45.45), (-75.40, 45.45), (-75.40, 45.85), (-76.10, 45.85)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -75.72, "centroid_lat": 45.62},
        "new brunswick|york": {"geom": Polygon([(-67.40, 45.65), (-66.45, 45.65), (-66.45, 46.35), (-67.40, 46.35)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -66.85, "centroid_lat": 45.98},
        "nova scotia|halifax": {"geom": Polygon([(-64.20, 44.45), (-62.80, 44.45), (-62.80, 45.15), (-64.20, 45.15)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -63.55, "centroid_lat": 44.82},
        "manitoba|winnipeg": {"geom": Polygon([(-97.35, 49.72), (-96.95, 49.72), (-96.95, 49.99), (-97.35, 49.99)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -97.14, "centroid_lat": 49.88},
    }
    return {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
        "counties": offline_counties,
        "counties_nospace": {k.replace(" ", ""): v for k, v in offline_counties.items()},
        "use_live_county_intersection": True,
    }


def _build_offline_canadian_cec_l3_index() -> Dict[str, Any]:
    """Deterministic CEC North American Level III polygons for Canadian coordinates when offline."""
    can_polys = [
        ("Lake_Erie_Lowland", Polygon([(-83.5, 41.6), (-78.8, 42.7), (-79.2, 44.2), (-81.9, 44.6), (-83.5, 43.2)])),
        ("Eastern_Great_Lakes_Lowlands", Polygon([(-79.8, 43.3), (-75.0, 44.2), (-71.2, 46.8), (-73.8, 47.0), (-76.8, 45.7), (-79.8, 44.6)])),
        ("Maritime_Lowlands", Polygon([(-65.6, 45.6), (-62.0, 45.8), (-61.8, 47.2), (-65.6, 47.8)])),
        ("Acadian_Plains_and_Hills", Polygon([(-67.8, 43.4), (-60.8, 45.2), (-64.2, 47.8), (-67.8, 46.2)])),
        ("Northern_Appalachian_and_Atlantic_Maritime_Highlands", Polygon([(-73.0, 44.5), (-64.0, 44.0), (-59.8, 46.8), (-64.5, 49.2), (-71.2, 46.8)])),
        ("Northern_Lakes_and_Forests", Polygon([(-92.5, 45.0), (-84.0, 45.2), (-78.5, 46.0), (-80.5, 48.5), (-90.0, 48.8)])),
        ("Boreal_Shield", Polygon([(-95.5, 48.8), (-80.5, 48.5), (-52.5, 46.5), (-52.5, 56.5), (-95.5, 56.5)])),
        ("Aspen_Parkland_Northern_Glaciated_Plains", Polygon([(-114.0, 49.0), (-96.0, 49.0), (-96.5, 53.0), (-114.0, 53.5)])),
        ("Boreal_Plain", Polygon([(-120.0, 53.5), (-98.0, 53.0), (-100.0, 58.5), (-120.0, 58.5)])),
        ("Canadian_Rockies", Polygon([(-121.0, 49.0), (-114.0, 49.0), (-115.0, 54.5), (-121.0, 54.5)])),
    ]
    eco_labels = [lbl for lbl, _ in can_polys]
    eco_geoms = [geom for _, geom in can_polys]
    return {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
    }


def _ensure_cec_na_level3_cache_file(cec_geojson_path: Path) -> bool:
    """Downloads official CEC North American Level III Ecoregions GeoJSON atomically if absent."""
    import scripts.core.epa_ecoregions as _epa

    if cec_geojson_path.exists():
        return True
    if _epa._is_epa_download_mocked():
        return False
    cec_geojson_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_d:
        tmp_dir = Path(tmp_d)
        zip_p = tmp_dir / "na_cec_eco_l3.zip"
        try:
            urllib.request.urlretrieve(_epa.CEC_NA_L3_OFFICIAL_ZIP_URL, str(zip_p))
            with zipfile.ZipFile(zip_p, "r") as zf:
                zf.extractall(tmp_dir)
        except Exception:
            return False
        shp_candidates = list(tmp_dir.rglob("*.shp"))
        if not shp_candidates:
            return False
        shp_path = shp_candidates[0]
        r_code = (
            "suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE);"
            f"eco <- st_read({json.dumps(str(shp_path))}, quiet=TRUE);"
            "nm_col <- intersect(c('NA_L3NAME', 'US_L3NAME', 'L3_KEY'), names(eco))[1];"
            "cd_col <- intersect(c('NA_L3CODE', 'US_L3CODE'), names(eco))[1];"
            "eco$ecoregion_label <- gsub('^_+|_+$', '', gsub('[^A-Za-z0-9]+', '_', trimws(as.character(eco[[nm_col]]))));"
            "eco$NA_L3NAME <- as.character(eco[[nm_col]]);"
            "eco$NA_L3CODE <- if (!is.na(cd_col)) as.character(eco[[cd_col]]) else NA_character_;"
            "eco <- eco[!is.na(eco$ecoregion_label) & nzchar(eco$ecoregion_label) & tolower(eco$ecoregion_label) != 'water', c('NA_L3CODE', 'NA_L3NAME', 'ecoregion_label')];"
            "eco_s <- if (isTRUE(st_is_longlat(eco))) st_simplify(eco, preserveTopology=TRUE, dTolerance=0.002) else st_simplify(eco, preserveTopology=TRUE, dTolerance=200);"
            "eco_w <- st_make_valid(st_transform(eco_s, 4326));"
            f"te <- {json.dumps(str(cec_geojson_path) + '.tmp.geojson')};"
            "st_write(eco_w, te, driver='GeoJSON', quiet=TRUE, delete_dsn=TRUE); file.rename(te, sub('\\.tmp\\.geojson$', '', te));"
        )
        try:
            subprocess.run(["Rscript", "-e", r_code], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            return False
    return cec_geojson_path.exists()


def load_cec_na_level3_spatial_index(allow_download: bool = True) -> Dict[str, Any]:
    """Loads CEC North American Level III spatial index for non-CONUS (e.g. Canadian) coordinates."""
    import scripts.core.epa_ecoregions as _epa

    cec_p = _epa.DEFAULT_CEC_NA_L3_GEOJSON
    if not cec_p.exists() and allow_download and not _epa._is_epa_download_mocked():
        _epa._ensure_cec_na_level3_cache_file(cec_p)
    if not cec_p.exists():
        return _build_offline_canadian_cec_l3_index()
    cec_mtime = cec_p.stat().st_mtime
    cache_key = (str(cec_p.resolve()), cec_mtime)
    if cache_key in _epa._CEC_NA_L3_INDEX_CACHE:
        return _epa._CEC_NA_L3_INDEX_CACHE[cache_key]
    with open(cec_p, "r", encoding="utf-8") as fh:
        cec_fc = json.load(fh)
    eco_geoms: List[Any] = []
    eco_labels: List[str] = []
    for feat in cec_fc.get("features", []):
        props = feat.get("properties") or {}
        raw_lbl = props.get("ecoregion_label") or props.get("NA_L3NAME") or props.get("US_L3NAME") or props.get("name")
        lbl = _epa.standardize_epa_l3_label(raw_lbl)
        geom_dict = feat.get("geometry")
        if lbl and geom_dict:
            geom_obj = shape(geom_dict)
            if not geom_obj.is_empty:
                eco_geoms.append(geom_obj)
                eco_labels.append(lbl)
    if not eco_geoms:
        return _build_offline_canadian_cec_l3_index()
    bundle = {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
    }
    _epa._CEC_NA_L3_INDEX_CACHE[cache_key] = bundle
    return bundle


def _point_in_epa_l3(lon_f: float, lat_f: float, idx_bundle: Dict[str, Any]) -> Optional[str]:
    pt = Point(lon_f, lat_f)
    eco_tree = idx_bundle["eco_tree"]
    eco_geoms = idx_bundle["eco_geoms"]
    eco_labels = idx_bundle["eco_labels"]
    eco_id_map = idx_bundle["eco_id_map"]
    int_types = (int, np.integer) if np is not None else (int,)
    for item in eco_tree.query(pt):
        idx = int(item) if isinstance(item, int_types) else eco_id_map.get(id(item), -1)
        if idx >= 0 and eco_geoms[idx].covers(pt):
            return eco_labels[idx]
    return None
