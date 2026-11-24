#!/usr/bin/env python3
"""USGS SGMC ScienceBase/GPKG and AWS Skadi / USGS 3DEP 30 m DEM helpers."""

from __future__ import annotations

import concurrent.futures
import contextlib
import csv
import gzip
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import subprocess
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
logger = logging.getLogger("packera.env_layers")

USER_AGENT = "PackeraMorphometricsPipeline/1.0 (UNC-Chapel-Hill-Herbarium-NCU)"
MIN_VALID_FILE_BYTES = 1024

SGMC_SCIENCEBASE_ITEM_ID = "5888bf4fe4b05ccb964bab9d"
SGMC_SCIENCEBASE_JSON_URL = f"https://www.sciencebase.gov/catalog/item/{SGMC_SCIENCEBASE_ITEM_ID}?format=json"
SGMC_SCIENCEBASE_CHILDREN_URL = f"https://www.sciencebase.gov/catalog/items?parentId={SGMC_SCIENCEBASE_ITEM_ID}&format=json"
SGMC_TARGET_GPKG = "sgmc_conus_lithology.gpkg"

SKADI_BASE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/skadi"
USGS_3DEP_EXPORT_URL = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage"


def format_dem_tile_id(lat_floor: int, lon_floor: int) -> str:
    """Format integer SW-corner (lat_floor, lon_floor) as SRTM/Skadi tile ID (e.g. N35W079)."""
    return f"{'N' if lat_floor >= 0 else 'S'}{abs(int(lat_floor)):02d}{'E' if lon_floor >= 0 else 'W'}{abs(int(lon_floor)):03d}"


def _is_nonempty_file(path: Path, min_bytes: int = MIN_VALID_FILE_BYTES) -> bool:
    try:
        return path.is_file() and path.stat().st_size >= min_bytes
    except OSError:
        return False


def _is_valid_tiff(path: Path) -> bool:
    if not _is_nonempty_file(path, min_bytes=256):
        return False
    try:
        with path.open("rb") as fh:
            return fh.read(4) in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")
    except OSError:
        return False


def _atomic_download(
    url: str, dest_path: Path, timeout: float = 120.0, retries: int = 3,
    min_bytes: int = MIN_VALID_FILE_BYTES, urlopen_fn: Optional[Callable[..., Any]] = None,
) -> bool:
    """Download URL atomically to dest_path with exponential backoff."""
    opener = urlopen_fn or urllib.request.urlopen
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dest_path.with_name(f".{dest_path.name}.part.{os.getpid()}")
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            with opener(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=timeout) as resp:
                status = getattr(resp, "status", 200)
                if status and int(status) >= 400:
                    raise urllib.error.HTTPError(url, int(status), "HTTP error", resp.headers, None)
                with tmp_path.open("wb") as out_f:
                    shutil.copyfileobj(resp, out_f, length=1024 * 256)
                    out_f.flush()
                    os.fsync(out_f.fileno())
            if not _is_nonempty_file(tmp_path, min_bytes=min_bytes):
                raise IOError(f"Downloaded file from {url} was smaller than {min_bytes} bytes")
            tmp_path.replace(dest_path)
            return True
        except urllib.error.HTTPError as exc:
            last_err = exc
            tmp_path.unlink(missing_ok=True)
            if exc.code == 404:
                return False
        except Exception as exc:
            last_err = exc
            tmp_path.unlink(missing_ok=True)
        if attempt < retries:
            time.sleep(min(6.0, 1.2 * (2 ** (attempt - 1))))
    if last_err is not None:
        logger.warning(f"Download failed for {url}: {last_err}")
    return False


def _fetch_json(url: str, timeout: float = 45.0, urlopen_fn: Optional[Callable[..., Any]] = None) -> Optional[Dict[str, Any]]:
    opener = urlopen_fn or urllib.request.urlopen
    try:
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        with opener(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
        return data if isinstance(data, dict) else None
    except Exception as exc:
        logger.warning(f"Metadata JSON query failed for {url}: {exc}")
        return None


def inspect_sgmc_cache(sgmc_dir: Path, env_dir: Optional[Path] = None) -> Optional[Path]:
    """Return path to an existing valid SGMC vector file if present, else None."""
    for d in [sgmc_dir] + ([env_dir] if env_dir is not None and env_dir != sgmc_dir else []):
        if d.is_dir():
            for p in sorted(d.iterdir()):
                if re.match(r"^sgmc.*\.(gpkg|shp|geojson)$", p.name, flags=re.I) and _is_nonempty_file(p):
                    return p
    return None


def resolve_sgmc_download_candidates(item_json: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Extract candidate (filename, download_url) pairs for SGMC zip bundles from ScienceBase JSON."""
    candidates: List[Tuple[str, str]] = []

    def _scan_files(entries: Any) -> None:
        for f_entry in entries or []:
            if isinstance(f_entry, dict):
                name = str(f_entry.get("name") or "").strip()
                uri = str(f_entry.get("downloadUri") or f_entry.get("url") or "").strip()
                if name.lower().endswith(".zip") and uri.startswith("http"):
                    candidates.append((name, uri))

    _scan_files(item_json.get("files"))
    for facet in item_json.get("facets") or []:
        if isinstance(facet, dict):
            _scan_files(facet.get("files"))
    for wl in item_json.get("webLinks") or []:
        if isinstance(wl, dict):
            uri = str(wl.get("uri") or "").strip()
            title = str(wl.get("title") or os.path.basename(urllib.parse.urlparse(uri).path)).strip()
            if uri.startswith("http") and (uri.lower().endswith(".zip") or "shapefile" in title.lower() or "gdb" in title.lower()):
                candidates.append((title or "sgmc_bundle.zip", uri))

    seen: Set[str] = set()
    uniq = [(n, u) for n, u in candidates if not (u in seen or seen.add(u))]
    uniq.sort(key=lambda it: (0 if ("shapefile" in it[0].lower() or "shp" in it[0].lower()) else (1 if ("geodatabase" in it[0].lower() or "gdb" in it[0].lower() or "sgmc" in it[0].lower()) else 2), it[0].lower()))
    return uniq


def _build_sgmc_join_r_script() -> str:
    """Return R script source that joins SGMC_Geology with SGMC_Units on UNIT_LINK and writes GPKG."""
    return (
        'if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)\n'
        'args <- commandArgs(trailingOnly = TRUE); staging_dir <- args[1]; out_gpkg <- args[2]; suppressPackageStartupMessages({ library(sf) })\n'
        'all_files <- list.files(staging_dir, recursive = TRUE, full.names = TRUE)\n'
        'gdb_dirs <- unique(dirname(all_files[grepl("\\\\.gdb(/|\\\\\\\\)", all_files, ignore.case = TRUE)]))\n'
        'if (length(gdb_dirs) == 0L) { cd <- list.dirs(staging_dir, recursive = TRUE, full.names = TRUE); gdb_dirs <- cd[grepl("\\\\.gdb$", cd, ignore.case = TRUE)] }\n'
        'shp_geol <- all_files[grepl("SGMC_Geology\\\\.shp$", all_files, ignore.case = TRUE)]; geol_sf <- NULL; units_df <- NULL\n'
        'if (length(shp_geol) > 0L) {\n'
        '  geol_sf <- sf::st_read(shp_geol[1L], quiet = TRUE, stringsAsFactors = FALSE); uc <- all_files[grepl("SGMC_Units\\\\.(dbf|csv)$|(^|/|\\\\\\\\)Units\\\\.(dbf|csv)$", all_files, ignore.case = TRUE)]\n'
        '  if (length(uc) > 0L) { up <- uc[1L]; units_df <- if (grepl("\\\\.csv$", up, ignore.case = TRUE)) utils::read.csv(up, stringsAsFactors = FALSE) else if (requireNamespace("foreign", quietly = TRUE)) foreign::read.dbf(up, as.is = TRUE) else as.data.frame(sf::st_read(up, quiet = TRUE, stringsAsFactors = FALSE)) }\n'
        '} else if (length(gdb_dirs) > 0L) {\n'
        '  gp <- gdb_dirs[1L]; lyrs <- sf::st_layers(gp)$name; gl <- lyrs[grepl("SGMC_Geology|Geology", lyrs, ignore.case = TRUE)][1L]; ul <- lyrs[grepl("^SGMC_Units$|^Units$", lyrs, ignore.case = TRUE)][1L]\n'
        '  if (!is.na(gl)) geol_sf <- sf::st_read(gp, layer = gl, quiet = TRUE, stringsAsFactors = FALSE); if (!is.na(ul)) units_df <- as.data.frame(sf::st_read(gp, layer = ul, quiet = TRUE, stringsAsFactors = FALSE))\n'
        '}\n'
        'if (is.null(geol_sf) || nrow(geol_sf) == 0L) stop("Failed to locate non-empty SGMC_Geology polygon layer.")\n'
        'pick_col <- function(df_obj, cands) { if (is.null(df_obj)) return(rep("", nrow(geol_sf))); nup <- toupper(names(df_obj)); for (c in toupper(cands)) { i <- match(c, nup); if (!is.na(i)) { v <- as.character(df_obj[[i]]); v[is.na(v)] <- ""; return(trimws(v)) } }; rep("", nrow(df_obj)) }\n'
        'g_link <- pick_col(geol_sf, c("UNIT_LINK", "UNITLINK")); g_gen <- pick_col(geol_sf, c("GENERALIZED_LITH", "GENERALIZE", "GEN_LITH")); g_orig <- pick_col(geol_sf, c("ORIG_LABEL", "SGMC_LABEL"))\n'
        'maj1 <- pick_col(geol_sf, c("MAJOR1", "ROCKTYPE1")); maj2 <- pick_col(geol_sf, c("MAJOR2", "ROCKTYPE2")); min1 <- pick_col(geol_sf, c("MINOR1")); unm <- pick_col(geol_sf, c("UNIT_NAME", "UNITNAME"))\n'
        'if (!is.null(units_df) && nrow(units_df) > 0L) {\n'
        '  u_link <- pick_col(units_df, c("UNIT_LINK", "UNITLINK")); ku <- nzchar(u_link) & !duplicated(u_link)\n'
        '  if (any(ku)) { uu <- units_df[ku, , drop = FALSE]; mi <- match(g_link, u_link[ku]); hit <- !is.na(mi); if (any(hit)) { um1 <- pick_col(uu, c("MAJOR1", "ROCKTYPE1"))[mi[hit]]; um2 <- pick_col(uu, c("MAJOR2", "ROCKTYPE2"))[mi[hit]]; un1 <- pick_col(uu, c("MINOR1"))[mi[hit]]; unm_u <- pick_col(uu, c("UNIT_NAME", "UNITNAME", "STRAT_UNIT"))[mi[hit]]; maj1[hit] <- ifelse(nzchar(um1), um1, maj1[hit]); maj2[hit] <- ifelse(nzchar(um2), um2, maj2[hit]); min1[hit] <- ifelse(nzchar(un1), un1, min1[hit]); unm[hit] <- ifelse(nzchar(unm_u), unm_u, unm[hit]) } }\n'
        '}\n'
        'out_sf <- sf::st_sf(UNIT_LINK = g_link, ORIG_LABEL = g_orig, UNIT_NAME = unm, MAJOR1 = maj1, MAJOR2 = maj2, ROCKTYPE1 = maj1, ROCKTYPE2 = ifelse(nzchar(maj2), maj2, min1), GENERALIZE = g_gen, geometry = sf::st_geometry(geol_sf))\n'
        'tmp_gpkg <- paste0(out_gpkg, ".tmp.gpkg"); if (file.exists(tmp_gpkg)) unlink(tmp_gpkg); sf::st_write(out_sf, dsn = tmp_gpkg, layer = "sgmc_lithology", driver = "GPKG", quiet = TRUE, append = FALSE)\n'
        'if (!file.exists(tmp_gpkg) || file.info(tmp_gpkg)$size < 1024) stop("Failed to write valid sgmc_conus_lithology.gpkg")\n'
        'if (file.exists(out_gpkg)) unlink(out_gpkg); if (!file.rename(tmp_gpkg, out_gpkg)) { file.copy(tmp_gpkg, out_gpkg, overwrite = TRUE); unlink(tmp_gpkg) }; cat(sprintf("SGMC_GPKG_OK rows=%d\\n", nrow(out_sf)))\n'
    )


def ensure_sgmc_layer(
    env_dir: Path, force: bool = False, rscript_bin: Optional[str] = None, timeout: float = 300.0,
    urlopen_fn: Optional[Callable[..., Any]] = None, run_fn: Optional[Callable[..., Any]] = None,
) -> Dict[str, Any]:
    """Ensure USGS SGMC lithology polygons (pre-joined with SGMC_Units) exist in env_dir/sgmc."""
    runner, sgmc_dir = run_fn or subprocess.run, env_dir / "sgmc"
    sgmc_dir.mkdir(parents=True, exist_ok=True)
    if not force and (existing := inspect_sgmc_cache(sgmc_dir, env_dir=env_dir)) is not None:
        return {"status": "cached", "path": str(existing)}
    if not (r_bin := rscript_bin or shutil.which("Rscript")):
        logger.warning("Rscript not found in PATH; skipping USGS SGMC GeoPackage build.")
        return {"status": "skipped_no_rscript", "dir": str(sgmc_dir)}

    item_json = _fetch_json(SGMC_SCIENCEBASE_JSON_URL, timeout=45.0, urlopen_fn=urlopen_fn)
    candidates: List[Tuple[str, str]] = resolve_sgmc_download_candidates(item_json) if item_json is not None else []
    if not candidates and (children_json := _fetch_json(SGMC_SCIENCEBASE_CHILDREN_URL, timeout=45.0, urlopen_fn=urlopen_fn)) is not None:
        candidates.extend(c for ch in (children_json.get("items") or []) if isinstance(ch, dict) for c in resolve_sgmc_download_candidates(ch))
    if not candidates:
        logger.warning(f"Could not resolve SGMC bundle download URL from ScienceBase ({SGMC_SCIENCEBASE_JSON_URL}).")
        return {"status": "failed_metadata", "dir": str(sgmc_dir)}

    target_gpkg = sgmc_dir / SGMC_TARGET_GPKG
    with tempfile.TemporaryDirectory(prefix="packera_sgmc_") as tmp_dir:
        tmp_root = Path(tmp_dir)
        zip_path, dl_ok, chosen_name = tmp_root / "sgmc_bundle.zip", False, ""
        for cand_name, cand_url in candidates:
            logger.info(f"Downloading USGS SGMC bundle ({cand_name})...")
            if _atomic_download(cand_url, zip_path, timeout=timeout, retries=3, min_bytes=1024 * 64, urlopen_fn=urlopen_fn):
                dl_ok, chosen_name = True, cand_name
                break
        if not dl_ok:
            return {"status": "failed_download", "dir": str(sgmc_dir)}

        extract_dir = tmp_root / "extracted"
        extract_dir.mkdir(parents=True, exist_ok=True)
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extract_dir)
        except Exception as exc:
            logger.warning(f"Failed to extract SGMC zip bundle ({chosen_name}): {exc}")
            return {"status": "failed_unzip", "dir": str(sgmc_dir)}

        r_script_file = tmp_root / "build_sgmc_gpkg.R"
        r_script_file.write_text(_build_sgmc_join_r_script(), encoding="utf-8")
        proc = runner([r_bin, str(r_script_file), str(extract_dir), str(target_gpkg)], capture_output=True, text=True, check=False)
        if proc.returncode != 0 or not _is_nonempty_file(target_gpkg):
            logger.warning(f"USGS SGMC relational join in R failed (exit {proc.returncode}): {proc.stderr.strip()}")
            return {"status": "failed_r_join", "dir": str(sgmc_dir)}

    logger.info(f"Built pre-joined USGS SGMC lithology layer: {target_gpkg}")
    return {"status": "downloaded", "path": str(target_gpkg)}


def collect_voucher_dem_tiles(
    vouchers_csv: Optional[Path] = None, workspace_root: Path = PROJECT_ROOT, db_path: Optional[Path] = None,
) -> Set[Tuple[int, int]]:
    """Collect unique 1x1-degree SW-corner (floor(lat), floor(lon)) cells across study vouchers."""
    tiles: Set[Tuple[int, int]] = set()

    def _add_coord(lat_raw: Any, lon_raw: Any) -> None:
        try:
            lat, lon = float(lat_raw), float(lon_raw)
        except (TypeError, ValueError):
            return
        if math.isfinite(lat) and math.isfinite(lon) and 14.0 <= lat <= 60.0 and -168.0 <= lon <= -50.0:
            tiles.add((int(math.floor(lat)), int(math.floor(lon))))

    canon_csv = workspace_root / "data" / "tables" / "curated_vouchers.csv"
    csv_cands = ([Path(vouchers_csv)] if vouchers_csv is not None else []) + ([canon_csv] if vouchers_csv is None or Path(vouchers_csv) != canon_csv else [])
    for c_path in csv_cands:
        if not c_path.is_file():
            continue
        try:
            with c_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
                reader = csv.DictReader(fh)
                flds = reader.fieldnames or []
                lat_col = next((c for c in ("decimalLatitude", "lat") if c in flds), None)
                lon_col = next((c for c in ("decimalLongitude", "lon") if c in flds), None)
                if lat_col and lon_col:
                    for row in reader:
                        _add_coord(row.get(lat_col), row.get(lon_col))
            if vouchers_csv is not None and c_path == Path(vouchers_csv) and tiles:
                return tiles
        except Exception as exc:
            logger.debug(f"Could not read coordinates from {c_path}: {exc}")

    if not tiles:
        sql_path = db_path or (workspace_root / "data" / "plantmatics.sqlite" if (workspace_root / "data" / "plantmatics.sqlite").is_file() else workspace_root / "data" / "packera_pipeline.sqlite")
        if sql_path.is_file():
            try:
                with contextlib.closing(sqlite3.connect(f"file:{sql_path}?mode=ro", uri=True, timeout=5.0)) as conn:
                    for lat_v, lon_v in conn.execute("SELECT decimalLatitude, decimalLongitude FROM vouchers WHERE decimalLatitude IS NOT NULL AND decimalLongitude IS NOT NULL"):
                        _add_coord(lat_v, lon_v)
            except Exception as exc:
                logger.debug(f"SQLite coordinate lookup skipped: {exc}")
    return tiles


def _load_void_tiles_memo(dem_dir: Path) -> Set[str]:
    memo_path = dem_dir / "ocean_void_tiles.json"
    try:
        data = json.loads(memo_path.read_text(encoding="utf-8")) if memo_path.is_file() else None
        return {str(x) for x in data} if isinstance(data, list) else set()
    except Exception:
        return set()


def _save_void_tiles_memo(dem_dir: Path, void_tiles: Set[str]) -> None:
    memo_path = dem_dir / "ocean_void_tiles.json"
    tmp_p = memo_path.with_suffix(".json.tmp")
    tmp_p.write_text(json.dumps(sorted(void_tiles), indent=2) + "\n", encoding="utf-8")
    tmp_p.replace(memo_path)


def _download_single_3dep_tile(
    lat_floor: int, lon_floor: int, staging_dir: Path, dem_dir: Path,
    timeout: float = 90.0, urlopen_fn: Optional[Callable[..., Any]] = None,
) -> Tuple[str, str, Optional[Path]]:
    """Fetch one 1x1-deg ~30 m DEM tile (Skadi .hgt.gz primary, USGS 3DEP REST fallback)."""
    opener, tile_id = urlopen_fn or urllib.request.urlopen, format_dem_tile_id(lat_floor, lon_floor)
    skadi_url = f"{SKADI_BASE_URL}/{'N' if lat_floor >= 0 else 'S'}{abs(int(lat_floor)):02d}/{tile_id}.hgt.gz"
    gz_path, hgt_path = staging_dir / f"{tile_id}.hgt.gz", staging_dir / f"{tile_id}.hgt"

    if _atomic_download(skadi_url, gz_path, timeout=timeout, retries=2, min_bytes=512, urlopen_fn=urlopen_fn):
        try:
            with gzip.open(gz_path, "rb") as gz_in, hgt_path.open("wb") as hgt_out:
                shutil.copyfileobj(gz_in, hgt_out, length=1024 * 256)
            gz_path.unlink(missing_ok=True)
            if _is_nonempty_file(hgt_path, min_bytes=2_000_000):
                return tile_id, "hgt", hgt_path
            hgt_path.unlink(missing_ok=True)
        except Exception:
            gz_path.unlink(missing_ok=True)
            hgt_path.unlink(missing_ok=True)

    if 24 <= lat_floor <= 50 and -126 <= lon_floor <= -65:
        params = urllib.parse.urlencode({"bbox": f"{lon_floor},{lat_floor},{lon_floor + 1},{lat_floor + 1}", "bboxSR": "4326", "imageSR": "4326", "size": "3600,3600", "format": "tiff", "pixelType": "F32", "noDataInterpretation": "esriNoDataMatchAny", "interpolation": "+RSP_BilinearInterpolation", "f": "image"})
        tif_dest = dem_dir / f"{tile_id}.tif"
        if _atomic_download(f"{USGS_3DEP_EXPORT_URL}?{params}", tif_dest, timeout=timeout, retries=2, min_bytes=4096, urlopen_fn=urlopen_fn):
            if _is_valid_tiff(tif_dest):
                return tile_id, "tif", tif_dest
            tif_dest.unlink(missing_ok=True)

    try:
        with opener(urllib.request.Request(skadi_url, method="HEAD", headers={"User-Agent": USER_AGENT}), timeout=15.0):
            pass
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return tile_id, "void", None
    except Exception:
        pass
    return tile_id, "failed", None


def _convert_hgt_batch_to_geotiff(
    hgt_paths: Sequence[Path], dem_dir: Path, rscript_bin: Optional[str] = None,
    run_fn: Optional[Callable[..., Any]] = None,
) -> int:
    """Convert staged .hgt tiles into DEFLATE-compressed EPSG:4326 GeoTIFFs via R terra (or keep .hgt)."""
    if not hgt_paths:
        return 0
    runner = run_fn or subprocess.run
    if (r_bin := rscript_bin or shutil.which("Rscript")):
        r_conv = (
            'if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)\n'
            'args <- commandArgs(trailingOnly = TRUE); dem_dir <- args[1]; hgt_list <- args[-1]; suppressPackageStartupMessages({ library(terra) }); n_ok <- 0L\n'
            'for (hp in hgt_list) { if (!file.exists(hp)) next; st <- sub("\\\\.hgt$", "", basename(hp), ignore.case = TRUE); ot <- file.path(dem_dir, paste0(st, ".tif")); tt <- file.path(dem_dir, paste0(".", st, ".tmp.tif"))\n'
            '  ok <- tryCatch({ r <- terra::rast(hp); if (!nzchar(terra::crs(r))) terra::crs(r) <- "EPSG:4326"; terra::NAflag(r) <- -32768; terra::writeRaster(r, tt, overwrite = TRUE, datatype = "INT2S", NAflag = -32768, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")); file.exists(tt) && file.info(tt)$size > 1024 }, error = function(e) FALSE)\n'
            '  if (isTRUE(ok)) { if (file.exists(ot)) unlink(ot); if (!file.rename(tt, ot)) { file.copy(tt, ot, overwrite = TRUE); unlink(tt) }; n_ok <- n_ok + 1L } else if (file.exists(tt)) unlink(tt) }\n'
            'cat(sprintf("CONVERTED_TIFS=%d\\n", n_ok))\n'
        )
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as rf:
            rf.write(r_conv)
            r_script_path = Path(rf.name)
        try:
            proc = runner([r_bin, str(r_script_path), str(dem_dir), *[str(p) for p in hgt_paths]], capture_output=True, text=True, check=False)
            if proc.returncode == 0 and (conv := sum(1 for p in hgt_paths if _is_valid_tiff(dem_dir / f"{p.stem}.tif"))) == len(hgt_paths):
                return conv
        finally:
            r_script_path.unlink(missing_ok=True)

    kept = 0
    for p in hgt_paths:
        if _is_valid_tiff(dem_dir / f"{p.stem}.tif"):
            kept += 1
        else:
            shutil.copy2(p, dst_hgt := dem_dir / p.name)
            kept += int(_is_nonempty_file(dst_hgt, min_bytes=2_000_000))
    return kept


def inspect_dem_3dep_cache(dem_dir: Path, required_tiles: Set[Tuple[int, int]]) -> Tuple[bool, List[Tuple[int, int]]]:
    """Return (all_ready, missing_tiles) for 30 m DEM coverage in dem_dir."""
    if not dem_dir.is_dir():
        return (len(required_tiles) == 0), sorted(required_tiles)
    (dem_dir / "dem_3dep_regional.tif").unlink(missing_ok=True)
    existing = [p for p in dem_dir.iterdir() if p.suffix.lower() in (".tif", ".tiff", ".hgt") and p.name.lower() != "dem_3dep_regional.tif" and _is_nonempty_file(p)]
    if not required_tiles:
        return (len(existing) > 0), []
    tile_names = {p.stem.upper() for p in existing if (p.suffix.lower() == ".hgt" or _is_valid_tiff(p))}
    if any(not re.match(r"^[NS]\d{2}[EW]\d{3}$", p.stem, flags=re.I) and _is_valid_tiff(p) for p in existing) and not tile_names:
        return True, []
    voids = _load_void_tiles_memo(dem_dir)
    missing = [(la, lo) for la, lo in sorted(required_tiles) if format_dem_tile_id(la, lo) not in tile_names and format_dem_tile_id(la, lo) not in voids]
    return (len(missing) == 0), missing


def ensure_dem_3dep_layers(
    env_dir: Path, vouchers_csv: Optional[Path] = None, workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None, force: bool = False, max_workers: int = 8,
    rscript_bin: Optional[str] = None, timeout: float = 90.0,
    urlopen_fn: Optional[Callable[..., Any]] = None, run_fn: Optional[Callable[..., Any]] = None,
) -> Dict[str, Any]:
    """Ensure 1-arc-second (~30 m) DEM tiles covering all study voucher coordinates exist."""
    dem_dir = env_dir / "dem_3dep_30m"
    dem_dir.mkdir(parents=True, exist_ok=True)
    required_tiles = collect_voucher_dem_tiles(vouchers_csv=vouchers_csv, workspace_root=workspace_root, db_path=db_path)
    if not required_tiles:
        ready_any, _ = inspect_dem_3dep_cache(dem_dir, set())
        return {"status": "cached" if ready_any else "deferred_no_vouchers", "downloaded": 0, "required_tiles": 0, "dir": str(dem_dir)}
    ready, missing_tiles = (False, sorted(required_tiles)) if force else inspect_dem_3dep_cache(dem_dir, required_tiles)
    if ready:
        return {"status": "cached", "downloaded": 0, "required_tiles": len(required_tiles), "dir": str(dem_dir)}

    logger.info(f"Downloading {len(missing_tiles)} missing 1x1-deg ~30 m DEM tile(s) (out of {len(required_tiles)} cohort cell(s)) into {dem_dir}...")
    void_tiles, staged_hgts, direct_tifs, failed_tiles = _load_void_tiles_memo(dem_dir), [], 0, []
    with tempfile.TemporaryDirectory(prefix="packera_3dep_") as tmp_dir:
        staging_path = Path(tmp_dir)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(int(max_workers), len(missing_tiles)))) as pool:
            for fut in concurrent.futures.as_completed([pool.submit(_download_single_3dep_tile, la, lo, staging_path, dem_dir, timeout, urlopen_fn) for la, lo in missing_tiles]):
                t_id, outcome, p_out = fut.result()
                if outcome == "hgt" and p_out is not None:
                    staged_hgts.append(p_out)
                elif outcome == "tif":
                    direct_tifs += 1
                elif outcome == "void":
                    void_tiles.add(t_id)
                else:
                    failed_tiles.append(t_id)
        converted = _convert_hgt_batch_to_geotiff(staged_hgts, dem_dir=dem_dir, rscript_bin=rscript_bin, run_fn=run_fn)

    if void_tiles:
        _save_void_tiles_memo(dem_dir, void_tiles)
    total_dl = converted + direct_tifs
    status = "downloaded" if not failed_tiles else ("partial" if total_dl > 0 else "failed")
    return {"status": status, "downloaded": total_dl, "required_tiles": len(required_tiles), "failed": sorted(failed_tiles), "dir": str(dem_dir)}
