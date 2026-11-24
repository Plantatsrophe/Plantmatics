#!/usr/bin/env python3
"""
===============================================================================
Script: download_environmental_layers.py
Project: Plantmatics Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Idempotent, atomic preflight provisioner for fine-grain environmental layers
    used in Step 05 ecological niche modeling and substrate contingency audits:
      1. USDA-NRCS SOLUS100 (100 m CONUS soil & lithic rasters, Nauman et al. 2024):
         - anylithicdpt_cm_p.tif (depth to lithic/paralithic bedrock, cm; scalar=1)
         - fragvol_{0,5,15}_cm_p.tif (% vol; scalar=1)
         - ph1to1h2o_{0,5,15}_cm_p.tif (pH*100; scalar=100)
         - sandtotal_{0,5,15}_cm_p.tif, claytotal_{0,5,15}_cm_p.tif,
           silttotal_{0,5,15}_cm_p.tif (% weight; scalar=1)
         - cec7_{0,5,15}_cm_p.tif (cmol+/kg * 10; scalar=10)
      2. USGS State Geologic Map Compilation (SGMC, Horton et al. 2017, DS-1052):
         - Resolves official bundle from USGS ScienceBase item 5888bf4fe4b05ccb964bab9d
         - Pre-joins SGMC_Geology polygons with SGMC_Units lithology attributes by UNIT_LINK
         - Atomically writes data/environmental/sgmc/sgmc_conus_lithology.gpkg
      3. 1-arc-second (~30 m) Digital Elevation Model tiles (AWS Terrain Skadi / USGS 3DEP):
         - Determines 1x1-degree tiles covering study voucher coordinates
         - Fetches missing 30 m tiles from AWS Open Data Terrain Tiles (skadi .hgt.gz)
           with automatic fallback to USGS 3DEP National Map ImageServer
         - Converts to DEFLATE-compressed GeoTIFFs in data/environmental/dem_3dep_30m/
===============================================================================
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.data_prep.env_sgmc_3dep_utils import (
    MIN_VALID_FILE_BYTES,
    SGMC_SCIENCEBASE_CHILDREN_URL,
    SGMC_SCIENCEBASE_ITEM_ID,
    SGMC_SCIENCEBASE_JSON_URL,
    SGMC_TARGET_GPKG,
    SKADI_BASE_URL,
    USER_AGENT,
    USGS_3DEP_EXPORT_URL,
    _atomic_download as _atomic_download_impl,
    _fetch_json as _fetch_json_impl,
    _is_nonempty_file,
    _is_valid_tiff,
    collect_voucher_dem_tiles,
    ensure_dem_3dep_layers as _ensure_dem_3dep_layers_impl,
    ensure_sgmc_layer as _ensure_sgmc_layer_impl,
    format_dem_tile_id,
    inspect_dem_3dep_cache,
    inspect_sgmc_cache,
    resolve_sgmc_download_candidates,
)

logger = logging.getLogger("packera.env_layers")
if not logger.handlers:
    _handler = logging.StreamHandler(sys.stdout)
    _handler.setFormatter(logging.Formatter("[Env Layers] %(message)s"))
    logger.addHandler(_handler)
    logger.setLevel(logging.INFO)

# Official USDA-NRCS SOLUS100 public bucket (Nauman et al. 2024; doi:10.15482/USDA.ADC/25033856)
SOLUS100_BASE_URL = "https://storage.googleapis.com/solus100pub"
SOLUS100_REQUIRED_LAYERS: Tuple[Dict[str, Any], ...] = (
    {"filename": "anylithicdpt_cm_p.tif", "property": "anylithicdpt", "scalar": 1, "units": "cm", "description": "Depth to lithic or paralithic bedrock contact"},
    {"filename": "fragvol_0_cm_p.tif", "property": "fragvol", "scalar": 1, "units": "percent_volume", "description": "Coarse rock fragments (>2 mm) at 0 cm"},
    {"filename": "fragvol_5_cm_p.tif", "property": "fragvol", "scalar": 1, "units": "percent_volume", "description": "Coarse rock fragments (>2 mm) at 5 cm"},
    {"filename": "fragvol_15_cm_p.tif", "property": "fragvol", "scalar": 1, "units": "percent_volume", "description": "Coarse rock fragments (>2 mm) at 15 cm"},
    {"filename": "ph1to1h2o_0_cm_p.tif", "property": "ph1to1h2o", "scalar": 100, "units": "pH_x100", "description": "Soil pH (1:1 H2O) at 0 cm (divide integer raster by 100)"},
    {"filename": "ph1to1h2o_5_cm_p.tif", "property": "ph1to1h2o", "scalar": 100, "units": "pH_x100", "description": "Soil pH (1:1 H2O) at 5 cm (divide integer raster by 100)"},
    {"filename": "ph1to1h2o_15_cm_p.tif", "property": "ph1to1h2o", "scalar": 100, "units": "pH_x100", "description": "Soil pH (1:1 H2O) at 15 cm (divide integer raster by 100)"},
    {"filename": "sandtotal_0_cm_p.tif", "property": "sandtotal", "scalar": 1, "units": "percent_weight", "description": "Total sand (0.05-2.0 mm, %) at 0 cm"},
    {"filename": "sandtotal_5_cm_p.tif", "property": "sandtotal", "scalar": 1, "units": "percent_weight", "description": "Total sand (0.05-2.0 mm, %) at 5 cm"},
    {"filename": "sandtotal_15_cm_p.tif", "property": "sandtotal", "scalar": 1, "units": "percent_weight", "description": "Total sand (0.05-2.0 mm, %) at 15 cm"},
    {"filename": "claytotal_0_cm_p.tif", "property": "claytotal", "scalar": 1, "units": "percent_weight", "description": "Total clay (<0.002 mm, %) at 0 cm"},
    {"filename": "claytotal_5_cm_p.tif", "property": "claytotal", "scalar": 1, "units": "percent_weight", "description": "Total clay (<0.002 mm, %) at 5 cm"},
    {"filename": "claytotal_15_cm_p.tif", "property": "claytotal", "scalar": 1, "units": "percent_weight", "description": "Total clay (<0.002 mm, %) at 15 cm"},
    {"filename": "silttotal_0_cm_p.tif", "property": "silttotal", "scalar": 1, "units": "percent_weight", "description": "Total silt (0.002-0.05 mm, %) at 0 cm"},
    {"filename": "silttotal_5_cm_p.tif", "property": "silttotal", "scalar": 1, "units": "percent_weight", "description": "Total silt (0.002-0.05 mm, %) at 5 cm"},
    {"filename": "silttotal_15_cm_p.tif", "property": "silttotal", "scalar": 1, "units": "percent_weight", "description": "Total silt (0.002-0.05 mm, %) at 15 cm"},
    {"filename": "cec7_0_cm_p.tif", "property": "cec7", "scalar": 10, "units": "cmol_kg_x10", "description": "Cation exchange capacity (pH 7, cmol+/kg x 10) at 0 cm"},
    {"filename": "cec7_5_cm_p.tif", "property": "cec7", "scalar": 10, "units": "cmol_kg_x10", "description": "Cation exchange capacity (pH 7, cmol+/kg x 10) at 5 cm"},
    {"filename": "cec7_15_cm_p.tif", "property": "cec7", "scalar": 10, "units": "cmol_kg_x10", "description": "Cation exchange capacity (pH 7, cmol+/kg x 10) at 15 cm"},
)


def _atomic_download(
    url: str,
    dest_path: Path,
    timeout: float = 120.0,
    retries: int = 3,
    min_bytes: int = MIN_VALID_FILE_BYTES,
) -> bool:
    return _atomic_download_impl(
        url,
        dest_path,
        timeout=timeout,
        retries=retries,
        min_bytes=min_bytes,
        urlopen_fn=urllib.request.urlopen,
    )


def _fetch_json(url: str, timeout: float = 45.0) -> Optional[Dict[str, Any]]:
    return _fetch_json_impl(url, timeout=timeout, urlopen_fn=urllib.request.urlopen)


def inspect_solus_cache(solus_dir: Path) -> Tuple[bool, List[Dict[str, Any]]]:
    """Return (all_ready, missing_specs) for SOLUS100 rasters in solus_dir."""
    if not solus_dir.is_dir():
        return False, list(SOLUS100_REQUIRED_LAYERS)

    missing = [spec for spec in SOLUS100_REQUIRED_LAYERS if not _is_valid_tiff(solus_dir / str(spec["filename"]))]
    if not missing:
        return True, []

    tifs = [p for p in solus_dir.glob("*.tif*") if _is_valid_tiff(p)]
    has_bed = any(re.search(r"anylithicdpt", p.name, re.I) for p in tifs)
    has_frg = any(re.search(r"(fragvol|cfvo)", p.name, re.I) for p in tifs)
    has_ph = any(re.search(r"^ph(_|1to1h2o)", p.name, re.I) for p in tifs)
    has_tex = (
        any(re.search(r"sandtotal", p.name, re.I) for p in tifs)
        and any(re.search(r"claytotal", p.name, re.I) for p in tifs)
        and any(re.search(r"silttotal", p.name, re.I) for p in tifs)
    )
    has_cec = any(re.search(r"cec7", p.name, re.I) for p in tifs)
    if has_bed and has_frg and has_ph and has_tex and has_cec:
        return True, []
    return False, missing


def ensure_solus100_layers(
    env_dir: Path,
    force: bool = False,
    max_workers: int = 4,
    timeout: float = 180.0,
) -> Dict[str, Any]:
    """Ensure required USDA-NRCS SOLUS100 100 m GeoTIFFs exist in env_dir/solus_100m."""
    solus_dir = env_dir / "solus_100m"
    solus_dir.mkdir(parents=True, exist_ok=True)
    if not force:
        ready, missing_specs = inspect_solus_cache(solus_dir)
        if ready:
            return {"status": "cached", "downloaded": 0, "dir": str(solus_dir)}
    else:
        missing_specs = list(SOLUS100_REQUIRED_LAYERS)

    logger.info(f"Downloading {len(missing_specs)} missing USDA-NRCS SOLUS100 100 m raster(s) into {solus_dir}...")

    def _dl_one(spec: Dict[str, Any]) -> Tuple[str, bool]:
        fname = str(spec["filename"])
        dest = solus_dir / fname
        ok = _atomic_download(f"{SOLUS100_BASE_URL}/{fname}", dest, timeout=timeout, retries=3, min_bytes=MIN_VALID_FILE_BYTES)
        if ok and not _is_valid_tiff(dest):
            dest.unlink(missing_ok=True)
            ok = False
        return fname, ok

    downloaded, failed = 0, []
    workers = max(1, min(int(max_workers), len(missing_specs)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for fname, ok in pool.map(_dl_one, missing_specs):
            if ok:
                downloaded += 1
            else:
                failed.append(fname)

    manifest_path = solus_dir / "solus100_layers_manifest.json"
    manifest_payload = {
        "source": "USDA-NRCS Soil Landscapes of the United States 100-meter (SOLUS100)",
        "doi": "https://doi.org/10.15482/USDA.ADC/25033856",
        "index_url": f"{SOLUS100_BASE_URL}/index.html",
        "layers": list(SOLUS100_REQUIRED_LAYERS),
    }
    tmp_manifest = manifest_path.with_suffix(".json.tmp")
    tmp_manifest.write_text(json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    tmp_manifest.replace(manifest_path)

    status = "downloaded" if not failed else ("partial" if downloaded > 0 else "failed")
    return {"status": status, "downloaded": downloaded, "failed": failed, "dir": str(solus_dir)}


def ensure_sgmc_layer(
    env_dir: Path,
    force: bool = False,
    rscript_bin: Optional[str] = None,
    timeout: float = 300.0,
) -> Dict[str, Any]:
    """Delegate USGS SGMC provisioning with runtime-bound urllib/subprocess callables."""
    return _ensure_sgmc_layer_impl(
        env_dir=env_dir,
        force=force,
        rscript_bin=rscript_bin,
        timeout=timeout,
        urlopen_fn=urllib.request.urlopen,
        run_fn=subprocess.run,
    )


def ensure_dem_3dep_layers(
    env_dir: Path,
    vouchers_csv: Optional[Path] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    force: bool = False,
    max_workers: int = 8,
    rscript_bin: Optional[str] = None,
    timeout: float = 90.0,
) -> Dict[str, Any]:
    """Delegate 3DEP/Skadi 30 m DEM tile provisioning with runtime-bound callables."""
    return _ensure_dem_3dep_layers_impl(
        env_dir=env_dir,
        vouchers_csv=vouchers_csv,
        workspace_root=workspace_root,
        db_path=db_path,
        force=force,
        max_workers=max_workers,
        rscript_bin=rscript_bin,
        timeout=timeout,
        urlopen_fn=urllib.request.urlopen,
        run_fn=subprocess.run,
    )


def ensure_environmental_layers(
    env_dir: Path | str = PROJECT_ROOT / "data" / "environmental",
    vouchers_csv: Optional[Path | str] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    force: bool = False,
    skip: bool = False,
    demo: bool = False,
    max_workers: int = 6,
    rscript_bin: Optional[str] = None,
) -> Dict[str, Any]:
    """Verify SOLUS100, USGS SGMC, and 30 m DEM caches; download only missing layers atomically."""
    if (
        demo
        or skip
        or os.environ.get("PLANTMATICS_DEMO_MODE") == "1"
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or os.environ.get("PLANTMATICS_SKIP_ENV_DOWNLOAD") == "1"
        or os.environ.get("PACKERA_SKIP_ENV_DOWNLOAD") == "1"
    ):
        return {"status": "skipped", "layers": {}}

    env_path = Path(env_dir)
    if not env_path.is_absolute():
        env_path = (workspace_root / env_path).resolve()
    vouchers_path = Path(vouchers_csv) if vouchers_csv is not None else None

    if not force:
        solus_ready, _ = inspect_solus_cache(env_path / "solus_100m")
        sgmc_existing = inspect_sgmc_cache(env_path / "sgmc", env_dir=env_path)
        req_tiles = collect_voucher_dem_tiles(vouchers_csv=vouchers_path, workspace_root=workspace_root, db_path=db_path)
        dem_ready, _ = inspect_dem_3dep_cache(env_path / "dem_3dep_30m", req_tiles)
        if solus_ready and sgmc_existing is not None and dem_ready:
            logger.info(
                f"Verified cached SOLUS100, USGS SGMC, and 3DEP 30 m DEM layers in {env_path} "
                f"( Cohort DEM tiles: {len(req_tiles)}; skipping download )."
            )
            return {
                "status": "cached",
                "layers": {
                    "solus_100m": {"status": "cached"},
                    "sgmc": {"status": "cached", "path": str(sgmc_existing)},
                    "dem_3dep_30m": {"status": "cached", "required_tiles": len(req_tiles)},
                },
            }

    env_path.mkdir(parents=True, exist_ok=True)
    solus_res = ensure_solus100_layers(env_path, force=force, max_workers=min(4, max(1, int(max_workers))))
    sgmc_res = ensure_sgmc_layer(env_path, force=force, rscript_bin=rscript_bin)
    dem_res = ensure_dem_3dep_layers(
        env_path,
        vouchers_csv=vouchers_path,
        workspace_root=workspace_root,
        db_path=db_path,
        force=force,
        max_workers=max(1, int(max_workers)),
        rscript_bin=rscript_bin,
    )

    statuses = {solus_res["status"], sgmc_res["status"], dem_res["status"]}
    overall = (
        "cached"
        if statuses == {"cached"}
        else ("ready" if statuses <= {"cached", "downloaded", "deferred_no_vouchers"} else "partial")
    )
    logger.info(
        f"Summary: SOLUS100={solus_res['status']}, "
        f"USGS_SGMC={sgmc_res['status']}, "
        f"DEM_3DEP_30M={dem_res['status']}."
    )
    return {
        "status": overall,
        "layers": {"solus_100m": solus_res, "sgmc": sgmc_res, "dem_3dep_30m": dem_res},
    }


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="download_environmental_layers.py",
        description="Idempotent downloader for SOLUS100 (100 m), USGS SGMC, and 3DEP/Skadi 30 m DEM layers.",
    )
    parser.add_argument("--env-dir", type=Path, default=PROJECT_ROOT / "data" / "environmental", help="Target environmental root directory (default: data/environmental).")
    parser.add_argument("--vouchers", type=Path, default=None, help="Optional curated_vouchers.csv path used to scope 1x1-degree 30 m DEM tiles.")
    parser.add_argument("--db", type=Path, default=None, help="Optional override path to SQLite database.")
    parser.add_argument("--force", action="store_true", default=False, help="Force re-download even when cached layers already exist.")
    parser.add_argument("--skip", action="store_true", default=False, help="Skip all downloads immediately (no-op).")
    parser.add_argument("--demo", action="store_true", default=False, help="Skip downloads in ephemeral --demo mode.")
    parser.add_argument("-j", "--workers", type=int, default=6, help="Concurrent download worker threads (default: 6).")
    parser.add_argument("--strict", action="store_true", default=False, help="Exit with non-zero status if any missing environmental layer fails to download.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_cli_parser().parse_args(argv)
    summary = ensure_environmental_layers(
        env_dir=args.env_dir,
        vouchers_csv=args.vouchers,
        workspace_root=PROJECT_ROOT,
        db_path=args.db,
        force=args.force,
        skip=args.skip,
        demo=args.demo,
        max_workers=args.workers,
    )
    if args.strict and summary["status"] not in ("cached", "ready", "skipped"):
        return 1
    return 0


__all__ = [
    "MIN_VALID_FILE_BYTES",
    "SGMC_SCIENCEBASE_CHILDREN_URL",
    "SGMC_SCIENCEBASE_ITEM_ID",
    "SGMC_SCIENCEBASE_JSON_URL",
    "SGMC_TARGET_GPKG",
    "SKADI_BASE_URL",
    "SOLUS100_BASE_URL",
    "SOLUS100_REQUIRED_LAYERS",
    "USER_AGENT",
    "USGS_3DEP_EXPORT_URL",
    "build_cli_parser",
    "collect_voucher_dem_tiles",
    "ensure_dem_3dep_layers",
    "ensure_environmental_layers",
    "ensure_sgmc_layer",
    "ensure_solus100_layers",
    "format_dem_tile_id",
    "inspect_dem_3dep_cache",
    "inspect_sgmc_cache",
    "inspect_solus_cache",
    "main",
    "resolve_sgmc_download_candidates",
]


if __name__ == "__main__":
    sys.exit(main())
