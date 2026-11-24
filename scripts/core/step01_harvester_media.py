"""Step 01 voucher image URL resolution, quality auditing, and asynchronous downloading."""

from __future__ import annotations

import asyncio
import logging
import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from scripts.core.config import (
    DEFAULT_MIN_FILE_SIZE_KB,
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN,
)

if TYPE_CHECKING:
    import aiohttp


def optimize_herbarium_image_url(url: str) -> str:
    """Transforms provider-specific URLs to request full-resolution original scans."""
    if not url or not isinstance(url, str):
        return ""

    optimized = url.strip()

    # Smithsonian NMNH: strip dimension clamp (e.g. &h=2000)
    if "collections.nmnh.si.edu/media/" in optimized:
        optimized = re.sub(r"[?&][hw]=\d+", "", optimized)
        if "?" not in optimized and "&" in optimized:
            optimized = optimized.replace("&", "?", 1)

    # Symbiota / SERNEC / SEINet / CCH portals: replace web/thumbnail with orig/large
    if any(k in optimized.lower() for k in ["symbiota", "sernec", "seinet", "cch2", "swbiodiversity"]):
        optimized = re.sub(r"/(?:web|tn|thumbnails?)/", "/orig/", optimized, flags=re.IGNORECASE)
        optimized = re.sub(r"_(?:tn|web|sm)\.(jpe?g|png)", r"_lg.\1", optimized, flags=re.IGNORECASE)

    # IIIF endpoints: replace constrained dimensions with /full/max/
    if "/full/!" in optimized or "/full/pct:" in optimized or re.search(r"/full/\d+,\d*/", optimized):
        optimized = re.sub(r"/full/(?:!?\d+,\d*|pct:\d+)/", "/full/max/", optimized)

    return optimized


def extract_high_res_image_url(media_list: Optional[List[Dict[str, Any]]]) -> Optional[str]:
    """Parses Darwin Core media records, scoring and prioritizing highest-quality specimen image."""
    if not media_list or not isinstance(media_list, list):
        return None

    scored_candidates: List[Tuple[float, str]] = []

    for item in media_list:
        if not isinstance(item, dict):
            continue

        m_type = str(item.get("type", ""))
        m_format = str(item.get("format", "")).lower()
        identifier = str(item.get("identifier", "")).strip()

        if not identifier or not identifier.startswith(("http://", "https://")):
            continue

        score = 0.0
        if m_type == "StillImage" or "image" in m_format:
            score += 50.0
        if identifier.lower().endswith((".jpg", ".jpeg", ".png", ".tif", ".tiff")):
            score += 20.0

        ident_lower = identifier.lower()
        if any(h in ident_lower for h in ["_lg", "_large", "original", "/orig/", "/master/", "/highres/", "/full/", "hires", "high_res"]):
            score += 100.0
        if "max" in ident_lower:
            score += 30.0

        if any(t in ident_lower for t in ["_tn", "_thumb", "thumbnail", "_sm", "_small", "preview", "icon", "mini"]):
            score -= 150.0
        if "detailimages" in ident_lower:
            score -= 25.0
        if re.search(r"[?&]h=(?:[1-9]\d{0,2}|1\d{3}|2000)\b", ident_lower):
            score -= 10.0

        optimized_url = optimize_herbarium_image_url(identifier)
        scored_candidates.append((score, optimized_url))

    if not scored_candidates:
        return None

    scored_candidates.sort(key=lambda x: x[0], reverse=True)
    return scored_candidates[0][1]


def validate_image_quality(
    image_path: Path,
    min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
    min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB,
    check_sharpness: bool = False,
    min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN,
) -> Tuple[bool, Dict[str, Any]]:
    """Evaluates image resolution, byte size, and optical sharpness metrics."""
    import scripts.core.harvester as _harv

    image_path = Path(image_path)
    if not image_path.exists():
        resolved_alt = _harv.resolve_voucher_image_path(image_path)
        if resolved_alt.is_file():
            image_path = resolved_alt
        else:
            return False, {"valid": False, "reason": "file_not_found"}

    file_size_kb = image_path.stat().st_size / 1024.0
    if file_size_kb < min_file_size_kb:
        return False, {
            "valid": False,
            "reason": "file_size_too_small",
            "file_size_kb": round(file_size_kb, 1),
            "min_file_size_kb": min_file_size_kb,
        }

    try:
        with _harv.Image.open(image_path) as img:
            w, h = img.size
            mp = round((w * h) / 1e6, 2)
    except Exception as e:
        return False, {
            "valid": False,
            "reason": f"corrupt_or_unreadable_image: {e}",
            "file_size_kb": round(file_size_kb, 1),
        }

    if mp < min_megapixels:
        return False, {
            "valid": False,
            "reason": "low_resolution",
            "megapixels": mp,
            "min_megapixels": min_megapixels,
            "width": w,
            "height": h,
            "file_size_kb": round(file_size_kb, 1),
        }

    sharpness_score = None
    if check_sharpness:
        try:
            import cv2
            img_cv = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
            if img_cv is not None:
                ch, cw = img_cv.shape[:2]
                y0, y1 = int(ch * 0.15), max(int(ch * 0.85), 1)
                x0, x1 = int(cw * 0.15), max(int(cw * 0.85), 1)
                roi_cv = img_cv[y0:y1, x0:x1] if (y1 > y0 and x1 > x0) else img_cv
                lap_roi = cv2.Laplacian(roi_cv, cv2.CV_64F)
                global_lap_var = float(lap_roi.var())
                rh, rw = lap_roi.shape[:2]
                tile_vars: List[float] = []
                if rh >= 8 and rw >= 8:
                    ys = _harv.np.linspace(0, rh, 9, dtype=int)
                    xs = _harv.np.linspace(0, rw, 9, dtype=int)
                    for ti in range(8):
                        for tj in range(8):
                            tile = lap_roi[ys[ti]:ys[ti + 1], xs[tj]:xs[tj + 1]]
                            if tile.size > 0:
                                tile_vars.append(float(tile.var()))
                tile_p90_lap_var = float(_harv.np.percentile(tile_vars, 90)) if tile_vars else global_lap_var
                sharpness_score = round(max(global_lap_var, tile_p90_lap_var), 2)
                if sharpness_score < min_sharpness:
                    return False, {
                        "valid": False,
                        "reason": "blurry_or_upscaled",
                        "sharpness": sharpness_score,
                        "min_sharpness": min_sharpness,
                        "megapixels": mp,
                        "width": w,
                        "height": h,
                        "file_size_kb": round(file_size_kb, 1),
                    }
        except Exception:
            pass

    return True, {
        "valid": True,
        "megapixels": mp,
        "width": w,
        "height": h,
        "file_size_kb": round(file_size_kb, 1),
        "sharpness": sharpness_score,
    }


async def download_single_image(
    session: aiohttp.ClientSession,
    image_url: str,
    destination_path: Path,
    semaphore: asyncio.Semaphore,
    max_retries: int = 3,
    force: bool = False,
    webp_quality: int = 92,
    logger: Optional[logging.Logger] = None,
) -> Tuple[bool, bool]:
    """Asynchronously downloads a single voucher image file to local storage with retry logic.

    Returns:
        Tuple[bool, bool]: (success, was_skipped)
    """
    import scripts.core.harvester as _harv

    destination_path = Path(destination_path)
    if not force:
        existing_path = (
            destination_path
            if (destination_path.exists() and os.path.getsize(destination_path) > 0)
            else _harv.resolve_voucher_image_path(destination_path)
        )
        if existing_path.is_file() and os.path.getsize(existing_path) > 0:
            if logger:
                logger.info(
                    f"Skipping download for existing voucher {existing_path.name} "
                    f"({os.path.getsize(existing_path)} bytes)."
                )
            return True, True

    async with semaphore:
        for attempt in range(1, max_retries + 1):
            temp_path: Optional[Path] = None
            try:
                if session is None:
                    req = _harv.urllib.request.Request(
                        image_url,
                        headers={"User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium)"},
                    )

                    def _fetch_urllib() -> Tuple[int, str, bytes]:
                        with _harv.urllib.request.urlopen(req, timeout=45) as resp:
                            st = int(getattr(resp, "status", 200))
                            ct = str(resp.headers.get("Content-Type", "")).lower()
                            return st, ct, resp.read()

                    status_code, content_type, content = await asyncio.to_thread(_fetch_urllib)
                else:
                    timeout = _harv.aiohttp.ClientTimeout(total=45, connect=15) if _harv.aiohttp is not None else 45
                    async with session.get(image_url, timeout=timeout) as response:
                        status_code = response.status
                        content_type = response.headers.get("Content-Type", "").lower() if response.status == 200 else ""
                        content = await response.read() if response.status == 200 else b""
                if status_code == 200:
                    if len(content) > 1024 and (
                        not content_type
                        or "image" in content_type
                        or "octet-stream" in content_type
                        or content[:3] == b"\xff\xd8\xff"
                    ):
                        payload_to_write = content
                        if destination_path.suffix.lower() == ".webp" and _harv.Image is not None:
                            try:
                                import io
                                with _harv.Image.open(io.BytesIO(content)) as img:
                                    rgb_img = img.convert("RGB")
                                    buf = io.BytesIO()
                                    rgb_img.save(buf, format="WEBP", quality=int(webp_quality), method=4)
                                    webp_bytes = buf.getvalue()
                                    if len(webp_bytes) > 0:
                                        payload_to_write = webp_bytes
                            except Exception:
                                payload_to_write = content

                        destination_path.parent.mkdir(parents=True, exist_ok=True)
                        temp_file = tempfile.NamedTemporaryFile(
                            dir=destination_path.parent,
                            prefix=f"{destination_path.stem}_",
                            suffix=".tmp",
                            delete=False,
                        )
                        temp_path = Path(temp_file.name)
                        temp_file.write(payload_to_write)
                        temp_file.flush()
                        os.fsync(temp_file.fileno())
                        temp_file.close()
                        temp_path.replace(destination_path)
                        return True, False
                elif status_code in {404, 410}:
                    return False, False
            except (*_harv._AIOHTTP_CLIENT_ERRORS, asyncio.TimeoutError, Exception):
                if temp_path and temp_path.exists():
                    try:
                        temp_path.unlink()
                    except Exception:
                        pass
                if attempt == max_retries:
                    return False, False
                await asyncio.sleep(1.0 * (2 ** (attempt - 1)))
        return False, False


async def download_all_voucher_images(
    records_to_download: List[Tuple[str, Path]],
    concurrency_limit: int = 15,
    force: bool = False,
    webp_quality: int = 92,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, int]:
    """Coordinates asynchronous batch downloading of voucher images with concurrency control."""
    import scripts.core.harvester as _harv

    semaphore = asyncio.Semaphore(concurrency_limit)
    headers = {
        "User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"
    }
    total_records = len(records_to_download)
    stats = {"success": 0, "skipped": 0, "failed": 0, "total": total_records}

    pending = []
    for url, dest in records_to_download:
        dest_p = Path(dest)
        existing_dest = (
            dest_p
            if (not force and dest_p.exists() and os.path.getsize(dest_p) > 0)
            else (_harv.resolve_voucher_image_path(dest_p) if not force else dest_p)
        )
        if not force and existing_dest.is_file() and os.path.getsize(existing_dest) > 0:
            if logger:
                logger.info(
                    f"Retaining existing voucher image: {existing_dest.name} "
                    f"({os.path.getsize(existing_dest)} bytes)."
                )
            stats["skipped"] += 1
        else:
            pending.append((url, dest))

    if not pending:
        if logger:
            logger.info(
                f"All {stats['skipped']} voucher images are already cached locally. "
                f"[Processed: 0 | Skipped: {stats['skipped']} | Total: {total_records}]"
            )
        return stats

    if logger:
        logger.info(
            f"Initiating asynchronous download of {len(pending)} pending images (Concurrency: {concurrency_limit})... "
            f"[Processed: 0 | Skipped: {stats['skipped']} | Total: {total_records}]"
        )

    aiohttp_mod = _harv.aiohttp
    tqdm_mod = _harv.async_tqdm
    if aiohttp_mod is not None:
        connector = aiohttp_mod.TCPConnector(limit=concurrency_limit, limit_per_host=5, ssl=False)
        async with aiohttp_mod.ClientSession(headers=headers, connector=connector) as session:
            tasks = [
                _harv.download_single_image(
                    session, url, dest, semaphore, force=force, webp_quality=webp_quality, logger=logger
                )
                for url, dest in pending
            ]
            results = (
                await tqdm_mod.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
                if tqdm_mod is not None
                else await asyncio.gather(*tasks)
            )
    else:
        tasks = [
            _harv.download_single_image(
                None, url, dest, semaphore, force=force, webp_quality=webp_quality, logger=logger
            )
            for url, dest in pending
        ]
        results = (
            await tqdm_mod.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
            if tqdm_mod is not None
            else await asyncio.gather(*tasks)
        )

    for success, was_skipped in results:
        if success:
            if was_skipped:
                stats["skipped"] += 1
            else:
                stats["success"] += 1
        else:
            stats["failed"] += 1

    if logger:
        logger.info(
            f"Voucher Download Summary: "
            f"[Processed: {stats['success']} | Skipped: {stats['skipped']} | Total: {total_records}]"
        )

    return stats
