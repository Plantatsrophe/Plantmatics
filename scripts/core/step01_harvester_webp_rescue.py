"""Step 01 WebP sheet migration and secondary duplicate media rescue."""

from __future__ import annotations

import asyncio
import concurrent.futures
import logging
import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Union

from scripts.core.step01_media_qc_cache import (
    audit_voucher_images_parallel,
    rescue_failed_exsiccatae_in_waves,
)

if TYPE_CHECKING:
    import pandas as pd


def migrate_existing_vouchers_to_webp(
    raw_dir: Union[str, Path],
    db_path: Optional[Union[str, Path]] = None,
    curated_csv: Optional[Union[str, Path]] = None,
    quality: int = 92,
    logger: Optional[logging.Logger] = None,
) -> int:
    """
    Scans `raw_dir` for any existing `.jpg`/`.jpeg` voucher images, converts each valid
    image to `.webp` (`quality=92, method=2` by default) via PIL, removes the old `.jpg`
    once the `.webp` is verified non-empty, and updates `.jpg`/`.jpeg` -> `.webp` in
    `curated_csv` and SQLite `vouchers.image_path` if present.
    """
    import scripts.core.harvester as _harv

    raw_dir_path = Path(raw_dir)
    migrated_count = 0
    img_mod = _harv.Image
    pd_mod = _harv.pd

    if raw_dir_path.exists() and raw_dir_path.is_dir() and img_mod is not None:
        jpg_files = sorted(
            p for p in raw_dir_path.iterdir()
            if p.is_file() and p.suffix.lower() in {".jpg", ".jpeg"}
        )
        total_jpgs = len(jpg_files)
        if total_jpgs > 0 and logger:
            logger.info(
                f"Found {total_jpgs} existing JPG voucher image(s) in {raw_dir_path}; "
                f"migrating to WebP (quality={quality})..."
            )

        def _convert_single_jpg(jpg_path: Path) -> int:
            webp_path = jpg_path.with_suffix(".webp")
            if webp_path.is_file() and webp_path.stat().st_size > 0:
                try:
                    jpg_path.unlink()
                    return 1
                except Exception:
                    return 0
            temp_path: Optional[Path] = None
            try:
                with img_mod.open(jpg_path) as img:
                    rgb_img = img.convert("RGB")
                    temp_file = tempfile.NamedTemporaryFile(
                        dir=raw_dir_path,
                        prefix=f"{jpg_path.stem}_",
                        suffix=".webp.tmp",
                        delete=False,
                    )
                    temp_path = Path(temp_file.name)
                    temp_file.close()
                    rgb_img.save(temp_path, format="WEBP", quality=int(quality), method=2)
                if temp_path.is_file() and temp_path.stat().st_size > 0:
                    temp_path.replace(webp_path)
                    if webp_path.is_file() and webp_path.stat().st_size > 0:
                        jpg_path.unlink()
                        return 1
                elif temp_path and temp_path.exists():
                    temp_path.unlink()
            except Exception as exc:
                if temp_path and temp_path.exists():
                    try:
                        temp_path.unlink()
                    except Exception:
                        pass
                if logger:
                    logger.warning(f"Skipping WebP migration for unreadable image {jpg_path.name}: {exc}")
            return 0

        if total_jpgs <= 4:
            for jpg_path in jpg_files:
                migrated_count += _convert_single_jpg(jpg_path)
        else:
            env_workers = (
                os.environ.get("PLANTMATICS_MAX_WORKERS", "").strip()
                or os.environ.get("PACKERA_MAX_WORKERS", "").strip()
            )
            max_workers = (
                max(1, int(env_workers))
                if env_workers.isdigit()
                else min(16, max(4, (os.cpu_count() or 4) * 2))
            )
            completed_jpgs = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
                for res in pool.map(_convert_single_jpg, jpg_files):
                    migrated_count += res
                    completed_jpgs += 1
                    if logger and (completed_jpgs % 250 == 0 or completed_jpgs == total_jpgs):
                        logger.info(
                            f"  [WebP Migration] Processed {completed_jpgs}/{total_jpgs} "
                            f"({100.0 * completed_jpgs / total_jpgs:.1f}%) JPG sheets..."
                        )

    if migrated_count > 0 and logger:
        logger.info(
            f"Migrated {migrated_count} existing JPG voucher image(s) in {raw_dir_path} to WebP (quality={quality})."
        )

    if curated_csv is not None and pd_mod is not None:
        csv_p = Path(curated_csv)
        if csv_p.is_file() and csv_p.stat().st_size > 0:
            try:
                df_csv = pd_mod.read_csv(csv_p)
                if "image_path" in df_csv.columns and not df_csv.empty:
                    old_paths = df_csv["image_path"].astype(str).tolist()
                    new_paths = [
                        re.sub(r"\.jpe?g$", ".webp", str(v), flags=re.IGNORECASE)
                        if not pd_mod.isna(v) and str(v).strip()
                        else v
                        for v in df_csv["image_path"]
                    ]
                    if [str(x) for x in new_paths] != old_paths:
                        df_csv["image_path"] = new_paths
                        tmp_csv = tempfile.NamedTemporaryFile(
                            mode="w",
                            delete=False,
                            dir=csv_p.parent,
                            suffix=".tmp",
                            encoding="utf-8",
                        )
                        tmp_csv_p = Path(tmp_csv.name)
                        try:
                            df_csv.to_csv(tmp_csv, index=False, encoding="utf-8")
                            tmp_csv.flush()
                            os.fsync(tmp_csv.fileno())
                            tmp_csv.close()
                            tmp_csv_p.replace(csv_p)
                        except Exception:
                            if tmp_csv_p.exists():
                                tmp_csv_p.unlink()
            except Exception as csv_exc:
                if logger:
                    logger.warning(f"Could not update image_path in {csv_p}: {csv_exc}")

    try:
        target_csv_ctx = Path(curated_csv) if curated_csv is not None else None
        if db_path is not None or (target_csv_ctx is not None and _harv.is_canonical_or_explicit_db_context(target_csv_ctx)):
            resolved_db = _harv.resolve_db_path(db_path=db_path)
            if resolved_db != ":memory:" and Path(resolved_db).is_file():
                conn = _harv.get_db_connection(db_path=resolved_db, init_schema=False)
                try:
                    conn.execute(
                        "UPDATE vouchers SET image_path = SUBSTR(image_path, 1, LENGTH(image_path) - 4) || '.webp' "
                        "WHERE image_path IS NOT NULL AND LOWER(SUBSTR(image_path, -4)) = '.jpg'"
                    )
                    conn.execute(
                        "UPDATE vouchers SET image_path = SUBSTR(image_path, 1, LENGTH(image_path) - 5) || '.webp' "
                        "WHERE image_path IS NOT NULL AND LOWER(SUBSTR(image_path, -5)) = '.jpeg'"
                    )
                    conn.commit()
                finally:
                    conn.close()
    except Exception as db_exc:
        if logger:
            logger.debug(f"SQLite WebP image_path migration skipped/warning: {db_exc}")

    return migrated_count


def download_and_validate_voucher_media(
    harvester_self: Any,
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Asynchronously downloads voucher specimen sheets, rescuing failed primaries via duplicate fallbacks."""
    import scripts.core.harvester as _harv

    if df.empty or "_image_url" not in df.columns:
        return df, {"success": 0, "skipped": 0, "failed": 0, "quality_rejected": 0}

    download_queue: List[Tuple[str, Path]] = []
    for _, row in df.iterrows():
        url = row["_image_url"]
        dest = harvester_self.workspace_dir / row["image_path"]
        download_queue.append((url, dest))

    harvester_self.logger.info(f"Starting asynchronous download of {len(download_queue)} voucher sheets...")
    stats: Dict[str, Any] = asyncio.run(
        _harv.download_all_voucher_images(
            records_to_download=download_queue,
            concurrency_limit=harvester_self.concurrency,
            force=harvester_self.force,
            webp_quality=harvester_self.webp_quality,
            logger=harvester_self.logger,
        )
    )

    harvester_self.logger.info(
        f"Auditing specimen image quality (min_megapixels={harvester_self.min_megapixels} MP, "
        f"min_file_size_kb={harvester_self.min_file_size_kb} KB)..."
    )
    raw_dir_attr = getattr(harvester_self, "raw_dir", None)
    qc_cache_path = (
        Path(raw_dir_attr).parent / ".voucher_image_qc_cache.json"
        if raw_dir_attr is not None and Path(raw_dir_attr).parent.is_dir()
        else None
    )
    primary_rows = [row for _, row in df.iterrows()]
    primary_paths = [harvester_self.workspace_dir / row["image_path"] for row in primary_rows]
    primary_qc = audit_voucher_images_parallel(
        image_paths=primary_paths,
        min_megapixels=harvester_self.min_megapixels,
        min_file_size_kb=harvester_self.min_file_size_kb,
        check_sharpness=harvester_self.check_sharpness,
        min_sharpness=harvester_self.min_sharpness,
        cache_path=qc_cache_path,
        logger=harvester_self.logger,
    )

    valid_rows: List[pd.Series] = []
    failed_exsiccatae_keys: List[str] = []
    quality_rejected = 0
    mp_values: List[float] = []

    for row, img_dest, (is_valid, q_metrics) in zip(primary_rows, primary_paths, primary_qc):
        if is_valid:
            valid_rows.append(row)
            if "megapixels" in q_metrics:
                mp_values.append(q_metrics["megapixels"])
        else:
            quality_rejected += 1
            ek = str(row.get("exsiccatae_key") or "").strip()
            if ek:
                failed_exsiccatae_keys.append(ek)
            for cand_unlink in (img_dest, _harv.resolve_voucher_image_path(img_dest)):
                if cand_unlink and cand_unlink.is_file():
                    try:
                        cand_unlink.unlink()
                    except Exception:
                        pass

    surviving_cats: Set[str] = {
        str(r.get("catalogNumber", "")).strip() for r in valid_rows if r.get("catalogNumber")
    }
    surviving_eks: Set[str] = {
        str(r.get("exsiccatae_key", "")).strip() for r in valid_rows if r.get("exsiccatae_key")
    }

    # Wave-based concurrent secondary duplicate rescue for failed collection events
    rescue_failed_exsiccatae_in_waves(
        harvester_self=harvester_self,
        failed_exsiccatae_keys=failed_exsiccatae_keys,
        valid_rows=valid_rows,
        surviving_cats=surviving_cats,
        surviving_eks=surviving_eks,
        mp_values=mp_values,
        stats=stats,
        qc_cache_path=qc_cache_path,
    )

    harvester_self.logger.info(
        f"Quality Audit Complete: Retained {len(valid_rows)} vouchers meeting quality standards "
        f"(Rejected {quality_rejected} substandard/low-res images)."
    )
    df_filtered = _harv.pd.DataFrame(valid_rows).reset_index(drop=True) if valid_rows else df.iloc[0:0].copy()
    if not df_filtered.empty:
        df_filtered = _harv.stratify_duplicates(
            df_filtered,
            workspace_dir=harvester_self.workspace_dir,
            logger=harvester_self.logger,
            home_herbarium=harvester_self.home_herbarium,
        )
        df_filtered["spatial_balance_weight"] = 1.0
    elif "spatial_balance_weight" in df.columns:
        df_filtered["spatial_balance_weight"] = _harv.pd.Series(dtype=float)
    stats["quality_rejected"] = quality_rejected
    if mp_values:
        stats["median_mp"] = float(_harv.np.median(mp_values))

    return df_filtered, stats
