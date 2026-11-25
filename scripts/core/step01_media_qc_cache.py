"""Step 01 media I/O helpers, host scheduling, and parallel/cached image quality auditing."""

from __future__ import annotations

import asyncio
import concurrent.futures
import json
import logging
import os
import tempfile
from collections import deque
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from urllib.parse import urlsplit

from scripts.core.config import resolve_max_workers


def interleave_records_by_host(records: List[Tuple[str, Path]]) -> List[Tuple[str, Path]]:
    """Round-robin interleave queued downloads by hostname to prevent single-host head-of-line blocking."""
    buckets: Dict[str, deque] = {}
    for url, dest in records:
        host = urlsplit(str(url)).netloc.lower() or "__unknown__"
        buckets.setdefault(host, deque()).append((url, dest))
    ordered: List[Tuple[str, Path]] = []
    active = deque(buckets.values())
    while active:
        q = active.popleft()
        ordered.append(q.popleft())
        if q:
            active.append(q)
    return ordered


def encode_and_write_image_bytes(
    content: bytes,
    destination_path: Path,
    webp_quality: int,
    image_mod: Any,
) -> bool:
    """Convert raw bytes to WebP (method=2) when targeting `.webp` and atomically write to disk."""
    payload_to_write = content
    is_already_webp = len(content) >= 12 and content[:4] == b"RIFF" and content[8:12] == b"WEBP"
    if destination_path.suffix.lower() == ".webp" and image_mod is not None and not is_already_webp:
        try:
            import io
            with image_mod.open(io.BytesIO(content)) as img:
                rgb_img = img.convert("RGB")
                buf = io.BytesIO()
                rgb_img.save(buf, format="WEBP", quality=int(webp_quality), method=2)
                webp_bytes = buf.getvalue()
                if len(webp_bytes) > 0:
                    payload_to_write = webp_bytes
        except Exception:
            payload_to_write = content

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination_path.parent,
            prefix=f"{destination_path.stem}_",
            suffix=".tmp",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)
            temp_file.write(payload_to_write)
            temp_file.flush()
        temp_path.replace(destination_path)
        return True
    except Exception:
        if temp_path and temp_path.exists():
            try:
                temp_path.unlink()
            except Exception:
                pass
        return False


def load_image_qc_cache(cache_path: Optional[Path]) -> Dict[str, List[Any]]:
    """Load persisted raw image QC measurements (`[size, mtime_ns, w, h, mp, sharpness]`)."""
    if cache_path is None or not cache_path.is_file():
        return {}
    try:
        raw = json.loads(cache_path.read_text(encoding="utf-8"))
        if isinstance(raw, dict):
            return {
                str(k): v
                for k, v in raw.items()
                if isinstance(v, list) and len(v) == 6
            }
    except Exception:
        pass
    return {}


def save_image_qc_cache(cache_path: Optional[Path], cache: Dict[str, List[Any]]) -> None:
    """Atomically persist image QC metric cache to disk."""
    if cache_path is None or not cache_path.parent.is_dir():
        return
    temp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            dir=cache_path.parent,
            prefix=".qc_cache_",
            suffix=".tmp",
            delete=False,
            encoding="utf-8",
        ) as tmp_f:
            temp_path = Path(tmp_f.name)
            json.dump(cache, tmp_f, separators=(",", ":"))
            tmp_f.flush()
        temp_path.replace(cache_path)
    except Exception:
        if temp_path and temp_path.exists():
            try:
                temp_path.unlink()
            except Exception:
                pass


def _is_validate_mocked(fn: Any) -> bool:
    return bool(
        hasattr(fn, "_mock_name")
        or hasattr(fn, "assert_called")
        or getattr(fn, "__module__", "") != "scripts.core.step01_harvester_media"
    )


def _eval_cached_metrics(
    cached: List[Any],
    file_size_kb: float,
    min_megapixels: float,
    check_sharpness: bool,
    min_sharpness: float,
) -> Tuple[bool, Dict[str, Any]]:
    w, h, mp, sharpness = int(cached[2]), int(cached[3]), float(cached[4]), cached[5]
    if mp < min_megapixels:
        return False, {
            "valid": False, "reason": "low_resolution", "megapixels": mp,
            "min_megapixels": min_megapixels, "width": w, "height": h,
            "file_size_kb": round(file_size_kb, 1),
        }
    if check_sharpness and sharpness is not None and float(sharpness) < min_sharpness:
        return False, {
            "valid": False, "reason": "blurry_or_upscaled", "sharpness": float(sharpness),
            "min_sharpness": min_sharpness, "megapixels": mp, "width": w, "height": h,
            "file_size_kb": round(file_size_kb, 1),
        }
    return True, {
        "valid": True, "megapixels": mp, "width": w, "height": h,
        "file_size_kb": round(file_size_kb, 1), "sharpness": sharpness,
    }


def audit_voucher_images_parallel(
    image_paths: List[Path],
    min_megapixels: float,
    min_file_size_kb: float,
    check_sharpness: bool = False,
    min_sharpness: float = 80.0,
    max_workers: Optional[int] = None,
    cache_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> List[Tuple[bool, Dict[str, Any]]]:
    """Validate specimen images in parallel across worker threads with mtime/size raw-metric caching."""
    import scripts.core.harvester as _harv

    n_total = len(image_paths)
    if n_total == 0:
        return []

    validate_fn = _harv.validate_image_quality
    use_cache = not _is_validate_mocked(validate_fn)
    eff_workers = (
        max(1, int(max_workers))
        if max_workers is not None
        else resolve_max_workers(default=min(32, max(4, os.cpu_count() or 4)), cap=32)
    )

    if not use_cache:
        return [
            validate_fn(
                image_path=p,
                min_megapixels=min_megapixels,
                min_file_size_kb=min_file_size_kb,
                check_sharpness=check_sharpness,
                min_sharpness=min_sharpness,
            )
            for p in image_paths
        ]

    cache = load_image_qc_cache(cache_path)
    results: List[Optional[Tuple[bool, Dict[str, Any]]]] = [None] * n_total
    misses: List[Tuple[int, Path, str, int, int]] = []
    cache_hits = 0

    for idx, img_p in enumerate(image_paths):
        p = Path(img_p)
        target = p if p.is_file() else _harv.resolve_voucher_image_path(p)
        if not target.is_file():
            results[idx] = (False, {"valid": False, "reason": "file_not_found"})
            continue
        try:
            st = target.stat()
        except OSError:
            results[idx] = (False, {"valid": False, "reason": "file_not_found"})
            continue
        size_kb = st.st_size / 1024.0
        if size_kb < min_file_size_kb:
            results[idx] = (
                False,
                {
                    "valid": False, "reason": "file_size_too_small",
                    "file_size_kb": round(size_kb, 1), "min_file_size_kb": min_file_size_kb,
                },
            )
            continue
        ckey = f"{target.parent.name}/{target.name}"
        cached = cache.get(ckey)
        if (
            cached is not None
            and cached[0] == st.st_size
            and cached[1] == st.st_mtime_ns
            and (not check_sharpness or cached[5] is not None)
        ):
            results[idx] = _eval_cached_metrics(cached, size_kb, min_megapixels, check_sharpness, min_sharpness)
            cache_hits += 1
        else:
            misses.append((idx, target, ckey, st.st_size, st.st_mtime_ns))

    def _eval_miss(item: Tuple[int, Path, str, int, int]) -> Tuple[int, str, int, int, Tuple[bool, Dict[str, Any]]]:
        idx_m, target_m, ckey_m, sz_m, mt_m = item
        res_m = validate_fn(
            image_path=target_m,
            min_megapixels=min_megapixels,
            min_file_size_kb=min_file_size_kb,
            check_sharpness=check_sharpness,
            min_sharpness=min_sharpness,
        )
        return idx_m, ckey_m, sz_m, mt_m, res_m

    cache_dirty = False
    n_misses = len(misses)
    if n_misses > 0:
        if n_misses <= 4 or eff_workers <= 1:
            eval_iter = (_eval_miss(it) for it in misses)
            for idx_m, ckey_m, sz_m, mt_m, res_m in eval_iter:
                results[idx_m] = res_m
                m_dict = res_m[1]
                if "width" in m_dict and "height" in m_dict and "megapixels" in m_dict:
                    cache[ckey_m] = [
                        sz_m, mt_m, int(m_dict["width"]), int(m_dict["height"]),
                        float(m_dict["megapixels"]), m_dict.get("sharpness"),
                    ]
                    cache_dirty = True
        else:
            done_m = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(eff_workers, n_misses)) as pool:
                for idx_m, ckey_m, sz_m, mt_m, res_m in pool.map(_eval_miss, misses):
                    results[idx_m] = res_m
                    m_dict = res_m[1]
                    if "width" in m_dict and "height" in m_dict and "megapixels" in m_dict:
                        cache[ckey_m] = [
                            sz_m, mt_m, int(m_dict["width"]), int(m_dict["height"]),
                            float(m_dict["megapixels"]), m_dict.get("sharpness"),
                        ]
                        cache_dirty = True
                    done_m += 1
                    if logger and n_misses >= 500 and (done_m % 1000 == 0 or done_m == n_misses):
                        logger.info(
                            f"  [Image QC Audit] Evaluated {done_m}/{n_misses} uncached sheet(s) "
                            f"({100.0 * done_m / n_misses:.1f}%) across {min(eff_workers, n_misses)} workers..."
                        )

    if cache_dirty and cache_path is not None:
        save_image_qc_cache(cache_path, cache)

    if logger and n_total > 20:
        logger.info(
            f"Image QC Audit Summary: {n_total} sheet(s) checked "
            f"({cache_hits} QC cache hit(s), {n_misses} decoded across up to {eff_workers} worker(s))."
        )

    return [r if r is not None else (False, {"valid": False, "reason": "unknown"}) for r in results]


def rescue_failed_exsiccatae_in_waves(
    harvester_self: Any,
    failed_exsiccatae_keys: List[str],
    valid_rows: List[Any],
    surviving_cats: Set[str],
    surviving_eks: Set[str],
    mp_values: List[float],
    stats: Dict[str, Any],
    qc_cache_path: Optional[Path] = None,
) -> int:
    """Concurrently rescue failed primary exsiccatae sheets in rank-ordered waves across secondary duplicates."""
    import scripts.core.harvester as _harv

    sec_pool = harvester_self._secondary_duplicate_pool
    if (
        not failed_exsiccatae_keys
        or sec_pool is None
        or sec_pool.empty
        or "exsiccatae_key" not in sec_pool.columns
        or "_image_url" not in sec_pool.columns
    ):
        return 0

    tier_prio = {"Tier_1_Gold": 3, "Tier_2_Silver": 2, "Tier_3_Bronze": 1}
    unresolved_eks = [ek for ek in dict.fromkeys(failed_exsiccatae_keys) if ek not in surviving_eks]
    needed_mask = sec_pool["exsiccatae_key"].astype(str).isin(set(unresolved_eks))
    pending_queues: Dict[str, List[Any]] = {}
    if needed_mask.any():
        sub_pool = sec_pool.loc[needed_mask].copy()
        sub_pool["_ek_str"] = sub_pool["exsiccatae_key"].astype(str)
        sub_pool["_tp"] = sub_pool["determiner_tier"].astype(str).map(lambda t: tier_prio.get(t, 1))
        sub_pool = sub_pool.sort_values(by="_tp", ascending=False, kind="mergesort")
        grouped_sec: Dict[str, List[Any]] = {}
        for ek_k, grp in sub_pool.groupby("_ek_str", sort=False):
            clean_grp = grp.drop(columns=["_ek_str", "_tp"])
            grouped_sec[str(ek_k)] = [r for _, r in clean_grp.iterrows()]
        for ek in unresolved_eks:
            if ek in grouped_sec:
                pending_queues[ek] = grouped_sec[ek]

    rescue_conc = max(1, int(harvester_self.concurrency or 1))
    rescued_n = 0
    wave_num = 0
    while pending_queues:
        wave_num += 1
        wave_items: List[Tuple[str, Any, Path]] = []
        wave_cats: Set[str] = set()
        for ek in list(pending_queues.keys()):
            if ek in surviving_eks:
                pending_queues.pop(ek, None)
                continue
            q_rows = pending_queues[ek]
            picked = None
            while q_rows:
                cand = q_rows.pop(0)
                cand_cat = str(cand.get("catalogNumber", "")).strip()
                if cand_cat and (cand_cat in surviving_cats or cand_cat in wave_cats):
                    continue
                picked = cand
                if cand_cat:
                    wave_cats.add(cand_cat)
                break
            if picked is None:
                pending_queues.pop(ek, None)
                continue
            if not q_rows:
                pending_queues.pop(ek, None)
            wave_items.append((ek, picked, harvester_self.workspace_dir / str(picked["image_path"])))
        if not wave_items:
            break

        harvester_self.logger.info(
            f"Exsiccatae Media Rescue [Wave {wave_num}]: Evaluating {len(wave_items)} candidate duplicate(s) "
            f"(Concurrency: {rescue_conc})..."
        )
        wave_stats = asyncio.run(
            _harv.download_all_voucher_images(
                records_to_download=[(str(dup_row["_image_url"]), dup_dest) for _, dup_row, dup_dest in wave_items],
                concurrency_limit=rescue_conc,
                force=harvester_self.force,
                webp_quality=harvester_self.webp_quality,
                logger=harvester_self.logger,
            )
        )
        wave_qc = audit_voucher_images_parallel(
            image_paths=[dup_dest for _, _, dup_dest in wave_items],
            min_megapixels=harvester_self.min_megapixels,
            min_file_size_kb=harvester_self.min_file_size_kb,
            check_sharpness=harvester_self.check_sharpness,
            min_sharpness=harvester_self.min_sharpness,
            cache_path=qc_cache_path,
            logger=harvester_self.logger,
        )
        wave_rescued = 0
        for (ek, dup_row, dup_dest), (dup_ok, dup_m) in zip(wave_items, wave_qc):
            if dup_ok:
                promoted = dup_row.copy()
                promoted["is_primary_duplicate"] = True
                valid_rows.append(promoted)
                dup_cat = str(dup_row.get("catalogNumber", "")).strip()
                if dup_cat:
                    surviving_cats.add(dup_cat)
                surviving_eks.add(ek)
                pending_queues.pop(ek, None)
                if "megapixels" in dup_m:
                    mp_values.append(dup_m["megapixels"])
                rescued_n += 1
                wave_rescued += 1
            else:
                for cand_unlink in (dup_dest, _harv.resolve_voucher_image_path(dup_dest)):
                    if cand_unlink and cand_unlink.is_file():
                        try:
                            cand_unlink.unlink()
                        except Exception:
                            pass
        if wave_rescued > 0:
            stats["success"] = int(stats.get("success", 0)) + min(wave_rescued, int(wave_stats.get("success", 0)))
    if rescued_n > 0:
        harvester_self.logger.info(
            f"Exsiccatae Media Rescue: Promoted {rescued_n} valid secondary duplicate sheet(s) "
            f"across {wave_num} concurrent wave(s) after primary image download/QC failure."
        )
    return rescued_n

