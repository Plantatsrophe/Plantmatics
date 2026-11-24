#!/usr/bin/env python3
"""Step 02 3-stage asynchronous GPU/CPU wave engine (prefetch, batched ONNX, parallel finalize, WAL commit)."""

from __future__ import annotations

import concurrent.futures
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import cv2
import numpy as np

from scripts.core.db import record_voucher_segmentation_transaction
from scripts.vision.lm3_data_loader import normalize_pixels_per_mm
from scripts.vision.lm3_inference_engines import _VoucherFinalizeBundle
from scripts.vision.step02_pipe_manifest_flush import (
    prefetch_voucher_data,
    retire_skipped_voucher,
    should_skip_voucher,
)
from scripts.vision.step02_pipe_stage1_candidates import _unpack_reduced_lm3_crop_preds

logger = logging.getLogger("SegmentAndExtract")


def _resolve_wave_worker_counts(pipe: Any) -> Tuple[int, int, int, int, int]:
    """Compute `(prefetch_workers, prep_workers, post_workers, wave_size, prefetch_window)`."""
    if pipe.threads is not None:
        t_val = max(1, int(pipe.threads))
        if t_val == 1:
            prefetch_workers, prep_workers, post_workers = 1, 1, 1
        else:
            prefetch_workers = min(t_val, max(4, min(64, t_val // 2)))
            prep_workers = min(t_val, max(4, min(64, (t_val * 3) // 4)))
            post_workers = min(t_val, 64)
    else:
        eff_cores = min(64, max(4, os.cpu_count() or 32))
        prefetch_workers = min(eff_cores, max(4, min(64, eff_cores // 2)))
        prep_workers = min(eff_cores, max(4, min(64, (eff_cores * 3) // 4)))
        post_workers = eff_cores
    if pipe.prefetch_workers_override is not None:
        prefetch_workers = max(1, int(pipe.prefetch_workers_override))
    if pipe.prep_workers_override is not None:
        prep_workers = max(1, int(pipe.prep_workers_override))
    if pipe.post_workers_override is not None:
        post_workers = max(1, int(pipe.post_workers_override))
    wave_size = max(1, int(pipe.stage1_wave_size))
    return prefetch_workers, prep_workers, post_workers, wave_size, max(wave_size * 3, prefetch_workers * 2)


def _build_wave_ordered_plan(
    pipe: Any, voucher_rows: Sequence[Tuple[Any, Any]], state: Dict[str, Any],
    manifest_path: Path, total_vouchers: int, project_root: Path,
) -> Tuple[List[Tuple[str, Any]], List[Tuple[int, str, Path, str]]]:
    """Build deterministic ordered retirement plan separating completed skips from active work."""
    ordered_plan: List[Tuple[str, Any]] = []
    active_tasks: List[Tuple[int, str, Path, str]] = []
    seen_active_cats: Set[str] = set()
    for _, row in voucher_rows:
        cat_num = str(row["catalogNumber"]).strip()
        raw_img_path = Path(str(row["image_path"]).strip())
        if not raw_img_path.is_absolute():
            raw_img_path = project_root / raw_img_path
        if not raw_img_path.exists():
            for ext in (".webp", ".jpg", ".jpeg", ".png"):
                cand_p = raw_img_path.with_suffix(ext)
                if cand_p.exists():
                    raw_img_path = cand_p
                    break
        if should_skip_voucher(pipe, cat_num, state, manifest_path, total_vouchers) or (
            not pipe.force and cat_num in seen_active_cats
        ):
            ordered_plan.append(("skip", cat_num))
        else:
            a_idx = len(active_tasks)
            taxon = str(row.get("scientificName", row.get("species_raw", ""))).strip()
            active_tasks.append((a_idx, cat_num, raw_img_path, taxon))
            ordered_plan.append(("active", a_idx))
            seen_active_cats.add(cat_num)
    return ordered_plan, active_tasks


def _commit_wave_bundles_to_db(db_conn: Any, wave_bundles: Sequence[_VoucherFinalizeBundle], output_dir: Path) -> None:
    """Persist a completed wave of voucher bundles inside one SQLite WAL transaction with SAVEPOINTs."""
    if db_conn is None or not wave_bundles:
        return
    try:
        with db_conn:
            for idx_b, b_item in enumerate(wave_bundles):
                sp_name = f"sp_v_{idx_b}"
                try:
                    db_conn.execute(f"SAVEPOINT {sp_name};")
                    record_voucher_segmentation_transaction(
                        conn=db_conn, catalog_number=b_item.catalog_number,
                        extracted_records=b_item.extracted_records, failure_record=b_item.failure_record,
                        repro_record=b_item.repro_record, pheno_record=b_item.pheno_record,
                        repro_by_pid=b_item.repro_by_pid, pheno_by_pid=b_item.pheno_by_pid,
                        crowns_by_pid=b_item.crowns_by_pid, rosette_centroids_by_pid=b_item.rosette_centroids_by_pid,
                        by_plant_petiole=b_item.by_plant_petiole, contour_coords_by_leaf=b_item.contour_coords_by_leaf,
                        output_dir=output_dir, image_path=b_item.image_path, commit=False,
                    )
                    db_conn.execute(f"RELEASE SAVEPOINT {sp_name};")
                except Exception as v_db_exc:
                    try:
                        db_conn.execute(f"ROLLBACK TO SAVEPOINT {sp_name};")
                        db_conn.execute(f"RELEASE SAVEPOINT {sp_name};")
                    except Exception:
                        pass
                    logger.warning(f"Per-voucher SQLite transaction failed for {b_item.catalog_number}: {v_db_exc}")
    except Exception as wave_db_exc:
        logger.warning(f"Wave SQLite WAL commit failed: {wave_db_exc}")


def _run_wave_stage_b_gpu(
    pipe: Any,
    wave_prefetched: Sequence[Tuple[int, str, Path, str, Optional[np.ndarray], Optional[np.ndarray]]],
    valid_indices: Sequence[int],
    s1_dets_by_local: Dict[int, List[Any]],
    prep_pool: concurrent.futures.ThreadPoolExecutor,
) -> Tuple[Dict[int, Tuple[Any, ...]], Dict[int, List[Any]], Dict[int, np.ndarray]]:
    """Run Stage B2 crop prep, Stage B2b whole-sheet specimen segmentation, and Stage B3 batched leaf/pose ONNX."""

    def _prep_s2_voucher(
        local_i: int,
    ) -> Tuple[int, Optional[Tuple[Any, ...]], List[Optional[Tuple[np.ndarray, int, int]]], bool]:
        _, cat_v, _, _, img_v, _ = wave_prefetched[local_i]
        dets_v = s1_dets_by_local.get(local_i, [])
        assert img_v is not None
        sh_v, sw_v = img_v.shape[:2]
        sf_v = max(sh_v, sw_v) / 5600.0 if max(sh_v, sw_v) > 1000 else 0.0
        hd_v = max(sh_v, sw_v) > 1000 and sf_v > 0.4
        if not dets_v:
            needs_ss = (not hd_v) or (pipe.stage1_detector.model is None)
            ss_pre = [pipe.model_engine.prepare_transformed_crop(img_v)] if (needs_ss and pipe.model_engine.model is not None) else []
            return local_i, None, ss_pre, needs_ss
        ppm_v: Optional[float] = normalize_pixels_per_mm(pipe.ruler_calibrations.get(cat_v))
        unpacked = pipe._unpack_stage1_candidates(
            catalog_number=cat_v, img_bgr=img_v, stage1_detections=dets_v, pixels_per_mm=ppm_v,
        )
        return local_i, unpacked, [pipe.model_engine.prepare_transformed_crop(c[3]) for c in unpacked[0]], False

    unpacked_by_local: Dict[int, Tuple[Any, ...]] = {}
    precrops_by_local: Dict[int, List[Optional[Tuple[np.ndarray, int, int]]]] = {}
    ss_fallback_locals: List[int] = []
    prep_iter = prep_pool.map(_prep_s2_voucher, valid_indices) if len(valid_indices) >= 2 else (_prep_s2_voucher(loc_i) for loc_i in valid_indices)
    for loc_i, unp, pre_c, needs_ss in prep_iter:
        if unp is not None:
            unpacked_by_local[loc_i], precrops_by_local[loc_i] = unp, pre_c
        elif needs_ss:
            ss_fallback_locals.append(loc_i)
            precrops_by_local[loc_i] = pre_c

    spec_mask_by_local: Dict[int, np.ndarray] = {}
    if getattr(pipe, "specimen_segmenter", None) is not None and getattr(pipe.specimen_segmenter, "model", None) is not None:
        for loc_i in valid_indices:
            img_v = wave_prefetched[loc_i][4]
            if img_v is None or img_v.ndim < 2:
                continue
            sh_v, sw_v = img_v.shape[:2]
            sf_v = max(sh_v, sw_v) / 5600.0 if max(sh_v, sw_v) > 1000 else 0.0
            if max(sh_v, sw_v) > 1000 and sf_v > 0.4:
                try:
                    sm_v = pipe.specimen_segmenter.predict_mask(img_v, return_sheet_size=True)
                except Exception:
                    sm_v = None
                spec_mask_by_local[loc_i] = sm_v if isinstance(sm_v, np.ndarray) else np.empty((0, 0), dtype=np.uint8)

    flat_crops: List[np.ndarray] = []
    flat_pre: List[Optional[Tuple[np.ndarray, int, int]]] = []
    flat_spec_crops: List[Optional[np.ndarray]] = []
    slice_ranges: List[Tuple[int, int, int]] = []
    for loc_i in valid_indices:
        if loc_i not in unpacked_by_local:
            continue
        cands = unpacked_by_local[loc_i][0]
        if not cands:
            slice_ranges.append((loc_i, 0, 0))
            continue
        s_idx = len(flat_crops)
        sm_loc, img_loc = spec_mask_by_local.get(loc_i), wave_prefetched[loc_i][4]
        sh_loc, sw_loc = img_loc.shape[:2] if isinstance(img_loc, np.ndarray) else (0, 0)
        has_sm_loc = isinstance(sm_loc, np.ndarray) and sm_loc.size > 0 and sm_loc.shape[:2] == (sh_loc, sw_loc)
        for c_item in cands:
            flat_crops.append(c_item[3])
            if has_sm_loc and sm_loc is not None:
                cy0_f, cx0_f, cy1_f, cx1_f = c_item[4]
                flat_spec_crops.append(sm_loc[cy0_f:cy1_f, cx0_f:cx1_f])
            else:
                flat_spec_crops.append(None)
        flat_pre.extend(precrops_by_local.get(loc_i, [None] * len(cands)))
        slice_ranges.append((loc_i, s_idx, len(flat_crops)))

    flat_preds = pipe.model_engine.predict_batch(flat_crops, batch_size=pipe.stage2_batch_size, pretransformed=flat_pre) if flat_crops else []
    flat_pose_preds: List[Optional[Any]] = [None] * len(flat_crops)
    if (
        flat_crops and flat_preds and bool(getattr(pipe, "do_apply_landmarks", True))
        and getattr(pipe, "pose_engine", None) is not None and getattr(pipe.pose_engine, "model", None) is not None
        and "predict_landmarks" not in getattr(pipe.pose_engine, "__dict__", {})
        and hasattr(pipe.pose_engine, "predict_batch")
    ):
        pose_indices: List[int] = []
        pose_crops: List[np.ndarray] = []
        pose_lamina_masks: List[Optional[np.ndarray]] = []
        pose_petiole_masks: List[Optional[np.ndarray]] = []
        for idx_f, (c_bgr, c_preds) in enumerate(zip(flat_crops, flat_preds)):
            if not c_preds or not any(int(p_cls) == 0 and np.any(pm_u8 > 0) for _bb, pm_u8, _p_sc, p_cls in c_preds):
                continue
            bm_f, _, pm_f, _ = _unpack_reduced_lm3_crop_preds(
                c_preds, c_bgr.shape[:2],
                specimen_paper_cleaned_mask=(flat_spec_crops[idx_f] if idx_f < len(flat_spec_crops) else None),
                subtract_holes=bool(getattr(pipe, "subtract_holes", True)),
            )
            pose_indices.append(idx_f)
            pose_crops.append(c_bgr)
            pose_lamina_masks.append(bm_f)
            pose_petiole_masks.append(pm_f)
        if pose_crops:
            pose_batch_out = pipe.pose_engine.predict_batch(
                pose_crops, lamina_masks=pose_lamina_masks, petiole_masks=pose_petiole_masks,
                batch_size=max(8, pipe.stage2_batch_size),
            )
            for idx_f, p_res in zip(pose_indices, pose_batch_out):
                flat_pose_preds[idx_f] = p_res

    s2_payload_by_local: Dict[int, Tuple[Any, ...]] = {}
    for loc_i, s_idx, e_idx in slice_ranges:
        unp = unpacked_by_local[loc_i]
        v_preds = list(flat_preds[s_idx:e_idx]) if e_idx > s_idx else []
        v_pose_preds = list(flat_pose_preds[s_idx:e_idx]) if e_idx > s_idx else []
        s2_payload_by_local[loc_i] = (*unp, v_preds, v_pose_preds)

    ss_payload_by_local: Dict[int, List[Any]] = {}
    if ss_fallback_locals:
        ss_imgs = [wave_prefetched[loc_i][4] for loc_i in ss_fallback_locals]
        ss_pre_list = [(precrops_by_local[loc_i][0] if precrops_by_local.get(loc_i) else None) for loc_i in ss_fallback_locals]
        ss_preds_batch = pipe.model_engine.predict_batch(
            ss_imgs, batch_size=pipe.stage2_batch_size, pretransformed=ss_pre_list,  # type: ignore[arg-type]
        )
        for loc_i, ss_raw in zip(ss_fallback_locals, ss_preds_batch):
            ss_payload_by_local[loc_i] = ss_raw
    return s2_payload_by_local, ss_payload_by_local, spec_mask_by_local


def execute_wave_engine(
    pipe: Any,
    voucher_rows: Sequence[Tuple[Any, Any]],
    state: Dict[str, Any],
    manifest_path: Path,
    total_vouchers: int,
    all_extracted_records: List[Dict[str, Any]],
    all_failed_records: List[Dict[str, Any]],
    project_root: Path,
) -> Tuple[int, int]:
    """Execute the 3-stage asynchronous GPU/CPU wave coordinator and return `(processed_count, skipped_count)`."""
    prefetch_workers, prep_workers, post_workers, wave_size, prefetch_window = _resolve_wave_worker_counts(pipe)
    prev_cv2_threads: Optional[int] = None
    prefetch_pool = prep_pool = post_pool = retire_pool = None
    processed_count = skipped_count = 0

    try:
        try:
            prev_cv2_threads = cv2.getNumThreads()
            cv2.setNumThreads(1 if (pipe.threads is not None and int(pipe.threads) == 1) else 2)
        except Exception:
            prev_cv2_threads = None

        ordered_plan, active_tasks = _build_wave_ordered_plan(
            pipe, voucher_rows, state, manifest_path, total_vouchers, project_root
        )
        plan_cursor = 0
        completed_bundles: Dict[int, _VoucherFinalizeBundle] = {}

        def _drain_ordered_plan() -> None:
            nonlocal plan_cursor, processed_count, skipped_count
            while plan_cursor < len(ordered_plan):
                kind, payload = ordered_plan[plan_cursor]
                if kind == "skip":
                    skipped_count = retire_skipped_voucher(
                        pipe, payload, state, all_extracted_records, all_failed_records,
                        processed_count, skipped_count, total_vouchers,
                    )
                    plan_cursor += 1
                    if (processed_count + skipped_count) == total_vouchers and processed_count > 0:
                        logger.info(
                            f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                            f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
                        )
                else:
                    if payload not in completed_bundles:
                        break
                    bundle = completed_bundles.pop(payload)
                    plan_cursor += 1
                    processed_count += 1
                    pipe._apply_finalize_bundle(bundle)
                    if bundle.extracted_records:
                        state["completed_on_disk"].add(bundle.catalog_number)
                    all_extracted_records.extend(bundle.extracted_records)
                    if bundle.failure_record:
                        all_failed_records.append(bundle.failure_record)
                    if processed_count % 25 == 0 or (processed_count + skipped_count) == total_vouchers:
                        logger.info(
                            f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                            f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
                        )

        _drain_ordered_plan()
        if active_tasks:
            logger.info(
                f"3-Stage Wave Engine active: {len(active_tasks)} vouchers queued "
                f"(wave_size={wave_size}, s2_batch={pipe.stage2_batch_size}, "
                f"prefetch_workers={prefetch_workers}, prep_workers={prep_workers}, post_workers={post_workers})."
            )
            prefetch_pool = concurrent.futures.ThreadPoolExecutor(max_workers=prefetch_workers)
            prep_pool = concurrent.futures.ThreadPoolExecutor(max_workers=prep_workers)
            post_pool = concurrent.futures.ThreadPoolExecutor(max_workers=post_workers)
            retire_pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
            pipe.stage1_detector._prep_executor = prep_pool  # type: ignore[attr-defined]
            pipe.model_engine._prep_executor = prep_pool  # type: ignore[attr-defined]
            if getattr(pipe, "pose_engine", None) is not None:
                pipe.pose_engine._prep_executor = prep_pool  # type: ignore[attr-defined]
            prefetch_futures: Dict[int, concurrent.futures.Future] = {}
            retire_futures: List[concurrent.futures.Future] = []
            use_s1_cache = bool(pipe.lm3_project_dir is not None and not pipe.force)

            def _retire_wave_batch(oldest_wave: List[Tuple[int, str, Path, concurrent.futures.Future]]) -> None:
                wave_bundles: List[_VoucherFinalizeBundle] = []
                for a_idx_w, cat_w, path_w, fut_w in oldest_wave:
                    try:
                        bundle = fut_w.result()
                    except Exception as exc_w:
                        logger.warning(f"Unexpected error processing voucher {cat_w} ({path_w}): {exc_w}")
                        bundle = _VoucherFinalizeBundle(
                            catalog_number=cat_w, image_path=str(path_w), extracted_records=[],
                            failure_record={
                                "catalogNumber": cat_w, "image_path": str(path_w),
                                "failure_reason": "SEGMENTATION_RUNTIME_ERROR",
                                "details": f"Runtime exception ({type(exc_w).__name__}) during voucher extraction: {exc_w}",
                            },
                        )
                    completed_bundles[a_idx_w] = bundle
                    wave_bundles.append(bundle)
                _drain_ordered_plan()
                _commit_wave_bundles_to_db(state["db_conn"], wave_bundles, pipe.output_dir)

            n_active = len(active_tasks)
            for w_start in range(0, n_active, wave_size):
                w_end = min(n_active, w_start + wave_size)
                for la_idx in range(w_start, min(n_active, w_end + prefetch_window)):
                    if la_idx not in prefetch_futures:
                        _, la_cat, la_path, _ = active_tasks[la_idx]
                        prefetch_futures[la_idx] = prefetch_pool.submit(prefetch_voucher_data, pipe, la_cat, la_path)

                wave_prefetched: List[Tuple[int, str, Path, str, Optional[np.ndarray], Optional[np.ndarray]]] = []
                valid_indices: List[int] = []
                valid_s1_inputs: List[Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]] = []
                for local_i, (a_idx, cat_num, raw_img_path, taxon) in enumerate(active_tasks[w_start:w_end]):
                    pre_img, pre_s1_lb = prefetch_futures.pop(a_idx).result()
                    wave_prefetched.append((a_idx, cat_num, raw_img_path, taxon, pre_img, pre_s1_lb))
                    if pre_img is not None:
                        valid_indices.append(local_i)
                        valid_s1_inputs.append((pre_img, cat_num, pre_s1_lb))

                s1_dets_by_local: Dict[int, List[Any]] = {}
                if valid_s1_inputs:
                    s1_batch_out = pipe.stage1_detector.predict_batch(
                        valid_s1_inputs, use_cache=use_s1_cache, max_batch_size=wave_size,
                    )
                    for local_i, dets in zip(valid_indices, s1_batch_out):
                        s1_dets_by_local[local_i] = dets

                s2_payload_by_local, ss_payload_by_local, spec_mask_by_local = _run_wave_stage_b_gpu(
                    pipe, wave_prefetched, valid_indices, s1_dets_by_local, prep_pool
                )
                wave_post_tasks: List[Tuple[int, str, Path, concurrent.futures.Future]] = []
                for local_i, (a_idx, cat_num, raw_img_path, taxon, pre_img, _pre_lb) in enumerate(wave_prefetched):
                    fut_post = post_pool.submit(
                        pipe.process_voucher, cat_num, raw_img_path, taxon, pre_img, None,
                        s1_dets_by_local.get(local_i, []) if pre_img is not None else None,
                        s2_payload_by_local.get(local_i), ss_payload_by_local.get(local_i),
                        True, spec_mask_by_local.get(local_i),
                    )
                    wave_post_tasks.append((a_idx, cat_num, raw_img_path, fut_post))

                retire_futures.append(retire_pool.submit(_retire_wave_batch, wave_post_tasks))
                while len(retire_futures) >= 4:
                    retire_futures.pop(0).result()

            while retire_futures:
                retire_futures.pop(0).result()
        _drain_ordered_plan()
    finally:
        pipe.stage1_detector._prep_executor = None  # type: ignore[attr-defined]
        pipe.model_engine._prep_executor = None  # type: ignore[attr-defined]
        if getattr(pipe, "pose_engine", None) is not None:
            pipe.pose_engine._prep_executor = None  # type: ignore[attr-defined]
        for pool_wait in (retire_pool, post_pool):
            if pool_wait is not None:
                pool_wait.shutdown(wait=True)
        for pool_nowait in (prep_pool, prefetch_pool):
            if pool_nowait is not None:
                pool_nowait.shutdown(wait=False, cancel_futures=True)
        if prev_cv2_threads is not None:
            try:
                cv2.setNumThreads(prev_cv2_threads)
            except Exception:
                pass
    return processed_count, skipped_count
