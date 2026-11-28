#!/usr/bin/env python3
"""
Fast Step 02 capitulum and peduncle micro-character rescorer (`rescore_step02_capitulum_microchars.py`).

Runs LeafMachine3's Stage 1 `LM3PlantDetector` (`yolo26x_det_1280` only, skipping
UNet++ v9 leaf segmentation and ECT) across existing `PASS` herbarium vouchers in
`data/plantmatics.sqlite`, extracts lateral involucre tip pigmentation, basal wool,
subtending peduncle wool, and calyculus flare per `(catalogNumber, plant_individual_id)`,
and updates `voucher_reproductive_phenology` plus `data/tables/voucher_reproductive_metrics.csv`.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import logging
import math
from pathlib import Path
import sqlite3
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import resolve_cpu_threads, resolve_gpu_batch_size  # noqa: E402
from scripts.core.db import get_db_connection  # noqa: E402
from scripts.vision.lm3_data_loader import LM3_REPRODUCTIVE_CLASS_IDS  # noqa: E402
from scripts.vision.lm3_inference_engines import LM3PlantDetector  # noqa: E402
from scripts.vision.step02_capitulum_microchars import (  # noqa: E402
    MICROCHAR_REPRO_KEYS,
    extract_single_capitulum_microchars,
    summarize_capitulum_microchars,
)
from scripts.vision.step02_photometry_corymb_utils import (  # noqa: E402
    estimate_sheet_paper_cie_l_ref_and_ab,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger("RescoreCapitulumMicrochars")


def _resolve_voucher_image(raw_path: Optional[str], cat: str, data_dir: Path) -> Optional[Path]:
    if raw_path:
        p = Path(raw_path)
        if p.is_file():
            return p
        p_ws = PROJECT_ROOT / raw_path
        if p_ws.is_file():
            return p_ws
    raw_dir = data_dir / "raw_vouchers"
    for ext in (".webp", ".jpg", ".jpeg", ".png"):
        cand = raw_dir / f"{cat}{ext}"
        if cand.is_file():
            return cand
    return None


def _query_target_vouchers(
    conn: sqlite3.Connection,
    taxa: Optional[Sequence[str]],
    missing_only: bool,
    data_dir: Path,
) -> Tuple[List[Tuple[str, Path, Optional[float], bool]], Dict[str, List[Tuple[int, Optional[Tuple[float, float]]]]]]:
    where_parts = ["v.qc_status = 'PASS'"]
    params: List[Any] = []
    if taxa:
        clean_t = [t.strip() for t in taxa if t and t.strip()]
        if clean_t:
            ph = ", ".join("?" for _ in clean_t)
            where_parts.append(f"v.species_standardized IN ({ph})")
            params.extend(clean_t)
    if missing_only:
        where_parts.append(
            "EXISTS (SELECT 1 FROM voucher_reproductive_phenology vrp "
            "WHERE vrp.catalogNumber = v.catalogNumber AND vrp.capitula_count > 0 "
            "AND vrp.phyllary_tip_darkening_delta IS NULL)"
        )
    sql = (
        "SELECT v.catalogNumber, v.image_path, v.pixels_per_mm, v.is_calibrated_mm "
        f"FROM vouchers v WHERE {' AND '.join(where_parts)} ORDER BY v.catalogNumber;"
    )
    rows = conn.execute(sql, params).fetchall()
    vouchers: List[Tuple[str, Path, Optional[float], bool]] = []
    cats: List[str] = []
    for cat_r, img_r, ppm_r, cal_r in rows:
        cat = str(cat_r).strip()
        img_p = _resolve_voucher_image(str(img_r) if img_r else None, cat, data_dir)
        if img_p is None:
            continue
        ppm = float(ppm_r) if (ppm_r is not None and np.isfinite(float(ppm_r)) and float(ppm_r) > 0) else None
        vouchers.append((cat, img_p, ppm, bool(cal_r)))
        cats.append(cat)

    pids_by_cat: Dict[str, List[Tuple[int, Optional[Tuple[float, float]]]]] = defaultdict(list)
    for i in range(0, len(cats), 500):
        chunk = cats[i : i + 500]
        ph = ", ".join("?" for _ in chunk)
        for cat_p, pid_p, crx, cry, rcx, rcy in conn.execute(
            f"SELECT catalogNumber, plant_individual_id, crown_x, crown_y, rosette_centroid_x, rosette_centroid_y "
            f"FROM plant_individuals WHERE catalogNumber IN ({ph}) ORDER BY catalogNumber, plant_individual_id;",
            chunk,
        ).fetchall():
            anchor: Optional[Tuple[float, float]] = None
            if crx is not None and cry is not None:
                anchor = (float(crx), float(cry))
            elif rcx is not None and rcy is not None:
                anchor = (float(rcx), float(rcy))
            pids_by_cat[str(cat_p)].append((int(pid_p or 0), anchor))
    return vouchers, pids_by_cat


def _score_voucher_detections(
    cat: str,
    img_bgr: np.ndarray,
    dets: Sequence[Tuple[Tuple[int, int, int, int], float, int]],
    ppm: Optional[float],
    is_cal: bool,
    pid_anchors: Sequence[Tuple[int, Optional[Tuple[float, float]]]],
) -> List[Tuple[Any, ...]]:
    sh, sw = img_bgr.shape[:2]
    repro_dets = [(bb, sc, int(cid)) for bb, sc, cid in dets if int(cid) in LM3_REPRODUCTIVE_CLASS_IDS]
    pids = list(pid_anchors) if pid_anchors else [(0, None)]
    if not repro_dets:
        return [(*([None] * len(MICROCHAR_REPRO_KEYS)), cat, int(pid)) for pid, _ in pids]

    paper_ref = estimate_sheet_paper_cie_l_ref_and_ab(img_bgr)
    items_by_pid: Dict[int, List[Dict[str, Any]]] = {int(pid): [] for pid, _ in pids}
    all_items: List[Dict[str, Any]] = []
    multi_pid = len(pids) > 1

    for bbox, score, cid in repro_dets:
        ymin, xmin, ymax, xmax = [int(v) for v in bbox]
        y0, y1, x0, x1 = max(0, ymin), min(sh, max(0, ymax)), max(0, xmin), min(sw, max(0, xmax))
        crop_bgr = img_bgr[y0:y1, x0:x1]
        if crop_bgr.size == 0:
            continue
        cx, cy = 0.5 * float(x0 + x1), 0.5 * float(y0 + y1)
        assigned_pid = int(pids[0][0])
        assigned_anchor = pids[0][1]
        if multi_pid:
            best_d = float("inf")
            for pid_cand, anc_cand in pids:
                if anc_cand is None:
                    continue
                d_c = math.hypot(cx - anc_cand[0], cy - anc_cand[1]) * (1.35 if cy > anc_cand[1] + 40.0 else 1.0)
                if d_c < best_d:
                    best_d, assigned_pid, assigned_anchor = d_c, int(pid_cand), anc_cand
        ped_dir = (cx - assigned_anchor[0], cy - assigned_anchor[1]) if assigned_anchor is not None else None
        mc = extract_single_capitulum_microchars(
            crop_bgr=crop_bgr, sheet_bgr=img_bgr, bbox=(y0, x0, y1, x1), mask=None,
            class_id=cid, peduncle_dir_xy=ped_dir, pixels_per_mm=ppm,
            is_calibrated_mm=is_cal, paper_lab_ref=paper_ref,
        )
        item = {"class_id": cid, "score": score, "microchars": mc}
        items_by_pid.setdefault(assigned_pid, []).append(item)
        all_items.append(item)

    updates: List[Tuple[Any, ...]] = []
    seen_pids = set()
    for pid, _ in pids:
        seen_pids.add(int(pid))
        summ = summarize_capitulum_microchars(items_by_pid.get(int(pid), []))
        vals = [float(summ[k]) if np.isfinite(summ[k]) else None for k in MICROCHAR_REPRO_KEYS]
        updates.append((*vals, cat, int(pid)))
    if 0 not in seen_pids and multi_pid:
        summ0 = summarize_capitulum_microchars(all_items)
        vals0 = [float(summ0[k]) if np.isfinite(summ0[k]) else None for k in MICROCHAR_REPRO_KEYS]
        updates.append((*vals0, cat, 0))
    return updates


def rescore_capitulum_microchars(
    db_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    taxa: Optional[Sequence[str]] = None,
    device: str = "cuda",
    threads: Optional[int] = None,
    batch_size: Optional[int] = None,
    wave_size: int = 64,
    missing_only: bool = False,
) -> int:
    eff_data = Path(data_dir) if data_dir is not None else (PROJECT_ROOT / "data")
    n_threads = resolve_cpu_threads(threads)
    eff_bs = resolve_gpu_batch_size(batch_size)
    conn = get_db_connection(db_path=db_path, init_schema=True)
    vouchers, pids_by_cat = _query_target_vouchers(conn, taxa=taxa, missing_only=missing_only, data_dir=eff_data)
    if not vouchers:
        logger.info("No matching PASS vouchers require capitulum micro-character extraction.")
        conn.close()
        return 0

    detector = LM3PlantDetector(device=device, conf_thres=0.25, iou_thres=0.95)
    logger.info(
        f"Rescoring capitulum & peduncle micro-characters on {len(vouchers)} sheets "
        f"(providers={detector.providers}, threads={n_threads}, batch={eff_bs})..."
    )

    set_clause = ", ".join(f"{k} = ?" for k in MICROCHAR_REPRO_KEYS)
    upd_sql = f"UPDATE voucher_reproductive_phenology SET {set_clause} WHERE catalogNumber = ? AND plant_individual_id = ?;"
    total_updated = 0

    def _load_item(v_tup: Tuple[str, Path, Optional[float], bool]) -> Tuple[str, Optional[np.ndarray], Optional[np.ndarray], Optional[float], bool]:
        cat_i, img_p, ppm_i, cal_i = v_tup
        img = cv2.imread(str(img_p))
        if img is None or img.size == 0:
            return (cat_i, None, None, ppm_i, cal_i)
        lb = detector.prepare_letterbox(img)
        return (cat_i, img, lb, ppm_i, cal_i)

    with ThreadPoolExecutor(max_workers=max(2, min(n_threads, 32))) as pool:
        for w0 in range(0, len(vouchers), wave_size):
            wave = vouchers[w0 : w0 + wave_size]
            loaded = list(pool.map(_load_item, wave))
            valid_wave = [(cat_i, img_i, lb_i, ppm_i, cal_i) for cat_i, img_i, lb_i, ppm_i, cal_i in loaded if img_i is not None]
            if not valid_wave:
                continue
            det_batches = detector.predict_batch(
                [(img_i, cat_i, lb_i) for cat_i, img_i, lb_i, _, _ in valid_wave],
                use_cache=False, max_batch_size=eff_bs, executor=pool,
            )
            post_tasks = [
                (cat_i, img_i, dets_i, ppm_i, cal_i, pids_by_cat.get(cat_i, []))
                for (cat_i, img_i, _lb_i, ppm_i, cal_i), dets_i in zip(valid_wave, det_batches)
            ]
            wave_updates_nested = list(pool.map(lambda t: _score_voucher_detections(*t), post_tasks))
            flat_updates = [u for sub in wave_updates_nested for u in sub]
            with conn:
                conn.executemany(upd_sql, flat_updates)
            total_updated += len(valid_wave)
            logger.info(f"Processed {min(w0 + wave_size, len(vouchers))} / {len(vouchers)} sheets...")

    repro_csv = eff_data / "tables" / "voucher_reproductive_metrics.csv"
    if repro_csv.parent.is_dir():
        cur = conn.execute(
            """
            SELECT vrp.catalogNumber, vrp.plant_individual_id, vrp.capitula_count, vrp.corymb_head_count_est,
                   vrp.involucre_height_px, vrp.involucre_width_px, vrp.involucre_height_mm, vrp.involucre_width_mm,
                   vrp.AR_inv, pi.plant_height_cm, pi.scape_to_basal_leaf_ratio, pi.is_bent_stem,
                   vrp.radiate_head_fraction, vrp.phyllary_tip_darkening_delta, vrp.phyllary_tip_anthocyanin_delta,
                   vrp.involucre_basal_cie_l_delta, vrp.peduncle_cie_l_sheet_delta, vrp.calyculus_flare_index
            FROM voucher_reproductive_phenology vrp
            LEFT JOIN plant_individuals pi
              ON pi.catalogNumber = vrp.catalogNumber AND pi.plant_individual_id = vrp.plant_individual_id
            ORDER BY vrp.catalogNumber, vrp.plant_individual_id;
            """
        )
        import csv

        cols = [d[0] for d in cur.description]
        tmp_p = repro_csv.with_suffix(".csv.tmp")
        with open(tmp_p, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(cols)
            writer.writerows(cur.fetchall())
        tmp_p.replace(repro_csv)

    conn.close()
    logger.info(f"Completed capitulum & peduncle micro-character extraction on {total_updated} vouchers.")
    return total_updated


def main() -> None:
    ap = argparse.ArgumentParser(description="Backfill Step 02 capitulum & peduncle micro-characters in SQLite.")
    ap.add_argument("--db-path", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=PROJECT_ROOT / "data")
    ap.add_argument("--taxon", action="append", dest="taxa", default=None, help="Repeatable species filter.")
    ap.add_argument("--species-list", type=str, default="", help="Comma-separated species filter.")
    ap.add_argument("--device", type=str, default="cuda")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--gpu-batch-size", type=int, default=None)
    ap.add_argument("--wave-size", type=int, default=64)
    ap.add_argument("--missing-only", action="store_true", default=False)
    args = ap.parse_args()
    taxa_list = list(args.taxa or [])
    if args.species_list:
        taxa_list.extend([s.strip() for s in args.species_list.split(",") if s.strip()])
    rescore_capitulum_microchars(
        db_path=args.db_path, data_dir=args.data_dir, taxa=taxa_list or None,
        device=args.device, threads=args.threads, batch_size=args.gpu_batch_size,
        wave_size=args.wave_size, missing_only=args.missing_only,
    )


if __name__ == "__main__":
    main()

