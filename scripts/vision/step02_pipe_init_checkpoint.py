#!/usr/bin/env python3
"""Step 02 pipeline initialization, LeafMachine3 engine setup, and checkpoint helpers."""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional, Set, Union

import pandas as pd
import yaml

from scripts.core.config import PipelineConfig, resolve_cpu_threads, resolve_gpu_batch_size
from scripts.core.db import is_canonical_or_explicit_db_context
from scripts.vision import lm3_inference_engines as _lm3_engines
from scripts.vision.configure_leafmachine3 import load_native_lm3_config
from scripts.vision.lm3_data_loader import RulerCalibrationMap, load_ruler_calibrations
from scripts.vision.lm3_inference_engines import _VoucherFinalizeBundle, is_stock_lm3_weights

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _resolve_env_int(*env_names: str) -> Optional[int]:
    """Return the first non-empty positive integer from the given environment variable names."""
    for name in env_names:
        val = os.environ.get(name, "").strip()
        if val:
            return max(1, int(val))
    return None


def _load_merged_lm3_settings(project_root: Path) -> Dict[str, Any]:
    """Load native LeafMachine3 config overlaid with `config/LM3_settings.yaml`."""
    try:
        lm3_settings: Dict[str, Any] = load_native_lm3_config(project_root=project_root)
    except Exception:
        lm3_settings = {}
    cfg_path = project_root / "config" / "LM3_settings.yaml"
    if not cfg_path.is_file():
        return lm3_settings
    try:
        loaded = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            for k_top, v_top in loaded.items():
                if isinstance(v_top, dict) and isinstance(lm3_settings.get(k_top), dict):
                    merged_sub = dict(lm3_settings[k_top])
                    for k_sub, v_sub in v_top.items():
                        if isinstance(v_sub, dict) and isinstance(merged_sub.get(k_sub), dict):
                            merged_sub[k_sub] = {**merged_sub[k_sub], **v_sub}
                        else:
                            merged_sub[k_sub] = v_sub
                    lm3_settings[k_top] = merged_sub
                else:
                    lm3_settings[k_top] = v_top
    except Exception:
        pass
    return lm3_settings


def _get_lm3_mod_section(
    lm3_settings: Dict[str, Any], name: str, legacy_alias: Optional[str] = None
) -> Dict[str, Any]:
    """Resolve a module configuration dict from `modules.<name>` or top-level keys."""
    mods_cfg = lm3_settings.get("modules") if isinstance(lm3_settings.get("modules"), dict) else {}
    for container in (mods_cfg, lm3_settings):
        sec = container.get(name)
        if isinstance(sec, dict):
            return sec
        if legacy_alias and isinstance(container.get(legacy_alias), dict):
            return container[legacy_alias]
    return {}


def init_segment_pipeline_state(
    pipe: Any,
    engine_ns: Optional[Dict[str, Any]],
    vouchers_csv: Union[str, Path],
    model_weights: Union[str, Path],
    output_dir: Union[str, Path],
    device: str = "cuda",
    score_thresh: float = 0.35,
    force: bool = False,
    config: Optional[PipelineConfig] = None,
    lm3_project_dir: Optional[Union[str, Path]] = None,
    stage1_weights: Optional[Union[str, Path]] = None,
    ruler_confidence_thresh: float = 0.50,
    db_path: Optional[Union[str, Path]] = None,
    threads: Optional[int] = None,
    gpu_batch_size: Optional[int] = None,
    wave_size: Optional[int] = None,
    prefetch_workers: Optional[int] = None,
    prep_workers: Optional[int] = None,
    post_workers: Optional[int] = None,
    save_debug_crops: Optional[bool] = None,
) -> None:
    """Initialize all attributes, directories, LeafMachine3 ONNX engines, and worker overrides on `pipe`."""
    ns = engine_ns or {}
    plant_det_cls = ns.get("LM3PlantDetector", _lm3_engines.LM3PlantDetector)
    leaf_seg_cls = ns.get("LM3LeafSegmenterEngine", _lm3_engines.LM3LeafSegmenterEngine)
    pose_cls = ns.get("LM3PoseLandmarkEngine", _lm3_engines.LM3PoseLandmarkEngine)
    scale_cls = ns.get("LM3ScaleCalibrator", _lm3_engines.LM3ScaleCalibrator)
    spec_cls = ns.get("LM3SpecimenSegmenter", _lm3_engines.LM3SpecimenSegmenter)

    pipe.db_path = db_path
    pipe.vouchers_csv = Path(vouchers_csv)
    pipe.model_weights = Path(model_weights)
    pipe.output_dir = Path(output_dir).resolve()
    if pipe.output_dir.name != "data":
        pipe.output_dir = pipe.output_dir / "data"
    pipe.device = device
    pipe.config = config or PipelineConfig.from_yaml()
    pipe.score_thresh = score_thresh
    pipe.threads = resolve_cpu_threads(threads, pipe.config)
    if pipe.threads is not None and pipe.device == "cpu":
        try:
            import torch as _torch_cpu

            _torch_cpu.set_num_threads(max(1, int(pipe.threads)))
        except Exception:
            pass
    pipe.gpu_batch_size = resolve_gpu_batch_size(
        gpu_batch_size, pipe.config, default_cuda=10, default_cpu=8, use_cuda=(pipe.device == "cuda"),
    )
    has_explicit_gpu_bs = (
        gpu_batch_size is not None
        or bool(os.environ.get("PLANTMATICS_GPU_BATCH_SIZE", "").strip() or os.environ.get("PACKERA_GPU_BATCH_SIZE", "").strip())
        or pipe.config.compute.gpu_batch_size is not None
    )
    pipe.force = force
    pipe.ruler_confidence_thresh = ruler_confidence_thresh
    pipe.lm3_project_dir = Path(lm3_project_dir) if lm3_project_dir is not None else None

    pipe.masks_dir = pipe.output_dir / "masks"
    pipe.contours_dir = pipe.output_dir / "contours"
    pipe.ect_dir = pipe.output_dir / "ect"
    pipe.tables_dir = pipe.output_dir / "tables"
    pipe.cropped_patches_dir = pipe.output_dir / "cropped_patches"
    pipe.capitula_patches_dir = pipe.cropped_patches_dir / "capitula"
    for d in (pipe.masks_dir, pipe.contours_dir, pipe.ect_dir, pipe.tables_dir):
        d.mkdir(parents=True, exist_ok=True)

    if not pipe.vouchers_csv.exists() and pipe.vouchers_csv.name == "curated_vouchers.csv":
        shared_vouchers = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
        if shared_vouchers.exists() and shared_vouchers.resolve() != pipe.vouchers_csv.resolve():
            pipe.vouchers_csv.parent.mkdir(parents=True, exist_ok=True)
            tmp_v = pipe.vouchers_csv.with_suffix(".csv.tmp")
            shutil.copy2(shared_vouchers, tmp_v)
            tmp_v.replace(pipe.vouchers_csv)

    pipe.lm3_settings = _load_merged_lm3_settings(PROJECT_ROOT)
    ect_cfg = _get_lm3_mod_section(pipe.lm3_settings, "ect")
    seg_cfg = getattr(pipe.config, "segmentation", None)
    pipe.save_debug_crops = (
        bool(save_debug_crops)
        if save_debug_crops is not None
        else bool(getattr(seg_cfg, "save_debug_crops", False) if seg_cfg is not None else False)
    )
    if pipe.save_debug_crops:
        pipe.capitula_patches_dir.mkdir(parents=True, exist_ok=True)
    cfg_ect_dirs = int(getattr(seg_cfg, "ect_num_dirs", 0) or 0) if seg_cfg is not None else 0
    pipe.ect_enabled = bool(ect_cfg.get("enabled", True))
    pipe.ect_num_dirs = max(8, cfg_ect_dirs if cfg_ect_dirs > 0 else int(ect_cfg.get("num_dirs", 360) or 360))
    pipe.ect_num_thetas = max(8, int(ect_cfg.get("num_thresh", ect_cfg.get("num_thetas", 32)) or 32))
    pipe.ect_include_petiole = bool(ect_cfg.get("include_petiole", False))
    pipe.ect_include_holes = bool(ect_cfg.get("include_holes", True))

    pd_cfg = _get_lm3_mod_section(pipe.lm3_settings, "plant_detector")
    ad_cfg = _get_lm3_mod_section(pipe.lm3_settings, "archival_detector")
    ss_cfg = _get_lm3_mod_section(pipe.lm3_settings, "specimen_segmenter")
    ls_cfg = _get_lm3_mod_section(pipe.lm3_settings, "leaf_segmenter")
    kp_cfg = _get_lm3_mod_section(pipe.lm3_settings, "landmark_detector", legacy_alias="pose_detector")

    s1_conf_cfg = float(pd_cfg.get("conf_thres", 0.25) or 0.25)
    s1_iou_cfg = float(pd_cfg.get("iou_thres", 0.55) or 0.55)
    s1_img_size = max(320, int(pd_cfg.get("img_size", 1280) or 1280))
    s1_work_dim = max(640, int(getattr(seg_cfg, "max_working_dim", 0) or pd_cfg.get("max_working_dim", 3200) or 3200))
    ad_conf_cfg = float(ad_cfg.get("conf_thres", 0.25) or 0.25)
    s2_min_area = max(1, int(ls_cfg.get("min_area_px", 64) or 64))
    s2_img_size = max(320, int(ls_cfg.get("img_size", 1024) or 1024))
    kp_conf_cfg = float(kp_cfg.get("conf_thres", 0.25) or 0.25)
    kp_img_size = max(256, int(kp_cfg.get("img_size", 640) or 640))
    kp_pad_frac = float(kp_cfg.get("white_pad_frac", 0.10) if kp_cfg.get("white_pad_frac") is not None else 0.10)
    ss_conf_cfg = float(ss_cfg.get("conf_thres", 0.50) or 0.50)
    ss_img_size = max(320, int(ss_cfg.get("img_size", 1024) or 1024))
    ss_work_dim = max(640, int(getattr(seg_cfg, "specimen_max_working_dim", 0) or ss_cfg.get("max_working_dim", 3200) or 3200))
    pipe.paperclean = bool(getattr(seg_cfg, "paperclean", ss_cfg.get("paperclean", True)) if seg_cfg is not None else ss_cfg.get("paperclean", True))
    pipe.subtract_holes = bool(getattr(seg_cfg, "subtract_holes", True)) if seg_cfg is not None else True
    pipe.do_apply_landmarks = bool(getattr(seg_cfg, "do_apply_landmarks", True)) if seg_cfg is not None else True
    pipe.use_cf_predicted_by_mp = bool(getattr(seg_cfg, "use_cf_predicted_by_mp", False)) if seg_cfg is not None else False

    pipe.ruler_calibrations = RulerCalibrationMap()
    if pipe.lm3_project_dir and pipe.lm3_project_dir.exists():
        pipe.ruler_calibrations.update(load_ruler_calibrations(pipe.lm3_project_dir, min_confidence=pipe.ruler_confidence_thresh))
    if pipe.output_dir.exists():
        pipe.ruler_calibrations.update(load_ruler_calibrations(pipe.output_dir, min_confidence=pipe.ruler_confidence_thresh))

    try:
        is_dummy_w = not pipe.model_weights.stem.lower().startswith("lm3_") and (
            not pipe.model_weights.is_file() or pipe.model_weights.stat().st_size < 1024
        )
    except OSError:
        is_dummy_w = True
    pipe._is_dummy_test_weights = bool(is_dummy_w)
    aux_default_weights = pipe.model_weights if is_dummy_w else None

    ls_conf_cfg = float(ls_cfg.get("conf_thres", 0.35) or 0.35)
    cfg_default_thresh = float(getattr(seg_cfg, "score_thresh", 0.35) or 0.35) if seg_cfg is not None else 0.35
    is_default_score = abs(float(pipe.score_thresh) - cfg_default_thresh) < 1e-6
    eff_s1_conf = s1_conf_cfg if (is_default_score and not is_dummy_w) else min(float(pipe.score_thresh), s1_conf_cfg)
    eff_s2_conf = max(float(pipe.score_thresh), ls_conf_cfg) if (is_default_score and not is_dummy_w) else float(pipe.score_thresh)

    # 1. Stage 1 Plant Detector (YOLO26x 1280px, 9 classes) & 2. Stage 2 Leaf Segmenter (YOLO26x-seg 1024px)
    pipe.stage1_detector = plant_det_cls(
        weights_path=stage1_weights if stage1_weights is not None else aux_default_weights,
        device=pipe.device, conf_thres=eff_s1_conf, iou_thres=s1_iou_cfg,
        img_size=s1_img_size, max_working_dim=s1_work_dim, lm3_project_dir=pipe.lm3_project_dir,
    )
    pipe.plant_detector = pipe.stage1_detector
    pipe.model_engine = leaf_seg_cls(
        weights_path=pipe.model_weights, device=pipe.device,
        score_thresh=eff_s2_conf, min_area_px=s2_min_area, img_size=s2_img_size,
    )
    pipe.leaf_segmenter = pipe.model_engine
    # 3. Pose Landmarks, 4. Ruler/MP Scale Calibrator, 5. Whole-Sheet Specimen Segmenter
    pipe.pose_engine = pose_cls(
        weights_path=aux_default_weights, device=pipe.device,
        conf_thres=kp_conf_cfg, img_size=kp_img_size, white_pad_frac=kp_pad_frac,
    )
    pipe.scale_calibrator = scale_cls(
        device=pipe.device, providers=([] if is_dummy_w else None),
        use_cf_predicted_by_mp=pipe.use_cf_predicted_by_mp, archival_conf_thres=ad_conf_cfg,
    )
    pipe.specimen_segmenter = spec_cls(
        weights_path=aux_default_weights, device=pipe.device, conf_thres=ss_conf_cfg,
        img_size=ss_img_size, max_working_dim=ss_work_dim, paperclean=pipe.paperclean,
        providers=([] if is_dummy_w else None),
    )

    pipe.is_stock_lm3 = is_stock_lm3_weights(pipe.model_weights)
    pipe.voucher_repro_records = {}
    pipe.voucher_pheno_records = {}
    pipe.voucher_repro_by_pid = {}
    pipe.voucher_pheno_by_pid = {}
    pipe.voucher_crowns_by_pid = {}
    pipe.voucher_rosette_centroids_by_pid = {}
    pipe.voucher_petiole_by_pid = {}
    pipe.voucher_contour_coords_by_leaf = {}
    pipe.stage2_batch_size = max(1, int(pipe.gpu_batch_size))

    env_wave_size = _resolve_env_int("PLANTMATICS_WAVE_SIZE", "PACKERA_WAVE_SIZE")
    if wave_size is not None:
        pipe.stage1_wave_size = max(1, int(wave_size))
    elif env_wave_size is not None:
        pipe.stage1_wave_size = env_wave_size
    elif has_explicit_gpu_bs:
        pipe.stage1_wave_size = max(1, min(64, int(pipe.gpu_batch_size)))
    else:
        pipe.stage1_wave_size = 16

    pipe.prefetch_workers_override = max(1, int(prefetch_workers)) if prefetch_workers is not None else _resolve_env_int("PLANTMATICS_PREFETCH_WORKERS", "PACKERA_PREFETCH_WORKERS")
    pipe.prep_workers_override = max(1, int(prep_workers)) if prep_workers is not None else _resolve_env_int("PLANTMATICS_PREP_WORKERS", "PACKERA_PREP_WORKERS")
    pipe.post_workers_override = max(1, int(post_workers)) if post_workers is not None else _resolve_env_int("PLANTMATICS_POST_WORKERS", "PACKERA_POST_WORKERS")

    pipe.vouchers_taxa = {}
    if pipe.vouchers_csv.exists() and pipe.vouchers_csv.stat().st_size > 0:
        try:
            vdf = pd.read_csv(pipe.vouchers_csv)
            if "catalogNumber" in vdf.columns:
                for _, r in vdf.iterrows():
                    cat = str(r["catalogNumber"]).strip()
                    tax = str(r.get("scientificName", r.get("species_raw", ""))).strip()
                    if tax:
                        pipe.vouchers_taxa[cat] = tax
        except Exception:
            pass


def apply_finalize_bundle(pipe: Any, bundle: _VoucherFinalizeBundle) -> None:
    """Merge a completed per-voucher bundle into pipeline state on the coordinator thread."""
    cat = bundle.catalog_number
    if bundle.repro_record is not None:
        pipe.voucher_repro_records[cat] = bundle.repro_record
    if bundle.pheno_record is not None:
        pipe.voucher_pheno_records[cat] = bundle.pheno_record
    if bundle.repro_by_pid:
        pipe.voucher_repro_by_pid[cat] = bundle.repro_by_pid
    if bundle.pheno_by_pid:
        pipe.voucher_pheno_by_pid[cat] = bundle.pheno_by_pid
    if bundle.crowns_by_pid:
        pipe.voucher_crowns_by_pid[cat] = bundle.crowns_by_pid
    if bundle.rosette_centroids_by_pid:
        pipe.voucher_rosette_centroids_by_pid[cat] = bundle.rosette_centroids_by_pid
    if bundle.by_plant_petiole:
        pipe.voucher_petiole_by_pid[cat] = bundle.by_plant_petiole
    if bundle.contour_coords_by_leaf:
        pipe.voucher_contour_coords_by_leaf[cat] = bundle.contour_coords_by_leaf


def index_completed_vouchers(pipe: Any) -> Set[str]:
    """Build an O(1) lookup set of catalogNumbers that already have mask or contour files."""
    completed: Set[str] = set()
    for d in (pipe.contours_dir, pipe.masks_dir):
        if not d.exists():
            continue
        try:
            with os.scandir(d) as it:
                for entry in it:
                    if not entry.is_file():
                        continue
                    name = entry.name
                    if "_leaf" in name and (name.endswith(".csv") or name.endswith(".png")):
                        try:
                            if entry.stat().st_size > 0:
                                completed.add(name.split("_leaf")[0])
                        except OSError:
                            pass
        except OSError:
            pass
    return completed


def check_voucher_completed(pipe: Any, catalog_number: str) -> bool:
    """Verify whether non-empty contour CSV or binary mask artifacts exist for `catalog_number`."""
    contours = list(pipe.contours_dir.glob(f"{catalog_number}_leaf*.csv"))
    if contours and any(c.is_file() and c.stat().st_size > 0 for c in contours):
        return True
    masks = list(pipe.masks_dir.glob(f"{catalog_number}_leaf*.png"))
    return bool(masks and any(m.is_file() and m.stat().st_size > 0 for m in masks))


def purge_voucher_disk_artifacts(pipe: Any, catalog_number: str) -> None:
    """Purge prior mask, contour, ECT, and capitula crop files for `catalog_number` before `--force`."""
    cat = str(catalog_number).strip()
    if not cat:
        return
    roots = [pipe.output_dir]
    if (
        pipe.force
        and os.environ.get("PLANTMATICS_DEMO_MODE") != "1"
        and os.environ.get("PACKERA_DEMO_MODE") != "1"
        and is_canonical_or_explicit_db_context(pipe.output_dir)
    ):
        canon_data = PROJECT_ROOT / "data"
        if canon_data.exists() and canon_data.resolve() != pipe.output_dir.resolve():
            roots.append(canon_data)

    for base_dir in roots:
        for d_path, pat in (
            (base_dir / "contours", f"{cat}_leaf*.csv"),
            (base_dir / "ect", f"{cat}_leaf*_ect.npy"),
            (base_dir / "masks", f"{cat}_leaf*.png"),
            (base_dir / "cropped_patches", f"{cat}_p*.*"),
            (base_dir / "cropped_patches" / "capitula", f"{cat}_cap*.jpg"),
        ):
            if not d_path.exists():
                continue
            for f_old in d_path.glob(pat):
                try:
                    if f_old.is_file() or f_old.is_symlink():
                        f_old.unlink(missing_ok=True)
                except OSError:
                    pass


def atomic_save_csv(df: pd.DataFrame, target_path: Path) -> None:
    """Atomically persist DataFrame to CSV via temporary file replacement."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    temp_file = tempfile.NamedTemporaryFile(
        mode="w", delete=False, dir=target_path.parent, suffix=".tmp", encoding="utf-8",
    )
    temp_path = Path(temp_file.name)
    try:
        df.to_csv(temp_file, index=False, encoding="utf-8")
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
        temp_path.replace(target_path)
    except Exception:
        if temp_path.exists():
            try:
                temp_path.unlink()
            except Exception:
                pass
        raise
