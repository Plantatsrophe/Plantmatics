#!/usr/bin/env python3
"""
Step 02 LeafMachine3 ONNX Runtime session management, CUDA allocator/library setup,
working-frame/letterbox preprocessing, and back-to-back GPU mini-batch execution.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np

from scripts.vision.lm3_data_loader import compute_work_scale

logger = logging.getLogger("SegmentAndExtract")

_INTER_LINEAR_EXACT: int = int(getattr(cv2, "INTER_LINEAR_EXACT", cv2.INTER_LINEAR))


def torch_cuda_available() -> bool:
    """Safely checks whether CUDA is available via ONNX Runtime or PyTorch."""
    try:
        import onnxruntime as ort
        if "CUDAExecutionProvider" in ort.get_available_providers():
            return True
    except Exception:
        pass
    try:
        import torch
        return bool(torch.cuda.is_available())
    except Exception:
        return False


def is_stock_lm3_weights(weights_path: Union[str, Path]) -> bool:
    """
    Returns True when using canonical stock LeafMachine3 ONNX weights rather than
    a taxon-specific custom fine-tuned checkpoint.
    """
    name = Path(weights_path).name.lower()
    if not name:
        return True
    return "finetuned" not in name and "packera" not in name


def resolve_onnx_providers(
    device: str = "cuda",
    requested_providers: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    Resolves the ordered ONNX Runtime execution provider priority list:
      ["CUDAExecutionProvider", "CoreMLExecutionProvider", "CPUExecutionProvider"]
    filtered against available runtime providers when device != 'cpu'.
    """
    import scripts.vision.lm3_inference_engines as _lm3_eng

    priority = list(requested_providers) if requested_providers else list(_lm3_eng.DEFAULT_ONNX_PROVIDERS)
    if str(device).strip().lower() == "cpu":
        return ["CPUExecutionProvider"]
    try:
        import onnxruntime as ort
        avail = set(ort.get_available_providers())
        if avail:
            filtered = [p for p in priority if p in avail]
            if "CPUExecutionProvider" not in filtered:
                filtered.append("CPUExecutionProvider")
            return filtered if filtered else list(_lm3_eng.DEFAULT_ONNX_PROVIDERS)
    except Exception:
        pass
    return priority


def _ensure_shared_cuda_arena() -> bool:
    """
    Registers a single process-wide shared CUDA BFCArena allocator in ONNX Runtime so all
    4 LeafMachine3 ONNX sessions (PlantDetector, LeafSegmenter, PoseLandmark, SpecimenSegmenter)
    reuse the same GPU activation pool instead of hoarding 4 isolated per-session arenas.
    """
    import scripts.vision.lm3_inference_engines as _lm3_eng

    if _lm3_eng._SHARED_CUDA_ARENA_REGISTERED:
        return True
    try:
        import onnxruntime as ort
        mem_info = ort.OrtMemoryInfo(
            "Cuda",
            ort.OrtAllocatorType.ORT_ARENA_ALLOCATOR,
            0,
            ort.OrtMemType.DEFAULT,
        )
        # OrtArenaCfg(max_mem=0 [auto], arena_extend_strategy=1 [kSameAsRequested], initial_chunk=-1, max_dead=-1)
        arena_cfg = ort.OrtArenaCfg(0, 1, -1, -1)
        if hasattr(ort, "create_and_register_allocator_v2"):
            ort.create_and_register_allocator_v2("CUDAExecutionProvider", mem_info, {}, arena_cfg)
        else:
            ort.create_and_register_allocator(mem_info, arena_cfg)
        _lm3_eng._SHARED_CUDA_ARENA_REGISTERED = True
        return True
    except Exception as exc:
        logger.debug(f"Shared CUDA arena registration skipped: {exc}")
        return False


def _preload_nvidia_cuda_libs() -> None:
    """
    Preloads bundled `nvidia-*` wheel CUDA 12 / cuDNN 9 shared libraries into the process
    symbol table (`RTLD_GLOBAL`) so `onnxruntime-gpu` binds `CUDAExecutionProvider` even when
    `LD_LIBRARY_PATH` was not exported prior to interpreter launch.
    """
    import scripts.vision.lm3_inference_engines as _lm3_eng

    if _lm3_eng._NVIDIA_LIBS_PRELOADED:
        return
    _lm3_eng._NVIDIA_LIBS_PRELOADED = True
    try:
        import ctypes
        import nvidia  # type: ignore
        roots = [Path(p) for p in getattr(nvidia, "__path__", ()) or ()]
        lib_dirs = sorted({d for root in roots for d in root.glob("*/lib") if d.is_dir()})
        so_files = [
            p
            for d in lib_dirs
            for p in sorted(d.glob("lib*.so*"))
            if p.is_file() and not p.is_symlink()
        ]
        for _ in range(2):
            for p in so_files:
                try:
                    ctypes.CDLL(str(p), mode=ctypes.RTLD_GLOBAL)
                except OSError:
                    pass
    except Exception:
        pass


def _load_onnx_session(
    weights_path: Optional[Union[str, Path]],
    default_model_name: str,
    device: str = "cuda",
    providers: Optional[Sequence[str]] = None,
) -> Tuple[Path, str, List[str], Optional[Any]]:
    """Loads or retrieves a cached ONNX Runtime InferenceSession with graceful fallback when absent."""
    import scripts.vision.lm3_inference_engines as _lm3_eng

    lookup_name = default_model_name
    if (
        (weights_path is None or not str(weights_path).strip())
        and default_model_name == "lm3_landmark_detector__yolo26x_pose_640"
    ):
        cand_canon = _lm3_eng.resolve_lm3_weights(
            "lm3_landmark_detector__yolo26x_pose_640_canonical", explicit_path=None
        )
        canon_path = Path(cand_canon[1]) if isinstance(cand_canon, tuple) else Path(cand_canon)
        try:
            if canon_path.is_file() and canon_path.stat().st_size >= 1024:
                lookup_name = "lm3_landmark_detector__yolo26x_pose_640_canonical"
        except OSError:
            pass
    resolved = _lm3_eng.resolve_lm3_weights(lookup_name, explicit_path=weights_path)
    if isinstance(resolved, tuple):
        m_name, w_path = str(resolved[0]), Path(resolved[1])
    else:
        w_path = Path(resolved)
        if weights_path is not None and str(weights_path).strip():
            wp_stem = Path(weights_path).stem
            m_name = wp_stem if (wp_stem and wp_stem != "model") else (Path(default_model_name).stem or default_model_name)
        else:
            m_name = (w_path.stem if (w_path.stem and w_path.stem != "model") else Path(default_model_name).stem) or default_model_name

    prov_list = _lm3_eng.resolve_onnx_providers(device=device, requested_providers=providers)
    try:
        has_valid_weights = w_path.is_file() and w_path.stat().st_size >= 1024
    except OSError:
        has_valid_weights = False

    if not has_valid_weights:
        logger.debug(f"LeafMachine3 ONNX weights not found at {w_path}; running in fallback/mock-ready mode.")
        return w_path, m_name, prov_list, None

    cache_key = (str(w_path.resolve()), tuple(prov_list))
    if cache_key in _lm3_eng._ONNX_SESSION_CACHE:
        return w_path, m_name, prov_list, _lm3_eng._ONNX_SESSION_CACHE[cache_key]

    try:
        import os
        import onnxruntime as ort
        if "CUDAExecutionProvider" in prov_list:
            _lm3_eng._preload_nvidia_cuda_libs()
        sess_opts = ort.SessionOptions()
        sess_opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        sess_opts.enable_mem_pattern = False
        sess_opts.enable_mem_reuse = True
        if "CUDAExecutionProvider" in prov_list and _lm3_eng._ensure_shared_cuda_arena():
            try:
                sess_opts.add_session_config_entry("session.use_env_allocators", "1")
            except Exception:
                pass
        sess_opts.intra_op_num_threads = 2
        sess_opts.inter_op_num_threads = 1
        try:
            sess_opts.log_severity_level = int(os.environ.get("LM3_ORT_LOG_SEVERITY", 3))
        except ValueError:
            sess_opts.log_severity_level = 3
        ort_providers: List[Any] = [
            (
                "CUDAExecutionProvider",
                {
                    "device_id": 0,
                    "arena_extend_strategy": "kSameAsRequested",
                    "cudnn_conv_use_max_workspace": "0",
                    "do_copy_in_default_stream": "1",
                },
            )
            if p == "CUDAExecutionProvider"
            else p
            for p in prov_list
        ]
        session = ort.InferenceSession(str(w_path), sess_options=sess_opts, providers=ort_providers)
        _lm3_eng._ONNX_SESSION_CACHE[cache_key] = session
        bound_provs = session.get_providers()
        logger.info(f"Loaded LeafMachine3 ONNX session ({m_name}) from {w_path} with providers={bound_provs}")
        return w_path, m_name, bound_provs, session
    except Exception as exc:
        logger.debug(f"ONNX Runtime session initialization skipped for {w_path}: {exc}")
        return w_path, m_name, prov_list, None


def _prepare_working_image(
    image_bgr: np.ndarray,
    max_working_dim: int = 3200,
) -> Tuple[np.ndarray, float]:
    """
    Downscales `image_bgr` with `cv2.INTER_AREA` when `max(orig_h, orig_w) > max_working_dim`
    (default 3200 px), returning `(work_img, work_scale)` matching LeafMachine3's working-image stage.
    """
    if image_bgr is None or image_bgr.size == 0:
        return image_bgr, 1.0
    orig_h, orig_w = int(image_bgr.shape[0]), int(image_bgr.shape[1])
    work_scale = compute_work_scale((orig_h, orig_w), max_working_dim=int(max_working_dim))
    if work_scale < 1.0:
        work_w = max(1, int(round(orig_w * work_scale)))
        work_h = max(1, int(round(orig_h * work_scale)))
        work_img = cv2.resize(image_bgr, (work_w, work_h), interpolation=cv2.INTER_AREA)
        return work_img, float(work_scale)
    return image_bgr, 1.0


def _letterbox_gain_pad_hw(
    img1_shape: Tuple[int, int],
    img0_shape: Tuple[int, int],
) -> Tuple[float, Tuple[float, float]]:
    """Computes LeafMachine3 `ultra_rep._letterbox_gain_pad` `(gain, (pad_x, pad_y))`."""
    try:
        from leafmachine3.inference import ultra_replacements as ultra_rep
        gain, pad_x, pad_y = ultra_rep._letterbox_gain_pad(img1_shape, img0_shape)
        return float(gain), (float(pad_x), float(pad_y))
    except Exception:
        h0, w0 = max(1, int(img0_shape[0])), max(1, int(img0_shape[1]))
        h1, w1 = int(img1_shape[0]), int(img1_shape[1])
        gain = min(float(h1) / float(h0), float(w1) / float(w0))
        pad_x = float(round((w1 - round(w0 * gain)) / 2.0 - 0.1))
        pad_y = float(round((h1 - round(h0 * gain)) / 2.0 - 0.1))
        return gain, (pad_x, pad_y)


def _letterbox_bgr(
    image_bgr: np.ndarray,
    target_size: int,
    pad_value: int = 114,
    stride: int = 32,
    auto: bool = True,
) -> Tuple[np.ndarray, float, Tuple[float, float]]:
    """
    Resizes `image_bgr` using LeafMachine3's native dynamic stride-32 letterboxing
    (`ultra_rep.letterbox` with `auto=True`, `stride=32`, and `cv2.INTER_LINEAR_EXACT`).
    Returns `(letterboxed_bgr, scale_ratio, (dw, dh))`.
    """
    h, w = image_bgr.shape[:2]
    if h <= 0 or w <= 0:
        empty = np.full((target_size, target_size, 3), pad_value, dtype=np.uint8)
        return empty, 1.0, (0.0, 0.0)
    try:
        from leafmachine3.inference import ultra_replacements as ultra_rep
        in_hw = ultra_rep.check_imgsz(int(target_size), stride=int(stride))
        padded = ultra_rep.letterbox(
            image_bgr,
            in_hw,
            stride=int(stride),
            auto=bool(auto),
            center=True,
        )
        gain, pad_x, pad_y = ultra_rep._letterbox_gain_pad(padded.shape[:2], (h, w))
        return padded, float(gain), (float(pad_x), float(pad_y))
    except Exception:
        r = min(float(target_size) / float(h), float(target_size) / float(w))
        new_w = max(1, int(round(w * r)))
        new_h = max(1, int(round(h * r)))
        dw = float(target_size - new_w)
        dh = float(target_size - new_h)
        if auto:
            dw = float(np.mod(dw, stride))
            dh = float(np.mod(dh, stride))
        dw *= 0.5
        dh *= 0.5
        top = int(round(dh - 0.1))
        bottom = int(round(dh + 0.1))
        left = int(round(dw - 0.1))
        right = int(round(dw + 0.1))
        if (new_w, new_h) != (w, h):
            resized = cv2.resize(image_bgr, (new_w, new_h), interpolation=cv2.INTER_LINEAR_EXACT)
        else:
            resized = image_bgr
        padded = cv2.copyMakeBorder(
            resized, top, bottom, left, right, cv2.BORDER_CONSTANT, value=(pad_value, pad_value, pad_value)
        )
        return padded, r, (float(left), float(top))


def _bgr_to_nchw_float32(letterboxed_bgr: np.ndarray) -> np.ndarray:
    """
    Converts a `(H, W, 3)` BGR uint8 array to `(1, 3, H, W)` RGB float32 in `[0, 1]`
    using LeafMachine3's reciprocal float32 scaling (`x *= float32(1.0) / float32(255.0)`).
    """
    im = letterboxed_bgr[None, ..., ::-1].transpose((0, 3, 1, 2))
    im_f32 = np.ascontiguousarray(im, dtype=np.float32)
    im_f32 *= np.float32(1.0) / np.float32(255.0)
    return im_f32


def _run_onnx_minibatches_back_to_back(
    active_sess: Any,
    input_name: str,
    valid_chw_list: Sequence[np.ndarray],
    initial_sub_bs: int,
) -> List[Tuple[int, int, List[np.ndarray]]]:
    """
    Executes GPU ONNX mini-batches back-to-back without interleaved CPU post-processing.
    Automatically halves the mini-batch size on transient BFCArena / VRAM allocation pressure
    and respects contiguous `(H_lb, W_lb)` shape boundaries.
    """
    import scripts.vision.lm3_inference_engines as _lm3_eng

    n_valid = len(valid_chw_list)
    if n_valid == 0:
        return []
    sub_bs = max(1, int(initial_sub_bs))
    batch_records: List[Tuple[int, int, List[np.ndarray]]] = []
    idx = 0
    while idx < n_valid:
        base_shape = valid_chw_list[idx].shape
        max_same_shape = 1
        while (
            idx + max_same_shape < n_valid
            and max_same_shape < sub_bs
            and valid_chw_list[idx + max_same_shape].shape == base_shape
        ):
            max_same_shape += 1
        cur_bs = min(sub_bs, max_same_shape)
        batch_tensor = np.stack(valid_chw_list[idx : idx + cur_bs], axis=0)
        try:
            with _lm3_eng._ONNX_GPU_LOCK:
                raw_outs = active_sess.run(None, {input_name: batch_tensor})
        except Exception as exc:
            if cur_bs > 1:
                sub_bs = max(1, cur_bs // 2)
                logger.debug(
                    f"ONNX mini-batch size {cur_bs} hit memory ceiling ({exc}); retrying with sub_bs={sub_bs}"
                )
                continue
            raise
        if raw_outs:
            batch_records.append((idx, idx + cur_bs, raw_outs))
        idx += cur_bs
    return batch_records
