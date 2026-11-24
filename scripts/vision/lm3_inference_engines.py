#!/usr/bin/env python3
"""
Step 02 LeafMachine3 ONNX Runtime inference facade re-exporting modular engines,
runtime state, pose/orientation helpers, scale calibration, and ECT descriptors.
"""

from __future__ import annotations

import logging
import sys
import threading
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import scripts.vision.lm3_compat  # noqa: F401

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
LM3_PACKAGE_ROOT = PROJECT_ROOT / "LeafMachine3"
if LM3_PACKAGE_ROOT.is_dir() and str(LM3_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(LM3_PACKAGE_ROOT))
_VENV_LM3_LIB = PROJECT_ROOT / ".venv_LM3" / "lib"


def _ensure_ect_module_loaded() -> None:
    if "ect" in sys.modules:
        return
    try:
        import ect as _ect_probe  # noqa: F401
        return
    except ImportError:
        pass
    if not _VENV_LM3_LIB.is_dir():
        return
    import importlib.util as _ilu
    import types as _types

    for _ect_init in sorted(_VENV_LM3_LIB.glob("python*/site-packages/ect/__init__.py")):
        if _ect_init.is_file():
            try:
                import numba as _numba_probe  # noqa: F401
            except ImportError:
                _numba_mod = _types.ModuleType("numba")
                _numba_mod.prange = range  # type: ignore[attr-defined]
                _numba_mod.njit = (  # type: ignore[attr-defined]
                    lambda *a, **kw: (a[0] if (a and callable(a[0])) else (lambda f: f))
                )
                _numba_typed = _types.ModuleType("numba.typed")
                _numba_typed.List = list  # type: ignore[attr-defined]
                _numba_mod.typed = _numba_typed  # type: ignore[attr-defined]
                sys.modules["numba"] = _numba_mod
                sys.modules["numba.typed"] = _numba_typed
            _ect_spec = _ilu.spec_from_file_location(
                "ect",
                str(_ect_init),
                submodule_search_locations=[str(_ect_init.parent)],
            )
            if _ect_spec is not None and _ect_spec.loader is not None:
                _ect_mod = _ilu.module_from_spec(_ect_spec)
                sys.modules["ect"] = _ect_mod
                _ect_spec.loader.exec_module(_ect_mod)
            break


try:
    from scripts.vision.configure_leafmachine3 import resolve_lm3_weights
except ImportError:
    def resolve_lm3_weights(
        default_filename: str,
        explicit_path: Optional[Union[str, Path]] = None,
        project_root: Optional[Path] = None,
    ) -> Path:
        root = Path(project_root) if project_root is not None else PROJECT_ROOT
        if explicit_path is not None and str(explicit_path).strip():
            p = Path(explicit_path)
            return p if p.is_absolute() else (root / p)
        stem = Path(default_filename).stem
        return root / "models" / "lm3" / f"{stem}.onnx"


from scripts.vision.lm3_data_loader import (  # noqa: E402
    LM3_LEAF_CLASS_IDS,
    LM3_PLANT_CLASSES,
    LM3_REPRODUCTIVE_CLASS_IDS,
    compute_work_scale,
    load_precomputed_plant_labels,
    unscale_working_bbox_to_sheet,
)

logger = logging.getLogger("SegmentAndExtract")

DEFAULT_ONNX_PROVIDERS: Tuple[str, ...] = (
    "CUDAExecutionProvider",
    "CoreMLExecutionProvider",
    "CPUExecutionProvider",
)

LM3_KEYPOINT_NAMES_31: Tuple[str, ...] = (
    "lamina_tip", "apex_1", "apex_2", "apex_3",
    "midvein_1", "midvein_2", "midvein_3", "midvein_4", "midvein_5",
    "midvein_6", "midvein_7", "midvein_8", "midvein_9", "midvein_10",
    "midvein_11", "midvein_12", "midvein_13", "midvein_14", "midvein_15",
    "base_1", "base_2", "base_3", "lamina_base",
    "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_5", "petiole_6",
    "width_1", "width_2",
)

_ONNX_SESSION_CACHE: Dict[Tuple[str, Tuple[str, ...]], Any] = {}
_NVIDIA_LIBS_PRELOADED: bool = False
_SHARED_CUDA_ARENA_REGISTERED: bool = False
_ONNX_GPU_LOCK: threading.Lock = threading.Lock()

from scripts.vision.step02_lm3_bbox_capitula import (  # noqa: E402
    DetectedInstance,
    _VoucherFinalizeBundle,
    compute_capitulum_metrics,
    crop_leaf_bbox_with_padding,
    reproject_crop_mask_to_sheet,
)
from scripts.vision.step02_lm3_crop_reducer import (  # noqa: E402
    _reduce_crop_predictions_native_lm3,
)
from scripts.vision.step02_lm3_damage_qc import (  # noqa: E402
    LM3SpecimenSegmenter,
    _count_half_width_lobes,
    evaluate_leaf_damage_and_fold_qc,
)
from scripts.vision.step02_lm3_lattice_ect import (  # noqa: E402
    _acorr_unbiased_1d,
    _parabolic_peak_lag,
    _scan_margin_strips_for_mm_lattice,
    compute_leaf_ect_descriptor,
    solve_ruler_tick_lattice_px_per_mm,
)
from scripts.vision.step02_lm3_leaf_orienter import (  # noqa: E402
    orient_and_score_leaf_lm3,
)
from scripts.vision.step02_lm3_leaf_segmenter import (  # noqa: E402
    LM3LeafSegmenterEngine,
)
from scripts.vision.step02_lm3_onnx_runtime import (  # noqa: E402
    _INTER_LINEAR_EXACT,
    _bgr_to_nchw_float32,
    _ensure_shared_cuda_arena,
    _letterbox_bgr,
    _letterbox_gain_pad_hw,
    _load_onnx_session,
    _preload_nvidia_cuda_libs,
    _prepare_working_image,
    _run_onnx_minibatches_back_to_back,
    is_stock_lm3_weights,
    resolve_onnx_providers,
    torch_cuda_available,
)
from scripts.vision.step02_lm3_petiole_kpts import (  # noqa: E402
    _LM3_CANONICAL_TO_NATIVE_KPT_MAP,
    _extract_native_lm3_petiole_anchor,
    normalize_lm3_keypoints_dict,
)
from scripts.vision.step02_lm3_plant_detector import (  # noqa: E402
    LM3PlantDetector,
)
from scripts.vision.step02_lm3_pose_engine import (  # noqa: E402
    LM3PoseLandmarkEngine,
)
from scripts.vision.step02_lm3_pose_fallback import (  # noqa: E402
    LM3PoseResult,
)
from scripts.vision.step02_lm3_scale_calibrator import (  # noqa: E402
    LM3ScaleCalibrator,
)

__all__ = [
    "DEFAULT_ONNX_PROVIDERS",
    "LM3_KEYPOINT_NAMES_31",
    "DetectedInstance",
    "_VoucherFinalizeBundle",
    "crop_leaf_bbox_with_padding",
    "reproject_crop_mask_to_sheet",
    "compute_capitulum_metrics",
    "resolve_onnx_providers",
    "LM3PlantDetector",
    "LM3SpecimenSegmenter",
    "LM3LeafSegmenterEngine",
    "LM3PoseResult",
    "LM3PoseLandmarkEngine",
    "LM3ScaleCalibrator",
    "normalize_lm3_keypoints_dict",
    "orient_and_score_leaf_lm3",
    "compute_leaf_ect_descriptor",
    "torch_cuda_available",
    "is_stock_lm3_weights",
    "_prepare_working_image",
    "_reduce_crop_predictions_native_lm3",
    "_extract_native_lm3_petiole_anchor",
    "evaluate_leaf_damage_and_fold_qc",
]
