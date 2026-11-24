#!/usr/bin/env python3
"""
Step 02 LeafMachine3 Scale Calibration Engine (`LM3ScaleCalibrator`):
archival ruler detection, classifier verification, lattice scale solving,
and megapixel conversion factor fallback.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

from scripts.vision.lm3_data_loader import compute_work_scale, unscale_working_bbox_to_sheet
from scripts.vision.step02_lm3_lattice_ect import _MP_PRIOR_SLOPE_PX_PER_CM
from scripts.vision.step02_lm3_onnx_runtime import _bgr_to_nchw_float32, _letterbox_bgr

logger = logging.getLogger("SegmentAndExtract")


class LM3ScaleCalibrator:
    """
    LeafMachine3 scale calibration engine combining:
      1. Archival ruler detection (`lm3_archival_detector__yolo26x_det_1280`, class `0: Ruler`)
      2. Ruler Verification Classifier (`RulerClassifierModel`) and Lattice Scale Solver (`RulerCFLattice`)
      3. Megapixel Conversion Factor prior (`MpConversionFactorModel`, `cf_px_per_cm = 27.112084 * sqrt(MP)`)
    """

    ARCHIVAL_MODEL_NAME: str = "lm3_archival_detector__yolo26x_det_1280"
    MP_PRIOR_SLOPE_PX_PER_CM: float = _MP_PRIOR_SLOPE_PX_PER_CM

    def __init__(
        self,
        device: str = "cuda",
        providers: Optional[Sequence[str]] = None,
        use_cf_predicted_by_mp: bool = True,
        archival_conf_thres: float = 0.25,
    ) -> None:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        self.device = device if (device == "cuda" and _lm3_eng.torch_cuda_available()) else "cpu"
        self.providers = providers
        self.use_cf_predicted_by_mp: bool = bool(use_cf_predicted_by_mp)
        self.archival_conf_thres: float = float(archival_conf_thres)
        self.archival_weights_path: Optional[Path] = None
        self.archival_session: Optional[Any] = None
        self.mp_model: Optional[Any] = None
        self.ruler_classifier: Optional[Any] = None
        self.ruler_lattice: Optional[Any] = None
        self._init_lm3_scale_models()
        self.archival_weights_path, _, _, self.archival_session = _lm3_eng._load_onnx_session(
            weights_path=None,
            default_model_name=self.ARCHIVAL_MODEL_NAME,
            device=self.device,
            providers=providers,
        )

    def _init_lm3_scale_models(self) -> None:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        lm3_root = _lm3_eng.LM3_PACKAGE_ROOT
        mp_json = lm3_root / "models" / "mp_conversion_factor" / "model.json"
        if mp_json.is_file():
            try:
                from leafmachine3.inference.mp_conversion_factor import MpConversionFactorModel
                self.mp_model = MpConversionFactorModel(mp_json)
            except Exception:
                self.mp_model = None

        cls_onnx = lm3_root / "models" / "ruler_classifier" / "model.onnx"
        if not cls_onnx.is_file():
            for _cls_stem in (
                "lm3_ruler_classifier__dinov3_vits16_tile_four_896_mlp",
                "lm3_ruler_classifier__effnet_b4",
            ):
                cand_cls = _lm3_eng.resolve_lm3_weights(_cls_stem, explicit_path=None)
                cand_path = Path(cand_cls[1]) if isinstance(cand_cls, tuple) else Path(cand_cls)
                if cand_path.is_file():
                    cls_onnx = cand_path
                    break
        if cls_onnx.is_file():
            try:
                from leafmachine3.inference.ruler_classifier import RulerClassifierModel
                self.ruler_classifier = RulerClassifierModel(cls_onnx)
            except Exception:
                self.ruler_classifier = None

        lat_onnx = lm3_root / "models" / "ruler_cf_lattice" / "model.onnx"
        if not lat_onnx.is_file():
            cand_lat = _lm3_eng.resolve_lm3_weights("lm3_ruler_cf_lattice__effb4_v11", explicit_path=None)
            lat_onnx = Path(cand_lat[1]) if isinstance(cand_lat, tuple) else Path(cand_lat)
        if lat_onnx.is_file():
            try:
                from leafmachine3.inference.ruler_lattice import RulerCFLattice
                self.ruler_lattice = RulerCFLattice(lat_onnx)
            except Exception:
                self.ruler_lattice = None

    def predict_megapixel_scale(self, height_px: int, width_px: int) -> Tuple[float, str]:
        """Predicts `pixels_per_mm` from native sheet dimensions `(height_px, width_px)` via `MpConversionFactorModel`."""
        h, w = max(1, int(height_px)), max(1, int(width_px))
        if self.mp_model is not None:
            try:
                pred = self.mp_model.predict_from_hw(h, w)
                px_per_mm = float(pred.cf_px_per_cm) / 10.0
                if np.isfinite(px_per_mm) and px_per_mm > 0.0:
                    return round(px_per_mm, 4), "predicted_from_megapixels"
            except Exception:
                pass
        mp = (float(h) * float(w)) / 1e6
        px_per_mm = (self.MP_PRIOR_SLOPE_PX_PER_CM * math.sqrt(max(mp, 1e-6))) / 10.0
        return round(px_per_mm, 4), "predicted_from_megapixels"

    def _detect_ruler_boxes(self, image_bgr: np.ndarray) -> List[Tuple[Tuple[int, int, int, int], float]]:
        """Runs `lm3_archival_detector__yolo26x_det_1280` and extracts class `0` (`Ruler`) bounding boxes."""
        import scripts.vision.lm3_inference_engines as _lm3_eng

        if self.archival_session is None or image_bgr is None or image_bgr.size == 0:
            return []
        try:
            sh, sw = image_bgr.shape[:2]
            w_scale = compute_work_scale((sh, sw), max_working_dim=3200)
            if w_scale < 1.0:
                work_w = max(1, int(round(sw * w_scale)))
                work_h = max(1, int(round(sh * w_scale)))
                work_img = cv2.resize(image_bgr, (work_w, work_h), interpolation=cv2.INTER_AREA)
            else:
                work_h, work_w = sh, sw
                work_img = image_bgr
            lb_bgr, r, (dw, dh) = _letterbox_bgr(work_img, 1280)
            inp_tensor = _bgr_to_nchw_float32(lb_bgr)
            input_name = self.archival_session.get_inputs()[0].name
            with _lm3_eng._ONNX_GPU_LOCK:
                raw_outs = self.archival_session.run(None, {input_name: inp_tensor})
            if not raw_outs:
                return []
            arr = np.asarray(raw_outs[0], dtype=np.float32)
            if arr.ndim == 3:
                arr = arr[0]
            if arr.size == 0:
                return []
            inv_r = 1.0 / max(float(r), 1e-6)
            conf_thres = float(getattr(self, "archival_conf_thres", 0.25))

            def _to_sheet_box(
                x1_l: float, y1_l: float, x2_l: float, y2_l: float
            ) -> Optional[Tuple[int, int, int, int]]:
                x1_w = max(0.0, min(float(work_w), (x1_l - dw) * inv_r))
                y1_w = max(0.0, min(float(work_h), (y1_l - dh) * inv_r))
                x2_w = max(0.0, min(float(work_w), (x2_l - dw) * inv_r))
                y2_w = max(0.0, min(float(work_h), (y2_l - dh) * inv_r))
                if y2_w <= y1_w or x2_w <= x1_w:
                    return None
                s_bbox = unscale_working_bbox_to_sheet(
                    (y1_w, x1_w, y2_w, x2_w), w_scale, sheet_shape=(sh, sw)
                )
                return s_bbox if (s_bbox[2] > s_bbox[0] and s_bbox[3] > s_bbox[1]) else None

            rulers: List[Tuple[Tuple[int, int, int, int], float]] = []
            if arr.ndim == 2 and arr.shape[1] in (6, 7):
                for row in arr:
                    sc = float(row[4])
                    cls_id = int(round(float(row[5])))
                    if cls_id != 0 or sc < conf_thres:
                        continue
                    s_bbox = _to_sheet_box(float(row[0]), float(row[1]), float(row[2]), float(row[3]))
                    if s_bbox is not None:
                        rulers.append((s_bbox, sc))
            else:
                if arr.ndim == 2 and arr.shape[0] < arr.shape[1] and arr.shape[0] >= 5:
                    arr = arr.T
                if arr.ndim != 2 or arr.shape[1] < 5:
                    return []
                boxes_xywh = arr[:, :4]
                cls_scores = arr[:, 4:]
                if cls_scores.size == 0:
                    return []
                best_cls = np.argmax(cls_scores, axis=1)
                best_scores = cls_scores[np.arange(len(arr)), best_cls]
                keep = (best_cls == 0) & (best_scores >= conf_thres)
                if not np.any(keep):
                    return []
                boxes_xywh = boxes_xywh[keep]
                best_scores = best_scores[keep]
                nms_boxes = [
                    [float(xc - 0.5 * bw), float(yc - 0.5 * bh), float(bw), float(bh)]
                    for xc, yc, bw, bh in boxes_xywh
                ]
                indices = cv2.dnn.NMSBoxes(
                    bboxes=nms_boxes,
                    scores=[float(s) for s in best_scores],
                    score_threshold=conf_thres,
                    nms_threshold=0.60,
                )
                if len(indices) == 0:
                    return []
                for idx_raw in np.asarray(indices).reshape(-1):
                    idx = int(idx_raw)
                    xc, yc, bw, bh = boxes_xywh[idx]
                    s_bbox = _to_sheet_box(
                        float(xc - 0.5 * bw),
                        float(yc - 0.5 * bh),
                        float(xc + 0.5 * bw),
                        float(yc + 0.5 * bh),
                    )
                    if s_bbox is not None:
                        rulers.append((s_bbox, float(best_scores[idx])))
            rulers.sort(key=lambda item: item[1], reverse=True)
            return rulers
        except Exception as exc:
            logger.debug(f"LM3ScaleCalibrator archival ruler detection fallback: {exc}")
            return []

    def calibrate_sheet(
        self,
        image_bgr: np.ndarray,
        catalog_number: str = "",
        use_cf_predicted_by_mp: Optional[bool] = None,
    ) -> Dict[str, Any]:
        """
        Resolves `pixels_per_mm`, `is_calibrated_mm`, `cf_source`, and `scale_source` for a herbarium sheet:
          1. Attempts direct ruler measurement via `lm3_archival_detector__yolo26x_det_1280` +
             `RulerClassifierModel` + `RulerCFLattice` (`measured_from_ruler`) or verified-crop CPU
             tick-lattice fallback (`ruler_calibrated`).
          2. Falls back to unboxed full-sheet outer-margin tick-lattice detection
             (`cf_source = 'tick_lattice_unboxed'`, `is_calibrated_mm = False`).
          3. When `use_cf_predicted_by_mp` is enabled, falls back to `MpConversionFactorModel`
             (`cf_source = 'predicted_from_megapixels'`, `is_calibrated_mm = False`); otherwise
             returns `uncalibrated` (`pixels_per_mm = None`, `is_calibrated_mm = False`).
        """
        import scripts.vision.lm3_inference_engines as _lm3_eng

        _ = catalog_number
        if image_bgr is None or image_bgr.size == 0:
            return {
                "pixels_per_mm": None,
                "is_calibrated_mm": False,
                "cf_source": "uncalibrated",
                "scale_source": "uncalibrated",
            }

        sh, sw = image_bgr.shape[:2]
        mp_prior_mm = (self.MP_PRIOR_SLOPE_PX_PER_CM * float(np.sqrt((sh * sw) / 1e6))) / 10.0
        if self.archival_session is not None:
            ruler_boxes = self._detect_ruler_boxes(image_bgr)
            min_det_sc = min(0.45, float(getattr(self, "archival_conf_thres", 0.25)))
            _IMPERIAL_TOKENS = ("INCH", "INCHES", "STD_IN", "IMPERIAL", "1/16", "1/8", "1/4")
            for (ymin, xmin, ymax, xmax), _det_sc in ruler_boxes:
                crop = image_bgr[ymin:ymax, xmin:xmax]
                if crop.size == 0 or min(crop.shape[:2]) < 16:
                    continue
                ruler_cls_ok = True
                ruler_cls_confirmed = False
                is_explicit_imperial = False
                lattice_cls_allowed = True
                corrob_tick_mm: Optional[float] = None
                corrob_tick_evaluated = False
                try:
                    if self.ruler_classifier is not None:
                        cls_pred = self.ruler_classifier.predict(crop)
                        if isinstance(cls_pred, dict):
                            lbl_u = str(cls_pred.get("label", "UNKNOWN")).strip().upper()
                            lbl_conf = float(cls_pred.get("conf", 0.0) or 0.0)
                            is_metric_lbl = (
                                lbl_u.startswith("METRIC")
                                or lbl_u in {"GRID_CM", "STD_IN8_CM", "FP"}
                            )
                            if is_metric_lbl and lbl_conf >= 0.35:
                                ruler_cls_confirmed = True
                            is_explicit_imperial = bool(
                                (not is_metric_lbl)
                                and any(tok in lbl_u for tok in _IMPERIAL_TOKENS)
                            )
                            if is_explicit_imperial and lbl_conf >= 0.45:
                                ruler_cls_ok = False
                            elif lbl_u == "MESSY" and lbl_conf >= 0.85 and float(_det_sc) < 0.50:
                                ruler_cls_ok = False
                            if is_explicit_imperial or ((not ruler_cls_confirmed) and lbl_conf >= 0.45):
                                lattice_cls_allowed = False
                        elif not getattr(cls_pred, "is_ruler", True):
                            ruler_cls_ok = False
                    if not ruler_cls_ok:
                        continue
                    if self.ruler_lattice is not None and lattice_cls_allowed:
                        crop_rgb = cv2.cvtColor(crop, cv2.COLOR_BGR2RGB)
                        lat_res = self.ruler_lattice.predict(crop_rgb)
                        px_per_cm = getattr(lat_res, "px_per_cm", None)
                        conf = float(getattr(lat_res, "confidence", 0.0) or 0.0)
                        if px_per_cm is not None and np.isfinite(float(px_per_cm)) and float(px_per_cm) > 10.0 and conf >= 0.45:
                            px_per_mm = round(float(px_per_cm) / 10.0, 4)
                            in_mp_bounds = (0.35 * mp_prior_mm) <= px_per_mm <= (2.25 * mp_prior_mm)
                            if in_mp_bounds:
                                return {
                                    "pixels_per_mm": px_per_mm,
                                    "is_calibrated_mm": True,
                                    "cf_source": "measured_from_ruler",
                                    "scale_source": "ruler_calibrated",
                                }
                            corrob_tick_mm = _lm3_eng.solve_ruler_tick_lattice_px_per_mm(
                                crop,
                                mp_prior_px_per_mm=mp_prior_mm,
                                verified_ruler_crop=True,
                            )
                            corrob_tick_evaluated = True
                            if (
                                corrob_tick_mm is not None
                                and np.isfinite(float(corrob_tick_mm))
                                and float(corrob_tick_mm) > 1.0
                                and abs(px_per_mm - float(corrob_tick_mm)) / max(float(corrob_tick_mm), 1e-6) <= 0.15
                            ):
                                return {
                                    "pixels_per_mm": px_per_mm,
                                    "is_calibrated_mm": True,
                                    "cf_source": "measured_from_ruler",
                                    "scale_source": "ruler_calibrated",
                                }
                except Exception:
                    pass
                if not ruler_cls_ok:
                    continue
                if ruler_cls_confirmed or float(_det_sc) >= min_det_sc:
                    crop_ruler_mm = (
                        corrob_tick_mm
                        if corrob_tick_evaluated
                        else _lm3_eng.solve_ruler_tick_lattice_px_per_mm(
                            crop,
                            mp_prior_px_per_mm=mp_prior_mm,
                            verified_ruler_crop=True,
                        )
                    )
                    if crop_ruler_mm is not None and np.isfinite(float(crop_ruler_mm)) and float(crop_ruler_mm) > 1.0:
                        return {
                            "pixels_per_mm": round(float(crop_ruler_mm), 4),
                            "is_calibrated_mm": True,
                            "cf_source": "ruler_calibrated",
                            "scale_source": "ruler_calibrated",
                        }

        cpu_ruler_mm = _lm3_eng.solve_ruler_tick_lattice_px_per_mm(
            image_bgr,
            mp_prior_px_per_mm=mp_prior_mm,
            verified_ruler_crop=False,
        )
        if cpu_ruler_mm is not None and np.isfinite(float(cpu_ruler_mm)) and float(cpu_ruler_mm) > 1.0:
            return {
                "pixels_per_mm": round(float(cpu_ruler_mm), 4),
                "is_calibrated_mm": False,
                "cf_source": "tick_lattice_unboxed",
                "scale_source": "tick_lattice_unboxed",
            }

        eff_use_mp = (
            self.use_cf_predicted_by_mp
            if use_cf_predicted_by_mp is None
            else bool(use_cf_predicted_by_mp)
        )
        if not eff_use_mp:
            return {
                "pixels_per_mm": None,
                "is_calibrated_mm": False,
                "cf_source": "uncalibrated",
                "scale_source": "uncalibrated",
            }
        px_per_mm, cf_source = self.predict_megapixel_scale(sh, sw)
        return {
            "pixels_per_mm": px_per_mm,
            "is_calibrated_mm": cf_source != "predicted_from_megapixels",
            "cf_source": cf_source,
            "scale_source": "mp_estimated" if cf_source == "predicted_from_megapixels" else "ruler_calibrated",
        }
