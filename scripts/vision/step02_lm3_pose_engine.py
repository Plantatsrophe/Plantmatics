#!/usr/bin/env python3
"""
Step 02 LeafMachine3 Stage 2 Pose Landmark Engine (YOLO26x-pose, 640 px, 31 keypoints).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np

from scripts.vision.step02_lm3_onnx_runtime import (
    _bgr_to_nchw_float32, _letterbox_bgr, _run_onnx_minibatches_back_to_back,
)
from scripts.vision.step02_lm3_pose_fallback import (
    _KEYPOINT_NAMES_31, LM3PoseResult, estimate_petiole_width_px, geometric_fallback_keypoints_31,
)

logger = logging.getLogger("SegmentAndExtract")

_PoseMeta = Tuple[int, int, int, int, float, float, float, float, float]


class LM3PoseLandmarkEngine:
    """
    LeafMachine3 Stage 2 ONNX Pose Landmark Engine (`lm3_landmark_detector__yolo26x_pose_640`).
    Runs 640 px YOLO26x-pose to predict the 31-keypoint `mid15_pet5` anatomical leaf skeleton
    (`lamina_tip`, `apex_1..3`, `midvein_1..15`, `base_1..3`, `lamina_base`, `petiole_1..6`,
    `width_1..2`) plus estimated `petiole_width_px`.
    """

    MODEL_NAME: str = "lm3_landmark_detector__yolo26x_pose_640"
    KEYPOINT_NAMES: Tuple[str, ...] = _KEYPOINT_NAMES_31
    WHITE_PAD_FRAC: float = 0.10

    def __init__(
        self, weights_path: Optional[Union[str, Path]] = None, device: str = "cuda",
        conf_thres: float = 0.25, img_size: int = 640, white_pad_frac: float = 0.10,
        providers: Optional[Sequence[str]] = None,
    ) -> None:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        self.device = device if (device == "cuda" and _lm3_eng.torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.img_size = int(img_size)
        self.WHITE_PAD_FRAC = float(white_pad_frac)
        self.weights_path, self.model_name, self.providers, self.session = _lm3_eng._load_onnx_session(
            weights_path=weights_path, default_model_name=self.MODEL_NAME, device=self.device, providers=providers,
        )
        self.model = self.session

    @staticmethod
    def _estimate_petiole_width_px(
        keypoints_31: Dict[str, Tuple[float, float, float]],
        petiole_mask: Optional[np.ndarray] = None, lamina_mask: Optional[np.ndarray] = None,
    ) -> Optional[float]:
        """Estimates petiole width in pixels by delegating to `estimate_petiole_width_px`."""
        return estimate_petiole_width_px(keypoints_31, petiole_mask=petiole_mask, lamina_mask=lamina_mask)

    @staticmethod
    def _geometric_fallback_keypoints_31(
        crop_bgr: Optional[np.ndarray], lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None, offset_xy: Tuple[float, float] = (0.0, 0.0),
        crown_hint_xy: Optional[Tuple[float, float]] = None,
    ) -> LM3PoseResult:
        """Derives the 31 canonical `mid15_pet5` keypoints deterministically from masks/crop."""
        return geometric_fallback_keypoints_31(
            crop_bgr=crop_bgr, lamina_mask=lamina_mask, petiole_mask=petiole_mask,
            offset_xy=offset_xy, crown_hint_xy=crown_hint_xy,
        )

    def _pad_and_letterbox(
        self, crop_bgr: np.ndarray, pad_f: float
    ) -> Tuple[np.ndarray, int, int, int, int, float, float, float]:
        crop_h, crop_w = crop_bgr.shape[:2]
        if pad_f > 0.0:
            px, py = int(round(crop_w * pad_f)), int(round(crop_h * pad_f))
            pose_in_bgr = cv2.copyMakeBorder(crop_bgr, py, py, px, px, cv2.BORDER_CONSTANT, value=(255, 255, 255))
        else:
            px, py = 0, 0
            pose_in_bgr = crop_bgr
        lb_bgr, r, (dw, dh) = _letterbox_bgr(pose_in_bgr, self.img_size, auto=True, stride=32)
        return lb_bgr, crop_h, crop_w, px, py, float(r), float(dw), float(dh)

    def _decode_pose_array(
        self,
        raw_arr: np.ndarray,
        meta: _PoseMeta,
        lamina_mask: Optional[np.ndarray],
        petiole_mask: Optional[np.ndarray],
    ) -> Optional[LM3PoseResult]:
        crop_h, crop_w, px, py, r, dw, dh, ox, oy = meta
        arr = np.asarray(raw_arr, dtype=np.float32)
        if arr.ndim == 3:
            arr = arr[0]
        if arr.ndim == 2 and arr.shape[0] < arr.shape[1] and arr.shape[0] >= 5 + 3 * 31:
            arr = arr.T
        if arr.ndim == 2 and arr.shape[0] > 0 and arr.shape[1] >= 5 + 3 * 31:
            best_row = arr[int(np.argmax(arr[:, 4]))]
            if float(best_row[4]) >= self.conf_thres:
                kp_offset = 6 if best_row.shape[0] >= 6 + 3 * 31 else 5
                kp_flat = best_row[kp_offset : kp_offset + 3 * 31].reshape(31, 3)
                inv_r = float(np.float32(1.0) / np.float32(max(float(r), 1e-6)))
                kps_out: Dict[str, Tuple[float, float, float]] = {}
                for idx_k, kp_name in enumerate(self.KEYPOINT_NAMES):
                    kx_pad = (float(kp_flat[idx_k, 0]) - dw) * inv_r
                    ky_pad = (float(kp_flat[idx_k, 1]) - dh) * inv_r
                    kx = max(0.0, min(float(crop_w), kx_pad - float(px))) + ox
                    ky = max(0.0, min(float(crop_h), ky_pad - float(py))) + oy
                    kc = float(np.clip(kp_flat[idx_k, 2], 0.0, 1.0))
                    kps_out[kp_name] = (round(kx, 2), round(ky, 2), round(kc, 4))
                pet_w = self._estimate_petiole_width_px(kps_out, petiole_mask=petiole_mask, lamina_mask=lamina_mask)
                return LM3PoseResult(kps_out, petiole_width_px=pet_w)
        return None

    def predict_landmarks(
        self,
        crop_bgr: np.ndarray,
        lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None,
        crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        """
        Runs 640 px YOLO26x-pose landmark detection on unmasked `crop_bgr`, returning an `LM3PoseResult`
        mapping all 31 canonical keypoint names to `(x, y, confidence)` plus `.petiole_width_px`.
        """
        import scripts.vision.lm3_inference_engines as _lm3_eng

        ox, oy = float(crop_offset[0]), float(crop_offset[1])
        active_sess = self.session if self.session is not None else self.model
        if active_sess is not None and crop_bgr is not None and crop_bgr.size > 0:
            try:
                pad_f = float(getattr(self, "WHITE_PAD_FRAC", 0.10))
                lb_bgr, crop_h, crop_w, px, py, r, dw, dh = self._pad_and_letterbox(crop_bgr, pad_f)
                inp_tensor = _bgr_to_nchw_float32(lb_bgr)
                input_name = active_sess.get_inputs()[0].name
                with _lm3_eng._ONNX_GPU_LOCK:
                    raw_outs = active_sess.run(None, {input_name: inp_tensor})
                if raw_outs:
                    decoded = self._decode_pose_array(
                        raw_outs[0], (crop_h, crop_w, px, py, r, dw, dh, ox, oy), lamina_mask, petiole_mask
                    )
                    if decoded is not None:
                        return decoded
            except Exception as exc:
                logger.debug(f"LM3PoseLandmarkEngine ONNX inference fallback: {exc}")

        return self._geometric_fallback_keypoints_31(
            crop_bgr=crop_bgr, lamina_mask=lamina_mask, petiole_mask=petiole_mask, offset_xy=(ox, oy),
        )

    def _predict_batch_sequential(
        self, crops_bgr: Sequence[np.ndarray],
        lamina_masks: Optional[Sequence[Optional[np.ndarray]]],
        petiole_masks: Optional[Sequence[Optional[np.ndarray]]],
        crop_offsets: Optional[Sequence[Tuple[float, float]]],
    ) -> List[LM3PoseResult]:
        return [
            self.predict_landmarks(
                crop_bgr=crops_bgr[i],
                lamina_mask=lamina_masks[i] if (lamina_masks and i < len(lamina_masks)) else None,
                petiole_mask=petiole_masks[i] if (petiole_masks and i < len(petiole_masks)) else None,
                crop_offset=crop_offsets[i] if (crop_offsets and i < len(crop_offsets)) else (0.0, 0.0),
            )
            for i in range(len(crops_bgr))
        ]

    def predict_batch(
        self, crops_bgr: Sequence[np.ndarray],
        lamina_masks: Optional[Sequence[Optional[np.ndarray]]] = None,
        petiole_masks: Optional[Sequence[Optional[np.ndarray]]] = None,
        crop_offsets: Optional[Sequence[Tuple[float, float]]] = None,
        batch_size: int = 8, executor: Optional[Any] = None,
    ) -> List[LM3PoseResult]:
        """Runs batched 640 px YOLO26x-pose landmark detection in shape-bucketed `(B, 3, H_lb, W_lb)` mini-batches."""
        if not crops_bgr:
            return []
        n_items = len(crops_bgr)
        active_sess = self.session if self.session is not None else self.model
        if (
            active_sess is None
            or "predict_landmarks" in self.__dict__
            or getattr(self.predict_landmarks, "__func__", None) is not LM3PoseLandmarkEngine.predict_landmarks
            or n_items <= 1
        ):
            return self._predict_batch_sequential(crops_bgr, lamina_masks, petiole_masks, crop_offsets)

        try:
            input_name = active_sess.get_inputs()[0].name
            sub_bs = min(16, max(1, int(batch_size)))
            pad_f = float(getattr(self, "WHITE_PAD_FRAC", 0.10))

            def _offset_at(idx: int) -> Tuple[float, float]:
                if crop_offsets and idx < len(crop_offsets):
                    return (float(crop_offsets[idx][0]), float(crop_offsets[idx][1]))
                return (0.0, 0.0)

            def _prep_pose_item(idx: int) -> Optional[Tuple[int, np.ndarray, _PoseMeta]]:
                crop_bgr = crops_bgr[idx]
                ox, oy = _offset_at(idx)
                if crop_bgr is None or crop_bgr.size == 0:
                    return None
                lb_bgr, crop_h, crop_w, px, py, r, dw, dh = self._pad_and_letterbox(crop_bgr, pad_f)
                rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
                chw_f32 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.float32)
                chw_f32 *= np.float32(1.0) / np.float32(255.0)
                return (idx, chw_f32, (crop_h, crop_w, px, py, r, dw, dh, ox, oy))

            all_idx = list(range(n_items))
            eff_exec = executor if executor is not None else getattr(self, "_prep_executor", None)
            prep_results = (
                list(eff_exec.map(_prep_pose_item, all_idx))
                if (eff_exec is not None and n_items >= 2)
                else [_prep_pose_item(i) for i in all_idx]
            )

            valid_indices: List[int] = []
            valid_chw_list: List[np.ndarray] = []
            valid_meta: List[_PoseMeta] = []
            for res in prep_results:
                if res is not None:
                    idx, chw_f32, meta = res
                    valid_indices.append(idx)
                    valid_chw_list.append(chw_f32)
                    valid_meta.append(meta)

            buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for valid_pos, chw_t in enumerate(valid_chw_list):
                buckets[(int(chw_t.shape[1]), int(chw_t.shape[2]))].append(valid_pos)

            raw_arr_by_idx: Dict[int, Tuple[np.ndarray, _PoseMeta]] = {}
            for _shape_hw, bucket_positions in buckets.items():
                bucket_tensors = [valid_chw_list[p] for p in bucket_positions]
                batch_records = _run_onnx_minibatches_back_to_back(active_sess, input_name, bucket_tensors, sub_bs)
                for start_i, end_i, raw_outs in batch_records:
                    out0 = np.asarray(raw_outs[0], dtype=np.float32)
                    for b_off, valid_pos in enumerate(bucket_positions[start_i:end_i]):
                        orig_idx = valid_indices[valid_pos]
                        arr = out0[b_off] if out0.ndim >= 3 else out0
                        raw_arr_by_idx[orig_idx] = (arr, valid_meta[valid_pos])

            def _decode_pose_item(idx: int) -> LM3PoseResult:
                lm_i = lamina_masks[idx] if (lamina_masks and idx < len(lamina_masks)) else None
                pm_i = petiole_masks[idx] if (petiole_masks and idx < len(petiole_masks)) else None
                ox, oy = _offset_at(idx)
                if idx in raw_arr_by_idx:
                    arr, meta_i = raw_arr_by_idx[idx]
                    decoded = self._decode_pose_array(arr, meta_i, lamina_mask=lm_i, petiole_mask=pm_i)
                    if decoded is not None:
                        return decoded
                return self._geometric_fallback_keypoints_31(
                    crop_bgr=crops_bgr[idx], lamina_mask=lm_i, petiole_mask=pm_i, offset_xy=(ox, oy),
                )

            if eff_exec is not None and n_items >= 2:
                return list(eff_exec.map(_decode_pose_item, all_idx))
            return [_decode_pose_item(i) for i in all_idx]
        except Exception as exc:
            logger.debug(f"LM3PoseLandmarkEngine batched inference fallback: {exc}")
            return self._predict_batch_sequential(crops_bgr, lamina_masks, petiole_masks, crop_offsets)

    def predict(
        self, crop_bgr: np.ndarray, lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None, crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        """Alias for `predict_landmarks`."""
        return self.predict_landmarks(
            crop_bgr=crop_bgr, lamina_mask=lamina_mask, petiole_mask=petiole_mask, crop_offset=crop_offset,
        )

    def __call__(
        self, crop_bgr: np.ndarray, lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None, crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        return self.predict_landmarks(
            crop_bgr=crop_bgr, lamina_mask=lamina_mask, petiole_mask=petiole_mask, crop_offset=crop_offset,
        )
