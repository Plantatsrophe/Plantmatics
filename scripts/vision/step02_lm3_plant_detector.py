#!/usr/bin/env python3
"""
Step 02 LeafMachine3 Stage 1 Plant Organ Detector (YOLO26x end-to-end NMS-free, 1280 px).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np

from scripts.vision.lm3_data_loader import (
    LM3_LEAF_CLASS_IDS,
    LM3_PLANT_CLASSES,
    LM3_REPRODUCTIVE_CLASS_IDS,
    compute_work_scale,
    load_precomputed_plant_labels,
    unscale_working_bbox_to_sheet,
)
from scripts.vision.step02_lm3_onnx_runtime import (
    _bgr_to_nchw_float32,
    _letterbox_bgr,
    _letterbox_gain_pad_hw,
    _prepare_working_image,
    _run_onnx_minibatches_back_to_back,
)

logger = logging.getLogger("SegmentAndExtract")

_DetTuple = Tuple[Tuple[int, int, int, int], float, int]
_S1Meta = Tuple[float, Tuple[float, float], Tuple[int, int], float, Tuple[int, int]]


class LM3PlantDetector:
    """
    LeafMachine3 Stage 1 ONNX Plant Organ Detector (`lm3_plant_detector__yolo26x_det_1280`).
    Runs 1280 px end-to-end NMS-free YOLO26x detection across 9 canonical plant classes:
      0: Leaf_WHOLE, 1: Leaf_PARTIAL, 2: Seed_Fruit_ONE, 3: Seed_Fruit_MANY,
      4: Flower_ONE, 5: Flower_MANY, 6: Bud, 7: Roots, 8: Wood.
    Scales detections from the working frame (`max_working_dim=3200`) back to native full-sheet
    coordinates via `work_scale` (`unscale_working_bbox_to_sheet`).
    """

    MODEL_NAME: str = "lm3_plant_detector__yolo26x_det_1280"
    CLASSES: Dict[int, str] = LM3_PLANT_CLASSES
    LEAF_CLASS_IDS: Tuple[int, ...] = LM3_LEAF_CLASS_IDS
    REPRODUCTIVE_CLASS_IDS: Tuple[int, ...] = LM3_REPRODUCTIVE_CLASS_IDS

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.30,
        iou_thres: float = 0.95,
        img_size: int = 1280,
        max_working_dim: int = 3200,
        lm3_project_dir: Optional[Union[str, Path]] = None,
        providers: Optional[Sequence[str]] = None,
        **kwargs: Any,
    ) -> None:
        import scripts.vision.lm3_inference_engines as _lm3_eng

        self.device = device if (device == "cuda" and _lm3_eng.torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.iou_thres = float(iou_thres)
        self.img_size = int(img_size)
        self.max_working_dim = int(max_working_dim)
        proj_dir = lm3_project_dir or kwargs.get("project_dir")
        self.lm3_project_dir = Path(proj_dir) if proj_dir else None
        self.stride = 32
        self.last_work_scale: float = 1.0
        self.weights_path, self.model_name, self.providers, self.session = _lm3_eng._load_onnx_session(
            weights_path=weights_path, default_model_name=self.MODEL_NAME, device=self.device, providers=providers,
        )
        self.model = self.session

    def prepare_letterbox(
        self,
        image_bgr: np.ndarray,
        max_working_dim: Optional[int] = None,
    ) -> Optional[np.ndarray]:
        """Precomputes the working-frame letterboxed BGR array on a background CPU worker thread."""
        if (self.session is None and self.model is None) or image_bgr is None or image_bgr.size == 0:
            return None
        work_dim = int(max_working_dim) if max_working_dim is not None else self.max_working_dim
        work_img, _ = _prepare_working_image(image_bgr, max_working_dim=work_dim)
        lb_bgr, _, _ = _letterbox_bgr(work_img, self.img_size, auto=True, stride=32)
        return lb_bgr

    _preprocess_letterbox = prepare_letterbox

    def _resolve_letterbox(
        self,
        image_bgr: np.ndarray,
        work_dim: int,
        precomputed_letterbox: Optional[np.ndarray],
    ) -> Tuple[np.ndarray, float, Tuple[float, float], Tuple[int, int], float]:
        sh, sw = image_bgr.shape[:2]
        if (
            precomputed_letterbox is not None
            and precomputed_letterbox.ndim == 3
            and max(precomputed_letterbox.shape[:2]) == self.img_size
            and precomputed_letterbox.shape[0] % 32 == 0
            and precomputed_letterbox.shape[1] % 32 == 0
        ):
            w_scale = compute_work_scale((sh, sw), max_working_dim=work_dim)
            if w_scale < 1.0:
                work_w, work_h = max(1, int(round(sw * w_scale))), max(1, int(round(sh * w_scale)))
            else:
                work_h, work_w = sh, sw
            lb_bgr = precomputed_letterbox
            r, pad_xy = _letterbox_gain_pad_hw(lb_bgr.shape[:2], (work_h, work_w))
            return lb_bgr, r, pad_xy, (work_h, work_w), w_scale
        work_img, w_scale = _prepare_working_image(image_bgr, max_working_dim=work_dim)
        work_h, work_w = work_img.shape[:2]
        lb_bgr, r, pad_xy = _letterbox_bgr(work_img, self.img_size, auto=True, stride=32)
        return lb_bgr, r, pad_xy, (work_h, work_w), w_scale

    def _decode_yolo26x_det_output(
        self,
        raw_output: np.ndarray,
        ratio: float,
        pad_xy: Tuple[float, float],
        work_shape: Tuple[int, int],
        work_scale: float,
        sheet_shape: Tuple[int, int],
    ) -> List[_DetTuple]:
        """
        Decodes end-to-end NMS-free YOLO26x detections `(B, N, 6)` `[x1, y1, x2, y2, score, cls]`
        or anchor-grid `(B, 4+C, A)` outputs, unletterboxing to the working frame and unscaling
        to native sheet coordinates via `unscale_working_bbox_to_sheet`.
        """
        arr = np.asarray(raw_output, dtype=np.float32)
        if arr.ndim == 3:
            arr = arr[0]
        if arr.size == 0:
            return []

        dw, dh = pad_xy
        work_h, work_w = int(work_shape[0]), int(work_shape[1])
        sheet_h, sheet_w = int(sheet_shape[0]), int(sheet_shape[1])
        inv_r = 1.0 / max(float(ratio), 1e-6)
        results: List[_DetTuple] = []

        def _to_sheet_bbox(x1_l: float, y1_l: float, x2_l: float, y2_l: float) -> Optional[Tuple[int, int, int, int]]:
            x1_w = max(0.0, min(float(work_w), (x1_l - dw) * inv_r))
            y1_w = max(0.0, min(float(work_h), (y1_l - dh) * inv_r))
            x2_w = max(0.0, min(float(work_w), (x2_l - dw) * inv_r))
            y2_w = max(0.0, min(float(work_h), (y2_l - dh) * inv_r))
            if y2_w <= y1_w or x2_w <= x1_w:
                return None
            s_bbox = unscale_working_bbox_to_sheet((y1_w, x1_w, y2_w, x2_w), work_scale, sheet_shape=(sheet_h, sheet_w))
            return s_bbox if (s_bbox[2] > s_bbox[0] and s_bbox[3] > s_bbox[1]) else None

        if arr.ndim == 2 and arr.shape[1] in (6, 7):
            # End-to-end NMS-free YOLO26x format: [x1, y1, x2, y2, score, class_id]
            for row in arr:
                score = float(row[4])
                if score < self.conf_thres:
                    continue
                cls_id = int(round(float(row[5])))
                if cls_id not in self.CLASSES:
                    continue
                sheet_bbox = _to_sheet_bbox(float(row[0]), float(row[1]), float(row[2]), float(row[3]))
                if sheet_bbox is not None:
                    results.append((sheet_bbox, score, cls_id))
            return self._deduplicate_sheet_boxes(results)

        # Transpose (4+C, A) -> (A, 4+C) if needed
        if arr.ndim == 2 and arr.shape[0] < arr.shape[1] and arr.shape[0] >= 5:
            arr = arr.T
        if arr.ndim != 2 or arr.shape[1] < 5:
            return []

        n_cls = len(self.CLASSES)
        boxes_xywh = arr[:, :4]
        cls_scores = arr[:, 4 : 4 + n_cls]
        if cls_scores.size == 0:
            return []
        best_cls = np.argmax(cls_scores, axis=1)
        best_scores = cls_scores[np.arange(len(arr)), best_cls]
        keep = best_scores >= self.conf_thres
        if not np.any(keep):
            return []

        boxes_xywh = boxes_xywh[keep]
        best_scores = best_scores[keep]
        best_cls = best_cls[keep]

        xyxy_boxes: List[List[float]] = [
            [float(xc - 0.5 * bw), float(yc - 0.5 * bh), float(bw), float(bh)]
            for (xc, yc, bw, bh) in boxes_xywh
        ]
        indices = cv2.dnn.NMSBoxes(
            bboxes=xyxy_boxes,
            scores=[float(s) for s in best_scores],
            score_threshold=self.conf_thres,
            nms_threshold=self.iou_thres,
        )
        if len(indices) == 0:
            return []

        for idx_raw in np.asarray(indices).reshape(-1):
            idx = int(idx_raw)
            xc, yc, bw, bh = boxes_xywh[idx]
            sheet_bbox = _to_sheet_bbox(
                float(xc - 0.5 * bw), float(yc - 0.5 * bh), float(xc + 0.5 * bw), float(yc + 0.5 * bh),
            )
            cls_id = int(best_cls[idx])
            if sheet_bbox is not None and cls_id in self.CLASSES:
                results.append((sheet_bbox, float(best_scores[idx]), cls_id))
        return self._deduplicate_sheet_boxes(results)

    _postprocess_yolo26x_det = _decode_yolo26x_det_output

    def _deduplicate_sheet_boxes(
        self,
        dets: List[_DetTuple],
        containment_thres: float = 0.95,
    ) -> List[_DetTuple]:
        """
        Applies LeafMachine3 native same-class duplicate box suppression
        (`suppress_duplicate_boxes(..., overlap=0.95, metric='min')`), sorting candidates in
        descending confidence order and suppressing a box only when
        `inter / min_area >= containment_thres` with a higher-confidence kept box.
        """
        if len(dets) <= 1:
            return dets
        eff_overlap = (
            max(float(containment_thres), float(self.iou_thres))
            if self.iou_thres > 0.85
            else float(containment_thres)
        )
        by_class: Dict[int, List[Tuple[int, Tuple[int, int, int, int], float, int]]] = defaultdict(list)
        for orig_idx, (bbox, score, cls_id) in enumerate(dets):
            by_class[int(cls_id)].append((orig_idx, bbox, float(score), int(cls_id)))

        kept_with_idx: List[Tuple[int, Tuple[int, int, int, int], float, int]] = []
        for _cid, group in by_class.items():
            group.sort(key=lambda item: item[2], reverse=True)
            accepted: List[Tuple[int, Tuple[int, int, int, int], float, int]] = []
            for cand in group:
                _, (y1_a, x1_a, y2_a, x2_a), _sc_a, _c_a = cand
                area_a = float(max(1, (y2_a - y1_a) * (x2_a - x1_a)))
                is_dup = False
                for _, (y1_b, x1_b, y2_b, x2_b), _sc_b, _c_b in accepted:
                    iy1, ix1 = max(y1_a, y1_b), max(x1_a, x1_b)
                    iy2, ix2 = min(y2_a, y2_b), min(x2_a, x2_b)
                    if iy2 <= iy1 or ix2 <= ix1:
                        continue
                    inter = float((iy2 - iy1) * (ix2 - ix1))
                    area_b = float(max(1, (y2_b - y1_b) * (x2_b - x1_b)))
                    if (inter / min(area_a, area_b)) >= eff_overlap:
                        is_dup = True
                        break
                if not is_dup:
                    accepted.append(cand)
            kept_with_idx.extend(accepted)

        kept_with_idx.sort(key=lambda item: item[0])
        return [(bbox, score, cls_id) for _, bbox, score, cls_id in kept_with_idx]

    def detect(
        self,
        image_bgr: np.ndarray,
        max_working_dim: int = 3200,
        catalog_number: Optional[str] = None,
        use_cache: bool = False,
        precomputed_letterbox: Optional[np.ndarray] = None,
    ) -> List[_DetTuple]:
        """
        Runs LeafMachine3 9-class full-sheet organ detection on `image_bgr`.
        Scales predictions from `max_working_dim` working coordinates back to native
        sheet dimensions `(sheet_h, sheet_w)` via `work_scale`.
        Returns list of `((ymin, xmin, ymax, xmax), score, class_id)`.
        """
        import scripts.vision.lm3_inference_engines as _lm3_eng

        if image_bgr is None or image_bgr.size == 0:
            return []
        sheet_h, sheet_w = image_bgr.shape[:2]
        work_dim = int(max_working_dim) if max_working_dim > 0 else self.max_working_dim
        work_scale = compute_work_scale((sheet_h, sheet_w), max_working_dim=work_dim)
        self.last_work_scale = work_scale

        active_sess = self.session if self.session is not None else self.model
        if use_cache and active_sess is None and self.lm3_project_dir and catalog_number:
            label_candidates = [
                self.lm3_project_dir / "Data" / "precomputed_annotations" / "Plant_Components" / "labels" / f"{catalog_number}.txt",
                self.lm3_project_dir / "Plant_Components" / "labels" / f"{catalog_number}.txt",
            ]
            for lf in label_candidates:
                if lf.exists():
                    dets = load_precomputed_plant_labels(lf, sheet_w, sheet_h, min_confidence=self.conf_thres)
                    if dets:
                        return dets

        if active_sess is None:
            return []

        try:
            lb_bgr, r, pad_xy, work_hw, work_scale = self._resolve_letterbox(
                image_bgr, work_dim, precomputed_letterbox
            )
            self.last_work_scale = work_scale
            inp_tensor = _bgr_to_nchw_float32(lb_bgr)
            input_name = active_sess.get_inputs()[0].name
            with _lm3_eng._ONNX_GPU_LOCK:
                raw_outs = active_sess.run(None, {input_name: inp_tensor})
            if not raw_outs:
                return []
            return self._decode_yolo26x_det_output(
                raw_outs[0], ratio=r, pad_xy=pad_xy, work_shape=work_hw,
                work_scale=work_scale, sheet_shape=(sheet_h, sheet_w),
            )
        except Exception as exc:
            logger.debug(f"LM3PlantDetector ONNX inference error: {exc}")
            return []

    def predict(
        self,
        image_bgr: np.ndarray,
        catalog_number: Optional[str] = None,
        use_cache: bool = False,
        precomputed_letterbox: Optional[np.ndarray] = None,
        max_working_dim: Optional[int] = None,
    ) -> List[_DetTuple]:
        """Alias for `.detect(...)` preserving Step 02 runner and unit test compatibility."""
        eff_dim = int(max_working_dim) if max_working_dim is not None else self.max_working_dim
        if "detect" in self.__dict__:
            return self.detect(
                image_bgr, max_working_dim=eff_dim, catalog_number=catalog_number,
                use_cache=use_cache, precomputed_letterbox=precomputed_letterbox,
            )
        return LM3PlantDetector.detect(
            self, image_bgr, max_working_dim=eff_dim, catalog_number=catalog_number,
            use_cache=use_cache, precomputed_letterbox=precomputed_letterbox,
        )

    def __call__(self, image_bgr: np.ndarray, max_working_dim: int = 3200, **kwargs: Any) -> List[_DetTuple]:
        if "predict" in self.__dict__ and "detect" not in self.__dict__:
            return self.predict(image_bgr, max_working_dim=max_working_dim, **kwargs)
        return self.detect(image_bgr, max_working_dim=max_working_dim, **kwargs)

    def predict_batch(
        self,
        items: Sequence[Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]],
        use_cache: bool = False,
        max_batch_size: int = 16,
        executor: Optional[Any] = None,
    ) -> List[List[_DetTuple]]:
        """Runs Stage 1 plant organ detection across a wave of herbarium sheets using shape-bucketed GPU tensors."""
        if not items:
            return []
        active_sess = self.session if self.session is not None else self.model
        if (
            active_sess is None
            or "detect" in self.__dict__
            or "predict" in self.__dict__
            or getattr(self.detect, "__func__", None) is not LM3PlantDetector.detect
            or getattr(self.predict, "__func__", None) is not LM3PlantDetector.predict
            or len(items) <= 1
        ):
            return [
                self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
                for img, cat, lb in items
            ]

        try:
            input_name = active_sess.get_inputs()[0].name
            sub_bs = min(8, max(1, int(max_batch_size)))
            all_outs: List[List[_DetTuple]] = [[] for _ in items]

            def _prep_s1_item(
                idx_item: Tuple[int, Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]]
            ) -> Optional[Tuple[int, np.ndarray, _S1Meta]]:
                idx, (img_bgr, _cat, lb_pre) = idx_item
                if img_bgr is None or img_bgr.size == 0:
                    return None
                sh, sw = img_bgr.shape[:2]
                lb_bgr, r, pad_xy, work_hw, w_scale = self._resolve_letterbox(
                    img_bgr, self.max_working_dim, lb_pre
                )
                rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
                chw_f32 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.float32)
                chw_f32 *= np.float32(1.0) / np.float32(255.0)
                return (idx, chw_f32, (r, pad_xy, work_hw, w_scale, (sh, sw)))

            indexed_items = list(enumerate(items))
            eff_exec = executor if executor is not None else getattr(self, "_prep_executor", None)
            if eff_exec is not None and len(indexed_items) >= 2:
                prep_results = list(eff_exec.map(_prep_s1_item, indexed_items))
            else:
                prep_results = [_prep_s1_item(it) for it in indexed_items]

            valid_indices: List[int] = []
            valid_chw_list: List[np.ndarray] = []
            valid_meta: List[_S1Meta] = []
            for res in prep_results:
                if res is not None:
                    idx, chw_f32, meta = res
                    self.last_work_scale = meta[3]
                    valid_indices.append(idx)
                    valid_chw_list.append(chw_f32)
                    valid_meta.append(meta)

            buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for valid_pos, chw_t in enumerate(valid_chw_list):
                buckets[(int(chw_t.shape[1]), int(chw_t.shape[2]))].append(valid_pos)

            decode_tasks: List[Tuple[int, np.ndarray, _S1Meta]] = []
            for _shape_hw, bucket_positions in buckets.items():
                bucket_tensors = [valid_chw_list[p] for p in bucket_positions]
                batch_records = _run_onnx_minibatches_back_to_back(active_sess, input_name, bucket_tensors, sub_bs)
                for start_i, end_i, raw_outs in batch_records:
                    out0 = np.asarray(raw_outs[0], dtype=np.float32)
                    for b_off, valid_pos in enumerate(bucket_positions[start_i:end_i]):
                        orig_idx = valid_indices[valid_pos]
                        slice_out = out0[b_off : b_off + 1] if out0.ndim >= 3 else out0
                        decode_tasks.append((orig_idx, slice_out, valid_meta[valid_pos]))

            def _decode_s1_item(task: Tuple[int, np.ndarray, _S1Meta]) -> Tuple[int, List[_DetTuple]]:
                orig_idx, slice_out, (r, pad_xy, work_hw, w_scale, sheet_hw) = task
                dets = self._decode_yolo26x_det_output(
                    slice_out, ratio=r, pad_xy=pad_xy, work_shape=work_hw,
                    work_scale=w_scale, sheet_shape=sheet_hw,
                )
                return orig_idx, dets

            if eff_exec is not None and len(decode_tasks) >= 2:
                decoded_pairs = list(eff_exec.map(_decode_s1_item, decode_tasks))
            else:
                decoded_pairs = [_decode_s1_item(t) for t in decode_tasks]

            for orig_idx, dets in decoded_pairs:
                all_outs[orig_idx] = dets
            return all_outs
        except Exception as exc:
            logger.debug(f"LM3PlantDetector batched inference fallback: {exc}")
            return [
                self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
                for img, cat, lb in items
            ]
