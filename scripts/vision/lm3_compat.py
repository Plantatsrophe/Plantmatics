#!/usr/bin/env python3
"""
===============================================================================
Script: lm3_compat.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline

Description:
    Optional dependency compatibility layer providing pure NumPy / OpenCV
    fallbacks for `sklearn.cluster.DBSCAN` and `leafmachine3` geometric helpers
    when running outside the dedicated `.venv_LM3` environment (e.g., unit tests,
    Step 01 harvesting, or lightweight CPU QA environments).
    When `leafmachine3` IS installed in the active interpreter, all calls
    delegate 100% unchanged to the canonical `leafmachine3` implementation.
===============================================================================
"""

from __future__ import annotations

import importlib.util
import sys
import types
from typing import Any

import numpy as np

from scripts.vision.step02_lm3_compat_ect_metrics import (
    _fb_compute_ect,
    _fb_compute_measurements,
    _fb_ect_product_for,
    _fb_measure_petiole,
    _fb_render_leaf_products,
    _fb_render_radial_ect_overlay,
    _fb_tip_up_directions,
    _fb_visual_log,
)
from scripts.vision.step02_lm3_compat_geometry import (
    _FbOrientedLeaf,
    _fb_apply_affine,
    _fb_components,
    _fb_compute_orientation,
    _fb_crop_to_box,
    _fb_extract_contour,
    _fb_largest_component,
    _fb_length_width_from_box,
    _fb_mask_bbox,
    _fb_measure_leaf,
    _fb_oriented_affine,
    _fb_outer_polygons,
    _fb_polygon_mask,
    _fb_polygon_morphology,
    _fb_prepare_mask,
    _fb_resample_closed,
    _fb_rotate_image,
    _fb_signed_area,
)

try:
    from sklearn.cluster import DBSCAN as _SklearnDBSCAN  # type: ignore[import-untyped]
except ImportError:
    _SklearnDBSCAN = None


class _NumpyDBSCAN:
    """Exact 2D/ND Euclidean DBSCAN fallback matching `sklearn.cluster.DBSCAN(eps, min_samples).fit(X)`."""

    def __init__(self, eps: float = 0.5, min_samples: int = 5, **_kwargs: Any) -> None:
        self.eps = float(eps)
        self.min_samples = int(min_samples)
        self.labels_: np.ndarray = np.empty(0, dtype=np.int64)

    def fit(self, X: Any) -> "_NumpyDBSCAN":
        arr = np.asarray(X, dtype=np.float64)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        n = int(arr.shape[0])
        if n == 0:
            self.labels_ = np.empty(0, dtype=np.int64)
            return self

        diff = arr[:, None, :] - arr[None, :, :]
        dist_mat = np.sqrt(np.sum(diff * diff, axis=-1))
        neighbors = [np.where(dist_mat[i] <= self.eps)[0] for i in range(n)]
        is_core = np.array([len(nb) >= self.min_samples for nb in neighbors], dtype=bool)

        labels = np.full(n, -1, dtype=np.int64)
        visited = np.zeros(n, dtype=bool)
        cluster_id = 0

        for i in range(n):
            if visited[i]:
                continue
            visited[i] = True
            if not is_core[i]:
                continue
            labels[i] = cluster_id
            queue = list(int(x) for x in neighbors[i] if int(x) != i)
            queued_set = set(queue)
            idx = 0
            while idx < len(queue):
                q = queue[idx]
                idx += 1
                if not visited[q]:
                    visited[q] = True
                    if is_core[q]:
                        for nb_pt in neighbors[q]:
                            nb_int = int(nb_pt)
                            if nb_int not in queued_set and labels[nb_int] == -1:
                                queue.append(nb_int)
                                queued_set.add(nb_int)
                if labels[q] == -1:
                    labels[q] = cluster_id
            cluster_id += 1

        self.labels_ = labels
        return self


DBSCAN = _SklearnDBSCAN if _SklearnDBSCAN is not None else _NumpyDBSCAN


def _register_virtual_leafmachine3_modules() -> None:
    """Register fallback `leafmachine3.*` submodules in `sys.modules` ONLY when `leafmachine3` is absent."""
    try:
        if importlib.util.find_spec("leafmachine3") is not None:
            return
    except Exception:
        pass
    if "leafmachine3" in sys.modules:
        return

    def _mk_pkg(name: str) -> types.ModuleType:
        m = types.ModuleType(name)
        m.__path__ = []  # type: ignore[attr-defined]
        sys.modules[name] = m
        return m

    def _mk_mod(parent: types.ModuleType, attr: str, full_name: str, **attrs: Any) -> types.ModuleType:
        m = types.ModuleType(full_name)
        for k, v in attrs.items():
            setattr(m, k, v)
        sys.modules[full_name] = m
        setattr(parent, attr, m)
        return m

    root_pkg = _mk_pkg("leafmachine3")
    core_pkg = _mk_pkg("leafmachine3.core")
    inf_pkg = _mk_pkg("leafmachine3.inference")
    rep_pkg = _mk_pkg("leafmachine3.reporting")
    root_pkg.core = core_pkg  # type: ignore[attr-defined]
    root_pkg.inference = inf_pkg  # type: ignore[attr-defined]
    root_pkg.reporting = rep_pkg  # type: ignore[attr-defined]

    _mk_mod(core_pkg, "leaf_outline", "leafmachine3.core.leaf_outline", extract_contour=_fb_extract_contour)
    _mk_mod(
        core_pkg,
        "momocs_prep",
        "leafmachine3.core.momocs_prep",
        prepare_mask=_fb_prepare_mask,
        resample_closed=_fb_resample_closed,
        signed_area=_fb_signed_area,
    )
    _mk_mod(core_pkg, "morphometrics", "leafmachine3.core.morphometrics", polygon_morphology=_fb_polygon_morphology)
    _mk_mod(
        core_pkg,
        "orientation",
        "leafmachine3.core.orientation",
        compute_orientation=_fb_compute_orientation,
        length_width_from_box=_fb_length_width_from_box,
    )
    _mk_mod(
        core_pkg,
        "imaging",
        "leafmachine3.core.imaging",
        polygon_mask=_fb_polygon_mask,
        mask_bbox=_fb_mask_bbox,
        crop_to_box=_fb_crop_to_box,
        rotate_image=_fb_rotate_image,
    )
    _mk_mod(
        core_pkg,
        "bilateral",
        "leafmachine3.core.bilateral",
        OrientedLeaf=_FbOrientedLeaf,
        oriented_affine=_fb_oriented_affine,
        apply_affine=_fb_apply_affine,
        measure_leaf=_fb_measure_leaf,
    )
    _mk_mod(
        core_pkg,
        "ect_compute",
        "leafmachine3.core.ect_compute",
        ect_product_for=_fb_ect_product_for,
        compute_ect=_fb_compute_ect,
        tip_up_directions=_fb_tip_up_directions,
    )
    _mk_mod(core_pkg, "petiole", "leafmachine3.core.petiole", measure_petiole=_fb_measure_petiole)
    _mk_mod(
        core_pkg,
        "landmark_metrics",
        "leafmachine3.core.landmark_metrics",
        compute_measurements=_fb_compute_measurements,
    )
    _mk_mod(
        inf_pkg,
        "leaf_segmenter",
        "leafmachine3.inference.leaf_segmenter",
        _components=_fb_components,
        _largest_component=_fb_largest_component,
        _outer_polygons=_fb_outer_polygons,
    )
    _mk_mod(
        rep_pkg,
        "leaf_products",
        "leafmachine3.reporting.leaf_products",
        render_leaf_products=_fb_render_leaf_products,
    )
    _mk_mod(
        rep_pkg,
        "ect_viz",
        "leafmachine3.reporting.ect_viz",
        visual_log=_fb_visual_log,
        render_radial_ect_overlay=_fb_render_radial_ect_overlay,
    )


_register_virtual_leafmachine3_modules()

from leafmachine3.core.bilateral import OrientedLeaf, apply_affine, measure_leaf, oriented_affine  # noqa: E402
from leafmachine3.core.ect_compute import compute_ect, ect_product_for, tip_up_directions  # noqa: E402
from leafmachine3.core.imaging import crop_to_box, mask_bbox, polygon_mask, rotate_image  # noqa: E402
from leafmachine3.core.landmark_metrics import compute_measurements  # noqa: E402
from leafmachine3.core.leaf_outline import extract_contour  # noqa: E402
from leafmachine3.core.momocs_prep import (  # noqa: E402
    resample_closed,
    signed_area as momocs_signed_area,
)
from leafmachine3.core.morphometrics import polygon_morphology  # noqa: E402
from leafmachine3.core.orientation import compute_orientation, length_width_from_box  # noqa: E402
from leafmachine3.core.petiole import measure_petiole  # noqa: E402
from leafmachine3.inference.leaf_segmenter import (  # noqa: E402
    _components,
    _largest_component,
    _outer_polygons,
)
from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log  # noqa: E402
from leafmachine3.reporting.leaf_products import render_leaf_products  # noqa: E402
