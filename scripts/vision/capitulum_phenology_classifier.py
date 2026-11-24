#!/usr/bin/env python3
"""
Module: capitulum_phenology_classifier.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline (NCU)

Vectorized "Crop-and-Classify" visual biomarker classifier for Asteraceae capitula
(LeafMachine3 reproductive Classes 2..6: Seed_Fruit_ONE, Seed_Fruit_MANY,
Flower_ONE, Flower_MANY, and Bud). Resolves HEAD_BUD, HEAD_ANTHESIS, and
HEAD_FRUIT at head level and aggregates voucher-level states via the Botanical
Determinate Cyme Rule (anthesis, fruit, sterile, bud).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from scripts.vision.step02_capitulum_biomarkers_early import (
    extract_biomarker_involucre_posture,
    extract_biomarker_pappus_plume,
    extract_biomarker_yellow_corollas,
    orient_capitulum_crop_upright,
)
from scripts.vision.step02_capitulum_biomarkers_late import (
    extract_biomarker_aged_ochre_corollas,
    extract_biomarker_radiate_ligules,
)

logger = logging.getLogger(__name__)

# State Constants
HEAD_ANTHESIS = "anthesis"
HEAD_FRUIT = "fruit"
HEAD_BUD = "bud"
HEAD_UNKNOWN = "unknown"

VOUCHER_ANTHESIS = "anthesis"
VOUCHER_FRUIT = "fruit"
VOUCHER_BUD = "bud"
VOUCHER_STERILE = "sterile"
VOUCHER_UNKNOWN = "unknown"


class CapitulumResult(dict):
    """Specialized dictionary returning head metrics while supporting string equality against head_state."""

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.get("head_state") == other
        return super().__eq__(other)

    def __hash__(self) -> int:
        return hash(self.get("head_state"))


HeadPhenologyResult = CapitulumResult


class VoucherPhenologyResult(dict):
    """Specialized dictionary returning voucher metrics while supporting string equality against voucher_phenological_state."""

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.get("voucher_phenological_state") == other
        return super().__eq__(other)

    def __hash__(self) -> int:
        return hash(self.get("voucher_phenological_state"))


def classify_single_capitulum(
    crop_rgb: np.ndarray,
    min_crop_dim: int = 5,
    min_pappus_luminance: float = 75.0,
    max_pappus_saturation: float = 0.18,
    min_pappus_area_ratio: float = 0.18,
    min_floret_yellow_saturation: float = 0.32,
    min_floret_yellow_ratio: float = 0.12,
    class_id: Optional[int] = None,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
) -> CapitulumResult:
    """
    Performs head-level phenological biomarker classification on a single cropped capitulum.

    Decision Logic:
      1. has_pappus_plume and not has_fresh_yellow_corollas -> HEAD_FRUIT.
      2. has_fresh_yellow_corollas and erect_involucre -> HEAD_ANTHESIS.
      3. Neither pappus nor fresh yellow, aged ochre passing class/posture gating -> HEAD_ANTHESIS.
      4. erect_involucre and low yellow and low white -> HEAD_BUD.
      5. Fallback: pappus -> HEAD_FRUIT, fresh yellow -> HEAD_ANTHESIS, reflexed -> HEAD_FRUIT, else HEAD_BUD.
    """
    if (
        crop_rgb is None
        or crop_rgb.size == 0
        or crop_rgb.ndim < 3
        or crop_rgb.shape[2] < 3
        or crop_rgb.shape[0] < min_crop_dim
        or crop_rgb.shape[1] < min_crop_dim
    ):
        return CapitulumResult({
            "head_state": HEAD_UNKNOWN,
            "has_pappus_plume": False,
            "has_fresh_yellow_corollas": False,
            "has_aged_ochre_corollas": False,
            "has_radiate_ligules": False,
            "erect_involucre": False,
            "reflexed_involucre": False,
            "white_fibrous_ratio": 0.0,
            "yellow_ratio": 0.0,
            "ochre_ratio": 0.0,
            "involucre_aspect_ratio": 0.0,
            "boundary_roughness": 1.0,
            "diagnostics": {},
        })

    oriented_rgb = orient_capitulum_crop_upright(crop_rgb, peduncle_dir_xy=peduncle_dir_xy)

    # Extract Biomarker A: Pappus Plume
    has_pappus_plume, white_fibrous_ratio, diag_a = extract_biomarker_pappus_plume(
        oriented_rgb,
        min_l_star=float(min_pappus_luminance),
        max_s=float(max_pappus_saturation),
        plume_threshold=float(min_pappus_area_ratio),
    )

    # Extract Biomarker B: Active Yellow Corollas
    has_fresh_yellow_corollas, yellow_ratio, diag_b = extract_biomarker_yellow_corollas(
        oriented_rgb,
        s_min=float(min_floret_yellow_saturation),
        corolla_threshold=float(min_floret_yellow_ratio),
    )

    # Extract Biomarker B2: Aged Ochre / Amber Corollas
    has_aged_ochre, ochre_ratio, diag_ochre = extract_biomarker_aged_ochre_corollas(
        oriented_rgb,
        s_max=float(min_floret_yellow_saturation),
    )

    # Extract Biomarker C: Involucre Posture
    erect_involucre, reflexed_involucre, diag_c = extract_biomarker_involucre_posture(oriented_rgb)
    inv_ar = float(diag_c.get("involucre_aspect_ratio", 1.0))

    aged_anthesis_ok = bool(
        (not has_pappus_plume)
        and (not has_fresh_yellow_corollas)
        and has_aged_ochre
        and (class_id != 6)
        and (
            (class_id in (4, 5))
            or (erect_involucre and inv_ar >= 0.95 and (not reflexed_involucre))
        )
    )

    # Head Decision Logic
    if has_pappus_plume and (not has_fresh_yellow_corollas):
        head_state = HEAD_FRUIT
    elif has_fresh_yellow_corollas and erect_involucre:
        head_state = HEAD_ANTHESIS
    elif aged_anthesis_ok:
        head_state = HEAD_ANTHESIS
    elif erect_involucre and (not has_fresh_yellow_corollas) and (not has_pappus_plume):
        head_state = HEAD_BUD
    else:
        # Resolving intermediate/ambiguous boundary cases
        if has_pappus_plume:
            head_state = HEAD_FRUIT
        elif has_fresh_yellow_corollas:
            head_state = HEAD_ANTHESIS
        elif reflexed_involucre:
            head_state = HEAD_FRUIT
        else:
            head_state = HEAD_BUD

    eval_radiate = bool(
        (class_id in (4, 5))
        or (class_id is None and (head_state == HEAD_ANTHESIS or has_fresh_yellow_corollas or aged_anthesis_ok))
    )
    if eval_radiate:
        has_radiate_ligules, diag_rad = extract_biomarker_radiate_ligules(
            oriented_rgb,
            min_floret_yellow_saturation=float(min_floret_yellow_saturation),
        )
    else:
        has_radiate_ligules = False
        diag_rad = {"outer_corolla_frac": 0.0, "corolla_to_cup_width_ratio": 0.0}

    return CapitulumResult({
        "head_state": head_state,
        "has_pappus_plume": has_pappus_plume,
        "has_fresh_yellow_corollas": has_fresh_yellow_corollas,
        "has_aged_ochre_corollas": bool(aged_anthesis_ok),
        "has_radiate_ligules": bool(has_radiate_ligules),
        "erect_involucre": erect_involucre,
        "reflexed_involucre": reflexed_involucre,
        "white_fibrous_ratio": white_fibrous_ratio,
        "yellow_ratio": yellow_ratio,
        "ochre_ratio": ochre_ratio,
        "involucre_aspect_ratio": diag_c.get("involucre_aspect_ratio", 1.0),
        "boundary_roughness": diag_c.get("boundary_roughness", 1.0),
        "diagnostics": {
            "pappus": diag_a,
            "corolla": diag_b,
            "aged_ochre": diag_ochre,
            "involucre": diag_c,
            "radiate": diag_rad,
        },
    })


def classify_voucher_phenology(
    catalog_number: str,
    head_classifications: List[Dict[str, Any]],
    min_anthesis_heads_for_flowering: int = 1,
) -> Dict[str, Any]:
    """
    Aggregates head-level classifications for a voucher into a single phenological state
    governed by the Botanical Determinate Cyme Rule with predominantly-fruiting corymb guarding:
      - N_anthesis >= 1 (with multi-head / proportion / fresh-yellow guard when N_fruit >= 3) -> "anthesis"
      - N_fruit > 0 and N_anthesis == 0 (or failing predominantly-fruiting guard) -> "fruit"
      - N_total == 0 -> "sterile"
      - N_bud > 0 and N_anthesis == 0 and N_fruit == 0 -> "bud"
    """
    total_capitula = len(head_classifications)
    n_anthesis = sum(1 for h in head_classifications if h.get("head_state") == HEAD_ANTHESIS)
    n_fruiting = sum(1 for h in head_classifications if h.get("head_state") == HEAD_FRUIT)
    n_bud = sum(1 for h in head_classifications if h.get("head_state") == HEAD_BUD)

    min_anth = max(1, int(min_anthesis_heads_for_flowering))
    has_strong_fresh_yellow = any(
        h.get("head_state") == HEAD_ANTHESIS
        and bool(h.get("has_fresh_yellow_corollas", False))
        and float(h.get("yellow_ratio", 0.0) or 0.0) >= 0.18
        for h in head_classifications
    )

    if n_anthesis >= min_anth:
        if n_fruiting >= 3:
            anth_frac = float(n_anthesis) / float(max(1, n_anthesis + n_fruiting))
            if (
                n_anthesis >= max(2, min_anth)
                or anth_frac >= 0.25
                or has_strong_fresh_yellow
            ):
                voucher_state = VOUCHER_ANTHESIS
            else:
                voucher_state = VOUCHER_FRUIT
        else:
            voucher_state = VOUCHER_ANTHESIS
    elif n_fruiting > 0 and n_anthesis == 0:
        voucher_state = VOUCHER_FRUIT
    elif total_capitula == 0:
        voucher_state = VOUCHER_STERILE
    elif n_bud > 0 and n_anthesis == 0 and n_fruiting == 0:
        voucher_state = VOUCHER_BUD
    else:
        voucher_state = VOUCHER_UNKNOWN

    pappus_ratios = [
        float(h.get("white_fibrous_ratio", 0.0))
        for h in head_classifications
        if "white_fibrous_ratio" in h
    ]
    dominant_pappus_ratio = round(max(pappus_ratios), 4) if pappus_ratios else 0.0

    flowering_eval_heads = [
        h
        for h in head_classifications
        if h.get("head_state") != HEAD_UNKNOWN
        and (
            h.get("head_state") == HEAD_ANTHESIS
            or (h.get("class_id") is not None and int(h.get("class_id")) in (4, 5))
        )
    ]
    if flowering_eval_heads:
        n_radiate = sum(1 for h in flowering_eval_heads if bool(h.get("has_radiate_ligules", False)))
        radiate_head_fraction = round(float(n_radiate) / float(len(flowering_eval_heads)), 4)
    else:
        radiate_head_fraction = float("nan")

    return VoucherPhenologyResult({
        "catalogNumber": str(catalog_number).strip(),
        "total_capitula": int(total_capitula),
        "n_anthesis": int(n_anthesis),
        "n_fruiting": int(n_fruiting),
        "n_bud": int(n_bud),
        "voucher_phenological_state": voucher_state,
        "dominant_pappus_ratio": float(dominant_pappus_ratio),
        "radiate_head_fraction": float(radiate_head_fraction),
    })
