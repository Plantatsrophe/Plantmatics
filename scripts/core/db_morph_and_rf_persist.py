"""
Steps 03-05 dynamic EFA harmonic, OTU morphometric, and Spatial RF SQLite persistence.
"""

from __future__ import annotations

import os
import sqlite3
from typing import Any, Dict, List, Mapping, Sequence, Tuple

from scripts.core.db_cohort_and_env_cache import (
    _normalize_semicolon_list_str,
    upsert_environmental_cache_records,
)
from scripts.core.db_column_constants import (
    HARMONIC_COLS, PC_COLS, SYM_HARMONIC_COLS, WC_BIO_COLUMNS,
    _to_bool_int, _to_float, _to_int, _to_str, format_lat_lon_key,
)
from scripts.core.db_conn_pragmas import _ensure_dynamic_table_migrated, upsert_rows_generic
from scripts.core.db_step02_loader import ensure_pipeline_run

_LEAF_EFA_BOOL_COLS: Tuple[str, ...] = (
    "is_dissected", "is_calibrated_mm", "is_contour_outlier", "mask_is_upright", "is_type_specimen",
    "is_promoted_bronze", "is_passive", "is_abaxial_pressed_leaf",
)
_LEAF_EFA_INT_COLS: Tuple[str, ...] = (
    "efd_harmonic_opt", "efd_k95_pcs", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
)
_LEAF_EFA_STR_COLS: Tuple[str, ...] = (
    "leaf_zone", "type_status", "anchor_tier_source", "scale_source",
)
_LEAF_EFA_FLOAT_COLS: Tuple[str, ...] = (
    "solidity", "dissection_index", "aspect_ratio", "area_px", "area_mm2", "pixels_per_mm",
    "authority_weight", "spatial_balance_weight", "petiole_width_mm", "blade_width_mm",
    "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "leaf_zone_confidence",
    "dist_from_crown_norm", "cy_rel_plant_bbox", "ECT_PC1", "ECT_PC2", "ECT_PC3",
    "ect_lobe_complexity", "ect_bilateral_asym", "blade_cie_l_norm", "blade_cie_l_sd",
    "blade_cie_l_sheet_delta", "blade_cie_a_delta", "blade_cie_b_delta", "base_angle_deg",
    "apex_angle_deg", "lamina_curvature", "lamina_trace_length_px", "lamina_trace_length_mm",
    "basal_cordate_index", "margin_serration_index", "margin_tooth_acuteness",
    "terminal_lobe_area_frac", "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta",
    "abaxial_cie_l_sheet_delta", *HARMONIC_COLS, "sym_var", "asym_var", "pct_asymmetry",
    "harmonic_roughness", "leaf_quality_weight", "centroid_size", "log_centsize", *PC_COLS,
)

def persist_leaf_efa_harmonics(
    conn: sqlite3.Connection,
    run_id: str,
    efa_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
) -> int:
    """
    Persist Step 03 leaf-level Elliptic Fourier harmonics (`A1..D24`), leaf PCs
    (`PC1..PC15`), and adaptive calibration metadata
    (`efd_harmonic_opt`, `efd_k95_pcs`) into `leaf_efa_harmonics`.
    """
    if demo or not efa_rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "leaf_efa_harmonics")
    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for idx, r in enumerate(efa_rows, start=1):
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        lid = _to_int(r.get("leaf_id"), default=idx) or idx
        d: Dict[str, Any] = {
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
            "leaf_id": lid,
        }
        for col in _LEAF_EFA_STR_COLS:
            if col in r:
                d[col] = _to_str(r.get(col))
        for col in _LEAF_EFA_BOOL_COLS:
            if col in r:
                d[col] = _to_bool_int(r.get(col))
        for col in _LEAF_EFA_INT_COLS:
            if col in r:
                d[col] = _to_int(r.get(col))
        for col in _LEAF_EFA_FLOAT_COLS:
            if col in r:
                d[col] = _to_float(r.get(col))
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="leaf_efa_harmonics",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id", "leaf_id"),
        preserve_existing_on_null=False,
    )

_OTU_MORPH_BOOL_COLS: Tuple[str, ...] = (
    "is_promoted_bronze", "is_passive", "is_calibrated_mm", "is_dissected", "is_contour_outlier",
    "morphometric_outlier_flag", "misidentification_flag", "is_type_specimen",
    "has_empirical_type_anchor", "is_typical_morphotype", "is_confident_label_noise",
    "is_abaxial_pressed_leaf",
)
_OTU_MORPH_INT_COLS: Tuple[str, ...] = (
    "leaf_count", "efd_harmonic_opt", "efd_k95_pcs", "basal_leaf_count", "cauline_leaf_count",
    "self_morph_df", "open_set_df",
)
_OTU_MORPH_STR_COLS: Tuple[str, ...] = (
    "species_standardized", "determiner_tier", "type_status", "scale_source",
    "reproductive_trait_status", "morph_tier_source", "gmm_cluster", "cda_predicted_taxon",
    "cda_second_taxon", "ecoregion_cda_predicted", "triage_priority", "discordance_reason",
    "nearest_type_morph_taxon", "nearest_type_basionym", "min_congener_taxon",
)
_OTU_MORPH_FLOAT_COLS: Tuple[str, ...] = (
    "authority_weight", "spatial_balance_weight", "pixels_per_mm", "area_mm2", "petiole_width_mm",
    "blade_width_mm", "solidity", "dissection_index", "q25_solidity", "dissected_leaf_fraction",
    "cauline_solidity_med", "cauline_dissection_index_med", "cauline_aspect_ratio_med",
    "cauline_blade_length_mm", "basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta",
    "cauline_cie_l_q85_sheet_delta", "basal_to_cauline_wool_delta",
    "heteroblasty_efa_distance", "Cauline_PC1", "Cauline_PC2", "Cauline_PC3", "ECT_PC1", "ECT_PC2",
    "ECT_PC3", "Cauline_ECT_PC1", "Cauline_ECT_PC2", "ect_lobe_complexity",
    "cauline_ect_lobe_complexity", "heteroblasty_ect_distance", "aspect_ratio", "area_px",
    "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "centroid_size",
    "log_centsize", "AR_inv", "AR_inv_std", *SYM_HARMONIC_COLS, "sym_var", "asym_var",
    "pct_asymmetry", "harmonic_roughness", "leaf_quality_weight", "foliar_variance",
    "within_rosette_efa_dispersion", *PC_COLS, "blade_cie_l_norm", "blade_cie_l_sd",
    "blade_cie_l_sheet_delta", "blade_cie_a_delta", "blade_cie_b_delta", "base_angle_deg",
    "apex_angle_deg", "lamina_curvature", "lamina_trace_length_mm", "basal_cordate_index",
    "margin_serration_index", "margin_tooth_acuteness", "terminal_lobe_area_frac",
    "blade_cie_l_q85_sheet_delta", "basal_zone_cie_l_delta", "abaxial_leaf_fraction",
    "abaxial_cie_l_sheet_delta", "log1p_capitulum_count", "log1p_capitulum_count_std",
    "gmm_uncertainty", "cda_posterior_prob", "cda_oof_posterior_prob", "cl_target_threshold",
    "self_morph_mahal_d2", "cda_second_posterior_prob", "cda_posterior_margin", "loocv_accuracy",
    "can1", "can2", "can3", "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2",
    "type_morph_mahal_d2", "nearest_type_basionym_d2", "type_morph_concordance",
    "min_congener_mahal_d2", "open_set_novelty_p",
)

def persist_otu_morphometrics(
    conn: sqlite3.Connection,
    run_id: str,
    otu_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    preserve_existing_on_null: bool = True,
) -> int:
    """
    Persist Step 03/04 specimen-level OTU medians (`A1..D24`, `PC1..PC15`,
    adaptive EFD/feature metadata (`efd_harmonic_opt`,
    `efd_k95_pcs`, `selected_morph_features`), GMM clusters, and CDA scores
    into `otu_morphometrics`.
    """
    if demo or not otu_rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "otu_morphometrics")
    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for r in otu_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        d: Dict[str, Any] = {
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
        }
        for col in _OTU_MORPH_STR_COLS:
            if col in r:
                d[col] = _to_str(r.get(col))
        if "selected_morph_features" in r:
            d["selected_morph_features"] = _normalize_semicolon_list_str(
                r.get("selected_morph_features")
            )
        for col in _OTU_MORPH_BOOL_COLS:
            if col in r:
                d[col] = _to_bool_int(r.get(col))
        for col in _OTU_MORPH_INT_COLS:
            if col in r:
                d[col] = _to_int(r.get(col))
        for col in _OTU_MORPH_FLOAT_COLS:
            if col in r:
                d[col] = _to_float(r.get(col))
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="otu_morphometrics",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null=preserve_existing_on_null,
    )

_SRF_BOOL_COLS: Tuple[str, ...] = (
    "substrate_specialist_flag", "multimodal_conflict_flag", "matches_type_lithology",
    "env_mcd_outlier_flag",
)
_SRF_INT_COLS: Tuple[str, ...] = (
    "n_observed_modalities",
)
_SRF_STR_COLS: Tuple[str, ...] = (
    "ecotype_regime", "rf_predicted_taxon", "rf_prediction_source", "edaphic_best_fit_taxon",
    "nearest_type_niche_taxon", "nearest_type_niche_basionym", "spatial_block_id",
    "pheno_timing_category", "anthesis_filter_status", "multimodal_concordance", "triage_category",
    "conflict_discordance_notes",
)
_SRF_LIST_STR_COLS: Tuple[str, ...] = (
    "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors",
)
_SRF_FLOAT_COLS: Tuple[str, ...] = (
    "lithic_outcrop_index", "summer_moisture_deficit_idx", "maritime_continentality_idx",
    "rf_predicted_prob", "rf_given_taxon_prob", "spatial_block_cv_acc", "spatial_autocorr_range_km",
    "doy", "expected_doy", "hopkins_expected", "delta_doy", "hopkins_adjusted_doy",
    "hopkins_adjusted_doy_resid", "anthesis_dominance_weight", "allochronic_bimodality_delta_bic",
    "edaphic_z_score", "climatic_z_score", "type_edaphic_z_score",
    "nearest_type_basionym_edaphic_z", "env_mahal_d2", *WC_BIO_COLUMNS,
)

def persist_spatial_rf_predictions(
    conn: sqlite3.Connection,
    run_id: str,
    srf_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    also_upsert_env_cache: bool = False,
) -> int:
    """
    Persist Step 05 Spatial Random Forest predictions, Hopkins phenology residuals,
    all 19 WorldClim bioclimatic variables,
    `climatic_z_score`, and `primary_*_predictors` into `spatial_rf_predictions`.
    """
    if demo or not srf_rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "spatial_rf_predictions")
    if also_upsert_env_cache:
        upsert_environmental_cache_records(conn, srf_rows, demo=demo)

    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for r in srf_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        lat = _to_float(r.get("latitude") if "latitude" in r else r.get("decimalLatitude"))
        lon = _to_float(r.get("longitude") if "longitude" in r else r.get("decimalLongitude"))
        ll_key = _to_str(r.get("lat_lon_key")) or format_lat_lon_key(lat, lon)
        d: Dict[str, Any] = {
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
        }
        if ll_key is not None or "lat_lon_key" in r:
            d["lat_lon_key"] = ll_key
        for col in _SRF_STR_COLS:
            if col in r:
                d[col] = _to_str(r.get(col))
        for col in _SRF_LIST_STR_COLS:
            if col in r:
                d[col] = _normalize_semicolon_list_str(r.get(col))
        for col in _SRF_BOOL_COLS:
            if col in r:
                d[col] = _to_bool_int(r.get(col))
        for col in _SRF_INT_COLS:
            if col in r:
                d[col] = _to_int(r.get(col))
        for col in _SRF_FLOAT_COLS:
            if col in r:
                d[col] = _to_float(r.get(col))
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="spatial_rf_predictions",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null=False,
    )

__all__ = (
    "_LEAF_EFA_BOOL_COLS", "_LEAF_EFA_INT_COLS", "_LEAF_EFA_STR_COLS", "_LEAF_EFA_FLOAT_COLS",
    "_OTU_MORPH_BOOL_COLS", "_OTU_MORPH_INT_COLS", "_OTU_MORPH_STR_COLS", "_OTU_MORPH_FLOAT_COLS",
    "_SRF_BOOL_COLS", "_SRF_INT_COLS", "_SRF_STR_COLS", "_SRF_LIST_STR_COLS", "_SRF_FLOAT_COLS",
    "persist_leaf_efa_harmonics", "persist_otu_morphometrics", "persist_spatial_rf_predictions",
)
