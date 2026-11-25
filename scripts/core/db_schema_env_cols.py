"""
Column-set constants for harmonic, PCA, WorldClim, and environmental cache SQLite schemas.
"""

from __future__ import annotations

from typing import Set, Tuple

HARMONIC_COLS: Tuple[str, ...] = tuple(
    f"{coef}{idx}" for idx in range(1, 25) for coef in ("A", "B", "C", "D")
)

SYM_HARMONIC_COLS: Tuple[str, ...] = tuple(
    f"{coef}{idx}" for idx in range(1, 25) for coef in ("A", "D")
)

PC_COLS: Tuple[str, ...] = tuple(f"PC{idx}" for idx in range(1, 16))

WC_BIO_COLUMNS: Tuple[str, ...] = (
    "bio1_temp_mean",
    "bio2_mean_diurnal_range",
    "bio3_isothermality",
    "bio4_temp_seasonality",
    "bio5_max_temp_warmest_month",
    "bio6_min_temp_coldest_month",
    "bio7_temp_annual_range",
    "bio8_mean_temp_wettest_qtr",
    "bio9_mean_temp_driest_qtr",
    "bio10_mean_temp_warmest_qtr",
    "bio11_mean_temp_coldest_qtr",
    "bio12_precip_annual",
    "bio13_precip_wettest_month",
    "bio14_precip_driest_month",
    "bio15_precip_seasonality",
    "bio16_precip_wettest_qtr",
    "bio17_precip_driest_qtr",
    "bio18_precip_warmest_qtr",
    "bio19_precip_coldest_qtr",
)

ENVIRONMENTAL_CACHE_COLUMNS: Tuple[str, ...] = (
    "lat_lon_key",
    "catalogNumber",
    "latitude",
    "longitude",
    "musym",
    "muname",
    "taxsuborder",
    "taxgreatgroup",
    "is_rock_outcrop",
    "outcrop_class",
    "ssurgo_rock_outcrop_pct",
    "ssurgo_lithic_series_pct",
    "ssurgo_bedrock_depth_cm",
    "ssurgo_sand_pct",
    "ssurgo_clay_pct",
    "ssurgo_ph",
    "ssurgo_om_pct",
    "ssurgo_cec",
    "ssurgo_aws_0_50cm",
    "polaris_ph_30m",
    "polaris_sand_30m",
    "polaris_clay_30m",
    "polaris_om_30m",
    "polaris_bd_30m",
    "polaris_ksat_30m",
    "solus_resdept_cm",
    "solus_fragvol_pct",
    "solus_sand_pct",
    "solus_clay_pct",
    "solus_silt_pct",
    "solus_cec",
    "topo_tpi_30m",
    "topo_slope_deg",
    "topo_hli_30m",
    "topo_twi_30m",
    "elevation_m",
    "fine_soil_ph",
    "fine_soil_sand",
    "fine_soil_clay",
    "soil_silt_pct",
    "ilr_texture_1",
    "ilr_texture_2",
    "fine_soil_cec",
    "fine_cfvo_pct",
    "fine_bedrock_depth_cm",
    "sgmc_lithology_class",
    "is_azonal_lithology",
    "buffer_radius_m",
    "high_coord_uncertainty_flag",
    "is_training_eligible_coord",
    "is_clim_eligible_coord",
    "is_edaphic_eligible_coord",
    "edaphic_coord_weight",
    "ph_q10_buffer",
    "ph_q90_buffer",
    "min_bedrock_depth_buffer_cm",
    "max_cfvo_buffer_pct",
    "azonal_lithology_presence_buffer",
    *WC_BIO_COLUMNS,
    "soil_ph",
    "soil_cec",
    "soil_sand",
    "soil_clay",
    "soil_silt",
    "soil_cfvo",
    "soil_bulk_density",
    "soil_ph_soilgrids",
    "soil_cec_soilgrids",
    "soil_sand_soilgrids",
    "soil_clay_soilgrids",
    "soil_silt_soilgrids",
    "soil_cfvo_soilgrids",
    "soil_bulk_density_soilgrids",
    "soil_ca_mg_molar_ratio",
    "soil_ultramafic_ni_cr_index",
    "soil_carbonate_c_pct",
    "is_nwi_wetland",
    "nwi_wetland_system",
    "nwi_wetland_dist_m",
)

ENVIRONMENTAL_TEXT_COLS: Set[str] = {
    "lat_lon_key",
    "catalogNumber",
    "musym",
    "muname",
    "taxsuborder",
    "taxgreatgroup",
    "outcrop_class",
    "sgmc_lithology_class",
    "nwi_wetland_system",
}

ENVIRONMENTAL_BOOL_COLS: Set[str] = {
    "is_rock_outcrop",
    "is_azonal_lithology",
    "high_coord_uncertainty_flag",
    "is_training_eligible_coord",
    "is_clim_eligible_coord",
    "is_edaphic_eligible_coord",
    "azonal_lithology_presence_buffer",
    "is_nwi_wetland",
}

__all__ = (
    "HARMONIC_COLS", "SYM_HARMONIC_COLS", "PC_COLS", "WC_BIO_COLUMNS",
    "ENVIRONMENTAL_CACHE_COLUMNS", "ENVIRONMENTAL_TEXT_COLS", "ENVIRONMENTAL_BOOL_COLS",
)
