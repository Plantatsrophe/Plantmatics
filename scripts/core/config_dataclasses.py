"""Strongly-typed configuration dataclasses for Plantmatics pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

from scripts.core.config_run_context import PROJECT_ROOT

@dataclass(frozen=True)
class PathsConfig:
    """Strongly-typed, project-root anchored filesystem paths."""
    workspace_root: Path = PROJECT_ROOT
    db_path: Path = PROJECT_ROOT / "data" / "plantmatics.sqlite"
    raw_vouchers_dir: Path = PROJECT_ROOT / "data" / "raw_vouchers"
    curated_vouchers_csv: Path = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
    failed_qc_vouchers_csv: Path = PROJECT_ROOT / "data" / "tables" / "failed_qc_vouchers.csv"
    extraction_manifest_csv: Path = PROJECT_ROOT / "data" / "tables" / "dataset_manifest.csv"
    contours_dir: Path = PROJECT_ROOT / "data" / "contours"
    masks_dir: Path = PROJECT_ROOT / "data" / "masks"
    leaf_efa_harmonics_csv: Path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics.csv"
    morphometric_flags_csv: Path = PROJECT_ROOT / "data" / "tables" / "morphometric_flags.csv"
    cda_biplot_png: Path = PROJECT_ROOT / "outputs" / "figures" / "fig07a_integrative_taxonomic_synthesis.png"
    revision_report_md: Path = PROJECT_ROOT / "outputs" / "reports" / "Packera_Taxonomic_Revision_Summary.md"
    gmm_report_csv: Path = PROJECT_ROOT / "outputs" / "reports" / "gmm_bayes_factors_summary.csv"
    summary_log: Path = PROJECT_ROOT / "outputs" / "reports" / "voucher_ingestion_summary.log"
    lm3_root: Path = PROJECT_ROOT / "LeafMachine3"
    models_dir: Path = PROJECT_ROOT / "models" / "lm3"
    model_weights: Path = (
        PROJECT_ROOT / "models" / "lm3" / "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
    )
    quarantine_dir: Path = PROJECT_ROOT / "data" / "raw_vouchers_quarantine"
    tables_dir: Path = PROJECT_ROOT / "data" / "tables"
    figures_dir: Path = PROJECT_ROOT / "outputs" / "figures"


@dataclass(frozen=True)
class TaxaConfig:
    """Target focal taxa, synonymy mappings, and outgroup query sets."""
    target_species: List[str] = field(default_factory=list)
    synonyms: Dict[str, List[str]] = field(default_factory=dict)
    outgroups: List[str] = field(default_factory=list)
    target_genus: str = "Packera"
    target_family: str = "Asteraceae"
    synonym_genera: List[str] = field(default_factory=lambda: ["Packera", "Senecio", "Cineraria"])
    eastern_resident_whitelist: List[str] = field(default_factory=list)
    transcontinental_boreal_taxa: List[str] = field(
        default_factory=lambda: ["Packera paupercula", "Packera indecora", "Packera pauciflora"]
    )

    @property
    def target_taxa(self) -> List[str]:
        return self.target_species


@dataclass(frozen=True)
class CapitulumPhenologyConfig:
    """Thresholds for capitulum phenology biomarker classification."""
    min_pappus_luminance: float = 75.0
    max_pappus_saturation: float = 0.18
    min_pappus_area_ratio: float = 0.18
    min_floret_yellow_saturation: float = 0.32
    min_floret_yellow_ratio: float = 0.12
    min_anthesis_heads_for_flowering: int = 1


VALID_GEOGRAPHIC_SCOPES: Set[str] = {"eastern_na", "north_america", "global"}


def normalize_countries_list(raw_countries: Any, default: Optional[List[str]] = None) -> List[str]:
    """Normalize country configuration values (list or comma/whitespace-separated string) to uppercase ISO-2 codes."""
    if raw_countries is None:
        return list(default) if default is not None else ["US", "CA"]
    items: List[str] = []
    if isinstance(raw_countries, str):
        raw_items = [tok for part in raw_countries.split(",") for tok in part.strip().split()]
    elif isinstance(raw_countries, (list, tuple, set)):
        raw_items = [
            tok
            for elem in raw_countries
            for part in str(elem).split(",")
            for tok in part.strip().split()
        ]
    else:
        raw_items = [str(raw_countries).strip()]
    for tok in raw_items:
        code = tok.strip().upper()
        if code and code not in items:
            items.append(code)
    return items


@dataclass(frozen=True)
class ThresholdsConfig:
    """Botanical optical quality and harvesting thresholds."""
    pcd_conf: float = 0.35
    minimum_confidence_threshold: float = 0.65
    min_megapixels: float = 8.0
    min_file_size_kb: float = 500.0
    min_sharpness_laplacian: float = 80.0
    max_uncertainty_meters: float = 5000.0
    western_longitude_threshold: float = -100.0
    geographic_scope: str = "eastern_na"
    countries: List[str] = field(default_factory=lambda: ["US", "CA"])
    capitulum_phenology: CapitulumPhenologyConfig = field(default_factory=CapitulumPhenologyConfig)


@dataclass(frozen=True)
class NormalizationConfig:
    """
    Homologized Elliptic Fourier contour invariants.
    Scale, rotation, and start-point homology are enforced via LeafMachine3's
    `lamina_base -> lamina_tip` anatomical frame (`coo_center() %>% coo_scale()`
    followed by `Momocs::efourier(..., norm = FALSE, start = FALSE)`), never via
    first-ellipse major-axis normalization (`efourier_norm_ellipse = False`).
    """
    align_major_axis: bool = True
    scale_invariant: bool = True
    rotation_invariant: bool = True
    start_point_invariant: bool = True
    homologous_base_to_tip: bool = True
    efourier_norm_ellipse: bool = False
    efourier_start_search: bool = False


@dataclass(frozen=True)
class MorphometricsConfig:
    """Elliptic Fourier Analysis and downstream clustering parameters."""
    cohort_adaptive_mode: bool = True
    min_anchor_n_for_empirical: int = 5
    nb_harmonics: int = 12
    harmonics: int = 12
    min_harmonics: int = 6
    max_harmonics: int = 24
    target_harmonic_power: float = 0.99
    target_pca_variance: float = 0.95
    num_pcs: int = 5
    min_export_pcs: int = 5
    max_export_pcs: int = 15
    max_k: int = 8
    random_seed: int = 42
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    correct_allometry: str = "auto"
    allometry_r2_threshold: float = 0.10
    symmetrize_efa_harmonics: bool = True


@dataclass(frozen=True)
class PhenologyConfig:
    """Phenological modeling, latitudinal spring progression, and anomaly detection."""
    min_flowering_doy: int = 30
    max_flowering_doy: int = 220
    robust_regression: bool = True


@dataclass(frozen=True)
class MacroecologyConfig:
    """Phase 6 Macroecology, spatial random forest, and phenological parameters."""
    phenology: PhenologyConfig = field(default_factory=PhenologyConfig)



@dataclass(frozen=True)
class HarvestingConfig:
    """Phase 1 GBIF query, download, and filtration hyperparameters."""
    max_records_per_taxon: int = 5000
    max_uncertainty_meters: float = 5000.0
    geographic_scope: str = "eastern_na"
    countries: List[str] = field(default_factory=lambda: ["US", "CA"])
    exclude_western: bool = True
    harvest_full_continental_range: bool = False
    min_eastern_residency_records: int = 10
    western_longitude_threshold: float = -100.0
    cis_mississippi_longitude_threshold: float = -90.5
    min_cis_mississippi_fraction: float = 0.05
    transcontinental_ecoregion_gate: bool = True
    transcontinental_longitude_threshold: float = -100.0
    home_herbarium: str = "NCU"
    min_megapixels: float = 8.0
    min_file_size_kb: float = 500.0
    min_sharpness_laplacian: float = 80.0
    check_sharpness: bool = True
    webp_quality: int = 92
    download_concurrency: int = 15
    usda_county_filter: bool = False
    usda_county_buffer_km: float = 50.0
    usda_exempt_gold: bool = True
    include_ungeoreferenced_state_vouchers: bool = True


@dataclass(frozen=True)
class SegmentationConfig:
    """Phase 2 LeafMachine3 ONNX segmentation and morphological extraction parameters."""
    lm3_root: Path = PROJECT_ROOT / "LeafMachine3"
    models_dir: Path = PROJECT_ROOT / "models" / "lm3"
    model_weights: Path = (
        PROJECT_ROOT / "models" / "lm3" / "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
    )
    plant_detector_model: str = "lm3_plant_detector__yolo26x_det_1280"
    archival_detector_model: str = "lm3_archival_detector__yolo26x_det_1280"
    specimen_segmenter_model: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
    leaf_segmenter_model: str = "lm3_leaf_segmenter__yolo26x_seg_1024"
    landmark_detector_model: str = "lm3_landmark_detector__yolo26x_pose_640"
    mp_conversion_model: str = "lm3_mp_conversion_factor__sqrt_fit"
    detector_version: str = "lm3_plant_detector__yolo26x_det_1280"
    max_working_dim: int = 3200
    specimen_max_working_dim: int = 3200
    paperclean: bool = True
    llr_hyst_keep_blobs: bool = True
    paperclean_llr_hyst_keepblobs: bool = True
    subtract_holes: bool = True
    use_cf_predicted_by_mp: bool = False
    ect_num_dirs: int = 360
    minimum_confidence_threshold: float = 0.35
    treat_leaflet_as_leaf_whole: bool = False
    use_efds_for_png_masks: bool = False
    do_apply_landmarks: bool = True
    segment_partial_leaves: bool = False
    rosette_plant_bbox_lower_fraction: float = 0.60
    score_thresh: float = 0.35
    save_debug_crops: bool = False
    device: str = "cuda"


@dataclass(frozen=True)
class ModelsConfig:
    """LeafMachine3 ONNX model identifiers, checkpoint paths, remote URLs, and integrity hashes."""
    lm3_root: Path = PROJECT_ROOT / "LeafMachine3"
    models_dir: Path = PROJECT_ROOT / "models" / "lm3"
    plant_detector_model: str = "lm3_plant_detector__yolo26x_det_1280"
    archival_detector_model: str = "lm3_archival_detector__yolo26x_det_1280"
    specimen_segmenter_model: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
    leaf_segmenter_model: str = "lm3_leaf_segmenter__yolo26x_seg_1024"
    landmark_detector_model: str = "lm3_landmark_detector__yolo26x_pose_640"
    mp_conversion_model: str = "lm3_mp_conversion_factor__sqrt_fit"
    pcd_weights_path: Path = (
        PROJECT_ROOT / "models" / "lm3" / "lm3_plant_detector__yolo26x_det_1280.onnx"
    )
    pcd_weights_url: str = ""
    pcd_weights_sha256: str = ""


@dataclass(frozen=True)
class ComputeConfig:
    """CPU multithreading and GPU inference batch size configuration."""
    threads: Optional[int] = None
    gpu_batch_size: Optional[int] = None
