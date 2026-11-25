"""Centralized `PipelineConfig` container and YAML loader for Plantmatics."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import os
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Union

import yaml

from scripts.core.config_dataclasses import (
    CapitulumPhenologyConfig,
    ComputeConfig,
    HarvestingConfig,
    MacroecologyConfig,
    ModelsConfig,
    MorphometricsConfig,
    NormalizationConfig,
    PathsConfig,
    PhenologyConfig,
    SegmentationConfig,
    TaxaConfig,
    ThresholdsConfig,
    VALID_GEOGRAPHIC_SCOPES,
    normalize_countries_list,
)
from scripts.core.config_run_context import (
    DEFAULT_CONFIG_PATH,
    PROJECT_ROOT,
    RUN_SCOPED_PATH_FIELDS,
)

@dataclass
class PipelineConfig:
    """Centralized, type-hinted configuration for Packera morphometrics pipeline."""
    paths: PathsConfig
    taxa: TaxaConfig
    thresholds: ThresholdsConfig
    morphometrics: MorphometricsConfig
    harvesting: HarvestingConfig
    segmentation: SegmentationConfig
    models: ModelsConfig = field(default_factory=ModelsConfig)
    macroecology: MacroecologyConfig = field(default_factory=MacroecologyConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    _raw_dict: Dict[str, Any] = field(default_factory=dict, repr=False)

    def __getitem__(self, key: str) -> Any:
        return self._raw_dict[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self._raw_dict.get(key, default)

    def to_dict(self) -> Dict[str, Any]:
        return dict(self._raw_dict)

    def with_run_dir(self, run_dir: Union[str, Path]) -> PipelineConfig:
        """
        Return a new PipelineConfig whose run-scoped output/table/figure paths
        point inside `run_dir`, while shared read-only caches (raw_vouchers_dir,
        quarantine_dir, models_dir, model_weights) remain anchored at PROJECT_ROOT.
        """
        rdir = Path(run_dir).resolve()
        new_paths = replace(
            self.paths,
            curated_vouchers_csv=rdir / "data" / "tables" / "curated_vouchers.csv",
            failed_qc_vouchers_csv=rdir / "data" / "tables" / "failed_qc_vouchers.csv",
            extraction_manifest_csv=rdir / "data" / "tables" / "extracted_leaf_manifest.csv",
            contours_dir=rdir / "data" / "contours",
            masks_dir=rdir / "data" / "masks",
            leaf_efa_harmonics_csv=rdir / "data" / "tables" / "leaf_efa_harmonics.csv",
            morphometric_flags_csv=rdir / "data" / "tables" / "morphometric_flags.csv",
            cda_biplot_png=rdir / "outputs" / "figures" / "fig07a_integrative_taxonomic_synthesis.png",
            revision_report_md=rdir / "outputs" / "reports" / "Packera_Taxonomic_Revision_Summary.md",
            gmm_report_csv=rdir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv",
            summary_log=rdir / "logs" / "pipeline.log",
            tables_dir=rdir / "data" / "tables",
            figures_dir=rdir / "outputs" / "figures",
        )
        new_raw = dict(self._raw_dict)
        raw_paths = dict(new_raw.get("paths", {}))
        for k in RUN_SCOPED_PATH_FIELDS:
            raw_paths[k] = str(getattr(new_paths, k))
        new_raw["paths"] = raw_paths
        new_raw["resolved_paths"] = {
            k: getattr(new_paths, k) for k in new_paths.__dataclass_fields__
        }
        return replace(self, paths=new_paths, _raw_dict=new_raw)

    @classmethod
    def from_yaml(
        cls,
        yaml_path: Optional[Union[str, Path]] = None,
        run_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineConfig:
        """Dynamically parse YAML and resolve paths relative to project root (or active run_dir if supplied)."""
        _cfg_mod = sys.modules.get("scripts.core.config")
        _default_cfg_path = getattr(_cfg_mod, "DEFAULT_CONFIG_PATH", DEFAULT_CONFIG_PATH)
        _proj_root = getattr(_cfg_mod, "PROJECT_ROOT", PROJECT_ROOT)
        target = Path(yaml_path) if yaml_path else _default_cfg_path
        target = target if target.is_absolute() else (_proj_root / target).resolve()

        if not target.exists():
            raise FileNotFoundError(f"Configuration file not found at: {target}")

        try:
            with open(target, "r", encoding="utf-8") as f:
                raw_cfg = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError(f"Invalid YAML syntax in configuration file '{target}': {exc}") from exc

        if not isinstance(raw_cfg, dict) or not raw_cfg:
            raise ValueError(f"Configuration file at '{target}' is empty or not a valid dictionary mapping.")

        def _resolve(rel_path: Union[str, Path]) -> Path:
            p = Path(rel_path)
            return p if p.is_absolute() else (_proj_root / p).resolve()

        p_kwargs = {
            k: _resolve(v) for k, v in raw_cfg.get("paths", {}).items()
            if k in PathsConfig.__dataclass_fields__
        }
        paths = PathsConfig(**p_kwargs)

        taxa_dict = raw_cfg.get("taxa", {}) or {}
        legacy_taxa = (raw_cfg.get("taxonomy") or {}).get("target_taxa", [])
        target_genus_val = (
            os.environ.get("PLANTMATICS_TARGET_GENUS", "").strip()
            or str(taxa_dict.get("target_genus", "") or "").strip()
            or "Packera"
        )
        target_family_val = (
            os.environ.get("PLANTMATICS_TARGET_FAMILY", "").strip()
            or str(taxa_dict.get("target_family", "") or "").strip()
            or "Asteraceae"
        )
        raw_syn_gen = [
            str(g).strip()
            for g in (taxa_dict.get("synonym_genera") or ("Packera", "Senecio", "Cineraria"))
            if str(g).strip()
        ]
        if target_genus_val and target_genus_val not in raw_syn_gen:
            raw_syn_gen.insert(0, target_genus_val)
        raw_wl = [str(s).strip() for s in (taxa_dict.get("eastern_resident_whitelist") or []) if str(s).strip()]
        raw_tc = [
            str(s).strip()
            for s in (
                taxa_dict.get("transcontinental_boreal_taxa")
                or ("Packera paupercula", "Packera indecora", "Packera pauciflora")
            )
            if str(s).strip()
        ]
        taxa = TaxaConfig(
            target_species=list(taxa_dict.get("target_species", legacy_taxa) or []),
            synonyms=dict(taxa_dict.get("synonyms", {}) or {}),
            outgroups=list(taxa_dict.get("outgroups", []) or []),
            target_genus=target_genus_val,
            target_family=target_family_val,
            synonym_genera=raw_syn_gen,
            eastern_resident_whitelist=raw_wl,
            transcontinental_boreal_taxa=raw_tc,
        )

        thresh_dict = raw_cfg.get("thresholds", {})
        pheno_dict = thresh_dict.get("capitulum_phenology")
        if isinstance(pheno_dict, dict):
            pheno_cfg = CapitulumPhenologyConfig(
                min_pappus_luminance=float(pheno_dict.get("min_pappus_luminance", 75.0)),
                max_pappus_saturation=float(pheno_dict.get("max_pappus_saturation", 0.18)),
                min_pappus_area_ratio=float(pheno_dict.get("min_pappus_area_ratio", 0.18)),
                min_floret_yellow_saturation=float(pheno_dict.get("min_floret_yellow_saturation", 0.32)),
                min_floret_yellow_ratio=float(pheno_dict.get("min_floret_yellow_ratio", 0.12)),
                min_anthesis_heads_for_flowering=int(pheno_dict.get("min_anthesis_heads_for_flowering", 1)),
            )
        else:
            pheno_cfg = CapitulumPhenologyConfig()

        harv_dict = raw_cfg.get("harvesting", {})
        non_float_thresh_fields = {"capitulum_phenology", "geographic_scope", "countries"}
        thresh_kwargs: Dict[str, Any] = {
            k: float(v)
            for k, v in thresh_dict.items()
            if k in ThresholdsConfig.__dataclass_fields__
            and k not in non_float_thresh_fields
            and not isinstance(v, (dict, list))
        }
        if "pcd_conf" not in thresh_dict:
            pcd_min_conf = (raw_cfg.get("plant_component_detector") or {}).get("minimum_confidence_threshold")
            if pcd_min_conf is not None:
                thresh_kwargs["pcd_conf"] = float(pcd_min_conf)
        raw_t_scope = str(thresh_dict.get("geographic_scope", harv_dict.get("geographic_scope", "eastern_na"))).strip().lower()
        thresh_kwargs["geographic_scope"] = raw_t_scope if raw_t_scope in VALID_GEOGRAPHIC_SCOPES else "eastern_na"
        thresh_kwargs["countries"] = normalize_countries_list(
            thresh_dict.get("countries", harv_dict.get("countries", ["US", "CA"])),
            default=["US", "CA"],
        )
        thresh_kwargs["capitulum_phenology"] = pheno_cfg

        thresholds = ThresholdsConfig(**thresh_kwargs)

        norm_dict = raw_cfg.get("morphometrics", {}).get("normalization", {})
        norm = NormalizationConfig(**{
            k: bool(v) for k, v in norm_dict.items() if k in NormalizationConfig.__dataclass_fields__
        })

        m_dict = raw_cfg.get("morphometrics", {}) or {}
        f_dict = m_dict.get("fourier", {}) if isinstance(m_dict.get("fourier"), dict) else {}
        nb_h = int(m_dict.get("nb_harmonics", m_dict.get("harmonics", 12)))
        morph = MorphometricsConfig(
            cohort_adaptive_mode=bool(f_dict.get("cohort_adaptive_mode", m_dict.get("cohort_adaptive_mode", True))),
            min_anchor_n_for_empirical=int(f_dict.get("min_anchor_n_for_empirical", m_dict.get("min_anchor_n_for_empirical", 5))),
            nb_harmonics=nb_h,
            harmonics=int(m_dict.get("harmonics", nb_h)),
            min_harmonics=int(f_dict.get("min_harmonics", 6)),
            max_harmonics=int(f_dict.get("max_harmonics", 24)),
            target_harmonic_power=float(f_dict.get("target_harmonic_power", 0.99)),
            target_pca_variance=float(f_dict.get("target_pca_variance", 0.95)),
            num_pcs=int(m_dict.get("num_pcs", 5)),
            min_export_pcs=int(f_dict.get("min_export_pcs", 5)),
            max_export_pcs=int(f_dict.get("max_export_pcs", 15)),
            max_k=int(m_dict.get("max_k", 8)),
            random_seed=int(m_dict.get("random_seed", 42)),
            normalization=norm,
            correct_allometry=str(m_dict.get("correct_allometry", "auto")),
            allometry_r2_threshold=float(m_dict.get("allometry_r2_threshold", 0.10)),
            symmetrize_efa_harmonics=bool(m_dict.get("symmetrize_efa_harmonics", True)),
        )

        harv_kwargs: Dict[str, Any] = {
            k: v for k, v in harv_dict.items() if k in HarvestingConfig.__dataclass_fields__
        }
        raw_h_scope = str(harv_kwargs.get("geographic_scope", thresh_kwargs["geographic_scope"])).strip().lower()
        harv_kwargs["geographic_scope"] = raw_h_scope if raw_h_scope in VALID_GEOGRAPHIC_SCOPES else "eastern_na"
        if "countries" in harv_kwargs:
            harv_kwargs["countries"] = normalize_countries_list(harv_kwargs["countries"], default=["US", "CA"])
        elif harv_kwargs["geographic_scope"] == "north_america":
            harv_kwargs["countries"] = ["US", "CA", "MX"]
        elif harv_kwargs["geographic_scope"] == "global":
            harv_kwargs["countries"] = []
        if "home_herbarium" in harv_kwargs:
            harv_kwargs["home_herbarium"] = str(harv_kwargs["home_herbarium"]).strip().upper() or "NCU"
        if "usda_county_filter" in harv_kwargs:
            harv_kwargs["usda_county_filter"] = bool(harv_kwargs["usda_county_filter"])
        if "usda_county_buffer_km" in harv_kwargs and harv_kwargs["usda_county_buffer_km"] is not None:
            harv_kwargs["usda_county_buffer_km"] = float(harv_kwargs["usda_county_buffer_km"])
        if "usda_exempt_gold" in harv_kwargs:
            harv_kwargs["usda_exempt_gold"] = bool(harv_kwargs["usda_exempt_gold"])
        if "transcontinental_ecoregion_gate" in harv_kwargs:
            harv_kwargs["transcontinental_ecoregion_gate"] = bool(harv_kwargs["transcontinental_ecoregion_gate"])
        if "include_ungeoreferenced_state_vouchers" in harv_kwargs:
            harv_kwargs["include_ungeoreferenced_state_vouchers"] = bool(harv_kwargs["include_ungeoreferenced_state_vouchers"])
        harvesting = HarvestingConfig(**harv_kwargs)

        vision_lm3_fallback = (raw_cfg.get("vision") or {}).get("leafmachine3", {}) or {}
        leaf_seg_fallback = raw_cfg.get("leaf_segmentation", {}) or {}
        seg_dict = {**vision_lm3_fallback, **leaf_seg_fallback, **(raw_cfg.get("segmentation", {}) or {})}
        seg_kwargs: Dict[str, Any] = {
            k: v for k, v in seg_dict.items() if k in SegmentationConfig.__dataclass_fields__
        }
        for path_k in ("model_weights", "lm3_root", "models_dir"):
            if path_k in seg_kwargs and seg_kwargs[path_k] is not None:
                seg_kwargs[path_k] = _resolve(seg_kwargs[path_k])
        for bool_k in (
            "paperclean",
            "llr_hyst_keep_blobs",
            "paperclean_llr_hyst_keepblobs",
            "subtract_holes",
            "use_cf_predicted_by_mp",
            "treat_leaflet_as_leaf_whole",
            "use_efds_for_png_masks",
            "do_apply_landmarks",
            "segment_partial_leaves",
            "save_debug_crops",
        ):
            if bool_k in seg_kwargs and seg_kwargs[bool_k] is not None:
                seg_kwargs[bool_k] = bool(seg_kwargs[bool_k])
        for int_k in ("max_working_dim", "specimen_max_working_dim", "ect_num_dirs"):
            if int_k in seg_kwargs and seg_kwargs[int_k] is not None:
                seg_kwargs[int_k] = int(seg_kwargs[int_k])
        for float_k in ("minimum_confidence_threshold", "rosette_plant_bbox_lower_fraction", "score_thresh"):
            if float_k in seg_kwargs and seg_kwargs[float_k] is not None:
                seg_kwargs[float_k] = float(seg_kwargs[float_k])
        segmentation = SegmentationConfig(**seg_kwargs)

        models_dict = raw_cfg.get("models", {})
        models_kwargs: Dict[str, Any] = {
            k: v for k, v in models_dict.items() if k in ModelsConfig.__dataclass_fields__
        }
        for path_k in ("pcd_weights_path", "lm3_root", "models_dir"):
            if path_k in models_kwargs and models_kwargs[path_k] is not None:
                models_kwargs[path_k] = _resolve(models_kwargs[path_k])
        for str_k in (
            "plant_detector_model",
            "archival_detector_model",
            "specimen_segmenter_model",
            "leaf_segmenter_model",
            "landmark_detector_model",
            "mp_conversion_model",
            "pcd_weights_url",
            "pcd_weights_sha256",
        ):
            if str_k in models_kwargs and models_kwargs[str_k] is not None:
                models_kwargs[str_k] = str(models_kwargs[str_k])
        models = ModelsConfig(**models_kwargs)

        macro_dict = raw_cfg.get("macroecology", {})
        pheno_dict = macro_dict.get("phenology", {})
        phenology_cfg = PhenologyConfig(
            min_flowering_doy=int(pheno_dict.get("min_flowering_doy", 30)),
            max_flowering_doy=int(pheno_dict.get("max_flowering_doy", 220)),
            robust_regression=bool(pheno_dict.get("robust_regression", True)),
        )
        macroecology = MacroecologyConfig(phenology=phenology_cfg)

        comp_dict = raw_cfg.get("compute", {}) or {}
        raw_threads = comp_dict.get("threads")
        raw_gpu_bs = comp_dict.get("gpu_batch_size")
        compute_cfg = ComputeConfig(
            threads=int(raw_threads) if raw_threads is not None and int(raw_threads) >= 1 else None,
            gpu_batch_size=int(raw_gpu_bs) if raw_gpu_bs is not None and int(raw_gpu_bs) >= 1 else None,
        )

        resolved_raw = dict(raw_cfg)
        resolved_raw["resolved_paths"] = {
            k: getattr(paths, k) for k in paths.__dataclass_fields__
        }

        instance = cls(
            paths=paths, taxa=taxa, thresholds=thresholds,
            morphometrics=morph, harvesting=harvesting,
            segmentation=segmentation, models=models,
            macroecology=macroecology, compute=compute_cfg,
            _raw_dict=resolved_raw,
        )
        if run_dir is not None:
            return instance.with_run_dir(run_dir)
        return instance
