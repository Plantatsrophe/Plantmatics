"""
Static hygiene specifications and per-module strict line-count ceilings
consumed by `scripts/tests/test_code_hygiene.py`.
"""

from __future__ import annotations

from typing import Tuple

ARCHIVE_MARKER = "_" + "archive"

EXCLUDED_DIRS = {
    ARCHIVE_MARKER,
    ".git",
    ".venv",
    ".venv_LM3",
    "LeafMachine3",
    "__pycache__",
}

FORBIDDEN_LEGACY_LM2_TOKENS: Tuple[str, ...] = (
    "LeafMachine2",
    "leafmachine2",
    ".venv_LM2",
    "setup_leafmachine2",
    "configure_leafmachine2",
    "lm2_data_loader",
    "lm2_inference_engines",
    "lm2_geometry_utils",
    "finetune_lm2_pcd",
    "PointRendInferenceEngine",
    "LM2Stage1Detector",
    "detectron2",
)

FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS: Tuple[str, ...] = (
    # SAM-2 & tuning tokens
    "sam2", "SAM2", "annotate_with_sam2", "sam2_annotator_utils", "SAM2BoundaryRefiner",
    "SAM2PromptBundle", "refine_mask_with_sam2", "tune_geometry_parameters", "run_calibrate_geometry",
    # Pruned helper tokens
    "detect_ruler_scale_hough", "trim_petiole_tail", "_sample_cross_sectional_widths",
    "_disambiguate_base_apex_endpoints", "homologize_contour_starting_point",
    "_align_and_normalize_contour", "align_contour_basal_apical", "is_non_lamina_petiole_ribbon",
    "extract_protologue_authors",
    # Pre-LeafMachine3 legacy Python QC, mask-bridging & Tier-2 reflection tokens
    "leaf_gatekeeper", "LeafRoutingResult", "TIER_1_DIRECT", "TIER_2_REFLECTED",
    "Tier_1_Direct", "Tier_2_Reflected", "Tier_1_Pure", "Tier_2_Fallback",
    "tier1_intact", "tier2_reflected", "assigned_tier", "tier_assignment",
    "ucs_score", "min_ucs", "min_solidity", "enable_geometric_qc",
    "is_botanical_dissection", "compute_geometric_metrics", "compute_classical_shape_descriptors",
    "_prepare_lamina_mask", "extract_tier1_pristine", "compute_sheet_adaptive_min_blade_area_cm2",
    "FoldDetectionConfig", "DissectionConfig", "n_tier1_leaves", "n_tier2_leaves",
    "is_folded", "reflection_applied", "render_tier1_vs_tier2_2d_pdf",
    "tier_symmetry_validation", "tier1_vs_tier2_density_overlay", "find_lamina_midrib_endpoints",
    # Purged Steps 1-4 twin columns & circular helpers (Zero-Trust Variable Overhaul)
    "determiner_raw", "pheno_sin", "pheno_cos", "calculate_circular_phenology",
    "capitulum_aspect_ratio", "blade_length_to_width_ratio", "blade_area_cm2",
    "x_norm", "y_norm", "Sym_PC", "is_promoted_anchor", "heteroblasty_index",
    "whiten_intraspecific_deme_plasticity", "compute_cluster_spatial_segregation_p",
    "residualize_within_deme_phenophase_plasticity", ".augment_foliar_architecture_features",
    ".compute_soft_geo_range_log_prior",
    # Purged Step 6 twin columns, obsolete vision hooks, and misnumbered step07/08 helpers
    "w_vision_kappa", "dinov2_predicted_taxon", "vision_conflict_flag", "final_curation_triage",
    "morphometrics_misidentification_flags", "probabilistic_anomaly_score", "composite_conflict_score",
    "integrative_consensus_prob", "n_active_modalities", "render_step07_synthesis_plate_pdf",
    "render_step08_complex_discovery_pdf", "load_step08_merged_inputs",
    "fig07_multimodal_discordance_summary", "fig07b_taxonomic_treatment_summary",
    # Retired Step 06 synthesis / Step 06b / Step 08 helper names and complex-table twin columns
    "06_triage_dashboard_synthesis", "step06b", "Step 06b", "fig06a_integrative_taxonomic_synthesis",
    "fig06b_taxon_complex_delimitation_network", "fig06c_", "ingest_genetic_data",
    "07_discover_taxon_complexes", "render_step07_complex_discovery_pdf",
    ".discover_step08_eco_cols", ".discover_step08_dynamic_morph_cols",
    "resolve_step08_morph_features", "resolve_step08_eco_predictors",
    ".parse_step08_semicolon", ".sort_step08_pcs", "n_Pure_Gold_Silver",
    "leiden_community", "ward_silhouette_cluster", "top_neighbor_1", "top_neighbor_2",
    "kappa_morph_weight", "kappa_eco_weight", "fig_empirical_taxon_complex_network",
    "vision_xai_audit", "dinov2_feature_cache", "cda_biplot_pdf",
    "generate_cda_biplot_pdf", "render_allometry_diagnostic_pdf", "render_cda_biplot_4panel_pdf",
    # Retired Step 02 mask-mutating & duplicated helper tokens
    "_subtract_interior_mounting_paper_pockets", "_build_native_lm3_crop_masks",
    "harmonize_upright_mask_sinus_bridges", "open_hairline_sinus_bridges_and_measure_voids",
    "bridge_already_opened", "compute_eroded_blade_cie_l_metrics",
    # Retired Step 03 dead & duplicated helper tokens
    "summarize_otu_dissection_metrics", "strip_deprecated_step03_cli_args",
    "compute_polygon_solidity",
    # Retired Step 04 dead & duplicated helper tokens
    "generate_allometry_diagnostic_plot", "run_gmm_k1_vs_k2_bootstrap_lrt",
    "generate_cda_biplot_png", ".project_shared_topology_rows_step04",
    # Retired Step 05 pedology, hardpan, & discontinuous composite tokens
    "resdept_all_cm_p", "dom_comp <- tab_res", "100.0 * pmax(0.0, 1.0 - bd_v",
    "0.65 * az_flag", "0.50 * wet_u",
    # LeafMachine2 tokens
    *FORBIDDEN_LEGACY_LM2_TOKENS,
)

RAW_STRICT_LINE_CEILINGS: Tuple[Tuple[str, int], ...] = (
    # Hygiene infrastructure
    ("scripts/tests/_hygiene_specs.py", 320), ("scripts/tests/test_code_hygiene.py", 280),
    # Step 01 core facades, modular helpers, & split test suites
    ("scripts/core/harvester.py", 450), ("scripts/core/taxonomy_authority.py", 260),
    ("scripts/core/usda_plants_ranges.py", 360), ("scripts/core/epa_ecoregions.py", 390),
    ("scripts/core/coord_hygiene_utils.py", 260), ("scripts/core/step01_biogeography_utils.py", 300),
    ("scripts/core/step01_concept_authority_utils.py", 380), ("scripts/data_prep/prune_eastern_cohort.py", 500),
    ("scripts/core/step01_harvester_filters.py", 260), ("scripts/core/step01_harvester_synonyms.py", 400),
    ("scripts/core/step01_harvester_taxon_std.py", 380), ("scripts/core/step01_harvester_prescreen.py", 220),
    ("scripts/core/step01_harvester_media.py", 440), ("scripts/core/step01_harvester_webp_rescue.py", 380),
    ("scripts/core/step01_harvester_export.py", 460), ("scripts/core/step01_harvester_record_ingest.py", 460),
    ("scripts/core/step01_harvester_pagination.py", 460), ("scripts/core/step01_exsiccatae_fingerprint.py", 290),
    ("scripts/core/step01_duplicate_stratifier.py", 460), ("scripts/core/step01_exsiccatae_reconcile.py", 250),
    ("scripts/core/step01_spatial_thinning.py", 460), ("scripts/core/step01_determiner_surnames.py", 410),
    ("scripts/core/step01_corpus_authority_discovery.py", 300), ("scripts/core/step01_corpus_authority_tiers.py", 440),
    ("scripts/core/step01_eastern_cohort_scope.py", 430), ("scripts/core/step01_outgroup_anchors.py", 260),
    ("scripts/core/step01_cec_ecoregion_utils.py", 280), ("scripts/core/step01_usda_symbol_resolver.py", 420),
    ("scripts/core/step01_usda_csv_parser.py", 300), ("scripts/core/step01_usda_distribution_io.py", 440),
    ("scripts/core/step01_county_spatial_graph.py", 320), ("scripts/core/step01_usda_natureserve_utils.py", 220),
    ("scripts/core/step01_usda_voucher_filter.py", 430), ("scripts/tests/test_step01_robustness_audit.py", 380),
    ("scripts/tests/test_step01_taxonomic_clarity.py", 400), ("scripts/tests/test_step01_harvester_geo_cli.py", 460),
    ("scripts/tests/test_step01_resumption_exsiccatae.py", 460), ("scripts/tests/test_step01_thinning_ecoregions.py", 300),
    ("scripts/tests/test_step01_dynamic_authority_synonyms.py", 460), ("scripts/tests/test_step01_genus_wide_harvest.py", 440),
    ("scripts/tests/test_step01_authority_conflict_types.py", 440), ("scripts/tests/test_step01_residency_specialists.py", 320),
    ("scripts/tests/test_step01_monographer_loan_hygiene.py", 260), ("scripts/tests/test_step01_usda_county_ranges.py", 400),
    ("scripts/tests/test_step01_systematics_backbone_hygiene.py", 460), ("scripts/tests/test_step01_synonym_unification_lag.py", 400),
    ("scripts/tests/test_step01_portal_metadata_lag.py", 220),
    # Step 02 vision facades, modular helpers, & tests
    ("scripts/vision/lm3_inference_engines.py", 260), ("scripts/vision/organ_crown_metrics.py", 180),
    ("scripts/vision/capitulum_phenology_classifier.py", 360), ("scripts/vision/lm3_compat.py", 280),
    ("scripts/vision/lm3_data_loader.py", 260), ("scripts/vision/step02_lm3_bbox_capitula.py", 410),
    ("scripts/vision/step02_lm3_onnx_runtime.py", 410), ("scripts/vision/step02_lm3_plant_detector.py", 470),
    ("scripts/vision/step02_lm3_crop_reducer.py", 230), ("scripts/vision/step02_lm3_leaf_segmenter.py", 470),
    ("scripts/vision/step02_lm3_pose_fallback.py", 360), ("scripts/vision/step02_lm3_pose_engine.py", 310),
    ("scripts/vision/step02_lm3_petiole_kpts.py", 450), ("scripts/vision/step02_lm3_damage_qc.py", 360),
    ("scripts/vision/step02_lm3_leaf_orienter.py", 460), ("scripts/vision/step02_lm3_lattice_ect.py", 390),
    ("scripts/vision/step02_lm3_scale_calibrator.py", 395), ("scripts/vision/step02_plant_cluster_attach.py", 320),
    ("scripts/vision/step02_plant_cluster_merge.py", 270), ("scripts/vision/step02_rosette_crown_pose.py", 410),
    ("scripts/vision/step02_rosette_crown_solver.py", 330), ("scripts/vision/step02_petiole_scape_metrics.py", 390),
    ("scripts/vision/step02_leaf_zone_geometry.py", 460), ("scripts/vision/step02_leaf_zone_classifier.py", 470),
    ("scripts/vision/step02_capitulum_biomarkers_early.py", 410), ("scripts/vision/step02_capitulum_biomarkers_late.py", 270),
    ("scripts/vision/step02_lm3_compat_geometry.py", 400), ("scripts/vision/step02_lm3_compat_ect_metrics.py", 350),
    ("scripts/vision/step02_ruler_calibration_io.py", 400),
    ("scripts/vision/step02_leaf_traits_utils.py", 360), ("scripts/vision/step02_photometry_corymb_utils.py", 360),
    ("scripts/vision/leaf_geometry_utils.py", 480), ("scripts/tests/test_step02_vision_hygiene.py", 380),
    ("scripts/tests/test_step02_taxonomic_clarity.py", 400),
    # Step 02B pipeline, rescorer, figure facades, submodules, & split test suites
    ("scripts/pipeline/02_segment_and_extract.py", 470), ("scripts/vision/rescore_step02_qc_and_zones.py", 260),
    ("scripts/vision/step02_figures.py", 270), ("scripts/vision/step02_pipe_init_checkpoint.py", 430),
    ("scripts/vision/step02_pipe_stage1_candidates.py", 380), ("scripts/vision/step02_pipe_stage2_inference.py", 430),
    ("scripts/vision/step02_pipe_repro_clustering.py", 350), ("scripts/vision/step02_pipe_leaf_readout_a.py", 455),
    ("scripts/vision/step02_pipe_leaf_readout_b.py", 475), ("scripts/vision/step02_pipe_manifest_flush.py", 380),
    ("scripts/vision/step02_pipe_wave_runner.py", 455), ("scripts/vision/step02_pipe_cli_args.py", 360),
    ("scripts/vision/step02_rescore_workers.py", 380), ("scripts/vision/step02_rescore_photometry.py", 380),
    ("scripts/vision/step02_rescore_sqlite_sync.py", 430), ("scripts/vision/step02_fig_style_assets.py", 460),
    ("scripts/vision/step02_fig_ect_exemplars.py", 390), ("scripts/vision/step02_fig_sheet_overlays.py", 380),
    ("scripts/vision/step02_fig_plate02a.py", 445), ("scripts/vision/step02_fig_plate02b.py", 360),
    ("scripts/tests/test_step02_core_metrics_resumption.py", 360), ("scripts/tests/test_step02_routing_repro_clustering.py", 360),
    ("scripts/tests/test_step02_wave_engine_sqlite.py", 380), ("scripts/tests/test_step02_letterbox_unet_cli.py", 460),
    ("scripts/tests/test_step02_lm3_delegation_disentangle.py", 440), ("scripts/tests/test_step02_poles_upright_masks.py", 460),
    ("scripts/tests/test_step02_audit_figures_12plate.py", 440), ("scripts/tests/test_step02_photometry_void_qc.py", 360),
    ("scripts/tests/test_step02_shoot_zonation_rescorer.py", 320),
    # Step 03 morphometrics & tests
    ("scripts/morphometrics/step03_config_scale_utils.R", 240), ("scripts/morphometrics/step03_ect_traits_utils.R", 360),
    ("scripts/morphometrics/03_fourier_extractor.R", 500), ("scripts/morphometrics/fourier_parallel_utils.R", 500),
    ("scripts/morphometrics/plasticity_mitigation_utils.R", 470), ("scripts/figures/step03_figures.R", 385),
    ("scripts/tests/test_step03_audit_contracts.py", 360), ("scripts/tests/test_step03_taxonomic_clarity.py", 380),
    # Step 04 morphometrics & tests
    ("scripts/morphometrics/04_gmm_morphotools.R", 500), ("scripts/morphometrics/step04_allometry_merger_utils.R", 380),
    ("scripts/morphometrics/step04_trait_projection_utils.R", 360), ("scripts/morphometrics/gmm_icl_bootstrap_utils.R", 280),
    ("scripts/morphometrics/pub_synthesis_plates.R", 120), ("scripts/figures/step04_figures.R", 500),
    ("scripts/tests/test_step04_audit_hygiene.py", 380), ("scripts/tests/test_step04_taxonomic_clarity.py", 380),
    ("scripts/tests/test_step04_05_model_selection.py", 400),
    # Step 05 macroecology, data prep, & tests
    ("scripts/analysis/05_multimodal_spatial_rf.R", 500), ("scripts/analysis/step05_env_raster_extract_utils.R", 420),
    ("scripts/analysis/step05_spatial_rf_cv_utils.R", 500), ("scripts/analysis/step05_consensus_niche_utils.R", 500),
    ("scripts/analysis/edaphic_fine_grain_utils.R", 500), ("scripts/analysis/edaphic_ssurgo_sda_utils.R", 480),
    ("scripts/analysis/edaphic_solus_3dep_utils.R", 480), ("scripts/analysis/edaphic_collinearity_contingency_utils.R", 480),
    ("scripts/analysis/edaphic_geochem_wetland_utils.R", 340), ("scripts/analysis/step05_circular_hopkins_models.R", 340),
    ("scripts/analysis/step05_phenology_ecotypes.R", 500), ("scripts/analysis/step05_allochronic_pheno_utils.R", 480),
    ("scripts/analysis/step05_type_niche_utils.R", 480), ("scripts/analysis/step05_config_geo_utils.R", 320),
    ("scripts/analysis/step05_layer_hygiene_utils.R", 280), ("scripts/analysis/step05_env_outlier_utils.R", 280),
    ("scripts/analysis/step05_sdm_maxent.R", 500), ("scripts/analysis/step05_sdm_engine_utils.R", 500),
    ("scripts/data_prep/download_environmental_layers.py", 480), ("scripts/data_prep/env_sgmc_3dep_utils.py", 480),
    ("scripts/data_prep/download_geochem_wetland_layers.py", 450), ("scripts/data_prep/download_taxonomic_backbone.py", 450),
    ("scripts/figures/step05_figures.R", 500), ("scripts/figures/step05_sdm_figures.R", 320),
    ("scripts/tests/test_step05_audit_hygiene.py", 380), ("scripts/tests/test_step05_taxonomic_clarity.py", 380),
    ("scripts/tests/test_step05_layer_hygiene.py", 380), ("scripts/tests/test_multimodal_spatial_rf.py", 500),
    ("scripts/tests/test_step05_edaphic_niche.py", 500), ("scripts/tests/test_step05_zero_trust_enm.py", 500),
    ("scripts/tests/test_external_data_sources.py", 440), ("scripts/tests/test_outlier_quarantine_gates.py", 430),
    # Step 06 genomics, cytology, & tests
    ("scripts/data_prep/download_public_packera_sequences.py", 500), ("scripts/data_prep/step06_dryad_backbone_utils.py", 240),
    ("scripts/data_prep/ena_paraloss_discovery.py", 340), ("scripts/data_prep/download_public_cytology.py", 450),
    ("scripts/data_prep/step06_cytology_parse_utils.py", 220), ("scripts/data_prep/organellar_cytonuclear_utils.py", 300),
    ("scripts/data_prep/skim_offtarget_organellar.py", 460), ("scripts/data_prep/06_ingest_genetic_data.py", 500),
    ("scripts/data_prep/step06_synonym_utils.py", 200), ("scripts/data_prep/step06_sample_metadata_utils.py", 460),
    ("scripts/data_prep/step06_voucher_linker_utils.py", 460), ("scripts/data_prep/step06_vcf_dosage_utils.py", 460),
    ("scripts/data_prep/step06_ploidy_hdplot_utils.py", 460), ("scripts/data_prep/step06_public_cohort_ingestor.py", 450),
    ("scripts/data_prep/step06_public_tree_qc_utils.py", 200), ("scripts/analysis/06_genomic_phylogenomics.R", 500),
    ("scripts/analysis/step06_dynamic_tier_engine.R", 470), ("scripts/analysis/step06_phylogenomics_utils.R", 480),
    ("scripts/analysis/step06_model_selection.R", 380), ("scripts/analysis/step06_public_backbone.R", 475),
    ("scripts/analysis/step06_organellar_containment.R", 200), ("scripts/analysis/step06_rogue_sequence_qc.R", 320),
    ("scripts/analysis/step06_reticulation_coancestry.R", 500), ("scripts/analysis/step06_fbranch_quartets.R", 380),
    ("scripts/analysis/step06_coalescent_models.py", 420), ("scripts/analysis/step06_coalescent_alignment_utils.py", 480),
    ("scripts/analysis/step06_coalescent_delimitation.py", 500), ("scripts/analysis/step06_spatial_landscape_genomics.R", 380),
    ("scripts/figures/step06_figures.R", 480), ("scripts/figures/step06_tier_figures.R", 480),
    ("scripts/figures/step06c_figures.R", 350), ("scripts/tests/test_public_cytology_integration.py", 420),
    ("scripts/tests/test_step06_model_selection.py", 420), ("scripts/tests/test_step06_audit_hygiene.py", 400),
    ("scripts/tests/test_step06_dynamic_tiers.py", 420), ("scripts/tests/test_step06_ingestion_hdplot_ploidy.py", 460),
    ("scripts/tests/test_step06_organellar_coalescent.py", 460), ("scripts/tests/test_step06_phylogenomics_r_pipeline.py", 460),
    ("scripts/tests/test_step06_public_backbone_dryad.py", 460), ("scripts/tests/test_step06_voucher_exsiccatae_linking.py", 460),
    # Step 07 synthesis, figures, & tests
    ("scripts/analysis/07_triage_dashboard_synthesis.R", 500), ("scripts/analysis/step07_hopkins_consensus_utils.R", 480),
    ("scripts/analysis/step07_ancova_merge_utils.R", 490), ("scripts/analysis/step07_reports_html_utils.R", 490),
    ("scripts/analysis/step07_biogeo_hybrid_utils.R", 490), ("scripts/analysis/step07_cluster_cohesion_utils.R", 460),
    ("scripts/analysis/step07_decision_stages_utils.R", 490), ("scripts/analysis/step07_postprocess_guards_utils.R", 460),
    ("scripts/analysis/step07_cytology_integration.R", 380), ("scripts/analysis/step07_deme_propagation_utils.R", 260),
    ("scripts/analysis/step07_genomic_propagation.R", 400), ("scripts/analysis/step07_organellar_triage_utils.R", 260),
    ("scripts/analysis/taxon_complex_discovery_utils.R", 480), ("scripts/analysis/step07_complex_affinity_utils.R", 450),
    ("scripts/figures/step07_figures.R", 420), ("scripts/figures/step07b_network_figures.R", 320),
    ("scripts/tests/test_step07_deme_propagation.py", 300), ("scripts/tests/test_step07_organellar_triage.py", 280),
    ("scripts/tests/test_step07_genomic_propagation_contract.py", 460), ("scripts/tests/test_step07_audit_hygiene.py", 400),
    ("scripts/tests/test_parallel_thread_scaling.py", 300),
)

EXTRA_500_LINE_BOUNDED_MODULES: Tuple[str, ...] = (
    "scripts/analysis/genomic_triage_utils.R",
    "scripts/pipeline/resolve_eastern_na_packera.py",
    "scripts/tests/test_eastern_na_monograph_cli.py",
)

