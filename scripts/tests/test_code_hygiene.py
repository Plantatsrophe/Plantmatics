#!/usr/bin/env python3
"""
===============================================================================
Script: test_code_hygiene.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Static AST and text/grep hygiene test suite verifying that all active .py,
    .R, and .sh files in the repository (excluding archived script folders,
    .git, .venv, .venv_LM3, and LeafMachine3) are 100% isolated from legacy
    archived code/data directories and contain zero legacy LeafMachine2 or
    Detectron2 PointRend remnants.
===============================================================================
"""

from __future__ import annotations

import re
import ast
import sys
import unittest
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

ARCHIVE_MARKER = "_" + "archive"
SELF_PATH = Path(__file__).resolve()

EXCLUDED_DIRS = {
    ARCHIVE_MARKER,
    ".git",
    ".venv",
    ".venv_LM3",
    "LeafMachine3",
    "__pycache__",
}

FORBIDDEN_LEGACY_LM2_TOKENS = (
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

FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS = (
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


def iter_active_source_files(include_shell: bool = False) -> List[Path]:
    """Return all active .py, .R (and optionally .sh) files under scripts/ and root, excluding archive and submodule/venv directories."""
    exts = ("*.py", "*.R", "*.sh") if include_shell else ("*.py", "*.R")
    files: List[Path] = []
    for ext in exts:
        for root_file in sorted(PROJECT_ROOT.glob(ext)):
            resolved = root_file.resolve()
            if resolved.is_file():
                files.append(resolved)

    scripts_dir = PROJECT_ROOT / "scripts"
    for ext in exts:
        for path in sorted(scripts_dir.rglob(ext)):
            resolved = path.resolve()
            if not resolved.is_file():
                continue
            rel_parts = set(resolved.relative_to(PROJECT_ROOT).parts)
            if rel_parts & EXCLUDED_DIRS:
                continue
            files.append(resolved)
    return files


class TestCodeHygiene(unittest.TestCase):
    """Static AST and grep hygiene checks for archive isolation and zero legacy LM2 remnants across active files."""

    def test_no_active_files_reference_legacy_folders_via_grep(self):
        """Scan all active .py, .R, and .sh files and fail if any references the archive folder."""
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            for lineno, line in enumerate(content.splitlines(), start=1):
                if ARCHIVE_MARKER in line:
                    rel = src_path.relative_to(PROJECT_ROOT)
                    violations.append(f"{rel}:{lineno}: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            f"Active source files must not reference '{ARCHIVE_MARKER}':\n" + "\n".join(violations),
        )

    def test_no_files_or_folders_named_lm2_or_leafmachine2(self):
        """Verify that zero files or folders in the workspace use 'lm2' or 'leafmachine2' in their naming."""
        skip_roots = {".git", ".venv", ".venv_LM3", "LeafMachine3", "data", "runs", "outputs"}
        bad_paths: List[str] = []
        for entry in sorted(PROJECT_ROOT.iterdir()):
            name_low = entry.name.lower()
            if "lm2" in name_low or "leafmachine2" in name_low:
                bad_paths.append(str(entry.relative_to(PROJECT_ROOT)))
            if entry.name in skip_roots or not entry.is_dir():
                continue
            for sub in entry.rglob("*"):
                sub_low = sub.name.lower()
                if "lm2" in sub_low or "leafmachine2" in sub_low:
                    bad_paths.append(str(sub.relative_to(PROJECT_ROOT)))

        self.assertEqual(
            bad_paths,
            [],
            "No files or folders in the workspace may use 'lm2' or 'leafmachine2' in their name:\n"
            + "\n".join(bad_paths),
        )

    def test_no_active_files_reference_obsolete_sam2_lm2_or_pruned_helpers(self):
        """Scan all active .py, .R, and .sh files (excluding test_code_hygiene.py) and assert zero forbidden SAM-2, LM2, or pruned helper tokens."""
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            if src_path == SELF_PATH:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                for token in FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS:
                    if token in ("x_norm", "y_norm", "ingest_genetic_data"):
                        hit = bool(re.search(rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])", line))
                    elif token == "fig06c_":
                        hit = bool(re.search(r"fig06c_(?!cytonuclear_paralog_cytotype)", line))
                    else:
                        hit = token in line
                    if hit:
                        violations.append(f"{rel}:{lineno} [{token}]: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            "Active source files must not contain any obsolete SAM-2, LeafMachine2, or pruned helper tokens:\n"
            + "\n".join(violations),
        )

    def test_new_public_genomics_modules_respect_500_line_ceiling(self):
        """Verify public Packera genomics, cytology, taxonomy, geochemistry, monograph, and hygiene modules stay within line ceilings."""
        raw_ceilings = (
            ("scripts/data_prep/download_public_packera_sequences.py", 500), ("scripts/data_prep/step06_dryad_backbone_utils.py", 240),
            ("scripts/data_prep/ena_paraloss_discovery.py", 340), ("scripts/data_prep/download_public_cytology.py", 450),
            ("scripts/data_prep/step06_cytology_parse_utils.py", 220), ("scripts/analysis/step07_cytology_integration.R", 380),
            ("scripts/tests/test_public_cytology_integration.py", 420), ("scripts/data_prep/download_taxonomic_backbone.py", 450),
            ("scripts/data_prep/download_geochem_wetland_layers.py", 450), ("scripts/analysis/edaphic_geochem_wetland_utils.R", 340),
            ("scripts/tests/test_external_data_sources.py", 440), ("scripts/figures/step06c_figures.R", 350),
            ("scripts/analysis/step06_coalescent_models.py", 420), ("scripts/analysis/step06_model_selection.R", 380),
            ("scripts/analysis/step06_spatial_landscape_genomics.R", 380), ("scripts/tests/test_step06_model_selection.py", 420),
            ("scripts/morphometrics/gmm_icl_bootstrap_utils.R", 280), ("scripts/analysis/step05_circular_hopkins_models.R", 340),
            ("scripts/tests/test_step04_05_model_selection.py", 400), ("scripts/core/coord_hygiene_utils.py", 260),
            ("scripts/analysis/step05_env_outlier_utils.R", 280), ("scripts/analysis/step06_rogue_sequence_qc.R", 320),
            ("scripts/tests/test_outlier_quarantine_gates.py", 430), ("scripts/analysis/step06_public_backbone.R", 475),
            ("scripts/analysis/step06_organellar_containment.R", 200), ("scripts/analysis/step07_deme_propagation_utils.R", 260),
            ("scripts/analysis/step07_genomic_propagation.R", 400), ("scripts/tests/test_step07_deme_propagation.py", 300),
            ("scripts/data_prep/organellar_cytonuclear_utils.py", 300), ("scripts/data_prep/skim_offtarget_organellar.py", 460),
            ("scripts/analysis/step07_organellar_triage_utils.R", 260), ("scripts/tests/test_step07_organellar_triage.py", 280),
            ("scripts/analysis/step05_layer_hygiene_utils.R", 280), ("scripts/tests/test_step05_layer_hygiene.py", 380),
            ("scripts/core/step01_biogeography_utils.py", 300), ("scripts/core/step01_concept_authority_utils.py", 380),
            ("scripts/data_prep/prune_eastern_cohort.py", 500), ("scripts/tests/test_step01_robustness_audit.py", 380),
            ("scripts/tests/test_step01_taxonomic_clarity.py", 400), ("scripts/vision/step02_leaf_traits_utils.py", 360),
            ("scripts/vision/step02_photometry_corymb_utils.py", 360), ("scripts/vision/leaf_geometry_utils.py", 480),
            ("scripts/tests/test_step02_vision_hygiene.py", 380), ("scripts/tests/test_step02_taxonomic_clarity.py", 400),
            ("scripts/morphometrics/step03_config_scale_utils.R", 240), ("scripts/morphometrics/step03_ect_traits_utils.R", 360),
            ("scripts/morphometrics/03_fourier_extractor.R", 500), ("scripts/morphometrics/fourier_parallel_utils.R", 500),
            ("scripts/morphometrics/plasticity_mitigation_utils.R", 470), ("scripts/figures/step03_figures.R", 385),
            ("scripts/tests/test_step03_audit_contracts.py", 360), ("scripts/tests/test_step03_taxonomic_clarity.py", 380),
            ("scripts/morphometrics/04_gmm_morphotools.R", 500), ("scripts/morphometrics/step04_allometry_merger_utils.R", 380),
            ("scripts/morphometrics/step04_trait_projection_utils.R", 360), ("scripts/morphometrics/pub_synthesis_plates.R", 120),
            ("scripts/figures/step04_figures.R", 500), ("scripts/tests/test_step04_audit_hygiene.py", 380),
            ("scripts/tests/test_step04_taxonomic_clarity.py", 380), ("scripts/analysis/05_multimodal_spatial_rf.R", 500),
            ("scripts/analysis/step05_env_raster_extract_utils.R", 420), ("scripts/analysis/step05_spatial_rf_cv_utils.R", 500),
            ("scripts/analysis/step05_consensus_niche_utils.R", 500), ("scripts/analysis/edaphic_fine_grain_utils.R", 500),
            ("scripts/analysis/edaphic_ssurgo_sda_utils.R", 480), ("scripts/analysis/edaphic_solus_3dep_utils.R", 480),
            ("scripts/analysis/edaphic_collinearity_contingency_utils.R", 480), ("scripts/analysis/step05_phenology_ecotypes.R", 500),
            ("scripts/analysis/step05_allochronic_pheno_utils.R", 480), ("scripts/analysis/step05_type_niche_utils.R", 480),
            ("scripts/analysis/step05_config_geo_utils.R", 320), ("scripts/analysis/step05_sdm_maxent.R", 500),
            ("scripts/analysis/step05_sdm_engine_utils.R", 500), ("scripts/data_prep/download_environmental_layers.py", 480),
            ("scripts/data_prep/env_sgmc_3dep_utils.py", 480), ("scripts/figures/step05_figures.R", 500),
            ("scripts/figures/step05_sdm_figures.R", 320), ("scripts/tests/test_step05_audit_hygiene.py", 380),
            ("scripts/tests/test_step05_taxonomic_clarity.py", 380), ("scripts/tests/test_multimodal_spatial_rf.py", 500),
            ("scripts/tests/test_step05_edaphic_niche.py", 500), ("scripts/tests/test_step05_zero_trust_enm.py", 500),
            ("scripts/data_prep/06_ingest_genetic_data.py", 500), ("scripts/data_prep/step06_synonym_utils.py", 200),
            ("scripts/data_prep/step06_sample_metadata_utils.py", 460), ("scripts/data_prep/step06_voucher_linker_utils.py", 460),
            ("scripts/data_prep/step06_vcf_dosage_utils.py", 460), ("scripts/data_prep/step06_ploidy_hdplot_utils.py", 460),
            ("scripts/data_prep/step06_public_cohort_ingestor.py", 450), ("scripts/data_prep/step06_public_tree_qc_utils.py", 200),
            ("scripts/analysis/06_genomic_phylogenomics.R", 500), ("scripts/analysis/step06_dynamic_tier_engine.R", 470),
            ("scripts/analysis/step06_phylogenomics_utils.R", 480), ("scripts/analysis/step06_coalescent_delimitation.py", 500),
            ("scripts/analysis/step06_coalescent_alignment_utils.py", 480), ("scripts/analysis/step06_reticulation_coancestry.R", 500),
            ("scripts/analysis/step06_fbranch_quartets.R", 380), ("scripts/figures/step06_figures.R", 480),
            ("scripts/figures/step06_tier_figures.R", 480), ("scripts/tests/test_step06_audit_hygiene.py", 400),
            ("scripts/tests/test_step06_dynamic_tiers.py", 420), ("scripts/tests/test_step06_ingestion_hdplot_ploidy.py", 460),
            ("scripts/tests/test_step06_organellar_coalescent.py", 460), ("scripts/tests/test_step06_phylogenomics_r_pipeline.py", 460),
            ("scripts/tests/test_step06_public_backbone_dryad.py", 460), ("scripts/tests/test_step06_voucher_exsiccatae_linking.py", 460),
            ("scripts/tests/test_step07_genomic_propagation_contract.py", 460), ("scripts/analysis/07_triage_dashboard_synthesis.R", 500),
            ("scripts/analysis/step07_hopkins_consensus_utils.R", 480), ("scripts/analysis/step07_ancova_merge_utils.R", 490),
            ("scripts/analysis/step07_reports_html_utils.R", 490), ("scripts/analysis/step07_biogeo_hybrid_utils.R", 490),
            ("scripts/analysis/step07_cluster_cohesion_utils.R", 460), ("scripts/analysis/step07_decision_stages_utils.R", 490),
            ("scripts/analysis/step07_postprocess_guards_utils.R", 460), ("scripts/analysis/taxon_complex_discovery_utils.R", 480),
            ("scripts/analysis/step07_complex_affinity_utils.R", 450), ("scripts/figures/step07_figures.R", 420),
            ("scripts/figures/step07b_network_figures.R", 320), ("scripts/tests/test_step07_audit_hygiene.py", 400),
            ("scripts/tests/test_parallel_thread_scaling.py", 300),
        )
        strict_ceilings = {PROJECT_ROOT / rel: lim for rel, lim in raw_ceilings}
        bounded_paths = (
            *strict_ceilings.keys(),
            PROJECT_ROOT / "scripts" / "analysis" / "genomic_triage_utils.R",
            PROJECT_ROOT / "scripts" / "pipeline" / "resolve_eastern_na_packera.py",
            PROJECT_ROOT / "scripts" / "tests" / "test_eastern_na_monograph_cli.py",
        )
        for mod_path in bounded_paths:
            self.assertTrue(mod_path.exists(), f"Missing expected bounded module: {mod_path}")
            line_count = len(mod_path.read_text(encoding="utf-8").splitlines())
            self.assertLess(line_count, 500, f"{mod_path.relative_to(PROJECT_ROOT)} exceeds 500-line ceiling ({line_count} lines)")
            if mod_path in strict_ceilings:
                ceiling = strict_ceilings[mod_path]
                self.assertLess(line_count, ceiling, f"{mod_path.relative_to(PROJECT_ROOT)} exceeds {ceiling}-line ceiling ({line_count} lines)")

    def test_step06_split_posix_and_mmrr_hygiene(self):
        """Assert legacy test monolith removal, pure POSIX regexes across Step 06 R files, and zero geographic abs(lat) IBE forgery."""
        legacy_monolith = PROJECT_ROOT / "scripts" / "tests" / "test_genomic_integration.py"
        self.assertFalse(legacy_monolith.exists(), "Legacy 2,426-line scripts/tests/test_genomic_integration.py must remain deleted")
        perl_re = re.compile(r"\\\\[sdwb]|\(\?:")
        step06_r_files = sorted((PROJECT_ROOT / "scripts" / "analysis").glob("*06*.R")) + sorted(
            (PROJECT_ROOT / "scripts" / "figures").glob("*06*.R")
        )
        self.assertTrue(step06_r_files, "Expected Step 06 R files to audit")
        perl_hits: List[str] = []
        for r_path in step06_r_files:
            rel = r_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(r_path.read_text(encoding="utf-8").splitlines(), start=1):
                if perl_re.search(line):
                    perl_hits.append(f"{rel}:{lineno}: {line.strip()}")
        self.assertEqual(perl_hits, [], "Step 06 R files must use pure POSIX ERE instead of Perl escapes:\n" + "\n".join(perl_hits))
        mmrr_text = (PROJECT_ROOT / "scripts" / "analysis" / "step06_spatial_landscape_genomics.R").read_text(encoding="utf-8")
        self.assertNotIn("abs(lat", mmrr_text, "step06_spatial_landscape_genomics.R must never substitute abs(lat into eco_mat")

    def test_data_cytology_directory_stays_flat(self):
        """Verify `data/cytology/` exists and contains zero nested subdirectories."""
        cytology_dir = PROJECT_ROOT / "data" / "cytology"
        self.assertTrue(cytology_dir.is_dir(), f"Missing expected cytology directory: {cytology_dir}")
        subdirs = [p.name for p in cytology_dir.iterdir() if p.is_dir()]
        self.assertEqual(subdirs, [], f"data/cytology/ must remain flat with zero subdirectories: {subdirs}")

    def test_data_taxonomy_directory_stays_flat(self):
        """Verify `data/taxonomy/` contains zero nested subdirectories when present."""
        taxonomy_dir = PROJECT_ROOT / "data" / "taxonomy"
        if taxonomy_dir.exists():
            self.assertTrue(taxonomy_dir.is_dir(), f"Expected directory at {taxonomy_dir}")
            subdirs = [p.name for p in taxonomy_dir.iterdir() if p.is_dir()]
            self.assertEqual(subdirs, [], f"data/taxonomy/ must remain flat with zero subdirectories: {subdirs}")

    def test_no_hardcoded_workstation_legacy_repo_path(self):
        """Scan all active .py, .R, and .sh files and assert zero hardcoded workstation legacy repo paths."""
        legacy_abs_path = "/" + "home/brandon/" + "Packera_dubia_morphometrics"
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            if src_path == SELF_PATH:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                if legacy_abs_path in line:
                    violations.append(f"{rel}:{lineno}: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            "Active source files must not hardcode the legacy workstation repository path:\n"
            + "\n".join(violations),
        )

    def test_no_duplicate_python_defs_in_same_scope(self):
        """Assert zero shadowed or duplicate def/async def/class definitions within any module or class body."""
        def _is_overload(node: ast.AST) -> bool:
            for dec in getattr(node, "decorator_list", []):
                if isinstance(dec, ast.Name) and dec.id == "overload":
                    return True
                if isinstance(dec, ast.Attribute) and dec.attr == "overload":
                    return True
            return False

        def _check_body(body: List[ast.stmt], scope_label: str, dupes: List[str]) -> None:
            seen: dict[str, int] = {}
            for stmt in body:
                if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    if _is_overload(stmt):
                        continue
                    if stmt.name in seen:
                        dupes.append(f"{scope_label}:{stmt.lineno} duplicates {stmt.name} (first at L{seen[stmt.name]})")
                    else:
                        seen[stmt.name] = stmt.lineno
                    if isinstance(stmt, ast.ClassDef):
                        _check_body(stmt.body, f"{scope_label}.{stmt.name}", dupes)

        dupes: List[str] = []
        py_files = [PROJECT_ROOT / "main.py"] + sorted(
            p for p in (PROJECT_ROOT / "scripts").rglob("*.py")
            if ARCHIVE_MARKER not in p.parts
        )
        for py_path in py_files:
            if not py_path.exists():
                continue
            tree = ast.parse(py_path.read_text(encoding="utf-8"), filename=str(py_path))
            _check_body(tree.body, str(py_path.relative_to(PROJECT_ROOT)), dupes)

        self.assertEqual(dupes, [], "Duplicate Python definitions found in same scope:\n" + "\n".join(dupes))

    def test_no_broken_symlinks(self):
        """Assert zero broken symlinks across models/, scripts/, and runs/latest*."""
        broken: List[str] = []
        for root_dir in (PROJECT_ROOT / "models", PROJECT_ROOT / "scripts"):
            if root_dir.exists():
                for entry in root_dir.rglob("*"):
                    if "__pycache__" in entry.parts:
                        continue
                    if entry.is_symlink() and not entry.exists():
                        broken.append(str(entry.relative_to(PROJECT_ROOT)))
        runs_dir = PROJECT_ROOT / "runs"
        if runs_dir.exists():
            for entry in runs_dir.glob("latest*"):
                if entry.is_symlink() and not entry.exists():
                    broken.append(str(entry.relative_to(PROJECT_ROOT)))
        self.assertEqual(broken, [], "Broken symlinks detected:\n" + "\n".join(broken))

    def test_no_pdf_figure_artifacts(self):
        """Assert zero Rplots.pdf or figure .pdf artifacts at repo root or figure directories."""
        found_pdfs: List[str] = []
        rplots = PROJECT_ROOT / "Rplots.pdf"
        if rplots.exists() or rplots.is_symlink():
            found_pdfs.append("Rplots.pdf")
        out_figs = PROJECT_ROOT / "outputs" / "figures"
        if out_figs.exists():
            for pdf_p in out_figs.rglob("*.pdf"):
                found_pdfs.append(str(pdf_p.relative_to(PROJECT_ROOT)))
        runs_dir = PROJECT_ROOT / "runs"
        if runs_dir.exists():
            for fig_dir in runs_dir.glob("*/outputs/figures"):
                if fig_dir.is_symlink():
                    continue
                for pdf_p in fig_dir.rglob("*.pdf"):
                    found_pdfs.append(str(pdf_p.relative_to(PROJECT_ROOT)))
        self.assertEqual(found_pdfs, [], "Unexpected PDF figure artifacts found:\n" + "\n".join(sorted(found_pdfs)))


if __name__ == "__main__":
    unittest.main()
