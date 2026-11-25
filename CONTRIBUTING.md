# Contributing and Developer Guide (`Plantmatics`)

**Standard Operating Procedure:** [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md) (`UNC-BOT-SOP-2026-04-REV5`)
**Architectural Contract:** [`AGENTS.md`](AGENTS.md) | **Schema Reference:** [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md)

I welcome bug fixes, botanical authority-table additions, regional environmental adapters, and statistical improvements that preserve the mathematical and curatorial rigor of Plantmatics. Because Python harvesting and ONNX computer vision stages hand off directly to R morphometric and spatial niche models through a shared SQLite Write-Ahead Logging (`WAL`) store (`data/plantmatics.sqlite`), every pull request must satisfy the architectural boundaries and automated audit checks below before merge.

## 1. Dual Virtual Environment and `LeafMachine3` Submodule Boundary

To keep compiled geospatial C++ libraries (`GDAL`, `GEOS`, `PROJ`) from conflicting with pinned ONNX Runtime and `CUDA` wheels, the repository maintains two strictly separated Python environments:

* **`.venv` (Primary Orchestration & Geospatial Python):** Created from `requirements.txt`. Runs `main.py`, `run_pipeline.sh` preflight helpers, Step 01 GBIF/iDigBio harvesting, EPA ecoregion and USDA county range lookups, `download_environmental_layers.py`, SQLite/DuckDB migrations, and the Python test runner.
* **`LeafMachine3/.venv` (`.venv_LM3` Vision Isolation):** Provisioned automatically via `uv` by `bash setup_leafmachine3.sh`. Executes Step 02 (`02_segment_and_extract.py`) ONNX inference, native `LeafMachine3` visual rendering, and $32 \times 360$ ECT calculation, writing to SQLite strictly through Python's standard-library `sqlite3` module.

Never modify, patch, or commit changes inside `LeafMachine3/`. All `Packera`-specific computer vision adaptations belong in our read-only adapter layer under `scripts/vision/` and `scripts/pipeline/02_segment_and_extract.py`, and all working-frame coordinates (`coord_work`) must divide by `work_scale` before entering SQLite or CSV outputs.

## 2. Module Ownership and Universal `< 500`-Line POSIX Hygiene

Keep functions cohesive within their designated directories (`scripts/core/`, `scripts/data_prep/`, `scripts/vision/`, `scripts/pipeline/`, `scripts/morphometrics/`, `scripts/figures/`, and `scripts/analysis/`), and never import or `source()` anything from `_archive` folders.

Every active Python and R module alongside every shell script obeys a universal `< 500`-line ceiling enforced automatically by `scripts/tests/test_code_hygiene.py`, paired with tighter per-module budgets in `scripts/tests/_hygiene_specs.py` and a repository-wide pure POSIX Extended Regular Expression gate across 100% of active `.R` files (zero Perl escapes such as `\s`, `\d`, `\w`, `\b`, or `(?:`):

| Entrypoint / Facade | Facade Ceiling | Modular Companion Files (All `< 500` Lines) |
| :--- | :---: | :--- |
| `main.py` | `< 360` lines | `scripts/pipeline/cli_*.py` (`cli_run_context.py`, `cli_env_doctor.py`, `cli_stage_*.py`, `cli_parser_subcommands_*.py`) |
| `run_pipeline.sh` | `< 460` lines | `scripts/pipeline/sh/pipe_*.sh` (`pipe_cli_args.sh`, `pipe_env_purge.sh`, `pipe_preflight_cohort_*.sh`, `pipe_steps_01_04.sh`) |
| `scripts/core/config.py` | `< 240` lines | `config_run_context.py`, `config_dataclasses.py`, `config_pipeline_cls.py`, `config_resolvers.py` |
| `scripts/core/db_schema.py` & `db.py` | `< 240` / `< 250` lines | Eight `db_schema_*.py` DDL/migration modules paired with thirteen `db_*.py` persistence/purge modules |
| `scripts/core/db_utils.R` | `< 150` lines | Thirteen auto-sourced `db_r_*.R` SQLite/DuckDB and stage-persistence modules |
| `scripts/morphometrics/03_fourier_extractor.R` | `< 500` lines | `step03_config_scale_utils.R`, `step03_ect_traits_utils.R`, `dissected_leaf_utils.R` (`dissect_*.R`), `fourier_parallel_utils.R`, `scripts/figures/step03_figures.R` |
| `scripts/morphometrics/04_gmm_morphotools.R` | `< 500` lines | `step04_allometry_merger_utils.R`, `step04_trait_projection_utils.R`, `gmm_icl_bootstrap_utils.R`, `authority_weighting_utils.R` (`auth_*.R`), `plasticity_mitigation_utils.R`, `pub_ecoregion_manova_utils.R` (`ecoregion_*.R`), `pub_figures_utils.R` (`pub_taxon_label_utils.R`), `scripts/figures/step04_figures.R` |
| `scripts/figures/theme_packera_pub.R` | `< 400` lines | `scripts/figures/basemap_exemplar_utils.R` (`< 240`), `scripts/figures/step03_figures.R` to `step07b_network_figures.R` |
| `scripts/analysis/05_multimodal_spatial_rf.R` | `< 500` lines | `step05_*_utils.R`, `edaphic_*_utils.R`, `step05_sdm_maxent.R`, `scripts/figures/step05*_figures.R` |
| `scripts/analysis/06_genomic_phylogenomics.R` | `< 500` lines | `step06_*.R`, `step06_coalescent_*.py`, `scripts/data_prep/06_ingest_genetic_data.py`, `scripts/data_prep/step06_*.py`, `scripts/figures/step06*_figures.R` |
| `scripts/analysis/07_triage_dashboard_synthesis.R` | `< 500` lines | `step07_*_utils.R`, `step07_cytology_integration.R`, `step07_genomic_propagation.R`, `genomic_triage_utils.R`, `taxon_complex_discovery_utils.R`, `scripts/figures/step07*_figures.R` |
| `scripts/pipeline/resolve_focal_species_complex.py` | `< 340` lines | `step07_focal_cohort_utils.py` (`< 450`), `step07_focal_stage_runners.py` (`< 450`), `resolve_eastern_na_packera.py` (`< 500`) |

When adding new morphometric diagnostics or plots, place shared plotting logic in `scripts/figures/theme_packera_pub.R` paired with `scripts/figures/basemap_exemplar_utils.R` or inside `scripts/figures/step0X_figures.R` so core entrypoints stay under their line budgets.

## 3. SQLite Schema (`20261127`), Composite Keys, and Demo Isolation

1. **Synchronized DDL Across Python and R:** Whenever you add or rename a table column or modify `v_integrative_triage_master`, update `scripts/core/db_schema_*.py` (re-exported by `db_schema.py`) alongside `scripts/core/db_*.py` (re-exported by `db.py`) plus `scripts/core/db_r_*.R` (auto-sourced by `db_utils.R`) and [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) in the same commit, and increment `PLANTMATICS_SCHEMA_USER_VERSION` (preserving `PACKERA_SCHEMA_USER_VERSION` alias compatibility) if breaking migrations are required.
2. **Composite Individual Keys:** Group and join all specimen-level data on `(catalogNumber, plant_individual_id)` (plus `run_id` in Tier B tables). Never deduplicate by `catalogNumber` alone in Steps 02 through 07, as multi-plant herbarium sheets carry distinct DBSCAN individuals (`plant_individual_id = 1, 2, ...`).
3. **Atomic Writes and Honest Missingness:** Write CSV, HTML, Markdown, and GeoTIFF deliverables to process-unique `.tmp` or `.part` files first and `rename()` atomically upon completion. Outside explicit `--demo` runs (`PLANTMATICS_DEMO_MODE=1`, with `PACKERA_DEMO_MODE=1` compatibility, which isolates writes inside `demo_ephemeral.sqlite`), never synthesize mock coordinates, random embeddings, or fabricated raster values when external layers are absent; return `NA` / `NULL` cleanly.
4. **RNG State Hygiene in R:** Any stochastic call inside helper routines (`igraph::cluster_leiden`, `kmeans`, bootstrap/permutation loops) must preserve and restore the caller's `.Random.seed` via `on.exit` (see `.with_preserved_rng_seed()` in `taxon_complex_discovery_utils.R`).

## 4. 300-DPI PNG-Only Figure Protocol

All pipeline visualizations export strictly as **300-DPI PNG** files (`0` `.pdf` outputs) through `save_packera_png()` in `scripts/figures/theme_packera_pub.R` (or Matplotlib `dpi=300` in `scripts/vision/step02_figures.py`):

* **Dynamic Taxon and Variable Scaling:** Figures must handle $N = 1$ intraspecific runs (falling back to unsupervised GMM morphotype clusters), $G = 2$ pairwise comparisons (1D `Can1` density + raincloud rather than degenerate 2D axes), and $60+$ species or infraspecific ranks (`subsp.`, `var.`, `f.`) using `build_packera_taxon_aesthetic_map()` (12 colorblind-safe hues $\times$ 5 unfilled hollow shapes) and dynamic legend canvas heights.
* **Botanical & Cartographic Standards:** Format all taxon labels through `format_packera_taxon_markdown()` / ICN expression labellers so genus and species/infraspecific epithets appear in italics while rank markers (`subsp.`, `var.`, `f.`, `s.s.`, `s.l.`, `sp. nov.`) remain Roman. Project all geographic maps to equal-area `EPSG:5070` (`NAD83 / Conus Albers`).

## 5. Documentation and Prose Style (`Writing DNA v5.2`)

All Markdown documentation (`README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, and `docs/*.md`) follows the prose rules in [`.agents/skills/write-in-my-voice/SKILL.md`](.agents/skills/write-in-my-voice/SKILL.md):

* Write in clear, active first-person or direct imperative sentences with varied paragraph lengths and explicit causal or contrastive connectors (`, as ...`, `, whereas ...`, `, so ...`).
* Enforce **zero em dashes**, **zero en dashes**, and **zero prose double-hyphens** (write numerical ranges as `"0 to 15 cm"` or `"Steps 01 to 07"`). CLI flags inside code spans (`--focal-taxon`) are fine.
* Cite academic literature using clean `doi: 10.xxxx/...` identifiers rather than raw publisher URLs, and describe active pipeline architecture in positive present-state terms.

## 6. Running the Pre-Commit Verification Suite

Before opening a pull request or generating a patch, run the full verification gate from the repository root:

```bash
# 1a. Run the complete Python, SQLite handshake, warm-R statistical, figure, and hygiene suite
.venv/bin/python scripts/tests/run_audit.py

# 1b. Optional rapid Python-only pre-check during active iteration (~25s)
.venv/bin/python scripts/tests/run_audit.py --fast

# 2. Verify shell syntax across all orchestrator and setup scripts
bash -n run_pipeline.sh scripts/pipeline/sh/*.sh setup_leafmachine3.sh resolve_packera_dubia_complex.sh resolve_eastern_na_packera.sh

# 3. Verify zero em/en dashes across all Markdown files alongside skill and citation files
grep -rnP "[\x{2013}\x{2014}]" README.md AGENTS.md CONTRIBUTING.md CHANGELOG.md CITATION.cff .agents/skills/pipeline-contract-audit/SKILL.md docs/*.md && exit 1 || echo "PASS: Zero em/en dashes"

# 4. Verify universal <500-line ceiling and pure POSIX R hygiene explicitly
.venv/bin/python -m unittest scripts/tests/test_code_hygiene.py -v

# 5. Verify git whitespace cleanliness
git diff --check
```
