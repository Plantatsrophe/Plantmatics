"""Tier A canonical shared SQLite table, trigger, R*Tree, and index DDL."""

from __future__ import annotations

from typing import List
from scripts.core.db_schema_env_cols import ENVIRONMENTAL_BOOL_COLS, ENVIRONMENTAL_CACHE_COLUMNS, ENVIRONMENTAL_TEXT_COLS

def _build_tier_a_schema_ddl() -> str:
    env_col_defs: List[str] = ["lat_lon_key TEXT PRIMARY KEY"]
    for col in ENVIRONMENTAL_CACHE_COLUMNS:
        if col == "lat_lon_key":
            continue
        if col in ENVIRONMENTAL_TEXT_COLS:
            env_col_defs.append(f"{col} TEXT")
        elif col in ENVIRONMENTAL_BOOL_COLS:
            env_col_defs.append(f"{col} INTEGER")
        else:
            env_col_defs.append(f"{col} REAL")
    env_ddl = ",\n        ".join(env_col_defs)
    return f"""
    -- =========================================================================
    -- TIER A: CANONICAL SHARED STORE
    -- =========================================================================

    CREATE TABLE IF NOT EXISTS taxon_synonyms (
        raw_name_lower TEXT PRIMARY KEY,
        canonical_packera TEXT NOT NULL,
        protologue_authors_json TEXT,
        updated_at TEXT
    );

    CREATE TABLE IF NOT EXISTS usda_county_ranges (
        canonical_packera TEXT NOT NULL,
        usda_symbol TEXT NOT NULL,
        state_norm TEXT NOT NULL,
        county_norm TEXT NOT NULL DEFAULT '',
        fips_code TEXT,
        native_status TEXT,
        is_state_only INTEGER NOT NULL DEFAULT 0,
        updated_at TEXT,
        PRIMARY KEY (canonical_packera, state_norm, county_norm)
    );

    CREATE INDEX IF NOT EXISTS idx_usda_county_ranges_taxon ON usda_county_ranges(canonical_packera);

    CREATE TABLE IF NOT EXISTS vouchers (
        catalogNumber TEXT PRIMARY KEY,
        institutionCode TEXT,
        scientificName TEXT,
        species_raw TEXT,
        recordedBy TEXT,
        recordNumber TEXT,
        identifiedBy TEXT,
        determiner_tier TEXT,
        type_status TEXT,
        county TEXT,
        stateProvince TEXT,
        locality TEXT,
        decimalLatitude REAL,
        decimalLongitude REAL,
        coordinateUncertainty REAL,
        year INTEGER,
        month INTEGER,
        day INTEGER,
        eventDate TEXT,
        doy REAL,
        regional_group TEXT,
        image_path TEXT,
        exsiccatae_key TEXT,
        is_primary_duplicate INTEGER,
        duplicate_count INTEGER,
        exsiccatae_taxonomic_conflict INTEGER DEFAULT 0,
        exsiccatae_gold_taxon TEXT,
        usda_symbol TEXT,
        usda_county_status TEXT,
        infraspecific_taxon TEXT,
        georef_qc_status TEXT,
        is_centroid_or_inst_coord INTEGER DEFAULT 0,
        coord_quarantine_reason TEXT,
        typified_basionym TEXT,
        verbatim_taxon_concept TEXT,
        concept_synonymy_status TEXT,
        identification_year INTEGER,
        is_open_or_hybrid_concept INTEGER DEFAULT 0,
        spatial_balance_weight REAL,
        species_standardized TEXT,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        qc_status TEXT DEFAULT 'PENDING' CHECK(qc_status IN ('PASS', 'FAILED_QC', 'PENDING')),
        failure_reason TEXT,
        failure_details TEXT
    );

    CREATE INDEX IF NOT EXISTS idx_vouchers_species_std ON vouchers(species_standardized);
    CREATE INDEX IF NOT EXISTS idx_vouchers_qc_status ON vouchers(qc_status);
    CREATE INDEX IF NOT EXISTS idx_vouchers_det_tier ON vouchers(determiner_tier);

    CREATE TABLE IF NOT EXISTS plant_individuals (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        rosette_centroid_x REAL,
        rosette_centroid_y REAL,
        crown_x REAL,
        crown_y REAL,
        plant_height_px REAL,
        plant_height_mm REAL,
        plant_height_cm REAL,
        is_bent_stem INTEGER,
        petiole_length_px REAL,
        petiole_length_mm REAL,
        petiole_width_px REAL,
        petiole_width_mm REAL,
        blade_length_px REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        scape_to_basal_leaf_ratio REAL,
        rosette_patch_path TEXT,
        rosette_mask_path TEXT,
        PRIMARY KEY (catalogNumber, plant_individual_id),
        FOREIGN KEY (catalogNumber) REFERENCES vouchers(catalogNumber) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS extracted_leaves (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        leaf_id INTEGER NOT NULL,
        shape_id TEXT UNIQUE,
        solidity REAL,
        circularity REAL,
        extent REAL,
        aspect_ratio REAL,
        area_px REAL,
        area_mm2 REAL,
        midrib_angle_deg REAL,
        petiole_length_px REAL,
        petiole_length_mm REAL,
        petiole_width_px REAL,
        petiole_width_mm REAL,
        blade_length_px REAL,
        blade_length_mm REAL,
        blade_width_px REAL,
        blade_width_mm REAL,
        petiole_to_blade_ratio REAL,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        scale_source TEXT,
        mask_path TEXT,
        mask_is_upright INTEGER DEFAULT 1,
        crop_path TEXT,
        contour_path TEXT,
        ect_path TEXT,
        contour_n_pts INTEGER DEFAULT 120,
        contour_xy_json TEXT NOT NULL,
        leaf_zone TEXT DEFAULT 'basal',
        leaf_zone_confidence REAL,
        dist_from_crown_norm REAL,
        cy_rel_plant_bbox REAL,
        archetype_score REAL,
        bilateral_dice REAL,
        si_a REAL,
        lm3_gates_pass INTEGER,
        n_lamina_components INTEGER DEFAULT 1,
        detector_class_id INTEGER DEFAULT 0,
        is_crop_clipped INTEGER DEFAULT 0,
        hole_area_frac REAL DEFAULT 0.0,
        enclosed_void_frac REAL DEFAULT 0.0,
        bilateral_lobe_symmetry INTEGER DEFAULT 0,
        efa_qc_pass INTEGER DEFAULT 1,
        qc_flag_reason TEXT DEFAULT 'PASS',
        blade_cie_l_norm REAL,
        blade_cie_l_sd REAL,
        blade_cie_l_sheet_delta REAL,
        blade_cie_a_delta REAL,
        blade_cie_b_delta REAL,
        bbox_x1 INTEGER,
        bbox_y1 INTEGER,
        bbox_x2 INTEGER,
        bbox_y2 INTEGER,
        base_angle_deg REAL,
        apex_angle_deg REAL,
        lamina_curvature REAL,
        lamina_trace_length_px REAL,
        basal_cordate_index REAL,
        margin_serration_index REAL,
        margin_tooth_acuteness REAL,
        terminal_lobe_area_frac REAL,
        blade_cie_l_q85_sheet_delta REAL,
        basal_zone_cie_l_delta REAL,
        is_abaxial_pressed_leaf INTEGER DEFAULT 0,
        abaxial_cie_l_sheet_delta REAL,
        PRIMARY KEY (catalogNumber, plant_individual_id, leaf_id),
        FOREIGN KEY (catalogNumber, plant_individual_id)
            REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
    );

    CREATE INDEX IF NOT EXISTS idx_extracted_leaves_cat ON extracted_leaves(catalogNumber);

    CREATE TABLE IF NOT EXISTS voucher_reproductive_phenology (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        capitula_count INTEGER,
        involucre_height_px REAL,
        involucre_width_px REAL,
        involucre_height_mm REAL,
        involucre_width_mm REAL,
        AR_inv REAL,
        total_capitula INTEGER,
        n_anthesis INTEGER,
        n_bud INTEGER,
        n_fruiting INTEGER,
        n_senescent INTEGER,
        voucher_phenological_state TEXT,
        is_verified_anthesis INTEGER,
        reproductive_stage_confidence REAL,
        dominant_pappus_ratio REAL,
        corymb_head_count_est INTEGER DEFAULT 0,
        radiate_head_fraction REAL,
        PRIMARY KEY (catalogNumber, plant_individual_id),
        FOREIGN KEY (catalogNumber, plant_individual_id)
            REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS environmental_cache (
        {env_ddl}
    );

    CREATE INDEX IF NOT EXISTS idx_environmental_cache_catalog
        ON environmental_cache(catalogNumber);

    CREATE VIRTUAL TABLE IF NOT EXISTS environmental_rtree USING rtree(
        id,
        min_lon, max_lon,
        min_lat, max_lat
    );

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_insert
    AFTER INSERT ON environmental_cache
    WHEN NEW.latitude IS NOT NULL AND NEW.longitude IS NOT NULL
    BEGIN
        INSERT OR REPLACE INTO environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)
        VALUES (NEW.rowid, NEW.longitude, NEW.longitude, NEW.latitude, NEW.latitude);
    END;

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_update
    AFTER UPDATE ON environmental_cache
    BEGIN
        DELETE FROM environmental_rtree WHERE id = OLD.rowid;
        INSERT INTO environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)
        SELECT NEW.rowid, NEW.longitude, NEW.longitude, NEW.latitude, NEW.latitude
        WHERE NEW.latitude IS NOT NULL AND NEW.longitude IS NOT NULL;
    END;

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_delete
    AFTER DELETE ON environmental_cache
    BEGIN
        DELETE FROM environmental_rtree WHERE id = OLD.rowid;
    END;

    CREATE TABLE IF NOT EXISTS genomic_samples (
        genetic_sample_id TEXT PRIMARY KEY,
        catalogNumber TEXT NULL,
        project_id TEXT,
        population_id TEXT,
        replicate_pair_id TEXT,
        taxon_raw TEXT,
        species_standardized TEXT,
        is_putative_novel_or_hybrid_query INTEGER DEFAULT 0,
        ploidy_level INTEGER,
        marker_platform TEXT,
        decimalLatitude REAL,
        decimalLongitude REAL,
        stateProvince TEXT,
        county TEXT,
        eventDate TEXT,
        hq_reads INTEGER,
        mapped_reads INTEGER,
        pct_repetitive REAL,
        n_called_loci INTEGER,
        missing_genotype_rate REAL,
        observed_het REAL,
        geno_quality_pass INTEGER,
        linkage_tier TEXT,
        hdplot_retained_loci INTEGER,
        inferred_cytotype TEXT,
        nquire_delta_logl REAL,
        plastid_clade TEXT,
        its_ribotype TEXT,
        cytonuclear_discordant INTEGER,
        public_plastome_match TEXT,
        public_its_match TEXT,
        waxy_lineage_match TEXT,
        exsiccatae_key TEXT,
        genbank_accession TEXT,
        coordinateUncertainty REAL,
        is_centroid_or_inst_coord INTEGER DEFAULT 0,
        coord_quarantine_reason TEXT,
        typified_basionym TEXT,
        concept_synonymy_status TEXT,
        is_open_or_hybrid_concept INTEGER DEFAULT 0
    );
    CREATE INDEX IF NOT EXISTS idx_genomic_samples_catalog ON genomic_samples(catalogNumber);
    CREATE INDEX IF NOT EXISTS idx_genomic_samples_pop ON genomic_samples(population_id);
    CREATE INDEX IF NOT EXISTS idx_genomic_samples_taxon ON genomic_samples(species_standardized);

    -- 8. Cytological records (canonical chromosome counts, C-values, breeding systems, voucher links)
    CREATE TABLE IF NOT EXISTS cytological_records (
        cytology_id TEXT PRIMARY KEY,
        source_db TEXT NOT NULL,
        source_record_id TEXT,
        original_taxon_name TEXT NOT NULL,
        basionym_name TEXT,
        species_standardized TEXT NOT NULL,
        catalogNumber TEXT,
        exsiccatae_key TEXT,
        cytology_linkage_tier TEXT,
        raw_count_str TEXT,
        gametophytic_n INTEGER,
        sporophytic_2n INTEGER,
        base_x INTEGER,
        ploidy_level INTEGER,
        cytotype_label TEXT,
        b_chromosomes INTEGER DEFAULT 0,
        is_aneuploid_or_odd_ploid INTEGER DEFAULT 0,
        is_approximate INTEGER DEFAULT 0,
        c_value_1c_pg REAL,
        c_value_2c_pg REAL,
        breeding_system TEXT,
        genbank_accession TEXT,
        collector_surname TEXT,
        record_number_clean TEXT,
        institution_code TEXT,
        voucher_citation TEXT,
        locality_text TEXT,
        reference_citation TEXT
    );
    CREATE INDEX IF NOT EXISTS idx_cytological_records_taxon ON cytological_records(species_standardized);
    CREATE INDEX IF NOT EXISTS idx_cytological_records_catalog ON cytological_records(catalogNumber);
    CREATE INDEX IF NOT EXISTS idx_cytological_records_exsiccatae ON cytological_records(exsiccatae_key);

    -- 9. Taxon-level cytological summary (aggregated modal ploidy, series, and C-value statistics)
    CREATE TABLE IF NOT EXISTS taxon_cytology_summary (
        species_standardized TEXT PRIMARY KEY,
        n_cytology_records INTEGER DEFAULT 0,
        n_direct_cytovouchers INTEGER DEFAULT 0,
        modal_sporophytic_2n INTEGER,
        modal_base_x INTEGER,
        modal_ploidy_level INTEGER,
        reported_cytotypes TEXT,
        is_multi_cytotype_taxon INTEGER DEFAULT 0,
        has_dysploid_series INTEGER DEFAULT 0,
        has_odd_polyploid_or_apomixis INTEGER DEFAULT 0,
        mean_2c_pg REAL,
        breeding_system_summary TEXT
    );

    -- 9. Taxonomic & Protologue Backbone (CoL/POWO/BHL + NatureServe Global Rank + NWPL Indicator)
    CREATE TABLE IF NOT EXISTS taxon_backbone_profiles (
        species_standardized TEXT PRIMARY KEY,
        target_genus TEXT NOT NULL DEFAULT 'Packera',
        accepted_scientific_name TEXT,
        accepted_authorship TEXT,
        basionym_name TEXT,
        basionym_authorship TEXT,
        homotypic_synonyms_json TEXT,
        heterotypic_synonyms_json TEXT,
        is_hybrid_taxon INTEGER DEFAULT 0,
        hybrid_formula TEXT,
        publication_citation TEXT,
        publication_year INTEGER,
        protologue_url TEXT,
        locus_classicus_text TEXT,
        tdwg_native_codes_json TEXT,
        tdwg_introduced_codes_json TEXT,
        natureserve_element_id TEXT,
        natureserve_global_rank TEXT,
        nwpl_indicator_code TEXT,
        nwpl_wetland_score REAL,
        backbone_source TEXT,
        updated_at TEXT
    );
    CREATE INDEX IF NOT EXISTS idx_tbp_target_genus ON taxon_backbone_profiles(target_genus);

    -- 10. NatureServe Subnational US State / Canadian Province Ranges
    CREATE TABLE IF NOT EXISTS natureserve_subnational_ranges (
        species_standardized TEXT NOT NULL,
        country_code TEXT NOT NULL,
        subnation_code TEXT NOT NULL,
        s_rank TEXT,
        is_native INTEGER DEFAULT 1,
        is_exotic INTEGER DEFAULT 0,
        is_historical_or_extirpated INTEGER DEFAULT 0,
        PRIMARY KEY (species_standardized, country_code, subnation_code)
    );
    CREATE INDEX IF NOT EXISTS idx_ns_ranges_taxon_sub ON natureserve_subnational_ranges(species_standardized, subnation_code);

"""

__all__ = ("_build_tier_a_schema_ddl",)
