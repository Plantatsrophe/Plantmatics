"""
Second half of `v_integrative_triage_master` SQLite view DDL (Step 06-07 joins) and assembly.
"""

from __future__ import annotations

from scripts.core.db_schema_master_view_a import _build_master_view_select_part_a

_MASTER_VIEW_PART_B_SQL: str = """\
        -- Step 07 Curation Triage
        ct.w_morph_kappa AS w_morph_kappa,
        ct.w_edaphic_kappa AS w_edaphic_kappa,
        ct.w_pheno_kappa AS w_pheno_kappa,
        ct.w_geno_kappa AS w_geno_kappa,
        ct.n_available_modalities AS n_available_modalities,
        ct.probabilistic_consensus_conflict AS probabilistic_consensus_conflict,
        ct.probabilistic_consensus_score AS probabilistic_consensus_score,
        COALESCE(srf.hopkins_adjusted_doy, ct.hopkins_adjusted_doy) AS hopkins_adjusted_doy,
        COALESCE(ct.delta_doy_resid, srf.hopkins_adjusted_doy_resid) AS delta_doy_resid,
        ct.taxonomic_status_call AS taxonomic_status_call,
        ct.recommended_determination AS recommended_determination,
        ct.taxonomic_treatment_broad AS taxonomic_treatment_broad,
        ct.taxonomic_treatment_subtype AS taxonomic_treatment_subtype,
        ct.synthesis_triage_priority AS synthesis_triage_priority,
        ct.synthesis_triage_action AS synthesis_triage_action,
        ct.synthesis_rationale AS synthesis_rationale,
        ct.cryptic_cluster_cohesion_pass AS cryptic_cluster_cohesion_pass,
        ct.cryptic_non_clinal_pass AS cryptic_non_clinal_pass,
        ct.resurrection_basionym AS resurrection_basionym,
        ct.hybrid_zone_cohesion_pass AS hybrid_zone_cohesion_pass,
        ct.hybrid_parent_1 AS hybrid_parent_1,
        ct.hybrid_parent_2 AS hybrid_parent_2,
        ct.is_completed_voucher AS is_completed_voucher,
        ct.genomic_propagation_tier AS genomic_propagation_tier,
        ct.tier0_calibrated_taxon AS tier0_calibrated_taxon,
        ct.tier0_calibrated_prob AS tier0_calibrated_prob,
        ct.tier0_abstention_flag AS tier0_abstention_flag,
        ct.deme_mahal_d2 AS deme_mahal_d2,
        ct.cluster_propagated_gdi AS cluster_propagated_gdi,
        ct.cytology_linkage_tier AS cytology_linkage_tier,
        ct.voucher_sporophytic_2n AS voucher_sporophytic_2n,
        ct.voucher_ploidy_level AS voucher_ploidy_level,
        ct.taxon_reported_cytotypes AS taxon_reported_cytotypes,
        ct.is_multi_cytotype_taxon AS is_multi_cytotype_taxon,
        ct.cytotype_reproductive_barrier_pass AS cytotype_reproductive_barrier_pass,
        ct.hydro_edaphic_niche_regime AS hydro_edaphic_niche_regime,
        ct.wetland_fidelity_concordance AS wetland_fidelity_concordance,
        ct.basionym_synonymy_type AS basionym_synonymy_type,
        ct.w_morph_simplex AS w_morph_simplex,
        ct.w_edaphic_simplex AS w_edaphic_simplex,
        ct.w_pheno_simplex AS w_pheno_simplex,
        ct.w_geno_simplex AS w_geno_simplex,
        ct.geno_deme_weight_scale AS geno_deme_weight_scale,
        ct.portal_metadata_lag_flag AS portal_metadata_lag_flag,
        ct.portal_lag_trigger_gate AS portal_lag_trigger_gate,
        ct.is_marginal_range_extension AS is_marginal_range_extension,
        -- Tier A Taxonomic Backbone & NWPL Wetland Indicator
        tbp.basionym_name AS backbone_basionym_name,
        tbp.locus_classicus_text AS locus_classicus_text,
        tbp.protologue_url AS protologue_url,
        tbp.natureserve_global_rank AS natureserve_global_rank,
        tbp.nwpl_indicator_code AS nwpl_indicator_code,
        tbp.nwpl_wetland_score AS nwpl_wetland_score,
        -- Step 06 Genomic Phylogenomics (1-to-1 per run_id, catalogNumber)
        geno_v.genetic_sample_id AS genetic_sample_id,
        geno_v.geno_population_id AS geno_population_id,
        geno_v.geno_replicate_pair_id AS geno_replicate_pair_id,
        geno_v.geno_taxon_raw AS geno_taxon_raw,
        geno_v.geno_species_standardized AS geno_species_standardized,
        geno_v.is_putative_novel_or_hybrid_query AS is_putative_novel_or_hybrid_query,
        geno_v.geno_ploidy_level AS geno_ploidy_level,
        geno_v.geno_marker_platform AS geno_marker_platform,
        geno_v.geno_hq_reads AS geno_hq_reads,
        geno_v.geno_missing_genotype_rate AS geno_missing_genotype_rate,
        geno_v.geno_observed_het AS geno_observed_het,
        geno_v.geno_quality_pass AS geno_quality_pass,
        geno_v.geno_linkage_tier AS geno_linkage_tier,
        geno_v.geno_predicted_taxon AS geno_predicted_taxon,
        geno_v.geno_posterior_prob AS geno_posterior_prob,
        geno_v.geno_margin AS geno_margin,
        geno_v.genetic_cluster AS genetic_cluster,
        geno_v.q_max AS q_max,
        geno_v.q_admixture_json AS q_admixture_json,
        geno_v.is_geno_admixed AS is_geno_admixed,
        geno_v.geno_parent_1 AS geno_parent_1,
        geno_v.geno_parent_2 AS geno_parent_2,
        geno_v.nearest_tree_sister_taxon AS nearest_tree_sister_taxon,
        geno_v.patristic_dist_to_dubia AS patristic_dist_to_dubia,
        geno_v.replicate_pair_concordance AS replicate_pair_concordance,
        geno_v.n_genomic_samples_linked AS n_genomic_samples_linked,
        geno_v.hdplot_retained_loci AS hdplot_retained_loci,
        geno_v.inferred_cytotype AS inferred_cytotype,
        geno_v.nquire_delta_logl AS nquire_delta_logl,
        geno_v.plastid_clade AS plastid_clade,
        geno_v.its_ribotype AS its_ribotype,
        geno_v.cytonuclear_discordant AS cytonuclear_discordant,
        geno_v.hybrid_index_h AS hybrid_index_h,
        geno_v.interspecific_ancestry_q12 AS interspecific_ancestry_q12,
        geno_v.hybrid_generation_class AS hybrid_generation_class,
        geno_v.coancestry_cluster AS coancestry_cluster,
        geno_v.f_branch_max AS f_branch_max,
        geno_v.f_branch_donor_taxon AS f_branch_donor_taxon,
        geno_v.snaq_gamma_minor AS snaq_gamma_minor,
        geno_v.coalescent_gdi AS coalescent_gdi,
        geno_v.bpp_introgression_bf AS bpp_introgression_bf,
        geno_v.mmrr_beta_ibd AS mmrr_beta_ibd,
        geno_v.mmrr_beta_ibe AS mmrr_beta_ibe,
        geno_v.mmrr_p_ibe AS mmrr_p_ibe,
        geno_v.public_plastome_match AS public_plastome_match,
        geno_v.public_its_match AS public_its_match,
        geno_v.waxy_lineage_match AS waxy_lineage_match,
        geno_v.public_backbone_clade AS public_backbone_clade,
        geno_v.public_phytop_ils_index AS public_phytop_ils_index,
        geno_v.public_phytop_ih_index AS public_phytop_ih_index
    FROM run_active_cohort rac
    LEFT JOIN vouchers v
        ON v.catalogNumber = rac.catalogNumber
    LEFT JOIN plant_individuals pi
        ON pi.catalogNumber = rac.catalogNumber
       AND pi.plant_individual_id = rac.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp
        ON vrp.catalogNumber = rac.catalogNumber
       AND vrp.plant_individual_id = rac.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp0
        ON vrp0.catalogNumber = rac.catalogNumber
       AND vrp0.plant_individual_id = 0
    LEFT JOIN otu_morphometrics om
        ON om.run_id = rac.run_id
       AND om.catalogNumber = rac.catalogNumber
       AND om.plant_individual_id = rac.plant_individual_id
    LEFT JOIN spatial_rf_predictions srf
        ON srf.run_id = rac.run_id
       AND srf.catalogNumber = rac.catalogNumber
       AND srf.plant_individual_id = rac.plant_individual_id
    LEFT JOIN environmental_cache ec
        ON ec.lat_lon_key = COALESCE(
            srf.lat_lon_key,
            CASE
                WHEN v.decimalLatitude IS NOT NULL
                 AND v.decimalLongitude IS NOT NULL
                THEN printf('%.5f_%.5f', v.decimalLatitude, v.decimalLongitude)
                ELSE NULL
            END,
            (SELECT ec2.lat_lon_key FROM environmental_cache ec2 WHERE ec2.catalogNumber = rac.catalogNumber LIMIT 1)
        )
    LEFT JOIN curation_triage ct
        ON ct.run_id = rac.run_id
       AND ct.catalogNumber = rac.catalogNumber
       AND ct.plant_individual_id = rac.plant_individual_id
    LEFT JOIN (
        SELECT
            run_id,
            linked_catalog_number AS catalogNumber,
            genetic_sample_id,
            geno_population_id,
            geno_replicate_pair_id,
            geno_taxon_raw,
            geno_species_standardized,
            is_putative_novel_or_hybrid_query,
            geno_ploidy_level,
            geno_marker_platform,
            geno_hq_reads,
            geno_missing_genotype_rate,
            geno_observed_het,
            geno_quality_pass,
            geno_linkage_tier,
            geno_predicted_taxon,
            geno_posterior_prob,
            geno_margin,
            genetic_cluster,
            q_max,
            q_admixture_json,
            is_geno_admixed,
            geno_parent_1,
            geno_parent_2,
            nearest_tree_sister_taxon,
            patristic_dist_to_dubia,
            replicate_pair_concordance,
            n_genomic_samples_linked,
            hdplot_retained_loci,
            inferred_cytotype,
            nquire_delta_logl,
            plastid_clade,
            its_ribotype,
            cytonuclear_discordant,
            hybrid_index_h,
            interspecific_ancestry_q12,
            hybrid_generation_class,
            coancestry_cluster,
            f_branch_max,
            f_branch_donor_taxon,
            snaq_gamma_minor,
            coalescent_gdi,
            bpp_introgression_bf,
            mmrr_beta_ibd,
            mmrr_beta_ibe,
            mmrr_p_ibe,
            public_plastome_match,
            public_its_match,
            waxy_lineage_match,
            public_backbone_clade,
            public_phytop_ils_index,
            public_phytop_ih_index
        FROM (
            SELECT
                gla.run_id AS run_id,
                COALESCE(gla.catalogNumber, gs.catalogNumber) AS linked_catalog_number,
                gla.genetic_sample_id AS genetic_sample_id,
                COALESCE(gla.population_id, gs.population_id) AS geno_population_id,
                gs.replicate_pair_id AS geno_replicate_pair_id,
                gs.taxon_raw AS geno_taxon_raw,
                gs.species_standardized AS geno_species_standardized,
                gs.is_putative_novel_or_hybrid_query AS is_putative_novel_or_hybrid_query,
                gs.ploidy_level AS geno_ploidy_level,
                gs.marker_platform AS geno_marker_platform,
                gs.hq_reads AS geno_hq_reads,
                gs.missing_genotype_rate AS geno_missing_genotype_rate,
                gs.observed_het AS geno_observed_het,
                gs.geno_quality_pass AS geno_quality_pass,
                gs.linkage_tier AS geno_linkage_tier,
                gla.geno_predicted_taxon AS geno_predicted_taxon,
                gla.geno_posterior_prob AS geno_posterior_prob,
                gla.geno_margin AS geno_margin,
                gla.genetic_cluster AS genetic_cluster,
                gla.q_max AS q_max,
                gla.q_admixture_json AS q_admixture_json,
                gla.is_geno_admixed AS is_geno_admixed,
                gla.geno_parent_1 AS geno_parent_1,
                gla.geno_parent_2 AS geno_parent_2,
                gla.nearest_tree_sister_taxon AS nearest_tree_sister_taxon,
                gla.patristic_dist_to_dubia AS patristic_dist_to_dubia,
                gla.replicate_pair_concordance AS replicate_pair_concordance,
                gs.hdplot_retained_loci AS hdplot_retained_loci,
                gs.inferred_cytotype AS inferred_cytotype,
                gs.nquire_delta_logl AS nquire_delta_logl,
                gs.plastid_clade AS plastid_clade,
                gs.its_ribotype AS its_ribotype,
                gs.cytonuclear_discordant AS cytonuclear_discordant,
                gla.hybrid_index_h AS hybrid_index_h,
                gla.interspecific_ancestry_q12 AS interspecific_ancestry_q12,
                gla.hybrid_generation_class AS hybrid_generation_class,
                gla.coancestry_cluster AS coancestry_cluster,
                gla.f_branch_max AS f_branch_max,
                gla.f_branch_donor_taxon AS f_branch_donor_taxon,
                gla.snaq_gamma_minor AS snaq_gamma_minor,
                gla.coalescent_gdi AS coalescent_gdi,
                gla.bpp_introgression_bf AS bpp_introgression_bf,
                gla.mmrr_beta_ibd AS mmrr_beta_ibd,
                gla.mmrr_beta_ibe AS mmrr_beta_ibe,
                gla.mmrr_p_ibe AS mmrr_p_ibe,
                gs.public_plastome_match AS public_plastome_match,
                gs.public_its_match AS public_its_match,
                gs.waxy_lineage_match AS waxy_lineage_match,
                gla.public_backbone_clade AS public_backbone_clade,
                gla.public_phytop_ils_index AS public_phytop_ils_index,
                gla.public_phytop_ih_index AS public_phytop_ih_index,
                COUNT(*) OVER (
                    PARTITION BY gla.run_id, COALESCE(gla.catalogNumber, gs.catalogNumber)
                ) AS n_genomic_samples_linked,
                ROW_NUMBER() OVER (
                    PARTITION BY gla.run_id, COALESCE(gla.catalogNumber, gs.catalogNumber)
                    ORDER BY
                        CASE gs.linkage_tier
                            WHEN 'G1_Direct_Voucher' THEN 1
                            WHEN 'G1_Exsiccatae_Isoduplicate' THEN 2
                            WHEN 'G2_Population_Voucher_Proxy' THEN 3
                            WHEN 'G2_Conditional_Deme_Hull' THEN 4
                            WHEN 'G3_Georeferenced_Unvouchered' THEN 5
                            ELSE 6
                        END ASC,
                        COALESCE(gs.geno_quality_pass, 0) DESC,
                        COALESCE(gs.hq_reads, 0) DESC,
                        COALESCE(gla.geno_posterior_prob, -1.0) DESC,
                        gla.genetic_sample_id ASC
                ) AS _geno_rn
            FROM genetic_lineage_assignments gla
            LEFT JOIN genomic_samples gs
                ON gs.genetic_sample_id = gla.genetic_sample_id
            WHERE COALESCE(gla.catalogNumber, gs.catalogNumber) IS NOT NULL
        ) geno_ranked
        WHERE _geno_rn = 1
    ) geno_v
        ON geno_v.run_id = rac.run_id
       AND geno_v.catalogNumber = rac.catalogNumber
    LEFT JOIN taxon_backbone_profiles tbp
        ON tbp.species_standardized = COALESCE(rac.species_standardized, om.species_standardized, v.species_standardized);
    """

def _build_master_view_ddl() -> str:
    return _build_master_view_select_part_a() + _MASTER_VIEW_PART_B_SQL

MASTER_VIEW_DDL: str = _build_master_view_ddl()

__all__ = (
    "_MASTER_VIEW_PART_B_SQL", "_build_master_view_ddl", "MASTER_VIEW_DDL",
)
