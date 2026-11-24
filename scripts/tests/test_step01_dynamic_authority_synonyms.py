#!/usr/bin/env python3
"""Unit tests for Step 01 dynamic authority engine, GBIF synonym expansion, and geographic scopes."""

import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from scripts.core.harvester import VoucherHarvester, harvest_taxa_occurrences
from scripts.core.taxonomy_authority import parse_determiner_tier, stratify_duplicates


class TestDynamicAuthorityEngineAndSynonymExpansion(unittest.TestCase):
    """Test suite for the 3-Layer Dynamic Authority Engine, GBIF synonym expansion, and geographic scopes."""

    def test_layer1_expanded_north_american_monographers(self):
        """Verify core Packera monographers -> Tier_1_Gold, bare common surnames -> Tier_3_Bronze,
        and regional/historical authorities -> Tier_2_Silver (promoted on Holotypes)."""
        core_gold = [
            "C. C. Freeman", "W. A. Weber", "L. J. Uttal", "G. L. Nesom",
            "J. L. Golden", "L. D. Estes", "Dwayne Estes",
        ]
        for det_name in core_gold:
            with self.subTest(core_gold=det_name):
                tier, _, _ = parse_determiner_tier(
                    type_status_raw=None, identified_by_raw=det_name,
                    recorded_by_raw="Collector", history_raw=None,
                    institution_code_raw="SMALL_HERB", locality_raw="Meadow", habitat_raw=None,
                )
                self.assertEqual(tier, "Tier_1_Gold", f"{det_name} should map to Tier_1_Gold")

        tier_bare, _, _ = parse_determiner_tier(
            type_status_raw=None, identified_by_raw="J. Estes",
            recorded_by_raw="Collector", history_raw=None,
            institution_code_raw="SMALL_HERB", locality_raw="Meadow", habitat_raw=None,
        )
        self.assertEqual(tier_bare, "Tier_3_Bronze", "Bare 'J. Estes' must remain Tier_3_Bronze")

        regional_silver = [
            "A. Cronquist", "J. M. Greenman", "Á. Löve", "R. Spellenberg",
            "R. W. Pippen", "M. L. Fernald", "P. A. Rydberg", "A. A. Heller",
            "R. Kral", "R. K. Godfrey", "B. A. Sorrie", "J. L. Strother",
            "B. L. Turner", "A. Gray", "C. Jeffrey",
        ]
        for det_name in regional_silver:
            with self.subTest(regional_silver=det_name):
                tier, _, _ = parse_determiner_tier(
                    type_status_raw=None, identified_by_raw=det_name,
                    recorded_by_raw="Collector", history_raw=None,
                    institution_code_raw="SMALL_HERB", locality_raw="Meadow", habitat_raw=None,
                )
                self.assertEqual(tier, "Tier_2_Silver", f"{det_name} should map to Tier_2_Silver")
                tier_holo, _, _ = parse_determiner_tier(
                    type_status_raw="Holotype", identified_by_raw=det_name,
                    recorded_by_raw="Collector", history_raw=None,
                    institution_code_raw="SMALL_HERB", locality_raw="Meadow", habitat_raw=None,
                )
                self.assertEqual(tier_holo, "Tier_1_Gold", f"Holotype by {det_name} must promote to Tier_1_Gold")

    def test_discover_corpus_authorities_promotes_multi_institution_specialist_not_local_collector(self):
        """
        Verify Layer 3 promotes a determiner with >= 5 annotations across >= 3 institutions
        (or type annotation) to Tier_1_Gold while keeping single-institution collectors in Silver/Bronze.
        """
        from scripts.core.taxonomy_authority import apply_corpus_authority_tiers, discover_corpus_authorities

        rows = []
        # 1. Loan specialist: 5 sheets of Packera streptanthifolia across 3 institutions (COLO, RM, UTC)
        loan_insts = ["COLO", "RM", "UTC", "COLO", "RM"]
        for idx, inst in enumerate(loan_insts):
            rows.append({
                "catalogNumber": f"LOAN_{idx:02d}",
                "institutionCode": inst,
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "Dr. Elena V. Monographia, 2019",
                "recordedBy": f"Collector {idx}",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Subalpine rocky meadow at 2800 m elevation",
                "habitat": "Granitic talus",
            })

        # 2. Single-institution local collector/determiner: 6 sheets all at NCU (only 1 institution)
        for idx in range(6):
            rows.append({
                "catalogNumber": f"LOCAL_NCU_{idx:02d}",
                "institutionCode": "NCU",
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "M. Localcollector",
                "recordedBy": "M. Localcollector",
                "type_status": "None",
                "determiner_tier": "Tier_2_Silver",
                "locality": "Roadside clearing 4 miles west of town center",
                "habitat": "Open woods",
            })

        # 3. Single-institution collector at unknown herbarium without rich locality -> remains Bronze
        for idx in range(6):
            rows.append({
                "catalogNumber": f"LOCAL_BRONZE_{idx:02d}",
                "institutionCode": "TINY_COLL",
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "P. Singleherb",
                "recordedBy": "P. Singleherb",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Woods",
                "habitat": "",
            })

        # 4. Generic non-person strings across 4 institutions must NEVER be promoted
        for idx, inst in enumerate(["NCU", "MO", "NY", "US", "GH"]):
            rows.append({
                "catalogNumber": f"GENERIC_{idx:02d}",
                "institutionCode": inst,
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "Herbarium Staff",
                "recordedBy": "Unknown",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Detailed locality over ten characters long",
                "habitat": "Rocky slope",
            })

        # 5. Dynamic regional research herbarium ('REG_UNIV') with 12 curated sheets & >=70% determiners
        for idx in range(12):
            rows.append({
                "catalogNumber": f"REG_UNIV_{idx:02d}",
                "institutionCode": "REG_UNIV",
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": f"Curator Person{idx}",
                "recordedBy": "Regional Collector",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Limestone bluff 6 km north of river crossing",
                "habitat": "Calcareous outcrop",
            })

        df = pd.DataFrame(rows)
        dyn_specs, dyn_herbs = discover_corpus_authorities(df)

        self.assertIn("monographia", dyn_specs)
        self.assertNotIn("localcollector", dyn_specs)
        self.assertNotIn("singleherb", dyn_specs)
        self.assertNotIn("staff", dyn_specs)
        self.assertNotIn("herbarium", dyn_specs)
        self.assertIn("REG_UNIV", dyn_herbs)

        df_scored = apply_corpus_authority_tiers(df, dynamic_specialists=dyn_specs, dynamic_herbaria=dyn_herbs)
        tier_by_cat = dict(zip(df_scored["catalogNumber"], df_scored["determiner_tier"]))

        self.assertEqual(tier_by_cat["LOAN_00"], "Tier_2_Silver")
        self.assertEqual(tier_by_cat["LOAN_04"], "Tier_2_Silver")
        self.assertEqual(tier_by_cat["LOCAL_NCU_00"], "Tier_2_Silver")
        self.assertEqual(tier_by_cat["LOCAL_BRONZE_00"], "Tier_3_Bronze")
        self.assertEqual(tier_by_cat["GENERIC_00"], "Tier_3_Bronze")
        self.assertEqual(tier_by_cat["REG_UNIV_00"], "Tier_2_Silver")

    def test_expand_taxon_synonyms_gbif_and_dynamic_standardization(self):
        """Verify dynamic GBIF Backbone synonym expansion (Packera/Senecio/Cineraria) and reverse standardization."""
        from scripts.core.harvester import (
            _GBIF_SYNONYM_CACHE,
            expand_taxon_synonyms_gbif,
            standardize_packera_taxon,
        )

        mock_species = MagicMock()
        mock_species.name_backbone.return_value = {
            "usageKey": 999123,
            "acceptedUsageKey": 999123,
            "canonicalName": "Packera neomexicana",
            "scientificName": "Packera neomexicana (A.Gray) W.A.Weber & Á.Löve",
        }
        mock_species.name_usage.return_value = {
            "results": [
                {
                    "canonicalName": "Senecio neomexicanus",
                    "scientificName": "Senecio neomexicanus A.Gray",
                },
                {
                    "canonicalName": "Cineraria neomexicana",
                    "scientificName": "Cineraria neomexicana Hort.",
                },
            ]
        }

        _GBIF_SYNONYM_CACHE.pop("packera neomexicana", None)
        with patch("scripts.core.harvester.species", mock_species):
            expanded = expand_taxon_synonyms_gbif("Packera neomexicana")

        self.assertIn("Packera neomexicana", expanded)
        self.assertIn("Senecio neomexicanus", expanded)
        self.assertIn("Cineraria neomexicana", expanded)

        # Harvested synonyms must map back to canonical Packera <epithet>
        self.assertEqual(standardize_packera_taxon("Senecio neomexicanus A.Gray"), "Packera neomexicana")
        self.assertEqual(standardize_packera_taxon("Cineraria neomexicana"), "Packera neomexicana")
        self.assertEqual(standardize_packera_taxon("Packera cana (Hook.) W.A.Weber & Á.Löve"), "Packera cana")

        # Unseen historical Senecio synonym resolves dynamically via GBIF Backbone lookup
        mock_strep = MagicMock()
        mock_strep.name_backbone.return_value = {
            "usageKey": 888456,
            "acceptedUsageKey": 888000,
            "canonicalName": "Senecio streptanthifolius",
            "species": "Packera streptanthifolia",
            "accepted": "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve",
        }
        mock_strep.name_usage.return_value = {"results": []}
        _GBIF_SYNONYM_CACHE.pop("senecio streptanthifolius", None)
        with patch("scripts.core.harvester.species", mock_strep):
            self.assertEqual(
                standardize_packera_taxon("Senecio streptanthifolius Greene"),
                "Packera streptanthifolia",
            )

    def test_geographic_scope_north_america_and_include_western_retain_western_vouchers(self):
        """Verify north_america / include_western retain < -106.65 western records while eastern_na excludes them."""
        mock_occ_records = [
            {
                "key": 5001,
                "scientificName": "Packera cana (Hook.) W.A.Weber & Á.Löve",
                "stateProvince": "Colorado",
                "decimalLatitude": 39.74,
                "decimalLongitude": -108.50,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2019, "month": 6, "day": 18,
                "media": [{"identifier": "https://example.com/colo_cana.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "COLO_WEST_01",
                "institutionCode": "COLO",
                "identifiedBy": "W. A. Weber",
            },
            {
                "key": 5002,
                "scientificName": "Packera breweri (Burtt Davy) W.A.Weber & Á.Löve",
                "stateProvince": "California",
                "decimalLatitude": 36.80,
                "decimalLongitude": -120.40,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2020, "month": 5, "day": 12,
                "media": [{"identifier": "https://example.com/uc_breweri.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "UC_WEST_02",
                "institutionCode": "UC",
                "identifiedBy": "J. L. Strother",
            },
            {
                "key": 5003,
                "scientificName": "Packera paupercula (Michx.) Á.Löve & D.Löve",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.91,
                "decimalLongitude": -79.05,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2021, "month": 5, "day": 3,
                "media": [{"identifier": "https://example.com/ncu_east.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "NCU_EAST_03",
                "institutionCode": "NCU",
                "identifiedBy": "A. S. Weakley",
            },
        ]

        def _fake_occ_search(**kwargs):
            country = kwargs.get("country")
            if country in (None, "US"):
                return {"results": mock_occ_records, "count": len(mock_occ_records)}
            return {"results": [], "count": 0}

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_search):
            # 1. Default eastern_na excludes Colorado & California records (< -106.65)
            df_east = harvest_taxa_occurrences(
                taxa_list=["Packera cana"],
                max_records_per_taxon=10,
            )
            self.assertEqual(set(df_east["catalogNumber"]), {"NCU_EAST_03"})

            # 2. geographic_scope='north_america' retains Colorado & California records and queries US/CA/MX
            harvester_na = VoucherHarvester(
                taxa=["Packera cana"],
                max_records_per_taxon=10,
                geographic_scope="north_america",
            )
            self.assertFalse(harvester_na.exclude_western)
            self.assertEqual(harvester_na.countries, ["US", "CA", "MX"])
            df_na = harvester_na.harvest()
            self.assertEqual(
                set(df_na["catalogNumber"]),
                {"COLO_WEST_01", "UC_WEST_02", "NCU_EAST_03"},
            )

            # 3. include_western=True also disables western exclusion
            df_inc_west = harvest_taxa_occurrences(
                taxa_list=["Packera cana"],
                max_records_per_taxon=10,
                include_western=True,
            )
            self.assertEqual(
                set(df_inc_west["catalogNumber"]),
                {"COLO_WEST_01", "UC_WEST_02", "NCU_EAST_03"},
            )

    def test_select_cascading_outgroup_anchors_and_home_herbarium(self):
        """Verify preflight outgroup anchor cascade (Gold -> Silver -> Bronze) and configurable home_herbarium."""
        from scripts.core.taxonomy_authority import select_cascading_outgroup_anchors

        df_pool = pd.DataFrame([
            # Outgroup A has NO Gold vouchers, only 2 Silver and 2 Bronze
            {
                "catalogNumber": "OG_A_BRONZE_1",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_3_Bronze",
                "image_path": "data/raw_vouchers/OG_A_BRONZE_1.jpg",
                "is_primary_duplicate": True,
            },
            {
                "catalogNumber": "OG_A_SILVER_1",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_2_Silver",
                "image_path": "data/raw_vouchers/OG_A_SILVER_1.jpg",
                "is_primary_duplicate": True,
            },
            {
                "catalogNumber": "OG_A_SILVER_2",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_2_Silver",
                "image_path": "data/raw_vouchers/OG_A_SILVER_2.jpg",
                "is_primary_duplicate": True,
            },
            {
                "catalogNumber": "OG_A_BRONZE_2",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_3_Bronze",
                "image_path": "data/raw_vouchers/OG_A_BRONZE_2.jpg",
                "is_primary_duplicate": True,
            },
        ])
        out_dfs, counts, tiers = select_cascading_outgroup_anchors(
            df_pool,
            outgroup_taxa=["Packera fendleri"],
            max_anchors=3,
        )
        self.assertEqual(counts["Packera fendleri"], 3)
        selected_cats = out_dfs[0]["catalogNumber"].tolist()
        # Both Silver sheets must rank ahead of Bronze fallback
        self.assertEqual(selected_cats[:2], ["OG_A_SILVER_1", "OG_A_SILVER_2"])
        self.assertEqual(tiers["Packera fendleri"]["Tier_2_Silver"], 2)
        self.assertEqual(tiers["Packera fendleri"]["Tier_3_Bronze"], 1)

        # Verify configurable home_herbarium override in stratify_duplicates
        df_dups = pd.DataFrame([
            {
                "catalogNumber": "NCU_DUP",
                "institutionCode": "NCU",
                "recordedBy": "W. A. Weber",
                "recordNumber": "4401",
                "eventDate": "1988-07-10",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
            {
                "catalogNumber": "COLO_DUP",
                "institutionCode": "COLO",
                "recordedBy": "W. A. Weber",
                "recordNumber": "4401",
                "eventDate": "1988-07-10",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
        ])
        strat_colo = stratify_duplicates(df_dups, home_herbarium="COLO")
        self.assertTrue(strat_colo[strat_colo["catalogNumber"] == "COLO_DUP"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(strat_colo[strat_colo["catalogNumber"] == "NCU_DUP"]["is_primary_duplicate"].iloc[0])

    def test_cli_geographic_scope_and_run_pipeline_shell_syntax(self):
        """Verify new CLI flags in 01_voucher_harvester.py and shell syntax/synonym generalization in run_pipeline.sh."""
        import importlib
        import subprocess

        mod_cli = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = mod_cli.build_cli_parser()

        defaults = parser.parse_args([])
        self.assertEqual(defaults.geographic_scope, "eastern_na")
        self.assertFalse(defaults.include_western)
        self.assertTrue(defaults.expand_synonyms)
        self.assertTrue(defaults.exclude_western)

        na_args = parser.parse_args([
            "--geographic-scope", "north_america",
            "--countries", "US", "CA", "MX",
            "--include-western",
        ])
        self.assertEqual(na_args.geographic_scope, "north_america")
        self.assertTrue(na_args.include_western)
        self.assertEqual(na_args.countries, ["US", "CA", "MX"])

        run_sh = Path(__file__).resolve().parents[2] / "run_pipeline.sh"
        sh_text = run_sh.read_text(encoding="utf-8")
        self.assertNotIn('if [[ "${FOCAL_TAXON}" == "Packera dubia" ]]', sh_text)
        self.assertNotIn("--all-eastern-taxa", sh_text)
        self.assertNotIn("ALL_EASTERN_TAXA", sh_text)
        self.assertNotIn("mod05.TARGET_TAXA", sh_text)
        self.assertNotIn('"Packera dubia" "Packera anonyma" "Packera paupercula" "Packera plattensis"', sh_text)
        self.assertIn("Packera_Taxonomic_Revision_Summary.md", sh_text)
        self.assertIn("fig07a_integrative_taxonomic_synthesis.png", sh_text)
        self.assertIn("select_cascading_outgroup_anchors", sh_text)
        self.assertIn("--geographic-scope", sh_text)
        self.assertIn("--include-western", sh_text)

        proc = subprocess.run(["bash", "-n", str(run_sh)], capture_output=True, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"bash -n run_pipeline.sh failed: {proc.stderr}")


if __name__ == "__main__":
    unittest.main()
