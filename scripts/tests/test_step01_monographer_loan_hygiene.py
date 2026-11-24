#!/usr/bin/env python3
"""Unit tests for Step 01 nomenclatural adjustment isolation, institutional batch demotion, and monographer loan protection."""

import unittest

import pandas as pd


class TestPhase1MonographerLoanAndBatchHygiene(unittest.TestCase):
    """Tests for nomenclatural adjustment isolation, institutional batch demotion, and monographer loan protection."""

    def test_nomenclatural_adjustment_isolation_and_institutional_batch_demotion(self):
        """
        Verify that:
          1. `identifiedBy = 'Nomenclatural Adjustment'` is not promoted to `Tier_1_Gold`
             even when `recordedBy` contains a specialist surname.
          2. Undated institutional batch determinations (`>= 25` sheets across `>= 10`
             collectors at one herbarium) are demoted to `Tier_2_Silver` for historical
             collectors while personal collections (`recordedBy == identifiedBy`) remain
             `Tier_1_Gold`.
        """
        from scripts.core.taxonomy_authority import (
            apply_corpus_authority_tiers,
            parse_determiner_tier,
        )

        tier_adj, _, _ = parse_determiner_tier(
            type_status_raw="None",
            identified_by_raw="Nomenclatural Adjustment",
            recorded_by_raw="D. K. Trock",
            history_raw="",
            institution_code_raw="NCU",
            locality_raw="Orange County, Chapel Hill, rich mesic slope",
            habitat_raw="Deciduous forest",
            dynamic_specialists={"trock", "weakley"},
            dynamic_herbaria={"NCU"},
        )
        self.assertNotEqual(tier_adj, "Tier_1_Gold")

        batch_rows = [
            {
                "catalogNumber": "NCU_NOM_ADJ_01",
                "institutionCode": "NCU",
                "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                "species_raw": "Packera dubia",
                "type_status": "None",
                "identifiedBy": "Nomenclatural Adjustment",
                "recordedBy": "A. S. Weakley",
                "dateIdentified": "",
                "locality": "Orange County, Chapel Hill, rich mesic slope",
                "habitat": "Deciduous forest",
            }
        ]
        for idx in range(26):
            batch_rows.append(
                {
                    "catalogNumber": f"NCU_BATCH_HIST_{idx:02d}",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species_raw": "Packera dubia",
                    "type_status": "None",
                    "identifiedBy": "D. K. Trock",
                    "recordedBy": f"J. HistoricalCollector{idx:02d}",
                    "dateIdentified": "",
                    "locality": "Orange County, North Carolina, alluvial woodland",
                    "habitat": "Bottomland forest",
                }
            )
        batch_rows.append(
            {
                "catalogNumber": "NCU_BATCH_PERSONAL_01",
                "institutionCode": "NCU",
                "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                "species_raw": "Packera dubia",
                "type_status": "None",
                "identifiedBy": "D. K. Trock",
                "recordedBy": "D. K. Trock",
                "dateIdentified": "",
                "locality": "Orange County, North Carolina, alluvial woodland",
                "habitat": "Bottomland forest",
            }
        )

        df_scored = apply_corpus_authority_tiers(pd.DataFrame(batch_rows))
        tier_by_cat = dict(zip(df_scored["catalogNumber"], df_scored["determiner_tier"]))

        self.assertNotEqual(tier_by_cat["NCU_NOM_ADJ_01"], "Tier_1_Gold")
        for idx in range(26):
            self.assertEqual(tier_by_cat[f"NCU_BATCH_HIST_{idx:02d}"], "Tier_2_Silver")
        self.assertEqual(tier_by_cat["NCU_BATCH_PERSONAL_01"], "Tier_1_Gold")

    def test_north_american_packera_specialists_and_monographer_loan_protection(self):
        """
        Verify that:
          1. Newly added North American Packera specialists are classified as Tier_1_Gold.
          2. apply_corpus_authority_tiers() preserves Tier_1_Gold for a specialist
             annotating >= 30 sheets across multiple species or herbaria even when
             dateIdentified is missing, while demoting an undated single-herbarium,
             single-species batch transfer (>= 25 sheets, >= 10 collectors, >= 85%
             in 1 species at 1 institution where recordedBy != identifiedBy) to Tier_2_Silver.
        """
        import inspect
        import scripts.core.taxonomy_authority as tax_auth
        from scripts.core.taxonomy_authority import (
            apply_corpus_authority_tiers,
            parse_determiner_tier,
        )

        classify_fn = getattr(tax_auth, "classify_determiner_tier", None)

        def _eval_tier(det_str: str) -> str:
            if callable(classify_fn):
                sig = inspect.signature(classify_fn)
                if "identified_by_raw" in sig.parameters:
                    res = classify_fn(
                        type_status_raw="None",
                        identified_by_raw=det_str,
                        recorded_by_raw="J. Collector",
                        history_raw="",
                        institution_code_raw="NCU",
                        locality_raw="Orange County, Chapel Hill, rich mesic slope",
                        habitat_raw="Deciduous forest",
                    )
                else:
                    res = classify_fn(det_str)
                return res[0] if isinstance(res, tuple) else str(res)
            tier_out, _, _ = parse_determiner_tier(
                type_status_raw="None",
                identified_by_raw=det_str,
                recorded_by_raw="J. Collector",
                history_raw="",
                institution_code_raw="NCU",
                locality_raw="Orange County, Chapel Hill, rich mesic slope",
                habitat_raw="Deciduous forest",
            )
            return tier_out

        specialists = [
            "A.M. Mahoney",
            "R.R. Kowal",
            "D.K. Trock",
            "T.M. Barkley 2002",
            "L.J. Uttal",
            "J.F. Bain",
            "J.L. Golden",
            "C.C. Freeman",
            "G.L. Nesom",
            "E.R. Moore",
            "A.S. Weakley",
            "D.E. Boufford",
        ]
        for spec_name in specialists:
            with self.subTest(specialist=spec_name):
                self.assertEqual(_eval_tier(spec_name), "Tier_1_Gold")
        self.assertEqual(_eval_tier("E.L. Greene"), "Tier_2_Silver")

        rows = []
        loan_species = [
            "Packera paupercula",
            "Packera plattensis",
            "Packera aurea",
            "Packera obovata",
        ]
        loan_herbaria = ["NCU", "WIS", "MO", "NY"]
        for idx in range(32):
            sp = loan_species[idx % len(loan_species)]
            inst = loan_herbaria[idx % len(loan_herbaria)]
            rows.append(
                {
                    "catalogNumber": f"LOAN_MAHONEY_{idx:02d}",
                    "institutionCode": inst,
                    "scientificName": sp,
                    "species_raw": sp,
                    "type_status": "None",
                    "identifiedBy": "A.M. Mahoney",
                    "recordedBy": f"M. LoanCollector{idx:02d}",
                    "dateIdentified": "",
                    "locality": "Calcareous bluff along river valley",
                    "habitat": "Limestone glade",
                }
            )

        for idx in range(28):
            rows.append(
                {
                    "catalogNumber": f"SINGLE_BATCH_{idx:02d}",
                    "institutionCode": "TENN",
                    "scientificName": "Packera anonyma (Wood) W.A. Weber & A. Love",
                    "species_raw": "Packera anonyma",
                    "type_status": "None",
                    "identifiedBy": "T.M. Barkley",
                    "recordedBy": f"H. BatchCollector{idx:02d}",
                    "dateIdentified": "",
                    "locality": "Knox County, Tennessee, roadside clay bank",
                    "habitat": "Disturbed meadow",
                }
            )

        df_out = apply_corpus_authority_tiers(pd.DataFrame(rows))
        tier_map = dict(zip(df_out["catalogNumber"], df_out["determiner_tier"]))

        for idx in range(32):
            self.assertEqual(tier_map[f"LOAN_MAHONEY_{idx:02d}"], "Tier_1_Gold")
        for idx in range(28):
            self.assertEqual(tier_map[f"SINGLE_BATCH_{idx:02d}"], "Tier_2_Silver")


if __name__ == "__main__":
    unittest.main()
