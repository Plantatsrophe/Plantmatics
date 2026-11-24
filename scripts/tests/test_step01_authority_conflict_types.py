#!/usr/bin/env python3
"""Unit tests for Step 01 authority leakage prevention, exsiccatae conflict auditing, and type priority."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from scripts.core.harvester import export_curated_table
from scripts.core.taxonomy_authority import (
    is_excluded_western_region,
    parse_determiner_tier,
    stratify_duplicates,
)


class TestPhase1AuthorityHygieneAndExsiccataeConflict(unittest.TestCase):
    """Tests for collector-as-determiner leakage prevention, corpus authority duplicate screening, and exsiccatae conflict auditing."""

    def test_collector_as_determiner_leakage_prevented_in_parse_determiner_tier(self):
        """Specialist in recordedBy must NOT promote a generalist/student identifiedBy sheet to Tier_1_Gold unless identifiedBy is blank or references collector."""
        # 1. Specialist collected, student determined at local herbarium -> Tier_3_Bronze (NOT Tier_1_Gold)
        tier_student_local, _, det_raw_1 = parse_determiner_tier(
            type_status_raw="",
            identified_by_raw="J. Undergrad Student",
            recorded_by_raw="T.M. Barkley",
            history_raw="",
            institution_code_raw="UNCC",
            locality_raw="Granite outcrop near Concord",
            habitat_raw="Shallow sandy soil",
        )
        self.assertEqual(tier_student_local, "Tier_3_Bronze")
        self.assertEqual(det_raw_1, "J. Undergrad Student")

        # 2. Specialist collected, generalist determined at major research herbarium (NCU) -> Tier_2_Silver (NOT Tier_1_Gold)
        tier_generalist_ncu, _, _ = parse_determiner_tier(
            type_status_raw="",
            identified_by_raw="J. Generalist",
            recorded_by_raw="T.M. Barkley",
            history_raw="",
            institution_code_raw="NCU",
            locality_raw="Wake Co., granite flatrock 3 km E of Rolesville",
            habitat_raw="Granitic outcrop seepage margin",
            dynamic_specialists={"mahoney", "kowal"},
        )
        self.assertEqual(tier_generalist_ncu, "Tier_2_Silver")

        # 3. Dynamic corpus specialist in recordedBy, student in identifiedBy -> NOT Tier_1_Gold
        tier_dyn_leak, _, _ = parse_determiner_tier(
            type_status_raw="",
            identified_by_raw="Field Intern",
            recorded_by_raw="A.M. Mahoney",
            history_raw="",
            institution_code_raw="UNCC",
            locality_raw="Roadside meadow",
            habitat_raw="",
            dynamic_specialists={"mahoney"},
        )
        self.assertNotEqual(tier_dyn_leak, "Tier_1_Gold")

        # 4. Specialist collected, identifiedBy is blank or explicitly references collector -> Tier_1_Gold
        for id_by in ["", "collector", "id. by collector", "same"]:
            with self.subTest(identified_by=id_by):
                tier_fb, _, _ = parse_determiner_tier(
                    type_status_raw="",
                    identified_by_raw=id_by,
                    recorded_by_raw="T.M. Barkley",
                    history_raw="",
                    institution_code_raw="UNCC",
                    locality_raw="Prairie slope",
                    habitat_raw="",
                )
                self.assertEqual(tier_fb, "Tier_1_Gold")

    def test_discover_corpus_authorities_disqualifies_high_duplicate_disagreement(self):
        """Candidate corpus specialists with >25% cross-duplicate taxonomic disagreement across shared exsiccatae_key groups are disqualified."""
        from scripts.core.taxonomy_authority import discover_corpus_authorities

        rows = []
        insts = ["NCU", "NY", "MO", "GH", "US", "BRIT"]
        # Consistent determiner: 6 determinations across 6 herbaria, 3 duplicate groups all concordant (0% disagreement)
        for i in range(6):
            grp_num = 100 + (i // 2)
            rows.append({
                "catalogNumber": f"CONS-{i:02d}",
                "institutionCode": insts[i],
                "scientificName": "Packera obovata",
                "identifiedBy": "E.R. Consistent",
                "recordedBy": "R. Collector",
                "recordNumber": str(grp_num),
                "eventDate": "2020-05-10",
            })

        # Noisy determiner: 6 determinations of Packera aurea across 6 herbaria,
        # participating in 3 duplicate groups where 2 of 3 (66.7% > 25%) have duplicate sheets determined as Packera obovata
        for i in range(6):
            grp_num = 200 + (i % 3)
            rows.append({
                "catalogNumber": f"NOISY-{i:02d}",
                "institutionCode": insts[i],
                "scientificName": "Packera aurea",
                "identifiedBy": "B.L. Noisy",
                "recordedBy": "R. Collector",
                "recordNumber": str(grp_num),
                "eventDate": "2021-04-18",
            })
        # Add conflicting duplicate partners for groups 200 and 201 (2 of 3 duplicate groups conflicted = 66.7%)
        rows.append({
            "catalogNumber": "DUP-PARTNER-200",
            "institutionCode": "TENN",
            "scientificName": "Packera obovata",
            "identifiedBy": "T.M. Barkley",
            "recordedBy": "R. Collector",
            "recordNumber": "200",
            "eventDate": "2021-04-18",
        })
        rows.append({
            "catalogNumber": "DUP-PARTNER-201",
            "institutionCode": "DUKE",
            "scientificName": "Packera obovata",
            "identifiedBy": "T.M. Barkley",
            "recordedBy": "R. Collector",
            "recordNumber": "201",
            "eventDate": "2021-04-18",
        })

        dyn_specs, _ = discover_corpus_authorities(pd.DataFrame(rows))
        self.assertIn("consistent", dyn_specs)
        self.assertNotIn("noisy", dyn_specs)

    def test_stratify_duplicates_conflict_elevation_sqlite_persistence_and_26_col_csv(self):
        """Conflicting Tier_3_Bronze home-herbarium primary duplicate elevates Tier_1_Gold duplicate determination, logs to SQLite, and preserves 26-column CSV."""
        import sqlite3
        from scripts.core.harvester import EXPORT_COLUMNS

        df_in = pd.DataFrame([
            {
                "catalogNumber": "NCU-DUP-1",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "species_raw": "Packera aurea",
                "recordedBy": "J. B. Fuller",
                "recordNumber": "1042",
                "identifiedBy": "Field Student",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "",
                "county": "Wake",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.85,
                "decimalLongitude": -78.55,
                "coordinateUncertainty": 50.0,
                "year": 2022,
                "month": 5,
                "day": 1,
                "eventDate": "2022-05-01",
                "regional_group": "Piedmont",
                "image_path": "",
            },
            {
                "catalogNumber": "MO-DUP-2",
                "institutionCode": "MO",
                "scientificName": "Packera obovata",
                "species_raw": "Packera obovata",
                "recordedBy": "J. B. Fuller",
                "recordNumber": "1042",
                "identifiedBy": "T.M. Barkley",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "",
                "county": "Wake",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.85,
                "decimalLongitude": -78.55,
                "coordinateUncertainty": 50.0,
                "year": 2022,
                "month": 5,
                "day": 1,
                "eventDate": "2022-05-01",
                "regional_group": "Piedmont",
                "image_path": "",
            },
        ])

        strat_df = stratify_duplicates(df_in, home_herbarium="NCU")
        self.assertTrue((strat_df["exsiccatae_gold_taxon"] == "Packera obovata").all())
        mo_gold = strat_df[strat_df["catalogNumber"] == "MO-DUP-2"].iloc[0]
        self.assertFalse(bool(mo_gold["exsiccatae_taxonomic_conflict"]))
        self.assertTrue(bool(mo_gold["is_primary_duplicate"]))

        brz_dup = strat_df[strat_df["catalogNumber"] == "NCU-DUP-1"].iloc[0]
        self.assertTrue(bool(brz_dup["exsiccatae_taxonomic_conflict"]))
        self.assertFalse(bool(brz_dup["is_primary_duplicate"]))
        # Conflicted Bronze sheet preserves original taxon, tier, and determiner without silent relabeling
        self.assertEqual(brz_dup["scientificName"], "Packera aurea")
        self.assertEqual(brz_dup["determiner_tier"], "Tier_3_Bronze")
        self.assertEqual(brz_dup["species_raw"], "Packera aurea")
        self.assertEqual(brz_dup["identifiedBy"], "Field Student")

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            out_csv = tmp_p / "curated_vouchers.csv"
            db_path = tmp_p / "test_packera.sqlite"
            export_curated_table(strat_df, out_csv, home_herbarium="NCU", db_path=db_path)

            csv_df = pd.read_csv(out_csv)
            from scripts.core.harvester import EXPORT_COLUMNS
            self.assertEqual(len(EXPORT_COLUMNS), 23)
            self.assertEqual(list(csv_df.columns[: len(EXPORT_COLUMNS)]), EXPORT_COLUMNS)
            self.assertIn("exsiccatae_taxonomic_conflict", csv_df.columns)
            self.assertIn("exsiccatae_gold_taxon", csv_df.columns)

            conn = sqlite3.connect(str(db_path))
            try:
                row_db = conn.execute(
                    "SELECT scientificName, determiner_tier, exsiccatae_taxonomic_conflict, exsiccatae_gold_taxon "
                    "FROM vouchers WHERE catalogNumber = 'NCU-DUP-1'"
                ).fetchone()
            finally:
                conn.close()
            self.assertIsNotNone(row_db)
            self.assertEqual(row_db[0], "Packera aurea")
            self.assertEqual(row_db[1], "Tier_3_Bronze")
            self.assertEqual(int(row_db[2]), 1)
            self.assertEqual(row_db[3], "Packera obovata")

    def test_classify_nomenclatural_type_rank_and_cascading_outgroup_type_priority(self):
        """Verify classify_nomenclatural_type_rank(), classify_nomenclatural_type_subrank(), type-annotator leakage guard, and cascading outgroup type priority."""
        from scripts.core.taxonomy_authority import (
            classify_nomenclatural_type_rank,
            classify_nomenclatural_type_subrank,
            discover_corpus_authorities,
            select_cascading_outgroup_anchors,
        )

        self.assertEqual(classify_nomenclatural_type_rank("Holotype of Senecio dubius")[0], 2)
        self.assertEqual(classify_nomenclatural_type_rank("Isolectotype")[0], 2)
        self.assertEqual(classify_nomenclatural_type_rank("Paratype")[0], 1)
        self.assertEqual(classify_nomenclatural_type_rank("Syntype")[0], 1)
        self.assertEqual(classify_nomenclatural_type_rank("Not a type"), (0, "None"))
        self.assertEqual(classify_nomenclatural_type_rank(None), (0, "None"))

        # Verify ICN Art. 9.6 sub-ranking: Holotype > Isotype > Syntype > Paratype > Non-type
        self.assertGreater(
            classify_nomenclatural_type_subrank("Holotype"),
            classify_nomenclatural_type_subrank("Isotype"),
        )
        self.assertGreater(
            classify_nomenclatural_type_subrank("Syntype"),
            classify_nomenclatural_type_subrank("Paratype"),
        )
        self.assertEqual(classify_nomenclatural_type_subrank("Not a type"), 0.0)

        # Verify Layer 3 type-annotator leakage guard: a single-sheet annotator on 1 type specimen
        # at 1 herbarium is NOT promoted to dynamic_specialists unless multi-institution criteria are met
        leak_rows = [
            {
                "catalogNumber": "TYPE-SINGLE-1",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "identifiedBy": "Q.Z. Singleannotator",
                "recordedBy": "A. Collector",
                "recordNumber": "10",
                "type_status": "Holotype",
                "eventDate": "1950-05-01",
            },
            {
                "catalogNumber": "TYPE-MULTI-1",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "identifiedBy": "M.K. Multiannotator",
                "recordedBy": "B. Collector",
                "recordNumber": "20",
                "type_status": "Holotype",
                "eventDate": "1960-05-01",
            },
            {
                "catalogNumber": "NONTYPE-MULTI-2",
                "institutionCode": "NY",
                "scientificName": "Packera aurea",
                "identifiedBy": "M.K. Multiannotator",
                "recordedBy": "C. Collector",
                "recordNumber": "21",
                "type_status": "None",
                "eventDate": "1962-05-01",
            },
        ]
        dyn_specs_leak, _ = discover_corpus_authorities(pd.DataFrame(leak_rows))
        self.assertNotIn("singleannotator", dyn_specs_leak)
        self.assertIn("multiannotator", dyn_specs_leak)

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            # Create local images where non-type Gold sheets have LARGER dimensions than the Holotype/Paratype
            holo_img = tmp_p / "holo.jpg"
            para_img = tmp_p / "para.jpg"
            Image.new("RGB", (2000, 2000), color=(120, 120, 120)).save(holo_img)
            Image.new("RGB", (2100, 2100), color=(120, 120, 120)).save(para_img)

            rows = [
                {
                    "catalogNumber": "HOLO-01",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "Holotype",
                    "is_primary_duplicate": True,
                    "image_path": str(holo_img),
                },
                {
                    "catalogNumber": "PARA-01",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "Paratype",
                    "is_primary_duplicate": True,
                    "image_path": str(para_img),
                },
            ]
            for idx in range(5):
                hi_img = tmp_p / f"gold_{idx}.jpg"
                Image.new("RGB", (3500, 3500), color=(120, 120, 120)).save(hi_img)
                rows.append({
                    "catalogNumber": f"GOLD-NONTYPE-{idx:02d}",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "None",
                    "is_primary_duplicate": True,
                    "image_path": str(hi_img),
                })

            df_in = pd.DataFrame(rows)
            outgroup_dfs, counts, breakdown = select_cascading_outgroup_anchors(
                df_in,
                outgroup_taxa=["Packera aurea"],
                max_anchors=3,
                workspace_dir=tmp_p,
            )
            self.assertEqual(counts["Packera aurea"], 3)
            selected_df = outgroup_dfs[0]
            self.assertEqual(len(selected_df), 3)
            sel_cats = list(selected_df["catalogNumber"])
            self.assertEqual(sel_cats[0], "HOLO-01")
            self.assertEqual(sel_cats[1], "PARA-01")
            self.assertEqual(breakdown["Packera aurea"]["Type_Specimens"], 2)
            self.assertEqual(breakdown["Packera aurea"]["Tier_0_Primary_Type"], 1)
            self.assertEqual(breakdown["Packera aurea"]["Tier_0_Secondary_Type"], 1)

    def test_canadian_western_provinces_and_ungeoreferenced_nontype_exclusion(self):
        """Verify Canadian western/northern provinces are excluded and un-georeferenced non-types without state are excluded while types are retained."""
        canadian_western = [
            "BC", "British Columbia", "AB", "Alberta", "SK", "Saskatchewan",
            "YT", "Yukon", "Yukon Territory", "NT", "Northwest Territories", "NU", "Nunavut",
        ]
        for prov in canadian_western:
            with self.subTest(province=prov):
                self.assertTrue(
                    is_excluded_western_region(prov, exclude_western=True),
                    f"Canadian western/northern province/territory '{prov}' should be excluded.",
                )

        # Un-georeferenced non-type lacking stateProvince is excluded when exclude_western=True
        self.assertTrue(
            is_excluded_western_region(None, lat=None, lon=None, exclude_western=True, is_type_specimen=False)
        )
        self.assertTrue(
            is_excluded_western_region("", lat=np.nan, lon=np.nan, exclude_western=True, is_type_specimen=False)
        )
        # Un-georeferenced type specimen without western state/lon is retained
        self.assertFalse(
            is_excluded_western_region(None, lat=None, lon=None, exclude_western=True, is_type_specimen=True)
        )
        # Un-georeferenced non-type WITH valid eastern state/province is retained
        self.assertFalse(
            is_excluded_western_region("North Carolina", lat=None, lon=None, exclude_western=True, is_type_specimen=False)
        )
        self.assertFalse(
            is_excluded_western_region("Ontario", lat=None, lon=None, exclude_western=True, is_type_specimen=False)
        )


if __name__ == "__main__":
    unittest.main()
