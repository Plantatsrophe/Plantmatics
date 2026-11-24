#!/usr/bin/env python3
"""Unit tests for Step 01 portal metadata-lag conflict detection, exsiccatae DB reconciliation, and outgroup balancing."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.core.taxonomy_authority import stratify_duplicates


class TestPortalMetadataLagPythonGates(unittest.TestCase):
    """Test Gates 1 & 2 Python metadata-lag conflict detection, exsiccatae DB reconciliation, and outgroup balancing."""

    def test_verbatim_conflict_detection_and_duplicate_reconciliation(self):
        from scripts.core.db import get_db_connection
        from scripts.core.harvester import _detect_verbatim_identification_conflict
        from scripts.core.taxonomy_authority import reconcile_cross_taxon_exsiccatae_with_db

        # 1. _detect_verbatim_identification_conflict detects stale vs updated identification
        self.assertEqual(
            _detect_verbatim_identification_conflict(
                {"verbatimIdentification": "Packera schweinitziana", "identifiedBy": "T.M. Barkley"},
                "Packera dubia",
                "Tier_1_Gold",
            ),
            (True, "Packera schweinitziana"),
        )
        self.assertEqual(
            _detect_verbatim_identification_conflict(
                {"identifications": [{"scientificName": "Packera anonyma", "identifiedBy": "T.M. Barkley", "isCurrent": True}]},
                "Packera dubia",
                "Tier_1_Gold",
            ),
            (True, "Packera anonyma"),
        )
        self.assertEqual(
            _detect_verbatim_identification_conflict(
                {"previousIdentifications": "Packera anonyma", "verbatimIdentification": "Packera anonyma"},
                "Packera dubia",
                "Tier_1_Gold",
            ),
            (False, ""),
        )

        # 2. stratify_duplicates prioritizes in-range Tier_1_Gold over extralimital Tier_1_Gold
        df = pd.DataFrame([
            {
                "catalogNumber": "DUP_NY_OOR",
                "species_standardized": "Packera dubia",
                "scientificName": "Packera dubia",
                "recordedBy": "House",
                "recordNumber": "26514",
                "eventDate": "1939-06-20",
                "identifiedBy": "T.M. Barkley",
                "determiner_tier": "Tier_1_Gold",
                "usda_county_status": "out_of_range",
                "exsiccatae_key": "house_26514_1939-06-20",
                "specimen_completeness_score": 0.9,
            },
            {
                "catalogNumber": "DUP_PA_INR",
                "species_standardized": "Packera schweinitziana",
                "scientificName": "Packera schweinitziana",
                "recordedBy": "House",
                "recordNumber": "26514",
                "eventDate": "1939-06-20",
                "identifiedBy": "T.M. Barkley",
                "determiner_tier": "Tier_1_Gold",
                "usda_county_status": "in_range",
                "exsiccatae_key": "house_26514_1939-06-20",
                "specimen_completeness_score": 0.85,
            },
        ])
        strat = stratify_duplicates(df)
        self.assertEqual(int(strat["exsiccatae_taxonomic_conflict"].iloc[0]), 1)
        self.assertEqual(int(strat["exsiccatae_taxonomic_conflict"].iloc[1]), 0)
        self.assertTrue(bool(strat["is_primary_duplicate"].iloc[1]))
        self.assertEqual(strat["exsiccatae_gold_taxon"].iloc[0], "Packera schweinitziana")
        self.assertEqual(strat["exsiccatae_gold_taxon"].iloc[1], "Packera schweinitziana")

        # Unanimous in-range disagreement vs stale duplicate
        df_unam = pd.DataFrame([
            {
                "catalogNumber": "U1",
                "species_standardized": "Packera dubia",
                "recordedBy": "Stone",
                "recordNumber": "10",
                "eventDate": "1930-06-01",
                "determiner_tier": "Tier_3_Bronze",
                "usda_county_status": "out_of_range",
                "exsiccatae_key": "stone_10_1930-06-01",
                "specimen_completeness_score": 0.8,
            },
            {
                "catalogNumber": "U2",
                "species_standardized": "Packera schweinitziana",
                "recordedBy": "Stone",
                "recordNumber": "10",
                "eventDate": "1930-06-01",
                "determiner_tier": "Tier_1_Gold",
                "usda_county_status": "in_range",
                "exsiccatae_key": "stone_10_1930-06-01",
                "specimen_completeness_score": 0.9,
            },
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            db_p = Path(tmpdir) / "plantmatics.sqlite"
            conn = get_db_connection(db_p, init_schema=True)
            with conn:
                conn.execute(
                    """
                    INSERT INTO vouchers (
                        catalogNumber, species_standardized, scientificName, determiner_tier,
                        usda_county_status, exsiccatae_key, exsiccatae_taxonomic_conflict, exsiccatae_gold_taxon
                    ) VALUES ('U2', 'Packera schweinitziana', 'Packera schweinitziana', 'Tier_1_Gold',
                              'in_range', 'stone_10_1930-06-01', 0, 'Packera schweinitziana');
                    """
                )
            conn.close()
            conn_r = get_db_connection(db_p, init_schema=False, read_only=True)
            try:
                rec = reconcile_cross_taxon_exsiccatae_with_db(df_unam.iloc[[0]].copy(), conn_r)
            finally:
                conn_r.close()
            self.assertEqual(int(rec["exsiccatae_taxonomic_conflict"].iloc[0]), 1)
            self.assertEqual(rec["exsiccatae_gold_taxon"].iloc[0], "Packera schweinitziana")

    def test_outgroup_anchor_prioritizes_focal_states_with_tier_cascade(self):
        from scripts.core.taxonomy_authority import select_cascading_outgroup_anchors

        focal_df = pd.DataFrame([
            {"catalogNumber": "F1", "species_standardized": "Packera dubia", "stateProvince": "New York"},
            {"catalogNumber": "F2", "species_standardized": "Packera dubia", "stateProvince": "North Carolina"},
        ])
        cands = []
        for k in range(4):
            cands.append({
                "catalogNumber": f"OG_NY_G_{k}",
                "species_standardized": "Packera schweinitziana",
                "determiner_tier": "Tier_1_Gold",
                "stateProvince": "New York",
                "image_url": f"https://example.org/ny_{k}.jpg",
                "specimen_completeness_score": 0.9,
            })
        for k in range(6):
            cands.append({
                "catalogNumber": f"OG_NH_G_{k}",
                "species_standardized": "Packera schweinitziana",
                "determiner_tier": "Tier_1_Gold",
                "stateProvince": "New Hampshire",
                "image_url": f"https://example.org/nh_{k}.jpg",
                "specimen_completeness_score": 0.9,
            })
        out_dfs, counts, _ = select_cascading_outgroup_anchors(
            pd.DataFrame(cands),
            outgroup_taxa=["Packera schweinitziana"],
            max_anchors=5,
            focal_df=focal_df,
        )
        self.assertEqual(counts["Packera schweinitziana"], 5)
        out_df = out_dfs[0]
        self.assertEqual(len(out_df), 5)
        # Round-robin across states guarantees both NY (focal overlap) and NH Gold anchors appear
        self.assertGreaterEqual(int((out_df["stateProvince"] == "New York").sum()), 2)


if __name__ == "__main__":
    unittest.main()
