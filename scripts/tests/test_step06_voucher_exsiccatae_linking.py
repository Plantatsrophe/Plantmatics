#!/usr/bin/env python3
"""
Unit and Warm-R integration tests for Step 06 voucher/exsiccatae primary duplicate
linking, Step 01 concept/quarantine propagation into `genomic_samples`, and
`G1_Exsiccatae_Isoduplicate` retention in `step06_dynamic_tier_engine.R`.
"""

from __future__ import annotations

import csv
import importlib
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection
from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

INGEST_MOD = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")


class TestStep06VoucherExsiccataeLinking(unittest.TestCase):
    """Verify exsiccatae primary-duplicate linking, concept/quarantine propagation, and dynamic-tier retention."""

    def test_link_sample_to_voucher_exsiccatae_primary_duplicate_priority(self) -> None:
        """Verify `_link_sample_to_voucher` selects the imaged primary duplicate across herbaria."""
        vouchers = [
            {
                "catalogNumber": "DUKE001", "recordNumber": "4421", "recordedBy": "A. E. Radford",
                "scientificName": "Packera dubia", "species_standardized": "Packera dubia",
                "decimalLatitude": 35.80, "decimalLongitude": -79.05, "stateProvince": "North Carolina",
                "county": "Orange", "eventDate": "1968-05-14", "determiner_tier": "Tier_1_Gold",
                "type_status": "none", "image_path": "", "exsiccatae_key": "radford|4421|1968-05-14",
                "is_primary_duplicate": 0,
            },
            {
                "catalogNumber": "US003", "recordNumber": "4421", "recordedBy": "Radford, A. E.",
                "scientificName": "Packera dubia", "species_standardized": "Packera dubia",
                "decimalLatitude": 35.80, "decimalLongitude": -79.05, "stateProvince": "North Carolina",
                "county": "Orange", "eventDate": "1968-05-14", "determiner_tier": "Tier_1_Gold",
                "type_status": "none", "image_path": "data/raw_vouchers/US003.webp",
                "exsiccatae_key": "radford|4421|1968-05-14", "is_primary_duplicate": 0,
            },
            {
                "catalogNumber": "NCU002", "recordNumber": "4421", "recordedBy": "J. R. Radford",
                "scientificName": "Packera dubia", "species_standardized": "Packera dubia",
                "decimalLatitude": 35.80, "decimalLongitude": -79.05, "stateProvince": "North Carolina",
                "county": "Orange", "eventDate": "1968-05-14", "determiner_tier": "Tier_2_Silver",
                "type_status": "none", "image_path": "data/raw_vouchers/NCU002.webp",
                "exsiccatae_key": "radford|4421|1968-05-14", "is_primary_duplicate": 1,
            },
        ]
        exs_index = INGEST_MOD._build_exsiccatae_index(vouchers)
        parsed = {"population_id": "OHS1", "taxon_code": "DUB", "species_standardized": "Packera dubia"}

        # 1. Pinning unimaged DUKE001 upgrades via sibling resolution to imaged primary NCU002
        res_pin = INGEST_MOD._link_sample_to_voucher("PACDUBOHS1M", parsed, {"catalogNumber": "DUKE001"}, None, vouchers, exsiccatae_index=exs_index)
        self.assertEqual((res_pin["catalogNumber"], res_pin["linkage_tier"], res_pin["exsiccatae_key"]), ("NCU002", "G1_Direct_Voucher", "radford|4421|1968-05-14"))

        # 2. Explicit exsiccatae metadata lookup resolves NCU002 over non-primary imaged US003
        res_meta = INGEST_MOD._link_sample_to_voucher("PACDUBOHS1M", parsed, {"exsiccatae_key": "radford|4421|1968-05-14"}, None, vouchers, exsiccatae_index=None)
        self.assertEqual((res_meta["catalogNumber"], res_meta["linkage_tier"]), ("NCU002", "G1_Direct_Voucher"))

        # 3. Record number fallback expands exsiccatae siblings and selects primary NCU002
        res_rec = INGEST_MOD._link_sample_to_voucher("4421-M", parsed, None, None, vouchers, exsiccatae_index=exs_index)
        self.assertEqual((res_rec["catalogNumber"], res_rec["linkage_tier"]), ("NCU002", "G1_Direct_Voucher"))

    def test_step01_concept_and_quarantine_propagation_to_genomic_samples(self) -> None:
        """Verify Step 01 concept/quarantine/exsiccatae/genbank columns propagate into CSV and SQLite genomic_samples."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir, db_path = tmp / "genomics", tmp / "parsed", tmp / "prop_vouchers.sqlite"
            geno_dir.mkdir(parents=True, exist_ok=True)

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.executemany(
                    """
                    INSERT INTO vouchers (
                        catalogNumber, recordNumber, recordedBy, scientificName, species_standardized,
                        decimalLatitude, decimalLongitude, coordinateUncertainty, stateProvince, county,
                        eventDate, determiner_tier, type_status, image_path, exsiccatae_key, is_primary_duplicate,
                        is_centroid_or_inst_coord, coord_quarantine_reason, typified_basionym,
                        verbatim_taxon_concept, concept_synonymy_status, identification_year, is_open_or_hybrid_concept
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                    """,
                    [
                        (
                            "NCU_ISO_1", "77", "B. Fuller", "Senecio tomentosus", "Packera dubia",
                            35.40, -79.30, 1250.0, "NC", "Moore", "2024-05-01", "Tier_0_Type", "ISOTYPE",
                            "data/raw_vouchers/NCU_ISO_1.webp", "fuller|77|2024-05-01", 1,
                            1, "county_centroid_match", "Senecio tomentosus",
                            "Packera cf. dubia", "heterotypic_synonym", 2024, 1,
                        ),
                        (
                            "NCU_CLN_2", "78", "B. Fuller", "Packera schweinitziana", "Packera schweinitziana",
                            35.82, -82.91, 45.0, "NC", "Haywood", "2024-06-15", "Tier_1_Gold", "none",
                            "data/raw_vouchers/NCU_CLN_2.webp", "fuller|78|2024-06-15", 1,
                            0, "", "Senecio robbinsii",
                            "Packera schweinitziana", "accepted", 2024, 0,
                        ),
                    ],
                )
                conn.commit()
                loaded_v = INGEST_MOD._load_db_vouchers(conn)
                by_cat = {v["catalogNumber"]: v for v in loaded_v}
                self.assertEqual(by_cat["NCU_ISO_1"]["typified_basionym"], "Senecio tomentosus")
                self.assertEqual(by_cat["NCU_ISO_1"]["concept_synonymy_status"], "heterotypic_synonym")
                self.assertEqual(int(by_cat["NCU_ISO_1"]["is_open_or_hybrid_concept"]), 1)
                self.assertEqual(int(by_cat["NCU_ISO_1"]["is_centroid_or_inst_coord"]), 1)
                self.assertEqual(by_cat["NCU_ISO_1"]["coord_quarantine_reason"], "county_centroid_match")
                self.assertAlmostEqual(float(by_cat["NCU_ISO_1"]["coordinateUncertainty"]), 1250.0, places=2)
            finally:
                conn.close()

            (geno_dir / "prop.vcf").write_text(
                "##fileformat=VCFv4.2\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS_ISO_1\tS_CLN_2\n"
                "c1\t10\tL1\tA\tG\t99\tPASS\t.\tGT:AD\t0/0:20,0\t1/1:0,20\n"
                "c1\t20\tL2\tA\tG\t99\tPASS\t.\tGT:AD\t0/1:15,15\t0/0:20,0\n",
                encoding="utf-8",
            )
            (geno_dir / "user_genetic_metadata.csv").write_text(
                "genetic_sample_id,catalogNumber,linkage_tier,genbank_accession,species_standardized\n"
                "S_ISO_1,NCU_ISO_1,G1_Exsiccatae_Isoduplicate,OQ998877,Packera dubia\n"
                "S_CLN_2,NCU_CLN_2,G1_Direct_Voucher,OQ998878,Packera schweinitziana\n",
                encoding="utf-8",
            )

            res = INGEST_MOD.run_genomic_ingestion(geno_dir, out_dir, db_path=db_path, min_hq_reads=0)
            self.assertEqual(res["status"], "ingested")

            with open(out_dir / "parsed_genomic_samples.csv", encoding="utf-8", newline="") as fh:
                csv_rows = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            r_iso = csv_rows["S_ISO_1"]
            self.assertEqual(r_iso["linkage_tier"], "G1_Exsiccatae_Isoduplicate")
            self.assertEqual(r_iso["exsiccatae_key"], "fuller|77|2024-05-01")
            self.assertEqual(r_iso["genbank_accession"], "OQ998877")
            self.assertAlmostEqual(float(r_iso["coordinateUncertainty"]), 1250.0, places=2)
            self.assertEqual(r_iso["is_centroid_or_inst_coord"], "1")
            self.assertEqual(r_iso["coord_quarantine_reason"], "county_centroid_match")
            self.assertEqual(r_iso["typified_basionym"], "Senecio tomentosus")
            self.assertEqual(r_iso["concept_synonymy_status"], "heterotypic_synonym")
            self.assertEqual(r_iso["is_open_or_hybrid_concept"], "1")

            with sqlite3.connect(str(db_path)) as conn_c:
                self.assertEqual(conn_c.execute("PRAGMA user_version;").fetchone()[0], PLANTMATICS_SCHEMA_USER_VERSION)
                db_row = conn_c.execute(
                    """
                    SELECT linkage_tier, exsiccatae_key, genbank_accession, coordinateUncertainty,
                           is_centroid_or_inst_coord, coord_quarantine_reason, typified_basionym,
                           concept_synonymy_status, is_open_or_hybrid_concept
                    FROM genomic_samples WHERE genetic_sample_id = 'S_ISO_1';
                    """
                ).fetchone()
                self.assertIsNotNone(db_row)
                self.assertEqual(
                    (db_row[0], db_row[1], db_row[2], float(db_row[3]), db_row[4], db_row[5], db_row[6], db_row[7], db_row[8]),
                    ("G1_Exsiccatae_Isoduplicate", "fuller|77|2024-05-01", "OQ998877", 1250.0, 1, "county_centroid_match", "Senecio tomentosus", "heterotypic_synonym", 1),
                )

    def test_dynamic_tier_engine_retains_g1_exsiccatae_isoduplicate(self) -> None:
        """Verify `step06_dynamic_tier_engine.R` retains and persists `G1_Exsiccatae_Isoduplicate` in Tier 3."""
        if not rscript_available():
            self.skipTest("Rscript unavailable or PLANTMATICS_FAST_AUDIT=1")

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            raw_dir, parsed_dir, out_dir = tmp / "raw_geno", tmp / "parsed", tmp / "out"
            tables_dir, figs_dir = tmp / "tables", tmp / "figs"
            for d in (raw_dir, parsed_dir, out_dir, tables_dir, figs_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path = tmp / "tier3_iso.sqlite"
            run_id = "run_t3_isodup"

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.execute("INSERT INTO pipeline_runs (run_id, run_dir, intraspecific_mode) VALUES (?, ?, 0);", (run_id, str(tmp)))
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, species_standardized, determiner_tier) VALUES (?, ?, 'Tier_1_Gold');",
                    [("NCU_ISO_T3", "Packera dubia"), ("NCU_SANG_T3", "Packera dubia")],
                )
                conn.executemany(
                    "INSERT INTO genomic_samples (genetic_sample_id, catalogNumber, linkage_tier, species_standardized, geno_quality_pass) VALUES (?, ?, ?, ?, 1);",
                    [
                        ("S_ISO_T3", "NCU_ISO_T3", "G1_Exsiccatae_Isoduplicate", "Packera dubia"),
                        ("S_SANG_T3", "NCU_SANG_T3", "G3_Public_Sanger_NonMutating", "Packera dubia"),
                    ],
                )
                conn.commit()
            finally:
                conn.close()

            (parsed_dir / "parsed_genomic_samples.csv").write_text(
                "genetic_sample_id,catalogNumber,species_standardized,marker_platform,linkage_tier,geno_linkage_tier,geno_quality_pass,geno_predicted_taxon,geno_posterior_prob,public_its_match,role_in_run\n"
                "S_ISO_T3,NCU_ISO_T3,Packera dubia,nextRAD,G1_Exsiccatae_Isoduplicate,G1_Exsiccatae_Isoduplicate,1,Packera dubia,0.96,Packera dubia,focal\n"
                "S_SANG_T3,NCU_SANG_T3,Packera dubia,Sanger_Barcode,G3_Public_Sanger_NonMutating,G3_Public_Sanger_NonMutating,1,,,Packera dubia,focal\n",
                encoding="utf-8",
            )
            tm = {
                "active_genomic_tier": "Tier_3_Public_MultiLocus",
                "focal_taxon": "Packera dubia",
                "intraspecific_mode": False,
                "complex_neighbors": ["Packera aurea"],
                "polarizing_sister_taxon": "Packera aurea",
                "outgroup_taxon": "Erechtites hieraciifolius",
            }
            (parsed_dir / "genomic_tier_manifest.json").write_text(json.dumps(tm), encoding="utf-8")

            r_file = tmp / "check_t3_isodup.R"
            r_file.write_text(
                'source("scripts/core/db_utils.R")\n'
                'source("scripts/analysis/step06_public_backbone.R")\n'
                'source("scripts/analysis/step06_dynamic_tier_engine.R")\n'
                f'tm <- jsonlite::fromJSON({json.dumps(str(parsed_dir / "genomic_tier_manifest.json"))}, simplifyVector = TRUE)\n'
                f'opts <- list("run-id" = {json.dumps(run_id)}, demo = FALSE, vouchers = "", "figure-a" = "", "figure-b" = "", "figure-c" = "")\n'
                '.resolve_first_nonempty <- function(...) { for (v in list(...)) if (!is.null(v) && nzchar(trimws(as.character(v)))) return(trimws(as.character(v))); "" }\n'
                f'res <- run_step06_dynamic_public_tier(tm, opts, {json.dumps(str(raw_dir))}, {json.dumps(str(parsed_dir))}, {json.dumps(str(out_dir))}, {json.dumps(str(tables_dir))}, {json.dumps(str(figs_dir))}, {json.dumps(str(db_path))})\n'
                'r_iso <- res[res$genetic_sample_id == "S_ISO_T3", , drop = FALSE]\n'
                'stopifnot(nrow(r_iso) == 1L, identical(r_iso$linkage_tier[1L], "G1_Exsiccatae_Isoduplicate"))\n'
                'cat("TIER3_ISODUP_OK\\n")\n',
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step06", cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Dynamic Tier 3 isoduplicate check failed:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("TIER3_ISODUP_OK", proc.stdout)

            with sqlite3.connect(str(db_path)) as conn_c:
                persisted = [r[0] for r in conn_c.execute("SELECT genetic_sample_id FROM genetic_lineage_assignments WHERE run_id = ?;", (run_id,)).fetchall()]
                self.assertIn("S_ISO_T3", persisted, "G1_Exsiccatae_Isoduplicate must be retained & persisted to SQLite in Tier 3")
                self.assertNotIn("S_SANG_T3", persisted, "G3_Public_Sanger_NonMutating must not mutate SQLite genetic_lineage_assignments")


if __name__ == "__main__":
    unittest.main()
