#!/usr/bin/env python3
"""
Unit and integration tests for Step 06 ingestion, VCF/dosage polyploid scaling,
2x-anchored HDplot + ParaLoss-1272 weights, nQuire cytotype guards, concept/type
holdouts, coordinate quarantine, and master-view replicate cardinality.
"""

from __future__ import annotations

import csv
import importlib
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main as root_main
from scripts.core.db import get_db_connection
from scripts.core.db_schema import PACKERA_SCHEMA_USER_VERSION, PLANTMATICS_SCHEMA_USER_VERSION

INGEST_MOD = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")


def _build_synthetic_vcf_text(
    samples: list[str], locus_specs: list[tuple[str, list[tuple[str, int, int, int]]]],
) -> str:
    """Build a minimal VCF v4.2 payload with `GT:DP:AD` calls per locus."""
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for pos, (lid, cells) in enumerate(locus_specs, start=1):
        fmt_cells = [f"{gt}:{dp}:{r},{a}" for gt, dp, r, a in cells]
        lines.append(f"contig_{pos}\t{100 + pos}\t{lid}\tA\tG\t99\tPASS\t.\tGT:DP:AD\t" + "\t".join(fmt_cells))
    return "\n".join(lines) + "\n"


class TestStep06IngestionHdplotPloidy(unittest.TestCase):
    """Step 06 ingestion, schema contract, HDplot/ParaLoss-1272, nQuire, and hygiene tests."""

    def test_gitignore_tracks_gitkeep_and_ignores_genomic_payloads(self) -> None:
        text = (PROJECT_ROOT / ".gitignore").read_text(encoding="utf-8")
        lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
        for required in (
            "data/genomics/*", "data/genomics/**", "!data/genomics/.gitkeep",
            "*.vcf", "*.vcf.gz", "*.fastq", "*.gz", "*.phylip", "*.nex", "*.bam", "*.tree", "*raxml*",
        ):
            self.assertIn(required, lines, f"Missing required .gitignore pattern: {required}")

        git_bin = shutil.which("git")
        if git_bin and (PROJECT_ROOT / ".git").exists():
            res_keep = subprocess.run([git_bin, "check-ignore", "data/genomics/.gitkeep"], cwd=PROJECT_ROOT, capture_output=True, text=True, check=False)
            self.assertNotEqual(res_keep.returncode, 0, "data/genomics/.gitkeep must remain unignored")
            for probe in (
                "data/genomics/synthetic_test.vcf", "data/genomics/snpsaurus_nextrad/synthetic_test.vcf",
                "data/genomics/nextrad_raw_reads/synthetic_reads_R1.fastq.gz", "data/genomics/snpsaurus_nextrad/RAxML_bestTree.synthetic",
            ):
                res_ign = subprocess.run([git_bin, "check-ignore", probe], cwd=PROJECT_ROOT, capture_output=True, text=True, check=False)
                self.assertEqual(res_ign.returncode, 0, f"Genomic payload not ignored by git: {probe}")

    def test_master_view_genomic_replicate_cardinality_invariance(self) -> None:
        self.assertEqual(PLANTMATICS_SCHEMA_USER_VERSION, 20261127)
        self.assertEqual(PACKERA_SCHEMA_USER_VERSION, 20261127)
        with tempfile.TemporaryDirectory() as td:
            ws = Path(td)
            data_dir, tables_dir, runs_dir = ws / "data", ws / "data" / "tables", ws / "runs" / "run_alpha"
            tables_dir.mkdir(parents=True, exist_ok=True)
            runs_dir.mkdir(parents=True, exist_ok=True)
            conn = get_db_connection(db_path=data_dir / "plantmatics.sqlite", init_schema=True, workspace_root=ws)
            try:
                gs_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genomic_samples);").fetchall()}
                for req_col in (
                    "exsiccatae_key", "genbank_accession", "coordinateUncertainty",
                    "is_centroid_or_inst_coord", "coord_quarantine_reason",
                    "typified_basionym", "concept_synonymy_status", "is_open_or_hybrid_concept",
                ):
                    self.assertIn(req_col, gs_cols, f"Missing genomic_samples column: {req_col}")
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, species_standardized, determiner_tier) VALUES (?, ?, ?);",
                    [("NCU001", "Packera dubia", "Tier_1_Gold"), ("NCU002", "Packera schweinitziana", "Tier_1_Gold")],
                )
                conn.executemany(
                    "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, rosette_centroid_x) VALUES (?, ?, ?);",
                    [("NCU001", 0, 100.0), ("NCU001", 1, 200.0), ("NCU002", 0, 150.0)],
                )
                conn.executemany(
                    "INSERT INTO genomic_samples (genetic_sample_id, catalogNumber, linkage_tier, species_standardized, geno_quality_pass, mapped_reads, inferred_cytotype, hdplot_retained_loci) VALUES (?, ?, ?, ?, ?, ?, ?, ?);",
                    [
                        ("PACTOMNCU1M", "NCU001", "G1_Direct_Voucher", "Packera dubia", 1, 500000, "4x", 1200),
                        ("PACTOMNCU1P", "NCU001", "G2_Population_Voucher_Proxy", "Packera dubia", 1, 420000, "4x", 1180),
                    ],
                )
                conn.execute("INSERT INTO pipeline_runs (run_id, run_dir, intraspecific_mode) VALUES ('run_alpha', ?, 0);", (str(runs_dir),))
                conn.executemany(
                    "INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier, authority_weight) VALUES (?, ?, ?, ?, ?, ?);",
                    [
                        ("run_alpha", "NCU001", 0, "Packera dubia", "Tier_1_Gold", 2.0),
                        ("run_alpha", "NCU001", 1, "Packera dubia", "Tier_1_Gold", 2.0),
                        ("run_alpha", "NCU002", 0, "Packera schweinitziana", "Tier_1_Gold", 2.0),
                    ],
                )
                conn.executemany(
                    "INSERT INTO genetic_lineage_assignments (run_id, genetic_sample_id, catalogNumber, geno_predicted_taxon, geno_posterior_prob, hybrid_index_h, interspecific_ancestry_q12, hybrid_generation_class, coalescent_gdi) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
                    [
                        ("run_alpha", "PACTOMNCU1M", "NCU001", "Packera dubia", 0.99, 0.04, 0.05, "Pure_Parent_1", 0.88),
                        ("run_alpha", "PACTOMNCU1P", "NCU001", "Packera dubia", 0.96, 0.06, 0.07, "Pure_Parent_1", 0.88),
                    ],
                )
                conn.execute(
                    "INSERT INTO curation_triage (run_id, catalogNumber, plant_individual_id, recommended_determination, synthesis_triage_priority, taxonomic_treatment_subtype) VALUES ('run_alpha', 'NCU001', 0, 'Packera dubia', 'RESOLVED', 'Confirmed_Species');"
                )
                conn.commit()
                self.assertEqual(conn.execute("SELECT COUNT(*) FROM v_integrative_triage_master;").fetchone()[0], 3)
            finally:
                conn.close()

    def test_empty_genomics_transient_and_flat_dir_guards(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir = tmp / "genomics", tmp / "parsed_out"
            run_dir, tables_dir, figs_dir = tmp / "runs" / "latest", tmp / "data_tables", tmp / "figures"
            for d in (geno_dir, run_dir, tables_dir, figs_dir):
                d.mkdir(parents=True, exist_ok=True)
            for fn, payload in (
                (".gitkeep", b""), ("reads.fastq.gz.crdownload", b"partial"), ("bundle.tar.gz.part", b"partial"),
                ("scratch.tmp", b"tmp"), (".hidden.vcf", b"##fileformat=VCFv4.2\n"),
                ("raw_sample_R1.fastq.gz", b"raw"), ("aligned.bam", b"bam"),
            ):
                (geno_dir / fn).write_bytes(payload)
            (geno_dir / "nextrad_raw_reads").mkdir(parents=True, exist_ok=True)
            (geno_dir / "nextrad_raw_reads" / "sample.F.fq.g").write_bytes(b"raw")

            self.assertFalse(INGEST_MOD.discover_genomic_inputs(geno_dir).has_genetic_payload)
            self.assertEqual(INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir)["status"], "skipped_empty")
            self.assertEqual(len(list(geno_dir.iterdir())), 8, "Empty no-op must never write into data/genomics/")
            with self.assertRaises(ValueError):
                INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=geno_dir / "parsed_subdir")

            cli_args = type("Ns", (), {"intraspecific_mode": False, "focal_taxon": None, "db_path": None, "demo": False, "genomics_dir": None})()
            with mock.patch.dict(os.environ, {"PLANTMATICS_GENOMICS_DIR": str(geno_dir)}), mock.patch("subprocess.run") as m_run:
                self.assertIsNone(root_main.run_genomics_stream(cli_args, run_dir, tables_dir, figs_dir, "Rscript"))
                m_run.assert_not_called()

            (geno_dir / "snpsaurus_nextrad").mkdir(parents=True, exist_ok=True)
            (geno_dir / "snpsaurus_nextrad" / "packera_test.vcf").write_text(
                _build_synthetic_vcf_text(["PACTOMSCG1M", "PACSCHNCB1M"], [("loc1", [("0/0", 20, 20, 0), ("1/1", 20, 0, 20)])]),
                encoding="utf-8",
            )
            with mock.patch.dict(os.environ, {"PLANTMATICS_GENOMICS_DIR": str(geno_dir)}), mock.patch("subprocess.run") as m_run:
                self.assertEqual(root_main.run_genomics_stream(cli_args, run_dir, tables_dir, figs_dir, "Rscript"), run_dir / "outputs" / "genomics")
                self.assertEqual(m_run.call_count, 2)
                cmd0 = " ".join(str(x) for x in m_run.call_args_list[0].args[0])
                cmd1 = " ".join(str(x) for x in m_run.call_args_list[1].args[0])
                self.assertIn("06_ingest_genetic_data.py", cmd0)
                self.assertIn("06_genomic_phylogenomics.R", cmd1)
                self.assertIn("fig06a_genomic_pcoa_tree_concordance.png", cmd1)
                self.assertIn("fig06b_reticulation_gdi_landscape.png", cmd1)
                self.assertTrue(any("fig06c_cytonuclear_paralog_cytotype_audit.png" in str(a) for a in m_run.call_args_list[1].args[0]))

    def test_hdplot_eb_dosage_cytotype_and_run5563_linker(self) -> None:
        self.assertEqual(INGEST_MOD._extract_cell_depths(["0/1", "40", "12,28"], {"GT": 0, "DP": 1, "AD": 2}, 0.5), (12.0, 28.0))
        self.assertEqual(INGEST_MOD._extract_cell_depths(["0/1", "40", "18"], {"GT": 0, "DP": 1, "AD": 2}, 0.5), (22.0, 18.0))

        n_loci = 40
        ref_2x, alt_2x, dos_2x = np.full(n_loci, 20.0, dtype=np.float32), np.full(n_loci, 20.0, dtype=np.float32), np.full(n_loci, 0.5, dtype=np.float32)
        self.assertEqual(INGEST_MOD._infer_sample_cytotype(ref_2x, alt_2x, dos_2x, None)[0], "2x")
        ref_4x = np.array([30.0 if i % 2 == 0 else 10.0 for i in range(n_loci)], dtype=np.float32)
        alt_4x = np.array([10.0 if i % 2 == 0 else 30.0 for i in range(n_loci)], dtype=np.float32)
        cyto_4x, dlog_4x = INGEST_MOD._infer_sample_cytotype(ref_4x, alt_4x, dos_2x, None)
        self.assertEqual(cyto_4x, "4x")
        self.assertIsNotNone(dlog_4x)
        self.assertGreater(float(dlog_4x), 2.0)

        hex_pairs = [(52.0, 8.0)] * 30 + [(8.0, 52.0)] * 30 + [(40.0, 20.0)] * 10 + [(20.0, 40.0)] * 10
        cyto_6x, dlog_6x = INGEST_MOD._infer_sample_cytotype(
            np.array([p[0] for p in hex_pairs], dtype=np.float32),
            np.array([p[1] for p in hex_pairs], dtype=np.float32),
            np.full(len(hex_pairs), 0.5, dtype=np.float32), None,
        )
        self.assertEqual(cyto_6x, "6x")
        self.assertIsNotNone(dlog_6x)
        self.assertGreater(float(dlog_6x), 2.0)
        self.assertEqual(INGEST_MOD._infer_sample_cytotype(ref_2x[:5], alt_2x[:5], dos_2x[:5], 6)[0], "6x")

        eb = INGEST_MOD._compute_eb_polyrad_dosage(
            np.array([[0.5, 0.0], [0.5, 1.0], [0.5, 0.5]], dtype=np.float32),
            np.array([[30.0, 30.0], [10.0, 0.0], [20.0, 20.0]], dtype=np.float32),
            np.array([[10.0, 0.0], [30.0, 30.0], [20.0, 20.0]], dtype=np.float32),
            ["4x", "4x", "2x"],
        )
        self.assertEqual(eb.shape, (3, 2))
        self.assertTrue(np.all(np.isfinite(eb)) and np.all((eb >= 0.0) & (eb <= 1.0)))
        self.assertLess(float(eb[0, 0]), 0.40)
        self.assertGreater(float(eb[1, 0]), 0.60)

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir, db_path = tmp / "genomics", tmp / "parsed", tmp / "packera_test.sqlite"
            geno_dir.mkdir(parents=True, exist_ok=True)
            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, recordNumber, recordedBy, scientificName, species_standardized, decimalLatitude, decimalLongitude, stateProvince, determiner_tier) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
                    [
                        ("NCU00445952", "PACSCH01-NCB1", "Fuller", "Packera schweinitziana", "Packera schweinitziana", 35.80, -82.90, "NC", "Tier_1_Gold"),
                        ("NCU00112233", "PACTOM-SCG1-01", "Fuller", "Packera tomentosa", "Packera dubia", 35.10, -79.20, "NC", "Tier_1_Gold"),
                    ],
                )
                conn.commit()
            finally:
                conn.close()

            samples = ["NCU00445952", "PACTOMSCG1M", "PACTOMSCG1P", "PACMILVAL1M", "PACPAUVAL1M", "PACXYZ99M"]
            locus_specs = [
                ("L_PARALOG_HET", [("0/1", 40, 20, 20) for _ in samples]),
                ("L_PARALOG_DBIAS", [("0/1", 40, 36, 4)] * 3 + [("0/0", 40, 40, 0)] * 3),
                ("L_MONOMORPHIC", [("0/0", 30, 30, 0) for _ in samples]),
                ("L_LOWMAF", [("0/0", 30, 30, 0)] * 5 + [("0/1", 30, 15, 15)]),
            ] + [
                (f"L_GOOD_{k + 1}", [
                    ("1/1", 30, 0, 30), ("0/0", 30, 30, 0), ("0/0", 30, 30, 0),
                    ("0/1", 40, 30, 10) if k % 2 == 0 else ("0/1", 40, 10, 30),
                    ("1/1", 30, 0, 30), ("0/1", 30, 15, 15),
                ])
                for k in range(20)
            ]
            (geno_dir / "packera.vcf").write_text(_build_synthetic_vcf_text(samples, locus_specs), encoding="utf-8")
            (geno_dir / "packera_key.txt").write_text("SampleID\tLibrary\tBarcode\n" + "\n".join(f"{s}\tlib1\tAA{i:02d}" for i, s in enumerate(samples)) + "\n", encoding="utf-8")
            (geno_dir / "packera_stats.txt").write_text(
                "Sample\tHQ reads\tHQ reads mapped to reference\tpercent reads mapped to reference\tpercent repetitive reads\t# ./. genotypes\t# 0/0 genotypes\t# 0/1 genotypes\t# 1/1 genotypes\n"
                + "\n".join(f"{s}\t400000\t360000\t90.0\t4.0\t0\t8\t4\t7" for s in samples) + "\n",
                encoding="utf-8",
            )
            (geno_dir / "user_genetic_metadata.csv").write_text(
                "genetic_sample_id,catalogNumber,decimalLatitude,decimalLongitude,taxon_raw\n"
                "PACMILVAL1M,,36.15,-81.65,Packera millefolium\nPACPAUVAL1M,,36.15,-81.65,Packera paupercula\n",
                encoding="utf-8",
            )

            res = INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir, db_path=db_path)
            self.assertEqual((res["status"], res["n_samples_total"], res["n_loci_raw"], res["n_loci_retained"]), ("ingested", 6, 24, 21))
            with open(out_dir / "hdplot_locus_diagnostics.tsv", encoding="utf-8", newline="") as fh:
                hd_rows = {row["locus_id"]: row for row in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual((hd_rows["L_PARALOG_HET"]["hdplot_keep"], hd_rows["L_PARALOG_HET"]["filter_reason"]), ("0", "HDplot_H>0.60"))
            self.assertEqual(hd_rows["L_PARALOG_DBIAS"]["hdplot_keep"], "0")
            self.assertIn("HDplot_|D|>", hd_rows["L_PARALOG_DBIAS"]["filter_reason"])
            self.assertEqual((hd_rows["L_MONOMORPHIC"]["hdplot_keep"], hd_rows["L_MONOMORPHIC"]["filter_reason"]), ("0", "monomorphic_or_sparse"))
            # Singleton private allele (L_LOWMAF) is 100% retained with zero MAF pruning
            self.assertEqual((hd_rows["L_LOWMAF"]["hdplot_keep"], hd_rows["L_LOWMAF"]["filter_reason"]), ("1", "PASS"))
            self.assertEqual((hd_rows["L_GOOD_1"]["hdplot_keep"], hd_rows["L_GOOD_1"]["filter_reason"]), ("1", "PASS"))

            with open(out_dir / "parsed_genomic_samples.csv", encoding="utf-8", newline="") as fh:
                samp_rows = {row["genetic_sample_id"]: row for row in csv.DictReader(fh)}
            self.assertEqual((samp_rows["NCU00445952"]["linkage_tier"], samp_rows["NCU00445952"]["catalogNumber"]), ("G1_Direct_Voucher", "NCU00445952"))
            self.assertEqual((samp_rows["PACTOMSCG1M"]["linkage_tier"], samp_rows["PACTOMSCG1M"]["catalogNumber"]), ("G2_Population_Voucher_Proxy", "NCU00112233"))
            self.assertEqual((samp_rows["PACMILVAL1M"]["linkage_tier"], samp_rows["PACMILVAL1M"]["inferred_cytotype"]), ("G3_Georeferenced_Unvouchered", "4x"))
            self.assertEqual(samp_rows["PACXYZ99M"]["linkage_tier"], "G4_Ungeoreferenced_Tip")
            with sqlite3.connect(str(db_path)) as conn_c:
                self.assertEqual(conn_c.execute("SELECT COUNT(*) FROM genomic_samples;").fetchone()[0], 6)

    def test_w4_post_hdplot_diploid_x23_guard_stays_2x(self) -> None:
        from scripts.data_prep import step06_ploidy_hdplot_utils as hd

        n_sc, n_para = 30, 50
        sub_pass = np.tile(np.array([[0.0], [0.5], [1.0]], dtype=np.float32), (2, n_sc + n_para))
        sub_ref = sub_pass.copy()
        ref_d, alt_d = np.full((6, n_sc + n_para), 20.0, dtype=np.float32), np.full((6, n_sc + n_para), 20.0, dtype=np.float32)
        sub_ref[:, n_sc:], ref_d[:, n_sc:], alt_d[:, n_sc:] = 0.25, 30.0, 10.0
        keep_mask, _ = hd._evaluate_hdplot_loci(sub_pass, sub_ref, ref_d, alt_d, [f"L{i}" for i in range(n_sc + n_para)])
        self.assertTrue(bool(np.all(keep_mask[:n_sc])))
        self.assertFalse(bool(np.any(keep_mask[n_sc:])))

        ref_v = np.concatenate([np.full(n_sc, 20.0), np.full(n_para, 30.0)])
        alt_v = np.concatenate([np.full(n_sc, 20.0), np.full(n_para, 10.0)])
        dos_v = np.full(n_sc + n_para, 0.5, dtype=np.float32)
        self.assertEqual(hd._infer_sample_cytotype(ref_v[keep_mask], alt_v[keep_mask], dos_v[keep_mask], None, species_modal_ploidy=2, is_multi_cytotype=False)[0], "2x")
        self.assertEqual(
            hd._infer_sample_cytotype(
                np.array([30.0, 10.0] * 20), np.array([10.0, 30.0] * 20), np.full(40, 0.5, dtype=np.float32), None,
                species_modal_ploidy=2, is_multi_cytotype=False, min_het_loci=60, strong_delta_bic=25.0, min_delta_logl_per_locus=0.15,
            ),
            ("2x", 0.0),
        )

    def test_polyploid_vcf_ds_and_rowwise_mixed_ploidy_dosage_scaling(self) -> None:
        """Verify 4x/6x VCF FORMAT/DS ploidy scaling and row-wise mixed-ploidy dosage table scaling."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            vcf_p = tmp / "poly.vcf"
            vcf_p.write_text(
                "##fileformat=VCFv4.2\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                '##FORMAT=<ID=DS,Number=1,Type=Float,Description="DS">\n'
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2X\tS4X\tS6X\n"
                "c1\t10\tL1\tA\tG\t99\tPASS\t.\tGT:DS:AD\t0/1:1.0:15,15\t0/0/1/1:2.0:20,20\t0/0/0/1/1/1:3.0:21,21\n"
                "c1\t20\tL2\tA\tG\t99\tPASS\t.\tGT:DS:AD\t1/1:2.0:0,30\t0/0/0/1:1.0:30,10\t0/0/0/0/0/1:1.0:30,6\n",
                encoding="utf-8",
            )
            sids, lids, dmat, _, _ = INGEST_MOD._parse_vcf_with_depths(vcf_p)
            self.assertEqual((sids, lids), (["S2X", "S4X", "S6X"], ["L1", "L2"]))
            np.testing.assert_allclose(dmat, [[0.50, 1.00], [0.50, 0.25], [0.50, 1.0 / 6.0]], atol=1e-3)

            dos_p = tmp / "mixed_ploidy_dosage.tsv"
            dos_p.write_text("genetic_sample_id\tL1\tL2\tL3\tL4\tL5\nS2X\t0\t1\t2\t1\t0\nS4X\t0\t1\t2\t3\t4\nS6X\t0\t1\t3\t5\t6\n", encoding="utf-8")
            ds_ids, _, ds_mat = INGEST_MOD._parse_dosage_table(dos_p)
            self.assertEqual(ds_ids, ["S2X", "S4X", "S6X"])
            np.testing.assert_allclose(ds_mat[0], [0.0, 0.5, 1.0, 0.5, 0.0], atol=1e-4)
            np.testing.assert_allclose(ds_mat[1], [0.0, 0.25, 0.5, 0.75, 1.0], atol=1e-4)
            np.testing.assert_allclose(ds_mat[2], [0.0, 1.0 / 6.0, 0.5, 5.0 / 6.0, 1.0], atol=1e-4)

    def test_2x_anchored_hdplot_paraloss_rescue_weights_concepts_and_rooftop_quarantine(self) -> None:
        """Verify 2x-anchored HDplot, 0.65 ParaLoss rescue, locus_weights.tsv, concept holdout, and vouchered rooftop quarantine."""
        self.assertFalse(INGEST_MOD.normalize_taxon_name("Packera \u00d7memmingeri")[2])
        self.assertFalse(INGEST_MOD.normalize_taxon_name("Packera x memmingeri")[2])
        for qry in ("Packera cf. aurea", "Packera aff. dubia", "Packera aurea x Packera obovata"):
            self.assertTrue(INGEST_MOD.normalize_taxon_name(qry)[2], f"Expected query=True for {qry!r}")

        q_res = INGEST_MOD.sanitize_step06_sample_coordinates(35.9108, -79.0505, state_prov="NC", county="Orange", locality="NCU Herbarium", is_vouchered=True)
        self.assertEqual((q_res[0], q_res[1], bool(q_res[2])), (None, None, True))

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir, db_path = tmp / "genomics", tmp / "parsed", tmp / "audit.sqlite"
            (geno_dir / "public_catalogs").mkdir(parents=True, exist_ok=True)
            (geno_dir / "public_catalogs" / "public_packera_paraloss1272_loci.csv").write_text(
                "locus_id,paralogy_class,locus_paralogy_weight\nL_SC_1,single_copy,1.5\nL_MC_1,multi_copy,0.5\nL_PL_RESCUE,single_copy,1.5\n",
                encoding="utf-8",
            )
            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, species_standardized, decimalLatitude, decimalLongitude, stateProvince, county, type_status, is_centroid_or_inst_coord, coord_quarantine_reason, concept_synonymy_status, is_open_or_hybrid_concept) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                    [
                        ("V_ROOF", "Packera dubia", 35.9108, -79.0505, "NC", "Orange", "none", 0, "", "accepted", 0),
                        ("V_HET", "Packera dubia", 35.40, -79.30, "NC", "Moore", "ISOTYPE", 0, "", "heterotypic_synonym", 0),
                        ("V_OPEN", "Packera cf. dubia", 35.42, -79.32, "NC", "Moore", "none", 0, "", "open_nomenclature", 1),
                    ],
                )
                conn.commit()
            finally:
                conn.close()

            sids = ["D1", "D2", "D3", "D4", "P1", "P2", "P3", "P4", "H_SYN", "H_OPEN"]
            loci = [
                ("L_SC_1", [("0/0", 30, 30, 0), ("0/0", 30, 30, 0), ("0/1", 40, 20, 20), ("1/1", 30, 0, 30)] + [("0/1", 40, 30, 10)] * 4 + [("0/1", 40, 36, 4)] * 2),
                ("L_UN_1", [("0/0", 30, 30, 0), ("1/1", 30, 0, 30), ("0/1", 40, 20, 20), ("0/0", 30, 30, 0)] + [("0/0", 30, 30, 0), ("1/1", 30, 0, 30), ("0/0", 30, 30, 0), ("1/1", 30, 0, 30)] + [("0/1", 40, 36, 4)] * 2),
                ("L_MC_1", [("0/0", 30, 30, 0), ("1/1", 30, 0, 30), ("0/0", 30, 30, 0), ("1/1", 30, 0, 30)] * 2 + [("0/1", 40, 36, 4)] * 2),
                ("L_DROP_1", [("0/1", 40, 20, 20)] * 10),
            ]
            (geno_dir / "cohort.vcf").write_text(_build_synthetic_vcf_text(sids, loci), encoding="utf-8")
            (geno_dir / "user_genetic_metadata.csv").write_text(
                "genetic_sample_id,catalogNumber,species_standardized,ploidy_level\n"
                "D1,V_ROOF,Packera dubia,2\nD2,,Packera dubia,2\nD3,,Packera dubia,2\nD4,,Packera dubia,2\n"
                "P1,,Packera dubia,4\nP2,,Packera dubia,4\nP3,,Packera dubia,4\nP4,,Packera dubia,4\n"
                "H_SYN,V_HET,Packera dubia,2\nH_OPEN,V_OPEN,Packera cf. dubia,2\n",
                encoding="utf-8",
            )

            res = INGEST_MOD.run_genomic_ingestion(geno_dir, out_dir, db_path=db_path, min_hq_reads=0)
            self.assertEqual(res["status"], "ingested")
            with open(out_dir / "allele_dosage_matrix.tsv", encoding="utf-8") as fh:
                dos_cols = next(csv.reader(fh, delimiter="\t"))[1:]
            with open(out_dir / "locus_weights.tsv", encoding="utf-8") as fh:
                w_rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual([r["locus_id"] for r in w_rows], dos_cols)
            w_by_id = {r["locus_id"]: r for r in w_rows}
            self.assertEqual(
                (w_by_id["L_SC_1"]["paraloss_status"], w_by_id["L_UN_1"]["paraloss_status"], w_by_id["L_MC_1"]["paraloss_status"]),
                ("single_copy", "unannotated", "multi_copy"),
            )
            self.assertAlmostEqual(float(w_by_id["L_SC_1"]["locus_paralogy_weight"]), 1.5, places=3)
            self.assertAlmostEqual(float(w_by_id["L_UN_1"]["locus_paralogy_weight"]), 1.0, places=3)
            self.assertAlmostEqual(float(w_by_id["L_MC_1"]["locus_paralogy_weight"]), 0.5, places=3)
            self.assertNotIn("L_DROP_1", w_by_id)
            with open(out_dir / "hdplot_locus_diagnostics.tsv", encoding="utf-8") as fh:
                hd_map = {r["locus_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(hd_map["L_SC_1"]["hdplot_keep"], "1")
            # Exactly 4 diploid anchors used: excludes 4x P1..P4 AND held-out 2x H_SYN (heterotypic type) & H_OPEN (cf.)
            self.assertEqual(int(hd_map["L_UN_1"]["n_called_ref"]), 4)
            with open(out_dir / "parsed_genomic_samples.csv", encoding="utf-8") as fh:
                sm = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            self.assertEqual(sm["D1"]["decimalLatitude"], "")
            self.assertEqual(sm["D1"]["is_centroid_or_inst_coord"], "1")

            geno_poly, out_poly = tmp / "geno_poly", tmp / "out_poly"
            (geno_poly / "public_catalogs").mkdir(parents=True, exist_ok=True)
            (geno_poly / "public_catalogs" / "public_packera_paraloss1272_loci.csv").write_text(
                "locus_id,paralogy_class,locus_paralogy_weight\nL_PL_RESCUE,single_copy,1.5\n", encoding="utf-8",
            )
            p8 = [f"TP{i}" for i in range(8)]
            p_loci = [("L_PL_RESCUE", [("0/1", 40, 20, 20)] * 5 + [("0/0", 30, 30, 0)] * 3), ("L_UNANN_HI", [("0/1", 40, 20, 20)] * 5 + [("0/0", 30, 30, 0)] * 3)]
            (geno_poly / "poly8.vcf").write_text(_build_synthetic_vcf_text(p8, p_loci), encoding="utf-8")
            (geno_poly / "user_genetic_metadata.csv").write_text(
                "genetic_sample_id,species_standardized,ploidy_level\n" + "".join(f"{s},Packera dubia,4\n" for s in p8), encoding="utf-8",
            )
            INGEST_MOD.run_genomic_ingestion(geno_poly, out_poly, min_hq_reads=0)
            with open(out_poly / "hdplot_locus_diagnostics.tsv", encoding="utf-8") as fh:
                hd_p = {r["locus_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual((hd_p["L_PL_RESCUE"]["hdplot_keep"], hd_p["L_UNANN_HI"]["hdplot_keep"]), ("1", "0"))


if __name__ == "__main__":
    unittest.main()
