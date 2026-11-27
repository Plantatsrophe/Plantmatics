#!/usr/bin/env python3
"""
Unit and integration tests for Step 06 off-target organellar/ITS/waxy skimming,
containment-gated cytonuclear reconciliation, IUPAC 4-state substitution model
selection, dual-allele multi-sample SVDquartets, and multi-block dosage gdi.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import random
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis import step06_coalescent_alignment_utils as aln_mod
from scripts.analysis import step06_coalescent_delimitation as coal_mod
from scripts.analysis import step06_coalescent_models as mod_sel
from scripts.data_prep import organellar_cytonuclear_utils as ocu_mod
from scripts.data_prep import skim_offtarget_organellar as skim_mod
from scripts.tests._r_session import rscript_available, run_warm_r_cmd


def _build_r1_fastq_gz_bytes(reads: list[tuple[str, str]]) -> bytes:
    """Compress synthetic FASTQ reads into `.fastq.gz` bytes."""
    buf = [f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in reads]
    return gzip.compress("".join(buf).encode("utf-8"))


class TestStep06OrganellarCoalescent(unittest.TestCase):
    """Tests for organellar skimming, containment gating, IUPAC models, SVDquartets, and gdi."""

    def test_offtarget_organellar_skim_and_coalescent_delimitation(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, raw_dir, out_dir = tmp / "genomics", tmp / "raw_reads", tmp / "coalescent_out"
            geno_dir.mkdir(parents=True, exist_ok=True)
            raw_dir.mkdir(parents=True, exist_ok=True)

            rbcl_seed = skim_mod._SEED_PLASTID_ANCHORS["Packera_rbcL_core"]
            its_seed = skim_mod._SEED_ITS_ANCHORS["Packera_ITS_58S_core"]
            rbcl_hap1 = rbcl_seed[:120] + ("A" if rbcl_seed[120] != "A" else "G") + rbcl_seed[121:]
            its_sch = its_seed
            its_tom = its_seed[:100] + ("C" if its_seed[100] != "C" else "T") + its_seed[101:]

            zip_path = raw_dir / "run_5563_fastqs.zip"
            with zipfile.ZipFile(zip_path, "w") as zf:
                for sid, p_seq, i_seq in (
                    ("PACSCHNCB1M", rbcl_hap1, its_sch),
                    ("PACSCHNCB2M", rbcl_hap1, its_sch),
                    ("PACTOMSCG1M", rbcl_hap1, its_tom),
                ):
                    reads = [(f"{sid}_p_{j}", p_seq) for j in range(5)] + [(f"{sid}_i_{j}", i_seq) for j in range(5)]
                    zf.writestr(f"fastqs/{sid}_R1.fastq.gz", _build_r1_fastq_gz_bytes(reads))

            skim_res = skim_mod.run_offtarget_organellar_skim(fastq_source=zip_path, genomics_dir=geno_dir)
            self.assertEqual(skim_res["status"], "skimmed")
            self.assertEqual([p.name for p in geno_dir.iterdir() if p.is_dir()], ["public_organellar_and_loci"])
            org_dir = geno_dir / "public_organellar_and_loci"
            for fn in ("skimmed_plastid_consensus.fasta", "skimmed_its_consensus.fasta", "cytonuclear_discordance.csv"):
                self.assertTrue((org_dir / fn).is_file())

            with open(org_dir / "cytonuclear_discordance.csv", "r", encoding="utf-8", newline="") as fh:
                cyto_rows = {row["genetic_sample_id"]: row for row in csv.DictReader(fh)}
            self.assertEqual(cyto_rows["PACTOMSCG1M"]["cytonuclear_discordant"], "1")
            self.assertEqual(cyto_rows["PACSCHNCB1M"]["cytonuclear_discordant"], "0")

            rng = np.random.default_rng(7)
            bases = np.array(list("ACGT"))
            anc = rng.choice(bases, size=732)
            seq_b_base = anc.copy()
            for idx in range(0, 732, 14):
                seq_b_base[idx] = "G" if anc[idx] != "G" else "T"

            def _mutate(seq_arr: np.ndarray, n_mut: int, seed_off: int) -> str:
                arr = seq_arr.copy()
                r_local = np.random.default_rng(100 + seed_off)
                if n_mut > 0:
                    for p in r_local.choice(len(arr), size=n_mut, replace=False):
                        arr[p] = "C" if arr[p] != "C" else "A"
                return "".join(arr.tolist())

            phy_records = [
                ("PACTOMSCG1M", _mutate(anc, 1, 1)), ("PACTOMSCG1P", _mutate(anc, 1, 2)),
                ("PACTOMVAL2M", _mutate(anc, 2, 3)), ("PACSCHNCB1M", _mutate(seq_b_base, 1, 4)),
                ("PACSCHNCB1P", _mutate(seq_b_base, 1, 5)), ("PACANOOHS1M", _mutate(anc, 0, 6)),
                ("PACANOOHS2M", _mutate(anc, 1, 7)),
            ]
            (geno_dir / "snpsaurus_nextrad").mkdir(parents=True, exist_ok=True)
            phy_lines = [f" {len(phy_records)} {len(phy_records[0][1])}"] + [f"{nm:<16}{sq}" for nm, sq in phy_records]
            (geno_dir / "snpsaurus_nextrad" / "packera_invariant.phylip").write_text("\n".join(phy_lines) + "\n", encoding="utf-8")

            parsed_csv = tmp / "parsed_samples.csv"
            with open(parsed_csv, "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["genetic_sample_id", "species_standardized", "replicate_pair_id", "geno_quality_pass"])
                for nm, _ in phy_records:
                    sp = "Packera dubia" if nm.startswith("PACTOM") else ("Packera schweinitziana" if nm.startswith("PACSCH") else "Packera anonyma")
                    wr.writerow([nm, sp, nm[:-1], "1"])

            coal_res = coal_mod.run_coalescent_delimitation(
                genomics_dir=geno_dir, out_dir=out_dir, parsed_samples_csv=parsed_csv, block_len=122, n_boot=25,
            )
            self.assertEqual(coal_res["status"], "completed_alignment")
            for req in (
                "svdquartets_support.csv", "bpp_msci_loci.phy", "bpp_msci.ctl",
                "bpp_imap.txt", "coalescent_gdi_delimitation.csv", "coalescent_gdi_samples.csv",
            ):
                self.assertTrue((out_dir / req).is_file(), f"Missing coalescent output: {req}")

            with open(out_dir / "coalescent_gdi_delimitation.csv", "r", encoding="utf-8", newline="") as fh:
                gdi_map = {row["taxon"]: row for row in csv.DictReader(fh)}
            self.assertIn("Packera schweinitziana", gdi_map)
            self.assertGreater(float(gdi_map["Packera schweinitziana"]["coalescent_gdi"]), 0.70)
            self.assertIn("Packera anonyma", gdi_map)
            self.assertLess(float(gdi_map["Packera anonyma"]["coalescent_gdi"]), 0.25)

    def test_skim_offtarget_organellar_baits_public_refs_including_waxy(self) -> None:
        def _seq(seed_int: int, length: int = 140) -> str:
            return "".join(np.random.default_rng(seed_int).choice(list("ACGT"), size=length).tolist())

        cp_ser, cp_aur = _seq(202601), _seq(202602)
        its_ser, its_aur = _seq(202603), _seq(202604)
        waxy_ser, waxy_aur = _seq(202605), _seq(202606)

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, raw_dir = tmp / "genomics", tmp / "raw_reads"
            org_dir = geno_dir / "public_organellar_and_loci"
            org_dir.mkdir(parents=True, exist_ok=True)
            raw_dir.mkdir(parents=True, exist_ok=True)
            (org_dir / "public_packera_18complete_plastomes.fasta").write_text(f">MK_CP_SER Packera serpenticola\n{cp_ser}\n>MK_CP_AUR Packera aurea\n{cp_aur}\n", encoding="utf-8")
            (org_dir / "public_packera_its_ets.fasta").write_text(f">MT_ITS_SER Packera serpenticola\n{its_ser}\n>MT_ITS_AUR Packera aurea\n{its_aur}\n", encoding="utf-8")
            (org_dir / "public_packera_waxy_nuclear.fasta").write_text(f">JX_WAXY_SER Packera serpenticola\n{waxy_ser}\n>JX_WAXY_AUR Packera aurea\n{waxy_aur}\n", encoding="utf-8")

            self.assertEqual(skim_mod._infer_sample_from_fq_name("EhieraciifoliusEREHIENCS1M_S76.F.fq.gz"), "EREHIENCS1M")
            zip_path = raw_dir / "run_5563_public_skim.zip"
            with zipfile.ZipFile(zip_path, "w") as zf:
                for fq_name, p_s, i_s, w_s in (
                    ("fastqs/PACSERNCS1M_R1.fastq.gz", cp_ser, its_ser, waxy_ser),
                    ("fastqs/PaureaPACAURFK1M_S2.F.fq.gz", cp_aur, its_aur, waxy_aur),
                ):
                    sid = fq_name.split("/")[-1].split("_")[0].replace("Paurea", "")
                    reads = [(f"{sid}_cp_{j}", p_s) for j in range(4)] + [(f"{sid}_its_{j}", i_s) for j in range(4)] + [(f"{sid}_waxy_{j}", w_s) for j in range(4)]
                    zf.writestr(fq_name, _build_r1_fastq_gz_bytes(reads))

            res = skim_mod.run_offtarget_organellar_skim(fastq_source=zip_path, genomics_dir=geno_dir)
            self.assertEqual(res["status"], "skimmed")
            self.assertEqual([p.name for p in geno_dir.iterdir() if p.is_dir()], ["public_organellar_and_loci"])
            for out_fa in ("skimmed_plastid_consensus.fasta", "skimmed_its_consensus.fasta", "skimmed_waxy_consensus.fasta"):
                p = org_dir / out_fa
                self.assertTrue(p.is_file(), f"Missing skimmed consensus FASTA: {out_fa}")
                txt = p.read_text(encoding="utf-8")
                self.assertIn(">PACSERNCS1M|", txt)
                self.assertIn(">PACAURFK1M|", txt)

            with open(org_dir / "cytonuclear_discordance.csv", "r", encoding="utf-8", newline="") as fh:
                cyto_rows = {row["genetic_sample_id"]: row for row in csv.DictReader(fh)}
            for sid, sp in (("PACSERNCS1M", "Packera serpenticola"), ("PACAURFK1M", "Packera aurea")):
                self.assertEqual(cyto_rows[sid]["public_plastome_match"], sp)
                self.assertEqual(cyto_rows[sid]["public_its_match"], sp)
                self.assertEqual(cyto_rows[sid]["waxy_lineage_match"], sp)

    def test_organellar_ets_filter_and_complete_plastome_stride1_matching(self) -> None:
        rng = random.Random(42)
        bases = "ACGT"
        its_body = "".join(rng.choice(bases) for _ in range(450))
        its_query = its_body[:260] + "".join(rng.choice(bases) for _ in range(40)) + its_body[300:]
        its_refs = {
            "GU818219.1 Packera porteri external transcribed spacer, partial sequence": its_query,
            "MK123456.1 Packera hartiana internal transcribed spacer 1, 5.8S, internal transcribed spacer 2": its_body,
        }
        self.assertEqual(ocu_mod._match_best_public_lineage(its_query, its_refs), "Packera hartiana")

        plastome = "".join(rng.choice(bases) for _ in range(130000))
        decoy = "".join(rng.choice(bases) for _ in range(3000))
        cp_refs = {
            "NC_099999.1 Packera subnuda complete chloroplast genome": plastome,
            "KX000001.1 Packera porteri psbA-trnH intergenic spacer": decoy,
        }
        kcache: dict = {}
        self.assertEqual(ocu_mod._match_best_public_lineage(plastome[20000:26000], cp_refs, _kcache=kcache), "Packera subnuda")
        self.assertEqual(ocu_mod._match_best_public_lineage(plastome[20001:26001], cp_refs, _kcache=kcache), "Packera subnuda")

    def test_coalescent_delimitation_collapses_technical_replicates_and_broadcasts_gdi(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            geno_dir, phy_dir, out_dir = root / "genomics", root / "genomics" / "snpsaurus_nextrad", root / "out"
            phy_dir.mkdir(parents=True, exist_ok=True)
            out_dir.mkdir(parents=True, exist_ok=True)
            seq_dub = "ACGT" * 24
            seq_ano = "TGCA" * 24
            phy_lines = [
                "6 96",
                f"TEST_TOMSC1M {'NNNN' + seq_dub[4:]}", f"TEST_TOMSC1P {seq_dub}", f"TEST_TOMGA1M {seq_dub}",
                f"TEST_ANOJE1M {seq_ano}", f"TEST_ANOJE1P {seq_ano}", f"TEST_ANONC1M {seq_ano}",
            ]
            (phy_dir / "packera_invariant.phylip").write_text("\n".join(phy_lines) + "\n", encoding="utf-8")
            parsed_csv = out_dir / "parsed_samples.csv"
            pd.DataFrame([
                {"genetic_sample_id": "TEST_TOMSC1M", "species_standardized": "Packera dubia", "replicate_pair_id": "REP_TOM_SC1", "hq_reads": 150000, "missing_genotype_rate": 0.15, "geno_quality_pass": 1, "is_putative_novel_or_hybrid_query": 0},
                {"genetic_sample_id": "TEST_TOMSC1P", "species_standardized": "Packera dubia", "replicate_pair_id": "REP_TOM_SC1", "hq_reads": 300000, "missing_genotype_rate": 0.05, "geno_quality_pass": 1, "is_putative_novel_or_hybrid_query": 0},
                {"genetic_sample_id": "TEST_TOMGA1M", "species_standardized": "Packera dubia", "replicate_pair_id": "REP_TOM_GA1", "hq_reads": 220000, "missing_genotype_rate": 0.08, "geno_quality_pass": 1, "is_putative_novel_or_hybrid_query": 0},
                {"genetic_sample_id": "TEST_ANOJE1M", "species_standardized": "Packera anonyma", "replicate_pair_id": "REP_ANO_JE1", "hq_reads": 180000, "missing_genotype_rate": 0.10, "geno_quality_pass": 1, "is_putative_novel_or_hybrid_query": 0},
                {"genetic_sample_id": "TEST_ANOJE1P", "species_standardized": "Packera anonyma", "replicate_pair_id": "REP_ANO_JE1", "hq_reads": 250000, "missing_genotype_rate": 0.06, "geno_quality_pass": 1, "is_putative_novel_or_hybrid_query": 0},
                {"genetic_sample_id": "TEST_ANONC1M", "species_standardized": "Packera anonyma", "replicate_pair_id": "REP_ANO_NC1", "hq_reads": 240000, "missing_genotype_rate": 0.07, "geno_quality_pass": 1, "is_putative_novel_or_hybrid_query": 0},
            ]).to_csv(parsed_csv, index=False)

            res = coal_mod.run_coalescent_delimitation(genomics_dir=geno_dir, out_dir=out_dir, parsed_samples_csv=parsed_csv, block_len=24, n_boot=10)
            self.assertEqual((res["status"], res["n_samples"]), ("completed_alignment", 4))
            samp_df = pd.read_csv(out_dir / "coalescent_gdi_samples.csv").set_index("genetic_sample_id")
            self.assertEqual(set(samp_df.index), {"TEST_TOMSC1M", "TEST_TOMSC1P", "TEST_TOMGA1M", "TEST_ANOJE1M", "TEST_ANOJE1P", "TEST_ANONC1M"})
            self.assertAlmostEqual(float(samp_df.loc["TEST_TOMSC1M", "coalescent_gdi"]), float(samp_df.loc["TEST_TOMSC1P", "coalescent_gdi"]), places=6)
            self.assertAlmostEqual(float(samp_df.loc["TEST_ANOJE1M", "coalescent_gdi"]), float(samp_df.loc["TEST_ANOJE1P", "coalescent_gdi"]), places=6)

    def test_iupac_heterozygote_partial_likelihood_svdquartets_and_multiblock_dosage_gdi(self) -> None:
        """Verify IUPAC 4-state masks/+ASC/+G4, dual-allele multi-sample SVDquartets, and multi-block dosage gdi."""
        sids4 = ["S1", "S2", "S3", "S4"]
        inv_tail = "ACGT" * 10
        het_col = "RYMKSW" * 3 + "RY"
        hom_col = "ACAGCACAGCACAGCACAGC"
        seqs_het = {"S1": hom_col + inv_tail, "S2": hom_col + inv_tail, "S3": hom_col + inv_tail, "S4": het_col + inv_tail}
        seqs_snp = {"S1": ("ACAGCA" * 10), "S2": ("ACAGCA" * 10), "S3": ("ACAGCA" * 10), "S4": ("RYMKSW" * 10)}
        seqs_inv = {s: hom_col + inv_tail for s in sids4}

        prev_skip = os.environ.get("PLANTMATICS_SKIP_IQTREE2")
        os.environ["PLANTMATICS_SKIP_IQTREE2"] = "1"
        try:
            res_het = mod_sel.select_nucleotide_substitution_model(seqs_het, sids4, alignment_name="het_iupac.phy")
            res_snp = mod_sel.select_nucleotide_substitution_model(seqs_snp, sids4, alignment_name="snp_iupac.phy")
            res_inv = mod_sel.select_nucleotide_substitution_model(seqs_inv, sids4, alignment_name="inv_iupac.phy")
        finally:
            if prev_skip is None:
                os.environ.pop("PLANTMATICS_SKIP_IQTREE2", None)
            else:
                os.environ["PLANTMATICS_SKIP_IQTREE2"] = prev_skip

        self.assertEqual(int(res_het["n_var_sites"]), 20)
        self.assertEqual(int(res_snp["n_var_sites"]), 60)
        self.assertTrue(bool(res_snp["asc_applied"]))
        self.assertIn("+ASC", str(res_snp["best_model"]))
        self.assertTrue(any(bool(r.get("has_gamma")) for r in res_snp["table_rows"]))
        ll_het = next(r["logL"] for r in res_het["table_rows"] if r["model_name"] == "JC69")
        ll_inv = next(r["logL"] for r in res_inv["table_rows"] if r["model_name"] == "JC69")
        self.assertNotAlmostEqual(float(ll_het), float(ll_inv), places=2)

        # Dual-allele (a1, a2) SVDquartets: a1 identical across taxa, tree ((T1,T2),(T3,T4)) exclusively in a2
        rng_sv = np.random.default_rng(202611)
        def _sim_quartet_tree(n_sites: int, order_pairs: Tuple[Tuple[int, int], Tuple[int, int]]) -> np.ndarray:
            def _mut(arr: np.ndarray, p: float) -> np.ndarray:
                flip = rng_sv.random(arr.size) < p
                return np.where(flip, (arr + rng_sv.integers(1, 4, size=arr.size, dtype=np.int8)) % 4, arr).astype(np.int8)
            root = rng_sv.integers(0, 4, size=n_sites, dtype=np.int8)
            clade_a, clade_b = _mut(root, 0.35), _mut(root, 0.35)
            out = np.zeros((4, n_sites), dtype=np.int8)
            for ti in order_pairs[0]:
                out[ti] = _mut(clade_a, 0.15)
            for ti in order_pairs[1]:
                out[ti] = _mut(clade_b, 0.15)
            return out

        n_s = 400
        a1_dual = np.zeros((8, n_s), dtype=np.int8)
        tree_12 = _sim_quartet_tree(n_s, ((0, 1), (2, 3)))
        a2_dual = np.repeat(tree_12, 2, axis=0)
        tx8 = ["T1", "T1", "T2", "T2", "T3", "T3", "T4", "T4"]
        id8 = [f"{t}_{k}" for t in ("T1", "T2", "T3", "T4") for k in ("a", "b")]
        sv_dual = aln_mod.compute_svdquartets_support(a1_dual, a2_dual, tx8, id8, n_boot=30, seed=202611)
        self.assertEqual(len(sv_dual), 1)
        self.assertEqual(sv_dual.loc[0, "best_split"], "T1,T2|T3,T4")
        self.assertLess(float(sv_dual.loc[0, "svd_score_best"]), float(sv_dual.loc[0, "svd_score_second"]))

        # Multi-sample SVDquartets: sample 0 per taxon is constant (400 sites); sample 1 has 380 sites encoding (T1,T3)|(T2,T4)
        a1_ms = np.zeros((8, n_s), dtype=np.int8)
        tree_13 = _sim_quartet_tree(n_s, ((0, 2), (1, 3)))
        for t_i in range(4):
            a1_ms[2 * t_i + 1, :380] = tree_13[t_i, :380]
            a1_ms[2 * t_i + 1, 380:] = -1
        sv_ms = aln_mod.compute_svdquartets_support(a1_ms, a1_ms.copy(), tx8, id8, n_boot=30, seed=202611)
        self.assertEqual(sv_ms.loc[0, "best_split"], "T1,T3|T2,T4")
        self.assertLess(float(sv_ms.loc[0, "svd_score_best"]), float(sv_ms.loc[0, "svd_score_second"]))

        # Multi-block gametic-mismatch dosage gdi: finite CI width > 0 and active bic_coalescent_mixture
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            geno_e, out_d = root / "geno_empty", root / "dos_out"
            geno_e.mkdir(); out_d.mkdir()
            dos_tsv, parsed_csv = root / "dosage.tsv", root / "parsed.csv"
            rng_d = np.random.default_rng(202611)
            n_loc = 120
            mat = np.zeros((6, n_loc), dtype=np.float64)
            # First 60 loci divergent between taxa; next 60 loci mixed/polymorphic with 0.5 hets -> heterogeneous blocks
            mat[:3, :60] = rng_d.choice([0.0, 0.0, 0.5], size=(3, 60))
            mat[3:, :60] = rng_d.choice([1.0, 1.0, 0.5], size=(3, 60))
            mat[:, 60:] = rng_d.choice([0.0, 0.5, 1.0], size=(6, 60))
            s6 = ["PACTOM1M", "PACTOM2M", "PACTOM3M", "PACSCH1M", "PACSCH2M", "PACSCH3M"]
            pd.DataFrame(mat, columns=[f"L{i}" for i in range(n_loc)]).assign(genetic_sample_id=s6)[["genetic_sample_id", *[f"L{i}" for i in range(n_loc)]]].to_csv(dos_tsv, sep="\t", index=False)
            pd.DataFrame([
                {"genetic_sample_id": s, "species_standardized": "Packera dubia" if s.startswith("PACTOM") else "Packera schweinitziana", "replicate_pair_id": s[:-1], "geno_quality_pass": 1}
                for s in s6
            ]).to_csv(parsed_csv, index=False)

            res_d = coal_mod.run_coalescent_delimitation(genomics_dir=geno_e, out_dir=out_d, parsed_samples_csv=parsed_csv, dosage_tsv=dos_tsv, n_boot=60)
            self.assertEqual(res_d["status"], "completed_dosage_fallback")
            gdi_df = pd.read_csv(out_d / "coalescent_gdi_delimitation.csv")
            adeq_df = pd.read_csv(out_d / "coalescent_model_adequacy.csv")
            self.assertGreaterEqual(int(gdi_df["n_loci_blocks"].min()), 4)
            self.assertTrue(any(float(r["gdi_ci_high"]) > float(r["gdi_ci_low"]) for _, r in gdi_df.iterrows()))
            self.assertTrue(all(str(m) == "bic_coalescent_mixture" for m in adeq_df["bf_method"]))
            self.assertTrue(any(math.isfinite(float(v)) and abs(float(v)) > 1e-6 for v in adeq_df["bic_m0"]))

            # Intraspecific deme gdi (.attrs-safe pd.concat + coalescent_gdi_samples.csv deme hydration)
            out_intra, cl_csv = root / "dos_intra_out", root / "clusters.csv"
            out_intra.mkdir()
            mat_i = np.vstack([
                rng_d.choice([0.0, 0.0, 0.0, 0.5], size=(2, n_loc)),
                rng_d.choice([1.0, 1.0, 1.0, 0.5], size=(2, n_loc)),
                rng_d.choice([0.0, 1.0], size=(2, n_loc)),
            ])
            s6_i = ["PACTOM1M", "PACTOM2M", "PACTOM3M", "PACTOM4M", "PACSCH1M", "PACSCH2M"]
            pd.DataFrame(mat_i, columns=[f"L{i}" for i in range(n_loc)]).assign(genetic_sample_id=s6_i)[["genetic_sample_id", *[f"L{i}" for i in range(n_loc)]]].to_csv(dos_tsv, sep="\t", index=False)
            pd.DataFrame([{"genetic_sample_id": s, "species_standardized": "Packera dubia" if s.startswith("PACTOM") else "Packera schweinitziana", "replicate_pair_id": s, "geno_quality_pass": 1} for s in s6_i]).to_csv(parsed_csv, index=False)
            pd.DataFrame([{"genetic_sample_id": s, "genetic_cluster": ("Geno_Cluster_1" if s in {"PACTOM1M", "PACTOM2M"} else ("Geno_Cluster_2" if s in {"PACTOM3M", "PACTOM4M"} else "Polarizing_Sister"))} for s in s6_i]).to_csv(cl_csv, index=False)
            res_in = coal_mod.run_coalescent_delimitation(genomics_dir=geno_e, out_dir=out_intra, parsed_samples_csv=parsed_csv, dosage_tsv=dos_tsv, focal_taxon="Packera dubia", intraspecific_mode=True, cluster_csv=cl_csv, n_boot=40)
            self.assertEqual(res_in["status"], "completed_dosage_fallback")
            gdi_in = pd.read_csv(out_intra / "coalescent_gdi_delimitation.csv")
            adeq_in = pd.read_csv(out_intra / "coalescent_model_adequacy.csv")
            samp_in = pd.read_csv(out_intra / "coalescent_gdi_samples.csv").set_index("genetic_sample_id")
            self.assertTrue({"Geno_Cluster_1", "Geno_Cluster_2"}.issubset(set(gdi_in["taxon"])))
            self.assertTrue({"Geno_Cluster_1", "Geno_Cluster_2"}.issubset(set(adeq_in["taxon"])))
            self.assertEqual((samp_in.loc["PACTOM1M", "lineage_taxon"], samp_in.loc["PACTOM3M", "lineage_taxon"]), ("Geno_Cluster_1", "Geno_Cluster_2"))
            self.assertGreater(float(samp_in.loc["PACTOM3M", "coalescent_gdi"]), 0.50)

    def test_organellar_containment_gate_and_discordance_preservation(self) -> None:
        """Verify >=0.80 containment AND >=80 shared 11-mers gate rejects noisy cont=0.25 hits and preserves discordance."""
        if not rscript_available():
            self.skipTest("Rscript unavailable or PLANTMATICS_FAST_AUDIT=1")

        rng = np.random.default_rng(202611)
        ref_aur = "".join(rng.choice(list("ACGT"), size=160).tolist())
        div_tail = "".join(rng.choice(list("ACGT"), size=120).tolist())
        # Noisy query: shares first 47 bp (~37 11-mers, containment ~0.25) with Packera aurea
        noisy_q = ref_aur[:47] + div_tail[:113]
        self.assertEqual(ocu_mod.match_best_public_lineage_with_metrics(noisy_q, {"MK1 Packera aurea": ref_aur})[0], "")

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir = tmp / "geno", tmp / "out"
            org_dir = geno_dir / "public_organellar_and_loci"
            org_dir.mkdir(parents=True, exist_ok=True)
            out_dir.mkdir(parents=True, exist_ok=True)
            (org_dir / "public_packera_18complete_plastomes.fasta").write_text(f">MK_AUR Packera aurea\n{ref_aur}\n", encoding="utf-8")
            (org_dir / "skimmed_plastid_consensus.fasta").write_text(
                f">S_HQ|clade=cp_1\n{ref_aur}\n>S_NOISY|clade=cp_2\n{noisy_q}\n", encoding="utf-8",
            )
            r_file = tmp / "check_containment_gate.R"
            r_file.write_text(
                f"source({str(PROJECT_ROOT / 'scripts' / 'analysis' / 'step06_public_backbone.R')!r})\n"
                f"source({str(PROJECT_ROOT / 'scripts' / 'analysis' / 'step06_phylogenomics_utils.R')!r})\n"
                "df <- data.frame(genetic_sample_id = c('S_HQ', 'S_NOISY'), species_standardized = c('Packera aurea', 'Packera dubia'), stringsAsFactors = FALSE)\n"
                f"ev <- evaluate_public_backbone_and_organellar(df, {str(geno_dir)!r}, {str(out_dir)!r})\n"
                "pl <- ev$placement_df\n"
                "stopifnot('containment_gate_pass' %in% names(pl))\n"
                "hq_r <- pl[pl$genetic_sample_id == 'S_HQ', , drop = FALSE]\n"
                "stopifnot(nrow(hq_r) == 1L, as.integer(hq_r$containment_gate_pass[1L]) == 1L, identical(hq_r$best_match_taxon[1L], 'Packera aurea'))\n"
                "stopifnot(identical(ev$public_plastome_match[1L], 'Packera aurea'))\n"
                "stopifnot(is.na(ev$public_plastome_match[2L]) || !nzchar(ev$public_plastome_match[2L]))\n"
                "ns_r <- pl[pl$genetic_sample_id == 'S_NOISY', , drop = FALSE]\n"
                "stopifnot(nrow(ns_r) == 1L, as.integer(ns_r$containment_gate_pass[1L]) == 0L, ns_r$kmer_containment[1L] < 0.80, ns_r$n_shared_11mers[1L] < 80L)\n"
                "# Noisy cont=0.25 placement must blank match and NEVER clear pre-existing cytonuclear_discordant=1\n"
                "aud <- data.frame(\n"
                "  genetic_sample_id = c('S_NOISY', 'S_CLEAN'), species_standardized = c('Packera dubia', 'Packera dubia'),\n"
                "  plastid_clade = c('cp_singleton_1', 'cp_singleton_2'), cytonuclear_discordant = c(1L, 0L),\n"
                "  public_plastome_match = c('Packera aurea', 'Packera dubia'), public_its_match = c('Packera dubia', 'Packera dubia'),\n"
                "  waxy_lineage_match = c('Packera dubia', 'Packera dubia'), geno_predicted_taxon = c('Packera dubia', 'Packera dubia'),\n"
                "  geno_posterior_prob = c(0.96, 0.96), is_geno_admixed = c(0L, 0L), stringsAsFactors = FALSE\n"
                ")\n"
                "noisy_pl <- data.frame(genetic_sample_id = c('S_NOISY', 'S_CLEAN'), locus = c('plastid', 'plastid'), best_match_taxon = c('Packera aurea', 'Packera aurea'), kmer_containment = c(0.25, 0.25), n_shared_11mers = c(30L, 30L), containment_gate_pass = c(0L, 0L), stringsAsFactors = FALSE)\n"
                "rec <- reconcile_step06_cytonuclear_discordance(aud, placement_df = noisy_pl)\n"
                "stopifnot(is.na(rec$public_plastome_match[rec$genetic_sample_id == 'S_NOISY']) || !nzchar(rec$public_plastome_match[rec$genetic_sample_id == 'S_NOISY']))\n"
                "stopifnot(as.integer(rec$cytonuclear_discordant[rec$genetic_sample_id == 'S_NOISY']) == 1L)\n"
                "stopifnot(as.integer(rec$cytonuclear_discordant[rec$genetic_sample_id == 'S_CLEAN']) == 0L)\n"
                "cat('ORGANELLAR_GATE_OK\\n')\n",
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step06", cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Organellar containment gate check failed:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("ORGANELLAR_GATE_OK", proc.stdout)

    def test_shared_site_gdi_eliminates_asymmetric_missingness_bias(self) -> None:
        """Verify shared-site gdi evaluates theta_A, theta_B, and d_XY on identical jointly callable loci."""
        # A0, B0, B1 come from one panmictic population: loci 0..39 are invariant (0.0) where A0 dropped out (NaN),
        # while loci 40..79 are polymorphic (0.5) across all 3 samples. Without shared-site masking, theta_B is
        # diluted by 2x across loci 0..79 while d_XY is evaluated only on loci 40..79, inflating gdi > 0.20.
        mat = np.zeros((3, 80), dtype=np.float64)
        mat[:, 40:] = 0.5
        mat[0, :40] = np.nan
        blk_cols = [cols for cols in np.array_split(np.arange(80), 8)]
        n_blk = len(blk_cols)
        diff_blk = np.zeros((n_blk, 3, 3), dtype=np.float64)
        call_blk = np.zeros((n_blk, 3, 3), dtype=np.float64)
        w_het = np.zeros((n_blk, 3), dtype=np.float64)
        w_call = np.zeros((n_blk, 3), dtype=np.float64)
        for bi, cols in enumerate(blk_cols):
            sub = mat[:, cols]
            ok = np.isfinite(sub)
            ok_f = ok.astype(np.float64)
            xs, omx = np.where(ok, sub, 0.0), np.where(ok, 1.0 - sub, 0.0)
            call_blk[bi] = (ok_f @ ok_f.T) * 50.0
            diff_blk[bi] = xs @ omx.T + omx @ xs.T
            w_call[bi] = np.sum(ok_f, axis=1) * 50.0
            w_het[bi] = np.sum(2.0 * xs * omx, axis=1)
        biased_df, _ = coal_mod.compute_coalescent_gdi_delimitation(
            ["A0", "B0", "B1"], ["Deme_A", "Deme_B", "Deme_B"], diff_blk, call_blk, w_het, w_call, n_boot=50
        )
        sh_ctx = {"mode": "dosage", "mat": mat, "ok_all": np.isfinite(mat), "blk_cols": blk_cols, "site_scale": 50.0}
        gdi_df, _ = coal_mod.compute_coalescent_gdi_delimitation(
            ["A0", "B0", "B1"], ["Deme_A", "Deme_B", "Deme_B"], diff_blk, call_blk, w_het, w_call, n_boot=50, shared_ctx=sh_ctx
        )
        self.assertGreater(float(biased_df[biased_df["taxon"] == "Deme_A"].iloc[0]["coalescent_gdi"]), 0.20)
        row_a = gdi_df[gdi_df["taxon"] == "Deme_A"].iloc[0]
        self.assertEqual(float(row_a["coalescent_gdi"]), 0.0)
        self.assertAlmostEqual(float(row_a["theta_lineage"]), float(row_a["theta_sister"]), places=6)


if __name__ == "__main__":
    unittest.main()
