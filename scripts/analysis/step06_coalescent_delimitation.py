#!/usr/bin/env python3
"""
Module: scripts/analysis/step06_coalescent_delimitation.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Invariant-site Multispecies Coalescent (MSC) & Reticulation (MSci) engine for
  Step 06 phylogenomics (`SVDquartets`, `BPP` control/PHYLIP export + conjugate
  Bayes factors, and locus-bootstrapped `gdi = 1 - exp(-2 * tau / theta)`).
"""

from __future__ import annotations

import argparse
import importlib
import math
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis.step06_coalescent_alignment_utils import (
    _BASE_TO_IDX,
    _IUPAC_ALLELES,
    _clean_id,
    _collapse_alignment_replicates,
    _collapse_dosage_replicates,
    _encode_haplotype_alleles,
    _expand_gdi_samples,
    _group_technical_replicates,
    _infer_taxon,
    _load_sample_metadata,
    _pairwise_site_divergence_blocks,
    _select_models_for_skimmed_fastas,
    compute_svdquartets_support,
    export_bpp_msci_files,
    parse_invariant_alignment,
)
from scripts.analysis.step06_coalescent_models import (
    _apply_evo_distance,
    _resolve_max_workers,
    compute_coalescent_introgression_bf,
    compute_model_corrected_divergence_blocks,
    select_nucleotide_substitution_model,
)
from scripts.core.config import resolve_genomics_config
from scripts.core.logger import setup_logging

logger = setup_logging(name="packera.genomics.coalescent")


def _calc_lineage_theta_and_dxy(
    diff_blk: np.ndarray,
    call_blk: np.ndarray,
    within_het: np.ndarray,
    within_call: np.ndarray,
    idx_a: Sequence[int],
    idx_b: Sequence[int],
    blk_indices: Optional[np.ndarray] = None,
) -> Tuple[float, float, float]:
    """Compute `(theta_A, theta_B, d_XY)` on identical per-site coalescent scale."""
    bi = np.arange(diff_blk.shape[0]) if blk_indices is None else blk_indices
    sub_d, sub_c, sub_wh, sub_wc = diff_blk[bi], call_blk[bi], within_het[bi], within_call[bi]

    def _within_pi(idxs: Sequence[int]) -> float:
        num = den = 0.0
        if len(idxs) == 1:
            num, den = float(np.sum(sub_wh[:, idxs[0]])), float(np.sum(sub_wc[:, idxs[0]]))
        else:
            for p_i in range(len(idxs)):
                for p_j in range(p_i + 1, len(idxs)):
                    ia, ib = idxs[p_i], idxs[p_j]
                    num += float(np.sum(sub_d[:, ia, ib]))
                    den += float(np.sum(sub_c[:, ia, ib]))
        return max(1e-4, num / den) if den > 0 else 1e-4

    theta_a, theta_b = _within_pi(idx_a), _within_pi(idx_b)
    d_num = sum(float(np.sum(sub_d[:, ia, ib])) for ia in idx_a for ib in idx_b)
    d_den = sum(float(np.sum(sub_c[:, ia, ib])) for ia in idx_a for ib in idx_b)
    return theta_a, theta_b, (d_num / d_den) if d_den > 0 else 0.0


def _build_shared_pair_block_stats(
    idx_a: Sequence[int],
    idx_b: Sequence[int],
    diff_blk: np.ndarray,
    call_blk: np.ndarray,
    shared_ctx: Optional[Dict[str, Any]],
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Evaluate per-block `(na, da, nb, db, nxy, dxy)` strictly over sites callable in both lineages."""
    if not shared_ctx:
        return None
    n_blk = diff_blk.shape[0]
    ix_r, ix_c = np.ix_(idx_a, idx_b)
    nxy = np.sum(diff_blk[:, ix_r, ix_c], axis=(1, 2))
    dxy = np.sum(call_blk[:, ix_r, ix_c], axis=(1, 2))
    na, da = np.zeros(n_blk, dtype=np.float64), np.zeros(n_blk, dtype=np.float64)
    nb, db = np.zeros(n_blk, dtype=np.float64), np.zeros(n_blk, dtype=np.float64)

    if shared_ctx.get("mode") == "alignment":
        a1, a2 = shared_ctx["a1"], shared_ctx["a2"]
        ok_all, slices, mspec = shared_ctx["ok_all"], shared_ctx["block_slices"], shared_ctx.get("model_spec")
        joint_ab = np.any(ok_all[idx_a], axis=0) & np.any(ok_all[idx_b], axis=0)

        def _aln_within(idxs: Sequence[int], out_n: np.ndarray, out_d: np.ndarray) -> None:
            for bi, (st, en) in enumerate(slices):
                jm = joint_ab[st:en]
                if not np.any(jm):
                    continue
                x1, x2 = a1[:, st:en], a2[:, st:en]
                ok = ok_all[:, st:en] & jm[None, :]
                if len(idxs) == 1:
                    i0, oki = idxs[0], ok[idxs[0]]
                    wc = float(np.sum(oki))
                    if wc > 0.0:
                        xd = np.where(oki, x1[i0] ^ x2[i0], 0)
                        p1, p2 = np.asarray([np.sum(oki & (xd == 2) & ((x1[i0] & 1) == 0)) / wc]), np.asarray([np.sum(oki & (xd == 2) & ((x1[i0] & 1) == 1)) / wc])
                        qs = np.asarray([np.sum(oki & ((xd & 1) == 1)) / wc])
                        out_n[bi], out_d[bi] = float(_apply_evo_distance(p1, p2, qs, mspec)[0]) * wc, wc
                else:
                    bn = bd = 0.0
                    for p_i, ia in enumerate(idxs):
                        if not np.any(ok[ia]):
                            continue
                        hx = np.stack([x1[ia], x1[ia], x2[ia], x2[ia]], axis=0)
                        for ib in idxs[p_i + 1 :]:
                            vij = ok[ia] & ok[ib]
                            if (cij := float(np.sum(vij))) <= 0.0:
                                continue
                            xor, den = np.where(vij[None, :], hx ^ np.stack([x1[ib], x2[ib], x1[ib], x2[ib]], axis=0), 0), 4.0 * cij
                            p1, p2 = np.asarray([np.sum(vij[None, :] & (xor == 2) & ((hx & 1) == 0)) / den]), np.asarray([np.sum(vij[None, :] & (xor == 2) & ((hx & 1) == 1)) / den])
                            qs = np.asarray([np.sum(vij[None, :] & ((xor & 1) == 1)) / den])
                            bn, bd = bn + float(_apply_evo_distance(p1, p2, qs, mspec)[0]) * cij, bd + cij
                    out_n[bi], out_d[bi] = bn, bd

        _aln_within(idx_a, na, da)
        _aln_within(idx_b, nb, db)
        return na, da, nb, db, nxy, dxy

    if shared_ctx.get("mode") == "dosage":
        mat, ok_all, blk_cols, sc = shared_ctx["mat"], shared_ctx["ok_all"], shared_ctx["blk_cols"], float(shared_ctx.get("site_scale", 50.0))
        joint_ab = np.any(ok_all[idx_a], axis=0) & np.any(ok_all[idx_b], axis=0)

        def _dos_within(idxs: Sequence[int], out_n: np.ndarray, out_d: np.ndarray) -> None:
            for bi, cols in enumerate(blk_cols):
                if cols.size == 0:
                    continue
                jm = joint_ab[cols]
                if not np.any(jm):
                    continue
                sub = np.clip(mat[:, cols], 0.0, 1.0)
                ok = ok_all[:, cols] & jm[None, :]
                if len(idxs) == 1:
                    i0 = idxs[0]
                    oki = ok[i0]
                    xs = np.where(oki, sub[i0], 0.0)
                    out_n[bi], out_d[bi] = float(np.sum(2.0 * xs * np.where(oki, 1.0 - sub[i0], 0.0))), float(np.sum(oki)) * sc
                else:
                    ok_f = ok[idxs].astype(np.float64)
                    xs = np.where(ok[idxs], sub[idxs], 0.0)
                    omx = np.where(ok[idxs], 1.0 - sub[idxs], 0.0)
                    c_m = (ok_f @ ok_f.T) * sc
                    d_m = xs @ omx.T + omx @ xs.T
                    iu = np.triu_indices(len(idxs), k=1)
                    out_n[bi], out_d[bi] = float(np.sum(d_m[iu])), float(np.sum(c_m[iu]))

        _dos_within(idx_a, na, da)
        _dos_within(idx_b, nb, db)
        return na, da, nb, db, nxy, dxy
    return None


def compute_coalescent_gdi_delimitation(
    sample_ids: Sequence[str],
    taxa_labels: Sequence[str],
    diff_blk: np.ndarray,
    call_blk: np.ndarray,
    within_het: np.ndarray,
    within_call: np.ndarray,
    n_boot: int = 500,
    seed: int = 202611,
    substitution_model: str = "JC69",
    bpp_out_path: Optional[Path] = None,
    shared_ctx: Optional[Dict[str, Any]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per-taxon and per-sample `gdi = 1 - exp(-2*tau/theta)` with locus-bootstrap 95% CIs."""
    by_taxon: Dict[str, List[int]] = defaultdict(list)
    for i, tx in enumerate(taxa_labels):
        by_taxon[tx].append(i)
    taxa = sorted(by_taxon.keys())
    n_blk = diff_blk.shape[0]
    rng = np.random.default_rng(seed)

    pair_dxy: Dict[Tuple[str, str], float] = {}
    for ta in taxa:
        for tb in taxa:
            if ta < tb:
                _, _, dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, by_taxon[ta], by_taxon[tb])
                pair_dxy[(ta, tb)] = pair_dxy[(tb, ta)] = dxy

    taxon_rows: List[Dict[str, Any]] = []
    adequacy_rows: List[Dict[str, Any]] = []
    taxon_gdi_map: Dict[str, Dict[str, Any]] = {}
    boot_draws = rng.integers(0, n_blk, size=(max(10, n_boot), n_blk)) if n_blk > 1 else np.zeros((10, 1), dtype=int)

    def _eval_tx_gdi(tx: str) -> Optional[Tuple[str, Dict[str, Any], Dict[str, Any]]]:
        others = [o for o in taxa if o != tx]
        if not others:
            return None
        sister = min(others, key=lambda o: pair_dxy.get((tx, o), 1.0))
        idx_a, idx_b = by_taxon[tx], by_taxon[sister]
        sh_vecs = _build_shared_pair_block_stats(idx_a, idx_b, diff_blk, call_blk, shared_ctx)
        if sh_vecs is not None:
            na_v, da_v, nb_v, db_v, nxy_v, dxy_v = sh_vecs
            th_a = max(1e-4, float(np.sum(na_v) / np.sum(da_v))) if float(np.sum(da_v)) > 0.0 else 1e-4
            th_b = max(1e-4, float(np.sum(nb_v) / np.sum(db_v))) if float(np.sum(db_v)) > 0.0 else 1e-4
            dxy = float(np.sum(nxy_v) / np.sum(dxy_v)) if float(np.sum(dxy_v)) > 0.0 else 0.0
            tau_net = max(0.0, 0.5 * (dxy - 0.5 * (th_a + th_b)))
            gdi_obs = float(1.0 - math.exp(-2.0 * tau_net / max(1e-4, th_a)))
            b_ta = np.where(np.sum(da_v[boot_draws], axis=1) > 0.0, np.maximum(1e-4, np.sum(na_v[boot_draws], axis=1) / np.maximum(1e-12, np.sum(da_v[boot_draws], axis=1))), 1e-4)
            b_tb = np.where(np.sum(db_v[boot_draws], axis=1) > 0.0, np.maximum(1e-4, np.sum(nb_v[boot_draws], axis=1) / np.maximum(1e-12, np.sum(db_v[boot_draws], axis=1))), 1e-4)
            b_dxy = np.where(np.sum(dxy_v[boot_draws], axis=1) > 0.0, np.sum(nxy_v[boot_draws], axis=1) / np.maximum(1e-12, np.sum(dxy_v[boot_draws], axis=1)), 0.0)
            boot_vals = 1.0 - np.exp(-2.0 * np.maximum(0.0, 0.5 * (b_dxy - 0.5 * (b_ta + b_tb))) / np.maximum(1e-4, b_ta))
        else:
            th_a, th_b, dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, idx_a, idx_b)
            tau_net = max(0.0, 0.5 * (dxy - 0.5 * (th_a + th_b)))
            gdi_obs = float(1.0 - math.exp(-2.0 * tau_net / max(1e-4, th_a)))
            boot_vals = [
                1.0 - math.exp(-2.0 * max(0.0, 0.5 * (b_dxy - 0.5 * (b_ta + b_tb))) / max(1e-4, b_ta))
                for b_row in boot_draws
                for b_ta, b_tb, b_dxy in (_calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, idx_a, idx_b, blk_indices=b_row),)
            ]
        ci_low, ci_high = float(np.percentile(boot_vals, 2.5)), float(np.percentile(boot_vals, 97.5))

        ix_r, ix_c = np.ix_(idx_a, idx_b)
        call_mean_vec = np.mean(call_blk[:, ix_r, ix_c], axis=(1, 2))
        diff_mean_vec = np.mean(diff_blk[:, ix_r, ix_c], axis=(1, 2))
        bf_info = compute_coalescent_introgression_bf(np.divide(diff_mean_vec, np.maximum(1.0, call_mean_vec)), call_mean_vec, bpp_out_path=bpp_out_path)
        bpp_bf = float(bf_info["bf_10"])
        gdi_cls = "Single_Lineage_Oversplit" if gdi_obs < 0.20 else ("Distinct_Coalescent_Lineage" if gdi_obs > 0.70 else "Ambiguous_Incomplete_Sorting")

        rec = {
            "taxon": tx, "nearest_sister_taxon": sister, "n_samples": len(idx_a), "n_loci_blocks": n_blk,
            "theta_lineage": round(th_a, 6), "theta_sister": round(th_b, 6), "d_xy": round(dxy, 6), "tau_net": round(tau_net, 6),
            "coalescent_gdi": round(gdi_obs, 4), "gdi_ci_low": round(ci_low, 4), "gdi_ci_high": round(ci_high, 4),
            "gdi_classification": gdi_cls, "bpp_introgression_bf": bpp_bf, "substitution_model": substitution_model, "bf_method": bf_info["bf_method"],
        }
        ad_row = {
            "taxon": tx, "nearest_sister_taxon": sister, "substitution_model": substitution_model, "bf_method": bf_info["bf_method"],
            "bic_m0": bf_info["bic_m0"], "bic_m1": bf_info["bic_m1"], "log_bf_10": bf_info["log_bf_10"], "bf_10": bf_info["bf_10"],
            "w_introgr": bf_info["w_introgr"], "tau_introgr": bf_info["tau_introgr"], "tau_sp": bf_info["tau_sp"],
            "coalescent_gdi": rec["coalescent_gdi"], "gdi_classification": gdi_cls,
        }
        return tx, rec, ad_row

    nw_tx = min(_resolve_max_workers(), len(taxa)) if len(taxa) >= 2 else 1
    if nw_tx > 1:
        with ThreadPoolExecutor(max_workers=nw_tx) as pool:
            tx_outs = list(pool.map(_eval_tx_gdi, taxa))
    else:
        tx_outs = [_eval_tx_gdi(tx) for tx in taxa]
    for item in tx_outs:
        if item is None:
            continue
        tx, rec, ad_row = item
        taxon_rows.append(rec)
        taxon_gdi_map[tx] = rec
        adequacy_rows.append(ad_row)

    sample_rows = [
        {
            "genetic_sample_id": sid, "lineage_taxon": tx, "nearest_sister_taxon": info.get("nearest_sister_taxon"),
            "coalescent_gdi": info.get("coalescent_gdi"), "gdi_ci_low": info.get("gdi_ci_low"),
            "gdi_ci_high": info.get("gdi_ci_high"), "gdi_classification": info.get("gdi_classification"),
            "bpp_introgression_bf": info.get("bpp_introgression_bf"),
        }
        for sid, tx in zip(sample_ids, taxa_labels)
        for info in (taxon_gdi_map.get(tx, {}),)
    ]
    taxon_df = pd.DataFrame(taxon_rows)
    taxon_df.attrs["adequacy_df"] = pd.DataFrame(adequacy_rows)
    return taxon_df, pd.DataFrame(sample_rows)


def run_coalescent_delimitation(
    genomics_dir: Path,
    out_dir: Path,
    parsed_samples_csv: Optional[Path] = None,
    dosage_tsv: Optional[Path] = None,
    block_len: Optional[int] = None,
    n_boot: Optional[int] = None,
    focal_taxon: Optional[str] = None,
    intraspecific_mode: bool = False,
    cluster_csv: Optional[Path] = None,
    gdi_source: str = "nextrad_msc",
    alignment_file: Optional[Path] = None,
) -> Dict[str, Any]:
    """Execute invariant-site SVDquartets, BPP export, and gdi delimitation into `out_dir`."""
    gcfg = resolve_genomics_config()
    eff_block_len = int(gcfg["coalescent_block_len"]) if block_len is None else int(block_len)
    eff_n_boot = int(gcfg["coalescent_n_boot"]) if n_boot is None else int(n_boot)

    discover_fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").discover_genomic_inputs
    assets = discover_fn(genomics_dir)
    align_cands = [alignment_file] if (alignment_file and alignment_file.is_file()) else sorted(
        assets.alignment_files, key=lambda p: (0 if p.name.lower().endswith((".phylip", ".phy", ".nex", ".nexus")) else 1, p.name)
    )
    if not align_cands and not (dosage_tsv and dosage_tsv.exists()):
        return {"status": "skipped_empty", "n_samples": 0}

    pass_set, sample_taxon_map, rep_pair_map, hq_reads_map, miss_rate_map = _load_sample_metadata(genomics_dir, parsed_samples_csv)
    cluster_map: Dict[str, str] = {}
    cl_src = cluster_csv if (cluster_csv and cluster_csv.is_file()) else parsed_samples_csv
    if cl_src and cl_src.is_file():
        cdf = pd.read_csv(cl_src, dtype=str, keep_default_na=False)
        for _, r in cdf.iterrows():
            if (sid_c := str(r.get("genetic_sample_id") or "").strip()) and (gc := str(r.get("genetic_cluster") or "").strip()):
                cluster_map[sid_c] = gc
    out_dir.mkdir(parents=True, exist_ok=True)
    gdi_csv, gdi_samp_csv = out_dir / "coalescent_gdi_delimitation.csv", out_dir / "coalescent_gdi_samples.csv"
    svdq_csv, model_sel_csv, adequacy_csv = out_dir / "svdquartets_support.csv", out_dir / "substitution_model_selection.csv", out_dir / "coalescent_model_adequacy.csv"

    def _apply_intraspecific_demes(
        gdi_df: pd.DataFrame, gdi_s_df: pd.DataFrame, adeq_df: pd.DataFrame,
        sids_in: Sequence[str], taxa_in: Sequence[str], diff_blk: np.ndarray, call_blk: np.ndarray,
        within_het: np.ndarray, within_call: np.ndarray, nb: int, smod: str, bpp_p: Optional[Path],
        sh_ctx: Optional[Dict[str, Any]] = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        if not intraspecific_mode or not focal_taxon or not cluster_map:
            return gdi_df, gdi_s_df, adeq_df
        f_low = focal_taxon.strip().lower()
        deme_lbls = [cluster_map.get(sid, tx) if tx.strip().lower() == f_low else tx for sid, tx in zip(sids_in, taxa_in)]
        focal_demes = {lbl for lbl, tx in zip(deme_lbls, taxa_in) if tx.strip().lower() == f_low and lbl.startswith("Geno_Cluster_")}
        if len(focal_demes) < 2:
            return gdi_df, gdi_s_df, adeq_df
        d_df, d_s_df = compute_coalescent_gdi_delimitation(sids_in, deme_lbls, diff_blk, call_blk, within_het, within_call, n_boot=nb, substitution_model=smod, bpp_out_path=bpp_p, shared_ctx=sh_ctx)
        d_adeq = d_df.attrs.pop("adequacy_df", pd.DataFrame())
        base_clean = gdi_df.copy(); base_clean.attrs.clear()
        d_sub = d_df[d_df["taxon"].isin(focal_demes)].copy() if not d_df.empty else pd.DataFrame()
        d_sub.attrs.clear()
        out_gdi = pd.concat([base_clean, d_sub], ignore_index=True) if not d_sub.empty else base_clean
        ad_clean = adeq_df.copy(); ad_clean.attrs.clear()
        ad_sub = d_adeq[d_adeq["taxon"].isin(focal_demes)].copy() if not d_adeq.empty else pd.DataFrame()
        ad_sub.attrs.clear()
        out_adeq = pd.concat([ad_clean, ad_sub], ignore_index=True) if not ad_sub.empty else ad_clean
        out_s = gdi_s_df.copy(); out_s.attrs.clear()
        if not d_s_df.empty and not out_s.empty:
            d_s_map = {str(r["genetic_sample_id"]): r for _, r in d_s_df.iterrows() if str(r.get("lineage_taxon") or "") in focal_demes}
            for r_i, s_id in enumerate(out_s["genetic_sample_id"].astype(str).tolist()):
                if (dr := d_s_map.get(s_id)) is not None:
                    for col_s in ("lineage_taxon", "nearest_sister_taxon", "coalescent_gdi", "gdi_ci_low", "gdi_ci_high", "gdi_classification", "bpp_introgression_bf"):
                        if col_s in out_s.columns and col_s in dr:
                            out_s.at[r_i, col_s] = dr[col_s]
        return out_gdi, out_s, out_adeq

    if align_cands:
        primary_path = align_cands[0]
        sids, seqs = parse_invariant_alignment(primary_path)
        if pass_set is not None and len(kept := [(sid, sq) for sid, sq in zip(sids, seqs) if sid in pass_set]) >= 4:
            sids, seqs = [x[0] for x in kept], [x[1] for x in kept]
        orig_sids = list(sids)
        taxa_labels = [_infer_taxon(sid, sample_taxon_map) for sid in sids]
        sids, seqs, taxa_labels, group_members = _collapse_alignment_replicates(sids, seqs, taxa_labels, rep_pair_map, hq_reads_map, miss_rate_map)
        n_sites = min(len(s) for s in seqs) if seqs else 0
        eff_blk = eff_block_len if n_sites >= eff_block_len * 2 else max(20, n_sites // 4 or n_sites)
        block_slices = [(st, min(n_sites, st + eff_blk)) for st in range(0, n_sites, eff_blk) if (min(n_sites, st + eff_blk) - st) >= 15]

        if paup_bin := shutil.which("paup"):
            logger.info("Detected external PAUP* binary (%s); running SVDquartets check alongside Python engine.", paup_bin)
        if bpp_bin := shutil.which("bpp"):
            logger.info("Detected external BPP binary (%s); generated control file for MSci.", bpp_bin)

        primary_spec = select_nucleotide_substitution_model(dict(zip(sids, seqs)), sids, alignment_name=primary_path.name)
        model_rows = list(primary_spec.get("table_rows") or []) + _select_models_for_skimmed_fastas(genomics_dir, [primary_path])
        pd.DataFrame(model_rows).to_csv(model_sel_csv, index=False)

        best_model = str(primary_spec.get("best_model") or "JC69")
        a1, a2 = _encode_haplotype_alleles(seqs)
        ok_raw = a1 >= 0
        nc_site = np.sum(ok_raw, axis=0)
        ok_non_para = ok_raw & (~((nc_site >= 4) & (np.sum(ok_raw & (a1 != a2), axis=0) > 0.60 * nc_site)))[None, :]
        aln_sh_ctx = {"mode": "alignment", "a1": a1, "a2": a2, "ok_all": ok_non_para, "block_slices": block_slices, "model_spec": primary_spec}
        diff_blk, call_blk, within_het, within_call = compute_model_corrected_divergence_blocks(a1, a2, block_slices, model_spec=primary_spec)
        svdq_lbls = list(taxa_labels)
        if len(set(svdq_lbls)) < 4 and intraspecific_mode and focal_taxon and cluster_map:
            f_low_q = focal_taxon.strip().lower()
            cand_q = [cluster_map.get(sid, tx) if tx.strip().lower() == f_low_q else tx for sid, tx in zip(sids, taxa_labels)]
            if len(set(cand_q)) >= 4:
                svdq_lbls = cand_q
        svdq_df = compute_svdquartets_support(a1, a2, svdq_lbls, sids, n_boot=200, seed=202611)
        export_bpp_msci_files(sids, taxa_labels, seqs, block_slices, out_dir, substitution_model=best_model)
        gdi_df, gdi_s_df = compute_coalescent_gdi_delimitation(
            sids, taxa_labels, diff_blk, call_blk, within_het, within_call,
            n_boot=eff_n_boot, substitution_model=best_model, bpp_out_path=out_dir / "bpp_msci_out.txt", shared_ctx=aln_sh_ctx,
        )
        adequacy_df = gdi_df.attrs.pop("adequacy_df", pd.DataFrame())
        gdi_df, gdi_s_df, adequacy_df = _apply_intraspecific_demes(gdi_df, gdi_s_df, adequacy_df, sids, taxa_labels, diff_blk, call_blk, within_het, within_call, eff_n_boot, best_model, out_dir / "bpp_msci_out.txt", sh_ctx=aln_sh_ctx)
        gdi_df["gdi_source"] = gdi_source
        gdi_s_df = _expand_gdi_samples(gdi_s_df, group_members, orig_sids)
        gdi_s_df["gdi_source"] = gdi_source
        svdq_df.to_csv(svdq_csv, index=False)
        gdi_df.to_csv(gdi_csv, index=False)
        gdi_s_df.to_csv(gdi_samp_csv, index=False)
        adequacy_df.to_csv(adequacy_csv, index=False)
        return {
            "status": "completed_alignment", "n_samples": len(sids), "n_blocks": len(block_slices),
            "best_substitution_model": best_model, "gdi_csv": str(gdi_csv), "svdq_csv": str(svdq_csv),
            "model_selection_csv": str(model_sel_csv), "adequacy_csv": str(adequacy_csv),
        }

    if dosage_tsv and dosage_tsv.exists():
        ddf = pd.read_csv(dosage_tsv, sep="\t")
        sids = [str(x) for x in ddf["genetic_sample_id"].tolist()]
        mat = ddf.drop(columns=["genetic_sample_id"]).to_numpy(dtype=np.float64)
        orig_sids = list(sids)
        taxa_labels = [_infer_taxon(sid, sample_taxon_map) for sid in sids]
        sids, mat, taxa_labels, group_members = _collapse_dosage_replicates(sids, mat, taxa_labels, rep_pair_map, hq_reads_map, miss_rate_map)
        n_s, n_l = mat.shape
        n_blk = min(32, max(4, n_l // 4)) if n_l >= 16 else (min(16, n_l) if n_l >= 4 else max(1, n_l))
        blk_cols = [cols for cols in np.array_split(np.arange(n_l), n_blk) if cols.size > 0] if n_l > 0 else [np.arange(0)]
        n_blk = len(blk_cols)
        site_scale = 50.0
        diff_blk = np.zeros((n_blk, n_s, n_s), dtype=np.float64)
        call_blk = np.zeros((n_blk, n_s, n_s), dtype=np.float64)
        w_het = np.zeros((n_blk, n_s), dtype=np.float64)
        w_call = np.zeros((n_blk, n_s), dtype=np.float64)
        for bi, cols in enumerate(blk_cols):
            if cols.size == 0:
                call_blk[bi], w_call[bi] = 1.0, 1.0
                continue
            sub = np.clip(mat[:, cols], 0.0, 1.0)
            ok = np.isfinite(sub)
            ok_f = ok.astype(np.float64)
            xs, omx = np.where(ok, sub, 0.0), np.where(ok, 1.0 - sub, 0.0)
            call_blk[bi] = (ok_f @ ok_f.T) * site_scale
            diff_blk[bi] = xs @ omx.T + omx @ xs.T
            w_call[bi] = np.sum(ok_f, axis=1) * site_scale
            w_het[bi] = np.sum(2.0 * xs * omx, axis=1)
        dos_sh_ctx = {"mode": "dosage", "mat": mat, "ok_all": np.isfinite(mat), "blk_cols": blk_cols, "site_scale": site_scale}
        gdi_df, gdi_s_df = compute_coalescent_gdi_delimitation(sids, taxa_labels, diff_blk, call_blk, w_het, w_call, n_boot=eff_n_boot, shared_ctx=dos_sh_ctx)
        adequacy_df = gdi_df.attrs.pop("adequacy_df", pd.DataFrame())
        gdi_df, gdi_s_df, adequacy_df = _apply_intraspecific_demes(gdi_df, gdi_s_df, adequacy_df, sids, taxa_labels, diff_blk, call_blk, w_het, w_call, eff_n_boot, "JC69", None, sh_ctx=dos_sh_ctx)
        gdi_df["gdi_source"] = gdi_source
        gdi_s_df = _expand_gdi_samples(gdi_s_df, group_members, orig_sids)
        gdi_s_df["gdi_source"] = gdi_source
        gdi_df.to_csv(gdi_csv, index=False); gdi_s_df.to_csv(gdi_samp_csv, index=False); adequacy_df.to_csv(adequacy_csv, index=False)
        pd.DataFrame(columns=["quartet_id", "best_split", "svd_score_best", "svd_score_second", "quartet_support", "quartet_boot_support", "quartet_p_value", "n_sites_used"]).to_csv(svdq_csv, index=False)
        return {"status": "completed_dosage_fallback", "n_samples": n_s, "n_blocks": n_blk, "gdi_csv": str(gdi_csv)}

    return {"status": "skipped_empty", "n_samples": 0}


def main(argv: Optional[Sequence[str]] = None) -> int:
    gcfg = resolve_genomics_config()
    parser = argparse.ArgumentParser(description="Invariant-site coalescent delimitation (SVDquartets, BPP, gdi).")
    parser.add_argument("--genomics-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--parsed-samples", type=Path, default=None)
    parser.add_argument("--dosage-tsv", type=Path, default=None)
    parser.add_argument("--alignment-file", type=Path, default=None)
    parser.add_argument("--focal-taxon", type=str, default=None)
    parser.add_argument("--intraspecific-mode", action="store_true", default=False)
    parser.add_argument("--cluster-csv", type=Path, default=None)
    parser.add_argument("--gdi-source", type=str, default="nextrad_msc")
    parser.add_argument("--block-len", type=int, default=int(gcfg["coalescent_block_len"]))
    parser.add_argument("--n-boot", type=int, default=int(gcfg["coalescent_n_boot"]))
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_coalescent_delimitation(
        args.genomics_dir, args.out_dir, parsed_samples_csv=args.parsed_samples, dosage_tsv=args.dosage_tsv,
        block_len=args.block_len, n_boot=args.n_boot, focal_taxon=args.focal_taxon,
        intraspecific_mode=bool(args.intraspecific_mode), cluster_csv=args.cluster_csv,
        gdi_source=args.gdi_source, alignment_file=args.alignment_file,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
