#!/usr/bin/env python3
"""Newick/NEXUS tree parsing with collision-free unique tip IDs and outgroup-exempt MAD rogue screening."""
from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Callable, Dict, List, Mapping, Optional, Set, Tuple, Union

import numpy as np

try:
    from scripts.data_prep.step06_synonym_utils import (
        CANONICAL_SENECIO_PACKERA_SYNONYMS,
        extract_species_from_seed_key,
        normalize_step06_taxon_name,
    )
except ImportError:
    from step06_synonym_utils import (  # type: ignore[no-redef]
        CANONICAL_SENECIO_PACKERA_SYNONYMS,
        extract_species_from_seed_key,
        normalize_step06_taxon_name,
    )

_NWK_TOK_RE = re.compile(r"([(),:;])")
_OUTGROUP_GENERA_RE = re.compile(
    r"(?:^|[^A-Za-z])("
    r"Erechtites|Jacobaea|Senecio|Tephroseris|Emilia|Pericallis|Cineraria|Doronicum|"
    r"Arnoglossum|Cacalia|Ligularia|Roldana|Barkleyanthus|Psacaliopsis|Psacalium|"
    r"Crassocephalum|Gynura|Arrhenechthites|Telanthophora|Digitacalia|Pippenalia|"
    r"Blennosperma|Petasites|Pittocaulon|Robinsonecio|Tussilago|Sinacalia|Cacaliopsis|"
    r"Brachyglottis|Cotula|Artemisia|Achillea|Anthemis|Elekmania|Rugelia|Luina"
    r")(?:[\s_]|$)",
    re.IGNORECASE,
)


def _std_taxon_tree(raw: Any) -> str:
    txt = str(raw or "").strip()
    if not txt or txt.upper() in {"NA", "NAN", "NONE", "NULL", "UNKNOWN"}:
        return ""
    txt_clean = re.sub(r"__tip\d+$", "", txt)
    txt_sp = re.sub(r"\s+", " ", txt_clean.replace("_", " ")).strip()
    low = txt_sp.lower()
    for syn, canon in CANONICAL_SENECIO_PACKERA_SYNONYMS.items():
        if re.search(rf"(?:^|[^a-z]){re.escape(syn)}(?:[^a-z]|$)", low):
            return canon
    try:
        from scripts.core.harvester import standardize_packera_taxon
        if std := standardize_packera_taxon(txt_sp):
            return std
    except Exception:
        pass
    if ex := extract_species_from_seed_key(txt_clean):
        return ex
    return normalize_step06_taxon_name(txt_sp)[1] or ""


def _is_outgroup_tip(
    tip_id: str, norm_taxon: Optional[str] = None, *, has_packera_ingroup: bool = False,
) -> bool:
    """Return True iff `tip_id`/`norm_taxon` represents a non-Packera outgroup exempted from ingroup MAD pruning."""
    tx = str(norm_taxon or "").strip() or _std_taxon_tree(tip_id)
    if re.match(r"^(?:Packera|P\.)(?:\s|_|$)", tx, re.IGNORECASE):
        return False
    raw_clean = re.sub(r"[|_]+", " ", re.sub(r"__tip\d+$", "", str(tip_id or ""))).strip()
    if _OUTGROUP_GENERA_RE.search(tx) or _OUTGROUP_GENERA_RE.search(raw_clean):
        return True
    return bool(has_packera_ingroup)


def _extract_nexus_tree_and_translate(text: str) -> Tuple[str, Dict[str, str]]:
    s = text.strip()
    tmap: Dict[str, str] = {}
    if "#nexus" in s[:300].lower():
        if tm := re.search(r"\btranslate\b(.*?);", s, re.I | re.S):
            for entry in tm.group(1).split(","):
                if parts := entry.strip().split(None, 1):
                    if len(parts) == 2:
                        tmap[parts[0].strip()] = parts[1].strip().strip("'\"")
        m = re.search(r"\btree\b[^=]*=\s*(?:\[[^\]]*\]\s*)?(\([^;]+;)", s, re.I | re.S)
        s = m.group(1) if m else ""
    return re.sub(r"\[[^\]]*\]", "", s), tmap


def _parse_newick_tips_and_dists(
    text: str, norm_fn: Optional[Callable[[Any], str]] = None,
) -> Tuple[Dict[str, float], Dict[str, float], Dict[Tuple[str, str], float], Dict[str, List[str]]]:
    """Parse Newick/NEXUS tree keying tips by unique tip IDs so conspecific accessions never overwrite each other."""
    s, tmap = _extract_nexus_tree_and_translate(text)
    if "(" not in s or ";" not in s:
        return {}, {}, {}, {}
    tokens = [t.strip() for t in _NWK_TOK_RE.split(s[: s.index(";") + 1]) if t.strip()]
    children: Dict[int, List[int]] = defaultdict(list)
    blen, label, supp, stack = {}, {}, {}, []
    nid, cur, i = 0, -1, 0
    while i < len(tokens):
        tk = tokens[i]
        if tk == "(":
            nid += 1; (children[stack[-1]].append(nid) if stack else None); stack.append(nid); cur = nid
        elif tk == ",": cur = -1
        elif tk == ")": cur = stack.pop() if stack else -1
        elif tk == ":":
            if i + 1 < len(tokens) and cur > 0:
                try: blen[cur] = max(0.0, float(tokens[i + 1]))
                except ValueError: blen[cur] = 0.0
                i += 1
        elif tk == ";": break
        elif cur in children and cur > 0:
            try: supp[cur] = float(tk)
            except ValueError: pass
        else:
            nid += 1; (children[stack[-1]].append(nid) if stack else None)
            raw_tk = tk.strip("'\""); label[nid] = tmap.get(raw_tk, raw_tk); cur = nid
        i += 1

    resolve_tx = norm_fn or _std_taxon_tree
    parent_of = {ch: pa for pa, chs in children.items() for ch in chs}
    tip_lens, tip_supp, tip_hdrs, tip_nodes_by_tx, node_anc = {}, {}, defaultdict(list), defaultdict(list), {}
    for nd, raw_lbl in label.items():
        if not (tx := resolve_tx(raw_lbl)):
            continue
        uid = raw_lbl if raw_lbl not in tip_lens else f"{raw_lbl}__tip{nd}"
        tip_lens[uid], tip_supp[uid] = blen.get(nd, 0.0), supp.get(parent_of.get(nd, 0), 0.90)
        tip_hdrs[tx].append(uid); tip_nodes_by_tx[tx].append(nd)
        dmap, acc, c = {}, 0.0, nd
        while c in parent_of:
            dmap[c] = acc; acc += blen.get(c, 0.01); c = parent_of[c]
        dmap[c] = acc; node_anc[nd] = dmap

    coph: Dict[Tuple[str, str], float] = {}
    taxa = list(tip_nodes_by_tx.keys())
    for ia, ta in enumerate(taxa):
        for tb in taxa[ia + 1 :]:
            dists = [
                min(node_anc[na][a] + node_anc[nb][a] for a in common)
                for na in tip_nodes_by_tx[ta] for nb in tip_nodes_by_tx[tb]
                if (common := [a for a in node_anc[na] if a in node_anc[nb]])
            ]
            if dists:
                coph[(ta, tb)] = coph[(tb, ta)] = min(dists)
    return tip_lens, tip_supp, coph, tip_hdrs


def _mad_screen_tips(
    tip_lens: Mapping[str, Union[float, Tuple[str, float]]],
    mad_z: float = 5.0,
    tip_taxa: Optional[Mapping[str, str]] = None,
) -> Set[str]:
    """Screen ingroup tips via median/MAD while unconditionally exempting non-Packera outgroups."""
    if not tip_lens:
        return set()
    parsed = [
        (uid, (str(v[0]) if isinstance(v, tuple) and len(v) == 2 else str((tip_taxa or {}).get(uid) or "")) or _std_taxon_tree(uid),
         float(v[1]) if isinstance(v, tuple) and len(v) == 2 else float(v))
        for uid, v in tip_lens.items()
    ]
    has_pk = any(re.match(r"^(?:Packera|P\.)(?:\s|_|$)", tx, re.IGNORECASE) for _, tx, _ in parsed)
    outgroup_uids = {uid for uid, tx, _ in parsed if _is_outgroup_tip(uid, tx, has_packera_ingroup=has_pk)}
    ingroup_items = [(uid, dv) for uid, tx, dv in parsed if uid not in outgroup_uids]
    if len(ingroup_items) < 6:
        return set(tip_lens.keys())
    vals = np.array([v for _, v in ingroup_items], dtype=np.float64)
    med = float(np.median(vals))
    mad = max(1e-6, 1.4826 * float(np.median(np.abs(vals - med))))
    return outgroup_uids | {uid for uid, v in ingroup_items if not (((v - med) / mad) > mad_z and v > 3.0 * max(med, 1e-4))}
