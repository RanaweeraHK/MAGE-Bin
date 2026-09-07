from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data

torch.serialization.add_safe_globals([Data])

# ---- tier config (VB_TIER=low|medium|high, default high; see README's tier
# table) ----
TIER = os.environ.get("VB_TIER", "high")
if TIER not in ("low", "medium", "high"):
    raise SystemExit(f"error: VB_TIER must be low|medium|high (got {TIER!r})")
TIER_SUFFIX = {"low": "_low", "medium": "_medium", "high": "_high"}[TIER]

SIM = os.environ.get("VB_SIM", "iss")
ASSEMBLER = os.environ.get("VB_ASM", "metaspades")
SIM_TAG = "" if SIM == "iss" else f"__{SIM}"
RUN_NAME = f"multibiome{TIER_SUFFIX}{SIM_TAG}__coasm__{ASSEMBLER}__noviral"


REPO_HOME = Path(__file__).resolve().parent.parent.parent
OUT_DIR = (REPO_HOME
           / f"Data/Processed_data/{RUN_NAME}")
POOL_DIR = REPO_HOME / f"Data/work/multibiome{TIER_SUFFIX}/genome_pool"

STRONG_MULTI_SPECIES_FRAC = 0.10
USABLE_MULTI_SPECIES_FRAC = 0.04


def tier(labelled_mask, labels, name):
    lab = labels[labelled_mask]
    uniq, cnt = np.unique(lab, return_counts=True)
    multi = {u for u, c in zip(uniq, cnt) if c >= 2}
    n_in_multi = sum(c for u, c in zip(uniq, cnt) if u in multi)
    # "species" is already plural; only "strain" takes the s.
    plural = name if name.endswith("s") else f"{name}s"
    print(f"\n=== {name} tier ===")
    print(f"  distinct {plural} recovered : {len(uniq)}")
    denom = len(uniq) or 1
    print(f"  MULTI-contig {plural}        : {len(multi)} ({100*len(multi)/denom:.0f}%)")
    lab_total = labelled_mask.sum() or 1
    print(f"  contigs in multi-{name}     : {n_in_multi} "
          f"({100*n_in_multi/lab_total:.0f}% of labelled)")
    return multi


def main() -> None:
    graph_path = OUT_DIR / "viral_graph.pt"
    meta_path = OUT_DIR / "contig_metadata.tsv"
    smap_path = POOL_DIR / "strain_map.tsv"
    for p in (graph_path, meta_path, smap_path):
        if not p.is_file():
            raise SystemExit(f"error: {p} not found - run the earlier phases first")

    d = torch.load(graph_path, weights_only=False)
    meta = pd.read_csv(meta_path, sep="\t").sort_values("node_index").reset_index(drop=True)
    smap = pd.read_csv(smap_path, sep="\t")
    strain2species = dict(zip(smap.strain_id, smap.species_id))
    strain2biome = dict(zip(smap.strain_id, smap.biome))

    n_pool_species = smap.species_id.nunique()
    n = len(meta)
    print(f"tier                     : {TIER}  ({n_pool_species} pool species, "
          f"{len(smap)} strains)")

    strain_lbl = meta["genome_label"].astype(str).to_numpy()
    labelled = np.array([s in strain2species for s in strain_lbl])
    species_lbl = np.array([strain2species.get(s, "NA") for s in strain_lbl])
    biome_lbl = np.array([strain2biome.get(s, "NA") for s in strain_lbl])
    lengths = meta["length"].to_numpy()

    print(f"contigs total            : {n}")
    print(f"labelled (map to genome) : {labelled.sum()}  ({100*labelled.mean():.0f}%)")
    print(f"contig length (bp)       : min={lengths.min():,} "
          f"median={int(np.median(lengths)):,} max={lengths.max():,}")
    for lo, hi, name in [(0, 5000, "<5kb"), (5000, 10000, "5-10kb"),
                          (10000, 10**9, ">10kb")]:
        c = ((lengths >= lo) & (lengths < hi)).sum()
        print(f"    {name:<7}: {c} ({100*c/n:.0f}%)")

    sp_multi = tier(labelled, species_lbl, "species")
    tier(labelled, strain_lbl, "strain")

    # GFA reachability at species tier.
    ei = d.edge_index.numpy()
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for a, b in ei.T:
        if a != b:
            g.add_edge(int(a), int(b))
    comp = np.full(n, -1)
    for c, cc in enumerate(nx.connected_components(g)):
        for i in cc:
            comp[i] = c
    idx_by_species: dict[str, list[int]] = {}
    for i in range(n):
        if labelled[i] and species_lbl[i] in sp_multi:
            idx_by_species.setdefault(species_lbl[i], []).append(i)
    sp = rr = 0
    for _sp_id, nodes in idx_by_species.items():
        for a, bb in combinations(nodes, 2):
            sp += 1
            rr += int(comp[a] == comp[bb])
    print("\n=== assembly graph (species tier) ===")
    print(f"  GFA edges                  : {g.number_of_edges()}")
    print(f"  same-species pairs (multi) : {sp}")
    print(f"  reachable in GFA           : {rr} ({100*rr/max(sp,1):.0f}% -- merge ceiling)")

    print("\n=== per-biome (species tier) ===")
    print(f"  {'biome':<11}{'contigs':>8}{'species':>8}{'multi-sp':>9}{'%contigs in multi':>18}")
    for b in sorted(set(biome_lbl[labelled])):
        mask = labelled & (biome_lbl == b)
        sub_sp = species_lbl[mask]
        uniq, cnt = np.unique(sub_sp, return_counts=True)
        multi = {u for u, c in zip(uniq, cnt) if c >= 2}
        n_in_multi = sum(c for u, c in zip(uniq, cnt) if u in multi)
        print(f"  {b:<11}{mask.sum():>8}{len(uniq):>8}{len(multi):>9}"
              f"{100*n_in_multi/max(mask.sum(),1):>17.0f}%")

    multi_frac = len(sp_multi) / max(n_pool_species, 1)
    verdict = ("STRONG benchmark" if multi_frac >= STRONG_MULTI_SPECIES_FRAC
               else ("usable" if multi_frac >= USABLE_MULTI_SPECIES_FRAC
                     else "too easy"))
    print(f"\nmulti-contig species     : {len(sp_multi)} / {n_pool_species} pool "
          f"species ({100*multi_frac:.0f}%)")
    print("VERDICT:", verdict)


if __name__ == "__main__":
    main()
