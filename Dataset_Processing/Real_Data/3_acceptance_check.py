#!/usr/bin/env python
# Phase 3 - is this real dataset a usable binning problem?
#
#   VB_DATASET=aloha_illumina python 3_acceptance_check.py
#

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data

torch.serialization.add_safe_globals([Data])

# ---- config (kept local so this script runs standalone) ----
# .parent.parent: this script lives in Dataset_Processing/Real_Data/, so the
# repo root is two levels up.
HERE = Path(__file__).resolve().parent
REPO_HOME = HERE.parent.parent
REGISTRY = HERE / "datasets.tsv"


COHORT_MIN_BP = 2000

STRONG_JOINABLE_FRAC = 0.10
USABLE_JOINABLE_FRAC = 0.04


def load_registry(path: Path) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(
            (line for line in fh if line.strip() and not line.startswith("#")),
            delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise SystemExit(f"error: {path} has no header row")
        for parts in reader:
            row = dict(zip(header, parts))
            rows[row["dataset"]] = row
    return rows


DATASET = os.environ.get("VB_DATASET", "")
_registry = load_registry(REGISTRY)
if DATASET not in _registry:
    raise SystemExit(
        f"error: VB_DATASET must name a row in {REGISTRY}\n"
        f"       got {DATASET!r}; known: {', '.join(sorted(_registry)) or '(none)'}")
_row = _registry[DATASET]
_declared_assembler = _row.get("assembler", "-")
ASSEMBLER = (_declared_assembler if _declared_assembler not in ("-", "")
             else {"short": "metaspades", "long": "metaflye"}[_row["read_type"]])
VIRAL_ID_TAG = "noviral" if _row["viral_enriched"] == "yes" else "genomad"
RUN_NAME = f"{DATASET}__coasm__{ASSEMBLER}__{VIRAL_ID_TAG}"
OUT_DIR = REPO_HOME / f"Data/Processed_data/{RUN_NAME}"


def main() -> None:
    graph_path = OUT_DIR / "viral_graph.pt"
    meta_path = OUT_DIR / "contig_metadata.tsv"
    config_path = OUT_DIR / "config.json"
    for p in (graph_path, meta_path, config_path):
        if not p.is_file():
            raise SystemExit(f"error: {p} not found - run phase 2 first")

    config = json.loads(config_path.read_text())
    data = torch.load(graph_path, weights_only=False)
    meta = (pd.read_csv(meta_path, sep="\t")
            .sort_values("node_index").reset_index(drop=True))
    n = len(meta)
    lengths = meta["length"].to_numpy()
    cohort = lengths >= COHORT_MIN_BP

    print("=" * 68)
    print(f"acceptance check  ::  {RUN_NAME}")
    print("=" * 68)
    print(f"dataset          : {DATASET} ({_row['read_type']} reads, {ASSEMBLER})")
    print(f"accessions       : {_row['accessions']}")
    print("ground truth     : none - reference-free checks only")

    # A real dataset must NOT claim ground truth: everything downstream keys off
    # this, so it is checked rather than assumed.
    if config.get("has_ground_truth", True):
        print("  [FAIL] config.json says has_ground_truth=true for a real dataset")
    # An all-empty genome_label column reads back as all-NaN float, so emptiness
    # is tested on the values rather than on their string form - `astype(str)`
    # of an all-NaN column is not the literal "nan" on every pandas version.
    label = meta["genome_label"]
    has_label = label.notna() & label.astype("string").fillna("").str.strip().ne("")
    if bool(has_label.any()):
        print(f"  [FAIL] contig_metadata.tsv carries genome labels for "
              f"{int(has_label.sum())} contig(s); a real dataset has none")
    if int((data.y >= 0).sum()) != 0:
        print("  [FAIL] viral_graph.pt has labelled nodes for a real dataset")

    print("\n=== contigs ===")
    print(f"  total                    : {n:,}")
    print(f"  >= {COHORT_MIN_BP:,} bp (cohort)      : {cohort.sum():,} "
          f"({100 * cohort.mean():.0f}%)")
    print(f"  total assembled bases    : {lengths.sum():,}")
    print(f"  contig length (bp)       : min={lengths.min():,} "
          f"median={int(np.median(lengths)):,} max={lengths.max():,}")
    order = np.sort(lengths)[::-1]
    n50 = order[np.searchsorted(np.cumsum(order), order.sum() / 2)]
    print(f"  N50                      : {n50:,}")
    for lo, hi, name in [(0, 5000, "<5kb"), (5000, 10000, "5-10kb"),
                         (10000, 10 ** 9, ">10kb")]:
        c = int(((lengths >= lo) & (lengths < hi)).sum())
        print(f"    {name:<7}: {c:,} ({100 * c / max(n, 1):.0f}%)")
    if "viral_identification" in json.loads(
            (OUT_DIR / "manifest.json").read_text()):
        report = json.loads((OUT_DIR / "manifest.json").read_text())["viral_identification"]
        print(f"  viral identification     : {report.get('tool')} "
              f"kept {report.get('kept'):,}/{report.get('input_contigs'):,}")

    # ---- 1. differential coverage -----------------------------------------
    tnf_end = config["coverage_feature_start"]
    cov_end = config["coverage_feature_end"]
    n_samples = config["n_samples"]
    cov = data.x.numpy()[:, tnf_end:cov_end]
    print(f"\n=== differential coverage ({n_samples} sample(s)) ===")
    if cov.shape[1] != n_samples:
        print(f"  [FAIL] coverage block is {cov.shape[1]} wide, config says {n_samples}")
    dead = [s for s in range(cov.shape[1]) if np.allclose(cov[:, s], cov[0, s])]
    if dead:
        print(f"  [WARN] sample(s) with no contig-to-contig variation: {dead}")
    # Spearman across samples: if every sample ranks contigs the same way there
    # is one abundance axis, not several, and the coverage branch is redundant.
    if cov.shape[1] >= 2 and cohort.sum() > 2:
        ranks = pd.DataFrame(cov[cohort]).rank()
        corr = ranks.corr().to_numpy()
        off = corr[np.triu_indices(len(corr), 1)]
        print(f"  between-sample rank corr : median {np.median(off):.3f} "
              f"(min {off.min():.3f}, max {off.max():.3f})")
        if np.median(off) > 0.95:
            print("  [WARN] samples rank contigs almost identically - little "
                  "differential-coverage signal to bin on")
    else:
        print("  [WARN] fewer than 2 samples: no differential coverage at all")

    # ---- 2. assembly graph -------------------------------------------------
    edge_index = data.edge_index.numpy()
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    for a, b in edge_index.T:
        if a != b:
            graph.add_edge(int(a), int(b))
    component = np.full(n, -1)
    for c, members in enumerate(nx.connected_components(graph)):
        for i in members:
            component[i] = c
    degree = np.asarray([graph.degree(i) for i in range(n)])
    multi_component = [len(c) for c in nx.connected_components(graph) if len(c) > 1]
    print("\n=== assembly graph ===")
    print(f"  GFA contig edges         : {graph.number_of_edges():,}")
    print(f"  contigs with >=1 edge    : {int((degree > 0).sum()):,} "
          f"({100 * (degree > 0).mean():.1f}%)")
    print(f"  multi-contig components  : {len(multi_component):,}"
          + (f" (largest {max(multi_component):,} contigs)" if multi_component else ""))

    # ---- 3. shared protein clusters ---------------------------------------
    clusters_path = OUT_DIR / "protein_clusters.tsv"
    shared_partner = np.zeros(n, dtype=bool)
    print("\n=== shared protein clusters ===")
    if clusters_path.is_file():
        clusters = pd.read_csv(clusters_path, sep="\t")
        index = {c: i for i, c in enumerate(meta.contig_id.astype(str))}
        by_cluster = clusters.groupby("cluster_id")["contig_id"].apply(
            lambda s: sorted(set(s.astype(str))))
        multi = [members for members in by_cluster if len(members) > 1]
        for members in multi:
            for c in members:
                if c in index:
                    shared_partner[index[c]] = True
        print(f"  protein clusters         : {len(by_cluster):,}")
        print(f"  multi-contig clusters    : {len(multi):,}")
        print(f"  contigs sharing a cluster: {int(shared_partner.sum()):,} "
              f"({100 * shared_partner.mean():.1f}%)")
    else:
        print(f"  [WARN] {clusters_path.name} not found - phase 2's viral-feature "
              "stage did not run")

    # ---- verdict -----------------------------------------------------------
    # "Joinable" = this contig has at least one reference-free reason to share a
    # bin with another: an assembly link, or a protein cluster in common. It is
    # a LOWER BOUND on the merges available, never an expected bin count.
    joinable = (degree > 0) | shared_partner
    cohort_joinable = joinable & cohort
    frac = cohort_joinable.sum() / max(cohort.sum(), 1)
    print("\n=== verdict ===")
    print(f"  cohort contigs (>= {COHORT_MIN_BP:,} bp) : {int(cohort.sum()):,}")
    print(f"  with a reference-free partner : {int(cohort_joinable.sum()):,} "
          f"({100 * frac:.1f}%)")
    verdict = ("STRONG - plenty to merge" if frac >= STRONG_JOINABLE_FRAC
               else ("usable" if frac >= USABLE_JOINABLE_FRAC
                     else "WEAK - little evidence any contigs co-bin; expect "
                          "near-singleton output from every method"))
    print("VERDICT:", verdict)
    print("\nNote: these are lower bounds from assembly and protein evidence, not "
          "truth.\n      Bin quality is measured after binning by "
          "4_checkv_evaluate.py.")


if __name__ == "__main__":
    main()
