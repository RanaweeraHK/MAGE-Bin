#!/usr/bin/env python
"""Phase 3 - label each contig with the reference genome it came from.

    VB_DATASET=phage_mock_illumina python 3_label_from_references.py
    VB_DATASET=dsmz_dsrna python 3_label_from_references.py --min-identity 0.92

A mock community was built from known genomes, so the truth the simulated
pipeline gets from a simulator, and the real pipeline cannot get at all, can be
recovered here by aligning contigs back to the reference set. Once
`genome_label` is filled, precision, recall, F1, ARI and NMI are defined on
real sequencing data.

    contigs_filt.fasta --minimap2 -x asm20--> reference_genomes.fasta
        -> per (contig, reference): matching bases, covered query span
        -> best reference, subject to --min-identity and --min-query-cov
        -> contig_metadata.tsv, truth_per_reference.tsv, truth_summary.json

Three rules decide what a later score means:

  Ambiguity. Genomes that no assembly can separate (this phage mock contains
  J1 and J2 at 100% identity) are merged into one truth genome, `J1|J2`,
  across the whole dataset. Picking one of them would invent a truth and
  penalise a binner for not reproducing it.

  Chimeras. A contig covered by two genomes on different parts has no single
  true genome, so it is left unlabelled and counted. Scoring it would measure
  the assembler.

  Multipartite genomes. `genome_label` holds the AGENT, so a binner is asked
  to co-bin a segmented virus; `molecule_label` keeps the segment beside it.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_HOME = HERE.parent.parent
REGISTRY = HERE / "datasets.tsv"

# Contigs come from the isolates in the reference set, so real hits are
# near-identical; these thresholds reject a contig that merely shares a gene.
MIN_IDENTITY = 0.90        # matching bases / alignment block length
MIN_QUERY_COV = 0.50       # fraction of the contig covered by that reference
AMBIGUITY_MARGIN = 0.02    # within 2% of the best score counts as a tie
CHIMERA_MIN_COV = 0.25     # rival genome covering this much of the rest of the contig


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _load_env_sh() -> None:
    env_sh = HERE.parent / "env.sh"
    if not env_sh.exists():
        return
    for line in env_sh.read_text().splitlines():
        line = line.strip()
        if not line.startswith("export "):
            continue
        key, _, value = line[len("export "):].partition("=")
        key, value = key.strip(), value.strip().strip('"').strip("'")
        if "${" in value or not key:
            # env.sh writes VAR="${VAR:-default}"; take the default.
            value = value.split(":-", 1)[-1].rstrip("}") if ":-" in value else value
        os.environ.setdefault(key, value)


_load_env_sh()


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


def resolve_minimap2() -> str:
    """Find minimap2 in asm_env, where the coverage stage finds it too."""
    for candidate in (Path(os.environ.get("VB_ASM_CONDA_BIN", "")) / "minimap2",
                      Path(os.environ.get("VB_CONDA_BIN", "")) / "minimap2"):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    found = shutil.which("minimap2")
    if found:
        return found
    raise SystemExit(
        "error: minimap2 not found.\n"
        "       Set VB_ASM_CONDA_BIN in Dataset_Processing/env.sh, or:\n"
        "       conda create -n asm_env -c bioconda -c conda-forge minimap2")


# ===========================================================================
# Alignment
# ===========================================================================
def run_minimap2(contigs: Path, references: Path, paf: Path, *,
                 threads: int, force: bool) -> Path:
    if paf.is_file() and paf.stat().st_size > 0 and not force:
        log(f"alignment already present: {paf.name} (skipping minimap2)")
        return paf
    binary = resolve_minimap2()
    # asm20 allows ~20% divergence between a reference and the isolate used.
    # -N/-p keep secondary hits, without which ties would be invisible.
    cmd = [binary, "-cx", "asm20", "-N", "50", "-p", "0.1",
           "-t", str(threads), str(references), str(contigs)]
    log("minimap2 " + " ".join(cmd[1:]))
    tmp = paf.with_suffix(".paf.part")
    with open(tmp, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    tmp.replace(paf)
    return paf


def merge_intervals(spans: list[tuple[int, int]]) -> int:
    """Total length covered by ``spans``, counting overlap once."""
    if not spans:
        return 0
    spans = sorted(spans)
    total, cur_start, cur_end = 0, *spans[0]
    for start, end in spans[1:]:
        if start > cur_end:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    return total + (cur_end - cur_start)


def parse_paf(paf: Path) -> tuple[dict, dict]:
    """contig -> reference -> stats, plus each reference's covered spans."""
    hits: dict[str, dict[str, dict]] = defaultdict(lambda: defaultdict(
        lambda: {"matches": 0, "block": 0, "qspans": [], "tspans": []}))
    ref_spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    ref_len: dict[str, int] = {}
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            qname, qstart, qend = f[0], int(f[2]), int(f[3])
            tname, tlen, tstart, tend = f[5], int(f[6]), int(f[7]), int(f[8])
            matches, block = int(f[9]), int(f[10])
            entry = hits[qname][tname]
            entry["matches"] += matches
            entry["block"] += block
            entry["qspans"].append((qstart, qend))
            entry["tspans"].append((tstart, tend))
            ref_spans[tname].append((tstart, tend))
            ref_len[tname] = tlen
    return hits, {"spans": ref_spans, "lengths": ref_len}


# ===========================================================================
# Assignment
# ===========================================================================
class UnionFind:
    """Merges genomes that no assembly in this dataset can tell apart."""

    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)

    def groups(self) -> dict[str, list[str]]:
        out: dict[str, list[str]] = defaultdict(list)
        for item in self.parent:
            out[self.find(item)].append(item)
        return out


def score_contig(per_ref: dict, ref_to_agent: dict, contig_len: int,
                 *, min_identity: float, min_query_cov: float) -> list[dict]:
    """One row per reference that clears both thresholds, best first."""
    scored = []
    for ref, stats in per_ref.items():
        if stats["block"] == 0:
            continue
        identity = stats["matches"] / stats["block"]
        covered = merge_intervals(stats["qspans"])
        query_cov = covered / contig_len if contig_len else 0.0
        if identity < min_identity or query_cov < min_query_cov:
            continue
        scored.append({"reference": ref, "agent": ref_to_agent.get(ref, ref),
                       "identity": identity, "query_cov": query_cov,
                       "covered": covered, "matches": stats["matches"],
                       "qspans": stats["qspans"]})
    scored.sort(key=lambda r: (-r["matches"], r["reference"]))
    return scored


def disjoint_fraction(spans_a: list[tuple[int, int]],
                      spans_b: list[tuple[int, int]], contig_len: int) -> float:
    """Fraction of the contig covered by B but not A - the chimera signal."""
    if not spans_b or not contig_len:
        return 0.0
    covered_a = set()
    for start, end in spans_a:
        covered_a.update(range(start, end))
    extra = set()
    for start, end in spans_b:
        extra.update(range(start, end))
    return len(extra - covered_a) / contig_len


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--min-identity", type=float, default=MIN_IDENTITY,
                        help=f"matching bases / block length (default {MIN_IDENTITY})")
    parser.add_argument("--min-query-cov", type=float, default=MIN_QUERY_COV,
                        help=f"fraction of contig covered (default {MIN_QUERY_COV})")
    parser.add_argument("--ambiguity-margin", type=float, default=AMBIGUITY_MARGIN,
                        help=f"tie window below the best hit (default {AMBIGUITY_MARGIN})")
    parser.add_argument("--level", choices=("agent", "molecule"), default="agent",
                        help="what genome_label holds (default agent)")
    parser.add_argument("--force", action="store_true",
                        help="re-run minimap2 even if the PAF is present")
    parser.add_argument("--threads", type=int,
                        default=int(os.environ.get("VB_THREADS", "8")))
    args = parser.parse_args()

    dataset = os.environ.get("VB_DATASET", "")
    registry = load_registry(REGISTRY)
    if dataset not in registry:
        raise SystemExit(
            f"error: VB_DATASET must name a row in {REGISTRY}\n"
            f"       got {dataset!r}; known: {', '.join(sorted(registry))}")
    row = registry[dataset]

    # Same run naming as phase 2 and run_all.sh; keep the three in step.
    declared = row.get("assembler", "-")
    assembler = (declared if declared not in ("-", "")
                 else {"short": "metaspades", "long": "metaflye"}[row["read_type"]])
    tag = "noviral" if row["viral_enriched"] == "yes" else "genomad"
    run_name = f"{dataset}__coasm__{assembler}__{tag}"

    mwork = REPO_HOME / "Data/work" / f"mock_{dataset}" / "standardize_work" / run_name
    mout = REPO_HOME / "Data/Processed_data" / run_name
    refdir = (REPO_HOME / "Data/Mock_dataset" / row["study"] / "references"
              / row["reference_set"])
    contigs = mwork / ("contigs_filt.fasta" if row["viral_enriched"] == "yes"
                       else "contigs_viral.fasta")
    references = refdir / "reference_genomes.fasta"
    refmap = refdir / "reference_map.tsv"
    metadata = mout / "contig_metadata.tsv"

    for path, hint in ((contigs, "run phase 2 first"),
                       (metadata, "run phase 2 first"),
                       (references, "run phase 1 first"),
                       (refmap, "run phase 1 first")):
        if not path.is_file():
            raise SystemExit(f"error: {path} not found - {hint}\n"
                             f"       VB_DATASET={dataset} ...")

    ref_to_agent, ref_to_molecule = {}, {}
    with open(refmap) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            ref_to_agent[rec["reference_id"]] = rec["agent_label"]
            ref_to_molecule[rec["reference_id"]] = rec["molecule_label"]
    log(f"reference set: {len(ref_to_agent)} molecules, "
        f"{len(set(ref_to_agent.values()))} agents")

    paf = run_minimap2(contigs, references, mwork / "contigs_vs_references.paf",
                       threads=args.threads, force=args.force)
    hits, ref_cov = parse_paf(paf)
    log(f"alignments: {len(hits):,} contigs have at least one hit")

    rows = list(csv.DictReader(open(metadata), delimiter="\t"))
    lengths = {r["contig_id"]: int(r["length"]) for r in rows}

    # ---- pass 1: score every contig, and learn which genomes are ties ----
    level_of = ref_to_agent if args.level == "agent" else ref_to_molecule
    scored_by_contig: dict[str, list[dict]] = {}
    ties = UnionFind()
    for contig, per_ref in hits.items():
        if contig not in lengths:
            continue                       # filtered out by phase 2
        scored = score_contig(per_ref, level_of, lengths[contig],
                              min_identity=args.min_identity,
                              min_query_cov=args.min_query_cov)
        if not scored:
            continue
        scored_by_contig[contig] = scored
        best = scored[0]["matches"]
        tied = {r["agent"] for r in scored
                if r["matches"] >= best * (1.0 - args.ambiguity_margin)}
        for other in sorted(tied)[1:]:
            ties.union(sorted(tied)[0], other)

    groups = ties.groups()
    group_label = {}
    for members in groups.values():
        if len(members) > 1:
            label = "|".join(sorted(members))
            for member in members:
                group_label[member] = label
    merged = sorted({v for v in group_label.values()})
    if merged:
        log(f"indistinguishable genomes merged into {len(merged)} truth "
            f"genome(s): {', '.join(merged)}")

    # ---- pass 2: assign ----
    labels, molecules, status, identities, covs = {}, {}, {}, {}, {}
    counts = defaultdict(int)
    for contig, scored in scored_by_contig.items():
        best = scored[0]
        truth = group_label.get(best["agent"], best["agent"])
        rival = next((r for r in scored
                      if group_label.get(r["agent"], r["agent"]) != truth
                      and disjoint_fraction(best["qspans"], r["qspans"],
                                            lengths[contig]) >= CHIMERA_MIN_COV),
                     None)
        if rival is not None:
            # Two genomes on different parts of one contig: no single truth.
            status[contig] = "chimeric"
            counts["chimeric"] += 1
            continue
        labels[contig] = truth
        molecules[contig] = ref_to_molecule.get(best["reference"], best["reference"])
        identities[contig] = best["identity"]
        covs[contig] = best["query_cov"]
        status[contig] = "merged" if "|" in truth else "labelled"
        counts[status[contig]] += 1

    for rec in rows:
        counts["total"] += 1
        if rec["contig_id"] not in status:
            counts["unlabelled"] += 1

    # ---- rewrite contig_metadata.tsv ----
    fieldnames = ["node_index", "contig_id", "source_sample", "original_contig_id",
                  "length", "genome_label", "molecule_label", "label_identity",
                  "label_query_cov", "label_status"]
    tmp = metadata.with_suffix(".tsv.part")
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for rec in rows:
            cid = rec["contig_id"]
            rec["genome_label"] = labels.get(cid, "")
            rec["molecule_label"] = molecules.get(cid, "")
            rec["label_identity"] = f"{identities[cid]:.4f}" if cid in identities else ""
            rec["label_query_cov"] = f"{covs[cid]:.4f}" if cid in covs else ""
            rec["label_status"] = status.get(cid, "unlabelled")
            writer.writerow(rec)
    tmp.replace(metadata)
    log(f"wrote {metadata}")

    # ---- per-reference recovery: was this genome even assembled? ----
    per_ref_out = mout / "truth_per_reference.tsv"
    contigs_per_label = defaultdict(int)
    bp_per_label = defaultdict(int)
    for cid, label in labels.items():
        contigs_per_label[label] += 1
        bp_per_label[label] += lengths[cid]
    with open(per_ref_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["reference_id", "molecule_label", "agent_label",
                         "truth_genome", "reference_bp", "reference_bp_covered",
                         "reference_fraction_covered", "contigs_assigned",
                         "contig_bp_assigned"])
        for ref in sorted(ref_to_agent):
            agent = ref_to_agent[ref]
            truth = group_label.get(agent, agent)
            rlen = ref_cov["lengths"].get(ref, 0)
            covered = merge_intervals(ref_cov["spans"].get(ref, []))
            writer.writerow([
                ref, ref_to_molecule[ref], agent, truth, rlen, covered,
                f"{covered / rlen:.4f}" if rlen else "0.0000",
                contigs_per_label.get(truth, 0), bp_per_label.get(truth, 0)])
    log(f"wrote {per_ref_out}")

    recovered = sum(1 for ref in ref_to_agent
                    if merge_intervals(ref_cov["spans"].get(ref, []))
                    >= 0.5 * ref_cov["lengths"].get(ref, 1))
    summary = {
        "run": run_name, "dataset": dataset, "level": args.level,
        "reference_set": row["reference_set"],
        "references": len(ref_to_agent),
        "agents": len(set(ref_to_agent.values())),
        "truth_genomes": len(set(labels.values())),
        "merged_truth_genomes": merged,
        "references_half_covered": recovered,
        "contigs_total": counts["total"],
        "contigs_labelled": counts["labelled"] + counts["merged"],
        "contigs_labelled_merged": counts["merged"],
        "contigs_chimeric": counts["chimeric"],
        "contigs_unlabelled": counts["unlabelled"],
        "thresholds": {"min_identity": args.min_identity,
                       "min_query_cov": args.min_query_cov,
                       "ambiguity_margin": args.ambiguity_margin,
                       "chimera_min_cov": CHIMERA_MIN_COV},
    }
    (mout / "truth_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    manifest_path = mout / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        manifest["has_ground_truth"] = True
        manifest["labelled"] = True
        manifest["evaluation"] = "reference_labels+checkv"
        manifest["truth"] = summary
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        log(f"updated {manifest_path} (has_ground_truth=True)")

    labelled = counts["labelled"] + counts["merged"]
    print()
    print("=" * 70)
    print(f"  truth for {run_name}")
    print("=" * 70)
    print(f"  contigs              {counts['total']:,}")
    print(f"  labelled             {labelled:,} "
          f"({100.0 * labelled / max(counts['total'], 1):.1f}%)")
    print(f"    of which merged    {counts['merged']:,} (indistinguishable genomes)")
    print(f"  chimeric             {counts['chimeric']:,}")
    print(f"  unlabelled           {counts['unlabelled']:,} "
          f"(no reference hit clearing the thresholds)")
    print(f"  truth genomes        {len(set(labels.values())):,} "
          f"of {len(set(ref_to_agent.values())):,} in the community")
    print(f"  references >=50% covered by the assembly   "
          f"{recovered:,}/{len(ref_to_agent):,}")
    print()
    print("  Unlabelled contigs are material the reference set does not")
    print("  explain, not binner errors: score on the labelled subset.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
