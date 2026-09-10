#!/usr/bin/env python

# METRICS REPORTED
# ----------------
#   Predicted vMAGs      bins evaluated (and, separately, multi-contig bins -
#                        a single-contig "bin" is the assembler's result, not
#                        the binner's, so both denominators are printed)
#   CheckV Complete      checkv_quality == "Complete"      (100% + DTR/ITR etc.)
#   CheckV HQ            checkv_quality == "High-quality"  (>90% complete)
#   CheckV MQ            checkv_quality == "Medium-quality"(50-90% complete)
#   MQ or better         Complete + HQ + MQ
#   Mean completeness    mean of CheckV's completeness estimate
#   Median completeness  median of the same
#   plus: low-quality / not-determined counts, mean+median contamination,
#         bins with contamination > 5% (the HQ MIUViG ceiling), and the
#         contig/length distribution of the bins - so a method that scores well
#         by emitting one enormous bin is visible rather than hidden.

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

def _load_env_sh() -> None:
    env_sh = Path(__file__).resolve().parent.parent / "env.sh"
    if not env_sh.exists():
        return
    for line in env_sh.read_text().splitlines():
        line = line.strip()
        if line.startswith("export "):
            line = line[len("export "):]
        key, sep, value = line.partition("=")
        key = key.strip()
        if not sep or key.startswith("#"):
            continue
        if not (key.startswith("VB_") or key in ("CHECKVDB", "GENOMAD_DB")):
            continue
        value = value.split("#")[0].strip().strip('"').strip("'")
        # Only the `${VAR:-default}` form the template uses is expanded; any
        # other shell syntax is left alone rather than half-interpreted.
        if value.startswith("${") and value.endswith("}") and ":-" in value:
            value = value[2:-1].split(":-", 1)[1]
        value = os.path.expandvars(os.path.expanduser(value))
        if value:
            os.environ.setdefault(key, value)

_load_env_sh()

# ---- config (kept local so this script runs standalone) ----
HERE = Path(__file__).resolve().parent
REPO_HOME = HERE.parent.parent
REGISTRY = HERE / "datasets.tsv"
CONDA_BIN = os.environ.get(
    "VB_CONDA_BIN", os.path.expanduser("~/miniforge3/envs/viralbin/bin"))
CHECKV_DB = os.environ.get(
    "CHECKVDB", os.path.expanduser("~/checkv-db/checkv-db-v1.5"))
THREADS = int(os.environ.get("VB_THREADS", 8))

os.environ["PATH"] = os.pathsep.join([CONDA_BIN, os.environ.get("PATH", "")])
CHECKV_EXTERNALS = ("hmmsearch", "diamond")

# vRhyme's default linker: 1,500 Ns between scaffolds of one bin.
LINKER_N = 1500
LINKER_CHAR = "N"
# CheckV's own quality tiers; no thresholds are re-derived here.
QUALITY_TIERS = ["Complete", "High-quality", "Medium-quality",
                 "Low-quality", "Not-determined"]


def resolve(name: str) -> str:
    candidate = os.path.join(CONDA_BIN, name)
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    found = shutil.which(name)
    if found:
        return found
    raise SystemExit(f"error: {name} not found (install: conda install -n "
                     f"viralbin -c bioconda {name})")


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


def read_fasta(path: Path) -> dict[str, str]:
    """{contig_id: sequence}; the id is the header up to the first whitespace."""
    sequences: dict[str, str] = {}
    name, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences

# Parse a bin assignment table -> {bin_id: [contig_id, ...]}
def read_bins(path: Path) -> dict[str, list[str]]:
    with open(path, encoding="utf-8") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh
                if line.strip() and not line.startswith("#")]
    if not rows:
        raise SystemExit(f"error: {path} is empty")
    if len(rows[0]) < 2:
        raise SystemExit(f"error: {path} is not a TSV with >=2 columns "
                         f"(first line: {rows[0]!r})")

    CONTIG_NAMES = {"contig", "contig_id", "contigname", "contig_name",
                    "sequence", "sequence_id", "scaffold", "scaffold_id",
                    "name", "seq_name"}
    BIN_NAMES = {"bin", "bin_id", "cluster", "cluster_id", "clustername",
                 "cluster_name", "binid", "group"}
    header = [h.strip().lower() for h in rows[0]]
    named_contig = next((i for i, h in enumerate(header) if h in CONTIG_NAMES), None)
    named_bin = next((i for i, h in enumerate(header) if h in BIN_NAMES), None)
    if named_contig is not None and named_bin is not None:
        contig_col, bin_col, body = named_contig, named_bin, rows[1:]
    else:
        body = rows[1:] if (named_contig is not None or named_bin is not None) else rows
        left = {r[0] for r in body if len(r) > 1}
        right = {r[1] for r in body if len(r) > 1}
        contig_col, bin_col = (0, 1) if len(left) >= len(right) else (1, 0)

    bins: dict[str, list[str]] = defaultdict(list)
    for row in body:
        if len(row) <= max(contig_col, bin_col):
            continue
        bins[row[bin_col].strip()].append(row[contig_col].strip())
    if not bins:
        raise SystemExit(f"error: no assignments parsed from {path}")
    return dict(bins)


def link_bins(bins: dict[str, list[str]], sequences: dict[str, str],
              out_fasta: Path) -> tuple[dict[str, dict], list[str]]:
    """Write one linked pseudo-genome per bin (vRhyme's linking recipe).

    Returns per-bin bookkeeping and the list of contigs that were assigned but
    are absent from the contig FASTA (reported, never silently dropped).
    """
    linker = LINKER_CHAR * LINKER_N
    info: dict[str, dict] = {}
    missing: list[str] = []
    with open(out_fasta, "w", encoding="utf-8", newline="\n") as fh:
        for bin_id in sorted(bins):
            members, parts, total = [], [], 0
            for contig in bins[bin_id]:
                seq = sequences.get(contig)
                if seq is None:
                    missing.append(contig)
                    continue
                members.append(contig)
                parts.append(seq)
                total += len(seq)
            if not parts:
                continue
            # Bin ids come from arbitrary tools; make them FASTA/CheckV safe and
            # keep the mapping so results can be traced back.
            safe = "bin_" + "".join(
                c if c.isalnum() or c in "._-" else "_" for c in str(bin_id))
            info[safe] = {"bin_id": bin_id, "contigs": members,
                          "n_contigs": len(members), "bases": total}
            fh.write(f">{safe}\n")
            linked = linker.join(parts)
            for i in range(0, len(linked), 60):
                fh.write(linked[i:i + 60] + "\n")
    return info, missing


def run_checkv(fasta: Path, outdir: Path, *, force: bool, threads: int = THREADS) -> Path:
    summary = outdir / "quality_summary.tsv"
    if summary.is_file() and not force:
        print(f"[checkv] reusing {summary} (--force to rerun)")
        return summary
    if not os.path.isdir(CHECKV_DB):
        raise SystemExit(
            f"error: CheckV database not found at {CHECKV_DB}\n"
            "  install: checkv download_database ~/checkv-db\n"
            "  or set CHECKVDB in Dataset_Processing/env.sh")
    # Resolved up front: CheckV reports a missing hmmsearch only after it has
    # already called every gene, which on a large dataset is a long wait for a
    # failure that was knowable at the start.
    for tool in CHECKV_EXTERNALS:
        resolve(tool)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [resolve("checkv"), "end_to_end", str(fasta), str(outdir),
           "-d", CHECKV_DB, "-t", str(threads), "--remove_tmp"]
    print("[checkv] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    if not summary.is_file():
        raise SystemExit(f"error: CheckV produced no {summary}")
    return summary


def summarise(summary_tsv: Path, info: dict[str, dict], method: str,
              dataset: str) -> dict:
    with open(summary_tsv, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    def number(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    quality = {tier: 0 for tier in QUALITY_TIERS}
    completeness, contamination = [], []
    per_bin = []
    for row in rows:
        tier = row.get("checkv_quality", "Not-determined")
        quality[tier] = quality.get(tier, 0) + 1
        comp = number(row.get("completeness"))
        cont = number(row.get("contamination"))
        if comp is not None:
            completeness.append(comp)
        if cont is not None:
            contamination.append(cont)
        meta = info.get(row["contig_id"], {})
        per_bin.append({
            "bin": meta.get("bin_id", row["contig_id"]),
            "linked_id": row["contig_id"],
            "n_contigs": meta.get("n_contigs"),
            "bases": meta.get("bases"),
            "checkv_quality": tier,
            "completeness": comp,
            "contamination": cont,
            "miuvig_quality": row.get("miuvig_quality"),
            "gene_count": row.get("gene_count"),
            "viral_genes": row.get("viral_genes"),
            "host_genes": row.get("host_genes"),
            "warnings": row.get("warnings", ""),
        })

    def mean(values):
        return round(sum(values) / len(values), 2) if values else None

    def median(values):
        if not values:
            return None
        ordered = sorted(values)
        mid = len(ordered) // 2
        return round(ordered[mid] if len(ordered) % 2
                     else (ordered[mid - 1] + ordered[mid]) / 2, 2)

    sizes = [b["n_contigs"] for b in per_bin if b["n_contigs"]]
    multi = [b for b in per_bin if (b["n_contigs"] or 0) > 1]
    mq_or_better = (quality["Complete"] + quality["High-quality"]
                    + quality["Medium-quality"])
    return {
        "dataset": dataset,
        "method": method,
        "predicted_vmags": len(per_bin),
        "multi_contig_vmags": len(multi),
        "contigs_binned": sum(sizes),
        "checkv_complete": quality["Complete"],
        "checkv_high_quality": quality["High-quality"],
        "checkv_medium_quality": quality["Medium-quality"],
        "mq_or_better": mq_or_better,
        "checkv_low_quality": quality["Low-quality"],
        "checkv_not_determined": quality["Not-determined"],
        "mean_completeness": mean(completeness),
        "median_completeness": median(completeness),
        "completeness_n": len(completeness),
        "mean_contamination": mean(contamination),
        "median_contamination": median(contamination),
        "contaminated_over_5pct": sum(1 for c in contamination if c > 5.0),
        "largest_bin_contigs": max(sizes) if sizes else 0,
        "median_bin_contigs": median([float(s) for s in sizes]),
        "checkv_db": CHECKV_DB,
        "linker_n": LINKER_N,
        "_per_bin": per_bin,
    }


def print_report(result: dict) -> None:
    print("\n" + "=" * 68)
    print(f"CheckV evaluation  ::  {result['method']}  on  {result['dataset']}")
    print("=" * 68)
    total = result["predicted_vmags"] or 1
    rows = [
        ("Predicted vMAGs", result["predicted_vmags"], ""),
        ("  multi-contig vMAGs", result["multi_contig_vmags"],
         f"{100 * result['multi_contig_vmags'] / total:.0f}% of bins"),
        ("  contigs binned", result["contigs_binned"], ""),
        ("CheckV Complete", result["checkv_complete"],
         f"{100 * result['checkv_complete'] / total:.1f}%"),
        ("CheckV HQ", result["checkv_high_quality"],
         f"{100 * result['checkv_high_quality'] / total:.1f}%"),
        ("CheckV MQ", result["checkv_medium_quality"],
         f"{100 * result['checkv_medium_quality'] / total:.1f}%"),
        ("MQ or better", result["mq_or_better"],
         f"{100 * result['mq_or_better'] / total:.1f}%"),
        ("CheckV LQ", result["checkv_low_quality"], ""),
        ("CheckV not-determined", result["checkv_not_determined"], ""),
    ]
    for label, value, note in rows:
        print(f"  {label:<24}{value:>8}   {note}")
    print(f"  {'Mean completeness':<24}{str(result['mean_completeness']):>8}   "
          f"over {result['completeness_n']} bin(s) CheckV could place")
    print(f"  {'Median completeness':<24}{str(result['median_completeness']):>8}")
    print(f"  {'Mean contamination':<24}{str(result['mean_contamination']):>8}")
    print(f"  {'Median contamination':<24}{str(result['median_contamination']):>8}")
    print(f"  {'Bins >5% contaminated':<24}{result['contaminated_over_5pct']:>8}")
    print(f"  {'Largest bin (contigs)':<24}{result['largest_bin_contigs']:>8}")
    print(f"  {'Median bin (contigs)':<24}{str(result['median_bin_contigs']):>8}")
    if result["checkv_not_determined"]:
        print(f"\n  note: {result['checkv_not_determined']} bin(s) are "
              "Not-determined and are excluded from the completeness averages.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="CheckV evaluation of viral bins (vRhyme-style linking).")
    parser.add_argument("--bins", nargs="+", required=True, type=Path,
                        help="Bin assignment TSV(s): contig_id + bin_id.")
    parser.add_argument("--method", nargs="+",
                        help="Method name per --bins file (default: file stem).")
    parser.add_argument("--dataset",
                        help="Registry dataset name; supplies the contig FASTA "
                             "and names the output. Omit with --contigs.")
    parser.add_argument("--contigs", type=Path,
                        help="Contig FASTA, if not taken from --dataset.")
    parser.add_argument("--outdir", type=Path,
                        help="Where to write results "
                             "(default Outputs/checkv/<dataset>).")
    parser.add_argument("--force", action="store_true",
                        help="Rerun CheckV even if a summary already exists.")
    parser.add_argument("--max-bin-contigs", type=int,
                        help="Exclude bins larger than this from CheckV. A bin of "
                             "thousands of contigs is not a genome - a phage "
                             "genome is 10-200 kb - so its completeness estimate "
                             "is meaningless, and linking it produces a single "
                             "sequence tens of Mbp long that exhausts memory in "
                             "hmmsearch and takes the whole evaluation down with "
                             "it. Excluding it keeps the numbers for every other "
                             "bin; what was dropped is printed and recorded.")
    parser.add_argument("--threads", type=int, default=THREADS,
                        help=f"Threads for CheckV (default {THREADS}). CheckV runs "
                             "this many hmmsearch processes at once, each holding "
                             "the query proteins, so on a large or badly-binned "
                             "input the default can exhaust memory - the symptom "
                             "is 'N hmmsearch tasks failed'. Lower it to trade "
                             "time for memory.")
    args = parser.parse_args()

    methods = args.method or [p.stem for p in args.bins]
    if len(methods) != len(args.bins):
        raise SystemExit(f"error: --method has {len(methods)} name(s) for "
                         f"{len(args.bins)} --bins file(s)")

    # Contigs come either from a registry dataset or from an explicit path, so
    # this works on a real dataset and on a simulated one (where it is a
    # reference-free second opinion alongside the truth-based metrics).
    if args.contigs:
        contigs_path = args.contigs
        dataset = args.dataset or contigs_path.stem
    else:
        if not args.dataset:
            raise SystemExit("error: pass --dataset (registry name) or --contigs")
        registry = load_registry(REGISTRY)
        if args.dataset not in registry:
            raise SystemExit(
                f"error: dataset {args.dataset!r} is not in {REGISTRY}; "
                f"known: {', '.join(sorted(registry))}\n"
                "       (for a dataset outside the registry, pass --contigs)")
        row = registry[args.dataset]
        assembler = {"short": "metaspades", "long": "metaflye"}[row["read_type"]]
        tag = "noviral" if row["viral_enriched"] == "yes" else "genomad"
        run_name = f"{args.dataset}__coasm__{assembler}__{tag}"
        config_path = REPO_HOME / f"Data/Processed_data/{run_name}/config.json"
        if not config_path.is_file():
            raise SystemExit(f"error: {config_path} not found - run phase 2 first")
        contigs_path = Path(json.loads(config_path.read_text())["combined_fasta"])
        dataset = run_name
    if not contigs_path.is_file():
        raise SystemExit(f"error: contig FASTA not found: {contigs_path}")

    outdir = args.outdir or (REPO_HOME / "Outputs" / "checkv" / dataset)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"contigs : {contigs_path}")
    print(f"outputs : {outdir}")
    sequences = read_fasta(contigs_path)
    print(f"          {len(sequences):,} contigs loaded")

    results = []
    for bins_path, method in zip(args.bins, methods):
        if not bins_path.is_file():
            raise SystemExit(f"error: bins file not found: {bins_path}")
        print(f"\n--- {method} ({bins_path}) ---")
        bins = read_bins(bins_path)
        excluded = {}
        if args.max_bin_contigs:
            excluded = {bin_id: members for bin_id, members in bins.items()
                        if len(members) > args.max_bin_contigs}
            bins = {bin_id: members for bin_id, members in bins.items()
                    if bin_id not in excluded}
            if excluded:
                dropped_contigs = sum(len(m) for m in excluded.values())
                print(f"[link] EXCLUDED {len(excluded)} bin(s) over "
                      f"{args.max_bin_contigs} contigs "
                      f"({dropped_contigs:,} contigs): "
                      + ", ".join(f"{bin_id}({len(members):,})"
                                  for bin_id, members in
                                  sorted(excluded.items(),
                                         key=lambda kv: -len(kv[1]))[:5]))
                print("       these are not scored; a bin this size has no "
                      "meaningful completeness")
            if not bins:
                raise SystemExit(
                    f"error: every bin exceeds --max-bin-contigs "
                    f"{args.max_bin_contigs}; nothing left to score")
        work = outdir / method
        work.mkdir(parents=True, exist_ok=True)
        linked_fasta = work / "linked_bins.fna"
        info, missing = link_bins(bins, sequences, linked_fasta)
        print(f"[link] {len(bins):,} bin(s) -> {len(info):,} linked sequence(s) "
              f"({LINKER_N} {LINKER_CHAR}s between contigs)")
        if missing:
            print(f"[link] WARNING: {len(missing)} assigned contig(s) are not in "
                  f"the FASTA and were dropped, e.g. {sorted(set(missing))[:3]}")
        if not info:
            raise SystemExit(f"error: no bin had any contig present in "
                             f"{contigs_path}")

        summary_tsv = run_checkv(linked_fasta, work / "checkv", force=args.force,
                                 threads=args.threads)
        result = summarise(summary_tsv, info, method, dataset)
        result["bins_excluded_oversized"] = len(excluded)
        result["contigs_excluded_oversized"] = sum(len(m) for m in excluded.values())
        print_report(result)

        per_bin = result.pop("_per_bin")
        columns = list(per_bin[0])
        with open(work / "per_bin.tsv", "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(columns) + "\n")
            for entry in sorted(per_bin, key=lambda e: -(e["completeness"] or -1)):
                fh.write("\t".join("" if entry[c] is None else str(entry[c])
                                   for c in columns) + "\n")
        print(f"  per-bin table -> {work / 'per_bin.tsv'}")
        results.append(result)

    # One row per method, so several methods on one dataset compare directly.
    columns = [c for c in results[0] if not c.startswith("_")]
    out_csv = outdir / "checkv_metrics.csv"
    existing = []
    if out_csv.is_file():
        with open(out_csv, encoding="utf-8") as fh:
            existing = [r for r in csv.DictReader(fh)
                        if r.get("method") not in methods]
    with open(out_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in existing:
            writer.writerow({c: row.get(c, "") for c in columns})
        for row in results:
            writer.writerow({c: row[c] for c in columns})
    print(f"\nwrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
