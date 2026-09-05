#!/usr/bin/env python3
"""Run the external-tool benchmark on REAL datasets and score it with CheckV.

The sibling of ``benchmark.py``.  Same tools, same wrappers, same >=2 kb cohort
rule - one difference that changes everything downstream:

    THERE IS NO GROUND TRUTH.

``benchmark.py`` scores every tool against the strain that generated each
contig.  A real dataset has no such reference, so they are undefined.  Bin quality is instead
measured after binning by CheckV, exactly as
``Dataset_Processing/Real_Data/4_checkv_evaluate.py`` does it: each bin's
contigs are concatenated into ONE sequence separated by 1,500 Ns (vRhyme's
linking recipe) so CheckV sees a bin as a single pseudo-genome, then CheckV
estimates its completeness and contamination against its reference database.

    Why the linking matters: CheckV scores one sequence at a time.  Feed it a
    bin as separate contigs and a genome correctly recovered as 5 contigs is
    reported as 5 fragments at ~20% completeness each - so a binner that
    reassembled it perfectly is indistinguishable from one that did nothing.

"""

from __future__ import annotations

import argparse
import csv
import fcntl
import importlib.util
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def _load_benchmark():
    """Import the sibling benchmark module for the shared tool-running code."""
    spec = importlib.util.spec_from_file_location(
        "vb_benchmark", HERE / "benchmark.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = _load_benchmark()

DEFAULT_DATA_ROOT = ROOT / "Data" / "Processed_data"
DEFAULT_OUTPUT = ROOT / "Outputs" / "comparison_real.csv"
DEFAULT_RUN_ROOT = ROOT / "Outputs" / "other_tools_real"
DEFAULT_CHECKV_ROOT = ROOT / "Outputs" / "checkv"
CHECKV_EVALUATOR = ROOT / "Dataset_Processing" / "Real_Data" / "4_checkv_evaluate.py"
COHORT_BP = B.COHORT_BP
DEFAULT_TOOLS = B.DEFAULT_TOOLS

# One row per dataset+method: there are no evaluation levels without truth.
CSV_FIELDS = (
    "dataset",
    "cohort_bp",
    "method",
    "method_class",
    "tool_version",
    "status",
    "effective_min_contig_bp",
    "parameter_policy",
    "metric_definition",
    # what the tool emitted, before CheckV sees it
    "n_bins",
    "multi_contig_bins",
    "singletons",
    "assigned_contigs",
    "total_contigs",
    "assigned_contig_fraction",
    "assigned_bases",
    "total_bases",
    "assigned_base_fraction",
    "largest_bin_contigs",
    # CheckV, on the linked bins
    "predicted_vmags",
    "checkv_complete",
    "checkv_high_quality",
    "checkv_medium_quality",
    "mq_or_better",
    "checkv_low_quality",
    "checkv_not_determined",
    "mean_completeness",
    "median_completeness",
    "completeness_n",
    "mean_contamination",
    "median_contamination",
    "contaminated_over_5pct",
    "checkv_bins_scored",
    "checkv_bins_excluded",
    "runtime_seconds",
    "checkv_seconds",
    "command",
    "notes",
    "run_utc",
)
KEY_FIELDS = ("dataset", "cohort_bp", "method")
CHECKV_COLUMNS = (
    "predicted_vmags", "checkv_complete", "checkv_high_quality",
    "checkv_medium_quality", "mq_or_better", "checkv_low_quality",
    "checkv_not_determined", "mean_completeness", "median_completeness",
    "completeness_n", "mean_contamination", "median_contamination",
    "contaminated_over_5pct",
)


# dataset discovery
def has_ground_truth(dataset_dir: Path) -> bool:
    config = json.loads((dataset_dir / "config.json").read_text(encoding="utf-8"))
    return bool(config.get("has_ground_truth", True))


def discover_real_datasets(root: Path):
    complete, skipped = B.discover_datasets(root)
    real, with_truth = [], []
    for path in complete:
        (real if not has_ground_truth(path) else with_truth).append(path)
    for path in with_truth:
        skipped.append((path.name, "has ground truth; use benchmark.py"))
    return real, skipped


# cohort preparation
def prepare_dataset(dataset_dir: Path, output_dir: Path, minimum: int = COHORT_BP):
    """The same tool inputs benchmark.py builds, minus the truth table.

    benchmark.py's cohort is "contigs >= minimum bp THAT CARRY A LABEL"; here it
    is simply "contigs >= minimum bp", because nothing carries a label.  Every
    other input - cohort FASTA, coverage TSV, CLMB abundance NPZ, CoCoBin
    alias FASTA / k-mer table / mapping - is built to the identical contract so
    the wrappers are unchanged.
    """
    config_path, config, metadata_path, fasta, coverage, graph = \
        B.locate_dataset_inputs(dataset_dir)
    if not fasta.is_file():
        raise FileNotFoundError(f"assembly FASTA not found: {fasta}")
    if not coverage.is_file():
        raise FileNotFoundError(f"coverage table not found: {coverage}")

    sequences = dict(B.read_fasta(fasta))
    lengths = {}
    with metadata_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lengths[row["contig_id"]] = int(row["length"])
    cohort = [name for name, seq in sequences.items() if len(seq) >= minimum]
    if not cohort:
        raise ValueError(f"{dataset_dir.name}: no contigs >= {minimum} bp")

    output_dir.mkdir(parents=True, exist_ok=True)
    cohort_fasta = output_dir / f"contigs.min{minimum}.fasta"
    B.write_fasta(cohort_fasta, ((name, sequences[name]) for name in cohort))

    with coverage.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.reader(source)
        coverage_header = next(reader)
        coverage_rows = {row[0]: row[1:] for row in reader if row}
    missing = [name for name in cohort if name not in coverage_rows]
    if missing:
        raise ValueError(
            f"{dataset_dir.name}: coverage is missing for {len(missing)} cohort contigs")
    abundance = np.asarray(
        [[float(value) for value in coverage_rows[name]] for name in cohort],
        dtype=np.float32)

    coverage_tsv = output_dir / f"coverage.min{minimum}.tsv"
    with coverage_tsv.open("w", newline="", encoding="utf-8") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["contigname", *coverage_header[1:]])
        writer.writerows(([name, *coverage_rows[name]] for name in cohort))

    clmb_rpkm = output_dir / f"clmb.rpkm.min{minimum}.npz"
    np.savez_compressed(clmb_rpkm, abundance)

    # CoCoBin parses `_length_` and `_cov_` out of the header, so it gets stable
    # aliases that are mapped back before anything is scored.
    tnf_path = dataset_dir / "tnf_counts.npz"
    if not tnf_path.is_file():
        raise FileNotFoundError(
            f"CoCoBin k-mer counts not found: {tnf_path}\\n"
            "  the model notebook writes this cache the first time it reads the "
            "dataset; run that cell once, or drop cocobin from --tools")
    with np.load(tnf_path) as archive:
        tnf_by_name = {str(name): row for name, row in
                       zip(archive["names"], np.asarray(archive["counts"], dtype=np.float64))}
    missing_tnf = [name for name in cohort if name not in tnf_by_name]
    if missing_tnf:
        raise ValueError(
            f"{dataset_dir.name}: k-mer counts are missing for {len(missing_tnf)} cohort contigs")
    aliases = [
        f"CB{index:07d}_length_{len(sequences[name])}_cov_{float(abundance[index - 1].mean()):.8g}"
        for index, name in enumerate(cohort, start=1)
    ]
    cocobin_fasta = output_dir / f"cocobin.min{minimum}.fasta"
    B.write_fasta(cocobin_fasta,
                  ((alias, sequences[name]) for alias, name in zip(aliases, cohort)))
    cocobin_mapping = output_dir / f"cocobin.mapping.min{minimum}.tsv"
    with cocobin_mapping.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["alias", "contig"])
        writer.writerows(zip(aliases, cohort))
    cocobin_kmer = output_dir / f"cocobin.kmer.min{minimum}.csv"
    with cocobin_kmer.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        for alias, name in zip(aliases, cohort):
            counts = tnf_by_name[name]
            total = float(counts.sum())
            features = counts / total if total else counts
            writer.writerow([alias, *(f"{value:.10g}" for value in features)])

    manifest = {
        "dataset": dataset_dir.name,
        "config": str(config_path.resolve()),
        "cohort_bp": minimum,
        "cohort_fasta": str(cohort_fasta.resolve()),
        "coverage_tsv": str(coverage_tsv.resolve()),
        "clmb_rpkm": str(clmb_rpkm.resolve()),
        "cocobin_fasta": str(cocobin_fasta.resolve()),
        "cocobin_kmer": str(cocobin_kmer.resolve()),
        "cocobin_mapping": str(cocobin_mapping.resolve()),
        "truth_tsv": None,                     # there is none, and nothing may invent one
        "assembly_graph": str(graph.resolve()) if graph else None,
        "n_samples": int(config.get("n_samples", 0)),
        "threads": int(config.get("cpus", 8)),
        "read_type": config.get("read_type", "unknown"),
        "cohort_contigs": len(cohort),
        "cohort_bases": sum(len(sequences[name]) for name in cohort),
        "has_ground_truth": False,
        "evaluation": "checkv",
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest, {name: lengths.get(name, len(sequences[name])) for name in cohort}


# ---------------------------------------------------------------------------
# CheckV
# ---------------------------------------------------------------------------
def describe_bins(assignments, lengths):
    """Reference-free description of what a tool emitted."""
    by_bin = {}
    for contig, bin_id in assignments.items():
        by_bin.setdefault(bin_id, []).append(contig)
    sizes = [len(members) for members in by_bin.values()]
    assigned_bases = sum(lengths.get(contig, 0) for contig in assignments)
    total_contigs, total_bases = len(lengths), sum(lengths.values())
    return {
        "n_bins": len(by_bin),
        "multi_contig_bins": sum(1 for size in sizes if size > 1),
        "singletons": sum(1 for size in sizes if size == 1),
        "assigned_contigs": len(assignments),
        "total_contigs": total_contigs,
        "assigned_contig_fraction": round(len(assignments) / total_contigs, 6)
        if total_contigs else "",
        "assigned_bases": assigned_bases,
        "total_bases": total_bases,
        "assigned_base_fraction": round(assigned_bases / total_bases, 6)
        if total_bases else "",
        "largest_bin_contigs": max(sizes) if sizes else 0,
    }, by_bin


def write_assignment_table(path: Path, assignments, lengths):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig_id", "bin_id", "length"])
        for contig, bin_id in sorted(assignments.items()):
            writer.writerow([contig, bin_id, lengths.get(contig, "")])
    return path


def run_checkv(*, contigs: Path, bins_tsv: Path, dataset: str, method: str,
               outdir: Path, threads: int, log: Path):
    """Invoke the real-data pipeline's CheckV evaluator; return its metric row."""
    command = [
        sys.executable, str(CHECKV_EVALUATOR),
        "--contigs", str(contigs),
        "--dataset", dataset,
        "--bins", str(bins_tsv),
        "--method", method,
        "--outdir", str(outdir),
        "--threads", str(threads),
    ]
    result, runtime = B.run_capture(command, log=log)
    metrics_csv = outdir / "checkv_metrics.csv"
    if result.returncode or not metrics_csv.is_file():
        return None, runtime, command
    with metrics_csv.open(newline="", encoding="utf-8-sig") as handle:
        rows = [row for row in csv.DictReader(handle) if row.get("method") == method]
    return (rows[-1] if rows else None), runtime, command


def cap_bins(assignments, limit: int):
    """Drop bins above `limit` contigs; return the rest and what was dropped.

    A degenerate bin - thousands of contigs linked into one pseudo-genome tens
    of Mbp long - has no meaningful completeness (a phage genome is 10-200 kb)
    and will exhaust memory in CheckV's hmmsearch stage, taking the whole
    evaluation down with it.  Excluding it keeps the numbers for every other bin
    and records what was removed instead of failing silently.
    """
    by_bin = {}
    for contig, bin_id in assignments.items():
        by_bin.setdefault(bin_id, []).append(contig)
    dropped = {bin_id for bin_id, members in by_bin.items() if len(members) > limit}
    kept = {contig: bin_id for contig, bin_id in assignments.items()
            if bin_id not in dropped}
    dropped_contigs = len(assignments) - len(kept)
    return kept, len(dropped), dropped_contigs


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
def empty_metrics(manifest):
    blank = {name: "" for name in (
        "n_bins", "multi_contig_bins", "singletons", "assigned_contigs",
        "assigned_contig_fraction", "assigned_bases", "assigned_base_fraction",
        "largest_bin_contigs", "checkv_bins_scored", "checkv_bins_excluded",
        "checkv_seconds", *CHECKV_COLUMNS)}
    blank["total_contigs"] = manifest["cohort_contigs"]
    blank["total_bases"] = manifest["cohort_bases"]
    return blank


def append_rows(path: Path, rows, *, replace_existing=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = path.with_suffix(path.suffix + ".lock")
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        existing = []
        if path.is_file() and path.stat().st_size:
            with path.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                if tuple(reader.fieldnames or ()) != CSV_FIELDS:
                    raise ValueError(
                        f"{path} has an incompatible schema; expected {CSV_FIELDS}")
                existing = list(reader)
        incoming_keys = {tuple(str(row[field]) for field in KEY_FIELDS) for row in rows}
        if replace_existing:
            existing = [row for row in existing
                        if tuple(row[field] for field in KEY_FIELDS) not in incoming_keys]
        else:
            seen = {tuple(row[field] for field in KEY_FIELDS) for row in existing}
            rows = [row for row in rows
                    if tuple(str(row[field]) for field in KEY_FIELDS) not in seen]
        if not rows:
            return 0
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(existing + list(rows))
        return len(rows)


def existing_result_keys(path: Path):
    if not path.is_file() or not path.stat().st_size:
        return set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != CSV_FIELDS:
            raise ValueError(f"{path} has an incompatible schema")
        return {tuple(row[field] for field in KEY_FIELDS) for row in reader
                if str(row.get("status", "")).startswith("completed")}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--checkv-root", type=Path, default=DEFAULT_CHECKV_ROOT)
    parser.add_argument("--dataset", action="append",
                        help="dataset directory name; repeatable")
    parser.add_argument("--tools", nargs="+", choices=DEFAULT_TOOLS,
                        default=list(DEFAULT_TOOLS))
    parser.add_argument("--threads", type=int)
    parser.add_argument("--checkv-threads", type=int, default=4,
                        help="threads for CheckV (default 4). CheckV runs this "
                             "many hmmsearch processes at once; the default of "
                             "8 can exhaust memory on a large or badly-binned "
                             "input - the symptom is 'N hmmsearch tasks failed'")
    parser.add_argument("--max-bin-contigs", type=int,
                        help="exclude bins larger than this from CheckV. A bin of "
                             "thousands of contigs is not a genome and will take "
                             "the evaluation down with it; excluding it keeps the "
                             "numbers for every other bin")
    parser.add_argument("--rerun", action="store_true",
                        help="rerun tools and replace matching CSV rows")
    parser.add_argument("--skip-checkv", action="store_true",
                        help="run the tools and write the bins, but do not score")
    parser.add_argument("--bam-dir", type=Path,
                        help="existing sorted cohort BAMs for CoCoNet/vRhyme; "
                             "their references must match the cohort FASTA")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    catalog = json.loads((HERE / "tools.json").read_text(encoding="utf-8"))
    datasets, skipped = discover_real_datasets(args.data_root.resolve())
    selected = set(args.dataset or ())
    if selected:
        datasets = [path for path in datasets if path.name in selected]
        missing = selected.difference(path.name for path in datasets)
        if missing:
            raise SystemExit(f"datasets unavailable, incomplete, or not real: "
                             f"{sorted(missing)}")
    if args.list:
        for path in datasets:
            print(f"ready\t{path.name}")
        for name, reason in skipped:
            print(f"skipped\t{name}\t{reason}")
        return 0
    if not datasets:
        raise SystemExit(
            "no ground-truth-free datasets found; build one with "
            "Dataset_Processing/Real_Data/run_all.sh")
    if not CHECKV_EVALUATOR.is_file() and not args.skip_checkv:
        raise SystemExit(f"CheckV evaluator not found: {CHECKV_EVALUATOR}")

    versions = {}
    recorded = set() if args.rerun else existing_result_keys(args.output.resolve())
    all_rows = []

    for dataset in datasets:
        print(f"\n=== {dataset.name}: >= {COHORT_BP} bp cohort, no ground truth ===",
              flush=True)
        work = args.run_root.resolve() / dataset.name
        manifest, lengths = prepare_dataset(dataset, work / "inputs", COHORT_BP)
        if args.threads:
            manifest["threads"] = args.threads
        metabat_depth = B.build_auxiliary_tables(manifest)
        print(f"cohort: {manifest['cohort_contigs']:,} contigs, "
              f"{manifest['cohort_bases'] / 1e6:,.1f} Mbp, "
              f"{manifest['n_samples']} samples", flush=True)

        bams, bam_note = None, ""
        if args.bam_dir:
            candidates = sorted(args.bam_dir.glob("*.bam"))
            bams = candidates or None
            if bams:
                # Supplied BAMs were produced by some earlier mapping run, not
                # by this pipeline. That is provenance a reader of the CSV needs:
                # same reads and reference or not, it is a different alignment
                # than the one behind the coverage matrix.
                bam_note = (f"; coverage BAMs supplied via --bam-dir "
                            f"{args.bam_dir} ({len(bams)} sample(s)), not "
                            "generated by this pipeline")
                print(f"using {len(bams)} BAM(s) from {args.bam_dir}", flush=True)
        elif not args.dry_run:
            bams = B.existing_bams(dataset, manifest)

        checkv_out = args.checkv_root.resolve() / dataset.name

        # ---- the external tools ------------------------------------------
        jobs = []
        for key in args.tools:
            spec = catalog[key]
            method = spec["method"]
            if (dataset.name, str(COHORT_BP), method) in recorded:
                print(f"{method}: already_recorded", flush=True)
                continue
            if key not in versions:
                versions[key] = B.tool_version(spec)
            tool_out = work / "results" / key
            tool_out.parent.mkdir(parents=True, exist_ok=True)
            log = work / "logs" / f"{key}.log"
            wrapper = HERE / "tools" / f"{key}.sh"
            source = B.assignment_source(key, tool_out)
            status, runtime = "completed", 0.0
            notes = spec["notes"]

            if args.dry_run:
                status = "dry_run"
            elif source is not None and not args.rerun:
                status = "completed_reused"
                runtime = B.logged_runtime(tool_out / "log.txt") or B.logged_runtime(log)
                notes += "; reused existing assignment output"
            elif spec.get("availability") == "no_public_implementation":
                status = "no_public_implementation"
            elif not Path(spec["version_command"][0]).is_file():
                status = "not_installed"
            elif key in {"coconet", "vrhyme"} and not bams:
                status = "missing_bam_input"
                notes += ("; this dataset has no reads on disk, so cohort BAMs "
                          "cannot be built; pass --bam-dir to supply them")
            elif tool_out.is_dir() and any(tool_out.iterdir()) and not args.rerun:
                status = "incomplete_existing_output"
                notes += ("; nonempty output has no assignment table; preserved "
                          "unchanged; inspect its log and use --rerun explicitly")
            else:
                if tool_out.exists() and args.rerun:
                    shutil.rmtree(tool_out)
                env = os.environ.copy()
                env.update({
                    "COHORT_FASTA": manifest["cohort_fasta"],
                    "COVERAGE_TSV": manifest["coverage_tsv"],
                    "CLMB_RPKM": manifest["clmb_rpkm"],
                    "COCOBIN_FASTA": manifest["cocobin_fasta"],
                    "COCOBIN_KMER": manifest["cocobin_kmer"],
                    "COCOBIN_MAPPING": manifest["cocobin_mapping"],
                    "METABAT_DEPTH": str(metabat_depth),
                    "TOOL_OUT": str(tool_out),
                    "THREADS": str(manifest["threads"]),
                    "ASSEMBLY_GRAPH": manifest.get("assembly_graph") or "",
                    "BAM_FILES": "\n".join(map(str, bams or ())),
                })
                result, runtime = B.run_capture(["bash", str(wrapper)], env=env, log=log)
                if result.returncode:
                    status = {4: "missing_required_input",
                              5: "not_installed"}.get(result.returncode, "failed")
                    notes += f"; runner exit {result.returncode}; see {log.relative_to(ROOT)}"
                source = B.assignment_source(key, tool_out)
                if status == "completed" and source is None:
                    if key == "metabat2" and (tool_out / "bins").is_dir():
                        source = tool_out / "bins"
                    else:
                        status = "failed_missing_assignments"
            if key in {"coconet", "vrhyme"} and bams and bam_note:
                notes += bam_note
            jobs.append((method, spec["method_class"], versions[key],
                         max(COHORT_BP, int(spec["hard_minimum_bp"])),
                         status, runtime, f"bash {wrapper}", notes, source))
            print(f"{method}: {status}", flush=True)

        # ---- score everything with CheckV --------------------------------
        for (method, method_class, version, effective_min, status, runtime,
             command_text, notes, source) in jobs:
            row = {
                "dataset": dataset.name,
                "cohort_bp": COHORT_BP,
                "method": method,
                "method_class": method_class,
                "tool_version": version,
                "status": status,
                "effective_min_contig_bp": effective_min,
                "parameter_policy": "common_2000bp",
                "metric_definition": "checkv_linked_bins_1500N",
                "runtime_seconds": round(runtime, 3),
                "command": command_text,
                "notes": notes,
                "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            }
            row.update(empty_metrics(manifest))
            if source is None or not str(status).startswith("completed"):
                all_rows.append(row)
                continue

            assignments = (
                B.assignments_from_fasta_dir(source, set(lengths))
                if source.is_dir()
                else B.assignments_from_table(source, set(lengths)))
            description, _ = describe_bins(assignments, lengths)
            row.update(description)

            if args.skip_checkv or args.dry_run:
                row["notes"] += "; CheckV not run"
                all_rows.append(row)
                continue

            scored = assignments
            excluded_bins = 0
            if args.max_bin_contigs:
                scored, excluded_bins, excluded_contigs = cap_bins(
                    assignments, args.max_bin_contigs)
                if excluded_bins:
                    row["notes"] += (
                        f"; {excluded_bins} bin(s) over {args.max_bin_contigs} "
                        f"contigs ({excluded_contigs} contigs) excluded from CheckV")
            if not scored:
                row["status"] = "checkv_no_scoreable_bins"
                all_rows.append(row)
                continue

            bins_tsv = write_assignment_table(
                work / "checkv_input" / f"{method.replace(' ', '_').replace('/', '_')}.tsv",
                scored, lengths)
            metrics, checkv_seconds, _ = run_checkv(
                contigs=Path(manifest["cohort_fasta"]), bins_tsv=bins_tsv,
                dataset=dataset.name, method=method, outdir=checkv_out,
                threads=args.checkv_threads,
                log=work / "logs" / f"checkv.{method.replace(' ', '_')}.log")
            row["checkv_seconds"] = round(checkv_seconds, 3)
            row["checkv_bins_excluded"] = excluded_bins
            if metrics is None:
                row["status"] = "checkv_failed"
                row["notes"] += ("; CheckV failed - a very large bin can exhaust "
                                 "memory in hmmsearch; retry with "
                                 "--max-bin-contigs or fewer --checkv-threads")
            else:
                row["checkv_bins_scored"] = metrics.get("predicted_vmags", "")
                for column in CHECKV_COLUMNS:
                    row[column] = metrics.get(column, "")
            print(f"  CheckV {method}: {row['status']}", flush=True)
            all_rows.append(row)

    if args.dry_run:
        print("\nDry run complete; comparison CSV was not modified.")
        return 0
    appended = append_rows(args.output.resolve(), all_rows, replace_existing=True)
    print(f"\nWrote {appended} rows to {args.output.resolve()}")
    print("Rows for other datasets and methods were preserved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
