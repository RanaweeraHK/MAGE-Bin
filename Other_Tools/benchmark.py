#!/usr/bin/env python3
"""Run the common 2-kb external-tool benchmark and append new result rows."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_DATA_ROOT = ROOT / "Data" / "Processed_data"
DEFAULT_OUTPUT = ROOT / "Outputs" / "comparison_simulated.csv"
DEFAULT_RUN_ROOT = ROOT / "Outputs" / "other_tools"

COHORT_BP = 2000
DEFAULT_TOOLS = (
    "coconet", "vrhyme", "vamb", "cocobin", "clmb",
    "metabat2", "concoct",
)
LEVELS = ("strain", "species")
CSV_FIELDS = (
    "dataset",
    "cohort_bp",
    "evaluation_level",
    "method",
    "method_class",
    "tool_version",
    "status",
    "effective_min_contig_bp",
    "parameter_policy",
    "metric_definition",
    "precision",
    "recall",
    "f1",
    "ari",
    "nmi",
    "hq_bins",
    "n_bins",
    "multi_contig_bins",
    "singletons",
    "assigned_contigs",
    "total_contigs",
    "assigned_contig_fraction",
    "assigned_bases",
    "total_bases",
    "assigned_base_fraction",
    "runtime_seconds",
    "command",
    "notes",
    "run_utc",
)
KEY_FIELDS = ("dataset", "cohort_bp", "evaluation_level", "method")
METHOD_ALIASES = {}
FASTA_SUFFIXES = {".fa", ".fna", ".fas", ".fasta"}

# reads FASTA file , allow both contigs.fasta and contigs.fasta.gz
def read_fasta(path: Path):
    name, sequence = None, []
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, mode="rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name, sequence = line[1:].split()[0], []
            elif name is None:
                raise ValueError(f"{path}: sequence before first FASTA header")
            else:
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence)

# writes sequences into a FASTA file
def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

# resolve the relative paths in contig.json
def resolve_config_path(config_path: Path, raw: str) -> Path:
    path = Path(raw).expanduser()
    if path.exists():
        return path.resolve()
    if not path.is_absolute():
        candidate = (config_path.parent / path).resolve()
        if candidate.exists():
            return candidate
    return path

# chceck the files used for each dataset. if even one file is missing, dataset is rejected.
def locate_dataset_inputs(dataset_dir: Path):
    config_path = dataset_dir / "config.json"
    metadata_path = dataset_dir / "contig_metadata.tsv"
    if not config_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError("config.json or contig_metadata.tsv is missing")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    fasta = resolve_config_path(config_path, config.get("combined_fasta", ""))
    work = resolve_config_path(config_path, config.get("work_dir", ""))
    if not fasta.is_file():
        fasta = next((p for p in (
            work / "contigs_filt.fasta",
            work / str(config.get("contig_mode", "")) / "contigs_filt.fasta",
        ) if p.is_file()), fasta)
    coverage = work / "coverage.csv"
    if not coverage.is_file():
        coverage = work / str(config.get("contig_mode", "")) / "coverage.csv"
    graph = next((p for p in (
        work / "assembly_graph.gfa",
        work / "assembly_graph_with_scaffolds.gfa",
        work / "final.contigs.gfa",
        work.parent / "assembly" / "assembly_graph.gfa",
    ) if p.is_file()), None)
    return config_path, config, metadata_path, fasta, coverage, graph

# load contig_metadata.tsv, -> it should contain contig_id, length , genome_label
def load_metadata(path: Path):
    rows = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig_id", "length", "genome_label"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in reader:
            genome = row["genome_label"].strip()
            if genome and genome.lower() not in {"nan", "none", "unmapped"}:
                rows[row["contig_id"]] = (int(row["length"]), genome)
    return rows

# converts starin labels into species label-> starin : virus_25_s1 - species : virus_25
def parent_species(genome: str) -> str:
    return re.sub(r"_s\d+$", "", genome)

### prepare all the input representations needed by the different tools
def prepare_dataset(dataset_dir: Path, output_dir: Path, minimum: int = 2000):
    config_path, config, metadata_path, fasta, coverage, graph = locate_dataset_inputs(dataset_dir)
    if not fasta.is_file():
        raise FileNotFoundError(f"assembly FASTA not found: {fasta}")
    if not coverage.is_file():
        raise FileNotFoundError(f"coverage table not found: {coverage}")
    sequences, metadata = dict(read_fasta(fasta)), load_metadata(metadata_path)
    cohort = [name for name, seq in sequences.items() if len(seq) >= minimum and name in metadata]
    if not cohort:
        raise ValueError(f"{dataset_dir.name}: no labelled contigs >= {minimum} bp")
    output_dir.mkdir(parents=True, exist_ok=True)
    cohort_fasta = output_dir / f"contigs.min{minimum}.fasta"
    write_fasta(cohort_fasta, ((name, sequences[name]) for name in cohort))
    truth_path = output_dir / f"truth.min{minimum}.tsv"
    with truth_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "length", "strain", "species"])
        for name in cohort:
            length, strain = metadata[name]
            writer.writerow([name, length, strain, parent_species(strain)])
    with coverage.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.reader(source)
        coverage_header = next(reader)
        coverage_rows = {row[0]: row[1:] for row in reader if row}
    missing_coverage = [name for name in cohort if name not in coverage_rows]
    if missing_coverage:
        raise ValueError(
            f"{dataset_dir.name}: coverage is missing for {len(missing_coverage)} cohort contigs"
        )
    abundance = np.asarray(
        [[float(value) for value in coverage_rows[name]] for name in cohort],
        dtype=np.float32,
    )
    coverage_tsv = output_dir / f"coverage.min{minimum}.tsv"
    with coverage_tsv.open("w", newline="", encoding="utf-8") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["contigname", *coverage_header[1:]])
        writer.writerows(([name, *coverage_rows[name]] for name in cohort))

    # CLMB's public CLI accepts a float32 abundance NPZ in FASTA order.
    clmb_rpkm = output_dir / f"clmb.rpkm.min{minimum}.npz"
    np.savez_compressed(clmb_rpkm, abundance)

    # CoCoBin's public script parses ``_length_`` and ``_cov_`` from headers.
    # Stable aliases avoid modifying source contig identifiers and are mapped
    # back before evaluation. Its k-mer input is the canonical 4-mer count
    # matrix already generated for this exact processed assembly.
    tnf_path = dataset_dir / "tnf_counts.npz"
    if not tnf_path.is_file():
        raise FileNotFoundError(f"CoCoBin k-mer counts not found: {tnf_path}")
    with np.load(tnf_path) as archive:
        tnf_names = [str(name) for name in archive["names"]]
        tnf_counts = np.asarray(archive["counts"], dtype=np.float64)
    tnf_by_name = {name: row for name, row in zip(tnf_names, tnf_counts)}
    missing_tnf = [name for name in cohort if name not in tnf_by_name]
    if missing_tnf:
        raise ValueError(
            f"{dataset_dir.name}: k-mer counts are missing for {len(missing_tnf)} cohort contigs"
        )
    aliases = [
        f"CB{index:07d}_length_{len(sequences[name])}_cov_{float(abundance[index - 1].mean()):.8g}"
        for index, name in enumerate(cohort, start=1)
    ]
    cocobin_fasta = output_dir / f"cocobin.min{minimum}.fasta"
    write_fasta(cocobin_fasta, (
        (alias, sequences[name]) for alias, name in zip(aliases, cohort)
    ))
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
        "truth_tsv": str(truth_path.resolve()),
        "assembly_graph": str(graph.resolve()) if graph else None,
        "n_samples": int(config.get("n_samples", 0)),
        "threads": int(config.get("cpus", 8)),
        "read_type": config.get("read_type", "unknown"),
        "cohort_contigs": len(cohort),
        "cohort_bases": sum(len(sequences[name]) for name in cohort),
        "strain_genomes": len({metadata[name][1] for name in cohort}),
        "species_genomes": len({parent_species(metadata[name][1]) for name in cohort}),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest

# load ground truth
def load_truth(path: Path, level: str):
    truth, lengths = {}, {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            truth[row["contig"]], lengths[row["contig"]] = row[level], int(row["length"])
    return truth, lengths

# Some tools output bins as FASTA files , Metabat2 -> converts this into 
# {
#     "contig1": "bin.1",
#     "contig5": "bin.1",
#     "contig8": "bin.2"
# } format
def assignments_from_fasta_dir(path: Path, valid: set[str]):
    assignments = {}
    fastas = sorted(
        p for p in path.rglob("*")
        if p.is_file() and (
            p.suffix.lower() in FASTA_SUFFIXES
            or any(p.name.lower().endswith(suffix + ".gz") for suffix in FASTA_SUFFIXES)
        )
    )
    for fasta in fastas:
        bin_name = str(fasta.relative_to(path).with_suffix(""))
        for contig, _ in read_fasta(fasta):
            if contig in valid:
                previous = assignments.setdefault(contig, bin_name)
                if previous != bin_name:
                    raise ValueError(f"{contig} occurs in two bins")
    return assignments

# Other tools output tables instead of FASTA bins ->  convert this into 
# {
#     "contig1": "bin.1",
#     "contig5": "bin.1",
#     "contig8": "bin.2"
# } format
def assignments_from_table(path: Path, valid: set[str]):
    assignments = {}
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.replace(",", "\t").split()
            positions = [i for i, value in enumerate(fields) if value in valid]
            if not positions:
                continue
            index, contig = positions[0], fields[positions[0]]
            candidates = [value for i, value in enumerate(fields) if i != index]
            if not candidates:
                continue
            previous = assignments.setdefault(contig, candidates[0])
            if previous != candidates[0]:
                raise ValueError(f"{contig} occurs in two bins")
    return assignments


def comb2(value: int) -> int:
    return value * (value - 1) // 2

# calculates ARI without using scikit-learn
def adjusted_rand_index(truth, predicted):
    n = len(truth)
    if n < 2:
        return 1.0
    cells = Counter(zip(truth, predicted))
    true_counts, pred_counts = Counter(truth), Counter(predicted)
    observed = sum(comb2(value) for value in cells.values())
    true_sum = sum(comb2(value) for value in true_counts.values())
    pred_sum = sum(comb2(value) for value in pred_counts.values())
    total = comb2(n)
    expected, maximum = true_sum * pred_sum / total, (true_sum + pred_sum) / 2
    return (observed - expected) / (maximum - expected) if maximum != expected else 1.0


def normalized_mutual_info(truth, predicted):
    n = len(truth)
    cells = Counter(zip(truth, predicted))
    true_counts, pred_counts = Counter(truth), Counter(predicted)
    mutual = sum(
        count / n * math.log(count * n / (true_counts[t] * pred_counts[p]))
        for (t, p), count in cells.items()
    )
    h_true = -sum(count / n * math.log(count / n) for count in true_counts.values())
    h_pred = -sum(count / n * math.log(count / n) for count in pred_counts.values())
    denominator = (h_true + h_pred) / 2
    return mutual / denominator if denominator else 1.0

# main evaluation function
def evaluate_assignments(truth, lengths, emitted):
    """GraphBin2 contig-count PRF plus ARI/NMI and truth-based HQ bins.

    GraphBin2 defines a K-by-S matrix a[k,s] of binned contigs, where K is
    predicted bins and S is truth species/genomes. Precision is the sum of the
    dominant truth count in every predicted bin divided by binned contigs.
    Recall is the sum of the best recovered-bin count for every truth group
    divided by all cohort contigs, including unclassified contigs.
    """
    contigs = list(truth)
    assigned = set(emitted).intersection(truth)
    matrix = Counter((emitted[c], truth[c]) for c in assigned)
    by_bin, by_truth = defaultdict(Counter), defaultdict(Counter)
    for (predicted_bin, truth_group), count in matrix.items():
        by_bin[predicted_bin][truth_group] += count
        by_truth[truth_group][predicted_bin] += count
    precision_numerator = sum(max(counts.values()) for counts in by_bin.values())
    recall_numerator = sum(
        max(by_truth[group].values(), default=0) for group in set(truth.values())
    )
    precision = precision_numerator / len(assigned) if assigned else 0.0
    recall = recall_numerator / len(contigs) if contigs else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0

    # ARI/NMI keep unclassified contigs as distinct singleton predictions so
    # discarded contigs cannot improve those secondary clustering metrics.
    predicted = {c: emitted.get(c, f"__unassigned_singleton__{c}") for c in contigs}
    pred_counts, true_counts = Counter(predicted.values()), Counter(truth.values())
    bin_members = defaultdict(Counter)
    for contig in assigned:
        bin_members[emitted[contig]][truth[contig]] += 1
    hq_bins = 0
    for counts in bin_members.values():
        dominant, hits = counts.most_common(1)[0]
        if hits / true_counts[dominant] >= 0.90 and 1.0 - hits / sum(counts.values()) <= 0.05:
            hq_bins += 1
    truth_labels, pred_labels = [truth[c] for c in contigs], [predicted[c] for c in contigs]
    assigned_bases, total_bases = sum(lengths[c] for c in assigned), sum(lengths.values())
    return {
        "precision": precision, "recall": recall, "f1": f1,
        "ari": adjusted_rand_index(truth_labels, pred_labels),
        "nmi": normalized_mutual_info(truth_labels, pred_labels),
        "hq_bins": hq_bins, "n_bins": len(by_bin),
        "multi_contig_bins": sum(sum(counts.values()) >= 2 for counts in by_bin.values()),
        "singletons": sum(sum(counts.values()) == 1 for counts in by_bin.values()),
        "assigned_contigs": len(assigned), "total_contigs": len(contigs),
        "assigned_contig_fraction": len(assigned) / len(contigs),
        "assigned_bases": assigned_bases, "total_bases": total_bases,
        "assigned_base_fraction": assigned_bases / total_bases,
    }

# searhes : Data/Processed_data
def discover_datasets(root: Path):
    complete, skipped = [], []
    if not root.is_dir():
        raise FileNotFoundError(root)
    for path in sorted(item for item in root.iterdir() if item.is_dir()):
        required = (path / "config.json", path / "contig_metadata.tsv")
        if not all(item.is_file() for item in required):
            skipped.append((path.name, "missing config.json or contig_metadata.tsv"))
            continue
        try:
            _, _, _, fasta, coverage, _ = locate_dataset_inputs(path)
            if not fasta.is_file() or not coverage.is_file():
                raise FileNotFoundError("assembly FASTA or coverage.csv is missing")
        except Exception as error:  # discovery should not abort other datasets
            skipped.append((path.name, str(error)))
            continue
        complete.append(path)
    return complete, skipped

# common functions to run external commands Ex: bash vamb.sh -> outputs at Outputs/other_tools/dataset/logs/vamb.log
def run_capture(command, *, env=None, log=None, timeout=None):
    started = time.monotonic()
    if log is None:
        result = subprocess.run(
            command, cwd=ROOT, env=env, text=True, capture_output=True, timeout=timeout
        )
    else:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w", encoding="utf-8") as handle:
            result = subprocess.run(
                command,
                cwd=ROOT,
                env=env,
                text=True,
                stdout=handle,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )
    return result, time.monotonic() - started

# for the reproducability , identifies tool version used from tools.json
def tool_version(spec):
    if spec.get("availability") == "no_public_implementation":
        return "no_public_implementation"
    command = spec["version_command"]
    executable = Path(command[0])
    if not executable.is_file():
        return "not_installed"
    try:
        env = os.environ.copy()
        env["PATH"] = str(executable.parent) + os.pathsep + env.get("PATH", "")
        result, _ = run_capture(command, env=env, timeout=20)
        text = (result.stdout + "\n" + result.stderr).strip()
    except Exception as error:
        return f"unavailable ({type(error).__name__})"
    commit = re.search(r"(?i)\bcommit\s+([0-9a-f]{7,40})\b", text)
    if commit:
        return f"git-{commit.group(1)[:12]}"
    match = re.search(r"(?i)(?:version\s*[:=]?\s*|v)(\d+(?:\.\d+)+)", text)
    if not match:
        match = re.search(r"\b(\d+(?:\.\d+)+)\b", text)
    return match.group(1) if match else "installed_version_not_reported"


def logged_runtime(path: Path) -> float:
    """Recover a tool-reported runtime when assignments are evaluated later."""
    if not path.is_file():
        return 0.0
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(r"Completed Vamb in ([0-9.]+) seconds", text)
    return float(matches[-1]) if matches else 0.0

#  model specific converage file to metabat2 
def build_auxiliary_tables(manifest):
    coverage = Path(manifest["coverage_tsv"])
    output_dir = coverage.parent
    metabat = output_dir / "coverage.metabat2.tsv"
    with coverage.open(newline="", encoding="utf-8") as source:
        rows = list(csv.reader(source, delimiter="\t"))
    with metabat.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    return metabat

### Some methods need BAM files rather than only a coverage table
def build_bams(dataset_dir: Path, manifest, *, force=False):
    """Create sorted cohort BAMs only for BAM-dependent tools."""
    config = json.loads((dataset_dir / "config.json").read_text(encoding="utf-8"))
    samples = config.get("samples", [])
    if not samples:
        raise ValueError("config has no samples for BAM generation")
    bam_dir = Path(manifest["cohort_fasta"]).parent / "bam"
    bam_dir.mkdir(parents=True, exist_ok=True)
    minimap2 = Path(os.environ.get(
        "MINIMAP2_CMD", "/home/hashini/miniforge3/envs/asm_env/bin/minimap2"
    ))
    samtools = Path(os.environ.get(
        "SAMTOOLS_CMD", "/home/hashini/miniforge3/envs/asm_env/bin/samtools"
    ))
    if not minimap2.is_file() or not samtools.is_file():
        raise FileNotFoundError("minimap2 or samtools is unavailable")
    preset = config.get("minimap_preset", "sr")
    fasta = manifest["cohort_fasta"]
    bams = []
    for index, sample in enumerate(samples):
        bam = bam_dir / f"sample{index}.bam"
        bams.append(bam)
        if bam.is_file() and Path(str(bam) + ".bai").is_file() and not force:
            continue
        reads = [sample.get("r1"), sample.get("r2")]
        reads = [str(Path(item).expanduser()) for item in reads if item]
        if not reads or not all(Path(item).is_file() for item in reads):
            raise FileNotFoundError(f"reads missing for sample {sample.get('name', index)}")
        # Avoid shell pipelines: stream minimap2 directly into samtools sort.
        mapper = subprocess.Popen(
            [str(minimap2), "-x", preset, "-a", "-t", str(manifest["threads"]), fasta, *reads],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        sort = subprocess.run(
            [str(samtools), "sort", "-@", str(manifest["threads"]), "-o", str(bam), "-"],
            cwd=ROOT,
            stdin=mapper.stdout,
        )
        if mapper.stdout:
            mapper.stdout.close()
        mapper_status = mapper.wait()
        if mapper_status or sort.returncode:
            raise RuntimeError(f"alignment failed for {sample.get('name', index)}")
        subprocess.run([str(samtools), "index", str(bam)], check=True, cwd=ROOT)
    return bams


def existing_bams(dataset_dir: Path, manifest):
    """Return the complete retained BAM set, never a partial sample set."""
    config = json.loads((dataset_dir / "config.json").read_text(encoding="utf-8"))
    expected = len(config.get("samples", []))
    bam_dir = Path(manifest["cohort_fasta"]).parent / "bam"
    bams = [bam_dir / f"sample{index}.bam" for index in range(expected)]
    if expected and all(
        bam.is_file() and Path(str(bam) + ".bai").is_file() for bam in bams
    ):
        return bams
    return None

# different tools has different assignments 
def assignment_source(tool: str, output: Path):
    if tool == "coconet":
        return next(iter(sorted(output.glob("bins_*.csv"))), None) or next(
            iter(sorted(output.glob("bins-*.csv"))), None
        )
    if tool == "vrhyme":
        return next(iter(sorted(output.glob("vRhyme_best_bins*.membership.tsv"))), None)
    if tool == "vamb":
        names = ("vae_clusters_unsplit.tsv", "vae_clusters_split.tsv", "clusters.tsv")
        for name in names:
            if (output / name).is_file():
                return output / name
        return None
    if tool == "cocobin":
        return output / "assignments.tsv" if (output / "assignments.tsv").is_file() else None
    if tool == "clmb":
        return output / "clusters.tsv" if (output / "clusters.tsv").is_file() else None
    if tool == "metabat2":
        return output / "bins" if (output / "bins").is_dir() else None
    if tool == "concoct":
        return next(iter(sorted(output.glob("clustering_*.csv"))), None)
    return None

# append results into csv
def append_rows(path: Path, rows, *, replace_existing=False):
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
                        f"{path} has an incompatible schema; expected {CSV_FIELDS}"
                    )
                existing = list(reader)
                for row in existing:
                    row["method"] = METHOD_ALIASES.get(row["method"], row["method"])
        existing_keys = {tuple(row[field] for field in KEY_FIELDS) for row in existing}
        incoming = [
            row for row in rows
            if replace_existing or tuple(str(row[field]) for field in KEY_FIELDS) not in existing_keys
        ]
        if replace_existing:
            incoming_keys = {tuple(str(row[field]) for field in KEY_FIELDS) for row in incoming}
            existing = [
                row for row in existing
                if tuple(row[field] for field in KEY_FIELDS) not in incoming_keys
            ]
        if not incoming:
            return 0
        mode = "a" if path.is_file() and path.stat().st_size and not replace_existing else "w"
        output_rows = incoming if mode == "a" else existing + incoming
        with path.open(mode, newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
            if mode == "w":
                writer.writeheader()
            writer.writerows(output_rows)
        return len(incoming)


def existing_result_keys(path: Path):
    if not path.is_file() or not path.stat().st_size:
        return set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != CSV_FIELDS:
            raise ValueError(f"{path} has an incompatible schema; expected {CSV_FIELDS}")
        return {
            tuple(row[field] for field in KEY_FIELDS)
            for row in reader
            if str(row.get("status", "")).startswith("completed")
        }


def empty_metrics(manifest):
    return {
        name: "" for name in (
            "precision", "recall", "f1", "ari", "nmi", "hq_bins", "n_bins",
            "multi_contig_bins", "singletons", "assigned_contigs", "assigned_contig_fraction",
            "assigned_bases", "assigned_base_fraction",
        )
    } | {"total_contigs": manifest["cohort_contigs"], "total_bases": manifest["cohort_bases"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--dataset", action="append", help="dataset directory name; repeatable")
    parser.add_argument("--tools", nargs="+", choices=DEFAULT_TOOLS, default=list(DEFAULT_TOOLS))
    parser.add_argument("--threads", type=int)
    parser.add_argument("--rerun", action="store_true", help="rerun tools and replace matching CSV rows")
    parser.add_argument(
        "--refresh-csv",
        action="store_true",
        help="re-evaluate existing raw assignments and replace matching CSV rows without rerunning tools",
    )
    parser.add_argument(
        "--prepare-bams", action="store_true",
        help="generate sorted BAMs needed by BAM-dependent tools",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    catalog = json.loads((HERE / "tools.json").read_text(encoding="utf-8"))
    datasets, skipped = discover_datasets(args.data_root.resolve())
    selected = set(args.dataset or ())
    if selected:
        datasets = [path for path in datasets if path.name in selected]
        missing = selected.difference(path.name for path in datasets)
        if missing:
            raise SystemExit(f"datasets unavailable or incomplete: {sorted(missing)}")
    if args.list:
        for path in datasets:
            print(f"ready\t{path.name}")
        for name, reason in skipped:
            print(f"skipped\t{name}\t{reason}")
        return 0
    if not datasets:
        raise SystemExit("no complete processed datasets found")
    versions = {}
    recorded = (
        set() if args.rerun or args.refresh_csv
        else existing_result_keys(args.output.resolve())
    )

    all_rows = []
    for dataset in datasets:
        print(f"\n=== {dataset.name}: common >= {COHORT_BP} bp cohort ===", flush=True)
        work = args.run_root.resolve() / dataset.name
        manifest = prepare_dataset(dataset, work / "inputs", COHORT_BP)
        if args.threads:
            manifest["threads"] = args.threads
        metabat_depth = build_auxiliary_tables(manifest)
        pending_bam_tools = {
            key for key in {"coconet", "vrhyme"}.intersection(args.tools)
            if not {
                (dataset.name, str(COHORT_BP), level, catalog[key]["method"])
                for level in LEVELS
            }.issubset(recorded)
        }
        bams = existing_bams(dataset, manifest)
        if pending_bam_tools and args.prepare_bams and not args.dry_run:
            bams = build_bams(dataset, manifest, force=False)

        for key in args.tools:
            spec = catalog[key]
            method = spec["method"]
            expected_keys = {
                (dataset.name, str(COHORT_BP), level, method) for level in LEVELS
            }
            if expected_keys.issubset(recorded):
                print(f"{method}: already_recorded", flush=True)
                continue
            if key not in versions:
                versions[key] = tool_version(spec)
            effective_min = max(COHORT_BP, int(spec["hard_minimum_bp"]))
            tool_out = work / "results" / key
            tool_out.parent.mkdir(parents=True, exist_ok=True)
            log = work / "logs" / f"{key}.log"
            wrapper = HERE / "tools" / f"{key}.sh"
            source = assignment_source(key, tool_out)
            status, runtime, command_text = "completed", 0.0, f"bash {wrapper}"
            notes = spec["notes"]
            provenance = tool_out / "PROVENANCE.txt"
            if provenance.is_file():
                notes += "; " + provenance.read_text(encoding="utf-8").strip()

            if args.dry_run:
                status = "dry_run"
            elif source is not None and not args.rerun:
                status = "completed_reused"
                runtime = logged_runtime(tool_out / "log.txt") or logged_runtime(log)
                notes += "; reused existing assignment output"
            elif spec.get("availability") == "no_public_implementation":
                status = "no_public_implementation"
            elif not Path(spec["version_command"][0]).is_file():
                status = "not_installed"
            elif key in {"coconet", "vrhyme"} and bams is None:
                status = "missing_bam_input"
                notes += "; rerun with --prepare-bams to generate sorted cohort BAMs; synthetic coverage variance is not used"
            elif args.refresh_csv:
                if tool_out.is_dir() and any(tool_out.iterdir()):
                    status = "incomplete_existing_output"
                    notes += "; retained output has no final assignment table; tool was not rerun"
                else:
                    status = "missing_existing_output"
                    notes += "; no retained assignment output exists; tool was not run"
            elif tool_out.is_dir() and any(tool_out.iterdir()) and not args.rerun:
                status = "incomplete_existing_output"
                notes += "; nonempty output has no assignment table; preserved unchanged; inspect its log and use --rerun explicitly"
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
                result, runtime = run_capture(["bash", str(wrapper)], env=env, log=log)
                if result.returncode:
                    status = {
                        4: "missing_required_input", 5: "not_installed"
                    }.get(result.returncode, "failed")
                    notes += f"; runner exit {result.returncode}; see {log.relative_to(ROOT)}"
                source = assignment_source(key, tool_out)
                if status == "completed" and source is None:
                    if key == "metabat2" and (tool_out / "bins").is_dir():
                        source = tool_out / "bins"
                    else:
                        status = "failed_missing_assignments"

            version = versions[key]
            for level in LEVELS:
                base = {
                    "dataset": dataset.name,
                    "cohort_bp": COHORT_BP,
                    "evaluation_level": level,
                    "method": method,
                    "method_class": spec["method_class"],
                    "tool_version": version,
                    "status": status,
                    "effective_min_contig_bp": effective_min,
                    "parameter_policy": "common_2000bp",
                    "metric_definition": "graphbin2_contig_count",
                    "runtime_seconds": round(runtime, 3),
                    "command": command_text,
                    "notes": notes,
                    "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                }
                if source is not None and status.startswith("completed"):
                    truth, lengths = load_truth(Path(manifest["truth_tsv"]), level)
                    emitted = (
                        assignments_from_fasta_dir(source, set(truth))
                        if source.is_dir()
                        else assignments_from_table(source, set(truth))
                    )
                    metrics = evaluate_assignments(truth, lengths, emitted)
                    base.update({key: round(value, 8) if isinstance(value, float) else value for key, value in metrics.items()})
                else:
                    base.update(empty_metrics(manifest))
                all_rows.append(base)
            print(f"{method}: {status}", flush=True)
    if args.dry_run:
        print("\nDry run complete; comparison CSV was not modified.")
        return 0
    appended = append_rows(
        args.output.resolve(),
        all_rows,
        replace_existing=True,
    )
    print(f"\nWrote {appended} rows to {args.output.resolve()}")
    print("Existing rows for other datasets were preserved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
