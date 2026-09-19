"""End-to-end MAGE-Bin inference and reproducible result writing."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .clustering import cluster_unknown_number_of_bins
from .completeness import detect_complete_contigs, read_fasta
from .config import MageBinConfig
from .dataset import load_binning_dataset
from .graph import filter_supported_assembly_edges, load_assembly_edges
from .model import train_mage_embedding


@dataclass(frozen=True, slots=True)
class BinningResult:
    """Paths and summary values produced by one MAGE-Bin run."""

    assignment_file: Path
    run_metadata_file: Path
    contig_count: int
    bin_count: int
    singleton_count: int
    bin_fasta_directory: Path | None = None


def _json_default(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialize {type(value).__name__}")


def _write_bin_fastas(
    fasta: str | Path, output: Path, contig_ids: list[str], labels: np.ndarray
) -> Path:
    """Write one FASTA per predicted bin using the original contig sequences."""

    source = Path(fasta).expanduser()
    if not source.is_file():
        raise FileNotFoundError(f"source contig FASTA not found: {source}")
    bins = output / "bins"
    bins.mkdir(parents=True, exist_ok=True)
    for previous in bins.glob("MAGE-Bin_*.fasta"):
        previous.unlink()
    assignments = dict(zip(contig_ids, labels, strict=True))
    found: set[str] = set()
    for name, sequence in read_fasta(source):
        label = assignments.get(name)
        if label is None:
            continue
        found.add(name)
        with (bins / f"MAGE-Bin_{int(label):06d}.fasta").open("a") as handle:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    missing = assignments.keys() - found
    if missing:
        example = ", ".join(sorted(missing)[:3])
        raise ValueError(
            f"source FASTA is missing {len(missing)} binned contigs: {example}"
        )
    return bins


def run_binning(
    dataset_directory: str | Path,
    output_directory: str | Path,
    *,
    config: MageBinConfig | None = None,
    coverage_file: str | Path | None = None,
    device: str = "auto",
) -> BinningResult:
    """Run MAGE-Bin and write assignments plus complete run diagnostics."""

    began = time.perf_counter()
    parameters = config or MageBinConfig()
    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    dataset = load_binning_dataset(
        dataset_directory,
        parameters,
        coverage_file=coverage_file,
    )

    assembly_candidates = load_assembly_edges(dataset)
    supported_assembly = filter_supported_assembly_edges(
        dataset.attributes,
        dataset.coverage,
        assembly_candidates,
        parameters,
    )
    training = train_mage_embedding(
        dataset,
        parameters,
        supported_assembly_edges=supported_assembly,
        assembly_candidate_count=len(assembly_candidates),
        device=device,
    )
    completeness_gate = None
    if parameters.isolate_complete_contigs:
        completeness_gate = detect_complete_contigs(
            dataset,
            parameters,
            cache_file=output / "terminal_repeats.tsv",
        )
    clustering = cluster_unknown_number_of_bins(
        training.embedding,
        dataset.coverage,
        parameters,
        supported_assembly_edges=supported_assembly,
        assembly_candidate_count=len(assembly_candidates),
        completeness_gate=completeness_gate,
    )

    labels = clustering.labels
    sizes = np.bincount(labels)
    assignment_file = output / "assignments.tsv"
    pd.DataFrame(
        {
            "contig_id": dataset.contig_ids,
            "bin_id": [f"MAGE-Bin_{label:06d}" for label in labels],
            "length": dataset.metadata["length"].astype(int),
        }
    ).to_csv(assignment_file, sep="\t", index=False)

    bin_fasta_directory = None
    source_fasta = dataset.manifest.get("combined_fasta")
    if source_fasta and Path(source_fasta).is_file():
        bin_fasta_directory = _write_bin_fastas(
            source_fasta, output, dataset.contig_ids, labels
        )

    metadata = {
        "tool": "MAGE-Bin",
        "dataset": dataset.name,
        "dataset_directory": dataset.directory,
        "supervised_labels_used": False,
        "checkv_used_for_selection": False,
        "parameters": parameters.to_dict(),
        "features": {
            "contigs": len(dataset),
            "dimensions": dataset.attributes.shape[1],
            "coverage_dimensions": dataset.coverage.shape[1],
            "feature_seconds": dataset.feature_seconds,
        },
        "assembly": {
            "candidate_edges": len(assembly_candidates),
            "supported_edges": len(supported_assembly),
        },
        "training": training.statistics,
        "training_graph": training.graph_statistics.to_dict(),
        "clustering": clustering.statistics,
        "inference_graph": clustering.graph_statistics.to_dict(),
        "completeness_gate": (
            {key: value for key, value in completeness_gate.items() if key != "locked"}
            if completeness_gate is not None
            else None
        ),
        "result": {
            "bins": len(sizes),
            "bin_fasta_directory": bin_fasta_directory,
            "singletons": int((sizes == 1).sum()),
            "multi_contig_bins": int((sizes > 1).sum()),
            "largest_bin_contigs": int(sizes.max()),
            "seconds": time.perf_counter() - began,
        },
    }
    run_metadata_file = output / "run.json"
    run_metadata_file.write_text(
        json.dumps(metadata, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    return BinningResult(
        assignment_file=assignment_file,
        run_metadata_file=run_metadata_file,
        contig_count=len(dataset),
        bin_count=len(sizes),
        singleton_count=int((sizes == 1).sum()),
        bin_fasta_directory=bin_fasta_directory,
    )
