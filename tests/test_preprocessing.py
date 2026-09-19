"""Integration checks for the user-facing FASTA/coverage/GFA path."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from magebin.cli import main
from magebin.completeness import read_fasta
from magebin.config import MageBinConfig
from magebin.dataset import load_binning_dataset
from magebin.graph import load_assembly_edges
from magebin.preprocessing import prepare_dataset


def _inputs(tmp_path):
    fasta = tmp_path / "contigs.fasta"
    names = [f"contig_{index}" for index in range(6)]
    sequences = ["ACGT" * 625] * 3 + ["AGCT" * 625] * 3
    fasta.write_text(
        "".join(
            f">{name} example\n{sequence}\n"
            for name, sequence in zip(names, sequences, strict=True)
        )
    )
    coverage = tmp_path / "coverage.tsv"
    pd.DataFrame(
        {
            "contig": names,
            "sample_a": [100, 95, 90, 2, 2, 1],
            "sample_b": [1, 2, 2, 90, 95, 100],
        }
    ).to_csv(coverage, sep="\t", index=False)
    graph = tmp_path / "assembly_graph.gfa"
    graph.write_text(
        "H\tVN:Z:1.0\n"
        + "".join(f"S\t{name}\t*\tLN:i:2500\n" for name in names)
        + "L\tcontig_0\t+\tcontig_1\t+\t0M\n"
        + "L\tcontig_3\t+\tcontig_4\t+\t0M\n"
    )
    return fasta, coverage, graph


def test_run_from_standard_inputs(tmp_path):
    fasta, coverage, graph = _inputs(tmp_path)
    output = tmp_path / "results"
    assert (
        main(
            [
                "run",
                "--contigs",
                str(fasta),
                "--coverage",
                str(coverage),
                "--graph",
                str(graph),
                "--output",
                str(output),
                "--device",
                "cpu",
                "--minimum-epochs",
                "2",
                "--maximum-epochs",
                "2",
                "--no-complete-contig-gate",
            ]
        )
        == 0
    )
    dataset = load_binning_dataset(output / "preprocessed", MageBinConfig())
    assert dataset.attributes.shape[0] == 6
    with np.load(output / "preprocessed" / "tnf_counts.npz") as saved:
        assert saved["counts"].shape == (6, 136)
    assert load_assembly_edges(dataset).tolist() == [[0, 1], [3, 4]]
    assignments = pd.read_csv(output / "assignments.tsv", sep="\t")
    assert assignments["contig_id"].tolist() == dataset.contig_ids
    expected_bins = dict(zip(assignments.contig_id, assignments.bin_id, strict=True))
    bin_fastas = sorted((output / "bins").glob("*.fasta"))
    assert len(bin_fastas) == assignments.bin_id.nunique()
    for bin_fasta in bin_fastas:
        for contig_id, sequence in read_fasta(bin_fasta):
            assert expected_bins[contig_id] == bin_fasta.stem
            assert len(sequence) == 2500
    assert (
        json.loads((output / "run.json").read_text())["assembly"]["candidate_edges"]
        == 2
    )


def test_gfa_path_projection_and_unmapped_error(tmp_path):
    fasta, coverage, graph = _inputs(tmp_path)
    graph.write_text(
        "S\ta\t*\nS\tb\t*\nP\tcontig_0\ta+\t*\nP\tcontig_1\tb+\t*\nL\ta\t+\tb\t+\t0M\n"
    )
    dataset_dir = prepare_dataset(
        fasta, coverage, graph, tmp_path / "mapped", min_contig_length=2000
    )
    dataset = load_binning_dataset(dataset_dir, MageBinConfig())
    assert load_assembly_edges(dataset).tolist() == [[0, 1]]

    graph.write_text("S\ta\t*\nS\tb\t*\nL\ta\t+\tb\t+\t0M\n")
    with pytest.raises(ValueError, match="do not map"):
        prepare_dataset(
            fasta, coverage, graph, tmp_path / "unmapped", min_contig_length=2000
        )

    (tmp_path / "contigs.paths").write_text("contig_0\na+\ncontig_1\nb+\n")
    dataset_dir = prepare_dataset(
        fasta, coverage, graph, tmp_path / "spades", min_contig_length=2000
    )
    dataset = load_binning_dataset(dataset_dir, MageBinConfig())
    assert load_assembly_edges(dataset).tolist() == [[0, 1]]


def test_missing_coverage_row_is_rejected(tmp_path):
    fasta, coverage, graph = _inputs(tmp_path)
    frame = pd.read_csv(coverage, sep="\t").iloc[:-1]
    frame.to_csv(coverage, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing 1 contigs"):
        prepare_dataset(
            fasta, coverage, graph, tmp_path / "invalid", min_contig_length=2000
        )


def test_run_requires_usable_gfa(tmp_path):
    fasta, coverage, graph = _inputs(tmp_path)
    with pytest.raises(SystemExit) as exit_info:
        main(
            [
                "run",
                "--contigs",
                str(fasta),
                "--coverage",
                str(coverage),
                "--output",
                str(tmp_path / "missing-graph"),
            ]
        )
    assert exit_info.value.code == 2

    graph.write_text("S\tcontig_0\t*\nS\tcontig_1\t*\n")
    with pytest.raises(ValueError, match="do not map"):
        prepare_dataset(
            fasta, coverage, graph, tmp_path / "empty-graph", min_contig_length=2000
        )
