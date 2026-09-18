import json
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch

from magebin.config import MageBinConfig
from magebin.pipeline import run_binning


def _write_tiny_dataset(root):
    dataset = root / "tiny"
    work = root / "work"
    dataset.mkdir()
    work.mkdir()
    names = np.array([f"contig_{index}" for index in range(6)])
    pd.DataFrame(
        {
            "node_index": np.arange(6),
            "contig_id": names,
            "length": np.full(6, 2_500),
        }
    ).to_csv(dataset / "contig_metadata.tsv", sep="\t", index=False)
    counts = np.array(
        [
            [50, 45, 1, 1],
            [48, 44, 1, 1],
            [46, 42, 1, 1],
            [1, 1, 45, 50],
            [1, 1, 44, 48],
            [1, 1, 42, 46],
        ],
        dtype=float,
    )
    np.savez(dataset / "tnf_counts.npz", names=names, counts=counts)
    pd.DataFrame(
        {
            "contig": names,
            "sample_a": [100, 95, 90, 2, 2, 1],
            "sample_b": [1, 2, 2, 90, 95, 100],
        }
    ).to_csv(work / "coverage.csv", index=False)
    (dataset / "config.json").write_text(json.dumps({"work_dir": str(work)}))
    graph = SimpleNamespace(
        contig_names=names.tolist(),
        edge_index=torch.tensor([[0, 1, 3, 4], [1, 2, 4, 5]]),
        edge_attr=torch.ones(4, 1),
    )
    torch.save(graph, dataset / "viral_graph.pt")
    return dataset


def test_tiny_end_to_end_run(tmp_path):
    dataset = _write_tiny_dataset(tmp_path)
    output = tmp_path / "output"
    config = replace(
        MageBinConfig(),
        minimum_epochs=2,
        maximum_epochs=2,
        graph_update_interval=1,
        null_pair_samples=2_000,
        latent_dimension=4,
        hidden_dimension=8,
        isolate_complete_contigs=False,
    )
    result = run_binning(dataset, output, config=config, device="cpu")
    assignments = pd.read_csv(result.assignment_file, sep="\t")
    metadata = json.loads(result.run_metadata_file.read_text())
    assert len(assignments) == 6
    assert assignments.columns.tolist() == ["contig_id", "bin_id", "length"]
    assert assignments["bin_id"].str.startswith("MAGE-Bin_").all()
    assert metadata["supervised_labels_used"] is False
    assert metadata["checkv_used_for_selection"] is False
    assert 0 < metadata["training"]["graph_gate_weight"] < 1

