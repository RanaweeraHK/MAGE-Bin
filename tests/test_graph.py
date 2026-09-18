from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch

from magebin.dataset import BinningDataset
from magebin.graph import (
    build_normalized_adjacency,
    load_assembly_edges,
    merge_supported_assembly_edges,
)


def test_assembly_union_preserves_statistical_edges():
    statistical = np.array([[0, 1], [1, 2]])
    assembly = np.array([[1, 2], [2, 3]])
    merged, stats = merge_supported_assembly_edges(
        statistical, assembly, assembly_candidate_count=3
    )
    assert {tuple(edge) for edge in merged.tolist()} == {(0, 1), (1, 2), (2, 3)}
    assert stats.statistical_edges == 2
    assert stats.assembly_supported == 2
    assert stats.assembly_edges_added == 1


def test_normalized_adjacency_is_symmetric_with_self_loops():
    adjacency = build_normalized_adjacency(
        3, np.array([[0, 1], [1, 2]]), torch.device("cpu")
    ).to_dense()
    assert torch.allclose(adjacency, adjacency.T)
    assert torch.all(torch.diag(adjacency) > 0)


def test_load_assembly_edges_maps_and_deduplicates(tmp_path):
    graph = SimpleNamespace(
        contig_names=["a", "b", "filtered"],
        edge_index=torch.tensor([[0, 1, 1, 0], [1, 0, 2, 1]]),
        edge_attr=torch.ones(4, 1),
    )
    torch.save(graph, tmp_path / "viral_graph.pt")
    dataset = BinningDataset(
        name="tiny",
        directory=tmp_path,
        manifest={},
        metadata=pd.DataFrame({"contig_id": ["a", "b"], "length": [2000, 2000]}),
        attributes=np.eye(2, dtype=np.float32),
        coverage=np.eye(2, dtype=np.float32),
        strain_labels=None,
        species_labels=None,
        feature_seconds=0.0,
    )
    assert load_assembly_edges(dataset).tolist() == [[0, 1]]
