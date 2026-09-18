"""Statistical and assembly-evidence graph construction."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import torch

from .config import MageBinConfig
from .dataset import BinningDataset
from .evidence import (
    build_positive_evidence_graph,
    collect_candidate_edges,
    compute_edge_evidence,
    estimate_positive_edge_prior,
    fit_edge_score_mixture,
    fit_robust_null,
    mixture_is_informative,
    normalize_rows,
    sample_null_similarity_scores,
)


@dataclass(frozen=True, slots=True)
class GraphMergeStats:
    """Counts describing statistical and assembly-supported graph fusion."""

    assembly_candidates: int
    assembly_supported: int
    statistical_edges: int
    assembly_edges_added: int
    evidence_graph_edges: int

    def to_dict(self) -> dict[str, int]:
        return {
            "assembly_candidates": self.assembly_candidates,
            "assembly_supported": self.assembly_supported,
            "statistical_edges": self.statistical_edges,
            "assembly_edges_added": self.assembly_edges_added,
            "evidence_graph_edges": self.evidence_graph_edges,
        }


def _canonical_edges(edges: np.ndarray) -> np.ndarray:
    array = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    if not len(array):
        return array
    array = np.sort(array, axis=1)
    array = array[array[:, 0] != array[:, 1]]
    if not len(array):
        return array.reshape(-1, 2)
    return np.unique(array, axis=0)


def load_assembly_edges(dataset: BinningDataset) -> np.ndarray:
    """Map assembly-graph links onto the retained MAGE-Bin contig indices."""

    graph_path = dataset.directory / "viral_graph.pt"
    if not graph_path.is_file():
        return np.empty((0, 2), dtype=np.int64)
    graph = torch.load(graph_path, map_location="cpu", weights_only=False)
    if not hasattr(graph, "contig_names") or not hasattr(graph, "edge_index"):
        raise ValueError("viral_graph.pt must provide contig_names and edge_index")

    graph_names = [str(name) for name in graph.contig_names]
    model_index = {name: index for index, name in enumerate(dataset.contig_ids)}
    edge_index = graph.edge_index.detach().cpu().numpy().T
    if edge_index.ndim != 2 or edge_index.shape[1] != 2:
        raise ValueError("assembly edge_index must have shape [2, edge_count]")

    edge_attributes = getattr(graph, "edge_attr", None)
    edge_types = getattr(graph, "edge_type", None)
    if edge_attributes is not None:
        values = edge_attributes.detach().cpu().numpy()
        if values.ndim == 1:
            is_assembly_link = values > 0
        else:
            is_assembly_link = np.any(values[:, : min(2, values.shape[1])] > 0, axis=1)
    elif edge_types is not None:
        values = edge_types.detach().cpu().numpy().reshape(-1).astype(np.int64)
        is_assembly_link = (values & 3) != 0
    else:
        is_assembly_link = np.ones(len(edge_index), dtype=bool)
    if len(is_assembly_link) != len(edge_index):
        raise ValueError("assembly edge attributes do not match edge_index")

    pairs: list[tuple[int, int]] = []
    for source, target in edge_index[is_assembly_link]:
        if not (0 <= source < len(graph_names) and 0 <= target < len(graph_names)):
            raise ValueError("assembly edge_index contains an invalid node index")
        left = model_index.get(graph_names[int(source)])
        right = model_index.get(graph_names[int(target)])
        if left is not None and right is not None and left != right:
            pairs.append((left, right))
    return _canonical_edges(np.asarray(pairs, dtype=np.int64).reshape(-1, 2))


def filter_supported_assembly_edges(
    attributes: np.ndarray,
    coverage: np.ndarray,
    assembly_edges: np.ndarray,
    config: MageBinConfig,
) -> np.ndarray:
    """Retain assembly links with positive sequence/coverage evidence.

    Assembly connectivity is treated as candidate biological evidence, not as
    an unconditional same-genome label.
    """

    assembly_edges = _canonical_edges(assembly_edges)
    if not len(assembly_edges):
        return assembly_edges
    rng = np.random.default_rng(config.seed)
    attribute_null = sample_null_similarity_scores(
        attributes, rng, config.null_pair_samples
    )
    coverage_null = sample_null_similarity_scores(
        coverage, rng, config.null_pair_samples
    )
    coverage_valid = fit_robust_null(coverage_null)[1] >= 1e-3
    prior = estimate_positive_edge_prior(
        coverage_null if coverage_valid else attribute_null
    )
    attribute_fit = fit_edge_score_mixture(attribute_null, prior)
    coverage_fit = fit_edge_score_mixture(coverage_null, prior)
    left, right = assembly_edges.T
    evidence = np.zeros(len(assembly_edges), dtype=np.float64)
    active = False
    if mixture_is_informative(attribute_fit):
        evidence += compute_edge_evidence(
            (attributes[left] * attributes[right]).sum(1), attribute_fit, prior
        )
        active = True
    if mixture_is_informative(coverage_fit):
        evidence += compute_edge_evidence(
            (coverage[left] * coverage[right]).sum(1), coverage_fit, prior
        )
        active = True
    keep = active & np.isfinite(evidence) & (evidence > 0.0)
    return assembly_edges[keep]


def limit_node_degree(
    reference: np.ndarray,
    coverage: np.ndarray,
    edges: np.ndarray,
    config: MageBinConfig,
) -> tuple[np.ndarray, int]:
    """Keep each endpoint's strongest edges under a dataset-sized ceiling."""

    edges = _canonical_edges(edges)
    node_count = len(reference)
    ceiling = max(8, min(20, math.ceil(math.log2(max(node_count, 2)))))
    if (
        node_count <= config.exact_knn_limit
        or not len(edges)
        or len(edges) <= node_count * ceiling // 2
    ):
        return edges, ceiling

    left, right = edges.T
    reference_score = (reference[left] * reference[right]).sum(1)
    coverage_score = (coverage[left] * coverage[right]).sum(1)
    coverage_nonzero = np.linalg.norm(coverage, axis=1) > 1e-8
    scores = reference_score + coverage_score * (
        coverage_nonzero[left] & coverage_nonzero[right]
    )
    incident: list[list[int]] = [[] for _ in range(node_count)]
    for edge_id, (left_node, right_node) in enumerate(edges):
        incident[int(left_node)].append(edge_id)
        incident[int(right_node)].append(edge_id)
    selected = [
        set(sorted(row, key=lambda edge_id: -scores[edge_id])[:ceiling])
        for row in incident
    ]
    keep = np.fromiter(
        (
            edge_id in selected[int(left_node)]
            and edge_id in selected[int(right_node)]
            for edge_id, (left_node, right_node) in enumerate(edges)
        ),
        dtype=bool,
        count=len(edges),
    )
    return edges[keep], ceiling


def merge_supported_assembly_edges(
    statistical_edges: np.ndarray,
    supported_assembly_edges: np.ndarray,
    *,
    assembly_candidate_count: int,
) -> tuple[np.ndarray, GraphMergeStats]:
    """Union supported assembly links with an unchanged statistical edge set."""

    statistical = {tuple(edge) for edge in _canonical_edges(statistical_edges).tolist()}
    assembly = {
        tuple(edge) for edge in _canonical_edges(supported_assembly_edges).tolist()
    }
    merged = np.asarray(sorted(statistical | assembly), dtype=np.int64).reshape(-1, 2)
    return merged, GraphMergeStats(
        assembly_candidates=assembly_candidate_count,
        assembly_supported=len(assembly),
        statistical_edges=len(statistical),
        assembly_edges_added=len(assembly - statistical),
        evidence_graph_edges=len(merged),
    )


def propose_training_graph(
    attributes: np.ndarray,
    coverage: np.ndarray,
    config: MageBinConfig,
    *,
    embedding: np.ndarray | None = None,
    supported_assembly_edges: np.ndarray | None = None,
    assembly_candidate_count: int = 0,
) -> tuple[np.ndarray, np.ndarray, float, int, GraphMergeStats]:
    """Construct a seed or learned graph and add supported assembly evidence."""

    assembly = (
        np.empty((0, 2), dtype=np.int64)
        if supported_assembly_edges is None
        else supported_assembly_edges
    )
    if embedding is None:
        reference = attributes
        statistical_edges = collect_candidate_edges(
            (reference, coverage),
            config,
            neighbors=config.initial_graph_neighbors,
        )
        node_count = len(reference)
        prior = float(
            np.clip(
                2.0 * len(statistical_edges) / max(node_count * (node_count - 1), 1),
                1e-9,
                0.5,
            )
        )
    else:
        reference = (
            attributes
            if len(attributes) <= config.exact_knn_limit
            else normalize_rows(
                np.concatenate(
                    (attributes, config.learned_view_weight * embedding), axis=1
                )
            )
        )
        statistical_edges, prior = build_positive_evidence_graph(
            reference, coverage, config
        )
    statistical_edges, ceiling = limit_node_degree(
        reference, coverage, statistical_edges, config
    )
    merged, merge_stats = merge_supported_assembly_edges(
        statistical_edges,
        assembly,
        assembly_candidate_count=assembly_candidate_count,
    )
    return reference, merged, prior, ceiling, merge_stats


def build_normalized_adjacency(
    node_count: int,
    edges: np.ndarray,
    device: torch.device,
    *,
    drop_rate: float = 0.0,
    rng: np.random.Generator | None = None,
) -> torch.Tensor:
    """Build symmetric D^-1/2 A D^-1/2 adjacency with self loops."""

    edges = _canonical_edges(edges)
    if drop_rate:
        if rng is None:
            raise ValueError("rng is required when drop_rate is non-zero")
        edges = edges[rng.random(len(edges)) >= drop_rate]
    if len(edges):
        row = np.r_[edges[:, 0], edges[:, 1], np.arange(node_count)]
        column = np.r_[edges[:, 1], edges[:, 0], np.arange(node_count)]
    else:
        row = column = np.arange(node_count)
    degree = np.bincount(row, minlength=node_count)
    values = 1.0 / np.sqrt(np.maximum(degree[row] * degree[column], 1))
    indices = torch.as_tensor(np.stack((row, column)), dtype=torch.long, device=device)
    weights = torch.as_tensor(values, dtype=torch.float32, device=device)
    return torch.sparse_coo_tensor(
        indices, weights, (node_count, node_count)
    ).coalesce()
