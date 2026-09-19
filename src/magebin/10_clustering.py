"""Unknown-K evidence-calibrated partitioning for MAGE-Bin."""

from __future__ import annotations

import math
from dataclasses import dataclass

import igraph as ig
import leidenalg
import numpy as np

from .completeness import isolate_complete_contigs
from .config import MageBinConfig
from .evidence import (
    collect_candidate_edges,
    compute_edge_evidence,
    estimate_positive_edge_prior,
    fit_edge_score_mixture,
    fit_robust_null,
    mixture_is_informative,
    sample_null_similarity_scores,
)
from .graph import (
    GraphMergeStats,
    limit_node_degree,
    merge_supported_assembly_edges,
)


@dataclass(slots=True)
class ClusteringResult:
    """Contig bin labels and diagnostics from unknown-K inference."""

    labels: np.ndarray
    statistics: dict[str, object]
    graph_statistics: GraphMergeStats


def estimate_partition_pair_prior(labels: np.ndarray) -> float:
    """Estimate the probability that two nodes share a predicted bin."""

    sizes = np.bincount(np.asarray(labels, dtype=np.int64))
    node_count = sizes.sum()
    pairs = (sizes * (sizes - 1) / 2).sum()
    possible = max(node_count * (node_count - 1) / 2, 1)
    return float(np.clip(pairs / possible, 1e-7, 0.5))


def partition_sparse_graph(
    node_count: int,
    pairs: np.ndarray,
    weights: np.ndarray,
    omitted_evidence: float,
    initial_labels: np.ndarray,
    *,
    coverage_active: bool,
    coverage_available: float,
    config: MageBinConfig,
) -> tuple[np.ndarray, float, float]:
    """Run sparse Leiden CPM with evidence-density-adaptive resolution."""

    correction = np.asarray(weights - omitted_evidence, dtype=np.float64)
    keep = np.isfinite(correction) & (correction > 0)
    if not np.any(keep):
        return np.arange(node_count), 0.0, 1.0
    average_degree = float(2 * keep.sum() / max(node_count, 1))
    if node_count <= config.exact_knn_limit:
        density_scale = config.small_graph_resolution_scale
    elif not coverage_active and coverage_available >= 0.98:
        density_scale = float(np.clip(average_degree / 10.0, 0.35, 0.85))
    elif not coverage_active:
        density_scale = 1.25
    else:
        density_scale = float(np.clip(average_degree / 8.0, 1.0, 1.75))
    density_scale *= 1.05
    graph = ig.Graph(n=node_count, edges=pairs[keep].tolist(), directed=False)
    partition = leidenalg.find_partition(
        graph,
        leidenalg.CPMVertexPartition,
        weights=correction[keep].tolist(),
        resolution_parameter=max(-omitted_evidence * density_scale, 1e-6),
        initial_membership=np.asarray(initial_labels, dtype=int).tolist(),
        n_iterations=3,
        seed=config.seed,
    )
    labels = np.unique(np.asarray(partition.membership), return_inverse=True)[1]
    return labels, average_degree, density_scale


def cluster_unknown_number_of_bins(
    embedding: np.ndarray,
    coverage: np.ndarray,
    config: MageBinConfig,
    *,
    supported_assembly_edges: np.ndarray | None = None,
    assembly_candidate_count: int = 0,
    completeness_gate: dict[str, object] | None = None,
) -> ClusteringResult:
    """Infer bins without labels, a specified K, or CheckV feedback."""

    node_count = len(embedding)
    assembly = (
        np.empty((0, 2), dtype=np.int64)
        if supported_assembly_edges is None
        else supported_assembly_edges
    )
    rng = np.random.default_rng(config.seed)
    embedding_null = sample_null_similarity_scores(
        embedding, rng, config.null_pair_samples
    )
    coverage_null = sample_null_similarity_scores(
        coverage, rng, config.null_pair_samples
    )
    coverage_null_valid = fit_robust_null(coverage_null)[1] >= 1e-3
    prior = estimate_positive_edge_prior(
        coverage_null if coverage_null_valid else embedding_null
    )

    statistical_pairs = collect_candidate_edges((embedding, coverage), config)
    statistical_pairs, degree_ceiling = limit_node_degree(
        embedding, coverage, statistical_pairs, config
    )
    pairs, graph_stats = merge_supported_assembly_edges(
        statistical_pairs,
        assembly,
        assembly_candidate_count=assembly_candidate_count,
    )
    withheld_edges = 0
    if completeness_gate is not None and len(pairs):
        free = ~np.asarray(completeness_gate["locked"], dtype=bool)
        retained = free[pairs[:, 0]] & free[pairs[:, 1]]
        withheld_edges = int((~retained).sum())
        pairs = pairs[retained]

    if not len(pairs):
        labels = np.arange(node_count)
        if completeness_gate is not None:
            labels = isolate_complete_contigs(labels, completeness_gate)
        return ClusteringResult(
            labels=labels,
            statistics={
                "pair_prior": prior,
                "rounds": 0,
                "inference_graph": "empty",
                "candidate_edges": 0,
                "degree_ceiling": degree_ceiling,
                "positive_average_degree": 0.0,
                "density_scale": 1.0,
                "omitted_evidence": math.nan,
                "active_views": "",
                "withheld_edges": withheld_edges,
                "locked_contigs": int(
                    completeness_gate["locked_contigs"]
                    if completeness_gate is not None
                    else 0
                ),
            },
            graph_statistics=graph_stats,
        )

    labels = np.arange(node_count)
    omitted_evidence = math.nan
    average_degree = 0.0
    density_scale = 1.0
    embedding_active = False
    coverage_active = False
    iteration = -1
    for iteration in range(config.prior_refinement_rounds):
        embedding_fit = fit_edge_score_mixture(embedding_null, prior)
        coverage_fit = fit_edge_score_mixture(coverage_null, prior)
        left, right = pairs.T
        embedding_evidence = compute_edge_evidence(
            (embedding[left] * embedding[right]).sum(1), embedding_fit, prior
        )
        coverage_evidence = compute_edge_evidence(
            (coverage[left] * coverage[right]).sum(1), coverage_fit, prior
        )
        embedding_active = mixture_is_informative(embedding_fit) and bool(
            np.any(embedding_evidence > 0)
        )
        coverage_active = mixture_is_informative(coverage_fit) and bool(
            np.any(coverage_evidence > 0)
        )
        if not embedding_active and not coverage_active:
            embedding_active = True
        weights = (
            embedding_evidence if embedding_active else 0.0
        ) + (coverage_evidence if coverage_active else 0.0)
        if embedding_active and coverage_active:
            weights -= config.view_disagreement_penalty * np.abs(
                embedding_evidence - coverage_evidence
            )
        null_terms = []
        if embedding_active:
            null_terms.append(
                np.median(compute_edge_evidence(embedding_null, embedding_fit, prior))
            )
        if coverage_active:
            null_terms.append(
                np.median(compute_edge_evidence(coverage_null, coverage_fit, prior))
            )
        omitted_evidence = float(min(-1e-6, sum(null_terms)))
        coverage_available = float(
            (np.linalg.norm(coverage, axis=1) > 1e-8).mean()
        )
        labels, average_degree, density_scale = partition_sparse_graph(
            node_count,
            pairs,
            weights,
            omitted_evidence,
            labels,
            coverage_active=coverage_active,
            coverage_available=coverage_available,
            config=config,
        )
        implied_prior = estimate_partition_pair_prior(labels)
        if abs(math.log(implied_prior) - math.log(prior)) < config.prior_tolerance:
            break
        prior = math.sqrt(prior * implied_prior)

    if completeness_gate is not None:
        labels = isolate_complete_contigs(labels, completeness_gate)
    return ClusteringResult(
        labels=labels,
        statistics={
            "locked_contigs": int(
                completeness_gate["locked_contigs"]
                if completeness_gate is not None
                else 0
            ),
            "withheld_edges": withheld_edges,
            "pair_prior": prior,
            "rounds": iteration + 1,
            "inference_graph": "dynamic-degree-controlled-knn-plus-assembly",
            "candidate_edges": len(pairs),
            "degree_ceiling": degree_ceiling,
            "positive_average_degree": average_degree,
            "density_scale": density_scale,
            "omitted_evidence": omitted_evidence,
            "active_views": "+".join(
                name
                for name, active in (
                    ("embedding", embedding_active),
                    ("coverage", coverage_active),
                )
                if active
            ),
        },
        graph_statistics=graph_stats,
    )

