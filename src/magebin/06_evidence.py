"""Empirical-null calibration and sparse candidate-edge construction."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
from scipy.stats import norm
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors

from .config import MageBinConfig

MixtureFit = tuple[float, float, float, float]


def normalize_rows(values: np.ndarray) -> np.ndarray:
    """Return finite, unit-length rows as float32 values."""

    matrix = np.nan_to_num(np.asarray(values, dtype=np.float64))
    if matrix.ndim != 2:
        raise ValueError("expected a two-dimensional matrix")
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return (matrix / np.maximum(norms, 1e-12)).astype(np.float32)


def cosine_similarity_to_logit(cosine: np.ndarray) -> np.ndarray:
    """Map cosine similarity from [-1, 1] to stable log-odds."""

    probability = np.clip((np.asarray(cosine, float) + 1.0) * 0.5, 1e-4, 1 - 1e-4)
    return np.log(probability) - np.log1p(-probability)


def log_normal_density(values: np.ndarray, mean: float, scale: float) -> np.ndarray:
    """Evaluate a univariate normal log-density."""

    return (
        -0.5 * ((values - mean) / scale) ** 2
        - math.log(scale)
        - 0.5 * math.log(2 * math.pi)
    )


def fit_robust_null(scores: np.ndarray) -> tuple[float, float]:
    """Estimate null location and scale using median and MAD."""

    scores = np.asarray(scores, dtype=float)
    if scores.size == 0:
        return 0.0, 1e-9
    center = float(np.median(scores))
    spread = float(1.4826 * np.median(np.abs(scores - center)) + 1e-9)
    return center, spread


def estimate_positive_edge_prior(scores: np.ndarray) -> float:
    """Estimate the excess upper-tail mass over a robust Gaussian null."""

    scores = np.asarray(scores, dtype=float)
    if scores.size == 0:
        return 1e-9
    center, spread = fit_robust_null(scores)
    quantiles = np.linspace(0.5, 1 - 1e-5, min(2_000, max(100, scores.size)))
    thresholds = np.quantile(scores, quantiles)
    excess = ((1 - quantiles) - norm.sf((thresholds - center) / spread)).max()
    return float(np.clip(excess, 1e-9, 0.5))


def fit_edge_score_mixture(scores: np.ndarray, prior: float) -> MixtureFit:
    """Fit a fixed-null two-component Gaussian mixture to similarity scores."""

    scores = np.asarray(scores, dtype=float)
    null_mean, null_scale = fit_robust_null(scores)
    alternative_mean = null_mean + 2 * null_scale
    alternative_scale = null_scale
    if scores.size == 0:
        return null_mean, null_scale, alternative_mean, alternative_scale
    prior = float(np.clip(prior, 1e-9, 0.5))
    for _ in range(200):
        positive = math.log(prior) + log_normal_density(
            scores, alternative_mean, alternative_scale
        )
        negative = math.log1p(-prior) + log_normal_density(
            scores, null_mean, null_scale
        )
        top = np.maximum(positive, negative)
        responsibility = np.exp(positive - top) / (
            np.exp(positive - top) + np.exp(negative - top)
        )
        mass = responsibility.sum()
        if mass < 10:
            break
        new_mean = max(float((responsibility * scores).sum() / mass), null_mean)
        new_scale = math.sqrt(
            max(
                float((responsibility * (scores - new_mean) ** 2).sum() / mass),
                1e-6,
            )
        )
        if abs(new_mean - alternative_mean) + abs(new_scale - alternative_scale) < 1e-9:
            break
        alternative_mean, alternative_scale = new_mean, new_scale
    return null_mean, null_scale, alternative_mean, alternative_scale


def mixture_is_informative(mixture: MixtureFit) -> bool:
    """Return whether a fitted view has finite, non-degenerate parameters."""

    return bool(np.all(np.isfinite(mixture)) and mixture[1] >= 1e-3)


def compute_edge_evidence(
    cosine: np.ndarray,
    mixture: MixtureFit,
    prior: float,
) -> np.ndarray:
    """Compute posterior log-odds that candidate pairs are true edges."""

    null_mean, null_scale, alternative_mean, alternative_scale = mixture
    score = cosine_similarity_to_logit(cosine)
    return (
        math.log(prior)
        - math.log1p(-prior)
        + log_normal_density(score, alternative_mean, alternative_scale)
        - log_normal_density(score, null_mean, null_scale)
    )


def sample_null_similarity_scores(
    view: np.ndarray,
    rng: np.random.Generator,
    maximum_samples: int,
) -> np.ndarray:
    """Sample non-self pair similarities for empirical-null estimation."""

    matrix = normalize_rows(view)
    node_count = len(matrix)
    if node_count < 2:
        return np.empty(0, dtype=float)
    sample_count = min(maximum_samples, max(1_000, node_count * 200))
    left = rng.integers(0, node_count, sample_count)
    right = rng.integers(0, node_count, sample_count)
    keep = left != right
    return cosine_similarity_to_logit((matrix[left[keep]] * matrix[right[keep]]).sum(1))


def find_reciprocal_knn_pairs(
    view: np.ndarray,
    *,
    neighbors: int,
    exact_limit: int,
) -> set[tuple[int, int]]:
    """Find reciprocal cosine-nearest-neighbor pairs without self edges."""

    matrix = normalize_rows(view)
    node_count = len(matrix)
    if node_count < 2:
        return set()
    neighbors = min(neighbors, node_count - 1)
    pairs: set[tuple[int, int]] = set()

    def add_bucket(nodes: np.ndarray) -> None:
        nodes = np.asarray(nodes, dtype=np.int64)
        if len(nodes) < 2:
            return
        query_count = min(neighbors + 1, len(nodes))
        nearest = NearestNeighbors(
            n_neighbors=query_count,
            metric="cosine",
            algorithm="brute",
        ).fit(matrix[nodes]).kneighbors(matrix[nodes], return_distance=False)
        memberships = [set(nodes[row].tolist()) for row in nearest]
        positions = {int(node): local for local, node in enumerate(nodes)}
        for local, row in enumerate(nearest):
            left = int(nodes[local])
            for right_value in nodes[row]:
                right = int(right_value)
                if left != right and left in memberships[positions[right]]:
                    pairs.add((min(left, right), max(left, right)))

    if node_count <= exact_limit:
        add_bucket(np.arange(node_count))
    else:
        anchor_count = min(512, max(32, math.ceil(math.sqrt(node_count))))
        for seed in (0, 1):
            groups = MiniBatchKMeans(
                n_clusters=anchor_count,
                batch_size=4_096,
                n_init=1,
                max_iter=60,
                random_state=seed,
            ).fit_predict(matrix)
            for anchor in range(anchor_count):
                add_bucket(np.flatnonzero(groups == anchor))
    return pairs


def collect_candidate_edges(
    views: Iterable[np.ndarray],
    config: MageBinConfig,
    *,
    neighbors: int | None = None,
) -> np.ndarray:
    """Return the sorted union of reciprocal-neighbor edges from all views."""

    pairs: set[tuple[int, int]] = set()
    for view in views:
        pairs.update(
            find_reciprocal_knn_pairs(
                view,
                neighbors=neighbors or config.sparse_neighbors,
                exact_limit=config.exact_knn_limit,
            )
        )
    return np.asarray(sorted(pairs), dtype=np.int64).reshape(-1, 2)


def build_positive_evidence_graph(
    attributes: np.ndarray,
    coverage: np.ndarray,
    config: MageBinConfig,
) -> tuple[np.ndarray, float]:
    """Keep candidate edges with positive combined calibrated evidence."""

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
    candidates = collect_candidate_edges((attributes, coverage), config)
    if not len(candidates):
        return candidates, prior
    left, right = candidates.T
    attribute_evidence = compute_edge_evidence(
        (attributes[left] * attributes[right]).sum(1), attribute_fit, prior
    )
    coverage_evidence = compute_edge_evidence(
        (coverage[left] * coverage[right]).sum(1), coverage_fit, prior
    )
    attribute_active = mixture_is_informative(attribute_fit) and bool(
        np.any(attribute_evidence > 0)
    )
    coverage_active = mixture_is_informative(coverage_fit) and bool(
        np.any(coverage_evidence > 0)
    )
    if not attribute_active and not coverage_active:
        return candidates[:0], prior
    evidence = (
        attribute_evidence if attribute_active else 0.0
    ) + (coverage_evidence if coverage_active else 0.0)
    return candidates[evidence > 0], prior
