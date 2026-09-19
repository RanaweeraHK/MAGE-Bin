"""Learnable-gate MAGE-Bin neural model and training loop."""

from __future__ import annotations

import copy
import math
import time
from dataclasses import dataclass

import numpy as np
import torch
from torch import nn
from torch.nn import functional

from .config import MageBinConfig
from .dataset import BinningDataset
from .evidence import normalize_rows
from .graph import GraphMergeStats, build_normalized_adjacency, propose_training_graph


class GraphMessageLayer(nn.Module):
    """Combine transformed self features with normalized neighbor messages."""

    def __init__(self, input_dimension: int, output_dimension: int) -> None:
        super().__init__()
        self.self_projection = nn.Linear(input_dimension, output_dimension)
        self.neighbor_projection = nn.Linear(
            input_dimension, output_dimension, bias=False
        )
        self.normalization = nn.BatchNorm1d(output_dimension)

    def forward(self, features: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        message = self.self_projection(features) + self.neighbor_projection(
            torch.sparse.mm(adjacency, features)
        )
        return functional.gelu(self.normalization(message))


class MAGEEncoder(nn.Module):
    """Fuse identity and graph representations through a learned convex gate."""

    def __init__(self, input_dimension: int, config: MageBinConfig) -> None:
        super().__init__()
        self.identity_projection = nn.Linear(
            input_dimension, config.latent_dimension, bias=False
        )
        self.first_graph_layer = GraphMessageLayer(
            input_dimension, config.hidden_dimension
        )
        self.second_graph_layer = GraphMessageLayer(
            config.hidden_dimension, config.latent_dimension
        )
        initial = config.initial_graph_weight
        self.graph_gate_logit = nn.Parameter(
            torch.tensor(
                math.log(initial) - math.log1p(-initial), dtype=torch.float32
            )
        )

    def fusion_weights(self) -> tuple[torch.Tensor, torch.Tensor]:
        """Return identity and graph weights that are positive and sum to one."""

        graph_weight = torch.sigmoid(self.graph_gate_logit)
        return 1.0 - graph_weight, graph_weight

    def forward(self, features: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        identity = self.identity_projection(features)
        graph = self.second_graph_layer(
            self.first_graph_layer(features, adjacency), adjacency
        )
        identity_weight, graph_weight = self.fusion_weights()
        return functional.normalize(
            identity_weight * identity + graph_weight * graph, dim=1
        )


class MaskedMAGEAutoencoder(nn.Module):
    """Reconstruct masked node features from the learned graph representation."""

    def __init__(self, input_dimension: int, config: MageBinConfig) -> None:
        super().__init__()
        self.mask_token = nn.Parameter(torch.zeros(input_dimension))
        self.encoder = MAGEEncoder(input_dimension, config)
        self.decoder = nn.Sequential(
            nn.Linear(config.latent_dimension, config.hidden_dimension),
            nn.GELU(),
            nn.Linear(config.hidden_dimension, input_dimension),
        )

    def encode(self, features: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        return self.encoder(features, adjacency)

    def forward(
        self,
        features: torch.Tensor,
        adjacency: torch.Tensor,
        masked_nodes: torch.Tensor,
    ) -> torch.Tensor:
        corrupted = features.clone()
        corrupted[masked_nodes] = self.mask_token
        return self.decoder(self.encode(corrupted, adjacency))


def masked_cosine_reconstruction_loss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    masked_nodes: torch.Tensor,
) -> torch.Tensor:
    """Return mean cosine reconstruction loss on the masked nodes."""

    if masked_nodes.numel() == 0:
        raise ValueError("at least one masked node is required")
    return 1.0 - functional.cosine_similarity(
        prediction[masked_nodes], target[masked_nodes]
    ).mean()


def normalized_effective_rank(embedding: np.ndarray) -> float:
    """Measure representation collapse; one denotes full effective rank."""

    sample = (
        embedding[np.linspace(0, len(embedding) - 1, 4_096, dtype=int)]
        if len(embedding) > 4_096
        else embedding
    )
    centered = sample - sample.mean(axis=0)
    eigenvalues = np.linalg.eigvalsh(centered.T @ centered)
    singular_values = np.sqrt(np.clip(eigenvalues, 0.0, None))
    probabilities = singular_values / max(singular_values.sum(), 1e-12)
    rank = math.exp(
        float(-(probabilities * np.log(probabilities + 1e-12)).sum())
    )
    return rank / max(min(sample.shape), 1)


@dataclass(slots=True)
class TrainingResult:
    """Embedding and diagnostics selected by held-out reconstruction."""

    embedding: np.ndarray
    statistics: dict[str, object]
    graph_statistics: GraphMergeStats


def train_mage_embedding(
    dataset: BinningDataset,
    config: MageBinConfig,
    *,
    supported_assembly_edges: np.ndarray | None = None,
    assembly_candidate_count: int = 0,
    device: str | torch.device = "auto",
) -> TrainingResult:
    """Train MAGE-Bin without using truth labels or external quality scores."""

    began = time.perf_counter()
    resolved_device = torch.device(
        "cuda" if device == "auto" and torch.cuda.is_available() else (
            "cpu" if device == "auto" else device
        )
    )
    attributes = dataset.attributes
    coverage = dataset.coverage
    node_count = len(attributes)
    if node_count < 2:
        raise ValueError("training requires at least two nodes")

    _reference, edges, pair_prior, degree_ceiling, graph_stats = propose_training_graph(
        attributes,
        coverage,
        config,
        supported_assembly_edges=supported_assembly_edges,
        assembly_candidate_count=assembly_candidate_count,
    )
    torch.manual_seed(config.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(config.seed)
    rng = np.random.default_rng(config.seed)
    features = torch.as_tensor(attributes, dtype=torch.float32, device=resolved_device)
    model = MaskedMAGEAutoencoder(attributes.shape[1], config).to(resolved_device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay,
    )

    validation_size = min(
        node_count - 1,
        max(1, round(node_count * config.validation_fraction)),
    )
    validation_nodes_array = rng.choice(
        node_count, validation_size, replace=False
    )
    training_pool = np.setdiff1d(np.arange(node_count), validation_nodes_array)
    validation_nodes = torch.as_tensor(
        validation_nodes_array, dtype=torch.long, device=resolved_device
    )

    best_loss = math.inf
    best_state: dict[str, torch.Tensor] | None = None
    best_edges = edges.copy()
    best_graph_stats = graph_stats
    best_epoch = 0
    stale_checks = 0
    previous_edges = set(map(tuple, edges.tolist()))
    edge_jaccard = 0.0
    effective_rank = 0.0
    stop_reason = "maximum_epochs"
    adjacency = build_normalized_adjacency(node_count, edges, resolved_device)

    final_epoch = 0
    for epoch in range(config.maximum_epochs):
        final_epoch = epoch + 1
        model.train()
        mask_size = min(
            len(training_pool),
            max(1, round(len(training_pool) * config.mask_fraction)),
        )
        masked_array = rng.choice(training_pool, mask_size, replace=False)
        masked_nodes = torch.as_tensor(
            masked_array, dtype=torch.long, device=resolved_device
        )
        prediction = model(features, adjacency, masked_nodes)
        loss = masked_cosine_reconstruction_loss(
            prediction, features, masked_nodes
        )
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if final_epoch % config.graph_update_interval:
            continue

        model.eval()
        with torch.no_grad():
            validation_prediction = model(features, adjacency, validation_nodes)
            validation_loss = float(
                masked_cosine_reconstruction_loss(
                    validation_prediction, features, validation_nodes
                ).cpu()
            )
            clean_embedding = model.encode(features, adjacency).cpu().numpy()

        (
            _reference,
            proposed_edges,
            pair_prior,
            degree_ceiling,
            proposed_stats,
        ) = propose_training_graph(
            attributes,
            coverage,
            config,
            embedding=clean_embedding,
            supported_assembly_edges=supported_assembly_edges,
            assembly_candidate_count=assembly_candidate_count,
        )
        proposed_set = set(map(tuple, proposed_edges.tolist()))
        edge_jaccard = len(previous_edges & proposed_set) / max(
            len(previous_edges | proposed_set), 1
        )
        effective_rank = normalized_effective_rank(clean_embedding)

        if validation_loss < best_loss - 1e-4:
            best_loss = validation_loss
            best_state = copy.deepcopy(model.state_dict())
            best_edges = proposed_edges.copy()
            best_graph_stats = proposed_stats
            best_epoch = final_epoch
            stale_checks = 0
        else:
            stale_checks += 1

        edges = proposed_edges
        graph_stats = proposed_stats
        previous_edges = proposed_set
        adjacency = build_normalized_adjacency(node_count, edges, resolved_device)
        if (
            final_epoch >= config.minimum_epochs
            and stale_checks >= config.validation_patience
            and effective_rank >= config.minimum_effective_rank
        ):
            stop_reason = (
                "validation_and_graph_converged"
                if edge_jaccard >= config.graph_jaccard_stop
                else "heldout_reconstruction_plateau"
            )
            break

    if best_state is not None:
        model.load_state_dict(best_state)
        edges = best_edges
        graph_stats = best_graph_stats
        adjacency = build_normalized_adjacency(node_count, edges, resolved_device)
    model.eval()
    with torch.no_grad():
        embedding = model.encode(features, adjacency).cpu().numpy()
        identity_weight, graph_weight = model.encoder.fusion_weights()
    degree = (
        np.bincount(edges.ravel(), minlength=node_count)
        if len(edges)
        else np.zeros(node_count)
    )
    return TrainingResult(
        embedding=normalize_rows(embedding),
        statistics={
            "device": str(resolved_device),
            "training_edges": len(edges),
            "graph_touched": int((degree > 0).sum()),
            "mean_degree": float(2 * len(edges) / max(node_count, 1)),
            "degree_ceiling": degree_ceiling,
            "selected_epoch": best_epoch or final_epoch,
            "stop_reason": stop_reason,
            "heldout_reconstruction": best_loss,
            "edge_jaccard": edge_jaccard,
            "effective_rank": effective_rank,
            "pair_prior": pair_prior,
            "identity_gate_weight": float(identity_weight.cpu()),
            "graph_gate_weight": float(graph_weight.cpu()),
            "training_seconds": time.perf_counter() - began,
        },
        graph_statistics=graph_stats,
    )
