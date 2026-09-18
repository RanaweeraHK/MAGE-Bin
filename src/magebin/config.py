"""Typed configuration for MAGE-Bin training and inference."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True, slots=True)
class MageBinConfig:
    """All reproducibility-sensitive MAGE-Bin parameters.

    Defaults reproduce the learnable-gate notebook experiment. The initial
    graph weight of ``1 / 11`` is exactly equivalent to the prototype's
    ``identity + 0.10 * graph`` ratio after convex normalization.
    """

    min_contig_length: int = 2_000
    seed: int = 0
    sparse_neighbors: int = 48
    exact_knn_limit: int = 4_000
    null_pair_samples: int = 600_000
    prior_refinement_rounds: int = 3
    prior_tolerance: float = 0.02

    latent_dimension: int = 64
    hidden_dimension: int = 128
    mask_fraction: float = 0.30
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    minimum_epochs: int = 40
    maximum_epochs: int = 100
    graph_update_interval: int = 10
    validation_fraction: float = 0.10
    validation_patience: int = 3
    graph_jaccard_stop: float = 0.90
    minimum_effective_rank: float = 0.20
    learned_view_weight: float = 0.15
    view_disagreement_penalty: float = 0.25
    initial_graph_neighbors: int = 8
    small_graph_resolution_scale: float = 1.05
    initial_graph_weight: float = 1.0 / 11.0

    isolate_complete_contigs: bool = True
    terminal_repeat_seed: int = 21
    terminal_repeat_minimum: int = 20
    terminal_repeat_max_fraction: float = 0.20

    def __post_init__(self) -> None:
        if self.min_contig_length < 1:
            raise ValueError("min_contig_length must be positive")
        if self.sparse_neighbors < 1 or self.initial_graph_neighbors < 1:
            raise ValueError("neighbor counts must be positive")
        if self.minimum_epochs < 1:
            raise ValueError("minimum_epochs must be positive")
        if self.maximum_epochs < self.minimum_epochs:
            raise ValueError("maximum_epochs must be >= minimum_epochs")
        if self.graph_update_interval < 1:
            raise ValueError("graph_update_interval must be positive")
        for name in ("mask_fraction", "validation_fraction", "initial_graph_weight"):
            value = getattr(self, name)
            if not 0.0 < value < 1.0:
                raise ValueError(f"{name} must lie strictly between 0 and 1")

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable parameter mapping."""

        return asdict(self)

