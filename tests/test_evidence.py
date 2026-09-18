import numpy as np
import pytest

from magebin.config import MageBinConfig
from magebin.evidence import (
    collect_candidate_edges,
    cosine_similarity_to_logit,
    normalize_rows,
)


def test_normalize_rows_is_finite_and_unit_length():
    values = normalize_rows(np.array([[3.0, 4.0], [0.0, 0.0], [np.nan, 2.0]]))
    assert np.isfinite(values).all()
    assert np.linalg.norm(values[0]) == pytest.approx(1.0)
    assert np.linalg.norm(values[1]) == pytest.approx(0.0)
    assert np.linalg.norm(values[2]) == pytest.approx(1.0)


def test_cosine_logit_is_monotonic_and_finite():
    transformed = cosine_similarity_to_logit(np.array([-1.0, 0.0, 1.0]))
    assert np.isfinite(transformed).all()
    assert np.all(np.diff(transformed) > 0)


def test_candidate_edges_have_no_duplicates_or_self_edges():
    view = np.array(
        [[1.0, 0.0], [0.99, 0.01], [0.0, 1.0], [0.01, 0.99]],
        dtype=float,
    )
    config = MageBinConfig(sparse_neighbors=2)
    edges = collect_candidate_edges((view, view), config)
    assert edges.shape[1] == 2
    assert np.all(edges[:, 0] < edges[:, 1])
    assert len(edges) == len({tuple(edge) for edge in edges.tolist()})

