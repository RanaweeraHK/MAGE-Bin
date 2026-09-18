import numpy as np
import pytest

from magebin.metrics import (
    compute_binning_precision_recall_f1,
    partitions_are_equivalent,
)


def test_partition_equivalence_ignores_label_names():
    assert partitions_are_equivalent(
        np.array([0, 0, 1, 2, 2]), np.array(["a", "a", "z", "q", "q"])
    )
    assert not partitions_are_equivalent(
        np.array([0, 0, 1]), np.array(["a", "b", "b"])
    )


def test_perfect_binning_metrics():
    precision, recall, f1 = compute_binning_precision_recall_f1(
        np.array([4, 4, 8, 8]), np.array(["a", "a", "b", "b"])
    )
    assert precision == pytest.approx(1.0)
    assert recall == pytest.approx(1.0)
    assert f1 == pytest.approx(1.0)

