"""Evaluation helpers kept separate from label-free MAGE-Bin inference."""

from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd


def compute_binning_precision_recall_f1(
    labels: np.ndarray,
    truth: np.ndarray,
) -> tuple[float, float, float]:
    """Compute contig-count binning precision, recall, and harmonic F1."""

    predicted = np.asarray(labels)
    actual = np.asarray(truth)
    if len(predicted) != len(actual):
        raise ValueError("labels and truth must have equal length")
    keep = pd.notna(actual)
    predicted = predicted[keep]
    actual = actual[keep]
    if not len(actual):
        return math.nan, math.nan, math.nan
    by_bin: defaultdict[object, Counter] = defaultdict(Counter)
    by_truth: defaultdict[object, Counter] = defaultdict(Counter)
    for predicted_label, actual_label in zip(predicted, actual, strict=True):
        by_bin[predicted_label][actual_label] += 1
        by_truth[actual_label][predicted_label] += 1
    precision = sum(max(counts.values()) for counts in by_bin.values()) / len(actual)
    recall = sum(max(counts.values()) for counts in by_truth.values()) / len(actual)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    return precision, recall, f1


def partitions_are_equivalent(
    left_labels: np.ndarray,
    right_labels: np.ndarray,
) -> bool:
    """Compare partitions while ignoring arbitrary cluster-label names."""

    left = np.asarray(left_labels)
    right = np.asarray(right_labels)
    if len(left) != len(right):
        return False
    left_to_right: defaultdict[object, set] = defaultdict(set)
    right_to_left: defaultdict[object, set] = defaultdict(set)
    for left_label, right_label in zip(left, right, strict=True):
        left_to_right[left_label].add(right_label)
        right_to_left[right_label].add(left_label)
    return all(len(values) == 1 for values in left_to_right.values()) and all(
        len(values) == 1 for values in right_to_left.values()
    )
