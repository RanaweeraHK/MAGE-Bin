from dataclasses import replace

import pytest

from magebin.config import MageBinConfig


def test_default_gate_matches_original_ratio():
    assert MageBinConfig().initial_graph_weight == pytest.approx(1 / 11)


def test_configuration_rejects_invalid_epoch_order():
    with pytest.raises(ValueError, match="maximum_epochs"):
        replace(MageBinConfig(), minimum_epochs=10, maximum_epochs=5)

