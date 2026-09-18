import numpy as np
import pytest
import torch

from magebin.config import MageBinConfig
from magebin.graph import build_normalized_adjacency
from magebin.model import MAGEEncoder, masked_cosine_reconstruction_loss


def test_encoder_gate_is_convex_and_uses_requested_initial_value():
    config = MageBinConfig(latent_dimension=3, hidden_dimension=4)
    encoder = MAGEEncoder(2, config)
    identity, graph = encoder.fusion_weights()
    assert float(identity + graph) == pytest.approx(1.0)
    assert float(graph) == pytest.approx(config.initial_graph_weight)


def test_encoder_output_is_normalized_and_gate_receives_gradient():
    torch.manual_seed(0)
    config = MageBinConfig(latent_dimension=3, hidden_dimension=4)
    encoder = MAGEEncoder(2, config)
    encoder.eval()
    features = torch.tensor([[1.0, 0.0], [0.8, 0.2], [0.0, 1.0]])
    adjacency = build_normalized_adjacency(
        3, np.array([[0, 1], [1, 2]]), torch.device("cpu")
    )
    output = encoder(features, adjacency)
    assert output.shape == (3, 3)
    assert torch.allclose(torch.linalg.norm(output, dim=1), torch.ones(3), atol=1e-5)
    output.sum().backward()
    assert encoder.graph_gate_logit.grad is not None


def test_masked_loss_requires_nodes():
    values = torch.eye(2)
    with pytest.raises(ValueError, match="masked node"):
        masked_cosine_reconstruction_loss(
            values, values, torch.tensor([], dtype=torch.long)
        )

