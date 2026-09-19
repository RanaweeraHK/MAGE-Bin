"""Public package interface for MAGE-Bin."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .config import MageBinConfig
    from .pipeline import BinningResult, run_binning

__all__ = ["BinningResult", "MageBinConfig", "run_binning"]
try:
    __version__ = version("magebin")
except PackageNotFoundError:
    __version__ = "unknown"


def __getattr__(name: str):
    """Load the training stack only when a public runtime object is requested."""

    if name == "MageBinConfig":
        from .config import MageBinConfig

        return MageBinConfig
    if name in {"BinningResult", "run_binning"}:
        from .pipeline import BinningResult, run_binning

        return {"BinningResult": BinningResult, "run_binning": run_binning}[name]
    raise AttributeError(name)
