"""Public package interface for MAGE-Bin."""

from __future__ import annotations

import sys
from importlib import import_module
from importlib.machinery import ModuleSpec
from importlib.metadata import PackageNotFoundError, version
from types import ModuleType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .config import MageBinConfig
    from .pipeline import BinningResult, run_binning

__all__ = ["BinningResult", "MageBinConfig", "run_binning"]


class _ModuleAlias(ModuleType):
    """Load a numbered module when an established import path is used."""

    def __init__(self, name: str, target: str) -> None:
        super().__init__(name)
        self.__package__ = __name__
        self.__spec__ = ModuleSpec(name, loader=None)
        self._target = target

    def _implementation(self) -> ModuleType:
        return import_module(f".{self._target}", __name__)

    def __getattr__(self, name: str):
        if name == "__all__":
            return [item for item in dir(self._implementation()) if not item.startswith("_")]
        if name.startswith("__"):
            raise AttributeError(name)
        return getattr(self._implementation(), name)

    def __dir__(self) -> list[str]:
        return sorted(set(super().__dir__()) | set(dir(self._implementation())))


_MODULE_ALIASES = {
    "config": "00_config",
    "dataset": "05_dataset",
    "evidence": "06_evidence",
    "graph": "07_graph",
    "model": "08_model",
    "completeness": "09_completeness",
    "clustering": "10_clustering",
    "pipeline": "11_pipeline",
    "metrics": "12_metrics",
    "cli": "13_cli",
}
for _name, _target in _MODULE_ALIASES.items():
    _alias = _ModuleAlias(f"{__name__}.{_name}", _target)
    sys.modules[_alias.__name__] = _alias
    globals()[_name] = _alias

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
