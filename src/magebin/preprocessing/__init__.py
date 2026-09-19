"""Prepare standard user inputs for the MAGE-Bin dataset loader.

The numbered modules correspond to the order of the preprocessing stages.
Python cannot use their names in regular import statements, so load them here.
"""

from __future__ import annotations

from importlib import import_module

prepare_dataset = import_module(f"{__name__}.4_dataset").prepare_dataset

__all__ = ["prepare_dataset"]
