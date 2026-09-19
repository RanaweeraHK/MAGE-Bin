#!/usr/bin/env python
"""Phase 2 - co-assembly, coverage and graph: runs ../Real_Data's phase 2.

    VB_DATASET=phage_mock_illumina python 2_assemble_and_build_graph.py
    VB_DATASET=dsmz_dsrna python 2_assemble_and_build_graph.py --dry-run

A mock is real sequencing, so every step from reads to a model-ready dataset is
the work the real-dataset pipeline already does, and the artifacts must be
identical for a model to read either unchanged. This runs that script with
three environment variables pointed at this pipeline's registry and
directories, rather than keeping a second copy of it that would drift.

Ground truth is added afterwards by phase 3, which fills the `genome_label`
column this leaves empty.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DELEGATE = HERE.parent / "Real_Data" / "2_assemble_and_build_graph.py"


def main() -> int:
    if not DELEGATE.is_file():
        raise SystemExit(f"error: {DELEGATE} not found")
    env = dict(os.environ)
    env["VB_REGISTRY"] = str(HERE / "datasets.tsv")
    env["VB_WORK_PREFIX"] = "mock"           # Data/work/mock_<dataset>/
    env["VB_DATASET_ROOT"] = "Mock_dataset"  # Data/Mock_dataset/<study>/
    env["VB_DATASET_TYPE"] = "mock_community"
    # has_ground_truth stays False in the manifest until phase 3 writes labels.
    os.execve(sys.executable, [sys.executable, str(DELEGATE), *sys.argv[1:]], env)


if __name__ == "__main__":
    sys.exit(main())
