#!/usr/bin/env python
"""Phase 4 - is this a usable binning problem? Runs ../Real_Data's phase 3.

    VB_DATASET=phage_mock_illumina python 4_acceptance_check.py

The question - is there differential coverage, an assembly graph and shared
protein content - is the same one a real dataset has to answer, so this runs
that check against this registry instead of copying it.

Read it together with phase 3's truth_summary.json. A mock can pass this and
still be a weak benchmark: phage_mock_illumina's three runs are replicates of
one community, so its coverage columns barely differ, and a reference set that
explains few contigs shows up in truth_summary.json rather than here.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DELEGATE = HERE.parent / "Real_Data" / "3_acceptance_check.py"


def main() -> int:
    if not DELEGATE.is_file():
        raise SystemExit(f"error: {DELEGATE} not found")
    env = dict(os.environ)
    env["VB_REGISTRY"] = str(HERE / "datasets.tsv")
    os.execve(sys.executable, [sys.executable, str(DELEGATE), *sys.argv[1:]], env)


if __name__ == "__main__":
    sys.exit(main())
