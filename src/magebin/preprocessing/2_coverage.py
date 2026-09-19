"""Stage 2: validate and normalize a user supplied coverage table."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def prepare_coverage(source: Path, dataset: Path, retained: list[str]) -> None:
    """Accept tab or comma separated coverage with contig and sample columns."""

    if not source.is_file():
        raise FileNotFoundError(f"coverage table not found: {source}")
    frame = pd.read_csv(source, sep=None, engine="python", dtype=str)
    if "contig" not in frame.columns:
        raise ValueError("coverage table needs a 'contig' column")
    samples = [column for column in frame.columns if column != "contig"]
    if not samples:
        raise ValueError("coverage table needs at least one sample column")
    if frame.columns.duplicated().any():
        raise ValueError("coverage table has duplicate column names")
    if frame["contig"].isna().any() or frame["contig"].duplicated().any():
        raise ValueError("coverage table has missing or duplicate contig identifiers")
    missing = set(retained) - set(frame["contig"])
    if missing:
        example = ", ".join(sorted(missing)[:3])
        raise ValueError(f"coverage table is missing {len(missing)} contigs: {example}")
    for sample in samples:
        try:
            frame[sample] = pd.to_numeric(frame[sample], errors="raise")
        except (ValueError, TypeError) as error:
            raise ValueError(f"coverage sample {sample!r} must be numeric") from error
    values = frame[samples].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("coverage values must be finite and nonnegative")
    frame.to_csv(dataset / "coverage.csv", index=False)
