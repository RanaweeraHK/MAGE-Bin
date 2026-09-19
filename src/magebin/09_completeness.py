"""Reference-free protection of apparently complete viral contigs."""

from __future__ import annotations

import gzip
import time
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd

from .config import MageBinConfig
from .dataset import BinningDataset

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""

    return sequence.translate(_COMPLEMENT)[::-1]


def find_terminal_repeats(
    sequence: str,
    *,
    seed_length: int,
) -> tuple[int, int]:
    """Return longest exact direct and inverted terminal-repeat lengths."""

    sequence = sequence.upper()
    size = len(sequence)
    if size < 4 * seed_length:
        return 0, 0
    direct = 0
    position = sequence.rfind(sequence[:seed_length])
    if position > 0:
        span = size - position
        if span <= size // 2 and sequence[position:] == sequence[:span]:
            direct = span

    inverted = 0
    if sequence.startswith(reverse_complement(sequence[-seed_length:])):
        inverted = seed_length
        while inverted < size // 2:
            prefix_base = sequence[inverted]
            suffix_base = sequence[size - inverted - 1]
            if prefix_base != reverse_complement(suffix_base):
                break
            inverted += 1
    return direct, inverted


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield identifier and sequence records from plain or gzipped FASTA."""

    fasta_path = Path(path)
    opener = gzip.open if fasta_path.suffix == ".gz" else open
    name: str | None = None
    parts: list[str] = []
    with opener(fasta_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def scan_terminal_repeats(
    dataset: BinningDataset,
    config: MageBinConfig,
    *,
    cache_file: str | Path | None = None,
) -> pd.DataFrame | None:
    """Scan cohort FASTA records, optionally caching the repeat table."""

    resolved_cache = Path(cache_file) if cache_file is not None else None
    if resolved_cache is not None and resolved_cache.is_file():
        return pd.read_csv(resolved_cache, sep="\t", dtype={"contig_id": str})
    fasta = dataset.manifest.get("combined_fasta")
    if not fasta or not Path(fasta).is_file():
        return None
    rows = []
    for name, sequence in read_fasta(fasta):
        direct, inverted = find_terminal_repeats(
            sequence, seed_length=config.terminal_repeat_seed
        )
        rows.append((name, len(sequence), direct, inverted))
    frame = pd.DataFrame(
        rows,
        columns=["contig_id", "length", "direct_repeat", "inverted_repeat"],
    )
    if resolved_cache is not None:
        resolved_cache.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(resolved_cache, sep="\t", index=False)
    return frame


def detect_complete_contigs(
    dataset: BinningDataset,
    config: MageBinConfig,
    *,
    cache_file: str | Path | None = None,
) -> dict[str, object] | None:
    """Identify terminal-repeat contigs that must remain singleton bins."""

    began = time.perf_counter()
    repeats = scan_terminal_repeats(dataset, config, cache_file=cache_file)
    if repeats is None:
        return None
    names = np.asarray(dataset.contig_ids)
    lengths = dataset.metadata["length"].to_numpy(float)
    aligned = (
        repeats.drop_duplicates("contig_id")
        .set_index("contig_id")
        .reindex(names)
        .fillna(0.0)
    )
    repeat_span = np.maximum(
        aligned["direct_repeat"].to_numpy(float),
        aligned["inverted_repeat"].to_numpy(float),
    )
    locked = (repeat_span >= config.terminal_repeat_minimum) & (
        repeat_span <= config.terminal_repeat_max_fraction * lengths
    )
    return {
        "locked": locked,
        "locked_contigs": int(locked.sum()),
        "locked_fraction": float(locked.mean()),
        "scan_seconds": time.perf_counter() - began,
    }


def isolate_complete_contigs(
    labels: np.ndarray,
    completeness_gate: dict[str, object],
) -> np.ndarray:
    """Move every detected complete contig into its own singleton bin."""

    result = np.asarray(labels, dtype=np.int64).copy()
    next_label = int(result.max()) + 1 if len(result) else 0
    for node in np.flatnonzero(np.asarray(completeness_gate["locked"], dtype=bool)):
        result[node] = next_label
        next_label += 1
    return np.unique(result, return_inverse=True)[1]

