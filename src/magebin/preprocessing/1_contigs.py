"""Stage 1: read contigs and calculate canonical tetranucleotide counts."""

from __future__ import annotations

import gzip
import itertools
from pathlib import Path

import numpy as np
import pandas as pd


def _kmer_columns() -> tuple[dict[str, int], int]:
    complement = str.maketrans("ACGT", "TGCA")
    columns: dict[str, int] = {}
    order: dict[str, int] = {}
    for bases in itertools.product("ACGT", repeat=4):
        kmer = "".join(bases)
        canonical = min(kmer, kmer.translate(complement)[::-1])
        if canonical not in order:
            order[canonical] = len(order)
        columns[kmer] = order[canonical]
    return columns, len(order)


def prepare_contigs(fasta: Path, dataset: Path, *, min_contig_length: int) -> list[str]:
    """Write contig metadata and TNF counts; return retained IDs."""

    if not fasta.is_file():
        raise FileNotFoundError(f"contig FASTA not found: {fasta}")
    columns, width = _kmer_columns()
    names: list[str] = []
    lengths: list[int] = []
    vectors: list[np.ndarray] = []
    seen: set[str] = set()
    opener = gzip.open if fasta.suffix == ".gz" else open
    name: str | None = None
    parts: list[str] = []

    def add_record() -> None:
        if name is None:
            return
        sequence = "".join(parts).upper()
        if not sequence:
            raise ValueError(f"contig {name!r} has no sequence")
        vector = np.zeros(width, dtype=np.int64)
        for index in range(len(sequence) - 3):
            column = columns.get(sequence[index : index + 4])
            if column is not None:
                vector[column] += 1
        names.append(name)
        lengths.append(len(sequence))
        vectors.append(vector)

    with opener(fasta, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                add_record()
                name = line[1:].split()[0] if line[1:].strip() else ""
                if not name:
                    raise ValueError("FASTA contains an empty contig identifier")
                if name in seen:
                    raise ValueError(f"duplicate contig identifier in FASTA: {name}")
                seen.add(name)
                parts = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence appears before its header")
                parts.append(line)
    add_record()
    retained = [
        contig
        for contig, length in zip(names, lengths, strict=True)
        if length >= min_contig_length
    ]
    if len(retained) < 2:
        raise ValueError(
            "MAGE-Bin needs at least two contigs at or above "
            f"{min_contig_length} bases; found {len(retained)}"
        )
    dataset.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {"node_index": range(len(names)), "contig_id": names, "length": lengths}
    ).to_csv(dataset / "contig_metadata.tsv", sep="\t", index=False)
    np.savez_compressed(
        dataset / "tnf_counts.npz",
        names=np.asarray(names, dtype=str),
        counts=np.stack(vectors),
    )
    return retained
