"""Stage 3: project GFA links onto FASTA contig identifiers."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import torch


def _segment(token: str) -> str:
    return token.strip().rstrip("+-")


def _spades_paths(path: Path, names: set[str]) -> dict[str, set[str]]:
    """Read SPAdes contigs.paths when GFA segment IDs differ from contig IDs."""

    mapped: dict[str, set[str]] = defaultdict(set)
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line in names:
                current = line
                continue
            if line.endswith("'") and line[:-1] in names:
                current = None
                continue
            if current is not None:
                mapped[current].update(
                    _segment(token)
                    for token in line.replace(";", ",").split(",")
                    if token.strip()
                )
    return mapped


def prepare_graph(
    source: Path,
    dataset: Path,
    names: list[str],
    *,
    paths: Path | None = None,
) -> int:
    """Write a contig-level graph compatible with load_assembly_edges."""

    if not source.is_file():
        raise FileNotFoundError(f"GFA graph not found: {source}")
    if paths is not None and not paths.is_file():
        raise FileNotFoundError(f"contig paths file not found: {paths}")
    index = {name: number for number, name in enumerate(names)}
    segment_to_contigs: dict[str, set[int]] = defaultdict(set)
    for name, number in index.items():
        segment_to_contigs[name].add(number)
    links: list[tuple[str, str]] = []
    with source.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if fields[0] == "P" and len(fields) >= 3 and fields[1] in index:
                for token in fields[2].split(","):
                    segment_to_contigs[_segment(token)].add(index[fields[1]])
            elif fields[0] == "L":
                if len(fields) < 6:
                    raise ValueError(f"invalid GFA L record on line {number}")
                links.append((fields[1], fields[3]))
    if paths is not None:
        for name, segments in _spades_paths(paths, set(index)).items():
            for segment in segments:
                segment_to_contigs[segment].add(index[name])

    edges: set[tuple[int, int]] = set()
    for left_segment, right_segment in links:
        left_contigs = segment_to_contigs.get(left_segment, set())
        right_contigs = segment_to_contigs.get(right_segment, set())
        for left in left_contigs:
            for right in right_contigs:
                if left != right:
                    edges.add((min(left, right), max(left, right)))
    if not edges:
        raise ValueError(
            "GFA links do not map to pairs of retained FASTA contigs. Supply a GFA with matching "
            "segment or P path names, or pass --paths for SPAdes contigs.paths."
        )
    ordered = sorted(edges)
    edge_index = torch.tensor(ordered, dtype=torch.long).T.contiguous()
    graph = SimpleNamespace(
        contig_names=names,
        edge_index=edge_index,
        edge_attr=torch.ones((len(ordered), 1), dtype=torch.float32),
    )
    torch.save(graph, dataset / "viral_graph.pt")
    return len(ordered)
