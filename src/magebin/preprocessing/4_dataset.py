"""Stage 4: assemble the model-ready dataset from user inputs."""

from __future__ import annotations

import json
from importlib import import_module
from pathlib import Path


def prepare_dataset(
    contigs: str | Path,
    coverage: str | Path,
    graph: str | Path,
    output: str | Path,
    *,
    min_contig_length: int,
    paths: str | Path | None = None,
) -> Path:
    """Generate the files consumed by the existing MAGE-Bin pipeline."""

    fasta = Path(contigs).expanduser().resolve()
    coverage_path = Path(coverage).expanduser().resolve()
    graph_path = Path(graph).expanduser().resolve()
    paths_path = Path(paths).expanduser().resolve() if paths is not None else None
    if paths_path is None:
        for candidate in (
            fasta.parent / "contigs.paths",
            graph_path.parent / "contigs.paths",
        ):
            if candidate.is_file():
                paths_path = candidate.resolve()
                break
    dataset = Path(output).expanduser().resolve() / "preprocessed"

    prepare_contigs = import_module(f"{__package__}.1_contigs").prepare_contigs
    prepare_coverage = import_module(f"{__package__}.2_coverage").prepare_coverage
    prepare_graph = import_module(f"{__package__}.3_graph").prepare_graph

    retained = prepare_contigs(fasta, dataset, min_contig_length=min_contig_length)
    prepare_coverage(coverage_path, dataset, retained)
    edge_count = prepare_graph(graph_path, dataset, retained, paths=paths_path)
    manifest = {
        "combined_fasta": str(fasta),
        "coverage_csv": str(dataset / "coverage.csv"),
        "source_contigs": str(fasta),
        "source_coverage": str(coverage_path),
        "source_graph": str(graph_path),
        "source_paths": str(paths_path) if paths_path is not None else None,
        "assembly_edges": edge_count,
    }
    (dataset / "config.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return dataset
