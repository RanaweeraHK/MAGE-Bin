"""Loading and validation of model-ready MAGE-Bin datasets."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import MageBinConfig
from .evidence import normalize_rows


@dataclass(slots=True)
class BinningDataset:
    """Validated features and metadata for one metagenomic cohort."""

    name: str
    directory: Path
    manifest: dict[str, Any]
    metadata: pd.DataFrame
    attributes: np.ndarray
    coverage: np.ndarray
    strain_labels: np.ndarray | None
    species_labels: np.ndarray | None
    feature_seconds: float

    @property
    def contig_ids(self) -> list[str]:
        return self.metadata["contig_id"].astype(str).tolist()

    def __len__(self) -> int:
        return len(self.metadata)


def _decode_name(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def _resolve_coverage_file(
    dataset_directory: Path,
    manifest: dict[str, Any],
    coverage_file: Path | None,
) -> Path:
    candidates = []
    if coverage_file is not None:
        candidates.append(coverage_file)
    candidates.append(dataset_directory / "coverage.csv")
    if manifest.get("coverage_csv"):
        candidates.append(Path(manifest["coverage_csv"]))
    if manifest.get("work_dir"):
        candidates.append(Path(manifest["work_dir"]) / "coverage.csv")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    attempted = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"coverage table not found; checked: {attempted}")


def load_binning_dataset(
    dataset_directory: str | Path,
    config: MageBinConfig,
    *,
    coverage_file: str | Path | None = None,
) -> BinningDataset:
    """Load standardized TNF and abundance features for one cohort.

    Ground-truth columns are retained only for optional evaluation. They are
    never included in the feature matrix or used by training and clustering.
    """

    began = time.perf_counter()
    directory = Path(dataset_directory).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"dataset directory does not exist: {directory}")

    manifest_path = directory / "config.json"
    metadata_path = directory / "contig_metadata.tsv"
    tnf_path = directory / "tnf_counts.npz"
    for required in (manifest_path, metadata_path, tnf_path):
        if not required.is_file():
            raise FileNotFoundError(f"required dataset file is missing: {required}")

    manifest = json.loads(manifest_path.read_text())
    metadata = pd.read_csv(metadata_path, sep="\t")
    required_columns = {"node_index", "contig_id", "length"}
    missing_columns = required_columns - set(metadata.columns)
    if missing_columns:
        names = ", ".join(sorted(missing_columns))
        raise ValueError(f"contig_metadata.tsv is missing columns: {names}")
    metadata = metadata.sort_values("node_index").reset_index(drop=True)
    metadata = metadata.loc[
        pd.to_numeric(metadata["length"], errors="coerce")
        >= config.min_contig_length
    ].reset_index(drop=True)
    if len(metadata) < 2:
        raise ValueError(
            "MAGE-Bin requires at least two contigs after length filtering"
        )
    if metadata["contig_id"].astype(str).duplicated().any():
        raise ValueError("contig IDs must be unique")
    contig_ids = metadata["contig_id"].astype(str).tolist()

    with np.load(tnf_path, allow_pickle=True) as saved:
        if not {"names", "counts"}.issubset(saved.files):
            raise ValueError("tnf_counts.npz must contain names and counts arrays")
        count_rows = {
            _decode_name(name): row
            for name, row in zip(saved["names"], saved["counts"], strict=True)
        }
    missing_tnf = [name for name in contig_ids if name not in count_rows]
    if missing_tnf:
        raise KeyError(f"{len(missing_tnf)} retained contigs have no TNF row")
    counts = np.stack([count_rows[name] for name in contig_ids]).astype(float)
    if counts.ndim != 2 or counts.shape[1] == 0:
        raise ValueError("TNF counts must be a non-empty two-dimensional matrix")
    counts += 0.5
    clr = np.log(counts) - np.log(counts).mean(axis=1, keepdims=True)
    clr = (clr - clr.mean(axis=0)) / (clr.std(axis=0) + 1e-5)

    resolved_coverage = _resolve_coverage_file(
        directory,
        manifest,
        Path(coverage_file).expanduser().resolve() if coverage_file else None,
    )
    coverage_frame = pd.read_csv(resolved_coverage)
    if "contig" not in coverage_frame.columns:
        raise ValueError("coverage table must contain a contig column")
    coverage_frame["contig"] = coverage_frame["contig"].astype(str)
    coverage_frame = coverage_frame.set_index("contig")
    numeric_coverage = coverage_frame.apply(pd.to_numeric, errors="coerce")
    coverage = numeric_coverage.reindex(contig_ids).fillna(0.0).to_numpy(float)
    if coverage.shape[1] == 0:
        raise ValueError("coverage table contains no sample columns")

    per_million = 1e6 / np.maximum(coverage.sum(axis=0), 1e-9)
    depth = coverage * per_million
    profile = np.log1p(depth)
    profile -= profile.mean(axis=1, keepdims=True)
    abundance = np.concatenate(
        (
            depth / np.maximum(depth.sum(axis=1, keepdims=True), 1e-9),
            profile,
            np.log1p(depth.sum(axis=1, keepdims=True)),
        ),
        axis=1,
    )
    abundance = (abundance - abundance.mean(axis=0)) / (
        abundance.std(axis=0) + 1e-5
    )
    raw_attributes = np.clip(np.concatenate((clr, abundance), axis=1), -8.0, 8.0)

    strain_labels = None
    if "genome_label" in metadata and metadata["genome_label"].notna().any():
        strain_labels = metadata["genome_label"].to_numpy()
    species_labels = None
    if strain_labels is not None and manifest.get("reference"):
        mapping_path = Path(manifest["reference"]).parent / "strain_map.tsv"
        if mapping_path.is_file():
            mapping = pd.read_csv(mapping_path, sep="\t")
            if {"strain_id", "species_id"}.issubset(mapping.columns):
                species_by_strain = dict(zip(mapping.strain_id, mapping.species_id))
                species_labels = metadata["genome_label"].map(species_by_strain).to_numpy()

    return BinningDataset(
        name=directory.name,
        directory=directory,
        manifest=manifest,
        metadata=metadata,
        attributes=normalize_rows(raw_attributes),
        coverage=normalize_rows(profile),
        strain_labels=strain_labels,
        species_labels=species_labels,
        feature_seconds=time.perf_counter() - began,
    )

