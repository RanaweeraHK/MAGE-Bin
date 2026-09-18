# MAGE-Bin

MAGE-Bin is a label-free viral metagenomic binner that combines standardized
tetranucleotide composition, multi-sample abundance, a masked graph
autoencoder, and biologically filtered assembly-graph links. Its encoder learns
a convex gate between the identity representation and graph messages. Neither
ground-truth labels nor CheckV scores are used for training, model selection, or
clustering.

The installable implementation in `src/magebin/` is the source of truth. The
notebooks under `Notebooks/` record the research process and previous
experiments; they are not imported at runtime.

## Installation

For development:

```bash
python -m pip install -e '.[test]'
magebin --version
magebin doctor
```

MAGE-Bin is distributed under the MIT License; see `LICENSE`.

## Input contract

The current CLI consumes the model-ready dataset contract produced by
`Dataset_Processing/`. A dataset directory must contain:

- `config.json`
- `contig_metadata.tsv`, including `node_index`, `contig_id`, and `length`
- `tnf_counts.npz`, including `names` and `counts` arrays
- `viral_graph.pt` when assembly links are available

Coverage is resolved from `--coverage`, `<dataset>/coverage.csv`, the
`coverage_csv` manifest field, or `<work_dir>/coverage.csv`, in that order.

## CLI

```bash
magebin bin Data/Processed_data/<dataset> --output Outputs/magebin/<dataset>
```

Useful reproducibility controls:

```bash
magebin bin Data/Processed_data/<dataset> \
  --output Outputs/magebin/<dataset> \
  --device cpu \
  --seed 0 \
  --minimum-epochs 40 \
  --maximum-epochs 100 \
  --graph-update-interval 10
```

Outputs are:

- `assignments.tsv`: `contig_id`, stable `bin_id`, and contig `length`
- `run.json`: parameters, learned fusion weights, graph statistics, timing,
  and output summary
- `terminal_repeats.tsv`: cached reference-free completeness evidence when a
  FASTA is available

## Testing

```bash
python -m ruff check src tests
python -m pytest
python -m build
python -m twine check --strict dist/*
```

Tests cover numerical utilities, graph invariants, assembly-edge fusion,
learnable-gate behavior, completeness gating, metrics, CLI behavior, and a tiny
end-to-end model run. Full biological benchmarks and CheckV evaluation remain
separate from pull-request tests because they require large datasets and
external databases.

## Design boundaries

- Assembly edges are candidates, not truth. They are retained only when their
  combined composition/coverage evidence is positive.
- Supported assembly edges are unioned with statistical edges without deleting
  or reweighting the statistical graph.
- The fusion gate is learned during masked reconstruction and starts at the
  notebook's original 0.10 graph-to-identity ratio.
- Complete-looking terminal-repeat contigs can be isolated from merging.
- CheckV is a downstream evaluation tool and never participates in fitting.

## Publishing

The GitHub workflows follow the pattern used by maintained metagenomic tools:
pull requests run linting, tests, CLI smoke checks, and package builds; a
semantic version tag such as `v0.1.0` builds the distributions, publishes them
to PyPI through Trusted Publishing, and creates a GitHub release. Configure a
protected `pypi` GitHub environment and the matching PyPI trusted publisher
before creating a release tag.

Bioconda submission comes after the first immutable PyPI release. Copy the
recipe template under `packaging/bioconda/` into a fork of
`bioconda-recipes`, replace its source checksum and maintainer placeholder, run
`bioconda-utils lint`, and open a pull request there. The template deliberately
cannot be submitted until this project has a released source archive.
