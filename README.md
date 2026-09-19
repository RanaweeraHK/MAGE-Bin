# MAGE-Bin

MAGE-Bin is a label-free method for viral metagenomic binning. It combines
tetranucleotide composition, multi-sample abundance, statistical graph evidence,
and biologically supported assembly-graph links to group viral contigs into bins.

MAGE-Bin uses self-supervised masked reconstruction to learn contig
representations. A learnable fusion gate controls how much graph information
contributes to the final representation. Ground-truth genome labels and CheckV
scores are not used during training, model selection, or clustering.

## Installation

MAGE-Bin requires Python 3.10 or later.

Install the current release from PyPI:

```bash
pip install magebin
```

Check the installation:

```bash
magebin --version
magebin doctor
```

### PyTorch

MAGE-Bin uses PyTorch and PyTorch Geometric. The standard installation above
allows `pip` to resolve these dependencies automatically.

If you require a particular CPU or CUDA build of PyTorch, install the appropriate
PyTorch build for your system before installing MAGE-Bin.

## Quick start

MAGE-Bin operates on a model-ready dataset directory.

A basic run is:

```bash
magebin bin /path/to/dataset \
    --output /path/to/results
```

For example, to force CPU execution:

```bash
magebin bin /path/to/dataset \
    --output /path/to/results \
    --device cpu \
    --seed 0
```

Run the following command to see all available options:

```bash
magebin bin --help
```

## Input

MAGE-Bin currently expects a preprocessed, model-ready dataset directory.

The directory must contain:

- `config.json`
- `contig_metadata.tsv`
- `tnf_counts.npz`
- `viral_graph.pt` when assembly-graph information is available

`contig_metadata.tsv` must include:

- `node_index`
- `contig_id`
- `length`

`tnf_counts.npz` must contain the `names` and `counts` arrays used for
tetranucleotide composition.

### Coverage

Coverage can be supplied explicitly:

```bash
magebin bin /path/to/dataset \
    --coverage /path/to/coverage.csv \
    --output /path/to/results
```

When `--coverage` is not provided, MAGE-Bin searches for coverage information
using the dataset configuration.

The preprocessing workflow used to construct model-ready datasets is available
in the `Dataset_Processing/` directory of the GitHub repository.

> **Current limitation:** this version does not yet take raw FASTA and read
> files directly through the `magebin bin` command. Input must first be
> converted to the model-ready dataset format.

## Output

MAGE-Bin writes its results to the directory specified with `--output`.

The main outputs are:

### `assignments.tsv`

Contains the final contig-to-bin assignments, including:

- `contig_id`
- `bin_id`
- contig `length`

### `run.json`

Records information about the run, including model parameters, learned fusion
weights, graph statistics, timing information, and output summary.

### `terminal_repeats.tsv`

Stores cached terminal-repeat evidence when the required sequence information is
available.

## Method overview

MAGE-Bin combines biological features and graph information in a self-supervised
binning framework.

The main steps are:

1. **Biological feature construction**  
   Tetranucleotide composition and multi-sample abundance are used to represent
   each viral contig.

2. **Statistical graph construction**  
   Candidate relationships between contigs are evaluated using biological
   similarity and empirical evidence calibration.

3. **Assembly-graph support**  
   Assembly relationships are treated as candidate links. An assembly link is
   retained only when it has positive independent biological evidence from the
   available composition and coverage information.

4. **Self-supervised representation learning**  
   A graph encoder is trained using masked feature reconstruction. No
   ground-truth genome labels are required.

5. **Learnable graph fusion**  
   The encoder learns a convex gate between the identity representation and
   graph-derived representation.

6. **Unknown-K clustering**  
   The learned representations and biological evidence are used to construct a
   sparse clustering graph. MAGE-Bin determines the bins without requiring the
   number of viral genomes to be specified beforehand.

7. **Completeness protection**  
   Terminal-repeat evidence can protect likely complete viral contigs from
   inappropriate merging.

Assembly edges are treated as biological evidence rather than ground truth.
Supported assembly edges are added to the statistical graph without deleting or
reweighting the existing statistical edges.

CheckV is used only for downstream biological evaluation. CheckV scores are not
used for training, model selection, graph construction, or clustering.

## Reproducibility

Important training parameters can be controlled from the command line.

For example:

```bash
magebin bin /path/to/dataset \
    --output /path/to/results \
    --device cpu \
    --seed 0 \
    --minimum-epochs 40 \
    --maximum-epochs 100 \
    --graph-update-interval 10
```

Using a fixed seed and recording `run.json` helps reproduce individual runs.

## Checking the installation

MAGE-Bin includes a diagnostic command:

```bash
magebin doctor
```

This reports the installed Python dependencies and the availability of optional
external tools.

## Development

Clone the repository:

```bash
git clone https://github.com/RanaweeraHK/MAGE-Bin.git
cd MAGE-Bin
```

Install MAGE-Bin with the development and testing dependencies:

```bash
python -m pip install -e '.[test]'
```

Run the automated checks:

```bash
python -m ruff check src tests
python -m pytest
```

Build and validate the distributions:

```bash
python -m build
python -m twine check --strict dist/*
```

The automated tests cover numerical utilities, graph construction,
assembly-edge integration, learnable fusion behavior, completeness gating,
metrics, CLI behavior, and a small end-to-end model run.

Large biological benchmarks and CheckV evaluations are kept separate from the
unit and integration tests because they require larger datasets and external
databases.

## Research code and notebooks

The installable MAGE-Bin implementation is maintained under `src/magebin/`.

The `Notebooks/` directory contains research experiments and development
notebooks. These notebooks are not required when MAGE-Bin is installed from
PyPI and are not imported by the runtime package.

Dataset preparation and benchmark-related scripts are maintained separately
from the installable model.

## Citation

If you use MAGE-Bin in research, please cite the MAGE-Bin paper.

The full citation and DOI will be added here when the paper becomes available.

## License

MAGE-Bin is distributed under the MIT License. See `LICENSE` for details.