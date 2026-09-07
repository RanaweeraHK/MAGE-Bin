# Multi-biome viral-binning benchmark — SIMULATED preprocessing pipeline

Turns the raw genome sources in `Data/Raw_dataset/multibiome_source_data/`
into three model-ready datasets under `Data/Processed_data/`, one per
complexity tier, each holding the same artifact set (`viral_graph.pt`,
`assembly_heterograph.pt`, `contig_metadata.tsv`, `protein_clusters.tsv`,
`viral_features.tsv`, `config.json`, `manifest.json`; plus
`assembly_graph.gfa`, `contigs.paths`, `contigs_filt.fasta`, `coverage.csv`,
`labels.csv` and `strain_map.tsv` left in that tier's work directory).


## Complexity tiers

The benchmark is built at three complexity levels. Complexity is two knobs:
how many genomes are in the community, and what share of them carry multiple
near-identical strains (the "strain complexity" that drives assembly
fragmentation and is the hard part of binning).

| tier | genomes | genomes/biome | strain complexity | strains | reference | reads/sample | build time | dataset directory |
|---|---|---|---|---|---|---|---|---|
| `low` | 50 | 10 | 5% | 54 | 3.18 Mb | 40,000 | ~6 min | `multibiome_low__coasm__metaspades__noviral/` |
| `medium` | 200 | 40 | 20% | 296 | 20.02 Mb | 250,000 | ~36 min | `multibiome_medium__coasm__metaspades__noviral/` |
| `high` | 998 | 200 | 40% | 1,995 | 128.66 Mb | 1,500,000 | ~4 h | `multibiome_high__coasm__metaspades__noviral/` |

What each tier came out as — re-measured against the datasets currently in
`Data/Processed_data/` by running `6_acceptance_check.py` at each tier, which
is how to reproduce this table after any rebuild:

| tier | contigs | labelled | species recovered | multi-contig species | GFA merge ceiling | verdict |
|---|---|---|---|---|---|---|
| `low` | 283 | 100% | 49 | 14 (28% of pool) | 1% | STRONG |
| `medium` | 814 | 100% | 193 | 57 (28% of pool) | 3% | STRONG |
| `high` | 4,046 | 100% | 950 | 347 (35% of pool) | 3% | STRONG |


`high` genomes total 998, not 1,000: two of the 200 freshwater draws are
exact duplicates of genomes already taken from another source and are
dropped by the MD5 guard. The low and medium draws hit no such collisions.

## Before you run anything

Setup is shared with the sibling real-data pipeline and lives one level up —
see [`../README.md`](../README.md) ("Setup — do this once"). In short:

```bash
cd /path/to/Viralbinning/Dataset_Processing
bash setup_environments.sh --only viralbin asm_env   # the two envs this pipeline needs
bash Simulated_Data/0_check_environment.sh           # preflight: seconds, not hours
```

`setup_environments.sh` also writes `Dataset_Processing/env.sh`, the one
gitignored file holding the machine-specific values every script here reads:
`VB_CONDA_BIN` and `VB_ASM_CONDA_BIN` (the two conda `bin/` directories),
`VB_THREADS` and `VB_MEMORY_GB`. Nothing in this directory hardcodes a path to
your machine, and every tool lookup falls back to `PATH`, so an activated
conda environment works even with no `env.sh` at all. Real environment
variables win over the file:

```bash
VB_THREADS=32 bash run_all.sh --tier low
```

The tier/simulator/assembler (`VB_TIER`, `VB_SIM`, `VB_ASM`) are *not* in
`env.sh` — they are per-run choices, passed on the command line as below.

## One command

Which datasets get built is a list at the top of `run_all.sh`. Edit it, then
run one command:

```bash
DATASETS=(
    "low:iss:metaspades"
    "medium:iss:metaspades"
    "high:iss:metaspades"
    "high:badread:metaflye"
)
```

```bash
cd /path/to/Viralbinning/Dataset_Processing/Simulated_Data
bash run_all.sh
```

Each entry is `tier:simulator:assembler`. Phases run in order, skipping any
phase whose output already exists, so it is safe to re-run after an
interruption. Useful flags:

```bash
bash run_all.sh --list                          # show what would be built, exit
bash run_all.sh --dataset high:badread:metaflye # just this one
bash run_all.sh --tier medium                   # filter the list by tier
bash run_all.sh --force                         # rebuild from scratch
bash run_all.sh --from 4                        # resume starting at phase 4
```

`--force` and `--from` apply to *every* dataset being built, so pair them
with `--dataset` when you only mean to redo one.

### Simulators and assemblers

| simulator | reads | assembler | minimap2 preset |
|---|---|---|---|
| `iss` (InSilicoSeq, hiseq) | Illumina paired-end 2x126 bp | `metaspades --only-assembler` | `-x sr` |
| `badread` (nanopore2023) | ONT single-end, 15 kb mean | `metaflye --meta --nano-hq` | `-x map-ont` |

The read type ties the two together, so only `iss+metaspades` and
`badread+metaflye` are valid. Any other pairing is rejected by `run_all.sh`
*before* any work starts rather than failing part-way through an assembly.

A long-read dataset is given the **same per-sample budget in bases** as the
short-read dataset at its tier (`reads/sample` x 126 bp), so per-genome depth
is comparable across simulators and only read length and error profile
differ. Phase 3's abundance tables are a share of each sample's *reads*;
Badread instead wants a per-sequence `depth=` and draws reads proportional to
depth x length, so phase 4 converts with `depth = abundance / length`, which
reproduces the iss read allocation exactly. Strains absent from a sample get
`depth=0` and yield no reads.

Every numbered script also takes its selection from the `VB_TIER`, `VB_SIM`
and `VB_ASM` environment variables (defaults `high`, `iss`, `metaspades`)
when run standalone:

```bash
VB_TIER=low python 1_build_genome_pool.py
VB_SIM=badread VB_ASM=metaflye python 5_assemble_and_build_graph.py
```

### Dataset directory naming

    multibiome<tier><sim>__coasm__<assembler>__noviral

where `<tier>` is `_low`/`_medium`/`_high` and `<sim>` is empty for `iss` and
`__<name>` otherwise — the `iss` segment being empty is what keeps the
short-read directories the plainer of the two names. The `DATASETS` list above
therefore produces:

    multibiome_low__coasm__metaspades__noviral
    multibiome_medium__coasm__metaspades__noviral
    multibiome_high__coasm__metaspades__noviral
    multibiome_high__badread__coasm__metaflye__noviral

**Phases 1-3 depend only on the tier**, so several simulators/assemblers at
one tier share one genome pool and one set of abundance tables and are built
once. Phase 4 writes per-simulator read directories (`reads/` for `iss`,
`reads_<sim>/` otherwise) and phases 5-7 write per-dataset output
directories, so datasets never overwrite each other.

The high tier's phase 4 (read simulation, ~2h) and phase 5 (assembly, ~1.5h)
are long — run `run_all.sh` under `nohup`/`tmux`/`screen` rather than in a
session that might get closed. Low (~6 min) and medium (~36 min) finish
first, so a full three-tier run has two datasets on disk within the hour.


## Scripts, in order

| Script | What it does | Needs |
|---|---|---|
| `0_check_environment.sh` | Preflight: verifies raw sources, tools (via `../env.sh`, then `PATH`), and the downstream pipeline scripts are all resolvable before burning hours in phases 1–5. Tier-independent, runs once. | — |
| `1_build_genome_pool.py` | Filters each of the 5 sources by its own quality flag, keeps genomes ≥10 kb, samples 10/40/200 per biome by tier (seeded), MD5 dedupes across sources. | python stdlib only |
| `2_expand_strains.py` | 5%/20%/40% of the pool genomes (by tier) become 2–5 near-identical strains at 98.0–99.5% ANI (real cause of assembly fragmentation). | python stdlib only |
| `3_generate_abundance.py` | Per-sample presence (p=0.70) + log-normal abundance for 15 samples → InSilicoSeq abundance files. | python stdlib only |
| `4_simulate_reads.sh` | Per sample, `iss generate` → paired-end Illumina reads, or `badread simulate` → ONT long reads, at the tier's depth-matched budget. Selected with `VB_SIM`. | `iss` or `badread` |
| `5_assemble_and_build_graph.py` | Co-assembly (metaSPAdes or metaFlye, `VB_ASM`) → length filter → coverage (minimap2) → labels (align to strains) → PyG graph + viral features (mmseqs). Self-contained (see below). | `metaspades.py` or `flye`, `minimap2`, `samtools`, `seqkit`, `mmseqs`, python w/ torch+PyG+pysam+pyrodigal-gv+sklearn |
| `6_acceptance_check.py` | Sanity check: is this actually a hard binning problem (multi-contig genomes, GFA reachability) at strain/species tiers, per biome. | python w/ torch+pandas+networkx |

`run_all.sh` is the only non-numbered file — it's the "one command" chain
runner, not a helper module each script depends on. Every numbered script
still runs standalone (e.g. `python 1_build_genome_pool.py`) with nothing
else to look up.


## Raw data provenance

| Biome / stratum | Source | Raw file(s) under `Data/Raw_dataset/multibiome_source_data/` | Genomes kept | Link |
|---|---|---|---|---|
| Reference (cultured) | NCBI RefSeq Viral | `reference/refseq_viral_raw.fasta` | 200 | https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ |
| Human gut | MGV (Metagenomic Gut Virus catalog) | `human_gut_mgv/*.fna` (280 files) | 200 | https://portal.nersc.gov/MGV/ |
| Freshwater | Freshwater viral genomes recovered by COBRA | `freshwater/freshwater_self_circular.fasta` | 198 | https://figshare.com/articles/dataset/viral_genomes_fasta/23282789 |
| Marine | IMG/VR v1 (2018-07-01_4) | `marine_imgvr/marine_imgvr.fasta` (already Marine + High-quality/Finished + ≥10kb extracted) | 200 | https://img.jgi.doe.gov/vr/ |
| Soil | Global Soil Virome (GSV) | `soil/soil_GSV.fasta` + `soil/soil_GSV_quality.csv` | 200 | https://zenodo.org/records/10463783 |
| **Total** | | | **998** | |

The "genomes kept" column is the `high` tier; `low` and `medium` draw 10 and
40 per biome from exactly the same candidate sets (50 and 200 total, no MD5
collisions at those sizes).

Verified by running `1_build_genome_pool.py` against this data: candidate
counts per biome (7,658 / 280 / 4,792 / 3,517 / 1,365) and high-tier kept
counts (200 / 200 / 198 / 200 / 200 → 998 total) match the original report
exactly.

## Where things land

`<S>` below is the tier suffix: `_low`, `_medium` or `_high`. Each
tier gets its own work and output directory and never touches another's.

```
Data/Raw_dataset/multibiome_source_data/   <- raw sources (input, untouched, shared by all tiers)
Data/work/multibiome<S>/                   <- scratch (phases 1-4)
    genome_pool/community_pool.fasta       <-   50 / 200 / 998 genomes
    genome_pool/genome_metadata.tsv        <-   genome_id, biome, source_id, length, quality
    genome_pool/expanded_genomes.fasta     <-   54 / 296 / 1,995 strain sequences
    genome_pool/strain_map.tsv             <-   strain_id -> species_id, biome, ani_to_base, ...
    abundance/sample0..14_abundance.txt    <-   per-sample strain proportions
    reads/sample0..14_R{1,2}.fastq.gz      <-   simulated Illumina reads
    standardize_work/<run>/                <-   phase 5's scratch; cleaned down after phase 5 to
                                                contigs_filt.fasta, coverage.csv, labels.csv,
                                                assembly_graph.gfa, contigs.paths (all still read
                                                by the model / analysis notebooks)
Data/Processed_data/
    multibiome<S>__coasm__metaspades__noviral/   <- FINAL dataset (kept)
        viral_graph.pt, assembly_heterograph.pt,
        contig_metadata.tsv, protein_clusters.tsv, viral_features.tsv,
        config.json, manifest.json
```

`Data/work/` and `Data/Processed_data/` are the only folders this pipeline
creates — nothing lands inside `Dataset_Processing/` itself.

Ground truth (`strain_map.tsv`) stays in the tier's `genome_pool/`, which is
where the notebooks read it from; it is not copied into the dataset
directory.

## Tool dependencies

`0_check_environment.sh` checks all of these and prints an install command
for anything missing. `iss` (system), `mmseqs` + `python` w/
torch/PyG/pysam/pyrodigal-gv/sklearn (`viralbin` conda env), and
`metaspades.py`/`minimap2`/`samtools`/`seqkit` (`asm_env` conda env, see
below) are expected.

`metaspades.py`/`minimap2`/`samtools`/`seqkit` live in their **own** env
(`asm_env`), not inside `viralbin`: bioconda's current `samtools`/`htslib`
builds need a newer `openssl`/`libdeflate` than what's already pinned inside
`viralbin` (by its `torch`/PyG stack), so `conda install -n viralbin ...`
fails to solve. A clean env has no such pins and resolves in seconds:

```bash
mamba create -n asm_env -c bioconda -c conda-forge \
  "spades>=3.15" "minimap2>=2.24" "samtools>=1.17" "seqkit>=2.5"
```

`setup_environments.sh` builds it, and `VB_ASM_CONDA_BIN` in `../env.sh`
points at it.
`5_assemble_and_build_graph.py` prepends it to `PATH` at import time, and its
`Toolbox.resolve()` falls back from `--conda-bin` to `PATH` (trying both
`name` and `name.py` either way), so the second env is picked up
automatically.

Phases 0–3 need nothing beyond the Python standard library and already run
against the real raw data (verified — see provenance table above).

