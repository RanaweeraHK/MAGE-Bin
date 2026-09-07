# Multi-biome viral-binning benchmark — SIMULATED preprocessing pipeline

*Real sequencing datasets are built by the sibling pipeline in
[`../Real_Data/`](../Real_Data/README.md), which has no ground truth and
is evaluated with CheckV instead. See [`../README.md`](../README.md) for
how the two relate.*

Turns the raw genome sources in `Data/Raw_dataset/multibiome_source_data/`
into three model-ready datasets under `Data/Processed_data/`, one per
complexity tier, each holding the same artifact set (`viral_graph.pt`,
`assembly_heterograph.pt`, `contig_metadata.tsv`, `protein_clusters.tsv`,
`viral_features.tsv`, `config.json`, `manifest.json`; plus
`assembly_graph.gfa`, `contigs.paths`, `contigs_filt.fasta`, `coverage.csv`,
`labels.csv` and `strain_map.tsv` left in that tier's work directory).

The recipe follows the original build's
`SImulated dataset/MULTIBIOME_BENCHMARK_REPORT.md`: complete viral genomes
across 5 biomes, expanded into strains, given differential abundance across
15 samples, sequenced in silico, and assembled — so that ground truth
(strain / species / biome) is known for every contig.

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

The multi-contig share lands at 28–35% across a 20× range in community size,
which is the depth scaling working: all three fragment the same way, so a
result that differs between tiers is a complexity effect, not a coverage one.
The `high` tier fragments somewhat more than the two smaller ones, which is
the strain fraction (40% vs 20% and 5%) showing through on top of size.

Everything except those knobs is held fixed across tiers: the same 15
samples, the same presence/log-normal abundance model, the same read
simulator and assembler, the same coverage/label/graph recipe, and the same
`MIN_CONTIG_LEN`. **`reads/sample` is scaled with the tier's total reference
size** so per-genome sequencing depth is roughly constant everywhere —
difficulty then varies by community complexity alone, not by how deeply each
genome happened to be sequenced. The high tier sits at 11,658 reads/Mb
(1.5M over 128.66 Mb); low and medium are rounded up to 40,000 and 250,000,
which is 12,579 and 12,488 reads/Mb — consistent with each other and ~7%
deeper than high, rather than matched to it exactly.
Keeping `N_SAMPLES=15` for all three also keeps the coverage feature block
one width, so a model reads all three tiers unchanged.

Tiers are drawn independently (each runs its own seeded sampling), so the
low pool is *not* a subset of the high pool. The `high` tier keeps the exact
parameters and seeds of the original single-tier build, so it reproduces the
same community as before, minus the two write-only files described under
"Write-only artifacts" below. Every tier carries its own suffix
(`Data/work/multibiome_high/`, `multibiome_high__coasm__metaspades__noviral/`),
so no tier's directory name is a prefix of another's.

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

**Everything lives flat in this one directory — no subfolders, and no shared
config file beyond `../env.sh`'s machine-specific values.** (The real-data pipeline is a sibling directory,
[`../Real_Data/`](../Real_Data/README.md), not a subfolder of this one; the two
share no code, only the output contract.) Every numbered script (`0_`–`6_`) is
fully self-contained:
its own small config block at the top (paths + parameters, including its own
copy of the three-line `VB_TIER` selector), no imports from a sibling
module. This includes phase 5's assemble → coverage → labels → graph engine,
which used to be a separate multi-file package (`standardize.py` + ~15
helper modules across `pipeline/`/`lib/` subfolders, itself once referenced
across two filesystems) — it is now one file,
`5_assemble_and_build_graph.py`, trimmed to exactly the recipe this pipeline
runs (see "Phase 5" below for what was dropped and why).

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

**To change a parameter** (sample count, a seed, ...), edit the small config
block at the top of the script that owns it — there is no central
`config.env`. The exception is the machine-specific handful (`CONDA_BIN`,
`ASM_CONDA_BIN`, `THREADS`, `MEMORY_GB`), which every script reads from
`../env.sh` / the environment with a built-in fallback, so those are edited in
one place. Anything else that is repeated verbatim across scripts (paths, the
`VB_TIER` selector itself) is duplicated by design, since each script must run
standalone; if you change one of those, grep for it across the directory and
update every copy.

**To change a tier** (or add a fourth), the per-tier values live in exactly
three places: `PER_BIOME` in `1_build_genome_pool.py`, `STRAIN_FRACTION` in
`2_expand_strains.py`, and `READS_PER_SAMPLE` in `4_simulate_reads.sh`.
Everything else only maps the tier name to a directory suffix
(`low` → `_low`, `medium` → `_medium`, `high` → `_high`). If you change a
genome count or strain fraction, re-derive `READS_PER_SAMPLE` from the new
reference size at roughly 12,000 reads/Mb, or the depth invariant breaks.

## Phase 5 — what it is and what got trimmed

`5_assemble_and_build_graph.py` is a from-scratch consolidation of what used
to be an external general-purpose tool (`standardize.py`) plus ~15 helper
modules (`tools.py`, `detect.py`, `assemble.py`, `coverage.py`, `labels.py`,
`graph.py`, `_coverage.py`, `_labels.py`, `_graph_utils.py`,
`build_graph.py`, `viral_features.py`, and others). That tool supported many
configurations (per-sample assembly, on-the-fly read simulation, viral-ID
filtering, interactive prompts, strobealign coverage, MEGAHIT/metaFlye/Raven
assembly) because it was shared across every dataset in the wider project.
This pipeline only ever invoked it one way — real pre-simulated reads,
metaSPAdes co-assembly (`--only-assembler`), no viral-ID filtering,
minimap2 coverage, reference-based labels, `--non-interactive` — so
everything else was dropped rather than carried along dead:

- **Kept, unchanged in substance:** tool resolution (`Toolbox`, tries
  `--conda-bin` then `PATH`, both `name` and `name.py`), FASTQ
  detection/pairing, metaSPAdes co-assembly + its resumable input-signature
  cache, seqkit length filter, minimap2 coverage mapping + pysam depth
  computation, reference-alignment labelling, TNF+coverage graph
  construction (canonical k-mer features, GFA→contig edge projection,
  the oriented assembly heterograph), mmseqs protein-cluster features,
  work-directory cleanup.
- **Merged in-process:** the three steps that used to shell out to a
  separate copied script (`coverage.py` → `_coverage.py`, `labels.py` →
  `_labels.py`, `graph.py` → `build_graph.py`) now call a plain Python
  function directly — no subprocess, no `config.json` round-trip through
  disk (a `config.json` is still written as a provenance record, since
  some downstream tooling expects it, but nothing re-reads it back).
- **Dropped (dead code on this recipe):** per-sample/individual assembly
  mode, on-the-fly read simulation, viral-ID tools (geNomad/VirSorter2/
  CheckV), interactive Q&A prompts, strobealign coverage, MEGAHIT/metaFlye/
  Raven assembly, sample-ID prefixing (only needed for individual mode),
  a few unused "compatibility view" helper functions.
- **Simplified:** no argparse config surface beyond `--force`/`--dry-run` —
  every other value is a constant at the top of the file, and the tier comes
  from `VB_TIER`, since this script only ever runs one recipe.

If you need one of the dropped modes back, the last version before this
consolidation is recoverable from git/session history; re-adding a mode is
easier by pulling in just that piece than by starting over.

### Write-only artifacts that are no longer produced

Two outputs were written by every run and read by nothing — not by a later
phase, not by the acceptance check, not by any notebook. They are gone:

- `label_map.json` (`<dataset>/`): the genome-label → integer-class mapping.
  The mapping is already recoverable from `contig_metadata.tsv`'s
  `genome_label` column alongside `viral_graph.pt`'s `y`, which is what
  every consumer actually uses.
- `scaffolds.paths` (work dir): SPAdes' scaffold-level path file. Only
  `contigs.paths` is read (`build_graph()` projects GFA links through it);
  the scaffold variant was copied out of the assembly directory and kept
  through cleanup for nothing.

Everything else phase 5 and 6 emit does have a consumer, including the two
that look like pure provenance: `config.json` supplies the TNF/coverage
feature-layout offsets to the model and analysis notebooks, and
`manifest.json` is read by `dataset_analysis.ipynb` and updated in place by
phase 6.

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

## Notes

- Every generative step is seeded (`POOL_SEED`, `STRAIN_SEED`,
  `ABUNDANCE_SEED` all 42, in scripts 1/2/3; reads use
  `READS_SEED_BASE + sample_index` in script 4), so a clean `--force` run
  is deterministic given identical inputs and tool versions. The seeds are
  shared across tiers but the draw sizes differ, so the three communities
  are independent samples rather than nested subsets.
- `5_assemble_and_build_graph.py`'s assembly step caches an input signature
  (`standardize_work/<run>/assembly_inputs.json`) and skips re-running
  metaSPAdes if it's unchanged — re-invoking after a crash in a *later*
  stage (coverage/labels/graph/viralfeat) does not redo the ~1.5h assembly.
  Coverage mapping has its own per-BAM freshness check, but note that the
  length-filter step always re-runs seqkit and rewrites
  `contigs_filt.fasta`, which bumps its mtime — so a resume after phase 5
  typically re-maps all samples (a few minutes with minimap2) even though
  it skips the assembly itself.
- `6_acceptance_check.py` needs the same `viralbin` python
  (torch/pandas/networkx) as phase 5 — run it as
  `"$VB_CONDA_BIN/python" 6_acceptance_check.py` (after `. ../env.sh`) if your
  default `python3` doesn't have those installed.
- `6_acceptance_check.py`'s verdict is a *share* of the tier's pool species
  that end up multi-contig (≥10% → "STRONG benchmark", ≥4% → "usable"), not
  the absolute counts the single-tier version used, so the same test means
  the same thing on a 50-genome and a 998-genome community.
- The `low` tier has only 2 multi-strain species by construction (5% of 50,
  rounded), so its strain tier is close to trivial. That is the intended
  meaning of "low complexity" here, not a bug — use `medium`/`high` when the
  question is specifically about strain separation.
