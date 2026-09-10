# Real-dataset preprocessing pipeline

Turns real sequencing runs into model-ready datasets under
`Data/Processed_data/`, holding the same artifact set as the simulated
benchmark (`viral_graph.pt`, `assembly_heterograph.pt`,
`contig_metadata.tsv`, `protein_clusters.tsv`, `viral_features.tsv`,
`config.json`, `manifest.json`) — with one difference that shapes
everything else: **there is no ground truth.**

No reference generated these reads, so no contig can be labelled with the
genome it came from. Precision, recall, F1, ARI and HQ-bin counts — every
metric the simulated benchmark reports — are undefined here. Bin quality is
instead measured *after* binning, by CheckV, in `4_checkv_evaluate.py`.


## Two kinds of dataset: `source = sra` or `source = assembly`

Not every dataset starts from reads. The registry's `source` column picks the
path, and it is the only thing you change:

| | `sra` | `assembly` |
|---|---|---|
| Input | SRA accessions | a directory holding `contigs.fasta` + GFA + path table |
| Phase 1 | `prefetch` + `fasterq-dump` → `reads/` | **imports** the assembly into the work dir |
| Phase 2 assembly | metaSPAdes / metaFlye | **skipped** — the assembly is given |
| Phase 2 coverage | minimap2 + pysam | adopted from `import_from`, or it fails |
| Phase 2 geNomad | run (bulk metagenome) | adopted from `import_from` if present |
| Phases 3, 4 | unchanged | unchanged |

### `source = assembly`, and what "adopting" means

A pre-built assembly carries **no reads**, and coverage cannot be conjured
without them. The `import_from` column therefore points at an earlier work
directory whose finished artifacts are taken as-is:

```
import_from/coverage.csv        -> the coverage matrix (defines n_samples)
import_from/contigs_filt.fasta  -> the contig set that build settled on
import_from/viral/genomad/      -> the whole geNomad output tree
```

## Before you run anything

Setup is shared with the sibling simulated pipeline and lives one level up —
see [`../README.md`](../README.md) ("Setup — do this once"), which lists every
environment and every `env.sh` variable. In short:

```bash
cd /path/to/Viralbinning/Dataset_Processing
bash setup_environments.sh              # all five envs + Dataset_Processing/env.sh
bash setup_environments.sh --databases  # optional: fetch the CheckV + geNomad DBs
bash Real_Data/0_check_environment.sh   # preflight: registry, tools, databases
```

This pipeline uses all five environments (`viralbin`, `asm_env`, `sra_env`,
and — only for `viral_enriched=no` and `assembler=megahit` rows —
`genomad_env` and `megahit_env`), plus two databases whose locations are
settings rather than downloads the pipeline makes for you: `CHECKVDB` and
`GENOMAD_DB`. All of it is `Dataset_Processing/env.sh`, which is gitignored,
so no path to your machine reaches the repo. Every script here reads that one
file — the shell phases source it, the Python phases parse it — and real
environment variables win over it:

```bash
VB_THREADS=32 bash run_all.sh --dataset lake_water
```

Which dataset gets built is *not* a machine setting: it comes from
`datasets.tsv` and the command line, as below.

## One command

```bash
cd /path/to/Viralbinning/Dataset_Processing/Real_Data
bash run_all.sh
```

Phases run in order, skipping any phase whose output already exists, so it
is safe to re-run after an interruption.

```bash
bash run_all.sh --list                      # show what would be built, exit
bash run_all.sh --dataset aloha_illumina    # just this one
bash run_all.sh --force                     # rebuild from scratch
bash run_all.sh --from 2                    # resume starting at phase 2
```

Phase 1 downloads tens of GB and phase 2's co-assembly is long — run under
`nohup`/`tmux`/`screen`, not in a session that might close.

## Scripts, in order

| Script | What it does | Needs |
|---|---|---|
| `0_check_environment.sh` | Preflight: registry validity, tools, CheckV/geNomad databases. Reports a missing geNomad as a *warning* while every registered dataset is a virome. | — |
| `1_get_inputs.sh` | `source=sra`: `prefetch` each accession into `Data/Real_dataset/<study>/`, `fasterq-dump --split-files` it, gzip into the dataset's `reads/`. `source=assembly`: import the assembly and adopt whatever `import_from` names. Resumable either way. | `sra-tools` (sra only) |
| `2_assemble_and_build_graph.py` | Co-assembly → length filter → coverage → viral ID → PyG graph + viral features. No labels. Any stage already in the work dir is `[adopted]` when `source=assembly`. | `metaspades.py`/`flye`, `minimap2`, `samtools`, `seqkit`, `mmseqs`, geNomad (bulk only), python w/ torch+PyG+pysam+pyrodigal-gv+sklearn |
| `3_acceptance_check.py` | Reference-free sanity: is there differential coverage, an assembly graph, and shared protein content — i.e. anything to merge? | python w/ torch+pandas+networkx |
| `4_checkv_evaluate.py` | **After binning**: bins → linked pseudo-genomes → CheckV → the metrics table. Not part of `run_all.sh`. | `checkv` + checkv-db |

Every script also runs standalone, selecting its dataset from `VB_DATASET`:

```bash
VB_DATASET=aloha_illumina bash 1_get_inputs.sh
VB_DATASET=lake_water     bash 1_get_inputs.sh
VB_DATASET=lake_water     python 2_assemble_and_build_graph.py --dry-run
```

## Phase 4 — CheckV evaluation

Not in `run_all.sh`, because it scores **bins**, which only exist once a
binner has run:

```bash
python 4_checkv_evaluate.py \
    --dataset aloha_illumina \
    --bins ../../Outputs/sigvib_assignments/aloha_illumina...tsv \
    --method SIGVIB-A
```

Several methods on one dataset compare directly (each is one row of
`Outputs/checkv/<run>/checkv_metrics.csv`):

### The linking step, and why it matters

**CheckV scores one sequence at a time.** Hand it a bin's contigs as separate
sequences and it returns one row per contig — a genome correctly recovered as
5 contigs is scored as 5 fragments at ~20% completeness each, and a binner
that reassembled it perfectly is indistinguishable from one that did nothing.

vRhyme's answer, in its `link_bin_sequences.py` helper, is to concatenate
each bin's scaffolds into a **single sequence separated by runs of N**, so
CheckV sees one pseudo-genome per bin and its completeness and contamination
describe the *bin*. This script does the same thing, at vRhyme's default
linker of **1,500 Ns**:

    bin_3  =  contigA + N×1500 + contigB + N×1500 + contigC

The linker must be long enough that gene callers and alignment do not read
through the junction and invent a gene spanning two contigs (vRhyme warns
below 1,000). Contig order within a bin is arbitrary and claims nothing about
the genome's real layout — CheckV's completeness comes from gene content and
AAI against its reference genomes, not from synteny, so the ordering does not
affect the estimate.

### Metrics

Per method, written to `Outputs/checkv/<run>/checkv_metrics.csv`, with a
`per_bin.tsv` beside it:

| Metric | Definition |
|---|---|
| Predicted vMAGs | bins evaluated |
| multi-contig vMAGs | bins with ≥2 contigs — a 1-contig "bin" is the assembler's result, not the binner's, so both denominators are reported |
| CheckV Complete | `checkv_quality == Complete` |
| CheckV HQ | `checkv_quality == High-quality` (>90% complete) |
| CheckV MQ | `checkv_quality == Medium-quality` (50–90%) |
| MQ or better | Complete + HQ + MQ |
| Mean / median completeness | over bins CheckV could place; `Not-determined` bins are counted separately and **excluded** from the averages |
| Mean / median contamination, bins >5% contaminated | 5% is the MIUViG high-quality ceiling |
| Largest / median bin size, contigs binned | so a method that scores well by emitting one enormous bin is visible rather than hidden |

Quality tiers are CheckV's own; no threshold is re-derived here.

## Where things land

```
Data/Real_dataset/<study>/<accession>/       <- prefetch's .sra cache (source=sra)
Data/Real_dataset/<assembly_dir>/            <- a pre-built assembly (source=assembly)
Data/work/real_<dataset>/
    reads/<accession>_{1,2}.fastq.gz         <- phase 1 output (source=sra)
    standardize_work/<run>/                  <- phase 1 imports here (source=assembly),
                                                phase 2 works here; cleaned down to
                                                contigs_filt.fasta, contigs_viral.fasta,
                                                coverage.csv, viral_contigs.txt,
                                                viral/genomad/, assembly_graph.gfa
                                                + path table
Data/Processed_data/<run>/                   <- FINAL dataset (kept)
Outputs/checkv/<run>/<method>/               <- phase 4: linked bins, CheckV, per_bin.tsv
Outputs/checkv/<run>/checkv_metrics.csv      <- one row per method
```

where `<run>` is `<dataset>__coasm__<assembler>__<noviral|genomad>`.




