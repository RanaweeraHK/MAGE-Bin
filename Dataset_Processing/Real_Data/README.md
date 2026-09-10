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

## What is different from `../Simulated_Data/`

| | Simulated | Real |
|---|---|---|
| Where reads come from | InSilicoSeq / Badread simulate them from a known genome pool | SRA accessions, downloaded (phase 1) |
| Contig labels | minimap2 contigs vs. the strain reference → `labels.csv` | **none** — `y = -1`, `genome_label` empty |
| Which contigs are viral | known by construction (the pool is viral) | virome library → all of them; bulk metagenome → geNomad |
| Evaluation | GraphBin2 P/R/F1, ARI/NMI, HQ bins vs. truth | CheckV completeness / contamination on linked bins |
| Selected by | tier + simulator + assembler (`low:iss:metaspades`) | a row in `datasets.tsv` |
| `config.json` says | `has_ground_truth` absent (truth exists) | `"has_ground_truth": false` |

Everything else — assembly, the length filter, minimap2 coverage, the
GFA→contig graph projection, the oriented assembly heterograph, ORF calling
and MMseqs2 protein clusters — is the **same code, copied verbatim** from
`Simulated_Data/5_assemble_and_build_graph.py`, so a real dataset and a
simulated one are never processed by two subtly different implementations.

## The datasets: `datasets.tsv`

One row per dataset. This is the only file you edit to add data, and the
only file shared between the scripts (each script parses it in a dozen lines
of stdlib and still runs standalone). Columns are documented in the file's
own header.

Currently registered — both from BioProject **PRJNA529454** / SRP189859,
*"Assembly Free Viral Genomes from Station ALOHA"* (HOT/SCOPE cruises,
DeLong Lab), the viral-TFF metagenomes CoCoNet was evaluated on:

| dataset | source | runs | read type | assembler | viral ID |
|---|---|---|---|---|---|
| `aloha_illumina` | `sra` | SRR8811962, SRR8811963, SRR10378148 | short (NextSeq 500, paired) | metaSPAdes | none (virome) |
| `aloha_nanopore` | `sra` | SRR8811960, SRR8811961, SRR8811964, SRR10881655, SRR10378147 | long (GridION, single) | metaFlye | none (virome) |
| `lake_water` | `assembly` | SRR15557808–SRR15557813 | short (HiSeq 2500, paired) | metaSPAdes (pre-built) | geNomad (adopted) |

The BioProject's 8 runs split cleanly into 3 Illumina and 5 Nanopore, which
is why there are two datasets rather than one: **one sample per accession,
one coverage column per sample**, and mixing read types in a single
co-assembly is not a thing either assembler does. Trim the accession lists
if you want a cheaper first pass — but note that the accession list *is* the
differential-coverage signal the binner sees, and 3 samples is already thin.

### On `badread`

`badread` is a long-read **simulator**, so it has no role on real data —
there is nothing to simulate. The short/long split it provides in the
simulated pipeline is provided here by the runs' own library type, and it
picks the assembler by the identical rule:

    short reads -> metaspades, minimap2 -x sr
    long  reads -> metaflye,   minimap2 -x map-ont

`badread` + `metaflye` therefore stays a simulated-pipeline pairing;
`aloha_nanopore` is the real-data counterpart to that dataset.

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

Phase 1 copies them into `Data/work/real_<dataset>/standardize_work/<run>/`
and derives `viral_contigs.txt` from geNomad's `*_virus_summary.tsv`, so the
adopted call is explicit rather than implied by a filename. Phase 2 then
prints `[adopted]` for every stage it skips, and records what was taken in the
dataset's `manifest.json`. **Nothing is silently reused:** a `source=sra`
dataset never adopts — it always recomputes, so an interrupted run cannot
leave a half-written artifact that a later run trusts.

geNomad's *entire* output tree is copied, not just the summary phase 2 reads.
It is the expensive artifact (hours on a 44k-contig assembly), and keeping the
annotation and provirus tables beside it means the call can be re-examined
later without the original run directory.

### The `lake_water` dataset, concretely

`Data/Real_dataset/Lake_water_metaSPAdes/` holds the metaSPAdes assembly
distributed with Phables — 344,873 contigs, a 309 MB GFA and `contigs.paths`.
It is a bulk lake metagenome, so `viral_enriched=no` and geNomad's calls define
the viral set. Running it is one command:

```bash
bash run_all.sh --dataset lake_water
```

which does:

```
phase 1  imports contigs.fasta / assembly_graph.gfa / contigs.paths from
         Data/Real_dataset/Lake_water_metaSPAdes/, and adopts
         coverage.csv (6 samples), contigs_filt.fasta and viral/genomad/
         from the earlier build named in import_from
phase 2  [adopted] assembly, filtered contigs, coverage, geNomad calls
         -> builds only the graph, contig metadata and protein clusters
phase 3  reference-free acceptance check
```

**Two honest caveats about this dataset**, both consequences of adopting
someone else's build rather than of the pipeline:

1. **Coverage was measured against the viral subset, not every contig.** The
   earlier build ran geNomad first and then mapped reads to the 18,716 viral
   contigs. This pipeline does it the other way round (see "Why coverage is
   computed before viral identification" below) precisely because mapping to a
   pre-filtered subset pushes reads off discarded contigs onto whichever viral
   contig they match next best. The adopted matrix carries that bias; the only
   fix is to re-fetch the six runs and re-map against all contigs.
2. **Read length is declared, not measured.** The reads are no longer on disk,
   so `mean_read_length` in the registry is the SRA-reported average for the six
   runs (595.3 bp per paired spot → 297.7 per read). Left blank, the model
   would default to 126 bp and overstate read counts by 2.4x, making its
   abundance posterior that much over-confident. `config.json` records
   `read_length_source: "registry"` so the distinction survives.

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

## Why coverage is computed before viral identification

Coverage is measured against **every** length-filtered contig, and geNomad's
viral subset is taken afterwards. Mapping to a pre-filtered subset would push
reads from discarded contigs onto whichever viral contig they matched next
best, inflating coverage exactly where the binner is most sensitive to it.
The full contig set stays in the work directory as `contigs_filt.fasta`; the
graph is built on `contigs_viral.fasta`, and `viral_contigs.txt` records the
decision.

For a virome (`viral_enriched=yes`) no predictor runs at all: every assembled
contig is taken as viral. Running a virus classifier over VLP-enriched
contigs would discard real but hard-to-call phage and silently shrink the
benchmark. The run name keeps the `__noviral` suffix, which means exactly
what it means in the simulated pipeline — *no viral-ID tool was applied*.

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

```bash
python 4_checkv_evaluate.py --dataset aloha_illumina \
    --bins sigvib.tsv vrhyme.tsv vamb.tsv \
    --method SIGVIB-A vRhyme VAMB
```

It accepts either a `contig_id`/`bin_id` header (what the model notebook
writes) or a bare two-column table (what most binners emit), and takes
`--contigs <fasta>` for anything outside the registry — including simulated
datasets, where it is a useful reference-free second opinion alongside the
truth-based metrics.

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

## Notes

- **The model notebook skips these datasets automatically.** Its dataset
  discovery filters on `config.json`'s `has_ground_truth`, because every
  metric it computes needs labels. Real datasets are for the CheckV route.
- geNomad's **database is already present** at `$GENOMAD_DB` (see `env.sh`);
  the binary is not installed yet. Neither is needed until a
  `viral_enriched=no` dataset is registered, which is why phase 0 reports it
  as a warning. Install with
  `mamba create -n genomad_env -c conda-forge -c bioconda genomad`.
- `1_get_inputs.sh` never deletes a `.sra.lock`/`.sra.tmp` pair — `prefetch`
  resumes from them, so an interrupted download costs only the run it died on.
- One accession = one sample = one coverage column. Sample order follows
  accession order and is stable across re-runs, so a rebuilt dataset's
  coverage columns line up with the previous build's.
