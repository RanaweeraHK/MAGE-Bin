# Mock-community preprocessing pipeline

Turns **viral mock communities** — real sequencing runs of a community that was
*built from known genomes* — into model-ready datasets under
`Data/Processed_data/`, holding the same artifact set as the other two
pipelines (`viral_graph.pt`, `assembly_heterograph.pt`, `contig_metadata.tsv`,
`protein_clusters.tsv`, `viral_features.tsv`, `config.json`, `manifest.json`).

It exists because of one gap. The simulated benchmark knows every contig's
source genome, but the reads came from a simulator. The real datasets have real
reads, but no contig can ever be labelled, so bins are scored after the fact
with CheckV. A mock has **both**: real sequencing, and a truth that can be
recovered by aligning contigs back to the genomes that went into the tube.

|  | [`Simulated_Data/`](../Simulated_Data/README.md) | [`Real_Data/`](../Real_Data/README.md) | `Mock_Data/` (here) |
|---|---|---|---|
| Reads | simulated | real | **real** |
| Ground truth | yes, from the simulator | none | **yes, from reference alignment** |
| Evaluation | precision/recall/F1, ARI/NMI, HQ bins | CheckV on linked bins | **both** |
| Unlabelled contigs | none | all | some — see [Reading a score](#reading-a-score) |

---

## The datasets

Four registry rows over two communities.

| dataset | community | libraries | reads | assembler | truth genomes | download |
|---|---|---|---|---|---|---|
| `dsmz_dsrna` | DSMZ plant-virus community | 22 | short | MEGAHIT | 115 molecules / 68 agents | 4.4 GB |
| `dsmz_vana` | DSMZ plant-virus community | 22 | short | MEGAHIT | 115 molecules / 68 agents | 4.4 GB |
| `phage_mock_illumina` | 15-phage mock | 3 | short | MEGAHIT | 15 | 1.1 GB |
| `phage_mock_nanopore` | 15-phage mock | 7 | long | metaFlye | 15 | 5.3 GB |

Short rows assemble with **MEGAHIT** rather than metaSPAdes: it is far lighter
on memory, which matters on a machine where the WSL disk shares a full `C:`.
The contig-level GFA comes from `megahit_toolkit contig2fastg` on the largest
k, the same path `../Real_Data`'s `lake_water_megahit` row uses.

### 1. DSMZ plant-virus synthetic community

A community assembled at DSMZ from quality-controlled, lyophilised infected
plant material, then sequenced by two virome protocols.

| | |
|---|---|
| Paper | *Benchmarking of virome metagenomic analysis approaches using a large, 60+ members, viral synthetic community*, J Virol 2023, [10.1128/jvi.01300-23](https://doi.org/10.1128/jvi.01300-23) |
| Reads | Recherche Data Gouv (INRAE), [doi:10.57745/EMM5EQ](https://doi.org/10.57745/EMM5EQ) — 47 files, 8.8 GB |
| References | [doi:10.57745/T4UYPC](https://doi.org/10.57745/T4UYPC) — `Reference genomes.fas` |
| Licence | Etalab Open Licence 2.0 (CC-BY compatible) |
| Composition | 68 viral agents over **115 molecules**, 21 families, 63 genera |
| Reference sizes | 230 bp – 15,469 bp, median 4,268 bp, 571,460 bp total |
| Read format | FASTA, unpaired, ~113 bp, already quality-trimmed, headers all `>No_name` |

The 44 libraries are **nested sub-communities** drawn from one pool, twice —
once per protocol:

| sub-community | libraries per protocol | viruses |
|---|---|---|
| P5 | 12 | 5 |
| P10 | 6 | 10 |
| P20 | 3 | 20 |
| P60 | 1 | 60 |

That nesting is why this is the stronger of the two datasets: each library is a
*different* subset of the pool, so the 22 coverage columns carry real
presence/absence structure instead of 22 replicates of one community. The three
control libraries in the deposit are not registered.

**What it tests.** 115 molecules over 68 agents means the truth asks a binner to
put the segments of a multipartite virus in one bin — the case coverage-plus-
composition binning is actually for, and one no other dataset in this repo can
score. `dsmz_dsrna` and `dsmz_vana` are the same communities through different
wet-lab chemistry, so the pair separates what the protocol makes recoverable
from what the binner does with it.

**Caveats.** Plant virus molecules are small: with phase 2's 1,000 bp length
filter, the satellites and small DNA components at the bottom of the size range
drop out before binning. And short genomes assemble into few contigs each, so
the assembly graph carries less than it does on a fragmented metagenome.

### 2. 15-phage mock community

Genomic DNA from 15 previously sequenced phages, combined at known copy numbers
and sequenced on three platforms.

| | |
|---|---|
| Paper | *The long and short of it: benchmarking viromics using Illumina, Nanopore and PacBio sequencing technologies*, Microb Genom 2024, [10.1099/mgen.0.001198](https://doi.org/10.1099/mgen.0.001198) |
| Reads | ENA study [PRJEB56639](https://www.ebi.ac.uk/ena/browser/view/PRJEB56639) |
| Assemblies | FigShare [10.25392/leicester.data.21346935](https://doi.org/10.25392/leicester.data.21346935) (not used here) |
| Runs | 3 Illumina MiSeq paired (2.4 Gbp), 7 ONT MinION (5.3 Gbp), 3 PacBio Sequel (1.1 Gbp, **not registered**) |
| Composition | 15 phages, 44.5–320 kb, 28.4–60.8 % GC |

The community, from Table S1 of the paper (copied into
`Data/Mock_dataset/phage_mock_15/references/`):

| phage | genome | length (bp) | GC % | genome copies |
|---|---|---|---|---|
| SLUR29 | Escherichia phage vB_Eco_SLUR29 | 48,593 | 44.7 | 684,329,546 |
| PARMAL1 | PARMAL1 | 44,565 | 57.6 | 200,736,667 |
| J2 | Escherichia phage vB_Eco_mar002J2 | 50,343 | 44.4 | 143,127,748 |
| HP1 | Bacteriophage HP1 | 32,355 | 40.0 | 119,757,670 |
| SM033 | Vibriophage vB_Vpa_sm033 | 320,253 | 43.2 | 82,540,671 |
| J3 | Escherichia phage vB_Eco_mar003J3 | 115,471 | 39.8 | 79,342,556 |
| SWAN | Escherichia phage vB_EcoS_swan01 | 50,865 | 44.7 | 76,036,616 |
| KUW1 | KUW1 | 44,509 | 60.8 | 72,995,152 |
| VP1 | Bacteriophage DSS3_PM1 | 70,044 | 47.3 | 62,567,273 |
| J1 | Escherichia phage vB_Eco_mar001J1 | 50,343 | 44.4 | 53,672,906 |
| PHAGE1 | Escherichia phage vB_Eco_mar005P1 | 167,773 | 37.7 | 53,450,592 |
| SM032 | Vibriophage vB_VpaS_sm032 | 79,660 | 45.7 | 52,465,265 |
| phix174 | Coliphage phi-X174 | 5,386 | 44.8 | 5,164,751 |
| SRSM4 | Synechococcus phage S-RSM4 | 194,454 | 41.1 | 465,530 |
| CDMH1 | Clostridium phage CDMH1 | 54,279 | 28.4 | 169,000 |

Copy numbers span **four orders of magnitude**, so recall across the abundance
range is measurable rather than assumed.

**What it tests.** The community was deliberately built with near-identical
members (Table S2):

| pair | similarity |
|---|---|
| J1 – J2 | 100 % |
| SLUR29 – SWAN | 87.0 % |
| J1/J2 – SWAN | 83.4 % |
| J1/J2 – SLUR29 | 81.8 % |
| KUW1 – PARMAL1 | 50.5 % |

That is strain-level microdiversity on real reads — the case the simulated
benchmark probes and no real dataset here can. J1 and J2 are *not separable*;
phase 3 merges them into one truth genome rather than inventing an answer (see
[Ambiguity](#ambiguity)).

**Caveats.** All runs of a platform sequence the **same** community, so the
coverage columns are replicates and carry almost no differential signal —
composition and the assembly graph do the work. And 15 genomes is a small
problem: read it as a purity and over-merging test, not as an F1 benchmark.
PacBio runs are excluded because they are CLR and would need Flye's
`--pacbio-raw`, while this pipeline uses `--nano-hq`.

---

## Before you run anything

Setup is shared with the sibling pipelines and lives one level up — see
[`../README.md`](../README.md). In short:

```bash
cd /path/to/Viralbinning/Dataset_Processing
bash setup_environments.sh --only viralbin asm_env megahit_env
bash Mock_Data/0_check_environment.sh
```

This pipeline needs **three** environments, not five: `viralbin` (python,
mmseqs, flye), `asm_env` (minimap2, samtools, seqkit) and `megahit_env`
(megahit, megahit_toolkit). No `sra_env` —
phase 1 pulls FASTQ straight off the ENA mirror. No geNomad and no CheckV
database are required, since every registered mock is pure viral material.

## One command

```bash
cd /path/to/Viralbinning/Dataset_Processing/Mock_Data
bash run_all.sh
```

```bash
bash run_all.sh --list                          # what would be built, then exit
bash run_all.sh --dataset phage_mock_illumina   # just this one
bash run_all.sh --force                         # rebuild from scratch
bash run_all.sh --from 3                        # resume at a given phase
bash run_all.sh --discard-downloads             # delete provider files once converted
```

Phases skip work already done, so re-running after an interruption is safe.
Downloads are tens of GB and co-assembly is long — run under `tmux`/`nohup`.

## Scripts, in order

| Script | What it does | Needs |
|---|---|---|
| `0_check_environment.sh` | Preflight: registry, tools, download endpoints, disk. | — |
| `1_get_inputs.py` | Downloads the libraries from ENA or the INRAE deposit, converts FASTA reads to FASTQ, and builds the reference set. | network |
| `2_assemble_and_build_graph.py` | Co-assembly → length filter → coverage → PyG graph + viral features. Runs `../Real_Data/2_assemble_and_build_graph.py` against this registry. | megahit/flye, minimap2, samtools, seqkit, mmseqs, torch+PyG |
| `3_label_from_references.py` | **The truth.** Aligns contigs to the reference set and fills `genome_label`. | minimap2 |
| `4_acceptance_check.py` | Is there anything to merge? Runs `../Real_Data/3_acceptance_check.py`. | torch, pandas, networkx |

Each also runs standalone:

```bash
VB_DATASET=phage_mock_illumina python 1_get_inputs.py --dry-run
VB_DATASET=dsmz_dsrna          python 3_label_from_references.py --min-identity 0.92
```

Phases 2 and 4 delegate to the real-dataset pipeline rather than copying it:
a mock is real sequencing, the artifacts must be byte-comparable, and a second
copy of a two-thousand-line script would drift. They set `VB_REGISTRY`,
`VB_WORK_PREFIX=mock` and `VB_DATASET_ROOT=Mock_dataset`; unset, those scripts
behave exactly as they always have.

## Phase 3 — how the truth is built

```
contigs_filt.fasta --minimap2 -cx asm20 -N 50 -p 0.1--> reference_genomes.fasta
```

For every (contig, reference) pair the alignments are summed into *matching
bases* and a *merged covered span* of the contig. The best reference wins, if it
clears both thresholds:

| threshold | default | meaning |
|---|---|---|
| `--min-identity` | 0.90 | matching bases / alignment block length |
| `--min-query-cov` | 0.50 | fraction of the contig covered by that reference |
| `--ambiguity-margin` | 0.02 | how close a rival must be to count as a tie |

Three rules then decide what a score can mean.

### Ambiguity

Genomes that no assembly can separate — J1 and J2 at 100 % identity — are
merged into **one** truth genome, `J1|J2`, consistently across the dataset
(union-find over every tie observed). Picking the higher-scoring reference
would invent a truth and then penalise a binner for not reproducing it. The
merge is printed and recorded in `truth_summary.json`, so the ceiling it implies
is visible.

### Chimeras

A contig covered by two different genomes on *different parts* of its length
(≥25 % disjoint) has no single true genome. It is left unlabelled and counted
separately: scoring it would measure the assembler, not the binner.

### Agent vs molecule

`genome_label` holds the **agent** — the virus — so a binner is asked to co-bin
a segmented genome. `molecule_label` keeps the segment beside it. For the phage
mock the two are identical; for the DSMZ community 115 molecules collapse onto
68 agents. `--level molecule` scores the finer question instead.

### Outputs

| File | Contents |
|---|---|
| `contig_metadata.tsv` | gains `genome_label`, `molecule_label`, `label_identity`, `label_query_cov`, `label_status` |
| `truth_per_reference.tsv` | per reference: how much of it the assembly covers, how many contigs were assigned |
| `truth_summary.json` | counts, merged genomes, thresholds used |
| `manifest.json` | `has_ground_truth` and `labelled` flip to true |

`label_status` is one of `labelled`, `merged`, `chimeric`, `unlabelled`.

## Reading a score

**Score on the labelled subset.** An unlabelled contig is material the
reference set does not explain — reagent contamination, host carry-over, or a
genome that drifted past the thresholds — not a binner error. This is the same
convention the simulated benchmark already uses when it ignores unlabelled
contigs.

Check `truth_summary.json` before trusting a number:

- **low labelled fraction** → the reference set explains little of the assembly;
  loosen `--min-identity`, or accept that the dataset is mostly background.
- **`merged_truth_genomes` non-empty** → some genomes are unseparable by
  construction, and no method can exceed that ceiling.
- **`references_half_covered` well below the community size** → genomes that
  were never assembled cannot be binned, so recall is capped by the assembler.

## CheckV as well

CheckV still works here and is worth running, because it scores bins the way
the real datasets are scored — and on a mock you can compare what CheckV says
to what the truth says:

```bash
VB_REGISTRY=$PWD/datasets.tsv python ../Real_Data/4_checkv_evaluate.py \
    --dataset phage_mock_illumina \
    --bins ../../Outputs/..._assignments/phage_mock_illumina...tsv \
    --method MAGE-Bin
```

## FASTA reads

The DSMZ libraries are distributed as quality-trimmed FASTA. Phase 1 rewrites
them as FASTQ with a flat placeholder quality, because phase 2 detects reads by
FASTQ extension. Nothing reads the placeholder: MEGAHIT does no error
correction, and minimap2 ignores qualities.

Both copies cost disk, so `--discard-downloads` deletes each provider file once
it has been converted:

```bash
bash run_all.sh --discard-downloads --dataset dsmz_dsrna
```

## Where things land

```
Data/Mock_dataset/<study>/downloads/            <- provider files as downloaded
Data/Mock_dataset/<study>/references/<set>/     <- reference_genomes.fasta
                                                   reference_map.tsv
Data/work/mock_<dataset>/reads/                 <- normalised FASTQ, one per library
Data/work/mock_<dataset>/standardize_work/<run>/ <- assembly scratch, + the
                                                   contigs_vs_references.paf
Data/Processed_data/<run>/                      <- FINAL dataset (kept)
Outputs/checkv/<run>/                           <- if CheckV is run
```

where `<run>` is `<dataset>__coasm__<assembler>__noviral`.

