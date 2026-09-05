# External-tool benchmark

This directory benchmarks viral binners and metagenomic binners against the
same **contigs >=2,000 bp** cohort. Methods that use
single-copy marker genes during binning are excluded.

The denominator is fixed before any tool runs. Tools that accept 2,000 bp are
configured at 2,000 bp. A stricter tool is not forced below its supported
limit: CoCoNet retains its `>2,048 bp` rule, its excluded 2,000--2,048 bp
contigs remain unassigned, and `effective_min_contig_bp=2048` records this in
the CSV. This avoids the misleading alternative of giving each method a
different evaluation denominator.

## Included methods

| Method | Category | Marker-gene policy | Effective minimum |
|---|---|---|---:|
| CoCoNet | viral binner | no universal single-copy markers | 2,048 bp (hard tool limit) |
| vRhyme | viral binner | protein redundancy, not universal single-copy markers | 2,000 bp |
| VAMB | metagenomic binner | core VAMB only; no marker reclustering | 2,000 bp |
| CoCoBin | graph metagenomic binner | no single-copy-marker stage | 2,000 bp |
| CLMB | contrastive-learning metagenomic binner | no single-copy-marker stage | 2,000 bp |
| MetaBAT2 | metagenomic binner | marker-free core | 2,000 bp |
| CONCOCT | metagenomic binner | marker-free core | 2,000 bp |


## One command

Run every complete dataset discovered under `Data/Processed_data/`, generating
the real read-alignment inputs required by viral binners when absent:

```bash
python Other_Tools/benchmark.py --prepare-bams
```

Run one new dataset or a subset of tools:

```bash
python Other_Tools/benchmark.py \
  --dataset multibiome_low__coasm__metaspades__noviral \
  --tools vamb cocobin clmb metabat2
```

CoCoNet and vRhyme need cohort BAMs (vRhyme also uses measured coverage
variance). The command above generates and retains the BAMs; later runs reuse
a complete set automatically. To run these tools, use:

```bash
python Other_Tools/benchmark.py \
  --tools coconet vrhyme --prepare-bams
```

The command discovers newly completed dataset directories automatically. It
appends only unseen `(dataset, cohort, level, method)` rows to
`Outputs/comparison_simulated.csv`; previous datasets are preserved and their
already-completed methods are skipped before execution. A failed or missing
status is retried on the next command and only its matching key is replaced;
unrelated rows remain intact. Use
`--rerun` only when deliberately replacing matching rows after rerunning tools.
Use `--refresh-csv` to recompute matching rows from retained raw assignments
without rerunning a tool.


## Metrics and CSV schema

Results are reported independently at strain and parent-virus species levels.
Precision, recall, and F1 use the **GraphBin2 paper's contig-count equations**.
For the bin-by-truth matrix `a[k,s]`, precision is
`sum_k max_s a[k,s] / binned_contigs`; recall is
`sum_s max_k a[k,s] / total_cohort_contigs`, where the denominator includes
unclassified/discarded contigs; F1 is their harmonic mean. Thus precision is
weighted bin purity and recall measures the best-bin recovery of every truth
group while penalizing unbinned contigs. These are not pairwise scores.
`hq_bins` requires at least 90% reference completeness and at most 5%
contamination. ARI/NMI remain secondary clustering metrics and represent every
unclassified contig as a distinct singleton so discarding cannot improve them.

The CSV records the tool version, status, effective cutoff, runtime, exact
wrapper command, cohort size, assigned fractions, GraphBin2 PRF, ARI/NMI, and bin
counts. A failed or unavailable tool still receives a status row, but metric
fields remain empty; failed runs are never presented as zero performance.

Intermediate inputs, logs, and raw tool results are retained under
`Outputs/other_tools/<dataset>/` for auditability.

Public sources: [CoCoBin preprint](https://doi.org/10.1101/2025.08.27.672549),
[CoCoBin code](https://github.com/cucpbioinfo/CoCoBin),
[CLMB paper](https://doi.org/10.1007/978-3-031-04749-7_24), and
[CLMB code](https://github.com/zpf0117b/CLMB).

## Real-data evaluation with CheckV — `benchmark_real.py`

Real metagenomes have no contig-to-genome truth labels, so they must not be
scored with the simulated-data precision/recall equations: without a reference,
precision, recall, F1, ARI and NMI are not approximate, they are **undefined**.
Bin quality is measured after binning, by CheckV.

`Other_Tools/benchmark_real.py` is the sibling of `benchmark.py` — same tools,
same wrappers, same >=2 kb cohort rule — and runs them all in one command:

```bash
python Other_Tools/benchmark_real.py
```

```bash
python Other_Tools/benchmark_real.py --list                 # what would run
python Other_Tools/benchmark_real.py --dry-run
python Other_Tools/benchmark_real.py --dataset lake_water__coasm__metaspades__genomad
python Other_Tools/benchmark_real.py --tools vamb metabat2 concoct
python Other_Tools/benchmark_real.py --max-bin-contigs 500 --checkv-threads 2
```


### Output

| path | contents |
|---|---|
| `Outputs/comparison_real.csv` | one row per dataset+method |
| `Outputs/other_tools_real/<dataset>/` | cohort inputs, per-tool results and logs |
| `Outputs/checkv/<dataset>/` | linked bins, CheckV output, per-bin tables |


