#!/usr/bin/env bash
set -euo pipefail
VAMB=${VAMB_CMD:-/home/hashini/miniforge3/envs/vamb_env/bin/vamb}
"$VAMB" bin default --outdir "$TOOL_OUT" --fasta "$COHORT_FASTA" \
  --abundance_tsv "$COVERAGE_TSV" -m 2000 -p "$THREADS" -o --seed 0
