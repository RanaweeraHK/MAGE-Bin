#!/usr/bin/env bash
set -euo pipefail
CONCOCT=${CONCOCT_CMD:-/home/hashini/miniforge3/envs/concoct_env/bin/concoct}
mkdir -p "$TOOL_OUT"
"$CONCOCT" --composition_file "$COHORT_FASTA" --coverage_file "$COVERAGE_TSV" \
  --basename "$TOOL_OUT/" --threads "$THREADS" --length_threshold 2000 --seed 1
