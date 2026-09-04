#!/usr/bin/env bash
set -euo pipefail
METABAT2=${METABAT2_CMD:-/home/hashini/miniforge3/envs/metabat2_env/bin/metabat2}
mkdir -p "$TOOL_OUT/bins"
"$METABAT2" -i "$COHORT_FASTA" -a "$METABAT_DEPTH" --cvExt \
  -o "$TOOL_OUT/bins/bin" -m 2000 -s 0 -t "$THREADS" --seed 1
