#!/usr/bin/env bash
set -euo pipefail
COCONET=${COCONET_CMD:-/home/hashini/miniforge3/envs/coconet_env/bin/coconet}
if [[ -z "${BAM_FILES:-}" ]]; then
  echo "CoCoNet requires BAM_FILES (sorted BAM files for the shared cohort)." >&2
  exit 4
fi
mapfile -t BAMS <<< "$BAM_FILES"
"$COCONET" run --fasta "$COHORT_FASTA" --bam "${BAMS[@]}" \
  --output "$TOOL_OUT" --min-ctg-len 2048 --min-prevalence 1 -t "$THREADS"
