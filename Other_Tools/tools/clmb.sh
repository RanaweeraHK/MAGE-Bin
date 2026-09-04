#!/usr/bin/env bash
set -euo pipefail

COMMIT=a92838924d8731769cc83a9bb436ee03a0fd79b6
if [[ "${1:-}" == "--version" ]]; then
  echo "CLMB commit $COMMIT"
  exit 0
fi

CLMB=${CLMB_CMD:-/home/hashini/miniforge3/envs/clmb_env/bin/vamb}
for required in "$CLMB" "$COHORT_FASTA" "$CLMB_RPKM"; do
  [[ -f "$required" ]] || { echo "Missing required CLMB input: $required" >&2; exit 4; }
done

"$CLMB" --outdir "$TOOL_OUT" --fasta "$COHORT_FASTA" --rpkm "$CLMB_RPKM" \
  --contrastive -m 2000 -p "$THREADS"
