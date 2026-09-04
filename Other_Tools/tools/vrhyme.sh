#!/usr/bin/env bash
set -euo pipefail
VRHYME_ENV_BIN=${VRHYME_ENV_BIN:-/home/hashini/miniforge3/envs/vrhyme_env/bin}
VRHYME=${VRHYME_CMD:-$VRHYME_ENV_BIN/vRhyme}
# Calling the entry point by absolute path does not activate its Conda
# environment. vRhyme also invokes mmseqs2, Prodigal, Mash, and related tools.
export PATH="$VRHYME_ENV_BIN:$PATH"
[[ -n "${BAM_FILES:-}" ]] || {
  echo "vRhyme requires BAM_FILES so coverage variance is measured, not fabricated." >&2
  exit 4
}
mapfile -t BAMS <<< "$BAM_FILES"
SAMTOOLS=${SAMTOOLS_CMD:-/home/hashini/miniforge3/envs/asm_env/bin/samtools}
[[ -x "$SAMTOOLS" ]] || { echo "samtools is required for vRhyme BAM filtering" >&2; exit 5; }

# vRhyme 1.1.0 divides NM by query_length for every BAM record. Long-read
# supplementary/secondary alignments can be hard-clipped with query_length=0,
# so retain only mapped primary records (4 + 256 + 2048 = flag mask 2308).
PRIMARY_BAM_DIR="${TOOL_OUT}.primary_bams"
mkdir -p "$PRIMARY_BAM_DIR"
PRIMARY_BAMS=()
for index in "${!BAMS[@]}"; do
  primary="$PRIMARY_BAM_DIR/sample${index}.primary.bam"
  if [[ ! -s "$primary" || ! -s "$primary.bai" ]]; then
    "$SAMTOOLS" view -@ "$THREADS" -F 2308 -b -o "$primary" "${BAMS[$index]}"
    "$SAMTOOLS" index -@ "$THREADS" "$primary"
  fi
  PRIMARY_BAMS+=("$primary")
done

"$VRHYME" -i "$COHORT_FASTA" -b "${PRIMARY_BAMS[@]}" -o "$TOOL_OUT" \
  -t "$THREADS" -l 2000 --verbose
