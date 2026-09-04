#!/usr/bin/env bash
set -euo pipefail

COMMIT=594fed94cdf3fd2d9e1420eb8c667f3f9eed32ca
if [[ "${1:-}" == "--version" ]]; then
  echo "CoCoBin commit $COMMIT"
  exit 0
fi

COCOBIN_PYTHON=${COCOBIN_PYTHON:-/home/hashini/miniforge3/envs/viralbin/bin/python}
COCOBIN_SCRIPT=${COCOBIN_SCRIPT:-/home/hashini/Viralbinning/Other_Tools/vendor/CoCoBin/Binning_project/Binning_Main.py}
for required in "$COCOBIN_PYTHON" "$COCOBIN_SCRIPT" "$COCOBIN_FASTA" "$COCOBIN_KMER" "$COCOBIN_MAPPING"; do
  [[ -f "$required" ]] || { echo "Missing required CoCoBin input: $required" >&2; exit 4; }
done

mkdir -p "$TOOL_OUT"
"$COCOBIN_PYTHON" "$COCOBIN_SCRIPT" "$COCOBIN_FASTA" "$COCOBIN_KMER" \
  -o "$TOOL_OUT/raw_assignments.csv"

awk -F',' 'BEGIN {OFS="\t"}
  NR == FNR {if (FNR > 1) {split($0, fields, "\t"); mapping[fields[1]] = fields[2]}; next}
  FNR == 1 {print "bin", "contig"; next}
  {gsub(/\r/, "", $2); if ($2 in mapping) print $1, mapping[$2]}
' "$COCOBIN_MAPPING" "$TOOL_OUT/raw_assignments.csv" > "$TOOL_OUT/assignments.tsv"
