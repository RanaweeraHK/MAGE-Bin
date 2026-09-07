#!/usr/bin/env bash

# run this file using bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Two dirname levels: this script lives in Dataset_Processing/Simulated_Data/.
REPO_HOME="$(dirname "$(dirname "$SCRIPT_DIR")")"
RAW_DIR="$REPO_HOME/Data/Raw_dataset/multibiome_source_data"

# ---- machine-specific settings: ../env.sh, shared by both pipelines --------
# Created by `bash ../setup_environments.sh` or if wants manually
[ -f "$SCRIPT_DIR/../env.sh" ] && . "$SCRIPT_DIR/../env.sh"
CONDA_BIN="${VB_CONDA_BIN:-$HOME/miniforge3/envs/viralbin/bin}"        # mmseqs, python (torch/PyG/pysam/pyrodigal-gv)
ASM_CONDA_BIN="${VB_ASM_CONDA_BIN:-$HOME/miniforge3/envs/asm_env/bin}" # metaspades.py, minimap2, samtools, seqkit
PY="$CONDA_BIN/python"
[ -x "$PY" ] || PY="python3"

missing=0

check_bin() {
    local name="$1" bin="$2" hint="$3"
    local found=""
    if [ -x "$CONDA_BIN/$bin" ]; then
        found="$CONDA_BIN/$bin"
    elif [ -x "$ASM_CONDA_BIN/$bin" ]; then
        found="$ASM_CONDA_BIN/$bin"
    elif command -v "$bin" >/dev/null 2>&1; then
        found="$(command -v "$bin")"
    fi
    if [ -n "$found" ]; then
        printf '  [ok]   %-12s -> %s\n' "$name" "$found"
    else
        printf '  [MISS] %-12s (looked in %s, %s and PATH)\n' "$name" "$CONDA_BIN" "$ASM_CONDA_BIN"
        printf '         install: %s\n' "$hint"
        missing=$((missing + 1))
    fi
}

# use to construct multibiome viral genome pool.
echo "=== raw data sources (RAW_DIR=$RAW_DIR) ==="
for f in \
    "reference/refseq_viral_raw.fasta" \
    "human_gut_mgv" \
    "freshwater/freshwater_self_circular.fasta" \
    "marine_imgvr/marine_imgvr.fasta" \
    "soil/soil_GSV.fasta" \
    "soil/soil_GSV_quality.csv"; do
    if [ -e "$RAW_DIR/$f" ]; then
        printf '  [ok]   %s\n' "$f"
    else
        printf '  [MISS] %s\n' "$f"
        missing=$((missing + 1))
    fi
done

echo
echo "=== tools (CONDA_BIN=$CONDA_BIN) ==="
REQUIRED_SIMS="${VB_REQUIRED_SIMS:-iss badread}"
REQUIRED_ASMS="${VB_REQUIRED_ASMS:-metaspades metaflye}"
echo "  (checking simulators:$REQUIRED_SIMS assemblers:$REQUIRED_ASMS)"

for sim in $REQUIRED_SIMS; do
    case "$sim" in
        iss) check_bin "iss" "iss" \
            "bash ../setup_environments.sh --only viralbin" ;;
        badread) check_bin "badread" "badread" \
            "bash ../setup_environments.sh --only viralbin" ;;
        *) printf '  [MISS] unknown simulator %s\n' "$sim"; missing=$((missing + 1)) ;;
    esac
done
for asm in $REQUIRED_ASMS; do
    case "$asm" in
        metaspades) check_bin "metaspades" "metaspades.py" \
            "bash ../setup_environments.sh --only asm_env" ;;
        metaflye) check_bin "metaflye" "flye" \
            "bash ../setup_environments.sh --only viralbin" ;;
        *) printf '  [MISS] unknown assembler %s\n' "$asm"; missing=$((missing + 1)) ;;
    esac
done

# always required tools
check_bin "minimap2"     "minimap2"      "bash ../setup_environments.sh --only asm_env"   # map reads back to contigs 
check_bin "samtools"     "samtools"      "bash ../setup_environments.sh --only asm_env"   # process BAM/SAM alignments
check_bin "seqkit"       "seqkit"        "bash ../setup_environments.sh --only asm_env"   # manipulate, filter FASTA/FASTQ
check_bin "mmseqs"       "mmseqs"        "bash ../setup_environments.sh --only viralbin"  # sequence similarity, clustering, genome pool processing

echo
echo "=== python ($PY) ==="
if [ -x "$PY" ]; then
    "$PY" - <<'PY' || missing=$((missing + 1))
import importlib, sys
mods = ["torch", "torch_geometric", "numpy", "pandas", "networkx", "pyrodigal_gv",
        "pysam", "sklearn"]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
        print(f"  [ok]   {m}")
    except ImportError:
        print(f"  [MISS] {m}")
        bad.append(m)
if bad:
    pip_name = {"sklearn": "scikit-learn"}
    pkgs = [pip_name.get(m, m) for m in bad]
    print("         install: bash ../setup_environments.sh --only viralbin")
    print("         (or, into the existing env: pip install " + " ".join(pkgs) + ")")
    sys.exit(1)
PY
else
    echo "  [MISS] no python found at $PY"
    missing=$((missing + 1))
fi

echo
echo "=== pipeline scripts ($SCRIPT_DIR) ==="
for script in \
    "5_assemble_and_build_graph.py" \
    "6_acceptance_check.py"; do
    if [ -f "$SCRIPT_DIR/$script" ]; then
        printf '  [ok]   %s\n' "$script"
    else
        printf '  [MISS] %s\n' "$script"
        missing=$((missing + 1))
    fi
done

echo
if [ "$missing" -gt 0 ]; then
    echo "preflight FAILED: $missing item(s) missing - install/fix them, then re-run this check."
    exit 1
fi
echo "preflight OK - all sources and tools resolved."
