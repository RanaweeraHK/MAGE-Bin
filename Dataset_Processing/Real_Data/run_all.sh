#!/usr/bin/env bash
# Build the REAL-dataset benchmarks, phases 0-3.
#
#   bash run_all.sh                          # build every dataset in datasets.tsv
#   bash run_all.sh --dataset aloha_illumina # just this one
#   bash run_all.sh --list                   # print what would be built, then exit
#   bash run_all.sh --force                  # rebuild every phase from scratch
#   bash run_all.sh --from 2                 # resume starting at a given phase
#
# WHICH DATASETS GET BUILT is `datasets.tsv`, not a list in this file - a real
# dataset is a set of SRA accessions plus its library type, which is data rather
# than a flag. Edit that file, then run this.
#
# Phase 4 (CheckV evaluation) is deliberately NOT part of this chain: it scores
# BINS, which only exist after a binner has run. See the README.
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"

# ---- config ----
# Two dirname levels: this script lives in Dataset_Processing/Real_Data/.
REPO_HOME="$(dirname "$(dirname "$D")")"
REGISTRY="$D/datasets.tsv"

[ -f "$D/../env.sh" ] && . "$D/../env.sh"
CONDA_BIN="${VB_CONDA_BIN:-$HOME/miniforge3/envs/viralbin/bin}"
PY="$CONDA_BIN/python"
[ -x "$PY" ] || PY="python3"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

FORCE=0
FROM=1
LIST_ONLY=0
SELECTED=()
while [ $# -gt 0 ]; do
    case "$1" in
        --force) FORCE=1; shift ;;
        --list) LIST_ONLY=1; shift ;;
        --from) FROM="$2"; shift 2 ;;
        --dataset) SELECTED+=("$2"); shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[ -f "$REGISTRY" ] || { echo "error: $REGISTRY not found" >&2; exit 2; }

# Every dataset name in the registry, in file order.
mapfile -t ALL < <(awk -F'\t' '!/^#/ && NF >= 11 && $1 != "dataset" {print $1}' "$REGISTRY")
if [ "${#SELECTED[@]}" -gt 0 ]; then
    BUILD=("${SELECTED[@]}")
else
    BUILD=("${ALL[@]}")
fi
if [ "${#BUILD[@]}" -eq 0 ]; then
    echo "error: no datasets selected (is $REGISTRY empty?)" >&2
    exit 2
fi

# field <dataset> <column-number>
field() {
    awk -F'\t' -v d="$1" -v c="$2" '!/^#/ && $1 == d {print $c; exit}' "$REGISTRY"
}

# run_name <dataset> -- mirrored in 2_assemble_and_build_graph.py and
# 3_acceptance_check.py; keep the three in step.
run_name() {
    local read_type viral asm tag declared
    read_type="$(field "$1" 4)"
    viral="$(field "$1" 5)"
    declared="$(field "$1" 9)"
    if [ -n "$declared" ] && [ "$declared" != "-" ]; then
        asm="$declared"
    else
        case "$read_type" in
            short) asm=metaspades ;;
            long)  asm=metaflye ;;
            *) echo "error: dataset '$1' has read_type '$read_type' (use short|long)" >&2; exit 2 ;;
        esac
    fi
    case "$viral" in
        yes) tag=noviral ;;
        no)  tag=genomad ;;
        *) echo "error: dataset '$1' has viral_enriched '$viral' (use yes|no)" >&2; exit 2 ;;
    esac
    printf '%s__coasm__%s__%s' "$1" "$asm" "$tag"
}

echo "datasets to build (${#BUILD[@]}):"
for name in "${BUILD[@]}"; do
    if [ -z "$(field "$name" 1)" ]; then
        echo "error: dataset '$name' is not in $REGISTRY" >&2
        echo "known: ${ALL[*]}" >&2
        exit 2
    fi
    printf '  %-14s %-9s %-6s %d run(s) -> Data/Processed_data/%s\n' \
        "$name" "$(field "$name" 6)" "$(field "$name" 4)" \
        "$(field "$name" 3 | tr ',' '\n' | grep -c .)" "$(run_name "$name")"
done
if [ "$LIST_ONLY" -eq 1 ]; then
    exit 0
fi

# skip_if <phase-number> <marker-file...>
skip_if() {
    local phase="$1"; shift
    [ "$phase" -lt "$FROM" ] && return 0        # earlier than --from: always skip
    [ "$FORCE" -eq 1 ] && return 1
    for f in "$@"; do
        [ -e "$f" ] || return 1
    done
    return 0
}

run_dataset() {
    local name="$1"
    export VB_DATASET="$name"
    local run; run="$(run_name "$name")"
    local reads_dir="$REPO_HOME/Data/work/real_$name/reads"
    local out_dir="$REPO_HOME/Data/Processed_data/$run"

    echo
    echo "##################################################################"
    echo "##########  $name  ->  $out_dir"
    echo "##################################################################"

    local source; source="$(field "$name" 6)"
    local mwork="$REPO_HOME/Data/work/real_$name/standardize_work/$run"

    echo "########## 1: get inputs ($source) ##########"
    if [ "$source" = "assembly" ]; then
        # A pre-built assembly: phase 1 imports it (and any adopted artifacts)
        # into the work directory. Cheap and idempotent, so it always runs -
        # it is what decides which later stages can be skipped.
        bash "$D/1_get_inputs.sh"
    else
        # One marker per accession: a partially fetched dataset must not be skipped.
        local ready=1
        for acc in $(field "$name" 3 | tr ',' ' '); do
            compgen -G "$reads_dir/${acc}*.fastq.gz" >/dev/null || ready=0
        done
        if [ "$ready" -eq 1 ] && skip_if 1 "$reads_dir"; then
            log "[skip] reads already fetched -> $reads_dir"
        else
            bash "$D/1_get_inputs.sh"
        fi
    fi

    echo "########## 2: assemble / coverage / viral ID / graph ##########"
    # For source=assembly, phase 2 adopts whatever phase 1 put in $mwork and
    # only computes what is missing; it prints [adopted] for each stage it skips.
    if skip_if 2 "$out_dir/viral_graph.pt"; then
        log "[skip] graph already built -> $out_dir/viral_graph.pt"
    else
        local extra=()
        [ "$FORCE" -eq 1 ] && extra+=(--force)
        "$PY" "$D/2_assemble_and_build_graph.py" "${extra[@]}"
    fi

    echo "########## 3: acceptance check ##########"
    "$PY" "$D/3_acceptance_check.py"
}

START=$(date +%s)
echo
echo "########## 0: environment check ##########"
# Only demand the assemblers the selected datasets actually need.
REQ_ASMS=""
for name in "${BUILD[@]}"; do
    # A pre-built assembly needs no assembler installed on this machine.
    [ "$(field "$name" 6)" = "assembly" ] && continue
    declared="$(field "$name" 9)"
    if [ -n "$declared" ] && [ "$declared" != "-" ]; then
        REQ_ASMS="$REQ_ASMS $declared"
    else
        case "$(field "$name" 4)" in
            short) REQ_ASMS="$REQ_ASMS metaspades" ;;
            long)  REQ_ASMS="$REQ_ASMS metaflye" ;;
        esac
    fi
done
VB_REQUIRED_ASMS="$(printf '%s\n' $REQ_ASMS | sort -u | tr '\n' ' ')" \
    bash "$D/0_check_environment.sh"

for name in "${BUILD[@]}"; do
    run_dataset "$name"
done

ELAPSED=$(( $(date +%s) - START ))
printf '\n########## DONE in %dh%02dm%02ds ##########\n' \
    $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
for name in "${BUILD[@]}"; do
    echo "  $name -> $REPO_HOME/Data/Processed_data/$(run_name "$name")"
done
echo
echo "Next: run a binner, then score its bins without a reference:"
echo "  python $D/4_checkv_evaluate.py --dataset <dataset> --bins <assignments.tsv>"
