#!/usr/bin/env bash
# Build the MOCK-COMMUNITY benchmarks, phases 0-4.
#
#   bash run_all.sh                             # build every dataset in datasets.tsv
#   bash run_all.sh --dataset phage_mock_illumina
#   bash run_all.sh --list                      # print what would be built, then exit
#   bash run_all.sh --force                     # rebuild every phase from scratch
#   bash run_all.sh --from 2                    # resume starting at a given phase
#   bash run_all.sh --discard-downloads         # delete provider files once converted
#
# Which datasets get built is `datasets.tsv`, as in the sibling pipelines.
#
# Phase 3 (labelling) IS part of this chain: it needs only the assembly. CheckV
# is not, because it scores bins, which exist only after a binner has run.
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"

# Two dirname levels: this script lives in Dataset_Processing/Mock_Data/.
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
DISCARD=0
SELECTED=()
while [ $# -gt 0 ]; do
    case "$1" in
        --force) FORCE=1; shift ;;
        --list) LIST_ONLY=1; shift ;;
        --from) FROM="$2"; shift 2 ;;
        --dataset) SELECTED+=("$2"); shift 2 ;;
        --discard-downloads) DISCARD=1; shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[ -f "$REGISTRY" ] || { echo "error: $REGISTRY not found" >&2; exit 2; }

mapfile -t ALL < <(awk -F'\t' '!/^#/ && NF >= 15 && $1 != "dataset" {print $1}' "$REGISTRY")
if [ "${#SELECTED[@]}" -gt 0 ]; then
    BUILD=("${SELECTED[@]}")
else
    BUILD=("${ALL[@]}")
fi
[ "${#BUILD[@]}" -gt 0 ] || { echo "error: no datasets selected (is $REGISTRY empty?)" >&2; exit 2; }

# field <dataset> <column-number>
field() {
    awk -F'\t' -v d="$1" -v c="$2" '!/^#/ && $1 == d {print $c; exit}' "$REGISTRY"
}

# run_name <dataset> -- also in 3_label_from_references.py and phase 2.
run_name() {
    local read_type viral asm declared
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
        yes) echo "$1__coasm__${asm}__noviral" ;;
        no)  echo "$1__coasm__${asm}__genomad" ;;
        *) echo "error: dataset '$1' has viral_enriched '$viral' (use yes|no)" >&2; exit 2 ;;
    esac
}

if [ "$LIST_ONLY" -eq 1 ]; then
    printf '%-24s %-18s %-6s %-10s %s\n' DATASET STUDY READS PROVIDER RUN
    for d in "${BUILD[@]}"; do
        printf '%-24s %-18s %-6s %-10s %s\n' \
            "$d" "$(field "$d" 2)" "$(field "$d" 4)" "$(field "$d" 12)" "$(run_name "$d")"
    done
    exit 0
fi

for DATASET in "${BUILD[@]}"; do
    RUN="$(run_name "$DATASET")"
    OUT="$REPO_HOME/Data/Processed_data/$RUN"
    export VB_DATASET="$DATASET"

    log "=== $DATASET  ->  $RUN"

    if [ "$FROM" -le 0 ]; then
        bash "$D/0_check_environment.sh"
    fi

    # Phase 1 always runs: it downloads only what is not already there.
    if [ "$FROM" -le 1 ]; then
        log "[1/4] inputs"
        if [ "$DISCARD" -eq 1 ]; then
            "$PY" "$D/1_get_inputs.py" --discard-downloads
        else
            "$PY" "$D/1_get_inputs.py"
        fi
    fi

    # Phase 2 skips a dataset that is already built, unless --force.
    if [ "$FROM" -le 2 ]; then
        if [ "$FORCE" -eq 0 ] && [ -f "$OUT/viral_graph.pt" ] && [ -f "$OUT/manifest.json" ]; then
            log "[2/4] already built: $OUT (use --force to rebuild)"
        else
            log "[2/4] assemble + graph"
            if [ "$FORCE" -eq 1 ]; then
                "$PY" "$D/2_assemble_and_build_graph.py" --force
            else
                "$PY" "$D/2_assemble_and_build_graph.py"
            fi
        fi
    fi

    # Phase 3 rewrites only the label columns, so it is cheap to re-run.
    if [ "$FROM" -le 3 ]; then
        log "[3/4] label against the reference set"
        if [ "$FORCE" -eq 1 ]; then
            "$PY" "$D/3_label_from_references.py" --force
        else
            "$PY" "$D/3_label_from_references.py"
        fi
    fi

    if [ "$FROM" -le 4 ]; then
        log "[4/4] acceptance check"
        "$PY" "$D/4_acceptance_check.py" || log "acceptance check reported problems (continuing)"
    fi

    log "=== $DATASET done"
done

log "all datasets complete"
