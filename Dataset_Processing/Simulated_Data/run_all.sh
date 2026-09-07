#!/usr/bin/env bash
# Build the multi-biome benchmark datasets, phases 0-6.
#
#   bash run_all.sh                              # build everything in DATASETS
#   bash run_all.sh --dataset high:badread:metaflye     # just this one
#   bash run_all.sh --tier low                   # every DATASETS entry at this tier
#   bash run_all.sh --force                      # rebuild every phase from scratch
#   bash run_all.sh --from 4                     # resume starting at a given phase
#   bash run_all.sh --list                       # print what would be built, then exit
#
# ============================================================================
# THE DATASETS TO BUILD - edit this list, then run `bash run_all.sh`
# ============================================================================
# One entry per dataset, "tier:simulator:assembler".
#
#   tier       low | medium | high      community complexity (see table below)
#   simulator  iss | badread            Illumina paired-end | ONT long reads
#   assembler  metaspades | metaflye    must match the simulator's read type
#
# Valid pairs are iss+metaspades (short) and badread+metaflye (long); anything
# else is rejected up front rather than failing an hour into an assembly.
DATASETS=(
    "low:iss:metaspades"
    "medium:iss:metaspades"
    "high:iss:metaspades"
    "high:badread:metaflye"
)
#
# Tiers (the complexity knob is genome count x strain complexity):
#
#   tier    genomes  strains  strain complexity  reads/sample (iss)
#   low         50       54          5%              40,000
#   medium     200      296         20%             250,000
#   high       998    1,995         40%          1,500,000
#
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"

# ---- config (see README's parameter table) ----
# Two dirname levels: this script lives in Dataset_Processing/Simulated_Data/.
REPO_HOME="$(dirname "$(dirname "$D")")"

[ -f "$D/../env.sh" ] && . "$D/../env.sh"
CONDA_BIN="${VB_CONDA_BIN:-$HOME/miniforge3/envs/viralbin/bin}"
PY="$CONDA_BIN/python"
[ -x "$PY" ] || PY="python3"
N_SAMPLES=15

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

FORCE=0
FROM=1
LIST_ONLY=0
TIER_FILTER=""
SELECTED=()
while [ $# -gt 0 ]; do
    case "$1" in
        --force) FORCE=1; shift ;;
        --list) LIST_ONLY=1; shift ;;
        --from) FROM="$2"; shift 2 ;;
        --dataset) SELECTED+=("$2"); shift 2 ;;
        --tier)
            case "$2" in
                low|medium|high) TIER_FILTER="$2" ;;
                all) TIER_FILTER="" ;;
                *) echo "unknown tier: $2 (use low|medium|high|all)" >&2; exit 2 ;;
            esac
            shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# --dataset overrides the built-in list; --tier filters whichever list is used.
if [ "${#SELECTED[@]}" -gt 0 ]; then
    BUILD=("${SELECTED[@]}")
else
    BUILD=("${DATASETS[@]}")
fi
if [ -n "$TIER_FILTER" ]; then
    filtered=()
    for spec in "${BUILD[@]}"; do
        [ "${spec%%:*}" = "$TIER_FILTER" ] && filtered+=("$spec")
    done
    BUILD=("${filtered[@]}")
fi
if [ "${#BUILD[@]}" -eq 0 ]; then
    echo "error: no datasets selected (DATASETS list empty, or --tier filtered everything out)" >&2
    exit 2
fi

# ---- naming + validation, applied to every entry BEFORE any work starts ----
# tier_suffix <tier>
tier_suffix() {
    case "$1" in
        low) printf '_low' ;;
        medium) printf '_medium' ;;
        high) printf '_high' ;;
    esac
}

# run_name <tier> <sim> <asm>  ->  the Data/Processed_data directory name.
# Mirrored in 5_assemble_and_build_graph.py; keep the two in step.
run_name() {
    local sim_tag=""
    [ "$2" != "iss" ] && sim_tag="__$2"
    printf 'multibiome%s%s__coasm__%s__noviral' "$(tier_suffix "$1")" "$sim_tag" "$3"
}

validate() {  # validate <tier> <sim> <asm>
    case "$1" in low|medium|high) ;; *)
        echo "error: bad tier '$1' in '$1:$2:$3' (use low|medium|high)" >&2; exit 2 ;;
    esac
    case "$2" in iss|badread) ;; *)
        echo "error: bad simulator '$2' in '$1:$2:$3' (use iss|badread)" >&2; exit 2 ;;
    esac
    case "$3" in metaspades|metaflye) ;; *)
        echo "error: bad assembler '$3' in '$1:$2:$3' (use metaspades|metaflye)" >&2; exit 2 ;;
    esac
    local sim_reads asm_reads
    case "$2" in iss) sim_reads=short ;; badread) sim_reads=long ;; esac
    case "$3" in metaspades) asm_reads=short ;; metaflye) asm_reads=long ;; esac
    if [ "$sim_reads" != "$asm_reads" ]; then
        echo "error: '$1:$2:$3' pairs a $sim_reads-read simulator with a" \
             "$asm_reads-read assembler." >&2
        echo "       Valid pairs: iss+metaspades (short), badread+metaflye (long)." >&2
        exit 2
    fi
}

echo "datasets to build (${#BUILD[@]}):"
for spec in "${BUILD[@]}"; do
    IFS=: read -r tier sim asm <<< "$spec"
    validate "$tier" "$sim" "$asm"
    printf '  %-26s -> Data/Processed_data/%s\n' "$spec" "$(run_name "$tier" "$sim" "$asm")"
done
if [ "$LIST_ONLY" -eq 1 ]; then
    exit 0
fi

# skip_if <phase-number> <marker-file...> -- true (skip) when not forced,
# phase >= FROM has already been reached, and every marker file exists.
skip_if() {
    local phase="$1"; shift
    [ "$phase" -lt "$FROM" ] && return 0        # earlier than --from: always skip
    [ "$FORCE" -eq 1 ] && return 1
    for f in "$@"; do
        [ -e "$f" ] || return 1
    done
    return 0
}

# run_dataset <tier> <sim> <asm> -- phases 1-7 for one dataset. Every numbered
# script reads VB_TIER/VB_SIM/VB_ASM for its own paths and parameters, so they
# are set once here.
run_dataset() {
    local tier="$1" sim="$2" asm="$3"
    export VB_TIER="$tier" VB_SIM="$sim" VB_ASM="$asm"
    local suffix; suffix="$(tier_suffix "$tier")"
    local name; name="$(run_name "$tier" "$sim" "$asm")"

    local work_dir="$REPO_HOME/Data/work/multibiome$suffix"
    local pool_dir="$work_dir/genome_pool"
    local abund_dir="$work_dir/abundance"
    local reads_dir="$work_dir/reads"
    [ "$sim" != "iss" ] && reads_dir="$work_dir/reads_$sim"
    local out_dir="$REPO_HOME/Data/Processed_data/$name"

    echo
    echo "##################################################################"
    echo "##########  $tier / $sim / $asm  ->  $out_dir"
    echo "##################################################################"

    # Phases 1-3 depend only on the tier and are shared by every simulator and
    # assembler at that tier, so a second dataset at the same tier reuses them.
    echo "########## 1: genome pool ##########"
    if skip_if 1 "$pool_dir/community_pool.fasta" "$pool_dir/genome_metadata.tsv"; then
        log "[skip] genome pool already built"
    else
        "$PY" "$D/1_build_genome_pool.py"
    fi

    echo "########## 2: strain expansion ##########"
    if skip_if 2 "$pool_dir/expanded_genomes.fasta" "$pool_dir/strain_map.tsv"; then
        log "[skip] strain expansion already built"
    else
        "$PY" "$D/2_expand_strains.py"
    fi

    echo "########## 3: differential abundance ##########"
    if skip_if 3 "$abund_dir/sample0_abundance.txt"; then
        log "[skip] abundance tables already built"
    else
        "$PY" "$D/3_generate_abundance.py"
    fi

    echo "########## 4: read simulation ($sim) ##########"
    local last_sample=$((N_SAMPLES - 1))
    if skip_if 4 "$reads_dir/sample${last_sample}.done"; then
        log "[skip] reads already simulated -> $reads_dir"
    else
        bash "$D/4_simulate_reads.sh"
    fi

    echo "########## 5: assemble ($asm) / coverage / labels / graph ##########"
    if skip_if 5 "$out_dir/viral_graph.pt"; then
        log "[skip] graph already built -> $out_dir/viral_graph.pt"
    else
        local extra=()
        [ "$FORCE" -eq 1 ] && extra+=(--force)
        "$PY" "$D/5_assemble_and_build_graph.py" "${extra[@]}"
    fi

    echo "########## 6: acceptance check ##########"
    "$PY" "$D/6_acceptance_check.py"
}

START=$(date +%s)
echo
echo "########## 0: environment check ##########"
# The check only demands the tools the selected datasets actually need.
VB_REQUIRED_SIMS="$(printf '%s\n' "${BUILD[@]}" | cut -d: -f2 | sort -u | tr '\n' ' ')" \
VB_REQUIRED_ASMS="$(printf '%s\n' "${BUILD[@]}" | cut -d: -f3 | sort -u | tr '\n' ' ')" \
    bash "$D/0_check_environment.sh"

for spec in "${BUILD[@]}"; do
    IFS=: read -r tier sim asm <<< "$spec"
    run_dataset "$tier" "$sim" "$asm"
done

ELAPSED=$(( $(date +%s) - START ))
printf '\n########## DONE in %dh%02dm%02ds ##########\n' \
    $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
for spec in "${BUILD[@]}"; do
    IFS=: read -r tier sim asm <<< "$spec"
    echo "  $spec -> $REPO_HOME/Data/Processed_data/$(run_name "$tier" "$sim" "$asm")"
done
