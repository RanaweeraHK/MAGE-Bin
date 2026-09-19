#!/usr/bin/env bash
# Build the conda environments both dataset pipelines need, then write env.sh.
#
#   bash setup_environments.sh              # create every missing env, write env.sh
#   bash setup_environments.sh --list       # what would be done, then exit
#   bash setup_environments.sh --dry-run    # print the commands, run nothing
#   bash setup_environments.sh --only viralbin asm_env
#   bash setup_environments.sh --force      # recreate envs that already exist
#   bash setup_environments.sh --databases  # also download CheckV + geNomad DBs

set -euo pipefail

D="$(cd "$(dirname "$0")" && pwd)"

# ---- the environments, in the order they are created --------------------
# name|used by|packages (conda)|packages (pip, after the conda solve)
ENVS=(

"viralbin|both pipelines: mmseqs, checkv, flye, badread, and the python holding torch/PyG|python=3.12 mmseqs2 checkv flye pyrodigal-gv pysam numpy>=2.4 pandas networkx scikit-learn|torch torch_geometric badread InSilicoSeq"
"asm_env|both pipelines: metaspades.py, minimap2, samtools, seqkit|spades>=3.15 minimap2>=2.24 samtools>=1.17 seqkit>=2.5 pigz|"
"sra_env|Real_Data phase 1: prefetch, fasterq-dump|sra-tools|"
"genomad_env|Real_Data: genomad (only for viral_enriched=no datasets)|genomad|"
"megahit_env|Real_Data: megahit (only for assembler=megahit datasets)|megahit|"
)

CHANNELS=(-c conda-forge -c bioconda)

LIST_ONLY=0
DRY_RUN=0
FORCE=0
DO_DBS=0
ONLY=()
while [ $# -gt 0 ]; do
    case "$1" in
        --list) LIST_ONLY=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --force) FORCE=1; shift ;;
        --databases) DO_DBS=1; shift ;;
        --only)
            shift
            while [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; do ONLY+=("$1"); shift; done ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "error: unknown flag '$1' (try --help)" >&2; exit 2 ;;
    esac
done

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

# ---- conda front-end: mamba is much faster, but plain conda works -------
CONDA_EXE=""
for c in mamba micromamba conda; do
    if command -v "$c" >/dev/null 2>&1; then CONDA_EXE="$c"; break; fi
done
if [ -z "$CONDA_EXE" ]; then
    echo "error: no mamba/micromamba/conda on PATH." >&2
    echo "  install miniforge first: https://github.com/conda-forge/miniforge" >&2
    exit 1
fi

# Where envs live, so the script can tell "already exists" from "to build" and
# can write the right paths into env.sh. `conda info --base` is the source of
# truth; CONDA_PREFIX_BASE falls back to a stock miniforge layout.
ENV_ROOT="$(conda info --base 2>/dev/null || true)"
[ -n "$ENV_ROOT" ] || ENV_ROOT="${MAMBA_ROOT_PREFIX:-$HOME/miniforge3}"
ENV_ROOT="$ENV_ROOT/envs"

wanted() {  # wanted <env> -- is this env in --only (or was --only not given)?
    [ ${#ONLY[@]} -eq 0 ] && return 0
    local e; for e in "${ONLY[@]}"; do [ "$e" = "$1" ] && return 0; done
    return 1
}

run() {  # run <cmd...> -- honour --dry-run
    if [ "$DRY_RUN" = 1 ]; then printf '  would run: %s\n' "$*"; return 0; fi
    "$@"
}

echo "conda front-end : $CONDA_EXE"
echo "environments in : $ENV_ROOT"
echo

created=0
skipped=0
for spec in "${ENVS[@]}"; do
    IFS='|' read -r name purpose conda_pkgs pip_pkgs <<< "$spec"
    wanted "$name" || continue

    if [ -x "$ENV_ROOT/$name/bin/python" ] || [ -d "$ENV_ROOT/$name/bin" ]; then
        if [ "$FORCE" = 0 ]; then
            printf '  [have] %-12s %s\n' "$name" "$purpose"
            skipped=$((skipped + 1))
            continue
        fi
        printf '  [redo] %-12s (--force: removing existing env)\n' "$name"
        [ "$LIST_ONLY" = 1 ] || run "$CONDA_EXE" env remove -n "$name" -y
    else
        printf '  [make] %-12s %s\n' "$name" "$purpose"
    fi
    created=$((created + 1))
    [ "$LIST_ONLY" = 1 ] && continue

    log "creating $name ..."
    # shellcheck disable=SC2086  # package list is intentionally word-split
    run "$CONDA_EXE" create -y -n "$name" "${CHANNELS[@]}" $conda_pkgs
    if [ -n "$pip_pkgs" ]; then
        # torch/PyG are installed with pip, not conda: the conda-forge builds
        # drag in a CUDA/BLAS stack this pipeline never uses, and pip's default
        # wheel is the CPU build that every phase here actually runs on.
        log "pip installing into $name: $pip_pkgs"
        # shellcheck disable=SC2086
        run "$ENV_ROOT/$name/bin/python" -m pip install --quiet $pip_pkgs
    fi
done

if [ "$LIST_ONLY" = 1 ]; then
    echo
    echo "$created to create, $skipped already present. Nothing was changed."
    exit 0
fi

# ---- env.sh -------------------------------------------------------------
echo
if [ -f "$D/env.sh" ]; then
    log "env.sh already exists - left untouched (delete it to regenerate)"
else
    if [ "$DRY_RUN" = 1 ]; then
        echo "  would write $D/env.sh from env.sh.example"
    else
        # Keep the $HOME form when the envs really are under $HOME - the file
        # then still reads as the template does. Only a root somewhere else
        # (a shared /opt install, a conda outside the home directory) gets
        # written out as an absolute path.
        case "$ENV_ROOT" in
            "$HOME"/*) replacement="\$HOME/${ENV_ROOT#"$HOME"/}" ;;
            *) replacement="$ENV_ROOT" ;;
        esac
        sed -e "s#\$HOME/miniforge3/envs#$replacement#g" \
            "$D/env.sh.example" > "$D/env.sh"
        log "wrote $D/env.sh (environments in $ENV_ROOT)"
    fi
fi

# ---- databases (opt-in: ~1.4 GB geNomad, several GB CheckV) -------------
if [ "$DO_DBS" = 1 ]; then
    # shellcheck source=/dev/null
    [ -f "$D/env.sh" ] && . "$D/env.sh"
    CHECKV_DB="${CHECKVDB:-$HOME/checkv-db/checkv-db-v1.5}"
    GDB="${GENOMAD_DB:-$HOME/db/genomad/genomad_db}"
    if [ -d "$CHECKV_DB/genome_db" ]; then
        log "checkv database already at $CHECKV_DB"
    else
        log "downloading CheckV database to $(dirname "$CHECKV_DB")"
        run "$ENV_ROOT/viralbin/bin/checkv" download_database "$(dirname "$CHECKV_DB")"
    fi
    if [ -d "$GDB" ]; then
        log "genomad database already at $GDB"
    else
        log "downloading geNomad database to $(dirname "$GDB")"
        run "$ENV_ROOT/genomad_env/bin/genomad" download-database "$(dirname "$GDB")"
    fi
fi

echo
echo "next:"
echo "  1. review $D/env.sh (threads, memory, database paths)"
echo "  2. bash Simulated_Data/0_check_environment.sh   # simulated pipeline preflight"
echo "  3. bash Real_Data/0_check_environment.sh        # real-data pipeline preflight"
echo "  Each prints where every tool resolved to; fix env.sh until both are clean."
