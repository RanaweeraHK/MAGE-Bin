#!/usr/bin/env bash
       
#   VB_TIER=low VB_SIM=iss bash 4_simulate_reads.sh
set -euo pipefail

_D="$(cd "$(dirname "$0")" && pwd)"
[ -f "$_D/../env.sh" ] && . "$_D/../env.sh"

TIER="${VB_TIER:-high}"
case "$TIER" in
    low)    TIER_SUFFIX="_low";    READS_PER_SAMPLE=40000 ;;
    medium) TIER_SUFFIX="_medium"; READS_PER_SAMPLE=250000 ;;
    high)   TIER_SUFFIX="_high";   READS_PER_SAMPLE=1500000 ;;
    *) echo "error: VB_TIER must be low|medium|high (got '$TIER')" >&2; exit 2 ;;
esac

# ---- simulator config (VB_SIM=iss|badread, default iss) ----
SIM="${VB_SIM:-iss}"
case "$SIM" in
    iss|badread) ;;
    *) echo "error: VB_SIM must be iss|badread (got '$SIM')" >&2; exit 2 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Two dirname levels: this script lives in Dataset_Processing/Simulated_Data/.
REPO_HOME="$(dirname "$(dirname "$SCRIPT_DIR")")"
WORK_DIR="$REPO_HOME/Data/work/multibiome$TIER_SUFFIX"
POOL_DIR="$WORK_DIR/genome_pool"
ABUND_DIR="$WORK_DIR/abundance"

if [ "$SIM" = "iss" ]; then
    READS_DIR="$WORK_DIR/reads"
else
    READS_DIR="$WORK_DIR/reads_$SIM"
fi
CONDA_BIN="${VB_CONDA_BIN:-$HOME/miniforge3/envs/viralbin/bin}"
THREADS="${VB_THREADS:-8}"
N_SAMPLES=15
READS_SEED_BASE=42

# iss: 2x126 bp HiSeq. The read length is also the conversion factor used to
# give badread an equal budget in bases.
ISS_MODEL=hiseq
ISS_READ_LENGTH=126
BASES_PER_SAMPLE=$((READS_PER_SAMPLE * ISS_READ_LENGTH))

# badread: Badread's own current ONT defaults, stated here rather than left
# implicit so the profile is visible and changeable in one place.
BADREAD_MODEL=nanopore2023
BADREAD_LENGTH="15000,13000"     # fragment length mean,stdev
BADREAD_IDENTITY="95,99,2.5"     # beta(mean,max,stdev) read identity

BADREAD_JOBS="${VB_BADREAD_JOBS:-$THREADS}"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
mkdir -p "$READS_DIR"

if [ ! -f "$POOL_DIR/expanded_genomes.fasta" ]; then
    echo "error: $POOL_DIR/expanded_genomes.fasta not found - run 2_expand_strains.py first" >&2
    exit 1
fi

resolve() {  # resolve <binary> -- prefer CONDA_BIN, fall back to PATH
    if [ -x "$CONDA_BIN/$1" ]; then printf '%s\n' "$CONDA_BIN/$1"
    elif command -v "$1" >/dev/null 2>&1; then command -v "$1"
    else return 1
    fi
}

PY="$(resolve python || true)"
[ -n "$PY" ] || PY=python3

# ============================================================================
# iss - Illumina paired-end
# ============================================================================
simulate_iss() {
    local j="$1"
    "$ISS" generate --genomes "$POOL_DIR/expanded_genomes.fasta" \
        --abundance_file "$ABUND_DIR/sample${j}_abundance.txt" \
        --n_reads "$READS_PER_SAMPLE" --model "$ISS_MODEL" --cpus "$THREADS" \
        --seed $((READS_SEED_BASE + j)) --compress \
        --output "$READS_DIR/sample${j}"
    # iss also leaves a copy of the abundance file + vcf/tmp files next to the
    # reads; drop those so READS_DIR only holds the two fastq.gz per sample.
    rm -f "$READS_DIR/sample${j}"_abundance.txt "$READS_DIR/sample${j}".iss.tmp.* 2>/dev/null || true
}

# ============================================================================
# badread - Oxford Nanopore single-end
#
# abundance -> depth.  Phase 3's tables give each strain a share of the sample's
# READS; iss consumes that directly. Badread instead reads a per-sequence
# `depth=` from the FASTA header and draws reads in proportion to depth x
# length. Setting depth = abundance / length therefore reproduces exactly the
# iss read allocation, and hence the same per-genome COVERAGE pattern that the
# existing short-read datasets have - which is what differential-coverage
# binning is being asked to recover. Strains absent from a sample get depth 0,
# which Badread honours by emitting no reads for them.
# ============================================================================
badread_reference() {  # badread_reference <sample-index> <destination fasta>
    VB_ABUNDANCE="$ABUND_DIR/sample${1}_abundance.txt" \
    VB_REFERENCE="$POOL_DIR/expanded_genomes.fasta" \
    VB_DESTINATION="$2" \
    "$PY" - <<'PY'
import os

abundance = {}
with open(os.environ["VB_ABUNDANCE"]) as fh:
    for line in fh:
        if line.strip():
            name, value = line.split()[:2]
            abundance[name] = float(value)

# Pass one: sequence lengths, so abundance (a share of reads) can be turned
# into a depth (a share of coverage).
lengths, name = {}, None
with open(os.environ["VB_REFERENCE"]) as fh:
    for line in fh:
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        elif name is not None:
            lengths[name] += len(line.strip())

missing = sorted(set(abundance) - set(lengths))
if missing:
    raise SystemExit(f"error: {len(missing)} abundance ids absent from the "
                     f"reference, first few: {missing[:5]}")

depth = {n: (abundance.get(n, 0.0) / L if L else 0.0) for n, L in lengths.items()}
# Badread uses only the ratios between depths; scale so the mean non-zero depth
# is 1.0 purely so the headers are readable.
positive = [d for d in depth.values() if d > 0]
scale = (len(positive) / sum(positive)) if positive else 1.0

with open(os.environ["VB_REFERENCE"]) as src, \
        open(os.environ["VB_DESTINATION"], "w") as dst:
    for line in src:
        if line.startswith(">"):
            n = line[1:].split()[0]
            dst.write(f">{n} depth={depth[n] * scale:.10g}\n")
        else:
            dst.write(line)
PY
}

simulate_badread() {
    local j="$1"
    local reference="$READS_DIR/sample${j}.depth.fasta"
    badread_reference "$j" "$reference"
    # Badread streams FASTQ on stdout and a per-read progress bar on stderr;
    # the progress goes to a log so the pipeline output stays readable.
    "$BADREAD" simulate --reference "$reference" \
        --quantity "$BASES_PER_SAMPLE" \
        --length "$BADREAD_LENGTH" --identity "$BADREAD_IDENTITY" \
        --error_model "$BADREAD_MODEL" --qscore_model "$BADREAD_MODEL" \
        --seed $((READS_SEED_BASE + j)) \
        2> "$READS_DIR/sample${j}.badread.log" \
        | gzip -c > "$READS_DIR/sample${j}.fastq.gz.partial"
    mv "$READS_DIR/sample${j}.fastq.gz.partial" "$READS_DIR/sample${j}.fastq.gz"
    rm -f "$reference"
}

# ============================================================================
# dispatch
# ============================================================================
case "$SIM" in
    iss)
        ISS="$(resolve iss || true)"
        if [ -z "$ISS" ]; then
            echo "error: iss not found (checked $CONDA_BIN and PATH)" >&2; exit 1
        fi
        log "simulator: iss $("$ISS" --version 2>&1 | head -1 || true) -> $ISS"
        log "tier=$TIER: $N_SAMPLES samples, $READS_PER_SAMPLE reads/sample, model=$ISS_MODEL"
        ;;
    badread)
        BADREAD="$(resolve badread || true)"
        if [ -z "$BADREAD" ]; then
            echo "error: badread not found (checked $CONDA_BIN and PATH)" >&2
            echo "       install: mamba install -n viralbin -c bioconda badread" >&2
            exit 1
        fi
        log "simulator: $("$BADREAD" --version 2>&1 | head -1 || true) -> $BADREAD"
        log "tier=$TIER: $N_SAMPLES samples, $BASES_PER_SAMPLE bases/sample" \
            "(= $READS_PER_SAMPLE x ${ISS_READ_LENGTH}bp, the iss budget), model=$BADREAD_MODEL"
        ;;
esac
log "reads -> $READS_DIR"


PENDING=()
for j in $(seq 0 $((N_SAMPLES - 1))); do
    if [ -f "$READS_DIR/sample${j}.done" ]; then
        log "[skip] sample${j} already done"
    else
        PENDING+=("$j")
    fi
done

if [ "${#PENDING[@]}" -eq 0 ]; then
    log "ALL READS DONE -> $READS_DIR"
    exit 0
fi

# one_sample <index> -- simulate and mark done. 
one_sample() {
    local j="$1"
    log "=== sample${j} start ==="
    # `set -e` (and pipefail, for badread's simulate | gzip) mean a failed
    # simulator aborts here rather than falling through to the marker - a
    # half-written sample must never look done to the next run or to phase 5.
    "simulate_$SIM" "$j"
    touch "$READS_DIR/sample${j}.done"
    log "[ok] sample${j} done"
}

if [ "$SIM" = "badread" ] && [ "$BADREAD_JOBS" -gt 1 ]; then
    log "running ${#PENDING[@]} samples, $BADREAD_JOBS at a time"
    failed=0
    for j in "${PENDING[@]}"; do
        while [ "$(jobs -rp | wc -l)" -ge "$BADREAD_JOBS" ]; do
            wait -n || failed=1
        done
        one_sample "$j" &
    done
    while [ -n "$(jobs -rp)" ]; do
        wait -n || failed=1
    done
    if [ "$failed" -ne 0 ]; then
        echo "error: at least one sample failed to simulate; see" \
             "$READS_DIR/sample*.badread.log" >&2
        exit 1
    fi
else
    for j in "${PENDING[@]}"; do
        one_sample "$j"
    done
fi

for j in $(seq 0 $((N_SAMPLES - 1))); do
    if [ ! -f "$READS_DIR/sample${j}.done" ]; then
        echo "error: sample${j} did not complete" >&2
        exit 1
    fi
done

log "ALL READS DONE -> $READS_DIR"
