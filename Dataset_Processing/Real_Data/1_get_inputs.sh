#!/usr/bin/env bash
# Phase 1 - put one dataset's inputs where phase 2 expects them.
#
#   VB_DATASET=aloha_illumina bash 1_get_inputs.sh    # source=sra      -> reads
#   VB_DATASET=lake_water     bash 1_get_inputs.sh    # source=assembly -> import
#
# Two paths, chosen by the registry's `source` column:
#
#   sra       prefetch each run into Data/Real_dataset/<study>/<accession>/ (the
#             layout already on disk), fasterq-dump it, gzip the result into
#             Data/work/real_<dataset>/reads/.  Phase 2 then assembles.
#
#             The FASTQ names are what phase 2's read detector pairs on:
#             `--split-files` writes <ACC>_1.fastq / <ACC>_2.fastq for a paired
#             run and <ACC>.fastq for a single-end one, which the detector's ENA
#             `_1/_2` pattern already handles - one sample per accession, one
#             coverage column per sample, ordered by accession so the column
#             order is stable across re-runs.
#
#   assembly  the assembly already exists.  Copy it into the run's work
#             directory under the names phase 2 reads, and - if `import_from`
#             names an earlier work directory - ADOPT its finished artifacts
#             (coverage matrix, filtered contigs, geNomad calls) rather than
#             recomputing them.  Nothing is assembled and no reads are needed.
#
# Every step is resumable: an input already in place is skipped, so an
# interrupted run costs only the file it died on.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Two dirname levels: this script lives in Dataset_Processing/Real_Data/.
REPO_HOME="$(dirname "$(dirname "$SCRIPT_DIR")")"
REGISTRY="$SCRIPT_DIR/datasets.tsv"

[ -f "$SCRIPT_DIR/../env.sh" ] && . "$SCRIPT_DIR/../env.sh"
SRA_CONDA_BIN="${VB_SRA_CONDA_BIN:-$HOME/miniforge3/envs/sra_env/bin}"
THREADS="${VB_THREADS:-8}"

PREFETCH_RETRIES=4
PREFETCH_BACKOFF=60
DATASET="${VB_DATASET:?set VB_DATASET to a dataset name from datasets.tsv}"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

# ---- read this dataset's row out of the registry ----
row="$(awk -F'\t' -v d="$DATASET" '!/^#/ && $1 == d {print; exit}' "$REGISTRY")"
if [ -z "$row" ]; then
    echo "error: dataset '$DATASET' is not in $REGISTRY" >&2
    echo "known datasets:" >&2
    awk -F'\t' '!/^#/ && $1 != "dataset" && NF >= 10 {print "  " $1}' "$REGISTRY" >&2
    exit 2
fi
STUDY="$(printf '%s' "$row" | cut -f2)"
ACCESSIONS="$(printf '%s' "$row" | cut -f3 | tr ',' ' ')"
READ_TYPE="$(printf '%s' "$row" | cut -f4)"
VIRAL_ENRICHED="$(printf '%s' "$row" | cut -f5)"
SOURCE="$(printf '%s' "$row" | cut -f6)"
ASSEMBLY_DIR="$(printf '%s' "$row" | cut -f7)"
IMPORT_FROM="$(printf '%s' "$row" | cut -f8)"
ASSEMBLER_COL="$(printf '%s' "$row" | cut -f9)"

# Each assembler names its three outputs differently. metaSPAdes writes
# contigs.fasta / contigs.paths / assembly_graph_with_scaffolds.gfa; metaFlye
# writes assembly.fasta / assembly_info.txt / assembly_graph.gfa. Phase 2 reads
# one set of names, so phase 1 maps whichever the assembler produced onto them.
case "$READ_TYPE" in
    short|long) ;;
    *) echo "error: read_type must be short|long (got '$READ_TYPE')" >&2; exit 2 ;;
esac
# read_type picks the assembler only when this pipeline does the assembling.
# For a pre-built assembly the registry may name it, because the two are
# independent: the assembly's format is whatever produced it, while read_type
# describes the reads being mapped for coverage.
if [ "$ASSEMBLER_COL" != "-" ] && [ -n "$ASSEMBLER_COL" ]; then
    ASSEMBLER="$ASSEMBLER_COL"
else
    case "$READ_TYPE" in
        short) ASSEMBLER=metaspades ;;
        long)  ASSEMBLER=metaflye ;;
    esac
fi

case "$ASSEMBLER" in
    metaspades) PATH_TABLE=contigs.paths
                CONTIG_CANDIDATES="contigs.fasta scaffolds.fasta" ;;
    metaflye)   PATH_TABLE=assembly_info.txt
                CONTIG_CANDIDATES="assembly.fasta contigs.fasta" ;;
    megahit)    PATH_TABLE=""
                CONTIG_CANDIDATES="final.contigs.fa final.contigs.fasta contigs.fasta" ;;
    *) echo "error: assembler must be metaspades|metaflye|megahit (got '$ASSEMBLER')" >&2; exit 2 ;;
esac
case "$VIRAL_ENRICHED" in
    yes) VIRAL_TAG=noviral ;;
    no)  VIRAL_TAG=genomad ;;
    *) echo "error: viral_enriched must be yes|no (got '$VIRAL_ENRICHED')" >&2; exit 2 ;;
esac
RUN_NAME="${DATASET}__coasm__${ASSEMBLER}__${VIRAL_TAG}"

WORK_DIR="$REPO_HOME/Data/work/real_$DATASET"
READS_DIR="$WORK_DIR/reads"
MWORK="$WORK_DIR/standardize_work/$RUN_NAME"

echo "======================================================================"
echo "phase 1  ::  $DATASET  ($SOURCE, $READ_TYPE reads, $ASSEMBLER)"
echo "======================================================================"

resolve() {
    if [ -x "$SRA_CONDA_BIN/$1" ]; then printf '%s' "$SRA_CONDA_BIN/$1"
    elif command -v "$1" >/dev/null 2>&1; then command -v "$1"
    else echo "error: $1 not found (install: mamba create -n sra_env -c bioconda sra-tools)" >&2; exit 3
    fi
}

# fetch_reads -- download every accession and lay it out as per-sample FASTQ.
# Used by source=sra (before assembly) and by source=assembly when the assembly
# was given but its reads still have to be mapped for coverage.
fetch_reads() {
    local PREFETCH FASTERQ SRA_DIR acc sra fq n_files n_acc
    PREFETCH="$(resolve prefetch)"
    FASTERQ="$(resolve fasterq-dump)"
    SRA_DIR="$REPO_HOME/Data/Real_dataset/$STUDY"
    mkdir -p "$SRA_DIR" "$READS_DIR"
    echo "  accessions : $ACCESSIONS"
    echo "  sra cache  : $SRA_DIR"
    echo "  reads      : $READS_DIR"
    echo

    for acc in $ACCESSIONS; do
        if compgen -G "$READS_DIR/${acc}*.fastq.gz" >/dev/null; then
            log "[skip] $acc: FASTQ already present"
            continue
        fi
        sra="$SRA_DIR/$acc/$acc.sra"
        if [ -f "$sra" ]; then
            log "[skip] $acc: .sra already downloaded ($(du -h "$sra" | cut -f1))"
        else
            # A killed prefetch leaves a .sra.lock behind
            attempt=1
            until [ "$attempt" -gt "$PREFETCH_RETRIES" ]; do
                if [ -f "$SRA_DIR/$acc/$acc.sra.lock" ]; then
                    log "$acc: clearing a stale lock from an interrupted download"
                    rm -f "$SRA_DIR/$acc/$acc.sra.lock" \
                          "$SRA_DIR/$acc/$acc.sra.tmp" \
                          "$SRA_DIR/$acc/$acc.sra.prf"
                fi
                log "$acc: prefetch (attempt $attempt/$PREFETCH_RETRIES)"
                if "$PREFETCH" --max-size u --output-directory "$SRA_DIR" "$acc"; then
                    break
                fi
                if [ "$attempt" -eq "$PREFETCH_RETRIES" ]; then
                    echo "error: prefetch failed $PREFETCH_RETRIES times for $acc" >&2
                    echo "       re-run this script to resume; finished runs are skipped" >&2
                    exit 4
                fi
                log "$acc: prefetch failed, retrying in ${PREFETCH_BACKOFF}s"
                sleep "$PREFETCH_BACKOFF"
                attempt=$((attempt + 1))
            done
        fi
        [ -f "$sra" ] || { echo "error: prefetch produced no $sra" >&2; exit 4; }

        log "$acc: fasterq-dump"
        "$FASTERQ" --split-files --threads "$THREADS" \
            --temp "$SRA_DIR/$acc" --outdir "$READS_DIR" "$sra"
        ZIP="gzip"
        for _d in "${VB_ASM_CONDA_BIN:-$HOME/miniforge3/envs/asm_env/bin}" "$SRA_CONDA_BIN"; do
            if [ -x "$_d/pigz" ]; then ZIP="$_d/pigz -p $THREADS"; break; fi
        done
        [ "$ZIP" = "gzip" ] && command -v pigz >/dev/null 2>&1 && ZIP="pigz -p $THREADS"
        for fq in "$READS_DIR/$acc"*.fastq; do
            [ -e "$fq" ] || continue
            log "$acc: compress $(basename "$fq") [${ZIP%% *}]"
            $ZIP -f "$fq"
        done

        # Drop the .sra once its FASTQ exists. 
        if compgen -G "$READS_DIR/${acc}*.fastq.gz" >/dev/null; then
            rm -rf "$SRA_DIR/$acc/fasterq.tmp."*
            if [ -f "$sra" ]; then
                log "$acc: removing .sra ($(du -h "$sra" | cut -f1)) - FASTQ is on disk"
                rm -f "$sra"
            fi
        else
            echo "error: no FASTQ produced for $acc; keeping $sra" >&2
            exit 4
        fi
    done

    echo
    echo "=== $DATASET reads ==="
    n_files=0
    for fq in "$READS_DIR"/*.fastq.gz; do
        [ -e "$fq" ] || continue
        printf '  %-28s %s\n' "$(basename "$fq")" "$(du -h "$fq" | cut -f1)"
        n_files=$((n_files + 1))
    done
    [ "$n_files" -gt 0 ] || { echo "error: no FASTQ produced under $READS_DIR" >&2; exit 5; }
    n_acc=$(printf '%s\n' $ACCESSIONS | wc -l)
    echo "  -> $n_files file(s) for $n_acc accession(s)"
    if [ "$READ_TYPE" = "short" ] && [ "$n_files" -lt $((n_acc * 2)) ]; then
        echo "  warning: read_type=short but some runs produced a single file;" >&2
        echo "           check whether those runs are really paired-end." >&2
    fi
    if [ "$READ_TYPE" = "long" ] && [ "$n_files" -gt "$n_acc" ]; then
        echo "  warning: read_type=long but some runs produced mate pairs;" >&2
        echo "           check whether those runs are really single-end." >&2
    fi
}

# =============================================================================
case "$SOURCE" in
sra)
# =============================================================================
    fetch_reads
    echo "phase 1 done -> $READS_DIR"
    ;;

# =============================================================================
assembly)
# =============================================================================
    # Resolve the assembly directory: absolute, or relative to Data/Real_dataset.
    [ "$ASSEMBLY_DIR" = "-" ] && {
        echo "error: source=assembly needs an assembly_dir in the registry" >&2; exit 2; }
    case "$ASSEMBLY_DIR" in
        /*) ASM="$ASSEMBLY_DIR" ;;
        *)  ASM="$REPO_HOME/Data/Real_dataset/$ASSEMBLY_DIR" ;;
    esac
    [ -d "$ASM" ] || { echo "error: assembly directory not found: $ASM" >&2; exit 2; }

    mkdir -p "$MWORK"
    echo "  assembly   : $ASM"
    echo "  work dir   : $MWORK"
    [ "$IMPORT_FROM" != "-" ] && echo "  adopting   : $IMPORT_FROM"
    echo

    # place <source-file> <destination-name> [required]
    # Copies only when absent or stale. These are hundreds of MB, so an
    # identical file already in place is left alone rather than recopied.
    place() {
        local src="$1" dst="$MWORK/$2" required="${3:-yes}"
        if [ ! -f "$src" ]; then
            [ "$required" = "yes" ] && {
                echo "error: required input not found: $src" >&2; exit 4; }
            return 0
        fi
        if [ -f "$dst" ] && [ "$dst" -nt "$src" ] \
           && [ "$(stat -c%s "$dst")" = "$(stat -c%s "$src")" ]; then
            log "[skip] $2 already in place"
            return 0
        fi
        log "copy $2  ($(du -h "$src" | cut -f1))"
        cp -f "$src" "$dst"
    }

    # ---- the assembly itself -------------------------------------------------
    # Assemblers name the GFA differently; take whichever this one wrote.
    gfa=""
    for candidate in assembly_graph_with_scaffolds.gfa assembly_graph.gfa \
                     assembly_graph_after_simplification.gfa final.contigs.gfa; do
        [ -f "$ASM/$candidate" ] && { gfa="$ASM/$candidate"; break; }
    done
    [ -n "$gfa" ] || { echo "error: no GFA found in $ASM" >&2; exit 4; }
    log "GFA: $(basename "$gfa")"

    contigs=""
    for candidate in $CONTIG_CANDIDATES; do
        [ -f "$ASM/$candidate" ] && { contigs="$ASM/$candidate"; break; }
    done
    if [ -z "$contigs" ]; then
        echo "error: no contig FASTA in $ASM" >&2
        echo "       looked for: $CONTIG_CANDIDATES" >&2
        exit 4
    fi
    log "contigs: $(basename "$contigs")"
    place "$contigs" "contigs.fasta"
    place "$gfa" "assembly_graph.gfa"
    if [ -n "$PATH_TABLE" ]; then
        place "$ASM/$PATH_TABLE" "$PATH_TABLE"
    else
        log "no path table for $ASSEMBLER; GFA links will be matched by contig name"
    fi

    # ---- adopt an earlier build's finished artifacts -------------------------
    if [ "$IMPORT_FROM" != "-" ]; then
        [ -d "$IMPORT_FROM" ] || {
            echo "error: import_from directory not found: $IMPORT_FROM" >&2; exit 4; }

        # The contig set the earlier build settled on, and its coverage matrix.
        # Both are optional: phase 2 recomputes whatever is missing.
        place "$IMPORT_FROM/contigs_filt.fasta" "contigs_filt.fasta" no
        place "$IMPORT_FROM/coverage.csv" "coverage.csv" no

        # geNomad. The whole output tree is copied into the work directory, not
        # just the summary phase 2 reads: it is the expensive artifact here
        # (hours on a 44k-contig assembly) and keeping the annotation and
        # provirus tables beside it means the call can be re-examined later
        # without the original run directory.
        if [ -d "$IMPORT_FROM/viral/genomad" ]; then
            if [ -d "$MWORK/viral/genomad" ]; then
                log "[skip] geNomad output already in place"
            else
                log "copy geNomad output ($(du -sh "$IMPORT_FROM/viral/genomad" | cut -f1))"
                mkdir -p "$MWORK/viral"
                cp -r "$IMPORT_FROM/viral/genomad" "$MWORK/viral/genomad"
            fi
            # Phase 2 keys off viral_contigs.txt; derive it from the summary so
            # the adopted call is explicit rather than implied by a file name.
            summary="$(find "$MWORK/viral/genomad" -name "*_virus_summary.tsv" | head -1)"
            if [ -n "$summary" ]; then
                # column 1 is seq_name; a provirus is "<contig>|provirus_x_y",
                # and the assembly graph only knows the host contig.
                awk -F'\t' 'NR > 1 {split($1, a, "|"); print a[1]}' "$summary" \
                    | awk '!seen[$0]++' > "$MWORK/viral_contigs.txt"
                log "viral_contigs.txt: $(wc -l < "$MWORK/viral_contigs.txt") contigs from $(basename "$summary")"
            else
                echo "  warning: no *_virus_summary.tsv under $MWORK/viral/genomad;" >&2
                echo "           phase 2 will run geNomad itself." >&2
            fi
        fi
    fi

    # Coverage needs reads mapped to THIS assembly.  If none was adopted and the
    # registry lists real accessions, fetch them now; phase 2 maps them to the
    # imported contigs.  Without either, phase 2 stops - a viral binner with no
    # differential coverage is a different experiment, not a weaker one.
    if [ ! -f "$MWORK/coverage.csv" ]; then
        if [ "$ACCESSIONS" != "-" ]; then
            echo
            echo "  no coverage adopted; fetching reads to map against this assembly"
            fetch_reads
        else
            echo
            echo "  NOTE: no coverage adopted and no accessions in the registry."
            echo "        The assembly is in place, but phase 2 will stop: coverage"
            echo "        needs reads mapped to it. Add the SRA runs to the"
            echo "        accessions column, or point import_from at a build that"
            echo "        already has coverage.csv."
        fi
    fi

    echo
    echo "=== $DATASET inputs ==="
    for f in contigs.fasta assembly_graph.gfa ${PATH_TABLE:+"$PATH_TABLE"} \
             contigs_filt.fasta coverage.csv viral_contigs.txt; do
        if [ -f "$MWORK/$f" ]; then
            printf '  [ok]   %-24s %s\n' "$f" "$(du -h "$MWORK/$f" | cut -f1)"
        else
            printf '  [--]   %-24s (phase 2 will build it)\n' "$f"
        fi
    done
    [ -d "$MWORK/viral/genomad" ] && \
        printf '  [ok]   %-24s %s\n' "viral/genomad/" "$(du -sh "$MWORK/viral/genomad" | cut -f1)"
    echo "phase 1 done -> $MWORK"
    ;;

*)
    echo "error: source must be sra|assembly (got '$SOURCE')" >&2
    exit 2
    ;;
esac
