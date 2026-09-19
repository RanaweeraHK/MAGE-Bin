#!/usr/bin/env bash
# Preflight: registry, tools, download endpoints, free disk.
#
#   bash 0_check_environment.sh
#
# Unlike ../Real_Data this needs no sra-tools (phase 1 uses the ENA mirror) and
# no geNomad (every registered mock is viral_enriched=yes). Short rows assemble
# with megahit, long rows with metaflye.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Two dirname levels: this script lives in Dataset_Processing/Mock_Data/.
REPO_HOME="$(dirname "$(dirname "$SCRIPT_DIR")")"
REGISTRY="$SCRIPT_DIR/datasets.tsv"

[ -f "$SCRIPT_DIR/../env.sh" ] && . "$SCRIPT_DIR/../env.sh"
CONDA_BIN="${VB_CONDA_BIN:-$HOME/miniforge3/envs/viralbin/bin}"          # mmseqs, flye, python
ASM_CONDA_BIN="${VB_ASM_CONDA_BIN:-$HOME/miniforge3/envs/asm_env/bin}"   # minimap2, samtools, seqkit
MEGAHIT_CONDA_BIN="${VB_MEGAHIT_CONDA_BIN:-$HOME/miniforge3/envs/megahit_env/bin}"  # megahit
PY="$CONDA_BIN/python"
[ -x "$PY" ] || PY="python3"

missing=0
warned=0

check_bin() {
    local name="$1" bin="$2" hint="$3"
    local found=""
    for dir in "$CONDA_BIN" "$ASM_CONDA_BIN" "$MEGAHIT_CONDA_BIN"; do
        [ -x "$dir/$bin" ] && { found="$dir/$bin"; break; }
    done
    [ -z "$found" ] && command -v "$bin" >/dev/null 2>&1 && found="$(command -v "$bin")"
    if [ -n "$found" ]; then
        printf '  [ok]   %-14s -> %s\n' "$name" "$found"
    else
        printf '  [MISS] %-14s (not in the search path above, nor on PATH)\n' "$name"
        printf '         install: %s\n' "$hint"
        missing=$((missing + 1))
    fi
}

echo "=== dataset registry ($REGISTRY) ==="
if [ -f "$REGISTRY" ]; then
    "$PY" - "$REGISTRY" <<'PY' || missing=$((missing + 1))
import sys

REQUIRED = ["dataset", "study", "accessions", "read_type", "viral_enriched",
            "source", "assembly_dir", "import_from", "assembler",
            "mean_read_length", "notes", "provider", "source_id",
            "reads_format", "reference_set"]
KNOWN_REFERENCE_SETS = {"dsmz_115_molecules", "phage_mock_15"}

header, rows = None, []
with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))

bad = 0
if header != REQUIRED:
    print(f"  [BAD]  header is {header}\n         expected {REQUIRED}")
    bad += 1
seen = set()
for r in rows:
    name = r.get("dataset", "?")
    problems = []
    if name in seen:
        problems.append("duplicate dataset name")
    seen.add(name)
    if r.get("read_type") not in ("short", "long"):
        problems.append(f"read_type={r.get('read_type')!r} (use short|long)")
    if r.get("viral_enriched") not in ("yes", "no"):
        problems.append(f"viral_enriched={r.get('viral_enriched')!r} (use yes|no)")
    if r.get("provider") not in ("ena", "dataverse"):
        problems.append(f"provider={r.get('provider')!r} (use ena|dataverse)")
    if r.get("reads_format") not in ("fastq", "fasta"):
        problems.append(f"reads_format={r.get('reads_format')!r} (use fastq|fasta)")
    if r.get("reference_set") not in KNOWN_REFERENCE_SETS:
        problems.append(f"reference_set={r.get('reference_set')!r} has no builder "
                        f"in 1_get_inputs.py (known: {sorted(KNOWN_REFERENCE_SETS)})")
    libraries = [a for a in r.get("accessions", "").split(",") if a]
    if not libraries:
        problems.append("no libraries in `accessions`")
    elif len(libraries) < 3:
        print(f"  [warn] {name}: {len(libraries)} librarie(s) - fewer than 3 "
              f"coverage columns is thin")
    if problems:
        bad += 1
        print(f"  [BAD]  {name}: " + "; ".join(problems))
    else:
        print(f"  [ok]   {name}: {len(libraries)} librarie(s), "
              f"{r['read_type']} reads from {r['provider']}, "
              f"truth = {r['reference_set']}")
sys.exit(1 if bad else 0)
PY
else
    echo "  [MISS] $REGISTRY not found"
    missing=$((missing + 1))
fi

echo
echo "=== tools ==="
check_bin "minimap2" minimap2 "conda create -n asm_env -c bioconda -c conda-forge minimap2 samtools seqkit"
check_bin "samtools" samtools "conda install -n asm_env -c bioconda samtools"
check_bin "seqkit"   seqkit   "conda install -n asm_env -c bioconda seqkit"
check_bin "megahit"  megahit  "conda create -n megahit_env -c bioconda megahit   # read_type=short rows"
check_bin "megahit_toolkit" megahit_toolkit "ships with megahit; without it no GFA is built"
check_bin "flye"     flye     "conda install -n viralbin -c bioconda flye            # read_type=long rows"
check_bin "mmseqs"   mmseqs   "conda install -n viralbin -c bioconda mmseqs2"

echo
echo "=== python (the viralbin env) ==="
"$PY" - <<'PY' || missing=$((missing + 1))
import importlib, sys
need = ["numpy", "pandas", "torch", "torch_geometric", "pysam", "sklearn",
        "networkx", "pyrodigal_gv"]
bad = 0
for module in need:
    try:
        importlib.import_module(module)
        print(f"  [ok]   {module}")
    except Exception as exc:                      # noqa: BLE001 - reporting only
        print(f"  [MISS] {module}: {exc}")
        bad += 1
print(f"  python -> {sys.executable}")
sys.exit(1 if bad else 0)
PY

echo
echo "=== download endpoints (phase 1) ==="
for url in "https://www.ebi.ac.uk/ena/portal/api/" \
           "https://entrepot.recherche.data.gouv.fr/api/info/version" \
           "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"; do
    # --retry: one refused connection is not a broken environment.
    if curl -sSf -m 30 --retry 2 --retry-delay 3 -o /dev/null "$url" 2>/dev/null; then
        printf '  [ok]   %s\n' "${url#https://}"
    else
        printf '  [warn] %s unreachable - phase 1 cannot download from it\n' "${url#https://}"
        warned=$((warned + 1))
    fi
done

echo
echo "=== disk ==="
# Reads are ~17 GB; assembly scratch is several times that.
avail_kb="$(df -Pk "$REPO_HOME" | awk 'NR==2 {print $4}')"
avail_gb=$((avail_kb / 1024 / 1024))
printf '  %s GB free under %s\n' "$avail_gb" "$REPO_HOME"
if [ "$avail_gb" -lt 100 ]; then
    printf '  [warn] under 100 GB - reads are ~17 GB and assembly scratch is several times that\n'
    warned=$((warned + 1))
fi

echo
if [ "$missing" -gt 0 ]; then
    echo "$missing problem(s) must be fixed before running run_all.sh"
    exit 1
fi
[ "$warned" -gt 0 ] && echo "$warned warning(s); nothing blocking"
echo "environment OK"
