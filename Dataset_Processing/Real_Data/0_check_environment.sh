#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Two dirname levels: this script lives in Dataset_Processing/Real_Data/.
REPO_HOME="$(dirname "$(dirname "$SCRIPT_DIR")")"
REGISTRY="$SCRIPT_DIR/datasets.tsv"

[ -f "$SCRIPT_DIR/../env.sh" ] && . "$SCRIPT_DIR/../env.sh"
CONDA_BIN="${VB_CONDA_BIN:-$HOME/miniforge3/envs/viralbin/bin}"            # mmseqs, checkv, python
ASM_CONDA_BIN="${VB_ASM_CONDA_BIN:-$HOME/miniforge3/envs/asm_env/bin}"     # metaspades.py, minimap2, samtools, seqkit
SRA_CONDA_BIN="${VB_SRA_CONDA_BIN:-$HOME/miniforge3/envs/sra_env/bin}"     # prefetch, fasterq-dump
GENOMAD_CONDA_BIN="${VB_GENOMAD_CONDA_BIN:-$HOME/miniforge3/envs/genomad_env/bin}"  # genomad
MEGAHIT_CONDA_BIN="${VB_MEGAHIT_CONDA_BIN:-$HOME/miniforge3/envs/megahit_env/bin}"  # megahit
PY="$CONDA_BIN/python"
[ -x "$PY" ] || PY="python3"

# CheckV needs its database; geNomad needs its own. 
CHECKV_DB="${CHECKVDB:-$HOME/checkv-db/checkv-db-v1.5}"
GENOMAD_DB="${GENOMAD_DB:-$HOME/db/genomad/genomad_db}"

missing=0
warned=0

check_bin() {
    local name="$1" bin="$2" hint="$3"
    local found=""
    for dir in "$CONDA_BIN" "$ASM_CONDA_BIN" "$SRA_CONDA_BIN" \
               "$GENOMAD_CONDA_BIN" "$MEGAHIT_CONDA_BIN"; do
        [ -x "$dir/$bin" ] && { found="$dir/$bin"; break; }
    done
    [ -z "$found" ] && command -v "$bin" >/dev/null 2>&1 && found="$(command -v "$bin")"
    if [ -n "$found" ]; then
        printf '  [ok]   %-12s -> %s\n' "$name" "$found"
    else
        printf '  [MISS] %-12s (not in the search path above, nor on PATH)\n' "$name"
        printf '         install: %s\n' "$hint"
        missing=$((missing + 1))
    fi
}

echo "=== dataset registry ($REGISTRY) ==="
if [ -f "$REGISTRY" ]; then
    "$PY" - "$REGISTRY" <<'PY' || missing=$((missing + 1))
import sys
rows, seen = [], set()
with open(sys.argv[1]) as fh:
    header = None
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))
bad = 0
import os
for r in rows:
    accs = [a for a in r["accessions"].split(",") if a]
    problems = []
    if r["read_type"] not in ("short", "long"):
        problems.append(f"read_type={r['read_type']!r} (use short|long)")
    if r["viral_enriched"] not in ("yes", "no"):
        problems.append(f"viral_enriched={r['viral_enriched']!r} (use yes|no)")
    source = r.get("source", "sra")
    if source not in ("sra", "assembly"):
        problems.append(f"source={source!r} (use sra|assembly)")
    if not accs:
        problems.append("no accessions")
    if r["dataset"] in seen:
        problems.append("duplicate dataset name")
    seen.add(r["dataset"])
    if source == "assembly":
        # The assembly must exist before anything else is worth checking.
        d = r.get("assembly_dir", "-")
        if d == "-":
            problems.append("source=assembly needs an assembly_dir")
        else:
            full = d if d.startswith("/") else os.path.join(
                os.path.dirname(os.path.dirname(os.path.dirname(
                    os.path.abspath(sys.argv[1])))), "Data/Real_dataset", d)
            if not os.path.isdir(full):
                problems.append(f"assembly_dir not found: {full}")
        imp = r.get("import_from", "-")
        if imp != "-" and not os.path.isdir(imp):
            problems.append(f"import_from not found: {imp}")
    # read_type picks the assembler only when this pipeline assembles; a
    # pre-built assembly may declare its own, because the assembly's format and
    # the reads' type are independent.
    declared = r.get("assembler", "-")
    asm = (declared if declared not in ("-", "")
           else ("metaspades" if r["read_type"] == "short" else "metaflye"))
    if asm not in ("metaspades", "metaflye", "megahit"):
        problems.append(f"assembler={asm!r} (use metaspades|metaflye|megahit|-)")
    tag = "noviral" if r["viral_enriched"] == "yes" else "genomad"
    if problems:
        print(f"  [BAD]  {r['dataset']}: " + "; ".join(problems))
        bad += 1
    else:
        print(f"  [ok]   {r['dataset']:16} {source:8} {len(accs)} run(s), "
              f"{r['read_type']:5} -> {r['dataset']}__coasm__{asm}__{tag}")
        if len(accs) < 3:
            print(f"         note: {len(accs)} sample(s) is thin differential "
                  "coverage; 3+ is preferable")
        if source == "assembly" and r.get("import_from", "-") != "-":
            print(f"         adopts artifacts from {r['import_from']}")
print(f"  {len(rows)} dataset(s) in registry")
sys.exit(1 if bad or not rows else 0)
PY
else
    printf '  [MISS] %s\n' "$REGISTRY"
    missing=$((missing + 1))
fi

echo
echo "=== tools ==="
for dir in "$CONDA_BIN" "$ASM_CONDA_BIN" "$SRA_CONDA_BIN" \
           "$GENOMAD_CONDA_BIN" "$MEGAHIT_CONDA_BIN"; do
    printf '  search: %s%s\n' "$dir" "$([ -d "$dir" ] || echo '   (missing)')"
done
# run_all.sh exports only the assemblers its selected datasets need, so a
# short-read-only run is not blocked by a missing long-read assembler.
REQUIRED_ASMS="${VB_REQUIRED_ASMS:-metaspades metaflye}"
echo "  (checking assemblers:$REQUIRED_ASMS)"
for asm in $REQUIRED_ASMS; do
    case "$asm" in
        metaspades) check_bin "metaspades" "metaspades.py" \
            "bash ../setup_environments.sh --only asm_env" ;;
        metaflye) check_bin "metaflye" "flye" \
            "bash ../setup_environments.sh --only viralbin" ;;
        megahit) check_bin "megahit" "megahit" \
            "bash ../setup_environments.sh --only megahit_env" ;;
        *) printf '  [MISS] unknown assembler %s\n' "$asm"; missing=$((missing + 1)) ;;
    esac
done
check_bin "prefetch"     "prefetch"      "bash ../setup_environments.sh --only sra_env"
check_bin "fasterq-dump" "fasterq-dump"  "bash ../setup_environments.sh --only sra_env"
check_bin "minimap2"     "minimap2"      "bash ../setup_environments.sh --only asm_env"
check_bin "samtools"     "samtools"      "bash ../setup_environments.sh --only asm_env"
check_bin "seqkit"       "seqkit"        "bash ../setup_environments.sh --only asm_env"
check_bin "mmseqs"       "mmseqs"        "bash ../setup_environments.sh --only viralbin"
check_bin "checkv"       "checkv"        "bash ../setup_environments.sh --only viralbin"

echo
echo "=== databases ==="
if [ -d "$CHECKV_DB" ] && [ -d "$CHECKV_DB/genome_db" ] && [ -d "$CHECKV_DB/hmm_db" ]; then
    printf '  [ok]   checkv-db    -> %s\n' "$CHECKV_DB"
else
    printf '  [MISS] checkv-db    (looked at %s)\n' "$CHECKV_DB"
    printf '         install: bash ../setup_environments.sh --databases\n'
    printf '         or set CHECKVDB in env.sh (see ../env.sh.example)\n'
    missing=$((missing + 1))
fi
# geNomad is only needed for bulk-metagenome datasets (viral_enriched=no), so
# both its binary and its database are WARNINGS rather than failures while every
# dataset in the registry is a virome. Phase 2 fails loudly if a dataset that
# does need geNomad is built without it.
GENOMAD_BIN=""
for dir in "$CONDA_BIN" "$GENOMAD_CONDA_BIN"; do
    [ -x "$dir/genomad" ] && { GENOMAD_BIN="$dir/genomad"; break; }
done
[ -z "$GENOMAD_BIN" ] && command -v genomad >/dev/null 2>&1 && GENOMAD_BIN="$(command -v genomad)"
if [ -n "$GENOMAD_BIN" ]; then
    printf '  [ok]   genomad      -> %s\n' "$GENOMAD_BIN"
else
    printf '  [warn] genomad      (not installed)\n'
    printf '         Only needed for datasets with viral_enriched=no. Install:\n'
    printf '           bash ../setup_environments.sh --only genomad_env\n'
    warned=$((warned + 1))
fi
if [ -d "$GENOMAD_DB" ]; then
    printf '  [ok]   genomad-db   -> %s\n' "$GENOMAD_DB"
else
    printf '  [warn] genomad-db   (looked at %s)\n' "$GENOMAD_DB"
    printf '         install: bash ../setup_environments.sh --databases\n'
    printf '         or set GENOMAD_DB in env.sh (see ../env.sh.example)\n'
    warned=$((warned + 1))
fi

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
    print("         install: bash ../setup_environments.sh --only viralbin")
    print("         (or, into the existing env: pip install "
          + " ".join(pip_name.get(m, m) for m in bad) + ")")
    sys.exit(1)
PY
else
    echo "  [MISS] no python found at $PY"
    missing=$((missing + 1))
fi

echo
echo "=== pipeline scripts ($SCRIPT_DIR) ==="
for script in \
    "1_get_inputs.sh" \
    "2_assemble_and_build_graph.py" \
    "3_acceptance_check.py" \
    "4_checkv_evaluate.py"; do
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
if [ "$warned" -gt 0 ]; then
    echo "preflight OK with $warned warning(s) - fine for the virome datasets in the registry."
else
    echo "preflight OK - registry, tools and databases all resolved."
fi
