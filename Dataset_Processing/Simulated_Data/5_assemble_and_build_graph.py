# ============================================================================
# DATASET PREPARATION PIPELINE
# ============================================================================
#
# Input sequencing reads
#        ↓
# [1] Detect samples
#     function: detect_reads()
#        ↓
# [2] Co-assemble all samples
#     function: assemble()
#       ├── _assemble_metaspades()   # short reads
#       └── _assemble_metaflye()     # long reads
#        ↓
# [3] Keep contigs ≥ 1,000 bp
#     function: length_filter()
#        ↓
# [4] Map reads back to contigs → calculate coverage
#     functions:
#       ├── map_samples()
#       ├── compute_coverage_matrix()
#       └── remove_bams()
#        ↓
# [5] Map contigs to known reference → ground-truth labels
#     function: label_from_reference()
#        ↓
# [6] Calculate TNF + coverage node features
#     function: build_graph()
#       ├── build_kmer_index()
#       ├── tnf_vector()
#       └── zscore()
#        ↓
# [7] Convert GFA segment graph → CONTIG-level graph
#     function: build_graph()
#       ├── parse_gfa_rich()
#       ├── parse_spades_paths_oriented() / parse_flye_paths_oriented()
#       └── contig_graph_from_paths()
#        ↓
# [8] Save PyTorch graph + contig metadata
#     functions:
#       ├── build_graph()            # saves viral_graph.pt
#       │                            # and assembly_heterograph.pt
#       └── write_contig_metadata() # saves contig_metadata.tsv
#        ↓
# [9] Predict ORFs + cluster proteins
#     function: extract_viral_features()
#       ├── _call_orfs()             # pyrodigal-gv
#       └── _cluster_proteins()      # MMseqs2
#        ↓
# Processed dataset for the viral-binning model
# ============================================================================


from __future__ import annotations

import argparse
import datetime
import gzip
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam
import torch
from sklearn.neighbors import kneighbors_graph
from torch_geometric.data import Data, HeteroData


def _load_env_sh() -> None:
    env_sh = Path(__file__).resolve().parent.parent / "env.sh"
    if not env_sh.exists():
        return
    for line in env_sh.read_text().splitlines():
        line = line.strip()
        if line.startswith("export "):
            line = line[len("export "):]
        key, sep, value = line.partition("=")
        key = key.strip()
        if not sep or key.startswith("#"):
            continue
        if not (key.startswith("VB_") or key in ("CHECKVDB", "GENOMAD_DB")):
            continue
        value = value.split("#")[0].strip().strip('"').strip("'")
        # Only the `${VAR:-default}` form the template uses is expanded; any
        # other shell syntax is left alone rather than half-interpreted.
        if value.startswith("${") and value.endswith("}") and ":-" in value:
            value = value[2:-1].split(":-", 1)[1]
        value = os.path.expandvars(os.path.expanduser(value))
        if value:
            os.environ.setdefault(key, value)

_load_env_sh()

TIER = os.environ.get("VB_TIER", "high")
if TIER not in ("low", "medium", "high"):
    raise SystemExit(f"error: VB_TIER must be low|medium|high (got {TIER!r})")
TIER_SUFFIX = {"low": "_low", "medium": "_medium", "high": "_high"}[TIER]

SIM = os.environ.get("VB_SIM", "iss")
ASSEMBLER = os.environ.get("VB_ASM", "metaspades")
SIM_READ_TYPE = {"iss": "short", "badread": "long"}
ASM_READ_TYPE = {"metaspades": "short", "metaflye": "long"}
if SIM not in SIM_READ_TYPE:
    raise SystemExit(f"error: VB_SIM must be one of "
                     f"{'|'.join(SIM_READ_TYPE)} (got {SIM!r})")
if ASSEMBLER not in ASM_READ_TYPE:
    raise SystemExit(f"error: VB_ASM must be one of "
                     f"{'|'.join(ASM_READ_TYPE)} (got {ASSEMBLER!r})")
READ_TYPE = SIM_READ_TYPE[SIM]
if ASM_READ_TYPE[ASSEMBLER] != READ_TYPE:
    raise SystemExit(
        f"error: {ASSEMBLER} assembles {ASM_READ_TYPE[ASSEMBLER]} reads but "
        f"{SIM} produces {READ_TYPE} reads. Valid pairs: "
        "iss+metaspades (short), badread+metaflye (long).")
# minimap2 preset for mapping reads back to contigs: short-read vs ONT.
MINIMAP_PRESET = {"short": "sr", "long": "map-ont"}[READ_TYPE]

HERE = Path(__file__).resolve().parent
REPO_HOME = HERE.parent.parent
WORK_DIR = REPO_HOME / f"Data/work/multibiome{TIER_SUFFIX}"
POOL_DIR = WORK_DIR / "genome_pool"
# The genome pool and abundance tables are simulator-independent and shared by
# every simulator for a tier; only the reads differ. iss keeps the unsuffixed
# directory so the three pre-existing datasets are untouched. (Phase 4 applies
# the identical rule - keep the two in step.)
READS_DIR = WORK_DIR / ("reads" if SIM == "iss" else f"reads_{SIM}")
OUT_ROOT = REPO_HOME / "Data/Processed_data"

SIM_TAG = "" if SIM == "iss" else f"__{SIM}"
RUN_NAME = f"multibiome{TIER_SUFFIX}{SIM_TAG}__coasm__{ASSEMBLER}__noviral"
REFERENCE_FASTA = POOL_DIR / "expanded_genomes.fasta"

CONDA_BIN = os.environ.get(
    "VB_CONDA_BIN", os.path.expanduser("~/miniforge3/envs/viralbin/bin"))
ASM_CONDA_BIN = os.environ.get(
    "VB_ASM_CONDA_BIN", os.path.expanduser("~/miniforge3/envs/asm_env/bin"))

os.environ["PATH"] = os.pathsep.join(
    [ASM_CONDA_BIN, CONDA_BIN, os.environ.get("PATH", "")])

THREADS = int(os.environ.get("VB_THREADS", 8))
MEMORY_GB = int(os.environ.get("VB_MEMORY_GB", 10))
MIN_CONTIG_LEN = 1000
TNF_K = 4                        # canonical tetranucleotide frequency width
KNN_K = 5                        # unused while USE_KNN_EDGES is False
USE_GFA_EDGES = True
USE_KNN_EDGES = False
INCLUDE_SHARED_SEGMENT_EDGES = False
MIN_MAPQ = 0


# always produces: coassembly + metaspades + no viral tools.
MWORK = WORK_DIR / "standardize_work" / RUN_NAME  # scratch for this run
MOUT = OUT_ROOT / RUN_NAME                        # final artifacts (kept)


# ============================================================================
# Tool discovery + subprocess helpers 
# ============================================================================
def _fmt_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m{s:02d}s"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


class _Heartbeat:
    """Print an 'elapsed' line every couple of minutes while a tool runs.

    External assemblers/aligners expose no percentage or ETA, so this at
    least proves the run is alive and shows how long it has been going.
    """

    def __init__(self, label: str, interval: float = 120.0):
        self.label = label
        self.interval = interval
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._start = time.monotonic()

    def _loop(self):
        while not self._stop.wait(self.interval):
            elapsed = _fmt_elapsed(time.monotonic() - self._start)
            print(f"    ... {self.label} still running - elapsed {elapsed}",
                  flush=True)

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()


class ToolError(RuntimeError):
    """Raised when a required external tool is missing or exits non-zero."""


@dataclass
class Toolbox:
    """Resolves an external binary: --conda-bin first, then PATH.

    Tries both ``name`` and ``name.py`` in each location, since SPAdes-family
    wrappers (metaspades.py, ...) install under the .py name.
    """

    conda_bin: str | None = None

    def resolve(self, name: str, required: bool = True) -> str | None:
        if self.conda_bin:
            for candidate in (name, f"{name}.py"):
                path = os.path.join(self.conda_bin, candidate)
                if os.path.isfile(path) and os.access(path, os.X_OK):
                    return path
        for candidate in (name, f"{name}.py"):
            found = shutil.which(candidate)
            if found:
                return found
        if required:
            hint = f" (looked in {self.conda_bin!r} and PATH)" if self.conda_bin else ""
            raise ToolError(f"Required tool not found: {name}{hint}")
        return None

    def python(self) -> str:
        if self.conda_bin:
            for candidate in ("python", "python3"):
                path = os.path.join(self.conda_bin, candidate)
                if os.path.isfile(path):
                    return path
        return sys.executable


def run(cmd, *, cwd=None, env=None, log_prefix="", heartbeat=True) -> None:
    """Run ``cmd`` (a list), streaming output; raise ToolError on failure."""
    printable = " ".join(str(c) for c in cmd)
    print(f"{log_prefix}$ {printable}", flush=True)
    label = (log_prefix.strip(" []") or os.path.basename(str(cmd[0]))).strip()
    start = time.monotonic()
    if heartbeat:
        with _Heartbeat(label):
            result = subprocess.run(cmd, cwd=cwd, env=env)
    else:
        result = subprocess.run(cmd, cwd=cwd, env=env)
    if result.returncode != 0:
        raise ToolError(f"Command failed ({result.returncode}): {printable}")
    print(f"{log_prefix}done in {_fmt_elapsed(time.monotonic() - start)}",
          flush=True)


def run_piped(producer, consumer, *, out_path=None, log_prefix="") -> None:
    """Run ``producer | consumer``, optionally redirecting stdout to a file.

    Used for ``minimap2 ... | samtools sort`` so the SAM stream is never
    materialised on disk.
    """
    print(f"{log_prefix}$ {' '.join(producer)} | {' '.join(consumer)}"
          + (f" > {out_path}" if out_path else ""), flush=True)
    p = subprocess.Popen(producer, stdout=subprocess.PIPE)
    out = open(out_path, "wb") if out_path else None
    try:
        c = subprocess.Popen(consumer, stdin=p.stdout, stdout=out)
        p.stdout.close()  # allow producer to receive SIGPIPE
        c_rc = c.wait()
        p_rc = p.wait()
    finally:
        if out:
            out.close()
    if p_rc != 0 or c_rc != 0:
        raise ToolError(
            f"Piped command failed (producer={p_rc}, consumer={c_rc})")

# measure actual average read length from FASTAQ
def measure_read_length(det: Detection, sample_reads: int = 20000) -> float:
    path = det.samples[0].r1
    opener = gzip.open if is_gzip(path) else open
    total = count = 0
    with opener(path, "rt") as fh:
        for index, line in enumerate(fh):
            if index % 4 == 1:                      # FASTQ sequence lines
                total += len(line.strip())
                count += 1
                if count >= sample_reads:
                    break
    return (total / count) if count else float("nan")


def is_gzip(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return path.endswith((".gz", ".bgz"))


def fasta_ids(path) -> list[str]:
    """First-token id of every FASTA header, in file order. Used throughout
    (coverage/labels/graph all need the contig id list)."""
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


# ============================================================================
# Read discovery + pairing (formerly detect.py)
# ============================================================================
_FASTQ_EXT = (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fastq.bgz", ".fq.bgz")
_MATE_PATTERNS = [
    re.compile(r"^(?P<stem>.+?)[._]R(?P<mate>[12])(?:_\d+)?$"),  # _R1 / .R1 / _R1_001
    re.compile(r"^(?P<stem>.+?)[._](?P<mate>[12])$"),            # _1 / .1  (ENA)
]


@dataclass
class Sample:
    """One sequencing sample: a stem plus its mate file(s)."""

    name: str                    # sanitised sample id, e.g. "sample0"
    stem: str                    # original shared prefix
    r1: str
    r2: str | None = None


@dataclass
class Detection:
    """Result of scanning the reads directory."""

    dataset_dir: str
    samples: list[Sample] = field(default_factory=list)
    single_end: bool = False

    @property
    def n_samples(self) -> int:
        return len(self.samples)


def _strip_ext(filename: str) -> tuple[str, str] | None:
    lower = filename.lower()
    for ext in _FASTQ_EXT:
        if lower.endswith(ext):
            return filename[: -len(ext)], filename[-len(ext):]
    return None


def _sanitise(stem: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", stem)


def detect_reads(dataset_dir) -> Detection:
    """Scan ``dataset_dir`` for FASTQs and pair R1/R2 by shared stem."""
    dataset_dir = str(dataset_dir)
    if not os.path.isdir(dataset_dir):
        raise FileNotFoundError(f"reads directory not found: {dataset_dir}")

    fastqs = []
    for dirpath, dirnames, filenames in os.walk(dataset_dir):
        dirnames.sort()
        dirnames[:] = [d for d in dirnames if d not in ("work", "bam", "assembly")]
        for entry in sorted(filenames):
            stripped = _strip_ext(entry)
            if stripped:
                fastqs.append((os.path.join(dirpath, entry), dirpath, stripped[0]))

    if not fastqs:
        raise FileNotFoundError(
            f"no FASTQ files found under {dataset_dir} (searched recursively)")

    groups: dict[tuple[str, str], dict[str, str]] = {}
    for full, parent, base in fastqs:
        matched = None
        for pattern in _MATE_PATTERNS:
            m = pattern.match(base)
            if m:
                matched = (m.group("stem"), m.group("mate"))
                break
        if matched:
            stem, mate = matched
            group = groups.setdefault((parent, stem), {})
            if mate in group:
                raise ValueError(
                    "multiple FASTQ files resolve to the same sample/mate: "
                    f"{group[mate]!r}, {full!r}")
            group[mate] = full
        else:
            groups.setdefault((parent, base), {})["1"] = full

    samples: list[Sample] = []
    any_single = False
    for i, key in enumerate(sorted(groups, key=lambda k: (k[1], k[0]))):
        parent, stem = key
        mates = groups[key]
        r1 = mates.get("1")
        r2 = mates.get("2")
        if r1 is None:
            r1, r2 = mates["2"], None
        if r2 is None:
            any_single = True
        samples.append(Sample(name=f"sample{i}", stem=_sanitise(stem), r1=r1, r2=r2))

    return Detection(dataset_dir=dataset_dir, samples=samples, single_end=any_single)


# ============================================================================
# Assembly - metaSPAdes co-assembly 
# ============================================================================
def _assembly_signature(det: Detection) -> dict:
    """Inputs/parameters whose change invalidates a cached assembly."""
    reads = []
    for sample in det.samples:
        for role, path in (("r1", sample.r1), ("r2", sample.r2)):
            if path:
                stat = os.stat(path)
                reads.append({"sample": sample.stem, "role": role,
                              "path": os.path.abspath(path), "size": stat.st_size,
                              "mtime_ns": stat.st_mtime_ns})
    return {"assembler": ASSEMBLER, "threads": THREADS, "memory_gb": MEMORY_GB,
            "only_assembler": ASSEMBLER == "metaspades", "read_type": READ_TYPE,
            "minimap_preset": MINIMAP_PRESET, "reads": reads}


def _load_json(path: str) -> dict | None:
    try:
        with open(path) as fh:
            value = json.load(fh)
        return value if isinstance(value, dict) else None
    except (OSError, ValueError):
        return None


def _concat_gz(sources: list[str], dst: str) -> None:
    """Pool FASTQ inputs into one valid gzip stream (no decompressed copy)."""
    if os.path.isfile(dst) and os.path.getsize(dst) > 0:
        newest_src = max(os.path.getmtime(s) for s in sources)
        if os.path.getmtime(dst) >= newest_src:
            return
    if all(is_gzip(src) for src in sources):
        with open(dst, "wb") as out:
            for src in sources:
                with open(src, "rb") as fh:
                    shutil.copyfileobj(fh, out, length=1 << 20)
        return
    with gzip.open(dst, "wb") as out:
        for src in sources:
            opener = gzip.open if is_gzip(src) else open
            with opener(src, "rb") as fh:
                shutil.copyfileobj(fh, out, length=1 << 20)


def _pool_reads(det: Detection, mwork: str) -> list[str]:
    """Pool all 15 samples into one metaSPAdes paired-end library."""
    paired = [s for s in det.samples if s.r2]
    single = [s for s in det.samples if not s.r2]
    args: list[str] = []
    if paired:
        r1 = os.path.join(mwork, "all_R1.fastq.gz")
        r2 = os.path.join(mwork, "all_R2.fastq.gz")
        _concat_gz([s.r1 for s in paired], r1)
        _concat_gz([s.r2 for s in paired], r2)
        args += ["-1", r1, "-2", r2]
    if single:
        pooled = os.path.join(mwork, "all_single.fastq.gz")
        _concat_gz([s.r1 for s in single], pooled)
        args += ["-s", pooled]
    return args


def _free_output_dir(path: str) -> None:
    """Guarantee ``path`` does not exist (metaSPAdes aborts if it does)."""
    if not os.path.exists(path):
        return
    trash = f"{path}.trash"
    if os.path.exists(trash):
        shutil.rmtree(trash, ignore_errors=True)
    os.rename(path, trash)
    shutil.rmtree(trash, ignore_errors=True)

# copy only the files contigs.fasta, assembly_graph.gfa, contigs.paths
def _surface_spades_outputs(asm: str, mwork: str) -> None:
    shutil.copy(os.path.join(asm, "contigs.fasta"),
               os.path.join(mwork, "contigs.fasta"))
    for g in ("assembly_graph_with_scaffolds.gfa",
              "assembly_graph_after_simplification.gfa"):
        src = os.path.join(asm, g)
        if os.path.isfile(src) and os.path.getsize(src) > 0:
            shutil.copy(src, os.path.join(mwork, "assembly_graph.gfa"))
            break
    # Only contigs.paths is surfaced: build_graph() projects GFA links through
    # it. scaffolds.paths is never read by anything downstream.
    src = os.path.join(asm, "contigs.paths")
    if os.path.isfile(src) and os.path.getsize(src) > 0:
        shutil.copy(src, os.path.join(mwork, "contigs.paths"))


def _surface_flye_outputs(asm: str, mwork: str) -> None:
    shutil.copy(os.path.join(asm, "assembly.fasta"),
                os.path.join(mwork, "contigs.fasta"))
    for name in ("assembly_graph.gfa", "assembly_info.txt"):
        src = os.path.join(asm, name)
        if os.path.isfile(src) and os.path.getsize(src) > 0:
            shutil.copy(src, os.path.join(mwork, name))


def _assemble_metaspades(det: Detection, mwork: str, asm: str, tb: Toolbox) -> None:
    binary = tb.resolve("metaspades")
    read_args = _pool_reads(det, mwork)
    _free_output_dir(asm)  # metaSPAdes refuses to start if the dir exists
    run([binary, *read_args, "-t", str(THREADS), "-m", str(MEMORY_GB),
         "-o", asm, "--only-assembler"], log_prefix="[assemble] ")
    _surface_spades_outputs(asm, mwork)


def _assemble_metaflye(det: Detection, mwork: str, asm: str, tb: Toolbox) -> None:
    binary = tb.resolve("flye")
    read_args = _pool_reads(det, mwork)
    # _pool_reads returns metaSPAdes-style flags; Flye takes bare read paths.
    reads = [a for a in read_args if not a.startswith("-")]
    _free_output_dir(asm)
    run([binary, "--meta", "--nano-hq", *reads, "--out-dir", asm,
         "--threads", str(THREADS)], log_prefix="[assemble] ")
    _surface_flye_outputs(asm, mwork)


_ASSEMBLERS = {"metaspades": _assemble_metaspades, "metaflye": _assemble_metaflye}
# Where each assembler leaves its contigs, used to detect a cached assembly.
_ASSEMBLY_CONTIGS = {"metaspades": "contigs.fasta", "metaflye": "assembly.fasta"}


def assemble(det: Detection, mwork: str, tb: Toolbox, *, force: bool = False) -> None:
    os.makedirs(mwork, exist_ok=True)
    signature_path = os.path.join(mwork, "assembly_inputs.json")
    signature = _assembly_signature(det)
    rebuild = force or _load_json(signature_path) != signature
    if rebuild:
        for name in ("contigs.fasta", "assembly_graph.gfa", "contigs.paths",
                     "assembly_info.txt", "all_R1.fastq.gz", "all_R2.fastq.gz",
                     "all_single.fastq.gz"):
            path = os.path.join(mwork, name)
            if os.path.exists(path):
                os.remove(path)

    asm = os.path.join(mwork, "assembly")
    contigs_out = os.path.join(asm, _ASSEMBLY_CONTIGS[ASSEMBLER])
    if rebuild:
        _free_output_dir(asm)
    if os.path.isfile(contigs_out) and os.path.getsize(contigs_out) > 0:
        print(f"assembly already present: {contigs_out} (skipping {ASSEMBLER})")
        (_surface_flye_outputs if ASSEMBLER == "metaflye"
         else _surface_spades_outputs)(asm, mwork)
    else:
        _ASSEMBLERS[ASSEMBLER](det, mwork, asm, tb)

    temporary = signature_path + ".tmp"
    with open(temporary, "w") as fh:
        json.dump(signature, fh, indent=2)
    os.replace(temporary, signature_path)


def length_filter(mwork: str, tb: Toolbox, *, min_len: int) -> str:
    """seqkit length filter -> mwork/contigs_filt.fasta."""
    seqkit = tb.resolve("seqkit")
    src = os.path.join(mwork, "contigs.fasta")
    dst = os.path.join(mwork, "contigs_filt.fasta")
    with open(dst, "w") as fh:
        subprocess.run([seqkit, "seq", "-m", str(min_len), src], stdout=fh, check=True)
    run([seqkit, "stats", dst], log_prefix="[coverage] ")
    return dst


# ============================================================================
# Coverage - minimap2 mapping + per-contig mean depth 
# ============================================================================
def _ok(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _fresh(path: str, inputs: list[str]) -> bool:
    """True when a non-empty output is at least as new as every input."""
    return (_ok(path) and all(os.path.isfile(src) for src in inputs)
            and os.path.getmtime(path) >= max(os.path.getmtime(src) for src in inputs))


def map_samples(det: Detection, mwork: str, tb: Toolbox, *, force: bool = False) -> None:
    """minimap2 map each sample -> sorted+indexed BAM in mwork/bam."""
    minimap2 = tb.resolve("minimap2")
    samtools = tb.resolve("samtools")
    ref = os.path.join(mwork, "contigs_filt.fasta")
    bamdir = os.path.join(mwork, "bam")
    os.makedirs(bamdir, exist_ok=True)

    for s in det.samples:
        idx = int(s.name.replace("sample", ""))
        bam = os.path.join(bamdir, f"sample{idx}.bam")
        bai = bam + ".bai"
        reads = [s.r1] + ([s.r2] if s.r2 else [])
        inputs = [ref, *reads]
        if not force and _fresh(bam, inputs) and _fresh(bai, [bam]):
            print(f"[coverage] {bam} is current, skipping")
            continue
        for old in (bam, bai):
            if os.path.exists(old):
                os.remove(old)
        producer = [minimap2, "-ax", MINIMAP_PRESET, "-t", str(THREADS),
                    ref, *reads]
        consumer = [samtools, "sort", "-@", str(THREADS), "-o", bam, "-"]
        run_piped(producer, consumer, log_prefix="[coverage] ")
        run([samtools, "index", bam], log_prefix="[coverage] ")


def compute_coverage_matrix(mwork: str, n_samples: int) -> str:
    """Per-contig mean depth for each sample -> mwork/coverage.csv."""
    contigs = fasta_ids(os.path.join(mwork, "contigs_filt.fasta"))
    if not contigs:
        raise ValueError("contigs_filt.fasta is empty")

    def usable(read):
        return (not read.is_unmapped and not read.is_secondary
                and not read.is_supplementary and not read.is_duplicate
                and not read.is_qcfail and read.mapping_quality >= MIN_MAPQ)

    bam0 = os.path.join(mwork, "bam", "sample0.bam")
    with pysam.AlignmentFile(bam0) as bf:
        if list(bf.references) != contigs:
            raise ValueError(
                "sample0 BAM reference does not match contigs_filt.fasta; "
                "delete/rebuild stale BAM files")
        lengths = dict(zip(bf.references, bf.lengths))

    cov = np.zeros((len(contigs), n_samples), dtype=float)
    for s in range(n_samples):
        bam = os.path.join(mwork, "bam", f"sample{s}.bam")
        with pysam.AlignmentFile(bam) as bf:
            if list(bf.references) != contigs:
                raise ValueError(
                    f"sample{s} BAM reference does not match contigs_filt.fasta; "
                    "delete/rebuild stale BAM files")
            for i, c in enumerate(contigs):
                counts = bf.count_coverage(c, quality_threshold=0, read_callback=usable)
                depth = np.sum(counts, axis=0)
                cov[i, s] = depth.mean() if lengths[c] > 0 else 0.0
        print(f"[coverage] sample{s}: coverage done")

    out = os.path.join(mwork, "coverage.csv")
    with open(out, "w") as fh:
        fh.write("contig," + ",".join(f"sample{s}" for s in range(n_samples)) + "\n")
        for i, c in enumerate(contigs):
            fh.write(c + "," + ",".join(f"{cov[i, s]:.4f}" for s in range(n_samples)) + "\n")
    print(f"[coverage] wrote {out} ({len(contigs)} contigs x {n_samples} samples)")
    return out


def remove_bams(mwork: str) -> None:
    """Remove BAM/index files after coverage.csv was written successfully."""
    bamdir = os.path.join(mwork, "bam")
    if not os.path.isdir(bamdir):
        return
    removed = 0
    for name in os.listdir(bamdir):
        if name.endswith((".bam", ".bam.bai", ".bai")):
            os.remove(os.path.join(bamdir, name))
            removed += 1
    if removed:
        print(f"[coverage] removed {removed} BAM/index files")


# ============================================================================
# Labels - align contigs to the strain reference, keep the best hit (formerly
# labels.py + the subprocess-invoked _labels.py, merged into one step)
# ============================================================================
def label_from_reference(mwork: str, reference: str, tb: Toolbox) -> None:
    """minimap2 (asm10) contigs vs reference, then best-match-by-bases labels."""
    if not os.path.isfile(reference):
        raise FileNotFoundError(f"reference not found: {reference}")
    minimap2 = tb.resolve("minimap2")
    filt_path = os.path.join(mwork, "contigs_filt.fasta")
    paf = os.path.join(mwork, "contigs_vs_ref.paf")
    with open(paf, "w") as fh:
        subprocess.run([minimap2, "-x", "asm10", "-t", str(THREADS), "--secondary=no",
                        reference, filt_path], stdout=fh, check=True)

    contigs = fasta_ids(filt_path)
    qlen: dict[str, int] = {}
    matched: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q, ql, tgt, nmatch = f[0], int(f[1]), f[5], int(f[9])
            qlen[q] = ql
            matched[q][tgt] += nmatch

    out = os.path.join(mwork, "labels.csv")
    n_assigned = 0
    with open(out, "w") as fh:
        fh.write("contig,genome_id,aligned_frac\n")
        for c in contigs:
            if c in matched and matched[c]:
                best = max(matched[c], key=matched[c].get)
                frac = matched[c][best] / max(qlen.get(c, 1), 1)
                fh.write(f"{c},{best},{frac:.4f}\n")
                n_assigned += 1
            else:
                fh.write(f"{c},unbinned,0.0000\n")
    print(f"[labels] wrote {out}: {n_assigned}/{len(contigs)} contigs labelled")


# ============================================================================
# Graph construction 
# ============================================================================
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _revcomp(kmer: str) -> str:
    return "".join(_COMP[b] for b in reversed(kmer))


def build_kmer_index(k: int = 4) -> tuple[dict, int]:
    """Map each k-mer over A/C/G/T to a canonical reverse-complement column."""
    kmers = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
    canon, order = {}, {}
    for km in kmers:
        c = min(km, _revcomp(km))
        if c not in order:
            order[c] = len(order)
        canon[km] = order[c]
    return canon, len(order)


def tnf_vector(seq: str, canon: dict, ncol: int, k: int = 4) -> np.ndarray:
    """Normalized canonical k-mer frequency vector for one sequence."""
    seq = seq.upper()
    vec = np.zeros(ncol, dtype=float)
    for i in range(len(seq) - k + 1):
        idx = canon.get(seq[i:i + k])
        if idx is not None:
            vec[idx] += 1.0
    total = vec.sum()
    return vec / total if total > 0 else vec


def _clean_segment_token(tok: str) -> str:
    """'10851-;' -> '10851', '9415+,' -> '9415', 'NODE_1+' -> 'NODE_1'."""
    return tok.strip().rstrip(";,").rstrip("+-").strip()


def parse_oriented_segment_token(token: str):
    """Return (segment_id, orientation), default orientation '+' if absent."""
    token = token.strip().rstrip(";,").strip()
    if not token:
        return None
    orientation = token[-1] if token[-1] in "+-" else "+"
    segment = (token[:-1] if token[-1] in "+-" else token).strip()
    return (segment, orientation) if segment else None


def _parse_gfa_tags(fields: list[str]) -> dict:
    """Parse optional GFA ``TAG:TYPE:VALUE`` fields into JSON-safe values."""
    tags = {}
    for field_ in fields:
        parts = field_.split(":", 2)
        if len(parts) != 3:
            continue
        tag, value_type, value = parts
        try:
            if value_type in {"i", "I"}:
                value = int(value)
            elif value_type == "f":
                value = float(value)
        except ValueError:
            pass
        tags[tag] = value
    return tags


def cigar_overlap_length(cigar: str) -> int:
    """Aligned overlap length represented by a GFA CIGAR string."""
    if not cigar or cigar == "*":
        return 0
    return int(sum(int(length) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
                   if op in {"M", "=", "X"}))


def parse_gfa_rich(path: str) -> tuple[dict, list]:
    """Read oriented, attributed GFA segment (S) and link (L) records."""
    segments, links = {}, []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S" and len(fields) >= 3:
                segment_id, sequence = fields[1], fields[2]
                tags = _parse_gfa_tags(fields[3:])
                length = len(sequence) if sequence != "*" else int(tags.get("LN", 0))
                segments[segment_id] = {"length": int(length),
                                        "depth": float(tags.get("DP", 0.0)),
                                        "kmer_count": int(tags.get("KC", 0)),
                                        "tags": tags}
            elif fields[0] == "L" and len(fields) >= 6:
                overlap_cigar = fields[5]
                links.append({"from_segment": fields[1], "from_orientation": fields[2],
                             "to_segment": fields[3], "to_orientation": fields[4],
                             "overlap_cigar": overlap_cigar,
                             "overlap_length": cigar_overlap_length(overlap_cigar),
                             "tags": _parse_gfa_tags(fields[6:])})
    return segments, links


def parse_gfa_edges(path: str, contig_names: list[str]) -> set:
    """Fallback contig adjacency by direct GFA segment/contig name matching.

    Only reliable when GFA segment names already match final contig names;
    used when no contigs.paths file exists (e.g. non-SPAdes assembler).
    """
    name_set = set(contig_names)
    idx = {n: i for i, n in enumerate(contig_names)}
    seg_to_contigs: dict[str, set] = {}
    links = []
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if not t:
                continue
            if t[0] == "P" and len(t) >= 3:
                pname = t[1]
                if pname in idx:
                    for seg in t[2].split(","):
                        seg = _clean_segment_token(seg)
                        if seg:
                            seg_to_contigs.setdefault(seg, set()).add(idx[pname])
            elif t[0] == "L" and len(t) >= 5:
                links.append((_clean_segment_token(t[1]), _clean_segment_token(t[3])))
            elif t[0] == "S" and len(t) >= 2 and t[1] in name_set:
                seg_to_contigs.setdefault(t[1], set()).add(idx[t[1]])

    edges = set()
    for a, b in links:
        for i in seg_to_contigs.get(a, set()):
            for j in seg_to_contigs.get(b, set()):
                if i != j:
                    edges.add((min(i, j), max(i, j)))
    return edges


def parse_spades_paths_oriented(path: str, contig_names: list[str]) -> dict:
    """Map final SPAdes contigs to ordered, oriented internal GFA walks."""
    name_set = set(contig_names)
    mapping: dict[str, list] = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            base_name = s[:-1] if s.endswith("'") else s
            if base_name in name_set:
                if s.endswith("'"):
                    cur = None
                else:
                    cur = base_name
                    mapping.setdefault(cur, [])
                continue
            if cur is not None:
                for tok in s.replace(";", ",").split(","):
                    parsed = parse_oriented_segment_token(tok)
                    if parsed is not None:
                        mapping[cur].append(parsed)
    return mapping


def parse_flye_paths_oriented(path: str, contig_names: list[str]) -> dict:
    name_set = set(contig_names)
    mapping: dict[str, list] = {}
    with open(path) as fh:
        header = None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [f.lstrip("#").strip() for f in fields]
                continue
            if header is None or len(fields) != len(header):
                continue
            row = dict(zip(header, fields))
            name = row.get("seq_name", "")
            if name not in name_set:
                continue
            walk = []
            for token in row.get("graph_path", "").split(","):
                token = token.strip()
                if not token or token in ("*", "??"):
                    continue
                orientation = "-" if token.startswith("-") else "+"
                segment = token.lstrip("+-").strip()
                if segment:
                    walk.append((f"edge_{segment}", orientation))
            mapping[name] = walk
    return mapping


def contig_graph_from_paths(gfa_path: str, paths: dict, contig_names: list[str],
                            include_shared_segment_edges: bool = False):
    idx = {name: i for i, name in enumerate(contig_names)}
    segments, links = parse_gfa_rich(gfa_path)

    seg_to_contigs: dict[str, set] = {}
    for name, walk in paths.items():
        contig_index = idx[name]
        for segment, _ in walk:
            seg_to_contigs.setdefault(segment, set()).add(contig_index)

    edge_attributes: dict[tuple, dict] = {}

    def attributes_for(left, right):
        pair = (min(int(left), int(right)), max(int(left), int(right)))
        return pair, edge_attributes.setdefault(
            pair, {"gfa_link_count": 0, "shared_segment_count": 0,
                  "max_overlap_length": 0})

    for link in links:
        left_contigs = seg_to_contigs.get(link["from_segment"], set())
        right_contigs = seg_to_contigs.get(link["to_segment"], set())
        for left in left_contigs:
            for right in right_contigs:
                if left == right:
                    continue
                _, attributes = attributes_for(left, right)
                attributes["gfa_link_count"] += 1
                attributes["max_overlap_length"] = max(
                    attributes["max_overlap_length"], int(link["overlap_length"]))

    shared_pair_counts = defaultdict(int)
    for contig_set in seg_to_contigs.values():
        ordered = sorted(contig_set)
        for a in range(len(ordered)):
            for b in range(a + 1, len(ordered)):
                shared_pair_counts[(ordered[a], ordered[b])] += 1

    for pair, count in shared_pair_counts.items():
        if include_shared_segment_edges:
            _, attributes = attributes_for(*pair)
            attributes["shared_segment_count"] = int(count)
        else:
            edge_attributes.pop(pair, None)

    path_segments = set(seg_to_contigs)
    shared_segments = set(segments) & path_segments
    edges = set(edge_attributes)
    diag = {
        "gfa_segments": len(segments), "gfa_links": len(links),
        "filtered_contigs": len(contig_names), "contigs_in_paths": len(paths),
        "path_segments": len(path_segments), "shared_segments": len(shared_segments),
        "gfa_link_edges": sum(a["gfa_link_count"] > 0 for a in edge_attributes.values()),
        "shared_segment_edges_detected": len(shared_pair_counts),
        "shared_segment_edges_included": (
            len(shared_pair_counts) if include_shared_segment_edges else 0),
        "shared_segment_edges_excluded": (
            0 if include_shared_segment_edges else len(shared_pair_counts)),
        "projected_edges": len(edges),
    }
    rich = {"segments": segments, "links": links, "paths": paths}
    return edges, edge_attributes, rich, diag


def read_fasta(path: str) -> tuple[list[str], list[str]]:
    names, seqs, cur, buf = [], [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    names.append(cur); seqs.append("".join(buf))
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur is not None:
        names.append(cur); seqs.append("".join(buf))
    return names, seqs


def zscore(a: np.ndarray) -> np.ndarray:
    mu, sd = a.mean(0), a.std(0)
    sd[sd == 0] = 1.0
    return (a - mu) / sd


def _flip_orientation(orientation: str) -> str:
    return "-" if orientation == "+" else "+"


def build_assembly_heterograph(contig_names, contig_features, labels, rich):
    """Lossless oriented contig/segment view of the SPAdes assembly graph.

    GFA links become bidirected traversals between oriented segment states;
    contig paths retain every repeated segment, its orientation, and its
    position in the walk.
    """
    graph = HeteroData()
    graph["contig"].x = torch.from_numpy(contig_features)
    graph["contig"].y = torch.from_numpy(labels)
    graph["contig"].num_nodes = len(contig_names)

    segments, links, paths = rich["segments"], rich["links"], rich["paths"]
    all_segment_ids = set(segments)
    for link in links:
        all_segment_ids.add(link["from_segment"])
        all_segment_ids.add(link["to_segment"])
    for walk in paths.values():
        all_segment_ids.update(segment for segment, _ in walk)
    segment_ids = sorted(all_segment_ids)

    oriented_states = [(segment, orientation) for segment in segment_ids
                       for orientation in ("+", "-")]
    oriented_index = {state: index for index, state in enumerate(oriented_states)}

    segment_features = []
    for segment, orientation in oriented_states:
        attributes = segments.get(segment, {})
        segment_features.append([
            np.log1p(float(attributes.get("length", 0))),
            np.log1p(max(float(attributes.get("depth", 0.0)), 0.0)),
            np.log1p(max(float(attributes.get("kmer_count", 0)), 0.0)),
            1.0 if orientation == "+" else -1.0,
        ])
    graph["oriented_segment"].x = torch.tensor(segment_features, dtype=torch.float32)
    graph["oriented_segment"].num_nodes = len(oriented_states)

    link_source, link_target, link_attributes = [], [], []
    for link in links:
        forward = (link["from_segment"], link["from_orientation"])
        target = (link["to_segment"], link["to_orientation"])
        if forward not in oriented_index or target not in oriented_index:
            continue
        link_source.append(oriented_index[forward])
        link_target.append(oriented_index[target])
        link_attributes.append([float(link["overlap_length"]), 0.0])
        reverse_source = (link["to_segment"], _flip_orientation(link["to_orientation"]))
        reverse_target = (link["from_segment"], _flip_orientation(link["from_orientation"]))
        link_source.append(oriented_index[reverse_source])
        link_target.append(oriented_index[reverse_target])
        link_attributes.append([float(link["overlap_length"]), 1.0])

    link_relation = ("oriented_segment", "gfa_link", "oriented_segment")
    graph[link_relation].edge_index = (
        torch.tensor([link_source, link_target], dtype=torch.long)
        if link_source else torch.empty((2, 0), dtype=torch.long))
    graph[link_relation].edge_attr = (
        torch.tensor(link_attributes, dtype=torch.float32)
        if link_attributes else torch.empty((0, 2), dtype=torch.float32))
    graph[link_relation].edge_attr_names = ["overlap_length", "reverse_complement_record"]

    contig_index = {name: index for index, name in enumerate(contig_names)}
    path_contigs, path_segments, path_attributes = [], [], []
    for contig_name, walk in paths.items():
        if contig_name not in contig_index:
            continue
        denominator = max(len(walk) - 1, 1)
        for position, state in enumerate(walk):
            if state not in oriented_index:
                continue
            path_contigs.append(contig_index[contig_name])
            path_segments.append(oriented_index[state])
            path_attributes.append(
                [float(position), float(position / denominator), float(len(walk))])

    traversal = ("contig", "traverses", "oriented_segment")
    reverse_traversal = ("oriented_segment", "part_of_path", "contig")
    traversal_index = (torch.tensor([path_contigs, path_segments], dtype=torch.long)
                       if path_contigs else torch.empty((2, 0), dtype=torch.long))
    traversal_attr = (torch.tensor(path_attributes, dtype=torch.float32)
                      if path_attributes else torch.empty((0, 3), dtype=torch.float32))
    graph[traversal].edge_index = traversal_index
    graph[traversal].edge_attr = traversal_attr
    graph[traversal].edge_attr_names = [
        "path_position", "normalized_path_position", "path_length"]
    graph[reverse_traversal].edge_index = traversal_index.flip(0)
    graph[reverse_traversal].edge_attr = traversal_attr.clone()
    graph[reverse_traversal].edge_attr_names = list(graph[traversal].edge_attr_names)

    graph.contig_names = list(contig_names)
    graph.segment_names = segment_ids
    graph.oriented_segment_names = [f"{s}{o}" for s, o in oriented_states]
    graph.contig_paths = paths
    graph.gfa_links = links
    graph.segment_feature_names = ["log_length", "log_depth", "log_kmer_count", "orientation"]
    return graph


def build_graph(mwork: str, mout: str, *, n_samples: int) -> Data:
    """Build viral_graph.pt (+ assembly_heterograph.pt).

    Nodes = contigs. Features = canonical TNF (z-scored) | log coverage
    (z-scored). Edges = real SPAdes contig links (GFA projected via
    contigs.paths) + optional kNN. Labels = source strain; unbinned -> -1.
    """
    os.makedirs(mout, exist_ok=True)
    names, seqs = read_fasta(os.path.join(mwork, "contigs_filt.fasta"))
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    print(f"[graph] {n} contigs")
    if n == 0:
        raise ValueError("contigs_filt.fasta has no contigs; cannot build graph")
    if len(idx) != n:
        raise ValueError("contigs_filt.fasta contains duplicate contig identifiers")

    canon, ncol = build_kmer_index(TNF_K)
    tnf = np.vstack([tnf_vector(s, canon, ncol, TNF_K) for s in seqs])

    cov = np.zeros((n, n_samples))
    covered = set()
    with open(os.path.join(mwork, "coverage.csv")) as fh:
        header = next(fh).rstrip("\n").split(",")
        expected = ["contig"] + [f"sample{i}" for i in range(n_samples)]
        if header != expected:
            raise ValueError(f"coverage.csv header mismatch: expected {expected}, got {header}")
        for line in fh:
            p = line.rstrip("\n").split(",")
            if p[0] in idx:
                if p[0] in covered:
                    raise ValueError(f"coverage.csv has a duplicate row for {p[0]!r}")
                if len(p) != n_samples + 1:
                    raise ValueError(f"coverage row for {p[0]!r} has {len(p)-1} values; "
                                     f"expected {n_samples}")
                cov[idx[p[0]]] = [float(x) for x in p[1:]]
                covered.add(p[0])
    missing = [name for name in names if name not in covered]
    if missing:
        raise ValueError(f"coverage.csv is missing contigs, e.g. {missing[:5]}")
    x = np.hstack([zscore(tnf), zscore(np.log1p(cov))]).astype(np.float32)

    lab = {}
    with open(os.path.join(mwork, "labels.csv")) as fh:
        next(fh)
        for line in fh:
            c, g, _ = line.rstrip("\n").split(",")
            lab[c] = g
    classes = sorted({g for g in lab.values() if g != "unbinned"})
    cls_id = {g: i for i, g in enumerate(classes)}
    y = np.array([cls_id.get(lab.get(c, "unbinned"), -1) for c in names], dtype=np.int64)
    print(f"[graph] {len(classes)} genome classes; "
          f"{(y >= 0).sum()} labelled / {(y < 0).sum()} unbinned")

    edges: set = set()
    edge_attributes: dict = {}
    rich_graph = None
    method = "none"
    gfa = os.path.join(mwork, "assembly_graph.gfa")
    # Each assembler states contig -> graph-walk in its own file; the
    # projection downstream of the parser is identical for both.
    path_table, path_parser = {
        "metaspades": ("contigs.paths", parse_spades_paths_oriented),
        "metaflye": ("assembly_info.txt", parse_flye_paths_oriented),
    }[ASSEMBLER]
    paths = os.path.join(mwork, path_table)
    if USE_GFA_EDGES and os.path.exists(gfa):
        if os.path.exists(paths):
            ge, edge_attributes, rich_graph, diag = contig_graph_from_paths(
                gfa, path_parser(paths, names), names,
                include_shared_segment_edges=INCLUDE_SHARED_SEGMENT_EDGES)
            method = "oriented-path-projection"
            print("[graph] --- assembly-graph diagnostics ---")
            print(f"[graph] raw GFA segments:   {diag['gfa_segments']}")
            print(f"[graph] raw GFA links:      {diag['gfa_links']}")
            print(f"[graph] filtered contigs:   {n}")
            print(f"[graph] contigs in paths:   {diag['contigs_in_paths']}")
            print(f"[graph] projected edges:    {diag['projected_edges']}")
        else:
            ge = parse_gfa_edges(gfa, names)
            method = "direct-gfa"
            edge_attributes = {e: {"gfa_link_count": 1, "shared_segment_count": 0,
                                   "max_overlap_length": 0} for e in ge}
            print(f"[graph] WARNING: no {path_table}; falling back to direct GFA matching")
        edges |= ge
    if USE_KNN_EDGES and n > KNN_K:
        A = kneighbors_graph(x, KNN_K, mode="connectivity").tocoo()
        for i, j in zip(A.row, A.col):
            if i != j:
                pair = (min(int(i), int(j)), max(int(i), int(j)))
                edges.add(pair)
                attributes = edge_attributes.setdefault(
                    pair, {"gfa_link_count": 0, "shared_segment_count": 0,
                          "max_overlap_length": 0})
                attributes["knn_edge"] = 1
        print(f"[graph] edges after kNN: {len(edges)}")

    if edges:
        ordered_edges = sorted(edges)
        undirected_edge_index = np.array(ordered_edges).T
        ei = np.hstack([undirected_edge_index, undirected_edge_index[::-1]])
        undirected_edge_attr = np.asarray([
            [edge_attributes.get(e, {}).get("gfa_link_count", 0),
             edge_attributes.get(e, {}).get("shared_segment_count", 0),
             edge_attributes.get(e, {}).get("max_overlap_length", 0),
             edge_attributes.get(e, {}).get("knn_edge", 0)]
            for e in ordered_edges], dtype=np.float32)
        edge_attr = np.vstack([undirected_edge_attr, undirected_edge_attr])
        edge_type = np.asarray([
            (1 if a[0] > 0 else 0) | (2 if a[1] > 0 else 0) | (4 if a[3] > 0 else 0)
            for a in edge_attr], dtype=np.int64)
    else:
        ei = np.zeros((2, 0), dtype=np.int64)
        edge_attr = np.zeros((0, 4), dtype=np.float32)
        edge_type = np.zeros((0,), dtype=np.int64)

    degree = np.zeros(n, dtype=int)
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    connected = int((degree > 0).sum())
    print(f"[graph] edge method:        {method}")
    print(f"[graph] contig-level edges: {len(edges)}")
    print(f"[graph] connected nodes:    {connected}  isolated: {n - connected}")

    data = Data(x=torch.from_numpy(x), edge_index=torch.from_numpy(ei.astype(np.int64)),
               edge_attr=torch.from_numpy(edge_attr), edge_type=torch.from_numpy(edge_type),
               y=torch.from_numpy(y))
    data.contig_names = names
    data.num_classes = len(classes)
    data.edge_attr_names = ["gfa_link_count", "shared_segment_count",
                            "max_overlap_length", "knn_edge"]
    data.edge_type_bitmask = {"gfa_link": 1, "shared_segment": 2, "feature_knn": 4}

    out = os.path.join(mout, "viral_graph.pt")
    torch.save(data, out)
    if rich_graph is not None:
        rich_out = os.path.join(mout, "assembly_heterograph.pt")
        heterograph = build_assembly_heterograph(names, x, y, rich_graph)
        heterograph["contig"].num_classes = len(classes)
        torch.save(heterograph, rich_out)
        data.assembly_heterograph_path = rich_out
        torch.save(data, out)  # re-save with the companion-graph path attached
        print(f"[graph] saved rich assembly graph: {rich_out}")
    print(f"[graph] saved {out}\n{data}")
    return data


# ============================================================================
# Contig metadata table
# ============================================================================
def write_contig_metadata(mwork: str, mout: str) -> str:
    """node_index, contig_id, source_sample, original_contig_id, length, genome_label.

    Node order matches build_graph() (reads contigs_filt.fasta in file order),
    so node_index here is the graph node index.
    """
    labels = {}
    labels_csv = os.path.join(mwork, "labels.csv")
    if _ok(labels_csv):
        with open(labels_csv) as fh:
            next(fh, None)
            for line in fh:
                p = line.rstrip("\n").split(",")
                if len(p) >= 2:
                    labels[p[0]] = p[1]

    combined = os.path.join(mwork, "contigs_filt.fasta")
    out = os.path.join(mout, "contig_metadata.tsv")
    node = 0
    with open(combined) as fin, open(out, "w") as fout:
        fout.write("node_index\tcontig_id\tsource_sample\toriginal_contig_id\t"
                   "length\tgenome_label\n")
        cur_id, cur_len = None, 0

        def flush(cid, clen):
            fout.write(f"{node}\t{cid}\tall\t{cid}\t{clen}\t{labels.get(cid, '')}\n")

        for line in fin:
            if line.startswith(">"):
                if cur_id is not None:
                    flush(cur_id, cur_len)
                    node += 1
                cur_id = line[1:].split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id is not None:
            flush(cur_id, cur_len)
    print(f"[graph] wrote node metadata: {out} ({node + 1 if cur_id else 0} rows)")
    return out


# ============================================================================
# Viral-specific features: ORF architecture + protein clusters (formerly
# viral_features.py).
# ============================================================================
MIN_SEQ_ID = 0.3      # mmseqs clustering stringency: "same protein family"
MIN_COVERAGE = 0.5
MIN_ORF_AA = 30


def _read_fasta_desc(path: str):
    """Yield (name, sequence); name is the header up to the first whitespace."""
    name, chunks = None, []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def _call_orfs(contigs_fasta: str, proteins_out: str) -> dict:
    """Call genes per contig with pyrodigal-gv; write proteins, return stats."""
    import pyrodigal_gv
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    stats = {}
    n_proteins = 0
    with open(proteins_out, "w", encoding="utf-8", newline="\n") as out:
        for contig_id, sequence in _read_fasta_desc(contigs_fasta):
            length = len(sequence)
            genes = finder.find_genes(sequence.encode())
            coding_bp, strands, lengths, gc_values, tables = 0, [], [], [], []
            index = 0
            for gene in genes:
                protein = gene.translate()
                if len(protein) < MIN_ORF_AA:
                    continue
                index += 1
                span = abs(gene.end - gene.begin) + 1
                coding_bp += span
                lengths.append(span)
                strands.append(gene.strand)
                gc_values.append(gene.gc_cont)
                tables.append(gene.translation_table)
                out.write(f">{contig_id}__orf{index}\n")
                for start in range(0, len(protein), 60):
                    out.write(protein[start:start + 60] + "\n")
                n_proteins += 1
            switches = sum(1 for a, b in zip(strands, strands[1:]) if a != b)
            stats[contig_id] = {
                "contig_id": contig_id, "length": length, "n_orfs": index,
                "coding_density": round(coding_bp / length, 4) if length else 0.0,
                "mean_orf_len": round(sum(lengths) / index, 1) if index else 0.0,
                "orfs_per_kb": round(1000.0 * index / length, 4) if length else 0.0,
                "strand_switch_rate": round(switches / (index - 1), 4) if index > 1 else 0.0,
                "gc_mean": round(sum(gc_values) / index, 4) if index else 0.0,
                "translation_table": (Counter(tables).most_common(1)[0][0] if tables else 11),
            }
    print(f"[viralfeat] called {n_proteins:,} proteins over {len(stats):,} contigs")
    return stats


def _cluster_proteins(proteins: str, workdir: str, tb: Toolbox) -> dict:
    """mmseqs2 easy-cluster; return protein_id -> cluster representative."""
    mmseqs = tb.resolve("mmseqs")
    prefix = os.path.join(workdir, "clu")
    tmp = os.path.join(workdir, "mmseqs_tmp")
    os.makedirs(workdir, exist_ok=True)
    run([mmseqs, "easy-cluster", proteins, prefix, tmp,
         "--min-seq-id", str(MIN_SEQ_ID), "-c", str(MIN_COVERAGE), "--cov-mode", "1",
         "--threads", str(THREADS)], log_prefix="[viralfeat] ")
    shutil.rmtree(tmp, ignore_errors=True)
    table = f"{prefix}_cluster.tsv"
    if not os.path.isfile(table):
        raise RuntimeError(f"mmseqs produced no cluster table at {table}")
    membership = {}
    with open(table, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                membership[parts[1]] = parts[0]
    print(f"[viralfeat] {len(set(membership.values())):,} protein clusters "
          f"over {len(membership):,} proteins")
    return membership


def extract_viral_features(mwork: str, mout: str, tb: Toolbox, *, force: bool = False) -> dict:
    """Build viral_features.tsv + protein_clusters.tsv; returns a report dict."""
    contigs = os.path.join(mwork, "contigs_filt.fasta")
    features_tsv = os.path.join(mout, "viral_features.tsv")
    clusters_tsv = os.path.join(mout, "protein_clusters.tsv")
    if not force and os.path.isfile(features_tsv) and os.path.isfile(clusters_tsv):
        print(f"[viralfeat] reusing existing {features_tsv} (--force to rebuild)")
        return {"status": "reused"}

    workdir = os.path.join(mwork, "viral_features")
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(mout, exist_ok=True)
    proteins = os.path.join(workdir, "proteins.faa")
    stats = _call_orfs(contigs, proteins)

    columns = ["contig_id", "length", "n_orfs", "coding_density", "mean_orf_len",
              "orfs_per_kb", "strand_switch_rate", "gc_mean", "translation_table"]
    with open(features_tsv, "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(columns) + "\n")
        for contig_id, row in stats.items():
            out.write("\t".join(str(row[c]) for c in columns) + "\n")
    print(f"[viralfeat] wrote {features_tsv}")

    if os.path.getsize(proteins) == 0:
        print("[viralfeat] no proteins called; writing an empty cluster table")
        with open(clusters_tsv, "w", encoding="utf-8", newline="\n") as out:
            out.write("contig_id\tcluster_id\tn_proteins\tcluster_size\tcluster_contigs\n")
        return {"status": "ok", "contigs": len(stats), "proteins": 0, "clusters": 0}

    membership = _cluster_proteins(proteins, workdir, tb)
    per_pair = Counter()
    cluster_contigs = defaultdict(set)
    cluster_size = Counter()
    for protein_id, representative in membership.items():
        contig_id = protein_id.rsplit("__orf", 1)[0]
        per_pair[(contig_id, representative)] += 1
        cluster_contigs[representative].add(contig_id)
        cluster_size[representative] += 1

    with open(clusters_tsv, "w", encoding="utf-8", newline="\n") as out:
        out.write("contig_id\tcluster_id\tn_proteins\tcluster_size\tcluster_contigs\n")
        for (contig_id, representative), count in sorted(per_pair.items()):
            out.write(f"{contig_id}\t{representative}\t{count}\t"
                      f"{cluster_size[representative]}\t{len(cluster_contigs[representative])}\n")
    print(f"[viralfeat] wrote {clusters_tsv}")

    singleton = sum(1 for c in cluster_contigs.values() if len(c) == 1)
    return {"status": "ok", "contigs": len(stats), "proteins": len(membership),
            "clusters": len(cluster_contigs), "single_contig_clusters": singleton}


# ============================================================================
# Work-directory cleanup 
# ============================================================================

FINAL_WORK_FILES = {
    "contigs_filt.fasta", "coverage.csv", "labels.csv",
    "assembly_graph.gfa", "contigs.paths", "assembly_info.txt",
}


def cleanup_intermediate_work(mwork: str, work_root: str, output_dir: str) -> list[str]:
    run_dir = os.path.realpath(mwork)
    root = os.path.realpath(work_root)
    output = os.path.realpath(output_dir)
    under_root = os.path.commonpath([run_dir, root]) == root
    output_under_run = os.path.commonpath([run_dir, output]) == run_dir
    if not under_root or run_dir == root:
        raise RuntimeError(f"refusing cleanup outside a dedicated run directory: {run_dir}")
    if output_under_run:
        raise RuntimeError(f"refusing cleanup: output directory is inside scratch: {output}")

    removed = []
    for entry in os.scandir(run_dir):
        if entry.is_file(follow_symlinks=False) and entry.name in FINAL_WORK_FILES:
            continue
        if entry.is_dir(follow_symlinks=False):
            shutil.rmtree(entry.path)
        else:
            os.remove(entry.path)
        removed.append(entry.name)
    return removed


# ============================================================================
# Orchestration (formerly standardize.py's main(), trimmed to this one recipe)
# ============================================================================
def check_tools(tb: Toolbox) -> None:
    """Resolve every binary this fixed recipe needs before starting."""
    tb.resolve("flye" if ASSEMBLER == "metaflye" else ASSEMBLER)
    tb.resolve("seqkit")
    tb.resolve("minimap2")
    tb.resolve("samtools")
    tb.resolve("mmseqs")
    import pyrodigal_gv  


def main() -> int:
    # Not read from __doc__: the module docstring is not guaranteed to survive
    # comment stripping, and argparse must not be what breaks the pipeline.
    parser = argparse.ArgumentParser(
        description="Phase 5 - assemble, filter, coverage, labels, graph.")
    parser.add_argument("--force", action="store_true",
                        help="Recompute every stage even if cached.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the plan and tool check, then exit.")
    args = parser.parse_args()

    tb = Toolbox(conda_bin=CONDA_BIN)

    print("=" * 68)
    print(f"multi-biome benchmark  ::  {RUN_NAME}")
    print("=" * 68)
    print(f"tier          : {TIER}")
    print(f"simulator     : {SIM} ({READ_TYPE} reads)")
    print(f"assembler     : {ASSEMBLER} (minimap2 -x {MINIMAP_PRESET})")
    print(f"reads         : {READS_DIR}")
    print(f"reference     : {REFERENCE_FASTA}")
    print(f"work dir      : {MWORK}")
    print(f"output dir    : {MOUT}")

    try:
        check_tools(tb)
        print("tools         : all required binaries resolved OK")
    except (ToolError, ImportError) as exc:
        print(f"tools         : MISSING -> {exc}")
        print("\nRun 0_check_environment.sh for install commands.")
        return 3

    det = detect_reads(READS_DIR)
    if det.n_samples <= 0:
        print("error: no read samples found under", READS_DIR, file=sys.stderr)
        return 2
    print(f"samples       : {det.n_samples} "
          f"({'single/mixed-end' if det.single_end else 'paired-end'})")
    for s in det.samples:
        mate = os.path.basename(s.r2) if s.r2 else "(none)"
        print(f"    {s.name}: R1={os.path.basename(s.r1)}  R2={mate}")

    os.makedirs(MWORK, exist_ok=True)
    os.makedirs(MOUT, exist_ok=True)

    # config.json is written as a provenance record; nothing in THIS script
    # reads it back, but downstream consumers (analysis notebooks) do -- in
    # particular the feature-layout triple below, which tells them where the
    # TNF block ends and the coverage block begins inside `x`.
    tnf_width = build_kmer_index(TNF_K)[1]   # canonical k-mer count (136 for k=4)
    observed_read_length = measure_read_length(det)
    print(f"read length   : {observed_read_length:,.1f} bp (measured)")
    config = {
        "_generated_by": "Dataset_Processing/Simulated_Data/5_assemble_and_build_graph.py",
        "tier": TIER, "dataset_source": str(READS_DIR),
        "work_dir": str(MWORK), "dataset_dir": str(MOUT), "contig_mode": RUN_NAME,
        "assembly_mode": "coassembly", "n_samples": det.n_samples,
        "min_contig_len": MIN_CONTIG_LEN, "tnf_k": TNF_K, "knn_k": KNN_K,
        # Feature layout of graph.x = [ canonical TNF | log coverage ].
        "tnf_feature_width": tnf_width,
        "coverage_feature_start": tnf_width,
        "coverage_feature_end": tnf_width + det.n_samples,
        "use_gfa_edges": USE_GFA_EDGES, "use_knn_edges": USE_KNN_EDGES,
        "include_shared_segment_edges": INCLUDE_SHARED_SEGMENT_EDGES,
        "cpus": THREADS, "memory_gb": MEMORY_GB, "dataset_type": "real",
        "read_type": READ_TYPE, "minimap_preset": MINIMAP_PRESET,
        "mean_read_length": round(observed_read_length, 2),
        "simulator": SIM, "assembler": ASSEMBLER,
        "only_assembler": ASSEMBLER == "metaspades", "coverage_method": "minimap2",
        "multimap_policy": "primary_only", "min_mapq": MIN_MAPQ,
        "viral_tools": [], "reference": str(REFERENCE_FASTA),
        "conda_bin": CONDA_BIN,
        "combined_fasta": os.path.join(MWORK, "contigs_filt.fasta"),
        "graph_path": os.path.join(MOUT, "viral_graph.pt"),
        "contig_metadata": os.path.join(MOUT, "contig_metadata.tsv"),
        "samples": [{"name": s.name, "stem": s.stem, "r1": s.r1, "r2": s.r2}
                    for s in det.samples],
    }
    config_path = os.path.join(MOUT, "config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"config        : {config_path}")

    if args.dry_run:
        print("\n[dry-run] plan looks good; exiting before execution.")
        return 0

    pipeline_start = time.monotonic()

    def stage(label):
        elapsed = time.strftime("%Hh%Mm%Ss", time.gmtime(time.monotonic() - pipeline_start))
        now = datetime.datetime.now().strftime("%H:%M:%S")
        print(f"\n=== {label}  ({now}, +{elapsed} total) ===", flush=True)

    stage(f"assemble ({ASSEMBLER})")
    assemble(det, str(MWORK), tb, force=args.force)

    stage(f"filter (contigs >= {MIN_CONTIG_LEN} bp)")
    length_filter(str(MWORK), tb, min_len=MIN_CONTIG_LEN)

    stage(f"coverage (map {det.n_samples} samples)")
    map_samples(det, str(MWORK), tb, force=args.force)
    compute_coverage_matrix(str(MWORK), det.n_samples)
    remove_bams(str(MWORK))

    stage("labels (align contigs to strain reference)")
    label_from_reference(str(MWORK), str(REFERENCE_FASTA), tb)

    stage("graph (build viral_graph.pt)")
    build_graph(str(MWORK), str(MOUT), n_samples=det.n_samples)
    write_contig_metadata(str(MWORK), str(MOUT))

    stage("viral features (ORF architecture + protein clusters)")
    viralfeat_report = extract_viral_features(str(MWORK), str(MOUT), tb, force=args.force)

    manifest = {
        "run": RUN_NAME, "dataset": f"multibiome{TIER_SUFFIX}", "tier": TIER,
        "dataset_type": "real",
        "assembly_mode": "coassembly", "read_type": READ_TYPE, "read_qc": "none",
        "simulator": SIM, "assembler": ASSEMBLER, "viral_tools": [],
        "coverage_method": "minimap2",
        "config": config_path, "graph": os.path.join(MOUT, "viral_graph.pt"),
        "contig_metadata": os.path.join(MOUT, "contig_metadata.tsv"),
        "n_samples": det.n_samples, "labelled": True,
        "viral_features": viralfeat_report,
    }
    with open(os.path.join(MOUT, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    graph_path = os.path.join(MOUT, "viral_graph.pt")
    if os.path.isfile(graph_path):
        removed = cleanup_intermediate_work(str(MWORK), str(WORK_DIR / "standardize_work"), str(MOUT))
        print(f"\n[cleanup] removed {len(removed)} intermediate entries from {MWORK}")
        if removed:
            print("[cleanup] removed: " + ", ".join(sorted(removed)))

    print("\n" + "=" * 68)
    print("DONE")
    print(f"  graph    -> {graph_path}")
    print(f"  manifest -> {os.path.join(MOUT, 'manifest.json')}")
    print("=" * 68)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
