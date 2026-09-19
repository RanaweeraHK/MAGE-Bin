#!/usr/bin/env python
"""Phase 1 - download one mock dataset's reads and reference set.

    VB_DATASET=phage_mock_illumina python 1_get_inputs.py
    VB_DATASET=dsmz_dsrna         python 1_get_inputs.py --dry-run

Three steps:
  1. Download the libraries in the registry's `accessions` from its `provider`
     (ena = FASTQ off the ENA mirror, md5-verified; dataverse = an INRAE
     deposit, addressed by DOI).
  2. Put them in Data/work/mock_<dataset>/reads/ under names phase 2 pairs on,
     converting FASTA reads to FASTQ on the way.
  3. Build the reference set - reference_genomes.fasta + reference_map.tsv -
     which is the ground truth phase 3 labels contigs against.

Resumable: anything already in place and intact is skipped.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_HOME = HERE.parent.parent
REGISTRY = HERE / "datasets.tsv"

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DATAVERSE_HOST = "https://entrepot.recherche.data.gouv.fr"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

NCBI_DELAY_S = 0.4          # NCBI allows 3 requests/second without an API key
PLACEHOLDER_QUAL = "I"      # quality for FASTA reads; nothing downstream reads it


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_registry(path: Path) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(
            (line for line in fh if line.strip() and not line.startswith("#")),
            delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise SystemExit(f"error: {path} has no header row")
        for parts in reader:
            row = dict(zip(header, parts))
            rows[row["dataset"]] = row
    return rows


def fetch(url: str, *, retries: int = 4, backoff: int = 5) -> bytes:
    last = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                return resp.read()
        except Exception as exc:                      # noqa: BLE001 - network
            last = exc
            log(f"  retry {attempt}/{retries} {url.split('/')[-1][:60]}: {exc}")
            time.sleep(backoff * attempt)
    raise SystemExit(f"error: could not fetch {url}: {last}")


def download_to(url: str, dest: Path, *, expect_md5: str = "",
                expect_bytes: int = 0, retries: int = 4) -> None:
    """Stream ``url`` to ``dest``, verified. An intact file already there is kept."""
    if dest.exists() and _intact(dest, expect_md5, expect_bytes):
        log(f"  have {dest.name}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        tmp = dest.with_suffix(dest.suffix + ".part")
        try:
            with urllib.request.urlopen(url, timeout=600) as resp, open(tmp, "wb") as out:
                shutil.copyfileobj(resp, out, length=1 << 20)
        except Exception as exc:                      # noqa: BLE001 - network
            log(f"  retry {attempt}/{retries} {dest.name}: {exc}")
            tmp.unlink(missing_ok=True)
            time.sleep(5 * attempt)
            continue
        if _intact(tmp, expect_md5, expect_bytes):
            tmp.replace(dest)
            log(f"  ok   {dest.name}")
            return
        log(f"  retry {attempt}/{retries} {dest.name}: checksum/size mismatch")
        tmp.unlink(missing_ok=True)
        time.sleep(5 * attempt)
    raise SystemExit(f"error: could not download {dest.name} intact from {url}")


def _intact(path: Path, expect_md5: str, expect_bytes: int) -> bool:
    if expect_bytes and path.stat().st_size != expect_bytes:
        return False
    if expect_md5:
        h = hashlib.md5()                             # noqa: S324 - ENA publishes md5
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != expect_md5:
            return False
    return bool(expect_md5 or expect_bytes) or path.stat().st_size > 0


# ===========================================================================
# Providers - where a library's bytes come from
# ===========================================================================
def ena_manifest(study: str, runs: list[str]) -> dict[str, dict]:
    """run -> {urls, md5s, bytes} from ENA's filereport for the whole study."""
    query = urllib.parse.urlencode({
        "accession": study, "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes,instrument_platform",
        "format": "tsv"})
    text = fetch(f"{ENA_FILEREPORT}?{query}").decode()
    wanted = set(runs)
    out: dict[str, dict] = {}
    for line in text.splitlines()[1:]:
        run, ftp, md5, size, platform = (line.split("\t") + [""] * 5)[:5]
        if run not in wanted:
            continue
        out[run] = {"urls": [u for u in ftp.split(";") if u],
                    "md5s": md5.split(";"), "bytes": [int(b) for b in size.split(";") if b],
                    "platform": platform}
    missing = wanted - set(out)
    if missing:
        raise SystemExit(f"error: {study} has no runs {sorted(missing)}")
    return out


def dataverse_manifest(doi: str, names: list[str]) -> dict[str, dict]:
    """file name -> {id, bytes} from the deposit's own file listing."""
    query = urllib.parse.urlencode({"persistentId": doi})
    payload = json.loads(fetch(f"{DATAVERSE_HOST}/api/datasets/:persistentId/?{query}"))
    listing = {f["dataFile"]["filename"]: f["dataFile"]
               for f in payload["data"]["latestVersion"]["files"]}
    out = {}
    for name in names:
        if name not in listing:
            raise SystemExit(
                f"error: {doi} has no file {name!r}\n"
                f"       it holds: {', '.join(sorted(listing)[:8])}...")
        out[name] = {"id": listing[name]["id"], "bytes": listing[name]["filesize"]}
    return out


# ===========================================================================
# Reads - download, then normalise into the names phase 2 pairs on
# ===========================================================================
def _open_maybe_gzip(path: Path):
    with open(path, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")


def fasta_to_fastq_gz(src: Path, dst: Path) -> int:
    """Rewrite a FASTA read file as gzipped FASTQ with a placeholder quality.

    The deposited reads are already quality-trimmed and carry no qualities.
    Phase 2 detects reads by FASTQ extension and metaSPAdes runs
    --only-assembler, so nothing downstream reads the placeholder.
    """
    n = 0
    with _open_maybe_gzip(src) as fin, gzip.open(dst, "wt", compresslevel=1) as fout:
        name, seq = None, []
        def flush():
            nonlocal n
            if name is None:
                return
            s = "".join(seq)
            if not s:
                return
            # The INRAE files name every read ">No_name", so number them.
            fout.write(f"@{name}_{n}\n{s}\n+\n{PLACEHOLDER_QUAL * len(s)}\n")
            n += 1
        for line in fin:
            if line.startswith(">"):
                flush()
                name, seq = (line[1:].split() or ["read"])[0], []
            else:
                seq.append(line.strip())
        flush()
    return n


def get_reads(row: dict, reads_dir: Path, study_dir: Path, *, dry_run: bool,
              discard_downloads: bool = False) -> None:
    libraries = [a for a in row["accessions"].split(",") if a]
    provider = row["provider"]
    cache = study_dir / "downloads"
    reads_dir.mkdir(parents=True, exist_ok=True)

    if provider == "ena":
        manifest = {} if dry_run else ena_manifest(row["source_id"], libraries)
        for run in libraries:
            if dry_run:
                log(f"  would fetch {run} from ENA {row['source_id']}")
                continue
            entry = manifest[run]
            for i, url in enumerate(entry["urls"]):
                # ENA's <run>_1/_2.fastq.gz naming is what phase 2 pairs on.
                dest = reads_dir / Path(url).name
                download_to(f"https://{url}", dest,
                            expect_md5=entry["md5s"][i] if i < len(entry["md5s"]) else "",
                            expect_bytes=entry["bytes"][i] if i < len(entry["bytes"]) else 0)
    elif provider == "dataverse":
        manifest = {} if dry_run else dataverse_manifest(row["source_id"], libraries)
        for name in libraries:
            stem = re.sub(r"\.(fa|fasta|fq|fastq)(\.gz)?$", "", name)
            final = reads_dir / f"{stem}.fastq.gz"
            if dry_run:
                log(f"  would fetch {name} from {row['source_id']} -> {final.name}")
                continue
            if final.exists() and final.stat().st_size > 0:
                log(f"  have {final.name}")
                continue
            entry = manifest[name]
            raw = cache / name
            download_to(f"{DATAVERSE_HOST}/api/access/datafile/{entry['id']}?format=original",
                        raw, expect_bytes=entry["bytes"])
            if row["reads_format"] == "fasta":
                tmp = final.with_suffix(".part.gz")
                count = fasta_to_fastq_gz(raw, tmp)
                tmp.replace(final)
                log(f"  conv {final.name} ({count:,} reads, placeholder quality)")
            else:
                shutil.copy(raw, final)
            if discard_downloads:
                raw.unlink()
                log(f"  rm   {raw.name} (--discard-downloads)")
    else:
        raise SystemExit(f"error: unknown provider {provider!r} (use ena|dataverse)")


# ===========================================================================
# Reference sets - the ground truth
# ===========================================================================
_ACCESSION = re.compile(r"^[A-Z]{1,2}_?\d{5,8}(\.\d+)?$")
_SEGMENT = re.compile(
    r"^(RNA\d+[A-Za-z]*|DNA-?[A-Za-z0-9]+|segment[_-]?\w+|S\d|M\d|L\d|clone_\d+)$",
    re.IGNORECASE)


def split_reference_name(header_id: str) -> tuple[str, str]:
    """Split a reference id into (agent, molecule).

    INRAE ids stack virus, segment and accession:
    `Alfalfa_mosaic_virus_RNA1_MZ405653`, `Banana_bunchy_top_virus_DNA-C`.
    The agent is the virus a binner should recover; the molecule is the
    segment a contig came from. Phase 3 writes both.
    """
    tokens = header_id.split("_")
    while tokens and _ACCESSION.match(tokens[-1]):
        tokens.pop()
    agent_tokens = list(tokens)
    if len(agent_tokens) > 1 and _SEGMENT.match(agent_tokens[-1]):
        agent_tokens.pop()
    agent = "_".join(agent_tokens) or header_id
    return agent, header_id


def write_reference_map(fasta: Path, out: Path) -> int:
    rows = []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            ref = line[1:].split()[0]
            agent, molecule = split_reference_name(ref)
            rows.append((ref, molecule, agent))
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reference_id", "molecule_label", "agent_label"])
        w.writerows(rows)
    return len(rows)


def reference_set_dsmz(refdir: Path, *, dry_run: bool) -> None:
    """The 115 molecules of the DSMZ/INRAE community, from the P60 deposit."""
    fasta = refdir / "reference_genomes.fasta"
    if dry_run:
        log("  would fetch Reference genomes.fas from doi:10.57745/T4UYPC")
        return
    if not (fasta.exists() and fasta.stat().st_size > 0):
        manifest = dataverse_manifest("doi:10.57745/T4UYPC", ["Reference genomes.fas"])
        entry = manifest["Reference genomes.fas"]
        download_to(
            f"{DATAVERSE_HOST}/api/access/datafile/{entry['id']}?format=original",
            fasta, expect_bytes=entry["bytes"])
    n = write_reference_map(fasta, refdir / "reference_map.tsv")
    agents = len({line.split("\t")[2] for line in
                  (refdir / "reference_map.tsv").read_text().splitlines()[1:]})
    log(f"  reference set: {n} molecules over {agents} agents")


# The 15 phages of PRJEB56639, from Table S1 of Microb Genom 2024
# 10.1099/mgen.0.001198. Names like "PHAGE1" and "J1" are not unique in NCBI,
# so a candidate is accepted only if its length and GC match this table.
PHAGE_MOCK_15 = [
    # label,      NCBI search term,                            length,  GC%
    ("SLUR29",    "vB_Eco_SLUR29",                              48593,  44.75),
    ("PARMAL1",   "PARMAL1",                                    44565,  57.58),
    ("J2",        "Escherichia phage vB_Eco_mar002J2",          50343,  44.38),
    ("HP1",       "Haemophilus phage HP1",                      32355,  40.01),
    ("SM033",     "vB_VpaM_sm033",                             320253,  43.23),
    ("J3",        "Escherichia phage vB_Eco_mar003J3",         115471,  39.78),
    ("SWAN",      "vB_Eco_swan01",                              50865,  44.74),
    ("KUW1",      "ParKuw1",                                    44509,  60.76),
    ("VP1",       "DSS3_PM1",                                   70044,  47.34),
    ("J1",        "Escherichia phage vB_Eco_mar001J1",          50343,  44.40),
    ("PHAGE1",    "vB_Eco_mar005P1",                           167773,  37.72),
    # Table S1 calls this one sm032; it is deposited as vB_VpaS_sm030, from the
    # same submission as sm033, at exactly the published length.
    ("SM032",     "vB_VpaS_sm030",                              79660,  45.67),
    ("phix174",   "NC_001422",                                   5386,  44.76),
    ("SRSM4",     "Synechococcus phage S-RSM4",                194454,  41.12),
    ("CDMH1",     "NC_024144",                                  54279,  28.45),
]
# The published length is the authors' own assembly, which can differ slightly
# from the deposited genome (SLUR29: 48,593 vs 48,466), so the tolerance is
# proportional. GC then separates candidates of the same size.
LENGTH_TOLERANCE_FRAC = 0.01
LENGTH_TOLERANCE_MIN_BP = 30
GC_TOLERANCE_PCT = 1.0


def _eutils(path: str, params: dict) -> bytes:
    time.sleep(NCBI_DELAY_S)
    return fetch(f"{NCBI_EUTILS}/{path}?{urllib.parse.urlencode(params)}")


def _norm(text: str) -> str:
    """Lowercase alphanumerics only, so vB_Eco_swan01 matches "vB_Eco_swan01"."""
    return re.sub(r"[^a-z0-9]", "", text.lower())


def _gc_percent(seq: str) -> float:
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    return 100.0 * (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0


def resolve_phage_genome(term: str, want_len: int, want_gc: float) -> tuple[str, str, str]:
    """Return (accession, header, sequence) for the candidate matching length+GC."""
    ids = re.findall(rb"<Id>(\d+)</Id>", _eutils("esearch.fcgi", {
        "db": "nuccore", "term": term, "retmax": "20"}))
    if not ids:
        raise LookupError(f"no nuccore hit for {term!r}")

    # Filter on length first: esummary is a few kB, while fetching twenty full
    # genomes to measure them is tens of MB.
    summary = json.loads(_eutils("esummary.fcgi", {
        "db": "nuccore", "id": ",".join(i.decode() for i in ids),
        "retmode": "json"}).decode())["result"]
    tolerance = max(LENGTH_TOLERANCE_MIN_BP,
                    int(want_len * LENGTH_TOLERANCE_FRAC))
    # Relatives in this mock differ by less than the tolerance - vB_Eco_Sip is
    # 50,809 bp against swan01's 50,865 - so candidates are RANKED, not taken in
    # the order esearch returns them: exact length first, and a title carrying
    # the phage's own name ahead of one that does not.
    hint = _norm(term.split("_")[-1])
    seen_lengths, ranked = [], []
    for uid in summary.get("uids", []):
        length = int(summary[uid].get("slen", 0))
        seen_lengths.append(length)
        if abs(length - want_len) > tolerance:
            continue
        named = hint in _norm(summary[uid].get("title", ""))
        ranked.append((0 if named else 1, abs(length - want_len), uid))
    candidates = [uid for _, _, uid in sorted(ranked)]
    if not candidates:
        found = ", ".join(f"{n:,} bp" for n in seen_lengths[:6]) or "nothing"
        raise LookupError(f"no candidate for {term!r} is {want_len:,} bp "
                          f"(found: {found})")

    fasta = _eutils("efetch.fcgi", {
        "db": "nuccore", "id": ",".join(candidates),
        "rettype": "fasta", "retmode": "text"}).decode()
    records, header, seq = [], None, []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        else:
            seq.append(line.strip())
    if header:
        records.append((header, "".join(seq)))

    for head, sequence in sorted(
            records, key=lambda r: (0 if hint in _norm(r[0]) else 1,
                                    abs(len(r[1]) - want_len))):
        if abs(len(sequence) - want_len) > tolerance:
            continue
        if abs(_gc_percent(sequence) - want_gc) > GC_TOLERANCE_PCT:
            continue
        return head.split()[0], head, sequence
    found = ", ".join(f"{_gc_percent(s):.1f}% GC" for _, s in records[:6])
    raise LookupError(f"{term!r}: right length, wrong GC "
                      f"(want {want_gc:.2f}%, found: {found})")


def reference_set_phage_mock(refdir: Path, *, dry_run: bool) -> None:
    """The 15 phage genomes, fetched from NCBI and verified against Table S1."""
    fasta = refdir / "reference_genomes.fasta"
    if dry_run:
        log(f"  would resolve {len(PHAGE_MOCK_15)} phage genomes from NCBI")
        return
    if fasta.exists() and fasta.read_text().count(">") == len(PHAGE_MOCK_15):
        log(f"  have {fasta.name} ({len(PHAGE_MOCK_15)} genomes)")
    else:
        resolved, failed = [], []
        for label, term, want_len, want_gc in PHAGE_MOCK_15:
            try:
                acc, head, seq = resolve_phage_genome(term, want_len, want_gc)
            except LookupError as exc:
                failed.append((label, str(exc)))
                log(f"  MISS {label}: {exc}")
                continue
            resolved.append((label, acc, head, seq))
            log(f"  ok   {label} -> {acc} ({len(seq):,} bp)")
        if failed:
            raise SystemExit(
                "error: could not verify " + ", ".join(l for l, _ in failed) + "\n"
                "       Every genome must be confirmed by length and GC before it can\n"
                "       define ground truth. Fix the search term in PHAGE_MOCK_15, or\n"
                "       drop the sequence in by hand as >"
                "<label> and re-run - this script keeps an existing file.")
        with open(fasta, "w") as fh:
            for label, acc, head, seq in resolved:
                # Use the mock's own label as the id; keep the accession beside it.
                fh.write(f">{label} {head}\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + "\n")
    # No multipartite genomes here, so agent == molecule == genome.
    with open(refdir / "reference_map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reference_id", "molecule_label", "agent_label"])
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                ref = line[1:].split()[0]
                w.writerow([ref, ref, ref])
    log(f"  reference set: {len(PHAGE_MOCK_15)} genomes, 1 molecule each")


REFERENCE_SETS = {
    "dsmz_115_molecules": reference_set_dsmz,
    "phage_mock_15": reference_set_phage_mock,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dry-run", action="store_true",
                        help="print what would be fetched, download nothing")
    parser.add_argument("--reads-only", action="store_true",
                        help="skip the reference set")
    parser.add_argument("--references-only", action="store_true",
                        help="skip the reads")
    parser.add_argument("--discard-downloads", action="store_true",
                        help="delete each provider file once converted; halves "
                             "peak disk at the cost of re-downloading on a rebuild")
    args = parser.parse_args()

    dataset = os.environ.get("VB_DATASET", "")
    registry = load_registry(REGISTRY)
    if dataset not in registry:
        raise SystemExit(
            f"error: VB_DATASET must name a row in {REGISTRY}\n"
            f"       got {dataset!r}; known: {', '.join(sorted(registry))}")
    row = registry[dataset]

    study_dir = REPO_HOME / "Data" / "Mock_dataset" / row["study"]
    reads_dir = REPO_HOME / "Data" / "work" / f"mock_{dataset}" / "reads"
    refdir = study_dir / "references" / row["reference_set"]
    refdir.mkdir(parents=True, exist_ok=True)

    log(f"dataset {dataset}  provider={row['provider']}  "
        f"libraries={len([a for a in row['accessions'].split(',') if a])}")
    log(f"  reads      -> {reads_dir}")
    log(f"  references -> {refdir}")

    if not args.references_only:
        get_reads(row, reads_dir, study_dir, dry_run=args.dry_run,
                  discard_downloads=args.discard_downloads)
    if not args.reads_only:
        builder = REFERENCE_SETS.get(row["reference_set"])
        if builder is None:
            raise SystemExit(
                f"error: no builder for reference_set {row['reference_set']!r}\n"
                f"       known: {', '.join(sorted(REFERENCE_SETS))}")
        builder(refdir, dry_run=args.dry_run)

    if not args.dry_run:
        log("phase 1 complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
