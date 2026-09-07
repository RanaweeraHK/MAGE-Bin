#!/usr/bin/env python
from __future__ import annotations

import csv
import glob
import gzip
import hashlib
import io
import os
import random
from pathlib import Path


# ---- tier config (VB_TIER=low|medium|high, default high; see README's tier table)
TIER = os.environ.get("VB_TIER", "high")
if TIER not in ("low", "medium", "high"):
    raise SystemExit(f"error: VB_TIER must be low|medium|high (got {TIER!r})")
TIER_SUFFIX = {"low": "_low", "medium": "_medium", "high": "_high"}[TIER]

# x5 biomes -> 50 / 200 / 998 genomes (high loses 2 to the MD5 dedupe guard).
PER_BIOME = {"low": 10, "medium": 40, "high": 200}[TIER]

# .parent.parent.parent: this script lives in Dataset_Processing/Simulated_Data/,
REPO_HOME = Path(__file__).resolve().parent.parent.parent
RAW = REPO_HOME / "Data/Raw_dataset/multibiome_source_data"
POOL_DIR = REPO_HOME / f"Data/work/multibiome{TIER_SUFFIX}/genome_pool"
MIN_GENOME_LEN = 10000
POOL_SEED = 42


def open_text(path):
    """Encoding-tolerant FASTA reader (RefSeq's raw dump is UTF-16 + BOM)."""
    raw = gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")
    head = raw.read(2)
    raw.seek(0)
    if head in (b"\xff\xfe", b"\xfe\xff"):
        return io.TextIOWrapper(raw, encoding="utf-16", newline="")
    return io.TextIOWrapper(raw, encoding="utf-8-sig", errors="replace", newline="")

# fasta parser
def iter_fasta(path):
    name, seq = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


# r[0] = source ID , r[1] = sequence , r[2] = quality
def sample(records, n, min_len):
    """records: list of (source_id, seq, quality). Length-stratified random."""
    records = [r for r in records if len(r[1]) >= min_len]
    if len(records) <= n:
        return records
    return random.sample(records, n)

# load different biomes, gether only complete genomes 
def gather_reference():
    recs = []
    for sid, seq in iter_fasta(RAW / "reference" / "refseq_viral_raw.fasta"):
        if len(seq) >= MIN_GENOME_LEN:
            recs.append((sid, seq, "RefSeq-complete"))
    return recs


def gather_gut():
    recs = []
    for f in sorted(glob.glob(str(RAW / "human_gut_mgv" / "*.fna"))):
        for sid, seq in iter_fasta(f):
            recs.append((sid, seq, "MGV-complete"))
    return recs


def gather_freshwater():
    return [(sid, seq, "self-circular")
            for sid, seq in iter_fasta(RAW / "freshwater" / "freshwater_self_circular.fasta")]


def gather_soil():
    keep = {}
    with open(RAW / "soil" / "soil_GSV_quality.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("checkv_quality") in ("Complete", "High-quality"):
                keep[row["id"]] = row["checkv_quality"]
    recs = []
    for sid, seq in iter_fasta(RAW / "soil" / "soil_GSV.fasta"):
        if sid in keep:
            recs.append((sid, seq, keep[sid]))
    return recs


def gather_marine():
    return [(sid, seq, "IMGVR-HQ-marine")
            for sid, seq in iter_fasta(RAW / "marine_imgvr" / "marine_imgvr.fasta")]


BIOMES = {
    "reference": gather_reference,
    "gut": gather_gut,
    "freshwater": gather_freshwater,
    "soil": gather_soil,
    "marine": gather_marine,
}


def main() -> None:
    random.seed(POOL_SEED)
    print(f"tier: {TIER}  ({PER_BIOME} genomes/biome)  ->  {POOL_DIR}")
    POOL_DIR.mkdir(parents=True, exist_ok=True)
    pool_fasta = POOL_DIR / "community_pool.fasta"
    meta_tsv = POOL_DIR / "genome_metadata.tsv"

    seen_hashes: set[str] = set()
    rows = []
    with open(pool_fasta, "w") as pf:
        for biome, fn in BIOMES.items():
            cands = fn()
            chosen = sample(cands, PER_BIOME, MIN_GENOME_LEN)
            kept = 0
            for sid, seq, qual in sorted(chosen, key=lambda r: r[0]):
                h = hashlib.md5(seq.encode()).hexdigest()
                if h in seen_hashes:  # exact-duplicate guard across sources
                    continue
                seen_hashes.add(h)
                kept += 1
                gid = f"{biome}_{kept:03d}"
                pf.write(f">{gid}\n")
                for j in range(0, len(seq), 70):
                    pf.write(seq[j:j + 70] + "\n")
                rows.append({"genome_id": gid, "biome": biome, "source_id": sid,
                             "length": len(seq), "quality": qual})
            print(f"{biome:<11}: {len(cands):>6} candidates "
                  f"(>= {MIN_GENOME_LEN} bp)  ->  kept {kept}")

    with open(meta_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome_id", "biome", "source_id",
                                            "length", "quality"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    lens = [r["length"] for r in rows]
    print(f"\nPOOL: {len(rows)} genomes -> {pool_fasta}")
    if lens:
        print(f"  length: min={min(lens):,}  median={sorted(lens)[len(lens)//2]:,}  "
              f"max={max(lens):,}  total={sum(lens)/1e6:.1f} Mb")
    print(f"  metadata -> {meta_tsv}")


if __name__ == "__main__":
    main()
