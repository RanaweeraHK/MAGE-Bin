#!/usr/bin/env python
from __future__ import annotations

import csv
import os
import random
from pathlib import Path

TIER = os.environ.get("VB_TIER", "high")
if TIER not in ("low", "medium", "high"):
    raise SystemExit(f"error: VB_TIER must be low|medium|high (got {TIER!r})")
TIER_SUFFIX = {"low": "_low", "medium": "_medium", "high": "_high"}[TIER]
STRAIN_FRACTION = {"low": 0.05, "medium": 0.20, "high": 0.40}[TIER]

# parameter sets
POOL_DIR = (Path(__file__).resolve().parent.parent.parent
            / f"Data/work/multibiome{TIER_SUFFIX}/genome_pool")
STRAINS_MIN = 2
STRAINS_MAX = 5
ANI_MIN = 98.0
ANI_MAX = 99.5
STRAIN_SEED = 42

BASES = "ACGT"


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def mutate(seq: str, ani: float) -> str:
    """Variant of seq at ~ani% nucleotide identity. it gets original sequence+ desired ANI"""
    sub_rate = (100 - ani) / 100.0
    out = []
    for base in seq:
        r = random.random()
        if r < sub_rate * 0.90:                 # substitution
            out.append(random.choice(BASES.replace(base, "") if base in BASES else BASES))
        elif r < sub_rate * 0.95:               # deletion (skip base)
            continue
        elif r < sub_rate * 1.00:               # insertion
            out.append(base)
            out.append(random.choice(BASES))
        else:
            out.append(base)
    return "".join(out)


def main() -> None:
    random.seed(STRAIN_SEED)
    in_fasta = POOL_DIR / "community_pool.fasta"
    in_meta = POOL_DIR / "genome_metadata.tsv"
    out_fasta = POOL_DIR / "expanded_genomes.fasta"
    out_map = POOL_DIR / "strain_map.tsv"

    if not in_fasta.is_file():
        raise SystemExit(f"error: {in_fasta} not found - run 1_build_genome_pool.py first")

    meta = {r["genome_id"]: r for r in csv.DictReader(open(in_meta), delimiter="\t")}
    genomes = list(read_fasta(in_fasta))
    n_multi = int(round(len(genomes) * STRAIN_FRACTION))
    multi_ids = set(random.sample([g for g, _ in genomes], n_multi))

    rows = []
    total_strains = 0
    with open(out_fasta, "w") as pf:
        for gid, seq in genomes:
            biome = meta[gid]["biome"]
            n = random.randint(STRAINS_MIN, STRAINS_MAX) if gid in multi_ids else 1
            for k in range(n):
                if k == 0:
                    variant, ani = seq, 100.0            # strain 0 = original
                else:
                    ani = round(random.uniform(ANI_MIN, ANI_MAX), 2)
                    variant = mutate(seq, ani)
                ssid = f"{gid}_s{k}"
                pf.write(f">{ssid}\n")
                for j in range(0, len(variant), 70):
                    pf.write(variant[j:j + 70] + "\n")
                rows.append({"strain_id": ssid, "species_id": gid, "biome": biome,
                             "ani_to_base": ani, "length": len(variant),
                             "n_strains_in_species": n})
                total_strains += 1

    with open(out_map, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["strain_id", "species_id", "biome",
                                            "ani_to_base", "length",
                                            "n_strains_in_species"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"tier                       : {TIER} "
          f"(strain complexity {100*STRAIN_FRACTION:.0f}%)")
    print(f"species (pool genomes)     : {len(genomes)}")
    print(f"  multi-strain species     : {n_multi} ({100*n_multi/len(genomes):.0f}%)")
    print(f"  single-strain species    : {len(genomes) - n_multi}")
    print(f"total strain sequences     : {total_strains}")
    print(f"  extra (variant) strains  : {total_strains - len(genomes)}")
    print(f"expanded genomes -> {out_fasta}")
    print(f"strain map       -> {out_map}")


if __name__ == "__main__":
    main()
