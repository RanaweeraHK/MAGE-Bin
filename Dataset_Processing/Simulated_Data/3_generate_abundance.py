from __future__ import annotations

import math
import os
import random
from pathlib import Path

TIER = os.environ.get("VB_TIER", "high")
if TIER not in ("low", "medium", "high"):
    raise SystemExit(f"error: VB_TIER must be low|medium|high (got {TIER!r})")
TIER_SUFFIX = {"low": "_low", "medium": "_medium", "high": "_high"}[TIER]

# .parent.parent.parent: this script lives in Dataset_Processing/Simulated_Data/,
# so the repo root is two levels up.
WORK_DIR = (Path(__file__).resolve().parent.parent.parent
            / f"Data/work/multibiome{TIER_SUFFIX}")
POOL_DIR = WORK_DIR / "genome_pool"
ABUND_DIR = WORK_DIR / "abundance"
N_SAMPLES = 15
PRESENCE_PROB = 0.70  # each strain has 70% probability to being presents in any sample. (15*0.7 = 10.5 samples contain that strain)
LOGNORM_SIGMA = 2.0 # log-normal distribution as abundance pattern
ABUNDANCE_SEED = 42


def main() -> None:
    random.seed(ABUNDANCE_SEED)
    strain_fasta = POOL_DIR / "expanded_genomes.fasta"
    if not strain_fasta.is_file():
        raise SystemExit(f"error: {strain_fasta} not found - run 2_expand_strains.py first")

    strain_ids = [line[1:].split()[0] for line in open(strain_fasta) if line.startswith(">")]
    ABUND_DIR.mkdir(parents=True, exist_ok=True)

    n = N_SAMPLES
    raw = {s: [0.0] * n for s in strain_ids}
    for s in strain_ids:
        present = [random.random() < PRESENCE_PROB for _ in range(n)]
        if not any(present):                       # guarantee each strain appears once
            present[random.randrange(n)] = True
        for j in range(n):
            if present[j]:
                raw[s][j] = math.exp(random.gauss(0.0, LOGNORM_SIGMA))

    present_counts = []
    for j in range(n):
        col_sum = sum(raw[s][j] for s in strain_ids) or 1.0
        n_present = sum(1 for s in strain_ids if raw[s][j] > 0)
        present_counts.append(n_present)
        with open(ABUND_DIR / f"sample{j}_abundance.txt", "w") as fh:
            for s in strain_ids:
                fh.write(f"{s}\t{raw[s][j]/col_sum:.10f}\n")

    print(f"tier               : {TIER}")
    print(f"strains            : {len(strain_ids)}")
    print(f"samples            : {n}")
    print(f"presence prob      : {PRESENCE_PROB}")
    print(f"strains per sample : min={min(present_counts)} "
          f"median={sorted(present_counts)[n//2]} max={max(present_counts)}")
    print(f"abundance files    -> {ABUND_DIR}")


if __name__ == "__main__":
    main()
