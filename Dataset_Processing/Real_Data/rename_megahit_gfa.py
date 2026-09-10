#!/usr/bin/env python3
"""Rename a MEGAHIT GFA's segments to the contig names used in its FASTA.

MEGAHIT publishes no GFA.  The usual recipe runs `megahit_toolkit contig2fastg`
and converts the result, and that tool SYNTHESISES SPAdes-style names:

    GFA    S  NODE_1_length_315_cov_1.0000_ID_1  GATCC...
    FASTA  >  k141_1325870 flag=0 multi=1.0000 len=315

The two share nothing.  Phase 2 has no path table for MEGAHIT, so it matches
GFA segment names against contig names directly - and against these names it
matches NOTHING, producing a graph with zero edges.  That failure is silent: it
looks exactly like a genuinely fragmented assembly, which this pipeline has
already seen for real once (wastewater, 10 edges over 7,842 contigs).

Segments are therefore renamed by SEQUENCE, checking both strands.  The
correspondence happens to be positional in the files seen so far, but position
is an assumption about tool internals while sequence identity is a fact.

    python rename_megahit_gfa.py final.contigs.fa final.contigs.gfa out.gfa
"""
import sys


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        yield name, "".join(chunks)


def main(fasta, gfa_in, gfa_out):
    by_seq = {}
    for name, seq in read_fasta(fasta):
        by_seq[seq] = name
    print(f"contigs in FASTA        : {len(by_seq):,}")

    # Pass 1: segment name -> contig name, via the segment's own sequence.
    rename, segments, unmatched = {}, 0, 0
    with open(gfa_in) as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue
            segments += 1
            _, sid, seq = line.rstrip("\n").split("\t")[:3]
            hit = by_seq.get(seq) or by_seq.get(revcomp(seq))
            if hit:
                rename[sid] = hit
            else:
                unmatched += 1
    print(f"segments in GFA         : {segments:,}")
    print(f"  matched to a contig   : {len(rename):,}")
    print(f"  unmatched (dropped)   : {unmatched:,}")

    # Pass 2: rewrite, dropping segments and links that have no contig.
    kept_s = kept_l = dropped_l = 0
    with open(gfa_in) as fh, open(gfa_out, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "S":
                new = rename.get(parts[1])
                if new:
                    parts[1] = new
                    out.write("\t".join(parts) + "\n")
                    kept_s += 1
            elif parts[0] == "L":
                a, b = rename.get(parts[1]), rename.get(parts[3])
                if a and b:
                    parts[1], parts[3] = a, b
                    out.write("\t".join(parts) + "\n")
                    kept_l += 1
                else:
                    dropped_l += 1
            else:
                out.write(line)
    print(f"written                 : {kept_s:,} S, {kept_l:,} L "
          f"({dropped_l:,} link(s) dropped)")
    covered = 100.0 * kept_s / len(by_seq) if by_seq else 0.0
    print(f"contig coverage         : {covered:.1f}%")
    if covered < 50.0:
        print("WARNING: fewer than half the contigs appear in the graph")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
