#!/usr/bin/env python3
"""
Track B, Step 3 (substituted) - extract the longest ORF per consensus genome as a
replicase/polymerase proxy marker.

palmscan (the standard RdRp palmprint tool) would not build on Apple Silicon, so
instead of a true RdRp we take the longest stop-free reading frame across all six
frames of each consensus. For these RNA viruses the longest ORF IS the replicase:
the picornavirus polyprotein, the coronavirus ORF1ab, and the paramyxo/pneumovirus
L gene all carry the RdRp. It keeps a single comparable protein per genome for the
ESM2 tier. Ns translate to X (not a stop), so partial genomes are not fragmented.
"""
from Bio import SeqIO
from Bio.Seq import Seq

RES = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest/results"
MIN_AA = 200

def longest_orf(seq):
    best = ""
    s = Seq(str(seq).upper().replace("-", "N"))
    for strand in (s, s.reverse_complement()):
        n = len(strand)
        for frame in range(3):
            aa = str(strand[frame: frame + (n - frame) // 3 * 3].translate(table=1))
            for piece in aa.split("*"):
                if len(piece) > len(best):
                    best = piece
    return best

recs = list(SeqIO.parse(f"{RES}/consensus_all.fasta", "fasta"))
rows, kept = [], []
for r in recs:
    orf = longest_orf(r.seq)
    rows.append((r.id, len(r.seq), len(orf)))
    if len(orf) >= MIN_AA:
        kept.append((r.id, orf))

with open(f"{RES}/longest_orf.faa", "w") as fh:
    for i, o in kept:
        fh.write(f">{i}\n{o}\n")

print(f"{'genome':45s} {'genome_bp':>9} {'longest_orf_aa':>14}")
for i, gl, ol in rows:
    print(f"{i:45s} {gl:9d} {ol:14d}")
print(f"\nkept {len(kept)}/{len(rows)} markers with ORF >= {MIN_AA} aa -> {RES}/longest_orf.faa")
