#!/usr/bin/env python3
"""Fetch a panel of reference rhinovirus genomes (multiple types per species) from
NCBI, so our swab consensuses can be placed among the known rhinovirus diversity."""
import time
from Bio import Entrez, SeqIO

Entrez.email = "zephyr.worktest@example.com"
OUT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest/refs/rhino_refs.fasta"
TERMS = {
    "RV-A": "Rhinovirus A[Organism] AND complete genome AND 6400:7600[SLEN]",
    "RV-B": "Rhinovirus B[Organism] AND complete genome AND 6400:7600[SLEN]",
    "RV-C": "Rhinovirus C[Organism] AND complete genome AND 6400:7600[SLEN]",
}
RETMAX = 16

written = 0
with open(OUT, "w") as out:
    for sp, term in TERMS.items():
        ids = Entrez.read(Entrez.esearch(db="nucleotide", term=term, retmax=RETMAX))["IdList"]
        time.sleep(0.5)
        if not ids:
            print(f"{sp}: no ids"); continue
        recs = list(SeqIO.parse(Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text"), "fasta"))
        time.sleep(0.5)
        for r in recs:
            acc = r.id.split(".")[0]
            out.write(f">{sp}__{acc}\n{str(r.seq)}\n")
            written += 1
        print(f"{sp}: {len(recs)} genomes")
print(f"\nWrote {written} reference genomes to {OUT}")
