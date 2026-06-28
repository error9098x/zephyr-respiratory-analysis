#!/usr/bin/env python3
"""
Track B, Step 1 - build one denoised consensus genome per selected detection.

For each kept (pool, virus), call samtools consensus over that virus's reference
using the Part 2 BAM. Positions with depth < 5 become N. We then QC each consensus
on length-vs-reference and N-fraction, drop the junk, and write the clean set to a
single FASTA (headers: virus__pool) for the downstream sketch / RdRp / clustering.
"""
import subprocess, os
import pandas as pd

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, REFS = f"{ROOT}/results", f"{ROOT}/refs"
os.makedirs(f"{RES}/consensus", exist_ok=True)

# virus -> reference accession (exact names as they appear in the BAM)
VREF = {
    "HCoV-229E": "NC_002645.1", "HCoV-HKU1": "NC_006577.2", "HPIV-4": "NC_021928.1",
    "RSV-A": "NC_038235.1", "RSV-B": "NC_001781.1",
    "Rhinovirus A": "NC_038311.1", "Rhinovirus B": "NC_038312.1", "Rhinovirus C": "NC_009996.1",
    "SARS-CoV-2": "NC_045512.2", "hMPV": "NC_039199.1",
}

# reference lengths
subprocess.run(["samtools", "faidx", f"{REFS}/respiratory_panel.fasta"], check=True)
reflen = {}
for ln in open(f"{REFS}/respiratory_panel.fasta.fai"):
    p = ln.split("\t"); reflen[p[0]] = int(p[1])

sel = pd.read_csv(f"{RES}/selected_detections.csv")
rows, seqs = [], {}
for _, r in sel.iterrows():
    ref = VREF[r.virus]
    label = f"{r.virus.replace(' ', '_')}__{r.pool}"
    bam = f"{RES}/{r.pool}.bam"
    out = subprocess.run(["samtools", "consensus", "-f", "fasta", "-r", ref, "-d", "5", bam],
                         capture_output=True, text=True, check=True).stdout
    seq = "".join(out.strip().split("\n")[1:])
    L = len(seq); Ns = seq.upper().count("N")
    nf = Ns / L if L else 1.0
    ratio = L / reflen[ref] if reflen.get(ref) else 0.0
    rows.append(dict(label=label, virus=r.virus, pool=r.pool, len=L, reflen=reflen.get(ref),
                     length_ratio=round(ratio, 2), N_frac=round(nf, 3)))
    seqs[label] = seq
    # also write individual file (handy for palmscan / inspection)
    with open(f"{RES}/consensus/{label}.fasta", "w") as fh:
        fh.write(f">{label}\n{seq}\n")

df = pd.DataFrame(rows)
df["keep"] = (df.N_frac <= 0.30) & (df.length_ratio >= 0.5)
print("=== consensus QC ===")
print(df[["label", "len", "reflen", "length_ratio", "N_frac", "keep"]].to_string(index=False))

kept = df[df.keep]
with open(f"{RES}/consensus_all.fasta", "w") as fh:
    for _, r in kept.iterrows():
        fh.write(f">{r.label}\n{seqs[r.label]}\n")

print(f"\nKept {len(kept)}/{len(df)} consensus genomes across {kept.virus.nunique()} viruses.")
if (~df.keep).any():
    print("Dropped:", ", ".join(df[~df.keep].label))
df.to_csv(f"{RES}/consensus_qc.csv", index=False)
print(f"Wrote {RES}/consensus_all.fasta and {RES}/consensus_qc.csv")
