#!/usr/bin/env python3
"""
Track B, Step 7 - skani ANI as a within-cluster validity check (NOT a cross-taxon
gradient). For each pair of consensus genomes we compute average nucleotide identity
with skani's small-genome (viral) preset. The point is twofold:
  (1) genomes the clustering put together (same virus) should sit at high ANI (>~95%),
      confirming they really are the same species and not a coincidental cluster;
  (2) genomes from different families should fall BELOW skani's 80% screen and report
      no value at all - i.e. ANI is undefined across distant taxa, so it cannot be used
      as a global similarity gradient. That undefined-ness is the expected, correct result.
"""
import os, subprocess, tempfile, numpy as np, pandas as pd
from Bio import SeqIO
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, FIG = f"{ROOT}/results", f"{ROOT}/figures"
SKANI = "/opt/homebrew/bin/skani"

recs = list(SeqIO.parse(f"{RES}/consensus_all.fasta", "fasta"))
labels = [r.id for r in recs]
virus = [l.split("__")[0].replace("_", " ") for l in labels]
order = sorted(range(len(labels)), key=lambda i: (virus[i], labels[i]))
labels = [labels[i] for i in order]; virus = [virus[i] for i in order]
recs = [recs[i] for i in order]

tmp = tempfile.mkdtemp()
files = []
for r in recs:
    p = os.path.join(tmp, r.id + ".fa")
    with open(p, "w") as fh:
        fh.write(f">{r.id}\n{str(r.seq)}\n")
    files.append(p)

out = os.path.join(tmp, "ani.txt")
cmd = [SKANI, "triangle", "--small-genomes", "--full-matrix", "--diagonal",
       "-s", "80", "-t", "4", "-o", out] + files
res = subprocess.run(cmd, capture_output=True, text=True)
print("skani stderr (tail):", res.stderr.strip().splitlines()[-2:] if res.stderr.strip() else "none")

# parse full matrix: first line = N, then rows "name v1 v2 ..."
lines = [l for l in open(out).read().splitlines() if l.strip()]
n = int(lines[0])
names, M = [], []
for ln in lines[1:]:
    parts = ln.split("\t") if "\t" in ln else ln.split()
    names.append(os.path.basename(parts[0]).rsplit(".fa", 1)[0])
    M.append([float(x) if x not in ("NA", "") else 0.0 for x in parts[1:]])
M = np.array(M)
idx = {nm: i for i, nm in enumerate(names)}
A = np.array([[M[idx[a]][idx[b]] for b in labels] for a in labels])  # reorder to labels

# summarise within- vs across-species
within, across_def = [], 0
for i in range(len(labels)):
    for j in range(i+1, len(labels)):
        v = A[i, j]
        if virus[i] == virus[j]:
            within.append((f"{virus[i]}", labels[i], labels[j], v))
        elif v > 0:
            across_def += 1

print("\nWITHIN-species genome pairs (should be high ANI):")
for vname, a, b, v in within:
    tag = "" if v >= 95 else ("  <-- below 95%" if v > 0 else "  <-- NO ANI (screened)")
    print(f"  {vname:14s} {v:6.2f}%   {a}  vs  {b}{tag}")
total_cross = sum(1 for i in range(len(labels)) for j in range(i+1, len(labels)) if virus[i] != virus[j])
print(f"\nACROSS-species pairs with a reported ANI (expect ~0): {across_def} / {total_cross} (rest screened out = undefined, as intended)")

# heatmap
fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(np.where(A > 0, A, np.nan), cmap="viridis", vmin=80, vmax=100)
ax.set_xticks(range(len(labels))); ax.set_yticks(range(len(labels)))
ax.set_xticklabels(virus, rotation=90, fontsize=7); ax.set_yticklabels(virus, fontsize=7)
ax.set_title("skani ANI (small-genome preset); grey = below 80% screen = undefined")
ax.set_facecolor("0.85")
fig.colorbar(im, ax=ax, fraction=0.046, label="ANI %")
fig.tight_layout(); fig.savefig(f"{FIG}/fig5_ani_heatmap.png", dpi=200); plt.close(fig)
pd.DataFrame(A, index=labels, columns=labels).to_csv(f"{RES}/ani_matrix.csv")
print(f"\nWrote {FIG}/fig5_ani_heatmap.png and {RES}/ani_matrix.csv")
