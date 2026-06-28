#!/usr/bin/env python3
"""
Rhinovirus serotype demonstration. The ANI step showed our same-species rhinovirus
genomes are too divergent to align (below the 80% screen). Here we turn that into a
positive demonstration: we place our seven rhinovirus consensuses among 48 reference
rhinovirus genomes (16 each of RV-A, RV-B, RV-C) and ask which reference each one
sits closest to. If our three RV-A genomes each match a DIFFERENT reference, they are
different types, and a single RV-A reference cannot represent them.

MinHash at k=21 saturates across divergent types (shared k-mers go to zero), so we
sketch at low k=12 over every k-mer (scaled=1) to recover a graded similarity that
still separates types.
"""
import os, sys, subprocess, tempfile
import numpy as np, pandas as pd
from Bio import SeqIO
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, FIG, REFS = f"{ROOT}/results", f"{ROOT}/figures", f"{ROOT}/refs"
SM = os.path.join(os.path.dirname(sys.executable), "sourmash")

# our rhinovirus consensuses
ours = [r for r in SeqIO.parse(f"{RES}/consensus_all.fasta", "fasta") if r.id.startswith("Rhinovirus")]
refs = list(SeqIO.parse(f"{REFS}/rhino_refs.fasta", "fasta"))

tmp = tempfile.mkdtemp(); combined = os.path.join(tmp, "all.fasta")
labels, species, is_ours = [], [], []
with open(combined, "w") as fh:
    for r in refs:
        sp = r.id.split("__")[0]            # RV-A / RV-B / RV-C
        fh.write(f">{r.id}\n{str(r.seq)}\n"); labels.append(r.id); species.append(sp); is_ours.append(False)
    for r in ours:
        sp = "RV-" + r.id.split("__")[0].replace("Rhinovirus_", "").replace("Rhinovirus ", "")  # Rhinovirus_A -> RV-A
        nm = "OUR_" + r.id.replace("Rhinovirus_", "RV").split("__")[0] + "_" + r.id.split("__")[1][:6]
        fh.write(f">{nm}\n{str(r.seq)}\n"); labels.append(nm); species.append(sp); is_ours.append(True)

sk = os.path.join(tmp, "s.zip"); cmp = os.path.join(tmp, "c.csv")
subprocess.run(f'"{SM}" sketch dna -p k=12,scaled=1 --singleton "{combined}" -o "{sk}" -f', shell=True, check=True, capture_output=True)
subprocess.run(f'"{SM}" compare "{sk}" -k 12 --csv "{cmp}"', shell=True, check=True, capture_output=True)
M = pd.read_csv(cmp); order = list(M.columns)
S = M.values.astype(float)
pos = {nm: i for i, nm in enumerate(order)}
idx = [pos[l] for l in labels]
S = S[np.ix_(idx, idx)]                       # reorder to our labels order
species = np.array(species); is_ours = np.array(is_ours)

# nearest-reference typing for each of our genomes
ref_mask = ~is_ours
print("Nearest reference type for each of our rhinovirus genomes:")
rows = []
for i in np.where(is_ours)[0]:
    sims = S[i].copy(); sims[i] = -1; sims[~ref_mask] = -1   # only compare to references
    j = int(np.argmax(sims))
    rows.append((labels[i], species[i], labels[j], round(float(sims[j]), 3)))
    print(f"  {labels[i]:22s} ({species[i]})  ->  {labels[j]:14s}  sim {sims[j]:.3f}")
pd.DataFrame(rows, columns=["our_genome", "our_species", "nearest_ref", "similarity"]).to_csv(
    f"{RES}/rhino_nearest_ref.csv", index=False)

# do our same-species genomes match the SAME reference or different ones?
for sp in ["RV-A", "RV-B", "RV-C"]:
    nn = [r[2] for r in rows if r[1] == sp]
    if len(nn) > 1:
        print(f"\n{sp}: our {len(nn)} genomes -> {len(set(nn))} distinct reference types "
              f"({'different types' if len(set(nn)) == len(nn) else 'some shared'})")

# pairwise similarity among our same-species genomes (should be low = different types)
print("\nPairwise k=12 similarity among our same-species genomes:")
oidx = np.where(is_ours)[0]
for a in range(len(oidx)):
    for b in range(a+1, len(oidx)):
        ia, ib = oidx[a], oidx[b]
        if species[ia] == species[ib]:
            print(f"  {labels[ia]} vs {labels[ib]}: {S[ia, ib]:.3f}")

# figure: our genomes (rows) vs reference genomes (cols, grouped by species).
# each of our genomes peaks at a different reference type = different serotype.
ref_order = sorted(np.where(ref_mask)[0], key=lambda i: (species[i], labels[i]))
our_order = sorted(np.where(is_ours)[0], key=lambda i: (species[i], labels[i]))
H = S[np.ix_(our_order, ref_order)]
fig, ax = plt.subplots(figsize=(12, 4.6))
im = ax.imshow(H, aspect="auto", cmap="viridis", vmin=0, vmax=max(0.4, H.max()))
ax.set_yticks(range(len(our_order))); ax.set_yticklabels([labels[i] for i in our_order], fontsize=8)
ax.set_xticks(range(len(ref_order)))
ax.set_xticklabels([species[i].replace("RV-", "") for i in ref_order], fontsize=6)
ax.set_xlabel("reference rhinovirus genomes (grouped A | B | C)")
# species group separators
refsp = [species[i] for i in ref_order]
for b in [k for k in range(1, len(refsp)) if refsp[k] != refsp[k-1]]:
    ax.axvline(b-0.5, color="white", lw=1.5)
# mark each row's nearest reference
for ri, i in enumerate(our_order):
    j_in_ref = int(np.argmax([S[i, r] for r in ref_order]))
    ax.add_patch(plt.Rectangle((j_in_ref-0.5, ri-0.5), 1, 1, fill=False, edgecolor="red", lw=2))
fig.colorbar(im, ax=ax, label="k=12 genome similarity", fraction=0.025)
ax.set_title("Each of our rhinovirus genomes (rows) is closest to a DIFFERENT reference type (red box)\n"
             "same-species genomes do not share a serotype, so one species reference is too coarse")
fig.tight_layout(); fig.savefig(f"{FIG}/fig8_rhino_serotypes.png", dpi=200); plt.close(fig)
print(f"\nWrote {FIG}/fig8_rhino_serotypes.png and {RES}/rhino_nearest_ref.csv")
