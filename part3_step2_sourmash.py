#!/usr/bin/env python3
"""
Track B, Step 2 - sourmash MinHash baseline.

Sketch each consensus genome, build an all-vs-all similarity matrix, cluster it,
and score the clusters against the Part 2 taxonomy. This is the simple,
alignment-free nucleotide baseline against which the ESM2/RdRp tier is compared.
"""
import subprocess, os, sys
import numpy as np, pandas as pd
SM = os.path.join(os.path.dirname(sys.executable), "sourmash")  # venv sourmash
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, FIG = f"{ROOT}/results", f"{ROOT}/figures"
KMER, SCALED = 21, 50   # small scaled because viral genomes are tiny

sh = lambda c: subprocess.run(c, shell=True, check=True)
sh(f'"{SM}" sketch dna -p k={KMER},scaled={SCALED} --singleton "{RES}/consensus_all.fasta" -o "{RES}/sketch.zip" -f 2>/dev/null')
sh(f'"{SM}" compare "{RES}/sketch.zip" -k {KMER} --csv "{RES}/sourmash_compare.csv" 2>/dev/null')

M = pd.read_csv(f"{RES}/sourmash_compare.csv")
labels = list(M.columns)
S = M.values.astype(float)
virus = [l.split("__")[0].replace("_", " ") for l in labels]
k = len(set(virus))

# distance + average-linkage hierarchical clustering at k = #viruses
D = 1.0 - S; np.fill_diagonal(D, 0.0); D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), method="average")
clusters = fcluster(Z, t=k, criterion="maxclust")

ari = adjusted_rand_score(virus, clusters)
nmi = normalized_mutual_info_score(virus, clusters)
print(f"n = {len(labels)} genomes, k = {k} viruses")
print(f"sourmash baseline:  ARI = {ari:.3f}   NMI = {nmi:.3f}\n")
print("contingency (virus x cluster):")
print(pd.crosstab(pd.Series(virus, name="virus"), pd.Series(clusters, name="cluster")).to_string())

# heatmap ordered by virus
order = np.argsort(virus)
So = S[np.ix_(order, order)]
vlab = [virus[i] for i in order]
fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(So, cmap="magma", vmin=0, vmax=1)
ax.set_xticks(range(len(vlab))); ax.set_xticklabels([labels[i] for i in order], rotation=90, fontsize=6)
ax.set_yticks(range(len(vlab))); ax.set_yticklabels([labels[i] for i in order], fontsize=6)
fig.colorbar(im, ax=ax, label="MinHash Jaccard similarity")
ax.set_title(f"sourmash genome similarity (k={KMER})  -  baseline ARI {ari:.2f}")
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_sourmash_heatmap.png", dpi=200); plt.close(fig)

pd.DataFrame({"label": labels, "virus": virus, "cluster": clusters}).to_csv(
    f"{RES}/sourmash_clusters.csv", index=False)
print(f"\nWrote {FIG}/fig3_sourmash_heatmap.png and {RES}/sourmash_clusters.csv")
