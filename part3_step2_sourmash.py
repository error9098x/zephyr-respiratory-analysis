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

# heatmap ordered by virus then cluster, with side strips for true virus + assigned cluster
varr = np.array(virus)
order = np.lexsort((clusters, varr))
So = S[np.ix_(order, order)]
vlab = [virus[i] for i in order]
clab = [int(clusters[i]) for i in order]
uniq = sorted(set(virus)); vcode = np.array([uniq.index(v) for v in vlab])
cuniq = sorted(set(clab));  ccode = np.array([cuniq.index(c) for c in clab])

fig = plt.figure(figsize=(11, 8.2))
gs = fig.add_gridspec(1, 4, width_ratios=[0.5, 0.5, 12, 0.35], wspace=0.05)
ax_v, ax_c, ax_h, cax = (fig.add_subplot(gs[0, i]) for i in range(4))

ax_v.imshow(vcode[:, None], aspect="auto", cmap="tab20")
ax_v.set_xticks([]); ax_v.set_yticks(range(len(vlab))); ax_v.set_yticklabels(vlab, fontsize=7)
ax_v.set_title("true", fontsize=8)

ax_c.imshow(ccode[:, None], aspect="auto", cmap="tab20b")
ax_c.set_xticks([]); ax_c.set_yticks([]); ax_c.set_title("clust", fontsize=8)
for i, c in enumerate(clab):
    ax_c.text(0, i, str(c), ha="center", va="center", fontsize=8, fontweight="bold", color="white")

im = ax_h.imshow(So, cmap="magma", vmin=0, vmax=1)
ax_h.set_yticks([])
ax_h.set_xticks(range(len(vlab))); ax_h.set_xticklabels(vlab, rotation=90, fontsize=7)
for b in [i for i in range(1, len(vlab)) if vlab[i] != vlab[i-1]]:
    ax_h.axhline(b-0.5, color="white", lw=0.6); ax_h.axvline(b-0.5, color="white", lw=0.6)
fig.colorbar(im, cax=cax, label="MinHash Jaccard similarity")
ax_h.set_title(f"sourmash similarity, ordered by virus  -  ARI {ari:.2f}\n"
               f"side strips: true virus vs assigned cluster (same number = same cluster)", fontsize=10)
fig.savefig(f"{FIG}/fig3_sourmash_heatmap.png", dpi=200, bbox_inches="tight"); plt.close(fig)

pd.DataFrame({"label": labels, "virus": virus, "cluster": clusters}).to_csv(
    f"{RES}/sourmash_clusters.csv", index=False)
print(f"\nWrote {FIG}/fig3_sourmash_heatmap.png and {RES}/sourmash_clusters.csv")
