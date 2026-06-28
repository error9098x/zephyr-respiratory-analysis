#!/usr/bin/env python3
"""
Track B, Steps 4-6 - embed the longest-ORF marker proteins with ESM2-8M, cluster,
and score against the Part 2 taxonomy (ARI/NMI), versus the sourmash baseline.

The markers (1.2k-4.5k aa) exceed ESM2's ~1022-residue context, so each protein is
split into 1000-aa windows, each window mean-pooled over residues, then the windows
are averaged into one whole-protein vector (L2-normalised). Clustering is k-means at
k = number of viruses (known from Part 2) plus Ward hierarchical, both compared to
the taxonomy. n is small (17), so we report it alongside the scores and lean on the
contingency table and plot, not the ARI number alone.
"""
import numpy as np, pandas as pd
import torch, esm
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, FIG = f"{ROOT}/results", f"{ROOT}/figures"
torch.manual_seed(0)

print("loading ESM2-8M ...")
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
model.eval()
bc = alphabet.get_batch_converter()
LAYER = 6

def embed(seq, win=1000):
    chunks = [seq[i:i+win] for i in range(0, len(seq), win)] or [seq]
    vecs = []
    with torch.no_grad():
        for j, ch in enumerate(chunks):
            _, _, toks = bc([(f"s{j}", ch)])
            rep = model(toks, repr_layers=[LAYER])["representations"][LAYER][0]
            vecs.append(rep[1:len(ch)+1].mean(0).numpy())   # exclude BOS/EOS
    v = np.mean(vecs, axis=0)
    return v / (np.linalg.norm(v) + 1e-9)

recs = list(SeqIO.parse(f"{RES}/longest_orf.faa", "fasta"))
labels = [r.id for r in recs]
virus = [l.split("__")[0].replace("_", " ") for l in labels]
X = np.vstack([embed(str(r.seq)) for r in recs])
print(f"embedded {X.shape[0]} markers, dim {X.shape[1]}")

k = len(set(virus))
P = PCA(n_components=min(10, X.shape[0]-1), random_state=0).fit_transform(X)
km = KMeans(n_clusters=k, n_init=20, random_state=0).fit_predict(P)
hc = fcluster(linkage(X, method="ward"), t=k, criterion="maxclust")

def score(name, cl):
    a, n = adjusted_rand_score(virus, cl), normalized_mutual_info_score(virus, cl)
    print(f"  {name:18s} ARI = {a:.3f}   NMI = {n:.3f}")
    return a, n

print(f"\nn = {len(labels)} markers, k = {k} viruses")
print("ESM2-8M (longest-ORF marker):")
ari_km, _ = score("k-means", km)
ari_hc, _ = score("ward-hier", hc)
print("\ncontingency (virus x k-means cluster):")
print(pd.crosstab(pd.Series(virus, name="virus"), pd.Series(km, name="cluster")).to_string())

# PCA 2D scatter coloured by virus
P2 = PCA(n_components=2, random_state=0).fit_transform(X)
fig, ax = plt.subplots(figsize=(9, 7))
vs = sorted(set(virus))
cmap = plt.get_cmap("tab20")
for i, v in enumerate(vs):
    m = [j for j, x in enumerate(virus) if x == v]
    ax.scatter(P2[m, 0], P2[m, 1], s=90, color=cmap(i % 20), label=v, edgecolor="k", linewidth=0.4)
ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
ax.set_title(f"ESM2-8M embeddings of longest-ORF markers (ARI {max(ari_km, ari_hc):.2f})")
ax.legend(fontsize=7, markerscale=0.8, ncol=2, loc="best")
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_esm_pca.png", dpi=200); plt.close(fig)

pd.DataFrame({"label": labels, "virus": virus, "kmeans": km, "ward": hc}).to_csv(
    f"{RES}/esm_clusters.csv", index=False)
print(f"\nWrote {FIG}/fig4_esm_pca.png and {RES}/esm_clusters.csv")
