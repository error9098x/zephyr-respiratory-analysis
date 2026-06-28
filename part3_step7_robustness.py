#!/usr/bin/env python3
"""
Track B robustness - is the clustering agreement real given only 17 genomes?
Two checks for both the sourmash and ESM2-8M arms:
  (1) Permutation null: shuffle the virus labels many times and recompute ARI against
      the fixed clustering. The observed ARI should sit far in the tail (small p).
  (2) Leave-one-out stability: drop each genome, recluster the rest, recompute ARI.
      A stable method gives a tight spread, not a number propped up by one or two points.
"""
import numpy as np, pandas as pd
import torch, esm
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, FIG = f"{ROOT}/results", f"{ROOT}/figures"
rng = np.random.default_rng(0); torch.manual_seed(0)
NPERM = 2000

# --- sourmash arm: load similarity matrix ---
M = pd.read_csv(f"{RES}/sourmash_compare.csv")
labels = list(M.columns); S = M.values.astype(float)
virus = np.array([l.split("__")[0].replace("_", " ") for l in labels])
n = len(labels)

def sourmash_cluster(idx):
    sub = S[np.ix_(idx, idx)]; v = virus[idx]; kk = len(set(v))
    D = 1.0 - sub; np.fill_diagonal(D, 0.0); D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method="average")
    return fcluster(Z, t=kk, criterion="maxclust")

# --- ESM arm: embed 8M once ---
recs = {r.id: str(r.seq) for r in SeqIO.parse(f"{RES}/longest_orf.faa", "fasta")}
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D(); model.eval(); bc = alphabet.get_batch_converter()
def embed(seq, win=1000):
    chunks = [seq[i:i+win] for i in range(0, len(seq), win)] or [seq]; vs = []
    with torch.no_grad():
        for j, ch in enumerate(chunks):
            _, _, t = bc([(f"s{j}", ch)])
            r = model(t, repr_layers=[6])["representations"][6][0]
            vs.append(r[1:len(ch)+1].mean(0).numpy())
    v = np.mean(vs, axis=0); return v / (np.linalg.norm(v) + 1e-9)
X = np.vstack([embed(recs[l]) for l in labels])

def esm_cluster(idx):
    sub = X[idx]; v = virus[idx]; kk = len(set(v))
    P = PCA(n_components=min(10, len(idx)-1), random_state=0).fit_transform(sub)
    return KMeans(n_clusters=kk, n_init=20, random_state=0).fit_predict(P)

allidx = np.arange(n)
arms = {"sourmash": sourmash_cluster, "ESM2-8M": esm_cluster}
summary = {}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (name, fn) in zip(axes, arms.items()):
    cl = fn(allidx)
    obs = adjusted_rand_score(virus, cl)
    null = np.array([adjusted_rand_score(rng.permutation(virus), cl) for _ in range(NPERM)])
    p = (np.sum(null >= obs) + 1) / (NPERM + 1)
    loo = np.array([adjusted_rand_score(virus[[i for i in range(n) if i != d]],
                    fn([i for i in range(n) if i != d])) for d in range(n)])
    summary[name] = dict(observed=round(obs, 3), perm_p=p,
                         null_mean=round(null.mean(), 3), null_95=round(np.quantile(null, 0.95), 3),
                         loo_mean=round(loo.mean(), 3), loo_min=round(loo.min(), 3), loo_max=round(loo.max(), 3))
    ax.hist(null, bins=30, color="0.7", edgecolor="none")
    ax.axvline(obs, color="crimson", lw=2, label=f"observed {obs:.2f}")
    ax.set_title(f"{name}: permutation null (p = {p:.4f})", fontsize=11)
    ax.set_xlabel("ARI under shuffled labels"); ax.set_ylabel("count"); ax.legend()
fig.suptitle("Is the clustering agreement real at n=17? Permutation null vs observed", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/fig7_robustness.png", dpi=200); plt.close(fig)

print(pd.DataFrame(summary).T.to_string())
pd.DataFrame(summary).T.to_csv(f"{RES}/robustness.csv")
print(f"\nWrote {FIG}/fig7_robustness.png and {RES}/robustness.csv")
print("\nLeave-one-out reading: spread = how much ARI moves when any single genome is removed.")
