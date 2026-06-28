#!/usr/bin/env python3
"""
Track B, Step 5 (extension) - does a larger protein language model resolve the
species that ESM2-8M could not? We embed the same longest-ORF markers with four
ESM2 sizes (8M, 35M, 150M, 650M), cluster each identically, and score against the
Part 2 taxonomy. We also report two targeted checks: whether RSV-A separates from
RSV-B, and a rhinovirus-only ARI (do the A/B/C species split apart).
"""
import gc, numpy as np, pandas as pd
import torch, esm
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
RES, FIG = f"{ROOT}/results", f"{ROOT}/figures"
torch.manual_seed(0); torch.set_num_threads(max(1, torch.get_num_threads()))

MODELS = [
    ("esm2_t6_8M_UR50D",   6,   "8M"),
    ("esm2_t12_35M_UR50D", 12,  "35M"),
    ("esm2_t30_150M_UR50D", 30, "150M"),
    ("esm2_t33_650M_UR50D", 33, "650M"),
]
PARAMS = {"8M": 8, "35M": 35, "150M": 150, "650M": 650}

recs = list(SeqIO.parse(f"{RES}/longest_orf.faa", "fasta"))
labels = [r.id for r in recs]
virus = [l.split("__")[0].replace("_", " ") for l in labels]
seqs = [str(r.seq) for r in recs]
k = len(set(virus))
rh = [i for i, v in enumerate(virus) if v.startswith("Rhinovirus")]

def embed_all(model_name, layer, win=1000):
    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval(); bc = alphabet.get_batch_converter()
    out = []
    for s in seqs:
        chunks = [s[i:i+win] for i in range(0, len(s), win)] or [s]
        vecs = []
        with torch.no_grad():
            for j, ch in enumerate(chunks):
                _, _, toks = bc([(f"s{j}", ch)])
                rep = model(toks, repr_layers=[layer])["representations"][layer][0]
                vecs.append(rep[1:len(ch)+1].mean(0).numpy())
        v = np.mean(vecs, axis=0)
        out.append(v / (np.linalg.norm(v) + 1e-9))
    del model; gc.collect()
    return np.vstack(out)

rows = []
for name, layer, tag in MODELS:
    print(f"\n=== {tag} ({name}) ===", flush=True)
    X = embed_all(name, layer)
    P = PCA(n_components=min(10, X.shape[0]-1), random_state=0).fit_transform(X)
    km = KMeans(n_clusters=k, n_init=20, random_state=0).fit_predict(P)
    ari = adjusted_rand_score(virus, km)
    nmi = normalized_mutual_info_score(virus, km)
    rsv_split = km[virus.index("RSV-A")] != km[virus.index("RSV-B")]
    rhino_ari = adjusted_rand_score([virus[i] for i in rh], [km[i] for i in rh])
    print(f"  dim {X.shape[1]:5d}  ARI {ari:.3f}  NMI {nmi:.3f}  RSV-A/B split: {rsv_split}  rhino-only ARI {rhino_ari:.3f}", flush=True)
    rows.append({"model": tag, "params_M": PARAMS[tag], "dim": X.shape[1],
                 "ARI": round(ari, 3), "NMI": round(nmi, 3),
                 "RSV_split": bool(rsv_split), "rhino_ARI": round(rhino_ari, 3)})

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/esm_modelsize_sweep.csv", index=False)
print("\n" + df.to_string(index=False))

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(df.params_M, df.ARI, "o-", label="overall ARI", lw=2)
ax.plot(df.params_M, df.NMI, "s--", label="overall NMI", lw=2)
ax.plot(df.params_M, df.rhino_ARI, "^:", label="rhinovirus-only ARI", lw=2)
ax.axhline(0.88, color="grey", ls="-", lw=1, label="sourmash baseline ARI")
ax.set_xscale("log"); ax.set_xticks(df.params_M); ax.set_xticklabels(df.model)
ax.set_xlabel("ESM2 model size (parameters)"); ax.set_ylabel("score")
ax.set_ylim(0, 1.02); ax.legend(); ax.grid(alpha=0.3)
ax.set_title("Does a larger protein model recover the taxonomy better?")
fig.tight_layout(); fig.savefig(f"{FIG}/fig6_esm_modelsize.png", dpi=200); plt.close(fig)
print(f"\nWrote {FIG}/fig6_esm_modelsize.png and {RES}/esm_modelsize_sweep.csv")
