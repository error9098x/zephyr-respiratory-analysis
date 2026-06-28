#!/usr/bin/env python3
"""
Part 2 - Taxonomic classification and genome coverage for Zephyr respiratory pools.

For each swab pool we competitively map the published viral reads against a curated
respiratory-virus reference panel with minimap2 (Nanopore preset), keep primary
alignments only, and summarise per-virus read counts and genome coverage
(breadth and depth) with samtools. Breadth of coverage (the percent of the genome
covered) is treated as the main "is it really present" signal, read count and
mapping quality as supporting evidence.

Reasoning: the reads are already pre-filtered to vertebrate-infecting viruses, so
the job is to assign them to specific viruses and measure coverage, not to profile
a whole community. Mapping to a reference panel does both; a k-mer classifier would
give abundance but not coverage. Competitive mapping plus primary-only counting
stops one read being counted against several closely related references.
"""
import subprocess, glob, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest"
DATA, REFS, RES, FIG = f"{ROOT}/data", f"{ROOT}/refs", f"{ROOT}/results", f"{ROOT}/figures"
PANEL = f"{REFS}/respiratory_panel.fasta"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# ---- accession -> virus name (groups the influenza segments) ----
FLU_A = {f"NC_0264{n}" for n in range(31, 39)}   # NC_026431..NC_026438
FLU_B = {f"NC_0022{n:02d}" for n in range(4, 12)} # NC_002204..NC_002211
FLU_C = {f"NC_0063{n:02d}" for n in range(7, 13)} # NC_006307..NC_006312
SINGLE = {
    "NC_045512": "SARS-CoV-2", "NC_038235": "RSV-A", "NC_001781": "RSV-B",
    "NC_002645": "HCoV-229E", "NC_006213": "HCoV-OC43", "NC_005831": "HCoV-NL63",
    "NC_006577": "HCoV-HKU1", "NC_039199": "hMPV",
    "NC_003461": "HPIV-1", "NC_003443": "HPIV-2", "NC_001796": "HPIV-3", "NC_021928": "HPIV-4",
    "NC_038311": "Rhinovirus A", "NC_038312": "Rhinovirus B", "NC_009996": "Rhinovirus C",
}
def virus_of(rname):
    a = rname.split(".")[0]
    if a in SINGLE: return SINGLE[a]
    if a in FLU_A: return "Influenza A"
    if a in FLU_B: return "Influenza B"
    if a in FLU_C: return "Influenza C"
    return a

def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)

# ---- build minimap2 index once ----
print("Indexing reference panel ...")
sh(f'minimap2 -x map-ont -d "{REFS}/panel.mmi" "{PANEL}" 2>/dev/null')

pools = sorted(glob.glob(f"{DATA}/*.respiratory.fasta.gz"))
print(f"Found {len(pools)} pools")

rows = []
for p in pools:
    name = os.path.basename(p).replace(".respiratory.fasta.gz", "")
    bam = f"{RES}/{name}.bam"
    # competitive map, Nanopore preset, no secondary alignments, -L for safe long CIGARs
    sh(f'minimap2 -ax map-ont -L --secondary=no "{REFS}/panel.mmi" "{p}" 2>/dev/null '
       f'| samtools sort -o "{bam}" - 2>/dev/null')
    sh(f'samtools index "{bam}"')
    # coverage, excluding unmapped/secondary/supplementary/qcfail/dup
    cov = subprocess.run(
        ["samtools", "coverage", "--ff", "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY", bam],
        capture_output=True, text=True, check=True).stdout
    for line in cov.strip().split("\n")[1:]:
        f = line.split("\t")
        # columns: rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
        rows.append(dict(pool=name, rname=f[0], reflen=int(f[2]), numreads=int(f[3]),
                         covbases=int(f[4]), breadth=float(f[5]), meandepth=float(f[6]),
                         meanmapq=float(f[8]), virus=virus_of(f[0])))
    print(f"  mapped {name}")

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/coverage_by_reference.csv", index=False)

# ---- aggregate to virus level per pool (sums influenza segments) ----
def agg_group(g):
    reads = int(g.numreads.sum())
    L = g.reflen.sum()
    breadth = 100.0 * g.covbases.sum() / L if L else 0.0
    depth = (g.meandepth * g.reflen).sum() / L if L else 0.0
    mapq = (g.meanmapq * g.numreads).sum() / reads if reads else 0.0
    return pd.Series(dict(reads=reads, breadth=round(breadth, 2),
                          meandepth=round(depth, 2), meanmapq=round(mapq, 1)))

vp = df.groupby(["pool", "virus"]).apply(agg_group, include_groups=False).reset_index()
vp = vp[vp.reads > 0].sort_values(["pool", "reads"], ascending=[True, False])
vp.to_csv(f"{RES}/virus_by_pool.csv", index=False)
print("\nVirus calls per pool (reads>0):")
print(vp.to_string(index=False))

# ---- Figure 1: breadth heatmap (virus x pool) ----
piv = vp.pivot_table(index="virus", columns="pool", values="breadth", fill_value=0)
# order viruses by total reads, pools by name (date-sorted)
order = vp.groupby("virus").reads.sum().sort_values(ascending=False).index
piv = piv.reindex(order)
fig, ax = plt.subplots(figsize=(max(8, 0.7*piv.shape[1]+3), 0.5*piv.shape[0]+2))
im = ax.imshow(piv.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=90, fontsize=7)
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=8)
cbar = fig.colorbar(im, ax=ax); cbar.set_label("Genome breadth of coverage (%)")
ax.set_title("Respiratory viruses across Zephyr swab pools (breadth of coverage)")
fig.tight_layout(); fig.savefig(f"{FIG}/fig1_breadth_heatmap.png", dpi=200); plt.close(fig)
print(f"\nWrote {FIG}/fig1_breadth_heatmap.png")

# ---- Figure 2: per-genome depth for the best-covered non-segmented viruses ----
def depth_array(bam, rname):
    out = subprocess.run(["samtools", "depth", "-a", "-r", rname, bam],
                         capture_output=True, text=True).stdout
    pos, dep = [], []
    for ln in out.strip().split("\n"):
        if not ln: continue
        _, p_, d_ = ln.split("\t"); pos.append(int(p_)); dep.append(int(d_))
    return np.array(pos), np.array(dep)

# choose up to 4 single-genome viruses with the highest breadth in any pool
single_df = df[df.virus.isin(SINGLE.values())]
best = (single_df.sort_values("breadth", ascending=False)
        .drop_duplicates("virus").head(4))
if len(best):
    n = len(best); fig, axes = plt.subplots(n, 1, figsize=(9, 2.2*n), squeeze=False)
    for ax, (_, r) in zip(axes[:, 0], best.iterrows()):
        pos, dep = depth_array(f"{RES}/{r.pool}.bam", r.rname)
        ax.fill_between(pos, dep, step="mid", color="#1F3A5F")
        ax.set_title(f"{r.virus}  ({r.pool}, breadth {r.breadth:.0f}%, mean depth {r.meandepth:.1f}x)",
                     fontsize=9)
        ax.set_xlabel("genome position (bp)", fontsize=8); ax.set_ylabel("depth", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig2_genome_coverage.png", dpi=200); plt.close(fig)
    print(f"Wrote {FIG}/fig2_genome_coverage.png")

print("\nPart 2 done.")
