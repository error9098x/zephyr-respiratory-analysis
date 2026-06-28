#!/usr/bin/env python3
"""
Track B, Step 0 - select which (pool, virus) detections are clean enough to
reconstruct a consensus genome from.

Selection thresholds (tunable): breadth of coverage, mean depth, read count.
We also flag low mean mapping quality, which marks ambiguous cross-mapping
(the rhinovirus-species issue from Part 2) where a consensus would be unreliable.
"""
import pandas as pd

RES = "/Users/aviral/Desktop/Projects/MATS/zephyr_worktest/results"

BREADTH_MIN = 80.0   # % of genome covered - decisive "really present" signal
DEPTH_MIN   = 20.0   # mean depth - supports consensus base calls
READS_MIN   = 200    # minimum reads
MAPQ_FLAG   = 30.0   # below this we flag as ambiguous mapping

df = pd.read_csv(f"{RES}/virus_by_pool.csv")
df["keep"]     = (df.breadth >= BREADTH_MIN) & (df.meandepth >= DEPTH_MIN) & (df.reads >= READS_MIN)
df["low_mapq"] = df.meanmapq < MAPQ_FLAG

kept = df[df.keep].copy().sort_values(["virus", "pool"])
dropped = df[~df.keep].copy()

print(f"Thresholds: breadth >= {BREADTH_MIN}%, depth >= {DEPTH_MIN}x, reads >= {READS_MIN}\n")
print("=== KEPT detections (consensus candidates) ===")
print(kept[["pool", "virus", "reads", "breadth", "meandepth", "meanmapq", "low_mapq"]].to_string(index=False))

n_genomes = len(kept)
n_viruses = kept.virus.nunique()
multi = kept.groupby("virus").pool.nunique()
multi_pool = multi[multi >= 2]
print("\n=== SELECTION SUMMARY ===")
print(f"genomes kept            : {n_genomes}")
print(f"distinct viruses        : {n_viruses}")
print(f"viruses in >=2 pools    : {len(multi_pool)}  -> {list(multi_pool.index)}")
print(f"flagged low-MAPQ in kept: {int(kept.low_mapq.sum())}")

print("\n=== viruses present but DROPPED (below thresholds) ===")
print(dropped.groupby("virus").agg(detections=("pool", "nunique"),
      max_breadth=("breadth", "max"), max_reads=("reads", "max")).round(1).to_string())

kept.to_csv(f"{RES}/selected_detections.csv", index=False)
print(f"\nWrote {RES}/selected_detections.csv")
