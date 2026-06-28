# Zephyr respiratory-virus analysis

Taxonomic classification, genome coverage, and embedding-based clustering of 12
Zephyr Oxford Nanopore pools.

## Data

There are two kinds of input, kept separate.

**Sample pools (the 12 Zephyr pools).** These are the Oxford Nanopore metagenomic
read sets we analysed, from the Zephyr respiratory-virus dataset in the work-test
brief. They are ~149 MB and are not committed; the scripts read them from `data/`.
The 12 pools span 4 sites (BC, BoDT, Copl, MBTA_Ha) and January 2025 to June 2026:

```
250129-BoDT-NAS        260109-BoDT-NAS-P1     260401-BC-NAS-P1
250812-MBTA_Ha-NAS     260115-Copl-NAS-P1     260415-BC-NAS-P1
251217-BoDT-NAS-P1     260304-BC-NAS-P1       260514-BoDT-NAS-P1
251219-MBTA_Ha-NAS-P1  260313-MBTA_Ha-NAS-P1  260613-MBTA_Ha-NAS-P1
```

Naming is `YYMMDD-site-NAS[-P1]`. Even without the raw pools, the per-pool
detections are committed in `results/virus_by_pool.csv`, so the findings are
inspectable directly. Source for the pools: the Zephyr respiratory viral reads at
https://data.securebio.org/zephyr/#respiratory-viral-reads (full read archive on
the SRA, BioProject PRJNA1379685).

**Reference panels (committed, in `refs/`).** `respiratory_panel.fasta` is the
panel of 37 RefSeq genomes used for the Part 2 taxonomy. `rhino_refs.fasta` is the
48 reference rhinovirus genomes (16 each of A, B, C) used for the serotype
demonstration. These are reference databases, not the sample pools.

## Layout

```
part2_taxonomy.py           competitive minimap2 mapping + breadth/depth coverage  -> fig1, fig2
part3_step0_select.py       select detections clean enough to reconstruct
part3_step1_consensus.py    build consensus genomes (samtools consensus)
part3_step2_sourmash.py     nucleotide MinHash baseline clustering                 -> fig3
part3_step3_longorf.py      longest-ORF replicase-proxy marker per genome
part3_step4_esm_cluster.py  ESM2-8M protein embedding + clustering                 -> fig4
part3_step5_esm_sweep.py    ESM2 model-size sweep (8M -> 650M)                     -> fig6
part3_step6_skani.py        skani ANI within-cluster validity check                -> fig5
part3_step7_robustness.py   permutation null + leave-one-out on the ARI            -> fig7
part3_fetch_rhino_refs.py   fetch 48 reference rhinovirus genomes from NCBI
part3_step8_rhino_serotypes.py  place our rhinoviruses among reference types       -> fig8
refs/respiratory_panel.fasta   37 RefSeq genomes (Part 2 taxonomy panel)
refs/rhino_refs.fasta          48 reference rhinovirus genomes (serotype demo)
data/                       the 12 Zephyr pools (not committed; see Data above)
results/                    coverage / consensus / cluster / ANI tables
figures/                    fig1..fig8
```

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# external tools: minimap2 2.31, samtools 1.23.1, skani 0.3.2 (brew install skani)

python part2_taxonomy.py
python part3_step0_select.py
python part3_step1_consensus.py
python part3_step2_sourmash.py
python part3_step3_longorf.py
python part3_step4_esm_cluster.py
python part3_step5_esm_sweep.py        # downloads ESM2 weights up to 650M
python part3_step6_skani.py
python part3_step7_robustness.py
python part3_fetch_rhino_refs.py       # needs internet (NCBI Entrez)
python part3_step8_rhino_serotypes.py
```

The raw pools and BAM files are not committed (see `.gitignore` and the Data
section); the scripts expect the pools under `data/`.

## Results at a glance

- Taxonomy: 10 viruses confidently detected across the 12 pools by breadth of
  coverage. Influenza A only partial, due to an out-of-date reference.
- Clustering baseline (sourmash MinHash): ARI 0.88, NMI 0.96.
- Learned tier (ESM2-8M on the replicase marker): ARI 0.52; recovers family/genus
  structure but is coarser at species level. Scaling ESM2 from 8M to 650M does not
  change this at all (flat ARI), so the limit is the marker, not model capacity.
- skani ANI: coronaviruses confirm same-species (96-99.9%); within-species
  rhinovirus pairs fall below the 80% screen, i.e. they are different serotypes.
  A single reference per rhinovirus species is too coarse for this data.
- Robustness: permutation p = 0.0005 (sketch) and 0.0015 (embedding); leave-one-out
  ARI stays 0.83-1.00 (sketch) and 0.39-0.65 (embedding). The agreement is not a
  small-n fluke.
- Serotype demonstration: against 48 reference rhinovirus genomes, each of our
  same-species genomes matches a different reference type, confirming the serotype
  diversity directly.
