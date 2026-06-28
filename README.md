# Zephyr respiratory-virus analysis

Taxonomic classification, genome coverage, and embedding-based clustering of twelve
Zephyr Oxford Nanopore pools.

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
refs/respiratory_panel.fasta   37 RefSeq genomes
results/                    coverage / consensus / cluster / ANI tables
figures/                    fig1..fig5
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
python part3_step6_skani.py
```

Raw pools (~149 MB) and BAM files are not committed (see `.gitignore`); the scripts
expect the pools under `data/`.

## Results at a glance

- Taxonomy: ten viruses confidently detected across the twelve pools by breadth of
  coverage. Influenza A only partial, due to an out-of-date reference.
- Clustering baseline (sourmash MinHash): ARI 0.88, NMI 0.96.
- Learned tier (ESM2-8M on the replicase marker): ARI 0.52; recovers family/genus
  structure but is coarser at species level. Scaling ESM2 from 8M to 650M does not
  change this at all (flat ARI), so the limit is the marker, not model capacity.
- skani ANI: coronaviruses confirm same-species (96-99.9%); within-species
  rhinovirus pairs fall below the 80% screen, i.e. they are different serotypes.
  One reference per rhinovirus species is too coarse for this data.
