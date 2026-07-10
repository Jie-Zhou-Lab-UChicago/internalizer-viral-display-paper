# Harnessing Endocytosis with Viral Display to Select Highly Developable, Internalizing Antibodies

Analysis code and data for:

> **Harnessing Endocytosis with Viral Display to Select Highly Developable,
> Internalizing Antibodies for Precision Therapy and Delivery**


This repository contains two independent analysis modules, one per figure panel.

| Folder | Manuscript figure | What it does | Language |
|---|---|---|---|
| [`tap_developability/`](tap_developability/) | Fig. 4 (developability) | TAP-style developability profiling (PSH, PPC, PNC) and aggregation score of phage- vs viral-selected VHHs from predicted structures | R | 
| [`gene_analysis/`](gene_analysis/) |  Suppl. (target selectivity) | Tumor selectivity of candidate internalizing arms (CDCP1/MET/EphA2) vs conventional eTPD/uptake receptors across TCGA + GTEx | Python |

Each module is self-contained with its own README, code, small derived data, and results.
Large inputs (predicted PDB/RSA structures; the ~8 GB Xena expression matrix) are **not**
committed — see each module's README for how to obtain them; the repo ships the extracted /
derived tables needed to reproduce the figures.

## Quick start

```bash
# TAP developability (R):
cd tap_developability && Rscript R/run_master_pipeline.R        # -> results/analysis.csv
Rscript analysis/phage_vs_viral_stats.R                          # -> stats + Prism bin tables

# Target selectivity (Python):
cd gene_analysis && pip install -r requirements.txt && python scripts/make_figures.py
```

## Citation & license

Please cite the manuscript (see `CITATION.cff`). Code released under the MIT License.
Primary public data: TCGA and GTEx via UCSC Xena. Predicted VHH structures generated in-house.
