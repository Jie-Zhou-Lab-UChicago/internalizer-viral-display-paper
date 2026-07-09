# Developability profiling of phage- vs viral-selected VHHs

Two structure-based developability analyses on predicted VHH models, comparing
antibodies enriched by phage vs lentiviral (viral) selection:

1. **TAP-style surface patches** — PSH (hydrophobicity), PPC / PNC (positive / negative charge),
   following the Therapeutic Antibody Profiler (Raybould et al., PNAS 2019).
2. **Aggregation propensity** — Aggrescan3D (A3D) total / mean / per-CDR / CDR-vicinity scores
   (Zambrano et al. 2015; Kuriata et al. 2019), run in static mode (FoldX not used).

## Inputs (provide locally; not committed)

- `data/VHH_to_analyze.csv` — 2,000 VHHs (1,000 viral + 1,000 phage): columns `ID`, `seq`. *(included)*
- `data/VHH_model_str/` — one sub-folder per ID with a predicted `*.pdb` + FreeSASA `*.rsa`. → PSH/PPC/PNC
- `data/VHH_AD3/` — one sub-folder per ID with the A3D `*.csv` + A3D `*input.pdb`. → aggregation
  (see the README inside each folder for the exact expected contents).

## Run

Open R with the working directory set to **this `tap_developability/` folder**, then:

```r
setwd("/path/to/.../tap_developability")     # this folder
install.packages("renv"); renv::restore()    # first time only: pinned package versions
source("run_all.R")                          # runs everything -> results/
```

`run_all.R` runs the four steps in order (all paths are relative to this folder,
so no need to change directories):

```r
source("R/run_master_pipeline.R")               # PSH/PPC/PNC  -> results/analysis.csv
source("R/aggrescan3d_pipeline.R")              # Aggrescan3D  -> results/Aggregation3D_analysis.csv
source("analysis/phage_vs_viral_stats.R")       # PSH/PPC/PNC stats + histograms + Prism bins
source("analysis/aggregation_phage_vs_viral.R") # aggregation stats + histograms + Prism bins
```

VHHs that have no structure / A3D sub-folder are skipped automatically (a note is printed),
so you can run on a subset without editing anything.

## Outputs (in `results/`)

- `analysis.csv` — master PSH/PPC/PNC table (WW / EM / Kyte-Doolittle PSH).
- `Aggregation3D_analysis.csv` — total / mean / CDR1-3 / vicinity A3D scores.
- `*_histogram_bins.csv` — GraphPad Prism bin tables (one per metric, common breaks).
- `hist_<metric>.pdf` + `aggregation_histograms_all.pdf` — annotated phage-vs-viral histograms.
- `*_stats*.csv` — KS, Wilcoxon, Welch t, Cohen's d per metric.

## Metric definitions & differences from reference TAP

PSH/PPC/PNC: per-residue Kyte-Doolittle hydrophobicity normalized to [1,2] (or formal charge,
His = +0.1); contribution `H(R1)*H(R2)/r12^2` (or `|Q1||Q2|/r12^2`) over residue pairs with
closest heavy-atom distance < 7.5 A; salt-bridged charges (N-O <= 3.2 A) neutralized; scored over
surface-exposed CDR residues (FreeSASA rel. side-chain SASA > 15%). Differences from reference
TAP: single-domain VHH; CDRs from the fixed IGHV3-23 scaffold numbering; FreeSASA exposure; scored
over exposed CDR residues (set `SCORE_SET="vicinity"` for the full vicinity). Absolute values are
self-consistent for ranking this library but not directly comparable to TAP clinical thresholds.

## Reproducibility

Exact package versions are pinned in `renv.lock` (R 4.4.1). To restore them:

```r
install.packages("renv"); renv::restore()
```

## Requirements

R: `bio3d`, `stringr`, `dplyr`, `tidyr`, `ggplot2`, `effectsize`. External: FreeSASA, Aggrescan3D.
