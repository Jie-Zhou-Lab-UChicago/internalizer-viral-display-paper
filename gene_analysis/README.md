# Internalizing-arm tumor selectivity (TCGA + GTEx)

Transcriptomic analysis accompanying:

> **Harnessing Endocytosis with Viral Display to Select Highly Developable,
> Internalizing Antibodies for Precision Therapy and Delivery**
> Mengbing Zou#, Baiyan Qian#, Yufan Zhu, Jenny Huo, Keyuan Ren, Rui Fu,
> Xiaokang Jin, Xiaoyu Zhang, James A. Wells, Yamuna Krishnan, Weixin Tang,
> Jie Zhou\*  (#equal contribution; \*corresponding: jiezhou6@uchicago.edu)
> Zhou Lab, University of Chicago.

## What this analysis shows

Internalizing arms for extracellular targeted protein degradation (eTPD) and
targeted delivery are conventionally built on a small set of **fast-internalizing
but broadly-expressed** receptors — TfR1 (TFRC), CXCR7 (ACKR3), CI-M6PR (IGF2R),
ASGPR (ASGR1), and related uptake receptors (LDLR, sortilin/SORT1, LRP1). This
analysis compares those conventional arms with the **candidate internalizing-arm
targets isolated in this study — CDCP1, MET, and EphA2** — for tumor selectivity
across TCGA + GTEx, with **HER2 (ERBB2)** and **EGFR** as validated ADC-target
controls.

**Result:** MET and CDCP1 are significantly more tumor-upregulated than every
conventional uptake receptor, and more than the validated ADC targets HER2/EGFR,
while the conventional effector receptors are broadly expressed in normal tissue
and flat-to-negative in tumor.

| gene | class | mean Δ (tumor−normal) | q (BH) | sig |
|---|---|---|---|---|
| MET | internalizing arm | +1.39 | 1.2e-05 | **** |
| CDCP1 | internalizing arm | +1.13 | 7.3e-04 | *** |
| TFRC (TfR1) | uptake receptor | +0.52 | 3.2e-02 | * |
| ERBB2 (HER2) | ADC control | +0.49 | 2.5e-02 | * |
| EphA2 | internalizing arm | +0.48 | 1.8e-01 | ns |
| IGF2R (CI-M6PR) | uptake receptor | +0.34 | 2.2e-02 | * |
| ACKR3 (CXCR7) | uptake receptor | +0.19 | 1.0e+00 | ns |
| SORT1 (sortilin) | uptake receptor | +0.14 | 1.0e+00 | ns |
| EGFR | ADC control | −0.05 | 1.0e+00 | ns |
| LDLR | uptake receptor | −0.19 | 1.0e+00 | ns |
| ASGR1 (ASGPR) | uptake receptor | −0.37 | 4.9e-01 | ns |
| LRP1 | uptake receptor | −0.38 | 3.1e-01 | ns |

Δ = median(tumor) − median(matched normal) per cancer, log2(TPM+1).
Significance = one-sample Wilcoxon of the 29 per-cancer Δ values vs 0, BH-FDR.

## Environment

```bash
pip install -r requirements.txt         # or:
conda env create -f environment.yml && conda activate internalizing-arm
```

## Reproduce (no download needed)

```bash
pip install -r requirements.txt
python scripts/make_figures.py        # or: make
```

Reads the included per-sample table and writes
`results/figures/InternalizingArm_tumor_upregulation_dotsbar.(png|svg)` + stats.

## Regenerate from source (optional)

See `data/raw/README.md` to download the two Xena files, then:

```bash
bash scripts/extract_genes.sh /path/to/TcgaTargetGtex_rsem_gene_tpm > data/raw/gene_subset.tsv
```

## Layout

```
scripts/
  make_figures.py            self-contained: figure + stats from extracted data
  extract_genes.sh           pull the 12 gene rows from the 8 GB Xena matrix
data/
  raw/README.md              how to obtain source matrices (not committed)
  extracted/InternalizingArm_RAW_persample_log2TPM.csv   19,131 samples x 12 genes
  derived/InternalizingArm_perCancer_delta.csv           29 cancers x 12 genes (Δ)
          InternalizingArm_median_tumor.csv / _median_normal.csv
          InternalizingArm_GTEx_tissue_median.csv        31 tissues x 12 genes
          InternalizingArm_stats_summary.csv             per-gene Wilcoxon q + stars
results/figures/             PNG (300 dpi) + editable SVG
prism/InternalizingArm_for_Prism.xlsx   GraphPad Prism tables (5 sheets, incl. stats)
```

## Methods

UCSC Xena TOIL RSEM TPM (log2(TPM+0.001)) → log2(TPM+1). Tumor = TCGA tumor
samples; Normal = TCGA adjacent-normal + matched GTEx tissue (GEPIA2/TIMER
convention). Panel of 12 genes defined by Ensembl ID in `scripts/extract_genes.sh`.
MESO/UVM excluded (no matched normal).

Primary data are public (TCGA, GTEx via UCSC Xena). Code licensed MIT; see
`CITATION.cff`.
