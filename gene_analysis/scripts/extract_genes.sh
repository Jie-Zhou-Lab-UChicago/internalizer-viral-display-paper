#!/usr/bin/env bash
# Extract the gene-panel rows from the UCSC Xena TOIL expression matrix.
# Source (download once, cohort "TCGA TARGET GTEx", https://xenabrowser.net):
#   TcgaTargetGtex_rsem_gene_tpm  (log2(TPM+0.001), genes x 19,131 samples, ~8 GB)
#   TcgaTargetGTEX_phenotype.txt
# Usage: bash scripts/extract_genes.sh /path/to/TcgaTargetGtex_rsem_gene_tpm > data/raw/gene_subset.tsv
set -euo pipefail
MATRIX="${1:?path to TcgaTargetGtex_rsem_gene_tpm}"
IDS=(
  ENSG00000163814  # CDCP1     (candidate arm)
  ENSG00000105976  # MET       (candidate arm)
  ENSG00000142627  # EPHA2     (candidate arm)
  ENSG00000141736  # ERBB2/HER2 (ADC control)
  ENSG00000146648  # EGFR       (ADC control)
  ENSG00000072274  # TFRC/TfR1  (conventional uptake receptor)
  ENSG00000144476  # ACKR3/CXCR7
  ENSG00000197081  # IGF2R/CI-M6PR
  ENSG00000130164  # LDLR
  ENSG00000134243  # SORT1/sortilin
  ENSG00000123384  # LRP1
  ENSG00000141505  # ASGR1/ASGPR
)
head -1 "$MATRIX"
for id in "${IDS[@]}"; do grep -m1 "^${id}" "$MATRIX"; done
