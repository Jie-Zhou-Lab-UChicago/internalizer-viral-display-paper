# Source data (download from UCSC Xena — not committed)

Cohort "TCGA TARGET GTEx" at https://xenabrowser.net/datapages/
- Expression: `TcgaTargetGtex_rsem_gene_tpm` (log2(TPM+0.001), genes x 19,131, ~8 GB)
- Phenotype : `TcgaTargetGTEX_phenotype.txt`

Extract the 12 gene rows:
```
bash scripts/extract_genes.sh /path/to/TcgaTargetGtex_rsem_gene_tpm > data/raw/gene_subset.tsv
```
The pre-extracted per-sample table in data/extracted/ lets the figures reproduce
without this download.
