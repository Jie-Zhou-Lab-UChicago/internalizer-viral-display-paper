.PHONY: all figures extract
all: figures
figures:
	python scripts/make_figures.py
extract:
	bash scripts/extract_genes.sh $(MATRIX) > data/raw/gene_subset.tsv
