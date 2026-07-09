# VHH_AD3  (not committed — provide before running)

Aggrescan3D (A3D) results, one sub-folder per VHH named exactly as its ID (e.g. `Ab_viral_1/`),
each containing:
- the A3D per-residue score table `*.csv` (columns include `residue`, `score`), and
- the A3D structure `*input.pdb` / `*folded.pdb` (used for the CDR-vicinity score).

Used by `R/aggrescan3d_analysis.R`. A3D run in static mode (FoldX not used).
