# VHH_model_str  (not committed — provide before running)

One sub-folder per VHH, named exactly as its ID (e.g. `Ab_viral_1/`), each containing:
- the predicted VHH structure `*.pdb` (Boltz-2 / chosen predictor), and
- a FreeSASA relative-accessibility file `*.rsa`.

Used by `R/run_master_pipeline.R` (PSH / PPC / PNC). Numbering must be sequential 1..N
(resno == sequence position).
