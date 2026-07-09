## Run the whole developability analysis in order.
## 1) In R, set the working directory to THIS folder (tap_developability):
##       setwd("/path/to/.../tap_developability")
## 2) First time only, restore pinned package versions:
##       install.packages("renv"); renv::restore()
## 3) Then:
##       source("run_all.R")
source("R/run_master_pipeline.R")          # PSH/PPC/PNC  -> results/analysis.csv
source("R/aggrescan3d_pipeline.R")         # Aggrescan3D  -> results/Aggregation3D_analysis.csv
source("analysis/phage_vs_viral_stats.R")  # PSH/PPC/PNC stats + histograms + Prism bins
source("analysis/aggregation_phage_vs_viral.R")  # aggregation stats + histograms + Prism bins
message("All done. Outputs are in results/.")
