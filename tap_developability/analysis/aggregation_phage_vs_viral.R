## =============================================================================
## Aggrescan3D: phage- vs viral-selected VHH comparison + figures + Prism export
## Input : ../results/Aggregation3D_analysis.csv  (from R/aggrescan3d_pipeline.R)
## Output: per-metric annotated histogram PDFs, a combined multi-page PDF,
##         Prism-ready bin tables, and a stats summary (KS / Wilcoxon / t / d).
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(effectsize)
})

IN   <- "results/Aggregation3D_analysis.csv"
OUT  <- "results"
NBIN <- 50
LAB  <- c(total_agg_score = "Total A3D aggregation score",
          mean_agg_score  = "Mean per-residue A3D score",
          cdr1_agg_score  = "CDR-H1 A3D score",
          cdr2_agg_score  = "CDR-H2 A3D score",
          cdr3_agg_score  = "CDR-H3 A3D score",
          vicinity_agg_score = "CDR-vicinity A3D score")
pal  <- c(Ab_phage = "#4c78a8", Ab_viral = "#d1495b")

d <- read.csv(IN, stringsAsFactors = FALSE)
d$group <- substr(d$ID, 1, 8)                    # "Ab_phage" / "Ab_viral"
metrics <- intersect(names(LAB), names(d))
stat_rows <- list(); plots <- list()

for (m in metrics) {
  sub <- data.frame(value = d[[m]], group = d$group)
  sub <- sub[is.finite(sub$value), ]
  ph <- sub$value[sub$group == "Ab_phage"]; vi <- sub$value[sub$group == "Ab_viral"]

  ks <- ks.test(ph, vi); wx <- wilcox.test(ph, vi); tt <- t.test(ph, vi)
  cd <- effectsize::cohens_d(ph, vi)
  stat_rows[[m]] <- data.frame(metric = m, mean_phage = round(mean(ph),3), mean_viral = round(mean(vi),3),
    KS_D = round(unname(ks$statistic),3), KS_p = signif(ks$p.value,3),
    Wilcoxon_p = signif(wx$p.value,3), t_p = signif(tt$p.value,3),
    cohens_d = round(cd$Cohens_d,3), d_low = round(cd$CI_low,3), d_high = round(cd$CI_high,3))

  subtitle <- sprintf("phage n=%d, viral n=%d  |  Wilcoxon p=%s  |  Cohen's d=%.2f",
                      length(ph), length(vi), signif(wx$p.value, 3), cd$Cohens_d)
  p <- ggplot(sub, aes(value, fill = group, color = group)) +
    geom_histogram(alpha = 0.5, bins = NBIN, position = "identity") +
    scale_fill_manual(values = pal) + scale_color_manual(values = pal) +
    labs(x = LAB[[m]], y = "Count", title = paste0(LAB[[m]], ": phage vs viral"),
         subtitle = subtitle, fill = NULL, color = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  ggsave(file.path(OUT, paste0("hist_", m, ".pdf")), p, width = 6, height = 4.2)
  plots[[m]] <- p

  brks <- seq(min(sub$value), max(sub$value), length.out = NBIN + 1)
  wide <- sub %>%
    mutate(bin = cut(value, breaks = brks, include.lowest = TRUE, right = FALSE)) %>%
    group_by(group, bin) %>% summarise(count = n(), .groups = "drop") %>%
    mutate(bin_mid = (as.numeric(sub("\\[([^,]*),.*", "\\1", bin)) +
                      as.numeric(sub("[^,]*,([^]]*)\\)", "\\1", bin)))/2) %>%
    select(group, bin_mid, count) %>%
    pivot_wider(names_from = group, values_from = count, values_fill = 0)
  write.csv(wide, file.path(OUT, paste0(m, "_histogram_bins.csv")), row.names = FALSE)
}

pdf(file.path(OUT, "aggregation_histograms_all.pdf"), width = 6, height = 4.2)
for (m in metrics) print(plots[[m]]); dev.off()

stats_tbl <- do.call(rbind, stat_rows)
write.csv(stats_tbl, file.path(OUT, "aggregation_phage_vs_viral_stats.csv"), row.names = FALSE)
print(stats_tbl)

## ---- optional: merge with PSH master table for a PSH-vs-A3D scatter ----------
## psh <- read.csv("../results/analysis.csv")[, c("ID","psh_15_ww","psh_15_kd")]
## x <- merge(psh, d[, c("ID","total_agg_score","mean_agg_score")], by = "ID", all = TRUE)
## write.csv(x, file.path(OUT, "psh_vs_agg.csv"), row.names = FALSE)
