## =============================================================================
## Phage- vs viral-selected VHHs: developability metric comparison + Prism export
## Input : ../results/analysis.csv  (master table from R/run_master_pipeline.R)
## Output: per-metric histograms, Prism-ready bin tables, and a stats summary CSV
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(effectsize)
})

IN   <- "results/analysis.csv"
OUT  <- "results"
NBIN <- 50
METRICS <- c(psh_15_ww = "Patches of surface hydrophobicity (Wimley-White)",
             psh_15_em = "Patches of surface hydrophobicity (Eisenberg-McLachlan)",
             psh_15_kd = "Patches of surface hydrophobicity (Kyte-Doolittle, TAP)",
             ppc       = "Patches of positive charge (PPC)",
             pnc       = "Patches of negative charge (PNC)")

d <- read.csv(IN, stringsAsFactors = FALSE)
d$group <- substr(d$ID, 1, 8)                      # "Ab_phage" / "Ab_viral"

summary_rows <- list()
for (m in names(METRICS)) {
  if (!m %in% names(d)) next
  sub <- d[, c(m, "group")]; colnames(sub)[1] <- "value"
  sub <- sub[is.finite(sub$value), ]
  ph <- sub$value[sub$group == "Ab_phage"]
  vi <- sub$value[sub$group == "Ab_viral"]

  ## --- overlaid histogram ---
  p <- ggplot(sub, aes(value, fill = group, color = group)) +
    geom_histogram(alpha = 0.5, bins = NBIN, position = "identity") +
    labs(x = METRICS[[m]], y = "Count", title = paste0(m, ": phage vs viral")) +
    theme_minimal()
  ggsave(file.path(OUT, paste0("hist_", m, ".pdf")), p, width = 6, height = 4)

  ## --- common-break bin table for Prism ---
  brks <- seq(min(sub$value), max(sub$value), length.out = NBIN + 1)
  wide <- sub %>%
    mutate(bin = cut(value, breaks = brks, include.lowest = TRUE, right = FALSE)) %>%
    group_by(group, bin) %>% summarise(count = n(), .groups = "drop") %>%
    mutate(bin_mid = (as.numeric(sub("\\[([^,]*),.*", "\\1", bin)) +
                      as.numeric(sub("[^,]*,([^]]*)\\)", "\\1", bin))) / 2) %>%
    select(group, bin_mid, count) %>%
    pivot_wider(names_from = group, values_from = count, values_fill = 0)
  write.csv(wide, file.path(OUT, paste0(m, "_histogram_bins.csv")), row.names = FALSE)

  ## --- statistics (KS, Wilcoxon, Welch t, Cohen's d) ---
  ks <- ks.test(ph, vi); wx <- wilcox.test(ph, vi); tt <- t.test(ph, vi)
  cd <- effectsize::cohens_d(ph, vi)
  summary_rows[[m]] <- data.frame(
    metric = m, n_phage = length(ph), n_viral = length(vi),
    mean_phage = round(mean(ph), 3), mean_viral = round(mean(vi), 3),
    KS_D = round(unname(ks$statistic), 3),  KS_p = signif(ks$p.value, 3),
    Wilcoxon_W = unname(wx$statistic),      Wilcoxon_p = signif(wx$p.value, 3),
    t = round(unname(tt$statistic), 3),     t_p = signif(tt$p.value, 3),
    cohens_d = round(cd$Cohens_d, 3),
    d_CI_low = round(cd$CI_low, 3), d_CI_high = round(cd$CI_high, 3))
}
stats_tbl <- do.call(rbind, summary_rows)
write.csv(stats_tbl, file.path(OUT, "phage_vs_viral_stats_summary.csv"), row.names = FALSE)
print(stats_tbl)
