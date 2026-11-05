# R/10_plot_counts.R
#
# Boxplots of CAZyme, effector and BGC counts per isolate, by host group.
#
# Input:  <tables>/counts_per_isolate.tsv
# Output: <out>/counts_by_host.png
#
# Run from the repo root:
#   Rscript R/10_plot_counts.R results/tables/counts_per_isolate.tsv results/enrichment "Alternaria isolates"

args <- commandArgs(TRUE)
counts <- read.delim(args[1])
out <- args[2]
title <- args[3]

png(file.path(out, "counts_by_host.png"), width = 2400, height = 800, res = 200)
par(mfrow = c(1, 3), mar = c(7, 4, 3, 1), oma = c(0, 0, 2, 0))
for (feature in c("CAZymes", "Effectors", "BGCs")) {
  boxplot(counts[[feature]] ~ counts$host, main = feature, xlab = "", ylab = "Count per isolate", las = 2)
  stripchart(counts[[feature]] ~ counts$host, vertical = TRUE, add = TRUE, pch = 16, method = "jitter")
}
mtext(title, outer = TRUE, font = 2)
dev.off()
