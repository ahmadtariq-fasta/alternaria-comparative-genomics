# R/09_host_enrichment.R
#
# Fisher's exact test on presence/absence of each feature: each host group vs the
# reference host. Benjamini-Hochberg FDR across all tests; candidates are raw p < 0.01 and OR > 5.
#
# Input:  samples.tsv, <tables>/{orthogroups,cazyme_families,bgc_classes,effector_classes}.tsv
# Output: <out>/all_tests.tsv, <out>/candidates.tsv
#
# Run from the repo root:
#   Rscript R/09_host_enrichment.R config/samples.tsv results/tables results/enrichment wild_tomato

args <- commandArgs(TRUE)
samples <- read.delim(args[1])
tables <- args[2]
out <- args[3]
reference <- args[4]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (!reference %in% samples$host) stop("Reference host not found in samples: ", reference)
ref_samples <- samples$sample[samples$host == reference]
focal_hosts <- setdiff(unique(samples$host), reference)

results <- list()
for (type in c("orthogroups", "cazyme_families", "bgc_classes", "effector_classes")) {
  d <- read.delim(file.path(tables, paste0(type, ".tsv")), row.names = 1, check.names = FALSE)
  if (nrow(d) == 0) { print(paste("No features in", type, "- skipped")); next }
  present <- d > 0

  for (host in focal_hosts) {
    host_samples <- samples$sample[samples$host == host]
    a <- rowSums(present[, host_samples, drop = FALSE])   # present in focal host
    b <- length(host_samples) - a                         # absent in focal host
    c <- rowSums(present[, ref_samples, drop = FALSE])    # present in reference
    e <- length(ref_samples) - c                          # absent in reference

    p <- numeric(nrow(d)); or <- numeric(nrow(d))
    for (i in seq_len(nrow(d))) {
      test <- fisher.test(matrix(c(a[i], c[i], b[i], e[i]), nrow = 2))
      p[i] <- test$p.value
      or[i] <- test$estimate
    }
    results[[paste(type, host)]] <- data.frame(
      feature_type = type, feature = rownames(d), host = host,
      present_host = a, n_host = length(host_samples),
      present_ref = c, n_ref = length(ref_samples),
      odds_ratio = or, p_value = p)
  }
}

all_tests <- do.call(rbind, results)
all_tests$fdr <- p.adjust(all_tests$p_value, method = "BH")
all_tests <- all_tests[order(all_tests$p_value), ]
candidates <- all_tests[all_tests$p_value < 0.01 & all_tests$odds_ratio > 5, ]

write.table(all_tests, file.path(out, "all_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(candidates, file.path(out, "candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sprintf("%d tests, %d pass FDR < 0.05, %d candidates (p < 0.01, OR > 5)",
              nrow(all_tests), sum(all_tests$fdr < 0.05), nrow(candidates)))
