# R/make_dummy_data.R
#
# Writes completely simulated tool outputs for 24 made-up isolates, so the analysis
# steps can be run without any real data. Nothing here comes from a real genome.
# Planted signals: DUMMY_OG00001-20 in potato, 21-40 in tomato, 41-60 in brassicaceae;
# CAZyme family PL1 in potato; BGC class terpene in tomato.
#
# Output: data/dummy/samples.tsv, data/dummy/Orthogroups.GeneCount.tsv,
#         data/dummy/results/<sample>/{dbcan,effectorp,antismash}/
#
# Run from the repo root: Rscript R/make_dummy_data.R

set.seed(42)
hosts <- rep(c("potato", "tomato", "brassicaceae", "wild_tomato"), each = 6)
ids <- sprintf("DUMMY_%02d", seq_along(hosts))
root <- "data/dummy"
unlink(file.path(root, "results"), recursive = TRUE)
dir.create(root, recursive = TRUE, showWarnings = FALSE)

write.table(data.frame(sample = ids, host = hosts, r1 = "none", r2 = "none"),
            file.path(root, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Orthogroups: 60 planted, the rest random with a per-gene-family presence rate
n_og <- 500
og <- matrix(0L, n_og, length(ids), dimnames = list(sprintf("DUMMY_OG%05d", 1:n_og), ids))
for (i in 61:n_og) og[i, ] <- rbinom(length(ids), 1, runif(1, 0.3, 1)) * rpois(length(ids), 1.5)
og[1:20, hosts == "potato"] <- 1L
og[21:40, hosts == "tomato"] <- 1L
og[41:60, hosts == "brassicaceae"] <- 1L
write.table(data.frame(Orthogroup = rownames(og), og, Total = rowSums(og), check.names = FALSE),
            file.path(root, "Orthogroups.GeneCount.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (k in seq_along(ids)) {
  s <- ids[k]; host <- hosts[k]
  base <- file.path(root, "results", s)
  for (d in c("dbcan", "effectorp", "antismash")) dir.create(file.path(base, d), recursive = TRUE)

  # dbCAN overview: last row of each pair is seen by one tool only, so must not count
  fams <- c("GH5_7", "GH28", "AA9", "CE1", "GT2", if (host == "potato") "PL1")
  rows <- c("Gene ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools")
  for (f in fams) rows <- c(rows, sprintf("gene_%s\t-\t%s(10-200)\t%s_e1\t%s\t3", f, f, f, f))
  rows <- c(rows, "gene_single\t-\tCBM1(5-40)\t-\t-\t1", "gene_none\t-\t-\t-\t-\t0")
  writeLines(rows, file.path(base, "dbcan", "overview.txt"))

  # EffectorP table
  n_e <- rpois(1, if (host == "wild_tomato") 3 else 6) + 1
  pred <- sample(c("Cytoplasmic effector", "Apoplastic effector", "Apoplastic/cytoplasmic effector", "Non-effector"),
                 n_e, TRUE)
  writeLines(c("# Identifier\tCytoplasmic effector\tApoplastic effector\tNon-effector\tPrediction",
               sprintf("prot_%d\t0.1\t0.8\t-\t%s", seq_len(n_e), pred)),
             file.path(base, "effectorp", "effectorp.txt"))

  # antiSMASH region files; the CDS product must not be read as a BGC class
  bgcs <- c("T1PKS", "NRPS", if (k %% 2 == 0) "T1PKS", if (host == "tomato") "terpene")
  for (j in seq_along(bgcs)) {
    writeLines(c("LOCUS       dummy_contig", "FEATURES             Location/Qualifiers",
                 "     region          1..1000", sprintf('                     /product="%s"', bgcs[j]),
                 "     CDS             10..200", '                     /product="hypothetical protein"',
                 "ORIGIN", "//"),
               file.path(base, "antismash", sprintf("contig.region%03d.gbk", j)))
  }
}
print("Done: dummy data in data/dummy")
