# R/08_build_feature_tables.R
#
# Turns the tool outputs into one feature-by-isolate count table per feature type.
#
# Input:  samples.tsv, <results>/<sample>/dbcan/overview.txt, effectorp/effectorp.txt,
#         antismash/*region*.gbk, <orthogroups file> (OrthoFinder Orthogroups.GeneCount.tsv)
# Output: <out>/orthogroups.tsv, cazyme_families.tsv, bgc_classes.tsv, effector_classes.tsv,
#         counts_per_isolate.tsv
#
# Run from the repo root:
#   Rscript R/08_build_feature_tables.R config/samples.tsv results \
#     results/proteins/OrthoFinder/Results_alternaria/Orthogroups/Orthogroups.GeneCount.tsv results/tables

args <- commandArgs(TRUE)
samples <- read.delim(args[1])
results <- args[2]
orthogroup_file <- args[3]
out <- args[4]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Makes a features x isolates count table from a list of feature names per isolate
count_table <- function(per_sample) {
  features <- sort(unique(unlist(per_sample)))
  if (length(features) == 0) return(matrix(0L, 0, length(per_sample), dimnames = list(NULL, names(per_sample))))
  tab <- sapply(per_sample, function(x) as.integer(table(factor(x, levels = features))))
  matrix(tab, nrow = length(features), dimnames = list(features, names(per_sample)))
}

# Families one gene gets from one dbCAN tool, e.g. "GH5_7(30-300)+CBM1(320-350)" -> GH5, CBM1
tool_families <- function(cell) {
  if (cell %in% c("-", "N", "")) return(character(0))
  fam <- strsplit(cell, "+", fixed = TRUE)[[1]]
  unique(sub("_.*", "", sub("\\(.*", "", fam)))
}

# CAZyme families: a family counts for a gene if at least 2 of the 3 dbCAN tools report it
cazyme_families <- function(file) {
  o <- read.delim(file, check.names = FALSE, colClasses = "character")
  keep <- c()
  for (i in seq_len(nrow(o))) {
    calls <- c(tool_families(o$HMMER[i]), tool_families(o$dbCAN_sub[i]), tool_families(o$DIAMOND[i]))
    calls <- table(calls)
    keep <- c(keep, names(calls)[calls >= 2])
  }
  keep
}

# BGC classes: the product(s) on the "region" feature of each antiSMASH region file
bgc_classes <- function(folder) {
  files <- list.files(folder, pattern = "region[0-9]+\\.gbk$", full.names = TRUE)
  classes <- c()
  for (f in files) {
    lines <- readLines(f)
    start <- grep("^     region +", lines)[1]
    ends <- grep("^     [a-zA-Z_]+ +[0-9<]", lines)
    end <- c(ends[ends > start], grep("^ORIGIN", lines), length(lines) + 1)[1]
    block <- lines[start:(end - 1)]
    products <- sub('.*/product="(.*)".*', "\\1", grep("/product=", block, value = TRUE))
    classes <- c(classes, paste(sort(products), collapse = "+"))
  }
  classes
}

# Effector classes: EffectorP's prediction, using the first-named class for dual predictions
effector_classes <- function(file) {
  lines <- grep("^#|^\\s*$", readLines(file), value = TRUE, invert = TRUE)
  prediction <- sapply(strsplit(lines, "\t"), function(x) x[length(x)])
  prediction <- prediction[grepl("effector", prediction) & !grepl("^Non", prediction)]
  sub("[/ ].*", "", prediction)
}

caz <- list(); bgc <- list(); eff <- list()
for (s in samples$sample) {
  caz[[s]] <- cazyme_families(file.path(results, s, "dbcan", "overview.txt"))
  bgc[[s]] <- bgc_classes(file.path(results, s, "antismash"))
  eff[[s]] <- effector_classes(file.path(results, s, "effectorp", "effectorp.txt"))
}

og <- read.delim(orthogroup_file, row.names = 1, check.names = FALSE)
og <- as.matrix(og[, samples$sample])

write_table <- function(x, name) {
  write.table(data.frame(feature = as.character(rownames(x)), x, check.names = FALSE),
              file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
write_table(og, "orthogroups.tsv")
write_table(count_table(caz), "cazyme_families.tsv")
write_table(count_table(bgc), "bgc_classes.tsv")
write_table(count_table(eff), "effector_classes.tsv")

counts <- data.frame(sample = samples$sample, host = samples$host,
                     CAZymes = lengths(caz), Effectors = lengths(eff), BGCs = lengths(bgc))
write.table(counts, file.path(out, "counts_per_isolate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print("Done: feature tables written")
