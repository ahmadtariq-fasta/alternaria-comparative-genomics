# tests/check_demo.R
#
# Fails if the demo no longer recovers the signals planted by R/make_dummy_data.R.
#
# Run from the repo root (after run_demo.sh): Rscript tests/check_demo.R

cand <- read.delim("demo_output/enrichment/candidates.tsv")
caz <- read.delim("demo_output/tables/cazyme_families.tsv", row.names = 1)
bgc <- read.delim("demo_output/tables/bgc_classes.tsv", row.names = 1)
eff <- read.delim("demo_output/tables/effector_classes.tsv", row.names = 1)

found <- function(host, features) all(features %in% cand$feature[cand$host == host])
planted <- sprintf("DUMMY_OG%05d", 1:60)
stopifnot(
  found("potato", planted[1:20]),
  found("tomato", planted[21:40]),
  found("brassicaceae", planted[41:60]),
  found("potato", "PL1"),
  found("tomato", "terpene"),
  !"CBM1" %in% rownames(caz),                 # seen by one dbCAN tool only
  all(caz["GH5", ] == 1),                     # subfamily collapsed to family
  !"hypothetical protein" %in% rownames(bgc), # CDS product is not a BGC class
  setequal(rownames(eff), c("Apoplastic", "Cytoplasmic"))  # Non-effector excluded
)
print("All checks passed")
