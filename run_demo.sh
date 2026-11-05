#!/bin/bash
# run_demo.sh
#
# Runs the analysis steps (08-10) on simulated data and checks the planted signals are found.
#
# Output: demo_output/
#
# Run from the repo root: bash run_demo.sh

set -euo pipefail
Rscript R/make_dummy_data.R
Rscript R/08_build_feature_tables.R data/dummy/samples.tsv data/dummy/results \
  data/dummy/Orthogroups.GeneCount.tsv demo_output/tables
Rscript R/09_host_enrichment.R data/dummy/samples.tsv demo_output/tables demo_output/enrichment wild_tomato
Rscript R/10_plot_counts.R demo_output/tables/counts_per_isolate.tsv demo_output/enrichment "DUMMY DATA: not real results"
Rscript tests/check_demo.R
