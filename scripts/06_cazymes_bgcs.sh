#!/bin/bash
# scripts/06_cazymes_bgcs.sh
#
# Annotate CAZymes with dbCAN3 and biosynthetic gene clusters with antiSMASH.
#
# Input:  $RESULTS/proteins/<sample>.faa, $RESULTS/<sample>/funannotate/predict_results/*.gbk
# Output: $RESULTS/<sample>/dbcan/overview.txt, $RESULTS/<sample>/antismash/
#
# Run from the repo root: sbatch scripts/06_cazymes_bgcs.sh
#SBATCH --job-name=cazymes_bgcs
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

for sample in $(tail -n +2 "$SAMPLES" | cut -f1); do
  out="$RESULTS/$sample"
  use_env dbcan
  run_dbcan "$RESULTS/proteins/$sample.faa" protein --out_dir "$out/dbcan" --db_dir "$DBCAN_DB" --tools all

  use_env antismash
  antismash --taxon fungi --cb-general --asf --pfam2go --genefinding-tool none \
    --cpus "$THREADS" --output-dir "$out/antismash" \
    "$out"/funannotate/predict_results/*.gbk
done
