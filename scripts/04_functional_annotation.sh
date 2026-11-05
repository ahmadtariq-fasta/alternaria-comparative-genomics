#!/bin/bash
# scripts/04_functional_annotation.sh
#
# Orthology-based (eggNOG-mapper) and domain-based (InterProScan) annotation.
#
# Input:  $RESULTS/proteins/<sample>.faa
# Output: $RESULTS/<sample>/eggnog/, $RESULTS/<sample>/interproscan/
#
# Run from the repo root: sbatch scripts/04_functional_annotation.sh
#SBATCH --job-name=functional_annotation
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

for sample in $(tail -n +2 "$SAMPLES" | cut -f1); do
  out="$RESULTS/$sample"
  # InterProScan rejects the stop-codon asterisks in predicted proteins
  mkdir -p "$out/eggnog" "$out/interproscan"
  sed 's/\*$//' "$RESULTS/proteins/$sample.faa" > "$out/interproscan/input.faa"

  emapper.py -i "$RESULTS/proteins/$sample.faa" -o "$sample" --output_dir "$out/eggnog" \
    --data_dir "$EMAPPER_DB" --cpu "$THREADS" --override
  interproscan.sh -i "$out/interproscan/input.faa" -d "$out/interproscan" -f tsv -goterms -cpu "$THREADS"
done
