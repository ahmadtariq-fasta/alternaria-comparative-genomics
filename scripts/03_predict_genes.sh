#!/bin/bash
# scripts/03_predict_genes.sh
#
# Clean, sort and mask each assembly, then predict genes with Funannotate.
#
# Input:  $RESULTS/<sample>/spades/scaffolds.fasta
# Output: $RESULTS/<sample>/funannotate/predict_results/, $RESULTS/proteins/<sample>.faa
#
# Run from the repo root: sbatch scripts/03_predict_genes.sh
#SBATCH --job-name=predict_genes
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

mkdir -p "$RESULTS/proteins"
use_env funannotate_env

for sample in $(tail -n +2 "$SAMPLES" | cut -f1); do
  out="$RESULTS/$sample"
  if [[ -s "$RESULTS/proteins/$sample.faa" ]]; then continue; fi   # finished in an earlier run
  mkdir -p "$out/funannotate_input"
  funannotate clean -i "$out/spades/scaffolds.fasta" -o "$out/funannotate_input/clean.fa"
  funannotate sort  -i "$out/funannotate_input/clean.fa" -o "$out/funannotate_input/sorted.fa"
  funannotate mask  -i "$out/funannotate_input/sorted.fa" -o "$out/funannotate_input/masked.fa" --cpus "$THREADS"

  funannotate predict -i "$out/funannotate_input/masked.fa" -o "$out/funannotate" \
    -s "$SPECIES" --isolate "$sample" -d "$FUNANNOTATE_DB" \
    --busco_db ascomycota --organism fungus --cpus "$THREADS"

  cp "$out"/funannotate/predict_results/*.proteins.fa "$RESULTS/proteins/$sample.faa"
done
