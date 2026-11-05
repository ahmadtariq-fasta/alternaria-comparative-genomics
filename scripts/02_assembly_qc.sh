#!/bin/bash
# scripts/02_assembly_qc.sh
#
# Assembly statistics (QUAST) and gene-space completeness (BUSCO).
#
# Input:  $RESULTS/<sample>/spades/scaffolds.fasta
# Output: $RESULTS/<sample>/quast/, $RESULTS/<sample>/busco/
#
# Run from the repo root: sbatch scripts/02_assembly_qc.sh
#SBATCH --job-name=assembly_qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

for sample in $(tail -n +2 "$SAMPLES" | cut -f1); do
  out="$RESULTS/$sample"
  use_env quast
  quast.py -t "$THREADS" -o "$out/quast" "$out/spades/scaffolds.fasta"
  use_env busco
  busco -i "$out/spades/scaffolds.fasta" -l "$BUSCO_LINEAGE" -m genome -c "$THREADS" \
    --offline --download_path "$BUSCO_DB" --out_path "$out" -o busco -f
done
