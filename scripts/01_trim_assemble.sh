#!/bin/bash
# scripts/01_trim_assemble.sh
#
# Quality-trim reads with Trimmomatic, then assemble each isolate with SPAdes.
#
# Input:  paired reads listed in config/samples.tsv
# Output: $RESULTS/<sample>/trimmed/, $RESULTS/<sample>/spades/scaffolds.fasta
#
# Run from the repo root: sbatch scripts/01_trim_assemble.sh
#SBATCH --job-name=trim_assemble
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

for sample in $(tail -n +2 "$SAMPLES" | cut -f1); do
  r1=$(awk -F'\t' -v s="$sample" '$1==s{print $3}' "$SAMPLES")
  r2=$(awk -F'\t' -v s="$sample" '$1==s{print $4}' "$SAMPLES")
  out="$RESULTS/$sample"
  mkdir -p "$out/trimmed"
  if [[ -s "$out/spades/scaffolds.fasta" ]]; then continue; fi   # finished in an earlier run

  use_env trim_quality
  trimmomatic PE -phred33 -threads "$THREADS" "$r1" "$r2" \
    "$out/trimmed/R1.fastq.gz" "$out/trimmed/R1_unpaired.fastq.gz" \
    "$out/trimmed/R2.fastq.gz" "$out/trimmed/R2_unpaired.fastq.gz" \
    "ILLUMINACLIP:$ADAPTERS:2:30:10" SLIDINGWINDOW:4:20 MINLEN:50

  use_env assembly
  spades.py --careful -t "$THREADS" \
    -1 "$out/trimmed/R1.fastq.gz" -2 "$out/trimmed/R2.fastq.gz" \
    -o "$out/spades"
done
