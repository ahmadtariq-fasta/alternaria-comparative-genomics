#!/bin/bash
# scripts/05_secretome_effectors.sh
#
# Find secreted proteins with SignalP, then predict effectors among them with EffectorP.
#
# Input:  $RESULTS/proteins/<sample>.faa
# Output: $RESULTS/<sample>/signalp/, $RESULTS/<sample>/effectorp/effectorp.txt
#
# Run from the repo root: sbatch scripts/05_secretome_effectors.sh
#SBATCH --job-name=secretome_effectors
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

for sample in $(tail -n +2 "$SAMPLES" | cut -f1); do
  out="$RESULTS/$sample"
  mkdir -p "$out/signalp" "$out/effectorp"

  # SignalP needs its own licence; without it the effector counts stay empty
  if ! command -v signalp6 > /dev/null; then
    echo "⚠️  signalp6 not found: no effector predictions for $sample" >&2
    : > "$out/effectorp/effectorp.txt"
    continue
  fi

  signalp6 --fastafile "$RESULTS/proteins/$sample.faa" --organism eukarya \
    --output_dir "$out/signalp" --format txt --mode slow

  # keep only the proteins SignalP calls secreted (Prediction = SP)
  awk -F'\t' '!/^#/ && $2=="SP"{split($1,a," "); print a[1]}' \
    "$out/signalp/prediction_results.txt" > "$out/signalp/secreted_ids.txt"
  awk 'NR==FNR{keep[$1]=1; next} /^>/{id=substr($1,2); show=(id in keep)} show' \
    "$out/signalp/secreted_ids.txt" "$RESULTS/proteins/$sample.faa" > "$out/signalp/secreted.faa"

  python "$EFFECTORP" -f -i "$out/signalp/secreted.faa" -o "$out/effectorp/effectorp.txt"
done
