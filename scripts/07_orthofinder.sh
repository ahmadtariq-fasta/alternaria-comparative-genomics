#!/bin/bash
# scripts/07_orthofinder.sh
#
# Infer orthogroups across all isolates with OrthoFinder (DIAMOND + MCL).
#
# Input:  $RESULTS/proteins/*.faa
# Output: $RESULTS/proteins/OrthoFinder/Results_alternaria/Orthogroups/Orthogroups.GeneCount.tsv
#
# Run from the repo root: sbatch scripts/07_orthofinder.sh
#SBATCH --job-name=orthofinder
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
source config/config.sh

use_env orthofinder_env
orthofinder -f "$RESULTS/proteins" -S diamond -n alternaria -t "$THREADS" -a "$THREADS"
