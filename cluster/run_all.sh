#!/bin/bash
# cluster/run_all.sh
#
# Runs the whole pipeline on public Alternaria genomes on the CAU cluster, separately from
# any other project. On a login node it installs the conda environments, downloads the reads
# and databases, then submits itself; compute nodes have no internet. Finished steps are
# skipped, so a job that hits its time limit can be resubmitted the same way.
#
# Input:  cluster/public_samples.tsv (run accessions from ENA)
# Output: $PROJECT_ROOT/alternaria_run/, results/public_data/
#
# Run from the repo root on a login node: bash cluster/run_all.sh
#SBATCH --job-name=alternaria_cg
#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=%x_%j.out

set -eo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/work_beegfs/$USER/portfolio}"
export ENV_ROOT="$PROJECT_ROOT/conda"
RUN="$PROJECT_ROOT/alternaria_run"
export RESULTS="$RUN/results"
export SAMPLES="$RUN/samples.tsv"
export FUNANNOTATE_DB="$PROJECT_ROOT/funannotate_db"
export BUSCO_DB="$PROJECT_ROOT/busco_downloads"
export DBCAN_DB="$PROJECT_ROOT/dbcan_db"
export REFERENCE_HOST=malus
export TARGET_BASES=1400000000     # about 40x of a 35 Mb genome

module load gcc12-env/12.1.0 miniconda3
eval "$(conda shell.bash hook)"

# make_env <name> <tool that must exist afterwards> <packages...>
make_env() {
  local name=$1 tool=$2; shift 2
  [[ -x "$ENV_ROOT/$name/bin/$tool" ]] && return 0
  rm -rf "$ENV_ROOT/$name"     # left over from an install that did not finish
  echo "Creating conda environment: $name"
  conda create -y -p "$ENV_ROOT/$name" --override-channels -c conda-forge -c bioconda "$@"
}

check_tool() {
  local env=$1; shift
  conda activate "$ENV_ROOT/$env"
  for tool in "$@"; do
    command -v "$tool" > /dev/null || { echo "❌ $tool is missing from $env" >&2; exit 1; }
  done
  conda deactivate
}

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  # ---------- login node: install, download, submit ----------
  mkdir -p "$ENV_ROOT" "$RUN/reads" "$FUNANNOTATE_DB" "$BUSCO_DB" "$DBCAN_DB"
  make_env download seqtk seqtk
  make_env trim_quality trimmomatic trimmomatic=0.39 fastqc
  make_env assembly spades.py spades
  make_env quast quast.py quast "python=3.11"   # QUAST imports distutils, removed in Python 3.12
  make_env busco busco "busco=5"
  make_env funannotate_env funannotate funannotate
  make_env dbcan run_dbcan "dbcan=4.1.4"
  make_env antismash antismash antismash
  make_env orthofinder_env orthofinder orthofinder
  make_env r_stats Rscript r-base

  check_tool download seqtk
  check_tool trim_quality trimmomatic
  check_tool assembly spades.py
  check_tool quast quast.py
  "$ENV_ROOT/quast/bin/quast.py" --version
  check_tool busco busco
  check_tool funannotate_env funannotate
  check_tool dbcan run_dbcan dbcan_build
  check_tool antismash antismash download-antismash-databases
  check_tool orthofinder_env orthofinder
  check_tool r_stats Rscript

  # sample sheet: sample, host, reads
  { echo -e "sample\thost\tr1\tr2"
    tail -n +2 cluster/public_samples.tsv | while IFS=$'\t' read -r id host; do
      echo -e "$id\t$host\t$RUN/reads/${id}_1.fq.gz\t$RUN/reads/${id}_2.fq.gz"
    done; } > "$SAMPLES"

  # reads from ENA, reduced to about TARGET_BASES per isolate
  conda activate "$ENV_ROOT/download"
  for ID in $(tail -n +2 cluster/public_samples.tsv | cut -f1); do
    [[ -s "$RUN/reads/${ID}_2.fq.gz" ]] && continue
    echo "Downloading $ID"
    info=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$ID&result=read_run&fields=base_count,fastq_ftp&format=tsv" | tail -1)
    bases=$(cut -f2 <<< "$info")
    urls=$(cut -f3 <<< "$info" | tr ';' '\n' | grep -E '_[12]\.fastq\.gz$')
    [[ $(wc -l <<< "$urls") -eq 2 ]] || { echo "❌ $ID does not have two FASTQ files" >&2; exit 1; }
    fraction=$(awk -v t="$TARGET_BASES" -v b="$bases" 'BEGIN { print (b > t) ? t / b : 1 }')
    n=1
    for url in $urls; do
      wget -q -c -O "$RUN/reads/${ID}_raw_$n.fastq.gz" "https://$url"
      seqtk sample -s 11 "$RUN/reads/${ID}_raw_$n.fastq.gz" "$fraction" | gzip > "$RUN/reads/${ID}_$n.part.gz"
      mv "$RUN/reads/${ID}_$n.part.gz" "$RUN/reads/${ID}_$n.fq.gz"
      rm "$RUN/reads/${ID}_raw_$n.fastq.gz"
      n=$((n + 1))
    done
  done
  conda deactivate

  # databases (each downloaded once)
  if [[ ! -f "$FUNANNOTATE_DB/.setup_done" ]]; then
    conda activate "$ENV_ROOT/funannotate_env"
    funannotate setup -i all -b ascomycota -d "$FUNANNOTATE_DB"
    conda deactivate
    touch "$FUNANNOTATE_DB/.setup_done"
  fi
  if [[ ! -f "$BUSCO_DB/.setup_done" ]]; then
    conda activate "$ENV_ROOT/busco"
    busco --download ascomycota_odb10 --download_path "$BUSCO_DB"
    conda deactivate
    touch "$BUSCO_DB/.setup_done"
  fi
  if [[ ! -f "$DBCAN_DB/.setup_done" ]]; then
    conda activate "$ENV_ROOT/dbcan"
    dbcan_build --cpus 8 --db-dir "$DBCAN_DB" --clean
    conda deactivate
    touch "$DBCAN_DB/.setup_done"
  fi
  if [[ ! -f "$PROJECT_ROOT/.antismash_db_done" ]]; then
    conda activate "$ENV_ROOT/antismash"
    download-antismash-databases
    conda deactivate
    touch "$PROJECT_ROOT/.antismash_db_done"
  fi

  sbatch cluster/run_all.sh
  echo "✅ Submitted. Follow it with: squeue -u $USER"
  exit 0
fi

# ---------- compute node: run every step in order ----------
cd "$SLURM_SUBMIT_DIR"
export THREADS="$SLURM_CPUS_PER_TASK"
export ADAPTERS
ADAPTERS=$(find "$ENV_ROOT/trim_quality/share" -name 'TruSeq3-PE-2.fa' | head -1)

# Step 4 (eggNOG-mapper, InterProScan) is left out: its databases need well over 100 GB
# and steps 8-10 do not use its output.
mkdir -p "$RUN/.done"
for step in 01_trim_assemble 02_assembly_qc 03_predict_genes 05_secretome_effectors 06_cazymes_bgcs 07_orthofinder; do
  [[ -f "$RUN/.done/$step" ]] && { echo "Skipping $step (already finished)"; continue; }
  echo "=== $step: $(date)"
  bash "scripts/$step.sh"
  touch "$RUN/.done/$step"
done

conda activate "$ENV_ROOT/r_stats"
Rscript R/08_build_feature_tables.R "$SAMPLES" "$RESULTS" \
  "$RESULTS/proteins/OrthoFinder/Results_alternaria/Orthogroups/Orthogroups.GeneCount.tsv" "$RESULTS/tables"
Rscript R/09_host_enrichment.R "$SAMPLES" "$RESULTS/tables" "$RESULTS/enrichment" "$REFERENCE_HOST"
Rscript R/10_plot_counts.R "$RESULTS/tables/counts_per_isolate.tsv" "$RESULTS/enrichment" \
  "Public Alternaria genomes: pipeline test, hosts confounded with study and species"

mkdir -p results/public_data
cp "$RESULTS"/enrichment/* "$RESULTS"/tables/counts_per_isolate.tsv results/public_data/
echo "✅ Finished: $(date). Results in results/public_data/"
