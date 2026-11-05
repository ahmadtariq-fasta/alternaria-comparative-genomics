# config/config.sh: settings shared by all pipeline steps. Edit the paths before running.
# Values already set in the environment (for example by cluster/run_all.sh) take priority.
SAMPLES="${SAMPLES:-config/samples.tsv}"
RESULTS="${RESULTS:-results}"
THREADS="${THREADS:-16}"
ADAPTERS="${ADAPTERS:-/path/to/TruSeq3-PE-2.fa}"       # Trimmomatic adapter file
BUSCO_LINEAGE=ascomycota_odb10
BUSCO_DB="${BUSCO_DB:-busco_downloads}"
SPECIES="Alternaria alternata"
DBCAN_DB="${DBCAN_DB:-/path/to/dbcan_db}"
EMAPPER_DB="${EMAPPER_DB:-/path/to/eggnog_db}"
EFFECTORP="${EFFECTORP:-/path/to/EffectorP.py}"
FUNANNOTATE_DB="${FUNANNOTATE_DB:-$HOME/.funannotate_db}"
REFERENCE_HOST="${REFERENCE_HOST:-wild_tomato}"

# Conda: each group of tools can live in its own environment. Set ENV_ROOT to the folder
# that holds them, or leave it empty if everything is already on PATH.
use_env() {
  if [[ -n "${ENV_ROOT:-}" ]]; then conda activate "$ENV_ROOT/$1"; fi
}
if [[ -n "${ENV_ROOT:-}" ]]; then eval "$(conda shell.bash hook)"; fi
