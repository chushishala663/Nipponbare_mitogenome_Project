#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -N hifisr-variants
#PBS -l nodes=1:ppn=32
#PBS -l walltime=168:00:00

# Estimation of mitochondrial variant frequencies from HiFi reads.
# This script uses get_variants_in_reads.py from HiFi-SR:
# https://github.com/zouyinstein/hifisr
#
# Required environment variables:
#   HIFISR_DIR     HiFi-SR repository directory
#   SOFT_PATHS     HiFi-SR tab-delimited software-path file
#   SAMPLE         Sample identifier
#   ASSEMBLY_FASTA Mitochondrial assembly FASTA
#   FILTERED_READS Filtered mitochondrial HiFi reads in FASTQ format
#
# Optional environment variables:
#   WORK_DIR       Output directory (default: current directory)
#   RUN_LABEL      HiFi-SR run label (default: run_4)
#   THREADS        Number of threads (default: PBS_NP or 32)

set -euo pipefail

: "${HIFISR_DIR:?Set HIFISR_DIR to the HiFi-SR repository directory}"
: "${SOFT_PATHS:?Set SOFT_PATHS to the HiFi-SR software-path file}"
: "${SAMPLE:?Set SAMPLE to the sample identifier}"
: "${ASSEMBLY_FASTA:?Set ASSEMBLY_FASTA to the mitochondrial assembly FASTA}"
: "${FILTERED_READS:?Set FILTERED_READS to the filtered mitochondrial HiFi reads}"

WORK_DIR="${WORK_DIR:-$PWD}"
RUN_LABEL="${RUN_LABEL:-run_4}"
THREADS="${THREADS:-${PBS_NP:-32}}"
GET_VARIANTS="${HIFISR_DIR}/analysis_scripts/get_variants_in_reads.py"

for file in \
    "$SOFT_PATHS" \
    "$ASSEMBLY_FASTA" \
    "$FILTERED_READS" \
    "$GET_VARIANTS"; do
    [[ -f "$file" ]] || { echo "Required file not found: $file" >&2; exit 1; }
done

command -v python >/dev/null 2>&1 || {
    echo "Required command not found in PATH: python" >&2
    exit 1
}

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

python "$GET_VARIANTS" \
    "$SOFT_PATHS" \
    "$SAMPLE" \
    mito \
    "$RUN_LABEL" \
    "$ASSEMBLY_FASTA" \
    "$FILTERED_READS" \
    "$THREADS"
