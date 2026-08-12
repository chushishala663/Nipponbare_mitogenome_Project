#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -N hifisr-mt
#PBS -l nodes=1:ppn=32
#PBS -l walltime=168:00:00

# Mitochondrial HiFi-read extraction, filtering, subsampling, and draft assembly.
# This script uses analysis scripts from HiFi-SR:
# https://github.com/zouyinstein/hifisr
#
# Required environment variables:
#   HIFISR_DIR        HiFi-SR repository directory
#   SOFT_PATHS        HiFi-SR tab-delimited software-path file
#   SAMPLE            Sample identifier
#   MITO_REFERENCE    Mitochondrial reference FASTA (may contain multiple records)
#   PLASTID_REFERENCE Plastid reference FASTA
#   HIFI_READS        HiFi reads in FASTQ or FASTQ.GZ format
#
# Optional environment variables:
#   WORK_DIR          Output directory (default: current directory)
#   THREADS           Number of threads (default: PBS_NP or 32)
#   MIN_READ_LENGTH   Minimum read length in bp (default: 10000)
#   MIN_READ_QUALITY  Minimum read-quality cutoff (default: 0)
#   SUBSAMPLE_SIZE    Number of mitochondrial reads sampled (default: 8000)
#   RANDOM_SEED       Seed used by seqtk sample (default: 42)
#   FILTERED_IDS      Read-ID file produced by filt_read_ids.py
#   FILTERED_READS    FASTQ file created from FILTERED_IDS

set -euo pipefail

: "${HIFISR_DIR:?Set HIFISR_DIR to the HiFi-SR repository directory}"
: "${SOFT_PATHS:?Set SOFT_PATHS to the HiFi-SR software-path file}"
: "${SAMPLE:?Set SAMPLE to the sample identifier}"
: "${MITO_REFERENCE:?Set MITO_REFERENCE to the mitochondrial reference FASTA}"
: "${PLASTID_REFERENCE:?Set PLASTID_REFERENCE to the plastid reference FASTA}"
: "${HIFI_READS:?Set HIFI_READS to the HiFi FASTQ/FASTQ.GZ file}"

WORK_DIR="${WORK_DIR:-$PWD}"
THREADS="${THREADS:-${PBS_NP:-32}}"
MIN_READ_LENGTH="${MIN_READ_LENGTH:-10000}"
MIN_READ_QUALITY="${MIN_READ_QUALITY:-0}"
SUBSAMPLE_SIZE="${SUBSAMPLE_SIZE:-8000}"
RANDOM_SEED="${RANDOM_SEED:-42}"

GET_READS="${HIFISR_DIR}/analysis_scripts/get_mtpt_reads.py"
FILTER_IDS="${HIFISR_DIR}/analysis_scripts/filt_read_ids.py"
GET_ASSEMBLY="${HIFISR_DIR}/analysis_scripts/get_draft_assembly.py"

for file in \
    "$SOFT_PATHS" \
    "$MITO_REFERENCE" \
    "$PLASTID_REFERENCE" \
    "$HIFI_READS" \
    "$GET_READS" \
    "$FILTER_IDS" \
    "$GET_ASSEMBLY"; do
    [[ -f "$file" ]] || { echo "Required file not found: $file" >&2; exit 1; }
done

for command_name in python seqkit seqtk; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "Required command not found in PATH: $command_name" >&2
        exit 1
    }
done

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

python "$GET_READS" \
    "$SOFT_PATHS" \
    "$SAMPLE" \
    "$MITO_REFERENCE" \
    "$PLASTID_REFERENCE" \
    "$HIFI_READS" \
    "$THREADS"

python "$FILTER_IDS" \
    "$SOFT_PATHS" \
    "$SAMPLE" \
    mito_id_length_qual.txt \
    plastid_id_length_qual.txt \
    "$MIN_READ_LENGTH" \
    "$MIN_READ_QUALITY"

MITO_READS="${SAMPLE}/reads/${SAMPLE}_mito.fastq"
FILTERED_IDS="${FILTERED_IDS:-${SAMPLE}/reads/filt_reads/filt_L10K_mito_ids.txt}"
FILTERED_READS="${FILTERED_READS:-${SAMPLE}/reads/filt_reads/mito_filt_L10K_mito.fastq}"
SAMPLED_READS="${SAMPLE}/reads/sample_reads/sample_${SUBSAMPLE_SIZE}_mito.fastq"
ASSEMBLY_SAMPLE="${SAMPLE}_${SUBSAMPLE_SIZE}"

mkdir -p "$(dirname "$SAMPLED_READS")"

seqkit grep \
    --threads "$THREADS" \
    --pattern-file "$FILTERED_IDS" \
    --out-file "$FILTERED_READS" \
    "$MITO_READS"

seqtk sample -s"$RANDOM_SEED" "$FILTERED_READS" "$SUBSAMPLE_SIZE" \
    > "$SAMPLED_READS"

python "$GET_ASSEMBLY" \
    "$SOFT_PATHS" \
    "$ASSEMBLY_SAMPLE" \
    mito \
    "$MITO_REFERENCE" \
    "$SAMPLED_READS" \
    "$THREADS"
