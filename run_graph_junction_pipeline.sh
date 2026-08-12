#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -N graph-junctions
#PBS -l nodes=1:ppn=32
#PBS -l walltime=168:00:00

# End-to-end HiFi/ONT assembly-graph junction analysis with NUMT-aware filtering.
# Long reads are recruited with minimap2 and aligned to the complete GFA with
# GraphAligner. The resulting traversal events are summarized by distinct read
# counts and traversal-event counts.
#
# Required environment variables:
#   GFA                 Complete mitochondrial assembly graph in GFA format
#   MITO_REFERENCE      Corrected mitochondrial sequence in FASTA format
#   NUCLEAR_REFERENCE   Nuclear reference genome in FASTA format
#   NUMT_BED            Zero-based BED file of structurally defined NUMT loci
#   HIFI_READS          PacBio HiFi reads in FASTQ/FASTQ.GZ format
#   ONT_READS           ONT reads in FASTQ/FASTQ.GZ format
#
# Optional environment variables:
#   WORK_DIR            Output directory (default: graph_junction_analysis)
#   THREADS             Number of threads (default: PBS_NP or 32)
#   JUNCTIONS_TSV       Curated junction definitions; generated from GFA if unset
#   BRANCHES_TSV        Curated two-way branch definitions; generated if unset
#   JUNCTION_ID_MAP     Optional original-to-publication junction ID crosswalk

set -euo pipefail

: "${GFA:?Set GFA to the complete mitochondrial assembly graph}"
: "${MITO_REFERENCE:?Set MITO_REFERENCE to the corrected mitochondrial FASTA}"
: "${NUCLEAR_REFERENCE:?Set NUCLEAR_REFERENCE to the nuclear reference FASTA}"
: "${NUMT_BED:?Set NUMT_BED to the structurally defined NUMT BED file}"
: "${HIFI_READS:?Set HIFI_READS to the PacBio HiFi FASTQ/FASTQ.GZ file}"
: "${ONT_READS:?Set ONT_READS to the ONT FASTQ/FASTQ.GZ file}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR="${WORK_DIR:-graph_junction_analysis}"
THREADS="${THREADS:-${PBS_NP:-32}}"

for command_name in minimap2 GraphAligner seqkit python3; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "Required command not found in PATH: $command_name" >&2
        exit 1
    }
done
for file in "$GFA" "$MITO_REFERENCE" "$NUCLEAR_REFERENCE" "$NUMT_BED" "$HIFI_READS" "$ONT_READS"; do
    [[ -f "$file" ]] || { echo "Required file not found: $file" >&2; exit 1; }
done

mkdir -p "$WORK_DIR/metadata" "$WORK_DIR/HIFI" "$WORK_DIR/ONT" "$WORK_DIR/summary" "$WORK_DIR/logs"

generated_junctions="$WORK_DIR/metadata/junctions.tsv"
generated_branches="$WORK_DIR/metadata/branches.tsv"
metadata_arguments=(
    --gfa "$GFA"
    --junctions-output "$generated_junctions"
    --branches-output "$generated_branches"
)
if [[ -n "${JUNCTION_ID_MAP:-}" ]]; then
    metadata_arguments+=(--id-map "$JUNCTION_ID_MAP")
fi
python3 "$SCRIPT_DIR/prepare_graph_junctions.py" "${metadata_arguments[@]}"
JUNCTIONS_TSV="${JUNCTIONS_TSV:-$generated_junctions}"
BRANCHES_TSV="${BRANCHES_TSV:-$generated_branches}"

run_sample() {
    local sample=$1
    local reads=$2
    local preset=$3
    local recruitment_identity=$4
    local graph_identity=$5
    local precise_clipping=$6
    local outdir="$WORK_DIR/$sample"

    minimap2 -x "$preset" -t "$THREADS" -c --secondary=yes -N 100 \
        "$MITO_REFERENCE" "$reads" \
        > "$outdir/reads_to_mitogenome.paf" \
        2> "$WORK_DIR/logs/${sample}.mitogenome.minimap2.log"

    python3 "$SCRIPT_DIR/recruit_mitochondrial_reads.py" \
        --paf "$outdir/reads_to_mitogenome.paf" \
        --read-ids "$outdir/candidate_read_ids.txt" \
        --summary "$outdir/candidate_recruitment.tsv" \
        --minimum-block-length 300 \
        --minimum-identity "$recruitment_identity" \
        --minimum-union-query-bases 2000

    seqkit grep --threads "$THREADS" --pattern-file "$outdir/candidate_read_ids.txt" \
        --out-file "$outdir/candidate_reads.fastq.gz" "$reads"

    GraphAligner -g "$GFA" -f "$outdir/candidate_reads.fastq.gz" \
        -a "$outdir/candidate_reads_to_graph.gaf" -x vg -t "$THREADS" \
        --precise-clipping "$precise_clipping" \
        2> "$WORK_DIR/logs/${sample}.graphaligner.log"

    minimap2 -x "$preset" -t "$THREADS" -c --secondary=yes -N 100 \
        "$NUCLEAR_REFERENCE" "$outdir/candidate_reads.fastq.gz" \
        > "$outdir/candidate_reads_to_nuclear.paf" \
        2> "$WORK_DIR/logs/${sample}.nuclear.minimap2.log"

    local event_arguments=(
        --sample "$sample" \
        --gaf "$outdir/candidate_reads_to_graph.gaf" \
        --nuclear-paf "$outdir/candidate_reads_to_nuclear.paf" \
        --numt-bed "$NUMT_BED" \
        --junctions "$JUNCTIONS_TSV" \
        --minimum-graph-identity "$graph_identity" \
        --minimum-nuclear-identity "$recruitment_identity" \
        --events-output "$outdir/traversal_events.tsv" \
        --provenance-output "$outdir/read_provenance.tsv"
    )
    python3 "$SCRIPT_DIR/extract_graph_junction_events.py" "${event_arguments[@]}"
}

run_sample HIFI "$HIFI_READS" map-hifi 0.90 0.90 0.90
run_sample ONT "$ONT_READS" map-ont 0.70 0.70 0.70

python3 - "$WORK_DIR/HIFI/traversal_events.tsv" "$WORK_DIR/ONT/traversal_events.tsv" "$WORK_DIR/traversal_events.tsv" <<'PY'
import sys
from pathlib import Path

inputs = [Path(sys.argv[1]), Path(sys.argv[2])]
output = Path(sys.argv[3])
with output.open("w", encoding="utf-8") as destination:
    for index, path in enumerate(inputs):
        with path.open(encoding="utf-8") as source:
            header = source.readline()
            if index == 0:
                destination.write(header)
            destination.writelines(source)
PY

python3 "$SCRIPT_DIR/summarize_graph_junctions.py" \
    --events "$WORK_DIR/traversal_events.tsv" \
    --branches "$BRANCHES_TSV" \
    --output-dir "$WORK_DIR/summary" \
    --minimum-graph-query-coverage 0.95 \
    --exclude-provenance nuclear_boundary_NUMT
