# Nipponbare_mitogenome_Project

This GitHub repository contains scripts used for analyses presented in the following paper:

**Reconstruction and Correction of the Rice Mitochondrial Reference Genome**

## File Overview

- `call_mitogenome_variants.py` - Identifies single-nucleotide variants and insertions/deletions between mitochondrial genome assemblies. BLAST-defined collinear blocks are refined with MAFFT, and variants are reported with coordinates in both assemblies.

- `plot_mitogenome_circos.R` - Generates a circular mitochondrial genome plot showing genome coordinates, GC content, annotated genes, repeat regions, and internal repeat-associated links.

- `summarize_graph_junctions.py` - Summarizes long-read support for assembly-graph junctions. It reports distinct supporting reads and total traversal events, calculates the relative support for competing junctions at two-way branches, and estimates read-cluster bootstrap confidence intervals.

- `prepare_graph_junctions.py` - Extracts oriented graph junctions and two-way branch definitions from a mitochondrial assembly graph in GFA format.

- `recruit_mitochondrial_reads.py` - Selects candidate mitochondrial reads from read-to-mitogenome PAF alignments using unioned query coverage.

- `extract_graph_junction_events.py` - Converts complete-graph GAF paths into oriented junction traversal events and classifies reads using nuclear-boundary and internal-NUMT evidence.

- `run_graph_junction_pipeline.sh` - Runs the NUMT-aware graph-junction analysis from HiFi and ONT FASTQ files through `traversal_events.tsv` and the final junction-support summaries.

- `assemble_mitogenome_with_hifisr.sh` - Extracts mitochondrial HiFi reads, filters reads by length and quality, performs reproducible random subsampling, and generates a draft mitochondrial assembly with [HiFi-SR](https://github.com/zouyinstein/hifisr).

- `estimate_hifisr_variant_frequencies.sh` - Estimates mitochondrial variant frequencies from filtered HiFi reads using [HiFi-SR](https://github.com/zouyinstein/hifisr).

- `README.md` - Project overview and descriptions of the scripts.

- `LICENSE` - MIT license for the scripts in this repository.
