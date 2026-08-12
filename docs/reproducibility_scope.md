# Reproducibility scope

This code package separates three analysis levels.

## Included

1. **Assembly-to-assembly sequence differences**
   - Input: two FASTA files and BLAST outfmt 6 collinear blocks.
   - Output: coordinate-aware SNV and indel table after MAFFT refinement.

2. **Circular coordinate visualization**
   - Input: prepared gene, repeat, GC, and internal-link tables.
   - Output: PDF, SVG, or PNG circular plot.

3. **Graph-junction summary statistics**
   - Input: one event-level row for every read traversal of a graph junction,
     including graph query coverage and NUMT provenance class.
   - Output: distinct-read counts, traversal-event counts, two-way branch
     composition, and whole-read cluster-bootstrap intervals.

## Not yet included

The current public-ready package does not generate the event-level traversal
table directly from raw HiFi/ONT reads. That upstream step requires the exact
complete-GFA alignment, split-GAF joining, and NUMT provenance-classification
workflow used for the final analysis. The summary script must not be described
as a raw-read-to-result pipeline until those upstream scripts and their required
reference inputs are added and validated.

## Interpretation

Graph-junction fractions are local read-level compositions within a sequencing
library. They are not biological-replicate estimates and do not establish the
stoichiometry of complete chromosome-sized mitochondrial molecules.
