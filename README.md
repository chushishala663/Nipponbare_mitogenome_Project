# Nipponbare mitochondrial genome analysis scripts

This repository contains reusable scripts associated with the study
**“Reconstruction and Correction of the Rice Mitochondrial Reference Genome.”**

Release status: pre-publication code package (`v0.1.0`).

Repository: <https://github.com/chushishala663/Nipponbare_mitogenome_Project>

The scripts are parameterized and contain no machine-specific paths. Large raw
sequencing data are not stored in the repository and can be obtained from the
Sequence Read Archive:

- PacBio HiFi: `SRR25241090`
- Oxford Nanopore ultra-long reads: `SRR25241091`
- Illumina paired-end reads: `SRR25241092`

## Contents

- `scripts/call_mitogenome_variants.py`: identifies SNVs and indels within
  BLAST-defined collinear blocks after MAFFT refinement.
- `scripts/plot_mitogenome_circos.R`: produces a circular mitochondrial genome
  plot containing coordinates, GC content, genes, repeats, and internal links.
- `scripts/summarize_graph_junctions.py`: summarizes event-level graph
  traversals as distinct-read support, traversal-event support, two-way branch
  composition, and whole-read cluster-bootstrap confidence intervals.
- `docs/input_formats.md`: exact input schemas.
- `docs/reproducibility_scope.md`: analyses covered by this release and the
  explicit boundary of the graph-junction workflow.
- `environment.yml`: Conda environment for both scripts.

## Installation

```bash
conda env create -f environment.yml
conda activate nipponbare-mitogenome
```

## Variant identification

First generate standard BLAST outfmt 6 blocks, then run:

```bash
python scripts/call_mitogenome_variants.py \
  --query query.fasta \
  --reference reference.fasta \
  --blast query_vs_reference.blast.tsv \
  --min-block-length 100 \
  --output results/query_vs_reference.variants.tsv
```

The query and reference FASTA identifiers must match the first two columns of
the BLAST file. For multi-record FASTA files, use `--query-id` and
`--reference-id`.

## Circular genome plot

```bash
Rscript scripts/plot_mitogenome_circos.R \
  --repeats repeats.tsv \
  --genes genes.tsv \
  --gc gc.tsv \
  --links internal_links.tsv \
  --genome-length 376041 \
  --output results/mitogenome_circos.pdf
```

PDF and SVG retain vector elements; PNG is exported at 300 dpi. See
`docs/input_formats.md` for table layouts.

## Interpretation and limitations

The variant script calls differences only within the supplied BLAST blocks. It
does not infer rearrangements from unmatched sequence and does not replace a
whole-genome structural-variant caller. Overlapping BLAST blocks can generate
duplicate calls; identical calls are collapsed by default.

The circular plot is a coordinate-based visualization and should not be
interpreted as evidence that the physical plant mitochondrial genome exists as
a single circular molecule.

## Graph-junction support

```bash
python scripts/summarize_graph_junctions.py \
  --events junction_events.tsv \
  --branches branches.tsv \
  --minimum-graph-query-coverage 0.95 \
  --exclude-provenance nuclear_boundary_NUMT \
  --bootstrap-replicates 2000 \
  --output-dir results/junctions
```

The input is an event-level export from a complete-GFA read-alignment workflow.
A read that traverses one junction more than once contributes one distinct read
but multiple traversal events. Bootstrap resampling uses whole read identifiers
as clusters so events from the same read remain correlated. The resulting
fractions describe local junction composition in each sequencing library, not
the abundance of complete chromosome-sized mitochondrial molecules.

## Citation

Please cite the associated article when using these scripts. The final article
citation and repository DOI will be added after publication.
