# Input formats

All coordinate tables are tab-separated, have no header, and use one-based
inclusive coordinates unless the source data are explicitly converted before
use.

## `plot_mitogenome_circos.R`

### Repeats (`--repeats`)

```text
sequence_id  start  end  value  color
```

`color` may be a hexadecimal color, an R color name, or `R,G,B` integers.

### Genes (`--genes`)

```text
sequence_id  start  end  color
```

### GC track (`--gc`)

```text
sequence_id  start  end  gc_fraction
```

### Internal links (`--links`)

```text
query_id  query_start  query_end  subject_id  subject_start  subject_end
```

Exact full-length self matches are removed automatically.

## `call_mitogenome_variants.py`

The BLAST input must contain the standard 12 columns produced by:

```bash
blastn -query query.fasta -subject reference.fasta \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
```

The output is a headered TSV. Insertions are represented between
`reference_start` and `reference_end` with an empty `ref` allele. Deletions use
the deleted reference interval and an empty `alt` allele. These are explicit
event coordinates rather than normalized VCF alleles.

## `summarize_graph_junctions.py`

The event table is headered and contains one row per observed traversal:

```text
sample  read_id  junction_id  graph_query_coverage  provenance_class
```

Repeated rows for the same `read_id` and `junction_id` represent multiple
traversal events carried by one read. `distinct_reads` counts the read once,
whereas `traversal_events` counts every occurrence. By default, rows classified
as `nuclear_boundary_NUMT` are excluded and internally NUMT-ambiguous rows are
retained.

The branch-definition table is headered:

```text
branch_id  oriented_endpoint  junction_1  junction_2
```
