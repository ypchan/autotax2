# AutoTax2

AutoTax2 is a GTDB-anchored, incremental SSU taxonomy construction framework.
It uses a user-built `gtdb_ssu.arb` reference with SINA for orientation and
anchor inference, and vsearch for efficient local hierarchical clustering of
novel lineages below the trusted anchor rank.

## Key ideas

- No SILVA dependency is required.
- GTDB-derived SSU reference is treated as a stable backbone.
- SINA aligns query sequences to `gtdb_ssu.arb` and provides orientation,
  alignment identity, alignment quality, and LCA taxonomy fields in FASTA headers.
- vsearch clusters only the novel sequence pools below the trusted GTDB anchor.
- Placeholder IDs are global within each rank, not within each source.
- Source prefixes indicate the first dataset that created the placeholder.
- Cluster membership is retained for every rank to support source-specific and
  shared-lineage analysis.

## Installation

```bash
pip install -e .
```

External tools required in `$PATH`:

```bash
sina --version
vsearch --version
```

Recommended versions:

- SINA 1.7.2
- vsearch 2.30.6_linux_x86_64

## Quick start

### 1. Check dependencies

```bash
autotax2 check --debug
```

### 2. Initialize a GTDB backbone database

```bash
autotax2 init \
  --ref-fa gtdb_ssu_tax.fa \
  --ref-tax gtdb_ssu.taxonomy.tsv \
  --ref-arb gtdb_ssu.arb \
  --db autotax2_db \
  --threads 64 \
  --debug
```

Required columns in `gtdb_ssu.taxonomy.tsv`:

```text
seq_id	domain	phylum	class	order	family	genus	species
```

Optional columns are preserved when present:

```text
evidence_level	source	length
```

### 3. Add a new dataset incrementally

```bash
autotax2 add \
  --db autotax2_db \
  --input MiDAS5.3.ssu.fa \
  --source midas \
  --prefix midas \
  --threads 64 \
  --mode incremental \
  --debug
```

### 4. Export annotation references

```bash
autotax2 export \
  --db autotax2_db \
  --format all \
  --out export/
```

Supported export formats:

- `sintax`: vsearch SINTAX FASTA
- `dada2`: DADA2-compatible FASTA
- `qiime2`: QIIME2 sequence FASTA + taxonomy TSV
- `taxonomy`: flat taxonomy TSV
- `all`: all of the above

## SINA header requirements

AutoTax2 expects SINA output FASTA headers to contain fields such as:

```text
[align_ident_slv=88.5] [align_quality_slv=90] [lca_tax_gtdb=d__Bacteria;p__...;] [turn=reversed and complemented]
```

Although the field names may include `slv`, AutoTax2 treats them as reference
metrics from your `gtdb_ssu.arb` run. It can also parse normalized fields:

```text
[align_ident_ref=88.5] [align_quality_ref=90] [lca_tax_ref=d__Bacteria;p__...;]
```

## Placeholder logic

Thresholds are percentage identities from SINA/vsearch:

| rank | cutoff |
|---|---:|
| species | 97.2 |
| genus | 90.1 |
| family | 80.1 |
| order | 72.9 |
| class | 72.2 |
| phylum | 69.6 |

If a new sequence anchors to a GTDB family at 85% identity, AutoTax2 inherits
family/order/class/phylum from GTDB and constructs novel genus/species clusters
below that family.

New placeholders are numbered globally per rank:

```text
s__midas_s000001
s__mfd_s000002
g__midas_g000001
g__hifimeta_g000002
```

## Development

```bash
pip install -e '.[dev]'
pytest
ruff check autotax2 tests
```
