# AutoTax2

AutoTax2 is a VSEARCH-based workflow for SILVA-backed rRNA sequence processing, taxonomy assignment, rank-wise clustering, reference extension, and source-overlap analysis.

Version: `0.1`

## Design goals

AutoTax2 is designed around a fixed SILVA backbone.

1. SILVA remains the backbone and is not re-clustered together with user sequences.
2. User datasets and extension references are first inserted into the SILVA framework.
3. Sequences that match SILVA above a rank-specific identity threshold inherit the corresponding SILVA-backed taxon assignment.
4. Sequences that do not match SILVA above a threshold are clustered within the extension dataset and reported as novel rank-like taxa.
5. VSEARCH performs the search and clustering work; Python handles file preparation, metadata parsing, workflow logic, and result summaries.
6. AutoTax2 does not download SILVA automatically. Users provide local SILVA FASTA and metadata files.

## Installation

AutoTax2 is a Python package. Install it into an existing Python environment.

```bash
python -m pip install -U pip
python -m pip install -e .
```

VSEARCH must be available on `PATH`, or you must provide its path with `--vsearch`.

Check the installation:

```bash
autotax2 --help
autotax2 check
vsearch --version
```

## Command overview

```bash
autotax2 --help
```

Commands:

```text
check              Check external dependencies and optional reference files.
prepare-silva      Prepare local SILVA FASTA and metadata files.
detect-intron      Detect intron-like insertions and write analysis FASTA files.
insert-backbone    Insert extension sequences into the SILVA backbone.
overlap-backbone   Summarize taxon overlap across backbone assignment files.
derep              Run VSEARCH dereplication.
cluster            Run VSEARCH clustering at one or more identity levels.
classify           Classify representative sequences using SINTAX plus SILVA/type-strain hits.
assign             Assign new sequences to old centroids or create new clusters.
provenance         Summarize source composition from UC clustering levels.
summarize          Summarize existing classify outputs.
run                Run dereplication, clustering and classification in one workflow.
```

Every command supports:

```bash
autotax2 <command> --help
autotax2 <command> --example
```

## Prepare SILVA

AutoTax2 expects local SILVA files, for example:

```text
SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
SILVA_138.2_SSURef.full_metadata.gz
```

Prepare a local reference directory:

```bash
autotax2 prepare-silva \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --silva-metadata SILVA_138.2_SSURef.full_metadata.gz \
  --out refdatabases \
  --make-udb \
  --threads auto
```

`prepare-silva` cleans the input FASTA before downstream processing. The default behavior is strict:

```text
U/u is converted to T
only A/C/G/T is allowed after conversion
records containing N, ambiguity codes, gaps, dots, or any other non-ACGT character are dropped
a summary report and a dropped-record report are generated
```

The output directory contains files such as:

```text
refdatabases/
  <prefix>.fasta
  <prefix>.udb
  <prefix>_sintax.fasta
  <prefix>_sintax.udb
  <prefix>_typestrains.fasta
  <prefix>_typestrains.udb
  <prefix>_raw_input.fasta
  <prefix>_fasta_cleaning_summary.tsv
  <prefix>_fasta_cleaning_dropped.tsv
  silva_taxonomy.tsv
  typestrains_accessionIDs.txt
  autotax2_ref_manifest.tsv
```

## Optional intron detection

Use this when long-read 16S sequences may contain intron-like insertions.

```bash
autotax2 detect-intron \
  --input hifimeta.original.fa \
  --db refdatabases/<prefix>.udb \
  --source-label hifimeta \
  --out hifimeta_intron \
  --search-id 0.70 \
  --rescue-id 0.987 \
  --min-intron-len 50 \
  --min-flank-len 150 \
  --threads auto
```

Important outputs:

```text
hifimeta_intron/analysis_sequences.fa
hifimeta_intron/sequence_version_map.tsv
```

The analysis FASTA is used for matching and clustering. The original FASTA can still be used for final centroid output.

## Insert sequences into the SILVA backbone

After intron detection:

```bash
autotax2 insert-backbone \
  --input hifimeta_intron/analysis_sequences.fa \
  --original-fasta hifimeta.original.fa \
  --version-map hifimeta_intron/sequence_version_map.tsv \
  --source-label hifimeta \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --rank-thresholds default \
  --db-format auto \
  --out hifimeta_inserted \
  --threads auto
```

Without intron detection:

```bash
autotax2 insert-backbone \
  --input hifimeta.original.fa \
  --original-fasta hifimeta.original.fa \
  --source-label hifimeta \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --rank-thresholds default \
  --db-format auto \
  --out hifimeta_inserted \
  --threads auto
```

Common outputs include:

```text
hifimeta_inserted/
  sequence_rank_assignment.tsv
  rank_taxa_summary.tsv
  rank_uc/
  rank_centroids_core/
  rank_centroids_original/
```

## Multi-reference overlap

Run `insert-backbone` for each reference dataset, then compare assignments:

```bash
autotax2 overlap-backbone \
  --assignments \
    ref2_inserted/sequence_rank_assignment.tsv \
    ref3_inserted/sequence_rank_assignment.tsv \
    ref4_inserted/sequence_rank_assignment.tsv \
    hifimeta_inserted/sequence_rank_assignment.tsv \
  --labels ref2,ref3,ref4,hifimeta \
  --out ref_overlap
```

Typical outputs:

```text
ref_overlap/
  taxon_presence_by_source.tsv
  taxon_count_by_source.tsv
  source_pairwise_overlap_by_rank.tsv
  source_unique_taxa_by_rank.tsv
```

## Classic helper workflow

Dereplicate:

```bash
autotax2 derep \
  --input input.fa \
  --out work \
  --sort \
  --threads auto
```

Cluster:

```bash
autotax2 cluster \
  --input work/derep_sorted.fa \
  --out clusters \
  --ids 0.99,0.97 \
  --threads auto
```

Classify:

```bash
autotax2 classify \
  --input clusters/otu099_centroids.fa \
  --ref-manifest refdatabases/autotax2_ref_manifest.tsv \
  --out classify \
  --threads auto
```

Or run the helper workflow end to end:

```bash
autotax2 run \
  --input input.fa \
  --ref-manifest refdatabases/autotax2_ref_manifest.tsv \
  --out autotax2_out \
  --ids 0.99,0.97,0.95,0.90 \
  --threads auto
```

## Input requirements

### SILVA FASTA

SILVA FASTA headers should contain an accession followed by semicolon-separated taxonomy:

```text
>AB000001.1.1500 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus acidophilus;
ACGT...
```

### SILVA metadata

The metadata table must contain at least:

```text
acc
flags
```

Accessions whose `flags` column contains `[T]` or `[t]` are extracted as type-strain references.

### User FASTA

The first token of each header is treated as the sequence ID:

```text
>seq000001 optional description
ACGT...
```

## Development

Install in editable mode and run tests:

```bash
python -m pip install -e .
python -m pip install pytest
pytest
```

Check syntax quickly:

```bash
python -m py_compile autotax2/cli.py autotax2/prepare.py
```
