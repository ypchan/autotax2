# Optional Real-Tool Integration Tests

These tests are optional. Normal `pytest` skips them. They require real tools
and real input files:

- SINA
- VSEARCH
- official SILVA NR99 taxonomy FASTA
- official SILVA `full_metadata`
- GTDB r232 `ar53_taxonomy_r232.tsv`
- GTDB r232 `bac120_taxonomy_r232.tsv`
- a dataset FASTA that already contains extracted SSU/16S sequences

Run through pytest:

```bash
AUTOTAX2_RUN_INTEGRATION=1 \
AUTOTAX2_INTEGRATION_SILVA_FASTA=/db/silva/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
AUTOTAX2_INTEGRATION_TYPE_STRAIN_METADATA=/db/silva/SILVA_138.2_SSURef_Nr99.full_metadata.gz \
AUTOTAX2_INTEGRATION_GTDB_AR53_TAXONOMY=/db/gtdb/ar53_taxonomy_r232.tsv \
AUTOTAX2_INTEGRATION_GTDB_BAC120_TAXONOMY=/db/gtdb/bac120_taxonomy_r232.tsv \
AUTOTAX2_INTEGRATION_DATASET_FASTA=/path/to/example.ssu.fa \
AUTOTAX2_INTEGRATION_OUTDIR=/tmp/autotax2_real_test \
AUTOTAX2_INTEGRATION_DOMAIN=Archaea \
AUTOTAX2_INTEGRATION_DATASET_NAME=digester2020_test \
AUTOTAX2_INTEGRATION_PREFIX=TST \
AUTOTAX2_INTEGRATION_THREADS=8 \
python -m pytest tests/integration -m integration
```

Optional environment variables:

```bash
AUTOTAX2_INTEGRATION_SINA_BIN=/path/to/sina
AUTOTAX2_INTEGRATION_VSEARCH_BIN=/path/to/vsearch
AUTOTAX2_INTEGRATION_SINA_REFERENCE=/path/to/sina_reference.ptdb
AUTOTAX2_INTEGRATION_SINA_SEARCH_DB=/path/to/sina_search_db
AUTOTAX2_INTEGRATION_SEARCH_CANDIDATES=0
AUTOTAX2_INTEGRATION_REQUIRE_SINA_CANDIDATES=1
AUTOTAX2_INTEGRATION_STRICT_VALIDATE=1
```

Run the shell script directly:

```bash
bash scripts/run_real_integration_test.sh \
  --silva-fasta /db/silva/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata /db/silva/SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --gtdb-ar53-taxonomy /db/gtdb/ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy /db/gtdb/bac120_taxonomy_r232.tsv \
  --dataset-fasta /path/to/example.ssu.fa \
  --outdir /tmp/autotax2_real_test \
  --domain Archaea \
  --dataset-name digester2020_test \
  --prefix TST \
  --threads 8 \
  --search-candidates \
  --search-min-sim 0.5 \
  --search-max-result 10
```

The script runs:

```text
init -> validate -> resolve -> validate -> prepare -> orient -> cluster -> place -> summarize -> export -> validate
```

It checks export format compatibility plus SILVA resolve evidence, placement
evidence, and SINA candidate diagnostics when candidate search is enabled.

Do not commit large real data files. Keep real SILVA, GTDB, and dataset files
outside the repository.
