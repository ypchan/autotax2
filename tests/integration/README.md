# Optional Real-Tool Integration Tests

These tests are optional. They require real external tools and real input files:

- SINA
- VSEARCH
- a SILVA FASTA file or a small SILVA-compatible fixture
- a dataset FASTA that already contains extracted SSU/16S sequences

Normal `pytest` skips these tests. They run only when:

```bash
AUTOTAX2_RUN_INTEGRATION=1
```

The Python integration test calls `scripts/run_real_integration_test.sh`.
Provide paths through environment variables:

```bash
AUTOTAX2_RUN_INTEGRATION=1 \
AUTOTAX2_INTEGRATION_SILVA_FASTA=/home/database/silva_138.2/silva_rep/SILVA_138.2_SSURef_NR99_tax_silva.dna.fasta.gz \
AUTOTAX2_INTEGRATION_DATASET_FASTA=/path/to/example.ssu.fa \
AUTOTAX2_INTEGRATION_OUTDIR=/tmp/autotax2_real_test \
AUTOTAX2_INTEGRATION_DOMAIN=Archaea \
AUTOTAX2_INTEGRATION_DATASET_NAME=digester2020_test \
AUTOTAX2_INTEGRATION_PREFIX=TST \
AUTOTAX2_INTEGRATION_THREADS=8 \
python -m pytest tests/integration -m integration
```

You can also run the shell script directly:

```bash
bash scripts/run_real_integration_test.sh \
  --silva-fasta /home/database/silva_138.2/silva_rep/SILVA_138.2_SSURef_NR99_tax_silva.dna.fasta.gz \
  --dataset-fasta /path/to/example.ssu.fa \
  --outdir /tmp/autotax2_real_test \
  --domain Archaea \
  --dataset-name digester2020_test \
  --prefix TST \
  --threads 8
```

Do not commit large real data files. Keep real SILVA files and real datasets
outside the repository.
