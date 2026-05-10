# autotax2 Demo Workflow

This is the full command sequence for a typical archaeal build. Replace file
names and dataset metadata with the real project inputs.

`digester2020.ssu.fa` is expected to contain already extracted SSU/16S
sequences. Run any locus extraction step before this workflow; autotax2 starts
from the extracted FASTA.

```bash
autotax2 init \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.dna.fasta.gz \
  --outdir autotax2_build \
  --domain Archaea

autotax2 resolve \
  --build autotax2_build \
  --threads 48

autotax2 prepare \
  --build autotax2_build \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea

autotax2 orient \
  --build autotax2_build \
  --dataset digester2020 \
  --threads 48

autotax2 cluster \
  --build autotax2_build \
  --dataset digester2020 \
  --threads 48

autotax2 place \
  --build autotax2_build \
  --dataset digester2020

autotax2 export all \
  --build autotax2_build \
  --gzip

autotax2 summarize \
  --build autotax2_build

autotax2 validate \
  --build autotax2_build
```

Main outputs:

- `registry/`: durable build state.
- `silva/`: SILVA named backbone and unresolved framework outputs.
- `datasets/01_digester2020/`: dataset-specific normalized FASTA, prepared
  SSU/16S FASTA, SINA, VSEARCH, placement, and summary files.
- `export/`: SINTAX, QIIME2, and DADA2 reference files.
- `export/export_validation.tsv`: automatic export format self-check report.
- `reports/`: global summaries and validation reports.
- `logs/`: dated command audit logs such as `export_dateYYYYMMDDHHMMSS.log`.

Run `autotax2 validate --build autotax2_build --strict` before a release to
turn selected reproducibility and export-compatibility warnings into failures.
