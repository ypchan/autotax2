# AutoTax2 usage examples

## Initialize

```bash
autotax2 init \
  --ref-fa gtdb_ssu_tax.fa \
  --ref-tax gtdb_ssu.taxonomy.tsv \
  --ref-arb gtdb_ssu.arb \
  --db autotax2_db \
  --threads 64 \
  --debug
```

## Add MiDAS

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

## Add MFD

```bash
autotax2 add \
  --db autotax2_db \
  --input MFD16S.ssu.fa \
  --source mfd \
  --prefix mfd \
  --threads 64
```

## Export

```bash
autotax2 export --db autotax2_db --format all --out export/
```
