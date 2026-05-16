# autotax2 Real-Tool Run Checklist

Date: 2026-05-16

This checklist is for the first real SINA/VSEARCH run after the mocked test
suite passes. It is intentionally operational: follow it before using a large
official SILVA build for publication work.

## 1. Inputs

Required official/backbone inputs:

```text
SILVA NR99 taxonomy FASTA, plain or .gz
SILVA full_metadata, plain or .gz
GTDB r232 ar53_taxonomy_r232.tsv
GTDB r232 bac120_taxonomy_r232.tsv
```

Required dataset input:

```text
dataset.ssu.fa
```

The dataset FASTA must already contain extracted SSU/16S sequences. autotax2
does not run barrnap, nhmmer, BLAST, or any other extractor internally.

## 2. Tool Preflight

Check that real tools are visible:

```bash
sina --version
vsearch --version
gzip --version
autotax2 --help
```

If a tool is not on `PATH`, pass it explicitly:

```bash
--sina-bin /path/to/sina
--vsearch-bin /path/to/vsearch
```

## 3. First Smoke Run

Start with a small dataset FASTA. Use the optional integration script so every
step and output check stays reproducible:

```bash
bash scripts/run_real_integration_test.sh \
  --silva-fasta /db/silva/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata /db/silva/SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --gtdb-ar53-taxonomy /db/gtdb/ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy /db/gtdb/bac120_taxonomy_r232.tsv \
  --dataset-fasta /path/to/smoke_dataset.ssu.fa \
  --outdir /tmp/autotax2_smoke \
  --domain Archaea \
  --dataset-name smoke_archaea \
  --prefix SMK \
  --threads 8 \
  --search-candidates \
  --search-min-sim 0.5 \
  --search-max-result 10
```

If SINA candidate search needs a separate database in your install, add:

```bash
--sina-search-db /path/to/sina_search_db
```

If the goal is only to verify the main registry flow first, disable candidate
search:

```bash
--no-search-candidates
```

## 4. Command Order

The script runs this sequence:

```text
init
validate --no-check-exports
resolve
validate --no-check-exports
prepare
orient
cluster
place
summarize --overwrite
export all --gzip
validate
```

Use `--strict-validate` on the script once the smoke run has no warnings that
you intentionally accept.

## 5. Critical Output Checks

After `init`:

```text
silva/silva_named_backbone.fa
silva/silva_unresolved.fa
silva/silva_rejected.tsv
registry/legal_name_catalog.tsv
registry/protected_taxa_snapshot.tsv
logs/init_start_date*.log
logs/init_date*.log
```

Check:

- `silva_rejected.tsv` is not dominated by unexpected `invalid_sequence_characters`.
- `legal_name_catalog.tsv` contains GTDB r232 names/placeholders after prefix removal.
- organelle lineages in `silva_named_backbone.tax.tsv` are standardized to seven ranks.

After `resolve`:

```text
silva/silva_unresolved_taxa.tsv
silva/silva_unresolved_members.tsv
silva/silva_unresolved_mapping.tsv
silva/silva_unresolved_evidence.tsv
silva/silva_unresolved_clusters/<rank>/*.uc
```

Check:

- every unresolved member has evidence rows from phylum through species;
- ambiguous ranks are explicit, not hidden as new placeholders;
- parent-local UC files are grouped by rank and parent hash;
- placeholder IDs use the `SILVA` prefix and are never reused.

After `orient`:

```text
datasets/NN_name/sina.oriented.fa
datasets/NN_name/sina.summary.tsv
datasets/NN_name/sina.candidates.tsv
datasets/NN_name/sina.log
```

Check:

- `sina.summary.tsv` reports plus/minus/unknown and any fallback;
- if `--search-candidates` was used, `sina.candidates.tsv` exists;
- SINA similarity is only advisory and is not used as a final rank threshold.

After `cluster`:

```text
datasets/NN_name/internal_clusters/
datasets/NN_name/vs_registry.raw.tsv
datasets/NN_name/vs_registry.filtered.tsv
datasets/NN_name/cluster_search_summary.tsv
datasets/NN_name/sina_candidate_diagnostics.tsv
registry/current_representatives.fa
registry/current_representatives.sina_candidates.fa
```

Check:

- `cluster_search_summary.tsv` records `iddef=2`;
- `sina_candidate_*` fields explain whether SINA candidates matched the current registry;
- `sina_candidate_diagnostics.tsv` shows per-query candidate subset or fallback behavior;
- if many sequences later become unplaced, inspect `vs_registry.raw.tsv` before the filtered file.

After `place`:

```text
datasets/NN_name/assignments.tsv
datasets/NN_name/near_best_consensus.tsv
datasets/NN_name/placement_evidence.tsv
datasets/NN_name/created_taxa.tsv
```

Check:

- each assigned sequence has placement evidence at phylum through species;
- species decisions require `>= 0.972`;
- if multiple species-level candidates are `>= 0.972`, inherit the
  highest-identity species lineage rather than creating a new species placeholder;
- if different species tie for highest identity, keep the species call ambiguous;
- genus decisions require `>= 0.901`;
- ambiguous near-best consensus remains explicit.

After `summarize` and `validate`:

```text
reports/global_summary.tsv
reports/dataset_delta_summary.tsv
reports/dataset_overlap_matrix.tsv
reports/rank_novelty_summary.tsv
reports/dataset_rank_overlap_detail.tsv
reports/dataset_rank_novelty_detail.tsv
reports/dataset_increment_audit.md
reports/validation_report.md
```

Check:

- `global_summary.tsv` includes SILVA resolve evidence counts;
- `dataset_delta_summary.tsv` includes placement evidence rows and SINA candidate counts;
- `dataset_overlap_matrix.tsv` answers rank overlap against existing sources;
- `dataset_rank_overlap_detail.tsv` traces each rank overlap/new call to taxon IDs and sequence IDs;
- `dataset_rank_novelty_detail.tsv` lists new placeholder taxa with parent and support;
- `dataset_increment_audit.md` gives the human-readable per-dataset audit summary;
- `validation_report.md` has no unexpected warnings before publication work.

## 6. Common Failure Modes

SINA candidate TSV is empty:

- confirm `--search-candidates` was used;
- inspect `sina.log`;
- pass `--sina-search-db` if your SINA install requires a separate search DB;
- rerun with `--no-search-candidates` to isolate orientation from candidate search.

SINA candidates do not match current representatives:

- inspect `sina_candidate_diagnostics.tsv`;
- check accession formatting, especially version suffixes and pipe-delimited labels;
- use normal fallback first, then try `--require-sina-candidates` only after accession matching is understood.

VSEARCH raw hits exist but filtered hits are empty:

- inspect query and target coverage in `vs_registry.tsv`;
- check `--min-query-cov` and `--min-target-cov`;
- confirm `--iddef 2` is recorded.

Validation fails after `resolve`:

- inspect `silva_unresolved_evidence.tsv`;
- check whether any unresolved member lacks all six rank evidence rows;
- inspect placeholder format and deprecated placeholder registry.

Validation fails after `place`:

- inspect `placement_evidence.tsv`;
- confirm duplicate sequences are accounted for;
- confirm active representatives point to valid sequence IDs.

## 7. Data Hygiene

Do not commit large real SILVA, GTDB, SINA, or dataset files. Keep them outside
the repository or in ignored local paths. Commit only scripts, documentation,
small fixtures, and reproducible tests.
