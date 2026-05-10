# autotax2

autotax2 is a fixed-backbone, rank-aware, incremental rRNA gene reference
builder. It preserves a stable named SILVA backbone while allowing controlled
additions from downstream SSU/16S datasets. The project is designed for
traceable taxonomy extension, reproducible classifier exports, and explicit
validation of registry invariants.

This documentation is intentionally text-only. No workflow diagrams are used.

## 1. Core Guarantees

1. The named SILVA backbone is immutable.
2. Named SILVA taxa with reliable taxonomy are marked as protected.
3. SILVA unresolved records may form a mutable placeholder framework.
4. Custom datasets can extend the framework, but they cannot overwrite named
   SILVA taxa.
5. Each custom dataset must use a frozen prefix, such as `D20`.
6. Input sequence IDs are remapped to internal IDs such as `D20_000001`.
7. Sequence MD5 values are recorded, and exact duplicate sequences are not
   exported repeatedly.
8. VSEARCH uses fixed `--iddef 2` by default.
9. `prepare` expects an externally processed SSU/16S FASTA file.
10. autotax2 does not run internal sequence extraction tools and does not
    record external extraction-tool versions.
11. SINA is used for orientation correction with loose behavior.
12. Placement uses near-best hit consensus, not a single best hit.
13. Exports include SINTAX, QIIME 2, DADA2 toGenus, and DADA2 assignSpecies
    formats.
14. Major CLI runs write dated audit logs under `logs/`.

## 2. Input Contracts

### 2.1 SILVA Input

Recommended SILVA inputs:

```text
SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
SILVA_138.2_SSURef_Nr99.full_metadata.gz
```

`SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz` provides the sequence records and
SILVA taxonomy strings. `SILVA_138.2_SSURef_Nr99.full_metadata.gz` is used to
identify type-strain or type-material evidence. When such evidence is available,
named SILVA representatives derived from type material are treated as the most
reliable representative source.

### 2.2 Custom Dataset Input

`autotax2 prepare` accepts only FASTA files that already contain
externally extracted SSU/16S sequences.

If the original data are full-length 16S records, contigs, MAG frameworks, or
records containing non-target flanking sequence, extract the target rRNA region
before running autotax2. You may use any external preprocessing workflow. For
example, barrnap can be used before `prepare`, but autotax2 does not
invoke barrnap, validate barrnap versions, or produce barrnap output files.

The FASTA file passed to `prepare` should look like this:

```text
>sample_001
ACGT...
>sample_002
ACGT...
```

Sequences should already represent the target SSU/16S locus. Orientation may be
handled later by `orient`.

## 3. Installation

After downloading the repository:

```bash
git clone <your-autotax2-repo-url>
cd autotax2
python -m pip install -e 
autotax2 --help
```


Run the test suite:

```bash
python -m pytest
```

## 4. First Complete Build

The following example builds an archaeal reference.

### 4.1 Initialize the SILVA Backbone

```bash
autotax2 init \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --outdir autotax2_build \
  --domain Archaea
```

This creates `autotax2_build/` and separates SILVA records into named,
unresolved, and rejected categories.

### 4.2 Resolve the SILVA unresolved framework

```bash
autotax2 resolve \
  --build autotax2_build \
  --threads 48
```

This resolves non-domain SILVA unresolved records into SILVA placeholder taxa.

### 4.3 Prepare a Custom Dataset

Assume that external preprocessing has already produced:

```text
digester2020.ssu.fa
```

Run:

```bash
autotax2 prepare \
  --build autotax2_build \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea
```

The prefix `D20` is frozen for this dataset. New placeholders created from this
dataset use names such as:

```text
g__D20g000001
s__D20s000001
```

### 4.4 Correct Orientation with SINA

```bash
autotax2 orient \
  --build autotax2_build \
  --dataset digester2020 \
  --threads 48
```

### 4.5 Cluster and Search with VSEARCH

```bash
autotax2 cluster \
  --build autotax2_build \
  --dataset digester2020 \
  --threads 48
```

### 4.6 Place Dataset Representatives

```bash
autotax2 place \
  --build autotax2_build \
  --dataset digester2020
```

### 4.7 Export Classifier References

```bash
autotax2 export all \
  --build autotax2_build \
  --gzip
```

### 4.8 Summarize and Validate

```bash
autotax2 summarize \
  --build autotax2_build

autotax2 validate \
  --build autotax2_build
```

Before release, use strict validation:

```bash
autotax2 validate \
  --build autotax2_build \
  --strict
```

## 5. Adding a New Full-Length 16S Dataset

Assume the original file is:

```text
new_full_length_16s.fa
```

Recommended workflow:

1. Confirm outside autotax2 that the records represent the target SSU/16S
   region.
2. If the file still contains full-length, contig, or flanking sequence, extract
   the SSU/16S region with an external workflow.
3. Produce `new_dataset.ssu.fa`.
4. Run the autotax2 dataset workflow.

If you choose to use barrnap, it belongs in step 2, before `prepare`.
From `prepare` onward, autotax2 accepts only the processed SSU/16S FASTA.

```bash
autotax2 prepare \
  --build autotax2_build \
  --name new_dataset \
  --prefix ND1 \
  --fasta new_dataset.ssu.fa \
  --domain Bacteria

autotax2 orient \
  --build autotax2_build \
  --dataset new_dataset \
  --threads 48

autotax2 cluster \
  --build autotax2_build \
  --dataset new_dataset \
  --threads 48

autotax2 place \
  --build autotax2_build \
  --dataset new_dataset
```

Then regenerate exports and reports:

```bash
autotax2 export all \
  --build autotax2_build \
  --gzip

autotax2 summarize \
  --build autotax2_build

autotax2 validate \
  --build autotax2_build
```

## 6. Adding Multiple Datasets Incrementally

Assume three externally processed SSU/16S FASTA files:

```text
dataset_a.ssu.fa
dataset_b.ssu.fa
dataset_c.ssu.fa
```

Assign one unique prefix per dataset:

```text
dataset_a -> A01
dataset_b -> B01
dataset_c -> C01
```

Add datasets sequentially. Each dataset should complete `prepare ->
orient -> cluster -> place` before the next dataset begins, because
later placement runs use the updated registry produced by earlier datasets.

```bash
autotax2 prepare --build autotax2_build --name dataset_a --prefix A01 --fasta dataset_a.ssu.fa --domain Archaea
autotax2 orient     --build autotax2_build --dataset dataset_a --threads 48
autotax2 cluster  --build autotax2_build --dataset dataset_a --threads 48
autotax2 place           --build autotax2_build --dataset dataset_a

autotax2 prepare --build autotax2_build --name dataset_b --prefix B01 --fasta dataset_b.ssu.fa --domain Archaea
autotax2 orient     --build autotax2_build --dataset dataset_b --threads 48
autotax2 cluster  --build autotax2_build --dataset dataset_b --threads 48
autotax2 place           --build autotax2_build --dataset dataset_b

autotax2 prepare --build autotax2_build --name dataset_c --prefix C01 --fasta dataset_c.ssu.fa --domain Archaea
autotax2 orient     --build autotax2_build --dataset dataset_c --threads 48
autotax2 cluster  --build autotax2_build --dataset dataset_c --threads 48
autotax2 place           --build autotax2_build --dataset dataset_c
```

After all additions:

```bash
autotax2 export all --build autotax2_build --gzip
autotax2 summarize --build autotax2_build
autotax2 validate --build autotax2_build --strict
```

## 7. Handling SILVA Records That Become Unidentified at a Rank

autotax2 parses SILVA taxonomy into seven ranks:

```text
domain;phylum;class;order;family;genus;species
```

Classification rules:

1. If `domain` is empty, missing, or contains an unresolved token such as
   `unidentified`, the SILVA record is rejected.
2. Rejected records are excluded from both the named backbone and unresolved
   framework.
3. Rejected records are written to:

```text
silva/silva_rejected.tsv
```

4. If `domain` is reliable but a lower rank contains unresolved terms such as
   `unidentified`, `uncultured`, `unknown`, `environmental`, `metagenome`, or
   `sp.`, the first unresolved rank and all lower ranks are treated as
   unresolved.
5. The record is written to:

```text
silva/silva_unresolved.fa
silva/silva_unresolved.tsv
```

6. `resolve` creates SILVA placeholders for supported unresolved ranks.

Example:

```text
Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;unidentified;unidentified
```

Result:

```text
lowest_reliable_rank = family
unresolved_ranks = genus,species
genus placeholder = g__SILVAg000001
species placeholder = s__SILVAs000001
```

If unresolved taxonomy begins at class, order, or family, all lower ranks are
also treated as unresolved. The current resolver primarily implements genus and
species clustering. Family, order, class, and floor threshold options are
reserved and must remain at their defaults; changing them raises an error so
users do not mistake reserved controls for completed behavior.

## 8. Type-Strain Metadata

When `SILVA_138.2_SSURef_Nr99.full_metadata.gz` is provided, autotax2 attempts
to identify type-strain or type-material evidence.

Algorithm:

```text
For each named SILVA species:
  collect all named SILVA sequences assigned to that species
  if one or more sequences have type-strain or type-material evidence:
      choose a type-evidence sequence as representative
      representative_reason = type_strain
      is_type_strain = true
  else:
      choose the first named SILVA sequence for that species
      representative_reason = first_silva_named_for_species
```

Implications:

1. Named SILVA taxonomy remains immutable.
2. Type-evidence named SILVA representatives have the highest representative
   priority.
3. Custom datasets cannot overwrite protected named SILVA taxa.
4. Exports preferentially use the most reliable available representative.

## 9. Command Reference

### 9.1 `autotax2 init`

Purpose: initialize the SILVA backbone.

Example:

```bash
autotax2 init \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --outdir autotax2_build \
  --domain Archaea
```

Parameters:

```text
--silva-fasta            SILVA taxonomy FASTA, plain or gzipped
--type-strain-metadata   Optional SILVA full_metadata file or autotax2 TSV
--outdir                 Output build directory
--domain                 Optional domain filter, for example Archaea or Bacteria
```

Steps:

```text
1. Read the SILVA FASTA.
2. Parse seq_id and seven-rank taxonomy from each header.
3. Reject records with empty or unresolved domain values into silva_rejected.tsv.
4. Apply --domain filtering when requested.
5. Classify each retained record as named or unresolved.
6. Write complete reliable records into the named backbone.
7. Write non-domain unresolved records into the unresolved framework.
8. Read type-strain metadata when provided.
9. Build protected taxon_nodes, sequence_registry, and representative_registry.
10. Initialize placeholder counters and the protected SILVA snapshot.
```

Primary outputs:

```text
silva/silva_named_backbone.fa
silva/silva_named_backbone.tax.tsv
silva/silva_unresolved.fa
silva/silva_unresolved.tsv
silva/silva_rejected.tsv
registry/taxon_nodes.tsv
registry/sequence_registry.tsv
registry/representative_registry.tsv
registry/placeholder_counters.tsv
registry/placeholder_counters.yaml
registry/protected_taxa_snapshot.tsv
```

### 9.2 `autotax2 resolve`

Purpose: resolve SILVA unresolved records into a SILVA placeholder framework.

Example:

```bash
autotax2 resolve \
  --build autotax2_build \
  --threads 48
```

Parameters:

```text
--build        Build directory
--threads      VSEARCH thread count
--species-id   Species threshold, default 0.987
--genus-id     Genus threshold, default 0.945
--vsearch-bin  VSEARCH executable
--iddef        VSEARCH --iddef value, default 2
--dry-run      Write proposals without updating the registry
```

Reserved parameters:

```text
--family-id
--order-id
--class-id
--floor-id
```

These parameters must currently remain at their defaults.

Algorithm:

```text
Input:
  silva_unresolved.fa
  silva_unresolved.tsv

For each unresolved SILVA record:
  read lowest_reliable_rank
  read unresolved_ranks

Cluster unresolved records:
  species clusters at species_id
  genus clusters at genus_id

For each unresolved rank:
  if rank is genus:
      create or reuse g__SILVAgNNNNNN by genus cluster
  if rank is species:
      create or reuse s__SILVAsNNNNNN by species cluster within genus context
  if rank is class, order, or family:
      create or reuse a placeholder using the current unresolved cluster context

Never reuse deprecated placeholder IDs.
Never mutate named SILVA nodes.
```

Primary outputs:

```text
silva/silva_unresolved_taxa.tsv
silva/silva_unresolved_members.tsv
silva/silva_unresolved_mapping.tsv
silva/silva_unresolved.resolved.fa
registry/taxon_nodes.tsv
registry/name_index.tsv
registry/cluster_to_taxon.tsv
```

### 9.3 `autotax2 prepare`

Purpose: register a custom dataset, freeze its prefix, remap sequence IDs, and
record MD5-based de-duplication membership.

Example:

```bash
autotax2 prepare \
  --build autotax2_build \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea
```

Parameters:

```text
--build                 Build directory
--name                  Dataset name
--prefix                Frozen dataset prefix, for example D20
--fasta                 Externally processed SSU/16S FASTA
--domain                Archaea or Bacteria
--min-ssu-len-archaea   Minimum archaeal SSU/16S length, default 900
--min-ssu-len-bacteria  Minimum bacterial SSU/16S length, default 1200
--reject-non-atgc       Reject sequences containing non-ATGC symbols, enabled by default
--no-reject-non-atgc    Allow non-ATGC symbols into the normalized input
```

Steps:

```text
1. Validate that the dataset name is not empty.
2. Validate that the prefix is legal and not assigned to another dataset.
3. Record the input FASTA MD5.
4. Read the input FASTA.
5. Remap original IDs to prefix_000001, prefix_000002, and so on.
6. Convert sequence letters to uppercase and replace U with T.
7. Mark non_atgc records when reject_non_atgc is enabled.
8. Mark short_ssu records when length is below the domain-specific threshold.
9. Write input.normalized.fa.
10. Write accepted records to prepared.ssu.fa.
11. Compute sequence MD5 values.
12. Write unique_sequence_registry.tsv and sequence_membership.tsv.
13. Write prepare_summary.tsv.
```

Important behavior:

```text
autotax2 does not run external sequence extraction tools.
autotax2 does not write extraction-tool GFF3 files.
autotax2 does not check external extraction-tool versions.
prepared.ssu.fa is the input for downstream SINA and VSEARCH steps.
```

Primary outputs:

```text
datasets/01_<dataset>/input.normalized.fa
datasets/01_<dataset>/prepared.ssu.fa
datasets/01_<dataset>/sequence_id_map.tsv
datasets/01_<dataset>/unique_sequence_registry.tsv
datasets/01_<dataset>/sequence_membership.tsv
datasets/01_<dataset>/prepare_summary.tsv
```

### 9.4 `autotax2 orient`

Purpose: correct prepared SSU/16S sequence orientation with SINA.

Example:

```bash
autotax2 orient \
  --build autotax2_build \
  --dataset digester2020 \
  --threads 48
```

Parameters:

```text
--build                       Build directory
--dataset                     Prepared dataset name
--threads                     SINA thread count
--sina-bin                    SINA executable
--reference                   Optional SINA reference/PTDB path
--strict-tool-version         Fail if SINA version cannot be detected
--allow-sina-failure          Allow fallback behavior when SINA fails
--no-allow-sina-failure       Fail immediately when SINA fails
--fallback-copy-original      Copy prepared FASTA records when SINA output is missing
--no-fallback-copy-original   Disable fallback copying
```

Steps:

```text
1. Read prepared.ssu.fa.
2. Run SINA.
3. Read sina.oriented.fa.
4. Compare original and SINA output records to infer plus, minus, or unknown strand.
5. Copy prepared.ssu.fa as fallback output when SINA fails and fallback is allowed.
6. Write sina.summary.tsv.
```

Primary outputs:

```text
datasets/01_<dataset>/sina.oriented.fa
datasets/01_<dataset>/sina.summary.tsv
datasets/01_<dataset>/tool_versions.tsv
```

### 9.5 `autotax2 cluster`

Purpose: cluster the dataset internally and search dataset representatives
against the current registry representatives.

Example:

```bash
autotax2 cluster \
  --build autotax2_build \
  --dataset digester2020 \
  --threads 48
```

Parameters:

```text
--build                 Build directory
--dataset               Prepared and oriented dataset name
--threads               VSEARCH thread count
--vsearch-bin           VSEARCH executable
--strict-tool-version   Fail if VSEARCH version cannot be detected
--iddef                 VSEARCH --iddef value, default 2
--species-id            Species threshold, default 0.987
--genus-id              Genus threshold, default 0.945
--family-id             Family threshold, default 0.865
--order-id              Order threshold, default 0.820
--class-id              Class threshold, default 0.785
--floor-id              Minimum registry search threshold, default 0.750
--min-query-cov         Minimum query coverage, default 0.80
--min-target-cov        Minimum target coverage, default 0.0
--maxaccepts            VSEARCH --maxaccepts value
--maxrejects            VSEARCH --maxrejects value
--near-best-delta       Near-best hit retention delta, default 0.005
--strand                plus or both
```

Algorithm:

```text
Input:
  sina.oriented.fa
  registry/representative_registry.tsv
  registry/sequence_registry.tsv

Build current representatives:
  prefer durable representative registry records
  include named SILVA representatives
  include active placeholder and custom representatives
  avoid stale cached representatives when durable sources are available

Internal clustering:
  run VSEARCH at species, genus, family, order, and class thresholds
  write .uc membership files

Registry search:
  search oriented dataset representatives against current representatives
  keep hits that pass identity and coverage filters
  retain near-best hits for placement consensus
```

Primary outputs:

```text
datasets/01_<dataset>/internal_clusters/
datasets/01_<dataset>/vs_registry.raw.tsv
datasets/01_<dataset>/vs_registry.filtered.tsv
datasets/01_<dataset>/cluster_search_summary.tsv
registry/current_representatives.fa
```

### 9.6 `autotax2 place`

Purpose: place dataset representatives into the current rank-aware registry.

Example:

```bash
autotax2 place \
  --build autotax2_build \
  --dataset digester2020
```

Parameters:

```text
--build                 Build directory
--dataset               clustered dataset name
--near-best-delta       Identity delta for retaining near-best hits, default 0.005
--rank-consensus        Minimum near-best agreement fraction, default 0.80
--species-id            Known-like species threshold, default 0.987
--genus-id              New species threshold, default 0.945
--family-id             New genus threshold, default 0.865
--order-id              New family threshold, default 0.820
--class-id              New order threshold, default 0.785
--floor-id              Minimum placement search threshold, default 0.750
--dry-run               Write proposed placements without updating the registry
--allow-ambiguous       Allow ambiguous records to be written
--no-allow-ambiguous    Fail when a record is ambiguous
```

Algorithm:

```text
For each query representative:
  if an exact MD5 duplicate already exists:
      final_status = duplicate
      do not export the duplicate sequence repeatedly
  else:
      best_identity = maximum identity among filtered hits
      near_best_hits = hits where identity >= best_identity - near_best_delta

For each rank:
  count taxon IDs among near_best_hits
  consensus_fraction = top_count / near_best_hit_count
  stable = consensus_fraction >= rank_consensus

Identity status:
  identity >= species_id -> known_like
  identity >= genus_id   -> new_species
  identity >= family_id  -> new_genus
  identity >= order_id   -> new_family
  identity >= class_id   -> new_order
  identity >= floor_id   -> new_class
  otherwise              -> unplaced

Placement decision:
  known_like requires stable species consensus
  new_species requires stable genus consensus
  new_genus requires stable family consensus
  higher novelty uses the nearest stable higher rank
  unstable required consensus may produce an ambiguous final status
```

Primary outputs:

```text
datasets/01_<dataset>/assignments.tsv
datasets/01_<dataset>/created_taxa.tsv
datasets/01_<dataset>/near_best_consensus.tsv
datasets/01_<dataset>/placement_summary.tsv
registry/taxon_nodes.tsv
registry/representative_registry.tsv
registry/cluster_to_taxon.tsv
```

### 9.7 `autotax2 export`

Purpose: export classifier-ready reference files.

Example:

```bash
autotax2 export all \
  --build autotax2_build \
  --gzip
```

Parameters:

```text
format_name             all, sintax, qiime2, or dada2
--build                 Build directory
--outdir                Output directory, default <build>/export
--gzip / --no-gzip      Whether FASTA outputs are gzipped
--representatives-only  Export active representatives only
--all-unique            Export all active unique MD5 sequences
--prefix                Output filename prefix, default autotax2
--force                 Overwrite existing output files
```

Outputs:

```text
export/sintax/autotax2.sintax.fa.gz
export/qiime2/reference_sequences.fasta.gz
export/qiime2/reference_taxonomy.tsv
export/dada2/autotax2_toGenus_trainset.fa.gz
export/dada2/autotax2_assignSpecies.fa.gz
export/export_manifest.tsv
export/export_validation.tsv
```

Format contract:

```text
SINTAX:
  headers contain ;tax=d:...,p:...,c:...,o:...,f:...,g:...,s:...
  tax values do not contain g__ or s__ prefixes

QIIME 2:
  reference_sequences.fasta(.gz)
  reference_taxonomy.tsv with header Feature ID<TAB>Taxon

DADA2 toGenus:
  FASTA headers contain semicolon-delimited taxonomy through genus

DADA2 assignSpecies:
  FASTA headers contain clean genus species names
  headers do not contain g__ or s__ prefixes
```

`export` automatically performs format self-checks after writing output files.
No additional command is required. The report is written to:

```text
export/export_validation.tsv
```

If a newly exported SINTAX, QIIME 2, or DADA2 file fails the format contract,
`export` exits with an error and records the failure in `export_validation.tsv`.

### 9.8 `autotax2 summarize`

Purpose: generate global summaries and dataset-level delta reports.

Example:

```bash
autotax2 summarize \
  --build autotax2_build
```

Parameters:

```text
--build       Build directory
--outdir      Report directory, default <build>/reports
--overwrite   Overwrite existing report files
```

Outputs:

```text
reports/global_summary.tsv
reports/dataset_delta_summary.tsv
reports/dataset_overlap_matrix.tsv
reports/source_contributions.tsv
reports/representative_summary.tsv
reports/deduplication_summary.tsv
```

### 9.9 `autotax2 validate`

Purpose: validate build invariants and export compatibility.

Example:

```bash
autotax2 validate \
  --build autotax2_build \
  --strict
```

Parameters:

```text
--build              Build directory
--strict             Treat selected warnings as failures
--check-exports      Validate existing export files when present
--no-check-exports   Skip export validation
--report             Markdown report path
```

Validation checks:

```text
1. Named SILVA protected taxa remain immutable.
2. Placeholder IDs are unique.
3. Deprecated placeholders are not reused.
4. Dataset prefixes are unique and frozen.
5. Sequence MD5 and duplicate accounting are internally consistent.
6. Representative registry rows point to valid taxa and sequences.
7. Export files satisfy the SINTAX, QIIME 2, and DADA2 format contracts.
8. SINA and VSEARCH metadata are traceable.
9. Existing export files pass the same export self-check logic used by export.
```

Outputs:

```text
reports/validation_report.md
reports/validation_report.tsv
```

### 9.10 `autotax2 add`

`add` is currently a placeholder command, not the recommended workflow entry
point. Use the explicit command sequence:

```text
prepare -> orient -> cluster -> place
```

## 10. Build Directory Layout

Typical build directory:

```text
autotax2_build/
  silva/
  registry/
  datasets/
    01_digester2020/
    02_nextdataset/
  export/
  reports/
  logs/
```

Major CLI commands write dated audit logs, for example:

```text
logs/prepare_date20260510153022.log
logs/place_date20260510153210.log
logs/export_date20260510153503.log
```

Each audit log is a simple `key=value` text file containing the command,
timestamp, dataset, output paths, and key counts.

Key registry files:

```text
registry/taxon_nodes.tsv
registry/sequence_registry.tsv
registry/representative_registry.tsv
registry/cluster_to_taxon.tsv
registry/name_index.tsv
registry/placeholder_counters.tsv
registry/protected_taxa_snapshot.tsv
```

## 11. Frequently Asked Questions

### 11.1 Do externally processed 16S datasets need to be extracted again?

No. autotax2 assumes that `prepare --fasta` already points to target
SSU/16S sequences.

### 11.2 Does autotax2 run barrnap?

No. barrnap or any other extraction workflow may be used before autotax2, but
autotax2 does not invoke, validate, or record that preprocessing step.

### 11.3 Why are SILVA records with empty domains rejected?

The domain rank is the root of the rank-aware backbone. Records with empty or
unresolved domains cannot be placed safely under Archaea or Bacteria, so they
are written to `silva_rejected.tsv`.

### 11.4 Can placeholder IDs be reused?

No. Placeholder counters are durable, and deprecated IDs are never reused.

### 11.5 What happens when identical sequences appear in multiple datasets?

Each sequence receives an MD5 digest. Exact duplicate sequences retain
membership records, but the same sequence MD5 is not exported repeatedly.

## 12. Developer Checks

```bash
python -m pytest
python -m compileall autotax2
autotax2 --help
autotax2 prepare --help
```
