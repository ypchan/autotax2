# autotax2

autotax2 is a fixed-backbone, rank-aware, incremental rRNA gene reference
builder. It preserves a stable named SILVA backbone while allowing controlled
additions from downstream SSU/16S datasets. The project is designed for
traceable taxonomy extension, reproducible classifier exports, and explicit
validation of registry invariants.

This documentation is intentionally text-only. No workflow diagrams are used.

For the upstream AutoTax reading notes and the design translation into
autotax2, see `docs/autotax_upstream_notes.md`.

For the first real SINA/VSEARCH smoke run, see
`docs/real_run_checklist.md`.

## End-to-End Tutorial

This section is the recommended read path for a new user. It starts from a
fresh checkout and ends with classifier-ready exports and validation reports.
The tutorial is intentionally unnumbered; the numbered sections after it are
reference material. The later command reference gives the full parameter details
for each step.

### What autotax2 Does

autotax2 builds an incremental SSU/16S reference in three layers:

```text
official SILVA NR99 named backbone
  -> resolved SILVA placeholder framework
  -> one or more user datasets appended in chronological order
```

The named SILVA backbone is fixed after `init`. Later datasets can reuse named
SILVA taxa, reuse SILVA-derived placeholders, reuse earlier dataset
placeholders, or create new dataset-specific placeholders. They cannot rename,
replace, or mutate protected named SILVA taxa.

### Required Software

Required:

```text
Python >= 3.10
SINA
VSEARCH
```

Check that the command-line tools are visible:

```bash
python --version
vsearch --version
sina --version
```

If `vsearch` or `sina` is not on `PATH`, pass the executable paths explicitly
with `--vsearch-bin` and `--sina-bin`.

### Install autotax2

On a Linux server:

```bash
git clone <your-autotax2-repo-url>
cd autotax2
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
autotax2 --help
python -m pytest
```

### Prepare the Official Input Files

`init` accepts official files directly, including gzip-compressed files.

Required backbone inputs:

```text
SILVA NR99 taxonomy FASTA:
  SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz

SILVA full metadata:
  SILVA_138.2_SSURef_Nr99.full_metadata.gz

GTDB r232 taxonomy tables:
  ar53_taxonomy_r232.tsv
  bac120_taxonomy_r232.tsv
```

The GTDB files are required only during `init`. autotax2 reads their second
column, removes rank prefixes such as `d__`, `p__`, `c__`, `o__`, `f__`, `g__`,
and `s__`, and records accepted names/placeholders in
`registry/legal_name_catalog.tsv`.

Custom dataset input:

```text
dataset.ssu.fa
```

This FASTA must already contain extracted SSU/16S sequences. autotax2 does not
extract rRNA genes from genomes, contigs, MAGs, or amplicon flanking regions.

### Choose a Build Directory and Threads

The build directory is the persistent database. Every later command reads from
this directory and writes new state back into it.

```bash
BUILD=autotax2_build
THREADS=48
```

After initialization, do not manually edit files under `registry/` unless you
are deliberately repairing a build and understand the registry invariants.

### Initialize the SILVA Backbone

```bash
autotax2 init \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --gtdb-ar53-taxonomy ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy bac120_taxonomy_r232.tsv \
  --outdir "$BUILD" \
  --threads "$THREADS"
```

What this does:

```text
1. Reads SILVA full_metadata and records type-material accessions.
2. Reads SILVA NR99 FASTA.
3. Uppercases sequences and converts U to T.
4. Keeps type-material sequences even when they contain ambiguous bases.
5. Rejects non-type sequences with non-ATGC symbols.
6. Builds the strict legal-name catalog from SILVA type material and GTDB.
7. Keeps organelle rRNA lineages as standardized seven-rank named records.
8. Splits retained SILVA records into named and unresolved records.
9. Writes the protected named SILVA registry snapshot.
10. Initializes placeholder counters.
```

Important outputs:

```text
$BUILD/silva/silva_named_backbone.fa
$BUILD/silva/silva_named_backbone.tax.tsv
$BUILD/silva/silva_unresolved.fa
$BUILD/silva/silva_unresolved.tsv
$BUILD/silva/silva_rejected.tsv
$BUILD/registry/taxon_nodes.tsv
$BUILD/registry/sequence_registry.tsv
$BUILD/registry/legal_name_catalog.tsv
$BUILD/registry/protected_taxa_snapshot.tsv
$BUILD/logs/init_start_date*.log
$BUILD/logs/init_date*.log
```

Inspect after `init`:

```bash
autotax2 validate --build "$BUILD" --no-check-exports
```

Look especially at:

```text
silva/silva_rejected.tsv
  Confirms why records were rejected.

silva/silva_unresolved.tsv
  Shows which rank became unresolved and why.

registry/legal_name_catalog.tsv
  Confirms SILVA type-material names and GTDB placeholders were captured.
```

### Resolve the SILVA Unresolved Framework

```bash
autotax2 resolve \
  --build "$BUILD" \
  --threads "$THREADS" \
  --vsearch-bin vsearch \
  --iddef 2
```

Default rank thresholds:

```text
species  >= 0.972
genus    >= 0.901
family   >= 0.801
order    >= 0.729
class    >= 0.722
phylum   >= 0.696
```

What this does:

```text
1. Processes unresolved SILVA records from high rank to low rank.
2. Groups each rank by its resolved immediate parent context.
3. Adds same-parent known anchors to the parent-local VSEARCH job.
4. Runs VSEARCH cluster_fast with the rank threshold and --iddef 2.
5. Assigns an unresolved record to a known same-parent child when exactly one
   anchor taxon supports it.
6. Marks the rank ambiguous when multiple same-parent child anchors support it.
7. Creates SILVA placeholders only when no same-parent anchor supports it.
8. Writes one evidence row per unresolved record per rank.
```

Important outputs:

```text
$BUILD/silva/silva_unresolved_taxa.tsv
$BUILD/silva/silva_unresolved_members.tsv
$BUILD/silva/silva_unresolved_mapping.tsv
$BUILD/silva/silva_unresolved.resolved.fa
$BUILD/silva/silva_unresolved_evidence.tsv
$BUILD/silva/silva_unresolved_clusters/<rank>/*.fa
$BUILD/silva/silva_unresolved_clusters/<rank>/*.uc
$BUILD/registry/taxon_nodes.tsv
$BUILD/registry/cluster_to_taxon.tsv
$BUILD/logs/resolve_date*.log
```

Validate again:

```bash
autotax2 validate --build "$BUILD" --no-check-exports
```

### Add One Dataset with the Recommended Single Command

Use a short, stable prefix for each dataset. The prefix becomes part of internal
sequence IDs and placeholder IDs. It cannot be changed later for the same
dataset.

Example:

```bash
autotax2 add \
  --build "$BUILD" \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea \
  --threads "$THREADS" \
  --search-candidates \
  --search-min-sim 0.5 \
  --search-max-result 10 \
  --iddef 2 \
  --near-best-delta 0.005 \
  --rank-consensus 0.80 \
  --export \
  --force-export
```

`add` runs:

```text
prepare -> orient -> cluster -> place -> summarize -> validate -> optional export
```

Default behavior:

```text
prepare:
  remaps input IDs to D20_000001, D20_000002, ...
  records sequence MD5
  rejects short or invalid sequences

orient:
  uses SINA loose orientation
  optionally writes SINA candidate hits for later VSEARCH rescoring

cluster:
  clusters dataset representatives with VSEARCH
  searches current registry representatives with VSEARCH
  uses SINA candidates to reduce the search database when available

place:
  uses near-best hit consensus
  creates new placeholders only when identity and parent-rank evidence support it

summarize:
  writes global and per-dataset reports

validate:
  checks registry invariants and evidence completeness

export:
  writes SINTAX, QIIME 2, and DADA2 reference files
```

Important dataset outputs:

```text
$BUILD/datasets/01_digester2020/prepared.ssu.fa
$BUILD/datasets/01_digester2020/sequence_id_map.tsv
$BUILD/datasets/01_digester2020/sina.oriented.fa
$BUILD/datasets/01_digester2020/sina.candidates.tsv
$BUILD/datasets/01_digester2020/sina_candidate_diagnostics.tsv
$BUILD/datasets/01_digester2020/internal_clusters/species_0.972.uc
$BUILD/datasets/01_digester2020/vs_registry.filtered.tsv
$BUILD/datasets/01_digester2020/near_best_consensus.tsv
$BUILD/datasets/01_digester2020/assignments.tsv
$BUILD/datasets/01_digester2020/created_taxa.tsv
$BUILD/datasets/01_digester2020/placement_evidence.tsv
```

Important reports and exports:

```text
$BUILD/reports/global_summary.tsv
$BUILD/reports/dataset_delta_summary.tsv
$BUILD/reports/dataset_rank_overlap_detail.tsv
$BUILD/reports/dataset_rank_novelty_detail.tsv
$BUILD/reports/dataset_increment_audit.md
$BUILD/reports/validation_report.md
$BUILD/export/
```

### Add More Datasets Incrementally

Use the same build directory and a new frozen prefix:

```bash
autotax2 add \
  --build "$BUILD" \
  --name wetland2021 \
  --prefix W21 \
  --fasta wetland2021.ssu.fa \
  --domain Bacteria \
  --threads "$THREADS" \
  --search-candidates \
  --search-min-sim 0.5 \
  --search-max-result 10 \
  --export \
  --force-export
```

The second dataset sees:

```text
1. protected named SILVA taxa
2. active SILVA placeholders from resolve
3. active placeholders created by earlier datasets
```

Exact duplicate sequences are recorded as membership evidence but are not
exported repeatedly.

### Manual Step-by-Step Mode

Use the manual commands when debugging a dataset:

```bash
autotax2 prepare \
  --build "$BUILD" \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea

autotax2 orient \
  --build "$BUILD" \
  --dataset digester2020 \
  --threads "$THREADS" \
  --search-candidates \
  --search-min-sim 0.5 \
  --search-max-result 10

autotax2 cluster \
  --build "$BUILD" \
  --dataset digester2020 \
  --threads "$THREADS" \
  --iddef 2 \
  --near-best-delta 0.005

autotax2 place \
  --build "$BUILD" \
  --dataset digester2020 \
  --near-best-delta 0.005 \
  --rank-consensus 0.80

autotax2 summarize --build "$BUILD" --overwrite
autotax2 export all --build "$BUILD" --force
autotax2 validate --build "$BUILD" --check-exports
```

Manual mode writes the same registry and report files as `add`; it just exposes
every intermediate checkpoint.

### What to Check When a Run Finishes

For every run, inspect these files first:

```text
reports/validation_report.md
  Start here. It tells whether protected SILVA taxa, placeholders, duplicates,
  evidence tables, and exports are internally consistent.

reports/dataset_increment_audit.md
  Human-readable summary of what each dataset added or reused.

datasets/<dataset>/assignments.tsv
  One row per representative sequence. This is the final placement table.

datasets/<dataset>/placement_evidence.tsv
  One evidence row per rank per representative. This is the main debug table.

datasets/<dataset>/created_taxa.tsv
  New placeholders created by the dataset.

datasets/<dataset>/sina_candidate_diagnostics.tsv
  Whether SINA candidates were found and whether they matched current registry
  representatives before VSEARCH rescoring.
```

### Recommended Defaults

Use these defaults unless you have a documented reason to change them:

```text
--iddef 2
--search-min-sim 0.5
--search-max-result 10
--search-kmer-candidates 1000
--near-best-delta 0.005
--rank-consensus 0.80
--min-query-cov 0.80
--min-target-cov 0.0
--strand plus
```

The most important scientific thresholds are:

```text
species 0.972
genus   0.901
family  0.801
order   0.729
class   0.722
phylum  0.696
```

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
ar53_taxonomy_r232.tsv
bac120_taxonomy_r232.tsv
```

`SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz` provides the sequence records and
SILVA taxonomy strings. `SILVA_138.2_SSURef_Nr99.full_metadata.gz` is used to
identify type-strain or type-material evidence. When such evidence is available,
named SILVA representatives derived from type material are treated as the most
reliable representative source.

The two GTDB r232 taxonomy files are required during `init`. autotax2 reads the
second column only, strips GTDB rank prefixes such as `p__`, `g__`, and `s__`,
and stores the resulting accepted taxa names/placeholders in
`registry/legal_name_catalog.tsv`. Placeholder-like names that could collide
with autotax2 IDs, such as `SILVAg000001` or `D20s000001`, are skipped.

Organelle rRNA lineages are treated as named SILVA backbone records outside the
strict cellular legal-name filter. When autotax2 sees `Chloroplast`,
`Mitochondria`/`Mitochondrion`, `Plastid`, `Apicoplast`, or `Cyanelle`, it keeps
the upstream lineage and standardizes that rank and every lower rank to the
organelle label so the output is always seven ranks.

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
python -m pip install -e .
autotax2 --help
```


Run the test suite:

```bash
python -m pytest
```

## 4. First Complete Build

The following example initializes the complete SILVA backbone, then adds an
archaeal dataset.

### 4.1 Initialize the SILVA Backbone

```bash
autotax2 init \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --gtdb-ar53-taxonomy ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy bac120_taxonomy_r232.tsv \
  --outdir autotax2_build \
  --threads 48
```

This creates `autotax2_build/` and separates SILVA records into named,
unresolved, and rejected categories.

Full SILVA gzip inputs are large. `init` writes
`logs/init_start_dateYYYYMMDDHHMMSS.log` immediately, then writes
`logs/init_dateYYYYMMDDHHMMSS.log` after successful completion. The terminal
also prints the start-log path at launch and the final audit-log path with the
record counts at completion.

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

Candidate environmental framework labels can remain reliable when they occupy
higher ranks. For example, `GW2011` and `AR15` are retained here:

```text
Archaea;Nanoarchaeota;Nanoarchaeia;Woesearchaeales;GW2011;AR15;archaeon
```

Result:

```text
lowest_reliable_rank = genus
unresolved_ranks = species
family = GW2011
genus = AR15
species = archaeon
```

Generic environmental descriptors are not accepted as named taxa at the rank
where they appear. Current unresolved signals include `uncultured`,
`unidentified`, `uncultivated`, `environmental`, `metagenome`, `metagenomic`,
`archaeon`, `bacterium`, `organism`, `clone`, `sample`, standalone
`Candidatus`, `Incertae`, `incertae sedis`, and species-level habitat labels
such as `gut`, `marine`, `sediment`, `soil`, `sludge`, `wastewater`,
`freshwater`, `rumen`, `oral`, `fecal`, `faecal`, `biofilm`, `thermal`, and
`hydrothermal`.

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

If unresolved taxonomy begins at phylum, class, order, or family, all lower
ranks are also treated as unresolved. The resolver handles these ranks with
parent-scoped hierarchical clustering and creates SILVA placeholders only within
the fixed immediate parent.

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
3. Other named SILVA sequences assigned to the same species are still retained
   as species evidence during placement. If any named sequence from a species
   passes the species threshold, that species can be supported even when the
   chosen type-strain representative alone is below the threshold.
4. Custom datasets cannot overwrite protected named SILVA taxa.
5. Exports preferentially use the most reliable available representative.

## 9. Command Reference

### 9.1 `autotax2 init`

Purpose: initialize the SILVA backbone.

Example:

```bash
autotax2 init \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --type-strain-metadata SILVA_138.2_SSURef_Nr99.full_metadata.gz \
  --gtdb-ar53-taxonomy ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy bac120_taxonomy_r232.tsv \
  --outdir autotax2_build \
  --threads 48
```

Parameters:

```text
--silva-fasta              Official SILVA NR99 taxonomy FASTA, plain or gzipped
--type-strain-metadata     Required official SILVA full_metadata TSV/TSV.gz
--gtdb-ar53-taxonomy       Required GTDB r232 ar53 taxonomy TSV/TSV.gz
--gtdb-bac120-taxonomy     Required GTDB r232 bac120 taxonomy TSV/TSV.gz
--outdir                   Output build directory
--threads                  Threads reserved for init and recorded in metadata
```

Detailed initialization flow:

1. Validate required initialization inputs.

   `init` requires the official SILVA NR99 taxonomy FASTA, the official SILVA
   `full_metadata` file, both GTDB r232 taxonomy files, and an output build
   directory. The inputs may be plain text or gzip-compressed where noted.
   GTDB taxonomy is intentionally required only at initialization time; after
   `init`, downstream commands read the parsed build registries.

2. Read SILVA `full_metadata` first.

   The metadata parser expects an accession-like column and at least one
   type-material/type-strain-like column. It keeps only records with explicit
   type evidence. These records are used for two things: retaining type-derived
   sequences even when they contain non-ATGC symbols, and building the
   type-material part of the legal-name catalog.

3. Read and normalize the SILVA FASTA.

   Each FASTA header is parsed as:

   ```text
   seq_id rank1;rank2;rank3;rank4;rank5;rank6;rank7
   ```

   The seven target ranks are domain, phylum, class, order, family, genus, and
   species. If a header has fewer than seven taxonomy fields, missing ranks are
   padded internally as `unresolved_<rank>`. Sequence text is uppercased and
   RNA `U` is converted to DNA `T`.

4. Reject records that cannot enter the build.

   Records with an empty domain or unresolved domain are written to
   `silva/silva_rejected.tsv`. Non-type SILVA sequences containing characters
   outside `A/T/G/C` are also rejected as `invalid_sequence_characters`.
   Type-material sequences are retained even if they contain ambiguity symbols,
   because they carry important naming evidence.

5. Build the strict legal-name catalog.

   The catalog has two sources:

   - SILVA type-material records contribute accepted names for all ranks along
     their lineage.
   - GTDB r232 contributes accepted taxa names and candidate placeholders from
     the second taxonomy column after removing prefixes such as `p__`, `g__`,
     and `s__`.

   Species names from type-material records are canonicalized to the species
   binomial; strain, isolate, clone, culture, subspecies, and serovar tails are
   removed. GTDB species placeholders such as `Lokiarchaeum sp000002` are kept
   as accepted GTDB names. Names that look like autotax2-generated placeholder
   IDs, such as `SILVAg000001`, `SILVAs000001`, or `D20s000001`, are skipped so
   GTDB cannot collide with autotax2 namespaces.

6. Classify retained SILVA records as named or unresolved.

   A record remains named only when all seven ranks are reliable under the
   current rules. A rank becomes unresolved if it is missing, contains generic
   environmental wording, is a standalone ambiguous value such as `Incertae
   Sedis` or `Candidatus`, contains tokens such as `uncultured`,
   `unidentified`, `metagenome`, `archaeon`, or `bacterium`, or is absent from
   the strict legal-name catalog. Once one rank is unresolved, all lower ranks
   are treated as unresolved too.

7. Standardize organelle rRNA lineages.

   Organelle records are kept as named SILVA backbone records even when their
   lineage does not follow cellular naming rules. If the lineage contains
   `Chloroplast`, `Mitochondria`/`Mitochondrion`, `Plastid`, `Apicoplast`, or
   `Cyanelle`, autotax2 keeps the upstream lineage and fills the organelle rank
   and every lower rank with the canonical organelle label.

   Example:

   ```text
   Eukaryota;Viridiplantae;Chloroplast
   -> Eukaryota;Viridiplantae;Chloroplast;Chloroplast;Chloroplast;Chloroplast;Chloroplast
   ```

8. Write named, unresolved, and rejected SILVA outputs.

   Named records are written to `silva/silva_named_backbone.fa` and
   `silva/silva_named_backbone.tax.tsv`. Unresolved-but-retained records are
   written to `silva/silva_unresolved.fa` and `silva/silva_unresolved.tsv` for
   later SILVA placeholder resolution. Rejected records stay in
   `silva/silva_rejected.tsv` and do not enter registries.

9. Build registry tables.

   `init` creates protected taxon nodes for named SILVA records, sequence MD5
   entries for all retained SILVA records, the parsed legal-name catalog,
   name-index rows, dataset metadata, and representative rows. For named SILVA
   species, a type-material sequence is preferred as the active representative;
   otherwise the first named SILVA sequence for that species is used.

10. Freeze initial build state.

    Placeholder counters are initialized for SILVA placeholders, and
    `registry/protected_taxa_snapshot.tsv` records the immutable named SILVA
    backbone. Later validation compares against this snapshot to detect
    accidental changes to protected SILVA taxa.

Key checks before running `resolve`:

- Inspect `silva/silva_rejected.tsv`. Domain-level rejection or many
  `invalid_sequence_characters` rows usually means the SILVA input or metadata
  pairing should be checked before continuing.
- Inspect `silva/silva_unresolved.tsv`. Confirm that `lowest_reliable_rank`,
  `unresolved_ranks`, and `unresolved_reason` match expectations for
  environmental tails such as `uncultured`, `metagenome`, or `archaeon`.
- Inspect `registry/legal_name_catalog.tsv`. Confirm GTDB names were read from
  the second column, rank prefixes were removed, and useful placeholders such
  as Asgard/Loki-style names are present.
- Inspect organelle rows in `silva/silva_named_backbone.tax.tsv`. They should
  be seven-rank named records, not unresolved records.
- Run `autotax2 validate --build <BUILD>` after `init` and before `resolve`.

Runtime logging:

```text
logs/init_start_dateYYYYMMDDHHMMSS.log
  Written immediately before large input parsing starts.
  Contains status=started and the input paths.

logs/init_dateYYYYMMDDHHMMSS.log
  Written only after successful completion.
  Contains status=completed plus records, named, unresolved, and rejected counts.
```

Primary outputs:

```text
silva/silva_named_backbone.fa
  FASTA of named SILVA records with reliable seven-rank taxonomy.

silva/silva_named_backbone.tax.tsv
  Taxonomy table for the named backbone.
  Columns: seq_id, taxonomy_7rank, domain, phylum, class, order, family, genus, species.

silva/silva_unresolved.fa
  FASTA of SILVA records whose domain is reliable but one or more lower ranks are unresolved.

silva/silva_unresolved.tsv
  Classification table for unresolved SILVA records.
  Includes lowest_reliable_rank, unresolved_ranks, and unresolved_reason.

silva/silva_rejected.tsv
  SILVA records excluded from the build because the domain is empty/unresolved
  or because a non-type sequence contains illegal sequence characters.
  Includes reject_reason.

registry/taxon_nodes.tsv
  Protected taxon tree for the named SILVA backbone.
  Named SILVA rows have protected=true and is_silva_named=true.

registry/sequence_registry.tsv
  Sequence-level registry for retained SILVA records.
  Stores sequence MD5, sequence length, named/unresolved flags, taxonomy, and type-strain fields.

registry/dataset_registry.tsv
  Source registry row for the immutable SILVA backbone.

registry/name_index.tsv
  Lookup table from taxon name and rank to taxon_id.

registry/cluster_to_taxon.tsv
  Empty at init; later populated by SILVA resolve and dataset placement.

registry/representative_registry.tsv
  Active representative sequence per named SILVA species.
  Type-strain/type-material evidence is preferred when available.

registry/placeholder_counters.tsv
  Durable next counters for placeholder IDs such as g__SILVAg000001 and s__SILVAs000001.

registry/placeholder_counters.yaml
  YAML copy of the placeholder counters used by registry update code.

registry/protected_taxa_snapshot.tsv
  Snapshot of protected named SILVA taxa used by validation to detect accidental mutation.
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
--species-id   Species threshold, default 0.972
--genus-id     Genus threshold, default 0.901
--family-id    Family threshold, default 0.801
--order-id     Order threshold, default 0.729
--class-id     Class threshold, default 0.722
--phylum-id    Phylum threshold, default 0.696
--vsearch-bin  VSEARCH executable
--iddef        VSEARCH --iddef value, default 2
--dry-run      Write proposals without updating the registry
```

Detailed resolve flow:

1. Read unresolved SILVA records from the initialized build.

   `resolve` reads `silva/silva_unresolved.fa` and
   `silva/silva_unresolved.tsv`. Only sequence IDs present in the unresolved TSV
   are resolved. If there are no unresolved records, autotax2 writes empty
   resolve output tables and exits successfully.

2. Use rank-specific thresholds.

   The strict resolver processes all unresolved ranks from high to low:

   ```text
   phylum  = 0.696
   class   = 0.722
   order   = 0.729
   family  = 0.801
   genus   = 0.901
   species = 0.972
   ```

3. Create or reuse parent-scoped unresolved-cluster UC files.

   The resolver writes cluster assignments under
   `silva/silva_unresolved_clusters/<rank>/`. Each UC file belongs to one
   immediate-parent group at one rank.

   ```text
   silva_unresolved_clusters/class/class_0.722_<parent-hash>.uc
   silva_unresolved_clusters/genus/genus_0.901_<parent-hash>.uc
   silva_unresolved_clusters/species/species_0.972_<parent-hash>.uc
   ```

   If a UC file already exists, autotax2 reuses it. If it does not exist,
   autotax2 writes a parent-local FASTA containing same-parent known anchors
   first, followed by unresolved candidates, then runs VSEARCH `--cluster_fast`
   with `--iddef 2` by default. If VSEARCH is unavailable or the clustering
   command fails, autotax2 writes deterministic singleton UC assignments so the
   run can still produce traceable placeholder proposals. If you later install
   or fix VSEARCH and want real clustering, delete the existing fallback UC
   files and rerun `resolve`.

4. Schedule parent-local jobs without exceeding `--threads`.

   Within a rank, parent-local jobs are independent after the parent taxon is
   fixed. autotax2 assigns each job a dynamic `job_threads` value from its
   candidate count, then runs jobs in batches whose active thread sum does not
   exceed the command-level `--threads`. Placeholder allocation and evidence
   writing are replayed afterward in deterministic parent order, so repeated
   runs keep stable placeholder IDs.

5. Parse parent-local rank cluster assignments.

   Every unresolved SILVA sequence gets one decision row for every rank from
   phylum through species. If a record is missing from a UC file, autotax2
   assigns it to a deterministic singleton cluster rather than dropping it.

6. Load existing placeholder state.

   `resolve` reads existing placeholders from `registry/taxon_nodes.tsv`,
   `registry/cluster_to_taxon.tsv`, and `registry/placeholder_counters.tsv`.
   Active cluster-to-taxon rows are reused on reruns. Deprecated or superseded
   placeholder IDs are recorded and are never allocated again.

7. Restrict anchor competition to the immediate parent.

   For each rank, known child taxa inside the same parent act as anchors.
   Known children outside that parent do not compete. If one same-parent anchor
   clusters with the unresolved sequence, the sequence is assigned to that
   existing child. If multiple same-parent anchors support it, the rank is
   marked ambiguous. If no same-parent anchor supports it, the residual cluster
   creates or reuses a SILVA placeholder.

8. Resolve each unresolved record rank by rank.

   Reliable ranks are carried forward as prefixed taxonomy labels such as
   `d__Archaea`, `p__Euryarchaeota`, or `g__Methanobacterium`. Unresolved ranks
   that can receive placeholders are replaced with SILVA placeholder names.

   Current placeholder behavior:

   - unresolved phylum uses a phylum parent-local cluster and receives a
     `p__SILVApNNNNNN` placeholder;
   - unresolved class uses a class parent-local cluster and receives a
     `c__SILVAcNNNNNN` placeholder;
   - unresolved order uses an order parent-local cluster and receives an
     `o__SILVAoNNNNNN` placeholder;
   - unresolved family uses a family parent-local cluster and receives an
     `f__SILVAfNNNNNN` placeholder;
   - unresolved genus uses a genus parent-local cluster and receives a
     `g__SILVAgNNNNNN` placeholder;
   - unresolved species uses a species parent-local cluster and receives an
     `s__SILVAsNNNNNN` placeholder.

   Species placeholders intentionally include genus context in their cluster
   key. This prevents the same species UC cluster number from merging species
   placeholders across different genus frameworks.

9. Write resolved SILVA outputs.

   Normal mode writes:

   ```text
   silva/silva_unresolved_taxa.tsv
   silva/silva_unresolved_members.tsv
   silva/silva_unresolved_mapping.tsv
   silva/silva_unresolved_evidence.tsv
   silva/silva_unresolved.resolved.fa
   ```

   `silva_unresolved_taxa.tsv` is the list of created or reused placeholder
   taxa. `silva_unresolved_members.tsv` is per-sequence membership, including
   genus/species placeholder names, cluster keys, and warnings.
   `silva_unresolved_mapping.tsv` maps original SILVA taxonomy to resolved
   autotax2 taxonomy. The resolved FASTA uses the resolved prefixed taxonomy in
   each header.
   `silva_unresolved_evidence.tsv` records the per-sequence, per-rank reasoning
   path.

10. Update registries, unless this is a dry run.

   In normal mode, `resolve` appends or updates placeholder rows in
   `registry/taxon_nodes.tsv`, adds names to `registry/name_index.tsv`, records
   stable cluster mappings in `registry/cluster_to_taxon.tsv`, and advances
   `registry/placeholder_counters.tsv` plus
   `registry/placeholder_counters.yaml`.

   With `--dry-run`, output files are written with `.dry_run` suffixes, but
   registries and placeholder counters are not modified.

10. Preserve the named SILVA backbone.

    `resolve` only acts on the SILVA unresolved framework. Named SILVA taxon
    nodes remain protected and should be unchanged before and after the run.
    This should be verified with `autotax2 validate --build <BUILD>`.

Key checks after `resolve`:

- Inspect `silva/silva_unresolved_taxa.tsv` for new placeholder names,
  parent_taxon_id values, cluster keys, representative sequence IDs, member
  counts, and warnings.
- Inspect `silva/silva_unresolved_members.tsv` to confirm each sequence received
  the expected genus/species placeholders.
- Inspect `silva/silva_unresolved_evidence.tsv` for rank-by-rank parent scope,
  anchor count, residual cluster key, decision, output taxon, and reason.
- Inspect `silva/silva_unresolved_mapping.tsv` to confirm no generic SILVA tail
  such as `unidentified`, `uncultured`, or `metagenome` was kept as a
  placeholder name.
- Inspect `registry/cluster_to_taxon.tsv` before rerunning. Active rows make
  placeholder assignment stable across repeated `resolve` runs.
- Inspect `registry/placeholder_counters.tsv` to confirm counters advanced and
  deprecated placeholder IDs were not reused.
- Run `autotax2 validate --build <BUILD>` after `resolve`.

Primary outputs:

```text
silva/silva_unresolved_taxa.tsv
  Placeholder taxa created or reused for unresolved SILVA groups.

silva/silva_unresolved_members.tsv
  Per-sequence resolved taxonomy, cluster keys, placeholder names, and warnings.

silva/silva_unresolved_mapping.tsv
  Original SILVA taxonomy to resolved autotax2 taxonomy mapping.

silva/silva_unresolved.resolved.fa
  FASTA of unresolved SILVA records with resolved prefixed taxonomy headers.

registry/taxon_nodes.tsv
  Updated with protected SILVA placeholder taxa; named SILVA rows remain unchanged.

registry/name_index.tsv
  Updated with placeholder name lookups.

registry/cluster_to_taxon.tsv
  Stable mapping from unresolved cluster keys to placeholder taxa.

registry/placeholder_counters.tsv
registry/placeholder_counters.yaml
  Advanced next ordinals for SILVA placeholder allocation.
```

Strict resolver behavior:

The resolver uses parent-scoped hierarchical logic for SILVA unresolved records.
The same behavior should also be used later when custom datasets add new taxa.

Ranks are processed from high to low:

```text
phylum  within domain parent      threshold 0.696
class   within phylum parent      threshold 0.722
order   within class parent       threshold 0.729
family  within order parent       threshold 0.801
genus   within family parent      threshold 0.901
species within genus parent       threshold 0.972
```

At each rank and for each immediate parent:

1. Keep the parent fixed.

   A class-level decision is made only inside its phylum parent, an order-level
   decision only inside its class parent, and so on. A known class under another
   phylum is not allowed to compete, even if 16S identity is high, because upper
   ranks may reflect genome, morphology, physiology, biochemical evidence, or
   later taxonomic revision rather than only marker-gene identity.

2. Use known children as anchors.

   Records that already have a reliable child rank inside the same parent are
   not used to create new placeholders at that rank. They are kept as
   known-rank anchors. Unresolved records are first compared against these
   same-parent anchors.

3. Assign to an existing child if support is unique inside the parent.

   For example, while resolving class inside phylum `P1`:

   ```text
   query vs P1/Class1 anchor >= 72.2%
   query vs every other known class under P1 < 72.2%
   => assign query to Class1
   ```

   Anchors under another phylum, such as `P2/ClassX`, are ignored for this
   class-level decision.

4. Mark parent-local ambiguity if multiple same-parent anchors pass.

   ```text
   query vs P1/Class1 anchor >= 72.2%
   query vs P1/Class2 anchor >= 72.2%
   => class decision is ambiguous inside P1
   ```

   The resolver must not choose one arbitrarily and must not create a novel
   class placeholder simply to hide the ambiguity.

5. Cluster only the residual unresolved records.

   Only records that fail to match any same-parent known anchor enter
   same-parent residual clustering. Residual clusters create placeholders under
   that parent. This prevents a sequence from receiving a new placeholder when
   it is already close enough to an accepted same-parent child taxon.

Per-sequence inference evidence contract:

Every rank decision for every unresolved SILVA record, and later every custom
dataset placement decision, should be traceable by a machine-readable evidence
table. The same evidence model should be used in both stages so debugging one
sequence follows the same logic everywhere.

Recommended evidence rows:

```text
source_stage            silva_resolve or dataset_place
dataset                 SILVA or custom dataset name
seq_id                  sequence being judged
rank                    phylum/class/order/family/genus/species
parent_rank             immediate parent rank used to limit candidates
parent_taxon            fixed parent taxonomy/taxon_id
input_rank_value        original rank value, if any
input_rank_status       known, unresolved, missing, environmental, illegal, etc.
threshold               rank identity threshold
candidate_scope         same-parent only
anchor_count            known same-parent anchors considered
best_anchor_taxon       best same-parent anchor taxon_id/name
best_anchor_identity    identity to best anchor
passing_anchor_count    number of same-parent anchors >= threshold
competing_anchor_taxa   all passing same-parent anchors when ambiguous
residual_cluster_key    cluster key if no anchor assignment was possible
residual_cluster_size   residual cluster member count
decision                keep_known, assign_existing, create_placeholder, reuse_placeholder, duplicate, ambiguous, unplaced
output_taxon            assigned existing taxon or new placeholder
reason                  short human-readable decision reason
job_size                unresolved candidate count in this parent-local job
job_threads             threads assigned to this parent-local job
```

For SILVA `resolve`, this evidence is written beside the current resolve
outputs as `silva/silva_unresolved_evidence.tsv`. For custom dataset addition,
the same kind of evidence is written under the dataset directory as
`datasets/NN_name/placement_evidence.tsv`. Existing
files such as `silva_unresolved_members.tsv`, `assignments.tsv`, and
`near_best_consensus.tsv` are useful summaries, but the evidence table should
show the actual rank-by-rank reasoning path.

Parent-scoped parallel scheduling:

Parent-scoped resolution creates many independent work units. The global
`--threads` value is the hard upper limit for the whole command, not a value to
blindly pass to every VSEARCH call.

The current scheduler:

1. Treat each same-parent, same-rank task as an independent job once its parent
   is fixed.
2. Process rank dependencies from high to low with a barrier between ranks.
3. Assign threads dynamically by candidate count.

   A small parent group should usually receive one thread. A large parent group
   can receive more threads.

4. Batch jobs so the sum of VSEARCH threads used by concurrently running jobs
   does not exceed command-level `--threads`.

5. Use candidate sequence count as the first scheduling weight.

   ```text
   job_weight = max(1, residual_sequence_count)
   job_threads = clamp(1, max_threads_per_job,
                       round(total_threads * job_weight / active_job_weight_sum))
   ```

6. Prefer keeping all cores busy with several small jobs instead of giving all
   threads to one tiny cluster. Conversely, allow a large residual cluster to
   use more threads when it would otherwise dominate runtime.

This matters for both SILVA `resolve` and later custom dataset addition: after
parent boundaries are fixed, different parent-local jobs are independent and can
be scheduled safely without changing taxonomy decisions.

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
--search-candidates           Ask SINA to write nearest_slv candidate hits
--search-db                   Optional SINA search database
--search-min-sim              Loose SINA candidate similarity, default 0.500
--search-max-result           Maximum SINA candidates per query, default 10
--search-kmer-candidates      SINA k-mer candidate pool size, default 1000
```

Steps:

```text
1. Read prepared.ssu.fa.
2. Run SINA.
3. Read sina.oriented.fa.
4. Compare original and SINA output records to infer plus, minus, or unknown strand.
5. Copy prepared.ssu.fa as fallback output when SINA fails and fallback is allowed.
6. Optionally parse SINA `nearest_slv` search results into `sina.candidates.tsv`.
7. Write sina.summary.tsv.

SINA candidate hits are advisory only. autotax2 does not use SINA similarity as
the rank identity. When `--search-candidates` is used, SINA narrows the candidate
target set, and the later `cluster` step still rescoring with VSEARCH `--iddef
2` decides whether each rank threshold is passed.
```

Primary outputs:

```text
datasets/01_<dataset>/sina.oriented.fa
datasets/01_<dataset>/sina.summary.tsv
datasets/01_<dataset>/sina.candidates.tsv
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
--species-id            Species threshold, default 0.972
--genus-id              Genus threshold, default 0.901
--family-id             Family threshold, default 0.801
--order-id              Order threshold, default 0.729
--class-id              Class threshold, default 0.722
--phylum-id             Phylum threshold and registry search floor, default 0.696
--min-query-cov         Minimum query coverage, default 0.80
--min-target-cov        Minimum target coverage, default 0.0
--maxaccepts            VSEARCH --maxaccepts value
--maxrejects            VSEARCH --maxrejects value
--near-best-delta       Near-best hit retention delta, default 0.005
--strand                plus or both
--sina-candidates       Optional SINA candidate TSV path
--require-sina-candidates
                        Fail if SINA candidates do not match current representatives
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
  run VSEARCH at species, genus, family, order, class, and phylum thresholds
  write .uc membership files

Registry search:
  if sina.candidates.tsv exists, build a candidate subset database
  write per-query SINA candidate diagnostics
  VSEARCH searches the candidate subset when any target matches
  if no SINA candidate target matches, default fallback is the full registry
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
datasets/01_<dataset>/sina_candidate_diagnostics.tsv
registry/current_representatives.fa
registry/current_representatives.sina_candidates.fa
```

`sina_candidate_diagnostics.tsv` is written when SINA candidate rows are
available. It records, per input query, the species centroid, whether the query
is the centroid searched by VSEARCH, candidate target IDs, matched current
representatives, and whether the command used a candidate subset, full-registry
fallback, or fatal no-match mode.

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
--species-id            Known-like species threshold, default 0.972
--genus-id              New species threshold, default 0.901
--family-id             New genus threshold, default 0.801
--order-id              New family threshold, default 0.729
--class-id              New order threshold, default 0.722
--phylum-id             Minimum placement threshold for new class calls, default 0.696
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
  first discard near_best_hits below that rank's identity threshold
  count distinct taxon IDs among the remaining rank-qualified hits
  consensus_fraction = top_taxon_count / distinct_rank_taxon_count
  stable = consensus_fraction >= rank_consensus

Identity status:
  identity >= species_id -> known_like
  identity >= genus_id   -> new_species
  identity >= family_id  -> new_genus
  identity >= order_id   -> new_family
  identity >= class_id   -> new_order
  identity >= phylum_id  -> new_class
  otherwise              -> unplaced

Placement decision:
  known_like requires stable species consensus
  if multiple species are rank-qualified at >= species_id and species consensus
      is not stable, inherit the highest-identity hit's species lineage
      if different species tie for highest identity, final_status = ambiguous
      do not create a new species placeholder for this case
  new_species requires stable genus consensus
  new_genus requires stable family consensus
  higher novelty uses the nearest stable higher rank
  unstable required consensus may produce an ambiguous final status

Evidence:
  placement_evidence.tsv uses the same per-rank evidence columns as
  silva_unresolved_evidence.tsv. For every query representative and every
  phylum-to-species rank, it records the threshold, near-best scope, best-hit
  rank taxon, best-hit identity, consensus support count, final rank decision,
  output taxon, and any residual placeholder cluster key.
```

Primary outputs:

```text
datasets/01_<dataset>/assignments.tsv
datasets/01_<dataset>/created_taxa.tsv
datasets/01_<dataset>/near_best_consensus.tsv
datasets/01_<dataset>/placement_evidence.tsv
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
reports/rank_novelty_summary.tsv
reports/dataset_rank_overlap_detail.tsv
reports/dataset_rank_novelty_detail.tsv
reports/source_contribution.tsv
reports/representative_summary.tsv
reports/sequence_dedup_summary.tsv
reports/dataset_increment_audit.md
```

`dataset_delta_summary.tsv` includes SINA candidate diagnostics
(`sina_candidate_queries`, `sina_candidate_targets`,
`sina_candidate_target_matches`) and `placement_evidence_rows`, so each dataset
addition can be audited from candidate discovery through final rank evidence.

`global_summary.tsv` includes `silva_unresolved_evidence_records` and
`silva_unresolved_evidence_rows`, so the resolved SILVA placeholder framework
can be checked for sequence-level and rank-level evidence coverage.

`dataset_rank_overlap_detail.tsv` is the rank-level audit table. For each
dataset, rank, and assigned taxon, it reports whether the sequences overlap
named SILVA, unresolved SILVA, previous custom datasets, or the current dataset.
It also lists supporting sequence IDs and placement evidence decisions.

`dataset_rank_novelty_detail.tsv` lists every newly created taxon with parent,
cluster key, representative sequence, and supporting sequence IDs.

`dataset_increment_audit.md` is the human-readable increment report. It
summarizes processing counts, placement outcomes, assigned source categories,
rank overlap, exact sequence overlap, created taxa, and trace files for each
added dataset.

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
9. SILVA resolve evidence covers unresolved SILVA members at phylum through species.
10. Placement evidence covers assignment rows at phylum through species.
11. Existing export files pass the same export self-check logic used by export.
```

Outputs:

```text
reports/validation_report.md
reports/validation_report.tsv
```

### 9.10 `autotax2 add`

Purpose: add one externally extracted SSU/16S dataset by orchestrating:

```text
prepare -> orient -> cluster -> place -> summarize -> validate
```

Example:

```bash
autotax2 add \
  --build autotax2_build \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea \
  --threads 48
```

Full example with the recommended scientific defaults shown explicitly:

```bash
autotax2 add \
  --build autotax2_build \
  --name digester2020 \
  --prefix D20 \
  --fasta digester2020.ssu.fa \
  --domain Archaea \
  --threads 48 \
  --min-ssu-len-archaea 900 \
  --min-ssu-len-bacteria 1200 \
  --reject-non-atgc \
  --sina-bin sina \
  --allow-sina-failure \
  --fallback-copy-original \
  --search-candidates \
  --search-min-sim 0.5 \
  --search-max-result 10 \
  --search-kmer-candidates 1000 \
  --vsearch-bin vsearch \
  --iddef 2 \
  --species-id 0.972 \
  --genus-id 0.901 \
  --family-id 0.801 \
  --order-id 0.729 \
  --class-id 0.722 \
  --phylum-id 0.696 \
  --min-query-cov 0.80 \
  --min-target-cov 0.0 \
  --maxaccepts 50 \
  --maxrejects 256 \
  --near-best-delta 0.005 \
  --rank-consensus 0.80 \
  --strand plus \
  --no-require-sina-candidates \
  --allow-ambiguous \
  --export \
  --gzip \
  --force-export \
  --validate
```

Parameters:

```text
--build                     Initialized and resolved autotax2 build directory
--name                      New dataset name; used as the dataset folder label
--prefix                    Frozen dataset prefix, for example D20 or W21
--fasta                     Input FASTA already containing extracted SSU/16S
--domain                    Archaea or Bacteria
--threads                   Threads passed to SINA and VSEARCH stages

--min-ssu-len-archaea       Minimum archaeal sequence length, default 900
--min-ssu-len-bacteria      Minimum bacterial sequence length, default 1200
--reject-non-atgc           Reject non-ATGC input sequences after normalization
--no-reject-non-atgc        Keep non-ATGC input sequences in the prepared FASTA

--sina-bin                  SINA executable, default sina
--sina-reference            Optional SINA reference/PTDB path
--allow-sina-failure        Permit fallback behavior if SINA fails
--no-allow-sina-failure     Fail immediately if SINA fails
--fallback-copy-original    Copy prepared FASTA to sina.oriented.fa if SINA output is missing
--no-fallback-copy-original Do not create oriented fallback output
--search-candidates         Ask SINA to write nearest_slv candidates, default on
--no-search-candidates      Skip SINA candidate search
--search-db                 Optional SINA search database path
--search-min-sim            Loose SINA candidate similarity, default 0.5
--search-max-result         Maximum SINA candidates per query, default 10
--search-kmer-candidates    SINA k-mer candidate pool size, default 1000

--vsearch-bin               VSEARCH executable, default vsearch
--strict-tool-version       Fail if SINA/VSEARCH versions cannot be detected
--iddef                     VSEARCH identity definition, default 2
--species-id                Species identity threshold, default 0.972
--genus-id                  Genus identity threshold, default 0.901
--family-id                 Family identity threshold, default 0.801
--order-id                  Order identity threshold, default 0.729
--class-id                  Class identity threshold, default 0.722
--phylum-id                 Phylum identity threshold and search floor, default 0.696
--min-query-cov             Minimum VSEARCH query coverage, default 0.80
--min-target-cov            Minimum VSEARCH target coverage, default 0.0
--maxaccepts                VSEARCH --maxaccepts, default 50
--maxrejects                VSEARCH --maxrejects, default 256
--near-best-delta           Near-best identity window, default 0.005
--rank-consensus            Stable rank consensus fraction, default 0.80
--strand                    VSEARCH search strand: plus or both, default plus
--require-sina-candidates   Fail if SINA candidates do not match representatives
--no-require-sina-candidates
                            Fall back to full representative search, default
--allow-ambiguous           Write ambiguous records instead of failing, default
--no-allow-ambiguous        Fail when ambiguous placements are encountered

--export                    Run export all after placement
--no-export                 Do not export after placement, default
--gzip                      Gzip FASTA exports when --export is used, default
--no-gzip                   Write uncompressed FASTA exports
--force-export              Overwrite existing export files
--validate                  Run validation after summarize/export, default
--no-validate               Skip post-add validation
--strict-validate           Treat selected validation warnings as failures
                            Omit this flag to keep validation warnings non-fatal
```

By default, `add` asks SINA to write loose `nearest_slv` candidates with
`--search-min-sim 0.5` and `--search-max-result 10`, then uses VSEARCH
`--iddef 2` for final registry rescoring. Disable candidate search with
`--no-search-candidates`.

`add` writes `logs/add_start_date*.log` and `logs/add_date*.log`. It does not
export classifier references by default; add `--export` to run `export all`
after placement. Use `--strict-validate` when the dataset should fail on
selected validation warnings.

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
logs/init_start_date20260510152800.log
logs/init_date20260510152944.log
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

## 11. Algorithm Reference and Formulas

This section documents the internal decision logic used by autotax2. It is meant
to make the output auditable: a user should be able to read the evidence tables
and reproduce why a record was assigned, marked ambiguous, or used to create a
new placeholder.

### 11.1 Seven-Rank Model

autotax2 uses exactly seven ranks:

```text
domain -> phylum -> class -> order -> family -> genus -> species
```

The rank list is fixed for every registry and export. Missing lower ranks are
not exported as empty fields; they are either resolved to an existing taxon,
resolved to a placeholder, marked ambiguous/unplaced in evidence, or rejected
before entering the registry.

### 11.2 Default Identity Thresholds

Rank identity thresholds are fixed by default:

```text
theta_species = 0.972
theta_genus   = 0.901
theta_family  = 0.801
theta_order   = 0.729
theta_class   = 0.722
theta_phylum  = 0.696
```

They are used in both SILVA unresolved resolution and downstream dataset
placement. Command-line overrides exist for sensitivity testing, but changing
thresholds changes the scientific meaning of the registry and should be
documented in the final report.

### 11.3 VSEARCH Identity

autotax2 passes `--iddef 2` to VSEARCH by default. In the VSEARCH manual,
`--iddef 2` is the edit-distance identity definition that excludes terminal
gaps. See the official VSEARCH manual for the external definition:
`https://github.com/torognes/vsearch/blob/master/man/vsearch.1`. Conceptually:

```text
I = identity reported by VSEARCH using --iddef 2
```

autotax2 stores `I` as a fraction:

```text
97.2% -> 0.972
90.1% -> 0.901
```

The same identity definition is used for:

```text
1. SILVA parent-scoped unresolved clustering
2. dataset representative clustering
3. dataset-vs-registry rescoring
4. final placement threshold comparisons
```

The important rule is consistency: the threshold table is meaningful only when
all identity values use the same VSEARCH identity definition.

### 11.4 Sequence Normalization and Duplicate Accounting

For each input sequence:

```text
sequence_norm = uppercase(sequence)
sequence_norm = sequence_norm with U converted to T
md5 = MD5(sequence_norm)
```

During `prepare`, input IDs are remapped deterministically:

```text
internal_id_i = <PREFIX>_<six-digit index>

Example:
  prefix = D20
  first sequence  -> D20_000001
  second sequence -> D20_000002
```

If the same MD5 already exists in the registry:

```text
final_status = duplicate
new_taxa_created = 0
export_duplicate_sequence = false
```

The duplicate is still recorded in membership and evidence tables so the
dataset history is not lost.

### 11.5 Legal-Name and Placeholder Catalog

During `init`, autotax2 builds a legal-name catalog from:

```text
1. SILVA type-material lineages
2. GTDB r232 taxonomy names/placeholders from the second TSV column
```

For GTDB:

```text
d__Archaea;p__Methanobacteriota;c__Methanobacteria;...
```

becomes:

```text
Archaea
Methanobacteriota
Methanobacteria
...
```

Names that look like autotax2-generated placeholders are skipped:

```text
SILVAg000001
SILVAs000001
D20g000001
D20s000001
```

Reason: GTDB can contribute accepted names and external placeholder-like labels,
but it must not collide with autotax2's own durable placeholder namespace.

### 11.6 SILVA Named vs Unresolved Classification

For each SILVA record and rank `r`, autotax2 asks:

```text
is rank value present?
is rank value not a generic environmental descriptor?
is rank value legal for this rank?
```

If all are true for all seven ranks:

```text
record_status = named_silva
protected = true
```

If domain is valid but some lower rank fails:

```text
record_status = silva_unresolved
lowest_reliable_rank = highest rank before first failure
unresolved_ranks = first_failed_rank and every lower rank
```

If domain is missing or unresolved:

```text
record_status = rejected
reject_reason = unresolved_domain
```

If a non-type sequence contains symbols outside `A/T/G/C` after normalization:

```text
record_status = rejected
reject_reason = invalid_sequence_characters
```

Type-material sequences are retained even with ambiguous sequence symbols
because they carry naming evidence.

### 11.7 Organelle Exception

Organelle rRNA records are not forced through the strict cellular legal-name
filter. If a SILVA lineage contains:

```text
Chloroplast
Mitochondria
Mitochondrion
Plastid
Apicoplast
Cyanelle
```

autotax2 keeps the upstream lineage and fills the organelle rank and every lower
rank with the organelle label:

```text
Eukaryota;Viridiplantae;Chloroplast
-> Eukaryota;Viridiplantae;Chloroplast;Chloroplast;Chloroplast;Chloroplast;Chloroplast
```

These records enter the named backbone rather than the unresolved placeholder
framework.

### 11.8 SILVA Resolve: Parent-Scoped Hierarchical Clustering

SILVA unresolved records are resolved from high rank to low rank:

```text
phylum -> class -> order -> family -> genus -> species
```

For rank `r`, define:

```text
parent(r) = all ranks above r
P = resolved parent context
C(r, P) = unresolved SILVA records needing rank r under parent P
A(r, P) = known same-parent anchors for rank r under parent P
J(r, P) = A(r, P) union C(r, P)
```

autotax2 runs one parent-local VSEARCH job:

```text
vsearch --cluster_fast J(r, P).fa --id theta_r --iddef 2 --uc J(r, P).uc
```

Decision rule for each unresolved candidate in a parent-local cluster:

```text
if cluster has exactly one distinct same-parent anchor taxon:
    assign existing anchor taxon

elif cluster has two or more distinct same-parent anchor taxa:
    mark this rank ambiguous
    block lower-rank placeholder creation for this record

else:
    create or reuse a SILVA placeholder for this rank and parent-local cluster
```

Important consequence:

```text
Known anchors outside parent P do not compete.
```

This protects autotax2 from assigning a record across a manually revised or
genome-defined higher-rank boundary.

### 11.9 SILVA Placeholder IDs

SILVA placeholders are allocated from durable counters:

```text
phylum  placeholder example: p__SILVAp000001
class   placeholder example: c__SILVAc000001
order   placeholder example: o__SILVAo000001
family  placeholder example: f__SILVAf000001
genus   placeholder example: g__SILVAg000001
species placeholder example: s__SILVAs000001
```

The cluster mapping is durable:

```text
cluster_key -> taxon_id
```

If a run is repeated and the same active cluster key already exists:

```text
reuse existing placeholder
do not allocate a new ID
```

Deprecated placeholder IDs are never reused.

### 11.10 SILVA Resolve Thread Allocation

Within one rank, parent groups are independent after the parent context is
fixed. autotax2 assigns threads to parent-local VSEARCH jobs according to job
size.

For job `j`:

```text
size_j = number of unresolved candidate records in parent group j
weight_j = max(1, size_j)
T = command-level --threads
```

Thread allocation:

```text
job_threads_j = max(1, min(T, round(T * weight_j / sum(weight_all_jobs))))
```

Jobs are then batched so the sum of `job_threads_j` in a batch does not exceed
`T`. This lets small parent groups finish in parallel while larger groups
receive more VSEARCH threads.

### 11.11 SINA Candidate Search

SINA is used as a loose candidate finder and orientation helper. It is not the
final taxonomic identity scorer.

With:

```text
--search-candidates
--search-min-sim 0.5
--search-max-result 10
```

autotax2 asks SINA for up to 10 loose candidates per query:

```text
C_sina(q) = top SINA candidates where SINA similarity >= 0.5
```

Then VSEARCH rescoring is run against the matching current registry
representatives:

```text
C_vsearch(q) = current registry representatives matched from C_sina(q)
```

If no SINA candidates match current representatives:

```text
default:
  fall back to full current representative search

with --require-sina-candidates:
  fail instead of falling back
```

The diagnostic file records which path was used:

```text
datasets/<dataset>/sina_candidate_diagnostics.tsv
```

### 11.12 Dataset VSEARCH Clustering and Registry Search

During `cluster`, autotax2 writes internal VSEARCH clusters and registry search
hits. The most important outputs are:

```text
internal_clusters/species_0.972.uc
vs_registry.filtered.tsv
```

The registry search keeps hits that pass coverage and search settings:

```text
query_coverage >= --min-query-cov
target_coverage >= --min-target-cov
strand follows --strand
```

The default coverage settings are:

```text
--min-query-cov 0.80
--min-target-cov 0.0
--strand plus
```

### 11.13 Near-Best Hit Set

For one query `q`, let VSEARCH hits be:

```text
H(q) = all filtered registry hits for q
I(h) = VSEARCH identity of hit h
I_best(q) = max I(h) over H(q)
delta = --near-best-delta
```

The near-best set is:

```text
NB(q) = { h in H(q) where I(h) >= I_best(q) - delta }
```

Default:

```text
delta = 0.005
```

So if the best hit is `0.980`, then hits down to `0.975` remain in the
near-best evidence set.

### 11.14 Rank-Specific Near-Best Consensus

For each rank `r`, autotax2 applies the rank threshold before consensus:

```text
NB_r(q) = { h in NB(q) where I(h) >= theta_r }
```

Then it collects distinct taxon IDs at rank `r`:

```text
T_r(q) = distinct rank-r taxa represented in NB_r(q)
```

Duplicate hits to the same rank-r taxon do not create extra votes. The current
vote model is taxon-presence based:

```text
support_r(t) = 1 if taxon t is present in T_r(q), else 0
fraction_r(t) = support_r(t) / sum support_r(all taxa in T_r(q))
```

A rank is stable when:

```text
max fraction_r(t) >= --rank-consensus
```

Default:

```text
--rank-consensus 0.80
```

With the default, a rank with exactly one passing taxon is stable. A rank with
two or more conflicting taxon IDs is normally unstable because each distinct
taxon contributes one presence vote.

### 11.15 Identity Status Formula

The best VSEARCH identity controls the initial status:

```text
I_best >= theta_species:
    identity_status = known_like

theta_genus <= I_best < theta_species:
    identity_status = new_species

theta_family <= I_best < theta_genus:
    identity_status = new_genus

theta_order <= I_best < theta_family:
    identity_status = new_family

theta_class <= I_best < theta_order:
    identity_status = new_order

theta_phylum <= I_best < theta_class:
    identity_status = new_class

I_best < theta_phylum:
    identity_status = unplaced
```

Identity status alone does not create a placeholder. A stable parent rank is
also required.

### 11.16 Placement Decision Rules

For `identity_status = known_like`:

```text
if species consensus is stable:
    assign existing species

else if the highest-identity species is unique and I_best >= theta_species:
    assign the unique highest-identity species lineage

else:
    final_status = ambiguous
```

This means:

```text
query vs Species_A = 0.980
query vs Species_B = 0.976
both >= 0.972
highest species is unique
=> assign Species_A
```

But:

```text
query vs Species_A = 0.980
query vs Species_B = 0.980
both >= 0.972
highest species tie across different species
=> ambiguous
```

For `identity_status = new_species`:

```text
if genus consensus is stable:
    create species placeholder under that genus
elif family consensus is stable:
    create genus placeholder and species placeholder under that family
else:
    ambiguous
```

For `identity_status = new_genus`:

```text
if family consensus is stable:
    create genus + species placeholders
elif order consensus is stable:
    create family + genus + species placeholders
else:
    ambiguous
```

For `identity_status = new_family`:

```text
if order consensus is stable:
    create family + genus + species placeholders
elif class consensus is stable:
    create order + family + genus + species placeholders
else:
    ambiguous
```

For `identity_status = new_order`:

```text
if class consensus is stable:
    create order + family + genus + species placeholders
elif phylum consensus is stable:
    create class + order + family + genus + species placeholders
else:
    ambiguous
```

For `identity_status = new_class`:

```text
if phylum consensus is stable:
    create class + order + family + genus + species placeholders
else:
    ambiguous
```

For `identity_status = unplaced`:

```text
final_status = unplaced
new_taxa_created = 0
```

### 11.17 Dataset Placeholder IDs

Dataset placeholders use the frozen dataset prefix:

```text
g__D20g000001
s__D20s000001
```

For a created rank:

```text
cluster_key = <PREFIX>|<rank>|<parent_taxon_id>|<threshold>|<query_id>
```

If that active cluster key already exists:

```text
reuse existing taxon_id
```

Otherwise:

```text
allocate next durable placeholder ID
write taxon_nodes.tsv
write cluster_to_taxon.tsv
write created_taxa.tsv
```

### 11.18 Type-Strain and Same-Species Evidence

Type-material sequences are preferred as representatives, but they are not a
special bypass around the thresholds.

Rule:

```text
If query vs type-material representative < theta_species,
but query vs another named SILVA sequence from the same species >= theta_species,
then that species can still be supported.
```

Reason:

```text
The species is represented by the named SILVA species set, not only by one
selected representative.
```

Only when all sequences supporting that species fall below `theta_species` can
the query become eligible for a new species placeholder under the placement
rules.

### 11.19 Evidence Tables

SILVA resolve evidence:

```text
silva/silva_unresolved_evidence.tsv
```

One row per unresolved SILVA sequence per rank. Important columns:

```text
rank
parent_rank
parent_taxon
threshold
candidate_scope
anchor_count
best_anchor_taxon
best_anchor_identity
passing_anchor_count
competing_anchor_taxa
residual_cluster_key
decision
output_taxon
reason
job_size
job_threads
```

Dataset placement evidence:

```text
datasets/<dataset>/placement_evidence.tsv
```

One row per dataset representative per rank. Important decisions:

```text
assign_existing
create_placeholder
duplicate
ambiguous
unplaced
```

Important reasons:

```text
stable_near_best_consensus
highest_identity_species_lineage
exact_sequence_md5_duplicate
no_registry_hit
multiple_same_parent_anchor_clusters
no_same_parent_anchor_above_threshold
```

### 11.20 Overlap and Novelty Reports

`dataset_rank_overlap_detail.tsv` reports whether each rank assignment overlaps
with previous registry state or is new:

```text
overlap_named_silva
overlap_unresolved_silva
overlap_previous_custom
new_current_dataset
ambiguous
unplaced
```

`dataset_rank_novelty_detail.tsv` lists newly created taxa:

```text
dataset
rank
taxon_id
name
parent_taxon_id
cluster_key
representative_seq_id
support_sequence_ids
```

Together with `placement_evidence.tsv`, these reports answer:

```text
Which ranks overlapped existing reference taxa?
Which ranks were new in this dataset?
Which query sequences supported each new placeholder?
Which competing taxa prevented a confident assignment?
```

### 11.21 Export Semantics

Exports are derived from the active registry:

```text
SINTAX reference
QIIME 2 reference sequences
QIIME 2 taxonomy table
DADA2 toGenus trainset
DADA2 assignSpecies file
```

By default:

```text
--representatives-only = true
```

This exports active representatives rather than every active unique MD5. Use
`--all-unique` when a downstream workflow needs all active unique sequences.

Exact duplicate MD5 sequences are not exported repeatedly.

### 11.22 Validation Invariants

`validate` checks:

```text
1. protected named SILVA taxa still match protected_taxa_snapshot.tsv
2. placeholder IDs are unique
3. placeholder counters have not moved backward
4. exact duplicate MD5 accounting is consistent
5. representative rows point to active taxa
6. cluster_to_taxon mappings point to active taxa
7. export files, when present, satisfy their format contracts
8. SILVA resolve evidence covers unresolved SILVA records
9. dataset placement evidence covers assignment rows from phylum to species
10. SINA and VSEARCH metadata are traceable
```

Run validation after every major step:

```bash
autotax2 validate --build "$BUILD" --no-check-exports
autotax2 validate --build "$BUILD" --check-exports
```

Use `--strict` when validation warnings should fail the command.

## 12. Frequently Asked Questions

### 12.1 Do externally processed 16S datasets need to be extracted again?

No. autotax2 assumes that `prepare --fasta` already points to target
SSU/16S sequences.

### 12.2 Does autotax2 run barrnap?

No. barrnap or any other extraction workflow may be used before autotax2, but
autotax2 does not invoke, validate, or record that preprocessing step.

### 12.3 Why are SILVA records with empty domains rejected?

The domain rank is the root of the rank-aware backbone. Records with empty or
unresolved domains cannot be placed safely under Archaea or Bacteria, so they
are written to `silva_rejected.tsv`.

### 12.4 Can placeholder IDs be reused?

No. Placeholder counters are durable, and deprecated IDs are never reused.

### 12.5 What happens when identical sequences appear in multiple datasets?

Each sequence receives an MD5 digest. Exact duplicate sequences retain
membership records, but the same sequence MD5 is not exported repeatedly.

## 13. Developer Checks

```bash
python -m pytest
python -m compileall autotax2
autotax2 --help
autotax2 prepare --help
```
