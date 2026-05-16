# autotax2 Debug Workflow

Date: 2026-05-15

This document is the working debug map for autotax2. It records the current
implementation flow, expected inputs and outputs, core algorithms, parameters,
debug checkpoints, and optimization items. Keep it updated whenever a rule or
file contract changes.

For the first real SINA/VSEARCH smoke run, use `docs/real_run_checklist.md`.

## 0. Global Contract

autotax2 is a fixed-backbone, rank-aware, incremental rRNA reference builder.
The current command path is:

```text
init -> resolve -> prepare -> orient -> cluster -> place -> export
                         \-> summarize
                         \-> validate
```

The `add` command exists in the CLI but is currently only a placeholder. Use the
explicit step-by-step commands for real runs.

Core invariants:

- Named SILVA taxa are protected and immutable.
- SILVA unresolved records can become mutable SILVA placeholder framework.
- Custom datasets get frozen prefixes during `prepare`, such as `D20`.
- Input sequence IDs are remapped to internal IDs, such as `D20_000001`.
- Exact sequence MD5 is tracked to avoid repeated exports.
- VSEARCH identity defaults to `--iddef 2`.
- Input dataset FASTA files are already extracted SSU/16S sequences.
- SINA is used for loose orientation, not internal SSU extraction.
- Placement uses near-best hit consensus, not best hit only.

Default identity thresholds:

| Rank boundary | Default |
|---|---:|
| species | 0.972 |
| genus | 0.901 |
| family | 0.801 |
| order | 0.729 |
| class | 0.722 |
| phylum | 0.696 |

## 1. Build Initialization

Command:

```text
autotax2 init \
  --silva-fasta SILVA.fa.gz \
  --type-strain-metadata SILVA.full_metadata.tsv.gz \
  --gtdb-ar53-taxonomy ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy bac120_taxonomy_r232.tsv \
  --outdir BUILD \
  --threads 8
```

Implementation:

- CLI: `autotax2/cli.py::init`
- Core: `autotax2/silva.py::initialize_silva_build`

Inputs:

- `--silva-fasta`: official SILVA NR99 taxonomy FASTA, plain or gzipped.
- `--type-strain-metadata`: required official SILVA full_metadata TSV/TSV.gz.
  The parser requires accession and type-material-like fields.
- `--gtdb-ar53-taxonomy`: required GTDB r232 archaeal taxonomy TSV/TSV.gz.
- `--gtdb-bac120-taxonomy`: required GTDB r232 bacterial taxonomy TSV/TSV.gz.
  For both GTDB files, autotax2 reads only the second column, splits the GTDB
  taxonomy string on semicolons, strips rank prefixes such as `p__`, `g__`, and
  `s__`, and treats the remaining names/placeholders as accepted legal names.
  Names that look like autotax2/SILVA placeholder IDs, such as
  `SILVAg000001` or `D20s000001`, are skipped to prevent namespace overlap.
- `--outdir`: build directory.
- `--threads`: recorded in build metadata and reserved for initialization
  steps.

Algorithm:

1. Read SILVA full_metadata first and retain only type-strain/type-material
   evidence by sequence accession.
2. Read SILVA FASTA and normalize sequence text.
   - `U` is converted to `T`.
   - Non-type sequences must contain only `A/T/G/C`.
   - Type-strain sequences are retained even if they contain illegal sequence
     symbols such as `N` or other ambiguity characters.
3. Reject invalid records into `silva/silva_rejected.tsv`.
   - Empty or unresolved domain -> rejected.
   - Non-type sequence with illegal characters -> `invalid_sequence_characters`.
4. Build a strict legal-name catalog.
   - SILVA type-material records contribute legal domain/phylum/class/order/
     family/genus/species names.
   - GTDB r232 taxonomy contributes accepted names and candidate placeholders
     such as Asgard/Loki-style names.
   - Species names from type material are canonicalized to the binomial, with
     strain/isolate/clone suffixes removed.
5. Classify records.
   - If a rank contains unresolved tokens such as `uncultured`, `unidentified`,
     `metagenome`, or standalone `Incertae Sedis`, it is unresolved from that
     rank downward.
   - If the legal-name catalog is available and a rank value is not accepted by
     it, it becomes unresolved from that rank downward.
   - Organelle lineages bypass strict cellular legal-name filtering. If a rank
     contains `Chloroplast`, `Mitochondria`/`Mitochondrion`, `Plastid`,
     `Apicoplast`, or `Cyanelle`, that rank and every lower rank are
     standardized to the organelle label, producing a named seven-rank backbone
     record.
   - Species values such as `Methanobacterium formicicum DSM 1535` are stored as
     `Methanobacterium formicicum`.
   - Species tails such as `uncultured archaeon SAGMA-B` are unresolved and the
     strain-like tail is not kept as a name.
6. Split records into named backbone and unresolved SILVA records.
7. Write registries, snapshots, and audit logs.

Outputs:

| Path | Meaning |
|---|---|
| `silva/silva_named_backbone.fa` | Protected named SILVA sequences |
| `silva/silva_named_backbone.tax.tsv` | Named SILVA taxonomy table |
| `silva/silva_unresolved.fa` | SILVA records needing placeholder resolution |
| `silva/silva_unresolved.tsv` | Unresolved rank reasons and reliable parent |
| `silva/silva_rejected.tsv` | Rejected records and reasons |
| `registry/taxon_nodes.tsv` | Protected named SILVA taxon tree |
| `registry/sequence_registry.tsv` | Sequence MD5, type-strain metadata, taxonomy |
| `registry/legal_name_catalog.tsv` | Parsed SILVA type-material + GTDB r232 legal-name catalog |
| `registry/name_index.tsv` | Name to taxon lookup |
| `registry/representative_registry.tsv` | Active species representatives |
| `registry/protected_taxa_snapshot.tsv` | SILVA immutability baseline |
| `registry/placeholder_counters.yaml` | Placeholder allocation counters |
| `logs/init_start_date*.log`, `logs/init_date*.log` | Audit logs |

Organelle standardization example:

```text
Input:  Eukaryota;Viridiplantae;Chloroplast
Output: Eukaryota;Viridiplantae;Chloroplast;Chloroplast;Chloroplast;Chloroplast;Chloroplast

Input:  Eukaryota;Opisthokonta;Metazoa;Mitochondrion;uncultured;metagenome
Output: Eukaryota;Opisthokonta;Metazoa;Mitochondria;Mitochondria;Mitochondria;Mitochondria
```

Important parameters:

- `--type-strain-metadata`: required official SILVA full_metadata.
- `--gtdb-ar53-taxonomy` and `--gtdb-bac120-taxonomy`: required once during
  initialization. Downstream commands use the initialized build directory.
- `--threads`: supported and recorded; currently mostly an initialization
  metadata/audit parameter.

Debug checklist:

- Check `silva_rejected.tsv` before any downstream step.
- Check `silva_unresolved.tsv` for `lowest_reliable_rank` and `unresolved_reason`.
- Check `sequence_registry.tsv` for `is_type_strain`, `type_species_name`, and
  `type_strain_id`.
- Check `registry/legal_name_catalog.tsv` to confirm GTDB names/placeholders were
  parsed after prefix removal and reserved autotax2 placeholder-like names were
  skipped.
- Check organelle records in `silva_named_backbone.tax.tsv`: downstream ranks
  from the organelle marker should be filled by the canonical organelle label,
  not by `unresolved_*`, `uncultured`, or environmental tails.
- Check that `protected_taxa_snapshot.tsv` is created and later remains stable.

Optimization items:

- Add explicit summary counts for invalid non-type SILVA sequences and retained
  invalid type-strain sequences.
- Add optional LPSN valid-name input if strict type-material + GTDB coverage is
  still too narrow for a target SILVA release.

## 2. Resolve SILVA Unresolved Records

Command:

```text
autotax2 resolve --build BUILD --threads 8 --vsearch-bin vsearch
```

Implementation:

- CLI: `autotax2/cli.py::resolve_silva`
- Core: `autotax2/silva.py::resolve_silva_unresolved_build`

Inputs:

- `silva/silva_unresolved.fa`
- `silva/silva_unresolved.tsv`
- `registry/placeholder_counters.yaml`
- Existing `registry/cluster_to_taxon.tsv`, if any.

Algorithm:

1. Read unresolved SILVA records.
2. Process unresolved records rank by rank with parent-scoped VSEARCH jobs:
   - phylum within domain at `0.696`
   - class within phylum at `0.722`
   - order within class at `0.729`
   - family within order at `0.801`
   - genus within family at `0.901`
   - species within genus at `0.972`
3. Reuse existing parent-local UC files if present. If a UC file is absent, run
   VSEARCH on a parent-local FASTA containing same-parent anchors followed by
   unresolved candidates. If VSEARCH fails, write singleton UC assignments so
   the resolver can still produce deterministic placeholder proposals. Delete
   fallback UC files before rerunning if real VSEARCH clustering should be
   regenerated.
4. Schedule same-rank parent-local jobs in batches. Per-job `job_threads` is
   assigned from candidate count, and the sum of threads in one active batch
   does not exceed command-level `--threads`. Placeholder allocation is applied
   after job completion in deterministic parent order.
5. For each unresolved rank:
   - keep reliable upstream ranks as prefixed labels, such as `p__Aenigmarchaeota`;
   - allocate SILVA placeholders for unresolved ranks, such as
     `g__SILVAg000001` or `s__SILVAs000001`;
   - species placeholder keys include genus context to avoid merging species
     placeholders across different genus frameworks.
6. Assign unresolved records to a unique same-parent known anchor when one
   exists; mark parent-local ambiguity when multiple anchors pass the threshold.
7. Reuse active cluster-to-taxon mappings where possible.
8. Write resolved taxonomy, per-rank evidence, and registry updates unless
   `--dry-run` is used.

Strict resolver behavior:

- Uses parent-scoped hierarchical resolution rather than global all-unresolved
  clustering.
- Process ranks from high to low:
  - phylum within domain at `0.696`
  - class within phylum at `0.722`
  - order within class at `0.729`
  - family within order at `0.801`
  - genus within family at `0.901`
  - species within genus at `0.972`
- At each rank, keep the immediate parent fixed. Same-rank anchors outside the
  fixed parent do not compete, even when marker-gene identity is high, because
  upper ranks may reflect genome, morphology, physiology, biochemical evidence,
  or later taxonomic revision.
- Use same-parent known children as anchors. They do not create new
  placeholders, but unresolved children must first be tested against them.
- If exactly one same-parent anchor is `>=` the rank threshold, assign the
  unresolved record to that existing child.
- If more than one same-parent anchor is `>=` the rank threshold, mark the rank
  ambiguous inside that parent. Do not pick one arbitrarily and do not create a
  novel placeholder to hide ambiguity.
- If all same-parent anchors are `<` the rank threshold, send the residual
  sequence into parent-local novel clustering and create/reuse placeholders from
  those residual clusters.

Inference evidence contract:

- Every unresolved SILVA record must get machine-readable per-rank evidence.
- The same evidence model must also be used for later custom dataset placement
  so sequence-level debugging follows one logic across SILVA resolve and
  incremental dataset addition.
- Recommended fields:
  - `source_stage`
  - `dataset`
  - `seq_id`
  - `rank`
  - `parent_rank`
  - `parent_taxon`
  - `input_rank_value`
  - `input_rank_status`
  - `threshold`
  - `candidate_scope`
  - `anchor_count`
  - `best_anchor_taxon`
  - `best_anchor_identity`
  - `passing_anchor_count`
  - `competing_anchor_taxa`
  - `residual_cluster_key`
  - `residual_cluster_size`
  - `decision`
  - `output_taxon`
  - `reason`
  - `job_size`
  - `job_threads`
- SILVA output: `silva/silva_unresolved_evidence.tsv`.
- Dataset output: `datasets/NN_name/placement_evidence.tsv`.

Parent-scoped parallel scheduling:

- Treat each same-parent, same-rank anchor-search/residual-clustering task as
  an independent job after its parent is fixed.
- The command-level `--threads` value is a global cap, not a value passed to
  every VSEARCH call.
- Allocate per-job threads dynamically by task size, such as residual sequence
  count or total bases.
- Small jobs should usually get one thread; large jobs can get more threads.
- Current implementation processes ranks with barriers from high to low, and
  batches same-rank jobs so the sum of active VSEARCH threads does not exceed
  `--threads`.
- A later DAG scheduler can unlock child jobs as soon as their parent decisions
  are complete.

Outputs:

| Path | Meaning |
|---|---|
| `silva/silva_unresolved_clusters/<rank>/*_<parent-hash>.uc` | Parent-local rank cluster assignments |
| `silva/silva_unresolved_taxa.tsv` | Created/reused SILVA placeholder taxa |
| `silva/silva_unresolved_members.tsv` | Per-sequence placeholder membership |
| `silva/silva_unresolved_mapping.tsv` | Original SILVA taxonomy to autotax2 taxonomy |
| `silva/silva_unresolved_evidence.tsv` | Per-sequence, per-rank reasoning evidence |
| `silva/silva_unresolved.resolved.fa` | Resolved SILVA unresolved FASTA |
| `registry/taxon_nodes.tsv` | Updated with SILVA placeholder taxa |
| `registry/cluster_to_taxon.tsv` | Stable cluster-to-placeholder mapping |
| `registry/placeholder_counters.*` | Advanced SILVA counters |

Important parameters:

- `--species-id`: default `0.972`.
- `--genus-id`: default `0.901`.
- `--family-id`: default `0.801`.
- `--order-id`: default `0.729`.
- `--class-id`: default `0.722`.
- `--phylum-id`: default `0.696`.
- `--iddef`: default `2`.
- `--dry-run`: write proposal files without registry mutation.

Debug checklist:

- Inspect `silva_unresolved_members.tsv` for placeholder names and warnings.
- Inspect `silva_unresolved_evidence.tsv` for parent scope, same-parent anchor
  decisions, ambiguity, residual cluster keys, and job thread allocation.
- Run `validate` after registry mutation.

Optimization items:

- Replace rank barriers with a DAG scheduler if parent decisions can safely
  unlock child jobs earlier.
- Add optional SINA candidate finding before parent-local VSEARCH scoring.

## 3. Prepare A Custom Dataset

Command:

```text
autotax2 prepare --build BUILD --name D20 --prefix D20 --fasta input.fa --domain Archaea
```

Implementation:

- CLI: `autotax2/cli.py::prepare_dataset_command`
- Core: `autotax2/dataset_prepare.py::prepare_dataset`

Inputs:

- Externally extracted SSU/16S FASTA.
- Build registry from `init` and optionally `resolve`.
- Dataset name, frozen prefix, and domain.

Algorithm:

1. Validate dataset prefix.
   - Prefix must start with a letter and contain only letters and digits.
   - `SILVA` is reserved for the backbone.
2. Register dataset order and output directory:
   - `datasets/01_<safe_dataset_name>/`
3. Normalize input sequences.
   - `U` -> `T`.
   - sequence MD5 is computed after normalization.
4. Assign internal IDs:
   - `D20_000001`, `D20_000002`, ...
5. Reject sequences:
   - non-ATGC if `--reject-non-atgc` is enabled;
   - shorter than the domain-specific minimum length.
6. Build per-dataset unique sequence registry and duplicate membership table.

Outputs:

| Path | Meaning |
|---|---|
| `datasets/NN_name/input.normalized.fa` | Normalized input except non-ATGC rejects |
| `datasets/NN_name/prepared.ssu.fa` | Final accepted SSU/16S sequences |
| `datasets/NN_name/sequence_id_map.tsv` | Original ID to internal ID mapping |
| `datasets/NN_name/unique_sequence_registry.tsv` | Dataset-local unique MD5 table |
| `datasets/NN_name/sequence_membership.tsv` | Internal sequence to unique MD5 membership |
| `datasets/NN_name/prepare_summary.tsv` | Counts and rejection summary |
| `registry/dataset_registry.tsv` | Frozen dataset prefix and add order |

Important parameters:

- `--prefix`: frozen dataset prefix.
- `--domain`: `Archaea` or `Bacteria`.
- `--min-ssu-len-archaea`: default `900`.
- `--min-ssu-len-bacteria`: default `1200`.
- `--reject-non-atgc/--no-reject-non-atgc`: default reject.

Debug checklist:

- Check `prepare_summary.tsv`.
- Check `sequence_id_map.tsv` for `rejected=true` rows.
- Check `sequence_membership.tsv` for duplicate MD5 behavior.
- Ensure the input FASTA was already SSU/16S extracted. autotax2 does not do
  extraction internally.

Optimization items:

- Add an option to record a rejected FASTA for rejected custom sequences.
- Add cross-dataset MD5 overlap during prepare, not only in reports.

## 4. Orient With SINA

Command:

```text
autotax2 orient --build BUILD --dataset D20 --threads 8 --sina-bin sina
```

Implementation:

- CLI: `autotax2/cli.py::orient_sina_command`
- Core: `autotax2/sina.py::orient_dataset_with_sina`

Inputs:

- `datasets/NN_name/prepared.ssu.fa`
- Optional SINA reference path via `--reference`.
- Optional SINA search database via `--search-db` when candidate search is used.

Algorithm:

1. Build a SINA command using loose orientation settings.
2. Run SINA and capture stderr into `sina.log`.
3. Compare each SINA output sequence with the input sequence.
   - Same sequence -> plus orientation.
   - Reverse complement -> minus orientation.
   - Other modifications -> retained but flagged as `sina_modified_sequence`.
4. If SINA fails or omits records:
   - default behavior is conservative fallback to original sequences.
   - fallback behavior is recorded in `sina.summary.tsv`.
5. If `--search-candidates` is enabled:
   - ask SINA to run loose `nearest_slv` candidate search;
   - parse the candidate CSV into `sina.candidates.tsv`;
   - do not use SINA similarity as an autotax2 threshold identity.
     VSEARCH `--iddef 2` remains the final threshold scorer.

Outputs:

| Path | Meaning |
|---|---|
| `datasets/NN_name/sina.oriented.fa` | Oriented sequences for VSEARCH |
| `datasets/NN_name/sina.summary.tsv` | Per-sequence orientation status |
| `datasets/NN_name/sina.candidates.tsv` | Optional SINA candidate target list for VSEARCH rescoring |
| `datasets/NN_name/sina.log` | SINA stderr/log text |
| `datasets/NN_name/tool_versions.tsv` | SINA version/status metadata |

Important parameters:

- `--allow-sina-failure/--no-allow-sina-failure`: default allow fallback.
- `--fallback-copy-original/--no-fallback-copy-original`: default copy original.
- `--strict-tool-version`: fail if SINA version cannot be parsed.
- `--min-sina-identity`, `--min-sina-score`: currently unsupported; leave `0.0`.
- `--search-candidates`: optional candidate discovery; default disabled.
- `--search-min-sim`: loose SINA search cutoff; default `0.500`.
- `--search-max-result`: maximum SINA candidates per query; default `10`.

Debug checklist:

- Check `sina.summary.tsv` for `fallback_used=true`.
- Check `sina.candidates.tsv` if candidate search was enabled.
- Check `warning=sina_modified_sequence`.
- Check `tool_versions.tsv` for SINA command and version.

Optimization items:

- Parse SINA alignment/search metrics if they will be used for filtering.

## 5. Cluster And Search Current Registry

Command:

```text
autotax2 cluster --build BUILD --dataset D20 --threads 8 --vsearch-bin vsearch --iddef 2
```

Implementation:

- CLI: `autotax2/cli.py::cluster_search_command`
- Core: `autotax2/vsearch.py::cluster_search_dataset`

Inputs:

- `datasets/NN_name/sina.oriented.fa`
- Optional `datasets/NN_name/sina.candidates.tsv`
- Current registry:
  - `registry/representative_registry.tsv`
  - `registry/sequence_registry.tsv`
  - `registry/taxon_nodes.tsv`
  - resolved SILVA FASTA files and dataset FASTA files for sequence lookup.

Algorithm:

1. Cluster oriented dataset sequences at all rank thresholds:
   - species, genus, family, order, class, phylum.
2. Write internal cluster membership tables.
3. Rebuild `registry/current_representatives.fa`.
   - Active representatives are included.
   - Named SILVA species evidence sequences are also included so non-type named
     same-species sequences can support a species call.
4. If `sina.candidates.tsv` exists, build
   `registry/current_representatives.sina_candidates.fa` from matching candidate
   targets.
   - If at least one candidate target matches, VSEARCH searches the subset.
   - If no candidate target matches, default behavior falls back to the full
     registry.
   - `--require-sina-candidates` makes unmatched candidates fatal.
5. Search species centroids against current or candidate-subset representatives
   using VSEARCH `--iddef 2`.
6. Filter registry hits by:
   - identity >= phylum floor;
   - query coverage >= `--min-query-cov`;
   - target coverage >= `--min-target-cov`.
7. Preserve all passing candidates for later near-best consensus.

Outputs:

| Path | Meaning |
|---|---|
| `datasets/NN_name/internal_clusters/*_*.uc` | Rank-specific VSEARCH UC files |
| `datasets/NN_name/internal_clusters/*_*.centroids.fa` | Rank centroids |
| `datasets/NN_name/internal_clusters/*_*.members.tsv` | Cluster membership tables |
| `registry/current_representatives.fa` | Current search database |
| `registry/current_representatives.sina_candidates.fa` | Candidate subset database when SINA candidates are used |
| `registry/current_representatives.tax.tsv` | Current representative taxonomy |
| `datasets/NN_name/vs_registry.raw.tsv` | Raw VSEARCH userout |
| `datasets/NN_name/vs_registry.tsv` | Parsed raw hits |
| `datasets/NN_name/vs_registry.filtered.tsv` | Hits used by placement |
| `datasets/NN_name/cluster_search_summary.tsv` | Cluster/search counts |
| `datasets/NN_name/sina_candidate_diagnostics.tsv` | Per-query SINA candidate matching diagnostics |
| `datasets/NN_name/tool_versions.tsv` | VSEARCH version/status metadata |

Important parameters:

- `--iddef`: default `2`.
- `--species-id`: default `0.972`.
- `--genus-id`: default `0.901`.
- `--family-id`: default `0.801`.
- `--order-id`: default `0.729`.
- `--class-id`: default `0.722`.
- `--phylum-id`: default `0.696`; also the current search floor.
- `--min-query-cov`: default `0.80`.
- `--min-target-cov`: default `0.0`.
- `--maxaccepts`: default `50`.
- `--maxrejects`: default `256`.
- `--near-best-delta`: default `0.005`, recorded for placement consistency.
- `--strand`: `plus` or `both`, default `plus`.
- `--sina-candidates`: explicit SINA candidate TSV; by default cluster auto-uses
  `datasets/NN_name/sina.candidates.tsv` when present.
- `--require-sina-candidates`: fail if no SINA targets match current
  representatives.

Debug checklist:

- Check `cluster_search_summary.tsv` for centroid counts and filtered hit count.
- Check `sina_candidate_*` columns in `cluster_search_summary.tsv` when SINA
  candidates are used.
- Check `sina_candidate_diagnostics.tsv` to see whether each input query had
  SINA targets, whether those targets matched current representatives, and
  whether VSEARCH searched a candidate subset or fell back to the full registry.
- Confirm `current_representatives.fa` is not stale.
- If many records are unplaced later, inspect `vs_registry.raw.tsv` before the
  filtered file to distinguish search failure from coverage/identity filtering.

Optimization items:

- Consider storing best-hit group diagnostics directly in the cluster summary.

## 6. Place Dataset Representatives

Command:

```text
autotax2 place --build BUILD --dataset D20
```

Implementation:

- CLI: `autotax2/cli.py::place_command`
- Core: `autotax2/placement.py::place_dataset`

Inputs:

- `datasets/NN_name/vs_registry.filtered.tsv`
- `datasets/NN_name/sequence_membership.tsv`
- `datasets/NN_name/sequence_id_map.tsv`
- Current registry taxon, sequence, cluster, and representative tables.

Algorithm:

1. Load filtered VSEARCH hits.
2. Determine query representatives.
3. First handle exact MD5 duplicates:
   - exact duplicates of existing sequences get `duplicate` status.
4. For each non-duplicate query:
   - sort hits by identity;
   - retain near-best hits with identity >= best identity - `near_best_delta`;
   - compute rank consensus from near-best hits after applying each rank's own
     identity threshold.
5. Consensus counts distinct taxon IDs per rank rather than raw hit count.
   - This prevents one species with many strains from dominating consensus only
     by redundancy.
6. Determine identity status:
   - `identity >= species` -> `known_like`
   - `genus <= identity < species` -> `new_species`
   - `family <= identity < genus` -> `new_genus`
   - `order <= identity < family` -> `new_family`
   - `class <= identity < order` -> `new_order`
   - `phylum <= identity < class` -> `new_class`
   - `< phylum` -> `unplaced`
7. Decide final placement from identity status plus stable rank consensus.
   - If two or more species-level candidates are `>= 0.972` and species
     consensus is not stable, autotax2 inherits the highest-identity hit's
     species lineage. It must not create a new species placeholder for this
     case.
   - If different species tie for the highest identity, there is no unique
     highest-identity lineage, so the query remains `ambiguous`.
   - Near-best hits below the species threshold can still inform higher ranks,
     but they do not compete in the species decision.
8. Create new placeholder lineage where needed.
9. Write `placement_evidence.tsv` with the same per-rank evidence columns used
   by SILVA resolve. Each query representative gets phylum-to-species decision
   rows with threshold, near-best scope, best-hit rank taxon, identity,
   consensus support count, final rank decision, output taxon, and placeholder
   cluster key where applicable.
10. Update registry unless `--dry-run` is used.

Outputs:

| Path | Meaning |
|---|---|
| `datasets/NN_name/assignments.tsv` | Final per-query placement decisions |
| `datasets/NN_name/created_taxa.tsv` | Newly created placeholder taxa |
| `datasets/NN_name/near_best_consensus.tsv` | Rank consensus evidence |
| `datasets/NN_name/placement_evidence.tsv` | Per-query, per-rank placement reasoning |
| `datasets/NN_name/representative_updates.tsv` | Representatives to add |
| `datasets/NN_name/placement_summary.tsv` | Counts by status |
| `registry/taxon_nodes.tsv` | Updated active taxa |
| `registry/sequence_registry.tsv` | Updated sequence registry |
| `registry/cluster_to_taxon.tsv` | Updated reusable placement mappings |
| `registry/representative_registry.tsv` | Updated active representatives |
| `registry/placeholder_counters.*` | Advanced custom dataset counters |

Important parameters:

- `--near-best-delta`: default `0.005`.
- `--rank-consensus`: default `0.80`.
- rank identity thresholds as above.
- `--dry-run`: write `*.dry_run.tsv` files only.
- `--allow-ambiguous/--no-allow-ambiguous`: default allow.

Debug checklist:

- Start with `placement_summary.tsv`.
- For any unexpected call, inspect the row in `assignments.tsv`.
- Check `near_best_consensus.tsv` for rank fraction and hit count.
- Check `placement_evidence.tsv` for threshold, best-hit identity, consensus
  support, final rank decision, and placeholder cluster key.
- Check `created_taxa.tsv` for placeholder parentage.
- If `known_like` was assigned despite unstable species consensus, check whether
  the warning is `species_consensus_unstable_best_identity_lineage`.

Optimization items:

- Add a specific report for type-strain versus non-type named same-species
  evidence.
- Add optional `--hit-finder sina-vsearch`.

## 7. Export References

Command:

```text
autotax2 export all --build BUILD --force
```

Implementation:

- CLI: `autotax2/cli.py::export_command`
- Core: `autotax2/export/__init__.py::export_references`

Inputs:

- Current registry.
- Sequence FASTA sources from SILVA and datasets.

Algorithm:

1. Build exportable records from active representatives by default.
2. De-duplicate records by sequence MD5.
3. Resolve each record to exactly 7 taxonomy ranks.
4. Export requested formats.
5. Run export self-validation.

Outputs:

| Path | Meaning |
|---|---|
| `export/sintax/<prefix>.sintax.fa(.gz)` | SINTAX reference FASTA |
| `export/qiime2/reference_sequences.fasta(.gz)` | QIIME2 sequences |
| `export/qiime2/reference_taxonomy.tsv` | QIIME2 taxonomy |
| `export/dada2/*` | DADA2 trainset and species files |
| `export/export_manifest.tsv` | Export manifest |
| `export/export_validation.tsv` | Export self-check report |

Important parameters:

- `format_name`: `all`, `sintax`, `qiime2`, or `dada2`.
- `--gzip/--no-gzip`: default gzip FASTA outputs.
- `--representatives-only/--all-unique`: default representatives only.
- `--prefix`: default `autotax2`.
- `--force`: overwrite existing export files.

Debug checklist:

- If export fails, read `export_validation.tsv`.
- Confirm `records_exported` matches expectation after MD5 de-duplication.
- Confirm all exported taxonomy strings have exactly seven ranks.

Optimization items:

- Add format-specific smoke commands or import recipes for QIIME2/DADA2.

## 8. Summarize Build

Command:

```text
autotax2 summarize --build BUILD --overwrite
```

Implementation:

- CLI: `autotax2/cli.py::summarize`
- Core: `autotax2/reports.py::summarize_build`

Outputs:

| Path | Meaning |
|---|---|
| `reports/global_summary.tsv` | Global build counts |
| `reports/dataset_delta_summary.tsv` | Per-dataset added/known/ambiguous counts |
| `reports/dataset_overlap_matrix.tsv` | Exact and rank overlap against sources |
| `reports/rank_novelty_summary.tsv` | New versus existing taxa by rank |
| `reports/dataset_rank_overlap_detail.tsv` | Rank/taxon overlap detail with supporting sequence IDs |
| `reports/dataset_rank_novelty_detail.tsv` | Created taxon detail with parent, representative, and support |
| `reports/source_contribution.tsv` | Source contribution by sequence/taxon |
| `reports/representative_summary.tsv` | Active representative details |
| `reports/sequence_dedup_summary.tsv` | MD5 de-duplication summary |
| `reports/dataset_increment_audit.md` | Human-readable per-dataset increment audit |

Debug checklist:

- Use `global_summary.tsv` to confirm SILVA resolve evidence coverage via
  `silva_unresolved_evidence_records` and `silva_unresolved_evidence_rows`.
- Use `dataset_delta_summary.tsv` to understand what each dataset changed.
  It includes SINA candidate counts and placement evidence row counts.
- Use `dataset_overlap_matrix.tsv` to answer overlap/new rank questions.
- Use `dataset_rank_overlap_detail.tsv` when a count needs to be traced to
  the exact rank taxon and sequence IDs.
- Use `dataset_rank_novelty_detail.tsv` to inspect every created placeholder
  taxon, its parent, representative, and support.
- Use `dataset_increment_audit.md` for the human-readable per-dataset summary
  before inspecting detailed TSV rows.
- Use `source_contribution.tsv` to audit SILVA versus custom contributions.

Optimization items:

- Add a dedicated type-strain/name-legality report for SILVA init.
- Add a per-run diff report comparing before and after a dataset addition.

## 9. Validate Build

Command:

```text
autotax2 validate --build BUILD --strict
```

Once the single-step commands are behaving as expected, a normal incremental
dataset addition can use:

```text
autotax2 add --build BUILD --name D20 --prefix D20 --fasta D20.ssu.fa --domain Archaea --threads 8
```

Implementation:

- CLI: `autotax2/cli.py::validate`
- Core: `autotax2/validate.py::validate_build`

Checks:

- Required registry files exist.
- Dataset prefixes are unique and valid.
- Placeholder format and uniqueness are valid.
- Deprecated placeholders are not reused.
- Taxon parent-child ranks are valid and acyclic.
- Active species resolve to exactly seven ranks.
- Named SILVA immutability snapshot is stable.
- Sequence MD5 and mapping files are complete.
- Active representatives point to existing sequences.
- SILVA resolve evidence covers unresolved SILVA members at phylum through species.
- Placement evidence covers assignment rows at phylum through species.
- SINA candidate usage is recorded in cluster summary when candidate files are present.
- Export files validate when present.
- SINA/VSEARCH tool metadata is recorded.

Outputs:

| Path | Meaning |
|---|---|
| `reports/validation_report.md` | Human-readable validation report |
| `reports/validation_report.tsv` | Machine-readable findings |

Debug checklist:

- Run after every mutating step: `init`, `resolve`, `place`, and `export`.
- Use `--strict` when preparing a release or publication build.

## 10. Minimal Debug Run Order

For a small fixture or new real run, use this order:

```text
autotax2 init \
  --silva-fasta SILVA.fa.gz \
  --type-strain-metadata SILVA.full_metadata.tsv.gz \
  --gtdb-ar53-taxonomy ar53_taxonomy_r232.tsv \
  --gtdb-bac120-taxonomy bac120_taxonomy_r232.tsv \
  --outdir BUILD \
  --threads 8
autotax2 validate --build BUILD

autotax2 resolve --build BUILD --threads 8
autotax2 validate --build BUILD

autotax2 prepare --build BUILD --name D20 --prefix D20 --fasta D20.ssu.fa --domain Archaea
autotax2 orient --build BUILD --dataset D20 --threads 8
autotax2 cluster --build BUILD --dataset D20 --threads 8
autotax2 place --build BUILD --dataset D20 --dry-run

# inspect dry-run outputs
autotax2 place --build BUILD --dataset D20
autotax2 summarize --build BUILD --overwrite
autotax2 export all --build BUILD --force
autotax2 validate --build BUILD --strict
```

For real SINA/VSEARCH smoke tests, prefer:

```text
bash scripts/run_real_integration_test.sh --help
```

The script requires the official SILVA FASTA, SILVA full_metadata, both GTDB
r232 taxonomy files, and an externally extracted SSU/16S dataset FASTA.

## 11. Current Known Gaps

| Area | Current state | Proposed optimization |
|---|---|---|
| `add` command | Orchestrates prepare -> orient -> cluster -> place -> summarize -> validate, with optional export | Add resume/skip-step controls after real runs show which steps need restart support |
| SILVA unresolved resolver | Parent-scoped hierarchical resolve implemented; same-rank parent jobs are batched in parallel | Replace rank barriers with DAG scheduling if needed |
| SINA hit finding | Implemented as optional SINA `nearest_slv` candidate finding plus VSEARCH `--iddef 2` rescoring, with per-query candidate diagnostics | Use real runs to tune accession matching and fallback reporting |
| Strict legal-name catalog | Built from SILVA full_metadata type material plus required GTDB r232 taxonomy | Add optional LPSN valid-name input if needed |
| Invalid SILVA sequences | Non-type invalid rejected; type invalid retained | Add counts and retained-invalid warning report |
| Ambiguity reports | Present in assignments/consensus | Add clearer type-strain and multi-species ambiguity report |

## 12. Debug Evidence From Current Test Suite

Current verification command:

```text
$env:PYTHONPATH='.pytest_deps'; C:\Users\cheny\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe -m pytest
```

Latest local result during this debug pass:

```text
152 passed, 1 skipped
```
