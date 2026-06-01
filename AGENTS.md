# AGENTS.md

This file records project rules for future agents working on autotax2.

## Project Goal

autotax2 is a GTDB-anchored, rank-aware, incremental SSU/16S reference
builder. It preserves a stable GTDB-derived named backbone while allowing
controlled incremental additions from downstream datasets.

The current workflow uses a user-built `gtdb_ssu.arb` reference with SINA for
query orientation and anchor inference, and vsearch for local hierarchical
clustering of novel lineages below the trusted GTDB anchor rank.

## Design Rules

1. The GTDB-derived named backbone is immutable after `autotax2 init`.
2. GTDB backbone taxonomy is supplied by `--ref-fa`, `--ref-tax`, and
   `--ref-arb`; these inputs define the stable reference state.
3. New datasets are added with `autotax2 add`; the dataset `--source` and
   placeholder `--prefix` are registered and must remain stable for that
   source.
4. Dataset FASTA inputs are already externally extracted SSU/16S sequences;
   autotax2 does not run sequence extraction internally.
5. SINA is used against `gtdb_ssu.arb` for orientation, alignment metrics, and
   GTDB LCA taxonomy fields.
6. SINA header fields ending in `_slv` may still be emitted by SINA; when SINA
   is run against the GTDB ARB, autotax2 treats those fields as reference
   metrics. Normalized `_ref` and `_gtdb` fields should also be accepted.
7. Anchor rank inference uses GTDB LCA taxonomy plus fixed identity thresholds:
   species 97.2, genus 90.1, family 80.1, order 72.9, class 72.2, phylum 69.6.
8. New lineages are clustered only below the trusted anchor rank.
9. vsearch identity uses fixed `--iddef 2` by default.
10. Placeholder IDs are global within each rank, not within each source.
11. Placeholder format is rank-prefixed and source-prefixed, for example:
    - `s__midas_s000001`
    - `g__midas_g000001`
    - `s__mfd_s000002`
    - `g__hifimeta_g000002`
12. Deprecated or retired placeholder IDs are never reused.
13. Input sequence MD5 values are stored; exact duplicate sequences should not
    be exported repeatedly.
14. Cluster membership should be retained for every rank to support
    source-specific and shared-lineage summaries.
15. Exports must include:
    - vsearch SINTAX reference
    - QIIME2 reference sequences and taxonomy
    - DADA2-compatible training FASTA
    - flat taxonomy TSV
16. Summary and validation reports must preserve named GTDB immutability,
    duplicate-sequence accounting, placeholder uniqueness, and export format
    compatibility.
17. Major CLI runs should leave audit logs under `logs/`.

## Current CLI Shape

- `autotax2 check` verifies external `sina` and `vsearch` executables.
- `autotax2 init` initializes a GTDB-backed database.
- `autotax2 add` incrementally adds a new dataset.
- `autotax2 rebuild` runs a full hierarchical vsearch rebuild.
- `autotax2 export` writes annotation references.
- `autotax2 summarize` reports source overlap by rank.

## Engineering Notes

- Keep the initial algorithmic modules small and explicit until the contracts
  become clear.
- Prefer deterministic IDs and durable registries over inferred state.
- Treat exported references as reproducible build artifacts.
- Add tests for behavior before expanding the algorithm.
- When changing user-facing behavior, update README examples and CLI help
  together.
- When changing registry schemas, preserve upgrade paths for existing database
  directories whenever practical.
