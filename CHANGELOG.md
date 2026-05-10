# Changelog

## 0.1.0 - Unreleased

- Create initial Python package skeleton.
- Add Typer CLI with placeholder commands.
- Add config, taxonomy, registry, placeholder, and export helper stubs.
- Add placeholder allocation, dataset prefix registry, and taxonomy parser basics.
- Add Phase 2 FASTA IO, internal sequence ID assignment, MD5 de-duplication,
  and TSV registry writers.
- Add Phase 3 `autotax2 init` for SILVA FASTA parsing, optional domain
  filtering, named/unresolved taxonomy split, protected named taxon nodes, and
  initial registry outputs.
- Add Phase 4 `autotax2 resolve` for mutable SILVA placeholder framework
  creation, VSEARCH .uc parsing, placeholder reuse/deprecation protection, dry
  runs, and unresolved mapping outputs.
- Add Phase 5 `autotax2 prepare` for dataset prefix freezing, internal
  sequence IDs, normalized/prepared SSU FASTA outputs, MD5 duplicate
  membership, length checks, and preparation summaries from externally
  extracted SSU/16S input FASTA.
- Add Phase 6 `autotax2 orient` for loose SINA orientation correction,
  version recording, failure fallback, missing-output fallback, strand
  detection, and orientation summaries.
- Add Phase 7 `autotax2 cluster` for VSEARCH command wrapping, version
  recording, independent rank-threshold clustering, .uc membership parsing,
  current representative FASTA construction, registry search, coverage
  filtering, and cluster search summaries.
- Add Phase 8 `autotax2 place` for near-best hit consensus, identity/status
  placement decisions, dataset-specific placeholder creation, duplicate MD5
  handling, registry/representative updates, dry-run outputs, and placement
  summaries.
- Add Phase 9 `autotax2 export` for SINTAX, QIIME2, and DADA2 reference files,
  representative-only and all-unique modes, gzip FASTA output, 7-rank taxonomy
  validation, duplicate MD5 suppression, deprecated taxon skipping, and export
  manifests.
- Add Phase 10 `autotax2 summarize` and `autotax2 validate` for global build
  summaries, dataset delta reports, overlap matrices, source contributions,
  representative and de-duplication reports, registry validation, placeholder
  validation, taxonomy/tree checks, export checks, tool metadata checks, and
  Markdown/TSV validation reports.
- Add Phase 11 tiny fixtures, mocked end-to-end pipeline tests, full CLI help
  smoke coverage, SILVA unresolved parent-link cleanup, README/demo workflow
  documentation, AGENTS.md cleanup, and release checklist.
- Add Phase 12 optional real-tool integration support with a user-provided data
  workflow script, skipped-by-default pytest integration test, real-tool
  integration README, SINA/VSEARCH checks, and post-run export/report format
  checks.
- Expand README into a detailed command-by-command algorithm guide with input
  formats, parameter behavior, output file schemas, placement formulas, export
  contracts, and an AI-generated project overview image.
- Tighten early optimization semantics: initialize placeholder counter and
  snapshot registry files, rebuild representative search FASTA from durable
  sources when possible, use explicit prepared-SSU length/non-ATGC rejection,
  and fail loudly for currently unsupported reserved threshold controls.
- Reject SILVA records with empty or unresolved domain during initialization,
  write `silva_rejected.tsv`, parse SILVA full metadata gzip files for
  type-strain/type-material evidence, and prefer type-strain named SILVA
  representatives.
- Rewrite README as a text-only command-by-command tutorial with first-run,
  single-dataset, and multi-dataset workflows.
- Remove internal sequence extraction from `prepare`; dataset FASTA
  inputs are now explicitly expected to be externally processed SSU/16S
  sequences, and no internal extraction output files are produced.
- Add dated command audit logs under `logs/` and automatic export format
  self-check reports at `export/export_validation.tsv`.
- Add initial pytest coverage.
