# AGENTS.md

This file records project rules for future agents working on autotax2.

## Project Goal

autotax2 is a fixed-backbone, rank-aware, incremental rRNA gene reference
builder. It must preserve a stable named SILVA backbone while allowing controlled
incremental additions from downstream datasets.

## Design Rules

1. SILVA named backbone is immutable.
2. SILVA unresolved records may form mutable placeholder framework.
3. Placeholder format:
   - `g__SILVAg000001`
   - `s__SILVAs000001`
   - `g__D20g000001`
   - `s__D20s000001`
4. Deprecated placeholder IDs are never reused.
5. Dataset prefix is supplied during `autotax2 prepare` and then frozen.
6. Input sequence IDs are remapped to internal IDs like `D20_000001`.
7. Sequence MD5 is stored; exact duplicate sequences are not exported repeatedly.
8. VSEARCH identity uses fixed `--iddef 2` by default.
9. Dataset FASTA inputs are already externally extracted SSU/16S sequences;
   autotax2 does not run sequence extraction internally.
10. SINA orientation should use loose settings.
11. Placement uses near-best hit consensus, not best hit only.
12. Exports must include:
    - SINTAX reference
    - QIIME2 reference sequences + taxonomy
    - DADA2 toGenus trainset
    - DADA2 assignSpecies file
13. Summary and validation reports must preserve named SILVA immutability,
    duplicate-sequence accounting, placeholder uniqueness, and export format
    compatibility.
14. Major CLI runs should leave dated audit logs under `logs/`.

## Engineering Notes

- Keep the initial algorithmic modules small and explicit until the contracts
  become clear.
- Prefer deterministic IDs and durable registries over inferred state.
- Treat exported references as reproducible build artifacts.
- Add tests for behavior before expanding the algorithm.
