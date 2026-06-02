# AutoTax2

AutoTax2 builds a GTDB-anchored, incremental SSU/16S taxonomy reference. It
keeps the GTDB-derived backbone stable, then adds downstream SSU datasets by
orienting sequences with SINA, assigning a trusted GTDB anchor rank, clustering
novel sequence pools with vsearch, and exporting classifier-ready references.

## Workflow Map

![AutoTax2 GTDB-anchored workflow](docs/images/autotax2-workflow.png)

The workflow figure shows the two main phases: `autotax2 init` builds the
GTDB-backed database, and `autotax2 add` aligns each new dataset with SINA,
infers trusted anchor ranks, clusters novel groups with vsearch, updates
registries, and writes exports.

## Core Ideas

- GTDB-derived SSU taxonomy is the named backbone.
- SINA runs against a user-built `gtdb_ssu.arb` file to orient query sequences
  and emit alignment, identity, orientation, and LCA metadata.
- vsearch clusters only the novel lineage space below the trusted GTDB anchor
  rank.
- Placeholder names are global within each rank and carry the first source
  prefix, for example `s__midas_s000001`.
- Dataset FASTA inputs must already be extracted SSU/16S sequences. AutoTax2
  does not extract rRNA genes from genomes or contigs.

## Installation

```bash
git clone https://github.com/ypchan/autotax2.git
cd autotax2
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Check the package and external tools:

```bash
autotax2 --help
sina --version
vsearch --version
pytest
ruff check autotax2 tests
```

Recommended external tool versions:

- SINA 1.7.2
- vsearch 2.30.6 Linux x86_64 or newer compatible vsearch 2.x

## SINA Command Requirements

AutoTax2 uses SINA for three things:

1. orientation correction;
2. reference identity and alignment quality calculation;
3. LCA taxonomy metadata written into FASTA headers.

A minimal alignment command is not enough. The command must request search,
identity calculation, orientation reporting, LCA fields, DNA FASTA output, and
header metadata. A representative command is:

```bash
sina \
  -i bac.16s_rRNA.single_copy.fa \
  --db /path/to/gtdb_ssu.arb \
  --search \
  --search-min-sim 0.5 \
  --search-max-result 10 \
  --lca-fields=tax_slv,tax_gtdb \
  --lca-quorum 0.7 \
  --show-conf \
  --threads 48 \
  --fasta-write-dna \
  --turn all \
  --calc-idty \
  --fs-req 1 \
  --fs-req-full 0 \
  --fs-msc 0.5 \
  --meta-fmt header \
  -o bac.16s_rRNA.single_copy.sina.fa
```

Expected FASTA header shape:

```text
>ERR14064717_s_ERR14064717.1773379 [align_cutoff_head_slv=701] [align_cutoff_tail_slv=5] [align_ident_slv=60.1063843] [align_quality_slv=66] [aligned_slv=2026-05-22 10:58:37] [lca_tax_gtdb=Unclassified;] [lca_tax_slv=Unclassified;] [turn=reversed and complemented]
```

SINA may keep `_slv` suffixes even when the database is your custom
`gtdb_ssu.arb`. In that case AutoTax2 treats `align_ident_slv`,
`align_quality_slv`, `align_cutoff_head_slv`, and `align_cutoff_tail_slv` as
reference metrics from the GTDB ARB run. Normalized fields such as
`align_ident_ref`, `align_quality_ref`, and `lca_tax_ref` are also accepted.

Important SINA options:

| option | why it matters |
|---|---|
| `--db gtdb_ssu.arb` | Uses the custom GTDB ARB reference database. |
| `--search` | Runs reference search so LCA fields can be calculated. |
| `--search-min-sim 0.5` | Keeps distant candidates for high-rank anchoring. |
| `--search-max-result 10` | Limits candidates used for LCA. |
| `--lca-fields=tax_slv,tax_gtdb` | Requests taxonomy fields from ARB entries. |
| `--lca-quorum 0.7` | Requires 70 percent agreement for reported LCA. |
| `--show-conf` | Includes confidence metadata when SINA can report it. |
| `--fasta-write-dna` | Writes DNA FASTA output. |
| `--turn all` | Allows reverse-complement detection and reports `[turn=...]`. |
| `--calc-idty` | Calculates `align_ident_*`, required for anchor rank inference. |
| `--fs-req 1` | Enables flexible search behavior. |
| `--fs-req-full 0` | Avoids requiring full-length search behavior. |
| `--fs-msc 0.5` | Sets flexible-search minimum similarity. |
| `--meta-fmt header` | Writes metadata into FASTA headers for AutoTax2 parsing. |

If identity, quality, cutoff, LCA, or orientation fields are missing, first
check that `--calc-idty`, `--turn all`, `--lca-fields=...`, and
`--meta-fmt header` were used.

## Custom GTDB ARB Taxonomy Fields

The custom `gtdb_ssu.arb` must contain per-sequence taxonomy fields if SINA is
expected to emit `lca_tax_gtdb` or `lca_tax_slv`. The ARB topology and aligned
sequences alone are not enough. If the ARB entries do not have taxonomy
annotation fields, SINA can still orient and align sequences, but LCA taxonomy
headers will be empty or `Unclassified`.

Recommended input table:

```text
seq_id    tax_gtdb
GCF_000001.1_16S    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__...;f__...;g__...;s__...
```

Recommended ARB fields:

```text
tax_gtdb = d__Bacteria;p__...;s__...
tax_slv  = d__Bacteria;p__...;s__...
```

`tax_gtdb` is the semantically correct field for AutoTax2. Mirroring the same
taxonomy into `tax_slv` is useful because SINA examples and older metadata names
often expect `tax_slv`.

After building or post-processing the ARB, run SINA with:

```bash
--lca-fields=tax_slv,tax_gtdb --meta-fmt header
```

Then confirm the output FASTA contains:

```text
[lca_tax_gtdb=d__Bacteria;p__...;]
```

The critical invariant is that ARB sequence IDs, reference FASTA IDs, and
taxonomy TSV `seq_id` values match exactly.

## Quick Start

### 1. Check Dependencies

```bash
autotax2 check --debug
```

| parameter | default | meaning |
|---|---:|---|
| `--sina-bin` | `sina` | SINA executable name or path. |
| `--vsearch-bin` | `vsearch` | vsearch executable name or path. |
| `--debug` | `False` | Enable debug logging. |

### 2. Initialize a GTDB Backbone Database

```bash
autotax2 init \
  --ref-fa gtdb_ssu_tax.fa \
  --ref-tax gtdb_ssu.taxonomy.tsv \
  --ref-arb gtdb_ssu.arb \
  --db autotax2_db \
  --threads 64 \
  --debug
```

| parameter | required | meaning |
|---|---:|---|
| `--ref-fa` | yes | GTDB-derived SSU reference FASTA. IDs must match `seq_id`. |
| `--ref-tax` | yes | GTDB taxonomy TSV with one row per reference sequence. |
| `--ref-arb` | yes | SINA-compatible ARB file with taxonomy fields. |
| `--db` | yes | Output AutoTax2 database directory. |
| `--threads`, `-t` | no | CPU thread budget recorded in `config.yaml`. |
| `--force` | no | Rebuild an existing non-empty database directory. |
| `--debug` | no | Enable debug logging. |

Required taxonomy TSV columns:

```text
seq_id    domain    phylum    class    order    family    genus    species
```

Important outputs:

```text
autotax2_db/config.yaml
autotax2_db/reference/current.fa
autotax2_db/reference/current.taxonomy.tsv
autotax2_db/registry/sequence_registry.tsv
autotax2_db/registry/cluster_registry.tsv
autotax2_db/registry/placeholder_counter.tsv
autotax2_db/registry/source_registry.tsv
autotax2_db/clusters/<rank>.membership.tsv
autotax2_db/logs/init.log
```

### 3. Add a Dataset

```bash
autotax2 add \
  --db autotax2_db \
  --input MiDAS5.3.ssu.fa \
  --source midas \
  --prefix midas \
  --threads 64 \
  --group-jobs 4 \
  --mode incremental \
  --debug
```

![AutoTax2 add workflow](docs/images/autotax2-add-workflow.png)

| parameter | required | meaning |
|---|---:|---|
| `--db` | yes | Existing AutoTax2 database directory. |
| `--input` | yes | New dataset FASTA containing externally extracted SSU/16S sequences. |
| `--source` | yes | Source name, for example `midas`, `mfd`, or `hifimeta`. |
| `--prefix` | yes | Stable placeholder prefix for new lineages created by this source. |
| `--threads`, `-t` | no | Total CPU thread budget for SINA, dereplication, and clustering. |
| `--group-jobs` | no | Independent anchor groups to cluster concurrently. |
| `--mode` | no | `incremental` or `full`; default is `incremental`. |
| `--keep-temp` | no | Keep intermediate files for debugging. |
| `--dry-run` | no | Print external commands without running downstream generation. |
| `--debug` | no | Enable detailed logs. |

Important dataset outputs:

```text
autotax2_db/versions/v001_midas/01_sina/new.sina.aligned.fa
autotax2_db/versions/v001_midas/01_sina/new.sina.corrected.fa
autotax2_db/versions/v001_midas/sina_annotation.tsv
autotax2_db/versions/v001_midas/02_derep/new.derep.fa
autotax2_db/versions/v001_midas/02_derep/new.derep.uc
autotax2_db/versions/v001_midas/03_clusters/group_*.fa
autotax2_db/versions/v001_midas/03_clusters/group_*/*.uc
autotax2_db/versions/v001_midas/03_clusters/group_*.membership.tsv
autotax2_db/versions/v001_midas/cluster_summary.tsv
autotax2_db/versions/v001_midas/provisional_taxonomy.tsv
autotax2_db/logs/add_midas.log
```

## Anchor and Placeholder Logic

AutoTax2 picks the finest trusted anchor rank from SINA identity and GTDB LCA
taxonomy.

| rank | identity cutoff |
|---|---:|
| species | 97.2 |
| genus | 90.1 |
| family | 80.1 |
| order | 72.9 |
| class | 72.2 |
| phylum | 69.6 |

Example:

```text
SINA LCA: d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__...;f__...
SINA identity: 85.0
trusted anchor: family
novel ranks: genus, species
```

Placeholder IDs are allocated globally by rank:

```text
s__midas_s000001
s__mfd_s000002
g__midas_g000001
g__hifimeta_g000002
```

## Export References

```bash
autotax2 export \
  --db autotax2_db \
  --format all \
  --out export/
```

| parameter | required | meaning |
|---|---:|---|
| `--db` | yes | Existing AutoTax2 database directory. |
| `--format` | no | `all`, `sintax`, `dada2`, `qiime2`, or `taxonomy`. |
| `--out` | yes | Output file or directory. |
| `--debug` | no | Enable export logs. |

Export outputs:

```text
export/autotax2.sintax.fa
export/autotax2.dada2.fa
export/autotax2.qiime2-seqs.fa
export/autotax2.qiime2-taxonomy.tsv
export/autotax2.taxonomy.tsv
```

## Threading and Performance

`--threads` is the total thread budget. SINA and dereplication usually run as
one external process and receive the full budget.

During `add`, independent anchor groups can be clustered concurrently because
their input FASTA files and output paths are separate. Rank steps inside one
group remain sequential because species centroids feed genus clustering, genus
centroids feed family clustering, and so on.

Recommended strategy:

```text
few huge anchor groups:
  --group-jobs 1 or 2

many small or medium anchor groups:
  --group-jobs 4 to 16 on large servers

shared server:
  lower --threads or --group-jobs to avoid CPU oversubscription
```

Budget rule:

```text
group_jobs <= threads
threads_per_group = max(1, floor(threads / group_jobs))
```

## Development

```bash
python -m pip install -e ".[dev]"
pytest
ruff check autotax2 tests
```
