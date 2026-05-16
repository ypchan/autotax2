# AutoTax Upstream Notes

These notes record how autotax2 maps the original AutoTax idea onto a
registry-driven incremental builder.

## Upstream AutoTax Model

AutoTax was designed to build ecosystem-specific full-length 16S rRNA gene
reference databases. Its public workflow is a monolithic Bash/R pipeline around
USEARCH, SINA, SILVA reference databases, and R data-frame joins.

The upstream flow can be read as:

1. Accept full-length 16S DNA sequences.
2. Optionally cluster input sequences at 99 percent before denovo taxonomy.
3. Append new FL-ASVs to a previously processed FASTA when `-d` is supplied.
4. Orient and align full-length sequences with SINA.
5. Cluster denovo representatives at descending rank thresholds.
6. Merge cluster tables from species toward phylum in R.
7. Fill unresolved SILVA taxonomy fields with denovo names.
8. Emit classifier-ready references, especially SINTAX and QIIME 2 formats.

Important upstream command-line options:

```text
-i  Input FASTA with full-length DNA sequences.
-c  Cluster generated FL-ASVs at 99 percent before denovo taxonomy and chimera filtering.
-d  Previously processed FL-ASV FASTA to append before rerunning denovo taxonomy.
-t  Maximum threads.
-b  Run BATS unit tests.
-v  Print version.
-h  Help.
```

The upstream design is useful, but its append mode reruns taxonomy over a
combined sequence set instead of preserving a durable per-addition registry.
It also relies primarily on USEARCH command behavior and R joins for the main
taxonomy build.

## autotax2 Interpretation

autotax2 keeps the biological idea but changes the state model:

1. SILVA named taxonomy becomes an immutable protected backbone.
2. SILVA unresolved records can become controlled placeholder taxa.
3. Each custom dataset gets a frozen prefix and internal sequence IDs.
4. Every addition writes dataset-local artifacts under `datasets/NN_name/`.
5. Registry files preserve all cumulative taxa, representatives, sequence MD5s,
   placeholder counters, and per-dataset membership.
6. Exact duplicate sequences are tracked but are not exported repeatedly.
7. VSEARCH replaces USEARCH for clustering and registry search.
8. Python owns orchestration, parsing, validation, reporting, and exports.
9. Placement uses near-best hit consensus rather than only the single best hit.
10. Reports describe each dataset's overlap with named SILVA, unresolved SILVA,
    previous custom additions, and newly created current-dataset taxa.

## autotax2 Default Rank Thresholds

```text
species: 0.972
genus:   0.901
family:  0.801
order:   0.729
class:   0.722
phylum:  0.696
```

During placement:

```text
identity >= species -> known_like
identity >= genus   -> new_species
identity >= family  -> new_genus
identity >= order   -> new_family
identity >= class   -> new_order
identity >= phylum  -> new_class
otherwise           -> unplaced
```
