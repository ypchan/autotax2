# Threading Optimization Notes

This note records the current threading bottlenecks and the proposed
accuracy-preserving optimization for AutoTax2.

## Current Bottlenecks

The expensive part of `autotax2 add` is the per-anchor-group vsearch
clustering step.

In the current MVP workflow:

1. SINA runs once for the whole dataset.
2. The corrected FASTA is parsed repeatedly while building group FASTA files.
3. Anchor groups are clustered one after another.
4. Each group may run several rank steps, for example species then genus.
5. Rank steps within one group are inherently sequential because centroids from
   the finer rank feed the next coarser rank.

The key observation is that different anchor groups are independent. Their
cluster input FASTA files are disjoint, their vsearch output paths are separate,
and their taxonomy anchor contexts are already fixed before clustering starts.

## Safe Parallel Boundary

Safe to parallelize:

- Different anchor groups during `autotax2 add`.
- Writing group-specific FASTA, UC, centroid, membership, and log files.

Not safe to parallelize without redesign:

- Rank steps inside one group, because `species.centroids.fa` feeds genus,
  genus feeds family, and so on.
- Placeholder allocation, because placeholder counters must remain deterministic
  and never reuse retired IDs.
- Updating global registries, because registry writes should stay ordered and
  auditable.

## Proposed CLI

Add one option to `autotax2 add`:

```bash
autotax2 add \
  --db autotax2_db \
  --input dataset.ssu.fa \
  --source midas \
  --prefix midas \
  --threads 64 \
  --group-jobs 4
```

`--threads` remains the total CPU budget. `--group-jobs` controls how many
independent anchor groups run concurrently.

Budget rule:

```text
workers = min(group_count, threads, group_jobs_or_auto)
threads_per_worker = max(1, floor(threads / workers))
```

This prevents the common mistake of launching many vsearch processes that each
request the full thread count.

## Accuracy Boundary

The optimization does not change:

- SINA inputs or outputs.
- identity thresholds.
- vsearch `--iddef 2`.
- anchor rank inference.
- rank order inside each group.
- placeholder allocation order.
- registry update order.
- export format semantics.

Only scheduling changes: independent groups are clustered concurrently, then
their summaries are sorted by numeric group index before writing.

## Patch

The concrete code patch is stored at:

```text
threading_optimization.patch
```

It changes:

- `autotax2/cli.py`
  - adds `--group-jobs`
  - passes it to `add_sequences`
- `autotax2/core.py`
  - parses the corrected FASTA once
  - builds an explicit list of cluster groups
  - runs independent groups with `ThreadPoolExecutor`
  - splits the thread budget across workers
  - writes per-group vsearch logs
  - returns deterministic group-sorted summaries

## Expected Speedup

Expected improvement depends on group distribution:

```text
one huge anchor group:
  small speedup; use vsearch threads inside the single group

many medium groups:
  large speedup; wall time can approach max(group time) rather than sum(group time)

many tiny groups:
  moderate speedup; process launch overhead becomes visible
```

Recommended starting values:

```text
--threads 32 --group-jobs 4
--threads 64 --group-jobs 4
--threads 128 --group-jobs 8
```

For a shared server, keep `threads * concurrent_autotax2_runs` below the actual
available cores.

## Verification Plan

Use a small fixture database and one synthetic dataset:

```bash
autotax2 add --db db --input test.fa --source test --prefix test --threads 8 --group-jobs 1
mv db/versions/v001_test db/versions/v001_test_serial

autotax2 add --db db2 --input test.fa --source test --prefix test --threads 8 --group-jobs 4
mv db2/versions/v001_test db2/versions/v001_test_parallel
```

Compare:

```bash
diff -u db/versions/v001_test_serial/sina_annotation.tsv \
        db2/versions/v001_test_parallel/sina_annotation.tsv

diff -u db/versions/v001_test_serial/provisional_taxonomy.tsv \
        db2/versions/v001_test_parallel/provisional_taxonomy.tsv

diff -u db/versions/v001_test_serial/cluster_summary.tsv \
        db2/versions/v001_test_parallel/cluster_summary.tsv
```

The exact vsearch UC file row order can vary between vsearch versions and
thread counts, so validation should focus on member-to-centroid equivalence and
final taxonomy equality rather than raw byte identity for every UC file.
