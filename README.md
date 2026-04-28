# AutoTax2

AutoTax2 是一个以 **VSEARCH** 为核心、以 **SILVA 为 backbone** 的 near/full-length rRNA 序列分类、分层聚类、参考库扩充和 overlap 分析工具。它的设计目标是替代原 AutoTax 中对 USEARCH 的依赖，同时避免把 SILVA 已有分类体系重新打散。

AutoTax2 的核心原则是：

```text
1. SILVA 是固定 backbone，不重新聚类、不打破原有分类体系。
2. 用户数据、ref2/ref3/ref4 等 extension datasets 先尝试归入 SILVA。
3. 能归入 SILVA 的序列使用 SILVA taxon_id。
4. 不能归入 SILVA 的 extension 序列才在对应 rank 阈值下聚类，形成 novel taxa。
5. VSEARCH 能完成的任务直接交给 VSEARCH，Python 只做格式转换、流程 glue logic 和结果汇总。
6. 用户自己提供本地参考数据；AutoTax2 不自动下载 SILVA。
7. 对含 intron 的 16S，可先生成 intron-free analysis sequence 参与计算，但正式 centroid 仍输出 original sequence。
```

---

## 目录

- [1. 适用场景](#1-适用场景)
- [2. 核心概念](#2-核心概念)
- [3. 仓库结构](#3-仓库结构)
- [4. 安装](#4-安装)
- [5. 输入数据准备](#5-输入数据准备)
- [6. 推荐主流程](#6-推荐主流程)
- [7. 命令总览](#7-命令总览)
- [8. `check`](#8-check)
- [9. `prepare-silva`](#9-prepare-silva)
- [10. `detect-intron`](#10-detect-intron)
- [11. `insert-backbone`](#11-insert-backbone)
- [12. `overlap-backbone`](#12-overlap-backbone)
- [13. 传统 VSEARCH 辅助命令](#13-传统-vsearch-辅助命令)
- [14. 输出文件详解](#14-输出文件详解)
- [15. 下游分析数据怎么用](#15-下游分析数据怎么用)
- [16. 针对 HiFiMeta 500 万 16S 的建议](#16-针对-hifimeta-500-万-16s-的建议)
- [17. 参数调优建议](#17-参数调优建议)
- [18. 已知限制](#18-已知限制)
- [19. 开发路线](#19-开发路线)

---

## 1. 适用场景

AutoTax2 特别适合以下场景：

```text
你有大量 near full-length / full-length 16S rRNA 序列；
你希望基于 SILVA 138.2 这样的成熟参考分类体系进行归类；
你不希望把 SILVA 和自己的序列全部混合后完全按序列相似性重新聚类；
你希望知道自己的数据在不同 rank-like 阈值下有多少已知类群和潜在新类群；
你有 ref2/ref3/ref4 等额外参考库，希望看它们相对于 SILVA 和自己的数据集有什么补充；
你担心部分 16S 中含有 intron，导致全局比对 identity 偏低，需要 intron-aware rescue。
```

一个典型应用是：

```text
1000 个 HiFiMeta 样本
  ↓
barrnap 提取 near/full-length 16S rRNA
  ↓
100% identity + 90% qcov 去冗余
  ↓
得到约 500 万条 16S 序列
  ↓
基于 SILVA backbone 进行 rank-wise insertion
  ↓
统计 species/genus/family/order/class/phylum-like 层面的 known/novel taxa
  ↓
与 ref2/ref3/ref4 进行 SILVA-backed overlap 分析
```

---

## 2. 核心概念

### 2.1 SILVA backbone

SILVA 已经有自己的 taxonomy。AutoTax2 不会把 SILVA 序列和用户序列混在一起重新聚类，因为这样会破坏已有分类体系。

AutoTax2 使用 SILVA 作为 backbone：

```text
输入序列 vs SILVA
  ├── 达到某 rank 阈值 → 归入已有 SILVA taxon
  └── 未达到阈值 → extension 序列内部聚类，生成 novel taxon
```

### 2.2 source-label

`source-label` 是用户自定义的数据来源标签。它不是必须代表样本，也可以代表数据库来源。

例如：

```text
hifimeta
ref2
ref3
ref4
```

运行时：

```bash
--source-label hifimeta
```

这表示当前输入 FASTA 中所有序列都标记为 `hifimeta`。

### 2.3 rank thresholds

默认 rank 阈值来自常用 16S novelty threshold：

| rank | identity threshold | 含义 |
|---|---:|---|
| species | 0.987 | 低于此阈值，可能是新 species-like taxon |
| genus | 0.945 | 低于此阈值，可能是新 genus-like taxon |
| family | 0.865 | 低于此阈值，可能是新 family-like taxon |
| order | 0.820 | 低于此阈值，可能是新 order-like taxon |
| class | 0.785 | 低于此阈值，可能是新 class-like taxon |
| phylum | 0.750 | 低于此阈值，可能是新 phylum-like taxon |

AutoTax2 输出中使用 `species-like`、`genus-like` 等标签，避免把 operational threshold 误解为严格分类学定义。

### 2.4 analysis sequence vs original sequence

为处理 intron，AutoTax2 区分两种序列：

| 类型 | 说明 |
|---|---|
| original sequence | 原始输入序列，正式 centroid 输出应使用这个 |
| analysis sequence | 计算用序列；如果检测到 intron，则为 intron-free sequence |

如果没有 intron：

```text
analysis sequence = original sequence
```

如果检测到 intron：

```text
analysis sequence = original sequence 去掉 intron 后的序列
```

后续：

```text
分类 / 聚类 / 插入 backbone 使用 analysis sequence
正式 centroid FASTA 使用 original sequence
```

### 2.5 core centroid vs original centroid

`insert-backbone` 输出两套 centroids：

```text
rank_centroids_core/
  计算用 centroid，可能是 intron-free sequence

rank_centroids_original/
  正式输出 centroid，使用 original sequence
```

---

## 3. 仓库结构

```text
autotax2/
  README.md
  pyproject.toml
  environment.yml
  LICENSE
  .gitignore

  autotax2/
    __init__.py
    cli.py
    utils.py
    logging.py
    threads.py
    dependencies.py
    vsearch.py
    prepare.py
    ranks.py
    backbone.py
    overlap.py
    intron.py
    assign.py
    provenance.py
    summarize.py

  scripts/
    detect_intron.py

  examples/
    data/
      input.fa
      new.fa
      old_centroids_97.fa
      source_map.tsv
    ref/
      silva_tiny.fasta
      silva_tiny.full_metadata

  docs/
    FORMAT.md

  tests/
    test_prepare.py
    test_assign.py
    test_provenance.py
    test_threads.py
    test_ranks.py
    test_intron.py

  .github/
    workflows/
      ci.yml
```

### 3.1 主要 Python 模块

| 文件 | 作用 |
|---|---|
| `cli.py` | 命令行入口 |
| `prepare.py` | 本地 SILVA FASTA / metadata 预处理 |
| `vsearch.py` | VSEARCH 命令封装 |
| `ranks.py` | rank 阈值体系 |
| `backbone.py` | SILVA backbone 插入式分类 / 聚类 |
| `overlap.py` | 基于 backbone assignment 的 ref overlap |
| `intron.py` | intron 检测与 intron-free analysis FASTA 生成 |
| `assign.py` | 新增序列分配到已有簇或生成新簇 |
| `provenance.py` | 来源组成和 UC overlap 统计 |
| `summarize.py` | SINTAX / SILVA / type strain 分类证据汇总 |
| `dependencies.py` | 依赖检查 |
| `threads.py` | 线程解析与 capped 策略 |
| `logging.py` | Rich 日志、表格、进度条封装 |

---

## 4. 安装

### 4.1 使用 conda / mamba

```bash
mamba env create -f environment.yml
conda activate autotax2
```

### 4.2 在已有环境中安装

```bash
mamba install -c conda-forge -c bioconda vsearch seqkit rich rich-argparse python
pip install -e .
```

### 4.3 检查安装

```bash
autotax2 --help
autotax2 check
vsearch --version
```

---

## 5. 输入数据准备

### 5.1 SILVA FASTA

用户需要提前下载 SILVA FASTA，例如：

```text
SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
```

FASTA header 应类似：

```text
>AB000001.1.1500 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus acidophilus;
ACGT...
```

### 5.2 SILVA metadata

用户需要提前下载 SILVA metadata，例如：

```text
SILVA_138.2_SSURef.full_metadata.gz
```

必须包含：

```text
acc
flags
```

`flags` 中带 `[T]` 或 `[t]` 的 accession 会被识别为 type strain。

### 5.3 用户输入 16S FASTA

例如：

```text
hifimeta.original.fa
```

推荐 FASTA header 第一列是唯一序列 ID：

```text
>seq000001
ACGT...
>seq000002
ACGT...
```

如果 header 是：

```text
>seq000001 some description
```

AutoTax2 使用第一个空格前的 `seq000001` 作为 ID。

### 5.4 source-label

如果整个 FASTA 都属于一个来源：

```bash
--source-label hifimeta
```

如果你要更细粒度到样本、宿主、批次，可以另外使用 `source_map.tsv` 给 `provenance` 或旧流程使用，但 backbone 主流程推荐用 `source-label` 表示 dataset-level source。

---

## 6. 推荐主流程

下面是当前最推荐的主流程。

### Step 1. 准备 SILVA backbone

```bash
autotax2 prepare-silva \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --silva-metadata SILVA_138.2_SSURef.full_metadata.gz \
  --out refdatabases \
  --make-udb \
  --threads auto
```

输出：

```text
refdatabases/
  SILVA_138.2_SSURef_NR99_tax_silva.fasta
  SILVA_138.2_SSURef_NR99_tax_silva.udb
  SILVA_138.2_SSURef_NR99_tax_silva_sintax.fasta
  SILVA_138.2_SSURef_NR99_tax_silva_sintax.udb
  SILVA_138.2_SSURef_NR99_tax_silva_typestrains.fasta
  SILVA_138.2_SSURef_NR99_tax_silva_typestrains.udb
  SILVA_138.2_SSURef.full_metadata
  typestrains_accessionIDs.txt
  silva_taxonomy.tsv
  autotax2_ref_manifest.tsv
```

### Step 2. 可选：检测 intron

如果你怀疑部分 16S 有 intron：

```bash
autotax2 detect-intron \
  --input hifimeta.original.fa \
  --db refdatabases/SILVA_138.2_SSURef_NR99_tax_silva.udb \
  --source-label hifimeta \
  --out hifimeta_intron \
  --search-id 0.70 \
  --rescue-id 0.987 \
  --min-intron-len 50 \
  --min-flank-len 150 \
  --threads auto
```

后续用：

```text
hifimeta_intron/analysis_sequences.fa
hifimeta_intron/sequence_version_map.tsv
```

### Step 3. SILVA backbone 插入

如果做了 intron 检测：

```bash
autotax2 insert-backbone \
  --input hifimeta_intron/analysis_sequences.fa \
  --original-fasta hifimeta.original.fa \
  --version-map hifimeta_intron/sequence_version_map.tsv \
  --source-label hifimeta \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --rank-thresholds default \
  --db-format auto \
  --out hifimeta_inserted \
  --threads auto
```

如果不做 intron 检测：

```bash
autotax2 insert-backbone \
  --input hifimeta.original.fa \
  --original-fasta hifimeta.original.fa \
  --source-label hifimeta \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --rank-thresholds default \
  --db-format auto \
  --out hifimeta_inserted \
  --threads auto
```

### Step 4. 多参考 overlap

分别处理 ref2/ref3/ref4：

```bash
autotax2 insert-backbone \
  --input ref2.fa \
  --source-label ref2 \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --out ref2_inserted \
  --threads auto

autotax2 insert-backbone \
  --input ref3.fa \
  --source-label ref3 \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --out ref3_inserted \
  --threads auto

autotax2 insert-backbone \
  --input ref4.fa \
  --source-label ref4 \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --out ref4_inserted \
  --threads auto
```

然后：

```bash
autotax2 overlap-backbone \
  --assignments \
    ref2_inserted/sequence_rank_assignment.tsv \
    ref3_inserted/sequence_rank_assignment.tsv \
    ref4_inserted/sequence_rank_assignment.tsv \
    hifimeta_inserted/sequence_rank_assignment.tsv \
  --labels ref2,ref3,ref4,hifimeta \
  --out ref_overlap
```

---

## 7. 命令总览

```bash
autotax2 check
autotax2 prepare-silva
autotax2 detect-intron
autotax2 insert-backbone
autotax2 overlap-backbone

# 辅助 / 兼容命令
autotax2 derep
autotax2 cluster
autotax2 classify
autotax2 assign
autotax2 provenance
autotax2 summarize
autotax2 run
```

每个命令都支持：

```bash
autotax2 <command> --example
```

例如：

```bash
autotax2 detect-intron --example
autotax2 insert-backbone --example
autotax2 overlap-backbone --example
```

---

## 8. `check`

检查依赖和参考文件。

### 命令

```bash
autotax2 check
```

检查参考库：

```bash
autotax2 check \
  --ref-manifest refdatabases/autotax2_ref_manifest.tsv
```

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--vsearch` | `vsearch` | VSEARCH 路径 |
| `--ref-manifest` | 无 | 检查 manifest 中的参考文件 |
| `--source-map` | 无 | 检查 source map 文件 |

---

## 9. `prepare-silva`

处理本地 SILVA 文件，不下载数据。

### 命令

```bash
autotax2 prepare-silva \
  --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --silva-metadata SILVA_138.2_SSURef.full_metadata.gz \
  --out refdatabases \
  --make-udb \
  --threads auto
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--silva-fasta` | 是 | 无 | 本地 SILVA FASTA，可 `.gz` |
| `--silva-metadata` | 是 | 无 | 本地 SILVA full_metadata，可 `.gz` |
| `--out` | 否 | `refdatabases` | 输出目录 |
| `--prefix` | 否 | 从 FASTA 文件名推断 | 输出文件前缀 |
| `--make-udb` | 否 | False | 是否用 VSEARCH 生成 UDB |
| `--vsearch` | 否 | `vsearch` | VSEARCH 路径 |
| `--threads` | 否 | `auto` | 线程数 |
| `--dry-run` | 否 | False | 只打印命令 |

### 输出

| 文件 | 说明 |
|---|---|
| `*.fasta` | 解压/复制后的 SILVA FASTA |
| `*.udb` | 可选，VSEARCH UDB |
| `*_sintax.fasta` | SINTAX 格式 FASTA |
| `*_sintax.udb` | 可选，SINTAX UDB |
| `*_typestrains.fasta` | type strain FASTA |
| `*_typestrains.udb` | 可选，type strain UDB |
| `typestrains_accessionIDs.txt` | 从 metadata 中提取的 type strain IDs |
| `silva_taxonomy.tsv` | 从 SILVA header 解析出的 taxonomy |
| `autotax2_ref_manifest.tsv` | 后续步骤引用的 manifest |

---

## 10. `detect-intron`

检测 query 相对于 SILVA top subject 的长插入，生成 intron-free analysis FASTA。

### 命令

```bash
autotax2 detect-intron \
  --input hifimeta.original.fa \
  --db refdatabases/SILVA_138.2_SSURef_NR99_tax_silva.udb \
  --source-label hifimeta \
  --out hifimeta_intron \
  --search-id 0.70 \
  --rescue-id 0.987 \
  --min-intron-len 50 \
  --min-flank-len 150 \
  --threads auto
```

### 独立脚本

```bash
python scripts/detect_intron.py \
  --input hifimeta.original.fa \
  --db refdatabases/SILVA_138.2_SSURef_NR99_tax_silva.udb \
  --source-label hifimeta \
  --out hifimeta_intron \
  --threads auto
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--input` | 是 | 无 | 原始 rRNA FASTA |
| `--db` | 是 | 无 | SILVA FASTA 或 UDB |
| `--out` | 是 | 无 | 输出目录 |
| `--source-label` | 否 | `query` | 数据来源标签 |
| `--search-id` | 否 | 0.70 | 初始 VSEARCH 搜索 identity |
| `--rescue-id` | 否 | 0.987 | 去掉 candidate intron 后要求的 core identity |
| `--min-intron-len` | 否 | 50 | 最小插入长度 |
| `--min-flank-len` | 否 | 150 | 插入左右两侧最小 core 长度 |
| `--strand` | 否 | `both` | VSEARCH strand |
| `--maxaccepts` | 否 | 1 | 每条 query 保留 top hit 数 |
| `--threads` | 否 | `auto` | 线程数 |
| `--dry-run` | 否 | False | 只打印命令 |

### 检测逻辑

`detect-intron` 要求 VSEARCH 输出：

```text
query target id alnlen mism opens qlo qhi tlo thi qstrand qrow trow
```

然后找：

```text
query 有连续碱基
target 对应位置是连续 gap
```

如果该 query insertion：

```text
长度 >= min-intron-len
去掉后 core identity >= rescue-id
左右 flank >= min-flank-len
```

则认为：

```text
status = intron_rescued
```

### 输出

| 文件 | 说明 |
|---|---|
| `vsearch_alignment_rows.tsv` | VSEARCH pairwise alignment row 输出 |
| `vsearch_alignment.uc` | VSEARCH UC |
| `analysis_sequences.fa` | 后续计算用 FASTA |
| `sequence_version_map.tsv` | original 与 analysis sequence 对应关系 |
| `intron_summary.tsv` | intron 检测汇总表 |
| `intron_sequences.fa` | intron 序列本身 |
| `intron_regions.bed` | intron 坐标，0-based BED |

### `sequence_version_map.tsv`

| 列 | 说明 |
|---|---|
| `sequence_id` | 原始序列 ID |
| `analysis_sequence_id` | 计算用序列 ID |
| `source_label` | 来源标签 |
| `has_intron` | yes/no |
| `original_length` | 原始长度 |
| `analysis_length` | 去 intron 后长度 |
| `centroid_sequence_allowed` | 当前固定为 `original_only` |

### `intron_summary.tsv`

| 列 | 说明 |
|---|---|
| `sequence_id` | 原始序列 ID |
| `source_label` | 来源 |
| `subject_id` | top subject |
| `original_identity` | 原始比对 identity |
| `original_length` | 原始长度 |
| `analysis_sequence_id` | analysis sequence ID |
| `analysis_length` | analysis sequence 长度 |
| `has_intron` | yes/no |
| `n_introns` | intron 数量 |
| `intron_positions_query` | query 上 1-based 坐标，如 `720-909` |
| `intron_lengths` | intron 长度 |
| `rescued_identity` | 去掉 intron 后 core identity |
| `status` | `intron_rescued`, `candidate_insertion_not_rescued`, `no_intron_detected`, `no_hit` |

---

## 11. `insert-backbone`

核心命令。将 extension dataset 插入固定 SILVA backbone。

### 命令

```bash
autotax2 insert-backbone \
  --input hifimeta_intron/analysis_sequences.fa \
  --original-fasta hifimeta.original.fa \
  --version-map hifimeta_intron/sequence_version_map.tsv \
  --source-label hifimeta \
  --silva-manifest refdatabases/autotax2_ref_manifest.tsv \
  --rank-thresholds default \
  --db-format auto \
  --out hifimeta_inserted \
  --threads auto
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--input` | 是 | 无 | analysis FASTA |
| `--original-fasta` | 否 | input | 原始 FASTA；正式 centroid 使用此文件 |
| `--version-map` | 否 | 无 | detect-intron 生成的 version map |
| `--source-label` | 是 | 无 | 数据来源标签 |
| `--silva-manifest` | 是 | 无 | prepare-silva 生成的 manifest |
| `--rank-thresholds` | 否 | default | rank 阈值 |
| `--db-format` | 否 | auto | `auto`, `udb`, `fasta` |
| `--maxaccepts` | 否 | 1 | 每条 query 对 SILVA 保留命中数 |
| `--strand` | 否 | both | VSEARCH strand |
| `--centroid-policy` | 否 | original_seed | 正式 centroid 策略 |
| `--threads` | 否 | auto | 线程 |
| `--dry-run` | 否 | False | 只打印命令 |

### `--rank-thresholds`

默认：

```bash
--rank-thresholds default
```

等价于：

```text
species:0.987,genus:0.945,family:0.865,order:0.820,class:0.785,phylum:0.750
```

也可以自定义：

```bash
--rank-thresholds species:0.987,genus:0.945,family:0.865
```

或者只跑部分 rank：

```bash
--rank-thresholds species,genus
```

### 逻辑

对每个 rank：

```text
输入序列 vs SILVA
  ├── identity >= rank threshold
  │     → silva_existing
  │     → taxon_id = silva|rank|taxonomy_path
  │
  └── identity < rank threshold
        → 进入 unmatched
        → unmatched extension sequences 内部按该阈值聚类
        → 生成 novel|rank|source_label_rank_N
```

SILVA 不参与重新聚类。

### 输出目录

```text
hifimeta_inserted/
  sequence_rank_assignment.tsv
  rank_taxa_summary.tsv

  rank_uc/
    species.silva.uc
    species.novel.uc
    genus.silva.uc
    genus.novel.uc
    ...

  rank_centroids_core/
    species.centroids.core.fa
    genus.centroids.core.fa
    ...

  rank_centroids_original/
    species.centroids.original.fa
    genus.centroids.original.fa
    ...

  rank_hits/
    species_vs_silva.tsv
    species_vs_silva.uc
    ...

  rank_unmatched/
    species.unmatched.fa
    genus.unmatched.fa
    ...
```

---

## 12. `overlap-backbone`

基于 `sequence_rank_assignment.tsv` 统计不同 source 在不同 rank 的 shared / unique taxa。

### 命令

```bash
autotax2 overlap-backbone \
  --assignments \
    ref2_inserted/sequence_rank_assignment.tsv \
    ref3_inserted/sequence_rank_assignment.tsv \
    ref4_inserted/sequence_rank_assignment.tsv \
    hifimeta_inserted/sequence_rank_assignment.tsv \
  --labels ref2,ref3,ref4,hifimeta \
  --out ref_overlap
```

### 参数

| 参数 | 必需 | 说明 |
|---|---:|---|
| `--assignments` | 是 | 多个 `sequence_rank_assignment.tsv` |
| `--labels` | 是 | 与 assignments 顺序对应的逗号分隔标签 |
| `--out` | 是 | 输出目录 |

### 输出

```text
ref_overlap/
  taxon_presence_by_source.tsv
  taxon_count_by_source.tsv
  source_pairwise_overlap_by_rank.tsv
  source_unique_taxa_by_rank.tsv
```

---

## 13. 传统 VSEARCH 辅助命令

这些命令保留，用于兼容或手动流程。

### 13.1 `derep`

```bash
autotax2 derep \
  --input input.fa \
  --out work \
  --sort \
  --threads auto
```

调用：

```text
vsearch --derep_fulllength
vsearch --sortbysize
```

### 13.2 `cluster`

```bash
autotax2 cluster \
  --input work/derep_sorted.fa \
  --out clusters \
  --ids 0.99,0.97 \
  --method cluster_size \
  --threads auto
```

### 13.3 `classify`

```bash
autotax2 classify \
  --input clusters/otu099_centroids.fa \
  --ref-manifest refdatabases/autotax2_ref_manifest.tsv \
  --out classify \
  --threads auto
```

输出：

```text
sintax.tsv
silva_hits.tsv
typestrain_hits.tsv
taxonomy_summary.tsv
```

### 13.4 `assign`

新增序列分配到已有 centroids 或生成新簇：

```bash
autotax2 assign \
  --new new.fa \
  --old-centroids centroids_97.fa \
  --id 0.97 \
  --out assign_97 \
  --threads auto
```

### 13.5 `provenance`

基于 UC 文件和 source map 统计来源 overlap：

```bash
autotax2 provenance \
  --source-map source_map.tsv \
  --derep-uc work/derep.uc \
  --level-uc clusters/otu099.uc clusters/otu097.uc \
  --level-labels 0.99,0.97 \
  --out provenance
```

### 13.6 `summarize`

汇总分类证据：

```bash
autotax2 summarize \
  --sintax classify/sintax.tsv \
  --silva-hits classify/silva_hits.tsv \
  --typestrain-hits classify/typestrain_hits.tsv \
  --out classify/taxonomy_summary.tsv
```

---

## 14. 输出文件详解

### 14.1 `sequence_rank_assignment.tsv`

由 `insert-backbone` 生成，是主输出。

| 列 | 说明 |
|---|---|
| `sequence_id` | analysis sequence ID |
| `analysis_sequence_id` | analysis sequence ID |
| `original_sequence_id` | 原始序列 ID |
| `source_label` | 来源标签 |
| `rank` | species/genus/family/order/class/phylum |
| `rank_label` | species-like 等 |
| `identity_threshold` | 当前 rank 阈值 |
| `assignment_type` | `silva_existing`, `new_cluster_seed`, `new_cluster_member` |
| `taxon_id` | SILVA taxon 或 novel taxon |
| `parent_taxon_id` | parent SILVA taxon |
| `cluster_id` | 与 taxon_id 一致 |
| `silva_target` | 命中的 SILVA target ID |
| `silva_identity` | 对 SILVA identity |
| `silva_qcov` | query coverage |
| `silva_tcov` | target coverage |
| `is_novel` | yes/no |
| `has_intron` | yes/no/unknown |

### 14.2 `rank_taxa_summary.tsv`

每个 rank 下每个 taxon 的序列数。

| 列 | 说明 |
|---|---|
| `rank` | rank |
| `rank_label` | rank-like label |
| `identity_threshold` | 阈值 |
| `taxon_id` | taxon ID |
| `source_label` | 来源 |
| `n_sequences` | 序列数 |
| `is_novel` | 是否 novel |
| `assignment_types_json` | assignment type 组成 |

### 14.3 `rank_uc/*.uc`

每个 rank 有两个 UC：

| 文件 | 说明 |
|---|---|
| `species.silva.uc` | 达到 species 阈值并归入 SILVA 的 VSEARCH UC |
| `species.novel.uc` | 未归入 SILVA 的 extension sequences 内部聚类 UC |
| `genus.silva.uc` | genus 阈值 SILVA assignment UC |
| `genus.novel.uc` | genus 阈值 novel clustering UC |

### 14.4 `rank_centroids_core/*.fa`

计算用 centroids。

如果输入是 intron-free analysis FASTA，core centroid 可能是 intron-free 序列。

### 14.5 `rank_centroids_original/*.fa`

正式输出 centroids，使用 original sequence。

这满足：

```text
intron-free sequence 参与聚类
但 intron-free sequence 不作为正式 centroid
```

### 14.6 overlap 输出

#### `taxon_presence_by_source.tsv`

| 列 | 说明 |
|---|---|
| `rank` | rank |
| `taxon_id` | taxon |
| `sources` | 出现在哪些 source |
| `n_sources` | source 数 |
| `is_source_specific` | 是否只在一个 source 出现 |
| `source_specific_to` | 特异 source |
| `ref2/ref3/ref4/hifimeta...` | yes/no |

#### `taxon_count_by_source.tsv`

| 列 | 说明 |
|---|---|
| `rank` | rank |
| `taxon_id` | taxon |
| `source_label` | source |
| `n_sequences` | 该 source 中序列数 |
| `present` | yes/no |

#### `source_pairwise_overlap_by_rank.tsv`

| 列 | 说明 |
|---|---|
| `rank` | rank |
| `source_a` | source A |
| `source_b` | source B |
| `shared_taxon_count` | 共享 taxon 数 |
| `source_a_taxon_count` | A 的 taxon 数 |
| `source_b_taxon_count` | B 的 taxon 数 |
| `union_taxon_count` | A 或 B 的 taxon 总数 |
| `jaccard_taxon_overlap` | shared / union |

#### `source_unique_taxa_by_rank.tsv`

| 列 | 说明 |
|---|---|
| `rank` | rank |
| `taxon_id` | taxon |
| `source_label` | 唯一来源 |

---

## 15. 下游分析数据怎么用

AutoTax2 不做富集分析，但会输出适合下游分析的干净表。

### 15.1 每个 rank 有多少类群

用：

```text
rank_taxa_summary.tsv
```

按 `rank` 统计 `taxon_id` 数量。

### 15.2 每个类群有多少序列

用：

```text
rank_taxa_summary.tsv
```

字段：

```text
rank
taxon_id
n_sequences
is_novel
```

### 15.3 每条序列属于哪个 taxon

用：

```text
sequence_rank_assignment.tsv
```

可以转换成：

```text
sequence_id × rank → taxon_id
```

### 15.4 已知 vs novel

用：

```text
sequence_rank_assignment.tsv
rank_taxa_summary.tsv
```

字段：

```text
is_novel
assignment_type
```

### 15.5 ref2/ref3/ref4/hifimeta overlap

用：

```text
ref_overlap/taxon_presence_by_source.tsv
ref_overlap/source_pairwise_overlap_by_rank.tsv
ref_overlap/source_unique_taxa_by_rank.tsv
```

可以回答：

```text
hifimeta 中哪些 genus-like taxa 被 SILVA 覆盖？
哪些 hifimeta taxa 也出现在 ref2/ref3/ref4？
ref2/ref3/ref4 对 SILVA backbone 有多少补充？
哪些 taxa 是 hifimeta 特有？
```

### 15.6 intron 相关分析

用：

```text
intron_summary.tsv
intron_sequences.fa
intron_regions.bed
```

可以分析：

```text
哪些序列含 intron？
intron 出现在 query 的哪个位置？
intron-bearing 序列归入哪些 taxon？
intron-bearing 序列是否集中在某些 source 或类群？
```

把 `intron_summary.tsv` 与 `sequence_rank_assignment.tsv` 按 `sequence_id/original_sequence_id` join 即可。

### 15.7 引物评估 / 扩增子注释提升

当前 AutoTax2 不内置 primer evaluation。推荐使用以下输出作为下游输入：

```text
rank_centroids_original/*.fa
rank_centroids_core/*.fa
rank_uc/*.uc
sequence_rank_assignment.tsv
rank_taxa_summary.tsv
```

这些足够支持你后续独立评估：

```text
不同引物覆盖哪些 rank-like taxa？
某些扩增子区域是否合并了多个 full-length taxa？
加入 ref2/ref3/ref4/hifimeta 后，对扩增子注释提升多少？
```

---

## 16. 针对 HiFiMeta 500 万 16S 的建议

### 16.1 强烈建议先生成 UDB

```bash
autotax2 prepare-silva ... --make-udb
```

500 万 query 对 SILVA FASTA 直接搜索会慢，UDB 更适合重复使用。

### 16.2 分批运行

对于超大 FASTA，建议按 chunk 拆分：

```text
hifimeta.part001.fa
hifimeta.part002.fa
...
```

每个 chunk 独立：

```bash
autotax2 detect-intron ...
autotax2 insert-backbone ...
```

最后合并 `sequence_rank_assignment.tsv` 后做 overlap 或统计。

### 16.3 线程

使用：

```bash
--threads auto
```

或明确指定：

```bash
--threads 64
```

AutoTax2 会 capped 到检测到的 CPU 数，避免过度超售。

### 16.4 intron 检测建议

不一定所有序列都需要 intron 检测。如果计算量太大，可以先筛选候选：

```text
长度明显偏长
初步 SILVA identity 低于预期
但局部片段命中很好
```

然后只对候选 FASTA 运行 `detect-intron`。

当前版本提供全量 detect-intron 能力，但大规模数据建议先做候选筛选以节省时间。

---

## 17. 参数调优建议

### 17.1 `--rank-thresholds`

默认：

```bash
--rank-thresholds default
```

自定义：

```bash
--rank-thresholds species:0.987,genus:0.945,family:0.865
```

### 17.2 `--db-format`

| 值 | 说明 |
|---|---|
| `auto` | 优先 UDB，缺失时用 FASTA |
| `udb` | 尝试使用 UDB |
| `fasta` | 强制使用 FASTA |

推荐：

```bash
--db-format auto
```

### 17.3 `detect-intron --rescue-id`

默认：

```bash
--rescue-id 0.987
```

如果你希望更宽松：

```bash
--rescue-id 0.98
```

但更宽松会增加误判风险。

### 17.4 `detect-intron --min-intron-len`

默认：

```bash
--min-intron-len 50
```

如果只关心较大 intron：

```bash
--min-intron-len 100
```

### 17.5 `detect-intron --min-flank-len`

默认：

```bash
--min-flank-len 150
```

这个参数用于避免很短 flank 导致的误判。

### 17.6 `--maxaccepts`

`insert-backbone` 默认：

```bash
--maxaccepts 1
```

因为 backbone assignment 通常只需要 top hit。若想人工检查多命中，可在额外分析中提高。

---

## 18. 已知限制

### 18.1 SILVA taxonomy parsing 依赖 header 格式

`prepare-silva` 假设 header 格式是：

```text
>ID taxonomy;taxonomy;taxonomy;
```

如果使用非 SILVA 格式，需要先转换 header。

### 18.2 `detect-intron` 不是通用嵌合体检测

它识别的是：

```text
query 相对于同一个 subject 的长插入
```

不等价于 chimera detection。复杂拼接、多个 subject 混合、组装错误需要单独判断。

### 18.3 novel taxa 命名是 operational naming

`novel|genus|hifimeta_genus_1` 表示在当前阈值和当前流程下的 novel genus-like cluster，不等价于正式分类学命名。

### 18.4 当前没有内置丰度表

AutoTax2 当前主要处理 FASTA 和 UC，不内置样本丰度矩阵生成。你可以用 `sequence_rank_assignment.tsv` 自行构建。

### 18.5 当前没有内置 primer evaluation

目前只输出后续 primer evaluation 需要的 UC、centroids 和 assignment 表。

---

## 19. 开发路线

### v0.5 已有

```text
SILVA 本地准备
UDB 生成
rank threshold 体系
SILVA backbone insertion
backbone-based overlap
detect-intron 子命令
独立 scripts/detect_intron.py
core/original centroid 分离
Rich CLI
依赖检查
线程 auto/capping
```

### 建议 v0.6

```text
大数据 chunk runner
候选 intron 预筛选
merge assignment tables
rank-level abundance matrix
更严格的 parent-aware novel clustering
```

### 建议 v0.7

```text
primer evaluation 外挂模块
HTML summary report
可选 Snakemake / Nextflow workflow
```

---

## 许可证

MIT License。
