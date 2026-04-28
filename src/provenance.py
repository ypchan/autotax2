from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


def seq_id_from_header(header: str) -> str:
    return header.split()[0].split(";")[0]


def read_source_map(path: str | Path) -> Dict[str, str]:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(f"source map not found or empty: {p}")

    mapping: Dict[str, str] = {}
    with open(p, encoding="utf-8") as f:
        first = f.readline().rstrip("\n")
        if not first:
            return mapping

        parts = first.split("\t")
        lower = [x.lower() for x in parts]
        has_header = "sequence_id" in lower and "source" in lower

        if has_header:
            seq_i = lower.index("sequence_id")
            src_i = lower.index("source")
        else:
            seq_i = 0
            src_i = 1 if len(parts) > 1 else 0
            if len(parts) >= 2:
                mapping[parts[seq_i]] = parts[src_i] or "unknown"

        for line in f:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(seq_i, src_i):
                continue
            seq = row[seq_i].strip()
            src = row[src_i].strip()
            if seq:
                mapping[seq] = src if src else "unknown"
    return mapping


def parse_uc_query_to_cluster(uc_path: str | Path) -> Dict[str, str]:
    p = Path(uc_path)
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(f"UC file not found or empty: {p}")

    mapping: Dict[str, str] = {}
    cluster_num_to_seed: Dict[str, str] = {}

    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            rec_type = parts[0]
            cluster_num = parts[1]
            query = seq_id_from_header(parts[8])
            if rec_type == "S":
                cluster_num_to_seed[cluster_num] = query
                mapping[query] = query
            elif rec_type == "H":
                target = seq_id_from_header(parts[9]) if len(parts) > 9 else cluster_num_to_seed.get(cluster_num, "")
                mapping[query] = target
    return mapping


def compose_original_to_level_cluster(derep_uc: str | Path, level_uc: str | Path) -> Dict[str, str]:
    original_to_derep = parse_uc_query_to_cluster(derep_uc)
    derep_to_cluster = parse_uc_query_to_cluster(level_uc)
    return {original: derep_to_cluster.get(derep, derep) for original, derep in original_to_derep.items()}


def summarize_cluster_sources(
    seq_to_cluster: Dict[str, str],
    source_map: Dict[str, str],
    identity_label: str,
) -> Tuple[List[Dict[str, str]], Dict[str, Set[str]], Dict[str, Counter]]:
    cluster_source_counts: Dict[str, Counter] = defaultdict(Counter)
    cluster_members: Dict[str, Set[str]] = defaultdict(set)

    for seq, cluster in seq_to_cluster.items():
        source = source_map.get(seq, "unknown")
        cluster_source_counts[cluster][source] += 1
        cluster_members[cluster].add(seq)

    rows: List[Dict[str, str]] = []
    for cluster in sorted(cluster_source_counts):
        counts = cluster_source_counts[cluster]
        total = sum(counts.values())
        sources = sorted(counts)
        dominant_source, dominant_count = counts.most_common(1)[0]
        is_source_specific = len(sources) == 1
        rows.append({
            "level": identity_label,
            "cluster_id": cluster,
            "total_sequences": str(total),
            "unique_sources": str(len(sources)),
            "sources": ",".join(sources),
            "source_counts_json": json.dumps(dict(counts), ensure_ascii=False, sort_keys=True),
            "dominant_source": dominant_source,
            "dominant_source_count": str(dominant_count),
            "dominant_source_fraction": f"{dominant_count / total:.6f}" if total else "0",
            "is_source_specific": "yes" if is_source_specific else "no",
            "source_specific_to": sources[0] if is_source_specific else "",
        })
    return rows, cluster_members, cluster_source_counts


def pairwise_source_overlap(cluster_source_counts: Dict[str, Counter], identity_label: str) -> List[Dict[str, str]]:
    all_sources: Set[str] = set()
    source_clusters: Dict[str, Set[str]] = defaultdict(set)

    for cluster, counts in cluster_source_counts.items():
        for source, n in counts.items():
            if n > 0:
                all_sources.add(source)
                source_clusters[source].add(cluster)

    rows: List[Dict[str, str]] = []
    for a, b in combinations(sorted(all_sources), 2):
        clusters_a = source_clusters[a]
        clusters_b = source_clusters[b]
        shared = clusters_a & clusters_b
        union = clusters_a | clusters_b
        shared_seq_count = 0
        for cluster in shared:
            counts = cluster_source_counts[cluster]
            shared_seq_count += counts.get(a, 0) + counts.get(b, 0)
        rows.append({
            "level": identity_label,
            "source_a": a,
            "source_b": b,
            "shared_cluster_count": str(len(shared)),
            "source_a_cluster_count": str(len(clusters_a)),
            "source_b_cluster_count": str(len(clusters_b)),
            "union_cluster_count": str(len(union)),
            "jaccard_cluster_overlap": f"{len(shared) / len(union):.6f}" if union else "0",
            "shared_sequence_count_in_shared_clusters": str(shared_seq_count),
        })
    return rows


def source_level_summary(cluster_source_counts: Dict[str, Counter], identity_label: str) -> List[Dict[str, str]]:
    source_cluster_count: Dict[str, int] = defaultdict(int)
    source_specific_cluster_count: Dict[str, int] = defaultdict(int)
    source_sequence_count: Dict[str, int] = defaultdict(int)
    shared_cluster_count: Dict[str, int] = defaultdict(int)

    for _cluster, counts in cluster_source_counts.items():
        sources = [s for s, n in counts.items() if n > 0]
        for source, n in counts.items():
            source_sequence_count[source] += n
            source_cluster_count[source] += 1
            if len(sources) == 1:
                source_specific_cluster_count[source] += 1
            else:
                shared_cluster_count[source] += 1

    rows: List[Dict[str, str]] = []
    for source in sorted(source_cluster_count):
        total_clusters = source_cluster_count[source]
        specific = source_specific_cluster_count[source]
        shared = shared_cluster_count[source]
        rows.append({
            "level": identity_label,
            "source": source,
            "sequence_count": str(source_sequence_count[source]),
            "cluster_count": str(total_clusters),
            "source_specific_cluster_count": str(specific),
            "shared_cluster_count": str(shared),
            "source_specific_fraction": f"{specific / total_clusters:.6f}" if total_clusters else "0",
        })
    return rows


def write_tsv(rows: List[Dict[str, str]], path: str | Path, fieldnames: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def provenance_from_uc_levels(
    source_map_path: str | Path,
    derep_uc: str | Path,
    level_ucs: List[str | Path],
    outdir: str | Path,
    level_labels: Optional[List[str]] = None,
) -> Dict[str, str]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    source_map = read_source_map(source_map_path)

    if level_labels is None:
        level_labels = [Path(x).stem for x in level_ucs]
    if len(level_labels) != len(level_ucs):
        raise ValueError("--level-labels must have the same number of entries as --level-uc")

    all_cluster_rows: List[Dict[str, str]] = []
    all_pairwise_rows: List[Dict[str, str]] = []
    all_source_rows: List[Dict[str, str]] = []
    all_membership_rows: List[Dict[str, str]] = []

    for label, uc in zip(level_labels, level_ucs):
        seq_to_cluster = compose_original_to_level_cluster(derep_uc, uc)
        cluster_rows, _members, cluster_source_counts = summarize_cluster_sources(seq_to_cluster, source_map, label)
        all_cluster_rows.extend(cluster_rows)
        all_pairwise_rows.extend(pairwise_source_overlap(cluster_source_counts, label))
        all_source_rows.extend(source_level_summary(cluster_source_counts, label))
        for seq, cluster in sorted(seq_to_cluster.items()):
            all_membership_rows.append({
                "level": label,
                "sequence_id": seq,
                "source": source_map.get(seq, "unknown"),
                "cluster_id": cluster,
            })

    cluster_summary = outdir / "cluster_source_summary.tsv"
    pairwise_overlap = outdir / "source_pairwise_overlap.tsv"
    source_summary = outdir / "source_level_summary.tsv"
    membership = outdir / "sequence_cluster_membership_by_level.tsv"

    write_tsv(all_cluster_rows, cluster_summary, [
        "level", "cluster_id", "total_sequences", "unique_sources", "sources",
        "source_counts_json", "dominant_source", "dominant_source_count",
        "dominant_source_fraction", "is_source_specific", "source_specific_to",
    ])
    write_tsv(all_pairwise_rows, pairwise_overlap, [
        "level", "source_a", "source_b", "shared_cluster_count",
        "source_a_cluster_count", "source_b_cluster_count", "union_cluster_count",
        "jaccard_cluster_overlap", "shared_sequence_count_in_shared_clusters",
    ])
    write_tsv(all_source_rows, source_summary, [
        "level", "source", "sequence_count", "cluster_count",
        "source_specific_cluster_count", "shared_cluster_count",
        "source_specific_fraction",
    ])
    write_tsv(all_membership_rows, membership, [
        "level", "sequence_id", "source", "cluster_id",
    ])

    return {
        "cluster_source_summary": str(cluster_summary),
        "source_pairwise_overlap": str(pairwise_overlap),
        "source_level_summary": str(source_summary),
        "sequence_cluster_membership_by_level": str(membership),
    }


def annotate_assignments_with_source(assignments_tsv: str | Path, source_map_path: str | Path, output_tsv: str | Path) -> None:
    source_map = read_source_map(source_map_path)
    with open(assignments_tsv, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    out_fields = ["query", "sequence_source", "cluster_id", "assignment_type", "identity", "qcov", "tcov", "assignment_source"]
    with open(output_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "query": row.get("query", ""),
                "sequence_source": source_map.get(row.get("query", ""), "unknown"),
                "cluster_id": row.get("cluster_id", ""),
                "assignment_type": row.get("assignment_type", ""),
                "identity": row.get("identity", ""),
                "qcov": row.get("qcov", ""),
                "tcov": row.get("tcov", ""),
                "assignment_source": row.get("source", ""),
            })
