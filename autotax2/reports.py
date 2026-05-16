"""Global summaries and dataset overlap reports."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
NOVELTY_RANKS = ("species", "genus", "family", "order", "class")
OVERLAP_RANKS = ("exact_sequence", "species", "genus", "family", "order", "class")
SOURCE_TYPES = ("named_silva", "unresolved_silva", "previous_custom", "current_dataset")


@dataclass(frozen=True)
class ReportSummary:
    """Summary returned by :func:`summarize_build`."""

    build: Path
    outdir: Path
    files_written: int


def summarize_counts(counts: dict[str, int]) -> str:
    """Render simple count summaries."""
    return "\n".join(f"{key}: {value}" for key, value in sorted(counts.items()))


def summarize_build(
    build: str | Path,
    outdir: str | Path | None = None,
    overwrite: bool = False,
) -> ReportSummary:
    """Write global summary and reporting TSVs for an autotax2 build."""
    build_dir = Path(build)
    report_dir = Path(outdir) if outdir is not None else build_dir / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    registry_dir = build_dir / "registry"

    taxon_rows = _read_tsv(registry_dir / "taxon_nodes.tsv")
    sequence_rows = _read_tsv(registry_dir / "sequence_registry.tsv")
    dataset_rows = _custom_dataset_rows(registry_dir)
    representative_rows = _active_rows(_read_tsv(registry_dir / "representative_registry.tsv"))
    taxon_by_id = {row.get("taxon_id", ""): row for row in taxon_rows if row.get("taxon_id")}

    report_payloads = {
        "global_summary.tsv": (
            [_global_summary_row(build_dir, taxon_rows, sequence_rows, dataset_rows, representative_rows)],
            [
                "build_dir",
                "autotax2_version",
                "registry_version_or_build_id",
                "silva_named_sequences",
                "silva_unresolved_sequences",
                "silva_unresolved_active_placeholders",
                "silva_unresolved_evidence_records",
                "silva_unresolved_evidence_rows",
                "custom_datasets",
                "custom_input_sequences",
                "custom_unique_sequences",
                "duplicate_sequences",
                "active_taxa_total",
                "active_species",
                "active_genera",
                "active_families",
                "active_orders",
                "active_classes",
                "active_representatives",
                "deprecated_taxa",
                "superseded_taxa",
            ],
        ),
        "dataset_delta_summary.tsv": (
            _dataset_delta_rows(build_dir, dataset_rows),
            [
                "dataset",
                "prefix",
                "add_order",
                "input_sequences",
                "normalized_sequences",
                "prepared_sequences",
                "oriented_sequences",
                "unique_md5_sequences",
                "duplicate_sequences",
                "sina_candidate_queries",
                "sina_candidate_targets",
                "sina_candidate_target_matches",
                "placement_evidence_rows",
                "assigned_named_silva",
                "assigned_unresolved_silva",
                "assigned_previous_custom",
                "new_current_dataset",
                "known_like",
                "new_species",
                "new_genus",
                "new_family",
                "new_order",
                "new_class",
                "ambiguous",
                "unplaced",
                "created_species",
                "created_genera",
                "created_families",
                "created_orders",
                "created_classes",
                "representatives_added",
            ],
        ),
        "dataset_overlap_matrix.tsv": (
            _dataset_overlap_rows(build_dir, dataset_rows, taxon_by_id, sequence_rows),
            [
                "query_dataset",
                "query_prefix",
                "compared_source",
                "compared_source_type",
                "rank",
                "overlap_count",
                "query_total",
                "overlap_fraction",
            ],
        ),
        "rank_novelty_summary.tsv": (
            _rank_novelty_rows(build_dir, dataset_rows, taxon_by_id),
            ["dataset", "prefix", "rank", "created_taxa", "assigned_existing_taxa", "ambiguous", "unplaced"],
        ),
        "dataset_rank_overlap_detail.tsv": (
            _dataset_rank_overlap_detail_rows(build_dir, dataset_rows, taxon_by_id),
            [
                "dataset",
                "prefix",
                "rank",
                "rank_taxon_id",
                "rank_taxon_name",
                "parent_taxon_id",
                "parent_taxon_name",
                "source_type",
                "relation",
                "sequence_count",
                "sequence_ids",
                "final_statuses",
                "evidence_decisions",
                "created_in_dataset",
            ],
        ),
        "dataset_rank_novelty_detail.tsv": (
            _dataset_rank_novelty_detail_rows(build_dir, dataset_rows, taxon_by_id),
            [
                "dataset",
                "prefix",
                "rank",
                "taxon_id",
                "name",
                "parent_taxon_id",
                "parent_taxon_name",
                "cluster_key",
                "status",
                "is_placeholder",
                "representative_seq_id",
                "supporting_sequence_count",
                "supporting_sequence_ids",
            ],
        ),
        "source_contribution.tsv": (
            _source_contribution_rows(taxon_rows, sequence_rows, representative_rows),
            [
                "source",
                "source_type",
                "sequences",
                "unique_sequences",
                "representatives",
                "active_taxa",
                "active_species",
                "active_genera",
                "active_families",
                "active_orders",
                "active_classes",
            ],
        ),
        "representative_summary.tsv": (
            _representative_summary_rows(build_dir, representative_rows, taxon_by_id, sequence_rows),
            [
                "taxon_id",
                "rank",
                "taxon_name",
                "representative_seq_id",
                "representative_source",
                "representative_dataset",
                "representative_reason",
                "is_type_strain",
                "protected",
                "sequence_length",
                "sequence_md5",
            ],
        ),
        "sequence_dedup_summary.tsv": (
            _sequence_dedup_rows(build_dir, sequence_rows, representative_rows),
            [
                "sequence_md5",
                "unique_seq_id",
                "representative_internal_seq_id",
                "first_seen_dataset",
                "occurrence_count",
                "datasets",
                "exported",
                "taxon_id",
            ],
        ),
    }

    files_written = 0
    for filename, (rows, fieldnames) in report_payloads.items():
        path = report_dir / filename
        _ensure_writable(path, overwrite)
        _write_tsv(rows, path, fieldnames)
        files_written += 1

    audit_path = report_dir / "dataset_increment_audit.md"
    _ensure_writable(audit_path, overwrite)
    audit_path.write_text(
        _dataset_increment_audit_markdown(build_dir, report_payloads),
        encoding="utf-8",
        newline="\n",
    )
    files_written += 1

    return ReportSummary(build=build_dir, outdir=report_dir, files_written=files_written)


def _global_summary_row(
    build_dir: Path,
    taxon_rows: list[dict[str, str]],
    sequence_rows: list[dict[str, str]],
    dataset_rows: list[dict[str, str]],
    representative_rows: list[dict[str, str]],
) -> dict[str, str]:
    active_taxa = [row for row in taxon_rows if _is_active(row)]
    status_counts = Counter(row.get("status", "active") for row in taxon_rows)
    active_rank_counts = Counter(row.get("rank", "") for row in active_taxa)
    dataset_dirs = [_dataset_dir(build_dir, row) for row in dataset_rows]
    membership_rows = [member for directory in dataset_dirs for member in _read_tsv(directory / "sequence_membership.tsv")]
    sequence_md5s = [row.get("sequence_md5", "") for row in membership_rows if row.get("sequence_md5")]
    global_md5s = [row.get("sequence_md5", "") for row in sequence_rows if row.get("sequence_md5")]
    duplicate_count = sum(count - 1 for count in Counter(global_md5s + sequence_md5s).values() if count > 1)
    duplicate_count += sum(1 for row in membership_rows if _is_true(row.get("is_duplicate_sequence", "false")))
    silva_evidence = _read_tsv(build_dir / "silva" / "silva_unresolved_evidence.tsv")

    return {
        "build_dir": str(build_dir),
        "autotax2_version": _autotax2_version(),
        "registry_version_or_build_id": build_dir.name or "build",
        "silva_named_sequences": str(sum(1 for row in sequence_rows if _is_true(row.get("is_silva_named", "false")))),
        "silva_unresolved_sequences": str(sum(1 for row in sequence_rows if _is_true(row.get("is_silva_unresolved", "false")))),
        "silva_unresolved_active_placeholders": str(
            sum(
                1
                for row in active_taxa
                if _is_true(row.get("is_silva_unresolved", "false"))
                or (row.get("source_prefix") == "SILVA" and _is_true(row.get("is_placeholder", "false")))
            )
        ),
        "silva_unresolved_evidence_records": str(len({row.get("seq_id", "") for row in silva_evidence if row.get("seq_id")})),
        "silva_unresolved_evidence_rows": str(len(silva_evidence)),
        "custom_datasets": str(len(dataset_rows)),
        "custom_input_sequences": str(sum(_int(_prepare_summary(_dataset_dir(build_dir, row)).get("input_sequences")) for row in dataset_rows)),
        "custom_unique_sequences": str(len(set(sequence_md5s))),
        "duplicate_sequences": str(duplicate_count),
        "active_taxa_total": str(len(active_taxa)),
        "active_species": str(active_rank_counts.get("species", 0)),
        "active_genera": str(active_rank_counts.get("genus", 0)),
        "active_families": str(active_rank_counts.get("family", 0)),
        "active_orders": str(active_rank_counts.get("order", 0)),
        "active_classes": str(active_rank_counts.get("class", 0)),
        "active_representatives": str(len(representative_rows)),
        "deprecated_taxa": str(status_counts.get("deprecated", 0)),
        "superseded_taxa": str(status_counts.get("superseded", 0)),
    }


def _dataset_delta_rows(build_dir: Path, dataset_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_row in dataset_rows:
        dataset_dir = _dataset_dir(build_dir, dataset_row)
        prepare = _prepare_summary(dataset_dir)
        placement = _first_row(dataset_dir / "placement_summary.tsv")
        cluster_summary = _first_row(dataset_dir / "cluster_search_summary.tsv")
        assignments = _read_tsv(dataset_dir / "assignments.tsv")
        members = _read_tsv(dataset_dir / "sequence_membership.tsv")
        created = _read_tsv(dataset_dir / "created_taxa.tsv")
        reps = _read_tsv(dataset_dir / "representative_updates.tsv")
        placement_evidence = _read_tsv(dataset_dir / "placement_evidence.tsv")
        assignment_status = Counter(row.get("final_status", "") for row in assignments)
        created_counts = Counter(row.get("rank", "") for row in created)
        rows.append(
            {
                "dataset": dataset_row.get("dataset_name", ""),
                "prefix": dataset_row.get("prefix", ""),
                "add_order": dataset_row.get("add_order", ""),
                "input_sequences": prepare.get("input_sequences", str(len(_read_tsv(dataset_dir / "sequence_id_map.tsv")))),
                "normalized_sequences": prepare.get("normalized_sequences", str(len(members))),
                "prepared_sequences": prepare.get("final_prepared_sequences", _count_fasta_records(dataset_dir / "prepared.ssu.fa")),
                "oriented_sequences": _count_fasta_records(dataset_dir / "sina.oriented.fa"),
                "unique_md5_sequences": str(len({row.get("sequence_md5", "") for row in members if row.get("sequence_md5")})),
                "duplicate_sequences": str(sum(1 for row in members if _is_true(row.get("is_duplicate_sequence", "false")))),
                "sina_candidate_queries": cluster_summary.get("sina_candidate_queries", "0"),
                "sina_candidate_targets": cluster_summary.get("sina_candidate_targets", "0"),
                "sina_candidate_target_matches": cluster_summary.get("sina_candidate_target_matches", "0"),
                "placement_evidence_rows": str(len(placement_evidence)),
                "assigned_named_silva": placement.get("assigned_named_silva", "0"),
                "assigned_unresolved_silva": placement.get("assigned_unresolved_silva", "0"),
                "assigned_previous_custom": placement.get("assigned_previous_custom", "0"),
                "new_current_dataset": placement.get("new_current_dataset", "0"),
                "known_like": str(assignment_status.get("known_like", _int(placement.get("known_like")))),
                "new_species": str(assignment_status.get("new_species", _int(placement.get("new_species")))),
                "new_genus": str(assignment_status.get("new_genus", _int(placement.get("new_genus")))),
                "new_family": str(assignment_status.get("new_family", _int(placement.get("new_family")))),
                "new_order": str(assignment_status.get("new_order", _int(placement.get("new_order")))),
                "new_class": str(assignment_status.get("new_class", _int(placement.get("new_class")))),
                "ambiguous": str(assignment_status.get("ambiguous", _int(placement.get("ambiguous")))),
                "unplaced": str(assignment_status.get("unplaced", _int(placement.get("unplaced")))),
                "created_species": str(created_counts.get("species", _int(placement.get("created_species")))),
                "created_genera": str(created_counts.get("genus", _int(placement.get("created_genera")))),
                "created_families": str(created_counts.get("family", _int(placement.get("created_families")))),
                "created_orders": str(created_counts.get("order", _int(placement.get("created_orders")))),
                "created_classes": str(created_counts.get("class", _int(placement.get("created_classes")))),
                "representatives_added": str(len(reps)),
            }
        )
    return rows


def _dataset_overlap_rows(
    build_dir: Path,
    dataset_rows: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
    sequence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    md5_sources: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for sequence in sequence_rows:
        digest = sequence.get("sequence_md5", "")
        if not digest:
            continue
        source_type = _sequence_source_type(sequence, taxon_by_id)
        source = sequence.get("dataset") or sequence.get("source") or source_type
        md5_sources[digest].add((source, source_type))

    for dataset_row in dataset_rows:
        dataset = dataset_row.get("dataset_name", "")
        prefix = dataset_row.get("prefix", "")
        dataset_dir = _dataset_dir(build_dir, dataset_row)
        members = _read_tsv(dataset_dir / "sequence_membership.tsv")
        assignments = _read_tsv(dataset_dir / "assignments.tsv")
        query_md5s = {row.get("sequence_md5", "") for row in members if row.get("sequence_md5")}
        query_total_md5 = len(query_md5s)
        exact_counts: Counter[tuple[str, str]] = Counter()
        for digest in query_md5s:
            for source, source_type in md5_sources.get(digest, set()):
                if source != dataset:
                    exact_counts[(source, source_type)] += 1
        for source_type in SOURCE_TYPES:
            source = source_type
            count = sum(value for (src, typ), value in exact_counts.items() if typ == source_type for source in [src])
            rows.append(_overlap_row(dataset, prefix, source, source_type, "exact_sequence", count, query_total_md5))

        for rank in OVERLAP_RANKS[1:]:
            query_total = sum(1 for row in assignments if row.get("assigned_taxon_id"))
            counts: Counter[str] = Counter()
            for assignment in assignments:
                taxon_id = assignment.get("assigned_taxon_id", "")
                lineage = _lineage_by_rank(taxon_id, taxon_by_id)
                rank_taxon = lineage.get(rank, "")
                if not rank_taxon:
                    continue
                counts[_taxon_source_type(rank_taxon, taxon_by_id, current_dataset=dataset)] += 1
            for source_type in SOURCE_TYPES:
                rows.append(_overlap_row(dataset, prefix, source_type, source_type, rank, counts.get(source_type, 0), query_total))
    return rows


def _rank_novelty_rows(
    build_dir: Path,
    dataset_rows: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_row in dataset_rows:
        dataset = dataset_row.get("dataset_name", "")
        dataset_dir = _dataset_dir(build_dir, dataset_row)
        created = Counter(row.get("rank", "") for row in _read_tsv(dataset_dir / "created_taxa.tsv"))
        assignments = _read_tsv(dataset_dir / "assignments.tsv")
        status_counts = Counter(row.get("final_status", "") for row in assignments)
        for rank in NOVELTY_RANKS:
            assigned_existing = 0
            for assignment in assignments:
                taxon_id = assignment.get("assigned_taxon_id", "")
                lineage = _lineage_by_rank(taxon_id, taxon_by_id)
                rank_taxon = lineage.get(rank, "")
                if rank_taxon and _taxon_source_type(rank_taxon, taxon_by_id, dataset) != "current_dataset":
                    assigned_existing += 1
            rows.append(
                {
                    "dataset": dataset,
                    "prefix": dataset_row.get("prefix", ""),
                    "rank": rank,
                    "created_taxa": str(created.get(rank, 0)),
                    "assigned_existing_taxa": str(assigned_existing),
                    "ambiguous": str(status_counts.get("ambiguous", 0)),
                    "unplaced": str(status_counts.get("unplaced", 0)),
                }
            )
    return rows


def _dataset_rank_overlap_detail_rows(
    build_dir: Path,
    dataset_rows: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_row in dataset_rows:
        dataset = dataset_row.get("dataset_name", "")
        prefix = dataset_row.get("prefix", "")
        dataset_dir = _dataset_dir(build_dir, dataset_row)
        assignments = _read_tsv(dataset_dir / "assignments.tsv")
        evidence = _placement_evidence_by_seq_rank(dataset_dir)
        groups: dict[tuple[str, str, str, str], dict[str, Any]] = {}

        for assignment in assignments:
            seq_id = assignment.get("internal_seq_id", "")
            final_status = assignment.get("final_status", "")
            lineage = _lineage_by_rank(assignment.get("assigned_taxon_id", ""), taxon_by_id)
            for rank in RANKS[1:]:
                rank_taxon_id = lineage.get(rank, "")
                source_type = _taxon_source_type(rank_taxon_id, taxon_by_id, dataset) if rank_taxon_id else ""
                relation = _rank_relation(rank_taxon_id, source_type, final_status)
                key = (rank, rank_taxon_id, source_type, relation)
                payload = groups.setdefault(
                    key,
                    {
                        "seq_ids": [],
                        "statuses": Counter(),
                        "evidence_decisions": Counter(),
                    },
                )
                if seq_id:
                    payload["seq_ids"].append(seq_id)
                if final_status:
                    payload["statuses"][final_status] += 1
                evidence_decision = evidence.get((seq_id, rank), {}).get("decision", "")
                if evidence_decision:
                    payload["evidence_decisions"][evidence_decision] += 1

        for (rank, rank_taxon_id, source_type, relation), payload in sorted(groups.items()):
            taxon = taxon_by_id.get(rank_taxon_id, {})
            parent_id = taxon.get("parent_taxon_id", "")
            parent = taxon_by_id.get(parent_id, {})
            rows.append(
                {
                    "dataset": dataset,
                    "prefix": prefix,
                    "rank": rank,
                    "rank_taxon_id": rank_taxon_id,
                    "rank_taxon_name": taxon.get("name", ""),
                    "parent_taxon_id": parent_id,
                    "parent_taxon_name": parent.get("name", ""),
                    "source_type": source_type,
                    "relation": relation,
                    "sequence_count": str(len(payload["seq_ids"])),
                    "sequence_ids": ",".join(sorted(payload["seq_ids"])),
                    "final_statuses": _counter_text(payload["statuses"]),
                    "evidence_decisions": _counter_text(payload["evidence_decisions"]),
                    "created_in_dataset": taxon.get("created_in_dataset", ""),
                }
            )
    return rows


def _dataset_rank_novelty_detail_rows(
    build_dir: Path,
    dataset_rows: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset_row in dataset_rows:
        dataset = dataset_row.get("dataset_name", "")
        prefix = dataset_row.get("prefix", "")
        dataset_dir = _dataset_dir(build_dir, dataset_row)
        assignments = _read_tsv(dataset_dir / "assignments.tsv")
        reps = _read_tsv(dataset_dir / "representative_updates.tsv")
        created = _read_tsv(dataset_dir / "created_taxa.tsv")
        support_by_taxon = _supporting_sequences_by_taxon(assignments, taxon_by_id)
        reps_by_taxon: dict[str, list[str]] = defaultdict(list)
        for rep in reps:
            taxon_id = rep.get("taxon_id", "")
            seq_id = rep.get("representative_seq_id", "")
            if taxon_id and seq_id:
                reps_by_taxon[taxon_id].append(seq_id)

        for created_row in created:
            taxon_id = created_row.get("taxon_id", "")
            taxon = {**taxon_by_id.get(taxon_id, {}), **created_row}
            parent_id = taxon.get("parent_taxon_id", "")
            parent = taxon_by_id.get(parent_id, {})
            supporting = sorted(support_by_taxon.get(taxon_id, set()))
            rows.append(
                {
                    "dataset": dataset,
                    "prefix": prefix,
                    "rank": taxon.get("rank", ""),
                    "taxon_id": taxon_id,
                    "name": taxon.get("name", ""),
                    "parent_taxon_id": parent_id,
                    "parent_taxon_name": parent.get("name", ""),
                    "cluster_key": taxon.get("cluster_key", ""),
                    "status": taxon.get("status", "active"),
                    "is_placeholder": _bool_text(_is_true(taxon.get("is_placeholder", "false")) or _looks_like_placeholder(taxon.get("name", ""))),
                    "representative_seq_id": ",".join(sorted(set(reps_by_taxon.get(taxon_id, [])))),
                    "supporting_sequence_count": str(len(supporting)),
                    "supporting_sequence_ids": ",".join(supporting),
                }
            )
    return rows


def _dataset_increment_audit_markdown(
    build_dir: Path,
    report_payloads: dict[str, tuple[list[dict[str, Any]], list[str]]],
) -> str:
    global_summary = _payload_rows(report_payloads, "global_summary.tsv")
    delta_rows = _payload_rows(report_payloads, "dataset_delta_summary.tsv")
    matrix_rows = _payload_rows(report_payloads, "dataset_overlap_matrix.tsv")
    overlap_rows = _payload_rows(report_payloads, "dataset_rank_overlap_detail.tsv")
    novelty_rows = _payload_rows(report_payloads, "dataset_rank_novelty_detail.tsv")
    global_row = global_summary[0] if global_summary else {}
    lines = [
        "# autotax2 Dataset Increment Audit",
        "",
        f"Build: `{build_dir}`",
        "",
        "## Global Summary",
        "",
        *_markdown_table(
            ["Metric", "Value"],
            [
                ["Custom datasets", global_row.get("custom_datasets", "0")],
                ["SILVA named sequences", global_row.get("silva_named_sequences", "0")],
                ["SILVA unresolved sequences", global_row.get("silva_unresolved_sequences", "0")],
                ["SILVA resolve evidence rows", global_row.get("silva_unresolved_evidence_rows", "0")],
                ["Active species", global_row.get("active_species", "0")],
                ["Active representatives", global_row.get("active_representatives", "0")],
                ["Duplicate sequences", global_row.get("duplicate_sequences", "0")],
            ],
        ),
    ]

    for delta in sorted(delta_rows, key=lambda row: (_int(row.get("add_order")), row.get("dataset", ""))):
        dataset = delta.get("dataset", "")
        prefix = delta.get("prefix", "")
        dataset_matrix = [row for row in matrix_rows if row.get("query_dataset") == dataset]
        dataset_overlap = [row for row in overlap_rows if row.get("dataset") == dataset]
        dataset_novelty = [row for row in novelty_rows if row.get("dataset") == dataset]
        exact_rows = [row for row in dataset_matrix if row.get("rank") == "exact_sequence" and _int(row.get("overlap_count")) > 0]

        lines.extend(
            [
                "",
                f"## Dataset {dataset} ({prefix})",
                "",
                "### Processing",
                "",
                *_markdown_table(
                    ["Metric", "Value"],
                    [
                        ["Add order", delta.get("add_order", "")],
                        ["Input sequences", delta.get("input_sequences", "0")],
                        ["Prepared sequences", delta.get("prepared_sequences", "0")],
                        ["Oriented sequences", delta.get("oriented_sequences", "0")],
                        ["Unique MD5 sequences", delta.get("unique_md5_sequences", "0")],
                        ["Duplicate sequences", delta.get("duplicate_sequences", "0")],
                        ["SINA candidate queries", delta.get("sina_candidate_queries", "0")],
                        ["SINA candidate targets", delta.get("sina_candidate_targets", "0")],
                        ["SINA candidate target matches", delta.get("sina_candidate_target_matches", "0")],
                        ["Placement evidence rows", delta.get("placement_evidence_rows", "0")],
                    ],
                ),
                "",
                "### Placement Outcomes",
                "",
                *_markdown_table(
                    ["Outcome", "Count"],
                    [
                        ["Known-like", delta.get("known_like", "0")],
                        ["New species", delta.get("new_species", "0")],
                        ["New genus", delta.get("new_genus", "0")],
                        ["New family", delta.get("new_family", "0")],
                        ["New order", delta.get("new_order", "0")],
                        ["New class", delta.get("new_class", "0")],
                        ["Ambiguous", delta.get("ambiguous", "0")],
                        ["Unplaced", delta.get("unplaced", "0")],
                    ],
                ),
                "",
                "### Assigned Source Summary",
                "",
                *_markdown_table(
                    ["Source", "Assigned"],
                    [
                        ["Named SILVA", delta.get("assigned_named_silva", "0")],
                        ["Unresolved SILVA", delta.get("assigned_unresolved_silva", "0")],
                        ["Previous custom", delta.get("assigned_previous_custom", "0")],
                        ["Current dataset", delta.get("new_current_dataset", "0")],
                    ],
                ),
                "",
                "### Rank Overlap Detail",
                "",
                *_rank_overlap_markdown_table(dataset_overlap),
                "",
                "### Exact Sequence Overlap",
                "",
                *_exact_overlap_markdown_table(exact_rows),
                "",
                "### Created Taxa",
                "",
                *_novelty_markdown_table(dataset_novelty),
                "",
                "### Trace Files",
                "",
                "- `dataset_delta_summary.tsv`",
                "- `dataset_overlap_matrix.tsv`",
                "- `dataset_rank_overlap_detail.tsv`",
                "- `dataset_rank_novelty_detail.tsv`",
                "- `placement_evidence.tsv`",
                "- `sina_candidate_diagnostics.tsv` when SINA candidate search was used",
            ]
        )
    lines.append("")
    return "\n".join(lines)


def _source_contribution_rows(
    taxon_rows: list[dict[str, str]],
    sequence_rows: list[dict[str, str]],
    representative_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], dict[str, Any]] = {}
    taxon_by_id = {row.get("taxon_id", ""): row for row in taxon_rows if row.get("taxon_id")}

    def group(source: str, source_type: str) -> dict[str, Any]:
        return groups.setdefault(
            (source, source_type),
            {"seq": [], "reps": 0, "taxa": []},
        )

    for row in sequence_rows:
        source_type = _sequence_source_type(row, taxon_by_id)
        source = row.get("dataset") or row.get("source") or source_type
        group(source, source_type)["seq"].append(row)
    for row in representative_rows:
        source_type = row.get("source_category") or _taxon_source_type(row.get("taxon_id", ""), taxon_by_id, "")
        source = row.get("dataset") or source_type
        group(source, source_type)["reps"] += 1
    for row in taxon_rows:
        if not _is_active(row):
            continue
        source_type = _taxon_source_type(row.get("taxon_id", ""), taxon_by_id, "")
        source = row.get("created_in_dataset") or row.get("source") or source_type
        group(source, source_type)["taxa"].append(row)

    rows = []
    for (source, source_type), payload in sorted(groups.items()):
        rank_counts = Counter(row.get("rank", "") for row in payload["taxa"])
        md5s = {row.get("sequence_md5", "") for row in payload["seq"] if row.get("sequence_md5")}
        rows.append(
            {
                "source": source,
                "source_type": source_type,
                "sequences": str(len(payload["seq"])),
                "unique_sequences": str(len(md5s)),
                "representatives": str(payload["reps"]),
                "active_taxa": str(len(payload["taxa"])),
                "active_species": str(rank_counts.get("species", 0)),
                "active_genera": str(rank_counts.get("genus", 0)),
                "active_families": str(rank_counts.get("family", 0)),
                "active_orders": str(rank_counts.get("order", 0)),
                "active_classes": str(rank_counts.get("class", 0)),
            }
        )
    return rows


def _representative_summary_rows(
    build_dir: Path,
    representative_rows: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
    sequence_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    seq_meta = {
        row.get("seq_id") or row.get("internal_seq_id", ""): row
        for row in sequence_rows
        if row.get("seq_id") or row.get("internal_seq_id")
    }
    lengths = _sequence_lengths(build_dir)
    rows = []
    for rep in representative_rows:
        seq_id = rep.get("representative_seq_id", "")
        taxon = taxon_by_id.get(rep.get("taxon_id", ""), {})
        meta = seq_meta.get(seq_id, {})
        rows.append(
            {
                "taxon_id": rep.get("taxon_id", ""),
                "rank": taxon.get("rank", ""),
                "taxon_name": taxon.get("name", ""),
                "representative_seq_id": seq_id,
                "representative_source": rep.get("source_category") or _taxon_source_type(rep.get("taxon_id", ""), taxon_by_id, ""),
                "representative_dataset": rep.get("dataset", ""),
                "representative_reason": rep.get("reason", ""),
                "is_type_strain": meta.get("is_type_strain", ""),
                "protected": taxon.get("protected", rep.get("protected", "")),
                "sequence_length": str(lengths.get(seq_id, _int(meta.get("sequence_length")))),
                "sequence_md5": meta.get("sequence_md5", ""),
            }
        )
    return rows


def _sequence_dedup_rows(
    build_dir: Path,
    sequence_rows: list[dict[str, str]],
    representative_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    all_rows = list(sequence_rows)
    for dataset_row in _custom_dataset_rows(build_dir / "registry"):
        all_rows.extend(_read_tsv(_dataset_dir(build_dir, dataset_row) / "sequence_membership.tsv"))
    representative_ids = {row.get("representative_seq_id", "") for row in representative_rows}
    exported_md5s = {
        row.get("sequence_md5", "")
        for row in sequence_rows
        if (row.get("seq_id") or row.get("internal_seq_id")) in representative_ids
    }
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in all_rows:
        digest = row.get("sequence_md5", "")
        if digest:
            groups[digest].append(row)
    rows = []
    for digest, members in sorted(groups.items()):
        first = members[0]
        datasets = sorted({row.get("dataset") or row.get("source") or "" for row in members if row.get("dataset") or row.get("source")})
        rows.append(
            {
                "sequence_md5": digest,
                "unique_seq_id": first.get("unique_seq_id", ""),
                "representative_internal_seq_id": first.get("representative_internal_seq_id") or first.get("seq_id") or first.get("internal_seq_id", ""),
                "first_seen_dataset": first.get("first_seen_dataset") or (datasets[0] if datasets else ""),
                "occurrence_count": str(len(members)),
                "datasets": ",".join(datasets),
                "exported": _bool_text(digest in exported_md5s),
                "taxon_id": first.get("taxon_id") or first.get("assigned_taxon_id", ""),
            }
        )
    return rows


def _overlap_row(
    dataset: str,
    prefix: str,
    source: str,
    source_type: str,
    rank: str,
    count: int,
    total: int,
) -> dict[str, str]:
    return {
        "query_dataset": dataset,
        "query_prefix": prefix,
        "compared_source": source,
        "compared_source_type": source_type,
        "rank": rank,
        "overlap_count": str(count),
        "query_total": str(total),
        "overlap_fraction": f"{(count / total if total else 0):.6f}",
    }


def _lineage_by_rank(taxon_id: str, taxon_by_id: dict[str, dict[str, str]]) -> dict[str, str]:
    lineage: dict[str, str] = {}
    seen: set[str] = set()
    current = taxon_id
    while current and current not in seen:
        seen.add(current)
        row = taxon_by_id.get(current, {})
        if not row:
            break
        if row.get("rank") in RANKS:
            lineage[row["rank"]] = current
        current = row.get("parent_taxon_id", "")
    return lineage


def _taxon_source_type(taxon_id: str, taxon_by_id: dict[str, dict[str, str]], current_dataset: str) -> str:
    row = taxon_by_id.get(taxon_id, {})
    if row.get("created_in_dataset") == current_dataset and current_dataset:
        return "current_dataset"
    if row.get("created_in_dataset"):
        return "previous_custom"
    if _is_true(row.get("is_silva_unresolved", "false")) or row.get("source_prefix") == "SILVA":
        return "unresolved_silva"
    if row.get("source", "").startswith("SILVA") or _is_true(row.get("is_silva_named", "false")) or _is_true(row.get("protected", "false")):
        return "named_silva"
    return "current_dataset" if row.get("source") == "custom_dataset" else "named_silva"


def _rank_relation(rank_taxon_id: str, source_type: str, final_status: str) -> str:
    if not rank_taxon_id:
        if final_status in {"ambiguous", "unplaced"}:
            return final_status
        return "missing_rank_call"
    if source_type == "current_dataset":
        return "new_current_dataset"
    if source_type == "previous_custom":
        return "overlap_previous_custom"
    if source_type == "unresolved_silva":
        return "overlap_unresolved_silva"
    if source_type == "named_silva":
        return "overlap_named_silva"
    return "overlap_existing"


def _placement_evidence_by_seq_rank(dataset_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    evidence: dict[tuple[str, str], dict[str, str]] = {}
    for row in _read_tsv(dataset_dir / "placement_evidence.tsv"):
        seq_id = row.get("seq_id", "")
        rank = row.get("rank", "")
        if seq_id and rank:
            evidence[(seq_id, rank)] = row
    return evidence


def _supporting_sequences_by_taxon(
    assignments: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
) -> dict[str, set[str]]:
    support: dict[str, set[str]] = defaultdict(set)
    for assignment in assignments:
        seq_id = assignment.get("internal_seq_id", "")
        if not seq_id:
            continue
        lineage = _lineage_by_rank(assignment.get("assigned_taxon_id", ""), taxon_by_id)
        for taxon_id in lineage.values():
            if taxon_id:
                support[taxon_id].add(seq_id)
    return support


def _counter_text(counter: Counter[str]) -> str:
    return ",".join(f"{key}:{value}" for key, value in sorted(counter.items()) if key)


def _looks_like_placeholder(name: str) -> bool:
    return len(name) >= 5 and name[:3] in {"p__", "c__", "o__", "f__", "g__", "s__"} and any(char.isdigit() for char in name)


def _payload_rows(
    report_payloads: dict[str, tuple[list[dict[str, Any]], list[str]]],
    filename: str,
) -> list[dict[str, Any]]:
    return report_payloads.get(filename, ([], []))[0]


def _rank_overlap_markdown_table(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return ["No rank overlap rows were available."]
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        rank = str(row.get("rank", ""))
        relation = str(row.get("relation", ""))
        if rank:
            counts[rank][relation] += _int(str(row.get("sequence_count", "0")))
    table_rows = []
    relation_columns = [
        "overlap_named_silva",
        "overlap_unresolved_silva",
        "overlap_previous_custom",
        "new_current_dataset",
        "ambiguous",
        "unplaced",
    ]
    for rank in RANKS[1:]:
        table_rows.append([rank, *[str(counts[rank].get(relation, 0)) for relation in relation_columns]])
    return _markdown_table(
        ["Rank", "Named SILVA", "Unresolved SILVA", "Previous custom", "New current", "Ambiguous", "Unplaced"],
        table_rows,
    )


def _exact_overlap_markdown_table(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return ["No exact sequence overlap with existing sources was reported."]
    return _markdown_table(
        ["Compared source", "Source type", "Overlap", "Query total", "Fraction"],
        [
            [
                str(row.get("compared_source", "")),
                str(row.get("compared_source_type", "")),
                str(row.get("overlap_count", "0")),
                str(row.get("query_total", "0")),
                str(row.get("overlap_fraction", "0.000000")),
            ]
            for row in rows
        ],
    )


def _novelty_markdown_table(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return ["No new taxa were created for this dataset."]
    return _markdown_table(
        ["Rank", "Taxon", "Parent", "Representative", "Support"],
        [
            [
                str(row.get("rank", "")),
                str(row.get("name") or row.get("taxon_id", "")),
                str(row.get("parent_taxon_name") or row.get("parent_taxon_id", "")),
                str(row.get("representative_seq_id", "")),
                str(row.get("supporting_sequence_ids", "")),
            ]
            for row in sorted(rows, key=lambda row: (RANKS.index(str(row.get("rank", "domain"))) if str(row.get("rank", "")) in RANKS else 99, str(row.get("taxon_id", ""))))
        ],
    )


def _markdown_table(headers: list[str], rows: list[list[str]]) -> list[str]:
    return [
        "| " + " | ".join(_md_cell(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *["| " + " | ".join(_md_cell(cell) for cell in row) + " |" for row in rows],
    ]


def _md_cell(value: str) -> str:
    clean = str(value).replace("\n", " ").replace("\r", " ")
    return clean.replace("|", "\\|") if clean else "0"


def _sequence_source_type(row: dict[str, str], taxon_by_id: dict[str, dict[str, str]]) -> str:
    if _is_true(row.get("is_silva_unresolved", "false")):
        return "unresolved_silva"
    if _is_true(row.get("is_silva_named", "false")) or row.get("source", "").startswith("SILVA"):
        return "named_silva"
    taxon_id = row.get("taxon_id") or row.get("assigned_taxon_id", "")
    if taxon_id:
        return _taxon_source_type(taxon_id, taxon_by_id, row.get("dataset", ""))
    return "previous_custom" if row.get("dataset") else "named_silva"


def _custom_dataset_rows(registry_dir: Path) -> list[dict[str, str]]:
    return [
        row
        for row in _read_tsv(registry_dir / "dataset_registry.tsv")
        if row.get("dataset_name") and row.get("prefix") and row.get("prefix") != "SILVA"
    ]


def _dataset_dir(build_dir: Path, dataset_row: dict[str, str]) -> Path:
    raw = dataset_row.get("dataset_dir", "")
    if raw:
        path = Path(raw)
        if path.exists():
            return path
        candidate = build_dir / raw
        if candidate.exists():
            return candidate
    dataset_name = dataset_row.get("dataset_name", "")
    add_order = dataset_row.get("add_order", "")
    name = f"{int(add_order):02d}_{dataset_name}" if add_order.isdigit() else dataset_name
    return build_dir / "datasets" / name


def _prepare_summary(dataset_dir: Path) -> dict[str, str]:
    return _first_row(dataset_dir / "prepare_summary.tsv")


def _first_row(path: Path) -> dict[str, str]:
    rows = _read_tsv(path)
    return rows[0] if rows else {}


def _active_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row.get("status", "active") in {"", "active"}]


def _is_active(row: dict[str, str]) -> bool:
    return row.get("status", "active") not in {"deprecated", "superseded"}


def _sequence_lengths(build_dir: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for path in _fasta_paths(build_dir):
        if not path.exists():
            continue
        seq_id = ""
        chunks: list[str] = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_id:
                        lengths.setdefault(seq_id, len("".join(chunks)))
                    seq_id = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if seq_id:
                lengths.setdefault(seq_id, len("".join(chunks)))
    return lengths


def _count_fasta_records(path: Path) -> str:
    if not path.exists():
        return "0"
    count = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return str(count)


def _fasta_paths(build_dir: Path) -> list[Path]:
    paths = [
        build_dir / "registry" / "current_representatives.fa",
        build_dir / "silva" / "silva_named_backbone.fa",
        build_dir / "silva" / "silva_unresolved.resolved.fa",
    ]
    datasets_dir = build_dir / "datasets"
    if datasets_dir.exists():
        for dataset_dir in sorted(path for path in datasets_dir.iterdir() if path.is_dir()):
            paths.extend(
                [
                    dataset_dir / "sina.oriented.fa",
                    dataset_dir / "prepared.ssu.fa",
                    dataset_dir / "input.normalized.fa",
                ]
            )
    return paths


def _ensure_writable(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Report already exists; use --overwrite: {path}")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))


def _write_tsv(rows: list[dict[str, Any]], path: Path, fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=chr(9),
            lineterminator=chr(10),
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _autotax2_version() -> str:
    try:
        return version("autotax2")
    except PackageNotFoundError:
        return "0.1.0"


def _int(value: str | None) -> int:
    try:
        return int(float(value or "0"))
    except ValueError:
        return 0


def _is_true(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "y"}


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
