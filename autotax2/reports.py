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
        assignments = _read_tsv(dataset_dir / "assignments.tsv")
        members = _read_tsv(dataset_dir / "sequence_membership.tsv")
        created = _read_tsv(dataset_dir / "created_taxa.tsv")
        reps = _read_tsv(dataset_dir / "representative_updates.tsv")
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
                if taxon_id and _taxon_source_type(taxon_id, taxon_by_id, dataset) != "current_dataset":
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
