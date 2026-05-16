"""Validation helpers and registry checks."""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from autotax2.export_validation import validate_export_dir

DATASET_PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
PLACEHOLDER_LIKE_RE = re.compile(r"^[pcofgs]__[A-Za-z][A-Za-z0-9]*[pcofgs][0-9]{6}$")
PLACEHOLDER_ID_RE = re.compile(
    r"^(?P<rank>[pcofgs])__(?P<prefix>[A-Za-z][A-Za-z0-9]*)(?P=rank)(?P<ordinal>[0-9]{6})$"
)
RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
RANK_INDEX = {rank: index for index, rank in enumerate(RANKS)}
REPORT_FIELDS = ["level", "check", "status", "message", "path"]


@dataclass(frozen=True)
class ValidationFinding:
    """One validation finding."""

    level: str
    check: str
    status: str
    message: str
    path: str
    strict_error: bool = False


@dataclass(frozen=True)
class ValidationSummary:
    """Summary returned by :func:`validate_build`."""

    build: Path
    report_md: Path
    report_tsv: Path
    errors: int
    warnings: int

    @property
    def failed(self) -> bool:
        """Return whether validation should fail."""
        return self.errors > 0


def validate_dataset_prefix(dataset_prefix: str) -> str:
    """Validate and normalize a dataset prefix."""
    normalized = dataset_prefix.strip()
    if not normalized:
        raise ValueError("Dataset prefix must not be empty.")
    if not DATASET_PREFIX_RE.fullmatch(normalized):
        raise ValueError(
            "Dataset prefix must start with a letter and contain only letters and digits."
        )
    return normalized


def validate_build(
    build: str | Path,
    strict: bool = False,
    check_exports: bool = True,
    report: str | Path | None = None,
) -> ValidationSummary:
    """Validate an autotax2 build and write Markdown/TSV validation reports."""
    build_dir = Path(build)
    report_md = Path(report) if report is not None else build_dir / "reports" / "validation_report.md"
    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_tsv = report_md.with_suffix(".tsv")
    findings: list[ValidationFinding] = []
    ctx = _ValidationContext(build_dir=build_dir, findings=findings, strict=strict)

    _check_registry_files(ctx)
    taxon_rows = _read_tsv(build_dir / "registry" / "taxon_nodes.tsv")
    taxon_by_id = {row.get("taxon_id", ""): row for row in taxon_rows if row.get("taxon_id")}
    _check_dataset_prefixes(ctx)
    _check_placeholders(ctx, taxon_rows)
    _check_taxonomy_and_tree(ctx, taxon_rows, taxon_by_id)
    _check_silva_immutability(ctx, taxon_rows)
    _check_sequence_mapping(ctx)
    _check_representatives(ctx, taxon_by_id)
    _check_silva_resolve_evidence(ctx)
    _check_dataset_evidence(ctx)
    if check_exports:
        _check_exports(ctx)
    _check_tool_metadata(ctx)

    _write_validation_reports(findings, report_md, report_tsv)
    return ValidationSummary(
        build=build_dir,
        report_md=report_md,
        report_tsv=report_tsv,
        errors=sum(1 for finding in findings if finding.level == "error"),
        warnings=sum(1 for finding in findings if finding.level == "warning"),
    )


@dataclass
class _ValidationContext:
    build_dir: Path
    findings: list[ValidationFinding]
    strict: bool

    def ok(self, check: str, message: str, path: str | Path = "") -> None:
        self.findings.append(ValidationFinding("info", check, "ok", message, str(path)))

    def warning(self, check: str, message: str, path: str | Path = "", strict_error: bool = False) -> None:
        level = "error" if self.strict and strict_error else "warning"
        self.findings.append(ValidationFinding(level, check, "warning", message, str(path), strict_error))

    def error(self, check: str, message: str, path: str | Path = "") -> None:
        self.findings.append(ValidationFinding("error", check, "error", message, str(path)))


def _check_registry_files(ctx: _ValidationContext) -> None:
    registry_dir = ctx.build_dir / "registry"
    required = [
        "taxon_nodes.tsv",
        "sequence_registry.tsv",
        "dataset_registry.tsv",
        "name_index.tsv",
        "placeholder_counters.yaml",
    ]
    for filename in required:
        path = registry_dir / filename
        if path.exists():
            ctx.ok("registry_files", f"Found {filename}.", path)
        else:
            ctx.error("registry_files", f"Missing required registry file: {filename}.", path)
    for filename in ("cluster_to_taxon.tsv", "representative_registry.tsv"):
        path = registry_dir / filename
        if path.exists():
            ctx.ok("registry_files", f"Found optional registry file: {filename}.", path)
        else:
            ctx.warning("registry_files", f"Optional registry file is absent: {filename}.", path)


def _check_dataset_prefixes(ctx: _ValidationContext) -> None:
    path = ctx.build_dir / "registry" / "dataset_registry.tsv"
    rows = _read_tsv(path)
    prefixes: dict[str, str] = {}
    datasets: dict[str, str] = {}
    for row in rows:
        dataset = row.get("dataset_name", "")
        prefix = row.get("prefix", "")
        if not dataset or not prefix:
            continue
        if prefix == "SILVA":
            ctx.error("dataset_prefixes", "Custom dataset cannot use reserved prefix SILVA.", path)
        if prefix in prefixes and prefixes[prefix] != dataset:
            ctx.error("dataset_prefixes", f"Prefix {prefix} is assigned to multiple datasets.", path)
        prefixes[prefix] = dataset
        if dataset in datasets and datasets[dataset] != prefix:
            ctx.error("dataset_prefixes", f"Dataset {dataset} has conflicting frozen prefixes.", path)
        datasets[dataset] = prefix
    if len(prefixes) == len([row for row in rows if row.get("dataset_name") and row.get("prefix")]):
        ctx.ok("dataset_prefixes", "Dataset prefixes are unique.", path)


def _check_placeholders(ctx: _ValidationContext, taxon_rows: list[dict[str, str]]) -> None:
    active_names_by_rank: set[tuple[str, str]] = set()
    active_placeholder_names: set[str] = set()
    inactive_placeholder_names: set[str] = set()
    for row in taxon_rows:
        name = row.get("name", "")
        rank = row.get("rank", "")
        status = row.get("status", "active")
        is_active = status not in {"deprecated", "superseded"}
        key = (rank, name)
        if is_active:
            if key in active_names_by_rank:
                ctx.error("placeholder_uniqueness", f"Duplicate active taxon name at rank {rank}: {name}.")
            active_names_by_rank.add(key)
        is_placeholder = _is_true(row.get("is_placeholder", "false")) or PLACEHOLDER_LIKE_RE.fullmatch(name) is not None
        if not is_placeholder:
            continue
        if PLACEHOLDER_ID_RE.fullmatch(name) is None:
            ctx.error("placeholder_format", f"Invalid placeholder format: {name}.")
            continue
        if is_active:
            if name in active_placeholder_names:
                ctx.error("placeholder_uniqueness", f"Duplicate active placeholder name: {name}.")
            active_placeholder_names.add(name)
        else:
            inactive_placeholder_names.add(name)
    reused = active_placeholder_names & inactive_placeholder_names
    for name in sorted(reused):
        ctx.error("placeholder_reuse", f"Deprecated or superseded placeholder is reused by active taxon: {name}.")
    if not reused:
        ctx.ok("placeholder_reuse", "Deprecated placeholders are not reused.")


def _check_taxonomy_and_tree(
    ctx: _ValidationContext,
    taxon_rows: list[dict[str, str]],
    taxon_by_id: dict[str, dict[str, str]],
) -> None:
    for row in taxon_rows:
        taxon_id = row.get("taxon_id", "")
        rank = row.get("rank", "")
        parent_id = row.get("parent_taxon_id", "")
        if not taxon_id or rank not in RANK_INDEX:
            continue
        if rank == "domain":
            continue
        if not parent_id:
            ctx.error("parent_child", f"Non-root taxon has no parent: {taxon_id}.")
            continue
        parent = taxon_by_id.get(parent_id)
        if parent is None:
            ctx.error("parent_child", f"Taxon {taxon_id} references missing parent {parent_id}.")
            continue
        parent_rank = parent.get("rank", "")
        if RANK_INDEX.get(parent_rank, -99) != RANK_INDEX[rank] - 1:
            ctx.error("parent_child", f"Parent rank for {taxon_id} is not immediately above {rank}.")

    for row in taxon_rows:
        taxon_id = row.get("taxon_id", "")
        if not taxon_id:
            continue
        cycle = _detect_cycle(taxon_id, taxon_by_id)
        if cycle:
            ctx.error("taxon_cycles", f"Taxon tree contains a cycle: {' -> '.join(cycle)}.")

    for row in taxon_rows:
        if row.get("rank") != "species" or row.get("status", "active") in {"deprecated", "superseded"}:
            continue
        lineage = _lineage_ranks(row.get("taxon_id", ""), taxon_by_id)
        if lineage != list(RANKS):
            ctx.error("taxonomy_paths", f"Active species taxon does not resolve to exactly 7 ranks: {row.get('taxon_id')}.")
    ctx.ok("taxonomy_paths", "Taxonomy path checks completed.")


def _check_silva_immutability(ctx: _ValidationContext, taxon_rows: list[dict[str, str]]) -> None:
    registry_dir = ctx.build_dir / "registry"
    snapshot_path = registry_dir / "protected_taxa_snapshot.tsv"
    protected = [
        {
            "taxon_id": row.get("taxon_id", ""),
            "rank": row.get("rank", ""),
            "name": row.get("name", ""),
            "parent_taxon_id": row.get("parent_taxon_id", ""),
        }
        for row in taxon_rows
        if _is_true(row.get("is_silva_named", "false"))
        or (
            _is_true(row.get("protected", "false"))
            and row.get("source", "").startswith("SILVA")
            and not _is_true(row.get("is_silva_unresolved", "false"))
        )
    ]
    for row in protected:
        source_row = next((item for item in taxon_rows if item.get("taxon_id") == row["taxon_id"]), {})
        if not _is_true(source_row.get("protected", "false")):
            ctx.error("silva_immutability", f"Named SILVA taxon is not protected: {row['taxon_id']}.")
        if PLACEHOLDER_LIKE_RE.fullmatch(row["name"]):
            ctx.error("silva_immutability", f"Named SILVA taxon was renamed to placeholder: {row['taxon_id']}.")
    if snapshot_path.exists():
        existing = _read_tsv(snapshot_path)
        if existing and existing != protected:
            ctx.error("silva_immutability", "Protected SILVA taxon snapshot differs from current registry.", snapshot_path)
    else:
        _write_tsv(protected, snapshot_path, ["taxon_id", "rank", "name", "parent_taxon_id"])
        ctx.warning("silva_immutability", "Created protected SILVA taxon snapshot for future comparisons.", snapshot_path)


def _check_sequence_mapping(ctx: _ValidationContext) -> None:
    sequence_path = ctx.build_dir / "registry" / "sequence_registry.tsv"
    for row in _read_tsv(sequence_path):
        seq_id = row.get("internal_seq_id") or row.get("seq_id", "")
        if not row.get("sequence_md5"):
            ctx.error("sequence_mapping", f"Sequence lacks sequence_md5: {seq_id}.", sequence_path)
        is_custom = bool(row.get("dataset")) and not row.get("source", "").startswith("SILVA")
        if is_custom and not row.get("original_seq_id"):
            ctx.error("sequence_mapping", f"Internal sequence lacks original_seq_id: {seq_id}.", sequence_path)

    for dataset_dir in _dataset_dirs(ctx.build_dir):
        id_map = _read_tsv(dataset_dir / "sequence_id_map.tsv")
        mapped_ids = {row.get("internal_seq_id", "") for row in id_map}
        for row in id_map:
            if not row.get("internal_seq_id") or not row.get("original_seq_id"):
                ctx.error("sequence_mapping", "Dataset sequence_id_map has incomplete ID mapping.", dataset_dir / "sequence_id_map.tsv")
            if not row.get("sequence_md5"):
                ctx.error("sequence_mapping", f"Dataset mapped sequence lacks sequence_md5: {row.get('internal_seq_id', '')}.", dataset_dir / "sequence_id_map.tsv")
        for member in _read_tsv(dataset_dir / "sequence_membership.tsv"):
            if member.get("internal_seq_id", "") not in mapped_ids:
                ctx.error("sequence_mapping", f"Membership lacks sequence_id_map row: {member.get('internal_seq_id', '')}.", dataset_dir / "sequence_membership.tsv")
            if _is_true(member.get("is_duplicate_sequence", "false")) and not member.get("unique_seq_id"):
                ctx.error("sequence_mapping", f"Duplicate MD5 membership lacks unique_seq_id: {member.get('internal_seq_id', '')}.")


def _check_representatives(ctx: _ValidationContext, taxon_by_id: dict[str, dict[str, str]]) -> None:
    rep_path = ctx.build_dir / "registry" / "representative_registry.tsv"
    reps = [row for row in _read_tsv(rep_path) if row.get("status", "active") in {"", "active"}]
    sequence_ids = _sequence_ids(ctx.build_dir)
    represented_taxa = {row.get("taxon_id", "") for row in reps}
    for rep in reps:
        seq_id = rep.get("representative_seq_id", "")
        if seq_id not in sequence_ids:
            ctx.error("representatives", f"Representative sequence does not exist: {seq_id}.", rep_path)
    for taxon_id, row in taxon_by_id.items():
        if row.get("rank") == "species" and row.get("status", "active") not in {"deprecated", "superseded"}:
            if taxon_id not in represented_taxa:
                ctx.warning("representatives", f"Active species lacks representative: {taxon_id}.", rep_path, strict_error=True)


def _check_exports(ctx: _ValidationContext) -> None:
    export_dir = ctx.build_dir / "export"
    if not export_dir.exists():
        ctx.warning("exports", "No export directory found; skipping export file validation.", export_dir)
        return
    for finding in validate_export_dir(export_dir):
        if finding.level == "error":
            ctx.error(finding.check, finding.message, finding.path)
        elif finding.level == "warning":
            ctx.warning(finding.check, finding.message, finding.path)
        else:
            ctx.ok(finding.check, finding.message, finding.path)


def _check_tool_metadata(ctx: _ValidationContext) -> None:
    for dataset_dir in _dataset_dirs(ctx.build_dir):
        tool_versions = _read_tsv(dataset_dir / "tool_versions.tsv")
        tools = {row.get("tool", "") for row in tool_versions}
        if (dataset_dir / "sina.oriented.fa").exists() and "sina" not in tools:
            ctx.warning("tool_metadata", f"SINA version/status was not recorded for {dataset_dir.name}.", dataset_dir / "tool_versions.tsv", strict_error=True)
        cluster_summary = _read_tsv(dataset_dir / "cluster_search_summary.tsv")
        if cluster_summary and not cluster_summary[0].get("iddef"):
            ctx.warning("tool_metadata", f"VSEARCH iddef was not recorded for {dataset_dir.name}.", dataset_dir / "cluster_search_summary.tsv", strict_error=True)
        if (dataset_dir / "sina.candidates.tsv").exists() and cluster_summary:
            summary_row = cluster_summary[0]
            if not summary_row.get("sina_candidate_source"):
                ctx.warning("tool_metadata", f"SINA candidate source was not recorded for {dataset_dir.name}.", dataset_dir / "cluster_search_summary.tsv", strict_error=True)
            if _int(summary_row.get("sina_candidate_targets")) > 0 and _int(summary_row.get("sina_candidate_target_matches")) == 0:
                ctx.warning("tool_metadata", f"SINA candidates had no matching current representatives for {dataset_dir.name}; VSEARCH likely used full registry fallback.", dataset_dir / "cluster_search_summary.tsv", strict_error=True)


def _check_dataset_evidence(ctx: _ValidationContext) -> None:
    for dataset_dir in _dataset_dirs(ctx.build_dir):
        assignments_path = dataset_dir / "assignments.tsv"
        assignments = _read_tsv(assignments_path)
        if not assignments:
            continue
        evidence_path = dataset_dir / "placement_evidence.tsv"
        evidence = _read_tsv(evidence_path)
        if not evidence:
            ctx.warning("placement_evidence", f"Missing placement evidence for {dataset_dir.name}.", evidence_path, strict_error=True)
            continue
        query_ids = {row.get("internal_seq_id", "") for row in assignments if row.get("internal_seq_id")}
        evidence_by_query: dict[str, set[str]] = {}
        for row in evidence:
            seq_id = row.get("seq_id", "")
            rank = row.get("rank", "")
            if seq_id and rank:
                evidence_by_query.setdefault(seq_id, set()).add(rank)
        missing_queries = sorted(query_id for query_id in query_ids if query_id not in evidence_by_query)
        for query_id in missing_queries:
            ctx.warning("placement_evidence", f"Assignment lacks placement evidence rows: {query_id}.", evidence_path, strict_error=True)
        for query_id, ranks in sorted(evidence_by_query.items()):
            missing_ranks = [rank for rank in RANKS[1:] if rank not in ranks]
            if missing_ranks:
                ctx.warning("placement_evidence", f"Placement evidence for {query_id} lacks ranks: {','.join(missing_ranks)}.", evidence_path, strict_error=True)
        if not missing_queries and all(all(rank in ranks for rank in RANKS[1:]) for ranks in evidence_by_query.values()):
            ctx.ok("placement_evidence", f"Placement evidence covers assignments for {dataset_dir.name}.", evidence_path)


def _check_silva_resolve_evidence(ctx: _ValidationContext) -> None:
    silva_dir = ctx.build_dir / "silva"
    members_path = silva_dir / "silva_unresolved_members.tsv"
    members = _read_tsv(members_path)
    if not members:
        return
    evidence_path = silva_dir / "silva_unresolved_evidence.tsv"
    evidence = _read_tsv(evidence_path)
    if not evidence:
        ctx.warning("silva_resolve_evidence", "Missing SILVA unresolved resolve evidence.", evidence_path, strict_error=True)
        return
    member_ids = {row.get("seq_id", "") for row in members if row.get("seq_id")}
    evidence_by_seq: dict[str, set[str]] = {}
    for row in evidence:
        seq_id = row.get("seq_id", "")
        rank = row.get("rank", "")
        if seq_id and rank:
            evidence_by_seq.setdefault(seq_id, set()).add(rank)
    missing_members = sorted(seq_id for seq_id in member_ids if seq_id not in evidence_by_seq)
    for seq_id in missing_members:
        ctx.warning("silva_resolve_evidence", f"SILVA unresolved member lacks evidence rows: {seq_id}.", evidence_path, strict_error=True)
    for seq_id, ranks in sorted(evidence_by_seq.items()):
        missing_ranks = [rank for rank in RANKS[1:] if rank not in ranks]
        if missing_ranks:
            ctx.warning("silva_resolve_evidence", f"SILVA resolve evidence for {seq_id} lacks ranks: {','.join(missing_ranks)}.", evidence_path, strict_error=True)
    if not missing_members and all(all(rank in ranks for rank in RANKS[1:]) for ranks in evidence_by_seq.values()):
        ctx.ok("silva_resolve_evidence", "SILVA unresolved resolve evidence covers all members.", evidence_path)


def _write_validation_reports(
    findings: list[ValidationFinding],
    report_md: Path,
    report_tsv: Path,
) -> None:
    error_count = sum(1 for finding in findings if finding.level == "error")
    warning_count = sum(1 for finding in findings if finding.level == "warning")
    sections = [
        "# autotax2 Validation Report",
        "",
        "## Summary",
        "",
        f"- Errors: {error_count}",
        f"- Warnings: {warning_count}",
        "",
        "## Errors",
        "",
        *_markdown_findings(findings, "error"),
        "",
        "## Warnings",
        "",
        *_markdown_findings(findings, "warning"),
        "",
    ]
    for title, prefix in [
        ("Registry checks", "registry"),
        ("Placeholder checks", "placeholder"),
        ("Taxonomy checks", "taxonomy"),
        ("Sequence checks", "sequence"),
        ("Export checks", "export"),
        ("Tool metadata checks", "tool"),
    ]:
        sections.extend([f"## {title}", ""])
        matched = [finding for finding in findings if finding.check.startswith(prefix) or prefix in finding.check]
        sections.extend(_markdown_finding_lines(matched) or ["- No findings."])
        sections.append("")
    report_md.write_text("\n".join(sections), encoding="utf-8")
    _write_tsv(
        [
            {
                "level": finding.level,
                "check": finding.check,
                "status": finding.status,
                "message": finding.message,
                "path": finding.path.replace("\\", "/"),
            }
            for finding in findings
        ],
        report_tsv,
        REPORT_FIELDS,
    )


def _markdown_findings(findings: list[ValidationFinding], level: str) -> list[str]:
    matched = [finding for finding in findings if finding.level == level]
    return _markdown_finding_lines(matched) or ["- None."]


def _markdown_finding_lines(findings: list[ValidationFinding]) -> list[str]:
    return [
        f"- **{finding.check}**: {finding.message}"
        + (f" (`{finding.path}`)" if finding.path else "")
        for finding in findings
    ]


def _detect_cycle(taxon_id: str, taxon_by_id: dict[str, dict[str, str]]) -> list[str]:
    seen: list[str] = []
    current = taxon_id
    while current:
        if current in seen:
            return seen[seen.index(current) :] + [current]
        seen.append(current)
        current = taxon_by_id.get(current, {}).get("parent_taxon_id", "")
    return []


def _lineage_ranks(taxon_id: str, taxon_by_id: dict[str, dict[str, str]]) -> list[str]:
    ranks: list[str] = []
    seen: set[str] = set()
    current = taxon_id
    while current and current not in seen:
        seen.add(current)
        row = taxon_by_id.get(current, {})
        if not row:
            break
        ranks.append(row.get("rank", ""))
        current = row.get("parent_taxon_id", "")
    return list(reversed(ranks))


def _dataset_dirs(build_dir: Path) -> list[Path]:
    datasets_dir = build_dir / "datasets"
    if not datasets_dir.exists():
        return []
    return sorted(path for path in datasets_dir.iterdir() if path.is_dir())


def _sequence_ids(build_dir: Path) -> set[str]:
    ids: set[str] = set()
    for path in _fasta_paths(build_dir):
        if not path.exists():
            continue
        ids.update(_fasta_header_ids(path))
    for row in _read_tsv(build_dir / "registry" / "sequence_registry.tsv"):
        ids.add(row.get("seq_id") or row.get("internal_seq_id", ""))
    return {seq_id for seq_id in ids if seq_id}


def _fasta_paths(build_dir: Path) -> list[Path]:
    paths = [
        build_dir / "registry" / "current_representatives.fa",
        build_dir / "silva" / "silva_named_backbone.fa",
        build_dir / "silva" / "silva_unresolved.resolved.fa",
    ]
    for dataset_dir in _dataset_dirs(build_dir):
        paths.extend([dataset_dir / "sina.oriented.fa", dataset_dir / "prepared.ssu.fa", dataset_dir / "input.normalized.fa"])
    export_dir = build_dir / "export"
    if export_dir.exists():
        paths.extend(export_dir.glob("**/*.fa"))
        paths.extend(export_dir.glob("**/*.fa.gz"))
        paths.extend(export_dir.glob("**/*.fasta"))
        paths.extend(export_dir.glob("**/*.fasta.gz"))
    return paths


def _fasta_header_ids(path: Path) -> list[str]:
    return [header.split()[0].split(";tax=", maxsplit=1)[0] for header in _fasta_headers(path)]


def _fasta_headers(path: Path) -> list[str]:
    headers: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:  # type: ignore[arg-type]
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                headers.append(line[1:])
    return headers


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


def _is_true(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "y"}


def _int(value: str | None) -> int:
    try:
        return int(float(value or "0"))
    except ValueError:
        return 0
