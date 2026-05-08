"""Dataset preparation with mandatory barrnap recutting."""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from autotax2.io import FastaRecord, read_fasta, write_fasta
from autotax2.registry import (
    AssignedSequence,
    SequenceRegistry,
    assign_internal_ids,
    sequence_md5,
)
from autotax2.validate import validate_dataset_prefix


BARRNAP_VERSION = "1.10.5"
TOOL_VERSION_FIELDS = [
    "tool",
    "expected_version",
    "detected_version",
    "status",
    "command",
]
DATASET_REGISTRY_FIELDS = [
    "dataset_name",
    "prefix",
    "add_order",
    "input_fasta",
    "input_md5",
    "domain",
    "status",
    "dataset_dir",
]
SEQUENCE_ID_MAP_FIELDS = [
    "internal_seq_id",
    "original_seq_id",
    "original_header",
    "dataset",
    "prefix",
    "input_order",
    "sequence_md5",
    "input_length",
    "normalized_length",
    "rejected",
    "reject_reason",
]
BARRNAP_SUMMARY_FIELDS = [
    "internal_seq_id",
    "extracted_seq_id",
    "original_seq_id",
    "dataset",
    "domain",
    "input_length",
    "hit_count",
    "selected_hit_index",
    "barrnap_type",
    "barrnap_start",
    "barrnap_end",
    "barrnap_strand",
    "barrnap_score",
    "extracted_start",
    "extracted_end",
    "extracted_length",
    "reverse_complemented",
    "status",
    "warning",
]
PREPARE_SUMMARY_FIELDS = [
    "dataset",
    "prefix",
    "input_sequences",
    "normalized_sequences",
    "rejected_non_atgc",
    "duplicate_md5_sequences",
    "barrnap_extracted",
    "no_barrnap_hit",
    "multiple_rrna_hits",
    "rejected_short_rrna",
    "final_prepared_sequences",
]
ATGC_RE = re.compile(r"^[ATGC]*$")


@dataclass(frozen=True)
class BarrnapHit:
    """One rRNA feature parsed from barrnap GFF3."""

    seq_id: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: dict[str, str]

    @property
    def length(self) -> int:
        """Return 1-based inclusive feature length."""
        return self.end - self.start + 1

    @property
    def numeric_score(self) -> float | None:
        """Return numeric score if present."""
        try:
            return float(self.score)
        except ValueError:
            return None


@dataclass(frozen=True)
class PrepareDatasetSummary:
    """Summary from dataset preparation."""

    dataset: str
    prefix: str
    dataset_dir: Path
    final_prepared_sequences: int


def recutting_required() -> bool:
    """Return whether barrnap recutting is required."""
    return True


def check_barrnap_version(barrnap_bin: str = "barrnap") -> str:
    """Return the detected barrnap version string, or an empty string if unavailable."""
    try:
        completed = subprocess.run(
            [barrnap_bin, "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return ""
    output = f"{completed.stdout}\n{completed.stderr}".strip()
    match = re.search(r"barrnap\s+v?(\d+(?:\.\d+)+)", output, re.I)
    if match:
        return match.group(1)
    match = re.search(r"(\d+(?:\.\d+)+)", output)
    return match.group(1) if match else ""


def prepare_dataset(
    build: str | Path,
    name: str,
    prefix: str,
    fasta: str | Path,
    domain: str,
    threads: int = 4,
    barrnap_bin: str = "barrnap",
    barrnap_kingdom: str | None = None,
    strict_tool_version: bool = False,
    multi_rrna_policy: str = "longest",
    min_rrna_len_archaea: int = 900,
    min_rrna_len_bacteria: int = 1200,
    flank: int = 0,
    allow_partial: bool = False,
    reject_non_atgc: bool = True,
) -> PrepareDatasetSummary:
    """Prepare a custom dataset by normalizing input and recutting rRNA with barrnap."""
    del allow_partial

    build_dir = Path(build)
    registry_dir = build_dir / "registry"
    registry_dir.mkdir(parents=True, exist_ok=True)
    normalized_name = _validate_dataset_name(name)
    normalized_prefix = validate_dataset_prefix(prefix)
    normalized_domain = _validate_domain(domain)
    policy = _validate_policy(multi_rrna_policy)
    kingdom = barrnap_kingdom or _default_barrnap_kingdom(normalized_domain)
    input_fasta = Path(fasta)
    input_md5 = _file_md5(input_fasta)
    registry_path = registry_dir / "dataset_registry.tsv"
    dataset_rows = _read_tsv(registry_path)
    add_order = _resolve_dataset_registration(
        rows=dataset_rows,
        name=normalized_name,
        prefix=normalized_prefix,
    )
    dataset_dir = build_dir / "datasets" / f"{add_order:02d}_{_safe_dataset_dir_name(normalized_name)}"
    if dataset_dir.exists():
        raise FileExistsError(
            f"Dataset output directory already exists and would be overwritten: {dataset_dir}"
        )
    dataset_dir.mkdir(parents=True)
    detected_barrnap_version = check_barrnap_version(barrnap_bin)
    if strict_tool_version and detected_barrnap_version != BARRNAP_VERSION:
        raise RuntimeError(
            f"Expected barrnap {BARRNAP_VERSION}, detected "
            f"{detected_barrnap_version or 'unknown'} from {barrnap_bin}."
        )

    source_records = read_fasta(input_fasta)
    assigned = assign_internal_ids(source_records, normalized_name, normalized_prefix)
    sequence_map_rows: list[dict[str, str]] = []
    accepted_assigned: list[AssignedSequence] = []
    accepted_records: list[FastaRecord] = []
    rejected_non_atgc = 0

    for input_order, (source_record, assigned_record) in enumerate(
        zip(source_records, assigned, strict=True),
        start=1,
    ):
        reject_reason = ""
        rejected = False
        if reject_non_atgc and not ATGC_RE.fullmatch(assigned_record.sequence):
            rejected = True
            reject_reason = "non_atgc"
            rejected_non_atgc += 1
        sequence_map_rows.append(
            {
                "internal_seq_id": assigned_record.internal_seq_id,
                "original_seq_id": assigned_record.original_seq_id,
                "original_header": assigned_record.original_header,
                "dataset": assigned_record.dataset,
                "prefix": assigned_record.prefix,
                "input_order": str(input_order),
                "sequence_md5": assigned_record.sequence_md5,
                "input_length": str(len(source_record.sequence)),
                "normalized_length": str(assigned_record.sequence_length),
                "rejected": _bool_text(rejected),
                "reject_reason": reject_reason,
            }
        )
        if not rejected:
            accepted_assigned.append(assigned_record)
            accepted_records.append(
                FastaRecord(
                    seq_id=assigned_record.internal_seq_id,
                    header=assigned_record.internal_seq_id,
                    sequence=assigned_record.sequence,
                )
            )

    write_fasta(accepted_records, dataset_dir / "input.normalized.fa")
    _write_tsv(sequence_map_rows, dataset_dir / "sequence_id_map.tsv", SEQUENCE_ID_MAP_FIELDS)

    sequence_registry = SequenceRegistry()
    memberships = sequence_registry.add_many(accepted_assigned)
    _write_tsv(
        [record.__dict__ for record in sequence_registry.unique_records()],
        dataset_dir / "unique_sequence_registry.tsv",
        [
            "unique_seq_id",
            "sequence_md5",
            "representative_internal_seq_id",
            "sequence",
            "length",
            "first_seen_dataset",
        ],
    )
    _write_tsv(
        [
            {
                "internal_seq_id": record.internal_seq_id,
                "original_seq_id": record.original_seq_id,
                "dataset": record.dataset,
                "prefix": record.prefix,
                "sequence_md5": record.sequence_md5,
                "unique_seq_id": record.unique_seq_id,
                "is_duplicate_sequence": _bool_text(record.is_duplicate_sequence),
            }
            for record in sequence_registry.membership_records()
        ],
        dataset_dir / "sequence_membership.tsv",
        [
            "internal_seq_id",
            "original_seq_id",
            "dataset",
            "prefix",
            "sequence_md5",
            "unique_seq_id",
            "is_duplicate_sequence",
        ],
    )

    barrnap_gff = dataset_dir / "barrnap.gff3"
    barrnap_log = dataset_dir / "barrnap.log"
    _run_barrnap_if_needed(
        input_fasta=dataset_dir / "input.normalized.fa",
        output_gff=barrnap_gff,
        log_path=barrnap_log,
        kingdom=kingdom,
        threads=threads,
        barrnap_bin=barrnap_bin,
    )
    _write_tool_version(
        dataset_dir / "tool_versions.tsv",
        detected_version=detected_barrnap_version,
        status="ok" if detected_barrnap_version == BARRNAP_VERSION else "warning_unparsed",
        command=[
            barrnap_bin,
            "--kingdom",
            kingdom,
            "--threads",
            str(threads),
            str(dataset_dir / "input.normalized.fa"),
        ],
    )

    hits_by_seq = _hits_by_seq_id(parse_barrnap_gff3(barrnap_gff))
    min_rrna_len = min_rrna_len_archaea if normalized_domain == "Archaea" else min_rrna_len_bacteria
    barrnap_summary_rows: list[dict[str, str]] = []
    extracted_records: list[FastaRecord] = []
    assigned_by_internal_id = {record.internal_seq_id: record for record in accepted_assigned}
    original_id_by_internal_id = {
        row["internal_seq_id"]: row["original_seq_id"] for row in sequence_map_rows
    }
    normalized_by_internal_id = {record.seq_id: record for record in accepted_records}

    for row in sequence_map_rows:
        if row["rejected"] == "true":
            barrnap_summary_rows.append(
                _barrnap_summary_row(
                    internal_seq_id=row["internal_seq_id"],
                    original_seq_id=row["original_seq_id"],
                    dataset=normalized_name,
                    domain=normalized_domain,
                    input_length=row["normalized_length"],
                    status="rejected_non_atgc",
                )
            )

    for record in accepted_records:
        assigned_record = assigned_by_internal_id[record.seq_id]
        original_seq_id = original_id_by_internal_id[record.seq_id]
        hits = hits_by_seq.get(record.seq_id, [])
        selected_hits = _select_hits(hits, policy)
        if not hits:
            barrnap_summary_rows.append(
                _barrnap_summary_row(
                    internal_seq_id=record.seq_id,
                    original_seq_id=original_seq_id,
                    dataset=normalized_name,
                    domain=normalized_domain,
                    input_length=str(assigned_record.sequence_length),
                    hit_count=0,
                    status="no_barrnap_hit",
                )
            )
            continue
        if len(hits) > 1 and policy == "fail":
            barrnap_summary_rows.append(
                _barrnap_summary_row(
                    internal_seq_id=record.seq_id,
                    original_seq_id=original_seq_id,
                    dataset=normalized_name,
                    domain=normalized_domain,
                    input_length=str(assigned_record.sequence_length),
                    hit_count=len(hits),
                    status="rejected_multiple_rrna_hits",
                    warning="multiple_rrna_hits",
                )
            )
            continue

        for selected_index, hit in selected_hits:
            extracted = _extract_hit(record.sequence, hit, flank)
            extracted_seq_id = (
                record.seq_id
                if len(selected_hits) == 1
                else f"{record.seq_id}_rrna{selected_index + 1}"
            )
            status = "extracted"
            warning = "multiple_rrna_hits" if len(hits) > 1 else ""
            if len(extracted["sequence"]) < min_rrna_len:
                status = "rejected_short_rrna"
                warning = _append_warning(warning, "short_rrna")
            else:
                extracted_records.append(
                    FastaRecord(
                        seq_id=extracted_seq_id,
                        header=extracted_seq_id,
                        sequence=extracted["sequence"],
                    )
                )
            barrnap_summary_rows.append(
                _barrnap_summary_row(
                    internal_seq_id=record.seq_id,
                    extracted_seq_id=extracted_seq_id,
                    original_seq_id=original_seq_id,
                    dataset=normalized_name,
                    domain=normalized_domain,
                    input_length=str(assigned_record.sequence_length),
                    hit_count=len(hits),
                    selected_hit_index=selected_index + 1,
                    hit=hit,
                    extracted_start=extracted["start"],
                    extracted_end=extracted["end"],
                    extracted_length=len(extracted["sequence"]),
                    reverse_complemented=hit.strand == "-",
                    status=status,
                    warning=warning,
                )
            )

    write_fasta(extracted_records, dataset_dir / "barrnap.extracted.fa")
    _write_tsv(
        barrnap_summary_rows,
        dataset_dir / "barrnap.summary.tsv",
        BARRNAP_SUMMARY_FIELDS,
    )
    duplicate_md5_sequences = sum(1 for membership in memberships if membership.is_duplicate_sequence)
    prepare_summary = {
        "dataset": normalized_name,
        "prefix": normalized_prefix,
        "input_sequences": str(len(source_records)),
        "normalized_sequences": str(len(accepted_records)),
        "rejected_non_atgc": str(rejected_non_atgc),
        "duplicate_md5_sequences": str(duplicate_md5_sequences),
        "barrnap_extracted": str(_count_status(barrnap_summary_rows, "extracted")),
        "no_barrnap_hit": str(_count_status(barrnap_summary_rows, "no_barrnap_hit")),
        "multiple_rrna_hits": str(
            sum(1 for row in barrnap_summary_rows if "multiple_rrna_hits" in row["warning"])
        ),
        "rejected_short_rrna": str(_count_status(barrnap_summary_rows, "rejected_short_rrna")),
        "final_prepared_sequences": str(len(extracted_records)),
    }
    _write_tsv([prepare_summary], dataset_dir / "prepare_summary.tsv", PREPARE_SUMMARY_FIELDS)
    _write_dataset_registry(
        registry_path=registry_path,
        rows=dataset_rows,
        dataset_row={
            "dataset_name": normalized_name,
            "prefix": normalized_prefix,
            "add_order": str(add_order),
            "input_fasta": str(input_fasta),
            "input_md5": input_md5,
            "domain": normalized_domain,
            "status": "prepared",
            "dataset_dir": str(dataset_dir),
        },
    )

    return PrepareDatasetSummary(
        dataset=normalized_name,
        prefix=normalized_prefix,
        dataset_dir=dataset_dir,
        final_prepared_sequences=len(extracted_records),
    )


def parse_barrnap_gff3(path: str | Path) -> list[BarrnapHit]:
    """Parse barrnap GFF3 and return flexible SSU/rRNA hits."""
    hits: list[BarrnapHit] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            attributes = _parse_gff3_attributes(fields[8])
            hit = BarrnapHit(
                seq_id=fields[0],
                source=fields[1],
                feature_type=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                score=fields[5],
                strand=fields[6],
                phase=fields[7],
                attributes=attributes,
            )
            if _is_ssu_rrna_hit(hit):
                hits.append(hit)
    return hits


def _run_barrnap_if_needed(
    input_fasta: Path,
    output_gff: Path,
    log_path: Path,
    kingdom: str,
    threads: int,
    barrnap_bin: str,
) -> None:
    if output_gff.exists():
        log_path.write_text("Using existing mocked barrnap.gff3\n", encoding="utf-8")
        return
    command = [
        barrnap_bin,
        "--kingdom",
        kingdom,
        "--threads",
        str(threads),
        str(input_fasta),
    ]
    try:
        with output_gff.open("w", encoding="utf-8", newline="") as stdout_handle:
            completed = subprocess.run(
                command,
                check=True,
                stdout=stdout_handle,
                stderr=subprocess.PIPE,
                text=True,
            )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"barrnap executable not found: {barrnap_bin}. "
            f"Install barrnap {BARRNAP_VERSION} or pass --barrnap-bin."
        ) from exc
    except subprocess.CalledProcessError as exc:
        log_path.write_text(exc.stderr or "", encoding="utf-8")
        raise RuntimeError(f"barrnap failed with exit code {exc.returncode}: {barrnap_bin}") from exc
    log_path.write_text(completed.stderr or "", encoding="utf-8")


def _resolve_dataset_registration(
    rows: list[dict[str, str]],
    name: str,
    prefix: str,
) -> int:
    existing_orders: list[int] = []
    for row in rows:
        row_name = row.get("dataset_name", "")
        row_prefix = row.get("prefix", "")
        row_order = row.get("add_order", "")
        if row_order.isdigit():
            existing_orders.append(int(row_order))
        if row_name == name and row_prefix and row_prefix != prefix:
            raise ValueError(f"Dataset name {name!r} is already registered with prefix {row_prefix!r}.")
        if row_prefix == prefix and row_name and row_name != name:
            raise ValueError(f"Dataset prefix {prefix!r} is already registered to {row_name!r}.")
        if row_name == name and row_prefix == prefix and row_order.isdigit():
            return int(row_order)
    return (max(existing_orders) + 1) if existing_orders else 1


def _write_dataset_registry(
    registry_path: Path,
    rows: list[dict[str, str]],
    dataset_row: dict[str, str],
) -> None:
    updated = False
    output_rows: list[dict[str, str]] = []
    for row in rows:
        if row.get("dataset_name") == dataset_row["dataset_name"]:
            output_rows.append(dataset_row)
            updated = True
        else:
            output_rows.append(row)
    if not updated:
        output_rows.append(dataset_row)
    fieldnames = _merge_fieldnames(DATASET_REGISTRY_FIELDS, output_rows)
    _write_tsv(output_rows, registry_path, fieldnames)


def _write_tool_version(
    tool_versions_path: Path,
    detected_version: str,
    status: str,
    command: list[str],
) -> None:
    rows = [
        row
        for row in _read_tsv(tool_versions_path)
        if row.get("tool") != "barrnap"
    ]
    rows.append(
        {
            "tool": "barrnap",
            "expected_version": BARRNAP_VERSION,
            "detected_version": detected_version,
            "status": status,
            "command": " ".join(command),
        }
    )
    _write_tsv(rows, tool_versions_path, TOOL_VERSION_FIELDS)


def _select_hits(hits: list[BarrnapHit], policy: str) -> list[tuple[int, BarrnapHit]]:
    indexed = list(enumerate(hits))
    if not hits:
        return []
    if policy == "all":
        return indexed
    if policy == "fail" and len(hits) > 1:
        return []
    if policy == "best":
        numeric_hits = [(index, hit) for index, hit in indexed if hit.numeric_score is not None]
        if numeric_hits:
            return [max(numeric_hits, key=lambda item: (item[1].numeric_score or 0.0, item[1].length))]
    return [max(indexed, key=lambda item: item[1].length)]


def _extract_hit(sequence: str, hit: BarrnapHit, flank: int) -> dict[str, Any]:
    seq_len = len(sequence)
    cut_start = max(1, hit.start - flank)
    cut_end = min(seq_len, hit.end + flank)
    extracted = sequence[cut_start - 1 : cut_end]
    if hit.strand == "-":
        extracted = _reverse_complement(extracted)
    return {"start": cut_start, "end": cut_end, "sequence": extracted}


def _barrnap_summary_row(
    internal_seq_id: str,
    original_seq_id: str,
    dataset: str,
    domain: str,
    input_length: str,
    status: str,
    extracted_seq_id: str = "",
    hit_count: int = 0,
    selected_hit_index: int | str = "",
    hit: BarrnapHit | None = None,
    extracted_start: int | str = "",
    extracted_end: int | str = "",
    extracted_length: int | str = "",
    reverse_complemented: bool = False,
    warning: str = "",
) -> dict[str, str]:
    return {
        "internal_seq_id": internal_seq_id,
        "extracted_seq_id": extracted_seq_id,
        "original_seq_id": original_seq_id,
        "dataset": dataset,
        "domain": domain,
        "input_length": str(input_length),
        "hit_count": str(hit_count),
        "selected_hit_index": str(selected_hit_index),
        "barrnap_type": hit.feature_type if hit else "",
        "barrnap_start": str(hit.start) if hit else "",
        "barrnap_end": str(hit.end) if hit else "",
        "barrnap_strand": hit.strand if hit else "",
        "barrnap_score": hit.score if hit else "",
        "extracted_start": str(extracted_start),
        "extracted_end": str(extracted_end),
        "extracted_length": str(extracted_length),
        "reverse_complemented": _bool_text(reverse_complemented),
        "status": status,
        "warning": warning,
    }


def _hits_by_seq_id(hits: Iterable[BarrnapHit]) -> dict[str, list[BarrnapHit]]:
    by_seq_id: dict[str, list[BarrnapHit]] = {}
    for hit in hits:
        by_seq_id.setdefault(hit.seq_id, []).append(hit)
    return by_seq_id


def _is_ssu_rrna_hit(hit: BarrnapHit) -> bool:
    haystack = " ".join(
        [
            hit.feature_type,
            hit.attributes.get("product", ""),
            hit.attributes.get("Name", ""),
            hit.attributes.get("name", ""),
            hit.attributes.get("note", ""),
        ]
    ).lower()
    return (
        "rrna" in haystack
        and (
            "16s" in haystack
            or "ssu" in haystack
            or "small subunit" in haystack
            or hit.feature_type.lower() == "rrna"
        )
    )


def _parse_gff3_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for chunk in text.split(";"):
        if not chunk:
            continue
        if "=" in chunk:
            key, value = chunk.split("=", maxsplit=1)
        elif " " in chunk:
            key, value = chunk.split(" ", maxsplit=1)
        else:
            key, value = chunk, ""
        attributes[key.strip()] = value.strip()
    return attributes


def _default_barrnap_kingdom(domain: str) -> str:
    return {"Archaea": "arc", "Bacteria": "bac"}[domain]


def _validate_domain(domain: str) -> str:
    if domain not in {"Archaea", "Bacteria"}:
        raise ValueError("Domain must be Archaea or Bacteria.")
    return domain


def _validate_policy(policy: str) -> str:
    if policy not in {"longest", "best", "all", "fail"}:
        raise ValueError("multi-rrna-policy must be one of: longest, best, all, fail.")
    return policy


def _validate_dataset_name(name: str) -> str:
    normalized = name.strip()
    if not normalized:
        raise ValueError("Dataset name must not be empty.")
    return normalized


def _safe_dataset_dir_name(name: str) -> str:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return safe_name or "dataset"


def _file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ATGC", "TACG"))[::-1]


def _append_warning(existing: str, warning: str) -> str:
    if not existing:
        return warning
    if warning in existing.split(","):
        return existing
    return f"{existing},{warning}"


def _count_status(rows: Iterable[dict[str, str]], status: str) -> int:
    return sum(1 for row in rows if row.get("status") == status)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))


def _write_tsv(rows: Iterable[dict[str, Any]], path: Path, fieldnames: list[str]) -> None:
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


def _merge_fieldnames(preferred: list[str], rows: Iterable[dict[str, str]]) -> list[str]:
    fieldnames = list(preferred)
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    return fieldnames


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
