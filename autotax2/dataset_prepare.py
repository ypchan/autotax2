"""Dataset preparation for externally processed SSU/16S sequences."""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from autotax2.io import FastaRecord, read_fasta, write_fasta
from autotax2.registry import (
    AssignedSequence,
    SequenceRegistry,
    assign_internal_ids,
)
from autotax2.validate import validate_dataset_prefix


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
PREPARE_SUMMARY_FIELDS = [
    "dataset",
    "prefix",
    "input_sequences",
    "normalized_sequences",
    "rejected_non_atgc",
    "rejected_short_ssu",
    "duplicate_md5_sequences",
    "final_prepared_sequences",
    "input_contract",
]
ATGC_RE = re.compile(r"^[ATGC]*$")


@dataclass(frozen=True)
class PrepareDatasetSummary:
    """Summary from dataset preparation."""

    dataset: str
    prefix: str
    dataset_dir: Path
    final_prepared_sequences: int


def prepare_dataset(
    build: str | Path,
    name: str,
    prefix: str,
    fasta: str | Path,
    domain: str,
    min_ssu_len_archaea: int = 900,
    min_ssu_len_bacteria: int = 1200,
    reject_non_atgc: bool = True,
) -> PrepareDatasetSummary:
    """Prepare a dataset whose input FASTA already contains SSU/16S sequences."""
    build_dir = Path(build)
    registry_dir = build_dir / "registry"
    registry_dir.mkdir(parents=True, exist_ok=True)
    normalized_name = _validate_dataset_name(name)
    normalized_prefix = validate_dataset_prefix(prefix)
    normalized_domain = _validate_domain(domain)
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

    source_records = read_fasta(input_fasta)
    assigned = assign_internal_ids(source_records, normalized_name, normalized_prefix)
    min_ssu_len = min_ssu_len_archaea if normalized_domain == "Archaea" else min_ssu_len_bacteria
    sequence_map_rows: list[dict[str, str]] = []
    normalized_records: list[FastaRecord] = []
    final_assigned: list[AssignedSequence] = []
    final_records: list[FastaRecord] = []
    rejected_non_atgc = 0
    rejected_short_ssu = 0

    for input_order, (source_record, assigned_record) in enumerate(
        zip(source_records, assigned, strict=True),
        start=1,
    ):
        rejected = False
        reject_reason = ""
        if reject_non_atgc and not ATGC_RE.fullmatch(assigned_record.sequence):
            rejected = True
            reject_reason = "non_atgc"
            rejected_non_atgc += 1
        elif assigned_record.sequence_length < min_ssu_len:
            rejected = True
            reject_reason = "short_ssu"
            rejected_short_ssu += 1

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

        if reject_reason != "non_atgc":
            normalized_records.append(
                FastaRecord(
                    seq_id=assigned_record.internal_seq_id,
                    header=assigned_record.internal_seq_id,
                    sequence=assigned_record.sequence,
                )
            )
        if not rejected:
            final_assigned.append(assigned_record)
            final_records.append(
                FastaRecord(
                    seq_id=assigned_record.internal_seq_id,
                    header=assigned_record.internal_seq_id,
                    sequence=assigned_record.sequence,
                )
            )

    write_fasta(normalized_records, dataset_dir / "input.normalized.fa")
    write_fasta(final_records, dataset_dir / "prepared.ssu.fa")
    _write_tsv(sequence_map_rows, dataset_dir / "sequence_id_map.tsv", SEQUENCE_ID_MAP_FIELDS)

    sequence_registry = SequenceRegistry()
    memberships = sequence_registry.add_many(final_assigned)
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

    duplicate_md5_sequences = sum(1 for membership in memberships if membership.is_duplicate_sequence)
    prepare_summary = {
        "dataset": normalized_name,
        "prefix": normalized_prefix,
        "input_sequences": str(len(source_records)),
        "normalized_sequences": str(len(normalized_records)),
        "rejected_non_atgc": str(rejected_non_atgc),
        "rejected_short_ssu": str(rejected_short_ssu),
        "duplicate_md5_sequences": str(duplicate_md5_sequences),
        "final_prepared_sequences": str(len(final_records)),
        "input_contract": "externally_processed_ssu",
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
        final_prepared_sequences=len(final_records),
    )


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


def _validate_domain(domain: str) -> str:
    if domain not in {"Archaea", "Bacteria"}:
        raise ValueError("Domain must be Archaea or Bacteria.")
    return domain


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
