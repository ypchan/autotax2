"""Input/output helpers for autotax2."""

from __future__ import annotations

import csv
import gzip
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO


FASTA_WRAP_WIDTH = 80


@dataclass(frozen=True)
class FastaRecord:
    """A normalized FASTA record."""

    seq_id: str
    header: str
    sequence: str


def read_text(path: str | Path) -> str:
    """Read UTF-8 text from a file."""
    return Path(path).read_text(encoding="utf-8")


def write_text(path: str | Path, content: str) -> None:
    """Write UTF-8 text to a file."""
    Path(path).write_text(content, encoding="utf-8")


def normalize_sequence(sequence: str) -> str:
    """Normalize a biological sequence for storage and exact MD5 comparison."""
    return re.sub(r"\s+", "", sequence).upper().replace("U", "T")


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read plain or gzipped FASTA into normalized records."""
    records: list[FastaRecord] = []
    current_header: str | None = None
    current_chunks: list[str] = []

    with _open_text(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append(_build_fasta_record(current_header, current_chunks))
                current_header = line[1:].strip()
                if not current_header:
                    raise ValueError(f"Empty FASTA header at line {line_number}.")
                current_chunks = []
                continue
            if current_header is None:
                raise ValueError(f"FASTA sequence encountered before header at line {line_number}.")
            current_chunks.append(line)

    if current_header is not None:
        records.append(_build_fasta_record(current_header, current_chunks))

    return records


def write_fasta(records: Iterable[FastaRecord], path: str | Path) -> None:
    """Write FASTA records, gzipping output when ``path`` ends with ``.gz``."""
    with _open_text(path, "wt") as handle:
        for record in records:
            header = record.header.strip() or record.seq_id
            print(f">{header}", file=handle)
            sequence = normalize_sequence(record.sequence)
            for start in range(0, len(sequence), FASTA_WRAP_WIDTH):
                print(sequence[start : start + FASTA_WRAP_WIDTH], file=handle)


def iter_fasta_records(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield FASTA records as ``(identifier, sequence)`` pairs."""
    for record in read_fasta(path):
        yield record.seq_id, record.sequence


def write_sequence_id_map(records: Iterable[Any], path: str | Path) -> None:
    """Write ``sequence_id_map.tsv`` records."""
    _write_tsv(
        records,
        path,
        [
            "internal_seq_id",
            "original_seq_id",
            "original_header",
            "dataset",
            "prefix",
            "sequence_md5",
            "sequence_length",
        ],
    )


def write_unique_sequence_registry(records: Iterable[Any], path: str | Path) -> None:
    """Write ``unique_sequence_registry.tsv`` records."""
    _write_tsv(
        records,
        path,
        [
            "unique_seq_id",
            "sequence_md5",
            "representative_internal_seq_id",
            "sequence",
            "length",
            "first_seen_dataset",
        ],
    )


def write_sequence_membership(records: Iterable[Any], path: str | Path) -> None:
    """Write ``sequence_membership.tsv`` records."""
    _write_tsv(
        records,
        path,
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


def _build_fasta_record(header: str, chunks: list[str]) -> FastaRecord:
    seq_id = header.split(maxsplit=1)[0]
    return FastaRecord(
        seq_id=seq_id,
        header=header,
        sequence=normalize_sequence("".join(chunks)),
    )


def _open_text(path: str | Path, mode: str) -> TextIO:
    fasta_path = Path(path)
    if fasta_path.suffix == ".gz":
        return gzip.open(fasta_path, mode, encoding="utf-8", newline="")
    return fasta_path.open(mode, encoding="utf-8", newline="")


def _write_tsv(records: Iterable[Any], path: str | Path, fieldnames: list[str]) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=chr(9),
            lineterminator=chr(10),
            extrasaction="ignore",
        )
        writer.writeheader()
        for record in records:
            writer.writerow({field: _field_value(record, field) for field in fieldnames})


def _field_value(record: Any, field: str) -> Any:
    if isinstance(record, dict):
        return record[field]
    return getattr(record, field)
