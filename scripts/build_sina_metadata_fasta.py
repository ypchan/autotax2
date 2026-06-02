#!/usr/bin/env python3
"""Build an aligned FASTA with SINA/ARB metadata fields.

Input:
  1. an already aligned FASTA file;
  2. a tab-delimited metadata table whose ID column matches the FASTA IDs.

Output:
  an aligned FASTA whose headers contain fields such as:

      [tax_gtdb=d__Bacteria;p__...;s__...]
      [tax_slv=d__Bacteria;p__...;s__...]

The output can be converted to ARB with:

    sina -i gtdb.with_metadata.aligned.fa --prealigned --outtype arb -o gtdb_ssu.arb
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable


DEFAULT_EXTRA_COLUMNS = (
    "original_id",
    "genome_id",
    "contig_id",
    "copy_no",
    "location",
    "ssu_len",
    "contig_len",
    "seq_len",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add tax_gtdb/tax_slv metadata to an already aligned FASTA."
    )
    parser.add_argument("-i", "--input-fasta", required=True, type=Path)
    parser.add_argument("-m", "--metadata-tsv", required=True, type=Path)
    parser.add_argument("-o", "--output-fasta", required=True, type=Path)
    parser.add_argument(
        "--id-column",
        default="arb_id",
        help="Metadata column matching FASTA IDs. Default: arb_id.",
    )
    parser.add_argument(
        "--taxonomy-column",
        default="gtdb_taxonomy",
        help="Metadata column containing full GTDB lineage. Default: gtdb_taxonomy.",
    )
    parser.add_argument(
        "--tax-gtdb-field",
        default="tax_gtdb",
        help="Header field used for GTDB taxonomy. Default: tax_gtdb.",
    )
    parser.add_argument(
        "--tax-slv-field",
        default="tax_slv",
        help="Header field used as SINA-compatible taxonomy mirror. Default: tax_slv.",
    )
    parser.add_argument(
        "--extra-columns",
        default=",".join(DEFAULT_EXTRA_COLUMNS),
        help="Comma-separated metadata columns to copy into headers. Use '' to disable.",
    )
    parser.add_argument("--wrap", type=int, default=80, help="Sequence line width. Default: 80.")
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Keep FASTA records without metadata instead of failing.",
    )
    parser.add_argument(
        "--strict-extra-columns",
        action="store_true",
        help="Fail if any requested extra metadata column is absent.",
    )
    return parser.parse_args()


def sanitize_field_name(name: str) -> str:
    """Return a conservative metadata field name for SINA/ARB headers."""
    cleaned = "".join(char if char.isalnum() or char == "_" else "_" for char in name.strip())
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid metadata field name: {name!r}")
    if cleaned[0].isdigit():
        cleaned = f"field_{cleaned}"
    return cleaned


def clean_field_value(value: str | None) -> str:
    """Make a metadata value safe for a single SINA-style header field."""
    if value is None:
        return ""
    value = str(value).strip().replace("]", ")")
    return " ".join(value.split())


def read_metadata(path: Path, id_column: str) -> tuple[dict[str, dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Metadata file has no header: {path}")
        if id_column not in reader.fieldnames:
            raise ValueError(f"ID column {id_column!r} not found in {path}")

        rows: dict[str, dict[str, str]] = {}
        duplicates: list[str] = []
        for row in reader:
            seq_id = clean_field_value(row.get(id_column))
            if not seq_id:
                continue
            if seq_id in rows:
                duplicates.append(seq_id)
                continue
            rows[seq_id] = row

    if duplicates:
        raise ValueError(f"Duplicate metadata IDs found, examples: {', '.join(duplicates[:10])}")
    return rows, list(reader.fieldnames)


def parse_fasta(path: Path) -> Iterable[tuple[str, str]]:
    """Yield (seq_id, sequence) from FASTA."""
    seq_id: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if seq_id is not None:
        yield seq_id, "".join(chunks)


def build_header(
    seq_id: str,
    row: dict[str, str],
    taxonomy_column: str,
    tax_gtdb_field: str,
    tax_slv_field: str,
    extra_columns: list[str],
) -> str:
    taxonomy = clean_field_value(row.get(taxonomy_column))
    fields = [
        (tax_gtdb_field, taxonomy),
        (tax_slv_field, taxonomy),
    ]
    for column in extra_columns:
        value = clean_field_value(row.get(column))
        if value:
            fields.append((column, value))
    metadata = " ".join(
        f"[{sanitize_field_name(key)}={value}]" for key, value in fields if value
    )
    return f"{seq_id} {metadata}".rstrip()


def write_sequence(handle, sequence: str, wrap: int) -> None:
    if wrap <= 0:
        handle.write(sequence + "\n")
        return
    for start in range(0, len(sequence), wrap):
        handle.write(sequence[start : start + wrap] + "\n")


def main() -> int:
    args = parse_args()
    metadata, fieldnames = read_metadata(args.metadata_tsv, args.id_column)

    if args.taxonomy_column not in fieldnames:
        raise ValueError(f"Taxonomy column {args.taxonomy_column!r} not found")

    extra_columns = [item.strip() for item in args.extra_columns.split(",") if item.strip()]
    missing_extra = [column for column in extra_columns if column not in fieldnames]
    if missing_extra and args.strict_extra_columns:
        raise ValueError(f"Requested extra columns not found: {missing_extra}")
    if missing_extra:
        print(
            f"Warning: skipping absent extra columns: {', '.join(missing_extra)}",
            file=sys.stderr,
        )
        extra_columns = [column for column in extra_columns if column in fieldnames]

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    missing_metadata: list[str] = []
    seen_fasta_ids: set[str] = set()
    written = 0

    with args.output_fasta.open("w", encoding="utf-8", newline="") as out:
        for seq_id, sequence in parse_fasta(args.input_fasta):
            seen_fasta_ids.add(seq_id)
            row = metadata.get(seq_id)
            if row is None:
                missing_metadata.append(seq_id)
                if not args.allow_missing:
                    continue
                header = seq_id
            else:
                header = build_header(
                    seq_id=seq_id,
                    row=row,
                    taxonomy_column=args.taxonomy_column,
                    tax_gtdb_field=args.tax_gtdb_field,
                    tax_slv_field=args.tax_slv_field,
                    extra_columns=extra_columns,
                )
            out.write(f">{header}\n")
            write_sequence(out, sequence, args.wrap)
            written += 1

    metadata_only = sorted(set(metadata) - seen_fasta_ids)
    if missing_metadata:
        print(
            f"Warning: {len(missing_metadata)} FASTA IDs had no metadata row. "
            f"Examples: {', '.join(missing_metadata[:10])}",
            file=sys.stderr,
        )
    if metadata_only:
        print(
            f"Warning: {len(metadata_only)} metadata IDs were not found in FASTA. "
            f"Examples: {', '.join(metadata_only[:10])}",
            file=sys.stderr,
        )
    if missing_metadata and not args.allow_missing:
        print("Error: missing metadata rows were found.", file=sys.stderr)
        print("Use --allow-missing only if this is intentional.", file=sys.stderr)
        return 2

    print(f"Wrote {written} records to {args.output_fasta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
