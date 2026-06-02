#!/usr/bin/env python3
"""Build an aligned FASTA with SINA/ARB metadata fields, with progress output.

Example:

    python3 script/build_sina_metadata_fasta_parallel.py \
      --input-fasta gtdbr232_type_or_complete_genome.ssu.sina.align.fa \
      --metadata-tsv gtdb_type_or_complete_genome.ssu.metadata.tsv \
      --output-fasta gtdbr232_type_or_complete_genome.ssu.sina.align_withmetadata.fa \
      --workers 8

The output can be converted to ARB with:

    sina -i aligned_withmetadata.fa --prealigned --outtype arb -o gtdb_ssu.arb
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import sys
import time
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

METADATA_BY_ID: dict[str, dict[str, str]] = {}
CONFIG: dict[str, object] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add tax_gtdb/tax_slv SINA metadata to an already aligned FASTA."
    )
    parser.add_argument("-i", "--input-fasta", required=True, type=Path)
    parser.add_argument("-m", "--metadata-tsv", required=True, type=Path)
    parser.add_argument("-o", "--output-fasta", required=True, type=Path)
    parser.add_argument("--id-column", default="arb_id")
    parser.add_argument("--taxonomy-column", default="gtdb_taxonomy")
    parser.add_argument("--tax-gtdb-field", default="tax_gtdb")
    parser.add_argument("--tax-slv-field", default="tax_slv")
    parser.add_argument("--extra-columns", default=",".join(DEFAULT_EXTRA_COLUMNS))
    parser.add_argument("--wrap", type=int, default=80)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=5000)
    parser.add_argument("--progress-every", type=int, default=10000)
    parser.add_argument("--no-count", action="store_true")
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument("--strict-extra-columns", action="store_true")
    return parser.parse_args()


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def sanitize_field_name(name: str) -> str:
    """Return a conservative metadata field name for SINA/ARB headers."""
    cleaned = "".join(char if char.isalnum() or char == "_" else "_" for char in name.strip())
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Invalid metadata field name: {name!r}")
    if cleaned[0].isdigit():
        cleaned = f"field_{cleaned}"
    return cleaned


def clean_value(value: str | None) -> str:
    """Make a metadata value safe for a single FASTA header line."""
    if value is None:
        return ""
    return " ".join(str(value).strip().replace("]", ")").split())


def read_metadata(path: Path, id_column: str) -> tuple[dict[str, dict[str, str]], list[str]]:
    started = time.time()
    log(f"[metadata] reading {path}")
    rows: dict[str, dict[str, str]] = {}
    duplicates: list[str] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Metadata file has no header: {path}")
        if id_column not in reader.fieldnames:
            raise ValueError(f"ID column {id_column!r} not found in {path}")
        fieldnames = list(reader.fieldnames)
        for index, row in enumerate(reader, start=1):
            seq_id = clean_value(row.get(id_column))
            if not seq_id:
                continue
            if seq_id in rows:
                duplicates.append(seq_id)
                continue
            rows[seq_id] = row
            if index % 1_000_000 == 0:
                log(f"[metadata] loaded {index:,} rows")
    if duplicates:
        raise ValueError(f"Duplicate metadata IDs found, examples: {', '.join(duplicates[:10])}")
    log(f"[metadata] loaded {len(rows):,} unique rows in {time.time() - started:.1f}s")
    return rows, fieldnames


def count_fasta_records(path: Path) -> int:
    started = time.time()
    log(f"[fasta] counting records in {path}")
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
                if count % 1_000_000 == 0:
                    log(f"[fasta] counted {count:,} records")
    log(f"[fasta] counted {count:,} records in {time.time() - started:.1f}s")
    return count


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


def iter_batches(records: Iterable[tuple[str, str]], batch_size: int) -> Iterable[list[tuple[str, str]]]:
    batch: list[tuple[str, str]] = []
    for record in records:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def init_worker(metadata: dict[str, dict[str, str]], config: dict[str, object]) -> None:
    global METADATA_BY_ID, CONFIG
    METADATA_BY_ID = metadata
    CONFIG = config


def wrap_sequence(sequence: str, wrap: int) -> str:
    if wrap <= 0:
        return sequence + "\n"
    return "".join(sequence[start : start + wrap] + "\n" for start in range(0, len(sequence), wrap))


def build_header(seq_id: str, row: dict[str, str]) -> str:
    taxonomy = clean_value(row.get(str(CONFIG["taxonomy_column"])))
    fields = [
        (str(CONFIG["tax_gtdb_field"]), taxonomy),
        (str(CONFIG["tax_slv_field"]), taxonomy),
    ]
    for column in list(CONFIG["extra_columns"]):
        value = clean_value(row.get(str(column)))
        if value:
            fields.append((str(column), value))
    metadata = " ".join(
        f"[{sanitize_field_name(key)}={value}]" for key, value in fields if value
    )
    return f"{seq_id} {metadata}".rstrip()


def format_batch(batch: list[tuple[str, str]]) -> tuple[int, int, list[str], list[str]]:
    wrap = int(CONFIG["wrap"])
    allow_missing = bool(CONFIG["allow_missing"])
    missing: list[str] = []
    output: list[str] = []
    for seq_id, sequence in batch:
        row = METADATA_BY_ID.get(seq_id)
        if row is None:
            missing.append(seq_id)
            if not allow_missing:
                continue
            header = seq_id
        else:
            header = build_header(seq_id, row)
        output.append(f">{header}\n{wrap_sequence(sequence, wrap)}")
    return len(batch), len(output), output, missing


def progress(done: int, total: int | None, started: float) -> None:
    elapsed = max(time.time() - started, 0.001)
    rate = done / elapsed
    if total:
        pct = done / total * 100
        eta = (total - done) / rate if rate > 0 else 0
        log(f"[progress] {done:,}/{total:,} ({pct:.1f}%), {rate:,.0f} rec/s, ETA {eta/60:.1f} min")
    else:
        log(f"[progress] {done:,} records, {rate:,.0f} rec/s")


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be >= 1")
    if args.batch_size < 1:
        raise ValueError("--batch-size must be >= 1")

    metadata, fieldnames = read_metadata(args.metadata_tsv, args.id_column)
    if args.taxonomy_column not in fieldnames:
        raise ValueError(f"Taxonomy column {args.taxonomy_column!r} not found")

    extra_columns = [item.strip() for item in args.extra_columns.split(",") if item.strip()]
    missing_extra = [column for column in extra_columns if column not in fieldnames]
    if missing_extra and args.strict_extra_columns:
        raise ValueError(f"Requested extra columns not found: {missing_extra}")
    if missing_extra:
        log(f"Warning: skipping absent extra columns: {', '.join(missing_extra)}")
        extra_columns = [column for column in extra_columns if column in fieldnames]

    total = None if args.no_count else count_fasta_records(args.input_fasta)
    config = {
        "taxonomy_column": args.taxonomy_column,
        "tax_gtdb_field": args.tax_gtdb_field,
        "tax_slv_field": args.tax_slv_field,
        "extra_columns": extra_columns,
        "wrap": args.wrap,
        "allow_missing": args.allow_missing,
    }

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    log(f"[run] workers={args.workers}, batch_size={args.batch_size}, output={args.output_fasta}")
    started = time.time()
    processed = 0
    written = 0
    next_report = args.progress_every
    missing_metadata: list[str] = []
    batch_iter = iter_batches(parse_fasta(args.input_fasta), args.batch_size)

    with args.output_fasta.open("w", encoding="utf-8", newline="") as out:
        if args.workers == 1:
            init_worker(metadata, config)
            result_iter = map(format_batch, batch_iter)
        else:
            pool = ProcessPoolExecutor(
                max_workers=args.workers,
                initializer=init_worker,
                initargs=(metadata, config),
            )
            result_iter = pool.map(format_batch, batch_iter)

        try:
            for batch_count, written_count, formatted, missing in result_iter:
                out.writelines(formatted)
                processed += batch_count
                written += written_count
                missing_metadata.extend(missing)
                if processed >= next_report:
                    progress(processed, total, started)
                    next_report += args.progress_every
        finally:
            if args.workers > 1:
                pool.shutdown()

    progress(processed, total, started)
    if missing_metadata:
        log(
            f"Warning: {len(missing_metadata):,} FASTA IDs had no metadata row. "
            f"Examples: {', '.join(missing_metadata[:10])}"
        )
    if missing_metadata and not args.allow_missing:
        log("Error: missing metadata rows were found.")
        log("Use --allow-missing only if this is intentional.")
        return 2
    log(f"[done] wrote {written:,} records to {args.output_fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
