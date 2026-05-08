"""Export helpers for downstream classifier formats."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from autotax2.export.dada2 import DADA2_SPECIES_FORMAT, DADA2_TOGENUS_FORMAT
from autotax2.export.formats import (
    format_taxonomy_dada2_genus,
    format_taxonomy_dada2_species,
    format_taxonomy_qiime2,
    format_taxonomy_sintax,
    strip_rank_prefix,
    validate_taxonomy_7rank,
)
from autotax2.export.qiime2 import QIIME2_FORMATS, QIIME2_REFERENCE_SEQUENCES_FORMAT, QIIME2_TAXONOMY_FORMAT
from autotax2.export.sintax import SINTAX_FORMAT, format_sintax_header
from autotax2.io import FastaRecord, read_fasta, write_fasta
from autotax2.registry import sequence_md5


REQUIRED_EXPORT_FORMATS = (
    SINTAX_FORMAT,
    *QIIME2_FORMATS,
    DADA2_TOGENUS_FORMAT,
    DADA2_SPECIES_FORMAT,
)
EXPORT_COMMANDS = ("sintax", "qiime2", "dada2", "all")
MANIFEST_FIELDS = [
    "format",
    "path",
    "gzip",
    "representatives_only",
    "records_exported",
    "unique_md5_exported",
    "created_at",
    "autotax2_version",
    "registry_version_or_build_id",
]


@dataclass(frozen=True)
class ExportRecord:
    """One sequence/taxonomy pair ready for classifier export."""

    seq_id: str
    sequence: str
    taxon_id: str
    taxonomy_7rank: tuple[str, ...]
    source: str
    sequence_md5: str
    is_representative: bool
    is_duplicate_sequence: bool


@dataclass(frozen=True)
class ExportSummary:
    """Summary returned by :func:`export_references`."""

    build: Path
    outdir: Path
    formats: tuple[str, ...]
    records_exported: int
    manifest_path: Path


def export_references(
    build: str | Path,
    format_name: str,
    outdir: str | Path | None = None,
    gzip_output: bool = True,
    representatives_only: bool = True,
    output_prefix: str = "autotax2",
    force: bool = False,
) -> ExportSummary:
    """Export classifier reference files from an autotax2 build registry."""
    normalized_format = format_name.strip().lower()
    if normalized_format not in EXPORT_COMMANDS:
        raise ValueError(f"Export format must be one of: {', '.join(EXPORT_COMMANDS)}")

    build_dir = Path(build)
    export_dir = Path(outdir) if outdir is not None else build_dir / "export"
    export_dir.mkdir(parents=True, exist_ok=True)
    records = build_exportable_records(
        build_dir,
        representatives_only=representatives_only,
    )
    formats = ("sintax", "qiime2", "dada2") if normalized_format == "all" else (normalized_format,)
    manifest_rows: list[dict[str, str]] = []
    for fmt in formats:
        if fmt == "sintax":
            manifest_rows.extend(
                _export_sintax(
                    records=records,
                    outdir=export_dir,
                    gzip_output=gzip_output,
                    representatives_only=representatives_only,
                    output_prefix=output_prefix,
                    force=force,
                    build_dir=build_dir,
                )
            )
        elif fmt == "qiime2":
            manifest_rows.extend(
                _export_qiime2(
                    records=records,
                    outdir=export_dir,
                    gzip_output=gzip_output,
                    representatives_only=representatives_only,
                    force=force,
                    build_dir=build_dir,
                )
            )
        elif fmt == "dada2":
            manifest_rows.extend(
                _export_dada2(
                    records=records,
                    outdir=export_dir,
                    gzip_output=gzip_output,
                    representatives_only=representatives_only,
                    output_prefix=output_prefix,
                    force=force,
                    build_dir=build_dir,
                )
            )

    manifest_path = export_dir / "export_manifest.tsv"
    _ensure_writable(manifest_path, force=force)
    _write_tsv(manifest_rows, manifest_path, MANIFEST_FIELDS)
    return ExportSummary(
        build=build_dir,
        outdir=export_dir,
        formats=formats,
        records_exported=len(records),
        manifest_path=manifest_path,
    )


def build_exportable_records(
    build_dir: str | Path,
    representatives_only: bool = True,
) -> list[ExportRecord]:
    """Build de-duplicated exportable records from the current registry layout."""
    build_path = Path(build_dir)
    registry_dir = build_path / "registry"
    taxon_by_id = _taxon_by_id(registry_dir)
    sequences = _sequence_lookup(build_path)
    representative_rows = _active_representative_rows(registry_dir)
    sequence_rows = _read_tsv(registry_dir / "sequence_registry.tsv")
    sequence_rows_by_id = {
        _first_value(row, "seq_id", "internal_seq_id", "representative_seq_id"): row
        for row in sequence_rows
        if _first_value(row, "seq_id", "internal_seq_id", "representative_seq_id")
    }

    candidates: list[tuple[str, dict[str, str], bool]] = []
    if representatives_only:
        for row in representative_rows:
            seq_id = row.get("representative_seq_id", "")
            if seq_id:
                candidates.append((seq_id, row, True))
        for row in sequence_rows:
            if row.get("taxonomy_7rank") and _first_value(row, "seq_id", "internal_seq_id") not in {
                candidate[0] for candidate in candidates
            }:
                candidates.append((_first_value(row, "seq_id", "internal_seq_id"), row, True))
    else:
        for row in sequence_rows:
            seq_id = _first_value(row, "seq_id", "internal_seq_id", "representative_seq_id")
            if seq_id:
                candidates.append((seq_id, row, False))
        existing_ids = {seq_id for seq_id, _, _ in candidates}
        for row in representative_rows:
            seq_id = row.get("representative_seq_id", "")
            if seq_id and seq_id not in existing_ids:
                candidates.append((seq_id, row, True))

    records: list[ExportRecord] = []
    seen_md5: set[str] = set()
    for seq_id, candidate, is_representative in candidates:
        sequence = sequences.get(seq_id, "")
        if not sequence:
            raise ValueError(f"Missing sequence for export record: {seq_id}")
        digest = candidate.get("sequence_md5") or sequence_rows_by_id.get(seq_id, {}).get("sequence_md5", "")
        digest = digest or sequence_md5(sequence)
        if digest in seen_md5:
            continue

        taxon_id = _first_value(candidate, "taxon_id", "assigned_taxon_id", "species_taxon_id")
        taxonomy = candidate.get("taxonomy_7rank") or sequence_rows_by_id.get(seq_id, {}).get("taxonomy_7rank", "")
        if taxon_id:
            if not _taxon_is_active(taxon_id, taxon_by_id):
                continue
            taxonomy_values = _taxonomy_for_taxon(taxon_id, taxon_by_id)
        elif taxonomy:
            taxonomy_values = validate_taxonomy_7rank(taxonomy)
        else:
            raise ValueError(f"Missing 7-rank taxonomy for export record: {seq_id}")

        taxonomy_7rank = validate_taxonomy_7rank(taxonomy_values)
        source = (
            candidate.get("source")
            or sequence_rows_by_id.get(seq_id, {}).get("source", "")
            or taxon_by_id.get(taxon_id, {}).get("source", "")
        )
        records.append(
            ExportRecord(
                seq_id=seq_id,
                sequence=sequence,
                taxon_id=taxon_id,
                taxonomy_7rank=taxonomy_7rank,
                source=source,
                sequence_md5=digest,
                is_representative=is_representative,
                is_duplicate_sequence=_is_true(candidate.get("is_duplicate_sequence", "false")),
            )
        )
        seen_md5.add(digest)

    return records


def _export_sintax(
    records: list[ExportRecord],
    outdir: Path,
    gzip_output: bool,
    representatives_only: bool,
    output_prefix: str,
    force: bool,
    build_dir: Path,
) -> list[dict[str, str]]:
    sintax_dir = outdir / "sintax"
    sintax_dir.mkdir(parents=True, exist_ok=True)
    path = sintax_dir / f"{output_prefix}.sintax.fa{'.gz' if gzip_output else ''}"
    _ensure_writable(path, force=force)
    write_fasta(
        [
            FastaRecord(
                seq_id=record.seq_id,
                header=format_sintax_header(record.seq_id, record.taxonomy_7rank),
                sequence=record.sequence,
            )
            for record in records
        ],
        path,
    )
    return [
        _manifest_row(
            fmt="sintax",
            path=path,
            gzip_output=gzip_output,
            representatives_only=representatives_only,
            records=records,
            build_dir=build_dir,
        )
    ]


def _export_qiime2(
    records: list[ExportRecord],
    outdir: Path,
    gzip_output: bool,
    representatives_only: bool,
    force: bool,
    build_dir: Path,
) -> list[dict[str, str]]:
    qiime2_dir = outdir / "qiime2"
    qiime2_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = qiime2_dir / f"reference_sequences.fasta{'.gz' if gzip_output else ''}"
    taxonomy_path = qiime2_dir / "reference_taxonomy.tsv"
    _ensure_writable(fasta_path, force=force)
    _ensure_writable(taxonomy_path, force=force)
    write_fasta(
        [
            FastaRecord(seq_id=record.seq_id, header=record.seq_id, sequence=record.sequence)
            for record in records
        ],
        fasta_path,
    )
    _write_tsv(
        [
            {
                "Feature ID": record.seq_id,
                "Taxon": format_taxonomy_qiime2(record.taxonomy_7rank),
            }
            for record in records
        ],
        taxonomy_path,
        ["Feature ID", "Taxon"],
    )
    return [
        _manifest_row(
            fmt="qiime2-reference-sequences",
            path=fasta_path,
            gzip_output=gzip_output,
            representatives_only=representatives_only,
            records=records,
            build_dir=build_dir,
        ),
        _manifest_row(
            fmt="qiime2-taxonomy",
            path=taxonomy_path,
            gzip_output=False,
            representatives_only=representatives_only,
            records=records,
            build_dir=build_dir,
        ),
    ]


def _export_dada2(
    records: list[ExportRecord],
    outdir: Path,
    gzip_output: bool,
    representatives_only: bool,
    output_prefix: str,
    force: bool,
    build_dir: Path,
) -> list[dict[str, str]]:
    dada2_dir = outdir / "dada2"
    dada2_dir.mkdir(parents=True, exist_ok=True)
    genus_path = dada2_dir / f"{output_prefix}_toGenus_trainset.fa{'.gz' if gzip_output else ''}"
    species_path = dada2_dir / f"{output_prefix}_assignSpecies.fa{'.gz' if gzip_output else ''}"
    _ensure_writable(genus_path, force=force)
    _ensure_writable(species_path, force=force)
    write_fasta(
        [
            FastaRecord(
                seq_id=record.seq_id,
                header=f"{record.seq_id} {format_taxonomy_dada2_genus(record.taxonomy_7rank)}",
                sequence=record.sequence,
            )
            for record in records
        ],
        genus_path,
    )
    write_fasta(
        [
            FastaRecord(
                seq_id=record.seq_id,
                header=f"{record.seq_id} {format_taxonomy_dada2_species(record.taxonomy_7rank)}",
                sequence=record.sequence,
            )
            for record in records
        ],
        species_path,
    )
    return [
        _manifest_row(
            fmt="dada2-togenus",
            path=genus_path,
            gzip_output=gzip_output,
            representatives_only=representatives_only,
            records=records,
            build_dir=build_dir,
        ),
        _manifest_row(
            fmt="dada2-assignspecies",
            path=species_path,
            gzip_output=gzip_output,
            representatives_only=representatives_only,
            records=records,
            build_dir=build_dir,
        ),
    ]


def _sequence_lookup(build_dir: Path) -> dict[str, str]:
    paths = [
        build_dir / "registry" / "current_representatives.fa",
        build_dir / "silva" / "silva_named_backbone.fa",
        build_dir / "silva" / "silva_unresolved.resolved.fa",
        build_dir / "silva" / "silva_unresolved.fa",
    ]
    datasets_dir = build_dir / "datasets"
    if datasets_dir.exists():
        for dataset_dir in sorted(path for path in datasets_dir.iterdir() if path.is_dir()):
            paths.extend(
                [
                    dataset_dir / "sina.oriented.fa",
                    dataset_dir / "barrnap.extracted.fa",
                    dataset_dir / "input.normalized.fa",
                ]
            )
    sequences: dict[str, str] = {}
    for path in paths:
        if not path.exists():
            continue
        for record in read_fasta(path):
            sequences.setdefault(record.seq_id, record.sequence)
    return sequences


def _taxon_by_id(registry_dir: Path) -> dict[str, dict[str, str]]:
    return {
        row.get("taxon_id", ""): row
        for row in _read_tsv(registry_dir / "taxon_nodes.tsv")
        if row.get("taxon_id")
    }


def _active_representative_rows(registry_dir: Path) -> list[dict[str, str]]:
    rows = []
    for row in _read_tsv(registry_dir / "representative_registry.tsv"):
        if row.get("status", "active") in {"", "active"}:
            rows.append(row)
    return rows


def _taxonomy_for_taxon(taxon_id: str, taxon_by_id: dict[str, dict[str, str]]) -> tuple[str, ...]:
    labels: list[str] = []
    seen: set[str] = set()
    current = taxon_id
    while current and current not in seen:
        seen.add(current)
        row = taxon_by_id.get(current, {})
        if not row:
            break
        labels.append(row.get("name", current))
        current = row.get("parent_taxon_id", "")
    if len(labels) != 7:
        raise ValueError(f"Taxon {taxon_id} does not resolve to a complete 7-rank taxonomy.")
    return tuple(reversed(labels))


def _taxon_is_active(taxon_id: str, taxon_by_id: dict[str, dict[str, str]]) -> bool:
    row = taxon_by_id.get(taxon_id, {})
    return row.get("status", "active") not in {"deprecated", "superseded"}


def _manifest_row(
    fmt: str,
    path: Path,
    gzip_output: bool,
    representatives_only: bool,
    records: list[ExportRecord],
    build_dir: Path,
) -> dict[str, str]:
    return {
        "format": fmt,
        "path": str(path).replace("\\", "/"),
        "gzip": _bool_text(gzip_output),
        "representatives_only": _bool_text(representatives_only),
        "records_exported": str(len(records)),
        "unique_md5_exported": str(len({record.sequence_md5 for record in records})),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "autotax2_version": _autotax2_version(),
        "registry_version_or_build_id": build_dir.name or "build",
    }


def _autotax2_version() -> str:
    try:
        return version("autotax2")
    except PackageNotFoundError:
        return "0.1.0"


def _ensure_writable(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"Export output already exists; use --force to overwrite: {path}")


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


def _first_value(row: dict[str, str], *fields: str) -> str:
    for field in fields:
        value = row.get(field, "")
        if value:
            return value
    return ""


def _is_true(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "y"}


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


__all__ = [
    "DADA2_SPECIES_FORMAT",
    "DADA2_TOGENUS_FORMAT",
    "EXPORT_COMMANDS",
    "ExportRecord",
    "ExportSummary",
    "QIIME2_FORMATS",
    "REQUIRED_EXPORT_FORMATS",
    "SINTAX_FORMAT",
    "build_exportable_records",
    "export_references",
    "format_sintax_header",
    "format_taxonomy_dada2_genus",
    "format_taxonomy_dada2_species",
    "format_taxonomy_qiime2",
    "format_taxonomy_sintax",
    "strip_rank_prefix",
    "validate_taxonomy_7rank",
]
