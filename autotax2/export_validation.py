"""Self-checks for exported classifier reference files."""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


EXPORT_VALIDATION_FIELDS = ["level", "check", "status", "message", "path"]


@dataclass(frozen=True)
class ExportValidationFinding:
    """One export format validation finding."""

    level: str
    check: str
    status: str
    message: str
    path: str


def validate_export_dir(
    export_dir: str | Path,
    formats: Iterable[str] | None = None,
) -> list[ExportValidationFinding]:
    """Validate existing export files, optionally limited to selected formats."""
    root = Path(export_dir)
    selected = {fmt.lower() for fmt in formats} if formats is not None else set()
    findings: list[ExportValidationFinding] = []

    if "sintax" in selected or (not selected and (root / "sintax").exists()):
        _check_sintax(root, findings, require=bool(selected))
    if "qiime2" in selected or (not selected and (root / "qiime2").exists()):
        _check_qiime2(root, findings, require=bool(selected))
    if "dada2" in selected or (not selected and (root / "dada2").exists()):
        _check_dada2(root, findings, require=bool(selected))

    if not findings:
        findings.append(_finding("info", "export_formats", "ok", "No export files were selected for validation.", root))
    return findings


def write_export_validation_report(
    findings: Iterable[ExportValidationFinding],
    path: str | Path,
) -> Path:
    """Write export validation findings as TSV."""
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
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
        report_path,
        EXPORT_VALIDATION_FIELDS,
    )
    return report_path


def _check_sintax(root: Path, findings: list[ExportValidationFinding], require: bool) -> None:
    directory = root / "sintax"
    paths = _fasta_files(directory)
    if require and not paths:
        findings.append(_finding("error", "export_sintax", "error", "SINTAX export file is missing.", directory))
        return
    checked_headers = 0
    for path in paths:
        headers = _fasta_headers(path)
        if not headers:
            findings.append(_finding("error", "export_sintax", "error", "SINTAX FASTA has no records.", path))
            continue
        for header in headers:
            checked_headers += 1
            if ";tax=" not in header:
                findings.append(_finding("error", "export_sintax", "error", f"SINTAX header lacks ;tax=: {header}", path))
                continue
            tax = header.split(";tax=", maxsplit=1)[1].rstrip(";")
            for token in ("d:", "p:", "c:", "o:", "f:", "g:", "s:"):
                if token not in tax:
                    findings.append(_finding("error", "export_sintax", "error", f"SINTAX tax string lacks {token}: {header}", path))
            if "__" in tax:
                findings.append(_finding("error", "export_sintax", "error", f"SINTAX tax values contain rank prefixes: {header}", path))
    if checked_headers:
        findings.append(_finding("info", "export_sintax", "ok", f"Validated {checked_headers} SINTAX headers.", directory))


def _check_qiime2(root: Path, findings: list[ExportValidationFinding], require: bool) -> None:
    directory = root / "qiime2"
    fasta_paths = _fasta_files(directory)
    taxonomy_path = directory / "reference_taxonomy.tsv"
    if require and not fasta_paths:
        findings.append(_finding("error", "export_qiime2", "error", "QIIME2 reference FASTA is missing.", directory))
    if require and not taxonomy_path.exists():
        findings.append(_finding("error", "export_qiime2", "error", "QIIME2 taxonomy TSV is missing.", taxonomy_path))
        return
    if not taxonomy_path.exists():
        return

    rows = _read_tsv(taxonomy_path)
    if not rows:
        findings.append(_finding("error", "export_qiime2", "error", "QIIME2 taxonomy file is empty.", taxonomy_path))
        return
    taxonomy_ids: set[str] = set()
    for row in rows:
        if "Feature ID" not in row or "Taxon" not in row:
            findings.append(_finding("error", "export_qiime2", "error", "QIIME2 taxonomy TSV lacks Feature ID and Taxon columns.", taxonomy_path))
            break
        feature_id = row.get("Feature ID", "")
        taxon = row.get("Taxon", "")
        taxonomy_ids.add(feature_id)
        if "d__" not in taxon:
            findings.append(_finding("error", "export_qiime2", "error", "QIIME2 taxonomy lacks d__ prefix.", taxonomy_path))
        if len([part for part in taxon.split(";") if part.strip()]) != 7:
            findings.append(_finding("error", "export_qiime2", "error", f"QIIME2 taxonomy is not 7-rank: {feature_id}", taxonomy_path))

    fasta_ids = set()
    for path in fasta_paths:
        fasta_ids.update(_fasta_header_ids(path))
    if fasta_ids and taxonomy_ids and fasta_ids != taxonomy_ids:
        missing_tax = sorted(fasta_ids - taxonomy_ids)
        missing_seq = sorted(taxonomy_ids - fasta_ids)
        if missing_tax:
            findings.append(_finding("error", "export_qiime2", "error", f"QIIME2 FASTA IDs lack taxonomy rows: {','.join(missing_tax[:5])}", taxonomy_path))
        if missing_seq:
            findings.append(_finding("error", "export_qiime2", "error", f"QIIME2 taxonomy IDs lack FASTA records: {','.join(missing_seq[:5])}", taxonomy_path))
    findings.append(_finding("info", "export_qiime2", "ok", f"Validated {len(rows)} QIIME2 taxonomy rows.", taxonomy_path))


def _check_dada2(root: Path, findings: list[ExportValidationFinding], require: bool) -> None:
    directory = root / "dada2"
    genus_paths = list(directory.glob("*toGenus*.fa*")) if directory.exists() else []
    species_paths = list(directory.glob("*assignSpecies*.fa*")) if directory.exists() else []
    if require and not genus_paths:
        findings.append(_finding("error", "export_dada2", "error", "DADA2 toGenus export is missing.", directory))
    if require and not species_paths:
        findings.append(_finding("error", "export_dada2", "error", "DADA2 assignSpecies export is missing.", directory))

    genus_headers = 0
    for path in genus_paths:
        headers = _fasta_headers(path)
        if not headers:
            findings.append(_finding("error", "export_dada2", "error", "DADA2 toGenus FASTA has no records.", path))
        for header in headers:
            genus_headers += 1
            taxonomy = header.split(maxsplit=1)[1] if len(header.split(maxsplit=1)) > 1 else ""
            if taxonomy.count(";") < 5:
                findings.append(_finding("error", "export_dada2", "error", f"DADA2 toGenus header does not contain taxonomy to genus: {header}", path))
            if "g__" in header or "s__" in header:
                findings.append(_finding("error", "export_dada2", "error", f"DADA2 toGenus header contains rank prefixes: {header}", path))
    species_headers = 0
    for path in species_paths:
        headers = _fasta_headers(path)
        if not headers:
            findings.append(_finding("error", "export_dada2", "error", "DADA2 assignSpecies FASTA has no records.", path))
        for header in headers:
            species_headers += 1
            parts = header.split()
            if len(parts) < 3:
                findings.append(_finding("error", "export_dada2", "error", f"DADA2 assignSpecies header lacks seq_id genus species: {header}", path))
            if ";" in header:
                findings.append(_finding("error", "export_dada2", "error", f"DADA2 assignSpecies header contains semicolon taxonomy: {header}", path))
            if "g__" in header or "s__" in header:
                findings.append(_finding("error", "export_dada2", "error", f"DADA2 assignSpecies header contains rank prefixes: {header}", path))
    if genus_headers or species_headers:
        findings.append(_finding("info", "export_dada2", "ok", f"Validated DADA2 headers: toGenus={genus_headers}, assignSpecies={species_headers}.", directory))


def _fasta_files(directory: Path) -> list[Path]:
    if not directory.exists():
        return []
    paths: list[Path] = []
    for pattern in ("*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz"):
        paths.extend(directory.glob(pattern))
    return sorted(set(paths))


def _fasta_headers(path: Path) -> list[str]:
    headers: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:  # type: ignore[arg-type]
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                headers.append(line[1:])
    return headers


def _fasta_header_ids(path: Path) -> list[str]:
    return [header.split()[0].split(";tax=", maxsplit=1)[0] for header in _fasta_headers(path)]


def _read_tsv(path: Path) -> list[dict[str, str]]:
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


def _finding(
    level: str,
    check: str,
    status: str,
    message: str,
    path: str | Path,
) -> ExportValidationFinding:
    return ExportValidationFinding(level, check, status, message, str(path))
