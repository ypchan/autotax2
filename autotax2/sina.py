"""SINA orientation correction with loose settings."""

from __future__ import annotations

import csv
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from autotax2.io import FastaRecord, read_fasta, write_fasta


DEFAULT_ORIENTATION_MODE = "loose"
EXPECTED_SINA_VERSION = "not_pinned"
SINA_SUMMARY_FIELDS = [
    "internal_seq_id",
    "dataset",
    "input_length",
    "output_length",
    "sina_status",
    "strand",
    "orientation_confidence",
    "sequence_changed",
    "fallback_used",
    "warning",
]
TOOL_VERSION_FIELDS = [
    "tool",
    "expected_version",
    "detected_version",
    "status",
    "command",
]
SINA_CANDIDATE_FIELDS = [
    "query",
    "target",
    "sina_score",
    "source_field",
    "rank",
]


@dataclass(frozen=True)
class SinaRunSummary:
    """Summary from SINA orientation correction."""

    dataset: str
    dataset_dir: Path
    oriented_records: int
    fallback_used: bool


def default_orientation_mode() -> str:
    """Return the default SINA orientation mode."""
    return DEFAULT_ORIENTATION_MODE


def build_sina_command(
    input_fasta: str | Path,
    output_fasta: str | Path,
    threads: int = 4,
    sina_bin: str = "sina",
    reference: str | Path | None = None,
    search_candidates: bool = False,
    search_output_csv: str | Path | None = None,
    search_db: str | Path | None = None,
    search_min_sim: float = 0.5,
    search_max_result: int = 10,
    search_kmer_candidates: int = 1000,
) -> list[str]:
    """Build a loose SINA command for orientation/alignment support."""
    command = [
        sina_bin,
        "-i",
        str(input_fasta),
        "-o",
        str(output_fasta),
    ]
    if search_candidates and search_output_csv is not None:
        command.append(str(search_output_csv))
    command.extend(["--threads", str(threads)])
    if reference is not None:
        command.extend(["--ptdb", str(reference)])
    if search_candidates:
        command.extend(
            [
                "--search",
                "--search-min-sim",
                f"{search_min_sim:.3f}",
                "--search-max-result",
                str(search_max_result),
                "--search-kmer-candidates",
                str(search_kmer_candidates),
                "--fields",
                "name,nearest_slv",
            ]
        )
        if search_db is not None:
            command.extend(["--search-db", str(search_db)])
    return command


def check_sina_version(sina_bin: str = "sina") -> str:
    """Return the detected SINA version string, or an empty string if unavailable."""
    try:
        completed = subprocess.run(
            [sina_bin, "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return ""

    output = f"{completed.stdout}\n{completed.stderr}".strip()
    match = re.search(r"(\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?)", output)
    return match.group(1) if match else ""


def parse_sina_candidate_csv(path: str | Path) -> list[dict[str, str]]:
    """Parse SINA CSV search candidates into a compact TSV-ready row list."""
    csv_path = Path(path)
    if not csv_path.exists():
        return []
    rows: list[dict[str, str]] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = ","
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
        except csv.Error:
            pass
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            query = _first_value(row, "name", "query", "seq_id", "internal_seq_id")
            nearest = _first_value(row, "nearest_slv", "nearest", "search", "search_results")
            if not query or not nearest:
                continue
            for rank, candidate in enumerate(_split_sina_nearest_entries(nearest), start=1):
                target, score = _parse_sina_nearest_entry(candidate)
                if not target:
                    continue
                rows.append(
                    {
                        "query": query,
                        "target": target,
                        "sina_score": score,
                        "source_field": "nearest_slv",
                        "rank": str(rank),
                    }
                )
    return rows


def _split_sina_nearest_entries(value: str) -> list[str]:
    clean = value.strip().strip("[]")
    if not clean:
        return []
    if "|" in clean:
        return [part.strip() for part in clean.split("|") if part.strip()]
    if ";" in clean:
        return [part.strip() for part in clean.split(";") if part.strip()]
    if "," in clean:
        return [part.strip() for part in clean.split(",") if part.strip()]
    return [part.strip() for part in clean.split() if part.strip()]


def _parse_sina_nearest_entry(value: str) -> tuple[str, str]:
    clean = value.strip().strip("'\"")
    if not clean:
        return "", ""
    score = ""
    if "~" in clean:
        clean, score = clean.rsplit("~", maxsplit=1)
    elif ":" in clean and re.search(r":\d+(?:\.\d+)?$", clean):
        clean, score = clean.rsplit(":", maxsplit=1)
    target = clean.strip()
    return target, score.strip()


def _first_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = row.get(key, "")
        if value:
            return value
    return ""


def orient_dataset_with_sina(
    build: str | Path,
    dataset: str,
    threads: int = 4,
    sina_bin: str = "sina",
    reference: str | Path | None = None,
    strict_tool_version: bool = False,
    allow_sina_failure: bool = True,
    fallback_copy_original: bool = True,
    min_sina_identity: float = 0.0,
    min_sina_score: float = 0.0,
    search_candidates: bool = False,
    search_db: str | Path | None = None,
    search_min_sim: float = 0.5,
    search_max_result: int = 10,
    search_kmer_candidates: int = 1000,
) -> SinaRunSummary:
    """Orient a prepared dataset with SINA, falling back conservatively by default."""
    if min_sina_identity > 0.0 or min_sina_score > 0.0:
        raise NotImplementedError(
            "SINA identity/score threshold filtering is not implemented because "
            "the current parser does not read those metrics from SINA output. "
            "Leave --min-sina-identity and --min-sina-score at 0.0."
        )

    build_dir = Path(build)
    dataset_dir = _find_dataset_dir(build_dir, dataset)
    input_fasta = dataset_dir / "prepared.ssu.fa"
    output_fasta = dataset_dir / "sina.oriented.fa"
    candidates_csv = dataset_dir / "sina.candidates.csv"
    candidates_tsv = dataset_dir / "sina.candidates.tsv"
    summary_path = dataset_dir / "sina.summary.tsv"
    log_path = dataset_dir / "sina.log"
    tool_versions_path = dataset_dir / "tool_versions.tsv"
    input_records = read_fasta(input_fasta)
    input_by_id = {record.seq_id: record for record in input_records}
    command = build_sina_command(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        threads=threads,
        sina_bin=sina_bin,
        reference=reference,
        search_candidates=search_candidates,
        search_output_csv=candidates_csv if search_candidates else None,
        search_db=search_db,
        search_min_sim=search_min_sim,
        search_max_result=search_max_result,
        search_kmer_candidates=search_kmer_candidates,
    )
    detected_version = check_sina_version(sina_bin)
    version_status = "ok" if detected_version else "warning_unparsed"
    if strict_tool_version and not detected_version:
        raise RuntimeError(f"Could not parse SINA version from {sina_bin} --version.")

    command_failed = False
    stderr_text = ""
    try:
        completed = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stderr_text = completed.stderr or ""
    except (FileNotFoundError, OSError, subprocess.CalledProcessError) as exc:
        command_failed = True
        stderr_text = getattr(exc, "stderr", "") or str(exc)
        if not allow_sina_failure:
            log_path.write_text(stderr_text, encoding="utf-8")
            _write_tool_version(
                tool_versions_path=tool_versions_path,
                detected_version=detected_version,
                status="failed",
                command=command,
            )
            raise RuntimeError(f"SINA command failed and failures are not allowed: {sina_bin}") from exc

    log_path.write_text(stderr_text, encoding="utf-8")

    if command_failed:
        if not fallback_copy_original:
            _write_tsv(
                [
                    _summary_row(
                        record,
                        dataset,
                        status="failed",
                        strand="unknown",
                        confidence="unknown",
                        output_length=0,
                        sequence_changed=False,
                        fallback_used=False,
                        warning="sina_failed",
                    )
                    for record in input_records
                ],
                summary_path,
                SINA_SUMMARY_FIELDS,
            )
            _write_tool_version(tool_versions_path, detected_version, "failed", command)
            if search_candidates:
                _write_tsv([], candidates_tsv, SINA_CANDIDATE_FIELDS)
            return SinaRunSummary(dataset=dataset, dataset_dir=dataset_dir, oriented_records=0, fallback_used=False)
        write_fasta(input_records, output_fasta)
        rows = [
            _summary_row(
                record,
                dataset,
                status="sina_failed_fallback_original",
                strand="unknown",
                confidence="unknown",
                output_length=len(record.sequence),
                sequence_changed=False,
                fallback_used=True,
                warning="sina_failed",
            )
            for record in input_records
        ]
        _write_tsv(rows, summary_path, SINA_SUMMARY_FIELDS)
        _write_tool_version(tool_versions_path, detected_version, "failed_fallback", command)
        if search_candidates:
            _write_tsv([], candidates_tsv, SINA_CANDIDATE_FIELDS)
        return SinaRunSummary(
            dataset=dataset,
            dataset_dir=dataset_dir,
            oriented_records=len(input_records),
            fallback_used=True,
        )

    sina_output_records = read_fasta(output_fasta) if output_fasta.exists() else []
    output_by_id = {record.seq_id: record for record in sina_output_records}
    final_records: list[FastaRecord] = []
    summary_rows: list[dict[str, str]] = []
    used_fallback = False

    for input_record in input_records:
        output_record = output_by_id.get(input_record.seq_id)
        if output_record is None:
            if fallback_copy_original:
                used_fallback = True
                final_records.append(input_record)
                summary_rows.append(
                    _summary_row(
                        input_record,
                        dataset,
                        status="sina_missing_output",
                        strand="unknown",
                        confidence="unknown",
                        output_length=len(input_record.sequence),
                        sequence_changed=False,
                        fallback_used=True,
                        warning="missing_sina_output",
                    )
                )
                continue
            summary_rows.append(
                _summary_row(
                    input_record,
                    dataset,
                    status="failed",
                    strand="unknown",
                    confidence="unknown",
                    output_length=0,
                    sequence_changed=False,
                    fallback_used=False,
                    warning="missing_sina_output",
                )
            )
            continue

        oriented_record, summary = _classify_sina_output(input_record, output_record, dataset)
        final_records.append(oriented_record)
        summary_rows.append(summary)

    write_fasta(final_records, output_fasta)
    _write_tsv(summary_rows, summary_path, SINA_SUMMARY_FIELDS)
    if search_candidates:
        _write_tsv(parse_sina_candidate_csv(candidates_csv), candidates_tsv, SINA_CANDIDATE_FIELDS)
    _write_tool_version(tool_versions_path, detected_version, version_status, command)

    return SinaRunSummary(
        dataset=dataset,
        dataset_dir=dataset_dir,
        oriented_records=len(final_records),
        fallback_used=used_fallback,
    )


def reverse_complement(sequence: str) -> str:
    """Return the DNA reverse complement."""
    return sequence.translate(str.maketrans("ATGC", "TACG"))[::-1]


def _classify_sina_output(
    input_record: FastaRecord,
    output_record: FastaRecord,
    dataset: str,
) -> tuple[FastaRecord, dict[str, str]]:
    oriented_record = FastaRecord(
        seq_id=input_record.seq_id,
        header=input_record.seq_id,
        sequence=output_record.sequence,
    )
    if output_record.sequence == input_record.sequence:
        return oriented_record, _summary_row(
            input_record,
            dataset,
            status="oriented",
            strand="plus",
            confidence="high",
            output_length=len(output_record.sequence),
            sequence_changed=False,
            fallback_used=False,
            warning="",
        )
    if output_record.sequence == reverse_complement(input_record.sequence):
        return oriented_record, _summary_row(
            input_record,
            dataset,
            status="oriented",
            strand="minus",
            confidence="high",
            output_length=len(output_record.sequence),
            sequence_changed=True,
            fallback_used=False,
            warning="",
        )
    return oriented_record, _summary_row(
        input_record,
        dataset,
        status="oriented_modified",
        strand="unknown_or_aligned_modified",
        confidence="low",
        output_length=len(output_record.sequence),
        sequence_changed=True,
        fallback_used=False,
        warning="sina_modified_sequence",
    )


def _summary_row(
    input_record: FastaRecord,
    dataset: str,
    status: str,
    strand: str,
    confidence: str,
    output_length: int,
    sequence_changed: bool,
    fallback_used: bool,
    warning: str,
) -> dict[str, str]:
    return {
        "internal_seq_id": input_record.seq_id,
        "dataset": dataset,
        "input_length": str(len(input_record.sequence)),
        "output_length": str(output_length),
        "sina_status": status,
        "strand": strand,
        "orientation_confidence": confidence,
        "sequence_changed": _bool_text(sequence_changed),
        "fallback_used": _bool_text(fallback_used),
        "warning": warning,
    }


def _find_dataset_dir(build_dir: Path, dataset: str) -> Path:
    datasets_dir = build_dir / "datasets"
    if not datasets_dir.exists():
        raise FileNotFoundError(f"Dataset directory root not found: {datasets_dir}")

    matches = sorted(
        path
        for path in datasets_dir.iterdir()
        if path.is_dir() and path.name.split("_", maxsplit=1)[-1] == dataset
    )
    if not matches:
        raise FileNotFoundError(f"Prepared dataset not found: {dataset}")
    if len(matches) > 1:
        raise ValueError(f"Multiple prepared dataset directories match {dataset!r}.")
    return matches[0]


def _write_tool_version(
    tool_versions_path: Path,
    detected_version: str,
    status: str,
    command: list[str],
) -> None:
    rows = [
        row
        for row in _read_tsv(tool_versions_path)
        if row.get("tool") != "sina"
    ]
    rows.append(
        {
            "tool": "sina",
            "expected_version": EXPECTED_SINA_VERSION,
            "detected_version": detected_version,
            "status": status,
            "command": " ".join(command),
        }
    )
    _write_tsv(rows, tool_versions_path, TOOL_VERSION_FIELDS)


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


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
