from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .external import run_command
from .logging import step, success, warning
from .threads import validate_threads
from .utils import ensure_dir, ensure_file, read_fasta


ORIENTATION_REPORT_FIELDS = [
    "sequence_id",
    "length",
    "method",
    "decision",
    "reverse_complemented",
    "forward_score",
    "reverse_complement_score",
    "sina_sequence_found",
    "sina_aligned_length",
    "sina_ungapped_length",
    "reason",
]


DNA_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "u": "a",
        "r": "y",
        "y": "r",
        "s": "s",
        "w": "w",
        "k": "m",
        "m": "k",
        "b": "v",
        "d": "h",
        "h": "d",
        "v": "b",
        "n": "n",
    }
)


@dataclass(frozen=True)
class OrientationDecision:
    sequence_id: str
    header: str
    original_sequence: str
    oriented_sequence: str
    method: str
    decision: str
    reverse_complemented: bool
    forward_score: float
    reverse_complement_score: float
    sina_sequence_found: bool
    sina_aligned_length: int
    sina_ungapped_length: int
    reason: str


def sequence_id_from_header(header: str) -> str:
    return header.split()[0].split(";")[0]


def normalize_sequence(sequence: str) -> str:
    return (
        "".join(str(sequence).split())
        .replace("-", "")
        .replace(".", "")
        .upper()
        .replace("U", "T")
    )


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1].upper().replace("U", "T")


def fasta_to_dict(path: str | Path) -> Dict[str, Tuple[str, str]]:
    records: Dict[str, Tuple[str, str]] = {}

    for header, sequence in read_fasta(path):
        seq_id = sequence_id_from_header(header)
        records[seq_id] = (header, normalize_sequence(sequence))

    return records


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")

    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def write_report(path: str | Path, decisions: Iterable[OrientationDecision]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ORIENTATION_REPORT_FIELDS,
            delimiter="\t",
        )
        writer.writeheader()

        for item in decisions:
            writer.writerow(
                {
                    "sequence_id": item.sequence_id,
                    "length": len(item.original_sequence),
                    "method": item.method,
                    "decision": item.decision,
                    "reverse_complemented": "yes" if item.reverse_complemented else "no",
                    "forward_score": f"{item.forward_score:.6f}",
                    "reverse_complement_score": f"{item.reverse_complement_score:.6f}",
                    "sina_sequence_found": "yes" if item.sina_sequence_found else "no",
                    "sina_aligned_length": item.sina_aligned_length,
                    "sina_ungapped_length": item.sina_ungapped_length,
                    "reason": item.reason,
                }
            )


def kmer_set(sequence: str, k: int) -> set[str]:
    sequence = normalize_sequence(sequence)

    if len(sequence) < k:
        return {sequence} if sequence else set()

    return {
        sequence[i:i + k]
        for i in range(0, len(sequence) - k + 1)
        if "N" not in sequence[i:i + k]
    }


def kmer_jaccard(a: str, b: str, *, k: int = 15) -> float:
    kmers_a = kmer_set(a, k)
    kmers_b = kmer_set(b, k)

    if not kmers_a or not kmers_b:
        return 0.0

    intersection = len(kmers_a & kmers_b)
    union = len(kmers_a | kmers_b)

    if union == 0:
        return 0.0

    return intersection / union


def run_sina(
    input_fasta: str | Path,
    output_fasta: str | Path,
    *,
    sina: str | Path,
    sina_ref: str | Path,
    threads: int,
    turn: str = "all",
    extra_args: Optional[Sequence[str]] = None,
) -> str:
    """Run SINA to orient 16S sequences.

    The SINA output is used only to infer orientation. Downstream AutoTax2
    steps use the original ungapped sequence, reverse-complemented when needed.
    """

    threads = validate_threads(threads)

    input_fasta = ensure_file(input_fasta, "input FASTA")
    sina_ref = ensure_file(sina_ref, "SINA reference database")

    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    command: List[str | Path] = [
        sina,
        "-i",
        input_fasta,
        "-o",
        output_fasta,
        "--db",
        sina_ref,
        "--turn",
        turn,
        "--threads",
        str(threads),
    ]

    if extra_args:
        command.extend(extra_args)

    step("Running SINA orientation check")
    run_command(command)

    return str(output_fasta)


def infer_decisions_from_sina_output(
    input_records: Dict[str, Tuple[str, str]],
    sina_records: Dict[str, Tuple[str, str]],
    *,
    score_margin: float = 0.05,
    min_score: float = 0.20,
    kmer_size: int = 15,
) -> List[OrientationDecision]:
    decisions: List[OrientationDecision] = []

    for seq_id, (header, original_sequence) in input_records.items():
        sina_record = sina_records.get(seq_id)

        if sina_record is None:
            decisions.append(
                OrientationDecision(
                    sequence_id=seq_id,
                    header=header,
                    original_sequence=original_sequence,
                    oriented_sequence=original_sequence,
                    method="sina",
                    decision="ambiguous",
                    reverse_complemented=False,
                    forward_score=0.0,
                    reverse_complement_score=0.0,
                    sina_sequence_found=False,
                    sina_aligned_length=0,
                    sina_ungapped_length=0,
                    reason="sequence_not_found_in_sina_output",
                )
            )
            continue

        _sina_header, sina_sequence = sina_record
        sina_ungapped = normalize_sequence(sina_sequence)

        rc_sequence = reverse_complement(original_sequence)

        forward_score = kmer_jaccard(
            sina_ungapped,
            original_sequence,
            k=kmer_size,
        )
        reverse_score = kmer_jaccard(
            sina_ungapped,
            rc_sequence,
            k=kmer_size,
        )

        if reverse_score >= min_score and reverse_score > forward_score + score_margin:
            decision = "reverse_complement"
            reverse_complemented = True
            oriented_sequence = rc_sequence
            reason = "sina_output_matches_reverse_complement"

        elif forward_score >= min_score and forward_score >= reverse_score + score_margin:
            decision = "keep"
            reverse_complemented = False
            oriented_sequence = original_sequence
            reason = "sina_output_matches_original_orientation"

        else:
            decision = "ambiguous"
            reverse_complemented = False
            oriented_sequence = original_sequence
            reason = "sina_orientation_evidence_ambiguous"

        decisions.append(
            OrientationDecision(
                sequence_id=seq_id,
                header=header,
                original_sequence=original_sequence,
                oriented_sequence=oriented_sequence,
                method="sina",
                decision=decision,
                reverse_complemented=reverse_complemented,
                forward_score=forward_score,
                reverse_complement_score=reverse_score,
                sina_sequence_found=True,
                sina_aligned_length=len(sina_sequence),
                sina_ungapped_length=len(sina_ungapped),
                reason=reason,
            )
        )

    return decisions


def write_orientation_outputs(
    decisions: List[OrientationDecision],
    outdir: str | Path,
    *,
    output_fasta_name: str = "oriented_input.fasta",
) -> Dict[str, str]:
    outdir = ensure_dir(outdir)

    oriented_fasta = outdir / output_fasta_name
    report_tsv = outdir / "orientation_report.tsv"
    reverse_ids = outdir / "reverse_complemented_ids.txt"
    ambiguous_ids = outdir / "ambiguous_orientation_ids.txt"

    with oriented_fasta.open("w", encoding="utf-8") as handle:
        for item in decisions:
            write_fasta_record(handle, item.header, item.oriented_sequence)

    write_report(report_tsv, decisions)

    with reverse_ids.open("w", encoding="utf-8") as handle:
        for item in decisions:
            if item.reverse_complemented:
                handle.write(item.sequence_id + "\n")

    with ambiguous_ids.open("w", encoding="utf-8") as handle:
        for item in decisions:
            if item.decision == "ambiguous":
                handle.write(item.sequence_id + "\n")

    return {
        "oriented_fasta": str(oriented_fasta),
        "orientation_report": str(report_tsv),
        "reverse_complemented_ids": str(reverse_ids),
        "ambiguous_orientation_ids": str(ambiguous_ids),
    }


def orient_with_sina(
    input_fasta: str | Path,
    outdir: str | Path,
    *,
    sina: str | Path,
    sina_ref: str | Path,
    threads: int = 4,
    score_margin: float = 0.05,
    min_score: float = 0.20,
    kmer_size: int = 15,
    sina_extra_args: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    """Orient near-full-length or full-length 16S sequences with SINA.

    SINA is used to infer orientation. The final oriented FASTA contains
    ungapped original sequences, reverse-complemented only when SINA evidence
    supports that decision.
    """

    threads = validate_threads(threads)

    input_fasta = ensure_file(input_fasta, "input FASTA")
    outdir = ensure_dir(outdir)

    sina_aligned = outdir / "sina_oriented_aligned.fasta"

    input_records = fasta_to_dict(input_fasta)

    if not input_records:
        raise ValueError(f"No FASTA records found in {input_fasta}")

    run_sina(
        input_fasta=input_fasta,
        output_fasta=sina_aligned,
        sina=sina,
        sina_ref=sina_ref,
        threads=threads,
        turn="all",
        extra_args=sina_extra_args,
    )

    sina_records = fasta_to_dict(sina_aligned)

    decisions = infer_decisions_from_sina_output(
        input_records,
        sina_records,
        score_margin=score_margin,
        min_score=min_score,
        kmer_size=kmer_size,
    )

    outputs = write_orientation_outputs(decisions, outdir)

    outputs["sina_aligned_fasta"] = str(sina_aligned)

    n_reverse = sum(1 for item in decisions if item.reverse_complemented)
    n_ambiguous = sum(1 for item in decisions if item.decision == "ambiguous")

    success(
        f"Orientation finished: {len(decisions)} sequences, "
        f"{n_reverse} reverse-complemented, {n_ambiguous} ambiguous"
    )

    return outputs


def copy_without_orientation(
    input_fasta: str | Path,
    outdir: str | Path,
    *,
    output_fasta_name: str = "oriented_input.fasta",
) -> Dict[str, str]:
    """Copy input FASTA and write a report marking orientation as skipped."""

    input_fasta = ensure_file(input_fasta, "input FASTA")
    outdir = ensure_dir(outdir)

    output_fasta = outdir / output_fasta_name
    shutil.copyfile(input_fasta, output_fasta)

    decisions: List[OrientationDecision] = []

    for header, sequence in read_fasta(input_fasta):
        seq_id = sequence_id_from_header(header)
        seq = normalize_sequence(sequence)

        decisions.append(
            OrientationDecision(
                sequence_id=seq_id,
                header=header,
                original_sequence=seq,
                oriented_sequence=seq,
                method="none",
                decision="skipped",
                reverse_complemented=False,
                forward_score=0.0,
                reverse_complement_score=0.0,
                sina_sequence_found=False,
                sina_aligned_length=0,
                sina_ungapped_length=0,
                reason="orientation_skipped",
            )
        )

    write_report(outdir / "orientation_report.tsv", decisions)

    reverse_ids = outdir / "reverse_complemented_ids.txt"
    ambiguous_ids = outdir / "ambiguous_orientation_ids.txt"

    reverse_ids.write_text("", encoding="utf-8")
    ambiguous_ids.write_text("", encoding="utf-8")

    return {
        "oriented_fasta": str(output_fasta),
        "orientation_report": str(outdir / "orientation_report.tsv"),
        "reverse_complemented_ids": str(reverse_ids),
        "ambiguous_orientation_ids": str(ambiguous_ids),
    }