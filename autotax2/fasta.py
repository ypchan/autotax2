"""FASTA parsing and writing utilities."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

from .utils import md5_text, strip_gaps


@dataclass
class FastaRecord:
    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def md5(self) -> str:
        return md5_text(self.sequence.upper())


def parse_fasta(path: str | Path) -> Iterator[FastaRecord]:
    """Parse FASTA without loading the entire file."""
    seq_id = None
    desc = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield FastaRecord(seq_id, desc or seq_id, "".join(chunks))
                desc = line[1:].strip()
                seq_id = desc.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if seq_id is not None:
            yield FastaRecord(seq_id, desc or seq_id, "".join(chunks))


def write_fasta(records: Iterable[FastaRecord], path: str | Path, wrap: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for rec in records:
            out.write(f">{rec.description}\n")
            seq = rec.sequence
            for i in range(0, len(seq), wrap):
                out.write(seq[i : i + wrap] + "\n")


def strip_gaps_fasta(input_fa: str | Path, output_fa: str | Path) -> None:
    records = (
        FastaRecord(rec.id, rec.description, strip_gaps(rec.sequence)) for rec in parse_fasta(input_fa)
    )
    write_fasta(records, output_fa)


def fasta_lengths(path: str | Path) -> dict[str, int]:
    return {rec.id: rec.length for rec in parse_fasta(path)}
