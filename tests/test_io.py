from __future__ import annotations

import gzip
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest

from autotax2.io import (
    FastaRecord,
    read_fasta,
    write_fasta,
    write_sequence_id_map,
    write_sequence_membership,
    write_unique_sequence_registry,
)
from autotax2.registry import SequenceRegistry, assign_internal_ids


@pytest.fixture
def io_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_io_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_read_plain_fasta_normalizes_and_ignores_blank_lines(io_tmp_dir: Path) -> None:
    fasta_path = io_tmp_dir / "input.fasta"
    fasta_path.write_text(
        ">seq1 some header\n"
        "acgu\n"
        "\n"
        "tt aa\n"
        ">seq2\n"
        "\n"
        "uucc\n",
        encoding="utf-8",
    )

    records = read_fasta(fasta_path)

    assert records == [
        FastaRecord(seq_id="seq1", header="seq1 some header", sequence="ACGTTTAA"),
        FastaRecord(seq_id="seq2", header="seq2", sequence="TTCC"),
    ]


def test_read_gzipped_fasta(io_tmp_dir: Path) -> None:
    fasta_path = io_tmp_dir / "input.fasta.gz"
    with gzip.open(fasta_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(">seq1\nacgu\n")

    records = read_fasta(fasta_path)

    assert records == [FastaRecord(seq_id="seq1", header="seq1", sequence="ACGT")]


def test_write_fasta_and_read_back(io_tmp_dir: Path) -> None:
    fasta_path = io_tmp_dir / "output.fasta"
    records = [
        FastaRecord(seq_id="seq1", header="seq1 description", sequence="acgu"),
        FastaRecord(seq_id="seq2", header="seq2", sequence="tttt"),
    ]

    write_fasta(records, fasta_path)

    assert read_fasta(fasta_path) == [
        FastaRecord(seq_id="seq1", header="seq1 description", sequence="ACGT"),
        FastaRecord(seq_id="seq2", header="seq2", sequence="TTTT"),
    ]


def test_write_fasta_wraps_sequence_at_80_characters(io_tmp_dir: Path) -> None:
    fasta_path = io_tmp_dir / "wrapped.fasta"
    write_fasta(
        [FastaRecord(seq_id="seq1", header="seq1", sequence="A" * 81)],
        fasta_path,
    )

    lines = fasta_path.read_text(encoding="utf-8").splitlines()

    assert lines == [">seq1", "A" * 80, "A"]


def test_write_gzipped_fasta(io_tmp_dir: Path) -> None:
    fasta_path = io_tmp_dir / "output.fasta.gz"

    write_fasta([FastaRecord(seq_id="seq1", header="seq1", sequence="acgu")], fasta_path)

    assert read_fasta(fasta_path) == [
        FastaRecord(seq_id="seq1", header="seq1", sequence="ACGT")
    ]


def test_tsv_writers_use_real_tabs_and_newlines(io_tmp_dir: Path) -> None:
    assigned = assign_internal_ids(
        [
            FastaRecord(seq_id="seq1", header="seq1 first", sequence="ACGT"),
            FastaRecord(seq_id="seq2", header="seq2 second", sequence="acgu"),
        ],
        dataset_name="dataset-20",
        prefix="D20",
    )
    registry = SequenceRegistry()
    registry.add_many(assigned)

    sequence_map_path = io_tmp_dir / "sequence_id_map.tsv"
    unique_path = io_tmp_dir / "unique_sequence_registry.tsv"
    membership_path = io_tmp_dir / "sequence_membership.tsv"

    write_sequence_id_map(assigned, sequence_map_path)
    write_unique_sequence_registry(registry.unique_records(), unique_path)
    write_sequence_membership(registry.membership_records(), membership_path)

    sequence_map = sequence_map_path.read_text(encoding="utf-8")
    unique_registry = unique_path.read_text(encoding="utf-8")
    membership = membership_path.read_text(encoding="utf-8")

    assert "\t" in sequence_map
    assert "\n" in sequence_map
    assert "\\t" not in sequence_map
    assert "\\n" not in sequence_map
    assert "unique_seq_id\tsequence_md5" in unique_registry
    assert "internal_seq_id\toriginal_seq_id" in membership
