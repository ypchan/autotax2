from __future__ import annotations

import csv
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
from typer.testing import CliRunner

from autotax2.cli import app
from autotax2.io import read_fasta


runner = CliRunner()


@pytest.fixture
def prepare_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_prepare_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_prepare_dataset_ids_mapping_rejection_and_duplicates(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(
        ">seqA first header\n"
        "ACGTACGT\n"
        ">seqA duplicate original id\n"
        "acguacgu\n"
        ">seqBad non atgc\n"
        "ACGTN\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "prepare-dataset",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(fasta),
            "--domain",
            "Archaea",
            "--min-ssu-len-archaea",
            "4",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    id_rows = _read_tsv(dataset_dir / "sequence_id_map.tsv")
    membership_rows = _read_tsv(dataset_dir / "sequence_membership.tsv")
    unique_rows = _read_tsv(dataset_dir / "unique_sequence_registry.tsv")
    summary_rows = _read_tsv(dataset_dir / "prepare_summary.tsv")

    assert [row["internal_seq_id"] for row in id_rows] == [
        "D20_000001",
        "D20_000002",
        "D20_000003",
    ]
    assert id_rows[0]["original_seq_id"] == "seqA"
    assert id_rows[0]["original_header"] == "seqA first header"
    assert id_rows[2]["rejected"] == "true"
    assert id_rows[2]["reject_reason"] == "non_atgc"
    assert len(unique_rows) == 1
    assert membership_rows[0]["is_duplicate_sequence"] == "false"
    assert membership_rows[1]["is_duplicate_sequence"] == "true"
    assert summary_rows[0]["duplicate_md5_sequences"] == "1"
    assert summary_rows[0]["rejected_non_atgc"] == "1"
    assert summary_rows[0]["input_contract"] == "externally_processed_ssu"


def test_prepare_dataset_writes_prepared_ssu_without_internal_extraction_outputs(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(
        ">plus\n"
        "AAAACCCCGGGGTTTT\n"
        ">second\n"
        "TTTTCCCCAAAAGGGG\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "prepare-dataset",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(fasta),
            "--domain",
            "Archaea",
            "--min-ssu-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    records = {record.seq_id: record.sequence for record in read_fasta(dataset_dir / "prepared.ssu.fa")}

    assert records == {
        "D20_000001": "AAAACCCCGGGGTTTT",
        "D20_000002": "TTTTCCCCAAAAGGGG",
    }
    assert not list(dataset_dir.glob("*.gff3"))
    assert not list(dataset_dir.glob("*.extracted.fa"))


def test_prepare_dataset_short_ssu_is_rejected_by_domain_length(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGTACGT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "prepare-dataset",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(fasta),
            "--domain",
            "Bacteria",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    row = _read_tsv(dataset_dir / "sequence_id_map.tsv")[0]
    summary = _read_tsv(dataset_dir / "prepare_summary.tsv")[0]

    assert row["rejected"] == "true"
    assert row["reject_reason"] == "short_ssu"
    assert summary["rejected_short_ssu"] == "1"
    assert read_fasta(dataset_dir / "prepared.ssu.fa") == []


def test_prepare_dataset_outputs_real_newlines_tabs_and_valid_headers(
    prepare_tmp_dir: Path,
) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGTACGT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "prepare-dataset",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(fasta),
            "--domain",
            "Archaea",
            "--min-ssu-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    fasta_content = (dataset_dir / "prepared.ssu.fa").read_text(encoding="utf-8")
    tsv_content = (dataset_dir / "prepare_summary.tsv").read_text(encoding="utf-8")

    assert fasta_content.startswith(">D20_000001\n")
    assert "\n" in fasta_content
    assert "\t" in tsv_content
    assert "\n" in tsv_content
    assert "\\t" not in tsv_content
    assert "\\n" not in tsv_content


def test_prepare_dataset_prefix_conflict_raises_error(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    registry = build / "registry"
    registry.mkdir(parents=True)
    (registry / "dataset_registry.tsv").write_text(
        "dataset_name\tprefix\tadd_order\tinput_fasta\tinput_md5\tdomain\tstatus\tdataset_dir\n"
        "other\tD20\t1\tinput.fa\tabc\tArchaea\tprepared\tdatasets/01_other\n",
        encoding="utf-8",
    )
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "prepare-dataset",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(fasta),
            "--domain",
            "Archaea",
        ],
    )

    assert result.exit_code != 0
    assert "already registered" in str(result.exception)


def test_prepare_dataset_name_with_different_prefix_raises_error(
    prepare_tmp_dir: Path,
) -> None:
    build = prepare_tmp_dir / "build"
    registry = build / "registry"
    registry.mkdir(parents=True)
    (registry / "dataset_registry.tsv").write_text(
        "dataset_name\tprefix\tadd_order\tinput_fasta\tinput_md5\tdomain\tstatus\tdataset_dir\n"
        "digester2020\tD19\t1\tinput.fa\tabc\tArchaea\tprepared\tdatasets/01_digester2020\n",
        encoding="utf-8",
    )
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "prepare-dataset",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(fasta),
            "--domain",
            "Archaea",
        ],
    )

    assert result.exit_code != 0
    assert "already registered" in str(result.exception)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
