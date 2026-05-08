from __future__ import annotations

import csv
import os
import shutil
import stat
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
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t1\t8\t42\t+\t.\tName=16S_rRNA;product=16S_rRNA\n"
        "D20_000002\tbarrnap\trRNA\t1\t8\t42\t+\t.\tName=16S_rRNA;product=16S_rRNA\n",
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
            "--threads",
            "8",
            "--barrnap-bin",
            str(barrnap_bin),
            "--min-rrna-len-archaea",
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


def test_prepare_dataset_extracts_gff_interval_minus_strand_and_flank(
    prepare_tmp_dir: Path,
) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(
        ">plus\n"
        "AAAACCCCGGGGTTTT\n"
        ">minus\n"
        "AAAACCCCGGGGTTTT\n",
        encoding="utf-8",
    )
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t5\t8\t42\t+\t.\tName=16S_rRNA;product=16S_rRNA\n"
        "D20_000002\tbarrnap\trRNA\t5\t8\t42\t-\t.\tName=16S_rRNA;product=16S_rRNA\n",
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
            "--barrnap-bin",
            str(barrnap_bin),
            "--flank",
            "1",
            "--min-rrna-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    records = {record.seq_id: record.sequence for record in read_fasta(dataset_dir / "barrnap.extracted.fa")}
    summary = {row["internal_seq_id"]: row for row in _read_tsv(dataset_dir / "barrnap.summary.tsv")}

    assert records["D20_000001"] == "ACCCCG"
    assert records["D20_000002"] == "CGGGGT"
    assert summary["D20_000001"]["extracted_start"] == "4"
    assert summary["D20_000001"]["extracted_end"] == "9"
    assert summary["D20_000002"]["reverse_complemented"] == "true"


def test_prepare_dataset_flank_clips_at_sequence_boundaries(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGTACGT\n", encoding="utf-8")
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t1\t4\t42\t+\t.\tproduct=16S ribosomal RNA\n",
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
            "--barrnap-bin",
            str(barrnap_bin),
            "--flank",
            "10",
            "--min-rrna-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    row = _read_tsv(build / "datasets" / "01_digester2020" / "barrnap.summary.tsv")[0]

    assert row["extracted_start"] == "1"
    assert row["extracted_end"] == "8"
    assert row["extracted_length"] == "8"


def test_prepare_dataset_longest_policy_chooses_longest_hit(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\n" + "A" * 30 + "\n", encoding="utf-8")
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t2\t6\t10\t+\t.\tproduct=SSU_rRNA\n"
        "D20_000001\tbarrnap\trRNA\t2\t20\t5\t+\t.\tproduct=SSU_rRNA\n",
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
            "--barrnap-bin",
            str(barrnap_bin),
            "--multi-rrna-policy",
            "longest",
            "--min-rrna-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    row = _read_tsv(build / "datasets" / "01_digester2020" / "barrnap.summary.tsv")[0]

    assert row["selected_hit_index"] == "2"
    assert row["extracted_length"] == "19"
    assert row["warning"] == "multiple_rrna_hits"


def test_prepare_dataset_fail_policy_rejects_multiple_hits(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\n" + "A" * 30 + "\n", encoding="utf-8")
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t2\t6\t10\t+\t.\tproduct=SSU_rRNA\n"
        "D20_000001\tbarrnap\trRNA\t2\t20\t5\t+\t.\tproduct=SSU_rRNA\n",
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
            "--barrnap-bin",
            str(barrnap_bin),
            "--multi-rrna-policy",
            "fail",
            "--min-rrna-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    row = _read_tsv(dataset_dir / "barrnap.summary.tsv")[0]

    assert row["status"] == "rejected_multiple_rrna_hits"
    assert read_fasta(dataset_dir / "barrnap.extracted.fa") == []


def test_prepare_dataset_short_rrna_is_rejected_by_domain_length(prepare_tmp_dir: Path) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGTACGT\n", encoding="utf-8")
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t1\t8\t42\t+\t.\tproduct=16S_rRNA\n",
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
            "Bacteria",
            "--barrnap-bin",
            str(barrnap_bin),
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    row = _read_tsv(dataset_dir / "barrnap.summary.tsv")[0]

    assert row["status"] == "rejected_short_rrna"
    assert row["warning"] == "short_rrna"
    assert read_fasta(dataset_dir / "barrnap.extracted.fa") == []


def test_prepare_dataset_outputs_real_newlines_tabs_and_valid_headers(
    prepare_tmp_dir: Path,
) -> None:
    build = prepare_tmp_dir / "build"
    fasta = prepare_tmp_dir / "input.fa"
    fasta.write_text(">seq1\nACGTACGT\n", encoding="utf-8")
    barrnap_bin = _fake_barrnap(
        prepare_tmp_dir,
        "D20_000001\tbarrnap\trRNA\t1\t8\t42\t+\t.\tproduct=16S_rRNA\n",
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
            "--barrnap-bin",
            str(barrnap_bin),
            "--min-rrna-len-archaea",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    fasta_content = (dataset_dir / "barrnap.extracted.fa").read_text(encoding="utf-8")
    tsv_content = (dataset_dir / "barrnap.summary.tsv").read_text(encoding="utf-8")

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


def _fake_barrnap(tmp_dir: Path, gff_text: str) -> Path:
    mock_gff = tmp_dir / f"mock_{uuid.uuid4().hex}.gff3"
    mock_gff.write_text(gff_text, encoding="utf-8")
    if os.name == "nt":
        script = tmp_dir / f"fake_barrnap_{uuid.uuid4().hex}.cmd"
        script.write_text(f"@echo off\r\ntype \"{mock_gff}\"\r\n", encoding="utf-8")
    else:
        script = tmp_dir / f"fake_barrnap_{uuid.uuid4().hex}.sh"
        script.write_text(f"#!/bin/sh\ncat '{mock_gff}'\n", encoding="utf-8")
        script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
