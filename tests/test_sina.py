from __future__ import annotations

import csv
import shutil
import subprocess
import uuid
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

from autotax2.io import read_fasta
from autotax2.sina import (
    build_sina_command,
    check_sina_version,
    orient_dataset_with_sina,
    parse_sina_candidate_csv,
    reverse_complement,
)


@pytest.fixture
def sina_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_sina_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_sina_command_builder_includes_input_output_threads() -> None:
    command = build_sina_command(
        input_fasta="in.fa",
        output_fasta="out.fa",
        threads=8,
        sina_bin="sina",
    )

    assert command == ["sina", "-i", "in.fa", "-o", "out.fa", "--threads", "8"]


def test_sina_command_builder_includes_reference() -> None:
    command = build_sina_command("in.fa", "out.fa", reference="ref.ptdb")

    assert "--ptdb" in command
    assert "ref.ptdb" in command


def test_sina_command_builder_can_request_search_candidates() -> None:
    command = build_sina_command(
        "in.fa",
        "out.fa",
        search_candidates=True,
        search_output_csv="candidates.csv",
        search_db="search.arb",
        search_min_sim=0.650,
        search_max_result=100,
    )

    assert command[command.index("-o") + 1 : command.index("--threads")] == ["out.fa", "candidates.csv"]
    assert "--search" in command
    assert "--search-db" in command
    assert "search.arb" in command
    assert "--fields" in command
    assert command[command.index("--fields") + 1] == "name,nearest_slv"


def test_sina_candidate_search_defaults_are_loose() -> None:
    command = build_sina_command(
        "in.fa",
        "out.fa",
        search_candidates=True,
        search_output_csv="candidates.csv",
    )

    assert command[command.index("--search-min-sim") + 1] == "0.500"
    assert command[command.index("--search-max-result") + 1] == "10"


def test_parse_sina_candidate_csv_reads_nearest_slv(sina_tmp_dir: Path) -> None:
    path = sina_tmp_dir / "sina.candidates.csv"
    path.write_text(
        "name,nearest_slv\n"
        'D20_000001,"REF_A~0.98;REF_B~0.97"\n',
        encoding="utf-8",
    )

    rows = parse_sina_candidate_csv(path)

    assert rows == [
        {"query": "D20_000001", "target": "REF_A", "sina_score": "0.98", "source_field": "nearest_slv", "rank": "1"},
        {"query": "D20_000001", "target": "REF_B", "sina_score": "0.97", "source_field": "nearest_slv", "rank": "2"},
    ]


def test_version_parser_handles_typical_version_string(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args[0], 0, stdout="SINA 1.7.2\n", stderr="")

    monkeypatch.setattr("autotax2.sina.subprocess.run", fake_run)

    assert check_sina_version("sina") == "1.7.2"


def test_orient_sina_output_equal_input_is_plus_high_confidence(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTACGT"})
    _mock_sina_run(monkeypatch, lambda: ">D20_000001\nACGTACGT\n")

    orient_dataset_with_sina(sina_tmp_dir / "build", "digester2020", threads=8)

    rows = _read_tsv(dataset_dir / "sina.summary.tsv")
    records = read_fasta(dataset_dir / "sina.oriented.fa")

    assert rows[0]["sina_status"] == "oriented"
    assert rows[0]["strand"] == "plus"
    assert rows[0]["orientation_confidence"] == "high"
    assert rows[0]["sequence_changed"] == "false"
    assert [record.seq_id for record in records] == ["D20_000001"]


def test_orient_sina_search_candidates_writes_tsv(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTACGT"})
    _mock_sina_run_with_candidates(
        monkeypatch,
        fasta_text=">D20_000001\nACGTACGT\n",
        candidate_csv='name,nearest_slv\nD20_000001,"REF_A~0.98;REF_B~0.97"\n',
    )

    orient_dataset_with_sina(
        sina_tmp_dir / "build",
        "digester2020",
        search_candidates=True,
    )

    rows = _read_tsv(dataset_dir / "sina.candidates.tsv")

    assert [row["target"] for row in rows] == ["REF_A", "REF_B"]
    assert rows[0]["query"] == "D20_000001"


def test_orient_sina_reverse_complement_is_minus_high_confidence(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTCC"})
    reversed_sequence = reverse_complement("ACGTCC")
    _mock_sina_run(monkeypatch, lambda: f">D20_000001\n{reversed_sequence}\n")

    orient_dataset_with_sina(sina_tmp_dir / "build", "digester2020")

    rows = _read_tsv(dataset_dir / "sina.summary.tsv")
    records = read_fasta(dataset_dir / "sina.oriented.fa")

    assert rows[0]["sina_status"] == "oriented"
    assert rows[0]["strand"] == "minus"
    assert rows[0]["orientation_confidence"] == "high"
    assert rows[0]["sequence_changed"] == "true"
    assert records[0].sequence == reversed_sequence


def test_orient_sina_modified_output_is_low_confidence(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTACGT"})
    _mock_sina_run(monkeypatch, lambda: ">D20_000001\nACGTACGA\n")

    orient_dataset_with_sina(sina_tmp_dir / "build", "digester2020")

    row = _read_tsv(dataset_dir / "sina.summary.tsv")[0]

    assert row["sina_status"] == "oriented_modified"
    assert row["strand"] == "unknown_or_aligned_modified"
    assert row["orientation_confidence"] == "low"
    assert row["warning"] == "sina_modified_sequence"


def test_orient_sina_missing_output_falls_back_to_original(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(
        sina_tmp_dir,
        {"D20_000001": "ACGTACGT", "D20_000002": "TTTTCCCC"},
    )
    _mock_sina_run(monkeypatch, lambda: ">D20_000001\nACGTACGT\n")

    orient_dataset_with_sina(sina_tmp_dir / "build", "digester2020")

    rows = {row["internal_seq_id"]: row for row in _read_tsv(dataset_dir / "sina.summary.tsv")}
    records = {record.seq_id: record.sequence for record in read_fasta(dataset_dir / "sina.oriented.fa")}

    assert rows["D20_000002"]["sina_status"] == "sina_missing_output"
    assert rows["D20_000002"]["fallback_used"] == "true"
    assert records["D20_000002"] == "TTTTCCCC"


def test_orient_sina_failure_fallback_copies_all_records(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(
        sina_tmp_dir,
        {"D20_000001": "ACGTACGT", "D20_000002": "TTTTCCCC"},
    )
    _mock_sina_failure(monkeypatch)

    summary = orient_dataset_with_sina(sina_tmp_dir / "build", "digester2020")

    rows = _read_tsv(dataset_dir / "sina.summary.tsv")
    records = {record.seq_id: record.sequence for record in read_fasta(dataset_dir / "sina.oriented.fa")}

    assert summary.fallback_used is True
    assert {row["sina_status"] for row in rows} == {"sina_failed_fallback_original"}
    assert records == {"D20_000001": "ACGTACGT", "D20_000002": "TTTTCCCC"}


def test_orient_sina_failure_can_be_fatal(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTACGT"})
    _mock_sina_failure(monkeypatch)

    with pytest.raises(RuntimeError):
        orient_dataset_with_sina(
            sina_tmp_dir / "build",
            "digester2020",
            allow_sina_failure=False,
        )


def test_orient_sina_rejects_unimplemented_metric_thresholds(sina_tmp_dir: Path) -> None:
    _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTACGT"})

    with pytest.raises(NotImplementedError, match="threshold filtering is not implemented"):
        orient_dataset_with_sina(
            sina_tmp_dir / "build",
            "digester2020",
            min_sina_identity=0.9,
        )


def test_sina_summary_tsv_has_real_tabs_and_newlines(
    sina_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset_dir = _make_dataset(sina_tmp_dir, {"D20_000001": "ACGTACGT"})
    _mock_sina_run(monkeypatch, lambda: ">D20_000001\nACGTACGT\n")

    orient_dataset_with_sina(sina_tmp_dir / "build", "digester2020")

    content = (dataset_dir / "sina.summary.tsv").read_text(encoding="utf-8")

    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content


def _make_dataset(tmp_dir: Path, records: dict[str, str]) -> Path:
    dataset_dir = tmp_dir / "build" / "datasets" / "01_digester2020"
    dataset_dir.mkdir(parents=True)
    fasta_text = "".join(f">{seq_id}\n{sequence}\n" for seq_id, sequence in records.items())
    (dataset_dir / "prepared.ssu.fa").write_text(fasta_text, encoding="utf-8")
    return dataset_dir


def _mock_sina_run(
    monkeypatch: pytest.MonkeyPatch,
    output_factory: Callable[[], str],
) -> None:
    def fake_run(command, *args, **kwargs):
        if command[-1] == "--version":
            return subprocess.CompletedProcess(command, 0, stdout="SINA 1.7.2\n", stderr="")
        output_path = Path(command[command.index("-o") + 1])
        output_path.write_text(output_factory(), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="mock sina\n")

    monkeypatch.setattr("autotax2.sina.subprocess.run", fake_run)


def _mock_sina_run_with_candidates(
    monkeypatch: pytest.MonkeyPatch,
    fasta_text: str,
    candidate_csv: str,
) -> None:
    def fake_run(command, *args, **kwargs):
        if command[-1] == "--version":
            return subprocess.CompletedProcess(command, 0, stdout="SINA 1.7.2\n", stderr="")
        output_index = command.index("-o") + 1
        output_path = Path(command[output_index])
        output_path.write_text(fasta_text, encoding="utf-8")
        for value in command[output_index + 1 : command.index("--threads")]:
            extra_output = Path(value)
            if extra_output.suffix == ".csv":
                extra_output.write_text(candidate_csv, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="mock sina\n")

    monkeypatch.setattr("autotax2.sina.subprocess.run", fake_run)


def _mock_sina_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(command, *args, **kwargs):
        if command[-1] == "--version":
            return subprocess.CompletedProcess(command, 0, stdout="SINA 1.7.2\n", stderr="")
        raise subprocess.CalledProcessError(2, command, stderr="boom")

    monkeypatch.setattr("autotax2.sina.subprocess.run", fake_run)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
