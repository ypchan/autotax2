from __future__ import annotations

import csv
import shutil
import subprocess
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest

from autotax2.vsearch import (
    build_current_representatives,
    build_cluster_fast_command,
    build_usearch_global_command,
    check_vsearch_version,
    cluster_memberships_from_uc,
    cluster_search_dataset,
    filter_registry_hits,
    parse_registry_hits,
    parse_uc_records,
)


@pytest.fixture
def vsearch_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_vsearch_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_cluster_fast_command_includes_iddef_2() -> None:
    command = build_cluster_fast_command(
        input_fasta="in.fa",
        uc_path="out.uc",
        centroids_path="centroids.fa",
        identity=0.987,
        threads=16,
        vsearch_bin="vsearch",
        iddef=2,
    )

    assert "--iddef" in command
    assert command[command.index("--iddef") + 1] == "2"
    assert command[command.index("--uc") + 1] == "out.uc"
    assert command[command.index("--centroids") + 1] == "centroids.fa"


def test_usearch_global_command_includes_search_controls() -> None:
    command = build_usearch_global_command(
        query_fasta="query.fa",
        db_fasta="db.fa",
        userout_path="hits.tsv",
        identity=0.750,
        maxaccepts=50,
        maxrejects=256,
        strand="both",
    )

    assert "--maxaccepts" in command
    assert command[command.index("--maxaccepts") + 1] == "50"
    assert "--maxrejects" in command
    assert command[command.index("--maxrejects") + 1] == "256"
    assert "--strand" in command
    assert command[command.index("--strand") + 1] == "both"


def test_vsearch_version_parser_handles_typical_output(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(command, *args, **kwargs):
        return subprocess.CompletedProcess(command, 0, stdout="vsearch v2.29.3_linux_x86_64\n", stderr="")

    monkeypatch.setattr("autotax2.vsearch.subprocess.run", fake_run)

    assert check_vsearch_version("vsearch") == "2.29.3"


def test_uc_parser_identifies_centroid_and_hit_records(vsearch_tmp_dir: Path) -> None:
    uc_path = vsearch_tmp_dir / "species.uc"
    uc_path.write_text(
        "S\t0\t100\t*\t*\t*\t*\t*\tD20_000001\t*\n"
        "H\t0\t100\t99.1\t+\t0\t0\t100M\tD20_000002\tD20_000001\n",
        encoding="utf-8",
    )

    records = parse_uc_records(uc_path)
    memberships = cluster_memberships_from_uc(uc_path, "species_0.987")

    assert records[0].record_type == "S"
    assert records[1].record_type == "H"
    assert memberships[0].centroid_id == "D20_000001"
    assert memberships[0].member_id == "D20_000001"
    assert memberships[1].centroid_id == "D20_000001"
    assert memberships[1].member_id == "D20_000002"
    assert memberships[1].identity_to_centroid == "99.1"


def test_registry_hit_parser_computes_coverage(vsearch_tmp_dir: Path) -> None:
    hits_path = vsearch_tmp_dir / "hits.tsv"
    hits_path.write_text(
        "D20_000001\tREF1\t98.7\t80\t1\t80\t5\t84\t100\t160\t42\n",
        encoding="utf-8",
    )

    hits = parse_registry_hits(hits_path)

    assert hits[0].query_coverage == pytest.approx(0.8)
    assert hits[0].target_coverage == pytest.approx(0.5)


def test_registry_hit_filter_keeps_multiple_passing_hits(vsearch_tmp_dir: Path) -> None:
    hits_path = vsearch_tmp_dir / "hits.tsv"
    hits_path.write_text(
        "D20_000001\tREF1\t98.7\t90\t1\t90\t5\t94\t100\t160\t42\n"
        "D20_000001\tREF2\t98.3\t90\t1\t90\t5\t94\t100\t160\t41\n"
        "D20_000001\tREF3\t99.0\t60\t1\t60\t5\t64\t100\t160\t45\n",
        encoding="utf-8",
    )

    filtered = filter_registry_hits(
        parse_registry_hits(hits_path),
        floor_id=0.750,
        min_query_cov=0.80,
        min_target_cov=0.0,
    )

    assert [hit.target for hit in filtered] == ["REF1", "REF2"]


def test_cluster_search_writes_outputs_and_summary(
    vsearch_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    build = _make_cluster_search_build(vsearch_tmp_dir)
    _mock_vsearch_run(monkeypatch)

    summary = cluster_search_dataset(
        build=build,
        dataset="digester2020",
        threads=16,
        vsearch_bin="vsearch",
    )

    dataset_dir = build / "datasets" / "01_digester2020"
    cluster_dir = dataset_dir / "internal_clusters"
    filtered_content = (dataset_dir / "vs_registry.filtered.tsv").read_text(encoding="utf-8")
    filtered_rows = _read_tsv(dataset_dir / "vs_registry.filtered.tsv")
    summary_rows = _read_tsv(dataset_dir / "cluster_search_summary.tsv")
    tool_rows = _read_tsv(dataset_dir / "tool_versions.tsv")

    assert summary.registry_hits_filtered == 2
    assert (cluster_dir / "species_0.987.uc").exists()
    assert (cluster_dir / "species_0.987.centroids.fa").exists()
    assert (cluster_dir / "species_0.987.members.tsv").exists()
    assert [row["target"] for row in filtered_rows] == ["REF1", "REF2"]
    assert "\t" in filtered_content
    assert "\n" in filtered_content
    assert "\\t" not in filtered_content
    assert "\\n" not in filtered_content
    assert summary_rows[0]["dataset"] == "digester2020"
    assert summary_rows[0]["prefix"] == "D20"
    assert summary_rows[0]["input_sequences"] == "2"
    assert summary_rows[0]["species_centroids"] == "1"
    assert summary_rows[0]["registry_representatives"] == "2"
    assert summary_rows[0]["registry_hits_raw"] == "3"
    assert summary_rows[0]["registry_hits_filtered"] == "2"
    assert summary_rows[0]["iddef"] == "2"
    assert tool_rows[0]["tool"] == "vsearch"
    assert "--iddef 2" in tool_rows[0]["command"]
    assert "--maxaccepts 50" in tool_rows[0]["command"]


def test_build_current_representatives_rebuilds_existing_cache(vsearch_tmp_dir: Path) -> None:
    build = vsearch_tmp_dir / "build"
    registry = build / "registry"
    silva = build / "silva"
    registry.mkdir(parents=True)
    silva.mkdir(parents=True)
    (registry / "current_representatives.fa").write_text(">OLD\nAAAA\n", encoding="utf-8")
    (silva / "silva_named_backbone.fa").write_text(
        ">REF_NEW Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )

    count = build_current_representatives(build)
    content = (registry / "current_representatives.fa").read_text(encoding="utf-8")

    assert count == 1
    assert ">REF_NEW" in content
    assert ">OLD" not in content


def _make_cluster_search_build(tmp_dir: Path) -> Path:
    build = tmp_dir / "build"
    dataset_dir = build / "datasets" / "01_digester2020"
    registry_dir = build / "registry"
    dataset_dir.mkdir(parents=True)
    registry_dir.mkdir(parents=True)
    (dataset_dir / "sina.oriented.fa").write_text(
        ">D20_000001\nACGTACGTACGT\n"
        ">D20_000002\nACGTACGTACGA\n",
        encoding="utf-8",
    )
    (registry_dir / "dataset_registry.tsv").write_text(
        "dataset_name\tprefix\tadd_order\tinput_fasta\tinput_md5\tdomain\tstatus\tdataset_dir\n"
        "digester2020\tD20\t1\tinput.fa\tabc\tArchaea\tprepared\tdatasets/01_digester2020\n",
        encoding="utf-8",
    )
    (registry_dir / "current_representatives.fa").write_text(
        ">REF1\nACGTACGTACGT\n"
        ">REF2\nACGTACGTACGA\n",
        encoding="utf-8",
    )
    return build


def _mock_vsearch_run(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(command, *args, **kwargs):
        if command[-1] == "--version":
            return subprocess.CompletedProcess(command, 0, stdout="vsearch v2.29.3\n", stderr="")
        if "--cluster_fast" in command:
            uc_path = Path(command[command.index("--uc") + 1])
            centroids_path = Path(command[command.index("--centroids") + 1])
            uc_path.write_text(
                "S\t0\t100\t*\t*\t*\t*\t*\tD20_000001\t*\n"
                "H\t0\t100\t99.2\t+\t0\t0\t100M\tD20_000002\tD20_000001\n",
                encoding="utf-8",
            )
            centroids_path.write_text(">D20_000001\nACGTACGTACGT\n", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
        if "--usearch_global" in command:
            userout_path = Path(command[command.index("--userout") + 1])
            userout_path.write_text(
                "D20_000001\tREF1\t98.7\t90\t1\t90\t1\t90\t100\t120\t42\n"
                "D20_000001\tREF2\t98.2\t90\t1\t90\t1\t90\t100\t120\t41\n"
                "D20_000001\tREF3\t99.0\t60\t1\t60\t1\t60\t100\t120\t45\n",
                encoding="utf-8",
            )
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr("autotax2.vsearch.subprocess.run", fake_run)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
