from __future__ import annotations

import csv
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
from typer.testing import CliRunner

from autotax2.cli import app


runner = CliRunner()


@pytest.fixture
def resolve_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_resolve_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_resolve_silva_creates_placeholder_scaffold(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    named_taxa_before = _named_taxon_signature(build)
    _write_cluster_uc(build)

    result = runner.invoke(app, ["resolve-silva", "--build", str(build), "--threads", "8"])

    assert result.exit_code == 0, result.output
    taxa_rows = _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv")
    member_rows = _read_tsv(build / "silva" / "silva_unresolved_members.tsv")
    mapping_rows = _read_tsv(build / "silva" / "silva_unresolved_mapping.tsv")

    names = [row["name"] for row in taxa_rows]
    assert "g__SILVAg000001" in names
    assert "s__SILVAs000001" in names
    assert len(names) == len(set(names))

    members = {row["seq_id"]: row for row in member_rows}
    assert members["U2"]["genus_placeholder"] == "g__SILVAg000001"
    assert members["U2"]["species_placeholder"] == "s__SILVAs000001"
    assert members["U1"]["genus_placeholder"] == ""
    assert members["U1"]["species_placeholder"] == "s__SILVAs000002"
    assert len(members["U2"]["resolved_taxonomy"].split(";")) == 7
    assert "g__SILVAg000001" in members["U2"]["resolved_taxonomy"]
    assert "s__SILVAs000001" in members["U2"]["resolved_taxonomy"]
    assert "g__Methanobacterium" in members["U1"]["resolved_taxonomy"]
    assert "s__SILVAs000002" in members["U1"]["resolved_taxonomy"]

    assert members["U2"]["warning"] == "mixed_silva_parent_unresolved_cluster"
    assert members["U3"]["warning"] == "mixed_silva_parent_unresolved_cluster"
    assert mapping_rows
    assert all("unidentified" not in row["placeholder_taxa"] for row in mapping_rows)
    assert _named_taxon_signature(build) == named_taxa_before

    content = (build / "silva" / "silva_unresolved_taxa.tsv").read_text(encoding="utf-8")
    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content


def test_resolve_silva_rerun_reuses_existing_cluster_taxa(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)

    first = runner.invoke(app, ["resolve-silva", "--build", str(build)])
    second = runner.invoke(app, ["resolve-silva", "--build", str(build)])

    assert first.exit_code == 0, first.output
    assert second.exit_code == 0, second.output
    taxa_rows = _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv")
    cluster_rows = _read_tsv(build / "registry" / "cluster_to_taxon.tsv")

    assert [row["name"] for row in taxa_rows].count("g__SILVAg000001") == 1
    assert [row["name"] for row in cluster_rows].count("g__SILVAg000001") == 1
    counters = {
        row["rank"]: row["next_ordinal"]
        for row in _read_tsv(build / "registry" / "placeholder_counters.tsv")
    }
    assert counters["genus"] == "2"
    assert counters["species"] == "3"
    yaml_text = (build / "registry" / "placeholder_counters.yaml").read_text(encoding="utf-8")
    assert "SILVA:" in yaml_text
    assert "genus: 2" in yaml_text
    assert "species: 3" in yaml_text


def test_resolve_silva_does_not_reuse_deprecated_placeholder(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)
    _append_deprecated_taxon(build, "g__SILVAg000001", "genus")

    result = runner.invoke(app, ["resolve-silva", "--build", str(build)])

    assert result.exit_code == 0, result.output
    taxa_rows = _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv")
    names = [row["name"] for row in taxa_rows]

    assert "g__SILVAg000001" not in names
    assert "g__SILVAg000002" in names


def test_resolve_silva_dry_run_does_not_update_counters(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)
    counters_path = build / "registry" / "placeholder_counters.tsv"
    counters_path.write_text(
        "source_prefix\trank\tnext_ordinal\n"
        "SILVA\tclass\t1\n"
        "SILVA\torder\t1\n"
        "SILVA\tfamily\t1\n"
        "SILVA\tgenus\t7\n"
        "SILVA\tspecies\t9\n",
        encoding="utf-8",
    )
    before = counters_path.read_text(encoding="utf-8")

    result = runner.invoke(app, ["resolve-silva", "--build", str(build), "--dry-run"])

    assert result.exit_code == 0, result.output
    assert counters_path.read_text(encoding="utf-8") == before
    assert (build / "silva" / "silva_unresolved_taxa.dry_run.tsv").exists()
    assert _read_tsv(build / "registry" / "cluster_to_taxon.tsv") == []


def test_resolve_silva_rejects_unimplemented_rank_threshold_overrides(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)

    result = runner.invoke(app, ["resolve-silva", "--build", str(build), "--family-id", "0.900"])

    assert result.exit_code != 0
    assert "currently clusters unresolved SILVA records at genus and species" in str(result.exception)


def test_resolve_silva_with_no_unresolved_records_exits_cleanly(resolve_tmp_dir: Path) -> None:
    build = resolve_tmp_dir / "build"
    fasta = resolve_tmp_dir / "named.fa"
    fasta.write_text(
        ">N1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )
    init_result = runner.invoke(
        app,
        ["init", "--silva-fasta", str(fasta), "--outdir", str(build)],
    )

    result = runner.invoke(app, ["resolve-silva", "--build", str(build)])

    assert init_result.exit_code == 0, init_result.output
    assert result.exit_code == 0, result.output
    assert _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv") == []
    assert _read_tsv(build / "silva" / "silva_unresolved_members.tsv") == []
    assert _read_tsv(build / "silva" / "silva_unresolved_mapping.tsv") == []


def _init_resolve_fixture(tmp_dir: Path) -> Path:
    build = tmp_dir / "build"
    fasta = tmp_dir / "silva.fa"
    fasta.write_text(
        ">U2 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;unidentified;unidentified\n"
        "ACGTACGT\n"
        ">U3 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Otheraceae;unidentified;unidentified\n"
        "ACGTACGA\n"
        ">U1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium sp.\n"
        "ACGTACGC\n"
        ">N1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGTACGG\n",
        encoding="utf-8",
    )
    result = runner.invoke(
        app,
        ["init", "--silva-fasta", str(fasta), "--outdir", str(build), "--domain", "Archaea"],
    )
    assert result.exit_code == 0, result.output
    return build


def _write_cluster_uc(build: Path) -> None:
    cluster_dir = build / "silva" / "silva_unresolved_clusters"
    cluster_dir.mkdir(parents=True, exist_ok=True)
    (cluster_dir / "genus_0.945.uc").write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tU2\t*\n"
        "H\t0\t*\t99.0\t+\t0\t0\t0\tU3\tU2\n"
        "S\t1\t*\t*\t*\t*\t*\t*\tU1\t*\n",
        encoding="utf-8",
    )
    (cluster_dir / "species_0.987.uc").write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tU2\t*\n"
        "H\t0\t*\t99.0\t+\t0\t0\t0\tU3\tU2\n"
        "S\t1\t*\t*\t*\t*\t*\t*\tU1\t*\n",
        encoding="utf-8",
    )


def _named_taxon_signature(build: Path) -> list[tuple[str, str, str, str]]:
    return [
        (row["taxon_id"], row["rank"], row["name"], row["source"])
        for row in _read_tsv(build / "registry" / "taxon_nodes.tsv")
        if row.get("is_silva_named") == "true"
    ]


def _append_deprecated_taxon(build: Path, name: str, rank: str) -> None:
    path = build / "registry" / "taxon_nodes.tsv"
    rows = _read_tsv(path)
    fieldnames = list(rows[0].keys()) if rows else ["taxon_id", "rank", "name"]
    for field in ["status", "source_prefix"]:
        if field not in fieldnames:
            fieldnames.append(field)
    rows.append(
        {
            "taxon_id": name,
            "rank": rank,
            "name": name,
            "status": "deprecated",
            "source_prefix": "SILVA",
        }
    )
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


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
