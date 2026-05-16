from __future__ import annotations

import csv
import hashlib
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
from typer.testing import CliRunner

from autotax2.cli import app
from autotax2.silva import (
    ResolveParentJob,
    _candidate_anchor_identities,
    _make_resolve_anchor,
    _parent_job_batches,
)


runner = CliRunner()
FIXTURES = Path("tests") / "fixtures"
GTDB_AR53 = FIXTURES / "gtdb_ar53_taxonomy_r232.tsv"
GTDB_BAC120 = FIXTURES / "gtdb_bac120_taxonomy_r232.tsv"


@pytest.fixture
def resolve_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_resolve_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_resolve_silva_creates_placeholder_framework(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    named_taxa_before = _named_taxon_signature(build)
    _write_cluster_uc(build)

    result = runner.invoke(app, ["resolve", "--build", str(build), "--threads", "8"])

    assert result.exit_code == 0, result.output
    taxa_rows = _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv")
    member_rows = _read_tsv(build / "silva" / "silva_unresolved_members.tsv")
    mapping_rows = _read_tsv(build / "silva" / "silva_unresolved_mapping.tsv")
    evidence_rows = _read_tsv(build / "silva" / "silva_unresolved_evidence.tsv")

    names = [row["name"] for row in taxa_rows]
    assert "g__SILVAg000001" in names
    assert "g__SILVAg000002" in names
    assert "s__SILVAs000001" in names
    assert len(names) == len(set(names))

    members = {row["seq_id"]: row for row in member_rows}
    assert members["U2"]["genus_placeholder"] == "g__SILVAg000001"
    assert members["U2"]["species_placeholder"] == "s__SILVAs000001"
    assert members["U3"]["genus_placeholder"] == "g__SILVAg000002"
    assert members["U1"]["genus_placeholder"] == ""
    assert members["U1"]["species_placeholder"] == "s__SILVAs000003"
    assert len(members["U2"]["resolved_taxonomy"].split(";")) == 7
    assert "g__SILVAg000001" in members["U2"]["resolved_taxonomy"]
    assert "s__SILVAs000001" in members["U2"]["resolved_taxonomy"]
    assert "g__Methanobacterium" in members["U1"]["resolved_taxonomy"]
    assert "s__SILVAs000003" in members["U1"]["resolved_taxonomy"]

    assert members["U2"]["warning"] == ""
    assert members["U3"]["warning"] == ""
    assert mapping_rows
    assert all("unidentified" not in row["placeholder_taxa"] for row in mapping_rows)
    assert _named_taxon_signature(build) == named_taxa_before
    assert evidence_rows
    assert {
        row["decision"]
        for row in evidence_rows
        if row["seq_id"] == "U2" and row["rank"] == "genus"
    } == {"create_placeholder"}
    assert all(row["candidate_scope"] == "same-parent only" for row in evidence_rows)

    content = (build / "silva" / "silva_unresolved_taxa.tsv").read_text(encoding="utf-8")
    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content


def test_resolve_silva_rerun_reuses_existing_cluster_taxa(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)

    first = runner.invoke(app, ["resolve", "--build", str(build)])
    second = runner.invoke(app, ["resolve", "--build", str(build)])

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
    assert counters["genus"] == "3"
    assert counters["species"] == "4"
    yaml_text = (build / "registry" / "placeholder_counters.yaml").read_text(encoding="utf-8")
    assert "SILVA:" in yaml_text
    assert "genus: 3" in yaml_text
    assert "species: 4" in yaml_text


def test_resolve_silva_does_not_reuse_deprecated_placeholder(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)
    _append_deprecated_taxon(build, "g__SILVAg000001", "genus")

    result = runner.invoke(app, ["resolve", "--build", str(build)])

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

    result = runner.invoke(app, ["resolve", "--build", str(build), "--dry-run"])

    assert result.exit_code == 0, result.output
    assert counters_path.read_text(encoding="utf-8") == before
    assert (build / "silva" / "silva_unresolved_taxa.dry_run.tsv").exists()
    assert _read_tsv(build / "registry" / "cluster_to_taxon.tsv") == []


def test_resolve_silva_accepts_all_rank_threshold_overrides(resolve_tmp_dir: Path) -> None:
    build = _init_resolve_fixture(resolve_tmp_dir)
    _write_cluster_uc(build)

    result = runner.invoke(app, ["resolve", "--build", str(build), "--family-id", "0.900"])

    assert result.exit_code == 0, result.output


def test_resolve_silva_assigns_unknown_rank_to_same_parent_anchor_only(
    resolve_tmp_dir: Path,
) -> None:
    build = resolve_tmp_dir / "build"
    fasta = resolve_tmp_dir / "silva.fa"
    metadata = resolve_tmp_dir / "metadata.tsv"
    fasta.write_text(
        ">N1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGTACGT\n"
        ">U_EURY Archaea;Euryarchaeota;unidentified;unidentified;"
        "unidentified;unidentified;unidentified\n"
        "ACGTACGA\n"
        ">U_ASGARD Archaea;Asgardarchaeota;unidentified;unidentified;"
        "unidentified;unidentified;unidentified\n"
        "ACGTACGA\n",
        encoding="utf-8",
    )
    metadata.write_text(
        "primaryAccession\torganism_name\tstrain\ttype_material\n"
        "N1\tMethanobacterium formicicum\tDSM 1535\ttype strain\n",
        encoding="utf-8",
    )
    init_result = runner.invoke(app, _init_args(fasta, build, metadata))
    assert init_result.exit_code == 0, init_result.output

    class_uc = _parent_job_uc_path(
        build,
        "class",
        0.722,
        ("d__Archaea", "p__Euryarchaeota"),
    )
    class_uc.parent.mkdir(parents=True, exist_ok=True)
    class_uc.write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tANCHOR_000001\t*\n"
        "H\t0\t*\t99.0\t+\t0\t0\t0\tU_EURY\tANCHOR_000001\n",
        encoding="utf-8",
    )

    result = runner.invoke(app, ["resolve", "--build", str(build)])

    assert result.exit_code == 0, result.output
    members = {row["seq_id"]: row for row in _read_tsv(build / "silva" / "silva_unresolved_members.tsv")}
    evidence = {
        (row["seq_id"], row["rank"]): row
        for row in _read_tsv(build / "silva" / "silva_unresolved_evidence.tsv")
    }

    assert "c__Methanobacteria" in members["U_EURY"]["resolved_taxonomy"]
    assert evidence[("U_EURY", "class")]["decision"] == "assign_existing"
    assert evidence[("U_EURY", "class")]["output_taxon"] == "c__Methanobacteria"
    assert evidence[("U_EURY", "class")]["parent_taxon"] == "p__Euryarchaeota"
    assert evidence[("U_EURY", "class")]["best_anchor_identity"] == "0.990"

    assert evidence[("U_ASGARD", "class")]["decision"] == "create_placeholder"
    assert evidence[("U_ASGARD", "class")]["parent_taxon"] == "p__Asgardarchaeota"
    assert "c__Methanobacteria" not in members["U_ASGARD"]["resolved_taxonomy"]


def test_parent_job_batches_respect_global_thread_cap(resolve_tmp_dir: Path) -> None:
    jobs = [
        ResolveParentJob(
            rank="genus",
            parent_key=(f"p__P{index}",),
            candidates=[],
            anchors=[],
            threshold=0.901,
            cluster_dir=resolve_tmp_dir,
            job_threads=threads,
            vsearch_bin="vsearch",
            iddef=2,
        )
        for index, threads in enumerate([3, 2, 1, 4, 1])
    ]

    batches = _parent_job_batches(jobs, total_threads=4)

    assert [[job.job_threads for job in batch] for batch in batches] == [[3], [2, 1], [4], [1]]
    assert all(sum(job.job_threads for job in batch) <= 4 for batch in batches)


def test_candidate_anchor_identities_use_uc_hit_identity(resolve_tmp_dir: Path) -> None:
    anchor = _make_resolve_anchor(
        rank="class",
        taxon_name="Methanobacteria",
        seq_id="N1",
        sequence="ACGT",
        source="named_silva",
        index=1,
    )
    uc_path = resolve_tmp_dir / "anchor_identity.uc"
    uc_path.write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tANCHOR_000001\t*\n"
        "H\t0\t*\t99.0\t+\t0\t0\t0\tU_EURY\tANCHOR_000001\n",
        encoding="utf-8",
    )

    identities = _candidate_anchor_identities(
        uc_path=uc_path,
        candidate_ids={"U_EURY"},
        job_id_to_anchor={anchor.job_id: anchor},
    )

    assert identities[("U_EURY", "c__Methanobacteria")] == "0.990"


def test_resolve_silva_with_no_unresolved_records_exits_cleanly(resolve_tmp_dir: Path) -> None:
    build = resolve_tmp_dir / "build"
    fasta = resolve_tmp_dir / "named.fa"
    metadata = resolve_tmp_dir / "metadata.tsv"
    fasta.write_text(
        ">N1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata)
    init_result = runner.invoke(
        app,
        _init_args(fasta, build, metadata),
    )

    result = runner.invoke(app, ["resolve", "--build", str(build)])

    assert init_result.exit_code == 0, init_result.output
    assert result.exit_code == 0, result.output
    assert _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv") == []
    assert _read_tsv(build / "silva" / "silva_unresolved_members.tsv") == []
    assert _read_tsv(build / "silva" / "silva_unresolved_mapping.tsv") == []


def _init_resolve_fixture(tmp_dir: Path) -> Path:
    build = tmp_dir / "build"
    fasta = tmp_dir / "silva.fa"
    metadata = tmp_dir / "metadata.tsv"
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
    _write_empty_metadata(metadata)
    result = runner.invoke(
        app,
        _init_args(fasta, build, metadata),
    )
    assert result.exit_code == 0, result.output
    return build


def _write_empty_metadata(path: Path) -> None:
    path.write_text(
        "primaryAccession\torganism_name\tstrain\ttype_material\n",
        encoding="utf-8",
    )


def _init_args(fasta: Path, build: Path, metadata: Path) -> list[str]:
    return [
        "init",
        "--silva-fasta",
        str(fasta),
        "--outdir",
        str(build),
        "--type-strain-metadata",
        str(metadata),
        "--gtdb-ar53-taxonomy",
        str(GTDB_AR53),
        "--gtdb-bac120-taxonomy",
        str(GTDB_BAC120),
    ]


def _write_cluster_uc(build: Path) -> None:
    cluster_dir = build / "silva" / "silva_unresolved_clusters"
    cluster_dir.mkdir(parents=True, exist_ok=True)
    (cluster_dir / "genus_0.901.uc").write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tU2\t*\n"
        "H\t0\t*\t99.0\t+\t0\t0\t0\tU3\tU2\n"
        "S\t1\t*\t*\t*\t*\t*\t*\tU1\t*\n",
        encoding="utf-8",
    )
    (cluster_dir / "species_0.972.uc").write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tU2\t*\n"
        "H\t0\t*\t99.0\t+\t0\t0\t0\tU3\tU2\n"
        "S\t1\t*\t*\t*\t*\t*\t*\tU1\t*\n",
        encoding="utf-8",
    )


def _parent_job_uc_path(
    build: Path,
    rank: str,
    threshold: float,
    parent_key: tuple[str, ...],
) -> Path:
    parent_hash = hashlib.sha1(";".join(parent_key).encode("utf-8")).hexdigest()[:12]
    return build / "silva" / "silva_unresolved_clusters" / rank / f"{rank}_{threshold:.3f}_{parent_hash}.uc"


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
