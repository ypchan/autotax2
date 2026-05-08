from __future__ import annotations

import csv
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
import yaml

from autotax2.placement import place_dataset


@pytest.fixture
def placement_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_placement_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_place_known_like_with_stable_species(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    assignments = _read_tsv(build / "datasets" / "01_digester2020" / "assignments.tsv")
    row = _by_id(assignments, "D20_000001")
    assert row["identity_status"] == "known_like"
    assert row["final_status"] == "known_like"
    assert row["assigned_taxon_id"] == "S1"
    assert row["species_consensus"] == "s__Methanobacterium formicicum"


def test_place_new_species_under_stable_genus(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    dataset_dir = build / "datasets" / "01_digester2020"
    assignments = _read_tsv(dataset_dir / "assignments.tsv")
    created = _read_tsv(dataset_dir / "created_taxa.tsv")
    row = _by_id(assignments, "D20_000002")

    assert row["identity_status"] == "new_species"
    assert row["final_status"] == "new_species"
    assert row["assigned_taxon_id"] == "s__D20s000001"
    assert any(taxon["name"] == "s__D20s000001" for taxon in created)


def test_place_genus_unstable_family_stable_becomes_new_genus(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    assignments = _read_tsv(build / "datasets" / "01_digester2020" / "assignments.tsv")
    row = _by_id(assignments, "D20_000003")

    assert row["identity_status"] == "new_species"
    assert row["final_status"] == "new_genus"
    assert row["warning"] == "genus_consensus_unstable"
    assert row["family_consensus_fraction"] == "1.000000"
    assert row["genus_consensus_fraction"] == "0.500000"


def test_place_new_genus_from_family_identity(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    dataset_dir = build / "datasets" / "01_digester2020"
    assignments = _read_tsv(dataset_dir / "assignments.tsv")
    created = _read_tsv(dataset_dir / "created_taxa.tsv")
    row = _by_id(assignments, "D20_000004")

    assert row["identity_status"] == "new_genus"
    assert row["final_status"] == "new_genus"
    assert any(taxon["name"] == "g__D20g000002" for taxon in created)
    assert any(taxon["rank"] == "species" for taxon in created)


def test_place_no_hit_above_floor_is_unplaced(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    assignments = _read_tsv(build / "datasets" / "01_digester2020" / "assignments.tsv")
    row = _by_id(assignments, "D20_000005")
    assert row["identity_status"] == "unplaced"
    assert row["final_status"] == "unplaced"
    assert row["assigned_taxon_id"] == ""


def test_place_duplicate_md5_does_not_create_taxon(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    dataset_dir = build / "datasets" / "01_digester2020"
    assignments = _read_tsv(dataset_dir / "assignments.tsv")
    created = _read_tsv(dataset_dir / "created_taxa.tsv")
    row = _by_id(assignments, "D20_000006")

    assert row["final_status"] == "duplicate"
    assert row["assigned_taxon_id"] == "S1"
    assert "D20_000006" not in {taxon["representative_seq_id"] for taxon in created}


def test_new_genus_creates_unique_genus_and_species_placeholders(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    created = _read_tsv(build / "datasets" / "01_digester2020" / "created_taxa.tsv")
    names = [taxon["name"] for taxon in created]

    assert "g__D20g000001" in names
    assert "s__D20s000002" in names
    assert len(names) == len(set(names))


def test_cluster_key_reuse_returns_same_taxon(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)
    dataset_dir = build / "datasets" / "01_digester2020"

    place_dataset(build, "digester2020")
    first = _by_id(_read_tsv(dataset_dir / "assignments.tsv"), "D20_000002")["assigned_taxon_id"]

    place_dataset(build, "digester2020")
    second = _by_id(_read_tsv(dataset_dir / "assignments.tsv"), "D20_000002")["assigned_taxon_id"]
    cluster_rows = _read_tsv(build / "registry" / "cluster_to_taxon.tsv")
    cluster_keys = [row["cluster_key"] for row in cluster_rows]

    assert second == first == "s__D20s000001"
    assert len(cluster_keys) == len(set(cluster_keys))


def test_silva_named_representative_is_not_replaced(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)
    registry_path = build / "registry" / "representative_registry.tsv"
    before = _read_tsv(registry_path)

    place_dataset(build, "digester2020")

    after = _read_tsv(registry_path)
    ref_a_before = [row for row in before if row["representative_seq_id"] == "REF_A"]
    ref_a_after = [row for row in after if row["representative_seq_id"] == "REF_A"]

    assert len(ref_a_after) == len(ref_a_before) == 1
    assert ref_a_after[0]["representative_seq_id"] == ref_a_before[0]["representative_seq_id"]
    assert ref_a_after[0]["taxon_id"] == ref_a_before[0]["taxon_id"]
    assert ref_a_after[0]["status"] == ref_a_before[0]["status"]


def test_assignments_tsv_has_real_tabs_and_newlines(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)

    place_dataset(build, "digester2020")

    content = (build / "datasets" / "01_digester2020" / "assignments.tsv").read_text(
        encoding="utf-8"
    )
    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content


def test_dry_run_does_not_update_counters(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)
    counters_path = build / "registry" / "placeholder_counters.yaml"
    before = counters_path.read_text(encoding="utf-8")

    place_dataset(build, "digester2020", dry_run=True)

    after = counters_path.read_text(encoding="utf-8")
    dataset_dir = build / "datasets" / "01_digester2020"

    assert after == before
    assert (dataset_dir / "assignments.dry_run.tsv").exists()
    assert not (dataset_dir / "assignments.tsv").exists()


def test_named_silva_backbone_rows_remain_protected(placement_tmp_dir: Path) -> None:
    build = _make_place_build(placement_tmp_dir)
    taxon_path = build / "registry" / "taxon_nodes.tsv"
    before = _named_taxon_signature(_read_tsv(taxon_path))

    place_dataset(build, "digester2020")

    after = _named_taxon_signature(_read_tsv(taxon_path))
    assert after == before


def _make_place_build(tmp_path: Path) -> Path:
    build = tmp_path / "build"
    registry_dir = build / "registry"
    dataset_dir = build / "datasets" / "01_digester2020"
    cluster_dir = dataset_dir / "internal_clusters"
    registry_dir.mkdir(parents=True)
    cluster_dir.mkdir(parents=True)

    _write_tsv(
        registry_dir / "taxon_nodes.tsv",
        [
            _taxon("D1", "domain", "d__Archaea", ""),
            _taxon("P1", "phylum", "p__Euryarchaeota", "D1"),
            _taxon("C1", "class", "c__Methanobacteria", "P1"),
            _taxon("O1", "order", "o__Methanobacteriales", "C1"),
            _taxon("F1", "family", "f__Methanobacteriaceae", "O1"),
            _taxon("G1", "genus", "g__Methanobacterium", "F1"),
            _taxon("S1", "species", "s__Methanobacterium formicicum", "G1"),
            _taxon("G2", "genus", "g__Methanobrevibacter", "F1"),
            _taxon("S2", "species", "s__Methanobrevibacter smithii", "G2"),
        ],
        [
            "taxon_id",
            "rank",
            "name",
            "parent_taxon_id",
            "protected",
            "is_silva_named",
            "source",
        ],
    )
    _write_tsv(
        registry_dir / "representative_registry.tsv",
        [
            {"representative_seq_id": "REF_A", "taxon_id": "S1", "status": "active"},
            {"representative_seq_id": "REF_A2", "taxon_id": "S1", "status": "active"},
            {"representative_seq_id": "REF_B", "taxon_id": "S2", "status": "active"},
        ],
        ["representative_seq_id", "taxon_id", "status"],
    )
    _write_tsv(
        registry_dir / "sequence_registry.tsv",
        [
            {
                "seq_id": "OLD_DUP",
                "dataset": "older",
                "sequence_md5": "md5dup",
                "taxon_id": "S1",
            }
        ],
        ["seq_id", "dataset", "sequence_md5", "taxon_id"],
    )
    _write_tsv(
        registry_dir / "dataset_registry.tsv",
        [
            {
                "dataset_name": "digester2020",
                "prefix": "D20",
                "add_order": "1",
                "input_fasta": "input.fa",
                "input_md5": "abc",
                "domain": "Archaea",
                "status": "prepared",
                "dataset_dir": str(dataset_dir),
            }
        ],
        [
            "dataset_name",
            "prefix",
            "add_order",
            "input_fasta",
            "input_md5",
            "domain",
            "status",
            "dataset_dir",
        ],
    )
    (registry_dir / "placeholder_counters.yaml").write_text(
        yaml.safe_dump(
            {
                "D20": {
                    "class": 1,
                    "order": 1,
                    "family": 1,
                    "genus": 1,
                    "species": 1,
                }
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    _write_tsv(
        registry_dir / "cluster_to_taxon.tsv",
        [],
        ["cluster_key", "taxon_id", "rank", "name", "status", "source_prefix"],
    )

    memberships = [
        _membership("D20_000001", "orig1", "md501", "false"),
        _membership("D20_000002", "orig2", "md502", "false"),
        _membership("D20_000003", "orig3", "md503", "false"),
        _membership("D20_000004", "orig4", "md504", "false"),
        _membership("D20_000005", "orig5", "md505", "false"),
        _membership("D20_000006", "orig6", "md5dup", "false"),
    ]
    _write_tsv(
        dataset_dir / "sequence_membership.tsv",
        memberships,
        [
            "internal_seq_id",
            "original_seq_id",
            "dataset",
            "prefix",
            "sequence_md5",
            "unique_seq_id",
            "is_duplicate_sequence",
        ],
    )
    _write_tsv(
        dataset_dir / "sequence_id_map.tsv",
        [
            {
                "internal_seq_id": row["internal_seq_id"],
                "original_seq_id": row["original_seq_id"],
                "sequence_md5": row["sequence_md5"],
            }
            for row in memberships
        ],
        ["internal_seq_id", "original_seq_id", "sequence_md5"],
    )
    _write_tsv(
        dataset_dir / "cluster_search_summary.tsv",
        [{"dataset": "digester2020", "prefix": "D20"}],
        ["dataset", "prefix"],
    )
    (cluster_dir / "species_0.987.uc").write_text(
        "".join(
            f"S\t{index}\t100\t*\t*\t*\t*\t*\tD20_00000{index + 1}\t*\n"
            for index in range(6)
        ),
        encoding="utf-8",
    )
    _write_tsv(
        dataset_dir / "vs_registry.filtered.tsv",
        [
            {"query": "D20_000001", "target": "REF_A", "identity": "0.990"},
            {"query": "D20_000002", "target": "REF_A", "identity": "0.960"},
            {"query": "D20_000003", "target": "REF_A", "identity": "0.960"},
            {"query": "D20_000003", "target": "REF_B", "identity": "0.959"},
            {"query": "D20_000004", "target": "REF_A", "identity": "0.900"},
            {"query": "D20_000005", "target": "REF_A", "identity": "0.700"},
            {"query": "D20_000006", "target": "REF_A", "identity": "0.990"},
        ],
        ["query", "target", "identity"],
    )
    return build


def _taxon(taxon_id: str, rank: str, name: str, parent: str) -> dict[str, str]:
    return {
        "taxon_id": taxon_id,
        "rank": rank,
        "name": name,
        "parent_taxon_id": parent,
        "protected": "true",
        "is_silva_named": "true",
        "source": "SILVA138.2_NR99",
    }


def _membership(
    internal_seq_id: str,
    original_seq_id: str,
    sequence_md5: str,
    is_duplicate: str,
) -> dict[str, str]:
    return {
        "internal_seq_id": internal_seq_id,
        "original_seq_id": original_seq_id,
        "dataset": "digester2020",
        "prefix": "D20",
        "sequence_md5": sequence_md5,
        "unique_seq_id": f"U_{internal_seq_id[-6:]}",
        "is_duplicate_sequence": is_duplicate,
    }


def _by_id(rows: list[dict[str, str]], internal_seq_id: str) -> dict[str, str]:
    return next(row for row in rows if row["internal_seq_id"] == internal_seq_id)


def _named_taxon_signature(rows: list[dict[str, str]]) -> list[tuple[str, str, str, str, str]]:
    return [
        (
            row["taxon_id"],
            row["rank"],
            row["name"],
            row["parent_taxon_id"],
            row["protected"],
        )
        for row in rows
        if row.get("is_silva_named") == "true"
    ]


def _write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
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
