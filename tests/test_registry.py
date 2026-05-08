from __future__ import annotations

import pytest

from autotax2.io import FastaRecord
from autotax2.registry import (
    DatasetPrefixRegistry,
    Registry,
    SequenceRegistry,
    assign_internal_ids,
    remap_sequence_id,
    sequence_md5,
)


def test_remap_sequence_id() -> None:
    assert remap_sequence_id("D20", 1) == "D20_000001"


def test_sequence_md5_normalizes_case_and_whitespace() -> None:
    assert sequence_md5("ac gu") == sequence_md5("ACGU")


def test_sequence_md5_normalizes_u_to_t() -> None:
    assert sequence_md5("ACGU") == sequence_md5("ACGT")


def test_registry_exports_exact_duplicate_sequences_once() -> None:
    registry = Registry()

    first = registry.add_sequence("D20", 1, "seq-a", "ACGU")
    second = registry.add_sequence("D20", 2, "seq-b", "ac gu")

    assert first.sequence_md5 == second.sequence_md5
    assert registry.exported_internal_ids() == ["D20_000001"]


def test_dataset_prefix_registry_freezes_prefix_assignment() -> None:
    prefixes = DatasetPrefixRegistry()

    record = prefixes.register("D20", "dataset-20")

    assert record.prefix == "D20"
    assert record.frozen is True
    assert prefixes.register("D20", "dataset-20") is record


def test_dataset_prefix_registry_rejects_reassignment() -> None:
    prefixes = DatasetPrefixRegistry()
    prefixes.register("D20", "dataset-20")

    with pytest.raises(ValueError):
        prefixes.register("D20", "other-dataset")


def test_registry_auto_assigns_sequence_ordinals() -> None:
    registry = Registry()

    first = registry.add_sequence_auto("D20", "seq-a", "ACGU")
    second = registry.add_sequence_auto("D20", "seq-b", "UGCA")

    assert first.internal_id == "D20_000001"
    assert second.internal_id == "D20_000002"


def test_registry_rejects_duplicate_internal_id() -> None:
    registry = Registry()
    registry.add_sequence("D20", 1, "seq-a", "ACGU")

    with pytest.raises(ValueError):
        registry.add_sequence("D20", 1, "seq-b", "UGCA")


def test_assign_internal_ids_generates_stable_ids() -> None:
    assigned = assign_internal_ids(
        [
            FastaRecord(seq_id="seq-a", header="seq-a first", sequence="ACGT"),
            FastaRecord(seq_id="seq-b", header="seq-b second", sequence="TGCA"),
        ],
        dataset_name="dataset-20",
        prefix="D20",
    )

    assert [record.internal_seq_id for record in assigned] == [
        "D20_000001",
        "D20_000002",
    ]
    assert assigned[0].original_seq_id == "seq-a"
    assert assigned[0].original_header == "seq-a first"
    assert assigned[0].dataset == "dataset-20"
    assert assigned[0].prefix == "D20"
    assert assigned[0].sequence_length == 4


def test_duplicate_original_fasta_ids_get_unique_internal_ids() -> None:
    assigned = assign_internal_ids(
        [
            FastaRecord(seq_id="seq-a", header="seq-a first", sequence="ACGT"),
            FastaRecord(seq_id="seq-a", header="seq-a second", sequence="TGCA"),
        ],
        dataset_name="dataset-20",
        prefix="D20",
    )

    assert [record.original_seq_id for record in assigned] == ["seq-a", "seq-a"]
    assert [record.internal_seq_id for record in assigned] == [
        "D20_000001",
        "D20_000002",
    ]


def test_sequence_registry_marks_duplicate_sequence() -> None:
    assigned = assign_internal_ids(
        [
            FastaRecord(seq_id="seq-a", header="seq-a", sequence="ACGT"),
            FastaRecord(seq_id="seq-b", header="seq-b", sequence="acgu"),
        ],
        dataset_name="dataset-20",
        prefix="D20",
    )
    registry = SequenceRegistry()

    memberships = registry.add_many(assigned)

    assert memberships[0].is_duplicate_sequence is False
    assert memberships[1].is_duplicate_sequence is True
    assert memberships[0].unique_seq_id == memberships[1].unique_seq_id


def test_unique_sequence_registry_has_one_record_for_identical_sequences() -> None:
    assigned = assign_internal_ids(
        [
            FastaRecord(seq_id="seq-a", header="seq-a", sequence="ACGT"),
            FastaRecord(seq_id="seq-b", header="seq-b", sequence="acgu"),
        ],
        dataset_name="dataset-20",
        prefix="D20",
    )
    registry = SequenceRegistry()
    registry.add_many(assigned)

    unique_records = registry.unique_records()

    assert len(unique_records) == 1
    assert unique_records[0].representative_internal_seq_id == "D20_000001"
    assert unique_records[0].first_seen_dataset == "dataset-20"
    assert registry.representative_internal_ids() == ["D20_000001"]
