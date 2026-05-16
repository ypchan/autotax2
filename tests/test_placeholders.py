from __future__ import annotations

import pytest

from autotax2.placeholders import (
    PlaceholderAllocator,
    PlaceholderRank,
    make_placeholder_id,
    parse_placeholder_id,
)


@pytest.mark.parametrize(
    ("rank", "prefix", "expected"),
    [
        (PlaceholderRank.PHYLUM, "SILVA", "p__SILVAp000001"),
        (PlaceholderRank.CLASS, "SILVA", "c__SILVAc000001"),
        (PlaceholderRank.ORDER, "SILVA", "o__SILVAo000001"),
        (PlaceholderRank.FAMILY, "SILVA", "f__SILVAf000001"),
        (PlaceholderRank.GENUS, "SILVA", "g__SILVAg000001"),
        (PlaceholderRank.SPECIES, "SILVA", "s__SILVAs000001"),
        (PlaceholderRank.GENUS, "D20", "g__D20g000001"),
        (PlaceholderRank.SPECIES, "D20", "s__D20s000001"),
    ],
)
def test_make_placeholder_id(rank: PlaceholderRank, prefix: str, expected: str) -> None:
    assert make_placeholder_id(rank, prefix, 1) == expected


def test_parse_placeholder_id() -> None:
    parsed = parse_placeholder_id("g__D20g000042")

    assert parsed.rank is PlaceholderRank.GENUS
    assert parsed.dataset_prefix == "D20"
    assert parsed.ordinal == 42


def test_placeholder_ordinal_must_be_positive() -> None:
    with pytest.raises(ValueError):
        make_placeholder_id(PlaceholderRank.GENUS, "D20", 0)


def test_allocator_allocates_by_rank_and_prefix() -> None:
    allocator = PlaceholderAllocator()

    assert allocator.allocate(PlaceholderRank.GENUS, "D20") == "g__D20g000001"
    assert allocator.allocate(PlaceholderRank.GENUS, "D20") == "g__D20g000002"
    assert allocator.allocate(PlaceholderRank.SPECIES, "D20") == "s__D20s000001"


def test_allocator_never_reuses_deprecated_ids() -> None:
    allocator = PlaceholderAllocator()

    first = allocator.allocate(PlaceholderRank.GENUS, "D20")
    allocator.deprecate(first)

    assert allocator.is_deprecated("g__D20g000001")
    assert allocator.allocate(PlaceholderRank.GENUS, "D20") == "g__D20g000002"


def test_allocator_advances_past_existing_ids() -> None:
    allocator = PlaceholderAllocator(
        issued_ids={"g__D20g000003"},
        deprecated_ids={"g__D20g000004"},
    )

    assert allocator.allocate(PlaceholderRank.GENUS, "D20") == "g__D20g000005"
