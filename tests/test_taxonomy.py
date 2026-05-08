from __future__ import annotations

import pytest

from autotax2.taxonomy import (
    Rank,
    Taxon,
    parse_taxon_label,
    parse_taxonomy_path,
    rank_prefix,
    split_taxonomy_path,
)


def test_taxon_prefixed_name() -> None:
    assert Taxon(Rank.GENUS, "Example").prefixed_name == "g__Example"


def test_split_taxonomy_path_ignores_empty_parts() -> None:
    assert split_taxonomy_path("d__Bacteria; p__Firmicutes; ; g__Example") == [
        "d__Bacteria",
        "p__Firmicutes",
        "g__Example",
    ]


def test_rank_prefix_requires_prefixed_label() -> None:
    with pytest.raises(ValueError):
        rank_prefix("Bacteria")


def test_parse_taxon_label() -> None:
    taxon = parse_taxon_label(" g__Example ")

    assert taxon.rank is Rank.GENUS
    assert taxon.name == "Example"
    assert taxon.prefixed_name == "g__Example"


def test_parse_taxonomy_path_is_rank_aware() -> None:
    path = parse_taxonomy_path(
        "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Exampleaceae;"
        "g__Example;s__Example species;"
    )

    assert path.labels == (
        "d__Bacteria",
        "p__Firmicutes",
        "c__Bacilli",
        "o__Lactobacillales",
        "f__Exampleaceae",
        "g__Example",
        "s__Example species",
    )
    assert path.require(Rank.GENUS).name == "Example"
    assert path.to_semicolon_string(trailing_semicolon=True).endswith(";")


def test_parse_taxonomy_path_rejects_duplicate_ranks() -> None:
    with pytest.raises(ValueError):
        parse_taxonomy_path("d__Bacteria;g__One;g__Two")


def test_parse_taxonomy_path_rejects_out_of_order_ranks_by_default() -> None:
    with pytest.raises(ValueError):
        parse_taxonomy_path("g__Example;f__Exampleaceae")
