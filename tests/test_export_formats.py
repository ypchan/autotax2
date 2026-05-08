from __future__ import annotations

import csv
import gzip
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest

from autotax2.export import (
    REQUIRED_EXPORT_FORMATS,
    export_references,
    format_taxonomy_dada2_genus,
    format_taxonomy_dada2_species,
    format_taxonomy_qiime2,
    format_taxonomy_sintax,
)
from autotax2.export.dada2 import DADA2_SPECIES_FORMAT, DADA2_TOGENUS_FORMAT
from autotax2.export.qiime2 import QIIME2_REFERENCE_SEQUENCES_FORMAT, QIIME2_TAXONOMY_FORMAT
from autotax2.export.sintax import SINTAX_FORMAT, format_sintax_header
from autotax2.registry import sequence_md5


@pytest.fixture
def export_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_export_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_required_export_formats_are_declared() -> None:
    assert REQUIRED_EXPORT_FORMATS == (
        SINTAX_FORMAT,
        QIIME2_REFERENCE_SEQUENCES_FORMAT,
        QIIME2_TAXONOMY_FORMAT,
        DADA2_TOGENUS_FORMAT,
        DADA2_SPECIES_FORMAT,
    )


def test_sintax_header_uses_rank_codes_without_double_underscore() -> None:
    taxonomy = _placeholder_taxonomy()
    header = format_sintax_header("REF_PLACEHOLDER", taxonomy)

    assert header == (
        "REF_PLACEHOLDER;tax=d:Archaea,p:Thermoproteota,c:Nitrososphaeria,"
        "o:Nitrososphaerales,f:Nitrososphaeraceae,g:SILVAg000001,s:SILVAs000001;"
    )
    assert "d__" not in header
    assert "g__" not in header


def test_qiime2_taxonomy_retains_rank_prefixes() -> None:
    formatted = format_taxonomy_qiime2(_placeholder_taxonomy())

    assert formatted.startswith("d__Archaea; p__Thermoproteota")
    assert "g__SILVAg000001" in formatted


def test_dada2_to_genus_strips_prefixes_and_omits_species() -> None:
    formatted = format_taxonomy_dada2_genus(_placeholder_taxonomy())

    assert formatted.endswith("Nitrososphaeraceae;SILVAg000001")
    assert "SILVAs000001" not in formatted
    assert "g__" not in formatted


def test_dada2_assignspecies_formats_binomial_and_placeholder_species() -> None:
    assert format_taxonomy_dada2_species(_placeholder_taxonomy()) == "SILVAg000001 SILVAs000001"
    assert format_taxonomy_dada2_species(_vibrio_taxonomy()) == "Vibrio halioticoli"


def test_export_sintax_headers_and_gzip_are_valid(export_tmp_dir: Path) -> None:
    build = _make_export_build(export_tmp_dir)

    export_references(build, "sintax")

    content = _read_gzip_text(build / "export" / "sintax" / "autotax2.sintax.fa.gz")
    headers = [line for line in content.splitlines() if line.startswith(">")]

    assert headers[0].startswith(">REF_PLACEHOLDER;tax=d:Archaea,p:Thermoproteota")
    assert "g:SILVAg000001" in headers[0]
    assert "g__" not in headers[0]
    assert "s__" not in headers[0]
    assert len(headers) == 2


def test_export_qiime2_outputs_taxonomy_header_and_prefixes(export_tmp_dir: Path) -> None:
    build = _make_export_build(export_tmp_dir)

    export_references(build, "qiime2")

    taxonomy_path = build / "export" / "qiime2" / "reference_taxonomy.tsv"
    content = taxonomy_path.read_text(encoding="utf-8")
    rows = _read_tsv(taxonomy_path)

    assert content.splitlines()[0] == "Feature ID\tTaxon"
    assert rows[0]["Taxon"].startswith("d__Archaea; p__Thermoproteota")
    assert "g__SILVAg000001" in rows[0]["Taxon"]
    assert (build / "export" / "qiime2" / "reference_sequences.fasta.gz").exists()


def test_export_dada2_togenus_and_assignspecies_headers(export_tmp_dir: Path) -> None:
    build = _make_export_build(export_tmp_dir)

    export_references(build, "dada2")

    genus_content = _read_gzip_text(build / "export" / "dada2" / "autotax2_toGenus_trainset.fa.gz")
    species_content = _read_gzip_text(build / "export" / "dada2" / "autotax2_assignSpecies.fa.gz")
    genus_header = genus_content.splitlines()[0]
    species_headers = [line for line in species_content.splitlines() if line.startswith(">")]

    assert genus_header == (
        ">REF_PLACEHOLDER "
        "Archaea;Thermoproteota;Nitrososphaeria;Nitrososphaerales;"
        "Nitrososphaeraceae;SILVAg000001"
    )
    assert "SILVAs000001" not in genus_header
    assert "g__" not in genus_header
    assert ">REF_PLACEHOLDER SILVAg000001 SILVAs000001" in species_headers
    assert ">REF_VIBRIO Vibrio halioticoli" in species_headers
    assert all(";" not in header for header in species_headers)
    assert all("g__" not in header and "s__" not in header for header in species_headers)


def test_duplicate_md5_and_deprecated_taxa_are_not_exported(export_tmp_dir: Path) -> None:
    build = _make_export_build(export_tmp_dir)

    export_references(build, "all")

    sintax_content = _read_gzip_text(build / "export" / "sintax" / "autotax2.sintax.fa.gz")
    headers = [line for line in sintax_content.splitlines() if line.startswith(">")]

    assert len(headers) == 2
    assert not any("REF_DUP" in header for header in headers)
    assert not any("REF_DEPRECATED" in header for header in headers)


def test_missing_7rank_taxonomy_raises_validation_error(export_tmp_dir: Path) -> None:
    build = _make_export_build(export_tmp_dir, incomplete=True)

    with pytest.raises(ValueError, match="7-rank"):
        export_references(build, "sintax")


def test_export_manifest_has_real_tabs_and_newlines(export_tmp_dir: Path) -> None:
    build = _make_export_build(export_tmp_dir)

    export_references(build, "all")

    manifest = build / "export" / "export_manifest.tsv"
    content = manifest.read_text(encoding="utf-8")
    rows = _read_tsv(manifest)

    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content
    assert {row["format"] for row in rows} == {
        "sintax",
        "qiime2-reference-sequences",
        "qiime2-taxonomy",
        "dada2-togenus",
        "dada2-assignspecies",
    }
    assert all(row["records_exported"] == "2" for row in rows)


def _make_export_build(tmp_dir: Path, incomplete: bool = False) -> Path:
    build = tmp_dir / "build"
    registry_dir = build / "registry"
    registry_dir.mkdir(parents=True)
    taxon_rows = [
        _taxon("D1", "domain", "d__Archaea", ""),
        _taxon("P1", "phylum", "p__Thermoproteota", "D1"),
        _taxon("C1", "class", "c__Nitrososphaeria", "P1"),
        _taxon("O1", "order", "o__Nitrososphaerales", "C1"),
        _taxon("F1", "family", "f__Nitrososphaeraceae", "O1"),
        _taxon("G1", "genus", "g__SILVAg000001", "F1"),
        _taxon("S1", "species", "s__SILVAs000001", "G1"),
        _taxon("D2", "domain", "d__Bacteria", ""),
        _taxon("P2", "phylum", "p__Pseudomonadota", "D2"),
        _taxon("C2", "class", "c__Gammaproteobacteria", "P2"),
        _taxon("O2", "order", "o__Enterobacterales", "C2"),
        _taxon("F2", "family", "f__Vibrionaceae", "O2"),
        _taxon("G2", "genus", "g__Vibrio", "F2"),
        _taxon("S2", "species", "s__Vibrio halioticoli", "G2"),
        _taxon("SDEP", "species", "s__Deprecated species", "G1", status="deprecated"),
    ]
    if incomplete:
        taxon_rows.append(_taxon("BAD", "species", "s__Missing parent", ""))
    _write_tsv(
        registry_dir / "taxon_nodes.tsv",
        taxon_rows,
        ["taxon_id", "rank", "name", "parent_taxon_id", "status", "source"],
    )
    reps = [
        {"representative_seq_id": "REF_PLACEHOLDER", "taxon_id": "S1", "status": "active"},
        {"representative_seq_id": "REF_VIBRIO", "taxon_id": "S2", "status": "active"},
        {"representative_seq_id": "REF_DUP", "taxon_id": "S1", "status": "active"},
        {"representative_seq_id": "REF_DEPRECATED", "taxon_id": "SDEP", "status": "active"},
    ]
    if incomplete:
        reps = [{"representative_seq_id": "REF_BAD", "taxon_id": "BAD", "status": "active"}]
    _write_tsv(
        registry_dir / "representative_registry.tsv",
        reps,
        ["representative_seq_id", "taxon_id", "status"],
    )
    sequences = {
        "REF_PLACEHOLDER": "ACGTACGTACGT",
        "REF_VIBRIO": "TGCATGCATGCA",
        "REF_DUP": "ACGTACGTACGT",
        "REF_DEPRECATED": "CCCCCCCCCCCC",
        "REF_BAD": "GGGGGGGGGGGG",
    }
    (registry_dir / "current_representatives.fa").write_text(
        "".join(f">{seq_id}\n{sequence}\n" for seq_id, sequence in sequences.items()),
        encoding="utf-8",
    )
    _write_tsv(
        registry_dir / "sequence_registry.tsv",
        [
            {
                "seq_id": seq_id,
                "sequence_md5": sequence_md5(sequence),
                "taxon_id": next((row["taxon_id"] for row in reps if row["representative_seq_id"] == seq_id), ""),
            }
            for seq_id, sequence in sequences.items()
        ],
        ["seq_id", "sequence_md5", "taxon_id"],
    )
    return build


def _placeholder_taxonomy() -> tuple[str, ...]:
    return (
        "d__Archaea",
        "p__Thermoproteota",
        "c__Nitrososphaeria",
        "o__Nitrososphaerales",
        "f__Nitrososphaeraceae",
        "g__SILVAg000001",
        "s__SILVAs000001",
    )


def _vibrio_taxonomy() -> tuple[str, ...]:
    return (
        "d__Bacteria",
        "p__Pseudomonadota",
        "c__Gammaproteobacteria",
        "o__Enterobacterales",
        "f__Vibrionaceae",
        "g__Vibrio",
        "s__Vibrio halioticoli",
    )


def _taxon(
    taxon_id: str,
    rank: str,
    name: str,
    parent: str,
    status: str = "active",
) -> dict[str, str]:
    return {
        "taxon_id": taxon_id,
        "rank": rank,
        "name": name,
        "parent_taxon_id": parent,
        "status": status,
        "source": "SILVA138.2_NR99",
    }


def _read_gzip_text(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return handle.read()


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
