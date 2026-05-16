from __future__ import annotations

import csv
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
from typer.testing import CliRunner

from autotax2.cli import app
from autotax2.io import read_fasta
from autotax2.silva import classify_unresolved, parse_silva_header


runner = CliRunner()
FIXTURES = Path("tests") / "fixtures"
GTDB_AR53 = FIXTURES / "gtdb_ar53_taxonomy_r232.tsv"
GTDB_BAC120 = FIXTURES / "gtdb_bac120_taxonomy_r232.tsv"


@pytest.fixture
def silva_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_silva_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_init_includes_all_valid_domains_and_converts_rna_to_dna(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "metadata.tsv"
    outdir = silva_tmp_dir / "build"
    named_taxonomy = (
        "Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum"
    )
    bacterial_taxonomy = (
        "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
        "Vibrionaceae;Vibrio;Vibrio halioticoli"
    )
    fasta_path.write_text(
        f">AR1 {named_taxonomy}\n"
        "ACGUACGU\n"
        ">AR2 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium sp.\n"
        "ACGUACGA\n"
        ">AR3 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;unidentified;unidentified\n"
        "ACGUACGC\n"
        f">BA1 {bacterial_taxonomy}\n"
        "ACGUACGG\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata_path)

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "registry" / "taxon_nodes.tsv").exists()
    assert (outdir / "registry" / "sequence_registry.tsv").exists()
    assert (outdir / "registry" / "dataset_registry.tsv").exists()
    assert (outdir / "registry" / "name_index.tsv").exists()
    assert (outdir / "registry" / "cluster_to_taxon.tsv").exists()
    assert (outdir / "registry" / "representative_registry.tsv").exists()
    assert (outdir / "registry" / "placeholder_counters.yaml").exists()
    assert (outdir / "registry" / "placeholder_counters.tsv").exists()
    assert (outdir / "registry" / "protected_taxa_snapshot.tsv").exists()
    assert (outdir / "silva" / "silva_named_backbone.fa").exists()
    assert (outdir / "silva" / "silva_named_backbone.tax.tsv").exists()
    assert (outdir / "silva" / "silva_unresolved.fa").exists()
    assert (outdir / "silva" / "silva_unresolved.tsv").exists()
    assert (outdir / "silva" / "silva_rejected.tsv").exists()
    assert (outdir / "logs").exists()
    log_names = {path.name for path in (outdir / "logs").glob("*.log")}
    assert any(name.startswith("init_start_date") for name in log_names)
    assert any(name.startswith("init_date") for name in log_names)
    assert "start_log=" in result.output
    assert "audit_log=" in result.output

    named_fasta = read_fasta(outdir / "silva" / "silva_named_backbone.fa")
    unresolved_fasta = read_fasta(outdir / "silva" / "silva_unresolved.fa")
    named_rows = _read_tsv(outdir / "silva" / "silva_named_backbone.tax.tsv")
    unresolved_rows = _read_tsv(outdir / "silva" / "silva_unresolved.tsv")
    sequence_rows = _read_tsv(outdir / "registry" / "sequence_registry.tsv")
    taxon_rows = _read_tsv(outdir / "registry" / "taxon_nodes.tsv")
    representative_rows = _read_tsv(outdir / "registry" / "representative_registry.tsv")

    assert [record.seq_id for record in named_fasta] == ["AR1", "BA1"]
    assert sorted(record.seq_id for record in unresolved_fasta) == ["AR2", "AR3"]
    assert "BA1" in {row["seq_id"] for row in sequence_rows}
    assert all("U" not in record.sequence for record in named_fasta + unresolved_fasta)

    assert named_rows[0]["seq_id"] == "AR1"
    assert named_rows[0]["taxonomy_7rank"] == named_taxonomy
    assert named_rows[0]["genus"] == "Methanobacterium"
    assert named_rows[0]["species"] == "Methanobacterium formicicum"
    assert named_rows[1]["seq_id"] == "BA1"
    assert named_rows[1]["taxonomy_7rank"] == bacterial_taxonomy

    unresolved_by_id = {row["seq_id"]: row for row in unresolved_rows}
    assert unresolved_by_id["AR2"]["lowest_reliable_rank"] == "genus"
    assert unresolved_by_id["AR2"]["unresolved_ranks"] == "species"
    assert unresolved_by_id["AR2"]["genus"] == "Methanobacterium"
    assert unresolved_by_id["AR2"]["species"] == "Methanobacterium sp."

    assert unresolved_by_id["AR3"]["lowest_reliable_rank"] == "family"
    assert unresolved_by_id["AR3"]["unresolved_ranks"] == "genus,species"
    assert unresolved_by_id["AR3"]["genus"] == "unidentified"
    assert "g__SILVA" not in (outdir / "silva" / "silva_unresolved.tsv").read_text(
        encoding="utf-8"
    )
    assert "s__SILVA" not in (outdir / "silva" / "silva_unresolved.tsv").read_text(
        encoding="utf-8"
    )

    assert taxon_rows
    assert {row["protected"] for row in taxon_rows} == {"true"}
    assert {row["is_silva_named"] for row in taxon_rows} == {"true"}

    named_sequence = next(row for row in sequence_rows if row["seq_id"] == "AR1")
    unresolved_sequence = next(row for row in sequence_rows if row["seq_id"] == "AR2")
    assert named_sequence["protected"] == "true"
    assert named_sequence["is_silva_named"] == "true"
    assert named_sequence["is_silva_unresolved"] == "false"
    assert named_sequence["original_taxonomy"] == named_taxonomy
    assert named_sequence["taxon_id"]
    assert unresolved_sequence["protected"] == "false"
    assert unresolved_sequence["is_silva_unresolved"] == "true"
    assert representative_rows[0]["representative_seq_id"] == "AR1"
    assert representative_rows[0]["source_category"] == "named_silva"
    assert representative_rows[0]["representative_reason"] == "first_silva_named_for_species"


def test_init_rejects_records_with_empty_domain(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "metadata.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">GOOD Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n"
        ">EMPTY_DOMAIN ;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n"
        ">NO_TAXONOMY\n"
        "ACGT\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata_path)

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    sequence_ids = {row["seq_id"] for row in _read_tsv(outdir / "registry" / "sequence_registry.tsv")}
    rejected = {row["seq_id"]: row["reject_reason"] for row in _read_tsv(outdir / "silva" / "silva_rejected.tsv")}

    assert sequence_ids == {"GOOD"}
    assert rejected == {"EMPTY_DOMAIN": "empty_domain", "NO_TAXONOMY": "empty_domain"}


def test_init_rejects_invalid_silva_sequences_except_type_strains(
    silva_tmp_dir: Path,
) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "type.tsv"
    outdir = silva_tmp_dir / "build"
    taxonomy = (
        "Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum"
    )
    fasta_path.write_text(
        f">TYPE_BAD {taxonomy}\n"
        "ACGTN\n"
        f">NON_TYPE_BAD {taxonomy}\n"
        "ACGTR\n"
        f">NON_TYPE_GOOD {taxonomy}\n"
        "ACGTU\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "primaryAccession\torganism_name\tstrain\ttype_material\n"
        "TYPE_BAD\tMethanobacterium formicicum\tDSM 1535\ttype strain\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    named_fasta = {record.seq_id: record.sequence for record in read_fasta(outdir / "silva" / "silva_named_backbone.fa")}
    rejected = {row["seq_id"]: row["reject_reason"] for row in _read_tsv(outdir / "silva" / "silva_rejected.tsv")}
    sequence_ids = {row["seq_id"] for row in _read_tsv(outdir / "registry" / "sequence_registry.tsv")}

    assert named_fasta["TYPE_BAD"] == "ACGTN"
    assert named_fasta["NON_TYPE_GOOD"] == "ACGTT"
    assert "NON_TYPE_BAD" not in sequence_ids
    assert rejected == {"NON_TYPE_BAD": "invalid_sequence_characters"}


def test_environmental_tail_preserves_candidate_placeholder_framework() -> None:
    examples = [
        ("archaeon", "contains archaeon"),
        ("metagenome", "contains metagenome"),
        ("uncultured", "contains uncultured"),
        ("unidentified", "contains unidentified"),
        ("marine", "environmental descriptor marine"),
    ]
    for species_tail, expected_reason in examples:
        _, taxonomy = parse_silva_header(
            "SEQ Archaea;Nanoarchaeota;Nanoarchaeia;Woesearchaeales;GW2011;AR15;"
            f"{species_tail}"
        )

        classification = classify_unresolved(taxonomy)

        assert classification.lowest_reliable_rank == "genus"
        assert classification.unresolved_ranks == ("species",)
        assert expected_reason in classification.unresolved_reason


def test_incertae_and_standalone_candidatus_start_unresolved_at_their_rank() -> None:
    cases = [
        (
            "SEQ Archaea;Hydrothermarchaeota;Hydrothermarchaeia;Hydrothermarchaeales;Incertae",
            "order",
            ("family", "genus", "species"),
        ),
        (
            "SEQ Archaea;Iainarchaeota;Iainarchaeia;Iainarchaeales;Iainarchaeaceae;Candidatus",
            "family",
            ("genus", "species"),
        ),
        (
            "SEQ Archaea;Incertae",
            "domain",
            ("phylum", "class", "order", "family", "genus", "species"),
        ),
    ]
    for header, lowest_rank, unresolved_ranks in cases:
        _, taxonomy = parse_silva_header(header)

        classification = classify_unresolved(taxonomy)

        assert classification.lowest_reliable_rank == lowest_rank
        assert classification.unresolved_ranks == unresolved_ranks


def test_init_outputs_real_tabs_and_newlines(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "metadata.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">AR1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata_path)

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    content = (outdir / "silva" / "silva_named_backbone.tax.tsv").read_text(
        encoding="utf-8"
    )
    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content


def test_init_attaches_type_strain_metadata(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "type.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">AR1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "primaryAccession\torganism_name\tstrain\ttype_material\n"
        "AR1\tMethanobacterium formicicum\tDSM 1535\ttype strain\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    sequence_rows = _read_tsv(outdir / "registry" / "sequence_registry.tsv")
    representative_rows = _read_tsv(outdir / "registry" / "representative_registry.tsv")

    assert sequence_rows[0]["is_type_strain"] == "true"
    assert sequence_rows[0]["type_species_name"] == "Methanobacterium formicicum"
    assert sequence_rows[0]["type_strain_id"] == "DSM 1535"
    assert sequence_rows[0]["type_source"] == "SILVA full_metadata"
    assert sequence_rows[0]["type_evidence"] == "type_material=type strain"
    assert representative_rows[0]["representative_seq_id"] == "AR1"
    assert representative_rows[0]["is_type_strain"] == "true"
    assert representative_rows[0]["representative_reason"] == "type_strain"


def test_init_uses_type_strain_catalog_and_strips_species_strain_suffix(
    silva_tmp_dir: Path,
) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "type.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">AR1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum DSM 1535\n"
        "ACGT\n"
        ">AR2 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum JCM 101\n"
        "ACGA\n"
        ">AR3 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;uncultured archaeon SAGMA-B\n"
        "ACGC\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "primaryAccession\torganism_name\tstrain\ttype_material\n"
        "AR1\tMethanobacterium formicicum\tDSM 1535\ttype strain\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    named_rows = _read_tsv(outdir / "silva" / "silva_named_backbone.tax.tsv")
    unresolved_rows = _read_tsv(outdir / "silva" / "silva_unresolved.tsv")

    assert [row["seq_id"] for row in named_rows] == ["AR1", "AR2"]
    assert {row["species"] for row in named_rows} == {"Methanobacterium formicicum"}
    assert all("DSM 1535" not in row["taxonomy_7rank"] for row in named_rows)
    assert all("JCM 101" not in row["taxonomy_7rank"] for row in named_rows)

    assert unresolved_rows[0]["seq_id"] == "AR3"
    assert unresolved_rows[0]["lowest_reliable_rank"] == "genus"
    assert unresolved_rows[0]["unresolved_ranks"] == "species"
    assert "SAGMA-B" in unresolved_rows[0]["species"]


def test_init_preserves_gtdb_candidate_placeholder_lineages(
    silva_tmp_dir: Path,
) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "metadata.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">LOKI Archaea;Asgardarchaeota;Lokiarchaeia;Lokiarchaeales;"
        "Lokiarchaeaceae;Lokiarchaeum;Lokiarchaeum sp000002\n"
        "ACGT\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata_path)

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    named_rows = _read_tsv(outdir / "silva" / "silva_named_backbone.tax.tsv")
    unresolved_rows = _read_tsv(outdir / "silva" / "silva_unresolved.tsv")

    assert unresolved_rows == []
    assert named_rows[0]["phylum"] == "Asgardarchaeota"
    assert named_rows[0]["class"] == "Lokiarchaeia"
    assert named_rows[0]["genus"] == "Lokiarchaeum"
    assert named_rows[0]["species"] == "Lokiarchaeum sp000002"

    legal_rows = _read_tsv(outdir / "registry" / "legal_name_catalog.tsv")
    assert {
        (row["rank"], row["name"], row["source"])
        for row in legal_rows
        if row["name"].startswith("Loki")
    } == {
        ("class", "Lokiarchaeia", "GTDB r232 taxonomy"),
        ("order", "Lokiarchaeales", "GTDB r232 taxonomy"),
        ("family", "Lokiarchaeaceae", "GTDB r232 taxonomy"),
        ("genus", "Lokiarchaeum", "GTDB r232 taxonomy"),
        ("species", "Lokiarchaeum sp000002", "GTDB r232 taxonomy"),
    }
    assert "SILVAg000001" not in {row["name"] for row in legal_rows}
    assert "SILVAs000001" not in {row["name"] for row in legal_rows}


def test_init_standardizes_organelle_lineages_as_named_backbone(
    silva_tmp_dir: Path,
) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "metadata.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">CHL Eukaryota;Viridiplantae;Chloroplast\n"
        "ACGT\n"
        ">MIT Eukaryota;Opisthokonta;Metazoa;Mitochondrion;uncultured;metagenome\n"
        "ACGA\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata_path)

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    named_rows = {row["seq_id"]: row for row in _read_tsv(outdir / "silva" / "silva_named_backbone.tax.tsv")}
    unresolved_rows = _read_tsv(outdir / "silva" / "silva_unresolved.tsv")

    assert unresolved_rows == []
    assert named_rows["CHL"]["taxonomy_7rank"] == (
        "Eukaryota;Viridiplantae;Chloroplast;Chloroplast;"
        "Chloroplast;Chloroplast;Chloroplast"
    )
    assert named_rows["MIT"]["taxonomy_7rank"] == (
        "Eukaryota;Opisthokonta;Metazoa;Mitochondria;"
        "Mitochondria;Mitochondria;Mitochondria"
    )


def test_init_reads_gzipped_silva_full_metadata_type_material(silva_tmp_dir: Path) -> None:
    import gzip

    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "full_metadata.tsv.gz"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">AR1.1.100 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n"
        ">AR2.1.100 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGA\n",
        encoding="utf-8",
    )
    with gzip.open(metadata_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "primaryAccession\torganism_name\tstrain\ttype_material\n"
            "AR2\tMethanobacterium formicicum\tDSM 1535\ttype strain\n"
        )

    result = runner.invoke(
        app,
        _init_args(fasta_path, outdir, metadata_path),
    )

    assert result.exit_code == 0, result.output
    representative_rows = _read_tsv(outdir / "registry" / "representative_registry.tsv")
    sequence_rows = {row["seq_id"]: row for row in _read_tsv(outdir / "registry" / "sequence_registry.tsv")}

    assert representative_rows[0]["representative_seq_id"] == "AR2.1.100"
    assert representative_rows[0]["representative_reason"] == "type_strain"
    assert sequence_rows["AR2.1.100"]["is_type_strain"] == "true"
    assert sequence_rows["AR2.1.100"]["type_source"] == "SILVA full_metadata"


def test_init_requires_type_strain_metadata(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">AR1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        ["init", "--silva-fasta", str(fasta_path), "--outdir", str(outdir)],
    )

    assert result.exit_code != 0
    assert "type-strain-metadata" in result.output


def test_init_requires_gtdb_taxonomies(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    metadata_path = silva_tmp_dir / "metadata.tsv"
    outdir = silva_tmp_dir / "build"
    fasta_path.write_text(
        ">AR1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum\n"
        "ACGT\n",
        encoding="utf-8",
    )
    _write_empty_metadata(metadata_path)

    result = runner.invoke(
        app,
        [
            "init",
            "--silva-fasta",
            str(fasta_path),
            "--outdir",
            str(outdir),
            "--type-strain-metadata",
            str(metadata_path),
        ],
    )

    assert result.exit_code != 0
    assert "gtdb-ar53-taxonomy" in result.output

    result = runner.invoke(
        app,
        [
            "init",
            "--silva-fasta",
            str(fasta_path),
            "--outdir",
            str(outdir),
            "--type-strain-metadata",
            str(metadata_path),
            "--gtdb-ar53-taxonomy",
            str(GTDB_AR53),
        ],
    )

    assert result.exit_code != 0
    assert "gtdb-bac120-taxonomy" in result.output


def _write_empty_metadata(path: Path) -> None:
    path.write_text(
        "primaryAccession\torganism_name\tstrain\ttype_material\n",
        encoding="utf-8",
    )


def _init_args(fasta_path: Path, outdir: Path, metadata_path: Path) -> list[str]:
    return [
        "init",
        "--silva-fasta",
        str(fasta_path),
        "--outdir",
        str(outdir),
        "--type-strain-metadata",
        str(metadata_path),
        "--gtdb-ar53-taxonomy",
        str(GTDB_AR53),
        "--gtdb-bac120-taxonomy",
        str(GTDB_BAC120),
    ]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
