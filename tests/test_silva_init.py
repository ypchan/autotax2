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


runner = CliRunner()


@pytest.fixture
def silva_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_silva_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_init_domain_filter_and_silva_split(silva_tmp_dir: Path) -> None:
    fasta_path = silva_tmp_dir / "silva.fa"
    outdir = silva_tmp_dir / "build"
    named_taxonomy = (
        "Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;"
        "Methanobacteriaceae;Methanobacterium;Methanobacterium formicicum"
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
        ">BA1 Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
        "Vibrionaceae;Vibrio;Vibrio halioticoli\n"
        "ACGUACGG\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "init",
            "--silva-fasta",
            str(fasta_path),
            "--outdir",
            str(outdir),
            "--domain",
            "Archaea",
        ],
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

    named_fasta = read_fasta(outdir / "silva" / "silva_named_backbone.fa")
    unresolved_fasta = read_fasta(outdir / "silva" / "silva_unresolved.fa")
    named_rows = _read_tsv(outdir / "silva" / "silva_named_backbone.tax.tsv")
    unresolved_rows = _read_tsv(outdir / "silva" / "silva_unresolved.tsv")
    sequence_rows = _read_tsv(outdir / "registry" / "sequence_registry.tsv")
    taxon_rows = _read_tsv(outdir / "registry" / "taxon_nodes.tsv")
    representative_rows = _read_tsv(outdir / "registry" / "representative_registry.tsv")

    assert [record.seq_id for record in named_fasta] == ["AR1"]
    assert sorted(record.seq_id for record in unresolved_fasta) == ["AR2", "AR3"]
    assert "BA1" not in {row["seq_id"] for row in sequence_rows}

    assert named_rows[0]["seq_id"] == "AR1"
    assert named_rows[0]["taxonomy_7rank"] == named_taxonomy
    assert named_rows[0]["genus"] == "Methanobacterium"
    assert named_rows[0]["species"] == "Methanobacterium formicicum"

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

    result = runner.invoke(
        app,
        ["init", "--silva-fasta", str(fasta_path), "--outdir", str(outdir)],
    )

    assert result.exit_code == 0, result.output
    sequence_ids = {row["seq_id"] for row in _read_tsv(outdir / "registry" / "sequence_registry.tsv")}
    rejected = {row["seq_id"]: row["reject_reason"] for row in _read_tsv(outdir / "silva" / "silva_rejected.tsv")}

    assert sequence_ids == {"GOOD"}
    assert rejected == {"EMPTY_DOMAIN": "empty_domain", "NO_TAXONOMY": "empty_domain"}


def test_init_outputs_real_tabs_and_newlines(silva_tmp_dir: Path) -> None:
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
        "seq_id\tis_type_strain\tspecies_name\tstrain_id\tsource\tevidence\n"
        "AR1\ttrue\tMethanobacterium formicicum\tDSM 1535\tmanual\tfixture\n",
        encoding="utf-8",
    )

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

    assert result.exit_code == 0, result.output
    sequence_rows = _read_tsv(outdir / "registry" / "sequence_registry.tsv")
    representative_rows = _read_tsv(outdir / "registry" / "representative_registry.tsv")

    assert sequence_rows[0]["is_type_strain"] == "true"
    assert sequence_rows[0]["type_species_name"] == "Methanobacterium formicicum"
    assert sequence_rows[0]["type_strain_id"] == "DSM 1535"
    assert sequence_rows[0]["type_source"] == "manual"
    assert sequence_rows[0]["type_evidence"] == "fixture"
    assert representative_rows[0]["representative_seq_id"] == "AR1"
    assert representative_rows[0]["is_type_strain"] == "true"
    assert representative_rows[0]["representative_reason"] == "type_strain"


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

    assert result.exit_code == 0, result.output
    representative_rows = _read_tsv(outdir / "registry" / "representative_registry.tsv")
    sequence_rows = {row["seq_id"]: row for row in _read_tsv(outdir / "registry" / "sequence_registry.tsv")}

    assert representative_rows[0]["representative_seq_id"] == "AR2.1.100"
    assert representative_rows[0]["representative_reason"] == "type_strain"
    assert sequence_rows["AR2.1.100"]["is_type_strain"] == "true"
    assert sequence_rows["AR2.1.100"]["type_source"] == "SILVA full_metadata"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
