from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
from typer.testing import CliRunner

from autotax2.io import FastaRecord, read_fasta, write_fasta
from autotax2.vsearch import parse_uc_records
from autotax2.cli import app


FIXTURES = Path("tests") / "fixtures"
runner = CliRunner()


@pytest.fixture
def e2e_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_e2e_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_tiny_mocked_end_to_end_pipeline(
    e2e_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    build = e2e_tmp_dir / "autotax2_build"
    _mock_external_tools(monkeypatch)

    result = runner.invoke(
        app,
        [
            "init",
            "--silva-fasta",
            str(FIXTURES / "tiny_silva.fa"),
            "--outdir",
            str(build),
            "--type-strain-metadata",
            str(FIXTURES / "type_strains.tsv"),
            "--gtdb-ar53-taxonomy",
            str(FIXTURES / "gtdb_ar53_taxonomy_r232.tsv"),
            "--gtdb-bac120-taxonomy",
            str(FIXTURES / "gtdb_bac120_taxonomy_r232.tsv"),
        ],
    )
    assert result.exit_code == 0, result.output
    named_fasta = build / "silva" / "silva_named_backbone.fa"
    unresolved_fasta = build / "silva" / "silva_unresolved.fa"
    named_taxa_before = _named_taxon_signature(_read_tsv(build / "registry" / "taxon_nodes.tsv"))
    assert [record.seq_id for record in read_fasta(named_fasta)] == ["ARCH_NAMED", "BACT_NAMED"]
    assert {record.seq_id for record in read_fasta(unresolved_fasta)} == {"ARCH_SP", "ARCH_UNID"}
    assert "BACT_NAMED" in (build / "silva" / "silva_named_backbone.tax.tsv").read_text(encoding="utf-8")
    assert all(row["protected"] == "true" for row in _read_tsv(build / "registry" / "taxon_nodes.tsv"))

    result = runner.invoke(app, ["resolve", "--build", str(build)])
    assert result.exit_code == 0, result.output
    unresolved_taxa = _read_tsv(build / "silva" / "silva_unresolved_taxa.tsv")
    unresolved_members = _read_tsv(build / "silva" / "silva_unresolved_members.tsv")
    assert any(row["name"] == "g__SILVAg000001" for row in unresolved_taxa)
    assert any(row["name"] == "s__SILVAs000001" for row in unresolved_taxa)
    species_only = next(row for row in unresolved_members if row["seq_id"] == "ARCH_SP")
    genus_unresolved = next(row for row in unresolved_members if row["seq_id"] == "ARCH_UNID")
    assert species_only["genus_placeholder"] == ""
    assert species_only["species_placeholder"] == "s__SILVAs000001"
    assert genus_unresolved["genus_placeholder"] == "g__SILVAg000001"
    assert all(len(row["resolved_taxonomy"].split(";")) == 7 for row in unresolved_members)
    assert _named_taxon_signature(_read_tsv(build / "registry" / "taxon_nodes.tsv")) == named_taxa_before

    result = runner.invoke(
        app,
        [
            "prepare",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(FIXTURES / "tiny_dataset_d20.fa"),
            "--domain",
            "Archaea",
        ],
    )
    assert result.exit_code == 0, result.output
    dataset_dir = build / "datasets" / "01_digester2020"
    id_rows = _read_tsv(dataset_dir / "sequence_id_map.tsv")
    membership_rows = _read_tsv(dataset_dir / "sequence_membership.tsv")
    unique_rows = _read_tsv(dataset_dir / "unique_sequence_registry.tsv")
    assert [row["internal_seq_id"] for row in id_rows[:2]] == ["D20_000001", "D20_000002"]
    assert id_rows[0]["original_seq_id"] == "dup_silva"
    assert [row["original_seq_id"] for row in id_rows[-2:]] == ["dup_orig", "dup_orig"]
    assert id_rows[-2]["internal_seq_id"] != id_rows[-1]["internal_seq_id"]
    assert all(row["sequence_md5"] for row in membership_rows)
    assert len(unique_rows) == len({row["sequence_md5"] for row in membership_rows})
    assert read_fasta(dataset_dir / "prepared.ssu.fa")

    result = runner.invoke(app, ["orient", "--build", str(build), "--dataset", "digester2020"])
    assert result.exit_code == 0, result.output
    oriented = read_fasta(dataset_dir / "sina.oriented.fa")
    sina_rows = _read_tsv(dataset_dir / "sina.summary.tsv")
    assert [record.seq_id for record in oriented] == [f"D20_{index:06d}" for index in range(1, 6)]
    assert {row["strand"] for row in sina_rows} >= {"plus", "minus", "unknown"}
    assert any(row["sina_status"] == "sina_missing_output" for row in sina_rows)

    result = runner.invoke(app, ["cluster", "--build", str(build), "--dataset", "digester2020"])
    assert result.exit_code == 0, result.output
    assert (dataset_dir / "internal_clusters" / "species_0.972.uc").exists()
    assert (dataset_dir / "internal_clusters" / "genus_0.901.uc").exists()
    hits = _read_tsv(dataset_dir / "vs_registry.filtered.tsv")
    assert len([row for row in hits if row["query"] == "D20_000004"]) == 2

    result = runner.invoke(app, ["place", "--build", str(build), "--dataset", "digester2020"])
    assert result.exit_code == 0, result.output
    assignments = {row["internal_seq_id"]: row for row in _read_tsv(dataset_dir / "assignments.tsv")}
    created_taxa = _read_tsv(dataset_dir / "created_taxa.tsv")
    created_names = [row["name"] for row in created_taxa]
    assert assignments["D20_000001"]["final_status"] == "duplicate"
    assert assignments["D20_000002"]["final_status"] == "known_like"
    assert assignments["D20_000003"]["final_status"] == "new_species"
    assert assignments["D20_000004"]["final_status"] == "new_genus"
    assert "s__D20s000001" in created_names
    assert "g__D20g000001" in created_names
    assert len(created_names) == len(set(created_names))
    representative_rows = _read_tsv(build / "registry" / "representative_registry.tsv")
    named_rep = next(row for row in representative_rows if row.get("representative_seq_id") == "ARCH_NAMED")
    assert named_rep["source_category"] == "named_silva"
    assert named_rep["representative_reason"] == "type_strain"

    result = runner.invoke(app, ["export", "all", "--build", str(build)])
    assert result.exit_code == 0, result.output
    sintax_path = build / "export" / "sintax" / "autotax2.sintax.fa.gz"
    qiime_fasta = build / "export" / "qiime2" / "reference_sequences.fasta.gz"
    qiime_tax = build / "export" / "qiime2" / "reference_taxonomy.tsv"
    dada_genus = build / "export" / "dada2" / "autotax2_toGenus_trainset.fa.gz"
    dada_species = build / "export" / "dada2" / "autotax2_assignSpecies.fa.gz"
    export_validation = build / "export" / "export_validation.tsv"
    for path in (sintax_path, qiime_fasta, qiime_tax, dada_genus, dada_species):
        assert path.exists()
    assert export_validation.exists()
    sintax_headers = _gzip_headers(sintax_path)
    assert all(";tax=d:" in header and ",g:" in header and ",s:" in header for header in sintax_headers)
    assert all("g__" not in header and "s__" not in header for header in sintax_headers)
    qiime_content = qiime_tax.read_text(encoding="utf-8")
    assert "d__" in qiime_content and "g__" in qiime_content
    dada_species_headers = _gzip_headers(dada_species)
    assert all(len(header.split()) >= 3 for header in dada_species_headers)
    assert all(";" not in header and "g__" not in header and "s__" not in header for header in dada_species_headers)

    result = runner.invoke(app, ["summarize", "--build", str(build)])
    assert result.exit_code == 0, result.output
    report_dir = build / "reports"
    assert (report_dir / "global_summary.tsv").exists()
    assert (report_dir / "dataset_delta_summary.tsv").exists()
    assert (report_dir / "dataset_overlap_matrix.tsv").exists()
    delta = _read_tsv(report_dir / "dataset_delta_summary.tsv")[0]
    global_summary = _read_tsv(report_dir / "global_summary.tsv")[0]
    assert delta["assigned_named_silva"] == "2"
    assert delta["assigned_unresolved_silva"] == "1"
    assert int(global_summary["duplicate_sequences"]) >= 1

    result = runner.invoke(app, ["validate", "--build", str(build)])
    assert result.exit_code == 0, result.output
    assert (report_dir / "validation_report.md").exists()
    assert (report_dir / "validation_report.tsv").exists()


def test_add_command_orchestrates_dataset_workflow(
    e2e_tmp_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    build = e2e_tmp_dir / "autotax2_build"
    _mock_external_tools(monkeypatch)

    result = runner.invoke(
        app,
        [
            "init",
            "--silva-fasta",
            str(FIXTURES / "tiny_silva.fa"),
            "--outdir",
            str(build),
            "--type-strain-metadata",
            str(FIXTURES / "type_strains.tsv"),
            "--gtdb-ar53-taxonomy",
            str(FIXTURES / "gtdb_ar53_taxonomy_r232.tsv"),
            "--gtdb-bac120-taxonomy",
            str(FIXTURES / "gtdb_bac120_taxonomy_r232.tsv"),
        ],
    )
    assert result.exit_code == 0, result.output

    result = runner.invoke(app, ["resolve", "--build", str(build)])
    assert result.exit_code == 0, result.output

    result = runner.invoke(
        app,
        [
            "add",
            "--build",
            str(build),
            "--name",
            "digester2020",
            "--prefix",
            "D20",
            "--fasta",
            str(FIXTURES / "tiny_dataset_d20.fa"),
            "--domain",
            "Archaea",
            "--threads",
            "4",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Completed dataset add workflow" in result.output

    dataset_dir = build / "datasets" / "01_digester2020"
    assert (dataset_dir / "sina.oriented.fa").exists()
    assert (dataset_dir / "internal_clusters" / "species_0.972.uc").exists()
    assert (dataset_dir / "assignments.tsv").exists()
    assert (dataset_dir / "placement_evidence.tsv").exists()
    assert (dataset_dir / "sina_candidate_diagnostics.tsv").exists()
    assert (build / "reports" / "global_summary.tsv").exists()
    assert (build / "reports" / "validation_report.tsv").exists()
    assert any(path.name.startswith("add_date") for path in (build / "logs").glob("*.log"))


def test_all_cli_help_commands_work() -> None:
    commands = [
        [],
        ["init"],
        ["resolve"],
        ["prepare"],
        ["orient"],
        ["cluster"],
        ["place"],
        ["export"],
        ["summarize"],
        ["validate"],
    ]
    for command in commands:
        result = runner.invoke(app, [*command, "--help"])
        assert result.exit_code == 0, result.output


def _mock_external_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_subprocess_run(command, *args, **kwargs):
        command = [str(part) for part in command]
        executable = Path(command[0]).name.lower()
        if command[-1] == "--version":
            if "sina" in executable:
                return subprocess.CompletedProcess(command, 0, stdout="SINA 1.7.2\n", stderr="")
            if "vsearch" in executable:
                return subprocess.CompletedProcess(command, 0, stdout="vsearch v2.29.3\n", stderr="")
        if "sina" in executable:
            output_path = Path(command[command.index("-o") + 1])
            shutil.copyfile(FIXTURES / "fake_sina_output.fa", output_path)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="sina mock\n")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    def fake_silva_cluster(fasta_path, uc_path, identity, *args, **kwargs):
        del fasta_path, identity, args, kwargs
        Path(uc_path).write_text(
            "S\t0\t1000\t*\t*\t*\t*\t*\tARCH_SP\t*\n"
            "S\t1\t1000\t*\t*\t*\t*\t*\tARCH_UNID\t*\n",
            encoding="utf-8",
        )
        return ["vsearch", "--mock-silva"]

    def fake_dataset_cluster(fasta_path, uc_path, identity, *args, **kwargs):
        del args
        uc_path = Path(uc_path)
        if "species_" in uc_path.name:
            shutil.copyfile(FIXTURES / "fake_species_0.972.uc", uc_path)
        elif "genus_" in uc_path.name:
            shutil.copyfile(FIXTURES / "fake_genus_0.901.uc", uc_path)
        else:
            ids = [record.seq_id for record in read_fasta(fasta_path)]
            uc_path.write_text(
                "".join(f"S\t{index}\t1000\t*\t*\t*\t*\t*\t{seq_id}\t*\n" for index, seq_id in enumerate(ids)),
                encoding="utf-8",
            )
        centroids_path = Path(kwargs.get("centroids_path") or uc_path.with_suffix(".centroids.fa"))
        records_by_id = {record.seq_id: record for record in read_fasta(fasta_path)}
        centroid_records = [
            FastaRecord(seq_id=record.query_label, header=record.query_label, sequence=records_by_id[record.query_label].sequence)
            for record in parse_uc_records(uc_path)
            if record.record_type == "S" and record.query_label in records_by_id
        ]
        write_fasta(centroid_records, centroids_path)
        return ["vsearch", "--cluster_fast", str(fasta_path), "--id", f"{identity:.3f}"]

    def fake_registry_search(command):
        userout_path = Path(command[command.index("--userout") + 1])
        shutil.copyfile(FIXTURES / "fake_vs_registry.tsv", userout_path)

    monkeypatch.setattr("subprocess.run", fake_subprocess_run)
    monkeypatch.setattr("autotax2.silva.run_vsearch_cluster", fake_silva_cluster)
    monkeypatch.setattr("autotax2.vsearch.run_vsearch_cluster", fake_dataset_cluster)
    monkeypatch.setattr("autotax2.vsearch.run_vsearch_search", fake_registry_search)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))


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


def _gzip_headers(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return [line.strip()[1:] for line in handle if line.startswith(">")]
