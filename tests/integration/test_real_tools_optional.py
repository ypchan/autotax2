from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


pytestmark = pytest.mark.integration

if os.environ.get("AUTOTAX2_RUN_INTEGRATION") != "1":
    pytest.skip(
        "Set AUTOTAX2_RUN_INTEGRATION=1 to run optional real-tool integration tests.",
        allow_module_level=True,
    )


def test_real_tools_optional_integration(tmp_path: Path) -> None:
    """Run the optional shell integration workflow against user-provided files."""
    silva_fasta = os.environ.get("AUTOTAX2_INTEGRATION_SILVA_FASTA")
    type_strain_metadata = os.environ.get("AUTOTAX2_INTEGRATION_TYPE_STRAIN_METADATA")
    gtdb_ar53_taxonomy = os.environ.get("AUTOTAX2_INTEGRATION_GTDB_AR53_TAXONOMY")
    gtdb_bac120_taxonomy = os.environ.get("AUTOTAX2_INTEGRATION_GTDB_BAC120_TAXONOMY")
    dataset_fasta = os.environ.get("AUTOTAX2_INTEGRATION_DATASET_FASTA")
    if not (
        silva_fasta
        and type_strain_metadata
        and gtdb_ar53_taxonomy
        and gtdb_bac120_taxonomy
        and dataset_fasta
    ):
        pytest.skip(
            "Set AUTOTAX2_INTEGRATION_SILVA_FASTA, "
            "AUTOTAX2_INTEGRATION_TYPE_STRAIN_METADATA, "
            "AUTOTAX2_INTEGRATION_GTDB_AR53_TAXONOMY, "
            "AUTOTAX2_INTEGRATION_GTDB_BAC120_TAXONOMY, and "
            "AUTOTAX2_INTEGRATION_DATASET_FASTA to run this test."
        )

    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash is required to run scripts/run_real_integration_test.sh")

    outdir = Path(os.environ.get("AUTOTAX2_INTEGRATION_OUTDIR", str(tmp_path / "autotax2_real_test")))
    command = [
        bash,
        "scripts/run_real_integration_test.sh",
        "--silva-fasta",
        silva_fasta,
        "--type-strain-metadata",
        type_strain_metadata,
        "--gtdb-ar53-taxonomy",
        gtdb_ar53_taxonomy,
        "--gtdb-bac120-taxonomy",
        gtdb_bac120_taxonomy,
        "--dataset-fasta",
        dataset_fasta,
        "--outdir",
        str(outdir),
        "--domain",
        os.environ.get("AUTOTAX2_INTEGRATION_DOMAIN", "Archaea"),
        "--dataset-name",
        os.environ.get("AUTOTAX2_INTEGRATION_DATASET_NAME", "testdataset"),
        "--prefix",
        os.environ.get("AUTOTAX2_INTEGRATION_PREFIX", "TST"),
        "--threads",
        os.environ.get("AUTOTAX2_INTEGRATION_THREADS", "4"),
    ]
    if os.environ.get("AUTOTAX2_INTEGRATION_SINA_BIN"):
        command.extend(["--sina-bin", os.environ["AUTOTAX2_INTEGRATION_SINA_BIN"]])
    if os.environ.get("AUTOTAX2_INTEGRATION_VSEARCH_BIN"):
        command.extend(["--vsearch-bin", os.environ["AUTOTAX2_INTEGRATION_VSEARCH_BIN"]])
    if os.environ.get("AUTOTAX2_INTEGRATION_SINA_REFERENCE"):
        command.extend(["--sina-reference", os.environ["AUTOTAX2_INTEGRATION_SINA_REFERENCE"]])
    if os.environ.get("AUTOTAX2_INTEGRATION_SINA_SEARCH_DB"):
        command.extend(["--sina-search-db", os.environ["AUTOTAX2_INTEGRATION_SINA_SEARCH_DB"]])
    if os.environ.get("AUTOTAX2_INTEGRATION_SEARCH_CANDIDATES", "1") == "0":
        command.append("--no-search-candidates")
    if os.environ.get("AUTOTAX2_INTEGRATION_REQUIRE_SINA_CANDIDATES") == "1":
        command.append("--require-sina-candidates")
    if os.environ.get("AUTOTAX2_INTEGRATION_STRICT_VALIDATE") == "1":
        command.append("--strict-validate")
    completed = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert completed.returncode == 0, (
        "Optional real-tool integration failed.\n"
        f"STDOUT:\n{completed.stdout}\n"
        f"STDERR:\n{completed.stderr}\n"
    )

    expected_outputs = [
        outdir / "export" / "sintax" / "autotax2.sintax.fa.gz",
        outdir / "export" / "dada2" / "autotax2_toGenus_trainset.fa.gz",
        outdir / "export" / "dada2" / "autotax2_assignSpecies.fa.gz",
        outdir / "export" / "qiime2" / "reference_sequences.fasta.gz",
        outdir / "export" / "qiime2" / "reference_taxonomy.tsv",
        outdir / "reports" / "global_summary.tsv",
        outdir / "reports" / "dataset_delta_summary.tsv",
        outdir / "reports" / "validation_report.md",
        outdir / "export" / "export_validation.tsv",
        outdir / "silva" / "silva_unresolved_evidence.tsv",
    ]
    for path in expected_outputs:
        assert path.exists() and path.stat().st_size > 0, path

    dataset_dirs = list((outdir / "datasets").glob(f"*_{os.environ.get('AUTOTAX2_INTEGRATION_DATASET_NAME', 'testdataset')}"))
    assert dataset_dirs
    dataset_dir = dataset_dirs[0]
    assert (dataset_dir / "placement_evidence.tsv").exists()
    if os.environ.get("AUTOTAX2_INTEGRATION_SEARCH_CANDIDATES", "1") != "0":
        assert (dataset_dir / "sina_candidate_diagnostics.tsv").exists()
