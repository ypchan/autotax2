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
    dataset_fasta = os.environ.get("AUTOTAX2_INTEGRATION_DATASET_FASTA")
    if not silva_fasta or not dataset_fasta:
        pytest.skip(
            "Set AUTOTAX2_INTEGRATION_SILVA_FASTA and "
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
    if os.environ.get("AUTOTAX2_INTEGRATION_STRICT_TOOL_VERSION") == "1":
        command.append("--strict-tool-version")

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
    ]
    for path in expected_outputs:
        assert path.exists() and path.stat().st_size > 0, path
