from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .logging import print_command, step
from .threads import validate_threads
from .utils import ensure_dir


class VsearchError(RuntimeError):
    """Raised when a VSEARCH command fails."""


def _run(command: Sequence[str | Path]) -> None:
    """Run an external command and raise a readable error on failure."""

    print_command(command)

    completed = subprocess.run(
        [str(part) for part in command],
        check=False,
        text=True,
    )

    if completed.returncode != 0:
        raise VsearchError(
            "VSEARCH command failed with exit code "
            f"{completed.returncode}: "
            + " ".join(str(part) for part in command)
        )


def _identity_label(identity: float) -> str:
    """Convert identity threshold to a compact output label.

    Examples
    --------
    0.99 -> otu099
    0.97 -> otu097
    0.9  -> otu090
    """

    if not 0 < identity <= 1:
        raise ValueError("Identity must be > 0 and <= 1.")

    return f"otu{int(round(identity * 100)):03d}"


def _require_method(method: str) -> str:
    allowed = {"cluster_size", "cluster_fast", "cluster_smallmem"}

    if method not in allowed:
        raise ValueError(
            "Invalid VSEARCH clustering method. "
            f"Expected one of {sorted(allowed)}, got {method!r}."
        )

    return method


def _require_strand(strand: str) -> str:
    allowed = {"both", "plus"}

    if strand not in allowed:
        raise ValueError(
            f"Invalid strand value. Expected 'both' or 'plus', got {strand!r}."
        )

    return strand


def dereplicate(
    input_fasta: str | Path,
    outdir: str | Path,
    vsearch: str | Path,
    threads: int,
) -> str:
    """Run VSEARCH full-length dereplication.

    Outputs
    -------
    derep.fasta
        Dereplicated FASTA with size annotations.
    derep.uc
        UC mapping file from original sequences to dereplicated representatives.
    """

    threads = validate_threads(threads)
    outdir = ensure_dir(outdir)

    output_fasta = outdir / "derep.fasta"
    output_uc = outdir / "derep.uc"

    step("VSEARCH dereplication")

    command = [
        vsearch,
        "--derep_fulllength",
        input_fasta,
        "--output",
        output_fasta,
        "--uc",
        output_uc,
        "--sizeout",
        "--threads",
        str(threads),
    ]

    _run(command)
    return str(output_fasta)


def sort_by_size(
    input_fasta: str | Path,
    outdir: str | Path,
    vsearch: str | Path,
    minsize: Optional[int] = None,
) -> str:
    """Sort FASTA records by abundance using VSEARCH sortbysize."""

    outdir = ensure_dir(outdir)

    output_fasta = outdir / "derep_sorted.fasta"

    step("VSEARCH sortbysize")

    command: List[str | Path] = [
        vsearch,
        "--sortbysize",
        input_fasta,
        "--output",
        output_fasta,
        "--sizeout",
    ]

    if minsize is not None:
        if minsize < 1:
            raise ValueError("minsize must be at least 1.")
        command.extend(["--minsize", str(minsize)])

    _run(command)
    return str(output_fasta)


def cluster(
    input_fasta: str | Path,
    outdir: str | Path,
    identity: float,
    method: str,
    vsearch: str | Path,
    threads: int,
    relabel: Optional[str] = None,
) -> Dict[str, str]:
    """Cluster sequences with VSEARCH at one identity threshold."""

    threads = validate_threads(threads)
    method = _require_method(method)

    outdir = ensure_dir(outdir)
    label = _identity_label(identity)

    centroids = outdir / f"{label}_centroids.fasta"
    uc = outdir / f"{label}.uc"

    step(f"VSEARCH {method} at identity {identity}")

    command: List[str | Path] = [
        vsearch,
        f"--{method}",
        input_fasta,
        "--id",
        str(identity),
        "--centroids",
        centroids,
        "--uc",
        uc,
        "--sizein",
        "--sizeout",
        "--threads",
        str(threads),
    ]

    if relabel:
        command.extend(["--relabel", relabel])

    _run(command)

    return {
        "centroids": str(centroids),
        "uc": str(uc),
    }


def sintax(
    input_fasta: str | Path,
    db: str | Path,
    outdir: str | Path,
    cutoff: float,
    vsearch: str | Path,
    threads: int,
) -> str:
    """Classify sequences with VSEARCH SINTAX."""

    threads = validate_threads(threads)

    if not 0 <= cutoff <= 1:
        raise ValueError("SINTAX cutoff must be between 0 and 1.")

    outdir = ensure_dir(outdir)

    output_tsv = outdir / "sintax.tsv"

    step("VSEARCH SINTAX classification")

    command = [
        vsearch,
        "--sintax",
        input_fasta,
        "--db",
        db,
        "--tabbedout",
        output_tsv,
        "--sintax_cutoff",
        str(cutoff),
        "--threads",
        str(threads),
    ]

    _run(command)
    return str(output_tsv)


def usearch_global(
    input_fasta: str | Path,
    db: str | Path,
    outdir: str | Path,
    output_prefix: str,
    min_id: float,
    maxaccepts: int,
    strand: str,
    vsearch: str | Path,
    threads: int,
    use_udb: bool = False,
) -> Dict[str, str]:
    """Run VSEARCH global search.

    Outputs
    -------
    <output_prefix>.userout.tsv
        Custom tabular output.
    <output_prefix>.blast6.tsv
        BLAST6-like output.
    """

    threads = validate_threads(threads)
    strand = _require_strand(strand)

    if not 0 < min_id <= 1:
        raise ValueError("min_id must be > 0 and <= 1.")

    if maxaccepts < 1:
        raise ValueError("maxaccepts must be at least 1.")

    outdir = ensure_dir(outdir)

    userout = outdir / f"{output_prefix}.userout.tsv"
    blast6 = outdir / f"{output_prefix}.blast6.tsv"

    step(f"VSEARCH global search: {output_prefix}")

    command: List[str | Path] = [
        vsearch,
        "--usearch_global",
        input_fasta,
        "--id",
        str(min_id),
        "--maxaccepts",
        str(maxaccepts),
        "--strand",
        strand,
        "--userout",
        userout,
        "--blast6out",
        blast6,
        "--userfields",
        "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl",
        "--threads",
        str(threads),
    ]

    if use_udb:
        command.extend(["--db", db])
    else:
        command.extend(["--db", db])

    _run(command)

    return {
        "userout": str(userout),
        "blast6out": str(blast6),
    }


def make_udb(
    input_fasta: str | Path,
    output_udb: str | Path,
    vsearch: str | Path,
    threads: int,
) -> str:
    """Build a VSEARCH UDB database."""

    threads = validate_threads(threads)

    output_udb = Path(output_udb)
    output_udb.parent.mkdir(parents=True, exist_ok=True)

    step(f"Building VSEARCH UDB: {output_udb}")

    command = [
        vsearch,
        "--makeudb_usearch",
        input_fasta,
        "--output",
        output_udb,
        "--threads",
        str(threads),
    ]

    _run(command)
    return str(output_udb)


def vsearch_version(vsearch: str | Path) -> str:
    """Return the first line of `vsearch --version`."""

    completed = subprocess.run(
        [str(vsearch), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )

    output = (completed.stdout or completed.stderr or "").strip()

    if completed.returncode != 0:
        raise VsearchError(f"Failed to run {vsearch} --version")

    return output.splitlines()[0] if output else ""