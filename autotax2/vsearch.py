from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from .utils import ensure_dir, run_cmd, which_or_error


def vsearch_path(vsearch: str = "vsearch") -> str:
    if "/" in vsearch:
        return vsearch
    return which_or_error(vsearch)


def dereplicate(input_fasta: str, outdir: str, vsearch: str = "vsearch", threads: int = 8, dry_run: bool = False) -> str:
    ensure_dir(outdir)
    output = str(Path(outdir) / "derep.fa")
    uc = str(Path(outdir) / "derep.uc")
    cmd = [
        vsearch_path(vsearch),
        "--derep_fulllength", input_fasta,
        "--output", output,
        "--sizeout",
        "--uc", uc,
        "--threads", str(threads),
    ]
    run_cmd(cmd, dry_run=dry_run)
    return output


def sort_by_size(input_fasta: str, outdir: str, vsearch: str = "vsearch", minsize: Optional[int] = None, dry_run: bool = False) -> str:
    ensure_dir(outdir)
    output = str(Path(outdir) / "derep_sorted.fa")
    cmd = [vsearch_path(vsearch), "--sortbysize", input_fasta, "--output", output]
    if minsize is not None:
        cmd += ["--minsize", str(minsize)]
    run_cmd(cmd, dry_run=dry_run)
    return output


def cluster(
    input_fasta: str,
    outdir: str,
    identity: float,
    method: str = "cluster_size",
    vsearch: str = "vsearch",
    threads: int = 8,
    relabel: Optional[str] = None,
    dry_run: bool = False,
) -> dict:
    ensure_dir(outdir)
    label = str(identity).replace(".", "")
    centroids = str(Path(outdir) / f"otu{label}_centroids.fa")
    uc = str(Path(outdir) / f"otu{label}.uc")

    if method not in {"cluster_size", "cluster_fast", "cluster_smallmem"}:
        raise ValueError("--method must be cluster_size, cluster_fast, or cluster_smallmem")

    cmd = [
        vsearch_path(vsearch),
        f"--{method}", input_fasta,
        "--id", str(identity),
        "--centroids", centroids,
        "--uc", uc,
        "--threads", str(threads),
    ]
    if relabel:
        cmd += ["--relabel", relabel]
    run_cmd(cmd, dry_run=dry_run)
    return {"centroids": centroids, "uc": uc, "identity": identity, "method": method}


def sintax(
    query_fasta: str,
    sintax_db: str,
    outdir: str,
    cutoff: float = 0.8,
    vsearch: str = "vsearch",
    threads: int = 8,
    dry_run: bool = False,
) -> str:
    ensure_dir(outdir)
    out = str(Path(outdir) / "sintax.tsv")
    cmd = [
        vsearch_path(vsearch),
        "--sintax", query_fasta,
        "--db", sintax_db,
        "--sintax_cutoff", str(cutoff),
        "--tabbedout", out,
        "--threads", str(threads),
    ]
    run_cmd(cmd, dry_run=dry_run)
    return out


def usearch_global(
    query_fasta: str,
    db_fasta: str,
    outdir: str,
    prefix: str,
    identity: float = 0.70,
    maxaccepts: int = 10,
    strand: str = "both",
    vsearch: str = "vsearch",
    threads: int = 8,
    notmatched: bool = False,
    dry_run: bool = False,
) -> dict:
    ensure_dir(outdir)
    userout = str(Path(outdir) / f"{prefix}.tsv")
    uc = str(Path(outdir) / f"{prefix}.uc")
    cmd = [
        vsearch_path(vsearch),
        "--usearch_global", query_fasta,
        "--db", db_fasta,
        "--id", str(identity),
        "--strand", strand,
        "--maxaccepts", str(maxaccepts),
        "--maxrejects", "0",
        "--uc", uc,
        "--userout", userout,
        "--userfields", "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov+tcov",
        "--threads", str(threads),
    ]
    out = {"userout": userout, "uc": uc}
    if notmatched:
        nm = str(Path(outdir) / f"{prefix}_notmatched.fa")
        cmd += ["--notmatched", nm]
        out["notmatched"] = nm
    run_cmd(cmd, dry_run=dry_run)
    return out



def make_udb(input_fasta: str, output_udb: str, vsearch: str = "vsearch", dry_run: bool = False) -> str:
    """Build a USEARCH-compatible UDB with VSEARCH."""
    cmd = [vsearch_path(vsearch), "--makeudb_usearch", input_fasta, "--output", output_udb]
    run_cmd(cmd, dry_run=dry_run)
    return output_udb
