"""vsearch wrappers and UC parsing."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable

import pandas as pd

from .config import DEFAULT_THRESHOLDS_FRACTION, RANKS_COARSE_TO_FINE
from .utils import ensure_dir, run_command

logger = logging.getLogger("autotax2.vsearch")


@dataclass
class ClusterResult:
    rank: str
    input_fa: Path
    centroids_fa: Path
    uc_file: Path
    threshold: float


def derep_fulllength(
    input_fa: str | Path,
    output_fa: str | Path,
    uc_file: str | Path,
    threads: int,
    vsearch_bin: str = "vsearch",
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> tuple[Path, Path]:
    output_fa = Path(output_fa)
    uc_file = Path(uc_file)
    ensure_dir(output_fa.parent)
    cmd = [
        vsearch_bin,
        "--derep_fulllength",
        str(input_fa),
        "--output",
        str(output_fa),
        "--uc",
        str(uc_file),
        "--sizeout",
        "--threads",
        str(threads),
    ]
    run_command(cmd, log_path=log_file, dry_run=dry_run)
    return output_fa, uc_file


def cluster_fast(
    input_fa: str | Path,
    rank: str,
    identity_fraction: float,
    output_dir: str | Path,
    threads: int,
    iddef: int = 2,
    vsearch_bin: str = "vsearch",
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> ClusterResult:
    output_dir = ensure_dir(output_dir)
    centroids = output_dir / f"{rank}.centroids.fa"
    uc = output_dir / f"{rank}.uc"
    cmd = [
        vsearch_bin,
        "--cluster_fast",
        str(input_fa),
        "--id",
        str(identity_fraction),
        "--iddef",
        str(iddef),
        "--strand",
        "plus",
        "--centroids",
        str(centroids),
        "--uc",
        str(uc),
        "--threads",
        str(threads),
    ]
    run_command(cmd, log_path=log_file, dry_run=dry_run)
    return ClusterResult(rank, Path(input_fa), centroids, uc, identity_fraction)


def hierarchical_cluster(
    input_fa: str | Path,
    output_dir: str | Path,
    ranks: list[str] | None,
    thresholds_fraction: Dict[str, float] | None,
    threads: int,
    iddef: int = 2,
    vsearch_bin: str = "vsearch",
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> list[ClusterResult]:
    """Run strict hierarchical clustering.

    For a full hierarchy, input is all sequences for species, then species
    centroids for genus, genus centroids for family, and so on.
    """
    thresholds_fraction = thresholds_fraction or DEFAULT_THRESHOLDS_FRACTION
    ranks = ranks or ["species", "genus", "family", "order", "class", "phylum"]
    current_input = Path(input_fa)
    results: list[ClusterResult] = []
    for rank in ranks:
        res = cluster_fast(
            current_input,
            rank=rank,
            identity_fraction=thresholds_fraction[rank],
            output_dir=output_dir,
            threads=threads,
            iddef=iddef,
            vsearch_bin=vsearch_bin,
            log_file=log_file,
            dry_run=dry_run,
        )
        results.append(res)
        current_input = res.centroids_fa
    return results


def read_uc(path: str | Path) -> pd.DataFrame:
    """Read vsearch/usearch UC file.

    Columns follow UC format: type, cluster, size, pctid, strand, qlo, tlo,
    aln, query, target.
    """
    cols = ["record_type", "cluster", "size", "pctid", "strand", "qlo", "tlo", "aln", "query", "target"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, comment="#", dtype=str)
    return df


def uc_membership(path: str | Path, rank: str) -> pd.DataFrame:
    df = read_uc(path)
    df = df[df["record_type"].isin(["S", "H"])].copy()
    df["centroid"] = df.apply(
        lambda r: r["query"] if r["record_type"] == "S" or r["target"] == "*" else r["target"], axis=1
    )
    df["member"] = df["query"]
    df["rank"] = rank
    return df[["rank", "cluster", "member", "centroid", "record_type", "pctid"]]
