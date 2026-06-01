"""Core AutoTax2 workflows."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import os
from pathlib import Path
from typing import Dict

import pandas as pd

from .config import DEFAULT_THRESHOLDS_FRACTION, DEFAULT_THRESHOLDS_PERCENT, RANKS_COARSE_TO_FINE
from .fasta import FastaRecord, parse_fasta, write_fasta
from .registry import (
    AutoTaxDB,
    allocate_placeholder,
    create_initial_registries,
    load_table,
    register_source,
    save_table,
)
from .sina import parse_sina_fasta_to_table, run_sina, strip_sina_alignment
from .taxonomy import infer_anchor_rank, inherited_ranks, novel_ranks_below_anchor
from .utils import ensure_dir, require_executable, require_file
from .vsearch import derep_fulllength, hierarchical_cluster, uc_membership

logger = logging.getLogger("autotax2.core")


def check_tools(sina_bin: str = "sina", vsearch_bin: str = "vsearch") -> dict[str, str]:
    sina = require_executable(sina_bin)
    vsearch = require_executable(vsearch_bin)
    logger.info("Found SINA: %s", sina)
    logger.info("Found vsearch: %s", vsearch)
    return {"sina": sina, "vsearch": vsearch}


def init_db(
    ref_fa: str | Path,
    ref_tax: str | Path,
    ref_arb: str | Path,
    db: str | Path,
    threads: int,
    force: bool = False,
) -> None:
    ref_fa = require_file(ref_fa, "reference FASTA")
    ref_tax = require_file(ref_tax, "reference taxonomy TSV")
    ref_arb = require_file(ref_arb, "reference ARB")
    atdb = AutoTaxDB(Path(db))
    create_initial_registries(atdb, ref_fa, ref_tax, ref_arb, threads=threads, force=force)


def add_sequences(
    db: str | Path,
    input_fa: str | Path,
    source: str,
    prefix: str,
    threads: int,
    group_jobs: int | None = None,
    mode: str = "incremental",
    debug: bool = False,
    keep_temp: bool = False,
    dry_run: bool = False,
) -> None:
    """Add new sequences into an AutoTax2 database.

    SINA anchors the whole dataset, then independent anchor groups are clustered
    concurrently within the user's total thread budget.
    """
    if mode not in {"incremental", "full"}:
        raise ValueError("mode must be 'incremental' or 'full'")
    atdb = AutoTaxDB(Path(db))
    cfg = atdb.read_config()
    ref_arb = Path(cfg["reference"]["ref_arb"])
    input_fa = require_file(input_fa, "input FASTA")
    version = atdb.next_version_name(source)
    version_dir = ensure_dir(atdb.versions_dir / version)
    ensure_dir(version_dir / "01_sina")
    ensure_dir(version_dir / "02_derep")
    ensure_dir(version_dir / "03_clusters")
    register_source(atdb, source=source, prefix=prefix, version=version)

    logger.info("Adding dataset source=%s prefix=%s version=%s mode=%s", source, prefix, version, mode)

    sina_aln = version_dir / "01_sina" / "new.sina.aligned.fa"
    sina_log = version_dir / "01_sina" / "sina.log"
    run_sina(
        input_fa=input_fa,
        ref_arb=ref_arb,
        output_fa=sina_aln,
        threads=threads,
        sina_bin=cfg.get("tools", {}).get("sina", "sina"),
        log_file=sina_log,
        dry_run=dry_run,
    )
    corrected = version_dir / "01_sina" / "new.sina.corrected.fa"
    if not dry_run:
        strip_sina_alignment(sina_aln, corrected)
        sina_table = parse_sina_fasta_to_table(sina_aln)
        thresholds = cfg.get("thresholds_percent", DEFAULT_THRESHOLDS_PERCENT)
        sina_table["anchor_rank"] = sina_table.apply(
            lambda r: infer_anchor_rank(
                r.get("align_ident"),
                {rank: r.get(rank) for rank in ["domain"] + RANKS_COARSE_TO_FINE},
                thresholds,
            )
            or "",
            axis=1,
        )
        sina_table["novel_ranks"] = sina_table["anchor_rank"].apply(
            lambda x: ",".join(novel_ranks_below_anchor(x or None))
        )
        save_table(sina_table, version_dir / "sina_annotation.tsv")
    else:
        logger.info("Dry run: skipping downstream table generation.")
        return

    derep_fulllength(
        corrected,
        version_dir / "02_derep" / "new.derep.fa",
        version_dir / "02_derep" / "new.derep.uc",
        threads=threads,
        vsearch_bin=cfg.get("tools", {}).get("vsearch", "vsearch"),
        log_file=version_dir / "02_derep" / "vsearch.derep.log",
        dry_run=dry_run,
    )

    corrected_records = list(parse_fasta(corrected))
    seq_registry_new = build_sequence_registry_rows(
        records=corrected_records,
        sina_table=sina_table,
        source=source,
        version=version,
    )

    cluster_groups = build_cluster_groups(sina_table)
    cluster_summaries = []
    if cluster_groups:
        workers, threads_per_worker = resolve_group_parallelism(
            group_count=len(cluster_groups), threads=threads, group_jobs=group_jobs
        )
        logger.info(
            "Clustering %d anchor groups with %d worker(s), %d thread(s) per worker",
            len(cluster_groups),
            workers,
            threads_per_worker,
        )
        cluster_summaries = run_cluster_groups(
            groups=cluster_groups,
            records=corrected_records,
            version_dir=version_dir,
            cfg=cfg,
            threads_per_worker=threads_per_worker,
            workers=workers,
            dry_run=dry_run,
        )

    if seq_registry_new:
        old = load_table(atdb.sequence_registry_path)
        seq_new_df = pd.DataFrame(seq_registry_new)
        save_table(pd.concat([old, seq_new_df], ignore_index=True).fillna(""), atdb.sequence_registry_path)

    if cluster_summaries:
        save_table(pd.DataFrame(cluster_summaries), version_dir / "cluster_summary.tsv")

    provisional = build_provisional_taxonomy(atdb, sina_table, source, prefix, version)
    save_table(provisional, version_dir / "provisional_taxonomy.tsv")
    update_current_taxonomy(atdb, provisional)
    logger.info("Finished add workflow. Version directory: %s", version_dir)


def build_sequence_registry_rows(
    records: list[FastaRecord], sina_table: pd.DataFrame, source: str, version: str
) -> list[dict]:
    rows = []
    by_seq_id = {row["seq_id"]: row for _, row in sina_table.iterrows()}
    for rec in records:
        r = by_seq_id.get(rec.id)
        if r is None:
            continue
        rows.append(
            {
                "seq_id": rec.id,
                "seq_md5": rec.md5,
                "source": source,
                "source_seq_id": rec.id,
                "length": rec.length,
                "added_version": version,
                "status": "active",
                "align_ident": r.get("align_ident", ""),
                "align_quality": r.get("align_quality", ""),
                "anchor_rank": r.get("anchor_rank", ""),
            }
        )
    return rows


def build_cluster_groups(sina_table: pd.DataFrame) -> list[dict]:
    groups = []
    grouped = sina_table.groupby(
        ["anchor_rank", "phylum", "class", "order", "family", "genus", "species"],
        dropna=False,
    )
    for group_key, group_df in grouped:
        anchor_rank = group_key[0] or None
        novel_ranks = novel_ranks_below_anchor(anchor_rank)
        novel_ranks_fine = [
            r for r in ["species", "genus", "family", "order", "class", "phylum"] if r in novel_ranks
        ]
        if not novel_ranks_fine:
            continue
        groups.append(
            {
                "group_index": len(groups) + 1,
                "anchor_rank": anchor_rank,
                "seq_ids": set(group_df["seq_id"]),
                "novel_ranks_fine": novel_ranks_fine,
            }
        )
    return groups


def resolve_group_parallelism(
    group_count: int, threads: int, group_jobs: int | None = None
) -> tuple[int, int]:
    """Choose worker count without exceeding the requested CPU budget."""
    if group_count <= 0:
        return 0, 0
    total_threads = max(1, int(threads or 1))
    if group_jobs is None:
        cpu_count = os.cpu_count() or total_threads
        workers = min(group_count, total_threads, cpu_count)
    else:
        workers = min(group_count, total_threads, max(1, int(group_jobs)))
    threads_per_worker = max(1, total_threads // workers)
    return workers, threads_per_worker


def run_cluster_groups(
    groups: list[dict],
    records: list[FastaRecord],
    version_dir: Path,
    cfg: dict,
    threads_per_worker: int,
    workers: int,
    dry_run: bool,
) -> list[dict]:
    """Cluster independent anchor groups concurrently and return deterministic summaries."""
    if workers <= 1:
        summaries = []
        for group in groups:
            summaries.extend(
                cluster_one_group(
                    group=group,
                    records=records,
                    version_dir=version_dir,
                    cfg=cfg,
                    threads=threads_per_worker,
                    dry_run=dry_run,
                )
            )
        return summaries

    results: dict[int, list[dict]] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                cluster_one_group,
                group=group,
                records=records,
                version_dir=version_dir,
                cfg=cfg,
                threads=threads_per_worker,
                dry_run=dry_run,
            ): group["group_index"]
            for group in groups
        }
        for future in as_completed(futures):
            group_index = futures[future]
            results[group_index] = future.result()

    summaries = []
    for group_index in sorted(results):
        summaries.extend(results[group_index])
    return summaries


def cluster_one_group(
    group: dict,
    records: list[FastaRecord],
    version_dir: Path,
    cfg: dict,
    threads: int,
    dry_run: bool,
) -> list[dict]:
    group_index = int(group["group_index"])
    group_ids = group["seq_ids"]
    cluster_dir = version_dir / "03_clusters"
    group_fa = cluster_dir / f"group_{group_index}.fa"
    group_records = [rec for rec in records if rec.id in group_ids]
    write_fasta(group_records, group_fa)

    results = hierarchical_cluster(
        group_fa,
        output_dir=cluster_dir / f"group_{group_index}",
        ranks=group["novel_ranks_fine"],
        thresholds_fraction=cfg.get("thresholds_fraction", DEFAULT_THRESHOLDS_FRACTION),
        threads=threads,
        iddef=cfg.get("vsearch", {}).get("iddef", 2),
        vsearch_bin=cfg.get("tools", {}).get("vsearch", "vsearch"),
        log_file=cluster_dir / f"group_{group_index}.vsearch.cluster.log",
        dry_run=dry_run,
    )

    summaries = []
    for cr in results:
        mem = uc_membership(cr.uc_file, cr.rank)
        mem_out = cluster_dir / f"group_{group_index}.{cr.rank}.membership.tsv"
        save_table(mem, mem_out)
        summaries.append(
            {
                "group": str(group_index),
                "anchor_rank": group["anchor_rank"] or "",
                "rank": cr.rank,
                "threshold": cr.threshold,
                "uc_file": str(cr.uc_file),
                "centroids_fa": str(cr.centroids_fa),
                "n_members": len(mem),
            }
        )
    return summaries


def build_provisional_taxonomy(
    atdb: AutoTaxDB, sina_table: pd.DataFrame, source: str, prefix: str, version: str
) -> pd.DataFrame:
    """Build a conservative sequence-level taxonomy table."""
    rows = []
    for _, r in sina_table.iterrows():
        anchor_rank = r.get("anchor_rank") or None
        tax = {rank: (r.get(rank) if isinstance(r.get(rank), str) else "") for rank in ["domain"] + RANKS_COARSE_TO_FINE}
        inherited = set(inherited_ranks(anchor_rank))
        for rank in RANKS_COARSE_TO_FINE:
            if rank not in inherited:
                placeholder, num = allocate_placeholder(atdb, rank, prefix)
                tax[rank] = placeholder
        rows.append(
            {
                "seq_id": r["seq_id"],
                "source": source,
                "added_version": version,
                "anchor_rank": anchor_rank or "",
                "align_ident": r.get("align_ident", ""),
                "align_quality": r.get("align_quality", ""),
                **tax,
            }
        )
    return pd.DataFrame(rows)


def update_current_taxonomy(atdb: AutoTaxDB, new_tax: pd.DataFrame) -> None:
    current = load_table(atdb.current_taxonomy)
    merged = pd.concat([current, new_tax], ignore_index=True).fillna("")
    save_table(merged, atdb.current_taxonomy)


def rebuild_db(db: str | Path, threads: int, dry_run: bool = False) -> None:
    atdb = AutoTaxDB(Path(db))
    cfg = atdb.read_config()
    logger.info("Running full hierarchical rebuild for %s", atdb.root)
    results = hierarchical_cluster(
        atdb.current_fa,
        output_dir=atdb.centroids_dir,
        ranks=["species", "genus", "family", "order", "class", "phylum"],
        thresholds_fraction=DEFAULT_THRESHOLDS_FRACTION,
        threads=threads,
        iddef=cfg.get("vsearch", {}).get("iddef", 2),
        vsearch_bin=cfg.get("tools", {}).get("vsearch", "vsearch"),
        log_file=atdb.logs_dir / "rebuild.vsearch.log",
        dry_run=dry_run,
    )
    logger.info("Rebuild completed with %d rank clustering steps", len(results))
