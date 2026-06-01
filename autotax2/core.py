"""Core AutoTax2 workflows."""
from __future__ import annotations

import logging
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
from .taxonomy import infer_anchor_rank, inherited_ranks, novel_ranks_below_anchor, parse_tax_string
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
    mode: str = "incremental",
    debug: bool = False,
    keep_temp: bool = False,
    dry_run: bool = False,
) -> None:
    """Add new sequences into an AutoTax2 database.

    This is an MVP implementation of the designed workflow:
    SINA -> parse anchors -> derep -> group by anchor -> local hierarchical clustering ->
    placeholder assignment -> registry update.
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

    derep_fa, derep_uc = derep_fulllength(
        corrected,
        version_dir / "02_derep" / "new.derep.fa",
        version_dir / "02_derep" / "new.derep.uc",
        threads=threads,
        vsearch_bin=cfg.get("tools", {}).get("vsearch", "vsearch"),
        log_file=version_dir / "02_derep" / "vsearch.derep.log",
        dry_run=dry_run,
    )

    # MVP grouping: group by anchor_rank + anchor taxonomy string; cluster all required novel ranks per group.
    # For production-scale runs, groups can be processed independently in parallel.
    cluster_summaries = []
    seq_registry_new = []
    seq_lengths = {rec.id: rec.length for rec in parse_fasta(corrected)}
    for rec in parse_fasta(corrected):
        row = sina_table[sina_table["seq_id"] == rec.id]
        if row.empty:
            continue
        r = row.iloc[0]
        seq_registry_new.append(
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

    # Create placeholders for novel ranks per sequence cluster group. This MVP clusters each broad group;
    # later optimization can batch groups and perform inheritance overlap checks.
    grouped = sina_table.groupby(["anchor_rank", "phylum", "class", "order", "family", "genus", "species"], dropna=False)
    for group_key, group_df in grouped:
        anchor_rank = group_key[0] or None
        novel_ranks = novel_ranks_below_anchor(anchor_rank)
        if not novel_ranks:
            continue
        group_ids = set(group_df["seq_id"])
        group_fa = version_dir / "03_clusters" / f"group_{len(cluster_summaries)+1}.fa"
        records = [rec for rec in parse_fasta(corrected) if rec.id in group_ids]
        write_fasta(records, group_fa)
        # Cluster fine-to-coarse novel ranks only: species -> genus -> ... up to rank below anchor.
        novel_ranks_fine = [r for r in ["species", "genus", "family", "order", "class", "phylum"] if r in novel_ranks]
        if novel_ranks_fine:
            res = hierarchical_cluster(
                group_fa,
                output_dir=version_dir / "03_clusters" / f"group_{len(cluster_summaries)+1}",
                ranks=novel_ranks_fine,
                thresholds_fraction=cfg.get("thresholds_fraction", DEFAULT_THRESHOLDS_FRACTION),
                threads=threads,
                iddef=cfg.get("vsearch", {}).get("iddef", 2),
                vsearch_bin=cfg.get("tools", {}).get("vsearch", "vsearch"),
                log_file=version_dir / "03_clusters" / "vsearch.cluster.log",
                dry_run=dry_run,
            )
            for cr in res:
                mem = uc_membership(cr.uc_file, cr.rank)
                mem_out = version_dir / "03_clusters" / f"group_{len(cluster_summaries)+1}.{cr.rank}.membership.tsv"
                save_table(mem, mem_out)
                cluster_summaries.append(
                    {
                        "group": str(len(cluster_summaries) + 1),
                        "anchor_rank": anchor_rank or "",
                        "rank": cr.rank,
                        "threshold": cr.threshold,
                        "uc_file": str(cr.uc_file),
                        "centroids_fa": str(cr.centroids_fa),
                        "n_members": len(mem),
                    }
                )

    if seq_registry_new:
        old = load_table(atdb.sequence_registry_path)
        seq_new_df = pd.DataFrame(seq_registry_new)
        save_table(pd.concat([old, seq_new_df], ignore_index=True).fillna(""), atdb.sequence_registry_path)

    if cluster_summaries:
        save_table(pd.DataFrame(cluster_summaries), version_dir / "cluster_summary.tsv")

    # Write sequence-level provisional taxonomy table using SINA anchor + placeholder allocation per missing rank.
    provisional = build_provisional_taxonomy(atdb, sina_table, source, prefix, version)
    save_table(provisional, version_dir / "provisional_taxonomy.tsv")
    update_current_taxonomy(atdb, provisional)
    logger.info("Finished add workflow. Version directory: %s", version_dir)


def build_provisional_taxonomy(
    atdb: AutoTaxDB, sina_table: pd.DataFrame, source: str, prefix: str, version: str
) -> pd.DataFrame:
    """Build a conservative sequence-level taxonomy table.

    MVP behavior: sequence-level placeholders are allocated for each missing novel rank.
    Cluster-level placeholder consolidation is represented in cluster_summary and can be
    expanded in later versions.
    """
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
