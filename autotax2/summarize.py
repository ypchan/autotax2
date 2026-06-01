"""Summary utilities for source overlap."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from .config import RANKS_FINE_TO_COARSE
from .registry import AutoTaxDB, load_table, save_table


def summarize_sources(db: str | Path, rank: str, out: str | Path | None = None) -> pd.DataFrame:
    if rank not in RANKS_FINE_TO_COARSE:
        raise ValueError(f"rank must be one of {RANKS_FINE_TO_COARSE}")
    atdb = AutoTaxDB(Path(db))
    mem = load_table(atdb.clusters_dir / f"{rank}.membership.tsv")
    if mem.empty:
        tax = load_table(atdb.current_taxonomy)
        if tax.empty:
            raise ValueError("No membership or taxonomy table available")
        if rank not in tax.columns:
            raise ValueError(f"Rank column not found: {rank}")
        grp = tax.groupby(rank)
        rows = []
        for taxon, sub in grp:
            sources = sorted(set(sub.get("source", pd.Series(dtype=str))))
            rows.append(
                {
                    "rank": rank,
                    "taxon_name": taxon,
                    "n_sequences": len(sub),
                    "n_sources": len(sources),
                    "sources": ",".join(sources),
                    "source_pattern": "+".join(sources),
                    "is_source_specific": "yes" if len(sources) == 1 else "no",
                }
            )
        summary = pd.DataFrame(rows)
    else:
        rows = []
        for taxon, sub in mem.groupby("taxon_name"):
            sources = sorted(set(sub["source"]))
            rows.append(
                {
                    "rank": rank,
                    "taxon_name": taxon,
                    "n_sequences": len(sub),
                    "n_sources": len(sources),
                    "sources": ",".join(sources),
                    "source_pattern": "+".join(sources),
                    "is_source_specific": "yes" if len(sources) == 1 else "no",
                }
            )
        summary = pd.DataFrame(rows)
    if out:
        save_table(summary, Path(out))
    return summary
