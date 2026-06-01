"""Export annotation references for DADA2, QIIME2, and vsearch SINTAX."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict

import pandas as pd

from .config import RANK_PREFIX, RANKS_COARSE_TO_FINE
from .fasta import FastaRecord, parse_fasta, write_fasta
from .registry import AutoTaxDB, load_table
from .utils import ensure_dir

logger = logging.getLogger("autotax2.export")


def _taxonomy_string(row: pd.Series, sep: str = ";", qiime_spaces: bool = False) -> str:
    parts = [row.get("domain", "")] + [row.get(rank, "") for rank in RANKS_COARSE_TO_FINE]
    if qiime_spaces:
        return "; ".join(parts)
    return sep.join(parts)


def _sintax_taxonomy(row: pd.Series) -> str:
    pairs = []
    for rank in ["domain"] + RANKS_COARSE_TO_FINE:
        taxon = str(row.get(rank, "") or "")
        if not taxon:
            continue
        pref = RANK_PREFIX[rank]
        clean = taxon
        if "__" in clean:
            clean = clean.split("__", 1)[1]
        pairs.append(f"{pref}:{clean}")
    return ",".join(pairs)


def export_database(db: str | Path, out: str | Path, fmt: str = "all") -> None:
    atdb = AutoTaxDB(Path(db))
    out = ensure_dir(out) if fmt == "all" or Path(out).suffix == "" else Path(out)
    tax = load_table(atdb.current_taxonomy)
    if tax.empty:
        raise ValueError(f"No current taxonomy found in {atdb.current_taxonomy}")
    seqs = {rec.id: rec for rec in parse_fasta(atdb.current_fa)} if atdb.current_fa.exists() else {}

    formats = ["sintax", "dada2", "qiime2", "taxonomy"] if fmt == "all" else [fmt]
    for f in formats:
        if f == "taxonomy":
            path = (out / "autotax2.taxonomy.tsv") if isinstance(out, Path) and out.suffix == "" else out
            tax.to_csv(path, sep="\t", index=False)
            logger.info("Wrote taxonomy table: %s", path)
        elif f == "sintax":
            path = out / "autotax2.sintax.fa" if out.is_dir() else out
            export_sintax(tax, seqs, path)
        elif f == "dada2":
            path = out / "autotax2.dada2.fa" if out.is_dir() else out
            export_dada2(tax, seqs, path)
        elif f == "qiime2":
            out_dir = ensure_dir(out if out.is_dir() else out.parent)
            export_qiime2(tax, seqs, out_dir)
        else:
            raise ValueError(f"Unknown export format: {f}")


def export_sintax(tax: pd.DataFrame, seqs: Dict[str, FastaRecord], path: str | Path) -> None:
    records = []
    for _, row in tax.iterrows():
        sid = row["seq_id"]
        if sid not in seqs:
            continue
        sintax = _sintax_taxonomy(row)
        records.append(FastaRecord(sid, f"{sid};tax={sintax};", seqs[sid].sequence))
    write_fasta(records, path)
    logger.info("Wrote SINTAX FASTA: %s", path)


def export_dada2(tax: pd.DataFrame, seqs: Dict[str, FastaRecord], path: str | Path) -> None:
    records = []
    for _, row in tax.iterrows():
        sid = row["seq_id"]
        if sid not in seqs:
            continue
        header = _taxonomy_string(row, sep=";")
        records.append(FastaRecord(sid, header, seqs[sid].sequence))
    write_fasta(records, path)
    logger.info("Wrote DADA2 FASTA: %s", path)


def export_qiime2(tax: pd.DataFrame, seqs: Dict[str, FastaRecord], out_dir: Path) -> None:
    seq_path = out_dir / "autotax2.qiime2-seqs.fa"
    tax_path = out_dir / "autotax2.qiime2-taxonomy.tsv"
    records = [seqs[sid] for sid in tax["seq_id"] if sid in seqs]
    write_fasta(records, seq_path)
    qiime = pd.DataFrame(
        {
            "Feature ID": tax["seq_id"],
            "Taxon": tax.apply(lambda r: _taxonomy_string(r, qiime_spaces=True), axis=1),
        }
    )
    qiime.to_csv(tax_path, sep="\t", index=False)
    logger.info("Wrote QIIME2 sequence FASTA: %s", seq_path)
    logger.info("Wrote QIIME2 taxonomy TSV: %s", tax_path)
