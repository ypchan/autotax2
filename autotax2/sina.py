"""SINA integration and header parsing."""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable

import pandas as pd

from .fasta import FastaRecord, parse_fasta, strip_gaps_fasta
from .taxonomy import parse_tax_string
from .utils import ensure_dir, run_command

logger = logging.getLogger("autotax2.sina")

HEADER_FIELD_RE = re.compile(r"\[([^=\]]+)=([^\]]*)\]")


@dataclass
class SinaAnnotation:
    seq_id: str
    align_ident: float | None
    align_quality: float | None
    align_cutoff_head: int | None
    align_cutoff_tail: int | None
    lca_tax_ref: str | None
    lca_tax_gtdb: str | None
    lca_tax_slv: str | None
    turn: str | None
    raw_header: str


def _to_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _to_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def parse_sina_header(header: str) -> SinaAnnotation:
    """Parse SINA FASTA header fields.

    AutoTax2 treats `*_slv` fields as reference fields when SINA is run against
    user-provided `gtdb_ssu.arb`. Normalized `*_ref` fields are also accepted.
    """
    seq_id = header.split()[0]
    fields = {m.group(1): m.group(2) for m in HEADER_FIELD_RE.finditer(header)}

    align_ident = _to_float(fields.get("align_ident_ref") or fields.get("align_ident_slv"))
    align_quality = _to_float(fields.get("align_quality_ref") or fields.get("align_quality_slv"))
    head = _to_int(fields.get("align_cutoff_head_ref") or fields.get("align_cutoff_head_slv"))
    tail = _to_int(fields.get("align_cutoff_tail_ref") or fields.get("align_cutoff_tail_slv"))
    lca_ref = fields.get("lca_tax_ref") or fields.get("lca_tax_gtdb") or fields.get("lca_tax_slv")
    return SinaAnnotation(
        seq_id=seq_id,
        align_ident=align_ident,
        align_quality=align_quality,
        align_cutoff_head=head,
        align_cutoff_tail=tail,
        lca_tax_ref=lca_ref,
        lca_tax_gtdb=fields.get("lca_tax_gtdb"),
        lca_tax_slv=fields.get("lca_tax_slv"),
        turn=fields.get("turn"),
        raw_header=header,
    )


def run_sina(
    input_fa: str | Path,
    ref_arb: str | Path,
    output_fa: str | Path,
    threads: int,
    sina_bin: str = "sina",
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> Path:
    output_fa = Path(output_fa)
    ensure_dir(output_fa.parent)
    cmd = [
        sina_bin,
        "-i",
        str(input_fa),
        "-o",
        str(output_fa),
        "-r",
        str(ref_arb),
        "--threads",
        str(threads),
    ]
    if log_file:
        cmd += ["--log-file", str(log_file)]
    run_command(cmd, log_path=log_file, dry_run=dry_run)
    return output_fa


def parse_sina_fasta_to_table(sina_fa: str | Path) -> pd.DataFrame:
    rows = []
    for rec in parse_fasta(sina_fa):
        ann = parse_sina_header(rec.description)
        tax = parse_tax_string(ann.lca_tax_ref)
        rows.append(
            {
                "seq_id": ann.seq_id,
                "align_ident": ann.align_ident,
                "align_quality": ann.align_quality,
                "align_cutoff_head": ann.align_cutoff_head,
                "align_cutoff_tail": ann.align_cutoff_tail,
                "lca_tax_ref": ann.lca_tax_ref,
                "lca_tax_gtdb": ann.lca_tax_gtdb,
                "lca_tax_slv": ann.lca_tax_slv,
                "turn": ann.turn,
                "domain": tax.get("domain"),
                "phylum": tax.get("phylum"),
                "class": tax.get("class"),
                "order": tax.get("order"),
                "family": tax.get("family"),
                "genus": tax.get("genus"),
                "species": tax.get("species"),
                "raw_header": ann.raw_header,
            }
        )
    return pd.DataFrame(rows)


def strip_sina_alignment(sina_fa: str | Path, output_fa: str | Path) -> Path:
    strip_gaps_fasta(sina_fa, output_fa)
    return Path(output_fa)
