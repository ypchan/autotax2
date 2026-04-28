from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, List, Optional


SINTAX_RANKS = ["d", "p", "c", "o", "f", "g", "s"]


def parse_sintax_tax(tax: str) -> Dict[str, str]:
    """Parse VSEARCH SINTAX output taxonomy field.

    VSEARCH often writes tokens like d:Bacteria(1.0000),p:Firmicutes(0.9800).
    """
    out: Dict[str, str] = {}
    for token in tax.split(","):
        token = token.strip()
        if not token or ":" not in token:
            continue
        rank, value = token.split(":", 1)
        value = re.sub(r"\([0-9.]+\)$", "", value)
        out[rank] = value
    return out


def read_sintax(path: str | Path) -> Dict[str, Dict[str, str]]:
    data: Dict[str, Dict[str, str]] = {}
    if not path or not Path(path).exists():
        return data
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            query = parts[0]
            tax = parts[1]
            data[query] = parse_sintax_tax(tax)
    return data


def read_top_hit(path: str | Path) -> Dict[str, Dict[str, str]]:
    """Read first hit per query from --userout query+target+id+..."""
    data: Dict[str, Dict[str, str]] = {}
    if not path or not Path(path).exists():
        return data
    fields = ["query", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", "qcov", "tcov"]
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(fields, vals))
            q = row.get("query")
            if q and q not in data:
                data[q] = row
    return data


def format_taxonomy(rank_map: Dict[str, str]) -> str:
    parts = []
    for rank in SINTAX_RANKS:
        if rank in rank_map and rank_map[rank]:
            parts.append(f"{rank}:{rank_map[rank]}")
    return ",".join(parts)


def summarize(
    sintax_tsv: str | Path,
    silva_hits_tsv: str | Path,
    typestrain_hits_tsv: str | Path,
    output_tsv: str | Path,
) -> None:
    sintax = read_sintax(sintax_tsv)
    silva = read_top_hit(silva_hits_tsv)
    types = read_top_hit(typestrain_hits_tsv)

    queries = sorted(set(sintax) | set(silva) | set(types))
    with open(output_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "query",
            "sintax_taxonomy",
            "silva_top_target",
            "silva_top_identity",
            "silva_qcov",
            "silva_tcov",
            "typestrain_top_target",
            "typestrain_top_identity",
            "typestrain_qcov",
            "typestrain_tcov",
            "proposed_taxonomy",
            "proposed_source",
        ])
        for q in queries:
            sintax_tax = format_taxonomy(sintax.get(q, {}))
            sh = silva.get(q, {})
            th = types.get(q, {})

            # v0.1 conservative proposal:
            # prefer SINTAX taxonomy when available; hit data remains evidence columns.
            if sintax_tax:
                proposed = sintax_tax
                source = "sintax"
            elif sh:
                proposed = sh.get("target", "")
                source = "silva_top_hit"
            else:
                proposed = "unclassified"
                source = "none"

            writer.writerow([
                q,
                sintax_tax,
                sh.get("target", ""),
                sh.get("id", ""),
                sh.get("qcov", ""),
                sh.get("tcov", ""),
                th.get("target", ""),
                th.get("id", ""),
                th.get("qcov", ""),
                th.get("tcov", ""),
                proposed,
                source,
            ])
