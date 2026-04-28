#!/usr/bin/env python3
from __future__ import annotations

import argparse

from autotax2.intron import detect_introns
from autotax2.threads import resolve_threads


def main() -> None:
    p = argparse.ArgumentParser(description="Standalone intron detector for near/full-length rRNA sequences.")
    p.add_argument("--input", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--source-label", default="query")
    p.add_argument("--vsearch", default="vsearch")
    p.add_argument("--threads", default="auto")
    p.add_argument("--search-id", type=float, default=0.70)
    p.add_argument("--rescue-id", type=float, default=0.987)
    p.add_argument("--min-intron-len", type=int, default=50)
    p.add_argument("--min-flank-len", type=int, default=150)
    p.add_argument("--strand", default="both", choices=["both", "plus"])
    p.add_argument("--maxaccepts", type=int, default=1)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    outputs = detect_introns(
        input_fasta=args.input,
        db=args.db,
        outdir=args.out,
        source_label=args.source_label,
        vsearch=args.vsearch,
        threads=resolve_threads(args.threads),
        search_id=args.search_id,
        rescue_id=args.rescue_id,
        min_intron_len=args.min_intron_len,
        min_flank_len=args.min_flank_len,
        strand=args.strand,
        maxaccepts=args.maxaccepts,
        dry_run=args.dry_run,
    )
    for k, v in outputs.items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
