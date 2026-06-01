"""Command-line interface for AutoTax2."""
from __future__ import annotations

from pathlib import Path
import logging

import typer
from rich.console import Console

from .core import add_sequences, check_tools, init_db, rebuild_db
from .export import export_database
from .logging_config import setup_logging
from .summarize import summarize_sources

app = typer.Typer(
    help=(
        "AutoTax2: GTDB-anchored incremental SSU taxonomy construction with "
        "SINA and vsearch."
    ),
    no_args_is_help=True,
)
console = Console()


@app.command()
def check(
    sina_bin: str = typer.Option("sina", help="SINA executable name or path."),
    vsearch_bin: str = typer.Option("vsearch", help="vsearch executable name or path."),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging."),
):
    """Check external software availability.

    Example:
        autotax2 check --debug
    """
    setup_logging(debug=debug)
    tools = check_tools(sina_bin=sina_bin, vsearch_bin=vsearch_bin)
    console.print("[green]AutoTax2 dependency check passed.[/green]")
    for name, path in tools.items():
        console.print(f"  {name}: {path}")


@app.command()
def init(
    ref_fa: Path = typer.Option(..., "--ref-fa", help="Required. GTDB SSU reference FASTA."),
    ref_tax: Path = typer.Option(
        ..., "--ref-tax", help="Required. GTDB taxonomy TSV with seq_id and 7 ranks."
    ),
    ref_arb: Path = typer.Option(..., "--ref-arb", help="Required. GTDB SSU ARB file."),
    db: Path = typer.Option(..., "--db", help="Required. AutoTax2 database directory."),
    threads: int = typer.Option(1, "--threads", "-t", help="Optional. CPU threads."),
    force: bool = typer.Option(False, "--force", help="Overwrite an existing database directory."),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging."),
):
    """Initialize an AutoTax2 database from a GTDB-derived SSU backbone.

    Required:
      --ref-fa, --ref-tax, --ref-arb, --db

    Optional:
      --threads, --force, --debug

    Example:
      autotax2 init \
        --ref-fa gtdb_ssu_tax.fa \
        --ref-tax gtdb_ssu.taxonomy.tsv \
        --ref-arb gtdb_ssu.arb \
        --db autotax2_db \
        --threads 64
    """
    setup_logging(log_file=db / "logs" / "init.log", debug=debug)
    init_db(ref_fa=ref_fa, ref_tax=ref_tax, ref_arb=ref_arb, db=db, threads=threads, force=force)
    console.print(f"[green]Initialized AutoTax2 database:[/green] {db}")


@app.command()
def add(
    db: Path = typer.Option(..., "--db", help="Required. Existing AutoTax2 database."),
    input_fa: Path = typer.Option(..., "--input", help="Required. New SSU FASTA."),
    source: str = typer.Option(..., "--source", help="Required. Source name, e.g. midas."),
    prefix: str = typer.Option(..., "--prefix", help="Required. Placeholder prefix for new lineages."),
    threads: int = typer.Option(1, "--threads", "-t", help="Optional. CPU threads."),
    mode: str = typer.Option(
        "incremental", "--mode", help="Optional. incremental or full. Default: incremental."
    ),
    keep_temp: bool = typer.Option(False, "--keep-temp", help="Keep intermediate files."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate and print commands without running."),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging."),
):
    """Add a new SSU dataset into an AutoTax2 database.

    Required:
      --db, --input, --source, --prefix

    Optional:
      --threads, --mode, --keep-temp, --dry-run, --debug

    Example:
      autotax2 add \
        --db autotax2_db \
        --input MiDAS5.3.ssu.fa \
        --source midas \
        --prefix midas \
        --threads 64 \
        --mode incremental
    """
    setup_logging(log_file=db / "logs" / f"add_{source}.log", debug=debug)
    add_sequences(
        db=db,
        input_fa=input_fa,
        source=source,
        prefix=prefix,
        threads=threads,
        mode=mode,
        keep_temp=keep_temp,
        dry_run=dry_run,
        debug=debug,
    )
    console.print(f"[green]Added source {source} into database:[/green] {db}")


@app.command()
def rebuild(
    db: Path = typer.Option(..., "--db", help="Required. AutoTax2 database."),
    threads: int = typer.Option(1, "--threads", "-t", help="Optional. CPU threads."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print commands without running."),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging."),
):
    """Run a full hierarchical vsearch rebuild.

    Example:
      autotax2 rebuild --db autotax2_db --threads 128
    """
    setup_logging(log_file=db / "logs" / "rebuild.log", debug=debug)
    rebuild_db(db=db, threads=threads, dry_run=dry_run)
    console.print(f"[green]Rebuild finished for database:[/green] {db}")


@app.command("export")
def export_cmd(
    db: Path = typer.Option(..., "--db", help="Required. AutoTax2 database."),
    fmt: str = typer.Option(
        "all", "--format", help="Required/optional. all, sintax, dada2, qiime2, or taxonomy."
    ),
    out: Path = typer.Option(..., "--out", help="Required. Output file or directory."),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging."),
):
    """Export annotation references.

    Required:
      --db, --format, --out

    Export formats:
      all, sintax, dada2, qiime2, taxonomy

    Example:
      autotax2 export --db autotax2_db --format all --out export/
    """
    setup_logging(log_file=db / "logs" / "export.log", debug=debug)
    export_database(db=db, out=out, fmt=fmt)
    console.print(f"[green]Export finished:[/green] {out}")


@app.command()
def summarize(
    db: Path = typer.Option(..., "--db", help="Required. AutoTax2 database."),
    rank: str = typer.Option("species", "--rank", help="Rank to summarize."),
    out: Path | None = typer.Option(None, "--out", help="Optional output TSV."),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging."),
):
    """Summarize source overlap for a rank.

    Example:
      autotax2 summarize --db autotax2_db --rank species --out species.summary.tsv
    """
    setup_logging(log_file=db / "logs" / "summarize.log", debug=debug)
    summary = summarize_sources(db=db, rank=rank, out=out)
    console.print(summary.head(20).to_string(index=False))
    if out:
        console.print(f"[green]Wrote summary:[/green] {out}")


if __name__ == "__main__":
    app()
