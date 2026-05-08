"""Command-line interface for autotax2."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from autotax2.barrnap import prepare_dataset
from autotax2.config import load_config
from autotax2.export import export_references
from autotax2.placement import place_dataset
from autotax2.reports import summarize_build
from autotax2.sina import orient_dataset_with_sina
from autotax2.silva import initialize_silva_build, resolve_silva_unresolved_build
from autotax2.validate import validate_build
from autotax2.vsearch import cluster_search_dataset

app = typer.Typer(
    name="autotax2",
    help="Fixed-backbone, rank-aware, incremental rRNA gene reference builder.",
    no_args_is_help=True,
)
console = Console()


class DomainChoice(str, Enum):
    """Supported dataset domains."""

    ARCHAEA = "Archaea"
    BACTERIA = "Bacteria"


class MultiRrnaPolicy(str, Enum):
    """Supported policies for multiple barrnap rRNA hits."""

    LONGEST = "longest"
    BEST = "best"
    ALL = "all"
    FAIL = "fail"


class StrandChoice(str, Enum):
    """Supported VSEARCH strand modes."""

    PLUS = "plus"
    BOTH = "both"


@app.command()
def init(
    silva_fasta: Path = typer.Option(
        ...,
        "--silva-fasta",
        help="SILVA taxonomy FASTA file, plain or gzipped.",
    ),
    outdir: Path = typer.Option(
        ...,
        "--outdir",
        help="Output build directory to initialize.",
    ),
    domain: Optional[str] = typer.Option(
        None,
        "--domain",
        help="Optional SILVA domain filter, for example Archaea or Bacteria.",
    ),
    type_strain_metadata: Optional[Path] = typer.Option(
        None,
        "--type-strain-metadata",
        help="Optional type-strain metadata TSV keyed by seq_id.",
    ),
) -> None:
    """Initialize a build directory from a SILVA FASTA backbone."""
    summary = initialize_silva_build(
        silva_fasta=silva_fasta,
        outdir=outdir,
        domain=domain,
        type_strain_metadata=type_strain_metadata,
    )
    console.print(
        "[green]Initialized autotax2 build:[/green] "
        f"{outdir} records={summary['records']} "
        f"named={summary['named']} unresolved={summary['unresolved']}"
    )


@app.command("resolve-silva")
def resolve_silva(
    build: Path = typer.Option(
        ...,
        "--build",
        help="Initialized autotax2 build directory.",
    ),
    threads: int = typer.Option(4, "--threads", help="Threads for VSEARCH clustering."),
    species_id: float = typer.Option(0.987, "--species-id", help="Species cluster identity."),
    genus_id: float = typer.Option(0.945, "--genus-id", help="Genus cluster identity."),
    family_id: float = typer.Option(0.865, "--family-id", help="Family cluster identity."),
    order_id: float = typer.Option(0.820, "--order-id", help="Order cluster identity."),
    class_id: float = typer.Option(0.785, "--class-id", help="Class cluster identity."),
    floor_id: float = typer.Option(0.750, "--floor-id", help="Minimum clustering identity floor."),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin", help="VSEARCH executable."),
    iddef: int = typer.Option(2, "--iddef", help="VSEARCH --iddef value."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Write proposals without registry updates."),
) -> None:
    """Resolve SILVA unresolved records into a placeholder scaffold."""
    summary = resolve_silva_unresolved_build(
        build=build,
        threads=threads,
        species_id=species_id,
        genus_id=genus_id,
        family_id=family_id,
        order_id=order_id,
        class_id=class_id,
        floor_id=floor_id,
        vsearch_bin=vsearch_bin,
        iddef=iddef,
        dry_run=dry_run,
    )
    mode = "dry-run " if summary.dry_run else ""
    console.print(
        f"[green]Completed {mode}SILVA unresolved resolution:[/green] "
        f"{build} unresolved={summary.unresolved_records} "
        f"placeholder_taxa={summary.placeholder_taxa}"
    )


@app.command("prepare-dataset")
def prepare_dataset_command(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    name: str = typer.Option(..., "--name", help="Dataset name."),
    prefix: str = typer.Option(..., "--prefix", help="Frozen dataset prefix."),
    fasta: Path = typer.Option(..., "--fasta", help="Input intron-free FASTA file."),
    domain: DomainChoice = typer.Option(..., "--domain", help="Dataset domain."),
    threads: int = typer.Option(4, "--threads", help="Threads for barrnap."),
    barrnap_bin: str = typer.Option("barrnap", "--barrnap-bin", help="barrnap executable."),
    barrnap_kingdom: Optional[str] = typer.Option(
        None,
        "--barrnap-kingdom",
        help="Optional barrnap kingdom override: arc or bac.",
    ),
    strict_tool_version: bool = typer.Option(
        False,
        "--strict-tool-version",
        help="Fail unless barrnap 1.10.5 is detected.",
    ),
    multi_rrna_policy: MultiRrnaPolicy = typer.Option(
        MultiRrnaPolicy.LONGEST,
        "--multi-rrna-policy",
        help="Policy for sequences with multiple rRNA hits.",
    ),
    min_rrna_len_archaea: int = typer.Option(
        900,
        "--min-rrna-len-archaea",
        help="Minimum extracted archaeal rRNA length.",
    ),
    min_rrna_len_bacteria: int = typer.Option(
        1200,
        "--min-rrna-len-bacteria",
        help="Minimum extracted bacterial rRNA length.",
    ),
    flank: int = typer.Option(0, "--flank", help="Extra bases to include around barrnap hit."),
    allow_partial: bool = typer.Option(False, "--allow-partial", help="Allow partial barrnap hits."),
    reject_non_atgc: bool = typer.Option(
        True,
        "--reject-non-atgc/--no-reject-non-atgc",
        help="Reject normalized input sequences containing non-ATGC symbols.",
    ),
) -> None:
    """Prepare a custom dataset with mandatory barrnap recutting."""
    summary = prepare_dataset(
        build=build,
        name=name,
        prefix=prefix,
        fasta=fasta,
        domain=domain.value,
        threads=threads,
        barrnap_bin=barrnap_bin,
        barrnap_kingdom=barrnap_kingdom,
        strict_tool_version=strict_tool_version,
        multi_rrna_policy=multi_rrna_policy.value,
        min_rrna_len_archaea=min_rrna_len_archaea,
        min_rrna_len_bacteria=min_rrna_len_bacteria,
        flank=flank,
        allow_partial=allow_partial,
        reject_non_atgc=reject_non_atgc,
    )
    console.print(
        "[green]Prepared dataset:[/green] "
        f"{summary.dataset} prefix={summary.prefix} "
        f"final_prepared_sequences={summary.final_prepared_sequences} "
        f"dir={summary.dataset_dir}"
    )
    console.print(
        "[cyan]Next:[/cyan] "
        f"autotax2 orient-sina --build {build} --dataset {summary.dataset}"
    )


@app.command("orient-sina")
def orient_sina_command(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    dataset: str = typer.Option(..., "--dataset", help="Prepared dataset name."),
    threads: int = typer.Option(4, "--threads", help="Threads for SINA."),
    sina_bin: str = typer.Option("sina", "--sina-bin", help="SINA executable."),
    reference: Optional[Path] = typer.Option(
        None,
        "--reference",
        help="Optional SINA reference/PTDB path.",
    ),
    strict_tool_version: bool = typer.Option(
        False,
        "--strict-tool-version",
        help="Fail if SINA version cannot be detected.",
    ),
    allow_sina_failure: bool = typer.Option(
        True,
        "--allow-sina-failure/--no-allow-sina-failure",
        help="Allow fallback behavior when SINA fails.",
    ),
    fallback_copy_original: bool = typer.Option(
        True,
        "--fallback-copy-original/--no-fallback-copy-original",
        help="Copy original barrnap FASTA records when SINA output is missing or failed.",
    ),
    min_sina_identity: float = typer.Option(
        0.0,
        "--min-sina-identity",
        help="Reserved loose identity confidence threshold.",
    ),
    min_sina_score: float = typer.Option(
        0.0,
        "--min-sina-score",
        help="Reserved loose score confidence threshold.",
    ),
) -> None:
    """Orient prepared rRNA sequences with SINA loose settings."""
    summary = orient_dataset_with_sina(
        build=build,
        dataset=dataset,
        threads=threads,
        sina_bin=sina_bin,
        reference=reference,
        strict_tool_version=strict_tool_version,
        allow_sina_failure=allow_sina_failure,
        fallback_copy_original=fallback_copy_original,
        min_sina_identity=min_sina_identity,
        min_sina_score=min_sina_score,
    )
    fallback = " fallback_used=true" if summary.fallback_used else ""
    console.print(
        "[green]Completed SINA orientation:[/green] "
        f"{summary.dataset} records={summary.oriented_records}{fallback}"
    )


@app.command("cluster-search")
def cluster_search_command(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    dataset: str = typer.Option(..., "--dataset", help="Prepared and oriented dataset name."),
    threads: int = typer.Option(4, "--threads", help="Threads for VSEARCH."),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin", help="VSEARCH executable."),
    strict_tool_version: bool = typer.Option(
        False,
        "--strict-tool-version",
        help="Fail if VSEARCH version cannot be detected.",
    ),
    iddef: int = typer.Option(2, "--iddef", help="VSEARCH --iddef value."),
    species_id: float = typer.Option(0.987, "--species-id", help="Species cluster identity."),
    genus_id: float = typer.Option(0.945, "--genus-id", help="Genus cluster identity."),
    family_id: float = typer.Option(0.865, "--family-id", help="Family cluster identity."),
    order_id: float = typer.Option(0.820, "--order-id", help="Order cluster identity."),
    class_id: float = typer.Option(0.785, "--class-id", help="Class cluster identity."),
    floor_id: float = typer.Option(0.750, "--floor-id", help="Registry search identity floor."),
    min_query_cov: float = typer.Option(0.80, "--min-query-cov", help="Minimum query coverage."),
    min_target_cov: float = typer.Option(0.0, "--min-target-cov", help="Minimum target coverage."),
    maxaccepts: int = typer.Option(50, "--maxaccepts", help="VSEARCH --maxaccepts."),
    maxrejects: int = typer.Option(256, "--maxrejects", help="VSEARCH --maxrejects."),
    near_best_delta: float = typer.Option(
        0.005,
        "--near-best-delta",
        help="Near-best hit retention delta for later consensus.",
    ),
    strand: StrandChoice = typer.Option(StrandChoice.PLUS, "--strand", help="Search strand."),
) -> None:
    """Cluster oriented dataset sequences and search current registry representatives."""
    summary = cluster_search_dataset(
        build=build,
        dataset=dataset,
        threads=threads,
        vsearch_bin=vsearch_bin,
        strict_tool_version=strict_tool_version,
        iddef=iddef,
        species_id=species_id,
        genus_id=genus_id,
        family_id=family_id,
        order_id=order_id,
        class_id=class_id,
        floor_id=floor_id,
        min_query_cov=min_query_cov,
        min_target_cov=min_target_cov,
        maxaccepts=maxaccepts,
        maxrejects=maxrejects,
        near_best_delta=near_best_delta,
        strand=strand.value,
    )
    console.print(
        "[green]Completed VSEARCH cluster/search:[/green] "
        f"{summary.dataset} filtered_hits={summary.registry_hits_filtered}"
    )


@app.command("place")
def place_command(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    dataset: str = typer.Option(..., "--dataset", help="Cluster-searched dataset name."),
    near_best_delta: float = typer.Option(
        0.005,
        "--near-best-delta",
        help="Identity delta for retaining near-best hits.",
    ),
    rank_consensus: float = typer.Option(
        0.80,
        "--rank-consensus",
        help="Minimum near-best agreement fraction for stable rank consensus.",
    ),
    species_id: float = typer.Option(0.987, "--species-id", help="Known-like species identity."),
    genus_id: float = typer.Option(0.945, "--genus-id", help="New species identity boundary."),
    family_id: float = typer.Option(0.865, "--family-id", help="New genus identity boundary."),
    order_id: float = typer.Option(0.820, "--order-id", help="New family identity boundary."),
    class_id: float = typer.Option(0.785, "--class-id", help="New order identity boundary."),
    floor_id: float = typer.Option(0.750, "--floor-id", help="Minimum placement search identity."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Write proposed placements only."),
    allow_ambiguous: bool = typer.Option(
        True,
        "--allow-ambiguous/--no-allow-ambiguous",
        help="Allow ambiguous records to be written instead of failing.",
    ),
) -> None:
    """Place dataset representatives into the current rank-aware registry."""
    summary = place_dataset(
        build=build,
        dataset=dataset,
        near_best_delta=near_best_delta,
        rank_consensus=rank_consensus,
        species_id=species_id,
        genus_id=genus_id,
        family_id=family_id,
        order_id=order_id,
        class_id=class_id,
        floor_id=floor_id,
        dry_run=dry_run,
        allow_ambiguous=allow_ambiguous,
    )
    mode = "dry-run " if summary.dry_run else ""
    console.print(
        f"[green]Completed {mode}placement:[/green] "
        f"{summary.dataset} assignments={summary.assignments} "
        f"created_taxa={summary.created_taxa}"
    )


@app.command()
def add(
    dataset_prefix: str = typer.Argument(
        ...,
        help="Dataset prefix to freeze for this addition, for example D20.",
    ),
    config: Path = typer.Option(
        Path("autotax2.yaml"),
        "--config",
        "-c",
        help="Project configuration file.",
    ),
    fasta: Optional[Path] = typer.Option(None, "--fasta", help="Input FASTA file."),
) -> None:
    """Add a dataset to the incremental reference scaffold."""
    load_config(config)
    console.print(
        "[yellow]add is a placeholder command.[/yellow] "
        f"dataset_prefix={dataset_prefix} fasta={fasta}"
    )


@app.command("export")
def export_command(
    format_name: str = typer.Argument(
        "all",
        help="Export format: all, sintax, qiime2, or dada2.",
    ),
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    outdir: Optional[Path] = typer.Option(
        None,
        "--outdir",
        help="Export output directory. Defaults to <build>/export.",
    ),
    gzip_output: bool = typer.Option(
        True,
        "--gzip/--no-gzip",
        help="Gzip FASTA outputs.",
    ),
    representatives_only: bool = typer.Option(
        True,
        "--representatives-only/--all-unique",
        help="Export active representatives only, or all active unique MD5 sequences.",
    ),
    prefix: str = typer.Option(
        "autotax2",
        "--prefix",
        help="Output filename prefix.",
    ),
    force: bool = typer.Option(False, "--force", help="Overwrite existing export files."),
) -> None:
    """Export reference files for downstream classifiers."""
    summary = export_references(
        build=build,
        format_name=format_name,
        outdir=outdir,
        gzip_output=gzip_output,
        representatives_only=representatives_only,
        output_prefix=prefix,
        force=force,
    )
    console.print(
        "[green]Completed export:[/green] "
        f"formats={','.join(summary.formats)} records={summary.records_exported} "
        f"outdir={summary.outdir}"
    )


@app.command()
def summarize(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    outdir: Optional[Path] = typer.Option(
        None,
        "--outdir",
        help="Report output directory. Defaults to <build>/reports.",
    ),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing report files."),
) -> None:
    """Summarize the current reference build state."""
    summary = summarize_build(build=build, outdir=outdir, overwrite=overwrite)
    console.print(
        "[green]Wrote reports:[/green] "
        f"files={summary.files_written} outdir={summary.outdir}"
    )


@app.command()
def validate(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    strict: bool = typer.Option(False, "--strict", help="Treat selected warnings as failures."),
    check_exports: bool = typer.Option(
        True,
        "--check-exports/--no-check-exports",
        help="Validate existing export files when present.",
    ),
    report: Optional[Path] = typer.Option(
        None,
        "--report",
        help="Validation Markdown report path. Defaults to <build>/reports/validation_report.md.",
    ),
) -> None:
    """Validate configuration and registry invariants."""
    summary = validate_build(
        build=build,
        strict=strict,
        check_exports=check_exports,
        report=report,
    )
    level = "red" if summary.failed else "green"
    console.print(
        f"[{level}]Validation complete:[/{level}] "
        f"errors={summary.errors} warnings={summary.warnings} report={summary.report_md}"
    )
    if summary.failed:
        raise typer.Exit(1)


def main() -> None:
    """Console-script entry point."""
    app()


if __name__ == "__main__":
    main()
