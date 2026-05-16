"""Command-line interface for autotax2."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from autotax2.audit import write_event_log
from autotax2.dataset_prepare import prepare_dataset
from autotax2.export import export_references
from autotax2.placement import place_dataset
from autotax2.reports import summarize_build
from autotax2.sina import orient_dataset_with_sina
from autotax2.silva import initialize_silva_build, resolve_silva_unresolved_build
from autotax2.thresholds import (
    DEFAULT_CLASS_ID,
    DEFAULT_FAMILY_ID,
    DEFAULT_GENUS_ID,
    DEFAULT_ORDER_ID,
    DEFAULT_PHYLUM_ID,
    DEFAULT_SPECIES_ID,
)
from autotax2.validate import validate_build
from autotax2.vsearch import cluster_search_dataset

HELP_CONTEXT = {"help_option_names": ["-h", "--help"]}

app = typer.Typer(
    name="autotax2",
    help="Fixed-backbone, rank-aware, incremental rRNA gene reference builder.",
    no_args_is_help=True,
    add_completion=False,
    context_settings=HELP_CONTEXT,
)
console = Console()


class DomainChoice(str, Enum):
    """Supported dataset domains."""

    ARCHAEA = "Archaea"
    BACTERIA = "Bacteria"


class StrandChoice(str, Enum):
    """Supported VSEARCH strand modes."""

    PLUS = "plus"
    BOTH = "both"


@app.command(context_settings=HELP_CONTEXT)
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
    type_strain_metadata: Path = typer.Option(
        ...,
        "--type-strain-metadata",
        help="Required SILVA official full_metadata TSV/TSV.gz.",
    ),
    gtdb_ar53_taxonomy: Path = typer.Option(
        ...,
        "--gtdb-ar53-taxonomy",
        help="Required GTDB r232 ar53 taxonomy TSV/TSV.gz.",
    ),
    gtdb_bac120_taxonomy: Path = typer.Option(
        ...,
        "--gtdb-bac120-taxonomy",
        help="Required GTDB r232 bac120 taxonomy TSV/TSV.gz.",
    ),
    threads: int = typer.Option(
        1,
        "--threads",
        help="Threads reserved for initialization steps and recorded in build metadata.",
    ),
) -> None:
    """Initialize a build directory from a SILVA FASTA backbone."""
    start_log = write_event_log(
        outdir,
        "init_start",
        {
            "status": "started",
            "silva_fasta": silva_fasta,
            "type_strain_metadata": type_strain_metadata,
            "gtdb_ar53_taxonomy": gtdb_ar53_taxonomy,
            "gtdb_bac120_taxonomy": gtdb_bac120_taxonomy,
            "threads": threads,
        },
    )
    console.print(
        "[cyan]Starting SILVA initialization.[/cyan] "
        "Full SILVA gzipped inputs can take several minutes; "
        f"start_log={start_log}"
    )
    summary = initialize_silva_build(
        silva_fasta=silva_fasta,
        outdir=outdir,
        type_strain_metadata=type_strain_metadata,
        gtdb_taxonomies=(gtdb_ar53_taxonomy, gtdb_bac120_taxonomy),
        threads=threads,
    )
    audit_log = write_event_log(
        outdir,
        "init",
        {
            "status": "completed",
            "silva_fasta": silva_fasta,
            "type_strain_metadata": type_strain_metadata,
            "gtdb_ar53_taxonomy": gtdb_ar53_taxonomy,
            "gtdb_bac120_taxonomy": gtdb_bac120_taxonomy,
            "records": summary["records"],
            "named": summary["named"],
            "unresolved": summary["unresolved"],
            "rejected": summary.get("rejected", 0),
            "threads": threads,
        },
    )
    console.print(
        "[green]Initialized autotax2 build:[/green] "
        f"{outdir} records={summary['records']} "
        f"named={summary['named']} unresolved={summary['unresolved']} "
        f"rejected={summary.get('rejected', 0)} audit_log={audit_log}"
    )


@app.command("resolve", context_settings=HELP_CONTEXT)
def resolve_silva(
    build: Path = typer.Option(
        ...,
        "--build",
        help="Initialized autotax2 build directory.",
    ),
    threads: int = typer.Option(4, "--threads", help="Threads for VSEARCH clustering."),
    species_id: float = typer.Option(DEFAULT_SPECIES_ID, "--species-id", help="Species-level parent-scoped resolve identity."),
    genus_id: float = typer.Option(DEFAULT_GENUS_ID, "--genus-id", help="Genus-level parent-scoped resolve identity."),
    family_id: float = typer.Option(
        DEFAULT_FAMILY_ID,
        "--family-id",
        help="Family-level parent-scoped resolve identity.",
    ),
    order_id: float = typer.Option(
        DEFAULT_ORDER_ID,
        "--order-id",
        help="Order-level parent-scoped resolve identity.",
    ),
    class_id: float = typer.Option(
        DEFAULT_CLASS_ID,
        "--class-id",
        help="Class-level parent-scoped resolve identity.",
    ),
    phylum_id: float = typer.Option(
        DEFAULT_PHYLUM_ID,
        "--phylum-id",
        help="Phylum-level parent-scoped resolve identity.",
    ),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin", help="VSEARCH executable."),
    iddef: int = typer.Option(2, "--iddef", help="VSEARCH --iddef value."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Write proposals without registry updates."),
) -> None:
    """Resolve SILVA unresolved records into a placeholder framework."""
    summary = resolve_silva_unresolved_build(
        build=build,
        threads=threads,
        species_id=species_id,
        genus_id=genus_id,
        family_id=family_id,
        order_id=order_id,
        class_id=class_id,
        phylum_id=phylum_id,
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
    write_event_log(
        build,
        "resolve",
        {
            "threads": threads,
            "species_id": species_id,
            "genus_id": genus_id,
            "dry_run": summary.dry_run,
            "unresolved_records": summary.unresolved_records,
            "placeholder_taxa": summary.placeholder_taxa,
        },
    )


@app.command("prepare", context_settings=HELP_CONTEXT)
def prepare_dataset_command(
    build: Path = typer.Option(..., "--build", help="Initialized autotax2 build directory."),
    name: str = typer.Option(..., "--name", help="Dataset name."),
    prefix: str = typer.Option(..., "--prefix", help="Frozen dataset prefix."),
    fasta: Path = typer.Option(
        ...,
        "--fasta",
        help="Input FASTA already processed to SSU/16S sequences outside autotax2.",
    ),
    domain: DomainChoice = typer.Option(..., "--domain", help="Dataset domain."),
    min_ssu_len_archaea: int = typer.Option(
        900,
        "--min-ssu-len-archaea",
        help="Minimum prepared archaeal SSU/16S length.",
    ),
    min_ssu_len_bacteria: int = typer.Option(
        1200,
        "--min-ssu-len-bacteria",
        help="Minimum prepared bacterial SSU/16S length.",
    ),
    reject_non_atgc: bool = typer.Option(
        True,
        "--reject-non-atgc/--no-reject-non-atgc",
        help="Reject normalized input sequences containing non-ATGC symbols.",
    ),
) -> None:
    """Prepare a custom dataset from externally processed SSU/16S sequences."""
    summary = prepare_dataset(
        build=build,
        name=name,
        prefix=prefix,
        fasta=fasta,
        domain=domain.value,
        min_ssu_len_archaea=min_ssu_len_archaea,
        min_ssu_len_bacteria=min_ssu_len_bacteria,
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
        f"autotax2 orient --build {build} --dataset {summary.dataset}"
    )
    write_event_log(
        build,
        "prepare",
        {
            "dataset": summary.dataset,
            "prefix": summary.prefix,
            "input_fasta": fasta,
            "domain": domain.value,
            "final_prepared_sequences": summary.final_prepared_sequences,
            "dataset_dir": summary.dataset_dir,
        },
    )


@app.command("orient", context_settings=HELP_CONTEXT)
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
        help="Copy prepared SSU/16S FASTA records when SINA output is missing or failed.",
    ),
    min_sina_identity: float = typer.Option(
        0.0,
        "--min-sina-identity",
        help="Currently unsupported SINA metric filter; leave at 0.0.",
    ),
    min_sina_score: float = typer.Option(
        0.0,
        "--min-sina-score",
        help="Currently unsupported SINA metric filter; leave at 0.0.",
    ),
    search_candidates: bool = typer.Option(
        False,
        "--search-candidates/--no-search-candidates",
        help="Ask SINA to write nearest_slv candidate hits for later VSEARCH rescoring.",
    ),
    search_db: Optional[Path] = typer.Option(
        None,
        "--search-db",
        help="Optional SINA search database path.",
    ),
    search_min_sim: float = typer.Option(
        0.5,
        "--search-min-sim",
        help="Loose SINA candidate search similarity; VSEARCH still performs final threshold scoring.",
    ),
    search_max_result: int = typer.Option(
        10,
        "--search-max-result",
        help="Maximum SINA nearest_slv candidates retained per query.",
    ),
    search_kmer_candidates: int = typer.Option(
        1000,
        "--search-kmer-candidates",
        help="SINA k-mer candidate pool size for candidate search.",
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
        search_candidates=search_candidates,
        search_db=search_db,
        search_min_sim=search_min_sim,
        search_max_result=search_max_result,
        search_kmer_candidates=search_kmer_candidates,
    )
    fallback = " fallback_used=true" if summary.fallback_used else ""
    console.print(
        "[green]Completed SINA orientation:[/green] "
        f"{summary.dataset} records={summary.oriented_records}{fallback}"
    )
    write_event_log(
        build,
        "orient",
        {
            "dataset": summary.dataset,
            "threads": threads,
            "records": summary.oriented_records,
            "fallback_used": summary.fallback_used,
            "search_candidates": search_candidates,
            "dataset_dir": summary.dataset_dir,
        },
    )


@app.command("cluster", context_settings=HELP_CONTEXT)
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
    species_id: float = typer.Option(DEFAULT_SPECIES_ID, "--species-id", help="Species cluster identity."),
    genus_id: float = typer.Option(DEFAULT_GENUS_ID, "--genus-id", help="Genus cluster identity."),
    family_id: float = typer.Option(DEFAULT_FAMILY_ID, "--family-id", help="Family cluster identity."),
    order_id: float = typer.Option(DEFAULT_ORDER_ID, "--order-id", help="Order cluster identity."),
    class_id: float = typer.Option(DEFAULT_CLASS_ID, "--class-id", help="Class cluster identity."),
    phylum_id: float = typer.Option(DEFAULT_PHYLUM_ID, "--phylum-id", help="Phylum cluster and registry search identity floor."),
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
    sina_candidates: Optional[Path] = typer.Option(
        None,
        "--sina-candidates",
        help="Optional SINA candidate TSV; defaults to dataset/sina.candidates.tsv when present.",
    ),
    require_sina_candidates: bool = typer.Option(
        False,
        "--require-sina-candidates/--no-require-sina-candidates",
        help="Fail if the SINA candidate TSV has no target matching current representatives.",
    ),
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
        phylum_id=phylum_id,
        floor_id=phylum_id,
        min_query_cov=min_query_cov,
        min_target_cov=min_target_cov,
        maxaccepts=maxaccepts,
        maxrejects=maxrejects,
        near_best_delta=near_best_delta,
        strand=strand.value,
        sina_candidates=sina_candidates,
        require_sina_candidates=require_sina_candidates,
    )
    console.print(
        "[green]Completed VSEARCH cluster/search:[/green] "
        f"{summary.dataset} filtered_hits={summary.registry_hits_filtered}"
    )
    write_event_log(
        build,
        "cluster",
        {
            "dataset": summary.dataset,
            "threads": threads,
            "iddef": iddef,
            "sina_candidates": sina_candidates,
            "require_sina_candidates": require_sina_candidates,
            "registry_hits_filtered": summary.registry_hits_filtered,
        },
    )


@app.command("place", context_settings=HELP_CONTEXT)
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
    species_id: float = typer.Option(DEFAULT_SPECIES_ID, "--species-id", help="Known-like species identity."),
    genus_id: float = typer.Option(DEFAULT_GENUS_ID, "--genus-id", help="New species identity boundary."),
    family_id: float = typer.Option(DEFAULT_FAMILY_ID, "--family-id", help="New genus identity boundary."),
    order_id: float = typer.Option(DEFAULT_ORDER_ID, "--order-id", help="New family identity boundary."),
    class_id: float = typer.Option(DEFAULT_CLASS_ID, "--class-id", help="New order identity boundary."),
    phylum_id: float = typer.Option(DEFAULT_PHYLUM_ID, "--phylum-id", help="Minimum placement identity for new class calls."),
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
        phylum_id=phylum_id,
        dry_run=dry_run,
        allow_ambiguous=allow_ambiguous,
    )
    mode = "dry-run " if summary.dry_run else ""
    console.print(
        f"[green]Completed {mode}placement:[/green] "
        f"{summary.dataset} assignments={summary.assignments} "
        f"created_taxa={summary.created_taxa}"
    )
    write_event_log(
        build,
        "place",
        {
            "dataset": summary.dataset,
            "assignments": summary.assignments,
            "created_taxa": summary.created_taxa,
            "dry_run": summary.dry_run,
        },
    )


@app.command(context_settings=HELP_CONTEXT)
def add(
    build: Path = typer.Option(..., "--build", help="Initialized and resolved autotax2 build directory."),
    name: str = typer.Option(..., "--name", help="Dataset name."),
    prefix: str = typer.Option(..., "--prefix", help="Frozen dataset prefix."),
    fasta: Path = typer.Option(
        ...,
        "--fasta",
        help="Input FASTA already processed to SSU/16S sequences outside autotax2.",
    ),
    domain: DomainChoice = typer.Option(..., "--domain", help="Dataset domain."),
    threads: int = typer.Option(4, "--threads", help="Threads for SINA and VSEARCH steps."),
    min_ssu_len_archaea: int = typer.Option(
        900,
        "--min-ssu-len-archaea",
        help="Minimum prepared archaeal SSU/16S length.",
    ),
    min_ssu_len_bacteria: int = typer.Option(
        1200,
        "--min-ssu-len-bacteria",
        help="Minimum prepared bacterial SSU/16S length.",
    ),
    reject_non_atgc: bool = typer.Option(
        True,
        "--reject-non-atgc/--no-reject-non-atgc",
        help="Reject normalized input sequences containing non-ATGC symbols.",
    ),
    sina_bin: str = typer.Option("sina", "--sina-bin", help="SINA executable."),
    sina_reference: Optional[Path] = typer.Option(
        None,
        "--sina-reference",
        help="Optional SINA reference/PTDB path.",
    ),
    allow_sina_failure: bool = typer.Option(
        True,
        "--allow-sina-failure/--no-allow-sina-failure",
        help="Allow fallback behavior when SINA fails.",
    ),
    fallback_copy_original: bool = typer.Option(
        True,
        "--fallback-copy-original/--no-fallback-copy-original",
        help="Copy prepared SSU/16S FASTA records when SINA output is missing or failed.",
    ),
    search_candidates: bool = typer.Option(
        True,
        "--search-candidates/--no-search-candidates",
        help="Ask SINA to write nearest_slv candidate hits for later VSEARCH rescoring.",
    ),
    search_db: Optional[Path] = typer.Option(
        None,
        "--search-db",
        help="Optional SINA search database path.",
    ),
    search_min_sim: float = typer.Option(
        0.5,
        "--search-min-sim",
        help="Loose SINA candidate search similarity; VSEARCH still performs final threshold scoring.",
    ),
    search_max_result: int = typer.Option(
        10,
        "--search-max-result",
        help="Maximum SINA nearest_slv candidates retained per query.",
    ),
    search_kmer_candidates: int = typer.Option(
        1000,
        "--search-kmer-candidates",
        help="SINA k-mer candidate pool size for candidate search.",
    ),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin", help="VSEARCH executable."),
    strict_tool_version: bool = typer.Option(
        False,
        "--strict-tool-version",
        help="Fail if SINA or VSEARCH version cannot be detected.",
    ),
    iddef: int = typer.Option(2, "--iddef", help="VSEARCH --iddef value."),
    species_id: float = typer.Option(DEFAULT_SPECIES_ID, "--species-id", help="Species identity threshold."),
    genus_id: float = typer.Option(DEFAULT_GENUS_ID, "--genus-id", help="Genus identity threshold."),
    family_id: float = typer.Option(DEFAULT_FAMILY_ID, "--family-id", help="Family identity threshold."),
    order_id: float = typer.Option(DEFAULT_ORDER_ID, "--order-id", help="Order identity threshold."),
    class_id: float = typer.Option(DEFAULT_CLASS_ID, "--class-id", help="Class identity threshold."),
    phylum_id: float = typer.Option(DEFAULT_PHYLUM_ID, "--phylum-id", help="Phylum identity threshold and search floor."),
    min_query_cov: float = typer.Option(0.80, "--min-query-cov", help="Minimum VSEARCH query coverage."),
    min_target_cov: float = typer.Option(0.0, "--min-target-cov", help="Minimum VSEARCH target coverage."),
    maxaccepts: int = typer.Option(50, "--maxaccepts", help="VSEARCH --maxaccepts."),
    maxrejects: int = typer.Option(256, "--maxrejects", help="VSEARCH --maxrejects."),
    near_best_delta: float = typer.Option(
        0.005,
        "--near-best-delta",
        help="Near-best hit retention delta for placement consensus.",
    ),
    rank_consensus: float = typer.Option(
        0.80,
        "--rank-consensus",
        help="Minimum near-best agreement fraction for stable rank consensus.",
    ),
    strand: StrandChoice = typer.Option(StrandChoice.PLUS, "--strand", help="VSEARCH search strand."),
    require_sina_candidates: bool = typer.Option(
        False,
        "--require-sina-candidates/--no-require-sina-candidates",
        help="Fail if the SINA candidate TSV has no target matching current representatives.",
    ),
    allow_ambiguous: bool = typer.Option(
        True,
        "--allow-ambiguous/--no-allow-ambiguous",
        help="Allow ambiguous records to be written instead of failing.",
    ),
    export_after: bool = typer.Option(
        False,
        "--export/--no-export",
        help="Export classifier references after placement.",
    ),
    gzip_output: bool = typer.Option(
        True,
        "--gzip/--no-gzip",
        help="Gzip FASTA outputs when --export is used.",
    ),
    force_export: bool = typer.Option(
        False,
        "--force-export",
        help="Overwrite existing export files when --export is used.",
    ),
    validate_after: bool = typer.Option(
        True,
        "--validate/--no-validate",
        help="Run validation after summarize/export.",
    ),
    strict_validate: bool = typer.Option(
        False,
        "--strict-validate",
        help="Treat selected validation warnings as failures.",
    ),
) -> None:
    """Run prepare, orient, cluster, place, summarize, and optional export for one dataset."""
    start_log = write_event_log(
        build,
        "add_start",
        {
            "status": "started",
            "dataset": name,
            "prefix": prefix,
            "fasta": fasta,
            "domain": domain.value,
            "threads": threads,
            "search_candidates": search_candidates,
            "export_after": export_after,
            "validate_after": validate_after,
        },
    )
    console.print(
        "[cyan]Starting dataset add workflow.[/cyan] "
        f"dataset={name} prefix={prefix} start_log={start_log}"
    )
    prepared = prepare_dataset(
        build=build,
        name=name,
        prefix=prefix,
        fasta=fasta,
        domain=domain.value,
        min_ssu_len_archaea=min_ssu_len_archaea,
        min_ssu_len_bacteria=min_ssu_len_bacteria,
        reject_non_atgc=reject_non_atgc,
    )
    oriented = orient_dataset_with_sina(
        build=build,
        dataset=name,
        threads=threads,
        sina_bin=sina_bin,
        reference=sina_reference,
        strict_tool_version=strict_tool_version,
        allow_sina_failure=allow_sina_failure,
        fallback_copy_original=fallback_copy_original,
        search_candidates=search_candidates,
        search_db=search_db,
        search_min_sim=search_min_sim,
        search_max_result=search_max_result,
        search_kmer_candidates=search_kmer_candidates,
    )
    clustered = cluster_search_dataset(
        build=build,
        dataset=name,
        threads=threads,
        vsearch_bin=vsearch_bin,
        strict_tool_version=strict_tool_version,
        iddef=iddef,
        species_id=species_id,
        genus_id=genus_id,
        family_id=family_id,
        order_id=order_id,
        class_id=class_id,
        phylum_id=phylum_id,
        floor_id=phylum_id,
        min_query_cov=min_query_cov,
        min_target_cov=min_target_cov,
        maxaccepts=maxaccepts,
        maxrejects=maxrejects,
        near_best_delta=near_best_delta,
        strand=strand.value,
        require_sina_candidates=require_sina_candidates,
    )
    placed = place_dataset(
        build=build,
        dataset=name,
        near_best_delta=near_best_delta,
        rank_consensus=rank_consensus,
        species_id=species_id,
        genus_id=genus_id,
        family_id=family_id,
        order_id=order_id,
        class_id=class_id,
        phylum_id=phylum_id,
        allow_ambiguous=allow_ambiguous,
    )
    reports = summarize_build(build=build, overwrite=True)
    export_summary = None
    if export_after:
        export_summary = export_references(
            build=build,
            format_name="all",
            gzip_output=gzip_output,
            force=force_export,
        )
    validation = None
    if validate_after:
        validation = validate_build(
            build=build,
            strict=strict_validate,
            check_exports=export_after,
        )
    audit_log = write_event_log(
        build,
        "add",
        {
            "status": "completed",
            "dataset": name,
            "prefix": prefix,
            "dataset_dir": prepared.dataset_dir,
            "prepared_sequences": prepared.final_prepared_sequences,
            "oriented_records": oriented.oriented_records,
            "sina_fallback_used": oriented.fallback_used,
            "registry_hits_filtered": clustered.registry_hits_filtered,
            "assignments": placed.assignments,
            "created_taxa": placed.created_taxa,
            "reports_written": reports.files_written,
            "export_records": export_summary.records_exported if export_summary else 0,
            "validation_errors": validation.errors if validation else "",
            "validation_warnings": validation.warnings if validation else "",
            "strict_validate": strict_validate,
        },
    )
    console.print(
        "[green]Completed dataset add workflow:[/green] "
        f"{name} prepared={prepared.final_prepared_sequences} "
        f"oriented={oriented.oriented_records} "
        f"filtered_hits={clustered.registry_hits_filtered} "
        f"assignments={placed.assignments} "
        f"created_taxa={placed.created_taxa} "
        f"reports={reports.outdir} audit_log={audit_log}"
    )
    if export_summary is not None:
        console.print(
            "[green]Exported references:[/green] "
            f"formats={','.join(export_summary.formats)} records={export_summary.records_exported}"
        )
    if validation is not None:
        level = "red" if validation.failed else "green"
        console.print(
            f"[{level}]Validation after add:[/{level}] "
            f"errors={validation.errors} warnings={validation.warnings} report={validation.report_md}"
        )
        if validation.failed:
            raise typer.Exit(1)


@app.command("export", context_settings=HELP_CONTEXT)
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
        f"outdir={summary.outdir} validation={summary.validation_path}"
    )
    write_event_log(
        build,
        "export",
        {
            "formats": ",".join(summary.formats),
            "records_exported": summary.records_exported,
            "outdir": summary.outdir,
            "manifest": summary.manifest_path,
            "validation": summary.validation_path,
            "validation_errors": summary.validation_errors,
        },
    )


@app.command(context_settings=HELP_CONTEXT)
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
    write_event_log(
        build,
        "summarize",
        {
            "files_written": summary.files_written,
            "outdir": summary.outdir,
        },
    )


@app.command(context_settings=HELP_CONTEXT)
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
    write_event_log(
        build,
        "validate",
        {
            "errors": summary.errors,
            "warnings": summary.warnings,
            "report_md": summary.report_md,
            "report_tsv": summary.report_tsv,
            "strict": strict,
            "check_exports": check_exports,
        },
    )
    if summary.failed:
        raise typer.Exit(1)


def main() -> None:
    """Console-script entry point."""
    app()


if __name__ == "__main__":
    main()
