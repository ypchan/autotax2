from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


REQUIRED_SOFTWARE = [
    "vsearch",
    "blastn",
    "makeblastdb",
]

OPTIONAL_SOFTWARE = [
    "sina",
    "mafft",
    "meme",
    "streme",
    "fimo",
    "rnafold",
    "rnaalifold",
    "barrnap",
]

SOFTWARE_EXECUTABLE_NAMES = {
    "vsearch": "vsearch",
    "blastn": "blastn",
    "makeblastdb": "makeblastdb",
    "sina": "sina",
    "mafft": "mafft",
    "meme": "meme",
    "streme": "streme",
    "fimo": "fimo",
    "rnafold": "RNAfold",
    "rnaalifold": "RNAalifold",
    "barrnap": "barrnap",
}

PYTHON_PACKAGES = [
    "typer",
]

CORE_REQUIRED_MANIFEST_KEYS = [
    "silva_fasta",
    "silva_sintax_fasta",
    "silva_typestrains_fasta",
]

CORE_OPTIONAL_MANIFEST_PATH_KEYS = [
    "silva_metadata",
    "typestrains_accession_ids",
    "silva_taxonomy_tsv",
    "silva_udb",
    "silva_sintax_udb",
    "silva_typestrains_udb",
]

INTRON_OPTIONAL_MANIFEST_PATH_KEYS = [
    "species_rep1_fasta",
    "species_rep1_taxonomy",
    "species_rep1_clean_fasta",
    "species_rep1_clean_taxonomy",
    "species_rep10_fasta",
    "species_rep10_taxonomy",
    "species_rep10_raw_fasta",
    "species_rep10_raw_taxonomy",
    "species_rep10_raw_selection",
    "reference_local_hsps",
    "reference_intron_candidates",
    "reference_intron_audit",
    "reference_intron_blacklist",
    "reference_barrnap_gff",
    "reference_barrnap_audit",
    "coordinate_ref",
    "sina_ref",
]

INTRON_OPTIONAL_MANIFEST_BLASTDB_KEYS = [
    "species_rep1_blastdb",
    "species_rep1_clean_blastdb",
    "species_rep10_blastdb",
    "species_rep10_raw_blastdb",
]


@dataclass(frozen=True)
class CheckResult:
    section: str
    key: str
    ok: bool
    value: str = ""
    required: bool = True
    message: str = ""


def normalized_key(key: str) -> str:
    return key.strip().replace("-", "_").lower()


def status_label(result: CheckResult) -> str:
    if result.ok:
        return "OK"
    if result.required:
        return "ERROR"
    return "SKIP"


def status_style(result: CheckResult) -> str:
    if result.ok:
        return "bold green"
    if result.required:
        return "bold red"
    return "yellow"


def required_label(required: bool) -> str:
    return "required" if required else "optional"


def expand_path(value: str | Path | None) -> Optional[Path]:
    if value is None:
        return None

    text = str(value).strip()
    if not text:
        return None

    return Path(text).expanduser()


def executable_exists(value: str | Path | None) -> bool:
    if value is None:
        return False

    text = str(value).strip()
    if not text:
        return False

    path = Path(text).expanduser()

    if path.exists() and os.access(path, os.X_OK):
        return True

    return shutil.which(text) is not None


def resolve_executable(value: str | Path | None) -> Optional[str]:
    if value is None:
        return None

    text = str(value).strip()
    if not text:
        return None

    path = Path(text).expanduser()

    if path.exists() and os.access(path, os.X_OK):
        return str(path)

    detected = shutil.which(text)
    if detected:
        return detected

    return None


def executable_name_for_key(key: str) -> str:
    return SOFTWARE_EXECUTABLE_NAMES.get(key.lower(), key)


def executable_version(executable: str | Path) -> str:
    resolved = resolve_executable(executable)

    if not resolved:
        return ""

    version_commands = [
        [resolved, "--version"],
        [resolved, "-version"],
        [resolved, "-h"],
    ]

    for command in version_commands:
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue

        output = (completed.stdout or completed.stderr or "").strip()
        if output:
            return output.splitlines()[0].strip()

    return ""


def check_software(
    key: str,
    configured_path: str | Path | None = None,
    *,
    required: bool,
) -> CheckResult:
    normalized = normalized_key(key)
    executable_name = executable_name_for_key(normalized)

    candidates: List[str | Path] = []

    if configured_path:
        candidates.append(configured_path)

    candidates.append(executable_name)

    resolved: Optional[str] = None

    for candidate in candidates:
        resolved = resolve_executable(candidate)
        if resolved:
            break

    if resolved:
        version = executable_version(resolved)
        return CheckResult(
            section="Software",
            key=normalized,
            ok=True,
            value=resolved,
            required=required,
            message=version if version else "found",
        )

    return CheckResult(
        section="Software",
        key=normalized,
        ok=False,
        value=str(configured_path or ""),
        required=required,
        message=f"missing executable: {executable_name}",
    )


def python_package_available(package: str) -> bool:
    return importlib.util.find_spec(package) is not None


def check_python_package(package: str) -> CheckResult:
    ok = python_package_available(package)

    return CheckResult(
        section="Python packages",
        key=package,
        ok=ok,
        value="importable" if ok else "",
        required=True,
        message="importable" if ok else "not importable",
    )


def check_python_packages(
    packages: Sequence[str] = PYTHON_PACKAGES,
) -> List[CheckResult]:
    return [check_python_package(package) for package in packages]


def path_exists(path: str | Path | None) -> bool:
    expanded = expand_path(path)
    return bool(expanded and expanded.exists())


def check_path(
    key: str,
    path: str | Path | None,
    *,
    required: bool,
    section: str = "Reference data",
) -> CheckResult:
    expanded = expand_path(path)

    if expanded is None:
        return CheckResult(
            section=section,
            key=key,
            ok=False,
            value="",
            required=required,
            message="not configured",
        )

    ok = expanded.exists()

    return CheckResult(
        section=section,
        key=key,
        ok=ok,
        value=str(expanded),
        required=required,
        message="found" if ok else "missing file or directory",
    )


def blastdb_exists(prefix: str | Path | None) -> bool:
    expanded = expand_path(prefix)

    if expanded is None:
        return False

    candidates = [
        Path(str(expanded) + ".nhr"),
        Path(str(expanded) + ".nin"),
        Path(str(expanded) + ".nsq"),
        Path(str(expanded) + ".ndb"),
        Path(str(expanded) + ".not"),
        Path(str(expanded) + ".ntf"),
        Path(str(expanded) + ".nto"),
        Path(str(expanded) + ".00.nhr"),
        Path(str(expanded) + ".00.nin"),
        Path(str(expanded) + ".00.nsq"),
        Path(str(expanded) + ".00.ndb"),
        Path(str(expanded) + ".00.not"),
        Path(str(expanded) + ".00.ntf"),
        Path(str(expanded) + ".00.nto"),
    ]

    return any(path.exists() for path in candidates)


def check_blastdb(
    key: str,
    prefix: str | Path | None,
    *,
    required: bool = False,
    section: str = "Reference data",
) -> CheckResult:
    expanded = expand_path(prefix)

    if expanded is None:
        return CheckResult(
            section=section,
            key=key,
            ok=False,
            value="",
            required=required,
            message="not configured",
        )

    ok = blastdb_exists(expanded)

    return CheckResult(
        section=section,
        key=key,
        ok=ok,
        value=str(expanded),
        required=required,
        message="found" if ok else "missing BLAST database files",
    )


def read_manifest(manifest_path: str | Path) -> Dict[str, str]:
    manifest = Path(manifest_path).expanduser()

    if not manifest.exists():
        raise FileNotFoundError(f"Reference manifest not found: {manifest}")

    entries: Dict[str, str] = {}

    with manifest.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) < 2:
                continue

            key = parts[0].strip()
            value = parts[1].strip()

            if key:
                entries[key] = value

    return entries


def check_reference_manifest(
    manifest_path: str | Path | None,
) -> List[CheckResult]:
    results: List[CheckResult] = []

    manifest_check = check_path(
        "manifest",
        manifest_path,
        required=True,
        section="Reference manifest",
    )
    results.append(manifest_check)

    if not manifest_check.ok or manifest_path is None:
        return results

    try:
        manifest = read_manifest(manifest_path)
    except OSError as exc:
        results.append(
            CheckResult(
                section="Reference manifest",
                key="manifest_parse",
                ok=False,
                value=str(manifest_path),
                required=True,
                message=str(exc),
            )
        )
        return results

    for key in CORE_REQUIRED_MANIFEST_KEYS:
        results.append(
            check_path(
                key,
                manifest.get(key),
                required=True,
                section="Reference manifest",
            )
        )

    for key in CORE_OPTIONAL_MANIFEST_PATH_KEYS + INTRON_OPTIONAL_MANIFEST_PATH_KEYS:
        results.append(
            check_path(
                key,
                manifest.get(key),
                required=False,
                section="Reference manifest",
            )
        )

    for key in INTRON_OPTIONAL_MANIFEST_BLASTDB_KEYS:
        results.append(
            check_blastdb(
                key,
                manifest.get(key),
                required=False,
                section="Reference manifest",
            )
        )

    return results


def check_from_config(config_path: str | Path | None = None) -> List[CheckResult]:
    from .config import check_config, config_path_message

    results: List[CheckResult] = []

    results.append(
        check_path(
            "config",
            config_path_message(config_path),
            required=False,
            section="Configuration",
        )
    )

    for item in check_config(config_path):
        results.append(
            CheckResult(
                section=item.section,
                key=item.key,
                ok=item.ok,
                value=item.value,
                required=item.required,
                message=item.message,
            )
        )

    results.extend(check_python_packages())

    return results


def required_checks_passed(results: Iterable[CheckResult]) -> bool:
    return all(result.ok for result in results if result.required)


def clean_message(result: CheckResult) -> str:
    message = result.message.strip()

    if result.ok:
        return message or "found"

    if not result.required and (
        not message
        or message.startswith("optional;")
        or message in {"missing", "not configured or missing"}
    ):
        return "not configured"

    if message == "missing":
        return "missing file or directory"

    return message or "not configured"


def purpose_for_key(result: CheckResult) -> str:
    key = normalized_key(result.key)
    section = result.section.lower()

    purposes = {
        "config": "Persistent AutoTax2 settings.",
        "threads": "Default runtime thread count.",
        "vsearch": "Dereplication, clustering, SINTAX, and global searches.",
        "blastn": "Local HSP search for intron detection and coordinate mapping.",
        "makeblastdb": "Build temporary or permanent BLAST nucleotide databases.",
        "sina": "Optional orientation normalization for full-length 16S.",
        "sina_ref": "Reference database used by SINA orientation normalization.",
        "mafft": "Optional intron multiple sequence alignment.",
        "meme": "Optional intron motif discovery with MEME.",
        "streme": "Optional intron motif discovery with STREME.",
        "fimo": "Optional motif scanning against intron sequences.",
        "rnafold": "Optional per-intron secondary-structure prediction.",
        "rnaalifold": "Optional consensus structure prediction from aligned introns.",
        "barrnap": "Optional rRNA annotation audit for reference sequences.",
        "ref_manifest": "Main AutoTax2 reference manifest generated by prepare-silva.",
        "manifest": "Reference manifest supplied with --ref-manifest.",
        "silva_fasta": "Prepared SILVA FASTA reference.",
        "silva_sintax_fasta": "Prepared SILVA SINTAX FASTA reference.",
        "silva_typestrains_fasta": "Prepared type-strain FASTA reference.",
        "silva_metadata": "Original SILVA full_metadata file or copied metadata.",
        "silva_taxonomy_tsv": "Parsed SILVA taxonomy table; needed for build-intron-ref.",
        "typestrains_accession_ids": "Type-strain accession list extracted from SILVA metadata.",
        "silva_udb": "Optional VSEARCH UDB version of SILVA FASTA.",
        "silva_sintax_udb": "Optional VSEARCH UDB version of SINTAX FASTA.",
        "silva_typestrains_udb": "Optional VSEARCH UDB version of type-strain FASTA.",
        "coordinate_ref": "Optional coordinate reference for intron site mapping.",
        "source_map": "Optional sequence-to-source mapping for provenance.",
        "typer": "Python CLI framework used by AutoTax2.",
    }

    if key in purposes:
        return purposes[key]

    if key.startswith("species_rep1_clean"):
        return "Audited clean species reference used by detect-intron."
    if key.startswith("species_rep1"):
        return "Species-level representative reference for intron detection."
    if key.startswith("species_rep10_raw"):
        return "Raw species representatives used for reference intron audit."
    if key.startswith("species_rep10"):
        return "Species-level multi-representative reference for sensitivity checks."
    if key.startswith("reference_intron") or key == "reference_local_hsps":
        return "Reference-intron audit output from build-intron-ref."
    if key.startswith("reference_barrnap"):
        return "barrnap audit output for reference sequences."

    if "python" in section:
        return "Python runtime dependency."

    return "AutoTax2 input, reference, or runtime dependency."


def action_for_key(result: CheckResult) -> str:
    key = normalized_key(result.key)
    section = result.section.lower()

    if result.ok:
        return "ready"

    if key == "config":
        return "Run: autotax2 config init"

    if section == "software" or key in SOFTWARE_EXECUTABLE_NAMES:
        exe = executable_name_for_key(key)
        if key in {"meme", "streme", "fimo"}:
            return f"Optional for motif analysis; install MEME Suite or run: autotax2 config set {key} /path/to/{exe}"
        if key in {"rnafold", "rnaalifold"}:
            return f"Optional for structure analysis; install ViennaRNA or run: autotax2 config set {key} /path/to/{exe}"
        if key == "sina":
            return "Optional for orientation; install SINA or run: autotax2 config set sina /path/to/sina"
        if key == "barrnap":
            return "Optional for reference audit; install barrnap or run: autotax2 config set barrnap /path/to/barrnap"
        if key == "mafft":
            return "Optional for intron alignment; install MAFFT or run: autotax2 config set mafft /path/to/mafft"
        return f"Install {exe}, add it to PATH, or run: autotax2 config set {key} /path/to/{exe}"

    if section == "python packages":
        return f"Install Python package: pip install {result.key}"

    if key in {"ref_manifest", "manifest"}:
        return "Run prepare-silva, or set an existing manifest: autotax2 config set ref-manifest refdatabases/autotax2_ref_manifest.tsv"

    if key in {"silva_fasta", "silva_sintax_fasta", "silva_typestrains_fasta"}:
        return "Run: autotax2 prepare-silva --silva-fasta <SILVA.fasta.gz> --silva-metadata <SILVA.full_metadata.gz> --out refdatabases"

    if key in CORE_OPTIONAL_MANIFEST_PATH_KEYS:
        return "Normally produced by prepare-silva; rerun prepare-silva if this file is needed."

    if key == "sina_ref":
        return "Optional for orient/--orient; run: autotax2 config set sina-ref /path/to/SINA_reference.arb"

    if key == "coordinate_ref":
        return "Optional for map-intron-coordinates; run: autotax2 config set coordinate-ref /path/to/coordinate_reference.fasta"

    if key.startswith("species_rep") or key.startswith("reference_intron") or key.startswith("reference_local_hsps"):
        return "Optional for intron workflow; run: autotax2 build-intron-ref --silva-fasta refdatabases/silva.fasta --silva-taxonomy refdatabases/silva_taxonomy.tsv --out refdatabases/intron_ref"

    if key.startswith("reference_barrnap"):
        return "Optional audit output; rerun build-intron-ref with --run-barrnap if needed."

    if key.endswith("blastdb"):
        return "Run build-intron-ref with --blastdb, or set the BLAST DB prefix in config."

    if key == "source_map":
        return "Optional provenance input; provide --source-map when needed."

    return "Configure this path with autotax2 config set <key> <path>, or generate it with the related workflow command."


def detail_for_result(result: CheckResult) -> str:
    message = clean_message(result)

    if result.value and message:
        return f"{result.value} | {message}"
    if result.value:
        return result.value
    return message


def result_keys(results: Iterable[CheckResult], *, only_missing: bool = True) -> Set[str]:
    keys: Set[str] = set()
    for result in results:
        if only_missing and result.ok:
            continue
        keys.add(normalized_key(result.key))
    return keys


def build_recommended_next_steps(results: Sequence[CheckResult]) -> List[str]:
    missing_all = result_keys(results)
    missing_required = result_keys((r for r in results if r.required))
    missing_optional = result_keys((r for r in results if not r.required))

    steps: List[str] = []

    if "config" in missing_optional:
        steps.append(
            "Create the default config file:\n"
            "  autotax2 config init"
        )

    missing_required_software = missing_required & {"vsearch", "blastn", "makeblastdb"}
    if missing_required_software:
        tools = ", ".join(sorted(missing_required_software))
        steps.append(
            f"Install or configure required software ({tools}). Example:\n"
            "  autotax2 config set vsearch /path/to/vsearch\n"
            "  autotax2 config set blastn /path/to/blastn\n"
            "  autotax2 config set makeblastdb /path/to/makeblastdb"
        )

    core_reference_missing = missing_required & {
        "ref_manifest",
        "manifest",
        "silva_fasta",
        "silva_sintax_fasta",
        "silva_typestrains_fasta",
    }
    if core_reference_missing:
        steps.append(
            "Prepare the core SILVA reference files:\n"
            "  autotax2 prepare-silva \\\n"
            "    --silva-fasta /path/to/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \\\n"
            "    --silva-metadata /path/to/SILVA_138.2_SSURef.full_metadata.gz \\\n"
            "    --out refdatabases"
        )
        steps.append(
            "If the references already exist, point AutoTax2 to the manifest:\n"
            "  autotax2 config set ref-manifest refdatabases/autotax2_ref_manifest.tsv"
        )

    intron_reference_keys = {
        "species_rep1_fasta",
        "species_rep1_taxonomy",
        "species_rep1_clean_fasta",
        "species_rep1_clean_taxonomy",
        "species_rep1_blastdb",
        "species_rep1_clean_blastdb",
        "species_rep10_fasta",
        "species_rep10_taxonomy",
        "species_rep10_raw_fasta",
        "species_rep10_raw_taxonomy",
        "species_rep10_raw_selection",
        "species_rep10_blastdb",
        "species_rep10_raw_blastdb",
        "reference_local_hsps",
        "reference_intron_candidates",
        "reference_intron_audit",
        "reference_intron_blacklist",
    }
    if missing_optional & intron_reference_keys:
        steps.append(
            "Optional intron workflow references are missing. Build them after prepare-silva:\n"
            "  autotax2 build-intron-ref \\\n"
            "    --silva-fasta refdatabases/silva.fasta \\\n"
            "    --silva-taxonomy refdatabases/silva_taxonomy.tsv \\\n"
            "    --out refdatabases/intron_ref \\\n"
            "    --audit-reference-introns \\\n"
            "    --threads 4"
        )

    if "sina" in missing_optional or "sina_ref" in missing_all:
        steps.append(
            "Optional orientation correction with SINA is not fully configured:\n"
            "  autotax2 config set sina /path/to/sina\n"
            "  autotax2 config set sina-ref /path/to/SINA_reference.arb"
        )

    if "coordinate_ref" in missing_optional:
        steps.append(
            "Optional coordinate mapping reference is missing:\n"
            "  autotax2 config set coordinate-ref /path/to/coordinate_reference.fasta"
        )

    motif_missing = missing_optional & {"meme", "streme", "fimo"}
    if motif_missing:
        steps.append(
            "Optional motif analysis tools are missing. Install MEME Suite, then configure paths if needed:\n"
            "  autotax2 config set meme /path/to/meme\n"
            "  autotax2 config set streme /path/to/streme\n"
            "  autotax2 config set fimo /path/to/fimo"
        )

    structure_missing = missing_optional & {"rnafold", "rnaalifold"}
    if structure_missing:
        steps.append(
            "Optional RNA structure tools are missing. Install ViennaRNA, then configure paths if needed:\n"
            "  autotax2 config set rnafold /path/to/RNAfold\n"
            "  autotax2 config set rnaalifold /path/to/RNAalifold"
        )

    if "barrnap" in missing_optional:
        steps.append(
            "Optional barrnap audit is unavailable. Install barrnap or configure it:\n"
            "  autotax2 config set barrnap /path/to/barrnap"
        )

    if "typer" in missing_required:
        steps.append(
            "Install missing Python dependency:\n"
            "  pip install typer"
        )

    return steps


def print_plain_check_report(results: Sequence[CheckResult], title: str, ok: bool) -> None:
    grouped: Dict[str, List[CheckResult]] = {}

    for result in results:
        grouped.setdefault(result.section, []).append(result)

    print(title)
    print()

    for section, section_results in grouped.items():
        print(section)

        for result in section_results:
            status = status_label(result)
            required = required_label(result.required)
            detail = detail_for_result(result)
            action = action_for_key(result)

            print(f"  {result.key:<34} {status:<8} {required:<8} {detail}")
            if not result.ok:
                print(f"  {'':<34} {'':<8} {'':<8} -> {action}")

        print()

    print("Result")
    if ok:
        print("  OK: required checks passed.")
    else:
        print("  ERROR: one or more required checks failed.")

    steps = build_recommended_next_steps(results)
    if steps:
        print()
        print("Recommended next steps")
        for index, step in enumerate(steps, start=1):
            print(f"{index}. {step}")


def print_rich_check_report(results: Sequence[CheckResult], title: str, ok: bool) -> None:
    console = Console()
    grouped: Dict[str, List[CheckResult]] = {}

    for result in results:
        grouped.setdefault(result.section, []).append(result)

    panel_style = "green" if ok else "red"
    console.print(Panel.fit(f"[bold]{title}[/bold]", border_style=panel_style))

    for section, section_results in grouped.items():
        table = Table(
            title=section,
            title_style="bold cyan",
            box=box.SIMPLE_HEAVY,
            show_lines=False,
            expand=True,
            padding=(0, 1),
        )
        table.add_column("Key", style="bold", no_wrap=True)
        table.add_column("Status", no_wrap=True)
        table.add_column("Need", no_wrap=True)
        table.add_column("Detail", overflow="fold")
        table.add_column("Direction", overflow="fold")

        for result in section_results:
            status = Text(status_label(result), style=status_style(result))
            need = "required" if result.required else "optional"
            need_style = "bold red" if result.required and not result.ok else "dim"
            detail = detail_for_result(result)
            direction = "ready" if result.ok else action_for_key(result)

            if not result.ok and not result.required:
                detail = f"[dim]{detail}[/dim]"

            table.add_row(
                result.key,
                status,
                Text(need, style=need_style),
                detail,
                direction,
            )

        console.print(table)
        console.print()

    if ok:
        console.print("[bold green]Result[/bold green]")
        console.print("  [green]OK: required checks passed.[/green]")
    else:
        console.print("[bold red]Result[/bold red]")
        console.print("  [red]ERROR: one or more required checks failed.[/red]")

    steps = build_recommended_next_steps(results)
    if steps:
        console.print()
        console.print("[bold cyan]Recommended next steps[/bold cyan]")
        for index, step in enumerate(steps, start=1):
            console.print(f"[bold]{index}.[/bold] {step}")


def print_check_report(
    results: Iterable[CheckResult],
    title: str = "AutoTax2 check",
) -> bool:
    result_list = list(results)
    ok = required_checks_passed(result_list)

    if os.environ.get("NO_COLOR") or os.environ.get("AUTOTAX2_PLAIN_CHECK"):
        print_plain_check_report(result_list, title, ok)
    else:
        print_rich_check_report(result_list, title, ok)

    return ok


def check_all(
    *,
    ref_manifest: str | Path | None = None,
    source_map: str | Path | None = None,
    use_config: bool = True,
    config_path: str | Path | None = None,
) -> bool:
    results: List[CheckResult] = []

    if use_config:
        results.extend(check_from_config(config_path))
    else:
        for software in REQUIRED_SOFTWARE:
            results.append(
                check_software(
                    software,
                    configured_path=None,
                    required=True,
                )
            )

        for software in OPTIONAL_SOFTWARE:
            results.append(
                check_software(
                    software,
                    configured_path=None,
                    required=False,
                )
            )

        results.extend(check_python_packages())

    if ref_manifest:
        results.extend(check_reference_manifest(ref_manifest))

    if source_map:
        results.append(
            check_path(
                "source_map",
                source_map,
                required=False,
                section="Input files",
            )
        )

    return print_check_report(results)
