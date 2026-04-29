from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


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
    "sina_ref"
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


def status_label(ok: bool) -> str:
    return "OK" if ok else "MISSING"


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
    normalized_key = key.lower()
    executable_name = executable_name_for_key(normalized_key)

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
            key=normalized_key,
            ok=True,
            value=resolved,
            required=required,
            message=version if version else "found",
        )

    return CheckResult(
        section="Software",
        key=normalized_key,
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
        message="found" if ok else "missing",
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


def print_check_report(
    results: Iterable[CheckResult],
    title: str = "AutoTax2 check",
) -> bool:
    grouped: Dict[str, List[CheckResult]] = {}

    for result in results:
        grouped.setdefault(result.section, []).append(result)

    print(title)
    print()

    for section, section_results in grouped.items():
        print(section)

        for result in section_results:
            status = status_label(result.ok)
            required = required_label(result.required)

            if result.value and result.message:
                detail = f"{result.value} | {result.message}"
            elif result.value:
                detail = result.value
            else:
                detail = result.message

            print(f"  {result.key:<34} {status:<8} {required:<8} {detail}")

        print()

    ok = required_checks_passed(
        result
        for section_results in grouped.values()
        for result in section_results
    )

    print("Result")
    if ok:
        print("  OK: required checks passed.")
    else:
        print("  ERROR: one or more required checks failed.")

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
