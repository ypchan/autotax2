from __future__ import annotations

import configparser
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .threads import DEFAULT_THREADS, ThreadValueError, parse_threads


APP_NAME = "autotax2"
ENV_CONFIG = "AUTOTAX2_CONFIG"

DEFAULT_CONFIG_PATH = Path.home() / ".config" / APP_NAME / "config.ini"

SOFTWARE_SECTION = "software"
REFERENCE_SECTION = "reference"
RUNTIME_SECTION = "runtime"


REQUIRED_SOFTWARE_KEYS = [
    "vsearch",
    "blastn",
    "makeblastdb",
]

OPTIONAL_SOFTWARE_KEYS = [
    "sina",
    "mafft",
    "meme",
    "streme",
    "fimo",
    "rnafold",
    "rnaalifold",
    "barrnap",
]

SOFTWARE_KEYS = REQUIRED_SOFTWARE_KEYS + OPTIONAL_SOFTWARE_KEYS

SOFTWARE_EXECUTABLE_NAMES = {
    "sina": "sina", 
    "vsearch": "vsearch",
    "blastn": "blastn",
    "makeblastdb": "makeblastdb",
    "mafft": "mafft",
    "meme": "meme",
    "streme": "streme",
    "fimo": "fimo",
    "rnafold": "RNAfold",
    "rnaalifold": "RNAalifold",
    "barrnap": "barrnap",
}


KEY_ALIASES: Dict[str, Tuple[str, str]] = {
    # Software
    "vsearch": (SOFTWARE_SECTION, "vsearch"),
    "sina": (SOFTWARE_SECTION, "sina"),
    "blastn": (SOFTWARE_SECTION, "blastn"),
    "makeblastdb": (SOFTWARE_SECTION, "makeblastdb"),
    "mafft": (SOFTWARE_SECTION, "mafft"),
    "meme": (SOFTWARE_SECTION, "meme"),
    "streme": (SOFTWARE_SECTION, "streme"),
    "fimo": (SOFTWARE_SECTION, "fimo"),
    "rnafold": (SOFTWARE_SECTION, "rnafold"),
    "rna-fold": (SOFTWARE_SECTION, "rnafold"),
    "rnaalifold": (SOFTWARE_SECTION, "rnaalifold"),
    "rna-alifold": (SOFTWARE_SECTION, "rnaalifold"),
    "barrnap": (SOFTWARE_SECTION, "barrnap"),

    # Core reference data
    "ref-manifest": (REFERENCE_SECTION, "ref_manifest"),
    "ref_manifest": (REFERENCE_SECTION, "ref_manifest"),
    "manifest": (REFERENCE_SECTION, "ref_manifest"),

    "silva-fasta": (REFERENCE_SECTION, "silva_fasta"),
    "silva_fasta": (REFERENCE_SECTION, "silva_fasta"),

    "silva-sintax": (REFERENCE_SECTION, "silva_sintax_fasta"),
    "silva-sintax-fasta": (REFERENCE_SECTION, "silva_sintax_fasta"),
    "silva_sintax": (REFERENCE_SECTION, "silva_sintax_fasta"),
    "silva_sintax_fasta": (REFERENCE_SECTION, "silva_sintax_fasta"),

    "typestrains": (REFERENCE_SECTION, "silva_typestrains_fasta"),
    "typestrain-fasta": (REFERENCE_SECTION, "silva_typestrains_fasta"),
    "typestrains-fasta": (REFERENCE_SECTION, "silva_typestrains_fasta"),
    "silva-typestrains-fasta": (REFERENCE_SECTION, "silva_typestrains_fasta"),
    "silva_typestrains_fasta": (REFERENCE_SECTION, "silva_typestrains_fasta"),

    "metadata": (REFERENCE_SECTION, "silva_metadata"),
    "silva-metadata": (REFERENCE_SECTION, "silva_metadata"),
    "silva_metadata": (REFERENCE_SECTION, "silva_metadata"),

    "silva-udb": (REFERENCE_SECTION, "silva_udb"),
    "silva_udb": (REFERENCE_SECTION, "silva_udb"),

    "silva-sintax-udb": (REFERENCE_SECTION, "silva_sintax_udb"),
    "silva_sintax_udb": (REFERENCE_SECTION, "silva_sintax_udb"),

    "typestrains-udb": (REFERENCE_SECTION, "silva_typestrains_udb"),
    "silva-typestrains-udb": (REFERENCE_SECTION, "silva_typestrains_udb"),
    "silva_typestrains_udb": (REFERENCE_SECTION, "silva_typestrains_udb"),

    "typestrains-accession-ids": (REFERENCE_SECTION, "typestrains_accession_ids"),
    "typestrains_accession_ids": (REFERENCE_SECTION, "typestrains_accession_ids"),

    "silva-taxonomy-tsv": (REFERENCE_SECTION, "silva_taxonomy_tsv"),
    "silva_taxonomy_tsv": (REFERENCE_SECTION, "silva_taxonomy_tsv"),

    # Intron-detection reference data: clean production reference
    "species-rep1-fasta": (REFERENCE_SECTION, "species_rep1_fasta"),
    "species_rep1_fasta": (REFERENCE_SECTION, "species_rep1_fasta"),
    "species-rep1-taxonomy": (REFERENCE_SECTION, "species_rep1_taxonomy"),
    "species_rep1_taxonomy": (REFERENCE_SECTION, "species_rep1_taxonomy"),
    "species-rep1-blastdb": (REFERENCE_SECTION, "species_rep1_blastdb"),
    "species_rep1_blastdb": (REFERENCE_SECTION, "species_rep1_blastdb"),

    "species-rep1-clean-fasta": (REFERENCE_SECTION, "species_rep1_clean_fasta"),
    "species_rep1_clean_fasta": (REFERENCE_SECTION, "species_rep1_clean_fasta"),
    "species-rep1-clean-taxonomy": (REFERENCE_SECTION, "species_rep1_clean_taxonomy"),
    "species_rep1_clean_taxonomy": (REFERENCE_SECTION, "species_rep1_clean_taxonomy"),
    "species-rep1-clean-blastdb": (REFERENCE_SECTION, "species_rep1_clean_blastdb"),
    "species_rep1_clean_blastdb": (REFERENCE_SECTION, "species_rep1_clean_blastdb"),

    # Intron-detection reference data: raw audit reference
    "species-rep10-fasta": (REFERENCE_SECTION, "species_rep10_fasta"),
    "species_rep10_fasta": (REFERENCE_SECTION, "species_rep10_fasta"),
    "species-rep10-taxonomy": (REFERENCE_SECTION, "species_rep10_taxonomy"),
    "species_rep10_taxonomy": (REFERENCE_SECTION, "species_rep10_taxonomy"),
    "species-rep10-blastdb": (REFERENCE_SECTION, "species_rep10_blastdb"),
    "species_rep10_blastdb": (REFERENCE_SECTION, "species_rep10_blastdb"),

    "species-rep10-raw-fasta": (REFERENCE_SECTION, "species_rep10_raw_fasta"),
    "species_rep10_raw_fasta": (REFERENCE_SECTION, "species_rep10_raw_fasta"),
    "species-rep10-raw-taxonomy": (REFERENCE_SECTION, "species_rep10_raw_taxonomy"),
    "species_rep10_raw_taxonomy": (REFERENCE_SECTION, "species_rep10_raw_taxonomy"),
    "species-rep10-raw-blastdb": (REFERENCE_SECTION, "species_rep10_raw_blastdb"),
    "species_rep10_raw_blastdb": (REFERENCE_SECTION, "species_rep10_raw_blastdb"),
    "species-rep10-raw-selection": (REFERENCE_SECTION, "species_rep10_raw_selection"),
    "species_rep10_raw_selection": (REFERENCE_SECTION, "species_rep10_raw_selection"),

    # Reference intron audit outputs
    "reference-intron-audit": (REFERENCE_SECTION, "reference_intron_audit"),
    "reference_intron_audit": (REFERENCE_SECTION, "reference_intron_audit"),
    "reference-intron-candidates": (REFERENCE_SECTION, "reference_intron_candidates"),
    "reference_intron_candidates": (REFERENCE_SECTION, "reference_intron_candidates"),
    "reference-intron-blacklist": (REFERENCE_SECTION, "reference_intron_blacklist"),
    "reference_intron_blacklist": (REFERENCE_SECTION, "reference_intron_blacklist"),
    "reference-local-hsps": (REFERENCE_SECTION, "reference_local_hsps"),
    "reference_local_hsps": (REFERENCE_SECTION, "reference_local_hsps"),
    "reference-barrnap-gff": (REFERENCE_SECTION, "reference_barrnap_gff"),
    "reference_barrnap_gff": (REFERENCE_SECTION, "reference_barrnap_gff"),
    "reference-barrnap-audit": (REFERENCE_SECTION, "reference_barrnap_audit"),
    "reference_barrnap_audit": (REFERENCE_SECTION, "reference_barrnap_audit"),

    "coordinate-ref": (REFERENCE_SECTION, "coordinate_ref"),
    "coordinate_ref": (REFERENCE_SECTION, "coordinate_ref"),

    "sina-ref": (REFERENCE_SECTION, "sina_ref"),
    "sina_ref": (REFERENCE_SECTION, "sina_ref"),
    "sina-reference": (REFERENCE_SECTION, "sina_ref"),
    "sina_reference": (REFERENCE_SECTION, "sina_ref"),      

    "threads": (RUNTIME_SECTION, "threads"),
}


MANIFEST_TO_CONFIG_KEYS: Dict[str, Tuple[str, str]] = {
    "silva_fasta": (REFERENCE_SECTION, "silva_fasta"),
    "silva_sintax_fasta": (REFERENCE_SECTION, "silva_sintax_fasta"),
    "silva_typestrains_fasta": (REFERENCE_SECTION, "silva_typestrains_fasta"),
    "silva_metadata": (REFERENCE_SECTION, "silva_metadata"),
    "typestrains_accession_ids": (REFERENCE_SECTION, "typestrains_accession_ids"),
    "silva_taxonomy_tsv": (REFERENCE_SECTION, "silva_taxonomy_tsv"),
    "silva_udb": (REFERENCE_SECTION, "silva_udb"),
    "silva_sintax_udb": (REFERENCE_SECTION, "silva_sintax_udb"),
    "silva_typestrains_udb": (REFERENCE_SECTION, "silva_typestrains_udb"),

    "species_rep1_fasta": (REFERENCE_SECTION, "species_rep1_fasta"),
    "species_rep1_taxonomy": (REFERENCE_SECTION, "species_rep1_taxonomy"),
    "species_rep1_blastdb": (REFERENCE_SECTION, "species_rep1_blastdb"),
    "species_rep1_clean_fasta": (REFERENCE_SECTION, "species_rep1_clean_fasta"),
    "species_rep1_clean_taxonomy": (REFERENCE_SECTION, "species_rep1_clean_taxonomy"),
    "species_rep1_clean_blastdb": (REFERENCE_SECTION, "species_rep1_clean_blastdb"),

    "species_rep10_fasta": (REFERENCE_SECTION, "species_rep10_fasta"),
    "species_rep10_taxonomy": (REFERENCE_SECTION, "species_rep10_taxonomy"),
    "species_rep10_blastdb": (REFERENCE_SECTION, "species_rep10_blastdb"),
    "species_rep10_raw_fasta": (REFERENCE_SECTION, "species_rep10_raw_fasta"),
    "species_rep10_raw_taxonomy": (REFERENCE_SECTION, "species_rep10_raw_taxonomy"),
    "species_rep10_raw_blastdb": (REFERENCE_SECTION, "species_rep10_raw_blastdb"),
    "species_rep10_raw_selection": (REFERENCE_SECTION, "species_rep10_raw_selection"),

    "reference_local_hsps": (REFERENCE_SECTION, "reference_local_hsps"),
    "reference_intron_candidates": (REFERENCE_SECTION, "reference_intron_candidates"),
    "reference_intron_audit": (REFERENCE_SECTION, "reference_intron_audit"),
    "reference_intron_blacklist": (REFERENCE_SECTION, "reference_intron_blacklist"),
    "reference_barrnap_gff": (REFERENCE_SECTION, "reference_barrnap_gff"),
    "reference_barrnap_audit": (REFERENCE_SECTION, "reference_barrnap_audit"),

    "coordinate_ref": (REFERENCE_SECTION, "coordinate_ref"),

    "sina_ref": (REFERENCE_SECTION, "sina_ref"),
}


REFERENCE_DEFAULTS: Dict[str, str] = {
    "ref_manifest": "",
    "silva_fasta": "",
    "silva_sintax_fasta": "",
    "silva_typestrains_fasta": "",
    "silva_metadata": "",
    "typestrains_accession_ids": "",
    "silva_taxonomy_tsv": "",
    "silva_udb": "",
    "silva_sintax_udb": "",
    "silva_typestrains_udb": "",
    "species_rep1_fasta": "",
    "species_rep1_taxonomy": "",
    "species_rep1_blastdb": "",
    "species_rep1_clean_fasta": "",
    "species_rep1_clean_taxonomy": "",
    "species_rep1_clean_blastdb": "",
    "species_rep10_fasta": "",
    "species_rep10_taxonomy": "",
    "species_rep10_blastdb": "",
    "species_rep10_raw_fasta": "",
    "species_rep10_raw_taxonomy": "",
    "species_rep10_raw_blastdb": "",
    "species_rep10_raw_selection": "",
    "reference_local_hsps": "",
    "reference_intron_candidates": "",
    "reference_intron_audit": "",
    "reference_intron_blacklist": "",
    "reference_barrnap_gff": "",
    "reference_barrnap_audit": "",
    "coordinate_ref": "",
    "sina_ref": "",
}


@dataclass(frozen=True)
class ConfigCheck:
    section: str
    key: str
    value: str
    required: bool
    ok: bool
    message: str


def get_config_path(path: Optional[str | Path] = None) -> Path:
    if path is not None:
        return Path(path).expanduser()

    env_path = os.environ.get(ENV_CONFIG)
    if env_path:
        return Path(env_path).expanduser()

    return DEFAULT_CONFIG_PATH


def new_config() -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config[SOFTWARE_SECTION] = {key: "" for key in SOFTWARE_KEYS}
    config[REFERENCE_SECTION] = dict(REFERENCE_DEFAULTS)
    config[RUNTIME_SECTION] = {"threads": str(DEFAULT_THREADS)}
    return config


def ensure_sections(config: configparser.ConfigParser) -> configparser.ConfigParser:
    defaults = new_config()

    for section in defaults.sections():
        if not config.has_section(section):
            config.add_section(section)

        for key, value in defaults[section].items():
            if key not in config[section]:
                config[section][key] = value

    return config


def load_config(path: Optional[str | Path] = None) -> configparser.ConfigParser:
    config_path = get_config_path(path)
    config = new_config()

    if config_path.exists():
        config.read(config_path, encoding="utf-8")

    return ensure_sections(config)


def save_config(
    config: configparser.ConfigParser,
    path: Optional[str | Path] = None,
) -> Path:
    config_path = get_config_path(path)
    config_path.parent.mkdir(parents=True, exist_ok=True)

    with config_path.open("w", encoding="utf-8") as handle:
        config.write(handle)

    return config_path


def init_config(
    path: Optional[str | Path] = None,
    overwrite: bool = False,
) -> Path:
    config_path = get_config_path(path)

    if config_path.exists() and not overwrite:
        return config_path

    config = new_config()

    for key in SOFTWARE_KEYS:
        executable = SOFTWARE_EXECUTABLE_NAMES.get(key, key)
        detected = shutil.which(executable)
        if detected:
            config[SOFTWARE_SECTION][key] = detected

    return save_config(config, config_path)


def normalize_key(key: str) -> Tuple[str, str]:
    cleaned = key.strip().replace("_", "-")
    lowered = cleaned.lower()

    if lowered not in KEY_ALIASES:
        allowed = ", ".join(sorted(KEY_ALIASES))
        raise KeyError(f"Unknown config key: {key}. Allowed keys: {allowed}")

    return KEY_ALIASES[lowered]


def set_config_value(
    key: str,
    value: str | Path,
    path: Optional[str | Path] = None,
) -> Path:
    config = load_config(path)
    section, option = normalize_key(key)

    if section == RUNTIME_SECTION and option == "threads":
        try:
            parse_threads(str(value))
        except ThreadValueError as exc:
            raise ValueError(str(exc)) from exc

    if not config.has_section(section):
        config.add_section(section)

    config[section][option] = str(value)
    return save_config(config, path)


def get_config_value(
    key: str,
    path: Optional[str | Path] = None,
    default: Optional[str] = None,
) -> Optional[str]:
    config = load_config(path)
    section, option = normalize_key(key)

    value = config.get(section, option, fallback="")
    if value:
        return value

    return default


def get_threads(
    path: Optional[str | Path] = None,
    default: int = DEFAULT_THREADS,
) -> int:
    value = get_config_value("threads", path=path, default=str(default))

    if value is None:
        return default

    return parse_threads(value)


def executable_exists(value: str) -> bool:
    if not value:
        return False

    expanded = Path(value).expanduser()

    if expanded.exists() and os.access(expanded, os.X_OK):
        return True

    return shutil.which(value) is not None


def get_software(
    key: str,
    path: Optional[str | Path] = None,
    required: bool = True,
) -> str:
    _, option = normalize_key(key)

    if option not in SOFTWARE_KEYS:
        raise KeyError(f"{key!r} is not a software configuration key.")

    configured = get_config_value(option, path=path)

    if configured:
        return configured

    executable = SOFTWARE_EXECUTABLE_NAMES.get(option, option)
    detected = shutil.which(executable)

    if detected:
        return detected

    if required:
        raise FileNotFoundError(
            f"{key} is not configured and was not found in PATH. "
            f"Run: autotax2 config set {key} /path/to/{executable}"
        )

    return executable


def get_vsearch(
    path: Optional[str | Path] = None,
    required: bool = True,
) -> str:
    return get_software("vsearch", path=path, required=required)


def get_ref_manifest(
    path: Optional[str | Path] = None,
    required: bool = True,
) -> Optional[Path]:
    value = get_config_value("ref-manifest", path=path)

    if value:
        return Path(value).expanduser()

    if required:
        raise FileNotFoundError(
            "Reference manifest is not configured. "
            "Run prepare-silva first or set it with: "
            "autotax2 config set ref-manifest refdatabases/autotax2_ref_manifest.tsv"
        )

    return None


def get_reference_path(
    key: str,
    path: Optional[str | Path] = None,
    required: bool = True,
) -> Optional[Path]:
    value = get_config_value(key, path=path)

    if value:
        return Path(value).expanduser()

    if required:
        raise FileNotFoundError(
            f"Reference path '{key}' is not configured. "
            "Run autotax2 config show to inspect the current configuration."
        )

    return None


def parse_manifest(manifest_path: str | Path) -> Dict[str, str]:
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


def update_reference_from_manifest(
    manifest_path: str | Path,
    path: Optional[str | Path] = None,
) -> Path:
    manifest_path = Path(manifest_path).expanduser()
    manifest_entries = parse_manifest(manifest_path)

    config = load_config(path)
    config[REFERENCE_SECTION]["ref_manifest"] = str(manifest_path)

    for manifest_key, value in manifest_entries.items():
        if manifest_key not in MANIFEST_TO_CONFIG_KEYS:
            continue

        section, option = MANIFEST_TO_CONFIG_KEYS[manifest_key]
        config[section][option] = value

    return save_config(config, path)


def iter_config_rows(
    path: Optional[str | Path] = None,
) -> Iterable[Tuple[str, str, str]]:
    config = load_config(path)

    for section in (SOFTWARE_SECTION, REFERENCE_SECTION, RUNTIME_SECTION):
        if not config.has_section(section):
            continue

        for key, value in config[section].items():
            yield section, key, value


def path_exists(value: str) -> bool:
    if not value:
        return False

    return Path(value).expanduser().exists()


def blastdb_exists(prefix: str) -> bool:
    if not prefix:
        return False

    prefix_path = Path(prefix).expanduser()
    candidates = [
        Path(str(prefix_path) + ".nhr"),
        Path(str(prefix_path) + ".nin"),
        Path(str(prefix_path) + ".nsq"),
        Path(str(prefix_path) + ".ndb"),
        Path(str(prefix_path) + ".not"),
        Path(str(prefix_path) + ".ntf"),
        Path(str(prefix_path) + ".nto"),
        Path(str(prefix_path) + ".00.nhr"),
        Path(str(prefix_path) + ".00.nin"),
        Path(str(prefix_path) + ".00.nsq"),
        Path(str(prefix_path) + ".00.ndb"),
        Path(str(prefix_path) + ".00.not"),
        Path(str(prefix_path) + ".00.ntf"),
        Path(str(prefix_path) + ".00.nto"),
    ]

    return any(path.exists() for path in candidates)


def _add_software_checks(config: configparser.ConfigParser, checks: List[ConfigCheck]) -> None:
    for key in REQUIRED_SOFTWARE_KEYS:
        configured = config.get(SOFTWARE_SECTION, key, fallback="")
        executable = SOFTWARE_EXECUTABLE_NAMES.get(key, key)
        detected = shutil.which(executable)
        ok = executable_exists(configured) or detected is not None

        checks.append(
            ConfigCheck(
                section=SOFTWARE_SECTION,
                key=key,
                value=configured or (detected or ""),
                required=True,
                ok=ok,
                message="found" if ok else "missing",
            )
        )

    for key in OPTIONAL_SOFTWARE_KEYS:
        configured = config.get(SOFTWARE_SECTION, key, fallback="")
        executable = SOFTWARE_EXECUTABLE_NAMES.get(key, key)
        detected = shutil.which(executable)
        ok = executable_exists(configured) or detected is not None

        checks.append(
            ConfigCheck(
                section=SOFTWARE_SECTION,
                key=key,
                value=configured or (detected or ""),
                required=False,
                ok=ok,
                message="found" if ok else "optional; not configured or missing",
            )
        )


def _add_path_checks(
    config: configparser.ConfigParser,
    checks: List[ConfigCheck],
    keys: Iterable[str],
    *,
    required: bool,
) -> None:
    for key in keys:
        value = config.get(REFERENCE_SECTION, key, fallback="")
        ok = path_exists(value) if value else False

        checks.append(
            ConfigCheck(
                section=REFERENCE_SECTION,
                key=key,
                value=value,
                required=required,
                ok=ok,
                message="found" if ok else ("missing" if required else "optional; not configured or missing"),
            )
        )


def _add_blastdb_checks(
    config: configparser.ConfigParser,
    checks: List[ConfigCheck],
    keys: Iterable[str],
    *,
    required: bool,
) -> None:
    for key in keys:
        value = config.get(REFERENCE_SECTION, key, fallback="")
        ok = blastdb_exists(value) if value else False

        checks.append(
            ConfigCheck(
                section=REFERENCE_SECTION,
                key=key,
                value=value,
                required=required,
                ok=ok,
                message="found" if ok else ("missing BLAST database files" if required else "optional; BLAST DB not configured or missing"),
            )
        )


def check_config(path: Optional[str | Path] = None) -> List[ConfigCheck]:
    config = load_config(path)
    checks: List[ConfigCheck] = []

    _add_software_checks(config, checks)

    required_reference_keys = [
        "ref_manifest",
        "silva_fasta",
        "silva_sintax_fasta",
        "silva_typestrains_fasta",
    ]

    optional_reference_keys = [
        "silva_metadata",
        "typestrains_accession_ids",
        "silva_taxonomy_tsv",
        "silva_udb",
        "silva_sintax_udb",
        "silva_typestrains_udb",
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

    optional_blastdb_keys = [
        "species_rep1_blastdb",
        "species_rep1_clean_blastdb",
        "species_rep10_blastdb",
        "species_rep10_raw_blastdb",
    ]

    _add_path_checks(config, checks, required_reference_keys, required=True)
    _add_path_checks(config, checks, optional_reference_keys, required=False)
    _add_blastdb_checks(config, checks, optional_blastdb_keys, required=False)

    threads_value = config.get(RUNTIME_SECTION, "threads", fallback=str(DEFAULT_THREADS))

    try:
        parse_threads(threads_value)
        threads_ok = True
        threads_message = "valid"
    except ThreadValueError as exc:
        threads_ok = False
        threads_message = str(exc)

    checks.append(
        ConfigCheck(
            section=RUNTIME_SECTION,
            key="threads",
            value=threads_value,
            required=True,
            ok=threads_ok,
            message=threads_message,
        )
    )

    return checks


def required_checks_passed(checks: Iterable[ConfigCheck]) -> bool:
    return all(check.ok for check in checks if check.required)


def config_path_message(path: Optional[str | Path] = None) -> str:
    return str(get_config_path(path))
