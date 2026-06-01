"""Registry management for AutoTax2 databases."""
from __future__ import annotations

import json
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable

import pandas as pd
import yaml

from .config import DEFAULT_THRESHOLDS_PERCENT, RANKS_COARSE_TO_FINE, RANKS_FINE_TO_COARSE
from .fasta import parse_fasta
from .taxonomy import make_placeholder, normalize_taxon
from .utils import ensure_dir, require_file

logger = logging.getLogger("autotax2.registry")


@dataclass
class AutoTaxDB:
    root: Path

    @property
    def reference_dir(self) -> Path:
        return self.root / "reference"

    @property
    def registry_dir(self) -> Path:
        return self.root / "registry"

    @property
    def clusters_dir(self) -> Path:
        return self.root / "clusters"

    @property
    def centroids_dir(self) -> Path:
        return self.root / "centroids"

    @property
    def versions_dir(self) -> Path:
        return self.root / "versions"

    @property
    def logs_dir(self) -> Path:
        return self.root / "logs"

    @property
    def export_dir(self) -> Path:
        return self.root / "export"

    @property
    def config_path(self) -> Path:
        return self.root / "config.yaml"

    @property
    def sequence_registry_path(self) -> Path:
        return self.registry_dir / "sequence_registry.tsv"

    @property
    def cluster_registry_path(self) -> Path:
        return self.registry_dir / "cluster_registry.tsv"

    @property
    def placeholder_counter_path(self) -> Path:
        return self.registry_dir / "placeholder_counter.tsv"

    @property
    def source_registry_path(self) -> Path:
        return self.registry_dir / "source_registry.tsv"

    @property
    def current_fa(self) -> Path:
        return self.reference_dir / "current.fa"

    @property
    def current_taxonomy(self) -> Path:
        return self.reference_dir / "current.taxonomy.tsv"

    def create_dirs(self) -> None:
        for d in [
            self.root,
            self.reference_dir,
            self.registry_dir,
            self.clusters_dir,
            self.centroids_dir,
            self.versions_dir,
            self.logs_dir,
            self.export_dir,
        ]:
            ensure_dir(d)

    def exists(self) -> bool:
        return self.config_path.exists()

    def write_config(self, config: dict) -> None:
        self.create_dirs()
        with open(self.config_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(config, fh, sort_keys=False, allow_unicode=True)

    def read_config(self) -> dict:
        require_file(self.config_path, "AutoTax2 config")
        with open(self.config_path, "r", encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}

    def next_version_name(self, source: str) -> str:
        self.versions_dir.mkdir(parents=True, exist_ok=True)
        existing = sorted(p.name for p in self.versions_dir.glob("v[0-9][0-9][0-9]_*"))
        idx = len(existing)
        safe_source = source.replace("/", "_").replace(" ", "_")
        return f"v{idx:03d}_{safe_source}"


def initialize_placeholder_counter(path: Path) -> pd.DataFrame:
    df = pd.DataFrame({"rank": RANKS_FINE_TO_COARSE, "next_number": [1] * len(RANKS_FINE_TO_COARSE)})
    df.to_csv(path, sep="\t", index=False)
    return df


def load_placeholder_counter(path: Path) -> pd.DataFrame:
    if not path.exists():
        return initialize_placeholder_counter(path)
    return pd.read_csv(path, sep="\t")


def allocate_placeholder(db: AutoTaxDB, rank: str, prefix: str, digits: int = 6) -> tuple[str, int]:
    counter = load_placeholder_counter(db.placeholder_counter_path)
    row = counter[counter["rank"] == rank]
    if row.empty:
        raise ValueError(f"Unknown rank for placeholder allocation: {rank}")
    num = int(row.iloc[0]["next_number"])
    taxon = make_placeholder(rank, prefix, num, digits=digits)
    counter.loc[counter["rank"] == rank, "next_number"] = num + 1
    counter.to_csv(db.placeholder_counter_path, sep="\t", index=False)
    logger.debug("Allocated placeholder %s for rank=%s prefix=%s", taxon, rank, prefix)
    return taxon, num


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def save_table(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False)


def init_source_registry(db: AutoTaxDB, source: str = "gtdb", prefix: str = "gtdb") -> None:
    df = pd.DataFrame(
        [
            {
                "source": source,
                "prefix": prefix,
                "first_version": "v000_gtdb",
                "status": "active",
            }
        ]
    )
    save_table(df, db.source_registry_path)


def register_source(db: AutoTaxDB, source: str, prefix: str, version: str) -> None:
    df = load_table(db.source_registry_path)
    if df.empty or source not in set(df.get("source", [])):
        new = pd.DataFrame(
            [{"source": source, "prefix": prefix, "first_version": version, "status": "active"}]
        )
        df = pd.concat([df, new], ignore_index=True)
        save_table(df, db.source_registry_path)


def create_initial_registries(
    db: AutoTaxDB,
    ref_fa: Path,
    ref_tax: Path,
    ref_arb: Path,
    threads: int,
    force: bool = False,
) -> None:
    if db.root.exists() and any(db.root.iterdir()) and not force:
        raise FileExistsError(f"Database directory is not empty: {db.root}. Use --force to overwrite.")
    if force and db.root.exists():
        shutil.rmtree(db.root)
    db.create_dirs()

    shutil.copy2(ref_fa, db.reference_dir / "gtdb_ssu_tax.fa")
    shutil.copy2(ref_arb, db.reference_dir / "gtdb_ssu.arb")
    shutil.copy2(ref_tax, db.reference_dir / "gtdb_ssu.taxonomy.tsv")
    shutil.copy2(ref_fa, db.current_fa)

    tax = pd.read_csv(ref_tax, sep="\t", dtype=str).fillna("")
    required = ["seq_id", "domain", "phylum", "class", "order", "family", "genus", "species"]
    missing = [c for c in required if c not in tax.columns]
    if missing:
        raise ValueError(f"Missing required taxonomy columns: {missing}")

    seq_meta = []
    fasta_ids = set()
    for rec in parse_fasta(ref_fa):
        fasta_ids.add(rec.id)
        seq_meta.append(
            {
                "seq_id": rec.id,
                "seq_md5": rec.md5,
                "source": "gtdb",
                "source_seq_id": rec.id,
                "length": rec.length,
                "added_version": "v000_gtdb",
                "status": "active",
            }
        )
    seq_df = pd.DataFrame(seq_meta)
    seq_df = seq_df.merge(tax, on="seq_id", how="left")
    if seq_df["domain"].isna().any():
        missing_ids = seq_df.loc[seq_df["domain"].isna(), "seq_id"].head(10).tolist()
        raise ValueError(f"Some FASTA seq_ids are missing in taxonomy TSV, examples: {missing_ids}")
    save_table(seq_df.fillna(""), db.sequence_registry_path)

    cluster_rows = []
    membership_tables: Dict[str, list[dict]] = {rank: [] for rank in RANKS_FINE_TO_COARSE}
    for rank in RANKS_FINE_TO_COARSE:
        for taxon, sub in tax.groupby(rank, dropna=False):
            taxon = normalize_taxon(rank, taxon)
            if not taxon:
                continue
            cluster_id = f"{rank[0].upper()}{len([r for r in cluster_rows if r['rank'] == rank]) + 1:06d}"
            parent_rank = None
            parent_taxon = None
            if rank != "phylum":
                coarse = ["phylum", "class", "order", "family", "genus", "species"]
                idx = coarse.index(rank)
                parent_rank = coarse[idx - 1]
                parent_taxon = sub.iloc[0][parent_rank]
            rep_seq = choose_initial_representative(sub, seq_df)
            cluster_rows.append(
                {
                    "cluster_id": cluster_id,
                    "rank": rank,
                    "taxon_name": taxon,
                    "is_placeholder": "no",
                    "created_source": "gtdb",
                    "created_prefix": "gtdb",
                    "placeholder_number": "",
                    "parent_rank": parent_rank or "",
                    "parent_taxon": parent_taxon or "",
                    "representative_seq_id": rep_seq,
                    "created_version": "v000_gtdb",
                    "status": "active",
                    "merged_into": "",
                }
            )
            for sid in sub["seq_id"]:
                membership_tables[rank].append(
                    {
                        "seq_id": sid,
                        "source": "gtdb",
                        "rank": rank,
                        "cluster_id": cluster_id,
                        "taxon_name": taxon,
                        "is_representative": "yes" if sid == rep_seq else "no",
                        "added_version": "v000_gtdb",
                    }
                )
    save_table(pd.DataFrame(cluster_rows), db.cluster_registry_path)
    for rank, rows in membership_tables.items():
        save_table(pd.DataFrame(rows), db.clusters_dir / f"{rank}.membership.tsv")

    initialize_placeholder_counter(db.placeholder_counter_path)
    init_source_registry(db)
    save_table(seq_df[required], db.current_taxonomy)
    config = {
        "database": {"name": "AutoTax2-GTDB-SSU", "version": "v000_gtdb"},
        "reference": {
            "ref_fa": str(db.reference_dir / "gtdb_ssu_tax.fa"),
            "ref_tax": str(db.reference_dir / "gtdb_ssu.taxonomy.tsv"),
            "ref_arb": str(db.reference_dir / "gtdb_ssu.arb"),
        },
        "tools": {"sina": "sina", "vsearch": "vsearch"},
        "thresholds_percent": DEFAULT_THRESHOLDS_PERCENT,
        "vsearch": {"iddef": 2, "strand": "plus"},
        "placeholder": {"digits": 6, "numbering": "rank_global"},
        "threads_last_init": threads,
    }
    db.write_config(config)
    logger.info("Initialized AutoTax2 database at %s", db.root)


def choose_initial_representative(tax_subset: pd.DataFrame, seq_df: pd.DataFrame) -> str:
    sub = tax_subset.merge(seq_df[["seq_id", "length"]], on="seq_id", how="left")
    if "evidence_level" in sub.columns:
        priority = {
            "gtdb_type_complete_genome": 100,
            "gtdb_type_genome": 90,
            "gtdb_complete_genome": 80,
            "gtdb_representative_genome": 70,
            "gtdb_other_genome": 50,
        }
        sub["_priority"] = sub["evidence_level"].map(priority).fillna(0)
    else:
        sub["_priority"] = 0
    sub["length"] = pd.to_numeric(sub["length"], errors="coerce").fillna(0)
    sub = sub.sort_values(["_priority", "length", "seq_id"], ascending=[False, False, True])
    return str(sub.iloc[0]["seq_id"])
