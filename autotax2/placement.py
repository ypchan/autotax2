"""Near-best hit consensus and rank-aware placement."""

from __future__ import annotations

import csv
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from autotax2.placeholders import (
    PlaceholderAllocator,
    PlaceholderRank,
    make_placeholder_id,
    parse_placeholder_id,
)
from autotax2.vsearch import parse_uc_records


RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
LOW_TO_HIGH_RANKS = tuple(reversed(RANKS))
PLACEHOLDER_RANKS = ("class", "order", "family", "genus", "species")
PLACEHOLDER_RANK_BY_NAME = {
    "class": PlaceholderRank.CLASS,
    "order": PlaceholderRank.ORDER,
    "family": PlaceholderRank.FAMILY,
    "genus": PlaceholderRank.GENUS,
    "species": PlaceholderRank.SPECIES,
}
THRESHOLD_BY_STATUS = {
    "new_species": "species",
    "new_genus": "genus",
    "new_family": "family",
    "new_order": "order",
    "new_class": "class",
}
ASSIGNMENT_FIELDS = [
    "internal_seq_id",
    "original_seq_id",
    "dataset",
    "prefix",
    "sequence_md5",
    "is_duplicate_sequence",
    "best_hit_id",
    "best_hit_identity",
    "best_hit_taxon_id",
    "best_hit_taxonomy",
    "best_hit_source_category",
    "near_best_hit_count",
    "domain_consensus",
    "domain_consensus_fraction",
    "phylum_consensus",
    "phylum_consensus_fraction",
    "class_consensus",
    "class_consensus_fraction",
    "order_consensus",
    "order_consensus_fraction",
    "family_consensus",
    "family_consensus_fraction",
    "genus_consensus",
    "genus_consensus_fraction",
    "species_consensus",
    "species_consensus_fraction",
    "lowest_stable_rank",
    "identity_status",
    "final_status",
    "assigned_taxon_id",
    "assigned_taxonomy",
    "created_taxon_ids",
    "warning",
]
CREATED_TAXA_FIELDS = [
    "taxon_id",
    "rank",
    "name",
    "parent_taxon_id",
    "source",
    "source_prefix",
    "created_in_dataset",
    "cluster_key",
    "status",
    "representative_seq_id",
]
PLACEMENT_SUMMARY_FIELDS = [
    "dataset",
    "prefix",
    "input_representatives",
    "duplicate",
    "known_like",
    "new_species",
    "new_genus",
    "new_family",
    "new_order",
    "new_class",
    "ambiguous",
    "unplaced",
    "assigned_named_silva",
    "assigned_unresolved_silva",
    "assigned_previous_custom",
    "new_current_dataset",
    "created_species",
    "created_genera",
    "created_families",
    "created_orders",
    "created_classes",
]
CONSENSUS_FIELDS = [
    "internal_seq_id",
    "best_identity",
    "near_best_delta",
    "near_best_hit_count",
    "rank",
    "consensus_taxon_id",
    "consensus_name",
    "consensus_fraction",
    "hit_count",
    "stable",
]
REPRESENTATIVE_UPDATE_FIELDS = [
    "representative_seq_id",
    "taxon_id",
    "dataset",
    "action",
    "reason",
]


@dataclass(frozen=True)
class PlacementSummary:
    """Summary returned by :func:`place_dataset`."""

    dataset: str
    dataset_dir: Path
    assignments: int
    created_taxa: int
    dry_run: bool


@dataclass(frozen=True)
class PlacementHit:
    """One registry hit used by the placement engine."""

    query: str
    target: str
    identity: float
    taxon_id: str
    taxonomy: str


@dataclass(frozen=True)
class RankConsensus:
    """Consensus call at one rank for one query."""

    rank: str
    taxon_id: str
    name: str
    fraction: float
    hit_count: int
    stable: bool


@dataclass
class PlacementState:
    """Mutable placement state for one command run."""

    build_dir: Path
    registry_dir: Path
    dataset_dir: Path
    dataset: str
    prefix: str
    thresholds: dict[str, float]
    allocator: PlaceholderAllocator
    active_cluster_to_taxon: dict[str, str]
    taxon_rows: list[dict[str, str]]
    taxon_by_id: dict[str, dict[str, str]]
    created_taxa: list[dict[str, str]]
    representative_updates: list[dict[str, str]]


def placement_strategy() -> str:
    """Return the placement strategy."""
    return "near-best-hit-consensus"


def place_dataset(
    build: str | Path,
    dataset: str,
    near_best_delta: float = 0.005,
    rank_consensus: float = 0.80,
    species_id: float = 0.987,
    genus_id: float = 0.945,
    family_id: float = 0.865,
    order_id: float = 0.820,
    class_id: float = 0.785,
    floor_id: float = 0.750,
    dry_run: bool = False,
    allow_ambiguous: bool = True,
) -> PlacementSummary:
    """Place dataset representatives using near-best hit rank consensus."""
    build_dir = Path(build)
    registry_dir = build_dir / "registry"
    dataset_dir = _find_dataset_dir(build_dir, dataset)
    prefix = _dataset_prefix(build_dir, dataset)
    if not prefix:
        raise ValueError(f"Dataset prefix not found in registry for {dataset!r}.")

    thresholds = {
        "species": species_id,
        "genus": genus_id,
        "family": family_id,
        "order": order_id,
        "class": class_id,
        "floor": floor_id,
    }
    taxon_rows = _read_tsv(registry_dir / "taxon_nodes.tsv")
    state = PlacementState(
        build_dir=build_dir,
        registry_dir=registry_dir,
        dataset_dir=dataset_dir,
        dataset=dataset,
        prefix=prefix,
        thresholds=thresholds,
        allocator=_load_placeholder_allocator(registry_dir),
        active_cluster_to_taxon=_load_active_cluster_to_taxon(registry_dir),
        taxon_rows=taxon_rows,
        taxon_by_id={row.get("taxon_id", ""): dict(row) for row in taxon_rows if row.get("taxon_id")},
        created_taxa=[],
        representative_updates=[],
    )

    memberships = _membership_by_seq_id(dataset_dir)
    id_map = _id_map_by_seq_id(dataset_dir)
    duplicate_md5 = _existing_md5_to_taxon(registry_dir, current_dataset=dataset)
    representative_taxa = _representative_taxa(registry_dir)
    hits_by_query = _load_hits(dataset_dir / "vs_registry.filtered.tsv", representative_taxa, state)
    query_ids = _query_representatives(dataset_dir, memberships, hits_by_query, species_id)

    assignment_rows: list[dict[str, str]] = []
    consensus_rows: list[dict[str, str]] = []

    for query_id in query_ids:
        membership = memberships.get(query_id, {})
        id_row = id_map.get(query_id, {})
        md5 = membership.get("sequence_md5") or id_row.get("sequence_md5", "")
        is_local_duplicate = _is_true(membership.get("is_duplicate_sequence", "false"))
        duplicate_exists = md5 in duplicate_md5 if md5 else False
        duplicate_taxon_id = duplicate_md5.get(md5, "") if md5 else ""
        placement = _place_one_query(
            query_id=query_id,
            membership=membership,
            id_row=id_row,
            hits=hits_by_query.get(query_id, []),
            duplicate_taxon_id=duplicate_taxon_id,
            duplicate_exists=duplicate_exists,
            is_local_duplicate=is_local_duplicate,
            near_best_delta=near_best_delta,
            rank_consensus=rank_consensus,
            state=state,
        )
        if placement["final_status"] == "ambiguous" and not allow_ambiguous:
            raise RuntimeError(f"Ambiguous placement is not allowed: {query_id}")
        assignment_rows.append(placement)
        consensus_rows.extend(placement.pop("_consensus_rows"))

    suffix = ".dry_run.tsv" if dry_run else ".tsv"
    _write_tsv(assignment_rows, dataset_dir / f"assignments{suffix}", ASSIGNMENT_FIELDS)
    _write_tsv(state.created_taxa, dataset_dir / f"created_taxa{suffix}", CREATED_TAXA_FIELDS)
    _write_tsv(
        consensus_rows,
        dataset_dir / f"near_best_consensus{suffix}",
        CONSENSUS_FIELDS,
    )
    _write_tsv(
        state.representative_updates,
        dataset_dir / f"representative_updates{suffix}",
        REPRESENTATIVE_UPDATE_FIELDS,
    )
    _write_tsv(
        [_placement_summary_row(state, assignment_rows, query_ids)],
        dataset_dir / f"placement_summary{suffix}",
        PLACEMENT_SUMMARY_FIELDS,
    )

    if not dry_run:
        _update_registry(state, assignment_rows)

    return PlacementSummary(
        dataset=dataset,
        dataset_dir=dataset_dir,
        assignments=len(assignment_rows),
        created_taxa=len(state.created_taxa),
        dry_run=dry_run,
    )


def _place_one_query(
    query_id: str,
    membership: dict[str, str],
    id_row: dict[str, str],
    hits: list[PlacementHit],
    duplicate_taxon_id: str,
    duplicate_exists: bool,
    is_local_duplicate: bool,
    near_best_delta: float,
    rank_consensus: float,
    state: PlacementState,
) -> dict[str, Any]:
    sequence_md5 = membership.get("sequence_md5") or id_row.get("sequence_md5", "")
    original_seq_id = membership.get("original_seq_id") or id_row.get("original_seq_id", "")
    is_duplicate = is_local_duplicate or duplicate_exists
    if is_duplicate:
        assigned_taxon_id = duplicate_taxon_id
        assigned_taxonomy = _taxonomy_for_taxon(assigned_taxon_id, state) if assigned_taxon_id else ""
        row = _base_assignment_row(
            query_id=query_id,
            original_seq_id=original_seq_id,
            sequence_md5=sequence_md5,
            is_duplicate=True,
            state=state,
        )
        row.update(
            {
                "identity_status": "duplicate",
                "final_status": "duplicate",
                "assigned_taxon_id": assigned_taxon_id,
                "assigned_taxonomy": assigned_taxonomy,
            }
        )
        row["_consensus_rows"] = []
        return row

    if not hits:
        row = _base_assignment_row(
            query_id=query_id,
            original_seq_id=original_seq_id,
            sequence_md5=sequence_md5,
            is_duplicate=False,
            state=state,
        )
        row.update({"identity_status": "unplaced", "final_status": "unplaced"})
        row["_consensus_rows"] = []
        return row

    hits = sorted(hits, key=lambda hit: hit.identity, reverse=True)
    best_hit = hits[0]
    best_identity = best_hit.identity
    near_best_hits = [
        hit
        for hit in hits
        if hit.identity >= best_identity - near_best_delta
    ]
    consensus = _rank_consensus(near_best_hits, state, rank_consensus)
    identity_status = _identity_status(best_identity, state.thresholds)
    decision = _decide_placement(
        identity_status=identity_status,
        consensus=consensus,
        query_id=query_id,
        representative_seq_id=query_id,
        state=state,
    )
    row = _base_assignment_row(
        query_id=query_id,
        original_seq_id=original_seq_id,
        sequence_md5=sequence_md5,
        is_duplicate=False,
        state=state,
    )
    assigned_taxon_id = decision["assigned_taxon_id"]
    assigned_taxonomy = _taxonomy_for_taxon(assigned_taxon_id, state) if assigned_taxon_id else ""
    row.update(
        {
            "best_hit_id": best_hit.target,
            "best_hit_identity": f"{best_identity:.6f}",
            "best_hit_taxon_id": best_hit.taxon_id,
            "best_hit_taxonomy": _taxonomy_for_taxon(best_hit.taxon_id, state),
            "best_hit_source_category": _source_category(best_hit.taxon_id, state),
            "near_best_hit_count": str(len(near_best_hits)),
            "lowest_stable_rank": _lowest_stable_rank(consensus),
            "identity_status": identity_status,
            "final_status": decision["final_status"],
            "assigned_taxon_id": assigned_taxon_id,
            "assigned_taxonomy": assigned_taxonomy,
            "created_taxon_ids": ",".join(decision["created_taxon_ids"]),
            "warning": decision["warning"],
        }
    )
    for rank in RANKS:
        rank_call = consensus.get(rank)
        if rank_call is None:
            continue
        row[f"{rank}_consensus"] = rank_call.name
        row[f"{rank}_consensus_fraction"] = f"{rank_call.fraction:.6f}"

    row["_consensus_rows"] = [
        {
            "internal_seq_id": query_id,
            "best_identity": f"{best_identity:.6f}",
            "near_best_delta": f"{near_best_delta:.6f}",
            "near_best_hit_count": str(len(near_best_hits)),
            "rank": rank_call.rank,
            "consensus_taxon_id": rank_call.taxon_id,
            "consensus_name": rank_call.name,
            "consensus_fraction": f"{rank_call.fraction:.6f}",
            "hit_count": str(rank_call.hit_count),
            "stable": _bool_text(rank_call.stable),
        }
        for rank_call in consensus.values()
    ]
    return row


def _decide_placement(
    identity_status: str,
    consensus: dict[str, RankConsensus],
    query_id: str,
    representative_seq_id: str,
    state: PlacementState,
) -> dict[str, Any]:
    if identity_status == "unplaced":
        return _placement_decision("unplaced", "", [], "")

    if identity_status == "known_like":
        species = _stable_taxon(consensus, "species")
        if species:
            return _placement_decision("known_like", species, [], "")
        genus = _stable_taxon(consensus, "genus")
        if genus:
            created = _create_lineage(
                parent_taxon_id=genus,
                ranks=("species",),
                final_status="new_species",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision(
                "new_species",
                created[-1],
                created,
                "species_consensus_unstable_despite_high_identity",
            )
        family = _stable_taxon(consensus, "family")
        if family:
            created = _create_lineage(
                parent_taxon_id=family,
                ranks=("genus", "species"),
                final_status="new_genus",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision(
                "new_genus",
                created[-1],
                created,
                "species_and_genus_consensus_unstable",
            )
        return _placement_decision("ambiguous", "", [], "species_consensus_unstable")

    if identity_status == "new_species":
        genus = _stable_taxon(consensus, "genus")
        if genus:
            created = _create_lineage(
                parent_taxon_id=genus,
                ranks=("species",),
                final_status="new_species",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_species", created[-1], created, "")
        family = _stable_taxon(consensus, "family")
        if family:
            created = _create_lineage(
                parent_taxon_id=family,
                ranks=("genus", "species"),
                final_status="new_genus",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_genus", created[-1], created, "genus_consensus_unstable")
        return _placement_decision("ambiguous", "", [], "genus_consensus_unstable")

    if identity_status == "new_genus":
        family = _stable_taxon(consensus, "family")
        if family:
            created = _create_lineage(
                parent_taxon_id=family,
                ranks=("genus", "species"),
                final_status="new_genus",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_genus", created[-1], created, "")
        order = _stable_taxon(consensus, "order")
        if order:
            created = _create_lineage(
                parent_taxon_id=order,
                ranks=("family", "genus", "species"),
                final_status="new_family",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_family", created[-1], created, "family_consensus_unstable")
        return _placement_decision("ambiguous", "", [], "family_consensus_unstable")

    if identity_status == "new_family":
        order = _stable_taxon(consensus, "order")
        if order:
            created = _create_lineage(
                parent_taxon_id=order,
                ranks=("family", "genus", "species"),
                final_status="new_family",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_family", created[-1], created, "")
        class_taxon = _stable_taxon(consensus, "class")
        if class_taxon:
            created = _create_lineage(
                parent_taxon_id=class_taxon,
                ranks=("order", "family", "genus", "species"),
                final_status="new_order",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_order", created[-1], created, "order_consensus_unstable")
        return _placement_decision("ambiguous", "", [], "order_consensus_unstable")

    if identity_status == "new_order":
        class_taxon = _stable_taxon(consensus, "class")
        if class_taxon:
            created = _create_lineage(
                parent_taxon_id=class_taxon,
                ranks=("order", "family", "genus", "species"),
                final_status="new_order",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_order", created[-1], created, "")
        phylum = _stable_taxon(consensus, "phylum")
        if phylum:
            created = _create_lineage(
                parent_taxon_id=phylum,
                ranks=("class", "order", "family", "genus", "species"),
                final_status="new_class",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_class", created[-1], created, "class_consensus_unstable")
        return _placement_decision("ambiguous", "", [], "class_consensus_unstable")

    if identity_status == "new_class":
        phylum = _stable_taxon(consensus, "phylum")
        if phylum:
            created = _create_lineage(
                parent_taxon_id=phylum,
                ranks=("class", "order", "family", "genus", "species"),
                final_status="new_class",
                query_id=query_id,
                representative_seq_id=representative_seq_id,
                state=state,
            )
            return _placement_decision("new_class", created[-1], created, "")
        return _placement_decision("ambiguous", "", [], "phylum_consensus_unstable")

    return _placement_decision("ambiguous", "", [], "unsupported_identity_status")


def _create_lineage(
    parent_taxon_id: str,
    ranks: tuple[str, ...],
    final_status: str,
    query_id: str,
    representative_seq_id: str,
    state: PlacementState,
) -> list[str]:
    created_or_reused: list[str] = []
    current_parent = parent_taxon_id
    for rank in ranks:
        threshold_rank = rank if rank in state.thresholds else THRESHOLD_BY_STATUS[final_status]
        threshold = state.thresholds.get(threshold_rank, state.thresholds["floor"])
        cluster_key = f"{state.prefix}|{rank}|{current_parent}|{threshold:.3f}|{query_id}"
        existing_taxon_id = state.active_cluster_to_taxon.get(cluster_key)
        if existing_taxon_id:
            current_parent = existing_taxon_id
            created_or_reused.append(existing_taxon_id)
            continue

        name = state.allocator.allocate(PLACEHOLDER_RANK_BY_NAME[rank], state.prefix)
        taxon_id = name
        row = {
            "taxon_id": taxon_id,
            "rank": rank,
            "name": name,
            "parent_taxon_id": current_parent,
            "source": "custom_dataset",
            "source_prefix": state.prefix,
            "created_in_dataset": state.dataset,
            "cluster_key": cluster_key,
            "status": "active",
            "representative_seq_id": representative_seq_id if rank == "species" else "",
        }
        state.created_taxa.append(row)
        state.taxon_by_id[taxon_id] = {
            **row,
            "is_placeholder": "true",
            "protected": "false",
        }
        state.active_cluster_to_taxon[cluster_key] = taxon_id
        current_parent = taxon_id
        created_or_reused.append(taxon_id)
        if rank == "species":
            state.representative_updates.append(
                {
                    "representative_seq_id": representative_seq_id,
                    "taxon_id": taxon_id,
                    "dataset": state.dataset,
                    "action": "add",
                    "reason": final_status,
                }
            )
    return created_or_reused


def _rank_consensus(
    hits: list[PlacementHit],
    state: PlacementState,
    rank_consensus: float,
) -> dict[str, RankConsensus]:
    consensus: dict[str, RankConsensus] = {}
    hit_count = len(hits)
    for rank in RANKS:
        counts: Counter[str] = Counter()
        for hit in hits:
            lineage = _lineage_by_rank(hit.taxon_id, state)
            taxon_id = lineage.get(rank, "")
            if taxon_id:
                counts[taxon_id] += 1
        if not counts:
            consensus[rank] = RankConsensus(rank, "", "", 0.0, hit_count, False)
            continue
        taxon_id, count = counts.most_common(1)[0]
        fraction = count / hit_count if hit_count else 0.0
        consensus[rank] = RankConsensus(
            rank=rank,
            taxon_id=taxon_id,
            name=state.taxon_by_id.get(taxon_id, {}).get("name", taxon_id),
            fraction=fraction,
            hit_count=hit_count,
            stable=fraction >= rank_consensus,
        )
    return consensus


def _identity_status(identity: float, thresholds: dict[str, float]) -> str:
    if identity >= thresholds["species"]:
        return "known_like"
    if identity >= thresholds["genus"]:
        return "new_species"
    if identity >= thresholds["family"]:
        return "new_genus"
    if identity >= thresholds["order"]:
        return "new_family"
    if identity >= thresholds["class"]:
        return "new_order"
    if identity >= thresholds["floor"]:
        return "new_class"
    return "unplaced"


def _base_assignment_row(
    query_id: str,
    original_seq_id: str,
    sequence_md5: str,
    is_duplicate: bool,
    state: PlacementState,
) -> dict[str, Any]:
    row: dict[str, Any] = {field: "" for field in ASSIGNMENT_FIELDS}
    row.update(
        {
            "internal_seq_id": query_id,
            "original_seq_id": original_seq_id,
            "dataset": state.dataset,
            "prefix": state.prefix,
            "sequence_md5": sequence_md5,
            "is_duplicate_sequence": _bool_text(is_duplicate),
            "near_best_hit_count": "0",
        }
    )
    for rank in RANKS:
        row[f"{rank}_consensus_fraction"] = "0.000000"
    return row


def _placement_decision(
    final_status: str,
    assigned_taxon_id: str,
    created_taxon_ids: list[str],
    warning: str,
) -> dict[str, Any]:
    return {
        "final_status": final_status,
        "assigned_taxon_id": assigned_taxon_id,
        "created_taxon_ids": created_taxon_ids,
        "warning": warning,
    }


def _stable_taxon(consensus: dict[str, RankConsensus], rank: str) -> str:
    rank_call = consensus.get(rank)
    if rank_call is None or not rank_call.stable:
        return ""
    return rank_call.taxon_id


def _lowest_stable_rank(consensus: dict[str, RankConsensus]) -> str:
    for rank in LOW_TO_HIGH_RANKS:
        rank_call = consensus.get(rank)
        if rank_call is not None and rank_call.stable:
            return rank
    return ""


def _taxonomy_for_taxon(taxon_id: str, state: PlacementState) -> str:
    if not taxon_id:
        return ""
    labels: list[str] = []
    seen: set[str] = set()
    current = taxon_id
    while current and current not in seen:
        seen.add(current)
        row = state.taxon_by_id.get(current, {})
        if not row:
            labels.append(current)
            break
        labels.append(row.get("name", current))
        current = row.get("parent_taxon_id", "")
    return ";".join(reversed(labels))


def _lineage_by_rank(taxon_id: str, state: PlacementState) -> dict[str, str]:
    lineage: dict[str, str] = {}
    seen: set[str] = set()
    current = taxon_id
    while current and current not in seen:
        seen.add(current)
        row = state.taxon_by_id.get(current, {})
        if not row:
            break
        rank = row.get("rank", "")
        if rank in RANKS:
            lineage[rank] = current
        current = row.get("parent_taxon_id", "")
    return lineage


def _load_hits(
    path: Path,
    representative_taxa: dict[str, str],
    state: PlacementState,
) -> dict[str, list[PlacementHit]]:
    hits_by_query: dict[str, list[PlacementHit]] = {}
    for row in _read_tsv(path):
        query = row.get("query", "")
        target = row.get("target", "")
        if not query or not target:
            continue
        taxon_id = (
            row.get("target_taxon_id")
            or row.get("taxon_id")
            or representative_taxa.get(target, "")
        )
        if not taxon_id and target in state.taxon_by_id:
            taxon_id = target
        if not taxon_id:
            continue
        hit = PlacementHit(
            query=query,
            target=target,
            identity=_identity_fraction(row.get("identity", "0")),
            taxon_id=taxon_id,
            taxonomy=_taxonomy_for_taxon(taxon_id, state),
        )
        hits_by_query.setdefault(query, []).append(hit)
    return hits_by_query


def _identity_fraction(value: str) -> float:
    try:
        identity = float(value)
    except (TypeError, ValueError):
        return 0.0
    if identity > 1.0:
        return identity / 100.0
    return identity


def _query_representatives(
    dataset_dir: Path,
    memberships: dict[str, dict[str, str]],
    hits_by_query: dict[str, list[PlacementHit]],
    species_id: float,
) -> list[str]:
    cluster_path = dataset_dir / "internal_clusters" / f"species_{species_id:.3f}.uc"
    query_ids: set[str] = set()
    if cluster_path.exists():
        for record in parse_uc_records(cluster_path):
            if record.record_type == "S":
                query_ids.add(record.query_label)
    if not query_ids:
        query_ids.update(memberships)
    query_ids.update(hits_by_query)
    return sorted(query_ids)


def _membership_by_seq_id(dataset_dir: Path) -> dict[str, dict[str, str]]:
    return {
        row.get("internal_seq_id", ""): row
        for row in _read_tsv(dataset_dir / "sequence_membership.tsv")
        if row.get("internal_seq_id")
    }


def _id_map_by_seq_id(dataset_dir: Path) -> dict[str, dict[str, str]]:
    return {
        row.get("internal_seq_id", ""): row
        for row in _read_tsv(dataset_dir / "sequence_id_map.tsv")
        if row.get("internal_seq_id")
    }


def _representative_taxa(registry_dir: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for path in (
        registry_dir / "representative_registry.tsv",
        registry_dir / "representative_history.tsv",
        registry_dir / "sequence_registry.tsv",
    ):
        for row in _read_tsv(path):
            representative_seq_id = _first_value(
                row,
                "representative_seq_id",
                "seq_id",
                "internal_seq_id",
                "target",
            )
            taxon_id = _first_value(
                row,
                "taxon_id",
                "assigned_taxon_id",
                "species_taxon_id",
                "representative_taxon_id",
            )
            if representative_seq_id and taxon_id:
                mapping[representative_seq_id] = taxon_id
    return mapping


def _existing_md5_to_taxon(registry_dir: Path, current_dataset: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in _read_tsv(registry_dir / "sequence_registry.tsv"):
        row_dataset = row.get("dataset") or row.get("created_in_dataset", "")
        if row_dataset == current_dataset:
            continue
        digest = row.get("sequence_md5", "")
        if not digest:
            continue
        mapping[digest] = _first_value(row, "taxon_id", "assigned_taxon_id", "species_taxon_id")
    return mapping


def _source_category(taxon_id: str, state: PlacementState) -> str:
    row = state.taxon_by_id.get(taxon_id, {})
    if not row:
        return "none"
    created_in_dataset = row.get("created_in_dataset", "")
    if created_in_dataset == state.dataset:
        return "current_dataset"
    if created_in_dataset:
        return "previous_custom"
    if _is_true(row.get("is_silva_unresolved", "false")):
        return "unresolved_silva"
    if (
        _is_true(row.get("is_silva_named", "false"))
        or row.get("source", "").startswith("SILVA")
        or _is_true(row.get("protected", "false"))
    ):
        return "named_silva"
    return "none"


def _placement_summary_row(
    state: PlacementState,
    assignments: list[dict[str, str]],
    query_ids: list[str],
) -> dict[str, str]:
    status_counts = Counter(row.get("final_status", "") for row in assignments)
    source_counts = Counter(
        _source_category(row.get("assigned_taxon_id", ""), state)
        for row in assignments
        if row.get("assigned_taxon_id")
    )
    created_counts = Counter(row.get("rank", "") for row in state.created_taxa)
    return {
        "dataset": state.dataset,
        "prefix": state.prefix,
        "input_representatives": str(len(query_ids)),
        "duplicate": str(status_counts.get("duplicate", 0)),
        "known_like": str(status_counts.get("known_like", 0)),
        "new_species": str(status_counts.get("new_species", 0)),
        "new_genus": str(status_counts.get("new_genus", 0)),
        "new_family": str(status_counts.get("new_family", 0)),
        "new_order": str(status_counts.get("new_order", 0)),
        "new_class": str(status_counts.get("new_class", 0)),
        "ambiguous": str(status_counts.get("ambiguous", 0)),
        "unplaced": str(status_counts.get("unplaced", 0)),
        "assigned_named_silva": str(source_counts.get("named_silva", 0)),
        "assigned_unresolved_silva": str(source_counts.get("unresolved_silva", 0)),
        "assigned_previous_custom": str(source_counts.get("previous_custom", 0)),
        "new_current_dataset": str(source_counts.get("current_dataset", 0)),
        "created_species": str(created_counts.get("species", 0)),
        "created_genera": str(created_counts.get("genus", 0)),
        "created_families": str(created_counts.get("family", 0)),
        "created_orders": str(created_counts.get("order", 0)),
        "created_classes": str(created_counts.get("class", 0)),
    }


def _update_registry(state: PlacementState, assignment_rows: list[dict[str, str]]) -> None:
    _update_taxon_nodes(state)
    _update_name_index(state)
    _update_cluster_to_taxon(state)
    _write_placeholder_counters_yaml(state)
    _update_representative_registry(state)
    _update_representative_history(state)
    _update_sequence_registry(state, assignment_rows)


def _update_taxon_nodes(state: PlacementState) -> None:
    path = state.registry_dir / "taxon_nodes.tsv"
    existing = _read_tsv(path)
    rows_by_taxon_id = {row.get("taxon_id", ""): dict(row) for row in existing if row.get("taxon_id")}
    for taxon in state.created_taxa:
        row = dict(taxon)
        row.update(
            {
                "taxonomy_path": _taxonomy_for_taxon(taxon["taxon_id"], state),
                "protected": "false",
                "is_silva_named": "false",
                "is_placeholder": "true",
                "is_silva_unresolved": "false",
            }
        )
        rows_by_taxon_id[row["taxon_id"]] = row
    fieldnames = _merge_fieldnames(
        [
            "taxon_id",
            "rank",
            "name",
            "parent_taxon_id",
            "taxonomy_path",
            "protected",
            "is_silva_named",
            "source",
            "source_prefix",
            "created_in_dataset",
            "is_placeholder",
            "is_silva_unresolved",
            "status",
            "cluster_key",
            "representative_seq_id",
        ],
        list(rows_by_taxon_id.values()),
    )
    _write_tsv(list(rows_by_taxon_id.values()), path, fieldnames)


def _update_name_index(state: PlacementState) -> None:
    path = state.registry_dir / "name_index.tsv"
    rows = _read_tsv(path)
    keys = {(row.get("name", ""), row.get("rank", "")) for row in rows}
    for taxon in state.created_taxa:
        key = (taxon["name"], taxon["rank"])
        if key in keys:
            continue
        rows.append(
            {
                "name": taxon["name"],
                "rank": taxon["rank"],
                "taxon_id": taxon["taxon_id"],
                "protected": "false",
                "source": "custom_dataset",
                "source_prefix": state.prefix,
                "status": "active",
            }
        )
        keys.add(key)
    fieldnames = _merge_fieldnames(
        ["name", "rank", "taxon_id", "protected", "source", "source_prefix", "status"],
        rows,
    )
    _write_tsv(rows, path, fieldnames)


def _update_cluster_to_taxon(state: PlacementState) -> None:
    path = state.registry_dir / "cluster_to_taxon.tsv"
    rows = _read_tsv(path)
    rows_by_key = {row.get("cluster_key", ""): dict(row) for row in rows if row.get("cluster_key")}
    for taxon in state.created_taxa:
        rows_by_key[taxon["cluster_key"]] = {
            "cluster_key": taxon["cluster_key"],
            "taxon_id": taxon["taxon_id"],
            "rank": taxon["rank"],
            "name": taxon["name"],
            "status": "active",
            "source_prefix": state.prefix,
            "created_in_dataset": state.dataset,
        }
    fieldnames = _merge_fieldnames(
        ["cluster_key", "taxon_id", "rank", "name", "status", "source_prefix", "created_in_dataset"],
        list(rows_by_key.values()),
    )
    _write_tsv(list(rows_by_key.values()), path, fieldnames)


def _write_placeholder_counters_yaml(state: PlacementState) -> None:
    path = state.registry_dir / "placeholder_counters.yaml"
    counters = _read_yaml_mapping(path)
    prefix_counters = counters.setdefault(state.prefix, {})
    for rank in PLACEHOLDER_RANKS:
        prefix_counters[rank] = state.allocator.next_ordinal(
            PLACEHOLDER_RANK_BY_NAME[rank],
            state.prefix,
        )
    path.write_text(
        yaml.safe_dump(counters, sort_keys=True),
        encoding="utf-8",
    )


def _update_representative_registry(state: PlacementState) -> None:
    path = state.registry_dir / "representative_registry.tsv"
    rows = _read_tsv(path)
    keys = {(row.get("representative_seq_id", ""), row.get("taxon_id", "")) for row in rows}
    for update in state.representative_updates:
        key = (update["representative_seq_id"], update["taxon_id"])
        if key in keys:
            continue
        rows.append(
            {
                "representative_seq_id": update["representative_seq_id"],
                "taxon_id": update["taxon_id"],
                "dataset": state.dataset,
                "source_category": "current_dataset",
                "status": "active",
                "protected": "false",
            }
        )
        keys.add(key)
    fieldnames = _merge_fieldnames(
        [
            "representative_seq_id",
            "taxon_id",
            "dataset",
            "source_category",
            "status",
            "protected",
        ],
        rows,
    )
    _write_tsv(rows, path, fieldnames)


def _update_representative_history(state: PlacementState) -> None:
    path = state.registry_dir / "representative_history.tsv"
    rows = _read_tsv(path)
    rows.extend(state.representative_updates)
    fieldnames = _merge_fieldnames(REPRESENTATIVE_UPDATE_FIELDS, rows)
    _write_tsv(rows, path, fieldnames)


def _update_sequence_registry(state: PlacementState, assignment_rows: list[dict[str, str]]) -> None:
    path = state.registry_dir / "sequence_registry.tsv"
    rows = _read_tsv(path)
    keys = {(row.get("seq_id", "") or row.get("internal_seq_id", ""), row.get("dataset", "")) for row in rows}
    for assignment in assignment_rows:
        key = (assignment["internal_seq_id"], state.dataset)
        if key in keys:
            continue
        rows.append(
            {
                "seq_id": assignment["internal_seq_id"],
                "internal_seq_id": assignment["internal_seq_id"],
                "original_seq_id": assignment["original_seq_id"],
                "dataset": state.dataset,
                "prefix": state.prefix,
                "sequence_md5": assignment["sequence_md5"],
                "taxon_id": assignment["assigned_taxon_id"],
                "assigned_taxon_id": assignment["assigned_taxon_id"],
                "final_status": assignment["final_status"],
                "is_duplicate_sequence": assignment["is_duplicate_sequence"],
            }
        )
        keys.add(key)
    fieldnames = _merge_fieldnames(
        [
            "seq_id",
            "internal_seq_id",
            "original_seq_id",
            "dataset",
            "prefix",
            "sequence_md5",
            "taxon_id",
            "assigned_taxon_id",
            "final_status",
            "is_duplicate_sequence",
        ],
        rows,
    )
    _write_tsv(rows, path, fieldnames)


def _load_placeholder_allocator(registry_dir: Path) -> PlaceholderAllocator:
    issued: set[str] = set()
    deprecated: set[str] = set()
    for row in _read_tsv(registry_dir / "taxon_nodes.tsv"):
        name = row.get("name", "")
        try:
            parse_placeholder_id(name)
        except ValueError:
            continue
        if row.get("status", "active") in {"deprecated", "superseded"}:
            deprecated.add(name)
        else:
            issued.add(name)
    for row in _read_tsv(registry_dir / "cluster_to_taxon.tsv"):
        name = row.get("name", "")
        try:
            parse_placeholder_id(name)
        except ValueError:
            continue
        if row.get("status", "active") in {"deprecated", "superseded"}:
            deprecated.add(name)
        else:
            issued.add(name)

    counters = _read_yaml_mapping(registry_dir / "placeholder_counters.yaml")
    for prefix, rank_values in counters.items():
        if not isinstance(rank_values, dict):
            continue
        for rank, next_ordinal in rank_values.items():
            if rank not in PLACEHOLDER_RANK_BY_NAME:
                continue
            try:
                ordinal = int(next_ordinal)
            except (TypeError, ValueError):
                continue
            if ordinal > 1:
                issued.add(make_placeholder_id(PLACEHOLDER_RANK_BY_NAME[rank], prefix, ordinal - 1))

    for row in _read_tsv(registry_dir / "placeholder_counters.tsv"):
        rank = row.get("rank", "")
        next_ordinal = row.get("next_ordinal", "")
        prefix = row.get("source_prefix") or row.get("prefix") or row.get("dataset_prefix", "")
        if rank not in PLACEHOLDER_RANK_BY_NAME or not prefix or not next_ordinal.isdigit():
            continue
        ordinal = int(next_ordinal)
        if ordinal > 1:
            issued.add(make_placeholder_id(PLACEHOLDER_RANK_BY_NAME[rank], prefix, ordinal - 1))

    return PlaceholderAllocator(issued_ids=issued, deprecated_ids=deprecated)


def _load_active_cluster_to_taxon(registry_dir: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in _read_tsv(registry_dir / "cluster_to_taxon.tsv"):
        if row.get("status", "active") != "active":
            continue
        cluster_key = row.get("cluster_key", "")
        taxon_id = row.get("taxon_id", "")
        if cluster_key and taxon_id:
            mapping[cluster_key] = taxon_id
    return mapping


def _find_dataset_dir(build_dir: Path, dataset: str) -> Path:
    datasets_dir = build_dir / "datasets"
    if not datasets_dir.exists():
        raise FileNotFoundError(f"Dataset directory root not found: {datasets_dir}")
    matches = sorted(
        path
        for path in datasets_dir.iterdir()
        if path.is_dir() and path.name.split("_", maxsplit=1)[-1] == dataset
    )
    if not matches:
        raise FileNotFoundError(f"Prepared dataset not found: {dataset}")
    if len(matches) > 1:
        raise ValueError(f"Multiple prepared dataset directories match {dataset!r}.")
    return matches[0]


def _dataset_prefix(build_dir: Path, dataset: str) -> str:
    for row in _read_tsv(build_dir / "registry" / "dataset_registry.tsv"):
        if row.get("dataset_name") == dataset:
            return row.get("prefix", "")
    return ""


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))


def _write_tsv(rows: Iterable[dict[str, Any]], path: Path, fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=chr(9),
            lineterminator=chr(10),
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _merge_fieldnames(preferred: list[str], rows: list[dict[str, str]]) -> list[str]:
    fieldnames = list(preferred)
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    return fieldnames


def _read_yaml_mapping(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


def _first_value(row: dict[str, str], *fields: str) -> str:
    for field in fields:
        value = row.get(field, "")
        if value:
            return value
    return ""


def _is_true(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "y"}


def _bool_text(value: bool) -> str:
    return "true" if value else "false"
