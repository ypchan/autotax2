from __future__ import annotations

import csv
import shutil
import uuid
from collections.abc import Iterator
from pathlib import Path

import pytest
import yaml

from autotax2.reports import summarize_build
from autotax2.registry import sequence_md5
from autotax2.validate import validate_build


@pytest.fixture
def report_tmp_dir() -> Iterator[Path]:
    tmp_dir = Path("tests") / f"tmp_reports_{uuid.uuid4().hex}"
    tmp_dir.mkdir(parents=True)
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_summarize_writes_all_report_files(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summary = summarize_build(build)

    expected = {
        "global_summary.tsv",
        "dataset_delta_summary.tsv",
        "dataset_overlap_matrix.tsv",
        "rank_novelty_summary.tsv",
        "dataset_rank_overlap_detail.tsv",
        "dataset_rank_novelty_detail.tsv",
        "source_contribution.tsv",
        "representative_summary.tsv",
        "sequence_dedup_summary.tsv",
    }
    assert summary.files_written == 10
    assert expected == {path.name for path in (build / "reports").glob("*.tsv")}
    assert (build / "reports" / "dataset_increment_audit.md").exists()


def test_global_summary_has_expected_counts(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summarize_build(build)

    row = _read_tsv(build / "reports" / "global_summary.tsv")[0]
    assert row["custom_datasets"] == "1"
    assert row["active_species"] == "3"
    assert row["active_representatives"] == "3"
    assert row["silva_unresolved_evidence_records"] == "1"
    assert row["silva_unresolved_evidence_rows"] == "6"
    assert int(row["duplicate_sequences"]) >= 1


def test_dataset_delta_separates_named_and_unresolved_sources(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summarize_build(build)

    row = _read_tsv(build / "reports" / "dataset_delta_summary.tsv")[0]
    assert row["assigned_named_silva"] == "1"
    assert row["assigned_unresolved_silva"] == "1"
    assert row["new_current_dataset"] == "1"
    assert row["sina_candidate_queries"] == "2"
    assert row["sina_candidate_targets"] == "2"
    assert row["sina_candidate_target_matches"] == "2"
    assert row["placement_evidence_rows"] == "18"


def test_dataset_overlap_matrix_includes_exact_sequence_and_species(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summarize_build(build)

    rows = _read_tsv(build / "reports" / "dataset_overlap_matrix.tsv")
    ranks = {row["rank"] for row in rows}
    assert "exact_sequence" in ranks
    assert "species" in ranks
    assert any(row["compared_source_type"] == "named_silva" for row in rows)


def test_dataset_rank_overlap_detail_lists_rank_taxa_and_relations(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summarize_build(build)

    rows = _read_tsv(build / "reports" / "dataset_rank_overlap_detail.tsv")
    species_rows = [row for row in rows if row["rank"] == "species"]

    assert any(row["rank_taxon_id"] == "S1" and row["relation"] == "overlap_named_silva" for row in species_rows)
    assert any(row["rank_taxon_id"] == "s__SILVAs000001" and row["relation"] == "overlap_unresolved_silva" for row in species_rows)
    assert any(row["rank_taxon_id"] == "s__D20s000001" and row["relation"] == "new_current_dataset" for row in species_rows)
    new_species = next(row for row in species_rows if row["rank_taxon_id"] == "s__D20s000001")
    assert new_species["sequence_ids"] == "D20_000003"
    assert new_species["evidence_decisions"] == "assign_existing:1"


def test_dataset_rank_novelty_detail_lists_supporting_sequences(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summarize_build(build)

    rows = _read_tsv(build / "reports" / "dataset_rank_novelty_detail.tsv")
    by_taxon = {row["taxon_id"]: row for row in rows}

    assert by_taxon["g__D20g000001"]["rank"] == "genus"
    assert by_taxon["g__D20g000001"]["supporting_sequence_ids"] == "D20_000003"
    assert by_taxon["s__D20s000001"]["rank"] == "species"
    assert by_taxon["s__D20s000001"]["representative_seq_id"] == "D20_000003"
    assert by_taxon["s__D20s000001"]["is_placeholder"] == "true"


def test_dataset_increment_audit_markdown_summarizes_increment(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summarize_build(build)

    content = (build / "reports" / "dataset_increment_audit.md").read_text(encoding="utf-8")
    assert "# autotax2 Dataset Increment Audit" in content
    assert "## Dataset digester2020 (D20)" in content
    assert "overlap_named_silva" not in content
    assert "| species | 1 | 1 | 0 | 1 | 0 | 0 |" in content
    assert "s__D20s000001" in content
    assert "D20_000003" in content
    assert "`placement_evidence.tsv`" in content


def test_duplicate_placeholder_detection_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _append_tsv(
        build / "registry" / "taxon_nodes.tsv",
        _taxon("DUP_S", "species", "s__D20s000001", "g__D20g000001", is_placeholder="true", created="digester2020"),
    )

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "Duplicate active")


def test_reused_deprecated_placeholder_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _append_tsv(
        build / "registry" / "taxon_nodes.tsv",
        _taxon("OLD_DEP", "species", "s__D20s000009", "g__D20g000001", status="deprecated", is_placeholder="true", created="digester2020"),
    )
    _append_tsv(
        build / "registry" / "taxon_nodes.tsv",
        _taxon("OLD_ACTIVE", "species", "s__D20s000009", "g__D20g000001", is_placeholder="true", created="digester2020"),
    )

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "reused")


def test_duplicate_dataset_prefix_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _append_tsv(
        build / "registry" / "dataset_registry.tsv",
        {
            "dataset_name": "other",
            "prefix": "D20",
            "add_order": "2",
            "dataset_dir": str(build / "datasets" / "02_other"),
        },
    )

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "Prefix D20")


def test_invalid_placeholder_format_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _append_tsv(
        build / "registry" / "taxon_nodes.tsv",
        _taxon("BADP", "genus", "g__D20x000001", "F1", is_placeholder="true", created="digester2020"),
    )

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "Invalid placeholder format")


def test_missing_parent_taxon_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    rows = _read_tsv(build / "registry" / "taxon_nodes.tsv")
    for row in rows:
        if row["taxon_id"] == "S1":
            row["parent_taxon_id"] = "NO_SUCH_PARENT"
    _write_tsv(build / "registry" / "taxon_nodes.tsv", rows, list(rows[0]))

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "missing parent")


def test_taxon_cycle_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    rows = _read_tsv(build / "registry" / "taxon_nodes.tsv")
    for row in rows:
        if row["taxon_id"] == "G1":
            row["parent_taxon_id"] = "S1"
    _write_tsv(build / "registry" / "taxon_nodes.tsv", rows, list(rows[0]))

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "cycle")


def test_missing_original_id_mapping_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    rows = _read_tsv(build / "datasets" / "01_digester2020" / "sequence_id_map.tsv")
    rows[0]["original_seq_id"] = ""
    _write_tsv(build / "datasets" / "01_digester2020" / "sequence_id_map.tsv", rows, list(rows[0]))

    summary = validate_build(build, check_exports=False)

    assert summary.failed
    assert _has_error(build, "incomplete ID mapping")


def test_sintax_export_with_rank_prefix_values_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    path = build / "export" / "sintax"
    path.mkdir(parents=True)
    (path / "bad.fa").write_text(">seq1;tax=d:Archaea,g:g__Bad,s:s__Bad;\nACGT\n", encoding="utf-8")

    summary = validate_build(build)

    assert summary.failed
    assert _has_error(build, "SINTAX tax values contain rank prefixes")


def test_dada2_assignspecies_semicolon_fails_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    path = build / "export" / "dada2"
    path.mkdir(parents=True)
    (path / "bad_assignSpecies.fa").write_text(">seq1 Genus;Species\nACGT\n", encoding="utf-8")

    summary = validate_build(build)

    assert summary.failed
    assert _has_error(build, "contains semicolon")


def test_non_strict_succeeds_when_only_selected_warnings_exist(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _remove_representative(build, "s__D20s000001")

    summary = validate_build(build, strict=False, check_exports=False)

    assert not summary.failed
    assert summary.warnings > 0


def test_strict_mode_fails_selected_warnings(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _remove_representative(build, "s__D20s000001")

    summary = validate_build(build, strict=True, check_exports=False)

    assert summary.failed
    assert _has_error(build, "lacks representative")


def test_validation_reports_are_written(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)

    summary = validate_build(build, check_exports=False)

    assert summary.report_md.exists()
    assert summary.report_tsv.exists()
    content = summary.report_tsv.read_text(encoding="utf-8")
    assert "\t" in content
    assert "\n" in content
    assert "\\t" not in content
    assert "\\n" not in content


def test_missing_placement_evidence_warns_in_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    (build / "datasets" / "01_digester2020" / "placement_evidence.tsv").unlink()

    summary = validate_build(build, check_exports=False)

    assert not summary.failed
    assert _has_warning(build, "Missing placement evidence")


def test_missing_silva_resolve_evidence_warns_in_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    (build / "silva" / "silva_unresolved_evidence.tsv").unlink()

    summary = validate_build(build, check_exports=False)

    assert not summary.failed
    assert _has_warning(build, "Missing SILVA unresolved resolve evidence")


def test_strict_mode_fails_missing_silva_resolve_evidence(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    (build / "silva" / "silva_unresolved_evidence.tsv").unlink()

    summary = validate_build(build, strict=True, check_exports=False)

    assert summary.failed
    assert _has_error(build, "Missing SILVA unresolved resolve evidence")


def test_strict_mode_fails_missing_placement_evidence(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    (build / "datasets" / "01_digester2020" / "placement_evidence.tsv").unlink()

    summary = validate_build(build, strict=True, check_exports=False)

    assert summary.failed
    assert _has_error(build, "Missing placement evidence")


def test_sina_candidate_unmatched_warns_in_validation(report_tmp_dir: Path) -> None:
    build = _make_report_build(report_tmp_dir)
    _write_tsv(
        build / "datasets" / "01_digester2020" / "cluster_search_summary.tsv",
        [
            {
                "dataset": "digester2020",
                "iddef": "2",
                "sina_candidate_source": "sina.candidates.tsv",
                "sina_candidate_queries": "1",
                "sina_candidate_targets": "1",
                "sina_candidate_target_matches": "0",
            }
        ],
        [
            "dataset",
            "iddef",
            "sina_candidate_source",
            "sina_candidate_queries",
            "sina_candidate_targets",
            "sina_candidate_target_matches",
        ],
    )

    summary = validate_build(build, check_exports=False)

    assert not summary.failed
    assert _has_warning(build, "SINA candidates had no matching")


def _make_report_build(tmp_dir: Path) -> Path:
    build = tmp_dir / "build"
    registry = build / "registry"
    silva = build / "silva"
    dataset = build / "datasets" / "01_digester2020"
    registry.mkdir(parents=True)
    silva.mkdir(parents=True)
    dataset.mkdir(parents=True)
    rows = [
        _taxon("D1", "domain", "d__Archaea", "", protected="true", silva_named="true"),
        _taxon("P1", "phylum", "p__Euryarchaeota", "D1", protected="true", silva_named="true"),
        _taxon("C1", "class", "c__Methanobacteria", "P1", protected="true", silva_named="true"),
        _taxon("O1", "order", "o__Methanobacteriales", "C1", protected="true", silva_named="true"),
        _taxon("F1", "family", "f__Methanobacteriaceae", "O1", protected="true", silva_named="true"),
        _taxon("G1", "genus", "g__Methanobacterium", "F1", protected="true", silva_named="true"),
        _taxon("S1", "species", "s__Methanobacterium formicicum", "G1", protected="true", silva_named="true"),
        _taxon("g__SILVAg000001", "genus", "g__SILVAg000001", "F1", protected="true", silva_unresolved="true", is_placeholder="true", source_prefix="SILVA"),
        _taxon("s__SILVAs000001", "species", "s__SILVAs000001", "g__SILVAg000001", protected="true", silva_unresolved="true", is_placeholder="true", source_prefix="SILVA"),
        _taxon("g__D20g000001", "genus", "g__D20g000001", "F1", is_placeholder="true", created="digester2020"),
        _taxon("s__D20s000001", "species", "s__D20s000001", "g__D20g000001", is_placeholder="true", created="digester2020"),
    ]
    _write_tsv(registry / "taxon_nodes.tsv", rows, list(rows[0]))
    _write_tsv(
        registry / "dataset_registry.tsv",
        [
            {
                "dataset_name": "digester2020",
                "prefix": "D20",
                "add_order": "1",
                "input_fasta": "input.fa",
                "input_md5": "abc",
                "domain": "Archaea",
                "status": "prepared",
                "dataset_dir": str(dataset),
            }
        ],
        ["dataset_name", "prefix", "add_order", "input_fasta", "input_md5", "domain", "status", "dataset_dir"],
    )
    _write_tsv(
        registry / "name_index.tsv",
        [{"name": row["name"], "rank": row["rank"], "taxon_id": row["taxon_id"]} for row in rows],
        ["name", "rank", "taxon_id"],
    )
    (registry / "placeholder_counters.yaml").write_text(
        yaml.safe_dump({"D20": {"genus": 2, "species": 2}, "SILVA": {"genus": 2, "species": 2}}),
        encoding="utf-8",
    )
    _write_tsv(
        registry / "cluster_to_taxon.tsv",
        [
            {"cluster_key": "D20|genus|F1|0.901|D20_000003", "taxon_id": "g__D20g000001", "rank": "genus", "name": "g__D20g000001", "status": "active", "source_prefix": "D20"},
            {"cluster_key": "D20|species|g__D20g000001|0.972|D20_000003", "taxon_id": "s__D20s000001", "rank": "species", "name": "s__D20s000001", "status": "active", "source_prefix": "D20"},
        ],
        ["cluster_key", "taxon_id", "rank", "name", "status", "source_prefix"],
    )
    sequences = {
        "REF_NAMED": "ACGTACGT",
        "REF_UNRES": "TGCATGCA",
        "D20_000001": "ACGTACGT",
        "D20_000002": "TGCATGCA",
        "D20_000003": "CCCCGGGG",
    }
    _write_tsv(
        registry / "representative_registry.tsv",
        [
            {"representative_seq_id": "REF_NAMED", "taxon_id": "S1", "status": "active", "source_category": "named_silva"},
            {"representative_seq_id": "REF_UNRES", "taxon_id": "s__SILVAs000001", "status": "active", "source_category": "unresolved_silva"},
            {"representative_seq_id": "D20_000003", "taxon_id": "s__D20s000001", "status": "active", "source_category": "current_dataset", "dataset": "digester2020"},
        ],
        ["representative_seq_id", "taxon_id", "status", "source_category", "dataset"],
    )
    _write_tsv(
        registry / "sequence_registry.tsv",
        [
            {"seq_id": "REF_NAMED", "source": "SILVA138.2_NR99", "sequence_md5": sequence_md5(sequences["REF_NAMED"]), "taxon_id": "S1", "is_silva_named": "true"},
            {"seq_id": "REF_UNRES", "source": "SILVA138.2_NR99", "sequence_md5": sequence_md5(sequences["REF_UNRES"]), "taxon_id": "s__SILVAs000001", "is_silva_unresolved": "true"},
            {"seq_id": "D20_000001", "internal_seq_id": "D20_000001", "original_seq_id": "orig1", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000001"]), "taxon_id": "S1"},
            {"seq_id": "D20_000002", "internal_seq_id": "D20_000002", "original_seq_id": "orig2", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000002"]), "taxon_id": "s__SILVAs000001", "is_duplicate_sequence": "true"},
            {"seq_id": "D20_000003", "internal_seq_id": "D20_000003", "original_seq_id": "orig3", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000003"]), "taxon_id": "s__D20s000001"},
        ],
        ["seq_id", "internal_seq_id", "original_seq_id", "dataset", "prefix", "source", "sequence_md5", "taxon_id", "is_silva_named", "is_silva_unresolved", "is_duplicate_sequence"],
    )
    (registry / "current_representatives.fa").write_text(
        "".join(f">{seq_id}\n{seq}\n" for seq_id, seq in sequences.items()),
        encoding="utf-8",
    )
    _write_tsv(
        silva / "silva_unresolved_members.tsv",
        [{"seq_id": "REF_UNRES", "resolved_taxonomy": "d__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__SILVAg000001;s__SILVAs000001"}],
        ["seq_id", "resolved_taxonomy"],
    )
    _write_tsv(
        silva / "silva_unresolved_evidence.tsv",
        _resolve_evidence_rows(["REF_UNRES"]),
        ["source_stage", "dataset", "seq_id", "rank", "decision", "output_taxon"],
    )
    (dataset / "sina.oriented.fa").write_text(
        "".join(f">{seq_id}\n{sequences[seq_id]}\n" for seq_id in ["D20_000001", "D20_000002", "D20_000003"]),
        encoding="utf-8",
    )
    _write_tsv(
        dataset / "prepare_summary.tsv",
        [{"dataset": "digester2020", "prefix": "D20", "input_sequences": "3", "normalized_sequences": "3", "final_prepared_sequences": "3"}],
        ["dataset", "prefix", "input_sequences", "normalized_sequences", "final_prepared_sequences"],
    )
    _write_tsv(
        dataset / "sequence_id_map.tsv",
        [
            {"internal_seq_id": "D20_000001", "original_seq_id": "orig1", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000001"])},
            {"internal_seq_id": "D20_000002", "original_seq_id": "orig2", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000002"])},
            {"internal_seq_id": "D20_000003", "original_seq_id": "orig3", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000003"])},
        ],
        ["internal_seq_id", "original_seq_id", "dataset", "prefix", "sequence_md5"],
    )
    _write_tsv(
        dataset / "sequence_membership.tsv",
        [
            {"internal_seq_id": "D20_000001", "original_seq_id": "orig1", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000001"]), "unique_seq_id": "U_000001", "is_duplicate_sequence": "false"},
            {"internal_seq_id": "D20_000002", "original_seq_id": "orig2", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000002"]), "unique_seq_id": "U_000002", "is_duplicate_sequence": "true"},
            {"internal_seq_id": "D20_000003", "original_seq_id": "orig3", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000003"]), "unique_seq_id": "U_000003", "is_duplicate_sequence": "false"},
        ],
        ["internal_seq_id", "original_seq_id", "dataset", "prefix", "sequence_md5", "unique_seq_id", "is_duplicate_sequence"],
    )
    _write_tsv(
        dataset / "assignments.tsv",
        [
            {"internal_seq_id": "D20_000001", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000001"]), "final_status": "known_like", "assigned_taxon_id": "S1", "best_hit_source_category": "named_silva"},
            {"internal_seq_id": "D20_000002", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000002"]), "final_status": "known_like", "assigned_taxon_id": "s__SILVAs000001", "best_hit_source_category": "unresolved_silva"},
            {"internal_seq_id": "D20_000003", "dataset": "digester2020", "prefix": "D20", "sequence_md5": sequence_md5(sequences["D20_000003"]), "final_status": "new_species", "assigned_taxon_id": "s__D20s000001", "best_hit_source_category": "named_silva"},
        ],
        ["internal_seq_id", "dataset", "prefix", "sequence_md5", "final_status", "assigned_taxon_id", "best_hit_source_category"],
    )
    _write_tsv(
        dataset / "placement_summary.tsv",
        [{"dataset": "digester2020", "prefix": "D20", "assigned_named_silva": "1", "assigned_unresolved_silva": "1", "assigned_previous_custom": "0", "new_current_dataset": "1", "known_like": "2", "new_species": "1", "created_species": "1", "created_genera": "1"}],
        ["dataset", "prefix", "assigned_named_silva", "assigned_unresolved_silva", "assigned_previous_custom", "new_current_dataset", "known_like", "new_species", "created_species", "created_genera"],
    )
    _write_tsv(
        dataset / "created_taxa.tsv",
        [
            {"taxon_id": "g__D20g000001", "rank": "genus", "name": "g__D20g000001"},
            {"taxon_id": "s__D20s000001", "rank": "species", "name": "s__D20s000001"},
        ],
        ["taxon_id", "rank", "name"],
    )
    _write_tsv(
        dataset / "representative_updates.tsv",
        [{"representative_seq_id": "D20_000003", "taxon_id": "s__D20s000001", "dataset": "digester2020", "action": "add", "reason": "new_species"}],
        ["representative_seq_id", "taxon_id", "dataset", "action", "reason"],
    )
    _write_tsv(
        dataset / "placement_evidence.tsv",
        _placement_evidence_rows(["D20_000001", "D20_000002", "D20_000003"]),
        ["source_stage", "dataset", "seq_id", "rank", "decision", "output_taxon"],
    )
    _write_tsv(
        dataset / "sina.candidates.tsv",
        [
            {"query": "D20_000001", "target": "REF_NAMED", "sina_score": "0.98", "source_field": "nearest_slv", "rank": "1"},
            {"query": "D20_000002", "target": "REF_UNRES", "sina_score": "0.97", "source_field": "nearest_slv", "rank": "1"},
        ],
        ["query", "target", "sina_score", "source_field", "rank"],
    )
    _write_tsv(
        dataset / "cluster_search_summary.tsv",
        [
            {
                "dataset": "digester2020",
                "iddef": "2",
                "sina_candidate_source": "sina.candidates.tsv",
                "sina_candidate_queries": "2",
                "sina_candidate_targets": "2",
                "sina_candidate_target_matches": "2",
            }
        ],
        ["dataset", "iddef", "sina_candidate_source", "sina_candidate_queries", "sina_candidate_targets", "sina_candidate_target_matches"],
    )
    return build


def _taxon(
    taxon_id: str,
    rank: str,
    name: str,
    parent: str,
    status: str = "active",
    protected: str = "false",
    silva_named: str = "false",
    silva_unresolved: str = "false",
    is_placeholder: str = "false",
    source_prefix: str = "",
    created: str = "",
) -> dict[str, str]:
    return {
        "taxon_id": taxon_id,
        "rank": rank,
        "name": name,
        "parent_taxon_id": parent,
        "status": status,
        "protected": protected,
        "is_silva_named": silva_named,
        "is_silva_unresolved": silva_unresolved,
        "is_placeholder": is_placeholder,
        "source_prefix": source_prefix,
        "created_in_dataset": created,
        "source": "custom_dataset" if created else "SILVA138.2_NR99",
    }


def _remove_representative(build: Path, taxon_id: str) -> None:
    path = build / "registry" / "representative_registry.tsv"
    rows = [row for row in _read_tsv(path) if row["taxon_id"] != taxon_id]
    _write_tsv(path, rows, ["representative_seq_id", "taxon_id", "status", "source_category", "dataset"])


def _has_error(build: Path, text: str) -> bool:
    rows = _read_tsv(build / "reports" / "validation_report.tsv")
    return any(row["level"] == "error" and text.lower() in row["message"].lower() for row in rows)


def _has_warning(build: Path, text: str) -> bool:
    rows = _read_tsv(build / "reports" / "validation_report.tsv")
    return any(row["level"] == "warning" and text.lower() in row["message"].lower() for row in rows)


def _placement_evidence_rows(seq_ids: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for seq_id in seq_ids:
        for rank in ["phylum", "class", "order", "family", "genus", "species"]:
            rows.append(
                {
                    "source_stage": "dataset_place",
                    "dataset": "digester2020",
                    "seq_id": seq_id,
                    "rank": rank,
                    "decision": "assign_existing",
                    "output_taxon": rank,
                }
            )
    return rows


def _resolve_evidence_rows(seq_ids: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for seq_id in seq_ids:
        for rank in ["phylum", "class", "order", "family", "genus", "species"]:
            rows.append(
                {
                    "source_stage": "silva_resolve",
                    "dataset": "SILVA138.2_NR99",
                    "seq_id": seq_id,
                    "rank": rank,
                    "decision": "assign_existing",
                    "output_taxon": rank,
                }
            )
    return rows


def _append_tsv(path: Path, row: dict[str, str]) -> None:
    rows = _read_tsv(path)
    fieldnames = list(rows[0])
    rows.append(row)
    _write_tsv(path, rows, fieldnames)


def _write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
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


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=chr(9)))
