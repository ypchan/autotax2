from src.assign import parse_vsearch_uc_for_new_clusters, uc_assignments_from_old_hits


def test_old_assignments(tmp_path):
    p = tmp_path / "hits.tsv"
    p.write_text("q1\tc1\t97.5\t100\t0\t0\t1\t100\t1\t100\t1.0\t1.0\n")
    rows = uc_assignments_from_old_hits(p)
    assert rows[0]["query"] == "q1"
    assert rows[0]["cluster_id"] == "c1"
    assert rows[0]["assignment_type"] == "old_cluster"


def test_new_uc_parse(tmp_path):
    p = tmp_path / "new.uc"
    p.write_text(
        "S\t0\t100\t*\t*\t*\t*\t*\tseqA\t*\n"
        "H\t0\t100\t98.0\t+\t0\t0\t100M\tseqB\tseqA\n"
    )
    rows = parse_vsearch_uc_for_new_clusters(p, {"seqA": "new_1"})
    assert rows[0]["cluster_id"] == "new_1"
    assert rows[1]["cluster_id"] == "new_1"
