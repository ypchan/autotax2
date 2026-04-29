from src.provenance import pairwise_source_overlap, read_source_map, summarize_cluster_sources


def test_read_source_map(tmp_path):
    p = tmp_path / "source.tsv"
    p.write_text("sequence_id\tsource\nseq1\tA\nseq2\tB\n")
    m = read_source_map(p)
    assert m["seq1"] == "A"
    assert m["seq2"] == "B"


def test_cluster_source_summary_and_overlap():
    seq_to_cluster = {"s1": "c1", "s2": "c1", "s3": "c2"}
    source_map = {"s1": "A", "s2": "B", "s3": "A"}
    rows, members, counts = summarize_cluster_sources(seq_to_cluster, source_map, "0.97")
    c1 = [r for r in rows if r["cluster_id"] == "c1"][0]
    assert c1["unique_sources"] == "2"
    pair = pairwise_source_overlap(counts, "0.97")[0]
    assert pair["shared_cluster_count"] == "1"
