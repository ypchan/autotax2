from autotax2.sina import parse_sina_header, run_sina


def test_parse_sina_header_with_identity_cutoff_lca_and_turn():
    header = (
        "ERR14064717_s_ERR14064717.1773379 "
        "[align_cutoff_head_slv=701] "
        "[align_cutoff_tail_slv=5] "
        "[align_ident_slv=60.1063843] "
        "[align_quality_slv=66] "
        "[lca_tax_gtdb=Unclassified;] "
        "[lca_tax_slv=Unclassified;] "
        "[turn=reversed and complemented]"
    )

    ann = parse_sina_header(header)

    assert ann.seq_id == "ERR14064717_s_ERR14064717.1773379"
    assert ann.align_cutoff_head == 701
    assert ann.align_cutoff_tail == 5
    assert ann.align_ident == 60.1063843
    assert ann.align_quality == 66
    assert ann.lca_tax_ref == "Unclassified;"
    assert ann.turn == "reversed and complemented"


def test_run_sina_requests_header_metadata(monkeypatch, tmp_path):
    captured = {}

    def fake_run_command(cmd, log_path=None, dry_run=False):
        captured["cmd"] = cmd
        captured["log_path"] = log_path
        captured["dry_run"] = dry_run

    monkeypatch.setattr("autotax2.sina.run_command", fake_run_command)

    out = tmp_path / "out.fa"
    run_sina(
        input_fa=tmp_path / "in.fa",
        ref_arb=tmp_path / "gtdb_ssu.arb",
        output_fa=out,
        threads=48,
        dry_run=True,
    )

    cmd = captured["cmd"]
    assert "--db" in cmd
    assert "--search" in cmd
    assert "--search-min-sim" in cmd
    assert "--search-max-result" in cmd
    assert "--lca-fields=tax_slv,tax_gtdb" in cmd
    assert "--lca-quorum" in cmd
    assert "--show-conf" in cmd
    assert "--fasta-write-dna" in cmd
    assert "--turn" in cmd
    assert "all" in cmd
    assert "--calc-idty" in cmd
    assert "--meta-fmt" in cmd
    assert "header" in cmd
    assert str(out) in cmd
