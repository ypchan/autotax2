from pathlib import Path
from autotax2.vsearch import uc_membership


def test_uc_membership(tmp_path: Path):
    uc = tmp_path / "x.uc"
    uc.write_text("S\t0\t100\t*\t*\t*\t*\t*\tseq1\t*\nH\t0\t100\t99.0\t+\t0\t0\t100M\tseq2\tseq1\n")
    mem = uc_membership(uc, "species")
    assert len(mem) == 2
    assert set(mem["centroid"]) == {"seq1"}
