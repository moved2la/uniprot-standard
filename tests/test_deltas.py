"""isoform_processing_deltas: arithmetic on residue fractions."""
from pipeline import isoform_processing_deltas as ipd
from pipeline import common


def test_fractions_sum_to_one_and_count():
    f = ipd.fractions("GGGAAAK")
    assert abs(sum(f.values()) - 1.0) < 1e-12
    assert abs(f["G"] - 3 / 7) < 1e-12 and f["W"] == 0.0


def test_delta_and_bound():
    a, b = ipd.fractions("GGGAAAK"), ipd.fractions("GGGAAAKK")
    d = ipd.delta(b, a)
    m, aa = ipd.max_abs(d)
    assert aa == "K" and abs(m - (2 / 8 - 1 / 7)) < 1e-12
    assert abs(sum(d.values())) < 1e-12  # fractions of both sum to one


def test_master_sequence_uses_true_segments_only():
    cp = common.new_ini()
    cp["X.n_terminal_removed"] = {"start": "1", "end": "1", "in_master_molecule": "false"}
    cp["X.master"] = {"start": "2", "end": "7", "in_master_molecule": "true"}
    cp["Y.master"] = {"start": "1", "end": "3", "in_master_molecule": "true"}
    assert ipd.master_sequence("X", "MGGGAAAK"[:7], cp) == "GGGAAA"
    assert ipd.master_sequence("Z", "MGG", cp) == ""
