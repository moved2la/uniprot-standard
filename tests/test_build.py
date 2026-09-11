"""build_protein_set: the Chain rule (R2a/b/c/d), segment generation, and flag behaviour."""
import configparser

from pipeline import build_protein_set as bps
from pipeline import common


def _sec(length, chains, extra=()):
    d = {"length": str(length)}
    feats = list(chains) + list(extra)
    for i, (t, s, e, desc, fid) in enumerate(feats):
        d[f"feature.{i}.type"] = t
        d[f"feature.{i}.start"] = str(s)
        d[f"feature.{i}.end"] = str(e)
        d[f"feature.{i}.description"] = desc
        d[f"feature.{i}.id"] = fid
    cp = common.new_ini(); cp["X"] = d
    return cp["X"]


def _dec(**closures):
    cp = configparser.ConfigParser(); cp.optionxform = str
    for k, v in closures.items():
        cp[k] = v
    return cp


def test_one_chain_gives_master_and_outside_segments():
    sec = _sec(13, [("Chain", 2, 13, "p", "PRO_1")], [("Initiator methionine", 1, 1, "Removed", "")])
    rule, ranges, note, flag, chains = bps.chain_rule("X", sec, _dec())
    assert (rule, ranges, flag) == ("R2a", [(2, 13)], None)
    segs = bps.segments_for("X", sec, ranges, rule, note)
    assert [n for n, _ in segs] == ["X.n_terminal_removed", "X.master"]
    assert segs[0][1]["in_master_molecule"] == "false" and "Initiator methionine 1-1" in segs[0][1]["note"]
    assert segs[1][1]["start"] == "2" and segs[1][1]["end"] == "13"


def test_no_chain_uses_whole_sequence():
    sec = _sec(10, [])
    rule, ranges, note, flag, _ = bps.chain_rule("X", sec, _dec())
    assert (rule, ranges, flag) == ("R2b", [(1, 10)], None)
    assert [n for n, _ in bps.segments_for("X", sec, ranges, rule, note)] == ["X.master"]


def test_multiple_chains_union_overlapping():
    # precursor 2-30 recorded as intermediate 2-30 plus final 3-30 -> union 2-30, no flag
    sec = _sec(30, [("Chain", 2, 30, "intermediate", "PRO_1"), ("Chain", 3, 30, "final", "PRO_2")])
    rule, ranges, note, flag, chains = bps.chain_rule("X", sec, _dec())
    assert (rule, ranges, flag) == ("R2d", [(2, 30)], None) and len(chains) == 2 and "D24" in note
    assert [n for n, _ in bps.segments_for("X", sec, ranges, rule, note)] == ["X.n_terminal_removed", "X.master"]


def test_multiple_chains_union_adjacent_split_precursor():
    sec = _sec(40, [("Chain", 5, 20, "alpha", "PRO_1"), ("Chain", 21, 40, "beta", "PRO_2")])
    rule, ranges, *_ = bps.chain_rule("X", sec, _dec())
    assert (rule, ranges) == ("R2d", [(5, 40)])


def test_multiple_chains_union_with_internal_gap():
    sec = _sec(40, [("Chain", 5, 20, "a", "PRO_1"), ("Chain", 26, 38, "b", "PRO_2")],
               [("Propeptide", 21, 25, "linker", "PRO_3")])
    rule, ranges, note, flag, _ = bps.chain_rule("X", sec, _dec())
    assert (rule, ranges) == ("R2d", [(5, 20), (26, 38)])
    segs = dict(bps.segments_for("X", sec, ranges, rule, note))
    assert list(segs) == ["X.n_terminal_removed", "X.master", "X.internal_removed", "X.master_2", "X.c_terminal_removed"]
    assert "Propeptide 21-25" in segs["X.internal_removed"]["note"]
    assert segs["X.master_2"]["in_master_molecule"] == "true"


def test_uncertain_chain_position_flags_until_closed():
    sec = _sec(30, [("Chain", "2(UNSURE)", 30, "a", "PRO_1")])
    rule, ranges, note, flag, _ = bps.chain_rule("X", sec, _dec())
    assert rule == "R2c" and ranges == [] and flag["rule"] == "chain_position_uncertain"
    closed = _dec(**{"flag.X": {"rule": "chain_position_uncertain", "decision": "D99",
                                "master_start": "2", "master_end": "30"}})
    rule, ranges, note, flag, _ = bps.chain_rule("X", sec, closed)
    assert (rule, ranges) == ("R2c-closed", [(2, 30)]) and "D99" in note
    assert flag["rule"] == "chain_position_uncertain"  # still recorded in flags.tsv, just closed


def test_c_terminal_removed_segment():
    sec = _sec(20, [("Chain", 1, 15, "p", "PRO_1")], [("Propeptide", 16, 20, "Removed", "PRO_2")])
    rule, ranges, note, flag, _ = bps.chain_rule("X", sec, _dec())
    segs = dict(bps.segments_for("X", sec, ranges, rule, note))
    assert set(segs) == {"X.master", "X.c_terminal_removed"}
    assert "Propeptide 16-20" in segs["X.c_terminal_removed"]["note"]


def test_union_helper():
    assert bps._union([(3, 5), (1, 2), (10, 12), (6, 6)]) == [(1, 6), (10, 12)]
