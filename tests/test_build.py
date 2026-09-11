"""build_protein_set: the Chain rule, segment generation, and flag behaviour."""
import configparser
import pytest

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
    rule, s, e, note, flag = bps.chain_rule("X", sec, _dec())
    assert (rule, s, e, flag) == ("R2a", 2, 13, None)
    segs = bps.segments_for("X", sec, s, e, rule, note)
    names = [n for n, _ in segs]
    assert names == ["X.n_terminal_removed", "X.master"]
    assert segs[0][1]["in_master_molecule"] == "false" and "Initiator methionine 1-1" in segs[0][1]["note"]
    assert segs[1][1]["start"] == "2" and segs[1][1]["end"] == "13"


def test_no_chain_uses_whole_sequence():
    sec = _sec(10, [])
    rule, s, e, note, flag = bps.chain_rule("X", sec, _dec())
    assert (rule, s, e, flag) == ("R2b", 1, 10, None)
    assert [n for n, _ in bps.segments_for("X", sec, s, e, rule, note)] == ["X.master"]


def test_multiple_chains_flag_until_closed():
    sec = _sec(30, [("Chain", 2, 15, "a", "PRO_1"), ("Chain", 16, 30, "b", "PRO_2")])
    rule, s, e, note, flag = bps.chain_rule("X", sec, _dec())
    assert rule == "R2c" and flag["rule"] == "multiple_chain_features"
    closed = _dec(**{"flag.X": {"rule": "multiple_chain_features", "decision": "D99",
                                "master_start": "2", "master_end": "30"}})
    rule, s, e, note, flag = bps.chain_rule("X", sec, closed)
    assert (rule, s, e) == ("R2c-closed", 2, 30) and "D99" in note
    assert flag["rule"] == "multiple_chain_features"  # still recorded in flags.tsv, just closed


def test_uncertain_chain_position_flags():
    sec = _sec(30, [("Chain", "2(UNSURE)", 30, "a", "PRO_1")])
    rule, s, e, note, flag = bps.chain_rule("X", sec, _dec())
    assert rule == "R2c" and flag["rule"] == "chain_position_uncertain"


def test_c_terminal_removed_segment():
    sec = _sec(20, [("Chain", 1, 15, "p", "PRO_1")], [("Propeptide", 16, 20, "Removed", "PRO_2")])
    rule, s, e, note, flag = bps.chain_rule("X", sec, _dec())
    segs = dict(bps.segments_for("X", sec, s, e, rule, note))
    assert set(segs) == {"X.master", "X.c_terminal_removed"}
    assert "Propeptide 16-20" in segs["X.c_terminal_removed"]["note"]
