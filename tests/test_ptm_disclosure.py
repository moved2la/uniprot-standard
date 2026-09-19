"""fetch_ptmlist.parse_ptmlist and ptm_disclosure.price_feature, offline, synthetic fixtures.

What these tests check, in plain terms:
- the vocabulary parser splits the flat file into records keyed by ID and keeps
  the two-letter-coded lines, ignoring the file's own header;
- a feature whose description's first clause matches a vocabulary record with an
  MA line gets that mass; a record without MA is counted only; an unmatched
  description is counted only; a disulfide bond within the chain is priced as
  minus H2 and an interchain one is counted only;
- a vocabulary ID containing "..." matches as a wildcard pattern, a leading
  "(Microbial infection) " qualifier is set aside and recorded, and the rule that
  matched plus the vocabulary's correction formula are recorded per site;
- a site outside the master range is marked in_master = false.
"""
import json
import pytest
from conftest import FIXTURES

from pipeline import ptm_disclosure as ptm
from pipeline.fetch_ptmlist import parse_ptmlist


@pytest.fixture
def vocab():
    return parse_ptmlist((FIXTURES / "ptmlist_synthetic.txt").read_text())


@pytest.fixture
def entry():
    return json.loads((FIXTURES / "entry_ptm.json").read_text())


def test_parse_ptmlist(vocab):
    assert set(vocab) == {"Fixturyl-lysine", "Fixturyl-serine no mass", "Fixture-linked (Xyz...) asparagine",
                          "Fixtural lysine isopeptide (Lys-Fix) (interchain with F-...)"}
    assert vocab["Fixturyl-lysine"]["MA"] == ["14.03"] and vocab["Fixturyl-lysine"]["AC"] == ["PTM-9001"]
    assert "MA" not in vocab["Fixturyl-serine no mass"]


def test_parse_ptmlist_rejects_duplicate_id():
    text = (FIXTURES / "ptmlist_synthetic.txt").read_text()
    with pytest.raises(ValueError):
        parse_ptmlist(text + text)


def test_lookup_key_takes_first_clause_and_sets_qualifier_aside():
    assert ptm.lookup_key("Fixturyl-lysine; by FIXK") == ("Fixturyl-lysine", "")
    assert ptm.lookup_key("Plain") == ("Plain", "")
    assert ptm.lookup_key("(Microbial infection) Fixturyl-lysine; by X") == ("Fixturyl-lysine", "(Microbial infection)")


def test_find_record_exact_then_wildcard(vocab):
    pats = ptm.wildcard_patterns(vocab)
    assert [v for _, v in pats] == ["Fixture-linked (Xyz...) asparagine",
                                    "Fixtural lysine isopeptide (Lys-Fix) (interchain with F-...)"]
    rec, vid, rule = ptm.find_record("Fixturyl-lysine", vocab, pats)
    assert rule == "exact" and vid == "Fixturyl-lysine" and rec["AC"] == ["PTM-9001"]
    rec, vid, rule = ptm.find_record("Fixtural lysine isopeptide (Lys-Fix) (interchain with F-Cter in FIXUB)", vocab, pats)
    assert rule == "wildcard" and vid.endswith("(interchain with F-...)") and rec["MA"] == ["-18.02"]
    rec, vid, rule = ptm.find_record("Nothing like this", vocab, pats)
    assert rec is None and rule == "none"
    # the wildcard must not match a different linkage
    rec, vid, rule = ptm.find_record("Fixtural lysine isopeptide (Lys-Fix) (interchain with Q-Cter in X)", vocab, pats)
    assert rec is None and rule == "none"


def test_price_feature_rules(vocab, entry):
    master = [(2, 8)]
    h2 = 2.0
    feats = [f for f in entry["features"] if f["type"] in ptm.PTM_FEATURE_TYPES]
    priced = [ptm.price_feature(f, vocab, h2, master) for f in feats]
    by_desc = {p["description"]: p for p in priced}

    p = by_desc["Fixturyl-lysine; by FIXK"]
    assert p["mass_status"] == "vocabulary" and p["mass_delta"] == pytest.approx(14.03) and p["in_master"]
    assert p["vocabulary_ac"] == "PTM-9001"
    assert by_desc["Fixturyl-serine no mass"]["mass_status"] == "no_mass_in_vocabulary"
    assert by_desc["Fixturyl-serine no mass"]["mass_delta"] is None
    assert by_desc["Fixture-linked (Xyz...) asparagine"]["mass_status"] == "no_mass_in_vocabulary"
    assert by_desc["Unknown fixture mod"]["mass_status"] == "no_vocabulary_match"
    # site at position 1 is outside the master range 2-8
    outside = [p for p in priced if p["start"] == 1][0]
    assert outside["in_master"] is False and outside["mass_status"] == "vocabulary"
    # disulfides
    intra = [p for p in priced if p["feature_type"] == "Disulfide bond" and p["description"] == ""][0]
    assert intra["mass_status"] == "disulfide_arithmetic" and intra["mass_delta"] == pytest.approx(-2.0)
    inter = by_desc["Interchain (with fixture partner)"]
    assert inter["mass_status"] == "interchain_count_only" and inter["mass_delta"] is None
    # wildcard cross-link: priced from the vocabulary's MA, rule and CF recorded
    xl = [p for p in priced if p["feature_type"] == "Cross-link"][0]
    assert xl["match_rule"] == "wildcard" and xl["mass_delta"] == pytest.approx(-18.02)
    assert xl["correction_formula"] == "H-2 O-1" and xl["vocabulary_ac"] == "PTM-9004"
    # qualifier set aside, still priced exactly
    q = [p for p in priced if p["qualifier"]][0]
    assert q["qualifier"] == "(Microbial infection)" and q["lookup_key"] == "Fixturyl-lysine"
    assert q["match_rule"] == "exact" and q["mass_delta"] == pytest.approx(14.03)
    # the Region feature was not a PTM type and was never priced
    assert len(priced) == 9


def test_alternatives_at_one_position_contribute_once():
    """C3b (D118): priced features of one type at one position are alternatives -- the largest
    |delta| counts once; different positions or types still add."""
    from pipeline.ptm_disclosure import collapse_alternatives
    priced = [("Glycosylation", 180, 180, 1216.4), ("Glycosylation", 180, 180, 1622.6), ("Glycosylation", 180, 180, 203.2),
              ("Glycosylation", 300, 300, 203.2), ("Modified residue", 180, 180, 42.0), ("Disulfide bond", 22, 96, -2.016)]
    total, collapsed = collapse_alternatives(priced)
    assert collapsed == 2
    assert abs(total - (1622.6 + 203.2 + 42.0 - 2.016)) < 1e-9
    assert collapse_alternatives([]) == (0.0, 0)
