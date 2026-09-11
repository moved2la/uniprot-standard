"""resolve_ontology_term: parsing, exact-name matching, subtree and ancestor walks."""
import pytest
from conftest import FIXTURES

from pipeline import resolve_ontology_term as rot


@pytest.fixture(scope="module")
def onto():
    return rot.parse_obo((FIXTURES / "ontology.obo").read_text())


def test_header_and_term_count(onto):
    assert onto.header["data-version"] == "releases/2099-01-01"
    assert onto.header["format-version"] == "1.2"
    assert len(onto.terms) == 10


def test_exact_match_is_exact(onto):
    assert [t.id for t in rot.exact_name_matches(onto, "widget")] == ["TEST:0000002"]
    assert rot.exact_name_matches(onto, "Widget") == []          # case matters
    assert rot.exact_name_matches(onto, "widget ") == []         # whitespace matters
    assert rot.exact_name_matches(onto, "widge") == []           # no substring search


def test_duplicate_and_obsolete_names_are_visible(onto):
    dup = rot.exact_name_matches(onto, "duplicate name")
    assert {t.id for t in dup} == {"TEST:0000007", "TEST:0000008"}
    ret = rot.exact_name_matches(onto, "retired thing")
    assert len(ret) == 1 and ret[0].obsolete


def test_subtree_follows_only_is_a_and_part_of(onto):
    sub = rot.descendants(onto, "TEST:0000002")
    ids = [tid for tid, *_ in sub]
    assert ids[0] == "TEST:0000002"
    assert set(ids) == {"TEST:0000002", "TEST:0000003", "TEST:0000004", "TEST:0000005"}
    assert "TEST:0000010" not in ids  # linked by 'regulates', not followed
    rels = {tid: rel for tid, _, rel, _ in sub}
    assert rels["TEST:0000003"] == "part_of" and rels["TEST:0000004"] == "is_a"


def test_has_part_is_not_followed(onto):
    # widget shell has_part widget core; that must not make core a child of shell
    sub = rot.descendants(onto, "TEST:0000005")
    assert [tid for tid, *_ in sub] == ["TEST:0000005"]


def test_ancestors_and_neighbours(onto):
    anc = rot.ancestors(onto, "TEST:0000004")
    assert [tid for tid, *_ in anc] == ["TEST:0000003", "TEST:0000002", "TEST:0000001"]
    nb = rot.neighbours(onto, "TEST:0000002")
    assert ("broader", "is_a", "TEST:0000001") in nb
    assert ("narrower", "part_of", "TEST:0000003") in nb
    assert ("narrower", "part_of", "TEST:0000005") in nb
    assert not any(tid == "TEST:0000010" for _, _, tid in nb)


def test_resolve_stops_on_zero_or_many(onto, tmp_path):
    import configparser, logging
    log = logging.getLogger("t")
    for name in ("no such term", "duplicate name", "retired thing"):
        cp = configparser.ConfigParser(); cp.optionxform = str
        cp["tier.1"] = {"seed_word": "x", "lookup_name": name}
        with pytest.raises(SystemExit):
            rot.resolve(onto, cp, log, tmp_path, tmp_path / "flags.tsv")
        assert (tmp_path / "flags.tsv").exists()


def test_resolve_writes_tables(onto, tmp_path):
    import configparser, logging
    cp = configparser.ConfigParser(); cp.optionxform = str
    cp["tier.1"] = {"seed_word": "x", "lookup_name": "widget"}
    out = rot.resolve(onto, cp, logging.getLogger("t"), tmp_path, tmp_path / "flags.tsv")
    assert out == {"1": "TEST:0000002"}
    assert (tmp_path / "tier1_subtree.tsv").read_text().count("\n") == 5  # header + 4 terms
    assert not (tmp_path / "flags.tsv").exists()
