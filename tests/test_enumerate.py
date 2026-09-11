"""enumerate_pool: query construction, union-vs-term-level check, neighbour sensitivity.
Uses a fake client; no network."""
import logging

from pipeline import common
from pipeline import enumerate_pool as ep


class FakeClient:
    """Maps a query string to rows. Term-level queries expand descendants except
    where the fixture deliberately omits one, so the check has something to catch."""
    def __init__(self, table):
        self.table = table
        self.release = "2099_01"

    def search_tsv(self, query, fields):
        return [{"Entry": a, "Gene Names (primary)": g, "Protein names": "p", "Length": "1",
                 "Reviewed": "reviewed", "Organism (ID)": "9606"} for a, g in self.table.get(query, [])]


def test_term_query_syntax():
    assert ep.term_query("TEST:0000002") == "(organism_id:9606) AND (reviewed:true) AND (go:0000002)"


def test_enumerate_tier_records_union_check_and_neighbours(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "DATA_DIR", tmp_path)
    monkeypatch.setattr(common, "FLAGS_TSV", tmp_path / "flags.tsv")
    q = ep.term_query
    table = {
        q("TEST:0000002"): [("A1", "ga"), ("A2", "gb")],           # term-level: misses A3
        q("TEST:0000003"): [("A2", "gb"), ("A3", "gc")],           # subtree term
        q("TEST:0000001"): [("A1", "ga"), ("A2", "gb"), ("A3", "gc"), ("B9", "gz")],  # broader
    }
    subtree = [{"term_id": "TEST:0000002", "name": "widget", "relation_to_parent": "", "depth": "0"},
               {"term_id": "TEST:0000003", "name": "widget core", "relation_to_parent": "part_of", "depth": "1"}]
    neighbours = [{"direction": "broader", "relation": "is_a", "term_id": "TEST:0000001", "name": "root compartment"},
                  {"direction": "narrower", "relation": "part_of", "term_id": "TEST:0000003", "name": "widget core"}]
    pool, n_top, ok = ep.enumerate_tier(FakeClient(table), "1", "TEST:0000002", "widget",
                                        subtree, neighbours, logging.getLogger("t"))
    assert pool == {"A1", "A2", "A3"} and n_top == 2 and ok is False
    rows = common.read_tsv(tmp_path / "tier1_pool.tsv")
    assert [r["accession"] for r in rows] == ["A1", "A2", "A3"]
    assert dict((r["accession"], r["in_term_level_query"]) for r in rows) == {"A1": "true", "A2": "true", "A3": "false"}
    assert [r for r in rows if r["accession"] == "A2"][0]["subtree_terms_returning_entry"] == "TEST:0000002;TEST:0000003"
    sens = common.read_tsv(tmp_path / "tier1_pool_sensitivity.tsv")
    broader = [r for r in sens if r["query_label"] == "broader"][0]
    assert broader["n_not_in_pool"] == "1" and broader["accessions_not_in_pool"] == "B9"
    flags = common.read_tsv(tmp_path / "flags.tsv")
    assert flags[0]["rule"] == "descendant_expansion_check"
