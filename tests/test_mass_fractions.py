"""mass_fractions.py rules on synthetic inputs. No real biology; no network; no files
outside tmp_path. Gene symbols and accessions here are invented."""
import logging
import math

from pipeline import mass_fractions as mf

LOG = logging.getLogger("t")


def _row(idx, gene, v=(10.0, 20.0, 30.0), sd=(1.0, 2.0, 3.0), valid=(5, 5, 5)):
    r = {"row_index": idx, "gene_cell": gene}
    for ft, a, b, c in zip(mf.FIBER_TYPES, v, sd, valid):
        r[f"median_{ft}"], r[f"sd_{ft}"], r[f"valid_{ft}"] = a, b, c
    return r


def _pool(*items):
    return {acc: {"accession": acc, "gene": g, "tier": t, "tiers": t.split(";")} for acc, g, t in items}


def test_single_row_matches_by_exact_symbol():
    pool = _pool(("X00001", "GA", "1"), ("X00002", "GB", "1"))
    j = mf.join([_row(4, "GA"), _row(5, "GZ")], pool, LOG)
    assert j["measured"]["X00001"]["match_rule"] == "single"
    assert j["measured"]["X00002"]["match_rule"] == "none"     # in the pool, no row -> w = 0, listed
    assert j["measured"]["X00001"]["v_I"] == 10.0
    assert [r["gene"] for r in j["outside"]] == ["GZ"]           # the unmatched row lands in the completeness list


def test_two_rows_one_gene_are_summed_with_minor_share():
    pool = _pool(("X00001", "GA", "1"))
    j = mf.join([_row(4, "GA", v=(90.0, 90.0, 90.0)), _row(5, "GA", v=(10.0, 10.0, 10.0))], pool, LOG)
    m = j["measured"]["X00001"]
    assert m["match_rule"] == "summed:2" and m["v_I"] == 100.0
    assert math.isclose(m["minor_share_I"], 0.1)
    assert math.isclose(m["sd_I"], math.sqrt(1.0 ** 2 + 1.0 ** 2))


def test_nan_with_valid_zero_is_zero_and_nan_with_valid_positive_stops():
    pool = _pool(("X00001", "GA", "1"), ("X00002", "GB", "1"))
    ok = _row(4, "GA", v=(float("nan"), 5.0, 5.0), valid=(0, 5, 5))
    j = mf.join([ok], pool, LOG)
    assert j["measured"]["X00001"]["v_I"] == 0.0 and not j["stops"]
    bad = _row(5, "GB", v=(float("nan"), 5.0, 5.0), valid=(3, 5, 5))
    assert mf.join([bad], pool, LOG)["stops"]


def test_shared_cell_one_member_goes_to_it_two_members_split_by_own_rows():
    pool = _pool(("X00001", "GA", "1"), ("X00002", "GB", "1"), ("X00003", "GC", "1"))
    rows = [_row(4, "GA", v=(30.0, 30.0, 30.0)), _row(5, "GB", v=(10.0, 10.0, 10.0)),
            _row(6, "NOTINPOOL;GC", v=(8.0, 8.0, 8.0)),          # one member in pool
            _row(7, "GA;GB", v=(40.0, 40.0, 40.0))]              # two members: 3:1 by own rows
    j = mf.join(rows, pool, LOG)
    m = j["measured"]
    assert m["X00003"]["match_rule"] == "shared_group" and m["X00003"]["v_I"] == 8.0
    assert math.isclose(m["X00001"]["v_I"], 30.0 + 30.0) and math.isclose(m["X00002"]["v_I"], 10.0 + 10.0)
    assert m["X00001"]["shared_whole_I"] == 40.0 and m["X00002"]["shared_whole_I"] == 40.0   # whole row -> w_high
    assert m["X00001"]["match_rule"] == "single+shared_split"


def test_weights_sum_to_one_per_tier_and_low_high_use_fixed_denominator():
    pool = _pool(("X00001", "GA", "1"), ("X00002", "GB", "1;2"), ("X00003", "GC", "2"))
    rows = [_row(4, "GA", v=(10.0, 10.0, 10.0), sd=(5.0, 5.0, 5.0)),
            _row(5, "GB", v=(10.0, 10.0, 10.0), sd=(0.0, 0.0, 0.0)),
            _row(6, "GC", v=(0.0, 0.0, 0.0))]
    j = mf.join(rows, pool, LOG)
    mw = {"X00001": 1000.0, "X00002": 3000.0, "X00003": 500.0}
    mf.weights(j["measured"], mw, LOG)
    m = j["measured"]
    for t in ("1", "2"):
        for ft in mf.FIBER_TYPES:
            s = sum(x[f"w_{ft}_tier{t}"] or 0.0 for x in m.values())
            assert math.isclose(s, 1.0)
    assert math.isclose(m["X00001"]["w_I_tier1"], 10000 / 40000)
    assert math.isclose(m["X00001"]["w_I_tier1_low"], 5000 / 40000)     # (v - sd) * mw / same denominator
    assert math.isclose(m["X00001"]["w_I_tier1_high"], 15000 / 40000)
    assert m["X00001"]["w_I_tier2"] is None                               # not in tier 2
    assert m["X00002"]["w_I_tier2"] == 1.0                                # GC has zero: GB is all of tier 2


def test_families_union_find_respects_cutoff(tmp_path):
    p = tmp_path / "pairs.tsv"
    p.write_text("# comment\naccession_a\tgene_a\taccession_b\tgene_b\tn_shared_peptides\tfrac_of_a\tfrac_of_b\n"
                 "X00001\tGA\tX00002\tGB\t5\t0.1\t0.1\nX00002\tGB\tX00003\tGC\t1\t0.1\t0.1\n", encoding="utf-8")
    f1 = mf.families_from_pairs(1, p)
    f2 = mf.families_from_pairs(2, p)
    assert f1["X00001"] == {"X00001", "X00002", "X00003"}
    assert f2["X00001"] == {"X00001", "X00002"} and "X00003" not in f2


def test_accession_cell_matches_pool_with_isoform_suffix_dropped():
    pool = _pool(("X00001", "GA", "1"), ("X00002", "GB", "1"))
    assert mf.pool_accessions_in_cell("X00001-3;Z99999;X00001", pool) == ["X00001"]
    assert mf.pool_accessions_in_cell("X00002;X00001-2", pool) == ["X00001", "X00002"]
    assert mf.pool_accessions_in_cell("Z99999", pool) == []
