"""mass_fractions.py rules on synthetic inputs. No real biology; no network; no files
outside tmp_path. Gene symbols and accessions here are invented."""
import configparser
import logging
import math

import pytest

from pipeline import common
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


def test_per_entry_table_has_every_weight_and_sums_to_one():
    """D61: the one-table form carries within-tier and combined weights with low/high, one
    row per pool accession, provenance in the header; every weight column sums to one."""
    pool = _pool(("X00001", "GA", "1"), ("X00002", "GB", "1"), ("X00003", "GC", "2"))
    j = mf.join([_row(4, "GA", v=(30.0, 30.0, 30.0)), _row(5, "GB", v=(10.0, 10.0, 10.0)), _row(6, "GC", v=(20.0, 20.0, 20.0))], pool, LOG)
    measured = j["measured"]
    mw = {"X00001": 1000.0, "X00002": 2000.0, "X00003": 500.0}
    w = mf.weights(measured, mw, LOG)
    meta = {"source_line": "S", "retrieved": "T", "sheet": "Sheet1", "header_row": 3,
            "columns": {f"{k}_{ft}": f"{k} {ft}" for ft in mf.FIBER_TYPES for k in ("median", "sd", "valid")}}
    text = mf.per_entry_table(measured, w["tiers"], meta, ["generated   = now"])
    assert text.startswith("# GENERATED by mass_fractions.py")
    assert "BETWEEN FIBERS" in text                    # the D45 reading of the interval, stated in the header
    lines = [l for l in text.splitlines() if not l.startswith("#")]
    hdr = lines[0].split("\t")
    rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]
    assert [r["accession"] for r in rows] == ["X00001", "X00002", "X00003"]
    for col in ("w_I_tier1", "w_IIa_tier2", "w_I_combined", "w_IIx_combined"):
        assert col in hdr and f"{col}_low" in hdr and f"{col}_high" in hdr
        assert math.isclose(sum(float(r[col]) for r in rows if r[col] != ""), 1.0)
    assert rows[2]["w_I_tier1"] == "" and rows[2]["w_I_tier2"] != ""      # a tier-2 entry has no tier-1 weight
    assert math.isclose(float(rows[0]["w_I_combined"]), 30.0 * 1000 / (30.0 * 1000 + 10.0 * 2000 + 20.0 * 500))


# --------------------------------------------------------------------------- category scoping
#
# config/literature_sources.ini is shared by every category, so `role = primary` is a statement
# about a source FOR a category: blood's two primaries are as legitimate as muscle's one. Before
# these tests the stage scanned the whole file and stopped with three primaries. (Step 7b)

def _sources(*rows) -> configparser.ConfigParser:
    """A synthetic literature_sources.ini: (source_id, role, used_for[, extra keys])."""
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    for sid, role, used_for, *extra in rows:
        cp[f"source.{sid}"] = {"role": role, "used_for": used_for, **(extra[0] if extra else {})}
    return cp


def _decisions(category: str | None) -> configparser.ConfigParser:
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    if category is not None:
        cp["scope"] = {"category": category}
    return cp


def test_primary_is_the_one_used_for_this_stages_category():
    src = _sources(("muscle_source", "primary", "skeletal_muscle"),
                   ("blood_red", "primary", "blood"),
                   ("blood_plasma", "primary", "blood"))
    assert mf.sources_with_role(src, "primary", "skeletal_muscle") == ["muscle_source"]
    assert mf.sources_with_role(src, "primary", "blood") == ["blood_red", "blood_plasma"]
    assert mf.sources_with_role(src, "primary", "liver") == []


def test_a_source_used_for_several_categories_is_found_by_each():
    src = _sources(("shared", "method_citation", "original, skeletal_muscle, non_protein_metabolites"))
    for cat in ("original", "skeletal_muscle", "non_protein_metabolites"):
        assert mf.sources_with_role(src, "method_citation", cat) == ["shared"]
    assert mf.sources_with_role(src, "method_citation", "blood") == []


def test_role_holders_hint_names_every_primary_and_what_it_is_used_for():
    """The STOP says which categories the primaries belong to, so a typo in [meta] category
    reads as a typo rather than as a missing source."""
    src = _sources(("muscle_source", "primary", "skeletal_muscle"), ("blood_red", "primary", "blood"))
    hint = mf.role_holders_by_category(src, "primary")
    assert "muscle_source (used_for = skeletal_muscle)" in hint and "blood_red (used_for = blood)" in hint


def test_stage_category_stops_when_it_is_missing_or_a_placeholder():
    assert mf.stage_category(_decisions("skeletal_muscle")) == "skeletal_muscle"
    for dec in (_decisions(None), _decisions("___"), _decisions("")):
        with pytest.raises(SystemExit) as e:
            mf.stage_category(dec)
        assert "[scope]" in str(e.value)


# --------------------------------------------------------------------------- the real config

def test_this_repository_declares_a_category_with_exactly_one_primary():
    dec = common.read_ini(common.CONFIG_DIR / "mass_fraction_decisions.ini")
    src = common.read_ini(common.LITERATURE_SOURCES_INI)
    category = mf.stage_category(dec)
    declared = [c.strip() for c in src["categories"].get("columns", "").split(",") if c.strip()]
    assert category in declared, \
        f"[meta] category = {category} is not a column of [categories] in literature_sources.ini: {declared}"
    assert len(mf.sources_with_role(src, "primary", category)) == 1, \
        f"category {category} must have exactly one source with role = primary"
    