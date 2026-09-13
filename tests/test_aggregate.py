"""aggregate.py rules on synthetic inputs: invented sequences, weights, and symbols; no real
biology; no files outside tmp_path; no network."""
import math

import pytest

from pipeline import common, aggregate as ag

AA = ag.AA


def _masses():
    # invented masses: free = 10 + index, residue = free - 2
    return {a: {"free": 10.0 + i, "residue": 8.0 + i, "three": a * 3, "name": a} for i, a in enumerate(AA)}


def _counts(spec):
    """spec: {acc: {letter: count}} -> the shape load_counts returns."""
    out = {}
    for acc, c in spec.items():
        count = {a: c.get(a, 0) for a in AA}
        n = sum(count.values())
        free = {a: count[a] * (10 + AA.index(a)) for a in AA}
        res = {a: count[a] * (8 + AA.index(a)) for a in AA}
        mw = sum(res.values()) + 2.0                              # residue mass sum + one "water"
        out[acc] = {"count": count, "mw": mw,
                    "free_frac": {a: free[a] / sum(free.values()) for a in AA},
                    "residue_frac": {a: res[a] / sum(res.values()) for a in AA}}
    return out


def test_molar_profile_of_one_protein_equals_its_own_composition():
    masses = _masses()
    counts = _counts({"X00001": {"A": 3, "G": 1}})
    p = ag.molar_profile({"X00001": 1.0}, counts, masses, "free")
    assert math.isclose(sum(p.values()), 1.0)
    assert math.isclose(p["A"], 3 * 10 / (3 * 10 + 1 * (10 + AA.index("G"))))
    assert math.isclose(p["G"], counts["X00001"]["free_frac"]["G"])
    r = ag.molar_profile({"X00001": 1.0}, counts, masses, "residue")
    assert math.isclose(r["A"], 3 * 8 / (3 * 8 + 1 * (8 + AA.index("G"))))


def test_molar_mixture_is_by_moles_not_by_fraction():
    """Two proteins of different mean residue mass: A1 (molar) and A1c (fraction mixture) differ,
    and A1 is the physical mixture — moles of each residue per gram of protein, summed."""
    masses = _masses()
    counts = _counts({"X00001": {"A": 10}, "X00002": {"Y": 10}})
    w = {"X00001": 0.5, "X00002": 0.5}
    p = ag.molar_profile(w, counts, masses, "residue")
    c = ag.fraction_mixture(w, counts, "residue")
    # by mass, half a gram of each: moles of A = 0.5/mw1 * 10, moles of Y = 0.5/mw2 * 10
    mw1, mw2 = counts["X00001"]["mw"], counts["X00002"]["mw"]
    nA, nY = 0.5 / mw1 * 10, 0.5 / mw2 * 10
    assert math.isclose(p["A"], nA * 8 / (nA * 8 + nY * (8 + AA.index("Y"))))
    assert math.isclose(c["A"], 0.5) and math.isclose(c["Y"], 0.5)     # the fraction mixture ignores the water
    assert p["A"] != c["A"]


def test_g_per_100g_residue_sums_to_100_less_terminal_water():
    masses = _masses()
    counts = _counts({"X00001": {"A": 4, "C": 6}})
    g = ag.g_per_100g_protein({"X00001": 1.0}, counts, masses, "residue")
    mw = counts["X00001"]["mw"]
    assert math.isclose(sum(g.values()), 100.0 * (mw - 2.0) / mw)


def test_rescaled_weights_hits_target_and_renormalises():
    w = {"M1": 0.4, "M2": 0.2, "A1": 0.2, "Z": 0.2}
    new, r0 = ag.rescaled_weights(w, {"M1", "M2"}, {"A1"}, 1.5, "scale_mhc")
    assert math.isclose(r0, 3.0)
    assert math.isclose(sum(new.values()), 1.0)
    assert math.isclose((new["M1"] + new["M2"]) / new["A1"], 1.5)
    assert math.isclose(new["M1"] / new["M2"], 2.0)                   # members keep their proportions
    new2, _ = ag.rescaled_weights(w, {"M1", "M2"}, {"A1"}, 1.5, "scale_actin")
    assert math.isclose((new2["M1"] + new2["M2"]) / new2["A1"], 1.5)
    assert math.isclose(new2["Z"] / new2["M1"], 0.5)                   # untouched entries keep their ratio to the fixed band
    same, _ = ag.rescaled_weights(w, {"M1", "M2"}, {"A1"}, 3.0, "scale_mhc")
    assert all(math.isclose(same[k], w[k]) for k in w)                 # at its own ratio nothing changes


def _write_symbols(path):
    cp = common.new_ini()
    for i, a in enumerate(AA):
        cp.add_section(a)
        cp[a]["trivial_name"] = f"name{i}"
        cp[a]["three_letter"] = f"{a}{a.lower()}{a.lower()}"          # invented three-letter symbols
    common.write_ini(cp, path, ["synthetic"])


def _write_eaa(path, headings, groups, filled=True):
    lines = ["[meta]", "source_id = x", "decision = D62", f"transcribed_by = {'t' if filled else '___'}", "transcribed = t", "",
             "[indispensable_amino_acids]", f"headings = {headings}", "location = tbl", "as_reported = hdr", ""]
    for g, members in groups.items():
        lines += [f"[group.{g}]", f"members = {members}", "location = tbl", "as_reported = key", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def test_load_eaa_resolves_headings_and_groups(tmp_path, monkeypatch):
    sym = tmp_path / "symbols.ini"
    _write_symbols(sym)
    eaa = tmp_path / "eaa.ini"
    _write_eaa(eaa, "Aaa, Ggg, GRP", {"GRP": "Sss, Ttt"})
    monkeypatch.setitem(ag.INPUTS, "symbols", sym)
    monkeypatch.setitem(ag.INPUTS, "eaa", eaa)
    out = ag.load_eaa()
    assert out["headings"] == [("Aaa", ["A"], "single"), ("Ggg", ["G"], "single"), ("GRP", ["S", "T"], "group")]


def test_load_eaa_stops_on_placeholder_or_unknown_symbol(tmp_path, monkeypatch):
    sym = tmp_path / "symbols.ini"
    _write_symbols(sym)
    monkeypatch.setitem(ag.INPUTS, "symbols", sym)
    eaa = tmp_path / "eaa.ini"
    _write_eaa(eaa, "___", {})
    monkeypatch.setitem(ag.INPUTS, "eaa", eaa)
    with pytest.raises(SystemExit):
        ag.load_eaa()
    _write_eaa(eaa, "Aaa;Zzz", {})                                      # Zzz is not in the synthetic table
    with pytest.raises(SystemExit):
        ag.load_eaa()


def test_tsv_text_marks_generated():
    t = ag.tsv_text(["h"], ["a", "b"], [[1, 2]])
    assert t.startswith("# GENERATED by aggregate.py") and t.endswith("a\tb\n1\t2\n")


def test_split_list_accepts_commas_and_semicolons():
    assert ag.split_list("a, b;c ,d") == ["a", "b", "c", "d"]
