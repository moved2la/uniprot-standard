"""bound_pools.py rules on synthetic inputs: invented compositions, invented pool values, invented
masses; no real biology; no files outside tmp_path; no network."""
import math

import pytest

from pipeline import aggregate as ag, bound_pools as bp

AA = bp.AA


def _g100(letter_share: float, letter: str, total: float = 116.0) -> dict[str, float]:
    """An invented g-per-100-g-protein vector: `letter` gets letter_share of `total`, the rest is spread evenly."""
    rest = (total - letter_share * total) / 19
    return {a: (letter_share * total if a == letter else rest) for a in AA}


def test_adjusted_profile_adds_the_pool_to_one_amino_acid_only():
    g100 = _g100(0.025, "H")
    r = bp.adjusted_profile(g100, protein_g_per_kg=200.0, pool_mmol_per_kg=10.0, letter="H", n_per_mole=1.0, free_mass=155.0)
    # F_a = P x g100 / 100 ; C = c/1000 x n x m
    assert math.isclose(r["F"]["H"], 200.0 * g100["H"] / 100.0)
    assert math.isclose(r["C"], 10.0 / 1000.0 * 155.0)
    assert math.isclose(r["ratio"], r["C"] / r["F"]["H"])
    assert math.isclose(sum(r["fractions"].values()), 1.0)
    tot = sum(r["F"].values()) + r["C"]
    assert math.isclose(r["fractions"]["H"], (r["F"]["H"] + r["C"]) / tot)
    assert math.isclose(r["fractions"]["A"], r["F"]["A"] / tot)
    # every other amino acid moves down by the same factor
    before = g100["A"] / sum(g100.values())
    assert r["fractions"]["A"] < before
    assert math.isclose(r["fractions"]["A"] / before, r["fractions"]["G"] / (g100["G"] / sum(g100.values())))


def test_zero_pool_reproduces_the_protein_only_profile():
    g100 = _g100(0.03, "H")
    r = bp.adjusted_profile(g100, 180.0, 0.0, "H", 1.0, 155.0)
    tot = sum(g100.values())
    for a in AA:
        assert math.isclose(r["fractions"][a], g100[a] / tot)


def test_stoichiometry_scales_the_pool():
    g100 = _g100(0.03, "H")
    one = bp.adjusted_profile(g100, 180.0, 5.0, "H", 1.0, 155.0)
    two = bp.adjusted_profile(g100, 180.0, 5.0, "H", 2.0, 155.0)
    assert math.isclose(two["C"], 2 * one["C"])


def test_lognormal_bounds_bracket_the_mean_and_are_positive():
    lo, hi = bp.lognormal_bounds(20.0, 15.0)
    assert 0 < lo < 20.0 < hi
    lo2, hi2 = bp.lognormal_bounds(20.0, 1.0)
    assert hi2 - lo2 < hi - lo                        # a smaller SD gives a narrower band


def test_rate_conversions():
    assert math.isclose(bp.k_per_day_from_half_life_weeks(1.0), math.log(2) / 7)
    assert math.isclose(bp.k_per_day_from_fsr_percent_per_hour(0.05), 0.012)


# ----------------------------------------------------------------------------
# end to end on synthetic files
# ----------------------------------------------------------------------------

def _ini_masses(path):
    lines = ["# synthetic"]
    for i, a in enumerate(AA):
        lines += [f"[{a}]", f"trivial_name = name{a}", f"three_letter = {a}{a.lower()}x", f"free_mass = {100 + i}.00", f"residue_mass = {82 + i}.00", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def _ini_symbols(path):
    lines = ["# synthetic"]
    for a in AA:
        lines += [f"[{a}]", f"three_letter = {a}{a.lower()}x", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def _tsv(path, columns, rows):
    path.write_text("# GENERATED (synthetic)\n" + "\t".join(columns) + "\n" + "\n".join("\t".join(str(v) for v in r) for r in rows) + "\n", encoding="utf-8")


def _standard_files(tmp_path, letter="H"):
    cols = ["amino_acid", "three_letter", "name"] + [f"{pf}_{ft}" for pf, ft in bp.COL_ORDER] + ["standard"]
    g_rows, s_rows = [], []
    g100 = {ft: _g100(0.02 + 0.005 * i, letter) for i, ft in enumerate(bp.FIBER_TYPES)}
    for a in AA:
        g = [a, a * 3, a]
        s = [a, a * 3, a]
        for pf, ft in bp.COL_ORDER:
            g.append(f"{g100[ft][a]:.4f}")
            s.append(f"{100 * g100[ft][a] / sum(g100[ft].values()):.4f}")
        g.append("")
        s.append(f"{sum(100 * g100[ft][a] / sum(g100[ft].values()) for ft in bp.FIBER_TYPES) / 3:.4f}")
        g_rows.append(g)
        s_rows.append(s)
    g_rows.append(["sum", "", ""] + [f"{sum(g100[ft].values()):.4f}" for pf, ft in bp.COL_ORDER] + [""])
    s_rows.append(["sum", "", ""] + ["100.0000"] * 9 + ["100.0000"])
    _tsv(tmp_path / "amino_acid_g_per_100g_protein_free.tsv", cols, g_rows)
    _tsv(tmp_path / "_calculated_amino_acid_standard.tsv", cols, s_rows)
    return g100


def _mix(path):
    path.write_text("[meta]\nsource_id = synthetic\n\n[shares]\nI = 1\nIIa = 1\nIIx = 1\nbasis = synthetic\nlocation = synthetic\n", encoding="utf-8")


POOL_INI = """
[meta]
decisions = D0
version_label = synthetic

[pool.testpool]
amino_acid = Hhx
moles_amino_acid_per_mole = 1

[testpool.single_fiber.slow]
value = 10
spread = 5
spread_kind = synthetic
basis = dry
applies_to = I
source_id = synthetic

[testpool.single_fiber.fast]
value = 20
spread = 10
spread_kind = synthetic
basis = dry
applies_to = IIa, IIx
source_id = synthetic

[testpool.whole_muscle]
value = 4
basis = wet
source_id = synthetic

[testpool.whole_muscle.groupA]
value = 5
basis = wet
source_id = synthetic

[protein_content]
value = 200
basis = wet
source_id = synthetic

[water_content]
value = 750
source_id = synthetic

[turnover.testpool.a]
t_half_weeks = 5
source_id = synthetic

[turnover.muscle_protein.low]
value = 0.03
source_id = synthetic

[turnover.muscle_protein.high]
value = 0.09
source_id = synthetic
"""


class _Log:
    def info(self, *a, **k): pass
    def warning(self, *a, **k): pass
    def error(self, *a, **k): pass


@pytest.fixture
def synthetic_repo(tmp_path, monkeypatch):
    out = tmp_path / "outputs"
    out.mkdir()
    g100 = _standard_files(out)
    masses = tmp_path / "masses.ini"
    symbols = tmp_path / "symbols.ini"
    mix = tmp_path / "mix.ini"
    pools = tmp_path / "pools.ini"
    _ini_masses(masses)
    _ini_symbols(symbols)
    _mix(mix)
    pools.write_text(POOL_INI, encoding="utf-8")
    manifest = tmp_path / "manifest.ini"
    manifest.write_text("[file.synthetic.file.1]\nsource_id = synthetic\nsha256 = 00ff\n", encoding="utf-8")
    monkeypatch.setattr(bp, "INPUTS", {"pools": pools, "g100_free": out / "amino_acid_g_per_100g_protein_free.tsv",
                                       "standard": out / "_calculated_amino_acid_standard.tsv", "fiber_type_mix": mix,
                                       "masses": masses, "symbols": symbols, "manifest": manifest})
    monkeypatch.setattr(bp, "OUT_DIR", out)
    monkeypatch.setitem(ag.INPUTS, "masses", masses)
    monkeypatch.setitem(ag.INPUTS, "fiber_type_mix", mix)
    monkeypatch.setattr(bp.common, "REPO_ROOT", tmp_path)
    return {"tmp": tmp_path, "out": out, "g100": g100, "pools": pools, "manifest": manifest}


def test_build_writes_the_adjusted_table_and_the_sensitivities(synthetic_repo):
    files = bp.build(_Log())
    assert bp.TABLE in files
    for name in ("bound_pool_adjustment_per_amino_acid.tsv", "bound_pool_amounts_per_kg_muscle.tsv",
                 "sensitivity_bound_pool_turnover_frame.tsv", "sensitivity_bound_pool_basis.tsv",
                 "sensitivity_bound_pool_sex.tsv", "sensitivity_bound_pool_spread.tsv", "bound_pools_summary.ini"):
        assert name in files, name
    body = [ln.split("\t") for ln in files[bp.TABLE].splitlines() if not ln.startswith("#")]
    header, rows = body[0], body[1:]
    by_aa = {r[0]: dict(zip(header, r)) for r in rows}
    # every column sums to 100 (the sum row is computed, not typed); the synthetic fixture's nineteen equal
    # values round in the same direction, so the tolerance here is looser than the real file's 5e-4
    for col in header[3:]:
        assert abs(sum(float(by_aa[a][col]) for a in AA) - 100.0) < 2e-3, col
    # contractile and builders copied unchanged; total moved on the pool's amino acid
    std = {r["amino_acid"]: r for r in ag.read_tsv_skip_comments(synthetic_repo["out"] / "_calculated_amino_acid_standard.tsv")}
    for ft in bp.FIBER_TYPES:
        assert by_aa["A"][f"contractile_{ft}"] == std["A"][f"contractile_{ft}"]
        assert by_aa["A"][f"builders_{ft}"] == std["A"][f"builders_{ft}"]
        assert float(by_aa["H"][f"total_{ft}"]) > float(std["H"][f"total_{ft}"])
        assert float(by_aa["A"][f"total_{ft}"]) < float(std["A"][f"total_{ft}"])
    # the fast value (20) is applied to IIa and IIx alike, the slow value (10) to I: the ratio doubles
    per_kg = [ln.split("\t") for ln in files["bound_pool_amounts_per_kg_muscle.tsv"].splitlines() if not ln.startswith("#")]
    h, rs = per_kg[0], per_kg[1:]
    ratio = {r[0]: float(r[h.index("ratio_pool_to_protein")]) for r in rs}
    pool_g = {r[0]: float(r[h.index("pool_Hhx_g_per_kg")]) for r in rs}
    assert math.isclose(pool_g["IIa"], pool_g["IIx"]) and math.isclose(pool_g["IIa"], 2 * pool_g["I"])
    # protein 200 g/kg wet with 750 g water -> 800 g/kg dry; the slow ratio = (10/1000 x 107) / (800 x g100_H/100)
    g_h = synthetic_repo["g100"]["I"]["H"]
    assert math.isclose(ratio["I"], (10 / 1000 * (100 + AA.index("H"))) / (800.0 * g_h / 100.0), rel_tol=1e-9)


def test_placeholders_leave_the_table_unwritten_and_say_so(synthetic_repo):
    synthetic_repo["pools"].write_text(POOL_INI.replace("[protein_content]\nvalue = 200", "[protein_content]\nvalue = ___"), encoding="utf-8")
    files = bp.build(_Log())
    assert bp.TABLE not in files
    assert "NOT computed" in files["bound_pools_summary.ini"]
    assert "protein_content" in files["bound_pools_summary.ini"]


def test_a_missing_sensitivity_input_skips_only_that_sensitivity(synthetic_repo):
    synthetic_repo["pools"].write_text(POOL_INI.replace("t_half_weeks = 5", "t_half_weeks = ___"), encoding="utf-8")
    files = bp.build(_Log())
    assert bp.TABLE in files
    assert "sensitivity_bound_pool_turnover_frame.tsv" not in files
    assert "sensitivity_bound_pool_basis.tsv" in files
    assert "sensitivity_turnover_frame = NOT computed" in files["bound_pools_summary.ini"]


def test_check_mode_is_current_after_a_build(synthetic_repo):
    assert bp.main([]) == 0
    assert bp.main(["--check"]) == 0
    (synthetic_repo["out"] / bp.TABLE).write_text("tampered\n", encoding="utf-8")
    assert bp.main(["--check"]) == 1


def test_green_light_refuses_a_source_without_a_hashed_file(synthetic_repo):
    synthetic_repo["manifest"].write_text("[file.synthetic.file.1]\nsource_id = synthetic\nsha256 =\n", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        bp.build(_Log())
    assert "green light refused" in str(e.value.code) and "synthetic" in str(e.value.code)


def test_green_light_refuses_when_the_manifest_is_absent(synthetic_repo):
    synthetic_repo["manifest"].unlink()
    with pytest.raises(SystemExit) as e:
        bp.build(_Log())
    assert "green light" in str(e.value.code)


def test_missing_whole_muscle_sections_skip_only_those_sensitivities(synthetic_repo):
    ini = POOL_INI
    i = ini.index("[testpool.whole_muscle]"); j = ini.index("[protein_content]")
    synthetic_repo["pools"].write_text(ini[:i] + ini[j:], encoding="utf-8")
    files = bp.build(_Log())
    assert bp.TABLE in files
    assert "sensitivity_bound_pool_basis.tsv" not in files and "sensitivity_bound_pool_sex.tsv" not in files
    assert "sensitivity_bound_pool_spread.tsv" in files and "sensitivity_bound_pool_turnover_frame.tsv" in files
    line = [ln for ln in files["bound_pools_summary.ini"].splitlines() if ln.startswith("sensitivity_basis")][0]
    assert "NOT computed" in line
