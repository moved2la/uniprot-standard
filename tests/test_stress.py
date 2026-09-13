"""stress.py on synthetic inputs: invented entries, counts, and settings; no real biology."""
import math

from pipeline import stress as stt, aggregate as ag

AA = ag.AA


def _counts(n):
    c = {}
    for i in range(n):
        acc = f"X{i:05d}"
        count = {a: (i * 3 + j) % 4 + 1 for j, a in enumerate(AA)}
        c[acc] = {"count": count, "mw": 1000.0 * (i + 1), "free_frac": {a: 0.05 for a in AA}, "residue_frac": {a: 0.05 for a in AA}}
    return c


def _masses():
    return {a: {"free": 10.0 + i, "residue": 8.0 + i} for i, a in enumerate(AA)}


def test_tier_ratio_at_the_measured_share_reproduces_the_standard_and_knockout_renormalises():
    counts, masses = _counts(4), _masses()
    w0 = {"X00000": 0.4, "X00001": 0.3, "X00002": 0.2, "X00003": 0.1}
    tiers = {"X00000": "1", "X00001": "1", "X00002": "2", "X00003": "2"}
    base = ag.molar_profile(w0, counts, masses, "free")
    share1 = 0.7
    wk = {a: w * (share1 / 0.7 if tiers[a] == "1" else (1 - share1) / 0.3) for a, w in w0.items()}
    p = ag.molar_profile(wk, counts, masses, "free")
    assert all(math.isclose(p[a], base[a]) for a in AA)
    drop = {"X00000"}
    wk = {a: w for a, w in w0.items() if a not in drop}
    tot = sum(wk.values())
    wk = {a: w / tot for a, w in wk.items()}
    assert math.isclose(sum(wk.values()), 1.0)
    p = ag.molar_profile(wk, counts, masses, "free")
    assert math.isclose(sum(p.values()), 1.0) and any(abs(p[a] - base[a]) > 1e-6 for a in AA)


def test_max_shift_reports_the_amino_acid_and_percent():
    base = {a: 0.05 for a in AA}
    stressed = dict(base)
    stressed["K"] = 0.06
    m, a, pct = stt.max_shift(base, stressed)
    assert math.isclose(m, 0.01) and a == "K" and math.isclose(pct, 20.0)


def test_settings_lists_parse_with_semicolons_and_stop_on_placeholder(tmp_path, monkeypatch):
    ini = tmp_path / "agg.ini"
    ini.write_text("[stress]\ncomposition_distance_top_entries = 3\nconvergence_top_k = 2;3\nknockout_top_k = 1\n"
                   "tier1_share = 0.5;0.9\nsize_tilt_alpha = -0.5;0.5\nrandom_abuse_factor = 2\nrandom_abuse_draws = 5\nseed = 1\n", encoding="utf-8")
    monkeypatch.setattr(stt.common, "AGGREGATION_DECISIONS_INI", ini)
    st = stt.load_settings()
    assert st["convergence_k"] == [2, 3] and st["alphas"] == [-0.5, 0.5] and st["random_draws"] == 5
    ini.write_text("[stress]\ncomposition_distance_top_entries = ___\n", encoding="utf-8")
    try:
        stt.load_settings()
        raise AssertionError("did not stop")
    except SystemExit:
        pass
