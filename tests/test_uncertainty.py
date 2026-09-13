"""uncertainty.py on synthetic inputs: invented entries and counts; no network; no real biology."""
import logging
import math

import numpy as np

from pipeline import uncertainty as un, aggregate as ag

AA = ag.AA
LOG = logging.getLogger("t")


def test_lognormal_sigma_matches_the_published_sd_exactly():
    for m, s in ((100.0, 10.0), (100.0, 100.0), (30.0, 96.0), (5.0, 0.5)):
        sig = un.lognormal_sigma(m, s)
        var = (math.exp(sig * sig) - 1.0) * m * m * math.exp(sig * sig)      # variance of LogNormal(ln m, sig)
        assert math.isclose(math.sqrt(var), s, rel_tol=1e-12), (m, s)
    assert un.lognormal_sigma(100.0, 0.0) == 0.0
    assert un.lognormal_sigma(0.0, 5.0) == 0.0


def _entries(n, sd_factor):
    e = {}
    for i in range(n):
        acc = f"X{i:05d}"
        e[acc] = {"gene": f"G{i}", "tiers": ["1"] if i % 2 == 0 else ["2"], "mw": 1000.0 + i}
        for ft in ag.FIBER_TYPES:
            e[acc][f"v_{ft}"] = 10.0 + i
            e[acc][f"sd_{ft}"] = (10.0 + i) * sd_factor
            e[acc][f"valid_{ft}"] = 19
    return e


def _counts(n):
    c = {}
    for i in range(n):
        acc = f"X{i:05d}"
        count = {a: (i + j) % 5 for j, a in enumerate(AA)}
        c[acc] = {"count": count, "mw": 1000.0 + i, "free_frac": {a: 0.05 for a in AA}, "residue_frac": {a: 0.05 for a in AA}}
    return c


def _masses():
    return {a: {"free": 10.0 + i, "residue": 8.0 + i} for i, a in enumerate(AA)}


def test_simulation_is_deterministic_and_every_draw_is_a_composition():
    e, c, m = _entries(6, 0.5), _counts(6), _masses()
    a = un.simulate(e, c, m, draws=40, seed=7, log=LOG)
    b = un.simulate(e, c, m, draws=40, seed=7, log=LOG)
    for k in a:
        assert np.array_equal(a[k], b[k])
        assert a[k].shape == (40, 20)
        assert np.allclose(a[k].sum(axis=1), 1.0)
    assert not np.array_equal(a[("total", "I", "free", "between_fiber_spread")],
                              un.simulate(e, c, m, draws=40, seed=8, log=LOG)[("total", "I", "free", "between_fiber_spread")])


def test_zero_sd_reproduces_the_standard_in_every_draw():
    e, c, m = _entries(6, 0.0), _counts(6), _masses()
    sims = un.simulate(e, c, m, draws=10, seed=1, log=LOG)
    w = {acc: e[acc]["v_I"] * e[acc]["mw"] for acc in e}
    tot = sum(w.values())
    w = {k: v / tot for k, v in w.items()}
    p = ag.molar_profile(w, c, m, "free")
    P = sims[("total", "I", "free", "median_uncertainty")]
    assert np.allclose(P, np.array([[p[a] for a in AA]] * 10))


def test_median_uncertainty_is_narrower_than_between_fiber_spread():
    e, c, m = _entries(6, 0.8), _counts(6), _masses()
    sims = un.simulate(e, c, m, draws=400, seed=3, log=LOG)
    spread = sims[("total", "I", "free", "between_fiber_spread")].std(axis=0)
    med = sims[("total", "I", "free", "median_uncertainty")].std(axis=0)
    assert (med <= spread).all()
