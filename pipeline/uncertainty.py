#!/usr/bin/env python3
"""
uncertainty.py — offline stage: the standard's uncertainty, propagated from Layer B, and every
uncertainty term on one scale.

Two different quantities come out of the between-fiber data, and both are reported (D63):
  between_fiber_spread    what a random pure fiber of the type looks like — wide; biology plus
                          single-fiber measurement noise. The D45 reading.
  median_uncertainty      how well the fiber type's median profile — the standard — is known
                          from the n fibers measured. Narrow. The one to state beside the number.

Both are Monte Carlo in log space (D63, A3):
  per entry and fiber type, a log-normal is fitted to the published median m and standard
  deviation s exactly (two numbers, two parameters): sigma^2 = ln x with x = (1 + sqrt(1 + 4 s^2/m^2)) / 2,
  mu = ln m. That is the one assumption — the between-fiber distribution of an intensity is
  log-normal, the field's standard model — and it is stated in every header.
  between_fiber_spread: v_i ~ LogNormal(mu_i, sigma_i), independent per entry.
  median_uncertainty:   log v_i ~ Normal(mu_i, sigma_i / sqrt(n_i)), n_i = the entry's valid-value
                        count for the fiber type (the sampling error of a log-median).
  Every draw recomputes the whole profile (A1), so each draw is a composition. Entries with m = 0
  or s = 0 are fixed. Draws and seed come from config/aggregation_decisions.ini [monte_carlo].
  The profile of a draw depends on v and the residue counts only: n_a proportional to
  sum_i v_i * count_{i,a} (the MW in the weight cancels against the MW in the molar conversion).

Reads
  config/mass_fractions_per_entry.tsv, outputs/composition/amino_acid_composition_per_protein.tsv,
  data/pubchem/amino_acid_masses.ini, config/aggregation_decisions.ini [monte_carlo],
  outputs/standard/amino_acid_profiles.tsv, sensitivity_mhc_actin_spread.tsv, bounds_after_weighting.tsv,
  completeness_sensitivity.tsv (from aggregate.py)

Writes (outputs/standard/)
  uncertainty_intervals.tsv         per profile, fiber type, convention, term, amino acid: the standard,
                                    the 2.5th / 50th / 97.5th percentiles over draws, the half-width
  uncertainty_per_amino_acid.tsv    every term on one scale — the maximum absolute shift from the standard
                                    (Monte Carlo: the interval's half-width; bounds: their worst-case sum;
                                    sensitivity: max_abs_shift) — and the widest term named (provisional scale)
  sd_to_median_ratio.tsv            per entry and fiber type: m, s, s/m, the fitted sigma, n — the regime table
  uncertainty_summary.ini
  logs/uncertainty_<UTC>.log
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

from pipeline import common
from pipeline.mass_fractions import read_tsv_skip_comments, fnum, sha256_path
from pipeline import aggregate as ag

OUT_DIR = common.STANDARD_DIR
AA = ag.AA
FIBER_TYPES = ag.FIBER_TYPES
PROFILES = ag.PROFILES
CONVENTIONS = ag.CONVENTIONS
TERMS_MC = ("between_fiber_spread", "median_uncertainty")
CHUNK = 500


def lognormal_sigma(m: float, s: float) -> float:
    """sigma of the log-normal with median m and standard deviation s (exact moment match)."""
    if m <= 0 or s <= 0:
        return 0.0
    x = (1.0 + math.sqrt(1.0 + 4.0 * (s * s) / (m * m))) / 2.0
    return math.sqrt(math.log(x))


def load_monte_carlo_settings() -> tuple[int, int]:
    cp = common.read_ini(common.AGGREGATION_DECISIONS_INI)
    if "monte_carlo" not in cp:
        raise SystemExit(f"[STOP] {common.AGGREGATION_DECISIONS_INI.name} has no [monte_carlo] section (draws, seed; D63)")
    sec = cp["monte_carlo"]
    for k in ("draws", "seed"):
        if k not in sec or sec[k].strip() in ("", "___"):
            raise SystemExit(f"[STOP] [monte_carlo] {k} is not set (D63)")
    return int(sec["draws"]), int(sec["seed"])


def read_standard_profiles() -> dict[tuple, dict[str, float]]:
    out = {}
    for r in read_tsv_skip_comments(OUT_DIR / "amino_acid_profiles.tsv"):
        out[(r["profile"], r["fiber_type"], r["convention"])] = {a: float(r[f"frac_{a}"]) for a in AA}
    if not out:
        raise SystemExit("[STOP] outputs/standard/amino_acid_profiles.tsv is empty; run aggregate first")
    return out


def simulate(entries, counts, masses, draws: int, seed: int, log) -> dict[tuple, np.ndarray]:
    """Returns, per (profile, fiber type, convention, term), an array (draws, 20) of profile fractions."""
    accs = sorted(entries)
    C = np.array([[counts[a]["count"][x] for x in AA] for a in accs], dtype=float)          # entries x 20
    mass = {conv: np.array([masses[x][conv] for x in AA]) for conv in CONVENTIONS}
    tier1 = np.array(["1" in entries[a]["tiers"] for a in accs])
    rng = np.random.default_rng(seed)
    results: dict[tuple, list[np.ndarray]] = {}
    for ft in FIBER_TYPES:
        m = np.array([entries[a][f"v_{ft}"] for a in accs])
        s = np.array([entries[a][f"sd_{ft}"] for a in accs])
        n = np.array([max(entries[a][f"valid_{ft}"], 1) for a in accs], dtype=float)
        sigma = np.array([lognormal_sigma(mi, si) for mi, si in zip(m, s)])
        mu = np.where(m > 0, np.log(np.where(m > 0, m, 1.0)), 0.0)
        for term in TERMS_MC:
            sig = sigma if term == "between_fiber_spread" else sigma / np.sqrt(n)
            done = 0
            while done < draws:
                k = min(CHUNK, draws - done)
                z = rng.standard_normal((k, len(accs)))
                V = np.where(m > 0, np.exp(mu + sig * z), 0.0)                           # k x entries
                for profile in PROFILES:
                    Vp = V if profile == "combined" else V * tier1
                    N = Vp @ C                                                        # k x 20: molar amounts (up to a constant)
                    for conv in CONVENTIONS:
                        M = N * mass[conv]
                        P = M / M.sum(axis=1, keepdims=True)
                        results.setdefault((profile, ft, conv, term), []).append(P)
                done += k
            log.info("type %s %s: %d draws", ft, term, draws)
    return {k: np.concatenate(v, axis=0) for k, v in results.items()}


def build(log) -> tuple[dict[str, str], dict]:
    draws, seed = load_monte_carlo_settings()
    masses = ag.load_masses()
    entries, _ = ag.load_entries()
    counts = ag.load_counts(entries)
    standard = read_standard_profiles()
    inputs = {"weights": common.MASS_FRACTIONS_TSV, "composition": common.COMPOSITION_TSV, "masses": common.AMINO_ACID_MASSES_INI,
              "settings": common.AGGREGATION_DECISIONS_INI, "profiles": OUT_DIR / "amino_acid_profiles.tsv",
              "mhc_actin": OUT_DIR / "sensitivity_mhc_actin_spread.tsv", "bounds": OUT_DIR / "bounds_after_weighting.tsv",
              "completeness": OUT_DIR / "completeness_sensitivity.tsv"}
    hashes = {k: f"{p.relative_to(common.REPO_ROOT).as_posix()} sha256 {sha256_path(p)}" for k, p in inputs.items()}
    header_common = [f"generated = {common.iso_now()}"] + [f"input.{k} = {v}" for k, v in hashes.items()] + [
        f"A3: Monte Carlo in log space, {draws} draws, seed {seed}, independent per entry, the whole profile recomputed per draw (D63).",
        "Assumption: the between-fiber distribution of each entry's value is log-normal, fitted exactly to the published median and SD.",
        "between_fiber_spread = a random pure fiber of the type (the D45 reading), NOT the uncertainty of the standard;",
        "median_uncertainty = the sampling uncertainty of the fiber type's median profile from n fibers (sigma/sqrt(n) in log space): the interval to state beside the standard.",
        "The dataset's anchor protein has SD 0 by construction; its fiber-to-fiber variation is inside every other entry's SD.",
    ]
    sims = simulate(entries, counts, masses, draws, seed, log)

    files: dict[str, str] = {}
    # intervals
    rows = []
    halfwidth: dict[tuple, dict[str, float]] = {}
    for (profile, ft, conv, term), P in sorted(sims.items()):
        lo, med, hi = np.percentile(P, [2.5, 50, 97.5], axis=0)
        std = standard[(profile, ft, conv)]
        hw = {}
        for i, a in enumerate(AA):
            hw[a] = max(abs(hi[i] - std[a]), abs(lo[i] - std[a]))
            rows.append([profile, ft, conv, term, a, fnum(std[a]), fnum(lo[i]), fnum(med[i]), fnum(hi[i]), fnum(hi[i] - lo[i]), fnum(hw[a])])
        halfwidth[(profile, ft, conv, term)] = hw
    files["uncertainty_intervals.tsv"] = ag.tsv_text(tool="uncertainty.py", header_lines=header_common + ["percentiles over draws of each amino acid's fraction; half_width = max(|p97.5 - standard|, |p2.5 - standard|)"],
        columns=["profile", "fiber_type", "convention", "term", "amino_acid", "standard", "p2_5", "p50", "p97_5", "width", "half_width"], rows=rows)

    # one scale: every term as a maximum absolute shift from the standard
    mhc = {}
    for r in read_tsv_skip_comments(inputs["mhc_actin"]):
        mhc[(r["profile"], r["fiber_type"], r["convention"], r["amino_acid"])] = float(r["max_abs_shift"])
    bounds: dict[tuple, float] = {}
    for r in read_tsv_skip_comments(inputs["bounds"]):
        bounds[(r["profile"], r["fiber_type"], r["bound"], r["amino_acid"])] = float(r["worst_case_sum"])
    comp = {}
    for r in read_tsv_skip_comments(inputs["completeness"]):
        if r["scope"] == "all":
            comp[r["fiber_type"]] = float(r["bound_on_any_fraction"])
    bound_kinds = sorted({k[2] for k in bounds if not k[2].startswith("shared_row_excess_")})
    term_cols = list(TERMS_MC) + ["mhc_actin_range"] + bound_kinds + ["shared_row_excess", "completeness_all_outside_genes"]
    rows = []
    widest_count: dict[tuple, dict[str, int]] = {}
    for profile in PROFILES:
        for ft in FIBER_TYPES:
            for conv in CONVENTIONS:
                std = standard[(profile, ft, conv)]
                for a in AA:
                    vals = {}
                    for term in TERMS_MC:
                        vals[term] = halfwidth[(profile, ft, conv, term)][a]
                    vals["mhc_actin_range"] = mhc.get((profile, ft, conv, a), float("nan"))
                    for kind in bound_kinds:
                        vals[kind] = bounds.get((profile, ft, kind, a), float("nan"))
                    vals["shared_row_excess"] = bounds.get((profile, ft, f"shared_row_excess_{conv}", a), float("nan"))
                    vals["completeness_all_outside_genes"] = comp.get(ft, float("nan"))
                    finite = {k: v for k, v in vals.items() if not math.isnan(v) and k != "glycosylated_entries_sum_w"}
                    # widest_term: among the terms that describe THIS amino acid's uncertainty of the standard —
                    # excluding between_fiber_spread (not an uncertainty of the standard) and completeness (one
                    # global mass share, the same for every amino acid; its own column).
                    per_aa = {k: v for k, v in finite.items() if k not in ("between_fiber_spread", "completeness_all_outside_genes")}
                    widest = max(per_aa, key=lambda k: per_aa[k]) if per_aa else ""
                    widest_all = max(finite, key=lambda k: finite[k]) if finite else ""
                    rows.append([profile, ft, conv, a, fnum(std[a])] + [fnum(vals[t]) for t in term_cols] + [widest, widest_all])
                    widest_count.setdefault((profile, ft, conv), {}).setdefault(widest, 0)
                    widest_count[(profile, ft, conv)][widest] += 1
    files["uncertainty_per_amino_acid.tsv"] = ag.tsv_text(tool="uncertainty.py", header_lines=header_common + ["one scale (provisional): every term as the maximum absolute shift of the fraction from the standard — Monte Carlo terms: the half-width of the 95 % interval; mhc_actin_range: max |profile - standard| over the measured ratios (A6);",
                         "bound terms: the worst-case sum over entries or families (bounds_after_weighting.tsv); shared_row_excess: D48; completeness: the mass share of every outside gene (A7).",
                         "between_fiber_spread is not an uncertainty of the standard and is listed for comparison only; glycosylated_entries_sum_w is a mass share, not a shift, and is excluded from the widest-term choice.",
                         "widest_term = the largest term describing this amino acid's uncertainty of the standard (between_fiber_spread and the global completeness share excluded); widest_term_including_global = the largest of everything but the glycosylation share."],
        columns=["profile", "fiber_type", "convention", "amino_acid", "standard"] + term_cols + ["widest_term", "widest_term_including_global"], rows=rows)

    # regime table
    rows = []
    regime = {ft: {"entries_with_value": 0, "sd_ge_median": 0, "sd_ge_half_median": 0, "sd_zero": 0} for ft in FIBER_TYPES}
    for acc in sorted(entries):
        e = entries[acc]
        for ft in FIBER_TYPES:
            m, s, n = e[f"v_{ft}"], e[f"sd_{ft}"], e[f"valid_{ft}"]
            if m <= 0:
                continue
            regime[ft]["entries_with_value"] += 1
            if s == 0:
                regime[ft]["sd_zero"] += 1
            if s >= m:
                regime[ft]["sd_ge_median"] += 1
            if s >= 0.5 * m:
                regime[ft]["sd_ge_half_median"] += 1
            sig = lognormal_sigma(m, s)
            rows.append([acc, e["gene"], ";".join(e["tiers"]), ft, fnum(m), fnum(s), fnum(s / m), fnum(sig), n,
                         fnum(math.exp(sig)), fnum(math.exp(sig / math.sqrt(max(n, 1)))), fnum(e[f"w_{ft}_combined"])])
    files["sd_to_median_ratio.tsv"] = ag.tsv_text(tool="uncertainty.py", header_lines=header_common + ["per entry and fiber type with a value: the published median m and SD s, s/m, the fitted log-normal sigma, n = valid values;",
                         "spread_factor = e^sigma (one-sigma multiplicative spread between fibers); median_factor = e^(sigma/sqrt(n)) (one-sigma multiplicative uncertainty of the median)"],
        columns=["accession", "gene", "tier", "fiber_type", "median", "sd", "sd_over_median", "sigma_log", "n_valid", "spread_factor", "median_factor", "w_combined"], rows=rows)

    summ = common.new_ini()
    summ.add_section("monte_carlo")
    summ["monte_carlo"]["draws"] = str(draws)
    summ["monte_carlo"]["seed"] = str(seed)
    summ["monte_carlo"]["scheme"] = "log-normal per entry, exact moment match to the published median and SD; median_uncertainty uses sigma/sqrt(n_valid)"
    for ft in FIBER_TYPES:
        sec = f"regime.type_{ft}"
        summ.add_section(sec)
        for k, v in regime[ft].items():
            summ[sec][k] = str(v)
    for profile in PROFILES:
        for ft in FIBER_TYPES:
            for conv in CONVENTIONS:
                sec = f"{profile}.type_{ft}.{conv}"
                summ.add_section(sec)
                for term in TERMS_MC:
                    hw = halfwidth[(profile, ft, conv, term)]
                    a_max = max(AA, key=lambda a: hw[a])
                    summ[sec][f"{term}_max_half_width"] = f"{fnum(hw[a_max])} ({a_max})"
                wc = widest_count[(profile, ft, conv)]
                summ[sec]["widest_term_counts"] = "; ".join(f"{k}={v}" for k, v in sorted(wc.items(), key=lambda kv: -kv[1]))
    files["uncertainty_summary.ini"] = common.render_ini(summ, ["GENERATED by uncertainty.py. DO NOT EDIT BY HAND."] + header_common)
    return files, {"draws": draws, "seed": seed}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("uncertainty_check" if args.check else "uncertainty")
    files, _ = build(log)
    if args.check:
        stale = [rel for rel, content in files.items()
                 if not (OUT_DIR / rel).exists() or ag.strip_generated_line((OUT_DIR / rel).read_text(encoding="utf-8")) != ag.strip_generated_line(content)]
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), stale)
            return 1
        log.info("check: uncertainty outputs are current")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for rel, content in files.items():
        (OUT_DIR / rel).write_text(content, encoding="utf-8", newline="\n")
        log.info("wrote %s", (OUT_DIR / rel).relative_to(common.REPO_ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
