#!/usr/bin/env python3
"""
stress.py — offline stage: how far the standard moves when the weights are pushed hard.

The uncertainty stage says how well the standard is known from its source. This stage asks
the question a reader will actually have: what would it take to move it? Every scenario is
"perturb the weights, recompute the profile (rule A1), report the shift"; every magnitude is an
author-set number in config/stress_test_settings.ini [stress]; nothing biological is named —
scenarios are defined by rank, by tier, by band family, or by a factor.

Scenarios (rule A9, D65)
  composition_distance   why it might not move: the free-convention composition of the largest
                         entries against each other and against the standard (per amino acid and
                         as a single distance); no weights perturbed
  convergence            the profile from the top-k entries only (renormalised) against the full
                         standard, for the k listed in config — how many proteins carry the answer
  knockout               remove the largest entry, the top-k entries, or a whole band family;
                         renormalise
  influence              per entry: zero it, double it; the largest shift of any amino acid —
                         sorted, the entries that can move the standard at all
  tier_ratio             set the contractile : builders share to each listed value (measured ~72:28)
  size_tilt              w_i x MW_i^alpha for each listed alpha, renormalised — a systematic bias
                         for or against large proteins (the shape of the iBAQ-vs-other-methods question)
  tpa_weighting          the other quantification convention, exactly: w_i proportional to
                         v_i x N_i (theoretical tryptic peptides, from the digest) instead of
                         v_i x MW_i (D31) — the profile TPA weighting would give from the same values
  random_abuse           every weight x an independent log-normal factor with the listed sigma
                         (given as a factor, e.g. 3 means x/÷ 3 at one sigma), many draws,
                         renormalised — the distribution of shifts if every abundance were wrong by
                         up to that much

Reads   config/mass_fractions_per_entry.tsv, the composition counts, the fetched masses,
        config/stress_test_settings.ini [stress], outputs/mass_fractions/band_families.tsv,
        outputs/digest/theoretical_peptides.tsv, outputs/standard/amino_acid_profiles.tsv
Writes  outputs/standard/stress_shifts.tsv            every scenario x profile x fiber type x convention x amino acid
        outputs/standard/stress_summary_per_scenario.tsv  the largest shift per scenario, in fraction and in percent of the fraction
        outputs/standard/stress_influence_per_entry.tsv   per entry: largest shift from zeroing and from doubling it
        outputs/standard/composition_distance_top_entries.tsv
        outputs/standard/stress_summary.ini
        logs/stress_<UTC>.log
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
BAND_CUTOFF = ag.BAND_CUTOFF


def load_settings() -> dict:
    cp = common.read_ini(common.STRESS_SETTINGS_INI)
    if "stress" not in cp:
        raise SystemExit(f"[STOP] {common.STRESS_SETTINGS_INI.name} has no [stress] section (D65)")
    s = cp["stress"]

    def nums(key, cast=float):
        if key not in s or s[key].strip() in ("", "___"):
            raise SystemExit(f"[STOP] [stress] {key} is not set (D65)")
        return [cast(x) for x in s[key].replace(",", ";").split(";") if x.strip()]

    return {
        "top_k_distance": nums("composition_distance_top_entries", int),
        "convergence_k": nums("convergence_top_k", int),
        "knockout_k": nums("knockout_top_k", int),
        "tier_ratios": nums("contractile_share"),
        "alphas": nums("size_tilt_alpha"),
        "random_factors": nums("random_abuse_factor"),
        "random_draws": int(s.get("random_abuse_draws", "0") or 0),
        "seed": int(s.get("seed", "0") or 0),
    }


def profile_from_weights(w, counts, masses, conv):
    return ag.molar_profile(w, counts, masses, conv)


def shift_rows(scenario, param, profile, ft, base, stressed, rows):
    for conv in CONVENTIONS:
        for a in AA:
            b, s = base[conv][a], stressed[conv][a]
            rows.append([scenario, param, profile, ft, conv, a, fnum(b), fnum(s), fnum(s - b), fnum(100.0 * (s - b) / b) if b > 0 else "nan"])


def max_shift(base, stressed) -> tuple[float, str, float]:
    best, where, pct = 0.0, "", 0.0
    for a in AA:
        d = abs(stressed[a] - base[a])
        if d > best:
            best, where, pct = d, a, (100.0 * d / base[a] if base[a] > 0 else float("nan"))
    return best, where, pct


def build(log) -> dict[str, str]:
    st = load_settings()
    masses = ag.load_masses()
    entries, _ = ag.load_entries()
    counts = ag.load_counts(entries)
    bands = ag.load_band_families()
    peptides = {}
    pep_path = common.DIGEST_DIR / "theoretical_peptides.tsv"
    for r in read_tsv_skip_comments(pep_path):
        peptides[r["accession"]] = int(float(r["n_peptides"]))
    inputs = {"weights": common.MASS_FRACTIONS_TSV, "composition": common.COMPOSITION_TSV, "masses": common.AMINO_ACID_MASSES_INI,
              "settings": common.STRESS_SETTINGS_INI, "band_families": ag.INPUTS["band_families"], "peptides": pep_path}
    hashes = {k: f"{p.relative_to(common.REPO_ROOT).as_posix()} sha256 {sha256_path(p)}" for k, p in inputs.items()}
    header_common = [f"generated = {common.iso_now()}"] + [f"input.{k} = {v}" for k, v in hashes.items()] + [
        "A9 (D65): every scenario perturbs the weights as stated, renormalises, recomputes the profile by rule A1, and reports the shift; magnitudes are the author's numbers in [stress]; nothing biological is named in code.",
        "shift = stressed - standard (fraction); shift_pct = 100 x shift / standard."]

    rows: list[list] = []
    summary_rows: list[list] = []
    influence_rows: list[list] = []
    distance_rows: list[list] = []
    summ_max: dict[str, tuple] = {}

    def record(scenario, param, profile, ft, base, stressed):
        shift_rows(scenario, param, profile, ft, base, stressed, rows)
        for conv in CONVENTIONS:
            m, a, pct = max_shift(base[conv], stressed[conv])
            summary_rows.append([scenario, param, profile, ft, conv, fnum(m), a, fnum(pct)])
            key = (scenario, conv)
            if key not in summ_max or m > summ_max[key][0]:
                summ_max[key] = (m, a, param, profile, ft)

    for profile in PROFILES:
        for ft in FIBER_TYPES:
            w0 = ag.weight_vector(entries, profile, ft)
            base = {conv: profile_from_weights(w0, counts, masses, conv) for conv in CONVENTIONS}
            order = sorted(w0, key=lambda acc: -w0[acc])
            entry_lookup = {acc: entries[acc] for acc in w0}

            # --- composition distance among the largest entries (free convention), no perturbation
            if profile == "total":
                for k in st["top_k_distance"]:
                    top = order[:k]
                    for acc in top:
                        f = counts[acc]["free_frac"]
                        d_std = math.sqrt(sum((f[a] - base["free"][a]) ** 2 for a in AA))
                        distance_rows.append([ft, k, acc, entries[acc]["gene"], common.TIER_NAMES.get(entries[acc]["tiers"][0], entries[acc]["tiers"][0]), fnum(w0[acc]), "standard", fnum(d_std),
                                              max(AA, key=lambda a: abs(f[a] - base["free"][a]))] + [fnum(f[a] - base["free"][a]) for a in AA])
                    for i, x in enumerate(top):
                        for y in top[i + 1:]:
                            fx, fy = counts[x]["free_frac"], counts[y]["free_frac"]
                            d = math.sqrt(sum((fx[a] - fy[a]) ** 2 for a in AA))
                            distance_rows.append([ft, k, x, entries[x]["gene"], common.TIER_NAMES.get(entries[x]["tiers"][0], entries[x]["tiers"][0]), fnum(w0[x]), f"{y} {entries[y]['gene']}", fnum(d),
                                                  max(AA, key=lambda a: abs(fx[a] - fy[a]))] + [fnum(fx[a] - fy[a]) for a in AA])

            # --- convergence: top-k only
            for k in st["convergence_k"]:
                top = order[:k]
                tot = sum(w0[a] for a in top)
                wk = {a: w0[a] / tot for a in top}
                record("convergence", f"top_{k}", profile, ft, base, {conv: profile_from_weights(wk, counts, masses, conv) for conv in CONVENTIONS})

            # --- knockouts: top-k removed; each band family removed
            for k in st["knockout_k"]:
                drop = set(order[:k])
                wk = {a: w for a, w in w0.items() if a not in drop}
                tot = sum(wk.values())
                wk = {a: w / tot for a, w in wk.items()}
                record("knockout", f"top_{k}_removed", profile, ft, base, {conv: profile_from_weights(wk, counts, masses, conv) for conv in CONVENTIONS})
            for band, members in bands.items():
                wk = {a: w for a, w in w0.items() if a not in members}
                tot = sum(wk.values())
                if tot <= 0:
                    continue
                wk = {a: w / tot for a, w in wk.items()}
                record("knockout", f"band_{band}_removed", profile, ft, base, {conv: profile_from_weights(wk, counts, masses, conv) for conv in CONVENTIONS})

            # --- influence per entry: zero it, double it (free convention shift, largest amino acid)
            # done with linear algebra: N = sum_i (w_i/MW_i) c_i; zeroing i removes its term, doubling adds it
            accs = [a for a in order if w0[a] > 0]
            C = np.array([[counts[a]["count"][x] for x in AA] for a in accs], dtype=float)
            per_mol = np.array([w0[a] / counts[a]["mw"] for a in accs])
            N = per_mol @ C
            mfree = np.array([masses[x]["free"] for x in AA])
            P0 = N * mfree / (N * mfree).sum()
            for i, acc in enumerate(accs):
                term = per_mol[i] * C[i]
                for label, Nk in (("zeroed", N - term), ("doubled", N + term)):
                    Pk = Nk * mfree / (Nk * mfree).sum()
                    d = np.abs(Pk - P0)
                    j = int(d.argmax())
                    influence_rows.append([profile, ft, acc, entries[acc]["gene"], common.TIER_NAMES.get(entries[acc]["tiers"][0], entries[acc]["tiers"][0]), fnum(w0[acc]), label, fnum(float(d[j])), AA[j],
                                           fnum(100.0 * float(d[j]) / P0[j]) if P0[j] > 0 else "nan"])

            # --- tier ratio (combined only): set the contractile share, keep within-tier proportions
            if profile == "total":
                share1 = sum(w0[a] for a in w0 if "1" in entry_lookup[a]["tiers"])
                for target in st["tier_ratios"]:
                    wk = {}
                    for a, w in w0.items():
                        in1 = "1" in entry_lookup[a]["tiers"]
                        wk[a] = w * (target / share1 if in1 else (1 - target) / (1 - share1))
                    record("tier_ratio", f"contractile_share_{target:g}", profile, ft, base, {conv: profile_from_weights(wk, counts, masses, conv) for conv in CONVENTIONS})

            # --- size tilt
            for alpha in st["alphas"]:
                wk = {a: w * (counts[a]["mw"] ** alpha) for a, w in w0.items()}
                tot = sum(wk.values())
                wk = {a: w / tot for a, w in wk.items()}
                record("size_tilt", f"alpha_{alpha:g}", profile, ft, base, {conv: profile_from_weights(wk, counts, masses, conv) for conv in CONVENTIONS})

            # --- TPA weighting: w proportional to v x N_peptides instead of v x MW
            wk, missing = {}, 0
            for a, w in w0.items():
                n = peptides.get(a)
                if n is None:
                    missing += 1
                    continue
                v = entry_lookup[a][f"v_{ft}"]
                wk[a] = v * n
            tot = sum(wk.values())
            if tot > 0:
                wk = {a: w / tot for a, w in wk.items()}
                record("tpa_weighting", f"v_x_n_peptides (entries without a digest row: {missing})", profile, ft, base,
                       {conv: profile_from_weights(wk, counts, masses, conv) for conv in CONVENTIONS})

            # --- random abuse
            if st["random_draws"] > 0:
                rng = np.random.default_rng(st["seed"] + 10 * PROFILES.index(profile) + FIBER_TYPES.index(ft))   # deterministic per profile and type
                accs = [a for a in order if w0[a] > 0]
                C = np.array([[counts[a]["count"][x] for x in AA] for a in accs], dtype=float)
                per_mol = np.array([w0[a] / counts[a]["mw"] for a in accs])
                mass = {conv: np.array([masses[x][conv] for x in AA]) for conv in CONVENTIONS}
                for factor in st["random_factors"]:
                    sigma = math.log(factor)
                    F = np.exp(sigma * rng.standard_normal((st["random_draws"], len(accs))))
                    Nk = (F * per_mol) @ C
                    for conv in CONVENTIONS:
                        Pk = Nk * mass[conv]
                        Pk = Pk / Pk.sum(axis=1, keepdims=True)
                        b = np.array([base[conv][a] for a in AA])
                        lo, hi = np.percentile(Pk, [2.5, 97.5], axis=0)
                        worst = np.abs(Pk - b).max(axis=0)
                        for j, a in enumerate(AA):
                            rows.append(["random_abuse", f"factor_{factor:g}_p2.5", profile, ft, conv, a, fnum(b[j]), fnum(lo[j]), fnum(lo[j] - b[j]), fnum(100 * (lo[j] - b[j]) / b[j]) if b[j] > 0 else "nan"])
                            rows.append(["random_abuse", f"factor_{factor:g}_p97.5", profile, ft, conv, a, fnum(b[j]), fnum(hi[j]), fnum(hi[j] - b[j]), fnum(100 * (hi[j] - b[j]) / b[j]) if b[j] > 0 else "nan"])
                        j = int(worst.argmax())
                        summary_rows.append(["random_abuse", f"factor_{factor:g}_worst_of_{st['random_draws']}_draws", profile, ft, conv, fnum(float(worst[j])), AA[j], fnum(100 * float(worst[j]) / b[j])])
                        j95 = int(np.maximum(hi - b, b - lo).argmax())
                        summary_rows.append(["random_abuse", f"factor_{factor:g}_95pct_half_width", profile, ft, conv, fnum(float(max(hi[j95] - b[j95], b[j95] - lo[j95]))), AA[j95], fnum(100 * float(max(hi[j95] - b[j95], b[j95] - lo[j95])) / b[j95])])
                        key = ("random_abuse", conv)
                        m = float(worst[j])
                        if key not in summ_max or m > summ_max[key][0]:
                            summ_max[key] = (m, AA[j], f"factor_{factor:g}", profile, ft)

    influence_rows.sort(key=lambda r: -float(r[7]))
    files = {
        "stress_shifts.tsv": ag.tsv_text(header_common, ["scenario", "parameter", "profile", "fiber_type", "convention", "amino_acid", "standard", "stressed", "shift", "shift_pct"], rows, tool="stress.py"),
        "stress_summary_per_scenario.tsv": ag.tsv_text(header_common + ["per scenario and parameter: the largest |shift| of any amino acid"],
                                                       ["scenario", "parameter", "profile", "fiber_type", "convention", "max_abs_shift", "amino_acid", "max_abs_shift_pct"], summary_rows, tool="stress.py"),
        "stress_influence_per_entry.tsv": ag.tsv_text(header_common + ["per entry with weight: the largest shift of any amino acid (free convention) if the entry were removed (zeroed) or doubled, everything renormalised; sorted by shift"],
                                                      ["profile", "fiber_type", "accession", "gene", "tier", "w", "perturbation", "max_abs_shift", "amino_acid", "max_abs_shift_pct"], influence_rows, tool="stress.py"),
        "composition_distance_top_entries.tsv": ag.tsv_text(header_common + ["free-convention composition of the largest entries of the total: each against the standard and against each other; distance = Euclidean over the twenty fractions; delta_<a> = row entry minus comparator"],
                                                            ["fiber_type", "top_k", "accession", "gene", "tier", "w", "against", "distance", "amino_acid_at_max"] + [f"delta_{a}" for a in AA], distance_rows, tool="stress.py"),
    }
    summ = common.new_ini()
    summ.add_section("settings")
    for k, v in st.items():
        summ["settings"][k] = str(v)
    summ.add_section("largest_shift_per_scenario")
    for (scenario, conv), (m, a, param, profile, ft) in sorted(summ_max.items()):
        summ["largest_shift_per_scenario"][f"{scenario}_{conv}"] = f"{fnum(m)} ({a}; {param}; {profile} {ft})"
    files["stress_summary.ini"] = common.render_ini(summ, ["GENERATED by stress.py. DO NOT EDIT BY HAND."] + header_common)
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    log = common.make_logger("stress_check" if args.check else "stress")
    files = build(log)
    if args.check:
        stale = [rel for rel, content in files.items()
                 if not common.standard_path(rel).exists() or ag.strip_generated_line(common.standard_path(rel).read_text(encoding="utf-8")) != ag.strip_generated_line(content)]
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), stale)
            return 1
        log.info("check: stress outputs are current")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for rel, content in files.items():
        common.write_text_file(common.standard_path(rel), content)
        log.info("wrote %s", common.standard_path(rel).relative_to(common.REPO_ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
