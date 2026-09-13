#!/usr/bin/env python3
"""
plots.py — offline stage: figures of the standard, from the tables aggregate.py and
uncertainty.py wrote. Nothing is computed here that is not already in a table; every figure
names the table it was drawn from in its title line.

Reads   outputs/standard/amino_acid_profiles.tsv, uncertainty_intervals.tsv, uncertainty_per_amino_acid.tsv,
        sensitivity_mhc_actin_spread.tsv, profile_differences.tsv, eaa_subset.tsv (if present)
Writes  outputs/standard/plots/<name>.png and .svg
        logs/plots_<UTC>.log

Figures
  profile_<type>.png                 total, contractile, and builders profiles, free convention, with the
                                     median-uncertainty interval as error bars
  profile_fiber_types_total.png      the three fiber types side by side, total profile, free convention
  uncertainty_terms_<type>.png       every uncertainty term per amino acid on the one scale (log axis)
  sensitivity_mhc_actin.png          the profile's range over the measured MHC:actin ratios, types I and IIa
  stress_largest_shift_per_scenario.png  the largest shift per stress scenario (if the stress stage ran)
  eaa_subset.png                     the indispensable amino acids per profile and fiber type (if computed)

Plots are not covered by the currency check: image bytes depend on the plotting library's
version, and the tables they are drawn from are checked instead.
"""

from __future__ import annotations

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

from pipeline import common
from pipeline.mass_fractions import read_tsv_skip_comments
from pipeline import aggregate as ag

OUT_DIR = common.STANDARD_DIR
PLOT_DIR = OUT_DIR / "plots"
AA = ag.AA
FIBER_TYPES = ag.FIBER_TYPES
PROFILES = ag.PROFILES


def _save(fig, name: str, log) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        p = PLOT_DIR / f"{name}.{ext}"
        fig.savefig(p, dpi=150 if ext == "png" else None, bbox_inches="tight", metadata={"Software": None} if ext == "png" else {"Creator": None, "Date": None})
        log.info("wrote %s", p.relative_to(common.REPO_ROOT).as_posix())
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    log = common.make_logger("plots")
    prof = {(r["profile"], r["fiber_type"], r["convention"]): {a: float(r[f"frac_{a}"]) for a in AA}
            for r in read_tsv_skip_comments(OUT_DIR / "amino_acid_profiles.tsv")}
    intervals = {}
    for r in read_tsv_skip_comments(OUT_DIR / "uncertainty_intervals.tsv"):
        intervals[(r["profile"], r["fiber_type"], r["convention"], r["term"], r["amino_acid"])] = (float(r["p2_5"]), float(r["p97_5"]))
    x = range(len(AA))

    # 1. combined vs contractile per fiber type, free convention, median-uncertainty error bars
    for ft in FIBER_TYPES:
        fig, ax = plt.subplots(figsize=(11, 4))
        width = 0.27
        for j, profile in enumerate(PROFILES):
            p = prof[(profile, ft, "free")]
            vals = [100 * p[a] for a in AA]
            lo = [100 * (p[a] - intervals[(profile, ft, "free", "median_uncertainty", a)][0]) for a in AA]
            hi = [100 * (intervals[(profile, ft, "free", "median_uncertainty", a)][1] - p[a]) for a in AA]
            ax.bar([i + (j - 1) * width for i in x], vals, width, yerr=[lo, hi], capsize=2, label=profile)
        ax.set_xticks(list(x), AA)
        ax.set_ylabel("% of amino acid mass (free convention)")
        ax.set_title(f"Fiber type {ft}: total (D58), contractile, builders; error bars = 95 % median uncertainty (amino_acid_profiles.tsv, uncertainty_intervals.tsv)", fontsize=9)
        ax.legend()
        _save(fig, f"profile_{ft}", log)

    # 2. fiber types side by side, combined
    fig, ax = plt.subplots(figsize=(11, 4))
    width = 0.27
    for j, ft in enumerate(FIBER_TYPES):
        p = prof[("total", ft, "free")]
        ax.bar([i + (j - 1) * width for i in x], [100 * p[a] for a in AA], width, label=f"fiber type {ft}")
    ax.set_xticks(list(x), AA)
    ax.set_ylabel("% of amino acid mass (free convention)")
    ax.set_title("Total profile per fiber type (amino_acid_profiles.tsv)", fontsize=9)
    ax.legend()
    _save(fig, "profile_fiber_types_total", log)

    # 3. uncertainty terms on one scale
    terms_rows = read_tsv_skip_comments(OUT_DIR / "uncertainty_per_amino_acid.tsv")
    if terms_rows:
        skip = {"profile", "fiber_type", "convention", "amino_acid", "standard", "widest_term", "widest_term_including_global", "glycosylated_entries_sum_w"}
        term_cols = [c for c in terms_rows[0] if c not in skip]
        for ft in FIBER_TYPES:
            fig, ax = plt.subplots(figsize=(11, 4.5))
            sel = [r for r in terms_rows if r["profile"] == "total" and r["fiber_type"] == ft and r["convention"] == "free"]
            for term in term_cols:
                vals = []
                for a in AA:
                    r = next((r for r in sel if r["amino_acid"] == a), None)
                    v = float(r[term]) if r and r[term] not in ("", "nan") else float("nan")
                    vals.append(max(v, 1e-7))
                ax.plot(list(x), vals, marker="o", markersize=3, linewidth=1, label=term)
            ax.set_yscale("log")
            ax.set_xticks(list(x), AA)
            ax.set_ylabel("max absolute shift of the fraction")
            ax.set_title(f"Fiber type {ft}, total, free convention: every uncertainty term on one scale (uncertainty_per_amino_acid.tsv)", fontsize=9)
            ax.legend(fontsize=7, ncol=2)
            _save(fig, f"uncertainty_terms_{ft}", log)

    # 4. MHC:actin sensitivity
    spread = read_tsv_skip_comments(OUT_DIR / "sensitivity_mhc_actin_spread.tsv")
    if spread:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
        for ax, ft in zip(axes, ag.GEL_TYPES):
            sel = [r for r in spread if r["profile"] == "total" and r["fiber_type"] == ft and r["convention"] == "free"]
            std = [100 * float(next(r for r in sel if r["amino_acid"] == a)["standard"]) for a in AA]
            lo = [100 * float(next(r for r in sel if r["amino_acid"] == a)["min_over_variants"]) for a in AA]
            hi = [100 * float(next(r for r in sel if r["amino_acid"] == a)["max_over_variants"]) for a in AA]
            ax.bar(list(x), std, 0.6, color="lightgray", label="standard")
            ax.errorbar(list(x), std, yerr=[[s - l for s, l in zip(std, lo)], [h - s for s, h in zip(std, hi)]], fmt="none", ecolor="black", capsize=2, label="range over measured MHC:actin ratios")
            ax.set_xticks(list(x), AA)
            ax.set_title(f"Fiber type {ft}, total, free convention", fontsize=9)
        axes[0].set_ylabel("% of amino acid mass")
        axes[0].legend(fontsize=8)
        fig.suptitle("Sensitivity of the profile to the MHC:actin ratio across the measured range (sensitivity_mhc_actin_spread.tsv)", fontsize=9)
        _save(fig, "sensitivity_mhc_actin", log)

    # 5. EAA subset, if computed
    eaa_path = OUT_DIR / "eaa_subset.tsv"
    if eaa_path.exists():
        rows = [r for r in read_tsv_skip_comments(eaa_path) if r["kind"] != "sum"]
        headings = []
        for r in rows:
            if r["heading"] not in headings:
                headings.append(r["heading"])
        fig, ax = plt.subplots(figsize=(10, 4))
        width = 0.09
        combos = [(p, ft) for p in PROFILES for ft in FIBER_TYPES]
        for j, (p, ft) in enumerate(combos):
            vals = [float(next(r["mg_free_amino_acid_per_g_protein"] for r in rows if r["profile"] == p and r["fiber_type"] == ft and r["heading"] == h)) for h in headings]
            ax.bar([i + (j - 4) * width for i in range(len(headings))], vals, width, label=f"{p} {ft}")
        ax.set_xticks(range(len(headings)), headings)
        ax.set_ylabel("mg free amino acid per g protein")
        ax.set_title("Indispensable amino acids per profile and fiber type (eaa_subset.tsv; headings per D62)", fontsize=9)
        ax.legend(fontsize=7, ncol=3)
        _save(fig, "eaa_subset", log)
    else:
        log.info("eaa_subset.tsv not present (D62 transcription pending); no EAA figure")
    return 0


if __name__ == "__main__":
    sys.exit(main())
