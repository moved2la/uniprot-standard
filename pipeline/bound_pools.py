#!/usr/bin/env python3
"""
bound_pools.py — offline stage: the bound-metabolite adjustment of the skeletal muscle standard.

The protein-only standard counts the amino acids held as protein. Muscle also holds one amino
acid in a stable non-protein form: histidine, as the dipeptide carnosine. This stage adds that
histidine to the protein-bound histidine of the same kilogram of muscle and writes the adjusted
standard beside the protein-only one, which it never changes (D76).

Reads
  config/bound_metabolite_pools.ini                          the author's transcription: pool size per fiber type and
                                                             whole muscle, protein and water per kg muscle, turnover rates
  outputs/standard/amino_acid_g_per_100g_protein_free.tsv    g of free amino acid per 100 g protein, per set and fiber type
  outputs/standard/_calculated_amino_acid_standard.tsv       the protein-only standard (its columns are copied, not recomputed)
  config/fiber_type_mix.ini                                  the shares that make the final column (as the protein-only standard uses them)
  data/pubchem/amino_acid_masses.ini                         the free mass of the pool's amino acid
  data/iupac/amino_acid_symbols.ini                          three-letter -> one-letter

Writes (outputs/standard/, every file with a header naming inputs and hashes)
  _calculated_amino_acid_standard_with_bound_pools.tsv   the ten columns of the protein-only standard with the pool folded
                                                         into total_<fiber type> and standard; contractile and builders
                                                         copied unchanged; percent of amino acid mass, free convention,
                                                         summing to 100
  bound_pool_adjustment_per_amino_acid.tsv               per column and amino acid: protein-only, with pools, the difference
  bound_pool_amounts_per_kg_muscle.tsv                 per fiber type: the pool, the protein-bound amount, their ratio, in grams
                                                         per kg muscle on the pool's basis
  sensitivity_bound_pool_turnover_frame.tsv              the resupply frame: pool x (k_pool / k_protein) over the cited ranges
  sensitivity_bound_pool_basis.tsv                       whole-muscle value on every column vs the single-fiber values
  sensitivity_bound_pool_sex.tsv                         the men's and women's whole-muscle values through the same arithmetic
  sensitivity_bound_pool_spread.tsv                      the single-fiber spread carried through (log-normal fitted to mean and SD)
  bound_pools_summary.ini
  logs/bound_pools_<UTC>.log

Rules (docs/conventions.md, "Aggregation rules")
  A10 Bound pools (D73–D77): per fiber type ft, on the pool's mass basis, F_a = P x g100_ft,a / 100 is the free mass of
      amino acid a from protein per kg muscle (P = protein per kg on that basis; g100 = g per 100 g protein, free
      convention); C = c_ft x n x m / 1000 is the mass of the pool's amino acid per kg (c in mmol/kg, n moles of the
      amino acid per mole of pool, m its free mass). The adjusted profile has F_a for every amino acid but the pool's,
      F_a + C for that one, all renormalised to 100. Contractile and builders columns are protein-only by definition
      and are copied. The final column mixes the adjusted totals by the fiber-type shares, as the protein-only standard does.
      Protein content moves between bases only through the cited water content: P_dry = P_wet / (1 - W/1000).
  The turnover frame (sensitivity): C is replaced by C x k_pool / k_protein with k_pool = ln 2 / (7 x t_half) per day
      and k_protein = 24 x FSR(%/h) / 100 per day, over every cited pair.
  Nothing is named in code: the amino acid, its stoichiometry, every value, every basis, and the fiber-type columns
      a value fills come from the config.
  Green light (the author's rule, 2026-09-17): before anything is computed, every source_id the config names must have
      a hashed file in data/literature/manifest.ini; a missing or unhashed source is a [STOP], not a warning. No
      calculation runs on a source that is not on disk.

Usage
  python pipeline/bound_pools.py            # build
  python pipeline/bound_pools.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from pipeline import common
from pipeline.aggregate import (FIBER_TYPES, load_fiber_type_mix, load_masses, split_list, strip_generated_line,
                                tsv_text)
from pipeline.mass_fractions import read_tsv_skip_comments, fnum, sha256_path

OUT_DIR = common.STANDARD_DIR
AA = list(common.AMINO_ACIDS)
PLACEHOLDER = "___"
COL_ORDER = [(pf, ft) for pf in ("contractile", "builders", "total") for ft in FIBER_TYPES]
TABLE = "_calculated_amino_acid_standard_with_bound_pools.tsv"
Z_95 = 1.959963984540054   # the 2.5 / 97.5 percentile of a standard normal

INPUTS = {
    "pools": common.BOUND_POOLS_INI,
    "g100_free": OUT_DIR / "amino_acid_g_per_100g_protein_free.tsv",
    "standard": OUT_DIR / "_calculated_amino_acid_standard.tsv",
    "fiber_type_mix": common.FIBER_TYPE_MIX_INI,
    "masses": common.AMINO_ACID_MASSES_INI,
    "symbols": common.AMINO_ACID_SYMBOLS_INI,
    "manifest": common.DATA_DIR / "literature" / "manifest.ini",
}


class NotFilled(Exception):
    """A value the calculation needs is still `___`."""


# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------

def _val(cp, section: str, key: str) -> str:
    if section not in cp or key not in cp[section]:
        raise NotFilled(f"[{section}] {key} missing")
    v = cp[section][key].strip()
    if not v or v == PLACEHOLDER:
        raise NotFilled(f"[{section}] {key} is `{PLACEHOLDER}`")
    return v


def _num(cp, section: str, key: str) -> float:
    v = _val(cp, section, key)
    try:
        return float(v)
    except ValueError:
        raise NotFilled(f"[{section}] {key} = {v!r} is not a number")


def _basis(cp, section: str) -> str:
    b = _val(cp, section, "basis").lower()
    if b not in ("dry", "wet"):
        raise SystemExit(f"[STOP] [{section}] basis must be dry or wet, not {b!r}")
    return b


def _cite(cp, section: str) -> str:
    """The source line for a header: source id, location, as_reported, population — whatever is filled."""
    parts = []
    for k in ("source_id", "location", "population", "spread_kind", "as_reported"):
        v = cp[section].get(k, "").strip() if section in cp else ""
        if v and v != PLACEHOLDER:
            parts.append(f"{k}: {v}")
    return "; ".join(parts)


def green_light(cp, manifest_path: Path) -> list[str]:
    """Every source_id named in a filled config line must have at least one hashed file in the manifest.
    Returns the lines it checked (for the header); raises SystemExit naming every source that is not on disk."""
    named: dict[str, list[str]] = {}
    for sec in cp.sections():
        sid = cp[sec].get("source_id", "").strip()
        if sid and sid != PLACEHOLDER:
            named.setdefault(sid, []).append(sec)
    if not named:
        return []
    if not manifest_path.exists():
        raise SystemExit(f"[STOP] {manifest_path.relative_to(common.REPO_ROOT).as_posix()} does not exist: run the literature fetch first (green light)")
    man = common.read_ini(manifest_path)
    hashed = set()
    for msec in man.sections():
        if man[msec].get("sha256", "").strip():
            hashed.add(man[msec].get("source_id", "").strip())
    missing = sorted(sid for sid in named if sid not in hashed)
    if missing:
        detail = "; ".join(f"{sid} (named by {', '.join(named[sid])})" for sid in missing)
        raise SystemExit(f"[STOP] green light refused: no hashed file in the literature manifest for {detail}. "
                         f"Fetch or hand-obtain the file, run `python run.py mass-fractions --stop-after fetch_literature`, and rerun.")
    return [f"green light: {sid} -> hashed file in manifest ({', '.join(named[sid])})" for sid in sorted(named)]


def load_pool(cp, symbols_path: Path, masses: dict) -> dict:
    """The one pool carried (D75): its amino acid resolved to a letter, its stoichiometry, its free mass."""
    pools = [s for s in cp.sections() if s.startswith("pool.")]
    if len(pools) != 1:
        raise SystemExit(f"[STOP] exactly one [pool.<name>] section is expected (D75); found {pools}")
    sec = pools[0]
    name = sec[len("pool."):]
    three = _val(cp, sec, "amino_acid")
    sym = common.read_ini(symbols_path)
    three_to_one = {sym[a]["three_letter"].strip().lower(): a for a in sym.sections() if a in AA}
    if three.lower() not in three_to_one:
        raise SystemExit(f"[STOP] [{sec}] amino_acid {three!r} is not a three-letter symbol in the IUPAC table")
    letter = three_to_one[three.lower()]
    n = _num(cp, sec, "moles_amino_acid_per_mole")
    return {"name": name, "section": sec, "three": three, "letter": letter, "n_per_mole": n,
            "free_mass": masses[letter]["free"], "cite": _cite(cp, sec)}


def protein_per_kg(cp, basis: str) -> tuple[float, str]:
    """Protein content in g per kg muscle on `basis`, converted through the cited water content if the
    source's basis differs. Returns the value and a line saying how it was obtained."""
    p = _num(cp, "protein_content", "value")
    pb = _basis(cp, "protein_content")
    if pb == basis:
        return p, f"protein {p:g} g/kg {basis} muscle as cited ({_cite(cp, 'protein_content')})"
    w = _num(cp, "water_content", "value")          # g water per kg wet
    if not 0 < w < 1000:
        raise SystemExit(f"[STOP] [water_content] value {w} must be between 0 and 1000 g per kg wet")
    if pb == "wet" and basis == "dry":
        conv = p / (1.0 - w / 1000.0)
        how = f"P_dry = P_wet / (1 - W/1000) = {p:g} / (1 - {w:g}/1000)"
    else:
        conv = p * (1.0 - w / 1000.0)
        how = f"P_wet = P_dry x (1 - W/1000) = {p:g} x (1 - {w:g}/1000)"
    return conv, (f"protein {conv:.4f} g/kg {basis} muscle: {how} "
                  f"(protein: {_cite(cp, 'protein_content')}; water: {_cite(cp, 'water_content')})")


def single_fiber_values(cp, pool_name: str) -> dict[str, dict]:
    """Per fiber type: the single-fiber pool value, basis, spread, and the section it came from (via applies_to).
    Sections are [<pool>.single_fiber.<anything>]; which fiber-type columns each fills is its applies_to list."""
    out: dict[str, dict] = {}
    for sec in cp.sections():
        if not sec.startswith(f"{pool_name}.single_fiber."):
            continue
        for ft in split_list(_val(cp, sec, "applies_to")):
            if ft not in FIBER_TYPES:
                raise SystemExit(f"[STOP] [{sec}] applies_to names {ft!r}, not one of {FIBER_TYPES}")
            if ft in out:
                raise SystemExit(f"[STOP] fiber type {ft} is filled by two sections: {out[ft]['section']} and {sec}")
            out[ft] = {"section": sec, "value": _num(cp, sec, "value"), "basis": _basis(cp, sec),
                       "spread": cp[sec].get("spread", "").strip(), "cite": _cite(cp, sec)}
    missing = [ft for ft in FIBER_TYPES if ft not in out]
    if missing:
        raise NotFilled(f"no single-fiber section applies to fiber type(s) {missing}")
    bases = {v["basis"] for v in out.values()}
    if len(bases) != 1:
        raise SystemExit(f"[STOP] single-fiber values are on different bases: {bases}")
    return out


# ----------------------------------------------------------------------------
# inputs from the protein-only standard
# ----------------------------------------------------------------------------

def load_g100() -> dict[tuple[str, str], dict[str, float]]:
    """(set, fiber type) -> g of free amino acid per 100 g protein, per amino acid."""
    rows = read_tsv_skip_comments(INPUTS["g100_free"])
    out: dict[tuple[str, str], dict[str, float]] = {}
    body = [r for r in rows if r["amino_acid"] in AA]
    if len(body) != 20:
        raise SystemExit(f"[STOP] {INPUTS['g100_free'].name}: expected twenty amino acid rows, found {len(body)}")
    for pf, ft in COL_ORDER:
        col = f"{pf}_{ft}"
        out[(pf, ft)] = {r["amino_acid"]: float(r[col]) for r in body}
    return out


def load_protein_only_standard() -> tuple[list[dict], list[str]]:
    rows = read_tsv_skip_comments(INPUTS["standard"])
    header = [ln for ln in INPUTS["standard"].read_text(encoding="utf-8").splitlines() if ln.startswith("#")]
    return rows, header


# ----------------------------------------------------------------------------
# the arithmetic (A10)
# ----------------------------------------------------------------------------

def adjusted_profile(g100: dict[str, float], protein_g_per_kg: float, pool_mmol_per_kg: float,
                     letter: str, n_per_mole: float, free_mass: float) -> dict:
    """A10 for one column. Returns the masses per kg (protein-derived F_a, the pool C), the ratio C / F_letter,
    and the adjusted fractions summing to one."""
    f = {a: protein_g_per_kg * g100[a] / 100.0 for a in AA}
    c = pool_mmol_per_kg / 1000.0 * n_per_mole * free_mass
    adj = dict(f)
    adj[letter] = f[letter] + c
    tot = sum(adj.values())
    return {"F": f, "C": c, "ratio": (c / f[letter] if f[letter] > 0 else float("nan")),
            "fractions": {a: adj[a] / tot for a in AA}, "protein_free_total": sum(f.values())}


def lognormal_bounds(mean: float, sd: float) -> tuple[float, float]:
    """2.5 and 97.5 percentiles of a log-normal with the given mean and SD (moment match; the D63 pattern)."""
    if mean <= 0:
        return float("nan"), float("nan")
    s2 = math.log(1.0 + (sd / mean) ** 2)
    mu = math.log(mean) - s2 / 2.0
    return math.exp(mu - Z_95 * math.sqrt(s2)), math.exp(mu + Z_95 * math.sqrt(s2))


def k_per_day_from_half_life_weeks(t_half_weeks: float) -> float:
    return math.log(2.0) / (7.0 * t_half_weeks)


def k_per_day_from_fsr_percent_per_hour(fsr: float) -> float:
    return 24.0 * fsr / 100.0


# ----------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------

def build(log) -> dict[str, str]:
    cp = common.read_ini(INPUTS["pools"])
    light = green_light(cp, INPUTS["manifest"])
    for ln in light:
        log.info(ln)
    hashes = {k: f"{p.relative_to(common.REPO_ROOT).as_posix()} sha256 {sha256_path(p)}" for k, p in INPUTS.items() if p.exists()}
    meta = dict(cp["meta"]) if "meta" in cp else {}
    decisions = meta.get("decisions", "").strip()
    label = meta.get("version_label", "").strip()
    header_common = [f"generated = {common.iso_now()}"] + [f"input.{k} = {v}" for k, v in hashes.items()] + [
        f"version label: ref = {label} ({decisions})",
    ] + light + [
        "A10: on the pool's mass basis, F_a = P x g100_a / 100 = free mass of amino acid a from protein per kg muscle; "
        "C = c x n x m / 1000 = mass of the pool's amino acid per kg; that amino acid becomes F + C, every other stays F, all renormalised to 100. "
        "Contractile and builders are protein-only by definition and copied; the final column mixes the adjusted totals by the fiber-type shares.",
    ]
    files: dict[str, str] = {}
    summ = common.new_ini()
    summ.add_section("status")
    summ.add_section("inputs")
    for k, v in hashes.items():
        summ["inputs"][k] = v

    # ---- the pool and the protein-only standard
    masses = load_masses()
    g100 = load_g100()
    std_rows, std_header = load_protein_only_standard()
    mix = load_fiber_type_mix()
    try:
        pool = load_pool(cp, INPUTS["symbols"], masses)
        sf = single_fiber_values(cp, pool["name"])
        basis = next(iter(sf.values()))["basis"]
        protein, protein_line = protein_per_kg(cp, basis)
    except NotFilled as e:
        summ["status"]["adjusted_standard"] = f"NOT computed: config/bound_metabolite_pools.ini is not filled in — {e}"
        log.warning("adjusted standard NOT written: %s", e)
        files["bound_pools_summary.ini"] = common.render_ini(summ, ["GENERATED by bound_pools.py. DO NOT EDIT BY HAND."] + header_common)
        return files
    letter = pool["letter"]
    log.info("pool %s -> %s (%s), %g mol per mol, free mass %g; basis %s; %s", pool["name"], pool["three"], letter,
             pool["n_per_mole"], pool["free_mass"], basis, protein_line)

    # ---- 1. the adjusted standard
    adj_total: dict[str, dict] = {}
    for ft in FIBER_TYPES:
        adj_total[ft] = adjusted_profile(g100[("total", ft)], protein, sf[ft]["value"], letter, pool["n_per_mole"], pool["free_mass"])
        log.info("fiber type %s: pool %g mmol/kg %s -> %s %.4f g/kg; protein-bound %s %.4f g/kg; ratio %.4f; %s %.4f -> %.4f %%",
                 ft, sf[ft]["value"], basis, letter, adj_total[ft]["C"], letter, adj_total[ft]["F"][letter], adj_total[ft]["ratio"],
                 letter, 100 * g100[("total", ft)][letter] / sum(g100[("total", ft)].values()), 100 * adj_total[ft]["fractions"][letter])
    final = None
    if mix is not None:
        final = {a: sum(mix["shares"][ft] * adj_total[ft]["fractions"][a] for ft in FIBER_TYPES) for a in AA}

    protein_only = {r["amino_acid"]: r for r in std_rows}
    cols = ["amino_acid", "three_letter", "name"] + [f"{pf}_{ft}" for pf, ft in COL_ORDER] + ["standard"]
    body = []
    for a in AA:
        row = [a, masses[a]["three"], masses[a]["name"]]
        for pf, ft in COL_ORDER:
            row.append(f"{100 * adj_total[ft]['fractions'][a]:.4f}" if pf == "total" else protein_only[a][f"{pf}_{ft}"])
        row.append(f"{100 * final[a]:.4f}" if final else "")
        body.append(row)
    sums = ["sum", "", ""]
    for pf, ft in COL_ORDER:
        sums.append(f"{100 * sum(adj_total[ft]['fractions'].values()):.4f}" if pf == "total" else protein_only.get("sum", {}).get(f"{pf}_{ft}", ""))
    sums.append(f"{100 * sum(final.values()):.4f}" if final else "")
    body.append(sums)
    pool_lines = [f"pool {pool['name']} -> {pool['three']}, {pool['n_per_mole']:g} mol per mol, free mass {pool['free_mass']} g/mol from PubChem ({pool['cite']})"]
    for ft in FIBER_TYPES:
        pool_lines.append(f"{pool['name']} fiber type {ft} = {sf[ft]['value']:g} mmol/kg {basis} muscle ({sf[ft]['cite']}) -> "
                          f"{letter} from the pool {adj_total[ft]['C']:.4f} g/kg vs {letter} from protein {adj_total[ft]['F'][letter]:.4f} g/kg; ratio {adj_total[ft]['ratio']:.4f}")
    mix_line = (f"standard = " + " + ".join(f"{mix['shares'][ft]:.4f} x total_{ft}" for ft in FIBER_TYPES) + f" (config/fiber_type_mix.ini, as the protein-only standard mixes)"
                if mix else "standard = BLANK until config/fiber_type_mix.ini is filled")
    header = header_common + [
        f"THE CALCULATED AMINO ACID STANDARD of human skeletal muscle protein WITH BOUND METABOLITE POOLS — percent of total amino acid mass, free convention, summing to 100.",
        "the protein-only standard (_calculated_amino_acid_standard.tsv) is unchanged; here total_<fiber type> and standard include the amino acid held as the pool; contractile and builders are copied from it.",
        protein_line, *pool_lines, mix_line,
        "the per-amino-acid difference is in bound_pool_adjustment_per_amino_acid.tsv; the sensitivities in sensitivity_bound_pool_*.tsv",
    ]
    files[TABLE] = tsv_text(header, cols, body, tool="bound_pools.py")

    # ---- 2. the difference
    rows = []
    for pf, ft in COL_ORDER:
        if pf != "total":
            continue
        tot = sum(g100[("total", ft)].values())
        for a in AA:
            before = g100[("total", ft)][a] / tot
            after = adj_total[ft]["fractions"][a]
            rows.append([f"total_{ft}", a, fnum(before), fnum(after), fnum(after - before), fnum((after - before) / before if before > 0 else float("nan"))])
    if final:
        for a in AA:
            before = float(protein_only[a]["standard"]) / 100.0 if protein_only[a]["standard"] else float("nan")
            after = final[a]
            rows.append(["standard", a, fnum(before), fnum(after), fnum(after - before), fnum((after - before) / before if before > 0 else float("nan"))])
    files["bound_pool_adjustment_per_amino_acid.tsv"] = tsv_text(
        header_common + ["per column and amino acid: the protein-only fraction, the fraction with the pool folded in, their difference, and the difference as a share of the protein-only fraction (free convention)"],
        ["column", "amino_acid", "protein_only", "with_bound_pools", "difference", "relative_difference"], rows, tool="bound_pools.py")

    # ---- 3. the amounts per kg
    rows = []
    for ft in FIBER_TYPES:
        r = adj_total[ft]
        rows.append([ft, basis, fnum(sf[ft]["value"]), sf[ft]["section"], fnum(r["C"]), fnum(protein), fnum(r["protein_free_total"]),
                     fnum(r["F"][letter]), fnum(r["ratio"]), fnum(100 * g100[("total", ft)][letter] / sum(g100[("total", ft)].values())), fnum(100 * r["fractions"][letter])])
    files["bound_pool_amounts_per_kg_muscle.tsv"] = tsv_text(
        header_common + [f"per fiber type, per kg muscle on the pool's basis: the pool in mmol and as grams of {pool['three']}; protein per kg; total free amino acid mass from that protein; {pool['three']} from protein; the ratio pool : protein"],
        ["fiber_type", "basis", "pool_mmol_per_kg", "pool_section", f"pool_{pool['three']}_g_per_kg", "protein_g_per_kg", "free_amino_acid_from_protein_g_per_kg",
         f"protein_{pool['three']}_g_per_kg", "ratio_pool_to_protein", f"protein_only_{pool['three']}_percent", f"with_pool_{pool['three']}_percent"], rows, tool="bound_pools.py")

    # ---- 4. sensitivities (each skipped, and said so, when its inputs are not filled)
    # 4a. turnover frame
    try:
        halves = []
        for sec in cp.sections():
            if sec.startswith(f"turnover.{pool['name']}."):
                halves.append((sec, _num(cp, sec, "t_half_weeks")))
        fsrs = [(sec, _num(cp, sec, "value")) for sec in cp.sections() if sec.startswith("turnover.muscle_protein.")]
        if not halves or not fsrs:
            raise NotFilled("no turnover sections")
        rows = []
        for ft in FIBER_TYPES:
            base = adj_total[ft]
            for hsec, th in halves:
                k_pool = k_per_day_from_half_life_weeks(th)
                for fsec, fsr in fsrs:
                    k_prot = k_per_day_from_fsr_percent_per_hour(fsr)
                    factor = k_pool / k_prot
                    r = adjusted_profile(g100[("total", ft)], protein, sf[ft]["value"] * factor, letter, pool["n_per_mole"], pool["free_mass"])
                    rows.append([ft, hsec, fnum(th), fnum(k_pool), fsec, fnum(fsr), fnum(k_prot), fnum(factor), fnum(base["ratio"]), fnum(r["ratio"]),
                                 fnum(100 * base["fractions"][letter]), fnum(100 * r["fractions"][letter])])
        files["sensitivity_bound_pool_turnover_frame.tsv"] = tsv_text(
            header_common + ["the resupply (turnover-weighted) frame, a sensitivity only (D73): the pool counted as C x k_pool / k_protein, k_pool = ln 2 / (7 x t_half) per day, k_protein = 24 x FSR / 100 per day, over every cited pair;",
                             "ratio_pool_frame = the inventory ratio of the standard; ratio_turnover_frame = the same ratio in this frame; the percent columns are the pool's amino acid in total_<fiber type> under each frame"],
            ["fiber_type", "pool_half_life_section", "t_half_weeks", "k_pool_per_day", "protein_fsr_section", "fsr_percent_per_hour", "k_protein_per_day", "k_pool_over_k_protein",
             "ratio_pool_frame", "ratio_turnover_frame", f"{pool['three']}_percent_pool_frame", f"{pool['three']}_percent_turnover_frame"], rows, tool="bound_pools.py")
        summ["status"]["sensitivity_turnover_frame"] = f"written: {len(halves)} half-lives x {len(fsrs)} protein rates"
    except NotFilled as e:
        summ["status"]["sensitivity_turnover_frame"] = f"NOT computed: {e}"
        log.warning("turnover-frame sensitivity NOT written: %s", e)

    # 4b. basis (whole-muscle value on every column) and 4c. sex
    def whole_muscle_profiles(section: str) -> dict[str, dict]:
        wm_val = _num(cp, section, "value")
        wm_basis = _basis(cp, section)
        p_wm, _ = protein_per_kg(cp, wm_basis)
        return {ft: adjusted_profile(g100[("total", ft)], p_wm, wm_val, letter, pool["n_per_mole"], pool["free_mass"]) for ft in FIBER_TYPES} | {"_value": wm_val, "_basis": wm_basis, "_protein": p_wm}

    try:
        wm = whole_muscle_profiles(f"{pool['name']}.whole_muscle")
        rows = []
        for ft in FIBER_TYPES:
            rows.append([f"total_{ft}", "single_fiber", basis, fnum(sf[ft]["value"]), fnum(protein), fnum(adj_total[ft]["ratio"]), fnum(100 * adj_total[ft]["fractions"][letter]),
                         "whole_muscle", wm["_basis"], fnum(wm["_value"]), fnum(wm["_protein"]), fnum(wm[ft]["ratio"]), fnum(100 * wm[ft]["fractions"][letter]),
                         fnum(100 * (wm[ft]["fractions"][letter] - adj_total[ft]["fractions"][letter]))])
        if final:
            fin_wm = sum(mix["shares"][ft] * wm[ft]["fractions"][letter] for ft in FIBER_TYPES)
            rows.append(["standard", "single_fiber", basis, "", fnum(protein), "", fnum(100 * final[letter]),
                         "whole_muscle", wm["_basis"], fnum(wm["_value"]), fnum(wm["_protein"]), "", fnum(100 * fin_wm), fnum(100 * (fin_wm - final[letter]))])
        files["sensitivity_bound_pool_basis.tsv"] = tsv_text(
            header_common + ["two measurement options for the pool, no verdict (the D54 pattern): the single-fiber values that fill the columns (D74) against one whole-muscle value applied to every column; the pool's amino acid in percent under each, and the difference"],
            ["column", "option_a", "basis_a", "pool_mmol_per_kg_a", "protein_g_per_kg_a", "ratio_a", f"{pool['three']}_percent_a",
             "option_b", "basis_b", "pool_mmol_per_kg_b", "protein_g_per_kg_b", "ratio_b", f"{pool['three']}_percent_b", "percent_b_minus_a"], rows, tool="bound_pools.py")
        summ["status"]["sensitivity_basis"] = "written"
    except NotFilled as e:
        summ["status"]["sensitivity_basis"] = f"NOT computed: {e}"
        log.warning("basis sensitivity NOT written: %s", e)
        wm = None

    try:
        if wm is None:
            raise NotFilled("whole-muscle value not filled")
        groups = {}
        for sec in cp.sections():
            if sec.startswith(f"{pool['name']}.whole_muscle."):
                groups[sec.rsplit(".", 1)[1]] = whole_muscle_profiles(sec)
        if not groups:
            raise NotFilled("no whole-muscle subgroup sections")
        rows = []
        for gname, g in groups.items():
            for ft in FIBER_TYPES:
                rows.append([gname, f"total_{ft}", g["_basis"], fnum(g["_value"]), fnum(g[ft]["ratio"]), fnum(100 * g[ft]["fractions"][letter]),
                             fnum(wm["_value"]), fnum(100 * wm[ft]["fractions"][letter]), fnum(100 * (g[ft]["fractions"][letter] - wm[ft]["fractions"][letter]))])
            if final:
                fin_g = sum(mix["shares"][ft] * g[ft]["fractions"][letter] for ft in FIBER_TYPES)
                fin_wm = sum(mix["shares"][ft] * wm[ft]["fractions"][letter] for ft in FIBER_TYPES)
                rows.append([gname, "standard", g["_basis"], fnum(g["_value"]), "", fnum(100 * fin_g), fnum(wm["_value"]), fnum(100 * fin_wm), fnum(100 * (fin_g - fin_wm))])
        files["sensitivity_bound_pool_sex.tsv"] = tsv_text(
            header_common + ["each whole-muscle subgroup's value through the same arithmetic, against the all-subjects whole-muscle value; the pool's amino acid in percent and the difference"],
            ["group", "column", "basis", "pool_mmol_per_kg", "ratio", f"{pool['three']}_percent", "pool_mmol_per_kg_all_subjects", f"{pool['three']}_percent_all_subjects", "percent_group_minus_all"], rows, tool="bound_pools.py")
        summ["status"]["sensitivity_sex"] = f"written: {', '.join(groups)}"
    except NotFilled as e:
        summ["status"]["sensitivity_sex"] = f"NOT computed: {e}"
        log.warning("sex sensitivity NOT written: %s", e)

    # 4d. the single-fiber spread carried through
    try:
        rows = []
        for ft in FIBER_TYPES:
            s = sf[ft]
            if not s["spread"] or s["spread"] == PLACEHOLDER:
                raise NotFilled(f"[{s['section']}] spread")
            lo, hi = lognormal_bounds(s["value"], float(s["spread"]))
            r_lo = adjusted_profile(g100[("total", ft)], protein, lo, letter, pool["n_per_mole"], pool["free_mass"])
            r_hi = adjusted_profile(g100[("total", ft)], protein, hi, letter, pool["n_per_mole"], pool["free_mass"])
            rows.append([ft, s["section"], fnum(s["value"]), s["spread"], cp[s["section"]].get("spread_kind", ""), fnum(lo), fnum(hi),
                         fnum(100 * adj_total[ft]["fractions"][letter]), fnum(100 * r_lo["fractions"][letter]), fnum(100 * r_hi["fractions"][letter])])
        files["sensitivity_bound_pool_spread.tsv"] = tsv_text(
            header_common + ["the published single-fiber spread carried through: a log-normal fitted to the mean and SD (the D63 pattern; a between-fiber spread, not the uncertainty of the mean), its 2.5 and 97.5 percentiles, and the pool's amino acid in percent at each"],
            ["fiber_type", "section", "mean_mmol_per_kg", "spread", "spread_kind", "p2_5_mmol_per_kg", "p97_5_mmol_per_kg", f"{pool['three']}_percent_at_mean", f"{pool['three']}_percent_at_p2_5", f"{pool['three']}_percent_at_p97_5"], rows, tool="bound_pools.py")
        summ["status"]["sensitivity_spread"] = "written"
    except NotFilled as e:
        summ["status"]["sensitivity_spread"] = f"NOT computed: {e}"
        log.warning("spread sensitivity NOT written: %s", e)

    # ---- summary
    summ["status"]["adjusted_standard"] = f"written: outputs/standard/{TABLE}"
    summ.add_section("the_pool")
    summ["the_pool"]["name"] = pool["name"]
    summ["the_pool"]["amino_acid"] = f"{pool['three']} ({letter}); {pool['n_per_mole']:g} mol per mol; free mass {pool['free_mass']} g/mol"
    summ["the_pool"]["basis"] = basis
    summ["the_pool"]["protein_per_kg"] = protein_line
    for ft in FIBER_TYPES:
        r = adj_total[ft]
        summ["the_pool"][f"fiber_type_{ft}"] = (f"pool {sf[ft]['value']:g} mmol/kg -> {letter} {r['C']:.4f} g/kg; protein-bound {letter} {r['F'][letter]:.4f} g/kg; ratio {r['ratio']:.4f}; "
                                               f"{letter} {100 * g100[('total', ft)][letter] / sum(g100[('total', ft)].values()):.4f} -> {100 * r['fractions'][letter]:.4f} %")
    if final:
        summ["the_pool"]["standard"] = f"{letter} {protein_only[letter]['standard']} -> {100 * final[letter]:.4f} %"
    files["bound_pools_summary.ini"] = common.render_ini(summ, ["GENERATED by bound_pools.py. DO NOT EDIT BY HAND."] + header_common)
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("bound_pools_check" if args.check else "bound_pools")
    files = build(log)
    if args.check:
        stale = [rel for rel, content in files.items()
                 if not (OUT_DIR / rel).exists() or strip_generated_line((OUT_DIR / rel).read_text(encoding="utf-8")) != strip_generated_line(content)]
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), stale)
            return 1
        log.info("check: bound-pool outputs are current")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for rel, content in files.items():
        (OUT_DIR / rel).write_text(content, encoding="utf-8", newline="\n")
        log.info("wrote %s", (OUT_DIR / rel).relative_to(common.REPO_ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
