#!/usr/bin/env python3
"""
aggregate.py — offline stage: the standard. The amino acid profile of human skeletal muscle
protein per fiber type, in both mass conventions, from the per-protein compositions (Layer A)
and the per-protein weights (Layer B), and nothing else.

Reads
  config/mass_fractions_per_entry.tsv                        v, sd, valid, mw, weights per entry and fiber type (D61)
  outputs/composition/amino_acid_composition_per_protein.tsv residue counts per entry, segment_set = master
  data/pubchem/amino_acid_masses.ini                         free and residue mass per letter; three-letter symbol; name
  data/iupac/amino_acid_symbols.ini                          three-letter -> one-letter (for the EAA config)
  config/aggregation_decisions.ini                           the indispensable amino acids and their groupings (D62)
  outputs/mass_fractions/classical_check_carroll_2004.tsv    the MHC:actin ratio as measured by each method
  outputs/mass_fractions/band_families.tsv                   the entries in each band's family (cutoff 2, D52)
  outputs/mass_fractions/family_bounds.tsv                   within-tier family bounds (D55)
  outputs/isoform_deltas.tsv, outputs/composition/processing_mass_deltas.tsv,
  outputs/composition/ptm_mass_deltas.tsv                    per-entry bounds
  outputs/mass_fractions/dataset_rows_outside_pool.tsv       the completeness list (molar shares)

Writes (outputs/standard/, every file with a header naming inputs and hashes)
  amino_acid_profiles.tsv              one row per profile x fiber type x convention: twenty fractions summing
                                       to one (A1), the fraction-mixing check (A1c), entries with weight
  amino_acid_profiles_free_wide.tsv    g per 100 g protein, one row per amino acid, one column per profile and
  amino_acid_profiles_residue_wide.tsv fiber type — the human-readable tables
  profile_differences.tsv              per amino acid: combined - contractile per fiber type; between fiber
                                       types; the fiber-type bracket (max - min over the three types, A5)
  sensitivity_mhc_actin_profiles.tsv   the profile at each measured MHC:actin ratio, both single-band ways (A6)
  sensitivity_mhc_actin_spread.tsv     per amino acid: spread and max absolute shift across those profiles
  bounds_after_weighting.tsv           isoform, processing, PTM, family, shared-row bounds x weight, per profile,
                                       fiber type, and amino acid: largest single entry and sum over entries
  completeness_sensitivity.tsv         the bound on the profile from the genes outside the pool, by rank (A7)
  eaa_subset.tsv                       the indispensable amino acids (D62) per profile and fiber type
  standard_summary.ini
  logs/aggregate_<UTC>.log

Rules (docs/conventions.md, "Aggregation rules")
  A0 every input file is hashed into every output header.
  A1 the profile is the molar mixture: n_a = sum_i (w_i / MW_i) * count_{i,a} = sum_i v_i * count_{i,a} / D,
     then p_a = n_a * m_a / sum_b n_b * m_b with m the residue mass (residue convention) or the free mass
     (free convention). A1c: the fraction mixture sum_i w_i * f_{i,a}, renormalised, is written beside it
     as a check; the difference is disclosed.
  A2 sets: combined = every weighted entry (D58); contractile = tier 1. No other set is named.
  A5 the bracket: max - min of each amino acid over the three pure fiber types bounds every mix.
  A6 MHC:actin sensitivity (D54b as reframed): at each ratio the three methods measured, the band family is
     rescaled two ways (MHC family scaled, actin fixed; actin family scaled, MHC fixed), all weights
     renormalised; no method is a reference; types I and IIa only; no IIx factor.
  A7 completeness: a mass share s of unknown composition moves any fraction by at most s; the molar share of
     the outside genes is converted with the pool's molar-mean molecular weight, stated as the assumption.
  Nothing is named in code: the band anchors, the EAA list, and the sources come from config and from files.

Usage
  python pipeline/aggregate.py            # build
  python pipeline/aggregate.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

from pipeline import common
from pipeline.mass_fractions import read_tsv_skip_comments, fnum, sha256_path

OUT_DIR = common.STANDARD_DIR
MF_DIR = common.MASS_FRACTIONS_DIR
FIBER_TYPES = ("I", "IIa", "IIx")
PROFILES = ("combined", "contractile")
CONVENTIONS = ("free", "residue")
AA = list(common.AMINO_ACIDS)
GEL_TYPES = ("I", "IIa")            # the fiber types the gel measured; no IIx factor is invented (A6)
BAND_CUTOFF = "2"                   # the D52 cutoff
COMPLETENESS_RANKS = (1, 5, 10, 30)

INPUTS = {
    "weights": common.MASS_FRACTIONS_TSV,
    "composition": common.COMPOSITION_TSV,
    "masses": common.AMINO_ACID_MASSES_INI,
    "symbols": common.AMINO_ACID_SYMBOLS_INI,
    "eaa": common.AGGREGATION_DECISIONS_INI,
    "mhc_actin": MF_DIR / "classical_check_carroll_2004.tsv",
    "band_families": MF_DIR / "band_families.tsv",
    "family_bounds": MF_DIR / "family_bounds.tsv",
    "isoform_deltas": common.OUTPUTS_DIR / "isoform_deltas.tsv",
    "processing_deltas": common.COMPOSITION_DIR / "processing_mass_deltas.tsv",
    "ptm_deltas": common.COMPOSITION_DIR / "ptm_mass_deltas.tsv",
    "outside_pool": MF_DIR / "dataset_rows_outside_pool.tsv",
}


# ----------------------------------------------------------------------------
# inputs
# ----------------------------------------------------------------------------

def load_masses() -> dict[str, dict]:
    cp = common.read_ini(INPUTS["masses"])
    out = {}
    for a in AA:
        if a not in cp:
            raise SystemExit(f"[STOP] amino acid {a} missing from {INPUTS['masses'].name}")
        s = cp[a]
        out[a] = {"free": float(s["free_mass"]), "residue": float(s["residue_mass"]),
                  "three": s.get("three_letter", ""), "name": s.get("trivial_name", "")}
    return out


def load_entries() -> tuple[dict[str, dict], list[str]]:
    """Per entry: tier list, mw, and per fiber type v, sd, valid, the combined weight, the tier-1
    weight, and the D48 shared-row excess (recovered from w_high). Plus the header lines of the
    weights file, so its provenance can be carried."""
    header = [ln for ln in INPUTS["weights"].read_text(encoding="utf-8").splitlines() if ln.startswith("#")]
    rows = read_tsv_skip_comments(INPUTS["weights"])
    entries = {}
    for r in rows:
        e = {"gene": r["gene"], "tiers": r["tier"].split(";"), "mw": float(r["mw_master"]), "match_rule": r["match_rule"]}
        for ft in FIBER_TYPES:
            e[f"v_{ft}"] = float(r[f"v_{ft}"])
            e[f"sd_{ft}"] = float(r[f"sd_{ft}"])
            e[f"valid_{ft}"] = int(float(r[f"valid_{ft}"])) if r[f"valid_{ft}"] else 0
            e[f"w_{ft}_combined"] = float(r[f"w_{ft}_combined"])
            e[f"w_{ft}_combined_high"] = float(r[f"w_{ft}_combined_high"])
            t1 = r.get(f"w_{ft}_tier1", "")
            e[f"w_{ft}_tier1"] = float(t1) if t1 not in ("", None) else None
        entries[r["accession"]] = e
    if not entries:
        raise SystemExit(f"[STOP] no rows in {INPUTS['weights']}")
    for ft in FIBER_TYPES:
        s = sum(e[f"w_{ft}_combined"] for e in entries.values())
        if abs(s - 1.0) > 1e-6:
            raise SystemExit(f"[STOP] combined weights for type {ft} sum to {s}, not 1 (B8)")
    # the denominator D per fiber type = sum v * mw; the shared-row excess per entry from w_high (B5)
    for ft in FIBER_TYPES:
        d = sum(e[f"v_{ft}"] * e["mw"] for e in entries.values())
        for e in entries.values():
            excess = e[f"w_{ft}_combined_high"] * d / e["mw"] - e[f"v_{ft}"] - e[f"sd_{ft}"] if e["mw"] > 0 else 0.0
            e[f"shared_excess_{ft}"] = max(0.0, excess)
        for e in entries.values():
            e[f"D_{ft}"] = d
    return entries, header


def load_counts(entries: dict[str, dict]) -> dict[str, dict]:
    """Residue counts of the master molecule, plus the free-fraction vector for the A1c check."""
    counts, missing = {}, []
    for r in common.read_tsv(INPUTS["composition"]):
        if r["segment_set"] != "master":
            continue
        acc = r["accession"]
        if acc not in entries:
            continue
        counts[acc] = {"count": {a: int(r[f"count_{a}"]) for a in AA},
                       "free_frac": {a: float(r[f"free_frac_{a}"]) for a in AA},
                       "residue_frac": {a: float(r[f"residue_frac_{a}"]) for a in AA},
                       "mw": float(r["mw"])}
    for acc in entries:
        if acc not in counts:
            missing.append(acc)
    if missing:
        raise SystemExit(f"[STOP] {len(missing)} weighted entries have no master composition row, e.g. {missing[:5]}")
    return counts


# ----------------------------------------------------------------------------
# the profile (A1)
# ----------------------------------------------------------------------------

def weight_vector(entries, profile: str, ft: str) -> dict[str, float]:
    """The weights of a profile's set for one fiber type (A2). Tier-1 weights sum to one over
    the tier; combined over everything."""
    if profile == "combined":
        return {acc: e[f"w_{ft}_combined"] for acc, e in entries.items()}
    if profile == "contractile":
        return {acc: e[f"w_{ft}_tier1"] for acc, e in entries.items() if e[f"w_{ft}_tier1"] is not None}
    raise ValueError(profile)


def molar_profile(weights: dict[str, float], counts: dict[str, dict], masses: dict[str, dict], convention: str) -> dict[str, float]:
    """A1: n_a = sum_i (w_i / MW_i) * count_{i,a}; p_a = n_a * m_a / sum. Returns fractions summing to one."""
    n = {a: 0.0 for a in AA}
    for acc, w in weights.items():
        if w == 0.0:
            continue
        c = counts[acc]
        per_mol = w / c["mw"]
        for a in AA:
            n[a] += per_mol * c["count"][a]
    mass = {a: n[a] * masses[a][convention] for a in AA}
    tot = sum(mass.values())
    return {a: (mass[a] / tot if tot > 0 else 0.0) for a in AA}


def molar_amounts(weights: dict[str, float], counts: dict[str, dict]) -> dict[str, float]:
    """mol of residue a per gram of protein (the w_i are g/g, MW g/mol)."""
    n = {a: 0.0 for a in AA}
    for acc, w in weights.items():
        if w == 0.0:
            continue
        c = counts[acc]
        per_mol = w / c["mw"]
        for a in AA:
            n[a] += per_mol * c["count"][a]
    return n


def fraction_mixture(weights: dict[str, float], counts: dict[str, dict], convention: str) -> dict[str, float]:
    """A1c: sum_i w_i * f_{i,a}, renormalised — the check against A1."""
    key = f"{convention}_frac"
    p = {a: 0.0 for a in AA}
    for acc, w in weights.items():
        if w == 0.0:
            continue
        for a in AA:
            p[a] += w * counts[acc][key][a]
    tot = sum(p.values())
    return {a: (p[a] / tot if tot > 0 else 0.0) for a in AA}


def g_per_100g_protein(weights, counts, masses, convention) -> dict[str, float]:
    """Grams of amino acid (residue or free) per 100 g protein. Residue sums to 100; free sums to
    more, by the water added on hydrolysis."""
    n = molar_amounts(weights, counts)
    return {a: n[a] * masses[a][convention] * 100.0 for a in AA}


# ----------------------------------------------------------------------------
# MHC:actin sensitivity (A6)
# ----------------------------------------------------------------------------

def load_band_families() -> dict[str, set[str]]:
    fams: dict[str, set[str]] = defaultdict(set)
    for r in read_tsv_skip_comments(INPUTS["band_families"]):
        if r["cutoff"] == BAND_CUTOFF:
            fams[r["band"]].add(r["accession"])
    if not fams:
        raise SystemExit(f"[STOP] no band families at cutoff {BAND_CUTOFF} in {INPUTS['band_families'].name}")
    return dict(fams)


def load_measured_ratios() -> dict[str, list[tuple[str, float]]]:
    """Per gel fiber type: [(method label, MHC:actin ratio)] for every method in the table at the
    D52 cutoff. The gel's own ratio is one of them. Labels are the table's, not typed here."""
    out: dict[str, list[tuple[str, float]]] = {ft: [] for ft in GEL_TYPES}
    seen_gel = set()
    for r in read_tsv_skip_comments(INPUTS["mhc_actin"]):
        if r["cutoff"] != BAND_CUTOFF:
            continue
        ftcell = r["fiber_type"]
        # rows for other sources are labelled "<pool> (vs gel type <T>)"; map them to that gel type
        gel_ft = None
        for ft in GEL_TYPES:
            if ftcell == ft or ftcell.endswith(f"(vs gel type {ft})"):
                gel_ft = ft
        if gel_ft is None:
            continue
        label = r["method"] if ftcell == gel_ft else f"{r['method']} [{ftcell}]"
        out[gel_ft].append((label, float(r["mhc_to_actin_method"])))
        if gel_ft not in seen_gel:
            out[gel_ft].append(("gel", float(r["mhc_to_actin_gel"])))
            seen_gel.add(gel_ft)
    for ft in GEL_TYPES:
        if len(out[ft]) < 2:
            raise SystemExit(f"[STOP] fewer than two MHC:actin measurements for type {ft} in {INPUTS['mhc_actin'].name}")
    return out


def rescaled_weights(weights: dict[str, float], mhc: set[str], actin: set[str], target: float, way: str) -> tuple[dict[str, float], float]:
    """Scale one band family so MHC:actin (sum over families) equals `target`, then renormalise.
    Returns the new weights and the ratio before scaling."""
    wm = sum(weights.get(a, 0.0) for a in mhc)
    wa = sum(weights.get(a, 0.0) for a in actin)
    if wa <= 0 or wm <= 0:
        raise SystemExit("[STOP] a band family has zero weight; the ratio is undefined (A6)")
    r0 = wm / wa
    factor = target / r0
    new = dict(weights)
    if way == "scale_mhc":
        for a in mhc:
            if a in new:
                new[a] *= factor
    elif way == "scale_actin":
        for a in actin:
            if a in new:
                new[a] /= factor
    else:
        raise ValueError(way)
    tot = sum(new.values())
    return {a: w / tot for a, w in new.items()}, r0


# ----------------------------------------------------------------------------
# bounds after weighting
# ----------------------------------------------------------------------------

def load_per_entry_bounds() -> dict[str, dict[str, dict[str, float]]]:
    """Per kind, per entry, per amino acid: the unweighted bound. isoform: largest |delta| over the
    entry's isoforms; processing: |delta| mature vs full; ptm: the entry's PTM mass fraction of MW,
    applied to every amino acid (a mass bound, not a per-amino-acid one)."""
    iso: dict[str, dict[str, float]] = defaultdict(lambda: {a: 0.0 for a in AA})
    if INPUTS["isoform_deltas"].exists():
        for r in common.read_tsv(INPUTS["isoform_deltas"]):
            for a in AA:
                iso[r["accession"]][a] = max(iso[r["accession"]][a], abs(float(r[f"delta_{a}"])))
    proc: dict[str, dict[str, float]] = {}
    for r in common.read_tsv(INPUTS["processing_deltas"]):
        proc[r["accession"]] = {a: abs(float(r[f"delta_free_frac_{a}"])) for a in AA}
    ptm: dict[str, dict[str, float]] = {}
    glyc: dict[str, dict[str, float]] = {}
    for r in common.read_tsv(INPUTS["ptm_deltas"]):
        f = abs(float(r["ptm_mass_delta_fraction_of_mw"]))
        ptm[r["accession"]] = {a: f for a in AA}
        if int(float(r.get("n_glycosylation", "0") or 0)) > 0:
            glyc[r["accession"]] = {a: 1.0 for a in AA}
    return {"isoform": dict(iso), "processing": proc, "ptm_mass": ptm, "glycosylated_entries_sum_w": glyc}


def load_family_bounds() -> dict[str, list[dict]]:
    """Within-tier family bounds (D55) per (tier, fiber type)."""
    out: dict[str, list[dict]] = defaultdict(list)
    for r in read_tsv_skip_comments(INPUTS["family_bounds"]):
        out[f"{r['tier']}/{r['fiber_type']}"].append(r)
    return dict(out)


# ----------------------------------------------------------------------------
# EAA (D62)
# ----------------------------------------------------------------------------

def load_eaa() -> dict:
    """The indispensable amino acids and their groupings as the author transcribed them, resolved
    to one-letter symbols through the IUPAC table. Any `___` is a stop."""
    cp = common.read_ini(INPUTS["eaa"])
    sym = common.read_ini(INPUTS["symbols"])
    three_to_one = {sym[a]["three_letter"].strip().lower(): a for a in sym.sections() if a in AA}
    if len(three_to_one) != 20:
        raise SystemExit(f"[STOP] {INPUTS['symbols'].name} does not give twenty three-letter symbols")

    def need(section, key):
        if section not in cp or key not in cp[section]:
            raise SystemExit(f"[STOP] {INPUTS['eaa'].name}: [{section}] {key} missing")
        val = cp[section][key].strip()
        if val == "___" or not val:
            raise SystemExit(f"[STOP] {INPUTS['eaa'].name}: [{section}] {key} is not filled in (D62; the author transcribes it from the pinned PDF)")
        return val

    headings = [h.strip() for h in need("indispensable_amino_acids", "headings").split(";") if h.strip()]
    groups = {}
    for sec in cp.sections():
        if sec.startswith("group."):
            gname = sec[len("group."):]
            if gname == "___":
                raise SystemExit(f"[STOP] {INPUTS['eaa'].name}: a [group.___] template section is still present — name it or delete it (D62)")
            members = [m.strip() for m in need(sec, "members").split(";") if m.strip()]
            groups[gname] = {"members": members, "location": cp[sec].get("location", ""), "as_reported": cp[sec].get("as_reported", "")}
    resolved = []   # (heading, [one-letter symbols])
    for h in headings:
        if h in groups:
            letters = []
            for m in groups[h]["members"]:
                if m.lower() not in three_to_one:
                    raise SystemExit(f"[STOP] group {h}: member {m!r} is not a three-letter symbol in the IUPAC table")
                letters.append(three_to_one[m.lower()])
            resolved.append((h, letters, "group"))
        elif h.lower() in three_to_one:
            resolved.append((h, [three_to_one[h.lower()]], "single"))
        else:
            raise SystemExit(f"[STOP] heading {h!r} is neither a three-letter symbol in the IUPAC table nor a [group.{h}] section")
    return {"headings": resolved, "location": need("indispensable_amino_acids", "location"),
            "as_reported": cp["indispensable_amino_acids"].get("as_reported", ""),
            "meta": dict(cp["meta"]) if "meta" in cp else {}, "groups": groups}


# ----------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------

def tsv_text(header_lines: list[str], columns: list[str], rows: list[list], tool: str = "aggregate.py") -> str:
    out = [f"# GENERATED by {tool}. DO NOT EDIT BY HAND."] + [f"# {h}" for h in header_lines]
    out.append("\t".join(columns))
    for r in rows:
        out.append("\t".join("" if v is None else str(v) for v in r))
    return "\n".join(out) + "\n"


def build(log) -> dict[str, str]:
    masses = load_masses()
    entries, weights_header = load_entries()
    counts = load_counts(entries)
    hashes = {k: f"{p.relative_to(common.REPO_ROOT).as_posix()} sha256 {sha256_path(p)}" for k, p in INPUTS.items() if p.exists()}
    header_common = [f"generated = {common.iso_now()}"] + [f"input.{k} = {v}" for k, v in hashes.items()] + [
        "A1: profile = molar mixture n_a = sum_i (w_i / MW_i) count_{i,a}; p_a = n_a m_a / sum (m = residue mass or free mass). "
        "A1c: fraction mixture sum_i w_i f_{i,a}, renormalised, written as a check.",
        "A2: combined = every weighted entry (D58); contractile = tier 1 (within-tier weights, D31).",
    ]
    files: dict[str, str] = {}

    # ---- 1. profiles
    prof: dict[tuple, dict[str, float]] = {}
    g100: dict[tuple, dict[str, float]] = {}
    check: dict[tuple, dict[str, float]] = {}
    nweighted: dict[tuple, int] = {}
    rows = []
    for profile in PROFILES:
        for ft in FIBER_TYPES:
            w = weight_vector(entries, profile, ft)
            nweighted[(profile, ft)] = sum(1 for x in w.values() if x > 0)
            wsum = sum(w.values())
            for conv in CONVENTIONS:
                p = molar_profile(w, counts, masses, conv)
                c = fraction_mixture(w, counts, conv)
                prof[(profile, ft, conv)] = p
                check[(profile, ft, conv)] = c
                g100[(profile, ft, conv)] = g_per_100g_protein(w, counts, masses, conv)
                maxdiff = max(abs(p[a] - c[a]) for a in AA)
                rows.append([profile, ft, conv, nweighted[(profile, ft)], fnum(wsum), fnum(sum(p.values())),
                             fnum(maxdiff), max(AA, key=lambda a: abs(p[a] - c[a]))]
                            + [fnum(p[a]) for a in AA] + [fnum(c[a]) for a in AA])
                log.info("profile %s type %s %s: %d entries with weight; A1 vs A1c max |diff| %.3g", profile, ft, conv, nweighted[(profile, ft)], maxdiff)
    files["amino_acid_profiles.tsv"] = tsv_text(
        header_common + ["fractions of total amino acid mass in the named convention, summing to one (A1); a1c_<a> = the fraction-mixture check"],
        ["profile", "fiber_type", "convention", "entries_with_weight", "sum_of_weights", "sum_of_fractions", "max_abs_a1_minus_a1c", "amino_acid_at_max"]
        + [f"frac_{a}" for a in AA] + [f"a1c_{a}" for a in AA], rows)

    # wide tables: g per 100 g protein
    for conv in CONVENTIONS:
        cols = ["amino_acid", "three_letter", "name"] + [f"{p}_{ft}" for p in PROFILES for ft in FIBER_TYPES]
        body = []
        for a in AA:
            body.append([a, masses[a]["three"], masses[a]["name"]] + [f"{g100[(p, ft, conv)][a]:.4f}" for p in PROFILES for ft in FIBER_TYPES])
        body.append(["sum", "", ""] + [f"{sum(g100[(p, ft, conv)].values()):.4f}" for p in PROFILES for ft in FIBER_TYPES])
        note = ("residue convention: in-chain residue mass; the column sums to 100 g per 100 g protein less each chain's one terminal water (the MW includes it)"
                if conv == "residue" else
                "free convention: free amino acid mass (the USDA / laboratory convention, D7); the column sums to more than 100 g per 100 g protein by the water added on hydrolysis")
        files[f"amino_acid_profiles_{conv}_wide.tsv"] = tsv_text(header_common + [f"g of amino acid per 100 g protein; {note}"], cols, body)

    # ---- 2. differences (A5)
    rows = []
    for conv in CONVENTIONS:
        for ft in FIBER_TYPES:
            for a in AA:
                rows.append([conv, "combined_minus_contractile", ft, a, fnum(prof[("combined", ft, conv)][a] - prof[("contractile", ft, conv)][a])])
        for profile in PROFILES:
            for x, y in (("I", "IIa"), ("I", "IIx"), ("IIa", "IIx")):
                for a in AA:
                    rows.append([conv, f"{x}_minus_{y}", profile, a, fnum(prof[(profile, x, conv)][a] - prof[(profile, y, conv)][a])])
            for a in AA:
                vals = [prof[(profile, ft, conv)][a] for ft in FIBER_TYPES]
                rows.append([conv, "fiber_type_bracket_max_minus_min", profile, a, fnum(max(vals) - min(vals))])
    files["profile_differences.tsv"] = tsv_text(
        header_common + ["difference of fractions (convention as named); fiber_type_bracket = max - min over the three pure types, which bounds every mix of them (A5)"],
        ["convention", "difference", "set", "amino_acid", "value"], rows)
    bracket = {}
    for conv in CONVENTIONS:
        for profile in PROFILES:
            bracket[(profile, conv)] = max(max(prof[(profile, ft, conv)][a] for ft in FIBER_TYPES) - min(prof[(profile, ft, conv)][a] for ft in FIBER_TYPES) for a in AA)

    # ---- 3. MHC:actin sensitivity (A6)
    bands = load_band_families()
    if "myosin_heavy_chain" not in bands or "actin" not in bands:
        raise SystemExit(f"[STOP] band families must name myosin_heavy_chain and actin (config/carroll_classical_fractionation.ini): {sorted(bands)}")
    ratios = load_measured_ratios()
    sens_rows, spread_rows = [], []
    mhc_shift: dict[tuple, dict[str, float]] = {}
    identity: dict[tuple, float] = {}
    for profile in PROFILES:
        for ft in GEL_TYPES:
            base_w = weight_vector(entries, profile, ft)
            variants: list[tuple[str, float, str, dict]] = []
            for label, target in ratios[ft]:
                for way in ("scale_mhc", "scale_actin"):
                    w2, r0 = rescaled_weights(base_w, bands["myosin_heavy_chain"], bands["actin"], target, way)
                    variants.append((label, target, way, w2))
            for conv in CONVENTIONS:
                base = prof[(profile, ft, conv)]
                per_variant = []
                for label, target, way, w2 in variants:
                    p2 = molar_profile(w2, counts, masses, conv)
                    per_variant.append(p2)
                    sens_rows.append([profile, ft, conv, label, fnum(target), fnum(r0), way] + [fnum(p2[a]) for a in AA])
                # identity check: the variants at the ratio nearest this repository's own reproduce the standard
                nearest = min((abs(target - r0) for (label, target, way, w2) in variants))
                own = [p2 for (label, target, way, w2), p2 in zip(variants, per_variant) if abs(target - r0) == nearest]
                identity[(profile, ft, conv)] = max((abs(p2[a] - base[a]) for p2 in own for a in AA), default=float("nan"))
                shifts = {}
                for a in AA:
                    vals = [p2[a] for p2 in per_variant]
                    shifts[a] = max(abs(v - base[a]) for v in vals)
                    spread_rows.append([profile, ft, conv, a, fnum(base[a]), fnum(min(vals)), fnum(max(vals)), fnum(max(vals) - min(vals)), fnum(shifts[a])])
                mhc_shift[(profile, ft, conv)] = shifts
    files["sensitivity_mhc_actin_profiles.tsv"] = tsv_text(
        header_common + ["A6: the profile when the MHC band family (D52, cutoff 2) or the actin band family is rescaled so MHC:actin equals each ratio the methods measured; the ratio labelled ibaq_x_mw_within_tier1 is this repository's own and reproduces the standard;",
                         "no method is a reference for another; the other sources' slow/fast pools are placed against the gel's fiber types as classical_check_carroll_2004.tsv places them; types I and IIa only — no factor is invented for IIx"],
        ["profile", "fiber_type", "convention", "ratio_source", "mhc_to_actin_target", "mhc_to_actin_before_scaling", "way"] + [f"frac_{a}" for a in AA], sens_rows)
    files["sensitivity_mhc_actin_spread.tsv"] = tsv_text(
        header_common + ["per amino acid: the profile's range over every measured ratio and both single-band rescalings; max_abs_shift = the largest |profile - standard| (A6)"],
        ["profile", "fiber_type", "convention", "amino_acid", "standard", "min_over_variants", "max_over_variants", "spread", "max_abs_shift"], spread_rows)

    # ---- 4. bounds after weighting
    per_entry = load_per_entry_bounds()
    fam = load_family_bounds()
    bound_rows = []
    bound_table: dict[tuple, dict[str, dict[str, float]]] = {}   # (profile, ft) -> kind -> aa -> worst-case sum
    for profile in PROFILES:
        for ft in FIBER_TYPES:
            w = weight_vector(entries, profile, ft)
            bound_table[(profile, ft)] = {}
            for kind, table in per_entry.items():
                sums, best, where = {a: 0.0 for a in AA}, {a: 0.0 for a in AA}, {a: "" for a in AA}
                for acc, ww in w.items():
                    if ww == 0.0:
                        continue
                    b = table.get(acc)
                    if not b:
                        continue
                    for a in AA:
                        d = b[a] * ww
                        sums[a] += d
                        if d > best[a]:
                            best[a], where[a] = d, acc
                bound_table[(profile, ft)][kind] = sums
                for a in AA:
                    bound_rows.append([profile, ft, kind, a, fnum(best[a]), where[a], fnum(sums[a])])
            # family bounds: within-tier rows scaled to the profile's denominator
            fam_sum, fam_best, fam_where = {a: 0.0 for a in AA}, {a: 0.0 for a in AA}, {a: "" for a in AA}
            if profile == "contractile":
                scale = {"1": 1.0}
            else:
                scale = {}
                for t in ("1", "2"):
                    scale[t] = sum(e[f"w_{ft}_combined"] for e in entries.values() if t in e["tiers"])
            for t, s in scale.items():
                for r in fam.get(f"{t}/{ft}", []):
                    for a in AA:
                        d = float(r[f"bound_{a}"]) * s
                        fam_sum[a] += d
                        if d > fam_best[a]:
                            fam_best[a], fam_where[a] = d, r["genes"]
            bound_table[(profile, ft)]["family"] = fam_sum
            for a in AA:
                bound_rows.append([profile, ft, "family", a, fnum(fam_best[a]), fam_where[a], fnum(fam_sum[a])])
            # shared-row excess (D48): the profile if every shared row belonged wholly to each member
            for conv in CONVENTIONS:
                base = prof[(profile, ft, conv)]
                w_plus = {}
                if profile == "combined":
                    d = next(iter(entries.values()))[f"D_{ft}"]
                    for acc, e in entries.items():
                        w_plus[acc] = (e[f"v_{ft}"] + e[f"shared_excess_{ft}"]) * e["mw"] / d
                else:
                    members = {acc: e for acc, e in entries.items() if "1" in e["tiers"]}
                    d1 = sum(e[f"v_{ft}"] * e["mw"] for e in members.values())
                    for acc, e in members.items():
                        w_plus[acc] = (e[f"v_{ft}"] + e[f"shared_excess_{ft}"]) * e["mw"] / d1
                tot = sum(w_plus.values())
                w_plus = {k: v / tot for k, v in w_plus.items()}
                p_plus = molar_profile(w_plus, counts, masses, conv)
                shifts = {a: abs(p_plus[a] - base[a]) for a in AA}
                bound_table[(profile, ft)][f"shared_row_{conv}"] = shifts
                for a in AA:
                    bound_rows.append([profile, ft, f"shared_row_excess_{conv}", a, fnum(shifts[a]), "", fnum(shifts[a])])
    files["bounds_after_weighting.tsv"] = tsv_text(
        header_common + ["per-entry bound (isoform: largest |delta| over the entry's isoforms; processing: |delta| mature vs full; ptm_mass: PTM mass / MW applied to every amino acid; glycosylated_entries_sum_w: 1 per glycosylated entry) x the entry's weight;",
                         "family: the D55 within-tier bounds scaled by the tier's share of the profile; shared_row_excess: |profile with every shared row given wholly to each member - standard| (D48).",
                         "largest_single = the biggest single entry or family; worst_case_sum = summed over every entry or family (all errors in one direction). Both are bounds on a fraction, not estimates."],
        ["profile", "fiber_type", "bound", "amino_acid", "largest_single", "where", "worst_case_sum"], bound_rows)

    # ---- 5. completeness (A7)
    outside = read_tsv_skip_comments(INPUTS["outside_pool"])
    comp_rows = []
    comp_bound: dict[str, float] = {}
    for ft in FIBER_TYPES:
        d = next(iter(entries.values()))[f"D_{ft}"]
        sum_v = sum(e[f"v_{ft}"] for e in entries.values())
        mw_molar_mean = d / sum_v if sum_v > 0 else float("nan")
        ranked = sorted(outside, key=lambda r: -float(r[f"v_{ft}"]))
        ks = [k for k in COMPLETENESS_RANKS if k <= len(ranked)] + [len(ranked)]
        for k in sorted(set(ks)):
            v_sum = sum(float(r[f"v_{ft}"]) for r in ranked[:k])
            molar_share = v_sum / (sum_v + v_sum) if sum_v + v_sum > 0 else 0.0
            vm = v_sum * mw_molar_mean
            mass_share = vm / (d + vm) if d + vm > 0 else 0.0
            comp_rows.append([ft, k, "all" if k == len(ranked) else "top", fnum(v_sum), fnum(molar_share), f"{mw_molar_mean:.1f}", fnum(mass_share), fnum(mass_share)])
            if k == len(ranked):
                comp_bound[ft] = mass_share
    files["completeness_sensitivity.tsv"] = tsv_text(
        header_common + ["A7: the genes the primary dataset quantified that map to no pool entry, ranked by their median value per fiber type; molar_share = their share of all quantified signal if included;",
                         "mass_share assumes each has the pool's molar-mean molecular weight (sum v*MW / sum v) — an assumption, stated, under which the mass share equals the molar share; bound_on_any_fraction = mass_share: a share of unknown composition can move no amino acid's fraction by more than itself"],
        ["fiber_type", "n_top_genes", "scope", "sum_v", "molar_share_if_included", "assumed_mw_g_per_mol", "mass_share_if_included", "bound_on_any_fraction"], comp_rows)

    # ---- 6. EAA (D62) — last, so everything above is written even if the transcription is not done
    try:
        eaa = load_eaa()
    except SystemExit as e:
        files["__stop__"] = str(e.code)
        eaa = None
    if eaa is None:
        files["standard_summary.ini"] = summary_text(hashes, header_common, nweighted, prof, check, g100, bound_table, mhc_shift, bracket, comp_bound, None, identity)
        return files
    eaa_rows = []
    for profile in PROFILES:
        for ft in FIBER_TYPES:
            w = weight_vector(entries, profile, ft)
            n = molar_amounts(w, counts)
            tot_free = sum(n[a] * masses[a]["free"] for a in AA)
            for heading, letters, kind in eaa["headings"]:
                mg_per_g_protein = sum(n[a] * masses[a]["free"] for a in letters) * 1000.0
                frac_free = sum(n[a] * masses[a]["free"] for a in letters) / tot_free
                eaa_rows.append([profile, ft, heading, kind, ";".join(letters), f"{mg_per_g_protein:.2f}", fnum(frac_free)])
            all_letters = sorted({l for _, ls, _ in eaa["headings"] for l in ls})
            eaa_rows.append([profile, ft, "sum_of_headings", "sum", ";".join(all_letters),
                             f"{sum(n[a] * masses[a]['free'] for a in all_letters) * 1000.0:.2f}",
                             fnum(sum(n[a] * masses[a]['free'] for a in all_letters) / tot_free)])
    files["eaa_subset.tsv"] = tsv_text(
        header_common + [f"D62: headings and groupings as transcribed by the author from the FAO 2013 report ({eaa['location']}); as_reported: {eaa['as_reported']}",
                         "mg_free_amino_acid_per_g_protein: free amino acid mass per gram of protein (residue-mass basis) — the unit of the report's scoring patterns; fraction_of_free_mass: share of the free-convention profile"],
        ["profile", "fiber_type", "heading", "kind", "letters", "mg_free_amino_acid_per_g_protein", "fraction_of_free_mass"], eaa_rows)

    files["standard_summary.ini"] = summary_text(hashes, header_common, nweighted, prof, check, g100, bound_table, mhc_shift, bracket, comp_bound, eaa, identity)
    return files


def summary_text(hashes, header_common, nweighted, prof, check, g100, bound_table, mhc_shift, bracket, comp_bound, eaa, identity=None) -> str:
    identity = identity or {}
    summ = common.new_ini()
    summ.add_section("inputs")
    for k, v in hashes.items():
        summ["inputs"][k] = v
    for profile in PROFILES:
        for ft in FIBER_TYPES:
            sec = f"{profile}.type_{ft}"
            summ.add_section(sec)
            summ[sec]["entries_with_weight"] = str(nweighted[(profile, ft)])
            for conv in CONVENTIONS:
                p = prof[(profile, ft, conv)]
                top = sorted(AA, key=lambda a: -p[a])[:3]
                summ[sec][f"{conv}_three_largest"] = "; ".join(f"{a} {fnum(p[a])}" for a in top)
                summ[sec][f"{conv}_max_abs_a1_minus_a1c"] = fnum(max(abs(p[a] - check[(profile, ft, conv)][a]) for a in AA))
                summ[sec][f"{conv}_g_per_100g_protein_sum"] = f"{sum(g100[(profile, ft, conv)].values()):.4f}"
            for kind, sums in bound_table[(profile, ft)].items():
                a_max = max(AA, key=lambda a: sums[a])
                if kind == "glycosylated_entries_sum_w":
                    summ[sec]["glycosylated_entries_sum_w"] = fnum(sums[a_max])      # a mass share, not a shift
                else:
                    summ[sec][f"bound_{kind}_worst_case_max"] = f"{fnum(sums[a_max])} ({a_max})"
            if ft in GEL_TYPES:
                for conv in CONVENTIONS:
                    s = mhc_shift[(profile, ft, conv)]
                    a_max = max(AA, key=lambda a: s[a])
                    summ[sec][f"mhc_actin_max_abs_shift_{conv}"] = f"{fnum(s[a_max])} ({a_max})"
                    summ[sec][f"mhc_actin_identity_check_{conv}"] = fnum(identity.get((profile, ft, conv), float("nan")))
    summ.add_section("differences")
    for conv in CONVENTIONS:
        for ft in FIBER_TYPES:
            diffs = {a: abs(prof[("combined", ft, conv)][a] - prof[("contractile", ft, conv)][a]) for a in AA}
            a_max = max(AA, key=lambda a: diffs[a])
            summ["differences"][f"{conv}_combined_minus_contractile_max_abs_{ft}"] = f"{fnum(diffs[a_max])} ({a_max})"
        for profile in PROFILES:
            summ["differences"][f"{conv}_fiber_type_bracket_max_{profile}"] = fnum(bracket[(profile, conv)])
    summ.add_section("completeness")
    for ft in FIBER_TYPES:
        summ["completeness"][f"bound_all_outside_genes_{ft}"] = fnum(comp_bound[ft])
    summ.add_section("eaa")
    if eaa is None:
        summ["eaa"]["status"] = "not computed: config/aggregation_decisions.ini is not filled in (D62)"
    else:
        summ["eaa"]["headings"] = "; ".join(f"{h}={''.join(ls)}" for h, ls, _ in eaa["headings"])
        summ["eaa"]["location"] = eaa["location"]
    return common.render_ini(summ, ["GENERATED by aggregate.py. DO NOT EDIT BY HAND."] + header_common)


def strip_generated_line(text: str) -> str:
    return "\n".join(l for l in text.splitlines() if not l.startswith("# generated"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("aggregate_check" if args.check else "aggregate")
    files = build(log)
    stop = files.pop("__stop__", None)
    if args.check:
        stale = [rel for rel, content in files.items()
                 if not (OUT_DIR / rel).exists() or strip_generated_line((OUT_DIR / rel).read_text(encoding="utf-8")) != strip_generated_line(content)]
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), stale)
            return 1
        log.info("check: standard outputs are current")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for rel, content in files.items():
        (OUT_DIR / rel).write_text(content, encoding="utf-8", newline="\n")
        log.info("wrote %s", (OUT_DIR / rel).relative_to(common.REPO_ROOT).as_posix())
    if stop:
        # Not a stop for the stage: every other output is written and the later stages can run.
        # The unfilled transcription is recorded in standard_summary.ini [eaa] and a test fails
        # until it is filled (the flag pattern: unsettled is visible, never defaulted).
        log.warning("EAA subset NOT written: %s", stop)
    return 0


if __name__ == "__main__":
    sys.exit(main())
