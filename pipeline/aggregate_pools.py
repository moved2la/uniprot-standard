"""The amino acid standard of a category whose pools are mixed by a MEASURED SPLIT (Step 7b, blood).

Muscle's aggregate mixes one dataset's fiber-type columns by a fiber-type mix. Blood is two pools
with their own weights (config/<category>/mass_fractions_per_entry.tsv), each a compartment, mixed
by the compartments' protein masses (D100, config/<category>/compartment_mix.ini). The weights are
mass shares (D105), so every profile here is a mass-weighted mean of per-protein mass fractions:

    profile_pool[a]  = 100 * sum_i w_i * frac_i[a]         (g of amino acid a per 100 g protein)
    standard[a]      = f_plasma * profile_plasma[a] + (1 - f_plasma) * profile_erythrocytes[a]
    f_plasma         = plasma protein / (plasma protein + erythrocyte protein)
    plasma protein   = plasma volume x plasma protein concentration
    erythrocyte protein = haemoglobin mass / haemoglobin's share of erythrocyte protein
    haemoglobin mass = red-cell volume x MCHC                  (D119; blood Hb x blood volume as a cross-check)

Both conventions (free amino acids; residues) are written, as for muscle.

The split and the dominant proteins (D100, D102, D117)
  Haemoglobin's share of erythrocyte protein is not a config number: the base standard uses the
  share the weights measure (the summed w of the haemoglobin entries, ~0.53), and the sensitivity
  recomputes everything at each published share (0.70, 0.95) -- the erythrocyte PROFILE (the
  haemoglobin entries scaled to s, every other entry to 1 - s: the two-piece compartment) AND the
  compartment's protein mass. Plasma has the same shape for albumin (measured ~0.24 by LFQ; 0.51 by
  clinical concentration), with the plasma protein mass held at its measured concentration.
  The standard is reported at every combination; the base is as measured.

Other sensitivities   the per-donor profiles of the replicate-share pool (D108); the profile under
                      each donor alone, under the mean (the weight) and under the median
Uncertainty (D63)     Monte Carlo in log space per pool from v, sd, n: a log-normal fitted to (v, sd)
                      exactly; `spread` draws each entry's value, `mean_uncertainty` draws with
                      sigma / sqrt(n). Each draw recomputes the whole profile. The standard's draws
                      combine the two pools at the base split.
Bounds                the weighted isoform / processing / PTM bounds, the D113 immunoglobulin bound,
                      the completeness gap, all on the one scale of uncertainty_per_amino_acid.tsv
Drivers               which entries carry each amino acid, and histidine in particular: the entries,
                      the compartments, and each sensitivity's effect on it, in one place

Reads   config/<category>/mass_fractions_per_entry.tsv, compartment_mix.ini, protein_set_decisions.ini,
        outputs/<category>/intermediate/composition/amino_acid_composition_per_protein.tsv,
        outputs/<category>/intermediate/mass_fractions/{weighted_bounds, immunoglobulin_variable_domain_bound_per_chain}.tsv,
        mass_fractions_summary.ini, config/fao_2013_indispensable_amino_acids.ini, data/iupac/amino_acid_symbols.ini,
        config/uncertainty_settings.ini
Writes  outputs/<category>/standard/*  (tables, standard_summary.ini, plots/)
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path

import numpy as np

from pipeline import common
from pipeline.mass_fractions_pools import fnum, tsv_text, strip_generated_line

STAGE = "aggregate_pools"
AA = list(common.AMINO_ACIDS)
CONVENTIONS = ("free", "residue")
FRAC = {"free": "free_frac", "residue": "residue_frac"}
CHUNK = 500
MC_TERMS = ("spread", "mean_uncertainty")


def read_commented_tsv(path: Path) -> list[dict[str, str]]:
    lines = [l for l in path.read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l]


def split_list(value: str) -> list[str]:
    return [x.strip() for x in value.replace(";", ",").split(",") if x.strip()]


def need(cp, section: str, key: str, path: Path) -> str:
    if section not in cp or key not in cp[section] or cp[section][key].strip() in ("", "___"):
        raise SystemExit(f"[STOP] {path.name}: [{section}] {key} is not filled in")
    return cp[section][key].strip()


# ----------------------------------------------------------------------------- inputs

def load_weights(files, pools: list[str]) -> dict[str, dict]:
    rows = read_commented_tsv(files["mass_fractions_tsv"])
    out = {}
    for r in rows:
        e = {"gene": r["gene"], "w": {}, "v": {}, "sd": {}, "n": {}, "reps": {}}
        for p in pools:
            if r.get(f"w_{p}", "") == "":
                continue
            e["w"][p] = float(r[f"w_{p}"]); e["v"][p] = float(r[f"v_{p}"]); e["sd"][p] = float(r[f"sd_{p}"]); e["n"][p] = int(r[f"n_{p}"])
            reps = []
            k = 1
            while f"v_{p}_rep{k}" in r:
                reps.append(float(r[f"v_{p}_rep{k}"]) if r[f"v_{p}_rep{k}"] else 0.0)
                k += 1
            if reps:
                e["reps"][p] = reps
        out[r["accession"]] = e
    return out


def load_composition(files) -> dict[str, dict]:
    out = {}
    for r in common.read_tsv(files["composition_tsv"]):
        if r["segment_set"] == "master":
            out[r["accession"]] = {conv: np.array([float(r[f"{FRAC[conv]}_{a}"]) for a in AA]) for conv in CONVENTIONS}
    return out


def load_eaa() -> dict:
    """The FAO 2013 indispensable amino acids as the author transcribed them (D62), resolved to
    one-letter symbols through the IUPAC table; paths read at call time."""
    path, sym_path = common.EAA_INI, common.AMINO_ACID_SYMBOLS_INI
    cp, sym = common.read_ini(path), common.read_ini(sym_path)
    three_to_one = {sym[a]["three_letter"].strip().lower(): a for a in sym.sections() if a in AA}
    headings = split_list(need(cp, "indispensable_amino_acids", "headings", path))
    groups = {s[len("group."):]: split_list(need(cp, s, "members", path)) for s in cp.sections() if s.startswith("group.")}
    resolved = []
    for h in headings:
        if h in groups:
            resolved.append((h, [three_to_one[m.lower()] for m in groups[h]]))
        elif h.lower() in three_to_one:
            resolved.append((h, [three_to_one[h.lower()]]))
        else:
            raise SystemExit(f"[STOP] {path.name}: heading {h!r} is neither a three-letter symbol nor a [group.{h}] section")
    return {"headings": resolved, "location": need(cp, "indispensable_amino_acids", "location", path)}


def load_mix(category: str) -> dict:
    """config/<category>/compartment_mix.ini: the D100 numbers with their locations, the dominant
    proteins of each compartment with their published shares (D102, D117), the cross-checks."""
    path = common.category_config_dir(category) / "compartment_mix.ini"
    cp = common.read_ini(path)
    mix = {"path": path, "sex": need(cp, "reference", "sex", path), "compartments": {}, "dominant": {}, "cross_checks": {}}
    for s in cp.sections():
        if s.startswith("compartment."):
            name = s[len("compartment."):]
            sec = cp[s]
            mix["compartments"][name] = {k: sec[k].strip() for k in sec}
            for k in ("volume_ml", "volume_location"):
                need(cp, s, k, path)
        elif s.startswith("dominant."):
            pool = s[len("dominant."):]
            sec = cp[s]
            shares = []
            for k in sec:
                if k.startswith("share.") and not k.endswith(".location"):
                    val = sec[k].strip()
                    shares.append((k[len("share."):], None if val == "computed" else float(val), sec.get(f"{k}.location", "")))
            mix["dominant"][pool] = {"genes": split_list(need(cp, s, "genes", path)), "genes_location": sec.get("genes_location", ""),
                                     "shares": shares}
        elif s.startswith("cross_check."):
            mix["cross_checks"][s[len("cross_check."):]] = {k: cp[s][k].strip() for k in cp[s]}
    mix["hb_shares"] = mix["dominant"]["erythrocytes"]["shares"] if "erythrocytes" in mix["dominant"] else []
    return mix


# ----------------------------------------------------------------------------- profiles and the split

def profile(weights: dict[str, float], comp, conv: str) -> np.ndarray:
    """g of each amino acid per 100 g protein: 100 x the w-weighted mean of the per-protein mass fractions."""
    tot = sum(weights.values())
    acc = np.zeros(len(AA))
    if tot <= 0:
        return acc
    for a, w in weights.items():
        if w > 0 and a in comp:
            acc += (w / tot) * comp[a][conv]
    return 100.0 * acc


def rescale_dominant(weights: dict[str, float], dominant: set[str], share: float) -> dict[str, float]:
    """The two-piece compartment: the dominant entries scaled to `share` of the pool, every other
    entry to 1 - share, each group keeping its internal proportions."""
    s_dom = sum(w for a, w in weights.items() if a in dominant)
    s_rest = sum(w for a, w in weights.items() if a not in dominant)
    out = {}
    for a, w in weights.items():
        if a in dominant:
            out[a] = w * share / s_dom if s_dom > 0 else 0.0
        else:
            out[a] = w * (1 - share) / s_rest if s_rest > 0 else 0.0
    return out


def split(mix: dict, hb_share: float, female: bool = False) -> dict:
    """D100 / D119: the compartments' protein masses and the plasma fraction at one haemoglobin share.
    `female` uses the *_female volumes and blood haemoglobin where the config gives them (reported, not the base)."""
    c = mix["compartments"]
    pl, er = c["plasma"], c["erythrocytes"]
    sfx = "_female" if female else ""
    plasma_l = float(pl.get(f"volume_ml{sfx}", pl["volume_ml"])) / 1000
    rbc_l = float(er.get(f"volume_ml{sfx}", er["volume_ml"])) / 1000
    plasma_protein = plasma_l * float(pl["protein_g_per_l"])
    hb_mass = rbc_l * float(er["mchc_g_per_l"])
    hb_alt = er.get(f"blood_hb_g_per_l{sfx}", er.get("blood_hb_g_per_l"))
    hb_mass_alt = (plasma_l + rbc_l) * float(hb_alt) if hb_alt else None
    ery_protein = hb_mass / hb_share
    f_plasma = plasma_protein / (plasma_protein + ery_protein)
    return {"plasma_protein_g": plasma_protein, "hb_mass_g": hb_mass, "hb_mass_g_blood_route": hb_mass_alt,
            "erythrocyte_protein_g": ery_protein, "whole_blood_protein_g": plasma_protein + ery_protein, "f_plasma": f_plasma}


# ----------------------------------------------------------------------------- Monte Carlo (D63)

def lognormal_sigma(m: float, s: float) -> float:
    if m <= 0 or s <= 0:
        return 0.0
    x = (1.0 + math.sqrt(1.0 + 4.0 * (s * s) / (m * m))) / 2.0
    return math.sqrt(math.log(x))


def simulate(pools: list[str], entries: dict[str, dict], comp, draws: int, seed: int, f_plasma: float, log) -> dict[tuple, np.ndarray]:
    """Per (pool or 'standard', convention, term): an array (draws, 20) of g/100 g profiles."""
    rng = np.random.default_rng(seed)
    results: dict[tuple, list[np.ndarray]] = {}
    per_pool_draws: dict[tuple, np.ndarray] = {}
    for p in pools:
        accs = sorted(a for a in entries if p in entries[a]["w"] and a in comp)
        F = {conv: np.array([comp[a][conv] for a in accs]) for conv in CONVENTIONS}
        m = np.array([entries[a]["v"][p] for a in accs])
        s = np.array([entries[a]["sd"][p] for a in accs])
        n = np.array([max(entries[a]["n"][p], 1) for a in accs], dtype=float)
        sigma = np.array([lognormal_sigma(mi, si) for mi, si in zip(m, s)])
        mu = np.where(m > 0, np.log(np.where(m > 0, m, 1.0)), 0.0)
        for term in MC_TERMS:
            sig = sigma if term == "spread" else sigma / np.sqrt(n)
            chunks = {conv: [] for conv in CONVENTIONS}
            done = 0
            while done < draws:
                k = min(CHUNK, draws - done)
                z = rng.standard_normal((k, len(accs)))
                V = np.where(m > 0, np.exp(mu + sig * z), 0.0)
                W = V / V.sum(axis=1, keepdims=True)
                for conv in CONVENTIONS:
                    chunks[conv].append(100.0 * (W @ F[conv]))
                done += k
            for conv in CONVENTIONS:
                arr = np.concatenate(chunks[conv], axis=0)
                results[(p, conv, term)] = arr
                per_pool_draws[(p, conv, term)] = arr
            log.info("pool %s %s: %d draws", p, term, draws)
    if {"plasma", "erythrocytes"} <= set(pools):
        for conv in CONVENTIONS:
            for term in MC_TERMS:
                results[("standard", conv, term)] = f_plasma * per_pool_draws[("plasma", conv, term)] + (1 - f_plasma) * per_pool_draws[("erythrocytes", conv, term)]
    return results


# ----------------------------------------------------------------------------- build

def build(category: str, log, make_plots: bool = True) -> tuple[dict[Path, str], dict[Path, bytes]]:
    files = common.category_files(category)
    out = common.category_outputs(category)["standard"]
    pools = sorted(common.pool_names(category))
    entries = load_weights(files, pools)
    comp = load_composition(files)
    mix = load_mix(category)
    eaa = load_eaa()
    mf_dir = files["mass_fractions_out"]
    mf_summary = common.read_ini(mf_dir / "mass_fractions_summary.ini")
    sha = {"weights": common.sha256_file(files["mass_fractions_tsv"]), "composition": common.sha256_file(files["composition_tsv"]),
           "mix": common.sha256_file(mix["path"])}
    header = ["GENERATED by aggregate_pools.py. DO NOT EDIT BY HAND.", f"generated   = {common.iso_now()}", f"category    = {category}",
              f"weights     = {files['mass_fractions_tsv'].relative_to(common.REPO_ROOT).as_posix()} sha256 {sha['weights']}",
              f"composition = {files['composition_tsv'].relative_to(common.REPO_ROOT).as_posix()} sha256 {sha['composition']}",
              f"mix         = {mix['path'].relative_to(common.REPO_ROOT).as_posix()} sha256 {sha['mix']} (reference {mix['sex']})",
              "profile[a] = 100 x sum_i w_i x mass fraction of a in entry i (D105 weights are mass shares); standard = f_plasma x plasma + (1 - f_plasma) x erythrocytes (D100)"]

    W = {p: {a: e["w"][p] for a, e in entries.items() if p in e["w"]} for p in pools}
    P = {(p, conv): profile(W[p], comp, conv) for p in pools for conv in CONVENTIONS}

    # the dominant proteins as measured (D102 / D117): the summed weight of the named genes
    dom_sets = {p: {a for a, e in entries.items() if p in e["w"] and e["gene"] in set(mix["dominant"][p]["genes"])} for p in mix["dominant"]}
    measured_share = {p: sum(W[p][a] for a in dom_sets[p]) for p in dom_sets}
    for p in dom_sets:
        log.info("pool %s: dominant entries %s carry %.4f of the pool as measured", p, ", ".join(sorted(entries[a]["gene"] for a in dom_sets[p])), measured_share[p])
    hb_base = measured_share.get("erythrocytes", 1.0)
    base_split = split(mix, hb_base)
    cc = mix["cross_checks"]
    if "whole_blood_protein" in cc:                                  # the share ICRP's whole-blood protein implies (a cited-number computation)
        ref = float(cc["whole_blood_protein"]["blood_mass_g"]) * float(cc["whole_blood_protein"]["protein_fraction"])
        implied = base_split["hb_mass_g"] / (ref - base_split["plasma_protein_g"])
        mix["hb_shares"] = [(l, (implied if v is None else v), loc + (f" -> {implied:.4f}" if v is None else "")) for l, v, loc in mix["hb_shares"]]
        log.info("haemoglobin share implied by ICRP whole-blood protein %.0f g: %.4f", ref, implied)
    else:
        mix["hb_shares"] = [(l, v, loc) for l, v, loc in mix["hb_shares"] if v is not None]
    mix["dominant"]["erythrocytes"]["shares"] = mix["hb_shares"]
    f_pl = base_split["f_plasma"]
    S = {conv: f_pl * P[("plasma", conv)] + (1 - f_pl) * P[("erythrocytes", conv)] for conv in CONVENTIONS}
    log.info("base split at haemoglobin share %.4f (as measured): plasma protein %.1f g, erythrocyte protein %.1f g, f_plasma %.4f",
             hb_base, base_split["plasma_protein_g"], base_split["erythrocyte_protein_g"], f_pl)

    outputs: dict[Path, str] = {}
    images: dict[Path, bytes] = {}

    # -- the standard and the profiles
    for conv in CONVENTIONS:
        rows = [[a] + [fnum(P[(p, conv)][i]) for p in pools] + [fnum(S[conv][i])] for i, a in enumerate(AA)]
        name = "_calculated_amino_acid_standard.tsv" if conv == "free" else "_calculated_amino_acid_standard_residue_convention.tsv"
        outputs[out / name] = tsv_text(header + [f"g per 100 g protein, {conv} amino acid convention; columns per pool (D97: total_<pool>) and the standard at the base split"],
                                       ["amino_acid"] + [f"total_{p}" for p in pools] + ["standard"], rows)
    prof_rows = [[p, conv, a, fnum(P[(p, conv)][i])] for p in pools for conv in CONVENTIONS for i, a in enumerate(AA)]
    prof_rows += [["standard", conv, a, fnum(S[conv][i])] for conv in CONVENTIONS for i, a in enumerate(AA)]
    outputs[out / "amino_acid_profiles.tsv"] = tsv_text(header, ["profile", "convention", "amino_acid", "g_per_100g_protein"], prof_rows)

    # -- EAA subset
    eaa_rows = []
    for heading, letters in eaa["headings"]:
        eaa_rows.append([heading, "+".join(letters)] + [fnum(sum(P[(p, "free")][AA.index(x)] for x in letters)) for p in pools]
                        + [fnum(sum(S["free"][AA.index(x)] for x in letters))])
    outputs[out / "eaa_subset.tsv"] = tsv_text(header + [f"FAO 2013 indispensable amino acids as transcribed (D62; {eaa['location']}); free convention"],
                                               ["heading", "amino_acids"] + [f"total_{p}" for p in pools] + ["standard"], eaa_rows)

    # -- the split table (D100, D119): every published haemoglobin share, the base marked
    split_rows = []
    shares = [("as_measured", hb_base, "the summed weight of the haemoglobin entries in the erythrocyte pool")] + list(mix["hb_shares"])
    for label, s_val, loc in shares:
        sp = split(mix, s_val)
        split_rows.append([label, fnum(s_val), fnum(sp["plasma_protein_g"]), fnum(sp["hb_mass_g"]), fnum(sp["erythrocyte_protein_g"]),
                           fnum(sp["whole_blood_protein_g"]), fnum(sp["f_plasma"]), fnum(1 - sp["f_plasma"]), loc])
    if any(k.endswith("_female") for k in mix["compartments"]["plasma"]) and any(k.endswith("_female") for k in mix["compartments"]["erythrocytes"]):
        for label, s_val, loc in shares:
            sp = split(mix, s_val, female=True)
            split_rows.append([f"female {label}", fnum(s_val), fnum(sp["plasma_protein_g"]), fnum(sp["hb_mass_g"]), fnum(sp["erythrocyte_protein_g"]),
                               fnum(sp["whole_blood_protein_g"]), fnum(sp["f_plasma"]), fnum(1 - sp["f_plasma"]), "female volumes and blood Hb from the config; reported, not the base"])
    cross = [f"haemoglobin mass by the blood route (blood Hb x blood volume) = {fnum(base_split['hb_mass_g_blood_route'])} g against {fnum(base_split['hb_mass_g'])} g by red-cell volume x MCHC (D119)"]
    if "whole_blood_protein" in cc:
        ref = float(cc["whole_blood_protein"]["blood_mass_g"]) * float(cc["whole_blood_protein"]["protein_fraction"])
        cross.append(f"ICRP whole-blood protein = {fnum(ref)} g ({cc['whole_blood_protein'].get('location', '')}); the split's whole-blood protein per haemoglobin share is in the table")
    outputs[out / "compartment_split.tsv"] = tsv_text(header + cross,
                                                       ["haemoglobin_share", "share", "plasma_protein_g", "haemoglobin_g", "erythrocyte_protein_g", "whole_blood_protein_g", "f_plasma", "f_erythrocytes", "share_source"], split_rows)

    # -- D117: the dominant-protein sensitivity, every combination; base = as measured
    ery_shares = [("as_measured", hb_base)] + [(l, v) for l, v, _ in mix["hb_shares"]]
    pl_shares = [("as_measured", measured_share.get("plasma", 0.0))] + [(l, v) for l, v, _ in mix["dominant"].get("plasma", {"shares": []})["shares"]]
    sens_rows, sens_profiles = [], {}
    for el, es in ery_shares:
        sp = split(mix, es)
        W_e = rescale_dominant(W["erythrocytes"], dom_sets["erythrocytes"], es) if el != "as_measured" else W["erythrocytes"]
        for pl, ps in pl_shares:
            W_p = rescale_dominant(W["plasma"], dom_sets["plasma"], ps) if pl != "as_measured" else W["plasma"]
            for conv in CONVENTIONS:
                pe, pp = profile(W_e, comp, conv), profile(W_p, comp, conv)
                st = sp["f_plasma"] * pp + (1 - sp["f_plasma"]) * pe
                sens_profiles[(el, pl, conv)] = (pp, pe, st)
                for i, a in enumerate(AA):
                    sens_rows.append([el, fnum(es), pl, fnum(ps), fnum(sp["f_plasma"]), conv, a, fnum(pp[i]), fnum(pe[i]), fnum(st[i]),
                                      fnum(st[i] - S[conv][i])])
    outputs[out / "sensitivity_dominant_protein.tsv"] = tsv_text(header + [
        "D117: the erythrocyte profile and protein mass at each published haemoglobin share (the haemoglobin entries scaled to the share, every other entry to 1 - share); the plasma profile with albumin likewise, plasma mass fixed; the standard at every combination; base = as measured (D102)"],
        ["hb_share_label", "hb_share", "albumin_share_label", "albumin_share", "f_plasma", "convention", "amino_acid", "profile_plasma", "profile_erythrocytes", "standard", "standard_minus_base"], sens_rows)
    spread_rows = []
    for conv in CONVENTIONS:
        for i, a in enumerate(AA):
            vals = [sens_profiles[k][2][i] for k in sens_profiles if k[2] == conv]
            ery_only = [sens_profiles[(el, "as_measured", conv)][2][i] for el, _ in ery_shares]
            pl_only = [sens_profiles[("as_measured", pl, conv)][2][i] for pl, _ in pl_shares]
            spread_rows.append([conv, a, fnum(S[conv][i]), fnum(min(vals)), fnum(max(vals)), fnum(max(vals) - min(vals)),
                                fnum(max(ery_only) - min(ery_only)), fnum(max(pl_only) - min(pl_only))])
    outputs[out / "sensitivity_dominant_protein_spread.tsv"] = tsv_text(header + ["the standard's range over every D117 combination; haemoglobin_only / albumin_only vary one compartment with the other as measured"],
                                                                        ["convention", "amino_acid", "standard_base", "min", "max", "spread_all", "spread_haemoglobin_only", "spread_albumin_only"], spread_rows)

    # -- D108: per-donor profiles of the replicate-share pool, and the standard under each
    donor_rows, donor_spread = [], {}
    rep_pools = [p for p in pools if any(p in e["reps"] for e in entries.values())]
    for p in rep_pools:
        accs = [a for a in entries if p in entries[a]["reps"]]
        n_rep = len(entries[accs[0]]["reps"][p])
        variants = {f"donor_{k + 1}": {a: entries[a]["reps"][p][k] for a in accs} for k in range(n_rep)}
        variants["mean_of_donors (the weight)"] = {a: entries[a]["w"][p] for a in accs}
        variants["median_of_donors"] = {a: statistics.median(entries[a]["reps"][p]) for a in accs}
        for conv in CONVENTIONS:
            profs = {name: profile(wts, comp, conv) for name, wts in variants.items()}
            other = "plasma" if p == "erythrocytes" else "erythrocytes"
            for i, a in enumerate(AA):
                vals = [profs[name][i] for name in variants]
                for name in variants:
                    st = f_pl * P[(other, conv)][i] + (1 - f_pl) * profs[name][i] if p == "erythrocytes" else f_pl * profs[name][i] + (1 - f_pl) * P[(other, conv)][i]
                    donor_rows.append([p, name, conv, a, fnum(profs[name][i]), fnum(st)])
                donor_spread[(p, conv, a)] = (max(vals) - min(vals), max(vals[:n_rep]) - min(vals[:n_rep]))
    outputs[out / "sensitivity_per_donor.tsv"] = tsv_text(header + ["D108: the pool's profile under each donor's shares alone, under their mean (the weight) and under their median; the standard at the base split with that pool replaced"],
                                                          ["pool", "weights", "convention", "amino_acid", "profile_pool", "standard"], donor_rows)
    outputs[out / "sensitivity_per_donor_spread.tsv"] = tsv_text(header, ["pool", "convention", "amino_acid", "spread_all_variants", "spread_across_donors"],
                                                                 [[p, conv, a, fnum(v[0]), fnum(v[1])] for (p, conv, a), v in donor_spread.items()])

    # -- uncertainty (D63)
    cpu = common.read_ini(common.UNCERTAINTY_SETTINGS_INI)
    draws, seed = int(need(cpu, "monte_carlo", "draws", common.UNCERTAINTY_SETTINGS_INI)), int(need(cpu, "monte_carlo", "seed", common.UNCERTAINTY_SETTINGS_INI))
    mc = simulate(pools, entries, comp, draws, seed, f_pl, log)
    mc_rows, half = [], {}
    for (prof, conv, term), arr in sorted(mc.items()):
        base = S[conv] if prof == "standard" else P[(prof, conv)]
        lo, med, hi = np.percentile(arr, [2.5, 50, 97.5], axis=0)
        for i, a in enumerate(AA):
            hw = max(hi[i] - base[i], base[i] - lo[i])
            half[(prof, conv, term, a)] = hw
            mc_rows.append([prof, conv, term, a, fnum(base[i]), fnum(lo[i]), fnum(med[i]), fnum(hi[i]), fnum(hw)])
    outputs[out / "uncertainty_intervals.tsv"] = tsv_text(header + [f"D63: log-normal Monte Carlo, {draws} draws, seed {seed}; spread = each entry's value drawn; mean_uncertainty = sigma / sqrt(n); the standard combines the pools' draws at the base split"],
                                                          ["profile", "convention", "term", "amino_acid", "value", "p2.5", "p50", "p97.5", "half_width"], mc_rows)

    # -- bounds on one scale
    wb = read_commented_tsv(mf_dir / "weighted_bounds.tsv")
    ig = {(r["pool"], r["chain"]): float(r["variable_mass_bound"]) for r in read_commented_tsv(mf_dir / "immunoglobulin_variable_domain_bound_per_chain.tsv")}
    gap = {p: float(mf_summary[f"pool.{p}"]["completeness_gap_share_of_listing"]) for p in pools}
    terms_rows, widest = [], {}
    for prof in pools + ["standard"]:
        f = {"plasma": f_pl, "erythrocytes": 1 - f_pl} if prof == "standard" else {prof: 1.0}
        for conv in CONVENTIONS:
            for i, a in enumerate(AA):
                terms = {}
                for term in MC_TERMS:
                    terms[f"monte_carlo_{term}"] = half[(prof, conv, term, a)]
                for kind in ("isoform", "processing"):
                    terms[f"bound_{kind}"] = 100 * sum(f[p] * max((float(r["max_weighted_delta"]) for r in wb if r["pool"] == p and r["kind"] == kind and r["amino_acid"] == a), default=0.0) for p in f)
                terms["bound_ptm_mass"] = 100 * sum(f[p] * max((float(r["max_weighted_delta"]) for r in wb if r["pool"] == p and r["kind"] == "ptm_mass"), default=0.0) for p in f)
                terms["bound_ig_variable_domain"] = 100 * sum(f[p] * ig.get((p, "all"), 0.0) for p in f)
                terms["bound_completeness_gap"] = 100 * sum(f[p] * gap[p] for p in f)
                terms["sensitivity_dominant_protein"] = next(float(r[5]) for r in spread_rows if r[0] == conv and r[1] == a) if prof == "standard" else 0.0
                terms["sensitivity_per_donor"] = donor_spread.get((prof, conv, a), (0.0, 0.0))[1] if prof in rep_pools else 0.0
                base = S[conv][i] if prof == "standard" else P[(prof, conv)][i]
                top = max(terms, key=terms.get)
                widest[(prof, conv, a)] = top
                terms_rows.append([prof, conv, a, fnum(base)] + [fnum(terms[k]) for k in sorted(terms)] + [top])
    term_names = sorted(terms)
    outputs[out / "uncertainty_per_amino_acid.tsv"] = tsv_text(header + [
        "every term as a maximum absolute shift of the g/100 g value: Monte Carlo half-widths; bounds = the pool's weighted delta (x 100, mixed at the split); the completeness gap and the immunoglobulin bound as shares of the profile (A7: a share of unknown composition moves no value by more than itself); the two sensitivities as their spread"],
        ["profile", "convention", "amino_acid", "value"] + term_names + ["widest_term"], terms_rows)

    # -- drivers: which entries carry each amino acid; histidine in full
    driver_rows, his_rows = [], []
    hi = AA.index("H")
    for prof in pools + ["standard"]:
        f = {"plasma": f_pl, "erythrocytes": 1 - f_pl} if prof == "standard" else {prof: 1.0}
        contrib = {}
        for p, fp in f.items():
            for a, w in W[p].items():
                if w > 0 and a in comp:
                    contrib[a] = contrib.get(a, np.zeros(len(AA))) + 100 * fp * w * comp[a]["free"]
        for i, aa in enumerate(AA):
            total = sum(c[i] for c in contrib.values())
            top = sorted(contrib.items(), key=lambda kv: -kv[1][i])[:5]
            driver_rows.append([prof, aa, fnum(total)] + [f"{entries[a]['gene']} {fnum(c[i])} ({fnum(c[i] / total) if total else ''})" for a, c in top])
        total_h = sum(c[hi] for c in contrib.values())
        cum = 0.0
        for a, c in sorted(contrib.items(), key=lambda kv: -kv[1][hi])[:30]:
            cum += c[hi]
            his_rows.append([prof, a, entries[a]["gene"], fnum(sum(f[p] * W[p].get(a, 0.0) for p in f)), fnum(100 * comp[a]["free"][hi]), fnum(c[hi]),
                             fnum(c[hi] / total_h if total_h else 0.0), fnum(cum / total_h if total_h else 0.0)])
    outputs[out / "drivers_top_entries_per_amino_acid.tsv"] = tsv_text(header + ["free convention; per profile and amino acid, the five entries contributing most g/100 g: gene, contribution, share of the value"],
                                                                      ["profile", "amino_acid", "value", "top_1", "top_2", "top_3", "top_4", "top_5"], driver_rows)
    outputs[out / "histidine_drivers.tsv"] = tsv_text(header + ["free convention; histidine: the thirty entries contributing most, with their weight in the profile, their own histidine content and their share of the profile's histidine"],
                                                     ["profile", "accession", "gene", "weight_in_profile", "histidine_g_per_100g_of_entry", "contribution_g_per_100g", "share_of_histidine", "cumulative_share"], his_rows)

    # -- summary
    his_ctx = {}
    muscle_std = common.standard_dir() / "_calculated_amino_acid_standard.tsv"
    if muscle_std.exists():
        for r in read_commented_tsv(muscle_std):
            if r.get("amino_acid") == "H":
                his_ctx.update({f"muscle_{k}": v for k, v in r.items() if k != "amino_acid"})
    gor = common.COMPARISON_OUT_DIR / "calculated_vs_gorissen_2018.tsv"
    if gor.exists():
        for r in read_commented_tsv(gor):
            if r.get("amino_acid") == "H" or r.get(list(r)[0]) == "H":
                his_ctx["muscle_vs_gorissen_row"] = "\t".join(r.values())
    summary = common.new_ini()
    summary["inputs"] = {**sha, "draws": str(draws), "seed": str(seed), "reference": mix["sex"]}
    summary["split"] = {"haemoglobin_share_as_measured": fnum(hb_base), "albumin_share_as_measured": fnum(measured_share.get("plasma", 0.0)),
                        "plasma_protein_g": fnum(base_split["plasma_protein_g"]), "haemoglobin_g": fnum(base_split["hb_mass_g"]),
                        "haemoglobin_g_blood_route": fnum(base_split["hb_mass_g_blood_route"]), "erythrocyte_protein_g": fnum(base_split["erythrocyte_protein_g"]),
                        "whole_blood_protein_g": fnum(base_split["whole_blood_protein_g"]), "f_plasma": fnum(f_pl), "f_erythrocytes": fnum(1 - f_pl)}
    for label, s_val, _ in mix["hb_shares"]:
        sp = split(mix, s_val)
        summary["split"][f"f_plasma_at_hb_share_{label}"] = f"{fnum(sp['f_plasma'])} (whole-blood protein {fnum(sp['whole_blood_protein_g'])} g)"
    summary["standard_free_convention"] = {a: fnum(S["free"][i]) for i, a in enumerate(AA)}
    summary["histidine"] = {
        **{f"{p}_g_per_100g": fnum(P[(p, "free")][hi]) for p in pools}, "standard_g_per_100g": fnum(S["free"][hi]),
        "spread_dominant_protein_all": next(r[5] for r in spread_rows if r[0] == "free" and r[1] == "H"),
        "spread_haemoglobin_only": next(r[6] for r in spread_rows if r[0] == "free" and r[1] == "H"),
        "spread_albumin_only": next(r[7] for r in spread_rows if r[0] == "free" and r[1] == "H"),
        "spread_across_donors_erythrocytes": fnum(donor_spread.get(("erythrocytes", "free", "H"), (0.0, 0.0))[1]),
        "monte_carlo_half_width_standard_spread": fnum(half[("standard", "free", "spread", "H")]),
        "monte_carlo_half_width_standard_mean": fnum(half[("standard", "free", "mean_uncertainty", "H")]),
        "widest_term_standard": widest[("standard", "free", "H")],
        "haemoglobin_entries_share_of_standard_histidine": fnum(sum(float(r[6]) for r in his_rows if r[0] == "standard" and r[1] in dom_sets["erythrocytes"])),
        **his_ctx,
    }
    outputs[out / "standard_summary.ini"] = common.render_ini(summary, header)

    if make_plots:
        images.update(plot_all(out, pools, P, S, sens_profiles, ery_shares, pl_shares, terms_rows, term_names, hb_base))
    return outputs, images


# ----------------------------------------------------------------------------- plots

def plot_all(out: Path, pools, P, S, sens_profiles, ery_shares, pl_shares, terms_rows, term_names, hb_base) -> dict[Path, bytes]:
    import io
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    images = {}
    x = np.arange(len(AA))

    def png(fig) -> bytes:
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=150, bbox_inches="tight"); plt.close(fig); return buf.getvalue()

    fig, ax = plt.subplots(figsize=(12, 4.5))
    n = len(pools) + 1
    wdt = 0.8 / n
    for k, p in enumerate(pools):
        ax.bar(x + (k - n / 2 + 0.5) * wdt, P[(p, "free")], wdt, label=p)
    ax.bar(x + (len(pools) - n / 2 + 0.5) * wdt, S["free"], wdt, label="standard", color="black")
    ax.set_xticks(x); ax.set_xticklabels(AA); ax.set_ylabel("g / 100 g protein (free)"); ax.set_title("blood: pool profiles and the standard at the base split"); ax.legend()
    images[out / "plots" / "profile_pools_and_standard.png"] = png(fig)

    fig, ax = plt.subplots(figsize=(12, 4.5))
    base = sens_profiles[("as_measured", "as_measured", "free")][2]
    for el, es in ery_shares:
        for pl, ps in pl_shares:
            st = sens_profiles[(el, pl, "free")][2]
            ax.plot(x, st - base, marker="o", ms=3, lw=1, label=f"Hb {el} ({es:.2f}) / ALB {pl} ({ps:.2f})")
    ax.axhline(0, color="black", lw=0.8); ax.set_xticks(x); ax.set_xticklabels(AA); ax.set_ylabel("standard minus base, g / 100 g")
    ax.set_title("D117: the standard under each dominant-protein share"); ax.legend(fontsize=7, ncol=2)
    images[out / "plots" / "sensitivity_dominant_protein.png"] = png(fig)

    fig, ax = plt.subplots(figsize=(12, 4.5))
    rows = [r for r in terms_rows if r[0] == "standard" and r[1] == "free"]
    bottom = np.zeros(len(AA))
    for j, t in enumerate(term_names):
        vals = np.array([float(r[4 + j]) for r in rows])
        ax.bar(x, vals, 0.7, bottom=bottom, label=t)
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(AA); ax.set_ylabel("max absolute shift, g / 100 g (terms stacked)")
    ax.set_title("uncertainty terms of the blood standard, one scale"); ax.legend(fontsize=7, ncol=3)
    images[out / "plots" / "uncertainty_terms_standard.png"] = png(fig)
    return images


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--category", required=True)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare the tables with the files on disk (plots not compared)")
    args = ap.parse_args(argv)
    log = common.make_logger(f"{STAGE}_{args.category}" + ("_check" if args.check else ""))
    outputs, images = build(args.category, log, make_plots=not args.check)
    if args.check:
        stale = [p for p, t in outputs.items() if not p.exists() or strip_generated_line(p.read_text(encoding="utf-8")) != strip_generated_line(t)]
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), [p.name for p in stale])
            return 1
        log.info("check: the standard's outputs are current")
        return 0
    for p, t in outputs.items():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t, encoding="utf-8", newline="\n")
    for p, b in images.items():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b)
    log.info("wrote %d tables and %d plots under %s", len(outputs), len(images), common.category_outputs(args.category)["standard"].relative_to(common.REPO_ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
