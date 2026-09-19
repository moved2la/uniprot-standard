#!/usr/bin/env python3
"""
composite.py — offline stage: the composite standard, the category standards mixed by protein mass (Layer C).

Every category standard is a percent of that category's amino acid mass. The composite puts the
categories on one denominator: each profile is weighted by the protein the category holds in the
reference person (config/tissue_mass_fractions.ini, D67 / D121), and the weighted profiles are
summed. With muscle and blood it is two categories of the reference body, not the whole body; the
header of every output lists which categories are in it (D121).

    composite[a] = Σ_c f_c × profile_c[a],   f_c = P_c / Σ P   (L1)

Reads
  config/tissue_mass_fractions.ini                 hand-written: one [category.<name>] per category standard, the
                                                   protein mass of each and where it comes from (Layer C)
  the category standards it names                  outputs/standard/..., outputs/blood/standard/... (free convention)
  the residue-convention tables it names           for the amino-acid-mass line (L4)
  config/non_protein_metabolite_pools.ini          [protein_content] (D77) when a category says protein_content_from
  outputs/blood/standard/compartment_split.tsv     a category whose own stage computed its protein mass (blood, D100 / D119)
  outputs/blood/standard/sensitivity_dominant_protein.tsv   the blood profile at every published share (D117), for L3
  outputs/standard/non_protein_metabolite_pool_amounts_per_kg_muscle.tsv, config/fiber_type_mix.ini   the pool's mass (L4)
  data/iupac/amino_acid_symbols.ini, data/pubchem/amino_acid_masses.ini

Writes (outputs/composite/)
  _calculated_amino_acid_standard_composite.tsv   the deliverable: one column per category and the composite, free convention
  category_protein_masses.tsv                     Layer C: per category and sex, the tissue mass, the protein mass, how, the share
  amino_acid_profiles.tsv                         long form: every profile written above, and the composite at the female masses
  category_contributions.tsv                      per amino acid, the grams each category contributes per 100 g of composite
  eaa_subset.tsv                                  the FAO 2013 indispensable amino acids (D62) per category and composite
  sensitivity_category_split.tsv                  L3: the composite at every published share of a category's split (blood, D117)
  sensitivity_category_split_spread.tsv           its range per amino acid (male rows)
  sensitivity_weighting_basis.tsv                 L4: the composite weighted by amino acid mass beside the protein-mass composite
  composite_summary.ini                           masses, shares, the composite, a [histidine] section, the L4 factors, hashes

Rules (docs/conventions.md, "Composite rules")
  L1 Weights are protein masses, from a cited tissue mass × a cited protein content, or from the category's own split.
  L2 Reference male; the female masses recorded and reported as a line, never the base (D100 / D119 pattern).
  L3 A category with a published-share split is reported at every row of that split, its mass moving with the share (D122).
  L4 Mixing by protein mass assumes one amino-acid mass per gram of protein across categories; the stage computes each
     category's factor from its residue-convention table and reports the composite under amino-acid-mass weights beside it.

Usage
  python pipeline/composite.py            # build
  python pipeline/composite.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline import common
from pipeline.aggregate_pools import load_eaa, read_commented_tsv
from pipeline.mass_fractions_pools import fnum, strip_generated_line, tsv_text

AA = list(common.AMINO_ACIDS)
CONFIG_INI = common.TISSUE_MASS_FRACTIONS_INI
OUT_DIR = common.COMPOSITE_OUT_DIR
TOOL = "composite.py"
SUM_TOLERANCE = 0.01          # a category table written at four decimals sums to 100 within this


class Stop(SystemExit):
    def __init__(self, msg: str):
        super().__init__(f"[STOP] {msg}")


def _rel(p: Path) -> str:
    try:
        return p.relative_to(common.REPO_ROOT).as_posix()
    except ValueError:
        return p.name


def _path(rel: str) -> Path:
    p = Path(rel.strip())
    return p if p.is_absolute() else common.REPO_ROOT / p


def _need(sec, key: str, where: str) -> str:
    v = sec.get(key, "").strip()
    if v in ("", "___"):
        raise Stop(f"{where} {key} is not filled in")
    return v


# --------------------------------------------------------------------------- inputs

def read_profile(path: Path, column: str) -> dict[str, float]:
    """{letter: g per 100 g} from a standard table's named column; the twenty letters, summing to 100."""
    if not path.exists():
        raise Stop(f"{_rel(path)} not found — its category's standard has not been built")
    rows = read_commented_tsv(path)
    if not rows or "amino_acid" not in rows[0] or column not in rows[0]:
        raise Stop(f"{_rel(path)} has no column {column!r} beside amino_acid")
    out = {r["amino_acid"]: float(r[column]) for r in rows if r["amino_acid"] in AA and r[column].strip()}
    missing = [a for a in AA if a not in out]
    if missing:
        raise Stop(f"{_rel(path)} column {column} is missing {', '.join(missing)}")
    s = sum(out.values())
    if abs(s - 100.0) > SUM_TOLERANCE:
        raise Stop(f"{_rel(path)} column {column} sums to {s:.4f}, not 100")
    return out


def read_symbols() -> dict[str, tuple[str, str]]:
    cp = common.read_ini(common.AMINO_ACID_SYMBOLS_INI)
    return {a: (cp[a]["three_letter"].strip(), cp[a]["trivial_name"].strip()) for a in AA}


def read_masses() -> dict[str, dict[str, float]]:
    cp = common.read_ini(common.AMINO_ACID_MASSES_INI)
    return {a: {"free": float(cp[a]["free_mass"]), "residue": float(cp[a]["residue_mass"])} for a in AA}


def protein_content(spec: str) -> dict:
    """`<config path> [<section>]` -> the section's value (g/kg), spread, basis, source and location."""
    spec = spec.strip()
    if "[" not in spec or not spec.endswith("]"):
        raise Stop(f"protein_content_from must read `<path> [<section>]`, not {spec!r}")
    rel, section = spec[:spec.index("[")].strip(), spec[spec.index("[") + 1:-1].strip()
    path = _path(rel)
    if not path.exists():
        raise Stop(f"{rel} not found")
    cp = common.read_ini(path)
    if section not in cp:
        raise Stop(f"{rel} has no [{section}]")
    s = cp[section]
    where = f"{rel} [{section}]"
    return {"path": path, "section": section, "g_per_kg": float(_need(s, "value", where)),
            "spread": float(s.get("spread", "0") or 0), "basis": _need(s, "basis", where),
            "source_id": s.get("source_id", "").strip(), "location": s.get("location", "").strip(),
            "spread_kind": s.get("spread_kind", "").strip()}


def read_split(path: Path, column: str) -> list[dict[str, str]]:
    if not path.exists():
        raise Stop(f"{_rel(path)} not found — its category's standard has not been built")
    rows = read_commented_tsv(path)
    if not rows or column not in rows[0]:
        raise Stop(f"{_rel(path)} has no column {column!r}")
    return rows


def fiber_mix() -> dict[str, float]:
    cp = common.read_ini(common.FIBER_TYPE_MIX_INI)
    if "shares" not in cp:
        raise Stop(f"{_rel(common.FIBER_TYPE_MIX_INI)} has no [shares]")
    out = {}
    for k, v in cp["shares"].items():
        try:
            out[k] = float(v)
        except ValueError:
            continue                      # basis, location: text keys beside the shares
    if abs(sum(out.values()) - 1.0) > 1e-6:
        raise Stop(f"{_rel(common.FIBER_TYPE_MIX_INI)} [shares] sum to {sum(out.values())}, not 1")
    return out


def pool_g_per_g_protein(amounts_path: Path, mix_path: Path) -> tuple[float, str]:
    """The non-protein metabolite pool's amino acid mass per gram of protein, per fiber type from the
    amounts table, mixed by the fiber-type shares (A10's inputs, read back; basis cancels in the ratio)."""
    if not amounts_path.exists():
        raise Stop(f"{_rel(amounts_path)} not found — run `python run.py standard` first")
    rows = read_commented_tsv(amounts_path)
    need = {"fiber_type", "pool_His_g_per_kg", "protein_g_per_kg"}
    if not rows or not need <= set(rows[0]):
        raise Stop(f"{_rel(amounts_path)} lacks {sorted(need - set(rows[0] if rows else []))}")
    mix = fiber_mix()
    total, parts = 0.0, []
    for r in rows:
        ft = r["fiber_type"]
        if ft not in mix:
            raise Stop(f"{_rel(amounts_path)}: fiber type {ft!r} has no share in {_rel(mix_path)}")
        ratio = float(r["pool_His_g_per_kg"]) / float(r["protein_g_per_kg"])
        total += mix[ft] * ratio
        parts.append(f"{ft} {ratio:.6f} x {mix[ft]:g}")
    return total, "; ".join(parts)


def free_mass_factor(residue_profile: dict[str, float], masses: dict) -> float:
    """Grams of free amino acids per gram of protein, from the residue-convention profile: each residue
    becomes its free amino acid by taking back a water (C1: residue mass = free mass − water)."""
    return sum(residue_profile[a] / 100.0 * masses[a]["free"] / masses[a]["residue"] for a in AA)


# --------------------------------------------------------------------------- categories

def load_categories(cp) -> list[dict]:
    cats = []
    for sec in cp.sections():
        if not sec.startswith("category."):
            continue
        name, s, where = sec[len("category."):], cp[sec], f"[{sec}]"
        c = {"name": name, "label": _need(s, "label", where),
             "standard_file": _path(_need(s, "standard_file", where)), "standard_column": _need(s, "standard_column", where),
             "residue_file": _path(_need(s, "residue_file", where)), "residue_column": _need(s, "residue_column", where),
             "tissue_mass_g": float(_need(s, "tissue_mass_g", where)),
             "tissue_mass_g_female": float(s["tissue_mass_g_female"]) if s.get("tissue_mass_g_female", "").strip() else None,
             "tissue_mass_location": _need(s, "tissue_mass_location", where),
             "how": _need(s, "protein_mass", where), "inputs": {}}
        c["profile"] = read_profile(c["standard_file"], c["standard_column"])
        c["residue_profile"] = read_profile(c["residue_file"], c["residue_column"])
        c["inputs"][f"standard.{name}"] = c["standard_file"]
        c["inputs"][f"residue.{name}"] = c["residue_file"]
        if c["how"] == "tissue_mass_x_protein_content":
            pc = protein_content(_need(s, "protein_content_from", where))
            basis = _need(s, "protein_content_basis", where)
            if pc["basis"] != basis:
                raise Stop(f"{where}: protein_content_basis = {basis} but {pc['path'].name} [{pc['section']}] basis = {pc['basis']}")
            c["protein_content"] = pc
            c["inputs"][f"protein_content.{name}"] = pc["path"]
            c["protein_g"] = c["tissue_mass_g"] / 1000.0 * pc["g_per_kg"]
            c["protein_g_female"] = (c["tissue_mass_g_female"] / 1000.0 * pc["g_per_kg"]) if c["tissue_mass_g_female"] is not None else None
            c["mass_how"] = (f"{c['tissue_mass_g']:g} g x {pc['g_per_kg']:g} g/kg ({pc['basis']}) / 1000")
            c["mass_location"] = f"{pc['source_id']}: {pc['location']}"
        elif c["how"] == "category_split":
            split_path, col = _path(_need(s, "split_file", where)), _need(s, "split_column", where)
            base_row, fem = _need(s, "split_base_row", where), s.get("split_female_prefix", "female").strip()
            rows = read_split(split_path, col)
            key = list(rows[0].keys())[0]                                   # the row-label column
            by_label = {r[key]: r for r in rows}
            if base_row not in by_label:
                raise Stop(f"{_rel(split_path)} has no row {base_row!r} in column {key}")
            c["split"] = {"path": split_path, "column": col, "label_column": key, "rows": rows, "base_row": base_row, "female_prefix": fem}
            c["inputs"][f"split.{name}"] = split_path
            c["protein_g"] = float(by_label[base_row][col])
            f_row = f"{fem} {base_row}"
            c["protein_g_female"] = float(by_label[f_row][col]) if f_row in by_label else None
            c["mass_how"] = f"{_rel(split_path)} row {base_row}, column {col}"
            c["mass_location"] = by_label[base_row].get("share_source", "")
            if s.get("sensitivity_file", "").strip():
                sp = _path(s["sensitivity_file"])
                if not sp.exists():
                    raise Stop(f"{_rel(sp)} not found")
                c["sensitivity"] = {"path": sp, "row_column": _need(s, "sensitivity_row_column", where),
                                    "variant_column": _need(s, "sensitivity_variant_column", where),
                                    "profile_column": _need(s, "sensitivity_profile_column", where),
                                    "row_share_column": _need(s, "sensitivity_row_share_column", where),
                                    "variant_share_column": _need(s, "sensitivity_variant_share_column", where),
                                    "note": s.get("sensitivity_note", "").strip()}
                c["inputs"][f"sensitivity.{name}"] = sp
        else:
            raise Stop(f"{where} protein_mass must be tissue_mass_x_protein_content or category_split, not {c['how']!r}")
        if s.get("non_protein_metabolite_pool_amounts", "").strip():
            ap, mp = _path(s["non_protein_metabolite_pool_amounts"]), _path(_need(s, "fiber_type_mix", where))
            c["pool_g_per_g_protein"], c["pool_how"] = pool_g_per_g_protein(ap, mp)
            c["inputs"][f"pool_amounts.{name}"] = ap
            c["inputs"][f"fiber_type_mix.{name}"] = mp
        else:
            c["pool_g_per_g_protein"], c["pool_how"] = 0.0, ""
        cats.append(c)
    if len(cats) < 2:
        raise Stop(f"{CONFIG_INI.name} names {len(cats)} [category.*]; a composite needs at least two")
    return cats


def mix(cats: list[dict], masses_g: dict[str, float]) -> tuple[dict[str, float], dict[str, float]]:
    """L1: shares from the masses given, composite = Σ share × profile."""
    total = sum(masses_g.values())
    if total <= 0:
        raise Stop("the category protein masses sum to zero")
    f = {c["name"]: masses_g[c["name"]] / total for c in cats}
    comp = {a: sum(f[c["name"]] * c["profile"][a] for c in cats) for a in AA}
    return f, comp


# --------------------------------------------------------------------------- build

def build(log) -> dict[str, str]:
    if not CONFIG_INI.exists():
        raise Stop(f"{_rel(CONFIG_INI)} not found")
    cp = common.read_ini(CONFIG_INI)
    sex = cp["meta"].get("reference_sex", "male").strip() if "meta" in cp else "male"
    cats = load_categories(cp)
    sym, masses = read_symbols(), read_masses()
    eaa = load_eaa()
    names = [c["name"] for c in cats]

    P = {c["name"]: c["protein_g"] for c in cats}
    f, comp = mix(cats, P)
    for c in cats:
        log.info("%s: protein %.1f g (%s); share of composite %.4f", c["name"], c["protein_g"], c["mass_how"], f[c["name"]])
    female_ok = all(c["protein_g_female"] is not None for c in cats)
    P_f = {c["name"]: c["protein_g_female"] for c in cats} if female_ok else None
    f_f, comp_f = mix(cats, P_f) if female_ok else (None, None)

    # -- provenance header: what was read, hashed; the categories with their masses (D121)
    inputs = {"config": CONFIG_INI, "symbols": common.AMINO_ACID_SYMBOLS_INI, "masses": common.AMINO_ACID_MASSES_INI,
              "eaa": common.EAA_INI}
    for c in cats:
        inputs.update(c["inputs"])
    seen, hash_lines = set(), []
    for k, p in inputs.items():
        if p.exists() and p not in seen:
            seen.add(p)
            hash_lines.append(f"input.{k} = {_rel(p)} sha256 {common.sha256_file(p)}")
    header = ([f"GENERATED by {TOOL}. DO NOT EDIT BY HAND.", f"generated = {common.iso_now()}",
               f"categories = {', '.join(names)} (reference {sex}; L1: mixed by protein mass, config/tissue_mass_fractions.ini)"]
              + hash_lines
              + [f"category.{c['name']} = {c['label']} | {_rel(c['standard_file'])} column {c['standard_column']} | protein {fnum(c['protein_g'])} g"
                 f" = {c['mass_how']} | share {fnum(f[c['name']])}" for c in cats])

    files: dict[str, str] = {}

    # ---- 1. the deliverable
    rows = [[a, sym[a][0], sym[a][1]] + [fnum(c["profile"][a]) for c in cats] + [fnum(comp[a])] for a in AA]
    files["_calculated_amino_acid_standard_composite.tsv"] = tsv_text(
        header + ["g per 100 g amino acid mass, free convention; one column per category standard as read, and the composite (L1)"],
        ["amino_acid", "three_letter", "name"] + names + ["composite"], rows)

    # ---- 2. Layer C: the masses
    mass_rows = []
    for c in cats:
        mass_rows.append([c["name"], sex, fnum(c["tissue_mass_g"]), fnum(c["protein_g"]), fnum(f[c["name"]]), c["how"], c["mass_how"],
                          c["tissue_mass_location"], c["mass_location"]])
    if female_ok:
        for c in cats:
            mass_rows.append([c["name"], "female", fnum(c["tissue_mass_g_female"]), fnum(c["protein_g_female"]), fnum(f_f[c["name"]]), c["how"],
                              "female tissue mass, or the female row of the split; reported, not the base (L2)", c["tissue_mass_location"], c["mass_location"]])
    files["category_protein_masses.tsv"] = tsv_text(
        header + ["Layer C (D67, D121): per category and sex, the cited tissue mass, the protein mass the composite weights by, and its share"],
        ["category", "sex", "tissue_mass_g", "protein_g", "share_of_composite", "protein_mass_rule", "protein_mass_how", "tissue_mass_location", "protein_mass_location"],
        mass_rows)

    # ---- 3. long form
    prof_rows = [[c["name"], a, fnum(c["profile"][a])] for c in cats for a in AA] + [["composite", a, fnum(comp[a])] for a in AA]
    if female_ok:
        prof_rows += [["composite_female_masses", a, fnum(comp_f[a])] for a in AA]
    files["amino_acid_profiles.tsv"] = tsv_text(header + ["free convention; composite_female_masses = the same profiles mixed by the female masses (L2)"],
                                                ["profile", "amino_acid", "g_per_100g"], prof_rows)

    # ---- 4. contributions: the grams each category puts into 100 g of composite, per amino acid
    contrib_rows = []
    for a in AA:
        parts = [f[c["name"]] * c["profile"][a] for c in cats]
        contrib_rows.append([a, sym[a][0]] + [fnum(x) for x in parts] + [fnum(comp[a])]
                            + [fnum(x / comp[a]) if comp[a] else "" for x in parts])
    files["category_contributions.tsv"] = tsv_text(
        header + ["per amino acid: g contributed by each category per 100 g of composite (share x profile), their sum, and each category's fraction of the amino acid"],
        ["amino_acid", "three_letter"] + [f"g_from_{n}" for n in names] + ["composite"] + [f"fraction_from_{n}" for n in names], contrib_rows)

    # ---- 5. EAA subset (D62)
    eaa_rows = []
    for heading, letters in eaa["headings"]:
        eaa_rows.append([heading, "+".join(letters)] + [fnum(sum(c["profile"][x] for x in letters)) for c in cats] + [fnum(sum(comp[x] for x in letters))])
    eaa_rows.append(["sum", ""] + [fnum(sum(float(r[i + 2]) for r in eaa_rows)) for i in range(len(cats) + 1)])
    files["eaa_subset.tsv"] = tsv_text(header + [f"FAO 2013 indispensable amino acids as transcribed (D62; {eaa['location']}); free convention"],
                                       ["heading", "amino_acids"] + names + ["composite"], eaa_rows)

    # ---- 6. L3: the composite at every published share of a category's split (blood, D117)
    sens_rows, sens_vals = [], {}
    for c in cats:
        if "sensitivity" not in c:
            continue
        sp, sv = c["split"], c["sensitivity"]
        srows = read_commented_tsv(sv["path"])
        need = {sv["row_column"], sv["variant_column"], sv["profile_column"], sv["row_share_column"], sv["variant_share_column"], "amino_acid", "convention"}
        if not srows or not need <= set(srows[0]):
            raise Stop(f"{_rel(sv['path'])} lacks {sorted(need - set(srows[0] if srows else []))}")
        profiles: dict[tuple[str, str], dict[str, float]] = {}
        shares: dict[tuple[str, str], tuple[str, str]] = {}
        for r in srows:
            if r["convention"] != "free":
                continue
            k = (r[sv["row_column"]], r[sv["variant_column"]])
            profiles.setdefault(k, {})[r["amino_acid"]] = float(r[sv["profile_column"]])
            shares[k] = (r[sv["row_share_column"]], r[sv["variant_share_column"]])
        if not profiles:
            raise Stop(f"{_rel(sv['path'])} has no free-convention rows")
        for split_row in sp["rows"]:
            label = split_row[sp["label_column"]]
            female = label.startswith(sp["female_prefix"] + " ")
            row_key = label[len(sp["female_prefix"]) + 1:] if female else label
            if female and not female_ok:
                continue
            masses_here = dict(P_f if female else P)
            masses_here[c["name"]] = float(split_row[sp["column"]])
            for (rk, vk), prof in sorted(profiles.items()):
                if rk != row_key:
                    continue
                if any(a not in prof for a in AA):
                    raise Stop(f"{_rel(sv['path'])}: row {rk} / {vk} lacks an amino acid")
                total = sum(masses_here.values())
                f_here = {n: masses_here[n] / total for n in masses_here}
                for a in AA:
                    val = sum(f_here[o["name"]] * (prof[a] if o["name"] == c["name"] else o["profile"][a]) for o in cats)
                    base = comp_f[a] if female else comp[a]
                    sens_rows.append([c["name"], "female" if female else sex, rk, shares[(rk, vk)][0], vk, shares[(rk, vk)][1],
                                      fnum(masses_here[c["name"]]), fnum(f_here[c["name"]]), a, fnum(prof[a]), fnum(val), fnum(val - base)])
                    if not female:
                        sens_vals.setdefault(a, []).append(val)
    files["sensitivity_category_split.tsv"] = tsv_text(
        header + ["L3 (D117, D122): a category's profile at each published share of its split, its protein mass at that share, the composite recomputed; composite_minus_base against the base of the same sex"],
        ["category", "sex", "split_row", "split_share", "variant", "variant_share", "category_protein_g", "category_share_of_composite",
         "amino_acid", "profile_category", "composite", "composite_minus_base"], sens_rows)
    spread_rows = [[a, fnum(comp[a]), fnum(min(v)), fnum(max(v)), fnum(max(v) - min(v))] for a, v in ((a, sens_vals.get(a)) for a in AA) if v]
    files["sensitivity_category_split_spread.tsv"] = tsv_text(header + ["the composite's range over every L3 row of the reference sex"],
                                                              ["amino_acid", "composite_base", "min", "max", "spread"], spread_rows)

    # ---- 7. L4: weighting by amino acid mass instead of protein mass
    factors, aa_mass = {}, {}
    for c in cats:
        factors[c["name"]] = free_mass_factor(c["residue_profile"], masses)
        aa_mass[c["name"]] = c["protein_g"] * (factors[c["name"]] + c["pool_g_per_g_protein"])
    f_aa, comp_aa = mix(cats, aa_mass)
    wb_rows = [[a, fnum(comp[a]), fnum(comp_aa[a]), fnum(comp_aa[a] - comp[a])] for a in AA]
    files["sensitivity_weighting_basis.tsv"] = tsv_text(
        header + ["L4: the composite with each category weighted by its amino acid mass (protein x free g per g protein from its residue-convention table, plus its non-protein metabolite pool) beside the protein-mass composite; "
                  + "; ".join(f"{n} factor {factors[n]:.4f}" + (f" + pool {c['pool_g_per_g_protein']:.4f} g/g" if c["pool_g_per_g_protein"] else "") + f" -> {fnum(aa_mass[n])} g, share {fnum(f_aa[n])}"
                              for n, c in zip(names, cats))],
        ["amino_acid", "composite_by_protein_mass", "composite_by_amino_acid_mass", "difference"], wb_rows)

    # ---- 8. summary
    s = common.new_ini()
    s.add_section("run")
    s["run"]["categories"] = ", ".join(names)
    s["run"]["reference_sex"] = sex
    s["run"]["weighting"] = "protein mass (L1)"
    s.add_section("masses")
    for c in cats:
        s["masses"][f"{c['name']}_tissue_g"] = fnum(c["tissue_mass_g"])
        s["masses"][f"{c['name']}_protein_g"] = fnum(c["protein_g"])
        s["masses"][f"{c['name']}_share"] = fnum(f[c["name"]])
    s["masses"]["total_protein_g"] = fnum(sum(P.values()))
    if female_ok:
        for c in cats:
            s["masses"][f"{c['name']}_protein_g_female"] = fnum(c["protein_g_female"])
            s["masses"][f"{c['name']}_share_female"] = fnum(f_f[c["name"]])
    s.add_section("composite_free_convention")
    for a in AA:
        s["composite_free_convention"][a] = fnum(comp[a])
    s.add_section("histidine")
    for c in cats:
        s["histidine"][f"{c['name']}_g_per_100g"] = fnum(c["profile"]["H"])
        s["histidine"][f"{c['name']}_g_in_100g_composite"] = fnum(f[c["name"]] * c["profile"]["H"])
        s["histidine"][f"{c['name']}_fraction_of_composite_histidine"] = fnum(f[c["name"]] * c["profile"]["H"] / comp["H"])
    s["histidine"]["composite_g_per_100g"] = fnum(comp["H"])
    if female_ok:
        s["histidine"]["composite_g_per_100g_female_masses"] = fnum(comp_f["H"])
    if sens_vals.get("H"):
        s["histidine"]["spread_across_split_rows"] = fnum(max(sens_vals["H"]) - min(sens_vals["H"]))
        s["histidine"]["min_across_split_rows"] = fnum(min(sens_vals["H"]))
        s["histidine"]["max_across_split_rows"] = fnum(max(sens_vals["H"]))
    s["histidine"]["difference_by_amino_acid_mass_weighting"] = fnum(comp_aa["H"] - comp["H"])
    s.add_section("phenylalanine")
    for c in cats:
        s["phenylalanine"][f"{c['name']}_g_per_100g"] = fnum(c["profile"]["F"])
    s["phenylalanine"]["composite_g_per_100g"] = fnum(comp["F"])
    if sens_vals.get("F"):
        s["phenylalanine"]["spread_across_split_rows"] = fnum(max(sens_vals["F"]) - min(sens_vals["F"]))
    s.add_section("weighting_basis")
    for c in cats:
        s["weighting_basis"][f"{c['name']}_free_g_per_g_protein"] = fnum(factors[c["name"]])
        s["weighting_basis"][f"{c['name']}_pool_g_per_g_protein"] = fnum(c["pool_g_per_g_protein"]) + (f"  ; {c['pool_how']}" if c["pool_how"] else "")
        s["weighting_basis"][f"{c['name']}_amino_acid_mass_g"] = fnum(aa_mass[c["name"]])
        s["weighting_basis"][f"{c['name']}_share_by_amino_acid_mass"] = fnum(f_aa[c["name"]])
    s["weighting_basis"]["largest_absolute_difference"] = fnum(max(abs(comp_aa[a] - comp[a]) for a in AA)) + "  ; g per 100 g, over the twenty amino acids"
    s.add_section("inputs")
    for line in hash_lines:
        k, v = line.split(" = ", 1)
        s["inputs"][k[len("input."):]] = v
    files["composite_summary.ini"] = common.render_ini(s, header)
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("composite_check" if args.check else "composite")
    files = build(log)

    if args.check:
        differ = [n for n, t in files.items()
                  if not (OUT_DIR / n).exists()
                  or strip_generated_line((OUT_DIR / n).read_text(encoding="utf-8")) != strip_generated_line(t)]
        if differ:
            log.error("check: %d file(s) differ from a fresh build: %s", len(differ), differ)
            return 1
        log.info("check: the composite outputs are current")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        common.write_text_file(OUT_DIR / name, text)
        log.info("wrote %s", _rel(OUT_DIR / name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
