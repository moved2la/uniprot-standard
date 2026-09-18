#!/usr/bin/env python3
"""
comparison.py — offline stage: the calculated standard beside an independent measurement.

Gorissen et al. 2018 measured the amino acid composition of human skeletal muscle protein
(m. vastus lateralis, n = 10, acid hydrolysis, UPLC-MS/MS) and published it as one column of
Table 1. It is the reference the original Match Rate was calibrated against, and it is the only
measurement of the same thing this repository derives from sequence. This stage puts the two
side by side and computes the difference. It reaches no verdict: neither measurement is treated
as the truth the other is scored against (the D54 pattern).

Two things make the comparison narrower than twenty amino acids, and both come from the paper's
own method rather than from a choice made here:

  - Tryptophan is decomposed by acid hydrolysis, and aspartic acid is listed as not measured.
    Asparagine converts into aspartic acid, so it is not measurable either. Three of the twenty
    (W, D, N) are therefore absent from the measured side and are excluded from BOTH sides.
  - Two further rows print as 0.0 in the human muscle column (proline, cysteine) and are declared
    not measured in the config, with the reason (D85, rule X5); they too are excluded from both.
  - Glutamine converts into glutamic acid during hydrolysis, so the paper's glutamic acid row
    holds both. The calculated standard's E and Q are summed into one row to match it.

The set is whatever the config declares under those rules (fourteen rows across fifteen amino
acids for Gorissen 2018). Each side is then renormalised over that set, so the two columns are
percentages of the same thing. Comparing
the paper's own "% of total protein" figures with this repository's would not be like for like:
the paper's denominator is a protein mass that includes the amino acids it did not measure.

Reads
  config/gorissen_2018_comparison.ini                    the transcribed Table 1 column (hand-written, D79/D84)
  config/literature_sources.ini                          the [source.gorissen_2018] citation
  data/literature/manifest.ini                           the green light: the PDF must be on disk and hashed (D78)
  data/iupac/amino_acid_symbols.ini                      trivial name -> letter: nothing is named in code
  outputs/standard/_calculated_amino_acid_standard.tsv                              protein only
  outputs/standard/_calculated_amino_acid_standard_with_non_protein_metabolite_pools.tsv

Writes (outputs/comparison/)
  gorissen_2018_human_muscle.tsv          the transcribed column as letters, with the location of every value
  calculated_vs_gorissen_2018.tsv         the comparison: both sides renormalised over the common set, difference and ratio
  gorissen_2018_mass_balance.tsv          whether the published values reconcile with a stated protein content (X6)
  comparison_summary.ini                  what was excluded and why, the transcription checks, the largest differences

Rules (docs/conventions.md, "Comparison rules", X1-X6)
  X1 Both sides renormalised over the SAME set, and the set is decided by the measured side's own
     stated method, never by which amino acids happen to agree.
  X2 An amino acid the method converts into another is compared as the sum of the two, on both
     sides. An amino acid the method destroys or does not report is excluded from both.
  X3 No verdict. The stage computes difference and ratio; it does not call either side correct,
     does not fit, scale, or adjust either side toward the other, and writes no "error".
  X4 Transcription is checked against the paper's own printed sums before anything is computed;
     a mismatch is a [STOP], not a flag.
  X5 A zero is never compared as a measurement by default. A measured row of 0.0 must be either
     declared in [not_measured] with its reason, or explicitly supported as a measurement in
     [zero_is_a_measurement]; an undeclared zero is a [STOP]. A zero that reaches the comparison
     unexamined is the failure this rule exists to prevent: it sits in one side's denominator and
     not the other's, and biases every other row.
  X6 The mass balance is computed, not assumed. Hydrolysis adds a water per peptide bond, so a gram
     of protein yields more than a gram of free amino acids; the stage computes that factor from
     this repository's own composition and the PubChem masses, and reports what share of the
     expected yield the published values account for, under each protein content offered. It draws
     no conclusion from the result.

Usage
  python pipeline/comparison.py            # build
  python pipeline/comparison.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline import common
from pipeline.aggregate import strip_generated_line, tsv_text
from pipeline.mass_fractions import sha256_path

AA = list(common.AMINO_ACIDS)
PLACEHOLDER = "___"
OUT_DIR = common.COMPARISON_OUT_DIR
SOURCE_ID = "gorissen_2018"

INPUTS = {
    "comparison": common.CONFIG_DIR / "gorissen_2018_comparison.ini",
    "sources": common.LITERATURE_SOURCES_INI,
    "symbols": common.AMINO_ACID_SYMBOLS_INI,
    "standard": common.standard_path("_calculated_amino_acid_standard.tsv"),
    "standard_npmp": common.standard_path("_calculated_amino_acid_standard_with_non_protein_metabolite_pools.tsv"),
    "masses": common.AMINO_ACID_MASSES_INI,
    "pools": common.NON_PROTEIN_METABOLITE_POOLS_INI,
}
MANIFEST = common.literature_manifest()


class Stop(SystemExit):
    def __init__(self, msg: str):
        super().__init__(f"[STOP] {msg}")


def _rel(p: Path) -> str:
    try:
        return p.relative_to(common.REPO_ROOT).as_posix()
    except ValueError:
        return p.name


# --------------------------------------------------------------------------- inputs

def trivial_to_letter() -> dict[str, str]:
    """{trivial name lowercased: letter} from the IUPAC table (rule U3's join, reused)."""
    cp = common.read_ini(INPUTS["symbols"])
    out = {}
    for a in AA:
        if a not in cp:
            raise Stop(f"amino acid {a} missing from {INPUTS['symbols'].name}")
        name = cp[a].get("trivial_name", "").strip().lower()
        if not name:
            raise Stop(f"{INPUTS['symbols'].name} [{a}] has no trivial_name")
        out[name] = a
    return out


def green_light(log) -> dict:
    """D78: the cited PDF must be on disk and hashed, or nothing is computed."""
    if not MANIFEST.exists():
        raise Stop(f"{_rel(MANIFEST)} not found — run the literature hashing before the comparison (D78)")
    man = common.read_ini(MANIFEST)
    hits = [s for s in man.sections() if man[s].get("source_id", "").strip() == SOURCE_ID]
    if not hits:
        raise Stop(f"no file for source {SOURCE_ID} in {_rel(MANIFEST)}; "
                   f"add [source.{SOURCE_ID}] to config/literature_sources.ini and run "
                   f"`python pipeline/fetch_literature.py --only {SOURCE_ID}` (D78)")
    sec = man[hits[0]]
    path = common.REPO_ROOT / sec.get("path", "").strip()
    if not path.exists():
        raise Stop(f"{sec.get('path')} is in the manifest but not on disk (D78)")
    got = sha256_path(path)
    if got != sec.get("sha256", "").strip().lower():
        raise Stop(f"{sec.get('path')} sha256 {got} != manifest {sec.get('sha256')} (B0)")
    log.info("green light: %s %s sha256 %s ok", SOURCE_ID, sec.get("path"), got[:12])
    return {"path": sec.get("path", "").strip(), "sha256": got,
            "published_name": sec.get("published_name", "").strip()}


def citation() -> str:
    cp = common.read_ini(INPUTS["sources"])
    sec = f"source.{SOURCE_ID}"
    if sec not in cp:
        raise Stop(f"[{sec}] not in {INPUTS['sources'].name}")
    c, doi = cp[sec].get("citation", "").strip(), cp[sec].get("doi", "").strip()
    return f"{cp[sec].get('id_short', '').strip()} {c} DOI {doi}"


def read_standard(path: Path, value_col: str) -> dict[str, float]:
    if not path.exists():
        raise Stop(f"{_rel(path)} not found — run `python run.py standard` first")
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    cols = lines[0].split("\t")
    if value_col not in cols:
        raise Stop(f"{_rel(path)} has no column {value_col!r}; found {cols}")
    i_aa, i_v = cols.index("amino_acid"), cols.index(value_col)
    out = {}
    for ln in lines[1:]:
        parts = ln.split("\t")
        if parts[i_aa] in AA and parts[i_v].strip():
            out[parts[i_aa]] = float(parts[i_v])
    missing = [a for a in AA if a not in out]
    if missing:
        raise Stop(f"{_rel(path)} is missing {', '.join(missing)}")
    return out


# --------------------------------------------------------------------------- the comparison set

def build_sets(cp, name_to_letter: dict[str, str]) -> tuple[dict, dict, dict]:
    """Rule X2. Returns (measured {row key: (letters, value)}, excluded {letter: why}, labels)."""
    if "values.human_muscle" not in cp:
        raise Stop(f"{INPUTS['comparison'].name} has no [values.human_muscle]")

    conversions: dict[str, list[str]] = {}
    if "hydrolysis_conversions" in cp:
        for key, val in cp["hydrolysis_conversions"].items():
            if key in ("conversion_basis", "location"):
                continue
            names = [n.strip() for n in val.split("+")]
            conversions[key.strip().lower()] = names

    supported_zero = set()
    if "zero_is_a_measurement" in cp:
        supported_zero = {k.strip().lower() for k in cp["zero_is_a_measurement"]}

    measured: dict[str, dict] = {}
    for row_label, raw in cp["values.human_muscle"].items():
        label = row_label.strip()
        if raw.strip() == PLACEHOLDER or raw.strip() == "":
            raise Stop(f"[values.human_muscle] {label} is not filled in")
        if float(raw) == 0.0 and label.lower() not in supported_zero:
            raise Stop(f"X5: [values.human_muscle] {label} is 0.0. A zero is not compared as a "
                       "measurement unless it is declared: move it to [not_measured] with the reason "
                       "it is not a measurement, or add it to [zero_is_a_measurement] with the "
                       "evidence that it is one. An undeclared zero biases every other row, because "
                       "it sits in one side's denominator and not the other's.")
        names = conversions.get(label.lower(), [label])
        letters = []
        for n in names:
            letter = name_to_letter.get(n.strip().lower())
            if letter is None:
                raise Stop(f"[values.human_muscle] row {label!r} -> {n!r} is not an IUPAC trivial name "
                           f"in {INPUTS['symbols'].name}")
            letters.append(letter)
        key = "+".join(letters)
        measured[key] = {"label": label, "letters": letters, "value": float(raw)}

    excluded: dict[str, str] = {}
    if "not_measured" in cp:
        for name, why in cp["not_measured"].items():
            letter = name_to_letter.get(name.strip().lower())
            if letter is None:
                raise Stop(f"[not_measured] {name!r} is not an IUPAC trivial name")
            excluded[letter] = why.strip()

    covered = {a for m in measured.values() for a in m["letters"]}
    unaccounted = [a for a in AA if a not in covered and a not in excluded]
    if unaccounted:
        raise Stop(f"{INPUTS['comparison'].name} accounts for neither a measured row nor "
                   f"[not_measured] for: {', '.join(unaccounted)}. Every one of the twenty must be "
                   "either compared or excluded with a stated reason (X2).")
    both = sorted(covered & set(excluded))
    if both:
        raise Stop(f"{', '.join(both)} is both measured and listed as not measured")
    return measured, excluded, conversions


def check_printed_sums(cp, measured: dict, log) -> list[str]:
    """Rule X4: the transcription must reproduce the paper's own printed sums."""
    loc = cp["values.human_muscle.location"] if "values.human_muscle.location" in cp else {}
    stated = loc.get("sums_as_printed", "")
    if not stated:
        return ["no printed sums given in config; transcription not cross-checked"]
    import re
    nums = [float(x) for x in re.findall(r"(\d+\.\d+)", stated)]
    if len(nums) != 2:
        return [f"could not read two printed sums from {stated!r}; transcription not cross-checked"]
    rows = loc.get("essential_rows", "")
    if not rows.strip():
        raise Stop("X4: [values.human_muscle.location] essential_rows is missing — the paper's own statement of which "
                   "rows its essential sum covers is needed to check the printed sums; code names no amino acid")
    eaa_names = {r.strip().lower() for r in rows.replace(";", ",").split(",") if r.strip()}
    unknown = [r for r in eaa_names if r not in {m["label"].lower() for m in measured.values()}]
    if unknown:
        raise Stop(f"X4: essential_rows names {', '.join(sorted(unknown))}, which is not a row of [values.human_muscle]")
    got_eaa = sum(m["value"] for m in measured.values() if m["label"].lower() in eaa_names)
    got_neaa = sum(m["value"] for m in measured.values() if m["label"].lower() not in eaa_names)
    notes = []
    for what, got, want in (("essential", got_eaa, nums[0]), ("non-essential", got_neaa, nums[1])):
        if abs(got - want) > 0.05:
            raise Stop(f"X4: the transcribed {what} values sum to {got:.1f}, the paper prints {want:.1f}. "
                       "Re-check the transcription against Table 1 before anything is computed.")
        notes.append(f"{what} sum {got:.1f} = printed {want:.1f}")
        log.info("X4 %s sum %.1f matches printed %.1f", what, got, want)
    return notes


def mass_balance(cp, measured: dict, ours: dict[str, float], log) -> tuple[list[list], float]:
    """X6. How much of the expected free-amino-acid yield the published values account for.

    Acid hydrolysis breaks every peptide bond and adds a water to each, so a gram of protein
    yields MORE than a gram of free amino acids. The factor is computed here from this
    repository's own composition and the PubChem free and residue masses — nothing about it comes
    from the paper. Each protein content named in config is then tested: expected yield, the share
    of it the amino acids actually measurable could carry, and what fraction of that the paper's
    printed values reach. No conclusion is drawn from the number (X3).
    """
    mcp = common.read_ini(INPUTS["masses"])
    free = {a: float(mcp[a]["free_mass"]) for a in AA}
    residue = {a: float(mcp[a]["residue_mass"]) for a in AA}
    moles = {a: ours[a] / free[a] for a in AA}
    protein_mass = sum(moles[a] * residue[a] for a in AA)
    free_mass = sum(moles[a] * free[a] for a in AA)
    gain = free_mass / protein_mass                       # free amino acid g per g of protein

    measurable = {a for m in measured.values() for a in m["letters"]}
    share_measurable = sum(ours[a] for a in measurable) / sum(ours.values())
    reported = sum(m["value"] for m in measured.values())

    contents = []
    g_pc = cp["protein_content_of_the_muscle_sample"]
    contents.append(("gorissen_nitrogen_derived", float(g_pc["value"]),
                     f"{INPUTS['comparison'].name} [protein_content_of_the_muscle_sample]: "
                     f"{g_pc.get('as_reported', '').strip()} ({g_pc.get('location', '').strip()})"))
    pcp = common.read_ini(INPUTS["pools"])
    if "protein_content" in pcp and "water_content" in pcp:
        wet_protein, water = float(pcp["protein_content"]["value"]), float(pcp["water_content"]["value"])
        dry = 1000.0 - water
        if dry > 0:
            contents.append(("mingrone_2001_measured", 100.0 * wet_protein / dry,
                             f"{INPUTS['pools'].name} [protein_content] {wet_protein:g} g/kg wet and "
                             f"[water_content] {water:g} g/kg wet, so {100 * wet_protein / dry:.1f} % of dry mass "
                             f"(source id {pcp['protein_content'].get('source_id', '').strip()}, D77)"))

    rows = []
    for name, pct, where in contents:
        expected = pct * gain
        carriable = expected * share_measurable
        rows.append([name, f"{pct:.1f}", f"{gain:.4f}", f"{expected:.1f}", f"{100 * share_measurable:.2f}",
                     f"{carriable:.1f}", f"{reported:.1f}", f"{100 * reported / carriable:.1f}" if carriable else "",
                     where])
        log.info("X6 %s: protein %.1f %%, expected %.1f g, measurable %.1f g, reported %.1f g, recovery %.0f %%",
                 name, pct, expected, carriable, reported, 100 * reported / carriable if carriable else 0)
    return rows, gain


# --------------------------------------------------------------------------- build

def build(log) -> dict[str, str]:
    cp = common.read_ini(INPUTS["comparison"])
    name_to_letter = trivial_to_letter()
    lit = green_light(log)
    cite = citation()

    measured, excluded, conversions = build_sets(cp, name_to_letter)
    sum_notes = check_printed_sums(cp, measured, log)

    ours = {"protein_only": read_standard(INPUTS["standard"], "standard"),
            "with_non_protein_metabolite_pools": read_standard(
                INPUTS["standard_npmp"], "standard_with_non_protein_metabolite_pools")}

    g_total = sum(m["value"] for m in measured.values())
    if g_total <= 0:
        raise Stop("the transcribed values sum to zero")
    folded = {k: {c: sum(ours[c][a] for a in m["letters"]) for c in ours} for k, m in measured.items()}
    our_total = {c: sum(folded[k][c] for k in measured) for c in ours}

    mb_rows, gain = mass_balance(cp, measured, ours["with_non_protein_metabolite_pools"], log)

    hashes = [f"input.{k} = {_rel(p)} sha256 {sha256_path(p)}" for k, p in INPUTS.items() if p.exists()]
    header_common = ([f"generated = {common.iso_now()}"] + hashes +
                     [f"source = {cite}", f"source file = {lit['path']} sha256 {lit['sha256']}"])

    caveats = {k.strip().lower(): v.strip() for k, v in cp["caveats"].items()} if "caveats" in cp else {}

    # ---- the transcription, as letters
    t_rows = []
    for key in sorted(measured, key=lambda k: -measured[k]["value"]):
        m = measured[key]
        t_rows.append([key, m["label"], ", ".join(m["letters"]), f"{m['value']:.1f}",
                       f"{100 * m['value'] / g_total:.4f}",
                       caveats.get(m["label"].lower(), "")])
    t_rows.append(["sum", "", "", f"{g_total:.1f}", "100.0000", ""])
    files = {"gorissen_2018_human_muscle.tsv": tsv_text(
        header_common + [
            "Gorissen et al. 2018, Table 1, column \"Human muscle\" (p. 1690), transcribed by hand and resolved to "
            "one-letter symbols through the IUPAC-IUBMB trivial names; no amino acid is named in code.",
            "g_per_100g_raw_material is the paper's unit: grams per 100 g of freeze-dried tissue, NOT per 100 g of "
            "protein. The same tissue is stated as 84 % protein (p. 1688).",
            "percent_of_measured renormalises over the rows compared (the paper's reported rows less the declared "
            "not-measured zeros, D85), which is the only denominator the two measurements share (X1).",
            f"transcription cross-check (X4): {'; '.join(sum_notes)}.",
        ],
        ["row", "paper_row_label", "amino_acids", "g_per_100g_raw_material", "percent_of_measured", "caveat"],
        t_rows, tool="comparison.py")}

    # ---- the comparison
    c_rows, biggest = [], []
    for key in sorted(measured, key=lambda k: -measured[k]["value"]):
        m = measured[key]
        g_pct = 100 * m["value"] / g_total
        row = [key, m["label"], "compared", f"{g_pct:.4f}"]
        for c in ("protein_only", "with_non_protein_metabolite_pools"):
            o_pct = 100 * folded[key][c] / our_total[c]
            row += [f"{o_pct:.4f}", f"{o_pct - g_pct:+.4f}",
                    f"{o_pct / g_pct:.4f}" if g_pct else ""]
            if c == "with_non_protein_metabolite_pools":
                biggest.append((abs(o_pct - g_pct), key, o_pct - g_pct, g_pct, o_pct))
        row.append(caveats.get(m["label"].lower(), ""))
        c_rows.append(row)
    for letter in sorted(excluded):
        c_rows.append([letter, "", "not measured", "", "", "", "", "", "", "", excluded[letter]])

    cols = ["row", "paper_row_label", "status", "gorissen_percent_of_measured"]
    for c in ("protein_only", "with_non_protein_metabolite_pools"):
        cols += [f"calculated_{c}_percent", f"difference_pp_{c}", f"ratio_calculated_over_gorissen_{c}"]
    cols.append("caveat")
    files["calculated_vs_gorissen_2018.tsv"] = tsv_text(
        header_common + [
            "THE CALCULATED STANDARD BESIDE AN INDEPENDENT MEASUREMENT. Both sides are percentages of the SAME set: "
            "the amino acids Gorissen et al. 2018 report (X1). Neither side is treated as the truth, neither is fitted "
            "or scaled toward the other, and no error term is computed (X3).",
            "The set is decided by the paper's method, not by agreement (X2). Excluded from BOTH sides: "
            + ", ".join(sorted(excluded))
            + " — each one's reason is on its 'not measured' row at the foot of this table. A combined row holds amino "
            "acids the method converts into one another (E+Q: glutamine converts to glutamic acid during hydrolysis), "
            "summed on both sides.",
            f"the comparison set is {our_total['with_non_protein_metabolite_pools']:.2f} % of the calculated standard's "
            f"twenty amino acids; the three excluded carry the rest.",
            "difference_pp is calculated minus Gorissen in percentage points of the comparison set; ratio is "
            "calculated / Gorissen. A row whose Gorissen value is 0.0 has no ratio — see its caveat.",
        ], cols, c_rows, tool="comparison.py")

    # ---- the mass balance (X6)
    files["gorissen_2018_mass_balance.tsv"] = tsv_text(
        header_common + [
            "X6: can the published values be reconciled with a stated protein content? Acid hydrolysis breaks every "
            "peptide bond and adds a water to each, so a gram of protein yields MORE than a gram of free amino acids. "
            f"For this composition the factor is {gain:.3f}, computed from the PubChem free and residue masses and this "
            "repository's own standard — nothing in it comes from the paper.",
            "expected_free_amino_acids_g = protein_content_percent x the factor. measurable_g is the share of that the "
            "amino acids the method can report would carry, on this repository's composition. recovery_percent is the "
            "published total against it.",
            "Two protein contents are tested and NEITHER is chosen (X3). The paper's own is nitrogen x 6.25 on whole "
            "freeze-dried tissue, which counts non-protein nitrogen — creatine, carnosine, free amino acids, "
            "nucleotides, urea. The other is a direct measurement of muscle protein already cited in this repository.",
        ],
        ["protein_content_source", "protein_content_percent", "free_per_protein_factor",
         "expected_free_amino_acids_g", "measurable_share_percent", "measurable_g", "reported_g",
         "recovery_percent", "where_the_protein_content_comes_from"],
        mb_rows, tool="comparison.py")

    # ---- the summary
    s = common.new_ini()
    s.add_section("source")
    s["source"].update({"citation": cite, "file": lit["path"], "sha256": lit["sha256"],
                        "tissue": cp["sample"].get("tissue", ""), "n": cp["sample"].get("n", ""),
                        "protein_content_percent": cp["protein_content_of_the_muscle_sample"].get("value", "")})
    s.add_section("comparison_set")
    s["comparison_set"].update({
        "rows": str(len(measured)),
        "amino_acids_compared": "".join(sorted({a for m in measured.values() for a in m["letters"]})),
        "amino_acids_excluded": "".join(sorted(excluded)),
        "combined_rows": "; ".join(f"{k} (paper row {measured[k]['label']!r})" for k in measured if len(measured[k]["letters"]) > 1) or "none",
        "gorissen_sum_g_per_100g_raw": f"{g_total:.1f}",
        "share_of_calculated_standard_percent": f"{our_total['with_non_protein_metabolite_pools']:.2f}",
        "transcription_check": "; ".join(sum_notes),
    })
    s.add_section("mass_balance")
    for r in mb_rows:
        s["mass_balance"][r[0]] = (f"protein {r[1]} %, expected {r[3]} g, measurable {r[5]} g, "
                                   f"reported {r[6]} g, recovery {r[7]} %")
    s.add_section("largest_differences_with_non_protein_metabolite_pools")
    for i, (_, key, d, g, o) in enumerate(sorted(biggest, reverse=True)[:5], 1):
        s["largest_differences_with_non_protein_metabolite_pools"][f"d{i}"] = (
            f"{key}: Gorissen {g:.2f} %, calculated {o:.2f} %, difference {d:+.2f} pp")
    files["comparison_summary.ini"] = common.render_ini(
        s, ["GENERATED by comparison.py. DO NOT EDIT BY HAND."] + header_common)
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("comparison_check" if args.check else "comparison")
    files = build(log)

    if args.check:
        differ = [n for n, t in files.items()
                  if not (OUT_DIR / n).exists()
                  or strip_generated_line((OUT_DIR / n).read_text(encoding="utf-8")) != strip_generated_line(t)]
        if differ:
            log.error("check: %d file(s) differ from a fresh build: %s", len(differ), differ)
            return 1
        log.info("check: the comparison outputs are current")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (OUT_DIR / name).write_text(text, encoding="utf-8", newline="\n")
        log.info("wrote %s", _rel(OUT_DIR / name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
