"""Disclose the mass of post-translational modifications not modelled by composition.py.

PTM = post-translational modification: a chemical change made to a protein after it
has been synthesised from its sequence -- a phosphate, methyl, acetyl, sugar, or lipid
group attached to a residue, or a disulfide bond formed between two cysteines. Cells
use them to regulate a protein (switch its activity on or off, change where it goes,
how long it lasts, what it binds). For this repository only the mechanical fact
matters: the atoms a PTM adds or removes are not in the sequence letters.

composition.py counts sequence letters. A post-translational modification (PTM)
adds or removes atoms after the chain is made; the letter stays the same, so the
engine counts the unmodified residue. This stage decides nothing about that. It
computes how much mass the annotated PTMs of each entry would add to or remove
from the mature chain, so the methods section can state the size of what is not
modelled as a number.

Per entry it reads the raw UniProt JSON (data/uniprot_raw/<accession>.json) and
takes every feature of the PTM-class types

    Modified residue, Lipidation, Glycosylation, Cross-link, Disulfide bond

For each feature:
  - the description's first clause (the text before the first ';') is looked up
    as an ID in UniProt's PTM vocabulary (data/uniprot_ptmlist/ptmlist.txt). The
    lookup is exact first; if that fails, vocabulary IDs that contain "..." are
    treated as the wildcard patterns they are in the file (e.g. the vocabulary ID
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-...)" covers the feature
    description "... (interchain with G-Cter in ubiquitin)"). A leading
    "(Microbial infection) " qualifier on the description is set aside for the
    lookup and recorded in its own column. The rule that matched is recorded per site;
  - the vocabulary's correction formula (CF) is recorded per site, so the reader can
    see what the mass delta represents -- for an interchain cross-link it is the
    linkage (e.g. loss of one water), not the partner molecule;
  - if the record has an MA line, that average mass difference is the site's mass
    delta (mass_status = vocabulary);
  - a Disulfide bond has no vocabulary record: within one chain it removes two
    hydrogens, so its delta is minus the molecular weight of H2 from
    data/pubchem/amino_acid_masses.ini (mass_status = disulfide_arithmetic);
    a bond whose description contains "Interchain" is counted only;
  - otherwise the site is counted but carries no mass (mass_status = no_mass_in_vocabulary
    or no_vocabulary_match), and the description is listed in ptm_summary.ini so the
    reader can see exactly what was not priced;
  - a site whose position lies outside the master molecule (config/segments.ini)
    is recorded with in_master = false and excluded from the sums;
  - C3b (D118): priced features of the same type at the same position are ALTERNATIVES
    (UniProt lists the glycoforms seen at one glycosylation site as separate features),
    not additions: one position contributes the largest |delta| among them to the entry's
    sum, and the number of features so collapsed is recorded per entry.

Outputs  outputs/composition/ptm_sites.tsv        one row per PTM feature
         outputs/composition/ptm_mass_deltas.tsv  one row per entry: counts, summed delta, fraction of MW
         outputs/composition/ptm_summary.ini      set-wide maxima, totals, unpriced descriptions

Inputs   data/uniprot_raw/*.json, data/uniprot_ptmlist/ptmlist.txt,
         data/pubchem/amino_acid_masses.ini, config/accessions.ini, config/segments.ini,
         outputs/composition/amino_acid_composition_per_protein.tsv (for mw of the master chain)
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from pipeline import common
from pipeline.composition import segment_ranges
from pipeline.fetch_ptmlist import parse_ptmlist

STAGE = "ptm_disclosure"
PTM_FEATURE_TYPES = ("Modified residue", "Lipidation", "Glycosylation", "Cross-link", "Disulfide bond")
FRAC_DIGITS = 10


QUALIFIER_PREFIX = "(Microbial infection) "


def lookup_key(description: str) -> tuple[str, str]:
    """(key, qualifier): the description's first clause, with a leading qualifier set aside."""
    key = description.split(";", 1)[0].strip()
    qualifier = ""
    if key.startswith(QUALIFIER_PREFIX):
        qualifier = QUALIFIER_PREFIX.strip()
        key = key[len(QUALIFIER_PREFIX):].strip()
    return key, qualifier


def wildcard_patterns(vocab: dict) -> list[tuple[re.Pattern, str]]:
    """Vocabulary IDs containing '...' as compiled full-match patterns, in file order."""
    out = []
    for vid in vocab:
        if "..." in vid:
            out.append((re.compile(".*".join(re.escape(part) for part in vid.split("..."))), vid))
    return out


def find_record(key: str, vocab: dict, patterns: list[tuple[re.Pattern, str]]) -> tuple[dict | None, str, str]:
    """Return (record, vocabulary_id, match_rule): exact match first, then wildcard."""
    if key in vocab:
        return vocab[key], key, "exact"
    hits = [(vid, pat) for pat, vid in patterns if pat.fullmatch(key)]
    if len(hits) == 1:
        return vocab[hits[0][0]], hits[0][0], "wildcard"
    if len(hits) > 1:
        return None, ";".join(v for v, _ in hits), "ambiguous_wildcard"
    return None, "", "none"


def _loc(feature: dict) -> tuple[int | None, int | None, str, str]:
    loc = feature.get("location", {})
    s, e = loc.get("start", {}), loc.get("end", {})
    return s.get("value"), e.get("value"), s.get("modifier", ""), e.get("modifier", "")


def in_ranges(pos: int | None, ranges: list[tuple[int, int]]) -> bool:
    return pos is not None and any(a <= pos <= b for a, b in ranges)


def price_feature(feature: dict, vocab: dict, h2_mass: float, master_ranges: list[tuple[int, int]],
                  patterns: list[tuple[re.Pattern, str]] | None = None) -> dict:
    """Classify one PTM feature and assign its average mass delta, if determinable."""
    if patterns is None:
        patterns = wildcard_patterns(vocab)
    ftype = feature.get("type", "")
    desc = feature.get("description", "")
    start, end, smod, emod = _loc(feature)
    key, qualifier = lookup_key(desc)
    row = {
        "feature_type": ftype, "start": start, "end": end,
        "start_modifier": smod, "end_modifier": emod,
        "description": desc, "qualifier": qualifier, "lookup_key": key,
        "vocabulary_id": "", "vocabulary_ac": "", "match_rule": "", "correction_formula": "",
        "mass_delta": None, "mass_status": "",
        "in_master": in_ranges(start, master_ranges) and in_ranges(end, master_ranges),
    }
    if ftype == "Disulfide bond":
        row["match_rule"] = "not_in_vocabulary"
        if "interchain" in desc.lower():
            row["mass_status"] = "interchain_count_only"
        else:
            row["mass_delta"] = -h2_mass
            row["mass_status"] = "disulfide_arithmetic"
            row["correction_formula"] = "H-2 (arithmetic)"
        return row
    rec, vid, rule = find_record(key, vocab, patterns)
    row["match_rule"], row["vocabulary_id"] = rule, vid
    if rec is None:
        row["mass_status"] = "no_vocabulary_match"
        return row
    row["vocabulary_ac"] = ";".join(rec.get("AC", []))
    row["correction_formula"] = ";".join(rec.get("CF", []))
    ma = rec.get("MA", [])
    if not ma:
        row["mass_status"] = "no_mass_in_vocabulary"
        return row
    row["mass_delta"] = float(ma[0])
    row["mass_status"] = "vocabulary"
    return row


def collapse_alternatives(priced: list[tuple[str, int | None, int | None, float]]) -> tuple[float, int]:
    """C3b: (the entry's summed delta with one contribution per (type, start, end) -- the largest
    |delta| among alternatives -- and the number of features collapsed away)."""
    best: dict[tuple, float] = {}
    for ftype, start, end, delta in priced:
        key = (ftype, start, end)
        if key not in best or abs(delta) > abs(best[key]):
            best[key] = delta
    return sum(best.values()), len(priced) - len(best)


def _f(x: float) -> str:
    return f"{x:.{FRAC_DIGITS}g}"


def master_mw_by_accession(all_tsv: Path) -> dict[str, float]:
    out = {}
    for r in common.read_tsv(all_tsv):
        if r["segment_set"] == "master":
            out[r["accession"]] = float(r["mw"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", default=None, help="folder of <accession>.json; default: data/uniprot_raw/ and every category subfolder in it (D110)")
    ap.add_argument("--category", default=None, help="a non-muscle category: its accessions and segments under config/<category>/, tables under outputs/<category>/intermediate/composition/ (Step 7b)")
    args = ap.parse_args(argv)
    log = common.make_logger(STAGE if args.category is None else f"{STAGE}_{args.category}")
    files = common.category_files(args.category)
    out_dir = files["composition_dir"]

    raw_dir = Path(args.raw_dir) if args.raw_dir else common.UNIPROT_RAW_DIR
    if not raw_dir.exists():
        log.error("%s not found: this stage needs the raw entry JSON saved by fetch_sequences.py", raw_dir)
        return 1

    def raw_entry(acc: str) -> Path | None:
        if args.raw_dir:
            p = raw_dir / f"{acc}.json"
            return p if p.exists() else None
        return common.uniprot_raw_entry(acc)
    accessions = common.read_ini(files["accessions"])
    segments = common.read_ini(files["segments"])
    vocab = parse_ptmlist(common.PTMLIST_TXT.read_text(encoding="utf-8", errors="replace"))
    patterns = wildcard_patterns(vocab)
    masses = common.read_ini(common.AMINO_ACID_MASSES_INI)
    h2_mass = float(masses["hydrogen"]["molecular_weight"])
    mw_master = master_mw_by_accession(files["composition_tsv"])
    ptm_src = common.read_ini(common.PTMLIST_TXT.parent / "source.ini")["ptmlist"]

    site_rows: list[list] = []
    entry_rows: list[list] = []
    unpriced: dict[str, int] = {}
    status_totals: dict[str, int] = {}
    rule_totals: dict[str, int] = {}
    best = (0.0, "")
    n_sites_total = 0
    n_collapsed_total = 0
    for acc in accessions.sections():
        entry = accessions[acc]
        if entry.get("flag_open") == "true":
            continue
        path = raw_entry(acc)
        if path is None:
            log.error("%s: raw JSON missing under %s (root or any category subfolder)", acc, raw_dir)
            return 1
        data = json.loads(path.read_text(encoding="utf-8"))
        gene, tier = entry.get("gene", ""), entry.get("tier", "")
        master_ranges = segment_ranges(segments, acc, master_only=True)
        by_type = {t: 0 for t in PTM_FEATURE_TYPES}
        n_in_master = n_with_mass = n_count_only = 0
        priced: list[tuple[str, int | None, int | None, float]] = []
        for f in data.get("features", []):
            if f.get("type") not in PTM_FEATURE_TYPES:
                continue
            p = price_feature(f, vocab, h2_mass, master_ranges, patterns)
            n_sites_total += 1
            by_type[p["feature_type"]] += 1
            status_totals[p["mass_status"]] = status_totals.get(p["mass_status"], 0) + 1
            rule_totals[p["match_rule"]] = rule_totals.get(p["match_rule"], 0) + 1
            site_rows.append([acc, gene, tier, p["feature_type"], p["start"], p["end"],
                              p["start_modifier"], p["end_modifier"], p["description"], p["qualifier"],
                              p["lookup_key"], p["vocabulary_id"], p["vocabulary_ac"], p["match_rule"],
                              p["correction_formula"], "" if p["mass_delta"] is None else _f(p["mass_delta"]),
                              p["mass_status"], str(p["in_master"]).lower()])
            if not p["in_master"]:
                continue
            n_in_master += 1
            if p["mass_delta"] is None:
                n_count_only += 1
                label = f"{p['feature_type']}: {p['lookup_key']} [{p['mass_status']}]"
                unpriced[label] = unpriced.get(label, 0) + 1
            else:
                n_with_mass += 1
                priced.append((p["feature_type"], p["start"], p["end"], p["mass_delta"]))
        delta_sum, n_collapsed = collapse_alternatives(priced)                     # C3b (D118)
        n_collapsed_total += n_collapsed
        mw = mw_master.get(acc)
        if mw is None:
            log.error("%s: no master row in %s -- run composition first", acc, files["composition_tsv"].name)
            return 1
        frac = delta_sum / mw if mw else 0.0
        if abs(frac) > best[0]:
            best = (abs(frac), acc)
        entry_rows.append([acc, gene, tier, sum(by_type.values()),
                           n_in_master, n_with_mass, n_count_only, *[by_type[t] for t in PTM_FEATURE_TYPES],
                           _f(delta_sum), f"{mw:.4f}", _f(frac), n_collapsed])

    out_dir.mkdir(parents=True, exist_ok=True)
    common.write_tsv(out_dir / "ptm_sites.tsv",
                     ["accession", "gene", "tier", "feature_type", "start", "end", "start_modifier", "end_modifier",
                      "description", "qualifier", "lookup_key", "vocabulary_id", "vocabulary_ac", "match_rule",
                      "correction_formula", "average_mass_delta_g_per_mol", "mass_status", "in_master"],
                     site_rows)
    type_cols = [f"n_{t.lower().replace(' ', '_').replace('-', '_')}" for t in PTM_FEATURE_TYPES]
    common.write_tsv(out_dir / "ptm_mass_deltas.tsv",
                     ["accession", "gene", "tier", "n_sites_total", "n_sites_in_master", "n_sites_with_mass",
                      "n_sites_count_only", *type_cols, "sum_mass_delta_g_per_mol", "mw_master",
                      "ptm_mass_delta_fraction_of_mw", "n_alternative_features_collapsed"],
                     entry_rows)
    summary = common.new_ini()
    summary["ptm_disclosure"] = {
        "ptmlist_url": ptm_src.get("url", ""),
        "ptmlist_sha256": ptm_src.get("sha256", ""),
        "ptmlist_retrieved": ptm_src.get("retrieved", ""),
        "feature_types": "; ".join(PTM_FEATURE_TYPES),
        "entries": str(len(entry_rows)),
        "entries_with_any_ptm_site": str(sum(1 for r in entry_rows if r[3] > 0)),
        "sites_total": str(n_sites_total),
        "sites_by_mass_status": "; ".join(f"{k}={v}" for k, v in sorted(status_totals.items())),
        "sites_by_match_rule": "; ".join(f"{k}={v}" for k, v in sorted(rule_totals.items())),
        "alternative_features_collapsed_c3b": str(n_collapsed_total),
        "largest_abs_ptm_mass_fraction_of_mw": _f(best[0]),
        "largest_abs_ptm_mass_fraction_accession": best[1],
        "h2_mass_used_for_disulfide": str(h2_mass),
    }
    summary["unpriced_sites_in_master"] = {f"site.{i}": f"{k} x{v}" for i, (k, v) in enumerate(sorted(unpriced.items()))} \
        or {"none": "every in-master site received a mass delta"}
    common.write_ini(summary, out_dir / "ptm_summary.ini",
                     ["ptm_summary.ini -- GENERATED by ptm_disclosure.py. DO NOT EDIT BY HAND.",
                      "Mass of annotated post-translational modifications not modelled by composition.py."])
    log.info("%d PTM sites across %d entries; by status: %s; largest |delta|/MW = %s (%s); %d unpriced descriptions",
             n_sites_total, len(entry_rows), summary["ptm_disclosure"]["sites_by_mass_status"],
             _f(best[0]), best[1], len(unpriced))
    return 0


if __name__ == "__main__":
    sys.exit(main())
