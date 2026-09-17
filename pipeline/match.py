#!/usr/bin/env python3
"""
match.py — offline stage: Match Rate, every food against every reference.

The formula is the author's spreadsheet formula, transcribed as it is (D87) and written out by
hand in docs/formula.md, "Match Rate". In one sentence:

    Match Rate = the smallest of (food share / reference share) over the scored amino acids,
                 both shares taken over the scored set only.

The spreadsheet reaches the same number by a longer road — it scales the food to the reference's
essential-amino-acid total, divides each amino acid by the limiting one, sums the "need to
consume" and calls the excess "wasted" — and docs/formula.md shows that every step of that road
cancels back to the minimum ratio. The per-amino-acid ratios are written beside the score so
the surplus each amino acid carries is visible; they are a display, not a second score (M5).

Nothing about the food's protein content enters the score (M1): the score is a property of the
proportion alone, so a food reported per 100 g of food, per 100 g of protein, or per serving
scores the same.

Reads
  config/match_rate.ini                                 hand-written: references, scored sets, food tables, version
  config/fao_2013_indispensable_amino_acids.ini         the nine headings (D62), reduced per D86
  config/gorissen_2018_comparison.ini                   the original reference, and the paper's own essential rows
  data/iupac/amino_acid_symbols.ini                     three-letter / trivial name -> letter: nothing is named in code
  outputs/standard/_calculated_amino_acid_standard*.tsv the calculated references (one column each)
  outputs/usda/amino_acids_per_food.tsv                 the USDA foods (grams per 100 g food, as published)
  config/food_amino_acids_other_sources.csv             the author's hand-maintained foods (M4)

Writes (outputs/match/)
  match_rate_per_food.tsv           one row per food per reference: the score, the limiting amino acid, every ratio
  match_rate_by_reference.tsv       one row per food, the score under each reference side by side (old beside new)
  match_rate_ranked_<reference>.tsv every scored food under that reference, best first
  limiting_amino_acid_counts.tsv    per reference, how many foods each amino acid limits
  foods_not_scored.tsv              food x reference pairs that could not be scored, and why
  match_summary.ini                 references, scored sets as resolved, counts, medians, hashes

Rules (docs/conventions.md, "Match Rate rules")
  M1 The score is the minimum ratio, shares over the scored set; protein content never enters.
  M2 Every reference names its scored set; a scored set is resolved to letters through config, never in code.
  M3 A food missing a scored amino acid is listed, not scored (U5 carried). A published zero is scored as
     published: the ratio is zero, the score is zero, and the limiting amino acid is named.
  M4 Foods come from the USDA table and from the author's other-sources file; no food is dropped or
     preferred across sources (U6 carried), and every row names its source.
  M5 The per-amino-acid ratios are written beside the score. They are the spreadsheet's surplus
     arithmetic, which is the same number, not a second score.
  M6 This is the public single-food scorer. The blend / fortification script is separate (D88) and
     calls match_rate() from here, so the arithmetic lives in one place.

Usage
  python pipeline/match.py            # build
  python pipeline/match.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

from pipeline import common
from pipeline.aggregate import split_list, strip_generated_line, tsv_text
from pipeline.mass_fractions import sha256_path

AA = list(common.AMINO_ACIDS)
OUT_DIR = common.MATCH_OUT_DIR
CONFIG_INI = common.MATCH_RATE_INI


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


# --------------------------------------------------------------------------- the formula (M1)

def match_rate(food: dict[str, float], reference: dict[str, float], scored: list[str]) -> dict:
    """The transcribed formula (D87), on one food against one reference.

    food, reference: {letter: amount} in any consistent unit; only the scored letters are read.
    scored: the letters of the scored set.
    Returns {"score": fraction 0..1, "limiting": [letters], "ratio": {letter: ratio},
             "food_share": {...}, "reference_share": {...}} or {"score": None, "reason": str}.
    """
    missing = [a for a in scored if food.get(a) is None]
    if missing:
        return {"score": None, "reason": f"food does not report {', '.join(missing)}"}
    ref_missing = [a for a in scored if reference.get(a) is None]
    if ref_missing:
        raise Stop(f"the reference does not carry the scored amino acid(s) {', '.join(ref_missing)}")
    f_sum = sum(food[a] for a in scored)
    r_sum = sum(reference[a] for a in scored)
    if r_sum <= 0:
        raise Stop("the reference's scored amino acids sum to zero")
    if f_sum <= 0:
        return {"score": None, "reason": "every scored amino acid is zero"}
    f_share = {a: food[a] / f_sum for a in scored}
    r_share = {a: reference[a] / r_sum for a in scored}
    for a in scored:
        if r_share[a] <= 0:
            raise Stop(f"the reference's share of {a} is zero; the ratio is undefined")
    ratio = {a: f_share[a] / r_share[a] for a in scored}
    k = min(ratio.values())
    limiting = [a for a in scored if ratio[a] == k]
    return {"score": k, "limiting": limiting, "ratio": ratio, "food_share": f_share, "reference_share": r_share}


# --------------------------------------------------------------------------- symbols

def symbol_tables() -> tuple[dict[str, str], dict[str, str]]:
    """({three-letter lowercased: letter}, {trivial name lowercased: letter}) from the IUPAC table."""
    cp = common.read_ini(common.AMINO_ACID_SYMBOLS_INI)
    three, names = {}, {}
    for a in AA:
        if a not in cp:
            raise Stop(f"amino acid {a} missing from {common.AMINO_ACID_SYMBOLS_INI.name}")
        t = cp[a].get("three_letter", "").strip().lower()
        n = cp[a].get("trivial_name", "").strip().lower()
        if not t or not n:
            raise Stop(f"{common.AMINO_ACID_SYMBOLS_INI.name} [{a}] lacks three_letter or trivial_name")
        three[t] = a
        names[n] = a
    return three, names


def _to_letter(token: str, three: dict, names: dict, where: str) -> str:
    t = token.strip().lower()
    if t in three:
        return three[t]
    if t in names:
        return names[t]
    if token.strip().upper() in AA and len(token.strip()) == 1:
        return token.strip().upper()
    raise Stop(f"{where}: {token!r} is neither a three-letter symbol nor a trivial name in the IUPAC table")


# --------------------------------------------------------------------------- scored sets (M2)

def load_scored_sets(cp, three: dict, names: dict) -> dict[str, dict]:
    """{name: {"letters": [...], "how": str, "decision": str}} for every [scored_set.<name>]."""
    out = {}
    for sec in cp.sections():
        if not sec.startswith("scored_set."):
            continue
        name = sec[len("scored_set."):]
        s = cp[sec]
        src = _path(s.get("from_file", ""))
        if not s.get("from_file", "").strip():
            raise Stop(f"[{sec}] has no from_file")
        if not src.exists():
            raise Stop(f"[{sec}] from_file {_rel(src)} not found")
        scp = common.read_ini(src)
        kind = s.get("kind", "").strip()
        letters: list[str] = []
        if kind == "fao_headings_reduced":
            # the nine headings as transcribed (D62); each group scored as the member the author names (D86)
            headings = split_list(scp["indispensable_amino_acids"]["headings"])
            for h in headings:
                gsec = f"group.{h}"
                if gsec in scp:
                    key = f"group_scored_as.{h}"
                    if key not in s:
                        raise Stop(f"[{sec}]: heading {h!r} is a group in {src.name}; say which member is scored "
                                   f"with {key} = <three-letter symbol> (D86)")
                    member = s[key].strip()
                    members = split_list(scp[gsec]["members"])
                    if member.lower() not in [m.lower() for m in members]:
                        raise Stop(f"[{sec}] {key} = {member!r} is not a member of [{gsec}] ({', '.join(members)})")
                    letters.append(_to_letter(member, three, names, f"[{sec}] {key}"))
                else:
                    letters.append(_to_letter(h, three, names, f"[{sec}] heading"))
            how = (f"the {len(headings)} headings of {_rel(src)} [indispensable_amino_acids], each group scored as "
                   + ", ".join(f"{h} -> {s[f'group_scored_as.{h}'].strip()}" for h in headings if f"group.{h}" in scp))
        elif kind == "rows_named_in_file":
            section, key = s.get("from_section", "").strip(), s.get("from_key", "").strip()
            if section not in scp or key not in scp[section]:
                raise Stop(f"[{sec}]: {_rel(src)} has no [{section}] {key}")
            rows = split_list(scp[section][key])
            letters = [_to_letter(r, three, names, f"[{sec}] {key}") for r in rows]
            how = f"the rows {_rel(src)} [{section}] {key} names: {', '.join(rows)}"
        else:
            raise Stop(f"[{sec}] kind must be fao_headings_reduced or rows_named_in_file, not {kind!r}")
        if len(set(letters)) != len(letters):
            raise Stop(f"[{sec}] names an amino acid twice: {letters}")
        out[name] = {"letters": letters, "how": how, "decision": s.get("decision", "").strip(),
                     "note": s.get("note", "").strip()}
    if not out:
        raise Stop(f"{CONFIG_INI.name} defines no [scored_set.*]")
    return out


# --------------------------------------------------------------------------- references (M2)

def read_standard_column(path: Path, column: str) -> dict[str, float]:
    if not path.exists():
        raise Stop(f"{_rel(path)} not found — run `python run.py standard` first")
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    cols = lines[0].split("\t")
    if column not in cols or "amino_acid" not in cols:
        raise Stop(f"{_rel(path)} has no column {column!r} (found {cols})")
    i_aa, i_v = cols.index("amino_acid"), cols.index(column)
    out = {}
    for ln in lines[1:]:
        parts = ln.split("\t")
        if parts[i_aa] in AA and parts[i_v].strip():
            out[parts[i_aa]] = float(parts[i_v])
    missing = [a for a in AA if a not in out]
    if missing:
        raise Stop(f"{_rel(path)} is missing {', '.join(missing)}")
    return out


def read_config_reference(path: Path, section: str, three: dict, names: dict) -> dict[str, float]:
    """A reference transcribed in a config section: {row label: value}, labels resolved through the IUPAC table."""
    if not path.exists():
        raise Stop(f"{_rel(path)} not found")
    cp = common.read_ini(path)
    if section not in cp:
        raise Stop(f"{_rel(path)} has no [{section}]")
    out = {}
    for label, raw in cp[section].items():
        if raw.strip() in ("", "___"):
            raise Stop(f"{_rel(path)} [{section}] {label} is not filled in")
        out[_to_letter(label, three, names, f"{_rel(path)} [{section}]")] = float(raw)
    return out


def load_references(cp, scored_sets: dict, three: dict, names: dict) -> list[dict]:
    refs = []
    for sec in cp.sections():
        if not sec.startswith("reference."):
            continue
        name = sec[len("reference."):]
        s = cp[sec]
        ss = s.get("scored_set", "").strip()
        if ss not in scored_sets:
            raise Stop(f"[{sec}] scored_set {ss!r} is not a [scored_set.*] in {CONFIG_INI.name}")
        if s.get("file", "").strip():
            path = _path(s["file"])
            column = s.get("column", "").strip()
            if not column:
                raise Stop(f"[{sec}] names a file but no column")
            values = read_standard_column(path, column)
            where = f"{_rel(path)} column {column}"
        elif s.get("config", "").strip():
            path = _path(s["config"])
            section = s.get("section", "").strip()
            values = read_config_reference(path, section, three, names)
            where = f"{_rel(path)} [{section}]"
        else:
            raise Stop(f"[{sec}] names neither file+column nor config+section")
        label = s.get("label", "").strip()
        if not label:
            raise Stop(f"[{sec}] has no label")
        refs.append({"name": name, "label": label, "scored_set": ss, "letters": scored_sets[ss]["letters"],
                     "values": values, "where": where, "path": path, "decisions": s.get("decisions", "").strip(),
                     "note": s.get("note", "").strip()})
    if not refs:
        raise Stop(f"{CONFIG_INI.name} defines no [reference.*]")
    return refs


# --------------------------------------------------------------------------- foods (M4)

def _num(x: str) -> float | None:
    x = (x or "").strip()
    return None if x == "" else float(x)


def read_usda_foods(path: Path) -> list[dict]:
    if not path.exists():
        raise Stop(f"{_rel(path)} not found — run `python run.py usda` first")
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    cols = lines[0].split("\t")
    need = ["fdc_id", "data_type", "release", "description", "protein_g_per_100g"] + [f"{a}_g_per_100g" for a in AA]
    absent = [c for c in need if c not in cols]
    if absent:
        raise Stop(f"{_rel(path)} lacks the column(s) {', '.join(absent)}")
    foods = []
    for ln in lines[1:]:
        r = dict(zip(cols, ln.split("\t")))
        foods.append({
            "source": "usda", "source_detail": f"{r['data_type']} {r['release']}", "food_id": r["fdc_id"],
            "description": r["description"], "category": r.get("food_category", ""),
            "protein": _num(r["protein_g_per_100g"]), "min_data_points": r.get("min_data_points", ""),
            "values": {a: _num(r[f"{a}_g_per_100g"]) for a in AA},
            "citation": "",
        })
    return foods


OTHER_REQUIRED = ["label", "source", "location", "note", "protein_g_per_100g"] + [f"{a}_g_per_100g" for a in AA]


def read_other_foods(path: Path) -> list[dict]:
    """The author's hand-maintained foods (M4). A header-only file is an empty list, not an error."""
    if not path.exists():
        raise Stop(f"{_rel(path)} not found (the hand-maintained other-sources file; it may be header-only)")
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        absent = [c for c in OTHER_REQUIRED if c not in cols]
        if absent:
            raise Stop(f"{_rel(path)} lacks the column(s) {', '.join(absent)}")
        foods, seen = [], set()
        for i, r in enumerate(reader, 2):
            label = (r.get("label") or "").strip()
            if not label:
                continue
            if label in seen:
                raise Stop(f"{_rel(path)} line {i}: label {label!r} appears twice; labels are the food's identity here")
            seen.add(label)
            if not (r.get("source") or "").strip():
                raise Stop(f"{_rel(path)} line {i} ({label}): source is empty; every hand-entered food cites where its values come from")
            try:
                values = {a: _num(r[f"{a}_g_per_100g"]) for a in AA}
                protein = _num(r["protein_g_per_100g"])
            except ValueError as e:
                raise Stop(f"{_rel(path)} line {i} ({label}): {e}")
            foods.append({"source": "other", "source_detail": (r.get("source") or "").strip(), "food_id": label,
                          "description": label, "category": "", "protein": protein, "min_data_points": "",
                          "values": values,
                          "citation": "; ".join(x for x in ((r.get("source") or "").strip(),
                                                           (r.get("location") or "").strip()) if x)})
    return foods


# --------------------------------------------------------------------------- build

def _pct(x: float, d: int) -> str:
    return f"{100 * x:.{d}f}"


def build(log) -> dict[str, str]:
    if not CONFIG_INI.exists():
        raise Stop(f"{_rel(CONFIG_INI)} not found")
    cp = common.read_ini(CONFIG_INI)
    three, names = symbol_tables()
    scored_sets = load_scored_sets(cp, three, names)
    refs = load_references(cp, scored_sets, three, names)
    primary = cp["formula"].get("primary_reference", "").strip() if "formula" in cp else ""
    if primary and primary not in [r["name"] for r in refs]:
        raise Stop(f"[formula] primary_reference {primary!r} is not a [reference.*]")
    version = cp["meta"].get("version", "").strip() if "meta" in cp else ""
    decimals = int(cp["output"].get("decimals", "4")) if "output" in cp else 4

    food_paths = {}
    for sec in cp.sections():
        if sec.startswith("foods."):
            food_paths[sec[len("foods."):]] = _path(cp[sec]["file"])
    if "usda" not in food_paths or "other_sources" not in food_paths:
        raise Stop(f"{CONFIG_INI.name} must name [foods.usda] and [foods.other_sources]")
    foods = read_usda_foods(food_paths["usda"]) + read_other_foods(food_paths["other_sources"])
    log.info("%d foods (%d USDA, %d other sources); %d references",
             len(foods), sum(f["source"] == "usda" for f in foods), sum(f["source"] == "other" for f in foods), len(refs))

    # ---- score
    all_letters = sorted({a for r in refs for a in r["letters"]}, key=AA.index)
    per_food_rows, not_scored, scores = [], [], {}       # scores[(source, food_id)][ref] = result
    for f in foods:
        key = (f["source"], f["food_id"])
        scores[key] = {}
        for r in refs:
            res = match_rate(f["values"], r["values"], r["letters"])
            scores[key][r["name"]] = res
            if res["score"] is None:
                not_scored.append([f["source"], f["source_detail"], f["food_id"], f["description"], r["name"], res["reason"]])
                continue
            scored_sum = sum(f["values"][a] for a in r["letters"])
            row = [f["source"], f["source_detail"], f["food_id"], f["description"], f["category"],
                   r["name"], r["label"], r["scored_set"], "".join(r["letters"]),
                   _pct(res["score"], decimals), "+".join(res["limiting"]),
                   f"{scored_sum:.{decimals}f}", "" if f["protein"] is None else f"{f['protein']:g}", f["min_data_points"]]
            for a in all_letters:
                row.append(f"{res['ratio'][a]:.{decimals}f}" if a in res["ratio"] else "")
            row.append(f["citation"])
            per_food_rows.append(row)

    # ---- headers
    inputs = {"config": CONFIG_INI, "symbols": common.AMINO_ACID_SYMBOLS_INI,
              "usda_foods": food_paths["usda"], "other_sources": food_paths["other_sources"]}
    for r in refs:
        inputs[f"reference.{r['name']}"] = r["path"]
    for ss in [s for s in cp.sections() if s.startswith("scored_set.")]:
        inputs[ss] = _path(cp[ss]["from_file"])
    seen_hash, hash_lines = set(), []
    for k, p in inputs.items():
        if p.exists() and p not in seen_hash:
            seen_hash.add(p)
            hash_lines.append(f"input.{k} = {_rel(p)} sha256 {sha256_path(p)}")
    header_common = ([f"generated = {common.iso_now()}"] + hash_lines +
                     [f"Match Rate version = {version} (config/match_rate.ini [meta] version)" if version else "Match Rate version = (unlabelled)"] +
                     [f"reference.{r['name']} = {r['label']} | {r['where']} | scored set {r['scored_set']} = {''.join(r['letters'])} | {r['decisions']}"
                      for r in refs] +
                     [f"primary reference = {primary}" if primary else "primary reference = (none named)"] +
                     ["M1 (D87): match_rate_percent = 100 x min over the scored set of (food share / reference share), both shares "
                      "taken over the scored set only. Protein content never enters: the score is a property of the proportion. "
                      "docs/formula.md 'Match Rate' shows the spreadsheet's surplus arithmetic reduces to this number.",
                      "M3: a food that does not report a scored amino acid is in foods_not_scored.tsv, not here; a published zero "
                      "scores zero with the limiting amino acid named. M4: usda rows are outputs/usda/amino_acids_per_food.tsv as "
                      "published; other rows are the author's config/food_amino_acids_other_sources.csv, each citing its source. "
                      "A food in two archives, or in USDA and the other-sources file, is two rows (U6)."])

    files: dict[str, str] = {}
    cols = (["source", "source_detail", "food_id", "description", "food_category",
             "reference", "reference_label", "scored_set", "scored_amino_acids",
             "match_rate_percent", "limiting_amino_acid", "scored_sum_g_per_100g", "protein_g_per_100g", "min_data_points"]
            + [f"ratio_{a}" for a in all_letters] + ["citation"])
    files["match_rate_per_food.tsv"] = tsv_text(
        header_common + ["ratio_<letter> = food share / reference share for that amino acid (M5): 1.0000 is the limiting "
                         "amino acid, anything above it is surplus relative to the reference; blank where the amino acid is not "
                         "in that reference's scored set."],
        cols, per_food_rows, tool="match.py")

    # ---- wide: one row per food, every reference side by side
    ref_names = [r["name"] for r in refs]
    wide_rows = []
    for f in foods:
        key = (f["source"], f["food_id"])
        row = [f["source"], f["source_detail"], f["food_id"], f["description"], "" if f["protein"] is None else f"{f['protein']:g}"]
        for rn in ref_names:
            res = scores[key][rn]
            row += [_pct(res["score"], decimals) if res["score"] is not None else "",
                    "+".join(res["limiting"]) if res["score"] is not None else ""]
        if primary and len(ref_names) > 1:
            p = scores[key][primary]
            for rn in ref_names:
                if rn == primary:
                    continue
                o = scores[key][rn]
                row.append(_pct(p["score"] - o["score"], decimals) if (p["score"] is not None and o["score"] is not None) else "")
        wide_rows.append(row)
    wcols = ["source", "source_detail", "food_id", "description", "protein_g_per_100g"]
    for rn in ref_names:
        wcols += [f"match_{rn}_percent", f"limiting_{rn}"]
    if primary and len(ref_names) > 1:
        wcols += [f"difference_pp_{primary}_minus_{rn}" for rn in ref_names if rn != primary]
    files["match_rate_by_reference.tsv"] = tsv_text(
        header_common + ["One row per food; the score and limiting amino acid under every reference, side by side. difference_pp "
                         "is the primary reference's score minus the other's, in percentage points; blank where either is unscored. "
                         "Two references with different scored sets (eight vs nine amino acids) score different questions: the "
                         "difference is reported, not explained (X3 pattern)."],
        wcols, wide_rows, tool="match.py")

    # ---- ranked, per reference
    for r in refs:
        rows = []
        for f in foods:
            res = scores[(f["source"], f["food_id"])][r["name"]]
            if res["score"] is None:
                continue
            rows.append((res["score"], f, res))
        rows.sort(key=lambda t: (-t[0], t[1]["description"].lower(), t[1]["source"], t[1]["food_id"]))
        files[f"match_rate_ranked_{r['name']}.tsv"] = tsv_text(
            header_common + [f"Every scored food under reference {r['name']} ({r['label']}), best first; ties by description. "
                             f"Scored set {r['scored_set']} = {''.join(r['letters'])}."],
            ["rank", "source", "source_detail", "food_id", "description", "match_rate_percent", "limiting_amino_acid", "protein_g_per_100g"],
            [[i, f["source"], f["source_detail"], f["food_id"], f["description"], _pct(s, decimals), "+".join(res["limiting"]),
              "" if f["protein"] is None else f"{f['protein']:g}"] for i, (s, f, res) in enumerate(rows, 1)],
            tool="match.py")

    # ---- limiting counts
    lc_rows = []
    for r in refs:
        counts = {a: 0 for a in r["letters"]}
        n = 0
        for f in foods:
            res = scores[(f["source"], f["food_id"])][r["name"]]
            if res["score"] is None:
                continue
            n += 1
            for a in res["limiting"]:
                counts[a] += 1
        for a in r["letters"]:
            lc_rows.append([r["name"], a, counts[a], n, f"{100 * counts[a] / n:.2f}" if n else ""])
    files["limiting_amino_acid_counts.tsv"] = tsv_text(
        header_common + ["Per reference: how many scored foods each amino acid limits (a tie counts for each tied amino acid). "
                         "This is the table that answers 'where do scores move between references'."],
        ["reference", "amino_acid", "foods_limited", "foods_scored", "percent_of_scored"], lc_rows, tool="match.py")

    files["foods_not_scored.tsv"] = tsv_text(
        header_common + ["Food x reference pairs the formula could not score, with the reason (M3). Nothing here is a judgement "
                         "of the food; it is the record of what the source did not report."],
        ["source", "source_detail", "food_id", "description", "reference", "reason"], not_scored, tool="match.py")

    # ---- summary
    s = common.new_ini()
    s.add_section("run")
    s["run"]["match_rate_version"] = version
    s["run"]["primary_reference"] = primary
    s["run"]["foods_read"] = str(len(foods))
    s["run"]["foods_usda"] = str(sum(f["source"] == "usda" for f in foods))
    s["run"]["foods_other_sources"] = str(sum(f["source"] == "other" for f in foods))
    s["run"]["references"] = ", ".join(ref_names)
    for ss_name, ss in scored_sets.items():
        s.add_section(f"scored_set.{ss_name}")
        s[f"scored_set.{ss_name}"]["amino_acids"] = "".join(ss["letters"])
        s[f"scored_set.{ss_name}"]["how"] = ss["how"]
        s[f"scored_set.{ss_name}"]["decision"] = ss["decision"]
    for r in refs:
        sec = f"reference.{r['name']}"
        s.add_section(sec)
        vals = [scores[(f["source"], f["food_id"])][r["name"]]["score"] for f in foods]
        vals = [v for v in vals if v is not None]
        s[sec]["label"] = r["label"]
        s[sec]["read_from"] = r["where"]
        s[sec]["scored_set"] = f"{r['scored_set']} = {''.join(r['letters'])}"
        s[sec]["reference_shares_percent"] = ", ".join(
            f"{a} {100 * r['values'][a] / sum(r['values'][b] for b in r['letters']):.{decimals}f}" for a in r["letters"])
        s[sec]["foods_scored"] = str(len(vals))
        s[sec]["foods_not_scored"] = str(len(foods) - len(vals))
        s[sec]["median_match_rate_percent"] = _pct(statistics.median(vals), decimals) if vals else ""
        s[sec]["foods_at_zero"] = str(sum(1 for v in vals if v == 0))
    files["match_summary.ini"] = common.render_ini(s, ["GENERATED by match.py. DO NOT EDIT BY HAND."] + header_common)
    return files


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("match_check" if args.check else "match")
    files = build(log)

    if args.check:
        differ = [n for n, t in files.items()
                  if not (OUT_DIR / n).exists()
                  or strip_generated_line((OUT_DIR / n).read_text(encoding="utf-8")) != strip_generated_line(t)]
        if differ:
            log.error("check: %d file(s) differ from a fresh build: %s", len(differ), differ)
            return 1
        log.info("check: the match outputs are current")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (OUT_DIR / name).write_text(text, encoding="utf-8", newline="\n")
        log.info("wrote %s", _rel(OUT_DIR / name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
