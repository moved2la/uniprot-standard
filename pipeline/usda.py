#!/usr/bin/env python3
"""
usda.py — offline stage: the USDA FoodData Central food tables, as one row per food.

Match Rate scores a food's amino acid profile against the standard. The food side comes from
USDA FoodData Central, which publishes its data as a long table — one row per food per nutrient,
per 100 g of food, with vitamins, minerals, macronutrients and amino acids all in the same file,
distinguished only by a nutrient id. This stage pivots that long table to one row per food with a
column per amino acid, and nothing else: it makes no scoring decision, drops no food for being
incomplete, and converts no value.

The archives are downloaded by hand from the FoodData Central download page (D80 — fetch by code
is retired; a static file is obtained once and hashed) and dropped unrenamed into data/usda/.
They are gitignored; data/usda/manifest.ini is the record and is committed.

Reads
  config/usda_food_data.ini                which archives, which release, which extra nutrients (hand-written)
  data/usda/*.zip                          the archives as downloaded
  data/iupac/amino_acid_symbols.ini        trivial name -> one-letter symbol: the join key for the twenty

Writes (outputs/usda/, every file with a header naming its inputs and their hashes)
  amino_acids_per_food.tsv         one row per food that carries at least one amino acid: fdc_id, data type,
                                   release, description, food category, protein and nitrogen where given, then
                                   <letter>_g_per_100g for each of the twenty the archive reports, the extra
                                   carried nutrients, the count present and the letters absent
  foods_without_amino_acids.tsv    the archive's own foods that carry no amino acid row at all, so that what is
                                   NOT in the table above is on the record too
  usda_nutrient_map.tsv            the audit trail of rule U3: every nutrient id carried, the name the archive
                                   gives it, its unit, and what it was carried as
  usda_archive_summary.ini         per archive: file, hash, size, counts; and the run's flags
  data/usda/manifest.ini           per archive: file name, SHA-256, bytes, data type, release, first seen
  logs/usda_<UTC>.log

Rules (docs/conventions.md, "USDA rules")
  U1 Archives. Only the data types config/usda_food_data.ini names. A zip is identified by reading its own
     food.csv, never by its file name; the file may be named anything. The SHA-256 is the archive's identity
     and goes in the manifest and in every output header. Two zips of the same data type is a [STOP].
  U2 Membership. A food counts when it is listed in the archive's own membership file <data_type>.csv. That
     file is the archive's statement of which rows are current: FoodData Central ships superseded earlier
     versions of a food in food.csv, and they are not in the membership file. If the membership file is
     absent, the data_type column of food.csv is the fallback and the log says which rule was used.
  U3 The name join. A nutrient is amino acid X when its name in the archive's nutrient.csv equals, ignoring
     case and surrounding space, the IUPAC-IUBMB trivial name of X in data/iupac/amino_acid_symbols.ini. No
     amino acid is named in code or in config. Nutrients carried that are not one of the twenty are named in
     config by the name the archive uses, and carried under their own column.
  U4 Values. The archive's `amount` is the published value, per 100 g of food, as reported; min, max, median
     and data_points are carried where the archive has them. Nothing is converted, scaled or filled in. A
     carried nutrient whose unit is not the one config expects is a flag and its values are not carried.
  U5 Completeness is stated, not decided. USDA reports no asparagine and no glutamine — acid hydrolysis
     converts them — so no food can carry all twenty. Each row states how many of the twenty it carries and
     which letters it lacks; which amino acids a score requires is the match step's decision, not this one.
  U6 No cross-archive dedup. A food measured in two archives is two rows with two fdc_ids and its archive
     named. Which to prefer is a match-step decision; `precedence` in config is recorded, not applied here.
  U7 Nothing in data/ or outputs/ is deleted, renamed or moved by this stage (the standing rule, D61).

Usage
  python pipeline/usda.py            # build
  python pipeline/usda.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import zipfile
from pathlib import Path

from pipeline import common
from pipeline.aggregate import strip_generated_line, tsv_text
from pipeline.mass_fractions import sha256_path

AA = list(common.AMINO_ACIDS)
PLACEHOLDER = "___"
OUT_DIR = common.USDA_OUT_DIR
ARCHIVE_DIR = common.USDA_DIR
MANIFEST = ARCHIVE_DIR / "manifest.ini"
MAIN_TABLE = "amino_acids_per_food.tsv"

INPUTS = {
    "usda_config": common.USDA_CONFIG_INI,
    "symbols": common.AMINO_ACID_SYMBOLS_INI,
}


def tsv_rows(path: Path) -> list[dict[str, str]]:
    """Rows of one of this stage's TSVs, skipping the '#' header block."""
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    cols = lines[0].split("\t")
    return [dict(zip(cols, ln.split("\t"))) for ln in lines[1:]]


def _rel(p: Path) -> str:
    """Path relative to the repository root where it is under it, else the file name (tests)."""
    try:
        return p.relative_to(common.REPO_ROOT).as_posix()
    except ValueError:
        return p.name


class Stop(SystemExit):
    def __init__(self, msg: str):
        super().__init__(f"[STOP] {msg}")


# --------------------------------------------------------------------------- reading an archive

def _csv_rows(zf: zipfile.ZipFile, member: str):
    """Rows of one CSV inside the archive, as dicts. FoodData Central ships its CSVs inside a
    single top-level folder; members are matched on the file name alone."""
    with zf.open(member) as raw:
        text = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
        for row in csv.DictReader(text):
            yield row


def _member(zf: zipfile.ZipFile, name: str) -> str | None:
    """The archive member whose file name is `name`, or None."""
    for m in zf.namelist():
        if m.endswith("/"):
            continue
        if Path(m).name.lower() == name.lower():
            return m
    return None


def archive_data_types(zf: zipfile.ZipFile) -> dict[str, int]:
    """The data_type values in the archive's food.csv, with row counts (rule U1)."""
    member = _member(zf, "food.csv")
    if member is None:
        raise Stop("an archive in data/usda/ has no food.csv; is it a FoodData Central CSV download?")
    counts: dict[str, int] = {}
    for row in _csv_rows(zf, member):
        dt = (row.get("data_type") or "").strip()
        counts[dt] = counts.get(dt, 0) + 1
    return counts


def identify_archives(log, wanted: dict[str, dict]) -> dict[str, dict]:
    """Match every zip in data/usda/ to a configured data type by reading it (rule U1)."""
    if not ARCHIVE_DIR.exists():
        raise Stop(f"{_rel(ARCHIVE_DIR)} does not exist; "
                   "download the CSV archives from https://fdc.nal.usda.gov/download-datasets and put the zips there")
    zips = sorted(p for p in ARCHIVE_DIR.glob("*.zip"))
    if not zips:
        raise Stop(f"no .zip in {_rel(ARCHIVE_DIR)}; "
                   f"expected the CSV downloads for: {', '.join(sorted(wanted))}")
    found: dict[str, dict] = {}
    for path in zips:
        try:
            with zipfile.ZipFile(path) as zf:
                counts = archive_data_types(zf)
        except zipfile.BadZipFile:
            log.warning("%s is not a readable zip; skipped", path.name)
            continue
        mine = {dt: n for dt, n in counts.items() if dt in wanted}
        if not mine:
            log.info("%s holds data types %s — none configured; skipped", path.name, ", ".join(sorted(counts)))
            continue
        dt = max(mine, key=lambda d: mine[d])
        if dt in found:
            raise Stop(f"two archives in data/usda/ hold data type {dt}: "
                       f"{found[dt]['file']} and {path.name}. Keep one.")
        found[dt] = {"path": path, "file": path.name, "sha256": sha256_path(path),
                     "bytes": path.stat().st_size, "food_rows": counts[dt]}
        log.info("%s -> %s (%d rows in food.csv)", path.name, dt, counts[dt])
    missing = [dt for dt in wanted if dt not in found]
    if missing:
        have = ", ".join(p.name for p in zips) or "nothing"
        raise Stop("no archive in data/usda/ holds data type(s): " + ", ".join(missing) +
                   f". Found: {have}. Download the CSV (not JSON) archive from "
                   "https://fdc.nal.usda.gov/download-datasets and drop the zip in unrenamed.")
    return found


def membership(zf: zipfile.ZipFile, data_type: str, log) -> tuple[set[str], str]:
    """The fdc_ids that count for this archive, and the rule used (U2)."""
    member = _member(zf, f"{data_type}.csv")
    if member is not None:
        ids = {(r.get("fdc_id") or "").strip() for r in _csv_rows(zf, member)}
        ids.discard("")
        log.info("%s: membership from %s.csv — %d foods", data_type, data_type, len(ids))
        return ids, f"{data_type}.csv"
    ids = {(r.get("fdc_id") or "").strip() for r in _csv_rows(zf, _member(zf, "food.csv"))
           if (r.get("data_type") or "").strip() == data_type}
    ids.discard("")
    log.warning("%s: no %s.csv in the archive; fell back to the data_type column of food.csv — %d foods",
                data_type, data_type, len(ids))
    return ids, "food.csv data_type column (fallback)"


# --------------------------------------------------------------------------- the nutrient join

def load_symbols() -> dict[str, dict[str, str]]:
    """{letter: {trivial_name, three_letter}} from the IUPAC table parsed in the composition step."""
    cp = common.read_ini(INPUTS["symbols"])
    out = {}
    for a in AA:
        if a not in cp:
            raise Stop(f"amino acid {a} missing from {INPUTS['symbols'].name}")
        out[a] = {"trivial_name": cp[a].get("trivial_name", "").strip(),
                  "three_letter": cp[a].get("three_letter", "").strip()}
    return out


def nutrient_map(zf: zipfile.ZipFile, symbols: dict, extras: dict[str, str], units: dict[str, str],
                 log) -> tuple[dict[str, dict], list[list]]:
    """Rule U3: {nutrient_id: {carried_as, name, unit}} plus the rows of the audit table.

    carried_as is a one-letter symbol for one of the twenty, or the config key for an extra."""
    member = _member(zf, "nutrient.csv")
    if member is None:
        raise Stop("an archive has no nutrient.csv")
    by_name: dict[str, list[dict]] = {}
    for row in _csv_rows(zf, member):
        by_name.setdefault((row.get("name") or "").strip().lower(), []).append(row)

    carried: dict[str, dict] = {}
    audit: list[list] = []
    flags: list[str] = []

    def take(key: str, wanted_name: str, kind: str, expect_unit: str):
        rows = by_name.get(wanted_name.strip().lower(), [])
        if not rows:
            audit.append(["", wanted_name, "", kind, key, "absent from this archive"])
            return
        for row in rows:
            nid, unit = (row.get("id") or "").strip(), (row.get("unit_name") or "").strip()
            if expect_unit and unit.upper() != expect_unit.upper():
                flags.append(f"nutrient {nid} {row.get('name')!r} is in {unit}, config expects {expect_unit}; not carried")
                audit.append([nid, row.get("name", ""), unit, kind, key, f"NOT carried: unit is {unit}, expected {expect_unit}"])
                continue
            carried[nid] = {"carried_as": key, "kind": kind, "name": row.get("name", ""), "unit": unit}
            audit.append([nid, row.get("name", ""), unit, kind, key, "carried"])

    for letter in AA:
        take(letter, symbols[letter]["trivial_name"], "amino_acid", units.get("amino_acids", ""))
    for key, name in extras.items():
        take(key, name, "extra", units.get(key, ""))
    for f in flags:
        log.warning(f)
    return carried, audit, flags


# --------------------------------------------------------------------------- the pivot

VALUE_EXTRAS = ["data_points", "min", "max", "median"]


def food_rows(zf: zipfile.ZipFile, ids: set[str], carried: dict[str, dict], log) -> dict[str, dict]:
    """Rule U4: {fdc_id: {carried_as: {amount, data_points, min, max, median}}}, streamed."""
    member = _member(zf, "food_nutrient.csv")
    if member is None:
        raise Stop("an archive has no food_nutrient.csv")
    with zf.open(member) as raw:
        text = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
        reader = csv.DictReader(text)
        cols = reader.fieldnames or []
        log.info("food_nutrient.csv columns: %s", ", ".join(cols))
        have_extra = [c for c in VALUE_EXTRAS if c in cols]
        out: dict[str, dict] = {}
        for row in reader:
            nid = (row.get("nutrient_id") or "").strip()
            spec = carried.get(nid)
            if spec is None:
                continue
            fid = (row.get("fdc_id") or "").strip()
            if fid not in ids:
                continue
            amount = (row.get("amount") or "").strip()
            if amount == "":
                continue
            rec = {"amount": amount}
            for c in have_extra:
                v = (row.get(c) or "").strip()
                if v:
                    rec[c] = v
            out.setdefault(fid, {})[spec["carried_as"]] = rec
    return out


def food_descriptions(zf: zipfile.ZipFile, ids: set[str]) -> dict[str, dict]:
    out = {}
    for row in _csv_rows(zf, _member(zf, "food.csv")):
        fid = (row.get("fdc_id") or "").strip()
        if fid in ids:
            out[fid] = {"description": (row.get("description") or "").strip(),
                        "food_category_id": (row.get("food_category_id") or "").strip(),
                        "publication_date": (row.get("publication_date") or "").strip()}
    return out


def food_categories(zf: zipfile.ZipFile) -> dict[str, str]:
    member = _member(zf, "food_category.csv")
    if member is None:
        return {}
    return {(r.get("id") or "").strip(): (r.get("description") or "").strip() for r in _csv_rows(zf, member)}


# --------------------------------------------------------------------------- build

def load_config() -> tuple[dict[str, dict], dict[str, str], dict[str, str], dict]:
    cp = common.read_ini(INPUTS["usda_config"])
    archives: dict[str, dict] = {}
    for sec in cp.sections():
        if not sec.startswith("archive."):
            continue
        s = cp[sec]
        dt = s.get("data_type", "").strip()
        if not dt or dt == PLACEHOLDER:
            raise Stop(f"{INPUTS['usda_config'].name} [{sec}] data_type is not filled in")
        archives[dt] = {"section": sec, "release": s.get("release", "").strip(),
                        "published_name": s.get("published_name", "").strip(),
                        "url": s.get("url", "").strip(), "obtained": s.get("obtained", "").strip(),
                        "precedence": s.get("precedence", "").strip(), "note": s.get("note", "").strip()}
        for k in ("release", "published_name"):
            if archives[dt][k] == PLACEHOLDER or not archives[dt][k]:
                raise Stop(f"{INPUTS['usda_config'].name} [{sec}] {k} is not filled in")
    if not archives:
        raise Stop(f"{INPUTS['usda_config'].name} names no [archive.*]")
    extras = {k: v.strip() for k, v in cp["carried_nutrients"].items()} if "carried_nutrients" in cp else {}
    units = {k: v.strip() for k, v in cp["units"].items()} if "units" in cp else {}
    meta = dict(cp["meta"]) if "meta" in cp else {}
    return archives, extras, units, meta


def read_manifest() -> "common.configparser.ConfigParser":
    cp = common.new_ini()
    if MANIFEST.exists():
        cp.read(MANIFEST, encoding="utf-8")
    return cp


def build(log) -> tuple[dict[str, str], str]:
    """Returns {file name: text} for outputs/usda/, and the manifest text."""
    archives, extras, units, meta = load_config()
    symbols = load_symbols()
    found = identify_archives(log, archives)

    hashes = {k: f"{_rel(p)} sha256 {sha256_path(p)}" for k, p in INPUTS.items() if p.exists()}
    archive_lines = []
    for dt in sorted(found):
        f, a = found[dt], archives[dt]
        archive_lines.append(f"archive.{dt} = {f['file']} sha256 {f['sha256']} ({a['published_name']}, "
                             f"release {a['release']}, obtained {a['obtained']}, {f['bytes']} bytes)")
    header_common = [f"generated = {common.iso_now()}"] + [f"input.{k} = {v}" for k, v in hashes.items()] + archive_lines

    manifest = read_manifest()
    rows: list[list] = []
    empty_rows: list[list] = []
    audit_rows: list[list] = []
    flags: list[str] = []
    summary = common.new_ini()
    summary.add_section("run")
    per_letter_counts = {a: 0 for a in AA}
    extra_cols = list(extras)

    for dt in sorted(found, key=lambda d: int(archives[d]["precedence"] or 99)):
        f, a = found[dt], archives[dt]
        with zipfile.ZipFile(f["path"]) as zf:
            ids, rule = membership(zf, dt, log)
            carried, audit, unit_flags = nutrient_map(zf, symbols, extras, units, log)
            flags += [f"{dt}: {x}" for x in unit_flags]
            for r in audit:
                audit_rows.append([dt, a["release"]] + r)
            values = food_rows(zf, ids, carried, log)
            desc = food_descriptions(zf, ids)
            cats = food_categories(zf)

        n_with_aa = 0
        for fid in sorted(ids, key=lambda x: int(x) if x.isdigit() else 0):
            v = values.get(fid, {})
            present = [a2 for a2 in AA if a2 in v]
            d = desc.get(fid, {})
            base = [fid, dt, a["release"], d.get("description", ""),
                    cats.get(d.get("food_category_id", ""), ""), d.get("publication_date", "")]
            if not present:
                empty_rows.append(base[:5] + [str(len(v))])
                continue
            n_with_aa += 1
            for a2 in present:
                per_letter_counts[a2] += 1
            row = list(base)
            for key in ("protein", "nitrogen"):
                row.append(v.get(key, {}).get("amount", "") if key in extras else "")
            for a2 in AA:
                row.append(v.get(a2, {}).get("amount", ""))
            for key in extra_cols:
                if key in ("protein", "nitrogen"):
                    continue
                row.append(v.get(key, {}).get("amount", ""))
            row.append(str(len(present)))
            row.append("".join(a2 for a2 in AA if a2 not in v))
            # the per-value detail USDA gives where it gives it: sample count of the amino acid rows
            pts = {a2: v[a2].get("data_points", "") for a2 in present if v[a2].get("data_points")}
            row.append(min(pts.values(), key=lambda x: int(x)) if pts and all(x.isdigit() for x in pts.values()) else "")
            rows.append(row)

        summary.add_section(f"archive.{dt}")
        summary[f"archive.{dt}"].update({
            "file": f["file"], "sha256": f["sha256"], "bytes": str(f["bytes"]),
            "published_name": a["published_name"], "release": a["release"], "url": a["url"],
            "obtained": a["obtained"], "precedence": a["precedence"] or "",
            "membership_rule": rule, "foods_in_membership": str(len(ids)),
            "foods_with_amino_acids": str(n_with_aa),
            "foods_without_amino_acids": str(len(ids) - n_with_aa),
            "amino_acids_carried": str(sum(1 for c in carried.values() if c["kind"] == "amino_acid")),
        })
        sec = f"archive.{dt}"
        if not manifest.has_section(sec):
            manifest.add_section(sec)
        first_seen = manifest[sec].get("first_seen", "") if manifest.has_section(sec) else ""
        same = manifest[sec].get("sha256", "") == f["sha256"]
        manifest[sec].update({
            "file": f["file"], "sha256": f["sha256"], "bytes": str(f["bytes"]), "data_type": dt,
            "release": a["release"], "published_name": a["published_name"], "url": a["url"],
            "obtained": a["obtained"],
            "first_seen": first_seen if (same and first_seen) else common.iso_now(),
        })
        log.info("%s: %d foods in membership, %d with amino acids, %d without",
                 dt, len(ids), n_with_aa, len(ids) - n_with_aa)

    # ---- the tables
    cols = (["fdc_id", "data_type", "release", "description", "food_category", "publication_date",
             "protein_g_per_100g", "nitrogen_g_per_100g"]
            + [f"{a2}_g_per_100g" for a2 in AA]
            + [f"{k}_g_per_100g" for k in extra_cols if k not in ("protein", "nitrogen")]
            + ["amino_acids_present", "amino_acids_absent", "min_data_points"])
    files: dict[str, str] = {}
    files[MAIN_TABLE] = tsv_text(
        header_common + [
            "USDA FoodData Central, one row per food: the amino acids the archive reports, in grams per 100 g of FOOD "
            "(not per 100 g protein), as published. Nothing converted, scaled or filled in (U4).",
            "U5: USDA reports no asparagine (N) and no glutamine (Q) — acid hydrolysis converts them — so no food "
            "carries all twenty; amino_acids_absent lists the letters this food does not carry. Which amino acids a "
            "score requires is the match step's decision.",
            "U6: a food measured in two archives appears twice, with its data_type named; no food is dropped or "
            "preferred here. cystine and cysteine: USDA reports one or the other per food, never both; whether a "
            "cystine value converts to cysteine is an open decision (see config/usda_food_data.ini).",
            "min_data_points: the smallest sample count USDA gives among this food's amino acid rows, where it gives any.",
        ], cols, rows, tool="usda.py")

    files["foods_without_amino_acids.tsv"] = tsv_text(
        header_common + ["foods in each archive's own membership list that carry no amino acid row at all, so that "
                         "what is not in amino_acids_per_food.tsv is on the record (U2, U5). "
                         "other_nutrients_present counts the carried non-amino-acid nutrients (protein, nitrogen, ...)."],
        ["fdc_id", "data_type", "release", "description", "food_category", "other_nutrients_present"],
        empty_rows, tool="usda.py")

    files["usda_nutrient_map.tsv"] = tsv_text(
        header_common + ["rule U3, the audit trail of the name join: every nutrient this stage looked for, the name and "
                         "unit the archive gives it, and what it was carried as. The twenty are matched to the "
                         "IUPAC-IUBMB trivial names in data/iupac/amino_acid_symbols.ini; nothing is mapped by hand."],
        ["data_type", "release", "nutrient_id", "nutrient_name", "unit", "kind", "carried_as", "outcome"],
        audit_rows, tool="usda.py")

    summary["run"]["foods_with_amino_acids"] = str(len(rows))
    summary["run"]["foods_without_amino_acids"] = str(len(empty_rows))
    summary["run"]["flags"] = str(len(flags))
    summary.add_section("amino_acid_coverage")
    for a2 in AA:
        summary["amino_acid_coverage"][a2] = str(per_letter_counts[a2])
    if flags:
        summary.add_section("flag")
        for i, f2 in enumerate(flags, 1):
            summary["flag"][f"f{i}"] = f2
    files["usda_archive_summary.ini"] = common.render_ini(
        summary, ["GENERATED by usda.py. DO NOT EDIT BY HAND."] + header_common)

    manifest_text = common.render_ini(
        manifest,
        ["manifest.ini — GENERATED by usda.py. DO NOT EDIT BY HAND.",
         "The record of which USDA FoodData Central archives this repository was run against.",
         "The .zip files themselves are gitignored — this repository is not a mirror of USDA's data.",
         "Download the same release from https://fdc.nal.usda.gov/download-datasets (CSV) and the hash must match."])
    return files, manifest_text


# --------------------------------------------------------------------------- main

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("usda_check" if args.check else "usda")

    files, manifest_text = build(log)

    if args.check:
        differ = []
        for name, text in files.items():
            path = OUT_DIR / name
            if not path.exists() or strip_generated_line(path.read_text(encoding="utf-8")) != strip_generated_line(text):
                differ.append(name)
        if differ:
            log.error("check: %d file(s) differ from a fresh build: %s", len(differ), differ)
            return 1
        log.info("check: USDA outputs are current")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (OUT_DIR / name).write_text(text, encoding="utf-8", newline="\n")
        log.info("wrote %s", _rel(OUT_DIR / name))
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(manifest_text, encoding="utf-8", newline="\n")
    log.info("wrote %s", _rel(MANIFEST))
    return 0


if __name__ == "__main__":
    sys.exit(main())
