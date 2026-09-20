#!/usr/bin/env python3
"""
score.py — tooling: Match Rate for a workbook of foods entered by hand.

NOT a pipeline stage. It is not in `common.COMMAND_ORDER`, nothing downstream reads what it
writes, and no currency test rebuilds it. It is the ad-hoc path: transcribe a label, paste a
supplier's analysis, enter a finished blend, and get the Match Rate back without editing config,
without rerunning `match` over every USDA food, and without any of it entering the repository's
record. `scoring/` is gitignored for that reason (M7).

    1. put a workbook in  scoring/
    2. python run.py score
    3. read  scoring/<name>_scored.xlsx

The arithmetic is `match.match_rate()` — the same function `pipeline/match.py` scores USDA with
(M6). Nothing about the formula is reimplemented here; what is here is the workbook in, the
workbook out, and the imputation rule for an amino acid a label does not print (M8).

The reference is whatever `config/match_rate.ini [formula] summary_reference` names — the most
complete standard built so far (D128). As categories 3 through 11 land, that key moves and this
tool follows it; no reference is named here or in the workbook.

THE SHEET IS TRANSPOSED (M9, D131), the way the author's own spreadsheets are laid out: the field
names run DOWN column A, and each column from B rightwards is one product.

    A                B                C
    label            Whey isolate     Pea isolate
    basis            product          protein
    grams_basis      100              100
    protein_g        90               100
    source           supplier COA     label
    Ala              4.5              4.1
    Arg              2.2              8.4
    ...              ...              ...

Reads
  scoring/*.xlsx                                the workbooks you put there (`_scored` files skipped)
  config/match_rate.ini                         the reference, its scored set, the version label
  config/fao_2013_indispensable_amino_acids.ini the nine headings (D62), reduced per D86
  data/iupac/amino_acid_symbols.ini             three-letter / trivial name -> letter: nothing is named in code
  the standard table the reference names        its values

Writes (scoring/, beside the input)
  <name>_scored.xlsx   sheet `match_rate`: YOUR sheet, every cell exactly as entered and in your order,
                                           with the calculation rows appended at the BOTTOM
                       sheet `steps`:      the spreadsheet's own walk (M5), transposed the same way —
                                           the step names down column A, one column per product
                       sheet `reference`:  what was scored against, where it was read, its hash, and its
                                           values — the provenance, as rows and columns, not as comments

The input workbook is never written to.

Rows of the input sheet (the name in column A)
  label         required. The product's identity, one per column. A column with no label is carried
                through but not scored.
  basis         required. `product`, `protein` or `eaa` — what the amounts are per.
  grams_basis   required. How many grams of that basis the amounts belong to (100, 250, 31.5 ...).
  protein_g     optional. Protein in that same amount. Blank leaves the density rows blank.
  <Xxx>         the amino acids, three-letter symbol or trivial name (His, Ile, Trp / Histidine ...),
                in the same unit as grams_basis. BLANK means not reported and is imputed (M8);
                0 means measured as zero and scores zero (M3).
  anything else carried through untouched: source, location, note, a supplier's code.

Rules M7-M10 in docs/conventions.md.

Usage
  python run.py score          # the whole scoring/ folder
  python pipeline/score.py     # the same
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline import common
from pipeline.match import (Stop, load_references, load_scored_sets, match_rate,
                            sheet_steps, symbol_tables)
from pipeline.mass_fractions import sha256_path

AA = list(common.AMINO_ACIDS)
CONFIG_INI = common.MATCH_RATE_INI
SCORED_SUFFIX = "_scored"
GENERATED_BY = "score.py"

REQUIRED_FIELDS = ["label", "basis", "grams_basis"]
OPTIONAL_FIELDS = ["protein_g", "source", "location", "note"]

# What the amounts are per. The score never sees the denominator (M1); the basis decides which
# density rows can be computed, and records what the numbers mean.
BASES = ("product", "protein", "eaa")

# The calculation rows, appended below the author's own. Order is the reading order: the score
# first, then what it was scored on, then the totals.
CALC_ROWS = ["match_rate_percent", "limiting_amino_acid", "scored_amino_acids",
             "imputed_at_reference_share", "not_scored_reason",
             "eaa_total", "eaa_percent_of_protein", "eaa_g_per_100g_product",
             "protein_percent_of_product"]


def _rel(p: Path) -> str:
    try:
        return p.relative_to(common.REPO_ROOT).as_posix()
    except ValueError:
        return p.name


def _col(index0: int) -> str:
    """Excel's name for a zero-based column index, so a message points at a cell the author sees."""
    out, n = "", index0 + 1
    while n > 0:
        n, r = divmod(n - 1, 26)
        out = chr(65 + r) + out
    return out


# --------------------------------------------------------------------------- the reference (one, D128)

def load_summary_reference() -> dict:
    """The reference `[formula] summary_reference` names, with its scored set resolved.

    The same loaders match.py uses, so a reference that scores in the pipeline scores here and a
    rename in config reaches both at once.
    """
    if not CONFIG_INI.exists():
        raise Stop(f"{_rel(CONFIG_INI)} not found")
    cp = common.read_ini(CONFIG_INI)
    try:
        three, names = symbol_tables()
    except FileNotFoundError:
        raise Stop(f"{_rel(common.AMINO_ACID_SYMBOLS_INI)} is not on this machine — it is the IUPAC-IUBMB "
                   f"table every stage resolves amino acid names through. Run `python run.py composition` "
                   f"first (docs/run_order.md).")
    scored_sets = load_scored_sets(cp, three, names)
    refs = load_references(cp, scored_sets, three, names)
    name = cp["formula"].get("summary_reference", "").strip() if "formula" in cp else ""
    if not name:
        raise Stop("[formula] summary_reference is not set in config/match_rate.ini. It names the most "
                   "complete standard built so far (D128), and it is what this tool scores against.")
    hit = [r for r in refs if r["name"] == name]
    if not hit:
        raise Stop(f"[formula] summary_reference {name!r} is not a [reference.*] in {CONFIG_INI.name}")
    ref = dict(hit[0])
    ref["version"] = cp["meta"].get("version", "").strip() if "meta" in cp else ""
    ref["decimals"] = int(cp["output"].get("decimals", "4")) if "output" in cp else 4
    ref["how"] = scored_sets[ref["scored_set"]]["how"]
    ref["symbol_of"] = {v: k for k, v in three.items()}          # letter -> three-letter, lowercased
    ref["three"], ref["names"] = three, names
    return ref


def three_letter(ref: dict, letter: str) -> str:
    return ref["symbol_of"].get(letter, letter).capitalize()


def alphabetical_letters(ref: dict) -> list[str]:
    """The twenty letters in alphabetical order of their IUPAC trivial name (the order labels and
    supplier sheets print: Ala, Arg, Asn, Asp, Cys, Glu, Gln, Gly, His ...). Read from the table,
    never listed here."""
    by_name = {v: k for k, v in ref["names"].items()}
    return sorted(AA, key=lambda a: by_name[a])


# --------------------------------------------------------------------------- the workbook in

def _cell(v) -> str:
    return "" if v is None else str(v).strip()


def read_sheet(path: Path) -> list[list]:
    """The grid of the sheet named `foods`, else the workbook's first sheet, padded to one width.

    Values come back as openpyxl read them — a number stays a number — so the output can write
    the author's own cells back unchanged.
    """
    from openpyxl import load_workbook
    try:
        wb = load_workbook(path, data_only=True, read_only=True)
    except Exception as e:
        raise Stop(f"{_rel(path)} could not be opened as a workbook: {e}")
    try:
        ws = wb["foods"] if "foods" in wb.sheetnames else wb.worksheets[0]
        rows = [list(r) for r in ws.iter_rows(values_only=True)]
    finally:
        wb.close()
    rows = [r for r in rows if any(_cell(c) for c in r)]
    if not rows:
        return []
    width = max(len(r) for r in rows)
    while width > 1 and not any(_cell(r[width - 1]) for r in rows if len(r) >= width):
        width -= 1
    # A sheet with nothing but column A — the unfilled blank sheet — is not an error: it is
    # skipped with a line in the log, so one unfilled template cannot block the folder.
    return [list(r[:width]) + [None] * max(0, width - len(r)) for r in rows]


def field_index(grid: list[list], ref: dict, where: str) -> tuple[dict[str, int], dict[str, int]]:
    """({field name: row index}, {letter: row index}) from column A.

    A row whose column A resolves through the IUPAC-IUBMB table is an amino acid row; every other
    row is a field, carried by whatever the author called it. Nothing is named in code.
    """
    fields: dict[str, int] = {}
    acids: dict[str, int] = {}
    for i, row in enumerate(grid):
        name = _cell(row[0])
        if not name:
            continue
        letter = ref["three"].get(name.lower()) or ref["names"].get(name.lower())
        if letter:
            if letter in acids:
                raise Stop(f"{where}: rows {acids[letter] + 1} and {i + 1} are both "
                           f"{three_letter(ref, letter)}")
            acids[letter] = i
        else:
            key = name.lower()
            if key in fields:
                raise Stop(f"{where}: rows {fields[key] + 1} and {i + 1} are both {name!r}")
            fields[key] = i
    return fields, acids


def _number(raw, where: str, field: str) -> float | None:
    if raw is None or str(raw).strip() == "":
        return None
    try:
        return float(str(raw).strip())
    except ValueError:
        raise Stop(f"{where}: {field} = {raw!r} is not a number")


# --------------------------------------------------------------------------- the imputation (M8)

def impute_unreported(values: dict[str, float | None], reference: dict[str, float],
                      scored: list[str]) -> tuple[dict[str, float], list[str], str]:
    """Fill a scored amino acid the product does not report at the reference's own share (ratio 1.0).

    Returns (the vector to score, the letters imputed, "" or the reason it cannot be scored).

    Arithmetically this is scoring over the amino acids the product DOES report, with both sides
    renormalised over that subset — the same number either way. Doing it as an imputation rather
    than a shrunken scored set is what keeps the walk in `steps` on the full scored set, so the
    sheet reads the same for a product with nine and a product with eight.

    An imputed amino acid can never become the limiting one on its own: the shares of both sides
    sum to one over the scored set, so at least one ratio is at or below 1.0 and it belongs to a
    reported amino acid unless every ratio is exactly 1.0.
    """
    reported = [a for a in scored if values.get(a) is not None]
    missing = [a for a in scored if values.get(a) is None]
    if not reported:
        return {}, missing, "the product reports none of the scored amino acids"
    r_total = sum(reference[a] for a in scored)
    missing_share = sum(reference[a] for a in missing) / r_total
    reported_sum = sum(values[a] for a in reported)
    if reported_sum <= 0:
        return {}, missing, "every scored amino acid the product reports is zero"
    scored_total = reported_sum / (1.0 - missing_share)
    out = {a: float(values[a]) for a in reported}
    for a in missing:
        out[a] = reference[a] / r_total * scored_total
    return out, missing, ""


# --------------------------------------------------------------------------- scoring one workbook

def score_workbook(path: Path, ref: dict) -> dict:
    """{grid, columns, scored} for one input workbook; `scored` is one dict per product column."""
    grid = read_sheet(path)
    if not grid:
        return {"grid": [], "columns": [], "scored": []}
    where = _rel(path)
    fields, acids = field_index(grid, ref, where)
    absent = [f for f in REQUIRED_FIELDS if f not in fields]
    if absent:
        raise Stop(f"{where}: column A has no row named {', '.join(absent)}. Required rows: "
                   f"{', '.join(REQUIRED_FIELDS)}; optional: protein_g; then the amino acids by "
                   f"three-letter symbol, one per row.")
    if not acids:
        raise Stop(f"{where}: no row of column A is an amino acid. They must be three-letter symbols "
                   f"(His, Ile, Leu ...) or trivial names (Histidine ...), one per row.")

    width = len(grid[0])
    label_row = grid[fields["label"]]
    columns = [j for j in range(1, width) if _cell(label_row[j])]
    if not columns:
        return {"grid": grid, "columns": [], "scored": []}

    out = []
    for j in columns:
        label = _cell(label_row[j])
        cw = f"{where} column {_col(j)} ({label})"
        basis = _cell(grid[fields["basis"]][j]).lower()
        if basis not in BASES:
            raise Stop(f"{cw}: basis = {_cell(grid[fields['basis']][j])!r}; it must be "
                       f"{', '.join(BASES[:-1])} or {BASES[-1]} — what the amounts are per")
        grams = _number(grid[fields["grams_basis"]][j], cw, "grams_basis")
        if grams is None or grams <= 0:
            raise Stop(f"{cw}: grams_basis must be a number above zero — how many grams of "
                       f"{basis} the amounts belong to")
        protein = _number(grid[fields["protein_g"]][j], cw, "protein_g") if "protein_g" in fields else None
        values = {**{a: None for a in AA},
                  **{a: _number(grid[i][j], cw, three_letter(ref, a)) for a, i in acids.items()}}

        vector, imputed, reason = impute_unreported(values, ref["values"], ref["letters"])
        entry = {"column": j, "label": label, "basis": basis, "grams_basis": grams, "protein": protein,
                 "imputed": imputed, "reason": reason, "result": None, "steps": None, "vector": vector}
        if not reason:
            res = match_rate(vector, ref["values"], ref["letters"])
            if res["score"] is None:
                entry["reason"] = res["reason"]
            else:
                entry["result"] = res
                entry["steps"] = sheet_steps(vector, ref["values"], ref["letters"], res)
        out.append(entry)
    return {"grid": grid, "columns": columns, "scored": out}


# --------------------------------------------------------------------------- the workbook out

def _round(x, d: int):
    return None if x is None else round(float(x), d)


def calculation_cells(e: dict, ref: dict) -> dict[str, object]:
    """{calculation row name: this product's value}."""
    d = ref["decimals"]
    res = e["result"]
    cells: dict[str, object] = {name: None for name in CALC_ROWS}
    cells["scored_amino_acids"] = ", ".join(three_letter(ref, a) for a in ref["letters"])
    cells["imputed_at_reference_share"] = ", ".join(three_letter(ref, a) for a in e["imputed"])
    if res is None:
        cells["limiting_amino_acid"] = ""
        cells["not_scored_reason"] = e["reason"]
        return cells
    eaa = sum(e["vector"][a] for a in ref["letters"])
    protein, grams = e["protein"], e["grams_basis"]
    per_product = e["basis"] == "product"
    cells["match_rate_percent"] = _round(100 * res["score"], d)
    cells["limiting_amino_acid"] = "+".join(three_letter(ref, a) for a in res["limiting"])
    cells["not_scored_reason"] = ""
    cells["eaa_total"] = _round(eaa, d)
    cells["eaa_percent_of_protein"] = _round(100 * eaa / protein, d) if protein else None
    cells["eaa_g_per_100g_product"] = _round(100 * eaa / grams, d) if per_product else None
    cells["protein_percent_of_product"] = _round(100 * protein / grams, d) \
        if (per_product and protein is not None) else None
    return cells


def match_rate_sheet(book: dict, ref: dict) -> list[list]:
    """The author's grid unchanged, with the calculation rows appended at the bottom (M9)."""
    width = len(book["grid"][0])
    rows = [list(r) for r in book["grid"]]
    cells = {e["column"]: calculation_cells(e, ref) for e in book["scored"]}
    rows.append([None] * width)                               # one blank row, then the calculations
    for name in CALC_ROWS:
        row: list = [name] + [None] * (width - 1)
        for col, c in cells.items():
            row[col] = c[name]
        rows.append(row)
    return rows


def steps_sheet(book: dict, ref: dict) -> list[list]:
    """The spreadsheet's walk (M5), transposed the same way: the step names down column A, one
    column per product — the shape of the author's own Base_Match_Rate_Formula.xlsx."""
    d, letters = ref["decimals"], ref["letters"]
    sym = [three_letter(ref, a) for a in letters]
    names = (["label", "match_rate_percent"]
             + [f"{s}_amount" for s in sym]                                 # rows 12-20
             + ["eaa_total", "reference_eaa_total", "percent_of_reference_total"]   # 22, A22, 23
             + [f"step3_scaled_{s}" for s in sym]                           # rows 30-37
             + [f"step7_percent_of_reference_{s}" for s in sym]             # rows 100-107
             + ["step7_min", "limiting_amino_acid"]                         # row 108
             + [f"step10b_{s}" for s in sym]                                # rows 145-152
             + ["step10b_total_need_to_consume", "percent_wasted", "percent_utilized"])   # 153-156
    columns: list[dict[str, object]] = []
    for e in book["scored"]:
        if e["result"] is None:
            continue
        s, res = e["steps"], e["result"]
        c: dict[str, object] = {"label": e["label"],
                                "match_rate_percent": _round(100 * res["score"], d)}
        for a, symbol in zip(letters, sym):
            c[f"{symbol}_amount"] = _round(e["vector"][a], d)
            c[f"step3_scaled_{symbol}"] = _round(s["scaled"][a], d)
            c[f"step7_percent_of_reference_{symbol}"] = _round(res["ratio"][a], d)
            c[f"step10b_{symbol}"] = None if s["final"] is None else _round(s["final"][a], d)
        c["eaa_total"] = _round(s["food_total"], d)
        c["reference_eaa_total"] = _round(s["ref_total"], d)
        c["percent_of_reference_total"] = _round(s["pct_of_ref_total"], d)
        c["step7_min"] = _round(res["score"], d)
        c["limiting_amino_acid"] = "+".join(three_letter(ref, a) for a in res["limiting"])
        c["step10b_total_need_to_consume"] = None if s["final"] is None else _round(s["need"], d)
        c["percent_wasted"] = None if s["final"] is None else _round(s["wasted"], d)
        c["percent_utilized"] = _round(s["utilized"], d)
        columns.append(c)
    return [[name] + [c.get(name) for c in columns] for name in names]


def reference_sheet(ref: dict, source: Path) -> list[list]:
    """The provenance and the reference's own values, as rows and columns. No comment lines."""
    total = sum(ref["values"][a] for a in ref["letters"])
    fields = [
        ("generated", common.iso_now()),
        ("generated_by", GENERATED_BY),
        ("match_rate_version", ref["version"] or "(unlabelled)"),
        ("input_workbook", source.name),
        ("input_workbook_sha256", sha256_path(source)),
        ("reference", ref["name"]),
        ("reference_label", ref["label"]),
        ("reference_read_from", ref["where"]),
        ("reference_sha256", sha256_path(ref["path"]) if ref["path"].exists() else ""),
        ("reference_decisions", ref["decisions"]),
        ("scored_set", ref["scored_set"]),
        ("scored_set_how", ref["how"]),
        ("scored_amino_acids", ", ".join(three_letter(ref, a) for a in ref["letters"])),
        ("config", f"{_rel(CONFIG_INI)} sha256 {sha256_path(CONFIG_INI)}"),
        ("amino_acid_symbols", f"{_rel(common.AMINO_ACID_SYMBOLS_INI)} sha256 "
                               f"{sha256_path(common.AMINO_ACID_SYMBOLS_INI)}"),
        ("formula", "Match Rate = the smallest of (food share / reference share) over the scored "
                    "amino acids, both shares over the scored set only (M1, D87)"),
        ("unreported_amino_acids", "imputed at the reference's own share, so their ratio is 1.0 "
                                   "(M8); the same number as scoring over the reported subset"),
    ]
    rows: list[list] = [["field", "value"]] + [[k, v] for k, v in fields]
    rows.append([])
    rows.append(["amino_acid", "symbol", "reference_amount", "reference_share_percent", "scored"])
    for a in alphabetical_letters(ref):
        value = ref["values"].get(a)          # a transcribed reference carries only its own rows
        share = 100 * value / total if (a in ref["letters"] and value is not None) else None
        rows.append([a, three_letter(ref, a), _round(value, 4),
                     _round(share, 4), "yes" if a in ref["letters"] else ""])
    return rows


def _fit(sheet) -> None:
    """Column A holds the field names, so it is the wide one; the product columns are numbers."""
    from openpyxl.utils import get_column_letter
    for i in range(1, sheet.max_column + 1):
        width = max((len(str(sheet.cell(row=r, column=i).value or ""))
                     for r in range(1, min(sheet.max_row, 80) + 1)), default=10)
        sheet.column_dimensions[get_column_letter(i)].width = min(max(width + 2, 10), 34 if i == 1 else 26)


def write_workbook(out_path: Path, book: dict, ref: dict, source: Path) -> None:
    from openpyxl import Workbook

    wb = Workbook()
    ws = wb.active
    ws.title = "match_rate"
    for r in match_rate_sheet(book, ref):
        ws.append(r)
    ws.freeze_panes = "B1"

    st = wb.create_sheet("steps")
    for r in steps_sheet(book, ref):
        st.append(r)
    st.freeze_panes = "B1"

    rs = wb.create_sheet("reference")
    for r in reference_sheet(ref, source):
        rs.append(r)

    for sheet in (ws, st, rs):
        _fit(sheet)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


def _is_ours(path: Path) -> bool:
    """True when `path` is a workbook this tool wrote — checked before overwriting it."""
    from openpyxl import load_workbook
    wb = None
    try:
        wb = load_workbook(path, data_only=True, read_only=True)
        if "reference" not in wb.sheetnames:
            return False
        for row in wb["reference"].iter_rows(values_only=True):
            if row and _cell(row[0]) == "generated_by":
                return _cell(row[1] if len(row) > 1 else "") == GENERATED_BY
        return False
    except Exception:
        return False
    finally:
        if wb is not None:
            try:
                wb.close()
            except Exception:
                pass


# --------------------------------------------------------------------------- the blank sheet

def write_blank_sheet(path: Path, ref: dict) -> None:
    """A workbook with nothing but column A, written when scoring/ holds no input."""
    from openpyxl import Workbook
    wb = Workbook()
    ws = wb.active
    ws.title = "foods"
    for name in (REQUIRED_FIELDS + OPTIONAL_FIELDS
                 + [three_letter(ref, a) for a in alphabetical_letters(ref)]):
        ws.append([name])
    ws.freeze_panes = "B1"
    _fit(ws)
    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)


# --------------------------------------------------------------------------- run

def input_workbooks(folder: Path) -> list[Path]:
    return sorted(p for p in folder.glob("*.xlsx")
                  if not p.name.startswith("~$") and not p.stem.endswith(SCORED_SUFFIX))


def run(log) -> list[Path]:
    folder = common.scoring_dir()
    folder.mkdir(parents=True, exist_ok=True)
    ref = load_summary_reference()
    log.info("reference %s (%s), scored on %s", ref["name"], ref["where"],
             ", ".join(three_letter(ref, a) for a in ref["letters"]))

    books = input_workbooks(folder)
    if not books:
        blank = folder / "blank_scoring_sheet.xlsx"
        if not blank.exists():
            write_blank_sheet(blank, ref)
            log.info("no workbook in %s; wrote %s with the field names down column A",
                     _rel(folder), _rel(blank))
        else:
            log.info("no workbook to score in %s", _rel(folder))
        return []

    written = []
    for src in books:
        book = score_workbook(src, ref)
        if not book["scored"]:
            log.info("%s: no product column (column A is the field names; B rightwards are the "
                     "products, each needing a label); nothing scored", _rel(src))
            continue
        out = src.with_name(f"{src.stem}{SCORED_SUFFIX}.xlsx")
        if out.exists() and not _is_ours(out):
            raise Stop(f"{_rel(out)} exists and was not written by {GENERATED_BY}. Move or rename it; "
                       f"this tool will not overwrite a workbook it did not write.")
        write_workbook(out, book, ref, src)
        ok = [e for e in book["scored"] if e["result"] is not None]
        imputed = sum(1 for e in ok if e["imputed"])
        log.info("%s -> %s: %d product(s), %d scored, %d with an imputed amino acid",
                 _rel(src), _rel(out), len(book["scored"]), len(ok), imputed)
        for e in book["scored"]:
            if e["result"] is None:
                log.warning("  column %s %s: not scored — %s", _col(e["column"]), e["label"], e["reason"])
            else:
                log.info("  column %s %s: %.2f %%, limiting %s%s", _col(e["column"]), e["label"],
                         100 * e["result"]["score"],
                         "+".join(three_letter(ref, a) for a in e["result"]["limiting"]),
                         (" (imputed " + ", ".join(three_letter(ref, a) for a in e["imputed"]) + ")")
                         if e["imputed"] else "")
        written.append(out)
    return written


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    log = common.make_logger("score")
    run(log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
