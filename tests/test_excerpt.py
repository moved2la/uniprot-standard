"""excerpt.py on synthetic workbooks. No real literature; no network; no files outside tmp_path.

The tool is what a column map in config is written against, so what it must get right is the
correspondence between a spreadsheet column and its letter: a note that says "columns P-S" has
to land on the same columns the excerpt shows.
"""
import pytest

from pipeline import common
from pipeline import excerpt as ex

openpyxl = pytest.importorskip("openpyxl")


def _book(path, sheets: dict[str, list[list]]):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, rows in sheets.items():
        ws = wb.create_sheet(name)
        for row in rows:
            ws.append(row)
    wb.save(path)
    return path


# a two-row header with a merged-looking group label, as the blood tables have
_S3 = [
    ["Table S3. identified proteins"],                                   # row 1: a title above the header
    ["", "", "Whole erythrocyte", None, "White ghosts"],                  # row 2: group labels, only in
    ["Protein IDs", "Gene names", "Donor 1", "Donor 2", "Donor 1", "Donor 2"],   # the first cell, and the
    ["P68871", "HBB", 40.5, 41.0, 1.2, 1.1],                             # row is SHORTER than row 3
    ["P69905", "HBA1", 30.0, 31.0, 0.9, 0.8],
    ["P00915", "CA1", 2.0, 2.5, 0.1, 0.2],
]


def test_column_letters_round_trip():
    for i, letter in ((1, "A"), (16, "P"), (26, "Z"), (27, "AA")):
        assert ex.col_letter(i) == letter and ex.col_index(letter) == i - 1


def test_whole_number_floats_lose_their_decimal_point():
    assert ex.cell_text(153.0) == "153" and ex.cell_text(40.5) == "40.5"
    assert ex.cell_text(None) == "" and ex.cell_text("  x  ") == "x"


def test_group_label_reaches_the_columns_it_spans(tmp_path):
    path = _book(tmp_path / "s.xlsx", {"Table S3": _S3})
    comments, header, rows, sheet = ex.sheet_table(path, {"sheet": "Table S3", "header_rows": [2, 3]})
    assert sheet == "Table S3"
    # C and D are the whole-erythrocyte donors, E and F the ghosts: the label carried across
    assert header[2] == "Whole erythrocyte | Donor 1" and header[3] == "Whole erythrocyte | Donor 2"
    assert header[4] == "White ghosts | Donor 1"
    # F is past the end of row 2: the label still reaches it, or the column looks unlabelled
    assert header[5] == "White ghosts | Donor 2"
    assert header[0] == "Protein IDs"                       # no group label above it, so just the name
    assert [r[1] for r in rows] == ["HBB", "HBA1", "CA1"]    # the header rows are not in the body
    assert any(c.startswith("# sheet row 1: Table S3.") for c in comments)   # the title travels
    assert any("C = Whole erythrocyte | Donor 1" in c for c in comments)     # the letter map


def test_top_by_column_letter_reads_that_column(tmp_path):
    path = _book(tmp_path / "s.xlsx", {"Table S3": _S3})
    _, header, rows, _ = ex.sheet_table(path, {"sheet": "Table S3", "header_rows": [2, 3]})
    cut = {"kind": "top", "name": "top2_by_col_C", "n": 2, "by_col": "C"}
    kept, rule = ex.apply_cut(cut, header, rows, {})
    assert [r[1] for r in kept] == ["HBB", "HBA1"]           # 40.5 then 30.0, not the ghost columns
    assert "column C (Whole erythrocyte | Donor 1)" in rule


def test_a_sheet_name_that_is_not_there_names_the_ones_that_are(tmp_path):
    path = _book(tmp_path / "s.xlsx", {"Table S1": [["a"]], "Table S3": _S3})
    with pytest.raises(ex.SheetNotFound) as e:
        ex.read_sheet(path, "Table S9")
    assert "'Table S1'" in str(e.value) and "'Table S3'" in str(e.value)
    assert ex.read_sheet(path, "table s3")[0] == "Table S3"   # case and spacing are forgiven


def test_sheet_list_reports_every_sheet(tmp_path):
    path = _book(tmp_path / "s.xlsx", {"Table S1": [["a", "b"]], "Table S3": _S3})
    assert [(n, r) for n, r, _ in ex.workbook_sheets(path)] == [("Table S1", 1), ("Table S3", 6)]


def test_the_sheet_list_says_its_counts_are_the_used_range(tmp_path, monkeypatch):
    """A used range can exceed the content — geyer_2016's mmc3 reports 1289 rows for a 323-row
    table — so the file must not read as a row count."""
    root = tmp_path / "repo"
    (root / "data").mkdir(parents=True)
    _book(root / "data" / "s.xlsx", {"Table S3": _S3})
    monkeypatch.setattr(common, "REPO_ROOT", root)
    monkeypatch.setattr(ex, "EXCERPTS_DIR", root / "excerpts")
    monkeypatch.setattr(ex, "SPEC", [("data/s.xlsx", [{"kind": "sheet_list", "name": "sheets"}])])
    assert ex.main([]) == 0
    stamp = [p for p in (root / "excerpts").iterdir() if p.is_dir()][0]
    text = (stamp / "data" / "s__sheets.tsv").read_text(encoding="utf-8")
    assert "used range" in text and "rows_used_range\tcolumns_used_range" in text


def test_trailing_padding_does_not_travel(tmp_path):
    path = _book(tmp_path / "s.xlsx", {"S": [["a", "b", None, None], ["1", "2", None, None], [None, None, None, None]]})
    _, rows = ex.read_sheet(path, "S")
    assert rows == [["a", "b"], ["1", "2"]]


def test_the_whole_run_writes_workbook_cuts_and_an_index(tmp_path, monkeypatch):
    """main() end to end on a synthetic tree: a workbook cut, a sheet list, a TSV cut, a
    conditional output that is absent, and the INDEX row for each."""
    root = tmp_path / "repo"
    (root / "data").mkdir(parents=True)
    _book(root / "data" / "s.xlsx", {"Table S3": _S3})
    (root / "data" / "t.tsv").write_text("# a header comment\nname\tv\na\t2\nb\t9\n", encoding="utf-8")
    monkeypatch.setattr(common, "REPO_ROOT", root)
    monkeypatch.setattr(ex, "EXCERPTS_DIR", root / "excerpts")
    monkeypatch.setattr(ex, "SPEC", [
        ("data/s.xlsx", [{"kind": "sheet_list", "name": "sheets"},
                         {"kind": "head", "name": "first2", "n": 2, "sheet": "Table S3", "header_rows": [2, 3]},
                         {"kind": "sheet", "name": "all", "sheet": "Table S9"}]),
        ("data/t.tsv", [{"kind": "top", "name": "top1", "n": 1, "by": "v"}]),
        ("data/never_written.tsv", [{"kind": "whole", "optional": True}]),
    ])
    assert ex.main([]) == 0
    out = sorted(p.relative_to(root / "excerpts").as_posix() for p in (root / "excerpts").rglob("*") if p.is_file())
    stamp = [p for p in (root / "excerpts").iterdir() if p.is_dir()][0]
    written = sorted(p.name for p in stamp.rglob("*") if p.is_file())
    assert written == ["INDEX.md", "s__first2.tsv", "s__sheets.tsv", "t__top1.tsv"]
    assert any(name.endswith(".tar.gz") for name in out)

    cut = (stamp / "data" / "s__first2.tsv").read_text(encoding="utf-8").splitlines()
    assert cut[0].startswith("# EXCERPT of data/s.xlsx: sheet 'Table S3', first 2 rows")
    assert cut[-1].split("\t")[1] == "HBA1"                  # two data rows, header rows excluded

    index = (stamp / "INDEX.md").read_text(encoding="utf-8")
    assert "has no sheet 'Table S9'" in index                # named, not crashed
    assert "not written by this run (conditional output)" in index


def test_the_one_xls_in_the_project_opens(tmp_path):
    """hortin_2008's supplementary table is the only .xls here, and it goes through xlrd rather
    than openpyxl. There is no synthetic substitute worth writing: what has to work is that
    file. The test runs where the file is on disk and skips where it is not."""
    pytest.importorskip("xlrd", reason="reading .xls needs xlrd (requirements.txt)")
    path = common.REPO_ROOT / "data/literature/hortin_2008/clinchem.2008.108175-2.xls"
    if not path.exists():
        pytest.skip("hortin_2008's .xls is not on this machine")
    sheets = ex.workbook_sheets(path)
    assert sheets, "the workbook reports no sheets"
    name, rows = ex.read_sheet(path, sheets[0][0])
    assert name == sheets[0][0] and rows and all(isinstance(c, str) for c in rows[0])