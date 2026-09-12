"""literature_inventory.py on synthetic files: a workbook, a text table, and a zip
holding both. No real biology; no network. The .rar path needs an external
extractor and is not exercised here — it is exercised by the stage itself on the
workstation, where a missing extractor is a stop, not a silent skip."""
import io
import logging
import zipfile

from openpyxl import Workbook

from pipeline import literature_inventory as inv


class Opts:
    preview_rows = 3
    preview_cols = 4
    cell_chars = 8


def _workbook_bytes() -> bytes:
    wb = Workbook()
    ws = wb.active
    ws.title = "First"
    ws.append(["title row", None, None])
    ws.append([None])
    ws.append(["Name ", "Median X ", "Valid  Y", "a-very-long-header-string"])
    ws.append(["g1", 1.5, 3, "x"])
    ws.append(["g2", "NaN", 0, "y"])
    ws2 = wb.create_sheet("Second")
    ws2.append(["only", "two"])
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def _ctx(**over):
    base = {"source_id": "s", "file_key": "s.1", "path": "data/literature/s/f", "member": "",
            "bytes": 0, "sha256": "", "match": "true"}
    base.update(over)
    return base


def test_workbook_counts_and_preview():
    fr, hr = inv.inspect_workbook(_workbook_bytes(), _ctx(), Opts)
    by_sheet = {r[8]: r for r in fr}
    assert by_sheet["First"][9] == 5 and by_sheet["First"][10] == 4      # exact rows, widest row
    assert by_sheet["Second"][9] == 1 and by_sheet["Second"][10] == 2
    first = [r for r in hr if r[3] == "First"]
    assert [r[4] for r in first] == [1, 2, 3]                           # preview rows only
    assert first[2][5] == "Name  | Median … | Valid  Y | a-very-…"  # whitespace kept, cells over 8 chars truncated


def test_text_table_counts():
    data = b"a\tb\tc\n1\t2\t3\n4\t5\n"
    fr, hr = inv.inspect_text_table(data, _ctx(member="t.tsv"), Opts)
    assert fr[0][4] == "text_table" and fr[0][9] == 3 and fr[0][10] == 3
    assert hr[0][5] == "a | b | c"


def test_zip_members_hashed_and_tabular_members_inspected(tmp_path):
    zpath = tmp_path / "bundle.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("inner/book.xlsx", _workbook_bytes())
        zf.writestr("notes.pdf", b"%PDF-1.4 synthetic")
    fr, hr = inv.inspect_archive(zpath, _ctx(path="bundle.zip"), Opts, None, logging.getLogger("t"))
    members = [r for r in fr if r[4] == "archive_member"]
    assert {r[3] for r in members} == {"inner/book.xlsx", "notes.pdf"}
    assert all(len(r[6]) == 64 for r in members)                          # every member hashed
    sheets = [r for r in fr if r[4] == "sheet"]
    assert {r[8] for r in sheets} == {"First", "Second"} and all(r[3] == "inner/book.xlsx" for r in sheets)
    assert any(r[2] == "inner/book.xlsx" for r in hr)


def test_rar_without_tool_is_a_stop(tmp_path):
    rpath = tmp_path / "x.rar"
    rpath.write_bytes(b"Rar!\x1a\x07\x00synthetic")
    try:
        inv.inspect_archive(rpath, _ctx(path="x.rar"), Opts, None, logging.getLogger("t"))
    except SystemExit as e:
        assert "STOP" in str(e)
    else:
        raise AssertionError("a .rar with no extractor must stop the stage")
