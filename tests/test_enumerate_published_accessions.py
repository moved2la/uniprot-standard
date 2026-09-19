"""enumerate_published_pool on synthetic inputs: a workbook written by the test, a fake UniProt
lookup, and a decisions file with one pool. Accessions, genes and keyword ids are invented
(rules P1-P6, D107; keratin rule D109)."""
import logging

from openpyxl import Workbook

from pipeline import common, enumerate_pool as ep

KW = "KW-0416"


class FakeClient:
    release = "TEST"

    def __init__(self, entries):
        self.entries = entries          # accession -> dict of the lookup fields
        self.queries = []

    def search_tsv(self, query, fields):
        self.queries.append(query)
        wanted = [t.split(":", 1)[1].rstrip(")") for t in query.split(" OR ")]
        wanted = [w.lstrip("(") for w in wanted]
        return [dict(self.entries[a], accession=a) for a in wanted if a in self.entries]


def _workbook(path):
    wb = Workbook()
    ws = wb.active
    ws.title = "T"
    ws.append(["title row"])
    ws.append(["Protein IDs", "Protein names", "Gene names", "x", "x", "x", "x", "x", "x", "x", "x", "x", "x", "x", "x",
               "Donor 1", "Donor 2", "Donor 3", "Donor 4"])
    rows = [
        ("X00001;X00001-2;T00001", "one", "GA", 40, 42, 38, 40),               # member, first token
        ("X00002-3;X00002", "two-iso", "GB", 10, 10, 10, 10),                    # isoform suffix stripped -> member
        ("T00003;X00003", "three", "GC", 5, 5, 5, 5),                            # first token unreviewed -> member via token 2
        ("CON__X00004;X00004", "four", "GD", 1, 1, 1, 1),                        # contaminant group
        ("X00005", "five", "KRTX", 3, 0.1, 0.1, 0.1),                            # keratin keyword
        ("X00006", "six", "GF;ALBX", 2, 2, 2, 2),                                # named contaminant via gene cell
        ("T00007;T00008", "seven", "GG", 1, 1, 1, 1),                            # no reviewed human entry
        ("X00001;Q99999", "one again", "GA", 0.5, 0.5, 0.5, 0.5),                # duplicate entry
        ("X00009", "nine", "GH", 0.2, 0.2, 0.2, 0.2),                            # not returned by UniProt
        ("M00010", "mouse", "GI", 0.1, 0.1, 0.1, 0.1),                           # reviewed but not human
        ("", "", "", "", "", "", ""),                                             # empty identity -> skipped
    ]
    for r in rows:
        ws.append(list(r[:3]) + [""] * 12 + list(r[3:]))
    wb.save(path)


def _setup(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(common, "CONFIG_DIR", tmp_path / "config")
    monkeypatch.setattr(common, "DATA_DIR", tmp_path / "data")
    monkeypatch.setattr(common, "OUTPUTS_DIR", tmp_path / "outputs")
    (tmp_path / "config" / "cat").mkdir(parents=True)
    (tmp_path / "data" / "literature" / "src").mkdir(parents=True)
    xlsx = tmp_path / "data" / "literature" / "src" / "t.xlsx"
    _workbook(xlsx)
    man = common.new_ini()
    man["file.src.file.1"] = {"path": "data/literature/src/t.xlsx", "sha256": common.sha256_file(xlsx)}
    common.write_ini(man, tmp_path / "data" / "literature" / "manifest.ini", ["test"])
    dec = common.new_ini()
    dec["scope"] = {"category": "cat"}
    dec["columns.src.file.1"] = {"sheet": "T", "header_row": "2", "identity_column": "A", "identity_kind": "accession_list",
                                 "gene_column": "C", "share_columns": "P, Q, R, S", "share_kind": "mass_share"}
    common.write_ini(dec, tmp_path / "config" / "cat" / "mass_fraction_decisions.ini", ["test"])
    psd = common.new_ini()
    psd["scope"] = {"category": "cat"}
    psd["pool.red"] = {"dataset_source": "src", "dataset_file": "file.1", "named_contaminant_genes": "ALBX, ZZZ",
                       "pool_decision": "D107", "contaminant_decision": "D98"}
    psd["contaminant_rules"] = {"keratin_keyword_id": KW}
    common.write_ini(psd, tmp_path / "config" / "cat" / "protein_set_decisions.ini", ["test"])
    return psd


def _entries():
    rev = {"reviewed": "reviewed", "organism_id": "9606", "gene_primary": "G", "protein_name": "p", "length": "10", "keywordid": ""}
    unrev = dict(rev, reviewed="unreviewed")
    return {
        "X00001": rev, "X00002": rev, "X00003": rev, "X00004": rev, "X00006": rev, "Q99999": rev,
        "X00005": dict(rev, keywordid="KW-0001;" + KW),
        "T00001": unrev, "T00003": unrev, "T00007": unrev, "T00008": unrev,
        "M00010": dict(rev, organism_id="10090"),
        # X00009 is deliberately absent: UniProt returns nothing for it
    }


def test_token_classification():
    assert ep.classify_token("P11277-2") == ("accession", "P11277")
    assert ep.classify_token("CON__P35908") == ("contaminant_marker", "")
    assert ep.classify_token("REV__P1") == ("decoy", "")
    assert ep.classify_token(" ") == ("blank", "")


def test_published_pool_outcomes(tmp_path, monkeypatch):
    psd = _setup(tmp_path, monkeypatch)
    log = logging.getLogger("t")
    client = FakeClient(_entries())
    members, record = ep.enumerate_published_pool(client, "cat", "red", psd["pool.red"], psd["contaminant_rules"], log)
    assert members == {"X00001", "X00002", "X00003"}
    rows = {r["identity_cell"]: r for r in common.read_tsv(tmp_path / "data" / "cat" / "red_dataset_rows.tsv")}
    assert rows["X00001;X00001-2;T00001"]["outcome"] == "member" and rows["X00001;X00001-2;T00001"]["entry_token_rank"] == "1"
    assert rows["X00002-3;X00002"]["outcome"] == "member" and rows["X00002-3;X00002"]["entry"] == "X00002"
    assert rows["T00003;X00003"]["outcome"] == "member_via_later_token" and rows["T00003;X00003"]["entry_token_rank"] == "2"
    assert rows["CON__X00004;X00004"]["outcome"] == "contaminant_group"
    assert rows["X00005"]["outcome"] == "keratin"
    assert rows["X00006"]["outcome"] == "named_contaminant"
    assert rows["T00007;T00008"]["outcome"] == "no_reviewed_human_entry"
    assert rows["X00001;Q99999"]["outcome"] == "member_duplicate_entry"
    assert rows["X00009"]["outcome"] == "no_reviewed_human_entry" and "not_returned" in rows["X00009"]["tokens_checked"]
    assert rows["M00010"]["outcome"] == "no_reviewed_human_entry" and "/10090" in rows["M00010"]["tokens_checked"]
    assert "" not in rows                                     # the empty row is skipped, not listed
    # the excluded listing is sorted by share, largest first, and carries the share
    exc = common.read_tsv(tmp_path / "data" / "cat" / "red_rows_excluded_with_share.tsv")
    assert [r["outcome"] for r in exc][:2] == ["named_contaminant", "contaminant_group"]
    assert all(r["published_share_percent"] for r in exc)
    pool = common.read_tsv(tmp_path / "data" / "cat" / "pool_red.tsv")
    assert {r["accession"] for r in pool} == members
    dup = next(r for r in pool if r["accession"] == "X00001")
    assert dup["dataset_rows"].count(";") == 1                # both rows of X00001 are recorded on its entry
    assert record["n_member"] == "2" and record["n_member_via_later_token"] == "1" and record["pool_size"] == "3"
    # two-pass lookup: the first pass asks about first tokens only, the second about the rest
    assert any("T00001" not in q for q in client.queries) and any("X00003" in q for q in client.queries)


def test_share_uses_mean_of_share_columns_and_log10_abundance():
    cols = {"share_columns": "P, Q", "abundance_column": ""}
    assert ep.published_share({"share_0": 1.0, "share_1": 3.0}, cols) == 2.0
    cols = {"share_columns": "", "abundance_column": "D", "abundance_scale": "log10"}
    assert ep.published_share({"abundance": 2}, cols) == 100.0
    assert ep.published_share({"abundance": ""}, cols) is None
