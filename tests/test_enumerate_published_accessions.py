"""enumerate_published_pool on synthetic inputs: a workbook written by the test, a fake UniProt
lookup, and a decisions file with one pool. Accessions, genes and keyword ids are invented
(rules P1-P6, D107; keratin rule D109)."""
import logging

from openpyxl import Workbook

from pipeline import common, enumerate_pool as ep

KW = "KW-0416"


class FakeClient:
    release = "TEST"

    def __init__(self, entries, genes=None):
        self.entries = entries          # accession -> dict of the lookup fields (may carry "sec_acc")
        self.genes = genes or {}        # gene symbol -> [accessions whose primary symbol or synonym it is]
        self.queries = []

    def search_tsv(self, query, fields):
        self.queries.append(query)
        if query.startswith("(organism_id:9606) AND (reviewed:true)"):          # the M2 gene query (P4b)
            wanted = [t.split('"')[1] for t in query.split(" OR ") if 'gene_exact:"' in t]
            rows = []
            for a, e in self.entries.items():                                  # every entry whose primary symbols or synonyms name the gene
                names = {x.strip() for x in e.get("gene_primary", "").split(";")} | set(e.get("gene_synonym", "").split())
                if e.get("reviewed") == "reviewed" and e.get("organism_id") == "9606" and names & set(wanted):
                    rows.append(dict(e, accession=a))
            return rows
        terms = [t.strip("() ") for t in query.split(" OR ")]
        by_acc = [t.split(":", 1)[1] for t in terms if t.startswith("accession:")]
        by_sec = [t.split(":", 1)[1] for t in terms if t.startswith("sec_acc:")]
        rows = [dict(self.entries[a], accession=a) for a in by_acc if a in self.entries]   # includes inactive stubs
        for a, e in self.entries.items():                                      # the sec_acc query field
            if any(w in e.get("sec_acc", "").split(";") for w in by_sec) and e.get("reviewed"):
                rows.append(dict(e, accession=a))
        if "sec_acc" not in fields:
            rows = [{k: v for k, v in r.items() if k != "sec_acc"} for r in rows]
        return rows


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
        ("S00011;S00011-2", "eleven", "GK", 0.3, 0.3, 0.3, 0.3),                  # secondary of exactly one entry -> X00011
        ("S00012", "twelve", "GL;GM", 0.4, 0.4, 0.4, 0.4),                        # demerged: secondary of X00012 and X00013 -> shared
        ("T00014", "fourteen", "GN14", 0.2, 0.2, 0.2, 0.2),                       # TrEMBL only; gene resolves to X00014 (whose primary field lists two genes)
        ("T00015", "fifteen", "GN15", 0.1, 0.1, 0.1, 0.1),                        # TrEMBL only; gene ambiguous -> excluded
        ("T00016", "sixteen", "GP16;GP17", 0.1, 0.1, 0.1, 0.1),                   # two genes, each one entry -> shared row (P4b)
        ("T00018", "eighteen", "GS18", 0.1, 0.1, 0.1, 0.1),                       # two entries answer, one by primary symbol -> M3b
        ("T00020", "twenty", "GP16;GN15", 0.1, 0.1, 0.1, 0.1),                    # one gene ambiguous -> the row is excluded
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
        "X00011": dict(rev, sec_acc="S00011;S00099"),
        "S00012": {"reviewed": "", "organism_id": "", "gene_primary": "", "protein_name": "", "length": "", "keywordid": ""},  # an INACTIVE stub, as UniProt answers a demerged accession
        "X00012": dict(rev, sec_acc="S00012"), "X00013": dict(rev, sec_acc="S00012"),
        "X00014": dict(rev, gene_primary="GN14; GN14B"), "X00015": dict(rev, gene_primary="GN15"), "X00016": dict(rev, gene_primary="GP16", gene_synonym="GN15"),
        "X00017": dict(rev, gene_primary="GP17"), "X00021": dict(rev, gene_primary="GN15"),      # GN15: two primary-symbol entries -> ambiguous
        "X00018": dict(rev, gene_primary="GS18"), "X00019": dict(rev, gene_primary="OTHER", gene_synonym="GS18 ALSO"),
        "T00014": unrev, "T00015": unrev, "T00016": unrev, "T00018": unrev, "T00020": unrev,
        # X00009 is deliberately absent: UniProt returns nothing for it, as primary or secondary
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
    assert members == {"X00001", "X00002", "X00003", "X00011", "X00012", "X00013", "X00014", "X00016", "X00017", "X00018"}
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
    assert "gene:GH:unmapped" in rows["X00009"]["tokens_checked"]
    assert rows["M00010"]["outcome"] == "no_reviewed_human_entry" and "/10090" in rows["M00010"]["tokens_checked"]
    # P3b: a secondary accession resolves to its current entry; a demerged one to both, as a shared row
    assert rows["S00011;S00011-2"]["outcome"] == "member_via_secondary_accession" and rows["S00011;S00011-2"]["entry"] == "X00011"
    assert rows["S00012"]["outcome"] == "member_via_secondary_accession_shared" and rows["S00012"]["entry"] == "X00012;X00013"
    assert "1:S00012:inactive" in rows["S00012"]["tokens_checked"] and "secondary_of:X00012;X00013" in rows["S00012"]["tokens_checked"]
    # P4b: a group with no usable token resolves gene by gene; an entry's ';'-joined primary field counts each symbol
    assert rows["T00014"]["outcome"] == "member_via_gene" and rows["T00014"]["entry"] == "X00014" and rows["T00014"]["entry_token_rank"] == "gene"
    assert rows["T00015"]["outcome"] == "no_reviewed_human_entry" and "GN15:ambiguous:X00015;X00016;X00021" in rows["T00015"]["tokens_checked"]
    assert rows["T00016"]["outcome"] == "member_via_gene_shared" and rows["T00016"]["entry"] == "X00016;X00017"
    assert rows["T00018"]["outcome"] == "member_via_gene" and rows["T00018"]["entry"] == "X00018"          # M3b: primary over synonym
    assert "primary over synonyms X00019" in rows["T00018"]["tokens_checked"]
    assert rows["T00020"]["outcome"] == "no_reviewed_human_entry"                                          # one gene of two unresolved
    assert "" not in rows                                     # the empty row is skipped, not listed
    # the excluded listing is sorted by share, largest first, and carries the share
    exc = common.read_tsv(tmp_path / "data" / "cat" / "red_rows_excluded_with_share.tsv")
    assert [r["outcome"] for r in exc][:2] == ["named_contaminant", "contaminant_group"]
    assert all(r["published_share_percent"] for r in exc)
    pool = common.read_tsv(tmp_path / "data" / "cat" / "pool_red.tsv")
    assert {r["accession"] for r in pool} == members
    dup = next(r for r in pool if r["accession"] == "X00001")
    assert dup["dataset_rows"].count(";") == 1                # both rows of X00001 are recorded on its entry
    assert record["n_member"] == "2" and record["n_member_via_later_token"] == "1" and record["pool_size"] == "10"
    assert record["n_member_via_secondary_accession"] == "1" and record["n_member_via_secondary_accession_shared"] == "1"
    assert record["n_member_via_gene"] == "2" and record["n_member_via_gene_shared"] == "1"
    assert record["tokens_resolved_as_secondary"] == "2" and record["tokens_inactive_in_uniprot"] == "1"
    shared = [r for r in pool if r["accession"] in ("X00012", "X00013")]
    assert all(r["dataset_rows"] == shared[0]["dataset_rows"] for r in shared)   # both carry the one shared row
    # two-pass lookup: the first pass asks about first tokens only, the second about the rest
    assert any("T00001" not in q for q in client.queries) and any("X00003" in q for q in client.queries)


def test_share_uses_mean_of_share_columns_and_log10_abundance():
    cols = {"share_columns": "P, Q", "abundance_column": ""}
    assert ep.published_share({"share_0": 1.0, "share_1": 3.0}, cols) == 2.0
    cols = {"share_columns": "", "abundance_column": "D", "abundance_scale": "log10"}
    assert ep.published_share({"abundance": 2}, cols) == 100.0
    assert ep.published_share({"abundance": ""}, cols) is None
