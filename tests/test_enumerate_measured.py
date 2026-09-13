"""enumerate_measured_tier on synthetic inputs: a fake UniProt client, a fake dataset gene
list, and a fake contaminant list. Gene symbols and accessions are invented (rules M1-M5)."""
import logging

from pipeline import common, enumerate_pool as ep


class FakeClient:
    release = "TEST"

    def __init__(self, rows):
        self.rows = rows

    def search_tsv(self, query, fields):
        names = lambda r: [r["gene_primary"]] + r["gene_synonym"].split()
        return [r for r in self.rows if any(f'gene_exact:"{n}"' in query for n in names(r))]


def test_measured_tier_outcomes(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "DATA_DIR", tmp_path)
    monkeypatch.setattr(ep, "dataset_gene_cells", lambda src, key, log: ["GA", "GB;GC", "GD", "GE", "GF", "GOLD"])
    monkeypatch.setattr(ep, "contaminant_accessions", lambda src, log: {"X00005"})
    rows = [
        {"accession": "X00001", "gene_primary": "GA", "gene_synonym": "", "protein_name": "a", "length": "10", "reviewed": "reviewed", "organism_id": "9606"},
        {"accession": "X00002", "gene_primary": "GB", "gene_synonym": "", "protein_name": "b", "length": "10", "reviewed": "reviewed", "organism_id": "9606"},
        {"accession": "X00003", "gene_primary": "GC", "gene_synonym": "", "protein_name": "c", "length": "10", "reviewed": "reviewed", "organism_id": "9606"},
        {"accession": "X00004", "gene_primary": "GC", "gene_synonym": "", "protein_name": "c2", "length": "10", "reviewed": "reviewed", "organism_id": "9606"},  # GC ambiguous
        {"accession": "X00005", "gene_primary": "GE", "gene_synonym": "", "protein_name": "e", "length": "10", "reviewed": "reviewed", "organism_id": "9606"},  # contaminant
        {"accession": "X00006", "gene_primary": "GNEW", "gene_synonym": "GOLD", "protein_name": "f", "length": "10", "reviewed": "reviewed", "organism_id": "9606"},  # synonym match
    ]
    sec = {"dataset_source": "s", "dataset_file": "file.1", "contaminant_source": "c"}
    pool, counts = ep.enumerate_measured_tier(FakeClient(rows), "2", sec, {"1": {"X00002"}}, logging.getLogger("t"))
    assert pool == {"X00001", "X00006"}                      # GA member; GOLD -> X00006 via synonym
    assert counts["unmapped"] == 2 and counts["ambiguous"] == 1 and counts["contaminant"] == 1 and counts["in_tier1"] == 1
    mapping = {r["gene"]: r for r in common.read_tsv(tmp_path / "tier2_gene_mapping.tsv")}
    assert mapping["GB"]["outcome"] == "in_tier1" and mapping["GC"]["outcome"] == "ambiguous"
    assert mapping["GD"]["outcome"] == "unmapped" and mapping["GE"]["outcome"] == "contaminant"
    assert mapping["GOLD"]["matched_by"] == "synonym" and mapping["GOLD"]["uniprot_primary_symbol"] == "GNEW"
    pool_rows = common.read_tsv(tmp_path / "tier2_pool.tsv")
    assert {r["accession"] for r in pool_rows} == {"X00001", "X00006"}
    assert all(r["subtree_terms_returning_entry"] == "" for r in pool_rows)


def test_contaminant_header_parsing(tmp_path, monkeypatch):
    import configparser, gzip
    lit = tmp_path / "literature"; (lit / "c").mkdir(parents=True)
    fasta = b">CON__X00005 something\nMKV\n>sp|X00007|NAME_HUMAN desc\nMKV\n>X00008-2 iso\nMKV\n"
    p = lit / "c" / "x.fasta.gz"; p.write_bytes(gzip.compress(fasta))
    man = configparser.ConfigParser(interpolation=None); man.optionxform = str
    man.add_section("file.c.file.1"); man["file.c.file.1"]["path"] = "literature/c/x.fasta.gz"
    man["file.c.file.1"]["sha256"] = common.sha256_file(p)
    monkeypatch.setattr(common, "DATA_DIR", tmp_path)
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    with open(lit / "manifest.ini", "w") as fh:
        man.write(fh)
    accs = ep.contaminant_accessions("c", logging.getLogger("t"))
    assert accs == {"X00005", "X00007", "X00008"}
