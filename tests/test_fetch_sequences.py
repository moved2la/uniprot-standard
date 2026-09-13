"""fetch_sequences: JSON extraction and sequence verification, offline."""
import json
import pytest
from conftest import FIXTURES

from pipeline import fetch_sequences as va
from pipeline import common


@pytest.fixture
def entry():
    return json.loads((FIXTURES / "entry_one_chain.json").read_text())


def test_extract_fields(entry):
    rec = va.extract(entry)
    assert rec["accession"] == "X99001" and rec["gene"] == "FIXA"
    assert rec["length"] == "13" and rec["seq_version"] == "2"
    assert [f["type"] for f in rec["features"]] == ["Initiator methionine", "Chain", "Alternative sequence"]
    assert rec["n_features_total"] == "4"
    assert [i["status"] for i in rec["isoforms"]] == ["Displayed", "Described"]
    assert rec["isoforms"][1]["note"].startswith("A free-text note")
    assert [a["term_id"] for a in rec["annotations"]] == ["TEST:0000003", "TEST:0000006"]
    assert rec["annotations"][0]["evidence"] == "IDA:UniProtKB"
    assert rec["tissue_specificity"] == "Recorded, never parsed."


def test_md5_verification_passes_and_fails(entry):
    rec = va.extract(entry)
    assert va.check_sequence(rec, rec["sequence"]) == []
    rec2 = va.extract(entry)
    problems = va.check_sequence(rec2, rec2["sequence"][:-1] + "A")
    assert any("MD5" in p for p in problems) and any("differs" in p for p in problems)


def test_nonstandard_letters_are_recorded_not_failures(entry):
    # R5 (D59): the letters are recorded on the record; build_protein_set excludes and lists the entry
    rec = va.extract(entry)
    rec["sequence"] = "MGGGAAAKKKWWXU"
    rec["md5_uniprot"] = common.md5_text(rec["sequence"])
    problems = va.check_sequence(rec, rec["sequence"])
    assert problems == [] and rec["non_standard_letters"] == "UX"


def test_md5_matches_uniprot_convention():
    # UniProt publishes the uppercase hex MD5 of the raw sequence string
    assert common.md5_text("MGGGAAAKKKWWY") == "B8BAB4B84AC5F378A79AC603BD2410F4"


def test_section_roundtrip(entry, tmp_path):
    rec = va.extract(entry)
    va.check_sequence(rec, rec["sequence"])
    rec["canonical_isoform_id"] = "X99001-1"
    cp = common.new_ini()
    cp["X99001"] = va.to_section(rec, ["1"])
    p = tmp_path / "v.ini"
    common.write_ini(cp, p, ["test"])
    back = common.read_ini(p)["X99001"]
    assert back["md5_match"] == "true" and back["feature.1.type"] == "Chain"
    assert back["annotation.1.evidence"] == "IEA:Ensembl"
