"""fetch_amino_acid_masses and composition.load_masses, offline, on synthetic fixtures.

What these tests check, in plain terms:
- the mass table loader refuses a file that is missing any of the twenty amino
  acids, water, or hydrogen, and for a complete file every residue mass equals
  free mass minus water;
- the PubChem response parser accepts exactly one compound and rejects zero or
  several (a name that resolves to several compounds must stop the run, not pick one);
- the symbol-table parser pulls "<trivial name> <three-letter> <one-letter>" rows
  out of tag-stripped page text, keeps only the standard letters, keeps two-word
  names ("... acid"), and refuses a letter that appears twice;
- the hand-written decisions file holds only URLs and the rule -- no name, symbol,
  or mass.
"""
import json
import pytest
from conftest import FIXTURES

from pipeline import common
from pipeline import composition
from pipeline import fetch_amino_acid_masses as fam


def test_load_masses_complete_and_consistent():
    m = composition.load_masses(FIXTURES / "masses_synthetic.ini")
    assert sorted(m["free"]) == sorted(common.AMINO_ACIDS)
    assert m["water"] == 10.0 and m["hydrogen"] == 2.0
    for a in common.AMINO_ACIDS:
        assert m["residue"][a] == pytest.approx(m["free"][a] - m["water"])


def test_load_masses_rejects_incomplete(tmp_path):
    text = (FIXTURES / "masses_synthetic.ini").read_text()
    broken = text.replace("[W]\n", "[W_missing]\n")
    p = tmp_path / "m.ini"
    p.write_text(broken)
    with pytest.raises(ValueError):
        composition.load_masses(p)
    p.write_text(text.replace("[hydrogen]\n", "[h]\n"))
    with pytest.raises(ValueError):
        composition.load_masses(p)


def test_pubchem_parser_exactly_one_compound():
    one = json.loads((FIXTURES / "pubchem_one.json").read_text())
    row = fam.parse_pubchem_properties(one)
    assert row["CID"] == 999001 and row["MolecularWeight"] == "123.45"
    two = json.loads((FIXTURES / "pubchem_two.json").read_text())
    with pytest.raises(ValueError):
        fam.parse_pubchem_properties(two)
    with pytest.raises(ValueError):
        fam.parse_pubchem_properties({"PropertyTable": {"Properties": []}})
    with pytest.raises(ValueError):
        fam.parse_pubchem_properties({"Fault": {"Message": "no such compound"}})


def test_symbol_table_parser_on_synthetic_page():
    page = (FIXTURES / "iupac_table_synthetic.html").read_text()
    text = fam.page_text(page)
    assert "<" not in text and "&nbsp;" not in text and "  " not in text
    found = fam.parse_symbol_table(text, letters="AB")
    assert found == {"A": {"trivial_name": "Fixturine", "three_letter": "Fix"},
                     "B": {"trivial_name": "Fixturic acid", "three_letter": "Fic"}}
    # a symbol outside the requested letters ('Unspecified Xaa X') is ignored
    assert "X" not in fam.parse_symbol_table(text, letters="AB")


def test_symbol_table_parser_refuses_duplicate_letter():
    with pytest.raises(ValueError):
        fam.parse_symbol_table("Fixturine Fix A 1-fixture Fixturone Fox A 2-fixture ", letters="A")


def test_decisions_file_holds_only_urls_and_rule():
    dec = common.read_ini(common.COMPOSITION_DECISIONS_INI)
    assert set(dec.sections()) == {"amino_acid_symbols", "mass_source", "query_name_rule", "extra_compounds"}
    assert dec["amino_acid_symbols"]["url"].startswith("https://")
    assert dec["mass_source"]["base_url"].startswith("https://")
    assert dec["query_name_rule"]["chiral_prefix"]
    for s in dec.sections():
        assert dec[s]["decision"].startswith("D")
        for k, v in dec[s].items():
            assert not any(ch.isdigit() for ch in v) or k in ("url", "base_url", "decision"), (s, k, v)
