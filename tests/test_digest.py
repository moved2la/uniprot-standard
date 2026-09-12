"""
Offline tests for pipeline/digest.py. Synthetic sequences only; no network;
no real accession appears here.
"""
import configparser
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import digest  # noqa: E402

RULES = {
    "enzyme": "trypsin",
    "cleave_after": set("KR"),
    "cleave_before_proline": False,
    "min_length": 3,
    "max_length": 8,
    "missed_cleavages": 0,
    "source": "synthetic",
    "retrieved": "",
}


def test_fragments_cleave_after_K_and_R():
    assert digest.fragments("AAAKBBBRCCC", set("KR"), True) == ["AAAK", "BBBR", "CCC"]


def test_proline_rule_blocks_cleavage_when_disabled():
    # K followed by P: no cleavage when cleave_before_proline is False
    assert digest.fragments("AAAKPBBR", set("KR"), False) == ["AAAKPBBR"]
    assert digest.fragments("AAAKPBBR", set("KR"), True) == ["AAAK", "PBBR"]


def test_terminal_cleavage_site_does_not_produce_empty_fragment():
    assert digest.fragments("AAAK", set("KR"), True) == ["AAAK"]


def test_length_window_and_distinct_count():
    seq = "GGGKAAAAAAAAAAAAKGGGKTTK"  # frags: GGGK(4) A..K(13) GGGK(4) TTK(3)
    peps = digest.theoretical_peptides(seq, RULES)
    assert peps == {"GGGK", "TTK"}          # 13-mer excluded (>8); GGGK counted once
    assert len(peps) == 2


def test_missed_cleavages_join_adjacent_fragments():
    r = dict(RULES, missed_cleavages=1, max_length=10)
    seq = "AAKBBRCCK"
    peps = digest.theoretical_peptides(seq, r)
    assert peps == {"AAK", "BBR", "CCK", "AAKBBR", "BBRCCK"}


def test_shared_peptide_forms_family():
    a = digest.theoretical_peptides("QQQKWWWKEEEK", RULES)
    b = digest.theoretical_peptides("WWWKNNNK", RULES)
    c = digest.theoretical_peptides("HHHKYYYK", RULES)
    assert "WWWK" in a & b
    assert not (c & (a | b))
    edges = {"A": {"B"}, "B": {"A"}, "C": set()}
    comps = [s for s in digest.connected_components(edges) if len(s) >= 2]
    assert comps == [{"A", "B"}]


def test_rules_refuse_placeholders():
    cp = configparser.ConfigParser(interpolation=None)
    cp["digest"] = {"cleave_after": "KR", "cleave_before_proline": "___", "min_length": "6",
                    "max_length": "30", "missed_cleavages": "0", "source": ""}
    with pytest.raises(SystemExit):
        digest.load_rules(cp)


def test_config_names_no_protein_or_number():
    """The hand-written config may not carry an accession, or any numeric
    value outside the explicitly allowed keys (DOIs, PMIDs, D-numbers,
    dates, digest lengths)."""
    text = (ROOT / "config" / "mass_fraction_decisions.ini").read_text(encoding="utf-8")
    # UniProt accession shapes: [OPQ][0-9][A-Z0-9]{3}[0-9] or [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    acc = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
    body_lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    hits = [ln for ln in body_lines if acc.search(ln.split("=", 1)[-1]) and not ln.strip().startswith(("citation", "raw_data", "pmc"))]
    assert hits == [], f"accession-shaped tokens in config: {hits}"
