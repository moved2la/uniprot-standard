"""The literature stage hashes what is on disk and never downloads (D80).

Two kinds of test here. The synthetic ones build a small config and a small data/literature
tree in a temp folder and check the rules. The currency ones read the repository's own
config and assert the two config files still agree with each other.
"""
import configparser
import hashlib

import pytest

from pipeline import common, fetch_literature as fl


class _Log:
    def __init__(self): self.lines = []
    def _add(self, m, *a): self.lines.append(m % a if a else m)
    info = warning = error = _add


PDF = b"%PDF-1.4\nsynthetic\n"
HTML = b"<!DOCTYPE html><html><head><title>Access denied</title></head></html>"


def _sha(b): return hashlib.sha256(b).hexdigest()


@pytest.fixture
def lit_repo(tmp_path, monkeypatch):
    """A config with four sources: pinned, unpinned, cited-only, and one file missing."""
    cfg = tmp_path / "config" / "literature_sources.ini"
    lit = tmp_path / "data" / "literature"
    cfg.parent.mkdir(parents=True)
    (lit / "alpha_2001").mkdir(parents=True)
    (lit / "alpha_2001" / "paper.pdf").write_bytes(PDF)
    (lit / "beta_2002").mkdir(parents=True)
    (lit / "beta_2002" / "table.pdf").write_bytes(PDF + b"beta")
    cfg.write_text(f"""[meta]
step = test

[categories]
columns = original, skeletal_muscle, blood, other
unused_marker = unused

[source.alpha_2001]
id_short = A1
citation = Alpha A. A paper. Journal 2001.
role = primary
role_decision = D1
used_for = skeletal_muscle
file.1.published_name = the paper
file.1.url = https://example.org/files/paper.pdf
file.1.sha256 = {_sha(PDF)}

[source.beta_2002]
id_short = B1
citation = Beta B. Another. Journal 2002.
role = cross_check
role_decision = D2
used_for = skeletal_muscle, blood
file.1.published_name = a table
file.1.url = https://example.org/files/table.pdf
file.1.sha256 =

[source.gamma_2003]
id_short = C1
citation = Gamma G. A method. Journal 2003.
role = method_citation
role_decision = D3
used_for = blood

[source.delta_2004]
id_short = D1
citation = Delta D. Missing. Journal 2004.
role = cross_check
role_decision = D4
used_for = skeletal_muscle
file.1.published_name = never saved
file.1.url = https://example.org/files/absent.pdf
file.1.sha256 =
""", encoding="utf-8")
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(fl, "ROOT", tmp_path)
    monkeypatch.setattr(fl, "CONFIG", cfg)
    monkeypatch.setattr(fl, "LIT_DIR", lit)
    monkeypatch.setattr(fl, "MANIFEST_INI", lit / "manifest.ini")
    monkeypatch.setattr(fl, "MANIFEST_MAP_TSV", lit / "manifest_map.tsv")
    return {"tmp": tmp_path, "cfg": cfg, "lit": lit}


def test_a_pinned_file_on_disk_is_hashed_and_lands_in_the_manifest(lit_repo):
    _, manifest, flags, pending, n_files, n_cited = fl.build(_Log())
    assert manifest["file.alpha_2001.file.1"]["sha256"] == _sha(PDF)
    assert manifest["file.alpha_2001.file.1"]["bytes"] == str(len(PDF))
    assert n_files == 2


def test_a_source_with_no_files_is_cited_only_and_never_flagged(lit_repo):
    _, manifest, flags, _, _, n_cited = fl.build(_Log())
    assert n_cited == 1
    assert not any("gamma_2003" in f for f in flags)
    assert not any("gamma" in s for s in manifest.sections())


def test_a_missing_file_is_a_flag(lit_repo):
    _, _, flags, _, _, _ = fl.build(_Log())
    assert [f for f in flags if f.startswith("delta_2004.file.1") and "not on disk" in f]


def test_a_blank_pin_is_recorded_not_guessed(lit_repo):
    _, _, _, pending, _, _ = fl.build(_Log())
    assert ("beta_2002", 1, _sha(PDF + b"beta")) in pending
    fl.record_hashes(pending)
    cp = common.read_ini(lit_repo["cfg"])
    assert cp["source.beta_2002"]["file.1.sha256"] == _sha(PDF + b"beta")
    # the pinned source is untouched
    assert cp["source.alpha_2001"]["file.1.sha256"] == _sha(PDF)


def test_a_changed_file_fails_its_pin(lit_repo):
    (lit_repo["lit"] / "alpha_2001" / "paper.pdf").write_bytes(PDF + b"tampered")
    _, _, flags, _, _, _ = fl.build(_Log())
    assert [f for f in flags if "alpha_2001.file.1" in f and "sha256 mismatch" in f]


def test_a_saved_paywall_page_is_not_a_pdf(lit_repo):
    (lit_repo["lit"] / "alpha_2001" / "paper.pdf").write_bytes(HTML)
    _, _, flags, _, _, _ = fl.build(_Log())
    assert [f for f in flags if "alpha_2001.file.1" in f and "HTML" in f]


def test_a_placeholder_url_is_a_flag(lit_repo):
    text = lit_repo["cfg"].read_text(encoding="utf-8").replace(
        "file.1.url = https://example.org/files/absent.pdf", "file.1.url = ___")
    lit_repo["cfg"].write_text(text, encoding="utf-8")
    _, _, flags, _, _, _ = fl.build(_Log())
    assert [f for f in flags if "delta_2004.file.1" in f and "placeholder" in f.lower() or "___" in f]


def test_the_map_columns_come_from_config_in_config_order(lit_repo):
    """Not first-appearance order: a person laid the columns out and the map keeps that."""
    cp = common.read_ini(lit_repo["cfg"])
    header, body, unknown = fl.manifest_map_rows(cp)
    assert header == ["source", "original", "skeletal_muscle", "blood", "other"]
    assert unknown == []
    rows = {r[0]: r[1:] for r in body}
    assert len(rows) == 4
    assert rows["alpha_2001"] == ["", "X", "", ""]
    assert rows["beta_2002"] == ["", "X", "X", ""]                # a source can serve several
    assert rows["gamma_2003"] == ["", "", "X", ""]                # cited-only still appears


def test_a_source_that_feeds_no_calculation_says_so(lit_repo):
    """A row of blanks reads as an unfilled line; the marker says it was decided."""
    text = lit_repo["cfg"].read_text(encoding="utf-8").replace(
        "role_decision = D4\nused_for = skeletal_muscle", "role_decision = D4\nused_for = unused")
    lit_repo["cfg"].write_text(text, encoding="utf-8")
    cp = common.read_ini(lit_repo["cfg"])
    header, body, unknown = fl.manifest_map_rows(cp)
    rows = {r[0]: r[1:] for r in body}
    assert rows["delta_2004"] == ["", "", "", "unused"]            # the marker sits in `other`
    assert unknown == []


def test_a_used_for_naming_a_column_that_does_not_exist_is_a_flag(lit_repo):
    text = lit_repo["cfg"].read_text(encoding="utf-8").replace(
        "used_for = skeletal_muscle, blood", "used_for = skeletal_muscle, liver")
    lit_repo["cfg"].write_text(text, encoding="utf-8")
    cp = common.read_ini(lit_repo["cfg"])
    _, _, unknown = fl.manifest_map_rows(cp)
    assert [u for u in unknown if "liver" in u]


def test_a_source_without_used_for_is_a_flag(lit_repo):
    text = lit_repo["cfg"].read_text(encoding="utf-8").replace("used_for = blood\n", "", 1)
    lit_repo["cfg"].write_text(text, encoding="utf-8")
    _, _, flags, _, _, _ = fl.build(_Log())
    assert [f for f in flags if "no used_for" in f]


def test_a_file_without_an_extension_is_found_by_its_magic_bytes(lit_repo):
    """F7b: the URL has no extension, the browser saved it with one."""
    (lit_repo["lit"] / "eps_2005").mkdir()
    # the URL's last segment is "abcd" with no extension; the browser saved it as abcd.pdf,
    # exactly as harris_1998's tandfonline abstract is stored as 026404198366443.pdf
    (lit_repo["lit"] / "eps_2005" / "abcd.pdf").write_bytes(PDF)
    text = lit_repo["cfg"].read_text(encoding="utf-8") + """
[source.eps_2005]
id_short = E1
citation = Eps E. No extension. Journal 2005.
role = cross_check
role_decision = D5
used_for = skeletal_muscle
file.1.published_name = article
file.1.url = https://example.org/doi/10.1234/abcd
file.1.sha256 =
"""
    lit_repo["cfg"].write_text(text, encoding="utf-8")
    _, manifest, flags, _, _, _ = fl.build(_Log())
    assert not [f for f in flags if "eps_2005" in f]
    assert manifest["file.eps_2005.file.1"]["path"].endswith("abcd.pdf")
    assert manifest["file.eps_2005.file.1"]["served_name"] == "abcd"


# --------------------------------------------------------------------------- the real config

def test_every_role_named_in_mass_fraction_decisions_resolves_to_a_source():
    """The Step 7a split: the column maps stayed, the sources left. The two must still agree."""
    dec = common.read_ini(common.CONFIG_DIR / "mass_fraction_decisions.ini")
    src = common.read_ini(common.LITERATURE_SOURCES_INI)
    known = {s[len("source."):] for s in src.sections() if s.startswith("source.")}
    assert known, "no [source.*] sections in literature_sources.ini"
    assert not [s for s in dec.sections() if s.startswith("source.")], \
        "sources belong in literature_sources.ini, not mass_fraction_decisions.ini"
    for name in dec.sections():
        if not name.startswith("columns."):
            continue
        sid = name[len("columns."):].split(".file.")[0]
        assert sid in known, f"[{name}] maps columns for a source not in literature_sources.ini"


def test_every_source_declares_what_it_is_used_for():
    src = common.read_ini(common.LITERATURE_SOURCES_INI)
    missing = [s for s in src.sections()
               if s.startswith("source.") and not src[s].get("used_for", "").strip()]
    assert not missing, f"sources without used_for: {missing}"


def test_the_real_map_places_every_source_in_a_declared_column():
    src = common.read_ini(common.LITERATURE_SOURCES_INI)
    header, body, unknown = fl.manifest_map_rows(src)
    assert not unknown, unknown
    assert header == ["source", "original", "skeletal_muscle", "non_protein_metabolites", "other"]
    blank = [r[0] for r in body if not any(c for c in r[1:])]
    assert not blank, f"sources with no column and no unused marker: {blank}"


def test_no_placeholder_hash_anywhere_in_config():
    """F2b: an unknown hash is blank. `___` is only ever a value a person types."""
    bad = []
    for path in sorted(common.CONFIG_DIR.rglob("*.ini")):
        cp = common.read_ini(path)
        for sec in cp.sections():
            for key, value in cp[sec].items():
                if key.endswith("sha256") and value.strip() == "___":
                    bad.append(f"{path.name} [{sec}] {key}")
    assert not bad, f"placeholder in a sha256 key: {bad}"
