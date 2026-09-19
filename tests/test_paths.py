"""Every generated path must follow a redirected root.

This is the rule two bugs were paid for. A path computed at import time from OUTPUTS_DIR or
DATA_DIR cannot follow a caller that moves the root — a test writing to a temp folder, or a
category writing under outputs/<category>/. Both times the failure landed far from the cause:
a stage silently wrote into the real repository while a test had built its own tree.

These tests fail the moment a path goes back to being frozen at import.
"""
import pytest

from pipeline import common


OUTPUT_PATHS = [
    "intermediate_dir", "protein_set_out_dir", "composition_dir", "digest_dir",
    "literature_inventory_dir", "mass_fractions_dir", "standard_dir",
    "composition_tsv", "composition_summary_ini",
]
DATA_PATHS = ["literature_dir", "literature_manifest", "literature_manifest_map"]


@pytest.mark.parametrize("name", OUTPUT_PATHS)
def test_an_outputs_path_follows_a_redirected_outputs_dir(name, tmp_path, monkeypatch):
    monkeypatch.setattr(common, "OUTPUTS_DIR", tmp_path / "elsewhere")
    got = getattr(common, name)()
    assert got.is_relative_to(tmp_path), f"common.{name}() ignored the redirect: {got}"


@pytest.mark.parametrize("name", DATA_PATHS)
def test_a_data_path_follows_a_redirected_data_dir(name, tmp_path, monkeypatch):
    monkeypatch.setattr(common, "DATA_DIR", tmp_path / "elsewhere")
    got = getattr(common, name)()
    assert got.is_relative_to(tmp_path), f"common.{name}() ignored the redirect: {got}"


def test_standard_path_follows_a_redirected_outputs_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "OUTPUTS_DIR", tmp_path / "elsewhere")
    for name in ("_calculated_amino_acid_standard.tsv", "uncertainty_intervals.tsv",
                 "stress_shifts.tsv", "completeness_sensitivity.tsv"):
        assert common.standard_path(name).is_relative_to(tmp_path), name
    for name, ext in (("profile_I", "png"), ("uncertainty_terms_I", "svg")):
        assert common.standard_plot_path(name, ext).is_relative_to(tmp_path), name


def test_a_category_gets_the_same_tree_under_its_own_root():
    """What Step 7b needs: blood's standard laid out exactly like muscle's."""
    blood = common.category_outputs("blood")
    assert blood["standard"] == common.OUTPUTS_DIR / "blood" / "standard"
    assert blood["intermediate"] == common.OUTPUTS_DIR / "blood" / "intermediate"
    got = common.standard_path("uncertainty_intervals.tsv", blood["standard"])
    assert got == common.OUTPUTS_DIR / "blood" / "standard" / "uncertainty" / "uncertainty_intervals.tsv"
    got = common.standard_plot_path("profile_I", "png", blood["standard"])
    assert got == common.OUTPUTS_DIR / "blood" / "standard" / "plots" / "profile_I.png"


def test_the_placement_map_holds_subfolder_names_not_absolute_paths():
    """The reason standard_path can take a base at all."""
    for name, sub in common.STANDARD_SUBFOLDERS.items():
        assert isinstance(sub, str), f"{name} maps to {sub!r}, not a subfolder name"
    for name, sub in common.STANDARD_PLOT_SUBFOLDERS.items():
        assert isinstance(sub, str), f"{name} maps to {sub!r}, not a subfolder name"


def test_category_files_follow_a_redirected_root(tmp_path, monkeypatch):
    """Step 7b: a category's decisions, config, pools, flags and tables all sit under the category."""
    monkeypatch.setattr(common, "CONFIG_DIR", tmp_path / "config")
    monkeypatch.setattr(common, "DATA_DIR", tmp_path / "data")
    monkeypatch.setattr(common, "OUTPUTS_DIR", tmp_path / "outputs")
    f = common.category_files("blood")
    assert f["decisions"] == tmp_path / "config" / "blood" / "protein_set_decisions.ini"
    assert f["accessions"] == tmp_path / "config" / "blood" / "accessions.ini"
    assert f["pool_queries"] == tmp_path / "data" / "blood" / "pool_queries.ini"
    assert f["flags"] == tmp_path / "outputs" / "blood" / "flags.tsv"
    assert f["protein_set_out"] == tmp_path / "outputs" / "blood" / "intermediate" / "protein_set"
    assert common.pool_file("blood", "plasma") == tmp_path / "data" / "blood" / "pool_plasma.tsv"
    assert common.pool_file(None, "2") == tmp_path / "data" / "tier2_pool.tsv"


def test_muscle_files_are_the_flat_constants():
    f = common.category_files(None)
    assert f["decisions"] == common.DECISIONS_INI and f["accessions"] == common.ACCESSIONS_INI
    assert f["flags"] == common.FLAGS_TSV and f["pool_dir"] == common.DATA_DIR


def test_raw_entry_is_found_in_the_root_or_any_category_folder(tmp_path, monkeypatch):
    """D110: the raw folder is nested per category; readers look everywhere, nothing is deleted."""
    raw = tmp_path / "uniprot_raw"
    monkeypatch.setattr(common, "UNIPROT_RAW_DIR", raw)
    (raw / "skeletal_muscle").mkdir(parents=True); (raw / "blood").mkdir()
    (raw / "P00001.json").write_text("{}"); (raw / "skeletal_muscle" / "P00002.json").write_text("{}")
    (raw / "blood" / "P00003.json").write_text("{}")
    assert common.uniprot_raw_entry("P00001") == raw / "P00001.json"
    assert common.uniprot_raw_entry("P00002") == raw / "skeletal_muscle" / "P00002.json"
    assert common.uniprot_raw_entry("P00003") == raw / "blood" / "P00003.json"
    assert common.uniprot_raw_entry("P00004") is None
    assert common.uniprot_raw_dir("blood") == raw / "blood"
    assert common.uniprot_raw_dir(None) == raw / common.MUSCLE_CATEGORY
