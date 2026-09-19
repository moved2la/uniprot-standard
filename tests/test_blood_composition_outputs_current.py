"""The blood composition outputs must equal a fresh build from the data on disk (Step 7b).

Skipped until `python run.py blood-composition` has run on this machine, and while an earlier
command is running (blood-protein-set rebuilds the set this reads).
"""
import pytest

from pipeline import common, composition

FILES = common.category_files("blood")

needs_run = pytest.mark.skipif(
    not FILES["composition_tsv"].exists() or not FILES["accessions"].exists()
    or common.running_command_is_before("blood-composition"),
    reason="blood-composition not run yet on this machine, or an earlier command is running")


@needs_run
def test_blood_composition_outputs_are_current():
    assert composition.main(["--category", "blood", "--check"]) == 0


@needs_run
def test_blood_composition_covers_every_blood_accession():
    """Every entry of config/blood/accessions.ini without an open flag has a master row."""
    acc = common.read_ini(FILES["accessions"])
    rows = common.read_tsv(FILES["composition_tsv"])
    have = {r["accession"] for r in rows if r["segment_set"] == "master"}
    want = {a for a in acc.sections() if acc[a].get("flag_open") != "true"}
    assert want <= have


@needs_run
def test_blood_ptm_disclosure_present():
    for name in ("ptm_sites.tsv", "ptm_mass_deltas.tsv", "ptm_summary.ini"):
        assert (FILES["composition_dir"] / name).exists(), name
