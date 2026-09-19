"""The blood standard's tables must equal a fresh build (Step 7b).

Skipped until `python run.py blood-standard` has run on this machine, and while an earlier
command is running.
"""
import pytest

from pipeline import common, aggregate_pools

OUT = common.category_outputs("blood")["standard"]

needs_run = pytest.mark.skipif(
    not (OUT / "standard_summary.ini").exists() or common.running_command_is_before("blood-standard"),
    reason="blood-standard not run yet on this machine, or an earlier command is running")


@needs_run
def test_blood_standard_outputs_are_current():
    assert aggregate_pools.main(["--category", "blood", "--check"]) == 0


@needs_run
def test_blood_standard_sums_to_one_hundred():
    for name in ("_calculated_amino_acid_standard.tsv", "_calculated_amino_acid_standard_residue_convention.tsv"):
        rows = aggregate_pools.read_commented_tsv(OUT / name)
        for col in [c for c in rows[0] if c != "amino_acid"]:
            s = sum(float(r[col]) for r in rows)
            assert abs(s - 100.0) < 1e-6, (name, col, s)


@needs_run
def test_blood_plots_exist():
    for name in ("profile_pools_and_standard.png", "sensitivity_dominant_protein.png", "uncertainty_terms_standard.png"):
        assert (OUT / "plots" / name).exists(), name
