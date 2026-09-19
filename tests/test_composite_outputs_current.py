"""The composite outputs must equal a fresh mix of the category standards on disk (Layer C, D121).

Skipped until `python run.py composite` has run on this machine, and while an earlier command is
running (the category standards it reads are rebuilt first).
"""
import pytest

from pipeline import common, composite

OUT = common.COMPOSITE_OUT_DIR

needs_run = pytest.mark.skipif(
    not (OUT / "composite_summary.ini").exists() or common.running_command_is_before("composite"),
    reason="composite not run yet on this machine, or an earlier command is running")


@needs_run
def test_composite_outputs_are_current():
    assert composite.main(["--check"]) == 0


@needs_run
def test_composite_sums_to_one_hundred_and_names_every_category():
    rows = composite.read_commented_tsv(OUT / "_calculated_amino_acid_standard_composite.tsv")
    cp = common.read_ini(common.TISSUE_MASS_FRACTIONS_INI)
    names = [s[len("category."):] for s in cp.sections() if s.startswith("category.")]
    for col in names + ["composite"]:
        assert col in rows[0], col
        s = sum(float(r[col]) for r in rows)
        assert abs(s - 100.0) < 1e-2, (col, s)


@needs_run
def test_the_shares_are_the_masses_over_their_sum():
    cp = common.read_ini(OUT / "composite_summary.ini")
    m = cp["masses"]
    names = [k[:-len("_protein_g")] for k in m if k.endswith("_protein_g") and not k.startswith("total")]
    total = sum(float(m[f"{n}_protein_g"]) for n in names)
    assert abs(total - float(m["total_protein_g"])) < 1e-6
    for n in names:
        assert abs(float(m[f"{n}_share"]) - float(m[f"{n}_protein_g"]) / total) < 1e-9
