"""The standard's outputs must equal a fresh build from the files on disk; the EAA
transcription must be filled in; the code names no accession."""
import re

import pytest

from pipeline import common, aggregate as ag, uncertainty as un, stress as stt, non_protein_metabolite_pools as npmp

needs_data = pytest.mark.skipif(not common.standard_path("amino_acid_profiles.tsv").exists(),
                                reason="standard not built on this machine")
upstream_running = pytest.mark.skipif(common.running_command_is_before("standard"),
                                      reason="an earlier command is running; the standard is rebuilt by `python run.py standard`")


@needs_data
@upstream_running
def test_aggregate_outputs_are_current():
    assert ag.main(["--check"]) == 0


@needs_data
@upstream_running
def test_metabolite_outputs_are_current():
    if not common.standard_path("non_protein_metabolite_pools_summary.ini").exists():
        pytest.skip("non_protein_metabolite_pools stage not run on this machine")
    assert npmp.main(["--check"]) == 0


@needs_data
@upstream_running
def test_uncertainty_outputs_are_current():
    assert un.main(["--check"]) == 0


@needs_data
@upstream_running
def test_stress_outputs_are_current():
    if not common.standard_path("stress_summary.ini").exists():
        pytest.skip("stress stage not run on this machine")
    assert stt.main(["--check"]) == 0


@needs_data
@upstream_running
def test_calculated_standard_columns_sum_to_100():
    rows = ag.read_tsv_skip_comments(common.standard_path("_calculated_amino_acid_standard.tsv"))
    body = [r for r in rows if r["amino_acid"] != "sum"]
    assert len(body) == 20
    for col in [c for c in rows[0] if c not in ("amino_acid", "three_letter", "name")]:
        if all(r[col] == "" for r in body):
            continue                                                    # the final column may be blank until the mix is filled
        assert abs(sum(float(r[col]) for r in body) - 100.0) < 5e-4, col


@needs_data
@upstream_running
def test_calculated_standard_with_non_protein_metabolite_pools_columns_sum_to_100():
    path = common.standard_path(npmp.TABLE)
    if not path.exists():
        pytest.skip("the adjusted standard is not written until config/non_protein_metabolite_pools.ini is filled")
    rows = ag.read_tsv_skip_comments(path)
    body = [r for r in rows if r["amino_acid"] != "sum"]
    assert len(body) == 20
    std = {r["amino_acid"]: r for r in ag.read_tsv_skip_comments(common.standard_path("_calculated_amino_acid_standard.tsv"))}
    for col in [c for c in rows[0] if c not in ("amino_acid", "three_letter", "name")]:
        if all(r[col] == "" for r in body):
            continue
        assert abs(sum(float(r[col]) for r in body) - 100.0) < 5e-4, col
        if col.startswith(("contractile_", "builders_")):
            assert all(r[col] == std[r["amino_acid"]][col] for r in body), f"{col} must be copied unchanged (D76)"


@needs_data
@upstream_running
def test_profiles_sum_to_one():
    for r in ag.read_tsv_skip_comments(common.standard_path("amino_acid_profiles.tsv")):
        s = sum(float(r[f"frac_{a}"]) for a in ag.AA)
        assert abs(s - 1.0) < 1e-9, (r["profile"], r["fiber_type"], r["convention"], s)


def _placeholders(path):
    body = [ln for ln in path.read_text(encoding="utf-8").splitlines() if not ln.lstrip().startswith("#")]
    return [ln for ln in body if "___" in ln]


@upstream_running
def test_eaa_transcription_is_filled_in():
    """D62: the indispensable-amino-acid file is the author's transcription; while any `___`
    remains the EAA subset is not computed, and this test says so."""
    left = _placeholders(common.EAA_INI)
    assert left == [], f"{common.EAA_INI.name} still has placeholders (D62): {left}"


@upstream_running
def test_fiber_type_mix_is_filled_in():
    """The final column of _calculated_amino_acid_standard.tsv needs the fiber-type shares from a
    cited source; while config/fiber_type_mix.ini holds `___` the column is blank, and this test says so."""
    left = _placeholders(common.FIBER_TYPE_MIX_INI)
    assert left == [], f"{common.FIBER_TYPE_MIX_INI.name} still has placeholders: {left}"


@upstream_running
def test_non_protein_metabolite_pools_is_filled_in():
    """D73–D77: the bound-pool file is the author's transcription; while any `___` remains the
    adjusted standard is not computed (or a sensitivity is skipped), and this test says so."""
    left = _placeholders(common.NON_PROTEIN_METABOLITE_POOLS_INI)
    assert left == [], f"{common.NON_PROTEIN_METABOLITE_POOLS_INI.name} still has placeholders: {left}"


def test_code_names_no_accession():
    acc = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
    for name in ("aggregate.py", "non_protein_metabolite_pools.py", "uncertainty.py", "stress.py", "plots.py"):
        text = (common.REPO_ROOT / "pipeline" / name).read_text(encoding="utf-8")
        hits = [ln for ln in text.splitlines() if acc.search(ln)]
        assert hits == [], (name, hits)
