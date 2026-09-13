"""The standard's outputs must equal a fresh build from the files on disk; the EAA
transcription must be filled in; the code names no accession."""
import re

import pytest

from pipeline import common, aggregate as ag, uncertainty as un, stress as stt

needs_data = pytest.mark.skipif(not (common.STANDARD_DIR / "amino_acid_profiles.tsv").exists(),
                                reason="standard not built on this machine")
upstream_running = pytest.mark.skipif(common.running_command_is_before("standard"),
                                      reason="an earlier command is running; the standard is rebuilt by `python run.py standard`")


@needs_data
@upstream_running
def test_aggregate_outputs_are_current():
    assert ag.main(["--check"]) == 0


@needs_data
@upstream_running
def test_uncertainty_outputs_are_current():
    assert un.main(["--check"]) == 0


@needs_data
@upstream_running
def test_stress_outputs_are_current():
    if not (common.STANDARD_DIR / "stress_summary.ini").exists():
        pytest.skip("stress stage not run on this machine")
    assert stt.main(["--check"]) == 0


@needs_data
@upstream_running
def test_profiles_sum_to_one():
    for r in ag.read_tsv_skip_comments(common.STANDARD_DIR / "amino_acid_profiles.tsv"):
        s = sum(float(r[f"frac_{a}"]) for a in ag.AA)
        assert abs(s - 1.0) < 1e-9, (r["profile"], r["fiber_type"], r["convention"], s)


@upstream_running
def test_eaa_transcription_is_filled_in():
    """D62: the indispensable-amino-acid config is the author's transcription; while any `___`
    remains the EAA subset is not computed, and this test says so."""
    text = common.AGGREGATION_DECISIONS_INI.read_text(encoding="utf-8")
    body = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    left = [ln for ln in body if "___" in ln]
    assert left == [], f"config/aggregation_decisions.ini still has placeholders (D62): {left}"


def test_code_names_no_accession():
    acc = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
    for name in ("aggregate.py", "uncertainty.py", "stress.py", "plots.py"):
        text = (common.REPO_ROOT / "pipeline" / name).read_text(encoding="utf-8")
        hits = [ln for ln in text.splitlines() if acc.search(ln)]
        assert hits == [], (name, hits)
