"""The match outputs must equal a fresh read of the references and the food tables on disk."""
import pytest

from pipeline import common, match

needs_data = pytest.mark.skipif(not (common.MATCH_OUT_DIR / "match_rate_per_food.tsv").exists(),
                                reason="Match Rate has not been run on this machine")
upstream_running = pytest.mark.skipif(common.running_command_is_before("match"),
                                      reason="an earlier command is running; the match tables are rebuilt by `python run.py match`")


@needs_data
@upstream_running
def test_match_outputs_are_current():
    assert match.main(["--check"]) == 0


@needs_data
@upstream_running
def test_every_step_row_has_a_score_between_zero_and_one_hundred_and_a_limiting_amino_acid():
    text = (common.MATCH_OUT_DIR / "match_rate_steps.tsv").read_text(encoding="utf-8")
    body = [ln.split("\t") for ln in text.splitlines() if not ln.startswith("#")]
    rows = [dict(zip(body[0], r)) for r in body[1:]]
    assert rows, "no food was scored"
    for r in rows:
        assert 0.0 <= float(r["match_rate_percent"]) <= 100.0 + 1e-9
        assert r["limiting_amino_acid"] and r["reference"] and r["source"]
        if r["percent_wasted"]:                                                              # blank on a zero score
            assert abs(float(r["percent_utilized"]) + float(r["percent_wasted"]) - 1.0) < 1e-3   # rows 155 + 156 = 1


@needs_data
@upstream_running
def test_the_scored_sets_resolve_to_the_expected_sizes():
    """The FAO-derived set is nine letters (D86); the original reference's own essential set is eight (no tryptophan)."""
    cp = common.read_ini(common.MATCH_OUT_DIR / "match_summary.ini")
    sizes = {s[len("scored_set."):]: len(cp[s]["amino_acids"].strip()) for s in cp.sections() if s.startswith("scored_set.")}
    assert sizes.get("fao_2013_indispensable") == 9, sizes
    assert sizes.get("gorissen_2018_essential") == 8, sizes
