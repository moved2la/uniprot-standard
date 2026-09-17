"""The USDA outputs must equal a fresh read of the archives on disk."""
import pytest

from pipeline import common, usda

needs_data = pytest.mark.skipif(not (common.USDA_OUT_DIR / usda.MAIN_TABLE).exists(),
                                reason="the USDA archives have not been read on this machine")
upstream_running = pytest.mark.skipif(common.running_command_is_before("usda"),
                                      reason="an earlier command is running; the USDA tables are rebuilt by `python run.py usda`")


@needs_data
@upstream_running
def test_usda_outputs_are_current():
    try:
        assert usda.main(["--check"]) == 0
    except SystemExit as stop:                      # an archive named in config is not on disk right now
        pytest.skip(str(stop))


@needs_data
@upstream_running
def test_every_carried_food_has_at_least_one_amino_acid():
    """The main table is the foods WITH amino acids; the others are listed in their own file (U5)."""
    rows = usda.tsv_rows(common.USDA_OUT_DIR / usda.MAIN_TABLE)
    assert rows, "no foods carried"
    for r in rows:
        assert int(r["amino_acids_present"]) >= 1
        assert len(r["amino_acids_absent"]) == 20 - int(r["amino_acids_present"])


@needs_data
@upstream_running
def test_the_config_names_no_amino_acid():
    """Rule U3: the twenty are joined by the IUPAC trivial name, never named in config or code."""
    text = common.USDA_CONFIG_INI.read_text(encoding="utf-8").lower()
    body = "\n".join(ln for ln in text.splitlines() if not ln.strip().startswith("#"))
    for name in ("histidine", "leucine", "lysine", "valine", "threonine", "tryptophan",
                 "isoleucine", "methionine", "phenylalanine"):
        assert name not in body, f"{name} is named in config/usda_food_data.ini; the join is by trivial name (U3)"


@needs_data
@upstream_running
def test_the_manifest_records_every_archive_that_was_read():
    summary = common.read_ini(common.USDA_OUT_DIR / "usda_archive_summary.ini")
    manifest = common.read_ini(usda.MANIFEST)
    for sec in [s for s in summary.sections() if s.startswith("archive.")]:
        assert manifest.has_section(sec), sec
        assert manifest[sec]["sha256"] == summary[sec]["sha256"]
