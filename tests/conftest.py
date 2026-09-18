"""pytest configuration: repository root on sys.path; every test run logged to logs/."""
import datetime as dt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
FIXTURES = Path(__file__).parent / "fixtures"


def pytest_configure(config):
    logs = ROOT / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    config.option.log_file = str(logs / f"pytest_{stamp}.log")
    config.option.log_file_level = "INFO"


# --------------------------------------------------------------------------- global-state guard
#
# A test that patches a module global must patch it in a way that is undone when the test ends.
# A direct assignment is not: it leaks into every later test in the session, and the failures then
# surface far from the cause — a stage resolving its paths against another test's temp directory,
# a literature file "missing on disk", a --check reporting that outputs differ. This fixture makes
# that failure land on the test that caused it.
#
# Add a name to _GUARDED when a new module-level global is worth protecting.

import pytest

from pipeline import common as _common

_GUARDED = ("REPO_ROOT", "CONFIG_DIR", "DATA_DIR", "OUTPUTS_DIR", "INTERMEDIATE_DIR",
            "PROTEIN_SET_OUT_DIR", "COMPOSITION_DIR", "DIGEST_DIR", "LITERATURE_INVENTORY_DIR",
            "MASS_FRACTIONS_DIR", "STANDARD_DIR", "STANDARD_PLOTS_DIR", "UNCERTAINTY_DIR",
            "STRESS_DIR", "SENSITIVITY_DIR", "MATCH_OUT_DIR", "COMPARISON_OUT_DIR", "USDA_OUT_DIR")


@pytest.fixture(autouse=True)
def _no_global_leaks():
    before = {name: getattr(_common, name) for name in _GUARDED if hasattr(_common, name)}
    yield
    changed = {name: (was, getattr(_common, name)) for name, was in before.items()
               if getattr(_common, name) != was}
    assert not changed, (
        "this test changed pipeline.common globals and did not restore them: "
        + "; ".join(f"{n}: {w} -> {g}" for n, (w, g) in changed.items())
        + ". Use monkeypatch.setattr(common, name, value) rather than assigning directly.")