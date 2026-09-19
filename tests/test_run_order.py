"""run.py's commands and common.COMMAND_ORDER (docs/run_order.md) must name the same commands: a
command missing from the order would never be recognised as "earlier" by the currency tests."""
import importlib.util
from pathlib import Path

from pipeline import common

ROOT = Path(__file__).resolve().parent.parent


def _run_module():
    spec = importlib.util.spec_from_file_location("run_py", ROOT / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_command_is_in_the_run_order_once():
    run = _run_module()
    assert sorted(run.COMMANDS) == sorted(common.COMMAND_ORDER)
    assert len(set(common.COMMAND_ORDER)) == len(common.COMMAND_ORDER)


def test_dependency_order():
    o = common.COMMAND_ORDER
    assert o.index("standard") < o.index("blood-protein-set")          # the sequence store is additive on top of muscle's (S1)
    assert o.index("blood-standard") < o.index("composite")            # the composite mixes the two standards
    assert o.index("standard") < o.index("composite")
    assert o.index("composite") < o.index("match")                     # match reads the composite as a reference
    assert o.index("usda") < o.index("match")
