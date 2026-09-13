#!/usr/bin/env python3
"""Single entry point for the repository.

    python run.py protein-set              fetch literature -> resolve -> enumerate -> fetch sequences -> build -> deltas -> tests
    python run.py protein-set --offline    skip the four network stages; rebuild from data/ on disk
    python run.py composition              fetch masses -> fetch ptmlist -> composition -> ptm disclosure -> tests
    python run.py composition --offline    skip the two network stages; recompute from data/ on disk
    python run.py mass-fractions           fetch literature -> digest -> literature inventory -> mass fractions -> tests
    python run.py mass-fractions --offline skip the fetch; recompute from data/ on disk
    python run.py <command> --stop-after <stage>
    python run.py test                     tests only

Every stage writes a timestamped log to logs/, and everything run.py itself prints —
banners, stop messages, crash tracebacks, the full test report — goes to
logs/run_<command>_<stamp>.log as well, so the screen and the logs never differ.
A non-zero exit from any stage stops the run; the log says why. Tests at the end of a
command skip the currency checks of later commands (their outputs are rebuilt next).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

COMMANDS = {
    # command: (network stages, offline stages). Each stage is pipeline/<name>.py with a main(argv).
    "protein-set": (["fetch_literature", "resolve_ontology_term", "enumerate_pool", "fetch_sequences"],
                    ["build_protein_set", "isoform_processing_deltas"]),
    "composition": (["fetch_amino_acid_masses", "fetch_ptmlist"],
                    ["composition", "ptm_disclosure"]),
    "mass-fractions": (["fetch_literature"],
                       ["digest", "literature_inventory", "mass_fractions"]),
}


RUN_LOG: Path | None = None
COMMAND_ORDER = ["protein-set", "composition", "mass-fractions"]


def say(text: str, err: bool = False) -> None:
    """Print AND append to the run log: nothing shown on screen is absent from logs/."""
    print(text, file=sys.stderr if err else sys.stdout, flush=True)
    if RUN_LOG is not None:
        with open(RUN_LOG, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")


def run_stage(name: str) -> int:
    say(f"\n=== {name} ===")
    mod = __import__(f"pipeline.{name}", fromlist=["main"])
    try:
        rc = mod.main([])
    except SystemExit as e:
        rc = int(e.code) if isinstance(e.code, int) else 1
        if e.code and not isinstance(e.code, int):
            say(str(e.code), err=True)          # a [STOP] message: on screen and in the run log
    except Exception:  # catch-all: nothing leaves a stage without a record on disk
        import traceback
        crash = ROOT / "logs" / f"{name}_crash_{STAMP}.log"
        crash.write_text(f"stage {name} crashed {STAMP}\n\n" + traceback.format_exc(), encoding="utf-8")
        say(traceback.format_exc(), err=True)
        say(f"CRASH: {name}; traceback written to {crash.relative_to(ROOT).as_posix()}", err=True)
        rc = 1
    return rc or 0


def run_tests(command: str | None = None) -> int:
    """Run the tests; the full console report is written to logs/pytest_report_<stamp>.log and to the run
    log. When called at the end of a command, the currency tests of LATER commands are skipped
    (their inputs were just regenerated and they are rebuilt by the next command); `python run.py
    test` runs everything."""
    say("\n=== tests ===")
    import os
    env = dict(os.environ)
    if command:
        env["UNIPROT_STANDARD_COMMAND"] = command
    proc = subprocess.run([sys.executable, "-m", "pytest", "-q", "-rs", str(ROOT / "tests")],
                          cwd=ROOT, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
    report = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
    (ROOT / "logs" / f"pytest_report_{STAMP}.log").write_text(report, encoding="utf-8")
    say(report)
    return proc.returncode


def main() -> int:
    global RUN_LOG, STAMP
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, (network, offline) in COMMANDS.items():
        cp = sub.add_parser(name)
        cp.add_argument("--offline", action="store_true")
        cp.add_argument("--stop-after", choices=network + offline)
    sub.add_parser("test")
    args = ap.parse_args()

    import datetime as _dt
    STAMP = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (ROOT / "logs").mkdir(parents=True, exist_ok=True)
    RUN_LOG = ROOT / "logs" / f"run_{args.cmd}_{STAMP}.log"
    say(f"run {args.cmd} {STAMP}  (run log: {RUN_LOG.relative_to(ROOT).as_posix()})")

    if args.cmd == "test":
        return run_tests()

    network, offline = COMMANDS[args.cmd]
    stages = (offline if args.offline else network + offline)
    for stage in stages:
        rc = run_stage(stage)
        if rc:
            say(f"\nSTOPPED: {stage} exited {rc}. See logs/ and outputs/flags.tsv.", err=True)
            return rc
        if args.stop_after == stage:
            return 0
    rc = run_tests(args.cmd)
    say(f"\n{'DONE' if rc == 0 else 'TESTS FAILED'}: {args.cmd}; run log logs/{RUN_LOG.name}", err=(rc != 0))
    return rc


if __name__ == "__main__":
    sys.exit(main())
