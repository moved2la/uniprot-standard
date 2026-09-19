#!/usr/bin/env python3
"""Single entry point for the repository.

    python run.py fetch-literature         hash the literature files on disk -> manifest, map, README (never downloads, D80)
    python run.py protein-set              resolve -> enumerate -> fetch sequences -> build -> deltas -> tests
    python run.py protein-set --offline    skip the three network stages; rebuild from data/ on disk
    python run.py composition              fetch masses -> fetch ptmlist -> composition -> ptm disclosure -> tests
    python run.py composition --offline    skip the two network stages; recompute from data/ on disk
    python run.py mass-fractions           digest -> literature inventory -> mass fractions -> tests (offline)
    python run.py standard                 aggregate -> non_protein_metabolite_pools -> uncertainty -> stress -> plots -> tests (all offline)
    python run.py <command> --stop-after <stage>
    python run.py usda                     read the USDA FoodData Central archives in data/usda/ -> tests (offline)
    python run.py comparison               the calculated standard beside Gorissen 2018's measurement (offline)
    python run.py match                    Match Rate: every food against every reference (offline)
    python run.py blood-protein-set        the blood category (Step 7b): enumerate the two pools by published
                                           accession -> fetch the sequences the store lacks (additive) ->
                                           build config/blood/ -> deltas -> tests; --offline skips the two network stages
    python run.py blood-composition        composition and PTM disclosure for the blood set (offline; the masses,
                                           symbols and ptmlist fetched by `composition` are shared)
    python run.py test                     tests only
    python run.py excerpt                  tooling: bounded excerpts of the large generated files
                                           into excerpts/<stamp>/ and a tarball (no tests; not the record)

Every stage writes its own timestamped log to logs/. In addition, logs/run_<command>_<stamp>.log
is the MASTER LOG: sys.stdout and sys.stderr are tee'd into it for the whole run, so it holds
everything that appeared on screen — stage log lines, banners, stop messages, crash tracebacks,
the full test report — in the exact order it appeared. Look there first.
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
    # fetch-literature is its own command (D80, D81): it hashes what is on disk, never downloads,
    # and is run when a literature file is added or replaced — not on every protein-set run.
    "fetch-literature": ([],
                         ["fetch_literature"]),
    "protein-set": (["resolve_ontology_term", "enumerate_pool", "fetch_sequences"],
                    ["build_protein_set", "isoform_processing_deltas"]),
    "composition": (["fetch_amino_acid_masses", "fetch_ptmlist"],
                    ["composition", "ptm_disclosure"]),
    "mass-fractions": ([],
                       ["digest", "literature_inventory", "mass_fractions"]),
    "standard": ([],
                 ["aggregate", "non_protein_metabolite_pools", "uncertainty", "stress", "plots"]),
    "usda": ([],
             ["usda"]),
    "comparison": ([],
                   ["comparison"]),
    "match": ([],
              ["match"]),
    # Step 7b: the blood category, built with the muscle stages passed --category blood. Each blood
    # command mirrors the muscle command it adapts (the adaptation list is kept per command); Step 7c
    # folds them into `<command> <category>`. Placement in the run order is provisional (after match).
    "blood-protein-set": (["enumerate_pool", "fetch_sequences"],
                          ["build_protein_set", "isoform_processing_deltas"]),
    # the masses, symbols and ptmlist are shared data already fetched by `composition`; blood has no network stage here
    "blood-composition": ([],
                          ["composition", "ptm_disclosure"]),
}

# The category a command builds; stages of these commands are called with --category <name>.
COMMAND_CATEGORY = {"blood-protein-set": "blood", "blood-composition": "blood"}


RUN_LOG: Path | None = None
COMMAND_ORDER = ["fetch-literature", "protein-set", "composition", "mass-fractions", "standard",
                 "usda", "comparison", "match", "blood-protein-set", "blood-composition"]


class _Tee:
    """A stream that writes to the screen and appends the same bytes, in the same order, to the
    master log. Installed over sys.stdout and sys.stderr for the whole run, so every stage's
    log lines, every print, every traceback, and the test report land in one file exactly as
    they appeared on screen."""

    def __init__(self, screen, path: Path):
        self.screen = screen
        self.path = path

    def write(self, s: str) -> int:
        self.screen.write(s)
        self.screen.flush()
        with open(self.path, "a", encoding="utf-8", errors="replace") as fh:
            fh.write(s)
        return len(s)

    def flush(self) -> None:
        self.screen.flush()

    def isatty(self) -> bool:
        return False

    def fileno(self):
        return self.screen.fileno()


def say(text: str, err: bool = False) -> None:
    """Print to the screen; the tee carries it into the master log."""
    print(text, file=sys.stderr if err else sys.stdout, flush=True)


def run_stage(name: str, argv: list[str] | None = None) -> int:
    argv = argv or []
    say(f"\n=== {name}{' ' + ' '.join(argv) if argv else ''} ===")
    mod = __import__(f"pipeline.{name}", fromlist=["main"])
    try:
        rc = mod.main(list(argv))
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
    sub.add_parser("excerpt")
    args = ap.parse_args()

    import datetime as _dt
    STAMP = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (ROOT / "logs").mkdir(parents=True, exist_ok=True)
    RUN_LOG = ROOT / "logs" / f"run_{args.cmd}_{STAMP}.log"
    sys.stdout = _Tee(sys.stdout, RUN_LOG)     # from here on, screen and master log are the same stream
    sys.stderr = _Tee(sys.stderr, RUN_LOG)
    say(f"run {args.cmd} {STAMP}  (master log: {RUN_LOG.relative_to(ROOT).as_posix()})")

    if args.cmd == "test":
        return run_tests()
    if args.cmd == "excerpt":
        return run_stage("excerpt")

    network, offline = COMMANDS[args.cmd]
    stages = (offline if args.offline else network + offline)
    category = COMMAND_CATEGORY.get(args.cmd)
    stage_argv = ["--category", category] if category else []
    for stage in stages:
        rc = run_stage(stage, stage_argv)
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
