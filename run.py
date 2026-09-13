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

Every stage writes a timestamped log to logs/. A non-zero exit from any
stage stops the run; the log says why.
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


def run_stage(name: str) -> int:
    print(f"\n=== {name} ===", flush=True)
    mod = __import__(f"pipeline.{name}", fromlist=["main"])
    try:
        rc = mod.main([])
    except SystemExit as e:
        rc = int(e.code) if isinstance(e.code, int) else 1
        if e.code and not isinstance(e.code, int):
            print(e.code, file=sys.stderr)
    except Exception:  # catch-all: nothing leaves a stage without a record on disk
        import datetime as _dt
        import traceback
        logs = ROOT / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        crash = logs / f"{name}_crash_{stamp}.log"
        crash.write_text(f"stage {name} crashed {stamp}\n\n" + traceback.format_exc(), encoding="utf-8")
        traceback.print_exc()
        print(f"\nCRASH: {name}; traceback written to {crash.relative_to(ROOT).as_posix()}", file=sys.stderr)
        rc = 1
    return rc or 0


def run_tests() -> int:
    print("\n=== tests ===", flush=True)
    return subprocess.call([sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")], cwd=ROOT)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, (network, offline) in COMMANDS.items():
        cp = sub.add_parser(name)
        cp.add_argument("--offline", action="store_true")
        cp.add_argument("--stop-after", choices=network + offline)
    sub.add_parser("test")
    args = ap.parse_args()

    if args.cmd == "test":
        return run_tests()

    network, offline = COMMANDS[args.cmd]
    stages = (offline if args.offline else network + offline)
    for stage in stages:
        rc = run_stage(stage)
        if rc:
            print(f"\nSTOPPED: {stage} exited {rc}. See logs/ and outputs/flags.tsv.", file=sys.stderr)
            return rc
        if args.stop_after == stage:
            return 0
    return run_tests()


if __name__ == "__main__":
    sys.exit(main())
