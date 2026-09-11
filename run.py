#!/usr/bin/env python3
"""Single entry point for the repository.

    python run.py protein-set              full run: resolve -> enumerate -> verify -> build -> deltas -> tests
    python run.py protein-set --offline    skip the three network stages; rebuild from data/ on disk
    python run.py protein-set --stop-after <stage>
    python run.py test                     tests only

Every stage writes a timestamped log to outputs/logs/. A non-zero exit from any
stage stops the run; the log says why.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

NETWORK_STAGES = ["resolve_ontology_term", "enumerate_pool", "verify_accessions"]
OFFLINE_STAGES = ["build_protein_set", "isoform_processing_deltas"]


def run_stage(name: str) -> int:
    print(f"\n=== {name} ===", flush=True)
    mod = __import__(f"pipeline.{name}", fromlist=["main"])
    try:
        rc = mod.main([])
    except SystemExit as e:
        rc = int(e.code) if isinstance(e.code, int) else 1
        if e.code and not isinstance(e.code, int):
            print(e.code, file=sys.stderr)
    return rc or 0


def run_tests() -> int:
    print("\n=== tests ===", flush=True)
    return subprocess.call([sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")], cwd=ROOT)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    ps = sub.add_parser("protein-set")
    ps.add_argument("--offline", action="store_true")
    ps.add_argument("--stop-after", choices=NETWORK_STAGES + OFFLINE_STAGES)
    sub.add_parser("test")
    args = ap.parse_args()

    if args.cmd == "test":
        return run_tests()

    stages = (OFFLINE_STAGES if args.offline else NETWORK_STAGES + OFFLINE_STAGES)
    for stage in stages:
        rc = run_stage(stage)
        if rc:
            print(f"\nSTOPPED: {stage} exited {rc}. See outputs/logs/ and outputs/flags.tsv.", file=sys.stderr)
            return rc
        if args.stop_after == stage:
            return 0
    return run_tests()


if __name__ == "__main__":
    sys.exit(main())
