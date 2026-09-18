"""Compare two outputs/ trees file by file, ignoring the header.

Step 7a moves generated files to new folders and changes nothing else. This script is the
acceptance check for that claim: every file that existed before must exist after, with an
identical body. It is tooling for one step, not a pipeline stage — it reads nothing from
config, writes nothing into outputs/, and is deleted when the step closes.

Why the body and not the whole file: a generated header names the paths it hashed and the time
it ran, so byte-identity legitimately stops at the header. "Body" here means every line that is
not a leading comment line and not the `# generated = ...` stamp.

Matching is by file NAME, not by path — that is the point of the run. Every name in the tree is
unique (checked, and asserted below), so a file that moved is still matched to its old self and
a file that changed name or vanished is reported rather than quietly passing.

Usage, one line:

    python compare_output_trees.py ../outputs_original outputs

Report goes to logs/compare_output_trees_<stamp>.log and to the screen. Exit code 0 when every file
matches, 1 when anything differs, is missing, or is unexpected.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pipeline import common  # noqa: E402

# Files the old tree holds that the current pipeline no longer produces. Each needs a reason:
# an entry here is a deliberate statement that its absence is correct, not an unnoticed loss.
NOT_REGENERATED = {
    "match_rate.tsv": "superseded in Step 6 by match_rate_per_food.tsv and match_rate_by_reference.tsv",
    "match_rate_filtered.xlsx": "not a pipeline output — the author's own working file",
}

BINARY_SUFFIXES = {".png", ".svg", ".xlsx", ".zip", ".gz", ".rar", ".pdf"}


def body(path: Path) -> list[str]:
    """Every line that is not part of the generated header."""
    if path.suffix.lower() in BINARY_SUFFIXES:
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    i = 0
    while i < len(lines) and (lines[i].startswith("#") or not lines[i].strip()):
        i += 1
    return [l for l in lines[i:] if not l.startswith("# generated")]


def index(root: Path) -> dict[str, Path]:
    files = [p for p in root.rglob("*") if p.is_file()]
    dupes = [n for n, c in Counter(p.name for p in files).items() if c > 1]
    if dupes:
        raise SystemExit(f"{root}: names are not unique, cannot match by name: {sorted(dupes)}")
    return {p.name: p for p in files}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old", type=Path, help="the outputs/ tree before the move (e.g. ../outputs_original)")
    ap.add_argument("new", type=Path, help="the outputs/ tree after the run (outputs)")
    ap.add_argument("--headers-only", action="store_true",
                    help="expect headers to differ too: report WHICH files differ, not that they do")
    args = ap.parse_args(argv)

    for p in (args.old, args.new):
        if not p.is_dir():
            raise SystemExit(f"not a directory: {p}")

    log = common.make_logger("compare_output_trees")
    old, new = index(args.old), index(args.new)
    log.info("old tree %s: %d files", args.old, len(old))
    log.info("new tree %s: %d files", args.new, len(new))

    same, moved, differ, missing, added, binary, expected_gone = [], [], [], [], [], [], []

    for name, opath in sorted(old.items()):
        if name not in new:
            (expected_gone if name in NOT_REGENERATED else missing).append(name)
            continue
        npath = new[name]
        orel = opath.relative_to(args.old).as_posix()
        nrel = npath.relative_to(args.new).as_posix()
        if npath.suffix.lower() in BINARY_SUFFIXES:
            binary.append((name, orel, nrel))
        elif body(opath) == body(npath):
            (moved if orel != nrel else same).append((name, orel, nrel))
        else:
            ob, nb = body(opath), body(npath)
            first = next((i for i, (a, b) in enumerate(zip(ob, nb)) if a != b), min(len(ob), len(nb)))
            differ.append((name, orel, nrel, len(ob), len(nb), first))

    for name in sorted(new):
        if name not in old:
            added.append(new[name].relative_to(args.new).as_posix())

    log.info("")
    log.info("identical body, path unchanged : %d", len(same))
    log.info("identical body, moved          : %d", len(moved))
    log.info("binary, not compared           : %d", len(binary))
    log.info("BODY DIFFERS                   : %d", len(differ))
    log.info("MISSING from the new tree      : %d", len(missing))
    log.info("new files not in the old tree  : %d", len(added))
    log.info("expected not-regenerated       : %d", len(expected_gone))

    if moved:
        log.info("")
        log.info("--- moved (body identical) ---")
        for name, o, n in moved:
            log.info("  %s  ->  %s", o, n)
    if binary:
        log.info("")
        log.info("--- binary, compared by presence only ---")
        for name, o, n in binary:
            log.info("  %s  ->  %s", o, n)
    if expected_gone:
        log.info("")
        log.info("--- not regenerated, by design ---")
        for name in expected_gone:
            log.info("  %s  (%s)", name, NOT_REGENERATED[name])
    if differ:
        log.error("")
        log.error("--- BODY DIFFERS ---")
        for name, o, n, lo, ln, first in differ:
            log.error("  %s -> %s : old %d body lines, new %d, first difference at body line %d",
                      o, n, lo, ln, first + 1)
    if missing:
        log.error("")
        log.error("--- MISSING (in the old tree, not produced by the run) ---")
        for name in missing:
            log.error("  %s", name)
    if added:
        log.info("")
        log.info("--- new files (not in the old tree) ---")
        for rel in added:
            log.info("  %s", rel)

    ok = not differ and not missing
    log.info("")
    log.info("RESULT: %s", "every file matches" if ok else "DIFFERENCES FOUND — see above")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
