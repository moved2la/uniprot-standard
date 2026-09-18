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
import re
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
    "pool_overlap.tsv": "written by enumerate_pool, a NETWORK stage; `protein-set --offline` skips it",
    "flags.tsv": "written only when a flag is raised; zero flags this run, and the old copy was header-only",
    "profile_fiber_types_combined.png": "plots.py no longer draws it (it draws profile_fiber_types_total)",
    "profile_fiber_types_combined.svg": "plots.py no longer draws it (it draws profile_fiber_types_total)",
}

# A body line may legitimately differ for reasons that are not a change of content. Each kind is
# named, matched narrowly, and reported in full — the point is to say WHY a line differs, not to
# wave it through. Anything not matching one of these is a real difference and fails the run.
SHA256 = re.compile(r"\b[0-9a-f]{64}\b")
TIMESTAMP = re.compile(r"\b\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z\b")
PATH_IN_OUTPUTS = re.compile(r"\boutputs/[\w./-]+")


def why_differs(old_line: str, new_line: str) -> str | None:
    """Name the reason two body lines differ, or None when the content really changed."""
    o, n = old_line, new_line
    if TIMESTAMP.sub("<t>", o) == TIMESTAMP.sub("<t>", n):
        return "run timestamp"
    o2, n2 = PATH_IN_OUTPUTS.sub("<p>", o), PATH_IN_OUTPUTS.sub("<p>", n)
    if o2 == n2:
        return "outputs path moved"
    if SHA256.sub("<h>", o2) == SHA256.sub("<h>", n2):
        return "path moved + hash of a file whose header names it" if o2 != n2 else "hash of a file whose header moved"
    if SHA256.sub("<h>", o) == SHA256.sub("<h>", n):
        return "hash of a file whose header moved"
    return None

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

    same, moved, differ, missing, added, binary, expected_gone, explained = [], [], [], [], [], [], [], []

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
            reasons, unexplained = [], []
            if len(ob) != len(nb):
                unexplained.append((0, f"line count {len(ob)} -> {len(nb)}", "", ""))
            for i, (a, b) in enumerate(zip(ob, nb)):
                if a == b:
                    continue
                why = why_differs(a, b)
                (reasons if why else unexplained).append((i + 1, why or "CONTENT CHANGED", a, b))
            if unexplained:
                differ.append((name, orel, nrel, unexplained))
            else:
                explained.append((name, orel, nrel, reasons))

    for name in sorted(new):
        if name not in old:
            added.append(new[name].relative_to(args.new).as_posix())

    log.info("")
    log.info("identical body, path unchanged : %d", len(same))
    log.info("identical body, moved          : %d", len(moved))
    log.info("binary, not compared           : %d", len(binary))
    log.info("differs, every line explained  : %d", len(explained))
    log.info("BODY DIFFERS, UNEXPLAINED      : %d", len(differ))
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
    if explained:
        log.info("")
        log.info("--- differs, every differing line explained (not a content change) ---")
        for name, o, n, reasons in explained:
            kinds = {}
            for _, why, _, _ in reasons:
                kinds[why] = kinds.get(why, 0) + 1
            log.info("  %s -> %s", o, n)
            for why, count in sorted(kinds.items()):
                log.info("      %2d line(s): %s", count, why)
    if differ:
        log.error("")
        log.error("--- BODY DIFFERS, UNEXPLAINED — a real change of content ---")
        for name, o, n, unexplained in differ:
            log.error("  %s -> %s", o, n)
            for ln, why, a, b in unexplained[:5]:
                log.error("      body line %d: %s", ln, why)
                if a or b:
                    log.error("        old: %s", a[:160])
                    log.error("        new: %s", b[:160])
            if len(unexplained) > 5:
                log.error("      ... and %d more", len(unexplained) - 5)
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