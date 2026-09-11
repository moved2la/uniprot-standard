"""Flags: what the rules could not settle from a database field.

A flag is a row in outputs/flags.tsv. It is closed only by a section in
config/protein_set_decisions.ini of the form

    [flag.<subject>]
    rule     = <rule name from the flag row>
    decision = D<n>
    <resolution keys specific to the rule>

Code never closes a flag on its own.
"""

from __future__ import annotations

from pathlib import Path

from pipeline import common

FLAG_HEADER = ["raised", "stage", "subject", "rule", "detail"]


def append_flag(path: Path, *, stage: str, subject: str, rule: str, detail: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists()
    with open(path, "a", encoding="utf-8") as fh:
        if new:
            fh.write("\t".join(FLAG_HEADER) + "\n")
        fh.write("\t".join([common.iso_now(), stage, subject, rule,
                            detail.replace("\t", " ").replace("\n", " ")]) + "\n")


def reset_flags(path: Path, stage: str) -> None:
    """Remove rows previously raised by `stage` so a rerun does not duplicate them."""
    if not path.exists():
        return
    rows = common.read_tsv(path)
    keep = [r for r in rows if r["stage"] != stage]
    common.write_tsv(path, FLAG_HEADER, [[r[h] for h in FLAG_HEADER] for r in keep])


def open_flags(path: Path, decisions) -> list[dict[str, str]]:
    """Flags in the file that have no matching [flag.<subject>] closure."""
    if not path.exists():
        return []
    rows = common.read_tsv(path)
    return [r for r in rows if not decisions.has_section(f"flag.{r['subject']}")]


def closure(decisions, subject: str):
    sec = f"flag.{subject}"
    return decisions[sec] if decisions.has_section(sec) else None
