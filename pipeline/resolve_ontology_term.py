"""Resolve each tier's lookup name to exactly one Gene Ontology term, by code.

Inputs
  config/protein_set_decisions.ini
      [ontology_source] url            -- where the ontology file is downloaded from
      [tier.N]          lookup_name    -- the exact term name to match (author decision)

What it does
  1. Downloads the ontology file unchanged to data/gene-ontology/ and records the
     resolved URL, the file's own version header lines, retrieval time, and SHA-256.
  2. Parses [Term] stanzas: id, name, namespace, definition, is_a, part_of, obsolete.
  3. For each tier, matches lookup_name against term names by exact string equality.
     Zero matches, or more than one match, is a hard stop -- the run exits non-zero
     and the condition is written to outputs/flags.tsv. No related-term search is
     performed.
  4. Writes the matched term's full descendant subtree (via is_a and part_of),
     its ancestors to the root, and its immediate neighbours (direct parents = next
     broader, direct children = next narrower).

Outputs
  data/gene-ontology/<file as downloaded>
  data/gene-ontology/source.ini
  data/gene-ontology/resolved_terms.ini          -- consumed by enumerate_pool.py
  data/gene-ontology/tier<N>_subtree.tsv
  data/gene-ontology/tier<N>_ancestors.tsv
  data/gene-ontology/tier<N>_neighbours.tsv
"""

from __future__ import annotations

import argparse
import sys
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

from pipeline import common
from pipeline.flags import append_flag

RELATIONS_FOLLOWED = ("is_a", "part_of")


@dataclass
class Term:
    id: str
    name: str = ""
    namespace: str = ""
    definition: str = ""
    obsolete: bool = False
    # (relation, parent_id) pairs; only RELATIONS_FOLLOWED are stored
    parents: list[tuple[str, str]] = field(default_factory=list)


@dataclass
class Ontology:
    terms: dict[str, Term]
    header: dict[str, str]  # format-version, data-version, date, ...
    children: dict[str, list[tuple[str, str]]]  # parent_id -> [(relation, child_id)]


def parse_obo(text: str) -> Ontology:
    """Parse the OBO 1.2/1.4 flat format. Only [Term] stanzas are kept."""
    header: dict[str, str] = {}
    terms: dict[str, Term] = {}
    current: Term | None = None
    in_term = False
    seen_stanza = False
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("["):
            seen_stanza = True
            if current is not None and in_term:
                terms[current.id] = current
            in_term = line == "[Term]"
            current = Term(id="") if in_term else None
            continue
        if current is None:
            if not seen_stanza and ":" in line:
                key, _, val = line.partition(":")
                header.setdefault(key.strip(), val.strip())
            continue
        key, _, val = line.partition(":")
        key = key.strip()
        val = val.strip()
        if key == "id":
            current.id = val
        elif key == "name":
            current.name = val
        elif key == "namespace":
            current.namespace = val
        elif key == "def":
            current.definition = val
        elif key == "is_obsolete":
            current.obsolete = val.lower().startswith("true")
        elif key == "is_a":
            current.parents.append(("is_a", val.split("!")[0].strip()))
        elif key == "relationship":
            parts = val.split("!")[0].split()
            if len(parts) >= 2 and parts[0] in RELATIONS_FOLLOWED:
                current.parents.append((parts[0], parts[1]))
    if current is not None and in_term and current.id:
        terms[current.id] = current
    children: dict[str, list[tuple[str, str]]] = {}
    for t in terms.values():
        for rel, pid in t.parents:
            children.setdefault(pid, []).append((rel, t.id))
    return Ontology(terms=terms, header=header, children=children)


def exact_name_matches(onto: Ontology, name: str) -> list[Term]:
    """Exact, case-sensitive equality on the term name. Obsolete terms are reported
    separately by the caller; here they are returned too so nothing is hidden."""
    return [t for t in onto.terms.values() if t.name == name]


def descendants(onto: Ontology, root_id: str) -> list[tuple[str, int, str, str]]:
    """Breadth-first over children via RELATIONS_FOLLOWED.
    Returns (term_id, depth, relation_to_parent, parent_id), root included at depth 0."""
    out = [(root_id, 0, "", "")]
    seen = {root_id}
    q = deque([(root_id, 0)])
    while q:
        tid, depth = q.popleft()
        for rel, cid in sorted(onto.children.get(tid, [])):
            if cid in seen:
                continue
            seen.add(cid)
            out.append((cid, depth + 1, rel, tid))
            q.append((cid, depth + 1))
    return out


def ancestors(onto: Ontology, start_id: str) -> list[tuple[str, int, str, str]]:
    """Breadth-first over parents. Returns (term_id, distance, relation, child_id)."""
    out = []
    seen = {start_id}
    q = deque([(start_id, 0)])
    while q:
        tid, dist = q.popleft()
        for rel, pid in sorted(onto.terms[tid].parents):
            if pid in seen:
                continue
            seen.add(pid)
            out.append((pid, dist + 1, rel, tid))
            q.append((pid, dist + 1))
    return out


def neighbours(onto: Ontology, tid: str) -> list[tuple[str, str, str]]:
    """Immediate parents (next broader) and children (next narrower).
    Returns (direction, relation, term_id)."""
    out = [("broader", rel, pid) for rel, pid in sorted(onto.terms[tid].parents)]
    out += [("narrower", rel, cid) for rel, cid in sorted(onto.children.get(tid, []))]
    return out


def _term_row(onto: Ontology, tid: str) -> tuple[str, str, str, str]:
    t = onto.terms.get(tid)
    if t is None:  # parent outside the parsed file (should not happen for go-basic)
        return (tid, "", "", "")
    return (t.name, t.namespace, t.definition, "true" if t.obsolete else "false")


def download(url: str, dest_dir: Path, log) -> tuple[Path, str]:
    """Download the ontology file unchanged. Returns (path, final_url)."""
    import requests

    dest_dir.mkdir(parents=True, exist_ok=True)
    log.info("downloading %s", url)
    resp = requests.get(url, timeout=120, allow_redirects=True)
    resp.raise_for_status()
    final_url = resp.url
    filename = final_url.rsplit("/", 1)[-1] or "ontology.obo"
    path = dest_dir / filename
    path.write_bytes(resp.content)
    log.info("saved %s (%d bytes) from %s", path, len(resp.content), final_url)
    return path, final_url


def resolve(onto: Ontology, decisions, log, onto_dir: Path, flags_path: Path) -> dict[str, str]:
    """Match every tier's lookup_name and write the per-tier tables.
    Returns {tier: term_id}. Raises SystemExit on any unresolvable tier."""
    resolved: dict[str, str] = {}
    hard_stop = False
    for tier in common.tiers_from_decisions(decisions):
        sec = decisions[f"tier.{tier}"]
        name = sec["lookup_name"]
        matches = exact_name_matches(onto, name)
        live = [m for m in matches if not m.obsolete]
        log.info("tier %s: lookup_name=%r exact matches=%d (non-obsolete=%d)",
                 tier, name, len(matches), len(live))
        if len(live) != 1:
            hard_stop = True
            detail = "; ".join(f"{m.id} [{m.namespace}]{' OBSOLETE' if m.obsolete else ''}"
                               for m in matches) or "no term has this exact name"
            append_flag(flags_path, stage="resolve_ontology_term", subject=f"tier.{tier}",
                        rule="exact_name_match", detail=f"lookup_name={name!r}: {detail}")
            log.error("tier %s: cannot resolve -- flagged", tier)
            continue
        term = live[0]
        resolved[tier] = term.id
        log.info("tier %s: resolved to %s (%s) namespace=%s", tier, term.id, term.name, term.namespace)

        sub = descendants(onto, term.id)
        common.write_tsv(
            onto_dir / f"tier{tier}_subtree.tsv",
            ["term_id", "depth", "relation_to_parent", "parent_id", "name", "namespace", "definition", "obsolete"],
            [[tid, d, rel, pid, *_term_row(onto, tid)] for tid, d, rel, pid in sub],
        )
        anc = ancestors(onto, term.id)
        common.write_tsv(
            onto_dir / f"tier{tier}_ancestors.tsv",
            ["term_id", "distance", "relation", "child_id", "name", "namespace", "definition", "obsolete"],
            [[tid, d, rel, cid, *_term_row(onto, tid)] for tid, d, rel, cid in anc],
        )
        nb = neighbours(onto, term.id)
        common.write_tsv(
            onto_dir / f"tier{tier}_neighbours.tsv",
            ["direction", "relation", "term_id", "name", "namespace", "definition", "obsolete"],
            [[dirn, rel, tid, *_term_row(onto, tid)] for dirn, rel, tid in nb],
        )
        log.info("tier %s: subtree=%d terms (incl. root), ancestors=%d, neighbours=%d",
                 tier, len(sub), len(anc), len(nb))
    if hard_stop:
        raise SystemExit("one or more tiers could not be resolved by exact name; see outputs/flags.tsv")
    return resolved


def write_resolved(resolved: dict[str, str], onto: Ontology, source: dict[str, str],
                   decisions, path: Path) -> None:
    cp = common.new_ini()
    cp["ontology"] = {
        "file": source["file"],
        "url_requested": source["url_requested"],
        "url_final": source["url_final"],
        "sha256": source["sha256"],
        "retrieved": source["retrieved"],
        "format_version": onto.header.get("format-version", ""),
        "data_version": onto.header.get("data-version", ""),
        "date": onto.header.get("date", ""),
        "relations_followed": ",".join(RELATIONS_FOLLOWED),
        "term_count": str(len(onto.terms)),
    }
    for tier, tid in resolved.items():
        t = onto.terms[tid]
        sec = decisions[f"tier.{tier}"]
        cp[f"tier.{tier}"] = {
            "seed_word": sec.get("seed_word", ""),
            "seed_decision": sec.get("seed_decision", ""),
            "lookup_name": sec["lookup_name"],
            "lookup_decision": sec.get("lookup_decision", ""),
            "term_id": tid,
            "term_name": t.name,
            "namespace": t.namespace,
            "definition": t.definition,
            "subtree_file": f"tier{tier}_subtree.tsv",
            "subtree_size": str(len(descendants(onto, tid))),
        }
    common.write_ini(cp, path, [
        "resolved_terms.ini -- GENERATED by resolve_ontology_term.py. DO NOT EDIT BY HAND.",
        "One Gene Ontology term per tier, matched by exact name equality from",
        "config/protein_set_decisions.ini. Consumed by enumerate_pool.py.",
    ])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ontology-file", type=Path, default=None,
                    help="use an already-downloaded ontology file instead of downloading")
    args = ap.parse_args(argv)

    log = common.make_logger("resolve_ontology_term")
    decisions = common.read_ini(common.DECISIONS_INI)
    url = decisions["ontology_source"]["url"]
    onto_dir = common.ONTOLOGY_DIR
    onto_dir.mkdir(parents=True, exist_ok=True)

    if args.ontology_file:
        path, final_url = args.ontology_file, f"(local file) {args.ontology_file}"
        log.info("using local ontology file %s", path)
    else:
        path, final_url = download(url, onto_dir, log)
    retrieved = common.iso_now()
    sha = common.sha256_file(path)
    log.info("sha256 %s", sha)

    onto = parse_obo(path.read_text(encoding="utf-8"))
    log.info("parsed %d terms; header: %s", len(onto.terms),
             {k: onto.header.get(k, "") for k in ("format-version", "data-version", "date")})

    source = {"file": path.name, "url_requested": url, "url_final": final_url,
              "sha256": sha, "retrieved": retrieved}
    src = common.new_ini()
    src["source"] = {**source, "format_version": onto.header.get("format-version", ""),
                     "data_version": onto.header.get("data-version", ""),
                     "date": onto.header.get("date", "")}
    common.write_ini(src, onto_dir / "source.ini",
                     ["source.ini -- GENERATED by resolve_ontology_term.py. Provenance of the ontology file."])

    resolved = resolve(onto, decisions, log, onto_dir, common.FLAGS_TSV)
    write_resolved(resolved, onto, source, decisions, common.RESOLVED_TERMS_INI)
    log.info("wrote %s", common.RESOLVED_TERMS_INI)
    return 0


if __name__ == "__main__":
    sys.exit(main())
