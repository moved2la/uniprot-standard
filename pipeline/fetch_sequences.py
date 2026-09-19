"""Fetch every pool member from UniProt and verify it, by code.

Inputs
  data/tier<N>_pool.tsv      (from enumerate_pool.py; the union over tiers is fetched)

Per accession
  - entry JSON  : rest.uniprot.org/uniprotkb/<acc>.json  (saved unchanged to data/uniprot_raw/)
  - FASTA       : rest.uniprot.org/uniprotkb/<acc>.fasta
  - verification: MD5 of the FASTA sequence == the MD5 UniProt publishes in the JSON;
                  FASTA sequence == JSON sequence; length agrees; only standard letters.
                  Any failure is a hard failure (exit non-zero).
  - captured    : processing-relevant features (see CAPTURED_FEATURE_TYPES), the
                  alternative-products isoform table with each isoform's sequence
                  fetched and MD5'd, the Gene Ontology cross-references with their
                  evidence codes, tissue-specificity text (recorded, never parsed).

Output
  data/uniprot_sequences.ini              -- everything build_protein_set.py reads; ONE store, shared
                                             by every category, written additively (D110)
  data/uniprot_raw/<category>/<acc>.json  -- raw entries, in the folder of the category whose fetch
                                             added them; nothing is ever deleted here (D61)

Additive store (Step 7b, D110)
  The store is rewritten as a MERGE: the sections fetched in this run are rendered; every section
  already in the file and not fetched now is carried over VERBATIM, in one sorted order. A
  category run (`--category <name>`) fetches only the accessions of its pools that the store does
  not already hold, so muscle's entries are never touched by blood's fetch and a muscle refetch
  cannot drop blood's. The header keeps the `fetched` / `uniprot_release` lines of the run that
  built the store from a full pool and adds `added.<category>.*` lines for each additive run; a
  section added by a category carries its own `fetched` and `uniprot_release`.
"""

from __future__ import annotations

import argparse
import json
import sys

from pipeline import common
from pipeline.uniprot_client import UniProtClient

STAGE = "fetch_sequences"
CAPTURED_FEATURE_TYPES = ("Chain", "Initiator methionine", "Signal", "Propeptide",
                          "Transit peptide", "Peptide", "Alternative sequence")
STANDARD_LETTERS = set(common.AMINO_ACIDS)


def pool_accessions(category: str | None = None) -> dict[str, list[str]]:
    """accession -> [pools it belongs to]: every data/tier<N>_pool.tsv for muscle, or the
    data/<category>/pool_<name>.tsv of the pools the category's decisions file declares."""
    out: dict[str, list[str]] = {}
    if category is None:
        for path in sorted(common.DATA_DIR.glob("tier*_pool.tsv")):
            tier = path.name[len("tier"):].split("_", 1)[0]
            for row in common.read_tsv(path):
                out.setdefault(row["accession"], []).append(tier)
        if not out:
            raise FileNotFoundError("no data/tier<N>_pool.tsv files found; run enumerate_pool first")
        return out
    for pool in common.pool_names(category):
        for row in common.read_tsv(common.pool_file(category, pool)):
            out.setdefault(row["accession"], []).append(pool)
    if not out:
        raise FileNotFoundError(f"no pool tables under data/{category}/; run enumerate_pool --category {category} first")
    return out


# ----------------------------------------------------------------------------- the store as text
#
# The store is read back as TEXT, not through configparser: a value re-read and re-rendered could
# lose an inline ';' or '#', and a section carried over must be the bytes it was. Blocks are
# keyed by section name; the header is the leading '# ' lines.

def read_store_text(path) -> tuple[list[str], dict[str, str]]:
    """(header lines without the '# ', {section name: the section's text block, no trailing blank line})."""
    if not path.exists():
        return [], {}
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")
    header: list[str] = []
    i = 0
    while i < len(lines) and lines[i].startswith("#"):
        header.append(lines[i][1:].lstrip(" ") if lines[i] != "#" else "")
        i += 1
    blocks: dict[str, str] = {}
    name, body = None, []
    for line in lines[i:]:
        if line.startswith("[") and line.rstrip().endswith("]"):
            if name is not None:
                blocks[name] = "\n".join(body).rstrip("\n")
            name, body = line.strip()[1:-1], [line.rstrip("\n")]
        elif name is not None:
            body.append(line)
    if name is not None:
        blocks[name] = "\n".join(body).rstrip("\n")
    return header, blocks


def render_section(name: str, items: dict[str, str]) -> str:
    """One section exactly as common.render_ini renders it (keys aligned within the section)."""
    width = max((len(k) for k in items), default=0)
    return "\n".join([f"[{name}]"] + [f"{k.ljust(width)} = {v}" for k, v in items.items()])


def merge_header(existing: list[str], category: str | None, release: str, n_new: int, n_total: int, n_failures: int) -> list[str]:
    """The header of the merged store. A full (muscle) run writes its own lines and keeps any
    `added.<category>.*` lines; an additive run keeps the existing lines and (re)writes its own
    `added.<category>.*` lines. `accessions` always counts the whole store."""
    keep = [l for l in existing if l and "=" in l and not l.split("=", 1)[0].strip().startswith(("added.", "accessions", "failures"))]
    added = [l for l in existing if l and "=" in l and l.split("=", 1)[0].strip().startswith("added.")
             and not (category and l.split("=", 1)[0].strip().startswith(f"added.{category}."))]
    if category is None:
        keep = [
            f"fetched         = {common.iso_now()}",
            f"uniprot_release = {release}",
            "source          = https://rest.uniprot.org/uniprotkb/",
            f"captured_feature_types = {'; '.join(CAPTURED_FEATURE_TYPES)}",
        ]
    else:
        added.append(f"added.{category}.fetched = {common.iso_now()}")
        added.append(f"added.{category}.uniprot_release = {release}")
        added.append(f"added.{category}.accessions = {n_new}")
    first = "uniprot_sequences.ini -- GENERATED by fetch_sequences.py. DO NOT EDIT BY HAND."
    keep = [l for l in keep if l != first]
    return [first] + keep + added + [f"accessions      = {n_total}", f"failures        = {n_failures}"]


def merged_store_text(header: list[str], fetched: dict[str, dict[str, str]], carried: dict[str, str]) -> str:
    """Header, then every section in one sorted order: fetched sections rendered, carried ones verbatim."""
    lines = [f"# {l}" if l else "#" for l in header]
    lines.append("")
    names = sorted(set(fetched) | set(carried))
    for name in names:
        lines.append(render_section(name, fetched[name]) if name in fetched else carried[name])
        lines.append("")
    return "\n".join(lines)


def _loc(feature: dict) -> tuple[str, str]:
    loc = feature.get("location", {})
    s, e = loc.get("start", {}), loc.get("end", {})

    def fmt(p):
        v = p.get("value")
        mod = p.get("modifier", "EXACT")
        v = "" if v is None else str(v)
        return v if mod == "EXACT" else f"{v}({mod})"
    return fmt(s), fmt(e)


def extract(entry: dict) -> dict:
    """Pull the fields we record from an entry JSON. Pure function, testable offline."""
    seq = entry.get("sequence", {})
    audit = entry.get("entryAudit", {})
    genes = entry.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
    pname = (entry.get("proteinDescription", {}).get("recommendedName", {})
             .get("fullName", {}).get("value", ""))
    out = {
        "accession": entry.get("primaryAccession", ""),
        "entry_name": entry.get("uniProtkbId", ""),
        "entry_type": entry.get("entryType", ""),
        "gene": gene,
        "protein_name": pname,
        "organism_id": str(entry.get("organism", {}).get("taxonId", "")),
        "sequence": seq.get("value", ""),
        "length": str(seq.get("length", "")),
        "md5_uniprot": str(seq.get("md5", "")).upper(),
        "seq_version": str(audit.get("sequenceVersion", "")),
        "entry_version": str(audit.get("entryVersion", "")),
        "last_seq_update": audit.get("lastSequenceUpdateDate", ""),
        "features": [],
        "n_features_total": str(len(entry.get("features", []))),
        "isoforms": [],
        "annotations": [],
        "tissue_specificity": "",
    }
    for f in entry.get("features", []):
        if f.get("type") in CAPTURED_FEATURE_TYPES:
            s, e = _loc(f)
            out["features"].append({
                "type": f.get("type", ""), "start": s, "end": e,
                "description": f.get("description", ""), "feature_id": f.get("featureId", ""),
            })
    for c in entry.get("comments", []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms", []):
                out["isoforms"].append({
                    "ids": ";".join(iso.get("isoformIds", [])),
                    "name": iso.get("name", {}).get("value", ""),
                    "status": iso.get("isoformSequenceStatus", ""),
                    "sequence_ids": ";".join(iso.get("sequenceIds", [])),
                    "note": " ".join(t.get("value", "") for t in iso.get("note", {}).get("texts", [])),
                })
        elif c.get("commentType") == "TISSUE SPECIFICITY":
            out["tissue_specificity"] = " ".join(t.get("value", "") for t in c.get("texts", []))
    for x in entry.get("uniProtKBCrossReferences", []):
        if x.get("database") == "GO":  # UniProt's literal cross-reference database name for the Gene Ontology
            props = {p.get("key"): p.get("value") for p in x.get("properties", [])}
            out["annotations"].append({"term_id": x.get("id", ""),
                                       "term": props.get("GoTerm", ""),
                                       "evidence": props.get("GoEvidenceType", "")})
    return out


def check_sequence(rec: dict, fasta_seq: str) -> list[str]:
    """Return a list of verification failures (empty = pass)."""
    problems = []
    if fasta_seq != rec["sequence"]:
        problems.append("FASTA sequence differs from JSON sequence")
    if str(len(fasta_seq)) != rec["length"]:
        problems.append(f"length mismatch: FASTA {len(fasta_seq)} vs JSON {rec['length']}")
    rec["non_standard_letters"] = "".join(sorted(set(fasta_seq) - STANDARD_LETTERS))   # recorded; R5 applies it in build (D59)
    computed = common.md5_text(fasta_seq)
    rec["md5_computed"] = computed
    if computed != rec["md5_uniprot"]:
        problems.append(f"MD5 mismatch: computed {computed} vs UniProt {rec['md5_uniprot']}")
    return problems


def to_section(rec: dict, tiers: list[str]) -> dict[str, str]:
    sec = {
        "tiers": ";".join(sorted(tiers)),
        "entry_name": rec["entry_name"], "entry_type": rec["entry_type"],
        "gene": rec["gene"], "protein_name": rec["protein_name"],
        "organism_id": rec["organism_id"],
        "length": rec["length"], "seq_version": rec["seq_version"],
        "entry_version": rec["entry_version"], "last_seq_update": rec["last_seq_update"],
        "md5_uniprot": rec["md5_uniprot"], "md5_computed": rec["md5_computed"],
        "md5_match": "true" if rec["md5_uniprot"] == rec["md5_computed"] else "false",
        "sequence": rec["sequence"],
        "non_standard_letters": rec.get("non_standard_letters", ""),
        "n_features_total": rec["n_features_total"],
        "n_features_captured": str(len(rec["features"])),
        "canonical_isoform_id": rec.get("canonical_isoform_id", ""),
        "n_isoforms": str(len(rec["isoforms"])),
        "tissue_specificity": rec["tissue_specificity"],
    }
    for i, f in enumerate(rec["features"]):
        p = f"feature.{i}."
        sec[p + "type"] = f["type"]
        sec[p + "start"] = f["start"]
        sec[p + "end"] = f["end"]
        sec[p + "description"] = f["description"]
        sec[p + "id"] = f["feature_id"]
    for i, iso in enumerate(rec["isoforms"]):
        p = f"isoform.{i}."
        for k in ("ids", "name", "status", "sequence_ids", "note"):
            sec[p + k] = iso.get(k, "")
        for k in ("fetched_id", "length", "md5_computed", "sequence", "non_standard_letters"):
            if k in iso:
                sec[p + k] = iso[k]
    for i, a in enumerate(rec["annotations"]):
        p = f"annotation.{i}."
        sec[p + "term_id"] = a["term_id"]
        sec[p + "term"] = a["term"]
        sec[p + "evidence"] = a["evidence"]
    return sec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", metavar="ACC", help="fetch a single accession (debugging)")
    ap.add_argument("--category", default=None, help="a non-muscle category: fetch only its pools' accessions that the store lacks (D110)")
    args = ap.parse_args(argv)
    log = common.make_logger(STAGE if args.category is None else f"{STAGE}_{args.category}")

    pool = pool_accessions(args.category)
    if args.only:
        pool = {args.only: pool.get(args.only, ["?"])}
    header_before, carried = read_store_text(common.SEQUENCES_INI)
    if args.category is not None and not args.only:
        already = {acc for acc in pool if acc in carried}
        log.info("%s: %d pool accessions, %d already in the store (not fetched again), %d to fetch",
                 args.category, len(pool), len(already), len(pool) - len(already))
        pool = {acc: t for acc, t in pool.items() if acc not in already}
    log.info("%d accessions to fetch", len(pool))

    raw_dir = common.uniprot_raw_dir(args.category)      # nothing is deleted here, ever (D61, D110)
    raw_dir.mkdir(parents=True, exist_ok=True)

    client = UniProtClient(log)
    cp = common.new_ini()
    failures: dict[str, list[str]] = {}
    for n, acc in enumerate(sorted(pool), 1):
        entry = client.entry_json(acc)
        (raw_dir / f"{acc}.json").write_text(json.dumps(entry, indent=1), encoding="utf-8", newline="\n")
        rec = extract(entry)
        if rec["accession"] != acc:
            failures[acc] = [f"primaryAccession {rec['accession']!r} != requested {acc!r} (merged/demerged entry?)"]
            log.error("%s: %s", acc, failures[acc][0])
            continue
        _, fasta_seq = client.fasta(acc)
        problems = check_sequence(rec, fasta_seq)
        # isoforms: fetch each non-displayed isoform's sequence
        for iso in rec["isoforms"]:
            ids = iso["ids"].split(";") if iso["ids"] else []
            if iso["status"] == "Displayed":
                rec["canonical_isoform_id"] = ids[0] if ids else ""
                continue
            if iso["status"] == "External" or not ids:
                continue
            try:
                _, iseq = client.fasta(ids[0])
            except Exception as exc:  # noqa: BLE001 - recorded, not fatal
                iso["fetched_id"] = ids[0]
                iso["length"] = ""
                iso["md5_computed"] = f"FETCH FAILED: {exc}"
                log.warning("%s isoform %s: fetch failed: %s", acc, ids[0], exc)
                continue
            iso["non_standard_letters"] = "".join(sorted(set(iseq) - STANDARD_LETTERS))   # recorded, not a failure (D59)
            iso["fetched_id"] = ids[0]
            iso["length"] = str(len(iseq))
            iso["md5_computed"] = common.md5_text(iseq)
            iso["sequence"] = iseq
        if problems:
            failures[acc] = problems
            log.error("%s: %s", acc, "; ".join(problems))
        else:
            log.info("[%d/%d] %s %s len=%s v%s md5 ok; features=%s/%s isoforms=%d annotations=%d",
                     n, len(pool), acc, rec["gene"], rec["length"], rec["seq_version"],
                     len(rec["features"]), rec["n_features_total"], len(rec["isoforms"]),
                     len(rec["annotations"]))
        sec = to_section(rec, pool[acc])
        if args.category is not None:
            sec["fetched"] = common.iso_now()                       # an added section says when and at which release it came
            sec["uniprot_release"] = client.release or "(header not returned)"
        cp[acc] = sec

    fetched = {acc: dict(cp[acc].items()) for acc in cp.sections()}
    keep = {name: block for name, block in carried.items() if name not in fetched}
    if args.category is None and not args.only:
        # a full run rebuilds the muscle-era store from its pool; the sections a category ADDED (they carry
        # their own `fetched` line) are carried over, so a muscle refetch cannot drop them
        keep = {name: block for name, block in keep.items() if "\nfetched " in block or "\nfetched=" in block}
    header = merge_header(header_before, args.category, client.release or "(header not returned)",
                          len(fetched), len(fetched) + len(keep), len(failures))
    common.SEQUENCES_INI.parent.mkdir(parents=True, exist_ok=True)
    common.SEQUENCES_INI.write_text(merged_store_text(header, fetched, keep), encoding="utf-8", newline="\n")
    log.info("wrote %s (%d accessions: %d fetched now, %d carried over unchanged)",
             common.SEQUENCES_INI, len(fetched) + len(keep), len(fetched), len(keep))
    if failures:
        for acc, probs in failures.items():
            log.error("VERIFICATION FAILED %s: %s", acc, "; ".join(probs))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
