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
  data/uniprot_verification.ini    -- everything build_protein_set.py reads
  data/uniprot_raw/<acc>.json      -- raw entries (folder emptied at start of run)
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys

from pipeline import common
from pipeline.uniprot_client import UniProtClient

STAGE = "verify_accessions"
CAPTURED_FEATURE_TYPES = ("Chain", "Initiator methionine", "Signal", "Propeptide",
                          "Transit peptide", "Peptide", "Alternative sequence")
STANDARD_LETTERS = set(common.AMINO_ACIDS)


def pool_accessions() -> dict[str, list[str]]:
    """accession -> [tiers it belongs to], from every data/tier<N>_pool.tsv."""
    out: dict[str, list[str]] = {}
    for path in sorted(common.DATA_DIR.glob("tier*_pool.tsv")):
        tier = path.name[len("tier"):].split("_", 1)[0]
        for row in common.read_tsv(path):
            out.setdefault(row["accession"], []).append(tier)
    if not out:
        raise FileNotFoundError("no data/tier<N>_pool.tsv files found; run enumerate_pool first")
    return out


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
    bad = sorted(set(fasta_seq) - STANDARD_LETTERS)
    if bad:
        problems.append(f"non-standard letters in sequence: {bad}")
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
        for k in ("fetched_id", "length", "md5_computed", "sequence"):
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
    args = ap.parse_args(argv)
    log = common.make_logger(STAGE)

    pool = pool_accessions()
    if args.only:
        pool = {args.only: pool.get(args.only, ["?"])}
    log.info("%d accessions to fetch", len(pool))

    if not args.only and common.UNIPROT_RAW_DIR.exists():
        shutil.rmtree(common.UNIPROT_RAW_DIR)
    common.UNIPROT_RAW_DIR.mkdir(parents=True, exist_ok=True)

    client = UniProtClient(log)
    cp = common.new_ini()
    failures: dict[str, list[str]] = {}
    for n, acc in enumerate(sorted(pool), 1):
        entry = client.entry_json(acc)
        (common.UNIPROT_RAW_DIR / f"{acc}.json").write_text(json.dumps(entry, indent=1), encoding="utf-8", newline="\n")
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
            bad = sorted(set(iseq) - STANDARD_LETTERS)
            if bad:
                problems.append(f"isoform {ids[0]}: non-standard letters {bad}")
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
        cp[acc] = to_section(rec, pool[acc])

    common.write_ini(cp, common.VERIFICATION_INI, [
        "uniprot_verification.ini -- GENERATED by verify_accessions.py. DO NOT EDIT BY HAND.",
        f"fetched         = {common.iso_now()}",
        f"uniprot_release = {client.release or '(header not returned)'}",
        f"source          = https://rest.uniprot.org/uniprotkb/",
        f"captured_feature_types = {'; '.join(CAPTURED_FEATURE_TYPES)}",
        f"accessions      = {len(cp.sections())}",
        f"failures        = {len(failures)}",
    ])
    log.info("wrote %s (%d accessions)", common.VERIFICATION_INI, len(cp.sections()))
    if failures:
        for acc, probs in failures.items():
            log.error("VERIFICATION FAILED %s: %s", acc, "; ".join(probs))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
