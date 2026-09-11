"""Enumerate the candidate pool for each tier by querying UniProt, by code.

Inputs
  data/gene-ontology/resolved_terms.ini     (from resolve_ontology_term.py)
  data/gene-ontology/tier<N>_subtree.tsv
  data/gene-ontology/tier<N>_neighbours.tsv

Definition of the pool (per tier)
  Human (organism_id:9606), reviewed (reviewed:true) UniProtKB entries annotated to
  the tier's term or to any term in its descendant subtree. The pool is computed as
  the UNION of one query per subtree term, so its definition does not depend on how
  UniProt expands a term. UniProt's own term-level query is run as a check and the
  two are compared; a difference is recorded, never silently absorbed.

Sensitivity
  The same term-level query is run for each immediate neighbour of the tier's term
  (next broader, next narrower) and the set differences against the pool are
  recorded. This is the sensitivity of the pool to the category boundary.

Outputs
  data/tier<N>_pool.tsv
  data/tier<N>_pool_sensitivity.tsv
  data/pool_queries.ini
  outputs/pool_overlap.tsv         (accessions appearing in more than one tier)
"""

from __future__ import annotations

import argparse
import sys

from pipeline import common
from pipeline.flags import append_flag, reset_flags
from pipeline.uniprot_client import UniProtClient

FIELDS = ["accession", "gene_primary", "protein_name", "length", "reviewed", "organism_id"]
STAGE = "enumerate_pool"


def term_query(term_id: str) -> str:
    """UniProt query-field syntax uses the numeric part of a Gene Ontology ID."""
    num = term_id.split(":", 1)[1]
    return f"(organism_id:9606) AND (reviewed:true) AND (go:{num})"


def _col(row: dict[str, str], *candidates: str) -> str:
    for c in candidates:
        if c in row:
            return row[c]
    return ""


def _accession(row: dict[str, str]) -> str:
    return _col(row, "Entry", "accession") or next(iter(row.values()))


def enumerate_tier(client: UniProtClient, tier: str, term_id: str, term_name: str,
                   subtree: list[dict[str, str]], neighbours: list[dict[str, str]], log):
    # --- definitional pool: union of one query per subtree term ---------------
    members: dict[str, dict[str, str]] = {}
    returned_by: dict[str, list[str]] = {}
    per_term_counts: list[list] = []
    for t in subtree:
        tid = t["term_id"]
        rows = client.search_tsv(term_query(tid), FIELDS)
        per_term_counts.append(["subtree", tid, t["name"], t["relation_to_parent"], t["depth"], len(rows)])
        for r in rows:
            acc = _accession(r)
            members.setdefault(acc, r)
            returned_by.setdefault(acc, []).append(tid)
    pool = set(members)
    log.info("tier %s: union over %d subtree terms -> %d entries", tier, len(subtree), len(pool))

    # --- check: UniProt's own expansion of the tier term ------------------------
    top_rows = client.search_tsv(term_query(term_id), FIELDS)
    top = {_accession(r) for r in top_rows}
    only_in_top = sorted(top - pool)
    only_in_union = sorted(pool - top)
    expansion_matches = not only_in_top and not only_in_union
    log.info("tier %s: term-level query -> %d entries; matches union: %s "
             "(only_in_term_query=%d, only_in_union=%d)",
             tier, len(top), expansion_matches, len(only_in_top), len(only_in_union))
    if not expansion_matches:
        append_flag(common.FLAGS_TSV, stage=STAGE, subject=f"tier.{tier}.expansion",
                    rule="descendant_expansion_check",
                    detail=(f"term-level query and per-term union differ: "
                            f"only_in_term_query={only_in_top} only_in_union={only_in_union}"))

    # --- sensitivity: immediate neighbours ---------------------------------------
    sens_rows: list[list] = []
    sens_rows.append(["tier_term", term_id, term_name, "", "term-level query", len(top),
                      len(only_in_top), len(only_in_union), ";".join(only_in_top), ";".join(only_in_union)])
    for nb in neighbours:
        tid = nb["term_id"]
        rows = client.search_tsv(term_query(tid), FIELDS)
        s = {_accession(r) for r in rows}
        extra = sorted(s - pool)
        missing = sorted(pool - s)
        sens_rows.append([nb["direction"], tid, nb["name"], nb["relation"], "term-level query",
                          len(rows), len(extra), len(missing), ";".join(extra), ";".join(missing)])
        log.info("tier %s: %s neighbour %s (%s): %d entries, %d not in pool, %d of pool absent",
                 tier, nb["direction"], tid, nb["name"], len(rows), len(extra), len(missing))

    # --- write --------------------------------------------------------------------
    pool_rows = []
    for acc in sorted(pool):
        r = members[acc]
        pool_rows.append([
            acc,
            _col(r, "Gene Names (primary)", "gene_primary"),
            _col(r, "Protein names", "protein_name"),
            _col(r, "Length", "length"),
            _col(r, "Reviewed", "reviewed"),
            _col(r, "Organism (ID)", "organism_id"),
            ";".join(returned_by[acc]),
            "true" if acc in top else "false",
        ])
    common.write_tsv(common.DATA_DIR / f"tier{tier}_pool.tsv",
                     ["accession", "gene_primary", "protein_name", "length", "reviewed", "organism_id",
                      "subtree_terms_returning_entry", "in_term_level_query"],
                     pool_rows)
    common.write_tsv(common.DATA_DIR / f"tier{tier}_pool_sensitivity.tsv",
                     ["query_label", "term_id", "term_name", "relation", "query_kind", "n_entries",
                      "n_not_in_pool", "n_pool_not_returned", "accessions_not_in_pool", "pool_accessions_not_returned"],
                     sens_rows)
    common.write_tsv(common.DATA_DIR / f"tier{tier}_pool_per_term_counts.tsv",
                     ["kind", "term_id", "term_name", "relation_to_parent", "depth", "n_entries"],
                     per_term_counts)
    return pool, len(top), expansion_matches


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    log = common.make_logger(STAGE)
    reset_flags(common.FLAGS_TSV, STAGE)
    resolved = common.read_ini(common.RESOLVED_TERMS_INI)
    client = UniProtClient(log)

    queries = common.new_ini()
    queries["run"] = {"date": common.iso_now(), "query_template": term_query("GO:NNNNNNN"),
                      "fields": ",".join(FIELDS), "endpoint": "https://rest.uniprot.org/uniprotkb/search"}
    pools: dict[str, set[str]] = {}
    for section in resolved.sections():
        if not section.startswith("tier."):
            continue
        tier = section.split(".", 1)[1]
        term_id = resolved[section]["term_id"]
        term_name = resolved[section]["term_name"]
        subtree = common.read_tsv(common.ONTOLOGY_DIR / f"tier{tier}_subtree.tsv")
        neighbours = common.read_tsv(common.ONTOLOGY_DIR / f"tier{tier}_neighbours.tsv")
        pool, n_top, ok = enumerate_tier(client, tier, term_id, term_name, subtree, neighbours, log)
        pools[tier] = pool
        queries[section] = {"term_id": term_id, "term_name": term_name,
                            "n_subtree_terms": str(len(subtree)), "pool_size": str(len(pool)),
                            "term_level_query_size": str(n_top),
                            "descendant_expansion_matches_union": "true" if ok else "false",
                            "pool_file": f"tier{tier}_pool.tsv"}
    queries["run"]["uniprot_release"] = client.release or "(header not returned)"
    common.write_ini(queries, common.POOL_QUERIES_INI, [
        "pool_queries.ini -- GENERATED by enumerate_pool.py. DO NOT EDIT BY HAND.",
        "Record of the UniProt queries that produced data/tier<N>_pool.tsv.",
    ])

    # overlap between tiers
    overlap_rows = []
    tiers = sorted(pools)
    for i, a in enumerate(tiers):
        for b in tiers[i + 1:]:
            for acc in sorted(pools[a] & pools[b]):
                overlap_rows.append([acc, a, b])
    common.write_tsv(common.OUTPUTS_DIR / "pool_overlap.tsv", ["accession", "tier_a", "tier_b"], overlap_rows)
    log.info("tier overlap: %d accessions in more than one tier", len(overlap_rows))
    log.info("UniProt release: %s", client.release)
    return 0


if __name__ == "__main__":
    sys.exit(main())
