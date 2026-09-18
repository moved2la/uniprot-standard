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

Measured tier (D56)
  A tier whose decisions section says `definition = measured_remainder` is not a Gene
  Ontology category. Its pool is every gene quantified in the primary Layer B dataset
  (the [source.*] named in the decisions section, read through its [columns.*] map in
  config/mass_fraction_decisions.ini and data/literature/manifest.ini) that maps by
  exact gene symbol to exactly one reviewed human UniProt entry, is not already in an
  ontology tier's pool, and is not in the cited contaminant list. Every gene's outcome
  is written to data/tier<N>_gene_mapping.tsv. Rules M1-M5 (docs/conventions.md).
  The ontology tiers are enumerated first, so "already in an ontology tier" is known.

Outputs
  data/tier<N>_pool.tsv
  data/tier<N>_pool_sensitivity.tsv
  data/pool_queries.ini
  outputs/pool_overlap.tsv         (accessions appearing in more than one tier)
"""

from __future__ import annotations

import argparse
import sys

import requests

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


# ----------------------------------------------------------------------------
# measured remainder (D56)
# ----------------------------------------------------------------------------

MEASURED_FIELDS = ["accession", "gene_primary", "gene_synonym", "protein_name", "length", "reviewed", "organism_id"]
GENE_BATCH = 50


def dataset_gene_cells(source_id: str, file_key: str, log) -> list[str]:
    """Every value of the identity column of the primary dataset, as read (rules B0, I1)."""
    import io
    from openpyxl import load_workbook
    dec = common.read_ini(common.CONFIG_DIR / "mass_fraction_decisions.ini")
    manifest = common.read_ini(common.DATA_DIR / "literature" / "manifest.ini")
    cols = dec[f"columns.{source_id}.{file_key}"]
    m = manifest[f"file.{source_id}.{file_key}"]
    path = common.REPO_ROOT / m["path"]
    got = common.sha256_file(path)
    if got != m["sha256"].lower():
        raise SystemExit(f"[STOP] {m['path']} sha256 {got} != manifest {m['sha256']}")
    wb = load_workbook(io.BytesIO(path.read_bytes()), read_only=True, data_only=True)
    ws = wb[cols["sheet"]]
    header_row = int(cols["header_row"])
    cells = []
    for r_index, row in enumerate(ws.iter_rows(values_only=True), 1):
        if r_index < header_row:
            continue
        if r_index == header_row:
            headers = ["" if v is None else str(v).strip() for v in row]
            hits = [j for j, h in enumerate(headers) if h == cols["identity_column"].strip()]
            if len(hits) != 1:
                raise SystemExit(f"[STOP] identity column {cols['identity_column']!r} matched {len(hits)} headers")
            j = hits[0]
            continue
        v = row[j] if j < len(row) else None
        if v is not None and str(v).strip():
            cells.append(str(v).strip())
    wb.close()
    log.info("measured tier: %d identity cells read from %s sha256 %s", len(cells), m["path"], got)
    return cells


def contaminant_accessions(source_id: str, log) -> set[str]:
    """Accessions in the cited contaminant FASTA (gzip or plain). The header's first
    token, with a leading 'CON__' or 'sp|'/'tr|' prefix removed, is the accession."""
    import gzip
    manifest = common.read_ini(common.DATA_DIR / "literature" / "manifest.ini")
    m = manifest[f"file.{source_id}.file.1"]
    path = common.REPO_ROOT / m["path"]
    got = common.sha256_file(path)
    if got != m["sha256"].lower():
        raise SystemExit(f"[STOP] {m['path']} sha256 {got} != manifest {m['sha256']}")
    data = path.read_bytes()
    text = (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode("utf-8", errors="replace")
    accs = set()
    for line in text.splitlines():
        if not line.startswith(">"):
            continue
        tok = line[1:].split()[0] if line[1:].split() else ""
        tok = tok.replace("CON__", "")
        if "|" in tok:
            parts = tok.split("|")
            tok = parts[1] if len(parts) > 1 else parts[0]
        accs.add(tok.split("-")[0])
    log.info("measured tier: %d contaminant accessions parsed from %s sha256 %s", len(accs), m["path"], got)
    return accs


def enumerate_measured_tier(client: UniProtClient, tier: str, sec, ontology_pools: dict[str, set[str]], log):
    """Rules
    M1  genes = every identity cell of the primary dataset, split on ';', stripped, unique.
    M2  each gene is looked up with (organism_id:9606) AND (reviewed:true) AND (gene_exact:<gene>),
        in batches; a returned entry belongs to a gene when the gene equals its primary
        symbol or one of its synonyms (case-sensitive, as UniProt reports them).
    M3  exactly one entry -> candidate. Zero -> unmapped, listed. Several -> ambiguous, listed, excluded.
        A symbol UniProt's parser rejects even when quoted -> query_rejected, listed.
    M4  a candidate already in an ontology tier's pool -> in_tier<N>, not a member.
    M5  a candidate in the contaminant list -> contaminant, not a member. Everything else -> member.
    Nothing is named; nothing is filtered by abundance."""
    cells = dataset_gene_cells(sec["dataset_source"], sec["dataset_file"], log)
    genes = sorted({g.strip() for c in cells for g in c.split(";") if g.strip()})
    contaminants = contaminant_accessions(sec["contaminant_source"], log)
    in_ontology = {acc: t for t, accs in ontology_pools.items() for acc in accs}
    found: dict[str, list[dict[str, str]]] = {g: [] for g in genes}
    rejected: set[str] = set()

    def quoted(g: str) -> str:
        return 'gene_exact:"' + g.replace('"', '') + '"'      # quoted: symbols may carry ':' '.' '-' (M2)

    def attribute(rows, batch):
        for r in rows:
            primary = _col(r, "Gene Names (primary)", "gene_primary").strip()
            syn = _col(r, "Gene Names (synonym)", "gene_synonym")
            names = {primary} | {s.strip() for s in syn.split() if s.strip()}
            for g in batch:
                if g in names:
                    found[g].append(dict(r, _matched_by="primary" if g == primary else "synonym"))

    for i in range(0, len(genes), GENE_BATCH):
        batch = genes[i:i + GENE_BATCH]
        q = "(organism_id:9606) AND (reviewed:true) AND (" + " OR ".join(quoted(g) for g in batch) + ")"
        try:
            rows = client.search_tsv(q, MEASURED_FIELDS)
            attribute(rows, batch)
            n = len(rows)
        except requests.exceptions.HTTPError as exc:                # a symbol UniProt cannot parse: retry one by one
            log.warning("measured tier: batch %d-%d rejected (%s); querying each gene alone", i + 1, i + len(batch), exc)
            n = 0
            for g in batch:
                try:
                    rows = client.search_tsv("(organism_id:9606) AND (reviewed:true) AND " + quoted(g), MEASURED_FIELDS)
                    attribute(rows, [g])
                    n += len(rows)
                except requests.exceptions.HTTPError as exc2:
                    rejected.add(g)
                    log.warning("measured tier: query rejected for %r (%s)", g, exc2)
        log.info("measured tier: genes %d-%d of %d queried; %d rows", i + 1, min(i + GENE_BATCH, len(genes)), len(genes), n)
    members: dict[str, dict[str, str]] = {}
    mapping_rows = []
    counts = {"member": 0, "unmapped": 0, "ambiguous": 0, "contaminant": 0, "query_rejected": 0}
    for g in genes:
        rows = found[g]
        accs = sorted({_accession(r) for r in rows})
        if g in rejected:
            outcome = "query_rejected"
        elif not accs:
            outcome = "unmapped"
        elif len(accs) > 1:
            outcome = "ambiguous"
        else:
            acc = accs[0]
            if acc in in_ontology:
                outcome = f"in_tier{in_ontology[acc]}"
            elif acc in contaminants:
                outcome = "contaminant"
            else:
                outcome = "member"
                members[acc] = rows[0]
        counts[outcome] = counts.get(outcome, 0) + 1
        mapping_rows.append([g, outcome, ";".join(accs), rows[0]["_matched_by"] if len(accs) == 1 else "",
                             _col(rows[0], "Gene Names (primary)", "gene_primary") if len(accs) == 1 else ""])
    common.write_tsv(common.DATA_DIR / f"tier{tier}_gene_mapping.tsv",
                     ["gene", "outcome", "accessions", "matched_by", "uniprot_primary_symbol"], mapping_rows)
    pool_rows = []
    for acc in sorted(members):
        r = members[acc]
        pool_rows.append([acc, _col(r, "Gene Names (primary)", "gene_primary"), _col(r, "Protein names", "protein_name"),
                          _col(r, "Length", "length"), _col(r, "Reviewed", "reviewed"), _col(r, "Organism (ID)", "organism_id"),
                          "", ""])
    common.write_tsv(common.DATA_DIR / f"tier{tier}_pool.tsv",
                     ["accession", "gene_primary", "protein_name", "length", "reviewed", "organism_id",
                      "subtree_terms_returning_entry", "in_term_level_query"], pool_rows)
    log.info("measured tier %s: %d genes -> %s", tier, len(genes), counts)
    return set(members), counts


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
        if resolved[section].get("definition") == "measured_remainder":
            continue                                    # after every ontology tier (M4)
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
    for section in resolved.sections():
        if section.startswith("tier.") and resolved[section].get("definition") == "measured_remainder":
            tier = section.split(".", 1)[1]
            pool, counts = enumerate_measured_tier(client, tier, resolved[section], pools, log)
            pools[tier] = pool
            queries[section] = {"definition": "measured_remainder", "dataset_source": resolved[section]["dataset_source"],
                                "dataset_file": resolved[section]["dataset_file"],
                                "contaminant_source": resolved[section]["contaminant_source"],
                                "query_template": "(organism_id:9606) AND (reviewed:true) AND (gene_exact:G1 OR ... )",
                                "genes": str(sum(counts.values())), "pool_size": str(len(pool)),
                                **{f"n_{k}": str(v) for k, v in counts.items()},
                                "pool_file": f"tier{tier}_pool.tsv", "mapping_file": f"tier{tier}_gene_mapping.tsv"}
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
    common.write_tsv(common.PROTEIN_SET_OUT_DIR / "pool_overlap.tsv", ["accession", "tier_a", "tier_b"], overlap_rows)
    log.info("tier overlap: %d accessions in more than one tier", len(overlap_rows))
    log.info("UniProt release: %s", client.release)
    return 0


if __name__ == "__main__":
    sys.exit(main())
