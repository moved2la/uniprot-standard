# Methods

Written incrementally, one section per step. Every sentence either cites a retrievable
source or describes a computation performed in this repository (`PROVENANCE.md`).
Angle-bracket placeholders are filled from the run records named beside them.

## Protein set definition

### Scope

The primary standard is the myofibrillar protein fraction of human skeletal muscle
(author decision D3, `docs/decisions.md`); the sarcoplasmic fraction is computed
separately as a sensitivity analysis (D4). Citations for the reasoning behind D3 are
attached in the paper step.

### Category

The seed words were resolved to Gene Ontology terms by `pipeline/resolve_ontology_term.py`
as follows. The basic release of the ontology was downloaded from its permanent URL
(D15); the file's version line, retrieval time, and SHA-256 are in
`data/gene-ontology/source.ini` (`data-version: <ONTOLOGY_DATA_VERSION>`). The noun form
of each seed word (D13, D14) was matched against term names by exact string equality;
the run stops if the match is not exactly one non-obsolete term. The matches were:

| Tier | Seed word | Lookup name | Term ID | Term definition (verbatim from the ontology file) |
|---|---|---|---|---|
| 1 | myofibrillar (D3) | myofibril (D13) | `<TIER1_TERM_ID>` | `<TIER1_DEFINITION>` |
| 2 | sarcoplasmic (D4) | sarcoplasm (D14) | `<TIER2_TERM_ID>` | `<TIER2_DEFINITION>` |

(`data/gene-ontology/resolved_terms.ini`.) The descendant subtree of each term was
enumerated over the `is_a` and `part_of` relations and written to
`data/gene-ontology/tier<N>_subtree.tsv` (`<TIER1_SUBTREE_SIZE>` and
`<TIER2_SUBTREE_SIZE>` terms including the root). The term's ancestors to the root and
its immediate broader and narrower terms are in `tier<N>_ancestors.tsv` and
`tier<N>_neighbours.tsv`.

### Pool

For each subtree term, UniProtKB was queried by `pipeline/enumerate_pool.py` with
`(organism_id:9606) AND (reviewed:true) AND (go:<term>)` at UniProt release
`<UNIPROT_RELEASE>` (`data/pool_queries.ini`). The pool of a tier is the union of these
per-term results (D20): `<TIER1_POOL_SIZE>` entries for Tier 1 and `<TIER2_POOL_SIZE>`
for Tier 2 (`data/tier<N>_pool.tsv`). The term-level query alone returned
`<TIER1_TERM_LEVEL_SIZE>` and `<TIER2_TERM_LEVEL_SIZE>` entries; whether it equalled
the union is recorded per tier in `data/pool_queries.ini`
(`descendant_expansion_matches_union`), and any difference is listed in
`outputs/flags.tsv`. The same query run for each immediate neighbour of the tier term
gave the counts and set differences in `data/tier<N>_pool_sensitivity.tsv`.
`<N_OVERLAP>` entries belong to both pools (`outputs/pool_overlap.tsv`).

The protein set is therefore defined as: the canonical sequences of human, reviewed
UniProtKB entries annotated to `<TIER1_TERM_ID>` or its `is_a`/`part_of` descendants
(Tier 1), and likewise for `<TIER2_TERM_ID>` (Tier 2), at Gene Ontology release
`<ONTOLOGY_DATA_VERSION>` and UniProt release `<UNIPROT_RELEASE>`.

### Verification

Each pool entry was fetched by `pipeline/verify_accessions.py` as entry JSON and FASTA.
The MD5 of the FASTA sequence was compared with the MD5 published in the entry JSON,
the FASTA and JSON sequences were compared character-for-character, lengths were
compared, and only the twenty standard letters were accepted; any failure stops the
run. `<N_ACCESSIONS>` entries passed (`data/uniprot_verification.ini`, header). The
entry JSON was saved unchanged to `data/uniprot_raw/`.

### Sequence and processing

Every entry is represented by its canonical sequence (D16). The mature-chain range is
taken from the UniProt feature table, type `Chain`: exactly one such feature gives the
range; residues outside it are marked `in_master_molecule = false` in
`config/segments.ini` (D19). Entries with no `Chain` feature use the whole sequence
(`<N_R2B>` entries, listed in the build log). Entries with more than one `Chain`
feature or a non-exact position were flagged and closed by decision
(`<N_R2C>` entries; `outputs/flags.tsv`, `docs/decisions.md`).

### Evidence

The Gene Ontology evidence code of every annotation that placed an entry in its tier's
subtree is recorded in `config/accessions.ini` and summarised in
`outputs/evidence_summary.tsv` (D18). `<EVIDENCE_SUMMARY_SENTENCE>` — filled from that
table, e.g. "Of the Tier 1 pool, N entries carry at least one experimental-type code
and M carry only electronically inferred codes", with the code classes named as the
consortium names them.

### Disclosed bounds

`pipeline/isoform_processing_deltas.py` computed, from residue counts alone, the
difference in each amino acid's fraction between every non-canonical isoform whose
sequence UniProt provides and the canonical sequence, and between the mature-chain
range and the full canonical sequence. Across the set, the largest such difference for
any single amino acid is `<ISOFORM_BOUND>` (isoform choice; `outputs/isoform_bound.tsv`)
and `<PROCESSING_BOUND>` (chain processing; `outputs/processing_bound.tsv`). No
decision is attached to these numbers; they bound what a different isoform or
processing choice could do to the composition of any one entry before mass weighting.

### Limitation

The pool is an annotation set. Its completeness relative to measurement is checked in
Step 3 against the Layer B abundance ranking (D22).
