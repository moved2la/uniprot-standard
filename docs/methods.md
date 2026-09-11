# Methods

Written incrementally, one section per step. Every sentence either cites a retrievable
source or describes a computation performed in this repository (`PROVENANCE.md`).
Numbers below are from the run of 2026-09-11 (ontology retrieved 06:22 UTC, UniProt fetched 06:22–06:36 UTC, config rebuilt 15:10 UTC).

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
`data/gene-ontology/source.ini` (`data-version: releases/2026-07-26`). The noun form
of each seed word (D13, D14) was matched against term names by exact string equality;
the run stops if the match is not exactly one non-obsolete term. The matches were:

| Tier | Seed word | Lookup name | Term ID | Term definition (verbatim from the ontology file) |
|---|---|---|---|---|
| 1 | myofibrillar (D3) | myofibril (D13) | `GO:0030016` | `The contractile element of skeletal and cardiac muscle; a long, highly organized bundle of actin, myosin, and other proteins that contracts by a sliding filament mechanism. [ISBN:0815316194]` |
| 2 | sarcoplasmic (D4) | sarcoplasm (D14) | `GO:0016528` | `The cytoplasm of a muscle cell; includes the sarcoplasmic reticulum. [ISBN:0198547684]` |

(`data/gene-ontology/resolved_terms.ini`.) The descendant subtree of each term was
enumerated over the `is_a` and `part_of` relations and written to
`data/gene-ontology/tier<N>_subtree.tsv` (18 and 13 terms including the root). The term's ancestors to the root and
its immediate broader and narrower terms are in `tier<N>_ancestors.tsv` and
`tier<N>_neighbours.tsv`.

### Pool

For each subtree term, UniProtKB was queried by `pipeline/enumerate_pool.py` with
`(organism_id:9606) AND (reviewed:true) AND (go:<term>)` at UniProt release
2026_03 (`data/pool_queries.ini`). The pool of a tier is the union of these
per-term results (D20): 242 entries for Tier 1 and 87 for Tier 2 (`data/tier<N>_pool.tsv`). The term-level query alone returned
242 and 87 entries; it equalled the
union for both tiers (`descendant_expansion_matches_union = true` in
`data/pool_queries.ini`), so UniProt's expansion of a term to its descendants is
confirmed by the run rather than assumed. The same query run for each immediate neighbour of the tier term
gave the counts and set differences in `data/tier<N>_pool_sensitivity.tsv`: for Tier 1,
the next broader term (`contractile muscle fiber`, GO:0043292) returns 252 entries, 10
more than the pool; the largest child (`sarcomere`, GO:0030017) returns 221; the child
`skeletal muscle myofibril` (GO:0098723) returns 0 and `cardiac myofibril`
(GO:0097512) returns 7. For Tier 2, the next broader term (`cytoplasm`, GO:0005737)
returns 11,855; the child `sarcoplasmic reticulum` (GO:0016529) returns 77 of the 87.
26 entries belong to both pools (`outputs/pool_overlap.tsv`).

The protein set is therefore defined as: the canonical sequences of human, reviewed
UniProtKB entries annotated to `GO:0030016` or its `is_a`/`part_of` descendants
(Tier 1), and likewise for `GO:0016528` (Tier 2), at Gene Ontology release
`releases/2026-07-26` and UniProt release 2026_03.

### Verification

Each pool entry was fetched by `pipeline/verify_accessions.py` as entry JSON and FASTA.
The MD5 of the FASTA sequence was compared with the MD5 published in the entry JSON,
the FASTA and JSON sequences were compared character-for-character, lengths were
compared, and only the twenty standard letters were accepted; any failure stops the
run. 303 entries passed (`data/uniprot_verification.ini`, header). The
entry JSON was saved unchanged to `data/uniprot_raw/`.

### Sequence and processing

Every entry is represented by its canonical sequence (D16). The mature-chain range is
taken from the UniProt feature table, type `Chain`. Exactly one such feature gives the
range (rule R2a; 295 entries); residues outside it are marked `in_master_molecule =
false` in `config/segments.ini` (D19). One entry has no `Chain` feature and uses the
whole sequence (rule R2b, build log). Seven entries have more than one `Chain` feature;
for these the master molecule is the union of the `Chain` ranges (rule R2d, D24), and
they are listed with their chains in `outputs/multi_chain_entries.tsv`. No `Chain`
feature had a non-exact position, so no flag was raised (`outputs/flags.tsv` is empty).
In total 93 residue runs are marked `in_master_molecule = false` (91 N-terminal, 2
C-terminal) across 92 entries.

### Evidence

The Gene Ontology evidence code of every annotation that placed an entry in its tier's
subtree is recorded in `config/accessions.ini` and summarised in
`outputs/evidence_summary.tsv` (D18). Of the 242 Tier 1 entries, 165 carry at least one
subtree annotation with a code other than IEA (Inferred from Electronic Annotation) and
77 carry only IEA annotations within the subtree; of the 87 Tier 2 entries, 64 and 23
respectively (`outputs/accession_evidence.tsv`). The per-code counts (IBA, IC, IDA,
IEA, IMP, IPI, ISS, NAS, TAS for Tier 1; EXP, IBA, IDA, IEA, IMP, IPI, ISS, TAS for
Tier 2) are in `outputs/evidence_summary.tsv`. Code names are the Gene Ontology
consortium's.

### Disclosed bounds

`pipeline/isoform_processing_deltas.py` computed, from residue counts alone, the
difference in each amino acid's fraction between every non-canonical isoform whose
sequence UniProt provides and the canonical sequence, and between the mature-chain
range and the full canonical sequence. 198 of the 303 entries have at least one such isoform
(662 isoform rows). Across the set, the largest such difference for any single amino acid
is 0.138 for isoform choice (glycine; one short isoform of one entry;
`outputs/isoform_bound.tsv`) and 0.046 for chain processing (serine; one entry with a
long propeptide; `outputs/processing_bound.tsv`), both as fractions of residues. No
decision is attached to these numbers; they bound what a different isoform or
processing choice could do to the composition of any one entry before mass weighting.
The per-isoform and per-entry values are in `outputs/isoform_deltas.tsv` and
`outputs/processing_deltas.tsv`; the two isoform sequences UniProt did not serve
(HTTP 404, recorded in `data/uniprot_verification.ini`) are absent from the isoform table.

### Limitation

The pool is an annotation set. Its completeness relative to measurement is checked in
Step 3 against the Layer B abundance ranking (D22).
