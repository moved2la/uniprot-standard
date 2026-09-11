# Handoff — Step 1: Protein Set (second attempt)

**Project:** Sequence-Derived Amino Acid Standard for Human Skeletal Muscle
**Plan:** `docs/00_project_plan.md`
**Governing rule:** `PROVENANCE.md`
**Depends on:** nothing (fresh repository `uniprot-standard`)
**Feeds:** Step 2 (harvest and composition engine), Step 3 (mass fractions)
**Status:** COMPLETE 2026-09-11. Two runs: the first raised 7 multi-chain flags, closed by rule R2d (D24); the second rebuilt green with zero flags. Config regenerated independently on two machines was byte-identical.

---

## Scope (as amended 2026-09-10, D22)

Step 1 involves no literature research. It is exactly:

seed word from D3/D4 → exact-name Gene Ontology lookup by code, subtree and
neighbours recorded → UniProt query by code for human reviewed entries under that
term → fetch and verify → the field-based rules (canonical sequence; Chain-feature
processing; evidence recorded), flagging what they cannot settle → generated config →
tests → methods section describing only what the code did.

The protein set's definition is: "canonical sequences of human, reviewed UniProtKB
entries annotated to Gene Ontology term <ID> or its descendants, Gene Ontology release
<X>, UniProt release <Y>." The annotation-completeness limitation is disclosed and its
check is scheduled for Step 3.

Removed from this step:
- Citing D3's reasoning → paper step. D3 is an author decision and is admissible as such.
- Measurement cross-check and fractionation-study search → Step 3.
- Fiber-type field → dropped entirely (D17).
- Isoform selection → dropped entirely (D16); deltas are disclosed, not acted on.

## Inputs (only these)

| Input | Origin |
|---|---|
| Seed words **myofibrillar** (D3), **sarcoplasmic** (D4); lookup names (D13, D14); ontology file variant (D15) | `config/protein_set_decisions.ini`, `docs/decisions.md` |
| Gene Ontology basic release | fetched by code, unchanged, hashed |
| UniProtKB | fetched by code, release recorded |

No gene name, accession, or term ID is typed anywhere.

## How to run

```
pip install -r requirements.txt
python run.py protein-set            # resolve → enumerate → verify → build → deltas → tests
python run.py protein-set --offline  # rebuild from data/ without network (e.g. after closing flags)
```

Upload back: `data/` (all of it), `config/accessions.ini`, `config/segments.ini`,
`outputs/` (all of it). Expect the first run to stop at `build_protein_set` if any
flag is open; the log names the flags. Close them in
`config/protein_set_decisions.ini` with a D-number, rerun with `--offline`.

## Rules (see `docs/conventions.md` for the field each reads)

R0 category by exact name · R0a subtree by `is_a`/`part_of` · R0b pool by union of
per-term queries · R1 canonical sequence always · R2 Chain feature (one → range; none
→ whole; several → flag) · R3 evidence recorded, never filtered · R4 tier = pool
membership.

## Deliverables

- `PROVENANCE.md`, `README.md`, `requirements.txt`, `run.py`, `pipeline/`, `tests/`
- `config/protein_set_decisions.ini` (hand-written; seeds and closures only)
- Generated at run time: `data/gene-ontology/*`, `data/tier<N>_pool*.tsv`,
  `data/pool_queries.ini`, `data/uniprot_raw/`, `data/uniprot_verification.ini`,
  `config/accessions.ini`, `config/segments.ini`, `outputs/*.tsv`, `outputs/logs/`
- `docs/conventions.md`, `docs/decisions.md` (D1–D23), `docs/methods.md` (placeholders
  filled after the run)

## Acceptance criteria

- [ ] Each tier's term was matched by exact name; match, definition, release, subtree,
      ancestors, and neighbours are in `data/gene-ontology/`. No related-term search.
- [ ] Each pool was produced by recorded UniProt queries; term-level-vs-union check and
      neighbour sensitivity counts are recorded.
- [ ] Every pool member fetched and MD5-verified; run exits zero.
- [ ] Every processing choice traces to the `Chain` feature or to a closed flag with a
      D-number. No isoform choice exists.
- [ ] `config/protein_set_decisions.ini` contains no protein that is not a flag closure.
- [ ] `pytest` passes, including generated-config-is-current and no-open-flags.
- [ ] `docs/methods.md` "Protein set definition" has its placeholders filled from the
      run records and no uncited sentence.
- [ ] No mass fraction, abundance, or literature-derived number appears anywhere.

## Explicitly out of scope

Mass fractions and weighting (Step 3); composition beyond the isoform/processing
deltas (Step 2); collagen/ECM; literature of any kind.

## Deferred / raised

- Pool completeness check against Layer B abundance ranking — Step 3.
- Citations for D3's reasoning (myofibrillar-synthesis methodology; training-response
  literature) — paper step.
- Overlap between Tier 1 and Tier 2 pools (`outputs/pool_overlap.tsv`) — handling rule
  decided in Step 3 when the list exists (within-tier normalisation may make it moot).
- Isoform and processing bounds (`outputs/isoform_bound.tsv`,
  `outputs/processing_bound.tsv`) — revisit in Step 3 only if mass-weighting makes one
  material.
- Any pool member returned by the UniProt query but lacking a subtree annotation in its
  entry cross-references (`pool_member_without_subtree_annotation` flags) — indicates
  the search index and the entry JSON disagree; close by recording, or investigate
  in Step 3.

## Handoff to Step 2

- Pool sizes: Tier 1 = 242, Tier 2 = 87, both = 26, unique accessions = 303.
- Gene Ontology release `releases/2026-07-26` (sha256 in `data/gene-ontology/source.ini`);
  UniProt release 2026_03; fetched 2026-09-11.
- Every accession uses its canonical sequence (R1). No isoform IDs anywhere.
- `segments.ini`: 303 master runs (295 R2a, 1 R2b, 7 R2d); 93 runs with
  `in_master_molecule = false` (91 N-terminal, 2 C-terminal, 0 internal) across 92 entries.
- Flags: none open, none raised in the final build.
- **Sequences are already on disk, verified.** `data/uniprot_verification.ini` holds every
  canonical sequence with UniProt's MD5 and the computed MD5. Step 2's harvest is
  therefore already done for this set; the composition engine reads
  `uniprot_verification.ini` + `segments.ini` and needs no fetcher of its own. Re-fetching
  is only required if the pool is regenerated at a new release.
- Two isoform FASTAs returned HTTP 404 (recorded in the verification file); irrelevant to
  Step 2 since only canonical sequences are used.
