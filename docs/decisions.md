# Decisions Log

Every author decision, numbered, dated, with its reason. A decision may cite sources
for its reason; the decision itself is the author's (`PROVENANCE.md`, origin type 1).
Decisions D1–D12 are carried from `docs/00_project_plan.md` §4 (planning thread,
2026-09-09). Nothing else from any earlier repository is inherited.

| # | Date | Decision | Reason |
|---|------|----------|--------|
| D1 | 2026-09-09 | The standard is built first; the Match Rate module comes after. | The standard is the input to the match calculation. |
| D2 | 2026-09-09 | Two-layer model: Layer A exact (sequence-derived), Layer B cited (tissue mass fractions). | UniProt gives composition, not abundance. Naming the layers makes the uncertainty explicit. |
| D3 | 2026-09-09 | The primary standard is the **myofibrillar** fraction. | Original scoping concept was muscle fiber; muscle-protein-synthesis literature measures myofibrillar synthesis separately; resistance training primarily drives the myofibrillar fraction. *Citations for this reasoning are attached at the paper step, not here.* |
| D4 | 2026-09-09 | The **sarcoplasmic** fraction is computed and reported as a sensitivity analysis, not folded into the primary standard. | Needed to compare against whole-tissue laboratory standards and to justify the scope with a number rather than by assertion. |
| D5 | 2026-09-09 | *(superseded by D20; retained for history)* Certain named sarcomeric components were to be included in Tier 1 from the first run. | Under the current design the pool is computed from the category, so no protein is placed by hand; the named list is not carried into this repository. |
| D6 | 2026-09-09 | Fiber-type-specific entries are used, not proxies. | Satisfied by construction: the computed pool contains each paralog as its own entry. |
| D7 | 2026-09-09 | Residue counts are the stored ground truth; both mass conventions are derived. | Residue-in-chain mass ≠ free amino acid mass. Match Rate must use the free-AA convention. |
| D8 | 2026-09-09 | Segment-flag architecture (`in_master_molecule`) retained from the collagen pipeline. | Proven architecture; the collagen tier reuses it. |
| D9 | 2026-09-09 | Layer B primary source = quantitative human proteomics; classical fractionation is the cross-check. | Human, fiber-type-resolved abundance; agreement between methods is itself a validation result. |
| D10 | 2026-09-09 | Literature is researched in-thread where a step needs it, not deferred. | Public scrutiny standard. |
| D11 | 2026-09-09 | The accession set is settled before mass fractions. | Cannot weight what has not been enumerated. |
| D12 | 2026-09-09 | The starting 9-protein table is not an input to anything. | Provenance unknown. |
| D13 | 2026-09-10 | The Gene Ontology lookup name for the Tier 1 seed word "myofibrillar" (D3) is the noun **`myofibril`**, matched by exact term-name equality. | The seed word is an adjective; the ontology names compartments as nouns. The lookup is exact-match only — no related-term search — so that the category is derived from the seed, not chosen. |
| D14 | 2026-09-10 | The Gene Ontology lookup name for the Tier 2 seed word "sarcoplasmic" (D4) is the noun **`sarcoplasm`**, matched the same way. | As D13. |
| D15 | 2026-09-10 | The ontology file used is the **basic** release (`go-basic.obo`, permanent URL `http://purl.obolibrary.org/obo/go/go-basic.obo`). | The consortium's download page (https://geneontology.org/docs/download-ontology/, retrieved 2026-09-10) describes it as filtered to be acyclic with relations restricted to `is a`, `part of`, and the regulates family, and recommends it for annotation tools. The pipeline follows only `is_a` and `part_of`. The file's own `data-version` line, URL, and SHA-256 are recorded at every run. |
| D16 | 2026-09-10 | **Every entry uses its canonical UniProt sequence.** No isoform is selected by any rule. The isoform-deltas tool computes each isoform's amino acid fraction distance from canonical and the methods report the set-wide bound, with no decision attached. | The only structured answer UniProt gives to "which sequence" is the canonical designation. Free-text isoform notes are not a category. If Step 3 mass-weighting shows a bound matters, that is a Step 3 question answered from the abundance data. |
| D17 | 2026-09-10 | **No fiber-type field** in the protein set. | Fiber type is a Layer B quantity (per-fiber-type abundance, including ~0), obtained from the proteomics dataset directly. A field parsed from UniProt free text would be a hand-written string map producing redundant information. |
| D18 | 2026-09-10 | Gene Ontology **evidence codes are recorded, not filtered.** | The set definition is "annotated to the term or its descendants"; the evidence class is disclosed per accession and summarised per tier so a reader can see it. |
| D19 | 2026-09-10 | **Processing rule** reads the UniProt `Chain` feature: exactly one → master range; none → whole sequence, logged; several or non-exact → flag. | The Chain feature is the structured statement of the mature chain. "Last Chain" or any other default would be a silent choice. |
| D20 | 2026-09-10 | The pool is defined as the **union of one UniProt query per subtree term**; UniProt's own term-level query is run as a check and any difference is flagged. | The definition must not depend on how UniProt's search index expands a term. Running both and comparing settles the expansion behaviour empirically and leaves a record. |
| D21 | 2026-09-10 | A surprising member of a subtree or pool is a **flag closed in Step 3 against measurement**, never a reason to re-pick the term. | The category must remain derived. |
| D22 | 2026-09-10 | Step 1 involves **no literature research**. Citations for D3's reasoning move to the paper step; the pool-completeness check against Layer B abundance moves to Step 3. | D3 is admissible as an author decision. Abundance data lives in Step 3. |
| D23 | 2026-09-10 | Nothing from the first repository is inherited. All code is rewritten from the specification in the Step 1 handoff. | `PROVENANCE.md`: a previous version of the repository is not a source. |
