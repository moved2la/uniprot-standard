# Project Plan — Sequence-Derived Amino Acid Standard for the Reference Adult Human

**Status:** Steps 1, 2, 3a, 3b, 4, and 4b complete (2026-09-17). Step 6 (Match Rate against USDA) is next and is the near-term deliverable pre-EOM for the Army DEVCOM RFP response. Step 5 (lab validation) deferred.
**Revision:** 9 — see §10.
**Planning thread:** this document is the output of the planning thread. Each step below is executed in its own thread using its handoff doc in `docs/handoffs/`.
**Last updated:** 2026-09-17 (rev. 9)

---

## 1. Purpose

Replace the laboratory-sampled amino acid profile currently used as the NutriMatch™ "standard" with a **public, auditable, sequence-derived standard** — the amino acid composition of the reference adult human, tissue by tissue, protein pool by protein pool, computed from UniProt protein sequences and cited quantitative proteomics.

The current lab-sampled standard has three problems:

1. It is a black box — quasi-public, not reproducible by a third party.
2. It carries laboratory error: acid hydrolysis destroys tryptophan, deamidates Asn→Asp and Gln→Glu, partially degrades Cys/Met/Ser/Thr; sampling variance between donors, muscles, and biopsy-contamination levels is unquantified.
3. It conflates protein-bound amino acids with free amino acid pool contributions and — in bulk-tissue measurements — with residual blood and connective tissue.

The sequence-derived standard is exact at the protein level, cites its metabolite pools separately, and carries uncertainty only in the tissue weighting, which is explicit, sourced, and propagated.

### Scope (rev. 7, D67 + D68)

The reference is the amino acid composition of the reference adult human, at whole-body scale. It is built as a **mass-weighted composite of per-tissue standards** — each tissue's proteome quantified from its own primary proteomics dataset using the same pipeline machinery. The current published tissues are:

- **Skeletal muscle** (Step 4, done): contractile proteins (Tier 1, myofibril GO set) + support proteins (Tier 2, measured remainder), per fiber type, weighted by fiber-type mix. Cited to Murgia 2021 and Mpampoulis 2021.
- **Bound metabolite pools inside skeletal muscle** (Step 4b, done): carnosine as histidine, at its cited pool size per fiber type, folded into the muscle standard as an inventory (D73); anserine is not carried because human muscle holds none (D75).

Additional tissue categories (blood, liver, collagen tissues, and further tissues as bound-shift math justifies) are added per the [Per-Category Gameplan](per_category_gameplan_R2.md) as separate, self-contained builds. Each category is independently publishable as its own standard and can support its own category-specific Match Rate (e.g., "muscle Match Rate," "collagen Match Rate," "whole-body Match Rate"). The whole-body composite is the mass-weighted summation of the per-category standards.

**Fork framing (D67):** the skeletal muscle standard is a first-class deliverable in its own right, not merely an intermediate. Match Rate can be computed against any single-category standard, against the whole-body composite, or against a purpose-specific subset (e.g., "muscle + blood" for anabolic-response framing). Consumers pick the reference that matches their question.

**What is included and excluded (D68, as executed in D73):** bound non-protein amino acid pools (carnosine now; for non-EAA extensions creatine and glutathione) are included at cited pool sizes — an inventory — within their host tissue's standard; turnover rates are a sensitivity, not part of the standard. These are amino acids the body has committed to non-protein form but must eat to synthesize; excluding them undercounts real dietary demand. The **free amino acid pool** remains excluded by design — it is transient, dynamically regulated, and does not represent bound demand.

In plain terms: the reference is what the body **is made of** — proteins and bound pools together, tissue by tissue. Match Rate Plus (Step 6, taxes and tolls between mouth and cell) handles what makes it there from what is eaten. The two calculations are distinct and remain separate.

End goal: a published methods paper plus an open repository. The standard against which Match Rate is measured should not be a black box.

---

## 2. Design Principles

- **Primary sources only.** Sequences from UniProt REST with MD5 verification. Mass fractions from peer-reviewed supplementary data. Metabolite pool sizes and turnover rates from cited human measurements. Tissue mass fractions from a cited reference-body publication (ICRP 89 or equivalent). No AI-estimated biological constants anywhere in the pipeline.
- **Every literature number carries a citation.** Source, table/figure, and page or supplementary file, on the same row as the number.
- **Human-editable, auditable config.** `.ini` files that a reviewer can open and check line by line. Biology lives in config, not in code.
- **Exactness is layered, and the layers are named.**
  - **Layer A** — per-protein residue composition from sequence. Exact. Zero uncertainty.
  - **Layer B** — mass fraction of each protein in the tissue, by fiber type (for muscle) or by tissue subtype (for other categories). Literature-derived. Uncertainty is quantified and propagated.
  - **Layer C** *(rev. 7, added for whole-body composite)* — mass fraction of each tissue in the reference adult body. Literature-derived from a cited reference-body publication. Uncertainty propagated.
- **Categories are self-contained.** Each tissue/category is its own build following the [Per-Category Gameplan](per_category_gameplan_R2.md), producing its own publishable standard, uncertainty, and stress outputs. The whole-body composite consumes those outputs; it does not recompute them.
- **Rule-of-three abstraction (D71).** The tier-generic pipeline refactor lands *after* the first non-muscle category (Step 7a: blood) is built with muscle-code adaptations. Building blood ad-hoc reveals what the actual axis of variance is between categories; the refactor codifies what is learned rather than what is guessed.
- **Scope is tiered, and tiers are justified in writing** rather than by convention.
- **Corrections are welcome.** Anything found wrong gets fixed and logged in the decisions log (§8).
- **Provenance governs everything.** `PROVENANCE.md` in the repository is the rule this plan operates under: every value is an author decision, a citation, or a computation. An assistant's recollection is never a source.

---

## 3. Architecture

```
UniProt REST ──► fetch ──► sequence store (.ini, MD5-verified)
                                    │
per-category                        │  Layer A (reused across categories)
accessions.ini ─────────────────────┤
segments.ini ───────────────────────┤
                                    ▼
                          composition engine
                    (residue counts → mass vectors,
                     residue-mass AND free-AA-mass conventions)
                                    │
per-category                        │  Layer B (per category, cited per row)
mass_fractions/*.ini ───────────────┤
per-category                        │
metabolite_pools.ini ───────────────┤  (bound pools only, D68)
                                    ▼
                          aggregate  (molar mixture, D64; tier-generic post-D71)
                          uncertainty (log-space Monte Carlo, D63)
                          stress  (named perturbations, D65)
                          plots
                                    │
                                    ▼
              CATEGORY STANDARD (one per tissue/category, per fiber type where applicable)
                                    │
tissue_mass_fractions.ini ──────────┤  Layer C (ICRP 89 or equivalent, D67)
                                    ▼
                          whole-body aggregate
                                    │
                                    ▼
                          WHOLE-BODY STANDARD (composite, mass-weighted)
                                    │
USDA food table ────────────────────┤
match formula ──────────────────────┤
                                    ▼
                    Match Rate (per category OR whole-body, consumer's choice)
```

### Repository layout (target — see `docs/conventions.md` for the current state)

```
uniprot-standard/
├── README.md
├── pyproject.toml
├── config/
│   ├── skeletal_muscle/            # Step 1–4 — muscle category (accessions, segments, mass_fractions, fiber_mixes, metabolite_pools)
│   ├── blood/                      # Step 7a — first non-muscle category, ad-hoc adaptations
│   ├── liver/                      # Step 7c-liver — clean per-category build post-refactor
│   ├── collagen_tissues/           # Step 7c-collagen — clean per-category build post-refactor
│   └── tissue_mass_fractions.ini   # Layer C (D67)
├── data/
│   ├── sequences.ini               # auto-generated, MD5-verified
│   ├── literature/                 # supplementary tables as downloaded, with provenance
│   └── usda/                       # Step 6
├── src/muscle_aa/                  # (name retained for continuity; not muscle-only in scope)
│   ├── fetch.py
│   ├── composition.py
│   ├── aggregate.py                # tier-generic post-D71; consumes any category
│   ├── uncertainty.py
│   ├── stress.py
│   ├── whole_body.py               # NEW — Layer C aggregation
│   ├── validate.py
│   ├── match.py                    # Step 6
│   └── plots.py
├── tests/
├── outputs/
│   ├── skeletal_muscle/            # per-category outputs
│   ├── blood/
│   ├── liver/
│   ├── collagen_tissues/
│   └── whole_body/                 # composite outputs
└── docs/
    ├── 00_project_plan_R<n>.md
    ├── per_category_gameplan_R<n>.md
    ├── planning_summary_expanding_to_full_body_match_R<n>.md
    ├── gorissen_comparison.md      # Step 6 artifact
    ├── adding_a_category.md        # Step 7b artifact — developer-facing map
    ├── handoff_template.md
    ├── handoffs/                   # one per step, one per category execution
    ├── decisions.md
    └── methods.md                  # incremental; becomes the paper's Methods
```

---

## 4. Decisions Already Made (planning thread)

*Authoritative copy: `docs/decisions.md` in the repository. This table is the planning-thread seed (D1–D12); D13 onward live only there.*

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Standard is built first; Match Rate module comes after (Step 6). | The standard is the input to the match calculation; must exist before the calculation is formalized in code. |
| D2 | Two-layer model (Layer A exact / Layer B cited). *(extended by D67 — Layer C for tissue mass fractions in the whole-body composite)* | UniProt gives composition, not abundance. Naming the layers makes the uncertainty honest and publishable. |
| D3 | *(scope superseded by D58, further extended by D67 — the reference is now the whole-body composite, of which the muscle contractile+support standard is one tissue category)* Tier 1 = myofibrillar (contractile) proteins. | Original scoping concept was muscle fiber; MPS literature measures myofibrillar synthesis separately; resistance training primarily drives the myofibrillar fraction. |
| D4 | *(absorbed by D56/D58 — Tier 2 is the measured remainder and part of the primary muscle standard)* Tier 2 = the support proteins. | Originally a sensitivity set for comparison with whole-tissue lab standards; the combined profile is now that comparison. |
| D5 | *(superseded by D20 — the pool is computed from the category; no protein is placed by hand)* | — |
| D6 | Fiber-type-specific entries are used, not proxies. | Satisfied by construction: the computed pool contains each paralog as its own entry; fiber type itself is a Layer B quantity (D17). |
| D7 | Residue counts are the stored ground truth; both mass conventions are derived. | Residue-in-chain mass ≠ free amino acid mass (one H₂O). USDA food tables and lab AAA report free-AA g/100 g. Match Rate must use the free-AA convention to be consistent with the food data. |
| D8 | Keep the `in_master_molecule` / `in_triple_helix`-style segment-flag architecture from the collagen pipeline. | Mature vs precursor is a small effect for cytosolic muscle proteins but the architecture is already proven and the collagen tier reuses it. |
| D9 | Layer B primary source = quantitative human proteomics. Classical biochemical fractionation values are the cross-check. | Proteomics gives human, fiber-type-resolved abundance. Classical fractionation gives direct mass measurement. Agreement between the two is itself a validation result. |
| D10 | Literature is researched in-thread, not deferred as "verify later". | Public scrutiny standard. |
| D11 | Accession list settled before mass fractions. | Can't weight what hasn't been enumerated. |
| D12 | The starting 9-protein table is treated as unverified and is not an input to anything. | Provenance unknown; cannot be reconstructed. |

---

## 5. Steps

Each step = one thread + one handoff doc. A step is done when its acceptance criteria are met and its outputs are committed.

### Step 0 — Repo scaffold and conventions — **DONE** (folded into Step 1)
Directory layout above; `.ini` schema documented in `docs/conventions.md`; `pyproject.toml`; test harness; `decisions.md` seeded from §4.

### Step 1 — Protein set (muscle) — **DONE 2026-09-11**
- **Result:** Tier 1 = `myofibril` GO:0030016, 242 entries; Tier 2 = `sarcoplasm` GO:0016528, 87 entries; 303 unique, all MD5-verified. Decisions D13–D24.
- **Handoff:** `docs/handoffs/01_protein_set_handoff.md`; narrative in `docs/step1_report.md`.

### Step 2 — Composition (Layer A) — **DONE 2026-09-12**
- **Result:** `outputs/composition/amino_acid_composition_per_protein.tsv`, 606 rows; residue counts and both mass conventions. Two disclosures (processing, PTMs) with no decision attached. Decisions D25–D31. 56 tests.
- **Handoff:** `docs/handoffs/02_composition_handoff.md`; narrative in `docs/step2_report.md`.

### Step 3a — Mass fractions: sources and rules — **DONE 2026-09-12**
- **Result:** Primary dataset A1 (Murgia et al. 2021); cross-checks A5, A4, H1. Decisions D31–D49.
- **Handoff:** `docs/handoffs/03a_close_and_03b_handoff.md`.

### Step 3b — Mass fractions: weights — **DONE 2026-09-13**
- **Result:** Weighted set 3,788 entries (242 contractile, 3,546 support); combined standard (D58) contractile share ≈ 72 %. Decisions D50–D60. 80 tests.
- **Handoff:** `docs/handoffs/03b_close_handoff.md`.

### Step 4 — Aggregation, uncertainty, and stress (skeletal muscle standard) — **DONE 2026-09-13**
- **Result:** `python run.py standard` = `aggregate` → `uncertainty` → `stress` → `plots` → tests. First fully-populated `_calculated_amino_acid_standard.tsv` shipped with all ten columns (per-fiber-type contractile, builders, total; standard combined via Mpampoulis 2021 fiber-type mix). Rules A0–A9. Decisions D61–D66. Test tolerance corrected for 4-decimal output precision (from `1e-6` to `2e-4`; principled minimum for the observed drift).
- **Known open item, not blocking:** the sum row of `_calculated_amino_acid_standard.tsv` is computed from unrounded internals and prints `100.0000` even when the visible rows sum to `100.0002`. Either compute the sum from the rounded values, or largest-remainder-adjust the 20 values so they genuinely sum to 100.0000. Deferred.
- **Handoff:** `docs/handoffs/04_aggregation_handoff.md`; narrative in `docs/step4_report.md`.

### Step 4b — Muscle metabolite adjustment (bound pools) — **DONE 2026-09-17**
- **Result:** `python run.py standard` gains stage `bound_pools` (rule A10) between `aggregate` and `uncertainty`; `outputs/standard/_calculated_amino_acid_standard_with_bound_pools.tsv` beside the unchanged protein-only table, version label `ref = skeletal muscle protein + carnosine`. Histidine 2.227 → 2.568 % of amino acid mass in the standard (+15 %); fiber type I +9.5 %, IIa and IIx +21 %; every other amino acid −0.35 % of its value. Inputs: Harris, Dunnett & Greenhaff 1998 (single-fiber carnosine, abstract), Mingrone et al. 2001 (protein and water per kg wet muscle), Harris et al. 2012 (carnosine the sole dipeptide). Decisions D72–D81. 121 tests.
- **Not computed, stated:** the resupply-frame sensitivity (no muscle-protein turnover source), the whole-muscle cross-check (Vega 2022 dropped, gated). The single-fiber spread sensitivity is written.
- **Rules changed in this step:** transcription may be assistant-extracted with page-level locations, author-audited (D79); fetch-by-code retired, pinned for 7b (D80); retrieval times out of generated config (D81); fetcher patches F2b, F3, F4b, F7b.
- **Handoff:** `docs/handoffs/04b_metabolite_adjustment_handoff.md`; narrative in `docs/step4b_report.md`.

### Step 5 — Validation against laboratory data — **DEFERRED**
Moved to future work / grant-funded research. The core reason to keep Step 5 in R6 (validating the muscle standard against lab-sampled human muscle AAA) is partly satisfied by the Gorissen reconciliation carried into Step 6. A full compiled dataset of published human muscle AAA measurements with per-amino-acid deviation attribution is out of scope for the pre-paper arc.

### Step 6 — Match Rate against USDA — **NEXT; NEAR-TERM DELIVERABLE (pre-EOM for Army DEVCOM RFP)**
- **Goal:** Formalize the Match Rate calculation in code against the new standard; run against the USDA food database; publish version-labeled results.
- **Deliverables:** Formula written out in `docs/methods.md`; `pipeline/match.py`; USDA ingestion in the free-AA convention; Match Rate for every food in the database against each available reference (muscle standard as of Step 4b, and any additional category standards or whole-body composite that are available at run time); re-scored reference foods (whey, pea, egg, steak, current blend) under old vs new standards, version-labeled; **`docs/gorissen_comparison.md`** as a supplementary artifact reconciling the new standard against the Gorissen 2018 reference the original NutriMatch was calibrated against.
- **Fork point:** Match Rate v1.0 runs against whichever reference is ready first — muscle-with-metabolites at minimum, whole-body composite if Step 7 iterations have completed enough categories. Each Match Rate output carries its reference version label (e.g., "Match Rate v1.0, ref: skeletal muscle + carnosine, D67/D68"). Version labeling prevents cross-version comparisons from looking like methodology inconsistency.
- **Input needed from Anthony:** the current spreadsheet formula; confirmation of EAA-only scope vs all-AA (working assumption: EAA-only per this session).
- **Handoff:** `docs/handoffs/06_match_rate_handoff.md` (drafted 2026-09-17 at the Step 4b close).

### Step 7a — First category expansion: blood — **OPEN**
- **Goal:** build the first non-muscle category standard, using the current (muscle-adapted) pipeline code. This is the deliberate second occurrence of the pattern (rule-of-three): muscle invented it, blood reveals its actual axis of variance, then Step 7b codifies what is learned.
- **Approach:** follow the [Per-Category Gameplan](per_category_gameplan_R2.md), but with the understanding that reused-code steps (Steps 4, 5, 7, 8 of the gameplan) will require ad-hoc muscle-code adaptations rather than clean parameterization. Record the adaptations in the handoff — they are the input to Step 7b's refactor scope.
- **Deliverables:** blood category standard (`outputs/blood/_calculated_amino_acid_standard.tsv`), uncertainty and stress outputs, tissue mass fraction added to `config/tissue_mass_fractions.ini`. Handoff explicitly names what code paths required adaptation and what the muscle-specific assumptions were.
- **Rationale for going first (before Step 7b refactor):** categories 2+ are the ones that benefit from a genuinely modular pipeline. Refactoring before category 2 is guessing what modularity should look like; refactoring after category 2 is codifying what has been observed. The rule-of-three principle (D71).
- **Handoff:** `docs/handoffs/07a_blood_handoff.md` (to be written).

### Step 7b — Pipeline refactor for per-category modularity — **AFTER STEP 7a**
- **Goal (D71):** generalize the pipeline to be genuinely tier- and category-agnostic, informed by what Step 7a revealed about muscle-specific assumptions.
- **Deliverables:**
  - Aggregate stage reads tier names, subtype names, and column templates from a per-category manifest instead of hard-coding "contractile / builders / type_I/IIa/IIx"
  - Config directory layout matches R7 §3 (per-category subdirectories under `config/`; per-category output subdirectories under `outputs/`)
  - Subtype-mix loader generalized (fiber-type mix becomes one instance of a per-category subtype mix)
  - Output header rendering templated with category-level substitutions
  - Tests parametrized over category × tier × subtype rather than hard-coded
  - `run.py standard <category>` and `run.py standard --all` commands
  - `config/tissue_mass_fractions.ini` scaffolding laid in (actual aggregation code deferred to Step 8)
  - `docs/adding_a_category.md` — developer-facing map showing where a new category is added (complements the per-category gameplan)
  - **From Step 4b (D80, D81):** `fetch-literature` becomes its own command whose stage only hashes what is in `data/literature/` and writes the manifest — no download attempts, no re-fetching of static files; the source list leaves `mass_fraction_decisions.ini` for its own file; each stage hashes only the config sections it reads, so bookkeeping edits cannot make unrelated outputs read as stale
  - **From Step 4b:** `docs/run_order.md` — a one-screen cheat sheet of what depends on what and which command to run after which edit, kept separate from the README; the README itself trimmed
- **Acceptance criteria:**
  - Skeletal muscle standard rebuilds byte-identical to pre-refactor output (primary regression test)
  - Blood standard rebuilds through the refactored pipeline and matches the Step 7a output
  - All existing tests pass after parametrization updates
- **Size:** 1–2 days, informed by what Step 7a required.
- **Handoff:** `docs/handoffs/07b_refactor_handoff.md` (to be written).

### Step 7c — Recurring category expansions per gameplan — **LOOPS**
- **Goal:** add subsequent tissue/category standards to the whole-body composite following the [Per-Category Gameplan](per_category_gameplan_R2.md), using the fully modular post-refactor code.
- **Categories in EAA-priority order** (from `planning_summary_expanding_to_full_body_match_R3.md`):
  1. **Liver** (~5 % body EAA) — Kim 2014 or Wang 2019 or Human Protein Atlas as candidate primary sources.
  2. **Collagen tissues** (~10 % body EAA combined) — already has its own project track; folds in as one lumped tier.
  3. Additional tissues (GI, brain, adipose protein, heart, kidneys, lungs, hair/nails, minor) as bound-shift math justifies. Stopping criterion: **provable bound** on how much uncharacterized remainder can shift any EAA in the composite, not "diminishing returns."
- **Per-category handoffs:** `docs/handoffs/07c_<category>_handoff.md`.
- **Deliverables per category:** category standard TSV, category uncertainty, category stress outputs, category-specific Match Rate potential. Runs cleanly through modular pipeline; no code changes expected per category (only config additions).

### Step 8 — Whole-body composite — **OPEN, RUN WHEN CATEGORIES JUSTIFY**
- **Goal:** aggregate category standards into the whole-body composite using tissue mass fractions (Layer C). Compute uncertainty envelope and bound-shift sensitivity for the uncharacterized remainder.
- **Deliverables:** `_amino_acid_standard_whole_body.tsv`; whole-body Match Rate; sensitivity analysis showing worst-case shift under alternate lumped-remainder assumptions.
- **Runnable at any point** after at least one Step 7c category is complete (or possibly after Step 7a if the aggregation is written to accept blood alone). Each rerun improves coverage; the version label records which categories were included.
- **Handoff:** `docs/handoffs/08_whole_body_composite_handoff.md` (to be written).

### Step 9 — Paper — **WRITEABLE AT ANY POINT AFTER STEP 6**
- **Goal:** Methods paper.
- **Deliverables:** `docs/methods.md` assembled from per-step docs; figures from Steps 4/6/7/8; Gorissen reconciliation as a named subsection; supplementary data = the config files themselves.
- **Cut-off decision:** the paper can be written at any point after Step 6 with whatever category coverage is complete. Coverage is a scope statement, not a gate. Version-labeled numbers throughout.

---

## 6. Conventions (formalized in `docs/conventions.md`)

**Config format:** `.ini` via `configparser`. One section per entity. Free-text `note` and `justification` keys allowed. Comments with `#`.

**Citation schema** (every literature-derived number):
```
value     = 0.431
unit      = mass_fraction_of_tier
source    = Author et al. YEAR, Journal Vol:pages, DOI
location  = Supplementary Table S3, column "iBAQ_typeI"
retrieved = 2026-09-XX
note      = optional — conversions applied, caveats
```

**Mass conventions:**
- `residue_mass` — in-chain residue mass (free AA minus 18.015 Da). Sums to protein MW.
- `free_aa_mass` — free amino acid mass. Sums to > protein MW. This is the USDA / lab AAA convention and the one Match Rate uses.
- Both reported. Ground truth is residue count.

**Accession identity:** `ACCESSION` or `ACCESSION-N` for a specific isoform. Sequence version recorded. MD5 recorded and checked against UniProt's published checksum.

**Version labeling:** every published standard output carries a version label naming the categories and metabolite pools included, and referencing the D-numbers that scope it. Every Match Rate output carries the reference version it was computed against.

**Author verification for citations:** author lists must be verified by fetching the article, not from search snippets. (Standing rule reinforced this session.)

---

## 7. Open Questions (carried forward, owned by the step that resolves them)

| Question | Resolved in |
|----------|-------------|
| Carnosine pool size in human vastus lateralis | Resolved Step 4b: Harris, Dunnett & Greenhaff 1998 (single fibers, abstract); anserine not carried (Harris et al. 2012) |
| Muscle-protein turnover range for the resupply-frame sensitivity | Open, not blocking: the sensitivity is stated as not computed until the primary studies at each end of the published range are on disk |
| Whether to run Match Rate v1.0 against muscle-with-metabolites first, or wait for whole-body composite | Step 6 (fork decision; both muscle references now exist with version labels) |
| Reference-body citation for tissue mass fractions (ICRP 89 or Snyder ICRP 23 or equivalent) | Step 8 (author verification required) |
| Primary quantitative proteomics source for blood tier | Step 7a |
| Which muscle-specific assumptions actually create friction for other categories (input to Step 7b refactor scope) | Step 7a handoff — recorded as it comes up |
| Primary quantitative proteomics source for liver tier (Kim 2014 vs Wang 2019 vs Human Protein Atlas) | Step 7c-liver |
| Collagen tissue mass fractions (skin, tendon, ligament, bone matrix, fascia — lumped vs subdivided) | Step 7c-collagen (working assumption: lumped, per this session) |
| Match Rate scope: EAA-only vs all AA; normalization; basis | Step 6 (working assumption: EAA-only per this session) |
| Whether the muscle standard's sum-row display inconsistency gets fixed (largest-remainder adjustment vs computed-from-rounded) | Deferred; not blocking any downstream step |
| Fiber-type mix — accept Mpampoulis 2021 (female-only, %CSA) with disclosed limitation, or search for a mixed-sex source | Open; not blocking |
| Stopping criterion for adding tissues to whole-body composite | Step 8 sensitivity analysis; provable bound on remainder-shift |
| Isoform, processing, and PTM bounds — material after weighting? | Resolved rev. 5; carried to Step 4 as bounds |
| Glycosylation mass (unbounded from the entry) | Resolved rev. 5; stated as a limitation |
| Which lab AAA datasets of human muscle to validate against | Was Step 5; deferred with Step 5 |
| Match Rate delivery layer: proportion × absorption × burn (Match Rate Plus) | Step 6 or later; distinct from the reference itself |
| Momenzadeh 2023: which table (S2 / S3) matches the paper's text on MYH4 | Paper step (Step 9) |
| Structural titin : MHC anchor from sarcomere stoichiometry | Paper step (Step 9) |
| Per-fiber values (PRIDE PXD006182 deposit; Murgia 2017 tables) for a bootstrap of the two uncertainty quantities | Not blocking; upgrade opportunity |

---

## 8. Decisions Log

Maintained in `docs/decisions.md`. Seeded from §4. Every subsequent decision made in a step thread is appended there with date, step, and rationale.

**New decisions in R7 (to be logged in `docs/decisions.md` when the reframe is committed):**
- **D67** — Scope extended from skeletal muscle to whole-body composite. Muscle standard remains a first-class deliverable; the whole-body composite is the mass-weighted summation of per-tissue category standards. Each category is independently publishable. Layer C added to the layered-scope model (tissue mass fractions from a cited reference-body publication).
- **D68** — Bound metabolite pools (carnosine, anserine, and — for non-EAA extensions — creatine, glutathione) are included at cited pool sizes and turnover rates within their host tissue's standard. Free amino acid pool remains excluded by design. Reverses R6 §1 problem #3's implicit rejection of carnosine/anserine.
- **D69** — Project title updated from "for Human Skeletal Muscle" to "for the Reference Adult Human." The muscle standard is retitled as an intermediate deliverable ("Skeletal Muscle Protein Composition Standard").
- **D70** — Per-category gameplan template established. Twelve-step procedure in `docs/per_category_gameplan_R<n>.md`. New categories follow the template as their own handoff, without requiring a project-plan revision per category.

**New decisions in R8:**
- **D71** — Pipeline refactor for per-category modularity is sequenced AFTER the first non-muscle category (Step 7a: blood), not before. Rule-of-three abstraction: muscle invented the pattern (Step 4), blood built with ad-hoc adaptations reveals the actual axis of variance (Step 7a), refactor codifies what is observed (Step 7b). Building blood ad-hoc is intentional — it produces the specification for what the refactor needs to generalize. Subsequent categories (Step 7c+) drop into the modular pipeline cleanly.

**Logged in R9 (all in `docs/decisions.md`):** D67–D71 written into the log (they had been plan-only); **D72** the fiber-type mix (Mpampoulis et al. 2021, %CSA, female); **D73** inventory frame, resupply a sensitivity; **D74** per-fiber-type single-fiber values, no whole-muscle cross-check; **D75** carnosine only; **D76** protein-only table untouched, adjusted table separate; **D77** Mingrone et al. 2001 for protein and water per kg muscle, measured over computed; **D78** green light before calculation; **D79** assistant-extracted transcription with page-level locations, author-audited; **D80** fetch-by-code retired (7b); **D81** retrieval times out of generated config.

---

## 9. Thread Protocol

1. Open a new thread named for the step (e.g., "Step 4b — Muscle metabolite adjustment," "Step 7a — Blood category").
2. Paste the handoff doc as the first message.
3. Work only that step. Anything out of scope goes into the handoff's "Deferred / raised" section, not into the work.
4. On completion: update `decisions.md`, write the next step's handoff (or update the draft), tick this plan's step status and add a row to §10, commit.
5. Return to the planning thread only to re-plan, not to do work.
6. For Step 7a and Step 7c category expansions: follow the [Per-Category Gameplan](per_category_gameplan_R2.md) as the executing template. The gameplan itself is not re-revised for each category (though it may be revised after Step 7a+7b close, informed by what those steps taught).

---

## 10. Revision Log

| Rev | Date | Change |
|---|---|---|
| 1 | 2026-09-09 | Planning thread output. |
| 2 | 2026-09-11 | Step 0 and Step 1 marked done; Step 1 goal restated as executed (D22); Step 2 rewritten; §7 rows resolved or moved; `PROVENANCE.md` added; §3 layout deferred to `docs/conventions.md`. |
| 3 | 2026-09-12 | Step 2 marked done (D25–D30); Step 3 marked next; §7 updates. |
| 4 | 2026-09-12 | Step 3 split into 3a (done: D31–D49) and 3b (next). |
| 5 | 2026-09-13 | Step 3b marked done (D50–D60); §1 scope restated under D58; Step 4 marked next. |
| 6 | 2026-09-13 | Step 4 marked in progress and restated as executed: `standard` command with four stages incl. stress (D65); uncertainty log-space (D63); molar profile (D64); FAO (D62); §3 architecture node updated. |
| 7 | 2026-09-13 | Step 4 marked done. **Scope reframed (D67):** whole-body composite as the reference, muscle as one category. **Metabolite pools included (D68):** bound pools in, free pool out. **Title change (D69):** "for the Reference Adult Human." **Per-category gameplan (D70)** established as its own artifact. Step 4b (muscle metabolite adjustment) inserted. Step 5 (lab validation) deferred. Step 6 (Match Rate) elevated as near-term deliverable for Army DEVCOM RFP; Gorissen reconciliation folded in as `docs/gorissen_comparison.md`. New Steps 7 (category expansions, looping per gameplan) and 8 (whole-body composite). Step 9 (paper) writeable at any point after Step 6. §3 architecture updated for per-category modules and Layer C. §7 rows updated. |
| 9 | 2026-09-17 | **Step 4b marked done** (D72–D81): `bound_pools` stage, `_calculated_amino_acid_standard_with_bound_pools.tsv`, histidine +15 % in the standard. §1 scope wording updated (inventory, not turnover; anserine not carried). Step 6 marked next with its handoff drafted. §5 Step 7b gains the fetch retirement, the source-list split, per-stage hashing, and the run-order cheat sheet. §7 carnosine row resolved; turnover-range row added. |
| 8 | 2026-09-13 | **Pipeline refactor sequenced after first non-muscle category (D71):** Step 7 restructured into 7a (blood, first category using muscle-adapted code), 7b (pipeline refactor codifying what 7a taught), and 7c (subsequent categories using modular code). §2 design principles gain "rule-of-three abstraction" bullet. §3 architecture note that aggregate is tier-generic post-D71. §5 Step 7 subdivided with acceptance criteria for refactor (byte-identical muscle rebuild, blood rebuild matches Step 7a output). §7 open questions gain row for "which muscle-specific assumptions actually create friction" (input to 7b scope, recorded during 7a). §11 references gain `adding_a_category.md` as 7b artifact. Per-category gameplan referenced as R2 (updated in parallel to reflect the after-blood refactor timing). |

---

## 11. Referenced Artifacts

- [`per_category_gameplan_R2.md`](per_category_gameplan_R2.md) — the twelve-step template for adding any new tissue/category to the whole-body composite. Referenced by Steps 7a, 7c, and 9. May be revised (R3+) after Step 7a+7b close, informed by what those steps teach about the actual per-category friction points.
- [`planning_summary_expanding_to_full_body_match_R3.md`](planning_summary_expanding_to_full_body_match_R3.md) — the reframe planning output that motivated R7. Contains the Gorissen comparison table, carnosine math walk-through, collagen EAA composition, cumulative body composition sketch, and per-tissue EAA ranking. Referenced by Steps 4b, 6, and 9.
- `docs/gorissen_comparison.md` (to be written in Step 6) — supplementary artifact reconciling the new standard against Gorissen 2018.
- `docs/adding_a_category.md` (to be written in Step 7b) — developer-facing map showing where in the code a new category is added; complements the per-category gameplan (which is planning-facing) with the code-facing perspective.
