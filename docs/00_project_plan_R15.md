# Project Plan — Sequence-Derived Amino Acid Standard for the Reference Adult Human

**Status:** Steps 1, 2, 3a, 3b, 4, 4b, 6 and 7a complete (2026-09-18). Match Rate 1.0 runs against the skeletal muscle standard with its non-protein metabolite pools and scores 4,952 USDA foods. **Step 7b (blood) is in progress:** part 1 — scope, sources, decisions D95–D105, the F7 amendment — closed 2026-09-18; part 2 — the build through the run order — continues in its own thread. Then Step 7c (the abstraction), then Step 7d (the category loops). Step 5 (lab validation) deferred.
**Revision:** 15 — see 10.
**Planning thread:** this document is the output of the planning thread. Each step below is executed in its own thread using its handoff doc in `docs/handoffs/`.
**Last updated:** 2026-09-18 (rev. 15)

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
- **Rule-of-three abstraction (D71; letters renumbered in R12).** The tier-generic pipeline refactor (Step 7c) lands *after* the first non-muscle category (Step 7b: blood) is built with muscle-code adaptations. Building blood ad-hoc reveals what the actual axis of variance is between categories; the refactor codifies what is learned rather than what is guessed. The repository reorganization and the bookkeeping items pinned in Steps 4b and 6 need nothing from blood and come first (Step 7a, D89).
- **Scope is tiered, and tiers are justified in writing** rather than by convention.
- **Corrections are welcome.** Anything found wrong gets fixed and logged in the decisions log (8).
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

### Repository layout (as built, Step 7a)

```
uniprot-standard/
├── README.md
├── requirements.txt
├── run.py                          # the one entry point, eight commands
├── config/                         # muscle's files stay flat
│   ├── literature_sources.ini      # every source; used_for; [categories] declares the map's columns
│   ├── tissue_mass_fractions.ini   # Layer C (D67) — Step 8
│   ├── blood/                      # Step 7b — every new category gets its own subfolder
│   ├── liver/                      # Step 7d
│   └── collagen_tissues/           # Step 7d
├── data/
│   ├── uniprot_sequences.ini       # shared Layer A store, MD5-verified
│   ├── literature/<source>/        # one folder per source, not committed
│   │   ├── manifest.ini            # generated: url, path, sha256, bytes
│   │   ├── manifest_map.tsv        # generated: source × purpose, A-Z
│   │   └── README.md               # generated provenance
│   └── usda/
├── pipeline/                       # code, one module per stage
├── tests/
├── outputs/
│   ├── flags.tsv
│   ├── intermediate/               # feeds the standard; is not the standard
│   │   ├── protein_set/  composition/  digest/
│   │   └── literature_inventory/  mass_fractions/
│   ├── standard/                   # the deliverable tables and profiles at the top
│   │   ├── plots/                  # the plots of the standard itself
│   │   ├── uncertainty/ + plots/
│   │   ├── stress/
│   │   └── sensitivity/ + plots/
│   ├── usda/  comparison/
│   ├── match/                      # per-food tables at the top
│   │   └── <reference>/            # one folder per [reference.*] in config/match_rate.ini
│   ├── blood/                      # Step 7b — intermediate/ and standard/ inside
│   └── whole_body/                 # Step 8
└── docs/
    ├── 00_project_plan_R<n>.md
    ├── run_order.md                # what to rerun after an edit
    ├── pipeline_map.md
    ├── conventions.md  methods.md  decisions.md
    ├── formula.md  gorissen_comparison.md
    ├── step<N>_report.md
    ├── adding_a_category.md        # Step 7c
    └── handoffs/                   # one per step (not committed)
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
Directory layout above; `.ini` schema documented in `docs/conventions.md`; `pyproject.toml`; test harness; `decisions.md` seeded from 4.

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
- **Result:** `python run.py standard` gains stage `non_protein_metabolite_pools` (rule A10) between `aggregate` and `uncertainty`; `outputs/standard/_calculated_amino_acid_standard_with_bound_pools.tsv` beside the unchanged protein-only table, version label `ref = skeletal muscle protein + carnosine`. Histidine 2.227 → 2.568 % of amino acid mass in the standard (+15 %); fiber type I +9.5 %, IIa and IIx +21 %; every other amino acid −0.35 % of its value. Inputs: Harris, Dunnett & Greenhaff 1998 (single-fiber carnosine, abstract), Mingrone et al. 2001 (protein and water per kg wet muscle), Harris et al. 2012 (carnosine the sole dipeptide). Decisions D72–D81. 121 tests.
- **Not computed, stated:** the resupply-frame sensitivity (no muscle-protein turnover source), the whole-muscle cross-check (Vega 2022 dropped, gated). The single-fiber spread sensitivity is written.
- **Rules changed in this step:** transcription may be assistant-extracted with page-level locations, author-audited (D79); fetch-by-code retired, pinned for 7a (D80, moved there by D89); retrieval times out of generated config (D81); fetcher patches F2b, F3, F4b, F7b.
- **Handoff:** `docs/handoffs/04b_metabolite_adjustment_handoff.md`; narrative in `docs/step4b_report.md`.

### Step 5 — Validation against laboratory data — **DEFERRED**
Moved to future work / grant-funded research. The core reason to keep Step 5 in R6 (validating the muscle standard against lab-sampled human muscle AAA) is partly satisfied by the Gorissen reconciliation, **delivered in Step 6 as `docs/gorissen_comparison.md`** (D84, D85). A full compiled dataset of published human muscle AAA measurements with per-amino-acid deviation attribution is out of scope for the pre-paper arc.

### Step 6 — Match Rate against USDA — **DONE 2026-09-18**
- **Goal:** Formalize the Match Rate calculation in code against the new standard; run against the USDA food database; publish version-labeled results.
- **Result:** `python run.py usda` → `python run.py comparison` → `python run.py match`. **Match Rate 1.0** scores every USDA food with an amino acid panel (5,156 foods, D83) against three references — the skeletal muscle standard with its non-protein metabolite pools (the primary, D82), the protein-only standard, and Gorissen et al. 2018's measured human muscle column (the original NutriMatch reference, on the eight amino acids it measures). Formula transcribed from the author's spreadsheet (D87): both sides as shares of the scored set, the score the smallest of (food share ÷ reference share), its amino acid the limiting one; the spreadsheet's surplus arithmetic reduces to the same number exactly; protein content never enters. Scored set = FAO 2013's nine headings, each group reduced to its indispensable member (D86). Two scripts: the public single-food scorer here, the blend/fortification script separate and private at publication (D88). Rules U1–U7 (USDA), X1–X6 (comparison), M1–M6 (Match Rate). Decisions D82–D88. Six output tables per run, including `match_rate_steps.tsv` — the spreadsheet's own walk per food and reference, so any published score can be rebuilt in Excel from its row.
- **What the run shows** (2026-09-18): median 69.9 % (primary) against 74.7 % under the original reference; **tryptophan limits 49.7 % of scored foods under the new reference** and is the one scored amino acid acid hydrolysis destroys, so the measured comparator cannot check it — the number deciding the most scores has the least independent measurement behind it and is exact at Layer A; phenylalanine and histidine, which limit 62 % of foods under the original reference, limit under 1 % under the new one (the Gorissen comparison seen from the food side). Documented, not acted on.
- **Gorissen comparison** (D84, D85): fourteen rows across fifteen amino acids, both sides renormalised, no verdict; valine 1.00×, leucine 0.97×, alanine 0.96×, lysine 1.04×; carnosine closes 21 % of the protein-only histidine gap; the mass balance (X6) is consistent with a protein content between the paper's nitrogen-derived 84 % and Mingrone 2001's direct measurement, reported as a bound, not a test.
- **Pinned by the author (do nothing yet, all in the handoff):** a tryptophan-preserving cross-check source; a minimum protein content for the ranked tables (lean: a floor in `config/match_rate.ini` applied to the ranked view only); consolidation of USDA preparations whose panels differ within rounding; archive precedence; cystine → cysteine; the reference-foods table; the leucine correction factor (with the blend script); Match Rate Plus (proportion × absorption × burn) as a later layer.
- **Handoff:** `docs/handoffs/06_match_rate_handoff.md` (closed); narrative in `docs/step6_report.md`.

### Step 7a — Reorganization — **DONE 2026-09-18**
- **Result:** `outputs/` laid out by what a file is, not which stage wrote it; the literature sources in their own config file with a category map; `fetch-literature` as its own command that never downloads; every generated path resolved at call time; a stage hashes only what it reads. **No number changed** — the acceptance check over 111 files found every difference explained by a moved path, a run timestamp, or a hash cascading one hop from a header that names a moved path. Decisions D90–D94. 213 tests.
- **`outputs/` as built**
  - `intermediate/` — `protein_set/` 11, `composition/` 7, `digest/` 6, `literature_inventory/` 3, `mass_fractions/` 17
  - `standard/` — 13 deliverable and profile tables at the top; `plots/` 12; `uncertainty/` 5 + `plots/` 6; `stress/` 5; `sensitivity/` 4 + `plots/` 2
  - `match/` — 6 at the top, one folder per `[reference.*]` of `config/match_rate.ini`
  - `usda/` 4, `comparison/` 4, `flags.tsv` — unchanged
  - a new category gets `outputs/<category>/{intermediate,standard}/` via `common.category_outputs()`
- **`config/`** — unchanged but for one new file, `literature_sources.ini`: the 22 `[source.*]` sections moved out of `mass_fraction_decisions.ini` unchanged, each with `used_for`, plus `[categories]` declaring the columns of the map. `data/` unchanged but for the generated `manifest_map.tsv`.
- **Bookkeeping closed:** `fetch-literature` its own command (D94); per-stage hashing with no timestamps in a hash (D92); the placeholder rule (`___` is only ever typed by a person; an unknown hash is blank); one excerpt spec per command; `docs/run_order.md`; README trimmed.
- **Rules:** A0 amended (D92); F3, F4, F4b, F6 retired with the download code; F8 (cited but not stored) and F9 (`used_for`) added.
- **Not done, deliberately:** the header rewrite of the older stages. Bringing `protein-set`, `composition`, `mass-fractions`, `aggregate`, `uncertainty` and `stress` to provenance-only headers changes headers on purpose and is better done with the header templating of Step 7c than twice.
- **Open, recorded in the report:** `comparison` has no `--check` currency test; `excerpt.py` lists three sensitivity files the pipeline never writes; whether a zero-flag run should write a header-only `flags.tsv`.
- **Handoff:** `docs/handoffs/07a_layout_handoff.md` (closed); narrative in `docs/step7a_report.md`.

### Step 7b — First category expansion: blood — **IN PROGRESS (part 1 closed 2026-09-18)**
- **Goal:** build the first non-muscle category standard, using the current (muscle-adapted) pipeline code. This is the deliberate second occurrence of the pattern (rule-of-three): muscle invented it, blood reveals its actual axis of variance, then Step 7c codifies what is learned.
- **Part 1 — scope, sources, decisions (closed 2026-09-18, D95–D105).**
  - Blood is plasma + erythrocytes, two compartments each from its own proteome; leucocytes and platelets excluded and bounded (D95). Blood owns plasma proteins; myoglobin stays with muscle (D96). In code terms one tier × two subtypes, mixed by a cited protein-mass share; the muscle names are left in place and recorded for 7c (D97). Pools are what the datasets quantified, with a per-compartment contaminant rule (D98). No metabolite pool this step (D99).
  - **The split** (D100): measured reference values, no remainders — ICRP 89 volumes and blood mass; NORIP concentrations (rustad_2004 plasma total protein, nordin_2004 hemoglobin / haematocrit / MCHC); adult male, midpoint of the male reference interval, female beside; ICRP's whole-blood protein and Hortin's summed table as cross-checks.
  - **Sources** (D101): `bryk_2017` primary erythrocytes (Table S3 per-donor mass shares, four donors); `gautier_2018` cross-check; `geyer_2016` primary plasma (Table S2, one donor, fifteen replicates, spread stated as technical; Table S5 for the completeness bound); `hortin_2008` immunoassay cross-check (153 rows with accessions); `anderson_2002` reference only; `rustad_2004`, `nordin_2004` the split; `icrp_2002` now `used_for = blood`. Ten sources, 34 files hashed locally, `fetch-literature` green, map 29 × 5.
  - **Two findings that shape the build** (D102): label-free MS under-weights the dominant protein of both compartments — hemoglobin 53 % (bryk) / ~70 % (gautier) / ~95 % (classical); albumin 24 % (geyer) / >50 % (hortin). Weights are taken as measured; the dominant protein is a named sensitivity across the measured shares, no adjustment. Immunoglobulins by constant-region entries, variable domains as an A8 bound (D103). Each source declares its abundance kind, `mass_share` or `molar`, and the code applies it (D105).
  - **Process** (D104): F7 amended — the stored file is the downloaded name (`file.<n>.name`), URL provenance only; the source acquisition protocol written into the gameplan (R3 §2a) after the thread lost an afternoon to the assistant writing URLs and hashes. Delivery 1 R2 (`fetch_literature.py`, its tests, the F7 row).
- **Part 2 — the build (next thread; `docs/handoffs/07b_blood_handoff.md` R2):** excerpts of the four tables; protein set for the two compartments; composition; mass fractions with the D105 unit rule; aggregate as one tier × two subtypes; uncertainty; stress; the D102 sensitivity; `config/blood/` and `outputs/blood/`; the Layer C row; **the adaptation list as the acceptance criterion**; `docs/step7b_report.md` completed; the Step 7c handoff.
- **Deliverables:** blood category standard (`outputs/blood/standard/_calculated_amino_acid_standard.tsv`), uncertainty and stress outputs, tissue mass fraction added to `config/tissue_mass_fractions.ini`. Handoff explicitly names what code paths required adaptation and what the muscle-specific assumptions were.
- **Rationale for going before the Step 7c refactor:** categories 2+ are the ones that benefit from a genuinely modular pipeline. Refactoring before category 2 is guessing what modularity should look like; refactoring after category 2 is codifying what has been observed. The rule-of-three principle (D71).
- **Upgrades recorded, not fetched:** the ten-individual plasma set (PRIDE PXD002854 `proteinGroups.txt`); Gautier's per-protein table (PXD009258 `txt.zip`, 4 GB).

### Step 7c — Pipeline refactor for per-category modularity — **AFTER STEP 7b**
- **Goal (D71):** generalize the pipeline to be genuinely tier- and category-agnostic, informed by what Step 7b revealed about muscle-specific assumptions. The bookkeeping items this step carried in R9–R11 moved to Step 7a (D89); what remains is the abstraction.
- **Deliverables:**
  - Aggregate stage reads tier names, subtype names, and column templates from a per-category manifest instead of hard-coding "contractile / builders / type_I/IIa/IIx"
  - Subtype-mix loader generalized (fiber-type mix becomes one instance of a per-category subtype mix)
  - Output header rendering templated with category-level substitutions
  - Tests parametrized over category × tier × subtype rather than hard-coded
  - `run.py standard <category>` and `run.py standard --all` commands
  - `config/tissue_mass_fractions.ini` scaffolding laid in (actual aggregation code deferred to Step 8)
  - Whether muscle's own config and outputs move under `skeletal_muscle/` subfolders — a decision for this step, not a default
  - `docs/adding_a_category.md` — developer-facing map showing where a new category is added (complements the per-category gameplan)
  - *(Moved to Step 7a by D89: the fetch retirement and source-list split, per-stage hashing, the placeholder convention, per-stage excerpt specs, the header rule, `docs/run_order.md`.)*
- **Acceptance criteria:**
  - Skeletal muscle standard rebuilds byte-identical to pre-refactor output (primary regression test)
  - Blood standard rebuilds through the refactored pipeline and matches the Step 7b output
  - All existing tests pass after parametrization updates
- **Size:** 1–2 days, informed by what Step 7b required.
- **Handoff:** `docs/handoffs/07c_refactor_handoff.md` (to be written).

### Step 7d — Recurring category expansions per gameplan — **LOOPS**
- **Goal:** add subsequent tissue/category standards to the whole-body composite following the [Per-Category Gameplan](per_category_gameplan_R2.md), using the fully modular post-refactor code.
- **Categories in EAA-priority order** (from `planning_summary_expanding_to_full_body_match_R3.md`):
  1. **Liver** (~5 % body EAA) — Kim 2014 or Wang 2019 or Human Protein Atlas as candidate primary sources.
  2. **Collagen tissues** (~10 % body EAA combined) — already has its own project track; folds in as one lumped tier.
  3. Additional tissues (GI, brain, adipose protein, heart, kidneys, lungs, hair/nails, minor) as bound-shift math justifies. Stopping criterion: **provable bound** on how much uncharacterized remainder can shift any EAA in the composite, not "diminishing returns."
- **Per-category handoffs:** `docs/handoffs/07d_<category>_handoff.md`.
- **Deliverables per category:** category standard TSV, category uncertainty, category stress outputs, category-specific Match Rate potential. Runs cleanly through modular pipeline; no code changes expected per category (only config additions).

### Step 8 — Whole-body composite — **OPEN, RUN WHEN CATEGORIES JUSTIFY**
- **Goal:** aggregate category standards into the whole-body composite using tissue mass fractions (Layer C). Compute uncertainty envelope and bound-shift sensitivity for the uncharacterized remainder.
- **Deliverables:** `_amino_acid_standard_whole_body.tsv`; whole-body Match Rate; sensitivity analysis showing worst-case shift under alternate lumped-remainder assumptions.
- **Runnable at any point** after at least one Step 7d category is complete (or possibly after Step 7b if the aggregation is written to accept blood alone). Each rerun improves coverage; the version label records which categories were included.
- **Handoff:** `docs/handoffs/08_whole_body_composite_handoff.md` (to be written).

### Step 9 — Paper — **WRITEABLE AT ANY POINT AFTER STEP 6**
- **Goal:** Methods paper.
- **Deliverables:** `docs/methods.md` assembled from per-step docs; figures from Steps 4/6/7/8; Gorissen reconciliation as a named subsection; supplementary data = the config files themselves.
- **A claim Step 6 put on the table:** under the new reference tryptophan limits half of all foods, and it is the amino acid acid hydrolysis destroys — a sequence-derived standard carries exactly the value a measured one cannot. The paper states it with the caveat that tryptophan is also the one scored amino acid without a measured cross-check until one is transcribed.
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
| Whether to run Match Rate v1.0 against muscle-with-metabolites first, or wait for whole-body composite | **Resolved Step 6 (2026-09-17): muscle with non-protein metabolite pools is v1.0's reference (D82).** |
| Reference-body citation for tissue mass fractions (ICRP 89 or Snyder ICRP 23 or equivalent) | Step 8 (author verification required) |
| Primary quantitative proteomics source for blood tier | **Resolved Step 7b part 1 (D101):** bryk_2017 (erythrocytes), geyer_2016 Table S2 (plasma); cross-checks gautier_2018 and hortin_2008 |
| Which muscle-specific assumptions actually create friction for other categories (input to Step 7c refactor scope) | Step 7b part 2 — the adaptation list; already known from part 1: tier names and the tier/subtype meaning (D97), per-compartment contaminants (D98), abundance kinds (D105), the fiber-type mix loader, header text |
| Hemoglobin's share of erythrocyte protein (53 / ~70 / ~95 % by three methods) and albumin's share of plasma protein (24 % LFQ vs >50 % immunoassay) — whether the sensitivity's movement justifies a two-piece erythrocyte compartment | Step 7b part 2 sensitivity (D102); follow-up decision if the movement is large |
| Which files of `outputs/standard/` go to `uncertainty/`, `stress/`, `sensitivity/` | Step 7a — a placement table of every file, approved by the author before the code is changed |
| Primary quantitative proteomics source for liver tier (Kim 2014 vs Wang 2019 vs Human Protein Atlas) | Step 7d-liver |
| Collagen tissue mass fractions (skin, tendon, ligament, bone matrix, fascia — lumped vs subdivided) | Step 7d-collagen (working assumption: lumped, per this session) |
| Match Rate scope: EAA-only vs all AA | **Resolved Step 6 (D86): EAA-only, the nine FAO 2013 headings each reduced to its indispensable member; conditionally indispensable amino acids become a separate scored set for a later output.** |
| Whether the muscle standard's sum-row display inconsistency gets fixed | **Resolved: fixed in code; the R9 text was stale.** |
| Fiber-type mix — accept Mpampoulis 2021 (female-only, %CSA) with disclosed limitation, or search for a mixed-sex source | **Closed 2026-09-17: accepted as D72, limitation disclosed in methods.** |
| Stopping criterion for adding tissues to whole-body composite | Step 8 sensitivity analysis; provable bound on remainder-shift |
| Isoform, processing, and PTM bounds — material after weighting? | Resolved rev. 5; carried to Step 4 as bounds |
| Glycosylation mass (unbounded from the entry) | Resolved rev. 5; stated as a limitation |
| Which lab AAA datasets of human muscle to validate against | Was Step 5; deferred with Step 5 |
| Match Rate delivery layer: proportion × absorption × burn (Match Rate Plus) | After Step 6; distinct from the reference itself. The spreadsheet's 95 % / 99 % absorption rows belong here, and the blend script (D88) is where they would first be used |
| The Match Rate formula: basis, normalization, and how surplus above the limiting amino acid enters the score | **Resolved Step 6 (D87):** shares over the scored set; the smallest ratio is the score; the surplus arithmetic reduces to it exactly; protein content never enters. The four spreadsheet foods are the transcription test |
| Which archive wins when a food appears in both Foundation Foods and SR Legacy | Open, deferred by the author (2026-09-18) — U6 and M4 carry both rows; every row names its archive |
| A minimum protein content for the ranked tables | Open, deferred by the author (2026-09-18) — the score is a proportion, so peaches score 86 % on 0.9 g protein; lean: a floor in `config/match_rate.ini` applied to the ranked view only, the per-food tables complete. A published zero scores 0 % with the limiting amino acid named (M3), as it should |
| Consolidation of USDA preparations whose panels differ within rounding (five rabbits at 92.83–92.92 %; sunflower kernels at 54.5717 / 54.5717 / 54.5707 %) | Open, pinned by the author (2026-09-18) — a later Match Rate revision; the amino acids survive simple processing, so one food should appear once |
| A tryptophan-preserving measurement of human muscle, to cross-check the one scored amino acid the Gorissen comparison cannot | Open, pinned by the author (2026-09-18) — transcribed into the comparison stage with a D-number when found; no change to the standard |
| The reference-foods table (whey, pea isolate, egg, steak, the current blend; old score beside new) | Open — whey, egg and steak are USDA rows and already in `match_rate_by_reference.tsv`; pea isolate and the blend need rows in `config/food_amino_acids_other_sources.csv` with their source |
| The leucine correction factor (concentration density; taxes and tolls) | With the blend script (D88), after Step 6 |
| Whether USDA's cystine figure is already expressed as cysteine equivalent | Open, not blocking — needs the FoodData Central documentation, not an assumption. Cysteine appears for only 25 of 5,156 foods; cystine is how USDA reports the sulphur amino acid for essentially the whole database. Nothing in the indispensable set depends on it; any all-twenty output does |
| The 1.54 pp of histidine, and the phenylalanine, serine, tyrosine and threonine disagreements that remain after the Gorissen comparison | Open — candidates are the blood category, the protein set's coverage, and the measurement itself; recorded in `docs/gorissen_comparison.md`, no verdict |
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

**New in R10 (all in `docs/decisions.md`):** **D82** the metabolite tier is named `non_protein_metabolite_pools` and the name is carried into the output file and column names, with a three-state definition in methods — "bound" was rejected because it already means *protein*-bound in this literature; **D83** the food side is USDA FoodData Central, CSV, hand-downloaded, SR Legacy and Foundation Foods only, zips gitignored and the manifest committed, the twenty joined by IUPAC trivial name so no amino acid is named in code or config; **D84** the calculated standard is compared with Gorissen et al. 2018's measured human muscle composition and nothing else, as a stage rather than a document, with no verdict; **D85** the paper's printed zeros for proline and cysteine are not measurements and are excluded from both sides, rule X5 rewritten so an undeclared zero stops the stage, and a mass balance (X6) computed every run.

**New in R15 (all in `docs/decisions.md`):** **D95** blood = plasma + erythrocytes, leucocytes/platelets bounded; **D96** blood owns plasma proteins; **D97** one tier × two subtypes, muscle names left for 7c; **D98** pools as measured, per-compartment contaminants; **D99** no metabolite pool; **D100** the split from measured reference values (ICRP volumes, NORIP concentrations, male midpoints); **D101** sources and roles; **D102** dominant-protein sensitivity, no adjustment; **D103** immunoglobulins by constant region with an A8 bound; **D104** F7 amended (`file.<n>.name`) and the source acquisition protocol; **D105** abundance kind declared per source.

**New in R12 (to be logged in `docs/decisions.md`):** **D89** — the repository reorganization (`outputs/intermediate/`, the `outputs/standard/` subfolders, `outputs/match/<reference>/`, the literature `used_by` line and map) and the bookkeeping items pinned in Steps 4b and 6 are executed as Step 7a, before the blood category; the abstraction stays after blood as Step 7c. A scoping refinement of D71, not a reversal: the abstraction is what needs a second category first; the layout does not. Step letters renumbered: blood 7b, refactor 7c, category loops 7d. D71's own text keeps the R8 letters (7a blood, 7b refactor, 7c loops) as logged.

**New in R11 (all in `docs/decisions.md`):** **D86** the scored set is FAO 2013's nine headings, each group reduced to its indispensable member (Met for SAA, Phe for AAA), cysteine and tyrosine unscored; **D87** the Match Rate formula is the author's spreadsheet formula transcribed as it is — shares over the scored set, the smallest ratio the score, its amino acid the limiting one; the surplus arithmetic reduces to the same number; protein content never enters; the leucine-target variants not carried; the original reference (Gorissen 2018, eight amino acids) carried beside the new; **D88** two scripts — the public single-food scorer in this repository, the blend / fortification script separate, importing the scorer, and private at publication.

**Logged in R9 (all in `docs/decisions.md`):** D67–D71 written into the log (they had been plan-only); **D72** the fiber-type mix (Mpampoulis et al. 2021, %CSA, female); **D73** inventory frame, resupply a sensitivity; **D74** per-fiber-type single-fiber values, no whole-muscle cross-check; **D75** carnosine only; **D76** protein-only table untouched, adjusted table separate; **D77** Mingrone et al. 2001 for protein and water per kg muscle, measured over computed; **D78** green light before calculation; **D79** assistant-extracted transcription with page-level locations, author-audited; **D80** fetch-by-code retired (7b); **D81** retrieval times out of generated config.

---

## 9. Thread Protocol

1. Open a new thread named for the step (e.g., "Step 4b — Muscle metabolite adjustment," "Step 7b — Blood category").
2. Paste the handoff doc as the first message.
3. Work only that step. Anything out of scope goes into the handoff's "Deferred / raised" section, not into the work.
4. On completion: update `decisions.md`, write the next step's handoff (or update the draft), tick this plan's step status and add a row to 10, commit.
5. Return to the planning thread only to re-plan, not to do work.
6. For Step 7b and Step 7d category expansions: follow the [Per-Category Gameplan](per_category_gameplan_R2.md) as the executing template. The gameplan itself is not re-revised for each category (though it may be revised after Step 7b+7c close, informed by what those steps taught). The gameplan's text still uses the R8 letters (7a blood, 7b refactor, 7c loops) until R3; its step 2 line `data/literature/<category>/` is superseded — literature stays one folder per source with `used_by` (Step 7a).

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
| 9 | 2026-09-17 | **Step 4b marked done** (D72–D81): `non_protein_metabolite_pools` stage, `_calculated_amino_acid_standard_with_bound_pools.tsv`, histidine +15 % in the standard. §1 scope wording updated (inventory, not turnover; anserine not carried). Step 6 marked next with its handoff drafted. §5 Step 7b gains the fetch retirement, the source-list split, per-stage hashing, and the run-order cheat sheet. §7 carnosine row resolved; turnover-range row added. |
| 10 | 2026-09-17 | **Step 6 in progress** (D82–D85), three deliveries. Reference tier renamed `non_protein_metabolite_pools` with a three-state definition in methods (D82). USDA food side built: `pipeline/usda.py`, rules U1–U7, 5,156 foods with amino acid panels (D83). Gorissen comparison built as a stage with rules X1–X6 and written up in `docs/gorissen_comparison.md` (D84); its printed zeros for proline and cysteine excluded and a mass balance added after they were shown not to be measurements (D85). §5 Step 6 restated as executed with what remains; Step 5's Gorissen row marked delivered; Step 7b gains the placeholder-convention fix and a concrete case for per-stage hashing. §7: four rows closed, five opened (formula basis, archive precedence, minimum protein content, the cystine convention, the residual Gorissen disagreements). §11 marks `gorissen_comparison.md` written. |
| 12 | 2026-09-18 | **Step 7a (reorganization) inserted before blood (D89); step letters renumbered** — blood 7b, refactor 7c, category loops 7d (live sections only; §8 and §10 entries keep the letters they were written with). §3 tree updated: `outputs/intermediate/`, the `outputs/standard/` subfolders, `outputs/match/<reference>/`, `config/literature_sources.ini`, `data/literature/manifest_map.tsv`, `docs/run_order.md`; `pipeline/` and `requirements.txt` in place of the stale `src/muscle_aa/` and `pyproject.toml`. §5 Step 7a written as a nested list; Step 7c's bookkeeping items moved to it. §7 placement row. §8 D89. §9 gameplan letters noted. §11 the 7a handoff and `run_order.md`. Nothing in Steps 1–6 changes. |
| 15 | 2026-09-18 | **Step 7b part 1 closed** (D95–D105): scope, the ten blood sources on disk and hashed, the split from NORIP + ICRP, the dominant-protein finding in both compartments, immunoglobulins, the F7 amendment and the source acquisition protocol (gameplan R3 §2a). Delivery 1 R2 (`fetch_literature.py`). §5 Step 7b restated as part 1 / part 2; §7 rows resolved; §11 artifacts. |
| 14 | 2026-09-18 | **Step 7a marked done** (D90–D94). `outputs/` reorganized into `intermediate/` and a subfoldered `standard/`, with `match/<reference>/`; placement held as data in `common.STANDARD_SUBFOLDERS`. Generated paths resolved at call time, not at import, with `tests/test_paths.py` to keep them that way (D91). A stage hashes only the config sections and manifest entries it reads, timestamps excluded; A0 amended (D92). `config/literature_sources.ini` with `used_for` and `data/literature/manifest_map.tsv` (D93). `fetch-literature` its own command, never downloads; F3/F4/F4b/F6 retired, F8/F9 added (D94). §3 tree replaced with the tree as built. §5 Step 7a restated as executed; Step 7b marked next. §11 gains `run_order.md` and `step7a_report.md`. Nothing in Steps 1–6 changes; no number changed. The older stages' header rewrite moved to Step 7c, to be done once with the header templating rather than twice. |
| 13 | 2026-09-18 | **Errata only — no scope, status, or decision changes.** Two live sections still carried the pre-D89 step letters. §Status: "Step 7b (reorganization) is next, then Step 7c (blood), then Step 7d (the abstraction)" corrected to 7a reorganization → 7b blood → 7c abstraction → 7d loops, matching §5. §5 Step 4b, "Rules changed in this step": the fetch-by-code retirement (D80) said "pinned for 7c"; it was moved to Step 7a by D89, as §5 Step 7c already states. Issued at the opening of the Step 7a thread rather than at its close, so the thread runs against correct letters; the Step 7a close therefore issues R14. |
| 11 | 2026-09-18 | **Step 6 marked done** (D86–D88; deliveries 3 and 4 R2). The formula transcribed from the author's spreadsheet and shown to reduce to the limiting ratio (D87); the scored set (D86); the two-script decision (D88); `pipeline/match.py`, rules M1–M6, six output tables including the hand-validation table `match_rate_steps.tsv`. The run's findings recorded (tryptophan limits half of foods under the new reference; phenylalanine and histidine no longer limit). §7: formula and scope rows resolved; archive precedence and the protein floor deferred by the author; rows added for the tryptophan cross-check, the consolidation of near-identical USDA preparations, the reference-foods table and the leucine correction factor. §5 Step 7b gains the per-stage excerpt spec and the header rule; Step 9 gains the tryptophan claim. Standing rules from this step (in `docs/conventions.md`): generated headers are provenance only; every reported number is rebuildable by hand from its row. Step 7a marked next with its handoff written. |
| 8 | 2026-09-13 | **Pipeline refactor sequenced after first non-muscle category (D71):** Step 7 restructured into 7a (blood, first category using muscle-adapted code), 7b (pipeline refactor codifying what 7a taught), and 7c (subsequent categories using modular code). §2 design principles gain "rule-of-three abstraction" bullet. §3 architecture note that aggregate is tier-generic post-D71. §5 Step 7 subdivided with acceptance criteria for refactor (byte-identical muscle rebuild, blood rebuild matches Step 7a output). §7 open questions gain row for "which muscle-specific assumptions actually create friction" (input to 7b scope, recorded during 7a). §11 references gain `adding_a_category.md` as 7b artifact. Per-category gameplan referenced as R2 (updated in parallel to reflect the after-blood refactor timing). |

---

## 11. Referenced Artifacts

- [`per_category_gameplan_R3.md`](per_category_gameplan_R3.md) — **R3 2026-09-18:** adds §2a, the source acquisition protocol (who does what, in what order) and the `file.<n>.name` rule. The twelve-step template for adding any new tissue/category to the whole-body composite. Referenced by Steps 7b, 7d, and 9. May be revised (R3+) after Step 7b+7c close, informed by what those steps teach about the actual per-category friction points.
- [`planning_summary_expanding_to_full_body_match_R3.md`](planning_summary_expanding_to_full_body_match_R3.md) — the reframe planning output that motivated R7. Contains the Gorissen comparison table, carnosine math walk-through, collagen EAA composition, cumulative body composition sketch, and per-tissue EAA ranking. Referenced by Steps 4b, 6, and 9.
- [`docs/gorissen_comparison.md`](gorissen_comparison.md) — **written 2026-09-17 (D84, D85).** Reconciles the standard against Gorissen et al. 2018's measured human muscle composition. Every figure comes from `outputs/comparison/`; carries the mass balance that showed the paper's protein content to be inflated by non-protein nitrogen, and the finding that human muscle has the highest histidine of the sixteen protein sources in the paper's own table.
- `docs/step6_report.md` — **written 2026-09-18.** The Step 6 narrative: procedure, the run's results, the tryptophan finding, method lessons.
- `docs/formula.md`, "Match Rate" — the calculation by hand with a worked example (pea against the original reference) checkable cell by cell against the author's spreadsheet.
- `docs/handoffs/07a_layout_handoff.md` — **written 2026-09-18.** The reorganization as a nested list, the bookkeeping items, the acceptance criteria.
- [`docs/run_order.md`](run_order.md) — **written 2026-09-18.** One screen: the eight commands in order with what each reads and writes, and what to rerun after any edit.
- `docs/step7b_report.md` — **part 1 written 2026-09-18;** the build appends part 2.
- `docs/handoffs/07b_blood_handoff.md` — **R2 2026-09-18:** the build handoff for the new thread — everything settled in part 1, the run order for blood, the adaptation list started.
- `docs/step7a_report.md` — **written 2026-09-18.** The Step 7a narrative: the placement table, the acceptance check over 111 files, the three bugs that were the same bug, the open items.
- `docs/adding_a_category.md` (to be written in Step 7c) — developer-facing map showing where in the code a new category is added; complements the per-category gameplan (which is planning-facing) with the code-facing perspective.
