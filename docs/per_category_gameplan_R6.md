# Per-Category Gameplan

**Revision:** R6 (2026-09-20) — **step 11 corrected to the composite stage as built.** `pipeline/composite.py` shipped inside Step 7b (D121–D124) between R4 and R5, and R5's step 11 still described the composite as something Step 8 would pick up later; it is command 10 of the run and it runs every time. Step 11, the deliverables checklist and this line are the only changes. R5 (2026-09-20) — **first revision informed by an executed category.** Blood (Step 7b) ran the twelve steps end to end; what it taught is folded in: the steps that turned out to be new work rather than reused code, the two join patterns a category can have, the adaptation list consolidated in full as Step 7c's specification, and time estimates replaced with actuals. R4 (2026-09-19) — step letters renumbered to the project plan's (7b blood, 7c refactor, 7d loops). R3 (2026-09-18) — the source acquisition protocol (§2a), after the blood thread lost hours to the assistant writing file URLs it had not been given.
**Established by:** D70 (Project Plan R7, §8); refactor sequencing set by D71 (R8, §8).
**Referenced from:** Project Plan §5 Steps 7b and 7d, §9 Thread Protocol.
**Purpose:** Template for adding any new tissue/category to the whole-body composite. Each execution produces one publishable category standard that consumes pipeline infrastructure and feeds the whole-body composite. This gameplan is not re-revised for each category — the category itself is the deliverable, and its handoff document instantiates the template. R5 is the revision R2 anticipated: blood has executed, and the template now says what the first execution actually cost and where it diverged.

---

## When to use this

**Use this gameplan for:**
- Adding a new tissue/category tier (blood, liver, collagen tissues, brain, heart, etc.)
- Building a category standard that will feed the whole-body composite
- Any category ranked in the estimated-EAA-contribution table from `planning_summary_expanding_to_full_body_match_R3.md`

**Do NOT use this gameplan for:**
- Modifications to a category that has already been shipped (those get their own decisions and a handoff, not a template re-run)
- The whole-body composite itself (Step 8 in the project plan — the close-out that decides the categories are enough and bounds the remainder; it consumes category outputs and does not follow this template). The composite *stage* is not Step 8's to build: it exists, and step 11 below is how a category enters it.
- Metabolite pool adjustments to an existing category (small enough to be a handoff on their own; e.g., Step 4b muscle metabolite adjustment did not use this template)
- The pipeline refactor for per-category modularity (Step 7c in the project plan; it is infrastructure work, not a category build)

---

## Pipeline state to know about

The tier-generic pipeline refactor (Step 7c, D71) lands **after** the first non-muscle category (Step 7b: blood) is built. Executions of this gameplan therefore split into two eras:

- **Pre-refactor (Step 7b only — blood; executed 2026-09-19/20).** The pipeline assumed muscle-specific structure in more places than the pre-build list predicted: twenty-two of them, now enumerated below. Four of the twelve steps needed new modules or new rules rather than parameters. Blood's commands live beside muscle's (`blood-protein-set`, `blood-composition`, `blood-mass-fractions`, `blood-standard`) rather than inside them; muscle rebuilds byte-identical. This is intentional (rule-of-three abstraction, D71) and the twenty-two items are Step 7c's specification.
- **Post-refactor (Step 7d onward — liver, collagen, and beyond).** The pipeline is category-agnostic. The reused-code steps reuse code; a new category adds config, and its standard drops into `outputs/<category>/`.

**If a second category is built before 7c** (a real possibility if liver is wanted sooner), it reuses blood's two modules rather than muscle's when its shape matches: `mass_fractions_pools.py` if its dataset publishes a mass share per protein, `aggregate_pools.py` if its compartments are mixed by measured masses. What is genuinely per-category is config: the pool decisions, the column maps by letter, the mix. Expect a day, not blood's four.

### The two shapes a category can have

Blood proved the pipeline's real axis of variance is not tissue — it is **how the primary dataset identifies and quantifies a protein.** Fix these two before any code:

| | Muscle (Step 3b) | Blood (Step 7b) |
|---|---|---|
| **Identity** | gene symbol (Murgia published no accessions) → rules M1–M5 | published accession, as a MaxQuant protein group → rules P1–P6, with UniProt's secondary accessions (P3b) and the gene cell (P4b) as fallbacks |
| **Quantity** | iBAQ-like value × molecular weight → mass share (B2) | a mass share as published; no MW step (D105, `share_kind = mass_share`) |
| **Subtype** | one dataset with a column per fiber type, mixed by a fiber-type mix | one dataset per compartment, mixed by the compartments' protein masses (D100) |
| **Spread** | the dataset's own SD column per fiber type | per-replicate columns (mean and SD across donors, D108) or per-run columns (technical, D114) |

A new category answers those four rows in its scoping delivery. If both the identity and quantity rows match blood's, the existing modules run with new config only.

## Prerequisites (all true after muscle build shipped)

- UniProt fetch code with MD5 verification (Step 2 infrastructure)
- Composition engine producing both mass conventions (Step 2 infrastructure)
- Aggregate stage (category-generic post-Step 7c; muscle-adapted plus blood's `aggregate_pools.py` pre-Step 7c)
- Uncertainty stage (log-space Monte Carlo, D63; blood's is inside `aggregate_pools.py`)
- Stress stage (named perturbations, D65; a category whose sensitivities already cover its scenarios may skip it, D120)
- Additive sequence store and per-category raw entries (S1, D110) — a category's fetch never touches another's
- `docs/conventions.md` (config schema, citation schema, mass conventions)
- `PROVENANCE.md` standing rule enforced

---

## The Twelve Steps

Each step names its inputs, outputs, reused code, and where new work sits. Steps marked "reused code" run cleanly post-refactor (Step 7d onward) and require adaptation pre-refactor (Step 7b only).

### 1. Pick the category and its scope
**New work.** Decide what tissue(s) the category represents, what proteins are in-scope, what the boundary decisions are (contaminants, dual-role proteins between categories, tissue-vs-cell scope). One D-number for this scoping decision, referencing D67 (whole-body composite framing). Written down in the handoff's opening section.

**Deliverable:** scoping paragraph in the handoff. D-number logged in `docs/decisions.md`.

**Examples:**
- Blood → hemoglobin subunits + serum albumin + immunoglobulins + complement + minor plasma proteins. Boundary: is myoglobin in muscle-cell or blood category? (Muscle-cell, since it's synthesized and stored in the myocyte.)
- Liver → hepatocyte proteome. Boundary: are Kupffer cells and other non-hepatocyte cells included? (Working assumption: include all liver-resident cell types, weighted as measured by the primary dataset.)
- Collagen tissues → skin + tendon + ligament + bone matrix organic fraction + fascia, dominated by Type I collagen. Boundary: is cartilage in this category or its own? (Working assumption: separate if ever added — cartilage has Type II + aggrecan.)

### 2. Find the primary quantitative proteomics dataset
**New work.** Search for the peer-reviewed, quantitative (iBAQ × MW or equivalent) proteomics dataset with sufficient coverage of the category. Prefer human, prefer recent, prefer datasets used in other published references. If no single dataset covers, name the fallback: multiple sources with a stated join rule, or a well-characterized single-protein dominance (e.g., hemoglobin at ~65 % of blood protein).

**Author verification required:** fetch the article, verify the author list is not a search-snippet artifact.

**Deliverable:** cited primary dataset + any cross-check datasets. Referenced in the handoff. Downloaded supplementary tables in `data/literature/<category>/`, SHA-256 pinned.

**Candidate datasets by category (starting points, not commitments):**
- Blood: hematology-standard hemoglobin fraction + published plasma proteomics (Anderson & Anderson-type surveys)
- Liver: Kim et al. 2014 (Nature) draft human proteome, Wang et al. 2019 (Mol Syst Biol) quantitative tissue proteomes, or Human Protein Atlas tissue expression
- Collagen tissues: already in-project via the collagen track; extend with tissue-specific collagen-to-non-collagen ratios if needed
- Brain: Kim 2014, Wang 2019, or brain-specific proteomics (multiple candidates)

### 2a. Source acquisition protocol (R3) — who does what, in what order

The point of this protocol is that nothing about a source is ever typed from memory. Each step is done by the party who can verify it directly; the other party never does it for them.

| # | Step | Who | Verified against |
|---|---|---|---|
| 1 | Identify candidate sources; state the role each would fill (primary / cross-check / split / reference) and why | Claude | search results and pages fetched and quoted in the thread |
| 2 | Obtain the **article URL** (landing page or DOI), the citation, and the author list | Claude | the page or PDF itself; the author list is re-verified from the PDF once it is on disk |
| 3 | Open the article, locate every file to be stored (supplementary tables, article PDF, figures), **download each one, and record its URL exactly as the browser shows it** (download history or right-click → copy link) | Anthony | the browser |
| 4 | Place each file in `data/literature/<source_id>/` **under the name it was downloaded with**; no renaming, ever | Anthony | — |
| 5 | Fill `file.<n>.url`, `file.<n>.name` (the name as downloaded), `file.<n>.published_name`, `file.<n>.obtained = manual` in `config/literature_sources.ini`; leave `sha256` blank | Anthony (Claude may draft the section text with the URL and name fields left empty) | — |
| 6 | Run `python run.py fetch-literature` — it hashes what is on disk, records blank hashes, verifies pinned ones, and regenerates `manifest.ini`, `manifest_map.tsv`, `README.md` | Anthony | the run log; a flag is the only signal that something is off |
| 7 | Read the stored files and report every value to be used with its page / table / column location; propose roles and D-numbers | Claude | the files on disk |

What Claude does **not** do: write a file URL it was not handed; build a URL from a publisher's pattern; ask for a file to be renamed; pre-fill a hash; declare a file "probably" named as served. What Anthony does **not** do: transcribe values out of a paper (that is step 7, and the values come with locations).

**Rule change carried to `docs/conventions.md` (F7 amended, D104):** the stored filename is the name the file was downloaded with, recorded in config as `file.<n>.name`; `file.<n>.url` is provenance only and is never used to derive a name. Since D80 nothing downloads, so a name derived from a URL served no purpose and failed on real links (ACS supplement links end in `/`; browsers rename on save). `fetch-literature` finds the file by `file.<n>.name` and falls back to the URL rule only for the existing sources that predate the key.

**Naming a source:** `<first author>_<year>`, all lower case; a second paper by the same author in the same year gets a letter (`smith_2011b`).

**Time cost when the protocol is followed:** steps 1–2 an hour; steps 3–6 fifteen minutes per source; step 7 an hour per source. When it is not followed, the blood thread's cost was a working afternoon.

### 3. Find or write the accession list
**New work.** Turn the primary dataset's identifiers into a UniProt accession list. Which rules apply depends on the identity row of the table above.

- **Published accessions (blood's case, rules P1–P6, D107).** The dataset's identity cell is a protein group — accessions joined by `;` in the authors' order. The entry is the first token that is a reviewed human entry, isoform suffix stripped. Two fallbacks are mandatory, not optional: a published accession may since have been **merged or demerged** (P3b: ask UniProt by the `sec_acc` query field, one token per query, and confirm against the entry's own `secondaryAccessions`), and a group may list only unreviewed fragments while the reviewed entry sits outside it (P4b: resolve the row's gene cell through M2/M3, with M3b breaking a synonym collision in favour of the one primary-symbol entry). Skipping the fallbacks cost blood 4.1 % of plasma mass; with them, 0.35 %.
- **Gene symbols (muscle's case, rules M1–M5).** Unchanged. Note item 19: M2 compares the whole primary-symbol field, and an entry with several gene names carries them `;`-joined (`HBA1; HBA2`), so it never matches — a muscle-side fix for 7c.
- **Contaminants** are per pool and cited: a deterministic rule the code can apply to any dataset (blood used UniProt's keyword KW-0416 for keratins, D109, because the MaxQuant `CON__` marker differed between the two datasets' searches) plus the genes the paper itself names, matched on the dataset's own gene cell. Nothing is excluded silently: every row's outcome and every removed row's share are written.

**Sizes, for planning:** blood's two pools were 322 and 2,653 rows → 283 and 2,437 entries. "A small category has dozens of proteins" was wrong — a modern proteomics table has thousands of rows whatever the tissue, and the pool is what the dataset quantified (D98), not what the textbook lists.

**Deliverable:** `config/<category>/protein_set_decisions.ini` (hand-written: the pools, their datasets, the contaminant rules with page locations) → `config/<category>/accessions.ini` and `segments.ini` (generated, D61), plus the outcome and excluded-share tables under `data/<category>/`.

### 4. Harvest sequences from UniProt
**Reused code; additive (S1, D110).** The sequence store is one file shared by every category. A category's fetch asks only for the accessions the store lacks, renders those sections, and carries every other section over verbatim; the raw entry JSON goes to `data/uniprot_raw/<category>/`, named for the fetch that brought it in, and nothing is ever deleted. A muscle refetch cannot drop a blood entry and vice versa.

**Deliverable:** sections added to `data/uniprot_sequences.ini`; `data/uniprot_raw/<category>/`. Blood added ~1,400 entries in about 45 minutes.

### 5. Compute residue compositions
**Reused code, `--category <name>`.** Both mass conventions; the PTM disclosure runs beside it. One rule earned here: features of the same type at the same position are **alternatives**, not additions (C3b, D118) — UniProt records each glycoform seen at a site as its own feature, and summing them put IGHG1's annotated PTM mass at 113 % of its molecular weight.

**Deliverable:** `outputs/<category>/intermediate/composition/`.

### 6. Assign mass fractions from the primary dataset
**New work (config); reused code if the shape matches.** The join was already done by Step 3 — the enumeration's row-outcome table is the join table, and this step never re-joins. What it needs is the **column map by letter** (`[columns.<source>.<file>]`: sheet, header row, identity, gene, the value columns, the replicate columns): blood's tables carried duplicate and stray header names, so a header string cannot identify a column and a name-based map is a silent wrong answer waiting to happen.

Rules Q1–Q7 (`mass_fractions_pools.py`) cover a published-share pool: summed rows (D47), shared rows split in proportion to the members' own rows with the whole to each `_high` (D48), the weight within pool (D31), the spread from replicates (D108) or runs (D114), what the rules removed and what no entry answers to as two separate shares (Q4, A7), and the cross-checks the category has (Q7). A category whose dataset is iBAQ-like uses muscle's `mass_fractions.py` B-rules instead.

**Every category needs its own bound for what its dominant entries do not measure.** Blood's was the immunoglobulin variable domains (D113): a constant-region entry is quantified by its constant domains, and the molecules carry variable domains no row measures — bounded by w × L_V / L_C from the pool's own V-region entries, 6.1 % of plasma against a loose 15.6 %. The pattern generalizes: find the entries whose measurement covers part of a molecule, bound the rest from the pool's own numbers, never from an assumption.

**Deliverable:** `config/<category>/mass_fraction_decisions.ini` (hand-written, `[scope] category = <name>`, the maps by letter with page locations) → `config/<category>/mass_fractions_per_entry.tsv` (generated, D61) and `outputs/<category>/intermediate/mass_fractions/`.

### 7. Aggregate to category-level standard
**New module pre-refactor; reused post-refactor.** Rules T1–T8 (`aggregate_pools.py`) cover a category whose compartments are mixed by measured masses: the profile as the w-weighted mean of per-protein mass fractions (no molar step when the weights are already mass), the mix computed from cited reference quantities (T2), the dominant-protein and per-replicate sensitivities (T3, T4), the Monte Carlo (T5), every bound on one scale (T6), and the per-amino-acid drivers (T7).

**The mix is arithmetic on cited quantities, never a ratio someone published.** Blood's (D100, D119): plasma protein = plasma volume × plasma protein concentration; haemoglobin = red-cell volume × MCHC; erythrocyte protein = haemoglobin ÷ haemoglobin's share of erythrocyte protein; the second route (blood haemoglobin × blood volume) reported as a cross-check because it differs by 4 %. Write the arithmetic into `config/<category>/<mix>.ini` with a page location per number, and have the stage compute any share a cited number implies rather than typing a rounded one — a number without a location does not go into config (D119).

**Expect the dominant protein to be under-weighted.** Label-free MS compressed both of blood's: haemoglobin 53 % measured against 70 % (a second MS anchor) and ~98 % (implied by ICRP's whole-blood protein); albumin 24 % against 51 % (clinical concentrations). The rule is D54's — weights as measured, the dominant protein a named sensitivity across the measured shares, no adjustment — and the sensitivity must carry **both** consequences: the compartment's profile and the compartment's protein mass (D117).

**Deliverable:** `outputs/<category>/standard/_calculated_amino_acid_standard.tsv` (+ the residue convention), the profiles, the EAA subset, the mix table, the sensitivities, the drivers, plots.

### 8. Uncertainty and stress on the category standard
**Reused pattern; inside the aggregate module pre-refactor.** The D63 Monte Carlo is unchanged in kind: a log-normal fitted to each entry's (v, sd), one term drawing the spread and one drawing σ/√n, every draw recomputing the whole profile.

**Stress may be redundant, and saying so is the deliverable.** Muscle's scenarios (TPA weighting, completeness, knockouts) exist to ask "what if the weighting is wrong". Where the category already answers that with named sensitivities — blood's dominant-protein and per-replicate tables — a stress stage only renames them, and it is skipped with a stated reason (D120). Decide this explicitly; do not carry the stage because muscle had it.

**Put every term on one scale** (T6): Monte Carlo half-widths, the weighted isoform/processing/PTM bounds, the category's own unmeasured-mass bound, the completeness gap, and each sensitivity's spread, all as a maximum absolute shift of the g/100 g value, with the widest named per amino acid. That table is what a reader of the standard needs and the only place the terms are comparable.

**Deliverable:** `outputs/<category>/standard/uncertainty_intervals.tsv`, `uncertainty_per_amino_acid.tsv`, the sensitivity spreads, and — if stress is skipped — the reason in the summary.

### 9. Identify and cite category-specific metabolite pools (if any)
**New work.** Not every category has meaningful bound metabolite pools. For those that do, produce a `config/<category>/metabolite_pools.ini` following the D68 template (pool size + turnover rate + amino acid stoichiometry). Fold into the category standard as a tier adjustment.

**Deliverable:** `config/<category>/metabolite_pools.ini` (only if applicable); updated category standard with metabolite adjustment noted in header.

**Known metabolite pools by category:**
- **Muscle:** carnosine + anserine (→ His). Handled in Step 4b, not repeated per gameplan.
- **Blood:** none (D99, confirmed in the build — free amino acids in plasma are transient and excluded by D68). Glutathione is carried as an all-twenty extension, not an EAA-scope item.
- **Liver:** glutathione (γ-Glu-Cys-Gly) at high concentration — relevant for Cys, Gly, Glu (only Gly among non-essentials matters for EAA scope; Cys is non-EAA in Match Rate).
- **Brain:** neurotransmitter-derived pools (GABA from Glu, serotonin from Trp, dopamine from Tyr, glycine as neurotransmitter). Trp and Tyr contributions are EAA-relevant.
- **Heart, kidneys, other:** creatine + phosphocreatine (→ Gly) in tissues with high energy demand; small vs muscle's contribution.

### 10. Publish category standard as its own file
**Deliverable.** The category standard is a first-class deliverable (D67 fork framing), published under `outputs/<category>/standard/` with the header naming its inputs' hashes and the D-numbers that scope it.

**The category's own Match Rate is additive, by construction.** Since Step 7a the match stage reads references from `config/match_rate.ini`: a category becomes a `[reference.<name>]` section naming its standard table, scored into its own `match/<reference>/` folder with a column beside the existing ones. Nothing of an earlier category's output is overwritten, so before / after is a diff rather than an argument. Blood's reference landed in the thread that followed the build.

### 11. Add the category to the composite
**Part of the category build, not a later step.** `pipeline/composite.py` is command 10 of the run (`python run.py composite`, D123), between `blood-standard` and `usda`. It mixes every category standard named in `config/tissue_mass_fractions.ini` by the protein mass that category holds in the reference person (Layer C, D67 / D121; rules L1–L4 in `docs/conventions.md`) and writes `outputs/composite/`. It refuses to run on fewer than two categories. Adding a category means adding its section and rerunning from command 10 to the end.

**What to write.** One `[category.<n>]` section. The keys every category needs:

- `label` — what the category's standard contains, in words; it is printed in the composite's header
- `standard_file` + `standard_column` — the category's free-convention standard, the profile that gets mixed
- `residue_file` + `residue_column` — its residue-convention table, read only for the amino-acid-mass line (L4)
- `tissue_mass_g`, `tissue_mass_g_female`, `tissue_mass_location` — the reference body's mass for the tissue, both sexes, with the page and table (blood: 5,600 g M / 4,100 g F, ICRP 89 p. 18 Table 2.8, row "Blood"). The male is the base; the female is recorded and reported as a line, never the base (L2, the D100 / D119 pattern)
- `protein_mass` — **which of the two rules states the category's protein mass.** This is the row to settle before writing anything else:
  - `tissue_mass_x_protein_content` — tissue mass × a cited protein content per kg. Needs `protein_content_from` (`<path> [<section>]`, so the value is transcribed once and read, never re-typed) and `protein_content_basis`; the stage stops if the section's own `basis` disagrees. Muscle's case.
  - `category_split` — the category's own stage already computed the protein mass, and the profile and the mass must be taken at the same point. Needs `split_file`, `split_column`, `split_base_row`, `split_female_prefix`. Blood's case (D122): the whole-blood protein at the base haemoglobin share.

Optional, where the category has them:

- `non_protein_metabolite_pool_amounts` + `fiber_type_mix` — a category carrying a non-protein metabolite pool adds the pool's amino acid mass to its amino-acid-mass weight (L4). Muscle only, so far.
- `sensitivity_file` and its five column keys (`sensitivity_row_column`, `sensitivity_row_share_column`, `sensitivity_variant_column`, `sensitivity_variant_share_column`, `sensitivity_profile_column`) — a category whose split has published alternative shares is reported at every one of them, its mass moving with the share (L3, D122). Blood's D117 table.

**Also update `config/match_rate.ini`.** A category adds **two** `[reference.*]` sections, not one (D124): its own standard, and the composite. The composite's reference id names the categories in it (`skeletal_muscle_with_non_protein_metabolite_pools_and_blood`), so it changes when a category is added — the section name, the `[reference.*]` file pointer if the column moved, and any `difference_pp_` column built from it. `primary_reference` is not changed by a category build: it stays the muscle reference (D82, D124), so each new category's effect is readable against one unchanging baseline.

**Deliverables:** the `[category.<n>]` section; `outputs/composite/` rebuilt (the composite standard, `category_protein_masses.tsv`, the contributions, the EAA subset, the split and weighting-basis sensitivities, `composite_summary.ini`); the two match references. Verify the tissue mass against the reference-body PDF before the run — a wrong mass produces a plausible composite rather than an error.

### 12. Log decisions, update handoff, tick the plan
**Delivery.** Every category-specific decision goes to `docs/decisions.md` with date, category, and rationale. The next handoff (either the next category or Step 8 whole-body composite) is written or updated. Project plan §5 Step 7b or 7d gets a status update. Commit.

**Decisions go in `docs/decisions.md` and rules in `docs/conventions.md` — never in a handoff or a delivery note.** The blood thread drifted on this because `docs/` was not in the thread's first tarball; the fix is to ask for the file, not to ship rows elsewhere. A thread's first tarball carries `docs/`, `pipeline/`, `tests/`, `config/`, `run.py`.

**For a pre-refactor category:** the handoff also enumerates the muscle-code adaptations required. Blood's list is consolidated below and is Step 7c's specification; a second pre-refactor category appends to it.

---

## The adaptation list (Step 7b output — Step 7c's specification)

Every muscle assumption blood met, with where it was met and what was done. Items 1–3 and 6–10 were predicted before the build; 11–22 were found during it. A second pre-refactor category appends rows here.

| # | Muscle assumption | Where | What blood did | 7c |
|---|---|---|---|---|
| 1 | tiers are named `contractile` / `builders` | `common.TIER_NAMES`; `aggregate.py` `PROFILES` and its three tier dicts, header strings | blood has one tier and two subtypes; names never entered the muscle code (blood's own module) | names from a per-category manifest |
| 2 | subtypes are `I` / `IIa` / `IIx` fiber types | `aggregate.py` `FIBER_TYPES`, `GEL_TYPES`, `load_fiber_type_mix()`, `fiber_type_mix.ini` | `plasma` / `erythrocytes` in `compartment_mix.ini`, read by `aggregate_pools.py` | one subtype-mix loader that reads names and its own arithmetic |
| 3 | "tier" = same-dataset weighting, "subtype" = cross-dataset mixing | the whole of `aggregate.py` | stated in the decisions (D97); blood's pools are both | say it in the manifest, one vocabulary |
| 4 | Tier 1 seeded from Gene Ontology, Tier 2 the measured remainder by gene symbol | `enumerate_pool.py` `enumerate_published_pool`; `build_protein_set.py` `category_membership_text` | pools are the quantified set, joined by published accession (P1–P6); no `resolved_terms.ini`, R3 does not apply | the pool rule as a per-category choice, not a branch |
| 5 | one contaminant list per category (`contaminant_source`, M5) | `enumerate_pool.py` (P5); `config/blood/protein_set_decisions.ini` `[pool.*] named_contaminant_genes`, `[contaminant_rules]` | per-pool rule: a UniProt keyword on the entry plus the genes the paper names, on the dataset's gene cell | contaminant rules as a list of cited predicates |
| 6 | the abundance column is iBAQ (× MW) | `mass_fractions.py` B-rules | `share_kind = mass_share` per source, used as-is (D105); the `molar` branch is not built until a source needs it | the value rule declared per column map |
| 7 | the primary dataset has a per-subtype column set | the column maps | each subtype is its own dataset and its own map | maps keyed by (category, subtype) |
| 8 | the MHC:actin and band-family sensitivities | `stress_test_settings.ini`; `aggregate.py` A6 | the dominant-protein share sensitivity in both compartments (D102, D117) | a per-category sensitivity list |
| 9 | header text says "fiber", "muscle", "contractile" | every writer | blood's writers say pool and compartment | templated headers with category substitutions |
| 10 | `comparison` is Gorissen-only | `comparison.py` | no external blood profile exists; the Hortin comparison lives in the mass-fraction stage (D116, D120) | comparison sources per category, or absent |
| 11 | a source's `role` is global | `mass_fractions.py` scanned every source in the shared file | `role` scoped to a category (**D106**) — found by a *muscle* test when blood's two primaries stopped the muscle stage | scoping is already in; keep the test |
| 12 | the sequence store is rebuilt from one pool and the raw folder emptied every fetch | `fetch_sequences.py` `read_store_text`, `merge_header`, `main` | merge-write, additive, nothing deleted (S1, **D110**) | unchanged; this is the post-refactor behaviour |
| 13 | the store's per-entry `tiers` line is a muscle fact in a shared file | `fetch_sequences.py` `to_section` | left as is; membership lives in the category's `accessions.ini` | per-category membership out of the shared store |
| 14 | `accessions.ini` field `tier` means tier | `build_protein_set.py` | for a category it holds pool names; not renamed | rename to `pools`, one migration |
| 15 | every stage's paths are import-time constants | `common.category_files()`; each stage's `--category` | resolved at call time from the category | the category is an argument of every stage; muscle's constants fold into the same map |
| 16 | the raw entries live in one flat folder | `common.uniprot_raw_dir()` / `uniprot_raw_entry()`; `ptm_disclosure.py`; `build_protein_set._uniprot_mol_weight` | one folder per category; readers look in the root and every subfolder | unchanged |
| 17 | composition and PTM disclosure read the flat config and write the flat outputs | `common.category_files()`; `composition.py`, `ptm_disclosure.py` `--category` | the category's files in, `outputs/<category>/intermediate/composition/` out; masses, symbols and ptmlist stay shared | unchanged |
| 18 | `load_masses()` bound its path at import | `composition.py` | resolved at call time (D91) | audit for other import-time binds |
| 19 | M2 compares the whole primary-symbol field | `enumerate_pool.py` `enumerate_measured_tier` `attribute()` | **muscle bug, left untouched for byte-identity**: an entry with several gene names (`HBA1; HBA2`) never matches, which is why HBA1 sits outside muscle's measured pool; blood's P4b splits the field | fix in 7c with a muscle rerun and a diff of what enters the pool |
| 20 | `mass_fractions.py` is one dataset × fiber types × tiers with a MW step, digest families and the muscle cross-checks | `pipeline/mass_fractions_pools.py` (whole module) | a second module for pools from published shares; no digest for blood (its families are the published groups, D115) | unify: the join table as input, the value rule per map, the cross-check list per category, digest optional |
| 21 | `literature_inventory.py` runs only in the muscle chain and does not open `.xls` / `.doc` | — | still true; Hortin's files inventory as `unknown_not_opened`, and the mass-fraction stage reads the `.xls` through xlrd | inventory per category; xlrd path shared |
| 22 | `aggregate.py` mixes fiber-type columns by a fiber mix; `uncertainty.py`, `stress.py`, `plots.py` are fiber-shaped | `pipeline/aggregate_pools.py` (whole module) | one module for pools mixed by measured masses, with its own Monte Carlo and plots | unify: the mix as config with its own arithmetic, the sensitivity list per category, stress optional (D120) |

**Acceptance for 7c is unchanged** (project plan §5): muscle rebuilds byte-identical, blood rebuilds and matches this step's output, tests pass parametrized. Item 19 is the one exception — fixing it *changes* muscle's pool, so it needs its own decision and a before/after table rather than a byte-identity check.

---

## Time estimates (post-muscle, per category)

The muscle build (Steps 1–4) took approximately 20 hours of active work across 3 days, including inventing the pipeline pattern. **Blood's actuals replace the pre-build estimates for the pre-refactor column.**

| Step | Effort (post-refactor, estimate) | Blood, actual (pre-refactor) |
|---|---|---|
| 1. Scope decision | 1–2 hours | ~half a day (part 1, with the split and the two findings) |
| 2. Find primary dataset | 2–4 hours | a working day, an afternoon of it lost before §2a existed |
| 3. Accession list | Hours to a day | **a day and three runs** — the two fallbacks (P3b, P4b) were each found by a run that had already shipped |
| 4. Harvest sequences | Minutes (automated) | 45 minutes of network, once |
| 5. Compute residue compositions | Minutes (automated) | minutes, plus C3b |
| 6. Assign mass fractions | Hours to a day | ~a day: a new module, the column maps by letter, the Ig bound, three cross-checks |
| 7. Aggregate to category standard | Minutes (automated) | ~a day: a new module, the mix config, the two sensitivities |
| 8. Uncertainty and stress | Minutes (automated) | folded into step 7; stress skipped with a reason |
| 9. Metabolite pools (if any) | 2–4 hours | none (D99) |
| 10. Publish category standard | Minutes | minutes; the Match Rate reference in the following thread |
| 11. Add tissue mass fraction | Minutes | open at close |
| 12. Log decisions and update plan | 30 minutes | ~2 hours, including this list |

**Blood, total: about four days of active work across three threads** — against the 2–4 days estimated, so the estimate held, but the shape did not: steps 3, 6 and 7 were the cost, and steps 4, 5, 8 were free.

**Realistic total per category post-refactor: 1–3 days.**
**A second pre-refactor category that matches blood's shape: 1–2 days** (config, not modules).

Where the time actually goes, in order: reading the primary table well enough to write a column map by letter; the identity rules and their fallbacks; the category's own unmeasured-mass bound; the mix arithmetic and its citations. Where it does not go: fetching, composition, the Monte Carlo.

## Deliverables checklist per category

Paths are the as-built layout (Step 7a, D89–D94; blood, Step 7b).

- [ ] every literature file obtained, placed and hashed per §2a — no URL in config that was not taken from the browser, no file renamed
- [ ] `config/<category>/protein_set_decisions.ini` — the pools, their datasets, the contaminant rules, each with a page location
- [ ] `config/<category>/mass_fraction_decisions.ini` — `[scope] category`, the column maps **by letter**, the cross-checks, the category's unmeasured-mass bound
- [ ] `config/<category>/<mix>.ini` — the mix arithmetic, a location per number (if the category has subtypes)
- [ ] `config/<category>/accessions.ini`, `segments.ini`, `mass_fractions_per_entry.tsv` — generated, covered by a currency test (D61)
- [ ] `data/<category>/` — pool tables, every row's outcome, the removed rows with their shares
- [ ] `outputs/<category>/intermediate/{protein_set,composition,mass_fractions}/`
- [ ] `outputs/<category>/standard/_calculated_amino_acid_standard.tsv` (+ residue convention), profiles, EAA subset
- [ ] `outputs/<category>/standard/uncertainty_intervals.tsv`, `uncertainty_per_amino_acid.tsv` — every term on one scale
- [ ] the category's sensitivities (dominant protein, replicates) and its drivers table
- [ ] stress outputs **or** a stated reason they are redundant (D120)
- [ ] `config/tissue_mass_fractions.ini` — the `[category.<n>]` section: the `protein_mass` rule and its keys, the tissue mass both sexes with its page location
- [ ] `outputs/composite/` rebuilt with the category in it (command 10), its share of the composite checked against the masses
- [ ] **two** sections in `config/match_rate.ini` (D124) — the category's own reference, and the composite reference with its id updated — each scored additively into `match/<reference>/`
- [ ] an excerpt spec per command in `pipeline/excerpt.py`
- [ ] decisions in `docs/decisions.md`, rules in `docs/conventions.md` — nowhere else
- [ ] `docs/handoffs/<step>_<category>_handoff.md` updated at close; `docs/<step>_report.md` written
- [ ] project plan §5 status updated

**Additional for a pre-refactor category:**
- [ ] its rows appended to the adaptation list above

## Notes

**On category-specific Match Rate:** each category standard supports its own Match Rate variant — same math, different reference — and since Step 7a the mechanism is additive: a `[reference.<name>]` section, its own `match/<reference>/` folder, a column beside the existing ones. The paper (Step 9) can present multiple framings; the RFP response uses one clearly labeled version.

**On what a category is likely to teach.** Blood's three transferable findings: (1) a published identifier is a label at a point in time, and a 2016 accession list needs UniProt's merge history to read in 2026; (2) label-free MS compresses whatever dominates the sample, so a category with a dominant protein needs its share as a named sensitivity carrying both the profile and the mass; (3) the honest bound for mass a method does not measure comes from the pool's own entries, not from an assumption. Expect the next category to add a fourth.

**On the stopping criterion for adding tissues:** the goal is not to characterize every tissue. It is to characterize enough tissues that the uncharacterized remainder can be bounded — mathematically, under worst-case assumptions about its composition — to within a stated tolerance on the composite standard. Once the bound is provably small (e.g., <0.5 pp on any EAA, <10 % of that EAA's value), further category builds are optional. See Project Plan R8 §5 Step 8.

**On maintenance:** if a category standard needs revision after shipping (e.g., a better primary dataset appears, a metabolite pool citation changes), the revision does not re-run this gameplan. It is a modification to an existing category, handled with its own handoff and D-number.

**Revision log.** R2 (planning thread, before the blood build): the twelve steps and the two eras. R3 (2026-09-18, Step 7b thread): §2a source acquisition protocol added; the deliverables checklist gains its first line; the `file.<n>.name` rule recorded here pending its D-number. R4 (2026-09-19, Step 7b thread): the step letters throughout renumbered to D89's. **R5 (2026-09-20, after the blood build closed):** the first revision informed by an executed category — the two shapes a category can have (identity × quantity) added to §Pipeline state; steps 3–8 rewritten with what blood actually required, including the mandatory identity fallbacks, the column-map-by-letter rule, the category's own unmeasured-mass bound, the mix-as-arithmetic rule and the dominant-protein expectation; step 10 gains the additive Match Rate mechanism; step 12 gains the decisions-live-in-decisions.md rule; **the adaptation list consolidated in full (items 1–22) as Step 7c's specification**; time estimates replaced with blood's actuals; the deliverables checklist rewritten against the as-built layout; "blood: dozens" corrected (2,720 entries). **R6 (2026-09-20, Step 7b debug thread):** step 11 rewritten against `pipeline/composite.py` as built — the composite is command 10 of the run, not something Step 8 picks up later; the deliverable is a `[category.<n>]` section with one of two `protein_mass` rules, not a row; the match-reference line corrected to two sections per category (D124); the deliverables checklist updated to match. R5 folded in delivery 10's match references but not the composite stage from the same delivery. Also: §2a's F7 note said "pending its D-number" — the number is D104 (2026-09-18) and the line now says so.

**On this gameplan's own future revisions:** R5 folds in Step 7b. The next revision is due after **Step 7c** closes: the adaptation list becomes a record rather than a specification, the two-era framing collapses, and the pre-refactor columns in the time table come out. A second pre-refactor category before 7c appends to the adaptation list without a revision.
