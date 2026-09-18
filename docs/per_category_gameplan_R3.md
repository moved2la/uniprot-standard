# Per-Category Gameplan

**Revision:** R3 (2026-09-18) — adds the source acquisition protocol (who does what, in what order) after the blood thread lost hours to the assistant writing file URLs it had not been given. R2 text is otherwise unchanged.
**Established by:** D70 (Project Plan R7, §8); refactor sequencing set by D71 (R8, §8).
**Referenced from:** Project Plan R8 §5 Steps 7a and 7c, §9 Thread Protocol.
**Purpose:** Template for adding any new tissue/category to the whole-body composite. Each execution produces one publishable category standard that consumes pipeline infrastructure and feeds the whole-body composite. This gameplan is not re-revised for each category — the category itself is the deliverable, and its handoff document instantiates the template. The gameplan itself may be revised (R3+) after Step 7a and 7b close, informed by what those steps teach.

---

## When to use this

**Use this gameplan for:**
- Adding a new tissue/category tier (blood, liver, collagen tissues, brain, heart, etc.)
- Building a category standard that will feed the whole-body composite
- Any category ranked in the estimated-EAA-contribution table from `planning_summary_expanding_to_full_body_match_R3.md`

**Do NOT use this gameplan for:**
- Modifications to a category that has already been shipped (those get their own decisions and a handoff, not a template re-run)
- The whole-body composite itself (Step 8 in the project plan; it consumes category outputs, does not follow this template)
- Metabolite pool adjustments to an existing category (small enough to be a handoff on their own; e.g., Step 4b muscle metabolite adjustment did not use this template)
- The pipeline refactor for per-category modularity (Step 7b in the project plan; it is infrastructure work, not a category build)

---

## Pipeline state to know about

The tier-generic pipeline refactor (Step 7b, D71) lands **after** the first non-muscle category (Step 7a: blood) is built. Executions of this gameplan therefore split into two eras:

- **Pre-refactor (Step 7a only — blood).** The pipeline still assumes muscle-specific column names ("contractile," "builders," "type_I/IIa/IIx") in the aggregate stage and other places. Steps 4, 5, 7, 8 of this gameplan (marked "reused code" below) will require ad-hoc adaptations rather than clean parameterization. Record what needed adapting in the handoff — those observations are the input to Step 7b's refactor scope. This is intentional (rule-of-three abstraction, D71).
- **Post-refactor (Step 7c onward — liver, collagen, and beyond).** The pipeline is tier- and category-agnostic. Steps 4, 5, 7, 8 of this gameplan reuse code cleanly with no muscle-specific adaptation; the new category just adds config and its standard drops into `outputs/<category>/`.

The gameplan's twelve steps are the same in both eras. What differs is how much muscle-specific code the reused-code steps have to work around.

---

## Prerequisites (all true after muscle build shipped)

- UniProt fetch code with MD5 verification (Step 2 infrastructure)
- Composition engine producing both mass conventions (Step 2 infrastructure)
- Aggregate stage (tier-generic post-Step 7b; muscle-adapted pre-Step 7b)
- Uncertainty stage (log-space Monte Carlo, D63)
- Stress stage (named perturbations, D65)
- `docs/conventions.md` (config schema, citation schema, mass conventions)
- `PROVENANCE.md` standing rule enforced

---

## The Twelve Steps

Each step names its inputs, outputs, reused code, and where new work sits. Steps marked "reused code" run cleanly post-refactor (Step 7c onward) and require adaptation pre-refactor (Step 7a only).

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

**Rule change carried to `docs/conventions.md` (F7 amended, pending its D-number):** the stored filename is the name the file was downloaded with, recorded in config as `file.<n>.name`; `file.<n>.url` is provenance only and is never used to derive a name. Since D80 nothing downloads, so a name derived from a URL served no purpose and failed on real links (ACS supplement links end in `/`; browsers rename on save). `fetch-literature` finds the file by `file.<n>.name` and falls back to the URL rule only for the existing sources that predate the key.

**Naming a source:** `<first author>_<year>`, all lower case; a second paper by the same author in the same year gets a letter (`smith_2011b`).

**Time cost when the protocol is followed:** steps 1–2 an hour; steps 3–6 fifteen minutes per source; step 7 an hour per source. When it is not followed, the blood thread's cost was a working afternoon.

### 3. Find or write the accession list
**New work.** From the primary dataset's gene/protein identifiers, produce a UniProt accession list for the category. Small categories (blood: dozens) are compact. Large categories (liver: thousands) reuse the pool-computation pattern from Step 1 of the muscle build.

**Deliverable:** `config/<category>/accessions.ini`, following the same schema as `config/skeletal_muscle/accessions.ini`.

### 4. Harvest sequences from UniProt
**Reused code, no changes.** Run `fetch.py` against the new accession list. MD5 verification per D25.

**Deliverable:** sequences added to `data/sequences.ini` (shared across categories).

### 5. Compute residue compositions
**Reused code, no changes.** Run `composition.py` on the new accessions. Both mass conventions computed.

**Deliverable:** `outputs/<category>/amino_acid_composition_per_protein.tsv`.

### 6. Assign mass fractions from the primary dataset
**New work.** For each accession, cite the row of the primary dataset that gives its weight in the category. Follow the same rules as Step 3b (D31 within-tier, D47 summed rows, D48 shared rows).

**Deliverable:** `config/<category>/mass_fractions.ini` (or per-subtype files if the category has subtypes analogous to muscle's fiber types). Every value has source, location, retrieved date per citation schema.

### 7. Aggregate to category-level standard
**Reused code — clean post-Step 7b refactor; ad-hoc adaptations pre-refactor.** Run the aggregate stage on the new category. Produces per-subtype and combined profiles in both conventions.

**Pre-refactor note (Step 7a only):** the aggregate stage assumes muscle-specific column names. Expect to fork or adapt the aggregate module. Record what adaptations were required in the handoff — that is the input to Step 7b's refactor scope.

**Deliverable:** `outputs/<category>/_calculated_amino_acid_standard.tsv` — the category standard in the same shape as the muscle standard.

### 8. Uncertainty and stress on the category standard
**Reused code — clean post-refactor.** Run uncertainty (log-space Monte Carlo, D63) and stress (D65) stages on the new category. The stress test scenarios apply per-category with sensible parameter ranges — knockout of top-N, tier ratio swings, size tilt, TPA weighting, random abuse.

**Pre-refactor note (Step 7a only):** may require similar adaptations as Step 7 above if the uncertainty/stress modules also assume muscle-specific tier structure.

**Deliverable:** `outputs/<category>/uncertainty_intervals.tsv`, `outputs/<category>/stress_summary.ini` and associated files, following the muscle-standard pattern.

### 9. Identify and cite category-specific metabolite pools (if any)
**New work.** Not every category has meaningful bound metabolite pools. For those that do, produce a `config/<category>/metabolite_pools.ini` following the D68 template (pool size + turnover rate + amino acid stoichiometry). Fold into the category standard as a tier adjustment.

**Deliverable:** `config/<category>/metabolite_pools.ini` (only if applicable); updated category standard with metabolite adjustment noted in header.

**Known metabolite pools by category:**
- **Muscle:** carnosine + anserine (→ His). Handled in Step 4b, not repeated per gameplan.
- **Blood:** no meaningful bound metabolite pools (free amino acids in plasma are transient and excluded by D68).
- **Liver:** glutathione (γ-Glu-Cys-Gly) at high concentration — relevant for Cys, Gly, Glu (only Gly among non-essentials matters for EAA scope; Cys is non-EAA in Match Rate).
- **Brain:** neurotransmitter-derived pools (GABA from Glu, serotonin from Trp, dopamine from Tyr, glycine as neurotransmitter). Trp and Tyr contributions are EAA-relevant.
- **Heart, kidneys, other:** creatine + phosphocreatine (→ Gly) in tissues with high energy demand; small vs muscle's contribution.

### 10. Publish category standard as its own file
**Deliverable.** The category standard is a first-class deliverable (D67 fork framing). Publish in `outputs/<category>/` with version-labeled header naming the D-numbers that scope it (D67, D68, plus category-specific decisions). Category-specific Match Rate can be computed against this standalone standard.

### 11. Feed category standard into whole-body composite
**Consumer step; separate from the category build.** Once the category standard exists, add its tissue mass fraction to `config/tissue_mass_fractions.ini` (with citation to the reference-body publication, D67). The whole-body composite (Step 8 in the project plan) picks it up on next run.

**Deliverable:** row added to `config/tissue_mass_fractions.ini`. Not a rerun of the whole-body composite — that is Step 8.

### 12. Log decisions, update handoff, tick the plan
**Delivery.** Every category-specific decision goes to `docs/decisions.md` with date, category, and rationale. The next handoff (either the next category or Step 8 whole-body composite) is written or updated. Project plan §5 Step 7a or 7c gets a status update. Commit.

**For Step 7a specifically:** the handoff must also enumerate the muscle-code adaptations that were required for Steps 7 and 8 above. That list is Step 7b's input.

---

## Time estimates (post-muscle, per category)

The muscle build (Steps 1–4) took approximately 20 hours of active work across 3 days, including inventing the pipeline pattern. Per-category work reuses that pattern.

| Step | Effort (post-refactor) | Effort (Step 7a, pre-refactor) |
|---|---|---|
| 1. Scope decision | 1–2 hours | 1–2 hours |
| 2. Find primary dataset | 2–4 hours | 2–4 hours |
| 3. Accession list | Hours to a day | Hours to a day |
| 4. Harvest sequences | Minutes (automated) | Minutes (automated) |
| 5. Compute residue compositions | Minutes (automated) | Minutes (automated) |
| 6. Assign mass fractions | Hours to a day | Hours to a day |
| 7. Aggregate to category standard | Minutes (automated) | Hours (ad-hoc adaptations) |
| 8. Uncertainty and stress | Minutes (automated) | Hours (ad-hoc adaptations) |
| 9. Metabolite pools (if any) | 2–4 hours | 2–4 hours |
| 10. Publish category standard | Minutes | Minutes |
| 11. Add tissue mass fraction | Minutes | Minutes |
| 12. Log decisions and update plan | 30 minutes | 1 hour (also enumerating adaptations for 7b) |

**Realistic total per category post-refactor: 1–3 days of active work.**
**Realistic total for Step 7a (blood, pre-refactor): 2–4 days, plus the enumeration of adaptations for Step 7b.**

Blood is likely the fastest category on scope (small proteome, well-characterized, one dominant protein — hemoglobin) but the slowest in wall-clock time because it runs pre-refactor. Liver runs post-refactor and benefits from the cleaner pipeline. Collagen tissues reuse the standalone collagen track's work and mostly need tissue mass fractions + lumping decisions.

---

## Deliverables checklist per category

At the end of each category execution, the following should exist:

- [ ] every literature file obtained, placed and hashed per §2a — no URL in config that was not taken from the browser, no file renamed
- [ ] `config/<category>/accessions.ini` — cited membership rule
- [ ] `config/<category>/mass_fractions.ini` (or per-subtype files) — cited per row
- [ ] `config/<category>/metabolite_pools.ini` — only if applicable
- [ ] `outputs/<category>/amino_acid_composition_per_protein.tsv`
- [ ] `outputs/<category>/_calculated_amino_acid_standard.tsv` — the category standard
- [ ] `outputs/<category>/uncertainty_intervals.tsv` and associated files
- [ ] `outputs/<category>/stress_summary.ini` and associated files
- [ ] Row added to `config/tissue_mass_fractions.ini` — cited tissue mass fraction
- [ ] Decisions logged in `docs/decisions.md`
- [ ] `docs/handoffs/07<a-or-c>_<category>_handoff.md` — the executed handoff, updated at close
- [ ] Project plan §5 Step 7a or 7c status updated

**Additional for Step 7a (blood only):**
- [ ] Handoff enumerates muscle-specific assumptions encountered and adaptations required — this becomes the specification for Step 7b's refactor

---

## Notes

**On category-specific Match Rate:** each category standard can support its own Match Rate variant. "Match to muscle" scores foods against the muscle standard. "Match to collagen tissues" scores foods against the collagen standard. "Match to whole body" scores against the composite. Same math, different reference. The paper (Step 9) can present multiple Match Rate framings if useful; the RFP response (Step 6) uses one clearly labeled version.

**On the stopping criterion for adding tissues:** the goal is not to characterize every tissue. It is to characterize enough tissues that the uncharacterized remainder can be bounded — mathematically, under worst-case assumptions about its composition — to within a stated tolerance on the composite standard. Once the bound is provably small (e.g., <0.5 pp on any EAA, <10 % of that EAA's value), further category builds are optional. See Project Plan R8 §5 Step 8.

**On maintenance:** if a category standard needs revision after shipping (e.g., a better primary dataset appears, a metabolite pool citation changes), the revision does not re-run this gameplan. It is a modification to an existing category, handled with its own handoff and D-number.

**Revision log.** R2 (planning thread, before Step 7a): the twelve steps and the two eras. R3 (2026-09-18, Step 7b thread): §2a source acquisition protocol added; the deliverables checklist gains its first line; the `file.<n>.name` rule recorded here pending its D-number and the `docs/conventions.md` edit. The step letters in R2's text (7a blood, 7b refactor) predate D89's renumbering (7b blood, 7c refactor, 7d loops) and are left as written.

**On this gameplan's own future revisions:** R2 was written before the first non-muscle category (blood, Step 7a) executed. After Step 7a and 7b close, this gameplan may benefit from R3 revisions informed by what those steps actually taught — particularly around Step 7 (aggregate) and Step 8 (uncertainty/stress) of the gameplan, which are the steps most likely to reveal muscle-specific assumptions the pipeline needs to shed.
