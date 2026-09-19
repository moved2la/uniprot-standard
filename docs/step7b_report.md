# Step 7b report — Blood (part 1: scope, sources, decisions; part 2: the build)

**Thread:** 2026-09-18. **Decisions:** D95–D105. **Deliveries:** 1 (R2), 2 (documents).
**Result of part 1:** blood is scoped as plasma + erythrocytes; ten sources are on disk, hashed
locally, `fetch-literature` green (34 files, 29 sources × 5 columns); the plasma : erythrocyte
split is specified from measured reference values; the two findings that will shape the build are
on record; F7 is amended and the source acquisition protocol is written down.
**Part 2 (the build) is appended below. Part 3 (analysis, Match Rate, the Layer C row) follows from its own thread.**

## What the step set out to do

Build the first non-muscle category with the muscle-adapted code and keep a list of every place a
muscle assumption has to bend — the list is the input to Step 7c (D71). Part 1 is everything that
has to be decided before a line of the build runs.

## Procedure

1. **Scope first** (gameplan step 1). Six questions from the handoff plus one it did not list
   (immunoglobulins), each with a recommendation; the author decided all seven in one pass
   (D95–D99, D103). The tier / subtype mapping was drawn out as a table before it was agreed:
   in the code, tiers are same-dataset weighting and subtypes are cross-dataset mixing by a
   cited share — so blood is one tier × two subtypes, and "tier" never meant "kind of protein".
2. **Sources** (gameplan step 2). Candidates from search, verified on their pages; the author
   downloaded and placed every file and recorded every URL; the script hashed. Hortin 2008 was
   bought ($43) after Anderson 2002 turned out to carry its concentrations only as a figure.
   NORIP's two papers were added when ICRP 89 turned out to carry no plasma protein value.
3. **The tables read** (D79 pattern: values from files on disk, with locations). Bryk's Table S3
   is Murgia-shaped — a measured mass share per donor, four donors. Geyer's Table S3 is CVs only,
   so plasma's primary is Table S2 (one donor, fifteen replicates). Hortin's `.xls` is the
   paper's database: 153 rows with accession, gene, mean mg/L and a reference each.
4. **Decisions written** as a batch (D95–D105), two carrying a recommendation the author then
   confirmed (male midpoints; Geyer S2 as plasma primary).

## Results

**The split, all inputs measured** (D100). ICRP 89: red cells 2300 ml, plasma 3000 ml, blood
5600 g (male). NORIP: plasma total protein 64–79 g/L (rustad_2004, Table I p. 274, n = 877);
hemoglobin 134–170 g/L and haematocrit 0.395–0.500 in men, MCHC 317–357 g/L (nordin_2004,
Table V p. 394), pooled means Hb 142.9 g/L, EVF 0.424, MCHC 337 g/L (Table III p. 389). NORIP's
MCHC agrees with ICRP's 33.5 g/100 ml within one per cent, and its hemoglobin lands on ICRP's
red-cell route (770 g), not its whole-blood route (875 g) — which is what settled the two-figure
inconsistency found in ICRP.

**Finding 1 — the dominant protein is under-weighted by label-free MS in both compartments**
(D102). Hemoglobin's share of erythrocyte protein: 53 % by bryk_2017's TPA (median of four
donors), about 0.7 by gautier_2018's MS signal (Fig. 2A), about 95 % by classical chemistry.
Gautier says of Bryk outright (p. 2655) that the globin signals are underestimated. Albumin's
share of plasma protein: 24 % by geyer_2016's LFQ, above half by hortin_2008 (p. 1608; 51 % of
its summed table). Same shape twice. The rule is D54's: weights as measured, the dominant
protein a named sensitivity across the measured shares, no adjustment; the sensitivity table
will say how much the profile moves, and a two-piece erythrocyte compartment is the follow-up if
it moves a lot.

**Finding 2 — immunoglobulins are the first place Layer A is not exact** (D103). MaxQuant
quantifies them as constant-region protein groups because variable-region peptides cannot be
assembled; Hortin says the same of its concentrations (p. 1614). Constant regions as measured,
variable domains as an A8 bound — a few per cent of blood protein, under one per cent at
whole-body weighting.

**Finding 3 — F7 failed on a real link** (D104). ACS supplement URLs end in `/`, so "the last
URL segment" is empty; browsers rename on save. Since D80 nothing downloads, so the URL had no
business naming the file. `file.<n>.name` is now the stored name; the URL is provenance.

**Finding 4 — Anderson 2002 is not a cross-check table.** Its 289-protein list has no
concentration column; the values are bars in Figure 3. Hortin 2008 is the table.

**Three abundance kinds in one category** (D105): TPA fraction (mass share), LFQ intensity
(mass-proportional share), µmol/L × Mr and copies per cell (molar). Declared per source in config,
applied in code — the adaptation list's first entry that is a rule rather than a rename.

## Decisions and rules

D95 (what blood is), D96 (boundary with liver), D97 (one tier × two subtypes; muscle names left
for 7c), D98 (pools as measured; per-compartment contaminants), D99 (no metabolite pool),
D100 (the split), D101 (sources and roles), D102 (dominant-protein sensitivity), D103
(immunoglobulins), D104 (F7 amended; source acquisition protocol, gameplan R3 §2a), D105
(abundance kind per source).

## Method lessons

- **A deterministic pipeline is only as deterministic as its inputs' provenance.** The assistant
  wrote two file URLs from a publisher's pattern, pre-filled hashes it had computed on uploaded
  copies, and twice asked for files to be renamed to fit URLs it had guessed. Each looked like
  helpfulness; each broke the rule the project was founded on. The fix is a protocol with names
  in it (gameplan R3 §2a), not a resolution to be careful.
- **Measured over computed, again.** ICRP 89 has no plasma protein value; summing Hortin's
  153 means gives 85.7 g/L, above every clinical interval, and the paper's own "top dozen are
  95 %" does not hold in its own table (77 %). The measured number existed in NORIP and cost
  two downloads. D77's pattern held.
- **Read the supplement before the paper's number.** Bryk's copies-per-cell rest on an assumed
  90 fL and 20 % protein; the TPA fraction beside it is the measurement. Geyer's Table S3 title
  promises the ten individuals; the sheet holds CVs.
- **A delivery's caller list is a test's job.** Delivery 1 grew a tuple by one field and missed
  one of two callers; R2 added the end-to-end test that would have caught it.
- **Ask about the thing that is actually open** — the tier / subtype question was fuzzy until it
  was drawn as a two-by-two of muscle against blood.

## Open (documented, not acted on)

- The build (part 2) and the adaptation list it produces.
- Whether D102's sensitivity movement justifies a two-piece erythrocyte compartment.
- The ten-individual plasma set (PXD002854 `proteinGroups.txt`) and Gautier's per-protein table
  (PXD009258 `txt.zip`, 4 GB): upgrades, not fetched.
- Six of Hortin's 153 rows (immunoglobulin heavy chains) carry no accession; the cross-check
  joins the 147 that do.
- Glutathione for an all-twenty extension (D99).
- Everything carried from Steps 6 and 7a, untouched here.


---

# Part 2 — the build (2026-09-19/20)

**Thread:** build thread, part 2. **Decisions:** D106–D120. **Deliveries:** 6, 7, 7b, 7c, 8, 8b, 9.
**Result of part 2:** the four blood commands exist — `blood-protein-set`, `blood-composition`,
`blood-mass-fractions`, `blood-standard` — and run green from a clean start after the muscle
chain (commands 1–8, then the four, then `test`). The blood standard is on disk at
`outputs/blood/standard/_calculated_amino_acid_standard.tsv` with its uncertainty, sensitivity and
driver tables. Muscle rebuilds byte-identical throughout, except the PTM disclosure tables under a
decided rule change (D118). No analysis of the numbers is in this part (deferred to part 3).

## Procedure, delivery by delivery

| Delivery | Command / change | Decisions | What the run showed |
|---|---|---|---|
| 6 | `blood-protein-set`: pools by published accession (P1–P6), additive sequence store, raw entries nested per category, `config/blood/protein_set_decisions.ini` and `mass_fraction_decisions.ini` (column maps by letter) | D107–D110 | fetch ~45 min; 247 tests. Plasma 4.45 % of the listing excluded — mostly accessions UniProt had demerged or merged since 2016, and TrEMBL-only groups |
| 7 | `blood-composition` (the muscle stages with `--category`); D111's two fallbacks in the enumeration | D111 | the fallbacks did not fire: UniProt answers deleted / demerged accessions with inactive stubs; an entry's several gene names come `;`-joined in the primary field |
| 7b | P3b reads inactive stubs; P4b per gene; M3b; the multi-gene primary field split | D112 | crashed: `sec_acc` is a query field, not a return field |
| 7c | P3b one token per query by `sec_acc`, confirmed against the entry JSON; a refused field falls through to P4b | — | plasma exclusion 0.35 %, erythrocytes 0.79 %; 38 plasma tokens resolved as secondary accessions (the old immunoglobulin V-region entries had been *merged*, not deleted); 96 erythrocyte rows by gene |
| 8 | `blood-mass-fractions` (`mass_fractions_pools.py`): weights from published shares, bounds, the D113 immunoglobulin bound, the Hortin and deep-dataset cross-checks, the overlap line; `[cross_checks]`, `[immunoglobulin_bound]`, S5 and Hortin column maps | D113–D116 | ran after `blood-composition` had actually been run (it had been skipped twice); albumin 24 % vs 51 %; Ig bound 6.1 % vs loose 15.6 %; IGHG1's PTM sum 113 % of its mass |
| 8b | C3b in `ptm_disclosure.py` | D118 | shared code; muscle's two PTM tables move where an entry has alternatives at one site |
| 9 | `blood-standard` (`aggregate_pools.py`): profiles, the split, D117 in both compartments, per-donor, D63 Monte Carlo, one scale, drivers, plots; `config/blood/compartment_mix.ini` drafted | D117, D119, D120 | green; analysis deferred |

## The run record (transcribed from the excerpts; not analysed here)

**Pools** (`data/blood/pool_queries.ini`, delivery 7c run)

| | plasma (geyer_2016 S2) | erythrocytes (bryk_2017 S3) |
|---|---|---|
| rows | 322 (321 quantified) | 2,653 |
| entries | 283 | 2,437 |
| member rows: first token / later token / secondary accession / secondary shared / gene / gene shared / duplicate | 258 / 5 / 13 / 1 / 0 / 0 / 26 | 2,114 / 156 / 3 / 0 / 96 / 2 / 170 |
| excluded share of the listing | 0.35 % (HBB, HBA1, HBD, CA1 = 0.33 %; keratins 0.006 %; one orphon V-gene) | 0.79 % (keratins 0.78 %; ALB 0.016 %; four rows < 0.01 %) |
| tokens inactive in UniProt / asked as secondary / resolved | ~30 / 38 / 38 | — / 15 / 15 |
| entries in both pools | 62 | |

**Protein set** (`outputs/blood/intermediate/protein_set/`): 9 entries excluded under R5 (selenoproteins:
GPX1 0.05 % of erythrocytes, SELENOP in plasma, seven others); 126 multi-chain entries under D24
(complement C3 with twelve Chain features); the immunoglobulin constant-region entries carry
`Chain 1(OUTSIDE)–end`, resolved whole by R2e, consistent with D103; no open flags.

**Mass fractions** (`mass_fractions_summary.ini`)

| | plasma | erythrocytes |
|---|---|---|
| denominator as share of the listing | 0.9964 | 0.9914 |
| largest entry | ALB 24.0 % | HBB 26.8 % |
| ten largest | 68 % | 74 % |
| zero-weight members (ghost-only rows) | 0 | 547 |
| completeness gap (A7) | 0.02 % | 0.06 % |
| immunoglobulin variable-domain bound (D113) | 6.1 % (loose 15.6 %; L_V ≈ 97 residues from 13 / 13 / 11 V entries) | ~0 |
| Hortin cross-check | 118 accessions joined covering 85.7 % of the pool; albumin share 0.585 (Hortin, joined) vs 0.281 (Geyer, joined), ratio 0.48; APOA1 ×4.9, SERPINA1 ×4.3, A2M ×2.8, C3 ×2.2, fibrinogen ×2–2.7; largest profile difference C +0.015, then K, A +0.008, Q −0.008, G, S −0.007 | — |
| deep dataset (S5) signal absent from S2 | 0.6 % (657 rows) | — |
| overlap sum of w | 62.8 % (albumin and the large plasma proteins, present at ~0 in the erythrocyte table) | 5.5 % (PRDX2 and other erythrocyte proteins present at ~0 in plasma; not carry-over) |

**Standard** (`outputs/blood/standard/`): built and green; `standard_summary.ini` `[split]` and
`[histidine]`, `sensitivity_dominant_protein_spread.tsv`, `histidine_drivers.tsv` are the tables
part 3 reads first. The split's arithmetic on the D100 numbers, before the profiles are read:

| haemoglobin share of erythrocyte protein | erythrocyte protein (775 g Hb ÷ share) | plasma : erythrocyte protein | whole-blood protein |
|---|---|---|---|
| 0.53 (bryk_2017, as measured) | 1,462 g | 13 : 87 | 1,677 g |
| 0.70 (gautier_2018) | 1,107 g | 16 : 84 | 1,322 g |
| ~0.98 (implied by ICRP's 1,008 g whole-blood protein) | ~790 g | 21 : 79 | 1,008 g |

Plasma protein is 214.5 g in every row (3.0 L × 71.5 g/L).

## Decisions and rules

D106–D120 in `docs/decisions.md`; rules P1–P6, P3b, P4b, M3b, S1, C3b, Q1–Q7, T1–T8 in
`docs/conventions.md`; the R4 amendment (the `tier` field holds pool names for a category).

## Method lessons

- **UniProt is not static.** Of Geyer's 2016 accessions, ~30 were inactive by 2026: demerged
  (P0CG05 → IGLC2 + IGLC3), merged (the old immunoglobulin V-region entries into the current IGKV/IGLV
  entries), or TrEMBL fragments. The join key the dataset supplies still needs the "formerly known
  as" (P3b) and, failing that, the dataset's own gene cell (P4b). Both are deterministic; nothing was
  hand-mapped. Plasma's exclusion went from 4.45 % to 0.35 %.
- **Read the field list before writing the query.** `sec_acc` is a query field with no return
  field; a search for a deleted accession returns an inactive stub, not nothing; an entry's gene
  names come `;`-joined in one field. Each cost one run. Muscle's M2 has the same `;` blind spot
  (item 19), left as is for byte-identity.
- **The same compression twice.** Label-free MS under-weights the dominant protein of each
  compartment — haemoglobin (D102) and now albumin (Hortin: 51 % against LFQ's 24 %) — and ICRP's
  whole-blood protein reproduces only at a haemoglobin share near 0.98. Three independent lines say
  the as-measured share is low. The base stays as measured; the sensitivity carries both
  consequences (D117). The discussion is part 3's.
- **A number without a location does not go into config.** The "95 % by classical chemistry" of
  part 1 had none; the share ICRP's cited numbers imply took its place (D119).
- **Shared code vs the blood branch.** A delivery that touches `common.py`, `fetch_sequences.py`,
  `ptm_disclosure.py` or `run.py` reruns from the top; a delivery confined to a blood stage does not.
  Each delivery states which.
- **Run lists are numbered steps, one per line.** A prose "run" paragraph cost two skipped commands.
- **Decisions go in `docs/decisions.md`.** They drifted into handoff notes because `docs/` was not in
  the thread's tarball; the fix was to ask for the file, not to ship rows elsewhere.

## Errata carried

- D95 and D103 cite "(A8)" for the bound rule; A8 is the EAA rule, A7 the bound. Fix in `decisions.md`.
- Part 1's "about 95 % by classical chemistry" is uncited (D119 replaces it).
- The gameplan's "blood: dozens" (§3) is stale — blood is 2,653 + 322 rows.

## Open (documented, not acted on) — part 3

- The analysis of the standard: the split, histidine, whether the as-measured base stays.
- Blood in the Match Rate outputs as an additive reference (before / after comparable).
- The Layer C row: `config/tissue_mass_fractions.ini`, blood 5,600 g M / 4,100 g F.
- Run-order placement of the four blood commands (provisional: after `match`).
- The adaptation list (items 4–22, in `docs/handoffs/step7b_delivery*_notes.md`) consolidated
  into the gameplan for 7c.
- `docs/pipeline_map.md`, `docs/run_order.md` brought current.
- One `___` in `config/blood/mass_fraction_decisions.ini`: the article page defining Bryk's TPA
  fraction (file.3 is now in the manifest).
- The female split needs no more transcription (ICRP's female volumes are in the config); it is
  reported, not the base.
- The carry-over line in `pool_overlap_mass_share.tsv` is directional and should be restated as
  the plasma-side sum of entries at ~0 in the erythrocyte table.
- Everything carried from part 1, untouched here.
