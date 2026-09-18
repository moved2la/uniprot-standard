# Step 6 report — Match Rate against USDA

**Threads:** 2026-09-16 to 2026-09-18 (two: "building step 6" for the food side and the Gorissen comparison; "base match
rate formula" for the formula and the scoring). **Decisions:** D82–D88. **Deliveries:** 1, 2 (with 2b), 3, 4 R2.
**Result:** `python run.py match` scores every food in USDA FoodData Central that carries an amino acid panel against
three references — the skeletal muscle standard with its non-protein metabolite pools (the primary, D82), the
protein-only standard, and Gorissen et al. 2018's measured human muscle column (the original NutriMatch reference) —
and writes the score, the limiting amino acid, and the spreadsheet's own arithmetic per food so that any number can be
rebuilt by hand. Version label `Match Rate 1.0`, reference `skeletal muscle protein + non-protein metabolite pools
(carnosine), fiber-type mix`.

## What the step set out to do

Formalise the Match Rate calculation in code against the new standard, run it against the USDA food database, and
publish version-labelled results (plan R7 §5). Three things had to exist first: a food table with citable provenance,
a reconciliation of the calculated standard against a measured one so a reader knows how far apart they are, and the
formula itself — which lived in a spreadsheet and had never been written down.

## Procedure

1. **The reference tier named** (D82, delivery 1). The Step 4b tier was renamed `non_protein_metabolite_pools`
   after a literature check showed "bound" already means *protein*-bound in this field. The name is carried into
   every file, column, stage and label so that it survives being read cold.
2. **The food side** (D83, delivery 1). `pipeline/usda.py` reads two hand-downloaded FoodData Central archives —
   SR Legacy (April 2018) and Foundation Foods (April 2026) — as published, rules U1–U7: the archive is identified by
   reading its own `food.csv`, not its file name; superseded versions are dropped by the archive's own membership
   file; the twenty amino acids are joined to USDA's nutrient names through the IUPAC trivial names so that no amino
   acid is named in code or config; completeness is stated per food, never decided; nothing is deduplicated across
   archives. **5,156 foods carry an amino acid panel** — 5,102 of SR Legacy's 7,793 and 54 of Foundation's 395.
3. **The Gorissen comparison** (D84, D85, deliveries 2 and 2b). `pipeline/comparison.py` puts the calculated standard
   beside Gorissen et al. 2018's human vastus lateralis column, both renormalised over the fourteen rows the paper's
   method can measure (W, D, N excluded; P and C excluded because the paper prints 0.0 for them and a printed zero
   from a run that reports both in fifteen other sources is a missing cell, not a measurement; E and Q summed because
   hydrolysis converts glutamine). Rules X1–X6; no fit, no scaling, no verdict, no path back into the standard.
   `docs/gorissen_comparison.md` is written from the generated tables. The transcription is checked against the
   paper's own printed sums every run, with the paper's own statement of which rows it sums as essential read from
   the config rather than named in code.
4. **The formula transcribed** (D86, D87, delivery 3). The author's spreadsheet `Base_Match_Rate_Formula.xlsx` was
   read cell by cell before anything was discussed (the D50 pattern). What it computes: scale the food to the
   reference's essential-amino-acid total; per amino acid, food ÷ reference; take the minimum; divide everything by
   that minimum so the limiting amino acid meets the reference; sum the "need to consume"; call the excess above the
   reference total "wasted"; Match Rate = 1 − wasted. Written out algebraically, every step after the minimum cancels
   — wasted is exactly 1 − k — so the score is the smallest of (food share ÷ reference share) over the scored set,
   and protein content never enters. The four foods in the spreadsheet are the transcription test: `match.py` must
   reproduce their four scores to the last digit, and a literal re-run of the spreadsheet's Steps 3–10b must equal
   the minimum ratio, before any other test counts. The scored set is FAO 2013's nine headings with each group
   reduced to its indispensable member (methionine, phenylalanine), the author's decision (D86); the original
   reference is scored on the eight amino acids the paper sums (no tryptophan — acid hydrolysis destroys it).
5. **Two scripts** (D88). `pipeline/match.py` is the public single-food scorer. The blend and fortification script
   — a fixed blend of sources plus free amino acids; the inverse problem of which additions lift the limiting ratio
   to a target; reverse-engineering a nutrition facts panel from its ingredient list — is a separate script that
   imports `match_rate()` and is not part of the public repository at publication.
6. **The output layout** (delivery 4 R2), after the author's review of delivery 3's tables: every number must be
   rebuildable in Excel on the author's own spreadsheet from the row itself. So `match_rate_steps.tsv` carries the
   spreadsheet's walk per food and reference in the spreadsheet's order and names, each reference's values in the
   header (the spreadsheet's column A); `match_rate_per_food.tsv` carries the score, limiting amino acid, ΣEAA, ΣEAA
   as a percent of protein and the Step 7 ratio per amino acid; `match_rate_by_reference.tsv` puts the three scores
   side by side per food. Headers were cut to provenance only. Delivery 4 R1 had also shortened the reference ids
   and replaced two tables; the author rejected that as beyond what was asked, and R2 reverted every name.
7. **Runs** 2026-09-17 and 2026-09-18 on the workstation: every test passed; excerpts uploaded and read; twenty
   steps-table rows re-derived by hand from their own grams and the header's reference values, all agreeing to four
   decimals.

## Results

All figures from `outputs/match/` and `outputs/comparison/` (run of 2026-09-18T00:26Z).

| | primary (with pools) | protein only | Gorissen 2018 (original) |
|---|---|---|---|
| scored set | nine (H I L K M F T W V) | nine | eight (no W) |
| foods scored / not scored | 4,952 / 204 | 4,952 / 204 | 5,019 / 137 |
| median Match Rate | 69.86 % | 69.32 % | 74.74 % |
| foods at 0 % (a published zero in the scored set) | 20 | 20 | 6 |

**Where the limiting amino acid falls** (`limiting_amino_acid_counts.tsv`, percent of scored foods; ties count for each):

| amino acid | primary | protein only | Gorissen 2018 |
|---|---:|---:|---:|
| tryptophan | **49.7** | 49.7 | not scored |
| lysine | 27.0 | 27.0 | 23.2 |
| methionine | 15.0 | 15.1 | 9.9 |
| threonine | 3.8 | 3.9 | 0.9 |
| isoleucine | 2.6 | 2.6 | 0.6 |
| valine | 1.2 | 1.4 | 2.4 |
| histidine | 0.6 | 0.2 | 23.2 |
| leucine | 0.3 | 0.3 | 0.9 |
| phenylalanine | 0.1 | 0.1 | **39.1** |

Under the original reference, phenylalanine and histidine between them limit 62 % of foods; under the new reference
they limit under 1 %. That is the Gorissen comparison seen from the food side: the measured column asks for 6.25 %
phenylalanine and 4.61 % histidine where the calculated standard carries 4.83 and 3.07, and every food that was short
of those two is no longer short. Fish shows it most plainly — pollock scores 59.7 % (histidine-limited) against the
old reference and 89.2 % against the new one. The pools column moves scores by under a point except where histidine
limits (pollock 91.7 → 89.2 with the pools, which add histidine to the reference), as it should.

**The tryptophan finding.** Under the new reference, tryptophan is the limiting amino acid for half of all scored foods.
Half of all scores are therefore set by the food's tryptophan against the standard's — 1.541 % of amino acid mass in
the free convention, 1.98 g per 100 g protein — and by nothing else. Tryptophan is also the one scored amino acid the
Gorissen comparison cannot check, because acid hydrolysis destroys it and the paper does not report it; the comparison
excludes it from both sides. So the number that decides the most scores is the one with the least independent
measurement behind it — and it is exact at Layer A, because it comes from sequence, which is the argument for a
sequence-derived standard in one sentence. USDA's meats sit near 1.0–1.1 g tryptophan per 100 g protein; the standard
asks for more. Three candidate explanations, none separable with what is on disk: the sequence-derived value is right
and measured muscle tryptophan is low; USDA's tryptophan (largely older or calculated values) is low; Layer B weights
tryptophan-rich proteins too heavily. **Documented, not acted on** (author, 2026-09-18): a tryptophan-preserving muscle
measurement transcribed into the comparison stage would give tryptophan the treatment histidine got.

**The ranked tables** (`match_rate_ranked_*.tsv`). Under the new reference the top of the list is rabbit — five
preparations of it at 92.83–92.92 %, raw, wild, roasted and stewed within a tenth of a point of each other — then
scallop, sockeye salmon, abalone, duck. Under the original reference the top ten are led by fried chicken, meat sticks,
beef lungs, sardines and two entries of canned peaches at 0.5 and 0.9 g protein. Two things the ranked tables show
that the score does not decide: **the score is a proportion, so a food with almost no protein can score high** (peaches
86.2 % new / 90.5 % old on 0.91 g protein; almond milk 11 %; cantaloupe 39 %); and **USDA lists many preparations of
one food whose panels differ by a rounding error** (sunflower kernels oil-roasted with salt, without salt, and dried:
54.5717, 54.5717, 54.5707 %; the five rabbits). Both are display questions; both are documented and deferred.

**The transcription test.** The spreadsheet's four foods against its column A: rice 44.99 % (lysine), chia 67.38 %
(lysine), pea 51.31 % (methionine), chicken breast 82.38 % (phenylalanine) — cells G156:J156 to the last digit, and
the spreadsheet's Steps 3–10b re-run literally give the same four numbers (`tests/test_match.py`).

## Decisions and rules

D82 (the tier name), D83 (USDA as the food side, rules U1–U7), D84 (the Gorissen comparison as a stage, rules X1–X6),
D85 (the printed zeros excluded, the mass balance computed), D86 (the scored set), D87 (the formula as transcribed,
rules M1–M6), D88 (two scripts, split at publication). Standing rules added to `docs/conventions.md`: generated-file
headers are provenance only; every reported number must be rebuildable by hand from its row.

## Method lessons (all pinned in the handoff)

- **Transcribe first, litigate after** held again. The algebra showing the surplus arithmetic collapses to the
  limiting ratio was written before the spreadsheet arrived and confirmed cell by cell in it; the author's framing
  and the minimum ratio are the same quantity stated twice.
- **A mathematically equivalent shortcut is not the deliverable.** Delivery 3 wrote the Step 7 ratios and dropped the
  Step 10b frame, the grams, and the totals the author reads the sheet in; it was correct and useless for validation.
  The steps table is the fix, and the rule is now standing: the row carries the inputs and the sheet's own vocabulary.
- **Do what was asked and nothing beside it.** Delivery 4 R1 shortened the D82 reference ids and replaced two tables
  while executing a list of column removals. Renames and removals are proposed and approved before they are built.
- **Every changed package gets a new revision label.** A corrected package was twice re-issued under the same file
  name on 2026-09-17; the author ran the wrong one. R2 is the correction and the rule is restated.
- **Headers are provenance, not explanation.** Eighteen header lines of rule prose were in the way of reading the
  table; the hashes and the reference lines are what a reader needs, and the rest lives in `docs/`.
- **The mass balance is a bound, not a test.** Uncorrected hydrolysis losses mean recovery below 100 % is expected
  under any protein figure; the zeros decision stands on the proline argument alone. §3 and §5 of the comparison
  document and D85's rationale were reworded to report mode.
- **For the refactor (Step 7b):** the code names no amino acid, and where it briefly did (the essential-row grouping
  for X4) the fix was a config key with the paper's location; the pattern holds.

## Open (documented; the author's call on 2026-09-18 is to do nothing yet)

Tryptophan cross-check source; a minimum protein content for the ranked tables (a floor in `config/match_rate.ini`
applied to the ranked view only, the per-food tables complete — the author's lean); consolidation of USDA
preparations whose panels differ within rounding (rabbit, sunflower kernels); which archive wins when a food is in
both; cystine → cysteine; the reference-foods table (pea isolate and the current blend need rows in
`config/food_amino_acids_other_sources.csv`); the leucine correction factor, with the blend script; Mingrone 2001's
protein assay; the "v1.0" label against the no-versioning-before-release rule. All in
`docs/handoffs/06_match_rate_handoff.md`, "Open".
