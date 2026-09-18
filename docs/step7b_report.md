# Step 7b report — Blood (part 1: scope, sources, decisions)

**Thread:** 2026-09-18. **Decisions:** D95–D105. **Deliveries:** 1 (R2), 2 (documents).
**Result of part 1:** blood is scoped as plasma + erythrocytes; ten sources are on disk, hashed
locally, `fetch-literature` green (34 files, 29 sources × 5 columns); the plasma : erythrocyte
split is specified from measured reference values; the two findings that will shape the build are
on record; F7 is amended and the source acquisition protocol is written down.
**Part 2 (the build) appends to this report from its own thread.**

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
