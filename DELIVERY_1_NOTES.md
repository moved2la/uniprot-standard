# Step 6 — delivery 1 (revised)

**Date:** 2026-09-17 · **Decisions:** D82 (the tier's name), D83 (the USDA archives)
**Contents:** 20 files — 16 changed, 4 new. Extract over the repo root.
**Supersedes** the first cut of this delivery, which used `_with_metabolites`.

---

## 1. What to do, in order

1. **Extract** the tarball over the repo root.
2. **Delete three files by hand** — their replacements ship under new names, and code never deletes (D61):
   - `config/bound_metabolite_pools.ini`
   - `pipeline/bound_pools.py`
   - `tests/test_bound_pools.py`

   (`__pycache__` can go too if the old `.pyc` lingers.)
3. **Download both USDA archives** from https://fdc.nal.usda.gov/download-datasets — the **CSV** links, not JSON:
   Foundation Foods April 2026 (3.7 MB) and SR Legacy April 2018 (6.7 MB). Put both **zips** in `data/usda/`,
   **unrenamed**, not unzipped. They are gitignored; `data/usda/manifest.ini` is the committed record.
4. `python run.py usda`
5. `python run.py standard`
6. **Delete nine superseded files** from `outputs/standard/` (list in §2).
7. `python run.py excerpt` and upload — I have not seen the real adjusted table or any SR Legacy output.

If a zip is missing, the stage stops and names the data type. Deliberate: a table built from one archive would look
like the whole database.

---

## 2. The tier's name (D82)

**The tier is `non_protein_metabolite_pools`.** Long on purpose — it has to survive being read cold by a reviewer.

### Why not "bound pools"

"Bound" is already taken. In amino acid nutrition, free amino acids are *defined* as those not bound in protein, and
papers reporting "free and bound amino acids" in a tissue mean free versus in-protein. `with_bound_pools` would point
a reviewer at the protein-only state — the other file. That is a collision no definition repairs, because the reader
arrives with the prior.

### Why not just "metabolites"

Free amino acids are metabolites too, and they are exactly what D68 excludes. Under-specified rather than
misdirecting, but it leans on the definition to do all the work.

### Why this is coined at all

There is no established term for this class. The literature has **"non-protein nitrogen-containing compound"**
(how carnosine was described on discovery — but it also covers free amino acids and urea) and
**"histidine-containing dipeptides"** (carnosine's family only, no creatine or glutathione). Nobody groups exactly
what we group, so nobody named it. We coin, and we define.

### The definition, now in `docs/methods.md`

Three states, stated together so the exclusion is visible:

1. **Protein-bound** — in peptide chains. The protein-only standard.
2. **Free** — single amino acids in solution, turned over continuously. Excluded by design (D68).
3. **Non-protein metabolite pools** — held in stable non-protein compounds: carnosine (histidine), later creatine
   (glycine) and glutathione. Replaced over weeks. Eaten to be there, so a real demand a protein-only count misses.

A pointer sits in the Naming section of `docs/conventions.md`. Both note that **D68, D73 and D75 predate the term and
say "bound pool"** — they mean state 3. Logged decisions are not rewritten; D82 records the change.

### What carries the name

| | |
|---|---|
| table | `_calculated_amino_acid_standard_with_non_protein_metabolite_pools.tsv` |
| the four computed columns | `total_I_with_non_protein_metabolite_pools`, `total_IIa_…`, `total_IIx_…`, `standard_with_non_protein_metabolite_pools` |
| companions | `non_protein_metabolite_pool_adjustment_per_amino_acid.tsv`, `non_protein_metabolite_pool_amounts_per_kg_muscle.tsv`, `sensitivity_non_protein_metabolite_pool_*.tsv`, `eaa_subset_with_non_protein_metabolite_pools.tsv`, `non_protein_metabolite_pools_summary.ini` |
| config | `config/non_protein_metabolite_pools.ini` |
| stage / module | `non_protein_metabolite_pools` / `pipeline/non_protein_metabolite_pools.py` |
| rule | A10 "Non-protein metabolite pools" |
| source roles | `role = non_protein_metabolite_pool_single_fiber`, `…_turnover`, `…_scope_citation` |
| version label | `ref = skeletal muscle protein + non-protein metabolite pools (carnosine)` |

The six contractile / builders columns are copied from the protein-only standard and keep their names, as agreed.
The protein-only standard is untouched (D76).

**On the scope.** You asked for four columns; this changes far more than four. The reason it is right *now* and was
not right before: with the term settled, leaving `bound_pool_*` file names beside
`_with_non_protein_metabolite_pools` columns would leave both coined terms live in the repo, which is the confusion
we just spent four messages removing. One term, everywhere. If you would still rather cut it back to the four column
headers only, say so and it is a ten-minute patch.

### The nine files to delete after step 5

```
_calculated_amino_acid_standard_with_bound_pools.tsv
bound_pools_summary.ini
bound_pool_adjustment_per_amino_acid.tsv
bound_pool_amounts_per_kg_muscle.tsv
sensitivity_bound_pool_turnover_frame.tsv
sensitivity_bound_pool_basis.tsv
sensitivity_bound_pool_sex.tsv
sensitivity_bound_pool_spread.tsv
eaa_subset_with_bound_pools.tsv
```

Same numbers under the old names; nothing reads them.

---

## 3. The food side (D83)

**New:** `pipeline/usda.py`, `config/usda_food_data.ini`, `tests/test_usda.py`, `tests/test_usda_outputs_current.py`.
**Command:** `python run.py usda` (offline; nothing downloaded by code — D80).

FoodData Central publishes one long table: a row per food per nutrient, per 100 g of food, every nutrient class in it,
separated only by a nutrient id. The stage pivots that into one row per food and does nothing else — no conversion,
no dropped food, no preferred archive.

### Outputs (`outputs/usda/`)

| File | What it is |
|---|---|
| `amino_acids_per_food.tsv` | one row per food carrying at least one amino acid: id, data type, release, description, category, protein and nitrogen, `<letter>_g_per_100g` for the twenty, cystine and hydroxyproline, the count present, the letters absent, the smallest sample count |
| `foods_without_amino_acids.tsv` | the archive's own foods that carry none, so what is *not* in the table above is on the record |
| `usda_nutrient_map.tsv` | the audit trail of the name join: every nutrient looked for, the id and unit the archive gives it, what it was carried as |
| `usda_archive_summary.ini` | per archive: file, hash, size, membership rule, counts; per-letter coverage; any flags |
| `data/usda/manifest.ini` | name, SHA-256, size, data type, release, first seen — **committed**; the zips are not |

### The rules (U1–U7, in `docs/conventions.md`)

- **U1** A zip is identified by reading its own `food.csv`, never by its file name. Its SHA-256 is its identity.
- **U2** A food counts when the archive's *own* membership file lists it — what drops the superseded versions. In the
  April 2026 Foundation archive `food.csv` has 469 Foundation rows and `foundation_food.csv` lists 395; the other 74
  are older versions of the same foods.
- **U3** A nutrient is amino acid X when its name equals the IUPAC-IUBMB trivial name of X, from the table the
  composition step already parsed. **No amino acid is named in code or config**, and a test fails if one appears in
  the config file.
- **U4** The archive's `amount` as published, with min/max/median/sample count where given. A nutrient in an
  unexpected unit is flagged, not converted.
- **U5** Completeness is stated, not decided.
- **U6** No cross-archive dedup.
- **U7** Nothing in `data/` or `outputs/` is deleted, renamed or moved.

### What the sandbox run found (Foundation, April 2026)

395 foods in membership · 54 with amino acids · 341 without · archive hash `70457ee9…c655`.

The join matched all twenty letters, **asparagine** (id 1231) and **glutamine** (1233) included: those ids exist in
FoodData Central's nutrient list but no food populates them, because acid hydrolysis converts both. So every row's
`amino_acids_absent` contains at least `NQ`. SR Legacy is where the thousands of foods are; not run here.

---

## 4. Open, for the formula discussion

1. **The formula itself.** `match.py` waits for the spreadsheet.
2. **The scored set** — the nine FAO headings with each group reduced to its indispensable member (methionine alone,
   phenylalanine alone; cysteine and tyrosine unscored). Your reduction, not FAO's, so it needs a D-number.
3. **Which archive wins** when a food is in both Foundation and SR Legacy.
4. **A minimum protein content** below which a food is not scored — 18 of the 54 Foundation panels carry a zero
   amino acid value (butter, cantaloupe). A zero in the scored set scores the food 0 %.
5. **Cystine → cysteine.** Some foods report one, some the other, never both. Whether USDA's cystine figure is
   already cysteine-equivalent has to be read from the documentation, not assumed. Nothing in the indispensable set
   depends on it.

---

## 5. Tests

14 synthetic tests for the USDA stage, 4 currency tests, and the 12 existing tests of the renamed stage — all green
here. Your run is the one that counts; this sandbox has no network and no pytest, so I ran them through a minimal
stand-in.
