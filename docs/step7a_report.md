# Step 7a report — Reorganization

**Thread:** 2026-09-18. **Decisions:** D90–D94. **Deliveries:** 1 (R3), 2 (R3), 3 (R4), 4 (+ a
two-file fix), 5.
**Result:** `outputs/` is laid out by what a file is; every generated path resolves at call time;
a stage hashes only what it reads; the literature sources have their own config file and a
category map; `fetch-literature` is its own command and never downloads. No number changed.

## What the step set out to do

Put `outputs/` into a shape a reader can navigate, close the bookkeeping pinned in Steps 4b and
6, and change nothing else — before blood is built, so blood is built into its own folders from
its first file (D89). The abstraction stays after blood (Step 7c, D71).

## Procedure

1. **The placement table first** (delivery 1). The author posted the file list of `outputs/`;
   every one of 111 files was assigned a new location and approved before any code changed.
   Two amendments came out of that review: the protein-set tables go in
   `intermediate/protein_set/` rather than loose in `intermediate/`, and the ranked Match Rate
   tables keep their full names inside their reference folders — "if I'm looking at the
   directory, I have no idea" is the same argument that decided both.
2. **The move** (delivery 1). `outputs_original/` was moved outside the repository, `outputs/`
   emptied, and the whole pipeline rerun so the code wrote to the new paths.
   `compare_output_trees.py` — tooling for this step, deleted at its close — matched files by
   name across the two trees and compared bodies below the header.
3. **The literature split** (deliveries 2 and 3). Twenty-two `[source.*]` sections moved into
   `config/literature_sources.ini` byte-unchanged apart from one added `used_for` line each;
   `fetch_literature.py` was rewritten as a hash-only stage; `manifest_map.tsv` was added.
4. **The bookkeeping** (delivery 4). Paths resolved at call time, per-stage hashing, per-command
   excerpt specs, sorting, and the removal of a code default.
5. **The docs** (delivery 5). `run_order.md` written, `pipeline_map.md` and `conventions.md`
   brought to the tree as built, the README trimmed.

## Results

**The acceptance check.** 111 files before, 105 after. 25 identical with their path unchanged,
57 identical after moving, 18 binary (matplotlib does not write byte-reproducible PNGs, so
plots are checked for presence), 3 differing only in lines the check can explain, 2 differing
for a reason outside this step, 6 not regenerated.

The three explained files differ in a run timestamp or in a recorded hash of a file whose own
header names a moved path — the cascade the handoff predicted, one hop deep. `standard_summary.ini`
is the clearest case: 217 body lines, 1 timestamp, 16 path lines, 2 hash lines, **every amino
acid number identical**. The hashes that did *not* move are the better evidence:
`amino_acid_composition_per_protein.tsv`, `isoform_deltas.tsv`, `processing_mass_deltas.tsv`
and `ptm_mass_deltas.tsv` are byte-identical across the move, because those stages read only
`data/` and `config/`, nothing moved, and their headers are unchanged.

The two unexplained differences were real and not ours: `gorissen_2018` had been added to the
literature manifest after the old outputs were generated, so the inventory gained a row and its
summary moved 20 → 21 files and 11 → 12 documents-not-opened. **The check was left failing.** A
check widened every time it fires stops being worth running, and this one had just isolated a
genuine difference to two files and two lines while explaining the other seven.

Four files are not regenerated, each for a stated reason: `match_rate.tsv` (superseded in Step
6), `match_rate_filtered.xlsx` (the author's own working file), `pool_overlap.tsv` (written by
`enumerate_pool`, a network stage, skipped by `--offline`), and `profile_fiber_types_combined.{png,svg}`
(`plots.py` no longer draws it). `flags.tsv` is absent because the run raised zero flags and
`outputs/` started empty; the old copy was header-only, so absence and presence say the same
thing.

**The final state.** 213 tests pass with no skips. The tree:

```
outputs/  flags.tsv
          intermediate/  protein_set 11, composition 7, digest 6,
                         literature_inventory 3, mass_fractions 17
          standard/      13 at the top, plots 12,
                         uncertainty 5 + plots 6, stress 5, sensitivity 4 + plots 2
          match/         6 at the top, one folder per reference
          usda/ 4        comparison/ 4
```

**The literature map**, generated from `used_for`, sorted A–Z to read beside the folder listing:
22 sources across four columns — 10 skeletal muscle, 6 non-protein metabolites, 2 original,
1 (`fao_2013`) in three at once, 4 marked `unused`. 14 sources declare files (21 files), 8 are
cited but not stored.

## Decisions and rules

D90 (the layout, placement as data), D91 (paths resolved at call time), D92 (a stage hashes only
what it reads; no timestamps in a hash; A0 amended), D93 (`literature_sources.ini`, `used_for`,
the category map), D94 (fetch retired; `fetch-literature` its own command; F3/F4/F4b/F6 retired,
F8 and F9 added).

## Method lessons

- **Three bugs in this step were the same bug**: a value that looks like an input but is not.
  A path frozen at import cannot follow a redirected root — twice — and a retrieval timestamp
  inside a hash makes an unchanged file look changed. Each surfaced far from its cause: a
  `relative_to` ValueError, a missing config section, six outputs reading stale after a no-op
  re-hash. `tests/test_paths.py` and `VOLATILE_KEYS` exist so the fourth one fails immediately.
- **A delivery contains only the files it changed.** A config file was shipped in its pre-run
  state with a note in the delivery document saying not to overwrite the author's copy, which
  held three hashes his own run had recorded. Documenting a hazard is not removing it. Where a
  config file must ship, it is checked against the author's current state plus exactly the
  intended edits before packaging.
- **Do not widen a check because it fired.** See the two remaining differences above.
- **The acceptance check earned its cost before it was ever run in anger** — writing it forced
  the placement table to be explicit, and the first real run of the reorganized code caught a
  defect that would otherwise have appeared months later as a blood bug.
- **A name is read where it is found, not where it was written.** `scripts/` was rejected
  ("every python file is a script — it describes nothing"), the map is sorted A–Z because it is
  read beside the folder listing, and the ranked tables keep their long names because a file
  reviewed alone must say what it is. The same argument produced three different fixes.
- **Ask about the thing that is actually open.** A question about whether to create a config
  file was read three times as a request for the file. A decision question and a file request
  need to look different on the page.

## Open (documented, not acted on)

- **`pipeline/comparison.py` has no `--check` currency test.** Every other stage has one, which
  is how `match` reported itself stale after the standard was rewritten; `comparison` reads the
  same two tables and would have sat there quietly. It should have one.
- **`pipeline/excerpt.py` lists three files the pipeline never writes** —
  `sensitivity_non_protein_metabolite_pool_{turnover_frame,basis,sex}.tsv`, planned in Step 4b
  and stated rather than computed (no turnover source, no whole-muscle carnosine value). The
  stage warns about all three on every run. Either compute them or remove the entries.
- **Should a zero-flag run write a header-only `flags.tsv`?** Today "no flags raised" and "the
  stage never ran" are indistinguishable on disk.
- Everything carried from Step 6, untouched here: the tryptophan cross-check source, the protein
  floor for the ranked tables, consolidation of near-identical USDA preparations, archive
  precedence, cystine → cysteine, the reference-foods table, the leucine correction factor,
  Mingrone 2001's protein assay, the "v1.0" label.
