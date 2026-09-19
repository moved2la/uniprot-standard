# Run order

One screen. What to run, in what order, and how far back to start.

## The run

The commands in dependency order. How much of this list you run depends on why you are running
it — see "Where to start", below.

```
0. apply the files from any delivery
1. python run.py fetch-literature
2. python run.py protein-set --offline
3. python run.py composition --offline
4. python run.py mass-fractions
5. python run.py standard
6. python run.py blood-protein-set --offline
7. python run.py blood-composition
8. python run.py blood-mass-fractions
9. python run.py blood-standard
10. python run.py composite
11. python run.py usda
12. python run.py comparison
13. python run.py match
14. python run.py test
```

- If a command stops, fix it and rerun **that same command** until it is green, then continue.
- **Tests are step 14, not step 1.** A test confirms a state that has already been established; it
  cannot establish one. Run against a repo in an unknown state, it reports failures whose cause is
  ambiguous — was it the last change, or something three days old? Run against a tree that was
  never built, it passes and means nothing. Get it working, then test that it stays working.
- Drop `--offline` on 2, 3 and 6 only when a sequence genuinely needs refetching — about two
  hours for muscle, about 45 minutes for blood.
- **The order is dependency order (D123):** blood's sequence fetch is additive on top of muscle's
  store (S1), so blood follows `standard`; `composite` mixes the muscle and blood standards, so it
  follows `blood-standard`; `match` reads the composite as a reference, so it comes last. The one
  declaration of this list is `common.COMMAND_ORDER`; `tests/test_run_order.py` holds these assertions.
- **A blood-only change** (a file under `config/blood/`, `pipeline/mass_fractions_pools.py`,
  `pipeline/aggregate_pools.py`) does not require 1–5: it reruns from the blood command that reads
  what changed, then 10–14. A change to shared code (`common.py`, `run.py`, a stage muscle also
  runs) is the full run from 1.
- A new command (a new category's stage) slots into this list in dependency order. Nothing is
  appended to the end because it is new.

## Where to start

**Start at 1 — the whole list — when:**

- closing out a thread
- deploying a build for the first time: a full delivery applied, a fresh clone, a repo you have
  not run today
- the change touches shared code (`pipeline/common.py`, `run.py`, a stage that more than one
  command runs)
- the state of the tree is unknown, for any reason

**Start further down for debugging and minor adjustments.** Find the earliest command that reads
anything you changed, start there, and **carry through to the end of the list.**

- The shortcut is at the front, never at the back. Starting at the right place is the saving;
  stopping early is not. Every command after the one you started at either reads what that command
  wrote or records its hashes, so leaving them unrun leaves the tree stale in a way only the next
  full run will find.
- The question to answer is "what is the earliest command that reads anything I changed?" —
  the table under "What to rerun after an edit" answers it for the usual cases.
- Carrying to the end can be a single command. A patch to `pipeline/match.py` alone starts at 13,
  and 13 is the end of the list, so `python run.py match` is the whole run (2026-09-19).
- `python run.py test` is not needed as a separate step when you carry to 13: `match` is last in
  `COMMAND_ORDER`, so the test run at the end of it skips nothing and is the full suite.
- A repo you have not run today is not a minor adjustment. If in doubt about the state of the
  tree, that is the unknown-state case: start at 1.

Earned in Step 7b, 2026-09-19: a thread opened with `run.py test` on a repo in an unknown state
and then jumped to command 4 on the strength of what the test reported. Both failures were real,
but half of them were a day-old staleness nobody had run into yet, and the sequence produced two
days of "you need to do this first". That is the unknown-state case, and it is what starting at 1
prevents — not a reason to rebuild the whole tree for a one-file patch.

## The commands, in order

| # | Command | Reads | Writes |
|---|---|---|---|
| 1 | `python run.py fetch-literature` | `config/literature_sources.ini`, the files in `data/literature/<source>/` | `data/literature/manifest.ini`, `manifest_map.tsv`, `README.md` |
| 2 | `python run.py protein-set` | `config/protein_set_decisions.ini`, Gene Ontology, UniProt | `config/accessions.ini`, `segments.ini`, `data/uniprot_sequences.ini`, `outputs/intermediate/protein_set/` |
| 3 | `python run.py composition` | `config/composition_decisions.ini`, the sequence store, PubChem, UniProt PTM list | `outputs/intermediate/composition/` |
| 4 | `python run.py mass-fractions` | `config/mass_fraction_decisions.ini`, `literature_sources.ini`, the manifest, the literature files | `config/mass_fractions/`, `mass_fractions_per_entry.tsv`, `outputs/intermediate/{digest,literature_inventory,mass_fractions}/` |
| 5 | `python run.py standard` | the weights, the composition table, `fao_2013_…ini`, `fiber_type_mix.ini`, `non_protein_metabolite_pools.ini`, `uncertainty_settings.ini`, `stress_test_settings.ini` | `outputs/standard/` and its `uncertainty/`, `stress/`, `sensitivity/`, `plots/` |
| 6 | `python run.py blood-protein-set` | `config/blood/protein_set_decisions.ini`, the column maps in `config/blood/mass_fraction_decisions.ini`, the Bryk and Geyer files, UniProt | `data/blood/`, `config/blood/accessions.ini`, `segments.ini`, the sequence store (additive), `outputs/blood/intermediate/protein_set/` |
| 7 | `python run.py blood-composition` | the blood accessions and segments, the sequence store, the shared masses and PTM list | `outputs/blood/intermediate/composition/` |
| 8 | `python run.py blood-mass-fractions` | `config/blood/mass_fraction_decisions.ini` (`[columns.*]`, `[cross_checks]`, `[immunoglobulin_bound]`), the manifest, the literature files | `config/blood/mass_fractions_per_entry.tsv`, `outputs/blood/intermediate/mass_fractions/` |
| 9 | `python run.py blood-standard` | the blood weights and composition table, `config/blood/compartment_mix.ini`, `fao_2013_…ini`, `uncertainty_settings.ini` | `outputs/blood/standard/` and its `plots/` |
| 10 | `python run.py composite` | `config/tissue_mass_fractions.ini`, the muscle standard (with pools) and residue tables, the blood standard, `compartment_split.tsv`, `sensitivity_dominant_protein.tsv`, `non_protein_metabolite_pools.ini`, `fiber_type_mix.ini` | `outputs/composite/` |
| 11 | `python run.py usda` | `config/usda_food_data.ini`, the archives in `data/usda/` | `outputs/usda/` |
| 12 | `python run.py comparison` | `gorissen_2018_comparison.ini`, the standard tables, the manifest | `outputs/comparison/` |
| 13 | `python run.py match` | `config/match_rate.ini`, the standard tables, the composite, `outputs/usda/` | `outputs/match/` and `match/<reference>/` |

`python run.py test` runs the tests alone. `python run.py excerpt` is tooling, not a stage.

`--offline` skips the stages that touch the network. Only `protein-set`, `composition` and
`blood-protein-set` have any: nothing else in the pipeline downloads (D80).

## What to rerun after an edit

**The rule: rerun the command that reads what you changed, then everything after it in the table
above.** "Rerun from" means from that number to the end of the list, not that number alone.

| You changed | Rerun from |
|---|---|
| a literature file added, replaced, or re-obtained | 1, then 4 |
| `literature_sources.ini` (a source, a `used_for`, a hash) | 1, then 4 |
| `protein_set_decisions.ini` | 2 |
| `composition_decisions.ini` | 3 |
| `mass_fraction_decisions.ini` (column map, digest rule) | 4 |
| `fao_2013_…`, `fiber_type_mix`, `non_protein_metabolite_pools`, `uncertainty_settings`, `stress_test_settings` | 5 (`fao_2013_…` and `uncertainty_settings` are read by 9 and 10 too, so through 13) |
| `config/blood/protein_set_decisions.ini` or a column map in `config/blood/mass_fraction_decisions.ini` | 6 |
| `config/blood/mass_fraction_decisions.ini` (`[cross_checks]`, `[immunoglobulin_bound]`) | 8 |
| `config/blood/compartment_mix.ini` | 9 |
| `tissue_mass_fractions.ini` | 10 |
| `usda_food_data.ini` or a USDA archive | 11, then 13 |
| `gorissen_2018_comparison.ini` | 12 |
| `match_rate.ini` | 13 |
| a stage's code | that stage's command, then everything after it |
| shared code (`pipeline/common.py`, `run.py`, a stage more than one command runs) | 1 |
| a test file only | the command the test belongs to, for its test run |
| documents only (`docs/`, `README.md`, `PROVENANCE.md`) | nothing |

**The chain that catches people: 5 → 10 → 12 → 13, and 9 → 10 → 13.** `composite`,
`comparison` and `match` all record the hashes of the standard tables in their headers, so any
rerun of `standard` or `blood-standard` makes them stale. Running a standard alone and stopping
leaves `match` claiming inputs it no longer has.

Each command runs the tests at the end, and a currency test per stage checks that what is on
disk equals a fresh build. That is what tells you a rerun is owed — except for `comparison`, which
has a `--check` that no test runs, and `literature_inventory`, which has neither. Those two go
stale silently (see `docs/step7a_report.md`, open items).

## What does not force a rerun

- Re-running `fetch-literature` on unchanged files. A stage hashes only the manifest entries it
  reads, and retrieval times are excluded from that hash, so re-hashing a file that has not
  moved changes nothing downstream.
- A source added to `literature_sources.ini` for another category. A stage reads only the sources
  whose `used_for` names its own category (D106), and hashes only the manifest entries it reads,
  so another category's source changes nothing for it. **This held only after D106:** before it,
  blood's two `role = primary` sources stopped the muscle mass-fraction stage, which scanned the
  whole file and found three primaries.

## A full rebuild

"The run" at the top of this file, from command 1. There is no other kind of rebuild. It is what a
thread close, a first deployment, a shared-code change and an unknown tree all get; a known tree
with a known change gets the starting point "Where to start" and the rerun table give it, carried
through to the end.
