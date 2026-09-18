# Run order

One screen. What to run, in what order, and what to rerun after an edit.

## The commands, in order

| # | Command | Reads | Writes |
|---|---|---|---|
| 1 | `python run.py fetch-literature` | `config/literature_sources.ini`, the files in `data/literature/<source>/` | `data/literature/manifest.ini`, `manifest_map.tsv`, `README.md` |
| 2 | `python run.py protein-set` | `config/protein_set_decisions.ini`, Gene Ontology, UniProt | `config/accessions.ini`, `segments.ini`, `data/uniprot_sequences.ini`, `outputs/intermediate/protein_set/` |
| 3 | `python run.py composition` | `config/composition_decisions.ini`, the sequence store, PubChem, UniProt PTM list | `outputs/intermediate/composition/` |
| 4 | `python run.py mass-fractions` | `config/mass_fraction_decisions.ini`, `literature_sources.ini`, the manifest, the literature files | `config/mass_fractions/`, `mass_fractions_per_entry.tsv`, `outputs/intermediate/{digest,literature_inventory,mass_fractions}/` |
| 5 | `python run.py standard` | the weights, the composition table, `fao_2013_…ini`, `fiber_type_mix.ini`, `non_protein_metabolite_pools.ini`, `uncertainty_settings.ini`, `stress_test_settings.ini` | `outputs/standard/` and its `uncertainty/`, `stress/`, `sensitivity/`, `plots/` |
| 6 | `python run.py usda` | `config/usda_food_data.ini`, the archives in `data/usda/` | `outputs/usda/` |
| 7 | `python run.py comparison` | `gorissen_2018_comparison.ini`, the standard tables, the manifest | `outputs/comparison/` |
| 8 | `python run.py match` | `config/match_rate.ini`, the standard tables, `outputs/usda/` | `outputs/match/` and `match/<reference>/` |

`python run.py test` runs the tests alone. `python run.py excerpt` is tooling, not a stage.

`--offline` skips the stages that touch the network. Only `protein-set` and `composition`
have any: nothing else in the pipeline downloads (D80).

## What to rerun after an edit

**The rule: rerun the command that reads what you changed, then everything after it in the table above.**

| You changed | Rerun from |
|---|---|
| a literature file added, replaced, or re-obtained | 1, then 4 |
| `literature_sources.ini` (a source, a `used_for`, a hash) | 1, then 4 |
| `protein_set_decisions.ini` | 2 |
| `composition_decisions.ini` | 3 |
| `mass_fraction_decisions.ini` (column map, digest rule) | 4 |
| `fao_2013_…`, `fiber_type_mix`, `non_protein_metabolite_pools`, `uncertainty_settings`, `stress_test_settings` | 5 |
| `usda_food_data.ini` or a USDA archive | 6, then 8 |
| `gorissen_2018_comparison.ini` | 7 |
| `match_rate.ini` | 8 |
| a stage's code | that stage's command, then everything after it |
| documents only | nothing |

**The chain that catches people: 5 → 7 → 8.** `comparison` and `match` both record the hashes
of the standard's tables in their headers, so any rerun of `standard` makes both stale. Running
`standard` alone and stopping leaves `match` claiming inputs it no longer has.

Each command runs the tests at the end, and a currency test per stage checks that what is on
disk equals a fresh build. That is what tells you a rerun is owed — except for `comparison`,
which has no currency test yet (see `docs/step7a_report.md`, open items).

## What does not force a rerun

- Re-running `fetch-literature` on unchanged files. A stage hashes only the manifest entries it
  reads, and retrieval times are excluded from that hash, so re-hashing a file that has not
  moved changes nothing downstream.
- A source added to `literature_sources.ini` for another category. A stage that does not read it
  is not affected.

## A full rebuild

In order, 1 through 8. On a machine that already has `data/` populated, use `--offline` on 2 and
3 to skip roughly two hours of UniProt fetching.
