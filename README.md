# uniprot-standard

An open, sequence-derived amino acid standard for human skeletal muscle protein,
computed from UniProtKB sequences and literature-derived tissue mass fractions.

The result is `outputs/standard/_calculated_amino_acid_standard.tsv`; `docs/formula.md` is the calculation by hand.
Read `PROVENANCE.md` first, then `docs/pipeline_map.md` for what runs in what order. Every value in this repository is an author decision, a
citation, or a computation, and nothing else.

## Run

```
pip install -r requirements.txt
python run.py protein-set            # protein set (needs network for four stages, incl. the literature fetch)
python run.py protein-set --offline  # rebuild config from data/ on disk
python run.py composition            # per-protein composition (needs network for two stages)
python run.py composition --offline  # recompute from data/ on disk
python run.py mass-fractions         # literature fetch, digest, inventory, mass fractions (needs network for one stage)
python run.py mass-fractions --offline
python run.py standard               # the standard: aggregate -> non_protein_metabolite_pools -> uncertainty -> stress -> plots (offline)
python run.py usda                   # the food side of Match Rate: read the USDA archives in data/usda/ (offline)
python run.py comparison             # the standard beside Gorissen 2018's measured human muscle composition
python run.py match                  # Match Rate: every food against every reference (offline)
python run.py test
python run.py excerpt                # tooling: bounded excerpts of the large generated files, for review (not the record)
```

Every stage writes a timestamped log to `logs/` (debugging only; not committed). Cases the rules cannot settle
go to `outputs/flags.tsv` and are closed, with a D-number, in
`config/protein_set_decisions.ini`.

## Layout

| Path | Meaning |
|---|---|
| `config/` | Hand-written decisions and transcriptions (`protein_set_decisions.ini`, `composition_decisions.ini`, `mass_fraction_decisions.ini`, `carroll_classical_fractionation.ini`, `fao_2013_indispensable_amino_acids.ini`, `fiber_type_mix.ini`, `non_protein_metabolite_pools.ini`, `gorissen_2018_comparison.ini`, `usda_food_data.ini`, `match_rate.ini`, `food_amino_acids_other_sources.csv`, `uncertainty_settings.ini`, `stress_test_settings.ini`) and generated config (`accessions.ini`, `segments.ini`, `mass_fractions/<type>.ini`, `mass_fractions_per_entry.tsv`) |
| `data/` | What the public databases said |
| `outputs/` | Everything computed here that is not config |
| `logs/` | Stage and test logs, debugging only; not committed |
| `excerpts/` | Bounded cuts of large generated files for review in chat (`run.py excerpt`); not committed |
| `pipeline/` | Code, one module per stage |
| `docs/` | Plan, decisions log, methods, conventions, handoffs |
| `tests/` | Offline tests with synthetic fixtures |

See `docs/conventions.md` for the rules the code applies and the field each reads.
