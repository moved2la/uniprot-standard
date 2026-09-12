# uniprot-standard

An open, sequence-derived amino acid standard for human skeletal muscle protein,
computed from UniProtKB sequences and literature-derived tissue mass fractions.

Read `PROVENANCE.md` first. Every value in this repository is an author decision, a
citation, or a computation, and nothing else.

## Run

```
pip install -r requirements.txt
python run.py protein-set            # protein set (needs network for three stages)
python run.py protein-set --offline  # rebuild config from data/ on disk
python run.py composition            # per-protein composition (needs network for two stages)
python run.py composition --offline  # recompute from data/ on disk
python run.py test
```

Every stage writes a timestamped log to `logs/` (debugging only; not committed). Cases the rules cannot settle
go to `outputs/flags.tsv` and are closed, with a D-number, in
`config/protein_set_decisions.ini`.

## Layout

| Path | Meaning |
|---|---|
| `config/` | Hand-written decisions (`protein_set_decisions.ini`, `composition_decisions.ini`) and generated config |
| `data/` | What the public databases said |
| `outputs/` | Everything computed here that is not config |
| `logs/` | Stage and test logs, debugging only; not committed |
| `pipeline/` | Code, one module per stage |
| `docs/` | Plan, decisions log, methods, conventions, handoffs |
| `tests/` | Offline tests with synthetic fixtures |

See `docs/conventions.md` for the rules the code applies and the field each reads.
