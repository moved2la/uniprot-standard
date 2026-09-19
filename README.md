# uniprot-standard

An open, sequence-derived amino acid standard for human skeletal muscle protein, computed from
UniProtKB sequences and literature-derived tissue mass fractions.

The result is `outputs/standard/_calculated_amino_acid_standard.tsv`.

Every value in this repository is an author decision, a citation, or a computation, and nothing
else. Read `PROVENANCE.md` first.

## Install and run

```
pip install -r requirements.txt
python run.py fetch-literature
python run.py protein-set
python run.py composition
python run.py mass-fractions
python run.py standard
python run.py blood-protein-set
python run.py blood-composition
python run.py blood-mass-fractions
python run.py blood-standard
python run.py composite
python run.py usda
python run.py comparison
python run.py match
```

That is the full rebuild, in order (`docs/run_order.md`). `--offline` on `protein-set`,
`composition` and `blood-protein-set` skips the network and rebuilds from `data/` on disk; no
other command touches the network.

`python run.py test` runs the tests alone.

## Where to look

| | |
|---|---|
| **What to rerun after an edit** | `docs/run_order.md` — one screen |
| **What each stage reads and writes** | `docs/pipeline_map.md` — every stage, every file, the flowchart |
| **The rules the code applies, and the field each reads** | `docs/conventions.md` |
| **Why each choice was made** | `docs/decisions.md` — every decision, numbered, with its rationale |
| **The calculation by hand** | `docs/formula.md` |
| **How the standard compares with a measurement** | `docs/gorissen_comparison.md` |

Every stage writes a timestamped log to `logs/` (debugging only, not committed), and each
command writes one master log holding exactly what appeared on screen. Cases the rules cannot
settle go to `outputs/flags.tsv` and are closed, with a D-number, in config.
