# Conventions

## Hand-written files

Exactly these files are written by a person. Everything else in the repository is
fetched or generated, and `pytest` fails if a generated file differs from what its
generator produces from the data on disk.

| File | What it may contain |
|---|---|
| `config/protein_set_decisions.ini` | Per tier: the seed word and either the exact ontology lookup name (D13) or `definition = measured_remainder` with the dataset and contaminant-list source keys (D56), each with its D-number; the ontology source URL; closures of flags raised by the code. Never a protein, accession, or term ID chosen by a person. |
| `config/composition_decisions.ini` | The URL of the IUPAC-IUBMB table that defines the one-letter code, the PubChem endpoint, the chiral-prefix query rule, and the two extra compound names (water, hydrogen), each with its D-number (D26). Never an amino acid name, symbol, or mass — the code parses the names from the cited table and fetches every mass. |
| `config/mass_fraction_decisions.ini` | Literature sources (citation, DOI, role with its D-number, download URLs), the column-to-fiber-type map for the primary file (each value the exact header string as it appears in the sheet, cited by sheet and row), and the in-silico digest rules with their citation. Never a protein, an accession, an abundance number, or a non-human data source — a test asserts the first three. |
| `config/carroll_classical_fractionation.ini` | The measured values of Carroll, Carrithers & Trappe 2004, transcribed by the author with page and table per row (D50). Feed the classical cross-check only, never a weight. |
| `config/fao_2013_indispensable_amino_acids.ini` | Which amino acids are dietary indispensable and how the reference scoring pattern groups them, transcribed by the author from FAO 2013 with table and page per value (D62). The report's three-letter symbols as printed; resolved to one-letter symbols by code from the IUPAC-IUBMB table. Never an abundance, a weight, or a protein. |
| `config/fiber_type_mix.ini` | The share of a muscle's protein in each fiber type, transcribed by the author from a cited source with page and table, and its basis (area, count); the only input to the final column of the standard. `___` until filled; the column stays blank. |
| `config/non_protein_metabolite_pools.ini` | The bound metabolite pool folded into the muscle standard (D73–D77): its amino acid as the three-letter symbol the sources print (resolved by code), its stoichiometry, its cited size per fiber type and per whole muscle with basis (dry / wet), the cited protein and water content per kg muscle, and the cited turnover rates for the resupply-frame sensitivity — transcribed by the author with location per value (the D50 pattern). Never a weight, a protein, or a mass. `___` until filled; the adjusted table is not written and a test says so. |
| `config/gorissen_2018_comparison.ini` | The human muscle column of Gorissen et al. 2018's Table 1, transcribed with a page-level location on every value, plus the paper's own statements about what its method converts and cannot measure (D84). Values and row labels only — no amino acid is named in code, and the comparison set is derived from the [hydrolysis_conversions] and [not_measured] sections rather than chosen. |
| `config/usda_food_data.ini` | Which USDA FoodData Central archives are read, the release of each, and the name each extra nutrient (protein, nitrogen, cystine, hydroxyproline) goes by in the archive's own `nutrient.csv` (D83). Never an amino acid name, a nutrient id, a food, or a value — the twenty are joined through the IUPAC-IUBMB trivial names and a test fails if one is named here. |
| `config/match_rate.ini` | What Match Rate scores and against what: the references (a file and column, or a config section), each with its scored set and label; the scored sets (the FAO headings with each group's scored member named as a three-letter symbol, D86; or a list of rows another config file names); the two food tables; the version label (D86–D88). |
| `config/food_amino_acids_other_sources.csv` | Foods USDA does not carry, maintained by the author (M4). CSV, one row per food. Columns: `label` (the food's identity; unique), `source` (where the values come from — a supplier analysis, a label, a paper; required), `location` (the spot in the source), `note`, `protein_g_per_100g` (carried, not scored), then `A_g_per_100g` … `Y_g_per_100g` — the twenty amino acids by one-letter symbol, grams per 100 g of product, blank where not reported. A header-only file is valid. |
| `config/uncertainty_settings.ini` | Monte Carlo draws and seed (D63). Nothing biological. |
| `config/stress_test_settings.ini` | The magnitudes the stress stage pushes the weights by (D65). Nothing biological. |
| `docs/*.md` | Prose. Every sentence in `docs/methods.md` either carries a citation or describes a computation performed here (`PROVENANCE.md`). |
| `tests/**` | Code and synthetic fixtures. Fixtures contain no real biology. |
| `pipeline/**`, `run.py` | Code. Contains no gene name, accession, term ID, or category. |

**Location format for every transcribed value (standing rule from Step 4b, D79):** `location` is a page number plus
the exact spot — table row and column, or column and paragraph — so a person opens the PDF and lands on it
(`p. E368, Table 1, row "Proteins, %", column "Control Group (n = 16)"`; `p. 841, first column, first paragraph`).
The file on disk goes in `file`, the source's wording in `as_reported`, and anything else about the study in `note`;
`location` carries nothing but the place. Values may be read from the PDF by the assistant and audited by the author;
`transcribed_by` says so.

## Placeholders

`___` means the author fills this in. It is only ever a value a person types, never something
code writes and never a hash: an unknown hash is left blank. Any stage that meets `___` in a
field it needs stops and flags; it never defaults. A test walks every `.ini` under `config/` and
fails on `___` in any `sha256` key.

The same principle covers defaults in code: where the word printed or the column shown is a
decision, config says it or the stage stops. `[categories] unused_marker` has no code fallback
for that reason.

## Folder meanings

One meaning per folder. A new category (blood, liver) repeats this shape under its own name:
`config/<category>/` and `outputs/<category>/{intermediate,standard}/`.

| Folder | Meaning |
|---|---|
| `config/` | Decisions (hand-written) and the config generated from them. `literature_sources.ini` holds every source the repository cites, one section each, with `used_for` naming what it is used for. The Layer B weights are generated in two forms with identical values: `config/mass_fractions/<type>.ini` (one section per accession, the per-entry citable record) and `config/mass_fractions_per_entry.tsv` (one row per accession, every fiber type side by side, provenance once in the header — the form the aggregation stage reads; D61). |
| `data/` | What the public databases said, unchanged or minimally tabulated. `data/gene-ontology/`, `data/uniprot_raw/`, `data/uniprot_sequences.ini`, `data/iupac/`, `data/pubchem/`, `data/uniprot_ptmlist/`, `data/usda/`. |
| `data/literature/` | Not committed (`.gitignore`). The literature files as obtained, unchanged, one folder per source id. Generated beside them: `manifest.ini` (url, path, sha256, bytes, retrieval time), `manifest_map.tsv` (one row per source, one column per purpose, `X` where used), `README.md` (human-readable provenance). Nothing here is downloaded by code (D80); `python run.py fetch-literature` hashes what is on disk. |
| `outputs/intermediate/` | What feeds the standard but is not the standard, one folder per stage: `protein_set/`, `composition/`, `digest/`, `literature_inventory/`, `mass_fractions/`. |
| `outputs/standard/` | The standard. The deliverable tables and the profiles at the top; `uncertainty/`, `stress/` and `sensitivity/` hold the tables that support them, each with its own `plots/`; `plots/` at the top holds the plots of the standard itself. `docs/pipeline_map.md` says which file answers which question. |
| `outputs/match/` | Match Rate. The per-food and cross-reference tables at the top; `match/<reference>/` — one folder per `[reference.*]` section of `config/match_rate.ini` — holds that reference's ranked table. |
| `outputs/usda/`, `outputs/comparison/` | The food side; the standard beside Gorissen 2018's measurement. |
| `outputs/flags.tsv` | Cases the rules could not settle. Written only when a flag is raised. |
| `logs/` | Debugging only: one timestamped log per stage per run, plus one per test run, plus one master log per command holding exactly what appeared on screen. Not committed; not the record. The record is the header of each generated file plus `outputs/flags.tsv`. |
| `docs/` | `run_order.md` (what to rerun after an edit), `pipeline_map.md` (every stage, its inputs and outputs, the flowchart), plan, decisions, methods, conventions, handoffs. |

## Paths in code

Every generated path is resolved through `pipeline/common.py`. Paths derived from a root that
can move — `OUTPUTS_DIR`, `DATA_DIR` — are **functions, not constants**: a path computed at
import time cannot follow a test that redirects the root, or a category writing under
`outputs/<category>/`. `tests/test_paths.py` fails the moment one goes back to being frozen
(D91). Where a `standard/` file lives is data, not code: `common.STANDARD_SUBFOLDERS` maps a
file name to its subfolder and `standard_path(name, base)` resolves it against whichever
standard folder the caller is writing (D90).

## One command

Eight commands, in run order: `fetch-literature`, `protein-set`, `composition`,
`mass-fractions`, `standard`, `usda`, `comparison`, `match`. Each runs its stages in order and
then the tests. `--offline` skips the stages that touch the network; only `protein-set` and
`composition` have any, because nothing in the pipeline downloads (D80).

`python run.py test` runs the tests alone. `python run.py excerpt` is tooling, not a stage: it
writes bounded, labelled cuts of the large generated files into `excerpts/` (not committed, not
the record) for review, one spec per command.

**`docs/run_order.md` is the one screen that says what to rerun after an edit.** The rule is:
the command that reads what changed, then everything after it.

Nothing in the repository is named by a project-plan step number.

## Rules applied by the code

Each rule names the database field it reads. Where the field cannot settle the case,
the code raises a flag; it never defaults.

**No stage deletes, renames, or moves a file** in `config/`, `data/`, or `outputs/`. A stage writes
the files it generates and leaves everything else as it found it; when a generated file is
superseded, removing it is the author's own hand action, recorded in `decisions.md`, never code's (D61).

| Rule | Reads | Result |
|---|---|---|
| **R0 Category** | Gene Ontology term **name**, exact string equality with `lookup_name` | Exactly one non-obsolete match → the tier's term. Zero or several → stop and flag. No related-term search. |
| **R0a Subtree** | Ontology relations `is_a` and `part_of` | All descendants of the term; immediate parents and children recorded as neighbours. |
| **R0b Pool** | UniProtKB search fields `organism_id`, `reviewed`, `go` | Union over one query per subtree term of `(organism_id:9606) AND (reviewed:true) AND (go:<id>)`. The term-level query is run as a check; a difference is flagged. |
| **R1 Sequence** | The entry's canonical sequence (the sequence UniProt displays and publishes an MD5 for) | Always used. No isoform is ever selected. Isoform composition deltas are computed and disclosed, not acted on. |
| **R2 Processing** | UniProt feature table, type `Chain` | R2a exactly one → its range is the master molecule. R2b none → whole sequence, logged. R2d several → the union of their ranges (D24), listed in `outputs/multi_chain_entries.tsv`. R2c any Chain with a non-exact position → flag. |
| **R3 Evidence** | Gene Ontology cross-reference property `GoEvidenceType` | Recorded per annotation and summarised per tier. Never filters. |
| **R4 Tier** | Membership of `data/tier<N>_pool.tsv` | An entry in more than one pool carries every tier it is in. |

| **R2e Non-exact Chain** (D60) | `Chain` positions with a qualifier or `UNKNOWN` | Qualified positions used as stated; an `UNKNOWN` boundary extends to the sequence end; listed in `outputs/chain_positions_resolved.tsv`. Supersedes R2c's flag. |
| **R5 Alphabet** (D59) | the canonical sequence's letters against the twenty coded amino acids | Any other letter excludes the entry from the set, listed with UniProt's molecular weight in `outputs/excluded_non_standard_alphabet.tsv`; an isoform with another letter is skipped in the delta table and listed. |
| **M1–M5 Measured tier** (D56) | The identity column of the primary Layer B dataset; UniProt search fields `organism_id`, `reviewed`, `gene_exact`; the contaminant FASTA headers | M1 genes = every identity cell split on `;`. M2 each gene looked up by `gene_exact`; an entry belongs to the gene when the gene equals its primary symbol or a synonym. M3 exactly one entry → candidate; zero → `unmapped`; several → `ambiguous`, excluded. M4 a candidate in an ontology tier → `in_tier<N>`, not a member. M5 a candidate in the contaminant list → `contaminant`, not a member. Every outcome is written to `data/tier<N>_gene_mapping.tsv`. R3 does not apply to a measured tier. |

There is no fiber-type rule. Fiber type is a Layer B quantity (per-fiber-type
abundance) and does not appear in the protein set.

### Composition rules

| Rule | Reads | Result |
|---|---|---|
| **C0 Symbols** | The IUPAC-IUBMB Table 1 page text: rows of the form `<trivial name> <three-letter symbol> <one-letter symbol>` | Exactly the twenty standard letters, each once → `data/iupac/amino_acid_symbols.ini`. Anything else → stop. No name is typed. |
| **C1 Masses** | PubChem `MolecularWeight`, `MolecularFormula`, `CID` for the query name (chiral prefix + trivial name; bare trivial name if PubChem has no such compound, recorded per row) | Exactly one compound → its served values, unchanged. Zero or several → stop. residue mass = free mass − water. |
| **C2 Counts** | The canonical sequence and `config/segments.ini` | `master` = residues with `in_master_molecule = true`; `metabolic` = every residue. Counts per letter are the ground truth (D7); both mass vectors and both fraction vectors are arithmetic on them. A letter outside the twenty standard ones stops the run. |
| **C3 PTM disclosure** | UniProt features of type `Modified residue`, `Lipidation`, `Glycosylation`, `Cross-link`, `Disulfide bond`; `ptmlist.txt` field `MA` | Description's first clause = vocabulary `ID` with an `MA` line → that average mass delta. Intrachain `Disulfide bond` → −H₂. Anything else → counted, no mass, listed. Sites outside the master molecule → recorded, excluded from sums. Nothing is decided (D28). |

### Literature rules (`pipeline/fetch_literature.py`)

Nothing here downloads (D80). A literature file is obtained by hand once, saved under
`data/literature/<source_id>/`, and hashed once; every run re-hashes what is on disk and must
reproduce the pin. **F3, F4, F4b and F6 retired with the download code** in Step 7a.

| Rule | Reads | Result |
|---|---|---|
| **F1 Placeholder** | `file.<n>.url` | A url of `___` is a flag: recorded, skipped, non-zero exit at the end. The stored filename comes from the url, so a missing one cannot be worked around. |
| **F2 Hash pinned** | `file.<n>.sha256` | Blank → the computed hash is written back into `config/literature_sources.ini` (the only thing this stage writes into config). Present → the bytes on disk must hash to it; a mismatch is a flag. |
| **F2b Blank, not `___`** | `file.<n>.sha256` | An unknown hash is left blank. `___` is only ever a value a person types; a test fails on `___` in any `sha256` key under `config/`. |
| **F5 On disk** | `data/literature/<source_id>/` | Every file must already be there; it is hashed in place and the file on disk is the record. A file that is missing, or whose bytes are an HTML page while its name implies a document, is a flag — that is what catches a saved paywall page filed under a `.pdf` name. |
| **F7 Stored name** | `file.<n>.name` — the filename as downloaded (amended 2026-09-18, Step 7b; D-number pending) | The URL is provenance only. A section with no `file.<n>.name` keeps the old rule (URL basename, Windows-illegal characters → `_`). A URL whose path ends in `/` has no filename (ACS supplement links), so its file needs `file.<n>.name`. The name looked for is kept in the manifest as `served_name`, with `name_from = config` or `url`. Files are never renamed to fit a URL. |
| **F7b Extension** | the file's magic bytes | When the URL basename has no extension, the one the bytes imply is appended (`%PDF` → `.pdf`, `PK` → `.zip`, `Rar!` → `.rar`, gzip → `.gz`), so every stored file opens in the program a person would use; the file is found under the bare name or under the name with any of those extensions. |
| **F8 Cited, not stored** | the absence of `file.<n>.*` | A source with no files is cited but not stored — a method citation, a provenance reference, a paywalled article whose values were never transcribed. It appears in the map and the README and is never flagged for a missing file. |
| **F9 Used for** | `used_for` | Names what the source is used for, comma-separated, from the columns declared in `[categories]`. A value not in that list is a flag, so a typo cannot quietly drop a source out of the map. A source that feeds no calculation carries the `unused_marker`, which prints in the `other` column — a row of blanks would read as an unfilled line. |

### Digest rules (`pipeline/digest.py`)

| Rule | Reads | Result |
|---|---|---|
| **G1 Cited or stop** | `[digest]` in `mass_fraction_decisions.ini` | Every parameter filled and `source` non-empty, else stop. |
| **G2 Cleavage** | `cleave_after`, `cleave_before_proline` | Cleave after each listed residue; `both` computes with and without cleavage before proline, primary = cleave-before-proline, difference per entry written. |
| **G3 Window** | `min_length`, `max_length`, `missed_cleavages` | Distinct peptide sequences per entry inside the window, with up to the allowed missed cleavages. |
| **G4 Sequence** | the full canonical sequence; composition table row `segment_set = metabolic` for MW | The digest runs on what a search engine digests, not the master molecule. |
| **G5 Nothing named** | — | No accession in code; nothing filtered. |

### Literature inventory rules (`pipeline/literature_inventory.py`)

| Rule | Reads | Result |
|---|---|---|
| **I1 Hash** | `data/literature/manifest.ini` `sha256` | Every file re-hashed; mismatch or missing file stops the run (B0). |
| **I2 In memory** | archive members | Nothing is written under `data/literature/`; member hashes go to the output table and the log. |
| **I3 Extractor** | PATH, default Windows install folders, `--rar-tool` | A `.rar` with no extractor found is a stop, stated in the log. |
| **I4 Preview** | sheet cells | Exact row and column counts; the first rows as text, truncated per cell; nothing interpreted. |
| **I5 Nothing named** | — | Records only. |

### Mass-fraction rules (`pipeline/mass_fractions.py`)

| Rule | Reads | Result |
|---|---|---|
| **B0 Source** | `config/mass_fraction_decisions.ini` `[source.*]`; `data/literature/manifest.ini` | Every file read is first re-hashed and compared with the manifest; mismatch stops the run. |
| **B1 Identity** | `identity_kind = gene_name`; the `gene` column of `config/accessions.ini`; `data/tier<N>_gene_mapping.tsv` | Exact symbol match after stripping whitespace; B1b (D57): a dataset symbol UniProt resolves to a pool entry joins that entry, marked `via_mapping`. Several dataset rows per gene → summed (D47). Multi-gene cells → split (D48). A pool gene with no row → w = 0, listed. A dataset gene with no pool entry → completeness list. |
| **B2 Quantity** | the mapped median column; composition table `mw` on `segment_set = master` | w_i = v_i · MW_i / Σ_tier (v_j · MW_j), within tier (D31). |
| **B3 Fiber type** | the dataset's own columns and purity threshold | Never re-derived. |
| **B4 Missing** | the mapped valid-values column | `NaN` with valid values = 0 → 0, listed; never imputed. |
| **B5 Spread** | the mapped SD column per fiber type | Carried per row for the aggregation step (D45); w_low / w_high = (v ∓ SD) × MW over the unchanged denominator, clipped at zero; a shared row adds its whole value to each member's w_high (D48). Rows summed under D47 carry the root-sum-square of their SDs (D53). |
| **B9 Config** (D61) | every measured entry | The weights written in two forms with identical values: `config/mass_fractions/<type>.ini` (one section per accession, within-tier and combined weights with low / high) and `config/mass_fractions_per_entry.tsv` (one row per accession, every fiber type side by side, provenance once in the header). Both are covered by the currency test. |
| **B7 Families** | `outputs/digest/shared_pairs.tsv`; composition `free_frac_*` on master rows | A family at the D52 cutoff is a unit (D55); its bound per amino acid = (max − min free-mass fraction across members) × family weight. |
| **B8 Combined** (D58) | every measured entry of every tier | w_i = v_i · MW_i / Σ_all (v_j · MW_j), one denominator per fiber type; the primary standard. Within-tier weights are kept for the checks. |
| **B6 Bands** | `config/carroll_classical_fractionation.ini` `anchor_genes`; `outputs/digest/shared_pairs.tsv` | A gel band is the shared-peptide family of its anchor genes at the D52 cutoff; its weight is Σ w over the family. Families at cutoffs 1 and 2 are both written. |

### Aggregation rules (`pipeline/aggregate.py`, `uncertainty.py`, `plots.py`)

| Rule | Reads | Result |
|---|---|---|
| **A0 Inputs** | every file read | Hashed into every output header. Amended in Step 7a (D92): a stage hashes only the config sections and manifest entries it **reads**, not whole shared files — a source added to `data/literature/manifest.ini` for another category must not make a stage's header claim its inputs moved. Keys that record *when* a line was written rather than *what* it says (`retrieved`, `generated`) are excluded, so re-hashing an unchanged file does not mark everything downstream stale. A named section that is absent hashes as absent rather than being skipped. |
| **A1 Profile** (D64) | `config/mass_fractions_per_entry.tsv` weights; composition `count_*` on `segment_set = master`; the fetched masses | The molar mixture: n_a = Σ_i (w_i / MW_i) count_{i,a}; p_a = n_a m_a / Σ (m = residue or free mass). A1c: the fraction mixture, renormalised, written as a check. Wide tables: g per 100 g protein, unnormalised. |
| **A2 Sets** | `w_<type>_combined` (D58); `w_<type>_tier1`, `w_<type>_tier2` (D31) | total = every weighted entry under one denominator; contractile and builders = the two tiers under their within-tier weights. No other set is named. |
| **A3 Uncertainty** (D63) | `v`, `sd`, `valid` per entry and fiber type; `[monte_carlo]` draws and seed | Log-normal per entry, exact moment match; between-fiber spread (σ) and median uncertainty (σ/√n) as two terms; every draw a composition; 2.5 / 50 / 97.5 percentiles. |
| **A4 One scale** (provisional) | the tables above | Every term as the maximum absolute shift of the fraction from the standard: Monte Carlo half-width; bounds as their worst-case sum over entries or families (largest single also written); sensitivity as max |profile − standard|; completeness as the outside genes' mass share. The widest term is named per amino acid. |
| **A5 Bracket** | the three pure-type profiles | max − min per amino acid bounds every mix of the types; the IIa − IIx difference is written (D46). |
| **A6 MHC:actin sensitivity** (D54b) | `classical_check_carroll_2004.tsv` at cutoff 2; `band_families.tsv` at cutoff 2 | At each ratio the methods measured, the MHC family is rescaled with actin fixed and the actin family with MHC fixed, all weights renormalised; combined and contractile; types I and IIa only. No method is a reference; files named `sensitivity_mhc_actin_*`. |
| **A7 Completeness** | `dataset_rows_outside_pool.tsv` | The outside genes' molar share by rank, converted to a mass share with the pool's molar-mean MW (stated assumption, under which the two are equal): a share of unknown composition moves no fraction by more than itself. |
| **A9 Stress** (D65) | `config/stress_test_settings.ini`; weights, counts, masses; `band_families.tsv`; `outputs/digest/theoretical_peptides.tsv` | Perturb the weights as stated, renormalise, recompute by A1, report the shift as a fraction and as a percentage; scenarios by rank, tier, band family, or factor — nothing named. TPA weighting = v × N_peptides in place of v × MW. |
| **A10 Non-protein metabolite pools** (D73–D77; `pipeline/non_protein_metabolite_pools.py`) | `config/non_protein_metabolite_pools.ini`; `amino_acid_g_per_100g_protein_free.tsv`; `_calculated_amino_acid_standard.tsv`; `fiber_type_mix.ini`; the fetched masses and symbols | Per fiber type, on the pool's mass basis: F_a = P × g100_a / 100 (free mass of amino acid a from protein per kg muscle; P = protein per kg on that basis, moved between bases only through the cited water content, P_dry = P_wet / (1 − W/1000)); C = c × n × m / 1000 (the pool's amino acid per kg; c in mmol/kg, n mol per mol, m the free mass). The pool's amino acid becomes F + C, every other stays F, all renormalised to 100. Contractile and builders are protein-only by definition and copied unchanged; the final column mixes the adjusted totals by the fiber-type shares as the protein-only standard does. The protein-only table is never changed. Sensitivities: the resupply frame (C × k_pool / k_protein, k_pool = ln 2 / (7 t½) per day, k_protein = 24 × FSR / 100 per day, over every cited pair), whole-muscle vs single-fiber basis (two options, no verdict), subgroups of the whole-muscle value, and the single-fiber spread (log-normal fitted to mean and SD, 2.5 / 97.5 percentiles). Any `___` leaves the adjusted table (or one sensitivity) unwritten and says so in `non_protein_metabolite_pools_summary.ini`. **Green light (D78):** before computing, every `source_id` the config names must have a hashed file in `data/literature/manifest.ini`; otherwise the stage stops. |
| **A8 EAA** (D62) | `config/fao_2013_indispensable_amino_acids.ini`; `data/iupac/amino_acid_symbols.ini` | Headings resolved from three-letter symbols by code; group headings expanded from their `[group.*]` sections; a `___` leaves the EAA subset unwritten, recorded in `standard_summary.ini`, and a test fails until it is filled. |

### USDA rules (`pipeline/usda.py`)

| Rule | Reads | Result |
|---|---|---|
| **U1 Archives** | `config/usda_food_data.ini`; `data/usda/*.zip` | Only the data types config names. A zip is identified by reading its own `food.csv`, never by its file name. The SHA-256 is the archive's identity: it goes in `data/usda/manifest.ini` and in every output header. Two archives of one data type, or a named data type with no archive, is a `[STOP]`. |
| **U2 Membership** | `<data_type>.csv` inside the archive | A food counts when the archive's own membership file lists it. FoodData Central ships superseded earlier versions of a food in `food.csv`; they are not in the membership file. No membership file → the `data_type` column of `food.csv`, and the log says which rule was used. |
| **U3 The name join** | `nutrient.csv`; `data/iupac/amino_acid_symbols.ini` | A nutrient is amino acid X when its name equals, ignoring case and space, the IUPAC-IUBMB trivial name of X. No amino acid is named in code or config. Nutrients carried that are not one of the twenty are named in config by the archive's own name. `usda_nutrient_map.tsv` is the audit trail. |
| **U4 Values** | `food_nutrient.csv` | The archive's `amount` per 100 g of food, as published; `min`, `max`, `median`, `data_points` carried where present. Nothing converted, scaled, or filled in. A carried nutrient in a unit config does not expect is a flag and is not carried. |
| **U5 Completeness is stated, not decided** | the pivoted rows | USDA reports no asparagine and no glutamine — acid hydrolysis converts them — so no food carries all twenty. Each row states how many it carries and which letters it lacks; which amino acids a score requires is the match step's decision. Foods with no amino acid row at all are listed in `foods_without_amino_acids.tsv`, not silently dropped. |
| **U6 No cross-archive dedup** | — | A food measured in two archives is two rows with two `fdc_id`s and its `data_type` named. `precedence` in config is recorded, not applied. Which to prefer is a match-step decision. |
| **U7 Nothing is removed** | — | This stage deletes, renames, and moves nothing in `data/` or `outputs/` (D61). |

### Comparison rules (`pipeline/comparison.py`)

| Rule | Reads | Result |
|---|---|---|
| **X1 One denominator** | the transcription; both standards | Both sides are renormalised over the SAME set before anything is compared. The measured side's published percentages use a denominator that includes what it did not measure, so they are never compared with this repository's directly. |
| **X2 The set comes from the method** | `[hydrolysis_conversions]`, `[not_measured]` | An amino acid the method converts into another is compared as the sum of the two, on both sides. One it destroys or does not report is excluded from both. The set is never chosen by which amino acids happen to agree. Every one of the twenty must be either compared or excluded with a stated reason, or the stage stops. |
| **X3 No verdict** | — | Difference and ratio only. Neither side is called correct, neither is fitted, scaled or adjusted toward the other, and no error term is written. A test asserts the column names carry no verdict word. |
| **X4 The transcription is checked against the paper** | the paper's own printed sums; `[values.human_muscle.location] essential_rows` (the paper's own statement of which rows its essential sum covers) | Recomputed from the transcribed values before anything else. A mismatch is a `[STOP]`, not a flag: a mistyped digit must not reach a comparison. |
| **X5 A zero is never compared unexamined** | `[not_measured]`, `[zero_is_a_measurement]` | A measured value of 0.0 must be declared: either in `[not_measured]` with the reason it is not a measurement, or in `[zero_is_a_measurement]` with the evidence that it is one. An undeclared zero is a `[STOP]`. A zero that reaches the comparison sits in one side's denominator and not the other's, and biases every other row (D85). |
| **X6 The mass balance is computed** | the PubChem masses; the calculated standard; every protein content config offers | Hydrolysis adds a water per peptide bond, so a gram of protein yields more than a gram of free amino acids. The stage computes that factor from this repository's own composition, then reports what share of the expected yield the source's published values reach, under each protein content offered. It draws no conclusion from the result (X3). |

### Match Rate rules (`pipeline/match.py`)

| Rule | Reads | Result |
|---|---|---|
| **M1 The score** | the food's and the reference's scored values | Both expressed as shares of their scored total; ratio = food share / reference share per amino acid; Match Rate = the smallest ratio, its amino acid the limiting one (D87). Protein content and the reporting unit never enter: the score is a property of the proportion. |
| **M2 The scored set is named, not coded** | `config/match_rate.ini` `[scored_set.*]`, `[reference.*]` | Every reference names its scored set. A set is either the FAO 2013 headings with each group's scored member named as a three-letter symbol (D86), or a list of rows another config file names (the paper's own essential rows for the original reference). Resolved through the IUPAC-IUBMB table; a group without a named member is a `[STOP]`. |
| **M3 Missing is listed, zero is scored** | the food's scored values | A food that does not report a scored amino acid goes to `foods_not_scored.tsv` with the letters named (U5 carried). A published zero is scored as published: ratio zero, score zero, limiting amino acid named. The stage does not decide what a source measured. |
| **M4 Two food tables, no dedup** | `outputs/usda/amino_acids_per_food.tsv`; `config/food_amino_acids_other_sources.csv` | Every row names its source; the same food from two sources, or two archives, is two rows (U6 carried). An other-sources row without a `source` is a `[STOP]`; a duplicate `label` is a `[STOP]`. |
| **M5 Every score can be rebuilt by hand** | — | `match_rate_steps.tsv` writes the spreadsheet's own walk per food and reference — the grams (rows 12–20), ΣEAA and "% of Human" (22–23), Step 3 scaled values (30–37), Step 7 percent of reference and its MIN (100–108), Step 10b values (145–152), Total Need to Consume, Percent Wasted, Percent Utilized (153–156) — in the spreadsheet's order and names, with the reference's values in the header (its column A). It is the same number as the minimum ratio (`docs/formula.md`); no second score is written. |
| **M6 One scorer** | — | `match_rate()` in `pipeline/match.py` is the only implementation of the formula. The blend / fortification script (D88) imports it; nothing copies it. |

## Flags

`outputs/flags.tsv` lists every case a rule could not settle, with the stage, subject,
rule, and the field contents that caused it. A flag is closed only by a `[flag.<subject>]`
section in `config/protein_set_decisions.ini` carrying a D-number. Closed flags stay in
the file; `pytest` fails while any flag is open.

## Ini style

`configparser`, one section per entity, keys aligned, `#` comments. Generated files
begin with `# GENERATED by <tool>. DO NOT EDIT BY HAND.` and a provenance header
(ontology file and version, UniProt release, fetch time). **The header is provenance, not
explanation** (author, 2026-09-17): what was generated when, the hash of every input, the
reference or source each value came from — the lines a reader needs to trace the file. What a
column means and why a rule exists lives in `docs/`, not in the header; one line naming the
column-to-source mapping is the most a generated table carries beyond provenance.

## Naming

**Non-protein metabolite pools (D82).** Amino acids held in stable non-protein compounds — carnosine today; creatine and glutathione later — as distinct from protein-bound amino acids (the standard) and from the free amino acid pool (excluded, D68). Defined in full in `docs/methods.md` §"Non-protein metabolite pools". The word "bound" is not used for them: in this literature it means protein-bound. Decisions D68, D73 and D75 predate the term and say "bound pool".

**Tiers and fiber types (D66).** The two protein tiers are **contractile** and **builders** in
every output and document; never "Tier 1 / Tier 2" in prose or headers. **Fiber** is the cell:
write "fiber type I / IIa / IIx", never "type I" alone. The identifiers `1` and `2` survive only
as values of `tier` columns in the config tables; `common.TIER_NAMES` maps them when a stage
writes a name.

**Config files** are named by source and content (`fao_2013_indispensable_amino_acids.ini`,
`fiber_type_mix.ini`) or by what they set (`uncertainty_settings.ini`,
`stress_test_settings.ini`), never by a step or a stage.

**The deliverable** is `outputs/standard/_calculated_amino_acid_standard.tsv`; the leading
underscore sorts it first in its folder. Its companion with the bound pool folded in,
`_calculated_amino_acid_standard_with_non_protein_metabolite_pools.tsv`, sorts beside it; the protein-only table is
never changed by the bound-pool stage.

Write "Gene Ontology" in full. The two-letter abbreviation appears only inside term
IDs (`GO:0000000`) and in the UniProt query field name `go:`, where it is a literal.
