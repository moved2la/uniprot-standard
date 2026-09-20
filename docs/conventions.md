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
| `scoring/*.xlsx` | The author's own foods, scored by `python run.py score` and by nothing else (M7, D129). Gitignored: not part of the repository's record, so it is the one hand-written place a proprietary formulation may sit. One row per product; required `label`, `basis` (`product` or `protein`), `grams_basis`; optional `protein_g` and any columns of the author's own (source, location, note), carried through untouched; the amino acids by three-letter symbol or trivial name, in the same unit as `grams_basis`. Blank means not reported (M8); `0` means measured as zero (M3). |
| `config/tissue_mass_fractions.ini` | Layer C (D67, D121, D122): one `[category.<name>]` per category standard — the table and column the composite reads, the tissue mass transcribed from ICRP 89 with page, table, row and column, and how the protein mass is set: a cited protein content read from the config that already transcribes it (`protein_content_from`), or the category's own split table and base row. Never a profile, a weight, or a protein. |
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
| `outputs/composite/` | The composite (Layer C, D121): the category standards `config/tissue_mass_fractions.ini` names, mixed by protein mass; the deliverable `_calculated_amino_acid_standard_composite.tsv`, the masses, the contributions per category, the split sensitivity, the summary. Not the whole body until the categories justify the name. |
| `outputs/usda/`, `outputs/comparison/` | The food side; the standard beside Gorissen 2018's measurement. |
| `outputs/flags.tsv` | Cases the rules could not settle. Written only when a flag is raised. |
| `logs/` | Debugging only: one timestamped log per stage per run, plus one per test run, plus one master log per command holding exactly what appeared on screen. Not committed; not the record. The record is the header of each generated file plus `outputs/flags.tsv`. |
| `docs/` | `run_order.md` (what to rerun after an edit), `pipeline_map.md` (every stage, its inputs and outputs, the flowchart), plan, decisions, methods, conventions, handoffs. |

Since Step 7b: `config/<category>/` holds a category's decisions and generated config; `data/<category>/` its pool tables and the enumeration's row outcomes; `outputs/<category>/` its `intermediate/` and `standard/`; `data/uniprot_raw/<category>/` the raw entries that category's fetch added (S1). `data/uniprot_sequences.ini` is one store for every category.

## Paths in code

Every generated path is resolved through `pipeline/common.py`. Paths derived from a root that
can move — `OUTPUTS_DIR`, `DATA_DIR` — are **functions, not constants**: a path computed at
import time cannot follow a test that redirects the root, or a category writing under
`outputs/<category>/`. `tests/test_paths.py` fails the moment one goes back to being frozen
(D91). Where a `standard/` file lives is data, not code: `common.STANDARD_SUBFOLDERS` maps a
file name to its subfolder and `standard_path(name, base)` resolves it against whichever
standard folder the caller is writing (D90).

## One command

Thirteen commands, in run order (`common.COMMAND_ORDER`, D123): `fetch-literature`, `protein-set`,
`composition`, `mass-fractions`, `standard`, `blood-protein-set`, `blood-composition`,
`blood-mass-fractions`, `blood-standard`, `composite`, `usda`, `comparison`, `match`. Each runs
its stages in order and then the tests. `--offline` skips the stages that touch the network; only `protein-set` and
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
| **R4 Tier** | Membership of `data/tier<N>_pool.tsv`; for a category, of `data/<category>/pool_<name>.tsv` | An entry in more than one pool carries every tier it is in. For a category the `tier` field of `config/<category>/accessions.ini` holds the pool names (Step 7b; the field is not renamed until 7c). |

| **R2e Non-exact Chain** (D60) | `Chain` positions with a qualifier or `UNKNOWN` | Qualified positions used as stated; an `UNKNOWN` boundary extends to the sequence end; listed in `outputs/chain_positions_resolved.tsv`. Supersedes R2c's flag. |
| **R5 Alphabet** (D59) | the canonical sequence's letters against the twenty coded amino acids | Any other letter excludes the entry from the set, listed with UniProt's molecular weight in `outputs/excluded_non_standard_alphabet.tsv`; an isoform with another letter is skipped in the delta table and listed. |
| **M1–M5 Measured tier** (D56) | The identity column of the primary Layer B dataset; UniProt search fields `organism_id`, `reviewed`, `gene_exact`; the contaminant FASTA headers | M1 genes = every identity cell split on `;`. M2 each gene looked up by `gene_exact`; an entry belongs to the gene when the gene equals its primary symbol or a synonym. M3 exactly one entry → candidate; zero → `unmapped`; several → `ambiguous`, excluded. M4 a candidate in an ontology tier → `in_tier<N>`, not a member. M5 a candidate in the contaminant list → `contaminant`, not a member. Every outcome is written to `data/tier<N>_gene_mapping.tsv`. R3 does not apply to a measured tier. |
| **P1–P6 Pool by published accession** (D107, Step 7b) | a category's `[pool.<name>]` sections; the dataset's identity column by letter; UniProt search field `accession` with `reviewed`, `organism_id`, `keywordid` | P1 rows = every non-empty data row. P2 the group = the identity cell split on `;` in the published order; a token is a `CON__` marker, a `REV__` decoy, or an accession with its isoform suffix removed (R1). P3 every accession token looked up; a row belongs to a token when the token equals the entry's primary accession; a deleted or demerged accession comes back as an inactive stub (no reviewed status, no organism). P4 the row's entry = the first token in group order that is a reviewed human entry; first token a marker → `contaminant_group`. P5 an entry carrying the keratin keyword (`[contaminant_rules] keratin_keyword_id`, D109) → `keratin`; a gene-cell token among the pool's `named_contaminant_genes` (D98) → `named_contaminant`; both excluded and listed with their share. P6 members = the distinct entries of the surviving rows; a later row choosing an entry already chosen is `member_duplicate_entry` (summed by the join, D47). Every row's outcome → `data/<category>/<pool>_dataset_rows.tsv`; the removed rows → `<pool>_rows_excluded_with_share.tsv`. R3 does not apply. |
| **P3b Secondary accession** (D111) | for a token not returned as a primary, or returned inactive: one query per token by the `sec_acc` field, confirmed against the entry's `secondaryAccessions` in its JSON | Exactly one reviewed human entry lists the token → `member_via_secondary_accession`; several (demerged) → every one joins and the row is shared between them, split by D48 (`…_shared`); none → obsolete. A refused query field is logged and the row falls to P4b. |
| **P4b Gene fallback** (D111, D112) | a group with no usable token after P3b: its gene cell, gene by gene, through the M2 query and attribution with an entry's `;`-joined primary symbols counted one by one | M3 per gene; M3b breaks a synonym collision in favour of the one primary-symbol entry. One gene → `member_via_gene`; several genes each resolved → `member_via_gene_shared` (split by D48); any gene unresolved → `no_reviewed_human_entry`, listed with its share. |
| **M3b Primary over synonym** (D112) | the M2 attribution's `primary` / `synonym` mark | When several reviewed human entries answer a gene and exactly one answers by its primary symbol, that entry stands. Applied in P4b only; muscle's M3 unchanged. |
| **S1 Additive store** (D110) | `data/uniprot_sequences.ini`; `data/uniprot_raw/<category>/` | The store is merged, never rebuilt from one category's pool: sections fetched now are rendered, every other section carried over verbatim in one sorted order; a category fetches only what the store lacks; a section a category adds carries its own `fetched` and `uniprot_release`, the header gains `added.<category>.*` lines. The raw folder is nested one folder per category (`skeletal_muscle/`, `blood/`) named for the fetch that brought the entry in; readers look in the root and every subfolder. Nothing is deleted. |

There is no fiber-type rule. Fiber type is a Layer B quantity (per-fiber-type
abundance) and does not appear in the protein set.

### Composition rules

| Rule | Reads | Result |
|---|---|---|
| **C0 Symbols** | The IUPAC-IUBMB Table 1 page text: rows of the form `<trivial name> <three-letter symbol> <one-letter symbol>` | Exactly the twenty standard letters, each once → `data/iupac/amino_acid_symbols.ini`. Anything else → stop. No name is typed. |
| **C1 Masses** | PubChem `MolecularWeight`, `MolecularFormula`, `CID` for the query name (chiral prefix + trivial name; bare trivial name if PubChem has no such compound, recorded per row) | Exactly one compound → its served values, unchanged. Zero or several → stop. residue mass = free mass − water. |
| **C2 Counts** | The canonical sequence and `config/segments.ini` | `master` = residues with `in_master_molecule = true`; `metabolic` = every residue. Counts per letter are the ground truth (D7); both mass vectors and both fraction vectors are arithmetic on them. A letter outside the twenty standard ones stops the run. |
| **C3 PTM disclosure** | UniProt features of type `Modified residue`, `Lipidation`, `Glycosylation`, `Cross-link`, `Disulfide bond`; `ptmlist.txt` field `MA` | Description's first clause = vocabulary `ID` with an `MA` line → that average mass delta. Intrachain `Disulfide bond` → −H₂. Anything else → counted, no mass, listed. Sites outside the master molecule → recorded, excluded from sums. Nothing is decided (D28). |
| **C3b Alternatives** (D118) | priced PTM features of one entry: type, start, end | Features of the same type at the same position are alternatives (UniProt records the glycoforms observed at one site as separate features); one position contributes the largest |delta| among them; the number of features collapsed is recorded per entry (`n_alternative_features_collapsed`) and in the summary. Different positions or types add as before. |

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
| **F7 Stored name** | `file.<n>.name` — the filename as downloaded (amended 2026-09-18, Step 7b; D104) | The URL is provenance only. A section with no `file.<n>.name` keeps the old rule (URL basename, Windows-illegal characters → `_`). A URL whose path ends in `/` has no filename (ACS supplement links), so its file needs `file.<n>.name`. The name looked for is kept in the manifest as `served_name`, with `name_from = config` or `url`. Files are never renamed to fit a URL. |
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

### Mass fractions from published shares (`pipeline/mass_fractions_pools.py`, Step 7b)

For a category whose primary datasets publish a mass share per protein (blood: `share_kind = mass_share`, D105). The join was done by the pool enumeration; there is no molecular-weight step.

| Rule | Reads | Result |
|---|---|---|
| **Q1 Join** | `data/<category>/<pool>_dataset_rows.tsv` | A row joins the entry (or entries) the enumeration gave it. Rows of one entry summed (D47); a shared row split in proportion to the members' own rows, equally if none, the whole to each `_high` (D48); summed rows carry the root-sum-square of their spreads (D53a). |
| **Q2 Value** | the column map by letter (`[columns.<source>.<file>]`): `share_columns` or `abundance_column` + `abundance_scale` | `share_columns` → v = mean across replicates (D108), sd = sample SD across them, n = replicates, the per-replicate values carried as `v_<pool>_rep<k>`. `abundance_scale = log10` → v = 10^x; `replicate_columns` → sd of 10^run in linear space, n = runs present (D114), the published CV carried beside it. `share_kind` must be `mass_share`; anything else stops (D105: the molar kind is not built until a source needs it). |
| **Q3 Weight** | the pool's entries | w = v / Σ_pool v (D31 within pool: the listing renormalised after the contaminant rule); `_low` / `_high` = (v ∓ sd) / Σ, clipped at zero (B5). Written to `config/<category>/mass_fractions_per_entry.tsv` (D61). |
| **Q4 What is missing** | `<pool>_rows_excluded_with_share.tsv`; `pool_<pool>.tsv` against `accessions.ini` | `removed_by_rule` (D98, on purpose) and `completeness_gap` (rows no entry answers to; R5-excluded entries) as shares of the listing, in `excluded_rows_mass_share.tsv`; the gap is the A7 bound. |
| **Q5 No family bound** (D115) | — | A published-share pool's groups are the entry with its own isoforms and fragments; the multi-gene groups are shared rows under Q1. The summary carries the shared-row count. |
| **Q6 Immunoglobulin bound** (D113) | `[immunoglobulin_bound]`; the HGNC gene symbol; master lengths | A constant-region entry (locus prefix IGH / IGK / IGL, not a V/D/J segment gene) carries unmeasured variable-domain mass bounded by w × L_V / L_C, L_V the mean master length of the pool's V-region entries of the same chain; a chain with no V entry gets the loose bound (its whole share) and says so. V-region entries are measured; no bound. |
| **Q7 Cross-checks** (D116) | `[cross_checks]`; Hortin 2008's `.xls` by letter (xlrd); Geyer's deep dataset | Reported, never used: the ratio per protein and the plasma profile under Hortin's mg/L beside the pool's over the joined set, unjoined rows listed with the reason; the share of the deep dataset's signal absent from the primary; the mass share of entries present in more than one pool, per pool. |

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

### Standard from pools mixed by a measured split (`pipeline/aggregate_pools.py`, Step 7b)

| Rule | Reads | Result |
|---|---|---|
| **T1 Profile** | the pool weights (mass shares, D105); composition `free_frac_*` / `residue_frac_*` on master rows | profile[a] = 100 × Σ_i w_i × frac_i[a], both conventions; no molar step (the weights are already mass). |
| **T2 Split** (D100, D119) | `config/<category>/compartment_mix.ini` | plasma protein = plasma volume × plasma protein concentration; haemoglobin = red-cell volume × MCHC (the blood route reported as a cross-check); erythrocyte protein = haemoglobin ÷ its share; f_plasma from the two masses; standard = f_plasma × plasma + (1 − f_plasma) × erythrocytes. The base share is as measured — the summed w of the `[dominant.erythrocytes]` genes. Reference male; the female volumes reported as a line per share. The share ICRP's whole-blood protein implies is computed and listed (`share.<x> = computed`). |
| **T3 Dominant protein** (D102, D117) | `[dominant.<pool>]` genes and published shares | The dominant entries scaled to the share and every other entry to 1 − share, each group keeping its internal proportions (the two-piece compartment); the erythrocyte protein mass follows the share, the plasma mass stays; the standard at every combination, `standard_minus_base` beside it; the spread per amino acid split into haemoglobin-only and albumin-only. |
| **T4 Per replicate** (D108) | `v_<pool>_rep<k>` | The pool's profile under each replicate alone, under the mean (the weight) and under the median; the standard with that pool replaced; the spread across replicates. |
| **T5 Uncertainty** (D63) | v, sd, n per pool; `[monte_carlo]` draws and seed | Log-normal per entry, exact moment match; `spread` and `mean_uncertainty` (σ / √n); every draw recomputes the whole profile; the standard's draws combine the pools at the base split; 2.5 / 50 / 97.5 percentiles. |
| **T6 One scale** | the weighted bounds, the immunoglobulin bound, the completeness gap, the two sensitivities | Every term as the maximum absolute shift of the g / 100 g value: Monte Carlo half-widths; bounds and gaps mixed at the split (A7); the sensitivities as their spread; the widest term named per amino acid. |
| **T7 Drivers** | the weights and compositions | Per profile and amino acid, the five entries contributing most; histidine's thirty largest with cumulative share; a `[histidine]` section in `standard_summary.ini` with every number that bears on it. |
| **T8 No stress, no comparison** (D120) | — | The dominant-protein and per-replicate sensitivities are the category's stress scenarios; no external measured blood profile exists to compare with. |

### Composite rules (`pipeline/composite.py`, Layer C, D121)

| Rule | Reads | Result |
|---|---|---|
| **L1 Protein mass** (D122) | `config/tissue_mass_fractions.ini` `[category.<name>]`; the standard table and column each names | Every category profile is a percent of its amino acid mass; the composite weights each by the protein the category holds in the reference person: `tissue_mass_x_protein_content` = cited tissue mass × cited protein per kg (the protein section read from the config that already transcribes it, once); `category_split` = the mass the category's own stage computed, at the row named (blood: `compartment_split.tsv`, `as_measured`). composite[a] = Σ f_c × profile_c[a], f_c = P_c / Σ P. A profile that does not sum to 100, a protein basis that does not match the tissue mass's, or fewer than two categories is a `[STOP]`. |
| **L2 Reference male** (D100, D119) | `tissue_mass_g_female`; the split's `female` rows | The base is the male masses. The female masses are recorded beside them and the composite at the female masses is written as a line (`amino_acid_profiles.tsv`, `category_protein_masses.tsv`), never the base. |
| **L3 A split is reported at every row** (D117, D122) | a category's `sensitivity_file` (the D117 table: the profile at every published share × variant) and its split table | For every row of the split, the category's profile at that share, its protein mass at that share, and the composite recomputed with it; `composite_minus_base` against the base of the same sex; the range per amino acid over the male rows. ICRP's whole-blood protein is one of those rows, so the two consistent pairs (profile and mass at one share) stand in one table. |
| **L4 Weighting basis is stated and checked** | the residue-convention table each category names; the PubChem masses; a category's non-protein metabolite pool amounts and the fiber-type mix | Mixing by protein mass assumes one amino-acid mass per gram of protein across categories. The stage computes each category's factor from its residue-convention profile (Σ residue share × free mass / residue mass, C1) and adds a category's non-protein metabolite pool per gram of protein (the amounts table mixed by the fiber-type shares, A10's inputs read back), then writes the composite under amino-acid-mass weights beside the protein-mass composite (`sensitivity_weighting_basis.tsv`; the factors and the largest difference in `composite_summary.ini`). Reported, not used. |
| **L5 Headers name the categories** (D121) | — | Every output's header lists the categories in the composite and each one's protein mass, its source, and its share. The file is `_calculated_amino_acid_standard_composite.tsv`; the word "whole body" is not used until the categories justify it. |

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
| **M4 One food table, no dedup** | `outputs/usda/amino_acids_per_food.tsv` | USDA is the only food table the pipeline scores (D129). Every row names its source and its archive; the same food in two archives is two rows (U6 carried). A food entered by hand is scored by `python run.py score` (M7), which writes nothing into `outputs/match/`. |
| **M5 Every score can be rebuilt by hand** | — | `match_rate_steps.tsv` writes the spreadsheet's own walk per food and reference — the grams (rows 12–20), ΣEAA and "% of Human" (22–23), Step 3 scaled values (30–37), Step 7 percent of reference and its MIN (100–108), Step 10b values (145–152), Total Need to Consume, Percent Wasted, Percent Utilized (153–156) — in the spreadsheet's order and names, with the reference's values in the header (its column A). It is the same number as the minimum ratio (`docs/formula.md`); no second score is written. `match_rate_percent` (row 156 × 100) is carried to column 6, beside the food name and the reference, so the headline is readable without scrolling (D126); the walk itself is untouched. In `match_rate_by_reference.tsv` the difference columns read `<other>_minus_<primary>`, so a positive number means that reference scores the food closer to 100 than the primary does (D125). |
| **M6 One scorer** | — | `match_rate()` in `pipeline/match.py` is the only implementation of the formula. `pipeline/score.py` (M7) and the blend / fortification script (D88) import it; nothing copies it. |
| **M7 The workbook scorer is tooling, not a stage** | `scoring/*.xlsx` | `python run.py score` scores the author's own workbooks and writes `<name>_scored.xlsx` beside each. It is not in `common.COMMAND_ORDER`, writes nothing under `outputs/`, `config/` or `data/`, has no currency test, and nothing downstream reads it — so a proprietary formulation can be scored without entering the repository's record (D129). `scoring/` is gitignored. A workbook whose name ends `_scored` is never taken as an input, and a workbook this tool did not write is never overwritten. |
| **M8 Unreported is imputed at the reference's share; the row says which** | the row's amino acid cells | A scored amino acid the row leaves BLANK is filled at the reference's own share, so its ratio is 1.0 (arithmetically identical to scoring over the amino acids the row does report, both sides renormalised). Every imputed amino acid is named on the row. A `0` is a measurement and scores zero (M3): blank and zero are different statements. This rule is the scorer's alone — `match.py` still lists an incomplete USDA food rather than scoring it (M3), because the pipeline's tables must not mix a nine-amino-acid score with an eight. |
| **M9 The author's rows come back untouched** | the input workbook | The output's first sheet is the input sheet's columns and values exactly as entered, in the author's order, with the calculation columns appended to the right; no column is renamed, reordered or dropped, and the input workbook is never written to. Rows that could not be scored keep their place and carry the reason in a column rather than moving to a separate sheet. |
| **M10 Provenance is a sheet, not a comment** | — | The `reference` sheet carries the generated time, the tool, the reference and where it was read, the hashes of the reference table, the config, the symbol table and the input workbook, and the reference's own values — as rows and columns. No output sheet has comment rows (the header rule, applied to a workbook). |

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
