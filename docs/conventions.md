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
| `config/aggregation_decisions.ini` | Which amino acids are dietary indispensable and how the reference scoring pattern groups them, transcribed by the author from FAO 2013 with table and page per value (D62). The report's three-letter symbols as printed; resolved to one-letter symbols by code from the IUPAC-IUBMB table. Never an abundance, a weight, or a protein. |
| `docs/*.md` | Prose. Every sentence in `docs/methods.md` either carries a citation or describes a computation performed here (`PROVENANCE.md`). |
| `tests/**` | Code and synthetic fixtures. Fixtures contain no real biology. |
| `pipeline/**`, `run.py` | Code. Contains no gene name, accession, term ID, or category. |

## Folder meanings

| Folder | Meaning |
|---|---|
| `config/` | Decisions (hand-written) and the config generated from them. `config/mass_fractions_per_entry.tsv` is the generated Layer B weights table (D61): one row per pool accession, every weight with its interval, provenance in the header — the input the aggregation stage reads. |
| `data/` | What the public databases said, unchanged or minimally tabulated. `data/gene-ontology/` holds the ontology file and the resolved term tables; `data/uniprot_raw/` holds entry JSON as fetched; `data/uniprot_sequences.ini` holds every canonical sequence in tabulated form, with its UniProt MD5, computed MD5, and captured features; `data/iupac/` holds the downloaded IUPAC-IUBMB Table 1 page, its hash, and `amino_acid_symbols.ini` parsed from it; `data/pubchem/` holds the PubChem responses as fetched and `amino_acid_masses.ini` tabulated from them; `data/uniprot_ptmlist/` holds UniProt's PTM vocabulary as fetched and its hash. |
| `data/literature/` | Not committed (`.gitignore`). The literature files as downloaded or hand-obtained, unchanged, one folder per source id; `manifest.ini` (generated) records url, path, sha256, bytes, retrieval time, and whether the file was obtained by hand; `README.md` (generated) is the human-readable provenance. A reviewer re-creates the folder with `fetch_literature.py` plus the hand-obtained files named in the manifest. |
| `outputs/` | Everything computed here that is not config: flags, evidence summaries, delta tables. `outputs/composition/` holds `amino_acid_composition_per_protein.tsv` and the processing and PTM disclosures. `outputs/digest/` holds the in-silico digest tables. `outputs/literature_inventory/` holds the structure of every literature file (members, sheets, header rows) read by code from the hashed files. `outputs/mass_fractions/` holds the ranked weight tables, completeness, cross-checks, and summary (the weights themselves are generated config, D61) — see `docs/pipeline_map.md` for which file answers which question. |
| `logs/` | Debugging only: one timestamped log per stage per run, plus one per test run. Not committed; not part of the record. The record is the header of each generated file plus `outputs/flags.tsv`. |
| `docs/` | `pipeline_map.md` (every stage, its inputs and outputs, and the flowchart), plan, decisions, methods, conventions, handoffs. |

## One command

`python run.py protein-set` runs the protein-set stages in order and then the tests;
`python run.py composition` runs the composition stages in order and then the tests;
`python run.py mass-fractions` runs the literature fetch, the digest, the literature
inventory, and the mass-fraction stages in order and then the tests. The literature fetch
is also the first stage of `protein-set`, because the measured tier (D56) reads the primary
dataset.
`--offline` skips the stages that touch the network and rebuilds from `data/`.
`python run.py excerpt` is tooling, not a stage: it writes bounded, labelled cuts of the large
generated files into `excerpts/` (not committed, not the record) for review in chat.
Nothing in the repository is named by a project-plan step number.

## Rules applied by the code

Each rule names the database field it reads. Where the field cannot settle the case,
the code raises a flag; it never defaults.

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

| Rule | Reads | Result |
|---|---|---|
| **F1 Placeholder** | `file.<n>.url` | A url of `___` is a flag: recorded, skipped, non-zero exit at the end. |
| **F2 Hash pinned** | `file.<n>.sha256` | Blank → the computed hash is written back into config (the only thing the fetcher writes there). Present → the bytes must hash to it; mismatch is a flag and nothing is written. |
| **F3 Unchanged** | served bytes | Bytes on disk == bytes served. Nothing is converted or re-saved. |
| **F4 HTTP failure** | response status | A flag, not a retry loop that hides the failure. |
| **F5 Hand-obtained** | `file.<n>.obtained = manual` | The file must already be on disk under the URL's filename; it is hashed in place and recorded as hand-obtained; never downloaded. |
| **F6 Blocked page** | served body | HTML where a document (pdf/xlsx/zip/rar/docx) was expected is a flag; nothing written. |
| **F7 Stored name** | URL basename | Windows-illegal characters → `_`; the served name is kept in the manifest. |

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
| **B9 Config** (D61) | every measured entry | One generated table, `config/mass_fractions_per_entry.tsv`: per entry and fiber type v, sd, valid, within-tier and combined weights with low / high, match rule, row numbers; provenance once in the header. The stage removes the superseded per-fiber-type `.ini` files and `weights_per_pool_entry.tsv` if present. |
| **B7 Families** | `outputs/digest/shared_pairs.tsv`; composition `free_frac_*` on master rows | A family at the D52 cutoff is a unit (D55); its bound per amino acid = (max − min free-mass fraction across members) × family weight. |
| **B8 Combined** (D58) | every measured entry of every tier | w_i = v_i · MW_i / Σ_all (v_j · MW_j), one denominator per fiber type; the primary standard. Within-tier weights are kept for the checks. |
| **B6 Bands** | `config/carroll_classical_fractionation.ini` `anchor_genes`; `outputs/digest/shared_pairs.tsv` | A gel band is the shared-peptide family of its anchor genes at the D52 cutoff; its weight is Σ w over the family. Families at cutoffs 1 and 2 are both written. |

## Flags

`outputs/flags.tsv` lists every case a rule could not settle, with the stage, subject,
rule, and the field contents that caused it. A flag is closed only by a `[flag.<subject>]`
section in `config/protein_set_decisions.ini` carrying a D-number. Closed flags stay in
the file; `pytest` fails while any flag is open.

## Ini style

`configparser`, one section per entity, keys aligned, `#` comments. Generated files
begin with `# GENERATED by <tool>. DO NOT EDIT BY HAND.` and a provenance header
(ontology file and version, UniProt release, fetch time).

## Naming

Write "Gene Ontology" in full. The two-letter abbreviation appears only inside term
IDs (`GO:0000000`) and in the UniProt query field name `go:`, where it is a literal.
