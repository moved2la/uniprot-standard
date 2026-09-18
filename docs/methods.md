# Methods

Written incrementally, one section per step. Every sentence either cites a retrievable
source or describes a computation performed in this repository (`PROVENANCE.md`).
Numbers below are from the run of 2026-09-11 (ontology retrieved 06:22 UTC, UniProt fetched 06:22–06:36 UTC, config rebuilt 15:10 UTC).

## Protein set definition

### Scope

The primary standard is the myofibrillar protein fraction of human skeletal muscle
(author decision D3, `docs/decisions.md`); the sarcoplasmic fraction is computed
separately as a sensitivity analysis (D4). Citations for the reasoning behind D3 are
attached in the paper step.

### Category

The seed words were resolved to Gene Ontology terms by `pipeline/resolve_ontology_term.py`
as follows. The basic release of the ontology was downloaded from its permanent URL
(D15); the file's version line, retrieval time, and SHA-256 are in
`data/gene-ontology/source.ini` (`data-version: releases/2026-07-26`). The noun form
of each seed word (D13, D14) was matched against term names by exact string equality;
the run stops if the match is not exactly one non-obsolete term. The matches were:

| Tier | Seed word | Lookup name | Term ID | Term definition (verbatim from the ontology file) |
|---|---|---|---|---|
| 1 | myofibrillar (D3) | myofibril (D13) | `GO:0030016` | `The contractile element of skeletal and cardiac muscle; a long, highly organized bundle of actin, myosin, and other proteins that contracts by a sliding filament mechanism. [ISBN:0815316194]` |
| 2 | sarcoplasmic (D4) | sarcoplasm (D14) | `GO:0016528` | `The cytoplasm of a muscle cell; includes the sarcoplasmic reticulum. [ISBN:0198547684]` |

(`data/gene-ontology/resolved_terms.ini`.) The descendant subtree of each term was
enumerated over the `is_a` and `part_of` relations and written to
`data/gene-ontology/tier<N>_subtree.tsv` (18 and 13 terms including the root). The term's ancestors to the root and
its immediate broader and narrower terms are in `tier<N>_ancestors.tsv` and
`tier<N>_neighbours.tsv`.

### Pool

For each subtree term, UniProtKB was queried by `pipeline/enumerate_pool.py` with
`(organism_id:9606) AND (reviewed:true) AND (go:<term>)` at UniProt release
2026_03 (`data/pool_queries.ini`). The pool of a tier is the union of these
per-term results (D20): 242 entries for Tier 1 and 87 for Tier 2 (`data/tier<N>_pool.tsv`). The term-level query alone returned
242 and 87 entries; it equalled the
union for both tiers (`descendant_expansion_matches_union = true` in
`data/pool_queries.ini`), so UniProt's expansion of a term to its descendants is
confirmed by the run rather than assumed. The same query run for each immediate neighbour of the tier term
gave the counts and set differences in `data/tier<N>_pool_sensitivity.tsv`: for Tier 1,
the next broader term (`contractile muscle fiber`, GO:0043292) returns 252 entries, 10
more than the pool; the largest child (`sarcomere`, GO:0030017) returns 221; the child
`skeletal muscle myofibril` (GO:0098723) returns 0 and `cardiac myofibril`
(GO:0097512) returns 7. For Tier 2, the next broader term (`cytoplasm`, GO:0005737)
returns 11,855; the child `sarcoplasmic reticulum` (GO:0016529) returns 77 of the 87.
26 entries belong to both pools (`outputs/pool_overlap.tsv`).

The protein set is therefore defined as: the canonical sequences of human, reviewed
UniProtKB entries annotated to `GO:0030016` or its `is_a`/`part_of` descendants
(Tier 1), and likewise for `GO:0016528` (Tier 2), at Gene Ontology release
`releases/2026-07-26` and UniProt release 2026_03.

### Verification

Each pool entry was fetched by `pipeline/fetch_sequences.py` as entry JSON and FASTA.
The MD5 of the FASTA sequence was compared with the MD5 published in the entry JSON,
the FASTA and JSON sequences were compared character-for-character, lengths were
compared, and only the twenty standard letters were accepted; any failure stops the
run. 303 entries passed (`data/uniprot_sequences.ini`, header). The
entry JSON was saved unchanged to `data/uniprot_raw/`.

### Sequence and processing

Every entry is represented by its canonical sequence (D16). The mature-chain range is
taken from the UniProt feature table, type `Chain`. Exactly one such feature gives the
range (rule R2a; 295 entries); residues outside it are marked `in_master_molecule =
false` in `config/segments.ini` (D19). One entry has no `Chain` feature and uses the
whole sequence (rule R2b; `config/segments.ini`, key `rule`). Seven entries have more than one `Chain` feature;
for these the master molecule is the union of the `Chain` ranges (rule R2d, D24), and
they are listed with their chains in `outputs/multi_chain_entries.tsv`. No `Chain`
feature had a non-exact position, so no flag was raised (`outputs/flags.tsv` is empty).
In total 93 residue runs are marked `in_master_molecule = false` (91 N-terminal, 2
C-terminal) across 92 entries.

### Evidence

The Gene Ontology evidence code of every annotation that placed an entry in its tier's
subtree is recorded in `config/accessions.ini` and summarised in
`outputs/evidence_summary.tsv` (D18). Of the 242 Tier 1 entries, 165 carry at least one
subtree annotation with a code other than IEA (Inferred from Electronic Annotation) and
77 carry only IEA annotations within the subtree; of the 87 Tier 2 entries, 64 and 23
respectively (`outputs/accession_evidence.tsv`). The per-code counts (IBA, IC, IDA,
IEA, IMP, IPI, ISS, NAS, TAS for Tier 1; EXP, IBA, IDA, IEA, IMP, IPI, ISS, TAS for
Tier 2) are in `outputs/evidence_summary.tsv`. Code names are the Gene Ontology
consortium's.

### Disclosed bounds

`pipeline/isoform_processing_deltas.py` computed, from residue counts alone, the
difference in each amino acid's fraction between every non-canonical isoform whose
sequence UniProt provides and the canonical sequence, and between the mature-chain
range and the full canonical sequence. 198 of the 303 entries have at least one such isoform
(662 isoform rows). Across the set, the largest such difference for any single amino acid
is 0.138 for isoform choice (glycine; one short isoform of one entry;
`outputs/isoform_bound.tsv`) and 0.046 for chain processing (serine; one entry with a
long propeptide; `outputs/processing_bound.tsv`), both as fractions of residues. No
decision is attached to these numbers; they bound what a different isoform or
processing choice could do to the composition of any one entry before mass weighting.
The per-isoform and per-entry values are in `outputs/isoform_deltas.tsv` and
`outputs/processing_deltas.tsv`; the two isoform sequences UniProt did not serve
(HTTP 404, recorded in `data/uniprot_sequences.ini`) are absent from the isoform table.

### Limitation

The pool is an annotation set. Its completeness relative to measurement is checked in
Step 3 against the Layer B abundance ranking (D22).

## Composition (Layer A)

*Placeholders in ⟨angle brackets⟩ are filled from the run records after
`python run.py composition` has been run; every other sentence describes a
computation in this repository or carries its citation.*

### Amino acid masses

The masses are the only external constants in Layer A. The correspondence between
one-letter symbol, three-letter symbol, and trivial name for the twenty ribosomally
incorporated amino acids is defined in Table 1 of the IUPAC-IUB Joint Commission on
Biochemical Nomenclature recommendations (*Nomenclature and Symbolism for Amino Acids
and Peptides. Recommendations 1983.* Pure Appl. Chem. 1984, 56(5), 595–624). At every
run, `pipeline/fetch_amino_acid_masses.py` downloads the IUBMB web copy of that table
(`https://iupac.qmul.ac.uk/AminoAcid/tab1.html`, retrieved ⟨date⟩, SHA-256 ⟨sha256⟩;
`data/iupac/source.ini`), parses each row of the form trivial name, three-letter
symbol, one-letter symbol from its text, and stops unless exactly the twenty standard
letters are found once each; the parsed rows are `data/iupac/amino_acid_symbols.ini`
(D26). No name or symbol is typed anywhere in the repository.

Table 1 footnote a states that for the chiral amino acids only the L form is used in
protein biosynthesis, and section 3AA-14.5 that the symbols denote the L configuration.
Each trivial name was therefore sent to PubChem PUG REST prefixed with "L-"
(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<name>/property/MolecularFormula,MolecularWeight,IUPACName/JSON`,
retrieved ⟨date⟩), with the bare trivial name used where PubChem holds no compound
under the prefixed name (⟨n_fallback⟩ of the twenty; the name that resolved is
recorded per row). "Water" and "hydrogen" were sent the same way. A name that resolved
to other than exactly one compound would have stopped the run; none did. The compound
identifier, molecular formula, IUPAC name, and `MolecularWeight` are recorded as served
in `data/pubchem/amino_acid_masses.ini`, with the responses unchanged in
`data/pubchem/raw/`. PubChem serves `MolecularWeight` to ⟨served_decimals⟩ decimal
places, so each free mass is known to ±0.005 g/mol. The in-chain residue mass of each
amino acid is its free mass minus the molecular weight of water (⟨water⟩ g/mol,
CID ⟨cid⟩).

### Residue counts and the two mass conventions

For every entry in `config/accessions.ini`, `pipeline/composition.py` counted the
residues of each of the twenty amino acids in two segment sets: the master molecule
(residues marked `in_master_molecule = true` in `config/segments.ini`; rule R2) and
the full canonical sequence. The counts are the stored ground truth (D7). From them:

- residue-mass vector: count × residue mass; its sum plus one water is the molecular
  weight of the chain;
- free-amino-acid-mass vector: count × free mass; its sum equals the molecular weight
  plus (n − 1) × water, where n is the residue count.

The second identity was evaluated for every row; the largest absolute discrepancy
was ⟨water_identity_max_abs_error_g_per_mol⟩ g/mol (`outputs/composition/composition_summary.ini`).
Each vector is normalised to fractions summing to one. The free-amino-acid
convention is the one used by food composition tables and by laboratory amino acid
analysis, and the one the Match Rate calculation consumes (D7); cysteine is reported
as cysteine at its free mass, and cysteine and methionine are separate columns (D27).
Hydrolysis-sensitive residues are reported at full sequence value (D30). The results
are in `outputs/composition/amino_acid_composition_per_protein.tsv`: ⟨rows_all⟩ rows, one per entry per segment set
(D29).

### Processing disclosure

For each entry, the difference between the free-amino-acid fraction of the master
molecule and of the full canonical sequence was computed per amino acid, together with
the mass removed by processing as a fraction of the full gene product's molecular
weight (`outputs/composition/processing_mass_deltas.tsv`). ⟨entries_with_processing⟩
entries have a master molecule shorter than the canonical sequence;
⟨entries_master_equals_full⟩ do not. The largest single-amino-acid difference in the
set is ⟨largest_processing_delta_free_frac⟩ (⟨amino acid⟩, ⟨accession⟩;
`outputs/composition/processing_mass_bound.tsv`), and the largest fraction of mass
removed from any one entry is ⟨largest_mass_removed_fraction_of_full⟩ (⟨accession⟩).
No decision is attached to these numbers.

### Post-translational modification disclosure

The composition engine counts sequence letters; a residue that carries a
post-translational modification is counted as its unmodified amino acid (D28). To
state the size of what is thereby not modelled, `pipeline/ptm_disclosure.py` read every
feature of type `Modified residue`, `Lipidation`, `Glycosylation`, `Cross-link`, and
`Disulfide bond` from the entry JSON in `data/uniprot_raw/` and priced each site as
follows. The first clause of the feature description (the text before any ";") was
looked up as an identifier in UniProt's controlled vocabulary of post-translational
modifications (`ptmlist.txt`,
`https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt`,
release ⟨release⟩, retrieved ⟨date⟩, SHA-256 ⟨sha256⟩, ⟨records⟩ records of which
⟨records_with_MA⟩ carry an average mass difference). The lookup is exact first; where
no exact identifier exists, the vocabulary's own wildcard identifiers — those
containing "...", such as `Glycyl lysine isopeptide (Lys-Gly) (interchain with G-...)`
— are matched as patterns against the description, and the rule that matched is
recorded per site. A leading "(Microbial infection)" qualifier is set aside for the
lookup and recorded. Where the matched record gives an average mass difference (`MA`),
that is the site's mass delta; the record's correction formula (`CF`) is recorded
beside it, so that for an interchain cross-link the reader can see the delta is that of
the linkage (loss of one water) and not of the partner molecule. An intrachain
disulfide bond has no vocabulary record and was priced as minus the molecular weight of
H₂ from PubChem (⟨h2_mass⟩ g/mol). Every other site — a vocabulary record without a
mass (glycosylation records do not specify the glycan), an interchain disulfide, a
description matching no identifier — was counted and not priced, and is listed by
description in `outputs/composition/ptm_summary.ini`. Sites outside the master
molecule were recorded and excluded from the sums.

Across ⟨entries⟩ entries there are ⟨sites_total⟩ such sites (⟨sites_by_mass_status⟩).
The summed mass delta of the priced sites, as a fraction of the master molecule's
molecular weight, is largest for ⟨accession⟩ at ⟨largest_abs_ptm_mass_fraction_of_mw⟩
(`outputs/composition/ptm_mass_deltas.tsv`); per-site rows are in
`outputs/composition/ptm_sites.tsv`. No decision is attached to these numbers.

## Mass fraction derivation (Layer B)

_Drafted 2026-09-11; revised 2026-09-13 from the run records of `pipeline/mass_fractions.py` (run 2026-09-13T04:35Z, UniProt release 2026_03, Gene Ontology release 2026-07-26). Every sentence either carries a citation (sources as in `config/mass_fraction_decisions.ini`, referred to by author and year, D51) or names a computation in this repository and the file it wrote._

### Inputs

Layer B reads: the protein set and mature-chain molecular weights from the protein-set and composition stages (`config/accessions.ini`; `outputs/composition/amino_acid_composition_per_protein.tsv`, column `mw`, rows `segment_set = master`); the literature files fetched and hashed by `pipeline/fetch_literature.py` (`data/literature/manifest.ini`); the hand-written `config/mass_fraction_decisions.ini` (source citations, roles, column maps, digest rules — by test, no protein and no abundance number); and the hand-transcribed `config/carroll_classical_fractionation.ini` (D50).

### The protein set as weighted

Tier 1 (contractile) is the 242-entry Gene Ontology `myofibril` set of the protein-set stage. Tier 2 (support) is the measured remainder (D56): every gene name in the primary dataset that is not in Tier 1, mapped by `gene_exact` to exactly one reviewed human UniProt entry, excluding entries in MaxQuant's contaminant list (rules M1–M5). Run record (`data/pool_queries.ini`, `data/tier2_gene_mapping.tsv`): 3,925 dataset genes → 3,564 members, 183 already in Tier 1, 98 ambiguous (several reviewed entries for one symbol), 77 unmapped, 3 contaminants. Six entries whose canonical sequence contains selenocysteine were excluded under rule R5 (D59; `outputs/excluded_non_standard_alphabet.tsv`); 54 entries with a non-exact Chain position were resolved under rule R2e (D60; `outputs/chain_positions_resolved.tsv`). The weighted set is 3,788 entries: 242 contractile, 3,546 support (`config/accessions.ini`).

### Primary dataset

The primary dataset is Murgia et al. 2021 (role by D32). Its quantification is described in its Methods: "the summed intensity of the peptides of each protein was divided by the number of theoretically observable peptides [8], and the resulting values were normalized to the expression of α-skeletal actin (ACTA1)" (Murgia 2021, Methods, "Single-fiber proteomics"); reference [8] there is Cox et al. 2014. Dataset 1 (`data/literature/murgia_2021/13395_2021_279_MOESM1_ESM.xlsx`, SHA-256 `223f8753…10d2a`) gives per protein and fiber type a median value, a standard deviation, and the number of fibers in which the protein was quantified, for pure type 1, 2A, and 2X fibers (≥ 80 % of one MYH by unique-peptide intensity; 19 / 29 / 13 fibers). As read by `pipeline/mass_fractions.py`: 3,857 data rows; 3,607 single gene names of which 84 occur on more than one row; 139 cells listing several genes; 491 median cells `NaN`, every one with a valid-value count of 0; ACTA1 = 100 000 in every median column (`outputs/mass_fractions/mass_fractions_summary.ini`, `[dataset]`).

### Quantity used

Intensity ÷ theoretical peptide count is iBAQ (Schwanhäusser et al. 2011). Two weightings are computed (`pipeline/mass_fractions.py`):

    within tier (D31):   w_i = v_i · MW_i / Σ_{j ∈ T} ( v_j · MW_j )         Σ_T w = 1
    combined  (D58):     w_i = v_i · MW_i / Σ_{j ∈ all} ( v_j · MW_j )       Σ_all w = 1

where v_i is the median value from the mapped column and MW_i the mature-chain molecular weight. The ACTA1 normalisation in v_i cancels in both ratios. The combined weighting is the primary standard (D58); the within-tier weighting is kept for the family bounds, the gel comparison, and the tier shares. Per-row intervals `w_low`, `w_high` use (v ∓ SD) over the unchanged denominator, clipped at zero (B5, D45, D53a); a shared row adds its whole value to each member's `w_high` (D48). Every weight, with its interval, is written per fiber type to `config/mass_fractions/<type>.ini` (one section per entry, the citation on the same lines) and, in one table with every fiber type side by side, to `config/mass_fractions_per_entry.tsv` (D61), the dataset file, hash, sheet, column names, and retrieval time in its header. The SD is the dataset's spread between fibers of one type, so the interval is that of a random pure fiber of the type, not of the median (D45).

Run record (`[combined]`): the contractile tier holds 71.7 % / 72.8 % / 71.6 % of the combined standard in types I / IIa / IIx; 3,575 / 3,683 / 3,460 entries carry weight above zero; the ten largest entries hold 54.8 % / 52.9 % / 50.3 %. The largest entry in every fiber type is titin (20.1 % of the combined standard in type I; `outputs/mass_fractions/combined_entries_ranked.tsv`). 55 pool entries have no Dataset 1 row and carry w = 0 (`pool_entries_without_dataset_row.tsv`); 121 Dataset 1 genes match no pool entry and are listed with their mapping outcome (`dataset_rows_outside_pool.tsv`); the six R5-excluded entries would together have held less than 0.01 % of the combined standard (`excluded_entries_mass_share.tsv`).

Why not LFQ: MaxLFQ is defined for comparing a protein across samples, and Murgia 2021 used it for its statistics (Murgia 2021, Methods, citing Cox et al. 2014). Momenzadeh et al. 2023 publish MYH isoform fractions from the same 53 fibers under both quantities (Table S2, iBAQ; Table S3, MaxLFQ). As read from the downloaded files (`outputs/mass_fractions/momenzadeh_2023_myh_fractions_ibaq_vs_lfq.tsv`), the mean MYH4 fraction is 0.105 under Table S2 and 0.0037 under Table S3 — a factor of 28 between the two quantities for one isoform. The paper's text, as summarised by later reviews, states the iBAQ value as the smaller; the published tables give the reverse, and the article mentions an iBAQ table "prior to removal" of a peptide shared with MYH2. The repository states what the files contain; which table matches the text is resolved at the paper step from the article itself. Either way, the two quantities differ by more than an order of magnitude for one protein from one set of fibers, which is the point. Datasets publishing only LFQ-type values (Moreno-Justicia 2025, DIA-NN; Deshmukh 2021, MaxLFQ per its Methods) are therefore used for within-protein ratios between fiber types and for an intensity-share comparison, never for weights.

Why not TPA: the total-protein approach (Wiśniewski 2017) reads mass fraction as raw intensity ÷ summed intensity. Murgia 2021 publishes iBAQ-type values, not raw intensities, so TPA is not computable from the file. iBAQ × MW and TPA differ per protein by N_i / MW_i (theoretical peptides per unit mass), computed from sequence by `pipeline/digest.py` (`outputs/digest/density_ranked.tsv`).

### Identity

Dataset 1 identifies proteins by gene name only. A row joins a pool entry by exact gene symbol against `config/accessions.ini` (rule B1), or, where the symbol differs from the entry's, through the gene-mapping table written by the measured-tier enumeration — a symbol UniProt resolves to a pool entry, primary or synonym (rule B1b, D57): 202 rows joined that way in the run. Several rows for one gene are summed and the minor rows' share written (D47; 83 entries). A cell listing several genes is split (D48): 138 of the 139 such rows touch the pool; the signal arriving through them is 0.011 % / 0.017 % / 0.022 % of the Tier 1 denominators (`[weights]`). Every match rule is written per entry (`weights_per_pool_entry.tsv`, `config/mass_fractions_per_entry.tsv`, `shared_gene_rows.tsv`).

### In-silico digest (Layer A computation used by Layer B)

`pipeline/digest.py` digests every canonical sequence of the set under the rules in `[digest]` (trypsin after K/R, 6–30 residues, 0 missed cleavages: Schwanhäusser et al. 2011 as restated in Nature Sci Data 2018;5:180128; cleavage before proline computed both ways, D49) and writes per entry the theoretical peptide count, `peptides_per_kDa`, and the peptides shared with other entries (`outputs/digest/`). Density is ranked and disclosed; no threshold is applied (D43). Shared-peptide pairs define families at a cutoff of ≥ 2 shared peptides (D52, D55).

### Fiber type

Fiber types are Murgia 2021's own (D41, B3). Moreno-Justicia et al. 2025 report that pure type 2X fibers "could only be detected at the RNA but not at the protein level" (Moreno-Justicia 2025, Results); the IIx column of this standard rests on 13 fibers and whether the standard reports I / IIa / IIx or slow / fast is decided at aggregation (D46).

Within the myosin heavy chain family, Dataset 1's gene-level values are not a clean split by isoform: in fibers selected as ≥ 80 % MYH2 or MYH1, MYH7 carries about half its type-I value, and MYH1 in 2X fibers sits below MYH2 and MYH7 (`band_families.tsv`); MYH1 and MYH2 share 80 in-window tryptic peptides (`outputs/digest/shared_pairs.tsv`). Under D55 a family is carried as a unit and the spread of each amino acid's free-mass fraction across its members, times the family's weight, is written as the family bound (`family_bounds.tsv`). Run record: the largest Tier 1 family bound is 0.020 (lysine; the titin family, which at the cutoff includes hemicentin-2) and 0.015 (arginine; the myosin heavy chain family); in Tier 2 the largest is 0.002.

### Cross-checks and their questions

Each cross-check is admitted only with a named question (D42) and reported as the numbers each method gives, with their quotient, without attribution of any difference to either method (D54).

| Source | Question it answers | Output | Run record |
|--------|---------------------|--------|------------|
| Moreno-Justicia et al. 2025, Supplementary Dataset 10 (pseudobulked log2 abundance per participant, slow and fast) | Are Murgia 2021's slow/fast ratios (61 fibers, 4 donors) reproduced at 1,038 fibers, 5 donors, different instrument and search? | `ratio_check_moreno-justicia_2025.tsv` | 2,983 rows, 2,303 pool entries compared; median absolute log2 disagreement between Murgia's I/IIa and their slow/fast, per protein: 0.58 |
| Deshmukh et al. 2021, Supplementary Data 3 (median LFQ intensity per protein, type 1 and type 2 pools) | Same question, pooled fibers, third workflow; joined by accession | `ratio_check_deshmukh_2021.tsv` | 3,864 / 3,786 rows; 2,870 / 2,842 matched a pool entry |
| Carroll, Carrithers & Trappe 2004 | Does iBAQ × MW reproduce the MHC : actin ratio measured by SDS-PAGE against standards in human single fibers, types I and IIa? And, since D58: does the band's share of the combined standard reproduce the gel's fraction of total fiber protein? | `classical_check_carroll_2004.tsv` | below |

MHC : actin (D52 families, cutoff 2): gel 1.56 (I) and 2.21 (IIa); iBAQ × MW 2.55 and 3.59 (quotients 1.64, 1.63); Deshmukh 2021 intensity share 3.83 (slow) and 3.91 (fast; fast pools ~80 % IIa / 20 % IIx per its Fig. 1f), quotients 2.46 and 1.77. At cutoff 1 the myosin family absorbs single-peptide bridges across the 3,788-entry set and the iBAQ × MW quotient rises to 3.9 / 3.5 — the disclosure that motivates the cutoff. Absolute shares are written by the same stage as rows `share_of_combined_<band>` (the band's share of the combined standard against the gel's µg per mg ÷ 1000). No usable human weighed-mass source exists for titin or nebulin (Trappe et al. 2002 examined and rejected: uncalibrated densitometry); the giant-protein question rests on the disclosed density rank (D39, D43).

### Layer B uncertainty terms (inputs to aggregation)

1. Between-fiber spread within a fiber type: Dataset 1's `Standard deviation` per fiber type, carried per row as `w_low` / `w_high` (D45, D53a).
2. Shared-peptide family bounds, per family and amino acid (D55; `family_bounds.tsv`).
3. Between-dataset ratio disagreement per protein (Moreno-Justicia 2025; Deshmukh 2021).
4. iBAQ × MW versus the gel (Carroll 2004), for MHC and actin, types I and IIa: reported; the aggregation step builds the profile under the weights as they are and under weights rescaled to the gel's MHC : actin and reports the per-amino-acid difference (D54b).
5. iBAQ × MW versus TPA — per-protein factor `peptides_per_kDa`; applied only if raw intensities become available.

### Carried checks (from the protein-set and composition stages)

All computed by `pipeline/mass_fractions.py` and written to `outputs/mass_fractions/` (`weighted_bounds.tsv`, `branch_exclusive_mass.tsv`, summary sections). Run record, Tier 1: largest weighted isoform bound 0.0055 (glutamine, titin), largest weighted processing bound 0.0002 (methionine, ACTA1), largest weighted PTM mass bound 0.0014 (TPM3); Σ w over glycosylated entries 2.1 % / 2.2 % / 3.0 % of Tier 1 (I / IIa / IIx) and 17 % / 17 % / 19 % of Tier 2. The 26 entries that were in both ontology tiers are now Tier 1 only, since Tier 2 is defined as the remainder. Every Gene Ontology branch of the Tier 1 subtree, including the cardiac-myofibril branch, carries its returned and branch-exclusive mass in `branch_exclusive_mass.tsv`.

## Aggregation (the standard)

*Placeholders in ⟨angle brackets⟩ are filled from the run records after `python run.py standard`.*

### Inputs

The per-entry weights (`config/mass_fractions_per_entry.tsv`, D61), the residue counts of every entry's mature chain (`outputs/composition/amino_acid_composition_per_protein.tsv`, `segment_set = master`), and the fetched amino acid masses (`data/pubchem/amino_acid_masses.ini`). Every input's hash is in every output header (rule A0).

### Profile

The standard is the molar mixture of its proteins (rule A1, D64). For a set of entries with mass fractions w_i and molecular weights MW_i, the moles of residue a per gram of protein are

    n_a = Σ_i (w_i / MW_i) · count_{i,a}

and the profile in either convention is p_a = n_a · m_a / Σ_b n_b · m_b, with m_a the in-chain residue mass (residue convention) or the free amino acid mass (free convention, D7). Three sets are aggregated (rule A2): the total, every weighted entry under one denominator (`w_combined`, D58); the contractile proteins; and the builders — the two tiers under their within-tier weights (D31); tiers are named, never numbered (D66). The deliverable, `outputs/standard/_calculated_amino_acid_standard.tsv`, gives the three sets per fiber type in percent of amino acid mass and a final column, the fiber-type mix of the totals from `config/fiber_type_mix.ini`; the calculation is written out in `docs/formula.md`. The fraction mixture Σ_i w_i f_{i,a}, renormalised, is written beside the profile as a check; the largest difference between the two is ⟨max_abs_a1_minus_a1c⟩ (`standard_summary.ini`). The wide tables give g per 100 g protein without renormalisation: the residue column sums to ⟨residue_g_sum⟩ (100 less the terminal water of each chain), the free column to ⟨free_g_sum⟩.

Run record: the contractile proteins hold ⟨share_contractile⟩ of the total; the three largest fractions of the total fiber-type-I profile (free convention) are ⟨three_largest_I⟩. The largest difference between the total and the contractile profile is ⟨total_minus_contractile_max⟩, between the contractile proteins and the builders ⟨contractile_minus_builders_max⟩ (`profile_differences.tsv`).

### Fiber types

The standard is reported for types I, IIa, and IIx as the primary dataset defines them (rule B3; D46). Because a mix of fiber types is a convex combination of the pure types, the range of each amino acid over the three pure profiles bounds every mix (rule A5); that bracket is at most ⟨bracket_max⟩ in the total and the IIa − IIx difference at most ⟨IIa_minus_IIx_max⟩ (`profile_differences.tsv`). No fiber-type mix is applied.

### Uncertainty

Layer B's between-fiber term is propagated in log space (rule A3, D63). For each entry and fiber type a log-normal is fitted exactly to the published median m and standard deviation s (σ² = ln x, x = (1 + √(1 + 4 s²/m²)) / 2; μ = ln m). Two quantities are sampled, ⟨draws⟩ draws each with seed ⟨seed⟩, the whole profile recomputed per draw: the *between-fiber spread*, v ~ LogNormal(μ, σ), what a random pure fiber of the type looks like (the D45 reading), and the *median uncertainty*, log v ~ Normal(μ, σ/√n) with n the entry's valid-value count, the sampling uncertainty of the fiber type's median profile — the interval stated beside the standard. The one assumption, that an entry's between-fiber distribution is log-normal, is the field's standard model for intensities and is stated in every header. On the linear scale the published SD exceeds the median for ⟨sd_ge_median_I⟩ of ⟨entries_with_value_I⟩ type-I entries with a value (`sd_to_median_ratio.tsv`), which is what a right-skewed spread looks like when summarised as an SD; the anchor protein's SD is zero by construction and its variation lives in every other entry's.

Run record: the 95 % median-uncertainty interval of the total fiber-type-I profile has a half-width of at most ⟨median_hw_max_I⟩ (free convention; `uncertainty_intervals.tsv`); the between-fiber spread at most ⟨spread_hw_max_I⟩.

The other terms are bounds, not sampled (`bounds_after_weighting.tsv`): the isoform, processing, and PTM bounds of each entry times its weight (D16, D28), the shared-peptide family bounds scaled to the profile (D55), and the D48 shared-row excess; each is written as the largest single contribution and as the worst-case sum. The MHC : actin sensitivity (rule A6, D54b) rescales the myosin heavy chain family with actin fixed, and the actin family with myosin fixed, to each ratio the three methods measured (`classical_check_carroll_2004.tsv`), for types I and IIa, and reports the profile's range; no method is a reference, and its largest shift of any fraction is ⟨mhc_actin_max_shift⟩ (`sensitivity_mhc_actin_spread.tsv`). The completeness bound (rule A7) converts the outside genes' molar share to a mass share under the stated assumption that they have the pool's molar-mean molecular weight, at most ⟨completeness_bound⟩ of the total. Every term is placed on one scale — the maximum absolute shift of the fraction from the standard — and the widest term named per amino acid (`uncertainty_per_amino_acid.tsv`); this scale is provisional.

### Stress test

Beyond the uncertainty of the standard from its source, `pipeline/stress.py` (rule A9, D65) reports how far the standard moves when the weights are pushed by stated amounts: the composition distance among the ⟨k⟩ largest entries and against the standard; the profile from the top-k entries only; the top-k entries and each band family removed; every entry zeroed and doubled (`stress_influence_per_entry.tsv`); the contractile share forced to ⟨tier1_shares⟩; every weight multiplied by MWᵅ for α in ⟨alphas⟩; the TPA weighting v × N_peptides, computed exactly from the digest in place of v × MW; and every weight multiplied by an independent log-normal factor of ×/÷ ⟨factors⟩ over ⟨random_draws⟩ draws. Run record: the largest shift of any fraction under each scenario is in `stress_summary.ini`; ⟨stress_largest_line⟩.

### Indispensable amino acids

The indispensable amino acids and the groupings of the reference scoring pattern are those of FAO 2013 (Food and Nutrition Paper 92), transcribed by the author with table and page into `config/fao_2013_indispensable_amino_acids.ini` (D62) and resolved to one-letter symbols by code through the IUPAC-IUBMB table. `eaa_subset.tsv` reports each heading per profile and fiber type as mg of free amino acid per gram of protein and as its share of the free-convention profile.

### Non-protein metabolite pools

**Definition.** Skeletal muscle holds amino acids in three states, and this standard treats each differently.

1. **Protein-bound** — amino acids in peptide chains. This is what the protein-only standard counts, and it is what
   "bound" means everywhere else in the amino acid literature, where free amino acids are defined as those *not*
   bound in protein.
2. **Free** — single amino acids in solution in the cell, turned over continuously. Excluded by design (D68):
   transient, heavily regulated, and not a standing demand.
3. **Non-protein metabolite pools** — amino acids held in stable low-molecular-weight compounds that are not
   proteins: carnosine (histidine), and, for later non-indispensable extensions, creatine (glycine) and glutathione
   (cysteine, glycine, glutamate). Replaced on a timescale of weeks rather than minutes. The body had to eat these
   amino acids to put them there, so they are a real dietary demand that a protein-only count misses.

This third state is what the adjusted standard adds. The name is this repository's, not the field's: the literature
has "non-protein nitrogen-containing compound" (which also covers the free pool and urea) and
"histidine-containing dipeptides" (carnosine's family only), but no term for exactly this class, because the class is
rarely grouped. "Bound pool" was avoided because "bound" already means protein-bound in this literature and would
point a reader at state 1. Decisions D68, D73 and D75 were written before the term was settled and call these
"bound pools"; they mean state 3 (D82).

Muscle holds histidine in one stable non-protein form, the dipeptide carnosine (β-alanyl-L-histidine); in human muscle carnosine is the only histidine-containing dipeptide (Harris et al. 2012, p. 6, "Species"), so it is the only bound pool carried (D75). The free amino acid pool is excluded by design (D68). The adjustment is an inventory (D73): the histidine held as carnosine is added to the histidine held as protein in the same kilogram of muscle (rule A10, `pipeline/non_protein_metabolite_pools.py`), and the adjusted standard is written beside the protein-only one, which is unchanged (D76). Per fiber type the pool is the single-fiber carnosine content of myosin-ATPase-typed fibres of human vastus lateralis from four men (Harris, Dunnett & Greenhaff 1998, abstract, p. 639): type I 10.5 and type II 23.2 mmol per kg dry muscle, the type II value applied to fiber types IIa and IIx alike (D74). The protein content is the direct chemical measurement of Mingrone et al. 2001 in 16 normal-weight subjects, 155 g protein and 763 g water per kg wet muscle (Table 1, p. E368), moved to the dry basis of the pool through that water content: 654 g protein per kg dry muscle (D77); the muscle measured there is rectus abdominis, stated in every header. One mole of carnosine yields one mole of histidine (free mass 155.15 g/mol, PubChem); the β-alanine is not one of the twenty and is not counted. The histidine from the pool amounts to 0.097, 0.211, and 0.213 of the protein-bound histidine in fiber types I, IIa, and IIx (`non_protein_metabolite_pool_amounts_per_kg_muscle.tsv`: 1.63 and 3.60 g per kg dry muscle from the pool against 16.7–17.1 g from protein), and histidine's share of the total profile moves from 2.206 to 2.416 % in fiber type I, 2.250 to 2.712 % in IIa, 2.229 to 2.691 % in IIx, and from 2.227 to 2.568 % in the standard under the fiber-type mix (`non_protein_metabolite_pool_adjustment_per_amino_acid.tsv`); every other amino acid moves down by the same factor, 0.35 % of its value. Before any of this is computed the stage checks that every cited source has a hashed file in the literature manifest and stops otherwise (D78). The indispensable amino acids of the adjusted profile are written beside the protein-only ones (`eaa_subset_with_non_protein_metabolite_pools.tsv`).

One sensitivity is reported and it is not a correction. The single-fiber spread, carried through as a log-normal fitted to the published mean and ± value (the D63 pattern; the abstract does not say whether its ± is an SD or a standard error, and the file says so), puts the 2.5 and 97.5 percentiles at 2.4–30.4 mmol per kg dry muscle for type I fibres and 4.9–69.8 for type II, giving histidine 2.25–2.81 % in fiber type I and 2.35–3.63 % in IIa (`sensitivity_non_protein_metabolite_pool_spread.tsv`) — a between-fiber spread from four subjects, not the uncertainty of the mean in the table. The resupply (turnover-weighted) frame is not computed: carnosine's washout half-life is cited (4.6 weeks, 95 % CI 3.2–7.0, Yamaguchi et al. 2021; 5.8 weeks as the first-order figure of Baguet et al. 2009, whose observed washout was linear), but the muscle-protein turnover rate it would be divided by was not sourced for this step, and the published range of that rate spans roughly four-fold, which is why the inventory frame is the standard (D73). No whole-muscle cross-check was made: the in vivo MRS measurement found (Vega et al. 2022) is paywalled and was not used, and the second single-fiber source (Tallon et al. 2007) is gated with no values in its abstract; the per-fiber-type columns therefore rest on one study of four men, a stated limitation. The single-fiber source is male and the fiber-type mix is female-only (`fiber_type_mix.ini`); both are stated in the table header.

Values from these papers were read from the PDFs on disk by the assistant and audited by the author, each with its page and table or paragraph beside it in `config/non_protein_metabolite_pools.ini` (D79); the calculation contains no model.

## The food tables (USDA FoodData Central)

Match Rate compares a food's amino acid profile with the standard. The food side is USDA FoodData Central, from which two data types carry measured amino acid panels: SR Legacy (April 2018, the final release of Standard Reference, frozen) and Foundation Foods (April 2026 release, analytical values with per-value sample counts and ranges). Branded Foods carries manufacturer label data and no amino acids, and the Food and Nutrient Database for Dietary Studies compiles its values from the other types; neither is used (D83). The CSV archives are downloaded by hand from the FoodData Central download page and placed unrenamed in `data/usda/`, where `pipeline/usda.py` reads them (rules U1–U7). The archives themselves are not committed — this repository is not a mirror of USDA's data — but `data/usda/manifest.ini` records each file's name, size and SHA-256, so that the same release downloaded by anyone else can be checked against the run.

FoodData Central publishes its data as a long table: one row per food per nutrient, per 100 g of food, with every nutrient class in the same file and distinguished only by a nutrient id. The stage pivots that table into one row per food (`outputs/usda/amino_acids_per_food.tsv`) and does nothing else to it. Each archive is identified by reading its own `food.csv` rather than by its file name (U1). A food counts when the archive's own membership file lists it (U2): FoodData Central ships superseded earlier versions of a food alongside the current one — in the April 2026 Foundation archive, 469 rows of `food.csv` carry the Foundation data type while `foundation_food.csv` lists 395 — and the membership file is the archive's own statement of which is current. A nutrient is amino acid X when its name in the archive's `nutrient.csv` equals the IUPAC-IUBMB trivial name of X already parsed in the composition step (U3), so no amino acid is named by hand anywhere in the ingestion; the resulting map is written out as `usda_nutrient_map.tsv`. Values are carried as published, with the archive's own minimum, maximum, median and sample count where it gives them, and nothing is converted or filled in; a nutrient reported in a unit other than the expected one is flagged and not carried (U4).

Two properties of the source are stated rather than resolved. USDA reports no asparagine and no glutamine — acid hydrolysis converts them to aspartic and glutamic acid — so no food can carry all twenty, and each row states how many it carries and which letters it lacks rather than being judged complete or incomplete here (U5); which amino acids a score requires belongs to the scoring step. And USDA reports the sulphur amino acid as "Cysteine" for some foods and as "Cystine" for others, never both: the cysteine values join as amino acid C by the name rule, the cystine values are carried in their own column, and whether a cystine figure converts to cysteine is left open, since it turns on whether USDA's figure is already expressed as cysteine equivalent. Neither affects the indispensable amino acids.

## Match Rate

Match Rate scores a food's amino acid proportion against a reference proportion. The formula is the one the NutriMatch spreadsheet has used, transcribed as it is (D87) and reduced to its arithmetic: over a scored set of amino acids, both the food and the reference are expressed as shares of their scored total; each amino acid's ratio is food share over reference share; the score is the smallest ratio, and the amino acid that gives it is the limiting amino acid (`pipeline/match.py`, rule M1; the calculation by hand in `docs/formula.md`). Because both sides are shares, the score is a property of the proportion alone: the unit a food is reported in and its protein content do not enter it. The spreadsheet reaches the same number by scaling the food to the reference's total, dividing every amino acid by the limiting one, and calling the excess above the reference "wasted"; that arithmetic reduces to one minus the smallest ratio exactly, so the per-amino-acid ratios are written beside the score as the record of where the surplus lies, not as a second score (M5). Four foods from the spreadsheet are carried as test fixtures, and the stage must reproduce their spreadsheet scores to the last digit before any other test is considered (`tests/test_match.py`).

The scored set is FAO 2013's nine indispensable headings (D62), each group reduced to its indispensable member — methionine for the sulphur pair, phenylalanine for the aromatic pair — so that cysteine and tyrosine are not scored (D86). The reduction is the author's decision, recorded in `config/match_rate.ini` as the three-letter symbol of the member scored and resolved through the IUPAC-IUBMB table; no amino acid is otherwise named in the configuration or the code (M2). A conditionally indispensable set is a separate scored set for a later output.

Three references are scored, each carrying its own scored set and a label that every output row repeats (M2). The primary reference is the skeletal muscle standard with its non-protein metabolite pools (`standard_with_non_protein_metabolite_pools`, D82), in the free-amino-acid convention that USDA reports foods in (D7); the protein-only column is scored beside it so the effect of the pools on scores is visible. The third is the original reference — Gorissen et al. 2018's human muscle column (D84), scored on the eight amino acids that paper sums as essential, because acid hydrolysis destroys tryptophan; it is carried so that old and new scores stand side by side, and it is not a standard of this repository. A food scored on eight amino acids and the same food scored on nine answer different questions, and the difference between them is reported, not explained.

Foods come from two tables (M4): the USDA rows of `outputs/usda/amino_acids_per_food.tsv`, as published, and `config/food_amino_acids_other_sources.csv`, a hand-maintained file for foods USDA does not carry — a supplier's analysis, a label, a paper's table — each row citing its source. No food is dropped or preferred across the two, or across the two USDA archives (U6 carried): the same food from two sources is two rows. A food that does not report a scored amino acid is listed in `foods_not_scored.tsv` rather than scored (U5 carried, M3); a published zero is scored as published — the ratio is zero, the score is zero, and the limiting amino acid is named — since the stage does not decide what USDA measured. Whether a food with almost no protein should be scored at all is a display question the score does not answer: the protein content is carried in every row for the reader to filter on, and the score itself is the proportion.

The outputs are `match_rate_per_food.tsv` (one row per food and reference: the score, the limiting amino acid, the scored amino acids' sum, that sum as a percent of the food's protein, and the per-amino-acid ratio), `match_rate_by_reference.tsv` (one row per food with every reference's score and limiting amino acid side by side and the primary's difference from each other in percentage points), `match_rate_steps.tsv` (one row per food and reference carrying the spreadsheet's own walk — the scaled values, the per-amino-acid ratios and their minimum, the Step 10b values, the total to consume, the percent wasted and utilized — so that any published score can be rebuilt in the spreadsheet from its row, M5), one ranked table per reference, `limiting_amino_acid_counts.tsv` (how many foods each amino acid limits under each reference — the table that says where scores move between references), `foods_not_scored.tsv`, and `match_summary.ini`. Every header carries the Match Rate version label from `config/match_rate.ini`, the hash of every input, and each reference's values. The stage is the public single-food scorer; a blend and fortification script — the inverse problem, which free amino acids lift a fixed blend's limiting ratio — is a separate script that calls the same function, so the arithmetic lives in one place (D88).
