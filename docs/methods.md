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
