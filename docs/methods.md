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

_Drafted 2026-09-11, revised 2026-09-12 after Dataset 1 was read. Every sentence either carries a citation (sources as in `config/mass_fraction_decisions.ini`, referred to by author and year, D51) or names a computation in this repository. Angle-bracket placeholders ⟨…⟩ are filled from run records at the mass-fraction run._

### Inputs

Layer B reads three things: the protein set and full-product molecular weights from the protein-set and composition stages (`config/accessions.ini`; `outputs/composition/amino_acid_composition_per_protein.tsv`, column `mw`); the literature files fetched and hashed by `pipeline/fetch_literature.py` (`data/literature/manifest.ini`); and the hand-written `config/mass_fraction_decisions.ini`, which holds source citations, role decisions, column-to-fiber-type mappings, and digest rules — and, by test, no protein and no abundance number.

### Primary dataset

The primary dataset is Murgia et al. 2021 (role by D32). Its quantification is described in its Methods: "the summed intensity of the peptides of each protein was divided by the number of theoretically observable peptides [8], and the resulting values were normalized to the expression of α-skeletal actin (ACTA1)" (Murgia 2021, Methods, "Single-fiber proteomics"). The reference [8] in that sentence is Cox et al. 2014. Dataset 1 reports, per protein and fiber type, "the median expression values and the percent of the maximal value" and "the number of fibers in which each protein was identified (valid values)"; "proteins are indicated by gene name" (Murgia 2021, Results, first paragraph). Fibers were "relatively pure myofibers containing at least 80% of either MYH7 or MYH2 or MYH1. A total of 61 fibers (19 type 1, 29 type 2A, and 13 type 2X fibers)" from "four younger (aged 22–27 years)" donors; "MYH expression was quantified by the intensities of peptides unique for each isoform" (Murgia 2021, Methods). Raw data are PXD006182 (Murgia 2021, Availability of data), shared with Murgia et al. 2017, which is therefore not an independent source (D33).

### Quantity used

The quantity intensity ÷ theoretical peptide count is iBAQ (Schwanhäusser et al. 2011). Mass fraction of protein *i* within tier *T* is computed by `pipeline/mass_fractions.py` as

    w_i = v_i · MW_i / Σ_{j ∈ T} ( v_j · MW_j )

where v_i is the median value read from the column mapped in `[columns.murgia_2021.file.1]` and MW_i is the composition-stage molecular weight of the mature chain (`segment_set = master`). The ACTA1 normalisation in v_i is a per-fiber constant and cancels in the ratio. Sums are taken within tier: Σ_{Tier 1} w = 1 and Σ_{Tier 2} w = 1 (D31). A pool entry with no row in Dataset 1 receives w = 0 and is listed (`outputs/mass_fractions/unmatched_pool_entries.tsv`); a Dataset 1 row with no pool entry is listed with its v·MW share of the whole file (`outputs/mass_fractions/dataset_rows_outside_pool.tsv`) — this is the pool-completeness check carried from the protein set.

Why not LFQ: MaxLFQ is defined for comparing a protein across samples, and Murgia 2021 used it for its statistics (Murgia 2021, Methods, citing Cox et al. 2014); Momenzadeh et al. 2023 reports MYH isoform fractions computed both by iBAQ (its Table S2) and by MaxLFQ (its Table S3) from the same fibers, and the two disagree by ⟨factor, from the downloaded tables⟩ for MYH4. Datasets that publish only LFQ-type values (Moreno-Justicia 2025, DIA-NN protein-group LFQ; Deshmukh 2021, MaxLFQ per its Methods) are therefore used for within-protein ratios between fiber types, never for cross-protein weights (D34, D36).

Why not TPA: the total-protein approach (Wiśniewski 2017) reads mass fraction as raw intensity ÷ summed intensity. Murgia 2021 publishes iBAQ-type values, not raw intensities, so TPA is not computable from the published file. iBAQ × MW and TPA differ per protein exactly by the factor N_i / MW_i (theoretical peptides per unit mass); that factor is computed from sequence alone by `pipeline/digest.py` and reported as `peptides_per_kDa` (`outputs/digest/theoretical_peptides.tsv`).

### Identity

Dataset 1 identifies proteins by gene name only (`Gene names`, sheet `Sheet1`, header row 3; no accession column). A row is joined to a pool entry by exact gene symbol against the `gene` column of `config/accessions.ini` (rule B1, `pipeline/mass_fractions.py`). Where one gene name occurs on several rows — separate MaxQuant protein groups of one gene product; 85 such names — the rows are summed and the share contributed by the minor rows is written per gene (D47). Where one cell lists several genes separated by `;` — 139 rows, 3.0 % of summed medians across the file — the value is assigned to the single pool member if only one is in the pool, or divided among pool members in proportion to their own single-name rows with the whole row carried as each member's upper bound (D48); the signal share arriving this way is printed per tier and fiber type. A pool gene with no row receives w = 0 and is listed; a Dataset 1 gene with no pool entry is listed with its v·MW share (the pool-completeness check). Dataset 1 was searched against "the human UniProt FASTA reference proteomes version of 2016" (Murgia 2021, Methods); gene symbols renamed since then surface in the unmatched lists. ⟨Counts from the run.⟩

### In-silico digest (Layer A computation used by Layer B)

`pipeline/digest.py` digests every canonical sequence under the rules in `[digest]` (6–30 residues, 0 missed cleavages: Schwanhäusser et al. 2011 as restated in Nature Sci Data 2018;5:180128), on the full gene product, and writes per entry the theoretical peptide count N_i, `peptides_per_kDa` = N_i / (MW_full,i / 1000), and the set of peptides shared with other entries. Two outputs feed Layer B:

- **Density.** `outputs/digest/density_ranked.tsv` ranks entries by `peptides_per_kDa` with a robust z-score ((x − median) / 1.4826·MAD). No threshold is applied by code (D43). Murgia 2021 does not state whether cleavage before proline was allowed in the search; both counts are computed and their per-entry difference written (D49). Run 2026-09-12T19:22:30Z: 303 entries; peptides per kDa median 0.530, MAD 0.061, range 0–0.942 (the zero is a 35-residue entry with no in-window peptide); proline-rule relative difference median 1.1 %, ≥ 10 % for 31 entries, maximum 23.8 %.
- **Families.** `outputs/digest/families.tsv` lists connected components of entries sharing at least one in-window peptide. Within a family, a search engine's assignment of shared-peptide intensity can move signal between members (Cox et al. 2014, on protein-group quantification). For each family the Murgia 2021 fiber-type ratio is compared with Moreno-Justicia 2025's, which quantified "proteotypic peptides" only (Moreno-Justicia 2025, Methods, "MS data processing"); disagreement is recorded as a per-family Layer B uncertainty term (D44). Run 2026-09-12T19:22:30Z: 18,039 distinct theoretical peptides, 530 shared across 141 pairs, 47 of them single-peptide links; 34 connected components, the largest being the ten MYH entries plus two single-peptide attachments. Pairwise counts and shared fractions are in `outputs/digest/shared_pairs.tsv`; the family comparison at the mass-fraction run is made per pair, not per component.

### Fiber type

Fiber types are Murgia 2021's own: ≥ 80 % of one MYH by unique-peptide intensity (Murgia 2021, Methods). The per-type standard is defined on those pure fibers; hybrid fibers enter only through the aggregation step's fiber-mix step (D41). Moreno-Justicia et al. 2025 report that in 1038 vastus lateralis fibers "pure type 2X fibers could only be detected at the RNA but not at the protein level" and that MYH1-positive fibers did not separate from other fast fibers on whole-proteome clustering (Moreno-Justicia 2025, Results, "Type 2X is not a distinct fiber type"). The IIx column of this standard therefore rests on Murgia 2021's 13 fibers; whether the published standard reports I / IIa / IIx or slow / fast is decided at aggregation (D46).

### Cross-checks and their questions

Each cross-check is admitted only with a named question (D42):

| Source | Question it answers | Output |
|--------|---------------------|--------|
| Moreno-Justicia et al. 2025 | Are Murgia 2021's slow/fast ratios (n = 61, 4 donors) reproduced at n = 1038, 5 donors, different instrument and search? | `outputs/mass_fractions/ratio_check_moreno-justicia_2025.tsv` |
| Deshmukh et al. 2021 | Same question, pooled fibers, third workflow | `outputs/mass_fractions/ratio_check_deshmukh_2021.tsv` |
| Carroll, Carrithers & Trappe 2004 | Does iBAQ×MW reproduce the MHC:actin ratio measured by quantitative SDS-PAGE against standards in human single fibers, types I and IIa (Carroll 2004, abstract: 252 fibers for MHC, 160 for actin, five donors)? | `outputs/mass_fractions/classical_check_carroll_2004.tsv` |

Carroll 2004's values are transcribed into `config/carroll_classical_fractionation.ini` with page and table per row (D50); they feed the check only and never the weights (D39). Carroll et al. 2004 report MHC and actin per total fiber protein, so the MHC:actin ratio compares directly and the absolute fractions await the Tier 1 : total ratio (aggregation step). No usable human weighed-mass source exists for titin or nebulin (Trappe et al. 2002 was examined: uncalibrated densitometry, arbitrary units); the giant-protein question rests on the disclosed density rank alone (D39, D43).

### Layer B uncertainty terms (inputs to aggregation)

1. Between-dataset ratio disagreement, per protein (Murgia 2021 vs Deshmukh 2021 vs Moreno-Justicia 2025).
2. iBAQ×MW vs weighed values (Carroll 2004), for myosin heavy chain and actin in types I and IIa.
3. Shared-peptide families, per family (Murgia 2021 vs Moreno-Justicia 2025).
4. iBAQ×MW vs TPA — per-protein factor `peptides_per_kDa`; applied only if raw intensities become available.
5. Between-fiber spread within a fiber type — Dataset 1's `Standard deviation` column per fiber type (19 / 29 / 13 fibers), carried per row (D45). Caveats: SD paired with a median; the 2X value rests on 13 fibers.

Each term is written per entry as `w_low`, `w_high` with the source of the interval named, in `config/mass_fractions/<fiber_type>.ini` (generated).

### Carried checks (from the protein-set and composition stages)

Computed in `pipeline/mass_fractions.py --checks` after weights exist: (a) weighted isoform bound = per-entry bound × w_i, largest per amino acid; (b) weighted processing bound, same; (c) weighted PTM mass bound, same, and Σ w_i over glycosylated entries; (d) Σ w_i over entries reached only via the cardiac-myofibril branch (`data/tier1_pool.tsv`, column `subtree_terms_returning_entry`); (e) the 26 dual-tier entries: with within-tier normalisation each carries a separate w in each tier — no rule needed, recorded; (f) Tier 2 redefinition from measurement: Dataset 1 rows not in Tier 1, ranked by v·MW share — decision D⟨n⟩ once the list exists.
