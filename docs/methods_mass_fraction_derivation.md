# Mass fraction derivation (Layer B)

_Section of `docs/methods.md`. Drafted Step 3a, 2026-09-11, revised 2026-09-12 after Dataset 1 was read, under `PROVENANCE.md`. Every sentence below either carries a citation (source id from `config/mass_fraction_decisions.ini`) or names a computation in this repository. Angle-bracket placeholders ⟨…⟩ are filled from run records in Step 3b._

## Inputs

Layer B reads three things: the protein set and full-product molecular weights from Steps 1–2 (`config/accessions.ini`; `outputs/composition/amino_acid_composition_per_protein.tsv`, column `mw`); the literature files fetched and hashed by `pipeline/fetch_literature.py` (`data/literature/manifest.ini`); and the hand-written `config/mass_fraction_decisions.ini`, which holds source citations, role decisions, column-to-fiber-type mappings, and digest rules — and, by test, no protein and no abundance number.

## Primary dataset

The primary dataset is A1 (Murgia et al. 2021; role by D32). Its quantification is described in its Methods: "the summed intensity of the peptides of each protein was divided by the number of theoretically observable peptides [8], and the resulting values were normalized to the expression of α-skeletal actin (ACTA1)" (A1, Methods, "Single-fiber proteomics"). The reference [8] in that sentence is B5 (Cox et al. 2014). Dataset 1 reports, per protein and fiber type, "the median expression values and the percent of the maximal value" and "the number of fibers in which each protein was identified (valid values)"; "proteins are indicated by gene name" (A1, Results, first paragraph). Fibers were "relatively pure myofibers containing at least 80% of either MYH7 or MYH2 or MYH1. A total of 61 fibers (19 type 1, 29 type 2A, and 13 type 2X fibers)" from "four younger (aged 22–27 years)" donors; "MYH expression was quantified by the intensities of peptides unique for each isoform" (A1, Methods). Raw data are PXD006182 (A1, Availability of data), shared with A2 (Murgia et al. 2017), which is therefore not an independent source (D33).

## Quantity used

The quantity intensity ÷ theoretical peptide count is iBAQ (B4, Schwanhäusser et al. 2011). Mass fraction of protein *i* within tier *T* is computed by `pipeline/mass_fractions.py` as

    w_i = v_i · MW_i / Σ_{j ∈ T} ( v_j · MW_j )

where v_i is the median value read from the column mapped in `[columns.murgia_2021.file.1]` and MW_i is the Step 2 molecular weight of the mature chain (`segment_set = master`). The ACTA1 normalisation in v_i is a per-fiber constant and cancels in the ratio. Sums are taken within tier: Σ_{Tier 1} w = 1 and Σ_{Tier 2} w = 1 (D31). A pool entry with no row in Dataset 1 receives w = 0 and is listed (`outputs/mass_fractions/unmatched_pool_entries.tsv`); a Dataset 1 row with no pool entry is listed with its v·MW share of the whole file (`outputs/mass_fractions/dataset_rows_outside_pool.tsv`) — this is the pool-completeness check carried from Step 1.

Why not LFQ: MaxLFQ is defined for comparing a protein across samples, and A1 used it for its statistics (A1, Methods, citing B5); A3 (Momenzadeh et al. 2023) reports MYH isoform fractions computed both by iBAQ (its Table S2) and by MaxLFQ (its Table S3) from the same fibers, and the two disagree by ⟨factor, from the downloaded tables⟩ for MYH4. Datasets that publish only LFQ-type values (A5, DIA-NN protein-group LFQ; A4, MaxLFQ per its Methods) are therefore used for within-protein ratios between fiber types, never for cross-protein weights (D34, D36).

Why not TPA: the total-protein approach (B1, Wiśniewski 2017) reads mass fraction as raw intensity ÷ summed intensity. A1 publishes iBAQ-type values, not raw intensities, so TPA is not computable from the published file. iBAQ × MW and TPA differ per protein exactly by the factor N_i / MW_i (theoretical peptides per unit mass); that factor is computed from sequence alone by `pipeline/digest.py` and reported as `peptides_per_kDa` (`outputs/digest/theoretical_peptides.tsv`).

## Identity

Dataset 1 identifies proteins by gene name only (`Gene names`, sheet `Sheet1`, header row 3; no accession column). A row is joined to a pool entry by exact gene symbol against the `gene` column of `config/accessions.ini` (rule B1, `pipeline/mass_fractions.py`). Where one gene name occurs on several rows — separate MaxQuant protein groups of one gene product; 85 such names — the rows are summed and the share contributed by the minor rows is written per gene (D47). Where one cell lists several genes separated by `;` — 139 rows, 3.0 % of summed medians across the file — the value is assigned to the single pool member if only one is in the pool, or divided among pool members in proportion to their own single-name rows with the whole row carried as each member's upper bound (D48); the signal share arriving this way is printed per tier and fiber type. A pool gene with no row receives w = 0 and is listed; a Dataset 1 gene with no pool entry is listed with its v·MW share (the pool-completeness check). Dataset 1 was searched against "the human UniProt FASTA reference proteomes version of 2016" (A1, Methods); gene symbols renamed since then surface in the unmatched lists. ⟨Counts from the run.⟩

## In-silico digest (Layer A computation used by Layer B)

`pipeline/digest.py` digests every canonical sequence under the rules in `[digest]` (6–30 residues, 0 missed cleavages: B4 as restated in Nature Sci Data 2018;5:180128), on the full gene product, and writes per entry the theoretical peptide count N_i, `peptides_per_kDa` = N_i / (MW_full,i / 1000), and the set of peptides shared with other entries. Two outputs feed Layer B:

- **Density.** `outputs/digest/density_ranked.tsv` ranks entries by `peptides_per_kDa` with a robust z-score ((x − median) / 1.4826·MAD). No threshold is applied by code (D43). A1 does not state whether cleavage before proline was allowed in the search; both counts are computed and their per-entry difference written (D49). Run 2026-09-12T19:22:30Z: 303 entries; peptides per kDa median 0.530, MAD 0.061, range 0–0.942 (the zero is a 35-residue entry with no in-window peptide); proline-rule relative difference median 1.1 %, ≥ 10 % for 31 entries, maximum 23.8 %.
- **Families.** `outputs/digest/families.tsv` lists connected components of entries sharing at least one in-window peptide. Within a family, a search engine's assignment of shared-peptide intensity can move signal between members (B5, Cox et al. 2014, on protein-group quantification). For each family the A1 fiber-type ratio is compared with A5's, which quantified "proteotypic peptides" only (A5, Methods, "MS data processing"); disagreement is recorded as a per-family Layer B uncertainty term (D44). Run 2026-09-12T19:22:30Z: 18,039 distinct theoretical peptides, 530 shared across 141 pairs, 47 of them single-peptide links; 34 connected components, the largest being the ten MYH entries plus two single-peptide attachments. Pairwise counts and shared fractions are in `outputs/digest/shared_pairs.tsv`; the family comparison in Step 3b is made per pair, not per component.

## Fiber type

Fiber types are A1's own: ≥ 80 % of one MYH by unique-peptide intensity (A1, Methods). The per-type standard is defined on those pure fibers; hybrid fibers enter only through Step 4's fiber-mix step (D41). A5 reports that in 1038 vastus lateralis fibers "pure type 2X fibers could only be detected at the RNA but not at the protein level" and that MYH1-positive fibers did not separate from other fast fibers on whole-proteome clustering (A5, Results, "Type 2X is not a distinct fiber type"). The IIx column of this standard therefore rests on A1's 13 fibers; whether the published standard reports I / IIa / IIx or slow / fast is decided in Step 4 (D46).

## Cross-checks and their questions

Each cross-check is admitted only with a named question (D42):

| id | Source | Question it answers | Output |
|----|--------|---------------------|--------|
| A5 | Moreno-Justicia et al. 2025 | Are A1's slow/fast ratios (n = 61, 4 donors) reproduced at n = 1038, 5 donors, different instrument and search? | `outputs/mass_fractions/ratio_check_A5.tsv` |
| A4 | Deshmukh et al. 2021 | Same question, pooled fibers, third workflow | `outputs/mass_fractions/ratio_check_A4.tsv` |
| H1 | Carroll, Carrithers & Trappe 2004 | Does iBAQ×MW reproduce the MHC:actin ratio measured by quantitative SDS-PAGE against standards in human single fibers, types I and IIa (H1, abstract: 252 fibers for MHC, 160 for actin, five donors)? | `outputs/mass_fractions/classical_check_H1.tsv` |

H1's values are transcribed into `config/classical_fractionation.ini` with page and table per row; they feed the check only and never the weights (D39). H1 reports MHC and actin per total fiber protein, so the MHC:actin ratio compares directly and the absolute fractions await the Tier 1 : total ratio (Step 4). No usable human weighed-mass source exists for titin or nebulin (Trappe et al. 2002 was examined: uncalibrated densitometry, arbitrary units); the giant-protein question rests on the disclosed density rank alone (D39, D43).

## Layer B uncertainty terms (inputs to Step 4)

1. Between-dataset ratio disagreement, per protein (A1 vs A4 vs A5).
2. iBAQ×MW vs weighed values (H1), for MHC and actin in types I and IIa.
3. Shared-peptide families, per family (A1 vs A5).
4. iBAQ×MW vs TPA — per-protein factor `peptides_per_kDa`; applied only if raw intensities become available.
5. Between-fiber spread within a fiber type — Dataset 1's `Standard deviation` column per fiber type (19 / 29 / 13 fibers), carried per row (D45). Caveats: SD paired with a median; the 2X value rests on 13 fibers.

Each term is written per entry as `w_low`, `w_high` with the source of the interval named, in `config/mass_fractions/<fiber_type>.ini` (generated).

## Carried checks (from Steps 1–2)

Computed in `pipeline/mass_fractions.py --checks` after weights exist: (a) weighted isoform bound = per-entry bound × w_i, largest per amino acid; (b) weighted processing bound, same; (c) weighted PTM mass bound, same, and Σ w_i over glycosylated entries; (d) Σ w_i over entries reached only via the cardiac-myofibril branch (`data/tier1_pool.tsv`, column `subtree_terms_returning_entry`); (e) the 26 dual-tier entries: with within-tier normalisation each carries a separate w in each tier — no rule needed, recorded; (f) Tier 2 redefinition from measurement: Dataset 1 rows not in Tier 1, ranked by v·MW share — decision D⟨n⟩ in 3b once the list exists.
