# The calculation, by hand

This is the complete arithmetic behind `outputs/standard/_calculated_amino_acid_standard.tsv`,
written so that anyone with the input files and a calculator can reproduce any number in it.
Every quantity below is either read from a named file or computed from the line before.
The code does nothing else.

Two layers, then a mixture:

- **Layer A** — what each protein is made of. From its sequence alone.
- **Layer B** — how much of each protein the tissue holds. From the measurement.
- **The standard** — the amino acids of the mixture, added up.

Plain-language names used throughout: a **fiber** is one muscle cell; it comes in three
**fiber types** (I, IIa, IIx). Inside every fiber are two kinds of protein: the
**contractile** proteins (the machinery that pulls) and the **builders** (everything else the
cell makes and keeps — pumps, enzymes, the machinery that builds and repairs). A fiber is about
72 % contractile and 28 % builders by protein mass, in every fiber type (as measured, below).

## Layer A — one protein

**Step 1. Count the letters.** Take the canonical sequence of the entry from UniProt
(`data/uniprot_sequences.ini`), keep the mature chain (`config/segments.ini`, rule R2), and
count how many times each of the twenty amino acid letters appears. Call the counts
c_A, c_C, …, c_Y.

**Step 2. Look up the masses.** `data/pubchem/amino_acid_masses.ini` gives, for each amino
acid, its **free mass** (the whole molecule) and its **residue mass** (free mass minus one
water, 18.015 g/mol, because a residue inside a chain has given up a water to form its
peptide bond). Both were fetched from PubChem by code.

**Step 3. The protein's molecular weight.**

    MW = Σ (count × residue mass) + 18.015

The one water at the end is the chain's two free ends.

That is all Layer A is. The per-protein rows are in
`outputs/composition/amino_acid_composition_per_protein.tsv`.

## Layer B — how much of each protein

**Step 4. The measured value.** Murgia 2021, Dataset 1, gives for each protein a median
intensity per fiber type, *v* (an iBAQ-type value: intensity divided by the number of
theoretical peptides, which makes it proportional to the number of molecules). Rows sharing a
gene are summed (D47); multi-gene rows are split (D48). These are the `v_I`, `v_IIa`, `v_IIx`
columns of `config/mass_fractions_per_entry.tsv`.

**Step 5. From molecules to mass.** Molecules × molecular weight = mass. The **weight** of
protein *i* in a set of proteins is its share of the set's mass:

    w_i = v_i × MW_i ÷ Σ_j (v_j × MW_j)

summed over the set (rule B2, D31). Three sets are used: the contractile proteins alone, the
builders alone, and the **total** (both together under one denominator, D58). The weights of a
set sum to one. The total set's weights are the `w_<type>_combined` columns of
`config/mass_fractions_per_entry.tsv`; the per-tier weights are `w_<type>_tier1` and
`w_<type>_tier2`.

## The standard — adding up the mixture

**Step 6. Moles of each amino acid per gram of protein.** One gram of the set contains w_i
grams of protein *i*, which is w_i ÷ MW_i moles of it, each molecule carrying c_{i,a}
residues of amino acid *a*. So

    n_a = Σ_i (w_i ÷ MW_i) × c_{i,a}        [mol of amino acid a per gram of protein]

This is rule A1 (D64): the mixture is added up by molecules, not by percentages, because a
gram of a giant protein holds fewer molecules than a gram of a small one.

**Step 7. Grams per 100 g protein, and percent.**

    g_a per 100 g protein  = n_a × m_a × 100
    percent of amino acids = n_a × m_a ÷ Σ_b (n_b × m_b) × 100

with *m_a* the **free mass** (the standard; what laboratories and food tables report; the
basis Match Rate uses) or the **residue mass** (the companion file). In the free convention
the grams sum to about 116 per 100 g protein — the water taken up on hydrolysis — and the
percent column always sums to 100.

**Step 8. The final column.** The nine columns are one per set and fiber type. The final
column mixes the three totals by the share of a muscle's protein in each fiber type:

    standard = share_I × total_I + share_IIa × total_IIa + share_IIx × total_IIx

The shares are not in the dataset (its fibers were selected pure, not sampled from a muscle);
they come, with their citation, from `config/fiber_type_mix.ini`. Until that file is filled the
column is blank.

## The bound pool — one more thing the muscle holds

The steps above count the amino acids held as protein. Muscle also holds one amino acid in a
stable non-protein form: histidine, as the dipeptide carnosine. That histidine had to be eaten
to be there, so the standard with bound pools adds it (D73–D77). The protein-only table above is
not changed; the adjusted table sits beside it as
`_calculated_amino_acid_standard_with_non_protein_metabolite_pools.tsv`.

**Step 9. Put the pool and the protein on the same kilogram of muscle.** The pool is cited per
kilogram of muscle (`config/non_protein_metabolite_pools.ini`: c mmol/kg, on the source's basis — dry
or wet muscle); the protein content P (g/kg) is cited on its own basis. If the bases differ, the
cited water content W (g per kg wet muscle) moves one to the other:

    P_dry = P_wet ÷ (1 − W ÷ 1000)

**Step 10. Grams of each amino acid per kilogram, from protein.** Step 7 gave g_a per 100 g
protein (free convention). Then

    F_a = P × g_a ÷ 100          [g of free amino acid a from protein, per kg muscle]

**Step 11. Grams of the pool's amino acid per kilogram, from the pool.** With n the moles of
the amino acid released per mole of pool (one for carnosine) and m its free mass (PubChem):

    C = c ÷ 1000 × n × m         [g of the amino acid held as the pool, per kg muscle]

**Step 12. Add, then renormalise.** The pool's amino acid becomes F + C; every other amino acid
stays F; the percent column is each over the new sum:

    percent_a = F_a ÷ (Σ_b F_b + C) × 100   (for the pool's amino acid, (F_a + C) in the numerator)

This is done per fiber type with that fiber type's pool value, for the total column; the
contractile and builders columns are protein-only by definition and are copied; the final column
mixes the three adjusted totals by the same fiber-type shares as step 8.

**A worked line.** Suppose a fiber type's protein-only profile gives 2.50 g of free histidine
per 100 g protein and 116.0 g of free amino acids in all; the muscle holds P = 800 g protein per
kg dry muscle and c = 20 mmol carnosine per kg dry muscle; histidine's free mass is 155.16.
Then F_His = 800 × 2.50 ÷ 100 = 20.00 g/kg; Σ F = 800 × 116.0 ÷ 100 = 928.0 g/kg;
C = 20 ÷ 1000 × 1 × 155.16 = 3.10 g/kg; histidine goes from 20.00 ÷ 928.0 = 2.155 % to
(20.00 + 3.10) ÷ (928.0 + 3.10) = 2.481 %, and every other amino acid is multiplied by
928.0 ÷ 931.1. The ratio C ÷ F_His = 0.155 is written in `non_protein_metabolite_pool_amounts_per_kg_muscle.tsv`.
(These are illustrative numbers; the real ones are in the config and the output headers.)

**The frame.** This is an inventory: what the muscle holds. The resupply frame — how much of
each pool the body has to replace per day — would multiply C by k_pool ÷ k_protein, the two
replacement rates; the cited ranges for those rates are wide, so that frame is reported as a
sensitivity (`sensitivity_non_protein_metabolite_pool_turnover_frame.tsv`), not as the standard (D73).

## A worked example

Two invented proteins, real masses (from `data/pubchem/amino_acid_masses.ini`: glycine
75.07 / 57.055, alanine 89.09 / 71.075, valine 117.15 / 99.135, lysine 146.19 / 128.175 g/mol,
free / residue; water 18.015):

| | Protein P: sequence G-A-V | Protein Q: sequence A-A-K-K |
|---|---|---|
| Step 1, counts | G 1, A 1, V 1 | A 2, K 2 |
| Step 3, MW | 57.055 + 71.075 + 99.135 + 18.015 = **245.280** | 2 × 71.075 + 2 × 128.175 + 18.015 = **416.515** |
| Step 4, measured *v* (invented) | 200 | 50 |
| v × MW | 49,056.0 | 20,825.75 |
| Step 5, weight w | 49,056.0 ÷ 69,881.75 = **0.70199** | 20,825.75 ÷ 69,881.75 = **0.29801** |
| w ÷ MW (mol of protein per g) | 0.0028620 | 0.00071550 |

Step 6, moles of each amino acid per gram of the mixture:

| amino acid | n_a | from |
|---|---|---|
| G | 0.0028620 × 1 = 0.0028620 | P |
| A | 0.0028620 × 1 + 0.00071550 × 2 = 0.0042930 | P and Q |
| V | 0.0028620 × 1 = 0.0028620 | P |
| K | 0.00071550 × 2 = 0.0014310 | Q |

Step 7, free convention:

| amino acid | n_a × free mass × 100 (g per 100 g protein) | percent |
|---|---|---|
| G | 0.0028620 × 75.07 × 100 = 21.49 | 18.82 |
| A | 0.0042930 × 89.09 × 100 = 38.25 | 33.50 |
| V | 0.0028620 × 117.15 × 100 = 33.53 | 29.36 |
| K | 0.0014310 × 146.19 × 100 = 20.92 | 18.32 |
| sum | 114.18 | 100.00 |

Residue convention, same n_a: G 16.33, A 30.51, V 28.37, K 18.34 g per 100 g protein, sum
93.56 (100 less the two chains' terminal water, weighted); percent 17.45 / 32.61 / 30.33 / 19.61.

Note what the percentages do *not* depend on: because w_i ÷ MW_i = v_i ÷ Σ(v × MW), the
molecular weights cancel between steps 5 and 6, and n_a is proportional to Σ_i v_i × c_{i,a}
— the measured molecule counts times the letter counts. The molecular weight matters for the
mass shares (which protein is "big"), not for the amino acid profile of the mixture.

## Doing it for a real entry

1. `config/mass_fractions_per_entry.tsv` — the entry's row: `v_I`, `mw_master`, `w_I_combined`.
2. `outputs/composition/amino_acid_composition_per_protein.tsv`, the entry's `master` row —
   its twenty `count_*` values.
3. `data/pubchem/amino_acid_masses.ini` — the twenty free and residue masses.
4. Steps 6–7 above over every row of the weights table give `_calculated_amino_acid_standard.tsv`;
   the `sum_of_weights` and `max_abs_a1_minus_a1c` columns of `amino_acid_profiles.tsv` are
   the checks the code runs on itself.

Every file named here carries the SHA-256 of what it was computed from in its header, so the
chain can be followed back to the UniProt release and the dataset file.
