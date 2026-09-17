# The calculated standard beside an independent measurement

**Gorissen et al. 2018 · human skeletal muscle protein · what agrees, what does not, and what the two methods cannot be asked to agree about**

Every number in this document comes from `outputs/comparison/calculated_vs_gorissen_2018.tsv`,
written by `pipeline/comparison.py` (`python run.py comparison`). Nothing here is computed by
hand. The source is `data/literature/gorissen_2018/726_2018_Article_2640.pdf`, sha256
`6f75965f…215f`, transcribed in `config/gorissen_2018_comparison.ini` with a page-level location
on every value (D79, D84).

---

## 1. Why this paper

This repository derives the amino acid composition of human skeletal muscle protein from protein
sequences and measured protein abundances. That is a computation, and a computation of something
that can also be measured directly. Gorissen et al. 2018 measured it: human *m. vastus lateralis*
from ten volunteers, freeze-dried, acid-hydrolysed, quantified by UPLC–MS/MS, published as one
column of Table 1 (p. 1690). They included human muscle deliberately, as the reference protein
with an "ideal" amino acid composition for muscle protein synthesis (p. 1686).

It is also the column the original Match Rate spreadsheet was calibrated against. So this
comparison answers two questions at once: does the sequence-derived standard land near an
independent measurement, and where will a Match Rate computed against the new reference differ
from one computed against the old.

The paper is a good comparator for three reasons beyond convenience. The muscle sample went
through the same instrument and the same procedure as the fifteen food proteins in the same
table, so its relationship to whey, egg and soy is internally consistent. The method is stated in
enough detail to know what it can and cannot see. And it reports its own totals, which means the
transcription can be checked rather than trusted — see §6.

---

## 2. What the two methods can be asked to agree about

Five of the twenty amino acids are absent from the measured side, and one pair cannot be
separated. Three of the exclusions are the paper's own statements about its method. Two are this
repository's reading of two printed zeros, argued in §3 and decided in D85.

| Amino acid | Why it is not a comparable number | Basis |
|---|---|---|
| Tryptophan (W) | Decomposed during acid hydrolysis, which precludes detection. The paper's own essential-amino-acid sums exclude it for this reason. | p. 1687, right column, last paragraph |
| Aspartic acid (D) | Listed as not measured. | p. 1690, Table 1 footnote |
| Asparagine (N) | Converted into aspartic acid during hydrolysis — and aspartic acid was not measured, so neither survives. | p. 1687, right column, last paragraph |
| **Proline (P)** | Printed as 0.0. Not a measurement: proline is not destroyed by acid hydrolysis, is never named in the paper among the amino acids its method struggles with, and the fifteen other sources in the same table carry 1.8–8.8 g/100 g. | §3, D85 |
| **Cysteine (C)** | Printed as 0.0. Not a measurement: the paper names cysteine as at risk and says hydrolysis was stopped at 12 h to limit its loss. The mass balance does not admit a true zero either. | §3, D85 |
| Glutamine (Q) | Converted into glutamic acid during hydrolysis. The glutamic acid row therefore holds both, so E and Q are summed on both sides. | p. 1687, right column, last paragraph |

So the comparison covers **fourteen rows across fifteen amino acids**: thirteen compared
individually, plus one combined `E+Q` row. W, D, N, P and C are excluded from *both* sides — not
just from the measured one. That symmetry is the whole point: an amino acid present in one
denominator and not the other biases every row in the table, which is exactly the error the first
version of this document made (§3).

**Both sides are then renormalised over that set.** The paper's own headline figures — leucine
7.6 % of total protein, essential amino acids 38 % — use a denominator that includes what it did
not measure, so they cannot be compared with this repository's directly. The comparison set is
**83.69 %** of the calculated standard's twenty.

## 3. The two zeros, and why they are not measurements

Table 1 prints 0.0 for proline and 0.0 for cysteine in the human muscle column. The first version
of this comparison carried both as measured zeros. That was wrong, and it was wrong in a way that
touched every other row.

**A protein with no proline does not exist.** Proline is among the most stable amino acids under
6 M HCl — the amino acids this method loses are tryptophan (destroyed), cysteine and methionine
(oxidised), serine and threonine (partly), and asparagine and glutamine (deamidated). Proline is
on none of those lists, the paper never names it, and the fifteen other protein sources in the
same table — same instrument, same procedure, same run — report proline from 1.8 to 8.8 g/100 g.
The method measures proline perfectly well. Only this one cell is zero.

**The mass balance will not admit either zero.** Acid hydrolysis breaks every peptide bond and
adds a water to each, so a gram of protein yields *more* than a gram of free amino acids — for
this composition, a factor of 1.160, computed in `gorissen_2018_mass_balance.tsv` from the
PubChem masses and this repository's own standard, with nothing taken from the paper.

| Protein content used | Expected free amino acids | What the method could measure | Reported | Recovery |
|---|---:|---:|---:|---:|
| Gorissen's 84 % (nitrogen × 6.25) | 97.5 g | 81.6 g | 60.8 g | **75 %** |
| Mingrone 2001's 65.4 % (measured) | 75.9 g | 63.5 g | 60.8 g | **96 %** |

Their sixteen printed values total 60.8 g per 100 g of tissue against a stated 84 g of protein.
That is 36.7 g short of the expected yield — and treating proline, cysteine, aspartate,
asparagine and tryptophan as entirely absent accounts for only 15.9 g of it. A quarter of the
expected mass is simply not in the table.

The second row is this repository's own D77 source: Mingrone et al. 2001 measured muscle protein
directly, 155 g/kg wet against 763 g/kg water, which is 65.4 % of dry mass. Gorissen's 84 % is
nitrogen × 6.25 on whole freeze-dried muscle, and that factor counts non-protein nitrogen —
creatine, carnosine, free amino acids, nucleotides, urea. Use the directly measured protein
content and treat the two zeros as unmeasured, and the balance closes to **96 %**.

Two independent corrections landing together is not proof, and the stage draws no conclusion from
it (X3). But it is a second line of evidence for what Step 4b found from the other direction: the
muscle column measures *tissue*, not protein.

## 3b. The comparison

Percent of the fourteen-row comparison set. The calculated standard is shown in both versions:
protein only, and with non-protein metabolite pools (carnosine, D75/D82). Difference is
calculated minus measured, in percentage points.

| Row | Gorissen | Calculated (protein only) | Calculated (with pools) | Difference | Ratio |
|---|---:|---:|---:|---:|---:|
| E+Q | 21.55 | 19.08 | 19.01 | −2.54 | 0.88× |
| K | 10.86 | 11.33 | 11.28 | +0.43 | 1.04× |
| L | 10.36 | 10.08 | 10.04 | −0.32 | 0.97× |
| R | 7.24 | 8.27 | 8.23 | +1.00 | 1.14× |
| V | 7.07 | 7.11 | 7.08 | +0.01 | **1.00×** |
| A | 6.74 | 6.49 | 6.46 | −0.28 | 0.96× |
| F | 6.25 | 4.83 | 4.81 | −1.44 | 0.77× |
| I | 5.59 | 6.52 | 6.49 | +0.90 | 1.16× |
| G | 5.10 | 4.09 | 4.07 | −1.03 | 0.80× |
| T | 4.77 | 6.08 | 6.05 | +1.28 | 1.27× |
| H | 4.61 | 2.66 | 3.07 | −1.54 | 0.67× |
| S | 3.78 | 5.69 | 5.67 | +1.88 | 1.50× |
| Y | 3.29 | 4.58 | 4.56 | +1.27 | 1.39× |
| M | 2.80 | 3.19 | 3.17 | +0.38 | 1.14× |

**Neither column is treated as the truth.** Nothing is fitted, scaled or adjusted toward the other
side, and no error term is computed (X3). Where the two disagree, the disagreement is the finding.

### Where they agree

Valine lands at 1.00×, leucine at 0.97×, alanine at 0.96×, lysine at 1.04×. Four of the five
largest amino acids by mass agree within 4 %. For a sequence-derived computation and a wet-lab
hydrolysis of different tissue from different people, that is closer than it had to be.

### Where they do not

**Glutamate + glutamine, 0.88×.** The largest remaining difference, and the one row where the
measured side is a sum of two amino acids the method could not separate — so any systematic
recovery difference between the two lands in one place.

**Phenylalanine 0.77× and histidine 0.67×**, both rows where the calculated standard asks for less
than the measurement does. **Serine 1.50×, tyrosine 1.39×, threonine 1.27×, isoleucine 1.16×**,
where it asks for more. Serine and threonine are both partly lost to acid hydrolysis, which points
in the right direction for those two but is not quantified here. Phenylalanine, tyrosine and
isoleucine have no such explanation in either method. They are the honest open items.

## 4. Histidine, and what carnosine does to it

Histidine deserves its own section because it is the amino acid Step 4b was about.

| | Percent of the comparison set |
|---|---:|
| Gorissen, measured | 4.61 |
| Calculated, protein only | 2.66 |
| Calculated, with non-protein metabolite pools | 3.07 |

Carnosine closes **21 %** of the protein-only gap — 0.41 of 1.95 percentage points. The remaining
1.54 points is the second largest unexplained difference in the table, after E+Q.

One observation from the paper's own table supports the direction of the adjustment without
settling the size of it. **Human muscle has the highest histidine of all sixteen protein sources
Gorissen measured**: 2.8 g/100 g raw material, against 2.2 for calcium caseinate, 1.9 for milk,
and 0.7–1.7 for everything else. A protein whose sequences contain no more histidine than casein
does, measuring higher than casein, is what you would expect if the sample carried histidine in a
form other than protein. The sample was freeze-dried whole muscle tissue, hydrolysed — and
hydrolysis releases carnosine's histidine exactly as it releases protein histidine. The measured
column is therefore not a measurement of muscle *protein* histidine alone; it is a measurement of
the histidine in muscle tissue, which is the quantity the adjusted standard is trying to
represent.

That is an argument for why the adjusted column is the right one to compare, not evidence that
its value is correct. The gap that remains after carnosine is not explained here, and the
candidates — the blood category, other histidine-bearing compounds, the protein set's coverage,
and the measurement itself — are for later steps. No verdict (X3).

---

## 5. What a Match Rate user should take from this

1. **Use the adjusted standard.** The measured comparator counts tissue histidine, not protein
   histidine, so the version of the standard that includes the non-protein metabolite pools is
   the one that answers the same question the measurement does.
2. **Scores against the new reference will differ most where the reference differs most.** Of the
   nine indispensable amino acids, phenylalanine (0.77×) and histidine (0.67×) are the rows where
   the calculated standard asks for *less* than Gorissen's column does, and threonine (1.27×),
   isoleucine (1.16×) and methionine (1.14×) are where it asks for more. A food scored against the
   new reference is less likely to be limited by phenylalanine or histidine than under the old
   one, and more likely to be limited by threonine or isoleucine.
3. **Proline and cysteine are not usable from this source at all.** Neither is indispensable, so
   v1.0 is unaffected; but any future all-twenty output must not take these two cells as
   measurements, and must not silently renormalise around them either.
4. **The paper's headline percentages are not this table's percentages.** Leucine "7.6 % of total
   protein" in the paper and 10.04 % here are the same measurement over different denominators.
   Quoting one beside the other would be an error.

---

## 6. How the transcription was checked

Table 1 prints its own column sums: 31.8 for the essential amino acids and 29.0 for the
non-essential. The stage recomputes both from the transcribed values before it computes anything
else, and stops if either misses by more than 0.05 (rule X4). Both reproduce exactly:

```
essential sum 31.8 = printed 31.8; non-essential sum 29.0 = printed 29.0
```

This is a check on the transcription, not on the paper. It catches a mistyped digit, a
transposed row, or a value read from the wrong column — the failure modes of a hand
transcription — before any comparison is drawn from it.

A second guard was added after the first version of this document shipped with the two zeros
compared. **An undeclared zero now stops the stage** (X5). A measured value of 0.0 must be either
declared in `[not_measured]` with the reason it is not a measurement, or declared in
`[zero_is_a_measurement]` with the evidence that it is one. Neither this repository nor a later
reader can compare a zero without having answered that question. And the mass balance (X6) is
computed on every run rather than being an argument made once in prose.

Three further guards run before the table is written: the cited PDF must be on disk and hash to
its manifest entry or nothing is computed (D78); every one of the twenty amino acids must be
either compared or listed as not measured with a stated reason, so none can go missing unnoticed;
and each row label in the config is resolved to a one-letter symbol through the IUPAC–IUBMB
trivial names, so no amino acid is named in code.

---

## 7. What this comparison is not

It is not a validation. Two measurements of the same quantity that agree tell you they agree;
they do not tell you either is right, and this one is a single published column from ten
volunteers against a computation from a different set of people's proteomes.

It is not a calibration. Nothing in this repository is adjusted toward Gorissen's values, and the
stage is built so that it cannot be: it has no path that writes back into the standard.

And it is not complete. Five of the twenty cannot be compared, a quarter of the expected amino
acid mass is missing from the source table under its own protein figure, and six rows disagree by
more than 10 % for reasons neither method explains. Those are on the record here rather than
resolved.

---

## Files

| Path | What |
|---|---|
| `config/gorissen_2018_comparison.ini` | the transcription, with a page-level location on every value |
| `pipeline/comparison.py` | the stage; rules X1–X5 in `docs/conventions.md` |
| `outputs/comparison/gorissen_2018_human_muscle.tsv` | Table 1's human muscle column as letters |
| `outputs/comparison/calculated_vs_gorissen_2018.tsv` | the comparison behind §3b |
| `outputs/comparison/gorissen_2018_mass_balance.tsv` | the reconciliation behind §3 |
| `outputs/comparison/comparison_summary.ini` | the set, the exclusions, the transcription check, the largest differences |
| `tests/test_comparison.py` | twelve tests, one per rule plus the guards |
