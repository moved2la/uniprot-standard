# Step 4b report — the bound-metabolite adjustment of the skeletal muscle standard

**Thread:** 2026-09-16 to 2026-09-17. **Decisions:** D72–D81. **Deliveries:** 1–9.
**Result:** `outputs/standard/_calculated_amino_acid_standard_with_bound_pools.tsv` — the muscle standard with the
histidine held as carnosine folded in, beside the unchanged protein-only standard. Version label
`ref = skeletal muscle protein + carnosine (D73–D79)`.

## What the step set out to do

The protein-only standard (Step 4) counts the amino acids held as protein. Muscle also holds one amino acid in a
stable non-protein form — histidine, as carnosine — which the body had to eat to put there. D68 (plan R7) said such
bound pools belong in the reference; this step put the first one in, with cited numbers replacing the planning
summary's guesses (a ~+11 % turnover-weighted figure built on two estimated rates).

## Procedure

1. **Sources found and read** (2026-09-16/17): single-fiber carnosine in ATPase-typed fibres of human vastus
   lateralis (Harris, Dunnett & Greenhaff 1998 — paywalled; its abstract carries every value used and the abstract
   page saved as a PDF is the file on disk); the statement that carnosine is the only histidine-containing dipeptide
   in human muscle (Harris et al. 2012, open access); washout half-lives (Baguet et al. 2009; Yamaguchi et al. 2021,
   accepted manuscript); muscle protein and water per kilogram (Mingrone et al. 2001, direct chemical analysis);
   ICRP Publication 89 as the reference-body document, read and found to give elemental composition only.
2. **Decisions before build:** inventory frame in the standard, resupply frame a sensitivity (D73); per-fiber-type
   single-fiber values in the fiber-type columns, type II applied to IIa and IIx alike (D74); carnosine only (D75);
   protein-only table untouched, adjusted table separate (D76); Mingrone 2001 as the denominator, a measured value
   chosen over one computed from ICRP 89's nitrogen (D77).
3. **Build:** stage `bound_pools` between `aggregate` and `uncertainty` in `python run.py standard` (rule A10).
   The pool's amino acid, stoichiometry, every value, every basis, and which columns a value fills come from
   `config/bound_metabolite_pools.ini`; the code names nothing. Twelve synthetic tests.
4. **Green light** (D78): the stage reads `data/literature/manifest.ini` and stops if any source the config names
   has no hashed file. Every 4b source was hashed by the fetch of 2026-09-17T01:28Z.
5. **Transcription** (D79): values read from the PDFs on disk by the assistant, each with page and table or
   paragraph, audited by the author. The muscle-protein turnover rate for the resupply sensitivity was not sourced
   (author's call); that sensitivity is reported as not computed.
6. **Run** 2026-09-17T02:12Z: 121 tests passed; excerpts uploaded and read.

## Results

| Column | Histidine, protein-only | with carnosine | pool ÷ protein histidine |
|---|---|---|---|
| fiber type I | 2.206 % | 2.416 % | 0.097 |
| fiber type IIa | 2.250 % | 2.712 % | 0.211 |
| fiber type IIx | 2.229 % | 2.691 % | 0.213 |
| standard (D72 mix) | 2.227 % | 2.568 % | — |

Per kg dry muscle: 654 g protein (155 g/kg wet ÷ (1 − 0.763)) holding 16.7–17.1 g histidine, against 1.63 g (type I)
and 3.60 g (type II) held as carnosine. Every other amino acid moves down by 0.35 % of its value. The between-fiber
spread of the pool (log-normal on 10.5 ± 7.6 and 23.2 ± 17.8; the ± kind is not stated in the abstract) puts
histidine at 2.25–2.81 % in fiber type I and 2.35–3.63 % in IIa at the 2.5 / 97.5 percentiles — the scatter of
single fibers around the mean, not the uncertainty of the mean.

Against the earlier planning guess: the inventory frame from cited numbers gives +15 % on the mixed standard where
the turnover-frame guess had given +11 %. Renormalised to the sixteen amino acids Gorissen et al. 2018 measure, the
standard's histidine moves from 2.50 to 2.88 against Gorissen's 4.61; carnosine accounts for a little under a fifth of
that difference. What the rest is — the blood category, the biopsy's residual blood — is a question for Steps 7a and
9, reported as measurements, not verdicts.

## Method lessons (for the refactor, Step 7b — all pinned)

- **Fetch-by-code retired (D80).** Nine sources: three fetched, five obtained by hand, one dropped. The
  requirement to re-fetch static files every run cost hours across blocked publishers, HTML served as PDFs, 403s,
  a file locked by a PDF viewer, and naming rules; it bought nothing the hash on disk does not. The replacement is
  a stage that hashes what is in `data/literature/` and writes the manifest. Until then: F3 (unchanged files are not
  rewritten; a locked file is a flag), F4b (a failed re-fetch of a verified file is not a flag), F7b (a stored file
  always carries a real extension), F2b (a file keeps its first retrieval time).
- **The staleness trap (D81).** Every fetch stamped every manifest entry with "now", and the generated weight files
  copied that stamp, so a fetch for an unrelated step made the weights read as stale and a test fail with nothing
  changed. Retrieval times are out of generated config. The general fix is 7b's: the source list leaves
  `mass_fraction_decisions.ini`, and each stage hashes only the config sections it reads.
- **A run-order cheat sheet** (`docs/run_order.md`): what depends on what, and which command to run after which
  edit; kept separate from the README.
- **Transcription (D79):** reading values off a PDF is the assistant's job with the page beside each value; the
  author audits; the calculation stays deterministic. The author-only rule was an over-reach.
- **Measured over computed:** a literature value is a harvested measurement; a number derived from a source's
  tables under assumptions is a last resort, presented as such.
- **Money:** the budget does not buy papers. Where an abstract carries the values, the abstract page is the record.
- **Delivery numbering:** every tarball is an increment, patches included.

## Carried forward

- The resupply-frame sensitivity needs a cited muscle-protein turnover range (the primary studies at each end of the
  four-fold published range); not pursued this step.
- Harris 2012 "p. 6" for the sole-dipeptide statement: read from the open HTML, which carries no page numbers;
  to confirm in the PDF on disk.
- A second single-fiber source (Tallon et al. 2007; Hill et al. 2007 reports type I 17.8 / type II 29.6 mmol/kg dm
  before supplementation) would cost a paper each; the columns rest on one study of four men, stated.
- Creatine and glutathione (non-EAA pools) are a later extension under D68.
