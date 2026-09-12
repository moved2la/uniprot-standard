# Provenance Standard

This repository produces a published standard. Every value, every list, and every
sentence in it must be traceable to one of exactly three kinds of origin. Anything
that cannot be is not admitted.

## The three admissible origins

1. **A decision by the author.** A scope choice, a threshold, an answer to a flag.
   Recorded in `docs/decisions.md` with a D-number, a date, and the reason. A decision
   may cite sources for its reason; the decision itself is the author's.

2. **A citation.** A published, retrievable source: a peer-reviewed paper, a public
   database (UniProt, Gene Ontology), a supplementary table, a government or consortium
   document. Recorded with enough detail to retrieve the exact value: source, version or
   release, table or field, date retrieved, and where applicable a checksum of the file.
   The citation sits on the same line as the value it supports.

3. **A computation.** Output of code in this repository, run on inputs of type 1 or 2,
   with the run logged. The code is inspectable, and tests exist for what it computes.

## What is not an origin

- **An AI assistant's recollection is never a source.** Not for a protein list, a
  category, a definition, a number, a citation, or a sentence of prose. An assistant may
  write code, may search for and surface published sources for the author to inspect,
  and may draft text that carries the citations it rests on. Its unsupported statements
  have the status of a colleague's hunch: worth checking, never worth recording.
- **"Well known" is not a source.** If it is well known, it is published; cite it.
- **A previous version of this repository is not a source.** Content is re-derived, not
  inherited, unless its own provenance was already one of the three kinds above.

## How the rule is enforced

- **Hand-written inputs are few and named.** The repository declares, in
  `docs/conventions.md`, which files are hand-written. Every other file is fetched or
  generated, and tests fail if a generated file differs from what its generator produces.
- **Unsettled means flagged, not decided.** When a rule in code cannot resolve a case
  from a database field, the code emits a flag. A flag is closed only by a logged
  decision (type 1) or a citation (type 2). Silent resolution — by code defaulting, or
  by a person quietly picking — is a defect.
- **Every generated file carries its own record.** A fetched or generated file states
  in its header what produced it: source URL, database release, file hash, query
  string, and retrieval time. What a rule could not settle is in `outputs/flags.tsv`.
  Those committed files, not logs, are the record a reviewer works from. Stage logs
  in `logs/` are for debugging only and are not committed.
- **Prose is held to the same rule as numbers.** A sentence in `docs/methods.md` either
  carries its citation or describes a computation performed here. Sentences that do
  neither are removed, not softened.

## Why this exists

The standard this repository replaces is a black box: a laboratory profile whose
sampling, hydrolysis, and provenance cannot be reconstructed by a third party. The
value of the replacement lies entirely in a reviewer being able to start from the
cited sources and the public databases, run the code, and arrive at the same numbers.
Any step that cannot be reproduced that way — including a step where an assistant
"knew" the answer — reintroduces the black box.
