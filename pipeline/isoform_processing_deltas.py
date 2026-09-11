"""Disclose how much the isoform and processing choices could move amino acid fractions.

This tool decides nothing. Rule R1 fixes the canonical sequence; rule R2 fixes the
master range. This tool computes, from residue counts only (no masses):

  isoform deltas    : for every non-canonical isoform whose sequence was fetched,
                      fraction(isoform) - fraction(canonical) per amino acid
  processing deltas : fraction(master range) - fraction(full canonical) per amino acid

and the set-wide bound on each (the largest absolute delta seen for each amino acid,
and where it occurs). The bound is what the methods section reports.

Inputs   data/uniprot_verification.ini, config/segments.ini, config/accessions.ini
Outputs  outputs/isoform_deltas.tsv, outputs/isoform_bound.tsv,
         outputs/processing_deltas.tsv, outputs/processing_bound.tsv
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter

from pipeline import common

STAGE = "isoform_processing_deltas"
AA = common.AMINO_ACIDS


def fractions(seq: str) -> dict[str, float]:
    c = Counter(seq)
    n = len(seq)
    return {a: (c[a] / n if n else 0.0) for a in AA}


def delta(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    """a - b per amino acid."""
    return {x: a[x] - b[x] for x in AA}


def max_abs(d: dict[str, float]) -> tuple[float, str]:
    aa = max(AA, key=lambda x: abs(d[x]))
    return abs(d[aa]), aa


def master_sequence(acc: str, seq: str, segments) -> str:
    parts = []
    for section in segments.sections():
        if section.split(".", 1)[0] != acc:
            continue
        s = segments[section]
        if s.get("in_master_molecule") == "true":
            parts.append((int(s["start"]), int(s["end"])))
    if not parts:
        return ""
    return "".join(seq[a - 1:b] for a, b in sorted(parts))


def _fmt(x: float) -> str:
    return f"{x:.6f}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    log = common.make_logger(STAGE)
    ver = common.read_ini(common.VERIFICATION_INI)
    seg = common.read_ini(common.SEGMENTS_INI)
    acc_cfg = common.read_ini(common.ACCESSIONS_INI)

    iso_rows, proc_rows = [], []
    iso_bound = {a: (0.0, "") for a in AA}
    proc_bound = {a: (0.0, "") for a in AA}
    delta_cols = [f"delta_{a}" for a in AA]

    for acc in acc_cfg.sections():
        sec = ver[acc]
        seq = sec["sequence"]
        gene = sec.get("gene", "")
        tier = acc_cfg[acc]["tier"]
        base = fractions(seq)

        # isoforms
        i = 0
        while f"isoform.{i}.ids" in sec:
            p = f"isoform.{i}."
            iseq = sec.get(p + "sequence", "")
            if iseq:
                d = delta(fractions(iseq), base)
                m, at = max_abs(d)
                iso_rows.append([acc, gene, tier, sec.get(p + "fetched_id", ""), sec.get(p + "name", ""),
                                 sec.get(p + "status", ""), len(seq), len(iseq), _fmt(m), at,
                                 *[_fmt(d[a]) for a in AA]])
                for a in AA:
                    if abs(d[a]) > iso_bound[a][0]:
                        iso_bound[a] = (abs(d[a]), f"{acc} {sec.get(p + 'fetched_id', '')}")
            i += 1

        # processing
        ms = master_sequence(acc, seq, seg)
        if not ms:
            log.warning("%s: no master segment in segments.ini (open flag?) -- skipped", acc)
            continue
        d = delta(fractions(ms), base)
        m, at = max_abs(d)
        proc_rows.append([acc, gene, tier, len(seq), len(ms), len(seq) - len(ms), _fmt(m), at,
                          *[_fmt(d[a]) for a in AA]])
        for a in AA:
            if abs(d[a]) > proc_bound[a][0]:
                proc_bound[a] = (abs(d[a]), acc)

    common.write_tsv(common.OUTPUTS_DIR / "isoform_deltas.tsv",
                     ["accession", "gene", "tier", "isoform_id", "isoform_name", "status",
                      "length_canonical", "length_isoform", "max_abs_delta", "amino_acid_at_max", *delta_cols],
                     iso_rows)
    common.write_tsv(common.OUTPUTS_DIR / "isoform_bound.tsv",
                     ["amino_acid", "max_abs_delta_fraction", "where"],
                     [[a, _fmt(iso_bound[a][0]), iso_bound[a][1]] for a in AA])
    common.write_tsv(common.OUTPUTS_DIR / "processing_deltas.tsv",
                     ["accession", "gene", "tier", "length_full", "length_master", "n_residues_removed",
                      "max_abs_delta", "amino_acid_at_max", *delta_cols],
                     proc_rows)
    common.write_tsv(common.OUTPUTS_DIR / "processing_bound.tsv",
                     ["amino_acid", "max_abs_delta_fraction", "where"],
                     [[a, _fmt(proc_bound[a][0]), proc_bound[a][1]] for a in AA])

    overall_iso = max(iso_bound.values(), default=(0.0, ""))
    overall_proc = max(proc_bound.values(), default=(0.0, ""))
    log.info("isoform rows: %d; largest single-amino-acid isoform delta: %.6f (%s)",
             len(iso_rows), overall_iso[0], overall_iso[1])
    log.info("processing rows: %d; largest single-amino-acid processing delta: %.6f (%s)",
             len(proc_rows), overall_proc[0], overall_proc[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
