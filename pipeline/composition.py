"""Per-protein amino acid composition from sequence, by arithmetic (Layer A).

For every accession in config/accessions.ini and for two segment sets:

  master     residues with in_master_molecule = true in config/segments.ini (the
             mature chain, rule R2)
  metabolic  every residue of the canonical sequence (the full gene product)

this stage computes

  count[aa]          number of residues of each of the twenty amino acids
                     (the stored ground truth, D7)
  residue_mass_sum   sum(count[aa] * residue_mass[aa])
  mw                 residue_mass_sum + water   -- molecular weight of the chain
  free_mass_sum      sum(count[aa] * free_mass[aa])
                     = mw + (n_residues - 1) * water    -- the water identity
  residue_frac[aa]   count[aa] * residue_mass[aa] / residue_mass_sum
  free_frac[aa]      count[aa] * free_mass[aa]    / free_mass_sum

residue_frac is the in-chain convention (sums to the protein); free_frac is the
free-amino-acid convention used by USDA food tables and by laboratory amino acid
analysis, and the one Match Rate consumes (D7). Both sum to 1 per row; g/100 g is
100 x the fraction. Masses come from data/pubchem/amino_acid_masses.ini (fetched from PubChem; the
letter-to-name table behind it is parsed from the IUPAC-IUBMB source); nothing here
is typed.

Disclosure (no decision attached): outputs/composition/processing_mass_deltas.tsv
gives, per accession, the mass removed by chain processing as a fraction of the
full gene product's mass, and free_frac(master) - free_frac(metabolic) per amino
acid; processing_mass_bound.tsv gives the set-wide maximum per amino acid.

Inputs   data/uniprot_sequences.ini, config/segments.ini, config/accessions.ini,
         data/pubchem/amino_acid_masses.ini
Outputs  outputs/composition/amino_acid_composition_per_protein.tsv  (one row per accession per segment set)
         outputs/composition/processing_mass_deltas.tsv
         outputs/composition/processing_mass_bound.tsv
         outputs/composition/composition_summary.ini  (counts, identities checked, bounds)

--check rebuilds everything in memory and fails if the files on disk differ.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from pipeline import common

STAGE = "composition"
AA = common.AMINO_ACIDS
SEGMENT_SETS = ("master", "metabolic")
FRAC_DIGITS = 10   # fractions are written with this many significant digits
MASS_DIGITS = 4    # sums and molecular weights, decimal places


# ----------------------------------------------------------------------------
# Masses
# ----------------------------------------------------------------------------

def load_masses(path: Path = common.AMINO_ACID_MASSES_INI) -> dict:
    """Return {'free': {aa: g/mol}, 'residue': {aa: g/mol}, 'water': g/mol, 'hydrogen': g/mol}."""
    cp = common.read_ini(path)
    missing = [a for a in AA if not cp.has_section(a)]
    if missing:
        raise ValueError(f"{path}: no section for {missing}")
    for name in ("water", "hydrogen"):
        if not cp.has_section(name):
            raise ValueError(f"{path}: no [{name}] section")
    return {
        "free": {a: float(cp[a]["free_mass"]) for a in AA},
        "residue": {a: float(cp[a]["residue_mass"]) for a in AA},
        "water": float(cp["water"]["molecular_weight"]),
        "hydrogen": float(cp["hydrogen"]["molecular_weight"]),
    }


# ----------------------------------------------------------------------------
# Sequence slicing
# ----------------------------------------------------------------------------

def segment_ranges(segments, acc: str, master_only: bool) -> list[tuple[int, int]]:
    """1-indexed inclusive ranges for an accession from segments.ini."""
    out = []
    for section in segments.sections():
        if section.split(".", 1)[0] != acc:
            continue
        s = segments[section]
        if master_only and s.get("in_master_molecule") != "true":
            continue
        out.append((int(s["start"]), int(s["end"])))
    return sorted(out)


def slice_sequence(seq: str, ranges: list[tuple[int, int]]) -> str:
    return "".join(seq[a - 1:b] for a, b in ranges)


# ----------------------------------------------------------------------------
# Composition of one sequence
# ----------------------------------------------------------------------------

def compose(seq: str, masses: dict) -> dict:
    """Residue counts, both mass vectors, both fraction vectors, and the identities."""
    bad = set(seq) - set(AA)
    if bad:
        raise ValueError(f"non-standard letters in sequence: {sorted(bad)}")
    counts = Counter(seq)
    n = len(seq)
    residue_vec = {a: counts[a] * masses["residue"][a] for a in AA}
    free_vec = {a: counts[a] * masses["free"][a] for a in AA}
    residue_sum = sum(residue_vec.values())
    free_sum = sum(free_vec.values())
    mw = residue_sum + masses["water"] if n else 0.0
    return {
        "n_residues": n,
        "count": {a: counts[a] for a in AA},
        "residue_mass_sum": residue_sum,
        "free_mass_sum": free_sum,
        "mw": mw,
        "residue_frac": {a: (residue_vec[a] / residue_sum if residue_sum else 0.0) for a in AA},
        "free_frac": {a: (free_vec[a] / free_sum if free_sum else 0.0) for a in AA},
    }


def water_identity_error(c: dict, water: float) -> float:
    """|free_mass_sum - (mw + (n - 1) * water)|; zero up to floating-point rounding."""
    if c["n_residues"] == 0:
        return 0.0
    return abs(c["free_mass_sum"] - (c["mw"] + (c["n_residues"] - 1) * water))


# ----------------------------------------------------------------------------
# Formatting (deterministic across machines: fixed digit counts)
# ----------------------------------------------------------------------------

def _f(x: float) -> str:
    return f"{x:.{FRAC_DIGITS}g}"


def _m(x: float) -> str:
    return f"{x:.{MASS_DIGITS}f}"


ALL_COLUMNS = (["accession", "gene", "tier", "segment_set", "n_residues", "residue_mass_sum", "mw", "free_mass_sum"]
               + [f"count_{a}" for a in AA] + [f"residue_frac_{a}" for a in AA] + [f"free_frac_{a}" for a in AA])

DELTA_COLUMNS = (["accession", "gene", "tier", "n_residues_full", "n_residues_master", "n_residues_removed",
                  "mw_full", "mw_master", "mass_removed_fraction_of_full", "max_abs_delta_free_frac", "amino_acid_at_max"]
                 + [f"delta_free_frac_{a}" for a in AA])


# ----------------------------------------------------------------------------
# Build
# ----------------------------------------------------------------------------

def build(accessions, sequences, segments, masses, log) -> tuple[list[list], list[list], list[list], dict]:
    """Return (all_rows, delta_rows, bound_rows, summary)."""
    all_rows: list[list] = []
    delta_rows: list[list] = []
    bound = {a: (0.0, "") for a in AA}
    max_identity_error = 0.0
    n_master_eq_full = 0
    skipped = []

    for acc in accessions.sections():
        entry = accessions[acc]
        if entry.get("flag_open") == "true":
            skipped.append(acc)
            log.warning("%s: flag_open = true in accessions.ini -- skipped", acc)
            continue
        if not sequences.has_section(acc):
            raise RuntimeError(f"{acc} is in accessions.ini but not in {common.SEQUENCES_INI}")
        seq = sequences[acc]["sequence"]
        gene, tier = entry.get("gene", ""), entry.get("tier", "")
        master_ranges = segment_ranges(segments, acc, master_only=True)
        if not master_ranges:
            raise RuntimeError(f"{acc}: no in_master_molecule = true segment in {common.SEGMENTS_INI}")
        comp = {
            "master": compose(slice_sequence(seq, master_ranges), masses),
            "metabolic": compose(seq, masses),
        }
        for segment_set in SEGMENT_SETS:
            c = comp[segment_set]
            max_identity_error = max(max_identity_error, water_identity_error(c, masses["water"]))
            all_rows.append([acc, gene, tier, segment_set, c["n_residues"],
                             _m(c["residue_mass_sum"]), _m(c["mw"]), _m(c["free_mass_sum"]),
                             *[c["count"][a] for a in AA],
                             *[_f(c["residue_frac"][a]) for a in AA],
                             *[_f(c["free_frac"][a]) for a in AA]])
        m, full = comp["master"], comp["metabolic"]
        if m["n_residues"] == full["n_residues"]:
            n_master_eq_full += 1
        d = {a: m["free_frac"][a] - full["free_frac"][a] for a in AA}
        aa_max = max(AA, key=lambda a: abs(d[a]))
        removed_fraction = (full["mw"] - m["mw"]) / full["mw"] if full["mw"] else 0.0
        delta_rows.append([acc, gene, tier, full["n_residues"], m["n_residues"],
                           full["n_residues"] - m["n_residues"], _m(full["mw"]), _m(m["mw"]),
                           _f(removed_fraction), _f(abs(d[aa_max])), aa_max, *[_f(d[a]) for a in AA]])
        for a in AA:
            if abs(d[a]) > bound[a][0]:
                bound[a] = (abs(d[a]), acc)

    bound_rows = [[a, _f(bound[a][0]), bound[a][1]] for a in AA]
    largest = max(AA, key=lambda a: bound[a][0])
    largest_removed = max(delta_rows, key=lambda r: float(r[8]), default=None)
    summary = {
        "accessions": str(len(accessions.sections())),
        "accessions_skipped_open_flag": ";".join(skipped),
        "rows_all": str(len(all_rows)),
        "rows_processing_deltas": str(len(delta_rows)),
        "entries_master_equals_full": str(n_master_eq_full),
        "entries_with_processing": str(len(delta_rows) - n_master_eq_full),
        "water_identity_max_abs_error_g_per_mol": f"{max_identity_error:.3e}",
        "largest_processing_delta_free_frac": _f(bound[largest][0]),
        "largest_processing_delta_amino_acid": largest,
        "largest_processing_delta_accession": bound[largest][1],
        "largest_mass_removed_fraction_of_full": largest_removed[8] if largest_removed else "",
        "largest_mass_removed_accession": largest_removed[0] if largest_removed else "",
    }
    return all_rows, delta_rows, bound_rows, summary


def _tsv_text(header: list[str], rows: list[list]) -> str:
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join("" if v is None else str(v) for v in row))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare a fresh build with outputs on disk; write nothing")
    args = ap.parse_args(argv)
    log = common.make_logger(STAGE + ("_check" if args.check else ""))

    accessions = common.read_ini(common.ACCESSIONS_INI)
    sequences = common.read_ini(common.SEQUENCES_INI)
    segments = common.read_ini(common.SEGMENTS_INI)
    masses = load_masses()
    masses_cp = common.read_ini(common.AMINO_ACID_MASSES_INI)
    masses_header = _ini_header(common.AMINO_ACID_MASSES_INI)

    all_rows, delta_rows, bound_rows, summary = build(accessions, sequences, segments, masses, log)

    header_note = [
        "GENERATED by composition.py. DO NOT EDIT BY HAND.",
        f"masses            = {_rel(common.AMINO_ACID_MASSES_INI)} (PubChem, fetched {masses_header.get('fetched', '')}; MolecularWeight served to {_first_token(masses_header.get('served_decimals', '?'))} decimals)",
        f"water             = {masses['water']} g/mol (PubChem CID {masses_cp['water']['cid']})",
        f"sequences         = {_rel(common.SEQUENCES_INI)} (UniProt release {_ini_header(common.SEQUENCES_INI).get('uniprot_release', '')})",
        f"segments          = {_rel(common.SEGMENTS_INI)}",
        f"fraction_digits   = {FRAC_DIGITS} significant; mass_decimals = {MASS_DIGITS}",
    ]
    summary_cp = common.new_ini()
    summary_cp["composition"] = summary

    files = {
        common.COMPOSITION_TSV: _tsv_text(ALL_COLUMNS, all_rows),
        common.COMPOSITION_DIR / "processing_mass_deltas.tsv": _tsv_text(DELTA_COLUMNS, delta_rows),
        common.COMPOSITION_DIR / "processing_mass_bound.tsv":
            _tsv_text(["amino_acid", "max_abs_delta_free_frac", "where"], bound_rows),
        common.COMPOSITION_SUMMARY_INI: common.render_ini(summary_cp, header_note),
    }

    if args.check:
        ok = True
        for path, text in files.items():
            if not path.exists() or path.read_text(encoding="utf-8") != text:
                ok = False
                log.error("STALE: %s differs from a fresh build", path)
        log.info("check: %s", "composition outputs are current" if ok else "composition outputs are STALE")
        return 0 if ok else 1

    common.COMPOSITION_DIR.mkdir(parents=True, exist_ok=True)
    for path, text in files.items():
        path.write_text(text, encoding="utf-8", newline="\n")
    log.info("wrote %d rows to amino_acid_composition_per_protein.tsv, %d processing rows; water identity max error %s g/mol; "
             "largest processing delta %s (%s, %s); largest mass removed %s (%s)",
             len(all_rows), len(delta_rows), summary["water_identity_max_abs_error_g_per_mol"],
             summary["largest_processing_delta_free_frac"], summary["largest_processing_delta_amino_acid"],
             summary["largest_processing_delta_accession"], summary["largest_mass_removed_fraction_of_full"],
             summary["largest_mass_removed_accession"])
    return 0


def _rel(path: Path) -> str:
    """Repository-relative path with forward slashes on every platform (byte-identical headers)."""
    return path.relative_to(common.REPO_ROOT).as_posix()


def _first_token(value: str) -> str:
    return value.split()[0] if value.split() else value


def _ini_header(path: Path) -> dict[str, str]:
    """Parse '# key = value' comment lines at the top of a generated .ini."""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("#"):
            break
        body = line.lstrip("#").strip()
        if "=" in body:
            k, _, v = body.partition("=")
            out[k.strip()] = v.strip()
    return out


if __name__ == "__main__":
    sys.exit(main())
