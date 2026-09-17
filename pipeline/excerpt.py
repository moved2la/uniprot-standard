#!/usr/bin/env python3
"""
excerpt.py — tooling, not a pipeline stage: bounded excerpts of the large generated files.

The generated tables run to thousands of rows; a review in chat needs their headers and a
few hundred well-chosen rows, not the bodies. This tool writes, for a fixed list of files,
either the whole file (small ones) or a labelled cut of it — the top N rows by a column,
the first N rows of an already-ranked table, the rows matching a filter, or the rows whose
key appears in an earlier cut — keeping every `#` header line so the provenance travels
with the rows. Nothing is computed, changed, or interpreted; every excerpt states which
file it came from, the rule that selected its rows, and how many of the file's rows it holds.

    python run.py excerpt

Writes  excerpts/<UTC stamp>/<same relative path>[__<cut>].<ext>
        excerpts/<UTC stamp>/INDEX.md                 what each excerpt is
        excerpts/excerpts_<UTC stamp>.tar.gz          the folder, for upload
        logs/excerpt_<UTC>.log

`excerpts/` is not committed and is not part of the record (like docs/handoffs/). A file
that is missing on this machine is listed as missing in INDEX.md and the log; it does not
stop the run. The list of files and cuts is SPEC below — edit it there.
"""

from __future__ import annotations

import sys
import tarfile
from pathlib import Path

from pipeline import common

EXCERPTS_DIR = common.REPO_ROOT / "excerpts"

# ----------------------------------------------------------------------------
# What to excerpt. Each entry: (relative path, [cuts]). A cut is a dict:
#   {"kind": "whole"}
#   {"kind": "head",   "n": N}                                  first N rows of an already-ranked table
#   {"kind": "top",    "n": N, "by": <column>}                  N largest by a numeric column
#   {"kind": "equals", "column": c, "value": v}                 rows where column == value
#   {"kind": "contains", "column": c, "value": v}               rows where value is a substring of the cell
#   {"kind": "in", "column": c, "from": <path>, "from_cut": <name>, "from_column": c2}
#                                                              rows whose column is in the set of
#                                                              values of c2 in an earlier cut
# Every cut has a "name" used in the output filename; "whole" needs none.
# ----------------------------------------------------------------------------

MF = "outputs/mass_fractions"
CO = "outputs/composition"

SPEC: list[tuple[str, list[dict]]] = [
    # --- one-screen summaries and provenance (whole) ---
    (f"{MF}/mass_fractions_summary.ini", [{"kind": "whole"}]),
    (f"{CO}/composition_summary.ini", [{"kind": "whole"}]),
    (f"{CO}/ptm_summary.ini", [{"kind": "whole"}]),
    ("outputs/digest/digest_summary.ini", [{"kind": "whole"}]),
    ("outputs/literature_inventory/literature_inventory_summary.ini", [{"kind": "whole"}]),
    ("data/literature/manifest.ini", [{"kind": "whole"}]),
    ("data/gene-ontology/source.ini", [{"kind": "whole"}]),
    ("data/iupac/source.ini", [{"kind": "whole"}]),
    ("data/iupac/amino_acid_symbols.ini", [{"kind": "whole"}]),
    ("data/pubchem/amino_acid_masses.ini", [{"kind": "whole"}]),
    ("outputs/flags.tsv", [{"kind": "whole"}]),
    # --- the weights (generated config, D61) ---
    ("config/mass_fractions_per_entry.tsv", [
        {"kind": "top", "name": "top100_by_w_I_combined", "n": 100, "by": "w_I_combined"},
        {"kind": "equals", "name": "no_dataset_row", "column": "match_rule", "value": "none"},
        {"kind": "contains", "name": "shared_split", "column": "match_rule", "value": "shared_split"},
        {"kind": "contains", "name": "via_mapping", "column": "match_rule", "value": "via_mapping"},
    ]),
    # --- Layer A rows for the same top-100 accessions ---
    (f"{CO}/amino_acid_composition_per_protein.tsv", [
        {"kind": "in", "name": "master_rows_of_top100", "column": "accession",
         "from": "config/mass_fractions_per_entry.tsv", "from_cut": "top100_by_w_I_combined", "from_column": "accession",
         "also_equals": ("segment_set", "master")},
    ]),
    # --- ranked tables (already ordered; first rows) ---
    (f"{MF}/combined_entries_ranked.tsv", [{"kind": "head", "name": "first50", "n": 50}]),
    (f"{MF}/tier1_entries_ranked.tsv", [{"kind": "head", "name": "first30", "n": 30}]),
    (f"{MF}/tier2_entries_ranked.tsv", [{"kind": "head", "name": "first30", "n": 30}]),
    (f"{MF}/dataset_rows_outside_pool.tsv", [{"kind": "head", "name": "first30", "n": 30}]),
    # --- bounds and checks (small; whole) ---
    (f"{MF}/weighted_bounds.tsv", [{"kind": "whole"}]),
    (f"{MF}/family_bounds.tsv", [{"kind": "whole"}]),
    (f"{MF}/band_families.tsv", [{"kind": "whole"}]),
    (f"{MF}/classical_check_carroll_2004.tsv", [{"kind": "whole"}]),
    (f"{MF}/excluded_entries_mass_share.tsv", [{"kind": "whole"}]),
    (f"{MF}/branch_exclusive_mass.tsv", [{"kind": "whole"}]),
    (f"{MF}/pool_entries_without_dataset_row.tsv", [{"kind": "whole"}]),
    (f"{MF}/momenzadeh_2023_myh_fractions_ibaq_vs_lfq.tsv", [{"kind": "whole"}]),
    (f"{MF}/ratio_check_moreno-justicia_2025.tsv", [{"kind": "head", "name": "first20", "n": 20}]),
    (f"{MF}/ratio_check_deshmukh_2021.tsv", [{"kind": "head", "name": "first20", "n": 20}]),
    ("outputs/isoform_bound.tsv", [{"kind": "whole"}]),
    ("outputs/processing_bound.tsv", [{"kind": "whole"}]),
    (f"{CO}/processing_mass_bound.tsv", [{"kind": "whole"}]),
    ("outputs/excluded_non_standard_alphabet.tsv", [{"kind": "whole"}]),
    ("outputs/chain_positions_resolved.tsv", [{"kind": "whole"}]),
    ("outputs/digest/density_ranked.tsv", [{"kind": "head", "name": "first20", "n": 20}]),
    ("outputs/digest/families.tsv", [{"kind": "head", "name": "first40", "n": 40}]),
    # --- the standard (small tables whole; the per-entry regime table cut) ---
    ("outputs/standard/_calculated_amino_acid_standard.tsv", [{"kind": "whole"}]),
    ("outputs/standard/_calculated_amino_acid_standard_residue_convention.tsv", [{"kind": "whole"}]),
    ("outputs/standard/standard_summary.ini", [{"kind": "whole"}]),
    ("outputs/standard/_calculated_amino_acid_standard_with_non_protein_metabolite_pools.tsv", [{"kind": "whole"}]),
    ("outputs/standard/non_protein_metabolite_pools_summary.ini", [{"kind": "whole"}]),
    ("outputs/standard/non_protein_metabolite_pool_adjustment_per_amino_acid.tsv", [{"kind": "whole"}]),
    ("outputs/standard/non_protein_metabolite_pool_amounts_per_kg_muscle.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sensitivity_non_protein_metabolite_pool_turnover_frame.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sensitivity_non_protein_metabolite_pool_basis.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sensitivity_non_protein_metabolite_pool_sex.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sensitivity_non_protein_metabolite_pool_spread.tsv", [{"kind": "whole"}]),
    ("outputs/standard/eaa_subset_with_non_protein_metabolite_pools.tsv", [{"kind": "whole"}]),
    ("outputs/standard/stress_summary.ini", [{"kind": "whole"}]),
    ("outputs/standard/stress_summary_per_scenario.tsv", [{"kind": "whole"}]),
    ("outputs/standard/stress_influence_per_entry.tsv", [{"kind": "head", "name": "first60", "n": 60}]),
    ("outputs/standard/composition_distance_top_entries.tsv", [{"kind": "whole"}]),
    ("outputs/standard/uncertainty_summary.ini", [{"kind": "whole"}]),
    ("outputs/standard/amino_acid_profiles.tsv", [{"kind": "whole"}]),
    ("outputs/standard/amino_acid_g_per_100g_protein_free.tsv", [{"kind": "whole"}]),
    ("outputs/standard/amino_acid_g_per_100g_protein_residue.tsv", [{"kind": "whole"}]),
    ("outputs/standard/profile_differences.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sensitivity_mhc_actin_profiles.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sensitivity_mhc_actin_spread.tsv", [{"kind": "whole"}]),
    ("outputs/standard/bounds_after_weighting.tsv", [{"kind": "whole"}]),
    ("outputs/standard/completeness_sensitivity.tsv", [{"kind": "whole"}]),
    # --- the USDA food tables (the big one cut; the rest whole) ---
    ("outputs/usda/usda_archive_summary.ini", [{"kind": "whole"}]),
    ("outputs/usda/usda_nutrient_map.tsv", [{"kind": "whole"}]),
    ("outputs/usda/amino_acids_per_food.tsv", [{"kind": "head", "name": "first60", "n": 60}]),
    ("outputs/usda/foods_without_amino_acids.tsv", [{"kind": "head", "name": "first40", "n": 40}]),
    # --- the comparison with the measured composition ---
    ("outputs/comparison/gorissen_2018_human_muscle.tsv", [{"kind": "whole"}]),
    ("outputs/comparison/calculated_vs_gorissen_2018.tsv", [{"kind": "whole"}]),
    ("outputs/comparison/comparison_summary.ini", [{"kind": "whole"}]),
    ("outputs/standard/eaa_subset.tsv", [{"kind": "whole"}]),
    ("outputs/standard/uncertainty_intervals.tsv", [{"kind": "whole"}]),
    ("outputs/standard/uncertainty_per_amino_acid.tsv", [{"kind": "whole"}]),
    ("outputs/standard/sd_to_median_ratio.tsv", [
        {"kind": "top", "name": "top100_by_w_combined", "n": 100, "by": "w_combined"},
        {"kind": "top", "name": "top100_by_sd_over_median", "n": 100, "by": "sd_over_median"},
    ]),
]


# ----------------------------------------------------------------------------

def read_table(path: Path) -> tuple[list[str], list[str], list[list[str]]]:
    """(comment lines, header columns, rows) of a TSV whose first lines may be '#' comments."""
    comments, header, rows = [], None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if header is None and line.startswith("#"):
                comments.append(line)
                continue
            if not line:
                continue
            if header is None:
                header = line.split("\t")
                continue
            rows.append(line.split("\t"))
    return comments, header or [], rows


def as_float(x: str) -> float:
    try:
        return float(x)
    except ValueError:
        return float("-inf")


def apply_cut(cut: dict, header: list[str], rows: list[list[str]], collected: dict) -> tuple[list[list[str]], str]:
    """Return (rows kept, one-line description of the rule)."""
    kind = cut["kind"]
    if kind == "head":
        return rows[: cut["n"]], f"first {cut['n']} rows as ordered in the file"
    if kind == "top":
        i = header.index(cut["by"])
        return sorted(rows, key=lambda r: -as_float(r[i]))[: cut["n"]], f"the {cut['n']} largest by column {cut['by']} (descending)"
    if kind == "equals":
        i = header.index(cut["column"])
        return [r for r in rows if r[i] == cut["value"]], f"rows where {cut['column']} == {cut['value']!r}"
    if kind == "contains":
        i = header.index(cut["column"])
        return [r for r in rows if cut["value"] in r[i]], f"rows where {cut['column']} contains {cut['value']!r}"
    if kind == "in":
        keys = collected.get((cut["from"], cut["from_cut"], cut["from_column"]))
        if keys is None:
            raise KeyError(f"cut {cut['name']}: its source cut {cut['from']}::{cut['from_cut']} did not run first")
        i = header.index(cut["column"])
        kept = [r for r in rows if r[i] in keys]
        rule = f"rows whose {cut['column']} is one of the {len(keys)} values of {cut['from_column']} in excerpt {cut['from']}::{cut['from_cut']}"
        if "also_equals" in cut:
            c, v = cut["also_equals"]
            j = header.index(c)
            kept = [r for r in kept if r[j] == v]
            rule += f", and {c} == {v!r}"
        return kept, rule
    raise ValueError(f"unknown cut kind {kind!r}")


def main(argv: list[str] | None = None) -> int:
    log = common.make_logger("excerpt")
    stamp = common.timestamp()
    out_root = EXCERPTS_DIR / stamp
    out_root.mkdir(parents=True, exist_ok=True)
    index = ["# Excerpts", "",
             f"Written {common.iso_now()} by `pipeline/excerpt.py` from the files on this machine. "
             "Each excerpt keeps its source file's `#` header lines and adds one line naming the rule that "
             "selected its rows. Nothing is computed or changed here.", "",
             "| Excerpt | Source file | Rule | Rows kept / in file |", "|---|---|---|---|"]
    collected: dict[tuple, set[str]] = {}
    written = 0
    for rel, cuts in SPEC:
        src = common.REPO_ROOT / rel
        if not src.exists():
            log.warning("missing on this machine: %s", rel)
            index.append(f"| — | `{rel}` | missing on this machine | — |")
            continue
        for cut in cuts:
            if cut["kind"] == "whole":
                dest = out_root / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(src.read_bytes())
                n = sum(1 for _ in open(src, encoding="utf-8", errors="replace"))
                index.append(f"| `{rel}` | `{rel}` | whole file | {n} lines |")
                log.info("whole  %s (%d lines)", rel, n)
                written += 1
                continue
            comments, header, rows = read_table(src)
            kept, rule = apply_cut(cut, header, rows, collected)
            name = cut["name"]
            p = Path(rel)
            dest = out_root / p.parent / f"{p.stem}__{name}{p.suffix}"
            dest.parent.mkdir(parents=True, exist_ok=True)
            lines = [f"# EXCERPT of {rel}: {rule}; {len(kept)} of {len(rows)} rows; written {common.iso_now()} by pipeline/excerpt.py"]
            lines += comments + ["\t".join(header)] + ["\t".join(r) for r in kept]
            dest.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
            for col in header:                       # every column of every cut is collectable by later cuts
                i = header.index(col)
                collected[(rel, name, col)] = {r[i] for r in kept}
            index.append(f"| `{dest.relative_to(out_root).as_posix()}` | `{rel}` | {rule} | {len(kept)} / {len(rows)} |")
            log.info("cut    %s :: %s -> %d of %d rows", rel, name, len(kept), len(rows))
            written += 1
    (out_root / "INDEX.md").write_text("\n".join(index) + "\n", encoding="utf-8", newline="\n")
    tar_path = EXCERPTS_DIR / f"excerpts_{stamp}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tf:
        tf.add(out_root, arcname=stamp)
    log.info("%d excerpts in %s; tarball %s (%d bytes)", written, out_root.relative_to(common.REPO_ROOT).as_posix(),
             tar_path.relative_to(common.REPO_ROOT).as_posix(), tar_path.stat().st_size)
    print(f"\nUPLOAD: {tar_path.relative_to(common.REPO_ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
