#!/usr/bin/env python3
"""
mass_fractions.py — offline stage: Layer B weights from the primary dataset.

Reads
  config/literature_sources.ini           the [source.*] with role = primary (Step 7a)
  config/mass_fraction_decisions.ini      its [columns.*] map;
                                          the method-citation source with two tables (iBAQ vs LFQ)
  data/literature/manifest.ini            path + sha256 of every literature file (B0)
  config/accessions.ini                   the pool: accession, gene, tier
  outputs/composition/amino_acid_composition_per_protein.tsv   mw on segment_set = master (B2)
  data/tier<N>_pool.tsv, data/gene-ontology/tier<N>_subtree.tsv   branch table
  outputs/isoform_deltas.tsv, outputs/composition/processing_mass_deltas.tsv,
  outputs/composition/ptm_mass_deltas.tsv                       weighted bounds
  outputs/digest/shared_pairs.tsv                                gel-band families (D52)
  config/carroll_classical_fractionation.ini                    gel values (D50)

Writes (outputs/mass_fractions/ unless stated; every file with a header naming inputs and hashes)
  config/mass_fractions/{I,IIa,IIx}.ini  generated config; one section per accession: within-tier
                                       weight (D31) and combined weight (D58), each with low / high,
                                       dataset values, match rule, rows, source line — unchanged shape
  config/mass_fractions_per_entry.tsv  the same weights as ONE TABLE (D61, an added format; nothing
                                       replaced): one row per pool accession, every fiber type and
                                       weight side by side, the constant provenance once in the header
  weights_per_pool_entry.tsv           one row per pool accession: what the dataset says about
                                       it, how the row was matched, and its weight per tier and
                                       fiber type — the full table
  combined_entries_ranked.tsv          D58: every entry of every tier under one denominator per
                                       fiber type — the combined (contractile + support) standard
  excluded_entries_mass_share.tsv      D59: the R5-excluded entries and the share they would have held
  tier<N>_entries_ranked.tsv           the entries of one tier, ordered by type-I weight, with
                                       rank per fiber type and cumulative share down each order —
                                       the table to read, quote, and plot
  pool_entries_without_dataset_row.tsv pool accessions with no Dataset 1 row (w = 0)
  dataset_rows_outside_pool.tsv        Dataset 1 genes matching no pool entry, ranked by median
                                       value — the completeness check and the Tier 2 evidence
  shared_gene_rows.tsv                 every ';' row touching the pool with its D48 outcome
  branch_exclusive_mass.tsv            per ontology term: sum of w over entries reached only
                                       through that term's branch
  weighted_bounds.tsv                  isoform / processing / PTM bounds x w, largest per amino
                                       acid; sum of w over glycosylated entries
  classical_check_carroll_2004.tsv     MHC:actin as measured by each method — Carroll 2004's gel,
                                       iBAQ x MW (ours), and the intensity share in the second
                                       mass-spec source that carries accessions — and their
                                       quotients; no attribution (D54)
  band_families.tsv                    the entries each band family contains, at cutoffs 1 and 2,
                                       with their dataset values and weights per fiber type
  family_bounds.tsv                    D55: per shared-peptide family, the composition spread across
                                       members x the family's weight, per amino acid
  ratio_check_<source>.tsv             per pool entry, the slow/fast ratio in each cross-check source
                                       against the primary dataset's I/IIa and I/IIx (D34, D36)
  <source>_myh_fractions_ibaq_vs_lfq.tsv   per gene of the two tables, mean fraction per fiber under
                                       iBAQ and MaxLFQ, over the fibers present in both tables
  mass_fractions_summary.ini           counts per file; per tier and fiber type: entries with
                                       mass, share of the ten largest, largest entry, maximum
                                       weighted bound of each kind; the gel-check factors
  logs/mass_fractions_<UTC>.log

Rules (docs/conventions.md, "Mass-fraction rules")
  B0 every file read is re-hashed against the manifest; mismatch stops the run.
  B1 identity by exact gene symbol after stripping whitespace, against accessions.ini.
     Several rows per gene: summed (D47). Multi-gene cells: split (D48).
  B2 w_i = v_i * MW_i / sum_tier(v_j * MW_j), within tier (D31).
  B3 fiber types are the dataset's own columns; never re-derived.
  B4 a NaN median with valid values = 0 is zero, listed; a NaN with valid values > 0 is a stop.
  B5 the SD column is carried per row; w_low / w_high = (v -/+ SD) * MW over the unchanged
     denominator, clipped at zero; a shared row adds its whole value to each member's w_high.
  Nothing is named in code: anchor genes for the gel bands come from the Carroll config (D52).

Usage
  python pipeline/mass_fractions.py            # build
  python pipeline/mass_fractions.py --check    # rebuild in memory, compare with disk, exit 1 on difference
"""

from __future__ import annotations

import argparse
import hashlib
import io
import math
import sys
from collections import defaultdict
from pathlib import Path

from pipeline import common
from pipeline.literature_inventory import find_rar_tool, read_archive_members

DECISIONS_INI = common.CONFIG_DIR / "mass_fraction_decisions.ini"   # column maps, digest rules
SOURCES_INI = common.LITERATURE_SOURCES_INI                          # the [source.*] sections (Step 7a)
CARROLL_INI = common.CONFIG_DIR / "carroll_classical_fractionation.ini"
MANIFEST_INI = common.LITERATURE_MANIFEST_INI
MASS_FRACTIONS_TSV = common.MASS_FRACTIONS_TSV                 # generated config (D61)
OUT_DIR = common.MASS_FRACTIONS_DIR
DIGEST_DIR = common.DIGEST_DIR
ISOFORM_DELTAS_TSV = common.PROTEIN_SET_OUT_DIR / "isoform_deltas.tsv"
PROCESSING_DELTAS_TSV = common.COMPOSITION_DIR / "processing_mass_deltas.tsv"
PTM_DELTAS_TSV = common.COMPOSITION_DIR / "ptm_mass_deltas.tsv"
SHARED_PAIRS_TSV = DIGEST_DIR / "shared_pairs.tsv"

PRIMARY_FILE_KEY = "file.1"             # the primary source is the [source.*] section whose role = primary (D32)
FIBER_TYPES = ("I", "IIa", "IIx")
SIG = 10                                # significant digits for fractions (as in composition, D29)


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------

def read_tsv_skip_comments(path: Path) -> list[dict[str, str]]:
    """TSV whose first lines may be '# ' comments (digest and inventory outputs)."""
    if not path.exists():
        raise FileNotFoundError(f"required file missing: {path}")
    with open(path, encoding="utf-8") as fh:
        header = None
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            rows.append(dict(zip(header, line.split("\t"))))
    return rows


def fnum(x: float) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "nan"
    return f"{x:.{SIG}g}"


def sha256_path(p: Path) -> str:
    return common.sha256_file(p)


def manifest_file(manifest, source_id: str, file_key: str):
    sec = f"file.{source_id}.{file_key}"
    if not manifest.has_section(sec):
        raise SystemExit(f"[STOP] manifest has no section {sec}; run fetch_literature first")
    return manifest[sec]


def verify_file(manifest, source_id: str, file_key: str, log) -> tuple[Path, str]:
    """B0: re-hash a literature file against the manifest. Returns (path, sha256)."""
    m = manifest_file(manifest, source_id, file_key)
    path = common.REPO_ROOT / m["path"]
    if not path.exists():
        raise SystemExit(f"[STOP] B0: file missing on disk: {m['path']}")
    got = sha256_path(path)
    if got != m["sha256"].lower():
        raise SystemExit(f"[STOP] B0: {m['path']} sha256 {got} != manifest {m['sha256']}")
    log.info("B0 %s.%s %s sha256 %s ok", source_id, file_key, m["path"], got)
    return path, got


def parse_value(cell, missing_token: str) -> float:
    if cell is None:
        return float("nan")
    if isinstance(cell, (int, float)):
        return float(cell)
    s = str(cell).strip()
    if s == "" or s == missing_token:
        return float("nan")
    return float(s)


# ----------------------------------------------------------------------------
# inputs
# ----------------------------------------------------------------------------

def load_pool() -> dict[str, dict]:
    cp = common.read_ini(common.ACCESSIONS_INI)
    pool = {}
    for acc in cp.sections():
        s = cp[acc]
        pool[acc] = {"accession": acc, "gene": s["gene"].strip(), "tier": s["tier"],
                     "tiers": [t for t in s["tier"].split(";") if t]}
    return pool


def load_master_mw() -> dict[str, float]:
    mw = {}
    for r in common.read_tsv(common.COMPOSITION_TSV):
        if r["segment_set"] == "master":
            mw[r["accession"]] = float(r["mw"])
    return mw


def read_dataset1(path: Path, cols, log) -> tuple[list[dict], dict]:
    """Every data row of Dataset 1 as read: row index (spreadsheet numbering), gene cell,
    median / sd / valid values per fiber type. Header matched by exact string after
    stripping whitespace (config `header_whitespace`); the fact that stripping was
    needed is recorded."""
    from openpyxl import load_workbook
    wb = load_workbook(io.BytesIO(path.read_bytes()), read_only=True, data_only=True)
    if cols["sheet"] not in wb.sheetnames:
        raise SystemExit(f"[STOP] sheet {cols['sheet']!r} not in {path.name}: {wb.sheetnames}")
    ws = wb[cols["sheet"]]
    header_row = int(cols["header_row"])
    wanted = {"identity": cols["identity_column"]}
    for ft in FIBER_TYPES:
        wanted[f"median_{ft}"] = cols[f"median_type_{ft}"]
        wanted[f"sd_{ft}"] = cols[f"sd_type_{ft}"]
        wanted[f"valid_{ft}"] = cols[f"valid_values_type_{ft}"]
    colmap: dict[str, int] = {}
    stripped_needed = []
    rows_out = []
    n_nonempty = 0
    for r_index, row in enumerate(ws.iter_rows(values_only=True), 1):
        if r_index < header_row:
            continue
        if r_index == header_row:
            headers = ["" if v is None else str(v) for v in row]
            for key, want in wanted.items():
                hits = [j for j, h in enumerate(headers) if h.strip() == want.strip()]
                if len(hits) != 1:
                    raise SystemExit(f"[STOP] header {want!r} matched {len(hits)} columns in row {header_row}")
                colmap[key] = hits[0]
                if headers[hits[0]] != headers[hits[0]].strip() or want != want.strip():
                    stripped_needed.append(want)
            continue
        if all(v is None or str(v).strip() == "" for v in row):
            continue
        n_nonempty += 1
        gene_cell = row[colmap["identity"]]
        gene = "" if gene_cell is None else str(gene_cell).strip()
        rec = {"row_index": r_index, "gene_cell": gene}
        for ft in FIBER_TYPES:
            rec[f"median_{ft}"] = parse_value(row[colmap[f"median_{ft}"]], cols["missing_value_token"])
            rec[f"sd_{ft}"] = parse_value(row[colmap[f"sd_{ft}"]], cols["missing_value_token"])
            vv = row[colmap[f"valid_{ft}"]]
            rec[f"valid_{ft}"] = int(float(vv)) if vv is not None and str(vv).strip() != "" else 0
        rows_out.append(rec)
    wb.close()
    meta = {"header_row": header_row, "sheet": cols["sheet"], "columns": {k: v for k, v in wanted.items()},
            "column_index": colmap, "header_whitespace_stripped_for": stripped_needed,
            "n_nonempty_rows_after_header": n_nonempty}
    log.info("Dataset 1: %d non-empty data rows; header whitespace stripped for %d columns",
             n_nonempty, len(stripped_needed))
    return rows_out, meta


# ----------------------------------------------------------------------------
# the join (B1, B4, D47, D48)
# ----------------------------------------------------------------------------

def join(rows: list[dict], pool: dict[str, dict], log) -> dict:
    """Returns per-accession measured values and the bookkeeping tables."""
    gene_to_acc = {p["gene"]: acc for acc, p in pool.items()}
    via_mapping: dict[str, str] = {}
    for mp in sorted(common.DATA_DIR.glob("tier*_gene_mapping.tsv")):     # written by enumerate_pool (D56)
        for r in common.read_tsv(mp):
            if r["accessions"] and ";" not in r["accessions"] and r["accessions"] in pool and r["gene"] not in gene_to_acc:
                via_mapping[r["gene"]] = r["accessions"]              # a dataset symbol that UniProt resolves to a pool entry (B1b)
    gene_to_acc.update(via_mapping)
    log.info("B1: %d pool symbols; B1b: %d dataset symbols joined through the gene-mapping table", len(pool), len(via_mapping))
    stops: list[str] = []

    # B4: a NaN median with valid values > 0 cannot be settled
    for r in rows:
        for ft in FIBER_TYPES:
            if math.isnan(r[f"median_{ft}"]) and r[f"valid_{ft}"] > 0:
                stops.append(f"row {r['row_index']} {r['gene_cell']!r}: median {ft} is NaN but valid values = {r[f'valid_{ft}']} (B4)")

    single_rows: dict[str, list[dict]] = defaultdict(list)   # gene -> its single-name rows
    multi_rows: list[dict] = []
    blank_rows: list[dict] = []
    for r in rows:
        g = r["gene_cell"]
        if g == "":
            blank_rows.append(r)
        elif ";" in g:
            multi_rows.append(r)
        else:
            single_rows[g].append(r)

    def v(r, ft):  # B4: NaN with valid = 0 -> 0
        x = r[f"median_{ft}"]
        return 0.0 if math.isnan(x) else x

    def sd(r, ft):
        x = r[f"sd_{ft}"]
        return 0.0 if math.isnan(x) else x

    # D47: sum single-name rows per gene
    per_gene: dict[str, dict] = {}
    for g, rs in single_rows.items():
        rec = {"gene": g, "n_rows": len(rs), "row_indices": [r["row_index"] for r in rs]}
        for ft in FIBER_TYPES:
            vals = [v(r, ft) for r in rs]
            total = sum(vals)
            rec[f"v_{ft}"] = total
            rec[f"sd_{ft}"] = math.sqrt(sum(sd(r, ft) ** 2 for r in rs))   # independent-sum convention, disclosed
            rec[f"valid_{ft}"] = max(r[f"valid_{ft}"] for r in rs)
            rec[f"minor_share_{ft}"] = (total - max(vals)) / total if total > 0 else 0.0
        per_gene[g] = rec

    # measured values per pool accession
    measured: dict[str, dict] = {}
    for acc, p in pool.items():
        g = p["gene"]
        rec = {"accession": acc, "gene": g, "tier": p["tier"], "tiers": p["tiers"],
               "match_rule": "none", "row_indices": [], "n_rows": 0}
        for ft in FIBER_TYPES:
            rec[f"v_{ft}"] = 0.0
            rec[f"sd_{ft}"] = 0.0
            rec[f"valid_{ft}"] = 0
            rec[f"minor_share_{ft}"] = 0.0
            rec[f"shared_add_{ft}"] = 0.0     # value assigned from ';' rows (D48)
            rec[f"shared_whole_{ft}"] = 0.0   # whole ';' row values, for w_high (D48)
        mapped_symbols = [s for s, a in via_mapping.items() if a == acc]
        for s in mapped_symbols:
            if s in per_gene:
                g = s
        if g in per_gene:
            pg = per_gene[g]
            rec["match_rule"] = ("single" if pg["n_rows"] == 1 else f"summed:{pg['n_rows']}") + ("" if g == p["gene"] else f"+via_mapping:{g}")
            rec["row_indices"] = list(pg["row_indices"])
            rec["n_rows"] = pg["n_rows"]
            for ft in FIBER_TYPES:
                for k in ("v", "sd", "valid", "minor_share"):
                    rec[f"{k}_{ft}"] = pg[f"{k}_{ft}"]
        measured[acc] = rec

    # D48: multi-gene cells
    shared_table = []
    for r in multi_rows:
        members = [m.strip() for m in r["gene_cell"].split(";") if m.strip()]
        in_pool = [m for m in members if m in gene_to_acc]
        rec = {"row_index": r["row_index"], "gene_cell": r["gene_cell"], "n_members": len(members),
               "members_in_pool": ";".join(in_pool)}
        for ft in FIBER_TYPES:
            rec[f"v_{ft}"] = v(r, ft)
        if not in_pool:
            rec["outcome"] = "ignored:no_member_in_pool"
        elif len(in_pool) == 1:
            acc = gene_to_acc[in_pool[0]]
            m = measured[acc]
            m["match_rule"] = "shared_group" if m["match_rule"] == "none" else m["match_rule"] + "+shared_group"
            m["row_indices"].append(r["row_index"])
            for ft in FIBER_TYPES:
                m[f"shared_add_{ft}"] += v(r, ft)
                m[f"shared_whole_{ft}"] += v(r, ft)
            rec["outcome"] = f"shared_group:{in_pool[0]}"
        else:
            # divide in proportion to the members' own single-name values (equally if none)
            parts = []
            for ft in FIBER_TYPES:
                own = [per_gene[m][f"v_{ft}"] if m in per_gene else 0.0 for m in in_pool]
                tot = sum(own)
                shares = [o / tot for o in own] if tot > 0 else [1.0 / len(in_pool)] * len(in_pool)
                for m, s in zip(in_pool, shares):
                    acc = gene_to_acc[m]
                    measured[acc][f"shared_add_{ft}"] += v(r, ft) * s
                    measured[acc][f"shared_whole_{ft}"] += v(r, ft)
                parts.append(f"{ft}:" + ",".join(f"{m}={s:.3f}" for m, s in zip(in_pool, shares)))
            for m in in_pool:
                acc = gene_to_acc[m]
                mm = measured[acc]
                mm["match_rule"] = "shared_split" if mm["match_rule"] == "none" else mm["match_rule"] + "+shared_split"
                mm["row_indices"].append(r["row_index"])
            rec["outcome"] = "shared_split " + " ".join(parts)
        shared_table.append(rec)

    for m in measured.values():
        for ft in FIBER_TYPES:
            m[f"v_{ft}"] += m[f"shared_add_{ft}"]

    # completeness: single-name genes not in the pool, plus ';' rows with no pool member
    outside = []
    for g, pg in per_gene.items():
        if g not in gene_to_acc:
            outside.append({"gene": g, "kind": "single_name", "n_rows": pg["n_rows"],
                            **{f"v_{ft}": pg[f"v_{ft}"] for ft in FIBER_TYPES},
                            **{f"valid_{ft}": pg[f"valid_{ft}"] for ft in FIBER_TYPES}})
    for r in multi_rows:
        members = [m.strip() for m in r["gene_cell"].split(";") if m.strip()]
        if not any(m in gene_to_acc for m in members):
            outside.append({"gene": r["gene_cell"], "kind": "multi_gene_cell", "n_rows": 1,
                            **{f"v_{ft}": v(r, ft) for ft in FIBER_TYPES},
                            **{f"valid_{ft}": r[f"valid_{ft}"] for ft in FIBER_TYPES}})

    # signal share arriving through ';' rows, per tier and fiber type (D48 disclosure)
    counts = {"rows_total": len(rows), "rows_blank_gene": len(blank_rows), "rows_single_gene": sum(len(x) for x in single_rows.values()),
              "genes_single_name": len(per_gene), "genes_on_more_than_one_row": sum(1 for pg in per_gene.values() if pg["n_rows"] > 1),
              "rows_multi_gene_cell": len(multi_rows), "rows_nan_median_valid0": sum(
                  1 for r in rows for ft in FIBER_TYPES if math.isnan(r[f"median_{ft}"]) and r[f"valid_{ft}"] == 0)}
    log.info("join: %s", counts)
    return {"measured": measured, "shared_table": shared_table, "outside": outside, "counts": counts, "stops": stops}


# ----------------------------------------------------------------------------
# weights (B2, B5)
# ----------------------------------------------------------------------------

def weights(measured: dict[str, dict], mw: dict[str, float], log) -> dict:
    tiers = sorted({t for m in measured.values() for t in m["tiers"]})
    denom = {(t, ft): 0.0 for t in tiers for ft in FIBER_TYPES}
    for acc, m in measured.items():
        if acc not in mw:
            raise SystemExit(f"[STOP] no master-chain MW for {acc} in composition table (B2)")
        m["mw"] = mw[acc]
        for t in m["tiers"]:
            for ft in FIBER_TYPES:
                denom[(t, ft)] += m[f"v_{ft}"] * m["mw"]
    for acc, m in measured.items():
        for t in tiers:
            for ft in FIBER_TYPES:
                key = f"w_{ft}_tier{t}"
                if t not in m["tiers"]:
                    m[key] = None
                    m[f"{key}_low"] = None
                    m[f"{key}_high"] = None
                    continue
                d = denom[(t, ft)]
                w = m[f"v_{ft}"] * m["mw"] / d if d > 0 else 0.0
                lo = max(0.0, (m[f"v_{ft}"] - m[f"sd_{ft}"])) * m["mw"] / d if d > 0 else 0.0
                hi = (m[f"v_{ft}"] + m[f"sd_{ft}"] + (m[f"shared_whole_{ft}"] - m[f"shared_add_{ft}"])) * m["mw"] / d if d > 0 else 0.0
                m[key], m[f"{key}_low"], m[f"{key}_high"] = w, lo, hi
    # D58: the combined standard — every entry of every tier, one denominator
    denom_all = {ft: sum(m[f"v_{ft}"] * m["mw"] for m in measured.values()) for ft in FIBER_TYPES}
    for m in measured.values():
        for ft in FIBER_TYPES:
            d = denom_all[ft]
            m[f"w_{ft}_combined"] = m[f"v_{ft}"] * m["mw"] / d if d > 0 else 0.0
            m[f"w_{ft}_combined_low"] = max(0.0, m[f"v_{ft}"] - m[f"sd_{ft}"]) * m["mw"] / d if d > 0 else 0.0
            m[f"w_{ft}_combined_high"] = (m[f"v_{ft}"] + m[f"sd_{ft}"] + (m[f"shared_whole_{ft}"] - m[f"shared_add_{ft}"])) * m["mw"] / d if d > 0 else 0.0
    shared_share = {}
    for t in tiers:
        for ft in FIBER_TYPES:
            d = denom[(t, ft)]
            s = sum(m[f"shared_add_{ft}"] * m["mw"] for m in measured.values() if t in m["tiers"])
            shared_share[(t, ft)] = s / d if d > 0 else 0.0
    log.info("denominators (sum v*MW per tier and fiber type): %s", {f"{t}/{ft}": f"{v:.6g}" for (t, ft), v in denom.items()})
    return {"tiers": tiers, "denom": denom, "denom_all": denom_all, "shared_share": shared_share}


# ----------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------

def load_subtree(tier: str) -> tuple[dict[str, str], dict[str, set[str]]]:
    """term -> name; term -> set of descendants (incl. itself)."""
    rows = common.read_tsv(common.ONTOLOGY_DIR / f"tier{tier}_subtree.tsv")
    name = {r["term_id"]: r["name"] for r in rows}
    parent = {r["term_id"]: r["parent_id"] for r in rows}
    desc = {t: {t} for t in name}
    for t in name:
        cur = t
        while parent.get(cur):
            desc[parent[cur]].add(t)
            cur = parent[cur]
    return name, desc


def branch_table(measured, tiers) -> list[list]:
    out = []
    for t in tiers:
        if not (common.ONTOLOGY_DIR / f"tier{t}_subtree.tsv").exists():
            continue                                   # measured tier (D56): no ontology branch
        name, desc = load_subtree(t)
        pool_rows = common.read_tsv(common.DATA_DIR / f"tier{t}_pool.tsv")
        terms_of = {r["accession"]: set(r["subtree_terms_returning_entry"].split(";")) for r in pool_rows}
        for term in name:
            excl = [a for a, ts in terms_of.items() if ts and ts <= desc[term] and a in measured]
            anyv = [a for a, ts in terms_of.items() if term in ts and a in measured]
            row = [t, term, name[term], len(anyv), len(excl)]
            for ft in FIBER_TYPES:
                row.append(fnum(sum(measured[a][f"w_{ft}_tier{t}"] or 0.0 for a in anyv)))
            for ft in FIBER_TYPES:
                row.append(fnum(sum(measured[a][f"w_{ft}_tier{t}"] or 0.0 for a in excl)))
            row.append(";".join(sorted(excl)))
            out.append(row)
    return out


def weighted_bounds(measured, tiers) -> tuple[list[list], dict]:
    aas = list(common.AMINO_ACIDS)
    iso: dict[str, dict[str, float]] = defaultdict(lambda: {a: 0.0 for a in aas})
    if ISOFORM_DELTAS_TSV.exists():
        for r in common.read_tsv(ISOFORM_DELTAS_TSV):
            for a in aas:
                iso[r["accession"]][a] = max(iso[r["accession"]][a], abs(float(r[f"delta_{a}"])))
    proc: dict[str, dict[str, float]] = {}
    for r in common.read_tsv(PROCESSING_DELTAS_TSV):
        proc[r["accession"]] = {a: abs(float(r[f"delta_free_frac_{a}"])) for a in aas}
    ptm: dict[str, dict] = {r["accession"]: r for r in common.read_tsv(PTM_DELTAS_TSV)}
    rows, extra = [], {}
    for t in tiers:
        for ft in FIBER_TYPES:
            wk = f"w_{ft}_tier{t}"
            for kind, table in (("isoform", iso), ("processing", proc)):
                for a in aas:
                    best, where = 0.0, ""
                    for acc, m in measured.items():
                        w = m.get(wk)
                        if not w:
                            continue
                        d = table.get(acc, {}).get(a, 0.0) * w
                        if d > best:
                            best, where = d, acc
                    rows.append([t, ft, kind, a, fnum(best), where])
            best, where = 0.0, ""
            glyc = 0.0
            for acc, m in measured.items():
                w = m.get(wk)
                if not w:
                    continue
                p = ptm.get(acc)
                if p:
                    d = float(p["ptm_mass_delta_fraction_of_mw"]) * w
                    if d > best:
                        best, where = d, acc
                    if int(p["n_glycosylation"]) > 0:
                        glyc += w
            rows.append([t, ft, "ptm_mass", "all", fnum(best), where])
            rows.append([t, ft, "glycosylated_entries_sum_w", "all", fnum(glyc), ""])
            extra[(t, ft)] = glyc
    return rows, extra


def family_bounds(measured, tiers, cutoff: int) -> list[list]:
    """D55: for every shared-peptide family at the cutoff with >= 2 pool members, the spread of
    each amino acid's free-mass fraction across the members (max - min) times the family's
    weight — what the within-family split could cost the profile if it were entirely wrong."""
    aas = list(common.AMINO_ACIDS)
    frac = {}
    for r in common.read_tsv(common.COMPOSITION_TSV):
        if r["segment_set"] == "master":
            frac[r["accession"]] = {a: float(r[f"free_frac_{a}"]) for a in aas}
    fams = families_from_pairs(cutoff)
    groups: dict[str, set[str]] = {}
    for acc, members in fams.items():
        if acc in measured:
            groups[min(members)] = {m for m in members if m in measured}
    rows = []
    for fid, members in sorted(groups.items()):
        if len(members) < 2:
            continue
        genes = ";".join(sorted(measured[m]["gene"] for m in members))
        spread = {a: max(frac[m][a] for m in members) - min(frac[m][a] for m in members) for a in aas}
        for tier in tiers:
            for ft in FIBER_TYPES:
                fw = sum(measured[m].get(f"w_{ft}_tier{tier}") or 0.0 for m in members)
                if fw == 0.0:
                    continue
                bounds = {a: spread[a] * fw for a in aas}
                a_max = max(aas, key=lambda a: bounds[a])
                rows.append([cutoff, fid, len(members), genes, tier, ft, fnum(fw), fnum(bounds[a_max]), a_max]
                            + [fnum(bounds[a]) for a in aas])
    return rows


def slow_fast_ratio_check(dec, manifest, measured, log) -> tuple[list[list], dict]:
    """The ratio cross-check (D34): the [columns.*] map whose `use` names a slow/fast ratio and
    whose scale is log2 — per pool gene, 2^(mean over participants of slow - fast) against the
    primary dataset's I/IIa and I/IIx value ratios. Joined by gene name (B1)."""
    maps = [s for s in dec.sections() if s.startswith("columns.") and dec[s].get("slow_columns") and dec[s].get("identity_kind") == "gene_name"]
    if len(maps) != 1:
        raise SystemExit(f"[STOP] expected exactly one [columns.*] map with slow_columns and identity_kind = gene_name, found {maps}")
    sec = dec[maps[0]]
    _, source_id, file_key = maps[0].split(".", 2)
    path, sha = verify_file(manifest, source_id, file_key, log)
    from openpyxl import load_workbook
    wb = load_workbook(io.BytesIO(path.read_bytes()), read_only=True, data_only=True)
    ws = wb[sec["sheet"]]
    it = ws.iter_rows(values_only=True)
    for _ in range(int(sec["header_row"]) - 1):
        next(it)
    headers = [("" if v is None else str(v).strip()) for v in next(it)]
    def idx(name):
        hits = [j for j, h in enumerate(headers) if h == name.strip()]
        if len(hits) != 1:
            raise SystemExit(f"[STOP] header {name!r} matched {len(hits)} columns in {path.name}")
        return hits[0]
    j_gene = idx(sec["identity_column"])
    slow = [idx(c) for c in sec["slow_columns"].split(";") if c.strip()]
    fast = [idx(c) for c in sec["fast_columns"].split(";") if c.strip()]
    if len(slow) != len(fast):
        raise SystemExit("[STOP] slow_columns and fast_columns differ in length; participants must pair")
    by_gene: dict[str, list[float]] = {}
    n_rows = 0
    for row in it:
        g = row[j_gene]
        if g is None or str(g).strip() == "":
            continue
        n_rows += 1
        diffs = []
        for js, jf in zip(slow, fast):
            a, b = parse_value(row[js], "NA"), parse_value(row[jf], "NA")
            if not (math.isnan(a) or math.isnan(b)):
                diffs.append(a - b if sec["scale"] == "log2" else math.log2(a / b) if a > 0 and b > 0 else float("nan"))
        by_gene.setdefault(str(g).strip(), []).extend(d for d in diffs if not math.isnan(d))
    wb.close()
    gene_to_acc = {m["gene"]: acc for acc, m in measured.items()}
    rows, log2_disagreement = [], []
    for g, diffs in sorted(by_gene.items()):
        acc = gene_to_acc.get(g)
        if acc is None or not diffs:
            continue
        m = measured[acc]
        mean = sum(diffs) / len(diffs)
        sd = math.sqrt(sum((d - mean) ** 2 for d in diffs) / (len(diffs) - 1)) if len(diffs) > 1 else float("nan")
        r_sf = 2 ** mean
        r_ia = m["v_I"] / m["v_IIa"] if m["v_IIa"] > 0 else float("nan")
        r_ix = m["v_I"] / m["v_IIx"] if m["v_IIx"] > 0 else float("nan")
        q = r_ia / r_sf if r_sf > 0 and not math.isnan(r_ia) else float("nan")
        if not math.isnan(q) and q > 0:
            log2_disagreement.append(abs(math.log2(q)))
        rows.append([acc, g, m["tier"], fnum(m["v_I"]), fnum(m["v_IIa"]), fnum(m["v_IIx"]), fnum(r_ia), fnum(r_ix),
                     len(diffs), fnum(mean), fnum(sd), fnum(r_sf), fnum(q), fnum(r_ix / r_sf) if r_sf > 0 and not math.isnan(r_ix) else "nan"])
    log2_disagreement.sort()
    summ = {"source": source_id, "file": f"{path.name} sha256 {sha}", "dataset_rows": n_rows, "pool_entries_compared": len(rows),
            "median_abs_log2_ratio_disagreement_I_over_IIa": fnum(log2_disagreement[len(log2_disagreement) // 2]) if log2_disagreement else "nan"}
    log.info("slow/fast ratio check %s: %d rows, %d pool entries compared", source_id, n_rows, len(rows))
    return rows, summ


def families_from_pairs(min_shared: int, pairs_path: Path = SHARED_PAIRS_TSV) -> dict[str, set[str]]:
    """Union-find over shared_pairs.tsv with n_shared_peptides >= min_shared."""
    parent: dict[str, str] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r in read_tsv_skip_comments(pairs_path):
        if int(r["n_shared_peptides"]) >= min_shared:
            a, b = find(r["accession_a"]), find(r["accession_b"])
            if a != b:
                parent[a] = b
    fam: dict[str, set[str]] = defaultdict(set)
    for x in list(parent):
        fam[find(x)].add(x)
    return {m: members for members in fam.values() for m in members}


def carroll_check(measured, pool, log) -> tuple[list[list], list[list]]:
    cp = common.read_ini(CARROLL_INI)
    gene_to_acc = {p["gene"]: acc for acc, p in pool.items()}
    bands = {}
    for sec in cp.sections():
        if sec.startswith("band."):
            bands[sec[len("band."):]] = [g.strip() for g in cp[sec]["anchor_genes"].split(";") if g.strip()]
    fam_rows, check_rows = [], []
    band_w: dict[tuple[int, str, str], float] = {}
    for cutoff in (1, 2):
        fams = families_from_pairs(cutoff)
        for band, anchors in bands.items():
            accs = set()
            for g in anchors:
                if g not in gene_to_acc:
                    raise SystemExit(f"[STOP] anchor gene {g} of band {band} is not in the pool (D52)")
                a = gene_to_acc[g]
                accs |= fams.get(a, {a})
            for a in sorted(accs, key=lambda a: (-(measured[a].get("w_I_tier1") or 0.0), a)):   # ties by accession: reproducible order
                m = measured[a]
                fam_rows.append([cutoff, band, a, m["gene"], m["tier"], "anchor" if m["gene"] in anchors else "family",
                                 m["match_rule"]] + [fnum(m[f"v_{ft}"]) for ft in FIBER_TYPES]
                                + [fnum(m.get(f"w_{ft}_tier1") or 0.0) for ft in FIBER_TYPES])
            for ft in ("I", "IIa"):
                band_w[(cutoff, band, ft)] = sum(measured[a].get(f"w_{ft}_tier1") or 0.0 for a in accs)
    # gel values
    def gel(band, ft):
        s = cp[f"{band}.type_{ft}"]
        return float(s["value"]), float(s["spread"]), s["unit"]
    for cutoff in (1, 2):
        for ft in ("I", "IIa"):
            m_val, m_sd, unit = gel("myosin_heavy_chain", ft)
            a_val, a_sd, _ = gel("actin", ft)
            wm = band_w[(cutoff, "myosin_heavy_chain", ft)]
            wa = band_w[(cutoff, "actin", ft)]
            gel_ratio = m_val / a_val
            our_ratio = wm / wa if wa > 0 else float("nan")
            check_rows.append([cutoff, ft, "ibaq_x_mw_within_tier1", fnum(wm), fnum(wa), fnum(our_ratio),
                               f"{m_val}±{m_sd}", f"{a_val}±{a_sd}", unit, fnum(gel_ratio),
                               fnum(our_ratio / gel_ratio) if gel_ratio else "nan"])
            log.info("gel comparison cutoff %d type %s: MHC:actin iBAQxMW %.4f gel %.4f", cutoff, ft, our_ratio, gel_ratio)
    families = {cutoff: {band: set() for band in bands} for cutoff in (1, 2)}
    for cutoff in (1, 2):
        fams = families_from_pairs(cutoff)
        for band, anchors in bands.items():
            for g in anchors:
                a = gene_to_acc[g]
                families[cutoff][band] |= fams.get(a, {a})
    return fam_rows, check_rows, families, cp


def pool_accessions_in_cell(cell: str, pool) -> list[str]:
    """Pool accessions named in a ';'-separated accession cell; an isoform suffix (-N) is
    the same entry. Sorted, unique."""
    accs = {a.strip().split("-")[0] for a in cell.split(";") if a.strip()}
    return sorted(a for a in accs if a in pool)


def intensity_share_check(dec, manifest, families: dict, carroll_cp, pool, measured_ref, log) -> tuple[list[list], list[list]]:
    """MHC:actin from the intensity share in the mass-spec source whose [columns.<source>.<file>]
    map declares identity_kind = accession_list (Deshmukh 2021, Supplementary Data 3). Rows are
    joined to pool entries by accession: any accession in the row's identity cell, isoform suffix
    dropped, equal to a pool accession. Several rows for one entry are summed. Band = sum of
    the intensity over the D52 family. Nothing is named; the sheets, columns, and file come from
    the config map."""
    maps = [s for s in dec.sections() if s.startswith("columns.") and dec[s].get("identity_kind") == "accession_list"]
    if len(maps) != 1:
        raise SystemExit(f"[STOP] expected exactly one [columns.*] map with identity_kind = accession_list, found {maps}")
    sec = dec[maps[0]]
    _, source_id, file_key = maps[0].split(".", 2)
    path, sha = verify_file(manifest, source_id, file_key, log)
    rar_tool = find_rar_tool(None, log)
    members = {name: data for name, _, data in read_archive_members(path, rar_tool)} if path.suffix.lower() in (".zip", ".rar") else {path.name: path.read_bytes()}
    if sec["member"] not in members:
        raise SystemExit(f"[STOP] member {sec['member']!r} not in {path.name}: {sorted(members)}")
    from openpyxl import load_workbook
    wb = load_workbook(io.BytesIO(members[sec["member"]]), read_only=True, data_only=True)
    header_row = int(sec["header_row"])
    out = []
    per_acc_by_sheet: dict[str, dict[str, float]] = {}
    for sheet_key, ft_label, gel_ft in (("sheet_slow", "slow", "I"), ("sheet_fast", "fast", "IIa")):
        ws = wb[sec[sheet_key]]
        it = ws.iter_rows(values_only=True)
        for _ in range(header_row - 1):
            next(it)
        headers = [("" if v is None else str(v).strip()) for v in next(it)]
        col = {}
        for key in ("identity_column", "intensity_column", "gene_column"):
            hits = [j for j, h in enumerate(headers) if h == sec[key].strip()]
            if len(hits) != 1:
                raise SystemExit(f"[STOP] header {sec[key]!r} matched {len(hits)} columns in sheet {sec[sheet_key]!r}")
            col[key] = hits[0]
        per_acc: dict[str, float] = defaultdict(float)
        per_acc_by_sheet[ft_label] = per_acc
        n_rows = n_matched = 0
        for row in it:
            ids = row[col["identity_column"]]
            if ids is None:
                continue
            n_rows += 1
            val = parse_value(row[col["intensity_column"]], "NaN")
            if math.isnan(val):
                continue
            hit = pool_accessions_in_cell(str(ids), pool)
            if hit:
                n_matched += 1
                for a in hit:
                    per_acc[a] += val / len(hit)
        log.info("%s %s sheet %r: %d rows, %d matched a pool entry by accession", source_id, file_key, sec[sheet_key], n_rows, n_matched)
        for cutoff in (1, 2):
            mhc = sum(per_acc.get(a, 0.0) for a in families[cutoff]["myosin_heavy_chain"])
            act = sum(per_acc.get(a, 0.0) for a in families[cutoff]["actin"])
            ratio = mhc / act if act > 0 else float("nan")
            g = carroll_cp[f"myosin_heavy_chain.type_{gel_ft}"]; a_ = carroll_cp[f"actin.type_{gel_ft}"]
            gel_ratio = float(g["value"]) / float(a_["value"])
            out.append([cutoff, f"{ft_label} (vs gel type {gel_ft})", f"intensity_share_{source_id}", fnum(mhc), fnum(act), fnum(ratio),
                        f"{g['value']}±{g['spread']}", f"{a_['value']}±{a_['spread']}", g["unit"], fnum(gel_ratio),
                        fnum(ratio / gel_ratio) if gel_ratio else "nan"])
    wb.close()
    ratio_rows = []
    for acc in sorted(set(per_acc_by_sheet["slow"]) | set(per_acc_by_sheet["fast"])):
        s, f = per_acc_by_sheet["slow"].get(acc, 0.0), per_acc_by_sheet["fast"].get(acc, 0.0)
        m = measured_ref[acc]
        r_sf = s / f if f > 0 else float("nan")
        r_ia = m["v_I"] / m["v_IIa"] if m["v_IIa"] > 0 else float("nan")
        ratio_rows.append([acc, m["gene"], m["tier"], fnum(s), fnum(f), fnum(r_sf), fnum(r_ia),
                           fnum(r_ia / r_sf) if r_sf > 0 and not math.isnan(r_ia) else "nan"])
    return out, ratio_rows


def ibaq_vs_lfq_table(manifest, source_id: str, log) -> list[list]:
    """Per gene of the two tables: mean fraction over fibers, iBAQ (file.1) vs MaxLFQ (file.2), on the
    fibers present in both tables and on all fibers of each. Nothing is named."""
    from openpyxl import load_workbook
    tables = {}
    for key, label in (("file.1", "ibaq"), ("file.2", "lfq")):
        path, _ = verify_file(manifest, source_id, key, log)
        wb = load_workbook(io.BytesIO(path.read_bytes()), read_only=True, data_only=True)
        ws = wb.worksheets[0]
        it = ws.iter_rows(values_only=True)
        header = [("" if v is None else str(v).strip()) for v in next(it)]
        cols = header[1:]
        data = {}
        for row in it:
            if row[0] is None:
                continue
            data[str(row[0]).strip()] = {f: (float(v) if v is not None and str(v).strip() != "" else float("nan"))
                                         for f, v in zip(cols, row[1:])}
        wb.close()
        # a fraction column holds values in [0, 1]; a column that does not is a summary or
        # another quantity, not a fiber, and is excluded and listed (rule in the file header)
        fibers, excluded = [], []
        for c in cols:
            vals = [d[c] for d in data.values() if c in d and not math.isnan(d[c])]
            if c == "" or not vals or max(vals) > 1.0000001 or min(vals) < 0:
                excluded.append(c or "<blank header>")
            else:
                fibers.append(c)
        log.info("%s %s: %d columns, %d kept as fiber fraction columns, excluded: %s", source_id, key, len(cols), len(fibers), excluded)
        tables[label] = (fibers, data, excluded)
    fib_common = [f for f in tables["ibaq"][0] if f in set(tables["lfq"][0])]
    genes = sorted(set(tables["ibaq"][1]) | set(tables["lfq"][1]))
    rows = []

    def mean(vals):
        vals = [v for v in vals if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else float("nan")

    for g in genes:
        ib = tables["ibaq"][1].get(g, {})
        lf = tables["lfq"][1].get(g, {})
        mi_c = mean([ib.get(f, float("nan")) for f in fib_common])
        ml_c = mean([lf.get(f, float("nan")) for f in fib_common])
        rows.append([g, len(fib_common), fnum(mi_c), fnum(ml_c),
                     fnum(ml_c / mi_c) if mi_c and mi_c > 0 else "nan",
                     len(tables["ibaq"][0]), fnum(mean([ib.get(f, float("nan")) for f in tables["ibaq"][0]])),
                     len(tables["lfq"][0]), fnum(mean([lf.get(f, float("nan")) for f in tables["lfq"][0]]))])
    log.info("iBAQ-vs-LFQ tables: %d genes; %d fibers in both tables", len(genes), len(fib_common))
    note = (f"fiber columns: values in [0,1] with a non-blank header; excluded from file.1: {tables['ibaq'][2]}; "
            f"excluded from file.2: {tables['lfq'][2]}; fibers in both: {len(fib_common)}")
    return rows, note


# ----------------------------------------------------------------------------
# generated config
# ----------------------------------------------------------------------------

def render_fiber_config(ft: str, measured, tiers, meta: dict, header_common: list[str]) -> str:
    cp = common.new_ini()
    for acc in sorted(measured):
        m = measured[acc]
        cp.add_section(acc)
        s = cp[acc]
        s["gene"] = m["gene"]
        s["tier"] = m["tier"]
        for t in tiers:
            w = m.get(f"w_{ft}_tier{t}")
            if w is None:
                continue
            s[f"w_tier{t}"] = fnum(w)
            s[f"w_tier{t}_low"] = fnum(m[f"w_{ft}_tier{t}_low"])
            s[f"w_tier{t}_high"] = fnum(m[f"w_{ft}_tier{t}_high"])
        s["unit"] = "mass_fraction_of_tier"
        s["w_combined"] = fnum(m[f"w_{ft}_combined"])            # D58: one denominator over every tier
        s["w_combined_low"] = fnum(m[f"w_{ft}_combined_low"])
        s["w_combined_high"] = fnum(m[f"w_{ft}_combined_high"])
        s["unit_combined"] = "mass_fraction_of_combined_standard"
        s["v"] = fnum(m[f"v_{ft}"])
        s["sd"] = fnum(m[f"sd_{ft}"])
        s["valid_values"] = str(m[f"valid_{ft}"])
        s["mw_master"] = f"{m['mw']:.4f}"
        s["match_rule"] = m["match_rule"]
        s["rows"] = ";".join(str(i) for i in m["row_indices"]) if m["row_indices"] else "none"
        s["source"] = meta["source_line"]
        s["location"] = (f"sheet {meta['sheet']}, header row {meta['header_row']}, columns "
                         f"{meta['columns'][f'median_{ft}']!r} / {meta['columns'][f'sd_{ft}']!r} / {meta['columns'][f'valid_{ft}']!r}, rows as listed")
        s["retrieved"] = "see data/literature/manifest.ini (the file is identified by its hash; the fetch time is not repeated here)"
    header = [f"GENERATED by mass_fractions.py. DO NOT EDIT BY HAND.", f"fiber_type = {ft} (the dataset's own column, B3)"] + header_common
    return common.render_ini(cp, header)


def per_entry_table(measured, tiers, meta: dict, header_common: list[str]) -> str:
    """D61: the Layer B weights as one table — one row per pool accession, every fiber type
    and weight side by side, the constant provenance once in the header. An ADDED format for
    reading and for the aggregation stage; the per-fiber-type .ini files are written unchanged
    alongside it. This stage never deletes a file."""
    hdr = ["accession", "gene", "tier", "mw_master", "match_rule", "n_rows", "rows"]
    for ft in FIBER_TYPES:
        hdr += [f"v_{ft}", f"sd_{ft}", f"valid_{ft}", f"minor_row_share_{ft}", f"shared_row_value_{ft}"]
    for t in tiers:
        for ft in FIBER_TYPES:
            hdr += [f"w_{ft}_tier{t}", f"w_{ft}_tier{t}_low", f"w_{ft}_tier{t}_high"]
    for ft in FIBER_TYPES:
        hdr += [f"w_{ft}_combined", f"w_{ft}_combined_low", f"w_{ft}_combined_high"]
    body = []
    for acc in sorted(measured):
        m = measured[acc]
        row = [acc, m["gene"], m["tier"], f"{m['mw']:.4f}", m["match_rule"], m["n_rows"],
               ";".join(str(i) for i in m["row_indices"]) if m["row_indices"] else "none"]
        for ft in FIBER_TYPES:
            row += [fnum(m[f"v_{ft}"]), fnum(m[f"sd_{ft}"]), m[f"valid_{ft}"], fnum(m[f"minor_share_{ft}"]), fnum(m[f"shared_whole_{ft}"])]
        for t in tiers:
            for ft in FIBER_TYPES:
                w = m.get(f"w_{ft}_tier{t}")
                row += ["" if w is None else fnum(w), "" if w is None else fnum(m[f"w_{ft}_tier{t}_low"]),
                        "" if w is None else fnum(m[f"w_{ft}_tier{t}_high"])]
        for ft in FIBER_TYPES:
            row += [fnum(m[f"w_{ft}_combined"]), fnum(m[f"w_{ft}_combined_low"]), fnum(m[f"w_{ft}_combined_high"])]
        body.append(row)
    provenance = [
        "GENERATED CONFIG (D61): the Layer B weights as one table, one row per pool accession; the same values as config/mass_fractions/<type>.ini. Rules B0-B8.",
        f"source      = {meta['source_line']}",
        "retrieved   = see data/literature/manifest.ini (the file is identified by its hash above; the fetch time is not repeated here, so a re-fetch of an unchanged file does not change this table)",
        f"location    = sheet {meta['sheet']}, header row {meta['header_row']}; column `rows` = the dataset row numbers each entry was read from",
    ]
    for ft in FIBER_TYPES:
        provenance.append(f"columns_{ft} = median {meta['columns'][f'median_{ft}']!r} / sd {meta['columns'][f'sd_{ft}']!r} / valid values {meta['columns'][f'valid_{ft}']!r}  (the dataset's own fiber type, B3)")
    provenance += [
        "v = median value (summed over the entry's rows, D47; split shares of multi-gene rows, D48); sd = the published standard deviation (root-sum-square over summed rows, D53a);",
        "  valid = the dataset's valid-value count; a NaN median with valid = 0 is v = 0 (B4). mw_master = Step-2 mature-chain molecular weight.",
        "w_<type>_tier<N> = v x mw / sum over the tier (D31; unit: mass fraction of the tier). w_<type>_combined = v x mw / sum over every tier (D58; unit: mass fraction of the combined standard).",
        "_low / _high = (v -/+ sd) x mw over the unchanged denominator, clipped at zero (B5, D45); _high also carries the whole value of every multi-gene row the entry shares (D48).",
        "The sd is Dataset 1's spread BETWEEN FIBERS of one type (D45): the interval it gives is 'a random pure fiber of this type', not the uncertainty of the median.",
        "match_rule: single | summed:n | shared_group | shared_split | via_mapping:<symbol> | none (no dataset row; w = 0).",
    ]
    return tsv_text(provenance + header_common, hdr, body)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def build(log) -> dict:
    dec = common.read_ini(DECISIONS_INI)
    src_cp = common.read_ini(SOURCES_INI)
    manifest = common.read_ini(MANIFEST_INI)
    primaries = [s for s in src_cp.sections() if s.startswith("source.") and src_cp[s].get("role") == "primary"]
    if len(primaries) != 1:
        raise SystemExit(f"[STOP] expected exactly one [source.*] with role = primary, found {primaries}")
    primary = primaries[0][len("source."):]
    src = src_cp[f"source.{primary}"]
    cols = dec[f"columns.{primary}.{PRIMARY_FILE_KEY}"]
    method_with_tables = [s[len("source."):] for s in src_cp.sections() if s.startswith("source.")
                          and src_cp[s].get("role") == "method_citation" and src_cp[s].get("file.2.url")]
    if len(method_with_tables) != 1:
        raise SystemExit(f"[STOP] expected exactly one method-citation source with two tables (the iBAQ-vs-LFQ comparison), found {method_with_tables}")

    path, sha = verify_file(manifest, primary, PRIMARY_FILE_KEY, log)
    mf = manifest_file(manifest, primary, PRIMARY_FILE_KEY)
    pool = load_pool()
    mw = load_master_mw()
    rows, meta = read_dataset1(path, cols, log)
    meta["source_line"] = f"{src['id_short']} {src['citation']} DOI {src['doi']}; file {mf['path']} sha256 {sha}"
    meta["retrieved"] = mf.get("retrieved", "")

    j = join(rows, pool, log)
    if j["stops"]:
        for s in j["stops"]:
            log.error("STOP: %s", s)
        raise SystemExit("[STOP] " + "; ".join(j["stops"][:5]))
    measured = j["measured"]
    wres = weights(measured, mw, log)
    tiers = wres["tiers"]

    hashes = {
        "dataset": f"{mf['path']} sha256 {sha}",
        "accessions": f"{common.ACCESSIONS_INI.relative_to(common.REPO_ROOT).as_posix()} sha256 {sha256_path(common.ACCESSIONS_INI)}",
        "composition": f"{common.COMPOSITION_TSV.relative_to(common.REPO_ROOT).as_posix()} sha256 {sha256_path(common.COMPOSITION_TSV)}",
    }
    header_common = [f"generated   = {common.iso_now()}",
                     f"dataset     = {hashes['dataset']}",
                     f"accessions  = {hashes['accessions']}",
                     f"composition = {hashes['composition']}",
                     f"rules       = B0-B8, D31 (within-tier), D47 (summed rows), D48 (shared rows), D58 (combined)",
                     f"header_whitespace_stripped_for = {meta['header_whitespace_stripped_for']}"]

    files: dict[str, str] = {}   # relative path -> content
    # the generated Layer B config (D61)
    # weights_per_pool_entry.tsv
    hdr = ["accession", "gene", "tier", "mw_master", "match_rule", "n_rows", "rows"]
    for ft in FIBER_TYPES:
        hdr += [f"v_{ft}", f"sd_{ft}", f"valid_{ft}", f"minor_row_share_{ft}", f"shared_row_value_{ft}"]
    for t in tiers:
        for ft in FIBER_TYPES:
            hdr += [f"w_{ft}_tier{t}", f"w_{ft}_tier{t}_low", f"w_{ft}_tier{t}_high"]
    for ft in FIBER_TYPES:
        hdr += [f"w_{ft}_combined", f"w_{ft}_combined_low", f"w_{ft}_combined_high"]
    body = []
    for acc in sorted(measured):
        m = measured[acc]
        row = [acc, m["gene"], m["tier"], f"{m['mw']:.4f}", m["match_rule"], m["n_rows"],
               ";".join(str(i) for i in m["row_indices"])]
        for ft in FIBER_TYPES:
            row += [fnum(m[f"v_{ft}"]), fnum(m[f"sd_{ft}"]), m[f"valid_{ft}"], fnum(m[f"minor_share_{ft}"]), fnum(m[f"shared_whole_{ft}"])]
        for t in tiers:
            for ft in FIBER_TYPES:
                w = m.get(f"w_{ft}_tier{t}")
                row += ["" if w is None else fnum(w), "" if w is None else fnum(m[f"w_{ft}_tier{t}_low"]),
                        "" if w is None else fnum(m[f"w_{ft}_tier{t}_high"])]
        for ft in FIBER_TYPES:
            row += [fnum(m[f"w_{ft}_combined"]), fnum(m[f"w_{ft}_combined_low"]), fnum(m[f"w_{ft}_combined_high"])]
        body.append(row)
    files["weights_per_pool_entry.tsv"] = tsv_text(header_common, hdr, body)

    files["config/mass_fractions_per_entry.tsv"] = per_entry_table(measured, tiers, meta, header_common)

    # combined ranked table (D58): every entry, one denominator, with its tier
    rank_c, cum_c = {}, {}
    for ft in FIBER_TYPES:
        order = sorted(measured.values(), key=lambda m: (-m[f"w_{ft}_combined"], m["accession"]))
        running = 0.0
        for i, m in enumerate(order, 1):
            running += m[f"w_{ft}_combined"]
            rank_c[(m["accession"], ft)] = i
            cum_c[(m["accession"], ft)] = running
    hdr = ["rank_I", "accession", "gene", "tier", "match_rule"]
    for ft in FIBER_TYPES:
        hdr += [f"w_{ft}", f"rank_{ft}", f"cumulative_share_{ft}_in_{ft}_order", f"w_{ft}_low", f"w_{ft}_high"]
    body = []
    for m in sorted(measured.values(), key=lambda m: rank_c[(m["accession"], "I")]):
        row = [rank_c[(m["accession"], "I")], m["accession"], m["gene"], m["tier"], m["match_rule"]]
        for ft in FIBER_TYPES:
            row += [fnum(m[f"w_{ft}_combined"]), rank_c[(m["accession"], ft)], fnum(cum_c[(m["accession"], ft)]),
                    fnum(m[f"w_{ft}_combined_low"]), fnum(m[f"w_{ft}_combined_high"])]
        body.append(row)
    files["combined_entries_ranked.tsv"] = tsv_text(
        header_common + ["D58: the combined standard — contractile (tier 1) and support (tier 2) together, one denominator per fiber type; ordered by type-I weight"],
        hdr, body)

    # per-tier ranked tables: rank per fiber type, cumulative share down each fiber type's own order
    for tier in tiers:
        members = [m for m in measured.values() if tier in m["tiers"]]
        rank, cum = {}, {}
        for ft in FIBER_TYPES:
            order = sorted(members, key=lambda m: (-(m[f"w_{ft}_tier{tier}"] or 0.0), m["accession"]))
            running = 0.0
            for i, m in enumerate(order, 1):
                running += m[f"w_{ft}_tier{tier}"] or 0.0
                rank[(m["accession"], ft)] = i
                cum[(m["accession"], ft)] = running
        hdr = ["rank_I", "accession", "gene", "tier", "match_rule"]
        for ft in FIBER_TYPES:
            hdr += [f"w_{ft}", f"rank_{ft}", f"cumulative_share_{ft}_in_{ft}_order", f"w_{ft}_low", f"w_{ft}_high", f"valid_{ft}"]
        body = []
        for m in sorted(members, key=lambda m: rank[(m["accession"], "I")]):
            row = [rank[(m["accession"], "I")], m["accession"], m["gene"], m["tier"], m["match_rule"]]
            for ft in FIBER_TYPES:
                row += [fnum(m[f"w_{ft}_tier{tier}"]), rank[(m["accession"], ft)], fnum(cum[(m["accession"], ft)]),
                        fnum(m[f"w_{ft}_tier{tier}_low"]), fnum(m[f"w_{ft}_tier{tier}_high"]), m[f"valid_{ft}"]]
            body.append(row)
        files[f"tier{tier}_entries_ranked.tsv"] = tsv_text(
            header_common + [f"tier {tier}: {len(members)} entries, ordered by type-I weight; rank_<type> and cumulative_share_<type> follow that fiber type's own order"],
            hdr, body)

    # pool entries without a row
    body = [[a, measured[a]["gene"], measured[a]["tier"]] for a in sorted(measured) if measured[a]["match_rule"] == "none"]
    files["pool_entries_without_dataset_row.tsv"] = tsv_text(header_common, ["accession", "gene", "tier"], body)

    # dataset rows outside the pool, ranked by the largest fiber-type median (molar; MW unknown outside the pool)
    outside = sorted(j["outside"], key=lambda r: -max(r[f"v_{ft}"] for ft in FIBER_TYPES))
    tot_v = {ft: sum((0.0 if math.isnan(r[f"median_{ft}"]) else r[f"median_{ft}"]) for r in rows) for ft in FIBER_TYPES}
    hdr = ["rank", "gene", "kind", "n_rows"] + [f"v_{ft}" for ft in FIBER_TYPES] + [f"v_share_of_file_{ft}" for ft in FIBER_TYPES] + [f"valid_{ft}" for ft in FIBER_TYPES]
    body = []
    for i, r in enumerate(outside, 1):
        body.append([i, r["gene"], r["kind"], r["n_rows"]] + [fnum(r[f"v_{ft}"]) for ft in FIBER_TYPES]
                    + [fnum(r[f"v_{ft}"] / tot_v[ft] if tot_v[ft] else 0.0) for ft in FIBER_TYPES] + [r[f"valid_{ft}"] for ft in FIBER_TYPES])
    files["dataset_rows_outside_pool.tsv"] = tsv_text(header_common + ["ranking = by largest median value (molar, as published: normalised to the dataset's anchor protein); MW is known only for pool entries, so no v*MW share here"], hdr, body)

    # shared gene rows
    hdr = ["row_index", "gene_cell", "n_members", "members_in_pool"] + [f"v_{ft}" for ft in FIBER_TYPES] + ["outcome"]
    body = [[r["row_index"], r["gene_cell"], r["n_members"], r["members_in_pool"]] + [fnum(r[f"v_{ft}"]) for ft in FIBER_TYPES] + [r["outcome"]]
            for r in j["shared_table"]]
    files["shared_gene_rows.tsv"] = tsv_text(header_common, hdr, body)

    # R5-excluded entries (D59): the mass they would have carried, from their dataset value and UniProt's molecular weight
    excl_path = common.PROTEIN_SET_OUT_DIR / "excluded_non_standard_alphabet.tsv"
    excl_rows = []
    if excl_path.exists():
        outside_by_gene = {r["gene"]: r for r in j["outside"]}
        for e in common.read_tsv(excl_path):
            o = outside_by_gene.get(e["gene"])
            mw_e = float(e["uniprot_mol_weight_da"]) if e.get("uniprot_mol_weight_da") else float("nan")
            row = [e["accession"], e["gene"], e["tier"], e["non_standard_letters"], e.get("uniprot_mol_weight_da", "")]
            for ft in FIBER_TYPES:
                v_e = o[f"v_{ft}"] if o else 0.0
                vm = v_e * mw_e if not math.isnan(mw_e) else float("nan")
                shares = []
                for tier in tiers:
                    if tier in e["tier"].split(";"):
                        d = wres["denom"][(tier, ft)]
                        shares.append(f"tier{tier}:" + fnum(vm / (d + vm)) if not math.isnan(vm) and d + vm > 0 else f"tier{tier}:nan")
                d_all = wres["denom_all"][ft]
                row += [fnum(v_e), fnum(vm), ";".join(shares), fnum(vm / (d_all + vm)) if not math.isnan(vm) and d_all + vm > 0 else "nan"]
            excl_rows.append(row)
    hdr = ["accession", "gene", "tier", "non_standard_letters", "uniprot_mol_weight_da"]
    for ft in FIBER_TYPES:
        hdr += [f"v_{ft}", f"v_x_mw_{ft}", f"share_of_tier_if_included_{ft}", f"share_of_combined_if_included_{ft}"]
    files["excluded_entries_mass_share.tsv"] = tsv_text(
        header_common + ["R5 (D59): entries excluded for letters outside the twenty coded amino acids; v from the dataset row, MW from UniProt's sequence.molWeight; share = v*MW / (denominator + v*MW)"],
        hdr, excl_rows)

    # branch table
    hdr = ["tier", "term_id", "term_name", "n_entries_returned_by_term", "n_entries_only_via_branch"] + \
          [f"sum_w_{ft}_returned" for ft in FIBER_TYPES] + [f"sum_w_{ft}_only_via_branch" for ft in FIBER_TYPES] + ["entries_only_via_branch"]
    files["branch_exclusive_mass.tsv"] = tsv_text(header_common, hdr, branch_table(measured, tiers))

    # weighted bounds
    wb_rows, glyc = weighted_bounds(measured, tiers)
    files["weighted_bounds.tsv"] = tsv_text(header_common + ["bound x w: per-entry bound (isoform: largest |delta| over the entry's isoforms; processing: |delta| mature vs full; ptm: summed PTM mass / MW) multiplied by the entry's weight; largest per amino acid"],
                                            ["tier", "fiber_type", "bound", "amino_acid", "max_bound_times_w", "accession"], wb_rows)

    # gel cross-check (D39, D50, D52) and the third measurement (D54)
    fam_rows, check_rows, band_families, carroll_cp = carroll_check(measured, pool, log)
    share_rows, deshmukh_ratio_rows = intensity_share_check(dec, manifest, band_families, carroll_cp, pool, measured, log)
    check_rows += share_rows
    acc_src = [s for s in dec.sections() if s.startswith("columns.") and dec[s].get("identity_kind") == "accession_list"][0].split(".")[1]
    files[f"ratio_check_{acc_src}.tsv"] = tsv_text(
        header_common + ["per pool entry: intensity in the slow and fast sheets (joined by accession), slow/fast, the primary dataset's I/IIa, and their quotient; fast pools ~80/20 IIa/IIx per config"],
        ["accession", "gene", "tier", "intensity_slow", "intensity_fast", "slow_over_fast", "primary_I_over_IIa", "primary_over_source"],
        deshmukh_ratio_rows)

    # slow/fast ratio check from the log2 per-participant source (D34)
    sf_rows, sf_summ = slow_fast_ratio_check(dec, manifest, measured, log)
    files[f"ratio_check_{sf_summ['source']}.tsv"] = tsv_text(
        header_common + [f"source file {sf_summ['file']}", "per pool entry: mean and SD over participants of log2(slow) - log2(fast); slow_over_fast = 2^mean; quotients against the primary dataset's I/IIa and I/IIx"],
        ["accession", "gene", "tier", "primary_v_I", "primary_v_IIa", "primary_v_IIx", "primary_I_over_IIa", "primary_I_over_IIx",
         "n_participants", "log2_slow_minus_fast_mean", "log2_slow_minus_fast_sd", "slow_over_fast", "primary_I_over_IIa_over_source", "primary_I_over_IIx_over_source"],
        sf_rows)

    # family bounds (D55)
    fb_rows = family_bounds(measured, tiers, 2)
    files["family_bounds.tsv"] = tsv_text(
        header_common + ["D55: families from outputs/digest/shared_pairs.tsv at n_shared_peptides >= 2; spread = max - min of each amino acid's free-mass fraction across members; bound = spread x family weight"],
        ["cutoff", "family_id", "n_members", "genes", "tier", "fiber_type", "family_w", "max_bound", "amino_acid_at_max"] + [f"bound_{a}" for a in common.AMINO_ACIDS],
        fb_rows)
    files["band_families.tsv"] = tsv_text(header_common + ["families from outputs/digest/shared_pairs.tsv at n_shared_peptides >= cutoff, anchored per config/carroll_classical_fractionation.ini (D52)"],
                                          ["cutoff", "band", "accession", "gene", "tier", "role", "match_rule"] + [f"v_{ft}" for ft in FIBER_TYPES] + [f"w_{ft}_tier1" for ft in FIBER_TYPES], fam_rows)
    files["classical_check_carroll_2004.tsv"] = tsv_text(
        header_common + ["gel values: config/carroll_classical_fractionation.ini (D50); band = the D52 family",
                         "method = ibaq_x_mw_within_tier1: mhc and actin are sums of Tier 1 weights over the family;",
                         "method = intensity_share_<source>: mhc and actin are sums of that file's intensity column over the family, per its sheet (slow / fast);",
                         "quotient = mhc_to_actin_method / mhc_to_actin_gel. Reported per D54 without attribution of the difference to any method."],
        ["cutoff", "fiber_type", "method", "mhc", "actin", "mhc_to_actin_method", "gel_mhc", "gel_actin", "gel_unit",
         "mhc_to_actin_gel", "method_over_gel"],
        check_rows)

    # iBAQ vs MaxLFQ on the same fibers (the method comparison, D35)
    ibaq_rows, ibaq_note = ibaq_vs_lfq_table(manifest, method_with_tables[0], log)
    files[f"{method_with_tables[0]}_myh_fractions_ibaq_vs_lfq.tsv"] = tsv_text(
        header_common + ["file.1 = iBAQ fractions, file.2 = MaxLFQ fractions (as published); mean over fibers", ibaq_note],
        ["gene", "n_fibers_in_both", "mean_fraction_ibaq_common", "mean_fraction_lfq_common", "lfq_over_ibaq",
         "n_fibers_ibaq", "mean_fraction_ibaq_all", "n_fibers_lfq", "mean_fraction_lfq_all"],
        ibaq_rows)

    # generated config per fiber type (kept; the table above is an added format, D61)
    for ft in FIBER_TYPES:
        files[f"config/mass_fractions/{ft}.ini"] = render_fiber_config(ft, measured, tiers, meta, header_common)

    # summary
    summ = common.new_ini()
    summ.add_section("dataset")
    summ["dataset"]["file"] = hashes["dataset"]
    for k, v in j["counts"].items():
        summ["dataset"][k] = str(v)
    summ["dataset"]["header_whitespace_stripped_for"] = str(meta["header_whitespace_stripped_for"])
    summ.add_section("pool")
    summ["pool"]["entries"] = str(len(measured))
    summ["pool"]["entries_without_row"] = str(sum(1 for m in measured.values() if m["match_rule"] == "none"))
    summ["pool"]["entries_summed"] = str(sum(1 for m in measured.values() if m["match_rule"].startswith("summed")))
    summ["pool"]["entries_touched_by_shared_rows"] = str(sum(1 for m in measured.values() if "shared" in m["match_rule"]))
    summ["pool"]["dataset_genes_outside_pool"] = str(len(j["outside"]))
    summ.add_section("weights")
    for t in tiers:
        for ft in FIBER_TYPES:
            s = sum(m[f"w_{ft}_tier{t}"] or 0.0 for m in measured.values())
            summ["weights"][f"sum_w_{ft}_tier{t}"] = fnum(s)
            summ["weights"][f"shared_row_signal_share_{ft}_tier{t}"] = fnum(wres["shared_share"][(t, ft)])
            summ["weights"][f"glycosylated_entries_sum_w_{ft}_tier{t}"] = fnum(glyc[(t, ft)])
    for tier in tiers:
        members = [m for m in measured.values() if tier in m["tiers"]]
        for ft in FIBER_TYPES:
            sec = f"tier{tier}.type_{ft}"
            summ.add_section(sec)
            ws = sorted(((m[f"w_{ft}_tier{tier}"] or 0.0), m["gene"], m["accession"]) for m in members)[::-1]
            summ[sec]["entries"] = str(len(members))
            summ[sec]["entries_with_weight_above_zero"] = str(sum(1 for w, _, _ in ws if w > 0))
            summ[sec]["share_of_ten_largest"] = fnum(sum(w for w, _, _ in ws[:10]))
            summ[sec]["largest_entry"] = f"{ws[0][2]} {ws[0][1]} {fnum(ws[0][0])}" if ws else ""
            for kind in ("isoform", "processing", "ptm_mass"):
                best = max((r for r in wb_rows if r[0] == tier and r[1] == ft and r[2] == kind), key=lambda r: float(r[4]), default=None)
                if best:
                    summ[sec][f"max_weighted_{kind}_bound"] = f"{best[4]} ({best[3]}, {best[5]})"
            summ[sec]["glycosylated_entries_sum_w"] = fnum(glyc[(tier, ft)])
    for tier in tiers:
        for ft in FIBER_TYPES:
            sec = f"tier{tier}.type_{ft}"
            best = max((r for r in fb_rows if r[4] == tier and r[5] == ft), key=lambda r: float(r[7]), default=None)
            if best:
                summ[sec]["max_family_bound"] = f"{best[7]} ({best[8]}, {best[3]})"
    summ.add_section("combined")
    for ft in FIBER_TYPES:
        summ["combined"][f"sum_w_{ft}"] = fnum(sum(m[f"w_{ft}_combined"] for m in measured.values()))
        for tier in tiers:
            summ["combined"][f"share_tier{tier}_{ft}"] = fnum(sum(m[f"w_{ft}_combined"] for m in measured.values() if tier in m["tiers"]))
        ws = sorted((m[f"w_{ft}_combined"], m["accession"], m["gene"]) for m in measured.values())[::-1]
        summ["combined"][f"entries_with_weight_above_zero_{ft}"] = str(sum(1 for w, _, _ in ws if w > 0))
        summ["combined"][f"share_of_ten_largest_{ft}"] = fnum(sum(w for w, _, _ in ws[:10]))
    summ["combined"]["r5_excluded_entries"] = str(len(excl_rows))
    summ.add_section("slow_fast_ratio_check")
    for k, v in sf_summ.items():
        summ["slow_fast_ratio_check"][k] = str(v)
    summ.add_section("mhc_to_actin")
    for r in check_rows:
        summ["mhc_to_actin"][f"cutoff{r[0]}_{r[1].split(' ')[0]}_{r[2]}_over_gel"] = r[10]
    files["mass_fractions_summary.ini"] = common.render_ini(summ, ["GENERATED by mass_fractions.py. DO NOT EDIT BY HAND."] + header_common)
    return files


def tsv_text(header_lines: list[str], columns: list[str], rows: list[list]) -> str:
    out = ["# GENERATED by mass_fractions.py. DO NOT EDIT BY HAND."] + [f"# {h}" for h in header_lines]
    out.append("\t".join(columns))
    for r in rows:
        out.append("\t".join("" if v is None else str(v) for v in r))
    return "\n".join(out) + "\n"


def target_path(rel: str) -> Path:
    return (common.REPO_ROOT / rel) if rel.startswith("config/") else (OUT_DIR / rel)


def strip_generated_line(text: str) -> str:
    return "\n".join(l for l in text.splitlines() if not l.startswith("# generated"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger("mass_fractions_check" if args.check else "mass_fractions")
    files = build(log)
    if args.check:
        stale = []
        for rel, content in files.items():
            p = target_path(rel)
            if not p.exists() or strip_generated_line(p.read_text(encoding="utf-8")) != strip_generated_line(content):
                stale.append(rel)
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), stale)
            return 1
        log.info("check: mass-fraction outputs and config are current")
        return 0
    for rel, content in files.items():
        p = target_path(rel)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content, encoding="utf-8", newline="\n")
        log.info("wrote %s", p.relative_to(common.REPO_ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
