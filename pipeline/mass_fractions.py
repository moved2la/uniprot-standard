#!/usr/bin/env python3
"""
mass_fractions.py — offline stage: Layer B weights from the primary dataset.

Reads
  config/mass_fraction_decisions.ini      the [source.*] with role = primary and its [columns.*] map;
                                          the method-citation source with two tables (iBAQ vs LFQ)
  data/literature/manifest.ini            path + sha256 of every literature file (B0)
  config/accessions.ini                   the pool: accession, gene, tier
  outputs/composition/amino_acid_composition_per_protein.tsv   mw on segment_set = master (B2)
  data/tier<N>_pool.tsv, data/gene-ontology/tier<N>_subtree.tsv   branch table
  outputs/isoform_deltas.tsv, outputs/composition/processing_mass_deltas.tsv,
  outputs/composition/ptm_mass_deltas.tsv                       weighted bounds
  outputs/digest/shared_pairs.tsv                                gel-band families (D52)
  config/carroll_classical_fractionation.ini                    gel values (D50)

Writes (outputs/mass_fractions/, every file with a header naming inputs and hashes)
  weights_per_pool_entry.tsv           one row per pool accession: what the dataset says about
                                       it, how the row was matched, and its weight per tier and
                                       fiber type — the human-readable "Tier 1 as measured"
  pool_entries_without_dataset_row.tsv pool accessions with no Dataset 1 row (w = 0)
  dataset_rows_outside_pool.tsv        Dataset 1 genes matching no pool entry, ranked by median
                                       value — the completeness check and the Tier 2 evidence
  shared_gene_rows.tsv                 every ';' row touching the pool with its D48 outcome
  branch_exclusive_mass.tsv            per ontology term: sum of w over entries reached only
                                       through that term's branch
  weighted_bounds.tsv                  isoform / processing / PTM bounds x w, largest per amino
                                       acid; sum of w over glycosylated entries
  classical_check_carroll_2004.tsv     gel MHC:actin vs sum of w over the D52 band families
  band_families.tsv                    the entries each band family contains, at cutoffs 1 and 2
  <source>_myh_fractions_ibaq_vs_lfq.tsv   per gene of the two tables, mean fraction per fiber under
                                       iBAQ and MaxLFQ, over the fibers present in both tables
  mass_fractions_summary.ini
  config/mass_fractions/{I,IIa,IIx}.ini            generated; one section per accession
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

DECISIONS_INI = common.CONFIG_DIR / "mass_fraction_decisions.ini"
CARROLL_INI = common.CONFIG_DIR / "carroll_classical_fractionation.ini"
MANIFEST_INI = common.DATA_DIR / "literature" / "manifest.ini"
MASS_FRACTIONS_CONFIG_DIR = common.CONFIG_DIR / "mass_fractions"
OUT_DIR = common.OUTPUTS_DIR / "mass_fractions"
DIGEST_DIR = common.OUTPUTS_DIR / "digest"
ISOFORM_DELTAS_TSV = common.OUTPUTS_DIR / "isoform_deltas.tsv"
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
        if g in per_gene:
            pg = per_gene[g]
            rec["match_rule"] = "single" if pg["n_rows"] == 1 else f"summed:{pg['n_rows']}"
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
    shared_share = {}
    for t in tiers:
        for ft in FIBER_TYPES:
            d = denom[(t, ft)]
            s = sum(m[f"shared_add_{ft}"] * m["mw"] for m in measured.values() if t in m["tiers"])
            shared_share[(t, ft)] = s / d if d > 0 else 0.0
    log.info("denominators (sum v*MW per tier and fiber type): %s", {f"{t}/{ft}": f"{v:.6g}" for (t, ft), v in denom.items()})
    return {"tiers": tiers, "denom": denom, "shared_share": shared_share}


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
            for a in sorted(accs):
                m = measured[a]
                fam_rows.append([cutoff, band, a, m["gene"], m["tier"], "anchor" if m["gene"] in anchors else "family",
                                 fnum(m.get("w_I_tier1") or 0.0), fnum(m.get("w_IIa_tier1") or 0.0)])
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
            check_rows.append([cutoff, ft, fnum(wm), fnum(wa), fnum(our_ratio),
                               f"{m_val}±{m_sd}", f"{a_val}±{a_sd}", unit, fnum(gel_ratio),
                               fnum(our_ratio / gel_ratio) if gel_ratio else "nan",
                               fnum(m_val / 1000.0 / wm) if wm > 0 else "nan",
                               fnum(a_val / 1000.0 / wa) if wa > 0 else "nan"])
            log.info("Carroll cutoff %d type %s: MHC:actin ours %.4f gel %.4f", cutoff, ft, our_ratio, gel_ratio)
    return fam_rows, check_rows


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
        fibers = header[1:]
        data = {}
        for row in it:
            if row[0] is None:
                continue
            data[str(row[0]).strip()] = {f: (float(v) if v is not None and str(v).strip() != "" else float("nan"))
                                         for f, v in zip(fibers, row[1:])}
        wb.close()
        tables[label] = (fibers, data)
    fib_common = [f for f in tables["ibaq"][0] if f in set(tables["lfq"][0]) and f]
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
                     len(ib), fnum(mean(ib.values())), len(lf), fnum(mean(lf.values()))])
    log.info("iBAQ-vs-LFQ tables: %d genes; %d fibers in both tables", len(genes), len(fib_common))
    return rows


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
        s["v"] = fnum(m[f"v_{ft}"])
        s["sd"] = fnum(m[f"sd_{ft}"])
        s["valid_values"] = str(m[f"valid_{ft}"])
        s["mw_master"] = f"{m['mw']:.4f}"
        s["match_rule"] = m["match_rule"]
        s["rows"] = ";".join(str(i) for i in m["row_indices"]) if m["row_indices"] else "none"
        s["source"] = meta["source_line"]
        s["location"] = (f"sheet {meta['sheet']}, header row {meta['header_row']}, columns "
                         f"{meta['columns'][f'median_{ft}']!r} / {meta['columns'][f'sd_{ft}']!r} / {meta['columns'][f'valid_{ft}']!r}, rows as listed")
        s["retrieved"] = meta["retrieved"]
    header = [f"GENERATED by mass_fractions.py. DO NOT EDIT BY HAND.", f"fiber_type = {ft} (the dataset's own column, B3)"] + header_common
    return common.render_ini(cp, header)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def build(log) -> dict:
    dec = common.read_ini(DECISIONS_INI)
    manifest = common.read_ini(MANIFEST_INI)
    primaries = [s for s in dec.sections() if s.startswith("source.") and dec[s].get("role") == "primary"]
    if len(primaries) != 1:
        raise SystemExit(f"[STOP] expected exactly one [source.*] with role = primary, found {primaries}")
    primary = primaries[0][len("source."):]
    src = dec[f"source.{primary}"]
    cols = dec[f"columns.{primary}.{PRIMARY_FILE_KEY}"]
    method_with_tables = [s[len("source."):] for s in dec.sections() if s.startswith("source.")
                          and dec[s].get("role") == "method_citation" and dec[s].get("file.2.url")]
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
                     f"rules       = B0-B5, D31 (within-tier), D47 (summed rows), D48 (shared rows)",
                     f"header_whitespace_stripped_for = {meta['header_whitespace_stripped_for']}"]

    files: dict[str, str] = {}   # relative path -> content
    # weights_per_pool_entry.tsv
    hdr = ["accession", "gene", "tier", "mw_master", "match_rule", "n_rows", "rows"]
    for ft in FIBER_TYPES:
        hdr += [f"v_{ft}", f"sd_{ft}", f"valid_{ft}", f"minor_row_share_{ft}", f"shared_row_value_{ft}"]
    for t in tiers:
        for ft in FIBER_TYPES:
            hdr += [f"w_{ft}_tier{t}", f"w_{ft}_tier{t}_low", f"w_{ft}_tier{t}_high"]
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
        body.append(row)
    files["weights_per_pool_entry.tsv"] = tsv_text(header_common, hdr, body)

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

    # branch table
    hdr = ["tier", "term_id", "term_name", "n_entries_returned_by_term", "n_entries_only_via_branch"] + \
          [f"sum_w_{ft}_returned" for ft in FIBER_TYPES] + [f"sum_w_{ft}_only_via_branch" for ft in FIBER_TYPES] + ["entries_only_via_branch"]
    files["branch_exclusive_mass.tsv"] = tsv_text(header_common, hdr, branch_table(measured, tiers))

    # weighted bounds
    wb_rows, glyc = weighted_bounds(measured, tiers)
    files["weighted_bounds.tsv"] = tsv_text(header_common + ["bound x w: per-entry bound (isoform: largest |delta| over the entry's isoforms; processing: |delta| mature vs full; ptm: summed PTM mass / MW) multiplied by the entry's weight; largest per amino acid"],
                                            ["tier", "fiber_type", "bound", "amino_acid", "max_bound_times_w", "accession"], wb_rows)

    # gel cross-check (D39, D50, D52)
    fam_rows, check_rows = carroll_check(measured, pool, log)
    files["band_families.tsv"] = tsv_text(header_common + ["families from outputs/digest/shared_pairs.tsv at n_shared_peptides >= cutoff, anchored per config/carroll_classical_fractionation.ini (D52)"],
                                          ["cutoff", "band", "accession", "gene", "tier", "role", "w_I_tier1", "w_IIa_tier1"], fam_rows)
    files["classical_check_carroll_2004.tsv"] = tsv_text(
        header_common + ["gel values: config/carroll_classical_fractionation.ini (D50); band = sum of Tier 1 w over the D52 family",
                         "implied_tier1_fraction_of_total_protein = gel fraction / our within-tier weight, per band — the Tier 1 : total protein ratio the gel implies"],
        ["cutoff", "fiber_type", "sum_w_mhc_band", "sum_w_actin_band", "mhc_to_actin_ours", "gel_mhc", "gel_actin", "gel_unit",
         "mhc_to_actin_gel", "ours_over_gel", "implied_tier1_fraction_of_total_protein_from_mhc", "implied_tier1_fraction_of_total_protein_from_actin"],
        check_rows)

    # iBAQ vs MaxLFQ on the same fibers (the method comparison, D35)
    files[f"{method_with_tables[0]}_myh_fractions_ibaq_vs_lfq.tsv"] = tsv_text(
        header_common + ["file.1 = iBAQ fractions (Table S2), file.2 = MaxLFQ fractions (Table S3); mean over fibers"],
        ["gene", "n_fibers_in_both", "mean_fraction_ibaq_common", "mean_fraction_lfq_common", "lfq_over_ibaq",
         "n_fibers_ibaq", "mean_fraction_ibaq_all", "n_fibers_lfq", "mean_fraction_lfq_all"],
        ibaq_vs_lfq_table(manifest, method_with_tables[0], log))

    # generated config per fiber type
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
    summ.add_section("carroll_2004")
    for r in check_rows:
        summ["carroll_2004"][f"cutoff{r[0]}_type_{r[1]}_mhc_to_actin_ours_over_gel"] = r[9]
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
