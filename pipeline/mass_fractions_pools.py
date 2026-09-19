"""Layer B weights for a category whose pools are PUBLISHED MASS SHARES (Step 7b, blood).

The muscle stage (mass_fractions.py) joins one dataset of three fiber-type medians to the pool by
gene symbol and turns molar-like values into mass shares through the molecular weight. A blood
pool is different in every one of those respects: the join was done by the pool enumeration
(rules P1-P6, D107/D111 -- data/<category>/<pool>_dataset_rows.tsv is the join table), each pool
has its own dataset, and each dataset already publishes a mass share (D105: `mass_share`, used
as-is, no molecular-weight step). Hence a separate module for the ad-hoc form of Step 7b; the
adaptation list records that Step 7c unifies the two.

What it computes, per pool
  the value v of each pool entry from the column map (by letter) of its primary dataset:
    share_columns   -- one mass share per replicate (Bryk: four donors); v = the MEAN across
                       replicates (D108), sd = the sample standard deviation across them, n = count
    abundance_column, abundance_scale = log10 -- (Geyer) v = 10^x; the replicate_columns (log10 per
                       run) give sd in linear space, n = runs present; the published CV is carried
  rows of one entry are summed (D47); a row shared between entries (P3b/P4b shared, `;` in the
  entry cell) is split in proportion to the members' own rows, equally if they have none, and its
  whole value goes to each member's _high (D48); sd of summed rows is root-sum-square (D53a)
  w = v / sum over the pool's entries (within pool, D31): the pool's listing renormalised after the
  contaminant rule removed its rows; w_low / w_high = (v -/+ sd) over the unchanged denominator (B5)
  bounds: the weighted isoform / processing / PTM bounds (D55 pattern) from the category's own
  delta tables; the shares the rules removed (D98) and the shares no entry carries (A7)
  D113: the immunoglobulin variable-domain bound from the pool's own V-region entries
  cross-checks, reported not used: the Hortin 2008 ratio table and the plasma profile under Hortin's
  weights; the share of Geyer's deep dataset (Table S5) that Table S2 does not carry; the mass share
  of entries present in more than one pool (plasma carry-over into the erythrocyte preparation)

Reads   config/<category>/protein_set_decisions.ini, config/<category>/mass_fraction_decisions.ini,
        config/<category>/accessions.ini, data/<category>/<pool>_dataset_rows.tsv and
        <pool>_rows_excluded_with_share.tsv and pool_<pool>.tsv, the primary tables (hashed, B0),
        outputs/<category>/intermediate/composition/*, .../protein_set/isoform_deltas.tsv, pool_overlap.tsv
Writes  config/<category>/mass_fractions_per_entry.tsv (generated config, D61)
        outputs/<category>/intermediate/mass_fractions/*
"""

from __future__ import annotations

import argparse
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from pipeline import common
from pipeline.enumerate_pool import classify_token, gene_tokens, read_sheet_by_letter, split_letters

STAGE = "mass_fractions_pools"
SIG = 10
AA = list(common.AMINO_ACIDS)
MEMBER_OUTCOMES = ("member", "member_via_later_token", "member_via_secondary_accession",
                   "member_via_secondary_accession_shared", "member_via_gene", "member_via_gene_shared",
                   "member_duplicate_entry")
UNIPROT_ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def fnum(x: float | None) -> str:
    if x is None:
        return ""
    if x == 0:
        return "0"
    return f"{x:.{SIG}g}"


def tsv_text(header_lines: list[str], columns: list[str], rows: list[list]) -> str:
    out = [f"# {l}" for l in header_lines]
    out.append("\t".join(columns))
    for r in rows:
        out.append("\t".join("" if c is None else str(c) for c in r))
    return "\n".join(out) + "\n"


def strip_generated_line(text: str) -> str:
    return "\n".join(l for l in text.splitlines() if not l.startswith("# generated"))


def cell_float(v) -> float | None:
    if v is None:
        return None
    s = str(v).strip().replace(",", "")
    if s == "" or s.lower() in ("nan", "na", "n/a"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def read_xls_by_letter(path: Path, sheet: str, header_row: int, letters: dict[str, str], log) -> list[dict]:
    """The legacy-.xls twin of enumerate_pool.read_sheet_by_letter (xlrd; Hortin 2008)."""
    import xlrd
    from openpyxl.utils import column_index_from_string
    book = xlrd.open_workbook(str(path))
    if sheet not in book.sheet_names():
        raise SystemExit(f"[STOP] sheet {sheet!r} not in {path.name}: {book.sheet_names()}")
    ws = book.sheet_by_name(sheet)
    idx = {k: column_index_from_string(v.strip()) - 1 for k, v in letters.items()}
    rows = []
    for r in range(header_row, ws.nrows):
        vals = ws.row_values(r)
        if all(v is None or str(v).strip() == "" for v in vals):
            continue
        rec = {"row_index": r + 1}
        for k, j in idx.items():
            rec[k] = vals[j] if j < len(vals) else None
        rows.append(rec)
    log.info("%s sheet %r: %d non-empty rows below row %d (xlrd)", path.name, sheet, len(rows), header_row)
    return rows


def read_table(path: Path, sheet: str, header_row: int, letters: dict[str, str], log) -> list[dict]:
    if path.suffix.lower() == ".xls":
        return read_xls_by_letter(path, sheet, header_row, letters, log)
    return read_sheet_by_letter(path, sheet, header_row, letters, log)


def verify_file(manifest, source_id: str, file_key: str) -> tuple[Path, str]:
    """B0: the file is re-hashed against the manifest before it is read."""
    sec = f"file.{source_id}.{file_key}"
    if not manifest.has_section(sec):
        raise SystemExit(f"[STOP] {sec} not in the literature manifest")
    path = common.REPO_ROOT / manifest[sec]["path"]
    if not path.exists():
        raise SystemExit(f"[STOP] {path} missing on disk")
    got = common.sha256_file(path)
    if got != manifest[sec]["sha256"].lower():
        raise SystemExit(f"[STOP] {path.name}: sha256 {got} != manifest {manifest[sec]['sha256']} (B0)")
    return path, got


def sample_sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


# ----------------------------------------------------------------------------- inputs

def load_membership(files) -> dict[str, dict]:
    cp = common.read_ini(files["accessions"])
    out = {}
    for acc in cp.sections():
        s = cp[acc]
        out[acc] = {"gene": s.get("gene", "").strip(), "pools": [p for p in s.get("tier", "").split(";") if p],
                    "flag_open": s.get("flag_open") == "true"}
    return out


def load_composition(files) -> dict[str, dict]:
    out = {}
    for r in common.read_tsv(files["composition_tsv"]):
        if r["segment_set"] == "master":
            out[r["accession"]] = {"mw": float(r["mw"]), "n_residues": int(r["n_residues"]),
                                   "free_frac": {a: float(r[f"free_frac_{a}"]) for a in AA},
                                   "residue_frac": {a: float(r[f"residue_frac_{a}"]) for a in AA}}
    return out


# ----------------------------------------------------------------------------- one pool

def pool_values(category: str, pool: str, psd, dec, manifest, membership, log) -> dict:
    """The value of every entry of one pool from its dataset, joined through the enumeration's row table."""
    sec = psd[f"pool.{pool}"]
    src, fkey = sec["dataset_source"], sec["dataset_file"]
    cols = dec[f"columns.{src}.{fkey}"]
    path, sha = verify_file(manifest, src, fkey)
    letters = {}
    reps = split_letters(cols.get("share_columns", ""))
    for k, letter in enumerate(reps):
        letters[f"rep_{k}"] = letter
    if cols.get("abundance_column"):
        letters["abundance"] = cols["abundance_column"]
    runs = split_letters(cols.get("replicate_columns", ""))
    for k, letter in enumerate(runs):
        letters[f"run_{k}"] = letter
    if cols.get("cv_column"):
        letters["cv"] = cols["cv_column"]
    sheet_rows = {r["row_index"]: r for r in read_table(path, cols["sheet"], int(cols["header_row"]), letters, log)}
    scale = 0.01 if cols.get("share_scale", "").strip().lower() == "percent" else 1.0
    log10 = cols.get("abundance_scale", "").strip().lower() == "log10"
    if cols.get("share_kind", "").strip() != "mass_share":
        raise SystemExit(f"[STOP] [columns.{src}.{fkey}] share_kind must be mass_share (D105: the molar kind is not built until a source needs it)")

    def row_value(rec) -> tuple[float, list[float], list[float | None], float | None]:
        """(v, per-replicate shares or [], per-run linear values or [], published cv or None)."""
        if reps:
            vals = [(cell_float(rec.get(f"rep_{k}")) or 0.0) * scale for k in range(len(reps))]
            return sum(vals) / len(vals), vals, [], None
        x = cell_float(rec.get("abundance"))
        v = 0.0 if x is None else (10 ** x if log10 else x)
        run_vals = []
        for k in range(len(runs)):
            y = cell_float(rec.get(f"run_{k}"))
            run_vals.append(None if y is None else (10 ** y if log10 else y))
        return v, [], run_vals, cell_float(rec.get("cv"))

    join_rows = common.read_tsv(common.category_data_dir(category) / f"{pool}_dataset_rows.tsv")
    entries: dict[str, dict] = {}
    for acc, m in membership.items():
        if pool in m["pools"] and not m["flag_open"]:
            entries[acc] = {"rows": [], "own_rows": [], "shared_rows": [], "v": 0.0, "rep": [0.0] * len(reps),
                            "run": [0.0] * len(runs), "run_present": [True] * len(runs), "cv": None,
                            "shared_add": 0.0, "shared_whole": 0.0, "sd_sq": 0.0, "published_share": 0.0}
    shared_rows: list[tuple[dict, list[str]]] = []
    listing_total = 0.0
    for jr in join_rows:
        rec = sheet_rows.get(int(jr["row_index"]))
        if rec is None:
            raise SystemExit(f"[STOP] {pool}: dataset row {jr['row_index']} of the join table is not in the sheet (rerun blood-protein-set)")
        v, rep_vals, run_vals, cv = row_value(rec)
        listing_total += v
        if jr["outcome"] not in MEMBER_OUTCOMES or not jr["entry"]:
            continue
        accs = [a for a in jr["entry"].split(";") if a in entries]          # an R5-excluded entry is not in accessions.ini
        if not accs:
            continue
        if len(accs) > 1:
            shared_rows.append(({"v": v, "rep": rep_vals, "run": run_vals, "row_index": jr["row_index"]}, accs))
            continue
        e = entries[accs[0]]
        e["rows"].append(jr["row_index"]); e["own_rows"].append(jr["row_index"])
        e["v"] += v
        for k, x in enumerate(rep_vals):
            e["rep"][k] += x
        for k, y in enumerate(run_vals):
            if y is None:
                e["run_present"][k] = False
            else:
                e["run"][k] += y
        if cv is not None and e["cv"] is None:
            e["cv"] = cv
    # D48: shared rows split in proportion to the members' own values, equally if none; the whole value to _high
    shared_table = []
    for rec, accs in shared_rows:
        own = [entries[a]["v"] for a in accs]
        tot = sum(own)
        shares = [o / tot for o in own] if tot > 0 else [1.0 / len(accs)] * len(accs)
        for a, s in zip(accs, shares):
            e = entries[a]
            e["rows"].append(rec["row_index"]); e["shared_rows"].append(rec["row_index"])
            e["shared_add"] += rec["v"] * s
            e["shared_whole"] += rec["v"]
            for k, x in enumerate(rec["rep"]):
                e["rep"][k] += x * s
            for k, y in enumerate(rec["run"]):
                if y is None:
                    e["run_present"][k] = False
                else:
                    e["run"][k] += y * s
        shared_table.append([pool, rec["row_index"], ";".join(accs), fnum(rec["v"]), " ".join(f"{a}={s:.3f}" for a, s in zip(accs, shares))])
    # per entry: v, sd, n (D108 mean of replicates; sd over replicates; summed rows carry the replicate sums)
    for acc, e in entries.items():
        e["v"] += e["shared_add"]
        if reps:
            e["sd"] = sample_sd(e["rep"])
            e["n"] = len(e["rep"])
        elif runs:
            present = [x for x, ok in zip(e["run"], e["run_present"]) if ok]
            e["sd"] = sample_sd(present) if e["rows"] else 0.0
            e["n"] = len(present) if e["rows"] else 0
        else:
            e["sd"], e["n"] = 0.0, 0
        e["published_share"] = e["v"] / listing_total if listing_total else 0.0
        e["match_rule"] = ("none" if not e["rows"] else
                           ("single" if len(e["own_rows"]) == 1 else f"summed:{len(e['own_rows'])}" if e["own_rows"] else "none")
                           + ("+shared_split" if e["shared_rows"] else ""))
        e["match_rule"] = e["match_rule"].replace("none+", "")
    denom = sum(e["v"] for e in entries.values())
    for acc, e in entries.items():
        e["w"] = e["v"] / denom if denom > 0 else 0.0
        e["w_low"] = max(0.0, e["v"] - e["sd"]) / denom if denom > 0 else 0.0
        e["w_high"] = (e["v"] + e["sd"] + (e["shared_whole"] - e["shared_add"])) / denom if denom > 0 else 0.0
    log.info("pool %s: %d entries, %d with rows, %d zero; listing total %.6g, denominator %.6g (%.4f of the listing)",
             pool, len(entries), sum(1 for e in entries.values() if e["rows"]), sum(1 for e in entries.values() if e["v"] == 0),
             listing_total, denom, denom / listing_total if listing_total else 0)
    meta = {"source": src, "file_key": fkey, "path": manifest[f"file.{src}.{fkey}"]["path"], "sha256": sha,
            "sheet": cols["sheet"], "header_row": cols["header_row"], "reps": reps, "runs": runs, "cols": cols,
            "listing_total": listing_total, "denominator": denom, "n_reps": len(reps), "n_runs": len(runs)}
    return {"entries": entries, "shared_table": shared_table, "meta": meta}


# ----------------------------------------------------------------------------- bounds

def weighted_bounds(pools: dict[str, dict], files) -> list[list]:
    iso: dict[str, dict[str, float]] = defaultdict(lambda: {a: 0.0 for a in AA})
    p = files["protein_set_out"] / "isoform_deltas.tsv"
    if p.exists():
        for r in common.read_tsv(p):
            for a in AA:
                iso[r["accession"]][a] = max(iso[r["accession"]][a], abs(float(r[f"delta_{a}"])))
    proc = {r["accession"]: {a: abs(float(r[f"delta_free_frac_{a}"])) for a in AA}
            for r in common.read_tsv(files["composition_dir"] / "processing_mass_deltas.tsv")}
    ptm = {r["accession"]: r for r in common.read_tsv(files["composition_dir"] / "ptm_mass_deltas.tsv")}
    rows = []
    for pool, pd in pools.items():
        ent = pd["entries"]
        for kind, table in (("isoform", iso), ("processing", proc)):
            for a in AA:
                best, where = 0.0, ""
                for acc, e in ent.items():
                    if e["w"]:
                        d = table.get(acc, {}).get(a, 0.0) * e["w"]
                        if d > best:
                            best, where = d, acc
                rows.append([pool, kind, a, fnum(best), where])
        best, where, glyc = 0.0, "", 0.0
        for acc, e in ent.items():
            if e["w"] and acc in ptm:
                d = float(ptm[acc]["ptm_mass_delta_fraction_of_mw"]) * e["w"]
                if d > best:
                    best, where = d, acc
                if int(ptm[acc]["n_glycosylation"]) > 0:
                    glyc += e["w"]
        rows.append([pool, "ptm_mass", "all", fnum(best), where])
        rows.append([pool, "glycosylated_entries_sum_w", "all", fnum(glyc), ""])
    return rows


def excluded_shares(category: str, pool: str, pd: dict, membership, log) -> tuple[list[list], dict]:
    """What the listing carries that the weights do not: the rows the contaminant rule removed (D98,
    on purpose), the rows no entry answers to and the R5-excluded entries (A7: a completeness gap)."""
    dat = common.category_data_dir(category)
    rows, removed, gap = [], 0.0, 0.0
    for r in common.read_tsv(dat / f"{pool}_rows_excluded_with_share.tsv"):
        share = float(r["published_share_percent"]) / 100 if r["published_share_percent"] else 0.0
        kind = "completeness_gap" if r["outcome"] == "no_reviewed_human_entry" else "removed_by_rule"
        rows.append([pool, r["row_index"], r["identity_cell"], r["gene_cell"], r["entry"], r["outcome"], kind, fnum(share)])
        if kind == "removed_by_rule":
            removed += share
        else:
            gap += share
    in_pool_file = {r["accession"] for r in common.read_tsv(common.pool_file(category, pool))}
    absent = sorted(a for a in in_pool_file if a not in membership or membership[a]["flag_open"])
    if absent:
        shares = {}
        for jr in common.read_tsv(dat / f"{pool}_dataset_rows.tsv"):
            for a in jr["entry"].split(";"):
                if a in absent and jr["published_share_percent"]:
                    shares[a] = shares.get(a, 0.0) + float(jr["published_share_percent"]) / 100 / max(1, jr["entry"].count(";") + 1)
        for a in absent:
            rows.append([pool, "", a, membership.get(a, {}).get("gene", ""), a, "not_in_accessions_ini (R5 or open flag)", "completeness_gap", fnum(shares.get(a, 0.0))])
            gap += shares.get(a, 0.0)
    log.info("pool %s: removed by rule %.5f of the listing; completeness gap %.5f (%d entries not in accessions.ini)", pool, removed, gap, len(absent))
    return rows, {"removed_by_rule": removed, "completeness_gap": gap, "entries_not_in_accessions_ini": len(absent)}


def immunoglobulin_bound(pools: dict[str, dict], membership, comp, dec, log) -> tuple[list[list], list[list]]:
    """D113: an immunoglobulin constant-region entry carries, per molecule, a variable domain no row
    measures. Its unmeasured mass is bounded by w x L_V / L_C, where L_C is the entry's master length
    and L_V the mean master length of the pool's own V-region entries of the same chain. Chains and
    segments come from the HGNC gene symbol: prefix IGH / IGK / IGL, then a segment letter (V, D, J)
    followed by a digit for a variable-region segment; anything else with the prefix is constant."""
    if not dec.has_section("immunoglobulin_bound"):
        raise SystemExit("[STOP] config/<category>/mass_fraction_decisions.ini needs [immunoglobulin_bound] (D113)")
    ig = dec["immunoglobulin_bound"]
    prefixes = [p.strip() for p in ig["chain_prefixes"].split(",") if p.strip()]
    seg_letters = {s.strip() for s in ig["variable_segment_letters"].split(",") if s.strip()}
    per_entry, per_chain = [], []

    def classify(gene: str) -> tuple[str, str]:
        g = gene_tokens(gene)[0] if gene_tokens(gene) else ""
        for p in prefixes:
            if g.startswith(p) and len(g) > len(p):
                tail = g[len(p):]
                if tail[0] in seg_letters and len(tail) > 1 and tail[1].isdigit():
                    return p, "variable_segment"
                return p, "constant"
        return "", ""

    for pool, pd in pools.items():
        ent = pd["entries"]
        classes = {acc: classify(membership[acc]["gene"]) for acc in ent}
        lv = {}
        for p in prefixes:
            lengths = [comp[acc]["n_residues"] for acc, (chain, seg) in classes.items()
                       if chain == p and seg == "variable_segment" and membership[acc]["gene"][len(p):len(p) + 1] == "V" and acc in comp]
            lv[p] = (sum(lengths) / len(lengths), len(lengths)) if lengths else (None, 0)
        totals = {p: {"constant_w": 0.0, "bound": 0.0, "loose": 0.0, "variable_measured_w": 0.0, "n_constant": 0} for p in prefixes}
        for acc, (chain, seg) in sorted(classes.items()):
            if not chain:
                continue
            e = ent[acc]
            lc = comp[acc]["n_residues"] if acc in comp else None
            if seg == "constant":
                if lv[chain][0] is not None and lc:
                    bound, how = e["w"] * lv[chain][0] / lc, f"w x L_V/L_C (L_V = mean of {lv[chain][1]} V entries)"
                else:
                    bound, how = e["w"], "loose: no V-region entry of this chain in the pool; the whole share"
                totals[chain]["constant_w"] += e["w"]; totals[chain]["bound"] += bound; totals[chain]["loose"] += e["w"]; totals[chain]["n_constant"] += 1
            else:
                bound, how = 0.0, "measured variable-region entry (recovered under D111); no bound"
                totals[chain]["variable_measured_w"] += e["w"]
            per_entry.append([pool, acc, membership[acc]["gene"], chain, seg, fnum(e["w"]), lc or "",
                              fnum(lv[chain][0]) if lv[chain][0] is not None else "", fnum(bound), how])
        for p in prefixes:
            t = totals[p]
            per_chain.append([pool, p, t["n_constant"], fnum(t["constant_w"]), fnum(lv[p][0]) if lv[p][0] is not None else "", lv[p][1],
                              fnum(t["variable_measured_w"]), fnum(t["bound"]), fnum(t["loose"])])
        per_chain.append([pool, "all", sum(t["n_constant"] for t in totals.values()), fnum(sum(t["constant_w"] for t in totals.values())), "", "",
                          fnum(sum(t["variable_measured_w"] for t in totals.values())), fnum(sum(t["bound"] for t in totals.values())),
                          fnum(sum(t["loose"] for t in totals.values()))])
        log.info("pool %s: immunoglobulin variable-domain bound %.5f of the pool (loose: %.5f)", pool,
                 sum(t["bound"] for t in totals.values()), sum(t["loose"] for t in totals.values()))
    return per_entry, per_chain


# ----------------------------------------------------------------------------- cross-checks

def profile(weights: dict[str, float], comp, key: str = "free_frac") -> dict[str, float]:
    tot = sum(weights.values())
    out = {a: 0.0 for a in AA}
    if tot <= 0:
        return out
    for acc, w in weights.items():
        for a in AA:
            out[a] += w / tot * comp[acc][key][a]
    return out


def hortin_check(dec, manifest, pools, membership, comp, log) -> dict:
    """The Carroll pattern for plasma: Hortin 2008's mg/L per polypeptide joined to the plasma pool by
    canonical UniProt accession; the ratio per protein and the plasma profile under both weightings
    over the joined set. Reported, never used (D101)."""
    cc = dec["cross_checks"]
    src, fkey, pool = cc["hortin_source"], cc["hortin_file"], cc["hortin_pool"]
    cols = dec[f"columns.{src}.{fkey}"]
    path, sha = verify_file(manifest, src, fkey)
    rows = read_table(path, cols["sheet"], int(cols["header_row"]), {"identity": cols["identity_column"], "gene": cols["gene_column"],
                                                                   "value": cols["value_column"], "name": cols.get("name_column", cols["gene_column"])}, log)
    ent = pools[pool]["entries"]
    joined, unjoined, total_value = {}, [], 0.0
    for r in rows:
        val = cell_float(r.get("value"))
        if val is None:
            continue
        total_value += val
        raw = "" if r["identity"] is None else str(r["identity"]).strip()
        token = raw.split()[0] if raw else ""
        kind, acc = classify_token(token) if token else ("blank", "")
        if not acc or not UNIPROT_ACC.match(acc):
            unjoined.append([r["row_index"], str(r.get("name") or ""), raw, str(r.get("gene") or ""), fnum(val), "no UniProt accession in the cell" if not acc else "not a UniProt accession"])
        elif acc not in ent:
            unjoined.append([r["row_index"], str(r.get("name") or ""), raw, str(r.get("gene") or ""), fnum(val), f"{acc} not in the {pool} pool"])
        else:
            joined[acc] = joined.get(acc, 0.0) + val
    hort_tot = sum(joined.values())
    geyer_joined = {a: ent[a]["v"] for a in joined}
    gey_tot = sum(geyer_joined.values())
    ratio_rows = []
    for acc in sorted(joined, key=lambda a: -joined[a]):
        sh_h = joined[acc] / hort_tot if hort_tot else 0.0
        sh_g = geyer_joined[acc] / gey_tot if gey_tot else 0.0
        ratio_rows.append([acc, membership[acc]["gene"], fnum(joined[acc]), fnum(sh_h), fnum(ent[acc]["w"]), fnum(sh_g),
                           fnum(sh_g / sh_h) if sh_h else ""])
    prof_g = profile({a: ent[a]["w"] for a in ent if a in comp and ent[a]["w"] > 0}, comp)
    prof_gj = profile({a: geyer_joined[a] for a in geyer_joined if a in comp}, comp)
    prof_h = profile({a: joined[a] for a in joined if a in comp}, comp)
    prof_rows = [[a, fnum(prof_g[a]), fnum(prof_gj[a]), fnum(prof_h[a]), fnum(prof_h[a] - prof_gj[a]),
                  fnum((prof_h[a] - prof_gj[a]) / prof_gj[a]) if prof_gj[a] else ""] for a in AA]
    summary = {"hortin_file": f"{path.name} sha256 {sha}", "hortin_rows_with_value": str(sum(1 for r in rows if cell_float(r.get("value")) is not None)),
               "hortin_rows_joined": str(len(joined)), "hortin_rows_not_joined": str(len(unjoined)),
               "hortin_summed_value_all_rows": fnum(total_value), "hortin_summed_value_joined": fnum(hort_tot),
               "geyer_share_of_pool_in_joined_set": fnum(sum(ent[a]["w"] for a in joined)),
               "max_abs_profile_difference_hortin_minus_geyer": max((abs(prof_h[a] - prof_gj[a]), a) for a in AA)[1] + " " + fnum(max(abs(prof_h[a] - prof_gj[a]) for a in AA))}
    log.info("Hortin cross-check: %d rows joined, %d not; Geyer's pool share inside the joined set %.4f; largest profile difference %s",
             len(joined), len(unjoined), sum(ent[a]["w"] for a in joined), summary["max_abs_profile_difference_hortin_minus_geyer"])
    return {"ratio_rows": ratio_rows, "unjoined": unjoined, "profile_rows": prof_rows, "summary": summary, "value_unit": cols.get("value_unit", "")}


def deep_dataset_check(category: str, dec, manifest, pools, log) -> dict:
    """The share of Geyer's deep dataset (Table S5) that Table S2 does not carry: S5 rows none of
    whose accessions appear in any S2 protein group, as a share of S5's summed signal. A completeness
    number for the plasma listing; reported, not used."""
    cc = dec["cross_checks"]
    src, fkey, pool = cc["deep_source"], cc["deep_file"], cc["deep_pool"]
    cols = dec[f"columns.{src}.{fkey}"]
    path, sha = verify_file(manifest, src, fkey)
    rows = read_table(path, cols["sheet"], int(cols["header_row"]),
                      {"identity": cols["identity_column"], "gene": cols["gene_column"], "abundance": cols["abundance_column"]}, log)
    log10 = cols.get("abundance_scale", "").strip().lower() == "log10"
    universe: set[str] = set()
    for jr in common.read_tsv(common.category_data_dir(category) / f"{pool}_dataset_rows.tsv"):
        for t in jr["identity_cell"].split(";"):
            kind, acc = classify_token(t)
            if acc:
                universe.add(acc)
    absent_rows, total, absent_total, n_val = [], 0.0, 0.0, 0
    for r in rows:
        x = cell_float(r.get("abundance"))
        if x is None:
            continue
        v = 10 ** x if log10 else x
        total += v; n_val += 1
        toks = [classify_token(t)[1] for t in str(r["identity"] or "").split(";")]
        if not any(t in universe for t in toks if t):
            absent_total += v
            absent_rows.append([r["row_index"], str(r["identity"] or ""), str(r.get("gene") or ""), fnum(v)])
    for row in absent_rows:
        row.append(fnum(float(row[3]) / total) if total else "")
    absent_rows.sort(key=lambda r: -float(r[3]))
    summary = {"deep_file": f"{path.name} sha256 {sha}", "deep_rows_with_value": str(n_val), "deep_rows_absent_from_primary": str(len(absent_rows)),
               "deep_signal_share_absent_from_primary": fnum(absent_total / total if total else 0.0)}
    log.info("deep dataset check: %d of %d rows absent from %s's primary table, %.4f of the deep signal",
             len(absent_rows), n_val, pool, absent_total / total if total else 0)
    return {"rows": absent_rows, "summary": summary}


def overlap_shares(pools: dict[str, dict], membership) -> tuple[list[list], dict]:
    names = sorted(pools)
    rows, sums = [], {p: 0.0 for p in names}
    for acc in sorted(membership):
        ps = [p for p in names if acc in pools[p]["entries"]]
        if len(ps) > 1:
            rows.append([acc, membership[acc]["gene"]] + [fnum(pools[p]["entries"][acc]["w"]) if p in ps else "" for p in names])
            for p in ps:
                sums[p] += pools[p]["entries"][acc]["w"]
    rows.sort(key=lambda r: -max(float(x) for x in r[2:] if x))
    return rows, sums


# ----------------------------------------------------------------------------- render

def per_entry_table(category: str, pools: dict[str, dict], membership, comp, header: list[str]) -> str:
    names = sorted(pools)
    hdr = ["accession", "gene", "tier", "mw_master"]
    for p in names:
        hdr += [f"match_rule_{p}", f"n_rows_{p}", f"rows_{p}", f"v_{p}", f"sd_{p}", f"n_{p}", f"published_share_{p}", f"shared_row_value_{p}",
                f"w_{p}", f"w_{p}_low", f"w_{p}_high"]
        for k in range(pools[p]["meta"]["n_reps"]):
            hdr.append(f"v_{p}_rep{k + 1}")
        if pools[p]["meta"]["n_runs"]:
            hdr.append(f"cv_percent_published_{p}")
    body = []
    for acc in sorted(membership):
        if membership[acc]["flag_open"]:
            continue
        row = [acc, membership[acc]["gene"], ";".join(membership[acc]["pools"]), f"{comp[acc]['mw']:.4f}" if acc in comp else ""]
        for p in names:
            e = pools[p]["entries"].get(acc)
            if e is None:
                row += [""] * 11 + [""] * pools[p]["meta"]["n_reps"] + ([""] if pools[p]["meta"]["n_runs"] else [])
                continue
            row += [e["match_rule"], len(e["rows"]), ";".join(str(i) for i in e["rows"]) or "none", fnum(e["v"]), fnum(e["sd"]), e["n"],
                    fnum(e["published_share"]), fnum(e["shared_whole"]), fnum(e["w"]), fnum(e["w_low"]), fnum(e["w_high"])]
            row += [fnum(x) for x in e["rep"]]
            if pools[p]["meta"]["n_runs"]:
                row.append("" if e["cv"] is None else fnum(e["cv"]))
        body.append(row)
    return tsv_text(header, hdr, body)


def build(category: str, log) -> tuple[dict[Path, str], dict]:
    files = common.category_files(category)
    psd = common.read_ini(files["decisions"])
    dec = common.read_ini(common.category_config_dir(category) / "mass_fraction_decisions.ini")
    if dec["scope"]["category"].strip() != category:
        raise SystemExit(f"[STOP] [scope] category in mass_fraction_decisions.ini is not {category} (D106)")
    manifest = common.read_ini(common.literature_manifest())
    sources = common.read_ini(common.LITERATURE_SOURCES_INI)
    membership = load_membership(files)
    comp = load_composition(files)
    for acc, m in membership.items():
        if not m["flag_open"] and acc not in comp:
            raise SystemExit(f"[STOP] {acc} is in accessions.ini but has no master row in the composition table; run blood-composition")
    pools = {p: pool_values(category, p, psd, dec, manifest, membership, log) for p in common.pool_names(category)}

    bounds = weighted_bounds(pools, files)
    excluded_rows, excluded_summary = [], {}
    for p in pools:
        r, s = excluded_shares(category, p, pools[p], membership, log)
        excluded_rows += r
        excluded_summary[p] = s
    ig_entries, ig_chains = immunoglobulin_bound(pools, membership, comp, dec, log)
    hortin = hortin_check(dec, manifest, pools, membership, comp, log)
    deep = deep_dataset_check(category, dec, manifest, pools, log)
    overlap_rows, overlap_sums = overlap_shares(pools, membership)

    acc_sha = common.sha256_file(files["accessions"])
    comp_sha = common.sha256_file(files["composition_tsv"])
    header = ["GENERATED by mass_fractions_pools.py. DO NOT EDIT BY HAND.",
              f"generated   = {common.iso_now()}", f"category    = {category}"]
    for p in sorted(pools):
        m = pools[p]["meta"]
        cite = sources[f"source.{m['source']}"].get("citation", "") if sources.has_section(f"source.{m['source']}") else ""
        header.append(f"pool {p}: {m['source']} {m['file_key']} {m['path']} sha256 {m['sha256']}; sheet {m['sheet']!r}, header row {m['header_row']}; {cite}")
    header += [f"accessions  = {files['accessions'].relative_to(common.REPO_ROOT).as_posix()} sha256 {acc_sha}",
               f"composition = {files['composition_tsv'].relative_to(common.REPO_ROOT).as_posix()} sha256 {comp_sha}",
               "rules       = B0 (hash), D105 (published mass share, no MW step), D108 (mean of replicates), D47 (summed rows), D48 (shared rows), D53a, D31 (within pool), B5 (low/high), A7, D113 (immunoglobulin bound)"]
    names = sorted(pools)
    out = files["mass_fractions_out"]
    outputs: dict[Path, str] = {}
    prov = header + [
        "v = the entry's published mass share (mass_share, D105), summed over its rows (D47) plus its split of shared rows (D48); as a fraction of the listing (not renormalised).",
        "sd / n: share_columns -> sd over the replicates of v (the replicate sums, D53a), n = replicates; replicate_columns (log10 runs) -> sd of 10^run in linear space, n = runs present.",
        "w_<pool> = v / sum over the pool's entries (D31 within pool: the listing renormalised after the contaminant rule removed its rows); _low/_high = (v -/+ sd) over the same denominator, _high plus the whole of every shared row (B5, D48).",
        "v_<pool>_rep<k> = the entry's share in replicate k (the per-donor sensitivity's input); cv_percent_published = the dataset's own CV column, carried for the first row of the entry.",
        "match_rule: single | summed:n | shared_split | none (no dataset row, or every row zero-weighted; w = 0).",
    ]
    outputs[files["mass_fractions_tsv"]] = per_entry_table(category, pools, membership, comp, prov)
    for p in names:
        ent = pools[p]["entries"]
        ranked = sorted(ent.items(), key=lambda kv: -kv[1]["w"])
        cum, rows = 0.0, []
        for acc, e in ranked:
            cum += e["w"]
            rows.append([acc, membership[acc]["gene"], fnum(e["w"]), fnum(e["w_low"]), fnum(e["w_high"]), fnum(cum), e["match_rule"]])
        outputs[out / f"{p}_entries_ranked.tsv"] = tsv_text(header, ["accession", "gene", "w", "w_low", "w_high", "cumulative_w", "match_rule"], rows)
    outputs[out / "shared_rows.tsv"] = tsv_text(header, ["pool", "row_index", "entries", "row_value", "split"], sum((pools[p]["shared_table"] for p in names), []))
    outputs[out / "weighted_bounds.tsv"] = tsv_text(header + ["D55 pattern: the largest single-entry isoform / processing / PTM delta times its weight; glycosylated_entries_sum_w = the weight of entries with any glycosylation site"],
                                                    ["pool", "kind", "amino_acid", "max_weighted_delta", "where"], bounds)
    outputs[out / "excluded_rows_mass_share.tsv"] = tsv_text(header + ["removed_by_rule = D98 contaminants (removed on purpose); completeness_gap = rows no entry answers to and entries build_protein_set excluded (R5), A7"],
                                                             ["pool", "row_index", "identity_cell", "gene_cell", "entry", "outcome", "kind", "share_of_listing"], excluded_rows)
    outputs[out / "immunoglobulin_variable_domain_bound.tsv"] = tsv_text(header + ["D113: per immunoglobulin entry of the pool; bound = w x L_V / L_C for a constant-region entry"],
                                                                          ["pool", "accession", "gene", "chain", "segment_class", "w", "L_C_master_residues", "L_V_mean_residues", "variable_mass_bound", "how"], ig_entries)
    outputs[out / "immunoglobulin_variable_domain_bound_per_chain.tsv"] = tsv_text(header,
                                                                                    ["pool", "chain", "n_constant_entries", "constant_entries_w", "L_V_mean_residues", "n_V_entries", "variable_entries_w_measured", "variable_mass_bound", "loose_bound_whole_share"], ig_chains)
    outputs[out / "ratio_check_hortin_2008.tsv"] = tsv_text(header + [f"Hortin 2008 ({hortin['summary']['hortin_file']}), value = {hortin['value_unit']}; share_hortin over the joined set; share_geyer_joined = the pool's v renormalised over the joined set; reported, not used"],
                                                             ["accession", "gene", "hortin_value", "share_hortin", "w_pool", "share_geyer_joined", "ratio_geyer_over_hortin"], hortin["ratio_rows"])
    outputs[out / "hortin_2008_rows_not_joined.tsv"] = tsv_text(header, ["row_index", "polypeptide", "identity_cell", "gene", "value", "why"], hortin["unjoined"])
    outputs[out / "profile_plasma_geyer_vs_hortin.tsv"] = tsv_text(header + ["amino acid mass fractions (free convention) of the plasma pool: under the pool's weights (all entries), under the same weights over the Hortin-joined set only, and under Hortin's mg/L over that set; reported, not used"],
                                                                    ["amino_acid", "profile_pool_all_entries", "profile_pool_joined_set", "profile_hortin_joined_set", "hortin_minus_pool_joined", "relative_difference"], hortin["profile_rows"])
    outputs[out / "completeness_check_geyer_2016_deep_dataset.tsv"] = tsv_text(header + ["rows of the deep dataset (Table S5) none of whose accessions appear in a Table S2 protein group, with their share of the deep dataset's summed signal; reported, not used"],
                                                                                ["row_index", "identity_cell", "gene", "value", "share_of_deep_signal"], deep["rows"])
    outputs[out / "pool_overlap_mass_share.tsv"] = tsv_text(header + ["entries present in more than one pool with their weight in each: plasma carry-over into the erythrocyte preparation and erythrocyte protein in plasma, at the shares each dataset measured"],
                                                            ["accession", "gene"] + [f"w_{p}" for p in names], overlap_rows)

    summary = common.new_ini()
    summary["inputs"] = {"accessions_sha256": acc_sha, "composition_sha256": comp_sha}
    for p in names:
        m, ent = pools[p]["meta"], pools[p]["entries"]
        ranked = sorted(ent.values(), key=lambda e: -e["w"])
        top = max(ent.items(), key=lambda kv: kv[1]["w"])
        summary[f"pool.{p}"] = {
            "dataset": f"{m['source']} {m['file_key']} {m['path']} sha256 {m['sha256']}",
            "entries": str(len(ent)), "entries_with_rows": str(sum(1 for e in ent.values() if e["rows"])),
            "entries_with_weight_above_zero": str(sum(1 for e in ent.values() if e["w"] > 0)),
            "entries_summed": str(sum(1 for e in ent.values() if e["match_rule"].startswith("summed"))),
            "entries_touched_by_shared_rows": str(sum(1 for e in ent.values() if e["shared_rows"])),
            "shared_rows": str(len(pools[p]["shared_table"])),
            "listing_total": fnum(m["listing_total"]), "denominator": fnum(m["denominator"]),
            "denominator_share_of_listing": fnum(m["denominator"] / m["listing_total"] if m["listing_total"] else 0),
            "sum_w": fnum(sum(e["w"] for e in ent.values())),
            "share_of_ten_largest": fnum(sum(e["w"] for e in ranked[:10])),
            "largest_entry": f"{top[0]} {membership[top[0]]['gene']} {fnum(top[1]['w'])}",
            "replicates": str(m["n_reps"]) if m["n_reps"] else f"{m['n_runs']} runs (technical)",
            "removed_by_rule_share_of_listing": fnum(excluded_summary[p]["removed_by_rule"]),
            "completeness_gap_share_of_listing": fnum(excluded_summary[p]["completeness_gap"]),
            "entries_not_in_accessions_ini": str(excluded_summary[p]["entries_not_in_accessions_ini"]),
            "overlap_entries_sum_w": fnum(overlap_sums[p]),
        }
    summary["immunoglobulin_bound"] = {f"{r[0]}.{r[1]}": f"bound {r[7]} (loose {r[8]}; constant entries w {r[3]}; measured variable entries w {r[6]})" for r in ig_chains}
    summary["cross_check.hortin_2008"] = hortin["summary"]
    summary["cross_check.geyer_2016_deep_dataset"] = deep["summary"]
    summary["overlap"] = {"entries_in_more_than_one_pool": str(len(overlap_rows)), **{f"sum_w_{p}": fnum(overlap_sums[p]) for p in names}}
    outputs[out / "mass_fractions_summary.ini"] = common.render_ini(summary, header)
    return outputs, {"pools": pools}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--category", required=True, help="the category whose pools are published mass shares (blood)")
    ap.add_argument("--check", action="store_true", help="rebuild in memory and compare with the files on disk")
    args = ap.parse_args(argv)
    log = common.make_logger(f"{STAGE}_{args.category}" + ("_check" if args.check else ""))
    outputs, _ = build(args.category, log)
    if args.check:
        stale = [p for p, text in outputs.items()
                 if not p.exists() or strip_generated_line(p.read_text(encoding="utf-8")) != strip_generated_line(text)]
        if stale:
            log.error("check: %d file(s) differ from a fresh build: %s", len(stale), [p.name for p in stale])
            return 1
        log.info("check: mass-fraction outputs and config are current")
        return 0
    for p, text in outputs.items():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8", newline="\n")
        log.info("wrote %s", p.relative_to(common.REPO_ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
