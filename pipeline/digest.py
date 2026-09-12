#!/usr/bin/env python3
"""
digest.py — offline stage: in-silico tryptic digest of every canonical sequence. Pure Layer A.

For every accession in data/uniprot_sequences.ini, performs an in-silico
tryptic digest under the rules in config/mass_fraction_decisions.ini [digest]
and writes:

  outputs/digest/theoretical_peptides.tsv
      accession, gene, tier, length, mw_full, n_peptides, peptides_per_kDa
  outputs/digest/density_ranked.tsv
      the same rows ranked by peptides_per_kDa, with robust z
      (x − median) / (1.4826·MAD). No threshold is applied: the table
      discloses; the mass-fraction stage's checks read it.
  outputs/digest/shared_peptides.tsv
      peptide, n_accessions, accessions   (only peptides in ≥ 2 entries)
  outputs/digest/shared_pairs.tsv
      accession_a, gene_a, accession_b, gene_b, n_shared_peptides,
      frac_of_a, frac_of_b   (edge weights; a family joined by one peptide
      is visible here as n_shared_peptides = 1)
  outputs/digest/families.tsv
      family_id, n_members, accessions, genes, n_shared_peptides
      (gene symbols come from the composition table, i.e. accessions.ini)
      (connected components of the "shares ≥ 1 in-window peptide" graph)
  outputs/digest/digest_summary.ini
  logs/digest_<UTC>.log

Why this exists
  iBAQ = summed intensity / N_theoretical.  iBAQ×MW and TPA differ exactly by
  N_i versus MW_i, so peptides_per_kDa is the per-protein factor between the
  two conventions — computable from sequence alone (methods §"Quantity used").
  Shared peptides are where a razor-peptide assignment can move intensity
  between entries; the families table names those groups by computation.

Rules
  G1  Every [digest] parameter must be filled and `source` non-empty, else stop.
  G2  Cleavage after each residue in `cleave_after`; if cleave_before_proline
      is false, no cleavage when the next residue is P. If it is `both`, the
      count is computed under both settings; the primary is cleave-before-P
      and the relative difference per entry is written and summarised.
  G3  Peptides of min_length ≤ L ≤ max_length are counted, with up to
      `missed_cleavages` joined fragments. Distinct sequences per entry.
  G4  Digest runs on the FULL canonical sequence (what a search engine
      digests), not the master molecule. MW is taken from the composition
      table's full-product row (segment_set = metabolic).
  G5  No accession is named in code. Nothing is filtered.
"""

from __future__ import annotations

import configparser
import csv
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "mass_fraction_decisions.ini"
SEQ_FILE = ROOT / "data" / "uniprot_sequences.ini"
COMP_TSV = ROOT / "outputs" / "composition" / "amino_acid_composition_per_protein.tsv"
OUT_DIR = ROOT / "outputs" / "digest"
LOG_DIR = ROOT / "logs"
PLACEHOLDER = "___"
AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------- config

def load_rules(cp: configparser.ConfigParser) -> dict:
    d = cp["digest"]
    missing = [k for k in ("cleave_after", "cleave_before_proline", "min_length", "max_length", "missed_cleavages", "source")
               if d.get(k, "").strip() in ("", PLACEHOLDER)]
    if missing:
        raise SystemExit(f"[STOP] [digest] parameters unfilled or uncited: {', '.join(missing)} (G1)")
    cbp = d["cleave_before_proline"].strip().lower()
    if cbp not in ("true", "false", "both"):
        raise SystemExit("[STOP] [digest] cleave_before_proline must be true, false, or both (G1)")
    return {
        "enzyme": d.get("enzyme", "trypsin").strip(),
        "cleave_after": set(d["cleave_after"].strip().upper()),
        "cleave_before_proline": cbp,
        "min_length": int(d["min_length"]),
        "max_length": int(d["max_length"]),
        "missed_cleavages": int(d["missed_cleavages"]),
        "source": d["source"].strip(),
        "retrieved": d.get("retrieved", "").strip(),
    }


# ---------------------------------------------------------------- digest

def fragments(seq: str, cleave_after: set[str], cleave_before_proline: bool) -> list[str]:
    """Fully cleaved fragments (0 missed cleavages)."""
    out, start = [], 0
    for i, aa in enumerate(seq):
        if aa in cleave_after and i + 1 < len(seq):
            if not cleave_before_proline and seq[i + 1] == "P":
                continue
            out.append(seq[start:i + 1])
            start = i + 1
    out.append(seq[start:])
    return [f for f in out if f]


def theoretical_peptides(seq: str, r: dict, cleave_before_proline: bool | None = None) -> set[str]:
    cbp = r["cleave_before_proline"] if cleave_before_proline is None else cleave_before_proline
    if cbp == "both":
        raise ValueError("pass cleave_before_proline explicitly when rules say both")
    frags = fragments(seq, r["cleave_after"], cbp is True or cbp == "true")
    peps: set[str] = set()
    for i in range(len(frags)):
        joined = ""
        for j in range(i, min(i + r["missed_cleavages"] + 1, len(frags))):
            joined += frags[j]
            if r["min_length"] <= len(joined) <= r["max_length"]:
                peps.add(joined)
            if len(joined) > r["max_length"]:
                break
    return peps


# ---------------------------------------------------------------- inputs

def load_sequences() -> dict[str, str]:
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    with SEQ_FILE.open(encoding="utf-8") as fh:
        cp.read_file(fh)
    seqs: dict[str, str] = {}
    for sec in cp.sections():
        s = cp[sec]
        seq_key = next((k for k in ("sequence", "seq", "canonical_sequence") if k in s), None)
        if seq_key is None:
            raise SystemExit(f"[STOP] section [{sec}] in {SEQ_FILE.name} has no sequence key; keys present: {list(s.keys())}")
        acc = s.get("accession", sec).strip()
        seq = "".join(s[seq_key].split()).upper()
        bad = set(seq) - AA20
        if bad:
            raise SystemExit(f"[STOP] {acc}: non-standard letters {sorted(bad)}")
        seqs[acc] = seq
    return seqs


def load_composition() -> dict[str, dict]:
    """accession -> {tier, mw_full, gene} from the full-product row."""
    if not COMP_TSV.exists():
        raise SystemExit(f"[STOP] composition table not found: {COMP_TSV.relative_to(ROOT).as_posix()} — run `python run.py composition` first")
    out: dict[str, dict] = {}
    with COMP_TSV.open(encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"accession", "tier", "segment_set", "mw"}
        if not need <= set(rd.fieldnames or []):
            raise SystemExit(f"[STOP] composition table lacks columns {sorted(need - set(rd.fieldnames or []))}")
        for row in rd:
            if row["segment_set"] == "metabolic":
                out[row["accession"]] = {"tier": row["tier"], "mw_full": float(row["mw"]), "gene": row.get("gene", "")}
    return out


# ---------------------------------------------------------------- families

def connected_components(edges: dict[str, set[str]]) -> list[set[str]]:
    seen, comps = set(), []
    for node in edges:
        if node in seen:
            continue
        stack, comp = [node], set()
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.add(n)
            stack.extend(edges[n] - seen)
        comps.append(comp)
    return comps


# ---------------------------------------------------------------- main

def main(argv: list[str] | None = None) -> int:
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=(";",))
    cp.optionxform = str
    with CONFIG.open(encoding="utf-8") as fh:
        cp.read_file(fh)
    rules = load_rules(cp)
    seqs = load_sequences()
    comp = load_composition()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    log = [f"digest.py started {started}", f"rules: {rules}", f"sequences: {len(seqs)}", ""]
    flags: list[str] = []

    both = rules["cleave_before_proline"] == "both"
    primary_cbp = True if both else rules["cleave_before_proline"] == "true"
    per_acc: dict[str, set[str]] = {}
    rows = []
    for acc, seq in seqs.items():
        peps = theoretical_peptides(seq, rules, primary_cbp)
        per_acc[acc] = peps
        n_alt = len(theoretical_peptides(seq, rules, not primary_cbp)) if both else None
        c = comp.get(acc)
        if c is None:
            flags.append(f"{acc}: no full-product row in composition table (G4)")
            tier, mw, gene = "", float("nan"), ""
        else:
            tier, mw, gene = c["tier"], c["mw_full"], c["gene"]
        density = len(peps) / (mw / 1000.0) if mw == mw and mw > 0 else float("nan")
        row = {"accession": acc, "gene": gene, "tier": tier, "length": len(seq), "mw_full": mw,
               "n_peptides": len(peps), "peptides_per_kDa": density}
        if both:
            row["n_peptides_trypsin_P"] = len(peps)          # cleaves before proline
            row["n_peptides_trypsin"] = n_alt                 # does not
            row["proline_rule_rel_diff"] = (len(peps) - n_alt) / len(peps) if peps else float("nan")
        rows.append(row)

    # theoretical_peptides.tsv
    hdr = ["accession", "gene", "tier", "length", "mw_full", "n_peptides", "peptides_per_kDa"]
    if both:
        hdr += ["n_peptides_trypsin_P", "n_peptides_trypsin", "proline_rule_rel_diff"]
    write_tsv(OUT_DIR / "theoretical_peptides.tsv", hdr, rows, rules)

    # density_ranked.tsv
    dens = [r["peptides_per_kDa"] for r in rows if r["peptides_per_kDa"] == r["peptides_per_kDa"]]
    med = statistics.median(dens) if dens else float("nan")
    mad = statistics.median(abs(x - med) for x in dens) if dens else float("nan")
    scale = 1.4826 * mad if mad and mad > 0 else float("nan")
    ranked = sorted(rows, key=lambda r: (r["peptides_per_kDa"] != r["peptides_per_kDa"], r["peptides_per_kDa"]))
    for i, r in enumerate(ranked, 1):
        r["rank"] = i
        r["robust_z"] = (r["peptides_per_kDa"] - med) / scale if scale == scale and r["peptides_per_kDa"] == r["peptides_per_kDa"] else float("nan")
    write_tsv(OUT_DIR / "density_ranked.tsv", ["rank"] + hdr + ["robust_z"], ranked, rules)

    # shared peptides
    pep_to_acc: dict[str, set[str]] = defaultdict(set)
    for acc, peps in per_acc.items():
        for p in peps:
            pep_to_acc[p].add(acc)
    shared = {p: a for p, a in pep_to_acc.items() if len(a) >= 2}
    srows = [{"peptide": p, "n_accessions": len(a), "accessions": ";".join(sorted(a))}
             for p, a in sorted(shared.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    write_tsv(OUT_DIR / "shared_peptides.tsv", ["peptide", "n_accessions", "accessions"], srows, rules)

    # shared_pairs.tsv — edge weights: how many in-window peptides each pair shares
    pair_counts: dict[tuple[str, str], int] = defaultdict(int)
    for p, accs in shared.items():
        al = sorted(accs)
        for i in range(len(al)):
            for j in range(i + 1, len(al)):
                pair_counts[(al[i], al[j])] += 1
    prows = [{"accession_a": a, "gene_a": comp.get(a, {}).get("gene", ""), "accession_b": b, "gene_b": comp.get(b, {}).get("gene", ""),
              "n_shared_peptides": n, "frac_of_a": n / len(per_acc[a]) if per_acc[a] else float("nan"),
              "frac_of_b": n / len(per_acc[b]) if per_acc[b] else float("nan")}
             for (a, b), n in sorted(pair_counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    write_tsv(OUT_DIR / "shared_pairs.tsv", ["accession_a", "gene_a", "accession_b", "gene_b", "n_shared_peptides", "frac_of_a", "frac_of_b"], prows, rules)

    # families
    edges: dict[str, set[str]] = {acc: set() for acc in per_acc}
    for p, accs in shared.items():
        for a in accs:
            edges[a] |= accs - {a}
    comps = [c for c in connected_components(edges) if len(c) >= 2]
    comps.sort(key=lambda c: (-len(c), sorted(c)[0]))
    frows = []
    for i, c in enumerate(comps, 1):
        n_shared = sum(1 for p, a in shared.items() if a & c)
        members = sorted(c)
        frows.append({"family_id": f"F{i:03d}", "n_members": len(c), "accessions": ";".join(members),
                      "genes": ";".join(comp.get(a, {}).get("gene", "") or "?" for a in members), "n_shared_peptides": n_shared})
    write_tsv(OUT_DIR / "families.tsv", ["family_id", "n_members", "accessions", "genes", "n_shared_peptides"], frows, rules)

    # summary
    summ = configparser.ConfigParser(interpolation=None)
    summ.optionxform = str
    summ["digest"] = {k: str(v) if not isinstance(v, set) else "".join(sorted(v)) for k, v in rules.items()}
    summ["result"] = {
        "generated": utc_now(),
        "n_entries": str(len(rows)),
        "n_peptides_total_distinct": str(len(pep_to_acc)),
        "n_shared_peptides": str(len(shared)),
        "n_families": str(len(comps)),
        "n_sharing_pairs": str(len(pair_counts)),
        "largest_family": str(max((len(c) for c in comps), default=0)),
        "density_median_peptides_per_kDa": f"{med:.6g}",
        "density_mad": f"{mad:.6g}",
        "density_min": f"{min(dens):.6g}" if dens else "nan",
        "density_max": f"{max(dens):.6g}" if dens else "nan",
        "flags": str(len(flags)),
    }
    if both:
        diffs = [r["proline_rule_rel_diff"] for r in rows if r["proline_rule_rel_diff"] == r["proline_rule_rel_diff"]]
        summ["result"]["proline_rule"] = "both computed; primary = Trypsin/P (cleaves before P); n_peptides columns for both"
        summ["result"]["proline_rule_rel_diff_median"] = f"{statistics.median(diffs):.6g}" if diffs else "nan"
        summ["result"]["proline_rule_rel_diff_max"] = f"{max(diffs):.6g}" if diffs else "nan"
    with (OUT_DIR / "digest_summary.ini").open("w", encoding="utf-8") as fh:
        fh.write("# Generated by pipeline/digest.py — DO NOT EDIT BY HAND\n")
        summ.write(fh)

    log += [f"entries {len(rows)}; distinct peptides {len(pep_to_acc)}; shared {len(shared)}; families {len(comps)}",
            f"density median {med:.6g} MAD {mad:.6g}", "", f"flags: {len(flags)}"] + [f"  - {f}" for f in flags] + [f"finished {utc_now()}"]
    (LOG_DIR / f"digest_{started.replace(':', '')}.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    return 1 if flags else 0


def write_tsv(path: Path, header: list[str], rows: list[dict], rules: dict) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# Generated by pipeline/digest.py — DO NOT EDIT BY HAND\n")
        fh.write(f"# generated = {utc_now()}\n")
        fh.write(f"# digest = {rules['enzyme']}; cleave_after={''.join(sorted(rules['cleave_after']))}; "
                 f"cleave_before_proline={rules['cleave_before_proline']}; length {rules['min_length']}-{rules['max_length']}; "
                 f"missed_cleavages={rules['missed_cleavages']}; source={rules['source']}\n")
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    sys.exit(main())