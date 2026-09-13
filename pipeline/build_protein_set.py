"""Generate config/accessions.ini and config/segments.ini from fetched data, by code.

Runs offline. Inputs are all on disk:
  config/protein_set_decisions.ini     (hand-written: seeds, rules, closed flags)
  data/gene-ontology/resolved_terms.ini, tier<N>_subtree.tsv
  data/tier<N>_pool.tsv
  data/pool_queries.ini
  data/uniprot_sequences.ini

Rules applied (each names the field it reads; see docs/conventions.md)
  R1  Sequence   : the canonical sequence of the entry. Always. No isoform choice.
  R2  Processing : UniProt feature table, type "Chain".
        R2a exactly one Chain feature -> its range is the master molecule; residues
            outside it are in_master_molecule = false.
        R2b no Chain feature          -> whole sequence is the master molecule; logged.
        R2d several Chain features    -> union of their ranges (D24); listed in
            outputs/multi_chain_entries.tsv.
        R2c any Chain with a non-exact position -> FLAG.
  R3  Evidence   : Gene Ontology evidence codes are recorded per annotation and
        summarised; they never filter.
  R4  Tier       : membership of data/tier<N>_pool.tsv. An entry in several pools
        carries all its tiers.

A flag is closed only by a [flag.<accession>] section in the decisions file:
    rule = chain_position_uncertain
    decision = D<n>
    master_start = <int>
    master_end   = <int>
Open flags leave the accession in accessions.ini with flag_open = true and OUT of
segments.ini, and the build exits non-zero.

  --check   rebuild in memory and compare with the files on disk; exit non-zero on
            any difference (used by the tests: generated config must be current).
"""

from __future__ import annotations

import argparse
import configparser
import sys
from collections import Counter

from pipeline import common
from pipeline.flags import append_flag, reset_flags, closure

STAGE = "build_protein_set"


# ----------------------------------------------------------------------------- helpers
def _features(sec, ftype: str) -> list[dict[str, str]]:
    out = []
    i = 0
    while f"feature.{i}.type" in sec:
        if sec[f"feature.{i}.type"] == ftype:
            out.append({k: sec.get(f"feature.{i}.{k}", "") for k in ("type", "start", "end", "description", "id")})
        i += 1
    return out


def _all_features(sec) -> list[dict[str, str]]:
    out = []
    i = 0
    while f"feature.{i}.type" in sec:
        out.append({k: sec.get(f"feature.{i}.{k}", "") for k in ("type", "start", "end", "description", "id")})
        i += 1
    return out


def _annotations(sec) -> list[dict[str, str]]:
    out = []
    i = 0
    while f"annotation.{i}.term_id" in sec:
        out.append({k: sec.get(f"annotation.{i}.{k}", "") for k in ("term_id", "term", "evidence")})
        i += 1
    return out


def _is_exact(pos: str) -> bool:
    return pos.isdigit()


def _union(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent 1-indexed inclusive ranges into sorted disjoint runs."""
    out: list[tuple[int, int]] = []
    for a, b in sorted(ranges):
        if out and a <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def chain_rule(acc: str, sec, decisions):
    """Apply R2. Returns (rule_applied, master_ranges, note, flag_or_None, chains).

    R2a one exact Chain           -> its range
    R2b no Chain                  -> whole sequence
    R2d several exact Chains      -> union of their ranges (D24); listed in
                                     outputs/multi_chain_entries.tsv, not flagged
    R2c any Chain with a non-exact position -> FLAG chain_position_uncertain; closable
         with master_start/master_end in the decisions file (rule 'R2c-closed').
    The flag dict is returned whenever the condition exists, closed or not, so it is
    always recorded; the caller checks the closure to decide whether it is open.
    """
    length = int(sec["length"])
    chains = _features(sec, "Chain")
    if not chains:
        return ("R2b", [(1, length)],
                "no Chain feature in UniProt feature table; whole sequence is the master molecule", None, chains)
    inexact = [c for c in chains if not (_is_exact(c["start"]) and _is_exact(c["end"]))]
    if inexact:
        desc = "; ".join(f"{c['id']} '{c['description']}' {c['start']}-{c['end']}" for c in inexact)
        flag = {"rule": "chain_position_uncertain", "detail": f"Chain feature(s) with non-exact position: {desc}"}
        cl = closure(decisions, acc)
        if cl is not None and cl.get("rule") == flag["rule"]:
            s, e = int(cl["master_start"]), int(cl["master_end"])
            return ("R2c-closed", [(s, e)], f"flag {flag['rule']} closed by {cl['decision']}: master = {s}-{e}", flag, chains)
        return ("R2c", [], "", flag, chains)
    ranges = [(int(c["start"]), int(c["end"])) for c in chains]
    if len(chains) == 1:
        c = chains[0]
        return ("R2a", ranges, f"Chain {c['id']} '{c['description']}' {c['start']}-{c['end']}", None, chains)
    merged = _union(ranges)
    desc = "; ".join(f"{c['id']} '{c['description']}' {c['start']}-{c['end']}" for c in chains)
    return ("R2d", merged,
            f"union of {len(chains)} Chain features = {_fmt_ranges(merged)} (D24); chains: {desc}", None, chains)


def _fmt_ranges(ranges: list[tuple[int, int]]) -> str:
    return ",".join(f"{a}-{b}" for a, b in ranges)


def segments_for(acc: str, sec, ranges: list[tuple[int, int]], rule: str, note: str) -> list[tuple[str, dict[str, str]]]:
    """Segment sections for one accession: each master run plus every residue run outside them.
    Section names: .master (first run), .master_2.. for further runs; .n_terminal_removed,
    .internal_removed[_k], .c_terminal_removed for uncovered runs."""
    length = int(sec["length"])
    feats = _all_features(sec)

    def covering(a: int, b: int) -> str:
        hits = []
        for f in feats:
            if f["type"] in ("Chain", "Alternative sequence"):
                continue
            if _is_exact(f["start"]) and _is_exact(f["end"]):
                fs, fe = int(f["start"]), int(f["end"])
                if fs <= b and fe >= a:
                    hits.append(f"{f['type']} {fs}-{fe}")
        return "; ".join(hits) if hits else "no captured feature covers this span"

    out = []
    pos = 1
    n_master = 0
    n_internal = 0
    for a, b in ranges:
        if a > pos:
            if pos == 1:
                name = f"{acc}.n_terminal_removed"
            else:
                n_internal += 1
                name = f"{acc}.internal_removed" + ("" if n_internal == 1 else f"_{n_internal}")
            out.append((name, {"start": str(pos), "end": str(a - 1), "in_master_molecule": "false",
                               "note": f"outside every Chain range; UniProt features here: {covering(pos, a - 1)}"}))
        n_master += 1
        name = f"{acc}.master" + ("" if n_master == 1 else f"_{n_master}")
        out.append((name, {"start": str(a), "end": str(b), "in_master_molecule": "true", "rule": rule, "note": note}))
        pos = b + 1
    if pos <= length:
        out.append((f"{acc}.c_terminal_removed", {"start": str(pos), "end": str(length), "in_master_molecule": "false",
                    "note": f"outside every Chain range; UniProt features here: {covering(pos, length)}"}))
    return out


# ----------------------------------------------------------------------------- build
def build(decisions, resolved, sequences, pools: dict[str, list[dict[str, str]]],
          subtrees: dict[str, set[str]], queries, log, flags_path):
    """Returns (accessions_cp, segments_cp, header_lines, open_flags, evidence_rows, acc_evidence_rows)."""
    tiers_of: dict[str, list[str]] = {}
    for tier, rows in pools.items():
        for r in rows:
            tiers_of.setdefault(r["accession"], []).append(tier)

    ver_header = _ini_header(common.SEQUENCES_INI)
    acc_cp = common.new_ini()
    seg_cp = common.new_ini()
    open_flags: list[str] = []
    evidence_counter: dict[str, Counter] = {t: Counter() for t in pools}
    acc_evidence_rows: list[list] = []
    multi_chain_rows: list[list] = []

    excluded_rows: list[list] = []
    for acc in sorted(tiers_of):
        if not sequences.has_section(acc):
            raise RuntimeError(f"{acc} is in a pool but not in {common.SEQUENCES_INI}; rerun fetch_sequences")
        sec = sequences[acc]
        if sec.get("md5_match") != "true":
            raise RuntimeError(f"{acc}: md5_match is not true in data/uniprot_sequences.ini; build refuses to proceed")
        tiers = sorted(tiers_of[acc])

        # R5 (D59): the alphabet is the twenty coded amino acids with fetched masses; any other
        # letter in the canonical sequence excludes the entry from the set, listed with its UniProt mass
        letters = "".join(sorted(set(sec["sequence"]) - set(common.AMINO_ACIDS)))
        if letters:
            excluded_rows.append([acc, sec.get("gene", ""), ";".join(tiers), letters, sec.get("length", ""),
                                  _uniprot_mol_weight(acc), sec.get("protein_name", "")])
            log.info("%s %s: excluded under R5 (letters %s)", acc, sec.get("gene", ""), letters)
            continue

        # R3: annotations inside each tier's subtree, with evidence codes
        anns = _annotations(sec)
        subtree_terms_by_tier = {}
        for t in tiers:
            if t not in subtrees:                      # measured tier (D56): no ontology category, no R3
                subtree_terms_by_tier[t] = []
                evidence_counter[t]["(measured remainder; no ontology annotation applies)"] += 1
                acc_evidence_rows.append([acc, t, "", ""])
                continue
            hits = [a for a in anns if a["term_id"] in subtrees[t]]
            subtree_terms_by_tier[t] = hits
            codes = sorted({a["evidence"].split(":")[0] for a in hits if a["evidence"]})
            for code in codes:
                evidence_counter[t][code] += 1
            if not hits:
                evidence_counter[t]["(no subtree annotation in entry cross-references)"] += 1
                subject = f"{acc}.tier{t}"
                append_flag(flags_path, stage=STAGE, subject=subject,
                            rule="pool_member_without_subtree_annotation",
                            detail=(f"returned by the UniProt query for tier {t} but no Gene Ontology "
                                    f"cross-reference in the entry JSON lies in the tier's subtree"))
                cl = closure(decisions, subject)
                if cl is None or cl.get("rule") != "pool_member_without_subtree_annotation":
                    open_flags.append(subject)
            acc_evidence_rows.append([acc, t, ";".join(f"{a['term_id']}({a['evidence']})" for a in hits),
                                      ";".join(codes)])

        # R2: processing
        rule, ranges, note, flag, chains = chain_rule(acc, sec, decisions)
        if flag is not None:
            append_flag(flags_path, stage=STAGE, subject=acc, rule=flag["rule"], detail=flag["detail"])
        if rule == "R2c":
            open_flags.append(acc)
            log.warning("%s: FLAG %s -- %s", acc, flag["rule"], flag["detail"])
        else:
            if rule in ("R2b", "R2c-closed", "R2d"):
                log.info("%s: %s", acc, note)
            if rule == "R2d":
                multi_chain_rows.append([acc, sec.get("gene", ""), ";".join(tiers), len(chains),
                                         " | ".join(f"{c['id']} '{c['description']}' {c['start']}-{c['end']}" for c in chains),
                                         _fmt_ranges(ranges), sec["length"]])
            for name, body in segments_for(acc, sec, ranges, rule, note):
                seg_cp[name] = body

        # accessions.ini section
        terms_txt = " ; ".join(
            f"tier{t}: " + (", ".join(f"{a['term_id']} [{a['evidence']}]" for a in subtree_terms_by_tier[t]) if t in subtrees else "(measured remainder)")
            for t in tiers)
        acc_cp[acc] = {
            "accession": acc,
            "entry_name": sec.get("entry_name", ""),
            "gene": sec.get("gene", ""),
            "protein_name": sec.get("protein_name", ""),
            "tier": ";".join(tiers),
            "isoform": "canonical",
            "isoform_rule": "R1",
            "canonical_isoform_id": sec.get("canonical_isoform_id", ""),
            "n_isoforms_in_entry": sec.get("n_isoforms", "0"),
            "length": sec["length"],
            "seq_version": sec.get("seq_version", ""),
            "entry_version": sec.get("entry_version", ""),
            "md5": sec.get("md5_uniprot", ""),
            "reviewed": "true" if "reviewed" in sec.get("entry_type", "").lower() else "false",
            "organism_id": sec.get("organism_id", ""),
            "subtree_annotations": terms_txt,
            "master_rule": rule if rule != "R2c" else "FLAGGED",
            "flag_open": "true" if rule == "R2c" else "false",
            "pool_membership": "computed: " + "; ".join(
                (f"tier {t} = subtree of {resolved[f'tier.{t}']['term_id']} ({resolved[f'tier.{t}']['term_name']})"
                 if t in subtrees else f"tier {t} = measured remainder of {resolved[f'tier.{t}']['dataset_source']} ({resolved[f'tier.{t}']['definition_decision']})")
                for t in tiers),
            "uniprot_release": ver_header.get("uniprot_release", ""),
            "fetched": ver_header.get("fetched", ""),
        }

    header = [
        "GENERATED by build_protein_set.py. DO NOT EDIT BY HAND -- edit config/protein_set_decisions.ini and rebuild.",
        f"ontology_file         = {resolved['ontology']['file']} (data_version {resolved['ontology']['data_version']}, sha256 {resolved['ontology']['sha256']})",
        f"uniprot_release       = {ver_header.get('uniprot_release', '')}",
        f"uniprot_fetched       = {ver_header.get('fetched', '')}",
        f"pool_queries          = data/pool_queries.ini (run {queries['run']['date']})",
    ]
    for t in sorted(pools):
        r = resolved[f"tier.{t}"]
        if t in subtrees:
            header.append(f"tier {t}: seed_word={r['seed_word']} ({r['seed_decision']}) lookup_name={r['lookup_name']} "
                          f"({r['lookup_decision']}) -> {r['term_id']} {r['term_name']} ; pool size {len(pools[t])}")
        else:
            header.append(f"tier {t}: seed_word={r['seed_word']} ({r['seed_decision']}) definition={r['definition']} "
                          f"({r['definition_decision']}) from {r['dataset_source']} {r['dataset_file']} minus {r['contaminant_source']} ; pool size {len(pools[t])}")
    evidence_rows = [[t, code, n] for t in sorted(evidence_counter)
                     for code, n in sorted(evidence_counter[t].items())]
    header.append(f"R5 (D59): {len(excluded_rows)} entries excluded for letters outside the twenty coded amino acids; see outputs/excluded_non_standard_alphabet.tsv")
    return acc_cp, seg_cp, header, open_flags, evidence_rows, acc_evidence_rows, multi_chain_rows, excluded_rows


def _uniprot_mol_weight(acc: str) -> str:
    """UniProt's own molecular weight (sequence.molWeight, Da) from the raw entry JSON, if on disk;
    used only to state the mass an R5-excluded entry would have carried."""
    import json
    p = common.DATA_DIR / "uniprot_raw" / f"{acc}.json"
    if not p.exists():
        return ""
    try:
        return str(json.loads(p.read_text(encoding="utf-8")).get("sequence", {}).get("molWeight", ""))
    except Exception:  # noqa: BLE001 - disclosure only
        return ""


def _ini_header(path) -> dict[str, str]:
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


def _render(cp: configparser.ConfigParser, header: list[str]) -> str:
    return common.render_ini(cp, header)


def load_inputs():
    decisions = common.read_ini(common.DECISIONS_INI)
    resolved = common.read_ini(common.RESOLVED_TERMS_INI)
    sequences = common.read_ini(common.SEQUENCES_INI)
    queries = common.read_ini(common.POOL_QUERIES_INI)
    pools, subtrees = {}, {}
    for section in resolved.sections():
        if section.startswith("tier."):
            t = section.split(".", 1)[1]
            pools[t] = common.read_tsv(common.DATA_DIR / f"tier{t}_pool.tsv")
            if resolved[section].get("term_id"):
                subtrees[t] = {r["term_id"] for r in common.read_tsv(common.ONTOLOGY_DIR / f"tier{t}_subtree.tsv")}
    return decisions, resolved, sequences, pools, subtrees, queries


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare a fresh build with config on disk; write nothing")
    args = ap.parse_args(argv)
    log = common.make_logger(STAGE + ("_check" if args.check else ""))

    decisions, resolved, sequences, pools, subtrees, queries = load_inputs()
    flags_path = common.FLAGS_TSV if not args.check else common.OUTPUTS_DIR / "_flags_check.tmp"
    if args.check and flags_path.exists():
        flags_path.unlink()
    if not args.check:
        reset_flags(common.FLAGS_TSV, STAGE)

    acc_cp, seg_cp, header, open_flags, evidence_rows, acc_evidence_rows, multi_chain_rows, excluded_rows = build(
        decisions, resolved, sequences, pools, subtrees, queries, log, flags_path)
    seg_header = [
        "GENERATED by build_protein_set.py. DO NOT EDIT BY HAND.",
        "Segment flags per accession. in_master_molecule = true marks the mature-chain range",
        "(rule R2, UniProt 'Chain' feature). Positions are 1-indexed, inclusive.",
    ]
    acc_text = _render(acc_cp, header)
    seg_text = _render(seg_cp, seg_header)

    if args.check:
        if flags_path.exists():
            flags_path.unlink()
        ok = True
        for path, text in ((common.ACCESSIONS_INI, acc_text), (common.SEGMENTS_INI, seg_text)):
            if not path.exists() or path.read_text(encoding="utf-8") != text:
                ok = False
                log.error("STALE: %s differs from a fresh build", path)
        log.info("check: %s", "generated config is current" if ok else "generated config is STALE")
        return 0 if ok else 1

    common.ACCESSIONS_INI.write_text(acc_text, encoding="utf-8", newline="\n")
    common.SEGMENTS_INI.write_text(seg_text, encoding="utf-8", newline="\n")
    common.write_tsv(common.OUTPUTS_DIR / "evidence_summary.tsv", ["tier", "evidence_code", "n_accessions"], evidence_rows)
    common.write_tsv(common.OUTPUTS_DIR / "accession_evidence.tsv",
                     ["accession", "tier", "subtree_annotations", "evidence_codes"], acc_evidence_rows)
    common.write_tsv(common.OUTPUTS_DIR / "excluded_non_standard_alphabet.tsv",
                     ["accession", "gene", "tier", "non_standard_letters", "length", "uniprot_mol_weight_da", "protein_name"],
                     excluded_rows)
    common.write_tsv(common.OUTPUTS_DIR / "multi_chain_entries.tsv",
                     ["accession", "gene", "tier", "n_chain_features", "chain_features", "master_ranges_union", "length"],
                     multi_chain_rows)
    log.info("wrote %s (%d accessions) and %s (%d segments)",
             common.ACCESSIONS_INI, len(acc_cp.sections()), common.SEGMENTS_INI, len(seg_cp.sections()))
    if open_flags:
        log.error("%d open flag(s); close them in %s and rebuild: %s",
                  len(open_flags), common.DECISIONS_INI, ", ".join(open_flags))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
