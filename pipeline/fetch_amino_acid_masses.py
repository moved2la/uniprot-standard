"""Fetch the free amino acid masses (and water, hydrogen), by code.

The composition engine multiplies residue counts by masses. Those masses are the
only external constants in Layer A, and none of them is typed. This stage reads
config/composition_decisions.ini (two source URLs and one rule; D26) and:

  1. downloads the IUPAC-IUBMB Table 1 page that defines the one-letter code,
     saves and hashes it, and parses every "<trivial name> <three-letter> <one-letter>"
     row from its text. The run stops unless exactly the twenty standard letters
     were found, each exactly once. Result: data/iupac/amino_acid_symbols.ini.
  2. for each letter, queries PubChem PUG REST by compound name -- the chiral
     prefix from the decisions file joined to the trivial name; if PubChem has no
     compound under that name, the bare trivial name -- for CID, MolecularFormula,
     MolecularWeight, IUPACName. Exactly one compound must come back; zero (after
     the fallback) or several stops the run. Water and hydrogen are fetched the
     same way. Raw responses are saved unchanged.
  3. writes data/pubchem/amino_acid_masses.ini with, per letter, the served values,
     which query name resolved, and residue_mass = free_mass - water_mass.

PubChem serves MolecularWeight to two decimal places; the precision seen is
recorded in the table header and carried into the methods section.

Inputs   config/composition_decisions.ini
Outputs  data/iupac/tab1.html, data/iupac/source.ini, data/iupac/amino_acid_symbols.ini
         data/pubchem/raw/<letter or name>.json, data/pubchem/amino_acid_masses.ini
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import time

import requests

from pipeline import common

STAGE = "fetch_amino_acid_masses"
MIN_INTERVAL_S = 0.25   # PubChem asks for no more than 5 requests per second
RETRIES = 5
USER_AGENT = "uniprot-standard/composition (contact: repository README)"

_TAG = re.compile(r"<[^>]+>")
_WS = re.compile(r"\s+")
# "<Trivial name> <Xxx> <X> " -- trivial name is one capitalised word, optionally
# followed by " acid"; three-letter symbol is Capital + two lower-case; one-letter
# symbol is a single capital. Applied to the tag-stripped, whitespace-collapsed text.
_ROW = re.compile(r"\b([A-Z][a-z]+(?: acid)?) ([A-Z][a-z]{2}) ([A-Z]) ")


class _Http:
    def __init__(self, log):
        self.log = log
        self.session = requests.Session()
        self.session.headers["User-Agent"] = USER_AGENT
        self._last = 0.0

    def get(self, url: str) -> requests.Response:
        for attempt in range(1, RETRIES + 1):
            wait = max(0.0, MIN_INTERVAL_S - (time.monotonic() - self._last))
            if wait:
                time.sleep(wait)
            try:
                resp = self.session.get(url, timeout=60)
                self._last = time.monotonic()
            except requests.RequestException as exc:
                self.log.warning("attempt %d/%d: %s -> %s", attempt, RETRIES, url, exc)
                time.sleep(2 ** attempt)
                continue
            if resp.status_code in (429, 500, 502, 503, 504):
                self.log.warning("attempt %d/%d: HTTP %d for %s", attempt, RETRIES, resp.status_code, url)
                time.sleep(2 ** attempt)
                continue
            return resp
        raise RuntimeError(f"gave up after {RETRIES} attempts: {url}")


# ----------------------------------------------------------------------------
# Pure functions (tested offline)
# ----------------------------------------------------------------------------

def page_text(html_source: str) -> str:
    """Strip tags, unescape entities, collapse whitespace."""
    return _WS.sub(" ", html.unescape(_TAG.sub(" ", html_source))).strip()


def parse_symbol_table(text: str, letters: str = common.AMINO_ACIDS) -> dict[str, dict[str, str]]:
    """Return {one_letter: {trivial_name, three_letter}} for the standard letters.

    Every row matching the pattern is collected; rows whose one-letter symbol is not
    one of the standard letters (e.g. the 'Unspecified Xaa X' row) are ignored. The
    caller stops the run unless exactly the standard letters were found, once each.
    """
    found: dict[str, dict[str, str]] = {}
    duplicates = []
    for trivial, three, one in _ROW.findall(text):
        if one not in letters:
            continue
        if one in found:
            duplicates.append(one)
            continue
        found[one] = {"trivial_name": trivial, "three_letter": three}
    if duplicates:
        raise ValueError(f"one-letter symbol matched more than once in the table: {sorted(set(duplicates))}")
    return found


def parse_pubchem_properties(payload: dict) -> dict:
    """Return the single Properties row, or raise if there is not exactly one."""
    props = payload.get("PropertyTable", {}).get("Properties", [])
    if len(props) != 1:
        raise ValueError(f"expected exactly one compound, got {len(props)}")
    row = props[0]
    for key in ("CID", "MolecularFormula", "MolecularWeight"):
        if key not in row:
            raise ValueError(f"PubChem response lacks {key}")
    return row


def _decimals(value: str) -> int:
    return len(value.split(".", 1)[1]) if "." in value else 0


# ----------------------------------------------------------------------------
# Stage
# ----------------------------------------------------------------------------

def _pubchem_url(base: str, name: str, properties: str) -> str:
    return f"{base}/compound/name/{requests.utils.quote(name)}/property/{properties}/JSON"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    log = common.make_logger(STAGE)
    dec = common.read_ini(common.COMPOSITION_DECISIONS_INI)
    table_url = dec["amino_acid_symbols"]["url"]
    base, properties = dec["mass_source"]["base_url"], dec["mass_source"]["properties"]
    prefix = dec["query_name_rule"]["chiral_prefix"]
    extras = {k: v for k, v in dec["extra_compounds"].items() if k != "decision"}
    http = _Http(log)

    # 1. Symbol table from the IUPAC-IUBMB page.
    resp = http.get(table_url)
    resp.raise_for_status()
    common.IUPAC_DIR.mkdir(parents=True, exist_ok=True)
    page_path = common.IUPAC_DIR / "tab1.html"
    page_path.write_bytes(resp.content)
    sha = common.sha256_file(page_path)
    retrieved = common.iso_now()
    src = common.new_ini()
    src["iupac_table1"] = {"url": table_url, "retrieved": retrieved, "sha256": sha, "bytes": str(len(resp.content))}
    common.write_ini(src, common.IUPAC_DIR / "source.ini",
                     ["source.ini -- GENERATED by fetch_amino_acid_masses.py. DO NOT EDIT BY HAND.",
                      "Provenance of the downloaded IUPAC-IUBMB Table 1 page (config/composition_decisions.ini, D26)."])
    try:
        symbols = parse_symbol_table(page_text(resp.text))
    except ValueError as exc:
        log.error("symbol table parse: %s -- stopping", exc)
        return 1
    missing = sorted(set(common.AMINO_ACIDS) - set(symbols))
    if missing:
        log.error("symbol table parse found %d letters; missing %s -- page shape changed? stopping", len(symbols), missing)
        return 1
    sym = common.new_ini()
    for letter in common.AMINO_ACIDS:
        sym[letter] = symbols[letter]
    common.write_ini(sym, common.AMINO_ACID_SYMBOLS_INI, [
        "amino_acid_symbols.ini -- GENERATED by fetch_amino_acid_masses.py. DO NOT EDIT BY HAND.",
        f"source    = {table_url}",
        f"retrieved = {retrieved}",
        f"sha256    = {sha}",
        "Parsed rows: <trivial name> <three-letter symbol> <one-letter symbol>, for the twenty standard letters.",
    ])
    log.info("IUPAC Table 1: %d bytes, sha256 %s; parsed %d symbol rows", len(resp.content), sha, len(symbols))

    # 2. PubChem, one (or two) queries per section.
    raw_dir = common.PUBCHEM_DIR / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    served: dict[str, tuple[dict, str, str]] = {}   # section -> (row, query_name_used, query_rule)
    queries: list[tuple[str, list[tuple[str, str]]]] = []
    for letter in common.AMINO_ACIDS:
        name = symbols[letter]["trivial_name"]
        queries.append((letter, [(prefix + name, "chiral_prefix"), (name, "bare_trivial_name")]))
    for section, name in extras.items():
        queries.append((section, [(name, "extra_compound")]))

    for section, candidates in queries:
        resolved = None
        for name, rule in candidates:
            url = _pubchem_url(base, name, properties)
            resp = http.get(url)
            (raw_dir / f"{section}.{rule}.json").write_text(
                json.dumps({"url": url, "status": resp.status_code, "body": _json_or_text(resp)},
                           indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
            if resp.status_code == 404:
                log.info("%s: PubChem has no compound named %r (%s); trying next", section, name, rule)
                continue
            if resp.status_code != 200:
                log.error("%s (%r): HTTP %d -- stopping", section, name, resp.status_code)
                return 1
            try:
                row = parse_pubchem_properties(resp.json())
            except ValueError as exc:
                log.error("%s (%r): %s -- stopping", section, name, exc)
                return 1
            resolved = (row, name, rule)
            break
        if resolved is None:
            log.error("%s: no PubChem compound for any of %s -- stopping", section, [n for n, _ in candidates])
            return 1
        served[section] = resolved
        row, name, rule = resolved
        log.info("%s %r [%s] -> CID %s %s MW %s", section, name, rule, row["CID"], row["MolecularFormula"], row["MolecularWeight"])

    # 3. Table.
    water = float(served["water"][0]["MolecularWeight"])
    decimals_aa = max(_decimals(str(served[a][0]["MolecularWeight"])) for a in common.AMINO_ACIDS)
    decimals = max(_decimals(str(r[0]["MolecularWeight"])) for r in served.values())
    out = common.new_ini()
    for section, _ in queries:
        row, name, rule = served[section]
        s = {}
        if len(section) == 1:
            s.update(trivial_name=symbols[section]["trivial_name"], three_letter=symbols[section]["three_letter"])
        s.update(query_name=name, query_rule=rule, cid=str(row["CID"]), molecular_formula=row["MolecularFormula"],
                 iupac_name=row.get("IUPACName", ""), molecular_weight=str(row["MolecularWeight"]),
                 source_url=_pubchem_url(base, name, properties))
        if len(section) == 1:
            s["free_mass"] = str(row["MolecularWeight"])
            s["residue_mass"] = f"{float(row['MolecularWeight']) - water:.{decimals}f}"
        out[section] = s
    n_fallback = sum(1 for s in common.AMINO_ACIDS if served[s][2] == "bare_trivial_name")
    common.write_ini(out, common.AMINO_ACID_MASSES_INI, [
        "amino_acid_masses.ini -- GENERATED by fetch_amino_acid_masses.py. DO NOT EDIT BY HAND.",
        f"fetched            = {common.iso_now()}",
        f"source             = PubChem PUG REST, {base}/compound/name/<query_name>/property/{properties}/JSON",
        f"symbols            = data/iupac/amino_acid_symbols.ini (parsed from {table_url}, sha256 {sha})",
        f"query_rule         = '{prefix}' + trivial name; bare trivial name if PubChem has no such compound ({n_fallback} of 20 used the bare name)",
        "units              = g/mol, average molecular weight as served by PubChem (MolecularWeight)",
        f"served_decimals    = {decimals_aa}",
        f"served_decimals_extras = {decimals}",
        "(served_decimals: decimal places PubChem gave MolecularWeight to, for the twenty amino acids; _extras: including water and hydrogen)",
        "residue_mass       = free_mass - molecular_weight(water); the in-chain residue mass",
        "free_mass          = MolecularWeight of the free amino acid; the USDA / laboratory amino-acid-analysis convention",
    ])
    log.info("wrote %s (%d sections; water = %s; %d letters resolved by bare name)",
             common.AMINO_ACID_MASSES_INI, len(out.sections()), water, n_fallback)
    return 0


def _json_or_text(resp):
    try:
        return resp.json()
    except ValueError:
        return resp.text


if __name__ == "__main__":
    sys.exit(main())
