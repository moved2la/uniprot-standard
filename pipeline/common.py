"""Shared helpers for every pipeline stage.

Nothing in this module knows anything about biology. It provides:
- repository paths,
- a timestamped file logger (every run leaves a record in logs/),
- .ini reading/writing with a fixed style,
- SHA-256 / MD5 helpers.
"""

from __future__ import annotations

import configparser
import datetime as _dt
import hashlib
import logging
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = REPO_ROOT / "config"
DATA_DIR = REPO_ROOT / "data"
ONTOLOGY_DIR = DATA_DIR / "gene-ontology"
UNIPROT_RAW_DIR = DATA_DIR / "uniprot_raw"
OUTPUTS_DIR = REPO_ROOT / "outputs"
LOGS_DIR = REPO_ROOT / "logs"
DOCS_DIR = REPO_ROOT / "docs"

DECISIONS_INI = CONFIG_DIR / "protein_set_decisions.ini"
RESOLVED_TERMS_INI = ONTOLOGY_DIR / "resolved_terms.ini"
POOL_QUERIES_INI = DATA_DIR / "pool_queries.ini"
SEQUENCES_INI = DATA_DIR / "uniprot_sequences.ini"
ACCESSIONS_INI = CONFIG_DIR / "accessions.ini"
SEGMENTS_INI = CONFIG_DIR / "segments.ini"
FLAGS_TSV = OUTPUTS_DIR / "flags.tsv"

# Composition (Layer A): masses fetched from PubChem, PTM vocabulary from UniProt.
COMPOSITION_DECISIONS_INI = CONFIG_DIR / "composition_decisions.ini"   # hand-written: two URLs and one rule
IUPAC_DIR = DATA_DIR / "iupac"
AMINO_ACID_SYMBOLS_INI = IUPAC_DIR / "amino_acid_symbols.ini"           # parsed from the IUPAC-IUBMB table
PUBCHEM_DIR = DATA_DIR / "pubchem"
AMINO_ACID_MASSES_INI = PUBCHEM_DIR / "amino_acid_masses.ini"           # fetched from PubChem
PTMLIST_TXT = DATA_DIR / "uniprot_ptmlist" / "ptmlist.txt"               # fetched from UniProt
COMPOSITION_DIR = OUTPUTS_DIR / "composition"
COMPOSITION_TSV = COMPOSITION_DIR / "amino_acid_composition_per_protein.tsv"
COMPOSITION_SUMMARY_INI = COMPOSITION_DIR / "composition_summary.ini"

# Mass fractions (Layer B): the generated weights table (D61) that aggregation reads.
MASS_FRACTIONS_TSV = CONFIG_DIR / "mass_fractions_per_entry.tsv"
MASS_FRACTIONS_DIR = OUTPUTS_DIR / "mass_fractions"

# The standard (aggregation, uncertainty, stress).
EAA_INI = CONFIG_DIR / "fao_2013_indispensable_amino_acids.ini"    # hand-written transcription (D62)
UNCERTAINTY_SETTINGS_INI = CONFIG_DIR / "uncertainty_settings.ini"   # hand-written: draws and seed (D63)
STRESS_SETTINGS_INI = CONFIG_DIR / "stress_test_settings.ini"        # hand-written: magnitudes (D65)
FIBER_TYPE_MIX_INI = CONFIG_DIR / "fiber_type_mix.ini"               # hand-written transcription: the shares that make the final column
NON_PROTEIN_METABOLITE_POOLS_INI = CONFIG_DIR / "non_protein_metabolite_pools.ini"           # hand-written transcription: the non-protein metabolite pools folded into the standard (D73-D77, D82)
TIER_NAMES = {"1": "contractile", "2": "builders"}                     # naming rule: tiers are named, never numbered, in outputs (D66)
STANDARD_DIR = OUTPUTS_DIR / "standard"

# Match Rate (Step 6): the USDA food tables.
USDA_CONFIG_INI = CONFIG_DIR / "usda_food_data.ini"    # hand-written: which archives, which release, extra nutrients
USDA_DIR = DATA_DIR / "usda"                            # the downloaded archives (gitignored) and their manifest
USDA_OUT_DIR = OUTPUTS_DIR / "usda"

# The twenty standard amino acid letters, in the order used for every vector.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def utc_now() -> _dt.datetime:
    return _dt.datetime.now(_dt.timezone.utc)


def timestamp() -> str:
    """Filesystem-safe UTC timestamp."""
    return utc_now().strftime("%Y%m%dT%H%M%SZ")


def iso_now() -> str:
    return utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")


def today() -> str:
    return utc_now().strftime("%Y-%m-%d")


def make_logger(stage: str, logs_dir: Path = LOGS_DIR) -> logging.Logger:
    """Logger that writes to logs/<stage>_<timestamp>.log and to stderr."""
    logs_dir.mkdir(parents=True, exist_ok=True)
    path = logs_dir / f"{stage}_{timestamp()}.log"
    logger = logging.getLogger(f"pipeline.{stage}.{timestamp()}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    fmt = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s", "%Y-%m-%dT%H:%M:%S")
    fmt.converter = time.gmtime  # log lines and log filenames are both UTC
    fh = logging.StreamHandler(open(path, "a", encoding="utf-8", newline="\n"))
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.info("log file: %s", path)
    return logger


def read_ini(path: Path) -> configparser.ConfigParser:
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=(";", "#"))
    cp.optionxform = str  # keep key case
    if not path.exists():
        raise FileNotFoundError(f"required file missing: {path}")
    cp.read(path, encoding="utf-8")
    return cp


def new_ini() -> configparser.ConfigParser:
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    return cp


def render_ini(cp: configparser.ConfigParser, header_lines: list[str]) -> str:
    """Render an .ini with a comment header. Keys are aligned for readability."""
    lines = [f"# {line}" if line else "#" for line in header_lines]
    lines.append("")
    for section in cp.sections():
        lines.append(f"[{section}]")
        items = list(cp.items(section))
        width = max((len(k) for k, _ in items), default=0)
        for k, v in items:
            lines.append(f"{k.ljust(width)} = {v}")
        lines.append("")
    return "\n".join(lines)


def write_ini(cp: configparser.ConfigParser, path: Path, header_lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_ini(cp, header_lines), encoding="utf-8", newline="\n")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_text(text: str) -> str:
    return hashlib.md5(text.encode("ascii")).hexdigest().upper()


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required file missing: {path}")
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(dict(zip(header, line.split("\t"))))
    return rows


def tiers_from_decisions(cp: configparser.ConfigParser) -> list[str]:
    """Return tier numbers declared as [tier.N] sections, in order."""
    tiers = sorted(s.split(".", 1)[1] for s in cp.sections() if s.startswith("tier."))
    if not tiers:
        raise ValueError("no [tier.N] sections in protein_set_decisions.ini")
    return tiers


def tier_definition(cp: configparser.ConfigParser, tier: str) -> str:
    """'ontology' (lookup_name present) or the value of `definition` (e.g. 'measured_remainder')."""
    sec = cp[f"tier.{tier}"]
    return sec.get("definition", "ontology") if not sec.get("lookup_name") else "ontology"


def ontology_tiers(cp: configparser.ConfigParser) -> list[str]:
    return [t for t in tiers_from_decisions(cp) if tier_definition(cp, t) == "ontology"]


def measured_tiers(cp: configparser.ConfigParser) -> list[str]:
    return [t for t in tiers_from_decisions(cp) if tier_definition(cp, t) == "measured_remainder"]


def running_command_is_before(command: str) -> bool:
    """True when run.py is executing an EARLIER command than `command` (env UNIPROT_STANDARD_COMMAND),
    so `command`'s outputs are legitimately stale until it is run next. False for `python run.py test`."""
    import os
    order = ["protein-set", "composition", "mass-fractions", "standard", "usda"]
    current = os.environ.get("UNIPROT_STANDARD_COMMAND", "")
    return current in order and command in order and order.index(current) < order.index(command)
