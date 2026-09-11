"""Shared helpers for every pipeline stage.

Nothing in this module knows anything about biology. It provides:
- repository paths,
- a timestamped file logger (every run leaves a record in outputs/logs/),
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
LOGS_DIR = OUTPUTS_DIR / "logs"
DOCS_DIR = REPO_ROOT / "docs"

DECISIONS_INI = CONFIG_DIR / "protein_set_decisions.ini"
RESOLVED_TERMS_INI = ONTOLOGY_DIR / "resolved_terms.ini"
POOL_QUERIES_INI = DATA_DIR / "pool_queries.ini"
SEQUENCES_INI = DATA_DIR / "uniprot_sequences.ini"
ACCESSIONS_INI = CONFIG_DIR / "accessions.ini"
SEGMENTS_INI = CONFIG_DIR / "segments.ini"
FLAGS_TSV = OUTPUTS_DIR / "flags.tsv"

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
    """Logger that writes to outputs/logs/<stage>_<timestamp>.log and to stderr."""
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
