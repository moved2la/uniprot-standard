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

# Literature: the sources are declared in config, the files sit on disk, the manifest is
# generated from the two by `python run.py fetch-literature`. Nothing is downloaded (D80).
LITERATURE_SOURCES_INI = CONFIG_DIR / "literature_sources.ini"   # hand-written: one section per source


def literature_dir() -> Path:
    """data/literature/ — one folder per source.

    A function, not a constant: DATA_DIR is redirected by tests and, later, by a category,
    and a path computed at import time cannot follow. Anything derived from a root that can
    move is resolved when it is used.
    """
    return DATA_DIR / "literature"


def literature_manifest() -> Path:
    """The generated manifest: one section per stored file, hashes and all."""
    return literature_dir() / "manifest.ini"


def literature_manifest_map() -> Path:
    """The generated source x category map."""
    return literature_dir() / "manifest_map.tsv"

# --------------------------------------------------------------------------- outputs/ layout
#
# Three meanings, one folder each:
#   outputs/intermediate/  what feeds the standard but is not the standard
#   outputs/standard/      the deliverable tables, at the top; their companions in subfolders
#   outputs/<the rest>     the food side (usda), the measured comparison, the match scores, flags
#
# A new category (blood, liver) gets outputs/<category>/ with the same intermediate/ and
# standard/ inside; category_outputs() below builds those paths so no folder name is typed
# into a stage. Every stage resolves its output paths through the constants here — nothing
# joins "outputs" to a name of its own.

# These are FUNCTIONS, not constants. A path computed at import time from OUTPUTS_DIR cannot
# follow a caller that redirects the root — a test writing to a temp folder, or a category
# writing under outputs/<category>/. That is not hypothetical: it broke standard_path in
# delivery 1 and the literature manifest in delivery 2, both times as a failure far from the
# cause. Anything derived from a root that can move is resolved when it is used.


def intermediate_dir() -> Path:
    """outputs/intermediate/ — what feeds the standard but is not the standard."""
    return OUTPUTS_DIR / "intermediate"


def protein_set_out_dir() -> Path:
    return intermediate_dir() / "protein_set"


def composition_dir() -> Path:
    return intermediate_dir() / "composition"


def digest_dir() -> Path:
    return intermediate_dir() / "digest"


def literature_inventory_dir() -> Path:
    return intermediate_dir() / "literature_inventory"


def mass_fractions_dir() -> Path:
    return intermediate_dir() / "mass_fractions"


def standard_dir() -> Path:
    """outputs/standard/ — the deliverable tables at the top, companions in subfolders."""
    return OUTPUTS_DIR / "standard"


def category_outputs(category: str) -> dict[str, Path]:
    """The outputs/<category>/ tree for a non-muscle category, same shape as the muscle tree."""
    root = OUTPUTS_DIR / category
    return {"root": root, "intermediate": root / "intermediate", "standard": root / "standard"}


# --------------------------------------------------------------------------- per-category files (Step 7b)
#
# A non-muscle category keeps its decisions and generated config under config/<category>/, its
# pool tables under data/<category>/, and its outputs under outputs/<category>/. Muscle's files
# stay where they are (D89). category_files() is the one place that says where a category's
# protein-set files live; a stage passed --category resolves everything through it. With no
# category the muscle constants are returned, so the muscle path is unchanged.
#
# Adaptation list, Step 7b: this is the ad-hoc form. Step 7c makes the category an argument
# of every stage and folds the muscle constants into the same map.


def category_config_dir(category: str) -> Path:
    return CONFIG_DIR / category


def category_data_dir(category: str) -> Path:
    return DATA_DIR / category


def category_files(category: str | None) -> dict[str, Path]:
    """Where a category's protein-set files live. None = muscle (the flat layout)."""
    if category is None:
        return {
            "decisions": DECISIONS_INI,
            "accessions": ACCESSIONS_INI,
            "segments": SEGMENTS_INI,
            "pool_dir": DATA_DIR,
            "pool_queries": POOL_QUERIES_INI,
            "flags": FLAGS_TSV,
            "protein_set_out": protein_set_out_dir(),
            "composition_dir": composition_dir(),
            "composition_tsv": composition_tsv(),
            "composition_summary": composition_summary_ini(),
        }
    cfg = category_config_dir(category)
    dat = category_data_dir(category)
    out = category_outputs(category)
    return {
        "decisions": cfg / "protein_set_decisions.ini",
        "accessions": cfg / "accessions.ini",
        "segments": cfg / "segments.ini",
        "pool_dir": dat,
        "pool_queries": dat / "pool_queries.ini",
        "flags": out["root"] / "flags.tsv",
        "protein_set_out": out["intermediate"] / "protein_set",
        "composition_dir": out["intermediate"] / "composition",
        "composition_tsv": out["intermediate"] / "composition" / "amino_acid_composition_per_protein.tsv",
        "composition_summary": out["intermediate"] / "composition" / "composition_summary.ini",
    }


def pool_file(category: str | None, pool: str) -> Path:
    """The pool table of one pool: data/tier<N>_pool.tsv for muscle, data/<category>/pool_<name>.tsv otherwise."""
    if category is None:
        return DATA_DIR / f"tier{pool}_pool.tsv"
    return category_data_dir(category) / f"pool_{pool}.tsv"


def pool_names(category: str | None) -> list[str]:
    """The pools a category declares: muscle's [tier.N] sections, or a category's [pool.<name>] sections."""
    cp = read_ini(category_files(category)["decisions"])
    if category is None:
        return tiers_from_decisions(cp)
    names = [s.split(".", 1)[1] for s in cp.sections() if s.startswith("pool.")]
    if not names:
        raise ValueError(f"no [pool.<name>] sections in {category_files(category)['decisions']}")
    return names


# --------------------------------------------------------------------------- the raw UniProt entries
#
# data/uniprot_raw/ holds one JSON per fetched entry. Since Step 7b it is nested one folder per
# category: uniprot_raw/muscle/ (the author moved muscle's files there by hand) and
# uniprot_raw/<category>/ for the entries THAT CATEGORY'S FETCH ADDED. An entry is fetched once —
# an accession already in the store is not fetched again by a later category — so a file's
# folder says which fetch brought it in, not which categories use it. Readers therefore look in
# every folder (and in the root, so files that were never moved are still found). Nothing here
# is ever deleted by code (D61).


MUSCLE_CATEGORY = "skeletal_muscle"     # the muscle category's name, as [categories] in literature_sources.ini spells it


def uniprot_raw_dir(category: str | None = None) -> Path:
    """Where a fetch WRITES: uniprot_raw/<category>/; the muscle fetch writes under the muscle category's name."""
    return UNIPROT_RAW_DIR / (MUSCLE_CATEGORY if category is None else category)


def uniprot_raw_entry(accession: str) -> Path | None:
    """Where an entry's JSON IS: the root, else the first subfolder (A-Z) that holds it; None if nowhere."""
    root = UNIPROT_RAW_DIR
    p = root / f"{accession}.json"
    if p.exists():
        return p
    if root.exists():
        for sub in sorted(d for d in root.iterdir() if d.is_dir()):
            q = sub / f"{accession}.json"
            if q.exists():
                return q
    return None

# Composition (Layer A): masses fetched from PubChem, PTM vocabulary from UniProt.
COMPOSITION_DECISIONS_INI = CONFIG_DIR / "composition_decisions.ini"   # hand-written: two URLs and one rule
IUPAC_DIR = DATA_DIR / "iupac"
AMINO_ACID_SYMBOLS_INI = IUPAC_DIR / "amino_acid_symbols.ini"           # parsed from the IUPAC-IUBMB table
PUBCHEM_DIR = DATA_DIR / "pubchem"
AMINO_ACID_MASSES_INI = PUBCHEM_DIR / "amino_acid_masses.ini"           # fetched from PubChem
PTMLIST_TXT = DATA_DIR / "uniprot_ptmlist" / "ptmlist.txt"               # fetched from UniProt
def composition_tsv() -> Path:
    return composition_dir() / "amino_acid_composition_per_protein.tsv"


def composition_summary_ini() -> Path:
    return composition_dir() / "composition_summary.ini"

# Mass fractions (Layer B): the generated weights table (D61) that aggregation reads.
MASS_FRACTIONS_TSV = CONFIG_DIR / "mass_fractions_per_entry.tsv"

# The standard (aggregation, uncertainty, stress).
EAA_INI = CONFIG_DIR / "fao_2013_indispensable_amino_acids.ini"    # hand-written transcription (D62)
UNCERTAINTY_SETTINGS_INI = CONFIG_DIR / "uncertainty_settings.ini"   # hand-written: draws and seed (D63)
STRESS_SETTINGS_INI = CONFIG_DIR / "stress_test_settings.ini"        # hand-written: magnitudes (D65)
FIBER_TYPE_MIX_INI = CONFIG_DIR / "fiber_type_mix.ini"               # hand-written transcription: the shares that make the final column
NON_PROTEIN_METABOLITE_POOLS_INI = CONFIG_DIR / "non_protein_metabolite_pools.ini"           # hand-written transcription: the non-protein metabolite pools folded into the standard (D73-D77, D82)
TIER_NAMES = {"1": "contractile", "2": "builders"}                     # naming rule: tiers are named, never numbered, in outputs (D66)

# Match Rate (Step 6): the USDA food tables.
USDA_CONFIG_INI = CONFIG_DIR / "usda_food_data.ini"    # hand-written: which archives, which release, extra nutrients
USDA_DIR = DATA_DIR / "usda"                            # the downloaded archives (gitignored) and their manifest
USDA_OUT_DIR = OUTPUTS_DIR / "usda"

# Step 6: the calculated standard beside an independent measurement.
GORISSEN_COMPARISON_INI = CONFIG_DIR / "gorissen_2018_comparison.ini"   # hand-written transcription (D79, D84)
COMPARISON_OUT_DIR = OUTPUTS_DIR / "comparison"

# --- Match Rate (D86-D88) ---
MATCH_RATE_INI = CONFIG_DIR / "match_rate.ini"                                   # hand-written: references, scored sets, food tables
FOOD_OTHER_SOURCES_CSV = CONFIG_DIR / "food_amino_acids_other_sources.csv"      # hand-maintained: foods USDA does not carry
MATCH_OUT_DIR = OUTPUTS_DIR / "match"

# --------------------------------------------------------------------------- placement of standard/ files
#
# The deliverable tables and the profiles they are read beside stay at the top of
# outputs/standard/. Everything that supports them sits in the subfolder named for what it is.
# This map IS the placement table: one place to read, one place to change. A file not listed
# here stays at the top — so adding a deliverable needs no entry, and moving one needs one line.

STANDARD_SUBFOLDERS = {
    # uncertainty/ — the interval reported beside the standard, and the terms behind it (D63)
    "uncertainty_intervals.tsv": "uncertainty",
    "uncertainty_per_amino_acid.tsv": "uncertainty",
    "uncertainty_summary.ini": "uncertainty",
    "sd_to_median_ratio.tsv": "uncertainty",
    "bounds_after_weighting.tsv": "uncertainty",
    # stress/ — named perturbations and how far each moves the standard (A9, D65)
    "stress_shifts.tsv": "stress",
    "stress_influence_per_entry.tsv": "stress",
    "stress_summary.ini": "stress",
    "stress_summary_per_scenario.tsv": "stress",
    "composition_distance_top_entries.tsv": "stress",
    # sensitivity/ — what the standard does under a cited range of an input
    "sensitivity_mhc_actin_profiles.tsv": "sensitivity",
    "sensitivity_mhc_actin_spread.tsv": "sensitivity",
    "sensitivity_non_protein_metabolite_pool_spread.tsv": "sensitivity",
    "completeness_sensitivity.tsv": "sensitivity",
}

# A plot goes where its table goes.
STANDARD_PLOT_SUBFOLDERS = {
    "uncertainty_terms_I": "uncertainty",
    "uncertainty_terms_IIa": "uncertainty",
    "uncertainty_terms_IIx": "uncertainty",
    "sensitivity_mhc_actin": "sensitivity",
}


def standard_path(name: str, base: Path | None = None) -> Path:
    """Where a standard/ file lives: its subfolder if it has one, else the top.

    `base` is the standard/ folder to resolve against and defaults to this repository's.
    Every stage passes its own OUT_DIR, so a test that redirects a stage to a temp folder,
    and a category that writes to outputs/<category>/standard/, both land where they mean to.
    The map holds subfolder NAMES, not absolute paths, for the same reason: an absolute path
    baked in at import time ignores where the caller is writing.
    """
    base = standard_dir() if base is None else base
    sub = STANDARD_SUBFOLDERS.get(name)
    return (base / sub / name) if sub else (base / name)


def standard_plot_path(name: str, ext: str, base: Path | None = None) -> Path:
    """Where a standard/ plot lives: in its table's subfolder, else standard/plots/."""
    base = standard_dir() if base is None else base
    sub = STANDARD_PLOT_SUBFOLDERS.get(name)
    return ((base / sub / "plots") if sub else (base / "plots")) / f"{name}.{ext}"


def write_text_file(path: Path, text: str) -> None:
    """Write a generated file, creating its folder. The one writer for files in subfolders."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")

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


# Keys that record WHEN a generated line was written, not WHAT it says. They change on every
# run of the stage that writes them, so hashing them makes an unchanged input look changed —
# the D81 problem: re-hashing a literature file that has not moved would mark every downstream
# output stale. A hash answers "would this produce the same numbers", and a timestamp cannot.
VOLATILE_KEYS = ("retrieved", "generated")


def hash_ini_sections(path: Path, sections: list[str]) -> str:
    """SHA-256 of just the named sections of an .ini, rendered canonically.

    A stage's header records what the stage READ, so that a reader can tell whether a rerun
    would produce the same numbers. Hashing a whole shared file breaks that: adding a source to
    data/literature/manifest.ini changed the hash recorded by every stage that touches the
    manifest, including ones that never read the new source — the header then claimed an input
    had changed when nothing the stage uses had. (A0, amended in Step 7a.)

    Sections are sorted and keys sorted within them, so the hash depends on the content read and
    not on where the section sits in the file or what order the stage asked for them. Keys in
    VOLATILE_KEYS are left out: they say when a line was written, not what it says. A section
    named here that does not exist is recorded as absent rather than skipped, so a section
    disappearing changes the hash instead of quietly matching.
    """
    cp = read_ini(path)
    parts = []
    for name in sorted(sections):
        if not cp.has_section(name):
            parts.append(f"[{name}]\n<absent>\n")
            continue
        body = "".join(f"{k} = {v}\n" for k, v in sorted(cp[name].items())
                       if k not in VOLATILE_KEYS)
        parts.append(f"[{name}]\n{body}")
    return hashlib.sha256("".join(parts).encode("utf-8")).hexdigest()


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
    order = ["fetch-literature", "protein-set", "composition", "mass-fractions", "standard",
             "usda", "comparison", "match", "blood-protein-set", "blood-composition"]
    current = os.environ.get("UNIPROT_STANDARD_COMMAND", "")
    return current in order and command in order and order.index(current) < order.index(command)