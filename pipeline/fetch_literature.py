#!/usr/bin/env python3
"""
fetch_literature.py — hash every literature file on disk and write the manifest.

Nothing here downloads (D80). A literature file is obtained by hand once, saved under
data/literature/<source_id>/, and hashed once; every later run re-hashes what is on disk and
must reproduce the pin. Fetch-by-code was retired because re-downloading a static document on
every run cost time, failed intermittently, and proved nothing the hash does not prove.

Reads  config/literature_sources.ini
Writes data/literature/manifest.ini        one section per file: url, path, sha256, bytes, obtained
       data/literature/manifest_map.tsv    one row per source, one column per category, X where used
       data/literature/README.md           human-readable provenance, generated
       logs/fetch_literature_<UTC>.log

A source with no file.<n>.* keys is cited but not stored — a method citation, a provenance
reference, or a paywalled article whose values were never transcribed. It appears in the map and
the README and is never flagged for a missing file: there is no file to miss.

Rules (the fetch-mechanics rules F3, F4, F4b and F6 are retired with the download code)
  F1  A url of `___` is a flag: recorded, skipped, non-zero exit at the end.
  F2  If config carries a sha256, the bytes on disk must hash to it; a mismatch is a flag.
      If config's sha256 is blank, the computed hash is written back into
      config/literature_sources.ini — the only thing this stage writes into config.
  F2b An unknown hash is blank, never `___`. `___` is only ever a value a person types.
  F5  Every file is obtained by hand and must already be on disk under
      data/literature/<source_id>/. This stage hashes it; the file on disk is the record.
  F7  The stored filename is the name the file was downloaded with, recorded in config as
      file.<n>.name (amended 2026-09-18, Step 7b). The URL is provenance only. For a section
      with no file.<n>.name the old rule still applies: the filename in the URL, with
      characters Windows forbids (colon, quotes, angle brackets, pipe, ?, *, slashes)
      replaced by '_'. Either way the name looked for is recorded in the manifest as
      served_name. A URL whose path ends in '/' has no filename (ACS supplement links do
      this), so such a file needs file.<n>.name.
  F7b When the URL's filename has no extension, the extension the file's own magic bytes imply
      is appended (%PDF -> .pdf, PK -> .zip, Rar! -> .rar, gzip -> .gz), so every stored file
      opens in the program a person would use. A file is looked for under the bare name and
      under the name with any of those extensions.

A file whose bytes are an HTML page while its name implies a document is still reported. That
check was F6's, applied to a served body; there is no served body now, so it is applied to the
bytes on disk under F5. It is what catches a saved paywall page filed under a .pdf name.

Usage
  python run.py fetch-literature
  python pipeline/fetch_literature.py --only murgia_2021
  python pipeline/fetch_literature.py --check      # hash and compare, write nothing
"""

from __future__ import annotations

import argparse
import configparser
import hashlib
import re
import sys
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline import common  # noqa: E402

ROOT = common.REPO_ROOT
CONFIG = common.LITERATURE_SOURCES_INI
LIT_DIR = common.literature_dir()
PLACEHOLDER = "___"

MANIFEST_INI = LIT_DIR / "manifest.ini"
MANIFEST_MAP_TSV = LIT_DIR / "manifest_map.tsv"

MAGIC = {".pdf": b"%PDF", ".zip": b"PK", ".xlsx": b"PK", ".docx": b"PK",
         ".rar": b"Rar!", ".gz": b"\x1f\x8b"}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def source_sections(cp: configparser.ConfigParser, only: str | None = None):
    for name in cp.sections():
        if not name.startswith("source."):
            continue
        sid = name[len("source."):]
        if only and sid != only:
            continue
        yield sid, cp[name]


def file_entries(section):
    """(n, published_name, url, sha256, obtained, name) for every file.<n>.* in the section.
    name is file.<n>.name — the filename as downloaded (F7) — or '' when the section has none."""
    ns = sorted({int(m.group(1)) for k in section
                 for m in [re.match(r"file\.(\d+)\.", k)] if m})
    for n in ns:
        yield (n,
               section.get(f"file.{n}.published_name", ""),
               (section.get(f"file.{n}.url", "") or "").strip(),
               (section.get(f"file.{n}.sha256", "") or "").strip(),
               (section.get(f"file.{n}.obtained", "manual") or "manual").strip(),
               (section.get(f"file.{n}.name", "") or "").strip())


def categories_used_for(section) -> list[str]:
    return [c.strip() for c in (section.get("used_for", "") or "").split(",") if c.strip()]


def looks_like_html(data: bytes) -> bool:
    head = data[:400].lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<head>" in head


def magic_ok(data: bytes, ext: str) -> bool:
    want = MAGIC.get(ext.lower())
    return True if want is None else data.startswith(want)


def filename_from_url(url: str) -> str:
    name = urllib.parse.urlsplit(url).path.rsplit("/", 1)[-1]
    return urllib.parse.unquote(name) or "download"


def safe_name(served: str) -> str:
    return re.sub(r'[:"<>|?*\\/]', "_", served)


def find_file(dest_dir: Path, served: str) -> Path | None:
    """The stored file under its safe name, with or without an implied extension (F7, F7b)."""
    if not dest_dir.is_dir():
        return None
    base = safe_name(served)
    if (dest_dir / base).is_file():
        return dest_dir / base
    for ext in MAGIC:
        if (dest_dir / (base + ext)).is_file():
            return dest_dir / (base + ext)
    return None


def manifest_map_rows(cp: configparser.ConfigParser) -> tuple[list[str], list[list[str]], list[str]]:
    """One row per source, one column per category, X where the source is used.

    The columns and their order come from [categories] in config, not from the data — so the
    map reads the way a person laid it out, and a new category is one line in config plus the
    used_for of the sources it is used for. A source that feeds no calculation carries the
    unused_marker in its used_for; the marker is printed in the `other` column, since a row of
    four blanks would read as an unfilled line rather than a deliberate one.

    The first column is the source id exactly as config spells it: `murgia_2021` is already the
    author-year name, and a prettier label would have to be derived from the citation prose.

    Rows are sorted A-Z by that id, NOT left in config order: the map is read beside the
    `data/literature/` folder listing, which the file explorer sorts A-Z. Same order, same scan.
    """
    cats = [c.strip() for c in cp["categories"]["columns"].split(",") if c.strip()]
    if "unused_marker" not in cp["categories"]:
        raise SystemExit("[STOP] [categories] in literature_sources.ini has no unused_marker; "
                         "the word printed for a source that feeds no calculation is a decision, not a default")
    marker = cp["categories"]["unused_marker"].strip()
    header = ["source"] + cats
    body, unknown = [], []
    for sid, sec in sorted(source_sections(cp), key=lambda pair: pair[0]):
        used = categories_used_for(sec)
        if used == [marker]:
            body.append([sid] + [marker if c == cats[-1] else "" for c in cats])
            continue
        for c in used:
            if c not in cats:
                unknown.append(f"source.{sid}: used_for names `{c}`, which is not in [categories] columns")
        body.append([sid] + ["X" if c in used else "" for c in cats])
    return header, body, unknown


def record_hashes(pending: list[tuple[str, int, str]]) -> None:
    """Fill blank `file.<n>.sha256 =` lines in place. The file is edited as text so comments,
    spacing and order are untouched; only those lines change."""
    nl = "\r\n" if b"\r\n" in CONFIG.read_bytes() else "\n"
    lines = CONFIG.read_text(encoding="utf-8").split("\n")
    section = None
    for i, ln in enumerate(lines):
        m = re.match(r"\s*\[(.+?)\]\s*$", ln)
        if m:
            section = m.group(1)
            continue
        m = re.match(r"(\s*file\.(\d+)\.sha256\s*=)\s*$", ln.rstrip("\r"))
        if m and section and section.startswith("source."):
            sid, n = section[len("source."):], int(m.group(2))
            for psid, pn, sha in pending:
                if psid == sid and pn == n:
                    lines[i] = f"{m.group(1)} {sha}"
    CONFIG.write_text("\n".join(lines), encoding="utf-8", newline=nl)


def write_readme(cp: configparser.ConfigParser, manifest: configparser.ConfigParser) -> None:
    out = [
        "# data/literature/ — provenance",
        "",
        "_Generated by `pipeline/fetch_literature.py` from `config/literature_sources.ini` and "
        "`data/literature/manifest.ini`. Do not edit; edit the config and rerun._",
        "",
        f"Generated: {utc_now()}",
        "",
        "Every file here was obtained by hand and is stored exactly as served (D80); each run "
        "re-hashes it. A file without a SHA-256 in the manifest has not been hashed and is not "
        "an admissible input. A source with no file table is cited but not stored.",
        "",
        "`manifest_map.tsv` beside this file says which category uses which source.",
        "",
    ]
    for sid, sec in sorted(source_sections(cp), key=lambda pair: pair[0]):
        used = ", ".join(categories_used_for(sec)) or "—"
        out.append(f"## {sid}  — {sec.get('id_short', '')}, role: `{sec.get('role', '')}` "
                   f"({sec.get('role_decision', '')}), used by: {used}")
        out.append("")
        out.append(f"- **Citation:** {sec.get('citation', '')}")
        for k in ("doi", "pmid", "pmc", "raw_data", "correction_doi", "species", "article_url"):
            if sec.get(k):
                out.append(f"- **{k}:** {sec.get(k)}")
        files = list(file_entries(sec))
        if files:
            out.append("")
            out.append("| n | Published name | URL | Path | Hashed | SHA-256 |")
            out.append("|---|---|---|---|---|---|")
            for n, pub, url, _, _obtained, _name in files:
                m = f"file.{sid}.file.{n}"
                if manifest.has_section(m):
                    ms = manifest[m]
                    out.append(f"| {n} | {pub} | {url} | `{ms['path']}` | {ms['retrieved']} | `{ms['sha256']}` |")
                else:
                    out.append(f"| {n} | {pub} | {url or PLACEHOLDER} | — | — | **not on disk** |")
        else:
            out.append("- *Cited only; no file stored.*")
        out.append("")
    common.write_text_file(LIT_DIR / "README.md", "\n".join(out) + "\n")


def build(log, only: str | None = None):
    """Hash what is on disk. Returns (manifest, flags, pending hashes, counts)."""
    cp = common.read_ini(CONFIG)
    manifest = configparser.ConfigParser(interpolation=None)
    manifest.optionxform = str
    flags: list[str] = []
    pending: list[tuple[str, int, str]] = []
    n_files = n_cited_only = 0

    # A-Z, like the map and the data/literature/ folder listing: the log is read beside them.
    for sid, sec in sorted(source_sections(cp, only), key=lambda pair: pair[0]):
        if not categories_used_for(sec):
            flags.append(f"source.{sid}: no used_for — every source names what it is used for, "
                         f"or the unused marker")
        files = list(file_entries(sec))
        if not files:
            n_cited_only += 1
            log.info("cited only, no file: %-24s role %s", sid, sec.get("role", ""))
            continue
        for n, pub_name, url, want_sha, obtained, name in files:
            key = f"{sid}.file.{n}"
            if not url or url == PLACEHOLDER:
                flags.append(f"{key}: url is blank or `{PLACEHOLDER}` (F1)")
                log.error("[FLAG] %s url is a placeholder", key)
                continue
            if name == PLACEHOLDER:
                flags.append(f"{key}: name is `{PLACEHOLDER}` — it is the filename as downloaded, or absent (F7)")
                log.error("[FLAG] %s name is a placeholder", key)
                continue
            # F7: the name as downloaded when config records it; the URL's filename otherwise
            served = name or filename_from_url(url)
            found = find_file(LIT_DIR / sid, served)
            if found is None:
                want = (LIT_DIR / sid / safe_name(served)).relative_to(ROOT).as_posix()
                hint = "" if Path(served).suffix else " (or that name plus the extension the browser gave it)"
                if not name and served == "download":   # filename_from_url found nothing after the last '/'
                    hint = f" — the URL has no filename; add file.{n}.name = <the filename as downloaded> (F7)"
                flags.append(f"{key}: not on disk — expected {want}{hint} (F5)")
                log.error("[FLAG] %s missing on disk", key)
                continue
            data = found.read_bytes()
            ext = found.suffix.lower()
            if looks_like_html(data) or not magic_ok(data, ext):
                flags.append(f"{key}: {found.name} is an HTML page or has the wrong magic bytes for {ext} (F5)")
                log.error("[FLAG] %s wrong file type on disk", key)
                continue
            got = sha256_bytes(data)
            if want_sha and want_sha.lower() != got:
                flags.append(f"{key}: sha256 mismatch — config {want_sha[:12]}… disk {got[:12]}… (F2)")
                log.error("[FLAG] %s sha256 mismatch", key)
                continue
            if not want_sha:
                pending.append((sid, n, got))
                log.info("%s: no pin in config; recording %s…", key, got[:12])
            manifest[f"file.{sid}.file.{n}"] = {
                "source_id": sid,
                "id_short": sec.get("id_short", ""),
                "published_name": pub_name,
                "url": url,
                "path": found.relative_to(ROOT).as_posix(),
                "served_name": served,
                "name_from": "config" if name else "url",
                "sha256": got,
                "bytes": str(len(data)),
                "content_type": obtained,
                "obtained": obtained,
                "retrieved": utc_now(),
            }
            n_files += 1
            log.info("[OK]   %-32s %9d bytes  sha256 %s…", key, len(data), got[:12])
    return cp, manifest, flags, pending, n_files, n_cited_only


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="source id, e.g. murgia_2021")
    ap.add_argument("--check", action="store_true", help="hash and compare; write nothing")
    args = ap.parse_args(argv)

    log = common.make_logger("fetch_literature")
    LIT_DIR.mkdir(parents=True, exist_ok=True)
    cp, manifest, flags, pending, n_files, n_cited_only = build(log, args.only)

    if args.check:
        log.info("check: %d file(s) hashed, %d source(s) cited only, %d flag(s)",
                 n_files, n_cited_only, len(flags))
    else:
        if pending:
            record_hashes(pending)
            log.info("recorded %d new hash(es) into %s", len(pending), CONFIG.relative_to(ROOT).as_posix())
        common.write_ini(manifest, MANIFEST_INI, [
            "Generated by pipeline/fetch_literature.py — DO NOT EDIT BY HAND.",
            f"generated = {utc_now()}",
            "Nothing is downloaded (D80): every file was obtained by hand and is re-hashed each run.",
        ])
        header, body, unknown = manifest_map_rows(cp)
        common.write_tsv(MANIFEST_MAP_TSV, header, body)
        flags.extend(unknown)
        log.info("wrote %s: %d source(s) x %d categor%s", MANIFEST_MAP_TSV.relative_to(ROOT).as_posix(),
                 len(body), len(header) - 1, "y" if len(header) == 2 else "ies")
        write_readme(cp, manifest)
        log.info("manifest: %d file(s), %d source(s) cited only", n_files, n_cited_only)

    if flags:
        log.error("%d flag(s):", len(flags))
        for f in flags:
            log.error("  %s", f)
        return 1
    log.info("no flags")
    return 0


if __name__ == "__main__":
    sys.exit(main())
