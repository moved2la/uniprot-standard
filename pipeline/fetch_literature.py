#!/usr/bin/env python3
"""
fetch_literature.py — networked stage: download and hash every literature file.

Reads config/mass_fraction_decisions.ini, downloads every file listed under
[source.*] file.<n>.url, stores it UNCHANGED under data/literature/<source_id>/<filename from the
URL>, computes SHA-256, and
generates:

    data/literature/manifest.ini   — one section per file: url, path, sha256,
                                      retrieved (UTC), bytes, content-type
    data/literature/README.md      — human-readable provenance, generated
    logs/fetch_literature_<UTC>.log

Rules
  F1  A url of `___` is a flag: recorded, skipped, exit non-zero at the end.
  F2  If config carries a sha256 for a file, the downloaded bytes must hash to
      it; mismatch is a flag and nothing is written. If config sha256 is
      blank, the computed hash is written back into config (the ONLY thing
      this script writes into config).
  F3  Nothing is ever modified, converted, or re-saved. Bytes on disk == bytes
      served. A re-fetch whose bytes already equal the file on disk writes
      nothing (logged as unchanged); a file that cannot be written (locked by
      a viewer or a sync client) is a flag, not a crash.
  F4  HTTP failure is a flag, not a retry loop that hides the failure.
  F5  A file obtained by hand (blocked download, paywalled PDF) is declared
      with file.<n>.obtained = manual. It must already be saved as
      data/literature/<source_id>/<filename from the URL> — the name the
      browser gives it by default. The script hashes it and records
      obtained=manual. The file on disk is the record; the script never
      tries to download it.
  F7  Stored filenames are the filename in the URL, for downloaded and
      hand-obtained files alike, with characters Windows forbids (colon, quotes,
      angle brackets, pipe, ?, *, slashes) replaced by '_'. The exact served
      name is in the manifest.
  F7b When the URL's filename has no extension, the extension the file's own
      magic bytes imply is appended (%PDF -> .pdf, PK -> .zip, Rar! -> .rar,
      gzip -> .gz), so every stored file opens in the program a person would
      use. A hand-obtained file is looked for under the bare name and under
      the name with any of those extensions.
  F6  A served body that is HTML while the expected type is a document
      (pdf/xlsx/zip/rar/docx) is a blocked or paywalled page: flagged in the
      log, nothing written to disk.

Usage
  python pipeline/fetch_literature.py            # fetch everything with a url
  python pipeline/fetch_literature.py --only murgia_2021
  python pipeline/fetch_literature.py --dry-run  # list what would be fetched
"""

from __future__ import annotations

import argparse
import configparser
import hashlib
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "mass_fraction_decisions.ini"
LIT_DIR = ROOT / "data" / "literature"
LOG_DIR = ROOT / "logs"
USER_AGENT = "uniprot-standard fetch_literature.py (provenance fetch; contact via repository)"
PLACEHOLDER = "___"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def read_config() -> configparser.ConfigParser:
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=(";",))
    cp.optionxform = str  # keep key case
    with CONFIG.open(encoding="utf-8") as fh:
        cp.read_file(fh)
    return cp


def source_sections(cp: configparser.ConfigParser, only: str | None):
    for sec in cp.sections():
        if not sec.startswith("source."):
            continue
        sid = sec[len("source."):]
        if only and sid != only:
            continue
        yield sid, cp[sec]


def file_entries(section):
    """Yield (n, published_name, url, sha256, obtained) for file.<n>.* keys."""
    ns = sorted({int(m.group(1)) for k in section if (m := re.match(r"file\.(\d+)\.", k))})
    for n in ns:
        yield (
            n,
            section.get(f"file.{n}.published_name", "").strip(),
            section.get(f"file.{n}.url", "").strip(),
            section.get(f"file.{n}.sha256", "").strip(),
            section.get(f"file.{n}.obtained", "download").strip().lower(),
        )


DOC_EXT = {".pdf", ".xlsx", ".xls", ".zip", ".rar", ".docx", ".csv", ".tsv", ".txt"}
DOC_MAGIC = {
    ".pdf": (b"%PDF",),
    ".xlsx": (b"PK\x03\x04",), ".docx": (b"PK\x03\x04",), ".zip": (b"PK\x03\x04",),
    ".rar": (b"Rar!",),
}


def looks_like_html(data: bytes) -> bool:
    head = data[:512].lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<html" in head


def magic_ok(data: bytes, ext: str) -> bool:
    sigs = DOC_MAGIC.get(ext)
    return True if not sigs else any(data.startswith(sig) for sig in sigs)


def filename_from_url(url: str) -> str:
    """Filename as served (URL basename), percent-decoded."""
    path = urllib.parse.urlparse(url).path
    return urllib.parse.unquote(Path(path).name) or "download"


MAGIC_EXT = ((b"%PDF", ".pdf"), (b"PK\x03\x04", ".zip"), (b"Rar!", ".rar"), (b"\x1f\x8b", ".gz"))
KNOWN_EXT = {".pdf", ".xlsx", ".xls", ".zip", ".rar", ".docx", ".csv", ".tsv", ".txt", ".gz", ".html", ".htm", ".json", ".xml", ".fasta", ".obo", ".ini"}


def extension_from_bytes(data: bytes) -> str:
    """The extension a file's own magic bytes imply, or '' (F7b)."""
    for sig, ext in MAGIC_EXT:
        if data.startswith(sig):
            return ext
    return ""


def stored_name(served: str, data: bytes | None = None) -> str:
    """Served name with characters Windows forbids replaced by '_' (F7). When the served name has no
    extension and the bytes are given, the extension the magic bytes imply is appended (F7b), so a
    PDF served from an extensionless URL is stored as <basename>.pdf and opens like any other PDF."""
    name = re.sub(r'[<>:"/\\|?*]', "_", served)
    if data is not None and Path(name).suffix.lower() not in KNOWN_EXT:
        name += extension_from_bytes(data)          # ".2008" or ".E365" is part of the name, not an extension
    return name


def find_manual_file(dest_dir: Path, served: str) -> Path | None:
    """A hand-obtained file is looked for under the bare stored name and, when that name has no
    extension, under every extension the magic bytes could add (F7b) — the person saved it with the
    extension their browser gave it."""
    bare = dest_dir / stored_name(served)
    if bare.exists():
        return bare
    if bare.suffix.lower() not in KNOWN_EXT:
        for _, ext in MAGIC_EXT:
            cand = bare.with_name(bare.name + ext)
            if cand.exists():
                return cand
    return None


def fetch(url: str) -> tuple[bytes, str]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read(), resp.headers.get("Content-Type", "")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="source id, e.g. murgia_2021")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    cp = read_config()
    LIT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    log_path = LOG_DIR / f"fetch_literature_{started.replace(':', '')}.log"
    log_lines: list[str] = [f"fetch_literature.py  started {started}", f"config {CONFIG.relative_to(ROOT).as_posix()}", ""]
    flags: list[str] = []
    config_dirty = False
    pending_hashes: list[tuple[str, int, str]] = []   # (source_id, n, sha256)

    manifest = configparser.ConfigParser(interpolation=None)
    manifest.optionxform = str
    manifest_path = LIT_DIR / "manifest.ini"
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8") as fh:
            manifest.read_file(fh)

    for sid, sec in source_sections(cp, args.only):
        for n, pub_name, url, want_sha, obtained in file_entries(sec):
            key = f"{sid}.file.{n}"
            dest_dir = LIT_DIR / sid

            if obtained == "manual":
                if not url or url == PLACEHOLDER:
                    flags.append(f"{key}: obtained=manual needs a url so the filename is known (F7)")
                    log_lines.append(f"[FLAG] {key} manual without url")
                    continue
                served_name = filename_from_url(url)
                found = find_manual_file(dest_dir, served_name)
                if found is None:
                    want = dest_dir / stored_name(served_name)
                    hint = "" if want.suffix else " (or that name plus the extension the browser gave it, e.g. .pdf)"
                    flags.append(f"{key}: obtained=manual but {want.relative_to(ROOT).as_posix()}{hint} not found (F5)")
                    log_lines.append(f"[FLAG] {key} manual file missing")
                    continue
                dest = found
                data = dest.read_bytes()
                ext = dest.suffix.lower()
                if looks_like_html(data) or not magic_ok(data, ext):
                    flags.append(f"{key}: {dest.name} is not a {ext} document (HTML or wrong magic bytes) (F6)")
                    log_lines.append(f"[FLAG] {key} manual file wrong type")
                    continue
                got_sha = sha256_bytes(data)
                if want_sha and want_sha.lower() != got_sha:
                    flags.append(f"{key}: sha256 mismatch on manual file — config {want_sha[:12]}… got {got_sha[:12]}… (F2)")
                    log_lines.append(f"[FLAG] {key} manual sha256 mismatch")
                    continue
                if not want_sha:
                    pending_hashes.append((sid, n, got_sha))
                    config_dirty = True
                ctype = "manual"
                log_lines.append(f"[MAN]  {key} {len(data)} bytes sha256 {got_sha} <- {dest.relative_to(ROOT).as_posix()}")
            else:
                if not url or url == PLACEHOLDER:
                    flags.append(f"{key}: url is `{PLACEHOLDER}` — not fetched (F1)")
                    log_lines.append(f"[FLAG] {key} url placeholder")
                    continue
                if args.dry_run:
                    log_lines.append(f"[DRY] {key} -> {url}")
                    continue
                dest_dir.mkdir(parents=True, exist_ok=True)
                served_name = filename_from_url(url)
                try:
                    data, ctype = fetch(url)
                except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
                    flags.append(f"{key}: fetch failed — {e} (F4)")
                    log_lines.append(f"[FLAG] {key} fetch failed: {e}")
                    continue
                fname = stored_name(served_name, data)          # F7b: extension from the bytes when the URL has none
                dest = dest_dir / fname
                ext = dest.suffix.lower()
                if ext in DOC_EXT and (looks_like_html(data) or not magic_ok(data, ext)):
                    flags.append(f"{key}: server returned HTML/wrong type for {fname} — blocked or paywalled; nothing saved (F6). Download it in a browser into {dest_dir.relative_to(ROOT).as_posix()}/ and set file.{n}.obtained = manual.")
                    log_lines.append(f"[FLAG] {key} html instead of {ext}")
                    continue
                got_sha = sha256_bytes(data)
                if want_sha and want_sha.lower() != got_sha:
                    flags.append(f"{key}: sha256 mismatch — config {want_sha[:12]}… got {got_sha[:12]}… (F2); nothing saved; the served file has changed since it was pinned")
                    log_lines.append(f"[FLAG] {key} sha256 mismatch")
                    continue
                if dest.exists() and sha256_bytes(dest.read_bytes()) == got_sha:
                    # F3: the bytes on disk already equal the bytes served; nothing to write (and nothing to
                    # collide with a viewer or a sync client holding the file open).
                    log_lines.append(f"[OK]   {key} {len(data)} bytes sha256 {got_sha} == {dest.relative_to(ROOT).as_posix()} (unchanged on disk)")
                else:
                    try:
                        dest.write_bytes(data)
                    except OSError as e:
                        flags.append(f"{key}: cannot write {dest.relative_to(ROOT).as_posix()} — {e.__class__.__name__}: {e} (F4). "
                                     f"The file is locked or read-only (a PDF viewer or a sync client holding it open); close it and rerun.")
                        log_lines.append(f"[FLAG] {key} cannot write: {e}")
                        continue
                    log_lines.append(f"[OK]   {key} {len(data)} bytes sha256 {got_sha} -> {dest.relative_to(ROOT).as_posix()}")
                if not want_sha:
                    pending_hashes.append((sid, n, got_sha))
                    config_dirty = True

            m = f"file.{key}"
            if not manifest.has_section(m):
                manifest.add_section(m)
            manifest[m]["source_id"] = sid
            manifest[m]["id_short"] = sec.get("id_short", "")
            manifest[m]["published_name"] = pub_name
            manifest[m]["url"] = url
            manifest[m]["path"] = dest.relative_to(ROOT).as_posix()
            manifest[m]["served_name"] = served_name
            manifest[m]["sha256"] = got_sha
            manifest[m]["bytes"] = str(len(data))
            manifest[m]["content_type"] = ctype
            manifest[m]["obtained"] = obtained
            manifest[m]["retrieved"] = utc_now()

    if not args.dry_run:
        with manifest_path.open("w", encoding="utf-8") as fh:
            fh.write(f"# Generated by pipeline/fetch_literature.py — DO NOT EDIT BY HAND\n# generated = {utc_now()}\n\n")
            manifest.write(fh)
        write_readme(cp, manifest)
        if config_dirty:
            record_hashes(pending_hashes)
            log_lines.append("config updated: sha256 recorded for newly fetched files (F2); nothing else in the file touched")

    log_lines.append("")
    log_lines.append(f"flags: {len(flags)}")
    log_lines.extend(f"  - {f}" for f in flags)
    log_lines.append(f"finished {utc_now()}")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print("\n".join(log_lines))
    return 1 if flags else 0


def record_hashes(pending: list[tuple[str, int, str]]) -> None:
    """Fill blank `file.<n>.sha256 =` lines in place. The file is edited as text so
    comments, spacing and order are untouched; only those lines change."""
    lines = CONFIG.read_text(encoding="utf-8").splitlines(keepends=True)
    section = None
    for i, ln in enumerate(lines):
        m = re.match(r"\s*\[(.+?)\]\s*$", ln)
        if m:
            section = m.group(1)
            continue
        m = re.match(r"(\s*file\.(\d+)\.sha256\s*=)\s*$", ln)
        if m and section and section.startswith("source."):
            sid, n = section[len("source."):], int(m.group(2))
            for psid, pn, sha in pending:
                if psid == sid and pn == n:
                    lines[i] = f"{m.group(1)} {sha}\n"
    CONFIG.write_text("".join(lines), encoding="utf-8")


def write_readme(cp: configparser.ConfigParser, manifest: configparser.ConfigParser) -> None:
    out = [
        "# data/literature/ — provenance",
        "",
        "_Generated by `pipeline/fetch_literature.py` from `config/mass_fraction_decisions.ini` and `data/literature/manifest.ini`. Do not edit; edit the config and rerun._",
        "",
        f"Generated: {utc_now()}",
        "",
        "Every file here is stored exactly as served. A file without a SHA-256 in the manifest has not been fetched and is not an admissible input.",
        "",
    ]
    for sid, sec in source_sections(cp, None):
        out.append(f"## {sid}  — {sec.get('id_short', '')}, role: `{sec.get('role', '')}` ({sec.get('role_decision', '')})")
        out.append("")
        out.append(f"- **Citation:** {sec.get('citation', '')}")
        for k in ("doi", "pmid", "pmc", "raw_data", "correction_doi", "species", "article_url"):
            if sec.get(k):
                out.append(f"- **{k}:** {sec.get(k)}")
        files = list(file_entries(sec))
        if files:
            out.append("")
            out.append("| n | Published name | URL | Path | Retrieved | SHA-256 |")
            out.append("|---|---|---|---|---|---|")
            for n, pub, url, _, obtained in files:
                m = f"file.{sid}.file.{n}"
                if manifest.has_section(m):
                    ms = manifest[m]
                    how = " (obtained by hand)" if ms.get("obtained", "download") == "manual" else ""
                    out.append(f"| {n} | {pub} | {url} | `{ms['path']}`{how} | {ms['retrieved']} | `{ms['sha256']}` |")
                else:
                    out.append(f"| {n} | {pub} | {url or '___'} | — | — | **not fetched** |")
        out.append("")
    (LIT_DIR / "README.md").write_text("\n".join(out) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())