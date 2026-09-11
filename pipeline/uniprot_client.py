"""Thin client for rest.uniprot.org.

Handles retries with backoff, a polite request rate, cursor pagination on the
search endpoint, and capture of the X-UniProt-Release header. Contains no
biology and no query content -- callers pass the query strings.
"""

from __future__ import annotations

import re
import time

import requests

BASE = "https://rest.uniprot.org/uniprotkb"
PAGE_SIZE = 500
MIN_INTERVAL_S = 0.34  # ~3 requests/second
RETRIES = 5

_NEXT_LINK = re.compile(r'<([^>]+)>;\s*rel="next"')


class UniProtClient:
    def __init__(self, log, user_agent: str = "uniprot-standard/protein-set (contact: repository README)"):
        self.log = log
        self.session = requests.Session()
        self.session.headers["User-Agent"] = user_agent
        self.release: str | None = None
        self._last = 0.0

    def _get(self, url: str, params: dict | None = None) -> requests.Response:
        for attempt in range(1, RETRIES + 1):
            wait = max(0.0, MIN_INTERVAL_S - (time.monotonic() - self._last))
            if wait:
                time.sleep(wait)
            try:
                resp = self.session.get(url, params=params, timeout=120)
                self._last = time.monotonic()
            except requests.RequestException as exc:
                self.log.warning("attempt %d/%d: %s -> %s", attempt, RETRIES, url, exc)
                time.sleep(2 ** attempt)
                continue
            rel = resp.headers.get("X-UniProt-Release")
            if rel:
                if self.release and rel != self.release:
                    self.log.error("UniProt release changed mid-run: %s -> %s", self.release, rel)
                    raise RuntimeError("UniProt release changed during the run; rerun from the start")
                self.release = rel
            if resp.status_code in (429, 500, 502, 503, 504):
                self.log.warning("attempt %d/%d: HTTP %d for %s", attempt, RETRIES, resp.status_code, url)
                time.sleep(2 ** attempt)
                continue
            return resp
        raise RuntimeError(f"gave up after {RETRIES} attempts: {url}")

    def search_tsv(self, query: str, fields: list[str]) -> list[dict[str, str]]:
        """Run a search and return every row (follows cursor pagination)."""
        url = f"{BASE}/search"
        params = {"query": query, "format": "tsv", "fields": ",".join(fields), "size": PAGE_SIZE}
        rows: list[dict[str, str]] = []
        total_reported = None
        while url:
            resp = self._get(url, params)
            resp.raise_for_status()
            if total_reported is None:
                total_reported = resp.headers.get("X-Total-Results")
            lines = resp.text.splitlines()
            if lines:
                header = lines[0].split("\t")
                for line in lines[1:]:
                    if line.strip():
                        rows.append(dict(zip(header, line.split("\t"))))
            m = _NEXT_LINK.search(resp.headers.get("Link", ""))
            url = m.group(1) if m else None
            params = None  # the next link carries its own parameters
        self.log.info("query %r -> %d rows (X-Total-Results=%s)", query, len(rows), total_reported)
        if total_reported is not None and total_reported.isdigit() and int(total_reported) != len(rows):
            raise RuntimeError(f"pagination mismatch for {query!r}: header says {total_reported}, got {len(rows)}")
        return rows

    def entry_json(self, accession: str) -> dict:
        resp = self._get(f"{BASE}/{accession}.json")
        resp.raise_for_status()
        return resp.json()

    def fasta(self, identifier: str) -> tuple[str, str]:
        """Return (header_line, sequence) for an accession or an isoform id of the form <accession>-<n>."""
        resp = self._get(f"{BASE}/{identifier}.fasta")
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        if not lines or not lines[0].startswith(">"):
            raise RuntimeError(f"no FASTA record returned for {identifier}")
        return lines[0], "".join(l.strip() for l in lines[1:])
