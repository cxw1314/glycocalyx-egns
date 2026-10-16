#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lightweight GEO series-matrix probe (stdlib only).

Streams the *header* of a GSE series matrix .gz file over HTTP and prints
sample-level metadata (titles / characteristics) so that case-control
composition can be judged before downloading full expression matrices.

Usage:
    python geo_probe.py GSE43292 GSE100927 ...
"""
from __future__ import annotations

import re
import sys
import urllib.request
import zlib

HOST = "https://ftp.ncbi.nlm.nih.gov/geo/series"
UA = {"User-Agent": "glycocalyx_cvd-probe"}


def _open(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout)


def matrix_dir(gse: str) -> str:
    stem = gse[:-3] + "nnn" if len(gse) > 6 else gse
    return f"{HOST}/{stem}/{gse}/matrix/"


def find_matrix(gse: str) -> str | None:
    url = matrix_dir(gse)
    try:
        html = _open(url).read().decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        print(f"  !! cannot list {url}: {exc}")
        return None
    names = re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', html)
    if not names:
        return None
    # prefer the largest / canonical file
    return url + sorted(set(names))[0]


def probe(gse: str, max_bytes: int = 400_000) -> dict:
    murl = find_matrix(gse)
    if not murl:
        print(f"{gse}: no series matrix found")
        return {}
    print(f"{gse}: {murl}")
    try:
        with _open(murl) as resp:
            raw = resp.read(max_bytes)
        # decompress incrementally so a deliberately truncated stream is fine
        dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
        text = dec.decompress(raw).decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        print(f"  !! cannot read matrix header: {exc}")
        return {}
    info: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("!"):
            break
        key, _, val = line.partition("\t")
        info.setdefault(key.strip(), []).append(val.strip())
    return info


def summarize(gse: str, info: dict) -> None:
    if not info:
        return
    print("-" * 78)
    print(f"### {gse}")
    for k in ("!Series_title", "!Series_summary", "!Series_overall_design",
              "!Series_platform_id", "!Series_sample_id"):
        if k in info:
            v = " ".join(info[k])[:300]
            print(f"  {k[1:]}: {v}")
    titles = info.get("!Sample_title", [])
    print(f"  n_samples(header) = {len(titles)}")
    for k in sorted(info):
        if k.startswith("!Sample_characteristics"):
            vals = info[k]
            # show the distinct values of this characteristic field
            distinct = []
            seen = set()
            for v in vals:
                if v not in seen:
                    seen.add(v)
                    distinct.append(v)
            print(f"  {k[1:]} ({len(distinct)} distinct): {distinct[:14]}")
        elif k in ("!Sample_source_name_ch1", "!Sample_organism_ch1"):
            distinct = []
            seen = set()
            for v in info[k]:
                if v not in seen:
                    seen.add(v)
                    distinct.append(v)
            print(f"  {k[1:]} ({len(distinct)} distinct): {distinct[:14]}")
    print(f"  sample_titles[:12] = {titles[:12]}")


def main() -> int:
    for gse in sys.argv[1:]:
        gse = gse.strip().upper()
        try:
            info = probe(gse)
            summarize(gse, info)
        except Exception as exc:  # noqa: BLE001
            print(f"{gse}: ERROR {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
