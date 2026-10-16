#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
List GEO matrix / supplementary files (name + size) for a GSE, so download
feasibility can be judged before pulling multi-GB objects.

Usage:
    python geo_files.py GSE123342 GSE206787 ...
"""
from __future__ import annotations

import re
import sys
import urllib.request

HOST = "https://ftp.ncbi.nlm.nih.gov/geo/series"
UA = {"User-Agent": "glycocalyx_cvd-files"}

ROW = re.compile(
    r'href="(?P<href>[^"]+)"[^>]*>[^<]*</a>\s*(?P<date>\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2})?\s*(?P<size>[0-9.]+[KMGT]?)?',
)


def stem(gse: str) -> str:
    return gse[:-3] + "nnn" if len(gse) > 6 else gse


def list_dir(url: str) -> list[tuple[str, str]]:
    try:
        req = urllib.request.Request(url, headers=UA)
        html = urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        print(f"  !! {url} -> {exc}")
        return []
    out = []
    for m in ROW.finditer(html):
        name = m.group("href")
        if name in ("../", "./") or name.startswith("/") or name.startswith("http"):
            continue
        out.append((name, m.group("size") or "?"))
    return out


def main() -> int:
    for gse in sys.argv[1:]:
        gse = gse.strip().upper()
        base = f"{HOST}/{stem(gse)}/{gse}"
        print("=" * 78)
        print(f"### {gse}")
        for sub in ("matrix", "suppl"):
            print(f"--- {sub}/")
            for name, size in list_dir(f"{base}/{sub}/"):
                print(f"    {size:>8}  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
