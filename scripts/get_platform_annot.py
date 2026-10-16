#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extract a probe -> gene-symbol map from the head of a GEO GPL SOFT family file.

GPL17586 (Affymetrix HTA 2.0) family.soft.gz is ~7.6 GB, but the ^PLATFORM
section -- which contains the full platform table with gene_assignment -- sits
at the very start of the file. We therefore download only a growing prefix of
the gzip stream (in parallel byte ranges, since the link is per-connection
throttled), decompress it incrementally, and stop as soon as the
!platform_table_end marker appears.

Output: <outdir>/GPL17586_annot.tsv  (probe_id<TAB>gene_symbol)

Usage:
    python get_platform_annot.py [--outdir data/raw/GPL17586] [--step-mb 40] [--max-mb 400]
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pfetch  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17586/soft/"
       "GPL17586_family.soft.gz")


def download_prefix(url: str, nbytes: int, workdir: str, workers: int) -> str:
    """Fetch bytes [0, nbytes) into workdir/prefix.bin using cached ranges."""
    os.makedirs(workdir, exist_ok=True)
    chunk = 2 * 1024 * 1024
    ranges = [(s, min(s + chunk - 1, nbytes - 1))
              for s in range(0, nbytes, chunk)]
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(pfetch.fetch_range, url, a, b,
                          os.path.join(workdir, f"{a:012d}-{b:012d}.bin")): (a, b)
                for a, b in ranges}
        for fut in cf.as_completed(futs):
            fut.result()
    return os.path.join(workdir, "prefix.bin"), ranges


def assemble(workdir: str, ranges: list[tuple[int, int]]) -> str:
    out = os.path.join(workdir, "prefix.bin")
    with open(out, "wb") as fh:
        for a, b in ranges:
            p = os.path.join(workdir, f"{a:012d}-{b:012d}.bin")
            if not os.path.exists(p):
                continue
            with open(p, "rb") as src:
                while True:
                    buf = src.read(1 << 22)
                    if not buf:
                        break
                    fh.write(buf)
    return out


def extract_platform_table(path: str) -> tuple[dict[str, str], bool]:
    """Stream-decompress and collect the ^PLATFORM table rows.

    Returns (rows, complete) where `complete` is True only when the
    !platform_table_end marker was actually reached. Note that many HTA 2.0
    transcript clusters carry no gene assignment, so the number of annotated
    rows is well below the platform's row count -- completion must therefore
    be judged by the end marker, not by a row threshold.
    """
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    rows: dict[str, str] = {}
    header: list[str] | None = None
    id_col = sym_col = None
    in_table = False
    tail = b""
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(1 << 22)
            if not buf:
                break
            try:
                text = dec.decompress(buf)
            except zlib.error:
                break
            if not text:
                continue
            data = tail + text
            lines = data.split(b"\n")
            tail = lines.pop()
            for raw in lines:
                line = raw.decode("utf-8", "replace").rstrip("\r")
                if not in_table:
                    if line.startswith("!platform_table_begin"):
                        in_table = True
                    continue
                if line.startswith("!platform_table_end"):
                    return rows, True
                if header is None:
                    header = line.split("\t")
                    id_col = header.index("ID") if "ID" in header else 0
                    sym_col = None
                    for cand in ("gene_assignment", "Gene Symbol", "Symbol"):
                        if cand in header:
                            sym_col = header.index(cand)
                            break
                    continue
                parts = line.split("\t")
                if sym_col is None or len(parts) <= max(id_col, sym_col):
                    continue
                sym = _first_symbol(parts[sym_col])
                if sym:
                    rows[parts[id_col]] = sym
    return rows, False


def _first_symbol(raw: str) -> str:
    import re
    if not raw or raw in ("---", "null", ""):
        return ""
    for chunk in raw.split("///"):
        fields = [f.strip() for f in chunk.split("//")]
        if len(fields) >= 2 and fields[1] and fields[1] != "---":
            cand = fields[1].split("//")[0].strip()
            if re.fullmatch(r"[A-Za-z0-9._\-@]{1,40}", cand):
                return cand.upper()
    return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "raw", "GPL17586"))
    ap.add_argument("--step-mb", type=int, default=40)
    ap.add_argument("--max-mb", type=int, default=400)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    workdir = os.path.join(args.outdir, "soft_prefix")
    limit = args.step_mb * 1024 * 1024
    maxb = args.max_mb * 1024 * 1024
    rows: dict[str, str] = {}
    complete = False

    while limit <= maxb:
        print(f"[prefix] downloading first {limit/1e6:.0f} MB ...", flush=True)
        _, ranges = download_prefix(URL, limit, workdir, args.workers)
        pre = assemble(workdir, ranges)
        print(f"[prefix] assembled {os.path.getsize(pre)/1e6:.1f} MB, extracting ...",
              flush=True)
        rows, complete = extract_platform_table(pre)
        print(f"[prefix] extracted {len(rows)} annotated probe rows "
              f"(table_complete={complete})", flush=True)
        if complete:            # reached !platform_table_end -> we have the whole table
            break
        limit += args.step_mb * 1024 * 1024

    if not rows:
        print("[error] no annotation rows extracted")
        return 1
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, "GPL17586_annot.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        for k, v in rows.items():
            fh.write(f"{k}\t{v}\n")
    print(f"[ok] {len(rows)} rows -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
