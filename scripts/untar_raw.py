#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract the 10x .h5 matrices (and optionally Visium image tars) from a GEO RAW.tar."""
from __future__ import annotations

import argparse
import os
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tar", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--include-images", action="store_true")
    ap.add_argument("--pattern", default=None,
                    help="only extract members whose basename contains this substring")
    ap.add_argument("--exclude", default=None,
                    help="skip members whose basename contains this substring")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    n = 0
    with tarfile.open(args.tar) as t:
        for m in t.getmembers():
            if not m.isfile():
                continue
            name = os.path.basename(m.name)
            if args.pattern and args.pattern not in name:
                continue
            if args.exclude and args.exclude in name:
                continue
            if name.endswith(".h5"):
                pass
            elif name.endswith((".mtx.gz", ".tsv.gz")):
                pass
            elif args.include_images and (name.endswith(".tar.gz") or name.endswith(".csv.gz")):
                pass
            else:
                print(f"[skip] {name}")
                continue
            dst = os.path.join(args.outdir, name)
            if os.path.exists(dst) and os.path.getsize(dst) == m.size:
                print(f"[have] {name}")
                continue
            with t.extractfile(m) as src, open(dst, "wb") as out:
                while True:
                    buf = src.read(1 << 22)
                    if not buf:
                        break
                    out.write(buf)
            print(f"[ok]   {name} ({m.size/1e6:.1f} MB)")
            n += 1
    print(f"[done] extracted {n} files -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
