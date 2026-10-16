#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Robust resumable downloader (stdlib only).

python download.py <url> <dest> [<url> <dest> ...]
"""
from __future__ import annotations

import os
import sys
import urllib.request

UA = {"User-Agent": "glycocalyx_cvd-dl"}


def size_of(url: str) -> int | None:
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:  # noqa: BLE001
        return None


def download(url: str, dest: str) -> bool:
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    if os.path.exists(dest):
        exp = size_of(url)
        if exp is not None and os.path.getsize(dest) == exp:
            print(f"[skip] already complete: {dest}")
            return True
    part = dest + ".part"
    exp = size_of(url)
    pos = os.path.getsize(part) if os.path.exists(part) else 0
    req = urllib.request.Request(url, headers=UA)
    if pos:
        req.add_header("Range", f"bytes={pos}-")
    with urllib.request.urlopen(req, timeout=180) as r:
        if pos and r.status != 206:      # server ignored Range -> restart
            print(f"[warn] no range support, restarting {dest}")
            pos = 0
        mode = "ab" if pos else "wb"
        total = pos + int(r.headers.get("Content-Length") or 0)
        done = pos
        with open(part, mode) as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                done += len(chunk)
                if total:
                    pct = 100.0 * done / total
                    print(f"\r  {os.path.basename(dest)}: {done/1e6:.1f}/{total/1e6:.1f} MB ({pct:.1f}%)",
                          end="", flush=True)
    print()
    got = os.path.getsize(part)
    if exp is not None and got != exp:
        print(f"[ERROR] size mismatch {dest}: {got} != {exp}")
        return False
    os.replace(part, dest)
    print(f"[ok] {dest} ({got/1e6:.1f} MB)")
    return True


def main() -> int:
    args = sys.argv[1:]
    if len(args) % 2:
        print("usage: download.py <url> <dest> [<url> <dest> ...]")
        return 2
    ok = True
    for i in range(0, len(args), 2):
        try:
            ok &= download(args[i], args[i + 1])
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL] {args[i]} -> {exc}")
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
