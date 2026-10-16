#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parallel segmented downloader for throughput-limited links.

Single connections to NCBI are throttled to ~10-15 KB/s, but aggregate
bandwidth scales with the number of concurrent connections. Files are
therefore split into byte ranges that are fetched concurrently with curl
and concatenated, with per-range retry/resume and post-hoc integrity checks.

Usage:
    python pfetch.py --url <url> --out <path> [--workers 8] [--chunk-mb 4]
    python pfetch.py --manifest <file>          # lines of "url|dest"
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import os
import shutil
import subprocess
import sys
import threading
import time
from urllib.parse import urlparse

CURL = r"C:\Windows\system32\curl.exe"
PRINT_LOCK = threading.Lock()


def log(msg: str) -> None:
    with PRINT_LOCK:
        print(msg, flush=True)


def remote_size(url: str, attempts: int = 6) -> int | None:
    """Total object size in bytes, with a ranged-GET fallback if HEAD is blocked."""
    for _ in range(attempts):
        try:
            out = subprocess.run(
                [CURL, "-sIL", "--ssl-no-revoke", "--max-time", "90", url],
                capture_output=True, text=True, timeout=120,
            ).stdout
            sizes = []
            for line in out.splitlines():
                if line.lower().startswith("content-length:"):
                    v = line.split(":", 1)[1].strip()
                    if v.isdigit():
                        sizes.append(int(v))
            if sizes:
                return sizes[-1]
        except Exception:  # noqa: BLE001
            pass
        try:
            out = subprocess.run(
                [CURL, "-s", "-D", "-", "-o", "NUL", "--ssl-no-revoke",
                 "--max-time", "90", "-r", "0-0", url],
                capture_output=True, text=True, timeout=120,
            ).stdout
            for line in out.splitlines():
                if line.lower().startswith("content-range"):
                    tot = line.split("/")[-1].strip()
                    if tot.isdigit():
                        return int(tot)
        except Exception:  # noqa: BLE001
            pass
        time.sleep(3)
    return None


def fetch_range(url: str, start: int, end: int, part: str, attempts: int = 80) -> bool:
    """Fetch bytes [start, end] into `part`, resuming by appending.

    Only payload from a successful HTTP 206 (or a 200 that genuinely covers the
    requested range) is appended, so rate-limit / error pages can never corrupt
    a chunk.
    """
    want = end - start + 1
    tmp = part + ".tmp"
    for _ in range(attempts):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have >= want:
            if have > want:                      # trim overshoot
                with open(part, "r+b") as fh:
                    fh.truncate(want)
            return True
        cur = start + have
        rng = f"{cur}-{end}"
        try:
            res = subprocess.run(
                [CURL, "-s", "--ssl-no-revoke", "--max-time", "300",
                 "--connect-timeout", "20", "--speed-limit", "200",
                 "--speed-time", "45", "-r", rng, "-o", tmp,
                 "-w", "%{http_code}", url],
                capture_output=True, text=True, timeout=360,
            )
            code = (res.stdout or "").strip()
        except Exception:  # noqa: BLE001
            code = ""

        ok = (code == "206")
        if code == "200":
            got0 = os.path.getsize(tmp) if os.path.exists(tmp) else 0
            ok = (cur == start and got0 == want)
        if ok and os.path.exists(tmp) and os.path.getsize(tmp) > 0:
            with open(part, "ab") as dst, open(tmp, "rb") as src:
                shutil.copyfileobj(src, dst, 1 << 22)
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        if not ok:
            time.sleep(1.5)
    have = os.path.getsize(part) if os.path.exists(part) else 0
    return have >= want


def verify(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) < 128:
        return False
    try:
        if path.endswith((".gz", ".tgz")):
            with gzip.open(path, "rb") as fh:
                while fh.read(1 << 22):
                    pass
        else:
            with open(path, "rb") as fh:
                fh.read(4096)
        return True
    except Exception:  # noqa: BLE001
        return False


def fetch_file(url: str, dest: str, workers: int, chunk_mb: int) -> bool:
    if verify(dest):
        log(f"[skip] verified: {dest}")
        return True
    total = remote_size(url)
    if not total:
        log(f"[error] cannot determine size: {url}")
        return False
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    chunk = max(chunk_mb * 1024 * 1024, total // (workers * 2) + 1)
    ranges = []
    s = 0
    while s < total:
        e = min(s + chunk - 1, total - 1)
        ranges.append((s, e))
        s = e + 1
    log(f"[get] {os.path.basename(dest)}  {total/1e6:.1f} MB -> {len(ranges)} ranges x {workers} workers")

    tmpdir = dest + ".parts"
    os.makedirs(tmpdir, exist_ok=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_range, url, a, b,
                          os.path.join(tmpdir, f"{a:012d}-{b:012d}.bin")): (a, b)
                for a, b in ranges}
        for n, fut in enumerate(cf.as_completed(futs), 1):
            a, b = futs[fut]
            try:
                good = fut.result()
            except Exception as exc:  # noqa: BLE001
                log(f"    !! range {a}-{b} failed: {exc}")
                good = False
            if not good:
                ok = False
            if n % 5 == 0 or n == len(ranges):
                log(f"    ranges done: {n}/{len(ranges)}")
    if not ok:
        log(f"[INCOMPLETE] {dest} (parts kept for resume)")
        return False

    with open(dest, "wb") as out:
        for a, b in ranges:
            with open(os.path.join(tmpdir, f"{a:012d}-{b:012d}.bin"), "rb") as fh:
                while True:
                    buf = fh.read(1 << 22)
                    if not buf:
                        break
                    out.write(buf)
    got = os.path.getsize(dest)
    if got != total:
        log(f"[ERROR] size mismatch {dest}: {got} != {total}")
        return False
    if not verify(dest):
        log(f"[ERROR] integrity check failed: {dest}")
        return False
    for a, b in ranges:
        try:
            os.remove(os.path.join(tmpdir, f"{a:012d}-{b:012d}.bin"))
        except OSError:
            pass
    try:
        os.rmdir(tmpdir)
    except OSError:
        pass
    log(f"[ok] {dest} ({got/1e6:.1f} MB)")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url")
    ap.add_argument("--out")
    ap.add_argument("--manifest")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chunk-mb", type=int, default=3)
    args = ap.parse_args()

    items = []
    if args.manifest:
        with open(args.manifest, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                u, _, d = line.partition("|")
                items.append((u.strip(), d.strip().replace("/", os.sep)))
    elif args.url and args.out:
        items.append((args.url, args.out))
    else:
        ap.error("need --manifest or both --url and --out")

    all_ok = True
    for url, dest in items:
        try:
            all_ok &= fetch_file(url, dest, args.workers, args.chunk_mb)
        except Exception as exc:  # noqa: BLE001
            log(f"[FAIL] {url} -> {exc}")
            all_ok = False
    log("[fetch done]" if all_ok else "[fetch finished with failures]")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
