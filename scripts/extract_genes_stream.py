#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Streaming gene extractor for very large genes x cells CSV.GZ matrices.

GSE183852_Integrated_Counts.csv.gz is a 487 MB gzip whose header alone is
7.5 MB (hundreds of thousands of cell barcodes); the fully expanded matrix is
far too large to load. We therefore stream it once, keep only the requested
gene rows, and write a compact cells x genes table.

Cell identity is taken from the barcode prefix (e.g. "H_ZC-11-292_ACGT..." ->
sample "H_ZC-11-292"), which can then be mapped to donor vs DCM via the GEO
series matrix.

Usage:
    python extract_genes_stream.py --path <csv.gz> --genes A,B,C --out <dir> [--limit N]
"""
from __future__ import annotations

import argparse
import gzip
import os
import pickle
import time

import numpy as np
import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", required=True)
    ap.add_argument("--genes", default="", help="comma-separated symbols")
    ap.add_argument("--genes-file", default=None,
                    help="optional file with one symbol per line (merged with --genes)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0, help="stop after N gene rows (for timing)")
    args = ap.parse_args()

    wanted = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
    if args.genes_file and os.path.exists(args.genes_file):
        with open(args.genes_file, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    wanted.add(line.upper())
    print(f"[genes] {len(wanted)} requested")

    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    rows: list[tuple[str, np.ndarray]] = []
    cells: list[str] | None = None
    n_rows = 0
    n_kept = 0

    with gzip.open(args.path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        parts = header.split(",")
        cells = parts[1:]
        print(f"[header] {len(cells)} cells (first: {cells[0] if cells else 'n/a'})")
        for line in fh:
            n_rows += 1
            if args.limit and n_rows > args.limit:
                break
            if not line:
                continue
            sym, sep, rest = line.partition(",")
            if not sep:
                continue
            if sym.strip().upper() not in wanted:
                if n_rows % 2000 == 0:
                    print(f"  scanned {n_rows} rows, kept {n_kept}, "
                          f"{time.time()-t0:.0f}s", flush=True)
                continue
            vals = np.array(rest.rstrip("\n").split(","), dtype=np.float32)
            rows.append((sym.strip().upper(), vals))
            n_kept += 1
            print(f"  [keep] {sym}  ({n_kept}/{len(wanted)})  {time.time()-t0:.0f}s",
                  flush=True)
            if n_kept >= len(wanted):
                break
        print(f"[scan] rows scanned={n_rows}, kept={n_kept}, "
              f"elapsed={time.time()-t0:.0f}s")

    if not rows:
        print("[error] no requested gene found")
        return 1
    genes = [r[0] for r in rows]
    mat = np.vstack([r[1] for r in rows])
    df = pd.DataFrame(mat, index=genes, columns=cells)
    df.index.name = "gene"
    print(f"[matrix] {df.shape} (genes x cells)")

    with open(os.path.join(args.out, "counts_selected.pkl"), "wb") as fh:
        pickle.dump(df, fh)
    pd.Series(cells, name="barcode").to_csv(
        os.path.join(args.out, "barcodes.csv"), index=False, encoding="utf-8-sig")
    print(f"[ok] written to {args.out}")
    print(f"[missing] {len(wanted) - n_kept} requested genes not found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
