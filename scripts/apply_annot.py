#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Collapse a probe-level expression matrix (probes x samples) to gene symbols
using a probe<TAB>symbol map, keeping the highest-mean probe per symbol.

Usage:
    python apply_annot.py --expr <expr.pkl> --annot <annot.tsv> --out <expr_genes.pkl> [--backup <path>]
"""
from __future__ import annotations

import argparse
import pickle
import shutil

import numpy as np
import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--backup", default=None)
    args = ap.parse_args()

    with open(args.expr, "rb") as fh:
        expr = pickle.load(fh)
    print(f"[expr] {expr.shape} probes x samples")

    mp: dict[str, str] = {}
    with open(args.annot, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                mp[parts[0]] = parts[1]
    print(f"[annot] {len(mp)} probe->symbol rows")

    sym = pd.Series({i: mp.get(i, "") for i in expr.index})
    keep = sym[sym != ""]
    matched = int(keep.shape[0])
    sub = expr.loc[keep.index]
    sub = sub.assign(__sym=keep.values, __mean=sub.mean(axis=1).values)
    sub = sub.sort_values("__mean", ascending=False).drop_duplicates("__sym")
    sub = sub.drop(columns=["__mean"]).set_index("__sym")
    sub.index.name = "gene"
    print(f"[genes] matched probes={matched}, unique symbols={sub.shape[0]}")

    if args.backup:
        shutil.copy(args.expr, args.backup)
    with open(args.out, "wb") as fh:
        pickle.dump(sub.astype(np.float32), fh)
    print(f"[ok] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
