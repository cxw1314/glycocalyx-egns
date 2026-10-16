#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build a single Visium AnnData for GSE176092 so that spatial deconvolution can
be run on it.

Design
------
The deconvolution step needs BOTH
  * raw counts (to re-normalise identically to the scRNA reference), and
  * the per-spot module scores / EGNS in .obs.

sc.tl.score_genes requires normalised data, so scoring is done on a copy and
only the resulting .obs columns are carried back; the saved .X keeps the raw
counts. obsm['spatial'] carries the Visium pixel coordinates.

Usage:
    python build_visium_h5ad.py --indir data/raw/GSE176092/visium \
                                --out results/GSE176092/visium_spatial.h5ad
"""
from __future__ import annotations

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analysis_visium_spatial import load_visium, load_modules, resolve_genes, condition_of  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/raw/GSE176092/visium")
    ap.add_argument("--out", default="results/GSE176092/visium_spatial.h5ad")
    ap.add_argument("--min-counts", type=int, default=500)
    args = ap.parse_args()

    import glob
    import scanpy as sc
    import anndata as ad

    def rp(p):
        return p if os.path.isabs(p) else os.path.join(ROOT, p)

    indir = rp(args.indir)
    out = rp(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)

    prefixes = sorted(p[:-len("_matrix.mtx.gz")] for p in
                      glob.glob(os.path.join(indir, "*_matrix.mtx.gz")))
    if not prefixes:
        print(f"[error] no Visium matrices under {indir}")
        return 1
    print(f"[load] {len(prefixes)} sections")

    modules = load_modules()
    parts = []

    for p in prefixes:
        name = os.path.basename(p)
        a = load_visium(p)
        if a is None:
            print(f"  !! {name}: load failed")
            continue
        sc.pp.filter_cells(a, min_counts=args.min_counts)
        if a.n_obs < 50:
            print(f"  !! {name}: {a.n_obs} spots after QC, skipped")
            continue

        # score on a normalised copy, keep counts in the object we save
        b = a.copy()
        sc.pp.normalize_total(b, target_sum=1e4)
        sc.pp.log1p(b)
        for key, spec in modules.items():
            genes = resolve_genes(spec, b.var_names)
            if genes:
                sc.tl.score_genes(b, genes, score_name=f"mod_{key}")

        def zc(col):
            s = b.obs[col]
            sd = s.std(ddof=0)
            return (s - s.mean()) / sd if sd else s * 0.0

        need = {"mod_GDS", "mod_EAS", "mod_EIDS", "mod_GCI", "mod_EIS"}
        if not need <= set(b.obs.columns):
            print(f"  !! {name}: missing module scores")
            continue
        b.obs["EGNS"] = (zc("mod_GDS") + zc("mod_EAS") + zc("mod_EIDS")
                         - zc("mod_GCI") - zc("mod_EIS"))

        a = a.copy()
        for c in b.obs.columns:
            if c.startswith("mod_") or c == "EGNS":
                a.obs[c] = b.obs[c].values
        a.obs["sample"] = name
        a.obs["condition"] = condition_of(name)
        a.obsm["spatial"] = np.asarray(a.obsm["spatial"], dtype=float)
        a.var_names_make_unique()
        parts.append(a)
        print(f"  ok {name:42s} {a.n_obs:5d} spots  cond={condition_of(name)}")

    if not parts:
        print("[error] nothing loaded")
        return 1

    sp = ad.concat(parts, join="outer", index_unique="|", fill_value=0)
    sp.obs_names_make_unique()
    print(f"\n[concat] {sp.n_obs} spots x {sp.n_vars} genes")
    print(sp.obs.groupby("condition")["sample"].nunique().to_string())
    sp.write_h5ad(out)
    print(f"[saved] {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
