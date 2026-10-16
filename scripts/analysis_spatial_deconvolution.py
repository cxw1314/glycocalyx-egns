#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Spatial deconvolution of Visium spots before re-testing EGNS spatial coupling.

Motivation (reviewer-proofing): an EGNS score averaged over a Visium spot
reflects both the true endothelial state AND how many endothelial cells that
spot contains. If high-EGNS spots are simply endothelial-rich spots, the
"EGNS-fibrotic niche adjacency" result is confounded by composition.

This script:
  1. builds a cell-type signature matrix from the matched scRNA reference
     (GSE176092 non-cardiomyocyte, lineage labels);
  2. deconvolves each Visium spot by non-negative least squares (NNLS),
     the same regression family RCTD/cell2location use, but fully offline;
  3. recomputes the EGNS-fibrotic-niche neighbourhood correlation with the
     endothelial *fraction* regressed out, versus the raw correlation.

Output:
    results/spatial_deconvolution/tables/spot_celltype_fractions.csv
    results/spatial_deconvolution/tables/coupling_adjusted.csv
"""
from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def z(v):
    sd = np.nanstd(v)
    return (v - np.nanmean(v)) / sd if sd else v * 0.0


def nnls_deconvolve(sig: np.ndarray, obs: np.ndarray) -> np.ndarray:
    """Column-wise NNLS with an intercept, renormalised to fractions.

    The intercept is implemented by augmenting BOTH the design matrix and the
    observation vector with a unit entry (appending only to the design matrix
    leaves the dimensions inconsistent).
    """
    from scipy.optimize import nnls
    n_spots = obs.shape[1]
    W = np.zeros((sig.shape[1], n_spots), dtype=np.float32)
    A = np.vstack([sig, np.ones((1, sig.shape[1]), dtype=sig.dtype)])
    for i in range(n_spots):
        y = np.concatenate([obs[:, i], np.ones(1, dtype=obs.dtype)])
        w, _ = nnls(A, y)
        W[:, i] = w
    s = W.sum(axis=0, keepdims=True)
    s[s == 0] = 1.0
    return W / s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scrna", required=True, help="scRNA reference h5ad")
    ap.add_argument("--spatial", required=True, help="Visium h5ad")
    ap.add_argument("--outdir", default="results/spatial_deconvolution")
    ap.add_argument("--knn", type=int, default=6)
    ap.add_argument("--n-genes", type=int, default=3000,
                    help="unused legacy knob (kept for CLI compatibility)")
    ap.add_argument("--n-markers", type=int, default=150,
                    help="top-N most cell-type-specific genes per lineage")
    args = ap.parse_args()

    import anndata as ad

    def rp(p):
        return p if os.path.isabs(p) else os.path.join(ROOT, p)

    outdir = rp(args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    sc_ref = ad.read_h5ad(rp(args.scrna))
    sp = ad.read_h5ad(rp(args.spatial))
    print(f"[scRNA] {sc_ref.n_obs} x {sc_ref.n_vars}   [spatial] {sp.n_obs} x {sp.n_vars}")

    # shared gene space (upper-cased)
    sc_raw = sc_ref.raw if sc_ref.raw is not None else sc_ref
    sp_raw = sp.raw if sp.raw is not None else sp
    sc_lut = {v.upper(): v for v in sc_raw.var_names}
    sp_lut = {v.upper(): v for v in sp_raw.var_names}
    shared = sorted(set(sc_lut) & set(sp_lut))
    print(f"[shared genes] {len(shared)}")

    def cpm(X):
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        X = np.log1p(X / np.maximum(X.sum(axis=1, keepdims=True), 1e-9) * 1e4)
        return X

    # Signature genes must DISCRIMINATE cell types. Using the most variable
    # genes instead lets one highly-abundant, highly-correlated lineage (SMC)
    # absorb essentially the whole spot and collapses the fractions; the
    # standard fix is a specificity-ranked marker panel per lineage.
    lab = sc_ref.obs["lineage"].astype(str)
    keep = ~lab.str.contains("Cardiomyocyte", case=False, na=False)
    sc_X = cpm(sc_raw[keep.values][:, [sc_lut[g] for g in shared]].X)
    lab = lab[keep.values]

    types = [ct for ct in sorted(lab.unique()) if (lab == ct).sum() >= 20]
    means = np.vstack([sc_X[(lab == ct).values].mean(axis=0) for ct in types])
    spec = np.zeros_like(means)
    for i in range(len(types)):
        others = np.delete(means, i, axis=0).max(axis=0)
        spec[i] = means[i] - others          # specificity vs the best other type
    n_mark = int(args.n_markers)
    sel: set[int] = set()
    for i in range(len(types)):
        for t in np.argsort(spec[i])[::-1][:n_mark]:
            sel.add(int(t))
    sel = sorted(sel)
    names = types
    sig = means[:, sel].T                     # genes x cell types
    print(f"[signature] {sig.shape[0]} marker genes x {sig.shape[1]} cell types: "
          f"{', '.join(names)}")
    sel_genes = [shared[i] for i in sel]

    sp_X = cpm(sp_raw[:, [sp_lut[g] for g in sel_genes]].X).T
    frac = nnls_deconvolve(sig.astype(np.float32), sp_X.astype(np.float32)).T
    fr = pd.DataFrame(frac, index=sp.obs_names, columns=names)
    fr.to_csv(os.path.join(tabdir, "spot_celltype_fractions.csv"),
              encoding="utf-8-sig")
    print("\n=== mean deconvolved fraction per condition ===")
    cond = sp.obs["condition"] if "condition" in sp.obs else pd.Series("all", index=sp.obs_names)
    print(fr.groupby(cond.values).mean().round(4).to_string())

    # ---- spatial neighbourhood, with and without EC fraction regressed out --
    coords = sp.obsm.get("spatial", None)
    if coords is None:
        print("[warn] no obsm['spatial']; writing fractions only")
        return 0
    coords = np.asarray(coords)

    from sklearn.neighbors import NearestNeighbors
    nb = NearestNeighbors(n_neighbors=args.knn + 1).fit(coords)
    _, ind = nb.kneighbors(coords)

    def neigh_mean(v):
        return np.array([v[ind[i][1:]].mean() for i in range(len(v))])

    eg = sp.obs["EGNS"].values if "EGNS" in sp.obs else None
    if eg is None:
        print("[warn] no obs['EGNS']; writing fractions only")
        return 0

    ec_col = [c for c in names if "Endothelial" in c]
    rows = []
    for target in [c for c in names if any(k in c for k in
                   ("Fibroblast", "Macrophage"))]:
        y_raw = neigh_mean(fr[target].values)
        r_raw = float(np.corrcoef(eg, y_raw)[0, 1])
        if ec_col:
            ec = fr[ec_col[0]].values
            # residualise both EGNS and the neighbour signal on EC fraction
            def resid(v, c):
                A = np.vstack([c, np.ones_like(c)]).T
                beta, *_ = np.linalg.lstsq(A, v, rcond=None)
                return v - A @ beta
            r_adj = float(np.corrcoef(resid(eg, ec), resid(y_raw, ec))[0, 1])
        else:
            r_adj = np.nan
        rows.append({"neighbour_programme": target,
                     "rho_EGNS_vs_neighbour_raw": round(r_raw, 4),
                     "rho_after_regressing_out_EC_fraction": round(r_adj, 4),
                     "delta": round(r_adj - r_raw, 4)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(tabdir, "coupling_adjusted.csv"),
               index=False, encoding="utf-8-sig")
    print("\n=== EGNS-fibrotic/immune niche coupling, composition-adjusted ===")
    print(out.to_string(index=False))

    # ---- also adjust the MODULE-SCORE coupling (comparable to main text) --
    if "EGNS" in sp.obs and {"mod_FNS", "mod_MFS"} <= set(sp.obs.columns):
        rows2 = []
        ec = fr[ec_col[0]].values if ec_col else None
        for mod in ("mod_FNS", "mod_MFS", "mod_GDS"):
            y_raw = neigh_mean(sp.obs[mod].values.astype(float))
            r_raw = float(np.corrcoef(eg, y_raw)[0, 1])
            if ec is not None:
                def resid(v, c):
                    A = np.vstack([c, np.ones_like(c)]).T
                    b, *_ = np.linalg.lstsq(A, v, rcond=None)
                    return v - A @ b
                r_adj = float(np.corrcoef(resid(eg, ec), resid(y_raw, ec))[0, 1])
            else:
                r_adj = np.nan
            rows2.append({"neighbour_module": mod,
                          "rho_EGNS_vs_neighbour_raw": round(r_raw, 4),
                          "rho_after_regressing_out_EC_fraction": round(r_adj, 4),
                          "delta": round(r_adj - r_raw, 4)})
        out2 = pd.DataFrame(rows2)
        out2.to_csv(os.path.join(tabdir, "coupling_adjusted_module_scores.csv"),
                    index=False, encoding="utf-8-sig")
        print("\n=== EGNS vs neighbourhood MODULE scores, composition-adjusted ===")
        print(out2.to_string(index=False))

    print("\n[ok] spatial deconvolution written to " + tabdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
