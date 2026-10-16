#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
In silico (virtual) knockout validation of the EGNS framework.

Design notes
------------
* CellOracle is the canonical tool, but it needs a base GRN fetched from GitHub
  (unreachable on this network) and benefits from scATAC motifs we do not have.
  We therefore implement the same two core steps offline:
    1. kNN imputation to overcome scRNA-seq dropout (as CellOracle does), and
    2. a ridge-regression GRN, perturbed CellOracle-style:
           shift_target = -Z[cell, ko_tf] * B[target, ko_tf]
           Z_knockout   = Z + shift
  after which the perturbed expression is pushed back through the EGNS modules.
* Regulators are restricted to the TRANSCRIPTION FACTOR panel (rather than all
  156 panel genes). Including every co-expressed gene as a regulator makes the
  design matrix collinear and ridge shrinks each per-edge coefficient to ~1e-3,
  which is uninterpretable; a TF -> target design is the biologically correct
  formulation and is well posed.
* The null distribution is built by shuffling the rows of the coefficient
  matrix, which preserves the marginal edge-weight distribution while
  destroying target-specific wiring.

Guard-rails: the GRN is learned from the same observational data, the model is
linear and has no chromatin layer, so results are prioritisation predictions,
not causal proof.

Usage:
    python analysis_virtual_ko.py --h5ad results/GSE176092/gse176092_nonCM_processed.h5ad \
                                  --outdir results/virtual_ko
"""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")
SEED = 42

TF_PANEL = [
    "KLF2", "KLF4", "RELA", "NFKB1", "NFKB2", "RELB", "ERG", "FLI1",
    "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
    "STAT3", "HIF1A", "SMAD1", "SMAD3", "JUN", "JUNB", "FOS", "FOSL2",
    "SOX17", "SOX18", "PROX1", "EPAS1", "PPARG", "MEF2C", "GATA6", "NFATC1",
]


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        out[key] = [g.upper() for g in val.get("genes", [])]
    return out


def zscore(M: np.ndarray) -> np.ndarray:
    mu = M.mean(axis=0, keepdims=True)
    sd = M.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (M - mu) / sd


def eg_ns(ms: dict) -> np.ndarray:
    """EGNS = z(GDS) + z(EAS) + z(EIDS) - z(GCI) - z(EIS) across cells."""
    def z(v):
        sd = v.std()
        return (v - v.mean()) / sd if sd else v * 0.0
    return (z(ms["GDS"]) + z(ms["EAS"]) + z(ms["EIDS"])
            - z(ms["GCI"]) - z(ms["EIS"]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-null", type=int, default=300)
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--knn", type=int, default=30)
    args = ap.parse_args()

    import anndata as ad
    from sklearn.linear_model import Ridge
    from sklearn.neighbors import NearestNeighbors

    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    lut: dict[str, str] = {}
    for v in raw.var_names:
        lut.setdefault(v.upper(), v)
    print(f"[load] {adata.n_obs} cells x {raw.n_vars} genes")

    modules = load_modules()
    panel = sorted({g for gs in modules.values() for g in gs if g in lut})
    tfs = [t for t in TF_PANEL if t in lut]
    print(f"[panel] {len(panel)} module genes; {len(tfs)} TFs present: {tfs}")
    if len(tfs) < 5:
        print("[error] too few TFs present")
        return 1

    def dense(genes: list[str]) -> np.ndarray:
        M = raw[:, [lut[g] for g in genes]].X
        M = np.asarray(M.todense()) if hasattr(M, "todense") else np.asarray(M)
        rs = M.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        return np.log1p(M / rs * 1e4)

    Xp = dense(panel)          # targets
    Xt = dense(tfs)            # regulators
    pidx = {g: i for i, g in enumerate(panel)}
    tidx = {g: i for i, g in enumerate(tfs)}
    module_idx = {k: [pidx[g] for g in gs if g in pidx] for k, gs in modules.items()}

    rng = np.random.default_rng(SEED)
    rows, edges = [], []

    for celltype in ("Endothelial", "Macrophage"):
        mask = adata.obs["lineage"].astype(str).str.contains(celltype).values
        n = int(mask.sum())
        if n < 100:
            print(f"[skip] {celltype}: {n} cells")
            continue

        emb, used = None, "panel z-score"
        for key in ("X_pca", "X_umap"):
            if key in adata.obsm:
                emb, used = np.asarray(adata.obsm[key])[mask], key
                break
        if emb is None or emb.shape[1] == 0:
            emb = zscore(Xp[mask])
        k = min(args.knn, n - 1)
        nn = NearestNeighbors(n_neighbors=k).fit(emb)
        _, nb = nn.kneighbors(emb)

        # kNN-imputed log-expression, deliberately NOT z-scored: the knockout
        # readout is a mean over cells, and a z-scored perturbation is
        # mean-zero by construction, which would force every delta to ~0.
        P = Xp[mask][nb].mean(axis=1)      # targets
        T = Xt[mask][nb].mean(axis=1)      # TFs
        print(f"\n=== {celltype}: {n} cells (kNN k={k} in {used}) ===")

        B = np.zeros((len(panel), len(tfs)))          # target x TF
        r2 = []
        for t in range(len(panel)):
            y = P[:, t]
            if y.std() == 0:
                continue
            m = Ridge(alpha=args.alpha, fit_intercept=True).fit(T, y)
            B[t] = m.coef_
            yc = y - y.mean()
            pc = m.predict(T)
            pc = pc - pc.mean()
            r2.append(float(1 - ((yc - pc) ** 2).sum() / max((yc ** 2).sum(), 1e-9)))
        print(f"[GRN] {len(panel)} targets x {len(tfs)} TFs; "
              f"mean |coef|={np.abs(B).mean():.4f}, max |coef|={np.abs(B).max():.3f}, "
              f"mean centred R2={np.mean(r2):.3f}")

        for probe in ("ICAM1", "SPP1", "POSTN", "PECAM1"):
            if probe not in pidx:
                continue
            r = B[pidx[probe]]
            top = np.argsort(-np.abs(r))[:3]
            print(f"        top TFs for {probe}: "
                  + ", ".join(f"{tfs[i]}({r[i]:+.3f})" for i in top))
            for i in top:
                edges.append({"cell_type": celltype, "target": probe,
                              "tf": tfs[i], "coef": round(float(r[i]), 4)})

        def scores(M: np.ndarray) -> dict:
            return {kk: M[:, ix].mean(axis=1) for kk, ix in module_idx.items() if ix}

        base = scores(P)
        # EGNS is a z-score combination and therefore has mean 0 by construction.
        # To make a knockout readout meaningful we anchor the z-scoring to the
        # BASE population and apply those fixed parameters to the perturbed
        # scores, so that a shift in the distribution becomes visible.
        EG_PARTS = ("GDS", "EAS", "EIDS", "GCI", "EIS")
        if celltype == "Endothelial":
            zp = {k: (float(base[k].mean()), float(base[k].std() or 1.0))
                  for k in EG_PARTS}

            def egns(ms: dict) -> np.ndarray:
                z = lambda k: (ms[k] - zp[k][0]) / zp[k][1]   # noqa: E731
                return (z("GDS") + z("EAS") + z("EIDS") - z("GCI") - z("EIS"))

            base_read = float(egns(base).mean())

            def read(mss: dict) -> float:
                return float(egns(mss).mean())
        else:
            base_read = float(base["MFS"].mean())

            def read(mss: dict) -> float:
                return float(mss["MFS"].mean())

        # ---- null: row-shuffled coefficient matrix ----------------------
        null = []
        for _ in range(args.n_null):
            Bs = B[rng.permutation(B.shape[0]), :]
            for j in range(len(tfs)):
                Pk = P - np.outer(T[:, j], Bs[:, j])
                null.append(read(scores(Pk)) - base_read)
        null = np.asarray(null)
        print(f"[null] {len(null)} shuffled-edge perturbations; "
              f"sd={null.std():.4f}")

        for g in tfs:
            j = tidx[g]
            Pk = P - np.outer(T[:, j], B[:, j])
            ms = scores(Pk)
            delta = read(ms) - base_read
            p = float((np.sum(np.abs(null) >= abs(delta)) + 1) / (len(null) + 1))
            row = {"cell_type": celltype, "ko_tf": g,
                   "delta_readout": round(delta, 4),
                   "z_vs_null": round(float(delta / (null.std() or 1)), 2),
                   "p_empirical": round(p, 4),
                   "n_cells": n}
            for kk in [k for k in ("GDS", "EAS", "EIDS", "GCI", "EIS", "MFS")
                       if k in ms]:
                row[f"d_{kk}"] = round(float(ms[kk].mean() - base[kk].mean()), 4)
            rows.append(row)

        del P, T

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(tabdir, "virtual_ko_results.csv"), index=False,
               encoding="utf-8-sig")
    pd.DataFrame(edges).to_csv(os.path.join(tabdir, "grn_top_edges.csv"),
                               index=False, encoding="utf-8-sig")
    if res.empty:
        print("[error] no knockout results")
        return 1

    pd.set_option("display.width", 240)
    for ct in res["cell_type"].unique():
        sub = res[res["cell_type"] == ct].sort_values("delta_readout")
        print(f"\n=== virtual TF knockout: {ct} ===")
        print(sub.to_string(index=False))

    print("\n[NOTE] observational GRN + linear model + no chromatin layer ->")
    print("       prioritisation predictions, not causal proof.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
