#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
scTenifoldKnk-equivalent virtual knockout for NON-transcription-factor effectors
(HPSE, MMP9, PLAUR).

Why this instead of CellOracle
------------------------------
CellOracle is a TF-perturbation method: it needs a TF-motif base GRN, and its
propagation is defined on a TF-target regulatory graph. HPSE/MMP9/PLAUR are
sheddases, not TFs, so a TF-motif prior does not apply to them. scTenifoldKnk
instead builds a single-cell GRN, deletes the knockout gene's outgoing edges,
and compares the wild-type and knockout networks. This script reproduces that
logic with the same manifold-distance idea:

  1. build a gene regulatory coefficient matrix B by principal-component
     regression (k = 3 PCs of each gene's kNN neighbourhood), scTenifoldNet style;
  2. zero out the row of the knockout gene (its outgoing edges) -> B_ko;
  3. quantify the network-level distance between B_wt and B_ko for each
     candidate KO and rank the sheddases by how strongly they rewire the
     macrophage fibrotic / endothelial-glycocalyx programmes.

Output:
    results/sctenifold_sheddase_ko/tables/sheddase_ko_effects.csv
    results/sctenifold_sheddase_ko/tables/sheddase_ko_network_distance.csv
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

CANDIDATES = ["Hpse", "Mmp9", "Plaur", "Mmp2", "Mmp14", "Adam17",
              "Hyal2", "Ctsl", "Adamts1", "Adamts4"]
# everything the module panel covers, plus the sheddases themselves
PANEL_GENES = ["GDS", "MFS", "FNS", "GCI", "EAS", "GDS_EXTRA"]


def z(v):
    sd = np.std(v)
    return (v - v.mean()) / sd if sd else v * 0.0


def pc_regression_grn(X: np.ndarray, genes: list[str], k: int = 3,
                      max_genes: int = 400, seed: int = 0) -> np.ndarray:
    """B[i, j] = coefficient of gene i (as a target) on gene j (as a regulator).

    scTenifoldNet/sciPCNet style: for each gene, take its kNN neighbourhood,
    extract k principal components of that neighbourhood, regress the gene on
    those PCs, and project the coefficients back onto the original gene axes.
    """
    rng = np.random.default_rng(seed)
    n_genes = X.shape[1]
    # keep the highest-variance genes to bound compute
    var = X.var(axis=0)
    idx = np.argsort(var)[::-1][:min(max_genes, n_genes)]
    Xs = X[:, idx]
    gs = [genes[i] for i in idx]
    n = Xs.shape[1]

    # neighbourhood graph is built once, scTenifoldNet style, over the panel
    k_nb = max(10, int(0.1 * X.shape[0]))
    from sklearn.neighbors import NearestNeighbors
    nb = NearestNeighbors(n_neighbors=min(k_nb, Xs.shape[0] - 1)).fit(Xs)
    _, ALL_NB = nb.kneighbors(Xs)
    from sklearn.decomposition import PCA

    B = np.zeros((n, n), dtype=np.float32)
    for j in range(n):
        Y = Xs[ALL_NB[j]]                  # cells x genes, neighbourhood of j
        p = min(k, Y.shape[0] - 1, Y.shape[1])
        pca = PCA(n_components=p, random_state=seed)
        W = pca.fit_transform(Y)            # cells x p
        # regress only gene j on the neighbourhood PCs, then project the
        # coefficients back onto gene space -> column j of B holds j's regulators
        yj = Y[:, j]
        coef = np.linalg.lstsq(W, yj, rcond=None)[0]   # p x 1
        B[:, j] = pca.components_.T @ coef             # n_genes
    return B, gs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default="results/sctenifold_sheddase_ko")
    ap.add_argument("--lineage", default="Macrophage")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import anndata as ad

    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    mask = adata.obs["lineage"].astype(str).str.contains(args.lineage).values
    print(f"[{args.lineage}] {int(mask.sum())} cells")

    lut = {}
    for v in raw.var_names:
        lut.setdefault(v.upper(), v)

    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    modules = {k: [g.upper() for g in v.get("genes", [])]
               for k, v in cfg.items()
               if isinstance(v, dict) and "genes" in v}

    # gene universe: module genes + sheddase candidates + top variable genes
    universe = set()
    for gs in modules.values():
        universe.update(gs)
    universe.update(g.upper() for g in CANDIDATES)
    present = [lut[g] for g in universe if g in lut]

    M = raw[mask][:, present].X
    M = np.asarray(M.todense()) if hasattr(M, "todense") else np.asarray(M)
    M = np.log1p(M / np.maximum(M.sum(axis=1, keepdims=True), 1e-9) * 1e4)
    genes = [g.upper() for g in present]
    print(f"[universe] {len(present)} panel genes")

    B, gs = pc_regression_grn(M, genes, seed=args.seed)
    gidx = {g: i for i, g in enumerate(gs)}
    ko_present = [c for c in CANDIDATES if c.upper() in gidx]
    print(f"[grn] {B.shape[0]} x {B.shape[1]} coefficient matrix; "
          f"KO candidates present: {', '.join(ko_present)}")

    # module scores in this lineage, as the readout
    def module_score(geneset):
        gg = [g for g in geneset if g in gidx]
        if not gg:
            return None
        return M[:, [gidx[g] for g in gg]].mean(axis=1)

    reads = {}
    for mod in ("MFS", "FNS", "GDS", "EAS", "GCI"):
        s = module_score(modules.get(mod, []))
        if s is not None:
            reads[mod] = s

    rows, dist = [], []
    for ko in ko_present:
        j = gidx[ko.upper()]
        B_ko = B.copy()
        B_ko[j, :] = 0.0                       # delete outgoing edges
        # manifold distance: Frobenius norm of the change, normalised
        fro = float(np.linalg.norm(B - B_ko, ord="fro") / np.linalg.norm(B, ord="fro"))
        spearman = float(pd.Series(B[j, :]).corr(pd.Series(B_ko[j, :]),
                                                 method="spearman"))
        dist.append({"ko": ko, "relative_frobenius_change": round(fro, 6),
                     "outgoing_edges_deleted": int((B[j, :] != 0).sum())})
        # effect on the knock-in programme: mean |column change| weighted by
        # whether the column gene is in each readout
        colchg = np.abs(B - B_ko).sum(axis=0)
        for mod, s in reads.items():
            gg = [gidx[g] for g in modules.get(mod, []) if g in gidx]
            if not gg:
                continue
            eff = float(colchg[gg].mean())
            base = float(colchg.mean())
            rows.append({"ko": ko, "readout_module": mod,
                         "mean_|B_change|_in_module": round(eff, 6),
                         "genome_wide_mean": round(base, 6),
                         "effect_ratio": round(eff / base, 4) if base else np.nan})
    eff = pd.DataFrame(rows)
    eff.to_csv(os.path.join(tabdir, "sheddase_ko_effects.csv"),
               index=False, encoding="utf-8-sig")
    dd = pd.DataFrame(dist).sort_values("relative_frobenius_change", ascending=False)
    dd.to_csv(os.path.join(tabdir, "sheddase_ko_network_distance.csv"),
              index=False, encoding="utf-8-sig")

    print("\n=== network-level rewiring by sheddase KO (macrophage scGRN) ===")
    print(dd.to_string(index=False))

    print("\n=== programme-specific effect (effect_ratio > 1 = above genome-wide) ===")
    piv = eff.pivot_table(index="ko", columns="readout_module",
                          values="effect_ratio")
    order = [m for m in ("MFS", "FNS", "GDS", "EAS", "GCI") if m in piv.columns]
    piv = piv[order].round(3)
    print(piv.to_string())

    rank = piv.get("MFS", pd.Series(dtype=float)).sort_values(ascending=False)
    if len(rank):
        print("\n[ranking] macrophage-profibrotic programme (MFS), highest first:")
        for i, (k, v) in enumerate(rank.items(), 1):
            print(f"   {i}. {k:8s} effect_ratio = {v:.3f}")
    print("\n[ok] scTenifoldKnk-equivalent effector KO written to " + tabdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
