#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CellOracle-STYLE OFFLINE virtual perturbation (NOT official CellOracle).

Why offline
-----------
Genuine CellOracle cannot be installed on this machine: it imports
`velocyto.serialization / diffusion / estimation`, and velocyto publishes no
binary wheel for Python 3.12 while its source build requires pysam -> a C
compiler, which is absent here (no cl.exe/gcc/g++/clang). No stub is used,
because `colDeltaCor` and `Diffusion` are the core functions the perturbation
simulation actually depends on.

What this script implements (the same *quantities* CellOracle reports)
---------------------------------------------------------------------
  1. PCA embedding of the lineage's panel-gene expression.
  2. A reference vector field: for each cell, the direction in PCA space
     towards its higher-pseudotime neighbours (diffusion pseudotime, root =
     the least "non-resolved" cell of that lineage).
  3. A knockout shift: dX = -1 * X_centred[:, ko] * B[:, ko], with B the gene
     regulatory network (ridge TF->target for endothelium; PC-regression
     gene->gene for the sheddase layer). dX is projected into PCA space.
  4. An INNER PRODUCT SCORE: per cell, the cosine between the reference
     transition vector and the knockout shift vector; the reported score is
     its mean. A drop means the knockout pushes cells off their normal
     trajectory.
  5. Top perturbed genes (mean |dX| per gene) and hypergeometric enrichment of
     those genes against the eight predefined biological modules.

Explicit limitation: this is a CellOracle-STYLE reimplementation. The GRN is
data-driven with no chromatin-accessibility layer, and results are hypothesis
ranking, not causal proof. All outputs are labelled accordingly.

Usage:
    python analysis_celloracle_offline.py \
        --h5ad results/GSE176092/gse176092_nonCM_processed.h5ad \
        --outdir results/celloracle_offline
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

TF_PANEL = ["KLF2", "KLF4", "RELA", "NFKB1", "NFKB2", "RELB", "ERG", "FLI1",
            "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "STAT3",
            "HIF1A", "SMAD1", "SMAD3", "JUN", "JUNB", "FOS", "FOSL2", "SOX17",
            "SOX18", "PROX1", "EPAS1", "PPARG", "MEFC2", "MEF2C", "GATA6", "NFATC1"]
SHEDDASE = ["Hpse", "Mmp9", "Plaur", "Mmp2", "Mmp14", "Adam17", "Hyal2",
            "Ctsl", "Adamts1", "Adamts4"]


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    return {k: [g.upper() for g in v.get("genes", [])]
            for k, v in cfg.items()
            if isinstance(v, dict) and "genes" in v}


def zscore(M: np.ndarray) -> np.ndarray:
    mu = M.mean(axis=0, keepdims=True)
    sd = M.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (M - mu) / sd


def transition_field(emb: np.ndarray, pt: np.ndarray, k: int = 30) -> np.ndarray:
    """Direction in embedding space towards higher-pseudotime neighbours."""
    from sklearn.neighbors import NearestNeighbors
    kk = min(k + 1, len(emb))
    nn = NearestNeighbors(n_neighbors=kk).fit(emb)
    _, idx = nn.kneighbors(emb)
    v = np.zeros_like(emb, dtype=float)
    for i in range(len(emb)):
        nb = idx[i, 1:]
        if not len(nb):
            continue
        higher = nb[pt[nb] > pt[i]]
        tgt = higher if len(higher) else nb
        v[i] = emb[tgt].mean(axis=0) - emb[i]
    return v


def hypergeom_enrich(genes_hit: list[str], universe: list[str],
                     modules: dict) -> pd.DataFrame:
    from scipy.stats import hypergeom
    N = len(universe)
    uset = set(universe)
    rows = []
    hit = [g for g in genes_hit if g in uset]
    K = len(hit)
    if K == 0:
        return pd.DataFrame()
    for mod, gs in modules.items():
        mod_in = [g for g in gs if g in uset]
        n = len(mod_in)
        if n == 0:
            continue
        x = len(set(hit) & set(mod_in))
        p = float(hypergeom.sf(x - 1, N, n, K)) if x > 0 else 1.0
        rows.append({"module": mod, "n_hit_in_module": x, "n_module": n,
                     "n_top": K, "universe": N,
                     "fold": round((x / K) / (n / N), 3) if K and n else np.nan,
                     "p": p})
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values("p")
        m = len(out)
        out["q"] = np.minimum(1.0, out["p"] * m / np.maximum(
            np.arange(1, m + 1), 1))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default="results/celloracle_offline")
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--knn", type=int, default=30)
    ap.add_argument("--n-top", type=int, default=20)
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--ridge-alpha", type=float, default=0.1)
    args = ap.parse_args()

    import anndata as ad
    import scanpy as sc
    from sklearn.decomposition import PCA
    from sklearn.linear_model import Ridge

    def rp(p):
        return p if os.path.isabs(p) else os.path.join(ROOT, p)

    outdir = rp(args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(rp(args.h5ad))
    raw = adata.raw if adata.raw is not None else adata
    lut: dict[str, str] = {}
    for v in raw.var_names:
        lut.setdefault(v.upper(), v)
    modules = load_modules()
    panel = sorted({g for gs in modules.values() for g in gs if g in lut})
    pidx = {g: i for i, g in enumerate(panel)}
    print(f"[load] {adata.n_obs} cells; panel = {len(panel)} genes")

    def dense(genes):
        M = raw[:, [lut[g] for g in genes]].X
        M = np.asarray(M.todense()) if hasattr(M, "todense") else np.asarray(M)
        rs = M.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        return np.log1p(M / rs * 1e4)

    rng = np.random.default_rng(SEED)
    ip_rows, top_rows, enr_rows = [], [], []

    for lineage in ("Endothelial", "Macrophage"):
        mask = adata.obs["lineage"].astype(str).str.contains(lineage).values
        n = int(mask.sum())
        if n < 100:
            continue

        # ---- expression + embedding ------------------------------------
        P = dense(panel)[mask]
        kk = min(args.knn, n - 1)
        from sklearn.neighbors import NearestNeighbors
        nn0 = NearestNeighbors(n_neighbors=kk).fit(zscore(P))
        _, nb0 = nn0.kneighbors(zscore(P))
        Pk = P[nb0].mean(axis=1)                      # kNN-imputed targets
        Pz = zscore(Pk)
        ncomp = min(args.n_pcs, Pz.shape[1] - 1, n - 1)
        pca = PCA(n_components=ncomp, random_state=SEED).fit(Pz)
        emb = pca.transform(Pz)
        print(f"\n=== {lineage}: {n} cells, PCA {ncomp} comps, "
              f"EV={pca.explained_variance_ratio_.sum():.3f}")

        # ---- diffusion pseudotime --------------------------------------
        a = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32))
        a.obsm["X_pca"] = emb
        sc.pp.neighbors(a, n_neighbors=min(30, n - 1), use_rep="X_pca")
        if lineage == "Endothelial":
            gds = [pidx[g] for g in modules["GDS"] if g in pidx]
            gci = [pidx[g] for g in modules["GCI"] if g in pidx]
            # higher = more non-resolved. The root must be the LEAST
            # non-resolved cell so that pseudotime (and therefore the
            # reference vector field) points TOWARDS non-resolution.
            # A previous version negated this metric and so rooted at the
            # MOST non-resolved cell, silently inverting the field.
            root_metric = Pk[:, gds].mean(axis=1) - Pk[:, gci].mean(axis=1)
        else:
            mfs = [pidx[g] for g in modules["MFS"] if g in pidx]
            root_metric = Pk[:, mfs].mean(axis=1)      # higher = more pro-fibrotic
        a.uns["iroot"] = int(np.argmin(root_metric))
        sc.tl.diffmap(a, n_comps=15)
        sc.tl.dpt(a)
        pt = np.asarray(a.obs["dpt_pseudotime"].values, dtype=float)
        pt = np.nan_to_num(pt, nan=0.0)
        # --- sanity check: pseudotime MUST track the target state --------
        from scipy.stats import spearmanr
        rho_pt = float(spearmanr(pt, root_metric).correlation)
        print(f"[dpt] range {pt.min():.3f} - {pt.max():.3f}; "
              f"Spearman(pseudotime, target state) = {rho_pt:+.3f}")
        if rho_pt < 0.3:
            print("  [WARN] pseudotime does not track the target state; "
                  "the vector-field direction is NOT interpretable")

        # ---- reference transition field --------------------------------
        V = transition_field(emb, pt, k=args.knn)

        # ---- GRN --------------------------------------------------------
        if lineage == "Endothelial":
            tfs = [t for t in TF_PANEL if t in lut]
            T = dense(tfs)[mask]
            T = T[nb0].mean(axis=1)
            B = np.zeros((len(panel), len(tfs)))
            for t in range(len(panel)):
                y = Pk[:, t]
                if y.std() == 0:
                    continue
                B[t] = Ridge(alpha=args.ridge_alpha).fit(T, y).coef_
            reg_names = tfs
            Bg = B                                   # target x regulator
            Xreg = T
        else:
            # PC-regression gene->gene network (edge-ablation layer)
            var = Pk.var(axis=0)
            keep = np.argsort(var)[::-1][:min(141, len(panel))]
            Xs = Pk[:, keep]
            gs = [panel[i] for i in keep]
            k_nb = max(10, int(0.1 * n))
            nnb = NearestNeighbors(n_neighbors=min(k_nb, n - 1)).fit(Xs)
            _, ALL = nnb.kneighbors(Xs)
            Bg = np.zeros((len(gs), len(gs)), dtype=np.float32)
            for j in range(len(gs)):
                Y = Xs[ALL[j]]
                p = min(3, Y.shape[0] - 1, Y.shape[1])
                pc = PCA(n_components=p, random_state=SEED).fit(Y)
                W = pc.transform(Y)
                coef = np.linalg.lstsq(W, Y[:, j], rcond=None)[0]
                Bg[:, j] = pc.components_.T @ coef
            reg_names = gs
            Bg = Bg.T                                # target x regulator
            Xreg = Xs
            Pk = Xs
            panel_local = gs
        if lineage == "Endothelial":
            panel_local = panel
        ridx: dict[str, int] = {}
        for i, g in enumerate(reg_names):
            ridx.setdefault(g, i)
            ridx.setdefault(g.upper(), i)          # SHEDDASE list is Title-case
        pidx_local = {g: i for i, g in enumerate(panel_local)}

        Xc = Pk - Pk.mean(axis=0, keepdims=True)
        load = pca.components_.T                     # n_panel x n_comp
        if Pk.shape[1] != load.shape[0]:
            # PC-regression branch uses a gene subset; refit PCA on it
            pca2 = PCA(n_components=min(ncomp, Pk.shape[1] - 1),
                       random_state=SEED).fit(zscore(Pk))
            load = pca2.components_.T
            emb = pca2.transform(zscore(Pk))
            V = transition_field(emb, pt, k=args.knn)

        cand = [g for g in (TF_PANEL if lineage == "Endothelial" else SHEDDASE)
                if g in ridx or g.upper() in ridx]
        print(f"[KO candidates] {', '.join(cand)}")

        # ---- per-knockout inner product score --------------------------
        def score_for(reg: str):
            j = ridx[reg] if reg in ridx else ridx[reg.upper()]
            dX = -np.outer(Xc[:, j], Bg[:, j])       # cells x targets
            dPC = dX @ load
            num = (V * dPC).sum(axis=1)
            den = (np.linalg.norm(V, axis=1) * np.linalg.norm(dPC, axis=1))
            cos = np.where(den > 1e-12, num / np.maximum(den, 1e-12), np.nan)
            return float(np.nanmean(cos)), dX, cos

        all_regs = list(reg_names)
        null = []
        for _ in range(args.n_null):
            r = all_regs[int(rng.integers(len(all_regs)))]
            null.append(score_for(r)[0])
        null = np.asarray(null, dtype=float)
        null = null[np.isfinite(null)]
        sd = float(null.std()) if len(null) else np.nan
        print(f"[null] {len(null)} random-gene perturbations, sd = {sd:.4f}")

        ko_stats, vf_store = {}, {}
        for g in cand:
            s, dX, cos = score_for(g)
            ip_rows.append({"lineage": lineage, "ko": g,
                            "inner_product_score": round(s, 4),
                            "spearman_pt_vs_target_state": round(rho_pt, 3),
                            "z_vs_null": round((s - float(np.nanmean(null))) / (sd or 1), 3)
                            if sd else np.nan,
                            "p_empirical": round(float(
                                (np.sum(np.abs(null) >= abs(s)) + 1) / (len(null) + 1)), 4),
                            "null_mean": round(float(np.nanmean(null)), 4),
                            "null_sd": round(sd, 4), "n_cells": n})
            ko_stats[g] = s
            vf_store[g] = (dX @ load)[:, :2].mean(axis=0)
            mag = np.abs(dX).mean(axis=0)
            order = np.argsort(mag)[::-1][:args.n_top]
            hits = [panel_local[i] for i in order]
            for r, i in enumerate(order, 1):
                top_rows.append({"lineage": lineage, "ko": g, "rank": r,
                                 "gene": panel_local[i],
                                 "mean_abs_shift": float(mag[i])})
            enr = hypergeom_enrich(hits, panel_local, modules)
            if len(enr):
                enr.insert(0, "ko", g)
                enr.insert(0, "lineage", lineage)
                enr_rows.append(enr)

        np.savez(os.path.join(outdir, f"vector_field_{lineage}.npz"),
                 emb=emb[:, :2], V=V[:, :2], pt=pt,
                 **{f"shift_{g}": vf_store[g] for g in vf_store})
        print(f"[saved] vector field for {lineage}")

    pd.DataFrame(ip_rows).to_csv(
        os.path.join(tabdir, "inner_product_scores.csv"),
        index=False, encoding="utf-8-sig")
    pd.DataFrame(top_rows).to_csv(
        os.path.join(tabdir, "top_perturbed_genes.csv"),
        index=False, encoding="utf-8-sig")
    if enr_rows:
        pd.concat(enr_rows, ignore_index=True).to_csv(
            os.path.join(tabdir, "module_enrichment.csv"),
            index=False, encoding="utf-8-sig")

    ip = pd.DataFrame(ip_rows)
    print("\n=== inner product score (CellOracle-STYLE offline) ===")
    print(ip.to_string(index=False))
    print("\n[NOTE] CellOracle-style offline reimplementation, NOT official "
          "CellOracle; hypothesis ranking, not causal proof.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
