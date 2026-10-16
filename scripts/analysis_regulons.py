#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Regulatory-network analysis of high-EGNS endothelium (offline DoRothEA/SCENIC
equivalent).

DoRothEA and SCENIC are R packages and their regulon resources are fetched
online, neither of which is available here. We therefore use a curated set of
canonical TF regulons for the transcription factors named in the study design
(KLF2/KLF4, NF-kB/RELA, SMAD1/3, AP-1, YAP-TAZ/TEAD, plus STAT3, HIF1A and
ERG as context), score their activity with a rank-based signature score, and
test whether activity differs between high- and low-EGNS endothelium.

Usage:
    python analysis_regulons.py --h5ad results/GSE176092/gse176092_nonCM_processed.h5ad \
                               --outdir results/regulons
"""
from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 42

# Curated TF regulons (mouse symbols). The TF itself is retained, as in
# DoRothEA regulons, so the signature reflects the whole module.
REGULONS: dict[str, list[str]] = {
    "KLF2_KLF4": ["Klf2", "Klf4", "Nos3", "Thbd", "Procr", "S1pr1", "Gja4",
                  "Vwf", "Emcn", "Sox18", "Cdh5"],
    "NFKB_RELA": ["Rela", "Nfkb1", "Nfkb2", "Relb", "Nfkbia", "Nfkbiz",
                  "Tnfaip3", "Il6", "Ccl2", "Vcam1", "Icam1", "Cxcl1", "Cxcl2"],
    "SMAD1_3": ["Smad1", "Smad3", "Smad6", "Smad7", "Serpine1", "Junb",
                "Id1", "Id2", "Col1a1", "Ccn2", "Cyr61"],
    "AP1": ["Jun", "Junb", "Jund", "Fos", "Fosb", "Fosl1", "Fosl2", "Mmp9",
            "Mmp13", "Ccl2", "Il6", "Timp1"],
    "YAP_TAZ_TEAD": ["Yap1", "Wwtr1", "Tead1", "Tead2", "Tead3", "Tead4",
                     "Ccn1", "Ccn2", "Amotl2", "Ankrd1", "Cyr61"],
    "STAT3": ["Stat3", "Socs3", "Il6", "Mmp9", "Vegfa", "Bcl2l1", "Hif1a"],
    "HIF1A": ["Hif1a", "Vegfa", "Slc2a1", "Ldha", "Pdk1", "Bnip3"],
    "ERG_FLI1": ["Erg", "Fli1", "Cdh5", "Pecam1", "Klf2", "Tek", "Sox18", "Vwf"],
}

EGNS_PARTS = {"GDS": "mod_GDS", "EAS": "mod_EAS", "EIDS": "mod_EIDS",
              "GCI": "mod_GCI", "EIS": "mod_EIS"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    import anndata as ad
    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    var_names = list(raw.var_names)
    lookup: dict[str, str] = {}
    for v in var_names:
        lookup.setdefault(v.upper(), v)
    print(f"[load] {adata.n_obs} cells")

    X = raw.X

    def gene_vec(sym: str) -> np.ndarray | None:
        real = lookup.get(sym.upper())
        if real is None:
            return None
        j = var_names.index(real)
        col = X[:, j]
        return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()

    # library-size normalisation of the raw counts for fair comparison
    if hasattr(X, "sum"):
        lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    scale = 1e4 / lib

    regulon_scores = {}
    coverage = {}
    for tf, genes in REGULONS.items():
        vecs = []
        for g in genes:
            v = gene_vec(g)
            if v is not None:
                vecs.append(v * scale)
        coverage[tf] = len(vecs)
        if vecs:
            regulon_scores[tf] = np.mean(np.log1p(np.vstack(vecs)), axis=0)
    cov_df = pd.DataFrame([{"regulon": k, "n_genes_found": v} for k, v in coverage.items()])
    cov_df.to_csv(os.path.join(tabdir, "regulon_coverage.csv"), index=False,
                  encoding="utf-8-sig")
    print("[coverage]")
    print(cov_df.to_string(index=False))

    rdf = pd.DataFrame(regulon_scores, index=adata.obs_names)
    for c, v in EGNS_PARTS.items():
        if v in adata.obs.columns:
            rdf[f"mod_{c}"] = adata.obs[v].values
    rdf["lineage"] = adata.obs["lineage"].values
    rdf["condition"] = adata.obs["condition"].values
    rdf["sample"] = adata.obs["sample"].values
    rdf.to_csv(os.path.join(tabdir, "regulon_scores_all_cells.csv"), encoding="utf-8-sig")

    # ---- endothelial cells ---------------------------------------------
    ec_mask = adata.obs["lineage"].astype(str).str.contains("Endothelial").values
    print(f"\n[endothelium] {int(ec_mask.sum())} cells")
    if ec_mask.sum() < 50:
        print("[error] too few endothelial cells")
        return 1
    ec = rdf[ec_mask].copy()

    ec["EGNS"] = adata.obs.loc[ec.index, "EGNS"].values if "EGNS" in adata.obs else np.nan

    tfs = list(regulon_scores.keys())

    # correlation of regulon activity with EGNS within endothelium
    rows = []
    for tf in tfs:
        rho = stats.spearmanr(ec[tf], ec["EGNS"])
        rows.append({"regulon": tf, "rho_vs_EGNS": round(float(rho.statistic), 4),
                     "p": float(f"{rho.pvalue:.3g}")})
    corr = pd.DataFrame(rows).sort_values("rho_vs_EGNS", ascending=False)
    corr.to_csv(os.path.join(tabdir, "regulon_vs_EGNS_correlation.csv"), index=False,
                encoding="utf-8-sig")
    print("\n=== regulon activity vs EGNS in endothelium ===")
    print(corr.to_string(index=False))

    # high vs low EGNS quartiles
    q75, q25 = ec["EGNS"].quantile([0.75, 0.25])
    hi = ec["EGNS"] >= q75
    lo = ec["EGNS"] <= q25
    rows = []
    for tf in tfs:
        u = stats.mannwhitneyu(ec.loc[hi, tf], ec.loc[lo, tf])
        rows.append({"regulon": tf,
                     "mean_high_EGNS": round(float(ec.loc[hi, tf].mean()), 4),
                     "mean_low_EGNS": round(float(ec.loc[lo, tf].mean()), 4),
                     "delta": round(float(ec.loc[hi, tf].mean() - ec.loc[lo, tf].mean()), 4),
                     "p": float(f"{u.pvalue:.3g}")})
    hlv = pd.DataFrame(rows).sort_values("delta", ascending=False)
    hlv.to_csv(os.path.join(tabdir, "regulon_high_vs_low_EGNS.csv"), index=False,
               encoding="utf-8-sig")
    print(f"\n=== regulon activity, high-EGNS (n={int(hi.sum())}) vs low-EGNS "
          f"(n={int(lo.sum())}) endothelium ===")
    print(hlv.to_string(index=False))

    # per-animal pseudobulk time course
    tc = ec.groupby(["sample", "condition"])[tfs].mean()
    tc["n_cells"] = ec.groupby(["sample", "condition"]).size()
    tc.to_csv(os.path.join(tabdir, "regulon_pseudobulk_by_sample.csv"),
              encoding="utf-8-sig")
    agg = ec.groupby("condition")[tfs].mean()
    agg.to_csv(os.path.join(tabdir, "regulon_pseudobulk_by_condition.csv"),
               encoding="utf-8-sig")
    print("\n=== regulon activity by condition (endothelium, cell-level mean) ===")
    print(agg.round(4).to_string())

    print("\n[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
