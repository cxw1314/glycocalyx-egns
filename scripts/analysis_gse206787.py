#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE206787 -- mouse MI single-cell analysis for the EGNS framework.

Reads the 10x .h5 matrices, performs QC/clustering, annotates lineages,
computes the eight predefined module scores plus EGNS, and produces
sample-level pseudobulk summaries.

It additionally tests the prespecified candidate HPSE for CELL SPECIFICITY:
is it expressed by endothelium, or only by immune cells? (This is the one
criterion that whole blood cannot answer.)

IMPORTANT: from 3dPI and 42dPI there is only ONE animal each, so time-point
comparisons at those levels have no biological replication and are reported
as descriptive per-sample pseudobulk only.

Usage:
    python analysis_gse206787.py --indir data/raw/GSE206787/extracted \
                                 --outdir results/GSE206787
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")

LINEAGE_MARKERS = {
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Tek", "Emcn", "Cd36"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Kcnj8", "Abcc9", "Notch3", "Higd1b"],
    "SMC": ["Myh11", "Acta2", "Tagln", "Des"],
    "Fibroblast": ["Dcn", "Pdgfra", "Col1a1", "Postn", "Cthrc1"],
    "Macrophage": ["Cd68", "Adgre1", "Csf1r", "Lyz2", "Ccr2"],
    "Tcell": ["Cd3e", "Cd3d", "Cd3g", "Lck"],
    "Bcell": ["Cd79a", "Ms4a1", "Cd19"],
    "Nkcells": ["Nkg7", "Klrb1c", "Gzma"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg"],
    "LymphaticEC": ["Prox1", "Lyve1", "Pdpn"],
}

# prespecified candidate genes to inspect for cell specificity
CANDIDATES = ["Hpse", "Hpse2", "Hyal1", "Hyal2", "Mmp2", "Mmp9", "Mmp14",
              "Adam17", "Sulf1", "Sulf2", "Neu1", "PlauR", "Plaur"]


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        out[key] = {"genes": list(val.get("genes", [])),
                    "aliases": {k: list(v) for k, v in (val.get("aliases") or {}).items()}}
    return out


def resolve_genes(spec: dict, var_names) -> list[str]:
    """Case-insensitive symbol resolution so the human-symbol module
    definitions can be applied to mouse data (SDC1 -> Sdc1, HPSE -> Hpse)."""
    lookup: dict[str, str] = {}
    for v in var_names:
        lookup.setdefault(v.upper(), v)
    out: list[str] = []
    for g in spec["genes"]:
        hit = lookup.get(g.upper())
        if hit:
            out.append(hit)
    for canon, alts in spec["aliases"].items():
        if lookup.get(canon.upper()):
            continue
        for al in alts:
            hit = lookup.get(al.upper())
            if hit:
                out.append(hit)
                break
    return list(dict.fromkeys(out))


def sample_meta(path: str) -> tuple[str, str]:
    base = re.sub(r"_filtered_feature_bc_matrix\.h5$", "", os.path.basename(path))
    sample = base.split("_", 1)[0] + "_" + "_".join(base.split("_")[1:3])
    if re.search(r"42dPI", base):
        tp = "42dPI"
    elif re.search(r"7dPI", base):
        tp = "7dPI"
    elif re.search(r"3dPI", base):
        tp = "3dPI"
    else:
        tp = "undamaged"
    return sample, tp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=20.0)
    ap.add_argument("--max-cells-per-sample", type=int, default=25000)
    args = ap.parse_args()

    indir = args.indir if os.path.isabs(args.indir) else os.path.join(ROOT, args.indir)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir, figdir = os.path.join(outdir, "tables"), os.path.join(outdir, "figures")
    for d in (tabdir, figdir):
        os.makedirs(d, exist_ok=True)

    files = sorted(glob.glob(os.path.join(indir, "*filtered_feature_bc_matrix.h5")))
    files = [f for f in files if "Visium" not in os.path.basename(f)]
    if not files:
        print(f"[error] no 10x h5 files in {indir}")
        return 1
    print(f"[load] {len(files)} scRNA-seq samples")

    rng = np.random.default_rng(42)
    adatas = []
    for f in files:
        sample, tp = sample_meta(f)
        try:
            a = sc.read_10x_h5(f)
        except Exception as exc:  # noqa: BLE001
            print(f"  !! {os.path.basename(f)}: {exc}")
            continue
        a.var_names_make_unique()
        if a.n_obs > args.max_cells_per_sample:
            keep = rng.choice(a.n_obs, args.max_cells_per_sample, replace=False)
            a = a[np.sort(keep)].copy()
        a.obs["sample"] = sample
        a.obs["timepoint"] = tp
        a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
        adatas.append(a)
        print(f"  {sample:22s} {tp:10s} {a.n_obs:6d} cells x {a.n_vars} genes")

    if not adatas:
        print("[error] nothing loaded")
        return 1
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.var_names_make_unique()
    print(f"[concat] {adata.shape}")

    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    qc = adata.obs.groupby(["sample", "timepoint"])[
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"]].median()
    qc.to_csv(os.path.join(tabdir, "qc_medians_by_sample.csv"), encoding="utf-8-sig")
    print(qc.to_string())

    before = adata.n_obs
    adata = adata[(adata.obs["n_genes_by_counts"] >= args.min_genes) &
                  (adata.obs["pct_counts_mt"] <= args.max_mito)].copy()
    print(f"[qc] {before} -> {adata.n_obs} cells")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.8, key_added="cluster", flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(adata)
    print(f"[cluster] {adata.obs['cluster'].nunique()} clusters")

    for lin, genes in LINEAGE_MARKERS.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{lin}", use_raw=True)
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    lab = adata.obs.groupby("cluster")[score_cols].mean()
    lab["lineage"] = lab.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["lineage"] = adata.obs["cluster"].map(lab["lineage"]).astype(str)
    lab.to_csv(os.path.join(tabdir, "cluster_lineage_scores.csv"), encoding="utf-8-sig")
    print(lab["lineage"].value_counts().to_string())

    # ---- diagnostic: marker expression per cluster ----------------------
    KEY = ["Pecam1", "Cdh5", "Emcn", "Kdr", "Pdgfrb", "Rgs5", "Kcnj8",
           "Myh11", "Acta2", "Dcn", "Postn", "Cd68", "Ptprc", "Hpse"]
    present_key = [g for g in KEY if g in adata.raw.var_names]
    if present_key:
        m = adata.raw[:, present_key].X
        m = np.asarray(m.todense()) if hasattr(m, "todense") else np.asarray(m)
        mdf = pd.DataFrame(m, columns=present_key, index=adata.obs_names)
        mdf["cluster"] = adata.obs["cluster"].values
        cm = mdf.groupby("cluster", observed=True).mean()
        cm.to_csv(os.path.join(tabdir, "cluster_key_markers.csv"), encoding="utf-8-sig")
        print("\n[diagnostic] mean log-expression of key markers per cluster:")
        print(cm.round(2).to_string())
        print(f"\n[diagnostic] cells with Pecam1>0: "
              f"{int((mdf['Pecam1'] > 0).sum()) if 'Pecam1' in mdf else 'n/a'} "
              f"/ {len(mdf)}")

    modules = load_modules()
    for key, spec in modules.items():
        present = resolve_genes(spec, adata.raw.var_names)
        if present:
            sc.tl.score_genes(adata, present, score_name=f"mod_{key}", use_raw=True)
            print(f"  [module] {key}: {len(present)} genes matched")

    def zc(col):
        s = adata.obs[col]
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd else s * 0.0

    have = set(adata.obs.columns)
    if {"mod_GDS", "mod_EAS", "mod_EIDS", "mod_GCI", "mod_EIS"} <= have:
        adata.obs["EGNS"] = (zc("mod_GDS") + zc("mod_EAS") + zc("mod_EIDS")
                             - zc("mod_GCI") - zc("mod_EIS"))

    # ---- candidate cell specificity (the key output) --------------------
    ec_mask = adata.obs["lineage"].isin(["Endothelial", "LymphaticEC"])
    rows = []
    for g in CANDIDATES:
        if g not in adata.raw.var_names:
            continue
        vec = adata.raw[:, g].X
        vec = np.asarray(vec.todense()).ravel() if hasattr(vec, "todense") else np.asarray(vec).ravel()
        raw_counts = adata.layers["counts"][:, adata.var_names.get_loc(g)].toarray().ravel() \
            if g in adata.var_names else None
        for lin in sorted(adata.obs["lineage"].unique()):
            m = (adata.obs["lineage"] == lin).values
            if m.sum() < 20:
                continue
            rows.append({"gene": g, "lineage": lin, "n_cells": int(m.sum()),
                         "mean_log_expr": float(np.mean(vec[m])),
                         "pct_expressing": float(100.0 * np.mean(vec[m] > 0))})
    spec = pd.DataFrame(rows)
    if len(spec):
        spec.to_csv(os.path.join(tabdir, "candidate_cell_specificity.csv"),
                    index=False, encoding="utf-8-sig")
        piv = spec.pivot(index="gene", columns="lineage", values="pct_expressing")
        print("\n[candidate cell specificity] % of cells expressing:")
        print(piv.round(1).to_string())
        ec_rank = []
        for g in piv.index:
            v = piv.loc[g]
            ec_v = v.get("Endothelial", np.nan)
            others = v.drop(labels=[l for l in ["Endothelial", "LymphaticEC"] if l in v.index])
            top = others.idxmax() if len(others) else None
            ec_rank.append({"gene": g, "pct_endothelial": ec_v,
                            "pct_max_other_lineage": float(others.max()) if len(others) else np.nan,
                            "top_other_lineage": top,
                            "endothelial_rank": (int(v.rank(ascending=False)["Endothelial"])
                                                 if "Endothelial" in v.index else -1)})
        pd.DataFrame(ec_rank).to_csv(
            os.path.join(tabdir, "candidate_endothelial_ranking.csv"),
            index=False, encoding="utf-8-sig")

    # ---- endothelial pseudobulk by sample/timepoint ---------------------
    mod_cols = [c for c in adata.obs.columns if c.startswith("mod_")]
    ec = adata[adata.obs["lineage"].isin(["Endothelial", "LymphaticEC"])].copy()
    print(f"\n[endothelium] {ec.n_obs} cells")
    if ec.n_obs > 0:
        cols = mod_cols + (["EGNS"] if "EGNS" in ec.obs else [])
        pb = ec.obs.groupby(["sample", "timepoint"], observed=True)[cols].mean()
        pb["n_cells"] = ec.obs.groupby(["sample", "timepoint"], observed=True).size()
        pb.to_csv(os.path.join(tabdir, "endothelial_pseudobulk_scores.csv"),
                  encoding="utf-8-sig")
        print(pb.round(3).to_string())

    adata.obs[["sample", "timepoint", "cluster", "lineage"] + mod_cols +
              (["EGNS"] if "EGNS" in adata.obs else [])].to_csv(
        os.path.join(tabdir, "cell_scores.csv"), encoding="utf-8-sig")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sc.pl.umap(adata, color=["lineage", "timepoint"], show=False, wspace=0.4)
        plt.savefig(os.path.join(figdir, "umap_lineage_timepoint.png"), dpi=200,
                    bbox_inches="tight"); plt.close()
        cols = [c for c in ["mod_GCI", "mod_GDS", "mod_EIS", "EGNS"] if c in adata.obs]
        if cols:
            sc.pl.umap(adata, color=cols, show=False, cmap="viridis")
            plt.savefig(os.path.join(figdir, "umap_module_scores.png"), dpi=200,
                        bbox_inches="tight"); plt.close()
        print(f"[fig] {figdir}")
    except Exception as exc:  # noqa: BLE001
        print(f"[fig] skipped: {exc}")

    os.makedirs(outdir, exist_ok=True)
    adata.write_h5ad(os.path.join(outdir, "gse206787_processed.h5ad"))
    print(f"[ok] saved -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
