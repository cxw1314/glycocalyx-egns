#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE176092 (non-cardiomyocyte scRNA-seq) -- mouse MI analysis for EGNS.

This dataset retains endothelial cells (unlike the mural-cell-focused
GSE206787), so it is used to (a) compute the module scores and EGNS on
endothelium and (b) adjudicate the CELL SPECIFICITY of the prespecified
candidate genes.

Design: Sham (n=2), MI day1 (n=2), MI day7 (n=2) -> sample-level pseudobulk.

Usage:
    python analysis_gse176092.py --indir data/raw/GSE176092/extracted \
                                 --outdir results/GSE176092
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import tempfile
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
    "Macrophage": ["Cd68", "Adgre1", "Csf1r", "Lyz2", "Ccr2", "Cd163"],
    "Tcell": ["Cd3e", "Cd3d", "Cd3g", "Lck"],
    "Bcell": ["Cd79a", "Ms4a1", "Cd19"],
    "Nkcells": ["Nkg7", "Klrb1c", "Gzma"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg"],
    "LymphaticEC": ["Prox1", "Lyve1", "Pdpn"],
    "Adipocyte": ["Adipoq", "Plin1", "Fabp4"],
}

CANDIDATES = ["Hpse", "Hpse2", "Hyal1", "Hyal2", "Mmp2", "Mmp9", "Mmp14",
              "Adam17", "Sulf1", "Sulf2", "Neu1", "Plaur", "Ext1", "Gpc1",
              "Sdc1", "Has2", "Ndst1"]

KEY_MARKERS = ["Pecam1", "Cdh5", "Emcn", "Kdr", "Pdgfrb", "Rgs5", "Myh11",
               "Dcn", "Postn", "Cd68", "Ptprc", "Hpse", "Mmp2", "Mmp9"]


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
    base = os.path.basename(path)
    if "Sham" in base:
        cond = "Sham"
    elif "day1" in base:
        cond = "MI_day1"
    elif "day7" in base:
        cond = "MI_day7"
    else:
        cond = "unknown"
    return base.replace("_matrix.mtx.gz", ""), cond


def load_mtx(prefix: str) -> ad.AnnData:
    tmp = tempfile.mkdtemp()
    try:
        for kind, pat in (("barcodes", "_barcodes.tsv.gz"),
                          ("features", "_features.tsv.gz"),
                          ("matrix", "_matrix.mtx.gz")):
            src = prefix + pat
            if not os.path.exists(src):
                alt = glob.glob(prefix + "*" + pat.split("_")[-1])
                if not alt:
                    raise FileNotFoundError(src)
                src = alt[0]
            dst = os.path.join(tmp, {"barcodes": "barcodes.tsv.gz",
                                     "features": "features.tsv.gz",
                                     "matrix": "matrix.mtx.gz"}[kind])
            shutil.copy(src, dst)
        return sc.read_10x_mtx(tmp, var_names="gene_symbols", make_unique=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=20.0)
    args = ap.parse_args()

    indir = args.indir if os.path.isabs(args.indir) else os.path.join(ROOT, args.indir)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir, figdir = os.path.join(outdir, "tables"), os.path.join(outdir, "figures")
    for d in (tabdir, figdir):
        os.makedirs(d, exist_ok=True)

    prefixes = sorted(p[:-len("_matrix.mtx.gz")] for p in
                      glob.glob(os.path.join(indir, "*_matrix.mtx.gz")))
    if not prefixes:
        print(f"[error] no 10x matrices in {indir}")
        return 1
    print(f"[load] {len(prefixes)} samples")

    adatas = []
    for p in prefixes:
        sample, cond = sample_meta(p)
        try:
            a = load_mtx(p)
        except Exception as exc:  # noqa: BLE001
            print(f"  !! {sample}: {exc}")
            continue
        a.obs["sample"] = sample
        a.obs["condition"] = cond
        a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
        adatas.append(a)
        print(f"  {sample:30s} {cond:10s} {a.n_obs:6d} cells")

    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.var_names_make_unique()
    print(f"[concat] {adata.shape}")

    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    qc = adata.obs.groupby(["sample", "condition"])[
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
    print("\n[lineage composition]")
    print(adata.obs["lineage"].value_counts().to_string())

    present_key = [g for g in KEY_MARKERS if g in adata.raw.var_names]
    m = adata.raw[:, present_key].X
    m = np.asarray(m.todense()) if hasattr(m, "todense") else np.asarray(m)
    mdf = pd.DataFrame(m, columns=present_key, index=adata.obs_names)
    mdf["lineage"] = adata.obs["lineage"].values
    cm = mdf.groupby("lineage", observed=True).mean()
    cm.to_csv(os.path.join(tabdir, "lineage_key_markers.csv"), encoding="utf-8-sig")
    print("\n[markers] mean log-expression by lineage:")
    print(cm.round(2).to_string())

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

    # ---- candidate cell specificity (key output) ------------------------
    rows = []
    for g in CANDIDATES:
        if g not in adata.raw.var_names:
            continue
        vec = adata.raw[:, g].X
        vec = np.asarray(vec.todense()).ravel() if hasattr(vec, "todense") else np.asarray(vec).ravel()
        for lin in sorted(adata.obs["lineage"].unique()):
            sel = (adata.obs["lineage"] == lin).values
            if sel.sum() < 20:
                continue
            rows.append({"gene": g, "lineage": lin, "n_cells": int(sel.sum()),
                         "mean_log_expr": round(float(np.mean(vec[sel])), 3),
                         "pct_expressing": round(float(100.0 * np.mean(vec[sel] > 0)), 2)})
    spec = pd.DataFrame(rows)
    if len(spec):
        spec.to_csv(os.path.join(tabdir, "candidate_cell_specificity.csv"),
                    index=False, encoding="utf-8-sig")
        piv = spec.pivot(index="gene", columns="lineage", values="pct_expressing")
        print("\n[candidate cell specificity] % of cells expressing:")
        print(piv.round(1).to_string())
        rank = []
        for g in piv.index:
            v = piv.loc[g]
            ec = v.get("Endothelial", np.nan)
            others = v.drop(labels=[l for l in ["Endothelial", "LymphaticEC"] if l in v.index])
            rank.append({"gene": g, "pct_endothelial": ec,
                         "pct_max_other": float(others.max()) if len(others) else np.nan,
                         "top_other_lineage": others.idxmax() if len(others) else None,
                         "endothelial_is_top": bool(ec >= (others.max() if len(others) else np.inf))})
        pd.DataFrame(rank).to_csv(os.path.join(tabdir, "candidate_endothelial_ranking.csv"),
                                  index=False, encoding="utf-8-sig")
        print("\n[endothelial specificity ranking]")
        print(pd.DataFrame(rank).to_string(index=False))

    # ---- endothelial pseudobulk: Sham vs MI_day1 vs MI_day7 --------------
    mod_cols = [c for c in adata.obs.columns if c.startswith("mod_")]
    ec = adata[adata.obs["lineage"].isin(["Endothelial", "LymphaticEC"])].copy()
    print(f"\n[endothelium] {ec.n_obs} cells")
    if ec.n_obs > 0:
        cols = mod_cols + (["EGNS"] if "EGNS" in ec.obs else [])
        pb = ec.obs.groupby(["sample", "condition"], observed=True)[cols].mean()
        pb["n_cells"] = ec.obs.groupby(["sample", "condition"], observed=True).size()
        pb.to_csv(os.path.join(tabdir, "endothelial_pseudobulk_scores.csv"),
                  encoding="utf-8-sig")
        print(pb.round(3).to_string())

        from scipy import stats
        res = []
        for c in cols:
            sham = pb.loc[pb.index.get_level_values("condition") == "Sham", c].values
            d1 = pb.loc[pb.index.get_level_values("condition") == "MI_day1", c].values
            d7 = pb.loc[pb.index.get_level_values("condition") == "MI_day7", c].values
            if len(sham) >= 2 and len(d1) >= 2 and len(d7) >= 2:
                res.append({"score": c, "n": len(sham),
                            "sham": round(float(np.mean(sham)), 3),
                            "MI_day1": round(float(np.mean(d1)), 3),
                            "MI_day7": round(float(np.mean(d7)), 3),
                            "delta_d1": round(float(np.mean(d1) - np.mean(sham)), 3),
                            "delta_d7": round(float(np.mean(d7) - np.mean(sham)), 3),
                            "p_d1": round(float(stats.mannwhitneyu(d1, sham).pvalue), 4),
                            "p_d7": round(float(stats.mannwhitneyu(d7, sham).pvalue), 4)})
        pd.DataFrame(res).to_csv(os.path.join(tabdir, "endothelial_pseudobulk_stats.csv"),
                                 index=False, encoding="utf-8-sig")
        print("\n[endothelial pseudobulk stats]")
        print(pd.DataFrame(res).to_string(index=False))

    adata.obs[["sample", "condition", "cluster", "lineage"] + mod_cols +
              (["EGNS"] if "EGNS" in adata.obs else [])].to_csv(
        os.path.join(tabdir, "cell_scores.csv"), encoding="utf-8-sig")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sc.pl.umap(adata, color=["lineage", "condition"], show=False, wspace=0.4)
        plt.savefig(os.path.join(figdir, "umap_lineage_condition.png"), dpi=200,
                    bbox_inches="tight"); plt.close()
        cols = [c for c in ["mod_GCI", "mod_GDS", "mod_EIS", "EGNS"] if c in adata.obs]
        if cols:
            sc.pl.umap(adata, color=cols, show=False, cmap="viridis")
            plt.savefig(os.path.join(figdir, "umap_module_scores.png"), dpi=200,
                        bbox_inches="tight"); plt.close()
        print(f"[fig] {figdir}")
    except Exception as exc:  # noqa: BLE001
        print(f"[fig] skipped: {exc}")

    adata.write_h5ad(os.path.join(outdir, "gse176092_nonCM_processed.h5ad"))
    print(f"[ok] saved -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
