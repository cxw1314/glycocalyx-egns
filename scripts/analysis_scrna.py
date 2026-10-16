#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mouse MI non-cardiomyocyte single-cell analysis (GSE176092) for the EGNS framework.

Pipeline: load 10x matrices -> QC -> normalise -> integrate -> cluster ->
annotate lineages -> module scores (GCI/GDS/EIS/EAS/EIDS/EGNS) on endothelium ->
sample-level pseudobulk differential testing (animal as the statistical unit).

Usage:
    python analysis_scrna.py --indir data/raw/GSE176092/nonCM --outdir results/GSE176092
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
import re
import shutil
import tempfile

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")

# marker panels for lineage annotation (mouse symbols)
LINEAGE_MARKERS = {
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Tek", "Emcn"],
    "Fibroblast": ["Dcn", "Pdgfra", "Col1a1", "Postn", "Cthrc1"],
    "Macrophage": ["Cd68", "Adgre1", "Csf1r", "Lyz2", "Ccr2"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Kcnj8", "Abcc9", "Notch3"],
    "Tcell": ["Cd3e", "Cd3d", "Cd3g", "Lck"],
    "Bcell": ["Cd79a", "Ms4a1", "Cd19"],
    "Nkcells": ["Nkg7", "Klrb1c", "Gzma"],
    "Adipocyte": ["Adipoq", "Plin1", "Fabp4"],
    "Neuronal": ["Plp1", "Sox10", "Ngf"],
    "LymphaticEC": ["Prox1", "Lyve1", "Pdpn"],
}


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        genes = list(val.get("genes", []))
        aliases = {k: list(v) for k, v in (val.get("aliases") or {}).items()}
        out[key] = {"genes": genes, "aliases": aliases}
    return out


def sample_label(fname: str) -> str:
    m = re.search(r"(Sham|MI_day\d+|MI_d\d+)", fname)
    base = re.sub(r"^GSM\d+_", "", fname)
    base = re.sub(r"_(barcodes|features|genes|matrix)\..*$", "", base)
    return base


def load_10x_dir(prefix: str) -> ad.AnnData:
    """prefix like .../GSM5355657_nonCM_Sham_1"""
    tmp = tempfile.mkdtemp()
    try:
        for kind, pats in (("barcodes", ["*barcodes.tsv.gz"]),
                           ("features", ["*features.tsv.gz", "*genes.tsv.gz"]),
                           ("matrix", ["*matrix.mtx.gz"])):
            src = None
            for pat in pats:
                hits = glob.glob(prefix + pat.split("*")[-1])
                if hits:
                    src = hits[0]
                    break
            if src is None:
                raise FileNotFoundError(f"missing {kind} for {prefix}")
            shutil.copy(src, os.path.join(tmp, f"{kind}.tsv.gz" if kind != "matrix"
                                          else "matrix.mtx.gz"))
        a = sc.read_10x_mtx(tmp, var_names="gene_symbols", make_unique=True)
        return a
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
    os.makedirs(outdir, exist_ok=True)
    tabdir = os.path.join(outdir, "tables")
    figdir = os.path.join(outdir, "figures")
    os.makedirs(tabdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)

    prefixes = sorted({re.sub(r"_(barcodes|features|genes|matrix)\.(tsv|mtx)\.gz$", "", p)
                       for p in glob.glob(os.path.join(indir, "*matrix.mtx.gz"))})
    if not prefixes:
        print(f"[error] no 10x matrix files in {indir}")
        return 1
    print(f"[load] {len(prefixes)} samples")

    adatas = []
    for p in prefixes:
        lab = sample_label(os.path.basename(p))
        try:
            a = load_10x_dir(p)
        except Exception as exc:  # noqa: BLE001
            print(f"  !! {lab}: {exc}")
            continue
        a.obs["sample"] = lab
        a.obs_names = [f"{lab}_{b}" for b in a.obs_names]
        adatas.append(a)
        print(f"  {lab}: {a.n_obs} cells x {a.n_vars} genes")
    if not adatas:
        print("[error] nothing loaded")
        return 1

    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.var_names_make_unique()
    print(f"[concat] {adata.shape}")

    # ---- QC (per-sample reporting, one rule applied transparently) ------
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    qc = adata.obs.groupby("sample")[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].median()
    qc.to_csv(os.path.join(tabdir, "qc_medians_by_sample.csv"), encoding="utf-8-sig")
    print(qc.to_string())

    before = adata.n_obs
    adata = adata[(adata.obs["n_genes_by_counts"] >= args.min_genes) &
                  (adata.obs["pct_counts_mt"] <= args.max_mito)].copy()
    print(f"[qc] {before} -> {adata.n_obs} cells")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
    sc.external.pp.harmony_integrate(adata, "sample") if False else None
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.6, key_added="cluster", flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(adata)
    print(f"[cluster] {adata.obs['cluster'].nunique()} clusters")

    # ---- lineage annotation --------------------------------------------
    for lin, genes in LINEAGE_MARKERS.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{lin}",
                              use_raw=True)
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    lab = adata.obs.groupby("cluster")[score_cols].mean()
    lab["lineage"] = lab.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["lineage"] = adata.obs["cluster"].map(lab["lineage"]).astype(str)
    lab.to_csv(os.path.join(tabdir, "cluster_lineage_scores.csv"), encoding="utf-8-sig")
    print(lab["lineage"].value_counts().to_string())

    # ---- module scores --------------------------------------------------
    modules = load_modules()
    for key, spec in modules.items():
        present = [g for g in spec["genes"] if g in adata.raw.var_names]
        for canon, alts in spec["aliases"].items():
            if canon not in present:
                for a in alts:
                    if a in adata.raw.var_names:
                        present.append(a)
                        break
        if present:
            sc.tl.score_genes(adata, present, score_name=f"mod_{key}", use_raw=True)
    mod_cols = [c for c in adata.obs.columns if c.startswith("mod_")]

    def zc(col):
        s = adata.obs[col]
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd else s * 0.0
    if {"mod_GDS", "mod_EAS", "mod_EIDS", "mod_GCI", "mod_EIS"} <= set(mod_cols):
        adata.obs["EGNS"] = (zc("mod_GDS") + zc("mod_EAS") + zc("mod_EIDS")
                             - zc("mod_GCI") - zc("mod_EIS"))

    adata.obs[["sample", "cluster", "lineage"] + mod_cols +
              (["EGNS"] if "EGNS" in adata.obs else [])].to_csv(
        os.path.join(tabdir, "cell_scores.csv"), encoding="utf-8-sig")

    # ---- endothelial pseudobulk (animal as the unit) --------------------
    ec = adata[adata.obs["lineage"].isin(["Endothelial", "LymphaticEC"])].copy()
    print(f"[endothelium] {ec.n_obs} cells")
    if ec.n_obs > 0 and mod_cols:
        pb = ec.obs.groupby("sample")[mod_cols + (["EGNS"] if "EGNS" in ec.obs else [])].mean()
        pb["n_cells"] = ec.obs.groupby("sample").size()
        pb["condition"] = [re.sub(r"_?\d+$", "", s) for s in pb.index]
        pb.to_csv(os.path.join(tabdir, "endothelial_pseudobulk_scores.csv"),
                  encoding="utf-8-sig")
        print(pb.to_string())

        from scipy import stats
        res = []
        for col in mod_cols + (["EGNS"] if "EGNS" in ec.obs else []):
            sham = pb.loc[pb["condition"].str.contains("Sham", na=False), col].values
            mi = pb.loc[~pb["condition"].str.contains("Sham", na=False), col].values
            if len(sham) >= 2 and len(mi) >= 2:
                u = stats.mannwhitneyu(mi, sham, alternative="two-sided")
                res.append({"score": col, "n_sham": len(sham), "n_mi": len(mi),
                            "median_sham": float(np.median(sham)),
                            "median_mi": float(np.median(mi)),
                            "delta": float(np.median(mi) - np.median(sham)),
                            "p": u.pvalue})
        pd.DataFrame(res).to_csv(os.path.join(tabdir, "endothelial_pseudobulk_stats.csv"),
                                 index=False, encoding="utf-8-sig")
        print(pd.DataFrame(res).to_string(index=False))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sc.pl.umap(adata, color=["lineage", "cluster"], show=False, wspace=0.4)
        plt.savefig(os.path.join(figdir, "umap_lineage_cluster.png"), dpi=200,
                    bbox_inches="tight")
        plt.close()
        cols = [c for c in ["mod_GCI", "mod_GDS", "mod_EIS", "EGNS"] if c in adata.obs]
        if cols:
            sc.pl.umap(adata, color=cols, show=False, cmap="viridis")
            plt.savefig(os.path.join(figdir, "umap_module_scores.png"), dpi=200,
                        bbox_inches="tight")
            plt.close()
        print(f"[fig] written to {figdir}")
    except Exception as exc:  # noqa: BLE001
        print(f"[fig] skipped: {exc}")

    adata.write_h5ad(os.path.join(outdir, "nonCM_processed.h5ad"))
    print(f"[ok] saved to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
