#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE135310 -- cardiac immune time series after MI (design section 2.2 / 4.5 step 7).

Resolves the CCR2+ monocyte recruitment wave and the macrophage state
transition (inflammatory -> resolving -> SPP1/TREM2/GPNMB pro-fibrotic) across
steady state, day 1, day 3, day 5 and day 7, and computes the macrophage
fibrosis score (MFS) on the same axis.

IMPORTANT: each time point is a SINGLE animal, so the series is descriptive and
directional only; no inferential test is possible and none is reported.

Usage:
    python analysis_gse135310.py --indir data/raw/GSE135310/extracted \
                                 --outdir results/GSE135310
"""
from __future__ import annotations

import argparse
import glob
import json
import os
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

# file prefix -> (label, time point, sorted index)
SAMPLES = {
    "GSM4005123_F8":  ("F8",  "Steady state", 0),
    "GSM4005125_F11": ("F11", "Day 1", 1),
    "GSM4005126_F12": ("F12", "Day 3", 2),
    "GSM4005127_G1":  ("G1",  "Day 5", 3),
    "GSM4005124_F9":  ("F9",  "Day 7", 4),
}

LINEAGE_MARKERS = {
    "Monocyte":   ["Ly6c2", "Ccr2", "Plac8", "Chil3", "Ms4a6c", "Vcan"],
    "Macrophage": ["Adgre1", "Cd68", "C1qa", "C1qb", "C1qc", "Mrc1", "Cd163"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg", "Il1b"],
    "DC":         ["Cd209a", "Itgax", "H2-Ab1", "Flt3", "Batf3"],
    "Tcell":      ["Cd3e", "Cd3d", "Cd3g", "Lck"],
    "Bcell":      ["Cd79a", "Ms4a1", "Cd19"],
    "Nkcells":    ["Nkg7", "Klrb1c", "Gzma"],
    "Mast":       ["Cpa3", "Ms4a2", "Kit"],
    "pDC":        ["Siglech", "Bst2", "Irf7"],
}

# macrophage activation states along the repair axis
STATES = {
    "Mono_classical":  ["Ly6c2", "Ccr2", "Plac8", "Chil3", "Ms4a6c"],
    "Mac_inflammatory": ["Il1b", "Tnf", "Ccl2", "Cxcl1", "Cxcl2", "Nfkb1", "Nlrp3"],
    "Mac_resolving":   ["Mrc1", "Cd163", "Arg1", "Il10", "Tgfb2", "Cd36"],
    "Mac_profibrotic": ["Spp1", "Trem2", "Gpnmb", "Lgals3", "Tgfb1", "Mmp14", "Cd9", "Fabp5"],
    "Mac_antigenpresent": ["H2-Ab1", "H2-Aa", "Cd74", "Ly6c2"],
}


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        out[key] = sorted({g.upper() for g in val.get("genes", [])})
    return out


def load_mtx(prefix: str) -> ad.AnnData | None:
    """Robust 10x loader: the deposit mixes genes.tsv.gz and features.tsv.gz."""
    tmp = tempfile.mkdtemp()
    specs = [("barcodes.tsv.gz", ["_barcodes.tsv.gz"]),
             ("features.tsv.gz", ["_features.tsv.gz", "_genes.tsv.gz"]),
             ("matrix.mtx.gz", ["_matrix.mtx.gz"])]
    try:
        for dst, suffixes in specs:
            src = None
            for suf in suffixes:
                if os.path.exists(prefix + suf):
                    src = prefix + suf
                    break
            if src is None:
                for suf in suffixes:
                    hits = sorted(glob.glob(prefix + "*" + suf))
                    if hits:
                        src = hits[0]
                        break
            if src is None:
                print(f"    missing component for {os.path.basename(prefix)}")
                return None
            if dst == "features.tsv.gz":
                # scanpy expects three columns (id, symbol, feature_type) but this
                # deposit mixes Cell Ranger v2 genes.tsv.gz (2 columns) with v3
                # features.tsv.gz (3 columns); normalise to three.
                df = pd.read_csv(src, sep="\t", header=None, dtype=str)
                if df.shape[1] == 2:
                    df[2] = "Gene Expression"
                df.to_csv(os.path.join(tmp, dst), sep="\t", header=False, index=False)
            else:
                shutil.copy(src, os.path.join(tmp, dst))
        return sc.read_10x_mtx(tmp, var_names="gene_symbols", make_unique=True)
    except Exception as exc:  # noqa: BLE001
        print(f"    loader error: {exc}")
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def resolve(spec_genes: list[str], var_names) -> list[str]:
    lut: dict[str, str] = {}
    for v in var_names:
        lut.setdefault(v.upper(), v)
    return list(dict.fromkeys(lut[g] for g in spec_genes if g in lut))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=20.0)
    args = ap.parse_args()

    indir = args.indir if os.path.isabs(args.indir) else os.path.join(ROOT, args.indir)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adatas = []
    for stem, (label, tp, idx) in SAMPLES.items():
        pre = os.path.join(indir, stem + "_mm10")
        a = load_mtx(pre)
        if a is None:
            print(f"  !! cannot load {stem}")
            continue
        a.obs["sample"] = label
        a.obs["timepoint"] = tp
        a.obs["t_index"] = idx
        a.obs_names = [f"{label}_{b}" for b in a.obs_names]
        adatas.append(a)
        print(f"  {label:5s} {tp:13s} {a.n_obs:6d} cells")
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
    print(qc.round(1).to_string())

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
    print(f"[cluster] {adata.obs['cluster'].nunique()} clusters")

    for lin, gs in LINEAGE_MARKERS.items():
        present = resolve([g.upper() for g in gs], adata.raw.var_names)
        if present:
            sc.tl.score_genes(adata, present, score_name=f"lin_{lin}", use_raw=True)
    scols = [c for c in adata.obs.columns if c.startswith("lin_")]
    lab = adata.obs.groupby("cluster")[scols].mean()
    lab["lineage"] = lab.idxmax(axis=1).str.replace("lin_", "", regex=False)
    adata.obs["lineage"] = adata.obs["cluster"].map(lab["lineage"])
    lab.to_csv(os.path.join(tabdir, "cluster_lineage_scores.csv"), encoding="utf-8-sig")
    print("\n[lineage composition]")
    print(adata.obs["lineage"].value_counts().to_string())

    for key, gs in STATES.items():
        present = resolve([g.upper() for g in gs], adata.raw.var_names)
        if present:
            sc.tl.score_genes(adata, present, score_name=f"st_{key}", use_raw=True)

    modules = load_modules()
    for key, gs in modules.items():
        present = resolve(gs, adata.raw.var_names)
        if present:
            sc.tl.score_genes(adata, present, score_name=f"mod_{key}", use_raw=True)
            print(f"  [module] {key}: {len(present)} genes")

    stcols = [c for c in adata.obs.columns if c.startswith("st_")]
    modcols = [c for c in adata.obs.columns if c.startswith("mod_")]

    # ---- composition over time ------------------------------------------
    tp_order = ["Steady state", "Day 1", "Day 3", "Day 5", "Day 7"]
    comp = (adata.obs.groupby(["timepoint", "lineage"], observed=True).size()
            .unstack(fill_value=0).reindex(tp_order))
    frac = comp.div(comp.sum(axis=1), axis=0)
    comp.to_csv(os.path.join(tabdir, "lineage_counts_by_timepoint.csv"),
                encoding="utf-8-sig")
    frac.to_csv(os.path.join(tabdir, "lineage_fraction_by_timepoint.csv"),
                encoding="utf-8-sig")
    print("\n[lineage fraction by time point]")
    print((frac * 100).round(1).to_string())

    # ---- CCR2+ recruitment ----------------------------------------------
    ccr2 = "Ccr2" if "Ccr2" in adata.raw.var_names else None
    mono = adata.obs["lineage"].isin(["Monocyte", "Macrophage"]).values
    rows = []
    if ccr2:
        v = adata.raw[:, ccr2].X
        v = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.asarray(v).ravel()
        for tp in tp_order:
            m = adata.obs["timepoint"].values == tp
            rows.append({
                "timepoint": tp,
                "n_cells": int(m.sum()),
                "pct_Ccr2_pos_all": round(float(100 * np.mean(v[m] > 0)), 2),
                "pct_Ccr2_pos_myeloid": round(float(100 * np.mean(v[m & mono] > 0)), 2),
                "mean_Ccr2_myeloid": round(float(np.mean(v[m & mono])), 4),
            })
        pd.DataFrame(rows).to_csv(os.path.join(tabdir, "ccr2_recruitment.csv"),
                                  index=False, encoding="utf-8-sig")
        print("\n[CCR2+ recruitment]")
        print(pd.DataFrame(rows).to_string(index=False))

    # ---- macrophage state trajectory ------------------------------------
    mac = adata.obs["lineage"] == "Macrophage"
    rows = []
    for tp in tp_order:
        m = (adata.obs["timepoint"] == tp).values
        r = {"timepoint": tp, "n_macrophage": int((m & mac).sum())}
        for c in stcols:
            r[c.replace("st_", "")] = round(float(adata.obs.loc[m & mac, c].mean()), 4)
        for c in modcols:
            r[c.replace("mod_", "")] = round(float(adata.obs.loc[m & mac, c].mean()), 4)
        rows.append(r)
    traj = pd.DataFrame(rows)
    traj.to_csv(os.path.join(tabdir, "macrophage_state_trajectory.csv"),
                index=False, encoding="utf-8-sig")
    print("\n[macrophage state trajectory]")
    print(traj.to_string(index=False))

    adata.obs[["sample", "timepoint", "t_index", "cluster", "lineage"] +
              stcols + modcols].to_csv(os.path.join(tabdir, "cell_scores.csv"),
                                       encoding="utf-8-sig")
    adata.write_h5ad(os.path.join(outdir, "gse135310_processed.h5ad"))
    print(f"\n[ok] saved -> {outdir}")
    print("[NOTE] one animal per time point: descriptive only, no inference.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
