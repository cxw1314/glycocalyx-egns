#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visium spatial neighbourhood analysis (GSE176092 WT samples).

Tests the central spatial hypothesis of the framework: are high-EGNS
endothelial-niche spots spatially adjacent to pro-fibrotic macrophage (MFS)
and activated fibroblast / ECM (FNS) niches, rather than merely co-expressed
with them in the same spot?

Method
  * per section: QC -> normalise -> module scores -> EGNS
  * k-nearest-neighbour graph on the Visium pixel coordinates (k=6)
  * same-spot correlation vs neighbour-mean correlation
  * top-quartile EGNS spots vs the rest: neighbour MFS/FNS/GDS, with a
    label-permutation test (1000 permutations) for the neighbour statistic
  * Moran's I spatial autocorrelation for EGNS and the niche scores
  * MI vs Sham comparison of spot-level EGNS

Usage:
    python analysis_visium_spatial.py --indir data/raw/GSE176092/visium \
                                      --outdir results/visium_spatial
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
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")
SEED = 42


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
    lookup: dict[str, str] = {}
    for v in var_names:
        lookup.setdefault(v.upper(), v)
    out = []
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


def condition_of(name: str) -> str:
    if "Sham" in name:
        return "Sham"
    m = re.search(r"day(\d+)", name)
    return f"MI_day{m.group(1)}" if m else "unknown"


def load_visium(prefix: str) -> ad.AnnData | None:
    tmp = tempfile.mkdtemp()
    try:
        for kind, pat, dst in (("barcodes", "_barcodes.tsv.gz", "barcodes.tsv.gz"),
                               ("features", "_features.tsv.gz", "features.tsv.gz"),
                               ("matrix", "_matrix.mtx.gz", "matrix.mtx.gz")):
            src = prefix + pat
            if not os.path.exists(src):
                cand = glob.glob(prefix + "*" + pat.split("_")[-1])
                if not cand:
                    return None
                src = cand[0]
            shutil.copy(src, os.path.join(tmp, dst))
        a = sc.read_10x_mtx(tmp, var_names="gene_symbols", make_unique=True)

        pos_file = None
        for cand in (prefix + "_tissue_positions_list.csv.gz",
                     prefix + "_tissue_positions.csv.gz"):
            if os.path.exists(cand):
                pos_file = cand
                break
        if pos_file is None:
            hits = glob.glob(prefix + "*tissue_positions*.csv.gz")
            pos_file = hits[0] if hits else None
        if pos_file is None:
            return None

        pos = pd.read_csv(pos_file, header=None)
        if pos.shape[1] >= 6:
            pos.columns = ["barcode", "in_tissue", "array_row", "array_col",
                           "pxl_row", "pxl_col"][:pos.shape[1]]
        else:
            pos = pd.read_csv(pos_file)
            pos = pos.rename(columns={"barcode": "barcode", "pxl_row_in_fullres": "pxl_row",
                                      "pxl_col_in_fullres": "pxl_col"})
        pos = pos.set_index("barcode")
        common = [b for b in a.obs_names if b in pos.index]
        a = a[common].copy()
        a.obsm["spatial"] = pos.loc[common, ["pxl_col", "pxl_row"]].to_numpy(dtype=float)
        return a
    except Exception:  # noqa: BLE001
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def knn_indices(coords: np.ndarray, k: int = 6) -> np.ndarray:
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(coords))).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)
    return idx[:, 1:]          # drop self


def morans_i(values: np.ndarray, idx: np.ndarray) -> float:
    x = values - values.mean()
    denom = np.sum(x ** 2)
    if denom <= 0:
        return np.nan
    num = 0.0
    for i, nb in enumerate(idx):
        num += np.sum(x[nb]) * x[i]
    n = len(x)
    w = idx.shape[1]
    return float((n / (n * w)) * (num / denom)) if w else np.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--min-counts", type=int, default=500)
    ap.add_argument("--n-perm", type=int, default=1000)
    args = ap.parse_args()

    indir = args.indir if os.path.isabs(args.indir) else os.path.join(ROOT, args.indir)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    prefixes = sorted(p[:-len("_matrix.mtx.gz")] for p in
                      glob.glob(os.path.join(indir, "*_matrix.mtx.gz")))
    if not prefixes:
        print(f"[error] no Visium matrices in {indir}")
        return 1
    print(f"[load] {len(prefixes)} Visium sections")

    modules = load_modules()
    rng = np.random.default_rng(SEED)
    per_section, summary = [], []

    for p in prefixes:
        name = os.path.basename(p)
        cond = condition_of(name)
        a = load_visium(p)
        if a is None:
            print(f"  !! {name}: could not load (matrix or positions missing)")
            continue
        sc.pp.filter_cells(a, min_counts=args.min_counts)
        if a.n_obs < 50:
            print(f"  !! {name}: only {a.n_obs} spots after QC, skipped")
            continue
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)

        for key, spec in modules.items():
            genes = resolve_genes(spec, a.var_names)
            if genes:
                sc.tl.score_genes(a, genes, score_name=f"mod_{key}")

        def zc(col):
            s = a.obs[col]
            sd = s.std(ddof=0)
            return (s - s.mean()) / sd if sd else s * 0.0

        needed = {"mod_GDS", "mod_EAS", "mod_EIDS", "mod_GCI", "mod_EIS"}
        if not needed <= set(a.obs.columns):
            print(f"  !! {name}: missing module scores")
            continue
        a.obs["EGNS"] = (zc("mod_GDS") + zc("mod_EAS") + zc("mod_EIDS")
                         - zc("mod_GCI") - zc("mod_EIS"))

        coords = a.obsm["spatial"]
        idx = knn_indices(coords, args.k)

        scores = a.obs[[c for c in a.obs.columns if c.startswith("mod_")] + ["EGNS"]]
        vals = {c: scores[c].to_numpy(float) for c in scores.columns}
        nb_mean = {c: np.array([np.mean(v[nb]) for nb in idx]) for c, v in vals.items()}

        e = vals["EGNS"]
        thr = np.quantile(e, 0.75)
        hi = e >= thr
        row = {"section": name, "condition": cond, "n_spots": int(a.n_obs)}
        for c in ["mod_MFS", "mod_FNS", "mod_GDS", "mod_EIS", "mod_GCI"]:
            if c not in vals:
                continue
            row[f"moranI_{c.replace('mod_','')}"] = round(morans_i(vals[c], idx), 4)
            row[f"rho_same_{c.replace('mod_','')}"] = round(
                float(stats.spearmanr(e, vals[c]).statistic), 4)
            row[f"rho_nbr_{c.replace('mod_','')}"] = round(
                float(stats.spearmanr(e, nb_mean[c]).statistic), 4)
            obs_diff = float(np.mean(nb_mean[c][hi]) - np.mean(nb_mean[c][~hi]))
            perm = []
            for _ in range(args.n_perm):
                lab = rng.permutation(hi)
                perm.append(np.mean(nb_mean[c][lab]) - np.mean(nb_mean[c][~lab]))
            row[f"d_nbr_{c.replace('mod_','')}"] = round(obs_diff, 4)
            row[f"p_nbr_{c.replace('mod_','')}"] = round(
                float((np.sum(np.abs(perm) >= abs(obs_diff)) + 1) / (args.n_perm + 1)), 4)
        row["moranI_EGNS"] = round(morans_i(e, idx), 4)
        summary.append(row)

        per_section.append(pd.DataFrame({
            "section": name, "condition": cond,
            **{c: vals[c] for c in vals},
            **{f"nbr_{c}": nb_mean[c] for c in nb_mean},
        }))
        print(f"  {name:28s} {cond:10s} {a.n_obs:5d} spots  "
              f"MoranI(EGNS)={row['moranI_EGNS']:.3f}  "
              f"rhoNbr(MFS)={row.get('rho_nbr_MFS', float('nan')):.3f}  "
              f"rhoNbr(FNS)={row.get('rho_nbr_FNS', float('nan')):.3f}")

    if not summary:
        print("[error] no section analysed")
        return 1
    summ = pd.DataFrame(summary)
    summ.to_csv(os.path.join(tabdir, "visium_spatial_summary.csv"), index=False,
                encoding="utf-8-sig")
    pd.concat(per_section, ignore_index=True).to_csv(
        os.path.join(tabdir, "visium_spot_scores.csv"), index=False, encoding="utf-8-sig")

    pd.set_option("display.width", 250)
    print("\n=== spatial neighbourhood summary ===")
    print(summ.to_string(index=False))

    print("\n=== across sections (mean +/- sd) ===")
    for col in [c for c in summ.columns if c.startswith(("moranI_", "rho_", "d_nbr_", "p_nbr_"))]:
        v = summ[col].astype(float)
        print(f"  {col:24s} {v.mean():+.4f} +/- {v.std():.4f}")

    print("\n[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
