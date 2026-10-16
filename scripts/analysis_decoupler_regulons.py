#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Real regulon-activity inference with decoupler (replaces the earlier hand-built
"mean of 6-13 curated target genes" proxy).

Two inference layers are produced so the manuscript can state exactly what was
done:
  (A) CollecTRI / DoRothEA regulons scored by the ULM (univariate linear model)
      method of decoupler -- a published, resource-backed regulon activity.
  (B) the curated-target mean used in v1.1, recomputed here for a head-to-head
      concordance check.

Output:
    results/decoupler_regulons/tables/regulon_activity_cells.csv
    results/decoupler_regulons/tables/regulon_activity_pseudobulk.csv
    results/decoupler_regulons/tables/regulon_source_audit.csv
    results/decoupler_regulons/tables/collecTRI_vs_curated_concordance.csv
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

PANEL_TFS = ["KLF2", "KLF4", "ERG", "FLI1", "SOX17", "SOX18", "EPAS1", "PROX1",
             "RELA", "NFKB1", "NFKB2", "RELB", "YAP1", "WWTR1", "TEAD1",
             "TEAD2", "TEAD3", "TEAD4", "STAT3", "HIF1A", "SMAD1", "SMAD3",
             "JUN", "JUNB", "FOS", "FOSL2", "MEF2C", "PPARG", "GATA6", "NFATC1"]


def curated_targets(modules: dict) -> dict[str, list[str]]:
    """Regulons approximated in v1.1: a TF plus its home module's genes."""
    home = {
        "KLF2": "EIS", "KLF4": "EIS", "ERG": "EIS", "FLI1": "EIS",
        "SOX17": "EIS", "SOX18": "EIS", "EPAS1": "EIS", "PROX1": "EIS",
        "RELA": "EAS", "NFKB1": "EAS", "NFKB2": "EAS", "RELB": "EAS",
        "YAP1": "MRS", "WWTR1": "MRS", "TEAD1": "MRS", "TEAD2": "MRS",
        "TEAD3": "MRS", "TEAD4": "MRS", "STAT3": "MFS", "HIF1A": "MFS",
    }
    out = {}
    for tf, mod in home.items():
        gs = list(modules.get(mod, []))
        if gs:
            out[tf] = sorted({tf, *gs})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default="results/decoupler_regulons")
    ap.add_argument("--organism", default="mouse")
    ap.add_argument("--net", default="collectri", choices=["collectri", "dorothea"])
    args = ap.parse_args()

    import anndata as ad
    import decoupler as dc

    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    print(f"[input] {adata.n_obs} x {adata.n_vars}")

    # ---- (A) resource-backed regulons ----------------------------------
    if args.net == "collectri":
        net = dc.op.collectri(organism=args.organism)
    else:
        net = dc.op.dorothea(organism=args.organism, levels=["A", "B", "C"])
    net.columns = [str(c).lower() for c in net.columns]
    print(f"[{args.net}] columns={list(net.columns)}  head symbols="
          f"{list(net['source'].astype(str)[:5])}")

    # Map resource symbols onto the data's own gene symbols, case-insensitively,
    # so a mouse/human capitalisation mismatch cannot silently drop edges.
    upper2data = {}
    for v in raw.var_names:
        upper2data.setdefault(str(v).upper(), v)
    net["source"] = [upper2data.get(str(s).upper(), str(s)) for s in net["source"]]
    net["target"] = [upper2data.get(str(s).upper(), str(s)) for s in net["target"]]
    tgt_cov = len(set(net["target"]) & set(raw.var_names)) / max(1, len(set(net["target"])))
    src_cov = len(set(net["source"]) & set(raw.var_names)) / max(1, len(set(net["source"])))
    print(f"[{args.net}] {net['source'].nunique()} TFs, {len(net)} edges; "
          f"target coverage {tgt_cov:.3f}, source coverage {src_cov:.3f}")
    net.to_csv(os.path.join(tabdir, "regulon_source_audit.csv"),
               index=False, encoding="utf-8-sig")
    tfs_avail = set(net["source"].astype(str))
    panel_hit = []
    for t in PANEL_TFS:
        cand = t if t in tfs_avail else upper2data.get(t.upper(), None)
        if cand is not None and cand in tfs_avail:
            panel_hit.append(cand)
    print(f"[panel] {len(panel_hit)}/{len(PANEL_TFS)} panel TFs present in resource")
    audit = pd.DataFrame({
        "tf": PANEL_TFS,
        "in_resource": [t.upper() in {str(x).upper() for x in tfs_avail}
                        for t in PANEL_TFS],
        "n_targets_in_resource": [int((net["source"].astype(str).str.upper()
                                       == t.upper()).sum()) for t in PANEL_TFS],
    })
    audit.to_csv(os.path.join(tabdir, "panel_tf_coverage.csv"),
                 index=False, encoding="utf-8-sig")

    # counts -> CPM log for the linear model
    M = raw.X
    M = np.asarray(M.todense()) if hasattr(M, "todense") else np.asarray(M)
    df = pd.DataFrame(M, index=raw.obs_names, columns=raw.var_names)
    df = np.log1p(df.div(df.sum(axis=1).replace(0, np.nan), axis=0) * 1e4).fillna(0)

    # decoupler 2.2 wants a net with columns exactly source/target/weight; the
    # column names are NOT passed as kwargs to ulm.
    for col in ("source", "target", "weight"):
        if col not in net.columns:
            raise ValueError(f"regulon table is missing column '{col}': {list(net.columns)}")
    A = ad.AnnData(
        X=df.values.astype(np.float32),
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=df.columns.astype(str)),
    )
    dc.mt.ulm(A, net[["source", "target", "weight"]], verbose=False)
    keys = list(A.obsm.keys())
    print(f"[decoupler] obsm keys = {keys}")
    score_key = next((k for k in keys if k.startswith("score")), keys[0])
    raw_score = A.obsm[score_key]
    if isinstance(raw_score, pd.DataFrame):
        est = raw_score.copy()
        est.index = A.obs_names
    else:
        arr = np.asarray(raw_score)
        srcs = None
        for k in list(A.uns.keys()) + list(A.obsm.keys()):
            v = A.uns[k] if k in A.uns else A.obsm[k]
            if isinstance(v, (pd.Index,)) or (
                    isinstance(v, (list, np.ndarray)) and not np.isscalar(v)
                    and len(v) == arr.shape[1]):
                try:
                    if len(v) == arr.shape[1]:
                        srcs = [str(x) for x in v]
                        break
                except Exception:
                    continue
        if srcs is None:
            srcs = [f"tf_{i}" for i in range(arr.shape[1])]
        print(f"[decoupler] {arr.shape[1]} regulons; source names from "
              f"{'obsm/uns' if srcs[0] != 'tf_0' else 'fallback'}")
        est = pd.DataFrame(arr, index=A.obs_names, columns=srcs)
    est["lineage"] = A.obs["lineage"].values
    est["condition"] = A.obs["condition"].values
    est["sample"] = A.obs["sample"].values
    keep_tf = [c for c in est.columns if c in panel_hit]
    if not keep_tf:
        print("[warn] no panel TF names matched the regulon columns; "
              "keeping all regulons")
        keep_tf = [c for c in est.columns
                   if c not in ("lineage", "condition", "sample")]
    est = est[keep_tf + ["lineage", "condition", "sample"]]
    est.to_csv(os.path.join(tabdir, "regulon_activity_cells.csv"),
               encoding="utf-8-sig")
    print(f"\n[ULM] {est.shape[1] - 3} regulon activities for {len(est)} cells")

    endo = est[est["lineage"].astype(str).str.contains("Endothelial")]
    print(f"[endothelium] {len(endo)} cells")
    agg = endo.groupby("condition")[[c for c in panel_hit]].mean()
    agg.to_csv(os.path.join(tabdir, "regulon_activity_pseudobulk.csv"),
               encoding="utf-8-sig")
    print("\n=== endothelial regulon activity (ULM) by condition ===")
    print(agg.round(3).to_string())

    # ---- (B) curated-target mean, for concordance ------------------------
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    modules = {k: [g.upper() for g in v.get("genes", [])]
               for k, v in cfg.items()
               if isinstance(v, dict) and "genes" in v}
    cur = curated_targets(modules)
    lut = {v.upper(): v for v in raw.var_names}
    rows = []
    for tf, tgts in cur.items():
        present = [lut[g] for g in tgts if g in lut]
        if not present:
            continue
        sub = raw[endo.index][:, present].X
        sub = np.asarray(sub.todense()) if hasattr(sub, "todense") else np.asarray(sub)
        m = np.log1p(sub / np.maximum(sub.sum(axis=1, keepdims=True), 1e-9) * 1e4).mean(axis=1)
        s = pd.Series(m, index=endo.index)
        by_cond = s.groupby(endo["condition"]).mean()
        if tf in endo.columns and set(by_cond.index) >= {"Sham", "MI_day1"}:
            r = np.corrcoef(s.values, endo[tf].values)[0, 1]
        else:
            r = np.nan
        rows.append({"tf": tf, "n_curated_targets": len(present),
                     "curated_D1_minus_Sham": float(by_cond.get("MI_day1", np.nan)
                                                    - by_cond.get("Sham", np.nan)),
                     "pearson_r_curated_vs_collecTRI": float(r) if r == r else np.nan})
    conc = pd.DataFrame(rows)
    conc.to_csv(os.path.join(tabdir, "collecTRI_vs_curated_concordance.csv"),
                index=False, encoding="utf-8-sig")
    print("\n=== curated-target mean vs CollecTRI/ULM (endothelium) ===")
    print(conc.round(4).to_string(index=False))
    ok = conc["pearson_r_curated_vs_collecTRI"].dropna()
    if len(ok):
        print(f"\n[concordance] median r = {ok.median():.3f} "
              f"(n = {len(ok)} TFs)")
    print("\n[ok] regulon analysis written to " + tabdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
