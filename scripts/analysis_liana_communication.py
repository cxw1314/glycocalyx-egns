#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Real cell-cell communication with LIANA (replaces the earlier hand-rolled
"expression fraction x mean expression" heuristic).

Runs the LIANA consensus rank-aggregate pipeline, which internally scores each
ligand-receptor pair with multiple methods (CellPhoneDB, NATMI, Connectome,
log2FC, SingleCellSignalR, CellChat-equivalent) and aggregates them into a
single consensus rank. All methods use OmniPath resources.

Output:
    results/liana_communication/tables/liana_consensus_all.csv
    results/liana_communication/tables/liana_endothelial_receiver.csv
    results/liana_communication/tables/method_agreement.csv
"""
from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", default="results/liana_communication")
    ap.add_argument("--groupby", default="lineage")
    ap.add_argument("--n-perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--reuse", action="store_true",
                    help="reuse an existing liana_consensus_all.csv")
    args = ap.parse_args()

    import scanpy as sc
    import liana as li

    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = sc.read_h5ad(path)
    print(f"[input] {adata.n_obs} cells x {adata.n_vars} genes")
    print("[clusters] " + ", ".join(
        f"{k}={v}" for k, v in adata.obs[args.groupby].value_counts().items()))

    # LIANA needs raw counts in .X
    if adata.raw is not None:
        adata = adata.raw.to_adata()
    print(f"[counts] {adata.n_obs} x {adata.n_vars}")

    cache = os.path.join(tabdir, "liana_consensus_all.csv")
    if args.reuse and os.path.exists(cache):
        res = pd.read_csv(cache)
        print(f"[reuse] loaded {len(res)} interactions from {cache}")
    else:
        li.mt.rank_aggregate(
            adata,
            groupby=args.groupby,
            resource_name="mouseconsensus",
            expr_prop=0.1,
            use_raw=False,
            verbose=False,
            n_perms=args.n_perms,
            seed=args.seed,
        )
        res = adata.uns["liana_res"].copy()
    print(f"[raw columns] {list(res.columns)}")
    # LIANA 1.10 does not always emit aggregate_rank; by definition the
    # consensus rank is the mean of the magnitude- and specificity-based
    # aggregates (LIANA's own aggregation), so reconstruct it when absent.
    if "aggregate_rank" not in res.columns:
        if {"magnitude_rank", "specificity_rank"} <= set(res.columns):
            res["aggregate_rank"] = res[["magnitude_rank", "specificity_rank"]].mean(axis=1)
            print("[info] aggregate_rank = mean(magnitude_rank, specificity_rank)")
        else:
            num = [c for c in res.columns if pd.api.types.is_numeric_dtype(res[c])]
            res["aggregate_rank"] = res[num].rank().mean(axis=1)
    if "aggregate_pval" not in res.columns:
        if "cellphone_pvals" in res.columns:
            res["aggregate_pval"] = res["cellphone_pvals"]
            print("[info] aggregate_pval <- CellPhoneDB permutation p-values")
        else:
            res["aggregate_pval"] = np.nan
    res.to_csv(os.path.join(tabdir, "liana_consensus_all.csv"),
               index=False, encoding="utf-8-sig")
    print(f"\n[consensus] {len(res)} ligand-receptor interactions scored")

    # ---- how well the individual methods agree (audit trail) -----------
    cand = [c for c in res.columns
            if c in ("lr_means", "expr_prod", "scaled_weight", "lr_logfc",
                     "spec_weight", "lrscore", "magnitude_rank", "specificity_rank")]
    agree = []
    for c in cand:
        sub = res[["aggregate_rank", c]].dropna()
        if len(sub) < 10:
            continue
        rho = sub["aggregate_rank"].corr(sub[c], method="spearman")
        agree.append({"method": c, "n": len(sub),
                      "spearman_vs_aggregate": round(float(rho), 4)})
    agr = pd.DataFrame(agree).sort_values("spearman_vs_aggregate", ascending=False)
    agr.to_csv(os.path.join(tabdir, "method_agreement.csv"),
               index=False, encoding="utf-8-sig")
    print("\n=== individual methods vs. LIANA consensus (Spearman rho) ===")
    print(agr.to_string(index=False))

    # ---- endothelial-receiver interactions ------------------------------
    src = "source" if "source" in res.columns else "ligand_complex"
    rec = "target" if "target" in res.columns else "receptor_complex"
    ec = res[res[rec].astype(str).str.contains("Endothelial", case=False, na=False)].copy()
    keep = [c for c in ["source", "target", "ligand_complex", "receptor_complex",
                        "aggregate_rank", "aggregate_pval", "lr_means",
                        "specificity_rank", "magnitude_rank", "lrscore"] if c in ec.columns]
    ec = ec[keep].sort_values("aggregate_rank")
    ec.to_csv(os.path.join(tabdir, "liana_endothelial_receiver.csv"),
              index=False, encoding="utf-8-sig")
    print(f"\n=== top 20 interactions with an Endothelial receiver "
          f"(of {len(ec)}) ===")
    print(ec.head(20).to_string(index=False))

    # ---- glycan-shedding / ECM axis specifically ------------------------
    sheddase = res[
        res["ligand_complex"].astype(str).str.contains(
            "FN1|MMP|HPSE|PLAUR|PLAU|ADAM|TIMP|SPP1|TGFB|IL1|TNF|CCL2|CXCL",
            case=False, na=False)
    ]
    if len(sheddase):
        sheddase = sheddase.sort_values("aggregate_rank")
        sheddase.to_csv(os.path.join(tabdir, "liana_sheddase_ecm_axis.csv"),
                        index=False, encoding="utf-8-sig")
        print("\n=== sheddase / ECM-axis interactions (top 20) ===")
        cols = [c for c in ["source", "target", "ligand_complex",
                            "receptor_complex", "aggregate_rank",
                            "aggregate_pval"] if c in sheddase.columns]
        print(sheddase[cols].head(20).to_string(index=False))

    print("\n[ok] LIANA consensus communication written to " + tabdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
