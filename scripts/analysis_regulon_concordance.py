#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Recompute the concordance between the v1.1 curated-target proxy regulons and
the resource-backed CollecTRI/ULM regulons, using the already-saved
regulon_activity_cells.csv (no CollecTRI re-download).

Fixes the case-sensitivity bug that made the Pearson column NaN: the resource
returns mouse symbols (Klf2) while the panel is written in upper case (KLF2).
"""
from __future__ import annotations

import json
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")

HOME = {
    "KLF2": "EIS", "KLF4": "EIS", "ERG": "EIS", "FLI1": "EIS", "SOX17": "EIS",
    "SOX18": "EIS", "EPAS1": "EIS", "PROX1": "EIS", "RELA": "EAS",
    "NFKB1": "EAS", "NFKB2": "EAS", "RELB": "EAS", "YAP1": "MRS",
    "WWTR1": "MRS", "TEAD1": "MRS", "TEAD2": "MRS", "TEAD3": "MRS",
    "TEAD4": "MRS", "STAT3": "MFS", "HIF1A": "MFS",
}


def main() -> int:
    import anndata as ad

    cell_csv = os.path.join(ROOT, "results", "decoupler_regulons", "tables",
                            "regulon_activity_cells.csv")
    acts = pd.read_csv(cell_csv, index_col=0)
    acts["lineage"] = acts["lineage"].astype(str)
    endo = acts[acts["lineage"].str.contains("Endothelial")]
    colmap = {c.upper(): c for c in endo.columns
              if c not in ("lineage", "condition", "sample")}
    print(f"[regulons available] {len(colmap)}; endothelium {len(endo)} cells")

    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    modules = {k: [g.upper() for g in v.get("genes", [])]
               for k, v in cfg.items()
               if isinstance(v, dict) and "genes" in v}

    adata = ad.read_h5ad(os.path.join(ROOT, "results", "GSE176092",
                                      "gse176092_nonCM_processed.h5ad"))
    raw = adata.raw if adata.raw is not None else adata
    lut = {v.upper(): v for v in raw.var_names}
    mask = adata.obs["lineage"].astype(str).str.contains("Endothelial").values
    cond = adata.obs["condition"].values[mask]

    rows = []
    for tf, mod in HOME.items():
        tgt = sorted({tf, *modules.get(mod, [])})
        present = [lut[g] for g in tgt if g in lut]
        sub = raw[mask][:, present].X
        sub = np.asarray(sub.todense()) if hasattr(sub, "todense") else np.asarray(sub)
        m = np.log1p(sub / np.maximum(sub.sum(axis=1, keepdims=True), 1e-9) * 1e4).mean(axis=1)
        by_cond = pd.Series(m).groupby(cond).mean()
        col = colmap.get(tf.upper())
        r = np.nan
        if col is not None:
            r = float(np.corrcoef(m, endo[col].values)[0, 1])
        rows.append({"tf": tf,
                     "n_curated_targets": len(present),
                     "curated_D1_minus_Sham": float(by_cond.get("MI_day1", np.nan)
                                                    - by_cond.get("Sham", np.nan)),
                     "collecTRI_regulon_present": col is not None,
                     "pearson_r_curated_vs_collecTRI": r})
    out = pd.DataFrame(rows)
    dest = os.path.join(ROOT, "results", "decoupler_regulons", "tables",
                        "collecTRI_vs_curated_concordance.csv")
    out.to_csv(dest, index=False, encoding="utf-8-sig")
    print("\n=== curated-target proxy vs CollecTRI/ULM (endothelium) ===")
    print(out.round(4).to_string(index=False))
    ok = out["pearson_r_curated_vs_collecTRI"].dropna()
    if len(ok):
        print(f"\n[concordance] median r = {ok.median():.3f} "
              f"(range {ok.min():.2f} to {ok.max():.2f}; n = {len(ok)})")
    print(f"\n[ok] written to {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
