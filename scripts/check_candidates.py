#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Report GSE123342 differential-expression statistics for the prespecified
candidate sheddases / biosynthesis genes (HPSE and alternatives)."""
from __future__ import annotations

import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")

CANDIDATES = ["HPSE", "HPSE2", "HYAL1", "HYAL2", "HYAL3", "MMP2", "MMP9", "MMP14",
              "ADAM17", "ADAMTS1", "ADAMTS5", "SULF1", "SULF2", "NEU1", "PLAUR",
              "GPC1", "EXT1", "EXT2", "EXTL3", "NDST1", "HS2ST1", "HS6ST1",
              "HAS1", "HAS2", "HAS3", "SDC1", "SDC2", "SDC4", "CSGALNACT2"]


def main() -> int:
    deg = pd.read_csv(os.path.join(TAB, "GSE123342_panel_DE.csv"), encoding="utf-8-sig")
    sel = pd.read_csv(os.path.join(TAB, "GSE123342_feature_selection.csv"),
                      encoding="utf-8-sig")
    merged = deg.merge(sel[["gene", "in_LASSO", "in_RF", "in_SVM_RFE", "n_methods"]],
                       on="gene", how="left")
    sub = merged[merged["gene"].isin(CANDIDATES)].copy()
    sub["direction"] = sub["delta"].apply(lambda d: "up in acute MI" if d > 0 else "down")
    sub = sub.sort_values("p")
    pd.set_option("display.width", 200)
    print("=== candidate molecules in GSE123342 (AcuteMI vs StableCAD) ===")
    print(sub[["gene", "median_acuteMI", "median_stableCAD", "delta", "p", "q",
               "n_methods", "direction"]].to_string(index=False))
    sub.to_csv(os.path.join(TAB, "GSE123342_candidate_genes.csv"),
               index=False, encoding="utf-8-sig")
    print()
    print("=== which candidates were selected by >=1 ML method ===")
    hit = sub[sub["n_methods"].fillna(0) >= 1]
    print(hit[["gene", "n_methods", "q"]].to_string(index=False) if len(hit)
          else "(none)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
