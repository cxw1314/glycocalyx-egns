#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Print the virtual-knockout results table, grouped by cell type."""
from __future__ import annotations

import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> int:
    p = os.path.join(ROOT, "results", "virtual_ko", "tables", "virtual_ko_results.csv")
    if not os.path.exists(p):
        print("[missing]", p)
        return 1
    d = pd.read_csv(p, encoding="utf-8-sig")
    pd.set_option("display.width", 240)
    print(f"[rows] {len(d)}  columns: {list(d.columns)}\n")
    for ct in d["cell_type"].unique():
        sub = d[d["cell_type"] == ct].sort_values("delta_readout")
        cols = [c for c in sub.columns if c != "cell_type"]
        print(f"=== {ct} ===")
        print(sub[cols].to_string(index=False))
        print()
    print("=== smallest p per cell type ===")
    for ct in d["cell_type"].unique():
        sub = d[d["cell_type"] == ct].nsmallest(5, "p_empirical")
        print(f"-- {ct}")
        print(sub[["ko_tf", "delta_readout", "z_vs_null", "p_empirical"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
