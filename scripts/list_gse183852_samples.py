#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""List GSE183852 sample titles and disease states for barcode-prefix reconciliation."""
from __future__ import annotations

import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "data", "processed", "GSE183852", "meta.csv")


def main() -> int:
    d = pd.read_csv(META, index_col=0, encoding="utf-8-sig")
    d.columns = [str(c).strip().strip('"') for c in d.columns]
    for c in d.columns:
        d[c] = d[c].astype(str).str.strip().str.strip('"')
    dis_col = next(c for c in d.columns if "disease" in c.lower())
    title_col = "Sample_title" if "Sample_title" in d.columns else d.columns[0]
    out = d[[title_col, dis_col]].copy()
    out.columns = ["sample_title", "disease_state"]
    print(out.to_string())
    print()
    print("titles starting with H / T W C M / H D C M:")
    for t in out["sample_title"]:
        if t.upper().startswith(("HDCM", "H_", "TWCM-H")):
            print("   ", t)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
