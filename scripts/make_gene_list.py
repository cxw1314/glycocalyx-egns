#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the gene list to extract from a very large genes x cells matrix:
all predefined module genes plus the lineage marker panel."""
from __future__ import annotations

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")
OUT = os.path.join(ROOT, "config", "gene_list.txt")

LINEAGE_MARKERS = [
    # endothelial
    "PECAM1", "CDH5", "EMCN", "KDR", "TEK", "CLDN5", "GJA5",
    # pericyte / SMC
    "PDGFRB", "RGS5", "KCNJ8", "ABCC9", "NOTCH3", "MYH11", "ACTA2", "TAGLN",
    # fibroblast
    "DCN", "PDGFRA", "COL1A1", "POSTN", "CTHRC1", "LUM",
    # myeloid
    "CD68", "PTPRC", "CD163", "CCR2", "LYZ", "AIF1", "ITGAM", "C1QA",
    # lymphoid
    "CD3E", "CD3D", "CD79A", "MS4A1", "NKG7",
    # lymphatic / adipocyte / neuronal
    "PROX1", "LYVE1", "PDPN", "ADIPOQ", "PLIN1", "PLP1",
    # general
    "VIM", "COL1A2", "FN1", "SPP1", "TGFB1", "IL6", "CCL2",
]


def main() -> int:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    genes: list[str] = []
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        genes.extend(g.upper() for g in val.get("genes", []))
        for canon, alts in (val.get("aliases") or {}).items():
            genes.append(canon.upper())
            genes.extend(a.upper() for a in alts)
    genes.extend(LINEAGE_MARKERS)
    seen = list(dict.fromkeys(genes))
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(seen) + "\n")
    print(f"[ok] {len(seen)} genes -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
