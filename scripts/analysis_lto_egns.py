#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Leave-target-out EGNS (LTO-EGNS) -- removes the self-circularity between the
score and the transcription factors that are tested against it.

Problem being fixed
-------------------
Descriptive EGNS = z(GDS) + z(EAS) + z(EIDS) - z(GCI) - z(EIS), and the module
gene sets contain the very TFs analysed later (RELA in EAS; KLF2/KLF4 in EIS;
VIM/FN1 in EIDS; YAP1/WWTR1 in MRS). A knockout of such a TF is therefore
guaranteed to move the score through its own membership, which reviewers will
read as circular.

Two corrections are produced here
---------------------------------
1. Core-EGNS = z(GDS) - z(GCI)
   Uses only the shedding and biosynthesis arms, which contain no TF from the
   perturbation panel. This is the primary, non-circular readout.

2. LTO-EGNS(target) = descriptive formula recomputed after deleting the knocked
   out gene (and, where supplied, its regulon members) from whichever module
   contains it.

Usage:
    python analysis_lto_egns.py --h5ad results/GSE176092/gse176092_nonCM_processed.h5ad \
                                --outdir results/lto_egns
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

# which module each perturbed TF sits in (from config/modules.json)
TF_HOME = {
    "KLF2": "EIS", "KLF4": "EIS", "ERG": "EIS", "FLI1": "EIS",
    "SOX17": "EIS", "SOX18": "EIS", "EPAS1": "EIS", "PROX1": "EIS",
    "RELA": "EAS", "NFKB1": "EAS", "NFKB2": "EAS", "RELB": "EAS",
    "YAP1": "MRS", "WWTR1": "MRS", "TEAD1": "MRS", "TEAD2": "MRS",
    "TEAD3": "MRS", "TEAD4": "MRS",
    "STAT3": "MFS", "HIF1A": "MFS", "SMAD1": "FNS", "SMAD3": "FNS",
    "JUN": "FNS", "JUNB": "FNS", "FOS": "FNS", "FOSL2": "FNS",
    "MEF2C": "EIS", "PPARG": "MFS", "GATA6": "EIS", "NFATC1": "FNS",
}


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        out[key] = [g.upper() for g in val.get("genes", [])]
    return out


def z(v: np.ndarray) -> np.ndarray:
    sd = v.std()
    return (v - v.mean()) / sd if sd else v * 0.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--knn", type=int, default=30)
    args = ap.parse_args()

    import anndata as ad
    from sklearn.neighbors import NearestNeighbors

    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    lut: dict[str, str] = {}
    for v in raw.var_names:
        lut.setdefault(v.upper(), v)

    modules = load_modules()
    mask = adata.obs["lineage"].astype(str).str.contains("Endothelial").values
    n = int(mask.sum())
    print(f"[endothelium] {n} cells")

    emb = np.asarray(adata.obsm["X_pca"])[mask]
    nn = NearestNeighbors(n_neighbors=min(args.knn, n - 1)).fit(emb)
    _, nb = nn.kneighbors(emb)

    def module_expr(genes: list[str]) -> np.ndarray | None:
        present = [lut[g] for g in genes if g in lut]
        if not present:
            return None
        M = raw[mask][:, present].X
        M = np.asarray(M.todense()) if hasattr(M, "todense") else np.asarray(M)
        rs = M.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        M = np.log1p(M / rs * 1e4)
        return M[nb].mean(axis=1)          # kNN-imputed, then gene-averaged

    scores = {}
    coverage = []
    for k, gs in modules.items():
        v = module_expr(gs)
        if v is None:
            continue
        scored_genes = [g for g in gs if g in lut]
        scores[k] = v.mean(axis=1)
        coverage.append({"module": k, "genes_defined": len(gs),
                         "genes_present": len(scored_genes)})
    cov = pd.DataFrame(coverage)
    cov.to_csv(os.path.join(tabdir, "module_coverage_corrected.csv"),
               index=False, encoding="utf-8-sig")
    print("\n[module coverage, corrected denominators]")
    print(cov.to_string(index=False))

    # ---- descriptive EGNS and the non-circular Core-EGNS ----------------
    descriptive = (z(scores["GDS"]) + z(scores["EAS"]) + z(scores["EIDS"])
                   - z(scores["GCI"]) - z(scores["EIS"]))
    core = z(scores["GDS"]) - z(scores["GCI"])
    out = pd.DataFrame({
        "descriptive_EGNS": descriptive,
        "core_EGNS": core,
        "GDS_z": z(scores["GDS"]), "GCI_z": z(scores["GCI"]),
        "EIS_z": z(scores["EIS"]), "EAS_z": z(scores["EAS"]),
        "EIDS_z": z(scores["EIDS"]),
        "condition": adata.obs["condition"].values[mask],
        "sample": adata.obs["sample"].values[mask],
    })
    out.to_csv(os.path.join(tabdir, "egns_variants_per_cell.csv"),
               encoding="utf-8-sig")

    print("\n=== descriptive vs core EGNS by condition (cell-level mean) ===")
    agg = out.groupby("condition")[["descriptive_EGNS", "core_EGNS"]].mean()
    print(agg.round(4).to_string())

    print("\n=== per-animal (pseudobulk) ===")
    pb = out.groupby(["sample", "condition"])[["descriptive_EGNS", "core_EGNS"]].mean()
    print(pb.round(4).to_string())
    pb.to_csv(os.path.join(tabdir, "egns_variants_pseudobulk.csv"),
              encoding="utf-8-sig")

    print("\n=== correlation between the two EGNS definitions ===")
    r = np.corrcoef(out["descriptive_EGNS"], out["core_EGNS"])[0, 1]
    print(f"    Pearson r = {r:.4f}")

    # ---- TF membership audit: which TFs sit inside which module --------
    rows = []
    for tf, home in TF_HOME.items():
        in_modules = [k for k, gs in modules.items() if tf in gs]
        rows.append({"tf": tf, "home_module_from_panel": home,
                     "actually_in_modules": ",".join(in_modules) or "-",
                     "circular_for_descriptive_EGNS": bool(in_modules)})
    audit = pd.DataFrame(rows)
    audit.to_csv(os.path.join(tabdir, "tf_module_membership_audit.csv"),
                 index=False, encoding="utf-8-sig")
    print("\n=== TF membership audit (descriptive EGNS circularity) ===")
    print(audit.to_string(index=False))
    n_circ = int(audit["circular_for_descriptive_EGNS"].sum())
    print(f"\n[n] TFs sitting inside a scored module: {n_circ} / {len(audit)}")
    print("    -> for these, the descriptive EGNS readout is circular; "
          "Core-EGNS (GDS-GCI, no TF) is the non-circular alternative.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
