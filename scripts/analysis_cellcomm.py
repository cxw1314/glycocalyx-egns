#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Cell-cell communication analysis (offline CellChat/LIANA equivalent).

Why offline: the R CellChat/NicheNet stack is unavailable here, and liana-py
fetches its ligand-receptor resource from the OmniPath web service, which is
not reachable over this network. We therefore use a curated LR panel relevant
to post-MI repair and implement the standard label-permutation test ourselves.

The framework's "triple condition" for a communication claim requires
  (1) ligand AND receptor expression,
  (2) a downstream response in the receptor cell,
  (3) spatial proximity.
(3) was established by the Visium analysis; this script addresses (1) and (2).

Usage:
    python analysis_cellcomm.py --h5ad results/GSE176092/gse176092_nonCM_processed.h5ad \
                               --outdir results/cellcomm
"""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 42

# Curated ligand -> receptor panel for post-MI repair
# (mouse symbols; complexes are given as lists of subunits, all required)
LR_PAIRS: dict[str, tuple[str, list[str]]] = {
    "CCL2-CCR2":        ("Ccl2", ["Ccr2"]),
    "CCL5-CCR5":        ("Ccl5", ["Ccr5"]),
    "CXCL1-CXCR2":      ("Cxcl1", ["Cxcr2"]),
    "CXCL12-CXCR4":     ("Cxcl12", ["Cxcr4"]),
    "SPP1-ITGAV_ITGB1": ("Spp1", ["Itgav", "Itgb1"]),
    "SPP1-CD44":        ("Spp1", ["Cd44"]),
    "TGFB1-TGFBR1_TGFBR2": ("Tgfb1", ["Tgfbr1", "Tgfbr2"]),
    "IL1B-IL1R1":       ("Il1b", ["Il1r1"]),
    "IL6-IL6R_IL6ST":   ("Il6", ["Il6r", "Il6st"]),
    "TNF-TNFRSF1A":     ("Tnf", ["Tnfrsf1a"]),
    "CSF1-CSF1R":       ("Csf1", ["Csf1r"]),
    "VEGFA-KDR":        ("Vegfa", ["Kdr"]),
    "VEGFA-FLT1":       ("Vegfa", ["Flt1"]),
    "ANGPT2-TEK":       ("Angpt2", ["Tek"]),
    "PDGFB-PDGFRB":     ("Pdgfb", ["Pdgfrb"]),
    "JAG1-NOTCH3":      ("Jag1", ["Notch3"]),
    "DLL4-NOTCH1":      ("Dll4", ["Notch1"]),
    "FN1-ITGB1":        ("Fn1", ["Itgb1"]),
    "POSTN-ITGAV_ITGB1": ("Postn", ["Itgav", "Itgb1"]),
    "LGALS3-ITGB1":     ("Lgals3", ["Itgb1"]),
    "THBS1-CD47":       ("Thbs1", ["Cd47"]),
    "MDK-NCL":          ("Mdk", ["Ncl"]),
    "HBEGF-EGFR":       ("Hbegf", ["Egfr"]),
    "GAS6-AXL":         ("Gas6", ["Axl"]),
    "IGF1-IGF1R":       ("Igf1", ["Igf1r"]),
    "BMP2-BMPR2":       ("Bmp2", ["Bmpr2"]),
}

# downstream programs to check in receptor cells (gene sets)
PROGRAMS = {
    "Macrophage_inflammatory": ["Il1b", "Tnf", "Ccl2", "Cxcl1", "Cxcl2", "Nfkb1", "Rel", "Il6"],
    "Macrophage_profibrotic": ["Spp1", "Trem2", "Gpnmb", "Lgals3", "Tgfb1", "Mmp14", "Cd9", "Fabp5"],
    "Fibroblast_activated": ["Postn", "Cthrc1", "Acta2", "Col1a1", "Col3a1", "Fn1", "Lox", "Tnc"],
    "Endothelial_activated": ["Vcam1", "Icam1", "Sele", "Ccl2", "Cxcl1", "Il6", "Rela"],
}

# which program to test for each receptor lineage
RECEIVER_PROGRAM = {
    "Macrophage": "Macrophage_profibrotic",
    "Fibroblast": "Fibroblast_activated",
    "Endothelial": "Endothelial_activated",
}


def resolve(symbol: str, lookup: dict[str, str]) -> str | None:
    return lookup.get(symbol.upper())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--min-pct", type=float, default=5.0)
    args = ap.parse_args()

    import anndata as ad
    path = args.h5ad if os.path.isabs(args.h5ad) else os.path.join(ROOT, args.h5ad)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    adata = ad.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    var_names = list(raw.var_names)
    lookup = {}
    for v in var_names:
        lookup.setdefault(v.upper(), v)
    print(f"[load] {adata.n_obs} cells x {len(var_names)} genes (raw)")

    lineage = adata.obs["lineage"].astype(str).values
    cond = adata.obs["condition"].astype(str).values
    samples = adata.obs["sample"].astype(str).values
    lineages = sorted(set(lineage))
    print(f"[lineages] {lineages}")

    X = raw.X
    dense_needed = False

    def gene_vec(sym: str) -> np.ndarray | None:
        real = resolve(sym, lookup)
        if real is None:
            return None
        j = var_names.index(real)
        col = X[:, j]
        return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()

    # expression weight per gene: %expressing x mean log-expression (bounded to [0,1] scale)
    weights: dict[str, np.ndarray] = {}
    for pair, (lig, recs) in LR_PAIRS.items():
        for g in [lig] + recs:
            if g in weights:
                continue
            v = gene_vec(g)
            if v is None:
                continue
            frac = np.mean(v > 0)
            weights[g] = (frac, v)

    # ---- lineages x lineage communication scores ------------------------
    rowsum = {}
    rows = []
    for pair, (lig, recs) in LR_PAIRS.items():
        if lig not in weights or any(r not in weights for r in recs):
            missing = [g for g in [lig] + recs if g not in weights]
            rows.append({"LR_pair": pair, "status": f"gene(s) absent: {','.join(missing)}"})
            continue
        lig_frac, lig_vec = weights[lig]
        rec_frac = float(np.mean([weights[r][0] for r in recs]))
        rec_vec = np.mean([weights[r][1] for r in recs], axis=0)
        if lig_frac < args.min_pct / 100 or rec_frac < args.min_pct / 100:
            rows.append({"LR_pair": pair, "status":
                         f"below {args.min_pct}% detection (lig {lig_frac*100:.1f}%, rec {rec_frac*100:.1f}%)"})
            continue

        # integer lineage codes -> vectorised group means (bincount), so the
        # permutation test can afford thousands of iterations
        codes = np.array([lineages.index(x) for x in lineage])
        G = len(lineages)

        def group_means(v: np.ndarray, lab_codes: np.ndarray) -> np.ndarray:
            c = np.bincount(lab_codes, minlength=G).astype(float)
            c[c == 0] = 1.0
            return np.bincount(lab_codes, weights=v, minlength=G) / c

        obs = np.outer(group_means(lig_vec, codes), group_means(rec_vec, codes))
        obs_flat = obs.ravel()

        # permutation null: shuffle lineage labels, keeping the cell counts
        rng = np.random.default_rng(SEED)
        null_flat = np.empty((args.n_perm, obs.size))
        for k in range(args.n_perm):
            pc = rng.permutation(codes)
            null_flat[k] = np.outer(group_means(lig_vec, pc),
                                    group_means(rec_vec, pc)).ravel()

        for i, a in enumerate(lineages):
            for j, b in enumerate(lineages):
                s = obs[i, j]
                if s <= 0:
                    continue
                idx = i * len(lineages) + j
                p = (np.sum(null_flat[:, idx] >= s) + 1) / (args.n_perm + 1)
                rows.append({"LR_pair": pair, "status": "ok",
                             "sender": a, "receiver": b,
                             "ligand": lig, "receptors": "+".join(recs),
                             "score": round(float(s), 4),
                             "p_perm": round(float(p), 4),
                             "lig_pct": round(lig_frac * 100, 2),
                             "rec_pct": round(rec_frac * 100, 2)})

    comm = pd.DataFrame(rows)
    comm.to_csv(os.path.join(tabdir, "cellcomm_all.csv"), index=False, encoding="utf-8-sig")

    ok = comm[comm["status"] == "ok"].copy()
    if len(ok):
        top = ok.sort_values("score", ascending=False).head(40)
        print("\n=== top lineage-to-lineage communications (all pairs) ===")
        print(top[["LR_pair", "sender", "receiver", "score", "p_perm",
                   "lig_pct", "rec_pct"]].to_string(index=False))
        ok.to_csv(os.path.join(tabdir, "cellcomm_significant.csv"), index=False,
                  encoding="utf-8-sig")

        # focus: which axes target endothelium, and which emanate from endothelium
        print("\n=== axes with Endothelial as SENDER (top 12) ===")
        es = ok[ok["sender"].str.contains("Endothelial")].sort_values("score", ascending=False).head(12)
        print(es[["LR_pair", "sender", "receiver", "score", "p_perm"]].to_string(index=False))
        print("\n=== axes with Endothelial as RECEIVER (top 12) ===")
        er = ok[ok["receiver"].str.contains("Endothelial")].sort_values("score", ascending=False).head(12)
        print(er[["LR_pair", "sender", "receiver", "score", "p_perm"]].to_string(index=False))

    # ---- condition (2): downstream response in receptor cells ------------
    resp = []
    for lin, prog in RECEIVER_PROGRAM.items():
        genes = [resolve(g, lookup) for g in PROGRAMS[prog]]
        genes = [g for g in genes if g is not None]
        if not genes:
            continue
        sel = lineage == lin
        if sel.sum() < 50:
            continue
        mat = np.vstack([weights[g.upper()][1] if False else
                         (np.asarray(X[:, var_names.index(g)].todense()).ravel()
                          if hasattr(X[:, var_names.index(g)], "todense")
                          else np.asarray(X[:, var_names.index(g)]).ravel())
                         for g in genes])
        prog_score = mat.mean(axis=0)
        # receptor genes for this lineage
        recs_here = sorted({r for _, recs in LR_PAIRS.values() for r in recs})
        for r in recs_here:
            real = resolve(r, lookup)
            if real is None:
                continue
            rv = (np.asarray(X[:, var_names.index(real)].todense()).ravel()
                  if hasattr(X[:, var_names.index(real)], "todense")
                  else np.asarray(X[:, var_names.index(real)]).ravel())
            if np.mean(rv[sel] > 0) < 0.05:
                continue
            rho = stats.spearmanr(rv[sel], prog_score[sel])
            resp.append({"lineage": lin, "receptor": r,
                         "program": prog,
                         "pct_expressing": round(float(100 * np.mean(rv[sel] > 0)), 2),
                         "rho_vs_program": round(float(rho.statistic), 4),
                         "p": float(f"{rho.pvalue:.3g}")})
    resp_df = pd.DataFrame(resp).sort_values("rho_vs_program", ascending=False)
    resp_df.to_csv(os.path.join(tabdir, "receptor_downstream_response.csv"),
                   index=False, encoding="utf-8-sig")
    print("\n=== condition (2): receptor expression vs downstream program ===")
    print(resp_df.to_string(index=False) if len(resp_df) else "(none)")

    # ---- time course of the niche programs (section 3.4 material) --------
    tc = []
    for lin, prog in RECEIVER_PROGRAM.items():
        genes = [resolve(g, lookup) for g in PROGRAMS[prog]]
        genes = [g for g in genes if g is not None]
        sel0 = lineage == lin
        if sel0.sum() < 50 or not genes:
            continue
        mat = np.vstack([(np.asarray(X[:, var_names.index(g)].todense()).ravel()
                          if hasattr(X[:, var_names.index(g)], "todense")
                          else np.asarray(X[:, var_names.index(g)]).ravel()) for g in genes])
        ps = mat.mean(axis=0)
        for c in ["Sham", "MI_day1", "MI_day7"]:
            for s in sorted(set(samples[sel0 & (cond == c)])):
                m = sel0 & (cond == c) & (samples == s)
                if m.sum() == 0:
                    continue
                tc.append({"lineage": lin, "program": prog, "condition": c,
                           "sample": s, "n_cells": int(m.sum()),
                           "program_score": round(float(ps[m].mean()), 4)})
    tc_df = pd.DataFrame(tc)
    tc_df.to_csv(os.path.join(tabdir, "niche_program_timecourse.csv"),
                 index=False, encoding="utf-8-sig")
    print("\n=== niche program time course (per sample) ===")
    if len(tc_df):
        piv = tc_df.pivot_table(index=["lineage", "program"], columns="condition",
                                values="program_score", aggfunc="mean")
        print(piv.round(4).to_string())

    print("\n[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
