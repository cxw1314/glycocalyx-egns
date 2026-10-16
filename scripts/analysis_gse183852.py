#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE183852 cross-species consistency check (human heart failure).

Question (design section 3.6): does reduced glycocalyx integrity / loss of
endothelial identity persist into end-stage human heart failure?

Scope caveat enforced in the code and in the output: this cohort is
DILATED CARDIOMYOPATHY, not post-MI heart failure, so it can only serve as a
consistency check in end-stage failure -- never as causal evidence that MI
causes heart failure.

Inputs:
  data/processed/GSE183852/counts_selected.pkl  (genes x cells, streamed subset)
  data/processed/GSE183852/meta.csv             (sample -> disease state)

Usage:
    python analysis_gse183852.py --indir data/processed/GSE183852 \
                                 --outdir results/GSE183852
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import re
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")
SEED = 42

LINEAGE_MARKERS = {
    "Endothelial": ["PECAM1", "CDH5", "EMCN", "KDR", "TEK", "CLDN5", "GJA5"],
    "Pericyte": ["PDGFRB", "RGS5", "KCNJ8", "ABCC9", "NOTCH3"],
    "SMC": ["MYH11", "ACTA2", "TAGLN"],
    "Fibroblast": ["DCN", "PDGFRA", "COL1A1", "POSTN", "CTHRC1", "LUM"],
    "Macrophage": ["CD68", "CD163", "LYZ", "AIF1", "ITGAM", "C1QA"],
    "Tcell": ["CD3E", "CD3D"],
    "Bcell": ["CD79A", "MS4A1"],
    "Nkcells": ["NKG7"],
    "LymphaticEC": ["PROX1", "LYVE1", "PDPN"],
    "Adipocyte": ["ADIPOQ", "PLIN1"],
    "Neuronal": ["PLP1"],
}


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        genes = [g.upper() for g in val.get("genes", [])]
        for canon, alts in (val.get("aliases") or {}).items():
            genes.append(canon.upper())
            genes.extend(a.upper() for a in alts)
        out[key] = sorted(set(genes))
    return out


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    idx = np.argsort(p[ok])
    n = int(ok.sum())
    prev = 1.0
    out = np.empty(n)
    for rank, i in enumerate(idx[::-1], start=1):
        prev = min(prev, p[ok][i] * n / (n - rank + 1))
        out[i] = prev
    q[ok] = out
    return q


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    indir = args.indir if os.path.isabs(args.indir) else os.path.join(ROOT, args.indir)
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(tabdir, exist_ok=True)

    with open(os.path.join(indir, "counts_selected.pkl"), "rb") as fh:
        counts = pickle.load(fh)
    counts.index = [str(g).upper() for g in counts.index]
    counts = counts[~counts.index.duplicated()]
    print(f"[counts] {counts.shape} (genes x cells)")

    meta = pd.read_csv(os.path.join(indir, "meta.csv"), index_col=0, encoding="utf-8-sig")
    meta.columns = [str(c).strip().strip('"') for c in meta.columns]
    for c in meta.columns:
        meta[c] = meta[c].astype(str).str.strip().str.strip('"')

    title_col = "Sample_title" if "Sample_title" in meta.columns else meta.columns[0]
    dis_col = next((c for c in meta.columns if "disease" in c.lower()), None)
    tech_col = next((c for c in meta.columns if "technology" in c.lower()), None)
    if dis_col is None:
        print("[error] no disease-state column in metadata")
        return 1
    samp2dis = dict(zip(meta[title_col], meta[dis_col]))
    samp2tech = dict(zip(meta[title_col], meta[tech_col])) if tech_col else {}
    print(f"[meta] {len(samp2dis)} samples; disease states: "
          f"{pd.Series(list(samp2dis.values())).value_counts().to_dict()}")

    # robust sample matching: barcode prefixes in this matrix do not always
    # equal the Sample_title verbatim (e.g. 'H_ZC-LVAD' vs 'H_ZC-LVAD-1')
    def norm(s: str) -> str:
        return re.sub(r"[^A-Z0-9]", "", str(s).upper())

    title_norm = {norm(t): t for t in meta[title_col]}

    # The seven single-cell libraries carry barcode prefixes HDCM<n>, whereas
    # the GEO Sample_title for the same libraries is TWCM-H<n> (optionally with
    # a -lib suffix). Reconcile them explicitly so those donors are not lost.
    hdcm_alias: dict[str, str] = {}
    for t in meta[title_col]:
        m = re.fullmatch(r"TWCM-H(\d+)(?:-lib\d+)?", str(t).strip(), flags=re.I)
        if m:
            hdcm_alias[f"HDCM{m.group(1)}"] = t
    if hdcm_alias:
        print(f"[map] HDCM->title alias ({len(hdcm_alias)}): {hdcm_alias}")

    def match_sample(s: str) -> str | None:
        n = norm(s)
        if not n:
            return None
        if s.strip().upper() in hdcm_alias:
            return hdcm_alias[s.strip().upper()]
        if n in title_norm:
            return title_norm[n]
        best = None
        for tn, t in title_norm.items():
            if tn.startswith(n) or n.startswith(tn):
                if best is None or len(tn) > len(best[0]):
                    best = (tn, t)
        return best[1] if best else None

    # barcode -> sample (everything before the LAST underscore)
    barcodes = list(counts.columns)
    raw_samples = [b.rsplit("_", 1)[0] for b in barcodes]
    matched_titles = [match_sample(s) for s in raw_samples]
    matched = sum(1 for t in matched_titles if t is not None)
    print(f"[map] {matched}/{len(barcodes)} cells mapped to a known sample")
    missing = sorted({s for s, t in zip(raw_samples, matched_titles) if t is None})[:8]
    if missing:
        print(f"[map] unmatched sample prefixes: {missing}")

    cell_dis = np.array([samp2dis.get(t, "unmapped") for t in matched_titles])
    cell_tech = np.array([samp2tech.get(t, "unknown") for t in matched_titles])
    cell_samp = np.array([t if t is not None else s
                          for t, s in zip(matched_titles, raw_samples)])
    print(f"[map] cells by disease: {pd.Series(cell_dis).value_counts().to_dict()}")

    # ---- normalise ------------------------------------------------------
    X = counts.to_numpy(dtype=np.float32)
    lib = X.sum(axis=0)
    lib[lib == 0] = 1.0
    Xn = np.log1p(X / lib * 1e4)

    genes = list(counts.index)
    gidx = {g: i for i, g in enumerate(genes)}

    # ---- lineage annotation by marker score ------------------------------
    # NOTE: only a 165-gene subset was streamed out of the full matrix, so the
    # QC threshold must be relative to that subset, not to a whole-transcriptome
    # expectation.
    keep_mask = np.ones(Xn.shape[1], dtype=bool)
    min_genes = 25
    ngene = (X > 0).sum(axis=0)
    keep_mask &= ngene >= min_genes
    print(f"[qc] {int(keep_mask.sum())} / {len(keep_mask)} cells kept "
          f"(>= {min_genes} of {X.shape[0]} panel genes detected)")

    lin_scores = {}
    for lin, markers in LINEAGE_MARKERS.items():
        present = [gidx[g] for g in markers if g in gidx]
        if not present:
            continue
        lin_scores[lin] = Xn[present].mean(axis=0)
    ldf = pd.DataFrame(lin_scores, index=counts.columns)
    ldf = ldf[keep_mask]
    lineage = ldf.idxmax(axis=1).values
    print("[lineage composition]")
    print(pd.Series(lineage).value_counts().to_string())

    # ---- module scores --------------------------------------------------
    modules = load_modules()
    mrows = []
    mscores = {}
    for key, genes_wanted in modules.items():
        present = [gidx[g] for g in genes_wanted if g in gidx]
        mrows.append({"module": key, "n_genes_found": len(present),
                      "of": len(genes_wanted)})
        if present:
            mscores[key] = Xn[present][:, keep_mask].mean(axis=0)
    pd.DataFrame(mrows).to_csv(os.path.join(tabdir, "module_coverage.csv"),
                               index=False, encoding="utf-8-sig")
    print("\n[module coverage]")
    print(pd.DataFrame(mrows).to_string(index=False))

    sdf = pd.DataFrame(mscores, index=counts.columns[keep_mask])
    keys = list(mscores.keys())

    def zc(col):
        sd = sdf[col].std(ddof=0)
        return (sdf[col] - sdf[col].mean()) / sd if sd else sdf[col] * 0.0

    if {"GDS", "EAS", "EIDS", "GCI", "EIS"} <= set(keys):
        sdf["EGNS"] = (zc("GDS") + zc("EAS") + zc("EIDS") - zc("GCI") - zc("EIS"))
    sdf["lineage"] = lineage
    sdf["disease"] = cell_dis[keep_mask]
    sdf["technology"] = cell_tech[keep_mask]
    sdf["sample"] = cell_samp[keep_mask]
    sdf.to_csv(os.path.join(tabdir, "cell_scores.csv"), encoding="utf-8-sig")

    # ---- endothelial pseudobulk: DCM vs non-diseased ---------------------
    ec = sdf[sdf["lineage"] == "Endothelial"].copy()
    print(f"\n[endothelium] {len(ec)} cells across {ec['sample'].nunique()} samples")
    if len(ec) < 100:
        print("[error] too few endothelial cells")
        return 1

    score_cols = keys + (["EGNS"] if "EGNS" in sdf.columns else [])
    pb = ec.groupby(["sample", "disease", "technology"], observed=True)[score_cols].mean()
    pb["n_cells"] = ec.groupby(["sample", "disease", "technology"], observed=True).size()
    pb.to_csv(os.path.join(tabdir, "endothelial_pseudobulk.csv"), encoding="utf-8-sig")
    print(pb.round(4).to_string())

    res = []
    for c in score_cols:
        for tech in ["all"] + sorted(set(pb.index.get_level_values("technology"))):
            sub = pb if tech == "all" else pb[pb.index.get_level_values("technology") == tech]
            dcm = sub.loc[sub.index.get_level_values("disease")
                          .str.contains("Dilated", case=False), c].values
            ctl = sub.loc[~sub.index.get_level_values("disease")
                          .str.contains("Dilated", case=False), c].values
            if len(dcm) < 3 or len(ctl) < 3:
                continue
            u = stats.mannwhitneyu(dcm, ctl, alternative="two-sided")
            res.append({"score": c, "technology": tech,
                        "n_DCM": len(dcm), "n_control": len(ctl),
                        "mean_DCM": round(float(dcm.mean()), 4),
                        "mean_control": round(float(ctl.mean()), 4),
                        "delta": round(float(dcm.mean() - ctl.mean()), 4),
                        "p": u.pvalue})
    rdf = pd.DataFrame(res)
    if len(rdf):
        rdf["q_all"] = np.nan
        m = rdf["technology"] == "all"
        rdf.loc[m, "q_all"] = bh(rdf.loc[m, "p"].values)
        rdf.to_csv(os.path.join(tabdir, "endothelial_DCM_vs_control.csv"),
                   index=False, encoding="utf-8-sig")
        print("\n=== endothelial DCM vs non-diseased (pseudobulk, donor as unit) ===")
        print(rdf.to_string(index=False))

    print("\n[NOTE] this cohort is dilated cardiomyopathy, not post-MI heart failure:")
    print("       results are a consistency check only and cannot establish that")
    print("       myocardial infarction causes heart failure.")
    print("\n[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
