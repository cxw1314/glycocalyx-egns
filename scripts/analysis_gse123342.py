#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE123342 (human peripheral blood, n=192) glycocalyx proxy-score analysis.

Stage 1 --audit : dump the parsed sample metadata so grouping can be verified.
Stage 2 (default): compute GCI/GDS/EIS/EAS/EIDS module scores and EGNS,
                   run group statistics, and fit/evaluate a diagnostic model.

Usage:
    python analysis_gse123342.py --audit
    python analysis_gse123342.py
"""
from __future__ import annotations

import argparse
import json
import os
import pickle

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed", "GSE123342")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
MODULES = os.path.join(ROOT, "config", "modules.json")

# primary scores used for EGNS
PRIMARY = ["GCI", "GDS", "EIS", "EAS", "EIDS", "MFS", "FNS", "MRS"]


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        genes = [g.upper() for g in val.get("genes", [])]
        aliases = {k.upper(): [a.upper() for a in v]
                   for k, v in (val.get("aliases") or {}).items()}
        out[key] = {"genes": genes, "aliases": aliases, "score": val.get("score", key)}
    return out


def gene_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    """z-score each gene across samples (expression already log-scale)."""
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def module_scores(z: pd.DataFrame, modules: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores, coverage = {}, {}
    for key, spec in modules.items():
        present = [g for g in spec["genes"] if g in z.index]
        # alias rescue (e.g. CCN1 -> CYR61)
        for canon, alts in spec["aliases"].items():
            if canon not in present:
                for a in alts:
                    if a in z.index:
                        present.append(a)
                        break
        coverage[key] = len(present)
        scores[key] = z.loc[present].mean(axis=0) if present else np.nan
    return pd.DataFrame(scores), pd.DataFrame([coverage])


def add_egns(scores: pd.DataFrame) -> pd.DataFrame:
    def z(col):
        s = scores[col]
        sd = s.std(ddof=0)
        return (s - s.mean()) / sd if sd and not np.isnan(sd) else s * 0.0
    scores = scores.copy()
    scores["EGNS"] = z("GDS") + z("EAS") + z("EIDS") - z("GCI") - z("EIS")
    return scores


def _clean_series(s: pd.Series) -> pd.Series:
    return (s.astype(str).str.strip().str.strip('"').str.strip().str.strip("'"))


# analysis groups for GSE123342
GROUP_ORDER = ["StableCAD", "AcuteMI", "MI_30d", "MI_1y"]


def assign_groups(meta: pd.DataFrame) -> pd.Series:
    """Map samples using the parsed characteristics of GSE123342.

    ch_"patient diagnosis : MI | CAD
    ch_"time point        : D0 | D30 | Y1
    Sample_title          : ..._Techn_Repl_... marks a technical replicate
    """
    def col(name: str) -> pd.Series:
        if name in meta.columns:
            return _clean_series(meta[name])
        return pd.Series("", index=meta.index)

    diag = col('ch_"patient diagnosis').str.upper()
    time = col('ch_"time point').str.upper()
    title = col("Sample_title")

    grp = pd.Series("Unknown", index=meta.index, dtype=object)
    grp[title.str.contains("Techn_Repl", case=False, na=False)] = "Replicate"
    m = grp.eq("Unknown")
    grp[m & diag.eq("CAD")] = "StableCAD"
    grp[m & diag.eq("MI") & time.eq("D0")] = "AcuteMI"
    grp[m & diag.eq("MI") & time.eq("D30")] = "MI_30d"
    grp[m & diag.eq("MI") & time.eq("Y1")] = "MI_1y"
    return grp


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = np.sum(~np.isnan(p))
    order = np.argsort(np.where(np.isnan(p), np.inf, p))
    q = np.full_like(p, np.nan, dtype=float)
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        if np.isnan(p[idx]):
            continue
        val = min(prev, p[idx] * n / (n - rank + 1))
        q[idx] = val
        prev = val
    return list(q)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", action="store_true")
    args = ap.parse_args()
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)

    meta = pd.read_csv(os.path.join(PROC, "meta.csv"), index_col=0, encoding="utf-8-sig")
    # SOFT values arrive quote-wrapped (e.g. '"GSM3501241"') -> normalise everywhere
    meta.index = [str(i).strip().strip('"').strip("'") for i in meta.index]
    meta.columns = [str(c).strip().strip('"').strip("'") for c in meta.columns]
    for c in meta.columns:
        meta[c] = _clean_series(meta[c])
    if args.audit:
        print("=== meta columns ===")
        print(list(meta.columns))
        for c in meta.columns:
            vals = meta[c].astype(str).unique()
            print(f"\n-- {c} ({len(vals)} distinct)")
            for v in vals[:20]:
                print("   ", v)
        return 0

    with open(os.path.join(PROC, "expr.pkl"), "rb") as fh:
        expr = pickle.load(fh)
    expr.index = [str(i).upper() for i in expr.index]
    expr = expr[~expr.index.duplicated()]
    print(f"[expr] {expr.shape} (genes x samples)")

    modules = load_modules()
    z = gene_zscore(expr)
    scores, coverage = module_scores(z, modules)
    scores = add_egns(scores)
    scores.index.name = "sample"
    scores.to_csv(os.path.join(TAB, "GSE123342_module_scores.csv"), encoding="utf-8-sig")
    coverage.T.rename(columns={0: "n_genes_found"}).to_csv(
        os.path.join(TAB, "GSE123342_module_coverage.csv"), encoding="utf-8-sig")
    print("[scores] module coverage:")
    print(coverage.T.to_string())

    meta = meta.reindex(scores.index)
    meta["group"] = assign_groups(meta)
    meta.to_csv(os.path.join(TAB, "GSE123342_sample_groups.csv"), encoding="utf-8-sig")
    print("[groups]", meta["group"].value_counts().to_dict())

    # ---- group statistics per score -------------------------------------
    tested = PRIMARY + ["EGNS"]
    rows = []
    for score in tested:
        if score not in scores.columns:
            continue
        sub = pd.DataFrame({"v": scores[score], "g": meta["group"]})
        sub = sub[sub["g"].isin(GROUP_ORDER)]
        groups = [g["v"].values for _, g in sub.groupby("g") if len(g) > 2]
        if len(groups) < 2:
            continue
        kw = stats.kruskal(*groups)
        comps = [("AcuteMI", "StableCAD"), ("AcuteMI", "MI_30d"), ("AcuteMI", "MI_1y"),
                 ("MI_30d", "MI_1y"), ("StableCAD", "MI_1y")]
        prow = {"score": score, "kruskal_H": kw.statistic, "kruskal_p": kw.pvalue}
        pvals = []
        for a, b in comps:
            va = sub.loc[sub["g"] == a, "v"].values
            vb = sub.loc[sub["g"] == b, "v"].values
            if len(va) > 2 and len(vb) > 2:
                u = stats.mannwhitneyu(va, vb, alternative="two-sided")
                prow[f"{a}_vs_{b}_p"] = u.pvalue
                prow[f"{a}_vs_{b}_delta"] = float(np.median(va) - np.median(vb))
                pvals.append(u.pvalue)
            else:
                prow[f"{a}_vs_{b}_p"] = np.nan
                prow[f"{a}_vs_{b}_delta"] = np.nan
                pvals.append(np.nan)
        for (a, b), q in zip(comps, bh(pvals)):
            prow[f"{a}_vs_{b}_q"] = q
        rows.append(prow)
    stat_df = pd.DataFrame(rows)
    stat_df.to_csv(os.path.join(TAB, "GSE123342_group_stats.csv"), index=False,
                   encoding="utf-8-sig")
    show = [c for c in ["score", "kruskal_p", "AcuteMI_vs_StableCAD_p",
                        "AcuteMI_vs_StableCAD_delta"] if c in stat_df.columns]
    print(stat_df[show].to_string(index=False) if show else "[stats] no comparisons ran")

    # ---- diagnostic model (AcuteMI vs StableCAD) ------------------------
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline

    mask = meta["group"].isin(["AcuteMI", "StableCAD"])
    y = (meta.loc[mask, "group"] == "AcuteMI").astype(int).values
    X = scores.loc[mask]
    if len(np.unique(y)) == 2 and mask.sum() >= 20:
        rng = np.random.default_rng(42)
        results = []
        feature_sets = {"EGNS": ["EGNS"],
                        "GDS": ["GDS"],
                        "EGNS+GCI": ["EGNS", "GCI"],
                        "all_primary": PRIMARY}
        for name, cols in feature_sets.items():
            cols = [c for c in cols if c in X.columns]
            if not cols:
                continue
            clf = make_pipeline(StandardScaler(),
                                LogisticRegression(max_iter=2000, C=1.0))
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            prob = cross_val_predict(clf, X[cols].values, y, cv=cv,
                                     method="predict_proba")[:, 1]
            auc = roc_auc_score(y, prob)
            boots = []
            for _ in range(1000):
                idx = rng.integers(0, len(y), len(y))
                if len(np.unique(y[idx])) < 2:
                    continue
                boots.append(roc_auc_score(y[idx], prob[idx]))
            lo, hi = np.percentile(boots, [2.5, 97.5])
            fpr, tpr, _ = roc_curve(y, prob)
            results.append({"model": name, "n_features": len(cols), "auc_cv": auc,
                            "auc_ci_low": lo, "auc_ci_high": hi, "n": len(y),
                            "n_acute": int(y.sum())})
            np.save(os.path.join(TAB, f"roc_{name}.npy"),
                    np.vstack([fpr, tpr]))
            print(f"[model] {name}: AUC={auc:.3f} (95% CI {lo:.3f}-{hi:.3f})")
        pd.DataFrame(results).to_csv(
            os.path.join(TAB, "GSE123342_model_performance.csv"),
            index=False, encoding="utf-8-sig")
    else:
        print("[model] skipped: group labels or sample size insufficient")

    # ---- figures --------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        order = [g for g in GROUP_ORDER if g in set(meta["group"])]
        plot_df = pd.DataFrame(scores[["GDS", "GCI", "EGNS"]]).join(meta["group"])
        plot_df = plot_df[plot_df["group"].isin(order)]
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for ax, col in zip(axes, ["GDS", "GCI", "EGNS"]):
            sns.boxplot(data=plot_df, x="group", y=col, order=order, ax=ax,
                        color="#cbd5e1", fliersize=2)
            sns.stripplot(data=plot_df, x="group", y=col, order=order, ax=ax,
                          color="#0f172a", size=3, alpha=0.6)
            ax.set_title(col)
            ax.set_xlabel("")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, "GSE123342_scores_by_group.png"), dpi=200)
        plt.close()
        print(f"[fig] {os.path.join(FIG, 'GSE123342_scores_by_group.png')}")
    except Exception as exc:  # noqa: BLE001
        print(f"[fig] skipped: {exc}")

    print("[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
