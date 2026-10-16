#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE123342 -- glycocalyx/niche gene-panel machine learning.

Implements the prespecified three-way feature screening from the study design
(LASSO, random forest, SVM-RFE), builds a cross-validated diagnostic model and
compares it with the module-score EGNS.

Panel   : union of the 8 predefined modules' genes present in the matrix.
Groups  : AcuteMI (n=65) vs StableCAD (n=22); technical replicates excluded.

Usage:
    python analysis_gse123342_ml.py
"""
from __future__ import annotations

import json
import os
import pickle
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed", "GSE123342")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
MODULES = os.path.join(ROOT, "config", "modules.json")
SEED = 42


def load_modules() -> dict:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = {}
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        genes = [g.upper() for g in val.get("genes", [])]
        for canon, alts in (val.get("aliases") or {}).items():
            genes.extend(a.upper() for a in alts)
        out[key] = genes
    return out


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    idx = np.argsort(p[ok])
    n = ok.sum()
    prev = 1.0
    ordered = np.empty(n)
    for rank, i in enumerate(idx[::-1], start=1):
        prev = min(prev, p[ok][i] * n / (n - rank + 1))
        ordered[i] = prev
    q[ok] = ordered
    return q


def main() -> int:
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    rng = np.random.default_rng(SEED)

    with open(os.path.join(PROC, "expr.pkl"), "rb") as fh:
        expr = pickle.load(fh)
    expr.index = [str(i).upper() for i in expr.index]
    expr = expr[~expr.index.duplicated()]

    groups = pd.read_csv(os.path.join(TAB, "GSE123342_sample_groups.csv"),
                         index_col=0, encoding="utf-8-sig")
    gcol = "group" if "group" in groups.columns else groups.columns[-1]
    grp = groups[gcol]
    mask = grp.isin(["AcuteMI", "StableCAD"])
    y = (grp[mask] == "AcuteMI").astype(int).values
    samples = grp[mask].index
    print(f"[samples] {mask.sum()} (AcuteMI={int(y.sum())}, StableCAD={int((1-y).sum())})")

    modules = load_modules()
    panel = sorted({g for genes in modules.values() for g in genes if g in expr.index})
    print(f"[panel] {len(panel)} module genes present in the matrix")

    X = expr.loc[panel, samples].T          # samples x genes
    X = X.loc[:, X.std(axis=0) > 0]

    # ---- per-gene differential expression (AcuteMI vs StableCAD) --------
    rows = []
    for g in X.columns:
        a = X.loc[y == 1, g].values
        b = X.loc[y == 0, g].values
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"gene": g, "median_acuteMI": float(np.median(a)),
                     "median_stableCAD": float(np.median(b)),
                     "delta": float(np.median(a) - np.median(b)), "p": u.pvalue})
    deg = pd.DataFrame(rows)
    deg["q"] = bh(deg["p"].values)
    deg = deg.sort_values("p")
    deg.to_csv(os.path.join(TAB, "GSE123342_panel_DE.csv"), index=False,
               encoding="utf-8-sig")
    sig = deg[deg["q"] < 0.05]
    print(f"[DE] {len(sig)} / {len(deg)} panel genes significant at q<0.05")
    print(deg.head(12).to_string(index=False))

    # ---- three-way feature screening ------------------------------------
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.feature_selection import RFE
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import LinearSVC
    from sklearn.metrics import roc_auc_score, roc_curve

    Xv = X.values
    scaler = StandardScaler().fit(Xv)
    Xs = scaler.transform(Xv)

    # 1) LASSO (L1 logistic), C chosen by inner CV
    best_c, best_auc = 1.0, -1
    for c in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]:
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=c,
                                 max_iter=5000, random_state=SEED)
        cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
        prob = cross_val_predict(clf, Xs, y, cv=cv, method="predict_proba")[:, 1]
        a = roc_auc_score(y, prob)
        if a > best_auc:
            best_c, best_auc = c, a
    lasso = LogisticRegression(penalty="l1", solver="liblinear", C=best_c,
                               max_iter=5000, random_state=SEED).fit(Xs, y)
    lasso_genes = [g for g, co in zip(X.columns, lasso.coef_[0]) if abs(co) > 1e-8]
    print(f"[LASSO] C={best_c} (inner-CV AUC={best_auc:.3f}), {len(lasso_genes)} genes")

    # 2) Random forest importance
    rf = RandomForestClassifier(n_estimators=1000, random_state=SEED,
                                class_weight="balanced", n_jobs=-1).fit(Xv, y)
    imp = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)
    rf_genes = list(imp.head(10).index)
    print(f"[RF] top10: {rf_genes}")

    # 3) SVM-RFE (linear SVM)
    rfe = RFE(LinearSVC(dual="auto", max_iter=20000, random_state=SEED),
              n_features_to_select=10, step=0.1).fit(Xs, y)
    svm_genes = list(X.columns[rfe.support_])
    print(f"[SVM-RFE] top10: {svm_genes}")

    sel = pd.DataFrame({
        "gene": X.columns,
        "in_LASSO": [g in lasso_genes for g in X.columns],
        "in_RF": [g in rf_genes for g in X.columns],
        "in_SVM_RFE": [g in svm_genes for g in X.columns],
        "rf_importance": imp.reindex(X.columns).values,
        "lasso_coef": lasso.coef_[0],
    }).merge(deg, on="gene", how="left")
    sel["n_methods"] = sel[["in_LASSO", "in_RF", "in_SVM_RFE"]].sum(axis=1)
    sel = sel.sort_values(["n_methods", "rf_importance"], ascending=False)
    sel.to_csv(os.path.join(TAB, "GSE123342_feature_selection.csv"), index=False,
               encoding="utf-8-sig")

    hub = sel.loc[sel["n_methods"] >= 2, "gene"].tolist()
    if not hub:
        hub = sel.head(8)["gene"].tolist()
    print(f"[hub] {len(hub)} genes selected by >=2 methods: {hub}")

    # ---- models + comparison with EGNS ----------------------------------
    # NOTE ON METHODOLOGY: selecting features on the full dataset and then
    # cross-validating those fixed features leaks label information from the
    # validation folds and inflates the AUC severely. We therefore report BOTH
    # (a) the leaky "selection outside CV" number, for transparency, and
    # (b) honest nested-CV numbers where the selector is a pipeline step and
    # refits inside every training fold.
    from sklearn.feature_selection import RFE, SelectFromModel
    from sklearn.pipeline import Pipeline

    scores = pd.read_csv(os.path.join(TAB, "GSE123342_module_scores.csv"),
                         index_col=0, encoding="utf-8-sig")
    results, roc_store = [], {}

    def evaluate(name: str, Xm: np.ndarray, selector=None):
        steps = [("sc", StandardScaler())]
        if selector is not None:
            steps.append(("sel", selector))
        steps.append(("clf", LogisticRegression(max_iter=5000, random_state=SEED)))
        pipe = Pipeline(steps)
        cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
        prob = cross_val_predict(pipe, Xm, y, cv=cv, method="predict_proba")[:, 1]
        auc = roc_auc_score(y, prob)
        boots = []
        for _ in range(1000):
            idx = rng.integers(0, len(y), len(y))
            if len(np.unique(y[idx])) < 2:
                continue
            boots.append(roc_auc_score(y[idx], prob[idx]))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        fpr, tpr, _ = roc_curve(y, prob)
        roc_store[name] = np.vstack([fpr, tpr])
        results.append({"model": name, "n_features": Xm.shape[1], "auc_cv": auc,
                        "auc_ci_low": lo, "auc_ci_high": hi})
        print(f"[model] {name:38s} AUC={auc:.3f} (95% CI {lo:.3f}-{hi:.3f})")
        return auc

    print("[models] leaky (selection performed on the full dataset):")
    leak_hub = evaluate("hub genes (selection OUTSIDE CV)", X[hub].values)
    leak_lasso = evaluate("LASSO genes (selection OUTSIDE CV)",
                          X[lasso_genes].values if lasso_genes else X[hub].values)

    print("[models] honest (nested CV; selector refits inside each fold):")
    sel_lasso = SelectFromModel(
        LogisticRegression(penalty="l1", solver="liblinear", C=0.3,
                           max_iter=5000, random_state=SEED),
        max_features=20, threshold=-np.inf)
    sel_rf = SelectFromModel(
        RandomForestClassifier(n_estimators=500, random_state=SEED,
                               class_weight="balanced", n_jobs=-1),
        max_features=20, threshold=-np.inf)
    sel_rfe = RFE(LinearSVC(dual="auto", max_iter=20000, random_state=SEED),
                  n_features_to_select=15, step=0.1)

    evaluate("EGNS (1 module score)", scores.loc[samples, ["EGNS"]].values)
    evaluate("all panel genes (no selection)", Xv)
    evaluate("nested: LASSO-select", Xv, sel_lasso)
    evaluate("nested: RF-select", Xv, sel_rf)
    evaluate("nested: SVM-RFE", Xv, sel_rfe)

    perf = pd.DataFrame(results)
    perf["leakage_note"] = ["leaky" if "OUTSIDE" in m else "ok" for m in perf["model"]]
    perf.to_csv(os.path.join(TAB, "GSE123342_ML_performance.csv"),
                index=False, encoding="utf-8-sig")
    print(f"[leakage] AUC inflation from selecting features outside CV: "
          f"{leak_hub:.3f} vs nested SVM-RFE {results[-1]['auc_cv']:.3f}")
    np.savez(os.path.join(TAB, "GSE123342_ML_roc.npz"),
             **{k.replace(" ", "_").replace("(", "").replace(")", ""): v
                for k, v in roc_store.items()})

    # ---- figure ---------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        for name, arr in roc_store.items():
            a = results[[r["model"] for r in results].index(name)]["auc_cv"]
            axes[0].plot(arr[0], arr[1], lw=1.6, label=f"{name} (AUC={a:.3f})")
        axes[0].plot([0, 1], [0, 1], "--", color="#94a3b8", lw=1)
        axes[0].set_xlabel("1 - specificity"); axes[0].set_ylabel("sensitivity")
        axes[0].set_title("Cross-validated ROC"); axes[0].legend(fontsize=7)
        top = sel.head(15).sort_values("rf_importance")
        axes[1].barh(top["gene"], top["rf_importance"], color="#475569")
        axes[1].set_title("Random-forest importance (top 15)")
        axes[1].tick_params(labelsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, "GSE123342_ML_ROC_importance.png"), dpi=200)
        plt.close()
        print(f"[fig] {os.path.join(FIG, 'GSE123342_ML_ROC_importance.png')}")
    except Exception as exc:  # noqa: BLE001
        print(f"[fig] skipped: {exc}")

    print("[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
