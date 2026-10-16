#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Nature-style figure set for the EGNS manuscript.

Style follows Nature/Nature Medicine conventions:
  * Arial, small type, single-column 89 mm / double-column 183 mm
  * no top/right spines, thin ticks, no chartjunk
  * bold lowercase panel letters
  * colour-blind-safe palette, effect sizes with 95% CI
Outputs PDF (vector) + PNG (600 dpi) per figure.

Usage:
    python make_figures_nature.py --outdir results/figures_nature
"""
from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------- style ---
MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM

C = {
    "blue":   "#2166AC",
    "red":    "#B2182B",
    "orange": "#E08214",
    "teal":   "#01665E",
    "purple": "#762A83",
    "green":  "#4D9221",
    "grey":   "#767676",
    "lgrey":  "#D9D9D9",
    "dark":   "#1A1A1A",
}
GROUP_COLORS = {
    "StableCAD": C["blue"], "AcuteMI": C["red"], "MI_30d": C["orange"], "MI_1y": C["teal"],
    "Sham": C["blue"], "MI_day1": C["red"], "MI_day7": C["orange"],
    "Non-diseased": C["blue"], "Dilated cardiomyopathy": C["red"],
}
LIN = {
    "Endothelial": "#2166AC", "LymphaticEC": "#92C5DE", "Pericyte": "#8C510A",
    "SMC": "#D8B365", "Fibroblast": "#4D9221", "Macrophage": "#B2182B",
    "Neutrophil": "#E08214", "Tcell": "#762A83", "Bcell": "#C51B7D",
    "Nkcells": "#01665E", "Adipocyte": "#F1B6DA", "Neuronal": "#999999",
}


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.2,
        "axes.labelsize": 6.5,
        "axes.titlesize": 6.8,
        "axes.linewidth": 0.55,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "legend.fontsize": 5.6,
        "legend.frameon": False,
        "legend.handlelength": 1.1,
        "legend.handletextpad": 0.4,
        "legend.labelspacing": 0.25,
        "lines.linewidth": 0.9,
        "lines.markersize": 2.4,
        "savefig.bbox": None,
        "savefig.pad_inches": 0.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left")


def save(fig, outdir, name, width_mm=None):
    """Write PDF (vector) + PNG at exactly the declared figure size.

    bbox_inches is deliberately NOT set to 'tight': Nature specifies physical
    column widths (89 / 183 mm), so the canvas must not be re-scaled by
    auto-cropping around long tick labels.

    FINAL_NAME renumbers the figures so that they appear in order of first
    citation in the manuscript (Nature requirement), e.g. the mouse
    single-cell figure is cited in section 3.2 and therefore becomes Fig. 2,
    while the human-blood figure is cited in section 3.8 and becomes Fig. 9.
    """
    name = FINAL_NAME.get(name, name)
    os.makedirs(outdir, exist_ok=True)
    w_in = fig.get_size_inches()[0]
    if width_mm is not None:
        w_mm = w_in * 25.4
        if abs(w_mm - width_mm) > 1.5:
            print(f"  [warn] {name}: {w_mm:.0f} mm != target {width_mm} mm")
    for ext in ("pdf", "png"):
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p, dpi=600 if ext == "png" else None, format=ext)
    plt.close(fig)
    print(f"[fig] {os.path.join(outdir, name)}.pdf/.png  ({w_in*25.4:.0f} mm)")


def T(rel):
    return os.path.join(ROOT, rel)


# internal figure name -> final name, ordered by first citation in the manuscript
FINAL_NAME = {
    "Fig1_data_overview":            "Fig1_data_overview",
    "Fig3_mouse_singlecell":         "Fig2_mouse_singlecell",
    "Fig9_immune_time_series":       "Fig3_immune_time_series",
    "Fig5_spatial_neighbourhood":    "Fig4_spatial_neighbourhood",
    "Fig4_cell_communication":       "Fig5_cell_communication",
    "Fig6_cross_species":            "Fig6_cross_species",
    "Fig7_candidate_decision_matrix": "Fig7_candidate_decision_matrix",
    "Fig8_feature_selection":        "Fig8_feature_selection",
    "Fig2_human_blood":              "Fig9_human_blood",
}


# =========================================================== Figure 2 =====
def fig_human_blood(outdir):
    scores = pd.read_csv(T("results/tables/GSE123342_module_scores.csv"),
                         index_col=0, encoding="utf-8-sig")
    groups = pd.read_csv(T("results/tables/GSE123342_sample_groups.csv"),
                         index_col=0, encoding="utf-8-sig")
    gcol = "group" if "group" in groups.columns else groups.columns[-1]
    g = groups[gcol].reindex(scores.index)
    keep = g.isin(["StableCAD", "AcuteMI", "MI_30d", "MI_1y"])
    scores, g = scores[keep], g[keep]
    order = ["StableCAD", "AcuteMI", "MI_30d", "MI_1y"]

    stats = pd.read_csv(T("results/tables/GSE123342_group_stats.csv"),
                        encoding="utf-8-sig").set_index("score")
    cand = pd.read_csv(T("results/tables/GSE123342_candidate_genes.csv"),
                       encoding="utf-8-sig")
    ml = pd.read_csv(T("results/tables/GSE123342_ML_performance.csv"),
                     encoding="utf-8-sig")

    fig = plt.figure(figsize=(DOUBLE, 4.55))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.92],
                          hspace=0.62, wspace=0.52,
                          left=0.085, right=0.975, top=0.93, bottom=0.10)

    mods = ["GCI", "GDS", "EIS", "EAS", "EIDS", "MFS", "FNS", "MRS", "EGNS"]

    # a  module-score heatmap of medians by group
    ax = fig.add_subplot(gs[0, 0])
    med = pd.DataFrame({k: scores.loc[g == k, mods].median() for k in order})
    z = med.sub(med.mean(axis=1), axis=0).div(med.std(axis=1).replace(0, np.nan), axis=0)
    im = ax.imshow(z.values, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["CAD", "AMI", "30d", "1y"], rotation=45, ha="right")
    ax.set_yticks(range(len(mods)))
    ax.set_yticklabels(mods)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("row z-score", fontsize=5.6)
    cb.ax.tick_params(labelsize=5.2, width=0.5)
    cb.outline.set_linewidth(0.4)
    ax.set_title("Module score (median)", pad=3)
    panel(ax, "a", dx=-0.30, dy=1.18)

    # b  EGNS distribution across groups
    ax = fig.add_subplot(gs[0, 1])
    data = [scores.loc[g == k, "EGNS"].values for k in order]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", lw=0.8),
                    whiskerprops=dict(lw=0.55), capprops=dict(lw=0.55))
    for patch, k in zip(bp["boxes"], order):
        patch.set_facecolor(GROUP_COLORS[k]); patch.set_edgecolor("none")
        patch.set_alpha(0.92)
    for i, (d, k) in enumerate(zip(data, order), start=1):
        ax.scatter(np.random.default_rng(1).normal(i, 0.07, len(d)), d,
                   s=1.4, color=C["dark"], alpha=0.35, linewidths=0)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(["CAD", "AMI", "30d", "1y"])
    ax.set_ylabel("EGNS")
    q = stats.loc["EGNS", "AcuteMI_vs_StableCAD_q"]
    ax.set_title(f"AMI vs CAD  q = {q:.3g}", pad=3)
    panel(ax, "b", dx=-0.30, dy=1.18)

    # c  follow-up trajectory of key scores
    ax = fig.add_subplot(gs[0, 2])
    traj = pd.DataFrame({k: scores.loc[g == k, mods].mean() for k in order})
    zt = traj.sub(med.mean(axis=1), axis=0).div(med.std(axis=1).replace(0, np.nan), axis=0)
    for m, col in zip(["GDS", "EAS", "EIS", "GCI"], [C["red"], C["orange"], C["blue"], C["teal"]]):
        ax.plot(range(4), zt.loc[m], marker="o", color=col, label=m, lw=1.0)
    ax.set_xticks(range(4)); ax.set_xticklabels(["CAD", "AMI", "30d", "1y"])
    ax.set_ylabel("z-score (vs group median)")
    ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    ax.legend(ncol=2, loc="upper right", columnspacing=0.7)
    ax.set_title("Recovery trajectory", pad=3)
    panel(ax, "c", dx=-0.30, dy=1.18)

    # d  ROC curves
    ax = fig.add_subplot(gs[0, 3])
    rocs = {}
    for f in os.listdir(T("results/tables")):
        if f.startswith("roc_") and f.endswith(".npy"):
            rocs[f[4:-4]] = np.load(T(f"results/tables/{f}"))
    perf = pd.read_csv(T("results/tables/GSE123342_model_performance.csv"), encoding="utf-8-sig")
    auc = dict(zip(perf["model"], perf["auc_cv"]))
    cmap = {"EGNS": C["red"], "GDS": C["orange"], "EGNS+GCI": C["teal"], "all_primary": C["purple"]}
    for k in ["EGNS", "GDS", "EGNS+GCI", "all_primary"]:
        if k in rocs:
            ax.plot(rocs[k][0], rocs[k][1], color=cmap.get(k, C["grey"]),
                    label=f"{k} ({auc.get(k, float('nan')):.2f})", lw=1.0)
    ax.plot([0, 1], [0, 1], ls=(0, (2, 2)), color=C["lgrey"], lw=0.6)
    ax.set_xlabel("1 - specificity"); ax.set_ylabel("Sensitivity")
    ax.legend(loc="lower right")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_title("Cross-validated ROC", pad=3)
    panel(ax, "d", dx=-0.30, dy=1.18)

    # e  candidate effect sizes, human blood
    ax = fig.add_subplot(gs[1, :2])
    sub = cand[cand["gene"].isin(["MMP9", "PLAUR", "NDST1", "EXTL3", "ADAM17", "HPSE",
                                  "EXT1", "HS6ST1", "HAS2", "HYAL2", "GPC1", "SDC1"])].copy()
    sub = sub.sort_values("delta")
    colors = [C["red"] if q < 0.05 else C["lgrey"] for q in sub["q"]]
    ax.barh(range(len(sub)), sub["delta"], color=colors, height=0.62)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([f"{g}   q={q:.2g}" for g, q in zip(sub["gene"], sub["q"])])
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_xlabel("Median difference, acute MI − stable CAD (log$_2$)")
    ax.set_title("Candidate molecules in peripheral blood", pad=3)
    panel(ax, "e", dx=-0.24, dy=1.10)

    # f  leakage correction
    ax = fig.add_subplot(gs[1, 2:])
    order_m = ["hub genes (selection OUTSIDE CV)", "LASSO genes (selection OUTSIDE CV)",
               "EGNS (1 module score)", "all panel genes (no selection)",
               "nested: LASSO-select", "nested: RF-select", "nested: SVM-RFE"]
    lab = {"hub genes (selection OUTSIDE CV)": "hub genes\n(selection outside CV)",
           "LASSO genes (selection OUTSIDE CV)": "LASSO genes\n(selection outside CV)",
           "EGNS (1 module score)": "EGNS (1 score)",
           "all panel genes (no selection)": "all panel genes",
           "nested: LASSO-select": "nested LASSO",
           "nested: RF-select": "nested random forest",
           "nested: SVM-RFE": "nested SVM-RFE"}
    ml = ml.set_index("model").reindex(order_m).dropna(subset=["auc_cv"]).reset_index()
    cols = [C["grey"] if "OUTSIDE" in m else C["blue"] for m in ml["model"]]
    y = np.arange(len(ml))
    ax.barh(y, ml["auc_cv"], color=cols, height=0.6)
    for i, r in ml.iterrows():
        ax.plot([r["auc_ci_low"], r["auc_ci_high"]], [i, i], color=C["dark"], lw=0.7)
        ax.plot([r["auc_cv"]], [i], "o", color=C["dark"], ms=1.9)
        ax.text(max(r["auc_ci_high"], r["auc_cv"]) + 0.012, i, f"{r['auc_cv']:.3f}",
                va="center", fontsize=5.2)
    ax.set_yticks(y); ax.set_yticklabels([lab.get(m, m) for m in ml["model"]])
    ax.axvline(0.5, color=C["lgrey"], lw=0.6, ls=(0, (2, 2)))
    ax.set_xlim(0.45, 1.06)
    ax.set_xlabel("Cross-validated AUC (95% CI)")
    ax.set_title("Selection-leakage correction", pad=3)
    panel(ax, "f", dx=-0.24, dy=1.10)

    save(fig, outdir, "Fig2_human_blood")


# =========================================================== Figure 3 =====
def fig_mouse_sc(outdir):
    pb = pd.read_csv(T("results/GSE176092/tables/endothelial_pseudobulk_scores.csv"),
                     encoding="utf-8-sig")
    spec = pd.read_csv(T("results/GSE176092/tables/candidate_cell_specificity.csv"),
                       encoding="utf-8-sig")
    reg = pd.read_csv(T("results/regulons/tables/regulon_vs_EGNS_correlation.csv"),
                      encoding="utf-8-sig")
    hlv = pd.read_csv(T("results/regulons/tables/regulon_high_vs_low_EGNS.csv"),
                      encoding="utf-8-sig")
    tc = pd.read_csv(T("results/cellcomm/tables/niche_program_timecourse.csv"),
                     encoding="utf-8-sig")

    fig = plt.figure(figsize=(DOUBLE, 4.7))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.95], hspace=0.60, wspace=0.55,
                          left=0.10, right=0.975, top=0.93, bottom=0.16)

    conds = ["Sham", "MI_day1", "MI_day7"]
    mods = ["GCI", "GDS", "EIS", "EAS", "EIDS", "MFS", "FNS", "MRS", "EGNS"]
    # the mouse table prefixes module columns with 'mod_'; harmonise the names
    pb = pb.rename(columns={f"mod_{m}": m for m in mods if f"mod_{m}" in pb.columns})

    # a  endothelial module scores by condition (heatmap of means)
    ax = fig.add_subplot(gs[0, 0])
    m = pb.groupby("condition")[mods].mean().reindex(conds)
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)
    im = ax.imshow(z.T.values, cmap="RdBu_r", vmin=-1.4, vmax=1.4, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(["Sham", "D1", "D7"], rotation=45, ha="right")
    ax.set_yticks(range(len(mods))); ax.set_yticklabels(mods)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=5.2, width=0.5); cb.outline.set_linewidth(0.4)
    cb.set_label("row z-score", fontsize=5.6)
    ax.set_title("Endothelial modules", pad=3)
    panel(ax, "a", dx=-0.34, dy=1.16)

    # b  EGNS trajectory with per-animal points
    ax = fig.add_subplot(gs[0, 1])
    for i, c in enumerate(conds):
        v = pb.loc[pb["condition"] == c, "EGNS"].values
        ax.scatter(np.full(len(v), i), v, s=14, color=GROUP_COLORS[c],
                   edgecolor="white", linewidth=0.4, zorder=3)
        ax.plot([i - 0.16, i + 0.16], [v.mean()] * 2, color=C["dark"], lw=1.1, zorder=4)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Sham", "D1", "D7"])
    ax.set_ylabel("EGNS (endothelium)")
    ax.set_title("Per-animal EGNS", pad=3)
    panel(ax, "b", dx=-0.30, dy=1.16)

    # c  regulon activity vs EGNS (lollipop)
    ax = fig.add_subplot(gs[0, 2])
    reg = reg.sort_values("rho_vs_EGNS")
    cols = [C["red"] if r > 0 else C["blue"] for r in reg["rho_vs_EGNS"]]
    ax.hlines(range(len(reg)), 0, reg["rho_vs_EGNS"], color=cols, lw=1.0)
    ax.plot(reg["rho_vs_EGNS"], range(len(reg)), "o", color=C["dark"], ms=2.2)
    ax.set_yticks(range(len(reg))); ax.set_yticklabels(reg["regulon"])
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_xlabel("Spearman $\\rho$ vs EGNS")
    ax.set_title("Regulon activity", pad=3)
    panel(ax, "c", dx=-0.46, dy=1.16)

    # d  high vs low EGNS delta
    ax = fig.add_subplot(gs[0, 3])
    hlv = hlv.sort_values("delta")
    cols = [C["red"] if d > 0 else C["blue"] for d in hlv["delta"]]
    ax.barh(range(len(hlv)), hlv["delta"], color=cols, height=0.6)
    ax.set_yticks(range(len(hlv))); ax.set_yticklabels(hlv["regulon"])
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_xlabel("$\\Delta$ high − low EGNS")
    ax.set_title("High vs low EGNS", pad=3)
    panel(ax, "d", dx=-0.46, dy=1.16)

    # e  candidate cell specificity dot plot (THE key panel)
    ax = fig.add_subplot(gs[1, :2])
    piv = spec.pivot(index="gene", columns="lineage", values="pct_expressing")
    order_g = ["Hpse", "Hpse2", "Mmp9", "Mmp2", "Hyal2", "Hyal1", "Plaur", "Mmp14",
               "Adam17", "Sulf1", "Sulf2", "Neu1"]
    order_l = ["Endothelial", "LymphaticEC", "Pericyte", "SMC", "Fibroblast",
               "Macrophage", "Neutrophil", "Tcell", "Bcell"]
    order_g = [g for g in order_g if g in piv.index]
    order_l = [l for l in order_l if l in piv.columns]
    P = piv.loc[order_g, order_l]
    xs, ys, ss, cc = [], [], [], []
    for i, g in enumerate(order_g):
        for j, l in enumerate(order_l):
            v = P.iloc[i, j]
            xs.append(j); ys.append(i); ss.append(3 + v * 1.15); cc.append(v)
    sc = ax.scatter(xs, ys, s=ss, c=cc, cmap="YlOrRd", vmin=0, vmax=90,
                    edgecolor=C["grey"], linewidth=0.25)
    ax.set_xticks(range(len(order_l)))
    ax.set_xticklabels(order_l, rotation=45, ha="right")
    ax.set_yticks(range(len(order_g))); ax.set_yticklabels(order_g)
    ax.set_ylim(len(order_g) - 0.5, -0.5)
    ax.set_xlim(-0.5, len(order_l) - 0.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.02)
    cb.set_label("% cells expressing", fontsize=5.6)
    cb.ax.tick_params(labelsize=5.2, width=0.5); cb.outline.set_linewidth(0.4)
    ax.set_title("Candidate genes: cell specificity in mouse heart", pad=3)
    panel(ax, "e", dx=-0.13, dy=1.10)

    # f  niche program time course
    ax = fig.add_subplot(gs[1, 2:])
    nice = {"Macrophage_profibrotic": ("Macrophage pro-fibrotic", C["red"]),
            "Fibroblast_activated": ("Fibroblast activated", C["green"]),
            "Endothelial_activated": ("Endothelial activated", C["blue"])}
    for prog, (lab, col) in nice.items():
        s = tc[tc["program"] == prog]
        xs, ys = [], []
        for c in conds:
            v = s.loc[s["condition"] == c, "program_score"]
            xs.append(conds.index(c)); ys.append(v.mean())
        ax.plot(xs, ys, marker="o", color=col, label=lab, lw=1.0)
        for c in conds:
            v = s.loc[s["condition"] == c, "program_score"]
            ax.scatter(np.full(len(v), conds.index(c)), v, s=8, color=col,
                       alpha=0.55, linewidths=0)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Sham", "D1", "D7"])
    ax.set_ylabel("Programme score")
    ax.legend(loc="upper right", ncol=1)
    ax.set_title("Niche programme time course (per animal)", pad=3)
    panel(ax, "f", dx=-0.11, dy=1.10)

    save(fig, outdir, "Fig3_mouse_singlecell")


# =========================================================== Figure 4 =====
def fig_cellcomm(outdir):
    ok = pd.read_csv(T("results/cellcomm/tables/cellcomm_significant.csv"),
                     encoding="utf-8-sig")
    resp = pd.read_csv(T("results/cellcomm/tables/receptor_downstream_response.csv"),
                       encoding="utf-8-sig")

    fig = plt.figure(figsize=(DOUBLE, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.25], wspace=0.85,
                          left=0.165, right=0.975, top=0.88, bottom=0.22)

    order_l = ["Endothelial", "LymphaticEC", "Pericyte", "SMC", "Fibroblast",
               "Macrophage", "Neutrophil", "Tcell", "Bcell"]
    order_p = ["FN1-ITGB1", "LGALS3-ITGB1", "THBS1-CD47", "SPP1-CD44",
               "SPP1-ITGAV_ITGB1", "VEGFA-KDR", "CCL2-CCR2", "CXCL12-CXCR4",
               "PDGFB-PDGFRB", "IL1B-IL1R1", "GAS6-AXL", "DLL4-NOTCH1"]
    ok = ok[ok["sender"].isin(order_l) & ok["receiver"].isin(order_l)]

    # a  sender->endothelial communication
    ax = fig.add_subplot(gs[0, 0])
    endo = ok[ok["receiver"] == "Endothelial"].copy()
    endo = endo.sort_values("score").tail(10)
    cols = [C["red"] if p <= 0.01 else C["grey"] for p in endo["p_perm"]]
    ax.barh(range(len(endo)), endo["score"], color=cols, height=0.6)
    ax.set_yticks(range(len(endo)))
    ax.set_yticklabels([f"{p}  ({s[:9]})" for p, s in zip(endo["LR_pair"], endo["sender"])])
    ax.set_xlabel("Communication score")
    ax.set_title("Axes targeting endothelium", pad=3)
    panel(ax, "a", dx=-0.62, dy=1.10)

    # b  endothelial as sender
    ax = fig.add_subplot(gs[0, 1])
    send = ok[ok["sender"].str.contains("Endothelial")].copy()
    send = send.sort_values("score").tail(8)
    cols = [C["blue"] if p <= 0.01 else C["grey"] for p in send["p_perm"]]
    ax.barh(range(len(send)), send["score"], color=cols, height=0.6)
    ax.set_yticks(range(len(send)))
    ax.set_yticklabels([f"{p}  ({s[:9]})" for p, s in zip(send["LR_pair"], send["receiver"])])
    ax.set_xlabel("Communication score")
    ax.set_title("Axes emitted by endothelium", pad=3)
    panel(ax, "b", dx=-0.62, dy=1.10)

    # c  receptor -> downstream response
    ax = fig.add_subplot(gs[0, 2])
    resp = resp.sort_values("rho_vs_program", ascending=False).head(14)
    cols = [LIN.get(l, C["grey"]) for l in resp["lineage"]]
    ax.barh(range(len(resp)), resp["rho_vs_program"], color=cols, height=0.62)
    ax.set_yticks(range(len(resp)))
    ax.set_yticklabels([f"{r} · {l[:4]}" for r, l in zip(resp["receptor"], resp["lineage"])])
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_xlabel("Spearman $\\rho$ (receptor vs downstream programme)")
    handles = [Patch(facecolor=LIN[k], label=k) for k in
               ["Macrophage", "Fibroblast", "Endothelial"] if k in set(resp["lineage"])]
    ax.legend(handles=handles, loc="lower right")
    ax.set_title("Downstream receptor response", pad=3)
    panel(ax, "c", dx=-0.52, dy=1.10)

    save(fig, outdir, "Fig4_cell_communication")


# =========================================================== Figure 5 =====
def fig_spatial(outdir):
    s = pd.read_csv(T("results/visium_spatial/tables/visium_spatial_summary.csv"),
                    encoding="utf-8-sig")
    s["tp"] = s["section"].str.extract(r"(Sham|MI_day\d+)")

    fig = plt.figure(figsize=(DOUBLE, 2.55))
    gs = fig.add_gridspec(1, 4, wspace=0.60, left=0.075, right=0.975,
                          top=0.87, bottom=0.26)

    order = ["Sham", "MI_day1", "MI_day7", "MI_day14"]

    # a  Moran's I per section
    ax = fig.add_subplot(gs[0, 0])
    for i, tp in enumerate(order):
        v = s.loc[s["tp"] == tp, "moranI_EGNS"].values
        if not len(v):
            continue
        ax.scatter(np.full(len(v), i), v, s=16, color=GROUP_COLORS.get(tp, C["grey"]),
                   edgecolor="white", linewidth=0.4, zorder=3)
        ax.plot([i - 0.18, i + 0.18], [v.mean()] * 2, color=C["dark"], lw=1.1, zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["Sham", "D1", "D7", "D14"])
    ax.set_ylabel("Moran's $I$ (EGNS)")
    ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    ax.set_title("Spatial autocorrelation", pad=3)
    panel(ax, "a", dx=-0.34, dy=1.10)

    # b  neighbour enrichment
    ax = fig.add_subplot(gs[0, 1])
    keys = [("d_nbr_FNS", "p_nbr_FNS", "FNS"), ("d_nbr_MFS", "p_nbr_MFS", "MFS"),
            ("d_nbr_GDS", "p_nbr_GDS", "GDS"), ("d_nbr_EIS", "p_nbr_EIS", "EIS")]
    for i, (dk, pk, lab) in enumerate(keys):
        v = s[dk].astype(float).values
        p = s[pk].astype(float).values
        col = C["red"] if np.nanmean(p) < 0.05 else C["grey"]
        ax.bar(i, np.nanmean(v), color=col, width=0.6,
               yerr=np.nanstd(v), error_kw=dict(lw=0.6, capsize=1.4, ecolor=C["dark"]))
    ax.axhline(0, color=C["dark"], lw=0.55)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([k[2] for k in keys])
    ax.set_ylabel("$\\Delta$ neighbour mean\n(high − low EGNS)")
    ax.set_title("Neighbourhood enrichment", pad=3)
    panel(ax, "b", dx=-0.34, dy=1.10)

    # c  same-spot vs neighbour correlation
    ax = fig.add_subplot(gs[0, 2])
    pairs = [("rho_same_FNS", "rho_nbr_FNS", "FNS"),
             ("rho_same_MFS", "rho_nbr_MFS", "MFS"),
             ("rho_same_GDS", "rho_nbr_GDS", "GDS"),
             ("rho_same_EIS", "rho_nbr_EIS", "EIS")]
    w = 0.36
    for i, (a, b, lab) in enumerate(pairs):
        ax.bar(i - w / 2, s[a].astype(float).mean(), width=w, color=C["blue"],
               yerr=s[a].astype(float).std(), error_kw=dict(lw=0.6, capsize=1.2, ecolor=C["dark"]))
        ax.bar(i + w / 2, s[b].astype(float).mean(), width=w, color=C["orange"],
               yerr=s[b].astype(float).std(), error_kw=dict(lw=0.6, capsize=1.2, ecolor=C["dark"]))
    ax.axhline(0, color=C["dark"], lw=0.55)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels([p[2] for p in pairs])
    ax.set_ylabel("Spearman $\\rho$ with EGNS")
    handles = [Patch(facecolor=C["blue"], label="same spot"),
               Patch(facecolor=C["orange"], label="neighbour mean")]
    ax.legend(handles=handles, loc="lower left")
    ax.set_title("Same-spot vs spatial coupling", pad=3)
    panel(ax, "c", dx=-0.34, dy=1.10)

    # d  Sham vs MI contrast
    ax = fig.add_subplot(gs[0, 3])
    mi = s[s["tp"] != "Sham"]
    sh = s[s["tp"] == "Sham"]
    rows = [("rho_nbr_FNS", "$\\rho$ FNS"), ("rho_nbr_MFS", "$\\rho$ MFS"),
            ("d_nbr_FNS", "$\\Delta$FNS"), ("moranI_EGNS", "Moran $I$")]
    for i, (k, lab) in enumerate(rows):
        ax.bar(i - w / 2, sh[k].astype(float).mean(), width=w, color=C["blue"])
        ax.bar(i + w / 2, mi[k].astype(float).mean(), width=w, color=C["red"])
    ax.axhline(0, color=C["dark"], lw=0.55)
    ax.set_xticks(range(len(rows)))
    # shorter tick labels and a gentler rotation: at 30 deg with a ~40 mm panel
    # the rotated labels ran into the neighbour panel and were clipped.
    ax.set_xticklabels([r[1] for r in rows], rotation=18, ha="right", fontsize=4.6)
    ax.set_ylabel("Mean statistic")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.22)
    handles = [Patch(facecolor=C["blue"], label="Sham"),
               Patch(facecolor=C["red"], label="MI (all)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.28),
              ncol=2, fontsize=4.4, handlelength=1.2, handletextpad=0.35,
              columnspacing=0.8, borderpad=0.1, frameon=False)
    ax.set_title("Sham vs infarcted", pad=3)
    panel(ax, "d", dx=-0.34, dy=1.10)

    save(fig, outdir, "Fig5_spatial_neighbourhood")


# =========================================================== Figure 6 =====
def fig_crossspecies(outdir):
    d = pd.read_csv(T("results/GSE183852/tables/endothelial_DCM_vs_control.csv"),
                    encoding="utf-8-sig")
    mouse = {"GCI": -0.102, "GDS": 0.018, "EIS": 1.097, "EAS": -0.223,
             "EIDS": 0.115, "MFS": 0.075, "FNS": -0.343, "MRS": 0.394}
    # mouse direction: EGNS/others up in MI. Derive sign from the paper's tables.
    mouse_dir = {"GCI": -1, "GDS": +1, "EIS": -1, "EAS": +1, "EIDS": +1,
                 "MFS": +1, "FNS": +1, "MRS": +1, "EGNS": +1}

    fig = plt.figure(figsize=(DOUBLE, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.68,
                          left=0.085, right=0.975, top=0.87, bottom=0.26)

    allr = d[d["technology"] == "all"].set_index("score")
    snr = d[d["technology"] == "Single Nuclei"].set_index("score")

    # a  human DCM effect sizes with CI
    ax = fig.add_subplot(gs[0, 0])
    order = ["EAS", "EGNS", "EIDS", "MRS", "GDS", "EIS", "GCI", "MFS", "FNS"]
    order = [o for o in order if o in allr.index]
    ax.barh(range(len(order)), [allr.loc[o, "delta"] for o in order],
            color=[C["red"] if allr.loc[o, "q_all"] < 0.05 else C["grey"] for o in order],
            height=0.6)
    for i, o in enumerate(order):
        ax.text(allr.loc[o, "delta"] + (0.02 if allr.loc[o, "delta"] >= 0 else -0.02), i,
                f"q={allr.loc[o,'q_all']:.2g}", va="center",
                ha="left" if allr.loc[o, "delta"] >= 0 else "right", fontsize=5.2)
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    ax.set_xlabel("$\\Delta$ (DCM − non-diseased)")
    ax.set_xlim(-0.45, 0.85)
    ax.set_title("Human DCM endothelium (all samples)", pad=3)
    panel(ax, "a", dx=-0.30, dy=1.10)

    # b  direction concordance mouse vs human
    ax = fig.add_subplot(gs[0, 1])
    comps = ["GCI", "GDS", "EIS", "EAS", "EIDS", "MRS", "EGNS"]
    for i, k in enumerate(comps):
        m = mouse_dir.get(k, 0)
        h = np.sign(allr.loc[k, "delta"]) if k in allr.index else 0
        agree = (m == h)
        ax.plot([0, 1], [m, h], color=C["green"] if agree else C["red"],
                lw=1.2, marker="o", ms=3.4, zorder=3)
        ax.text(-0.06, m, k, ha="right", va="center", fontsize=5.6)
        ax.text(1.06, h, k, ha="left", va="center", fontsize=5.6)
    ax.set_xlim(-0.5, 1.5); ax.set_ylim(-1.5, 1.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["mouse MI", "human DCM"])
    ax.set_yticks([-1, 0, 1]); ax.set_yticklabels(["down", "0", "up"])
    ax.set_ylabel("Direction of change")
    handles = [Patch(facecolor=C["green"], label="concordant"),
               Patch(facecolor=C["red"], label="discordant")]
    ax.legend(handles=handles, loc="lower center", ncol=1)
    ax.set_title("Directional concordance", pad=3)
    panel(ax, "b", dx=-0.16, dy=1.10)

    # c  technology-stratified comparison
    ax = fig.add_subplot(gs[0, 2])
    comps2 = ["GCI", "GDS", "EIS", "EAS", "EIDS", "MRS", "EGNS"]
    w = 0.36
    for i, k in enumerate(comps2):
        ax.bar(i - w / 2, allr.loc[k, "delta"], width=w, color=C["blue"])
        ax.bar(i + w / 2, snr.loc[k, "delta"], width=w, color=C["orange"])
    ax.axhline(0, color=C["dark"], lw=0.55)
    ax.set_xticks(range(len(comps2)))
    ax.set_xticklabels(comps2, rotation=45, ha="right")
    ax.set_ylabel("$\\Delta$ (DCM − control)")
    handles = [Patch(facecolor=C["blue"], label="all samples"),
               Patch(facecolor=C["orange"], label="single-nucleus only")]
    ax.legend(handles=handles, loc="lower left")
    ax.set_title("Technology-stratified", pad=3)
    panel(ax, "c", dx=-0.30, dy=1.10)

    save(fig, outdir, "Fig6_cross_species")


# =========================================================== Figure 1 =====
def fig_overview(outdir):
    """Data inventory + workflow summary derived from Table 1 and section 3.1."""
    datasets = [
        ("GSE123342", "human blood", 192, "discovery: clinical score"),
        ("GSE176092", "mouse heart sc/sn + Visium", 8473, "primary: endothelial EGNS"),
        ("GSE206787", "mouse heart sc + Visium", 49568, "mural-cell axis"),
        ("GSE183852", "human heart sc/sn", 92770, "cross-species check"),
    ]
    fig = plt.figure(figsize=(DOUBLE, 2.15))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.15], wspace=0.68,
                          left=0.115, right=0.975, top=0.85, bottom=0.22)

    # a  datasset sizes
    ax = fig.add_subplot(gs[0, 0])
    names = [d[0] for d in datasets]
    vals = [d[2] for d in datasets]
    cols = [C["blue"], C["red"], C["teal"], C["purple"]]
    ax.barh(range(len(datasets)), vals, color=cols, height=0.62)
    for i, v in enumerate(vals):
        ax.text(v * 1.05, i, f"{v:,}", va="center", fontsize=5.4)
    ax.set_yticks(range(len(datasets))); ax.set_yticklabels(names)
    ax.set_xscale("log")
    ax.set_xlabel("Cells / samples analysed")
    ax.set_xlim(50, 4e5)
    ax.set_title("Datasets", pad=3)
    panel(ax, "a", dx=-0.34, dy=1.13)

    # b  module coverage
    ax = fig.add_subplot(gs[0, 1])
    cov = {"GCI": 30, "GDS": 19, "EIS": 17, "EAS": 15, "EIDS": 16,
           "MFS": 14, "FNS": 18, "MRS": 19}
    tot = {"GCI": 30, "GDS": 20, "EIS": 17, "EAS": 15, "EIDS": 16,
           "MFS": 15, "FNS": 18, "MRS": 21}
    ks = list(cov.keys())
    ax.bar(range(len(ks)), [cov[k] / tot[k] for k in ks],
           color=[C["blue"]] * len(ks), width=0.62)
    ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks, rotation=45, ha="right")
    ax.set_ylim(0, 1.12); ax.set_ylabel("Genes detected / defined")
    ax.set_title("Module coverage", pad=3)
    panel(ax, "b", dx=-0.30, dy=1.13)

    # c  headline results summary
    ax = fig.add_subplot(gs[0, 2])
    rows = [
        ("Human blood EGNS\nAMI vs CAD", 0.688, "AUC"),
        ("Mouse endothelial EGNS\nD1 vs Sham", 0.891, "$\\Delta$ z"),
        ("Spatial coupling\nEGNS vs FNS", 0.510, "$\\rho$ (nbr)"),
        ("HPSE endothelial\nspecificity", 0.013, "frac."),
    ]
    ax.barh(range(len(rows)), [r[1] for r in rows],
            color=[C["red"], C["red"], C["teal"], C["grey"]], height=0.6)
    for i, r in enumerate(rows):
        ax.text(r[1] + 0.02, i, f"{r[1]:.3f}", va="center", fontsize=5.4)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Effect size")
    ax.set_title("Headline findings", pad=3)
    panel(ax, "c", dx=-0.40, dy=1.13)

    save(fig, outdir, "Fig1_data_overview")


def fig_candidate_matrix(outdir):
    """Table 3 (candidate-axis lock-in criteria) rendered as an evidence matrix,
    using ONLY quantities actually measured in this study; criteria that were
    not assessed are shown as such rather than being guessed."""
    spec = pd.read_csv(T("results/GSE176092/tables/candidate_cell_specificity.csv"),
                       encoding="utf-8-sig")
    rank = pd.read_csv(T("results/GSE176092/tables/candidate_endothelial_ranking.csv"),
                       encoding="utf-8-sig")
    blood = pd.read_csv(T("results/tables/GSE123342_candidate_genes.csv"),
                        encoding="utf-8-sig")

    cands = ["Hpse", "Hyal2", "Mmp9", "Plaur", "Adam17", "Mmp14", "Ext1", "Gpc1", "Sdc1", "Has2"]
    blood_map = {"Hpse": "HPSE", "Hyal2": "HYAL2", "Mmp9": "MMP9", "Plaur": "PLAUR",
                 "Adam17": "ADAM17", "Mmp14": "MMP14", "Ext1": "EXT1", "Gpc1": "GPC1",
                 "Sdc1": "SDC1", "Has2": "HAS2"}

    fig = plt.figure(figsize=(DOUBLE, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.55,
                          left=0.085, right=0.965, top=0.86, bottom=0.24)

    # --- a  measured evidence heatmap -------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ec = spec[spec["lineage"].str.contains("Endothelial")].set_index("gene")["pct_expressing"]
    mx = spec[spec["lineage"].str.contains("Endothelial")].set_index("gene")["lineage"]
    top_other = {}
    for g in spec["gene"].unique():
        s = spec[spec["gene"] == g]
        o = s[~s["lineage"].str.contains("Endothelial")]
        top_other[g] = float(o["pct_expressing"].max()) if len(o) else np.nan

    rows = []
    for c in cands:
        b = blood_map[c]
        bq = blood.loc[blood["gene"] == b, "q"]
        rows.append({
            "candidate": c,
            f"Endothelial\n% positive": ec.get(c, np.nan),
            "Highest other\nlineage %": top_other.get(c, np.nan),
            "Human blood\n$\\Delta$ (AMI$-$CAD)":
                float(blood.loc[blood["gene"] == b, "delta"].iloc[0]) if len(bq) else np.nan,
            "Human blood\n$-\\log_{10}q$":
                float(-np.log10(max(bq.iloc[0], 1e-12))) if len(bq) else np.nan,
        })
    M = pd.DataFrame(rows).set_index("candidate")
    norm = M.copy()
    for c in norm.columns:
        v = norm[c]
        norm[c] = (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else 0.5
    im = ax.imshow(norm.values, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels(M.columns, fontsize=5.2)
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels(M.index)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for i in range(len(M)):
        for j, c in enumerate(M.columns):
            v = M.iloc[i, j]
            txt = "NA" if not np.isfinite(v) else (f"{v:.2f}" if abs(v) < 10 else f"{v:.0f}")
            ax.text(j, i, txt, ha="center", va="center", fontsize=4.9,
                    color="white" if norm.iloc[i, j] > 0.55 else C["dark"])
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
    cb.set_label("column-normalised", fontsize=5.4)
    cb.ax.tick_params(labelsize=5.0, width=0.5); cb.outline.set_linewidth(0.4)
    ax.set_title("Measured evidence per candidate", pad=3)
    panel(ax, "a", dx=-0.20, dy=1.14)

    # --- b  Table 3 scorecard --------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    criteria = ["Cell specificity", "Temporal precedence", "Spatial consistency",
                "Cross-dataset replication", "Pathway interpretability", "Druggability"]
    #  2 = criterion met, 1 = partial, 0 = failed, -1 = not assessed in this study
    S = {
        "Hpse":   [0, -1, -1, 0, 2, 2],
        "Mmp9":   [0, -1, -1, 1, 2, 2],
        "Plaur":  [0, -1, -1, 1, 2, 1],
        "Adam17": [0, -1, -1, 1, 2, 2],
        "Mmp14":  [0, -1, -1, 1, 2, 2],
        "Hyal2":  [0, -1, -1, 0, 2, 1],
        "Ext1":   [0, -1, -1, 0, 2, 1],
        "Gpc1":   [0, -1, -1, 0, 2, 1],
        "Sdc1":   [0, -1, -1, 0, 2, 1],
        "Has2":   [0, -1, -1, 0, 2, 1],
    }
    names = list(S.keys())
    arr = np.array([S[n] for n in names], dtype=float)
    cmap = mpl.colors.ListedColormap(["#D9D9D9", C["red"], "#F6E8B1", C["green"]])
    bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
    im = ax.imshow(arr, cmap=cmap, norm=mpl.colors.BoundaryNorm(bounds, cmap.N),
                   aspect="auto")
    ax.set_xticks(range(len(criteria)))
    ax.set_xticklabels(criteria, rotation=40, ha="right", fontsize=5.2)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Table 3 criteria scorecard", pad=3)
    handles = [Patch(facecolor="#D9D9D9", label="not assessed"),
               Patch(facecolor=C["red"], label="failed"),
               Patch(facecolor="#F6E8B1", label="partial"),
               Patch(facecolor=C["green"], label="met")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.42),
              ncol=2, fontsize=5.0)
    panel(ax, "b", dx=-0.20, dy=1.14)

    save(fig, outdir, "Fig7_candidate_decision_matrix")


def fig_feature_selection(outdir):
    """Section 3.9: three-way feature screening of the glycocalyx panel."""
    sel = pd.read_csv(T("results/tables/GSE123342_feature_selection.csv"),
                      encoding="utf-8-sig")
    sel = sel.sort_values(["n_methods", "rf_importance"], ascending=False)
    top = sel.head(20).copy()

    fig = plt.figure(figsize=(DOUBLE, 2.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.55,
                          left=0.085, right=0.975, top=0.87, bottom=0.28)

    # a  random-forest importance, coloured by how many methods agreed
    ax = fig.add_subplot(gs[0, 0])
    cmap = {3: C["red"], 2: C["orange"], 1: C["lgrey"], 0: C["lgrey"]}
    y = np.arange(len(top))[::-1]
    ax.barh(y, top["rf_importance"],
            color=[cmap.get(int(n), C["lgrey"]) for n in top["n_methods"]],
            height=0.66)
    ax.set_yticks(y); ax.set_yticklabels(top["gene"], fontsize=5.4)
    ax.set_xlabel("Random-forest importance")
    handles = [Patch(facecolor=cmap[k], label=f"{k} method(s)") for k in (3, 2, 1)]
    ax.legend(handles=handles, loc="lower right", fontsize=5.0)
    ax.set_title("Top 20 panel genes", pad=3)
    panel(ax, "a", dx=-0.34, dy=1.12)

    # b  method-membership matrix
    ax = fig.add_subplot(gs[0, 1])
    genes = list(top["gene"])
    methods = [("in_LASSO", "LASSO"), ("in_RF", "Random\nforest"), ("in_SVM_RFE", "SVM-RFE")]
    Z = np.array([[1.0 if bool(top.iloc[i][m]) else 0.0 for m, _ in methods]
                  for i in range(len(top))])
    im = ax.imshow(Z.T, cmap="Greys", vmin=-0.35, vmax=1.0, aspect="auto")
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels([lab for _, lab in methods], fontsize=5.4)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=5.0)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for gi in range(Z.shape[0]):
        for mi in range(Z.shape[1]):
            if Z[gi, mi] > 0:
                # draw a marker rather than a glyph: Arial has no check mark
                ax.plot(gi, mi, "o", color="white", ms=1.7, zorder=3)
    ax.set_title(f"Selected by $\\geq$2 methods: "
                 f"{int((top['n_methods'] >= 2).sum())} of top 20", pad=3)
    panel(ax, "b", dx=-0.16, dy=1.12)

    save(fig, outdir, "Fig8_feature_selection")


def fig_immune_series(outdir):
    """GSE135310 cardiac immune time series (one animal per time point)."""
    frac = pd.read_csv(T("results/GSE135310/tables/lineage_fraction_by_timepoint.csv"),
                       index_col=0, encoding="utf-8-sig")
    ccr2 = pd.read_csv(T("results/GSE135310/tables/ccr2_recruitment.csv"),
                       encoding="utf-8-sig")
    traj = pd.read_csv(T("results/GSE135310/tables/macrophage_state_trajectory.csv"),
                       encoding="utf-8-sig")
    order = ["Steady state", "Day 1", "Day 3", "Day 5", "Day 7"]
    xlabels = ["St.st.", "D1", "D3", "D5", "D7"]
    frac = frac.reindex(order)
    ccr2 = ccr2.set_index("timepoint").reindex(order)
    traj = traj.set_index("timepoint").reindex(order)

    fig = plt.figure(figsize=(DOUBLE, 2.9))
    gs = fig.add_gridspec(1, 4, wspace=0.62, left=0.07, right=0.98,
                          top=0.86, bottom=0.26)

    # a  lineage composition
    ax = fig.add_subplot(gs[0, 0])
    show = ["Neutrophil", "Monocyte", "Macrophage", "Bcell"]
    other = 1.0 - frac[show].sum(axis=1)
    x = np.arange(len(order))
    ax.stackplot(x, [frac[c].values * 100 for c in show] + [other.values * 100],
                 labels=show + ["other"],
                 colors=[LIN["Neutrophil"], LIN["Macrophage"], LIN["Pericyte"],
                         LIN["Tcell"], "#E8E8E8"], alpha=0.95, lw=0.3)
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    ax.set_ylabel("% of CD45$^+$ cells")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper left", fontsize=4.9, ncol=1)
    ax.set_title("Immune composition", pad=3)
    panel(ax, "a", dx=-0.32, dy=1.12)

    # b  CCR2+ recruitment
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(x, ccr2["pct_Ccr2_pos_myeloid"], marker="o", color=C["red"],
            label="myeloid cells")
    ax.plot(x, ccr2["pct_Ccr2_pos_all"], marker="s", color=C["blue"],
            label="all cells")
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    ax.set_ylabel("Ccr2$^+$ (%)")
    ax.legend(loc="lower right", fontsize=5.0)
    ax.set_title("CCR2$^+$ recruitment", pad=3)
    panel(ax, "b", dx=-0.32, dy=1.12)

    # c  macrophage state transition
    ax = fig.add_subplot(gs[0, 2])
    states = [("Mac_profibrotic", "pro-fibrotic", C["red"]),
              ("Mac_inflammatory", "inflammatory", C["orange"]),
              ("Mono_classical", "classical mono", C["blue"]),
              ("Mac_antigenpresent", "antigen-presenting", C["teal"])]
    for col, lab, c in states:
        ax.plot(x, traj[col], marker="o", color=c, label=lab)
    ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    ax.set_ylabel("Macrophage state score")
    ax.legend(loc="upper left", fontsize=5.0)
    ax.set_title("Macrophage states", pad=3)
    panel(ax, "c", dx=-0.32, dy=1.12)

    # d  MFS vs integrated timeline
    ax = fig.add_subplot(gs[0, 3])
    ax.plot(x, traj["MFS"], marker="o", color=C["red"], label="MFS")
    ax.plot(x, traj["Mac_profibrotic"], marker="^", color=C["orange"],
            label="pro-fibrotic")
    # endothelial EGNS from GSE176092, aligned on the shared D1/D7 time points
    ex = [1, 4]
    ey = [0.891, 0.333]
    ax.plot(ex, ey, marker="D", ls=(0, (3, 2)), color=C["blue"],
            label="EGNS (GSE176092)")
    ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(xlabels)
    ax.set_ylabel("Score / $\\Delta$ vs Sham")
    # short labels + legend placed BELOW the axes: a 3-entry legend inside a
    # ~40 mm panel previously overlapped the curves and the adjacent panel.
    ax.set_ylim(top=ax.get_ylim()[1] * 1.18)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=1,
              fontsize=4.4, handlelength=1.4, handletextpad=0.4,
              labelspacing=0.25, borderpad=0.1, frameon=False)
    ax.set_title("Cascade alignment", pad=3)
    panel(ax, "d", dx=-0.32, dy=1.12)

    save(fig, outdir, "Fig9_immune_time_series")


def fig_virtual_ko(outdir):
    """In silico TF knockout: predicted effect on EGNS (endothelium) and on the
    pro-fibrotic macrophage programme."""
    d = pd.read_csv(T("results/virtual_ko/tables/virtual_ko_results.csv"),
                    encoding="utf-8-sig")
    ec = d[d["cell_type"] == "Endothelial"].sort_values("delta_readout")
    mp = d[d["cell_type"] == "Macrophage"].sort_values("delta_readout")

    fig = plt.figure(figsize=(DOUBLE, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.0], wspace=0.75,
                          left=0.085, right=0.975, top=0.86, bottom=0.24)

    # a  endothelial dEGNS
    ax = fig.add_subplot(gs[0, 0])
    cols = [C["red"] if p < 0.05 else (C["orange"] if p < 0.10 else C["lgrey"])
            for p in ec["p_empirical"]]
    y = np.arange(len(ec))
    ax.barh(y, ec["delta_readout"], color=cols, height=0.72)
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_yticks(y); ax.set_yticklabels(ec["ko_tf"], fontsize=5.0)
    ax.set_xlabel("$\\Delta$EGNS (in silico KO)")
    handles = [Patch(facecolor=C["red"], label="P < 0.05"),
               Patch(facecolor=C["orange"], label="P < 0.10"),
               Patch(facecolor=C["lgrey"], label="ns")]
    ax.legend(handles=handles, loc="lower right", fontsize=5.0)
    ax.set_title("Endothelial EGNS", pad=3)
    panel(ax, "a", dx=-0.42, dy=1.10)

    # b  module decomposition for the six most impactful TFs
    ax = fig.add_subplot(gs[0, 1])
    sel = pd.concat([ec.head(3), ec.tail(3)])
    parts = [("d_GDS", "GDS", C["red"]), ("d_EAS", "EAS", C["orange"]),
             ("d_EIDS", "EIDS", C["purple"]), ("d_GCI", "GCI", C["teal"]),
             ("d_EIS", "EIS", C["blue"])]
    w = 0.16
    x = np.arange(len(sel))
    for i, (col, lab, c) in enumerate(parts):
        ax.bar(x + (i - 2) * w, sel[col], width=w, color=c, label=lab)
    ax.axhline(0, color=C["dark"], lw=0.55)
    ax.set_xticks(x); ax.set_xticklabels(sel["ko_tf"], rotation=45, ha="right",
                                         fontsize=5.0)
    ax.set_ylabel("$\\Delta$ module score")
    ax.legend(ncol=3, loc="lower left", fontsize=5.0)
    ax.set_title("Module decomposition (6 TFs)", pad=3)
    panel(ax, "b", dx=-0.30, dy=1.10)

    # c  macrophage dMFS
    ax = fig.add_subplot(gs[0, 2])
    cols = [C["red"] if p < 0.05 else (C["orange"] if p < 0.10 else C["lgrey"])
            for p in mp["p_empirical"]]
    y = np.arange(len(mp))
    ax.barh(y, mp["delta_readout"], color=cols, height=0.72)
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_yticks(y); ax.set_yticklabels(mp["ko_tf"], fontsize=5.0)
    ax.set_xlabel("$\\Delta$MFS (in silico KO)")
    ax.set_title("Macrophage fibrosis score", pad=3)
    panel(ax, "c", dx=-0.42, dy=1.10)

    save(fig, outdir, "Fig10_virtual_knockout")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "figures_nature"))
    args = ap.parse_args()
    # resolve a relative --outdir against the PROJECT root, not the CWD: passing
    # a bare relative path previously wrote every figure into a stray
    # <cwd>/results/figures_nature while the manuscript pointed at ROOT.
    if not os.path.isabs(args.outdir):
        args.outdir = os.path.join(ROOT, args.outdir)
    os.makedirs(args.outdir, exist_ok=True)
    set_style()
    fig_overview(args.outdir)
    fig_human_blood(args.outdir)
    fig_mouse_sc(args.outdir)
    fig_cellcomm(args.outdir)
    fig_spatial(args.outdir)
    fig_crossspecies(args.outdir)
    fig_candidate_matrix(args.outdir)
    fig_feature_selection(args.outdir)
    fig_immune_series(args.outdir)
    fig_virtual_ko(args.outdir)
    print("[done] all figures written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
