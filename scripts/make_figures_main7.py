#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Consolidated MAIN figure set (7 figures) for the EGNS manuscript.

Reviewer-driven restructure: the main text should carry only the most credible
evidence, not the algorithm-iteration history. This script builds the three
COMPOSITE figures that merge formerly separate ones; the other four main
figures are reused unchanged from the existing generators.

Final numbering
---------------
  Fig.1  study design: datasets + module coverage + headline effects
         + descriptive vs Core-EGNS (per cell, per animal, TF audit)   [composite]
  Fig.2  endothelial temporal state                         (Fig2_mouse_singlecell)
  Fig.3  regulatory network                                 (Fig11_decoupler_regulons)
  Fig.4  immune temporal cascade                            (Fig3_immune_time_series)
  Fig.5  spatial neighbourhood                              (Fig4_spatial_neighbourhood)
  Fig.6  communication + immune-sheddase effector model                  [composite]
  Fig.7  cross-species + virtual perturbation                            [composite]
Supplementary (not main text):
  FigS1 heuristic 26-pair communication, FigS2 feature screening,
  FigS3 human-blood clinical score.

Usage:
    python make_figures_main7.py --outdir results/figures_nature
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from make_figures_nature import (  # noqa: E402
    C, LIN, GROUP_COLORS, DOUBLE, MM, ROOT, T, set_style, panel, save,
)


def _rp(p):
    return p if os.path.isabs(p) else os.path.join(ROOT, p)


# =========================================================== Fig 1 ==========
def fig1_framework(outdir):
    """Datasets + module coverage + headline effects + EGNS/Core-EGNS framework."""
    per = _rp("results/lto_egns/tables/egns_variants_per_cell.csv")
    pb = _rp("results/lto_egns/tables/egns_variants_pseudobulk.csv")
    audit = _rp("results/lto_egns/tables/tf_module_membership_audit.csv")
    if not all(os.path.exists(p) for p in (per, pb, audit)):
        print("[skip] Fig1 -- run analysis_lto_egns.py first")
        return
    d = pd.read_csv(per)
    p = pd.read_csv(pb)
    a = pd.read_csv(audit)
    a = a[a["actually_in_modules"].astype(str) != "-"]
    cond_order = ["Sham", "MI_day1", "MI_day7"]
    cond_colors = {"Sham": C["blue"], "MI_day1": C["red"], "MI_day7": C["orange"]}

    fig = plt.figure(figsize=(DOUBLE, 4.6))
    gs = fig.add_gridspec(2, 4, wspace=0.72, hspace=0.85,
                          left=0.075, right=0.978, top=0.90, bottom=0.10)

    # ---- row 1: study design ------------------------------------------
    ax = fig.add_subplot(gs[0, 0:2])
    datasets = [("GSE123342", 192), ("GSE176092", 8473),
                ("GSE206787", 49568), ("GSE183852", 92770)]
    names = [x[0] for x in datasets]
    vals = [x[1] for x in datasets]
    cols = [C["blue"], C["red"], C["teal"], C["purple"]]
    ax.barh(range(len(datasets)), vals, color=cols, height=0.62)
    for i, v in enumerate(vals):
        ax.text(v * 1.08, i, f"{v:,}", va="center", fontsize=5.4)
    ax.set_yticks(range(len(datasets)))
    ax.set_yticklabels(names)
    ax.set_xscale("log")
    ax.set_xlim(50, 6e5)
    ax.set_xlabel("Cells / samples analysed")
    ax.set_title("Datasets", pad=3)
    panel(ax, "a", dx=-0.30, dy=1.10)

    ax = fig.add_subplot(gs[0, 2])
    cov = {"GCI": 30, "GDS": 19, "EIS": 17, "EAS": 14, "EIDS": 16,
           "MFS": 15, "FNS": 18, "MRS": 19}
    tot = {"GCI": 30, "GDS": 19, "EIS": 17, "EAS": 15, "EIDS": 16,
           "MFS": 15, "FNS": 18, "MRS": 19}
    ks = list(cov.keys())
    ax.bar(range(len(ks)), [cov[k] / tot[k] for k in ks], color=C["blue"], width=0.62)
    ax.set_xticks(range(len(ks)))
    ax.set_xticklabels(ks, rotation=45, ha="right", fontsize=5.2)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Detected / defined")
    ax.set_title("Module coverage", pad=3)
    panel(ax, "b", dx=-0.32, dy=1.10)

    ax = fig.add_subplot(gs[0, 3])
    rows = [("Blood EGNS\nAUC", 0.688), ("Endo EGNS\nD1 vs Sham", 0.891),
            ("Spatial coupling\nEGNS-FNS", 0.510), ("Hpse endothelial\nspecificity", 0.013)]
    ax.barh(range(len(rows)), [r[1] for r in rows],
            color=[C["red"], C["red"], C["teal"], C["grey"]], height=0.6)
    for i, r in enumerate(rows):
        ax.text(r[1] + 0.025, i, f"{r[1]:.3f}", va="center", fontsize=5.2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.2)
    ax.set_xlim(0, 1.12)
    ax.set_xlabel("Effect size")
    ax.set_title("Headline effects", pad=3)
    panel(ax, "c", dx=-0.42, dy=1.10)

    # ---- row 2: EGNS / Core-EGNS framework ----------------------------
    for k, (col, title) in enumerate([
            ("descriptive_EGNS", "Descriptive EGNS (contains TFs)"),
            ("core_EGNS", "Core-EGNS = $z$(GDS) \u2212 $z$(GCI)  (no TF)")]):
        ax = fig.add_subplot(gs[1, k])
        data = [d.loc[d["condition"] == c, col].values for c in cond_order]
        bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                        medianprops=dict(color=C["dark"], lw=0.8),
                        whiskerprops=dict(lw=0.55), capprops=dict(lw=0.55))
        for patch, c in zip(bp["boxes"], cond_order):
            patch.set_facecolor(cond_colors[c])
            patch.set_alpha(0.75)
            patch.set_edgecolor(C["dark"])
            patch.set_linewidth(0.5)
        for i, c in enumerate(cond_order):
            v = d.loc[d["condition"] == c, col].values
            x = np.random.default_rng(0).normal(i + 1, 0.055, len(v))
            ax.plot(x, v, ".", ms=0.8, color=C["dark"], alpha=0.16, zorder=0)
        ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
        ax.set_xticklabels([c.replace("MI_", "MI ") for c in cond_order],
                           rotation=20, ha="right", fontsize=5.2)
        ax.set_ylabel("score (z)")
        ax.set_title(title, pad=3, fontsize=5.8)
        panel(ax, "de"[k], dx=-0.30, dy=1.10)

    ax = fig.add_subplot(gs[1, 2])
    for j, (col, mk, lab) in enumerate([("descriptive_EGNS", "o", "descriptive"),
                                        ("core_EGNS", "s", "Core-EGNS")]):
        xs, vs = [], []
        for i, c in enumerate(cond_order):
            for _, r in p[p["condition"] == c].iterrows():
                xs.append(i + (j - 0.5) * 0.26)
                vs.append(r[col])
        ax.plot(xs, vs, mk, ms=3.4, mfc="none",
                mec=C["dark"] if j == 0 else C["teal"], mew=0.8, label=lab, ls="none")
    ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    ax.set_xticks(range(len(cond_order)))
    ax.set_xticklabels([c.replace("MI_", "MI ") for c in cond_order],
                       rotation=20, ha="right", fontsize=5.2)
    ax.set_ylabel("pseudobulk mean (z)")
    ax.legend(loc="upper left", fontsize=4.8, handletextpad=0.3)
    ax.set_title("Per animal (n = 6)", pad=3)
    panel(ax, "f", dx=-0.32, dy=1.10)

    ax = fig.add_subplot(gs[1, 3])
    cnt = a["actually_in_modules"].value_counts()
    bars = ax.bar(range(len(cnt)), cnt.values,
                  color=[C["purple"] if m == "EIS" else
                         (C["red"] if m == "EAS" else C["orange"]) for m in cnt.index],
                  edgecolor=C["dark"], lw=0.5, width=0.6)
    for b, v in zip(bars, cnt.values):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.1, str(v), ha="center", fontsize=5.4)
    ax.set_xticks(range(len(cnt)))
    ax.set_xticklabels(list(cnt.index), fontsize=5.4)
    ax.set_ylim(0, cnt.max() + 1.2)
    ax.set_ylabel("TFs inside a score")
    ax.set_title("TF circularity audit", pad=3)
    panel(ax, "g", dx=-0.36, dy=1.10)

    save(fig, outdir, "Fig1_study_framework")


# =========================================================== Fig 6 ==========
def fig6_immune_sheddase(outdir):
    """Real LIANA communication + effector ranking + candidate decision."""
    ep = _rp("results/liana_communication/tables/liana_endothelial_receiver.csv")
    ag = _rp("results/liana_communication/tables/method_agreement.csv")
    eff = _rp("results/sctenifold_sheddase_ko/tables/sheddase_ko_effects.csv")
    dis = _rp("results/sctenifold_sheddase_ko/tables/sheddase_ko_network_distance.csv")
    spec_f = _rp("results/GSE176092/tables/candidate_cell_specificity.csv")
    blood_f = _rp("results/tables/GSE123342_candidate_genes.csv")
    if not all(os.path.exists(x) for x in (ep, eff, dis)):
        print("[skip] Fig6 -- run analysis_liana_communication.py and "
              "analysis_sctenifold_sheddase.py first")
        return

    fig = plt.figure(figsize=(DOUBLE, 4.5))
    gs = fig.add_gridspec(2, 3, wspace=0.70, hspace=0.80,
                          left=0.10, right=0.975, top=0.90, bottom=0.11)

    # a  endothelial-receiver interactions
    ec = pd.read_csv(ep).sort_values("aggregate_rank").head(14).copy()
    ec["pair"] = (ec["ligand_complex"].astype(str) + " \u2192 "
                  + ec["receptor_complex"].astype(str))
    ax = fig.add_subplot(gs[0, 0:2])
    y = np.arange(len(ec))[::-1]
    # x-axis = LIANA consensus rank as a strength score, larger = stronger.
    # NOTE aggregate_pval is EXACTLY 0 for every top-ranked pair (whole column
    # min = 0), so the previous -log10(P) axis collapsed onto a single vertical
    # line at 3.0 and carried no information at all. aggregate_rank spans three
    # decades over these 14 pairs (3.7e-06 .. 2.0e-03), so it drives the axis.
    rk = pd.to_numeric(ec["aggregate_rank"], errors="coerce").fillna(1.0).clip(lower=1e-12)
    xval = -np.log10(rk)
    # The source lineage is encoded by marker COLOUR plus a legend, NOT by a
    # second text layer. An earlier version drew the source name at
    # transform=ax.get_yaxis_transform() with x=-0.015, which lands exactly on
    # top of the y tick labels (14 pairs overlapped at ~83%).
    SRC_COL = {"Neutrophil": "#E08214", "Fibroblast": "#2166AC",
               "Macrophage": "#B2182B", "Pericyte": "#01665E",
               "LymphaticEC": "#762A83", "Endothelial": "#4D9221",
               "SMC": "#8C510A", "Tcell": "#999999", "Bcell": "#666666",
               "Adipocyte": "#CCCCCC"}
    srcs = list(dict.fromkeys(ec["source"].astype(str)))
    cols = [SRC_COL.get(s, "#767676") for s in ec["source"].astype(str)]
    ax.scatter(xval, y, s=26, c=cols, edgecolor=C["dark"], lw=0.4, zorder=3)
    ax.set_xlim(xval.min() - 0.40, xval.max() + 0.60)
    ax.set_yticks(y)
    ax.set_yticklabels(ec["pair"], fontsize=4.6)
    ax.set_xlabel("$-\\log_{10}$ (LIANA consensus rank)  —  larger = stronger",
                  fontsize=5.2)
    handles = [Line2D([], [], marker="o", ls="none", ms=3.6,
                      mfc=SRC_COL.get(s, "#767676"), mec=C["dark"], mew=0.4,
                      label=s) for s in srcs]
    ax.legend(handles=handles, loc="lower right", fontsize=4.4, frameon=False,
              handlelength=0.9, handletextpad=0.25, labelspacing=0.22,
              borderpad=0.1)
    ax.set_title("Endothelial-receiver interactions (LIANA)", pad=3)
    panel(ax, "a", dx=-0.46, dy=1.10)

    ax = fig.add_subplot(gs[0, 2])
    if os.path.exists(ag):
        m = pd.read_csv(ag).sort_values("spearman_vs_aggregate")
        yy = np.arange(len(m))
        ax.barh(yy, m["spearman_vs_aggregate"], color=C["teal"],
                edgecolor=C["dark"], lw=0.5, height=0.6)
        ax.set_yticks(yy)
        ax.set_yticklabels(m["method"], fontsize=4.8)
        ax.set_xlim(0, 1)
        ax.axvline(0.5, color=C["lgrey"], lw=0.6, ls="--")
        ax.set_xlabel("$\\rho$ vs consensus")
        ax.set_title("Method agreement", pad=3)
    panel(ax, "b", dx=-0.44, dy=1.10)

    # c  programme-specific sheddase effect
    e = pd.read_csv(eff)
    piv = e.pivot_table(index="ko", columns="readout_module", values="effect_ratio")
    mods = [m for m in ("MFS", "FNS", "GDS", "EAS", "GCI") if m in piv.columns]
    piv = piv[mods]
    if "MFS" in piv.columns:
        piv = piv.loc[piv["MFS"].sort_values(ascending=False).index]
    ax = fig.add_subplot(gs[1, 0])
    v = piv.values
    im = ax.imshow(v, cmap="PuOr_r", aspect="auto",
                   vmin=np.nanmin(v), vmax=np.nanmax(v))
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels(piv.columns, fontsize=5.2)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels(piv.index, fontsize=4.8)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            if np.isfinite(v[i, j]):
                ax.text(j, i, f"{v[i, j]:.1f}", ha="center", va="center", fontsize=4.0)
    cb = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.02)
    cb.set_label("effect ratio", fontsize=5.0)
    cb.ax.tick_params(labelsize=4.6)
    ax.set_title("Sheddase KO by programme", pad=3)
    panel(ax, "c", dx=-0.36, dy=1.10)

    # d  network-level rewiring
    ax = fig.add_subplot(gs[1, 1])
    dd = pd.read_csv(dis).sort_values("relative_frobenius_change")
    yy = np.arange(len(dd))
    cols = [C["red"] if k in ("Mmp9", "Hpse") else C["orange"] for k in dd["ko"]]
    ax.barh(yy, dd["relative_frobenius_change"], color=cols,
            edgecolor=C["dark"], lw=0.5, height=0.6)
    ax.set_yticks(yy)
    ax.set_yticklabels(dd["ko"], fontsize=4.8)
    ax.set_xlabel("relative Frobenius change")
    ax.set_title("Network rewiring", pad=3)
    panel(ax, "d", dx=-0.44, dy=1.10)

    # e  candidate evidence + f  scorecard
    ax = fig.add_subplot(gs[1, 2])
    if os.path.exists(spec_f) and os.path.exists(blood_f):
        spec = pd.read_csv(spec_f)
        blood = pd.read_csv(blood_f)
        cands = ["Hpse", "Hyal2", "Mmp9", "Plaur", "Adam17", "Mmp14"]
        bmap = {"Hpse": "HPSE", "Hyal2": "HYAL2", "Mmp9": "MMP9", "Plaur": "PLAUR",
                "Adam17": "ADAM17", "Mmp14": "MMP14"}
        ecp = spec[spec["lineage"].str.contains("Endothelial")].set_index("gene")["pct_expressing"]
        rows = []
        for c in cands:
            b = bmap[c]
            sub = blood.loc[blood["gene"] == b]
            rows.append({
                "candidate": c,
                "Endo %": float(ecp.get(c, np.nan)),
                "Blood $\\Delta$": float(sub["delta"].iloc[0]) if len(sub) else np.nan,
                "Blood $-\\log_{10}q$":
                    float(-np.log10(max(sub["q"].iloc[0], 1e-12))) if len(sub) else np.nan,
            })
        M = pd.DataFrame(rows).set_index("candidate")
        norm = M.copy()
        for c in norm.columns:
            vv = norm[c]
            norm[c] = (vv - vv.min()) / (vv.max() - vv.min()) if vv.max() > vv.min() else 0.5
        im = ax.imshow(norm.values, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(M.columns)))
        ax.set_xticklabels(M.columns, fontsize=4.8)
        ax.set_yticks(range(len(M)))
        ax.set_yticklabels(M.index, fontsize=4.8)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        for i in range(len(M)):
            for j in range(len(M.columns)):
                t = "NA" if not np.isfinite(M.iloc[i, j]) else f"{M.iloc[i, j]:.2f}"
                ax.text(j, i, t, ha="center", va="center", fontsize=4.0,
                        color="white" if norm.iloc[i, j] > 0.55 else C["dark"])
        ax.set_title("Effector candidates", pad=3)
    panel(ax, "e", dx=-0.36, dy=1.10)

    save(fig, outdir, "Fig6_immune_sheddase_model")


# =========================================================== Fig 7 ==========
def _model_panel(ax):
    """Final mechanistic model: four layers with their supporting evidence."""
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    boxes = [
        (0.015, "Regulatory layer\n(endothelium)\n\nKLF2 / ERG identity $\\downarrow$\n"
                "NF-$\\kappa$B / HIF1A / AP-1 $\\uparrow$", "#FDECEA"),
        (0.265, "Effector layer\n(immune-derived)\n\nMMP9 / PLAUR / ADAM17\n"
                "macrophage, neutrophil", "#FEF3E2"),
        (0.515, "Endothelial state\n\nGlycocalyx\nnon-resolution\n"
                "(Core-EGNS $\\uparrow$)", "#E8F0F8"),
        (0.765, "Tissue phenotype\n\nFibrotic niche\n(MFS / FNS $\\uparrow$)", "#EDF3EA"),
    ]
    w = 0.22
    for x, label, fc in boxes:
        ax.add_patch(plt.Rectangle((x, 0.40), w, 0.44, transform=ax.transAxes,
                                   facecolor=fc, edgecolor=C["dark"], lw=0.6,
                                   clip_on=False))
        ax.text(x + w / 2, 0.62, label, ha="center", va="center", fontsize=4.9,
                transform=ax.transAxes)
    for x in (0.235, 0.485, 0.735):
        ax.annotate("", xy=(x + 0.030, 0.62), xytext=(x, 0.62),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C["grey"]))
    evid = ["CollecTRI / ULM activity (\u00a73.3.1)\n"
            "Core-EGNS, TF-free (\u00a73.10.1)",
            "expression specificity (\u00a73.7)\n"
            "network edge-ablation KO (\u00a73.10.2)",
            "neighbourhood coupling (\u00a73.5)\n"
            "composition-adjusted (\u00a73.5.2)"]
    for x, t in zip((0.235, 0.485, 0.735), evid):
        ax.text(x + 0.015, 0.355, t, ha="center", va="top", fontsize=4.2,
                color=C["grey"], transform=ax.transAxes, linespacing=1.5)
    ax.text(0.5, 0.955, "Working model: immune-derived sheddase attack on an "
            "identity-compromised endothelium", ha="center", va="center",
            fontsize=5.6, fontweight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.115, "HPSE, the prespecified candidate, is refuted as the "
            "endothelial primary axis and retained only as a refuted candidate",
            ha="center", va="center", fontsize=4.3, style="italic",
            color=C["grey"], transform=ax.transAxes)


def fig7_crossspecies_perturbation(outdir):
    """Cross-species concordance + the final mechanistic model.

    NOTE: the CellOracle-STYLE vector-field perturbation is deliberately NOT
    plotted -- it failed its own trajectory-validity check and was withdrawn
    (manuscript Section 3.10.3). Only the validated perturbation readouts
    (delta EGNS, network edge-ablation) are cited in the model panel.
    """
    cs = _rp("results/GSE183852/tables/endothelial_DCM_vs_control.csv")
    if not os.path.exists(cs):
        print("[skip] Fig7 -- missing cross-species results")
        return

    d = pd.read_csv(cs, encoding="utf-8-sig")
    allr = d[d["technology"] == "all"].set_index("score")
    snr = d[d["technology"] == "Single Nuclei"].set_index("score")
    mouse_dir = {"GCI": -1, "GDS": +1, "EIS": -1, "EAS": +1, "EIDS": +1,
                 "MFS": +1, "FNS": +1, "MRS": +1, "EGNS": +1}

    fig = plt.figure(figsize=(DOUBLE, 4.9))
    gs = fig.add_gridspec(2, 3, wspace=0.72, hspace=0.80,
                          left=0.085, right=0.968, top=0.885, bottom=0.135)

    # ---------------- row 1: cross-species (a-c) ------------------------
    ax = fig.add_subplot(gs[0, 0])
    order = [o for o in ["EAS", "EGNS", "EIDS", "MRS", "GDS", "EIS", "GCI", "MFS", "FNS"]
             if o in allr.index]
    ax.barh(range(len(order)), [allr.loc[o, "delta"] for o in order],
            color=[C["red"] if allr.loc[o, "q_all"] < 0.05 else C["grey"] for o in order],
            height=0.6)
    ax.axvline(0, color=C["dark"], lw=0.55)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=5.0)
    ax.set_xlabel("$\\Delta$ (DCM \u2212 non-diseased)", fontsize=5.2)
    ax.set_xlim(-0.45, 0.9)
    ax.set_title("Human DCM endothelium", pad=3, fontsize=5.8)
    panel(ax, "a", dx=-0.36, dy=1.12)

    ax = fig.add_subplot(gs[0, 1])
    comps = ["GCI", "GDS", "EIS", "EAS", "EIDS", "MRS", "EGNS"]
    Mdir = np.zeros((len(comps), 2), dtype=float)
    for i, k in enumerate(comps):
        Mdir[i, 0] = mouse_dir.get(k, 0)
        Mdir[i, 1] = np.sign(allr.loc[k, "delta"]) if k in allr.index else 0
    im = ax.imshow(Mdir, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["mouse MI", "human DCM"], fontsize=5.0)
    ax.set_yticks(range(len(comps)))
    ax.set_yticklabels(comps, fontsize=5.2)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for i in range(len(comps)):
        agree = Mdir[i, 0] == Mdir[i, 1]
        ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 2.0, 1.0, fill=False,
                                   edgecolor=C["green"] if agree else C["red"],
                                   lw=1.0, zorder=3))
    handles = [Patch(facecolor="none", edgecolor=C["green"], label="concordant"),
               Patch(facecolor="none", edgecolor=C["red"], label="discordant")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=2, fontsize=4.6, frameon=False, handlelength=1.2,
              handletextpad=0.35, columnspacing=0.8, borderpad=0.1)
    ax.set_title("Directional concordance", pad=3, fontsize=5.8)
    panel(ax, "b", dx=-0.40, dy=1.12)

    ax = fig.add_subplot(gs[0, 2])
    w = 0.36
    for i, k in enumerate(comps):
        ax.bar(i - w / 2, allr.loc[k, "delta"], width=w, color=C["blue"])
        ax.bar(i + w / 2, snr.loc[k, "delta"], width=w, color=C["orange"])
    ax.axhline(0, color=C["dark"], lw=0.55)
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels(comps, rotation=45, ha="right", fontsize=5.0)
    ax.set_ylabel("$\\Delta$ (DCM \u2212 control)", fontsize=5.2)
    handles = [Patch(facecolor=C["blue"], label="all"),
               Patch(facecolor=C["orange"], label="snRNA only")]
    ax.legend(handles=handles, loc="lower left", fontsize=4.8, frameon=False)
    ax.set_title("Technology-stratified", pad=3, fontsize=5.8)
    panel(ax, "c", dx=-0.36, dy=1.12)

    # ---------------- row 2: final mechanistic model --------------------
    ax = fig.add_subplot(gs[1, 0:3])
    _model_panel(ax)
    panel(ax, "d", dx=-0.030, dy=1.10)

    save(fig, outdir, "Fig7_crossspecies_perturbation")
    print("  [note] vector-field perturbation panels omitted: that analysis "
          "was withdrawn after failing its trajectory-validity check")

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "figures_nature"))
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(ROOT, args.outdir)
    os.makedirs(outdir, exist_ok=True)
    set_style()
    todo = {"1": fig1_framework, "6": fig6_immune_sheddase,
            "7": fig7_crossspecies_perturbation}
    for k, fn in todo.items():
        if args.only and k != args.only:
            continue
        print(f"[Fig{k}]")
        try:
            fn(outdir)
        except Exception as exc:  # noqa: BLE001
            print(f"  [error] Fig{k}: {type(exc).__name__}: {exc}")
            import traceback
            traceback.print_exc()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
