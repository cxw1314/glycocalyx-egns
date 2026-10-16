#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Nature-style figure set, v2: figures introduced by the methodological revision.

  Fig11  EGNS circularity resolution -- descriptive vs Core-EGNS (no TF),
         per cell and per animal, plus the TF-membership audit
  Fig12  real LIANA consensus cell-cell communication (replaces the heuristic)
  Fig13  real decoupler / CollecTRI regulon activity + concordance with the
         curated-target approximation used in v1.1
  Fig14  scTenifoldKnk-equivalent sheddase knockout in macrophages

Every panel reads the CSV emitted by the corresponding analysis script; a
figure whose input is absent is skipped with a notice rather than fabricated.
Outputs PDF (vector) + PNG (600 dpi) at Nature double-column width 183 mm.
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MM = 1 / 25.4
DOUBLE = 183 * MM

C = {"blue": "#2166AC", "red": "#B2182B", "orange": "#E08214",
     "teal": "#01665E", "purple": "#762A83", "green": "#4D9221",
     "grey": "#767676", "lgrey": "#D9D9D9", "dark": "#1A1A1A"}
COND_COLORS = {"Sham": C["blue"], "MI_day1": C["red"], "MI_day7": C["orange"]}
COND_ORDER = ["Sham", "MI_day1", "MI_day7"]
MOD_ORDER = ["MFS", "FNS", "GDS", "EAS", "GCI"]


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.2, "axes.labelsize": 6.5, "axes.titlesize": 6.8,
        "axes.linewidth": 0.55, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8, "legend.fontsize": 5.6,
        "legend.frameon": False, "xtick.major.width": 0.55,
        "ytick.major.width": 0.55, "xtick.major.size": 1.8,
        "ytick.major.size": 1.8, "figure.dpi": 150, "savefig.dpi": 600,
    })


def panel(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left")


def save(fig, outdir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches=None)
    plt.close(fig)
    print(f"  [saved] {name}.pdf / {name}.png")


def rp(p):
    return p if os.path.isabs(p) else os.path.join(ROOT, p)


def have(*paths):
    return all(os.path.exists(rp(p)) for p in paths)


# ------------------------------------------------------------------ Fig11 ---
def fig11(outdir):
    per = "results/lto_egns/tables/egns_variants_per_cell.csv"
    pb = "results/lto_egns/tables/egns_variants_pseudobulk.csv"
    audit = "results/lto_egns/tables/tf_module_membership_audit.csv"
    if not have(per, pb, audit):
        print("[skip] Fig11 -- run analysis_lto_egns.py first")
        return
    d = pd.read_csv(rp(per))
    p = pd.read_csv(rp(pb))
    a = pd.read_csv(rp(audit))
    a = a[a["actually_in_modules"].astype(str) != "-"]

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 118 * MM),
                             gridspec_kw={"width_ratios": [1, 1], "hspace": 0.55,
                                          "wspace": 0.30})

    # (a, b) per-cell distributions, descriptive vs core
    for k, (col, title) in enumerate([
            ("descriptive_EGNS", "Descriptive EGNS  (contains TFs)"),
            ("core_EGNS", "Core-EGNS  $=z$(GDS)$-z$(GCI)  (no TF)")]):
        ax = axes[0, k]
        data = [d.loc[d["condition"] == c, col].values for c in COND_ORDER]
        bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                        medianprops=dict(color=C["dark"], lw=0.8),
                        whiskerprops=dict(lw=0.55), capprops=dict(lw=0.55))
        for patch, c in zip(bp["boxes"], COND_ORDER):
            patch.set_facecolor(COND_COLORS[c])
            patch.set_alpha(0.75)
            patch.set_edgecolor(C["dark"])
            patch.set_linewidth(0.5)
        for i, c in enumerate(COND_ORDER):
            v = d.loc[d["condition"] == c, col].values
            x = np.random.default_rng(0).normal(i + 1, 0.055, len(v))
            ax.plot(x, v, ".", ms=0.8, color=C["dark"], alpha=0.16, zorder=0)
        ax.set_xticklabels([c.replace("MI_", "MI ") for c in COND_ORDER],
                           rotation=20, ha="right")
        ax.set_ylabel("score (z)")
        ax.set_title(title, pad=3)
        panel(ax, "ab"[k])
    axes[0, 0].axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    axes[0, 1].axhline(0, color=C["lgrey"], lw=0.5, zorder=0)

    # (c) per-animal pseudobulk
    ax = axes[1, 0]
    xpos = np.arange(len(COND_ORDER))
    for j, (col, mk, lab) in enumerate([
            ("descriptive_EGNS", "o", "descriptive (n = 6 animals)"),
            ("core_EGNS", "s", "Core-EGNS (n = 6 animals)")]):
        vals, xs = [], []
        for i, c in enumerate(COND_ORDER):
            sub = p[p["condition"] == c]
            for _, r in sub.iterrows():
                vals.append(r[col])
                xs.append(i + (j - 0.5) * 0.26)
        ax.plot(xs, vals, mk, ms=3.4, mfc="none", mec=C["dark"] if j == 0 else C["teal"],
                mew=0.8, label=lab, ls="none")
    ax.set_xticks(xpos)
    ax.set_xticklabels([c.replace("MI_", "MI ") for c in COND_ORDER],
                       rotation=20, ha="right")
    ax.set_ylabel("pseudobulk mean (z)")
    ax.axhline(0, color=C["lgrey"], lw=0.5, zorder=0)
    ax.legend(loc="upper left", handletextpad=0.3)
    ax.set_title("per-animal direction is preserved", pad=3)
    panel(ax, "c")
    ax.annotate("both definitions:\nSham < D1, D7 intermediate",
                xy=(0.98, 0.03), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=5.2, color=C["grey"])

    # (d) TF circularity audit
    ax = axes[1, 1]
    cnt = a["actually_in_modules"].value_counts()
    bars = ax.bar(range(len(cnt)), cnt.values,
                  color=[C["red"] if m != "EIS" else C["purple"] for m in cnt.index],
                  edgecolor=C["dark"], lw=0.5, width=0.6)
    for b, v in zip(bars, cnt.values):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.12, str(v),
                ha="center", fontsize=5.6)
    ax.set_xticks(range(len(cnt)))
    ax.set_xticklabels([f"{m}\n({cnt[m]} TF{'s' if cnt[m] > 1 else ''})"
                        for m in cnt.index], fontsize=5.4)
    ax.set_ylabel("TFs embedded in a scored module")
    ax.set_ylim(0, cnt.max() + 1.4)
    ax.set_title(f"{int(a['circular_for_descriptive_EGNS'].sum())} of 30 panel TFs "
                 "are inside the descriptive score", pad=3)
    panel(ax, "d")

    fig.suptitle("EGNS circularity: the score is not an artefact of TF membership",
                 x=0.02, ha="left", fontsize=7.4, fontweight="bold", y=1.005)
    save(fig, outdir, "Fig13_EGNS_circularity")


# ------------------------------------------------------------------ Fig12 ---
def fig12(outdir):
    ep = "results/liana_communication/tables/liana_endothelial_receiver.csv"
    ag = "results/liana_communication/tables/method_agreement.csv"
    if not have(ep):
        print("[skip] Fig12 -- run analysis_liana_communication.py first")
        return
    ec = pd.read_csv(rp(ep))
    ec = ec.sort_values("aggregate_rank").head(18)
    ec["pair"] = ec["ligand_complex"].astype(str) + " -> " + ec["receptor_complex"].astype(str)

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 62 * MM),
                             gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.42})
    ax = axes[0]
    y = np.arange(len(ec))[::-1]
    # CellPhoneDB permutation p-values can be exactly 0 with finite permutations;
    # floor at the resolution of the test (1/1000) so -log10 stays finite.
    pv = pd.to_numeric(ec["aggregate_pval"], errors="coerce").fillna(1.0).clip(lower=1e-3)
    sc = ax.scatter(-np.log10(pv), y,
                    s=26, c=ec["aggregate_rank"], cmap="viridis_r",
                    edgecolor=C["dark"], lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(ec["pair"], fontsize=5.0)
    for yi, src in zip(y, ec["source"]):
        ax.text(-0.02, yi, str(src)[:12], transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=4.8, color=C["grey"])
    ax.set_xlabel("$-\\log_{10}$ (LIANA consensus $P$, floored at $10^{-3}$)")
    ax.set_title("Endothelial-receiver interactions", pad=3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("consensus rank", fontsize=5.4)
    cb.ax.tick_params(labelsize=5)
    panel(ax, "a", dx=-0.42)

    ax = axes[1]
    if have(ag):
        m = pd.read_csv(rp(ag)).sort_values("spearman_vs_aggregate")
        yy = np.arange(len(m))
        ax.barh(yy, m["spearman_vs_aggregate"], color=C["teal"],
                edgecolor=C["dark"], lw=0.5, height=0.62)
        ax.set_yticks(yy)
        ax.set_yticklabels(m["method"], fontsize=5.2)
        ax.set_xlabel("Spearman $\\rho$ vs consensus")
        ax.set_xlim(0, 1)
        ax.axvline(0.5, color=C["lgrey"], lw=0.6, ls="--")
        ax.set_title("Method agreement", pad=3)
    panel(ax, "b", dx=-0.30)
    fig.suptitle("Real LIANA consensus communication (mouse OmniPath resource)",
                 x=0.02, ha="left", fontsize=7.4, fontweight="bold", y=1.03)
    save(fig, outdir, "Fig12_LIANA_communication")


# ------------------------------------------------------------------ Fig13 ---
def fig13(outdir):
    pb = "results/decoupler_regulons/tables/regulon_activity_pseudobulk.csv"
    cc = "results/decoupler_regulons/tables/collecTRI_vs_curated_concordance.csv"
    if not have(pb):
        print("[skip] Fig13 -- run analysis_decoupler_regulons.py first")
        return
    act = pd.read_csv(rp(pb), index_col=0)
    act = act.loc[[c for c in COND_ORDER if c in act.index]]
    order = act.abs().mean(axis=0).sort_values(ascending=False).index
    act = act[order]

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 68 * MM),
                             gridspec_kw={"width_ratios": [2.0, 1], "wspace": 0.34})
    ax = axes[0]
    v = act.values
    lim = np.nanmax(np.abs(v))
    im = ax.imshow(v, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(act.columns)))
    ax.set_xticklabels(act.columns, rotation=90, fontsize=5.0)
    ax.set_yticks(range(len(act.index)))
    ax.set_yticklabels([c.replace("MI_", "MI ") for c in act.index], fontsize=5.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.015)
    cb.set_label("regulon score (ULM)", fontsize=5.4)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("CollecTRI regulon activity, endothelial cells", pad=3)
    panel(ax, "a", dx=-0.16)

    ax = axes[1]
    if have(cc):
        ccdf = pd.read_csv(rp(cc)).dropna(subset=["pearson_r_curated_vs_collecTRI"])
        ax.scatter(ccdf["curated_D1_minus_Sham"], ccdf["pearson_r_curated_vs_collecTRI"],
                   s=20, c=C["purple"], edgecolor=C["dark"], lw=0.4, alpha=0.85)
        # 19 TF labels previously overlapped heavily; repel them with adjustText
        # and fall back to unlabelled points if it is unavailable.
        texts = []
        try:
            from adjustText import adjust_text
            for _, r in ccdf.iterrows():
                texts.append(ax.text(r["curated_D1_minus_Sham"],
                                     r["pearson_r_curated_vs_collecTRI"],
                                     r["tf"], fontsize=4.0, color=C["dark"]))
            adjust_text(texts, ax=ax, expand=(1.06, 1.15),
                        arrowprops=dict(arrowstyle="-", lw=0.28,
                                        color=C["grey"], shrinkA=1, shrinkB=1))
        except Exception:
            for t in texts:
                t.remove()
            print("  [warn] adjustText unavailable; TF labels omitted (panels stay legible)")
        ax.axhline(0.5, color=C["lgrey"], lw=0.6, ls="--")
        ax.set_xlabel("curated proxy: $\\Delta$(D1 $-$ Sham)")
        ax.set_ylabel("Pearson $r$ vs CollecTRI")
        ax.margins(x=0.16, y=0.18)
        ax.set_title("Approximation vs real regulons", pad=3)
    panel(ax, "b", dx=-0.30)
    fig.suptitle("Resource-backed regulon inference (decoupler / CollecTRI)",
                 x=0.02, ha="left", fontsize=7.4, fontweight="bold", y=1.03)
    save(fig, outdir, "Fig11_decoupler_regulons")


# ------------------------------------------------------------------ Fig14 ---
def fig14(outdir):
    eff = "results/sctenifold_sheddase_ko/tables/sheddase_ko_effects.csv"
    dis = "results/sctenifold_sheddase_ko/tables/sheddase_ko_network_distance.csv"
    if not have(eff):
        print("[skip] Fig14 -- run analysis_sctenifold_sheddase.py first")
        return
    e = pd.read_csv(rp(eff))
    piv = e.pivot_table(index="ko", columns="readout_module", values="effect_ratio")
    piv = piv[[m for m in MOD_ORDER if m in piv.columns]]
    if "MFS" in piv.columns:
        piv = piv.loc[piv["MFS"].sort_values(ascending=False).index]

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 62 * MM),
                             gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.34})
    ax = axes[0]
    v = piv.values
    im = ax.imshow(v, cmap="PuOr_r", vmin=np.nanmin(v), vmax=np.nanmax(v),
                   aspect="auto")
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels(piv.columns, fontsize=5.6)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels(piv.index, fontsize=5.6)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            if np.isfinite(v[i, j]):
                ax.text(j, i, f"{v[i, j]:.2f}", ha="center", va="center",
                        fontsize=4.6, color=C["dark"])
    cb = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.015)
    cb.set_label("effect ratio vs genome-wide", fontsize=5.4)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("Programme-specific rewiring", pad=3)
    panel(ax, "a", dx=-0.22)

    ax = axes[1]
    if have(dis):
        dd = pd.read_csv(rp(dis)).sort_values("relative_frobenius_change")
        yy = np.arange(len(dd))
        ax.barh(yy, dd["relative_frobenius_change"], color=C["orange"],
                edgecolor=C["dark"], lw=0.5, height=0.62)
        ax.set_yticks(yy)
        ax.set_yticklabels(dd["ko"], fontsize=5.4)
        ax.set_xlabel("relative Frobenius change of scGRN")
        ax.set_title("Network-level rewiring", pad=3)
    panel(ax, "b", dx=-0.30)
    fig.suptitle("scTenifoldKnk-equivalent sheddase knockout (macrophage scGRN)",
                 x=0.02, ha="left", fontsize=7.4, fontweight="bold", y=1.03)
    save(fig, outdir, "Fig14_sheddase_KO")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results/figures_nature")
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    outdir = rp(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    set_style()
    todo = {"11": fig11, "12": fig12, "13": fig13, "14": fig14}
    for k, fn in todo.items():
        if args.only and k != args.only:
            continue
        print(f"[Fig{k}]")
        try:
            fn(outdir)
        except Exception as e:
            print(f"  [error] Fig{k}: {type(e).__name__}: {e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
