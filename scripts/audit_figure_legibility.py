#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Programmatic figure-legibility audit.

Renders every main/supplementary figure and measures, for every visible Text
artist (tick labels, annotations, legend entries, titles, axis labels), its
true on-canvas bounding box. Reports:
  * pairwise OVERLAPS (intersection area > threshold of the smaller box)
  * texts that extend OUTSIDE the parent axes (a common cause of clipping and
    of "text runs into the neighbouring panel")

This gives an objective, reproducible check instead of eyeballing PNGs.

Usage:
    python audit_figure_legibility.py
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import make_figures_nature as N  # noqa: E402
import make_figures_main7 as M7  # noqa: E402

ROOT = N.ROOT
OUT = os.path.join(ROOT, "results", "figures_nature")

CAPTURED: list[tuple[str, object]] = []


def _install_recorders():
    """Replace save()/close() so figures are kept for measurement only."""
    def make_rec(name):
        def _save(fig, outdir, fname, width_mm=None):
            CAPTURED.append((fname, fig))
        return _save

    N.save = make_rec("nature")
    M7.save = make_rec("main7")
    plt.close = lambda *a, **k: None


def _texts_of(ax):
    out = []
    for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        if t.get_text().strip():
            out.append(("tick", t))
    for t in ax.texts:
        if t.get_text().strip():
            out.append(("text", t))
    lg = ax.get_legend()
    if lg is not None:
        for t in lg.get_texts():
            if t.get_text().strip():
                out.append(("legend", t))
    out.append(("title", ax.title))
    out.append(("xlabel", ax.xaxis.label))
    out.append(("ylabel", ax.yaxis.label))
    return [(k, t) for k, t in out if t.get_text().strip()]


def _overlap_area(a, b):
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx * dy if (dx > 0 and dy > 0) else 0.0


def audit(name, fig, tol=0.22):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fw, fh = fig.canvas.get_width_height()
    problems = []

    def _is_panel_letter(kind, t):
        # panel letters are intentionally drawn just outside the axes
        s = t.get_text().strip()
        return kind == "text" and len(s) <= 2 and (s.isalnum() or s.startswith("S"))

    for ai, ax in enumerate(fig.axes):
        entries = _texts_of(ax)
        boxes = []
        for kind, t in entries:
            try:
                bb = t.get_window_extent(renderer=r)
            except Exception:
                continue
            boxes.append((kind, t, bb))
        # --- pairwise overlaps -----------------------------------------
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                k1, t1, b1 = boxes[i]
                k2, t2, b2 = boxes[j]
                if k1 == k2 == "legend":
                    continue
                inter = _overlap_area(b1, b2)
                sm = min(b1.width * b1.height, b2.width * b2.height)
                if sm <= 0:
                    continue
                frac = inter / sm
                if frac > tol:
                    problems.append(
                        f"ax{ai} OVERLAP {frac:.0%}: "
                        f"[{k1}] '{t1.get_text()[:26]}' <-> "
                        f"[{k2}] '{t2.get_text()[:26]}'")
        # --- out-of-axes / out-of-canvas --------------------------------
        ab = ax.get_window_extent(renderer=r)
        for kind, t, bb in boxes:
            if _is_panel_letter(kind, t):
                continue
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fw + 1 or bb.y1 > fh + 1:
                problems.append(
                    f"ax{ai} OFF-CANVAS [{kind}] '{t.get_text()[:26]}'")
            elif kind == "text" and (
                    bb.x0 < ab.x0 - 2 or bb.x1 > ab.x1 + 2
                    or bb.y0 < ab.y0 - 2 or bb.y1 > ab.y1 + 2):
                problems.append(
                    f"ax{ai} OUTSIDE-AXES [{kind}] '{t.get_text()[:26]}'")
    print(f"\n=== {name}  (panels={len(fig.axes)}, issues={len(problems)}) ===")
    for p in problems:
        print("   ", p)
    return len(problems)


def main() -> int:
    _install_recorders()
    N.set_style()
    for fn in (N.fig_overview, N.fig_human_blood, N.fig_mouse_sc, N.fig_cellcomm,
               N.fig_spatial, N.fig_crossspecies, N.fig_candidate_matrix,
               N.fig_feature_selection, N.fig_immune_series, N.fig_virtual_ko):
        try:
            fn(OUT)
        except Exception as e:  # noqa: BLE001
            print(f"[skip] {fn.__name__}: {type(e).__name__}: {e}")
    M7.set_style()
    for fn in (M7.fig1_framework, M7.fig6_immune_sheddase,
               M7.fig7_crossspecies_perturbation):
        try:
            fn(OUT)
        except Exception as e:  # noqa: BLE001
            print(f"[skip] {fn.__name__}: {type(e).__name__}: {e}")

    total = 0
    for name, fig in CAPTURED:
        try:
            total += audit(name, fig)
        except Exception as e:  # noqa: BLE001
            print(f"[audit-error] {name}: {type(e).__name__}: {e}")
    print(f"\n[total issues] {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
