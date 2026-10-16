#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sanity-check generated figures: dimensions, ink coverage (non-blank), DPI."""
from __future__ import annotations

import glob
import os

import numpy as np
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> int:
    d = os.path.join(ROOT, "results", "figures_nature")
    files = sorted(glob.glob(os.path.join(d, "*.png")))
    if not files:
        print("[error] no figures found")
        return 1
    print(f"{'figure':34s} {'px':>13s} {'mm @600dpi':>14s} {'ink%':>7s} {'KB':>8s}")
    for f in files:
        im = Image.open(f).convert("L")
        a = np.asarray(im, dtype=np.uint8)
        ink = float((a < 240).mean() * 100)
        w_mm = im.width / 600 * 25.4
        h_mm = im.height / 600 * 25.4
        print(f"{os.path.basename(f):34s} {im.width:6d}x{im.height:<6d} "
              f"{w_mm:6.0f}x{h_mm:<6.0f} {ink:6.1f}% {os.path.getsize(f)/1024:7.0f}")
    print("\nNature widths: single column 89 mm, double column 183 mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
