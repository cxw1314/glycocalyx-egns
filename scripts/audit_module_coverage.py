#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Audit module gene coverage: reported 'n matched / n defined' can exceed the
defined count because alias genes are appended on top of the canonical list.
This script separates canonical matches from alias-rescued matches.
"""
from __future__ import annotations

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = os.path.join(ROOT, "config", "modules.json")


def main() -> int:
    with open(MODULES, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    print(f"{'module':8s} {'defined':>8s} {'aliases':>8s} {'max':>5s}")
    for key, val in cfg.items():
        if key.startswith("_") or key == "EGNS" or not isinstance(val, dict):
            continue
        genes = list(val.get("genes", []))
        alts = val.get("aliases") or {}
        n_alias = sum(len(v) for v in alts.values())
        print(f"{key:8s} {len(genes):8d} {n_alias:8d} {len(genes)+n_alias:5d}")
    print("\nNOTE: the earlier 'MRS 19/18' report counted canonical + alias "
          "matches against the canonical denominator, which is why the "
          "numerator can exceed the denominator.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
