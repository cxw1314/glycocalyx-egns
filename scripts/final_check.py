#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Final consistency check: every figure referenced by the manuscripts exists."""
from __future__ import annotations

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = os.path.join(ROOT, "manuscript")

DOCS = [("论文_中文初稿.md", "中文稿"), ("paper_english_v1.md", "English")]


def main() -> int:
    rc = 0
    for fname, label in DOCS:
        path = os.path.join(MAN, fname)
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
        refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
        missing = []
        for r in refs:
            p = os.path.normpath(os.path.join(MAN, r))
            if not os.path.exists(p):
                missing.append(r)
        n_sec = len(re.findall(r"^### ", text, flags=re.M))
        print(f"[{label}] {fname}")
        print(f"   figures referenced : {len(refs)}")
        print(f"   subsections (###)  : {n_sec}")
        print(f"   missing image files: {len(missing)}"
              + (f" -> {missing}" if missing else ""))
        rc |= 1 if missing else 0
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
