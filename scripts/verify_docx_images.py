#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Verify that a generated .docx embeds every referenced PNG (no placeholders)."""
from __future__ import annotations

import re
import sys
import zipfile


def check(docx: str, md: str) -> int:
    with open(md, "r", encoding="utf-8") as fh:
        text = fh.read()
    refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    with zipfile.ZipFile(docx) as z:
        names = z.namelist()
        media = [n for n in names if n.startswith("word/media/")]
        doc = z.read("word/document.xml").decode("utf-8", "ignore")
    missing = doc.count("Image not found")
    print(f"[{docx}]")
    print(f"  markdown image refs : {len(refs)}")
    print(f"  embedded media parts: {len(media)}")
    print(f"  'Image not found'   : {missing}")
    print(f"  total docx parts    : {len(names)}")
    ok = (missing == 0) and (len(media) >= len(refs))
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    rc = 0
    for pair in sys.argv[1:]:
        docx, md = pair.split("::")
        rc |= check(docx, md)
    raise SystemExit(rc)
