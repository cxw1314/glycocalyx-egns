#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Verify a generated .docx: list embedded media, count images placed in the
document body, and flag any '[Image not found' placeholders.
"""
from __future__ import annotations

import os
import re
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(path: str) -> int:
    if not os.path.exists(path):
        print(f"[missing] {path}")
        return 1
    print(f"=== {os.path.basename(path)} ({os.path.getsize(path)/1024:.0f} KB)")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        media = [n for n in names if n.startswith("word/media/")]
        print(f"    embedded media files : {len(media)}")
        for m in sorted(media)[:12]:
            print(f"      {m}  ({z.getinfo(m).file_size/1024:.0f} KB)")
        if len(media) > 12:
            print(f"      ... and {len(media)-12} more")

        doc = z.read("word/document.xml").decode("utf-8", "replace")
        drawings = len(re.findall(r"<w:drawing>", doc))
        tbls = len(re.findall(r"<w:tbl>", doc))
        paras = len(re.findall(r"<w:p[ >]", doc))
        print(f"    inline drawings      : {drawings}")
        print(f"    tables               : {tbls}")
        print(f"    paragraphs           : {paras}")

        text = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", doc))
        ph = text.count("Image not found")
        print(f"    '[Image not found]'  : {ph}")
        for probe in ("摘要", "Abstract" if "english" in path.lower() else "§3.10",
                      "虚拟", "virtual" if "english" in path.lower() else "EGNS"):
            print(f"    contains '{probe}': {probe in text}")
    return 0


def main() -> int:
    for name in sys.argv[1:]:
        p = name if os.path.isabs(name) else os.path.join(ROOT, name)
        check(p)
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
