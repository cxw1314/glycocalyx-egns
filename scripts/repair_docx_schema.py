#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Repair three PRE-EXISTING OOXML schema defects in the unpacked JMCC manuscript.

These defects are present in the author-supplied manuscript itself (verified by
validating the untouched upload and the edit snapshot: identical 3 errors, and
word/settings.xml is byte-identical before and after our scoped text edit).
They are unrelated to the code-availability change. This repair is therefore a
SEPARATELY SCOPED operation, declared as such.

Defects and minimal fixes
  1. <w:shd w:fill="..."/>            -> add required w:val="clear"
  2. <w:tblBorders> after <w:tblLook> -> move before <w:tblLayout>
     (CT_TblPrBase sequence: tblStyle, tblW, jc, tblBorders, tblLayout, tblLook)
  3. <w:zoom w:val="bestFit"/>        -> add required w:percent="100"

Usage:
    python repair_docx_schema.py <unpacked_dir>
"""
from __future__ import annotations

import os
import re
import sys

SHADING_RE = re.compile(r"<w:shd(?![^>]*\bw:val=)([^>]*?)/>")
TBLPR_RE = re.compile(r"<w:tblPr\b.*?</w:tblPr>", re.S)
TBLBORDERS_RE = re.compile(r"\s*<w:tblBorders\b.*?</w:tblBorders>", re.S)
TBLBORDERS_SELFCLOSE_RE = re.compile(r"\s*<w:tblBorders\b[^>]*/>", re.S)
ZOOM_RE = re.compile(r"<w:zoom(?![^>]*\bw:percent=)([^>]*?)/>")


def fix_shading(xml: str) -> tuple[str, int]:
    n = 0

    def sub(m):
        nonlocal n
        n += 1
        return f'<w:shd w:val="clear"{m.group(1)}/>'

    return SHADING_RE.sub(sub, xml), n


def fix_zoom(xml: str) -> tuple[str, int]:
    n = 0

    def sub(m):
        nonlocal n
        n += 1
        return f'<w:zoom m_placeholder{m.group(1)}/>'.replace("m_placeholder", 'w:percent="100" ')

    return ZOOM_RE.sub(sub, xml), n


def fix_tblborders(xml: str) -> tuple[str, int]:
    n = 0

    def fix_block(m):
        nonlocal n
        block = m.group(0)
        b = TBLBORDERS_RE.search(block) or TBLBORDERS_SELFCLOSE_RE.search(block)
        if not b:
            return block
        borders = b.group(0).strip()
        # already correctly placed (before tblLayout)?
        rest = block[: b.start()] + block[b.end():]
        if "<w:tblLayout" not in rest:
            return block
        n += 1
        rest = TBLBORDERS_RE.sub("", rest)
        rest = TBLBORDERS_SELFCLOSE_RE.sub("", rest)
        m2 = re.search(r"<w:tblLayout", rest)
        return rest[: m2.start()] + "  " + borders + "\n  " + rest[m2.start():]

    return TBLPR_RE.sub(fix_block, xml), n


def main() -> int:
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    doc = os.path.join(root, "word", "document.xml")
    setx = os.path.join(root, "word", "settings.xml")

    with open(doc, "r", encoding="utf-8") as fh:
        xml = fh.read()
    xml, n_shd = fix_shading(xml)
    xml, n_tbl = fix_tblborders(xml)
    with open(doc, "w", encoding="utf-8") as fh:
        fh.write(xml)

    with open(setx, "r", encoding="utf-8") as fh:
        st = fh.read()
    st, n_zoom = fix_zoom(st)
    with open(setx, "w", encoding="utf-8") as fh:
        fh.write(st)

    print(f"[repair] w:shd  + w:val  -> {n_shd}")
    print(f"[repair] w:tblBorders moved -> {n_tbl}")
    print(f"[repair] w:zoom + w:percent -> {n_zoom}")

    # residual check
    res = []
    if SHADING_RE.search(xml):
        res.append("shd still missing val")
    if TBLBORDERS_RE.search(xml) and any(
            TBLBORDERS_RE.search(b) and "<w:tblLayout" in b and
            (TBLBORDERS_RE.search(b).start() > b.find("<w:tblLayout"))
            for b in TBLPR_RE.findall(xml)):
        res.append("tblBorders still misordered")
    if ZOOM_RE.search(st):
        res.append("zoom still missing percent")
    print("[residual]", res or "none")
    return 1 if res else 0


if __name__ == "__main__":
    raise SystemExit(main())
