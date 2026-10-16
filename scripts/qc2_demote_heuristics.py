#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
QC-2: demote the v1.1 heuristic approximations out of the main text.

Moves
  * the curated-target regulon proxy block in Section 3.3, and
  * the 26-pair heuristic communication block in Section 3.5(5)
into a new "Supplementary" section at the end of the manuscript, leaving only
the resource-backed results (decoupler/CollecTRI, LIANA) in the main text.

Run:
    python qc2_demote_heuristics.py --md manuscript/论文_中文初稿.md --lang zh
    python qc2_demote_heuristics.py --md manuscript/paper_english_v1.md --lang en
"""
from __future__ import annotations

import argparse
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ZH = {
    "s33_start": "### 3.3 高 EGNS 内皮的调控网络",
    "s33_next": "### 3.3.1 真实调控子活性推断（decoupler / CollecTRI）",
    "s35_start": "**（5）配体—受体通讯：内皮受到的旁分泌信号主要来自巨噬细胞**",
    "s35_next": "### 3.5.1 真实 LIANA 通讯推断及其与启发式评分的差异",
    "supp_header": "## 补充材料（Supplementary）",
    "supp_note": ("> 本补充部分收录 v1.1 的启发式近似结果，仅作为**敏感性分析 / 方法演进记录**，"
                  "不构成正文结论依据。正文结论一律基于资源化真实工具（decoupler/CollecTRI、LIANA）。\n"),
    "s33_pointer": (
        "### 3.3.1 真实调控子活性推断（decoupler / CollecTRI）"),
    "s33_summary": (
        "> **方法演进说明**：v1.1 曾以\"人工整理 6–13 个靶基因取均值\"近似调控子活性"
        "（该近似结果见补充材料 S1）。经方法学审读，正文改用资源化真实推断（下述 §3.3.1）。"
        "两者的定量差异一并报告于 S1。\n"),
    "s35_summary": (
        "> **方法演进说明**：v1.1 曾以\"表达比例 × 平均表达量\"的启发式对 26 对配受体打分"
        "（该近似结果见补充材料 S2）。正文改用 LIANA 共识流程（下述 §3.5.1）；"
        "两者在内皮接收端领先轴上的差异一并报告于 S2。\n"),
    "supp_s1": "### S1 人工整理调控子代理（v1.1，敏感性分析）",
    "supp_s2": "### S2 26 对配受体启发式通讯（v1.1，敏感性分析）",
}

EN = {
    "s33_start": "### 3.3 Regulatory network of high-EGNS endothelium",
    "s33_next": "### 3.3.1 Resource-backed regulon inference (decoupler / CollecTRI)",
    "s35_start": ("**(5) Ligand-receptor communication: paracrine input to endothelium "
                  "comes mainly from macrophages**"),
    "s35_next": ("### 3.5.1 Real LIANA communication and how it differs "
                 "from the heuristic"),
    "supp_header": "## Supplementary",
    "supp_note": ("> This supplementary part collects the v1.1 heuristic approximations as "
                  "**sensitivity analysis / method-evolution records only**; they do not "
                  "support any main-text conclusion. All main-text conclusions rest on the "
                  "resource-backed tools (decoupler/CollecTRI, LIANA).\n"),
    "s33_summary": (
        "> **Method-evolution note**: v1.1 approximated regulon activity by the mean of "
        "6-13 curated target genes (results in Supplementary S1). Following methodological "
        "review, the main text now uses resource-backed inference (Section 3.3.1 below); "
        "the quantitative difference is reported in S1.\n"),
    "s35_summary": (
        "> **Method-evolution note**: v1.1 scored 26 ligand-receptor pairs with an "
        "\"expression fraction x mean expression\" heuristic (results in Supplementary S2). "
        "The main text now uses the LIANA consensus pipeline (Section 3.5.1 below); the "
        "difference in the leading endothelial-receiver axes is reported in S2.\n"),
    "supp_s1": "### S1 Curated-target regulon proxy (v1.1, sensitivity analysis)",
    "supp_s2": "### S2 26-pair heuristic communication (v1.1, sensitivity analysis)",
}


def find_block(text: str, start: str, end: str) -> tuple[int, int]:
    i = text.find(start)
    if i < 0:
        raise SystemExit(f"[error] start anchor not found: {start[:60]}")
    j = text.find(end, i + len(start))
    if j < 0:
        raise SystemExit(f"[error] end anchor not found: {end[:60]}")
    return i, j


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", required=True)
    ap.add_argument("--lang", choices=["zh", "en"], required=True)
    args = ap.parse_args()

    L = ZH if args.lang == "zh" else EN
    path = args.md if os.path.isabs(args.md) else os.path.join(ROOT, args.md)
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()

    # ---- cut block 1 (curated regulon proxy in 3.3) ---------------------
    a1, b1 = find_block(text, L["s33_start"], L["s33_next"])
    block1 = text[a1:b1].rstrip() + "\n"
    text = text[:a1] + L["s33_summary"] + "\n" + L["s33_next"] + text[b1 + len(L["s33_next"]):]

    # ---- cut block 2 (heuristic communication in 3.5(5)) ----------------
    a2, b2 = find_block(text, L["s35_start"], L["s35_next"])
    block2 = text[a2:b2].rstrip() + "\n"
    text = text[:a2] + L["s35_summary"] + "\n" + L["s35_next"] + text[b2 + len(L["s35_next"]):]

    # ---- append supplementary ------------------------------------------
    supp = ("\n\n---\n\n" + L["supp_header"] + "\n\n" + L["supp_note"] + "\n"
            + L["supp_s1"] + "\n\n" + block1 + "\n\n" + L["supp_s2"] + "\n\n" + block2 + "\n")
    text = text.rstrip() + supp

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)

    print(f"[ok] {path}")
    print(f"     moved block1: {len(block1)} chars  (curated regulon proxy)")
    print(f"     moved block2: {len(block2)} chars  (26-pair heuristic)")
    print(f"     supplementary appended: {len(supp)} chars")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
