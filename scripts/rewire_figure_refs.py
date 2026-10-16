#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rewire the manuscripts onto the consolidated figure set (7 main + 3 supplementary).

For every `![图 N](path)` + following `**图 N｜...**` caption block:
  * renumber and repoint (kept figures),
  * delete (figures whose panels were merged into a composite), or
  * move into the Supplementary section (demoted heuristic figures).

Figures whose panels were merged are removed from the main text; the composite
that now carries them is inserted at a defined anchor.
"""
from __future__ import annotations

import argparse
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------- mapping ---
# old number -> ("keep", new_number, new_file, new_caption) | ("del",) | ("supp", ...)
ZH = {
    "img_tpl": "![图 {n}](../results/figures_nature/{f}.png)",
    "cap_tpl": "**图 {n}｜{c}**",
    "kept": {
        1: (1, "Fig1_study_framework",
            "研究设计与 EGNS / Core-EGNS 框架。** **a**，各数据集分析到的细胞/样本数（对数轴）。**b**，八个预定义模块在数据中检出的基因比例。**c**，四项主要效应量，含 HPSE 的内皮特异性比例（0.013）。**d, e**，描述性 EGNS 与 Core-EGNS 的逐细胞分布。**f**，按动物（n = 6）的 pseudobulk 均值，两种定义方向一致。**g**，30 个 panel TF 中有 14 个嵌入被评分模块（EIS 6、EAS 2、MRS 6）。"),
        2: (2, "Fig2_mouse_singlecell", None),
        11: (3, "Fig3_regulatory_network",
             "真实调控子活性推断（decoupler / CollecTRI）。** **(a)** 内皮细胞各条件 29 个 panel TF 的 regulon 活性（ULM）热图；**(b)** 人工靶基因代理与真实 regulon 活性的相关性（每个点为一个 TF）。"),
        3: (4, "Fig4_immune_cascade", None),
        4: (5, "Fig5_spatial_neighbourhood", None),
        12: (6, "Fig6_immune_sheddase_model",
             "细胞通讯与免疫来源 sheddase 效应层。** **(a)** LIANA 内皮接收端前 14 位配受体交互（颜色为共识秩）；**(b)** 各内置方法与共识秩的一致性；**(c)** 各 sheddase 敲除对 MFS/FNS/GDS/EAS/GCI 程序的效应比（>1 高于全基因组平均）；**(d)** 网络级重构幅度；**(e)** 候选分子的内皮表达比例、人血效应量与显著性。"),
        6: (7, "Fig7_crossspecies_perturbation",
            "跨物种一致性与计算扰动。** **(a)** 人 DCM 内皮效应量（红 = q < 0.05）；**(b)** 小鼠心梗与人 DCM 的方向一致性（绿 = 同向）；**(c)** 技术分层比较；**(d)** 内皮 TF 敲除对 EGNS 的影响；**(e)** 主要 TF 的模块级分解；**(f)** 巨噬细胞 TF 敲除对 MFS 的影响。"),
    },
    "deleted": [7, 10, 13, 14],
    "supp": {5: "S1", 8: "S2", 9: "S3"},
    "supp_header": "## 补充材料（Supplementary）",
    "supp_note": ("\n### 补充图图注\n"),
}

EN = {
    "img_tpl": "![Fig. {n}](../results/figures_nature/{f}.png)",
    "cap_tpl": "**Fig. {n} | {c}**",
    "kept": {
        1: (1, "Fig1_study_framework",
            "Study design and the EGNS / Core-EGNS framework. ** **a**, Cells or samples analysed per dataset (log axis). **b**, Fraction of each predefined module's genes detected. **c**, Headline effect sizes, including the endothelial-specificity fraction for HPSE (0.013). **d, e**, Per-cell distributions of the descriptive EGNS and Core-EGNS. **f**, Per-animal (n = 6) pseudobulk means; both definitions agree in direction. **g**, 14 of the 30 panel TFs are embedded in a scored module (EIS 6, EAS 2, MRS 6)."),
        2: (2, "Fig2_mouse_singlecell", None),
        11: (3, "Fig3_regulatory_network",
             "Resource-backed regulon inference (decoupler / CollecTRI). ** **(a)** ULM regulon activity for the 29 panel TFs across conditions in endothelial cells; **(b)** correlation between the curated-target proxy and the real regulon activity (one point per TF)."),
        3: (4, "Fig4_immune_cascade", None),
        4: (5, "Fig5_spatial_neighbourhood", None),
        12: (6, "Fig6_immune_sheddase_model",
             "Cell communication and the immune-derived sheddase effector layer. ** **(a)** Top 14 ligand-receptor interactions with an endothelial receiver (colour, consensus rank); **(b)** agreement of each built-in method with the consensus rank; **(c)** effect ratio of each sheddase knockout on the MFS/FNS/GDS/EAS/GCI programmes (>1 above the genome-wide mean); **(d)** magnitude of network-level rewiring; **(e)** endothelial expression fraction, human-blood effect size and significance of each candidate."),
        6: (7, "Fig7_crossspecies_perturbation",
            "Cross-species concordance and computational perturbation. ** **(a)** Effect sizes in human DCM endothelium (red, q < 0.05); **(b)** directional concordance between mouse MI and human DCM (green, concordant); **(c)** technology-stratified comparison; **(d)** effect of endothelial TF knockout on EGNS; **(e)** module-level decomposition for the leading TFs; **(f)** effect of macrophage TF knockout on MFS."),
    },
    "deleted": [7, 10, 13, 14],
    "supp": {5: "S1", 8: "S2", 9: "S3"},
    "supp_header": "## Supplementary",
    "supp_note": ("\n### Supplementary figure legends\n"),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", required=True)
    ap.add_argument("--lang", choices=["zh", "en"], required=True)
    args = ap.parse_args()

    L = ZH if args.lang == "zh" else EN
    fig_re = re.compile(
        r"!\[(?:图|Fig\.)\s*(\d+)\]\(([^)]+)\)[ \t]*\n+[ \t]*"
        r"\*\*(?:图|Fig\.)\s*\1\s*[｜|][^\n]*",
        re.MULTILINE)

    path = args.md if os.path.isabs(args.md) else os.path.join(ROOT, args.md)
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()

    supp_blocks = []
    stats = {"kept": 0, "deleted": 0, "supp": 0}

    def repl(m):
        old = int(m.group(1))
        block = m.group(0)
        # the caption is the LAST line of the block; take everything after the
        # "图 N｜" / "Fig. N |" prefix, and drop a single trailing bold marker.
        cap_line = block.rstrip().split("\n")[-1].strip()
        if "｜" in cap_line:
            cap_body = cap_line.split("｜", 1)[1]
        elif "|" in cap_line:
            cap_body = cap_line.split("|", 1)[1]
        else:
            cap_body = cap_line
        cap_body = cap_body.strip()
        if cap_body.endswith("**"):
            cap_body = cap_body[:-2].strip()

        prefix = "图 " if args.lang == "zh" else "Fig. "
        sep = "｜" if args.lang == "zh" else " | "

        if old in L["kept"]:
            new_n, new_f, new_cap = L["kept"][old]
            stats["kept"] += 1
            body = (new_cap if new_cap is not None else cap_body).strip()
            while body.endswith("*"):
                body = body[:-1].rstrip()
            return (L["img_tpl"].format(n=new_n, f=new_f) + "\n\n"
                    + f"**{prefix}{new_n}{sep}{body}**")
        if old in L["supp"]:
            tag = L["supp"][old]
            stats["supp"] += 1
            supp_blocks.append(
                L["img_tpl"].format(n=tag, f={
                    "S1": "FigS1_comm_heuristic",
                    "S2": "FigS2_feature_screening",
                    "S3": "FigS3_human_blood"}[tag]) + "\n\n"
                + f"**{prefix}{tag}{sep}{cap_body}**")
            return ""
        if old in L["deleted"]:
            stats["deleted"] += 1
            return ""
        return block

    text = fig_re.sub(repl, text)

    # ---- append the supplementary figures (forced S1, S2, S3 order) -----
    if supp_blocks:
        supp_blocks.sort(key=lambda b: ("S1" not in b, "S2" not in b))
        if L["supp_header"] in text:
            text = text.rstrip() + "\n" + L["supp_note"] + "\n" + \
                "\n\n".join(supp_blocks) + "\n"
        else:
            text = (text.rstrip() + "\n\n---\n\n" + L["supp_header"] + "\n"
                    + L["supp_note"] + "\n" + "\n\n".join(supp_blocks) + "\n")

    # tidy up triple blank lines left behind by deletions
    text = re.sub(r"\n{4,}", "\n\n\n", text)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"[ok] {path}")
    print(f"     kept={stats['kept']}  deleted(merged)={stats['deleted']}  "
          f"-> supplementary={stats['supp']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
