#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reorganise the figure set into 7 main figures + 3 supplementary figures.

Order of operations matters: copy-to-new-name FIRST, delete-old-name SECOND.
(An earlier PowerShell attempt deleted sources before the copies succeeded.)

Final main set
    Fig1_study_framework          dataseries + EGNS/Core-EGNS framework   [composite]
    Fig2_mouse_singlecell         endothelial temporal state
    Fig3_regulatory_network       decoupler/CollecTRI regulons   (was Fig11)
    Fig4_immune_cascade           cardiac immune time series     (was Fig3)
    Fig5_spatial_neighbourhood    Visium neighbourhood           (was Fig4)
    Fig6_immune_sheddase_model    LIANA + effector ranking       [composite]
    Fig7_crossspecies_perturbation cross-species + in silico KO  [composite]
Supplementary
    FigS1_comm_heuristic  (was Fig5)   FigS2_feature_screening (was Fig8)
    FigS3_human_blood     (was Fig9)
"""
from __future__ import annotations

import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "results", "figures_nature")

RENAMES = {
    "Fig11_decoupler_regulons": "Fig3_regulatory_network",
    "Fig3_immune_time_series": "Fig4_immune_cascade",
    "Fig4_spatial_neighbourhood": "Fig5_spatial_neighbourhood",
    "Fig5_cell_communication": "FigS1_comm_heuristic",
    "Fig8_feature_selection": "FigS2_feature_screening",
    "Fig9_human_blood": "FigS3_human_blood",
}
KEEP = ["Fig1_study_framework", "Fig2_mouse_singlecell",
        "Fig6_immune_sheddase_model", "Fig7_crossspecies_perturbation"]
RETIRE = ["Fig1_data_overview", "Fig3_immune_time_series",
          "Fig4_spatial_neighbourhood", "Fig5_cell_communication",
          "Fig6_cross_species", "Fig7_candidate_decision_matrix",
          "Fig8_feature_selection", "Fig9_human_blood",
          "Fig10_virtual_knockout", "Fig11_decoupler_regulons",
          "Fig12_LIANA_communication", "Fig13_EGNS_circularity",
          "Fig14_sheddase_KO"]


def main() -> int:
    missing = []
    # ---- 1. copy to the new names -------------------------------------
    for src, dst in RENAMES.items():
        for ext in ("pdf", "png"):
            s = os.path.join(FIG, f"{src}.{ext}")
            d = os.path.join(FIG, f"{dst}.{ext}")
            if not os.path.exists(s):
                missing.append(f"{src}.{ext}")
                continue
            shutil.copy2(s, d)
    if missing:
        print("[error] sources missing, not deleting anything:", missing)
        return 1

    # ---- 2. verify the final set exists BEFORE deleting ----------------
    want = KEEP + list(RENAMES.values())
    absent = [n for n in want
              if not all(os.path.exists(os.path.join(FIG, f"{n}.{e}"))
                         for e in ("pdf", "png"))]
    if absent:
        print("[error] final set incomplete, not deleting anything:", absent)
        return 1

    # ---- 3. retire merged / superseded panels -------------------------
    for n in RETIRE:
        for ext in ("pdf", "png"):
            p = os.path.join(FIG, f"{n}.{ext}")
            if os.path.exists(p):
                os.remove(p)

    print(f"[ok] main figures  ({len(KEEP) + len(RENAMES) - 3} expected -> 7):")
    for n in KEEP + ["Fig3_regulatory_network", "Fig4_immune_cascade",
                     "Fig5_spatial_neighbourhood"]:
        print(f"     {n}")
    print("[ok] supplementary: FigS1_comm_heuristic, FigS2_feature_screening, "
          "FigS3_human_blood")
    print("\n[files present]")
    for f in sorted(os.listdir(FIG)):
        if f.endswith(".png"):
            print("   ", f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
