# Figure to script map

| Figure | Generator | Content |
|---|---|---|
| Figure 1 | `make_figures_main7.py fig1_framework` | Study design, dataset/module coverage, headline effects, EGNS vs Core-EGNS |
| Figure 2 | `make_figures_nature.py fig_mouse_sc` | Endothelial temporal state and candidate cell specificity |
| Figure 3 | `make_figures_v2.py fig13 -> Fig11_decoupler_regulons` | decoupler/CollecTRI regulon activity and proxy concordance |
| Figure 4 | `make_figures_nature.py fig_immune_series` | Cardiac immune time series and EGNS alignment |
| Figure 5 | `make_figures_nature.py fig_spatial` | Visium spatial autocorrelation and neighbourhood coupling |
| Figure 6 | `make_figures_main7.py fig6_immune_sheddase` | LIANA endothelial-receiver interactions, sheddase ranking, candidate evidence |
| Figure 7 | `make_figures_main7.py fig7_crossspecies_perturbation` | Cross-species concordance and the integrated mechanistic model |
| Figure S1 | `make_figures_nature.py fig_cellcomm` | 26-pair heuristic communication (sensitivity analysis) |
| Figure S2 | `make_figures_nature.py fig_feature_selection` | Three-method screening and selection-leakage audit |
| Figure S3 | `make_figures_nature.py fig_human_blood` | Whole-blood module scores, recovery trajectories, discrimination |

`scripts/reorganise_figures.py` performs the final file renaming (copy -> verify -> delete) and `scripts/rewire_figure_refs.py` renumbers the captions in the manuscript. `scripts/audit_figure_legibility.py` renders every figure and flags overlapping or clipped text.
