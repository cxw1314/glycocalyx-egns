# glycocalyx-egns

Analysis code for the manuscript:

> **A transcriptionally defined endothelial glycocalyx non-resolution state links
> immune proteolysis to post-infarction fibrosis**

This repository contains every script needed to reproduce the module scores,
statistical tests, spatial analyses, communication inference, perturbation
analyses and figures reported in the manuscript. All input data are public; no
new data were generated.

---

## 1. Environment

* Python 3.12.13
* Pinned dependencies: see `requirements.txt`
* No R environment is required (the scTenifoldKnk-inspired analysis is a
  documented Python re-implementation, see Sec. 5 below)

```bash
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## 2. Input data

All datasets are obtained from NCBI GEO. Nothing is redistributed here.

| Accession | Role in the study |
|---|---|
| GSE176092 / GSE214611 | Primary murine MI snRNA-seq and Visium spatial data |
| GSE135310 | Time-series cardiac/blood immune scRNA-seq after murine MI |
| GSE206787 | Independent murine MI single-cell / positional context dataset |
| GSE183852 | Human donor and DCM cardiac single-cell/single-nucleus atlas |
| GSE123342 | Whole-blood transcriptomics after acute MI and follow-up |

Helpers `download.py`, `geo_files.py` and `untar_raw.py` fetch and unpack the
supplementary archives; `parse_series_matrix.py` and `get_platform_annot.py`
handle the microarray cohort (GSE123342, platform GPL17586).

## 3. Reproduction order

Run the scripts in the order listed in `docs/ANALYSIS_ORDER.md`. Each script
writes its tables under `results/<dataset>/tables/` and takes a `--outdir`
argument. Every stochastic step uses a fixed seed (`SEED = 42` unless stated in
the script header).

## 4. Figure generation

`docs/FIGURE_MAP.md` maps each main and supplementary figure to the script that
produces it. Figures are written as vector PDF plus 600-dpi PNG at Nature
double-column width (183 mm).

## 5. Analyses that are deliberately NOT official tool runs

Two methods are **documented re-implementations**, not runs of the cited
official packages. This is stated explicitly in the manuscript, and the
distinction matters for interpretation:

1. **Network edge-ablation perturbation** (inspired by scTenifoldKnk). The
   official R package is not executed; the same "delete the target gene's
   outgoing edges and compare networks" logic is implemented in Python.
   See `analysis_sctenifold_sheddase.py`.
2. **CellOracle-inspired vector-field sensitivity analysis.** Official
   CellOracle is not executed. This analysis **failed its own prespecified
   trajectory-validity check** (Spearman correlation between pseudotime and the
   target state was -0.145 for endothelium and -0.124 for macrophages, below the
   0.3 threshold), so its readout is not interpretable and it is reported only
   as a negative result in the Supplementary note.
   See `analysis_celloracle_offline.py`, which prints the validity check on every
   run and warns when it fails.

NicheNet-style multivariate ligand-target inference was **not** performed; all
communication results come from the LIANA consensus pipeline.

## 6. Licence and citation

Released under the MIT Licence, see `LICENSE.md`. Machine-readable
citation metadata is in `CITATION.cff`; please cite the associated manuscript.

* Repository: https://github.com/cxw1314/glycocalyx-egns
* Archived release DOI: assigned on archiving (see CITATION.cff)

## 7. Contact

Corresponding author: Xingwang Cao
