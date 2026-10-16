#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Assemble a submission-ready code repository bundle for the EGNS manuscript.

Includes only the analysis / figure / audit scripts that actually produce the
manuscript's results. Internal diagnostics, debug probes and one-off
manuscript-rewriting helpers are excluded (see EXCLUDE_PREFIXES).

Outputs
    code_repository/glycocalyx-egns/     (browsable directory, git-init ready)
    code_repository/glycocalyx-egns.zip  (upload / Zenodo-ready archive)
"""
from __future__ import annotations

import json
import os
import shutil
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = ROOT
DEST_ROOT = os.path.join(ROOT, "code_repository")
REPO = os.path.join(DEST_ROOT, "glycocalyx-egns")

# --------------------------------------------------------------- release metadata
# Author-supplied values. They live here, not in the generated files, so that a
# rebuild cannot silently revert the released metadata back to placeholders.
REPO_URL = "https://github.com/cxw1314/glycocalyx-egns"
LICENSE_ID = "MIT"
LICENSE_YEAR = "2026"
AUTHORS = [("Xingwang", "Cao"), ("Ting", "Wang")]   # (given, family), in order
CORRESPONDING = "Xingwang Cao"
ZENODO_DOI = ""       # filled once the Zenodo archive mints a DOI
MANUSCRIPT_DOI = ""   # filled on acceptance

EXCLUDE_PREFIXES = ("_", "diag_", "debug_", "peek_")
EXCLUDE_NAMES = {
    "_diag_api.py", "_diag_cov.py", "_diag_liana.py", "_diag_liana_cols.py",
    "_diag_sets.py", "_final_state.py", "_inspect_jmcc.py", "_inspect_jmcc2.py",
    "_rewrite_fig7.py", "_trim_orphan.py", "_check_captions.py",
    "trim_manuscript.py",  # one-off manuscript rewriting helper
}

# main-text analysis pipeline, in execution order
PIPELINE_MAIN = [
    ("analysis_gse123342.py", "GSE123342 whole-blood module scoring and group testing (Sec. 3.1/3.9)"),
    ("analysis_gse123342_ml.py", "GSE123342 machine-learning screening and nested-CV leakage audit (Sec. 3.8/3.9)"),
    ("analysis_gse176092.py", "GSE176092 non-cardiomyocyte snRNA-seq: lineages, module scores, candidate specificity (Sec. 3.2/3.7)"),
    ("analysis_scrna.py", "Shared single-cell scoring / pseudobulk helpers"),
    ("analysis_gse135310.py", "GSE135310 cardiac immune time series and CCR2+ recruitment (Sec. 3.4)"),
    ("analysis_visium_spatial.py", "GSE176092 Visium per-spot scores, Moran's I and neighbourhood tests (Sec. 3.5)"),
    ("build_visium_h5ad.py", "Builds the concatenated Visium AnnData used by the deconvolution step"),
    ("analysis_spatial_deconvolution.py", "NNLS deconvolution and composition-adjusted coupling (Sec. 3.5.2)"),
    ("analysis_decoupler_regulons.py", "decoupler/CollecTRI ULM regulon activity (Sec. 3.3.1)"),
    ("analysis_regulon_concordance.py", "Concordance between heuristic proxy and real regulons (Sec. 3.3.1)"),
    ("analysis_liana_communication.py", "LIANA rank_aggregate consensus cell-cell communication (Sec. 3.5.1)"),
    ("analysis_lto_egns.py", "Core-EGNS / leave-target-out circularity analysis (Sec. 3.7)"),
    ("analysis_virtual_ko.py", "In silico TF knockout on the ridge GRN (Sec. 3.7)"),
    ("analysis_sctenifold_sheddase.py", "Network edge-ablation perturbation of non-TF sheddases (Sec. 3.7)"),
    ("analysis_celloracle_offline.py", "CellOracle-inspired vector-field sensitivity analysis, rejected by its own validity check (Supplementary note)"),
    ("analysis_gse183852.py", "GSE183852 human heart atlas cross-species analysis (Sec. 3.6)"),
    ("analysis_gse206787.py", "GSE206787 dataset audit and context analysis (Sec. 2.2)"),
]

# supplementary / heuristic sensitivity analyses
PIPELINE_SUPP = [
    ("analysis_regulons.py", "Heuristic curated-target regulon proxy (Supplementary S1)"),
    ("analysis_cellcomm.py", "26-pair heuristic ligand-receptor scoring (Supplementary S2)"),
    ("verify_docx.py", "Manuscript figure/table reference integrity check"),
]

# figure and audit tooling
FIG_TOOLS = [
    ("make_figures_nature.py", "Base Nature-style figure set"),
    ("make_figures_v2.py", "Regulon / circularity figure panels"),
    ("make_figures_main7.py", "Consolidated 7-figure main set (composite panels)"),
    ("reorganise_figures.py", "Copy-then-verify-then-delete figure renaming"),
    ("rewire_figure_refs.py", "Renumber / repoint manuscript figure references"),
    ("qc2_demote_heuristics.py", "Move heuristic analyses into the supplement"),
    ("audit_figure_legibility.py", "Programmatic text-overlap and clipping audit for every figure"),
    ("audit_module_coverage.py", "Module gene-coverage audit with corrected denominators"),
    ("final_check.py", "Final figure-reference integrity check"),
    ("verify_docx_images.py", "Verify every referenced image is embedded in the DOCX"),
]

DATA_HELPERS = [
    ("download.py", "Resumable downloader for GEO supplementary files"),
    ("geo_files.py", "List GEO supplementary files for an accession"),
    ("geo_search.py", "GEO series search helper"),
    ("geo_probe.py", "Probe a GEO accession's available files"),
    ("get_platform_annot.py", "Download and parse a GEO platform annotation table"),
    ("parse_series_matrix.py", "Parse GEO series matrix into an expression matrix"),
    ("pfetch.py", "Parallel fetch helper"),
    ("untar_raw.py", "Extract downloaded archives"),
    ("extract_genes_stream.py", "Streaming gene extraction from very large matrices"),
    ("list_gse183852_samples.py", "Map GSE183852 cells to samples"),
    ("check_candidates.py", "Candidate-gene availability check"),
    ("check_figures.py", "Figure file presence check"),
    ("make_gene_list.py", "Build the candidate gene list"),
    ("apply_annot.py", "Apply cluster annotations"),
    ("show_virtual_ko.py", "Print virtual-knockout summary tables"),
]

README = """# glycocalyx-egns

Analysis code for the manuscript:

> **A transcriptionally defined endothelial glycocalyx non-resolution state links
> immune proteolysis to post-infarction fibrosis**

This repository contains every script needed to reproduce the module scores,
statistical tests, spatial analyses, communication inference, perturbation
analyses and figures reported in the manuscript. All input data are public; no
new data were generated.

---

## 1. Environment

* Python {pyver}
* Pinned dependencies: see `requirements.txt`
* No R environment is required (the scTenifoldKnk-inspired analysis is a
  documented Python re-implementation, see Sec. 5 below)

```bash
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\\Scripts\\activate
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

Released under the {license} Licence, see `LICENSE.md`. Machine-readable
citation metadata is in `CITATION.cff`; please cite the associated manuscript.

* Repository: {repourl}
* Archived release DOI: {doi}

## 7. Contact

Corresponding author: {corresponding}
"""

MIT_LICENCE = """MIT License

Copyright (c) {year} {holders}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


def _filter_scripts(names: list[str]) -> list[str]:
    out = []
    for n in sorted(names):
        if not n.endswith(".py"):
            continue
        if n in EXCLUDE_NAMES:
            continue
        if any(n.startswith(p) for p in EXCLUDE_PREFIXES):
            continue
        out.append(n)
    return out


def main() -> int:
    import importlib.metadata as md
    import sys

    if os.path.isdir(DEST_ROOT):
        shutil.rmtree(DEST_ROOT)
    os.makedirs(REPO)

    # ---- scripts -------------------------------------------------------
    src_scripts = os.path.join(SRC, "scripts")
    available = set(os.listdir(src_scripts))
    keep = _filter_scripts(list(available))
    os.makedirs(os.path.join(REPO, "scripts"))
    for n in keep:
        shutil.copy2(os.path.join(src_scripts, n), os.path.join(REPO, "scripts", n))
    print(f"[scripts] copied {len(keep)} analysis scripts "
          f"(excluded {len(available) - len(keep)} diagnostics/helpers)")

    # ---- config --------------------------------------------------------
    src_cfg = os.path.join(SRC, "config")
    os.makedirs(os.path.join(REPO, "config"), exist_ok=True)
    ncfg = 0
    for n in sorted(os.listdir(src_cfg)):
        p = os.path.join(src_cfg, n)
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(REPO, "config", n))
            ncfg += 1
    print(f"[config] copied {ncfg} files")

    # ---- requirements --------------------------------------------------
    pkgs = ["numpy", "pandas", "scipy", "scikit-learn", "anndata", "scanpy",
            "matplotlib", "seaborn", "statsmodels", "h5py", "networkx",
            "liana", "decoupler", "adjustText"]
    lines = ["# Pinned to the environment used for the manuscript",
             f"# Python {sys.version.split()[0]}", ""]
    for p in pkgs:
        try:
            lines.append(f"{p}=={md.version(p)}")
        except Exception:
            lines.append(f"# {p}: not installed in the reference environment")
    with open(os.path.join(REPO, "requirements.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- docs ----------------------------------------------------------
    os.makedirs(os.path.join(REPO, "docs"))
    with open(os.path.join(REPO, "docs", "ANALYSIS_ORDER.md"), "w",
              encoding="utf-8") as fh:
        fh.write("# Analysis order\n\n## Main pipeline\n\n")
        for i, (s, d) in enumerate(PIPELINE_MAIN, 1):
            fh.write(f"{i}. `scripts/{s}` - {d}\n")
        fh.write("\n## Supplementary / sensitivity analyses\n\n")
        for i, (s, d) in enumerate(PIPELINE_SUPP, 1):
            fh.write(f"{i}. `scripts/{s}` - {d}\n")
        fh.write("\n## Figure and audit tooling\n\n")
        for s, d in FIG_TOOLS:
            fh.write(f"* `scripts/{s}` - {d}\n")
        fh.write("\n## Data acquisition helpers\n\n")
        for s, d in DATA_HELPERS:
            fh.write(f"* `scripts/{s}` - {d}\n")

    figmap = {
        "Figure 1": ("make_figures_main7.py fig1_framework",
                     "Study design, dataset/module coverage, headline effects, EGNS vs Core-EGNS"),
        "Figure 2": ("make_figures_nature.py fig_mouse_sc",
                     "Endothelial temporal state and candidate cell specificity"),
        "Figure 3": ("make_figures_v2.py fig13 -> Fig11_decoupler_regulons",
                     "decoupler/CollecTRI regulon activity and proxy concordance"),
        "Figure 4": ("make_figures_nature.py fig_immune_series",
                     "Cardiac immune time series and EGNS alignment"),
        "Figure 5": ("make_figures_nature.py fig_spatial",
                     "Visium spatial autocorrelation and neighbourhood coupling"),
        "Figure 6": ("make_figures_main7.py fig6_immune_sheddase",
                     "LIANA endothelial-receiver interactions, sheddase ranking, candidate evidence"),
        "Figure 7": ("make_figures_main7.py fig7_crossspecies_perturbation",
                     "Cross-species concordance and the integrated mechanistic model"),
        "Figure S1": ("make_figures_nature.py fig_cellcomm",
                      "26-pair heuristic communication (sensitivity analysis)"),
        "Figure S2": ("make_figures_nature.py fig_feature_selection",
                      "Three-method screening and selection-leakage audit"),
        "Figure S3": ("make_figures_nature.py fig_human_blood",
                      "Whole-blood module scores, recovery trajectories, discrimination"),
    }
    with open(os.path.join(REPO, "docs", "FIGURE_MAP.md"), "w",
              encoding="utf-8") as fh:
        fh.write("# Figure to script map\n\n| Figure | Generator | Content |\n|---|---|---|\n")
        for k, (gen, desc) in figmap.items():
            fh.write(f"| {k} | `{gen}` | {desc} |\n")
        fh.write("\n`scripts/reorganise_figures.py` performs the final "
                 "file renaming (copy -> verify -> delete) and "
                 "`scripts/rewire_figure_refs.py` renumbers the captions in the "
                 "manuscript. `scripts/audit_figure_legibility.py` renders every "
                 "figure and flags overlapping or clipped text.\n")

    # ---- repo metadata --------------------------------------------------
    with open(os.path.join(REPO, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(README
                 .replace("{pyver}", sys.version.split()[0])
                 .replace("{license}", LICENSE_ID)
                 .replace("{repourl}", REPO_URL)
                 .replace("{doi}", ZENODO_DOI or
                          "assigned on archiving (see CITATION.cff)")
                 .replace("{corresponding}", CORRESPONDING))

    holders = ", ".join(f"{g} {f}" for g, f in AUTHORS)
    with open(os.path.join(REPO, "LICENSE.md"), "w", encoding="utf-8") as fh:
        fh.write(MIT_LICENCE.format(year=LICENSE_YEAR, holders=holders))

    # CFF is schema-validated (GitHub rejects the file outright if a field format
    # is wrong), so never emit a placeholder string for `doi` or `license`:
    # omit the field until a real value exists.
    cff = [
        "cff-version: 1.2.0",
        'message: "If you use this code, please cite the manuscript below."',
        'title: "glycocalyx-egns"',
        "abstract: >-",
        "  Analysis code for a transcriptionally defined endothelial glycocalyx",
        "  non-resolution state linking immune proteolysis to post-infarction",
        "  fibrosis.",
        "type: software",
        'version: "1.0.0"',
        f'date-released: "{LICENSE_YEAR}-10-16"',
        f"license: {LICENSE_ID}",
        f'repository-code: "{REPO_URL}"',
    ]
    if ZENODO_DOI:
        cff.append(f'doi: "{ZENODO_DOI}"')
    cff.append("authors:")
    for g, f in AUTHORS:
        cff.append(f'  - family-names: "{f}"')
        cff.append(f'    given-names: "{g}"')
    cff += [
        "preferred-citation:",
        "  type: article",
        "  title: >-",
        "    A transcriptionally defined endothelial glycocalyx non-resolution state",
        "    links immune proteolysis to post-infarction fibrosis",
        '  journal: "Journal of Molecular and Cellular Cardiology"',
        "  year: 2026",
    ]
    if MANUSCRIPT_DOI:
        cff.append(f'  doi: "{MANUSCRIPT_DOI}"')
    with open(os.path.join(REPO, "CITATION.cff"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(cff) + "\n")

    with open(os.path.join(REPO, ".zenodo.json"), "w", encoding="utf-8") as fh:
        json.dump({
            "title": "glycocalyx-egns: analysis code for endothelial glycocalyx "
                     "non-resolution after myocardial infarction",
            "description": "Analysis scripts reproducing all module scores, "
                           "statistical tests, spatial analyses, communication "
                           "inference, perturbation analyses and figures of the "
                           "manuscript. Two methods are documented "
                           "re-implementations rather than official tool runs "
                           "(network edge-ablation inspired by scTenifoldKnk; a "
                           "CellOracle-inspired vector-field analysis that failed "
                           "its own trajectory-validity check and is reported only "
                           "as a negative result).",
            "upload_type": "software",
            "access_right": "open",
            "license": LICENSE_ID,
            "version": "1.0.0",
            "language": "eng",
            "creators": [{"name": f"{f}, {g}"} for g, f in AUTHORS],
            "keywords": ["endothelial glycocalyx", "myocardial infarction",
                         "single-cell RNA sequencing", "spatial transcriptomics",
                         "fibrosis", "MMP9"],
            "related_identifiers": (
                [{"identifier": MANUSCRIPT_DOI,
                  "relation": "isSupplementTo",
                  "scheme": "doi",
                  "resource_type": "publication-article"}]
                if MANUSCRIPT_DOI else []),
        }, fh, indent=2, ensure_ascii=False)

    with open(os.path.join(REPO, "VERSION"), "w", encoding="utf-8") as fh:
        fh.write("1.0.0\n")

    with open(os.path.join(REPO, ".gitignore"), "w", encoding="utf-8") as fh:
        fh.write("\n".join([
            "__pycache__/", "*.py[cod]", ".venv/", "venv/",
            "# large inputs / derived matrices are fetched, not committed",
            "data/", "results/", "*.h5ad", "*.h5", "*.mtx", "*.mtx.gz",
            "figures/", "code_repository/", ""]))

    # ---- zip -----------------------------------------------------------
    zpath = os.path.join(DEST_ROOT, "glycocalyx-egns.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for base, _dirs, files in os.walk(REPO):
            for f in files:
                full = os.path.join(base, f)
                z.write(full, os.path.join("glycocalyx-egns",
                                           os.path.relpath(full, REPO)))
    print(f"[zip] {zpath} ({os.path.getsize(zpath)/1024:.0f} KB)")

    print("\n=== repository contents ===")
    for base, _dirs, files in os.walk(REPO):
        rel = os.path.relpath(base, REPO)
        for f in sorted(files):
            print(f"  {os.path.join(rel, f) if rel != '.' else f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
