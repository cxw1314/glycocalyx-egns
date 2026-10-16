# Analysis order

## Main pipeline

1. `scripts/analysis_gse123342.py` - GSE123342 whole-blood module scoring and group testing (Sec. 3.1/3.9)
2. `scripts/analysis_gse123342_ml.py` - GSE123342 machine-learning screening and nested-CV leakage audit (Sec. 3.8/3.9)
3. `scripts/analysis_gse176092.py` - GSE176092 non-cardiomyocyte snRNA-seq: lineages, module scores, candidate specificity (Sec. 3.2/3.7)
4. `scripts/analysis_scrna.py` - Shared single-cell scoring / pseudobulk helpers
5. `scripts/analysis_gse135310.py` - GSE135310 cardiac immune time series and CCR2+ recruitment (Sec. 3.4)
6. `scripts/analysis_visium_spatial.py` - GSE176092 Visium per-spot scores, Moran's I and neighbourhood tests (Sec. 3.5)
7. `scripts/build_visium_h5ad.py` - Builds the concatenated Visium AnnData used by the deconvolution step
8. `scripts/analysis_spatial_deconvolution.py` - NNLS deconvolution and composition-adjusted coupling (Sec. 3.5.2)
9. `scripts/analysis_decoupler_regulons.py` - decoupler/CollecTRI ULM regulon activity (Sec. 3.3.1)
10. `scripts/analysis_regulon_concordance.py` - Concordance between heuristic proxy and real regulons (Sec. 3.3.1)
11. `scripts/analysis_liana_communication.py` - LIANA rank_aggregate consensus cell-cell communication (Sec. 3.5.1)
12. `scripts/analysis_lto_egns.py` - Core-EGNS / leave-target-out circularity analysis (Sec. 3.7)
13. `scripts/analysis_virtual_ko.py` - In silico TF knockout on the ridge GRN (Sec. 3.7)
14. `scripts/analysis_sctenifold_sheddase.py` - Network edge-ablation perturbation of non-TF sheddases (Sec. 3.7)
15. `scripts/analysis_celloracle_offline.py` - CellOracle-inspired vector-field sensitivity analysis, rejected by its own validity check (Supplementary note)
16. `scripts/analysis_gse183852.py` - GSE183852 human heart atlas cross-species analysis (Sec. 3.6)
17. `scripts/analysis_gse206787.py` - GSE206787 dataset audit and context analysis (Sec. 2.2)

## Supplementary / sensitivity analyses

1. `scripts/analysis_regulons.py` - Heuristic curated-target regulon proxy (Supplementary S1)
2. `scripts/analysis_cellcomm.py` - 26-pair heuristic ligand-receptor scoring (Supplementary S2)
3. `scripts/verify_docx.py` - Manuscript figure/table reference integrity check

## Figure and audit tooling

* `scripts/make_figures_nature.py` - Base Nature-style figure set
* `scripts/make_figures_v2.py` - Regulon / circularity figure panels
* `scripts/make_figures_main7.py` - Consolidated 7-figure main set (composite panels)
* `scripts/reorganise_figures.py` - Copy-then-verify-then-delete figure renaming
* `scripts/rewire_figure_refs.py` - Renumber / repoint manuscript figure references
* `scripts/qc2_demote_heuristics.py` - Move heuristic analyses into the supplement
* `scripts/audit_figure_legibility.py` - Programmatic text-overlap and clipping audit for every figure
* `scripts/audit_module_coverage.py` - Module gene-coverage audit with corrected denominators
* `scripts/final_check.py` - Final figure-reference integrity check
* `scripts/verify_docx_images.py` - Verify every referenced image is embedded in the DOCX

## Data acquisition helpers

* `scripts/download.py` - Resumable downloader for GEO supplementary files
* `scripts/geo_files.py` - List GEO supplementary files for an accession
* `scripts/geo_search.py` - GEO series search helper
* `scripts/geo_probe.py` - Probe a GEO accession's available files
* `scripts/get_platform_annot.py` - Download and parse a GEO platform annotation table
* `scripts/parse_series_matrix.py` - Parse GEO series matrix into an expression matrix
* `scripts/pfetch.py` - Parallel fetch helper
* `scripts/untar_raw.py` - Extract downloaded archives
* `scripts/extract_genes_stream.py` - Streaming gene extraction from very large matrices
* `scripts/list_gse183852_samples.py` - Map GSE183852 cells to samples
* `scripts/check_candidates.py` - Candidate-gene availability check
* `scripts/check_figures.py` - Figure file presence check
* `scripts/make_gene_list.py` - Build the candidate gene list
* `scripts/apply_annot.py` - Apply cluster annotations
* `scripts/show_virtual_ko.py` - Print virtual-knockout summary tables
